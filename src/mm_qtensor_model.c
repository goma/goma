/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 
/*
 * $Id: mm_qtensor_model.c,v 5.2 2009-05-20 15:31:33 hkmoffa Exp $
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_eh.h"
#include "el_elm_info.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_numjac.h"
#include "mm_viscosity.h"

#define GOMA_MM_QTENSOR_MODEL_C
#include "mm_qtensor_model.h"

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
*
*       NAME                    TYPE            CALLED_BY (almost certainly not up to date!)
*   ------------              ---------       --------------
* hydro_qtensor_flux             int          hydro_flux
* cross_really_simple_vectors   void          various in this file
* normalize_really_simple_vector dbl          various in this file
* compute_principle_directions  void          hydro_qtensor_flux
* diagonalize_tensor             int          hydro_qtensor_flux, comptue_VQVt_directly
* find_eigenvector              void          compute_principle_shear_directions
* compute_VQVt_directly         void          hydro_qtensor_flux
* please_work                   void          hydro_qtensor_flux
* bias_eigenvector              void          assemble_vorticity_direction
* get_nodal_qtensor		 int          assemble_new_qtensor
*    
******************************************************************************/

/* This file contains the routines to implement the Q-tensor diffusive
 * flux model.  It is a particular sub-model in the HYDRO Diffusivity
 * model, which is part of the HYDRODYNAMIC constituitive equation.
 * Its requirement/context overlaps substantially with the other HYDRO
 * submodels.  In particular, the Q-tensor diffusivity require both
 * the coefficients from the shear-rate and the viscosity diffusivity
 * submodels.
 *
 * It also requires the VORT_DIR{k} and VORT_LAMBDA
 * equations/variables to be active.  Since the vorticity principle
 * flow direction is always a 3-vector, there should be "3"'s all over
 * the place instead of DIM's.  It just plain doesn't make any sense
 * for anything other than 3-vectors.
 */

//static int bias_eigenvector_to(dbl *, dbl *);

static int get_local_qtensor
(double [][DIM]);

#define MAX_GAUSS_POINTS 12

/*
 *  Expanded Guass Point Quadrature Rules for Use in this
 *  Routine Only.
 *
 * Locations (x_i) for Gaussian quadrature,
 *    \int_{-1}^1 f(x) \approx \sum_{i} w_i f(x_i)
 */
static double gauss_point[MAX_GAUSS_POINTS][MAX_GAUSS_POINTS] = {
  {0.0},
  {-0.577350269189626, 0.577350269189626},
  {-0.774596669241483, 0.0, 0.774596669241483},
  {-0.861136311594053, -0.339981043584856, 0.339981043584856, 
   0.861136311594053},
  {-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 
   0.906179845938664},
  {-0.932469514203152, -0.661209386466265, -0.238619186083197, 
   0.238619186083197, 0.661209386466265, 0.932469514203152},
  {0.0},
  {0.0},
  {0.0},
  {0.0},
  {0.0},
  {-0.981560634246719, -0.904117256370475, -0.769902674194305,
   -0.587317954286617, -0.367831498998180, -0.125233408511469,
   0.125233408511469, 0.367831498998180, 0.587317954286617,
   0.769902674194305, 0.904117256370475, 0.981560634246719}
};

/* Weights (w_i) for Gaussian quadrature,
 *    \int_{-1}^1 f(x) \approx \sum_{i} w_i f(x_i)
 */
static double gauss_weight[MAX_GAUSS_POINTS][MAX_GAUSS_POINTS] = {
  {2.0},
  {1.0, 1.0},
  {0.555555555555556, 0.888888888888889, 0.555555555555556},
  {0.347854845137454, 0.652145154862546, 0.652145154862546, 
   0.347854845137454},
  {0.236926885056189, 0.478628670499366, 0.568888888888889, 
   0.478628670499366, 0.236926885056189},
  {0.171324492379170, 0.360761573048139, 0.467913934572691, 
   0.467913934572691, 0.360761573048139, 0.171324492379170},
  {0.0},
  {0.0},
  {0.0},
  {0.0},
  {0.0},
  {0.047175336386512, 0.106939325995318, 0.160078328543346,
   0.203167426723066, 0.233492536538355, 0.249147045813403,
   0.249147045813403, 0.233492536538355, 0.203167426723066, 
   0.160078328543346, 0.106939325995318, 0.047175336386512}
};

int MMH_ip = 0;


dbl vort_dir[MDE][DIM];		/* vorticity direction */
dbl qtensor[MDE][DIM][DIM];	/* I - 1/2 v^t v for each gauss point */
dbl div_qtensor[MDE][DIM];	/* div(I - v^t v) for each gauss point */


/* This routine populates div(V^tQV) = div(I-v^tv) for each gauss
 * point used for the particle species equation.
 */
void
assemble_qtensor(dbl *el_length) /* 2 x approximate element length scales */
{
  int print;
  int a, i, j, p;
  int ielem_type, ip_total, ip;
  dbl delta[DIM], xi[DIM], delta_xi[DIM], new_xi[DIM];
  dbl E[DIM][DIM];
  dbl vort_dir_local[DIM], vort_dir_delta[DIM], tmp;
  dbl d_vort_dir_d_x[DIM][DIM];
  dbl gd_local = 0.0;
  dbl gd_delta = 0.0;
  dbl *h;
  dbl *hq[DIM];
  dbl v1[DIM], v2[DIM], v3[DIM];

  memset(v1, 0, DIM * sizeof(double));
  memset(v2, 0, DIM * sizeof(double));
  memset(v3, 0, DIM * sizeof(double));

  print = 0;
  /*
  if(ei[pg->imtrx]->ielem == 15) print = 1;
  */

  for(i = 0; i < DIM; i++) {
     vort_dir_local[i] = 0.0;
     vort_dir_delta[i] = 0.0;
     delta_xi[i] = 0.0;
  }

/*   if(pd->e[pg->imtrx][R_MESH1]) EH(-1, "assemble_qtensor is not deformable mesh friendly."); */

  ielem_type = ei[pg->imtrx]->ielem_type;	/* element type */
  ip_total = elem_info(NQUAD, ielem_type); /* number of quadrature points */
  for(ip = 0; ip < ip_total; ip++)
    {
      /* First get the qtensor at the gauss point. */
      find_stu (ip, ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      load_basis_functions(xi, bfd);
      beer_belly();
      load_fv();
      load_bf_grad();
      load_fv_grads();
      h = fv->h;
      for(i = 0; i < DIM; i++)
	hq[i] = fv->hq[i];
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  E[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
      if(print)
	{
	  printf("BASE LOCATION:\n");
	  printf("    + % 10.4g % 10.4g % 10.4g +\n", E[0][0], E[0][1], E[0][2]);
	  printf("E = | % 10.4g % 10.4g % 10.4g |\n", E[1][0], E[1][1], E[1][2]);
	  printf("    + % 10.4g % 10.4g % 10.4g +\n", E[2][0], E[2][1], E[2][2]);
	}
      find_super_special_eigenvector(E, vort_dir_local, v1, v2, v3, &tmp, print);

      for(i = 0; i < DIM; i++) {
        memset( d_vort_dir_d_x[i], 0, DIM * sizeof(dbl) );
	for(j = 0; j < DIM; j++)
	  qtensor[ip][i][j] = (dbl)delta(i,j) -
	    0.5 * vort_dir_local[i] * vort_dir_local[j];
      }

      /*
      gd_local = 0.0;
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  gd_local += E[i][j] * E[j][i];
      gd_local = sqrt(0.5 * gd_local);
      */

      /* Compute d(vd)/dq(p) */
      for(p = 0; p < pd->Num_Dim; p++)
	{
	  delta[0] = delta[1] = delta[2] = 0.0;
	  delta[p] = FINITE_DELTA * (el_length[p]/2.0);
	  map_direction(delta_xi, delta);
	  for(i = 0; i < DIM; i++)
	    new_xi[i] = xi[i] + delta_xi[i];
	  load_basis_functions(new_xi, bfd);
	  beer_belly();
	  load_fv();
	  load_bf_grad();
	  load_fv_grads();
	  for(i = 0; i < DIM; i++)
	    for(j = 0; j < DIM; j++)
	      E[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];

          find_super_special_eigenvector(E, vort_dir_delta, v1, v2, v3, &tmp, print);

	  /*
	  gd_delta = 0.0;
	  for(i = 0; i < DIM; i++)
	    for(j = 0; j < DIM; j++)
	      gd_delta += E[i][j] * E[j][i];
	  gd_delta = sqrt(0.5 * gd_delta);
	  */

	  if(!bias_eigenvector_to(vort_dir_delta, vort_dir_local) || print)
	    {
	      printf("DELTA %d LOCATION:\n", p);
	      printf("             x = [% 10.4g, % 10.4g % 10.4g]\n",
	             fv->x[0], fv->x[1], fv->x[2]);
	      printf("                 / % 10.4g % 10.4g % 10.4g \\\n",
				             E[0][0], E[0][1], E[0][2]);
	      printf("             E = | % 10.4g % 10.4g % 10.4g |\n",
				             E[1][0], E[1][1], E[1][2]);
	      printf("                 \\ % 10.4g % 10.4g % 10.4g /\n",
				             E[2][0], E[2][1], E[2][2]);
	      printf("      vort_dir = [% 10.4g % 10.4g % 10.4g]\n",
		     vort_dir_local[0], vort_dir_local[1], vort_dir_local[2]);
	      printf("vort_dir_delta = [% 10.4g % 10.4g % 10.4g]\n",
		     vort_dir_delta[0], vort_dir_delta[1], vort_dir_delta[2]);
	      printf("|v| = % 10.4g\n", sqrt(fv->v[0]*fv->v[0] +
					     fv->v[1]*fv->v[1] +
					     fv->v[2]*fv->v[2]));
	      printf("sh = % 10.4g, gd = % 10.4g, gd_delta = % 10.4g\n",
		     fv->SH, gd_local, gd_delta);
	    }

	  for(i = 0; i < DIM; i++)
	        d_vort_dir_d_x[i][p] = (vort_dir_delta[p] - vort_dir_local[p])
			             / (delta[p]);
	  /*
	  d_gd_d_x[p] = (gd_delta - gd_local) / delta[p];
	  */
	}

      /* Now comptue div_qtensor[a] for this gauss point */
      div_qtensor[ip][0] = div_qtensor[ip][1] = div_qtensor[ip][2] = 0.0;
      for(a = 0; a < DIM; a++)
	{
	  for(i = 0; i < DIM; i++)
	    {
	      div_qtensor[ip][a] +=
		(hq[a][i] * vort_dir_local[a] * vort_dir_local[i] / h[a]
		 - hq[i][a] * vort_dir_local[i] * vort_dir_local[i] / h[a]
		 - hq[i][i] * vort_dir_local[i] * vort_dir_local[a] / h[i]) 
		/ h[i];
	      for(j = 0; j < DIM; j++)
		div_qtensor[ip][a] += hq[i][j] * vort_dir_local[j] * vort_dir_local[a] / h[i] / h[j];
	      div_qtensor[ip][a] += vort_dir_local[i] * d_vort_dir_d_x[a][i] / h[i];
	      div_qtensor[ip][a] += vort_dir_local[a] * d_vort_dir_d_x[i][i] / h[i];
	    }
	  div_qtensor[ip][a] *= -0.5;
	}
      
      /*
      for(a = 0; a < DIM; a++)
	grad_gd[ip][a] = d_gd_d_x[a] / h[a];

      for(i = 0; i < DIM; i++)
	vort_dir[ip][i] = vort_dir_local[i];
      gd[ip] = gd_local;
      */

      if(print)
	{
	  printf("    d/dq[0] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][0], d_vort_dir_d_x[1][0], d_vort_dir_d_x[2][0]);
	  printf("    d/dq[1] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][1], d_vort_dir_d_x[1][1], d_vort_dir_d_x[2][1]);
	  printf("    d/dq[2] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][2], d_vort_dir_d_x[1][2], d_vort_dir_d_x[2][2]);
	  printf("div_qtensor = [% 10.4g % 10.4g % 10.4g]\n",
		 div_qtensor[ip][0], div_qtensor[ip][1], div_qtensor[ip][2]);
	  /*
	  printf("    grad_gd = [% 10.4g % 10.4g % 10.4g]\n",
		 grad_gd[ip][0], grad_gd[ip][1], grad_gd[ip][2]);
	  printf("    grad_SH = [% 10.4g % 10.4g % 10.4g]\n",
		 fv->grad_SH[0], fv->grad_SH[1], fv->grad_SH[2]);
	  */
	  printf( "\n");
	}
    }

  return;
}

#ifdef NOT_USED
/* 
 *This routine calculates the Qtensor and its divergence using the vorticity
 * direction variable.
 */
void
assemble_qtensor_vort(dbl *el_length) /* 2 x approximate element length scales */
{
  int print;
  int a, i, j, p;
  int ielem_type, ip_total, ip;
  dbl delta[DIM], xi[DIM], delta_xi[DIM], new_xi[DIM];
  dbl E[DIM][DIM];
  dbl vort_dir_local[DIM], vort_dir_delta[DIM], tmp;
  dbl d_vort_dir_d_x[DIM][DIM];
  dbl gd_local = 0.0;
  dbl gd_delta = 0.0;
  dbl *h;
  dbl *hq[DIM];
  dbl v1[DIM], v2[DIM], v3[DIM];

  memset(v1, 0, DIM * sizeof(double));
  memset(v2, 0, DIM * sizeof(double));
  memset(v3, 0, DIM * sizeof(double));

  print = 0;

  for(i = 0; i < DIM; i++) {
     vort_dir_local[i] = 0.0;
     vort_dir_delta[i] = 0.0;
     delta_xi[i] = 0.0;
  }


  ielem_type = ei[pg->imtrx]->ielem_type;	/* element type */
  ip_total = elem_info(NQUAD, ielem_type); /* number of quadrature points */
  for(ip = 0; ip < ip_total; ip++)
    {
      /* First get the qtensor at the gauss point. */
      find_stu (ip, ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      load_basis_functions(xi, bfd);
      beer_belly();
      load_fv();
      load_bf_grad();
      load_fv_grads();
      h = fv->h;
      for(i = 0; i < DIM; i++)
	hq[i] = fv->hq[i];
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  E[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
      if(print)
	{
	  printf("BASE LOCATION:\n");
	  printf("    + % 10.4g % 10.4g % 10.4g +\n", E[0][0], E[0][1], E[0][2]);
	  printf("E = | % 10.4g % 10.4g % 10.4g |\n", E[1][0], E[1][1], E[1][2]);
	  printf("    + % 10.4g % 10.4g % 10.4g +\n", E[2][0], E[2][1], E[2][2]);
	}
      find_super_special_eigenvector(E, vort_dir_local, v1, v2, v3, &tmp, print);

      for(i = 0; i < DIM; i++) {
        memset( d_vort_dir_d_x[i], 0, DIM * sizeof(dbl) );
	for(j = 0; j < DIM; j++)
	  qtensor[ip][i][j] = (dbl)delta(i,j) -
	    0.5 * vort_dir_local[i] * vort_dir_local[j];
      }

      /*
      gd_local = 0.0;
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  gd_local += E[i][j] * E[j][i];
      gd_local = sqrt(0.5 * gd_local);
      */

      /* Compute d(vd)/dq(p) */
      for(p = 0; p < pd->Num_Dim; p++)
	{
	  delta[0] = delta[1] = delta[2] = 0.0;
	  delta[p] = FINITE_DELTA * (el_length[p]/2.0);
	  map_direction(delta_xi, delta);
	  for(i = 0; i < DIM; i++)
	    new_xi[i] = xi[i] + delta_xi[i];
	  load_basis_functions(new_xi, bfd);
	  beer_belly();
	  load_fv();
	  load_bf_grad();
	  load_fv_grads();
	  for(i = 0; i < DIM; i++)
	    for(j = 0; j < DIM; j++)
	      E[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];

          find_super_special_eigenvector(E, vort_dir_delta, v1, v2, v3, &tmp, print);

	  /*
	  gd_delta = 0.0;
	  for(i = 0; i < DIM; i++)
	    for(j = 0; j < DIM; j++)
	      gd_delta += E[i][j] * E[j][i];
	  gd_delta = sqrt(0.5 * gd_delta);
	  */

	  if(!bias_eigenvector_to(vort_dir_delta, vort_dir_local) || print)
	    {
	      printf("DELTA %d LOCATION:\n", p);
	      printf("             x = [% 10.4g, % 10.4g % 10.4g]\n",
	             fv->x[0], fv->x[1], fv->x[2]);
	      printf("                 / % 10.4g % 10.4g % 10.4g \\\n",
				             E[0][0], E[0][1], E[0][2]);
	      printf("             E = | % 10.4g % 10.4g % 10.4g |\n",
				             E[1][0], E[1][1], E[1][2]);
	      printf("                 \\ % 10.4g % 10.4g % 10.4g /\n",
				             E[2][0], E[2][1], E[2][2]);
	      printf("      vort_dir = [% 10.4g % 10.4g % 10.4g]\n",
		     vort_dir_local[0], vort_dir_local[1], vort_dir_local[2]);
	      printf("vort_dir_delta = [% 10.4g % 10.4g % 10.4g]\n",
		     vort_dir_delta[0], vort_dir_delta[1], vort_dir_delta[2]);
	      printf("|v| = % 10.4g\n", sqrt(fv->v[0]*fv->v[0] +
					     fv->v[1]*fv->v[1] +
					     fv->v[2]*fv->v[2]));
	      printf("sh = % 10.4g, gd = % 10.4g, gd_delta = % 10.4g\n",
		     fv->SH, gd_local, gd_delta);
	    }

	  for(i = 0; i < DIM; i++)
	        d_vort_dir_d_x[i][p] = (vort_dir_delta[p] - vort_dir_local[p])
			             / (delta[p]);
	  /*
	  d_gd_d_x[p] = (gd_delta - gd_local) / delta[p];
	  */
	}

      /* Now comptue div_qtensor[a] for this gauss point */
      div_qtensor[ip][0] = div_qtensor[ip][1] = div_qtensor[ip][2] = 0.0;
      for(a = 0; a < DIM; a++)
	{
	  for(i = 0; i < DIM; i++)
	    {
	      div_qtensor[ip][a] +=
		(hq[a][i] * vort_dir_local[a] * vort_dir_local[i] / h[a]
		 - hq[i][a] * vort_dir_local[i] * vort_dir_local[i] / h[a]
		 - hq[i][i] * vort_dir_local[i] * vort_dir_local[a] / h[i]) 
		/ h[i];
	      for(j = 0; j < DIM; j++)
		div_qtensor[ip][a] += hq[i][j] * vort_dir_local[j] * vort_dir_local[a] / h[i] / h[j];
	      div_qtensor[ip][a] += vort_dir_local[i] * d_vort_dir_d_x[a][i] / h[i];
	      div_qtensor[ip][a] += vort_dir_local[a] * d_vort_dir_d_x[i][i] / h[i];
	    }
	  div_qtensor[ip][a] *= -0.5;
	}
      
      /*
      for(a = 0; a < DIM; a++)
	grad_gd[ip][a] = d_gd_d_x[a] / h[a];

      for(i = 0; i < DIM; i++)
	vort_dir[ip][i] = vort_dir_local[i];
      gd[ip] = gd_local;
      */

      if(print)
	{
	  printf("    d/dq[0] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][0], d_vort_dir_d_x[1][0], d_vort_dir_d_x[2][0]);
	  printf("    d/dq[1] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][1], d_vort_dir_d_x[1][1], d_vort_dir_d_x[2][1]);
	  printf("    d/dq[2] = [% 10.4g % 10.4g % 10.4g]\n",
		  d_vort_dir_d_x[0][2], d_vort_dir_d_x[1][2], d_vort_dir_d_x[2][2]);
	  printf("div_qtensor = [% 10.4g % 10.4g % 10.4g]\n",
		 div_qtensor[ip][0], div_qtensor[ip][1], div_qtensor[ip][2]);
	  /*
	  printf("    grad_gd = [% 10.4g % 10.4g % 10.4g]\n",
		 grad_gd[ip][0], grad_gd[ip][1], grad_gd[ip][2]);
	  printf("    grad_SH = [% 10.4g % 10.4g % 10.4g]\n",
		 fv->grad_SH[0], fv->grad_SH[1], fv->grad_SH[2]);
	  */
	  printf( "\n");
	}
    }

  return;
}
#endif

/*
 * This is Ryan's qtensor assembly.
 */
void
assemble_new_qtensor(dbl *el_length) /* 2 x approximate element length scales */
{
  int err=0;
  int i, j, p, k1;
  int ielem_type, ip_total, ip;
  dbl dp;
  dbl delta[DIM], xi[DIM], delta_xi[DIM];
  dbl local_q[DIM][DIM], local_q2[DIM][DIM];
  dbl d_qtensor[DIM][DIM][DIM];


  /* Catch cases which are not available */
  if(pd->e[pg->imtrx][R_MESH1]) EH(-1, "assemble_qtensor is not deformable mesh friendly.");
  if (pd->Num_Dim == 3) EH(-1, "qtensor not ready for 3D problems yet!");
  if (pd->CoordinateSystem != CARTESIAN
   && pd->CoordinateSystem != PROJECTED_CARTESIAN
   && pd->CoordinateSystem != CARTESIAN_2pt5D)
    EH(-1, "Qtensor requires CARTESIAN coordinates for now!");

  /* Do some initializations */
  memset(qtensor, 0, MDE*DIM*DIM*sizeof(double));
  memset(div_qtensor, 0, MDE*DIM*sizeof(double));
  ielem_type = ei[pg->imtrx]->ielem_type;	/* element type */
  ip_total = elem_info(NQUAD, ielem_type); /* number of quadrature points */

  /* Loop over Gauss points */
  for (ip = 0; ip < ip_total; ip++)
    {
      /* Perform usual assembly sequence as if at a Gauss point */
      find_stu(ip, ielem_type, &xi[0], &xi[1], &xi[2]);
      err = load_basis_functions(xi, bfd);
      EH(err, "Problem from load_basis_functions!");
      err = beer_belly();
      EH(err, "Problem from beer_belly!");
      err = load_fv();
      EH(err, "Problem from load_fv!");
      err = load_bf_grad();
      EH(err, "Problem from load_bf_grad!");
      err = load_fv_grads();
      EH(err, "Problem from load_fv_grads!");

      /* Calculate Q-tensor at this Gauss point */
      err = get_local_qtensor(local_q);
      EH(err, "Problem getting local qtensor!");

      /* Fill these values into local array */
      for (i = 0; i < DIM; i++)
	{
          for (j = 0; j < DIM; j++)
	    {
	      qtensor[ip][i][j] = local_q[i][j];
	    }
	}

      /* Perturb the Gauss point location in each direction */
      memset(d_qtensor, 0, DIM*DIM*DIM*sizeof(double));
      for (p = 0; p < pd->Num_Dim; p++)
        {
          memset(delta, 0, DIM*sizeof(double));
          dp = fabs(FINITE_DELTA * el_length[p] / 2.0);
          if (dp < 0.001) dp = 0.001;
          for (k1 = 0; k1 < 2; k1++)
            {
              delta[p] = dp * ( (k1 == 0) ? 1.0 : -1.0);
              map_direction(delta_xi, delta);

              /* Perform assembly sequence up to load_fv_grads at this point */
              err = load_basis_functions(xi, bfd);
              EH(err, "Problem from load_basis_functions!");
              err = beer_belly();
              EH(err, "Problem from beer_belly!");
              err = load_fv();
              EH(err, "Problem from load_fv!");
              err = load_bf_grad();
              EH(err, "Problem from load_bf_grad!");
              err = load_fv_grads();
              EH(err, "Problem from load_fv_grads!");

              /* Evaluate qtensor at perturbed location */
              if (k1 == 0)
                {
                  err = get_local_qtensor(local_q);
                  EH(err, "Problem getting local qtensor!");
                }
              else
                {
                  err = get_local_qtensor(local_q2);
                  EH(err, "Problem getting local qtensor!");
                }
            }

          /* Get centered-difference approximations for d_qtensor */
          for (i = 0; i < DIM; i++)
            {
              for (j = 0; j < DIM; j++)
                {
                  d_qtensor[i][j][p] = (local_q[i][j] - local_q2[i][j])
                                     / (2.0 * dp);

                  /* Also sum up terms for div_qtensor */
                  if (i == p)
                    {
                      div_qtensor[ip][j] += d_qtensor[i][j][p];
                    }
                }
            }

        }  /* End of loop over directions (p) */

    }  /* End of loop over Gauss points (ip) */

  return;
}

static int
get_local_qtensor(double q[DIM][DIM])
/*
 * Calculate local q-tensor at one point using Ryan Miller's model. 
 * All quantities are loaded up in the usual places when called.
 * Currently limited to 2D fixed-grid problems in Cartesian coordinates.
 */
{
  int err=0, print=0;
  int i, j;
  dbl E[DIM][DIM], Erot[DIM][DIM];
  dbl mag_E, mag_Erot, rho_k;
  dbl theta, cc, ss;
  dbl stt, stc, scc;
  dbl ev[DIM], v_comp[DIM], v_vort[DIM], v_tens[DIM];


  /* Evaluate strain tensors E and Erot */
  for(i = 0; i < DIM; i++)
    {
      for(j = 0; j < DIM; j++)
        {
	  E[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
	  Erot[i][j] = fv->grad_v[i][j] - fv->grad_v[j][i];
        }
    }

  /* Get invariants of E and Erot */
  mag_E = 0.0;
  mag_Erot = 0.0;
  for(i = 0; i < DIM; i++)
    {
      for(j = 0; j < DIM; j++)
        {
          mag_E += E[i][j] * E[j][i];
          mag_Erot += Erot[i][j] * Erot[i][j];
        }
    }
  mag_E = sqrt(mag_E);
  mag_Erot = sqrt(mag_Erot);

  /* Set parameter rho_k */
  if (mag_Erot > mag_E)
    {
      rho_k = 1.0;
    }
  else if (mag_E == 0.0)
    {
      rho_k = 1.0;
    }
  else
    {
      rho_k = mag_Erot / mag_E;
    }

  /* Debugging lines */
  if (print)
    {
      printf("BASE LOCATION:\n");
      printf("    + % 10.4g % 10.4g % 10.4g +\n", E[0][0], E[0][1], E[0][2]);
      printf("E = | % 10.4g % 10.4g % 10.4g |\n", E[1][0], E[1][1], E[1][2]);
      printf("    + % 10.4g % 10.4g % 10.4g +\n", E[2][0], E[2][1], E[2][2]);
    }

  /* Solve for the three principal flow directions */
  diagonalize_symmetric_tensor(E, v_comp, v_vort, v_tens, ev, print);

  /* Get angle theta from tension eigenvector */
  if (fabs(v_tens[0]) > 1.0e-20)
    {
      theta = atan(v_tens[1] / v_tens[0]);
    }
  else
    {
      theta = 0.5 * M_PIE;
    }

  /* Sine and cosine of theta */
  cc = cos(-theta);
  ss = sin(-theta);

  /*
   * Initialize tensor components in rotated directions
   * defined by compression and tension axes
   * Note: assuming symmetric tensor!
   */
  scc = 0.5 * (rho_k + (1.0 - rho_k) * mp->Qtensor_Extension_P);
  stc = -mp->Qtensor_Nct;
  stt = scc;

  /* Transform qtensor into Cartesian coordinates (2D case only) */
  q[0][0] = stt* cc * cc + scc * ss * ss + 2.0 * stc * cc * ss;
  q[0][1] = (scc - stt) * cc * ss + stc * (cc * cc - ss * ss);
  q[1][0] = q[0][1];
  q[1][1] = stt * ss * ss + scc * cc * cc - 2.0 * stc * ss * cc;

  /* Done! */
  return err;
}

/* This routine computes the residual and jacobian values for the
 * vorticity principle flow direction equations.  BOTH the VORT_DIR{k}
 * and the VORT_LAMBDA equations are included here.  There is never a
 * case where only one or the other is being solved, so both are
 * always solved together.
 */
int
assemble_vorticity_direction(void)
{
  int a, b, i, j, k, p, q, r;		/* iteration counters */
  int eqn, var, peqn, pvar;	/* equation/variable indices */
  int print=0;	        	/* debugging print toggle. */


  dbl gamma_dot[DIM][DIM];
  dbl gamma_dot_pert[DIM][DIM];
  dbl det_J,			/* mapping jacobian */
    wt,				/* gauss point weight */
    h3;				/* volume scale */
  dbl	 phi_j;			/* bf[var]->phi[j] */
  int status=1;
  dbl advection, source;
  dbl vort_dir_local[DIM];
  dbl vort_dir[DIM];
  dbl vort_dir_pert[DIM];
  dbl tmp, dV;
  dbl wt_func;
  dbl dim;
  dbl R_old[DIM][MDE], R_new[DIM][MDE], eps;              /* For numerical Jacobians */
  dbl vv[DIM][MDE], grad_v[DIM][DIM];
  dbl dd[DIM][MDE], x[DIM], J[DIM][DIM], B[DIM][DIM] ;
  dbl grad_phi[MDE][DIM], grad_phi_e[MDE][DIM][DIM][DIM], d_phi[MDE][DIM];
  int dofs, index;
  int node, WIM;
  dbl xcoor[DIM][MDE], hh3, h[DIM], hq[DIM][DIM], radius, theta;
  dbl grad_e[DIM][DIM][DIM];
  dbl detJ1, det_J_inv;
  dbl v1[DIM], v2[DIM], v3[DIM];

  memset(v1, 0, DIM * sizeof(double));
  memset(v2, 0, DIM * sizeof(double));
  memset(v3, 0, DIM * sizeof(double));
  
  if ( ! pd->e[pg->imtrx][eqn = R_VORT_DIR1] )
    {
      return(status);
    }

  for (a=0; a < DIM; a++)
    {
      for (b=0; b < DIM; b++)
	{
	  gamma_dot[a][b] = 0.;
	  gamma_dot_pert[a][b] = 0.;
	}
    }
  
  dim = pd->Num_Dim;
  WIM  = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN ||
      pd->CoordinateSystem == CARTESIAN_2pt5D)
    WIM = WIM+1;

  wt = fv->wt;                 /* Numerical integration weight */
  
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  
  h3 = fv->h3;			/* Differential volume element (scales). */  
  
  dV = wt * det_J * h3;

  /* Super-duper special case for these two coordinate systems.  The
   * vorticity principle flow direction turns out to be in the
   * out-of-plane direction.  So, we are really "solving" the
   * following system:
   *
   *   v0 = 0
   *   v1 = 0
   *   v2 = 1
   *
   * The only other possibility is a regime where there is no shear
   * anywhere, and the rate of deformation tensor is identically 0.
   * In that case, it shouldn't matter what the vorticity principle
   * flow direction is... */


  /*
   * Residuals_________________________________________________________________
   */


  if(af->Assemble_Residual)
    {
      /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */
      
      for( a=0; a<VIM; a++ )
	{
	  for( b=0; b<VIM; b++)
	    {
	      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
	    }
	}
      
      /* Get local vorticity direction for least squares projection*/
      for ( a=0; a<DIM; a++)
	{
	  vort_dir_local[a] = 0.0;
	  vort_dir[a] = 0.;
	}
      vort_dir[0] = 1.;

      find_super_special_eigenvector(gamma_dot, vort_dir_local, v1, v2, v3, &tmp, print); 
      
      bias_eigenvector_to(vort_dir_local, vort_dir);
      memset(R_old, 0, sizeof(double)*DIM*MDE);
      
      for(a = 0; a < DIM; a++)
	{
	  eqn = R_VORT_DIR1 + a;
	  peqn = upd->ep[pg->imtrx][eqn];
	  for(i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
	    {
	      wt_func = bf[eqn]->phi[i];  

	      advection = 0.;
	      
	      if ( pd->e[pg->imtrx][eqn] & T_ADVECTION )
		{
		  advection = -vort_dir_local[a];
		  advection *= wt_func * dV;
		  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
		}
	      
	      /*
	       * Source term...
	       */
	      
	      source = 0;
	      
	      if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		{
		  source = fv->vd[a];    
		  source *= wt_func * dV;
		  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		}
	  
	      lec->R[peqn][i] += 
		advection  + source ;  
	      /* Save a convenient local copy for numerical Jacobians */
	      R_old[a][i] += advection  + source ;
	    }
	}
      
    } /* end if(af->Assemble_Residual) */

  if(af->Assemble_Jacobian)
    {
      /* J_vd_vd */
      for(a = 0; a < DIM; a++)
	{
	  eqn = R_VORT_DIR1 + a;
	  peqn = upd->ep[pg->imtrx][eqn];
	  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
	    {
	      wt_func = bf[eqn]->phi[i];
	      var = VORT_DIR1 + a;
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
	
		  source = 0.0;
		  if ( pd->e[pg->imtrx][eqn] & T_SOURCE )
		    {
		      source += phi_j;
		      source *= wt_func*dV;
		      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		    }		      
		  lec->J[peqn][pvar][i][j] +=source;
		}
	    }
	}

      /* J_vd_v */
      /* Must calculate numerical Jacobians since analytical are intractable! */
      memset(vv, 0, sizeof(double)*DIM*MDE);
      var = VELOCITY1;
      for(a = 0; a < VIM; a++)
	{
	  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      vv[a][j]=*esp->v[a][j];
	    }
	}
      for(b = 0; b < VIM; b++)
	{
	  var = VELOCITY1 + b;
	  pvar = upd->vp[pg->imtrx][var];
	  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	    {
	      phi_j = bf[var]->phi[j];
	      
	      if(*esp->v[b][j] != 0.)
		{
		  eps = *esp->v[b][j] * DELTA_UNKNOWN ;
		}
	      else
		{
		  eps = DELTA_UNKNOWN ;
		}
	      
	      vv[b][j]-=eps;
	      memset(grad_v, 0, sizeof(double)*DIM*DIM);
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      for ( r=0; r<VIM; r++)
			{
			  for ( k=0; k<ei[pg->imtrx]->dof[var]; k++)
			    {
			      grad_v[p][q] +=
				vv[r][k] * bf[var]->grad_phi_e[k][r] [p][q];
			    }
			}
		    }
		}
	      memset(gamma_dot_pert, 0, sizeof(double)*DIM*DIM);
	      for( p=0; p<VIM; p++ )
		{
		  for ( q=0; q<VIM; q++)
		    {
		      gamma_dot_pert[p][q] = grad_v[p][q] + grad_v[q][p];
		      
		    }
		}
	      
	      /* Get local vorticity direction for least squares projection*/
	      for ( p=0; p<VIM; p++ )
		vort_dir_pert[p] = vort_dir_local[p];
	      find_super_special_eigenvector(gamma_dot_pert, vort_dir_pert, v1, v2, v3, &tmp, print);
	      bias_eigenvector_to(vort_dir_pert, vort_dir);
	      
	      memset(R_new, 0, sizeof(double)*DIM*MDE);
	      for(a = 0; a < DIM; a++)
		{
		  eqn = R_VORT_DIR1 + a;
		  peqn = upd->ep[pg->imtrx][eqn];
		  
		  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		    {
		      wt_func = bf[eqn]->phi[i];
		      R_new[a][i] = (-vort_dir_pert[a]*pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]+
				     fv->vd[a]*pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)])*wt_func * dV;
		      lec->J[peqn][pvar][i][j] += (R_old[a][i]-R_new[a][i])/eps;
		      
		    }
		}
	      vv[b][j]=*esp->v[b][j];
	    }

	}


      /* J_vd_dm */
      /* Must calculate numerical Jacobians since analytical are intractable! */
      var = MESH_DISPLACEMENT1;
      
      if (pd->v[pg->imtrx][var]&& 0) {
	dofs = ei[pg->imtrx]->dof[var];
	memset(dd, 0, sizeof(double)*DIM*MDE);
	for( j=0; j<dofs; j++)
	  {
	    node  = ei[pg->imtrx]->dof_list[R_MESH1][j];
	    index = Proc_Elem_Connect[ Proc_Connect_Ptr[ei[pg->imtrx]->ielem] + node ];
	    for(a = 0; a < dim; a++)
	      {
		dd[a][j]=*esp->d[a][j];
		xcoor[a][j] = Coor[a][index];
	      }
	  }
	for(b = 0; b < dim; b++)
	  {
	    var = MESH_DISPLACEMENT1 + b;
	    pvar = upd->vp[pg->imtrx][var];
	    for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
	      {
		if(*esp->d[b][j] != 0.)
		  {
		    eps = *esp->d[b][j] * DELTA_UNKNOWN ;
		  }
		else
		  {
		    eps = DELTA_UNKNOWN ;
		  }
		dd[b][j] -= eps;
		for ( a=0; a< dim; a++)
		  {
		    x[a] = 0.;
		    for ( k=0; k<dofs; k++)
		      {
			x[a]     +=  ( xcoor[a][k] + dd[a][k]) * bf[R_MESH1]->phi[k];
		      }
		  }
		memset(J, 0, sizeof(double)*DIM*DIM);
		for ( a=0; a<dim; a++)
		  {
		    for ( p=0; p<dim; p++)
		      {
			for ( k=0; k<dofs; k++)
			  {
			    J[a][p] 
			      += ( xcoor[p][k] + dd[p][k]) * bf[R_MESH1]->dphidxi[k][a];
			  }
		      }
		  }
		if(dim == 2 )
		  {
		    detJ1    =  J[0][0] * J[1][1] - J[0][1] * J[1][0];
		    det_J_inv = 1/detJ1;
		    
		    B[0][0] =  J[1][1] * det_J_inv;
		    B[0][1] = -J[0][1] * det_J_inv;
		    B[1][0] = -J[1][0] * det_J_inv;
		    B[1][1] =  J[0][0] * det_J_inv;
		  }    
		else if(dim == 3 )
		  {
		    detJ1 = J[0][0] * ( J[1][1] * J[2][2] - J[1][2] * J[2][1])
	                  - J[0][1] * ( J[1][0] * J[2][2] - J[2][0] * J[1][2])
		          + J[0][2] * ( J[1][0] * J[2][1] - J[2][0] * J[1][1]);
		    det_J_inv = 1/detJ1;
		    
		    B[0][0] = (  J[1][1] *  J[2][2] - J[2][1] *  J[1][2])* det_J_inv;
		    B[0][1] =-(  J[0][1] *  J[2][2] - J[2][1] *  J[0][2])* det_J_inv;
		    B[0][2] = (  J[0][1] *  J[1][2] - J[1][1] *  J[0][2])* det_J_inv;
		    B[1][0] =-(  J[1][0] *  J[2][2] - J[2][0] *  J[1][2])* det_J_inv;
		    B[1][1] = (  J[0][0] *  J[2][2] - J[2][0] *  J[0][2])* det_J_inv;
		    B[1][2] =-(  J[0][0] *  J[1][2] - J[1][0] *  J[0][2])* det_J_inv;
		    B[2][0] = (  J[1][0] *  J[2][1] - J[1][1] *  J[2][0])* det_J_inv;
		    B[2][1] =-(  J[0][0] *  J[2][1] - J[2][0] *  J[0][1])* det_J_inv;
		    B[2][2] = (  J[0][0] *  J[1][1] - J[1][0] *  J[0][1])* det_J_inv;
		  }
		
		memset(d_phi,0,DIM*MDE*sizeof(double));
		for ( k=0; k<dofs; k++)
		  {
		    for ( p=0; p<dim; p++)
		      {
			for ( q=0; q<dim; q++)
			  {
			    d_phi[k][p] +=B[p][q]*bf[var]->dphidxi[k][q];
			  }
		      }
		  }
		
		for ( p=0; p<DIM; p++)
		  {
		    h[p] = 1.;
		  }
		memset(hq, 0,DIM*DIM*sizeof(double));
		if( pd->CoordinateSystem == CARTESIAN 
		    || pd->CoordinateSystem == PROJECTED_CARTESIAN 
		    || pd->CoordinateSystem == CARTESIAN_2pt5D )
		  { /* Do nothing. Everything stays 1. */
		  }
		else if ( pd->CoordinateSystem == CYLINDRICAL 
			  || pd->CoordinateSystem == SWIRLING )
		  {   
		    radius = x[1];
		    h[2]     = radius;
		    hq[2][1] = 1.;
		    if (radius == 0.)
		      {
			h[2] = 1.;
			hq[2][1] =  0.;
		      } 
		  }
		else if ( pd->CoordinateSystem == SPHERICAL )
		  {
		    radius     = x[0];
		    theta = x[1];

		    h[1]     = radius;
		    h[2]     = radius*sin(theta);

		    hq[1][0] = 1.;
		    hq[2][0] = sin(theta);
		    hq[2][1] = radius*cos(theta);
		  }

		hh3 = 1.;
		for ( i=0; i<VIM; i++)
		  {
		    hh3 *= h[i];
		  }
		memset(grad_phi,0,DIM*MDE*sizeof(double));
		for ( k=0; k<dofs; k++)
		  {
		    for ( p=0; p<VIM; p++ )
		      {
			grad_phi[k][p] = (d_phi[k][p])/h[p];
		      }
		  }
		memset(grad_phi_e,0,DIM*DIM*DIM*MDE*sizeof(double));
		
		for ( k=0; k<dofs; k++)
		  {
		    for ( p=0; p<VIM; p++)
		      {
			for ( q=0; q<VIM; q++)
			  {
			    grad_phi_e[k][q] [p][q] =  grad_phi[k][p];
			  }
		      }
		  }
		memset(grad_e, 0, sizeof(double)*DIM*DIM*DIM);
	
		if( pd->CoordinateSystem != CARTESIAN 
		    && pd->CoordinateSystem != PROJECTED_CARTESIAN 
		    && pd->CoordinateSystem != CARTESIAN_2pt5D )
		  {
		    for ( q=0; q<VIM; q++)
		      {
			for ( p=0; p<VIM; p++)
			  {
			    grad_e[q][p][p] 
			      = hq[p][q] / ( h[p] * h[q] );
			  }
		      }
		    for ( p=0; p<VIM; p++)
		      {
			for (q=0; q<VIM; q++)
			  {
			    grad_e[p][p][q] 
			      -= hq[p][q] / ( h[p] * h[q] );
			  }
		      }
		  }

		if( pd->CoordinateSystem != CARTESIAN 
		    && pd->CoordinateSystem != PROJECTED_CARTESIAN 
		    && pd->CoordinateSystem != CARTESIAN_2pt5D )
		  {
		    for ( k=0; k<dofs; k++)
		      {
			for ( r=0; r<WIM; r++)
			  {
			    for ( p=0; p<VIM; p++)
			      {
				for ( q=0; q<VIM; q++)
				  {
				    grad_phi_e[k][r] [p][q] += bf[var]->phi[k] * grad_e[r][p][q];
				  }
			      }
			  }
		      }
		    
		  }
		
		
		memset(grad_v, 0, sizeof(double)*DIM*DIM);
		for ( p=0; p<VIM; p++)
		  {
		    for ( q=0; q<VIM; q++)
		      {
			for ( r=0; r<VIM; r++)
			  {
			    for ( k=0; k<ei[pg->imtrx]->dof[var]; k++)
			      {
				grad_v[p][q] +=
				  *esp->v[r][k] * grad_phi_e[k][r] [p][q];
			      }
			  }
		      }
		  }
		memset(gamma_dot_pert, 0, sizeof(double)*DIM*DIM);
		for( p=0; p<VIM; p++ )
		  {
		    for ( q=0; q<VIM; q++)
		      {
			gamma_dot_pert[p][q] = grad_v[p][q] + grad_v[q][p];
			
		      }
		  }
		
		/* Get local vorticity direction for least squares projection*/
		for ( p=0; p<VIM; p++ )
		  vort_dir_pert[p] = vort_dir_local[p];
		find_super_special_eigenvector(gamma_dot_pert, vort_dir_pert, v1, v2, v3, &tmp, print);
		bias_eigenvector_to(vort_dir_pert, vort_dir);
		
		memset(R_new, 0, sizeof(double)*DIM*MDE);
		for(a = 0; a < DIM; a++)
		  {
		    eqn = R_VORT_DIR1 + a;
		    peqn = upd->ep[pg->imtrx][eqn];
		    
		    for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
		      {
			wt_func = bf[eqn]->phi[i];
			R_new[a][i] = (-vort_dir_pert[a]*pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)]+
				       fv->vd[a]*pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)])*wt_func 
			  * wt * hh3 * detJ1;
			lec->J[peqn][pvar][i][j] += (R_old[a][i]-R_new[a][i])/eps;
			
		      }
		  }
		dd[b][j]=*esp->d[b][j];
	      }
	  }
	
      }
      
      
      
    }
  
  return 1;
}

/* This routine calculates the Eigenvectors of the shear rate tensor given a velocity field */

/* EDW: not used as of 07/05/2006!
int
get_shearrate_eigenvector(dbl velocity[DIM][MDE], dbl vort_dir[DIM])
{
  return 0;
}
*/

/* This is the main function that is called from hydro_flux.  It
 * implements the Q-tensor diffusivity model.  Upon exit, the
 * st->diff_flux structure is filled.
 *
 * It is assumed that the VORT_DIR{k} and VORT_LAMBDA
 * equations/variables are included in the problem...  NOT ANYMORE
 *  
 * A lot of the structures here are named in terms of VQVt.  This is
 * the same as qtensor.
 *
 *  Author: M.M. Hopkins, 8/99 - 03/01
 */
int
hydro_qtensor_flux (struct Species_Conservation_Terms *st,
		    int w)                      /* species number */
{
  int a, b, i, j, var;
  int status=1;
  dbl gammadot, *grad_gammadot, gamma_dot[DIM][DIM];

/*   dbl qtensor_const[DIM][DIM]; */ /* for analytical qtensor hardwire */

  dbl del_rho;
  dbl mu, grad_mu[DIM];
  dbl grad_mu1[DIM];
  dbl Dc, Dmu, Dg;
  dbl *Y, (*grad_Y)[DIM], *dmu_dY, *d2mu_dY2;

  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl c_term, mu_term;

  dbl phi_j;
  dbl *grad_phi_j;

  dbl div_gdYVQVt[DIM]; /* div(gammadot * Y * (V Q Vt)) */
  dbl VQVt_grad_mu[DIM]; /* (V Q V^t) . grad(mu) */
  dbl grad_Y_VQVt[DIM]; /* grad(Y) . (V Q V^t) */
  dbl VQVt_grad_Y[DIM]; /* (V Q V^t) . grad(Y) */
  dbl grad_gd_VQVt[DIM]; /* grad(gammadot) . (V Q V^t) */
  dbl grad_phi_j_VQVt[DIM]; /* grad(phi_j) . (V Q V^t) */
  dbl VQVt_grad_phi_j[DIM]; /* (V Q V^t) . grad(phi_j) */

  if(!pd->e[pg->imtrx][VORT_DIR3] && 0)
    {
      EH(-1, "Cannot use QTENSOR without the VORT_DIR{1,2,3} equations/variables active!");
      exit(-1);
    }

  if(MMH_ip < 0)
    {
      printf("Entered with bad MMH_ip\n");
      MMH_ip = 1;
    }

  /* Set up some convenient local variables and pointers */
  Y = fv->c;
  grad_Y = fv->grad_c;
  grad_gammadot = fv->grad_SH;
  gammadot = fv->SH;
  /*
  grad_gammadot = grad_gd[MMH_ip];
  gammadot = gd[MMH_ip];
  */
  
  dmu_dY = &(mp->d_viscosity[MAX_VARIABLE_TYPES]);
  d2mu_dY2 = &(mp->d2_viscosity[MAX_VARIABLE_TYPES]);
  
  memset(gamma_dot, 0, DIM*DIM*sizeof(dbl) );

  for(a = 0; a < DIM; a++)
    for(b = 0; b < DIM; b++)
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
  
  /* get mu and grad_mu */
  mu = viscosity(gn, gamma_dot, d_mu);

  memset(grad_mu, 0, DIM*sizeof(dbl));
  memset(grad_mu1, 0, DIM*sizeof(dbl));

  /* Controversial as to whether the grad should
   * depend on all variables or just concentration.
   * For now we will just do concentration.
   *
   * Separate grad_mu into two parts, so that we can
   * use two different coefficients. The first part
   * gives the sensitivity wrt to concentration and
   * the second gives the sensitivity wrt to other
   * fields such as shear-rate etc.
   * This is to get a better match with the data.
   */
  for(a=0; a<DIM; a++)
    {
      grad_mu[a] += dmu_dY[w] * fv->grad_c[w][a];
      grad_mu1[a] += d_mu->gd * fv->grad_SH[a];
    }
  
  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */

  /* If non-neutrally bouyant suspension, compute density difference 
   *  or it defaults to zero 
   */
  del_rho = 0.;
  if( mp->DensityModel == SUSPENSION)
    {
      del_rho = (mp->u_density[2] -  mp->u_density[1]);
    }
  
  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */
  
  if(Y[w] > 0. )
    {
      if (mp->GamDiffType[w] == LINEAR) 
	{
	  Dc  = mp->u_gadiffusivity[w][0] * 1.4 * Y[w];
	}
      else if (mp->GamDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_gadiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gadiffusivity[w][2];

	  ls_transport_property (  mp->u_gadiffusivity[w][0],
				   mp->u_gadiffusivity[w][1],
				   width,
				   &Dc,
				   NULL );
	}

      else
	{
	  Dc = mp->gam_diffusivity[w];
	}
      
      if (mp->MuDiffType[w] == LINEAR) 
	{
	  Dmu  = mp->u_mdiffusivity[w][0] * 1.4 * Y[w];
	}
      else if (mp->MuDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_mdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_mdiffusivity[w][2];

	  ls_transport_property (  mp->u_mdiffusivity[w][0],
				   mp->u_mdiffusivity[w][1],
				   width,
				   &Dmu,
				   NULL );
	}
      else
	{
	  Dmu = mp->mu_diffusivity[w];
	}
      
      if (mp->GravDiffType[w] == BISECTION)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if(mp->GravDiffType[w] == RZBISECTION)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if(mp->GravDiffType[w] == RICHARDSON_ZAKI)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if (mp->GravDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_gdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gdiffusivity[w][2];

	  ls_transport_property (  mp->u_gdiffusivity[w][0],
				   mp->u_gdiffusivity[w][1],
				   width,
				   &Dg,
				   NULL );
	  Dg  *= del_rho;
	}
      else
	{
	  Dg  = mp->g_diffusivity[w]*del_rho;
	}
    }
  else
    {
      Dg  = 0.;
      Dmu = 0.;
      Dc = 0.;
    } 



  
  /* MMH: Set Q tensor components.  The Q tensor is really V Q_0 V^t,
   * where V = [v_flow v_norm v_vort].  V should be orthogonal (v_*
   * should be an orthonormal basis for R^3) so that V^t = V^-1.
   * Although V is heavily dependant on the fluid velocity, no
   * contribution to the Jacobian is calculated for them.  If it
   * doesn't converge, then I will go back and compute some Jacobian
   * entries (or maybe secant approximations).  Analytical Jacobian
   * entries would be very, very, very, very ugly.
   *
   * NOTE 1: The Q-tensor formulation is qualitatively different enough
   * from the flux-diffusion model to warrant a new construction (as
   * opposed to modifying the flux diffusion model).
   *
   * NOTE 2: I assume that Dc = Kc * a * a. (Only the CONSTANT model,
   * not the LINEAR model.)  Similarly for Dmu.
   */

  /* Compute grad(Y) . (V Q V^t)
   *         grad(gammadot) . (V Q V^t)
   *         (V Q V^t) . grad(mu)
   */
  memset(grad_Y_VQVt, 0, DIM*sizeof(dbl));
  memset(grad_gd_VQVt, 0, DIM*sizeof(dbl));
  memset(VQVt_grad_mu, 0, DIM*sizeof(dbl));
  memset(VQVt_grad_Y, 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    for(i = 0; i <DIM; i++)
      {
	grad_Y_VQVt[a] += grad_Y[w][i] * qtensor[MMH_ip][i][a];
	grad_gd_VQVt[a] += grad_gammadot[i] * qtensor[MMH_ip][i][a];
	VQVt_grad_mu[a] += qtensor[MMH_ip][a][i] * grad_mu[i];
	VQVt_grad_Y[a] += qtensor[MMH_ip][a][i] * grad_Y[w][i];
      }


  /* Compute div(gammadot * Y * (V Q V^t)) */
  memset(div_gdYVQVt, 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    {
      div_gdYVQVt[a] += Y[w] * grad_gd_VQVt[a];
      div_gdYVQVt[a] += gammadot * grad_Y_VQVt[a];
      div_gdYVQVt[a] += gammadot * Y[w] * div_qtensor[MMH_ip][a];
    }


  /* Assemble residual */
  memset(st->diff_flux[w], 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    {
      st->diff_flux[w][a] -= Dc * Y[w] * div_gdYVQVt[a];
      st->diff_flux[w][a] -= Dmu * Y[w] * Y[w] * gammadot * VQVt_grad_mu[a] / mu;
    }


  /* Assemble Jacobian */
  /* Currently no mesh displacement dependencies. */
  if(af->Assemble_Jacobian)
    {
      var = MASS_FRACTION;
      memset(st->d_diff_flux_dc, 0, MAX_CONC*DIM*MAX_CONC*MDE*sizeof(dbl));
      for(a = 0; a < DIM; a++)
	for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	  {
	    c_term = 0.0;
	    mu_term = 0.0;
	    phi_j = bf[var]->phi[j];
	    grad_phi_j = bf[var]->grad_phi[j];
		
	    /* Compute grad(phi) . VQVt.  This is different for
	     * different dependency functions (when
	     * bf[var]->grad_phi[j] changes).
	     */
	    memset(VQVt_grad_phi_j, 0, DIM * sizeof(dbl));
	    memset(grad_phi_j_VQVt, 0, DIM * sizeof(dbl));
	    for(b = 0; b < DIM; b++)
	      for(i = 0; i < DIM; i++)
		{
		  VQVt_grad_phi_j[b] += qtensor[MMH_ip][b][i] * grad_phi_j[i];
		  grad_phi_j_VQVt[b] += grad_phi_j[i] * qtensor[MMH_ip][i][b];
		}


	    c_term += phi_j * div_gdYVQVt[a];
	    c_term += Y[w] * gammadot * grad_phi_j_VQVt[a];
	    c_term += Y[w] * phi_j * (grad_gd_VQVt[a] + gammadot * div_qtensor[MMH_ip][a]);
	    c_term *= -Dc;

	    mu_term += 2.0 * phi_j * dmu_dY[w] * VQVt_grad_Y[a];
	    mu_term += (d2mu_dY2[w] - dmu_dY[w] * dmu_dY[w] / mu) * Y[w] * VQVt_grad_Y[a];
	    mu_term += Y[w] * dmu_dY[w] * VQVt_grad_phi_j[a];
	    mu_term *= -Dmu * gammadot * Y[w] / mu;

	    st->d_diff_flux_dc[w][a] [w][j] = c_term + mu_term;
#ifdef DEBUG_QTENSOR
	    if(print_info)
	      if(st->d_diff_flux_dc[w][a][w][j] == 0.0)
		printf( "d_diff_flux_dc[%d][%d][%d][%d] is zero!\n",
			w, a, w, j);
#endif
	  }

      var = SHEAR_RATE;
      memset(st->d_diff_flux_dSH[w], 0, DIM*MDE*sizeof(dbl));
      for(a = 0; a < DIM; a++)
	for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	  {
	    c_term = 0.0;
	    mu_term = 0.0;
	    phi_j = bf[var]->phi[j];
	    grad_phi_j = bf[var]->grad_phi[j];
		
	    /* Compute grad(phi) . VQVt.  This is different for
	     * different dependency functions (when
	     * bf[var]->grad_phi[j] changes).
	     */
	    memset(grad_phi_j_VQVt, 0, DIM * sizeof(dbl));
	    for(b = 0; b < DIM; b++)
	      for(i = 0; i < DIM; i++)
		grad_phi_j_VQVt[b] += grad_phi_j[i] * qtensor[MMH_ip][i][b];

	    c_term += Y[w] * grad_phi_j_VQVt[a];
	    c_term += phi_j * (grad_Y_VQVt[a] + Y[w] * div_qtensor[MMH_ip][a]);
	    c_term *= -Dc * Y[w];
		
	    mu_term = -Dmu * Y[w] * Y[w] * VQVt_grad_mu[a] / mu * phi_j;

	    st->d_diff_flux_dSH[w][a][j] = c_term + mu_term;
#ifdef DEBUG_QTENSOR
	    if(print_info)
	      if(st->d_diff_flux_dSH[w][a][j] == 0.0)
		printf( "d_diff_flux_dSH[%d][%d][%d] is zero!\n",
			w, a, j);
#endif
	  }
    }	      
  return(status);
}

/* This is the main function that is called from hydro_flux.  It
 * implements the Q-tensor diffusivity model.  Upon exit, the
 * st->diff_flux structure is filled.
 *
 * It is assumed that the VORT_DIR{k}
 * equations/variables are included in the problem...  
 *  
 * A lot of the structures here are named in terms of VQVt.  This is
 * the same as qtensor.
 *
 *  Author: M.M. Hopkins, 8/99 - 03/01
 *          R.R. Rao, 7/04
 */
int
hydro_qtensor_flux_new (struct Species_Conservation_Terms *st,
		    int w)                      /* species number */
{
  int a, b, i, j, var;
  int p, q, imtrx, vort_dir_on=0;
  int status=1;
  dbl gammadot, *grad_gammadot, gamma_dot[DIM][DIM];

/*   dbl qtensor_const[DIM][DIM]; */ /* for analytical qtensor hardwire */

  dbl del_rho;
  dbl mu, grad_mu[DIM];
  dbl grad_mu1[DIM];
  dbl Dc, Dmu, Dg;

  dbl Dd[DIM];
  dbl dDd_dy[DIM];
  dbl dDd_dgrady[DIM];

  dbl *Y, (*grad_Y)[DIM], *dmu_dY, *d2mu_dY2;

  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl c_term, mu_term, d_term;

  dbl phi_j = 0.0;

  dbl maxpack = 0.0;  /*  ERROR !!! maxpack is never set before use !!!! */

  dbl div_gdYVQVt[DIM]; /* div(gammadot * Y * (V Q Vt)) */
  dbl VQVt_grad_mu[DIM]; /* (V Q V^t) . grad(mu) */
  dbl grad_Y_VQVt[DIM]; /* grad(Y) . (V Q V^t) */
  dbl VQVt_grad_Y[DIM]; /* (V Q V^t) . grad(Y) */
  dbl grad_gd_VQVt[DIM]; /* grad(gammadot) . (V Q V^t) */
  dbl grad_phi_Q[DIM][MDE]; /* grad(phi_j) . (Q ) */
  dbl Q_grad_phi[DIM][MDE]; /* (Q) . grad(phi_j) */

  dbl qtensor_loc[DIM][DIM];
  dbl d_qtensor_dvd[DIM][DIM][DIM][MDE];
  
  dbl div_phi_j_e_b;
  dbl div_q[DIM];
  dbl d_div_q_dvd[DIM][DIM][MDE];
  dbl d_grad_Y_VQVt_dvd[DIM][DIM][MDE];
  dbl d_grad_gd_VQVt_dvd[DIM][DIM][MDE];
  dbl d_VQVt_grad_mu_dvd[DIM][DIM][MDE];
  dbl d_VQVt_grad_Y_dvd[DIM][DIM][MDE];

  dbl d_div_gdYVQVt_dvd[DIM][DIM][MDE];
  
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) 
    {
      if(pd->e[imtrx][VORT_DIR1])
	{
	  vort_dir_on = 1;
	}
    }


  if(!vort_dir_on)
    {
      EH(-1, "Cannot use this QTENSOR model without the VORT_DIR{1,2,3} equations/variables active!");
      exit(-1);
    }


  /* Set up some convenient local variables and pointers */
  Y = fv->c;
  grad_Y = fv->grad_c;
  grad_gammadot = fv->grad_SH;
  gammadot = fv->SH;
  
  dmu_dY = &(mp->d_viscosity[MAX_VARIABLE_TYPES]);
  d2mu_dY2 = &(mp->d2_viscosity[MAX_VARIABLE_TYPES]);
  
  memset(gamma_dot, 0, DIM*DIM*sizeof(dbl) );
  memset(dDd_dgrady, 0, DIM*sizeof(dbl));

  for(a = 0; a < DIM; a++)
    for(b = 0; b < DIM; b++)
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
  
  /* get mu and grad_mu */
  mu = viscosity(gn, gamma_dot, d_mu);

  memset(grad_mu, 0, DIM*sizeof(dbl));
  memset(grad_mu1, 0, DIM*sizeof(dbl));

  /* Controversial as to whether the grad should
   * depend on all variables or just concentration.
   * For now we will just do concentration.
   *
   * Separate grad_mu into two parts, so that we can
   * use two different coefficients. The first part
   * gives the sensitivity wrt to concentration and
   * the second gives the sensitivity wrt to other
   * fields such as shear-rate etc.
   * This is to get a better match with the data.
   */
  for(a=0; a<DIM; a++)
    {
      grad_mu[a] += dmu_dY[w] * fv->grad_c[w][a];
      grad_mu1[a] += d_mu->gd * fv->grad_SH[a];
    }
  
  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */

  /* If non-neutrally bouyant suspension, compute density difference 
   *  or it defaults to zero 
   */
  del_rho = 0.;
  if( mp->DensityModel == SUSPENSION)
    {
      del_rho = (mp->u_density[2] -  mp->u_density[1]);
    }
  
  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */
  
  if(Y[w] > 0. )
    {
      if (mp->GamDiffType[w] == LINEAR) 
	{
	  Dc  = mp->u_gadiffusivity[w][0] * 1.4 * Y[w];
	}
      else if (mp->GamDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_gadiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gadiffusivity[w][2];

	  ls_transport_property (  mp->u_gadiffusivity[w][0],
				   mp->u_gadiffusivity[w][1],
				   width,
				   &Dc,
				   NULL );
	}

      else
	{
	  Dc = mp->gam_diffusivity[w];
	}
      
      if (mp->MuDiffType[w] == LINEAR) 
	{
	  Dmu  = mp->u_mdiffusivity[w][0] * 1.4 * Y[w];
	}
      else if (mp->MuDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_mdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_mdiffusivity[w][2];

	  ls_transport_property (  mp->u_mdiffusivity[w][0],
				   mp->u_mdiffusivity[w][1],
				   width,
				   &Dmu,
				   NULL );
	}
      else
	{
	  Dmu = mp->mu_diffusivity[w];
	}
      
      if (mp->GravDiffType[w] == BISECTION)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if(mp->GravDiffType[w] == RZBISECTION)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if(mp->GravDiffType[w] == RICHARDSON_ZAKI)
	{
	  Dg  = mp->u_gdiffusivity[w][0]*del_rho;
	}
      else if (mp->GravDiffType[w] == LEVEL_SET ) 
	{
	  double width;
	  if ( ls == NULL ) EH(-1,"Need to activate to Level Set Interface Tracking to use this model.\n");

	  width = ( mp->u_gdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gdiffusivity[w][2];

	  ls_transport_property (  mp->u_gdiffusivity[w][0],
				   mp->u_gdiffusivity[w][1],
				   width,
				   &Dg,
				   NULL );
	  Dg  *= del_rho;
	}
      else
	{
	  Dg  = mp->g_diffusivity[w]*del_rho;
	}
    }
  else
    {
      Dg  = 0.;
      Dmu = 0.;
      Dc = 0.;
    } 

  
  /*  anisotropic diffusion coefficient set
   *  in mat file it is direction dependent
   */
  
  if(mp->FickDiffType[w] == ANISOTROPIC )
    {
      for ( a=0; a<DIM; a++)
	{
	  Dd[a]  = mp->u_fdiffusivity[w][a];
	  dDd_dy[a] = 0.;
	}
    }
  else if(mp->FickDiffType[w] == EXP_DECAY )
    {
      for ( a=0; a<DIM; a++)
	{
	  if(Y[w] >= 0. && Y[w] <= maxpack)
	    {
	      Dd[a]  = mp->u_fdiffusivity[w][0] * 
		(exp(-mp->u_fdiffusivity[w][1] * Y[w])
		 +exp(-mp->u_fdiffusivity[w][1] * fabs (maxpack -Y[w])));
	      dDd_dy[a] = -mp->u_fdiffusivity[w][0] * mp->u_fdiffusivity[w][1]*
		(exp(-mp->u_fdiffusivity[w][1] * Y[w])
		 -exp(-mp->u_fdiffusivity[w][1] * fabs (maxpack -Y[w])));
	    }
	  else if(Y[w] < 0.)
	    {
	      Dd[a]  = mp->u_fdiffusivity[w][0] ;
	      dDd_dy[a] = 0.;
	    }
	  else if(Y[w] > maxpack)
	    {
	      Dd[a]  = mp->u_fdiffusivity[w][0] ;
	      dDd_dy[a] = 0.;
	    }
	  
	}
    } 
  

  
  /* MMH: Set Q tensor components.  The Q tensor is really V Q_0 V^t,
   * where V = [v_flow v_norm v_vort].  V should be orthogonal (v_*
   * should be an orthonormal basis for R^3) so that V^t = V^-1.
   * Although V is heavily dependant on the fluid velocity, no
   * contribution to the Jacobian is calculated for them.  If it
   * doesn't converge, then I will go back and compute some Jacobian
   * entries (or maybe secant approximations).  Analytical Jacobian
   * entries would be very, very, very, very ugly.
   *
   * NOTE 1: The Q-tensor formulation is qualitatively different enough
   * from the flux-diffusion model to warrant a new construction (as
   * opposed to modifying the flux diffusion model).
   *
   * NOTE 2: I assume that Dc = Kc * a * a. (Only the CONSTANT model,
   * not the LINEAR model.)  Similarly for Dmu.
   */

  /* Compute grad(Y) . (V Q V^t)
   *         grad(gammadot) . (V Q V^t)
   *         (V Q V^t) . grad(mu)
   */


  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  qtensor_loc[a][b] = (dbl)delta(a,b) -
	    0.5 * fv->vd[a] * fv->vd[b];
	}
    }
 
  memset(div_q, 0, DIM*sizeof(dbl));
  for ( a=0; a<VIM; a++)
    {
      div_q[a] = -0.5*fv->vd[a]*fv->div_vd;
      for ( b=0; b<VIM; b++)
	{
	  div_q[a] -= 0.5*fv->vd[b] * fv->grad_vd[b][a];
	}
    }
 

  memset(grad_Y_VQVt, 0, DIM*sizeof(dbl));
  memset(grad_gd_VQVt, 0, DIM*sizeof(dbl));
  memset(VQVt_grad_mu, 0, DIM*sizeof(dbl));
  memset(VQVt_grad_Y, 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    for(i = 0; i <DIM; i++)
      {
	grad_Y_VQVt[a] += grad_Y[w][i] * qtensor_loc[i][a];
	grad_gd_VQVt[a] += grad_gammadot[i] * qtensor_loc[i][a];
	VQVt_grad_mu[a] += qtensor_loc[a][i] * grad_mu[i];
	VQVt_grad_Y[a] += qtensor_loc[a][i] * grad_Y[w][i];
      }



  /* Compute div(gammadot * Y * (V Q V^t)) */
  memset(div_gdYVQVt, 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    {
      div_gdYVQVt[a] = Y[w] * grad_gd_VQVt[a];
      div_gdYVQVt[a] += gammadot * grad_Y_VQVt[a];
      div_gdYVQVt[a] += gammadot * Y[w] * div_q[a];
    }


  /* Assemble residual */
  memset(st->diff_flux[w], 0, DIM*sizeof(dbl));
  for(a = 0; a < DIM; a++)
    {
      st->diff_flux[w][a] -= Dc * Y[w] * div_gdYVQVt[a];
      st->diff_flux[w][a] -= Dmu * Y[w] * Y[w] * gammadot * VQVt_grad_mu[a] / mu;
      st->diff_flux[w][a] -= Dd[a]*grad_Y[w][a];

    }


  /* Assemble Jacobian */
  /* Currently no mesh displacement dependencies. */
  if(af->Assemble_Jacobian)
    {



      var = MASS_FRACTION;
      memset(st->d_diff_flux_dc, 0, MAX_CONC*DIM*MAX_CONC*MDE*sizeof(dbl));

		
      /* Compute grad(phi) . VQVt.  This is different for
       * different dependency functions (when
       * bf[var]->grad_phi[j] changes).
       */
      memset(Q_grad_phi, 0, DIM * MDE *sizeof(dbl));
      memset(grad_phi_Q, 0, DIM * MDE * sizeof(dbl));
      
      for(b = 0; b < DIM; b++)
	{
	  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	    {
	      for(a = 0; a < DIM; a++)
		{
		  Q_grad_phi[b][j] += qtensor_loc[b][a] * bf[var]->grad_phi[j][a];
		  grad_phi_Q[b][j] += bf[var]->grad_phi[j][a] * qtensor_loc[a][b];
		}
	    }
	}

      for(a = 0; a < DIM; a++)
	for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	  {
	    c_term = 0.0;
	    mu_term = 0.0;
	    phi_j = bf[var]->phi[j];

	    c_term += phi_j * div_gdYVQVt[a];
	    c_term += Y[w] * gammadot * grad_phi_Q[a][j];
	    c_term += Y[w] * phi_j * (grad_gd_VQVt[a] + gammadot * div_q[a]);
	    c_term *= -Dc;

	    mu_term += 2.0 * phi_j * dmu_dY[w] * VQVt_grad_Y[a];
	    mu_term += (d2mu_dY2[w] - dmu_dY[w] * dmu_dY[w] / mu) * Y[w] * VQVt_grad_Y[a];
	    mu_term += Y[w] * dmu_dY[w] * Q_grad_phi[a][j];
	    mu_term *= -Dmu * gammadot * Y[w] / mu;


	    d_term = -bf[var]->grad_phi[j][a] * Dd[a] 
	      - grad_Y[w][a]*dDd_dy[a]*bf[var]->phi[j]
	      - grad_Y[w][a]*dDd_dgrady[a]*bf[var]->grad_phi[j][a];
	    
	    st->d_diff_flux_dc[w][a] [w][j] = c_term + mu_term + d_term;
	  }

      var = SHEAR_RATE;
      memset(st->d_diff_flux_dSH[w], 0, DIM*MDE*sizeof(dbl));

      /* Compute grad(phi) . VQVt.  This is different for
       * different dependency functions (when
       * bf[var]->grad_phi[j] changes).
       */
      memset(grad_phi_Q, 0, DIM *MDE* sizeof(dbl));
      for(b = 0; b < VIM; b++)
	{
	  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	    {
	      for(a= 0; a < VIM; a++)
		{
		  grad_phi_Q[b][j] += bf[var]->grad_phi[j][a]* qtensor_loc[a][b];
		}
	    }
	}
      for(a = 0; a < VIM; a++)
	{
	  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	    {
	      c_term = 0.0;
	      mu_term = 0.0;
	      phi_j = bf[var]->phi[j];
	      
	      c_term += Y[w] * grad_phi_Q[a][j];
	      c_term += phi_j * (grad_Y_VQVt[a] + Y[w] * div_q[a]);
	      c_term *= -Dc * Y[w];
	      
	      mu_term = -Dmu * Y[w] * Y[w] * VQVt_grad_mu[a] / mu * phi_j;
	      
	      st->d_diff_flux_dSH[w][a][j] = c_term + mu_term;
	    }
	}

      var = VORT_DIR1;
      memset(st->d_diff_flux_dvd[w], 0, DIM*DIM*MDE*sizeof(dbl));
      
      
      memset(d_qtensor_dvd, 0, DIM*DIM*DIM*MDE*sizeof(dbl));
      for(p = 0; p < DIM; p++)
	{
	  for(q = 0; q < DIM; q++)
	    {
	      for(b = 0; b < DIM; b++)
		{
		  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		    {
		      d_qtensor_dvd[p][q][b][j]=-0.5*phi_j*
			(fv->vd[q]*delta(p,b)+fv->vd[p]*delta(q,b));
		    }
		}
	    }
	}
      
      
      memset(d_div_q_dvd, 0, DIM*DIM*MDE*sizeof(dbl));
      for(a = 0; a < DIM; a++)
	{
	  for(b = 0; b < DIM; b++)
	    {
	      for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  div_phi_j_e_b = 0.;
		  for ( p=0; p<VIM; p++)
		    {
		      div_phi_j_e_b += 
			bf[var]->grad_phi_e[j][b] [p][p];
		    }
		  
		  d_div_q_dvd[a][b][j] = -0.5*(phi_j*delta(a,b)*fv->div_vd +
					       fv->vd[a]*div_phi_j_e_b);
		  for ( p=0; p<VIM; p++)
		    {
		      d_div_q_dvd[a][b][j] -= 0.5*phi_j* delta(p,b)*fv->grad_vd[p][a] +
			0.5*fv->vd[p] * bf[var]->grad_phi_e[j][b] [p][a];
		    }
		}
	    }
	}
      
      memset(d_grad_Y_VQVt_dvd,  0, DIM*DIM*MDE*sizeof(dbl));
      memset(d_grad_gd_VQVt_dvd, 0, DIM*DIM*MDE*sizeof(dbl));
      memset(d_VQVt_grad_mu_dvd, 0, DIM*DIM*MDE*sizeof(dbl));
      memset(d_VQVt_grad_Y_dvd,  0, DIM*DIM*MDE*sizeof(dbl));
      for(b = 0; b < DIM; b++)
	{
	  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	    {
	      for(a = 0; a < DIM; a++)
		{
		  for(i = 0; i <DIM; i++)
		    {
		      d_grad_Y_VQVt_dvd[a][b][j] += grad_Y[w][i] * d_qtensor_dvd[i][a][b][j];
		      d_grad_gd_VQVt_dvd[a][b][j] += grad_gammadot[i] * d_qtensor_dvd[i][a][b][j];
		      d_VQVt_grad_mu_dvd[a][b][j] += d_qtensor_dvd[a][i][b][j]* grad_mu[i];
		      d_VQVt_grad_Y_dvd[a][b][j] += d_qtensor_dvd[a][i][b][j] * grad_Y[w][i];
		    }
		}
	    }
	}
      
      
      /* Compute d_div(gammadot * Y * (V Q V^t))_dvd */
      memset(d_div_gdYVQVt_dvd, 0, DIM*DIM*MDE*sizeof(dbl));
      for(b = 0; b < DIM; b++)
	{
	  for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
	    {
	      for(a = 0; a < DIM; a++)
		{
		  d_div_gdYVQVt_dvd[a][b][j] += Y[w] * d_grad_gd_VQVt_dvd[a][b][j];
		  d_div_gdYVQVt_dvd[a][b][j] += gammadot * d_grad_Y_VQVt_dvd[a][b][j];
		  d_div_gdYVQVt_dvd[a][b][j] += gammadot * Y[w] * d_div_q_dvd[a][b][j];
		}
	    }
	}
      
      
      for(a = 0; a < DIM; a++)
	{
	  for(b = 0; b < DIM; b++)
	    {
	      for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  
		  c_term = -Dc * Y[w] * d_div_gdYVQVt_dvd[a][b][j];
		  mu_term = -Dmu * Y[w] * Y[w] * gammadot * d_VQVt_grad_mu_dvd[a][b][j] / mu;
		  
		  
		  st->d_diff_flux_dvd[w][a][b][j] = c_term + mu_term;
		}
	    }
	}
    }
  
 
  return(status);
}


/* Cross the first two really simple dbl[DIM] vectors and return their
 * cross product.
 */
void
cross_really_simple_vectors(const dbl *v0, /* v0 */
			    const dbl *v1, /* v1 */
			    dbl *v2) /* v2 = v0 x v1 */
{
  int i, j, k;

  memset(v2, 0, DIM * sizeof(dbl));
  for(i = 0; i < DIM; i++)
    for(j = 0; j < DIM; j++)
      for(k = 0; k < DIM; k++)
	v2[k] += permute(i,j,k) * v0[i] * v1[j];
}

/* normalize a really simple dbl[dim] vector.  If the norm is 0.0,
 * then the zero vector is left alone.
 */
dbl
normalize_really_simple_vector(dbl *v, int dim)
{
  int i;
  dbl norm;
  
  norm = really_simple_vector_magnitude ( v, dim );

  if(norm != 0.0)
    for(i = 0; i < dim; i++)
      v[i] /= norm;
  return norm;
}

dbl
really_simple_vector_magnitude( dbl *v, int dim )
{
  int i;
  dbl norm;
  
  for (i=0, norm=0.0; i<dim; i++, v++ ) norm += (*v)*(*v);

  return( sqrt(norm) );
}

/* This routine computes the three principle directions used in the Q
 * tensor diffusivity model.  These directions should end up making an
 * orthonormal basis for R^3.  The only exception is when some of the
 * directions are zero, then it shouldn't matter..
 *
 * Upon entrance, v_flow = fv->v, v_vort = approx to d(fv->v)/dt
 *
 * NOTE: Assumes a VIM=3 model.
 *
 * Author: MMH 9/3/99
 */
void
compute_principle_directions(dbl *v_flow,
			     dbl *v_norm,
			     dbl *v_vort,
			     int print)
{
  dbl vf_norm;
  int i,j,k;
  dbl t1;

#ifdef DEBUG_QTENSOR
  if(print)
    {
      printf( "\nUpon entrance, v_flow = [% 10.4g %10.4g %10.4g]\n",
	      v_flow[0], v_flow[1], v_flow[2]);
      printf( "               v_norm = [% 10.4g %10.4g %10.4g]\n",
	      v_norm[0], v_norm[1], v_norm[2]);
      printf( "               v_vort = [% 10.4g %10.4g %10.4g]\n",
	      v_vort[0], v_vort[1], v_vort[2]);
    }
#endif

  t1 = 0.0;
  for(i = 0; i < VIM; i++)
    t1 += v_vort[i] * v_flow[i];

  /* First get a unit vector in the flow direction. */
  vf_norm = normalize_really_simple_vector(v_flow, VIM);

  if(vf_norm != 0.0)
    t1 /= vf_norm;
  else
    t1 = 0.0;

#ifdef DEBUG_QTENSOR
  if(print)
    {
      printf( "Subtracting [% 10.4g, % 10.4g, % 10.4g] from\n",
	      t1*v_flow[0], t1*v_flow[1], t1*v_flow[2]);
      printf( "       from [% 10.4g, % 10.4g, % 10.4g]\n",
	      v_vort[0], v_vort[1], v_vort[2]);
    }
#endif

  for(i = 0; i < VIM; i++)
    v_vort[i] -= t1 * v_flow[i];

#ifdef DEBUG_QTENSOR
  if(print)
    {
      printf( "                    v_flow = [%10.4g, %10.4g, %10.4g]\n", 
	      v_flow[0], v_flow[1], v_flow[2]);
      printf( "Before normalizing, v_vort = [%10.4g, %10.4g, %10.4g]\n", 
	      v_vort[0], v_vort[1], v_vort[2]);
    }
#endif

  normalize_really_simple_vector(v_vort, VIM);

  for(i = 0; i < VIM; i++)
    for(j = 0; j < VIM; j++)
      for(k = 0; k < VIM; k++)
	v_norm[k] += permute(i,j,k) * v_flow[i] * v_vort[j];

#ifdef DEBUG_QTENSOR
  if(print)
    {
      printf( "Upon exit, %10s %10s %10s\n","v_flow","v_norm","v_vort");
      printf( "      (1)  % 10.4g % 10.4g % 10.4g\n",
	      v_flow[0], v_norm[0], v_vort[0]);
      printf( "      (2)  % 10.4g % 10.4g % 10.4g\n",
	      v_flow[1], v_norm[1], v_vort[1]);
      printf( "      (3)  % 10.4g % 10.4g % 10.4g\n",
	      v_flow[2], v_norm[2], v_vort[2]);
    }
#endif
}


/* This routine is used by compute_principle_shear_directions().  It
 * computes an eigenvector of the matrix A, using the eigenvalue
 * lambda.  This does straight Gaussian elimination with partial
 * pivoting.  Since this matrix is so small (3x3), I just do it
 * directly here w/o the outer loop...
 *
 * On entrance, A = 3 x 3 matrix,
 *              lambda = nonzero eigenvalue.
 *
 * On exit, v contains a unit-length eigenvector.
 *
 * Author: MMH 11/8/99
 */
void
find_eigenvector(dbl AA[3][3],
		 dbl lambda,
		 dbl *v,
		 int print)
{
  int i, p[3];
  dbl A[3][3], m12, m13, m23;
  dbl a11, a12, a13, a22, a23;
  dbl x, y, z;

  memcpy(&A[0][0], &AA[0][0], 9 * sizeof(dbl));
  x = y = z = 0.0;
  
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "In find_eigenvector(), memcpy'ed to get AA = \n");
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    printf( "% 16.14g ", A[i][j]);
	  printf( "\n");
	}
      fflush(stdout);
    }
#endif

  /* Initialize the permutation matrix and subtract the eigenvalue
   * from the diagonals.
   */
  for(i = 0; i < 3; i++)
    {
      p[i] = i;
      A[i][i] -= lambda;
    }

  /* Zero out below 0,0 */
#ifdef DEBUG_DIAGONALIZATION
  if(print) { printf( "A  - lambda*I = \n"); for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) printf( "% 10.4g ", A[p[i]][j]); printf( "\n"); }
  printf( "p = ["); for(i = 0; i < 3; i++) printf( "%d ",p[i]);
  printf( "]\n"); fflush(stdout); }  
#endif

  if(fabs(A[0][0]) > fabs(A[1][0]))
    if(fabs(A[0][0]) > fabs(A[2][0]))
      {
	/* Do nothing, the first step can proceed. */
      }
    else
      {
#ifdef DEBUG_DIAGONALIZATION
	if(print) printf( "Swapping first and third rows.\n");
#endif
	p[0] = 2;
	p[2] = 0;
      }
  else
    {
      if(fabs(A[1][0]) > fabs(A[2][0]))
	{
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Swapping first and second rows.\n");
#endif
	  p[0] = 1;
	  p[1] = 0;
	}
      else
	{
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Swapping first and third rows.\n");
#endif
	  p[0] = 2;
	  p[2] = 0;
	}
    }
  if(fabs(A[p[0]][0]) < QTENSOR_SMALL_DBL)
    {
      /* Uh-oh, first column is zero.  If the first minor, A11, is
       * singular we have an eigenspace of dimension >= 2.  If it is
       * not singular, then the eigenvector is (1,0,0).
       */
      if(fabs(A[p[1]][1] * A[p[2]][2] - A[p[1]][2] * A[p[2]][1]) <
	 QTENSOR_SMALL_DBL)
	x = 1.0;
      else
	{ 
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Echeck, double evector\n");
#endif
	}
    }
  else
    {
      m12 = A[p[1]][0] / A[p[0]][0];
#ifdef DEBUG_DIAGONALIZATION
      if(print) printf( "First multiplier is %g\n", m12);
#endif
      A[p[1]][0] = 0.0;
      for(i = 1; i < 3; i++)
	A[p[1]][i] -= m12 * A[p[0]][i];
      m13 = A[p[2]][0] / A[p[0]][0];
#ifdef DEBUG_DIAGONALIZATION
      if(print) printf( "Second multiplier is %g\n", m13);
#endif
      A[p[2]][0] = 0.0;
      for(i = 1; i < 3; i++)
	A[p[2]][i] -= m13 * A[p[0]][i];
      /* Zero out below 1,1 */
#ifdef DEBUG_DIAGONALIZATION
      if(print) { printf( "A - lambda*I = \n"); for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) printf( "% 10.4g ", A[p[i]][j]); printf( "\n"); }
      printf( "p = ["); for(i = 0; i < 3; i++) printf( "%d ",p[i]);
      printf( "]\n"); fflush(stdout); }  
      fflush(stdout);
#endif
      if(fabs(A[p[1]][1]) > fabs(A[p[2]][1]))
	{
	  /* Do nothing, the next elimination step can proceed. */
	}
      else
	{
	  i = p[2];
	  p[2] = p[1];
	  p[1] = i;
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Swapping second and third rows.\n");
#endif
	}
      if(fabs(A[p[1]][1]) < QTENSOR_SMALL_DBL)
	{
	  /* We have a singular matrix, at least one eigenvalue is zero.
	     But DON'T try to do another GE step. */
	}
      else
	{
	  m23 = A[p[2]][1] / A[p[1]][1];
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Third multiplier is %g = %g/%g\n", m23,
			    A[p[2]][1], A[p[1]][1]);
#endif
	  A[p[2]][1] = 0.0;
	  A[p[2]][2] -= m23 * A[p[1]][2];
	}
      /* Swap last row if the second row was the zero-row */
      /* This should never happen, and it uses the 1,2 element, which
       * may actually be zero for perfectly good reasons, so I removed
       * it. */
#if 0
      if(fabs(A[p[2]][2]) > fabs(A[p[1]][2]))
	{
#ifdef DEBUG_DIAGONALIZATION
	  if(print) printf( "Swapping second and third rows b/c |%g| > |%g|\n",
			    A[p[2]][2], A[p[1]][2]);
#endif
	  i = p[2];
	  p[2] = p[1];
	  p[1] = i;
	}
#endif

      /* A[p[2]][*] would be zero analytically, but this procedure has a
       * some precision error in it, so we artifically set the last
       * row to be zero.  
       *
       * This isn't exactly what's done any more.  Since we have an
       * upper triangular matrix, it's really only the 22, 23, and 32
       * elements that we need to investigate. */
      A[p[2]][0] = A[p[2]][1] = 0.0;
      if(fabs(A[p[1]][1]) < fabs(A[p[2]][2]))
	{
	  A[p[1]][1] = 0.0;
	  i = p[2];
	  p[2] = p[1];
	  p[1] = i;
	}
      else
	A[p[2]][2] = 0.0;
#ifdef DEBUG_DIAGONALIZATION
      if(print) { printf( "A - lambda*I = \n"); for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) printf( "% 10.4g ", A[p[i]][j]); printf( "\n"); }
      printf( "p = ["); for(i = 0; i < 3; i++) printf( "%d ",p[i]);
      printf( "]\n"); fflush(stdout); }  
      fflush(stdout);
#endif
      /* Now the matrix is in upper triangular form.  Time to compute the
       * eigenvector.
       */
      a11 = A[p[0]][0];
      a12 = A[p[0]][1];
      a13 = A[p[0]][2];
      a22 = A[p[1]][1];
      a23 = A[p[1]][2];

      /* All of the numerical zeros should have been taken care of
       * before we get here. */
      if(a23 != 0.0)
	if(a22 != 0.0)
	  if(a11 != 0.0)
	    {
	      z = -1;
	      y = a23 / a22;
	      x = (a13 - a12 * y) / a11;
#ifdef DEBUG_DIAGONALIZATION
	      if(print) printf( "Echeck, full evalue\n");
#endif
	    }
	  else
	    {
	      y = a23;
	      x = -a22;
#ifdef DEBUG_DIAGONALIZATION
	      if(print) printf( "Echeck, rel rat err = %g\n", fabs((a12/a13 - a22/a23)/(a22/a23)));
#endif
	    }
	else
	  if(a12 != 0.0)
	    if(a11 != 0.0)
	      {
		x = a12;
		y = -a11;
	      }
	    else
	      {
#ifdef DEBUG_DIAGONALIZATION
		if(print) printf( "Echeck, 0 evector\n");
#endif
	      }
	  else
	    if(a11 != 0.0)
	      y = 1.0;
	    else
	      {
#ifdef DEBUG_DIAGONALIZATION
		if(print) printf( "Echeck, double evector\n");
#endif
	      }
      else
	{
	  if(a22 != 0.0)
	    {
	      if(a13 != 0.0)
		{
		  if(a11 != 0.0)
		    {
		      x = a13;
		      z = -a11;
		    }
		  else
		    x = 1.0;
		}
	      else
		{
		  if(a11 != 0.0)
		    z = 1.0;
#ifdef DEBUG_DIAGONALIZATION
		  else
		    { if(print) printf( "Echeck, double evector\n"); }
#endif
		}
	    }
#ifdef DEBUG_DIAGONALIZATION
	  else
	    { if(print) printf( "Echeck, double evector\n"); }
#endif
	}
    }
  
#ifdef DEBUG_DIAGONALZIATION
  printf( "x=%g, y=%g, z=%g\n", x, y, z);
#endif

  v[0] = x;
  v[1] = y;
  v[2] = z;
  normalize_really_simple_vector(v, 3);

#ifdef DEBUG_QTENSOR
  if(print)
    printf( "EIGENPAIR: { %.8g, (%.8g, %.8g, %.8g) }\n", lambda,
	    v[0],v[1],v[2]);
#endif
}


/* This routine returns the eigenvalue of smallest absolute value and
 * its associated eigenvector.
 *
 * NOTE: Assumes a VIM=3 model.
 *
 * Author: MMH 11/99-3/01 */

void
find_super_special_eigenvector(dbl T[DIM][DIM],
			       dbl *v,
			       dbl *v1,
			       dbl *v2,
			       dbl *v3,
			       dbl *eigenvalue,
			       int print)
{
  int num_zero_eigenvalues;
  int i,j;
  dbl A[3][3];
  dbl a0, a2, a1;
  dbl q, r, d, m, theta, z1, z2, z3;
  dbl eig1=0., eig2=0., eig3=0.;

  memset(v, 0, DIM * sizeof(dbl));
  memset(v1, 0, DIM * sizeof(dbl));
  memset(v2, 0, DIM * sizeof(dbl));
  memset(v3, 0, DIM * sizeof(dbl));
 
  *eigenvalue = 0.0;

  memset(A, 0, 9 * sizeof(dbl));
  for (i=0; i < DIM; i++)
    {
      for (j=0; j < DIM; j++)
	{
	  A[i][j] = T[i][j];
	}
    }

  /* Now find the zeros of the cubic polynomial characteristic
   * equation.  This requires some complex number mumbo-jubmo, which I
   * have significantly simplified on some assumptions (symmetric
   * tensor, etc.).
   */

  /* Note that if T really is a rate of deformation tensor, then a2 =
   * 0 because it is just the trace, and the trace should be zero if
   * the flow is incompressible. */
  get_characteristic_eq_coeffs(A, &a0, &a1, &a2);

  q = (3.0 * a1 - a2 * a2) / 9.0;
  /* q will only be > 0 due to roundoff. */
  if(q > 0.0)
    q = 0.0;
  r = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 *a2 / 27.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "q = % 10.4g, r = % 10.4g ", q, r);
      fflush(stdout);
    }
#endif
  /* This will only fail if q and r are corrupted by roundoff. */
  d = -q*q*q - r*r;
  if(d < 0.0)
    d = 0.0;
  else
    d = sqrt(d);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "d = % 10.4g\n", d);
#endif
  theta = atan2(d, r) / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "theta = % 10.4g, ", theta);
      fflush(stdout);
    }
#endif
  if(q == 0.0)
    m = 0.0;
  else
    m = sqrt(-q);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "m = % 10.4g\n", m);
#endif
  z1 = 2.0 * m * cos(theta) - a2 / 3.0;
  z2 = -m * (cos(theta) + sqrt(3.0) * sin(theta)) - a2 / 3.0;
  z3 = -m * (cos(theta) - sqrt(3.0) * sin(theta)) - a2 / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "z1 = % 16.14g, z2 = % 16.14g, z3 = % 16.14g\n", z1, z2, z3);
  fflush(stdout);
#endif

  /* Try to catch some roundoff errors. */
  if(fabs(z1) < QTENSOR_SMALL_DBL) z1 = 0.0;
  if(fabs(z2) < QTENSOR_SMALL_DBL) z2 = 0.0;
  if(fabs(z3) < QTENSOR_SMALL_DBL) z3 = 0.0;

  num_zero_eigenvalues = (z1 == 0.0) + (z2 == 0.0) + (z3 == 0.0);
  if(print)
    {
#ifdef DEBUG_QTENSOR
      printf( "EIGENVALUES = %g, %g, %g\n", z1, z2, z3);
#endif
#ifdef DEBUG_DIAGONALIZATION
      printf( "Before entering find_eigenvector(), A = \n");
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    printf( "% 10.4g ", A[i][j]);
	  printf( "\n");
	}
      fflush(stdout);
#endif
    }

  /* Arrange eigenvalues as tension, compression, vorticity */
  /* eig(compression) < eig(vorticity) < eig(tension)       */
  
  if(num_zero_eigenvalues <= 1)
    {
      /* Good, there must be a "largest" and "smallest" eigenvalue. */
      if(z1 > z2)
	if(z2 > z3)
	  { 
	    /* z3 < z2 < z1 */
	    find_eigenvector(A, z2, v, print);
	    *eigenvalue = z2;
	    eig1 = z1;
	    eig2 = z3;
	    eig3 = z2;
	  }
	else
	  if(z3 > z1)
	    {
	      /* z2 < z1 < z3 */
	      find_eigenvector(A, z1, v, print);
	      *eigenvalue = z1;
	      eig1 = z3;
	      eig2 = z2;
	      eig3 = z1;
	    }
	  else
	    {
	      /* z2 < z3 < z1 */
	      find_eigenvector(A, z3, v, print);
	      *eigenvalue = z3;
	      eig1 = z1;
	      eig2 = z2;
	      eig3 = z3;
	    }
      else
	if(z2 > z3)
	  if(z3 > z1)
	    {
	      /* z1 < z3 < z2 */
	      find_eigenvector(A, z3, v, print);
	      *eigenvalue = z3;
	      eig1 = z2;
	      eig2 = z1;
	      eig3 = z3;
	    }
	  else
	    {
	      /* z3 < z1 < z2 */
	      find_eigenvector(A, z1, v, print);
	      *eigenvalue = z1;
	      eig1 = z2;
	      eig2 = z3;
	      eig3 = z1;
	    }
	else
	  {
	    /* z1 < z2 < z3 */
	    find_eigenvector(A, z2, v, print);
	    *eigenvalue = z2;
	    eig1 = z3;
	    eig2 = z1;
	    eig3 = z2;
	  }

      if(print)
	{
	  printf( "special eigenvalue = % 10.4g\n", *eigenvalue);
	  printf( "                 v = [% 10.4g % 10.4g % 10.4g]\n",
		  v[0], v[1], v[2]);
	}
    }
  else
    {
      /* At least 2 zeros => all 3 are zero (symmetry) */
    }
  find_eigenvector(A, eig1, v1, print);
  find_eigenvector(A, eig2, v2, print);
  find_eigenvector(A, eig3, v3, print);
  /* Try to catch some roundoff errors. */
  if(fabs(v[0]) < QTENSOR_SMALL_DBL) v[0] = 0.0;
  if(fabs(v[1]) < QTENSOR_SMALL_DBL) v[1] = 0.0;
  if(fabs(v[2]) < QTENSOR_SMALL_DBL) v[2] = 0.0;
}

void
find_eigenvalues_eigenvectors(dbl T[DIM][DIM],
			      dbl *e1, dbl *e2, dbl *e3,
			       dbl *v1,
			       dbl *v2,
			       dbl *v3)
{
  int num_zero_eigenvalues;
  int i,j, print=0;
  dbl A[3][3];
  dbl a0, a2, a1;
  dbl q, r, d, m, theta, z1, z2, z3;
  dbl eig1=0., eig2=0., eig3=0.;
  dbl v[3];

  memset(v, 0, DIM * sizeof(dbl));
  memset(v1, 0, DIM * sizeof(dbl));
  memset(v2, 0, DIM * sizeof(dbl));
  memset(v3, 0, DIM * sizeof(dbl));

  memset(A, 0, 9 * sizeof(dbl));
  for (i=0; i < DIM; i++)
    {
      for (j=0; j < DIM; j++)
	{
	  A[i][j] = T[i][j];
	}
    }

  /* Now find the zeros of the cubic polynomial characteristic
   * equation.  This requires some complex number mumbo-jubmo, which I
   * have significantly simplified on some assumptions (symmetric
   * tensor, etc.).
   */

  /* Note that if T really is a rate of deformation tensor, then a2 =
   * 0 because it is just the trace, and the trace should be zero if
   * the flow is incompressible. */
  get_characteristic_eq_coeffs(A, &a0, &a1, &a2);

  q = (3.0 * a1 - a2 * a2) / 9.0;
  /* q will only be > 0 due to roundoff. */
  if(q > 0.0)
    q = 0.0;
  r = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 *a2 / 27.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "q = % 10.4g, r = % 10.4g ", q, r);
      fflush(stdout);
    }
#endif
  /* This will only fail if q and r are corrupted by roundoff. */
  d = -q*q*q - r*r;
  if(d < 0.0)
    d = 0.0;
  else
    d = sqrt(d);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "d = % 10.4g\n", d);
#endif
  theta = atan2(d, r) / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "theta = % 10.4g, ", theta);
      fflush(stdout);
    }
#endif
  if(q == 0.0)
    m = 0.0;
  else
    m = sqrt(-q);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "m = % 10.4g\n", m);
#endif
  z1 = 2.0 * m * cos(theta) - a2 / 3.0;
  z2 = -m * (cos(theta) + sqrt(3.0) * sin(theta)) - a2 / 3.0;
  z3 = -m * (cos(theta) - sqrt(3.0) * sin(theta)) - a2 / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "z1 = % 16.14g, z2 = % 16.14g, z3 = % 16.14g\n", z1, z2, z3);
  fflush(stdout);
#endif

  /* Try to catch some roundoff errors. */
  if(fabs(z1) < QTENSOR_SMALL_DBL) z1 = 0.0;
  if(fabs(z2) < QTENSOR_SMALL_DBL) z2 = 0.0;
  if(fabs(z3) < QTENSOR_SMALL_DBL) z3 = 0.0;

  num_zero_eigenvalues = (z1 == 0.0) + (z2 == 0.0) + (z3 == 0.0);
  if(print)
    {
#ifdef DEBUG_QTENSOR
      printf( "EIGENVALUES = %g, %g, %g\n", z1, z2, z3);
#endif
#ifdef DEBUG_DIAGONALIZATION
      printf( "Before entering find_eigenvector(), A = \n");
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    printf( "% 10.4g ", A[i][j]);
	  printf( "\n");
	}
      fflush(stdout);
#endif
    }

  /* Arrange eigenvalues as tension, compression, vorticity */
  /* eig(compression) < eig(vorticity) < eig(tension)       */
  
  if(num_zero_eigenvalues <= 1)
    {
      /* Good, there must be a "largest" and "smallest" eigenvalue. */
      if(z1 > z2)
	if(z2 > z3)
	  { 
	    /* z3 < z2 < z1 */
	    find_eigenvector(A, z2, v, print);
	    eig1 = z1;
	    eig2 = z3;
	    eig3 = z2;
	  }
	else
	  if(z3 > z1)
	    {
	      /* z2 < z1 < z3 */
	      find_eigenvector(A, z1, v, print);
	      eig1 = z3;
	      eig2 = z2;
	      eig3 = z1;
	    }
	  else
	    {
	      /* z2 < z3 < z1 */
	      find_eigenvector(A, z3, v, print);
	      eig1 = z1;
	      eig2 = z2;
	      eig3 = z3;
	    }
      else
	if(z2 > z3)
	  if(z3 > z1)
	    {
	      /* z1 < z3 < z2 */
	      find_eigenvector(A, z3, v, print);
	      eig1 = z2;
	      eig2 = z1;
	      eig3 = z3;
	    }
	  else
	    {
	      /* z3 < z1 < z2 */
	      find_eigenvector(A, z1, v, print);
	      eig1 = z2;
	      eig2 = z3;
	      eig3 = z1;
	    }
	else
	  {
	    /* z1 < z2 < z3 */
	    find_eigenvector(A, z2, v, print);
	    eig1 = z3;
	    eig2 = z1;
	    eig3 = z2;
	  }
    }
  else
    {
      /* At least 2 zeros => all 3 are zero (symmetry) */
    }
  find_eigenvector(A, eig1, v1, print);
  find_eigenvector(A, eig2, v2, print);
  find_eigenvector(A, eig3, v3, print);
  /* Try to catch some roundoff errors. */
  if(fabs(v[0]) < QTENSOR_SMALL_DBL) v[0] = 0.0;
  if(fabs(v[1]) < QTENSOR_SMALL_DBL) v[1] = 0.0;
  if(fabs(v[2]) < QTENSOR_SMALL_DBL) v[2] = 0.0;

  *e1 = eig1;
  *e2 = eig2;
  *e3 = eig3;
}



/* This routine diagonalizes a symmetric 2-tensor of size 3 by 3.
 * Usually, this is the rate-of-deformation tensor.  In that case,
 * this routine computes the three principle shear directions used in
 * the Q tensor diffusivity model.  It returns the eigenvalues, e[0] <
 * e[1] < e[2], and their associated eigenvectors v0, v1, v2.  Note
 * that if the rate-of-deformation tensor is of rank 1 or less that
 * everything is just returned as zeros.  It copies the incoming
 * tensor so as not to change it.
 *
 * NOTE: Assumes a VIM=3 model.
 *
 * Author: MMH 11/99-12/99 */
void
diagonalize_symmetric_tensor(dbl T[DIM][DIM],
			     dbl *v0,
			     dbl *v1,
			     dbl *v2,
			     dbl *eigenvalues,
			     int print)
{
  int i, j, k;
  int num_zero_eigenvalues;
  dbl A[3][3];
  dbl a0, a2, a1;
  dbl q, r, d, m, theta, z1, z2, z3;
  dbl tmp1, vfnorm;

  memset(v0, 0, DIM * sizeof(dbl));
  memset(v1, 0, DIM * sizeof(dbl));
  memset(v2, 0, DIM * sizeof(dbl));
  memset(eigenvalues, 0, DIM * sizeof(dbl));

  /* First, get a copy of the shear rate tensor. */
  memcpy(A, T, DIM * DIM * sizeof(dbl));

  /* Now find the zeros of the cubic polynomial characteristic
   * equation.  This requires some complex number mumbo-jubmo, which I
   * have significantly simplified on some assumptions (symmetric
   * tensor, etc.).
   */

  /* Note that if T really is a rate of deformation tensor, then a2 =
   * 0 because it is just the trace, and the trace should be zero if
   * the flow is incompressible. */
  a2 = -A[0][0] - A[1][1] - A[2][2];
  a1 = A[1][1] * A[2][2] + A[0][0] * A[2][2] + A[0][0] * A[1][1] -
    A[0][1] * A[1][0] - A[0][2] * A[2][0] - A[1][2] * A[2][1];
  a0 = -A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[0][1] * A[2][2] +
    A[0][0] * A[1][2] * A[1][2] + A[0][2] * A[0][2] * A[1][1] -
    A[0][1] * A[2][0] * A[1][2] - A[0][2] * A[1][0] * A[2][1];

  q = (3.0 * a1 - a2 * a2) / 9.0;
  /* q will only be > 0 due to roundoff. */
  if(q > 0.0)
    q = 0.0;
  r = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 *a2 / 27.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "q = % 10.4g, r = % 10.4g ", q, r);
      fflush(stdout);
    }
#endif
  /* This will only fail if q and r are corrupted by roundoff. */
  d = -q*q*q - r*r;
  if(d < 0.0)
    d = 0.0;
  else
    d = sqrt(d);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "d = % 10.4g\n", d);
#endif
  theta = atan2(d, r) / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "theta = % 10.4g, ", theta);
      fflush(stdout);
    }
#endif
  if(q == 0.0)
    m = 0.0;
  else
    m = sqrt(-q);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "m = % 10.4g\n", m);
#endif
  z1 = 2.0 * m * cos(theta) - a2 / 3.0;
  z2 = -m * (cos(theta) + sqrt(3.0) * sin(theta)) - a2 / 3.0;
  z3 = -m * (cos(theta) - sqrt(3.0) * sin(theta)) - a2 / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "z1 = % 16.14g, z2 = % 16.14g, z3 = % 16.14g\n", z1, z2, z3);
  fflush(stdout);
#endif

  /* Try to catch some roundoff errors. */
  if(fabs(z1) < QTENSOR_SMALL_DBL) z1 = 0.0;
  if(fabs(z2) < QTENSOR_SMALL_DBL) z2 = 0.0;
  if(fabs(z3) < QTENSOR_SMALL_DBL) z3 = 0.0;

  num_zero_eigenvalues = (z1 == 0.0) + (z2 == 0.0) + (z3 == 0.0);
  if(print)
    {
#ifdef DEBUG_QTENSOR
      printf( "EIGENVALUES = %g, %g, %g\n", z1, z2, z3);
#endif
#ifdef DEBUG_DIAGONALIZATION
      printf( "Before entering find_eigenvector(), A = \n");
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    printf( "% 10.4g ", A[i][j]);
	  printf( "\n");
	}
      fflush(stdout);
#endif
    }

  if(num_zero_eigenvalues <= 1)
    {
      /* Good, there must be a "largest" and "smallest" eigenvalue. */
      if(z1 > z2)
	if(z2 > z3)
	  { 
	    /* z3 < z2 < z1 */
	    find_eigenvector(A, z3, v0, print);
	    find_eigenvector(A, z1, v2, print);
	    eigenvalues[0] = z3;
	    eigenvalues[1] = z2;
	    eigenvalues[2] = z1;
	  }
	else
	  if(z3 > z1)
	    {
	      /* z2 < z1 < z3 */
	      find_eigenvector(A, z2, v0, print);
	      find_eigenvector(A, z3, v2, print);
	      eigenvalues[0] = z2;
	      eigenvalues[1] = z1;
	      eigenvalues[2] = z3;
	    }
	  else
	    {
	      /* z2 < z3 < z1 */
	      find_eigenvector(A, z2, v0, print);
	      find_eigenvector(A, z1, v2, print);
	      eigenvalues[0] = z2;
	      eigenvalues[1] = z3;
	      eigenvalues[2] = z1;
	    }
      else
	if(z2 > z3)
	  if(z3 > z1)
	    {
	      /* z1 < z3 < z2 */
	      find_eigenvector(A, z1, v0, print);
	      find_eigenvector(A, z2, v2, print);
	      eigenvalues[0] = z1;
	      eigenvalues[1] = z3;
	      eigenvalues[2] = z2;
	    }
	  else
	    {
	      /* z3 < z1 < z2 */
	      find_eigenvector(A, z3, v0, print);
	      find_eigenvector(A, z2, v2, print);
	      eigenvalues[0] = z3;
	      eigenvalues[1] = z1;
	      eigenvalues[2] = z2;
	    }
	else
	  {
	    /* z1 < z2 < z3 */
	    find_eigenvector(A, z1, v0, print);
	    find_eigenvector(A, z3, v2, print);
	    eigenvalues[0] = z1;
	    eigenvalues[1] = z2;
	    eigenvalues[2] = z3;
	  }
      /*
      if(print)
	for(i = 0; i < DIM; i++)
	  {
	    printf( "v0[%d] = %g\n", i, v0[i]);
	    printf( "v2[%d] = %g\n", i, v2[i]);
	  }
      */

      /* tmp1 holds the value of v0 . v2.  It is used to subtract the
       * part of v2 that is in the v0 direction.  We dot first, then
       * convert v0 into a unit vector.  Once v2 is orthogonal to v0
       * we normalize v2.  The cross product of these two othronormal
       * directions is the third eigenvector (assuming the incoming
       * matrix was real symmetric).
       *
       * I'm assuming that I don't compute the eigenvector for the
       * minimum eigenvalue (in absolute value) directly because of
       * stability issues.  Check this sometime.
       */
      tmp1 = 0.0;
      for(i = 0; i < DIM; i++)
	tmp1 += v0[i] * v2[i];
      vfnorm = normalize_really_simple_vector(v0, DIM);
      if(vfnorm > 1.0e-16)
	tmp1 /= vfnorm;
      for(i = 0; i < DIM; i++)
	{
	  /* was (10/18/00):
	  v1[i] -= tmp1 * v0[i];
	  */
	  v2[i] -= tmp1 * v0[i];
	  /*
	  if(print)
	    printf( "Setting v2[%d] = %g\n", i, v2[i]);
	  */
	}
      normalize_really_simple_vector(v2, DIM);
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  for(k = 0; k < DIM; k++)
	    v1[k] += permute(i,j,k) * v0[i] * v2[j];
      /*
      if(print)
	for(i = 0; i < DIM; i++)
	  printf( "v1[%d] = %g\n", i, v1[i]);
      */
    }
  else
    {
      /* At least 2 zeros => all 3 are zero (symmetry) */
    }
}

/* This routine diagonalizes a rate of deformation tensor.  I've tried
 * to incorporate as many simplifying assumptions as possible
 * (symmetric, trace is zero b/c of incompressibility, etc., etc.,
 * etc.).  This routine computes the three principle shear directions
 * used in the Q tensor diffusivity model.  It returns the
 * eigenvalues, e[0] < e[1] < e[2], and their associated eigenvectors
 * v0, v1, v2.  Note that if the rate-of-deformation tensor is of rank
 * 1 or less that everything is just returned as zeros.  It copies the
 * incoming tensor so as not to change it.
 *
 * NOTE: Assumes a VIM=3 model.
 *
 * Author: MMH 11/99-2/01 */
void
diagonalize_rate_of_deformation_tensor(dbl T[3][3],
				       dbl *v0,
				       dbl *v1,
				       dbl *v2,
				       dbl *eigenvalues,
				       int print)
{
  int i, j, k;
  int num_zero_eigenvalues;
  dbl A[3][3];
  dbl a0, a1, a2;		/* a2 is ignored b/c tr(T) = 0 */
  dbl q, r, d, m, theta, costheta, sintheta, z1, z2, z3;
  dbl tmp1, vfnorm;

  memset(v0, 0, DIM * sizeof(dbl));
  memset(v1, 0, DIM * sizeof(dbl));
  memset(v2, 0, DIM * sizeof(dbl));
  memset(eigenvalues, 0, DIM * sizeof(dbl));

  /* First, get a copy of the shear rate tensor. */
  memcpy(A, T, DIM * DIM * sizeof(dbl));

  /* Now find the zeros of the cubic polynomial characteristic
   * equation.  This requires some complex number mumbo-jubmo, which I
   * have significantly simplified on some assumptions (symmetric
   * tensor, etc.).
   */
  get_characteristic_eq_coeffs(A, &a0, &a1, &a2);

  q = a1 / 3.0;
  /* q will only be > 0 due to roundoff. */
  if(q > 0.0)
    q = 0.0;
  r = -a0 / 2.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "q = % 10.4g, r = % 10.4g ", q, r);
      fflush(stdout);
    }
#endif
  /* This will only fail if q and r are corrupted by roundoff. */
  d = -q*q*q - r*r;
  if(d < 0.0)
    d = 0.0;
  else
    d = sqrt(d);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "d = % 10.4g\n", d);
#endif
  theta = atan2(d, r) / 3.0;
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    {
      printf( "theta = % 10.4g, ", theta);
      fflush(stdout);
    }
#endif
  if(q == 0.0)
    m = 0.0;
  else
    m = sqrt(-q);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "m = % 10.4g\n", m);
#endif
  costheta = cos(theta);
  sintheta = sin(theta);
  z1 = 2.0 * m * costheta;
  z2 = -m * (costheta + sqrt(3.0) * sintheta);
  z3 = -m * (costheta - sqrt(3.0) * sintheta);
#ifdef DEBUG_DIAGONALIZATION
  if(print)
    printf( "z1 = % 16.14g, z2 = % 16.14g, z3 = % 16.14g\n", z1, z2, z3);
  fflush(stdout);
#endif

  /* Try to catch some roundoff errors. */
  if(fabs(z1) < QTENSOR_SMALL_DBL) z1 = 0.0;
  if(fabs(z2) < QTENSOR_SMALL_DBL) z2 = 0.0;
  if(fabs(z3) < QTENSOR_SMALL_DBL) z3 = 0.0;

  num_zero_eigenvalues = (z1 == 0.0) + (z2 == 0.0) + (z3 == 0.0);
  if(print)
    {
#ifdef DEBUG_QTENSOR
      printf( "EIGENVALUES = %g, %g, %g\n", z1, z2, z3);
#endif
#ifdef DEBUG_DIAGONALIZATION
      printf( "Before entering find_eigenvector(), A = \n");
      for(i = 0; i < 3; i++)
	{
	  for(j = 0; j < 3; j++)
	    printf( "% 10.4g ", A[i][j]);
	  printf( "\n");
	}
      fflush(stdout);
#endif
    }

  if(num_zero_eigenvalues <= 1)
    {
      /* Good, there must be a "largest" and "smallest" eigenvalue. */
      if(z1 > z2)
	if(z2 > z3)
	  { 
	    /* z3 < z2 < z1 */
	    find_eigenvector(A, z3, v0, print);
	    find_eigenvector(A, z1, v2, print);
	    eigenvalues[0] = z3;
	    eigenvalues[1] = z2;
	    eigenvalues[2] = z1;
	  }
	else
	  if(z3 > z1)
	    {
	      /* z2 < z1 < z3 */
	      find_eigenvector(A, z2, v0, print);
	      find_eigenvector(A, z3, v2, print);
	      eigenvalues[0] = z2;
	      eigenvalues[1] = z1;
	      eigenvalues[2] = z3;
	    }
	  else
	    {
	      /* z2 < z3 < z1 */
	      find_eigenvector(A, z2, v0, print);
	      find_eigenvector(A, z1, v2, print);
	      eigenvalues[0] = z2;
	      eigenvalues[1] = z3;
	      eigenvalues[2] = z1;
	    }
      else
	if(z2 > z3)
	  if(z3 > z1)
	    {
	      /* z1 < z3 < z2 */
	      find_eigenvector(A, z1, v0, print);
	      find_eigenvector(A, z2, v2, print);
	      eigenvalues[0] = z1;
	      eigenvalues[1] = z3;
	      eigenvalues[2] = z2;
	    }
	  else
	    {
	      /* z3 < z1 < z2 */
	      find_eigenvector(A, z3, v0, print);
	      find_eigenvector(A, z2, v2, print);
	      eigenvalues[0] = z3;
	      eigenvalues[1] = z1;
	      eigenvalues[2] = z2;
	    }
	else
	  {
	    /* z1 < z2 < z3 */
	    find_eigenvector(A, z1, v0, print);
	    find_eigenvector(A, z3, v2, print);
	    eigenvalues[0] = z1;
	    eigenvalues[1] = z2;
	    eigenvalues[2] = z3;
	  }

      /* Before gramm-schmidt'ing them ... */
      /*
      if(print)
	for(i = 0; i < DIM; i++)
	  {
	    printf( "v0[%d] = %g\n", i, v0[i]);
	    printf( "v2[%d] = %g\n", i, v2[i]);
	  }
      */

      /* tmp1 holds the value of v0 . v2.  It is used to subtract the
       * part of v2 that is in the v0 direction.  We dot first, then
       * convert v0 into a unit vector.  Once v2 is orthogonal to v0
       * we normalize v2.  The cross product of these two othronormal
       * directions is the third eigenvector (assuming the incoming
       * matrix was real symmetric).
       *
       * I'm assuming that I don't compute the eigenvector for the
       * minimum eigenvalue (in absolute value) directly because of
       * stability issues.  Check this sometime.
       */
      tmp1 = 0.0;
      for(i = 0; i < DIM; i++)
	tmp1 += v0[i] * v2[i];
      vfnorm = normalize_really_simple_vector(v0, DIM);
      if(vfnorm > 1.0e-16)
	tmp1 /= vfnorm;
      for(i = 0; i < DIM; i++)
	{
	  /* was (10/18/00):
	  v1[i] -= tmp1 * v0[i];
	  */
	  v2[i] -= tmp1 * v0[i];
	  /*
	  if(print)
	    printf( "Setting v2[%d] = %g\n", i, v2[i]);
	  */
	}
      normalize_really_simple_vector(v2, DIM);
      for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	  for(k = 0; k < DIM; k++)
	    v1[k] += permute(i,j,k) * v0[i] * v2[j];

      if(print)
	{
	  printf( "l0, v0 = % 10.4g, [% 10.4g, % 10.4g, % 10.4g]\n",
		  eigenvalues[0], v0[0], v0[1], v0[2]);
	  printf( "l1, v1 = % 10.4g, [% 10.4g, % 10.4g, % 10.4g]\n",
		  eigenvalues[1], v1[0], v1[1], v1[2]);
	  printf( "l2, v2 = % 10.4g, [% 10.4g, % 10.4g, % 10.4g]\n",
		  eigenvalues[2], v2[0], v2[1], v2[2]);
	}

    }
  else
    {
      /* At least 2 zeros => all 3 are zero (symmetry) */
    }
}

/* Ill-fated attempt to evaluate Brady's triple integral formulation
   for the Q-tensor by Gaussian quadrature.  */
void compute_VQVt_directly(dbl T[3][3],
			   dbl VQVt[3][3],
			   int print)
{
  int i, j, k, l, m, g1, g2, g3, G1, G2, G3;
  int num_nonzero_eigenvalues;
  dbl a, b, sum, minisum, prod, jac;
  dbl low1, low2, low3, high1, high2, high3;
  dbl v0[3], v1[3], v2[3], lambda[3], x[3];
  dbl S[3][3], D[3][3];

  memset(VQVt, 0, DIM * DIM * sizeof(dbl));
  memcpy(D, T, DIM * DIM * sizeof(dbl));

  diagonalize_rate_of_deformation_tensor(D, v0, v1, v2, lambda, print);

  num_nonzero_eigenvalues = (lambda[0] == 0.0) + (lambda[1] == 0.0) + (lambda[2] == 0.0);

  if(num_nonzero_eigenvalues > 1)
    return;

  /* Get the orthogonal transformation that diagonalizes the
     rate-of-deformation tensor, S = [v0 v1 v2]. */
  for(i = 0; i < VIM; i++)
    {
      S[i][0] = v0[i];	
      S[i][1] = v1[i];
      S[i][2] = v2[i];
    }

  /* Compute major and minor axes of integral domain. */
  a = (lambda[1] + lambda[2]) / (2.0 * lambda[1] + lambda[2]);
  b = (lambda[1] + lambda[2]) / (lambda[1] + 2.0 * lambda[2]);

#ifdef DEBUG_QTENSOR
  if(print)
    printf( "a = %g, b = %g\n", a, b);
#endif
  
  /* Set number of gauss points for each integral. */
  G1 = G2 = G3 = 2;

  for(i = 0; i < VIM; i++)
    for(j = 0; j < VIM; j++)
      {
	sum = 0.0;
	high3 = sqrt(b);
	low3 = -high3;
#ifdef DEBUG_QTENSOR
	if(print) printf( "high3 = %g\n", high3);
#endif
	for(g3 = 0; g3 < G3; g3++)
	  {
	    x[2] = gauss_point[G3 - 1][g3] * high3;
	    high2 = sqrt(a - a / b * x[2] * x[2]);
	    low2 = -high2;
#ifdef DEBUG_QTENSOR
	    if(print) printf( "high2 = %g\n", high2);
#endif
	    for(g2 = 0; g2 < G2; g2++)
	      {
		x[1] = gauss_point[G2 - 1][g2] * high2;
		high1 = sqrt(1.0 - x[1] * x[1] - x[2] * x[2]);
		low1 = -high1;
#ifdef DEBUG_QTENSOR
		if(print) printf( "high1 = %g\n", high1);
#endif
		jac = (high3 - low3) * (high2 - low2) * (high1 - low1) / 8.0;
#ifdef DEBUG_QTENSOR
		if(print) printf( "jac = %g\n", jac);
#endif
		for(g1 = 0; g1 < G1; g1++)
		  {
		    x[0] = gauss_point[G1 - 1][g1] * high1;
#ifdef DEBUG_QTENSOR
		    if(print) printf( "x3 = %g, x2 = %g, x1 = %g\n", x[2], x[1], x[0]);
#endif
		    for(k = 0; k < VIM; k++)
		      for(l = 0; l < VIM; l++)
			{
			  minisum = 0.0;
			  for(m = 0; m < VIM; m++)
			    minisum += lambda[m] * x[m] * x[m];
			  prod = x[k] * x[l] * minisum * minisum;
			  prod *= S[k][i] * S[l][j];
			  sum += gauss_weight[G3 - 1][g3] *
			    gauss_weight[G2 - 1][g2] * 
			    gauss_weight[G1 - 1][g1] * prod * jac;
			}
		  }
		  
	      }
	  }
	VQVt[i][j] = sum;
      }
  diagonalize_symmetric_tensor(VQVt, v0, v1, v2, lambda, print);
#ifdef DEBUG_QTENSOR
  if(print && ei[pg->imtrx]->ielem == 0)
    {
      printf( "Evalues = %g (flow), %g (norm), %g (vort)\n", 
	      lambda[0], lambda[1], lambda[2]);
      printf( "Flow direction = [% 10.4g, %10.4g, %10.4g]\n", v0[0], v0[1], v0[2]);
      printf( "Norm direction = [% 10.4g, %10.4g, %10.4g]\n", v1[0], v1[1], v1[2]);
      printf( "Vort direction = [% 10.4g, %10.4g, %10.4g]\n", v2[0], v2[1], v2[2]);
    }
#endif
}
  

/* This is yet another attempt to find the normal direction to the
 * plane of shear. 
 *
 * Author: Matt Hopkins (9114) 12/16/99
 */
void
please_work(dbl T[3][3],	/* Rate of deformation tensor */
	    dbl *v_flow,	/* first direction in plane of shear */
	    dbl *v_norm,	/* direction normal to plane of shear */
	    dbl *v_vort,	/* second direction in plane of shear */
	    int print)		/* print toggle */
{
  int i;
  dbl v0[3], v2[3], evalues[3];

  for(i = 0; i < DIM; i++)
    v_flow[i] = fv->v[i];
  normalize_really_simple_vector(v_flow, DIM);

  /* v1 will be the eigenvector with smallest |eigenvalue| since the
   * tensor is symmetric.
   */
  diagonalize_rate_of_deformation_tensor(T, v0, v_vort, v2, evalues, print);
  /*
  cross_really_simple_vectors(v_flow, v_vort, v_norm);
  */
  cross_really_simple_vectors(v_vort, v_flow, v_norm);
#ifdef DEBUG_QTENSOR
  if(print)
    {
      printf( "\nLEAVING, v_flow = [% 10.4g %10.4g %10.4g]\n", 
	      v_flow[0], v_flow[1], v_flow[2]);
      printf( "         v_norm = [% 10.4g %10.4g %10.4g]\n",
	      v_norm[0], v_norm[1], v_norm[2]);
      printf( "         v_vort = [% 10.4g %10.4g %10.4g]\n",
	      v_vort[0], v_vort[1], v_vort[2]);
    }
#endif
}

/* I am making this short snippet a function so I cannot inadvertently
 * use two different (and thus wrong) computations for the
 * coefficients of the charactersitic equation.  BTW, the
 * characteristic equation for a rate of deformation tensor comes out
 * to be l^3 + a1 l + a0 = 0.  a2 = zero b/c it is equal to tr(E),
 * which is zero for incompressible fluids.  At least with a correct
 * solution...  */
void
get_characteristic_eq_coeffs(dbl E[3][3],
			     dbl *a0,
			     dbl *a1,
			     dbl *a2)
{
  *a2 = -(E[0][0] + E[1][1] + E[2][2]);
  *a1 = E[0][0] * E[1][1] + E[0][0] * E[2][2] + E[1][1] * E[2][2]
    - E[0][1] * E[1][0] - E[0][2] * E[2][0] - E[1][2] * E[2][1];
  *a0 = -(E[0][0] * (E[1][1] * E[2][2] - E[1][2] * E[2][1])
	  + E[0][1] * (E[1][2] * E[2][0] - E[1][0] * E[2][2])
	  + E[0][2] * (E[1][0] * E[2][1] - E[1][1] * E[2][0]));
}



int
bias_eigenvector_to(dbl *v, dbl *target)
{
  dbl dot_prod = 0.0;
  dbl v_mag = 0.0;
  int i, retval;

  retval = 1;
  dot_prod = v_mag = 0.0;
  for(i = 0; i < DIM; i++)
    {
      dot_prod += v[i] * target[i];
      v_mag += v[i] * v[i];
    }
  if(fabs(v_mag) < QTENSOR_SMALL_DBL)
    return retval;
  if(fabs(dot_prod) < BAD_MAGIC_VECTOR_RELATIVE_TOLERANCE)
    {
/*       printf("WARNING: Ni bu hao!  You make bad tensor magic! (dot_prod = %g, ei = %d)\n", dot_prod, ei[pg->imtrx]->ielem); */
      retval = 0;
    }
  if(dot_prod < 0.0)
    {
      v[0] *= -1.0;
      v[1] *= -1.0;
      v[2] *= -1.0;
    }
  return retval;
}

/*****************************************************************************/
/* END of file mm_qtensor_model.c */
/*****************************************************************************/
