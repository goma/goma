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
 *$Id: mm_fill_fill.c,v 5.5 2009-05-19 23:07:03 hkmoffa Exp $
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "std.h"		/* This needs to be here. */
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_bc.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_qtensor_model.h"
#include "mm_unknown_map.h"
#include "rf_node_const.h"

/*
 * Set default memory usage for Aztec's call to y12m higher than
 * the paltry 16 MB that is the default. If you set it on the compile
 * line with a C preprocessor definition, it should override this value,
 * i.e., "-DAZ_MAX_MEMORY_SIZE=1099511627776" for a terabyte. 
 */
#ifndef AZ_MAX_MEMORY_SIZE
#define AZ_MAX_MEMORY_SIZE   536870912	/* 512 MB (2^29) */
#endif


#ifdef PARALLEL
#ifndef MPI
#define MPI			/* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

/* goma include files (of course!) */

#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_shell_util.h"

#define GOMA_MM_FILL_FILL_C

/*
 * This variable is used for more than one linear solver package.
 */


static int neighbor_fill
(Exo_DB *,                /* exo                                       */
       dbl [],                  /* x                                         */
       int ,                    /* current_elem                              */
       int ,                    /* neighbor_elem                             */
       dbl [],                  /* F_neighbor                                */
       int ,                    /* num_local_nodes                           */
       int ,                    /* nodes_per_side                            */
       int [],                  /* local_elem_node_id                        */
       int ,                    /* ielem_type                                */
       int ,                    /* ielem_type_fill                           */
       int [],                  /* node_to_fill                              */
       dbl [][DIM]);           /* x_n                                       */

static int neighbor_species
(Exo_DB *,                /* exo - ptr to exodus ii mesh information   */
       dbl [],                  /* x                                         */
       int ,                    /* current_elem                              */
       int ,                    /* neighbor_elem                             */
       dbl [][MAX_CONC],        /* conc_neighbor                             */
       int ,                    /* num_local_nodes                           */
       int ,                    /* nodes_per_side                            */
       int [],                  /* local_elem_node_id                        */
       int ,                    /* ielem_type                                */
       int ,                    /* ielem_type_fill                           */
       dbl **);                /* x_n                                       */



/*
 * Prototype declarations of static functions.
 */


/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
*
*						-All routines in this
*				NOTE:		 file are called only by
*						 mm_fill.c: matrix_fill
*						 except possibly for static
*						 functions.
*
*       NAME			TYPE			CALLED BY
*  -----------------------------------------------------------------
*  assemble_fill         	void
*
*
******************************************************************************/

#ifdef COUPLED_FILL
int assemble_fill(double tt,
                  double dt,
                  const PG_DATA *pg_data,
                  const int applied_eqn,
                  double xi[3],
                  Exo_DB *exo,
                  double time,
                  struct LS_Mass_Lumped_Penalty *mass_lumped_penalty) {
/******************************************************************************
 *
 * assemble_fill -- integrate fill equation.  This routine assembles
 *                  the fill equations with the other physics equations.
 *
 * in:
 *      tt -- Time integration parameter.
 *      dt -- Current time step size.
 *
 * Created: Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised: 9/24/94 by RRR
 * Revised: 9/24/01 by PKN
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 ******************************************************************************/
  int eqn, var, peqn, pvar, dim, status;
  int i, j, a, b, c;

  dbl F_dot;				/* Fill derivative wrt time. */
  dbl *grad_F;  			/* Fill gradient. */
  dbl grad_II_F[DIM];                   /* Fill surface gradient */
  dbl d_grad_II_F_dmesh[DIM][DIM][MDE];

  dbl *v;                               /* Local velocity. */
  dbl *v_old;                           /* Old v[]. */

  dbl *xx;                              /* Nodal coordinates. */
  dbl *x_old;                           /* Old xx[]. */

  dbl x_dot[DIM];                       /* Time derivative of the mesh displacements. */
  /*dbl x_dot_old[DIM];*/               /* Old x_dot[]. */
  dbl *x_dot_old;                       /* Old x_dot[]. */

  dbl v_rel[DIM];                       /* Velocity relative to the mesh. */
  dbl v_rel_old[DIM];                   /* Old v_rel[]. */

  dbl zero[3] = {0.0, 0.0, 0.0};        /* An array of zeros, for convienience. */

  dbl phi_i;                            /* i-th basis function for the FILL equation. */
  dbl *grad_phi_i;                  /* Gradient of phi_i. */
  dbl grad_II_phi_i[DIM];
  dbl d_grad_II_phi_i_dmesh[DIM][DIM][MDE];

  dbl phi_j;                    /* j-th basis function of a field variable. */
  dbl *grad_phi_j;              /* Gradient of phi_j. */
  dbl grad_II_phi_j[DIM];
  dbl h3;                       /* Volume element (scale factors). */
  dbl det_J;                    /* Determinant of the Jacoabian of transformation. */
  dbl wt;                       /* Gauss point weight. */
  dbl rmp[MDE];                 /* Hold on to the integrands from the residuals. */
  dbl wfcn;                     /* The weight function. */
  dbl d_wfcn_du;                /* Deriv. of wfcn w.r.t. the fluid velocity. */
  dbl d_wfcn_dx;                /* Deriv. of wfcn w.r.t. the mesh displacement. */
  dbl mass = 0.0, advection = 0.0;
  dbl v_dot_Dphi[MDE];          /* v.grad(phi) */
  dbl vc_dot_Dphi[MDE];          /* vcent.grad(phi) */
  dbl v_dot_DF;                 /* v.grad(F) */
  dbl dtinv;                    /* = 1 / dt */
  int Fill_Weight_Fcn;          /* Fill weight function. */

  /* See get_supg_stuff() for a better description of these variables */
  dbl supg_term;                     /* Major term for SUPG -- see get_supg_stuff(). */
  dbl vcent[DIM];                    /* Element centroid velocity (get_supg_stuff(). */
  dbl d_vcent_du[DIM][MDE][DIM];     /* deriv. of vcent[] w.r.t. nodal velocities.   */
  dbl d_supg_term_du[MDE][DIM];      /* deriv. of supg_term w.r.t. nodal velocities. */
  dbl d_supg_term_dx[MDE][DIM];      /* deriv. of supg_term w.r.t. mesh coords.      */
  dbl d_vrel_d_x_rs[DIM][DIM][MDE];  /* deriv of solid rel velo w.r.t. real-solid displ */

  double vmag_old, tau_gls;
  double h_elem;

  /* Alternative newer SUPG style */
  SUPG_terms supg_terms;


  status = 0;
  eqn	 = applied_eqn;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim		  = pd->Num_Dim;                  /* Number of dimensions. */
  Fill_Weight_Fcn = tran->Fill_Weight_Fcn;        /* Which weight function to use */
  wt		  = fv->wt;                       /* Gauss point weight. */
  h3		  = fv->h3;                       /* Differential volume element. */
  det_J		  = bf[eqn]->detJ;                /* Really, ought to be mesh eqn. */
  dtinv           = 1.0 / dt;                     /* Ah, 1 / dt. */

  if(eqn == R_FILL)
	{
  	if (pd->TimeIntegration != STEADY)
    		{ F_dot = fv_dot->F; }
  	else
    		{ F_dot = 0.0; }
  	grad_F = fv->grad_F;
	}
  else
	{
  	if (pd->TimeIntegration != STEADY)
    		{ F_dot = fv_dot->pF[ls->var-PHASE1]; }
  	else
    		{ F_dot = 0.0; }
  	grad_F = fv->grad_pF[ls->var-PHASE1];
        Fill_Weight_Fcn = FILL_WEIGHT_EXPLICIT;
	}

  h_elem = 0.;
  for ( a=0; a<dim; a++) h_elem += pg_data->hsquared[a];
  /* This is the size of the element */
  h_elem = sqrt(h_elem/ ((double )dim));

  /* On 1/11/01, MMH changed these loops from dim's to VIM's.  This
   * was necessary for the PROJECTED_CARTESIAN coordinate system and
   * 3d stability of 2d flow.  The only coordinate system this could
   * possibly cause problems for is the CYLINDRICAL one.  In
   * CYLINDRICAL coordinates, VIM = 3, but there are not always 3
   * components to vectors (b/c the theta-velocity is assumed to be
   * zero).  It passed the test suite, though.
   */

  /*
   * Calculate lubrication velocity for direct integration
   */
  int lubon = 0;
  if ( pd->e[pg->imtrx][R_LUBP] ) {
     if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G  ||
         tran->Fill_Weight_Fcn == FILL_WEIGHT_TG ) {
// if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G  ) {
      lubon = 1;
 } else {
      WH(-1,"\n Multiphase lubrication should be run with Galerkin weighting to \n take advantage of direct velocity calculations.  \n Talk to SAR.");
      lubon = 0;
    }
  }
  int *n_dof = NULL;
  int dof_map[MDE];
  if ( pd->e[pg->imtrx][R_LUBP] ) {

    /* Set up shells */
    n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Calculate velocity */
    calculate_lub_q_v(R_LUBP, time, dt, xi, exo);
    calculate_lub_q_v_old (R_LUBP, tran->time_value_old, tran->delta_t_old, xi, exo );

    /* Set up weights */
    wt    = fv->wt;
    h3    = fv->h3;
    det_J = fv->sdet;

    Inn(grad_F, grad_II_F);
 }
 if ( pd->e[pg->imtrx][R_LUBP_2] ) {

   EH(-1," if you have a fill equation turned on in the R_LUBP_2 phase, you are in the wrong place");
    /* Set up shells */
    n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Calculate velocity */
    calculate_lub_q_v(R_LUBP_2, time, dt, xi, exo);
    calculate_lub_q_v_old (R_LUBP_2, tran->time_value_old, tran->delta_t_old, xi, exo );

    /* Set up weights */
    wt    = fv->wt;
    h3    = fv->h3;
    det_J = fv->sdet;

    Inn(grad_F, grad_II_F);
  }



  /* Use pointers unless we need to do algebra. */
  if (lubon) {
    v     = LubAux->v_avg;
    v_old = LubAux_old->v_avg;
    xx    = fv->x;
    x_old = fv_old->x;
  } else {
    v     = fv->v;
    v_old = fv_old->v;
    xx    = fv->x;
    x_old = fv_old->x;
  }
  if(eqn == R_PHASE1)
	{
         memset(v, 0, sizeof(double)*VIM);
         memset(v_old, 0, sizeof(double)*VIM);
	}

  v_dot_DF = 0.0;
  if ( pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
    {
      x_dot_old = fv_dot_old->x;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	  if (lubon) x_dot[a] = (1+2*tt)/dt * (xx[a] - x_old[a]);
	  v_rel[a]     = v[a] - x_dot[a];
	  v_rel_old[a] = v_old[a] - x_dot_old[a];
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_F[a];       /* v.gradII(F) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_F[a];       /* v.grad(F) */
            }
	}
      if (lubon) ShellRotate(grad_F, fv->d_grad_F_dmesh, grad_II_F, d_grad_II_F_dmesh, n_dof[MESH_DISPLACEMENT1]);
    }
  else if  (pd->TimeIntegration != STEADY &&
	    pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)] &&
	    pd->MeshMotion == TOTAL_ALE)    /*This is the Eulerian solid-mech case */
    {
      xx	 = fv->d_rs;
      x_old	 = fv_old->d_rs;
      x_dot_old =  fv_dot_old->d_rs;
      for(a=0; a < VIM; a++)
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	}
      for(a=0; a < VIM; a++)
	{
	  v_rel[a] = x_dot[a];
	  v_rel_old[a] = x_dot_old[a];
	  for(b=0; b < VIM; b++)
	    {
	      v_rel[a] -= x_dot[b]*fv->grad_d_rs[b][a];
	      v_rel_old[a] -= x_dot_old[b]*fv_old->grad_d_rs[b][a];
	    }
	}
      for ( a=0; a < VIM; a++ )
	{
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_F[a];       /* v.gradII(F) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_F[a];       /* v.grad(F) */
            }
	}

    }
  else
    {
      x_dot_old = zero;
      for ( a=0; a < VIM; a++ )
	{
          x_dot[a]     = 0.0;
	  v_rel[a]     = v[a];
	  v_rel_old[a] = v_old[a];
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_F[a];       /* v.gradII(F) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_F[a];       /* v.grad(F) */
            }
	}
    }

  /* Get the SUPG stuff, if necessary. */
  if ( Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
    {
      memset(vcent,          0, sizeof(double)*DIM);
      memset(d_vcent_du,     0, sizeof(double)*DIM*MDE*DIM);
      memset(d_supg_term_du, 0, sizeof(double)*MDE*DIM);
      memset(d_supg_term_dx, 0, sizeof(double)*MDE*DIM);
      memset(vc_dot_Dphi, 0, sizeof(double)*MDE);
      supg_term = 0.;
      get_supg_stuff(&supg_term, vcent, d_vcent_du, d_supg_term_du, d_supg_term_dx,
		     pd->e[pg->imtrx][R_MESH1]);
    }

  if (Fill_Weight_Fcn == FILL_WEIGHT_SUPG_SHAKIB ||
      Fill_Weight_Fcn == FILL_WEIGHT_SUPG_GP) {
    supg_tau(&supg_terms, dim, 0, pg_data, dt,
             Fill_Weight_Fcn == FILL_WEIGHT_SUPG_SHAKIB, R_FILL);
  }


  /* Compute and save v.grad(phi) and vcent.grad(phi). */
  memset(v_dot_Dphi, 0, sizeof(double)*MDE);


  for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
      grad_phi_i = bf[eqn]->grad_phi[i];

      if (lubon)
        {
         Inn(grad_phi_i, grad_II_phi_i);

         for ( a=0; a < VIM; a++ )
	    {
             v_dot_Dphi[i] += v_rel[a] * grad_II_phi_i[a];
            }
        }
      else
        {
         for ( a=0; a < VIM; a++ )
	    {
	     v_dot_Dphi[i] += v_rel[a] * grad_phi_i[a];
	     if ( Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
	       {
		 vc_dot_Dphi[i] += vcent[a] * grad_phi_i[a];
	       }
            }
        }
    }

  vmag_old = 0.;
  for ( a=0; a<dim; a++)
    {
      vmag_old += fv_old->v[a] * fv_old->v[a];
    }
  vmag_old = sqrt( vmag_old );
  tau_gls = 0.0;
  if (Fill_Weight_Fcn == FILL_WEIGHT_EXPLICIT)
    {
      tau_gls = 1./sqrt((2./dt)*(2./dt) + (2.*vmag_old/h_elem)*(2.*vmag_old/h_elem));
    }

  if (ls != NULL && ls->Semi_Implicit_Integration &&
      (Fill_Weight_Fcn == FILL_WEIGHT_SUPG_SHAKIB || Fill_Weight_Fcn == FILL_WEIGHT_SUPG_GP)) {
    F_dot = (fv->F - fv_old->F) / dt;
    v_dot_DF = 0;
    for (int a = 0; a < VIM; a++) {
      v_dot_DF += (0.5 * (fv->v[a] + fv_old->v[a])) * (0.5 * (fv->grad_F[a] + fv_old->grad_F[a]));
    }

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      v_dot_Dphi[i] = 0;
      for (int a = 0; a < VIM; a++) {
        v_dot_Dphi[i] += (0.5 * (fv->v[a] + fv_old->v[a])) * bf[eqn]->grad_phi[i][a];
      }
    }
  } else if (ls != NULL && !ls->Semi_Implicit_Integration &&
             (Fill_Weight_Fcn == FILL_WEIGHT_SUPG_SHAKIB ||
              Fill_Weight_Fcn == FILL_WEIGHT_SUPG_GP)) {
    F_dot = fv_dot->F;
    v_dot_DF = 0;
    for (int a = 0; a < VIM; a++) {
      v_dot_DF += fv->v[a] * fv->grad_F[a];
    }

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      v_dot_Dphi[i] = 0;
      for (int a = 0; a < VIM; a++) {
        v_dot_Dphi[i] += fv->v[a] * bf[eqn]->grad_phi[i][a];
      }
    }
  } else if (ls != NULL && ls->Semi_Implicit_Integration) {
    EH(-1, "Error Level Set Semi-Implicit Time Integration can only be used "
           "with SUPG_GP and SUPG_SHAKIB");
    return -1;
  }

  double p_e = 0;
  double d_p_e[MDE];
  double d_inv_F_norm[MDE];
  double inv_F_norm = 0;
  if (ls != NULL && ls->Toure_Penalty) {
    for (int i = 0; i < ei[pg->imtrx]->dof[R_FILL]; i++) {
      p_e += mass_lumped_penalty->penalty[i] * bf[R_FILL]->phi[i];
      if (ls->Semi_Implicit_Integration) {
        p_e += mass_lumped_penalty->penalty_old[i] * bf[R_FILL]->phi[i];
      }
    }
    if (ls->Semi_Implicit_Integration) {
      p_e *= 0.5;
    }

    for (int k = 0; k < ei[pg->imtrx]->dof[R_FILL]; k++) {
      d_p_e[k] = 0;
      for (int i = 0; i < ei[pg->imtrx]->dof[R_FILL]; i++) {
        d_p_e[k] += mass_lumped_penalty->d_penalty[i][k] * bf[R_FILL]->phi[i];
      }
      if (ls->Semi_Implicit_Integration) {
        d_p_e[k] *= 0.5;
      }
    }

    double F_norm = 0;
    for (int i = 0; i < dim; i++) {
      F_norm += fv->grad_F[i] * fv->grad_F[i];
    }
    F_norm = sqrt(F_norm);
    inv_F_norm = 1 / (F_norm + 1e-32);
    for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
      d_inv_F_norm[j] = 0;
      for (int i = 0; i < dim; i++) {
        d_inv_F_norm[j] += 2 * fv->grad_F[i] * bf[eqn]->grad_phi[j][i];
      }
      d_inv_F_norm[j] *= -0.5 * inv_F_norm * inv_F_norm * inv_F_norm;
    }
  }

  dbl k_dc = 0;
  // dbl d_k_dc[MDE] = {0};
  if (ls != NULL && ls->YZbeta != YZBETA_NONE) {
    dbl strong_residual = 0;
    strong_residual = fv_dot_old->F;
    for (int p = 0; p < VIM; p++) {
      strong_residual += fv->v[p] * fv_old->grad_F[p];
    }
    // strong_residual -= s_terms.MassSource[w];
    dbl h_elem = 0;
    for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      h_elem += pg_data->hsquared[a];
    }
    /* This is the size of the element */
    h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));

    dbl inner = 0;
    for (int i = 0; i < dim; i++) {
      inner += fv_old->grad_F[i] * fv_old->grad_F[i];
    }

    dbl yzbeta = 0;

    dbl inv_sqrt_inner = (1 / sqrt(inner + 1e-12));
    dbl dc1 = fabs(strong_residual) * inv_sqrt_inner * h_elem * 0.5;
    dbl dc2 = fabs(strong_residual) * h_elem * h_elem * 0.25;

    // dc1 = fmin(supg_terms.supg_tau,dc1);//0.5*(dc1 + dc2);
    yzbeta = fmin(supg_terms.supg_tau, 0.5 * (dc1 + dc2)); // 0.5*(dc1 + dc2);
    //    for (int k = 0; k <  ei[pg->imtrx]->dof[eqn]; k++) {
    //      d_k_dc[k] = 0;
    //    }

    k_dc = yzbeta;
  }

  /**********************************************************************
   **********************************************************************
   ** Residuals
   **********************************************************************
   **********************************************************************/

  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      if(ls == NULL) var  = FILL;
      	else var  = ls->var;
      for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{
	  phi_i = bf[eqn]->phi[i];

	  /************************************************************
	   * Assemble according to the weight function selected.
	   ************************************************************/
	  switch(Fill_Weight_Fcn)
	    {
	    case FILL_WEIGHT_TG:  /* Taylor-Galerkin */

	       mass = F_dot  * phi_i ;

	       advection = 0.;

               if (lubon)
                 {
	          for (a = 0; a < dim; a++)
		     {
		      advection += phi_i * 0.5 * (v_rel[a] + v_rel_old[a]) * grad_II_F[a];
		     }
                 }
               else
                 {
	          for (a = 0; a < dim; a++)
		     {
		      advection += phi_i * 0.5 * (v_rel[a] + v_rel_old[a]) * grad_F[a];
		     }
                 }
	       advection += v_dot_Dphi[i] * v_dot_DF * dt * 0.5;


	      break;
#if 0
	    case FILL_WEIGHT_EXPLICIT:

	       mass = F_dot  * phi_i ;

	       advection = 0.;
	       for (a = 0; a < dim; a++) advection += phi_i * v_rel_old[a] * fv_old->grad_F[a];

	      break;
#endif
#if 0
	    case FILL_WEIGHT_EXPLICIT:

	       wfcn = phi_i;
               for (a = 0; a < dim; a++) wfcn += tau_gls * fv_old->v[a] * bf[eqn]->grad_phi[i][a];

               mass = F_dot  * wfcn ;

	       advection = 0.;
	       for (a = 0; a < dim; a++) advection += wfcn * v_rel_old[a] * fv_old->grad_F[a];

	      break;
#endif
#if 1
	    case FILL_WEIGHT_EXPLICIT:

	       wfcn = phi_i;
               for (a = 0; a < dim; a++) wfcn += tau_gls * fv_old->v[a] * bf[eqn]->grad_phi[i][a];

               mass = F_dot  * wfcn ;

	       advection = 0.;
	       for (a = 0; a < dim; a++) advection += wfcn * v_rel_old[a] * fv->grad_F[a];
	       if (pd->e[pg->imtrx][R_EXT_VELOCITY] &&
			(pfd == NULL || eqn == R_PHASE1)) 
			advection += fv_old->ext_v * wfcn;

	      break;
#endif
	    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

	      mass = F_dot  * phi_i ;   
	      advection = v_dot_DF  * phi_i ;
	      if (pd->e[pg->imtrx][R_EXT_VELOCITY] &&
		  (pfd == NULL || eqn == R_PHASE1)) 
		advection += fv->ext_v * phi_i;

	      break;

	    case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

	      mass = F_dot  * (vc_dot_Dphi[i] * supg_term + phi_i);   
	      advection = v_dot_DF  * (vc_dot_Dphi[i] * supg_term + phi_i);

	      break;

            case FILL_WEIGHT_SUPG_GP:
            case FILL_WEIGHT_SUPG_SHAKIB: {
              dbl wt_func = bf[eqn]->phi[i] + supg_terms.supg_tau * v_dot_Dphi[i];

              mass = F_dot * wt_func;
              advection = v_dot_DF;

              if (ls != NULL && ls->Enable_Div_Term) {
                advection += fv->F * fv->div_v;
              }

              advection *= wt_func;
            } break;

            default:

	      EH(-1,"Unknown Fill_Weight_Fcn");

	    }
	  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

          dbl discontinuity_capturing = 0;
          for (int a = 0; a < dim; a++) {
            discontinuity_capturing += k_dc * fv->grad_F[a] * bf[eqn]->grad_phi[i][a];
          }

          dbl source = 0;
          if (ls != NULL && ls->Toure_Penalty) {
            for (int a = 0; a < dim; a++) {
              source += p_e * fv->grad_F[a] * inv_F_norm * bf[eqn]->grad_phi[i][a];
              //          double tmp = 1.0 / sqrt(fv->grad_F[0] * fv->grad_F[0] +
              //                                  fv->grad_F[1] * fv->grad_F[1]);
              //          source += p_e * fv->grad_F[a] * tmp *
              //          bf[eqn]->grad_phi[i][a];
            }
          }

          /* hang on to the integrand (without the "dV") for use below. */
          rmp[i] = mass + advection + source + discontinuity_capturing;

          lec->R[peqn][i] += (mass + advection + source + discontinuity_capturing) * wt * det_J * h3;
	}
    }

  /**********************************************************************
   **********************************************************************
   * Jacobian terms...
   **********************************************************************
   **********************************************************************/
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{

	  phi_i = bf[eqn]->phi[i];

	  /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
	  grad_phi_i = bf[eqn]->grad_phi[i];

          if (lubon)
            {
             ShellBF(eqn, i, &phi_i, grad_phi_i, grad_II_phi_i,
                     d_grad_II_phi_i_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
            }


	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  /* The weight function for Galerkin and SUPG */
	  wfcn = phi_i;
	  if ( Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
	    {
	      wfcn += vc_dot_Dphi[i] * supg_term;
	    }

	  /* The weight function for SUPG */
	  if ( Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
	    {
	      wfcn = 0.;

	      /* vcent "dot" grad_phi */
	      for(a=0; a<dim; a++)
		wfcn += vcent[a] * bf[eqn]->grad_phi[i][a];

	      wfcn *= supg_term;
	      wfcn += bf[eqn]->phi[i];
	    }
	  else if ( Fill_Weight_Fcn == FILL_WEIGHT_G )
	    {
	      wfcn = bf[eqn]->phi[i];
	    }

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to FILL variable
	   *
	   *************************************************************/

          if(ls == NULL) var  = FILL;
      	        else var  = ls->var;
	  if ( pd->v[pg->imtrx][var] )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  phi_j = bf[eqn]->phi[j];

		  /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		  grad_phi_j = bf[eqn]->grad_phi[j];

                  if (lubon)
                  Inn(grad_phi_j, grad_II_phi_j);

		  /*
		   * Use the selected weight function
		   */
		  switch(Fill_Weight_Fcn)
		    {
		    case FILL_WEIGHT_TG:  /* Taylor-Galerkin */

		      mass = phi_i * phi_j * (1. + 2. * tt) * dtinv;
		      advection = 0.;
		      for ( a=0; a < VIM; a++)
			{
                          if (lubon)
                            {
			     advection += 0.5 * (v_rel[a] + v_rel_old[a]) * grad_II_phi_j[a] * phi_i;
                            }
                          else
                            {
			     advection += 0.5 * (v_rel[a] + v_rel_old[a]) * grad_phi_j[a] * phi_i;
                            }
			}
		      advection += v_dot_Dphi[j] * v_dot_Dphi[i] * dt * 0.5;
		      break;
#if 0
                    case FILL_WEIGHT_EXPLICIT:

		      mass = phi_i * phi_j * (1. + 2. * tt) * dtinv;
		      advection = 0.;
		      break;
#endif
#if 0
                    case FILL_WEIGHT_EXPLICIT:
		      wfcn = phi_i;
                      for (a = 0; a < dim; a++) wfcn += tau_gls * fv_old->v[a] * bf[eqn]->grad_phi[i][a];

		      mass = wfcn * phi_j * (1. + 2. * tt) * dtinv;
		      advection = 0.;
		      break;
#endif
#if 1
                    case FILL_WEIGHT_EXPLICIT:
		      wfcn = phi_i;
                      for (a = 0; a < dim; a++) wfcn += tau_gls * fv_old->v[a] * bf[eqn]->grad_phi[i][a];

		      mass = wfcn * phi_j * (1. + 2. * tt) * dtinv;
		      advection = 0.;
                      for (a = 0; a < dim; a++) advection += wfcn * v_rel_old[a] * bf[eqn]->grad_phi[j][a];
		      break;
#endif
		    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */
		    case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

                      mass = ( phi_j * (1. + 2.*tt) * dtinv ) * wfcn;
		      advection = v_dot_Dphi[j] * wfcn;
		      if (lubon) {
			for (a = 0; a < dim; a++) advection += LubAux->dv_avg_df[a][j] * grad_II_F[a] * wfcn;
		      }

		      break;

                    case FILL_WEIGHT_SUPG_GP:
                    case FILL_WEIGHT_SUPG_SHAKIB: {
                      dbl wt_func = bf[eqn]->phi[i] + supg_terms.supg_tau * v_dot_Dphi[i];
                      if (ls != NULL && ls->Semi_Implicit_Integration) {
                        mass = (phi_j) / dt * wt_func;
                        advection = 0.5 * v_dot_Dphi[j] * wt_func;
                      } else {

                        mass = (phi_j * (1. + 2. * tt) * dtinv) * wt_func;
                        advection = 0;
                        for (int a = 0; a < dim; a++) {
                          advection += fv->v[a] * bf[var]->grad_phi[j][a];
                        }

                        if (ls != NULL && ls->Enable_Div_Term) {
                          advection += bf[var]->phi[j] * fv->div_v;
                        }
                        advection *= wt_func;
                      }
                    } break;

                    default:

		      EH(-1,"Unknown Fill_Weight_Fcn");

		    } /* switch(Fill_Weight_Fcn) */
                    dbl source = 0;
                    if (ls != NULL && ls->Toure_Penalty) {
                      for (int a = 0; a < dim; a++) {
                        source += p_e * fv->grad_F[a] * d_inv_F_norm[j] * bf[eqn]->grad_phi[i][a];
                        source +=
                            p_e * bf[eqn]->grad_phi[j][a] * inv_F_norm * bf[eqn]->grad_phi[i][a];
                        source += d_p_e[j] * fv->grad_F[a] * inv_F_norm * bf[eqn]->grad_phi[i][a];
                        //              source +=
                        //                  fv->grad_F[a] * d_inv_F_norm * bf[eqn]->grad_phi[i][a];
                      }
                    }

                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                    dbl discontinuity_capturing = 0;
                    for (int a = 0; a < dim; a++) {
                      discontinuity_capturing +=
                          k_dc * bf[var]->grad_phi[j][a] * bf[eqn]->grad_phi[i][a];
                    }

                    lec->J[peqn][pvar][i][j] +=
                        (mass + advection + source + discontinuity_capturing) * wt * h3 * det_J;

		} /* for: FILL DoFs */

	    } /* if: FILL exisits */

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to VELOCITY variables
	   *
	   *************************************************************/
	  for ( b=0; b < VIM; b++ )
	    {
	      var = VELOCITY1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{

		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
		    {

		      phi_j = bf[var]->phi[j];

		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_TG:  /* Taylor Galerkin */

			  mass = 0.;
			  advection  = grad_phi_i[b] * v_dot_DF;
			  advection += v_dot_Dphi[i] * grad_F[b];
			  advection *= phi_j * dt * 0.5;
			  advection += 0.5 * phi_j * phi_i * grad_F[b];
			  if (lubon) advection = 0.0;

			  break;

                        case FILL_WEIGHT_EXPLICIT:

			  mass = 0.;
			  advection = 0.;

			  break;

			case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

			  mass = 0.;
			  advection = phi_j * grad_F[b] * phi_i;
			  if (lubon) advection = 0.0;

			  break;

			case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

			  d_wfcn_du = d_supg_term_du[j][b] * vc_dot_Dphi[i];
			  for ( a=0; a < dim; a++ )
			    {
			      d_wfcn_du += supg_term * d_vcent_du[a][j][b] * grad_phi_i[a];
			    }

			  mass = F_dot * d_wfcn_du;
			  advection = v_dot_DF * d_wfcn_du;
			  advection += phi_j * grad_F[b] * wfcn;

			  break;
                        case FILL_WEIGHT_SUPG_GP:
                        case FILL_WEIGHT_SUPG_SHAKIB: {
                          dbl wt_func = bf[eqn]->phi[i];
                          for (int a = 0; a < dim; a++) {
                            wt_func += supg_terms.supg_tau * fv->v[a] * bf[eqn]->grad_phi[i][a];
                          }

                          dbl d_wt_func = 0;
                          for (int a = 0; a < dim; a++) {
                            d_wt_func +=
                                supg_terms.supg_tau * bf[eqn]->phi[j] * bf[eqn]->grad_phi[i][a] +
                                supg_terms.d_supg_tau_dv[b][j] * fv->v[a] * bf[eqn]->grad_phi[i][a];
                          }
                          mass = fv_dot->F * d_wt_func;
                          advection = 0;
                          for (int a = 0; a < dim; a++) {
                            advection += fv->v[a] * bf[eqn]->grad_phi[j][a];
                          }
                          advection *= d_wt_func;
                          advection += bf[var]->phi[j] * bf[eqn]->grad_phi[j][b] * wt_func;
                        } break;

                        default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */

		      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += (mass + advection) * wt * det_J * h3;

		    }
		}
	    }
	      var = EXT_VELOCITY;
	      if ( pd->v[pg->imtrx][var] )
		{

		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
		    {
		      phi_j = bf[var]->phi[j];
		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_TG:  /* Taylor Galerkin */
                        case FILL_WEIGHT_EXPLICIT:

			  mass = 0.;
			  advection = 0.;

			  break;

			case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

			  mass = 0.;
			  advection = 0.;
	       if (pd->e[pg->imtrx][R_EXT_VELOCITY] &&
			(pfd == NULL || eqn == R_PHASE1)) 
			advection += phi_j * phi_i;

			  break;

			case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

			  d_wfcn_du = d_supg_term_du[j][b] * vc_dot_Dphi[i];
			  for ( a=0; a < dim; a++ )
			    {
			      d_wfcn_du += supg_term * d_vcent_du[a][j][b] * grad_phi_i[a];
			    }

			  mass = F_dot * d_wfcn_du;
			  advection = v_dot_DF * d_wfcn_du;
			  advection += phi_j * grad_F[b] * wfcn;

			  break;

			default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */

		      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += (mass + advection) * wt * det_J * h3;

		    }
		}

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to MESH_DISPLACEMENT variables
	   *
	   *************************************************************/
	  for ( b=0; b < VIM; b++ )
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{

		  if ( lubon ) {

		    pvar = upd->vp[pg->imtrx][var];
		    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++) {
		      c = dof_map[j];
		      phi_j = bf[var]->phi[j];

		      switch(Fill_Weight_Fcn)
			{

			case FILL_WEIGHT_TG:  /* Taylor Galerkin */

			  mass = 0.;
		          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

		          advection = 0.0;
		          for ( a=0; a < dim; a++ ) {
                            advection += 0.5 * LubAux->dv_avg_dx[a][b][j] * grad_II_F[a] * phi_i;
			    advection -= 0.5 * (1+2*tt)/dt * phi_j * delta(a,b) * grad_II_F[a] * phi_i;
                            advection += 0.5 * (v_rel[a] + v_rel_old[a]) * d_grad_II_F_dmesh[a][b][j] * phi_i;

                            advection += LubAux->dv_avg_dx[a][b][j] * grad_II_phi_i[a] * v_dot_DF * dt * 0.5;
                            advection -= (1+2*tt)/dt * phi_j * delta(a,b) * grad_II_phi_i[a] * v_dot_DF * dt * 0.5;
                            advection += v_rel[a] * d_grad_II_phi_i_dmesh[a][b][j] * v_dot_DF * dt * 0.5;

                            advection += v_dot_Dphi[i] * LubAux->dv_avg_dx[a][b][j] * grad_II_F[a] * dt * 0.5;
                            advection -= v_dot_Dphi[i] * (1+2*tt)/dt * phi_j * delta(a,b) * grad_II_F[a] * dt * 0.5;
                            advection += v_dot_Dphi[i] * v_rel[a] * d_grad_II_F_dmesh[a][b][j] * dt * 0.5;

		          }
	                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          break;

			case FILL_WEIGHT_G:  /* Plain ol' Galerkin */

			  mass = 0.;
		          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

		          advection = 0.0;
		          for ( a=0; a < dim; a++ ) {
//			    advection += LubAux->dv_avg_dx[a][b][j] * grad_F[a] * phi_i;
//			    advection -= (1+2*tt)/dt * phi_j * delta(a,b) * grad_F[a] * phi_i;
//			    advection += v_rel[a] * fv->d_grad_F_dmesh[a][b][j] * phi_i;

			    advection += LubAux->dv_avg_dx[a][b][j] * grad_II_F[a] * phi_i;
			    advection -= (1+2*tt)/dt * phi_j * delta(a,b) * grad_II_F[a] * phi_i;
			    advection += v_rel[a] * d_grad_II_F_dmesh[a][b][j] * phi_i;
		          }
	                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          break;
                        }
		      lec->J[peqn][pvar][i][c] += (mass + advection) * wt * h3 * fv->sdet;
		      lec->J[peqn][pvar][i][c] += rmp[i] * wt * h3 * fv->dsurfdet_dx[b][c];
		    } /* j loop */

		  } else {

		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {

		      phi_j = bf[var]->phi[j];

		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[var]->grad_phi[j];

		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_TG:  /* Taylor Galerkin */

			  mass = 0.;
			  advection  = -0.5 * phi_j * (1.+2.*tt) * dtinv * phi_i * grad_F[b];
			  advection += -0.5 * dt * phi_j * grad_phi_j[b] * v_dot_DF;

			  break;

                        case FILL_WEIGHT_EXPLICIT:

			  mass = 0.;
			  advection = 0.;

			  break;

			case FILL_WEIGHT_G:  /* Plain ol' Galerkin */

			  mass = 0.;
			  advection = -phi_j * (1. + 2.*tt) * dtinv * grad_F[b];

			  break;

			case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

			  d_wfcn_dx = d_supg_term_dx[j][b] * vc_dot_Dphi[i];
			  for ( a=0; a < dim; a++ )
			    {
			      d_wfcn_dx += supg_term * vcent[a] * bf[eqn]->d_grad_phi_dmesh[i][a][b][j];
			    }

			  mass = F_dot * d_wfcn_dx;
			  advection  = v_dot_DF * d_wfcn_dx ;
			  advection += -phi_j * (1. + 2.*tt) * dtinv * grad_F[b] * wfcn;

			  break;

			default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */

		      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += (mass + advection) * wt * det_J * h3;

		      /* Derivatives of the dV part. */
		      /* rmp[i] holds the integrand without the "dV" part. */

		      /* lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * bf[var]->phi[j]; */
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * phi_j;
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * h3    * bf[eqn]->d_det_J_dm[b][j];

		    } /* for 'j': MESH DoFs */

		  }

		} /* if: MESH exists? */

	    } /* for 'b': MESH componenets */

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to SHELL NORMAL variables
	   *
	   *************************************************************/
	  for ( b=0; b < VIM; b++ )
	    {
	      var = SHELL_NORMAL1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{
		  
		  if ( lubon ) {

		    pvar = upd->vp[pg->imtrx][var];
		    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++) {
		      
		      advection = 0.0;
		      for ( a=0; a < dim; a++ ) {
			advection += LubAux->dv_avg_dnormal[a][b][j] * grad_F[a] * phi_i;
		      }
	              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += advection * wt * h3 * fv->sdet;
		    } /* j loop */

		  } 

		} /* if: SHELL NORMAL exists? */

	    } /* for 'b': SHELL NORMAL components */


	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to REAL_SOLID displacements
	   * for Eulerian solid mechanics
	   *
	   *************************************************************/

	  memset(d_vrel_d_x_rs, 0, sizeof(double)*DIM*DIM*MDE);
	  for (a = 0; a < VIM; a++)
	    {
	      for(b=0; b < VIM; b++)
		{
		  var = SOLID_DISPLACEMENT1 + b;
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      d_vrel_d_x_rs[a][b][j] += phi_j * (1. + 2.*tt) * dtinv *delta(a,b) ;

		      grad_phi_j = bf[var]->grad_phi[j];

		      d_vrel_d_x_rs[a][b][j] -= (1. + 2.*tt) * phi_j * dtinv * fv->grad_d_rs[b][a] ;

		      for(c=0; c < VIM; c++)
			{
			  /*grad_phi_e[i][j][k][l] = e_k e_l : grad(phi_i e_j )*/
			  d_vrel_d_x_rs[a][b][j] -= x_dot[c]* bf[var]->grad_phi_e[j][b][c][a];

			}
		    }
		}
	    }
	  for ( b=0; b < VIM; b++ )
	    {
	      var = SOLID_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{

		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {

		      phi_j = bf[var]->phi[j];

		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[var]->grad_phi[j];

		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_TG:  /* Taylor Galerkin */

			  for(a=0; a < VIM; a++)
			    {
			      advection += phi_i*0.5*d_vrel_d_x_rs[a][b][j]*grad_F[a];
			    }

			  for(a=0; a < VIM; a++)
			    {
			      advection += (v_dot_Dphi[i] * d_vrel_d_x_rs[a][b][j] * grad_F[a]
				            + d_vrel_d_x_rs[a][b][j] * grad_phi_i[a] * v_dot_DF)*dt*0.5;
			    }


			  break;

                        case FILL_WEIGHT_EXPLICIT:

			  advection = 0.;

			  break;

			case FILL_WEIGHT_G:  /* Plain ol' Galerkin */

			  advection  = phi_j * (1. + 2.*tt) * dtinv * grad_F[b]*phi_i;
			  for(a=0; a < VIM; a++)
			    {
			      EH(-1,"This is easy.  Need to copy essentials from TG case above");
			    }

			  break;

			case FILL_WEIGHT_SUPG:  /* Streamline Upwind Petrov Galerkin (SUPG) */

			  EH(-1,"Level set fill for Eulerian solid mech incompatible for FILL_WEIGHT_SUPG");

			  break;

			default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */

		      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += advection * wt * det_J * h3;


		    } /* for 'j': real-solid DoFs */

		} /* if: real-solid exists? */

	    } /* for 'b': real solid componenets */




	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to LUBP variable
	   *
	   *************************************************************/

          var = LUBP;
	  if ( pd->v[pg->imtrx][var] && lubon )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  phi_j = bf[eqn]->phi[j];
		  grad_phi_j = bf[eqn]->grad_phi[j];

		  Inn(grad_phi_j, grad_II_phi_j);

		  mass = 0.0;
		  advection = 0.0;

	          switch(Fill_Weight_Fcn)
                    {

	             case FILL_WEIGHT_TG:  /* Taylor-Galerkin */
		       for (a = 0; a < dim; a++) {
		         advection += 0.5 * LubAux->dv_avg_dp1[a][j] * grad_II_F[a] * phi_i * grad_II_phi_j[a];
		         advection += 0.5 * LubAux->dv_avg_dp2[a][j] * grad_II_F[a] * phi_i * phi_j;

		         advection += 0.5 * LubAux->dv_avg_dp1[a][j] * grad_II_phi_i[a] * grad_II_phi_j[a] * v_dot_DF * dt;
		         advection += 0.5 * LubAux->dv_avg_dp2[a][j] * grad_II_phi_i[a] * phi_j * v_dot_DF * dt;

		         advection += 0.5 * LubAux->dv_avg_dp1[a][j] * grad_II_F[a] * grad_II_phi_j[a] * v_dot_Dphi[i] * dt;
		         advection += 0.5 * LubAux->dv_avg_dp2[a][j] * grad_II_F[a] * phi_j * v_dot_Dphi[i] * dt;
		      }

                       break;

	             case FILL_WEIGHT_G:  /* Plain ol' Galerkin */
		       for (a = 0; a < dim; a++) {
		         advection += LubAux->dv_avg_dp1[a][j] * grad_II_F[a] * wfcn * grad_II_phi_j[a];
		         advection += LubAux->dv_avg_dp2[a][j] * grad_II_F[a] * wfcn * phi_j;
		      }

                       break;

                    }

		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		  lec->J[peqn][pvar][i][j] += (mass + advection) * wt * h3 * det_J;

		} /* for: LUBP DoFs */

	    } /* if: LUBP exisits */

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to SHELL_LUB_CURV variable
	   *
	   *************************************************************/

          var = SHELL_LUB_CURV;
	  if ( pd->v[pg->imtrx][var] && lubon )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
		{
		  phi_j = bf[eqn]->phi[j];
		  grad_phi_j = bf[eqn]->grad_phi[j];

		  mass = 0.0;
		  advection = 0.0;

	          switch(Fill_Weight_Fcn)
                    {

	             case FILL_WEIGHT_TG:  /* Taylor-Galerkin */
		       for (a = 0; a < dim; a++) {

                         advection += 0.5 * LubAux->dv_avg_dk[a][j] * grad_II_F[a] * phi_i * phi_j;
                         advection += 0.5 * LubAux->dv_avg_dk[a][j] * grad_II_phi_i[a] * phi_j * v_dot_DF * dt;
                         advection += 0.5 * LubAux->dv_avg_dk[a][j] * grad_II_F[a] * phi_j * v_dot_Dphi[i] * dt;
		      }

                       break;

	             case FILL_WEIGHT_G:  /* Plain ol' Galerkin */
		       for (a = 0; a < dim; a++)
		           advection += LubAux->dv_avg_dk[a][j] * grad_II_F[a] * phi_i * phi_j;

                       break;
                    }

		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		  lec->J[peqn][pvar][i][j] += (mass + advection) * wt * h3 * det_J;

		} /* for: SHELL_LUB_CURV DoFs */

	    } /* if: SHELL_LUB_CURV exisits */


	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to SHELL_DELTAH variable
	   *
	   *************************************************************/

          var = SHELL_DELTAH;
	  if ( pd->v[pg->imtrx][var] && lubon )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  phi_j = bf[eqn]->phi[j];
		  grad_phi_j = bf[eqn]->grad_phi[j];

		  mass = 0.0;
		  advection = 0.0;
		  for (a = 0; a < dim; a++)
		    advection += LubAux->dv_avg_ddh[a][j] * grad_F[a] * phi_i * phi_j;

		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		  lec->J[peqn][pvar][i][j] += (mass + advection) * wt * h3 * det_J;

		} /* for: SHELL_DELTAH DoFs */

	    } /* if: SHELL_DELTAH exisits */

	} /* for 'i': FILL DoFs */

    } /* if: af->Assemble_Jacobian */

  
  /* clean-up */
  fv->wt = wt;  /* load_neighbor_var_data screws this up */
  safe_free((void *) n_dof);
  
  return(status);

} /* end of assemble_fill */

int 
assemble_fill_ext_v(double tt, 
		    double dt, 
		    dbl hsquared[DIM], 
		    dbl hh[DIM][DIM], 
		    dbl dh_dxnode[DIM][MDE] )
{
/******************************************************************************
 *
 * assemble_fill_ext_v -- integrate fill equation when it uses the extension
 * velocity instead of the fluid velocity.  This routine assembles
 * the fill equation coupled with the other physics equations.
 *
 * in:
 *      tt -- Time integration parameter.
 *      dt -- Current time step size.
 *      hsquared[DIM] element length size square in the three direction
 *      hh[DIM][DIM]
 *      dh_dxnode derivative of hsquared wrt mesh equation
 *
 * Created: Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised: 9/24/94 by RRR
 * Revised: 9/24/01 by PKN
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 ******************************************************************************/
  int eqn, var, peqn, pvar, dim, status;
  int i, j, a, b;
  
  dbl F_dot;				/* Fill derivative wrt time. */
  dbl *grad_F;  			/* Fill gradient. */

  dbl *xx;                              /* Nodal coordinates. */
  dbl *x_old;                           /* Old xx[]. */

  dbl x_dot[DIM];                       /* Time derivative of the mesh displacements. */
  /*dbl x_dot_old[DIM];*/               /* Old x_dot[]. */
  dbl *x_dot_old;                       /* Old x_dot[]. */

  dbl v_rel_n;                       /* Normal velocity relative to the mesh. */
  dbl v_rel_n_old;                   /* Old v_rel[]. */

  dbl zero[3] = {0.0, 0.0, 0.0};        /* An array of zeros, for convienience. */

  dbl phi_i;                            /* i-th basis function for the FILL equation. */

  dbl phi_j;                    /* j-th basis function of a field variable. */
  dbl *grad_phi_i;          /* Gradient of phi_i. */
  dbl *grad_phi_j;          /* Gradient of phi_j. */
  dbl h3;                       /* Volume element (scale factors). */
  dbl det_J;                    /* Determinant of the Jacoabian of transformation. */
  dbl wt;                       /* Gauss point weight. */
  dbl tmp;                      /* A temporary variable. */
  dbl mass = 0.0, advection = 0.0, source = 0.0;
  dbl rmp[MDE];                 /* Hold on to the integrands from the residuals. */
  dbl wfcn;                     /* The weight function. */
  dbl d_wfcn_dext_v;            /* Deriv. of wfcn w.r.t. the extension velocity. */
  dbl d_wfcn_dF;                /* Deriv. of wfcn w.r.t. the Level Set F. */
  dbl v_dot_Dphi[MDE];          /* v.grad(phi) */
  dbl v_dot_DF;                 /* v.grad(F) */
  dbl v_dot_n;                  /* v.n */
  dbl v_old_dot_DF, v_dot_n_old;
  dbl dtinv;                    /* = 1 / dt */
  int Fill_Weight_Fcn;          /* Fill weight function. */

  /* Terms needed for GLS stabilization */

  dbl h_elem, tau_gls;
  int dofs;
  dbl ext_v, ext_v_old, ext_v_mag, ext_v_avg;
  dbl d_tau_ext_v[MDE];
  dbl gfmag_old;

  /* terms for shock capturing term */
  dbl d_visc_dF[MDE];
  dbl d_visc_dext_v[MDE];
  dbl visc = 0.0, visc_coeff;
  dbl gradF_gradphi[MDE];       /* gradF.gradphi */
  dbl gradphi_gradphi[MDE][MDE];/* gradphi.gradphi */

  dbl d_n_dot_dF_dF[MDE][MDE];
  dbl *n_dot_dF = 0;

  dbl x_dot_n; /* scalar dot product */
  dbl x_dot_n_old; /* scalar dot product */

  status = 0;
  eqn	 = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim		  = pd->Num_Dim;                  /* Number of dimensions. */
  Fill_Weight_Fcn = tran->Fill_Weight_Fcn;        /* Which weight function to use */
  wt		  = fv->wt;                       /* Gauss point weight. */
  h3		  = fv->h3;                       /* Differential volume element. */
  det_J		  = bf[eqn]->detJ;                /* Really, ought to be mesh eqn. */
  dtinv           = 1.0 / dt;                     /* Ah, 1 / dt. */

  if(pfd == NULL)
	{
  	if (pd->TimeIntegration != STEADY)
    		{ F_dot = fv_dot->F; }
  	else
    		{ F_dot = 0.0; }
  	grad_F = fv->grad_F;
	}
  else
	{

  	if (pd->TimeIntegration != STEADY)
    		{ F_dot = fv_dot->pF[ls->var-PHASE1]; }
  	else
    		{ F_dot = 0.0; }
  	grad_F = fv->grad_pF[ls->var-PHASE1];
	}

  h_elem = 0.;
  for ( a=0; a<dim; a++) h_elem += hsquared[a];
  /* This is the size of the element */
  h_elem = sqrt(h_elem/ ((double )dim));


  /* On 1/11/01, MMH changed these loops from dim's to VIM's.  This
   * was necessary for the PROJECTED_CARTESIAN coordinate system and
   * 3d stability of 2d flow.  The only coordinate system this could
   * possibly cause problems for is the CYLINDRICAL one.  In
   * CYLINDRICAL coordinates, VIM = 3, but there are not always 3
   * components to vectors (b/c the theta-velocity is assumed to be
   * zero).  It passed the test suite, though. 
   */

  /* Use pointers unless we need to do algebra. */
  xx	 = fv->x;
  x_old	 = fv_old->x;

  load_lsi( 0. );
  load_lsi_derivs();

  dofs     = ei[pg->imtrx]->dof[EXT_VELOCITY];

  /* extension velocity magnitude for stabilization term */  
  ext_v_avg = 0.;
  for (i = 0; i < dofs; i++) ext_v_avg += *esp->ext_v[i] * *esp->ext_v[i];
  ext_v_mag = sqrt(ext_v_avg/dofs);

      tau_gls = 0.5 * h_elem;
      for (i = 0; i < dofs; i++)
        {
          d_tau_ext_v[i] = 0.;
        }

  tau_gls = 1./sqrt((2./dt)*(2./dt) + (2.*ext_v_mag/h_elem)*(2.*ext_v_mag/h_elem));
  for (i = 0; i < dofs; i++)
    {
      d_tau_ext_v[i] = -4. * tau_gls*tau_gls*tau_gls / (h_elem*h_elem) * *esp->ext_v[i] / dofs;
    }

  v_dot_DF = 0.0;
  ext_v = fv->ext_v;
  ext_v_old = fv_old->ext_v;
  gfmag_old = 0.;
  for ( a=0; a<dim; a++)
    {
      gfmag_old += fv_old->grad_F[a] * fv_old->grad_F[a];
    }
  gfmag_old = sqrt( gfmag_old );

  if ( pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
    {
      x_dot_old = fv_dot_old->x;
      v_rel_n = 0.0;
      v_rel_n_old = 0.0;
      x_dot_n = 0.0;
      x_dot_n_old = 0.0;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	  x_dot_n += x_dot[a]* lsi->normal[a];
	  x_dot_n_old += x_dot_old[a]* lsi->normal[a];
	}
      v_rel_n      = ext_v - x_dot_n;
      v_rel_n_old  = ext_v_old - x_dot_n_old;
      v_dot_n      = v_rel_n;                   /* v.n */
      v_dot_n_old  = v_rel_n_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* v.n.grad(F)_mag */
      v_old_dot_DF  = v_rel_n_old * lsi->gfmag;
    }
  else
    {
      x_dot_old = zero;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = 0.0;
	}
      v_rel_n      = ext_v ;
      v_rel_n_old  = ext_v_old;
      v_dot_n      = v_rel_n;                   /* v.n */
      v_dot_n_old  = v_rel_n_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* v.n.grad(F)_mag */
      v_old_dot_DF  = v_rel_n_old * lsi->gfmag;
    }

  var = FILL;
  eqn = R_FILL;
  for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      gradF_gradphi[i] = 0.;
      grad_phi_i = bf[var]->grad_phi[i];
      for ( a=0; a < VIM; a++ )
	{
	  gradF_gradphi[i] += fv->grad_F[a] * grad_phi_i[a];
	}
    }
 
  for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      grad_phi_i = bf[var]->grad_phi[i];
      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
	{
	  grad_phi_j = bf[var]->grad_phi[j];

	  gradphi_gradphi[i][j] = 0.;
	  for ( a=0; a < VIM; a++ )
	    {
	      gradphi_gradphi[i][j] += grad_phi_j[a] * grad_phi_i[a];
	    }
	}
    }

  if ( Fill_Weight_Fcn == FILL_WEIGHT_SC )
    {
      double num, den;
      int sign_adv = 1, sign_ext_v = 1;

      /* very bad, hardcoded viscosity coefficient */
      visc_coeff = 0.1;
      
      num = F_dot + ext_v * lsi->gfmag;
      if ( num < 0. )
        {
          sign_adv = -1;
          num *= -1.;
        }
      den = fabs(ext_v) * lsi->gfmag + h_elem;
      if ( ext_v < 0. ) sign_ext_v = -1;

      visc = h_elem * visc_coeff * num / den;

      var = EXT_VELOCITY;
      for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
        {
          phi_j = bf[var]->phi[j];
          d_visc_dext_v[j] = h_elem * visc_coeff / den * phi_j *
                             ( sign_adv * lsi->gfmag -
                               num * sign_ext_v * lsi->gfmag / den );
        }
        
      var = FILL;
      for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
        {
          phi_j = bf[var]->phi[j];
          d_visc_dF[j] = h_elem * visc_coeff / den *
                         ( sign_adv * ( phi_j * (1. + 2.*tt) * dtinv +
                                        ext_v * lsi->d_gfmag_dF[j] ) -
                           num * fabs(ext_v) * lsi->d_gfmag_dF[j] / den ); 
        }
    }

  /* As seen below there is a strange weighting term in the
     GLS and SC forms.
     Specifically there is a term that is d_gfmag_dF[i].
     This is the same as n . grad_phi[i].
     This is used here to simplify the calculation of the
     derivative of this strange weighting term wrt the nodal F
   */
  if ( Fill_Weight_Fcn == FILL_WEIGHT_SC ||
       Fill_Weight_Fcn == FILL_WEIGHT_GLS )
    {
      n_dot_dF = lsi->d_gfmag_dF;
      var = FILL;
      for ( i=0; i < ei[pg->imtrx]->dof[var]; i++ )
	{
          for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
	    {
              d_n_dot_dF_dF[i][j] = 0.;
              for ( a=0; a<dim; a++)
                d_n_dot_dF_dF[i][j] += bf[var]->grad_phi[i][a] *
                                       lsi->d_normal_dF[a][j];
            }
        }
    }
    
  for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
      grad_phi_i = bf[eqn]->grad_phi[i];
      v_dot_Dphi[i] = 0.;
      for ( a=0; a < VIM; a++ )
	{
	  v_dot_Dphi[i]	+= v_dot_n * lsi->normal[a] * grad_phi_i[a];
	}
    }
      
		
  /**********************************************************************
   **********************************************************************
   ** Residuals
   **********************************************************************
   **********************************************************************/

  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{
	  phi_i = bf[eqn]->phi[i];
	  grad_phi_i = bf[eqn]->grad_phi[i];
	  
	  /* The weight function for Galerkin and GLS */
	  wfcn = phi_i;
	  /* The weight function for GLS */
	  if ( Fill_Weight_Fcn == FILL_WEIGHT_GLS )
	    {
	      wfcn += tau_gls*ext_v*n_dot_dF[i];
	      /*wfcn += tau_gls*n_dot_dF[i];*/
	    }
	  
	  /************************************************************
	   * Assemble according to the weight function selected.
	   ************************************************************/
	  switch(Fill_Weight_Fcn)
	    {
	    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

	      mass = F_dot * phi_i;
	      advection = v_dot_n * phi_i;
	      source = 0.;

	      break;
	      
	    case FILL_WEIGHT_TG:  /* Taylor-Galerkin */
	    
	      mass = F_dot * phi_i;
	      advection = phi_i * 0.5 * (v_dot_DF + v_old_dot_DF);
	      advection += v_dot_Dphi[i] * v_dot_DF * dt * 0.5;
	      source = 0.;
	      
	      break;
#if 0
	    case FILL_WEIGHT_EXPLICIT:

	      mass = F_dot * phi_i;
	      advection = v_old_dot_DF * phi_i;
	      source = 0.;

	      break;
#endif
#if 1
	    case FILL_WEIGHT_EXPLICIT:

	      mass = F_dot * phi_i;
	      advection = fv_old->ext_v * phi_i;
	      source = 0.;

	      break;
#endif

	    case FILL_WEIGHT_GLS:  /* Galerkin Least Squares */

	      mass = F_dot * wfcn;
	      advection = v_dot_DF * wfcn;
	      source = 0.;

	      break;

	    case FILL_WEIGHT_SC:  /* Shock Capturing */

	      mass = F_dot * phi_i;
	      advection = v_dot_DF * phi_i + visc * lsi->gfmag * lsi->d_gfmag_dF[i];
	      source = 0.;

	      break;

              
	    default:

	      EH(-1,"Unknown Fill_Weight_Fcn");

	    }

	  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
	  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		  
	  /* hang on to the integrand (without the "dV") for use below. */
	  rmp[i] = mass + advection;

	  lec->R[peqn][i] += (mass + advection + source) * wt * det_J * h3;
	  
	}
    }

  /**********************************************************************
   **********************************************************************
   * Jacobian terms...
   **********************************************************************
   **********************************************************************/
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{
	  
	  phi_i = bf[eqn]->phi[i];
	  
	  /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
          grad_phi_i = bf[eqn]->grad_phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  /* The weight function for Galerkin and GLS */
	  wfcn = phi_i;
	  /* The weight function for GLS */
	  if ( Fill_Weight_Fcn == FILL_WEIGHT_GLS )
	    {
	      wfcn += tau_gls*ext_v*n_dot_dF[i];
	      /*wfcn += tau_gls*n_dot_dF[i];*/
	    }
	  
	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to FILL variable
	   *
	   *************************************************************/
	  var = ls->var;
	  if ( pd->v[pg->imtrx][var] )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
		{
		  phi_j = bf[var]->phi[j];

		  /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		  grad_phi_j = bf[var]->grad_phi[j];
		  
		  /*
		   * Use the selected weight function
		   */
		  switch(Fill_Weight_Fcn)
		    {
		    case FILL_WEIGHT_GLS:  /* Galerkin Least Squares*/
                    
                      d_wfcn_dF = tau_gls*ext_v*d_n_dot_dF_dF[i][j];
		      /*d_wfcn_dF = tau_gls*d_n_dot_dF_dF[i][j];*/
                      
		      mass = wfcn * (phi_j * (1. + 2.*tt) * dtinv) +
		             d_wfcn_dF * F_dot;
		      advection = wfcn*( v_dot_n * lsi->d_gfmag_dF[j] ) +
                                  d_wfcn_dF*( v_dot_DF );
		      source = 0.;

		      break;

		    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      advection = 0.;
		      source = 0.;

		      break;
                      
		    case FILL_WEIGHT_SC:  /* Shock Capturing */

		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      advection = wfcn*( v_rel_n * lsi->d_gfmag_dF[j] ) +
                            lsi->d_gfmag_dF[i] * ( d_visc_dF[j] * lsi->gfmag +
                                                   visc * lsi->d_gfmag_dF[j] ) +
                            d_n_dot_dF_dF[i][j] * ( visc * lsi->gfmag );
		      source = 0.;
		      break;
		      
		    case FILL_WEIGHT_TG:  /* Taylor-Galerkin */
		      
		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      
		      tmp = 0.;
		      for ( a=0; a < VIM; a++ ) tmp += v_dot_n * lsi->d_normal_dF[a][j] * grad_phi_i[a];
		      
		      advection = phi_i * 0.5 * ( v_dot_n + v_dot_n_old ) * lsi->d_gfmag_dF[j];
		      advection += dt * 0.5 * (v_dot_Dphi[i] * v_dot_n * lsi->d_gfmag_dF[j] + 
		                               v_dot_DF * tmp);
		      source = 0.;
		      break;
#if 0
                    case FILL_WEIGHT_EXPLICIT:
 
		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      advection = 0.;
		      source = 0.;

		      break;
#endif
#if 0
                    case FILL_WEIGHT_EXPLICIT:
 
		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      advection = wfcn * v_dot_n_old * lsi->d_gfmag_dF[j];
		      source = 0.;

		      break;
#endif
#if 1
                    case FILL_WEIGHT_EXPLICIT:
 
		      mass = wfcn*( phi_j * (1. + 2.*tt) * dtinv );
		      advection = 0.;
		      source = 0.;

		      break;
#endif
                      
		    default:

		      EH(-1,"Unknown Fill_Weight_Fcn");

		    } /* switch(Fill_Weight_Fcn) */

		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
		  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		  
		  lec->J[peqn][pvar][i][j] += wt * h3 * det_J * (mass + advection + source);		  

		} /* for: FILL DoFs */

	    } /* if: FILL exisits */
	  
	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to EXT_VELOCITY variables
	   *
	   *************************************************************/
	  var = EXT_VELOCITY;
	  if ( pd->v[pg->imtrx][var] )
	    {
	      
	      pvar = upd->vp[pg->imtrx][var];
	      for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
		{
		  
		  phi_j = bf[var]->phi[j];	      
		  
		  switch(Fill_Weight_Fcn)
		    {
		    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */
		      
		      mass = 0.;
		      advection = phi_j * phi_i;
		      source = 0.;
		      
		      break;
		      
		    case FILL_WEIGHT_GLS:  /* Galerkin Least Squares (GLS) */
		      
                      d_wfcn_dext_v =  d_tau_ext_v[j]*ext_v*n_dot_dF[i] +
			               tau_gls*phi_j*n_dot_dF[i];
		      /*d_wfcn_dext_v =  0.;*/
		     
		      mass = F_dot * d_wfcn_dext_v;
		      advection = v_dot_DF * d_wfcn_dext_v;
		      advection += phi_j * lsi->gfmag * wfcn;
		      source = 0.;
		      
		      break;

 		    case FILL_WEIGHT_SC:  /* Shock Capturing */

		      mass = 0.;
		      advection = phi_j * lsi->gfmag * phi_i +
                                  lsi->gfmag * lsi->d_gfmag_dF[i] * d_visc_dext_v[j];
		      source = 0.;
                            
		      break;
		      
		    case FILL_WEIGHT_TG:  /* Taylor-Galerkin */
		      
		      mass = 0.;

		      tmp = 0.;
		      for ( a=0; a < VIM; a++ ) tmp += lsi->normal[a] * grad_phi_i[a];
		      
		      advection = 0.5 * phi_i * phi_j * lsi->gfmag;
		      advection += dt * 0.5 * ( v_dot_Dphi[i] * lsi->gfmag + tmp * v_dot_DF ) * phi_j;
		      source = 0.;
		      
		      break;
                      
                    case FILL_WEIGHT_EXPLICIT:
		      
		      mass = 0.;
                      advection = 0.;
		      source = 0.;
		      
		      break;
                       
		    default:
		      
		      EH(-1,"Unknown Fill_Weight_Fcn");
		      
		    } /* switch(Fill_Weight_Fcn) */
		  
		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
		  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		  
		  lec->J[peqn][pvar][i][j] += wt * det_J * h3 * (mass + advection + source);
		  
		}
	    }
	  
	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to MESH_DISPLACEMENT variables
	   *
	   * Ugh - haven't touched this yet!
	   *************************************************************/
	  for ( b=0; b < VIM; b++ )
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{
		  
		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      
		      phi_j = bf[var]->phi[j];	      
		      
		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[var]->grad_phi[j];

		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_G:  /* Plain ol' Galerkin */
			  
			  mass = 0.;
			  advection = -phi_j * (1. + 2.*tt) * dtinv * grad_F[b];

			  break;

			case FILL_WEIGHT_GLS:  /* Galerkin Least Squares */
			  /*
			  d_wfcn_dx = d_supg_term_dx[j][b] * vc_dot_Dphi[i];
			  for ( a=0; a < dim; a++ )
			    {
			      d_wfcn_dx += supg_term * vcent[a] * bf[eqn]->d_grad_phi_dmesh[i][a][b][j];
			    }

			  tmp  = ( F_dot + v_dot_DF ) * d_wfcn_dx ;
			  tmp += -phi_j * (1. + 2.*tt) * dtinv * grad_F[b] * wfcn; */
			  /* Need to fix above for GLS */

			  break;

			default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */
		      
		      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		      lec->J[peqn][pvar][i][j] += wt * det_J * h3 * (mass + advection);

		      /* Derivatives of the dV part. */
		      /* rmp[i] holds the integrand without the "dV" part. */

		      /* lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * bf[var]->phi[j]; */
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * phi_j;
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * h3    * bf[eqn]->d_det_J_dm[b][j];

		    } /* for 'j': MESH DoFs */

		} /* if: MESH exists? */

	    } /* for 'b': MESH componenets */

	} /* for 'i': FILL DoFs */

    } /* if: af->Assemble_Jacobian */

  return(status);

} /* end of assemble_fill_ext_v */

int 
assemble_fill_gradf(double tt, 
		    double dt, 
		    dbl hsquared[DIM], 
		    dbl hh[DIM][DIM], 
		    dbl dh_dxnode[DIM][MDE] )
{
/******************************************************************************
 *
 * assemble_fill_gradf -- integrate fill equation when it uses the 
 *         sign(f) grad f * grad f = sign(f)
 * instead of the fluid velocity advection equation.  This routine assembles
 * the fill equation coupled with the other physics equations.
 *
 * in:
 *      tt -- Time integration parameter.
 *      dt -- Current time step size.
 *      hsquared[DIM] element length size square in the three direction
 *      hh[DIM][DIM]
 *      dh_dxnode derivative of hsquared wrt mesh equation
 *
 * Created: Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 ******************************************************************************/
  int eqn, var, peqn, pvar, dim, status;
  int i, j, a, b;
  
  dbl *grad_F;  			/* Fill gradient. */
  dbl phi_i;                    /* i-th basis function for the FILL equation. */

  dbl phi_j;                    /* j-th basis function of a field variable. */
  dbl *grad_phi_j;              /* Gradient of phi_j. */
  dbl *grad_phi_i;              /* Gradient of phi_i. */
  dbl h3;                       /* Volume element (scale factors). */
  dbl det_J;                    /* Determinant of the Jacoabian of transformation. */
  dbl wt;                       /* Gauss point weight. */
  dbl tmp = 0.0;                      /* A temporary variable. */
  dbl rmp[MDE];                 /* Hold on to the integrands from the residuals. */
  dbl wfcn;                     /* The weight function. */
  dbl d_wfcn_dx;                /* Deriv. of wfcn w.r.t. the mesh displacement. */
  dbl rhs = 0.0;                      /* Temp. storage of the residual. */
  dbl grad_F_2;                 /* gradF.grad_f */
  dbl gradF_gradphi[MDE];       /* gradF.gradphi */
  dbl gradphi_gradphi[MDE][MDE];/* gradphi.gradphi */
  dbl d_gradF_dmesh_gradF;      /* d_dgradF_dmesh.gradF */
  dbl S;                        /* sign(F) */
  int Fill_Weight_Fcn;          /* Fill weight function. */

  /* Terms needed for GLS stabilization */
  dbl h_elem;

  /* terms for shock capturing term */
  dbl d_visc_dF[MDE];
  dbl visc, visc_coeff;

  status = 0;
  eqn	 = R_FILL;
  var    = FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim		  = pd->Num_Dim;                  /* Number of dimensions. */
  Fill_Weight_Fcn = tran->Fill_Weight_Fcn;        /* Which weight function to use */
  wt		  = fv->wt;                       /* Gauss point weight. */
  h3		  = fv->h3;                       /* Differential volume element. */
  det_J		  = bf[eqn]->detJ;                /* Really, ought to be mesh eqn. */

  h_elem = 0.;
  for ( a=0; a<dim; a++) h_elem += hsquared[a];
  /* This is the size of the element */
  h_elem = sqrt(h_elem/ ((double )dim));

  /* Use pointers unless we need to do algebra. */
  grad_F = fv->grad_F;

  load_lsi( ls->Length_Scale );
  load_lsi_derivs();
  

  /* Precompute vector dot products */
  grad_F_2 = 0.;
  for ( a=0; a < VIM; a++ )
    {
      grad_F_2 += grad_F[a] *grad_F[a];
    }

  eqn	 = R_FILL;
  for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      gradF_gradphi[i] = 0.;
      grad_phi_i = bf[eqn]->grad_phi[i];
      for ( a=0; a < VIM; a++ )
	{
	  gradF_gradphi[i] += grad_F[a] * grad_phi_i[a];
	}
    }

  var    = FILL;
  for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      grad_phi_i = bf[eqn]->grad_phi[i];
      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
	{
	  grad_phi_j = bf[var]->grad_phi[j];

	  gradphi_gradphi[i][j] = 0.;
	  for ( a=0; a < VIM; a++ )
	    {
	      gradphi_gradphi[i][j] += grad_phi_j[a] * grad_phi_i[a];
	    }
	}
    }
  
  /* sign of LS */
  S = 2. * lsi->H - 1.;
  
  /* artificial viscosity for stabilization */
#if 0
  /* best consistent form so far, but not nearly as successful as inconsistent form below */
  visc_coeff = 0.5;
  if ( ( lsi->gfmag - 1. ) > 0.25 )
    visc = h_elem * visc_coeff * 0.25;
  else if ( ( lsi->gfmag - 1. ) < -0.02 )
    visc = -h_elem * visc_coeff * 0.02;
  else
    visc = h_elem * visc_coeff * ( lsi->gfmag - 1. );
    
  var = FILL;
  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
    {
      if ( ( lsi->gfmag - 1. ) > 0.25 )
        d_visc_dF[j] = 0.;
      else if ( ( lsi->gfmag - 1. ) < -0.02 )
        d_visc_dF[j] = 0.;
      else
        d_visc_dF[j] = h_elem * visc_coeff * lsi->d_gfmag_dF[j];

    }
#endif
#if 0
  /* inconsistent form that works pretty well, except is overly diffusive */
  visc_coeff = 0.5;
  visc = h_elem * visc_coeff;
  
  var = FILL;
  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
    {
      d_visc_dF[j] = 0.;
    }
#endif
#if 0
  /* simple consistent form */
  visc_coeff = 0.5;
  visc = h_elem * visc_coeff * ( 1. - lsi->gfmaginv );
    
  var = FILL;
  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
    {
      d_visc_dF[j] = -h_elem * visc_coeff * lsi->d_gfmaginv_dF[j];
    }
#endif
#if 1
  /* experimental consistent form */
  visc_coeff = 0.5;
  if ( lsi->gfmag < 1. )
    visc = 0.;
  else
    visc = h_elem * visc_coeff * ( 1. - lsi->gfmaginv );
    
  var = FILL;
  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
    {
      if ( lsi->gfmag < 1. )
        d_visc_dF[j] = 0.;
      else
        d_visc_dF[j] = -h_elem * visc_coeff * lsi->d_gfmaginv_dF[j];
    }
#endif





  /**********************************************************************
   **********************************************************************
   ** Residuals
   **********************************************************************
   **********************************************************************/

  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      for( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{
	  phi_i = bf[eqn]->phi[i];
	  grad_phi_i = bf[eqn]->grad_phi[i];

          /*if ( sign_change( F, *esp->F[i] ) ) continue;*/
	  
	  
	  /************************************************************
	   * Assemble according to the weight function selected.
	   ************************************************************/
	  switch(Fill_Weight_Fcn)
	    {
	    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

	      rhs = (lsi->gfmag - 1.) * S * phi_i;
	      rhs *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

	      break;

	    case FILL_WEIGHT_GLS:  /* Galerkin Least Squares */

	      rhs = (lsi->gfmag - 1.) * ( S * phi_i  ) +
	            visc * gradF_gradphi[i];	
              rhs *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
	      break;

              
	    default:

	      EH(-1,"Unknown Fill_Weight_Fcn");

	    }

	  /* hang on to the integrand (without the "dV") for use below. */
	  rmp[i] = rhs;

	  lec->R[peqn][i] += rhs * wt * det_J * h3;
	  
	}
    }

  /**********************************************************************
   **********************************************************************
   * Jacobian terms...
   **********************************************************************
   **********************************************************************/
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
	{
	  
	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  /* The weight function for Galerkin */
	  wfcn = S * phi_i;
	  
	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to FILL variable
	   *
	   *************************************************************/
	  var = FILL;
	  if ( pd->v[pg->imtrx][var] )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
		{
		  phi_j = bf[var]->phi[j];

		  /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		  grad_phi_j = bf[var]->grad_phi[j];
		  
		  /*
		   * Use the selected weight function
		   */
		  switch(Fill_Weight_Fcn)
		    {
		    case FILL_WEIGHT_GLS:  /* Galerkin Least Squares*/
                      tmp = wfcn * lsi->d_gfmag_dF[j] +
			    (d_visc_dF[j] * gradF_gradphi[i] + visc * gradphi_gradphi[i][j]);
                      tmp *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		      break;

		    case FILL_WEIGHT_G:  /* Plain ol' Galerkin method */

		      tmp = wfcn * lsi->d_gfmag_dF[j];
                      tmp *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
		      break;
                      
		    default:

		      EH(-1,"Unknown Fill_Weight_Fcn");

		    } /* switch(Fill_Weight_Fcn) */

		  lec->J[peqn][pvar][i][j] += wt * h3 * det_J * tmp;		  

		} /* for: FILL DoFs */

	    } /* if: FILL exisits */
	  
	  
	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to MESH_DISPLACEMENT variables
	   *
	   * Ugh - haven't touched this yet!
	   *************************************************************/
	  for ( b=0; b < VIM; b++ )
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if ( pd->v[pg->imtrx][var] )
		{
		  
		  pvar = upd->vp[pg->imtrx][var];
		  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
		    {
		      d_gradF_dmesh_gradF = 0.;
		      for ( a=0; a < VIM; a++ )
			{
			  d_gradF_dmesh_gradF += fv->d_grad_F_dmesh[a] [b][j]*grad_F[a];
			}
		      phi_j = bf[var]->phi[j];	      
		      
		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[var]->grad_phi[j];

		      
		      switch(Fill_Weight_Fcn)
			{
			case FILL_WEIGHT_G:  /* Plain ol' Galerkin */
			  
			  tmp = 2.*wfcn * S *d_gradF_dmesh_gradF;
			  tmp *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

			  break;

			case FILL_WEIGHT_GLS:  /* Galerkin Least Squares */

			  /* Ignore derivative  d/dx of h_elem for now */
			  d_wfcn_dx = h_elem*lsi->d_gfmaginv_dmesh[b][j];
			 

			  tmp  =  wfcn*(2.*S*d_gradF_dmesh_gradF ) +
                            d_wfcn_dx*(grad_F_2 - 1.0);
			  tmp *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
			  

			  break;

			default:

			  EH(-1,"Unknown Fill_Weight_Fcn");

			} /* switch(Fill_Weight_Fcn) */
		      
		      /* Add in the dV */
		      tmp *= wt * det_J * h3;

		      lec->J[peqn][pvar][i][j] += tmp;

		      /* Derivatives of the dV part. */
		      /* rmp[i] holds the integrand without the "dV" part. */

		      /* lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * bf[var]->phi[j]; */
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * det_J * fv->dh3dq[b] * phi_j;
		      lec->J[peqn][pvar][i][j] += rmp[i] * wt * h3    * bf[eqn]->d_det_J_dm[b][j];

		    } /* for 'j': MESH DoFs */

		} /* if: MESH exists? */

	    } /* for 'b': MESH componenets */

	} /* for 'i': FILL DoFs */

    } /* if: af->Assemble_Jacobian */

  return(status);

} /* end of assemble_fill_gradf */

#else /* #ifdef COUPLED_FILL */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 * assemble_fill -- integrate explicit fill equation with Taylor-Galerkin
 *                  formulation.  This routine is for the de-coupled method.
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	9/24/94 by RRR
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int 
assemble_fill(double afill[],	/* Jacobian matrix for fill equation  */
	      int ijaf[],	/* pointer to nonzeros in Jacobian matrix   */
	      double rf[],	/* rhs vector   */
	      double dt,	/* current time step size */
	      double tt,	/* parameter to vary time integration from
				 * explicit (tt = 1) to implicit (tt = 0) */
	      int node_to_fill[],
	      double xi[DIM],
	      Exo_DB *exo,
	      double time)
{
  dbl wfcn;                          /* Holds the final weight function (G, TG, or SUPG) */

  /* See get_supg_stuff() for a better description of these variables */
  dbl supg_term;                     /* major term for SUPG -- see get_supg_stuff() */
  dbl vcent[DIM];                    /* Element centroid velocity (get_supg_stuff() */
  dbl d_vcent_du[DIM][MDE][DIM];     /* deriv. of vcent[] w.r.t. nodal velocities */
  dbl d_supg_term_du[MDE][DIM];      /* deriv. of supg_term w.r.t. nodal velocities */
  dbl d_supg_term_dx[MDE][DIM];      /* deriv. of supg_term w.r.t. mesh coords. */

  int eqn, dim, a, i, j, status, I, J, idof,jdof, ie, je, ja;
  int ki, kj, nvdofi, nvdofj, offseti;
  dbl d_grad_phi_i;

  dbl F;				/* Fill. */
  dbl F_dot;				/* Fill derivative wrt time. */

  dbl v[DIM];				/* Velocity field. */
  dbl v_old[DIM];			/* Velocity field at last time step. */
  dbl v_rel[DIM];			/* relative Velocity field wrt mesh motion */
  dbl v_rel_old[DIM];			/* relative Velocity field wrt mesh motion  
					 * at last time step */

  dbl xx[DIM];	        		/* position field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl x_dot_old[DIM];
  dbl x_old[DIM];			/* old position field. */
  dbl grad_F[DIM];			/* gradient of Fill. */
  dbl num_diff;

  /*
   * Galerkin weighting functions for i-th energy residuals 
   * and some of their derivatives...
   */
  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;
  dbl wt;
  double rhs;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  
  /* Get the SUPG stuff, if necessary. */
  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
    {
      memset(vcent, 0, sizeof(double)*DIM);
      memset(d_vcent_du, 0, sizeof(double)*DIM*MDE*DIM);
      memset(d_supg_term_du, 0, sizeof(double)*MDE*DIM);
      memset(d_supg_term_dx, 0, sizeof(double)*MDE*DIM);
      supg_term = 0.0;
      get_supg_stuff(&supg_term, vcent, d_vcent_du, d_supg_term_du, d_supg_term_dx,
		     pd->e[pg->imtrx][R_MESH1]);
    }
  /* static char yo[] = "assemble_fill"; */


  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;
  eqn   = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  wt	= fv->wt;			/* Gauss point weight.           */
  h3	= fv->h3;			/* Differential volume element.  */
  det_J	= bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  F     = fv->F;
  F_dot = fv_dot->F;
  for (a = 0; a < dim; a++) {
    grad_F[a] = fv->grad_F[a];
  }

  /*
   * HKM -> Hmm, at this point in the code, xx[] is equal to x_old[].
   *        Therefore, the time derivatives calculated below are wrong.
   *        I won't try to fix this. However, I suggest that xx[] and
   *        xdot[] should be changed to represent the predicted 
   *        values at t = time2.
   *        Then, these formulas might start to make sense.
   *        A related issue is the Taylor-Galerkin representation.
   *        The treatment below doesn't agree with the Zienkiewicz
   *        book. Also, it doesn't take into account that the velocity
   *        field may be nonsolenodal and time dependent. 
   */
  for ( a=0; a<dim; a++) {
    v[a] = fv->v[a];
    /* note, these are zero for steady calculations */
    v_old[a] = fv_old->v[a];
    if (((pd->TimeIntegration) != STEADY) && 
	(pd->v[pg->imtrx][MESH_DISPLACEMENT1+a])) {
      xx[a] = fv->x[a];
      x_old[a] = fv_old->x[a];
      x_dot_old[a] = fv_dot->x[a];
      x_dot[a] = (1.+2.*tt) * (xx[a] - x_old[a])/dt - 2. * tt * x_dot_old[a];
      v_rel[a] = v[a] - x_dot[a];
      v_rel_old[a] = v_old[a] - x_dot_old[a];
    } else {
      x_dot[a] = 0.0;
      x_dot_old[a] = 0.0;
      v_rel[a] = v[a];
      v_rel_old[a] = v_old[a];
    }
  }

  num_diff = 0.0000;
  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  if ( af->Assemble_Residual )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	{
	  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
	  node = Nodes[I];
	  nv = node->Nodal_Vars_Info[pg->imtrx];
	  /* check for multiple dofs */
	  nvdofi = get_nv_ndofs(nv, eqn);
	  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);
	
	  for (ki = 0; ki < nvdofi; ki++)	
	    {
	      /* set rhs to zero for Dirichlet boundary conditions */
	      if (node->DBC[pg->imtrx] &&
		  (node->DBC[pg->imtrx][offseti + ki] != -1) &&
		  Debug_Flag >= 0)
		{
		  rf[node_to_fill[I] + ki] = 0.0;
		}
	      else 
		{
		  double v_dot_Dphi_i = 0.0, v_dot_DF = 0.0;
  
		  /* check to make sure that unknowns are defined at this node,
		   * otherwise don't add anything to this node */
 		  idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

		  /* also convert from node number to dof number */
		  phi_i = bf[eqn]->phi[idof];
		  
		  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_TG ) 
		    {
		      rhs = F_dot  * phi_i ;   /* try implementing an explicit/implicit method */
		      /*      rhs = (F-F_old)/dt * phi_i ; */
		      
		      for (a = 0; a < dim; a++)
			{
			  v_dot_Dphi_i += v_rel[a]* bf[eqn]->grad_phi[idof] [a];
			  v_dot_DF     += v_rel[a]*grad_F[a];		      
			  rhs	       += phi_i * (v_rel[a] + v_rel_old[a])/2. *  grad_F[a]; 
			}
		      
		      rhs += v_dot_Dphi_i * v_dot_DF * dt / 2.0;  /* Taylor-Galerkin */
		    }
		  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ) 
		    {
		      rhs = F_dot;

		      for (a=0; a < dim; a++)
			{
			  rhs += v_rel[a] * grad_F[a];
			}

		      rhs *= bf[eqn]->phi[idof];
		    }
		  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_EXPLICIT ) 
		    {
		      rhs = F_dot;

		      for (a=0; a < dim; a++)
			{
			  rhs += v_rel_old[a] * grad_F[a];
			}

		      rhs *= bf[eqn]->phi[idof];
		    }
		  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_SUPG ) 
		    {
		      /* Make the SU/PG weight function */
		      /* wfcn = w + p = phi + supg_term * (vcent "dot" grad(phi) ) */
		      
		      /* vcent "dot" grad_phi */
		      wfcn = 0.;
		      for(a=0; a<dim; a++)
			{
			  wfcn += vcent[a] * bf[eqn]->grad_phi[idof][a];
			}
		      
		      wfcn *= supg_term;
		      wfcn += bf[eqn]->phi[idof];
		      
		      rhs = F_dot;
		      for(a=0; a<dim; a++)
			rhs += v_rel[a] * grad_F[a];
		      rhs *= wfcn;
		    }

		  /* This is common to all "tran->Fill_Weight_Fcn"s */
		  rf[node_to_fill[I] + ki] += rhs * wt * det_J * h3;
		  
		}
	    }
	}
    } 

  if ( af->Assemble_Jacobian )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	  {
		  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i]; 
		  
		  node = Nodes[I];
		  nv = node->Nodal_Vars_Info[pg->imtrx];
		  nvdofi = get_nv_ndofs(nv, eqn);
		  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);	  
		  for (ki = 0; ki < nvdofi; ki++)
		  {
			  ie = node_to_fill[I] + ki;
			  idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
			  phi_i = bf[eqn]->phi[idof];
			  
			  /* The weight function for SUPG */
			  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
			  {
				  wfcn = 0.;
				  
				  /* vcent "dot" grad_phi */
				  for(a=0; a<dim; a++)
					  wfcn += vcent[a] * bf[eqn]->grad_phi[idof][a];
				  
				  wfcn *= supg_term;
				  wfcn += bf[eqn]->phi[idof];
			  }
			  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ||
						tran->Fill_Weight_Fcn == FILL_WEIGHT_EXPLICIT )
			  {
				  wfcn = bf[eqn]->phi[idof];
			  }
			  
			  /* derivatives of fill equation wrt to fill variable */	  
			  for( j=0; j< ei[pg->imtrx]->num_local_nodes; j++) 
			  {
				  J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j]; 
				  nvdofj = Dolphin[pg->imtrx][J][eqn];
				  for ( kj=0; kj<nvdofj; kj++)
					  
				  {
					  
					  jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
					  phi_j = bf[eqn]->phi[jdof];
					  
					  je = node_to_fill[J] + kj;
					  ja     = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie+1], ijaf);
					  EH(ja, "Could not find vbl in sparse matrix.");
					  
					  /* set diagonal to one for Dirichlet boundary conditions */
					  if (node->DBC[pg->imtrx] && 
						  ((int) node->DBC[pg->imtrx][offseti + ki] != -1) &&
						  Debug_Flag >= 0)
					  {
						  afill[ja] = delta(i,j); 
					  }
					  else
					  {
						  double v_dot_Dphi_i= 0.0, v_dot_Dphi_j = 0.0;
						  /*  afill[ja] += wt * det_J * phi_i * phi_j/dt*h3; */
						  
						  /*
						   * Use the selected Finite Element Method
						   */
						  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_TG )
						  {
							  
							  afill[ja] += wt * h3 * det_J * phi_i * phi_j * (1. + 2. * tt)/dt;
							  for ( a=0; a<dim; a++)
							  {
								  afill[ja] += 
								  wt * h3 * det_J * phi_i * (v_rel[a] + v_rel_old[a])/2. 
								  * bf[eqn]->grad_phi[jdof] [a] ;/*    + */
								  /* wt * h3 * det_J * bf[eqn]->grad_phi[idof] [a]  */
								  /* * v_rel[a] * v_rel_old[a] * bf[eqn]->grad_phi[jdof] [a] *dt/2. */
								  /* +  num_diff * wt * h3 * det_J * bf[eqn]->grad_phi[idof] [a]  */
								  /* bf[eqn]->grad_phi[jdof] [a] ; */
								  
								  v_dot_Dphi_j += v_rel[a]*bf[eqn]->grad_phi[jdof][a];
								  v_dot_Dphi_i += v_rel[a]* bf[eqn]->grad_phi[idof] [a];
								  
							  }
							  afill[ja] += 	wt * h3 * det_J * v_dot_Dphi_j*v_dot_Dphi_i*dt/2.0;
						  }
						  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G )
						  {
							  /* The "mass" term */
							  afill[ja] += wt*h3*det_J*wfcn*phi_j*(1. + 2.*tt)/dt;
							  
							  /* The advection term */
							  for ( a=0; a<dim; a++)
							  {
								  afill[ja] += wt*h3*det_J*wfcn*v_rel[a]*bf[eqn]->grad_phi[jdof][a];
							  }
							  
						  }
						  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_EXPLICIT )
						  {
							  /* The "mass" term */
							  afill[ja] += wt*h3*det_J*wfcn*phi_j*(1. + 2.*tt)/dt;
							  
							  /* The advection term */
							  for ( a=0; a<dim; a++)
							  {
								  afill[ja] += wt*h3*det_J*wfcn*v_rel_old[a]*bf[eqn]->grad_phi[jdof][a];
							  }
							  
						  }
						  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_SUPG )
						  {
							  /* The "mass" term */
							  afill[ja] += wt*h3*det_J*wfcn*phi_j*(1. + 2.*tt)/dt;
							  
							  /* The advection term */
							  for ( a=0; a<dim; a++)
							  {
								  afill[ja] += wt*h3*det_J*wfcn*v_rel[a]*bf[eqn]->grad_phi[jdof][a];
							  }
							  
						  }
					  }
				  }
			  }
		  }
	  }
    }

  return(status);
} /* end of assemble_fill */

/*
 * assemble_fill_ext_v -- integrate explicit fill equation with Taylor-Galerkin
 *                  formulation.  This routine is for the de-coupled method.
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	9/24/94 by RRR
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int 
assemble_fill_ext_v(double afill[],	/* Jacobian matrix for fill equation  */
		    int ijaf[],	/* pointer to nonzeros in Jacobian matrix   */
		    double rf[],	/* rhs vector   */
		    double dt,	/* current time step size */
		    double tt,	/* parameter to vary time integration from
					 * explicit (tt = 1) to implicit (tt = 0) */
		    int node_to_fill[],   
		    double hsquared[DIM]) /* element size */
{
  dbl wfcn;                          /* Holds the final weight function (G, TG, or SUPG) */

  int eqn, dim, a, i, j, status, I, J, idof,jdof, ie, je, ja;
  int ki, kj, nvdofi, nvdofj, offseti;
  dbl d_grad_phi_i;
  dbl *grad_phi_i;                  /* Gradient of phi_i. */
  dbl *grad_phi_j;                  /* Gradient of phi_j. */

  dbl tmp;

  dbl F;				/* Fill. */
  dbl F_dot;				/* Fill derivative wrt time. */
  dbl *grad_F;  			/* Fill gradient. */

  dbl v_rel_n;                       /* Normal velocity relative to the mesh. */
  dbl v_rel_n_old;                   /* Old v_rel[]. */

  dbl zero[3] = {0.0, 0.0, 0.0};        /* An array of zeros, for convienience. */

  dbl d_wfcn_dF;                /* Deriv. of wfcn w.r.t. the Level Set F. */

  dbl v_dot_DF;                 /* v.grad(F) */
  dbl dtinv;                    /* = 1 / dt */

  dbl *xx;	        		/* position field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl *x_dot_old;
  dbl *x_old;		         	/* old position field. */
  
  dbl num_diff;

  /*
   * Galerkin weighting functions for i-th energy residuals 
   * and some of their derivatives...
   */
  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;
  dbl wt;
  double rhs;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /* Terms needed for GLS stabilization */

  dbl h_elem, tau_gls;
  int dofs;
  dbl ext_v, ext_v_old, ext_v_mag, ext_v_avg;
  dbl gfmag_inv;

  dbl gradF_dot_Dphi[MDE];

  dbl x_dot_n; /* scalar dot product */
  dbl x_dot_n_old; /* scalar dot product */
  
  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;
  eqn   = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  wt	= fv->wt;			/* Gauss point weight.           */
  h3	= fv->h3;			/* Differential volume element.  */
  det_J	= bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  dtinv = 1.0 / dt;                     /* Ah, 1 / dt. */

  F     = fv->F;
  F_dot = fv_dot->F;
  /* Use pointers unless we need to do algebra. */
  grad_F = fv->grad_F;
  xx	 = fv->x;
  x_old	 = fv_old->x;

  h_elem = 0.;
  for ( a=0; a<dim; a++) h_elem += hsquared[a];
  /* This is the size of the element */
  h_elem = sqrt(h_elem/ ((double )dim));

  load_lsi( 0. );
  /* find derivs, even for decoupled fill */
  load_lsi_derivs();

  gfmag_inv = 1./lsi->gfmag;
  
  dofs     = ei[pg->imtrx]->dof[EXT_VELOCITY];

  /* extension velocity magnitude for stabilization term */
  ext_v_avg = 0.;
  for (i = 0; i < dofs; i++)
    {
      ext_v_avg += *esp->ext_v[i] * *esp->ext_v[i];
    }
  ext_v_mag = sqrt(ext_v_avg/dofs);

  if(ext_v_mag == 0.)
    {
      tau_gls = 0.;
      ext_v_mag =1.;
    }

  v_dot_DF = 0.0;
  ext_v = fv->ext_v;
  ext_v_old = fv_old->ext_v;
  dtinv = 1./dt;
  if ( pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
    {
      x_dot_old = fv_dot_old->x;
      v_rel_n = 0.0;
      v_rel_n_old = 0.0;
      x_dot_n = 0.0;
      x_dot_n_old = 0.0;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	  x_dot_n += x_dot[a]* lsi->normal[a];
	  x_dot_n_old += x_dot_old[a]* lsi->normal[a];
	}
      v_rel_n      = ext_v - x_dot_n;
      v_rel_n_old  = ext_v_old - x_dot_n_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* n.v.grad(F)_mag */
 
    }
  else
    {
      x_dot_old = zero;
      for ( a=0; a < VIM; a++ ) x_dot[a] = 0.0;
      v_rel_n      = ext_v ;
      v_rel_n_old  = ext_v_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* n.v.grad(F)_mag */
   }
  tau_gls = h_elem/ext_v_mag;

  /* Compute and save v.grad(phi) and vcent.grad(phi). */
  memset(gradF_dot_Dphi, 0, sizeof(double)*MDE);
  for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
      grad_phi_i = bf[eqn]->grad_phi[i];
      for ( a=0; a < VIM; a++ ) gradF_dot_Dphi[i] += grad_phi_i[a]*grad_F[a];
    }
    
  num_diff = 0.0000;
  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  if ( af->Assemble_Residual )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	{
	  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
	  node = Nodes[I];
	  nv = node->Nodal_Vars_Info[pg->imtrx];
	  /* check for multiple dofs */
	  nvdofi = get_nv_ndofs(nv, eqn);
	  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);
	
	  for (ki = 0; ki < nvdofi; ki++)	
	    {
	      /* set rhs to zero for Dirichlet boundary conditions */
	      if (node->DBC[pg->imtrx] &&
		  (node->DBC[pg->imtrx][offseti + ki] != -1) &&
		  Debug_Flag >= 0)
		{
		  rf[node_to_fill[I] + ki] = 0.0;
		}
	      else 
		{
		  /* check to make sure that unknowns are defined at this node,
		   * otherwise don't add anything to this node */
 		  idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

		  /* also convert from node number to dof number */
		  phi_i = bf[eqn]->phi[idof];
		  
		  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ) 
		    {
		      rhs = (F_dot + v_dot_DF) * phi_i;

		    }
		  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS ) 
		    {
		      /* Make the GLS weight function */
		      /* wfcn = w + p = phi + tau_gls * (ext_v*gradF/grad_F_mag "dot" grad(phi) ) */
		      
		      grad_phi_i = bf[eqn]->grad_phi[idof];
		      
		      rhs = (F_dot + v_dot_DF) * (phi_i +
						  tau_gls*ext_v/lsi->gfmag*gradF_dot_Dphi[i]);
		    }
		  
		  /* This is common to all "tran->Fill_Weight_Fcn"s */
		  rf[node_to_fill[I] + ki] += rhs * wt * det_J * h3;
		}
	    }
	}
    } 

  if ( af->Assemble_Jacobian )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	{
	  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i]; 

	  node = Nodes[I];
	  nv = node->Nodal_Vars_Info[pg->imtrx];
	  nvdofi = get_nv_ndofs(nv, eqn);
	  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);	  
	  for (ki = 0; ki < nvdofi; ki++)
	    {
	      ie = node_to_fill[I] + ki;
	      idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
	      phi_i = bf[eqn]->phi[idof];
	      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
	      grad_phi_i = bf[eqn]->grad_phi[idof];

	      /* The weight function for Galerkin and GLS */
	      wfcn = phi_i; 
	      /* The weight function for GLS */
	      if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS )
		{
		   /* wfcn = w + p = phi + tau_gls * (ext_v*gradF/grad_F_mag "dot" grad(phi) ) */
                  wfcn += tau_gls*ext_v*gfmag_inv*gradF_dot_Dphi[i];
		}
	      
	      /* derivatives of fill equation wrt to fill variable */	  
	      for( j=0; j< ei[pg->imtrx]->num_local_nodes; j++) 
		{
		  J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j]; 
		  nvdofj = Dolphin[pg->imtrx][J][eqn];
		  for ( kj=0; kj<nvdofj; kj++)

		    {
		      
		      jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
		      phi_j = bf[eqn]->phi[jdof];
		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[eqn]->grad_phi[jdof];

		      
		      je = node_to_fill[J] + kj;
		      ja     = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie+1], ijaf);
		      EH(ja, "Could not find vbl in sparse matrix.");
		      
		      /* set diagonal to one for Dirichlet boundary conditions */
		      if (node->DBC[pg->imtrx] && 
			  ((int) node->DBC[pg->imtrx][offseti + ki] != -1) &&
			  Debug_Flag >= 0)
			{
			  afill[ja] = delta(i,j); 
			}
		      else
			{
			  /*  afill[ja] += wt * det_J * phi_i * phi_j/dt*h3; */

			  /*
			   * Use the selected Finite Element Method
			   */
			  /* Galerkin terms */

			  tmp = wfcn*(phi_j * (1. + 2.*tt) * dtinv + 
				      v_rel_n * lsi->d_gfmag_dF[j] );
			  
			  /* Extra terms for GLS derivatives */
			  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS )
			    {
                              d_wfcn_dF = -tau_gls*ext_v*gradF_dot_Dphi[i]*lsi->d_gfmag_dF[j]*
                                (gfmag_inv*gfmag_inv);
			      for ( a=0; a<VIM; a++)
				{
                                  d_wfcn_dF += tau_gls*ext_v*gfmag_inv*grad_phi_i[a]* grad_phi_j[a];
				}
			      tmp += d_wfcn_dF*(F_dot + v_dot_DF);
			      
			    }

			  afill[ja] += wt*h3*det_J*tmp;

			}
		    }
		}
	    }
	}
    }

  return(status);
} /* end of assemble_fill_ext_v */

/*
 * assemble_fill_ext_v -- integrate explicit fill equation with Taylor-Galerkin
 *                  formulation.  This routine is for the de-coupled method.
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	9/24/94 by RRR
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int 
assemble_fill_gradf(double afill[],	/* Jacobian matrix for fill equation  */
		    int ijaf[],	/* pointer to nonzeros in Jacobian matrix   */
		    double rf[],	/* rhs vector   */
		    double dt,	/* current time step size */
		    double tt,	/* parameter to vary time integration from
					 * explicit (tt = 1) to implicit (tt = 0) */
		    int node_to_fill[],   
		    double hsquared[DIM]) /* element size */
{
  dbl wfcn;                          /* Holds the final weight function (G, TG, or SUPG) */

  int eqn, dim, a, i, j, status, I, J, idof,jdof, ie, je, ja;
  int ki, kj, nvdofi, nvdofj, offseti;
  dbl d_grad_phi_i;
  dbl *grad_phi_i;                  /* Gradient of phi_i. */
  dbl *grad_phi_j;                  /* Gradient of phi_j. */

  dbl tmp;

  dbl F;				/* Fill. */
  dbl F_dot;				/* Fill derivative wrt time. */
  dbl *grad_F;  			/* Fill gradient. */

  dbl v_rel_n;                       /* Normal velocity relative to the mesh. */
  dbl v_rel_n_old;                   /* Old v_rel[]. */

  dbl zero[3] = {0.0, 0.0, 0.0};        /* An array of zeros, for convienience. */

  dbl d_wfcn_dF;                /* Deriv. of wfcn w.r.t. the Level Set F. */

  dbl v_dot_DF;                 /* v.grad(F) */
  dbl dtinv;                    /* = 1 / dt */

  dbl *xx;	        		/* position field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl *x_dot_old;
  dbl *x_old;		         	/* old position field. */
  
  dbl num_diff;

  /*
   * Galerkin weighting functions for i-th energy residuals 
   * and some of their derivatives...
   */
  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;
  dbl wt;
  double rhs;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /* Terms needed for GLS stabilization */

  dbl h_elem, tau_gls;
  int dofs;
  dbl ext_v, ext_v_old, ext_v_mag, ext_v_avg;
  dbl gfmag_inv;

  dbl gradF_dot_Dphi[MDE];

  dbl x_dot_n; /* scalar dot product */
  dbl x_dot_n_old; /* scalar dot product */
  
  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;
  eqn   = R_FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(status);
    }

  wt	= fv->wt;			/* Gauss point weight.           */
  h3	= fv->h3;			/* Differential volume element.  */
  det_J	= bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  dtinv = 1.0 / dt;                     /* Ah, 1 / dt. */

  F     = fv->F;
  F_dot = fv_dot->F;
  /* Use pointers unless we need to do algebra. */
  grad_F = fv->grad_F;
  xx	 = fv->x;
  x_old	 = fv_old->x;

  h_elem = 0.;
  for ( a=0; a<dim; a++) h_elem += hsquared[a];
  /* This is the size of the element */
  h_elem = sqrt(h_elem/ ((double )dim));

  load_lsi( 0. );
  /* find derivs, even for decoupled fill */
  load_lsi_derivs();

  gfmag_inv = 1./lsi->gfmag;
  
  dofs     = ei[pg->imtrx]->dof[EXT_VELOCITY];

  /* extension velocity magnitude for stabilization term */
  ext_v_avg = 0.;
  for (i = 0; i < dofs; i++)
    {
      ext_v_avg += *esp->ext_v[i] * *esp->ext_v[i];
    }
  ext_v_mag = sqrt(ext_v_avg/dofs);

  if(ext_v_mag == 0.)
    {
      tau_gls = 0.;
      ext_v_mag =1.;
    }

  v_dot_DF = 0.0;
  ext_v = fv->ext_v;
  ext_v_old = fv_old->ext_v;
  dtinv = 1./dt;
  if ( pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
    {
      x_dot_old = fv_dot_old->x;
      v_rel_n = 0.0;
      v_rel_n_old = 0.0;
      x_dot_n = 0.0;
      x_dot_n_old = 0.0;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	  x_dot_n += x_dot[a]* lsi->normal[a];
	  x_dot_n_old += x_dot_old[a]* lsi->normal[a];
	}
      v_rel_n      = ext_v - x_dot_n;
      v_rel_n_old  = ext_v_old - x_dot_n_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* n.v.grad(F)_mag */
 
    }
  else
    {
      x_dot_old = zero;
      for ( a=0; a < VIM; a++ ) x_dot[a] = 0.0;
      v_rel_n      = ext_v ;
      v_rel_n_old  = ext_v_old;
      v_dot_DF     = v_rel_n * lsi->gfmag;      /* n.v.grad(F)_mag */
   }
  tau_gls = h_elem/ext_v_mag;

  /* Compute and save v.grad(phi) and vcent.grad(phi). */
  memset(gradF_dot_Dphi, 0, sizeof(double)*MDE);
  for ( i=0; i < ei[pg->imtrx]->dof[eqn]; i++ )
    {
      /* So: grad_phi_i[a] == bf[var]->grad_phi[i][a] */
      grad_phi_i = bf[eqn]->grad_phi[i];
      for ( a=0; a < VIM; a++ ) gradF_dot_Dphi[i] += grad_phi_i[a]*grad_F[a];
    }
    
  num_diff = 0.0000;
  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  if ( af->Assemble_Residual )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	{
	  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
	  node = Nodes[I];
	  nv = node->Nodal_Vars_Info[pg->imtrx];
	  /* check for multiple dofs */
	  nvdofi = get_nv_ndofs(nv, eqn);
	  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);
	
	  for (ki = 0; ki < nvdofi; ki++)	
	    {
	      /* set rhs to zero for Dirichlet boundary conditions */
	      if (node->DBC[pg->imtrx] &&
		  (node->DBC[pg->imtrx][offseti + ki] != -1) &&
		  Debug_Flag >= 0)
		{
		  rf[node_to_fill[I] + ki] = 0.0;
		}
	      else 
		{
		  /* check to make sure that unknowns are defined at this node,
		   * otherwise don't add anything to this node */
 		  idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;

		  /* also convert from node number to dof number */
		  phi_i = bf[eqn]->phi[idof];
		  
		  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ) 
		    {
		      rhs = (F_dot + v_dot_DF) * phi_i;

		    }
		  else if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS ) 
		    {
		      /* Make the GLS weight function */
		      /* wfcn = w + p = phi + tau_gls * (ext_v*gradF/grad_F_mag "dot" grad(phi) ) */
		      
		      grad_phi_i = bf[eqn]->grad_phi[idof];
		      
		      rhs = (F_dot + v_dot_DF) * (phi_i +
						  tau_gls*ext_v/lsi->gfmag*gradF_dot_Dphi[i]);
		    }
		  
		  /* This is common to all "tran->Fill_Weight_Fcn"s */
		  rf[node_to_fill[I] + ki] += rhs * wt * det_J * h3;
		}
	    }
	}
    } 

  if ( af->Assemble_Jacobian )
    {
      for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
	{
	  I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i]; 

	  node = Nodes[I];
	  nv = node->Nodal_Vars_Info[pg->imtrx];
	  nvdofi = get_nv_ndofs(nv, eqn);
	  offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);	  
	  for (ki = 0; ki < nvdofi; ki++)
	    {
	      ie = node_to_fill[I] + ki;
	      idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
	      phi_i = bf[eqn]->phi[idof];
	      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
	      grad_phi_i = bf[eqn]->grad_phi[idof];

	      /* The weight function for Galerkin and GLS */
	      wfcn = phi_i; 
	      /* The weight function for GLS */
	      if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS )
		{
		   /* wfcn = w + p = phi + tau_gls * (ext_v*gradF/grad_F_mag "dot" grad(phi) ) */
                  wfcn += tau_gls*ext_v*gfmag_inv*gradF_dot_Dphi[i];
		}
	      
	      /* derivatives of fill equation wrt to fill variable */	  
	      for( j=0; j< ei[pg->imtrx]->num_local_nodes; j++) 
		{
		  J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j]; 
		  nvdofj = Dolphin[pg->imtrx][J][eqn];
		  for ( kj=0; kj<nvdofj; kj++)

		    {
		      
		      jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
		      phi_j = bf[eqn]->phi[jdof];
		      /* So: grad_phi_j[a] == bf[var]->grad_phi[j][a] */
		      grad_phi_j = bf[eqn]->grad_phi[jdof];

		      
		      je = node_to_fill[J] + kj;
		      ja     = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie+1], ijaf);
		      EH(ja, "Could not find vbl in sparse matrix.");
		      
		      /* set diagonal to one for Dirichlet boundary conditions */
		      if (node->DBC[pg->imtrx] && 
			  ((int) node->DBC[pg->imtrx][offseti + ki] != -1) &&
			  Debug_Flag >= 0)
			{
			  afill[ja] = delta(i,j); 
			}
		      else
			{
			  /*  afill[ja] += wt * det_J * phi_i * phi_j/dt*h3; */

			  /*
			   * Use the selected Finite Element Method
			   */
			  /* Galerkin terms */

			  tmp = wfcn*(phi_j * (1. + 2.*tt) * dtinv + 
				      v_rel_n * lsi->d_gfmag_dF[j] );
			  
			  /* Extra terms for GLS derivatives */
			  if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_GLS )
			    {
                              d_wfcn_dF = -tau_gls*ext_v*gradF_dot_Dphi[i]*lsi->d_gfmag_dF[j]*
                                (gfmag_inv*gfmag_inv);
			      for ( a=0; a<VIM; a++)
				{
                                  d_wfcn_dF += tau_gls*ext_v*gfmag_inv*grad_phi_i[a]* grad_phi_j[a];
				}
			      tmp += d_wfcn_dF*(F_dot + v_dot_DF);
			      
			    }

			  afill[ja] += wt*h3*det_J*tmp;

			}
		    }
		}
	    }
	}
    }

  return(status);
} /* end of assemble_fill_gradf */

#endif /* COUPLED_FILL */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


#ifndef COUPLED_FILL
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int
integrate_explicit_eqn(
     struct Aztec_Linear_Solver_System *ams,       /* cf "sl_util_structs.h" */
     double rf[],                                  /* residual for fill
						      equation only     */
     double xf[],                                  /* vector with fill at
						      nodes only     */
     double xf_old[],
     double xfdot[],
     double xfdot_old[],
     double x[],		                   /* Solution vector for
						      the current processor */
     double x_old[],                               /* Solution vector at last
						      time step            */
     double xdot[],                                /* time derivative of
						      solution vector */
     dbl delta_t,		                   /* time step size */
     dbl theta,	             /* parameter to vary time integration from 
			        explicit (theta = 1) to implicit (theta = 0) */
     dbl *time2,
     int eqntype,
     int node_to_fill[],                           /* this is a map from the 
						    * node number to the 
						    * unnown number for the
						    * fill equation */
     Exo_DB *exo,                                  /* ptr to exodus file */
     Dpi *dpi,					   /* distributed processing 
						    * info */
     Comm_Ex *cx )				   /* cx - array of structs, 
						    * one for each neighbor */

/******************************************************************************
  Function which integrates the explicit fill equation using a segregated
  vector/matrix that has fill variables in it only

  Author:          R. R. Rao
  Date:            5 October 1996
  Revised:         

*******************************************************************************/
{
  int 	        i;
  int 		inewton;		/* Newton iteration counter */

  int		return_value;	        /* non zero if things screw up ...  */

  double	Norm[2][3];		/* Global norms... */
					/*   [0][0] == residual, L_oo norm */
					/*   [0][1] == residual, L_1  norm */
					/*   [0][2] == residual, L_2  norm */
					/*   [1][0] == correction, L_oo norm */
					/*   [1][1] == correction, L_1  norm */
					/*   [1][2] == correction, L_2  norm */

  int		 prnt_lev_linsolve;	/* Level of printing for the linear */
					/* solve funcs: */
					/* 0 -- no printing */
					/* 1 -- print brief summary */
					/* 2 -- print line @ every iteration */

  int       printing = ( (eqntype != CORRECT)  );

  static char yo[] = "integrate_explicit_eqn"; /* My name to take blame... */

  /* 
   * The following vectors all have lengths equal to the number of local
   * unknowns - i.e., equations for which the current processor is responsible.
   */

  double 	*delta_x;		/* update */


  int		num_unk_r, num_unk_x; 
  char dofname_r[80];
  char dofname_x[80];

  int		num_total_nodes; 
  int		Max_Steps; 

  static int UMF_system_id;	/* Used to uniquely identify the
				 * explicit fill system to solve from
				 * the other UMF systems. */
  int		Factor_Flag;	/* Legend::
				 * --------
				 * 1 -- this is the first LU,
				 * 2 -- >1st LU, similar matrix structure,
				 * 3 -- new rhs, use already factored LU

				 UMFPACK USAGE:

				 0   : LOAD MATRIX AND ANALYSIS/DECOMPOSITION 
				       AND BACK SUBSTITUTION
				 1   : LOAD MATRIX AND DECOMPOSITION USING PAST ANALYSIS 
				       AND BACK SUBSTITUTION
				 > 1 : BACK SUBSTITUTION ONLY

				 */

  int matr_form;                /* 1: MSR FORMAT MATRIX FOR UMFPACK DRIVER */

  int       err, why;

  Spfrtn spfrtn;
  dbl		lits;		/* number of linear solver iterations taken */

  char		stringer[80];	/* holding format of num linear solve itns */

  dbl		a_start;	/* mark start of assembly */
  dbl		a_end;		/* mark end of assembly */

  char		ctod[80];	/* hold current time of day */
  dbl		s_start;	/* mark start of solve */
  dbl		s_end;		/* mark end of solve */


  int	linear_solver_blk;	/* count calls to AZ_solve() */
  int	linear_solver_itns;	/* count cumulative linearsolver iterations */
  int	num_linear_solve_blks;	/* one pass for now */
  int	matrix_solved;		/* boolean */
  int   converged;

  int local_order;		/* of the unknowns that this processor owns */
  int local_order_plus;		/* and if the external rows are included */
  int local_nnz;		/* number of nonzero matrix entries that
				 * are owned by this processor */
  int local_nnz_plus;		/* and if the external rows are included */

  int global_order;		/* order of the global system */
  int global_order_plus;	/* and if the external rows are overincluded */
  int global_nnz;		/* a sum of the number of nonzero matrix 
				 * entries owned by each processor */
  int global_nnz_plus;		/* a sum that overincludes the external rows */

  double *scale = NULL; 	/* scale vector for row sum scaling */


  /*
   * Pointers for convenience purposes. Now, the matrix system is
   * imported via the "ams" structure. These pointers help smooth
   * the transition for the old usage.
   */

  double *afill; 	/* matrix for fill variables only ("val") */

  int *ijaf; 	/* pntr to nonzeros in fill matrix ("bindx") */
  int *ijaf_save; 	/* pntr to external nonzeros in fill matrix ("bindx"),
			 * a save area while Aztec crunches the distributed 
			 * matrix */

  extern struct elem_side_bc_struct **First_Elem_Side_BC_Array;

/*
 * Begin executable statements...
 */

  afill     = ams->val;
  ijaf      = ams->bindx;
  ijaf_save = ams->belfry;
  return_value = 0;		/* OK to begin */
  prnt_lev_linsolve = Iout;

  /*
   * The number of nonzero matrix entries that this processor owns.
   * External rows are not counted, though ija[] does stretch out to 
   * include them too for the parallel case. In the parallel case,
   * watch out to use the bigger ija[] for assembly purposes, but the
   * smaller one for Aztec (after processing via hide_external()).
   */

  /* PAS-warning: this should be based on the number of unknowns now
   * and not the nodes, so that it will work for discontinuous variables
   * as well
   */

  local_order      = owned_fill_unknowns;
  local_order_plus = num_fill_unknowns;

  local_nnz        = ijaf[local_order]-1;
  local_nnz_plus   = ijaf[local_order_plus]-1;

  canine_chaos(local_order, local_order_plus, local_nnz, local_nnz_plus,
	       &global_order, &global_order_plus, 
	       &global_nnz, &global_nnz_plus);

  /*
   * Setup some legacy variables.
   */

  num_total_nodes     = dpi->num_universe_nodes;
  fill_zeros          = local_nnz_plus; /* big estimate for assembly*/
  


  asdv(&delta_x, num_fill_unknowns);
  asdv(&scale, num_fill_unknowns);
  
  /*
   * Initial conditions at beginning of Newton iteration loop...
   */
  
  converged = FALSE;
  inewton   = 0;
  Max_Steps = 6;

  if( eqntype == CORRECT ) Max_Steps = 6;

  /*
   *  HKM->
   *  We need to make the incoming estimate of xfdot[] consistent
   *  with the estimate of xf[] and the incoming values of
   *  xf_old[] and xfdot_old[]. Note, the current implementation
   *  assignes xf[] = xf_old[]. However, we don't necessarily have
   *  to be this crude in the future.
   *  xfdot[]'s value must be propogated to xdot[] for the values
   *  to be represented in the FILL matrix and residuals.
   */


  for (i = 0; i < num_fill_unknowns; i++) {
    xfdot[i] = (1.0 + 2.0*theta) / delta_t  * (xf[i] - xf_old[i]) -
      2.0 * theta * xfdot_old[i];
  }
  put_fill_vector(num_total_nodes, xdot, xfdot, node_to_fill);
  exchange_dof(cx, dpi, xdot);  

  if (printing) {
    DPRINTF(stderr, "\n\t   L_2 in R   L_2 in dx   lis \n");
    DPRINTF(stderr, "\t ---------  ---------   -----\n");
  }
  
  /*
   * Top of the Newton's loop 
   */ 
  while ( ( ! converged ) && ( inewton < Max_Steps ) )
    {      
      /*
       * Re-zero vectors that would otherwise accumulate junk over more than
       * one Newton iteration...
       */
      if(Debug_Flag == -4)
	{

	  af->Assemble_Residual = TRUE;
	  af->Assemble_Jacobian = TRUE;
	  af->Assemble_LSA_Jacobian_Matrix = FALSE;
	  af->Assemble_LSA_Mass_Matrix = FALSE;

	  numerical_jacobian_fill(ijaf, afill, xf, rf, delta_t, theta,
				  x, x_old, xdot, Debug_Flag, node_to_fill,
				  exo, dpi);
	  if(printing) DPRINTF(stderr,"%s: numerical Jacobian done....\n", yo);
	  exit(-1);
	}
      else
	{
	  af->Assemble_Residual = TRUE;
	  af->Assemble_Jacobian = TRUE;
	  af->Assemble_LSA_Jacobian_Matrix = FALSE;
	  af->Assemble_LSA_Mass_Matrix = FALSE;

	  /*
	   * Initialize the rhs and matrix for the FILL equation
	   * problem to zero.
	   */
	  init_vec_value(rf,    0.0, num_fill_unknowns);
	  init_vec_value(afill, 0.0, (fill_zeros + 1) );
	  
	  get_time(ctod);
	  if (printing) DPRINTF(stderr,"\t   [%d] ", inewton );
	  
	  /*
	   * Formulate the rhs and matrix for the FILL equation
	   * by calling special matrix fill routines
	   */
	  a_start = ut();
	  err = fill_matrix(afill, ijaf, rf, xf, x, x_old, xdot, delta_t,
			    theta, eqntype, node_to_fill, 
			    First_Elem_Side_BC_Array, exo, dpi);
	  a_end = ut();
	  EH( err, "fill_matrix");
          
          /* Scale matrix first to get rid of problems with 
             penalty parameter. */
          row_sum_scale_MSR(num_fill_unknowns, afill, ijaf, rf, scale);
          
	}
      
      /*
       * We already know the resid vector now, so let's print it out
       * now instead of after a big long matrix elimination...
       */
      
      Norm[0][0] = Loo_norm(rf, owned_fill_unknowns, &num_unk_r, dofname_r);
      Norm[0][1] = L1_norm (rf, owned_fill_unknowns);
      Norm[0][2] = L2_norm (rf, owned_fill_unknowns);
      
      if (printing) DPRINTF(stderr, "%10.2e ", Norm[0][2]);
      
      /*
       * Before starting the solve of the linear system, fix up the matrix
       * to exclude external rows. However, we'll save all the column name 
       * values so we can reincorporate them back in before the next assembly 
       * step.
       */

      hide_external(num_fill_unknowns, owned_fill_unknowns, ijaf, ijaf_save,
                    afill);

      s_start = ut();
      
      /*
       * SOLVE THE LINEAR SYSTEM using one of our various solvers
       */
      switch (Linear_Solver)
	{
	case UMFPACK2:
	  if (first_linear_solver_call)
	    {
	      Factor_Flag = 0;
	      UMF_system_id = -1;
	    }
	  if(Linear_Solver == UMFPACK2F)
	    Factor_Flag = 0;
	  /*  */
	  matr_form = 1;
	  UMF_system_id =  SL_UMF(UMF_system_id,
				  &first_linear_solver_call,
				  &Factor_Flag, &matr_form,
				  &num_fill_unknowns, &local_nnz,
				  ijaf, ijaf, afill, rf, delta_x);
	  /*  */
	  first_linear_solver_call = FALSE;
	  Factor_Flag = 1;
	  /*  */
	  strcpy(stringer, " 1 ");
	  break;

	case SPARSE13a:
	case FRONT:  /* That's right if your using FRONT, LU solves fill equation */
	case UMFPACK2F:  /* Also if you are using UMF2 */
	case AMESOS:     /* for the moment we shall default the FILL solve to LU when using AMESOS */
	  dcopy1( num_fill_unknowns, rf, delta_x);
	  if(inewton == 0)
	    {
	      luf(num_fill_unknowns, 0, fill_zeros, 
		  afill, ijaf, delta_x, (first_linear_solver_call?1:2));
	      first_linear_solver_call = FALSE;
	    }
	  else
	    {
	      luf(num_fill_unknowns, 0, fill_zeros, 
		  afill, ijaf, delta_x, 3);
	      first_linear_solver_call = FALSE;
	    }

	  strcpy(stringer, " 1 ");
	  break;
	  
	case AZTEC:

	  /* hardwire tolerance for now */

	  ams->params[AZ_tol] = 1.e-12;
	  ams->options[AZ_precond] = AZ_dom_decomp;
	  ams->options[AZ_subdomain_solve] = AZ_ilu;

	  ams->data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
	  ams->mat_type                 = AZ_MSR_MATRIX;
	  ams->data_org[AZ_N_int_blk]   = ams->N_update;
	  ams->data_org[AZ_N_bord_blk]  = 0;
	  ams->data_org[AZ_N_ext_blk]   = ams->N_update;  

	  if ( first_linear_solver_call )
	    {
	      ams->options[AZ_pre_calc] = AZ_calc;
	    }
	  else
	    {
	      if ( strcmp(Matrix_Factorization_Reuse, "calc") == 0 )
		{
		  /*
		   * Gonna start from scratch even though I've cooked a
		   * preconditioner in the kitchen all day? Well, then
		   * you won't need the leftover pieces from all my
		   * hard preparation last time around.
		   */

		  AZ_free_memory(ams->data_org[AZ_name]); 

		  ams->options[AZ_pre_calc] = AZ_calc;

		}
	      else if ( strcmp(Matrix_Factorization_Reuse, "recalc") == 0 )
		{
		  ams->options[AZ_pre_calc] = AZ_recalc;
		}
	      else if ( strcmp(Matrix_Factorization_Reuse, "reuse") == 0 )
		{
		  ams->options[AZ_pre_calc] = AZ_reuse;
		}
	      else
		{
		  EH(-1, "Unknown factorization reuse specification.");
		}
	    }
	  
	  linear_solver_blk     = 0; /* count calls to AZ_solve() */
	  num_linear_solve_blks = 1; /* upper limit to AZ_solve() calls */
	  linear_solver_itns    = 0; /* cumulative number of iterations */
	  matrix_solved         = FALSE; 
	  while ( ( ! matrix_solved                            ) && 
		  ( linear_solver_blk < num_linear_solve_blks  ) )
	    {
	      /* 
	       * Someday the user may want to do fancy heuristics based
	       * on all kinds of cost functions, artificial intelligence
	       * neural networks, etc.
	       *
	       * For the linear system "Ax=b", we have
	       *    A -- indx, bindx(ija), rpntr, cpntr, bpntr, val(a)
	       *    x -- delta_x, newton correction vector
	       *    b -- resid_vector, newton residual equation vector
	       */
	      AZ_solve(delta_x, rf, ams->options, ams->params, 
		       ams->indx, ams->bindx, ams->rpntr, ams->cpntr, 
		       ams->bpntr, ams->val, ams->data_org, ams->status, 
		       ams->proc_config);

	      first_linear_solver_call = FALSE;

	      if ( Debug_Flag > 0 )
		{
		  dump_aztec_status(ams->status);
		}
	      
              why = (int)ams->status[AZ_why];
              aztec_stringer(why, ams->status[AZ_its], &stringer[0]);
	      
	      matrix_solved = ( ams->status[AZ_why] == AZ_normal) ;
	      linear_solver_blk++;
	      linear_solver_itns += ams->status[AZ_its];
	      
	    }
	  AZ_free_memory(ams->data_org[AZ_name]); 

	  /*	  if (ams->options[AZ_pre_calc] == AZ_calc) { */
	    /* FREE the memory used in storing preconditioner info
	     *   - unless using the RE_USE option */
	    /* 	      AZ_clear(ams->data_org[AZ_name]); */
	  /*	    AZ_free_memory(ams->data_org[AZ_name]); */
	  /*	  } */

	  break;
	  
	case MA28:
	  /*
	   * sl_ma28 keeps interntal static variables to determine whether
	   * it is the first call or not.
	   */
#ifdef HARWELL	  
	  err = cmsr_ma28 (num_fill_unknowns, local_nnz, afill, ijaf, 
			   delta_x, rf);
#endif /* HARWELL */
#ifndef HARWELL
	  EH(-1, "That linear solver package is not implemented.");
#endif /* HARWELL */
	  strcpy(stringer, " 1 ");
	  break;
	  
	default:
	  EH(-1, "That linear solver package is not implemented.");
	  break;
	}
      
      s_end = ut();
      
      /*
       * Back to normal. Readjust the ija[] pointers to mimic the bloated
       * system for this processor that pretends like it really cares about
       * filling in matrix entries for external nodal fill degrees of freedom.
       */

      show_external(num_fill_unknowns, owned_fill_unknowns, ijaf, ijaf_save,
                    afill);

      Norm[1][0] = Loo_norm(delta_x, owned_fill_unknowns,&num_unk_x, dofname_x);
      Norm[1][1] = L1_norm (delta_x, owned_fill_unknowns);
      Norm[1][2] = L2_norm (delta_x, owned_fill_unknowns);
      
      if (printing) DPRINTF(stderr, "%10.2e   %s \n", 
			    Norm[1][2], stringer);
      
      /*
       *  Update the value of the FILL unknowns on this processor from
       *  the delta_x[] just computed in the newton iteration
       */
      for (i = 0; i < num_fill_unknowns; i++) {
	xf[i] -= delta_x[i];
	xfdot[i] = (1.0 + 2.0*theta) / delta_t  * (xf[i] - xf_old[i]) -
	    2.0 * theta * xfdot_old[i];
      }
     
      /*
       *  Update the main solution vector, x[] with the value of the
       *  FILL variables, exchange the values between owned nodes and
       *  ghost nodes and recompute xf[]
       */
      put_fill_vector(num_total_nodes, x, xf, node_to_fill);
      exchange_dof(cx, dpi, x);
      put_fill_vector(num_total_nodes, xdot, xfdot, node_to_fill);
      exchange_dof(cx, dpi, xdot);


      get_fill_vector(num_total_nodes, x, xf, node_to_fill);
      get_fill_vector(num_total_nodes, xdot, xfdot, node_to_fill);
      converged = (Norm[0][2] <  Epsilon[0]);
      inewton++;
    }
  
  return_value = inewton;
  
  if ( ! converged ) 
    {
      if ( Debug_Flag  && ProcID==0 ) 
	{ 
	  if(printing)  DPRINTF(stderr, "\n%s:  Newton iteration FAILED.\n", yo);
	}

      return_value *= -1;
    }
  else
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  if(printing)  DPRINTF(stderr,"\n%s:  Newton iteration CONVERGED.\n", yo);
	  
	  /*  writesolution...to soln.dat and to exodus ii db.... */
	  
	}
    }
  
  safe_free( (void *) delta_x);
  safe_free( (void *) scale);

  
  return(return_value);
/*****************************************************************************/
}   /*   END OF integrate_explicit_eqn                                       */

   

#endif /* not COUPLED_FILL */

/*****************************************************************************/
/******************************************************************************/
/******************************************************************************/
int
assemble_fill_fake( double tt, double dt)

    /*************************************************************************
     *
     * assemble_fill_fake()
     *
     *   Dummy routine that puts a 1 on the diagonal of the fill equation
     *   and sets the residual equal to zero. The actual calculation takes
     *   place before the main time transient calculation. We seek here only
     *   to create a situation where the solution storred in the fill unknowns
     *   in x[] remains unchanged.
     *************************************************************************/
{
  int eqn, var, peqn, pvar, i;

  /*
   * Unpack variables from structures for local convenience...
   */
  eqn  = ls->var;
  peqn = upd->ep[pg->imtrx][eqn];

  /*
   * Bail out fast if there's nothing to do...
   */
  if (! pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  /*
   * Put local contributions into global right-hand side
   * if it is not a right-hand side variable-it won't get added in (contribution
   * is zero)
   */
  var  = ls->var;
  pvar = upd->vp[pg->imtrx][var];
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    lec->R[peqn][i] = 0.;
    /*
     *    J_f_F 
     */
    lec->J[peqn][pvar][i][i] = 1.;	
  }
  return (0);
} /* end of assemble_fill_fake */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
#ifndef COUPLED_FILL

int
fill_matrix(double afill[],	/* matrix for fill variables only      */
	    int ijaf[],		/* pointer to nonzeros in fill matrix  */
	    double rf[],	/* residual for fill equation only     */
	    double xf[],	/* vector with fill at nodes only     */
	    double x[],		/* Solution vector for the current proc   */
	    double x_old[],	/* Solution vector at last time step        */
	    double xdot[],	/* time derivative of solution vector        */
	    dbl delta_t,
	    dbl theta,
	    int eqntype,        /* type of equation being integrated */
	    int node_to_fill[],
	    struct elem_side_bc_struct *first_elem_side_BC_array[],
	    Exo_DB *exo,	/* ptr to finite element database */
	    Dpi *dpi)			/* distributed processing info */

/******************************************************************************
  Function which has the element loop, volume integration and surface
  integration loops for the fill equation

  Author:          R. R. Rao
  Date:            5 October 1996
  Revised:         

   TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.

******************************************************************************/
{
  int ebi;			/* element block index */
  int err;                    /* temp variable to hold diagnostic flags. */
  int ip;                     /* ip is the local quadrature point index       */
  int ip_total;               /* ip_total is the total number of volume
                                 quadrature points in the element             */
  int ielem, iel;             /* Index for elements                           */
  int ielem_type;             /* Element type of the current element          */
  int ielem_type_fill = -1;	/* Element type of the fill element          */
  int ielem_dim;              /* Element physical dimension                   */
  int num_local_nodes;        /* Number of local basis functions in the
                                 current element                              */
  int num_total_nodes;        /* Number of nodes that each processor is
                                 responsible for                              */
  int iconnect_ptr;           /* Pointer to the beginning of the connectivity
                                 list for the current element                 */

  double s, t, u;             /* Gaussian-quadrature point locations          */
  
  double xi[DIM];               /* Local element coordinates of Gauss point. */
  
  double wt;                  /* Quadrature weights
                                 units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  dbl hsquared[DIM];          /* square of the element size variable for GLS */
  dbl hhv[DIM][DIM];
  dbl dhv_dxnode[DIM][MDE];

  int e_start, e_end, mn;
  int id_side;
  int discontinuous;
  const int dim = pd->Num_Dim;

  static int inlet_bc_list[MAX_INLET_BC] ;  /* list of FILL_INLET conditions BCs indexs */
  static int bc_list_built = 0;   /* flag to indicate the inlet_bc_list has been CONSTRUCTED */

/******************************************************************************/
/*                                BLOCK 0                                     */
/*                         INITIALIZATION FOR ROUTINE                         */
/*      INITIALIZATION THAT IS DEPENDENT ON THE CURRENT PROCESSOR NUMBER      */
/******************************************************************************/

  /* Find the total number of local nodes for the current processor.  Note that
   * this is also the maximum value of the local node for which this processor
   * is responsible for.  Nodes, whose local node number on that processor is
   * greater than num_total_nodes, are owned by other processors.            
   */

  /*
   * Let it slide and pretend as if we're responsible for those external
   * nodes. In the end, the data exchange will override our external unknowns
   * anyway.
   */
  num_total_nodes = dpi->num_universe_nodes;

  ei[pg->imtrx]->elem_blk_id         = -1234; /* Hardwired to ridiculous, for now. -pas */

/******************************************************************************/
/*                                BLOCK 1                                     */
/*          LOOP OVER THE ELEMENTS DEFINED ON THE CURRENT PROCESSOR           */
/*          INITIALIZATION THAT IS DEPENDENT ON THE ELEMENT NUMBER            */
/*          First calculate field variables using volume integrals            */
/*          Then calculate stream and flux functions (need fields of          */
/*          diffusion and conduction vectors for the flux functions           */
/******************************************************************************/

  /* 
   * Loop over all the elements, one element block (matl) at a time, calculating
   * the required interaction coefficients.
   */

  for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
    {
      ei[pg->imtrx]->elem_blk_id	  = exo->eb_id[ebi];
      mn                  = Matilda[ebi];
      
      /* 
       * Assign local pointer pd to appropriate material
       */
      
      pd  = pd_glob[mn];
      cr  = cr_glob[mn];
      elc = elc_glob[mn];
      elc_rs = elc_rs_glob[mn];
      gn  = gn_glob[mn];
      mp  = mp_glob[mn];
      
      e_start = exo->eb_ptr[ebi];
      e_end   = exo->eb_ptr[ebi+1];
      
      /* element type, assumed to be consistent in each block */
      ielem_type      = Elem_Type(exo, e_start);
      
      discontinuous = 0;
      if(pd->i[pg->imtrx][FILL]==I_P1)
	{
	  if (dim == 2) ielem_type_fill = P1_QUAD;
	  if (dim == 3) ielem_type_fill = P1_HEX;
	  discontinuous = 1;
	}
      else if(pd->i[pg->imtrx][FILL]==I_P0)
	{
	  if (dim == 2) ielem_type_fill = P0_QUAD;
	  if (dim == 3) ielem_type_fill = P0_HEX;
	  discontinuous = 1;
	}
      else if(pd->i[pg->imtrx][FILL]==I_PQ1)
	{
	  if (dim == 2) ielem_type_fill = C_BILINEAR_QUAD;
	  if (dim == 3) EH(-1,"Sorry PQ1 interpolation has not been implemented in 3D yet.");
	  discontinuous = 1;
	}
      else if(pd->i[pg->imtrx][FILL]==I_PQ2)
	{
	  if (dim == 2) ielem_type_fill = BIQUAD_QUAD;
	  if (dim == 3) EH(-1,"Sorry PQ2 interpolation has not been implemented in 3D yet.");
	  discontinuous = 1;
	}
      else
	{
	  ielem_type_fill = ielem_type;
	}

/******************************************************************************/
/*                                BLOCK 1.5                                    */
/*          Do some preliminary side set checking to identify                 */
/*          the  FILL_INLET BCs indecies                                      */
/******************************************************************************/

      if( discontinuous && !bc_list_built )
	{
	  int k=0, bc_input_id=0;
	  
	  for(bc_input_id = 0; bc_input_id < Num_BC; bc_input_id++ )
	    {
	      if( BC_Types[bc_input_id].BC_Name == FILL_INLET_BC ) 
		{
		  inlet_bc_list[k++] = bc_input_id;
		}
	      
	      if( k == MAX_INLET_BC )
		{
		  EH(-1,"Maximum FILL_INLET BCs exceeded.  Adjust MAX_INLET_BC in mm_fill_fill.c");
		}
	    }
	  
	  inlet_bc_list[k] = -8458912;  /* My office phone number marks the end of the list.  TAB */
	  bc_list_built = 1;
	}

	  
	  for (iel = e_start; iel < e_end; iel++)
	    {
	      int fill_inlet_side = -1, fill_inlet_bc=-1;
	      int k;


	      ielem = iel;

	      /*
	       * For each variable there are generally different degrees of
	       * freedom that they and their equations contribute to.
	       *
	       * For this element, load up arrays that tell, for each variable,
	       *
	       *        (i) how many degrees of freedom they contribute towords
	       *        (ii) the local node number associated with each degree of
	       *             freedom
	       *        (iii) pointers in the "esp" structure that tell where
	       *              things are located in the global scheme...
	       *                (a) nodal unknowns in this element...
	       *                (b) Residual equations receiving contributions in
	       *                    this element.
	       *                (c) where the Jacobian entries go in the global
	       *                    "a" matrix in its MSR format...
	       */
	      	  	      
	      err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot, x, 0);
	      EH(err, "load_elem_dofptr");

	      err = bf_mp_init(pd);
	      EH(err, "bf_mp_init");


	      ielem_type      = ei[pg->imtrx]->ielem_type;  /* element type */
	      ip_total        = elem_info(NQUAD, ielem_type); /* number of  quadrature points */
	      
	      num_local_nodes = ei[pg->imtrx]->num_local_nodes; /* number of local  basis functions */
	      
	      ielem_dim      = ei[pg->imtrx]->ielem_dim; /* physical dimension  of this element */
	      
	      iconnect_ptr    = ei[pg->imtrx]->iconnect_ptr; /* find pointer to beginning */
	      /* of this element's */
	      /* connectivity list */

	      if ( discontinuous )
		{
		  for( k=0; 
		       (fill_inlet_bc = inlet_bc_list[k]) != -8458912 ; k++)
		    {
		      if ( ( fill_inlet_side = 
			  elem_on_ss( exo, BC_Types[fill_inlet_bc].BC_ID, 
				      ielem ) ) ) break;
		    }
		}

	      /* Get element size information for GLS terms in assemble_fill_ext_v */
	      h_elem_siz(hsquared, hhv, dhv_dxnode, pd->e[pg->imtrx][R_MESH1]);

	      
      
/********************************************************************************/
/*                              BLOCK 1                                         */
/*                   START OF VOLUME INTEGRATION LOOP                           */
/*                LOOP OVER THE NUMBER OF QUADRATURE POINTS                     */
/********************************************************************************/
	      /* Loop over all the Volume Quadrature integration points */
	      
	      for (ip = 0; ip < ip_total; ip++)
		{
		  find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */

/* 		  if (eqntype == CORRECT )  find_nodal_stu (ip, ielem_type, &s, &t, &u); */

		  /* locations for current ip */
		  
		  /*
		   * Local element coordinates...
		   */
		  
		  xi[0] = s;
		  xi[1] = t;
		  xi[2] = u;
		  
		  fv->wt = wt = Gq_weight (ip, ielem_type); /* find quadrature weights for */
		  /* current ip */
		  
		  /*
		   * Load up basis function information for ea variable...
		   * Old usage: fill_shape
		   */
		  
		  err = load_basis_functions( xi, bfd);
		  EH( err, "problem from load_basis_functions");
		  
		  /*
		   * This has elemental Jacobian transformation and some 
		   * basic mesh derivatives...
		   * Old usage: calc_Jac, jelly_belly
		   */
		  
		  err = beer_belly();
		  EH( err, "beer_belly");
		  
		  /*
		   * Load up field variable values at this Gauss point.
		   */
		  err = load_fv();
		  EH( err, "load_fv");
		  
		  /*
		   * Here, load in the final part of the necessary basis function
		   * information derivatives in the physical space coordinates...
		   *
		   *			grad(phi_i)
		   *
		   *			grad(phi_i e_a)
		   * where:
		   *		phi_i is the basis function at i-th dof
		   *		e_a   is a unit vector in the coordinate system
		   *
		   * 		grad() operator depends on the coordinate system.
		   */
		  err = load_bf_grad();
		  EH( err, "load_bf_grad");
		  
		  /*
		   * Load up physical space gradients of field variables at this
		   * Gauss point.
		   */
		  err = load_fv_grads();
		  EH( err, "load_fv_grads");	  
		  
		  do_LSA_mods(LSA_VOLUME);
		  
		  switch ( eqntype )
		    {
		    case ADVECT:
		      /* This calls the decoupled version of the extenion
		       * velocity code
		       */
		      if(pd->e[pg->imtrx][R_EXT_VELOCITY] )
			{
			  err = assemble_fill_ext_v(afill, ijaf, rf, delta_t, 
						    theta, node_to_fill, hsquared);
			  EH( err, "assemble_fill_ext_v");
			}
		      else
			{
			  err = assemble_fill(afill, ijaf, rf, delta_t, 
					      theta, node_to_fill);
			}

		      EH( err, "assemble_fill");

		      break;
		    case CORRECT:
		      err = assemble_level_correct ( afill, ijaf, rf, delta_t, 
                                                     theta, node_to_fill);
		      EH( err, "assemble_level_correct");
		      break;
		    case PROJECT:
		      err = assemble_level_project ( afill, ijaf, rf, delta_t,
                                                     theta, node_to_fill);
		      EH( err, "assemble_level_project");
		      break;
		    default:
		      EH(-1,"Unknown equation in fill_matrix.\n");
		    }
		}
	      /* END  for (ip = 0; ip < ip_total; ip++)                               */

/********************************************************************************/
/*                              BLOCK 2                                         */
/*                   START OF SURFACE INTEGRATION LOOP                          */
/*                LOOP OVER THE NUMBER OF SURFACE QUADRATURE POINTS             */
/********************************************************************************/
	      /* Loop over all the surface Quadrature integration points */
	      
	      if(discontinuous)
		{

		  int neighbor;   /* element number of current elements neighbor  */
		  int index, face;
		  double F_inlet = 0.0;

		  for (face = 0; face < ei[pg->imtrx]->num_sides; face++)
		    {
		      index = exo->elem_elem_pntr[ielem] + face;
		      neighbor = exo->elem_elem_list[index];

		      id_side = face + 1; /* the exo/patran side number starts with 1 not zero */

		      if( neighbor == -1 )
			{
			  if ( fill_inlet_side == id_side )
			    {
			      F_inlet = BC_Types[fill_inlet_bc].BC_Data_Float[0];
			    }
			}
		      

		      err = assemble_surface(exo, x, afill, ijaf, rf, delta_t, 
                                             theta, node_to_fill, ielem_type, 
                                             ielem_type_fill, id_side, F_inlet,
					     neighbor,ielem, num_local_nodes);
		      EH( err, "assemble_surface");
		    }
		}

/********************************************************************************/
/*                              BLOCK 3                                         */
/*                   APPLICATION OF FILL GRADIENT (CONTACT ANGLE)               */
/*                   BOUNDARY CONDITIONS                                        */
/********************************************************************************/
	      if ( ( eqntype == CORRECT ) || ( eqntype == PROJECT ) )
		{

		  struct elem_side_bc_struct *elem_side_bc = first_elem_side_BC_array[iel];

		  if ( elem_side_bc != NULL )
		    {		      
		      do
			{
			  int ibc = 0, bc_input_id;
		      
			  while ( ( bc_input_id = (int) elem_side_bc->BC_input_id[ibc++]) != -1 )
			    {
 			      if ( BC_Types[bc_input_id].BC_Name == STRONG_FILL_CA_BC )
                                {
                                  err = apply_strong_fill_ca_bc (afill, ijaf, x, rf, delta_t, theta,
				        node_to_fill, ielem, ielem_type, num_local_nodes, ielem_dim,
				        iconnect_ptr, elem_side_bc, num_total_nodes,
				        BC_Types[bc_input_id].BC_Data_Float[0], exo);
	                          EH( err, " apply_strong_fill_ca_bc");
				}
			    }
			} 
		      while ( (elem_side_bc = elem_side_bc->next_side_bc) != NULL );
		    }
		}

/********************************************************************************/
/*                              BLOCK 4                                         */
/*                   APPLICATION OF LS INLET                                    */
/*                   BOUNDARY CONDITIONS                                        */
/********************************************************************************/

	      if ( ls != NULL )
		{

		  struct elem_side_bc_struct *elem_side_bc = first_elem_side_BC_array[iel];

		  if ( elem_side_bc != NULL )
		    {		      
		      do
			{
			  int ibc = 0, bc_input_id;
		      
			  while ( ( bc_input_id = (int) elem_side_bc->BC_input_id[ibc++]) != -1 )
			    {
 			      if ( BC_Types[bc_input_id].BC_Name == LS_INLET_BC )
				{
				  err = apply_ls_inlet_bc( afill, ijaf, x, rf, node_to_fill, 
							   elem_side_bc, exo );
				  EH(err,"apply_ls_inlet_bc");
				}
			    }
			}
		      while ( ( elem_side_bc = elem_side_bc->next_side_bc ) != NULL );
		    }
		} 
		  
	      
	    } /* END  for (iel = 0; iel < num_internal_elem; iel++)                   */

	}  /* END for (mn =.... */


  return(0);
}
#endif /* not COUPLED_FILL */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int
assemble_surface (Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
		  double x[],	/* global vector containing all unknowns  */
		  double afill[], /* Jacobian matrix for fill equation  */
		  int ijaf[],	/* pointer to nonzeros in Jacobian matrix   */
		  double rf[],	/* rhs vector   */
		  double delta_t, /* current time step size */
		  double theta,	/* parameter to vary time integration from
				 * explicit (theta = 1),implicit (theta = 0) */
		  int node_to_fill[],
		  int ielem_type, /* element type  */
		  int ielem_type_fill, /* element type for fill function */
		  int id_side,	/* id number of current side according to 
				 * EXODUS convention  */
		  double F_inlet, /* source concentrations for inlet boundarys */
		  int neighbor,	/* element neighboring this side */
		  int ielem,	/* current element */
		  int num_local_nodes)   /* number of nodes per element */
{
/*    TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98. */

  int ip, ip1, i, j, I, J;     /* counters */
  int a, idof, jdof, ie, je, ja;    /* more counters */
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];

  int eqn;
  int err;         /* status variable for functions */
  int ip_total, ip_total_fill;
  int nvdofi, nvdofj, ki, kj;
  int found_it, offseti;

  dbl x_neighbor[MAX_SURF_GP][DIM];
  dbl F_neighbor[MAX_SURF_GP], F_n = 1e12;
/*  dbl F_neighbor; */

  double phi_j, phi_i;
  dbl rhs;
  double s, t, u;			/* Gaussian quadrature point locations  */
  double xi[DIM];             /* Local element coordinates of Gauss point. */
  dbl vdotn, vdotn_avg = 1e12;
  dbl vdotn_norm;
  double wt;          /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K) */
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  /***************************************************************************/

  /***************************************************************************/
  /*   START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)          */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                      */
  /***************************************************************************/
  /* Find out the number of surface quadrature points 
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_SURF, ielem_type);
  
  eqn = FILL;

  /* Use one point surface quadrature integration 
     to get the sign of v*n */

  ip_total_fill = elem_info(NQUAD_SURF, ielem_type_fill);

  /*
   * Initialize gauss point coordinate array - (avoid UMR).
   */

  for ( i=0; i<MAX_SURF_GP; i++)
    {
      for ( j=0; j<DIM; j++)
	{
	  x_neighbor[i][j] = 0;
	}
    }

  /* Surface integration over element */
  
  for (ip = 0; ip < ip_total_fill; ip++) 
    {
      /* find the quadrature point locations for current ip */
      
      
      find_surf_st (ip, ielem_type_fill, id_side, pd->Num_Dim, xi, &s, &t, &u);
      
      /* find the quadrature weight for current ip */
      wt = Gq_surf_weight (ip, ielem_type_fill);    
      
      /* ****************************************/
      err = load_basis_functions( xi, bfd);
      EH( err, "problem from load_basis_functions");

      err = beer_belly();
      EH( err, "beer_belly");
            
      /* precalculate variables at  current integration pt.*/
      err = load_fv();
      EH( err, "load_fv");

     /* calculate the determinant of the surface jacobian and the normal to 
      * the surface all at one time 
      */

      err =  get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
	  
      EH( err, "get_side_info");
      
      surface_determinant_and_normal (ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes, 
				      ei[pg->imtrx]->ielem_dim - 1,  
				      id_side,
				      nodes_per_side,
				      local_elem_node_id );

      do_LSA_mods(LSA_SURFACE);

      vdotn_avg = 0.;
      vdotn_norm = 0.;
      for( a=0; a< pd->Num_Dim; a++)
	{
	  vdotn_avg += fv->v[a] * fv->snormal[a];
	  vdotn_norm += fv->v[a]*fv->v[a];
	}
      if (vdotn_avg <  0.)
	{
	  /*  use average F value from center Gauss point for
	   *  all more detailed flux calculation below
	   */
	  if(neighbor != -1  && vdotn_avg*vdotn_avg/vdotn_norm > 1.e-12 )
	    {

	      /*
	       * Fundamentally, this routine fills F_neighbor[] with
	       * values at surface Gauss points from the perspective
	       * of the neighboring element. Those Gauss points will
	       * be in some unusual ordering that is picked apart later
	       * by also load spatial coordinates that are compared
	       * to coordinates of the surface Gauss points in this
	       * element. A better way later...
	       */

	      err =  neighbor_fill(exo, x, ielem, neighbor, F_neighbor,
				   num_local_nodes, nodes_per_side,
				   local_elem_node_id, ielem_type, 
				   ielem_type_fill, node_to_fill, x_neighbor);
	      EH( err, "neighbor_fill");
	    }
	  else
	    {
	      /* if there is no neighbor, set this value to F_inlet
	       * I am assuming we will only get here for inflow
	       * boundaries.
	       */
	      F_n = F_inlet; 
	    }
	}
    }


  /* Surface integration over element */
  if( vdotn_avg <  0.) 
    {
      for (ip = 0; ip < ip_total; ip++) 
	{
	  /* find the quadrature point locations for current ip */
	  
	  
	  find_surf_st (ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
	  
	  /* find the quadrature weight for current ip */
	  wt = Gq_surf_weight (ip, ielem_type);    
	  
	  /* ****************************************/
	  err = load_basis_functions( xi, bfd);
	  EH( err, "problem from load_basis_functions");

          err = beer_belly();
          EH( err, "beer_belly");
            
	  /* precalculate variables at  current integration pt.*/
	  err = load_fv();
	  EH( err, "load_fv");
	  
	  /* calculate the determinant of the surface jacobian and the normal to 
	   * the surface all at one time */
	  
	  err =  get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
	  EH( err, "get_side_info");
	  
	  surface_determinant_and_normal (ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes, 
					  ei[pg->imtrx]->ielem_dim - 1,  
					  id_side,
					  nodes_per_side,
					  local_elem_node_id );

	  do_LSA_mods(LSA_SURFACE);

	  if(neighbor != -1)
	    {
	      found_it = 0; 
	      for (ip1 = 0; ip1 < ip_total && (!found_it); ip1++) 
		{
		  if(  (fabs(fv->x[0]-x_neighbor[ip1][0])<1.e-7)
		       &&(fabs(fv->x[1]-x_neighbor[ip1][1])<1.e-7)
		       &&(fabs(fv->x[2]-x_neighbor[ip1][2])<1.e-7))
		    {
		      F_n = F_neighbor[ip1];
		      found_it = 1;
		    }
		}
	    }

	  vdotn = 0.;
	  for( a=0; a< pd->Num_Dim; a++)
	    {
	      vdotn += fv->v[a] * fv->snormal[a];
	    }
	  vdotn_norm = sqrt(vdotn*vdotn);
	  
	  if (vdotn_norm > 1.e-7 && vdotn < 0 )
	    {
	      /*
	       * Put local contributions into global right-hand side
	       * if it is not a right-hand side variable-it won't get added in (contribution
	       * is zero)
	       */
	      if ( af->Assemble_Residual )
		{
		  rhs = 0.;
		  for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
		    {
		      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
		      node = Nodes[I];
		      nv = node->Nodal_Vars_Info[pg->imtrx];
		      nvdofi = get_nv_ndofs(nv, eqn);
		      offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);
		      for ( ki=0; ki<nvdofi; ki++)
			{
			  /* set rhs to zero for Dirichlet boundary conditions */
			  if (node->DBC[pg->imtrx] &&
			      (node->DBC[pg->imtrx][offseti + ki] != -1) &&
			      Debug_Flag >= 0)
			    {
			      rf[node_to_fill[I]+ki] = 0;
			    }
			  else 
			    {
			      /* check to make sure that unknowns are defined at this node,
				 otherwise don't add anything to this node */
			      idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki;
			      /* also convert from node number to dof number */
			      phi_i = bf[eqn]->phi[idof];
			      
			      rhs = phi_i * ( wt * fv->sdet * 
					      vdotn  * (fv->F - F_n));
			      
			      rf[node_to_fill[I]+ki] -= rhs;
			    }
			}
		    }
		}
	      
	      if ( af->Assemble_Jacobian )
		{
		  for ( i=0; i<ei[pg->imtrx]->num_local_nodes; i++)
		    {
		      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i]; 
		      node = Nodes[I];
		      nv = node->Nodal_Vars_Info[pg->imtrx];
		      nvdofi = get_nv_ndofs(nv, eqn);
		      offseti = get_nodal_unknown_offset(nv, eqn, -2, 0, &vd);	  
		      
		      for ( ki=0; ki<nvdofi; ki++)
			{
			  ie = node_to_fill[I] + ki;
			  idof = ei[pg->imtrx]->ln_to_first_dof[eqn][i] + ki; 
			  phi_i = bf[eqn]->phi[idof];
			  
			  /* derivatives of fill equation wrt to fill variable */	  
			  for( j=0; j< ei[pg->imtrx]->num_local_nodes; j++) 
			    {
			      J = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j]; 
			      nvdofj = Dolphin[pg->imtrx][J][eqn];
			      for ( kj=0; kj<nvdofj; kj++)
				{
				  jdof = ei[pg->imtrx]->ln_to_first_dof[eqn][j] + kj;
				  phi_j = bf[eqn]->phi[jdof];
				  
				  je = node_to_fill[J] + kj;
				  ja     = (ie == je) ? ie : in_list(je, ijaf[ie], ijaf[ie+1], ijaf);
				  EH(ja, "Could not find vbl in sparse matrix.");
				  
				  /* set diagonal to one for Dirichlet boundary conditions */
				  if(node->DBC[pg->imtrx] &&
				     (node->DBC[pg->imtrx][offseti + ki] != -1) &&
				     Debug_Flag >= 0)
				    {
				      afill[ja] = delta(i,j); 
				    }
				  else
				    {
				      afill[ja] -= wt * fv->sdet * 
					vdotn  * phi_j * phi_i;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

  return 0;
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

static int 
neighbor_fill(Exo_DB *exo,
	      dbl x[],
	      int current_elem,
	      int neighbor_elem,
	      dbl F_neighbor[],
	      int num_local_nodes,
	      int nodes_per_side,
	      int local_elem_node_id[],
	      int ielem_type,
	      int ielem_type_fill,
	      int node_to_fill[],
	      dbl x_n[][DIM])
/*   
 *   This function take the current element and side and 
 *   knowing who the neighboring element is, finds the value
 *   of the fill function at the current gauss point location
 *   in the neighboring element.
 *
 *   Currently, this will only work for 1 gauss point per side.

 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 */

{
  int gnn;
  int i,j;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ie;
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int p, dim;     /* counters */
  int ktype;
  int nvdof;
  int status = 0;
  int v;

  dbl phi[MDE], phi_map[MDE];
  dbl s, t, u;
  dbl xi[DIM];

  ielem_shape     = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);
  dim =  pd->Num_Dim;

  for(ip = 0; ip < ip_total; ip++) 
    {
      F_neighbor[ip] = 0.;
      for(p=0; p<DIM; p++)
	{
	  x_n[ip][p] = 0.;
	}
    }

  /* first get global node numbers for current element
     and side */
  iconnect_ptr    = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */

  for (i = 0; i < nodes_per_side; i++)
    {
      id    = local_elem_node_id[i];
      /* load up global node numbers on this side */
      inode[i] = Proc_Elem_Connect[iconnect_ptr + id];
    }

  /* find localside number for neighbor element from
     global node numbers of current element */
  
  id_side = find_id_side (neighbor_elem, nodes_per_side,
			  inode, id_local_elem_coord, exo);

  for (ip = 0; ip < ip_total; ip++) 
    {
      /* find the quadrature point locations for current ip */
      find_surf_st (ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
      
      /* 
       *  we are cheating here and hoping that the "ei[pg->imtrx]->dof[v]" for
       *  the current element will work on the neighbor element that
       *  we are trying to get information for 
       */
      v = FILL;
      /* first load phi for the fill function */      
      for (i = 0; i < ei[pg->imtrx]->dof[v]; i++)  
	{
	  phi[i] = newshape(xi, ielem_type, PSI, 
			    ei[pg->imtrx]->dof_list[v][i], 
			    ielem_shape, pd->i[pg->imtrx][v], i);
	}

      
      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */

      
      ktype = 0;
      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */

      v = pd->ShapeVar;  /* ShapeVar is the variable index where we find the element mapping interpolation */

      for ( i=0; i < ei[pg->imtrx]->dof[v] ; i++)
	{
	  phi_map[i] = newshape(xi, ielem_type, PSI, 
				ei[pg->imtrx]->dof_list[v][i], 
				ielem_shape, pd->i[pg->imtrx][v], i);

	  gnn = Proc_Elem_Connect[iconnect_ptr + i];

	  for ( p=0; p<dim; p++)
	    {
	      x_n[ip][p] +=  Coor[p][gnn] * phi_map[i];  
	    }
	}

      ktype = 0;
      
      for ( i=0; i<num_local_nodes; i++)
	{
	  gnn = Proc_Elem_Connect[iconnect_ptr + i];
	  nvdof = num_varType_at_node(gnn, FILL);
	  for ( j=0; j<nvdof; j++)
	    {
	      ie = Index_Solution(gnn, FILL, ktype, j, -1, pg->imtrx);
	      EH(ie, "Could not find vbl in sparse matrix.");
	      	      
	      F_neighbor[ip] += x[ie] * phi[j];
	      status = 1;
	    }
	}
    }
  
  return(status);
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

int 
get_side_info(const int ielem_type, 
	      const int id_side, /* 1-based EXODUS/PATRAN side number */
	      int *nodes_per_side, 
	      int local_elem_node_id[])
{
 /* TAB certifies that this routine is consistent with the exo/patran side/node number
  * conventions.  11/9/98 */

  switch( ielem_type ){                 /* select element */

  case LINEAR_TRI:                   /* linear triangle */
  case BILINEAR_TRISHELL:
    *nodes_per_side = 2;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 0;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;

   case QUAD_TRI:                   /* quadratic triangle */
   case QUAD6_TRI:                  /* quadratic triangle with 6th order quadrature */
    *nodes_per_side = 3;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      local_elem_node_id[2] = 3;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 4;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 0;
      local_elem_node_id[2] = 5;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;
    
  case BILINEAR_QUAD:                   /* bilinear quadrilateral */
  case C_BILINEAR_QUAD:                   /* bilinear quadrilateral with centroid node */
    // case BILINEAR_SHELL:  /*bilinear quad shell (DSB 7/13) */
    *nodes_per_side = 2;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 3;
      break;
    case 4:
      local_elem_node_id[0] = 3;
      local_elem_node_id[1] = 0;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;

  case S_BIQUAD_QUAD:                    /* biquadratic serendipity quadrilateral */
  case   BIQUAD_QUAD:                    /* biquadratic quadrilateral */
    *nodes_per_side = 3;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      local_elem_node_id[2] = 4;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 5;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 3;
      local_elem_node_id[2] = 6;
      break;
    case 4:
      local_elem_node_id[0] = 3;
      local_elem_node_id[1] = 0;
      local_elem_node_id[2] = 7;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;

  case TRILINEAR_HEX:                    /* trilinear hexahedron */
  case C_TRILINEAR_HEX:
    *nodes_per_side = 4;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      local_elem_node_id[2] = 5;
      local_elem_node_id[3] = 4;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 6;
      local_elem_node_id[3] = 5;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 3;
      local_elem_node_id[2] = 7;
      local_elem_node_id[3] = 6;
      break;
    case 4:
      local_elem_node_id[0] = 3;
      local_elem_node_id[1] = 0;
      local_elem_node_id[2] = 4;
      local_elem_node_id[3] = 7;
      break;
    case 5:
      local_elem_node_id[0] = 3;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 1;
      local_elem_node_id[3] = 0;
      break;
    case 6:
      local_elem_node_id[0] = 4;
      local_elem_node_id[1] = 5;
      local_elem_node_id[2] = 6;
      local_elem_node_id[3] = 7;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;

  case S_TRIQUAD_HEX:                     /* serendipity triquadratic hexahedron */
  case   TRIQUAD_HEX:                     /* triquadratric hexahedron */
    *nodes_per_side = 9;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 8;
      local_elem_node_id[2] = 1;
      local_elem_node_id[3] = 13;
      local_elem_node_id[4] = 5;
      local_elem_node_id[5] = 16;
      local_elem_node_id[6] = 4;
      local_elem_node_id[7] = 12;
      local_elem_node_id[8] = 25;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 9;
      local_elem_node_id[2] = 2;
      local_elem_node_id[3] = 14;
      local_elem_node_id[4] = 6;
      local_elem_node_id[5] = 17;
      local_elem_node_id[6] = 5;
      local_elem_node_id[7] = 13;
      local_elem_node_id[8] = 24;
      break;
    case 3:
      local_elem_node_id[0] = 2;
      local_elem_node_id[1] = 10;
      local_elem_node_id[2] = 3;
      local_elem_node_id[3] = 15;
      local_elem_node_id[4] = 7;
      local_elem_node_id[5] = 18;
      local_elem_node_id[6] = 6;
      local_elem_node_id[7] = 14;
      local_elem_node_id[8] = 26;
      break;
    case 4:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 12;
      local_elem_node_id[2] = 4;
      local_elem_node_id[3] = 19;
      local_elem_node_id[4] = 7;
      local_elem_node_id[5] = 15;
      local_elem_node_id[6] = 3;
      local_elem_node_id[7] = 11;
      local_elem_node_id[8] = 23;
      break;
    case 5:
      local_elem_node_id[0] = 3;
      local_elem_node_id[1] = 10;
      local_elem_node_id[2] = 2;
      local_elem_node_id[3] = 9;
      local_elem_node_id[4] = 1;
      local_elem_node_id[5] = 8;
      local_elem_node_id[6] = 0;
      local_elem_node_id[7] = 11;
      local_elem_node_id[8] = 21;
      break;
    case 6:
      local_elem_node_id[0] = 4;
      local_elem_node_id[1] = 16;
      local_elem_node_id[2] = 5;
      local_elem_node_id[3] = 17;
      local_elem_node_id[4] = 6;
      local_elem_node_id[5] = 18;
      local_elem_node_id[6] = 7;
      local_elem_node_id[7] = 19;
      local_elem_node_id[8] = 22;
      break;
    default:
      EH(-1,"Illegal side number for 2-D element");
      break;
    }
    break;

  case LINEAR_TET:                   /* tri-linear tetrahedron */
    *nodes_per_side = 3;
    switch (id_side) {
    case 1:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 1;
      local_elem_node_id[2] = 3;
      break;
    case 2:
      local_elem_node_id[0] = 1;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 3;
      break;
    case 3:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 3;
      local_elem_node_id[2] = 2;
      break;
    case 4:
      local_elem_node_id[0] = 0;
      local_elem_node_id[1] = 2;
      local_elem_node_id[2] = 1;
      break;
    default:
      EH(-1,"Illegal side number for 3-D tetrahedral element");
      break;
    }
    break;

  default:
    EH(-1, "Unknown or unimplemented element type.");
    break;

  }

  return 1;
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

int
assemble_surface_species (Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
			  double x[],
			  double delta_t, /* current time step size */
			  double theta,	/* parameter to vary time integration from
					 * explicit (theta = 1) to implicit (theta = 0) */
			  int ielem_type, /* element type  */
			  int ielem_type_fill, /* element type for fill function */
			  int id_side,	/* id number of current side according to 
					 * EXODUS convention  */
			  int neighbor,	/* element neighboring this side */
			  int ielem,	/* current element */
			  int num_local_nodes)   /* number of nodes per element */
{
/*    TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98. */

/* LOCAL VARIABLES */
  int ip, ip1, i, j, dim;     /* counters */
  int a, b, p;			/* more counters */
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];

  int eqn;
  int var, var1, pvar;
  int err;         /* status variable for functions */
  int ip_total, ip_total_fill;
  int found_it;
  int species, w;

  dbl conc_neighbor[MAX_SURF_GP][MAX_CONC], c_n[MAX_CONC];
  dbl **x_neighbor;

  dbl advection;

  double phi_j, phi_i;
  dbl rhs;
  double s, t, u;			/* Gaussian quadrature point locations  */
  double xi[DIM];             /* Local element coordinates of Gauss point. */
  dbl vdotn, vdotn_avg = 1e12;
  dbl vdotn_norm;
  double wt;                  /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */
  
  /***************************************************************************/

/* EXTERNAL FUNCTIONS and PROTOTYPES */
  
  /********************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)             */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                           */
  /********************************************************************************/
  /* Find out the number of surface quadrature points 
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_SURF, ielem_type);
  
  eqn = MASS_FRACTION;
  dim =  pd->Num_Dim;

  /* safety measure until MMH/Buyevich model is fixed */
  species = 4;

  /* allocate space for x_neighbor */

  x_neighbor = (double **) array_alloc(2, ip_total, DIM, sizeof(double));

  /* Use one point surface quadrature integration 
     to get the sign of v*n */

  ip_total_fill = elem_info(NQUAD_SURF, ielem_type_fill);

  /* Surface integration over element */
  
  for (ip = 0; ip < ip_total_fill; ip++) 
    {
      /* find the quadrature point locations for current ip */
      
      
      find_surf_st (ip, ielem_type_fill, id_side, pd->Num_Dim, xi, &s, &t, &u);
      
      /* find the quadrature weight for current ip */
      wt = Gq_surf_weight (ip, ielem_type_fill);    
      
      /* ****************************************/
      err = load_basis_functions( xi, bfd);
      EH( err, "problem from load_basis_functions");

      err = beer_belly();
      EH( err, "beer_belly");
            
      /* precalculate variables at  current integration pt.*/
 	      
      err = load_fv();
      EH( err, "load_fv");


     /* calculate the determinant of the surface jacobian and the normal to 
      * the surface all at one time 
      */

      err =  get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
	  
      EH( err, "get_side_info");
      
      surface_determinant_and_normal (ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes, 
				      ei[pg->imtrx]->ielem_dim - 1,  
				      id_side,
				      nodes_per_side,
				      local_elem_node_id );

      do_LSA_mods(LSA_SURFACE);

      vdotn_avg = 0.;
      vdotn_norm = 0.;
      for( a=0; a< pd->Num_Dim; a++)
	{
	  vdotn_avg += fv->v[a] * fv->snormal[a];
	  vdotn_norm += fv->v[a]* fv->v[a];
	}

      if (vdotn_avg <  0.)
	{
/* 	  if(neighbor != -1  && vdotn_avg*vdotn_avg/vdotn_norm > 1.e-12 ) */
	  if(neighbor != -1 )
	    {
	      err =  neighbor_species(exo, x, ielem, neighbor, conc_neighbor,
				      num_local_nodes, nodes_per_side,
				      local_elem_node_id, ielem_type, 
				      ielem_type_fill,  x_neighbor);
	      EH( err, "neighbor_species");
	    }
	  else
	    {
	      /* if there is no neighbor, set this value to one 
	       * I am assuming we will only get here for inflow
	       * boundaries.
	       */

	      /* this needs to generalize though the input deck
	       * for other concentrations on inflow boundaries 
	       */
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  c_n[w] = 1.; 
		}
	    }
	}
    }


  /* Surface integration over element */
  if(vdotn_avg <  0.) 
    {
      for (ip = 0; ip < ip_total; ip++) 
	{
	  /* find the quadrature point locations for current ip */
	  
	  
	  find_surf_st (ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
	  
	  /* find the quadrature weight for current ip */
	  wt = Gq_surf_weight (ip, ielem_type);    
	  
	  /* ****************************************/
	  err = load_basis_functions( xi, bfd);
	  EH( err, "problem from load_basis_functions");

          err = beer_belly();
          EH( err, "beer_belly");
            
	  /* precalculate variables at  current integration pt.*/
	  err = load_fv();
	  EH( err, "load_fv");

	  /* calculate the determinant of the surface jacobian and the normal to 
	   * the surface all at one time */
	  
	  err =  get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
	  EH( err, "get_side_info");
	  
	  surface_determinant_and_normal (ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes, 
					  ei[pg->imtrx]->ielem_dim - 1,  
					  id_side,
					  nodes_per_side,
					  local_elem_node_id );

	  do_LSA_mods(LSA_SURFACE);

	  if(neighbor != -1)
	    {
	      found_it = 0; 
	      for (ip1 = 0; ip1 < ip_total && (!found_it); ip1++) 
		{
		  if(  (fabs(fv->x0[0]-x_neighbor[ip1][0])<1.e-7)
		       &&(fabs(fv->x0[1]-x_neighbor[ip1][1])<1.e-7)
		       &&(fabs(fv->x0[2]-x_neighbor[ip1][2])<1.e-7))
		    {
		      for( w=0; w< pd->Num_Species_Eqn; w++)
			{
			  c_n[w] = conc_neighbor[ip1][w];
			  found_it = 1;
			}
		    }
		}
	    }
	  vdotn = 0.;
	  for( a=0; a< pd->Num_Dim; a++)
	    {
	      vdotn += fv->v[a]* fv->snormal[a];
	    }
	  vdotn_norm = sqrt(vdotn*vdotn);
	  
	  if ((vdotn <  0.) && (vdotn_norm > 1.e-7))
	    {
	      /*
	       * Put local contributions into global right-hand side
	       * if it is not a right-hand side variable-it won't get added in (contribution
	       * is zero)
	       */
	      if ( af->Assemble_Residual )
		{
		  for ( w=0; w<pd->Num_Species_Eqn; w++)
		    {
		      rhs = 0.;
		      for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
			{
			  phi_i = bf[eqn]->phi[i];		      
			  
			  rhs = phi_i * wt * fv->sdet * 
			    vdotn  * (fv->c[w]- c_n[w]);
			  
			  lec->R[MAX_PROB_VAR + w][i] += rhs;
			}
		    }
		}
	      
	      if ( af->Assemble_Jacobian )
		{
		  for ( w=0; w<pd->Num_Species_Eqn; w++)
		    {
		     
		      for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
			{
			  phi_i = bf[eqn]->phi[i];	
			  
			  var = MASS_FRACTION;
			  /* derivatives of conc equation wrt to conc variable */	  
			  for( j=0; j<ei[pg->imtrx]->dof[var]; j++) 
			    {
			      phi_j = bf[var]->phi[j];	
			      advection = wt * fv->sdet * vdotn  * phi_j * phi_i;
			      
			      lec->J[MAX_PROB_VAR + w][MAX_PROB_VAR + w][i][j] += 
				advection;
			    }

			  /*
			   * J_s_v  sensitivity of species equation w.r.t. velocity
			   */			  

			  if(mp->DensityModel == SUSPENSION_PM && w==species)
			    var1 = PVELOCITY1;
			  else
			    var1 = VELOCITY1;
			  for ( b=0; b<dim; b++)
			    {
			      var = var1 + b;
			      if ( pd->v[pg->imtrx][var] )
				{	      
				  pvar = upd->vp[pg->imtrx][var];
				  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];	
				      
				      advection =  phi_i * wt * fv->sdet *  
					phi_j * fv->snormal[b] *
					(fv->c[w] - c_n[w]);
				      
				      lec->J[MAX_PROB_VAR + w][pvar][i][j] += 
					advection;				
				    }
				}
			    }

			  /*
			   * J_s_d  sensitivity of species equation w.r.t. mesh displacement
			   */
			  for ( b=0; b<dim; b++)
			    {
			      var = MESH_DISPLACEMENT1+b;
			      if ( pd->v[pg->imtrx][var] )
				{
				  pvar = upd->vp[pg->imtrx][var];
				  for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];	      
				      
				      advection = 0.;
				      for (p=0; p<pd->Num_Dim; p++)
					{
					  advection += fv->sdet* fv->v[p] 
					    * fv->dsnormal_dx[p][b][j];
					}

				      advection += fv->dsurfdet_dx[b][j] * vdotn;

				      advection *= phi_i * wt * (fv->c[w] - c_n[w]);
				      
				      lec->J[MAX_PROB_VAR + w][pvar][i][j] += 
					advection;
				    }
				}
			    }
			}

		    }
		}
	    }
	}
    }
		      

  safe_free( (void *) x_neighbor);

  return 0;
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

static int 
neighbor_species(Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
		 dbl x[],
		 int current_elem,
		 int neighbor_elem,
		 dbl conc_neighbor[][MAX_CONC],
		 int num_local_nodes,
		 int nodes_per_side,
		 int local_elem_node_id[],
		 int ielem_type,
		 int ielem_type_fill,
		 dbl **x_n)
/*   
*   This function take the current element and side and 
 *   knowing who the neighboring element is, finds the value
 *   of the fill function at the current gauss point location
 *   in the neighboring element.
 *

 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 */

{
  int gnn;
  int i,j;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ie;
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int p, dim;     /* counters */
  int nvdof;
  int status = 0;
  int v;
  int w;

  dbl phi[MDE], phi_map[MDE];
  dbl s, t, u;
  dbl xi[DIM];
  NODE_INFO_STRUCT *node;
  ielem_shape     = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);
  dim =  pd->Num_Dim;

  for(ip = 0; ip < ip_total; ip++) 
    {
      for(p=0; p<DIM; p++)
	{
	  x_n[ip][p] = 0.;
	}
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  conc_neighbor[ip][w] = 0.;
	}
    }

  /* first get global node numbers for current element
     and side */
  iconnect_ptr    = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */

  for (i = 0; i < nodes_per_side; i++)
    {
      id    = local_elem_node_id[i];
      /* load up global node numbers on this side */
      inode[i] = Proc_Elem_Connect[iconnect_ptr + id];
    }

  /* find localside number for neighbor element from
     global node numbers of current element */
  
  id_side = find_id_side (neighbor_elem, nodes_per_side,
			  inode, id_local_elem_coord, exo);

  for (ip = 0; ip < ip_total; ip++) 
    {
      /* find the quadrature point locations for current ip */
      find_surf_st (ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
      
      /* 
       *  we are cheating here and hoping that the "ei[pg->imtrx]->dof[v]" for
       *  the current element will work on the neighbor element that
       *  we are trying to get information for 
       */
      v = MASS_FRACTION;
     
      /* first load phi for the fill function */      
      for (i = 0; i < ei[pg->imtrx]->dof[v]; i++)  
	{
	  phi[i] = newshape(xi, ielem_type, PSI, 
			    ei[pg->imtrx]->dof_list[v][i], 
			    ielem_shape, pd->i[pg->imtrx][v], i);
	}

      v = pd->ShapeVar;
      /* 
       *  we are cheating here and hoping that the "ei[pg->imtrx]->dof[v]" for
       *  the current element will work on the neighbor element that
       *  we are trying to get information for 
       */
      
      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */

      for (i = 0; i < ei[pg->imtrx]->dof[v]; i++)  
	{
	  phi_map[i] = newshape(xi, ielem_type, PSI, 
				ei[pg->imtrx]->dof_list[v][i], 
				ielem_shape, pd->i[pg->imtrx][v], i);
	}

      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */
      for ( i=0; i< ei[pg->imtrx]->dof[v]; i++)
	{
	  gnn = Proc_Elem_Connect[iconnect_ptr + i];

	  for ( p=0; p<dim; p++)
	    {
	      x_n[ip][p] +=  Coor[p][gnn] * phi_map[i];  
	    }
	}
      /*
       * HKM -> I don't know what's going on in this 
       *        next loop. It needs attention
       */
      v = MASS_FRACTION;
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  for ( i=0; i< num_local_nodes; i++)
	    {
	      gnn = Proc_Elem_Connect[iconnect_ptr + i];
	      node = Nodes[gnn];
	      nvdof = get_nv_ndofs_modMF(node->Nodal_Vars_Info[pg->imtrx], v);
	      for (j = 0; j < nvdof; j++)
		{
		  ie = Index_Solution(gnn, v, w, 0, -2, pg->imtrx) + j;
		  EH(ie, "Could not find vbl in sparse matrix.");
		  conc_neighbor[ip][w] += x[ie] * phi[j];
		}
	    }
	}
    }

  status = 1;
  return(status);
}



int
elem_on_ss( Exo_DB *exo,
	    int ss_id,
	    int elem )
{
  int iss,j,elem_index;

  if( (iss = in_list( ss_id, 0, exo->num_side_sets, exo->ss_id) ) == -1 )
    {
      fprintf(stderr, "Error in elem_on_ss.  Cannot find side set id %d",ss_id);
      exit(-1);
    }

  elem_index = exo->ss_elem_index[iss];

  if( ( j = in_list( elem, 0, exo->ss_num_sides[iss], &(exo->ss_elem_list[elem_index]) ) ) == -1 )
    {
      return 0;
    }
  else
    {
      return( exo->ss_side_list[elem_index + j ] );
    }
}


/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void
map_direction( double *d_xi, double *n)
{
  int a,b;
  double mag = 0.0;

  for( a=0; a< pd->Num_Dim; a++)
    {
      d_xi[a] = 0.0;
      
      for( b=0; b<pd->Num_Dim; b++)
	{
	  d_xi[a] += bf[pd->ShapeVar]->B[b][a]*n[b];
	}
      mag += d_xi[a]*d_xi[a];
    }
  
  mag = sqrt( mag );

/*    patches to avoid UMRS for SWIRL and axisymmetric problems */

  for( a=0; a< pd->Num_Dim; a++) 
    if ( fabs( d_xi[a] ) < 1.e-3*mag ) d_xi[a] = 0.0;  

  if( DIM > pd->Num_Dim ) for( a=pd->Num_Dim; a< DIM; a++) d_xi[a] = 0.0;  

}



int 
boundary_curvature( double func[DIM],
		    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		    double *CA )
{
  int j, p, a;
  int var;
  int dim = pd->Num_Dim;

  int status=0;

  double n[DIM], d_n_dF[DIM][MDE], mag_grad_F;

  double cos_CA = 0.0;

  if ( CA != NULL ) cos_CA = cos( M_PIE*(*CA)/180.0);

  for( a=0; a<dim; a++) n[a] = fv->grad_F[a];

  mag_grad_F = normalize_really_simple_vector( n, dim );


  func[0] = 0.0;

  if( CA == NULL )
    {
      for( a=0 ; a<dim; a++)
	func[0] -= n[a]*fv->snormal[a]; 

      if( af->Assemble_Jacobian )
	{

	  for( a=0; a<dim ; a++ )
	    {
	      var = MESH_DISPLACEMENT1 + a;

	      if( pd->v[pg->imtrx][var])
		  {
			  for( j=0 ; j<ei[pg->imtrx]->dof[var] ; j++ )
			  {
				  for( p=0 ; p<dim; p++)
				  {
					  d_func[0][var][j] -= n[p]*fv->dsnormal_dx[p][a][j];
				  }
			  }
		  }
		}

#ifdef COUPLED_FILL
	  var = LS;

	  if( pd->v[pg->imtrx][var] )
	    {
	      /* 
	       * Compute d_n_F
	       */

	      for( a = 0; a<dim; a++ )
		{
		  
		  for ( j = 0 ; j<ei[pg->imtrx]->dof[LS]; j++)
		    {
		      d_n_dF[a][j] = 0.0;
			  
		      for( p=0 ; p<dim ; p++)
			{
			  d_n_dF[a][j] += bf[var]->grad_phi[j][p]*n[p];
			}
			  
		      d_n_dF[a][j] *= -n[a]/mag_grad_F;
		      
		      d_n_dF[a][j] += bf[var]->grad_phi[j][a]/mag_grad_F;
		    }
		}

	      for( j=0 ; j<ei[pg->imtrx]->dof[var]; j++)
		{
		  for( p=0 ; p<dim; p++)
		    {
		      d_func[0][var][j] -= d_n_dF[p][j]*fv->snormal[p];
		    }
		}
	    }
#endif
	} 
    }
  else
    {
      func[0] = cos_CA;
    }

return(status);

}



int
curvature_momentum_source ( double f[DIM],
			    double dfdX[DIM][DIM][MDE],
			    double dfdF[DIM][MDE],
			    double dfdH[DIM][MDE])
{

  int a,b,j;
  int dim = pd->Num_Dim;
  int var;

  memset( dfdH, 0, DIM*MDE*sizeof(double));

#define CONST_CURVE 0

#if CONST_CURVE
  fv->H = 1.;
#endif

  for( a=0; a<dim; a++)
    f[a] += -mp->surface_tension*fv->H*lsi->delta*lsi->normal[a];

      
      /*
       *  dfdX
       */

  
  for( b = 0 ; b < dim; b++ )
    {
      var = MESH_DISPLACEMENT1 + b;
      
      if ( pd->v[pg->imtrx][var] )
	{
	  
	  for( j = 0 ; j<ei[pg->imtrx]->dof[var] ; j++)
	    {
	      for( a=0 ; a<dim ; a++)
		{
		  dfdX[a][b][j] +=  -mp->surface_tension*fv->H*lsi->delta*lsi->d_normal_dmesh[a][b][j];
		}
	    }
	}
    }

      /*
       * dfdH
       */
#if CONST_CURVE
#else
  var = CURVATURE;

  if(  pd->v[pg->imtrx][var] )
    {
      for( j=0 ; j<ei[pg->imtrx]->dof[var] ; j++ )
	{
	  for( a=0; a<dim ; a++ )
	    {
	      dfdH[a][j] = -mp->surface_tension * bf[var]->phi[j] * lsi->delta*lsi->normal[a];
	    }
	}
    }
#endif

#ifdef COUPLED_FILL

  var = LS;

  if( pd->v[pg->imtrx][var] )
    {
      for( j=0 ; j<ei[pg->imtrx]->dof[var] ; j++ )
	{
	  
	  for( a=0; a<dim; a++ )
	    {
	      dfdF[a][j] += -mp->surface_tension * fv->H *( lsi->d_delta_dF[j]*lsi->normal[a] +
							    lsi->delta* lsi->d_normal_dF[a][j] );
	    }
	}
    }
#endif

return(0);

}

/* 
 * Phase functions kernels 
 *  For convenience the phase function kernels are located in this file.  They should be moved elsewhere as the capability warrants it.
 */

int
assemble_phase_function ( double time_value,
			  double tt,
			  double dt,
			  double xi[DIM],
			  Exo_DB *exo)
{
  int i, j;
  int a, b, p;
  int eqn, var, peqn, pvar, ledof;
  int dim, wim;
  int status = -1;
  int Fill_Wt_Fcn = tran->Fill_Weight_Fcn;

  double det_J, h3, wt, phi_i, phi_j, wfcn = 0, *grad_phi_j, *grad_phi_i;
  double grad_II_phi_i[DIM], grad_II_phi_j[DIM];
  double d_wfcn_du;
  double pf_dot,  *grad_pf;
  double grad_II_pf[DIM];
  double *v, *v_old;
  dbl *xx;                              /* Nodal coordinates. */
  dbl *x_old;                           /* Old xx[]. */
  dbl x_dot[DIM];                       /* Time derivative of the mesh displacements. */
  dbl *x_dot_old;                 /* Old x_dot[]. */

  
  dbl v_dot_DF;                 /* v.grad(F) */
  dbl dtinv;                    /* = 1 / dt */
  double v_rel[DIM], v_rel_old[DIM];

  double v_dot_grad_pf, v_dot_Dphi[MDE], vc_dot_Dphi[MDE];

  dbl zero[3] = {0.0, 0.0, 0.0};        /* An array of zeros, for convienience. */

  /* See get_supg_stuff() for a better description of these variables */
  dbl supg_term, vcent[DIM], d_vcent_du[DIM][MDE][DIM], d_supg_term_du[MDE][DIM], d_supg_term_dx[MDE][DIM];


  double rhs = 0, mass=0.0, advection=0.0, tmp = 0.0;
  /* double pf, d_n_gradT_dT[MDE]; */
  double tau_gls,vmag_old;
  struct Level_Set_Data *ls_old;

  eqn = PHASE1;
  
  if ( !pd->e[pg->imtrx][eqn] || pfd == NULL )
    {
      return(status);
    }
  
  memset(grad_II_phi_i, 0, sizeof(double)*DIM);

  /*
   * Calculate lubrication velocity for direct integration
   * On Phase fields we associate the level set with LUBP_2 field.   
   */
  int lubon = 0;
  if ( pd->e[pg->imtrx][R_LUBP_2] ) {
     if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ||
    	 tran->Fill_Weight_Fcn == FILL_WEIGHT_TG) {
   // if ( tran->Fill_Weight_Fcn == FILL_WEIGHT_G ) {
      lubon = 1;
    } else {
      WH(-1,"\n Multiphase lubrication should be run with Galerkin weighting to \n take advantage of direct velocity calculations.  \n Talk to SAR.");
      lubon = 0;
    }
  }

  dim = pd->Num_Dim;
  wim = dim;

  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN ||
      pd->CoordinateSystem == CARTESIAN_2pt5D)
    wim = wim+1;

  wt = fv->wt;

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			/* Differential volume element (scales). */

  dtinv           = 1.0 / dt;

  int *n_dof = NULL;
  int dof_map[MDE];
  if ( pd->e[pg->imtrx][R_LUBP_2] ) {

    /* Set up shells */
    n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Calculate velocity */
    calculate_lub_q_v(R_LUBP_2, time_value, dt, xi, exo);
    calculate_lub_q_v_old (R_LUBP_2, tran->time_value_old, tran->delta_t_old, xi, exo );

    /* Set up weights */
    wt    = fv->wt;
    h3    = fv->h3;
    det_J = fv->sdet;

  }

  grad_pf = fv->grad_pF[0];
  if (lubon) Inn(grad_pf, grad_II_pf);

  /* Use pointers unless we need to do algebra. */
  if (lubon) {
    v     = LubAux->v_avg;
    v_old = LubAux_old->v_avg;
    xx    = fv->x;
    x_old = fv_old->x;
  } else {
    v     = fv->v;
    v_old = fv_old->v;
    xx    = fv->x;
    x_old = fv_old->x;
  }

  v_dot_DF = 0.0;
  if ( pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1] )
    {
      x_dot_old = fv_dot_old->x;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	  if (lubon) x_dot[a] = (1+2*tt)/dt * (xx[a] - x_old[a]);
	  v_rel[a]     = v[a] - x_dot[a];
	  v_rel_old[a] = v_old[a] - x_dot_old[a];
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_pf[a];       /* v.gradII(pF) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_pf[a];       /* v.grad(pF) */
            }
	}
    }
  else if  (pd->TimeIntegration != STEADY &&
	    pd->etm[pg->imtrx][R_SOLID1][(LOG2_MASS)] &&
	    pd->MeshMotion == TOTAL_ALE)    /*This is the Eulerian solid-mech case */
    {
      xx	 = fv->d_rs;
      x_old	 = fv_old->d_rs;
      x_dot_old =  fv_dot_old->d_rs;
      for(a=0; a < VIM; a++)
	{
	  x_dot[a]     = (1.+2.*tt) * (xx[a] - x_old[a]) * dtinv - 2. * tt * x_dot_old[a];
	}
      for(a=0; a < VIM; a++)
	{
	  v_rel[a] = x_dot[a];
	  v_rel_old[a] = x_dot_old[a];
	  for(b=0; b < VIM; b++)
	    {
	      v_rel[a] -= x_dot[b]*fv->grad_d_rs[b][a];
	      v_rel_old[a] -= x_dot_old[b]*fv_old->grad_d_rs[b][a];
	    }
	}
      for ( a=0; a < VIM; a++ )
	{
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_pf[a];       /* v.gradII(pF) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_pf[a];       /* v.grad(pF) */
            }
	}

    }
  else
    {
      x_dot_old = zero;
      for ( a=0; a < VIM; a++ )
	{
	  x_dot[a]     = 0.0;
	  v_rel[a]     = v[a];
	  v_rel_old[a] = v_old[a];
          if (lubon)
            {
	     v_dot_DF    += v_rel[a] * grad_II_pf[a];       /* v.gradII(pF) */
            }
          else
            {
	     v_dot_DF    += v_rel[a] * grad_pf[a];       /* v.grad(pF) */
            }
	}
    }

  if ( Fill_Wt_Fcn == FILL_WEIGHT_SUPG )
    {
      memset(vcent,          0, sizeof(double)*DIM);
      memset(d_vcent_du,     0, sizeof(double)*DIM*MDE*DIM);
      memset(d_supg_term_du, 0, sizeof(double)*MDE*DIM);
      memset(d_supg_term_dx, 0, sizeof(double)*MDE*DIM);
      memset(vc_dot_Dphi, 0, sizeof(double)*MDE);
      supg_term = 0.;
      get_supg_stuff(&supg_term, vcent, d_vcent_du, d_supg_term_du, d_supg_term_dx, pd->e[pg->imtrx][R_MESH1]);
    }

  memset( v_dot_Dphi, 0, sizeof(double)*MDE);
  for ( i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
    {
      grad_phi_i = bf[eqn]->grad_phi[i];

      if (lubon)
        {
         for( a=0; a<VIM; a++)
	    {
             v_dot_Dphi[i] += v_rel[a] * grad_II_phi_i[a];
            }
        }
      else
        {
         for( a=0; a<VIM; a++)
	    {
             v_dot_Dphi[i] += v_rel[a] * grad_phi_i[a];

	     if ( Fill_Wt_Fcn == FILL_WEIGHT_SUPG )
		vc_dot_Dphi[i] += vcent[a] * grad_phi_i[a];
	    }
        }
  }

  vmag_old = 0.;
  for ( a=0; a<dim; a++)
    {
      vmag_old += v[a] * v_old[a];
    }
  vmag_old = sqrt( vmag_old );
/*  tau_gls = 1./sqrt((2./dt)*(2./dt) + (2.*vmag_old/h_elem)*(2.*vmag_old/h_elem));  */
  tau_gls = 1./sqrt((2./dt)*(2./dt) );
  /*
   * Residuals_________________________________________________________________
   */

      /*
       * Assemble each component "a" of the phase function vector equations...
       */
  for ( a=0; a<pfd->num_phase_funcs; a++)
    {

      eqn = R_PHASE1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      /* pf = fv->pF[a]; */
      pf_dot = fv_dot->pF[a];
      grad_pf = fv->grad_pF[a];
      if (lubon) Inn(grad_pf, grad_II_pf);
      ls_old = ls;
      ls = pfd->ls[a];
      load_lsi( 0.0 );
      load_lsi_derivs();

      for( p=0, v_dot_grad_pf = 0.0;  p<wim; p++)
	{
          if (lubon)
            {
	     v_dot_grad_pf += v_rel[p]*grad_II_pf[p];
            }
          else
            {
	     v_dot_grad_pf += v_rel[p]*grad_pf[p];
            }
	}

      if ( af->Assemble_Residual )
	{

	  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
	    {
	      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

	      phi_i = bf[eqn]->phi[i];

	      if (ei[pg->imtrx]->active_interp_ledof[ledof]) 
		{
		  switch( Fill_Wt_Fcn )
		    {
		    case FILL_WEIGHT_TG:
		      rhs = pf_dot*phi_i;

		      for( p=0; p<dim; p++)
			{
                          if (lubon)
                            {
			     rhs += phi_i*0.5*( v_rel[p] + v_rel_old[p] )*grad_II_pf[p];
                            }
                          else
                            {
			     rhs += phi_i*0.5*( v_rel[p] + v_rel_old[p] )*grad_pf[p];
                            }
			}

		      rhs += v_dot_Dphi[i] * v_dot_grad_pf * dt * 0.5;
		      break;


		    case FILL_WEIGHT_G: /* Vanilla Galerkin weight on advection equation */

		      rhs = ( pf_dot + v_dot_grad_pf ) * phi_i;
		      break;

		    case FILL_WEIGHT_SUPG:

		      rhs = ( pf_dot + v_dot_grad_pf ) * ( vc_dot_Dphi[i] * supg_term + phi_i );

		      break;
		    case FILL_WEIGHT_EXPLICIT: /* Vanilla Galerkin weight on advection equation */
	       wfcn = phi_i;
               for (p = 0; p < dim; p++) wfcn += tau_gls * v_rel_old[p] * bf[eqn]->grad_phi[i][p];

		      rhs = ( pf_dot + v_dot_grad_pf ) * wfcn;
		      break;

		    default:
		      EH(-1, "That Fill Weight Function type not currently supported for phase function\n");
		      break;
		    }

		  lec->R[peqn][i] += rhs*wt*det_J*h3;

		}
	    }
	}

      /*
       * Jacobian terms...
       */

      if ( af->Assemble_Jacobian )
	{


	  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
	    {
	      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

	      phi_i = bf[eqn]->phi[i];
	      grad_phi_i = bf[eqn]->grad_phi[i];

	      if( Fill_Wt_Fcn == FILL_WEIGHT_G )
		{
		  wfcn = phi_i;
		}
	      else if ( Fill_Wt_Fcn == FILL_WEIGHT_SUPG )
		{
		  wfcn = 0.0;

		  for( p=0; p< dim ; p++ )
		    {
		      wfcn += vcent[p] * grad_phi_i[p];
		    }

		  wfcn *= supg_term;
		  wfcn += phi_i;
		}


	      if (ei[pg->imtrx]->active_interp_ledof[ledof]) 
		{

		  /* J_pf_pf */

		  var = PHASE1 + a;

		  if ( pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{

			  phi_j = bf[var]->phi[j];
			  grad_phi_j = bf[var]->grad_phi[j];

                          if(lubon) Inn(grad_phi_j, grad_II_phi_j);

			  switch( Fill_Wt_Fcn )
			    {
			    case FILL_WEIGHT_TG:
			      tmp = phi_i * phi_j * ( 1. + 2.0*tt) / dt;

			      for( p=0; p<VIM; p++ )
				{
                                  if (lubon)
                                    {
                                     tmp += 0.5 * LubAux->dv_avg_df[p][j] * grad_II_pf[p] * phi_i;
				     tmp += 0.5 * (v_rel[p] + v_rel_old[p] ) *grad_II_phi_j[p] * phi_i;
                                    }
                                  else
                                    {
				     tmp += 0.5 * (v_rel[p] + v_rel_old[p] ) *grad_phi_j[p] * phi_i;
                                    }
				}

			      tmp += v_dot_Dphi[j] * v_dot_Dphi[i] * dt * 0.5;
			      tmp *= wt*h3*det_J;

			      break;

			    case FILL_WEIGHT_G: /* Vanilla Galerkin weight on advection equation */

			     tmp = ( phi_j * (1. + 2.*tt) * dtinv ) * wfcn;
			     tmp += v_dot_Dphi[j] * wfcn;
			     if (lubon) {
			       for (p = 0; p < dim; p++) tmp += LubAux->dv_avg_df[p][j] * grad_pf[p] * wfcn;
			     }
			     tmp *= wt*h3*det_J;

				break;
			    case FILL_WEIGHT_SUPG:

			      tmp = ( phi_j*(1 + 2.0*tt)/dt + v_dot_Dphi[j] )*wt*h3*det_J*wfcn;


			      break;
			    case FILL_WEIGHT_EXPLICIT: /* Vanilla Galerkin weight on advection equation */
			     tmp = phi_j;
				break;

			    default:
			      EH(-1, "That Fill Weight Function type not currently supported for phase function\n");
			      break;
				}

			  lec->J[peqn][pvar][i][j] += tmp;

			} /* for ( j = 0 ... */
		    } /* if ( pd->v[pg->imtrx][var] ) */
		  /* T */

		  var = TEMPERATURE;

		  if ( pd->v[pg->imtrx][var] )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
			{

			  phi_j = bf[var]->phi[j];
			  grad_phi_j = bf[var]->grad_phi[j];

			  switch( Fill_Wt_Fcn )
			    {
			    case FILL_WEIGHT_EXPLICIT:
			     tmp = phi_j;
			     break;
			     //case 99: PRS:Here for tracking melting fronts?
			     // tmp = phi_i * ( pf*d_n_gradT_dT[j] - phi_j);
			     // tmp *= wt*h3*det_J;
			     // break;
			    default:
			      EH(-1, "That Fill Weight Function type not currently supported for phase function\n");
			      break;
			     }

			  lec->J[peqn][pvar][i][j] += tmp;

			} /* for ( j = 0 ... */
		    } /* if ( pd->v[pg->imtrx][var] ) */

		  /*
		   *   J_pf_v 
		   */

		  for( b=0; b<VIM; b++)
		    {
		      var = VELOCITY1 +b;

		      if( pd->v[pg->imtrx][var] )
			{
			  pvar = upd->vp[pg->imtrx][var];

			  for ( j=0; j < ei[pg->imtrx]->dof[var]; j++ )
			    {
			      phi_j = bf[var]->phi[j];

			      switch(Fill_Wt_Fcn)
				{
				case FILL_WEIGHT_TG:

				  tmp = grad_phi_i[b] * v_dot_grad_pf;
				  tmp += v_dot_Dphi[i] * grad_pf[b];
				  tmp *= phi_j * dt * 0.5;
				  tmp += 0.5 * phi_j * phi_i * grad_pf[b];

				  break;
				case FILL_WEIGHT_G:
				  tmp = phi_j * grad_pf[b] * phi_i;
				  if (lubon) tmp = 0.0;

				  break;

				case FILL_WEIGHT_SUPG:
				  d_wfcn_du = d_supg_term_du[j][b] * vc_dot_Dphi[i];

				  for( p=0 ; p<dim ; p++ )
				    {
				      d_wfcn_du += supg_term * d_vcent_du[p][j][b] * grad_phi_i[p];
				    }

				  tmp = ( pf_dot + v_dot_grad_pf ) * d_wfcn_du;

				  tmp += phi_j * grad_pf[b] * wfcn;
				  break;

				default:
				  break;
				}

			      lec->J[peqn][pvar][i][j] += tmp*wt*det_J*h3;
			    }
			}
		    }

		  /*************************************************************
		   *
		   * Derivatives of fill equation w.r.t. to LUBP_2 variable
		   *
		   *************************************************************/

		  var = LUBP_2;
		  if (pd->v[pg->imtrx][var] && lubon )
		    {
		      pvar = upd->vp[pg->imtrx][var];
		      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
			{
			  phi_j = bf[eqn]->phi[j];
			  grad_phi_j = bf[eqn]->grad_phi[j];

			  dbl grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
			  ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);

			  mass = 0.0;
			  advection = 0.0;

                          switch(Fill_Wt_Fcn)
                           {

			    case FILL_WEIGHT_TG:
			      for (p = 0; p < dim; p++) {
                                advection += 0.5 * LubAux->dv_avg_dp1[p][j] * grad_II_pf[p] * grad_II_phi_j[p] * phi_i;
                                advection += 0.5 * LubAux->dv_avg_dp2[p][j] * grad_II_pf[p] * phi_j * phi_i;

                                advection += 0.5 * LubAux->dv_avg_dp1[p][j] * grad_II_phi_i[p] * grad_II_phi_j[p]
                                                 * v_dot_grad_pf * dt * 0.5;
                                advection += 0.5 * LubAux->dv_avg_dp1[p][j] * grad_II_phi_i[p] * phi_j
                                                 * v_dot_grad_pf * dt * 0.5;

                                advection += 0.5 * LubAux->dv_avg_dp1[p][j] * grad_II_pf[p] * grad_II_phi_j[p]
                                                 * v_dot_Dphi[i] * dt * 0.5;
                                advection += 0.5 * LubAux->dv_avg_dp1[p][j] * grad_II_pf[p] * phi_j
                                                 * v_dot_Dphi[i] * dt * 0.5;

                              }

                              break;

			    case FILL_WEIGHT_G:
			      for (p = 0; p < dim; p++) {
			        advection += LubAux->dv_avg_dp1[p][j] * grad_II_pf[p] * wfcn * grad_II_phi_j[p];
			        advection += LubAux->dv_avg_dp2[p][j] * grad_II_pf[p] * wfcn * phi_j;
			      }

			      break;
                           }

			  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
			  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

			  lec->J[peqn][pvar][i][j] += (mass + advection) * wt * h3 * det_J;

			} /* for: LUBP DoFs */

		    } /* if: LUBP_2 exists */

	  /*************************************************************
	   *
	   * Derivatives of fill equation w.r.t. to SHELL_LUB_CURV variable
	   *
	   *************************************************************/

          var = SHELL_LUB_CURV;
	  if ( pd->v[pg->imtrx][var] && lubon )
	    {
	      pvar = upd->vp[pg->imtrx][var];
	      for( j=0; j < ei[pg->imtrx]->dof[var]; j++) 
		{
		  phi_j = bf[eqn]->phi[j];

		  mass = 0.0;
		  advection = 0.0;

                  switch( Fill_Wt_Fcn )
                    {

                     case FILL_WEIGHT_TG:
		       for (p = 0; p < dim; p++) {
                         advection += 0.5 * LubAux->dv_avg_dk[p][j] * grad_II_pf[p] * phi_i * phi_j;
                         advection += 0.5 * LubAux->dv_avg_dk[p][j] * grad_II_phi_i[p] * phi_j 
                                          * v_dot_grad_pf * dt * 0.5;
                         advection += 0.5 * LubAux->dv_avg_dk[p][j] * grad_II_pf[p] * phi_j  
                                          * v_dot_Dphi[i] * dt * 0.5; 
                       }

                       break;

                     case FILL_WEIGHT_G:
		       for (p = 0; p < dim; p++)
		         advection += LubAux->dv_avg_dk[p][j] * grad_II_pf[p] * phi_i * phi_j;

                       break;
                    }

		  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

		  lec->J[peqn][pvar][i][j] += (mass + advection) * wt * h3 * det_J;

		} /* for: SHELL_LUB_CURV DoFs */

	    } /* if: SHELL_LUB_CURV exisits */

		} /* if( active...) */
	    } /* for(i : */
	} /* af->Jacobian */
      ls = ls_old;
    }

   if ( pd->e[pg->imtrx][R_LUBP_2] )  safe_free((void *) n_dof);
  return (0);
}


int
assemble_pf_constraint( double delta_t,
						double *pf_constraint,
						double lambda,
						double *ptr_d_pf_lm,
						double *ptr_d_lm_pf )
{
	int i,j, ledof;
	int var, pvar, ie;
	double wt, det_J, h3;
	double sum, sum_old, phi_j;
	
	struct Level_Set_Interface pf_lsi[MAX_PHASE_FUNC], pf_lsi_old[MAX_PHASE_FUNC];
	struct Level_Set_Data *ls_save = ls;
	struct Level_Set_Interface *lsi_save = lsi;
	
	double d_pf_lm[MAX_PHASE_FUNC][MDE];
	double d_lm_pf[MAX_PHASE_FUNC][MDE];
	
	
	wt = fv->wt;
	
	det_J = bf[R_PHASE1]->detJ;		/* Really, ought to be mesh eqn. */
	
	h3 = fv->h3;			/* Differential volume element (scales). */
	
	
	for( i=0, sum_old=sum=0.0; i<pfd->num_phase_funcs; i++)
	{
		ls = pfd->ls[i];
		lsi = &(pf_lsi[i]);
		
		load_lsi( ls->Length_Scale );
		load_lsi_derivs( );
		
		sum += pf_lsi[i].H ;
		
		if(lsi->near)
		{
			pf_lsi_old[i].dH = 0.5 * (1.0/lsi->alpha) * (1. + cos(M_PIE * fv_old->pF[i] / lsi->alpha));
			sum_old += 0.5 * (1. + fv_old->pF[i] / lsi->alpha + sin(M_PIE * fv_old->pF[i] / lsi->alpha) / M_PIE);
		}
		else
		{
			sum_old += ( fv_old->pF[i] < 0. ) ? 0.0 : 1.0 ;
			pf_lsi_old[i].dH = 0.0;
		}
	}
	
	sum -= 1.0;
	sum_old -= 1.0;
	
	*pf_constraint += 0.5*wt*sum*sum*h3*det_J;

	memset(d_pf_lm, 0, sizeof(double)*MAX_PHASE_FUNC*MDE );
	memset(d_lm_pf, 0, sizeof(double)*MAX_PHASE_FUNC*MDE );
	
	
	for( i=0; i<pfd->num_phase_funcs; i++)
	{
		var = PHASE1 + i;
		
		pvar = upd->vp[pg->imtrx][var];
		
		if(pd->v[pg->imtrx][var] )
		{
			for( j=0; j<ei[pg->imtrx]->dof[var] ; j++ )
			{
				phi_j = bf[var]->phi[j];
				
				d_pf_lm[i][j] = wt*phi_j*pf_lsi_old[i].dH*sum_old*h3*det_J/delta_t;
				d_lm_pf[i][j] = wt*phi_j*pf_lsi[i].dH*sum*h3*det_J;
				
				lec->R[pvar][j] += lambda*d_pf_lm[i][j];
				
				ledof = ei[pg->imtrx]->lvdof_to_ledof[var][j];
				if( ei[pg->imtrx]->owned_ledof[ledof] )
				{
					ie = ei[pg->imtrx]->ieqn_ledof[ledof];
					ptr_d_pf_lm[ie] += d_pf_lm[i][j];
					ptr_d_lm_pf[ie] += d_lm_pf[i][j];
				}
				
			}
		}
	}
	
	ls = ls_save;
	lsi = lsi_save;
	
	return(TRUE);
}

void assemble_ls_mass_lumped_penalty(
    struct LS_Mass_Lumped_Penalty *mass_lumped_penalty, int ip_total,
    int ielem_type, PG_DATA *pg_data) {
  dbl s, t, u;
  dbl MM[MDE];
  dbl rhs[MDE];
  dbl d_rhs[MDE][MDE];
  dbl rhs_old[MDE];

  for (int i = 0; i < MDE; i++) {
    MM[i] = 0;
    rhs[i] = 0;
    rhs_old[i] = 0;
    for (int k = 0; k < MDE; k++) {
      d_rhs[i][k] = 0;
    }
  }
  for (int ip = 0; ip < ip_total; ip++) {
    dbl xi[DIM];
    int err;
    find_stu(ip, ielem_type, &s, &t, &u); /* find quadrature point */
    dbl wt = Gq_weight(ip, ielem_type);   /* find quadrature weights for */

    xi[0] = s;
    xi[1] = t;
    xi[2] = u;

    fv->wt = wt;
    err = load_basis_functions(xi, bfd);
    EH(err, "problem from load_basis_functions");

    err = beer_belly();
    EH(err, "beer_belly");
    if (neg_elem_volume)
      return;
    if (zero_detJ)
      return;

    err = load_fv();
    EH(err, "load_fv");

    err = load_bf_grad();
    EH(err, "load_bf_grad");

    err = load_fv_grads();
    EH(err, "load_fv_grads");

    int eqn = R_FILL;
    int dim = pd->Num_Dim;

    double F_norm = 0;
    for (int i = 0; i < dim; i++) {
      F_norm += fv->grad_F[i] * fv->grad_F[i];
    }
    F_norm = sqrt(F_norm);
    double inv_F_norm = 1 / F_norm;

    double F_norm_old = 0;
    for (int i = 0; i < dim; i++) {
      F_norm_old += fv_old->grad_F[i] * fv_old->grad_F[i];
    }
    F_norm_old = sqrt(F_norm_old);

    double d_F_norm[MDE];
    for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
      d_F_norm[j] = 0;
      for (int i = 0; i < dim; i++) {
        d_F_norm[j] += fv->grad_F[i] * bf[eqn]->grad_phi[j][i] * inv_F_norm;
      }
    }

    double h_elem = 0.;
    for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      h_elem += pg_data->hsquared[a];
    }
    /* This is the size of the element */
    h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));

    double vnorm = 0;
    for (int a = 0; a < dim; a++) {
      vnorm += fv->v[a] * fv->v[a];
    }
    vnorm = sqrt(vnorm);

    double lambda = 0.8 * h_elem * h_elem * vnorm * 0.5;

    for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      rhs[i] += fv->wt * bf[eqn]->detJ * fv->h3 * lambda * bf[eqn]->phi[i] *
                (F_norm - 1);

      rhs_old[i] += fv->wt * bf[eqn]->detJ * fv->h3 * lambda * bf[eqn]->phi[i] *
                (F_norm_old - 1);
      //        rhs[i] += fv->wt * bf[R_FILL]->detJ * bf[R_FILL]->phi[i] *
      //        (eikonal_norm - 1);
      //      rhs[i] +=
      //          fv->wt * bf[R_FILL]->detJ * bf[R_FILL]->phi[i] * lambda *
      //          (sqrt(fv->grad_F[0] * fv->grad_F[0] + fv->grad_F[1] *
      //          fv->grad_F[1]) -
      //           1);
      for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
        d_rhs[i][k] += fv->wt * bf[eqn]->detJ * fv->h3 * lambda *
                       bf[eqn]->phi[i] * (d_F_norm[k]);
        //          d_rhs[i][k] += fv->wt * bf[R_FILL]->detJ *
        //          bf[R_FILL]->phi[i] * d_eikonal_norm[k];
        //        d_rhs[i][k] += fv->wt * bf[R_FILL]->detJ * bf[R_FILL]->phi[i]
        //        * lambda *
        //                       (2 * bf[R_FILL]->grad_phi[k][0] * fv->grad_F[0]
        //                       +
        //                        2 * bf[R_FILL]->grad_phi[k][1] *
        //                        fv->grad_F[1]) /
        //                       (2 * sqrt(fv->grad_F[0] * fv->grad_F[0] +
        //                                 fv->grad_F[1] * fv->grad_F[1]));
        //      }
        for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          MM[i] += fv->wt * bf[R_FILL]->detJ * bf[R_FILL]->phi[i] *
                   bf[R_FILL]->phi[j] * fv->h3;
        }
      }
    }

    for (int i = 0; i < ei[pg->imtrx]->dof[R_FILL]; i++) {
      mass_lumped_penalty->penalty[i] = rhs[i] / MM[i];
      mass_lumped_penalty->penalty_old[i] = rhs_old[i] / MM[i];
      for (int k = 0; k < ei[pg->imtrx]->dof[R_FILL]; k++) {
        mass_lumped_penalty->d_penalty[i][k] = d_rhs[i][k] / MM[i];
      }
    }
  }
}


/***************************************************************************************/
/* end of mm_fill_fill.c */

