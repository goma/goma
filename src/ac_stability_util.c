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
 *$Id: ac_stability_util.c,v 5.6 2010-04-07 22:27:00 prschun Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: ac_stability_util.c,v 5.6 2010-04-07 22:27:00 prschun Exp $";
#endif

#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "el_elm.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "dpi.h"
#include "exo_struct.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "wr_dpi.h"
#include "wr_exo.h"

#define GOMA_AC_STABILITY_UTIL_C

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
*									      *
*									      *
*  void do_LSA_mods(int):  Called after basic bf/fv assembly sequence from    *
*			   several places (e.g. matrix_fill).		      *
*			   Int argument mod_type identifies whether assembly  *
*			   is for volume, surface, or edge integrals.	      *
*									      *
*  Function do_LSA_mods then calls (as needed):				      *
*	void modify_basis_and_weight_functions_for_LSA_3D_of_2D(void)	      *
*	void modify_bf_mesh_derivs_for_LSA_3D_of_2D(void)		      *
*	void modify_fv_mesh_derivs_for_LSA_3D_of_2D(void)		      *
*	void modify_normal_vector_for_LSA_3D_of_2D(void)		      *
*									      *
*  int create_eigen_outfiles(void ***):  Called from LOCA interface to	      *
*					 create the ExodusII files needed     *
*                                        for ARPACK eigenvector output.       *
*									      *
*  void close_eigen_outfiles(void ***):  Closes the aforesaid files.          *
*									      *
*******************************************************************************/

/*
 * EDW 9/28/2001:
 * This is a wrapper around all modifications to basis/weight functions,
 * gradient mesh derivatives, and vector mesh derivatives needed
 * to do 3D stability of 2D flows with or without mesh displacement.
 */

void
do_LSA_mods(int mod_type)
{

/* Proceed only if this is a 3D of 2D stability problem */
  if (Linear_Stability != LSA_3D_OF_2D
   && Linear_Stability != LSA_3D_OF_2D_SAVE) return;

/* All cases: modify all basis and weight functions */
  modify_basis_and_weight_functions_for_LSA_3D_of_2D();

/* If not using a deforming mesh:  that's all, folks! */
  if (!pd->e[pg->imtrx][R_MESH1]) return;

/* First, do modifications for field variables and basis/weight functions */
      modify_fv_mesh_derivs_for_LSA_3D_of_2D();
      modify_bf_mesh_derivs_for_LSA_3D_of_2D();

/* Do modifications for surface integral assembly */
  if (mod_type == LSA_SURFACE)
    {
      modify_normal_vector_for_LSA_3D_of_2D();
    }

/* Do modifications for edge integral assembly */
/* EDW 10/4/2001   WARNING: This has not yet been tested! */
  else if (mod_type == LSA_EDGE)
    {
      printf("WARNING: LSA edge vectors not yet tested!\n");
      modify_normal_vector_for_LSA_3D_of_2D();
    }

}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/* This function will modify the basis and weight functions for 3D
 * stability of a 2D flow.  We want to perform this modification AFTER
 * the steady-state variables have been defined in terms of the
 * regular basis functions.  The regular steady-state solution gets
 * constructed by the sequence (beer_belly(), load_fv(),
 * load_bf_grad(), load_bf_mesh_derivs(), load_fv_grads(),
 * load_fv_mesh_derivs()).  This applies to both internal and boundary
 * nodes.
 *
 * This code must be called at the end of the above sequence.  From
 * then on, any direct reference to phi's and grad_phi's are assumed
 * to be calls for a weight function, grad of a weight function, a
 * basis function, or grad of a basis function.  If bfd->phi[], etc.,
 * is used AFTER this call is made, then you WILL NOT get the
 * regular steady-state solution.  Put another way: don't construct
 * steady-state variables outside of the above sequance...
 *
 * Once that sequence is executed for the regular steady-state
 * solution, we can modify the basis and weight functions to include
 * the 0/1 multiplications I figured out.  This allows us to create
 * (in 2 steps) the correct Jacobian/mass matrices in terms of the
 * cos(N z) and sin(N z) multipliers.
 *
 * The day I put this in, the global structure "bfd" is always used in
 * load_basis_functions().  Although it is passed to
 * load_basis_functions(), it is otherwise assumed to be global.  I
 * shall do the same.
 *
 * I'm not putting in all the mods for the "raw" derivatives
 * (dphidxi[i][a]), because those "shouldn't" be used in the
 * assemble_* routines.  If they are, then there's problems, because
 * those are coupled across a's for a resultant (i.e., there would be
 * mixing of sin's and cos's).  Bad bad bad!
 *
 * Author: Matt Hopkins, 12/18/00.  */
void
modify_basis_and_weight_functions_for_LSA_3D_of_2D(void)
{
  int i, j, a, b, bf_index, var, mn;
  BASIS_FUNCTIONS_STRUCT *bf_ptr;
  dbl N, cosine_value, sine_value, orig_phi;

  if(!af->Assemble_LSA_Jacobian_Matrix &&
     !af->Assemble_LSA_Mass_Matrix)
    return;

  mn = ei[pg->imtrx]->mn;
  N = LSA_3D_of_2D_wave_number;
  
  if(LSA_3D_of_2D_pass == 1)
    {
      cosine_value = 1.0;
      sine_value = 0.0;
    }
  else
    {
      cosine_value = 0.0;
      sine_value = 1.0;
    }

  /* The "special" case here is when the wavenumber is 0.  In that
   * case, we would like the original system to be computed when we
   * add our two passes.  The only piece that would get added that
   * shouldn't is the second pass sine_value changes.  These should
   * all be zero... Note that N == 0 kills a lot of the unnecessary
   * components, already.
   *
   * Unfortunately, this causes additional consequences because now
   * the w-momentum equations are all zero!  This causes zero rows in
   * the Jacobian... What to do about that? */
  if(N == 0.0)
    sine_value = 0.0;

  /*
  if(ei[pg->imtrx]->ielem == 0)
    printf("pass = %d, ielem = %d, cv = %3.1g, sv = %3.1g, N = %g\n",
	   LSA_3D_of_2D_pass, ei[pg->imtrx]->ielem, cosine_value, sine_value, N);
  */

  for(bf_index = 0; bf_index < Num_Basis_Functions; bf_index++)
    {
      bf_ptr = bfd[bf_index];
      var = bf_ptr->Var_Type_MatID[mn];
      if(var != -1)
	{
	  /* First, modify all of the basis functions and weight
	   * functions that are NOT associated with the
	   * w-velocity/w-momentum variable/equation.  The latter
	   * should have been tagged with the I_Q2_LSA or I_Q2_D_LSA
	   * weighting/interpolating functions.  */
	  if(bf_ptr->interpolation != I_Q2_LSA &&
	     bf_ptr->interpolation != I_Q2_D_LSA)
	    {
	      for(i = 0; i < ei[pg->imtrx]->dof[var]; i++)
		{
		  orig_phi = bf_ptr->phi[i];
		  bf_ptr->phi[i] *= cosine_value;
		  bf_ptr->grad_phi[i][0] *= cosine_value;
		  bf_ptr->grad_phi[i][1] *= cosine_value;
		  bf_ptr->grad_phi[i][2] = -N * sine_value * orig_phi;
		  /* Luckily, in CARTESIAN (and therefore,
		   * PROJECTED_CARTESIAN), grad_phi_e[i][a][p][a] =
		   * grad_phi[i][p] */
		  for(a = 0; a < VIM; a++)
		    for(b = 0; b < VIM; b++)
		      bf_ptr->grad_phi_e[i][a][b][a] =
			bf_ptr->grad_phi[i][b];
		  /* These next derivatives are ONLY used in the
		   * assemble_* routines for the mesh equations
		   * (assemble_mesh).  Note that grad_phi[i][2] == 0,
		   * so we don't need to multiply by anything. */
		  if(pd_glob[mn]->e[pg->imtrx][R_MESH1])
		    {
		      for(b = 0; b < 2; b++)
			for(j = 0; j < ei[pg->imtrx]->dof[R_MESH1]; j++)
			  {
			    bf_ptr->d_grad_phi_dmesh[i][0][b][j] *=
			      cosine_value;
			    bf_ptr->d_grad_phi_dmesh[i][1][b][j] *=
			      cosine_value;
			    /* Note that multiplying by N sine_value
                               is almost certainly wrong.  The "d_phi"
                               it's multiplying by should probably
                               just be "phi".  */
			    bf_ptr->d_grad_phi_dmesh[i][2][b][j] *=
			      -N * sine_value;
			    bf_ptr->d_d_phi_dmesh[i][0][b][j] *=
			      cosine_value;
			    bf_ptr->d_d_phi_dmesh[i][1][b][j] *=
			      cosine_value;
			    bf_ptr->d_d_phi_dmesh[i][2][b][j] *=
			      -N * sine_value;
			    for(a = 0; a < VIM; a++)
			      {
				bf_ptr->d_grad_phi_e_dmesh[i][a][0][a][b][j] *=
				  cosine_value;
				bf_ptr->d_grad_phi_e_dmesh[i][a][1][a][b][j] *=
				  cosine_value;
				bf_ptr->d_grad_phi_e_dmesh[i][a][2][a][b][j] *=
				  -N * sine_value;
			      }
			  }
		    }
		}
	    }
	  /* If it's a w-velocity/w-momentum variable/equation, then
	   * we want to modify their weighting/interpolating functions
	   * appropriately. */
	  else
	    {
	      for(i = 0; i < ei[pg->imtrx]->dof[var]; i++)
		{
		  orig_phi = bf_ptr->phi[i];
		  bf_ptr->phi[i] *= sine_value;
		  bf_ptr->grad_phi[i][0] *= sine_value;
		  bf_ptr->grad_phi[i][1] *= sine_value;
		  bf_ptr->grad_phi[i][2] = N * cosine_value * orig_phi;
		  for(a = 0; a < VIM; a++)
		    for(b = 0; b < VIM; b++)
		      bf_ptr->grad_phi_e[i][a][b][a] =
			bf_ptr->grad_phi[i][b];
		  /* These next derivatives are ONLY used in the
		   * assemble_* routines for the mesh equations
		   * (assemble_mesh).  Since this is
		   * w-velocity/equation stuff, I'm pretty sure that
		   * this is wasted work... Better to be consistent,
		   * though! */
		  if(pd_glob[mn]->e[pg->imtrx][R_MESH1])
		    for(b = 0; b < 2; b++)
		      for(j = 0; j < ei[pg->imtrx]->dof[R_MESH1]; j++)
			{
			  bf_ptr->d_grad_phi_dmesh[i][0][b][j] *=
			    sine_value;
			  bf_ptr->d_grad_phi_dmesh[i][1][b][j] *=
			    sine_value;
			  bf_ptr->d_grad_phi_dmesh[i][2][b][j] *=
			    N * cosine_value;
			  bf_ptr->d_d_phi_dmesh[i][0][b][j] *=
			    sine_value;
			  bf_ptr->d_d_phi_dmesh[i][1][b][j] *=
			    sine_value;
			  bf_ptr->d_d_phi_dmesh[i][2][b][j] *=
			    N * cosine_value;
			  for(a = 0; a < DIM; a++)
			    {
			      bf_ptr->d_grad_phi_e_dmesh[i][a][0][a][b][j] *=
				sine_value;
			      bf_ptr->d_grad_phi_e_dmesh[i][a][1][a][b][j] *=
				sine_value;
			      bf_ptr->d_grad_phi_e_dmesh[i][a][2][a][b][j] *=
				N * cosine_value;
			    }
			}
		}
	    }
	}
    }
}

void
modify_bf_mesh_derivs_for_LSA_3D_of_2D(void)

/************************************************************************
 *									*
 * modify_bf_mesh_derivs_for_LSA_3D_of_2D -- ancillary function for	*
 * load_bf_mesh_derivs.							*
 *									*
 * This function provides the entries of the arrays			*
 * d_grad_phi_dmesh and d_grad_phi_e_dmesh which correspond to		*
 * mesh derivatives of the first-order z-components of the basis	*
 * function gradients computed in load_fv_mesh_derivs.			*
 * These derivatives are nonzero and important in some 			*
 * free-surface 3D of 2D stability problems.				*
 * This function is called just after the basis & weight functions	*
 * are modified for 3D of 2D, and before equations are assembled.	*
 *									*
 * The general format used here is:					*
 *									*
 *     bf[v]->d_grad_phi_dmesh[i][2] [b][j] =				*
 *		- bf[v]->grad_phi[i][b] * bf[x]->grad_phi[j][2]		*
 *									*
 * where bf[v] points to each defined basis function (in turn), and	*
 * bf[x] points to the basis functions assigned to mesh displacement. 	*
 * Note that the MODIFIED basis functions are (intentionally) used here *
 *									*
 * Input and output are the same as for load_bf_mesh_derivs().		*
 *									*
 * Created:   Thu Sep 13 13:20:00 MDT 2001 edwilke@sandia.gov		*
 *									*
 ************************************************************************/
{
  int b, i, j, k, q;
  int mdof, vdof, var;
  struct Basis_Functions *bfv, *bfx;

  bfx = bf[R_MESH1];
  mdof = ei[pg->imtrx]->dof[R_MESH1];
  for (k=0; k<Num_Basis_Functions; k++)
    {
      bfv = bf[k];
      var = bfv->Var_Type_MatID[ei[pg->imtrx]->mn];
      if (var != -1)
	{
	  vdof = ei[pg->imtrx]->dof[var];
	  for (i=0; i<vdof; i++)
	    {
	      for (b=0; b<ei[pg->imtrx]->ielem_dim; b++)
		{
		  for (j=0; j<mdof; j++)
		    {
		      bfv->d_grad_phi_dmesh[i][2] [b][j] =
		             - bfv->grad_phi[i][b] * bfx->grad_phi[j][2];
		      for (q=0; q<VIM; q++)
			{
			  bfv->d_grad_phi_e_dmesh[i][q] [2][q] [b][j] =
			    bfv->d_grad_phi_dmesh[i][2] [b][j];
			}
		    }
		}
	    }
	}
    }
}

void
modify_fv_mesh_derivs_for_LSA_3D_of_2D(void)

/************************************************************************
 *									*
 * modify_fv_mesh_derivs_for_LSA_3D_of_2D -- ancillary function for	*
 * load_fv_mesh_derivs.							*
 *									*
 * This function provides the entries of the arrays			*
 * d_grad_(f)_dmesh which correspond to mesh derivatives of the		*
 * first-order z-components of the field variable gradients		*
 * computed in load_fv_mesh_derivs. These derivatives are nonzero	*
 * and important in some free-surface 3D of 2D stability problems.	*
 * This function is called just after the basis & weight functions	*
 * are modified for 3D of 2D, and before equations are assembled.	*
 *									*
 * The general format used here is:					*
 *									*
 *     fv-> d_grad_(f)_dmesh[2][b][j] =					*
 *		- fv->grad_(f)[b] * bf[x]->grad_phi[j][2]		*
 *									*
 * where (f) is any of the field variable gradients and bf[x] points	*
 * to the basis functions assigned to mesh displacement. 		*
 * Note that the MODIFIED basis functions are (intentionally) used here *
 *									*
 * Input and output are the same as for load_fv_mesh_derivs().		*
 *									*
 * Created:   Thu Sep 13 11:00:00 MDT 2001 edwilke@sandia.gov		*
 *									*
 ************************************************************************/
{
  int m, p, q, r, b, j, v, w, dim;
  int mdof;
  struct Basis_Functions *bfx;

/* Bail out fast if there's nothing to do */
  if ( !(pd->e[pg->imtrx][R_MESH1]) ) return;

/* Initialize values which will be constant */
  dim = ei[pg->imtrx]->ielem_dim;
  p = 2;
  mdof = ei[pg->imtrx]->dof[R_MESH1];
  bfx = bf[R_MESH1];

/* Modify each gradient in sequence */

/* d(grad(T))/dmesh */
  v = TEMPERATURE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_T_dmesh[p][b][j] =
		  - fv->grad_T[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(V))/dmesh */
  v = VOLTAGE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_V_dmesh[p][b][j] =
		  - fv->grad_V[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(qs))/dmesh */
  v = SURF_CHARGE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_qs_dmesh[p][b][j] =
		  - fv->grad_qs[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(sh_J))/dmesh */
  v = SHELL_DIFF_FLUX;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
        {
          for (j=0; j<mdof; j++)
            {
              fv->d_grad_sh_J_dmesh[p][b][j] =
                  - fv->grad_sh_J[b] * bfx->grad_phi[j][2];
            }
        }
    }

/* d(grad(apr,api))/dmesh */
  v = ACOUS_PREAL;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_apr_dmesh[p][b][j] =
		  - fv->grad_apr[b] * bfx->grad_phi[j][2];
	    }
	}
    }
  v = ACOUS_PIMAG;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_api_dmesh[p][b][j] =
		  - fv->grad_api[b] * bfx->grad_phi[j][2];
	    }
	}
    }
  v = LIGHT_INTP;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_poynt_dmesh[0][p][b][j] =
		  - fv->grad_poynt[0][b] * bfx->grad_phi[j][2];
	    }
	}
    }
  v = LIGHT_INTM;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_poynt_dmesh[1][p][b][j] =
		  - fv->grad_poynt[1][b] * bfx->grad_phi[j][2];
	    }
	}
    }
  v = LIGHT_INTD;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_poynt_dmesh[2][p][b][j] =
		  - fv->grad_poynt[2][b] * bfx->grad_phi[j][2];
	    }
	}
    }

  v = RESTIME;
  if (pd->v[v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_restime_dmesh[p][b][j] =
		  - fv->grad_restime[b] * bfx->grad_phi[j][2];
	    }
	}
    }  

  v = ACOUS_REYN_STRESS;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_ars_dmesh[p][b][j] =
		  - fv->grad_ars[b] * bfx->grad_phi[j][2];
	    }
	}
    }

  v = SHELL_BDYVELO;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_sh_bv_dmesh[p][b][j] =
		  - fv->grad_sh_bv[b] * bfx->grad_phi[j][2];
	    }
	}
    }

  v = SHELL_LUBP;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_sh_p_dmesh[p][b][j] =
                - fv->grad_sh_p[b] * bfx->grad_phi[j][2];
          }
      }
    }

  v = LUBP;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_lubp_dmesh[p][b][j] =
                - fv->grad_lubp[b] * bfx->grad_phi[j][2];
          }
      }
    }

 v = LUBP_2;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_lubp_2_dmesh[p][b][j] =
                - fv->grad_lubp_2[b] * bfx->grad_phi[j][2];
          }
      }
    }

  v = SHELL_TEMPERATURE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_sh_t_dmesh[p][b][j] =
                - fv->grad_sh_t[b] * bfx->grad_phi[j][2];
          }
      }
    }

  v = SHELL_FILMP;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_sh_fp_dmesh[p][b][j] =
                - fv->grad_sh_fp[b] * bfx->grad_phi[j][2];
          }
      }
    }

  v = SHELL_FILMH;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_sh_fh_dmesh[p][b][j] =
                - fv->grad_sh_fh[b] * bfx->grad_phi[j][2];
          }
      }
    }

  v = SHELL_PARTC;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
      {
        for (j=0; j<mdof; j++)
          {
            fv->d_grad_sh_pc_dmesh[p][b][j] =
                - fv->grad_sh_pc[b] * bfx->grad_phi[j][2];
          }
      }
    }


/* d(grad(SH))/dmesh */
  v = SHEAR_RATE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_SH_dmesh[p][b][j] =
		  - fv->grad_SH[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(F))/dmesh */
  v = FILL;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_F_dmesh[p][b][j] =
		  - fv->grad_F[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(P))/dmesh */
  v = PRESSURE;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_P_dmesh[p][b][j] =
		  - fv->grad_P[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(nn))/dmesh */
  v = BOND_EVOLUTION;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_nn_dmesh[p][b][j] =
		  - fv->grad_nn[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(c_w))/dmesh */
  v = MASS_FRACTION;
  for ( w=0; w<pd->Num_Species_Eqn; w++)
    {
      if (pd->v[pg->imtrx][v])
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_c_dmesh[p][w] [b][j] =
		      - fv->grad_c[w][b] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(porous_media_variables))/dmesh */

  v = POR_LIQ_PRES;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_p_liq_dmesh[p][b][j] =
		  - fv->grad_p_liq[b] * bfx->grad_phi[j][2];
	    }
	}
    }

  v = POR_GAS_PRES;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_p_gas_dmesh[p][b][j] =
		  - fv->grad_p_gas[b] * bfx->grad_phi[j][2];
	    }
	}
    }

  v = POR_POROSITY;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_porosity_dmesh[p][b][j] =
		  - fv->grad_porosity[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(v))/dmesh */
  v = VELOCITY1;
  if (pd->v[pg->imtrx][v])
    {
      for (q=0; q<VIM; q++)
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_v_dmesh[p][q][b][j] =
		      - fv->grad_v[b][q] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(ext_v))/dmesh */
  v = EXT_VELOCITY;
  if (pd->v[pg->imtrx][v])
    {
      for (b=0; b<dim; b++)
	{
	  for (j=0; j<mdof; j++)
	    {
	      fv->d_grad_ext_v_dmesh[p][b][j] =
		- fv->grad_ext_v[b] * bfx->grad_phi[j][2];
	    }
	}
    }

/* d(grad(E_field))/dmesh */
  v = EFIELD1;
  if (pd->v[pg->imtrx][v])
    {
      for (q=0; q<VIM; q++)
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_E_field_dmesh[p][q][b][j] =
		      - fv->grad_E_field[b][q] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(pv))/dmesh */
  v = PVELOCITY1;
  if (pd->v[pg->imtrx][v])
    {
      for (q=0; q<VIM; q++)
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_pv_dmesh[p][q][b][j] =
		      - fv->grad_pv[b][q] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(d))/dmesh */
  v = MESH_DISPLACEMENT1;
  if (pd->v[pg->imtrx][v])
    {
      for (q=0; q<VIM; q++)
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_d_dmesh[p][q][b][j] =
		      - fv->grad_d[b][q] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(d_rs))/dmesh */
  v = SOLID_DISPLACEMENT1;
  if (pd->v[pg->imtrx][v])
    {
      for (q=0; q<VIM; q++)
	{
	  for (b=0; b<dim; b++)
	    {
	      for (j=0; j<mdof; j++)
		{
		  fv->d_grad_d_rs_dmesh[p][q][b][j] =
		      - fv->grad_d_rs[b][q] * bfx->grad_phi[j][2];
		}
	    }
	}
    }

/* d(grad(S))/dmesh */
  v = POLYMER_STRESS11;
  if ( pd->v[pg->imtrx][v] )
    {
      for (m=0; m<vn->modes; m++)
        {
	  for (q=0; q<VIM; q++)
	    {
	      for (r=0; r<VIM; r++)
		{
		  for (b=0; b<dim; b++)
		    {
		      for (j=0; j<mdof; j++)
			{
			  fv->d_grad_S_dmesh[m][p][q][r] [b][j] =
			      - fv->grad_S[m][b][q][r] * bfx->grad_phi[j][2];
			}
		    }
		}
	    }
	}
    }

/* d(grad(G))/dmesh */
  v = VELOCITY_GRADIENT11;
  if ( pd->v[pg->imtrx][v] )
    {
      for (q=0; q<VIM; q++)
	{
	  for (r=0; r<VIM; r++)
	    {
	      for (b=0; b<dim; b++)
		{
		  for (j=0; j<mdof; j++)
		    {
		      fv->d_grad_G_dmesh[p][q][r] [b][j] =
			  - fv->grad_G[b][q][r] * bfx->grad_phi[j][2];
		    }
		}
	    }
	}
    }
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void
modify_normal_vector_for_LSA_3D_of_2D(void)

/*
 *  Function to insert the mesh sensitivities of the third (z)
 *  component of the normal vector for 3D of 2D stability of
 *  free surface flow problems. These terms are added to the
 *  fv->dsnormal_dx array.
 *  Author:  Edward Wilkes (9233) 9/13/2001.
 */

{
  int j, q;

/*
 *  Proceed only if doing 3D of 2D stability with deforming mesh.
 *  Also include modified mesh derivatives of stangent[1] in case
 *  they are needed.  NOTE: Binormal (stangent[0]) is not affected.
 */

  if (LSA_3D_of_2D_pass == 0 || !(pd->e[pg->imtrx][R_MESH1]) ) return;

  for (j=0; j<ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++)
    {
      for (q=0; q<2; q++)
        {
          fv->dsnormal_dx[2][q][j] = - fv->snormal[q]
                                     * bf[MESH_DISPLACEMENT1+q]->grad_phi[j][2];
          fv->dstangent_dx[1][0][q][j] = - fv->snormal[0] * fv->dsnormal_dx[2][q][j];
          fv->dstangent_dx[1][1][q][j] = - fv->snormal[1] * fv->dsnormal_dx[2][q][j];
        }
    }
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int
create_eigen_outfiles(Exo_DB *exo, Dpi *dpi,
                      RESULTS_DESCRIPTION_STRUCT *rd)
/*
 * Routine to open an array of MxN eigenvector output files.
 * M = specified Eigen_Record_Modes
 * N = specified LSA_number_wave_numbers
 * Pointers to these files are returned in array **fp.
 */
{
  int m = eigen->Eigen_Record_Modes;
  int n = LSA_number_wave_numbers;

  char fname[MAX_FNL];

  int i, j, k, n1;

/* Determine if doing 3D of 2D LSA */
  if (Linear_Stability == LSA_3D_OF_2D
   || Linear_Stability == LSA_3D_OF_2D_SAVE)
    {
      n1 = n;
    }
  else
    {
      n1 = 1;
      n = 0;
    }

/* Record mode loop */
  for (j=0; j<m; j++)
    {

/* Wave number loop */
      for (k=0; k<n1; k++)
        {

/* Create name for this file */
          for (i=0; i<MAX_FNL; i++) 
            {
              fname[i] = '\0';
            }
          strcpy(fname, eigen->Eigen_Output_File);
          get_eigen_outfile_name(fname, j, k);

/* Now, open a new ExodusII file with this composite name
   and initialize with basic mesh information */
          one_base(exo);
          wr_mesh_exo(exo, fname, 0);
          wr_result_prelim_exo(rd, exo, fname, NULL);
          if (Num_Proc > 1)
            wr_dpi(dpi, fname);
          zero_base(exo);

        }  /* End of wave number loop */

    }  /* End of mode loop */

/* Successful completion */
  return 1;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void
get_eigen_outfile_name(char *in_name, int j, int k)
{
  char prefix[MAX_FNL];
  char suffix[MAX_FNL];
  char new_name[MAX_FNL];
  char mode_string[MAX_FNL];
  char wave_string[MAX_FNL];

  int i;
  int suffix_length;
  int prefix_length;

/* Initialize string pieces */
  for ( i=0; i<MAX_FNL; i++)
    {
      prefix[i] = '\0';
      suffix[i] = '\0';
      new_name[i] = '\0';
      mode_string[i] = '\0';
      wave_string[i] = '\0';
    }

/* Save beginning and end of specified file name */
  suffix_length = get_suffix(suffix, in_name);
  prefix_length = get_prefix(prefix, in_name);

/* Determine if doing 3D of 2D LSA; set k = -1 if not */
  if ( !(Linear_Stability == LSA_3D_OF_2D)
    && !(Linear_Stability == LSA_3D_OF_2D_SAVE)) k = -1;

/* Always construct mode string */
  sprintf(mode_string, "_mode%d", j);

/* Construct wave string for 3D of 2D LSA only */
  if (k != -1) sprintf(wave_string, "_wn=%g", LSA_wave_numbers[k]);

/* Assemble the pieces of the name string */
  if (prefix_length > 0) strcpy(new_name, prefix);
  strcat(new_name, mode_string);
  if (k != -1) strcat(new_name, wave_string);
  if (suffix_length > 0) strcat(new_name, suffix);

/* Multiname for this processor if parallel */
  if (Num_Proc > 1) multiname(new_name, ProcID, Num_Proc);

/* Return the name as in_name */
  strcpy(in_name, new_name);
  return;
}
