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

/* mm_fill_stress -- auxiliary routines helpful in matrix & residual assembly 
 * for polymer stress equations
 *
 *   Copyright (C) 1993, 1994  Sandia National Laboratories
 */

/* Standard include files */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>

/* GOMA include files */
#include "std.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_eh.h"
#include "mm_fill_stress.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"

#define _MM_FILL_STRESS_C
#include "goma.h"

extern struct Boundary_Condition *inlet_BC[MAX_VARIABLE_TYPES+MAX_CONC];

// direct call to a fortran LAPACK eigenvalue routine
extern FSUB_TYPE dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA,
			double *W, double *WORK, int *LWORK, int *INFO,
			int len_jobz, int len_uplo);

/*  _______________________________________________________________________  */

/* assemble_stress -- assemble terms (Residual &| Jacobian) for polymer stress eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 * 	ija -- vector of pointers into the a matrix
 * 	a  -- global Jacobian matrix
 * 	R  -- global residual vector
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Dec  8 14:03:06 MST 1993 pasacki@sandia.gov
 *
 * Revised:	Sun Feb 27 06:53:12 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int 
assemble_stress(dbl tt,		/* parameter to vary time integration from 
				   explicit (tt = 1) to implicit (tt = 0) */
		dbl dt,		/* current time step size */
		dbl h[DIM],
		dbl hh[DIM][DIM],
		dbl dh_dxnode[DIM][MDE],
		dbl vcent[DIM],	/* average element velocity, which is the
				 * centroid velocity for Q2 and the average of
				 * the vertices for Q1. It comes from the 
				 * routine "element_velocity."               */
		dbl dvc_dnode[DIM][MDE])
{
  int dim, p, q, r, a, b, w = -1;
  
  int mode; /*counter for viscoelastic modes */
  
  int eqn, var;
  int peqn, pvar;

  int i, j, status;

  dbl v[DIM];			        /* Velocity field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl h3;		        	/* Volume element (scale factors). */
  dbl dh3dmesh_pj;	        	/* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];                  /* Shear-rate tensor based on velocity */
  dbl det_J;                            /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj;			/* for specific (p,j) mesh dof */

  dbl mass;			        /* For terms and their derivatives */
  dbl mass_a, mass_b, mass_c;	
  dbl advection;	
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a, source_b, source_c;

  /*
   * 
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */


  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals 
   * and some of their derivatives...
   */

  dbl wt_func;


  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM]; 
  int v_s[MAX_MODES][DIM][DIM]; 
  int v_g[DIM][DIM]; 

  dbl s[DIM][DIM];         /* stress tensor */
  dbl s_dot[DIM][DIM];     /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];         /* velocity gradient tensor */
  dbl gt[DIM][DIM];        /* transpose of velocity gradient tensor */
  dbl g_dot[DIM][DIM];     /* velocity gradient tensor time derivative */
  dbl gt_dot[DIM][DIM];    /* transpose of velocity gradient tensor time derivative */

  /* dot product tensors */

  dbl g_dot_g[DIM][DIM];
  dbl g_dot_gt[DIM][DIM];
  dbl gt_dot_g[DIM][DIM]; 
  dbl gt_dot_gt[DIM][DIM];
  dbl s_dot_s[DIM][DIM]; 
  dbl s_dot_g[DIM][DIM]; 
  dbl s_dot_gt[DIM][DIM]; 
  dbl g_dot_s[DIM][DIM];
  dbl gt_dot_s[DIM][DIM]; 

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  dbl d_mup_dv_pj;
  dbl d_mup_dmesh_pj;

  /* constitutive equation parameters */
  dbl alpha;     /* This is the Geisekus mobility parameter */
  dbl lambda=0;    /* polymer relaxation constant */

  /* advective terms are precalculated */
  dbl v_dot_del_g[DIM][DIM];
  dbl x_dot_del_g[DIM][DIM];
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdelg_dm;
  dbl d_xdotdelgt_dm;
  dbl d_xdotdels_dm;

  dbl d_vdotdelg_dm;
  dbl d_vdotdelgt_dm;
  dbl d_vdotdels_dm;

/* SUPG variables */
  dbl h_elem=0, h_elem_inv=0, h_elem_deriv=0;
  dbl supg=0;

  status = 0;

  eqn   = R_STRESS11;			

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			/* Differential volume element (scales). */


  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32; 
  v_g[2][2] = VELOCITY_GRADIENT33; 


  /*
   * Field variables...
   */
  
  for ( a=0; a<dim; a++)
    {
      v[a] = fv->v[a];

      /* note, these are zero for steady calculations */
      if (  pd->TimeIntegration != STEADY &&  pd->v[MESH_DISPLACEMENT1+a] )
	{
	  x_dot[a] = fv_dot->x[a];
	}
      else
	{
	  x_dot[a] = 0.;
	}
    }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */
  
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  grad_v[a][b] = fv->grad_v[a][b];
	}
    }
  
  /* load up shearrate tensor based on velocity */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  g[a][b]  = fv->G[a][b];
	  gt[b][a] = g[a][b];
	  if (pd->TimeIntegration != STEADY) {
	    g_dot[a][b] = fv_dot->G[a][b];
	    gt_dot[b][a] = g_dot[a][b];
	  } else {
	    g_dot[a][b] = 0.0;
	    gt_dot[a][b] =  0.0;
	  }
	}
    }
  
  if( vn->wt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( vn->wt_funcModel == SUPG)
    {
      supg = vn->wt_func;
    }


  if(supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p]*vcent[p]*h[p];
	}
      h_elem = sqrt(h_elem)/2.;
      if(h_elem == 0.) 
	{
	  h_elem_inv=1.;
	}
      else
	{
	  h_elem_inv=1./h_elem;
	}
	
    }
/* end Petrov-Galerkin addition */

  /* get tensor dot products for future use */
  
  (void) tensor_dot(g,  g,  g_dot_g,   VIM);
  (void) tensor_dot(g,  gt, g_dot_gt,  VIM);
  (void) tensor_dot(gt, g,  gt_dot_g,  VIM);
  (void) tensor_dot(gt, gt, gt_dot_gt, VIM);

  (void) stress_eqn_pointer(v_s);
  (void) stress_eqn_pointer(R_s);

  /* Begin loop over modes */
  for ( mode=0; mode<vn->modes; mode++)
    {
      
      load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);
      
      /* precalculate advective terms of form (v dot del tensor)*/
      
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      v_dot_del_g[a][b] = 0.;
	      x_dot_del_g[a][b] = 0.;
	      v_dot_del_s[a][b] = 0.;
	      x_dot_del_s[a][b] = 0.; 
	      for ( q=0; q<dim; q++)
		{
		  v_dot_del_g[a][b] +=  v[q] * fv->grad_G[q][a][b];
		  x_dot_del_g[a][b] +=  x_dot[q] * fv->grad_G[q][a][b];
		  v_dot_del_s[a][b] +=  v[q] * grad_s[q][a][b];
		  x_dot_del_s[a][b] +=  x_dot[q] * grad_s[q][a][b];
		} 
	    }
	}
      
      
      
      /*
       * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
       */
      
      /* get polymer viscosity */
      mup = viscosity(ve[mode]->gn, gamma, d_mup);

      /* get Geisekus mobility parameter */
      alpha = ve[mode]->alpha;

      /* get time constant */
      if(ve[mode]->time_constModel == CONSTANT)
	{
	  lambda = ve[mode]->time_const;
	}
      else if(ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW)
	{
	  lambda = mup/ve[mode]->time_const;
	}
      
      
      /* get tensor dot products for future use */
      
      (void) tensor_dot(s, s, s_dot_s, VIM);
      (void) tensor_dot(s, g, s_dot_g, VIM);
      (void) tensor_dot(s, gt, s_dot_gt, VIM);
      (void) tensor_dot(gt, s, gt_dot_s, VIM);
      (void) tensor_dot(g, s, g_dot_s, VIM);
      
      /*
       * Residuals_________________________________________________________________
       */
      
      if ( af->Assemble_Residual )
	{
	  /*
	   * Assemble each component "ab" of the polymer stress equation...
	   */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      eqn = R_s[mode][a][b];
		      
		      peqn = upd->ep[eqn];
		      
		      /*
		       * In the element, there will be contributions to this many equations
		       * based on the number of degrees of freedom...
		       */
		      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  mass = 0.;
			  
			  if ( pd->TimeIntegration != STEADY )
			    {
			      if ( pd->e[eqn] & T_MASS )
				{
				  mass = s_dot[a][b] + mup *(g_dot[a][b] + gt_dot[a][b]);
				  mass *= wt_func * lambda * det_J * wt;
				  mass *= h3;
				  mass *= pd->etm[eqn][(LOG2_MASS)];
				}
			    }
			  
			  advection = 0.;
			  if ( pd->e[eqn] & T_ADVECTION )
			    {
			      if(DOUBLE_NONZERO(lambda))
				{
				  
				  advection += ( v_dot_del_s[a][b]  -  x_dot_del_s[a][b]
						 + mup*(v_dot_del_g[a][b] - x_dot_del_g[a][b] 
							+ v_dot_del_g[b][a] - x_dot_del_g[b][a] ));
				  advection -= ( mup*(2.*gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b] ) 
						 + gt_dot_s[a][b] + s_dot_g[a][b] );
				  
				  advection *= wt_func * lambda *det_J * wt * h3;
				  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
				}     
			    }
			  
			  diffusion = 0.;
			  if ( pd->e[eqn] & T_DIFFUSION )
			    {
			      /* add SU term in here when appropriate */
			      diffusion += 0.;
			      diffusion *= -wt_func * det_J * wt * h3;
			      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			    }
			  
			  /*
			   * Source term...
			   */
			  
			  source = 0.;
			  source1 = 0.;
			  
			  if ( pd->e[eqn] & T_SOURCE )
			    {
			      source +=  s[a][b];
			      if(DOUBLE_NONZERO(alpha))
				{
				  source1 += (s_dot_g[a][b] + s_dot_gt[a][b] 
					      + g_dot_s[a][b] + gt_dot_s[a][b] + s_dot_s[a][b]/mup 
					      + mup * ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
							+ gt_dot_g[a][b] + g_dot_gt[a][b]));
				  
				  source1 *= alpha * lambda;
				  source  += source1;
				}
			      
			      source *= wt_func * det_J * h3 * wt;
			      
			      source *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
			  
			  /*
			   * Add contributions to this residual (globally into Resid, and 
			   * locally into an accumulator)
			   */
			  
			  lec->R[peqn][i] += 
			    mass + advection + diffusion + source;
			}
		    }
		}
	    }
	}
      
      /*
       * Jacobian terms...
       */
      
      if ( af->Assemble_Jacobian )
	{
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      
		      eqn = R_s[mode][a][b];
		      peqn = upd->ep[eqn];
		      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  
			  /*
			   * Set up some preliminaries that are needed for the (a,i)
			   * equation for bunches of (b,j) column variables...
			   */
			  
			  /* 
			   * J_S_T
			   */
			  
			  var = TEMPERATURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  mass = 0.;
				  
				  if ( pd->TimeIntegration != STEADY )
				    {
				      if ( pd->e[eqn] & T_MASS )
					{
					  mass = d_mup->T[j] *(g_dot[a][b] + gt_dot[a][b]);
					  mass *= wt_func;
					  mass *= h3 * lambda * det_J * wt;
					  mass *= pd->etm[eqn][(LOG2_MASS)];
					}
				    }
				  
				  advection = 0.;
				  if ( pd->e[eqn] & T_ADVECTION )
				    {
				      if(DOUBLE_NONZERO(lambda))
					{
					  advection += ( v_dot_del_g[a][b] - x_dot_del_g[a][b] +
							 v_dot_del_g[b][a] - x_dot_del_g[b][a]) -
					    2.*gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
					  advection *= wt_func * d_mup->T[j]* lambda * det_J * wt * h3;
					  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					}
				    }
				  
				  
				  diffusion = 0.;
				  if ( pd->e[eqn] & T_DIFFUSION )
				    {
				      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
				    }
				  
				  source    = 0.;
				  source1    = 0.;
				  if ( pd->e[eqn] & T_SOURCE )
				    {
			              if(DOUBLE_NONZERO(alpha))
					{
					  source1 += ( -s_dot_s[a][b]/(mup*mup) +
						       ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
							 + gt_dot_g[a][b] + g_dot_gt[a][b]));
					  source1 *= lambda * alpha * d_mup->T[j];
					  source  += source1;
					}
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    mass + advection + diffusion + source;
				}
			    }
			  
			  /*
			   * J_S_v
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = VELOCITY1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_mup_dv_pj  = d_mup->v[p][j];
				      
				      mass = 0.;
				      
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      mass = d_mup_dv_pj *(g_dot[a][b] + gt_dot[a][b]) *wt_func;
					      
					      if(supg!=0.)
						{
						  mass_a = supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
						  
						  for(w=0;w<dim;w++)
						    {
						      mass_a += supg * vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
							*v[w]*bf[eqn]->grad_phi[i] [w];
						    }
						  
						  mass_a*=  s_dot[a][b] + mup *(g_dot[a][b] + gt_dot[a][b]);
						  mass += mass_a;
						}
					      
					      mass *= pd->etm[eqn][(LOG2_MASS)] * lambda * det_J * wt * h3;
					    }
					  
					}
				      
				      advection = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  if(DOUBLE_NONZERO(lambda))
					    {
					      advection_a = phi_j * 
						(grad_s[p][a][b] + 
						 mup*(fv->grad_G[p][a][b] + fv->grad_G[p][b][a]));
					      
					      advection_a *=  wt_func;
					      advection_b = ( v_dot_del_g[a][b]  -  x_dot_del_g[a][b] + 
							      v_dot_del_g[b][a]  -  x_dot_del_g[b][a]);
					      
					      advection_b -= 2.*gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b];
					      
					      advection_b *=   wt_func * d_mup_dv_pj; 
					      
					      advection_c = 0.;
					      /* Petrov-Galerkin term */
					      if(supg !=0.)
						{
						  advection_c =  supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
						  for( w =0; w<dim; w++ )
						    {
						      advection_c += supg*vcent[p]*h[p]*dvc_dnode[p][j]*h_elem_inv/4.
							*v[w]*bf[eqn]->grad_phi[i] [w];
						    }
						  
						  advection_c*= (v_dot_del_s[a][b]  -  x_dot_del_s[a][b]
								 + mup*(v_dot_del_g[a][b] - x_dot_del_g[a][b] 
									+ v_dot_del_g[b][a] - x_dot_del_g[b][a] )
								 - mup*(2.*gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b] ) 
								 - gt_dot_s[a][b] - s_dot_g[a][b]);
						}
					      
					      advection = advection_a + advection_b + advection_c;
					      advection *= lambda * det_J * wt *h3;
					      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					    }
					}
				      
				      diffusion = 0.;
				      
				      if ( pd->e[eqn] & T_DIFFUSION )
					{
					  /* add SU term in here when appropriate */
					  
					  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					}
				      
				      source    = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a = 0.;
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source_a = -s_dot_s[a][b]/(mup*mup) +  g_dot_g[a][b] 
						+  gt_dot_gt[a][b] 
						+ gt_dot_g[a][b] + g_dot_gt[a][b];
					      source_a *= wt_func * alpha * lambda * d_mup_dv_pj;
					    }
					  
					  source_b = 0.;
					  if(supg != 0.)
					    {
					      source_b = supg * h_elem* phi_j*bf[eqn]->grad_phi[i][p];
					      
					      for(w=0;w<dim;w++)
						{
						  source_b += supg*vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
						    *v[w]*bf[eqn]->grad_phi[i] [w];
						}
					      
					      source_b *= s[a][b] + (s_dot_g[a][b] + s_dot_gt[a][b] 
								     + g_dot_s[a][b] + gt_dot_s[a][b] + s_dot_s[a][b]/mup 
								     + mup * ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
									       + gt_dot_g[a][b] + g_dot_gt[a][b])) *
						alpha * lambda;
					    }
					  
					  source = source_a + source_b;
					  source *=  det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					}
				      
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + diffusion + source;
				    }
				}
			    }
			  
			  /*
			   * J_S_c
			   */
			  var = MASS_FRACTION;
			  if ( pd->v[var] )
			    {
			      pvar = MAX_PROB_VAR + w;
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  for ( w=0; w<pd->Num_Species_Eqn; w++)
				    {
				      
				      mass    = 0.;
				      
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      mass = d_mup->C[w][j] *(g_dot[a][b] + gt_dot[a][b]);
					      mass *= wt_func * lambda * det_J * wt;
					      mass *= h3;
					      mass *= pd->etm[eqn][(LOG2_MASS)];
					    }
					}
				      
				      
				      advection = 0.;
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  if(DOUBLE_NONZERO(lambda))
					    {
					      advection += ( v_dot_del_g[a][b] - x_dot_del_g[a][b] +
							     v_dot_del_g[b][a] - x_dot_del_g[b][a]) -
						2.*gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
					      advection *= d_mup->C[w][j] * wt_func * lambda * det_J * wt * h3;
					      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					    }
					}
				      
				      diffusion = 0.;
				      if ( pd->e[eqn] & T_DIFFUSION )
					{
					  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					}
				      
				      source    = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source += -s_dot_s[a][b]/(mup*mup) + g_dot_g[a][b] 
						+  gt_dot_gt[a][b] 
						+ gt_dot_g[a][b] + g_dot_gt[a][b];
					      source *= alpha * lambda * d_mup->C[w][j];
					    }
					  
					  source *= wt_func * det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					  
					}
				      
				      if ( w > 1 )
					{
					  EH(-1, "Need more arrays for each species.");
					}
				      
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + diffusion + source;
				    }
				}
			    }
			  
			  /*
			   * J_S_P
			   */
			  var = PRESSURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  mass    = 0.;
				  if ( pd->TimeIntegration != STEADY )
				    {
				      if ( pd->e[eqn] & T_MASS )
					{
					  mass =  d_mup->P[j] *(g_dot[a][b] + gt_dot[a][b]);
					  mass *= wt_func * lambda * det_J * wt;
					  mass *= h3;
					  mass *= pd->etm[eqn][(LOG2_MASS)];
					}
				    }
				  
				  advection = 0.;
				  if ( pd->e[eqn] & T_ADVECTION )
				    {
				      if(DOUBLE_NONZERO(lambda))
					{
					  
					  advection += ( v_dot_del_g[a][b] - x_dot_del_g[a][b] +
							 v_dot_del_g[b][a] - x_dot_del_g[b][a]) -
					    2.*gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
					  advection *= wt_func  * lambda * det_J * wt * h3 * d_mup->P[j];
					  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					}
				    }
				  
				  diffusion = 0.;
				  
				  if ( pd->e[eqn] & T_DIFFUSION )
				    {
				      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
				    }
				  
				  source    = 0.;
				  
				  if ( pd->e[eqn] & T_SOURCE )
				    {
			              if(DOUBLE_NONZERO(alpha))
					{
					  source += ( -s_dot_s[a][b]/(mup*mup) 
						      + ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
							  + gt_dot_g[a][b] + g_dot_gt[a][b]));
					  source *= d_mup->P[j] * alpha * lambda;
					}
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    mass + advection + diffusion + source;
				}		      
			    }
			  
			  /*
			   * J_S_d
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
				      dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
				      d_mup_dmesh_pj  = d_mup->X [p][j];
				      if(supg!=0.)
					{
					  h_elem_deriv = 0.;
					  for( q=0; q<dim; q++ )
					    {
					      h_elem_deriv += 
						hh[q][p]*vcent[q]*vcent[q]*dh_dxnode[q][j]*h_elem_inv/4.;
					    } 
					}
				      
				      mass = 0.;
				      mass_a = 0.;
				      mass_b = 0.;
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      mass_a = s_dot[a][b] + mup *(g_dot[a][b] + gt_dot[a][b]);
					      mass_a *= wt_func * ( d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					      
					      mass_b = d_mup_dmesh_pj *(g_dot[a][b] + gt_dot[a][b]);
					      mass_b *= wt_func * h3 * det_J;
					      
					      mass_c =0.;
					      if(supg != 0.)
						{
						  for( w=0; w<dim; w++ )
						    {
						      mass_c += supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
									+  h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						    }
						  mass_c*= (s_dot[a][b] + mup *(g_dot[a][b] + gt_dot[a][b])) * h3 * det_J;
						}
					      
					      mass = mass_a + mass_b + mass_c;
					      mass *= lambda * wt;
					      mass *= pd->etm[eqn][(LOG2_MASS)];
					    }
					}
				      
				      advection   = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  if(DOUBLE_NONZERO(lambda))
					    {
					      /*
					       * Four parts: 
					       *    advection_a = 
					       *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					       *
					       *    advection_b = 
					       *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
					       *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					       *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
					       *
					       * For unsteady problems, we have an 
					       * additional term
					       *
					       *    advection_c = 
					       *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
					       */
					      
					      advection_a = 0.;
					      
					      advection_a += ( v_dot_del_s[a][b]  -  x_dot_del_s[a][b]
							       + mup*(v_dot_del_g[a][b] - x_dot_del_g[a][b] 
								      + v_dot_del_g[b][a] - x_dot_del_g[b][a]) );
					      advection_a -= ( mup*(2.*gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b] ) 
							       + gt_dot_s[a][b] + s_dot_g[a][b] );
					      
					      
					      advection_a *= wt_func *(  d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					      
					      advection_b = 0.;		
					      d_vdotdelg_dm = 0.;
					      d_vdotdelgt_dm = 0.;
					      d_vdotdels_dm = 0.;
					      for ( q=0; q<dim; q++)
						{
						  d_vdotdels_dm += (v[q]-x_dot[q]) * d_grad_s_dmesh[q][a][b] [p][j];
						  d_vdotdelg_dm += (v[q]-x_dot[q]) * fv->d_grad_G_dmesh[q][a][b] [p][j];
						  d_vdotdelgt_dm += (v[q]-x_dot[q]) * fv->d_grad_G_dmesh[q][b][a] [p][j];
						}
					      
					      advection_b +=  (d_vdotdels_dm + mup*(d_vdotdelg_dm + d_vdotdelgt_dm));
					      advection_b += d_mup_dmesh_pj *
						(( v_dot_del_g[a][b] - x_dot_del_g[a][b] +
						   v_dot_del_g[b][a] - x_dot_del_g[b][a]) -
						 2.*gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b]);
					      advection_b *=  wt_func *det_J * h3;
					      
					      advection_c = 0.;	
					      if ( pd->TimeIntegration != STEADY )
						{
						  if ( pd->e[eqn] & T_MASS )
						    {
						      d_xdotdels_dm = (1.+2.*tt) * phi_j/dt 
							* grad_s[p][a][b];
						      d_xdotdelg_dm = (1.+2.*tt) * phi_j/dt 
							* fv->grad_G[p][a][b];
						      d_xdotdelgt_dm = (1.+2.*tt) * phi_j/dt 
							* fv->grad_G[p][b][a];
						      
						      advection_c -= (d_xdotdels_dm + mup*(d_xdotdelg_dm + d_xdotdelgt_dm)) ;
						      
						      advection_c *= wt_func * h3 * det_J;
						    }
						}
					      
					      advection_d = 0.;	
					      if(supg!=0.)
						{
						  for( w=0; w<dim; w++ )
						    {
						      advection_d+= supg *(h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
									   + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						    }
						  
						  advection_d *= (v_dot_del_s[a][b]  -  x_dot_del_s[a][b]
								  + mup*(v_dot_del_g[a][b] - x_dot_del_g[a][b] 
									 + v_dot_del_g[b][a] - x_dot_del_g[b][a] )
								  - mup*(2.*gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b] ) 
								  - gt_dot_s[a][b] - s_dot_g[a][b]) * det_J * h3;
						  
						  
						}
					      
					      advection = advection_a + advection_b + advection_c + advection_d;
					      
					      advection *=  wt * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
					      
					    }
					}
				      
				      /*
				       * Diffusion...
				       */
				      
				      diffusion = 0.;
				      
				      if ( pd->e[eqn] & T_DIFFUSION )
					{
					  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					}
				      
				      /*
				       * Source term...
				       */
				      
				      
				      source = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a =  s[a][b];
					  source_b = 0.;
					  
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source_a += (s_dot_g[a][b] + s_dot_gt[a][b] 
							   + g_dot_s[a][b] + gt_dot_s[a][b] + s_dot_s[a][b]/mup 
							   + mup * ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
								     + gt_dot_g[a][b] + g_dot_gt[a][b])) * alpha * lambda;
					      source_b += -s_dot_s[a][b]/(mup*mup) + ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
										       + gt_dot_g[a][b] + g_dot_gt[a][b]);
					      source_b *= d_mup_dmesh_pj * alpha * lambda;
					    }
					  
					  source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
					  
					  source_b *= wt_func *det_J * h3;
					  
					  source_c = 0.;
					  if(supg !=0.)
					    {
					      
					      for( w=0;w<dim;w++)
						{
						  source_c+= supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
								     + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] ); 
						}
					      source_c *= (s[a][b]+
							   alpha * lambda * 
							   (s_dot_g[a][b] + s_dot_gt[a][b] 
							    + g_dot_s[a][b] + gt_dot_s[a][b] + s_dot_s[a][b]/mup 
							    + mup * ( g_dot_g[a][b] +  gt_dot_gt[a][b] 
								      + gt_dot_g[a][b] + g_dot_gt[a][b]))) * det_J * h3;
					    }
					  
					  source  = source_a + source_b + source_c;
					  
					  source *=  wt * pd->etm[eqn][(LOG2_SOURCE)];
					  
					}
				      
				      lec->J[peqn][pvar][i][j] +=
					mass + advection  + diffusion + source;
				    }
				}
			    }
			  
			  /*
			   * J_S_G
			   */
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  var =  v_g[p][q];
				  
				  if ( pd->v[var] )
				    {
				      pvar = upd->vp[var];
				      for ( j=0; j<ei->dof[var]; j++)
					{
					  phi_j = bf[var]->phi[j];
					  mass = 0.;
					  if ( pd->TimeIntegration != STEADY )
					    {
					      if ( pd->e[eqn] & T_MASS )
						{
						  mass = mup *(1.+2.*tt) * phi_j/dt * ((double)delta(a,p) * (double)delta(b,q) 
										       + (double)delta(b,p) * (double)delta(a,q)); 
						  mass *= h3 * det_J;
						  
						  mass *= wt_func * lambda * wt * pd->etm[eqn][(LOG2_MASS)];
						}
					    }
					  
					  advection   = 0.;
					  advection_a   = 0.;
					  
					  if ( pd->e[eqn] & T_ADVECTION )
					    {
					      if(DOUBLE_NONZERO(lambda))
						{
						  
						  for( r=0; r<dim; r++)
						    {
						      advection_a +=  mup * (v[r]-x_dot[r])* bf[var]->grad_phi[j][r];
						    }
						  
						  advection += advection_a * 
						    ((double)delta(a,p) * (double)delta(b,q) + (double)delta(b,p) * (double)delta(a,q));
						  
						  advection -= mup* phi_j *
						    (2. * (g[p][b] * (double)delta(a,q) + gt[a][p] * (double)delta(b,q) )
						     +  (gt[p][b] * (double)delta(a,q) + gt[a][q] * (double)delta(b,p) )
						     +  ( g[q][b] * (double)delta(a,p) +  g[a][p] * (double)delta(b,q) ));
						  
						  advection -=  phi_j * (s[p][b] * (double)delta(a,q) + s[a][p] * (double)delta(b,q));
						  
						  advection *=  h3 * det_J ;
						  
						  advection *= wt_func * wt * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
						} 
					    }
					  
					  /*
					   * Diffusion...
					   */
					  
					  diffusion = 0.;
					  
					  if ( pd->e[eqn] & T_DIFFUSION )
					    {
					      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					    }
					  
					  /*
					   * Source term...
					   */
					  
					  source = 0.;		      
					  
					  if ( pd->e[eqn] & T_SOURCE )
					    {
					      
			                      if(DOUBLE_NONZERO(alpha))
						{
						  source +=  phi_j *( s[a][p] * (double)delta(b,q) 
								      + s[a][q] * (double)delta(b,p) 
								      + s[q][b] * (double)delta(a,p)
								      + s[p][b] * (double)delta(a,q))  
						    + mup * phi_j *( g[q][b] * (double)delta(a,p) +  g[a][p] * (double)delta(b,q) 
								     + gt[p][b] * (double)delta(a,q) + gt[a][q] * (double)delta(b,p) 
								     +  g[p][b] * (double)delta(a,q) + gt[a][p] * (double)delta(b,q) 
								     +  g[a][q] * (double)delta(b,p) + gt[q][b] * (double)delta(a,p));
						  source *= alpha * lambda;
						}
					      
					      source *= det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					      
					    }
					  
					  lec->J[peqn][pvar][i][j] +=
					    mass + advection  + diffusion + source;
					}
				    }
				}
			    }
			  
			  
			  /*
			   * J_S_S
			   */
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  var =  v_s[mode][p][q];
				  
				  if ( pd->v[var] )
				    {
				      pvar = upd->vp[var];
				      for ( j=0; j<ei->dof[var]; j++)
					{
					  phi_j = bf[var]->phi[j];
					  mass = 0.;
					  if ( pd->TimeIntegration != STEADY )
					    {
					      if ( pd->e[eqn] & T_MASS )
						{
						  mass = (1.+2.*tt) * phi_j/dt * (double)delta(a,p) * (double)delta(b,q); 
						  mass *= h3 * det_J;
						  mass *= wt_func * lambda * wt * pd->etm[eqn][(LOG2_MASS)];
						}
					    }
					  
					  advection   = 0.;
					  
					  if ( pd->e[eqn] & T_ADVECTION )
					    {
					      if(DOUBLE_NONZERO(lambda))
						{
						  if((a == p) && (b == q))
						    {
						      for( r=0; r<dim; r++)
							{
							  advection +=  (v[r]-x_dot[r])*  bf[var]->grad_phi[j][r];
							}
						    }
						  
						  advection -=  phi_j * (gt[a][p] * (double)delta(b,q) + g[q][b] * (double)delta(a,p));
						  
						  advection *=  h3 * det_J ;
						  
						  advection *= wt_func * wt * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
						}
					    }
					  
					  /*
					   * Diffusion...
					   */
					  
					  diffusion = 0.;
					  
					  if ( pd->e[eqn] & T_DIFFUSION )
					    {
					      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					    }
					  
					  /*
					   * Source term...
					   */
					  
					  source = 0.;		      
					  
					  if ( pd->e[eqn] & T_SOURCE )
					    {
					      source_a  =  phi_j * (double)delta(a,p) * (double)delta(b,q);
					      
					      source_b  = 0.;
			                      if(DOUBLE_NONZERO(alpha))
						{
						  source_b  =  phi_j *( g[q][b] * (double)delta(a,p) 
									+ gt[q][b] * (double)delta(a,p) 
									+ g[a][p] * (double)delta(b,q)
									+ gt[a][p] * (double)delta(b,q))
						    + phi_j * (s[q][b] * (double)delta(a,p) + s[a][p] * (double)delta(b,q))/mup;
						  
						  source_b *= alpha * lambda;
						}
					      
					      source  = source_a + source_b;
					      
					      source *= det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					      
					    }
					  
					  lec->J[peqn][pvar][i][j] +=
					    mass + advection + diffusion + source;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

  return(status);
}

/* this stress routine does the EVSS formulation according to Fortin, 1995
 *who uses the regular stress equation and converts stress in the momentum
 *equation by adding the divergence of (g + gT).
 */

int 
assemble_stress_fortin(dbl tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
		       dbl dt,	/* current time step size */
		       dbl h[DIM], /* coordinate scale factors */
		       dbl hh[DIM][DIM], /* coordinate scale factors */
		       dbl dh_dxnode[DIM][MDE],
		       dbl vcent[DIM], /* Average element velocity, which is 
					* the centroid velocity for Q2 and the
					* average of the vertices for Q1. It 
					* comes from the routine 
					* "element_velocity." */
		       dbl dvc_dnode[DIM][MDE])
{
  int dim, p, q, r, a, b, w, k;

  int eqn, var;
  int peqn, pvar;
  int evss_gradv=0;

  int i, j, status, mode;
  dbl v[DIM];			        /* Velocity field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl h3;		        	/* Volume element (scale factors). */
  dbl dh3dmesh_pj;	        	/* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];                  /* Shear-rate tensor based on velocity */
  dbl det_J;                            /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj;			/* for specific (p,j) mesh dof */

  dbl mass;			        /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;	
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a=0, source_b=0, source_c=0;

  /*
   * 
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals 
   * and some of their derivatives...
   */

  dbl wt_func;


  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM]; 
  int v_s[MAX_MODES][DIM][DIM]; 
  int v_g[DIM][DIM]; 

  dbl s[DIM][DIM];         /* stress tensor */
  dbl s_dot[DIM][DIM];     /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];         /* velocity gradient tensor */
  dbl gt[DIM][DIM];        /* transpose of velocity gradient tensor */


  /* dot product tensors */

  dbl s_dot_s[DIM][DIM]; 
  dbl s_dot_g[DIM][DIM];
  dbl g_dot_s[DIM][DIM];
  dbl s_dot_gt[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  SARAMITO_DEPENDENCE_STRUCT d_saramito_struct;
  SARAMITO_DEPENDENCE_STRUCT *d_saramito = &d_saramito_struct;

  // todo: will want to parse necessary parameters... for now hard code
  const bool saramitoEnabled = (vn->ConstitutiveEquation == SARAMITO_OLDROYDB ||
				vn->ConstitutiveEquation == SARAMITO_PTT      ||
				vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  dbl saramitoCoeff = 1.;


  dbl d_mup_dv_pj; 
  dbl d_mup_dmesh_pj; 

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* constitutive equation parameters */
  dbl alpha;     /* This is the Geisekus mobility parameter */
  dbl lambda=0;    /* polymer relaxation constant */
  dbl ucwt, lcwt; /* Upper convected derviative weight, Lower convected derivative weight */
  dbl eps;       /* This is the PTT elongation parameter */
  dbl Z=1.0;         /* This is the factor appearing in front of the stress tensor in PTT */
  dbl dZ_dtrace =0.0;

  /* advective terms are precalculated */
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

  dbl trace=0.0; /* trace of the stress tensor */

  /* SUPG variables */
  dbl h_elem=0, h_elem_inv=0, h_elem_deriv=0;
  dbl supg=0;

  if(vn->evssModel == EVSS_GRADV)
    {
      evss_gradv = 1;
    }
  
  status = 0;

  eqn   = R_STRESS11;			

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			/* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void) stress_eqn_pointer(v_s);
  (void) stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32; 
  v_g[2][2] = VELOCITY_GRADIENT33; 


  /*
   * Field variables...
   */
  
  for ( a=0; a<dim; a++)
    {
      v[a] = fv->v[a];

      /* note, these are zero for steady calculations */
      if (  pd->TimeIntegration != STEADY &&  pd->v[MESH_DISPLACEMENT1+a] )
	{
	  x_dot[a] = fv_dot->x[a];
	}
      else
	{
	  x_dot[a] = 0.;
	}
    }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  grad_v[a][b] = fv->grad_v[a][b];
	}
    }

  /* load up shearrate tensor based on velocity */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  if(evss_gradv)
	    {
	      g[a][b] = fv->grad_v[a][b];
	      gt[a][b] = fv->grad_v[b][a];
	    }
	  else
	    {
	      g[a][b] = fv->G[a][b];
	      gt[b][a]  = g[a][b];
	    }
	}
    }
  
  if( vn->wt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( vn->wt_funcModel == SUPG)
    {
      supg = vn->wt_func;
    }


  if(supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p]*vcent[p]*h[p];
	}
      h_elem = sqrt(h_elem)/2.;
      if(h_elem == 0.) 
	{
	  h_elem_inv=1.;
	}
      else
	{
	  h_elem_inv=1./h_elem;
	}
	
    }
  /* end Petrov-Galerkin addition */

  /*  shift factor  */
  if( pd->e[TEMPERATURE])
    {
      if(vn->shiftModel == CONSTANT)
	{
	  at = vn->shift[0];
	  for( j=0 ; j<ei->dof[TEMPERATURE] ; j++)
	    {
	      d_at_dT[j]=0.;
	    }
	}
      else if(vn->shiftModel == MODIFIED_WLF)
	{
	  wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
	  if(wlf_denom != 0.)
	    {
	      at=exp(vn->shift[0]*(mp->reference[TEMPERATURE]-fv->T)/wlf_denom);
	      for( j=0 ; j<ei->dof[TEMPERATURE] ; j++)
		{
		  d_at_dT[j]= -at*vn->shift[0]*vn->shift[1]
		    /(wlf_denom*wlf_denom)*bf[TEMPERATURE]->phi[j];
		}
	    }
	  else
	    { 
	      at = 1.;
	    } 
	  for( j=0 ; j<ei->dof[TEMPERATURE] ; j++)
	    {
	      d_at_dT[j]=0.;
	    }
	}
    }
  else
    {
      at = 1.;
    }

  /* Begin loop over modes */
  for ( mode=0; mode<vn->modes; mode++)
    {
      
      load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);
      
      /* precalculate advective terms of form (v dot del tensor)*/

      trace = 0.0;

      for ( a=0; a<VIM; a++)
	{
	  trace += s[a][a];
	  for ( b=0; b<VIM; b++)
	    {
	      v_dot_del_s[a][b] = 0.;
	      x_dot_del_s[a][b] = 0.; 
	      for ( q=0; q<dim; q++)
		{
		  v_dot_del_s[a][b] +=  v[q] * grad_s[q][a][b];
		  x_dot_del_s[a][b] +=  x_dot[q] * grad_s[q][a][b];
		} 
	    }
	}
      
      /*
       * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
       */
      
      /* get polymer viscosity */
      mup = viscosity(ve[mode]->gn, gamma, d_mup);

      if(saramitoEnabled == TRUE)
	{
	  saramitoCoeff = compute_saramito_model_terms(s, ve[mode]->gn->tau_y, ve[mode]->gn->fexp, d_saramito);
	}
      else
	{
	  saramitoCoeff = 1.;
	  d_saramito->tau_y = 0;
			
	  for(int i=0; i<VIM; ++i){
	    for(int j=0; j<VIM; ++j){
	      d_saramito->s[i][j] = 0;
	    }
	  }
	}

      /* get Geisekus mobility parameter */
      alpha = ve[mode]->alpha;
      
      /* get time constant */
      if(ve[mode]->time_constModel == CONSTANT)
	{
	  lambda = ve[mode]->time_const;
	}
      else if(ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW)
	{
	  lambda = mup/ve[mode]->time_const;
	}

      ucwt = 1.0 - ve[mode]->xi / 2.0 ;
      lcwt = ve[mode]->xi / 2.0 ;
      
      eps = ve[mode]->eps;

      Z = exp( eps*lambda*trace/mup ); dZ_dtrace = Z*eps*lambda/mup ;

      /* get tensor dot products for future use */
      
      if(DOUBLE_NONZERO(alpha)) (void) tensor_dot(s, s, s_dot_s, VIM);

      if( ucwt != 0. )
	{
	  (void) tensor_dot(s, g, s_dot_g, VIM);
	  (void) tensor_dot(gt, s, gt_dot_s, VIM);
	}

      if( lcwt != 0.) 
	{
	  (void) tensor_dot(s, gt, s_dot_gt, VIM);
	  (void) tensor_dot(g, s, g_dot_s, VIM);
	}
      /*
       * Residuals_________________________________________________________________
       */
      
      if ( af->Assemble_Residual )
	{
	  /*
	   * Assemble each component "ab" of the polymer stress equation...
	   */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      eqn = R_s[mode][a][b];
		      
		      /*
		       * In the element, there will be contributions to this many equations
		       * based on the number of degrees of freedom...
		       */
		      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  mass = 0.;
			  
			  if ( pd->TimeIntegration != STEADY )
			    {
			      if ( pd->e[eqn] & T_MASS )
				{
				  mass = s_dot[a][b];
				  mass *= wt_func * at * lambda * det_J * wt;
				  mass *= h3;
				  mass *= pd->etm[eqn][(LOG2_MASS)];
				}
			    }
			  
			  advection = 0.;
			  if ( pd->e[eqn] & T_ADVECTION )
 			    {
			      if(DOUBLE_NONZERO(lambda))
				{
				  
				  advection +=  v_dot_del_s[a][b]  -  x_dot_del_s[a][b];
				  if( ucwt != 0.) advection -= ucwt*(gt_dot_s[a][b] + s_dot_g[a][b]);
				  if( lcwt != 0.) advection += lcwt*(s_dot_gt[a][b] + g_dot_s[a][b]);

				  advection *= wt_func * at * lambda *det_J * wt * h3;
				  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
				}     
			    }
			  
			  diffusion = 0.;
			  if ( pd->e[eqn] & T_DIFFUSION )
			    {
			      /* add SU term in here when appropriate */
			      diffusion += 0.;
			      diffusion *= -wt_func * det_J * wt * h3;
			      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			    }
			  
			  /*
			   * Source term...
			   */
			  
			  source = 0.;
			  if ( pd->e[eqn] & T_SOURCE )
			    {
			      // consider whether saramitoCoeff should multiply here
			      source +=  saramitoCoeff * Z* s[a][b] - at * mup * ( g[a][b] +  gt[a][b]);
			      
			      if(DOUBLE_NONZERO(alpha))
				{
				  source1 = ( s_dot_s[a][b]/mup);
				  
				  source1 *= alpha * lambda * saramitoCoeff;
				  source  += source1;
				}

			      source *= wt_func * det_J * h3 * wt;
			      
			      source *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
			  
			  /*
			   * Add contributions to this residual (globally into Resid, and 
			   * locally into an accumulator)
			   */
			  
			  lec->R[upd->ep[eqn]][i] += 
			    mass + advection + diffusion + source;
			}
		    }
		}
	    }
	}
      
      /*
       * Jacobian terms...
       */
      
      if ( af->Assemble_Jacobian )
	{
	  dbl R_source, R_advection; /* Places to put the raw residual portions 
					instead of constantly recalcing them */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      eqn = R_s[mode][a][b];
		      peqn = upd->ep[eqn];

		      R_advection =  v_dot_del_s[a][b]  -  x_dot_del_s[a][b];
		      if( ucwt != 0.) R_advection -= ucwt*(gt_dot_s[a][b] + s_dot_g[a][b]);
		      if( lcwt != 0.) R_advection += lcwt*(s_dot_gt[a][b] + g_dot_s[a][b]);

		      R_source =   Z*s[a][b];
			      
		      if(DOUBLE_NONZERO(alpha)) R_source  +=  alpha * lambda*( s_dot_s[a][b]/mup);
		      R_source *= saramitoCoeff;
		      R_source += - at * mup * ( g[a][b] +  gt[a][b]);
		      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  
			  /*
			   * Set up some preliminaries that are needed for the (a,i)
			   * equation for bunches of (b,j) column variables...
			   */
			  
			  /* 
			   * J_S_T
			   */
			  
			  var = TEMPERATURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  mass = 0.;
			  
				  if ( pd->TimeIntegration != STEADY )
				    {
				      if ( pd->e[eqn] & T_MASS )
					{
				  	  mass = s_dot[a][b];
				  	  mass *= wt_func * d_at_dT[j] * lambda * det_J * wt;
				  	  mass *= h3;
				  	  mass *= pd->etm[eqn][(LOG2_MASS)];
					}
				    }
			  
				  advection = 0.;
				  if ( pd->e[eqn] & T_ADVECTION )
				    {
			      	      if(DOUBLE_NONZERO(lambda))
					{
				  
					  advection +=  v_dot_del_s[a][b]  -  x_dot_del_s[a][b];
					  if( ucwt != 0.) advection -= ucwt*(gt_dot_s[a][b] + s_dot_g[a][b]);
					  if( lcwt != 0.) advection += lcwt*(s_dot_gt[a][b] + g_dot_s[a][b]);

					  advection *= wt_func * d_at_dT[j] * lambda *det_J * wt * h3;
					  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					}     
				    }
			  
				  source    = 0.;
				  source1    = 0.;
				  if ( pd->e[eqn] & T_SOURCE )
				    {
				      source = -(g[a][b] +  gt[a][b])
					*(at*d_mup->T[j]+mup*d_at_dT[j]);
				      
			              if(DOUBLE_NONZERO(alpha))
					{
					  source1 -= s_dot_s[a][b]/(mup*mup)*d_mup->T[j];
					  source1 *= lambda * alpha * saramitoCoeff ;
					  source  += source1;
					}
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    mass + advection + source;
				}
			    }
			  
			  /*
			   * J_S_v
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = VELOCITY1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_mup_dv_pj  = d_mup->v[p][j];
				      
				      mass = 0.;

				      
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      if(supg!=0.)
						{
						  mass = supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
						  
						  for(w=0;w<dim;w++)
						    {
						      mass += supg * vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
							*v[w]*bf[eqn]->grad_phi[i] [w];
						    }
						  
						  mass *=  s_dot[a][b];
						}
					      
					      mass *= pd->etm[eqn][(LOG2_MASS)] * at * lambda * det_J * wt * h3;
					    }
					  
					}
				      
				      advection = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  if(DOUBLE_NONZERO(lambda))
					    {
					      advection_a = phi_j * 
						(grad_s[p][a][b]);
					      
					      advection_a *=  wt_func;
					      
					      advection_b = 0.;
					      /* Petrov-Galerkin term */
					      if(supg !=0.)
						{
						 
						  advection_b =  supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
						  for( w =0; w<dim; w++ )
						    {
						      advection_b += supg*vcent[p]*h[p]*dvc_dnode[p][j]*h_elem_inv/4.
							*v[w]*bf[eqn]->grad_phi[i] [w];
						    }
						  
						  advection_b *= R_advection;
						}
					      
					      advection_c = 0.;
					      if(evss_gradv)
						{
						  if(pd->CoordinateSystem != CYLINDRICAL)
						    {
						      if( ucwt != 0)
							{
							  for( k=0; k<VIM; k++)
							    {
							      advection_c -= ucwt*(bf[VELOCITY1+a]->grad_phi_e[j][p][k][a]*s[k][b] +
										   bf[VELOCITY1+b]->grad_phi_e[j][p][k][b]*s[a][k]);
							    }
							}
						      if( lcwt != 0.)
							{
							  for( k=0; k<VIM; k++)
							    {
							      advection_c += lcwt*(bf[VELOCITY1+b]->grad_phi_e[j][p][b][k]*s[a][k] +
										   bf[VELOCITY1+a]->grad_phi_e[j][p][a][k]*s[k][b]);
							    }
							}
						    }
						  else
						    {
						      if( ucwt != 0)
							{
							  for( k=0; k<VIM; k++)
							    {
							      advection_c -= ucwt*(bf[VELOCITY1]->grad_phi_e[j][p][k][a]*s[k][b] +
										   bf[VELOCITY1]->grad_phi_e[j][p][k][b]*s[a][k]);
							    }
							}
						      if( lcwt != 0.)
							{
							  for( k=0; k<VIM; k++)
							    {
							      advection_c += lcwt*(bf[VELOCITY1]->grad_phi_e[j][p][b][k]*s[a][k] +
										   bf[VELOCITY1]->grad_phi_e[j][p][a][k]*s[k][b]);
							    }
							}
						    }
						  advection_c *= wt_func;
						}
						     
					      advection = advection_a +  advection_b + advection_c;
					      advection *= at * lambda * det_J * wt *h3;
					      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					    }
					}
				      
				      diffusion = 0.;
				      
				      if ( pd->e[eqn] & T_DIFFUSION )
					{
					  /* add SU term in here when appropriate */
					  
					  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					}
				      
				      source    = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_c =  -at * d_mup_dv_pj * ( g[a][b] +  gt[a][b]);
                                          if(evss_gradv)
                                            {
					      if(pd->CoordinateSystem != CYLINDRICAL)
						{
						  source_c -= at * mup * (bf[VELOCITY1+a]->grad_phi_e[j][p][a][b] +
									  bf[VELOCITY1+b]->grad_phi_e[j][p][b][a]);
						}
					      else
						{
						  source_c -= at * mup * (bf[VELOCITY1]->grad_phi_e[j][p][a][b] +
									  bf[VELOCITY1]->grad_phi_e[j][p][b][a]);
						}
                                            }
					  source_c *= wt_func;
					  
					  source_a = 0.;
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source_a = -s_dot_s[a][b]/(mup*mup);
					      source_a *= wt_func * saramitoCoeff * alpha * lambda * d_mup_dv_pj;
					    }
					  
					  source_b = 0.;
					  if(supg != 0.)
					    {
					      source_b = supg * h_elem* phi_j*bf[eqn]->grad_phi[i][p];
					      
					      for(w=0;w<dim;w++)
						{
						  source_b += supg*vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
						    *v[w]*bf[eqn]->grad_phi[i] [w];
						}
					      
					      source_b *= R_source;
					    }
					  
					  source = source_a + source_b + source_c;
					  source *=  det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					}
				  
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + diffusion + source;
				    }
				}
			    }
			  
			  /*
			   * J_S_c
			   */
			  var = MASS_FRACTION;
			  if ( pd->v[var] )
			    {
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  for ( w=0; w<pd->Num_Species_Eqn; w++)
				    {
				      
				      source    = 0.;
			  
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a =   -at * d_mup->C[w][j] * ( g[a][b] +  gt[a][b]);
					  
					  source_b = 0.;
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source_b -= s_dot_s[a][b]/(mup*mup);
					      source_b *= alpha * lambda * saramitoCoeff * d_mup->C[w][j];
					    }
					  source = source_a + source_b;
					  source *= wt_func * det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					  
					}
				      
				      if ( w > 1 )
					{
					  EH(-1, "Need more arrays for each species.");
					}
				      
				      lec->J[peqn][MAX_PROB_VAR + w][i][j] +=
					source;
				    }
				}
			    }
			  
			  /*
			   * J_S_P
			   */
			  var = PRESSURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  source    = 0.;
				  if ( pd->e[eqn] & T_SOURCE )
				    {
				      source_a +=  -at * d_mup->P[j] * ( g[a][b] +  gt[a][b]);

				      source_b = 0.;
			              if(DOUBLE_NONZERO(alpha))
					{
					  source_b -= ( s_dot_s[a][b]/(mup*mup));
					  source_b *= d_mup->P[j] * alpha * lambda * saramitoCoeff;
					}
				      source  = source_a + source_b;
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    source;
				}		      
			    }
			  
			  /*
			   * J_S_d
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
				      dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
				      d_mup_dmesh_pj  = d_mup->X [p][j];
				      if(supg!=0.)
					{
					  h_elem_deriv = 0.;
					  for( q=0; q<dim; q++ )
					    {
					      h_elem_deriv += 
						hh[q][p]*vcent[q]*vcent[q]*dh_dxnode[q][j]*h_elem_inv/4.;
					    } 
					}
				      
				      mass = 0.;
				      mass_a = 0.;
				      mass_b = 0.;
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      mass_a = s_dot[a][b];
					      mass_a *= wt_func * ( d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					      
					      if(supg != 0.)
						{
						  for( w=0; w<dim; w++ )
						    {
						      mass_b += supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
									+ h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						    }
						  mass_b*= s_dot[a][b] * h3 * det_J;
						}
					      
					      mass = mass_a + mass_b;
					      mass *= at * lambda * wt * pd->etm[eqn][(LOG2_MASS)];
					    }
					}
				      
				      advection   = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  if(DOUBLE_NONZERO(lambda))
					    {
					      /*
					       * Four parts: 
					       *    advection_a = 
					       *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					       *
					       *    advection_b = 
					       *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
					       *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					       *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
					       *
					       * For unsteady problems, we have an 
					       * additional term
					       *
					       *    advection_c = 
					       *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
					       */
					      
					      advection_a =  R_advection;

					      advection_a *= wt_func *(  d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					      
					      d_vdotdels_dm = 0.;
					      for ( q=0; q<dim; q++)
						{
						  d_vdotdels_dm += (v[q]-x_dot[q]) * d_grad_s_dmesh[q][a][b] [p][j];
						}
					  
					      advection_b =  d_vdotdels_dm;
					      advection_b *=  wt_func *det_J * h3;
					  
					      advection_c = 0.;	
					      if ( pd->TimeIntegration != STEADY )
						{
						  if ( pd->e[eqn] & T_MASS )
						    {
						      d_xdotdels_dm = (1.+2.*tt) * phi_j/dt 
							* grad_s[p][a][b];
						  
						      advection_c -= d_xdotdels_dm;
						      
						      advection_c *= wt_func * h3 * det_J;
						    }
						}
					      
					      advection_d = 0.;	
					      if(supg!=0.)
						{
						  for( w=0; w<dim; w++ )
						    {
						      advection_d+= supg *(h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
									   + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						    }
						  
						  advection_d *= ( R_advection )* det_J * h3;
						}
					  
					      advection = advection_a + advection_b + advection_c + advection_d;
					      
					      advection *=  wt * at * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
					      
					    }
					}
				      
				      /*
				       * Diffusion...
				       */
				      
				      diffusion = 0.;
				      
				      if ( pd->e[eqn] & T_DIFFUSION )
					{
					  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					}
				      
				      /*
				       * Source term...
				       */
				      
				      source = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a =  R_source;
					  source_b = -at * (g[a][b] +  gt[a][b]);
					  
			                  if(DOUBLE_NONZERO(alpha))
					    {
					      source_b += -s_dot_s[a][b]/(mup*mup) * alpha * lambda * saramitoCoeff;
					    }
					  
					  source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
					  
					  source_b *= wt_func * det_J * h3 * d_mup_dmesh_pj ;
					  
					  source_c = 0.;
					  if(supg !=0.)
					    {
					      for( w=0;w<dim;w++)
						{
						  source_c+= supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
								     + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] ); 
						}
					      source_c *= R_source * det_J * h3;
					    }
					  
					  source  = source_a + source_b + source_c;
					  
					  source *=  wt * pd->etm[eqn][(LOG2_SOURCE)];
				      
					}
				      
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + diffusion + source;
				    }
				}
			    }
		      
			  /*
			   * J_S_G
			   */
			  if(evss_gradv == 0)
			    {
			      for ( p=0; p<VIM; p++)
				{
				  for ( q=0; q<VIM; q++)
				    {
				      var =  v_g[p][q];
				  
				      if ( pd->v[var] )
					{
					  pvar = upd->vp[var];
					  for ( j=0; j<ei->dof[var]; j++)
					    {
					      phi_j = bf[var]->phi[j];
					      advection   = 0.;
					      advection_a   = 0.;
					      if ( pd->e[eqn] & T_ADVECTION )
						{
						  if(DOUBLE_NONZERO(lambda))
						    {
						  
						      advection -=  ucwt * (s[p][b] * (double)delta(a,q) + s[a][p] * (double)delta(b,q));
						      advection +=  lcwt * (s[a][q] * (double)delta(p,b) + s[q][b] * (double)delta(a,p));
						  
						      advection *=  phi_j* h3 * det_J ;
						  
						      advection *= wt_func * wt * at * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
						    } 
						}
					  
					      /*
					       * Diffusion...
					       */
					  
					      diffusion = 0.;
					  
					      if ( pd->e[eqn] & T_DIFFUSION )
						{
						  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
						}
					  
					      /*
					       * Source term...
					       */
					  
					      source = 0.;		      
					  
					      if ( pd->e[eqn] & T_SOURCE )
						{
						  source =  -at * mup *  phi_j *( (double)delta(a,p)*(double)delta(b,q) +  (double)delta(b,p)*(double)delta(a,q));
						  source *= det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
						}
					  
					      lec->J[peqn][pvar][i][j] +=
						advection + diffusion + source;
					    }
					}
				    }
				}
			    }
		      
			  
			  /*
			   * J_S_S
			   */
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  var =  v_s[mode][p][q];
				  
				  if ( pd->v[var] )
				    {
				      pvar = upd->vp[var];
				      for ( j=0; j<ei->dof[var]; j++)
					{
					  phi_j = bf[var]->phi[j];
					  mass = 0.;
					  if ( pd->TimeIntegration != STEADY )
					    {
					      if ( pd->e[eqn] & T_MASS )
						{
						  mass = (1.+2.*tt) * phi_j/dt * (double)delta(a,p) * (double)delta(b,q); 
						  mass *= h3 * det_J;
						  mass *= wt_func * at * lambda * wt * pd->etm[eqn][(LOG2_MASS)];
						}
					    }
					  
					  advection   = 0.;
					  
					  if ( pd->e[eqn] & T_ADVECTION )
					    {
					      if(DOUBLE_NONZERO(lambda))
						{
						  if((a == p) && (b == q))
						    {
						      for( r=0; r<dim; r++)
							{
							  advection +=  (v[r]-x_dot[r])*  bf[var]->grad_phi[j][r];
							}
						    }
						  advection -=  phi_j*ucwt * (gt[a][p] * (double)delta(b,q) + g[q][b] * (double)delta(a,p));
						  advection +=  phi_j*lcwt * (gt[q][b] * (double)delta(p,a) + g[a][p] * (double)delta(q,b));
						  
						  advection *=  h3 * det_J ;
						  
						  advection *= wt_func * wt * at * lambda * pd->etm[eqn][(LOG2_ADVECTION)];
						}
					    }
					  
					  /*
					   * Diffusion...
					   */
					  
					  diffusion = 0.;
					  
					  if ( pd->e[eqn] & T_DIFFUSION )
					    {
					      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
					    }
					  
					  /*
					   * Source term...
					   */
					  
					  source = 0.;		      
					  
					  if ( pd->e[eqn] & T_SOURCE )
					    {
					      source_a  =  Z * (double)delta(a,p) * (double)delta(b,q);
 					      if( p == q) source_a +=  s[a][b] * dZ_dtrace; 
					      source_a *= saramitoCoeff;
					      // sensitivities for saramito model:
                                              if (p <= q) {
                                                source_a +=  d_saramito->s[p][q] * s[a][b] * Z;
                                              }
		      			
					      source_b  =0.;
			                      if(DOUBLE_NONZERO(alpha))
						{
						  source_b  = alpha * lambda * saramitoCoeff * 
						    (s[q][b] * (double)delta(a,p) + s[a][p] * (double)delta(b,q))/mup; 
                                                  if (p <= q) {
						    source_b += d_saramito->s[p][q] * alpha * lambda*( s_dot_s[a][b]/mup);
                                                  }
						}
					      
					      source  = source_a + source_b;
					      
					      source *= phi_j * det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					      
					    }
					  
					  lec->J[peqn][pvar][i][j] +=
					    mass + advection + diffusion + source;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	} /* End Assemble Jacobian */
    } /* End loop over modes */
  
  return(status);
}


/*
 * This routine assembles the stress with a log-conformation tensor formulation.
 */
int 
assemble_stress_log_conf(dbl tt,	 
		     dbl dt,	              
		     dbl h[DIM], 
		     dbl hh[DIM][DIM], 
		     dbl dh_dxnode[DIM][MDE],
		     dbl vcent[DIM], 
		     dbl dvc_dnode[DIM][MDE])
{
  int dim, p, q, a, b, w;
  int eqn, siz;

  int i, j, status, mode;
  int logc_gradv = 0;
  dbl v[DIM];			       
  dbl x_dot[DIM];	      
  dbl h3;		        

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];                
  dbl det_J;                            

  dbl mass;			        
  dbl advection;
  dbl source;
  dbl source_term1[DIM][DIM];

  dbl wt_func;
  dbl wt;
  dbl tmp1[DIM][DIM],tmp2[DIM][DIM],tmp3[DIM][DIM];
  dbl advection_term1[DIM][DIM];

  //Variables for stress, velocity gradient
  int R_s[MAX_MODES][DIM][DIM]; 
  int v_s[MAX_MODES][DIM][DIM]; 
  dbl s[DIM][DIM], exp_s[DIM][DIM];
  dbl s_dot[DIM][DIM];    
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl gt[DIM][DIM];

  //Polymer viscosity
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  //Temperature shift
  dbl at = 0.0;
  dbl wlf_denom;

  //Consitutive prameters
  dbl alpha;     
  dbl lambda=0;    
  dbl d_lambda;
  dbl eps;      
  dbl Z=1.0;        

  // Decomposition of velocity vector
  dbl M1[DIM][DIM];
  dbl eig_values[DIM];
  dbl R1[DIM][DIM];
  dbl R1_T[DIM][DIM];
  dbl Rt_dot_gradv[DIM][DIM];
  dbl D[DIM][DIM];
  dbl D_dot_D[DIM][DIM];

  //Advective terms
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  //Trace of stress
  dbl trace=0.0; 

  //SUPG terms
  dbl h_elem=0;
  dbl supg=0;

  status = 0;
  if(vn->evssModel == LOG_CONF_GRADV)
    {
      logc_gradv = 1;
    }
  
  eqn   = R_STRESS11;			
  //Check if we are actually needed
  if(!pd->e[eqn])
    {
      return(status);
    }

  dim   = pd->Num_Dim;
  wt = fv->wt;
  det_J = bf[eqn]->detJ;
  h3 = fv->h3;		      

  //Load pointers
  (void) stress_eqn_pointer(v_s);
  (void) stress_eqn_pointer(R_s);


  memset( s, 0, sizeof(double)*DIM*DIM);
  memset( exp_s, 0, sizeof(double)*DIM*DIM);
  
  //Load up field variables
  for(a=0; a<dim; a++)
    {
      //Velocity
      v[a] = fv->v[a];
      //
      if(pd->TimeIntegration!=STEADY && pd->v[MESH_DISPLACEMENT1+a])
	{
	  x_dot[a] = fv_dot->x[a];
	}
      else
	{
	  x_dot[a] = 0.0;
	}
    }

  //Velocity gradient
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  grad_v[a][b] = fv->grad_v[a][b];
	}
    }

  //Shear rate
  for(a=0; a<VIM; a++)
    {
      for(b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }

    // Velocity gradient projection
  for (a=0; a<VIM; a++)
    {
      for (b=0; b<VIM; b++)
	{
	  gt[a][b] = fv->G[b][a];
	}
    }

  if(vn->wt_funcModel == GALERKIN)
    {
      supg = 0.0;
    }
  else if( vn->wt_funcModel == SUPG)
    {
      supg = vn->wt_func;
    }


  if(supg!=0.0)
    {
      h_elem = 0.0;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p]*vcent[p]*h[p];
	}
      h_elem = sqrt(h_elem)/2.0;

    }

  //Shift factor
  if(pd->e[TEMPERATURE])
    {
      if(vn->shiftModel == CONSTANT)
	{
	  at = vn->shift[0];
	}
      else if(vn->shiftModel == MODIFIED_WLF)
	{
	  wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
	  if(wlf_denom!= 0.0)
	    {
	      at=exp(vn->shift[0]*(mp->reference[TEMPERATURE]-fv->T)/wlf_denom);
	    }
	  else
	    { 
	      at = 1.0;
	    } 

	}
    }
  else
    {
      at = 1.0;
    }

  //Loop over modes
  for(mode=0; mode<vn->modes; mode++)
    {
      //Load up constants and some pointers
      load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

      //Polymer viscosity
      mup = viscosity(ve[mode]->gn, gamma, d_mup);

      //Giesekus mobility parameter
      alpha = ve[mode]->alpha;
      
      //Polymer time constant
      if(ve[mode]->time_constModel == CONSTANT)
	{
	  lambda = ve[mode]->time_const;
	}
      else if(ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW)
	{
	  lambda = mup/ve[mode]->time_const;
	}

#ifdef ANALEIG_PLEASE
      analytical_exp_s(s, exp_s, eig_values, R1);
#else
      compute_exp_s(s, exp_s, eig_values, R1);
#endif

      /* Check to make sure eigenvalues are positive (negative eigenvalues will not
         work for log-conformation formulation). These eigenvalues are for the
         conformation tensor, not the log-conformation tensor. */
      if(eig_values[0] < 0. || eig_values[1] < 0. || (VIM > 2 && eig_values[2] < 0.))
	{
	  WH(-1, "Error: Negative eigenvalue for conformation tensor");
	  return -1;
	}
      
      memset(D, 0, sizeof(double)*DIM*DIM);
      D[0][0] = eig_values[0];
      D[1][1] = eig_values[1];
      if (VIM > 2) { D[2][2] = eig_values[2]; }
      (void) tensor_dot(D, D, D_dot_D, VIM);

      // Decompose velocity gradient

      memset(M1, 0, sizeof(double)*DIM*DIM);
      memset(R1_T, 0, sizeof(double)*DIM*DIM);

      for(i=0; i<VIM; i++)
        {
          for(j=0; j<VIM; j++)
            {
              R1_T[i][j] = R1[j][i];
            }
        }

      for(i=0; i<VIM; i++) 
        {    
          for(j=0; j<VIM; j++) 
            {
              Rt_dot_gradv[i][j] = 0.;  
              for(w=0; w<VIM; w++) 
                {
		  if(logc_gradv)
		    {
		      Rt_dot_gradv[i][j] += R1_T[i][w] * grad_v[j][w];
		    }
		  else
		    {
		      Rt_dot_gradv[i][j] += R1_T[i][w] * gt[w][j];
		    }
                }
            }
        }     

      for(i=0; i<VIM; i++) 
        {    
          for(j=0; j<VIM; j++) 
            {
              M1[i][j] = 0.;
              for(w=0; w<VIM; w++) 
                {
                  M1[i][j] += Rt_dot_gradv[i][w] * R1[w][j];
                }
            }
        }  

      //Predetermine advective terms
      trace = eig_values[0]+eig_values[1]; 
      if (VIM > 2) { trace += eig_values[2]; }
      
      for(a=0; a<VIM; a++)
      	{
      	  for(b=0; b<VIM; b++)
      	    {
      	      v_dot_del_s[a][b] = 0.0;
      	      x_dot_del_s[a][b] = 0.0;
      	      for(q=0; q<dim; q++)
      		{
      		  v_dot_del_s[a][b] +=  v[q]*grad_s[q][a][b];
      		  x_dot_del_s[a][b] +=  x_dot[q]*grad_s[q][a][b];
      		}
      	    }
      	}

      //PTT exponent
      eps  = ve[mode]->eps;

      //Exponential term for PTT
      Z = exp(eps*(trace - (double) VIM));

      siz = sizeof(double)*DIM*DIM;
      memset(tmp1, 0, siz);
      memset(tmp2, 0, siz);
      memset(tmp3, 0, siz);
      memset(advection_term1, 0, siz);
      memset(source_term1, 0, siz);

      for (a=0; a<VIM; a++)
        {
          for (b=0; b<VIM; b++)
            {
	      if ( a != b )
 	        {
                  d_lambda = eig_values[b]-eig_values[a];
                  if (DOUBLE_NONZERO(d_lambda))  
                    {
                      double eiglog_a = log(DBL_SMALL), eiglog_b = log(DBL_SMALL);
                      if (DOUBLE_NONZERO(eig_values[b]))
                         { eiglog_b = fmax(eiglog_b,log(eig_values[b]));}
                      if (DOUBLE_NONZERO(eig_values[a]))
                         { eiglog_a = fmax(eiglog_a,log(eig_values[a]));}
                      tmp1[a][b] += (eiglog_b - eiglog_a)/d_lambda;
                      tmp1[a][b] *= (eig_values[a]*M1[b][a] + eig_values[b]*M1[a][b]);
                    }
		  else
		    {
                      tmp1[a][b] += eig_values[b]*(M1[a][b] + M1[b][a]);
	            }
	        }
              if ( a == b )
                {
                  source_term1[a][b] += Z * (1.0 - D[a][a]) /lambda;
		  if(DOUBLE_NONZERO(alpha))
                    {
                      source_term1[a][b] += alpha*(2.0 * D[a][a] - 1.0 - D_dot_D[a][a])/lambda;
 	            }
                  source_term1[a][b] /= fmax(DBL_SMALL,eig_values[a]);
  	          source_term1[a][b] += 2.0*M1[a][a];
                }
            }
        }

      (void) tensor_dot(R1, tmp1, tmp2, VIM);
      (void) tensor_dot(tmp2, R1_T, advection_term1, VIM);
      (void) tensor_dot(R1, source_term1, tmp3, VIM);
      (void) tensor_dot(tmp3, R1_T, source_term1, VIM);

      if(af->Assemble_Residual)
	{
	  for(a=0; a<VIM; a++)
	    {
	      for(b=0; b<VIM; b++)
		{		  
		  if(a<=b)
		    {
		      eqn = R_s[mode][a][b];

		      for(i=0; i<ei->dof[eqn]; i++)
			{
			  wt_func = bf[eqn]->phi[i];
			  
			  //SUPG weighting, this is SUPG with s, not e^s
			  if(DOUBLE_NONZERO(supg))
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg*h_elem*v[w]*bf[eqn]->grad_phi[i][w];
				}
			    }

			  mass = 0.0;			  
			  if(pd->TimeIntegration!=STEADY)
			    {
			      if(pd->e[eqn] & T_MASS)
				{
				  mass  = s_dot[a][b];
				  mass *= wt_func*at*det_J*wt*h3;
				  mass *= pd->etm[eqn][(LOG2_MASS)];
				}
			    } 

			  advection = 0.;
			  if(pd->e[eqn] & T_ADVECTION)
 			    {  
			      advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
			      advection -= advection_term1[a][b];
			      advection *= wt_func*at*det_J*wt*h3;
			      advection *= pd->etm[eqn][(LOG2_ADVECTION)];			      
			    }

                          source = 0.0; 
                          if(pd->e[eqn] & T_SOURCE)
                            {
		              source -= source_term1[a][b];
                              source *= wt_func*det_J*h3*wt;     
                              source *= pd->etm[eqn][(LOG2_SOURCE)];
                            }
                          lec->R[upd->ep[eqn]][i] += mass + advection + source;
			}//i loop
		    }//if a<=b
		}// b loop
	    }//a loop
	}//if Residual

    }
  return(status);
}



/* this stress routine uses the viscoelastic equations to do a solid-fluid
* interaction problem in an Eulerian context with the level set denote the
* solid-fluid interface.
 */

int 
assemble_stress_level_set(dbl tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
			  dbl dt,	/* current time step size */
			  dbl h[DIM], /* coordinate scale factors */
			  dbl hh[DIM][DIM], /* coordinate scale factors */
			  dbl dh_dxnode[DIM][MDE],
			  dbl vcent[DIM], /* Average element velocity, which is 
					* the centroid velocity for Q2 and the
					* average of the vertices for Q1. It 
					* comes from the routine 
					* "element_velocity." */
			  dbl dvc_dnode[DIM][MDE])
{
  int dim, p, q, r, a, b, w;

  int eqn, var;
  int peqn, pvar;

  int i, j, mode;
  dbl v[DIM];			        /* Velocity field. */
  dbl x_dot[DIM];			/* current position field derivative wrt time. */
  dbl h3;		        	/* Volume element (scale factors). */
  dbl dh3dmesh_pj;	        	/* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];                  /* Shear-rate tensor based on velocity */
  dbl det_J;                            /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj;			/* for specific (p,j) mesh dof */

  dbl mass;			        /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;	
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl source;
  dbl source_a=0, source_b=0, source_c=0;

  /*
   * 
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals 
   * and some of their derivatives...
   */

  dbl wt_func;


  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM]; 
  int v_s[MAX_MODES][DIM][DIM]; 
  int v_g[DIM][DIM]; 

  dbl s[DIM][DIM];         /* stress tensor */
  dbl s_dot[DIM][DIM];     /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];         /* velocity gradient tensor */
  dbl gt[DIM][DIM];        /* transpose of velocity gradient tensor */


  /* dot product tensors */

  dbl s_dot_g[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  dbl d_mup_dv_pj; 
  dbl d_mup_dmesh_pj; 

  /* advective terms are precalculated */
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

/* SUPG variables */
  dbl h_elem=0, h_elem_inv=0, h_elem_deriv=0;
  dbl supg=0;

  int status = 0;

  double H_ls, delta_ls, normal_ls[MAX_PDIM];
  int near_ls;
  dbl Gmod;

  eqn   = R_STRESS11;			

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			/* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void) stress_eqn_pointer(v_s);
  (void) stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32; 
  v_g[2][2] = VELOCITY_GRADIENT33; 


  /*
   * Field variables...
   */
  
  for ( a=0; a<dim; a++)
    {
      v[a] = fv->v[a];

      /* note, these are zero for steady calculations */
      if (  pd->TimeIntegration != STEADY &&  pd->v[MESH_DISPLACEMENT1+a] )
	{
	  x_dot[a] = fv_dot->x[a];
	}
      else
	{
	  x_dot[a] = 0.;
	}
    }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  grad_v[a][b] = fv->grad_v[a][b];
	}
    }

/* load up shearrate tensor based on velocity */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  g[a][b]  = fv->G[a][b];
	  gt[b][a]  = g[a][b];
	}
    }
  
  if( vn->wt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( vn->wt_funcModel == SUPG)
    {
      supg = vn->wt_func;
    }


  if(supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p]*vcent[p]*h[p];
	}
      h_elem = sqrt(h_elem)/2.;
      if(h_elem == 0.) 
	{
	  h_elem_inv=1.;
	}
      else
	{
	  h_elem_inv=1./h_elem;
	}
	
    }
/* end Petrov-Galerkin addition */

  /* Get level set length scales, Heaviside function etc */
  level_set_interface(fv->F, fv->grad_F, ls->Length_Scale, 0,
                      &near_ls,
		      &H_ls, NULL, NULL, 
		      &delta_ls, NULL, NULL, 
		      normal_ls, NULL, NULL);


  /* Begin loop over modes */
  for ( mode=0; mode<vn->modes; mode++)
    {
      
      load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);
      
      /* precalculate advective terms of form (v dot del tensor)*/


      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      v_dot_del_s[a][b] = 0.;
	      x_dot_del_s[a][b] = 0.; 
	      for ( q=0; q<dim; q++)
		{
		  v_dot_del_s[a][b] +=  v[q] * grad_s[q][a][b];
		  x_dot_del_s[a][b] +=  x_dot[q] * grad_s[q][a][b];
		} 
	    }
	}
      
      /*
       * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
       */
      
      /* get polymer viscosity */
      mup = viscosity(ve[mode]->gn, gamma, d_mup);
      
      /* hardwire shear modulus for now */
      /* Gmod = 1000000.; */
      Gmod = mup/ve[mode]->time_const;


      /* get tensor dot products for future use */
      (void) tensor_dot(s, g, s_dot_g, VIM);
      (void) tensor_dot(gt, s, gt_dot_s, VIM);

      /* get tensor dot products for future use */
      

      /*
       * Residuals_________________________________________________________________
       */
      
      if ( af->Assemble_Residual )
	{
	  /*
	   * Assemble each component "ab" of the polymer stress equation...
	   */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      eqn = R_s[mode][a][b];
		      
		      /*
		       * In the element, there will be contributions to this many equations
		       * based on the number of degrees of freedom...
		       */
		      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  /* The mass and advection terms constitute the substantial 
			   * derivative and the upper convective derivative of the
			   * stress tensor. These terms are divided by the modulus, 
			   * G, and form the elastic solid response.
			   */
			  mass = 0.;
			  if ( pd->TimeIntegration != STEADY )
			    {
			      if ( pd->e[eqn] & T_MASS )
				{
				  mass = s_dot[a][b]* H_ls / Gmod;
				  mass *= wt_func * det_J * wt * h3;
				  mass *= pd->etm[eqn][(LOG2_MASS)];
				}
			    }
			  
			  advection = 0.;
			  if ( pd->e[eqn] & T_ADVECTION )
 			    {
			      advection +=  v_dot_del_s[a][b]  -  x_dot_del_s[a][b];
			      advection -= (gt_dot_s[a][b] + s_dot_g[a][b]);
			      
			      advection *= wt_func * H_ls / Gmod *det_J * wt * h3;
			      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
			    }
			  
			  
			  /*
			   * Source term...
			   */
			  
			  source = 0.;
			  if ( pd->e[eqn] & T_SOURCE )
			    {
			      /* This source term involves the Newtonian fluid
			       * viscosity expressed in the mixed formulation.
			       */
			      source =  s[a][b] * (1.-H_ls) / mup; 

			      /* This source term involves the rate of deformation
			       * tensor, which is used for both the solid and fluid
			       * formulations.
			       */

			      source -= ( g[a][b] +  gt[a][b]);
			      source *= wt_func * det_J * h3 * wt;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
			  
			  /*
			   * Add contributions to this residual (globally into Resid, and 
			   * locally into an accumulator)
			   */
			  
			  lec->R[upd->ep[eqn]][i] += 
			    mass + advection + source;
			}
		    }
		}
	    }
	}
      
      /*
       * Jacobian terms...
       */
      
      if ( af->Assemble_Jacobian )
	{
	  dbl R_source, R_advection; /* Places to put the raw residual portions 
					        instead of constantly recalcing them */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  if(a <= b) /* since the stress tensor is symmetric, only assemble the upper half */ 
		    {
		      eqn = R_s[mode][a][b];
		      peqn = upd->ep[eqn];

		      R_advection =  v_dot_del_s[a][b]  -  x_dot_del_s[a][b];
		      R_advection -= (gt_dot_s[a][b] + s_dot_g[a][b]);
		      R_advection *= H_ls / Gmod;
		     
		      R_source =   s[a][b] * (1.-H_ls) / mup -  ( g[a][b] +  gt[a][b] );
			      
		      for ( i=0; i<ei->dof[eqn]; i++)
			{
			  
			  wt_func = bf[eqn]->phi[i];
			  /* add Petrov-Galerkin terms as necessary */
			  if(supg!=0.)
			    {
			      for(w=0; w<dim; w++)
				{
				  wt_func += supg * h_elem*v[w]*bf[eqn]->grad_phi[i] [w];
				}
			    }
			  
			  
			  /*
			   * Set up some preliminaries that are needed for the (a,i)
			   * equation for bunches of (b,j) column variables...
			   */
			  
			  /* 
			   * J_S_T
			   */
			  
			  var = TEMPERATURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  source    = 0.;
				  if ( pd->e[eqn] & T_SOURCE )
				    {
				      source -=  s[a][b] * (1.-H_ls)/(mup*mup)* d_mup->T[j];
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    source;
				}
			    }
			  
			  /*
			   * J_S_v
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = VELOCITY1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_mup_dv_pj  = d_mup->v[p][j];
				      
				      mass = 0.;
				      
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      if(supg!=0.)
						{
						  mass = supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
						  
						  for(w=0;w<dim;w++)
						    {
						      mass += supg * vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
							*v[w]*bf[eqn]->grad_phi[i] [w];
						    }
						  
						  mass *=  s_dot[a][b] * H_ls / Gmod;
						}
					      
					      mass *= pd->etm[eqn][(LOG2_MASS)] * det_J * wt * h3;
					    }
					  
					}
				      
				      advection = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  advection_a = phi_j * 
					    (grad_s[p][a][b]) * H_ls / Gmod;
					  
					  advection_a *=  wt_func;
					  
					  advection_b = 0.;
					  /* Petrov-Galerkin term */
					  if(supg !=0.)
					    {
					      advection_b =  supg * h_elem*phi_j*bf[eqn]->grad_phi[i][p];
					      for( w =0; w<dim; w++ )
						{
						  advection_b += supg*vcent[p]*h[p]*dvc_dnode[p][j]*h_elem_inv/4.
						    *v[w]*bf[eqn]->grad_phi[i] [w];
						}
					      
					      advection_b *= R_advection;
					    }
					  
					  advection = advection_a +  advection_b;
					  advection *= det_J * wt *h3;
					  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					}
				      
				      source    = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a =  -d_mup_dv_pj * (s[a][b] * (1.-H_ls) / (mup*mup ));
					  source_a *= wt_func;
					  
					  source_b = 0.;
					  if(supg != 0.)
					    {
					      source_b = supg * h_elem* phi_j*bf[eqn]->grad_phi[i][p];
					      
					      for(w=0;w<dim;w++)
						{
						  source_b += supg*vcent[p]*dvc_dnode[p][j]*h[p]*h_elem_inv/4.
						    *v[w]*bf[eqn]->grad_phi[i] [w];
						}
					      
					      source_b *= R_source;
					    }
					  
					  source = source_a + source_b;
					  source *=  det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					}
				  
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + source;
				    }
				}
			    }
			  
			  /*
			   * J_S_c
			   */
			  var = MASS_FRACTION;
			  if ( pd->v[var] )
			    {
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  for ( w=0; w<pd->Num_Species_Eqn; w++)
				    {
				      
				      source    = 0.;
			  
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source =   -d_mup->C[w][j] * (s[a][b] * (1.-H_ls) / (mup*mup) );
					  source *= wt_func * det_J * wt * h3;
					  source *= pd->etm[eqn][(LOG2_SOURCE)];
					  
					}
				      
				      if ( w > 1 )
					{
					  EH(-1, "Need more arrays for each species.");
					}
				      
				      lec->J[peqn][MAX_PROB_VAR + w][i][j] +=
					source;
				    }
				}
			    }
			  
			  /*
			   * J_S_P
			   */
			  var = PRESSURE;
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  source    = 0.;
				  if ( pd->e[eqn] & T_SOURCE )
				    {
				      source =  -d_mup->P[j] * (s[a][b] * (1.-H_ls) / (mup*mup) );
				      source *= wt_func * det_J * wt * h3;
				      source *= pd->etm[eqn][(LOG2_SOURCE)];
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    source;
				}		      
			    }
			  
			  /*
			   * J_S_d
			   */
			  for ( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1+p;
			      if ( pd->v[var] )
				{
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++)
				    {
				      phi_j = bf[var]->phi[j];
				      d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
				      dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
				      d_mup_dmesh_pj  = d_mup->X [p][j];
				      if(supg!=0.)
					{
					  h_elem_deriv = 0.;
					  for( q=0; q<dim; q++ )
					    {
					      h_elem_deriv += 
						hh[q][p]*vcent[q]*vcent[q]*dh_dxnode[q][j]*h_elem_inv/4.;
					    } 
					}
				      
				      mass = 0.;
				      mass_a = 0.;
				      mass_b = 0.;
				      if ( pd->TimeIntegration != STEADY )
					{
					  if ( pd->e[eqn] & T_MASS )
					    {
					      mass_a = s_dot[a][b]* H_ls / Gmod;
					      mass_a *= wt_func * ( d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					      
					      if(supg != 0.)
						{
						  for( w=0; w<dim; w++ )
						    {
						      mass_b += supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
									+ h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						    }
						  mass_b*= s_dot[a][b]* H_ls / Gmod * h3 * det_J;
						}
					      
					      mass = mass_a + mass_b;
					      mass *= wt * pd->etm[eqn][(LOG2_MASS)];
					    }
					}
				      
				      advection   = 0.;
				      
				      if ( pd->e[eqn] & T_ADVECTION )
					{
					  /*
					   * Four parts: 
					   *    advection_a = 
					   *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					   *
					   *    advection_b = 
					   *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
					   *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
					   *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
					   *
					   * For unsteady problems, we have an 
					   * additional term
					   *
					   *    advection_c = 
					   *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
					   */
					  
					  advection_a =  R_advection;
					  
					  advection_a *= wt_func *(  d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
					  
					  d_vdotdels_dm = 0.;
					  for ( q=0; q<dim; q++)
					    {
					      d_vdotdels_dm += (v[q]-x_dot[q]) * d_grad_s_dmesh[q][a][b] [p][j];
					    }
					  
					  advection_b =  d_vdotdels_dm;
					  advection_b *=  H_ls / Gmod * wt_func *det_J * h3;
					  
					  advection_c = 0.;	
					  if ( pd->TimeIntegration != STEADY )
					    {
					      if ( pd->e[eqn] & T_MASS )
						{
						  d_xdotdels_dm = (1.+2.*tt) * phi_j/dt 
						    * grad_s[p][a][b];
						  
						  advection_c -= d_xdotdels_dm;
						  
						  advection_c *= H_ls / Gmod * wt_func * h3 * det_J;
						}
					    }
					  
					  advection_d = 0.;	
					  if(supg!=0.)
					    {
					      for( w=0; w<dim; w++ )
						{
						  advection_d+= supg *(h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
								       + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] );
						}
					      
					      advection_d *= ( R_advection )* det_J * h3;
					    }
					  
					  advection = advection_a + advection_b + advection_c + advection_d;
					  
					  advection *=  wt * pd->etm[eqn][(LOG2_ADVECTION)];
					  
					}
				      
				      /*
				       * Source term...
				       */
				      
				      source = 0.;
				      
				      if ( pd->e[eqn] & T_SOURCE )
					{
					  source_a =  R_source;
					  source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

					  source_b = - s[a][b] * (1.-H_ls) / (mup*mup);
					  source_b *= wt_func * det_J * h3 * d_mup_dmesh_pj ;
					  
					  source_c = 0.;
					  if(supg !=0.)
					    {
					      for( w=0;w<dim;w++)
						{
						  source_c+= supg * (h_elem*v[w]* bf[eqn]->d_grad_phi_dmesh[i][w] [p][j]
								     + h_elem_deriv * v[w]*bf[eqn]->grad_phi[i] [w] ); 
						}
					      source_c *= R_source * det_J * h3;
					    }
					  
					  source  = source_a + source_b + source_c;
					  
					  source *=  wt * pd->etm[eqn][(LOG2_SOURCE)];
				      
					}
				      
				      lec->J[peqn][pvar][i][j] +=
					mass + advection + source;
				    }
				}
			    }
		      
			  /*
			   * J_S_G
			   */
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  var =  v_g[p][q];
				  
				  if ( pd->v[var] )
				    {
				      pvar = upd->vp[var];
				      for ( j=0; j<ei->dof[var]; j++)
					{
					  phi_j = bf[var]->phi[j];
					  advection   = 0.;
					  if ( pd->e[eqn] & T_ADVECTION )
					    {
					      advection -=  (s[p][b] * (double)delta(a,q) + s[a][p] * (double)delta(b,q));
					      advection *=  phi_j* h3 * det_J ;
					      
					      advection *= wt_func * wt * H_ls / Gmod * pd->etm[eqn][(LOG2_ADVECTION)];
					    }

					  /*
					   * Source term...
					   */
					  
					  source = 0.;		      
					  
					  if ( pd->e[eqn] & T_SOURCE )
					    {
					      source =  -phi_j *( (double)delta(a,p)*(double)delta(b,q) +  (double)delta(b,p)*(double)delta(a,q));
					      source *= det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					    }
					  
					  lec->J[peqn][pvar][i][j] +=
					    advection + source;
					}
				    }
				}
			    }
		      
			  
			  /*
			   * J_S_S
			   */
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  var =  v_s[mode][p][q];
				  
				  if ( pd->v[var] )
				    {
				      pvar = upd->vp[var];
				      for ( j=0; j<ei->dof[var]; j++)
					{
					  phi_j = bf[var]->phi[j];
					  mass = 0.;
					  if ( pd->TimeIntegration != STEADY )
					    {
					      if ( pd->e[eqn] & T_MASS )
						{
						  mass = (1.+2.*tt) * phi_j/dt * (double)delta(a,p) * (double)delta(b,q); 
						  mass *= h3 * det_J;
						  mass *= wt_func * H_ls / Gmod * wt * pd->etm[eqn][(LOG2_MASS)];
						}
					    }
					  
					  advection   = 0.;
					  
					  if ( pd->e[eqn] & T_ADVECTION )
					    {
					      if((a == p) && (b == q))
						{
						  for( r=0; r<dim; r++)
						    {
						      advection +=  (v[r]-x_dot[r])*  bf[var]->grad_phi[j][r];
						    }
						}
					      
					      advection -=  phi_j*(gt[a][p] * (double)delta(b,q) + g[q][b] * (double)delta(a,p));
						  
					      advection *=  h3 * det_J ;
						  
					      advection *= wt_func * wt *  H_ls / Gmod * pd->etm[eqn][(LOG2_ADVECTION)];
					    }
					  
					  /*
					   * Source term...
					   */
					  
					  source = 0.;		      
					  
					  if ( pd->e[eqn] & T_SOURCE )
					    {
					      source  = phi_j * (double)delta(a,p) * (double)delta(b,q) * (1.-H_ls) / mup;  
		      
					      source *= det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					    }
					  
					  lec->J[peqn][pvar][i][j] +=
					    mass + advection + source;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  return(status);
}


int
assemble_gradient(dbl tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
		  dbl dt)	/* current time step size */
{
  int dim;
  int p, q, a, b;
  
  int eqn, var;
  int peqn, pvar;
  int i, j;
  int status;
  
  dbl h3;		        	/* Volume element (scale factors). */
  dbl dh3dmesh_pj;	        	/* Sensitivity to (p,j) mesh dof. */
  
  dbl grad_v[DIM][DIM];
  dbl g[DIM][DIM];                      /* velocity gradient tensor */
  
  dbl det_J;                            /* determinant of element Jacobian */
  
  dbl d_det_J_dmesh_pj;			/* for specific (p,j) mesh dof */
  
  dbl advection;	
  dbl advection_a, advection_b;
  dbl source;
  
  /*
   * 
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */
  
  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals 
   * and some of their derivatives...
   */
  
  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals 
   * and some of their derivatives...
   */
  
  dbl wt_func;
  
  
  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  
  dbl phi_j;
  
  dbl wt;
  
  /* Variables for stress */
  
  int R_g[DIM][DIM]; 
  int v_g[DIM][DIM]; 
  
  
  
  status = 0;
  
  /*
   * Unpack variables from structures for local convenience...
   */
  
  dim   = pd->Num_Dim;
  
  eqn   = R_GRADIENT11;			
  
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[eqn] )
    {
      return(status);
    }
  
  wt = fv->wt;
  
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  
  h3 = fv->h3;			/* Differential volume element (scales). */
  
  /* load eqn and variable number in tensor form */
  
  
  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32; 
  v_g[2][2] = VELOCITY_GRADIENT33; 
  
  
  R_g[0][0] = R_GRADIENT11;
  R_g[0][1] = R_GRADIENT12;
  R_g[1][0] = R_GRADIENT21;
  R_g[1][1] = R_GRADIENT22;
  R_g[0][2] = R_GRADIENT13;
  R_g[1][2] = R_GRADIENT23;
  R_g[2][0] = R_GRADIENT31;
  R_g[2][1] = R_GRADIENT32; 
  R_g[2][2] = R_GRADIENT33; 
  
  
  /*
   * Field variables...
   */
  
  
  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  grad_v[a][b] = fv->grad_v[a][b];
	}
    }
  
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  g[a][b]  = fv->G[a][b];
	}
    }
  
  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      /*
       * Assemble each component "ab" of the velocity gradient equation...
       */
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      eqn = R_g[a][b];
	      /*
	       * In the element, there will be contributions to this many equations
	       * based on the number of degrees of freedom...
	       */
	      
	      for ( i=0; i<ei->dof[eqn]; i++)
		{
		  
		  wt_func = bf[eqn]->phi[i];   /* add Petrov-Galerkin terms as necessary */
		  
		  advection = 0.;
		  
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      advection -= grad_v[a][b];
		      advection *= wt_func * det_J * wt * h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }
		  
		  /*
		   * Source term...
		   */
		  
		  source = 0;
		  
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source +=  g[a][b];    
		      source *= wt_func * det_J * h3 * wt;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }
		  
		  lec->R[upd->ep[eqn]][i] += 
		    advection  + source;      
		}
	    }
	}
    }  
  
  /*
   * Jacobian terms...
   */
  
  if ( af->Assemble_Jacobian )
    {
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      eqn = R_g[a][b];
	      peqn = upd->ep[eqn];
	      
	      for ( i=0; i<ei->dof[eqn]; i++)
		{
		  wt_func = bf[eqn]->phi[i];   /* add Petrov-Galerkin terms as necessary */
  
		  /*
		   * J_G_v
		   */
		  for ( p=0; p<dim; p++)
		    {
		      var = VELOCITY1+p;
		      if ( pd->v[var] )
			{
			  pvar = upd->vp[var];
			  for ( j=0; j<ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];
			      
			      advection = 0.;
			      
			      if ( pd->e[eqn] & T_ADVECTION )
				{
				  advection -= bf[var]->grad_phi_e[j][p][a][b];
				  advection *= wt_func * det_J * wt *h3;
				  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
				  
				}
			      
			      source    = 0.;
			      
			      lec->J[peqn][pvar][i][j] +=
				advection + source;
			    }
			}
		    }
		  
		  /*
		   * J_G_d
		   */
		  for ( p=0; p<dim; p++)
		    {
		      var = MESH_DISPLACEMENT1+p;
		      if ( pd->v[var] )
			{
			  pvar = upd->vp[var];
			  for ( j=0; j<ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];
			      
			      d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
			      
			      dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
			      
			      advection   = 0.;
			      
			      if ( pd->e[eqn] & T_ADVECTION )
				{
				  /*
				   * three parts: 
				   *    advection_a = 
				   *    	Int ( d(Vv)/dmesh h3 |Jv| )
				   *
				   *    advection_b = 
				   *  (i)	Int ( Vv h3 d(|Jv|)/dmesh )
				   *  (ii)      Int ( Vv dh3/dmesh |Jv|   )
				   */
				  
				  
				  advection_a =  -grad_v[a][b];
				  
				  advection_a *=  (  d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj );
				  
				  advection_b  = -fv->d_grad_v_dmesh[a][b] [p][j];
				  
				  advection_b *=  det_J * h3;
				  
				  advection = advection_a + advection_b;
				  
				  advection *= wt_func * wt * pd->etm[eqn][(LOG2_ADVECTION)];
				  
				}
			      
			      source = 0.;
			      
			      if ( pd->e[eqn] & T_SOURCE )
				{
				  source +=  g[a][b];
				  
				  source *=  d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj;
				  
				  source *= wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
				  
				}
			      
			      lec->J[peqn][pvar][i][j] +=
				advection + source;
			    }
			}
		    }
		  
		  /*
		   * J_G_G
		   */
		  
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{
			  var =  v_g[p][q];
			  
			  if ( pd->v[var] )
			    {
			      pvar = upd->vp[var];
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];
				  
				  source = 0.;		      
				  
				  if ( pd->e[eqn] & T_SOURCE )
				    {
				      if((a == p) && (b == q))
					{
					  source = phi_j  * det_J * h3 * wt_func * wt * pd->etm[eqn][(LOG2_SOURCE)];
					}
				    }
				  
				  lec->J[peqn][pvar][i][j] +=
				    source;
				}
			    }
			}
		      
		    }
		}  
	    }
	}
    }
  return(status);
}


int tensor_dot(dbl t1[DIM][DIM],
	       dbl t2[DIM][DIM],
	       dbl t1_dot_t2[DIM][DIM],
	       const int dim)
{
  int i,j,k;
  int status;
  dbl v1[DIM];
  dbl v2[DIM];

  for(k=0;k<dim;k++)
    {
      for(i=0;i<dim;i++)
	{
	  v1[i]=t1[k][i];
	}
      for(j=0;j<dim;j++)
	{
	  for(i=0;i<dim;i++)
	    {
	      v2[i]=t2[i][j];
	    }
	  t1_dot_t2[k][j]=vec_dot(dim,v1,v2);
	}
    }

  status=1;
  return(status);
}

dbl
vec_dot(const int n1,
	dbl *v1,
	dbl *v2)
{
  int i;
  dbl rc = 0.0;

  for (i = 0; i < n1; i++) {
    rc += *v1 * *v2;
    v1++; v2++;
  }
  return(rc);
}

void 
load_modal_pointers(int ve_mode, /* mode number */
		    dbl tt,
		    dbl dt,
		    dbl s[DIM][DIM], /* stress tensor for mode ve_mode */
		    dbl s_dot[DIM][DIM], /* stress tensor time derivative for mode ve_mode */
		    dbl grad_s[DIM][DIM][DIM], /* grad of stress tensor for mode ve_mode */
		    dbl d_grad_s_dm[DIM][DIM][DIM][DIM][MDE]) /* derivative of grad of stress tensor for mode ve_mode */

{
  int a,b,p,q;			/* indeces for dimensions */
  int j;			/* indeces for dofs */
  int var;
  int siz;

  siz = sizeof(double)*DIM*DIM*DIM*DIM*MDE;
  memset(d_grad_s_dm,0,siz);

  /* load up things we need in the assembly routine for each mode in turn*/

  /* put stress in a nice working array */
  
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      s[a][b] = fv->S[ve_mode][a][b];	  
      if (pd->TimeIntegration != STEADY)  {
	s_dot[a][b] = fv_dot->S[ve_mode][a][b];
      } else {
	s_dot[a][b] = 0.;
      }
    }
  } 
  
  for ( p=0; p<VIM; p++)
    {
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      grad_s[p][a][b] = fv->grad_S[ve_mode][p][a][b];
	    }
	}
    } 
  
  var = MESH_DISPLACEMENT1;
  if ( pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  for ( q=0; q<ei->ielem_dim; q++)
		    {
		      for ( j=0; j<ei->dof[var]; j++)
			{
			  d_grad_s_dm[p][a][b][q][j] = fv->d_grad_S_dmesh[ve_mode][p][a][b][q][j];
			}
		    }
		}
	    }
	}
    }

}


/******************************************************************************/
/* END of routine modal_esp_alloc */
/******************************************************************************/

int
assemble_surface_stress (Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
			 double x[],
			 struct Aztec_Linear_Solver_System *ams,
			 dbl x_update[],    /* last update for x vector */
			 double delta_t, /* current time step size */
			 double t_,	/* parameter to vary time integration from
					 * explicit (tt = 1) to implicit (tt = 0) */
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
  int a, b, p, q;    /* more counters */
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];

  int eqn, peqn;
  int var, pvar;
  int err;         /* status variable for functions */
  int ip_total, ip_total_fill;
  int found_it;

  /* Variables for stress */
  int R_s[MAX_MODES][DIM][DIM]; 
  int v_s[MAX_MODES][DIM][DIM];
  int S_map[MAX_MODES][DIM][DIM]; /* map var index to stress mode component */

  int ibc, ins, side_index;
  int table_ibc[MAX_MODES][DIM][DIM]; /* maps table boundary condition index for each stress mode */
  int num_zeros; 

  int doit = 1;

  int mode;

  dbl s[DIM][DIM], s_dot[DIM][DIM], grad_s[DIM][DIM][DIM];
  dbl stress_neighbor[MAX_SURF_GP][MAX_MODES][DIM][DIM];
  dbl stress_update_v[MAX_SURF_GP][MAX_MODES][DIM][DIM];
  dbl s_n[MAX_MODES][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl **x_neighbor;

  dbl *J_S_S_v[MAX_MODES][DIM][DIM][MDE][DIM][DIM][MDE];
  dbl phi_neighbor[MAX_SURF_GP][MDE];

  dbl advection;

  double phi_j, phi_i;
  dbl rhs;
  double ss, tt, uu;			/* Gaussian quadrature point locations  */
  double xi[DIM];             /* Local element coordinates of Gauss point. */
  dbl vdotn, vdotn_avg;
  dbl vdotn_norm;
  double wt;                  /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  dbl *phi_v=NULL;

  dbl alpha = 0.5;
    
  /***************************************************************************/

  /* load eqn and variable number in tensor form */
  err = stress_eqn_pointer(S_map);

  err = stress_eqn_pointer(v_s);
  err = stress_eqn_pointer(R_s);

  /* initialization of the neighbor stress pointers array */


  memset( J_S_S_v, 0,  MAX_MODES*DIM*DIM*MDE*DIM*DIM*MDE);

  /********************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)             */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                           */
  /********************************************************************************/
  /* Find out the number of surface quadrature points 
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_SURF, ielem_type);
  
  dim =  pd->Num_Dim;


  /* allocate space for x_neighbor */

  /*  x_neighbor = (double **) array_alloc(2, ip_total, DIM, sizeof(double)); */
  /*  manually allocate space to correct (seeming) misalignment for HPUX */

   x_neighbor = (double **) smalloc( ip_total*sizeof(double*) );
     for( i=0; i<DIM; i++ )
      {
        x_neighbor[i] = ( double * ) smalloc( DIM * sizeof(double) );
      }



/* If no neighbor element found, check for table boundary condition on inlet */

     for(mode=0; mode<vn->modes; mode++) 
       {
	 for(a=0; a<VIM; a++)
	   {
	     for(b=0; b<VIM; b++)
	       {
		 table_ibc[mode][a][b]=-1;
	       }
	   }
       }

     for(ibc=0; ibc<Num_BC; ibc++)
       { 
	 if(BC_Types[ibc].BC_Name == TABLE_BC)
	   {
	     /*Loop over all side sets to find a match */
	     for(ins=0; ins<exo->num_side_sets; ins++)
	       {
		 if(Proc_SS_Ids[ins]==BC_Types[ibc].BC_ID)
		   {
		     
		     
		     /* Does it contain the element? */
		     for(side_index=0;side_index<exo->ss_num_sides[ins];side_index++)
		       {
			 
			 if(ielem==exo->ss_elem_list[exo->ss_elem_index[ins]+side_index])
			   {
			     /* which variable is the table for? */
			     switch(BC_Types[ibc].table->f_index)
			       {
			       case POLYMER_STRESS11:
				 {
				   table_ibc[0][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12:
				 {
				   table_ibc[0][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22:
				 {
				   table_ibc[0][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_1:
				 {
				   table_ibc[1][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_1:
				 {
				   table_ibc[1][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_1:
				 {
				   table_ibc[1][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_2:
				 {
				   table_ibc[2][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_2:
				 {
				   table_ibc[2][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_2:
				 {
				   table_ibc[2][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_3:
				 {
				   table_ibc[3][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_3:
				 {
				   table_ibc[3][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_3:
				 {
				   table_ibc[3][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_4:
				 {
				   table_ibc[4][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_4:
				 {
				   table_ibc[4][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_4:
				 {
				   table_ibc[4][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_5:
				 {
				   table_ibc[5][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_5:
				 {
				   table_ibc[5][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_5:
				 {
				   table_ibc[5][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_6:
				 {
				   table_ibc[6][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_6:
				 {
				   table_ibc[6][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_6:
				 {
				   table_ibc[6][1][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS11_7:
				 {
				   table_ibc[7][0][0]=ibc;
				   break;
				 }
			       case POLYMER_STRESS12_7:
				 {
				   table_ibc[7][0][1]=ibc;
				   break;
				 }
			       case POLYMER_STRESS22_7:
				 {
				   table_ibc[7][1][1]=ibc;
				   break;
				 }
			       default:
				 {
				   break;
				 }
			       }
			   } /* end check for right element */
		       } /* end loop over sides of ins */
		   } /*close if loop for BC_ID check */
	       } /*close loop over all side sets */
	   } /* end if loop over tables */
       } /* end loop over all bc's */

     num_zeros = 0;

     for(mode=0; mode<vn->modes; mode++) {
       for(a=0; a<VIM; a++){
	 for(b=0; b<VIM; b++)
	   {
	     if((a <= b)&&(table_ibc[mode][a][b] == -1)) 
	       { num_zeros++; }
	   }
       }
     }

     /* if num_zeros == 0, then table bc's exist for all stress modes
	if num_zeros == vn->modes*(VIM!), then no table bc's exist
	if num_zeros in between - there is a problem...  */


  /* Use one point surface quadrature integration 
     to get the sign of v*n */

/*   ip_total_fill = elem_info(NQUAD_SURF, ielem_type_fill); */

  ip_total_fill =1;

  /* Surface integration over element */
  
  for (ip = 0; ip < ip_total_fill; ip++) 
    {
      /* find the quadrature point locations for current ip */
      
      
      find_surf_st (ip, P1_QUAD, id_side, dim, xi, &ss, &tt, &uu);
      
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
      
      surface_determinant_and_normal (ielem, ei->iconnect_ptr, num_local_nodes, 
				      ei->ielem_dim - 1,  
				      id_side,
				      nodes_per_side,
				      local_elem_node_id );

      do_LSA_mods(LSA_SURFACE);

      vdotn_avg = 0.;
      vdotn_norm = 0.;
      for( a=0; a< dim; a++)
	{
	  vdotn_avg += fv->v[a] * fv->snormal[a];
	  vdotn_norm += fv->v[a]* fv->v[a];
	}

      if (vdotn_avg <  0. && vdotn_avg*vdotn_avg/vdotn_norm > 1.e-12 )
	{
 	  if(neighbor != -1 )
	    {
	      err =  neighbor_stress(exo, x, x_update, ielem, neighbor, stress_neighbor, 
				     stress_update_v, phi_neighbor,
				     num_local_nodes, nodes_per_side,
				     local_elem_node_id, ielem_type, 
				     ielem_type_fill,  x_neighbor, S_map);
	      EH( err, "neighbor_stress");
	    }
 	  else if((neighbor==-1)&&(num_zeros==0))
	    {
	      /* inlet table boundary consitions exist for the stress components */

	      err =  neighbor_stress_table(exo, x, x_update, ielem, stress_neighbor,  
				     stress_update_v, phi_neighbor,
				     num_local_nodes, nodes_per_side,
				     local_elem_node_id, ielem_type, 
				     ielem_type_fill,  x_neighbor, S_map, table_ibc);
	      EH( err, "neighbor_stress_table");
	    }
 	  else
	    {
	      /* if there is no neighbor and no tables, set this value to zero 
	       * I am assuming we will only get here for inflow
	       * boundaries... later we will do something better.*/
	      for (a = 0; a < MAX_SURF_GP; a++) {
		for (b = 0; b < MDE; b++) {
		  phi_neighbor[a][b] = 0.0;
		}
	      }

	      phi_v = phi_neighbor[0];

	      for ( mode=0; mode<vn->modes; mode++)
		{
		  for ( a=0; a<VIM; a++)
		    {
		      for ( b=0; b<VIM; b++)
			{
			  s_n[mode][a][b] = 0.;
			}
		    }
		}
	    }
	}
    }


  /* Surface integration over element */
  if(vdotn_avg <  0. && vdotn_avg*vdotn_avg/vdotn_norm > 1.e-12 && ((neighbor!=-1)||(num_zeros==0)) ) 
    {
      for (ip = 0; ip < ip_total; ip++) 
	{
	  /* find the quadrature point locations for current ip */
	  
	  
	  find_surf_st (ip, ielem_type, id_side, dim, xi, &ss, &tt, &uu);
	  
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
	  
	  surface_determinant_and_normal (ielem, ei->iconnect_ptr, num_local_nodes, 
					  ei->ielem_dim - 1,  
					  id_side,
					  nodes_per_side,
					  local_elem_node_id );

	  do_LSA_mods(LSA_SURFACE);

	  if((neighbor != -1)||(num_zeros==0))
	    {
	      found_it = 0; 
	      for (ip1 = 0; ip1 < ip_total && (!found_it); ip1++) 
		{
		  if(  (fabs(fv->x0[0]-x_neighbor[ip1][0])<1.e-7)
		       &&(fabs(fv->x0[1]-x_neighbor[ip1][1])<1.e-7)
		       &&(fabs(fv->x0[2]-x_neighbor[ip1][2])<1.e-7))
		    {
		      found_it = 1;
		      phi_v = phi_neighbor[ip1];

		      for ( mode=0; mode<vn->modes; mode++)
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
				  /* since the stress tensor is symmetric, only assemble the upper half */ 
				  if(a <= b)
				    {
				      s_n[mode][a][b] = stress_neighbor[ip1][mode][a][b];
				    }
				}
			    }
			}
		    }
		}
	    }

	  vdotn = 0.;
	  for( a=0; a< dim; a++)
	    {
	      vdotn += fv->v[a]* fv->snormal[a];
	    }

	  /* add alpha for upwinding and possible stability */
	  vdotn *= alpha;

	  vdotn_norm = sqrt(vdotn*vdotn);
	  
	  if ((vdotn <  0.) && (vdotn_norm > 1.e-7))
	    {
	      /*
	       * Put local contributions into global right-hand side
	       * if it is not a right-hand side variable-it won't get added in (contribution
	       * is zero)
	       */
	      for ( mode=0; mode<vn->modes; mode++)
		{
		  load_modal_pointers(mode, t_, delta_t, s, s_dot, grad_s, d_grad_s_dmesh);

		  if( vn->dg_J_model == FULL_DG && Linear_Solver != FRONT)
		    {
		      load_neighbor_pointers( exo, ams, neighbor, ielem_type, mode, R_s, v_s, J_S_S_v[mode] );
		    }

		  if ( af->Assemble_Residual )
		    {
		      for ( a=0; a<VIM; a++)
			{
			  for ( b=0; b<VIM; b++)
			    {
			      /* since the stress tensor is symmetric, only assemble the upper half */ 
			      if(a <= b)
				{
				  eqn = R_s[mode][a][b];
				  if ( pd->e[eqn] )
				    {
				      peqn = upd->ep[eqn];
				      
				      for ( i=0; i<ei->dof[eqn]; i++)
					{
					  phi_i = bf[eqn]->phi[i];

					  rhs = phi_i * wt * fv->sdet * 
					    ve[mode]->time_const * 
					    vdotn  * (s[a][b]- s_n[mode][a][b]);

					  lec->R[peqn][i] -= rhs;
					}
				    }
				}
			    }
			}
		    }
	      
		  if ( af->Assemble_Jacobian && doit)
		    {
		      for ( a=0; a<VIM; a++)
			{
			  for ( b=0; b<VIM; b++)
			    {
			      /* since the stress tensor is symmetric, only assemble the upper half */ 
			      if(a <= b)
				{
				  eqn = R_s[mode][a][b];
				  if ( pd->e[eqn] )
				    {
				      peqn = upd->ep[eqn];
				      
				      for ( i=0; i<ei->dof[eqn]; i++)
					{
					  phi_i = bf[eqn]->phi[i];
					  
					  /*
					   * J_S_S
					   */
					  var = eqn;
					  pvar = upd->vp[var];
					  for ( j=0; j<ei->dof[var]; j++)
					    {
					      phi_j = bf[var]->phi[j];
					      
					      advection = wt * fv->sdet * vdotn *
						ve[mode]->time_const *
						phi_j * phi_i;
					      
					      /* Work better without this?????, see PRS concern above 
					         Or is this a correction ???*/
					       lec->J[peqn][pvar][i][j] -= 
						advection; 

					       advection = 0;
					      if ( vn->dg_J_model == FULL_DG && found_it)
						{
						    advection = wt * fv->sdet * vdotn *
						      ve[mode]->time_const *
						      phi_v[j] * phi_i;
						  
						  if (Linear_Solver != FRONT)
						    {
						      *J_S_S_v[mode][a][b][i][a][b][j] += advection;
						    }
						  else
						    {

						      /* Notice how this is diagonal only */
						      /* i.e., T_12_i is only depending on T_12_j on other face element*/
						      lec->J_stress_neighbor[id_side-1][i][peqn][j] += advection;
						    }
						}

					    }

				      
					  /*
					   * J_S_v  sensitivity of stress equation w.r.t. velocity
					   */			  
					  
					  for ( p=0; p<dim; p++)
					    {
					      var = VELOCITY1+p;
					      if ( pd->v[var] )
						{
						  pvar = upd->vp[var];
						  for ( j=0; j<ei->dof[var]; j++)
						    {
						      phi_j = bf[var]->phi[j];
						      
						      advection =  phi_i * wt * fv->sdet * 
							ve[mode]->time_const *
							phi_j * fv->snormal[p] * alpha *
							(s[a][b]- s_n[mode][a][b]);
						      
						      lec->J[peqn][pvar][i][j] -=
							advection;
						      
						    }
						}
					    }
					  
					  /*
					   * J_S_d  sensitivity of stress equation w.r.t. mesh displacement
					   */
					  
					  for ( p=0; p<dim; p++)
					    {
					      var = MESH_DISPLACEMENT1+p;
					      if ( pd->v[var] )
						{
						  pvar = upd->vp[var];
						  for ( j=0; j<ei->dof[var]; j++)
						    {
						      phi_j = bf[var]->phi[j];
						      
						      advection = 0.;
						      for (q=0; q<dim; q++)
							{
							  advection += fv->sdet* fv->v[q] 
							    * fv->dsnormal_dx[q][p][j];
							
							}

						      advection += fv->dsurfdet_dx[p][j] * vdotn;

						      advection *= phi_i * wt * 
							ve[mode]->time_const *
							(s[a][b]- s_n[mode][a][b]);
						      
						      lec->J[peqn][pvar][i][j] -=
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
		} /* end mode loop */

	    }
	}
    }
		      
  for (i=0; i<ip_total; i++) {
    safer_free(( void ** ) &x_neighbor[i] );
    }

  safer_free( (void **) &x_neighbor);

  return 0;
}
/*****************************************************************************/
/* END of routine assemble_surface_stress */
/*****************************************************************************/



int 
neighbor_stress(Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
		dbl x[],
		dbl x_update[],
		int current_elem,
		int neighbor_elem,
		dbl stress_neighbor[][MAX_MODES][DIM][DIM],
		dbl snv[][MAX_MODES][DIM][DIM],
		dbl phi_v[][MDE],
		int num_local_nodes,
		int nodes_per_side,
		int local_elem_node_id[],
		int ielem_type,
		int ielem_type_fill,
		dbl **x_n,
		int v_s[MAX_MODES][DIM][DIM])
/*   
 *   This function take the current element and side and 
 *   knowing who the neighboring element is, finds the value
 *   of the fill function at the current gauss point location
 *   in the neighboring element.
 *

 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 */

{
  int gnn, ledof, i, j, mode;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ie;
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int a,b, p;     /* counters */
  const int dim =  pd->Num_Dim;
  int nvdof;
  int status = 0;
  int v;

  dbl phi[MDE], phi_map[MDE], arg_j, s, t, u;
  dbl xi[DIM];

  ielem_shape     = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);


  for(ip = 0; ip < ip_total; ip++) 
    {
      for(p=0; p<DIM; p++)
	{
	  x_n[ip][p] = 0.;
	}

      for(j=0; j<MDE; j++)
	{
	  phi_v[ip][j] = 0.;
	}

      for ( mode=0; mode<vn->modes; mode++)
	{
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  {
		    stress_neighbor[ip][mode][a][b] = 0.;
		    snv[ip][mode][a][b] = 0.;
		  }
		}
	    }
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
  
  id_side = find_id_side(neighbor_elem, nodes_per_side,
			 inode, id_local_elem_coord, exo);

  for (ip = 0; ip < ip_total; ip++) 
    {

      /* find the quadrature point locations for current ip */
      find_surf_st(ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
      
      /* 
       *  we are cheating here and hoping that the "ei->dof[v]" for
       *  the current element will work on the neighbor element that
       *  we are trying to get information for 
       */
      v = POLYMER_STRESS11;
     
      /* first load phi for the fill function */      
      for (i = 0; i < ei->dof[v]; i++)  
	{
	  phi[i] = newshape(xi, ielem_type, PSI, 
			    ei->dof_list[v][i], 
			    ielem_shape, pd->i[v], i);

	  phi_v[ip][i] = phi[i];
	}

      v = pd->ShapeVar;
      /* 
       *  we are cheating here and hoping that the "ei->dof[v]" for
       *  the current element will work on the neighbor element that
       *  we are trying to get information for 
       */
      
      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */
      for (i = 0; i < ei->dof[v]; i++)  
	{
	  phi_map[i] = newshape(xi, ielem_type, PSI, 
				ei->dof_list[v][i], 
				ielem_shape, pd->i[v], i);
	}

      iconnect_ptr    = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */
      for ( i=0; i< ei->dof[v]; i++)
	{
	  gnn = Proc_Elem_Connect[iconnect_ptr + i];
	  
	  for ( p=0; p<dim; p++)
	    {
	      x_n[ip][p] +=  Coor[p][gnn] * phi_map[i];  
	    }
	}
      
      for (mode = 0; mode < vn->modes; mode++) {
	for (a = 0; a < VIM; a++) {
	  for (b = 0; b < VIM; b++) {
	    if (a <= b) {
	      v = v_s[mode][a][b];
	      if (pd->v[v]) {
		for (i = 0; i < num_local_nodes; i++) {
		  gnn = Proc_Elem_Connect[iconnect_ptr + i];
		  nvdof = Dolphin[gnn][v];
		  for (j = 0; j < nvdof; j++) {
		    ledof = ei->lvdof_to_ledof[v][j];
		    ie = Index_Solution(gnn, v, 0, j,
					ei->matID_ledof[ledof]);
		    EH(ie, "Could not find vbl in sparse matrix.");
		    if ( vn->dg_J_model == EXPLICIT_DG ||
			 vn->dg_J_model == SEGREGATED) {
		      arg_j =  x[ie] -  vn->dg_J_model_wt[0] * x_update[ie];
		    } else {
		      arg_j =  x[ie] ;
		    }
		    stress_neighbor[ip][mode][a][b] += arg_j * phi[j] ;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

  status = 1;
  return(status);
}
/*****************************************************************************/
/* END of routine neighbor_stress */
/*****************************************************************************/


int 
neighbor_stress_table(Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
		dbl x[],
		dbl x_update[],
		int current_elem,
		dbl stress_neighbor[][MAX_MODES][DIM][DIM],
		dbl snv[][MAX_MODES][DIM][DIM],
		dbl phi_v[][MDE],
		int num_local_nodes,
		int nodes_per_side,
		int local_elem_node_id[],
		int ielem_type,
		int ielem_type_fill,
		dbl **x_n,
		int v_s[MAX_MODES][DIM][DIM],
		int table_ibc[MAX_MODES][DIM][DIM])
/*   
 *   This function take the current element and side and 
 *   knowing where table boundary conditions are stored,
 *   finds the value of the fill function at the current
 *   gauss point location in the current element.
 *
 */

{
  int gnn, i, j, mode;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int a,b, p;     /* counters */
  const int dim =  pd->Num_Dim;
  int status = 0;
  int v;

  dbl phi_map[MDE], s, t, u;
  dbl xi[DIM];
  dbl x_ip[DIM], interp_val, slope;



  ielem_shape     = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);


  for(p=0; p<DIM; p++)
    {
      x_ip[p] = 0.;
    }

  for(ip = 0; ip < ip_total; ip++) 
    {
      for(p=0; p<DIM; p++)
	{
	  x_n[ip][p] = 0.;
	}

      for(j=0; j<MDE; j++)
	{
	  phi_v[ip][j] = 0.;
	}

      for ( mode=0; mode<vn->modes; mode++)
	{
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  {
		    stress_neighbor[ip][mode][a][b] = 0.;
		    snv[ip][mode][a][b] = 0.;
		  }
		}
	    }
	}
    }


  /* first get global node numbers for current element and side */
  iconnect_ptr    = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */
  
  for (i = 0; i < nodes_per_side; i++)
    {
      id    = local_elem_node_id[i];
      /* load up global node numbers on this side */
      inode[i] = Proc_Elem_Connect[iconnect_ptr + id];
    }
  
  /* find localside number for neighbor element from
     global node numbers of current element */
  
  id_side = find_id_side(current_elem, nodes_per_side,
			 inode, id_local_elem_coord, exo);


  for (ip = 0; ip < ip_total; ip++) 
    { 

      /* find the quadrature point locations for current ip */
      find_surf_st(ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);
      

      v = pd->ShapeVar;

      
      iconnect_ptr    = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */
      for (i = 0; i < ei->dof[v]; i++)  
	{
	  phi_map[i] = newshape(xi, ielem_type, PSI, 
				ei->dof_list[v][i], 
				ielem_shape, pd->i[v], i);
	}

      for ( i=0; i< ei->dof[v]; i++)
	{
	  gnn = Proc_Elem_Connect[iconnect_ptr + i];

	  for ( p=0; p<dim; p++) 
	    {
	      x_n[ip][p] +=   Coor[p][gnn] * phi_map[i];  
	    } 

     
	}  /* sets global positions for gaussian points */

      for (p=0; p<dim; p++)
	{
	  x_ip[p] = x_n[ip][p];  
	}


      for (mode = 0; mode < vn->modes; mode++) {
	for (a = 0; a < VIM; a++) {
	  for (b = 0; b < VIM; b++) {
	    if (a <= b) {

	      if (table_ibc[mode][a][b] != 0) 
	        {
		  /*       printf("index %d, x_ip %f,%f \n",BC_Types[table_ibc[mode][a][b]].table->t_index[0],x_ip[0],x_ip[1]); */

	  /*set first coordinate of x_ip to interpolation abscissa (assumes only one) */
		  if(BC_Types[table_ibc[mode][a][b]].table->t_index[0]==1)
		    {
		      x_ip[0]=x_ip[1];
		    }
                interp_val = interpolate_table(BC_Types[table_ibc[mode][a][b]].table,x_ip,&slope, NULL);
		stress_neighbor[ip][mode][a][b] = interp_val;
		
		/* printf("mode %d a %d b %d: x %f interp_val %f \n",mode,a,b,x_ip[0],interp_val); */

		} 
	      else 
		{
		stress_neighbor[ip][mode][a][b] = 0.;
		}
	      }
	    }
	  }
	} /* close mode loop*/
      } /* close ip loop */


  status = 1;
  return(status);
}

/*****************************************************************************/
/* END of routine neighbor_stress_table */
/*****************************************************************************/

void
load_neighbor_pointers( Exo_DB *exo,
			struct Aztec_Linear_Solver_System *ams,
		        int ielem, /* neighbor element */
			int etype, /* element type */
			int mode,  /* stress mode */
			int R_s[MAX_MODES][DIM][DIM],
                             /* Equation number for mode ve_mode */
			int v_s[MAX_MODES][DIM][DIM],
			     /* Variable number for mode ve_mode */
			dbl *J_S_S[DIM][DIM][MDE][DIM][DIM][MDE] )

    /**********************************************************************
     *
     * load_neighbor_pointers:
     *
     *  This routine calculates the pointer array, J_S_S, defined
     *  below:
     *
     * J_S_S: Pointer array. This is a array of pointers to 
     *                   d_R_S_ab,i/d_S_v_pq,j  
     *        where R_S_ab,i is the residual equation to the
     *        ith dof of the ab stress component of the current element
     *        and S_v_pq,j is the jth dof of the pq stress 
     *        component in the current upstream neighbor of the 
     *        current element.
     *********************************************************************/
{
  int iconnect_ptr, ldof, v, ln, nunks, nnodes;
  int *enl;
  int dof_list[MAX_VARIABLE_TYPES+MAX_CONC][MDE];
  int gun_list[MAX_VARIABLE_TYPES+MAX_CONC][MDE];
  int dof[MAX_VARIABLE_TYPES+MAX_CONC];
  int gnn, i, j, ie, je, ja,  meqn1, meqn2, mvar1, mvar2, eqn, var;
  int I, row_dofs, blk_row, J, K, blk_col;
  int *rpntr, *bpntr, *bindx, *indx, *ija;
  double *a;
  NODAL_VARS_STRUCT *nv;
  
  (void) memset((void *)dof, 0, sizeof(int)*(MAX_VARIABLE_TYPES+MAX_CONC));
  iconnect_ptr = exo->elem_node_pntr[ielem];
  enl = exo->elem_node_list + iconnect_ptr;
  nnodes = exo->elem_node_pntr[ielem+1] - iconnect_ptr;
  
  /*
   *  Formulate a list of the number of degrees of freedom, dof[varType],
   *  in element, ielem. Restrict the list to vartypes of type, stress
   *  mode.  Also form the following maps:
   *     dof_list[v][ldof] -> local dof to local node map
   *     gun_list[v][ldof] -> local dof to proc unknown index.
   */
  for (v = v_s[mode][0][0]; v <= v_s[mode][2][2] ; v++) {
    if (Num_Var_In_Type[v]) {
      ldof = 0;
      for (ln = 0; ln < nnodes ; ln++) {
	/*
	 * For this local node "ln", what is the global node number, "gnn"?
	 */
	gnn = *(enl + ln);
	      
	/*
	 * For this variable at this local node, how many dofs are needed?
	 * (according to this element) Note: this can be zero...
	 *
	 */
	nv = Nodes[gnn]->Nodal_Vars_Info;
	nunks = get_nv_ndofs_modMF(nv, v);
#ifdef DEBUG_HKM
	if (nunks != node_info(ln, etype, v, gnn)) {
	  fprintf(stderr,"load_neighbor_pointers ERROR P_%d:", ProcID);
	  fprintf(stderr,"old and new nunks differ: %d %d\n",
		  nunks, node_info(ln, etype, v, gnn));
	  EH(-1, "load_neighbor_pointers: nunks problem");
	}
#endif
	dof[v] += nunks;
	EH(nunks, "problem with nun for this var.");
	      
	for (i = 0; i < nunks; i++ ) {
	  dof_list[v][ldof] = ln;
	  gun_list[v][ldof] = Index_Solution(gnn, v, 0, i, -1);
	  ldof++;
	}
      }
    }
  }

  if (strcmp( Matrix_Format,"msr") == 0) {
    ija = ams->bindx;
    a   = ams->val;
    for (meqn1 = 0; meqn1 < VIM; meqn1++) {
      for (meqn2 = 0; meqn2 < VIM; meqn2++)  {
	if (meqn1 <= meqn2)  {
	  eqn = R_s[mode][meqn1][meqn2];
	  if (pd->e[eqn]) {
	    for (i = 0; i < ei->dof[eqn]; i++) {
	      ie = ei->gun_list[eqn][i];
	      EH(ie, "Bad eqn index.");
	      mvar1 = meqn1;
	      mvar2 = meqn2;
	      if (mvar1 <= mvar2)  {
		var = v_s[mode][mvar1][mvar2];
		if (pd->v[var]) {
		  for (j = 0; j < dof[var]; j++ ) {
		    je = gun_list[var][j];
		    EH(je, "Bad var index.");
		    ja = (ie == je) ? ie :
			in_list(je, ija[ie], ija[ie+1], ija);
		    EH(ja, "Could not find vbl in sparse matrix.");
		    J_S_S[meqn1][meqn2][i] [mvar1][mvar2][j] =
			a + ja;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } else if (strcmp( Matrix_Format, "vbr") == 0) {
    rpntr = ams->rpntr;
    bpntr = ams->bpntr;
    bindx = ams->bindx;
    indx  = ams->indx;
    a     = ams->val;
    for (meqn1 = 0; meqn1 < VIM; meqn1++) {
      for (meqn2 = 0; meqn2 < VIM; meqn2++) {
	if (meqn1 <= meqn2) {
	  eqn = R_s[mode][meqn1][meqn2];
	  if (pd->e[eqn]) {
	    for (i = 0; i < ei->dof[eqn]; i++) {
	      ie = ei->gun_list[eqn][i];
	      I  = ei->gnn_list[eqn][i];
	      row_dofs = rpntr[I+1] - rpntr[I];
	      blk_row = Local_Offset[I][eqn];
	      EH(ie, "Bad eqn index.");
	      mvar1 = meqn1;
	      mvar2 = meqn2;
	      if (mvar1 <= mvar2) {
		var = v_s[mode][mvar1][mvar2];
		if (pd->v[var]) {
		  for (j = 0; j < dof[var]; j++ ) {
		    J =  *( enl +  dof_list[var][j] ) ;
		    K = in_list(J, bpntr[I], bpntr[I+1], bindx);
		    blk_col = Local_Offset[J][var];
		    J_S_S[meqn1][meqn2][i] [mvar1][mvar2][j] =
			a + indx[K] + row_dofs * (blk_col + j)
			+ blk_row + i;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } else if (strcmp( Matrix_Format, "epetra") == 0) {
    EH(-1, "load_neighbor_pointers unsupported by epetra");
  }
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
int
segregate_stress_update ( double x_update[] )
{
  int a,b,ieqn,pvar, p, q;
  int eqn, var;
  int ldof_eqn, ldof_var;
  int ln;
  int R_s[DIM][DIM];
  double lump=0;
  int status = 0;


  if( vn->dg_J_model != SEGREGATED ) return (0) ;

  R_s[0][0] = POLYMER_STRESS11;
  R_s[0][1] = POLYMER_STRESS12;
  R_s[0][2] = POLYMER_STRESS13;
  R_s[1][0] = POLYMER_STRESS12;
  R_s[1][1] = POLYMER_STRESS22;
  R_s[1][2] = POLYMER_STRESS23;
  R_s[2][0] = POLYMER_STRESS13;
  R_s[2][1] = POLYMER_STRESS23;
  R_s[2][2] = POLYMER_STRESS33;


  for(ln=0; ln<ei->num_local_nodes; ln++)
    {
      for( a=0; a<DIM; a++)
	{
	  for( b=0; b<DIM; b++)
	    {
	      if( a <= b )
		{
		  eqn = R_s[a][b];
		  
		  if( pd->e[eqn] && ( ldof_eqn = ei->ln_to_first_dof[eqn][ln] != -1 ) )
		    {
		      ieqn = upd->ep[eqn];

		      while ( ldof_eqn <= ei->ln_to_dof[eqn][ln] )
			{
			  for( p=0; p<DIM; p++ )
			    {
			      for( q=0; q<DIM; q++)
				{
				  if( p <= q && delta(a,p) && delta(q,b) )
				    {
				      var = R_s[p][q];

				      lump = 0.0;

				      if( pd->v[var] && (ldof_var = ei->ln_to_first_dof[var][ln] != -1 ) )
					{
					  pvar = upd->vp[var];

					  while ( ldof_var <= ei->ln_to_dof[var][ln] )
					    {
					      lump += lec->J[ieqn][pvar][ldof_eqn][ldof_var];
					      ldof_var++;
					    }
					}
				    }
				}
			    }
			  x_update[ ei->gun_list[eqn][ldof_eqn] ] = lec->R[ieqn][ldof_eqn]/lump;
			  ldof_eqn++;
			}
		    }
		}
	    }
	}
    }
 return(status);
}
		  
/* This routine calculates the adaptive viscosity from Sun et al., 1999.
 * The adaptive viscosity term multiplies the continuous and discontinuous
 * shear-rate, so it should cancel out and not affect the
 * solution, other than increasing the stability of the
 * algorithm in areas of high shear and stress.
 */
dbl
numerical_viscosity(dbl s[DIM][DIM],                       /* total stress */
		    dbl gamma_cont[DIM][DIM],              /* continuous shear rate */
		    dbl d_mun_dS[MAX_MODES][DIM][DIM][MDE],/* derivative of mun wrt S */ 
		    dbl d_mun_dG[DIM][DIM][MDE])           /* derivative of mun wrt G */
{
  int a,b,j;
  int var, mode;
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];
  dbl s_dbl_dot_s;
  dbl g_dbl_dot_g;
  dbl eps2;
  dbl eps; /* should migrate this to input deck */

  dbl mun;

  /* load pointers into equation/variable number */
  (void) stress_eqn_pointer(v_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32; 
  v_g[2][2] = VELOCITY_GRADIENT33; 
  
  eps = vn->eps;
  
  eps2 = eps/2.;

  s_dbl_dot_s = 0.;
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  s_dbl_dot_s += s[a][b]*s[b][a];
	}
    }

  g_dbl_dot_g = 0.;
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  g_dbl_dot_g += gamma_cont[a][b]*gamma_cont[b][a];
	}
    }

  mun = (sqrt(1.+ eps2*s_dbl_dot_s))/sqrt(1.+eps2*g_dbl_dot_g); 


  for ( mode=0; mode<vn->modes; mode++)
    {
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      var = v_s[mode][a][b];
	      
	      for ( j=0; j<ei->dof[var]; j++)
		{
 		  d_mun_dS[mode][a][b][j] = eps2/(sqrt(1.+ eps2*s_dbl_dot_s)*sqrt(1.+eps2*g_dbl_dot_g))
 		    *s[a][b]*bf[var]->phi[j];
		}
	    }
	}
    }

  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  var = v_g[a][b];
	  
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      /* d_mun_dG[a][b][j] = -eps*(sqrt(1.+ eps2*s_dbl_dot_s))/pow((1.+eps2*g_dbl_dot_g),1.5) */
  	      /* 	*gamma_cont[a][b]*bf[var]->phi[j]; */
	      d_mun_dG[a][b][j] = -eps*mun/(1.+eps2*g_dbl_dot_g)

  		*gamma_cont[a][b]*bf[var]->phi[j];
	    }
	}
    }

  
  return(mun);
}
		  
/* This routine sets a handy pointer to give the
 * correct equation number for a given mode
 * and a and b component: v[mode][a][b].
 */
int
stress_eqn_pointer(int v_s[MAX_MODES][DIM][DIM])
{
  int status = 1;
  /* mode 0, polymer stress */
  v_s[0][0][0] = POLYMER_STRESS11;
  v_s[0][0][1] = POLYMER_STRESS12;
  v_s[0][0][2] = POLYMER_STRESS13;
  v_s[0][1][0] = POLYMER_STRESS12;
  v_s[0][1][1] = POLYMER_STRESS22;
  v_s[0][1][2] = POLYMER_STRESS23;
  v_s[0][2][0] = POLYMER_STRESS13;
  v_s[0][2][1] = POLYMER_STRESS23;
  v_s[0][2][2] = POLYMER_STRESS33;

  /* mode 1, polymer stress */
  v_s[1][0][0] = POLYMER_STRESS11_1;
  v_s[1][0][1] = POLYMER_STRESS12_1;
  v_s[1][0][2] = POLYMER_STRESS13_1;
  v_s[1][1][0] = POLYMER_STRESS12_1;
  v_s[1][1][1] = POLYMER_STRESS22_1;
  v_s[1][1][2] = POLYMER_STRESS23_1;
  v_s[1][2][0] = POLYMER_STRESS13_1;
  v_s[1][2][1] = POLYMER_STRESS23_1;
  v_s[1][2][2] = POLYMER_STRESS33_1;

  /* mode 2, polymer stress */
  v_s[2][0][0] = POLYMER_STRESS11_2;
  v_s[2][0][1] = POLYMER_STRESS12_2;
  v_s[2][0][2] = POLYMER_STRESS13_2;
  v_s[2][1][0] = POLYMER_STRESS12_2;
  v_s[2][1][1] = POLYMER_STRESS22_2;
  v_s[2][1][2] = POLYMER_STRESS23_2;
  v_s[2][2][0] = POLYMER_STRESS13_2;
  v_s[2][2][1] = POLYMER_STRESS23_2;
  v_s[2][2][2] = POLYMER_STRESS33_2;

  /* mode 3, polymer stress */
  v_s[3][0][0] = POLYMER_STRESS11_3;
  v_s[3][0][1] = POLYMER_STRESS12_3;
  v_s[3][0][2] = POLYMER_STRESS13_3;
  v_s[3][1][0] = POLYMER_STRESS12_3;
  v_s[3][1][1] = POLYMER_STRESS22_3;
  v_s[3][1][2] = POLYMER_STRESS23_3;
  v_s[3][2][0] = POLYMER_STRESS13_3;
  v_s[3][2][1] = POLYMER_STRESS23_3;
  v_s[3][2][2] = POLYMER_STRESS33_3;

  /* mode 4, polymer stress */
  v_s[4][0][0] = POLYMER_STRESS11_4;
  v_s[4][0][1] = POLYMER_STRESS12_4;
  v_s[4][0][2] = POLYMER_STRESS13_4;
  v_s[4][1][0] = POLYMER_STRESS12_4;
  v_s[4][1][1] = POLYMER_STRESS22_4;
  v_s[4][1][2] = POLYMER_STRESS23_4;
  v_s[4][2][0] = POLYMER_STRESS13_4;
  v_s[4][2][1] = POLYMER_STRESS23_4;
  v_s[4][2][2] = POLYMER_STRESS33_4;

  /* mode 5, polymer stress */
  v_s[5][0][0] = POLYMER_STRESS11_5;
  v_s[5][0][1] = POLYMER_STRESS12_5;
  v_s[5][0][2] = POLYMER_STRESS13_5;
  v_s[5][1][0] = POLYMER_STRESS12_5;
  v_s[5][1][1] = POLYMER_STRESS22_5;
  v_s[5][1][2] = POLYMER_STRESS23_5;
  v_s[5][2][0] = POLYMER_STRESS13_5;
  v_s[5][2][1] = POLYMER_STRESS23_5;
  v_s[5][2][2] = POLYMER_STRESS33_5;

  /* mode 6, polymer stress */
  v_s[6][0][0] = POLYMER_STRESS11_6;
  v_s[6][0][1] = POLYMER_STRESS12_6;
  v_s[6][0][2] = POLYMER_STRESS13_6;
  v_s[6][1][0] = POLYMER_STRESS12_6;
  v_s[6][1][1] = POLYMER_STRESS22_6;
  v_s[6][1][2] = POLYMER_STRESS23_6;
  v_s[6][2][0] = POLYMER_STRESS13_6;
  v_s[6][2][1] = POLYMER_STRESS23_6;
  v_s[6][2][2] = POLYMER_STRESS33_6;

  /* mode 7, polymer stress */
  v_s[7][0][0] = POLYMER_STRESS11_7;
  v_s[7][0][1] = POLYMER_STRESS12_7;
  v_s[7][0][2] = POLYMER_STRESS13_7;
  v_s[7][1][0] = POLYMER_STRESS12_7;
  v_s[7][1][1] = POLYMER_STRESS22_7;
  v_s[7][1][2] = POLYMER_STRESS23_7;
  v_s[7][2][0] = POLYMER_STRESS13_7;
  v_s[7][2][1] = POLYMER_STRESS23_7;
  v_s[7][2][2] = POLYMER_STRESS33_7;

  return(status);
}  
/*****************************************************************************/
/* END OF FILE mm_fill_stress.c */
/*****************************************************************************/


void
compute_exp_s(double s[DIM][DIM],
	      double exp_s[DIM][DIM],
              double eig_values[DIM],
              double R[DIM][DIM])
{

  int N = VIM;
  int LDA = N;
  int i,j,k;
  double tmp;

  int INFO;
  int LWORK = 20;
  double WORK[LWORK];

  double A[VIM*VIM];
  double D[DIM][DIM];
  double EIGEN_MAX = sqrt(sqrt(DBL_MAX));
  double eig_S[DIM];
  memset(A, 0.0, sizeof(double)*VIM*VIM);
  memset(D, 0.0, sizeof(double)*DIM*DIM);
  memset(eig_values, 0.0, sizeof(double)*DIM);
  memset(eig_S, 0.0, sizeof(double)*DIM);
  memset(WORK, 0, sizeof(double)*LWORK);

  // convert to column major
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      A[i*VIM + j] = s[j][i];
    }
  }


  // eig solver
  dsyev_("V", "U", &N, A, &LDA, eig_S, WORK, &LWORK, &INFO, 1, 1);


  // transpose (revert to row major)
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      R[i][j] = A[j*VIM + i];
    }
  }

  // exponentiate diagonal
  for (i = 0; i < VIM; i++) {
	eig_values[i] = MIN(exp(eig_S[i]),EIGEN_MAX); 
        }

  /* matrix multiplication, the slow way */
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      exp_s[i][j] = R[i][j] * eig_values[j];
      }
    }
  
  // multiply by transpose
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      tmp = 0.;
      for (k = 0; k < VIM; k++) {
	tmp += exp_s[i][k] * R[j][k];
        }
      D[i][j] = tmp;
      }
    }

  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
	exp_s[i][j] = D[i][j];
        }
      }
} // End compute_exp_s

void
analytical_exp_s(double s[DIM][DIM],
	      double exp_s[DIM][DIM],
              double eig_values[DIM],
              double R[DIM][DIM])
{

  double I_S, II_S, disc, off_diag, q, p2, r, p, phi;
  int i,j,k;
  double tmp;


  double B[DIM][DIM], D[DIM][DIM],eig_S[DIM],Q1[DIM][DIM],Q2[DIM][DIM];
  memset(D, 0.0, sizeof(double)*DIM*DIM);

  /* Use Eigenvalue algorithm from Wikipedia - https://en.wikipedia.org/wiki/
     Eigenvalue_algorithm#Normal%2C_Hermitian%2C_and_real-symmetric_matrices */

  /*  Make sure stress is symmetric  */
  s[1][0] = s[0][1];  s[2][0] = s[0][2];  s[2][1] = s[1][2];

  if ((VIM==2 || pd->CoordinateSystem == CYLINDRICAL))
    {
      eig_S[2] = s[2][2];
      R[2][2] = 1.0;
      I_S = s[0][0]+s[1][1];
      II_S = s[0][0]*s[1][1] - SQUARE(s[0][1]);
      disc = SQUARE(I_S) - 4*II_S;
      eig_S[0] = 0.5*(I_S + sqrt(disc));
      eig_S[1] = 0.5*(I_S - sqrt(disc));
      if( DOUBLE_NONZERO(disc))	{
         for (j = 0; j < pd->Num_Dim; j++) {
            for (i = 0; i < pd->Num_Dim; i++) {
               R[i][j] = s[i][j] - eig_S[1-j]*delta(i,j);
               }
            }
      /*  Normalize eigenvectors   */
         for (j = 0; j < pd->Num_Dim; j++) {
            tmp = sqrt(SQUARE(R[0][j]) + SQUARE(R[1][j]));
            if( DOUBLE_NONZERO(tmp))	{
               for (i = 0; i < pd->Num_Dim; i++) {
                  R[i][j] /= tmp;
                  }
               }
            }
         }
      else
         {
          for (j = 0; j < pd->Num_Dim; j++) {
             for (k = 0; k < pd->Num_Dim; k++) {
                R[j][k] = delta(j,k);
                }
             }
         }
    }
  else
    {
     off_diag = SQUARE(s[0][1]) + SQUARE(s[0][2]) + SQUARE(s[1][2]);
     I_S = s[0][0] + s[1][1] + s[2][2];
     II_S = s[0][0]*s[1][1]+s[0][0]*s[2][2]+s[1][1]*s[2][2]
              -(SQUARE(s[0][1])+SQUARE(s[0][2])+SQUARE(s[1][2]));
     if( DOUBLE_NONZERO(off_diag))	{
        q = I_S/3.;
        p2 = 2*off_diag;
        for(i = 0 ; i < VIM ; i++)	{p2 += SQUARE(s[i][i]-q);}
        p = sqrt(p2/6.);
        for (i = 0; i < VIM; i++) {
           for (j = 0; j < VIM; j++) {
              B[i][j] = (s[i][j]-q*delta(i,j))/p;
              }
           }
        r = B[0][0]*B[1][1]*B[2][2] + 2.*(B[0][1]*B[1][2]*B[0][2])
              -B[0][0]*SQUARE(B[1][2])-B[1][1]*SQUARE(B[0][2])-B[2][2]*SQUARE(B[0][1]);
        r *= 0.5;
        if ( r <= -1.)
           {phi = M_PIE/3.;}
        else if (r >= 1.)
           {phi =0.;}
        else
           {phi = acos(r)/3.;}
        for (i = 0; i < VIM; i++) {
           eig_S[i] = q + 2*p*cos(phi+i*2*M_PIE/3.);
           }
        }
      else
        {
        for(i = 0 ; i < DIM ; i++)	{eig_S[i] = s[i][i];}
        }

     if( DOUBLE_NONZERO(I_S) || DOUBLE_NONZERO(II_S))	{
        for (i = 0; i < VIM; i++) {
          for (j = 0; j < VIM; j++) {
            for (k = 0; k < VIM; k++) {
               Q1[j][k] = s[j][k] - eig_S[(i+1)%3]*delta(j,k);
               Q2[j][k] = s[j][k] - eig_S[(i+2)%3]*delta(j,k);
               }
            }
          for (j = 0; j < VIM; j++) {
            for (k = 0; k < VIM; k++) {
               R[j][i] = Q1[j][k]*Q2[k][j];
               }
            }
          }
      /*  Normalize eigenvectors   */
        for (j = 0; j < VIM; j++) {
          tmp = sqrt(SQUARE(R[0][j]) + SQUARE(R[1][j]) + SQUARE(R[2][j]));
          if( DOUBLE_NONZERO(tmp))	{
            for (i = 0; i < VIM; i++) {
               R[i][j] /= tmp;
               }
            }
         }
       }
     else
       {
        for (j = 0; j < VIM; j++) {
          for (k = 0; k < VIM; k++) {
             R[j][k] = delta(j,k);
             }
          }
       }
    }


  // exponentiate diagonal
  for (i = 0; i < DIM; i++) {
      eig_values[i] = exp(eig_S[i]);
      }

  /* matrix multiplication, the slow way */
  for (i = 0; i < VIM; i++) {
     for (j = 0; j < VIM; j++) {
         exp_s[i][j] = R[i][j] * eig_values[j];
         }
     }
  
  // multiply by transpose
  for (i = 0; i < VIM; i++) {
     for (j = 0; j < VIM; j++) {
         tmp = 0.;
         for (k = 0; k < VIM; k++) {
             tmp += exp_s[i][k] * R[j][k];
             }
         D[i][j] = tmp;
         }
     }
  for (i = 0; i < VIM; i++) {
     for (j = 0; j < VIM; j++) {
         exp_s[i][j] = D[i][j];
         }
    }
  
} // End analytical_exp_s

void
compute_d_exp_s_ds(dbl s[DIM][DIM],                   //s - stress
		   dbl exp_s[DIM][DIM],
		   dbl d_exp_s_ds[DIM][DIM][DIM][DIM])
{
  double s_p[DIM][DIM], s_n[DIM][DIM];
  double exp_s_p[DIM][DIM], exp_s_n[DIM][DIM];
  double eig_values[DIM];
  double R1[DIM][DIM];
  int i,j,p,q;
  double ds, ds_den, fd = FD_FACTOR;


  memset(d_exp_s_ds, 0, sizeof(double)*DIM*DIM*DIM*DIM);
 
#if 1
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
       s_p[i][j] = s[i][j];
       s_n[i][j] = s[i][j];
       }
    }


  for (i = 0; i < VIM; i++) {
    for (j = i; j < VIM; j++) {

      // perturb s
      ds = fd*s[i][j];
      ds = (fabs(ds) < fd ? fd : ds);
      s_p[i][j] += ds;
      s_n[i][j] -= ds;
      if( i != j) {
        s_p[j][i] = s_p[i][j];
        s_n[j][i] = s_n[i][j];
      }
      ds_den = 0.5/ds;

      // find exp_s at perturbed value
#ifdef ANALEIG_PLEASE
      analytical_exp_s(s_p, exp_s_p, eig_values, R1);
#else
      compute_exp_s(s_p, exp_s_p, eig_values, R1);
#endif
#ifdef ANALEIG_PLEASE
      analytical_exp_s(s_n, exp_s_n, eig_values, R1);
#else
      compute_exp_s(s_n, exp_s_n, eig_values, R1);
#endif

      // approximate derivative
      for (p = 0; p < VIM; p++) {
	for (q = 0; q < VIM; q++) {
	  d_exp_s_ds[p][q][i][j] = ds_den*(exp_s_p[p][q] - exp_s_n[p][q]);
	  if(i != j)d_exp_s_ds[p][q][j][i] = d_exp_s_ds[p][q][i][j];
	}
      }
      s_p[i][j] = s[i][j];
      s_n[i][j] = s[i][j];
      if( i != j) {
        s_p[j][i] = s[j][i];
        s_n[j][i] = s[j][i];
        } 
    }
  }
#else
#ifdef ANALEIG_PLEASE
      analytical_exp_s(s, exp_s, eig_values, R1);
#else
      compute_exp_s(s, exp_s, eig_values, R1);
#endif
  memset(d_I_dS,    0, sizeof(double)*DIM*DIM);
  memset(d_II_dS,    0, sizeof(double)*DIM*DIM);
  memset(d_III_dS,    0, sizeof(double)*DIM*DIM);
  memset(d_eig_ds,    0, sizeof(double)*DIM*DIM*DIM);

  I_S = s[0][0] + s[1][1] + s[2][2];
  II_S = s[0][0]*s[1][1]+s[0][0]*s[2][2]+s[1][1]*s[2][2]
              -(SQUARE(s[0][1])+SQUARE(s[0][2])+SQUARE(s[1][2]));
  III_S = s[0][0]*s[1][1]*s[2][2] + 2.*(s[0][1]*s[1][2]*s[0][2])
              -s[0][0]*SQUARE(s[1][2])-s[1][1]*SQUARE(s[0][2])-s[2][2]*SQUARE(s[0][1]);
  for (i = 0; i < DIM; i++) {
      d_I_dS[i][i] = 1.;
      }
  d_II_S[0][0] = s[1][1] + s[2][2];
  d_II_S[1][1] = s[0][0] + s[2][2];
  d_II_S[2][2] = s[0][0] + s[1][1];
  d_II_S[0][1] = d_II_S[1][0] = -2.*s[0][1];
  d_II_S[0][2] = d_II_S[2][0] = -2.*s[0][2];
  d_II_S[1][2] = d_II_S[2][1] = -2.*s[1][2];
  d_III_S[0][0] = s[1][1]*s[2][2]-SQUARE(s[1][2]);
  d_III_S[1][1] = s[0][0]*s[2][2]-SQUARE(s[0][2]);
  d_III_S[2][2] = s[1][1]*s[0][0]-SQUARE(s[0][1]);
  d_III_S[0][1] = d_III_S[1][0] = 2.*s[1][2]*s[0][2] - 2*s[2][2]*s[0][1];
  d_III_S[0][2] = d_III_S[2][0] = 2.*s[0][1]*s[1][2] - 2*s[1][1]*s[0][2];
  d_III_S[1][2] = d_III_S[2][1] = 2.*s[0][1]*s[1][2] - 2*s[1][1]*s[0][2];

  for (i = 0; i < DIM; i++) {
     eig_S[i] = log(eig_values[i]);
     denom = 3*SQUARE(eig_S[i]) -2*eig_S[i]*I_S + II_S; 
     if(DOUBLE_NONZERO(denom))  {
        for (p = 0; p < DIM; p++) {
           for (q = 0; q < DIM; p++) {
              d_eig_ds[i][p][q] = denom_inv*
                     (SQUARE(eig_S[i])*d_I_S[p][q]-eig_S[i]*d_II_S[p][q] + d_III_S[p][q]);
              }
           }
        }
  /* matrix multiplication, the slow way */
  for (p = 0; i < VIM; i++) {
     for (q = 0; j < VIM; j++) {
        for (i = 0; i < DIM; i++) {
           for (j = 0; j < DIM; j++) {
              d_exp_s_ds[p][q][i][j] = eig_values[p]*d_eig_ds[p][i][j]*R1[q][p];
              }
           }
        for (i = 0; i < VIM; i++) {
           for (j = 0; j < VIM; j++) {
              for (k = 0; k < VIM; k++) {
                 d_exp_s_ds[p][q][i][j] += D[i][k] * R1[j][k];
             }
         }
        }
      }
    }

  
 
   }
#endif

}

dbl
compute_saramito_model_terms(dbl stress[DIM][DIM],
			     dbl yieldStress,
				 dbl yieldExpon,
			     SARAMITO_DEPENDENCE_STRUCT* d_sCoeff)
{
  /* start by computing the norm of the deviatoric stress,  sqrt(J_2), 
   * J_2 = 1/(2*DIM)[
   (stress_{0,0} - stress_{1,1})**2 + 
   (stress_{0,0} - stress_{2,2})**2 + 
   (stress_{2,2} - stress_{1,1})**2
   ] +
   stress_{0,1}**2 + stress_{0,2}**2 + stress_{2,1}**2

   * see the following wikipedia page:
   * https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Invariants_of_the_stress_deviator_tensor
   */
  dbl invDenom = 1./(2.*VIM);

  // square of the deviatoric sress norm

  dbl normOfStressDSqr = pow(stress[0][0] - stress[1][1], 2)*invDenom + pow(stress[0][1], 2);
  if(VIM>2){
    normOfStressDSqr += pow(stress[0][0] - stress[2][2], 2)*invDenom + pow(stress[0][2], 2)
      + pow(stress[1][1] - stress[2][2], 2)*invDenom + pow(stress[1][2], 2);
  }

  const dbl normOfStressD = sqrt(normOfStressDSqr);

  const dbl sc = 1. - yieldStress/normOfStressD;

  dbl sCoeff;
  dbl expYSC=1;
  if(yieldExpon > 0){
	  expYSC = exp(sc/yieldExpon);
	  sCoeff = log( 1. + expYSC )*yieldExpon;
  }
  else{
	  sCoeff = fmax(0, sc);
  }

  // take care of indeterminate behavior for normOfStressD == 0
  if(normOfStressD == 0){ 
	  if(yieldStress <= 0){sCoeff = 1;}
	  else                {sCoeff = 0;} 
	}

  if (d_sCoeff != NULL) {
    // if normStress_d < yieldStress, set sensitivities to zero and return 
    if( (sCoeff) == 0 || yieldStress <= 0){
      for(int i=0; i<VIM; i++){
        for(int j=0; j<VIM; j++){
          d_sCoeff->s[i][j] = 0.;
        }
      }

      d_sCoeff->tau_y = 0.;
  }
  else{
    // otherwise, sensitivities need to be calculated
	d_sCoeff->tau_y = -1./(normOfStressD);
	dbl d_sCoeff_d_normOfStressD = yieldStress/(normOfStressDSqr);
	if(yieldExpon > 0){
		const dbl expYSCDerivativeTerm = expYSC/(1+expYSC);
		d_sCoeff->tau_y          *= expYSCDerivativeTerm;
		d_sCoeff_d_normOfStressD *= expYSCDerivativeTerm;
	}
    // first assign elements values of d(normOfStressDSqr)/d(stress);
    d_sCoeff->s[0][0] = 2.*invDenom*( stress[0][0] - stress[1][1]);
    d_sCoeff->s[0][1] = 2.*stress[0][1];
    d_sCoeff->s[1][0] = d_sCoeff->s[0][1];
    d_sCoeff->s[1][1] = -d_sCoeff->s[0][0];

    if(VIM>2){
      d_sCoeff->s[0][0] += 2.*invDenom*( stress[0][0] - stress[2][2]);
	  d_sCoeff->s[0][1] += 4.*invDenom*stress[0][1];
	  d_sCoeff->s[1][0] += 4.*invDenom*stress[1][0];
      d_sCoeff->s[1][1] += 2.*invDenom*(stress[1][1] - stress[2][2]);
      d_sCoeff->s[0][2] = 2.*stress[0][2];
	  d_sCoeff->s[2][0] = d_sCoeff->s[0][2];
      d_sCoeff->s[1][2] = 2.*stress[1][2];
	  d_sCoeff->s[2][1] = d_sCoeff->s[1][2];
      d_sCoeff->s[2][2] = 2.*invDenom*( 2*stress[2][2] - stress[0][0] - stress[1][1]);
    }

      /* use the chain rule to computute sCoeff sensitivies to stress components
       * d(sCoeff)/d(stress)  =   d(sCoeff)/d(normOfStressD)
       *                        * d(normOfStressD)/d(normOfStressDSqr)
       *                        * d(normOfStressDSqr)/d(stress)
       * 
       *                      =   d(sCoeff)/d(normOfStressD)
       *                        * 0.5/normOfStressD
       *                        * d(normOfStressDSqr)/d(stress) 
       */ 
                          
      // invDenom will be used as d(normOfStressD)/d(stress)
      invDenom = 0.5/normOfStressD*d_sCoeff_d_normOfStressD;

      for(int i=0; i<VIM; i++){
        for(int j=i; j<VIM; j++){
          d_sCoeff->s[i][j] *= invDenom;
          d_sCoeff->s[j][i] = d_sCoeff->s[i][j];
        }
      }
    }
  } // d_sCoeff != NULL
  // for debugging purposes.
  /*
    for(int i=0; i<VIM; i++){
    for(int j=0; j<VIM; j++){
    printf("\nstress[%d][%d] = %E", i, j, stress[i][j]);
    }
    }
    printf("\n");
    for(int i=0; i<VIM; i++){
    for(int j=0; j<VIM; j++){
    printf("\nd(sCoeff)/d(stress[%d][%d] = %E)", i, j, d_sCoeff->s[i][j]);
    }
    }
    printf("\n");
    printf("yield stress = %E", yieldStress);
    printf("\n");
    printf("|stress_d|**2 = %E", normOfStressDSqr);
    printf("\n-------------------------------------------------------------");
  */
  return sCoeff;
}

