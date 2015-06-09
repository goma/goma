/************************************************************************ *
 * Goma - Multiphysics finite element software                             *
 * Sandia National Laboratories                                            *
 *                                                                         *
 * Copyright (c) 2015 Sandia Corporation.                                  *
 *                                                                         *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
 * the U.S. Government retains certain rights in this software.            *
 *                                                                         *
 * This software is distributed under the GNU General Public License.      *
\************************************************************************/
 

// Standard include files

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// GOMA include files
#define _MM_FILL_TERMS_SEGREGATED_C
#include "std.h"
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
#include "rf_solver.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_std_models.h"
#include "mm_std_models_shell.h"
#include "mm_segregated_structs.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_fill_terms.h"
#include "mm_fill_terms_segregated.h"


#include "goma.h"
#include "mm_species.h"
#include "rf_allo.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"


void
load_splitb_esp(int ielem, Exo_DB *exo)
{
  int imtrx;
  int eqn;
  int ie;
  int gnn;
  int dofs;
  int i;
  int d;
  int iNdof;

  /* Loop over all matrices and always fill values */
  /* Assuming that these equations are only active in 1 matrix */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    eqn = R_AUX_MOMENTUM1;
    if (upd->ep[imtrx][eqn] >= 0) {
      for (d = 0; d < VIM; d++) {
        eqn = R_AUX_MOMENTUM1 + d;
        dofs = pg->element_dof_info[imtrx][ielem][eqn].dof;
        for (i = 0; i < dofs; i++) {
          gnn = pg->element_dof_info[imtrx][ielem][eqn].gnn[i];
          iNdof = pg->element_dof_info[imtrx][ielem][eqn].iNdof[i];
          ie = Index_Solution(gnn, eqn, 0, iNdof, -1, imtrx);
          (pg->sbesp).v_star[d][i] = pg->matrices[imtrx].x[ie];
        }
      }
    }

    eqn = R_PRESSURE_POISSON;
    if (upd->ep[imtrx][eqn] >= 0) {
      dofs = pg->element_dof_info[imtrx][ielem][eqn].dof;
      for (i = 0; i < dofs; i++) {
        gnn = pg->element_dof_info[imtrx][ielem][eqn].gnn[i];
        iNdof = pg->element_dof_info[imtrx][ielem][eqn].iNdof[i];
        ie = Index_Solution(gnn, eqn, 0, iNdof, -1, imtrx);
        (pg->sbesp).P_star[i] = pg->matrices[imtrx].x[ie];
      }
    }

    eqn = R_MOMENTUM1;
    if (upd->ep[imtrx][eqn] >= 0) {
      for (d = 0; d < VIM; d++) {
        eqn = R_MOMENTUM1 + d;
        dofs = pg->element_dof_info[imtrx][ielem][eqn].dof;
        for (i = 0; i < dofs; i++) {
          gnn = pg->element_dof_info[imtrx][ielem][eqn].gnn[i];
          iNdof = pg->element_dof_info[imtrx][ielem][eqn].iNdof[i];
          ie = Index_Solution(gnn, eqn, 0, iNdof, -1, imtrx);
          (pg->sbesp).v_old[d][i] = pg->matrices[imtrx].x_old[ie];
        }
      }
    }

    eqn = R_PRESSURE;
    if (upd->ep[imtrx][eqn] >= 0) {
      dofs = pg->element_dof_info[imtrx][ielem][eqn].dof;
      for (i = 0; i < dofs; i++) {
        gnn = pg->element_dof_info[imtrx][ielem][eqn].gnn[i];
        iNdof = pg->element_dof_info[imtrx][ielem][eqn].iNdof[i];
        ie = Index_Solution(gnn, eqn, 0, iNdof, -1, imtrx);
        (pg->sbesp).P_old[i] = pg->matrices[imtrx].x_old[ie];
      }
    }
  }
}

void
load_splitb_fv(int ielem)
{
  int v;
  int d;
  int i;
  int imtrx;
  int *pdv;
  int var;

  /* Similar to esp, load values from all matrices */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    pdv = pd->v[imtrx];

    var = AUX_VELOCITY1;
    if (pdv[var]) {
      for (d = 0; d < VIM; d++) {
        v = AUX_VELOCITY1 + d;
        (pg->sbcfv).v_star[d] = 0;
        for (i = 0; i < pg->element_dof_info[imtrx][ielem][v].dof; i++) {
          (pg->sbcfv).v_star[d] += (pg->sbesp).v_star[d][i] * bf[var]->phi[i];
        }
      }
    }

    var = AUX_PRESSURE;
    if (pdv[AUX_PRESSURE]) {
      v = AUX_PRESSURE;
      (pg->sbcfv).P_star = 0;
      for (i = 0; i < pg->element_dof_info[imtrx][ielem][v].dof; i++) {
        (pg->sbcfv).P_star = (pg->sbesp).P_star[i] * bf[var]->phi[i];
      }
    }

    var = VELOCITY1;

    if (pdv[VELOCITY1]) {
      for (d = 0; d < VIM; d++) {
        v = VELOCITY1 + d;
        (pg->sbcfv).v_old[d] = 0;
        for (i = 0; i < pg->element_dof_info[imtrx][ielem][v].dof; i++) {
          (pg->sbcfv).v_old[d] += (pg->sbesp).v_old[d][i] * bf[var]->phi[i];
        }
      }
    }

    var = PRESSURE;
    if (pdv[PRESSURE]) {
      v = PRESSURE;
      (pg->sbcfv).P_old = 0;
      for (i = 0; i < pg->element_dof_info[imtrx][ielem][v].dof; i++) {
        (pg->sbcfv).P_old = (pg->sbesp).P_old[i] * bf[var]->phi[i];
      }
    }
  }
}

void load_splitb_fv_grads(int ielem)
{
  int p;
  int d;
  int i;
  int r;
  int imtrx;
  int *pdv;
  int var;

  /* Similar to esp, load values from all matrices */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    pdv = pd->v[imtrx];

    var = VELOCITY1;
    if (pdv[var]) {
      for (p = 0; p < VIM; p++) {
        var = VELOCITY1 + p;
        for (d = 0; d < VIM; d++) {
          (pg->sbcfv).grad_v_old[p][d] = 0;
          for (r = 0; r < VIM; r++) {
            for (i = 0; i < pg->element_dof_info[imtrx][ielem][var].dof; i++) {
              (pg->sbcfv).grad_v_old[p][d] += (pg->sbesp).v_old[r][i] *  bf[var]->grad_phi_e[i][r][p][d];
            }
          }
        }
      }
    }

    var = AUX_PRESSURE;
    if (pdv[var]) {
      var = AUX_PRESSURE;
      for (d = 0; d < VIM; d++) {
        (pg->sbcfv).grad_P_star[d] = 0;
        for (i = 0; i < pg->element_dof_info[imtrx][ielem][var].dof; i++) {
          (pg->sbcfv).grad_P_star[d] = (pg->sbesp).P_star[i]
              *  bf[var]->grad_phi[i][d];
        }
      }
    }

    var = AUX_VELOCITY1;
    if (pdv[var]) {
      for (p = 0; p < VIM; p++) {
        var = AUX_VELOCITY1 + p;
        for (d = 0; d < VIM; d++) {
          (pg->sbcfv).grad_v_star[p][d] = 0;
          for (r = 0; r < VIM; r++) {
            for (i = 0; i < pg->element_dof_info[imtrx][ielem][var].dof; i++) {
              (pg->sbcfv).grad_v_star[p][d] += (pg->sbesp).v_star[r][i] *  bf[var]->grad_phi_e[i][r][p][d];
            }
          }
        }
      }
    }

    /* div(v_star) */

    if (pdv[AUX_VELOCITY1]) {
      (pg->sbcfv).div_v_star = 0;
      for (d = 0; d < VIM; d++) {
        (pg->sbcfv).div_v_star += (pg->sbcfv).grad_v_star[d][d];
      }
    }
  }
}

void set_bf(const double xi[], BASIS_FUNCTIONS_STRUCT *bf_ptr,
    int type, int dof)
{
  int ln;

  const double s     = xi[0];
  const double t     = xi[1];
  const double u     = xi[2];

/* Assume 2D */
  switch (type) {
  case I_Q1:
    for (ln = 0; ln < dof; ln++) {
      bf_ptr->phi[ln] = shape(s, t, u, BILINEAR_QUAD, PSI, ln);
      bf_ptr->dphidxi[ln][0] = shape(s, t, u, BILINEAR_QUAD, DPSI_S, ln);
      bf_ptr->dphidxi[ln][1] = shape(s, t, u, BILINEAR_QUAD, DPSI_T, ln);
    }
    break;
  case I_Q2:
    for (ln = 0; ln < dof; ln++) {

      bf_ptr->phi[ln] = shape(s, t, u, BIQUAD_QUAD, PSI, ln);
      bf_ptr->dphidxi[ln][0] = shape(s, t, u, BIQUAD_QUAD, DPSI_S, ln);
      bf_ptr->dphidxi[ln][1] = shape(s, t, u, BIQUAD_QUAD, DPSI_T, ln);
    }
    break;
  }
}

int
load_segregated_basis_functions(const double xi[],             /*  [DIM]               */
                     struct Basis_Functions **bfa) /* ptr to basis function *
                                                    * array of interest     */

     /************************************************************************
      *
      * load_basis_functions():
      *
      *    Calculates the values of all the basis functions active in the
      * current element. It also calculates the derivatives of the basis
      * functions wrt local element coordinates.
      *
      ************************************************************************/
{
  int b, i;
  int shape;

  for (b = 0; b < Num_Basis_Functions; b++) {
    shape = QUADRILATERAL;

    if (pd->i[0][AUX_VELOCITY1] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[AUX_VELOCITY1] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }

    if (pd->i[0][AUX_VELOCITY2] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[AUX_VELOCITY2] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }
    if (pd->i[0][AUX_VELOCITY3] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[AUX_VELOCITY3] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }

    if (pd->i[1][AUX_PRESSURE] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[AUX_PRESSURE] = bfd[b];
      set_bf(xi, bfd[b], I_Q1, 4);
    }

    if (pd->i[2][VELOCITY1] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[VELOCITY1] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }

    if (pd->i[2][VELOCITY2] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[VELOCITY1] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }

    if (pd->i[2][VELOCITY3] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[VELOCITY1] = bfd[b];
      set_bf(xi, bfd[b], I_Q2, 9);
    }

    if (pd->i[3][PRESSURE] == bfd[b]->interpolation
        && shape == bfd[b]->element_shape) {
      bf[PRESSURE] = bfd[b];
      set_bf(xi, bfd[b], I_Q1, 4);
    }
  }
  return (0);
} /* END of routine load_basis_functions */

/* 
 *This function assembles the first step of the CBS, split-B, quasi-implicit
 * method.  Here, we solve for an auxiliary velocity from some form of the 
 * momentum equation.  Due to the time discretization, the only non-linear term
 * is mu_star*grad_v_star if mu_star depends on v_star.
 */
int
assemble_aux_u(dbl time,   // Current time
	       dbl dt)     // Current time step
{
  int dim, wim, eqn, peqn, var, pvar, ledof, status=0;
  int a, b, i, ii, j, p, q;
  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];

  // Density and viscosity
  dbl rho, mu_star, mu_old;

  // Relevant field variable quantities
  dbl *v_star, *v_old, *grad_v_star[DIM], *grad_v_old[DIM], P_old;
  dbl gamma_star[DIM][DIM], gamma_old[DIM][DIM];

  // Body force
  dbl f[DIM];
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct;
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  // Residual contributions
  dbl advection, diffusion, diffusion_a, mass, source;
  
  // Residual and Jacobian
  dbl *R, *J=NULL;
 
  // Integration factors
  dbl h3, wt, det_J;
  
  // Basis functions
  struct Basis_Functions *bfm;
  dbl phi_i, *phi_i_vector, (*grad_phi_i_e_a)[DIM]=NULL;
  dbl phi_j, *phi_j_vector;
  
  // Equation term multipliers
  dbl advection_etm, diffusion_etm, mass_etm, source_etm;
  
  
  // Set equation index
  //eqn = R_MOMENTUM1;
  eqn = R_AUX_MOMENTUM1;
  
  // Leave if we are unwanted
  if(!pd->e[pg->imtrx][eqn])
    {
      return(status);
    }

  // Fill field variables and parameters
  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN)
    {
      wim = wim+1;
    }

  wt = fv->wt;
  det_J = bf[eqn]->detJ;	       
  h3 = fv->h3;
  
  //v_star = fv->v_star;
  v_star = fv->v_star;
  v_old = (pg->sbcfv).v_old;
  P_old = (pg->sbcfv).P_old;

  for(a=0; a<VIM; a++)
    {
      grad_v_star[a] = fv->grad_v_star[a];
      //grad_v_star[a] = fv->grad_v[a];
      grad_v_old[a] = (pg->sbcfv).grad_v_old[a];
    }

  for(a=0; a<VIM; a++)
    {
      for(b=0; b<VIM; b++)
	{
	  gamma_star[a][b] = grad_v_star[a][b] + grad_v_star[b][a];
	  gamma_old[a][b] = grad_v_old[a][b] + grad_v_old[b][a];
	}
    }

  rho = density(NULL, time);
  mu_old = viscosity(gn, gamma_old, NULL);
  mu_star = viscosity(gn, gamma_star, NULL);
  (void) momentum_source_term(f, df, time);


  // Start assembling that residual
  if(af->Assemble_Residual)
    {
      for(a=0; a<wim; a++)
	{
	  eqn = R_AUX_MOMENTUM1 + a;
	  //eqn  = R_MOMENTUM1 + a;
	  bfm  = bf[eqn];
	  phi_i_vector = bfm->phi;

	  peqn = upd->ep[pg->imtrx][eqn];
	  R = lec->R[peqn];

	  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	  advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
	  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
	  source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

	  for(i=0; i<ei->dof[eqn]; i++) 
	    {	    
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      if(ei->active_interp_ledof[ledof]) 
	      	{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];
		  
		  phi_i = phi_i_vector[i];
		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];
		  
		  /*
		   * Time derivative contribution, typically fv_dot is used
		   * However, our time discretization does not fit that of the
		   * backward Euler or Crank-Nicolson schemes, so we need to do 
		   * this separately
		   */
		  mass = 0.0;
		  if(pde[eqn] & T_MASS)
		    {
		      mass += rho/dt*(v_star[a]-v_old[a]);
		      mass *= phi_i*h3*wt*det_J;
		      mass *= mass_etm;
		    }
		  
		  // Advection contribution, uses old velocity only
		  advection = 0.0;
		  if(pde[eqn] & T_ADVECTION)
		    {
		      for(p=0; p<wim; p++)
			{
			  advection += v_old[p]*grad_v_old[p][a];
			}
		      advection *= rho*phi_i*h3*wt*det_J;
		      advection *= advection_etm;
		    }
		  
		  // Diffusion contribution, uses both old and aux velocity
		  diffusion = 0.0;
		  if(pde[eqn] & T_DIFFUSION)
		    {
		      for(p=0; p<VIM; p++) 
			{
			  for(q=0; q<VIM; q++) 
			    {
			      diffusion += mu_star/2.0*gamma_star[q][p]*grad_phi_i_e_a[p][q];
			      diffusion += mu_old/2.0*gamma_old[q][p]*grad_phi_i_e_a[p][q];
			      
			      if(p==q)
				{
				  diffusion -= P_old*grad_phi_i_e_a[p][q];
				}
			    }
			}		      
		      diffusion *= h3*wt*det_J;
		      diffusion *= diffusion_etm;
		    }		 
		  
		  // Source term contribution
		  source = 0.0;
		  if(pde[eqn] & T_SOURCE)
		    {
		      source += f[a];
		      source *= phi_i*h3*wt*det_J;
		      source *= source_etm;
		    }
		  
		  R[ii] += mass + advection + diffusion + source;
		} // if active_interp_ledof
	    }     // for i<dof[eqn]
	}         // for a<wim
    }             // if residual
  

  if(af->Assemble_Jacobian)
    {
      for(a=0; a<wim; a++)
	{
	  eqn = R_AUX_MOMENTUM1 + a;
	  //eqn  = R_MOMENTUM1 + a;
	  peqn = upd->ep[pg->imtrx][eqn];
	  bfm  = bf[eqn];	 
	  
	  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
	  source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
	  
	  phi_i_vector = bfm->phi;
	  
	  for(i=0; i<ei->dof[eqn]; i++) 
	    {		  
	      ii = ei->lvdof_to_row_lvdof[eqn][i];	      
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      
	      if(ei->active_interp_ledof[ledof]) 
		{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];		  
		  phi_i = phi_i_vector[i];		  
		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];	
		  
		  // J_v_star
		  for(b=0; b<wim; b++)
		    {
		      var = AUX_VELOCITY1 + b;
		      //var = VELOCITY1 + b;
		      if(pdv[var])
			{
			  pvar = upd->vp[pg->imtrx][var];			  
			  J = lec->J[peqn][pvar][ii];			  
			  phi_j_vector = bf[var]->phi;
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {			      
			      phi_j = phi_j_vector[j];			      
			      
			      mass = 0.0;			     				
			      if((pde[eqn] & T_MASS) && (a == b))
				{
				  mass += rho/dt*phi_j;
				  mass *= phi_i*h3*wt*det_J;
				  mass *= mass_etm;				
				}
			      
			      diffusion = 0.0;
			      if(pde[eqn] & T_DIFFUSION)
				{
				  for(p=0; p<VIM; p++)
				    {
				      for(q=0; q<VIM; q++)
					{
					  diffusion_a = 0.0;
					  diffusion_a += bf[AUX_VELOCITY1+p]->grad_phi_e[j][b][q][p];
					  diffusion_a += bf[AUX_VELOCITY1+q]->grad_phi_e[j][b][p][q];
					  //diffusion_a += bf[VELOCITY1+p]->grad_phi_e[j][b][q][p];
					  //diffusion_a += bf[VELOCITY1+q]->grad_phi_e[j][b][p][q];

					  diffusion += mu_star/2.0*diffusion_a*grad_phi_i_e_a[p][q];
					  //diffusion += d_mu->v[b][j]*gamma_star[p][q]*grad_phi_i_e_a[p][q];
					}
				    }
				  diffusion *= h3*wt*det_J;
				  diffusion *= diffusion_etm;
				}
			      
			      source = 0.0;
			      if(pde[eqn] & T_SOURCE)
				{
				  source += df->v[a][b][j];
				  source *= phi_i*h3*wt*det_J;
				  source *= source_etm;
				}			      
			      J[j] +=  mass + diffusion + source;
			    } 
			}     
		    }	      	  

		} // if active_interp_ledof
	    }     // for i<dof[eqn]
	}         // for a<wim
    }             // if jacobian
  
  return(status);
} // assemble_aux_u
  

/*
 * This function assembles the pressure-poisson equation for the CBS, split-B, quasi-implicit method
 * Here, we solve for an intermediate/correction pressure P_star
 * It is assumed grad(P_star) = 0 on the boundary which allows for the boundary conditions
 * from the auxiliary velocity step to translate directly through the algorithm
 */
int
assemble_press_poisson(dbl time,  // Current time
		       dbl dt)    // Current time step
{
  // Some indices
  int a, i, j, dim, wim, eqn, peqn, var, pvar, status=0;
  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];

  //Density
  dbl rho;
  
  // Relevant field variable quantities
  dbl div_v_star, *grad_P_star;

  // Residual contributions
  dbl diffusion, mass;

  // Residual and Jacobian
  dbl *R, *J;
 
  // Integration factors
  dbl h3, wt, det_J;
  
  // Basis functions
  dbl phi_i, (*grad_phi)[DIM];
  
  // Equation term multipliers
  dbl diffusion_etm, mass_etm;
  
  // Set equation index
  eqn = R_PRESSURE_POISSON;
  peqn = upd->ep[pg->imtrx][eqn];

  // Leave if we are unwanted
  if(!pd->e[pg->imtrx][eqn])
    {
      return(status);
    }

  // Fill field variables and parameters
  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN)
    {
      wim = wim+1;
    }

  R = lec->R[peqn];
  wt = fv->wt;
  det_J = bf[eqn]->detJ;	       
  h3 = fv->h3;
  grad_phi = bf[eqn]->grad_phi;

  //div_v_star = fv->div_v_star;
  //grad_P_star = fv->grad_P_star;
  div_v_star = (pg->sbcfv).div_v_star;
  grad_P_star = fv->grad_P_star;

  rho = density(NULL, time);

  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

  // Residual
  if(af->Assemble_Residual)
    {
      for(i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  
	  // Incompressibility correction
	  mass = 0.0;
	  if(pde[eqn] & T_MASS)
	    {
	      mass += rho*div_v_star;
	      mass *= phi_i*h3*wt*det_J;
	      mass *= mass_etm;
	    }

	  // Pressure laplacian contribution
	  diffusion = 0.0;
	  if(pde[eqn] & T_DIFFUSION)
	    {
	      for(a=0; a<wim; a++)
		{
		  diffusion += grad_P_star[a]*grad_phi[i][a];
		}
	      diffusion *= h3*wt*det_J;
	      diffusion *= diffusion_etm;
	    }

	  R[i] += mass + diffusion*dt;
	  
	} // for i<dof[eqn]
    }     // if residual

  // Jacobian
  if(af->Assemble_Jacobian)
    {
      for(i=0; i<ei->dof[eqn]; i++)
	{ 
	  phi_i = bf[eqn]->phi[i];

	  // J_P_star
	  //var = AUX_PRESSURE
	  var = AUX_PRESSURE;
	  if(pdv[var])
	    {
	      pvar = upd->vp[pg->imtrx][var];	      
	      J = lec->J[peqn][pvar][i];	      
	      for(j=0; j<ei->dof[var]; j++)
		{		  
		  diffusion = 0.0;
		  if(pde[eqn] & T_DIFFUSION)
		    { 
		      for(a=0; a<wim; a++)
			{
			  diffusion += grad_phi[i][a]*bf[var]->grad_phi[j][a];
			}		      
		      diffusion *= h3*wt*det_J;
		      diffusion *= diffusion_etm;
		    }		  		      		    
		  J[j] += diffusion*dt;
		}
	    }	  	
  
	} // for i<dof[eqn] 
    }     // if jacobian

  return(status);
} // assemble_poisson_press



/*
 * This function assembles the pressure projection step of the CBS, split-B, quasi-implicit method
 * Here, we obtain the velocity from the pressure projection of the auxiliary velocity
 */

int assemble_press_proj(dbl time,  // Current time
			dbl dt)    // Current time step
{
  // Some indices
  int dim, wim, eqn, peqn, var, pvar, ledof, status=0;
  int a, b, i, j, ii;
  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];
  

  // Density
  dbl rho;

  // Relevant field variable quantities
  dbl *v_star, *v, *grad_P_star;

  // Residual contributions
  dbl diffusion, mass;
  
  // Residual and Jacobian
  dbl *R, *J;
 
  // Integration factors
  dbl h3, wt, det_J;
  
  // Basis functions
  struct Basis_Functions *bfm;
  dbl phi_i, *phi_i_vector;
  dbl phi_j, *phi_j_vector;
  
  // Equation term multipliers
  dbl diffusion_etm, mass_etm;
  
  
  // Set equation index
  //eqn = R_PRESSURE_PROJECTION1
  eqn = R_MOMENTUM1;

  // Leave if we are unwanted
  if(!pd->e[pg->imtrx][eqn])
    {
      return(status);
    }

  // Fill field variables and parameters
  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN)
    {
      wim = wim+1;
    }

  wt = fv->wt;
  det_J = bf[eqn]->detJ;	       
  h3 = fv->h3;

  //v_star = fv->v_star;
  v_star = (pg->sbcfv).v_star;
  v = fv->v;
  //grad_P_star = fv->grad_P_star
  grad_P_star = (pg->sbcfv).grad_P_star;

  rho = density(NULL, time);


  // Start assembling that residual
  if(af->Assemble_Residual)
    {
      for(a=0; a<wim; a++)
	{
	  //eqn = R_PRESSURE_PROJECTION1 + a;
	  eqn  = R_MOMENTUM1 + a;
	  peqn = upd->ep[pg->imtrx][eqn];
	  bfm  = bf[eqn];

	  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
	  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

	  R = lec->R[peqn];
	  phi_i_vector = bfm->phi;

	  for(i=0; i<ei->dof[eqn]; i++) 
	    {	    
	       ledof = ei->lvdof_to_ledof[eqn][i];
	      if(ei->active_interp_ledof[ledof]) 
	      	{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];		  
		  phi_i = phi_i_vector[i];
		  
		  // Velocity contribution
		  mass = 0.0;
		  if(pde[eqn] & T_MASS)
		    {
		      mass += v[a] - v_star[a];
		      mass *= phi_i*h3*wt*det_J;
		      mass *= mass_etm;
		    }	      
		  
		  // Pressure contribution
		  diffusion = 0.0;
		  if(pde[eqn] & T_DIFFUSION)
		    {		    
		      diffusion += dt/rho*grad_P_star[a];
		      diffusion *= phi_i*h3*wt*det_J;
		      diffusion *= diffusion_etm;
		    }
		  
		  R[ii] += mass + diffusion;		  
		} // if active_interp_ledof
	    }     // for i<dof[eqn]
	}         // for a<wim
    }             // if residual

  // Jacobian
  if(af->Assemble_Jacobian)
    {      
      for(a=0; a<wim; a++)
	{
	  eqn  = R_MOMENTUM1 + a;
	  peqn = upd->ep[pg->imtrx][eqn];
	  bfm  = bf[eqn];	 
	  
	  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
	  
	  phi_i_vector = bfm->phi;
	  
	  for(i=0; i<ei->dof[eqn]; i++) 
	    {		  
	      ii = ei->lvdof_to_row_lvdof[eqn][i];	      
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      
	      if(ei->active_interp_ledof[ledof]) 
		{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];		  
		  phi_i = phi_i_vector[i];	
		  
		  // J_v
		  for(b=0; b<wim; b++)
		    {
		      var = VELOCITY1 + b;
		      if(pdv[var])
			{
			  pvar = upd->vp[pg->imtrx][var];			  
			  J = lec->J[peqn][pvar][ii];			  
			  phi_j_vector = bf[var]->phi;
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {			      
			      phi_j = phi_j_vector[j];
			      
			      mass = 0.0;			  			    
			      if((pde[eqn] & T_MASS) && (a == b))
				{
				  mass += phi_i*phi_j;
				  mass *= h3*wt*det_J;
				  mass *= mass_etm;
				}
			      J[j] += mass;
			    }
			}
		    }
		  
		} // if active_interp_ledof
	    }     // for i<dof[eqn]
	}         // for a<wim
    }             // if jacobian
  
  return(status);
} // assemble_press_proj




/*
 * This function assembles the pressure update in the CBS, split-B, quasi-implicit method
 * Here, we find the new pressure from the old pressure, auxiliary pressure, and
 * incompressibility
 */
int
assemble_press_update(void)
{
  // Some indices
  int a, b, i, j, dim, wim, eqn, peqn, var, pvar, status=0;
  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];

  // Viscosity
  dbl mu_star;
  
  // Relevant field variable quantities
  dbl div_v_star, P, P_star, P_old;
  dbl *grad_v_star[DIM], gamma_star[DIM][DIM];

  // Residual contributions
  dbl diffusion, mass;

  // Residual
  dbl *R, *J;
 
  // Integration factors
  dbl h3, wt, det_J;
  
  // Basis functions
  dbl phi_i, phi_j;
  
  // Equation term multipliers
  dbl diffusion_etm, mass_etm;
  
  // Set equation index
  //eqn = R_PRESSURE_UPDATE;
  eqn = R_PRESSURE;
  peqn = upd->ep[pg->imtrx][eqn];

  // Leave if we are unwanted
  if(!pd->e[pg->imtrx][eqn])
    {
      return(status);
    }

  // Fill field variables and parameters
  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN)
    {
      wim = wim+1;
    }

  R = lec->R[peqn];
  wt = fv->wt;
  det_J = bf[eqn]->detJ;	       
  h3 = fv->h3;

  //div_v_star = fv->div_v_star;
  //grad_v_star = fv->grad_v_star;
  div_v_star = (pg->sbcfv).div_v_star;
  P = fv->P;
  //P_star = fv->P_star;
  P_star = (pg->sbcfv).P_star;
  P_old = fv_old->P;

  for(a=0; a<VIM; a++)
    {
      //grad_v_star[a] = fv->grad_v_star[a];
      grad_v_star[a] = (pg->sbcfv).grad_v_star[a];
    }
  
  for(a=0; a<VIM; a++)
    {
      for(b=0; b<VIM; b++)
	{
	  gamma_star[a][b] = grad_v_star[a][b] + grad_v_star[b][a];
	}
    }

  mu_star = viscosity(gn, gamma_star, NULL);

  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

  // Residual
  if(af->Assemble_Residual)
    {
      for(i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];
	  
	  // Pressure contribution
	  mass = 0.0;
	  if(pde[eqn] & T_SOURCE)
	    {
	      mass += P - P_star - P_old;
	      mass *= phi_i*h3*wt*det_J;
	      mass *= mass_etm;
	    }

	  // Incompressibility contribution
	  diffusion = 0.0;
	  if(pde[eqn] & T_ADVECTION)
	    {
	      diffusion += mu_star/2.0*div_v_star;
	      diffusion *= h3*wt*det_J;
	      diffusion *= diffusion_etm;
	    }

	  R[i] += mass + diffusion;	  
	} // for i<dof[eqn]
    }     // if residual
  

  // Jacobian
  if(af->Assemble_Jacobian)
    {
      for(i=0; i<ei->dof[eqn]; i++)
	{ 
	  phi_i = bf[eqn]->phi[i];

	  // J_P
	  var = PRESSURE;
	  if(pdv[var])
	    {
	      pvar = upd->vp[pg->imtrx][var];	      
	      J = lec->J[peqn][pvar][i];	      
	      for(j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		
		  mass = 0.0;
		  if(pde[eqn] & T_SOURCE)
		    { 
		      mass += phi_i*phi_j;		      
		      mass *= h3*wt*det_J;
		      mass *= mass_etm;
		    }		  		      		    
		  J[j] += mass;
		}
	    }	  	
  
	} // for i<dof[eqn] 
    }     // if jacobian

  return(status);
} // assemble_press_update
