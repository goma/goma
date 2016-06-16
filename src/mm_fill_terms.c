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
 *$Id: mm_fill_terms.c,v 5.82 2010-07-28 14:46:55 sarober Exp $
 */

/*
 * Copyright (C) 1993, 1994  Sandia National Laboratories
 */

/* Standard include files */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* GOMA include files */
#define _MM_FILL_TERMS_C
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


#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_fill_terms.h"


#include "goma.h"
#include "mm_species.h"
#include "rf_allo.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"

// direct call to a fortran LAPACK eigenvalue routine
extern FSUB_TYPE dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA,
			double *W, double *WORK, int *LWORK, int *INFO, 
			int len_jobz, int len_uplo); 

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
*  assemble_mesh		int 	
*  assemble_energy		int 	
*  assemble_momentum		int
*  assemble_continuity          int
*  assemble_volume              int 
*  assemble_curvature
*  assemble_div_normals
*  assemble_LSvelocity          int 
*
*  assemble_normals
*  assemble_ls_momentum_source
*  apply_ls_momentum_source 
*  scalar_fv_fill
*  grad_scalar_fv_fill
*  load_fv                      int
*  load_fv_grads                int
*  load_fv_mesh_derivs          int
*
*  density                      double
*  conductivity
*  heat_capacity
*  ls_modulate_thermalconductivity 
*  ls_modulate_heatcapacity
*  acoustic_impedance
*  wave_number
*  acoustic_absorption
*  light_absorption
*  refractive_index
*  momentum_source_term         int
*
*  ls_modulate_momentumsource
*  apply_table_mp
*  interpolate_table
*  continuous_surface_tension
*  quad_isomap_invert
*  scalar_fv_fill_altmatrl
*  scalar_fv_fill_adjmatrl
*  load_matrl_statevector
*  FoamVolumeSource
*  REFVolumeSource
*  assemble_projection_stabilization
*  assemble_PPPS_generalized
*  apply_distributed_sources 
*  assemble_pf_capillary
*  assemble_csf_tensor
*  assemble_div_n_source
*  assemble_div_s_n_source 
*  assemble_curvature_with_normals_source
*  assemble_curvature_source
*  assemble_q_source
*  assemble_t_source  
*  assemble_qlaser_source
*  assemble_qvapor_source 
*  assemble_qrad_source 
*  assemble_ls_yflux_source
*  assemble_cont_vel_source
*  assemble_extv_kinematic
*  assemble_interface_extension_velocity_sic 
*  assemble_eik_kinematic
*  assemble_p_source 
*  assemble_precoil_source
*  assemble_uvw_source
*  assemble_extension_velocity_path_dependence
*  assemble_fill_path_dependence
*  assemble_energy_path_dependence
*  assemble_momentum_path_dependence
*  assemble_continuity_path_dependence
*  assemble_LM_source 
*  fluid_stress                 void       assemble_momentum()
*                                          assemble_p_source()
*                                          assemble_precoil_source() 
*                                          assemble_uvw_source()
*                                          assemble_momentum_path_dependence()
*  heat_flux                    void       assemble_energy()
*                                          assemble_t_source()
*                                          assemble_energy_path_dependence()
*  heat_source                 double
*  ls_modulate_heatsource 
*  calc_pspg
*  calc_cont_gls
*  assemble_ls_latent_heat_source
*  assemble_acoustic		int 
*  assemble_acoustic_reynolds_stress int
*  accoustic_flux              void
*  assemble_ars_source          int
*  visc_diss_acoustic_source    double
******************************************************************************/

/*  _______________________________________________________________________  */

/* assemble_mesh -- assemble terms (Residual &| Jacobian) for mesh stress eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *	lec-- pointer to Local Element Contribution	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	Sat Mar 19 16:07:51 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_mesh (double time, 
	       double tt,
	       double dt,
	       int ielem,              /* current element number */
	       int ip,                 /* current integration point */
	       int ip_total)           /* total gauss integration points */

{
  int eqn, peqn, var, pvar;
  const int dim   = pd->Num_Dim;
  int p, q, a, b;

  int w;

  int i, j;
  int status, err;

  dbl det_J;
  dbl d_det_J_dmeshbj;
  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl VV[DIM][DIM];		/* Inertia Tensor... */

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  dbl TT[DIM][DIM];		/* Mesh stress tensor... */

  dbl dTT_dx[DIM][DIM][DIM][MDE];	/* Sensitivity of stress tensor...
                           to nodal displacements */
  dbl dTT_dp[DIM][DIM][MDE];	/* Sensitivity of stress tensor...
                           to nodal pressure*/
  dbl dTT_dc[DIM][DIM][MAX_CONC][MDE];	/* Sensitivity of stress tensor...
                           to nodal concentration */
  dbl dTT_dp_liq[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				    to nodal porous liquid pressure*/
  dbl dTT_dp_gas[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				    to nodal porous gas pressure*/
  dbl dTT_dporosity[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				    to nodal porosity*/
  dbl dTT_dsink_mass[DIM][DIM][MDE];/* Sensitivity of stress tensor...
				       to nodal sink_mass*/
  dbl dTT_dT[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
                           to temperature*/
  dbl dTT_dmax_strain[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
                           to max_strain*/
  dbl dTT_dcur_strain[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
                           to cur_strain*/

  dbl mu;
  dbl lambda, rho;

  dbl g[DIM];			/* Mesh body force. */

  dbl diff_a, diff_b, diff_c;	/* Temporary variables hold partially */
				/* constructed diffusion terms... */
  dbl diffusion=0., advection=0., advect_a=0., advect_b=0., advect_c=0.;
  dbl source=0., mass=0.;
  dbl wt;

  /*
   * Galerkin weighting functions for i-th and a-th mesh stress residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl (* grad_phi_i_e_a)[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  /*
   * Mesh derivatives...
   */
  dbl (* dgradphi_i_e_a_dmesh)[DIM][DIM][MDE];	/* for specific (i, a, b, j) !!  */

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Basis_Functions *bfm;

  /*finally we have some time-derivative quantities for the Newmark scheme */
  dbl x_c[DIM];          /*base coordinates of deformed mesh */
  dbl x_old[DIM];        /*old value of base coordinates on deformed mesh */
  dbl x_dot_old[DIM];    /*Old value of mesh velocity */
  dbl x_dbl_dot_old[DIM];/*old value of acceleration */
  dbl x_dbl_dot[DIM];    /*current value of mesh acceleration */
  dbl newmark_beta = 0.0;       /*Newmark scheme beta value. */

  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on =0;
  int source_on =0;
  int diffusion_on =0;

  dbl mass_etm, advection_etm, diffusion_etm, source_etm;
  
  dbl d_area;
  
  dbl * _J;
  


  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */



  eqn   = R_MESH1;			/* Well, yes, there really are 3, */
					/* but for now just use the first */
					/* to tell whether there is anything */
					/* to do...*/
  wt    = fv->wt;
  h3 = fv->h3;			/* Differential volume element (scales). */

  /*
   * Bail out fast if there's nothing to do...
   * Also bail out if we are in an a shell equation. In that case Num_Dim
   *      is equal to the problem dimension, but ielem_dim is one less than the
   *      problem dimension.
   */
  if (! pd->e[eqn] || (ei->ielem_dim < pd->Num_Dim))
    {
      return(status);
    }

  det_J = bf[eqn]->detJ;
  d_area = det_J*wt*h3;

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in mm_fill.
   */
  rho	 = density(d_rho, time);
  mu	 = elc->lame_mu;
  lambda = elc->lame_lambda;

   if(mp->MeshSourceModel == CONSTANT)
     {
       for (a = 0; a < dim; a++)
	 {
	   g[a] = mp->mesh_source[a];
	 }
     }
   else
     {
       EH(-1, "Unrecognized MeshSourceModel");
     }

  /*
   * Get the deformation gradients and tensors if needed
   */
  err = belly_flop(mu);
  EH(err, "error in belly flop");
  if (err == 2) return(err);

  /*
   * If using volume change metric for element quality,
   * sum local contribution here.
   */
  if (nEQM>0 && eqm->do_vol && af->Assemble_Residual)
    {
      eqm->vol_sum += fv->volume_change;
      if (fv->volume_change < eqm->vol_low)
        {
          eqm->vol_low = fv->volume_change;
        }
      eqm->vol_count++;
    }

  /*
   * Total mesh stress tensor...
   */
  /* initialize some arrays */
  memset( TT, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian) {
    memset( dTT_dx,         0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset( dTT_dp,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dc,         0, sizeof(double)*DIM*DIM*MAX_CONC*MDE);
    memset( dTT_dp_liq,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dp_gas,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dporosity,  0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dsink_mass, 0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dT,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dmax_strain,0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dcur_strain,0, sizeof(double)*DIM*DIM*MDE);
  }

  err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq,
			   dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
			   mu, lambda, dt, ielem, ip, ip_total);

  /* Calculate inertia of mesh if required. PRS side note: no sensitivity
   * with respect to porous media variables here for single compont pore liquids.
   * eventually there may be sensitivity for diffusion in gas phase to p_gas */
  if ( cr->MeshMotion == LAGRANGIAN ||
       cr->MeshMotion == DYNAMIC_LAGRANGIAN)
    {
      err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    }
  else /* No inertia in an Arbitrary Mesh */
    {
      memset( vconv, 0, sizeof(double)*MAX_PDIM );
      if (pd->v[MESH_DISPLACEMENT1])
	memset(d_vconv->X, 0, DIM*DIM*MDE*sizeof(dbl));


      if (pd->v[VELOCITY1] || pd->v[POR_LIQ_PRES])
	memset(d_vconv->v, 0, DIM*DIM*MDE*sizeof(dbl));

      if (pd->v[MASS_FRACTION] || pd->v[POR_LIQ_PRES])
	memset(d_vconv->C, 0, DIM*MAX_CONC*MDE*sizeof(dbl));


      if (pd->v[TEMPERATURE])
	memset(d_vconv->T, 0, DIM*MDE*sizeof(dbl) );

    }
  for ( a=0; a<dim; a++)
    {
      for ( b=0; b<dim; b++)
	{
	  VV[a][b] = vconv[a] * vconv[b];
	}
    }

  if ( tran->solid_inertia)
    {
      newmark_beta = tran->newmark_beta;
      for (a=0; a<dim; a++)
	{
	  x_c[a]      = fv->x[a];
	  x_old[a]    = fv_old->x[a];
	  x_dot_old[a]    = fv_dot_old->x[a];
	  x_dbl_dot_old[a] = fv_dot_dot_old->x[a];

	  /*generalized Newmark equations here */
	  /*Notice Gama is not needed here, unless we have a damping term */
	  x_dbl_dot[a] =  (x_c[a] - x_old[a])/newmark_beta/dt/dt
	    -x_dot_old[a]*(1.0/newmark_beta/dt)
	    -x_dbl_dot_old[a]*(1-2.*newmark_beta)/2./newmark_beta;
	}
    }
  else
    {
      for (a=0; a<dim; a++)
	{
	  x_dbl_dot[a]    = 0.;
	}
    }
  
  /*
   * Residuals_________________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      /*
       * Assemble each component "a" of the momentum equation...
       */
      for ( a=0; a<dim; a++)
	{
	  eqn = R_MESH1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];
	  
	  mass_on = pd->e[eqn] & T_MASS;
	  advection_on =  pd->e[eqn] & T_ADVECTION ;
	  diffusion_on = pd->e[eqn] & T_DIFFUSION;
	  source_on = pd->e[eqn] & T_SOURCE;
	  
	  mass_etm = pd->etm[eqn][(LOG2_MASS)];
	  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
	  diffusion_etm = pd->etm[eqn][(LOG2_DIFFUSION)];
	  source_etm = pd->etm[eqn][(LOG2_SOURCE)];


	  for ( i=0; i<ei->dof[eqn]; i++)
	    {
	      phi_i = bfm->phi[i];

	      grad_phi_i_e_a = bfm->grad_phi_e[i][a] ;
	      
	      mass = 0.;
	      if ( transient_run &&
		   pd_glob[ei->mn]->MeshMotion == DYNAMIC_LAGRANGIAN)
		{
		  if ( pd->e[eqn] & T_MASS )
		    {
		      mass  = -x_dbl_dot[a];
		      mass *= phi_i * rho * d_area;
		      mass *= mass_etm;
		    }
		}	
      
	      diffusion = 0.;
	      if ( diffusion_on )
		{
		  if (cr->MeshMotion == ARBITRARY)
		    {
		      /*
		       * use pseudo cartesian arbitrary mesh motion 
		       */
		      for ( p=0; p<dim; p++)
			{
			  diffusion += bfm->d_phi[i][p] * TT[a][p];
			}
		      diffusion *= - det_J * wt;
		    } else {
 		      for ( q=0; q<VIM; q++) 
 			{ 
			  for ( p=0; p<VIM; p++) 
			    { 
			      diffusion += grad_phi_i_e_a[p][q] * TT[q][p]; 
			    } 
			}
		      diffusion *= - d_area;
		    }

		  diffusion *= diffusion_etm;
		}

	      /* add inertia of moving mesh */
	      advection = 0.;
	      if ( advection_on )
		{
		  WH(-1,"Warning: mesh advection unverified.\n");
		  for ( p=0; p<dim; p++)
		    {
		      for ( q=0; q<dim; q++)
			{

			  advection += grad_phi_i_e_a[p][q] * VV[q][p];
			}
		    }
		  advection *= rho * d_area;
		  /* Advection term only applies in Lagrangian mesh motion */
		  advection *= advection_etm;
		}

	      source = 0.;
	      if ( source_on )
		{
		  source = phi_i * g[a] * d_area;
		  /* Source term only applies in Lagrangian mesh motion */

		  source *= source_etm;
		}

	      /*
	       * porous term removed for mesh equation
	       *  - the additional effects due  to porosity are entered
	       *    into the consitutive equation for stress
	       */

	      lec->R[peqn][i] += mass + diffusion + source + advection;
	    }
	}
    } 

  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      for ( a=0; a<dim; a++)
	{
	  eqn = R_MESH1+a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  mass_on = pd->e[eqn] & T_MASS;
	  advection_on =  pd->e[eqn] & T_ADVECTION ;
	  diffusion_on = pd->e[eqn] & T_DIFFUSION;
	  source_on = pd->e[eqn] & T_SOURCE;
	  
	  mass_etm = pd->etm[eqn][(LOG2_MASS)];
	  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
	  diffusion_etm = pd->etm[eqn][(LOG2_DIFFUSION)];
	  source_etm = pd->etm[eqn][(LOG2_SOURCE)];


	  for ( i=0; i<ei->dof[eqn]; i++)
	    {
	      phi_i = bfm->phi[i];
	      grad_phi_i_e_a = bfm->grad_phi_e[i][a] ;
	      dgradphi_i_e_a_dmesh = bfm->d_grad_phi_e_dmesh[i][a];

	      /*
	       * Set up some preliminaries that are needed for the (a,i)
	       * equation for bunches of (b,j) column variables...
	       */

	      /*
	       * J_d_d
	       */
	      for ( b=0; b<dim; b++)
		{
		  var = MESH_DISPLACEMENT1+b;
		  if ( pd->v[var] )
		    {
		      pvar = upd->vp[var];
			  _J = lec->J[peqn][pvar][i];

		      for ( j=0; j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  
			  d_det_J_dmeshbj = bfm->d_det_J_dm[b][j];

			  dh3dmesh_bj = fv->dh3dq[b] * phi_j;

			  mass = 0.;

			  if ( transient_run &&
			       pd_glob[ei->mn]->MeshMotion == DYNAMIC_LAGRANGIAN )
			    {
			      if ( mass_on )
				{
				  mass = -(double)delta(a,b)*phi_j/newmark_beta/dt/dt;
				  mass *= phi_i * rho * d_area;
				  mass -= x_dbl_dot[a]*phi_i*rho*wt*
				    (d_det_J_dmeshbj * h3 + dh3dmesh_bj * det_J);
				  mass *= mass_etm;
				}
				  
			    }

			  diffusion = 0.;			  
			  /* Three parts:
			   *   diff_a = Int ( d(grad(phi_i e_a))/dmesh : Pi |Jv| )
			   *   diff_b = Int ( grad(phi_i e_a) : d(Pi)/dmesh |Jv| )
			   *   diff_c = Int ( grad(phi_i e_a) : Pi d(|Jv|)/dmesh )
			   */
			  if ( diffusion_on )
			    {
			      if (cr->MeshMotion == ARBITRARY)
				{
				  /*
				   * use pseudo cartesian arbitrary mesh motion 
				   */
				  diff_a = 0.;
				  diff_b = 0.;				  
				  diff_c = 0.;

				  for ( p=0; p<dim; p++)
				    {
				      diff_a += 
					bfm->d_d_phi_dmesh[i][p][b][j] * TT[a][p];

				      diff_b += 
					bfm->d_phi[i][p] * dTT_dx[a][p][b][j]; 

				      diff_c += 
					bfm->d_phi[i][p] * TT[a][p]; 
				    }

				  diff_a *= -det_J * wt;
				  diff_b *= -det_J * wt;
				  diff_c *= - d_det_J_dmeshbj * wt;

				} else {
				  /* LAGRANGIAN MESH - Solve in curvilinear coordinates */
					  
				  diff_a = 0.;
				  diff_c = 0.;
				  diff_b = 0.;		
		   
				  for ( p=0; p<VIM; p++) 
				    { 
				      for ( q=0; q<VIM; q++) 
					{ 
					  diff_a +=  
					    dgradphi_i_e_a_dmesh[p][q][b][j] * TT[q][p]; 

					  diff_b +=  
					    grad_phi_i_e_a[p][q] * dTT_dx[q][p][b][j];  

					  diff_c +=  
					    grad_phi_i_e_a[p][q] * TT[q][p];  
					}
				    }
				  diff_a *= -d_area;

				  diff_b *= -d_area;

				  diff_c *= - wt * ( d_det_J_dmeshbj * h3 +
						     det_J * dh3dmesh_bj );


				}

			      diffusion = diff_a + diff_b + diff_c;
			      
			      diffusion *= diffusion_etm;
			    }
							  
			  /* add inertia of moving mesh */
			  advect_a = 0.;
			  advect_b = 0.;
			  advect_c = 0.;
			  if ( advection_on )
			    {
			      /*
			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      dgradphi_i_e_a_dmeshbj[p][q] = 
					bfm->
					d_grad_phi_e_dmesh[i][a] [p][q] [b][j];
				    }
				}
				*/
			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      

				      advect_a += grad_phi_i_e_a[p][q] * VV[q][p];
				      advect_b += grad_phi_i_e_a[p][q] * (
					     d_vconv->X[q][b][j] * vconv[p] +
					     vconv[q] *            d_vconv->X[p][b][j] );
				      advect_c += dgradphi_i_e_a_dmesh[p][q][b][j] * VV[q][p];
				    }
				}

			      advect_a *=  rho * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
			      advect_b *=  rho * det_J * wt;
			      advect_c *=  rho * det_J * wt;

			      advection = advect_a + advect_b + advect_c;
			      advection *= advection_etm;
			    }

			  source = 0.;

			  if ( source_on )
			    {
			      source = phi_i * g[a] * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
			      source *= source_etm;
			    }
			  
			  _J[j] += mass + diffusion + source + advection;
			/*  lec->J[peqn][pvar][i][j] += mass + diffusion + source + advection;*/
			}
		    }
		}

/*
 * J_d_P
 */
/* For Lagrangian Mesh, add in the sensitivity to pressure  */
	      if (pd->MeshMotion == LAGRANGIAN ||
		  pd->MeshMotion == DYNAMIC_LAGRANGIAN)
		{
		  var = PRESSURE;
		  if (pd->v[var])
		    {
		      pvar = upd->vp[var];
			  _J = lec->J[peqn][pvar][i];

		      for ( j=0; j<ei->dof[var]; j++)
			{
			  diffusion = 0.;
			  phi_j = bf[var]->phi[j];
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{				  
				  diffusion += grad_phi_i_e_a[p][q] * dTT_dp[q][p][j];
				}
			    }
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;

			  _J[j] += diffusion;
			 /* lec->J[peqn][pvar][i][j] += diffusion;*/
			}
		    }
		}
/*
 * J_d_c
 */
	      var = MASS_FRACTION;
	      if (pd->v[var])
		{

		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      for ( w=0; w<pd->Num_Species_Eqn; w++)
			{
				_J = lec->J[peqn][MAX_PROB_VAR+w][i];

			  /* add inertia of moving mesh */
			  advection = 0.;
			  /* For Lagrangian Mesh with advection, add in the 
			     sensitivity to concentration  */
			  if ((cr->MeshMotion == LAGRANGIAN ||
			       cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
			          ( pd->e[eqn] & T_ADVECTION ) ) 
			    {
			      advect_a = 0.;
			      advect_b = 0.;

			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      
				      advect_a += grad_phi_i_e_a[p][q] * (
						d_vconv->C[q][w][j] * vconv[p] +
						vconv[q] * d_vconv->C[p][w][j] );
				    }
				}
			      advect_a *= rho * det_J * wt * h3;
			      advect_a *= pd->etm[eqn][LOG2_ADVECTION];

			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
				    }
				}
			      advect_b *= d_rho->C[w][j] * det_J * wt * h3;
			      advect_b *= pd->etm[eqn][LOG2_ADVECTION];
			      advection = advect_a + advect_b;
			    }

			  diffusion = 0.;
			  if ((cr->MeshMotion == LAGRANGIAN ||
			       cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
			          ( pd->e[eqn] & T_DIFFUSION ) )
			    { 
			      phi_j = bf[var]->phi[j];
			      for ( p=0; p<VIM; p++)
				{
				  for ( q=0; q<VIM; q++)
				    {
				      grad_phi_i_e_a[p][q] = 
					bfm->grad_phi_e[i][a] [p][q];
				      
				      diffusion += grad_phi_i_e_a[p][q] 
					* dTT_dc[q][p][w][j];
				    }
				}
			      diffusion *= - det_J * wt * h3;
			      diffusion *= pd->etm[eqn][LOG2_DIFFUSION];
			    }

			  _J[j] += advection + diffusion;
			  /*lec->J[peqn][MAX_PROB_VAR + w][i][j] += advection + diffusion;*/
			}
		    }
		}
	      
/*
 * J_d_p_liq
 */
		  
		if (pd->v[POR_LIQ_PRES] &&  
		    (cr->MeshMotion == LAGRANGIAN ||
		     cr->MeshMotion == DYNAMIC_LAGRANGIAN) )
		{
		  var = POR_LIQ_PRES;
		  pvar = upd->vp[var];
		  _J = lec->J[peqn][pvar][i];
		  
			for ( j=0; j<ei->dof[var]; j++)
		    {
		      advection = 0.;
		      /* For Lagrangian Mesh with advection, add in the 
			 sensitivity to concentration. N.B. this is zero
			 for now. The only way this will be nonzero is if
			 we add in diffusion-induced velocity in each porous
			 phase. This may happen in multicomponent cases, but not
			 now. PRS 5/7/01 */
		      if  ( advection_on) 
			{
			  advect_a = 0.;
			  advect_b = 0.;
			  
			/*  for ( p=0; p<dim; p++)
			    {
			      for ( q=0; q<dim; q++)
				{
				  
				  advect_a += grad_phi_i_e_a[p][q] * (
				    dvconv_d_p_gas[q][w][j] * vconv[p] +
				    vconv[q] * dvconv_d_p_gas[p][w][j] ); 
				}
			    }*/
			  /* advect_a *= rho * det_J * wt * h3;
			     advect_a *= pd->etm[eqn][LOG2_ADVECTION]; */

			 /*  for ( p=0; p<dim; p++)
			    {
			      for ( q=0; q<dim; q++)
				{
				  advect_b += grad_phi_i_e_a[p][q] * VV[q][p]; 
				}
			    }*/
			  /* advect_b *= d_rho->C[w][j] * det_J * wt * h3;
			     advect_b *= pd->etm[eqn][LOG2_ADVECTION];
			     advection = advect_a + advect_b;  */
			}

		      diffusion = 0.;
		      if (  diffusion_on )
			{ 
			  phi_j = bf[var]->phi[j];
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  grad_phi_i_e_a[p][q] = 
				    bfm->grad_phi_e[i][a] [p][q];
				  
				  diffusion += grad_phi_i_e_a[p][q] 
				    * dTT_dp_liq[q][p][j];
				}
			    }
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;
			}
		      
		      _J[j] += advection + diffusion;
		     /* lec->J[peqn][pvar][i][j] += advection + diffusion;*/
		    }
		}
/*
 * J_d_p_gas
 */
	      if (pd->v[POR_GAS_PRES] &&
		      (cr->MeshMotion == LAGRANGIAN ||
			   cr->MeshMotion == DYNAMIC_LAGRANGIAN) )
		{
	      var = POR_GAS_PRES;
	      pvar = upd->vp[var];
		  
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      advection = 0.;
		      /* For Lagrangian Mesh with advection, add in the 
			 sensitivity to concentration. N.B. this is zero
			 for now. The only way this will be nonzero is if
			 we add in diffusion-induced velocity in each porous
			 phase. This may happen in multicomponent cases, but not
			 now. PRS 5/7/01 */
		      if ( advection_on) 
			{
			  advect_a = 0.;
			  advect_b = 0.;
			  
			/*  for ( p=0; p<dim; p++)
			    {
			      for ( q=0; q<dim; q++)
				{
				  
				  advect_a += grad_phi_i_e_a[p][q] * (
				    dvconv_d_p_gas[q][w][j] * vconv[p] +
				    vconv[q] * dvconv_d_p_gas[p][w][j] ); 
				}
			    }*/
			  /* advect_a *= rho * det_J * wt * h3;
			     advect_a *= pd->etm[eqn][LOG2_ADVECTION]; */
			  
			/*  for ( p=0; p<dim; p++)
			    {
			      for ( q=0; q<dim; q++)
				{
				   advect_b += grad_phi_i_e_a[p][q] * VV[q][p]; 
				}
			    }*/
			  /* advect_b *= d_rho->C[w][j] * det_J * wt * h3;
			     advect_b *= advection_etm;
			     advection = advect_a + advect_b;  */
			}
		      
		      diffusion = 0.;
		      if ( diffusion_on )
			{ 
			  phi_j = bf[var]->phi[j];
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  grad_phi_i_e_a[p][q] = 
				    bfm->grad_phi_e[i][a] [p][q];
				  
				  diffusion += grad_phi_i_e_a[p][q] 
				    * dTT_dp_gas[q][p][j];
				}
			    }
			  diffusion *= - det_J * wt * h3;
			  diffusion *= diffusion_etm;
			}
		      
		      lec->J[peqn][pvar][i][j] += advection + diffusion;
		    }
		}
/*
 * J_d_porosity
 */
	      if (pd->v[POR_POROSITY])
		{
	      var = POR_POROSITY;
	      pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      diffusion = 0.;
		      if ((cr->MeshMotion == LAGRANGIAN ||
			   cr->MeshMotion == DYNAMIC_LAGRANGIAN) && 
			  ( pd->e[eqn] & T_DIFFUSION ) )
			{ 
			  phi_j = bf[var]->phi[j];
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  grad_phi_i_e_a[p][q] = 
				    bfm->grad_phi_e[i][a] [p][q];
				  
				  diffusion += grad_phi_i_e_a[p][q] 
				    * dTT_dporosity[q][p][j];
				}
			    }
			  diffusion *= - d_area;
			  diffusion *= pd->etm[eqn][LOG2_DIFFUSION];
			}
		      
		      lec->J[peqn][pvar][i][j] += diffusion;
		    }
		}
/*
 * J_d_sink_mass
 */
	      if (pd->v[POR_SINK_MASS])
		{
		  var = POR_SINK_MASS;
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      diffusion = 0.;
		      if ((cr->MeshMotion == LAGRANGIAN ||
			   cr->MeshMotion == DYNAMIC_LAGRANGIAN) && 
			  ( pd->e[eqn] & T_DIFFUSION ) )
			{ 
			  phi_j = bf[var]->phi[j];
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  grad_phi_i_e_a[p][q] = 
				    bfm->grad_phi_e[i][a] [p][q];
				  
				  diffusion += grad_phi_i_e_a[p][q] 
				    * dTT_dsink_mass[q][p][j];
				}
			    }
			  diffusion *= - d_area;
			  diffusion *= pd->etm[eqn][LOG2_DIFFUSION];
			}
		      
		      lec->J[peqn][pvar][i][j] +=  diffusion;
		    }
		}
/*
 * J_d_T
 */
 	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&  pd->e[eqn]  ) 
		{
		  var = TEMPERATURE;
		  if (pd->v[var])
		    {
		      pvar = upd->vp[var];

		      for ( j=0; j<ei->dof[var]; j++)
			{
 			  diffusion = 0.;
 			  phi_j = bf[var]->phi[j];
 			  for ( p=0; p<VIM; p++)
 			    {
 			      for ( q=0; q<VIM; q++)
 				{				  
				  diffusion += grad_phi_i_e_a[p][q] * dTT_dT[q][p][j];
 				}
 			    }
 			  diffusion *= - d_area;
 			  diffusion *= diffusion_etm;
 
			  /* add inertia of moving mesh */
			  advection = 0.;
			  if ( advection_on )
			    {
			      /*
			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      grad_phi_i_e_a[p][q] = 
					bfm->grad_phi_e[i][a] [p][q];

				    }
				}
				*/

			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      
				      advect_a += grad_phi_i_e_a[p][q] * (
						  d_vconv->T[q][j] * vconv[p] +
						  vconv[q] *         d_vconv->T[p][j] );
				    }
				}
			      advect_a *=  rho * d_area;
			      advect_a *= advection_etm;

			      for ( p=0; p<dim; p++)
				{
				  for ( q=0; q<dim; q++)
				    {
				      advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
				    }
				}
			      advect_b *= d_rho->T[j] * d_area;
			      advect_b *= advection_etm;
			     
			      advection = advect_a + advect_b;
			    }

 			  lec->J[peqn][pvar][i][j] += advection + diffusion;
			}
		    }
		}
/*
 * J_d_max_strain
 */
 	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&  pd->e[eqn]  ) 
		{
		  var = MAX_STRAIN;
		  if (pd->v[var])
		    {
		      pvar = upd->vp[var];

		      for ( j=0; j<ei->dof[var]; j++)
			{

 			  diffusion = 0.0;
			  for ( q=0; q<VIM; q++) 
			    { 
			      for ( p=0; p<VIM; p++) 
				{ 
				  diffusion += grad_phi_i_e_a[p][q] * dTT_dmax_strain[q][p][j]; 
				} 
			    }
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;

 			  lec->J[peqn][pvar][i][j] +=  diffusion;
			}
		    }
		}
/*
 * J_d_cur_strain
 */
 	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&  pd->e[eqn]  ) 
		{
		  var = CUR_STRAIN;
		  if (pd->v[var])
		    {
		      pvar = upd->vp[var];

		      for ( j=0; j<ei->dof[var]; j++)
			{

 			  diffusion = 0.0;
			  for ( q=0; q<VIM; q++) 
			    { 
			      for ( p=0; p<VIM; p++) 
				{ 
				  diffusion += grad_phi_i_e_a[p][q] * dTT_dcur_strain[q][p][j]; 
				} 
			    }
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;

 			  lec->J[peqn][pvar][i][j] +=  diffusion;
			}
		    }
		}

	    }  /* end of loop over equations i  */
	} /* end of loop over equation directions a */
    } /* end of if jacobian */

  return(status);
}

/*  _______________________________________________________________________  */

/* assemble_energy -- assemble terms (Residual &| Jacobian) for energy eqns
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
 * Revised:     10/98 by KSC to enable the calculation of source term due to electrode kinetics as in thermal batteries.
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_energy(double time,	/* present time value */
		double tt,	/* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0) */
		double dt,	/* current time step size */
		const PG_DATA *pg_data )
{
  int eqn, var, peqn, pvar, dim, p, a, b, qq, w, i, j, status;

  dbl T_dot;				/* Temperature derivative wrt time. */

  dbl q[DIM];				/* Heat flux vector. */
  HEAT_FLUX_DEPENDENCE_STRUCT d_q_struct; /* Heat flux dependence. */
  HEAT_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  dbl rho;				/* Density (no variations allowed
					   here, for now) */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Thermal conductivity dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl Cp;				/* Heat capacity. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT d_Cp_struct; /* Heat capacity dependence. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp = &d_Cp_struct;

  dbl h;				/* Heat source. */
  HEAT_SOURCE_DEPENDENCE_STRUCT d_h_struct; /* Heat source dependence. */
  HEAT_SOURCE_DEPENDENCE_STRUCT *d_h = &d_h_struct;

  int v_s[MAX_MODES][DIM][DIM];
  int mode;

  dbl mass;		         	/* For terms and their derivatives */

  dbl advection;			/* For terms and their derivatives */
  dbl advection_a;
  dbl advection_b;
  dbl advection_c;
  dbl advection_d;
  dbl advection_e;
  dbl advection_f;

  dbl diffusion;
  dbl diff_a, diff_b, diff_c, diff_d;
  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

/* SUPG variables */
  dbl h_elem=0, h_elem_inv=0, h_elem_deriv=0, h_elem_inv_deriv=0.;
  dbl supg, d_wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM];

  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM];		/* ditto.  */
  dbl wt;

  dbl *grad_T = fv->grad_T;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  int err;
  
  const double *hsquared = pg_data->hsquared ;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
	  centroid velocity for Q2 and the average
	  of the vertices for Q1. It comes from
	  the routine "element_velocity." */

  /* initialize grad_phi_j */
  for (i = 0; i < DIM; i++) {
    grad_phi_j[i] = 0;
  }

  /*   static char yo[] = "assemble_energy";*/

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;

  eqn   = R_ENERGY;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  wt = fv->wt;				/* Gauss point weight. */

  h3 = fv->h3;			/* Differential volume element. */

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  supg = 0.;

  if( mp->Ewt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Ewt_funcModel == SUPG )
    {
      if( !pd->e[R_MOMENTUM1])
        EH(-1, " must have momentum equation velocity field for energy equation upwinding. You may want to turn it off");
      supg = mp->Ewt_func;
    }

  (void) stress_eqn_pointer(v_s);

  if(supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p] * vcent[p] / hsquared[p];
	}
      h_elem = sqrt(h_elem)/2.;
      if(h_elem == 0.)
	{
	  h_elem_inv=0.;
	}
      else
	{
	  h_elem_inv=1./h_elem;
	}

    }
/* end Petrov-Galerkin addition */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are rho, Cp, k, and h.  For now we will
   * take rho as constant.   Cp, h, and k we will allow to vary with temperature,
   * spatial coordinates, and species concentration.
   */

  rho = density(d_rho, time);

  /* CHECK FOR REMOVAL */
  conductivity( d_k, time );

  Cp = heat_capacity( d_Cp, time );

  h = heat_source( d_h, time, tt, dt );

  if (pd->TimeIntegration != STEADY) {
    T_dot = fv_dot->T;
  } else {
    T_dot = 0.0;
  }

  heat_flux( q, d_q, time );

/* get the convection velocity (it's different for arbitrary and
   lagrangian meshes) */
  if( cr->MeshMotion == ARBITRARY ||
      cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN)
    {
      err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
      EH(err, "Error in calculating effective convection velocity");
    }
  else if ( cr->MeshMotion == TOTAL_ALE)
    {
      err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, dt, tt);
      EH(err, "Error in calculating effective convection velocity_rs");
    }


  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      eqn = R_ENERGY;
      peqn = upd->ep[eqn];
      var = TEMPERATURE;
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  mass = 0.;
	  if ( pd->TimeIntegration != STEADY )
	    {
	      if ( pd->e[eqn] & T_MASS )
		{
		  mass  = T_dot;
		  mass *= - phi_i * rho * Cp * det_J * wt;
		  mass *= h3;
		  mass *= pd->etm[eqn][(LOG2_MASS)];
		}
	    }

	  /* only use Petrov Galerkin on advective term - if required */
	  wt_func = bf[eqn]->phi[i];
	  /* add Petrov-Galerkin terms as necessary */
	  if(supg!=0.)
	    {
	      for(p=0; p<dim; p++)
		{
		  wt_func += supg * h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
		}
	    }

	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {

	      for ( p=0; p<VIM; p++)
		{
		  advection += vconv[p]     * grad_T[p];
		}

	      advection *= - wt_func * rho * Cp * det_J * wt;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }

	  diffusion = 0.;
	  if ( pd->e[eqn] & T_DIFFUSION )
	    {
	      for ( p=0; p<VIM; p++)
		{
		  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
		}

	      for ( p=0; p<VIM; p++)
		{
		  diffusion += grad_phi_i[p] * q[p];
		}
	      diffusion *= det_J * wt;
	      diffusion *= h3;
	      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  source = 0.;
	  if ( pd->e[eqn] & T_SOURCE )
	    {
	      source += phi_i * h * det_J * wt;
	      source *= h3;
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	    }

	  lec->R[peqn][i] +=
	    mass + advection +  diffusion + source;

	}
    }


  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      eqn   = R_ENERGY;
      peqn = upd->ep[eqn];
      for ( i=0; i<ei->dof[eqn]; i++)
	{
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  wt_func = bf[eqn]->phi[i];
	  /* add Petrov-Galerkin terms as necessary */
	  if(supg!=0.)
	    {
	      for(p=0; p<dim; p++)
		{
		  wt_func += supg * h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
		}
	    }

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  for ( p=0; p<VIM; p++)
	    {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
	    }

	  /*
	   * J_e_T
	   */
	  var = TEMPERATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  mass = 0.;
		  if ( pd->TimeIntegration != STEADY )
		    {
		      if ( pd->e[eqn] & T_MASS )
			{
			  mass  =  rho         * d_Cp->T[j] * T_dot
			    +      d_rho->T[j] * Cp       * T_dot
			    +      rho *  Cp * (1 + 2. * tt) * phi_j / dt;
			  mass *= - phi_i * det_J * wt;
			  mass *= h3;
			  mass *= pd->etm[eqn][(LOG2_MASS)];
			}
		    }

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  advection +=  rho * d_Cp->T[j]     *
			                vconv[p] * grad_T[p]
			              + d_rho->T[j] * Cp   *
			                vconv[p] * grad_T[p]
				      + rho * Cp           *
					vconv[p] * grad_phi_j[p]
				      + rho * Cp           *
					d_vconv->T[p][j] * grad_T[p];
			}
		      advection *= - wt_func * det_J * wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  diffusion = 0.;
		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += d_q->T[p][j] * grad_phi_i[p];
			}
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->T[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }

	  /*
	   * J_e_V
	   */
	  var = VOLTAGE;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<dim; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  mass = 0.;

		  advection = 0.;

		  diffusion = 0.;

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->V[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }

	  /*
	   * J_e_v
	   */
	  for ( b=0; b<VIM; b++)
	    {
	      var = VELOCITY1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      mass = 0.;

		      advection = 0.;
		      advection_a = 0.;
		      if ( pd->e[eqn] & T_ADVECTION )
			{
			  advection_a += wt_func * rho * Cp * d_vconv->v[b][b][j] * grad_T[b];
			  advection_a *= - det_J * wt;
			  advection_a *= h3;
			  advection_a *= pd->etm[eqn][(LOG2_ADVECTION)];
			  if(supg!=0.)
			    {
			      h_elem_deriv = 0.;
			      if(hsquared[b] != 0.)
				{
				  h_elem_deriv =
				    vcent[b] * pg_data->dv_dnode[b][j] * h_elem_inv/4./hsquared[b];
				}
			      if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;

			    }
			  advection_b = 0.;
			  if(supg!=0.)
			    {
			      d_wt_func = supg * h_elem_inv * d_vconv->v[b][b][j] * grad_phi_i[b]
				+ supg * h_elem_inv_deriv * vconv[b] * grad_phi_i[b];

			      for(p=0;p<dim;p++)
				{
				  advection_b +=  rho * Cp * vconv[p] * grad_T[p];
				}

			      advection_b *=  d_wt_func;
			      advection_b *= - det_J * wt;
			      advection_b *= h3;
			      advection_b *= pd->etm[eqn][(LOG2_ADVECTION)];
			    }
			  advection = advection_a + advection_b;
			}

		      diffusion = 0.;
                      source = 0.;

		      if ( pd->e[eqn] & T_SOURCE )
			{
			  source += phi_i * d_h->v[b][j] * det_J * wt;
			  source *= h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;

		    }
		}
	    }

	  /*
	   * J_e_d_rs
	   */
	  for ( b=0; b<dim; b++)
	    {
	      var = SOLID_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

		      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];


		      

		      if(supg!=0.)
			{
			  h_elem_deriv = 0.;
			  h_elem_inv_deriv=0.;
			  for( qq=0; qq<dim; qq++ )
			    {
			      if(pg_data->hhv[qq][b] != 0.)
				{
				  h_elem_deriv -=
				    vcent[qq] * vcent[qq] *
				    pg_data->dhv_dxnode[qq][j] * pg_data->hhv[qq][b] * h_elem_inv/4./hsquared[qq]/hsquared[qq];
				}
			    }
			  if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;
			  // h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */

	      }

		      mass = 0.;

		      advection = 0.;

		      if ( pd->e[eqn] & T_ADVECTION )
			{
			  /*
			   * one parts:
			   *	Here, by "v" we mean "vconv" which is
			   *    vconv = v_real_solid - vmesh,   vmesh ~ xdot ~ dd/dt
			   *
			   *	(a)	Int( Cp * dvconv/dreal_solid . grad(T) h3 detJ )
			   *
			   */

			  advection_a = 0.;
			  if ( pd->e[eqn] & T_MASS )
			    {
			      for ( p=0; p<dim; p++)
				{ 
				  advection_a +=
				    d_vconv->X[p][b][j] * grad_T[p];
				}
			      advection_a *= -wt_func * rho * Cp * h3 *
				det_J * wt;
			    }

			  advection = advection_a;

			  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
			}


		      diffusion = 0.;

		      source = 0.;

		      lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		    }
		}
	    }
	  /*
	   * J_e_d
	   */
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

		      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];


		      if(supg!=0.)
			{
			  h_elem_deriv = 0.;
			  h_elem_inv_deriv=0.;
			  for( qq=0; qq<dim; qq++ )
			    {
			      if(pg_data->hhv[qq][b] != 0.)
				{
				  h_elem_deriv -=
				    vcent[qq] * vcent[qq] *
				    pg_data->dhv_dxnode[qq][j] * pg_data->hhv[qq][b] * h_elem_inv/4./hsquared[qq]/hsquared[qq];
				}
			    }
			  if (h_elem != 0.) h_elem_inv_deriv = -h_elem_deriv/h_elem/h_elem;
			  // h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */

			}

		      mass = 0.;
		      if ( pd->TimeIntegration != STEADY )
			{
			  if ( pd->e[eqn] & T_MASS )
			    {
			      mass  = T_dot;
			      mass *= - phi_i * rho  *
				(  Cp*              h3          * d_det_J_dmeshbj
				 + Cp*              dh3dmesh_bj * det_J
				 + d_Cp->X[b][j]*     h3          * det_J)
				  * wt;
			      mass *= pd->etm[eqn][(LOG2_MASS)];
			    }
			}


		      advection = 0.;

		      if ( pd->e[eqn] & T_ADVECTION )
			{
			  /*
			   * Four parts:
			   *	Here, by "v" we mean "vconv" which is
			   *    vconv = v - vmesh,   vmesh ~ xdot ~ dd/dt
			   *
			   *	d/dmesh [ Int(... vconv.grad(T) h3 detJ ) ]
			   *
			   *	(a)	Int( Cp * vconv.d(grad(T))/dmesh h3 detJ )
			   *	(b)	Int( Cp * vconv.grad(T) h3 ddetJ/dmesh )
			   *	(c)	Int( Cp * dvconv/dmesh . grad(T) h3 detJ )
			   *	(d)	Int( Cp * vconv.grad(T) dh3/dmesh detJ )
			   *	(e)	Int( dCp/dmesh * vconv.grad(T) h3 detJ )
			   *
			   */

			  advection_a = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      advection_a += vconv[p] *
				fv->d_grad_T_dmesh[p] [b][j];
			    }
			  advection_a *= -wt_func * rho * Cp * h3 * det_J * wt;

			  advection_b = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      advection_b += vconv[p] * grad_T[p];
			    }
			  advection_b *= - wt_func * rho * Cp * h3 *
			    d_det_J_dmeshbj *wt;

			  advection_c = 0.;
			  if ( pd->TimeIntegration != STEADY )
			    {
			      if ( pd->e[eqn] & T_MASS )
				{
				  for ( p=0; p<dim; p++)
				    {
				      advection_c +=
					d_vconv->X[p][b][j] * grad_T[p];
				    }
				  advection_c *= -wt_func * rho * Cp * h3 *
				    det_J * wt;
				}
			    }

			  advection_d = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      advection_d += vconv[p] * grad_T[p];
			    }
			  advection_d *= - wt_func * rho * Cp * dh3dmesh_bj *
			    det_J *wt;

                          advection_e = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      advection_e += vconv[p] * grad_T[p];
			    }

			  advection_e *= - wt_func * rho * d_Cp->X[b][j] * h3 *
			    det_J *wt;

			  advection_f = 0.;
			  if(supg != 0.)
			    {
			      d_wt_func = 0.;
			      for( p=0; p<dim; p++ )
				{
				  d_wt_func += supg
				    * (h_elem_inv *fv->v[p]* bf[eqn]->d_grad_phi_dmesh[i][p] [b][j]
				       +  h_elem_inv_deriv * fv->v[p] *grad_phi_i[p]);


				  advection_f += vconv[p] * grad_T[p];
				}
				  advection_f *= -d_wt_func * h3 * det_J * wt;
			    }

			  advection = advection_a + advection_b + advection_c
			    + advection_d + advection_e + advection_f;

			  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
			}

			  /*
			   * multiple parts:
			   * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
			   *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
			   *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
			   *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
			   */
		      diffusion = 0.;
                      if ( pd->e[eqn] & T_DIFFUSION )
		        {
                          diff_a = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      dgrad_phi_i_dmesh[p]
				= bf[eqn]->d_grad_phi_dmesh[i][p] [b][j];

			      diff_a += dgrad_phi_i_dmesh[p] * q[p];
			    }
			  diff_a *= det_J * h3 * wt;

                          diff_b = 0.;
                          for ( p=0; p<VIM; p++)
			    {
			      diff_b += d_q->X[p][b][j] * grad_phi_i[p];
			    }
		          diff_b *= det_J * h3 * wt;

			  diff_c = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      diff_c += grad_phi_i[p] * q[p];
			    }
			  diff_c *= d_det_J_dmeshbj * h3 * wt;

			  diff_d = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      diff_d += grad_phi_i[p] * q[p];
			    }
			  diff_d *= det_J * dh3dmesh_bj * wt;

			  diffusion = diff_a + diff_b + diff_c + diff_d;

			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}


		      source = 0.;

		      if ( pd->e[eqn] & T_SOURCE )
			{
			  source = phi_i * (h *          d_det_J_dmeshbj * h3 +
					    h *          det_J           * dh3dmesh_bj +
					    d_h->X[b][j]*  det_J           * h3) * wt;

			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		    }
		}
	    }

	  /*
	   * J_e_c
	   */
	  var = MASS_FRACTION;
	  if ( pd->e[eqn] && pd->v[var] )
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      mass = 0.;
		      if ( pd->TimeIntegration != STEADY )
			{
			  if ( pd->e[eqn] & T_MASS )
			    {
			      mass  = T_dot * (rho * d_Cp->C[w][j]  + d_rho->C[w][j] * Cp );
			      mass *= - phi_i * det_J * wt;
			      mass *= h3;
			      mass *= pd->etm[eqn][(LOG2_MASS)];
			    }
			}

		      advection = 0.;
		      if ( pd->e[eqn] & T_ADVECTION )
			{
			  for ( p=0; p<dim; p++)
			    {
			      advection +=  rho * Cp        * d_vconv->C[p][w][j] * grad_T[p];
			      advection += (rho * d_Cp->C[w][j]  + d_rho->C[w][j] * Cp)
				             * vconv[p] * grad_T[p] ;
			    }
			  advection *= - wt_func * h3 * det_J * wt;
			  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
			}

		      diffusion = 0.;
		      if ( pd->e[eqn] & T_DIFFUSION )
			{
			  for ( p=0; p<dim; p++)
			    {
			      diffusion += grad_phi_i[p] * d_q->C[p][w][j];
			    }
			  diffusion *= det_J * wt;
			  diffusion *= h3;
			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}

		      source = 0.;
		      if ( pd->e[eqn] & T_SOURCE )
			{
			  source += phi_i * d_h->C[w][j] * det_J * wt;
			  source *= h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][MAX_PROB_VAR + w][i][j] += advection + mass + diffusion + source;
		    }
		}
	    }


#ifdef COUPLED_FILL

	  /*
	   * J_e_F
	   */
	  var = FILL;
	  if ( pd->e[eqn] && pd->v[var] )
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  
		  mass = 0.;
		  if ( pd->TimeIntegration != STEADY )
		    {
		      if ( pd->e[eqn] & T_MASS )
			{
			  mass  = T_dot * (rho * d_Cp->F[j]  + d_rho->F[j] * Cp );
			  mass *= - phi_i * det_J * wt;
			  mass *= h3;
			  mass *= pd->etm[eqn][(LOG2_MASS)];
			}
		    }
		  
		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      for ( p=0; p<dim; p++)
			{
			  advection += (rho * d_Cp->F[j]  + d_rho->F[j] * Cp)
			    * vconv[p] * grad_T[p] ;
			}
		      advection *= - wt_func * h3 * det_J * wt;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }
		  
		  diffusion = 0.;
		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<dim; p++)
			{
			  diffusion += grad_phi_i[p] * d_q->F[p][j];
			}
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }
		  
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->F[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }
		  
		  lec->J[peqn][var][i][j] += advection + mass + diffusion + source;
		}
	    }
#endif /* COUPLED_FILL */


/**    add ve stress terms  **/
	  /*
	   * J_e_S
	   */
           for ( mode=0; mode<vn->modes; mode++)
              {
               for ( a=0; a<VIM; a++)
                 {
                   for ( b=0; b<VIM; b++)
                     {
                      var = v_s[mode][a][b];
	  	      if ( pd->v[var] )
	    		{
	      		pvar = upd->vp[var];
	      		for ( j=0; j<ei->dof[var]; j++)
			  {
		  	   source = 0.;
		  	   if ( pd->e[eqn] & T_SOURCE )
		    	    {
		      		source += phi_i * d_h->S[mode][a][b][j] * det_J * wt;
		      		source *= h3;
		      		source *= pd->etm[eqn][(LOG2_SOURCE)];
		    	    }
		  	lec->J[peqn][pvar][i][j] += source;
			}
	    	      }

		    }
		}
	    }
	  /*
	   * J_e_apr and api
	   */
	  var = ACOUS_PREAL;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  mass = 0.;
		  advection = 0.;
		  diffusion = 0.;
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->APR[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }
	  var = ACOUS_PIMAG;

	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  mass = 0.;
		  advection = 0.;
		  diffusion = 0.;
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->API[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }

	  var = LIGHT_INTP;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->INT[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] +=  source;
		}
	     }


	  var = LIGHT_INTM;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[LIGHT_INTM];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->INT[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }
		  lec->J[peqn][pvar][i][j] += source;
                }
            }


	  var = LIGHT_INTD;
	  if ( pd->v[var] )
	    {
	     pvar = upd->vp[var];
	     for ( j=0; j<ei->dof[var]; j++)
	       {
		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source += phi_i * d_h->INT[j] * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }
		  lec->J[peqn][pvar][i][j] += source;
               }
            }
	}
    }

  return(status);
} /* end of assemble_energy */




/*  _______________________________________________________________________  */

/* assemble_momentum -- assemble terms (Residual &| Jacobian) for momentum eqns
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
assemble_momentum(dbl time,       /* current time */
		  dbl tt,	  /* parameter to vary time integration from
				     explicit (tt = 1) to implicit (tt = 0) */
		  dbl dt,	  /* current time step size */
		  dbl h_elem_avg, /* average global element size for PSPG*/
		  const PG_DATA *pg_data, 
		  double xi[DIM],    /* Local stu coordinates */
		  const Exo_DB *exo)
{
#ifdef DEBUG_MOMENTUM_JAC
  int adx;
#endif

  int dim;
  //! wim is the length of the velocity vector
  int wim;
  int i, j, jk, p, q, a, b, c;

  int ledof, eqn, var, ii, peqn, pvar, w;
  
  int *pde = pd->e;
  int *pdv = pd->v;


  int status;
  struct Basis_Functions *bfm;

  dbl zero[3] = {0.0, 0.0, 0.0};        /* A zero array, for convenience. */

  /*dbl v_dot[DIM];*/		         	/* time derivative of velocity field. */
  dbl *v_dot;		         	/* time derivative of velocity field. */

  /*dbl x_dot[DIM];*/			/* current position field derivative wrt time. */
  dbl *x_dot;			/* current position field derivative wrt time. */

  dbl h3;			        /* Volume element (scale factors). */
  dbl dh3dmesh_bj;		        /* Sensitivity to (b,j) mesh dof. */

  /* field variables */
  /*dbl grad_v[DIM][DIM];*/		        /* Gradient of v. */
  dbl *grad_v[DIM];
  dbl *v = fv->v;

  dbl rho;				/* Density. */

  dbl f[DIM];				/* Body force. */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct;  /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  dbl det_J;                            /* determinant of element Jacobian */

  dbl d_det_J_dmesh_bj;			/* for specific (b,j) mesh dof */
  
  dbl d_area;

  dbl mass;			        /* For terms and their derivatives */
  dbl advection;
  dbl advection_a, advection_b, advection_c;
  dbl porous;
  dbl diffusion;
  dbl diff_a, diff_b, diff_c;
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
  dbl phi_i;
  dbl (* grad_phi_i_e_a ) [DIM] = NULL;
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl (* d_grad_phi_i_e_a_dmesh ) [DIM][DIM][MDE] = NULL;
  
  double * phi_i_vector, *phi_j_vector;

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  dbl wt;

  /* coefficient variables for the Brinkman Equation
     in flows through porous media: KSC on 5/10/95 */
  dbl por;         /* porosity of porous media */
  dbl por2;        /* square of porosity */
  dbl per = 0.0;         /* permeability of porous media */
  /* derivative of permeability wrt concentration */
  dbl d_per_dc[MAX_CONC][MDE];

  dbl vis;         /* flowing-liquid viscosity */
  dbl dvis_dT[MDE];/* sensitivities of flowing-liquid viscosity */
  dbl sc;          /* inertial coefficient */
  dbl speed;       /* magnitude of the velocity vector */

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  int v_g[DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int mode;


  /* Variables used for the modified fluid momentum equations when looking
   * at the particle momentum model.  Note that pd->MomentumFluxModel ==
   * SUSPENSION_PM when the Buyevich particle momentum equations are active.
   */
  int particle_momentum_on;	/* Boolean for particle momentum eq.'s active */
  int species=-1;		/* Species number of particle phase */
  double p_vol_frac=0;		/* Particle volume fraction (phi) */
  double ompvf=1;		/* 1 - p_vol_frac "One Minus Particle */
				/*    Volume Fraction" */
  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on =0;
  int source_on =0;
  int diffusion_on =0;
  int porous_brinkman_on =0; 
  int mesh_disp_on = pd->v[MESH_DISPLACEMENT1];
  
  dbl mass_etm, advection_etm, diffusion_etm, source_etm, porous_brinkman_etm;
  
  double *J = NULL;
  double *R ;
  
  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /* SUPG variables */
  dbl h_elem=0, h_elem_inv=0;
  dbl supg, d_wt_func;

  const double *hsquared = pg_data->hsquared ;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
					   centroid velocity for Q2 and the average
					   of the vertices for Q1. It comes from
					   the routine "element_velocity." */

  //Continuity stabilization
  dbl continuity_stabilization;
  dbl cont_gls;
  CONT_GLS_DEPENDENCE_STRUCT d_cont_gls_struct;
  CONT_GLS_DEPENDENCE_STRUCT *d_cont_gls = &d_cont_gls_struct;

  int *n_dof = NULL;
  int dof_map[MDE];
  
  
  status = 0;
  
  
  /*
   * Unpack variables from structures for local convenience...
   */

  eqn   = R_MOMENTUM1;

  /*
   * Bail out fast if there's nothing to do...
   */
#ifdef DEBUG_HKM
  bfm  = bf[R_MOMENTUM1];
  phi_i_vector = bfm->phi;
  checkFinite(phi_i_vector[0]);
#endif
  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;

  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  h3 = fv->h3;			/* Differential volume element (scales). */
  
  d_area = det_J * wt * h3;

  supg = 0.;

  if( mp->Mwt_funcModel == GALERKIN)
    { supg = 0.; }
  else if( mp->Mwt_funcModel == SUPG)
    { supg = mp->Mwt_func; }

  if (supg!=0.)
    {
      h_elem = 0.;
      for ( p=0; p<dim; p++)
        {
          h_elem += vcent[p] * vcent[p] * hsquared[p];
        }
      h_elem = sqrt(h_elem)/2.;
      if(h_elem == 0.)
        {
          h_elem_inv=0.;
        }
      else
        {
          h_elem_inv=1./h_elem;
        }

    }
  /* end Petrov-Galerkin addition */
  
  if( pd->v[POLYMER_STRESS11] )
    {
      (void)stress_eqn_pointer(v_s);

      v_g[0][0] = VELOCITY_GRADIENT11;
      v_g[0][1] = VELOCITY_GRADIENT12;
      v_g[1][0] = VELOCITY_GRADIENT21;
      v_g[1][1] = VELOCITY_GRADIENT22;
      v_g[0][2] = VELOCITY_GRADIENT13;
      v_g[1][2] = VELOCITY_GRADIENT23;
      v_g[2][0] = VELOCITY_GRADIENT31;
      v_g[2][1] = VELOCITY_GRADIENT32;
      v_g[2][2] = VELOCITY_GRADIENT33;
    }

  /* Set up variables for particle/fluid momentum coupling.
   */
  if(pd->e[R_PMOMENTUM1])
    {
      particle_momentum_on = 1;
      /* This is the species number of the particle phase. */
      species = (int) mp->u_density[0];
      p_vol_frac = fv->c[species];
      ompvf = 1.0 - p_vol_frac;
      /* Uncomment this to check for when the particle volume fraction
       * becomes non-physical.  Beware, however, that the intermediate
       * solutions may, indeed, become negative while converging to a
       * physical solution.
       if(p_vol_frac<0.0 || p_vol_frac>1.0)
       {
       if(fabs(p_vol_frac)<1e-14)
       p_vol_frac=0.0;
       else
       {
       printf("assemble_momentum: p_vol_frac=%g, exiting\n",p_vol_frac);
       exit(0);
       }
       }
      */
    }
  else
    particle_momentum_on = 0;

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */

  /*** Density ***/

  rho = density(d_rho, time);

  if ( pd->e[eqn] & T_POROUS_BRINK )
    {
      /* Load up remaining parameters for the Brinkman Equation. */
      por = mp->porosity;
      vis = mp->FlowingLiquid_viscosity;
      /* Load variable FlowingLiquid_viscosity */
      sc  = mp->Inertia_coefficient;

      /*
       * Special lubrication case
       */
      if ( mp->FSIModel > 0 ) {

	/* Store proper Gauss weights */
	wt = fv->wt;
	h3 = fv->h3;

	/* Prepare geometry */
	n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
	lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

	/* Adjust FEM weights */
	dim = pd->Num_Dim;
	fv->wt = wt;		
	det_J = fv->sdet;   

	/* 
	 * Calculate the velocity from the velocity calculator, of course!
	 * But then why do I even need this equation?
	 * I DON'T!!!
	 */
	calculate_lub_q_v(R_LUBP, time, dt, xi, exo);
  
	fv->wt = wt; /*load_neighbor_var_data screws fv->wt up */
      }
      


      if (mp->PorousMediaType != POROUS_BRINKMAN)
	WH(-1, "Set Porous term multiplier in continuous medium");

      if(mp->FlowingLiquidViscosityModel == CONSTANT)
	{
	  /* Do nothing */
	}
      else if (mp->FlowingLiquidViscosityModel == MOLTEN_GLASS)
	{
	  (void) molten_glass_viscosity(&(mp->FlowingLiquid_viscosity),
					dvis_dT, mp->u_FlowingLiquid_viscosity);
	}
      else if (mp->FlowingLiquidViscosityModel == USER)
	{
	  (void) usr_FlowingLiquidViscosity(mp->u_FlowingLiquid_viscosity);
	  var = TEMPERATURE;

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      dvis_dT[j]= mp->d_FlowingLiquid_viscosity[var]*bf[var]->phi[j];
	    }
	}
      else
	{
	  EH(-1,"Don't recognize your FlowingLiquidViscosity model");
	}

      if (mp->PermeabilityModel == SOLIDIFICATION)
	{
	  /* This is a permeability model that slows down
	   * the flow with increasing solid fraction.
	   * It should be useful for solidification problems.
	   */
	  per = solidification_permeability(h_elem_avg,d_per_dc);
	}
      else if (mp->PermeabilityModel == CONSTANT)
	{
	  per = mp->permeability;
	}
      else
	{
	  EH(-1,"Unrecognizable Permeability model");
	}
    }
  else
    {
      por = 1.;
      per = 1.;
      vis = mp->viscosity;
      sc  = 0.;
    }

  /*   eqn = R_MOMENTUM1; */
  /*
   * Field variables...
   */

  if (  transient_run &&  mesh_disp_on )
    x_dot = fv_dot->x;
  else
    x_dot = zero;
  if ( transient_run )
    v_dot = fv_dot->v;
  else
    v_dot = zero;


  /* for porous media stuff */
#ifdef DO_NO_UNROLL
  speed = 0.0;
  for ( a=0; a<wim; a++)
    {
      speed += v[a]*v[a];
    }
  speed = sqrt(speed);

  for ( a=0; a<VIM; a++) grad_v[a] = fv->grad_v[a];
#else
  speed = 0.0;
  speed += v[0]*v[0] + v[1]*v[1] ;
  if( wim == 3 ) speed += v[2]*v[2];
  
  speed = sqrt(speed);
  
  grad_v[0] = fv->grad_v[0];
  grad_v[1] = fv->grad_v[1];
  if( VIM == 3 ) grad_v[2] = fv->grad_v[2];
#endif
  
#ifdef DEBUG_HKM
  bfm  = bf[R_MOMENTUM1];
  phi_i_vector = bfm->phi;
  checkFinite(phi_i_vector[0]);
#endif

  /*
   * Calculate the momentum stress tensor at the current gauss point
   */
  fluid_stress( Pi, d_Pi );

  (void) momentum_source_term(f, df, time);

  //Call continuity stabilization if desired
  if(Cont_GLS)
    {
      calc_cont_gls(&cont_gls, d_cont_gls, time, pg_data);
    }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual)
    {
      /*
       * Assemble each component "a" of the momentum equation...
       */
      for ( a=0; a<wim; a++)
	{
	  eqn  = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm  = bf[eqn];
	  
	  
	  mass_on = pde[eqn] & T_MASS;
	  advection_on =  pde[eqn] & T_ADVECTION ;
	  diffusion_on = pde[eqn] & T_DIFFUSION;
	  source_on = pde[eqn] & T_SOURCE;
	  porous_brinkman_on = pde[eqn] & T_POROUS_BRINK;
	  
	  mass_etm = pd->etm[eqn][(LOG2_MASS)];
	  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
	  diffusion_etm = pd->etm[eqn][(LOG2_DIFFUSION)];
	  porous_brinkman_etm = pd->etm[eqn][(LOG2_POROUS_BRINK)];
	  source_etm = pd->etm[eqn][(LOG2_SOURCE)];

	  /*
	   * In the element, there will be contributions to this many equations
	   * based on the number of degrees of freedom...
	   */
	  
	  R = lec->R[peqn];
	  phi_i_vector = bfm->phi;
#ifdef DEBUG_HKM
	  checkFinite(phi_i_vector[0]);
#endif

	  for (i = 0; i < ei->dof[eqn]; i++) {
	    ledof = ei->lvdof_to_ledof[eqn][i];
	    if (ei->active_interp_ledof[ledof]) {
	      /*
	       *  Here is where we figure out whether the row is to placed in
	       *  the normal spot (e.g., ii = i), or whether a boundary condition
	       *  require that the volumetric contribution be stuck in another
	       *  ldof pertaining to the same variable type.
	       */
	      ii = ei->lvdof_to_row_lvdof[eqn][i];

		  
	      phi_i = phi_i_vector[i];
	      grad_phi_i_e_a = bfm->grad_phi_e[i][a];
		  

              /* this is an optimization for xfem */
	      if ( xfem != NULL )
                {
		  int xfem_active, extended_dof, base_interp, base_dof;
		  xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                                  &xfem_active, &extended_dof, &base_interp, &base_dof );
		  if ( extended_dof && !xfem_active ) continue;
		}

		  
	      mass = 0.;
	      if ( transient_run )
		{
		  if ( mass_on )
		    {
		      mass = v_dot[a] * rho;
		      mass *= - phi_i*d_area;
		      mass *= mass_etm;
		    }
		      
		  if ( porous_brinkman_on )
		    {
		      mass /= por; 
		    }
		      
		  if (particle_momentum_on)
		    {
		      mass *= ompvf;
		    }
		}

	      /* only use Petrov Galerkin on advective term - if required */
	      wt_func = phi_i;
	      /* add Petrov-Galerkin terms as necessary */
	      if (supg!=0.)
		{
		  for (p=0; p<dim; p++)
		    {
		      wt_func += supg * h_elem * v[p] * bfm->grad_phi[i][p];
		    }
		}
		  
	      advection = 0.;
	      if (advection_on)
		{
#ifdef DO_NO_UNROLL
		  for ( p=0; p<wim; p++)
		    {
		      advection += (v[p] - x_dot[p]) * grad_v[p][a];
		    }
#else
		  advection += (v[0] - x_dot[0]) * grad_v[0][a];
		  advection += (v[1] - x_dot[1]) * grad_v[1][a];
		  if (wim == 3) advection += (v[2] - x_dot[2]) * grad_v[2][a];
#endif  
		  advection *= rho;
		  advection *= - wt_func*d_area;
		  advection *= advection_etm;
		      
		  if (porous_brinkman_on)
		    {
		      por2 = por*por; 
		      advection /= por2; 
		    } 

		  if (particle_momentum_on)
		    {
		      advection *= ompvf;
		    }
		}
		  
	      porous = 0.;
	      if ( porous_brinkman_on )
		{
		  if(vis != 0.)
		    {
		      porous = v[a]*(rho*sc*speed/sqrt(per)+vis/per);
		      porous    *= -phi_i*d_area;
		      porous    *= porous_brinkman_etm;
		    }
		  else if ( mp->FSIModel > 0 )
		    {	      
		      porous = v[a] - LubAux->v_avg[a];
		      porous *= phi_i * wt * fv->sdet * h3;
		      porous *= porous_brinkman_etm;
		    }
		  else if (vis == 0. && mp->viscosity == 0.)
		    {
		      EH(-1, "cannot have both flowing liquid viscosity and mp->viscosity equal to zero");
		    }
		}
		  
	      diffusion = 0.;
	      if ( diffusion_on )
		{
#ifdef DO_NO_UNROLL
		  for ( p=0; p<VIM; p++) {
		    for ( q=0; q<VIM; q++) {
		      diffusion += grad_phi_i_e_a[p][q] * Pi[q][p]; 
		    }
		  }
#else
		  diffusion += grad_phi_i_e_a[0][0] * Pi[0][0]; 
		  diffusion += grad_phi_i_e_a[1][1] * Pi[1][1]; 
		  diffusion += grad_phi_i_e_a[1][0] * Pi[0][1]; 
		  diffusion += grad_phi_i_e_a[0][1] * Pi[1][0];
		  if( VIM == 3 ) {
		    diffusion += grad_phi_i_e_a[2][2] * Pi[2][2]; 
		    diffusion += grad_phi_i_e_a[2][1] * Pi[1][2]; 
		    diffusion += grad_phi_i_e_a[2][0] * Pi[0][2]; 
		    diffusion += grad_phi_i_e_a[1][2] * Pi[2][1]; 
		    diffusion += grad_phi_i_e_a[0][2] * Pi[2][0]; 
		  }  
#endif
		      
		  diffusion *= - d_area;
		  diffusion *= diffusion_etm;
		}
		  
	      /*
	       * Source term...
	       */
		  
	      source = 0.0;
	      if (source_on)
		{
		  source += f[a];		      
		  source *= phi_i * d_area;
		  source *= source_etm;
		}
		  
	      /* MMH For massful particles. */
	      if(Particle_Dynamics &&
		 (Particle_Model == SWIMMER_EXPLICIT || Particle_Model == SWIMMER_IMPLICIT))
		{
		  if(a == pd->Num_Dim - 1)
		    /* These data values should hold the entire
		     * source term. */
		    source = element_particle_info[ei->ielem].source_term[i];
		}
	      
	      //Continuity residual
	      continuity_stabilization = 0.0;
              if (Cont_GLS)
		{
		  for (p=0; p<VIM; p++){
		    continuity_stabilization += grad_phi_i_e_a[p][p];
		    }
                  continuity_stabilization *= cont_gls*d_area;
	        }

	      /*
	       * Add contributions to this residual (globally into Resid, and 
	       * locally into an accumulator)
	       */
		  
	      /*lec->R[peqn][ii] += mass + advection + porous + diffusion + source;*/
              R[ii] += mass + advection + porous + diffusion + source + continuity_stabilization;

#ifdef DEBUG_MOMENTUM_RES
	      printf("R_m[%d][%d] += %10f %10f %10f %10f %10f\n",
		     a,i,mass,advection,porous,diffusion,source);
#endif /* DEBUG_MOMENTUM_RES */
		  
	    }  /*end if (active_dofs) */
	  } /* end of for (i=0,ei->dofs...) */
	}
    }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian)
    {
      for (a = 0; a < wim; a++)
	{
	  eqn  = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm  = bf[eqn];
	  
	  mass_on = pde[eqn] & T_MASS;
	  advection_on =  pde[eqn] & T_ADVECTION ;
	  diffusion_on = pde[eqn] & T_DIFFUSION;
	  source_on = pde[eqn] & T_SOURCE;
	  porous_brinkman_on = pde[eqn] & T_POROUS_BRINK;
	  
	  mass_etm = pd->etm[eqn][(LOG2_MASS)];
	  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
	  diffusion_etm = pd->etm[eqn][(LOG2_DIFFUSION)];
	  porous_brinkman_etm = pd->etm[eqn][(LOG2_POROUS_BRINK)];
	  source_etm = pd->etm[eqn][(LOG2_SOURCE)];
	  
	  phi_i_vector = bfm->phi;
#ifdef DEBUG_HKM
	  checkFinite(phi_i_vector[0]);
#endif
	  
	  for (i = 0; i < ei->dof[eqn]; i++) {		  
	    ii = ei->lvdof_to_row_lvdof[eqn][i];

	    ledof = ei->lvdof_to_ledof[eqn][i];
	    if (ei->active_interp_ledof[ledof]) {
	      /*
	       *  Here is where we figure out whether the row is to placed in
	       *  the normal spot (e.g., ii = i), or whether a boundary condition
	       *  require that the volumetric contribution be stuck in another
	       *  ldof pertaining to the same variable type.
	       */
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
		  
	      phi_i = phi_i_vector[i];
		  
	      /* Assign pointers into the bf structure */
		  
	      grad_phi_i_e_a = bfm->grad_phi_e[i][a];	      
                 
	      /* this is an optimization for xfem */
	      if ( xfem != NULL )
		{
		  int xfem_active, extended_dof, base_interp, base_dof;
		  xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
				  &xfem_active, &extended_dof, &base_interp, &base_dof );
		  if ( extended_dof && !xfem_active ) continue;
		}
                    
	      d_grad_phi_i_e_a_dmesh = bfm->d_grad_phi_e_dmesh[i][a];
		  
#ifdef DEBUG_MOMENTUM_JAC		  
	      if(ei->ielem == 0)
		{
		  printf("\nASSEMBLE_MOMENTUM, a = %d, dof = %d\n", a, i);
		  printf("\tphi_i = %g\n", phi_i);
		}
#endif /* DEBUG_MOMENTUM_JAC */
		  
	      wt_func = phi_i;
	      /* add Petrov-Galerkin terms as necessary */
	      if(supg!=0.)
		{porous += (v[a]) * fv->dsurfdet_dx[b][j] * h3;
			      porous += (v[a]) * fv->sdet  * dh3dmesh_bj;

		  for(p=0; p<dim; p++)
		    {
		      wt_func += supg * h_elem * v[p] * bfm->grad_phi[i][p];
		    }
		}
		  
	      /* 
	       * J_m_T
	       */
	      if ( pdv[TEMPERATURE] )
		{
		  var = TEMPERATURE;
		  pvar = upd->vp[var];
		  phi_j_vector = bf[var]->phi;
		  J = lec->J[peqn][pvar][ii];
		      
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = phi_j_vector[j];	      
			  
		      mass = 0.;
			  
		      if (transient_run )
			{
			  if ( mass_on)
			    {
			      mass = d_rho->T[j] * v_dot[a];
			      mass *= - phi_i * d_area;
			      mass *= mass_etm;
			    }
			      
			  /* if porous flow is considered. KSC on 5/10/95 */ 
			  if (porous_brinkman_on ) 
			    {
			      mass /= por; 
			    } 
			}
			  
		      advection = 0.;
		      if ( advection_on )
			{
#ifdef DO_NO_UNROLL
			  for ( p=0; p<wim; p++)
			    {
			      advection += (v[p] - x_dot[p]) * grad_v[p][a];
			    }
#else
			  advection += (v[0] - x_dot[0]) * grad_v[0][a];
			  advection += (v[1] - x_dot[1]) * grad_v[1][a];
			  if(wim==3) advection += (v[2] - x_dot[2]) * grad_v[2][a];
#endif
				  
			  advection *= - wt_func * d_rho->T[j] * d_area;
			  advection *= advection_etm;
			      
			  if (porous_brinkman_on) 
			    {
			      por2 = por*por; 
			      advection /= por2; 
			    }
			}
			  
		      porous = 0.;
		      if (porous_brinkman_on)
			{
			  if (vis != 0.)
			    {
			      porous = v[a]*(d_rho->T[j] *sc*speed/sqrt(per));
			      porous += v[a]*(dvis_dT[j]/per);
			      porous    *= - phi_i * d_area;
			      porous    *= porous_brinkman_etm;
			    }
			  else
			    {
			      porous = 0.;
			    }
			}
		      
		      diffusion = 0.;
		      if ( diffusion_on )
			{
#ifdef DO_NO_UNROLL
			  for ( p=0; p<VIM; p++) {
			    for ( q=0; q<VIM; q++) {
			      diffusion +=  grad_phi_i_e_a[p][q] * d_Pi->T[q][p][j];
			    }
			  }
#else
			  diffusion +=  grad_phi_i_e_a[0][0] * d_Pi->T[0][0][j];
			  diffusion +=  grad_phi_i_e_a[1][1] * d_Pi->T[1][1][j];
			  diffusion +=  grad_phi_i_e_a[1][0] * d_Pi->T[0][1][j];
			  diffusion +=  grad_phi_i_e_a[0][1] * d_Pi->T[1][0][j];
			  if(VIM == 3)
			    {
			      diffusion +=  grad_phi_i_e_a[2][0] * d_Pi->T[0][2][j];
			      diffusion +=  grad_phi_i_e_a[2][1] * d_Pi->T[1][2][j];
			      diffusion +=  grad_phi_i_e_a[2][2] * d_Pi->T[2][2][j];
			      diffusion +=  grad_phi_i_e_a[0][2] * d_Pi->T[2][0][j];
			      diffusion +=  grad_phi_i_e_a[1][2] * d_Pi->T[2][1][j];
			    }
#endif
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;
			}
			  
		      source    = 0.;
		      if ( source_on )
			{
			  source = phi_i * df->T[a][j] * d_area;
			  source *= source_etm;
			}
			  
		      if (particle_momentum_on)
			{
			  mass *= (1.0 - p_vol_frac);
			  advection *= (1.0 - p_vol_frac);
			}
			  
		      /*lec->J[peqn][pvar][ii][j] += mass + advection + porous + diffusion + source;*/
		      J[j] += mass + advection + porous + diffusion + source; 
#ifdef DEBUG_HKM
		      checkFinite(J[j]);
#endif
		    }
		}


	      /*
	       * J_m_nn
	       */
		  
		  
	      if ( pdv[BOND_EVOLUTION] )
		{
		  var = BOND_EVOLUTION;
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];	      
			  
		      diffusion = 0.;
		      if ( diffusion_on)
			{
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  diffusion +=  grad_phi_i_e_a[p][q] * 
				    d_Pi->nn[q][p][j];
				}
			    }
			  diffusion *= - d_area ;
			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}
			  
		      lec->J[peqn][pvar][ii][j] +=
			diffusion;
		    }
		}

#ifdef COUPLED_FILL
	      /* 
	       * J_m_F
	       */
		  
		  
	      if ( pdv[FILL] )
		{
		  var = FILL;
		  pvar = upd->vp[var];
			  
		  J = lec->J[peqn][pvar][ii];
			  
		  phi_j_vector = bf[var]->phi;
			  
		  for ( j=0; j < ei->dof[var]; j++)
		    {
		      phi_j = phi_j_vector[j];	      

				  
		      mass = 0.;
		      if ( transient_run )
			{
			  if ( mass_on )
			    {
			      mass = d_rho->F[j] * v_dot[a];
			      mass *= - phi_i * d_area;
			      mass *= mass_etm;
			    }
			      
			  /* if porous flow is considered. KSC on 5/10/95 */ 
			  if (porous_brinkman_on) 
			    {
			      if (vis != 0.)
				{
				  mass /= por; 
				}
			      else if (mp->viscosity ==0 )
				{
				  EH(-1, "Error in Porous Brinkman specs.  See PRS ");
				}
			    } 
			}
			  
		      advection = 0.;
		      if ( advection_on )
			{
#ifdef DO_NO_UNROLL
			  for ( p=0; p<wim; p++)
			    {
			      advection += (v[p] - x_dot[p]) * grad_v[p][a];
			    }
#else 
			  advection += (v[0] - x_dot[0]) * grad_v[0][a];
			  advection += (v[1] - x_dot[1]) * grad_v[1][a];
			  if(wim==3) advection += (v[2] - x_dot[2]) * grad_v[2][a];
#endif
			  advection *= - wt_func * d_rho->F[j] * d_area;
			  advection *= advection_etm;
			      
			  /* if porous flow is considered. KSC on 5/10/95 */ 
			  if ( porous_brinkman_on) 
			    {

			      if (vis != 0.)
				{
				  por2	     = por * por; 
				  advection /= por2;  
				}

			      
			    } 
			}
			  
		      porous = 0.;
		      if ( porous_brinkman_on )
			{
			  if (vis != 0)
			    {
			      porous = v[a]*(d_rho->F[j] *sc*speed/sqrt(per));
			      /* porous += v[a]*(dvis_dF[j]/per); */
			      porous    *= - phi_i * d_area;
			      porous    *= porous_brinkman_etm;
			    }
			  else if ( mp->FSIModel > 0 )
			    {
			      porous = 0.0;
			      porous -= LubAux->dv_avg_df[a][j] * phi_i;
			      porous *= wt * fv->sdet * h3;
			      porous *= porous_brinkman_etm;
			    }
			  else
			    {
			      porous = 0.;
			    }
			}
			  
		      diffusion = 0.;
		      if ( diffusion_on )
			{
#ifdef DO_NO_UNROLL
			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  diffusion +=  grad_phi_i_e_a[p][q] *d_Pi->F[q][p][j];
				}
			    }
#else
			  diffusion +=  grad_phi_i_e_a[0][0] * d_Pi->F[0][0][j];
			  diffusion +=  grad_phi_i_e_a[1][1] * d_Pi->F[1][1][j];
			  diffusion +=  grad_phi_i_e_a[1][0] * d_Pi->F[0][1][j];
			  diffusion +=  grad_phi_i_e_a[0][1] * d_Pi->F[1][0][j];
				  
			  if(VIM == 3)
			    {
			      diffusion +=  grad_phi_i_e_a[2][0] * d_Pi->F[0][2][j];
			      diffusion +=  grad_phi_i_e_a[2][1] * d_Pi->F[1][2][j];
			      diffusion +=  grad_phi_i_e_a[2][2] * d_Pi->F[2][2][j];
			      diffusion +=  grad_phi_i_e_a[0][2] * d_Pi->F[2][0][j];
			      diffusion +=  grad_phi_i_e_a[1][2] * d_Pi->F[2][1][j];
			    }
#endif				
			  diffusion *= - d_area;
			  diffusion *= diffusion_etm;
			}
			  
		      source = 0.;
			  
		      /* Stay away from evil Hessians. */
		      if ( source_on )
			{
			  source = phi_i * df->F[a][j] * d_area;
			  source *= source_etm;
			}
			  
		      if (particle_momentum_on)
			{
			  mass *= (1.0 - p_vol_frac);
			  advection *= (1.0 - p_vol_frac);
			}
		      J[j] +=  mass + advection + porous + diffusion + source;
		      /*lec->J[peqn][pvar][ii][j] +=mass + advection + porous + diffusion + source; */
		    }
		}
#endif /* COUPLED_FILL */

	  
	      /*
	       * J_m_v
	       */
	      for ( b=0; b<wim; b++)
		{
		  var = VELOCITY1+b;
		  if ( pdv[var] )
		    {
		      pvar = upd->vp[var];
			  
		      J = lec->J[peqn][pvar][ii];
			  
		      phi_j_vector = bf[var]->phi;
			  
		      for ( j=0; j<ei->dof[var]; j++)
			{
			      
			  phi_j = phi_j_vector[j];
				  			      
			  mass = 0.;
			  if ( transient_run)
			    {
			      if ( mass_on && (a == b ))
				{
				  /*mass = (1.+2.*tt) * phi_j/dt * (double)delta(a,b); */
				  mass = (1.+2.*tt) * phi_j/dt ;
				  mass *= - phi_i * rho * d_area;
				  mass *= mass_etm;
				}
				  
			      if (porous_brinkman_on) 
				{
				  if(vis != 0.)
				    {
				      mass /= por; 
				    }
				  else if (mp->viscosity == 0)
				    {
				      EH(-1, "incorrect viscosity settings on porous_brinkman");
				    }
				} 

			      if (particle_momentum_on)
				{
				  mass *= ompvf;
				}
			    }
			      
			      
			  porous = 0.0;
			  if ( porous_brinkman_on )
			    {
			      if(vis != 0)
				{
				  porous =   ((rho*sc/sqrt(per))*(2.*v[b])*v[a] +
					      (rho*sc*speed/sqrt(per) + vis/per)*(double)delta(a,b));
				  porous *= -phi_i*phi_j*d_area;
				  porous *= porous_brinkman_etm;
				}
			      else if (mp->viscosity !=0)
				{
				  porous = delta(a,b);
				  porous *= phi_i*phi_j* wt * fv->sdet * h3;
				  // porous *= phi_i*phi_j* wt * h3;
				  porous *= porous_brinkman_etm;
				}
			      else
				{
				  EH(-1, "incorrect viscosity settings on porous_brinkman");
				}
			    }


			  advection = 0.;
			  advection_a = 0.;
			  if ( advection_on )
			    {
			      advection_a += phi_j * grad_v[b][a];

#ifdef DO_NO_UNROLL
			      for ( p=0; p<wim; p++)
				{
				  advection_a += (v[p] - x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
				}
#else 
			      advection_a += (v[0] - x_dot[0]) * bf[var]->grad_phi_e[j][b][0][a];
			      advection_a += (v[1] - x_dot[1]) * bf[var]->grad_phi_e[j][b][1][a];
			      if(wim==3) advection_a += (v[2] - x_dot[2]) * bf[var]->grad_phi_e[j][b][2][a];
#endif
				  

			      advection_a *= rho;
			      advection_b = 0.;
			      if(supg!=0.)
                            	{
				  d_wt_func = supg * h_elem * phi_j*bfm->grad_phi[i][b];

				  for(p=0;p<dim;p++)
				    {
				      d_wt_func += supg * vcent[b] * 
                                           pg_data->dv_dnode[b][j] *
					hsquared[b] * h_elem_inv / 4. *
					v[p] * bfm->grad_phi[i][p];

				    }
				  advection_b +=  advection_a;

				  advection_b *=  d_wt_func;
				  advection_b *= - det_J * wt;
				  advection_b *= h3;
                                }
			      advection_a *= -wt_func * d_area;
			      advection = advection_a + advection_b;
			      advection *= advection_etm;


			      /* if porous flow is considered. KSC on 5/10/95 */
			      if ( porous_brinkman_on)
				{
				  if(vis !=0)
				    {
				      por2 = por*por;
				      advection /= por2;
				    }
				}

			      if (particle_momentum_on)
				{
				  advection *= ompvf;
				}
			    }

			  diffusion = 0.;
			  if ( diffusion_on)
			    {
#ifdef DO_NO_UNROLL
			      for ( p=0; p<VIM; p++)
				{
				  for ( q=0; q<VIM; q++)
				    {
				      diffusion += grad_phi_i_e_a[p][q] * d_Pi->v[q][p][b][j];
				    }
				}
#else
			      diffusion += grad_phi_i_e_a[0][0] * d_Pi->v[0][0][b][j];
			      diffusion += grad_phi_i_e_a[1][1] * d_Pi->v[1][1][b][j];
			      diffusion += grad_phi_i_e_a[1][0] * d_Pi->v[0][1][b][j];
			      diffusion += grad_phi_i_e_a[0][1] * d_Pi->v[1][0][b][j];
					  
			      if(VIM==3)
				{
				  diffusion += grad_phi_i_e_a[2][0] * d_Pi->v[0][2][b][j];
				  diffusion += grad_phi_i_e_a[2][1] * d_Pi->v[1][2][b][j];
				  diffusion += grad_phi_i_e_a[0][2] * d_Pi->v[2][0][b][j];
				  diffusion += grad_phi_i_e_a[1][2] * d_Pi->v[2][1][b][j];
				  diffusion += grad_phi_i_e_a[2][2] * d_Pi->v[2][2][b][j];
				}
#endif	  
					  
			      diffusion *= -d_area;
			      diffusion *= diffusion_etm;
			    }

			  source    = 0.;
			  if ( source_on )
			    {
			      source    = phi_i * df->v[a][b][j] * d_area;
			      source   *= source_etm;
			    }

			  continuity_stabilization = 0.0;
			  if (Cont_GLS)
			    {
			      for (p=0; p<VIM; p++)
				{
				  continuity_stabilization += grad_phi_i_e_a[p][p];
				}
			      continuity_stabilization *= d_cont_gls->v[b][j]*d_area;
			    }


			  J[j] +=  mass + advection + porous + diffusion + source + continuity_stabilization;

			  /*lec->J[peqn][pvar][ii][j] +=  mass + advection + porous + diffusion + source; */
			}
		    }
		}

	      /*
	       * J_m_vd
	       */
	      if( pdv[VORT_DIR1] ) {
		for ( b=0; b<DIM; b++)
		  {
		    var = VORT_DIR1+b;
		    if ( pdv[var] )
		      {
			pvar = upd->vp[var];
			for ( j=0; j<ei->dof[var]; j++)
			  {
						  
			    phi_j = bf[var]->phi[j];
						  
			    diffusion = 0.;
			    if ( pd->e[eqn] & T_DIFFUSION )
			      {
							  
				for ( p=0; p<VIM; p++)
				  {
				    for ( q=0; q<VIM; q++)
				      {
					diffusion += grad_phi_i_e_a[p][q] *
					  d_Pi->vd[q][p][b][j];
				      }
				  }
				diffusion *= -det_J * wt;
				diffusion *= h3;
				diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			      }
						  
			    lec->J[peqn][pvar][ii][j] +=
			      diffusion;
			  }
		      }
		  }
	      }

	      /*
	       * J_m_E (E_field from EHD source)
	       */
	      if ( pdv[EFIELD1] ) {
		for ( b=0; b<wim; b++)
		  {
		    var = EFIELD1+b;
		    if ( pdv[var] )
		      {
			pvar = upd->vp[var];
			for ( j=0; j<ei->dof[var]; j++)
			  {
						  
			    phi_j = bf[var]->phi[j];
						  
						  
			    source    = 0.;
			    if ( pd->e[eqn] & T_SOURCE )
			      {
				source    = phi_i * df->E[a][b][j] * det_J * h3 *wt;
				source   *= pd->etm[eqn][(LOG2_SOURCE)];
			      }
						  
			    lec->J[peqn][pvar][ii][j] += source;
			  }
		      }
		  }
	      }

	      /*
	       * J_m_lubp (lubrication pressure field for porous-brinkman)
	       */
	      if ( pdv[LUBP] ) {
		var = LUBP;
		if ( pdv[var] )
		  {
		    pvar = upd->vp[var];
		    for ( j=0; j<ei->dof[var]; j++)
		      {
						  
			porous    = 0.;
			if ( porous_brinkman_on)
			  {
			    /*This needs no insulation as the code shouldn't get here 
			     *for traditional Brinkman form */
             
			    /* Need a few more basis functions */
			    dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
			    ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
			    
			    /* Assemble */
			    porous -= LubAux->dv_avg_dp1[a][j] * grad_II_phi_j[a];
			    porous -= LubAux->dv_avg_dp2[a][j] * phi_j;
			    porous *= phi_i * fv->sdet * wt * h3;
			    // porous *= phi_i * wt * h3;
			    porous *= porous_brinkman_etm;
			  }
						  
			lec->J[peqn][pvar][ii][j] += porous;
		      }
		  }
	      }

	      /*
	       * J_m_d_delta_h (Delta_H equation for shell for porous-brinkman)
	       */
	      if ( pdv[SHELL_DELTAH] ) {
		var = SHELL_DELTAH;
		if ( pdv[var] )
		  {
                   
		    pvar = upd->vp[var];
		    for ( j=0; j<ei->dof[var]; j++)
		      {
						  
			phi_j = bf[var]->phi[j];
						  
						  
			porous    = 0.;
			if ( porous_brinkman_on)
			  {
			    porous -= LubAux->dv_avg_ddh[a][j]*phi_j;
			    porous *= phi_i * fv->sdet * wt * h3;
			    porous *= porous_brinkman_etm;
			  }
						  
			lec->J[peqn][pvar][ii][j] += porous;
		      }
		  }
	      }

	      /*
	       * J_m_d_rs (realsolid equation for shell)
	       */
	      if (upd->vp[SOLID_DISPLACEMENT1]) {
		    for ( b=0; b<dim; b++)
		      {
			var = SOLID_DISPLACEMENT1+b;
			if ( pdv[var] )
			  {	    
			    if ( mp->FSIModel == FSI_REALSOLID_CONTINUUM ) {
			      for ( j = 0; j < ei->dof[SOLID_DISPLACEMENT1]; j++ ) {
				jk=dof_map[j];
			    
				phi_j = bf[var]->phi[j];

				porous = 0.;
				if ( porous_brinkman_on ) {
				  porous -= LubAux->dv_avg_drs[a][b][j] * fv->sdet * h3;
				}
				porous *= phi_i * wt;


				J[jk] += porous;
    
			      }
			    }

			  }
		      }
		  }

	      /*
	       * J_m_d_curv (lubrication curvature equation for shell)
	       */
	      if ( pdv[SHELL_LUB_CURV] ) {
		var = SHELL_LUB_CURV;
		if ( pdv[var] )
		  {
                   
		    pvar = upd->vp[var];
		    for ( j=0; j<ei->dof[var]; j++)
		      {
						  
			phi_j = bf[var]->phi[j];
						  
						  
			porous    = 0.;
			if ( porous_brinkman_on)
			  {
			    porous -= LubAux->dv_avg_dk[a][j]*phi_j;
			    porous *= phi_i * fv->sdet * wt * h3;
			    porous *= porous_brinkman_etm;
			  }
						  
			lec->J[peqn][pvar][ii][j] += porous;
		      }
		  }
	      }

	      /*
	       * J_m_sh_fp (film pressure field for porous-brinkman)
	       */
	      if ( pdv[SHELL_FILMP] ) {
		var = SHELL_FILMP;
		if ( pdv[var] )
		  {
		    pvar = upd->vp[var];
		    for ( j=0; j<ei->dof[var]; j++)
		      {
						  
			porous    = 0.;
			if ( porous_brinkman_on)
			  {
			    /*This needs no insulation as the code shouldn't get here 
			     *for traditional Brinkman form */
             
			    /* Need a few more basis functions */
			    dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
			    ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);

			    /* Assemble */
			    porous -= LubAux->dv_avg_dp1[a][j] * grad_II_phi_j[a];
			    porous *= phi_i * fv->sdet * wt * h3;
			    porous *= porous_brinkman_etm;
			  }
						  
			lec->J[peqn][pvar][ii][j] += porous;
		      }
		  }
	      }

	      /*
	       * J_m_sh_fh (film thickness for porous-brinkman)
	       */
	      if ( pdv[SHELL_FILMH] ) {
		var = SHELL_FILMH;
		if ( pdv[var] )
		  {
		    pvar = upd->vp[var];
		    for ( j=0; j<ei->dof[var]; j++)
		      {
			porous    = 0.;
			if ( porous_brinkman_on)
			  {
			    /*This needs no insulation as the code shouldn't get here 
			     *for traditional Brinkman form */

			    dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
             	            ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh, n_dof[MESH_DISPLACEMENT1], dof_map);
			    
			    /* Assemble */
			    porous -= LubAux->dv_avg_dh1[a][j] * grad_II_phi_j[a];
			    porous -= LubAux->dv_avg_dh2[a][j] * phi_j;
			    porous *= phi_i * fv->sdet * wt * h3;
			    porous *= porous_brinkman_etm;
			  }
						  
			lec->J[peqn][pvar][ii][j] += porous;
		      }
		  }
	      }

	      /*
	       * J_m_c
	       */
		  
	      if ( pdv[MASS_FRACTION] )
		{
		  var = MASS_FRACTION;
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      for ( w=0; w<pd->Num_Species_Eqn; w++)
			{
			  mass = 0.;
			  if ( transient_run )
			    {
			      if ( pd->e[eqn] & T_MASS )
				{
				  /* For particle momentum stuff, d_rho->* = 0, since the density per phase
				   * is actually constant.  However, the mass term was multiplied by (1-phi),
				   * so that we get a J_m_c entry now.
				   */
				  if (particle_momentum_on && w == species)
				    mass += -rho * v_dot[a];
				  else
				    mass += d_rho->C[w][j] * v_dot[a];
				  mass *= - phi_i * det_J * wt* h3;
				  mass *= pd->etm[eqn][(LOG2_MASS)];
				}

			      if (pd->e[eqn] & T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
				{
				  mass /= por;
				}
			    }

			  advection = 0.;
			  if ( pd->e[eqn] & T_ADVECTION )
			    {
			      for ( p=0; p<wim; p++)
				{
				  advection += (v[p] - x_dot[p]) *
				    grad_v[p][a];
				}
			      if(particle_momentum_on &&
				 w == species)
				advection *= -phi_j * (-wt_func * rho * det_J * wt * h3);
			      else
				advection *= - wt_func * d_rho->C[w][j] * det_J * wt * h3;
			      advection *= pd->etm[eqn][(LOG2_ADVECTION)];

			      if (pd->e[eqn] & T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
				{
				  por2 = por*por;
				  advection /= por2;
				}
			    }

			  porous    = 0.;
			  if ( pd->e[eqn] & T_POROUS_BRINK )
			    {
			      porous = v[a]*(d_rho->C[w][j] *sc*speed/sqrt(per)
					     -0.5*rho*sc*speed/(per*sqrt(per))*d_per_dc[w][j]
					     -vis/(per*per)*d_per_dc[w][j]);
			      porous    *= - phi_i * det_J * wt * h3;
			      porous    *= pd->etm[eqn][(LOG2_POROUS_BRINK)];

			      /*Hack again here*/
			      porous = 0.;
			    }


			  diffusion = 0.;
			  if ( pd->e[eqn] & T_DIFFUSION )
			    {
			      for ( p=0; p<VIM; p++)
				{
				  for ( q=0; q<VIM; q++)
				    {
				      diffusion += grad_phi_i_e_a[p][q] *
					d_Pi->C[q][p][w][j];
				    }
				}
			      diffusion *= - det_J * wt;
			      diffusion *= h3;
			      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			    }

			  source = 0.;
			  if ( pd->e[eqn] & T_SOURCE )
			    {
			      /* df->C was calculated in mm_std_models.c */
			      source    = phi_i * df->C[a][w][j] * det_J * h3 *wt;
			      source   *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
		  

			  lec->J[peqn][MAX_PROB_VAR + w][ii][j] +=
			    mass + advection + porous + diffusion + source;
			}
		    }
		}

	      /*
	       * J_m_P
	       */
		  
	      if ( pdv[PRESSURE] ){
		var = PRESSURE;
		pvar = upd->vp[var];
			   
		J = lec->J[peqn][pvar][ii];
			   
		for ( j=0; j<ei->dof[var]; j++)
		  {

		    /*  phi_j = bf[var]->phi[j]; */
				   

				   
		    diffusion = 0.;
				   
		    if ( diffusion_on)
		      {
#ifdef DO_NO_UNROLL
			for ( p=0; p<VIM; p++)
			  {
			    for ( q=0; q<VIM; q++)
			      {
							   
				diffusion -=  grad_phi_i_e_a[p][q] * d_Pi->P[q][p][j];
							   
#ifdef DEBUG_MOMENTUM_JAC
				if(ei->ielem == 0)
				  printf("P diffusion, %d, %d  -= %g %g\n",
					 p, q,
					 grad_phi_i_e_a[p][q], d_Pi->P[j][q][p]);
#endif /* DEBUG_MOMENTUM_JAC */
			      }
			  }
#else
			diffusion -=  grad_phi_i_e_a[0][0] * d_Pi->P[0][0][j];
			diffusion -=  grad_phi_i_e_a[1][1] * d_Pi->P[1][1][j];
			diffusion -=  grad_phi_i_e_a[1][0] * d_Pi->P[0][1][j];
			diffusion -=  grad_phi_i_e_a[0][1] * d_Pi->P[1][0][j];
					   
			if(VIM==3)
			  {
			    diffusion -=  grad_phi_i_e_a[2][2] * d_Pi->P[2][2][j];
			    diffusion -=  grad_phi_i_e_a[2][1] * d_Pi->P[1][2][j];
			    diffusion -=  grad_phi_i_e_a[2][0] * d_Pi->P[0][2][j];
			    diffusion -=  grad_phi_i_e_a[0][2] * d_Pi->P[2][0][j];
			    diffusion -=  grad_phi_i_e_a[1][2] * d_Pi->P[2][1][j];
			  }
#endif
			diffusion *= d_area;
			diffusion *= diffusion_etm;
		      }
				   
				   
#ifdef DEBUG_MOMENTUM_JAC
		    if(ei->ielem == 0)
		      fprintf(stdout,
			      "\tJ_m_P[%d][%d][%d] += %10f %10f %10f %10f\n",
			      a,i,j,mass,advection,porous,diffusion);
#endif /* DEBUG_MOMENTUM_JAC */
				   
		    /*lec->J[peqn][pvar][ii][j] += diffusion ;  */
		    J[j] += diffusion;
		  }
	      }

	      /*
	       * J_m_S
	       */
	      if ( pdv[POLYMER_STRESS11]) {
		for ( mode=0; mode<vn->modes; mode++)
		  {
		    for ( b=0; b<VIM; b++)
		      {
			for ( c=0; c<VIM; c++)
			  {
			    var = v_s[mode][b][c];
						  
			    if ( pdv[var] )
			      {
							  
				pvar = upd->vp[var];
							  
				for ( j=0; j<ei->dof[var]; j++)
								  
				  {
				    phi_j = bf[var]->phi[j];
								  
				    diffusion = 0.;
								  
				    if ( pd->e[eqn] & T_DIFFUSION )
				      {
					for ( p=0; p<VIM; p++)
					  {
					    for ( q=0; q<VIM; q++)
					      {
						diffusion -=
						  grad_phi_i_e_a[p][q] * d_Pi->S[p][q][mode][b][c][j];
					      }
					  }
					diffusion *= det_J * wt * h3;
					diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
				      }
								  
				    lec->J[peqn][pvar][ii][j] +=
				      diffusion;
				  }
			      }
			  }
		      }
		  }
	      }
		  
	      /*
	       * J_m_G
	       */
		  
	      if ( pdv[POLYMER_STRESS11] && (vn->evssModel == EVSS_F) )
		{
		  for ( b=0; b<VIM; b++)
		    {
		      for ( c=0; c<VIM; c++)
			{
			  var = v_g[b][c];
					  
			  if ( pdv[var] )
			    {
						  
			      pvar = upd->vp[var];
						  
			      for ( j=0; j<ei->dof[var]; j++)
							  
				{
				  phi_j = bf[var]->phi[j];
							  
				  diffusion = 0.;
							  
				  if ( pd->e[eqn] & T_DIFFUSION )
				    {
				      for ( p=0; p<VIM; p++)
					{
					  for ( q=0; q<VIM; q++)
					    {
					      diffusion +=
						grad_phi_i_e_a[p][q] * d_Pi->g[q][p][b][c][j];
					    }
					}
				      diffusion *= det_J * wt *h3;
				      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
				    }
							  
				  lec->J[peqn][pvar][ii][j] +=diffusion;
				}
			    }
			}
		    }
		}

	      /*
	       * J_m_d
	       */
	      if (pdv[MESH_DISPLACEMENT1]) {
		for ( b=0; b<dim; b++)
		  {
		    var = MESH_DISPLACEMENT1+b;
		    if ( pdv[var] )
		      {
			pvar = upd->vp[var];
					  
			J = lec->J[peqn][pvar][ii];
					  
			for ( j=0; j<ei->dof[var]; j++)
			  {
			    phi_j = bf[var]->phi[j];
						  
			    d_det_J_dmesh_bj = bfm->d_det_J_dm[b][j];
						  
			    dh3dmesh_bj = fv->dh3dq[b] * phi_j;
						  
			    mass = 0.;
						  
			    if ( transient_run)
			      {
				if ( mass_on)
				  {
				    mass = v_dot[a];
				    mass *= - phi_i * rho *
				      ( d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj )
				      * wt;
				    mass *= mass_etm;
				  }
							  
				/* if porous flow is considered. KSC on 5/10/95 */ 
				if (porous_brinkman_on) 
				  {
				    mass /= por; 
				  } 
							  
			      }
						  
			    porous = 0.;
			    if (porous_brinkman_on )
			      {
				if(vis != 0.)
				  {
				    porous += v[a]*(rho*sc*speed/sqrt(per)+vis/per)*(double)delta(a,b);
				    porous *= -phi_i * wt * (  d_det_J_dmesh_bj * h3  + det_J * dh3dmesh_bj );
				    porous *= porous_brinkman_etm;
				    
				    /* hack again here */
				    porous = 0.;
				     
				  } 

				/* Case for sheet-only  */
				else if ( mp->FSIModel == FSI_SHELL_ONLY ) {
				  jk=dof_map[j];
				  porous += (v[a] - LubAux->v_avg[a]) * fv->dsurfdet_dx[b][jk] * h3;
				  porous += (v[a] - LubAux->v_avg[a]) * fv->sdet  * dh3dmesh_bj;
				  porous -= LubAux->dv_avg_dx[a][b][j] * fv->sdet * h3;
				  porous *= phi_i * wt;    
				}
				
				else if (vis == 0. && mp->viscosity == 0.)
				  {
				    EH(-1, "cannot have both flowing liquid viscosity and mp->viscosity equal to zero");
				  }

			      }
						  
						  
			    advection   = 0.;
						  
			    if ( advection_on )
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
#ifdef DO_NO_UNROLL
				advection_a = 0.;
				for ( p=0; p<wim; p++)
				  {
				    advection_a += 
				      (v[p]-x_dot[p]) 
				      * fv->d_grad_v_dmesh[p][a] [b][j];
				  }
				advection_a *= -wt_func * rho * d_area;
							  
				advection_b = 0.;		      
				for ( p=0; p<wim; p++)
				  {
				    advection_b += (v[p]    -x_dot[p]) * grad_v[p][a];
				  }
#else
				advection_a = 0.;
				advection_a += (v[0]-x_dot[0]) * fv->d_grad_v_dmesh[0][a] [b][j];
				advection_a += (v[1]-x_dot[1]) * fv->d_grad_v_dmesh[1][a] [b][j];
				if( wim == 3 ) advection_a += (v[2]-x_dot[2]) * fv->d_grad_v_dmesh[2][a] [b][j];
							  
				advection_a *= -wt_func * rho * d_area;
							  
				advection_b = 0.;		      
				advection_b += (v[0]    -x_dot[0]) * grad_v[0][a];				  
				advection_b += (v[1]    -x_dot[1]) * grad_v[1][a];				  
				if ( wim == 3 ) advection_b += (v[2]    -x_dot[2]) * grad_v[2][a];				  
#endif
				advection_b *= 
				  -wt_func * rho * wt * 
				  ( d_det_J_dmesh_bj * h3 
				    + det_J * dh3dmesh_bj );
							  
							  
				advection_c = 0.;	
				if ( transient_run )
				  {
				    if ( mass_on)
				      {
#ifdef DO_NO_UNROLL
					for ( p=0; p<wim; p++)
					  {
					    advection_c +=   (-(1.+2.*tt) * phi_j/dt * (double)delta(p,b)) 
					      * grad_v[p][a];
					  }
#else
					advection_c +=   (-(1.+2.*tt) * phi_j/dt * (double)delta(0,b)) * grad_v[0][a];						  
					advection_c +=   (-(1.+2.*tt) * phi_j/dt * (double)delta(1,b)) * grad_v[1][a];						  
					if( wim == 3 ) advection_c +=   (-(1.+2.*tt) * phi_j/dt * (double)delta(2,b)) * grad_v[2][a];						  
#endif
					advection_c *= -wt_func * rho * d_area;
				      }
				  }
							  
				advection = advection_a + advection_b + advection_c;
							  
				advection *= advection_etm;
				/* if porous flow is considered. KSC on 5/10/95 */ 
				if (porous_brinkman_on) 
				  {
				    por2 = por*por; 
				    advection /= por2; 
				  } 
			      }
						  
			    /*
			     * Diffusion...
			     */
						  
			    diffusion = 0.;
						  
			    if ( diffusion_on )
			      {
							  
				/* Three parts:
				 *   diff_a = 
				 *   Int ( d(grad(phi_i e_a))/dmesh : Pi h3 |Jv|)
				 *
				 *   diff_b = 
				 *   Int ( grad(phi_i e_a) : d(Pi)/dmesh h3 |Jv|)
				 *
				 *   diff_c = 
				 *   Int ( grad(phi_i e_a) : Pi d(h3|Jv|)/dmesh )
				 */
							  
				diff_a = 0.;
				diff_b = 0.;
				diff_c = 0.;
#ifdef DO_NO_UNROLL
				for ( p=0; p<VIM; p++)
				  {
				    for ( q=0; q<VIM; q++)
				      {
					diff_a += 
					  d_grad_phi_i_e_a_dmesh[p][q] [b][j] * Pi[q][p];
									  
					diff_b +=
					  grad_phi_i_e_a[p][q] * d_Pi->X[q][p][b][j];
									  
					diff_c += grad_phi_i_e_a[p][q] * Pi[q][p];
				      }
				  }
#else
				diff_a += d_grad_phi_i_e_a_dmesh[0][0] [b][j] * Pi[0][0];				  
				diff_a += d_grad_phi_i_e_a_dmesh[1][1] [b][j] * Pi[1][1];				  
				diff_a += d_grad_phi_i_e_a_dmesh[0][1] [b][j] * Pi[1][0];				  
				diff_a += d_grad_phi_i_e_a_dmesh[1][0] [b][j] * Pi[0][1];		
							  
				diff_b += grad_phi_i_e_a[0][0] * d_Pi->X[0][0][b][j];
				diff_b += grad_phi_i_e_a[1][1] * d_Pi->X[1][1][b][j];
				diff_b += grad_phi_i_e_a[0][1] * d_Pi->X[1][0][b][j];
				diff_b += grad_phi_i_e_a[1][0] * d_Pi->X[0][1][b][j];
							  
				diff_c += grad_phi_i_e_a[0][0] * Pi[0][0];					
				diff_c += grad_phi_i_e_a[1][1] * Pi[1][1];					
				diff_c += grad_phi_i_e_a[0][1] * Pi[1][0];					
				diff_c += grad_phi_i_e_a[1][0] * Pi[0][1];					
							  
				if( VIM == 3)
				  {
				    diff_a += d_grad_phi_i_e_a_dmesh[2][2] [b][j] * Pi[2][2];				  
				    diff_a += d_grad_phi_i_e_a_dmesh[2][1] [b][j] * Pi[1][2];				  
				    diff_a += d_grad_phi_i_e_a_dmesh[2][0] [b][j] * Pi[0][2];				  
				    diff_a += d_grad_phi_i_e_a_dmesh[1][2] [b][j] * Pi[2][1];				  
				    diff_a += d_grad_phi_i_e_a_dmesh[0][2] [b][j] * Pi[2][0];				
								  
				    diff_b += grad_phi_i_e_a[2][2] * d_Pi->X[2][2][b][j];
				    diff_b += grad_phi_i_e_a[2][1] * d_Pi->X[1][2][b][j];
				    diff_b += grad_phi_i_e_a[2][0] * d_Pi->X[0][2][b][j];
				    diff_b += grad_phi_i_e_a[1][2] * d_Pi->X[2][1][b][j];
				    diff_b += grad_phi_i_e_a[0][2] * d_Pi->X[2][0][b][j];
								  
				    diff_c += grad_phi_i_e_a[2][2] * Pi[2][2];					
				    diff_c += grad_phi_i_e_a[2][1] * Pi[1][2];					
				    diff_c += grad_phi_i_e_a[2][0] * Pi[0][2];					
				    diff_c += grad_phi_i_e_a[1][2] * Pi[2][1];					
				    diff_c += grad_phi_i_e_a[0][2] * Pi[2][0];					
				  }
							  
							  
#endif
				diff_a *= - d_area;
				diff_b *= -d_area ;
				diff_c *= - wt * ( d_det_J_dmesh_bj * h3 +
												 
						   det_J * dh3dmesh_bj );
				diffusion = diff_a + diff_b + diff_c;
							  
				diffusion *= diffusion_etm;
			      }
						  
			    /*
			     * Source term...
			     */
						  
			    source = 0.;	  
			    if ( source_on)
			      {
				source += phi_i * wt * 
				  ( f[a]   *        d_det_J_dmesh_bj * h3 +
				    f[a]   *        det_J *            dh3dmesh_bj + 
				    df->X[a][b][j] * det_J *            h3);
				source *= source_etm;
			      }
						  
			    if (particle_momentum_on)
			      {
				mass *= (1.0 - p_vol_frac);
				advection *= (1.0 - p_vol_frac);
			      }

			    continuity_stabilization = 0.0;
			    if (Cont_GLS)
			      {
				for(p=0; p<VIM; p++)
				  {
				    continuity_stabilization += d_grad_phi_i_e_a_dmesh[p][p][b][j] * cont_gls*d_area;
				    continuity_stabilization += grad_phi_i_e_a[p][p] * d_cont_gls->X[b][j]*d_area;
				    continuity_stabilization += grad_phi_i_e_a[p][p]*cont_gls *
				                                (d_det_J_dmesh_bj *h3 + det_J *dh3dmesh_bj)*wt;
				  }
				
			      }

			    J[j] += mass + advection + porous + diffusion + source + continuity_stabilization;
			  }

			/* Special Case for shell lubrication velocities with bounding continuum.  Note j-loop is over bulk element */			
			if ( mp->FSIModel == FSI_MESH_CONTINUUM || 
			     mp->FSIModel == FSI_REALSOLID_CONTINUUM) {
			  for ( j = 0; j < ei->dof[MESH_DISPLACEMENT1]; j++ ) {
			    jk=dof_map[j];
			    
			    phi_j = bf[var]->phi[j];
			    dh3dmesh_bj = fv->dh3dq[b] * phi_j;

			    porous = 0.;
			    if ( porous_brinkman_on ) {
			      porous += (v[a] - LubAux->v_avg[a]) * fv->dsurfdet_dx[b][jk] * h3;
			      porous += (v[a] - LubAux->v_avg[a]) * fv->sdet  * dh3dmesh_bj;
			      porous -= LubAux->dv_avg_dx[a][b][j] * fv->sdet * h3;
			    }
			    porous *= phi_i * wt;

			    J[jk] += porous;
			  }
			}

		      }
		  }
	      }
	    } /* end of if(active_dofs) */
	  } /* end of for(i=ei->dof*/ 
	}
    }
  safe_free((void *) n_dof);
  return(status);
}

/* assemble_continuity -- assemble Residual &| Jacobian for continuity eqns
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
 *	a  -- gets loaded up with proper contribution 
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Mar  2 09:27:30 MST 1994 pasacki@sandia.gov
 *
 * Revised:	Sun Mar 20 13:24:50 MST 1994 pasacki@sandia.gov
 */
int 
assemble_continuity(dbl time_value,   /* current time */
		    dbl tt,	/* parameter to vary time integration from
				   explicit (tt = 1) to implicit (tt = 0)    */
		    dbl dt,	/* current time step size                    */
		    const PG_DATA *pg_data )
{
  int dim, wim;
  int p, q, a, b;

  int eqn, var;
  int peqn, pvar;
  int w;

  int i, j;
  int status, err;

  dbl time = 0.0; /*  RSL 6/6/02  */

  dbl *v = fv->v;			/* Velocity field. */
  dbl div_v = fv->div_v;		/* Divergence of v. */

  dbl epsilon = 0.0, derivative, sum;		/*  RSL 7/24/00  */
  dbl sum1, sum2;			/*  RSL 8/15/00  */
  dbl sum_a, sum_b;                     /*  RSL 9/28/01  */
  int jj;				/*  RSL 7/25/00  */

  dbl advection;
  dbl source;
  dbl pressure_stabilization;

  dbl volsolvent=0;		/* volume fraction of solvent                */
  dbl initial_volsolvent=0;	/* initial solvent volume fraction
				 * (in stress-free state) input as source
				 * constant from input file                  */

  dbl det_J;
  dbl h3;
  dbl wt;
  dbl d_area;

  dbl d_h3detJ_dmesh_bj;		/* for specific (b,j) mesh dof */

  /*
   * Galerkin weighting functions...
   */

  dbl phi_i;
  dbl phi_j;
  dbl div_phi_j_e_b;
  dbl ( *grad_phi )[DIM];               /* weight-function for PSPG term */

  dbl div_v_dmesh;			/* for specific (b,j) mesh dof */

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];
  int mode;
  
  int *pde = pd->e;
  int *pdv = pd->v;

  dbl pspg[DIM];
  PSPG_DEPENDENCE_STRUCT d_pspg_struct;
  PSPG_DEPENDENCE_STRUCT *d_pspg = &d_pspg_struct;

  dbl mass, mass_a;
  dbl source_a;
  dbl sourceBase = 0.0;

  dbl rho=0;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Species_Conservation_Terms s_terms;
  dbl rhos=0, rhof=0;
  dbl h_flux=0;
  int w0=0;

  /* For particle momentum model.
   */
  int species;			/* species number for particle phase,  */
  dbl ompvf;			/* 1 - partical volume fraction */
  int particle_momentum_on; 	/* boolean. */

  /* Foaming model TAB */
  double dFVS_dv[DIM][MDE];
  double dFVS_dT[MDE];
  double dFVS_dx[DIM][MDE];
  double dFVS_dC[MAX_CONC][MDE];
  double dFVS_dF[MDE];

  int transient_run = pd->TimeIntegration != STEADY;
  int advection_on =0;
  int source_on =0;
  int ion_reactions_on=0, electrode_kinetics_on=0;
  int lagrangian_mesh_motion=0, total_ale_on=0;
  int hydromassflux_on=0, suspensionsource_on=0;
  int foam_volume_source_on = 0;
  
  dbl advection_etm,  source_etm;
  
  double *J;
  
  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn   = R_PRESSURE;
  peqn = upd->ep[eqn];

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;
  
  if (pd->v[POLYMER_STRESS11])
    {
      err = stress_eqn_pointer(v_s);
	  
      v_g[0][0] = VELOCITY_GRADIENT11;
      v_g[0][1] = VELOCITY_GRADIENT12;
      v_g[1][0] = VELOCITY_GRADIENT21;
      v_g[1][1] = VELOCITY_GRADIENT22;
      v_g[0][2] = VELOCITY_GRADIENT13;
      v_g[1][2] = VELOCITY_GRADIENT23;
      v_g[2][0] = VELOCITY_GRADIENT31;
      v_g[2][1] = VELOCITY_GRADIENT32;
      v_g[2][2] = VELOCITY_GRADIENT33; 
    }
  
  wt = fv->wt;
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */
  h3 = fv->h3;			/* Differential volume element (scales). */
  
  d_area = wt * det_J * h3;

  grad_phi = bf[eqn]->grad_phi;

  /*
   * Get the deformation gradients and tensors if needed
   */
   
  lagrangian_mesh_motion = (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN);
  electrode_kinetics_on = (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS);
  ion_reactions_on = ( mp->SpeciesSourceModel[0]  == ION_REACTIONS );
  total_ale_on = (cr->MeshMotion == TOTAL_ALE); 
  hydromassflux_on = (cr->MassFluxModel == HYDRODYNAMIC);
  suspensionsource_on = (mp->MomentumSourceModel == SUSPENSION );
  
  if ( lagrangian_mesh_motion && pde[R_MESH1])
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if (err == 2) return(err);
    }
  

  if (( total_ale_on && !pdv[VELOCITY1]) && pde[R_SOLID1])
    {
      err = belly_flop_rs(elc_rs->lame_mu);
      EH(err, "error in belly flop for real solid");
      if (err == 2) return(err);
    }

  if(pd->e[R_PMOMENTUM1])
    {
      particle_momentum_on = 1;
      species = (int) mp->u_density[0];
      ompvf = 1.0 - fv->c[species];
    }
  else
    {
      particle_momentum_on = 0;
      species = -1;
      ompvf = 1.0;
    }

  if(PSPG)
    {
      calc_pspg( pspg, d_pspg,
	         time_value, tt, dt,
	         pg_data);
    }


  if ((lagrangian_mesh_motion ||(total_ale_on && !pd->v[VELOCITY1]))
      && (mp->PorousMediaType == CONTINUOUS))
    {
      initial_volsolvent = elc->Strss_fr_sol_vol_frac;
      volsolvent = 0.;
      for (w=0; w<pd->Num_Species_Eqn; w++) volsolvent += fv->c[w];
      if(particle_momentum_on)
	volsolvent -= fv->c[species];
    }

  

  
  if ( electrode_kinetics_on ||  ion_reactions_on ) 
    {
      if (mp->PorosityModel == CONSTANT)
	{
	  epsilon = mp->porosity;
	}
      else if (mp->PorosityModel == THERMAL_BATTERY)
	{
	  epsilon = mp->u_porosity[0];
	}
      else
	{
	  EH(-1, "invalid porosity model");
	}
    }
  
  if ( mp->MomentumSourceModel == SUSPENSION_PM ||
       electrode_kinetics_on ) /*  RSL 7/25/00  */
    {
      err = get_continuous_species_terms(&s_terms, 0.0, tt, dt, pg_data->hsquared);
      EH(err,"problem in getting the species terms");
    }

  if ( ion_reactions_on ) /*  RSL 3/19/01 and 6/6/02  */
    {
      zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
      err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
      EH(err,"problem in getting the species terms");
    }

	  
  if( (hydromassflux_on) && ( mp->DensityModel == SUSPENSION )&& ( suspensionsource_on ))
    {
      /*
       * Compute hydrodynamic/sedimentation flux and sensitivities.
       */


      w0 = (int) mp->u_density[0]; /* This is the species number that is transported HYDRODYNAMICally  */

      hydro_flux( &s_terms, w0, tt, dt, pg_data->hsquared);

      rhof = mp->u_density[1];
      rhos = mp->u_density[2];

    }
  rho  = density(d_rho, time_value);
	  
  advection_on = pd->e[eqn] & T_ADVECTION ;
  source_on = pd->e[eqn] & T_SOURCE ;

	  
  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
  source_etm = pd->etm[eqn][(LOG2_SOURCE)];
	  

  if (af->Assemble_Residual)
    {
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	 
#if 1
	  /* this is an optimization for xfem */

	  if ( xfem != NULL )
	    {
	      int xfem_active, extended_dof, base_interp, base_dof;
			  
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
			      &xfem_active, &extended_dof, &base_interp, &base_dof );
			  
	      if ( extended_dof && !xfem_active ) continue;
	    }
#endif

	  phi_i      = bf[eqn]->phi[i];

	  /* 
	   *  Mass Terms: drhodt terms (usually though problem dependent)
	   */
	  mass = 0.0;
	  if (particle_momentum_on)
	    {
	      if (transient_run)
		{
		  mass = -s_terms.Y_dot[species];
		  mass *= phi_i * d_area;
		}
	    }
	  if (electrode_kinetics_on || ion_reactions_on) 
	    {
	      mass = 0.0;
	      if (transient_run)
		{
		  var = MASS_FRACTION;
		  for (w = 0; w < pd->Num_Species-1; w++)
		    {
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  if ( bf[var]->phi[j] > 0.0 ) break;
			}
		      derivative = d_rho->C[w][j]/bf[var]->phi[j];
		      mass += derivative * s_terms.Y_dot[w];
		    }
		  mass *= epsilon/rho;
		  mass *= phi_i * d_area;
		}
	    }
	  
	  /*
	   *  Advection: 
	   *    This term refers to the standard del dot v . 
	   *    
	   *    int (phi_i div_v d_omega)
	   *
	   *   Note density is not multiplied into this term normally
	   */
	  advection = 0.0;
	  if (advection_on)
	    {
	      if (pdv[VELOCITY1]) /* then must be solving fluid mechanics in this material */
		{

		  /*
		   * Standard incompressibility constraint means we have
		   * a solenoidal velocity field
		   */

		  advection = div_v;

		  /* We get a more complicated advection term because the
		   * particle phase is not incompressible.
		   */
		  if (particle_momentum_on) advection *= ompvf;

		  advection *= phi_i * d_area;
		  advection *= advection_etm;
		}
	      else if (lagrangian_mesh_motion || total_ale_on)
		/* use divergence of displacement for linear elasticity */
		{
		  advection = fv->volume_change;

		  if( particle_momentum_on ) advection *= ompvf;

		  advection *= phi_i * h3 * det_J * wt;
		  advection *= advection_etm;
		}
	      if (electrode_kinetics_on || ion_reactions_on)
		{
		  advection = div_v;
		  var = MASS_FRACTION;
		  for (w=0; w<pd->Num_Species-1; w++)
		    {
		      for ( j=0; j<ei->dof[var]; j++)
			{
			  if ( bf[var]->phi[j] > 0.0 ) break;
			}
		      derivative = d_rho->C[w][j]/bf[var]->phi[j];
		      sum = 0.;
		      for (p=0; p<dim; p++)
			{
			  sum += s_terms.conv_flux[w][p];
			}
		      advection += derivative * sum / rho;
		    }
		  advection *= phi_i *d_area;
		  advection *= advection_etm;
		}
	    }

	  source = 0.0;
	  sourceBase = 0.0;
	  if (source_on)
	    {
	      if (pdv[VELOCITY1])
		{
		  /* DRN (07/13/05):
		     This was previously:
		     source     =  P;
		     But this messes with level set problems that have a density source
		     over part of the domain and constant in other regions.  If someone was
		     counting on this behavior as a form of a penalty method to give a non-zero
		     diagonal entry, we should implement a new density model that accomplishes this.
		     I really don't know what you want for DENSITY_IDEAL_GAS, though?!?
		      
		     source     =  0.;
		     }*/
		  if ( mp->DensityModel == DENSITY_FOAM || 
		       mp->DensityModel == DENSITY_FOAM_CONC || 
		       mp->DensityModel == DENSITY_FOAM_TIME ||
		       mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
		    {
		      /* These density models locally permit a time and spatially varying
			 density.  Consequently, the Lagrangian derivative of the density
			 terms in the continuity equation are not zero and are
			 included here as a source term
		      */
		      source = FoamVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT,
						dFVS_dx, dFVS_dC, dFVS_dF);
		      sourceBase = source;
		      foam_volume_source_on = 1;
		    }
		  else if ( mp->DensityModel == REACTIVE_FOAM )
		    {
		      /* These density models locally permit a time and spatially varying
			 density.  Consequently, the Lagrangian derivative of the density
			 terms in the continuity equation are not zero and are
			 included here as a source term
		      */
		      source = REFVolumeSource( time_value,
						dt,
						tt,
						dFVS_dv,
						dFVS_dT,
						dFVS_dx,
						dFVS_dC );
		      sourceBase  = source;
		      foam_volume_source_on =  1;
		    }
			  
		  /*
		    else if
		    ( mp->DensityModel == SUSPENSION ||
		    mp->DensityModel == SUSPENSION_PM )
		    {
		  */
		  /* Although the suspension density models meet the definition of
		     a locally variable density model, the Lagrangian derivative
		     of their densities can be represented as a divergence of
		     mass flux.  This term must be integrated by parts and so is
		     included separately later on and is not include as part of the "source"
		     terms 
		     source = 0.0;
		     }  */
				

		  /* To include, or not to include, that is the question
		   * when considering the particle momentum coupled eq.'s...
		   */

		  /* We get this source term because the fluid phase is not
		   * incompressible.
		   */
		  if (particle_momentum_on)
		    {
		      source = 0.0;
 		      for(a = 0; a < wim; a++ )
			{
			  /* Cannot use s_terms.conv_flux[a] here because that
			   * is defined in terms of the particle phase
			   * velocities.
			   */
			  source -= fv->grad_c[species][a] * v[a];
			}
		      sourceBase = source;
		    }
		  source *= phi_i * d_area;
		  source *= source_etm;
		}

	      if ((lagrangian_mesh_motion || (total_ale_on && !pd->v[VELOCITY1]) ))
		/* add swelling as a source of volume */
		{
		  if ( mp->PorousMediaType == CONTINUOUS )
		    {
		      source = - (1. - initial_volsolvent) / (1. - volsolvent);
		      sourceBase = source;
		      source *= phi_i * d_area;
		      source *= source_etm;

		    }
		}
	      if (electrode_kinetics_on || ion_reactions_on)
	        {
		  source = 0.0;
		  for (j=0; j<pd->Num_Species; j++)
	            {
		      source -= s_terms.MassSource[j] * mp->molecular_weight[j];
	            }
		  source /= rho;
		  sourceBase = source;
		  source *= phi_i * d_area;
		  source *= source_etm;
	        }
	    }
	  /* add Pressure-Stabilized Petrov-Galerkin term
	   * if desired.
	   */

	  pressure_stabilization = 0.0;
	  if (PSPG)
	    {
	      for ( a=0; a<wim; a++)
		{
		  meqn = R_MOMENTUM1+a;
		  if( pd->e[meqn])
		    {
		      pressure_stabilization += grad_phi[i][a] * pspg[a];
		    }
		}
	      pressure_stabilization *= d_area;
	    }

	  h_flux = 0.0;
	  if ( (hydromassflux_on) && ( suspensionsource_on))
	    {	
	      /* add divergence of particle phase flux as source term */
		
	      /* The particle flux terms has been integrated by parts.
	       * No boundary integrals are included in this formulation
	       * so it is tacitly assumed that the particle phase
	       * relative mass flux over all boundaries is zero
	       */
		  
	      for( p=0; p<dim; p++)
		{
		  h_flux += grad_phi[i][p]*s_terms.diff_flux[w0][p];
		}
		      
	      h_flux *= (rhos - rhof )/rhof;
	      h_flux *= h3 * det_J * wt;
	      h_flux *= pd->etm[eqn][(LOG2_ADVECTION)];
	      /*  h_flux = 0.0; */
	    } 

#ifdef DEBUG_CONTINUITY_RES
	  printf("R_c[%d] += %10f + %10f + %10f + %10f + %10f\n",i,mass,advection,source,
		 pressure_stabilization,h_flux);
#endif
	  /*
	   *  Add up the individual contributions and sum them into the local element
	   *  contribution for the total continuity equation for the ith local unknown
	   */
	  lec->R[peqn][i] += mass + advection + source + pressure_stabilization + h_flux; 
	}
    }

  if (af->Assemble_Jacobian)
    {
      for (i = 0; i < ei->dof[eqn]; i++)
	{ 
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
              int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
              if ( extended_dof && !xfem_active ) continue;
            }
#endif

	  phi_i = bf[eqn]->phi[i];

	  /* 
	   * J_c_v NOTE that this is applied whenever velocity is a variable
	   */
	  for ( b=0; b<wim; b++)
	    {
	      var = VELOCITY1+b;
	      if ( pdv[var] )
		{
		  pvar = upd->vp[var];
			  
		  J = lec->J[peqn][pvar][i];
			  
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];	      
				  
		      h_flux = 0.;
		      if((hydromassflux_on) && ( suspensionsource_on) )
			{
			  for( p=0; p<dim; p++)
			    {
			      h_flux += grad_phi[i][p] * s_terms.d_diff_flux_dv[w0][p][b][j] ;
			    }
			  h_flux *= ( rhos - rhof )/rhof * d_area * advection_etm;	
					  
			}
				  
		      advection  = 0.;
				  
		      if ( advection_on )
			{
			  div_phi_j_e_b = 0.;
#ifdef DO_NO_UNROLL
			  for ( p=0; p<VIM; p++)
			    {
			      div_phi_j_e_b += bf[var]->grad_phi_e[j][b] [p][p];
			    }
#else
			  div_phi_j_e_b += bf[var]->grad_phi_e[j][b] [0][0] + bf[var]->grad_phi_e[j][b] [1][1];
			  if( VIM == 3 ) div_phi_j_e_b += bf[var]->grad_phi_e[j][b] [2][2];
#endif
					  
			  if ( electrode_kinetics_on || ion_reactions_on ) 
			    {
			      sum = 0.;
			      for (jj=0; jj<pd->Num_Species-1; jj++)
				{
				  for ( q=0; q<ei->dof[MASS_FRACTION]; q++)
				    {
				      if ( bf[MASS_FRACTION]->phi[q] > 0.0 ) break;
				    }
				  derivative = d_rho->C[jj][q]/bf[MASS_FRACTION]->phi[q];
				  sum += derivative * s_terms.d_conv_flux_dv[jj][b][b][j];
				}
			      div_phi_j_e_b += sum/rho;
			    }
					  
			  advection = phi_i * div_phi_j_e_b * d_area;
					  
			  if(particle_momentum_on)
			    advection *= ompvf;
					  
			  advection *= advection_etm;
			}
				  
				  
		      source = 0.;
				  
		      if( source_on)
			{
					  
			  if ( foam_volume_source_on )
			    {
			      source = dFVS_dv[b][j];
			    }
					  
			  if( particle_momentum_on  )
			    {
			      /* From residual calculation.
				 source -= fv->grad_c[species][a] * v[a];
			      */
			      source = -s_terms.grad_Y[species][b] * phi_j;
			    }
					  
			  source *= phi_i * d_area;
			  source *= source_etm;
			}
				  
		      /* add Pressure-Stabilized Petrov-Galerkin term
		       * if desired.
		       */
		      pressure_stabilization = 0.;  
		      if (PSPG)
			{
			  for ( a=0; a<wim; a++)
			    {
			      meqn = R_MOMENTUM1+a;
			      if( pd->e[meqn])
				{
				  pressure_stabilization += grad_phi[i][a] * d_pspg->v[a][b][j];
				}
			    }
			  pressure_stabilization *= d_area;
			}
				  
		      J[j] += advection + source + pressure_stabilization + h_flux; 
		      /* lec->J[peqn][pvar][i][j] += advection + source + pressure_stabilization + h_flux; */
#ifdef DEBUG_CONTINUITY_JAC
		      printf("J_c_v[%d] [%d][%d] += %10f + %10f + %10f + %10f\n",
			     i,b,j,advection,source,pressure_stabilization,
			     h_flux);
#endif
				  
				  
		    }
		}
	    }

	  /*
	   * J_c_T This term comes from the temperature dependency of the momentume source
	   * which comes from the Pressure-Stabilized Petrov-Galerkin term
	   */
	  var = TEMPERATURE;
	  if ( PSPG && pd->v[var] )
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  pvar = upd->vp[var];

		  phi_j = bf[var]->phi[j];
		  
		  /* add Pressure-Stabilized Petrov-Galerkin term
		   * if desired.
		   */
		  pressure_stabilization = 0.;  

		  for ( a=0; a<wim; a++)
		    {
		      meqn = R_MOMENTUM1+a;
		      if( pd->e[meqn])
			{
			  pressure_stabilization += grad_phi[i][a] * d_pspg->T[a][j];
			}
		    }
		  pressure_stabilization *= d_area;
		      
		  lec->J[peqn][pvar][i][j] += pressure_stabilization;
		}
	    }

	  if( source_on )
	    {
	      if( foam_volume_source_on )
		{
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      pvar = upd->vp[var];
		      
		      lec->J[peqn][pvar][i][j] += d_area*phi_i*dFVS_dT[j]* source_etm;
		    }
		}
	    }

	  if ( (hydromassflux_on ) &&  pdv[var] )
	    {
	      if ( suspensionsource_on )
		{
		  pvar = upd->vp[var];

		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      h_flux = 0.;

		      for ( p=0; p<dim ; p++)
			{
			  h_flux += grad_phi[i][p]*s_terms.d_diff_flux_dT[w0][p][j];
			}

		      h_flux *= d_area*( rhos - rhof )/rhof;

		      /*  h_flux = 0.0; */	  	  
	       
		      lec->J[peqn][pvar][i][j] += h_flux;
		    }
		}
	    }

	  /*
	   *  J_c_F this is primarily the foam volume source terms 
	   */

	  var = FILL;
	  if( source_on && pdv[var] && ls != NULL)
	    {
	      pvar = upd->vp[var];
	      
	      J = lec->J[peqn][pvar][i];
	      
	      for ( j=0; j<ei->dof[var]; j++)
		{			
		  source = 0.0;
		  
		  if ( foam_volume_source_on )
		    {
		      source = dFVS_dF[j];
		      source *= phi_i*d_area;
		    }
		  
		  J[j] += source;
		  
		  /*lec->J[peqn][pvar][i][j] += source;*/
		  
		}
	    }
			
		
	  /*
	   * J_c_P here species act as a volume source in continuous lagrangian mesh motion
	   */
	  var = PRESSURE;
	  if ( pdv[var] )
	    {
	      pvar = upd->vp[var];
	      
	      J = lec->J[peqn][pvar][i];
	      
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  
		  advection  = 0.;
		  
		  if (advection_on && (lagrangian_mesh_motion ||
				       (total_ale_on && !pd->v[VELOCITY1])))
		    {
                      /*Need to compute this for total ALE.  Not done yet */
		      advection = fv->d_volume_change_dp[j];
		      
		      advection *= phi_i * d_area;
		      
		      advection *= advection_etm;
		    }
		  
		  source = 0.;
		  
		  /* add Pressure-Stabilized Petrov-Galerkin term
		   * if desired.
		   */
		  pressure_stabilization = 0.;
		  if(PSPG)
		    {
		      for ( a=0; a<wim; a++)
			{
			  meqn = R_MOMENTUM1 + a;
			  if ( pd->e[meqn] & T_DIFFUSION )
			    { 
			      pressure_stabilization += grad_phi[i][a] * d_pspg->P[a][j];
			    }
			}
		      pressure_stabilization *= d_area;
		      
		    }
		  J[j] += advection  + source + pressure_stabilization;
		  /*lec->J[peqn][pvar][i][j] += advection  + source + pressure_stabilization; */
#ifdef DEBUG_CONTINUITY_JAC
		  printf("J_c_P[%d] [%d] += %10f %10f %10f\n",i,j,
			 advection,source,pressure_stabilization);
#endif
		}
	    }  
	  
	  /*
	   * J_c_S this term is only present for PSPG
	   */
	  var = POLYMER_STRESS11;
	  if ( PSPG && pdv[var] )
	    {
	      for ( mode=0; mode<vn->modes; mode++)
		{
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{
			  var = v_s[mode][p][q];
			  pvar = upd->vp[var];
			  if ( pd->v[var] )
			    {
			      for ( j=0; j<ei->dof[var]; j++)
				{
				  pressure_stabilization = 0.;
				  
				  for ( a=0; a<wim; a++)
				    {
				      meqn = R_MOMENTUM1+a;
				      if( pd->e[meqn])
					{
					  pressure_stabilization += grad_phi[i][a] * d_pspg->S[a][mode][p][q][j];
					}
				    }
				  
				  pressure_stabilization *=  h3 * det_J * wt;
			      
				  lec->J[peqn][pvar][i][j] += pressure_stabilization;
				}
			    }
			}
		    }
		}
	    }
	  
	  /*
	   * J_c_G this term is only present for PSPG
	   */
	  var = VELOCITY_GRADIENT11;
	  if ( PSPG && pdv[var] )
	    {
	      for ( p=0; p<VIM; p++)
		{
		  for ( q=0; q<VIM; q++)
		    {
		      var = v_g[p][q];
		      pvar = upd->vp[var];
		      if ( pd->v[var] )
			{
			  for ( j=0; j<ei->dof[var]; j++)
			    {
			      pressure_stabilization = 0.;
				  
			      for ( a=0; a<wim; a++)
				{
				  meqn = R_MOMENTUM1+a;
				  if( pd->e[meqn])
				    {
				      pressure_stabilization += grad_phi[i][a] * d_pspg->g[a][p][q][j];
				    }
				}
			      pressure_stabilization *=  h3 * det_J * wt;
			      
			      lec->J[peqn][pvar][i][j] += pressure_stabilization;   
			      
			    }
			}
		    }
		}
	    }


	  /*
	   * J_c_SH this term is only present for HYDRODYNAMIC mass flux and SUSPENSION
	   *        momentum source
	   */

	  var = SHEAR_RATE;

	  if ( hydromassflux_on &&  pdv[var] )
	    {
	      if ( suspensionsource_on )
		{
		  if ( pd->v[var] )
		    {
		      pvar = upd->vp[var];

		      for ( j=0; j<ei->dof[var]; j++)
			{
			  h_flux = 0.0;
			  for ( a=0; a<dim ; a++)
			    {
			      h_flux += grad_phi[i][a]*s_terms.d_diff_flux_dSH[w0][a][j];
			    }

			  h_flux *= h3*det_J*wt*( rhos - rhof )/rhof;

			  /*  h_flux = 0.0; */	  	  

			  lec->J[peqn][pvar][i][j] += h_flux;
			}
		    }  
		}
	    }

  
	  /*
 	   * J_c_d
	   */
	  
	  for (b = 0; b < dim; b++)
	    {
	      var = MESH_DISPLACEMENT1 + b;
	      if (pdv[var])
		{
		  pvar = upd->vp[var];
		  J = lec->J[peqn][pvar][i];
		  for (j = 0; j < ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      /* derivative of |J| with extra term for axisymmetry e.g. 
			 d/dmesh [ r|J| ] */  
		      d_h3detJ_dmesh_bj =(h3 * bf[eqn]->d_det_J_dm[b][j] +
					  det_J * fv->dh3dq[b] * bf[var]->phi[j]);

                      mass = 0.0;  
		      if (electrode_kinetics_on || ion_reactions_on ) 
			{
			  if ( transient_run)
			    {
			      for (w=0; w<pd->Num_Species-1; w++)
				{
				  for ( q=0; q<ei->dof[MASS_FRACTION]; q++)
				    {
				      if ( bf[MASS_FRACTION]->phi[q] > 0.0 ) break;
				    }
				  derivative = d_rho->C[w][q]/bf[MASS_FRACTION]->phi[q];
				  mass += derivative * s_terms.Y_dot[w];
				}
			      mass *= epsilon/rho;
			      mass *= phi_i * d_h3detJ_dmesh_bj * wt;
			    }
			}

		      advection  = 0.0;
		      if (advection_on)
			{
			  if (pdv[VELOCITY1])
			    {
			      h_flux = 0.0;
			      if((hydromassflux_on) && ( suspensionsource_on ) )
				{
				  for( p=0; p<dim; p++)
				    {
				      h_flux += grad_phi[i][p]*s_terms.diff_flux[w0][p]
					* d_h3detJ_dmesh_bj + grad_phi[i][p]* s_terms.d_diff_flux_dmesh[w0][p][b][j] 
					* det_J * h3 + bf[eqn]->d_grad_phi_dmesh[i][p] [b][j] * s_terms.diff_flux[w0][p] * det_J * h3;
				    }
				  h_flux *= ( rhos - rhof )/rhof * wt * pd->etm[eqn][(LOG2_ADVECTION)];	
	      			  
				}

			      div_v_dmesh = fv->d_div_v_dmesh[b][j];
		      
			      
			      advection+= div_v_dmesh * det_J * h3 + div_v * ( d_h3detJ_dmesh_bj );
				  
			      if (electrode_kinetics_on || ion_reactions_on ) /*  RSL  9/28/01  */
				{
				  sum_a = 0.;
				  sum_b = 0.;
				  for (w=0; w<pd->Num_Species-1; w++)
				    {
				      for ( q=0; q<ei->dof[MASS_FRACTION]; q++)
					{
					  if ( bf[MASS_FRACTION]->phi[q] > 0.0 ) break;
					}
				      derivative = d_rho->C[w][q]/bf[MASS_FRACTION]->phi[q];
				      sum1 = 0.;
				      sum2 = 0.;
				      for (p=0; p<dim; p++)
					{
					  sum1 += s_terms.conv_flux[w][p];
					  sum2 += s_terms.d_conv_flux_dmesh[w][p][b][j];
					}
				      sum_a += derivative * sum1;
				      sum_b += derivative * sum2;
				    }
				  sum_a /= rho;
				  sum_b /= rho;
				  advection += sum_b * det_J * h3 + sum_a * d_h3detJ_dmesh_bj;
				}
					  
			      advection*= phi_i * wt;
			    } 
			  else if (lagrangian_mesh_motion || ( total_ale_on && !pd->v[VELOCITY1]) )
			    {
			      advection += fv->volume_change 
				* ( d_h3detJ_dmesh_bj ); 
					  
			      advection += fv->d_volume_change_dx[b][j] * 
				h3 * det_J; 
					  
			      advection *= phi_i * wt;
			    }
			  
			  advection *= advection_etm;
			}
		      
		      source = 0.0;
		      if (source_on)
			{
			  if( mp->DensityModel == DENSITY_FOAM || 
			      mp->DensityModel == DENSITY_FOAM_CONC || 
			      mp->DensityModel == DENSITY_FOAM_TIME ||
			      mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
			    {
			      source = sourceBase * d_h3detJ_dmesh_bj * wt + dFVS_dx[b][j] * d_area;
			      source *= phi_i * source_etm;
			    }
			  else if (mp->DensityModel == REACTIVE_FOAM)
			    {
			      source = sourceBase * d_h3detJ_dmesh_bj * wt + dFVS_dx[b][j] * d_area;
			      source *= phi_i * source_etm;
			    }
			}
		      if ( lagrangian_mesh_motion || (total_ale_on  && !pdv[VELOCITY1])  )
			{
			  /* add swelling as a source of volume */
			  if ( mp->PorousMediaType == CONTINUOUS )
			    {
			      source = - phi_i * (d_h3detJ_dmesh_bj) * wt 
				* (1. - initial_volsolvent) / (1. - volsolvent) 
				* source_etm;
			    }
			}

		      if (electrode_kinetics_on || ion_reactions_on ) 
			{
			  sum1 = 0.;
			  sum2 = 0.;
			  for (q=0; q<pd->Num_Species; q++)
			    {
			      sum1 -= s_terms.MassSource[q] * mp->molecular_weight[q];
			      sum2 -= s_terms.d_MassSource_dmesh[q][b][j] * mp->molecular_weight[q];
			    }
			  sum1 /= rho;
			  sum2 /= rho;
			  source += sum2 * det_J * h3 + sum1 * d_h3detJ_dmesh_bj;
			  source *= phi_i * wt;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      /* add Pressure-Stabilized Petrov-Galerkin term
		       * if desired.
		       */
		      pressure_stabilization = 0.0;  
		      if (PSPG)
			{
			  for ( a=0; a<wim; a++)
			    {
			      meqn = R_MOMENTUM1+a;
			      if( pd->e[meqn])
				{
				  pressure_stabilization += grad_phi[i][a] * d_pspg->X[a][b][j] * h3 * det_J * wt
				    + grad_phi[i][a] * pspg[a] *  wt * d_h3detJ_dmesh_bj
				    + bf[eqn]->d_grad_phi_dmesh[i][a][b][j] * pspg[a] * wt  * h3 * det_J;
				}
			    }
			}

		      /*lec->J[peqn][pvar][i][j] += advection  + source + pressure_stabilization +  h_flux + mass ; */
		      J[j] += advection  + source + pressure_stabilization +  h_flux + mass ; 
		    }
		}
	    }

	  /*
 	   * J_c_d_rs
	   */

	  for ( b=0; b<dim; b++)
	    {
	      var = SOLID_DISPLACEMENT1+b;
	      if ( pdv[var] )
		{
		  pvar = upd->vp[var];

		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      advection  = 0.;

		      if ( pd->e[eqn] & T_ADVECTION )
			{
			  if (cr->MeshMotion == TOTAL_ALE && !pd->v[VELOCITY1])
			    {
			      advection += fv->d_volume_change_drs[b][j] * 
				h3 * det_J; 
			      
			      advection *= phi_i * wt;
			    }
			  
			  advection *= pd->etm[eqn][(LOG2_ADVECTION)];
			}
		      
		      source = 0.;

		      lec->J[peqn][pvar][i][j] += advection  + source; 

		    }
		}
	    }

	  /*
	   * J_c_c
	   */
	      
	  var = MASS_FRACTION;
	  if ( pdv[var])
	    {
	      pvar = upd->vp[var];

	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  
		  for (w=0; w<pd->Num_Species_Eqn; w++) 
		    {
		      /* add swelling as a source of volume */
		      source = 0.;
		      if ( ( mp->PorousMediaType == CONTINUOUS ) && 
			   ( cr->MeshMotion == LAGRANGIAN ||
			     cr->MeshMotion == DYNAMIC_LAGRANGIAN))
			{
			  source = - phi_j * phi_i * h3 * det_J * wt 
			    * (1. - initial_volsolvent) 
			    / (1. - volsolvent) / (1. - volsolvent)
			    * pd->etm[eqn][LOG2_SOURCE];
			}
		      if (pd->e[eqn] & T_SOURCE)
			{
			  /* Foaming volume source term */

			  if (  mp->DensityModel == REACTIVE_FOAM || mp->DensityModel == DENSITY_FOAM )
			    {
			      source += dFVS_dC[w][j];
			      source *= phi_i * h3 * det_J * wt* pd->etm[eqn][LOG2_SOURCE];
			    }

			}

		      /* add Pressure-Stabilized Petrov-Galerkin term
		       * if desired.
		       */
		      pressure_stabilization = 0.;
		      if(PSPG)
			{
			  for ( a=0; a<wim; a++)
			    {
			      meqn = R_MOMENTUM1 + a;
			      if ( pd->e[meqn] & T_DIFFUSION )
				{
				  pressure_stabilization += grad_phi[i][a] * d_pspg->C[a][w][j];
				}
			    }
			  pressure_stabilization *= h3 * det_J * wt;
			}
		      
		      /* The fluid phase is not incompressible in the
		       * SUSPENION_PM model.
		       */
		      mass_a = 0.0;
		      advection = 0.0;
		      if( particle_momentum_on && w == species )
			{
			  mass_a = -(1.0 + 2.0*tt) * phi_j / dt;
			  mass_a *= phi_i * h3 * det_J * wt;
			  
			  if( pd->e[eqn] & T_ADVECTION )
			    {
			      advection = -phi_j * div_v;
			      advection *= phi_i * det_J * h3 * wt;
			    }

			  if( pd->e[eqn] & T_SOURCE )
			    {
			      source_a = 0.0;
			      for( a=0; a<wim; a++)
				source_a -= grad_phi[j][a]*v[a];
			      source_a *= phi_i * det_J * h3 * wt;
			      
			      source += source_a;
			    }
			}

		      if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
			  mp->SpeciesSourceModel[0]  == ION_REACTIONS) /*  RSL 3/19/01  */
			{
			  sum = 0.;
			  for (jj=0; jj<pd->Num_Species-1; jj++)
			    {
			      for ( q=0; q<ei->dof[MASS_FRACTION]; q++)
				{
				  if ( bf[MASS_FRACTION]->phi[q] > 0.0 ) break;
				}
			      derivative = d_rho->C[jj][q]/bf[MASS_FRACTION]->phi[q];
			      sum += derivative * (s_terms.d_Y_dot_dc[jj][w][j] -
						   d_rho->C[w][j]*s_terms.Y_dot[jj]/rho);
			    }
			  mass_a = sum*epsilon/rho;
			  mass_a *= phi_i * h3 * det_J * wt;
			}

		      if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
			  mp->SpeciesSourceModel[0]  == ION_REACTIONS) /*  RSL 3/19/01  */
			{
			  sum = 0.;
			  for (jj=0; jj<pd->Num_Species-1; jj++)
			    {
			      for ( q=0; q<ei->dof[MASS_FRACTION]; q++)
				{
				  if ( bf[MASS_FRACTION]->phi[q] > 0.0 ) break;
				}
			      derivative = d_rho->C[jj][q]/bf[MASS_FRACTION]->phi[q];
			      sum1 = 0.;
			      sum2 = 0.;
			      for ( p=0; p<VIM; p++)
				{
				  sum1 += s_terms.d_conv_flux_dc[jj][p][w][j];
				  sum2 += s_terms.conv_flux[jj][p];
				}
			      sum += derivative * (sum1 - d_rho->C[w][j]*sum2/rho);
			    }
			  advection = sum/rho;
			  advection *= phi_i * h3 * det_J * wt;
			}

		      if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
			  mp->SpeciesSourceModel[0]  == ION_REACTIONS) /*  RSL 3/19/01  */
			{
			  source_a = 0.0;

			  for (jj=0; jj<pd->Num_Species; jj++)
			    {
			      source_a += mp->molecular_weight[jj] *
				(s_terms.d_MassSource_dc[jj][w][j] -
				 d_rho->C[w][j]*s_terms.MassSource[jj]/rho);
			    }
			  source_a /= (-rho);
			  source_a *= phi_i * h3 * det_J * wt;

			  source += source_a;
			}

		      lec->J[peqn][MAX_PROB_VAR + w][i][j] += mass_a + advection + source + pressure_stabilization; 
		    }
		}
	    }
	  
	  if ( cr->MassFluxModel == HYDRODYNAMIC &&  pd->v[var] )
	    {
	      if ( mp->MomentumSourceModel == SUSPENSION )
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		  
		      h_flux = 0.0;
			  
		      for( a=0; a<dim; a++)
			{
			  h_flux +=  grad_phi[i][a]*s_terms.d_diff_flux_dc[w0][a][w0][j];
			}
		      
		      h_flux *= h3*det_J*wt*( rhos - rhof )/rhof;
		      /*  h_flux = 0.0; */

		      lec->J[peqn][MAX_PROB_VAR + w0][i][j] += h_flux;
		    }
		}
	    }
	}
    }
  return(status);
}

/******************************************************************************/
/* assemble_volume 
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *	lec-- pointer to Local Element Contribution	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution 
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Amy Sun 1/22/99
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_volume(bool owner)
{
  int eqn;
  int var;
  int dim;
  int b;

  int j;
  int status;
  int iAC=0;
  int vj;
  int VC_mode;
  int spec_id;
  int ktype, gnn, nvdof, ledof;

  dbl det_J;
  dbl h3;			/* Volume element (scale factors). */
  dbl wt;
  dbl func;
  dbl phi_j;
  dbl d_det_J_dmeshbj;
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */
  dbl d_detJh3_dmeshbj;
  dbl d_func_dmeshbj;           /* Sensitivity of func to (b,j)mesh dof. */

  status = 0;
  func=1.;

  /*
   * Unpack variables from structures for local convenience...
   */
  eqn = pd->ShapeVar;
  dim = pd->Num_Dim;
  h3  = fv->h3;			/* Differential volume element (scales). */
  det_J = bf[eqn]->detJ;

  iAC = pd->VolumeIntegral;
 
  VC_mode = augc[iAC].VOLID ;
  spec_id = augc[iAC].COMPID;


  if (pd_glob[0] ->CoordinateSystem != CYLINDRICAL &&
      pd_glob[0] ->CoordinateSystem != SWIRLING)
    {
      wt  = fv->wt;
    }
  else
    {
      wt = 2.0 * M_PIE * fv->wt;
    }

  /*
   * Load up the function the volume multiplies to
   * 1= pure volume calc for the element
   * 2= rho * volume
   * 3= mass fraction of species w * volume
   * Add as needed
   */

  switch (VC_mode)
    {
    case 1: /* pure volume calc */
    case 11:
      break;
    case 2: /* density*volume */
    case 12:
      func = mp->density;
      break;
    case 3: /* massfraction*volume */
    case 13:
      if (pd->Num_Species_Eqn < 1)
	{
	  fprintf(stderr, "Need at least one species to use this option"); 
	} 
      else  
	{
	  func = fv->c[spec_id];
	}
      break;
    case 4:
    case 14:
      if(pd->v[VELOCITY1])
	{
	  func = fv->v[0];
	}
      else
	{
	  EH(-1," must have momentum equation on for this AC");
	}
      break;

    default:
      EH(-1, "assemble_volue() unknown integral calculation\n");
      break;
    }

  /* total integral */

  if (owner)
    {
      augc[iAC].evol += func * det_J * wt * h3;
    }
    
  /* Calculate sensitivities so that we don't need to compute it numerically
     in mm_sol_nonlinear.c  */
  
  for (b = 0; b < dim; b++)
    {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[var])
	{
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      eqn =  R_MESH1+b;

	      phi_j = bf[var]->phi[j];

	      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];
	      
	      dh3dmesh_bj = fv->dh3dq[b] * phi_j;

	      d_detJh3_dmeshbj = d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj;

	      /*assume no dependency of func on meshbj, may change in the future*/

	      d_func_dmeshbj = 0.;

	      /*grab Index_Solution for location of x. see mm_fill_ptrs.c */

	      vj = ei->gun_list[var][j];

	      augc[iAC].d_evol_dx[vj] += wt*(func*d_detJh3_dmeshbj
					     + d_func_dmeshbj*det_J*h3);
	    }
	}
    }
     
  if (VC_mode == 3) 
    {
      var = MASS_FRACTION ;
      eqn = R_MASS;
      if (pd->v[var])
	{
	  for (j = 0; j < ei->dof[var]; j++) 
	    {
	      phi_j = bf[var]->phi[j];
	      ktype = spec_id;
	      gnn   = ei->gnn_list[eqn][j];
	      nvdof = ei->Baby_Dolphin[eqn][j];
	      ledof =  ei->lvdof_to_ledof[eqn][j];
	      vj    = Index_Solution(gnn, eqn, ktype , nvdof,
				     ei->matID_ledof[ledof]);
	      augc[iAC].d_evol_dx[vj] += phi_j * det_J * wt * h3;
	    }
	}
    }
  if (VC_mode == 4) 
    {
      var = VELOCITY1 ;
      eqn = R_MOMENTUM1;
      if (pd->v[var])
	{
	  for (j = 0; j < ei->dof[var]; j++) 
	    {
	      phi_j = bf[var]->phi[j];
	      ktype = 0;
	      gnn   = ei->gnn_list[eqn][j];
	      nvdof = ei->Baby_Dolphin[eqn][j];
	      ledof =  ei->lvdof_to_ledof[eqn][j];
	      vj    = Index_Solution(gnn, eqn, ktype , nvdof,
				     ei->matID_ledof[ledof]);
	      augc[iAC].d_evol_dx[vj] += phi_j * det_J * wt * h3;
		}
	    }
	
	}
  
  return(status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_curvature--
 *   
 *   Assemble the residual and jacobian contributions of the
 *   curvature of the level set function:
 *           H = div( grad_F/|grad_F|)
 *
 *    
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Feb 12, 2002 tabaer@sandia.gov
 *
 * Revised:	
 *
 *
 *
 */

int
assemble_curvature ( ) /*  time step size      */
{
  int i,j,p,a;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status=0;

  double det_J;
  double h3;			/* Volume element (scale factors). */
  double wt_func;
  double wt;

  double source, diffusion = 0.0;
  


  dim   = pd->Num_Dim;
  
 
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[eqn = R_CURVATURE] )
    {
      return(status);
    }

  /* only relevant if level set field exists */

  if ( ls == NULL ) EH(-1," Curvature equation can only be used in conjunction with level set.");

  peqn = upd->ep[eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;
  
  load_lsi( ls->Length_Scale );
  load_lsi_derivs();

  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      if ( pd->e[eqn] )
	{
	  for( i=0; i<ei->dof[eqn]; i++)
	    {
	      
	      wt_func = bf[eqn]->phi[i];

	      
	      diffusion = 0.;
	      /* Now granted this really isn't a diffusion term
	       * but it does involve a divergence operator integrated by parts
	       */
	      if( pd->e[eqn] & T_DIFFUSION )
		{
		  
		  for(a=0; a<dim; a++)
		    diffusion += bf[eqn]->grad_phi[i][a]*lsi->normal[a];

		  diffusion *=  det_J * wt * h3;
		  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		}

	      source = 0.0;

	      if( pd->e[eqn] & T_SOURCE )
		{
		  
 		  source += fv->H; 

		  source *= wt_func * det_J * h3 * wt;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];
		}

	      lec->R[peqn][i] += source + diffusion;  
	    }
	}
    }
  /*
   * Jacobian terms_________________________________________________________________
   */
  
  if ( af->Assemble_Jacobian )
    {


      
      if( pd->e[eqn] )
	{
	  for( i=0; i<ei->dof[eqn]; i++ )
	    {
	      wt_func = bf[eqn]->phi[i];

	      /*
	       * J_H_H
	       */

	      var = CURVATURE;
	      
	      if( pd->v[var] )
		{
		  pvar = upd->vp[var];
		    
		  for( j=0; j<ei->dof[var]; j++ )
		    {
		      source = 0.0;

		      if ( pd->e[eqn] & T_SOURCE )
			{
			  
 			  source += bf[var]->phi[j]; 
			  source *= wt_func*det_J*wt*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += source;

		    }
		}
	      /* 
	       * J_H_mesh
	       */


	      for( a=0 ; a<dim ; a++ )
		{
		  var = MESH_DISPLACEMENT1 + a;
	      
		  if( pd->v[var] )
		    {
		      pvar = upd->vp[var];

		      for( j=0; j<ei->dof[var]; j++ )
			{
			  source = 0.0;

			  if( pd->e[eqn] & T_SOURCE )
			    {
			      
			      source += h3*bf[var]->d_det_J_dm[a][j] + det_J*fv->dh3dmesh[a][j];
			      source *= wt_func*fv->H*wt;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
			  
			  diffusion = 0.0;

			  if( pd->e[eqn] & T_DIFFUSION )
			    {
			      double diff1 = 0.0, diff2 = 0.0, diff3= 0.0;

			      for( p=0; p<dim; p++) 
				diff1 += bf[eqn]->d_grad_phi_dmesh[i][p][a][j]*lsi->normal[p];
				
			      diff1 *= wt*h3*det_J;

			      for (p=0 ; p<dim ; p++ )
				diff2 += bf[eqn]->grad_phi[i][p]*lsi->d_normal_dmesh[p][a][j];
				
			      diff2 *= wt*h3*det_J;

			      for( p=0; p<dim; p++ )
				diff3 += bf[eqn]->grad_phi[p][i]*lsi->normal[p];
			      
			      diff3 *= h3*bf[var]->d_det_J_dm[a][j] + det_J*fv->dh3dmesh[a][j];
			      diff3 *= wt;

			      diffusion += diff1 + diff2 + diff3;

			      diffusion *=  pd->etm[eqn][(LOG2_DIFFUSION)];
			    }
			  
			  lec->J[peqn][pvar][i][j] += source + diffusion;
				  
			}
		    }
		}
#ifdef COUPLED_FILL
	      /* 
	       * J_H_F
	       */


	      var = LS;

	      if( pd->v[var] )
		{

		  pvar = upd->vp[var];

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      
		      if( pd->e[eqn] & T_DIFFUSION )
			{

			  diffusion = 0.0;
			  
			  for( p=0; p<dim; p++ )
			    diffusion += bf[eqn]->grad_phi[i][p]*lsi->d_normal_dF[p][j];
			
	
			  diffusion *=  det_J * wt * h3;
			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}

		      lec->J[peqn][pvar][i][j] += diffusion;

		    }

		}
#endif
	    }
	}
    }
  return(0);
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_div_normals --
 *   
 *   Assemble the residual and jacobian contributions of the
 *   curvature of the level set function via divegence of the level set normal vector field, videlicit:
 *
 *           H = div( n )
 *
 *    
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Aug 27, 2003 tabaer@sandia.gov
 *
 * Revised:	
 *
 *
 *
 */

int
assemble_div_normals (  ) /*  time step size      */
{
  int i,j,p,a;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status=0;

  double det_J;
  double h3;			/* Volume element (scale factors). */
  double wt_func;
  double wt;

  double source, diffusion;

  dim   = pd->Num_Dim;
  
 
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[eqn = R_CURVATURE] )
    {
      return(status);
    }

  /* only relevant if level set field exists */

  if ( ls == NULL ) EH(-1," Curvature equation can only be used in conjunction with level set.");

  peqn = upd->ep[eqn];

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      if ( pd->e[eqn] )
	{
	  for( i=0; i<ei->dof[eqn]; i++)
	    {
	      wt_func = bf[eqn]->phi[i];
	      
	      diffusion = 0.;

	      /* Now granted this really isn't a diffusion term
	       * but it does involve a divergence operator 
	       */
	      if( pd->e[eqn] & T_DIFFUSION )
		{
		  diffusion = fv->div_n;
		  diffusion *= wt_func*det_J * wt * h3;
		  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		}

	      source = 0.0;
	      if( pd->e[eqn] & T_SOURCE )
		{
		  
 		  source += fv->H; 
		  source *= wt_func * det_J * h3 * wt;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];
		}

	      lec->R[peqn][i] += source + diffusion;  
	    }
	}
    }


  /*
   * Jacobian terms_________________________________________________________________
   */
  
  if ( af->Assemble_Jacobian )
    {
      if( pd->e[eqn] )
	{
	  for( i=0; i<ei->dof[eqn]; i++ )
	    {
	      wt_func = bf[eqn]->phi[i];

	      /*
	       * J_H_H
	       */

	      var = CURVATURE;
	      
	      if( pd->v[var] )
		{
		  pvar = upd->vp[var];
		    
		  for( j=0; j<ei->dof[var]; j++ )
		    {
		      source = 0.0;

		      if ( pd->e[eqn] & T_SOURCE )
			{
			  
 			  source += bf[var]->phi[j]; 
			  source *= wt_func*det_J*wt*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += source;

		    }
		}

	      /*
	       * J_H_mesh
	       */

	      /*
	       * J_H_n
	       */

	      for( a=0; a<dim ; a++)
		{
		  var = NORMAL1 + a;

		  if(pd->v[var])
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++ )
			{
			  diffusion = 0.0;

			  if(pd->e[eqn] & T_DIFFUSION ) {
			    
			    for(p=0; p<dim; p++) {
			      diffusion += bf[var]->grad_phi_e[j][a][p][p];
			      }

			    diffusion *= wt_func * wt*det_J*h3;
			    diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			  }
		    
			  lec->J[peqn][pvar][i][j] += diffusion;   
			}
		    }
		}
	    }
	}
    }
  return(0);
}
		      



/******************************************************************************/
/* assemble_LSvelocity 
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *
 * out:
 *	augc.lsvel      -- contribution to velocity of phase SIGN
 *	augc.d_lsvel_dx -- contribution to derivative of velocity wrt. v, H & x
 *	augc.d_lsvol_dx -- contribution to derivative of volume wrt. v, H & x
 *
 * Created:	Anne Grillet  2/4/2002
 *
 *
 */

int
assemble_LSvelocity( bool owner, int ielem )
{
  int eqn;
  int var;
  int dim;
  int vel_dim;
  int b;

  int j;
  int status=0;
  int iAC=0;
  int vj;

  dbl det_J;
  dbl h3;			/* Volume element (scale factors). */
  dbl wt;
  dbl phi_j;
  dbl d_det_J_dmeshbj;
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */
  dbl d_detJh3_dmeshbj;
  dbl d_H_dmeshbj;              /* Sensitivity of H to (b,j) mesh dof. */
  dbl d_H_dFj;                  /* Sensitivity of H to j fill function dof. */

  dbl alpha;                    /* Level set length scale */
  dbl H;                        /* Level set heaviside function  */

  /*
   * Unpack variables from structures for local convenience...
   */
  eqn = pd->ShapeVar;
  dim = pd->Num_Dim;
  h3  = fv->h3;			/* Differential volume element (scales). */
  det_J = bf[eqn]->detJ;

  iAC = pd->LSVelocityIntegral;
 
  alpha = ls -> Length_Scale;


  /* Added 3/15/01 */
  if( pd_glob[0] ->CoordinateSystem != CYLINDRICAL &&
      pd_glob[0] ->CoordinateSystem != SWIRLING )
    {
      wt  = fv->wt;
    }
  else
    {
      wt = 2.0 * M_PIE * fv->wt;
    }


  vel_dim = ( augc[iAC].DIR == I_NEG_VX ? 0 :
	    ( augc[iAC].DIR == I_NEG_VY ? 1 :
	    ( augc[iAC].DIR == I_NEG_VZ ? 2 : -1 ) ) );

  if ( vel_dim == -1 )
    {
      vel_dim = ( augc[iAC].DIR == I_POS_VX ? 0 :
		( augc[iAC].DIR == I_POS_VY ? 1 :
	        ( augc[iAC].DIR == I_POS_VZ ? 2 : -1 ) ) );
    }

  if ( vel_dim == -1 )
    EH(-1, "assemble_LSVelocity(): Couldn't identify velocity component and phase.\n");

  load_lsi(alpha);
  load_lsi_derivs();

  H = augc[iAC].LSPHASE == I_POS_FILL ? lsi->H : ( 1.0 - lsi->H ) ;

  /* total integral */

  if( owner ){ 
    augc[iAC].lsvel += wt * det_J * H *h3 * fv->v[vel_dim];
    augc[iAC].lsvol += wt * det_J * H *h3;
  } 

  /* Calculate sensitivities so that we don't need to compute it numerically
     in mm_sol_nonlinear.c  */

  /* sensitivity to velocity                                          */

  for ( b=0; b<dim; b++)
    {
      var = vel_dim;
      if ( pd->v[var] )
	{
	  for( j=0; j<ei->dof[var]; j++)
	    {
	      vj = ei->gun_list[var][j];

	      augc[iAC].d_lsvel_dx[ vj ] += wt * det_J * h3 * H * bf[var]->phi[j];
	      augc[iAC].d_lsvol_dx[ vj ] += 0.0;
	    }

	}
    }
     
  /* sensitivity to mesh                                  */

  for ( b=0; b<dim; b++)
    {
      var = MESH_DISPLACEMENT1 + b;
       if ( pd->v[var] )
	{
	  for( j=0; j<ei->dof[var]; j++)
	    {
	      
	      eqn =  R_MESH1+b;

	      phi_j = bf[var]->phi[j];

	      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];
	      
	      dh3dmesh_bj = fv->dh3dq[b] * phi_j;

              d_detJh3_dmeshbj = d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj;

	      /* now, H doesn't depend on mesh (load_lsi_derivs), but leave this in 
		 just in case that changes  */

	      d_H_dmeshbj = lsi->d_H_dmesh[b][j];

	      /*grab Index_Solution for location of x. see mm_fill_ptrs.c */

	      vj = ei->gun_list[var][j];

	      augc[iAC].d_lsvel_dx[ vj ] += wt * fv->v[vel_dim] *(H*d_detJh3_dmeshbj
								 + d_H_dmeshbj*det_J*h3);
	      augc[iAC].d_lsvol_dx[ vj ] += wt *(H*d_detJh3_dmeshbj + d_H_dmeshbj*det_J*h3);
	    }

	}
    }
     
  /* sensitivity to FILL                                  */

  var = FILL;
  if ( pd->v[var] )
    {
      for( j=0; j<ei->dof[var]; j++)
	{
	  
	  d_H_dFj = lsi->d_H_dF[j];

	  /*grab Index_Solution for location of x. see mm_fill_ptrs.c */
	  
	  vj = ei->gun_list[var][j];

	  augc[iAC].d_lsvel_dx[ vj ] += wt * fv->v[vel_dim]* d_H_dFj*det_J*h3;
	  augc[iAC].d_lsvol_dx[ vj ] += wt * d_H_dFj*det_J*h3;
	}
    }
     
  
  return(status);
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_normals--
 *   
 *   Assemble the residual and jacobian contributions of the
 *     of the level set normal vector projection equation:
 *
 *           n = grad(phi);
 *
 *    
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Aug 24, 2003 tabaer@sandia.gov
 *
 * Revised:	
 *
 *
 *
 */
int
assemble_normals( )
{

  int i,j,p,a,b;
  int ii;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status=0;
  struct Basis_Functions *bfm;
 

  double det_J;
  double h3;			/* Volume element (scale factors). */
  double wt_func;
  double phi_i,phi_j;
  double wt;

  double normal[DIM];
  double P[DIM][DIM];

  double mag_grad_F;

  double advection = 0.0, advection1 = 0.0, source = 0.0;


  dim   = pd->Num_Dim;

  
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[eqn = R_NORMAL1] )
    {
      return(status);
    }

  /* only relevant if level set field exists */

  if ( ls == NULL ) {
    WH(-1," Normals to level set can only be used in conjunction with level set.");
    return(status);
  }

  peqn = upd->ep[eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;

  for( a=0; a<dim ; a++ ) normal[a] = fv->grad_F[a];
  
  mag_grad_F =  normalize_really_simple_vector( normal, dim );

  for( a=0; a<dim; a++ ) {
    for( b=0; b<dim; b++) {
      P[a][b] = (double) delta(a,b) - normal[a]*normal[b];
    }
  }
				 
  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      for( a=0; a<dim; a++)
	{
	  eqn = R_NORMAL1+a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];
	  
	  if ( pd->e[eqn] )
	    {
	      for( i=0; i<ei->dof[eqn]; i++)
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  wt_func = bfm->phi[i];

		  source = 0.0;
		  
		  /*
		   * Not really a source per se, but trying to follow the GOMA paradigm, such as it is.
		   */ 
	      
		  if( pd->e[eqn] & T_SOURCE )
		    {
		      source = fv->n[a];
		      source *= wt*wt_func*h3*det_J;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  /*
		   * Ditto
		   */ 

		  advection = 0.0;

		  if( pd->e[eqn] & T_ADVECTION )
		    {
		      advection = -normal[a];
		      advection *= wt*wt_func*h3*det_J;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  lec->R[peqn][ii] += source + advection;
		}
	    }
	}
    }
  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      for (a = 0; a < dim; a++)
	{

	  int ledof;

	  eqn  = R_NORMAL1 + a;
	  peqn = upd->ep[eqn];
	  bfm  = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      if (ei->active_interp_ledof[ledof]) 
		{		
		  ii = ei->lvdof_to_row_lvdof[eqn][i];
		  
		  phi_i = bfm->phi[i];

		  /*
		   *  J_n_n
		   */

		  for( b=0; b< dim; b++)
		    {
		      var = NORMAL1 + b;
		      pvar = upd->vp[var];

		      if( eqn == var )
			{
			  for( j=0; j<ei->dof[var]; j++ )
			    {
			      phi_j = bf[var]->phi[j];
			      
			      if( pd->e[eqn] & T_SOURCE )
				{
				  source = 0.0;
				  source = phi_j;

				  source *= phi_i*wt*det_J*h3;
				  source *= pd->etm[eqn][(LOG2_SOURCE)];
				}

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}
		    }

		  /*
		   *  J_n_F
		   */
#ifdef COUPLED_FILL
		  var = LS;
		  pvar = upd->vp[var];

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      advection = 0.0;

		      if( pd->e[eqn] & T_ADVECTION )
			{
			  for( b=0; b< dim ; b++)
			    {
			      advection  += P[a][b]*bf[var]->grad_phi[j][b]/mag_grad_F;
			    }

			  advection *= -phi_i*wt*h3*det_J;
			  advection *= pd->etm[eqn][(LOG2_ADVECTION)];	
		    
			}
		      lec->J[peqn][pvar][ii][j] += advection;
		    }
#endif  

		  /*
		   * J_n_d
		   */

		  for( b=0; b<dim ; b++)
		    {
		      var = R_MESH1 + b;
		      pvar = upd->vp[var];

		      for( j=0; j<ei->dof[var]; j++)
			{
			  advection = 0.0;

			  if( pd->e[eqn] & T_ADVECTION )
			    {

			      for( p=0; p<dim; p++ )
				{
				  advection += P[a][p]*fv->d_grad_F_dmesh[p][b][j]/mag_grad_F;
				}

			      advection *= -phi_i*wt*det_J*h3;
			      advection *= pd->etm[eqn][(LOG2_ADVECTION)];

			      advection1 = -normal[a]*( bf[var]->d_det_J_dm[b][j]*h3 + fv->dh3dmesh[b][j]*det_J);
			      advection1 *= phi_i*wt;
			      advection1 *= pd->etm[eqn][(LOG2_ADVECTION)];

			    }
			  lec->J[peqn][pvar][ii][j] += advection + advection1;
			}
		    }
		}
	    }
	}
    }
  return(status);
}

int
assemble_ls_momentum_source(void)
{

  int i, j, p, q, a;
  int ii;
  int peqn, pvar;
  int var;

  int dim,wim;
  int eqn;
  int status=0;
  struct Basis_Functions *bfm;
  double (* grad_phi_i_e_a ) [DIM] = NULL;
 

  double det_J;
  double h3;			/* Volume element (scale factors). */
  double phi_i,phi_j;
  double wt;

  double csf[DIM][DIM];
  double d_csf_dF[DIM][DIM][MDE];

  double source = 0.0;


  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;
    
  /*
   * Bail out fast if there's nothing to do...
   */
  
  if ( ! pd->e[eqn = R_LAGR_MULT1] )
    {
      return(status);
    }

  /* only relevant if level set field exists */

  if ( ls == NULL ) {
    WH(-1," Normals to level set can only be used in conjunction with level set.");
    return(status);
  }

  peqn = upd->ep[eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;

  memset( csf,          0, sizeof(double)*DIM*DIM);
  memset( d_csf_dF,     0, sizeof(double)*DIM*DIM*MDE);

  /* Fetch the level set interface functions. */
  load_lsi(ls->Length_Scale);

  /* Calculate the CSF tensor. */
  for ( p=0; p < VIM; p++)
    {
      for ( q=0; q < VIM; q++)
	{
 	  csf[p][q] = mp->surface_tension *
	    ( (double)delta(p,q) - lsi->normal[p] * lsi->normal[q] );
	}
    }

#ifdef COUPLED_FILL
  load_lsi_derivs();

  /* Calculate the derivatives. */
  var = FILL;
  for ( j=0 ; j < ei->dof[var] ; j++ )
    {
      for ( p=0 ; p < VIM ; p++ )
	{
	  for ( q=0; q < VIM; q++ )
	    {
	      d_csf_dF[p][q][j] = 
		- mp->surface_tension * lsi->d_normal_dF[p][j] * lsi->normal[q]
		- mp->surface_tension * lsi->normal[p] * lsi->d_normal_dF[q][j];
	    }
	}
    }
#endif

  /*
   * Residuals_________________________________________________________________
   */
  
  if ( af->Assemble_Residual )
    {
      for( a=0; a<dim; a++)
	{
	  eqn = R_LAGR_MULT1+a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];
	  
	  if ( pd->e[eqn] )
	    {
	      for( i=0; i<ei->dof[eqn]; i++)
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];
	      
		  if( pd->e[eqn] )
		    {
		      source = fv->lm[a]*phi_i;

                       for ( p=0; p<VIM; p++)
		        {
		          for ( q=0; q<VIM; q++)
			    {
			      source += grad_phi_i_e_a[p][q] * csf[q][p]; 
			    }
			}
		      source *= wt*h3*det_J;
		      
		      lec->R[peqn][ii] += source;
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
      for (a = 0; a < dim; a++)
	{

	  int ledof;

	  eqn  = R_LAGR_MULT1 + a;
	  peqn = upd->ep[eqn];
	  bfm  = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      if (ei->active_interp_ledof[ledof]) 
		{		
		  ii = ei->lvdof_to_row_lvdof[eqn][i];
		  
		  phi_i = bfm->phi[i];
		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];

		  /*
		   *  J_n_n
		   */
                  var = eqn;
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++ )
		    {
		      phi_j = bf[var]->phi[j];
			      
		      if( pd->e[eqn] )
			{
			  source = phi_j*phi_i;

			  source *= wt*det_J*h3;
			}

		      lec->J[peqn][pvar][ii][j] += source;
		    }

		  /*
		   *  J_n_F
		   */
#ifdef COUPLED_FILL
		  var = LS;
		  pvar = upd->vp[var];

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      source = 0.0;

		      if( pd->e[eqn] )
			{
			  for ( p=0; p<VIM; p++)
		            {
		              for ( q=0; q<VIM; q++)
			        {
			          source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j]; 
			        }
			    }
			  source *= wt*det_J*h3;
		    
			}
		      lec->J[peqn][pvar][ii][j] += source;
		    }
#endif  
		}
	    }
	}
    }

  return(status);
}

int
apply_ls_momentum_source(void)
{
  int i,j,a, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;

  eqn = R_MOMENTUM1;	
  if ( ! pd->e[eqn] )
    {
      return(0);
    }
    
  wt = fv->wt;	
  h3 = fv->h3;	
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;


#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  
		  source = fv->lm[a] * lsi->delta;
		    
		  source *= phi_i * det_J * wt * h3;

		  /* J_m_F
		   */
                  var = LS ;
		  pvar = upd->vp[var];
		      
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  
		  source = fv->lm[a] * lsi->delta;
		    
		  source *= phi_i * det_J * wt * h3;

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }


  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];


		  /* J_m_n
		   */
		  var = LAGR_MULT1 + a;
		  if( pd->v[var] )
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  
			  source = phi_j * lsi->delta;

			  source *= phi_i * det_J * wt * h3;

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
		  /* J_m_F
		   */
		  var = LS;
		  if( pd->v[var] )
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  
			  source = fv->lm[a] * lsi->d_delta_dF[j];

			  source *= phi_i * det_J * wt * h3;

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
		}
	    }
	}
    }
	  
  return ( 1 );
}

/***************************************************************************/
/****************************************************************************/
/****************************************************************************/



void 
scalar_fv_fill(double **base, double **base_dot, double **base_old,
	       double *phiv, int dofs,
	       double *val, double *val_dot, double *val_old
	       )
/*
 *   scalar_fv_fill:   Kernal operation for calculating the value of a
 *                     scalar field variable at the quadrature point from
 *                     the basis functions.
 *                     Note for the optimal results, this function should
 *                     be inlined during the compile.
 */
{
  int i;
  *val = *val_dot = *val_old = 0.0;
  if (pd->TimeIntegration == STEADY )  {
    for (i = 0; i < dofs; i++)  {
      *val += *(base[i]) * phiv[i];
    }
  } else { 
    for (i = 0; i < dofs; i++)  {
      *val     += *(base[i])     * phiv[i];
      *val_dot += *(base_dot[i]) * phiv[i];
      *val_old += *(base_old[i]) * phiv[i];
    }
  }
}


/***************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
grad_scalar_fv_fill(double **base, double (*grad_phiv)[DIM], int dofs,
			   double *grad_val )
/*
*/
{
	int i;

	grad_val[0] = grad_val[1] = grad_val[2]=0.0;
	for (i = 0; i < dofs; i++)  {
		grad_val[0] += *(base[i]) * grad_phiv[i][0];
		grad_val[1] += *(base[i]) * grad_phiv[i][1];
		if( VIM == 3) grad_val[2] += *(base[i]) * grad_phiv[i][2];
	}

	return;
}

void
grad_vector_fv_fill ( double ***base, double (*grad_phiv)[DIM][DIM][DIM], int dofs, double ( *grad_val )[DIM] )
{
	int r, i;
	
	double base_off;
	double (*grad_phiv_off)[DIM];
	
	memset( grad_val, 0, DIM*DIM*sizeof(double) );
	
	for ( r=0; r<VIM; r++)
	{
		for ( i=0; i<dofs; i++)
		{
			base_off = *base[r][i];
			grad_phiv_off = grad_phiv[i][r];
			
			grad_val[0][0] += base_off * grad_phiv_off[0][0];	      
			grad_val[1][1] += base_off * grad_phiv_off[1][1];	      
			grad_val[0][1] += base_off * grad_phiv_off[0][1];	      
			grad_val[1][0] += base_off * grad_phiv_off[1][0];	
			
			if( VIM == 3 )
			{
				grad_val[2][2] += base_off * grad_phiv_off[2][2];	      
				grad_val[2][1] += base_off * grad_phiv_off[2][1];	      
				grad_val[2][0] += base_off * grad_phiv_off[2][0];	      
				grad_val[1][2] += base_off * grad_phiv_off[1][2];	      
				grad_val[0][2] += base_off * grad_phiv_off[0][2];	      

			}
			
		/*	for ( p=0; p<VIM; p++)
			{
				for ( q=0; q<VIM; q++)
				{
					(grad_val[p][q] += *base[r][i] * grad_phiv[i][r] [p][q];)	      
					grad_val[p][q] += base_off * grad_phiv_off[p][q];	      
				}
			}*/
		}
	}
	return;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int 
load_fv(void)

     /*******************************************************************************
      * load_fv() -- load up values of all relevant field variables at the
      *              current gauss pt
      *
      * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
      * ------
      *	
      *
      * output: ( this routine fills in the following parts of the fv structure)
      * -------
      *		fv->
      *			T -- temperature	(scalar)
      *			v -- velocity		(vector)
      *			d -- mesh displacement	(vector)
      *			c -- concentration	(multiple scalars)
      *		      por -- porous media       (multiple scalars)
      *			P -- pressure		(scalar)
      *			S -- polymer stress	(tensor)
      *			G -- velocity gradient	(tensor)
      *                     pv -- particle velocity  (vector)
      *                     pG -- particle velocity gradient (tensor)
      *              mp->StateVector[] is filled in as well, for 
      *                    pertinent entries that make up the specification of 
      *                    the state of the material.
      *
      * NOTE: To accommodate shell elements, this function has been modified
      *       so that fv variables are not zeroed out when they are active
      *       on an element block other than the current one.
      *       The check done for variable v is then:
      *          if ( pd->v[v] || upd->vp[v] == -1 )
      *       In many cases below, this conditional zeroing is done in
      *       a separate small loop before the main one.
      *
      *
      * Return values:
      *		0 -- if things went OK
      *	       -1 -- if a problem occurred
      *
      * Created:	Fri Mar 18 06:44:41 MST 1994 pasacki@sandia.gov
      *
      * Modified:	
      ***************************************************************************/
{
  int v;			/* variable type indicator */
  int i;			/* index */
  int p, q;			/* dimension indeces */
  int dim;
  int velodim;			/* Someday...we might have more velocity */
				/* components than we have spatial dimensions */
				/* these are the 2.5 dimensional problems.*/
  int dofs;			/* degrees of freedom for a var in the elem */
  int w;			/* concentration species and porous media vars counter */
  int node, index;
  int N;
  int mode;
  int status = 0;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];
  double rho, *stateVector = mp->StateVector;
  BASIS_FUNCTIONS_STRUCT *bfv;
  int transient_run = pd->TimeIntegration != STEADY ;
  int *pdv = pd->v;

  status = 0;


  /* load eqn and variable number in tensor form */
  if( pdv[POLYMER_STRESS11] ) {
    status = stress_eqn_pointer(v_s);
    EH(status, "stress_eqn_pointer(v_s)");
  }
  if( pdv[VELOCITY_GRADIENT11] ) {
    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32; 
    v_g[2][2] = VELOCITY_GRADIENT33; 
  }
  
  /*
   * Since it is possible to have a 1D element in a 2D problem, 
   * (i.e. shell elements), use the problem dimension instead.
   */
  dim = pd->Num_Dim;

  /*
   * Temperature...
   *    HKM -> Introduced a function to handle the core scalar fill
   *           operation; it should be inlined for maximum speed.
   *           Added the concept that the reference temperature should
   *           be used, when no temperature variable exists in the solution
   *           vector.
   */
  
  if (pdv[TEMPERATURE]) {
    v = TEMPERATURE;
    scalar_fv_fill(esp->T, esp_dot->T, esp_old->T, bf[v]->phi, ei->dof[v],
                   &(fv->T), &(fv_dot->T), &(fv_old->T));
    stateVector[TEMPERATURE] = fv->T;
  } /*else if (upd->vp[v] == -1) {
      fv->T     = mp->reference[TEMPERATURE];
      fv_old->T = 0.;
      fv_dot->T = 0.;
      } */
  

  /*
   * Fill...
   */
  
  if (pdv[FILL]) {
    v = FILL;
    scalar_fv_fill(esp->F, esp_dot->F, esp_old->F, bf[v]->phi, ei->dof[v],
                   &(fv->F), &(fv_dot->F), &(fv_old->F));
    stateVector[v] = fv->F;
  } /*else if (upd->vp[v] == -1) {
      fv->F = fv_old->F = fv_dot->F = 0.;
      } */


  /*
   * Suspension Temperature...
   */
  
  if (pdv[BOND_EVOLUTION]) {
    v = BOND_EVOLUTION;
    scalar_fv_fill(esp->nn, esp_dot->nn, esp_old->nn, bf[v]->phi, ei->dof[v],
                   &(fv->nn), &(fv_dot->nn), &(fv_old->nn));
    stateVector[v] = fv->nn;
  } /*else if (upd->vp[v] == -1)  {
      fv->nn = fv_old->nn = fv_dot->nn = 0.;
      } */


  /*
   * Voltage...
   */
  
  if ( pdv[VOLTAGE] ) {
    v = VOLTAGE;
    scalar_fv_fill(esp->V, esp_dot->V, esp_old->V, bf[v]->phi, ei->dof[v],
                   &(fv->V), &(fv_dot->V), &(fv_old->V));
    stateVector[VOLTAGE] = fv->V;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->V = fv_old->V = fv_dot->V = 0.;
      }*/

  /*
   * Surface charge density...
   */
  
  if ( pdv[SURF_CHARGE] ) {
    v = SURF_CHARGE;
    scalar_fv_fill(esp->qs, esp_dot->qs, esp_old->qs, bf[v]->phi, ei->dof[v],
                   &(fv->qs), &(fv_dot->qs), &(fv_old->qs));
    stateVector[SURF_CHARGE] = fv->qs;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->qs = fv_old->qs = fv_dot->qs = 0.;
      } */

  /*
   * Structural shell curvature
   */
  
  if ( pdv[SHELL_CURVATURE] ) {
    v = SHELL_CURVATURE;
    scalar_fv_fill(esp->sh_K, esp_dot->sh_K, esp_old->sh_K, bf[v]->phi, ei->dof[v],
                   &(fv->sh_K), &(fv_dot->sh_K), &(fv_old->sh_K));
    stateVector[SHELL_CURVATURE] = fv->sh_K;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->sh_K = fv_old->sh_K = fv_dot->sh_K = 0.;
      } */

  if ( pdv[SHELL_CURVATURE2] ) {
    v = SHELL_CURVATURE2;
    scalar_fv_fill(esp->sh_K2, esp_dot->sh_K2, esp_old->sh_K2, bf[v]->phi, ei->dof[v],
                   &(fv->sh_K2), &(fv_dot->sh_K2), &(fv_old->sh_K2));
    stateVector[SHELL_CURVATURE2] = fv->sh_K2;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->sh_K2 = fv_old->sh_K2 = fv_dot->sh_K2 = 0.;
      } */  

  /*
   * Structural shell tension
   */
  
  if ( pdv[SHELL_TENSION] ) {
    v = SHELL_TENSION;
    scalar_fv_fill(esp->sh_tens, esp_dot->sh_tens, esp_old->sh_tens, bf[v]->phi, ei->dof[v],
                   &(fv->sh_tens), &(fv_dot->sh_tens), &(fv_old->sh_tens));
    stateVector[SHELL_TENSION] = fv->sh_tens;
  } /*else if ( updvp[v] == -1 ) {
      fv->sh_tens = fv_old->sh_tens = fv_dot->sh_tens = 0.;
      }*/
  /*
   * Structural shell x coordinate
   */
  
  if ( pdv[SHELL_X] ) {
    v = SHELL_X      ;
    scalar_fv_fill(esp->sh_x, esp_dot->sh_x, esp_old->sh_x, bf[v]->phi, ei->dof[v],
                   &(fv->sh_x), &(fv_dot->sh_x), &(fv_old->sh_x));
    stateVector[SHELL_X] = fv->sh_x;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->sh_x = fv_old->sh_x = fv_dot->sh_x    = 0.;
      } */
  
  /*
   * Structural shell y coordinate
   */
  
  if ( pdv[SHELL_Y] ) {
    v = SHELL_Y      ;
    scalar_fv_fill(esp->sh_y, esp_dot->sh_y, esp_old->sh_y, bf[v]->phi, ei->dof[v],
                   &(fv->sh_y), &(fv_dot->sh_y), &(fv_old->sh_y));
    stateVector[SHELL_Y] = fv->sh_y;
  } /* else if ( upd->vp[v] == -1 ) {
       fv->sh_y = fv_old->sh_y = fv_dot->sh_y    = 0.;
       } */
  
  /*
   * Shell user
   */
  
  if ( pdv[SHELL_USER] ) {
    v = SHELL_USER   ;
    scalar_fv_fill(esp->sh_u, esp_dot->sh_u, esp_old->sh_u, bf[v]->phi, ei->dof[v],
                   &(fv->sh_u), &(fv_dot->sh_u), &(fv_old->sh_u));
    stateVector[SHELL_USER] = fv->sh_u;
  } /* else if ( upd->vp[v] == -1 ) {
       fv->sh_u = fv_old->sh_u = fv_dot->sh_u    = 0.;
       } */
 
  /*
   * Shear rate from second invariant of rate of strain tensor
   */

  

  if ( pdv[SHEAR_RATE] )
    {
      v = SHEAR_RATE;
      fv->SH = 0.0;
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv->SH += *esp->SH[i] * bf[v]->phi[i];
	}
    }
  /*  else if ( upd->vp[v] == -1)
      {
      fv->SH = 0.0;
      } */

  /* Square of the norm of the potential field, |E|^2 */
  

  if(pdv[ENORM])
    {
      v = ENORM;
      fv->Enorm = 0.0;
      fv_old->Enorm = 0.0;
      dofs = ei->dof[v];
      for(i = 0; i < dofs; i++)
	{
	  fv->Enorm += *esp->Enorm[i] * bf[v]->phi[i];
	  fv_old->Enorm += *esp_old->Enorm[i] * bf[v]->phi[i];
	}
    } /*  else if ( upd->vp[v] == -1) {
	  fv->Enorm = 0.0;
	  fv_old->Enorm = 0.0;
	  } */


  /* 
   * Curvature of level set function
   */

  

  if ( pdv[CURVATURE] )
    {
      v = CURVATURE;
      scalar_fv_fill(esp->H, esp_dot->H, esp_old->H, bf[v]->phi, ei->dof[v],
		     &(fv->H), &(fv_dot->H), &(fv_old->H));

    }
  /*  else if ( upd->vp[v] == -1)
      {
      fv->H = fv_dot->H = fv_old->H = 0.0;
      }
  */
  
  /* 
   *  Normal to level set function
   */

  
  if ( pdv[NORMAL1] )
    {
      v = NORMAL1;
      scalar_fv_fill(esp->n[0], esp_dot->n[0], esp_old->n[0], bf[v]->phi, ei->dof[v],
		     &(fv->n[0]), &(fv_dot->n[0]), &(fv_old->n[0]));

    }
  /*  else if ( upd->vp[v] == -1)
      {
      fv->n[0] = fv_dot->n[0] = fv_old->n[0] = 0.0;
      } */

  
  if ( pdv[NORMAL2] )
    {
      v = NORMAL2;
      scalar_fv_fill(esp->n[1], esp_dot->n[1], esp_old->n[1], bf[v]->phi, ei->dof[v],
		     &(fv->n[1]), &(fv_dot->n[1]), &(fv_old->n[1]));

    }
  /* else if ( upd->vp[v] == -1)
     {
     fv->n[1] = fv_dot->n[1] = fv_old->n[1] = 0.0;
     }*/


  
  if ( pdv[NORMAL3] )
    {
      v = NORMAL3;
      scalar_fv_fill(esp->n[2], esp_dot->n[2], esp_old->n[2], bf[v]->phi, ei->dof[v],
		     &(fv->n[2]), &(fv_dot->n[2]), &(fv_old->n[2]));

    }
  /*  else if ( upd->vp[v] == -1)
      {
      fv->n[2] = fv_dot->n[2] = fv_old->n[2] = 0.0;
      } */
    
  /* 
   *  shell element orientation angles
   */

  
  if ( pdv[SHELL_ANGLE1] )
    {
      v = SHELL_ANGLE1;
      fv->sh_ang[0] = 0.;
      for (i = 0; i < ei->dof[v]; i++) {
        if ( ( *esp->sh_ang[0][i] - *esp->sh_ang[0][0] ) > M_PIE ) 
          fv->sh_ang[0] += bf[v]->phi[i] * ( *esp->sh_ang[0][i] - 2.*M_PIE );
        else if ( ( *esp->sh_ang[0][i] - *esp->sh_ang[0][0] ) < -M_PIE ) 
          fv->sh_ang[0] += bf[v]->phi[i] * ( *esp->sh_ang[0][i] + 2.*M_PIE );
        else 
          fv->sh_ang[0] += bf[v]->phi[i] * *esp->sh_ang[0][i];
      }
      /* surely no one will be using these */
      fv_dot->sh_ang[0] = fv_old->sh_ang[0] = 0.0;
    }/*
       else if ( upd->vp[v] == -1)
       {
       fv->sh_ang[0] = fv_dot->sh_ang[0] = fv_old->sh_ang[0] = 0.0;
       }*/

  
  if ( pdv[SHELL_ANGLE2] )
    {
      v = SHELL_ANGLE2;
      fv->sh_ang[1] = 0.;
      for (i = 0; i < ei->dof[v]; i++) {
        if ( ( *esp->sh_ang[1][i] - *esp->sh_ang[1][0] ) > M_PIE ) 
	  fv->sh_ang[1] += bf[v]->phi[i] * ( *esp->sh_ang[1][i] - 2.*M_PIE );
        else if ( ( *esp->sh_ang[1][i] - *esp->sh_ang[1][0] ) < -M_PIE ) 
          fv->sh_ang[1] += bf[v]->phi[i] * ( *esp->sh_ang[1][i] + 2.*M_PIE );
        else 
	  fv->sh_ang[1] += bf[v]->phi[i] * *esp->sh_ang[1][i];
      }
      /* surely no one will be using these */
      fv_dot->sh_ang[1] = fv_old->sh_ang[1] = 0.0;
    }
  /*  else if ( upd->vp[v] == -1)
      {
      fv->sh_ang[1] = fv_dot->sh_ang[1] = fv_old->sh_ang[1] = 0.0;
      }*/

  /*
   * Surface Rheo shell piece
   */
  
  if (pdv[SHELL_SURF_DIV_V]) {
    v = SHELL_SURF_DIV_V     ;
    scalar_fv_fill(esp->div_s_v, esp_dot->div_s_v, esp_old->div_s_v, bf[v]->phi, ei->dof[v],
                   &(fv->div_s_v), &(fv_dot->div_s_v), &(fv_old->div_s_v));
    stateVector[SHELL_SURF_DIV_V] = fv->div_s_v;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->div_s_v = fv_old->div_s_v = fv_dot->div_s_v    = 0.;
      }*/

  
  if ( pdv[SHELL_SURF_CURV] ) {
    v = SHELL_SURF_CURV    ;
    scalar_fv_fill(esp->curv, esp_dot->curv, esp_old->curv, bf[v]->phi, ei->dof[v],
                   &(fv->curv), &(fv_dot->curv), &(fv_old->curv));
    stateVector[SHELL_SURF_CURV] = fv->curv;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->curv = fv_old->curv = fv_dot->curv    = 0.;
      }*/


  if ( pdv[N_DOT_CURL_V] ) {
    v = N_DOT_CURL_V    ;
    scalar_fv_fill(esp->n_dot_curl_s_v, esp_dot->n_dot_curl_s_v, esp_old->n_dot_curl_s_v, bf[v]->phi, ei->dof[v],
                   &(fv->n_dot_curl_s_v), &(fv_dot->n_dot_curl_s_v), &(fv_old->n_dot_curl_s_v));
    stateVector[N_DOT_CURL_V] = fv->n_dot_curl_s_v;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->n_dot_curl_s_v = fv_old->n_dot_curl_s_v = fv_dot->n_dot_curl_s_v    = 0.;
      }*/


  if ( pdv[GRAD_S_V_DOT_N1] ) {
    v = GRAD_S_V_DOT_N1    ;
    scalar_fv_fill(esp->grad_v_dot_n[0], esp_dot->grad_v_dot_n[0], esp_old->grad_v_dot_n[0], bf[v]->phi, ei->dof[v],
                   &(fv->grad_v_dot_n[0]), &(fv_dot->grad_v_dot_n[0]), &(fv_old->grad_v_dot_n[0]));
    stateVector[GRAD_S_V_DOT_N1] = fv->grad_v_dot_n[0];
  }/* else if ( upd->vp[v] == -1 ) {
      fv->grad_v_dot_n[0] = fv_old->grad_v_dot_n[0] = fv_dot->grad_v_dot_n[0]    = 0.;
      }*/


  if ( pdv[GRAD_S_V_DOT_N2] ) {
    v = GRAD_S_V_DOT_N2    ;
    scalar_fv_fill(esp->grad_v_dot_n[1], esp_dot->grad_v_dot_n[1], esp_old->grad_v_dot_n[1], bf[v]->phi, ei->dof[v],
                   &(fv->grad_v_dot_n[1]), &(fv_dot->grad_v_dot_n[1]), &(fv_old->grad_v_dot_n[1]));
    stateVector[GRAD_S_V_DOT_N2] = fv->grad_v_dot_n[1];
  }/* else if ( upd->vp[v] == -1 ) {
      fv->grad_v_dot_n[1] = fv_old->grad_v_dot_n[1] = fv_dot->grad_v_dot_n[1]    = 0.;


      } */
  
  if ( pdv[GRAD_S_V_DOT_N3] ) {
    v = GRAD_S_V_DOT_N3    ;
    scalar_fv_fill(esp->grad_v_dot_n[2], esp_dot->grad_v_dot_n[2], esp_old->grad_v_dot_n[2], bf[v]->phi, ei->dof[v],
                   &(fv->grad_v_dot_n[2]), &(fv_dot->grad_v_dot_n[2]), &(fv_old->grad_v_dot_n[2]));
    stateVector[GRAD_S_V_DOT_N3] = fv->grad_v_dot_n[2];
  } /*else if ( upd->vp[v] == -1 ) {
      fv->grad_v_dot_n[2] = fv_old->grad_v_dot_n[2] = fv_dot->grad_v_dot_n[2]    = 0.;
      }*/


  if ( pdv[SHELL_DIFF_FLUX] ) {
    v = SHELL_DIFF_FLUX  ;
    scalar_fv_fill(esp->sh_J, esp_dot->sh_J, esp_old->sh_J, bf[v]->phi, ei->dof[v],
                   &(fv->sh_J), &(fv_dot->sh_J), &(fv_old->sh_J));
    stateVector[SHELL_DIFF_FLUX] = fv->sh_J;
  } /*else if ( upd->vp[v] == -1 ) {
      fv->sh_J = fv_old->sh_J = fv_dot->sh_J = 0.;
      }*/
                                                                                 

  if ( pdv[SHELL_DIFF_CURVATURE] ) {
    v = SHELL_DIFF_CURVATURE;
    scalar_fv_fill(esp->sh_Kd, esp_dot->sh_Kd, esp_old->sh_Kd, bf[v]->phi, ei->dof[v],
                   &(fv->sh_Kd), &(fv_dot->sh_Kd), &(fv_old->sh_Kd));
    stateVector[SHELL_DIFF_CURVATURE] = fv->sh_Kd;
  }/* else if ( upd->vp[v] == -1 ) {
      fv->sh_Kd = fv_old->sh_Kd = fv_dot->sh_Kd = 0.;
      }*/
 

  if ( pdv[SHELL_NORMAL1] ) {
    v = SHELL_NORMAL1    ;
    scalar_fv_fill(esp->n[0], esp_dot->n[0], esp_old->n[0], bf[v]->phi, ei->dof[v],
                   &(fv->n[0]), &(fv_dot->n[0]), &(fv_old->n[0]));
    stateVector[SHELL_NORMAL1] = fv->n[0];
  } /*else if ( upd->vp[v] == -1 ) {
      fv->n[0] = fv_old->n[0] = fv_dot->n[0] = 0.;
      }*/
                                                                                 

  if ( pdv[SHELL_NORMAL2] ) {
    v = SHELL_NORMAL2    ;
    scalar_fv_fill(esp->n[1], esp_dot->n[1], esp_old->n[1], bf[v]->phi, ei->dof[v],
                   &(fv->n[1]), &(fv_dot->n[1]), &(fv_old->n[1]));
    stateVector[SHELL_NORMAL2] = fv->n[1];
  } /*else if ( upd->vp[v] == -1 ) {
      fv->n[1] = fv_old->n[1] = fv_dot->n[1] = 0.;
      }*/

  if ( pdv[SHELL_NORMAL3] ) {
    v = SHELL_NORMAL3    ;
    scalar_fv_fill(esp->n[2], esp_dot->n[2], esp_old->n[2], bf[v]->phi, ei->dof[v],
                   &(fv->n[2]), &(fv_dot->n[2]), &(fv_old->n[2]));
    stateVector[SHELL_NORMAL3] = fv->n[2];
  } /*else if ( upd->vp[v] == -1 ) {
      fv->n[2] = fv_old->n[2] = fv_dot->n[2] = 0.;
      }*/

  if ( pdv[SHELL_NORMAL1] && pdv[SHELL_NORMAL2] && pdv[SHELL_NORMAL3] )
    {
     memset (fv->d_n_dxi, 0.0,sizeof(double)*DIM*DIM);
     for (i = 0; i < ei->dof[SHELL_NORMAL1]; i++)
        {
         fv->d_n_dxi[0][0]+= *esp->n[0][i] * bf[SHELL_NORMAL1]->dphidxi[i][0];
         fv->d_n_dxi[1][0]+= *esp->n[1][i] * bf[SHELL_NORMAL2]->dphidxi[i][0];
         fv->d_n_dxi[2][0]+= *esp->n[2][i] * bf[SHELL_NORMAL3]->dphidxi[i][0];

         fv->d_n_dxi[0][1]+= *esp->n[0][i] * bf[SHELL_NORMAL1]->dphidxi[i][1];
         fv->d_n_dxi[1][1]+= *esp->n[1][i] * bf[SHELL_NORMAL2]->dphidxi[i][1];
         fv->d_n_dxi[2][1]+= *esp->n[2][i] * bf[SHELL_NORMAL3]->dphidxi[i][1];
        }
    }

  /*
   *	Acoustic Pressure
   */

  if (pdv[ACOUS_PREAL]) {
    v = ACOUS_PREAL;
    scalar_fv_fill(esp->apr, esp_dot->apr, esp_old->apr, bf[v]->phi, ei->dof[v],
                   &(fv->apr), &(fv_dot->apr), &(fv_old->apr));
    stateVector[ACOUS_PREAL] = fv->apr;
  } /*else if (upd->vp[v] == -1) {
      fv->apr = fv_old->apr = fv_dot->apr = 0.;
      }*/
  


  if (pdv[ACOUS_PIMAG]) {
    v = ACOUS_PIMAG;
    scalar_fv_fill(esp->api, esp_dot->api, esp_old->api, bf[v]->phi, ei->dof[v],
                   &(fv->api), &(fv_dot->api), &(fv_old->api));
    stateVector[ACOUS_PIMAG] = fv->api;
  } /*else if (upd->vp[v] == -1) {
      fv->api = fv_old->api = fv_dot->api = 0.;
      }*/
  


  if (pdv[ACOUS_REYN_STRESS]) {
    v = ACOUS_REYN_STRESS;
    scalar_fv_fill(esp->ars, esp_dot->ars, esp_old->ars, bf[v]->phi, ei->dof[v],
                   &(fv->ars), &(fv_dot->ars), &(fv_old->ars));
    stateVector[ACOUS_REYN_STRESS] = fv->ars;
  } /*else if (upd->vp[v] == -1) {
      fv->ars = fv_old->ars = fv_dot->ars = 0.;
      }*/
 

  if (pdv[SHELL_BDYVELO]) {
    v = SHELL_BDYVELO;
    scalar_fv_fill(esp->sh_bv, esp_dot->sh_bv, esp_old->sh_bv, bf[v]->phi, ei->dof[v],
                   &(fv->sh_bv), &(fv_dot->sh_bv), &(fv_old->sh_bv));
    stateVector[SHELL_BDYVELO] = fv->sh_bv;
  } /*else if (upd->vp[v] == -1) { */

  if (pdv[SHELL_LUBP]) {
    v = SHELL_LUBP;
    scalar_fv_fill(esp->sh_p, esp_dot->sh_p, esp_old->sh_p, bf[v]->phi, ei->dof[v],
                   &(fv->sh_p), &(fv_dot->sh_p), &(fv_old->sh_p));
    stateVector[SHELL_LUBP] = fv->sh_p;
  } 

  if (pdv[LUBP]) {
    v = LUBP;
    scalar_fv_fill(esp->lubp, esp_dot->lubp, esp_old->lubp, bf[v]->phi, ei->dof[v],
                   &(fv->lubp), &(fv_dot->lubp), &(fv_old->lubp));
    stateVector[LUBP] = fv->lubp;
  } 

  if (pdv[LUBP_2]) {
    v = LUBP_2;
    scalar_fv_fill(esp->lubp_2, esp_dot->lubp_2, esp_old->lubp_2, bf[v]->phi, ei->dof[v],
                   &(fv->lubp_2), &(fv_dot->lubp_2), &(fv_old->lubp_2));
    stateVector[LUBP_2] = fv->lubp_2;
  } 

 if (pdv[SHELL_FILMP]) {
    v = SHELL_FILMP;
    scalar_fv_fill(esp->sh_fp, esp_dot->sh_fp, esp_old->sh_fp, bf[v]->phi, ei->dof[v],
                   &(fv->sh_fp), &(fv_dot->sh_fp), &(fv_old->sh_fp));
    stateVector[SHELL_FILMP] = fv->sh_fp;
  } 

  if (pdv[SHELL_FILMH]) {
    v = SHELL_FILMH;
    scalar_fv_fill(esp->sh_fh, esp_dot->sh_fh, esp_old->sh_fh, bf[v]->phi, ei->dof[v],
                   &(fv->sh_fh), &(fv_dot->sh_fh), &(fv_old->sh_fh));
    stateVector[SHELL_FILMH] = fv->sh_fh;
  } 
  
  if (pdv[SHELL_PARTC]) {
    v = SHELL_PARTC;
    scalar_fv_fill(esp->sh_pc, esp_dot->sh_pc, esp_old->sh_pc, bf[v]->phi, ei->dof[v],
                   &(fv->sh_pc), &(fv_dot->sh_pc), &(fv_old->sh_pc));
    stateVector[SHELL_PARTC] = fv->sh_pc;
  } 


  if (pdv[SHELL_SAT_CLOSED]) {
    v = SHELL_SAT_CLOSED;
    scalar_fv_fill(esp->sh_sat_closed, esp_dot->sh_sat_closed, esp_old->sh_sat_closed, bf[v]->phi, ei->dof[v],
                   &(fv->sh_sat_closed), &(fv_dot->sh_sat_closed), &(fv_old->sh_sat_closed));
    stateVector[SHELL_SAT_CLOSED] = fv->sh_sat_closed;
  } 
  if (pdv[SHELL_PRESS_OPEN]) {
    v = SHELL_PRESS_OPEN;
    scalar_fv_fill(esp->sh_p_open, esp_dot->sh_p_open, esp_old->sh_p_open, bf[v]->phi, ei->dof[v],
                   &(fv->sh_p_open), &(fv_dot->sh_p_open), &(fv_old->sh_p_open));
    stateVector[SHELL_PRESS_OPEN] = fv->sh_p_open;
  } 
  if (pdv[SHELL_PRESS_OPEN_2]) {
    v = SHELL_PRESS_OPEN_2;
    scalar_fv_fill(esp->sh_p_open_2, esp_dot->sh_p_open_2, esp_old->sh_p_open_2, bf[v]->phi, ei->dof[v],
                   &(fv->sh_p_open_2), &(fv_dot->sh_p_open_2), &(fv_old->sh_p_open_2));
    stateVector[SHELL_PRESS_OPEN_2] = fv->sh_p_open_2;
  } 
  if (pdv[SHELL_TEMPERATURE]) {
    v = SHELL_TEMPERATURE;
    scalar_fv_fill(esp->sh_t, esp_dot->sh_t, esp_old->sh_t, bf[v]->phi, ei->dof[v],
                   &(fv->sh_t), &(fv_dot->sh_t), &(fv_old->sh_t));
    stateVector[SHELL_TEMPERATURE] = fv->sh_t;
  } 
  if (pdv[SHELL_DELTAH]) {
    v = SHELL_DELTAH;
    scalar_fv_fill(esp->sh_dh, esp_dot->sh_dh, esp_old->sh_dh, bf[v]->phi, ei->dof[v],
                   &(fv->sh_dh), &(fv_dot->sh_dh), &(fv_old->sh_dh));
    stateVector[SHELL_DELTAH] = fv->sh_dh;
  } 
  if (pdv[SHELL_LUB_CURV]) {
    v = SHELL_LUB_CURV;
    scalar_fv_fill(esp->sh_l_curv, esp_dot->sh_l_curv, esp_old->sh_l_curv, bf[v]->phi, ei->dof[v],
                   &(fv->sh_l_curv), &(fv_dot->sh_l_curv), &(fv_old->sh_l_curv));
    stateVector[SHELL_LUB_CURV] = fv->sh_l_curv;
  } 
  if (pdv[SHELL_LUB_CURV_2]) {
    v = SHELL_LUB_CURV_2;
    scalar_fv_fill(esp->sh_l_curv_2, esp_dot->sh_l_curv_2, esp_old->sh_l_curv_2, bf[v]->phi, ei->dof[v],
                   &(fv->sh_l_curv_2), &(fv_dot->sh_l_curv_2), &(fv_old->sh_l_curv_2));
    stateVector[SHELL_LUB_CURV] = fv->sh_l_curv;
  } 
  if (pdv[SHELL_SAT_GASN]) {
    v = SHELL_SAT_GASN;
    scalar_fv_fill(esp->sh_sat_gasn, esp_dot->sh_sat_gasn, esp_old->sh_sat_gasn, bf[v]->phi, ei->dof[v],
                   &(fv->sh_sat_gasn), &(fv_dot->sh_sat_gasn), &(fv_old->sh_sat_gasn));
    stateVector[SHELL_SAT_GASN] = fv->sh_sat_gasn;
  } 
  if (pdv[SHELL_SHEAR_TOP]) {
    v = SHELL_SHEAR_TOP;
    scalar_fv_fill(esp->sh_shear_top, esp_dot->sh_shear_top, esp_old->sh_shear_top, bf[v]->phi, ei->dof[v],
                   &(fv->sh_shear_top), &(fv_dot->sh_shear_top), &(fv_old->sh_shear_top));
    stateVector[SHELL_SHEAR_TOP] = fv->sh_shear_top;
  }
  if (pdv[SHELL_SHEAR_BOT]) {
    v = SHELL_SHEAR_BOT;
    scalar_fv_fill(esp->sh_shear_bot, esp_dot->sh_shear_bot, esp_old->sh_shear_bot, bf[v]->phi, ei->dof[v],
                   &(fv->sh_shear_bot), &(fv_dot->sh_shear_bot), &(fv_old->sh_shear_bot));
    stateVector[SHELL_SHEAR_BOT] = fv->sh_shear_bot;
  }
  if (pdv[SHELL_CROSS_SHEAR]) {
    v = SHELL_CROSS_SHEAR;
    scalar_fv_fill(esp->sh_cross_shear, esp_dot->sh_cross_shear, esp_old->sh_cross_shear, bf[v]->phi, ei->dof[v],
                   &(fv->sh_cross_shear), &(fv_dot->sh_cross_shear), &(fv_old->sh_cross_shear));
    stateVector[SHELL_CROSS_SHEAR] = fv->sh_cross_shear;
  }
  if (pdv[MAX_STRAIN]) {
    v = MAX_STRAIN;
    scalar_fv_fill(esp->max_strain, esp_dot->max_strain, esp_old->max_strain, bf[v]->phi, ei->dof[v],
                   &(fv->max_strain), &(fv_dot->max_strain), &(fv_old->max_strain));
    stateVector[MAX_STRAIN] = fv->max_strain;
  }
  if (pdv[CUR_STRAIN]) {
    v = CUR_STRAIN;
    scalar_fv_fill(esp->cur_strain, esp_dot->cur_strain, esp_old->cur_strain, bf[v]->phi, ei->dof[v],
                   &(fv->cur_strain), &(fv_dot->cur_strain), &(fv_old->cur_strain));
    stateVector[CUR_STRAIN] = fv->cur_strain;
  }
  if ( pdv[LIGHT_INTP] ) {
    v = LIGHT_INTP    ;
    scalar_fv_fill(esp->poynt[0], esp_dot->poynt[0], esp_old->poynt[0], bf[v]->phi, ei->dof[v],
                   &(fv->poynt[0]), &(fv_dot->poynt[0]), &(fv_old->poynt[0]));
    stateVector[LIGHT_INTP] = fv->poynt[0];
  } 
  if ( pdv[LIGHT_INTM] ) {
    v = LIGHT_INTM    ;
    scalar_fv_fill(esp->poynt[1], esp_dot->poynt[1], esp_old->poynt[1], bf[v]->phi, ei->dof[v],
                   &(fv->poynt[1]), &(fv_dot->poynt[1]), &(fv_old->poynt[1]));
    stateVector[LIGHT_INTM] = fv->poynt[1];
  } 
  if ( pdv[LIGHT_INTD] ) {
    v = LIGHT_INTD    ;
    scalar_fv_fill(esp->poynt[2], esp_dot->poynt[2], esp_old->poynt[2], bf[v]->phi, ei->dof[v],
                   &(fv->poynt[2]), &(fv_dot->poynt[2]), &(fv_old->poynt[2]));
    stateVector[LIGHT_INTD] = fv->poynt[2];
  } 
  /*
   *	Porous sink mass
   */
  
  if (pdv[POR_SINK_MASS]) {
    v = POR_SINK_MASS;
    scalar_fv_fill(esp->sink_mass, esp_dot->sink_mass, esp_old->sink_mass, bf[v]->phi, ei->dof[v],
                   &(fv->sink_mass), &(fv_dot->sink_mass), &(fv_old->sink_mass));
    stateVector[POR_SINK_MASS] = fv->sink_mass;
  } /*else if (upd->vp[v] == -1) {
      fv->sink_mass = fv_old->sink_mass = fv_dot->sink_mass = 0.;
      }*/
  

  /*
   * Pressure...
   */


  if (pdv[PRESSURE]) 
    {
      v = PRESSURE;
      scalar_fv_fill(esp->P, esp_dot->P, esp_old->P, bf[v]->phi, ei->dof[v],
		     &(fv->P), &(fv_dot->P), &(fv_old->P));
      stateVector[v] = fv->P + upd->Pressure_Datum;
    } 

  /*
   *  Set the state vector pressure datum to include an additional
   *  uniform pressure datum. 
   */
  
	      
  /*
   * Mesh displacement (vector)...
   * and positions (vector)
   */
  
  /*
   * Insure default value of 3rd coordinate is zero for low
   * dimensional problems...DIM=3, while dim can be less...
   */
  


#ifdef DEBUG_HKM
  if (ei->ielem == 165) {
    //printf("load_fv: we are here\n");
  }
#endif

  if (ei->deforming_mesh) {
    v = pd->ShapeVar;
    /*
     * ShapeVar will always be mesh displacement where it is defined.
     * Otherwise (e.g. for shell elements), it will be set correctly.
     */
    bfv = bf[v];
    dofs     = ei->dof[v];

    for (p = 0; p < dim; p++)
      {

	  
	fv->x0[p]        = 0.0;
	fv->x[p]         = 0.0;
	fv_old->x[p]     = 0.0;
	fv_dot->x[p]     = 0.0;
	fv_dot_old->x[p] = 0.0;

	fv->d[p]     = 0.0;
	fv_old->d[p] = 0.0;
	fv_dot->d[p] = 0.0;
	fv_dot_old->d[p] = 0.0;

	if (tran->solid_inertia)
	  {
	    fv_dot_dot->x[p] = 0;
	    fv_dot_dot->d[p] = 0;
	  }

	/*
	 * If this is a shell element, mesh displacements may not be
	 * defined on this element block even if the mesh is deforming.
	 * In this case, displacements are defined on a neighbor block
	 * and ei->deforming_mesh will be TRUE, so that the true
	 * displaced coordinates can be loaded here.
	 */
	for (i = 0; i < dofs; i++)
	  {
	    node = ei->dof_list[R_MESH1][i];
	    index = Proc_Elem_Connect[Proc_Connect_Ptr[ei->ielem] +node];
	    fv->d[p] += *esp->d[p][i] * bfv->phi[i];
	    fv->x[p] +=  ( Coor[p][index] + *esp->d[p][i] ) * bfv->phi[i];
	    fv->x0[p] +=  Coor[p][index] * bfv->phi[i];
 
	    /*
	     * HKM -> We calculate the old xdot's by doing a little pointer
	     *        arithmetic. Basically we use esp_dot as an offset to
	     *        the base and then correct for the base address.
	     */
	    if (transient_run)
	      {
		fv_old->d[p] += *esp_old->d[p][i] * bfv->phi[i];
		fv_old->x[p] +=  ( Coor[p][index] + *esp_old->d[p][i] ) *
		  bf[v]->phi[i];
		fv_dot->d[p] += *esp_dot->d[p][i] * bfv->phi[i];
		fv_dot->x[p] += *esp_dot->d[p][i] * bfv->phi[i];
		if (tran->solid_inertia)
		  {
		    fv_dot_dot->d[p] += *esp_dbl_dot->d[p][i] * bfv->phi[i];
		    fv_dot_dot->x[p] += *esp_dbl_dot->d[p][i] * bfv->phi[i];
		  }
		fv_dot_old->x[p] += *(xdot_old_static - xdot_static +
				      esp_dot->d[p][i]) * bfv->phi[i];
		fv_dot_old->d[p] += *(xdot_old_static - xdot_static +
				      esp_dot->d[p][i]) * bfv->phi[i];
		if (tran->solid_inertia)
		  {
		    fv_dot_dot_old->d[p] += *(x_dbl_dot_old_static - 
					      x_dbl_dot_static + 
					      esp_dbl_dot->d[p][i])*bfv->phi[i];
		  }
	      }
	    else
	      {	  
		fv_old->d[p] = fv->d[p];
		fv_old->x[p] = fv->x[p]; /* Fixed grid stays fixed thru time. */
	      }

	  }
      }
      /* Zero these only if not using mesh displacements: */
  } else if (upd->vp[MESH_DISPLACEMENT1] == -1) {
    v            = pd->ShapeVar;
    dofs         = ei->dof[v];
    for (p = 0; p < dim; p++)
      {
	fv->d[p]     = 0.0;
	fv_old->d[p] = 0.0;
	fv_dot->d[p] = 0.0;
	fv_dot_dot->d[p] = 0.0;

	fv->x[p]     = 0.;
	for (i = 0; i < dofs; i++)
	  {
	    node = ei->dof_list[v][i];
	    index = Proc_Elem_Connect[Proc_Connect_Ptr[ei->ielem] +node];
	    fv->x[p] +=  ( Coor[p][index] ) * bf[v]->phi[i];
	  }
	fv_old->x[p] = fv->x[p]; /* Fixed grid stays fixed thru time. */
	fv->x0[p] = fv->x[p];
	fv_dot->x[p] = 0.0;
	fv_dot_dot->x[p] = 0.0;
      }
  }
  
  /*
   * SOLID displacement (vector)...
   * and positions (vector)
   */
  for (p = 0; pdv[SOLID_DISPLACEMENT1] && p < dim; p++)
    {
      v = SOLID_DISPLACEMENT1 + p;
      fv->d_rs[p] = 0.;
      fv_old->d_rs[p] = 0.;
      fv_dot->d_rs[p] = 0.;
      fv_dot_dot->d_rs[p] = 0.;
      fv_dot_old->d_rs[p] = 0.;
      fv_dot_dot_old->d_rs[p] = 0.;
	  
      if ( pdv[v] )
	{
	  dofs     = ei->dof[v];

	  for (i = 0; i < dofs; i++)
	    {
	      node = ei->dof_list[R_SOLID1][i];
	      index = Proc_Elem_Connect[Proc_Connect_Ptr[ei->ielem] +node];
	      fv->d_rs[p] += *esp->d_rs[p][i] * bf[v]->phi[i];
 
	      if (pd->TimeIntegration != STEADY)
		{
		  fv_old->d_rs[p] += *esp_old->d_rs[p][i] * bf[v]->phi[i];
		  fv_dot->d_rs[p] += *esp_dot->d_rs[p][i] * bf[v]->phi[i];
		  if (tran->solid_inertia)
		    {
		      fv_dot_dot->d_rs[p] += *esp_dbl_dot->d_rs[p][i] * bf[v]->phi[i];
		    }
		  fv_dot_old->d_rs[p] += *(xdot_old_static - xdot_static +
					   esp_dot->d_rs[p][i]) * bf[v]->phi[i];
		  if (tran->solid_inertia)
		    {
		      fv_dot_dot_old->d_rs[p] += *(x_dbl_dot_old_static - 
						   x_dbl_dot_static + 
						   esp_dbl_dot->d_rs[p][i])*bf[v]->phi[i];
		    }
		}
	      else
		{
                  fv_old->d_rs[p] = fv->d_rs[p];
		}

	    }
	}

    }

  status = load_coordinate_scales(pd->CoordinateSystem, fv);
  EH(status, "load_coordinate_scales(fv)");

  /*
   * Velocity (vector)...
   */
  
  velodim = dim;		/* Later, this might include v_theta... */

  /* Change here in case the 3rd velocity component is non-zero
   * (non LSA situation). */
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN)
    velodim = dim + 1; /* Later is Now!  Woo!!! */

  /*
   * Default: all velocities are zero...
   */
  for ( p=0; p<velodim; p++)
    {
      v = VELOCITY1 + p;
      if ( pdv[v] )
	{
	  dofs     = ei->dof[v];
          fv->v[p]     = 0.;
          fv_old->v[p] = 0.;
          fv_dot->v[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->v[p] += *esp->v[p][i] * bf[v]->phi[i];
	      if ( pd->TimeIntegration != STEADY )
		{
		  fv_old->v[p] += *esp_old->v[p][i] * bf[v]->phi[i];
		  fv_dot->v[p] += *esp_dot->v[p][i] * bf[v]->phi[i];
		}
	    }
	}
      stateVector[VELOCITY1+p] = fv->v[p];
    }

  /* 
   * Particle velocity (vector)...
   */
  
  /*
   * Default: all velocities are zero...
   */
  /*  for ( p=velodim; p<DIM; p++)
      {
      v = PVELOCITY1 + p;
      if ( pd->v[v] || (upd->vp[v] == -1) )
      {
      fv->pv[p]     = 0.;
      fv_old->pv[p] = 0.;
      fv_dot->pv[p] = 0.;
      }
      } */

  for ( p=0; pdv[PVELOCITY1] && p<velodim; p++)
    {
      v = PVELOCITY1 + p;
      if ( pdv[v] )
	{
	  fv->pv[p] = 0.0; fv_old->pv[p]=0.0;fv_dot->pv[p]=0.0;
		  
	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv->pv[p] += *esp->pv[p][i] * bf[v]->phi[i];

	      if ( pd->TimeIntegration != STEADY )
		{
		  fv_old->pv[p] += *esp_old->pv[p][i] * bf[v]->phi[i];
		  fv_dot->pv[p] += *esp_dot->pv[p][i] * bf[v]->phi[i];
		}
	    }
	}
    }

  /* Extension velocity */

  v = EXT_VELOCITY;
  if ( pdv[EXT_VELOCITY] )
    {
      scalar_fv_fill(esp->ext_v, esp_dot->ext_v, esp_old->ext_v, bf[v]->phi, 
		     ei->dof[v], &(fv->ext_v), &(fv_dot->ext_v), &(fv_old->ext_v));
      stateVector[EXT_VELOCITY] = fv->ext_v;
    }

  /* Electric Field */

  /*
   * Default: all velocities are zero...
   */
  for ( p=0; pdv[EFIELD1] && p<velodim; p++)
    {
      v =  EFIELD1 + p;
      
      if ( pdv[v] )
	{
	  fv->E_field[p]     = 0.;

	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv->E_field[p] += *esp->E_field[p][i] * bf[v]->phi[i];
	    }
	}
    }

  /* Phase functions */

  for( p=0; pdv[PHASE1] && p<pfd->num_phase_funcs; p++)
    {
      v = PHASE1 + p;
      if ( pdv[v] )
	{
	  scalar_fv_fill(esp->pF[p], esp_dot->pF[p], esp_old->pF[p], bf[v]->phi, ei->dof[v],
			 &(fv->pF[p]), &(fv_dot->pF[p]), &(fv_old->pF[p]));
	}
    }

  /*
   * Polymer Stress (tensor)...
   */
  
  for ( mode=0; pdv[POLYMER_STRESS11] && mode<vn->modes; mode++)
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
              v = v_s[mode][p][q];
              if ( pdv[v] || (upd->vp[v] == -1) )
		{
		  /* good default behavior */
		  fv->S[mode][p][q] = 0.;
		  fv_old->S[mode][p][q] = 0.;
		  fv_dot->S[mode][p][q] = 0.;
                }
            }
        }
    }
  
  for ( mode=0; pdv[POLYMER_STRESS11] &&  mode<vn->modes; mode++)
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      if( p <= q)
		{
		  v = v_s[mode][p][q];
		  if ( pdv[v] )
		    {
		      dofs     = ei->dof[v];
		      for ( i=0; i<dofs; i++)
			{
			  fv->S[mode][p][q] += *esp->S[mode][p][q][i] * bf[v]->phi[i];
			  if ( pd->TimeIntegration != STEADY )
			    {
			      fv_old->S[mode][p][q] += *esp_old->S[mode][p][q][i] * bf[v]->phi[i];
			      fv_dot->S[mode][p][q] += *esp_dot->S[mode][p][q][i] * bf[v]->phi[i];
			    }
			}
		    }
		  /* form the entire symmetric stress matrix for the momentum equation */
		  fv->S[mode][q][p] = fv->S[mode][p][q];
		  fv_old->S[mode][q][p] = fv_old->S[mode][p][q];
		  fv_dot->S[mode][q][p] = fv_dot->S[mode][p][q];
		}
	    }
	}
    }


  /*
   * Velocity Gradient (tensor)...
   */
   
  for ( p=0; pdv[VELOCITY_GRADIENT11] && p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  v = v_g[p][q];
	  if ( pdv[v] )
	    {
	      fv->G[p][q] = fv_old->G[p][q] = fv_dot->G[p][q] = 0.0;
	      dofs     = ei->dof[v];
	      for ( i=0; i<dofs; i++)
		{
		  fv->G[p][q] += *esp->G[p][q][i] * bf[v]->phi[i];
		  if ( pd->TimeIntegration != STEADY )
		    {
		      fv_old->G[p][q] += *esp_old->G[p][q][i] * bf[v]->phi[i];
		      fv_dot->G[p][q] += *esp_dot->G[p][q][i] * bf[v]->phi[i];
		    }
		}
	    }
	}
    }
  
  /*
   * Species Unknown Variable
   *      Modifications for non-dilute systems:
   *      -> We need to calculate the mole fraction of the last species
   *         in the mechanism, even if there isn't an equation
   *         for it. This is done via the sum MF = 1
   *         constraint, here. Generalization to other equations
   *         of state should be added here in the future.
   */


  if (pdv[MASS_FRACTION])
    {
      v = MASS_FRACTION;
      if (pd->Num_Species_Eqn != pd->Num_Species) {
	N = pd->Num_Species - 1;

	fv->c[N] = fv_old->c[N] = fv_dot->c[N] = 0.;
	for (w = 0; w < pd->Num_Species_Eqn; w++) {
	  scalar_fv_fill(esp->c[w], esp_dot->c[w], esp_old->c[w], 
			 bf[v]->phi, ei->dof[v],
			 &(fv->c[w]), &(fv_dot->c[w]), &(fv_old->c[w]));
	  fv->c[N]     -= fv->c[w];
	  fv_old->c[N] -= fv_old->c[w];
	  fv_dot->c[N] -= fv_dot->c[w];
	  stateVector[SPECIES_UNK_0+w] = fv->c[w];
	}
	switch (mp->Species_Var_Type) {
	case SPECIES_UNDEFINED_FORM:
	case SPECIES_MOLE_FRACTION:
	case SPECIES_MASS_FRACTION:
	case SPECIES_VOL_FRACTION:
	  fv->c[N] += 1.0;
	  stateVector[SPECIES_UNK_0+N] = fv->c[N];
	  break;
	case SPECIES_CONCENTRATION:
	  /*
	   *  HKM:  These are currently problematic cases that
	   *  argue for including the last species as a dependent
	   *  variable in the solution vector. The reason is that
	   *  a small nonlinear iteration must be performed here
	   *  for calculation of the total species concentration
	   *  and the species concentration of the last species
	   *  in the mechanism. An equation of state has to
	   *  be assumed, and then inverted. 
	   *
	   *   C = Sum_i=1,N[C_i],
	   *
	   *  where C is determined from the equation of state
	   *  as a function of T, P, [other_state_variables], and
	   *  C_1 to C_N.
	   *
	   *  If the problem is
	   *  dilute, then the whole process is circumvented:
	   *     C = C_N.
	   *  or at least C isn't a function of the last species
	   *  in the mechanism:
	   *     C = C(T, P, C_1, ..., C_N-1)
	   *
	   *  The dilute case is carried out below. However, 
	   *  note that c[N] = 0 is assumed in the function
	   *  call. If c[N] is used, then this routine will
	   *  give the wrong answer, until a nonlinear iteration
	   *  loop is installed below.
	   *
	   * -> The same is true for the SPECIES_DENSITY
	   *    loop.
	   */
	  switch (mp->DensityModel) {
	  case DENSITY_CONSTANT_LAST_CONC:
	    fv->c[N] = mp->u_density[0];
	    break;
	  default:
	    rho = calc_concentration(mp, FALSE, NULL);
	    fv->c[N] += rho;
	  }
	  stateVector[SPECIES_UNK_0+N] = fv->c[N];
	  break;
	case SPECIES_DENSITY:
	  /* note, this won't work for time dependent density models, but I got tired ... RRR*/
	  rho = calc_density(mp, FALSE, NULL, 0.);
	  fv->c[N] += rho;
	  stateVector[SPECIES_UNK_0+N] = fv->c[N];
	  break;
	case SPECIES_CAP_PRESSURE:
	  break;
	default:
	  break;
	}
	mp->StateVector_speciesVT = upd->Species_Var_Type;
      } else {
	for (w = 0; w < pd->Num_Species_Eqn; w++) {
	  scalar_fv_fill(esp->c[w], esp_dot->c[w], esp_old->c[w], 
			 bf[v]->phi, ei->dof[v],
			 &(fv->c[w]), &(fv_dot->c[w]), &(fv_old->c[w]));
	  stateVector[SPECIES_UNK_0+w] = fv->c[w];
	}
      }
    }     /* HKM->
	   *  (we may want to insert reference concentrations here
	   *   if available)
	   */

  /*else if (upd->vp[v] == -1) {

  for (w = 0; w < pd->Num_Species; w++) {
  fv->c[w]     = 0.;
  fv_old->c[w] = 0.;
  fv_dot->c[w] = 0.;
  stateVector[SPECIES_UNK_0+w] = fv->c[w];
  }
  } */
	
  /*
   * Porous media Variables
   */


  /*  if (pd->v[v] || (upd->vp[v] == -1) )
      {
      fv->p_liq=0.;
      fv_old->p_liq = 0.;
      fv_dot->p_liq = 0.;
      fv_dot_old->p_liq = 0.0;
      } */

  if (pdv[POR_LIQ_PRES]) 
    {
      v = POR_LIQ_PRES;
      fv->p_liq=0.;
      fv_old->p_liq = 0.;
      fv_dot->p_liq = 0.;
      fv_dot_old->p_liq = 0.0;
      dofs     = ei->dof[v];
      for (i = 0; i < dofs; i++) {
	fv->p_liq += *esp->p_liq[i] * bf[v]->phi[i];
	if (pd->TimeIntegration != STEADY) {
	  fv_old->p_liq += *esp_old->p_liq[i] * bf[v]->phi[i];
	  fv_dot->p_liq += *esp_dot->p_liq[i] * bf[v]->phi[i];
	  fv_dot_old->p_liq += 
	    *(xdot_old_static - xdot_static + esp_dot->p_liq[i]) *
	    bf[v]->phi[i];
	}
      }
    }


  /*  if (pd->v[v] || (upd->vp[v] == -1) )
      {
      fv->p_gas=0.;
      fv_old->p_gas = 0.;
      fv_dot->p_gas = 0.;
      fv_dot_old->p_gas = 0.0;
      }
  */
  if (pdv[POR_GAS_PRES]) 
    {
      v = POR_GAS_PRES;
      fv->p_gas=0.;
      fv_old->p_gas = 0.;
      fv_dot->p_gas = 0.;
      fv_dot_old->p_gas = 0.0;
      dofs     = ei->dof[v];
      for (i = 0; i < dofs; i++) {
	fv->p_gas += *esp->p_gas[i] * bf[v]->phi[i];
	if (pd->TimeIntegration != STEADY) {
	  fv_old->p_gas += *esp_old->p_gas[i] * bf[v]->phi[i];
	  fv_dot->p_gas += *esp_dot->p_gas[i] * bf[v]->phi[i];
	  fv_dot_old->p_gas += *(xdot_old_static - xdot_static + 
				 esp_dot->p_gas[i]) * bf[v]->phi[i];
	}
      }
    }


  /*  if (pd->v[v] || (upd->vp[v] == -1) )
      {
      fv->porosity=0.;
      fv_old->porosity = 0.;
      fv_dot->porosity = 0.;
      fv_dot_old->porosity = 0.0;
      }
  */
  if (pdv[POR_POROSITY]) 
    {
      v = POR_POROSITY;
      fv->porosity=0.;
      fv_old->porosity = 0.;
      fv_dot->porosity = 0.;
      fv_dot_old->porosity = 0.0;
	  
      dofs     = ei->dof[v];
      for (i = 0; i < dofs; i++) {
	fv->porosity += *esp->porosity[i] * bf[v]->phi[i];
	if (pd->TimeIntegration != STEADY) {
	  fv_old->porosity += *esp_old->porosity[i] * bf[v]->phi[i];
	  fv_dot->porosity += *esp_dot->porosity[i] * bf[v]->phi[i];
	  fv_dot_old->porosity += *(xdot_old_static - xdot_static + 
				    esp_dot->porosity[i]) * bf[v]->phi[i];
	}
      }
    }


  if (pdv[POR_TEMP]) 
    {
      v = POR_TEMP;
      fv->T=0.;
      fv_old->T = 0.;
      fv_dot->T = 0.;
      fv_dot_old->T = 0.0;
	  
      dofs     = ei->dof[v];
      for (i = 0; i < dofs; i++) {
	fv->T += *esp->T[i] * bf[v]->phi[i];
	if (pd->TimeIntegration != STEADY) {
	  fv_old->T += *esp_old->T[i] * bf[v]->phi[i];
	  fv_dot->T += *esp_dot->T[i] * bf[v]->phi[i];
	  fv_dot_old->T += *(xdot_old_static - xdot_static + 
			     esp_dot->T[i]) * bf[v]->phi[i];
			  
	}
      }
    }

  /* 
   * Vorticity principle flow direction
   */
  for(p = 0; pdv[VORT_DIR1] && p < DIM; p++)
    {
      v = VORT_DIR1 + p;
      /*if (pd->v[v] || upd->vp[v] == -1) fv->vd[p] = 0.0; */

      if(pdv[v])
	{
	  dofs = ei->dof[v];
	  fv->vd[p] = 0.0;
	  for(i = 0; i < dofs; i++)
	    fv->vd[p] += *esp->vd[p][i] * bf[v]->phi[i];
	}
    }

  /* 
   * Lagrange Multiplier Field
   */
  for(p = 0; pdv[LAGR_MULT1] && p < DIM; p++)
    {
      v = LAGR_MULT1 + p;
      /*      if (pd->v[v] || upd->vp[v] == -1)
	      {
	      fv->lm[p] = 0.0;
	      fv_old->lm[p] = 0.0;
	      }
      */	  
      if(pdv[v])
	{
	  fv->lm[p] = 0.0;
          fv_old->lm[p] = 0.0;

	  dofs = ei->dof[v];
	  for(i = 0; i < dofs; i++)
	    {
	      fv->lm[p] += *esp->lm[p][i] * bf[v]->phi[i];
	      if (pd->TimeIntegration != STEADY) {
		fv_old->lm[p] += *esp_old->lm[p][i] * bf[v]->phi[i];
	      }
	    }
	}
    }

  /*
   * Eigenvalue associated with vd.
   */
  
  /*  if (pd->v[v] || upd->vp[v] == -1) fv->vlambda = 0.0; */

  if (pdv[VORT_LAMBDA])
    {
      v = VORT_LAMBDA;
      fv->vlambda = 0.0;
      dofs = ei->dof[v];
      for(i = 0; i < dofs; i++)
	fv->vlambda += *esp->vlambda[i] * bf[v]->phi[i];
    }

  /*
   * External...
   */
  if (efv->ev)
    {
      int table_id=0;
      v = EXTERNAL;
      for ( w=0; w<efv->Num_external_field; w++)
	{
	  dofs     = ei->dof_ext[w];
	  fv->external_field[w] = 0.;
	  fv_old->external_field[w] = 0.;
	  fv_dot->external_field[w] = 0.;
	  
 	  if( efv->i[w] != I_TABLE )
	    {
	      for ( i=0; i<dofs; i++)
		{
		  fv->external_field[w] += *evp->external_field[w][i] * bfex[w]->phi[i];

		  if ( pd->TimeIntegration != STEADY )
		    {
		      fv_old->external_field[w] += *evp->external_field[w][i] * bfex[w]->phi[i];
		      fv_dot->external_field[w] += 0.;
		    }
		}
	    }
	  else
	    {
	      double slope;
	      fv->external_field[w] = interpolate_table( ext_Tables[table_id],
 							 fv->x, &slope, NULL);
	      table_id++;
	    }

          /* 
           * If the variable name is velocity, and the momentum equations are not active,
           * load the external_fields into the velocity fv for use in Advection-diffusion analysis
           */
          if (strcmp(efv->name[w], "VX") == 0 && !pd->v[VELOCITY1])
	    {
	      fv->v[0] = fv->external_field[w];
	      fv_old->v[0] = fv_old->external_field[w];
	      fv_dot->v[0] = fv_dot->external_field[w];
	    }
          if (strcmp(efv->name[w], "VY") == 0 && !pd->v[VELOCITY2])
	    {
	      fv->v[1] = fv->external_field[w];
	      fv_old->v[1] = fv_old->external_field[w];
	      fv_dot->v[1] = fv_dot->external_field[w];
	    }
          if (strcmp(efv->name[w], "VZ") == 0 && !pd->v[VELOCITY3])
	    {
	      fv->v[2] = fv->external_field[w];
	      fv_old->v[2] = fv_old->external_field[w];
	      fv_dot->v[2] = fv_dot->external_field[w];
	    }
          /* 
           * If the variable name is mesh displacement, and the mesh 
	   * equations are not active, load the external_fields into 
	   * the mesh fv.  PRS NOTE: I don't think we will need these 
	   * for the decoupled JAS/GOMA, as the displacments are 
	   * swallowed into the mesh
           */
          if (strcmp(efv->name[w], "DMX") == 0 && !pd->v[MESH_DISPLACEMENT1])
	    {
	      fv->d[0] = fv->external_field[w];
	      fv_old->d[0] = fv_old->external_field[w];
	      fv_dot->d[0] = fv_dot->external_field[w];

	    }
          if (strcmp(efv->name[w], "DMY") == 0 && !pd->v[MESH_DISPLACEMENT2])
	    {
	      fv->d[1] = fv->external_field[w];
	      fv_old->d[1] = fv_old->external_field[w];
	      fv_dot->d[1] = fv_dot->external_field[w];
	    }
          if (strcmp(efv->name[w], "DMZ") == 0 && !pd->v[MESH_DISPLACEMENT3])
	    {
	      fv->d[2] = fv->external_field[w];
	      fv_old->d[2] = fv_old->external_field[w];
	      fv_dot->d[2] = fv_dot->external_field[w];
	    }
          /* 
           * If the variable name is porosity, and the porosity equation is not active,
           * load the external_fields into the porosity fv
           */
          if (strcmp(efv->name[w], "P_POR") == 0 && !pd->v[POR_POROSITY])
	    {
	      fv->porosity = fv->external_field[w];
	      fv_old->porosity = fv_old->external_field[w];
	      fv_dot->porosity = fv_dot->external_field[w];
	    }
	  if (strcmp(efv->name[w], "F") == 0 && !pd->v[FILL])
	    {
	      fv->F = fv->external_field[w];
	      fv_old->F = fv_old->external_field[w];
	      fv_dot->F = fv_dot->external_field[w];
	    }
	  if (strcmp(efv->name[w], "F1") == 0 && !pd->v[PHASE1])
	    {
	      fv->pF[0] = fv->external_field[w];
	      fv_old->pF[0] = fv_old->external_field[w];
	      fv_dot->pF[0] = fv_dot->external_field[w];
	    }

	}
    }
  else
    {
      for ( w=0; w<efv->Num_external_field; w++)
	{
	  fv->external_field[w]     = 0.;
	  fv_old->external_field[w] = 0.;
	}
    }

  /* initial displacements for TALE anneals. all for KINEMATIC DISPLACEMENT BC */
  if (efv->TALE)
    {
      for ( w=0; w<dim; w++)
	{
	  v = MESH_DISPLACEMENT1 + w;
	  dofs     = ei->dof[v];
	  fv->initial_displacements[w] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->initial_displacements[w] += *evp->initial_displacements[w][i] * bf[v]->phi[i];

	    }
	  v = SOLID_DISPLACEMENT1 + w;
	  dofs     = ei->dof[v];
	  fv->initial_displacements[w+DIM] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->initial_displacements[w+DIM] += *evp->initial_displacements[w+DIM][i] * bf[v]->phi[i];

	    }
	}
    }
  else
    {
      for ( w=0; w<dim; w++)
	{
	  fv->initial_displacements[w]     = 0.;
	  fv->initial_displacements[w + DIM]     = 0.;
	}
    }
	      
  return(status);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
int
load_fv_grads(void)
     /*******************************************************************************
      *
      * load_fv_grads() -- load relevant field variable gradients at this gauss pt
      *
      * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
      * ------
      *
      *
      * output: ( this routine fills in the following parts of the fv structure)
      * -------
      *		fv->
      *			grad_T -- temperature gradient
      *			grad_P -- pressure gradient
      *			grad_c -- species unknown gradient
      *			grad_F -- fill gradient
      *			grad_V -- voltage potential gradient
      *			div_v  -- divergence of velocity 
      *			grad_v -- velocity gradient tensor
      *                      curl_v -- curl of velocity, a.k.a. vorticity
      *			div_d  -- divergence of displacement (dilatation)
      *			grad_d -- gradient of displacement ("strain")
      *                      grad_X -- gradient of the Eulerian solid reference state.
      *			grad_S -- gradient of the polymer stress tensor
      *			grad_G -- gradient of the velocity gradient tensor
      *                      grad_pv -- gradient of particle velocity
      *                      grad_p_liq, grad_p_gas, grad_porosity
      *                             -- gradient of porous media variables
      *                      grad_n -- gradient of level set normal vector 
      *
      * NOTE: To accommodate shell elements, this function has been modified
      *       so that fv variables are not zeroed out when they are active
      *       on an element block other than the current one.
      *       The check done for variable v is then:
      *          if ( pd->v[v] || upd->vp[v] == -1 )
      *       In many cases below, this conditional zeroing is done in
      *       a separate small loop before the main one.
      *
      * Return values:
      *		0 -- if things went OK
      *	       -1 -- if a problem occurred
      *
      * Created:	Fri Mar 18 07:36:14 MST 1994 pasacki@sandia.gov
      *
      * Modified:	Tue Feb 21 11:08 MST 1995 pasacki@sandia.gov
      ***************************************************************************/
{
  int v;			/* variable type indicator */
  int i, a;			/* index */
  int p, q, r, s;		/* dimension index */
  int dofs;			/* degrees of freedom */
  int w;			/* concentration species */
  int dim = pd->Num_Dim;
  int mode;                     /* modal counter */
  int status = 0;
  int transient_run =  (pd->TimeIntegration != STEADY) ;
  BASIS_FUNCTIONS_STRUCT *bfn;

  int wim = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN) {
    wim = wim+1;
  }

  /* Use a static flag so unused grads are zero on first call, but are not zero subsequently
  *  This is for efficieny
  */  
  static int zero_unused_grads = TRUE;

  /*
   * grad(T)
   */
  if (pd->v[TEMPERATURE])
    {
      v = TEMPERATURE;
      for (p = 0; p < VIM; p++) fv->grad_T[p] = 0.0;
      dofs  = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_T[p] += *esp->T[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } 
  else if (zero_unused_grads && upd->vp[TEMPERATURE] == -1 )
    {
      for (p=0; p<VIM; p++) fv->grad_T[p] = 0.0;
    }
   
  /*
   * grad(P)
   */
  
  if ( pd->v[PRESSURE] )
    {
      v = PRESSURE;
      dofs  = ei->dof[v];
#ifdef DO_NO_UNROLL
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_P[p] = 0.0;
		  
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_P[p] += *esp->P[i] * bf[v]->grad_phi[i] [p];
	    }
	}
#else
      grad_scalar_fv_fill( esp->P, bf[v]->grad_phi, dofs, fv->grad_P);
#endif
	  
    } else if ( zero_unused_grads &&  upd->vp[PRESSURE] == -1 ) {
      for (p=0; p<VIM; p++) fv->grad_P[p] = 0.0;
    }
  
  /*
   * grad(nn)
   */
  
  if ( pd->v[BOND_EVOLUTION] )
    {
      v = BOND_EVOLUTION;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_nn[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_nn[p] += *esp->nn[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  
  /*
   * grad(F)
   */
  
  if ( pd->v[FILL] )
  {
	  v = FILL;
      dofs  = ei->dof[v];
#ifdef DO_NO_UNROLL
      for ( p=0; p<VIM; p++)
	  {
		 fv->grad_F[p] = 0.0;fv_old->grad_F[p] = 0.0; fv_dot->grad_F[p] = 0.0;
		  
		  for ( i=0; i<dofs; i++)
		  {
			 fv->grad_F[p] += *esp->F[i] * bf[v]->grad_phi[i] [p];  
			  if ( transient_run )
			  {
				  /* keep this only for VOF/Taylor-Galerkin stuff */
				  fv_old->grad_F[p] += 
				  *esp_old->F[i] * bf[v]->grad_phi[i][p];
				  
				  fv_dot->grad_F[p] += 
					  *esp_dot->F[i] * bf[v]->grad_phi[i][p];
			  }
		  }
	  }
#else
	  grad_scalar_fv_fill( esp->F, bf[v]->grad_phi, dofs, fv->grad_F ); 
	  
	  if( transient_run )
	  {
		  grad_scalar_fv_fill ( esp_old->F,  bf[v]->grad_phi, dofs, fv_old->grad_F);
		  grad_scalar_fv_fill ( esp_dot->F,  bf[v]->grad_phi, dofs, fv_dot->grad_F);		  
	  }
#endif
	  
  } else if ( zero_unused_grads &&  upd->vp[FILL] == -1 ) {
      for (p=0; p<VIM; p++)
	  {
		  fv->grad_F[p] = 0.0; fv_old->grad_F[p] = 0.0; fv_dot->grad_F[p] = 0.0;
	  }
  }
	  
  
  
  /*
   * grad(H)
   */
  if ( pd->v[CURVATURE] )
    {
      v = CURVATURE;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_H[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_H[p] += *esp->H[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  else if ( zero_unused_grads &&  upd->vp[CURVATURE] == -1 )
    {
      for (p=0; p<VIM; p++) 
	{
	  fv->grad_H[p] = 0.0;
	}
    }
  
  /*
   * grad(V)
   */ 
  if ( pd->v[VOLTAGE] )
  {
	  v = VOLTAGE;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	  {
		  fv->grad_V[p] = 0.0;
		  for ( i=0; i<dofs; i++)
		  {
			  fv->grad_V[p] += *esp->V[i] * bf[v]->grad_phi[i] [p];
		  }
	  }
  }
  else  if ( zero_unused_grads &&  upd->vp[VOLTAGE] == -1 ) {
      for (p=0; p<VIM; p++) fv->grad_V[p] = 0.0;
  } 
  
  /*
   * grad(Enorm) 
   */
  if (pd->v[ENORM])
    {
      v = ENORM;
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_Enorm[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_Enorm[p] += *esp->Enorm[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  
  /*
   * grad(qs)
   */
  if (pd->v[SURF_CHARGE])
    {
      v = SURF_CHARGE;
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_qs[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      // HKM -> I changed dphidxi to grad_phi below. There didn't seem
	      //        to be a case where the raw element derivative (dphidxi), which
	      //        doesn't even take into account of the size of the element, should
	      //        be used at this level. 
	      //   fv->grad_qs[p] += *esp->qs[i] * bf[v]->dphidxi[i][p];
	      fv->grad_qs[p] += *esp->qs[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  else if (zero_unused_grads &&  upd->vp[SURF_CHARGE] == -1)
    {
      for (p = 0; p < VIM-1; p++)
	{
	  fv->grad_qs[p] = 0.0;
	}
    }
  
  
  /*
   * grad(SH)
   */
  
  if ( pd->v[SHEAR_RATE] )
  {
	  v = SHEAR_RATE;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	  {
		  fv->grad_SH[p] = 0.0;
		  for ( i=0; i<dofs; i++)
		  {
			  fv->grad_SH[p] += *esp->SH[i] * bf[v]->grad_phi[i] [p];
		  }
	  }
  } else  if (zero_unused_grads &&  upd->vp[SHEAR_RATE] == -1 )   {
      for (p=0; p<VIM; p++) fv->grad_SH[p] = 0.0;
  }
  
  
  /*
   * grad(d)
   */

  if (pd->v[MESH_DISPLACEMENT1] && !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1])
    {
      v = MESH_DISPLACEMENT1;
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      fv->grad_d[p][q]     = 0.0;
	      fv_old->grad_d[p][q] = 0.0;
	      for (r = 0; r < VIM; r++)
		{
		  for (i = 0; i < dofs; i++)
		    {
		      fv->grad_d[p][q]     += *esp->d[r][i]     * bf[v]->grad_phi_e[i][r] [p][q];
		      fv_old->grad_d[p][q] += *esp_old->d[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    } else if (upd->vp[MESH_DISPLACEMENT1] == -1) 
      {
	for (p = 0; p < VIM; p++)
	  {
	    for (q = 0; q < VIM; q++)
	      {
		fv->grad_d[p][q] = 0.0;
		fv_old->grad_d[p][q] = 0.0;
		fv->grad_d_dot[p][q] = 0.0;
	      }
	  }
      }
  
  /*
   * grad(d_dot)
   */
  if (pd->v[MESH_DISPLACEMENT1] && !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1]) 
    {
      for (p = 0; p < VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      v = MESH_DISPLACEMENT1 + p;
	      if ( pd->v[v])
		{
		  dofs     = ei->dof[v];
		  fv->grad_d_dot[p][q] = 0.;
		  for (r = 0; r < dim; r++)
		    {
		      for ( i=0; i<dofs; i++)
			{
			  fv->grad_d_dot[p][q] += 
			    *esp_dot->d[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
			}
		    }
		}
	    }
	}
    }
  
  /*
   * grad(d_rs)
   */
  if( pd->v[SOLID_DISPLACEMENT1] )
    {
      for ( p=0;  p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      v = SOLID_DISPLACEMENT1 + p;
			  
	      if ( pd->v[v] )
		{
		  dofs     = ei->dof[v];
		  fv->grad_d_rs[p][q] = 0.0;
		  for ( r = 0; r < dim; r++)
		    {
		      for ( i=0; i<dofs; i++)
			{
			  fv->grad_d_rs[p][q] += 
			    *esp->d_rs[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
						  
			}
		    }
		} 
	    }
	}
    } else if ( zero_unused_grads &&  upd->vp[SOLID_DISPLACEMENT1] == -1 )
      {
	memset(fv->grad_d_rs, 0, sizeof(double) *VIM*VIM);
      }
	  
  /*
   * div(d)
   */

  
  /*
   * div(d_dot)
   */

  
  if (pd->v[MESH_DISPLACEMENT1] && !ShellElementParentElementCoverageForVariable[MESH_DISPLACEMENT1])
    {
      v = MESH_DISPLACEMENT1;
      fv->div_d_dot = 0.0;
      fv->div_d = 0.0;
	  
      for (p = 0; p < VIM; p++)
	{
          fv->div_d_dot += fv->grad_d_dot[p][p];
	  fv->div_d += fv->grad_d[p][p];
	}
    } else if (zero_unused_grads && upd->vp[MESH_DISPLACEMENT1] == -1) {
      fv->div_d_dot = 0.0;
      fv->div_d     = 0.0;
    }
  
  
  /*
   * div(d_rs)
   */
  
  if ( pd->v[SOLID_DISPLACEMENT1] )
  {
      for ( p=0; p<VIM; p++)
	  {
          fv->div_d_rs += fv->grad_d_rs[p][p];
	  }
  } else if ( zero_unused_grads &&  upd->vp[SOLID_DISPLACEMENT1] == -1 ) {
	  fv->div_d_rs = 0.;
  }
  
  /*
   * grad(v)
   */
  if (pd->v[VELOCITY1])
    {
      dofs  = ei->dof[VELOCITY1];
      v = VELOCITY1;  
#ifdef DO_NO_UNROLL	  
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      fv->grad_v[p][q] = 0.0;
	      for (r = 0; r < wim; r++)
		{
		  for (i = 0; i < dofs; i++)
		    {
		      fv->grad_v[p][q] += (*esp->v[r][i]) * bf[v]->grad_phi_e[i][r] [p][q];	      
		    }
		}
	    }
	}
#else
      grad_vector_fv_fill(esp->v, bf[v]->grad_phi_e, dofs, fv->grad_v);
#endif
    } 
  else if (zero_unused_grads && upd->vp[VELOCITY1] == -1) 
    {
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      fv->grad_v[p][q] = 0.0; 
	    }
	}
    }


  /*
   * grad(pv), particle velocity gradients.
   */
  if (pd->v[PVELOCITY1])
    {
      dofs = ei->dof[PVELOCITY1];
      v = PVELOCITY1;

      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      for (r = 0; r < wim; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv->grad_pv[p][q] += 
			*esp->pv[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    } else if ( zero_unused_grads && upd->vp[PVELOCITY1] == -1 ) {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      fv->grad_pv[p][q] = 0.0;
	    }
	}

    }
  
  /*
   * grad(ext_v), extension velocity gradients.
   */


  if( pd->v[EXT_VELOCITY] )
  {
	  dofs = ei->dof[EXT_VELOCITY];
	  v = EXT_VELOCITY;
	  
	  for ( p=0; p<VIM; p++)
	  {
		  fv->grad_ext_v[p] = 0.0;
		  for ( i=0; i<dofs; i++)
		  {
			  fv->grad_ext_v[p] += 
			  *esp->ext_v[i] * bf[v]->grad_phi[i][p];
		  }
	  }
  } else if (zero_unused_grads &&  upd->vp[EXT_VELOCITY] == -1 ) {
      for ( p=0; p<VIM; p++)
	  {
		  fv->grad_ext_v[p] = 0.0;
	  }
  }	  
  
  /* 
   *  grad( n ) .  Level set or shell normal vector gradient 
   */
  if (pd->v[NORMAL1] || pd->v[SHELL_NORMAL1])
    {
      if (pd->v[NORMAL1]) v = NORMAL1;
      if (pd->v[SHELL_NORMAL1]) v = SHELL_NORMAL1;
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
            {
	      fv->grad_n[p][q] = 0.0;
	      for (r = 0; r < dim; r++)
		{
		  for (i = 0; i < dofs; i++)
		    {
		      fv->grad_n[p][q] +=
			*(esp->n[r][i]) * bf[v]->grad_phi_e[i][r][p][q];
		    }
		}
	    }
	}
    }
  else if ( (zero_unused_grads && upd->vp[NORMAL1] == -1) || upd->vp[SHELL_NORMAL1] == -1)
    {
      for(p = 0; p < VIM; p++ )
	{
	  for (q = 0; q < VIM; q++ )
	    {
	      fv->grad_n[p][q] = 0.0;
	    }
	}
    }
  
  /*
   *  div(n)
   */
  if (pd->v[NORMAL1] || pd->v[SHELL_NORMAL1])
    {
      fv->div_n = 0.0;
      for (p = 0 ; p < VIM; p++)
	{
	  fv->div_n += fv->grad_n[p][p];
	}
    }
  
  /*
   * div_s(n)
   */
  if (pd->v[NORMAL1] || pd->v[SHELL_NORMAL1])
    {
      fv->div_s_n = 0.0;
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      fv->div_s_n += (delta(p,q) * fv->grad_n[p][q]
			      - fv->n[p] * fv->n[q] * fv->grad_n[q][p]);
	    }
	}
#ifdef DEBUG_HKM
      /*
	for (p = 0; p < VIM; p++)
	{
	for (q = 0; q < VIM; q++)
	{
	printf("load_fv_grads: fv->grad_n[%d][%d] = %g\n", p, q, fv->grad_n[p][q]);
	}
	}
	printf("load_fv_grads: fv->div_s_n = %g\n", fv->div_s_n);
      */
#endif 
    }

  // Calculation of the surface curvature dyadic
  //             = - del_s (n)  = - (I - n n ) grad n
  if (pd->v[SHELL_NORMAL1])
    {
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
            {
	      fv->surfCurvatureDyadic[p][q] = - fv->grad_n[p][q];
	      for (a = 0; a < VIM; a++) 
		{
		  fv->surfCurvatureDyadic[p][q] += fv->n[p] * fv->n[a] * fv->grad_n[a][q];
		}
	    }
	}
    }
  
  /*
   * grad(E_field), Electric field gradients.
   */  
  if (pd->v[EFIELD1])
    {
      dofs = ei->dof[EFIELD1];
      v = EFIELD1;
	  
      for (p = 0; p <VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      fv->grad_E_field[p][q] = 0.0;
	      for (r = 0; r < wim; r++)
		{
		  for (i = 0; i < dofs; i++)
		    {
		      fv->grad_E_field[p][q] += 
			*esp->E_field[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    } else if (zero_unused_grads && upd->vp[EFIELD1] != -1 )  {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      fv->grad_E_field[p][q] = 0.0;
	    }
	}
    }
  
  
  /* Phase Functions */
  
  if (pd->v[PHASE1])
    {
      for (r = 0; r < pfd->num_phase_funcs; r++ )
	{
	  v = PHASE1 + r;
	  dofs = ei->dof[v];
		  
	  for (p = 0; p < VIM; p++)
	    {
	      fv->grad_pF[r][p] = 0.0;
			  
	      for ( i=0; i<dofs; i++)
		{
		  fv->grad_pF[r][p] += *esp->pF[r][i] * bf[v]->grad_phi[i] [p];
		}
	    }
	}
    }
  
  
  
  /*
   * curl(v)
   */
  if (CURL_V != -1 && !InShellElementWithParentElementCoverage)
    {
      v = VELOCITY1;
      dofs = ei->dof[VELOCITY1];
      bfn = bf[v];
     
      if (pd->v[v] || upd->vp[v] != -1)
	{
	  for (p = 0; p < DIM; p++) fv->curl_v[p] = 0.0;
	}
	  
      if (pd->v[v])
	{
	  /* Always compute all three components of the vorticity
	   * vector.  If we are really in 2D, then the output routine will
	   * only output the 3rd component. DIM is equal to 3.
	   */
	  for (p = 0; p < DIM; p++)
	    {
	      for (i = 0; i < dofs; i++)
		{
		  for (a = 0; a < VIM; a++) /* VIM */
		    {
		      fv->curl_v[p] += *esp->v[a][i] * bfn->curl_phi_e[i][a][p];
		    }
		}
	    }
	}
    }
  
  
  /*
   * div(v)
   */
  
  if (pd->v[VELOCITY1])
    {
      fv->div_v = fv->grad_v[0][0] + fv->grad_v[1][1];
      if (VIM == 3)
	{
	  fv->div_v += fv->grad_v[2][2];
	}
    } 
  else if (zero_unused_grads && upd->vp[VELOCITY1] == -1)
    {
      fv->div_v = 0.0;
    }
  
  
  /*
   * div(pv)
   */
  
  if (pd->v[PVELOCITY1])
    {
      v = PVELOCITY1;
      fv->div_pv = 0.0;
      for ( p=0; p<VIM; p++)
	{
	  fv->div_pv += fv->grad_pv[p][p];
	}
    } else if (zero_unused_grads && upd->vp[PVELOCITY1] == -1) 
      {
	fv->div_pv = 0.0;
      }   
  
  
  /*
   * Species unknown ...
   */
  
  if( pd->v[MASS_FRACTION] )
  {
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	  {
		  v = MASS_FRACTION;
		  dofs     = ei->dof[v];
		  for ( p=0; p<VIM; p++ )
		  {
			  fv->grad_c[w][p] = 0.;
			  fv_old->grad_c[w][p] = 0.;
			  
			  for ( i=0; i<dofs; i++)
			  {
				  fv->grad_c[w][p] += 
				  *esp->c[w][i] * bf[v]->grad_phi[i][p];
				  if ( pd->TimeIntegration != STEADY )
				  {
					  /* keep this only for VOF/Taylor-Galerkin stuff */
					  fv_old->grad_c[w][p] += 
					  *esp_old->c[w][i] * bf[v]->grad_phi[i][p];
				  }
			  }
		  }
	  }
  } else if (zero_unused_grads && upd->vp[MASS_FRACTION] == -1) {
	  	  for ( w=0; w<pd->Num_Species_Eqn; w++)
		  {			  
			  for ( p=0; p<VIM; p++ )
			  {
				  fv->grad_c[w][p] = 0.;
				  fv_old->grad_c[w][p] = 0.;
			  }
		  }
  }
  
  /*
   * Porous media unknowns ...
   */
  


  if ( pd->v[POR_LIQ_PRES] )
  {
	  v = POR_LIQ_PRES;
      dofs     = ei->dof[v];
      for ( p=0; p<VIM; p++ )
	  {
		  fv->grad_p_liq[p] = 0.0;
		  fv_old->grad_p_liq[p] = 0.0;

		  for ( i=0; i<dofs; i++)
		  {
			  fv->grad_p_liq[p] +=
			  *esp->p_liq[i] * bf[v]->grad_phi[i][p];
			  if ( pd->TimeIntegration != STEADY )
			  {
				  /* keep this only for VOF/Taylor-Galerkin stuff */
				  fv_old->grad_p_liq[p] +=
				  *esp_old->p_liq[i] * bf[v]->grad_phi[i][p];
			  }
		  }
	  }
  } else if ( zero_unused_grads && upd->vp[POR_LIQ_PRES] == -1 ) {
      for ( p=0; p<VIM; p++ )
	  {
		  fv->grad_p_liq[p] = 0.0;
		  fv_old->grad_p_liq[p] = 0.0;
	  }
  }
	  
	  


  if ( pd->v[POR_GAS_PRES] )
  {
	    v = POR_GAS_PRES;
      dofs     = ei->dof[v];
      for ( p=0; p<VIM; p++ )
	  {
		  fv->grad_p_gas[p] = 0.0;
		  fv_old->grad_p_gas[p] = 0.0;
		  for ( i=0; i<dofs; i++)
		  {
			  fv->grad_p_gas[p] +=
			  *esp->p_gas[i] * bf[v]->grad_phi[i][p];
			  if ( pd->TimeIntegration != STEADY )
			  {
				  /* keep this only for VOF/Taylor-Galerkin stuff */
				  fv_old->grad_p_gas[p] +=
				  *esp_old->p_gas[i] * bf[v]->grad_phi[i][p];
			  }
		  }
	  }
  } else if (zero_unused_grads &&  upd->vp[POR_GAS_PRES] == -1 ) {
	  for ( p=0; p<VIM; p++ )
	  {
		  fv->grad_p_gas[p] = 0.0;
		  fv_old->grad_p_gas[p] = 0.0;
	  }
  }  


  if (pd->v[POR_POROSITY])
    {
      v = POR_POROSITY;
      dofs = ei->dof[v];  
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_porosity[p] = 0.0;
	  fv_old->grad_porosity[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_porosity[p] += *esp->porosity[i] * bf[v]->grad_phi[i][p];
	      if (pd->TimeIntegration != STEADY)
		{
		  /* keep this only for VOF/Taylor-Galerkin stuff */
		  fv_old->grad_porosity[p] += *esp_old->porosity[i] * bf[v]->grad_phi[i][p];
		}
	    }
	}
    } 
  else if (zero_unused_grads &&  upd->vp[POR_POROSITY] == -1) 
    {
      for (p = 0; p < VIM; p++ )
	{
	  fv->grad_porosity[p] = 0.0;
	  fv_old->grad_porosity[p] = 0.0;
	}
    }
  

  if ( pd->v[POR_TEMP] )
    {
      v = POR_TEMP;
      dofs     = ei->dof[v];
      for ( p=0; p<VIM; p++ )
	{
	  fv->grad_T[p] = 0.;
	  fv_old->grad_T[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_T[p] +=
		*esp->T[i] * bf[v]->grad_phi[i][p];
	      if ( pd->TimeIntegration != STEADY )
		{
		  /* keep this only for VOF/Taylor-Galerkin stuff */
		  fv_old->grad_T[p] +=
		    *esp_old->T[i] * bf[v]->grad_phi[i][p];
		}
	    }
	}
    }
  
  /*
   * grad(S), mode 0
   */
 
  
  if (pd->v[POLYMER_STRESS11])
    {  
      v = POLYMER_STRESS11;
      dofs = ei->dof[v];
      for ( mode=0; mode<vn->modes; mode++)
	{
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  for ( r=0; r<VIM; r++)
		    {
		      fv->grad_S[mode][r][p][q]=0.;
					  
		      for ( i=0; i<dofs; i++)
			{
			  if(p<=q)
			    {
			      fv->grad_S[mode][r][p][q] += 
				*esp->S[mode][p][q][i] * bf[v]->grad_phi[i][r] ;
			    }
			  else
			    {
			      fv->grad_S[mode][r][p][q] += 
				*esp->S[mode][q][p][i] * bf[v]->grad_phi[i][r] ;
			    }
			}
					  
		    }
		}
	    }
	}
	  
      /*
       * div(S) - this is a vector!
       */
      for ( mode=0; mode<vn->modes; mode++)
	  {
		  for ( r=0; r<dim; r++)
		  {
			 fv->div_S[mode][r] = 0.0;
			  
			  for ( q=0; q<dim; q++)
			  {
				  fv->div_S[mode][r] += 
				  fv->grad_S[mode][q][q][r];
			  }
		  }
	  }
      
      if ( pd->CoordinateSystem != CARTESIAN )
	  {
		  for ( mode=0; mode<vn->modes; mode++)
		  {
			  for ( s=0; s<VIM; s++)
			  {
				  for ( r=0; r<VIM; r++)
				  {
					  for ( p=0; p<VIM; p++)
					  {
						  fv->div_S[mode][s] += 
						  fv->S[mode][p][s]*fv->grad_e[p][r][s];
					  }
				  }
			  }
			  for ( s=0; s<VIM; s++)
			  {
				  for ( r=0; r<VIM; r++)
				  {
					  for ( q=0; q<VIM; q++)
					  { 
						  fv->div_S[mode][s] += 
						  fv->S[mode][r][q]* fv->grad_e[q][r][s] ;
					  }
				  }
			  }
		  }
	  }
	  
	  
  } else if ( zero_unused_grads && upd->vp[POLYMER_STRESS11] == -1 )  {  
      for ( mode=0; mode<vn->modes; mode++)
	  {
		  for ( p=0; p<VIM; p++)
		  {
              fv->div_S[mode][p] = 0.0;
			  for ( q=0; q<VIM; q++)
			  {
				  for ( r=0; r<VIM; r++)
				  {
					  fv->grad_S[mode][r][p][q]=0.;
				  }
			  }
		  }
	  }
  }
  
  
  /*
   * grad(G)
   */
  
  if ( pd->v[VELOCITY_GRADIENT11] )
  {  
	  v = VELOCITY_GRADIENT11;
      dofs = ei->dof[v];
      for ( p=0; p<VIM; p++)
	  {
		  for ( q=0; q<VIM; q++)
		  {
			  for ( r=0; r<VIM; r++)
			  {
				  fv->grad_G[r][p][q] = 0.0;
				  
				  for ( i=0; i<dofs; i++)
				  {
					  fv->grad_G[r][p][q] += 
					  *esp->G[p][q][i] * bf[v]->grad_phi[i][r];
				  }
			  }
		  }
	  }
      
      /*
       * div(G) - this is a vector!
       */
      for ( r=0; r<dim; r++)
	  {
		   fv->div_G[r]  = 0.0;
		  for ( q=0; q<dim; q++)
		  {
			  fv->div_G[r] += 
			  fv->grad_G[q][q][r];
		  }
	  }
      
      if ( pd->CoordinateSystem != CARTESIAN )
	  {
          for ( s=0; s<VIM; s++)
		  {
              for ( r=0; r<VIM; r++)
			  {
                  for ( p=0; p<VIM; p++)
				  {
                      fv->div_G[s] +=
					  fv->G[p][s]*fv->grad_e[p][r][s];
				  }
			  }
		  }
		  
          for ( s=0; s<VIM; s++)
		  {
              for ( r=0; r<VIM; r++)
			  {
                  for ( q=0; q<VIM; q++)
				  {
                      fv->div_G[s] +=
					  fv->G[r][q]* fv->grad_e[q][r][s] ;
				  }
			  }
		  }
	  }
	  
      /*
       * div(Gt) - this is a vector and the divergence of the
       * transpose of G.
       */
      for ( r=0; r<dim; r++)
	  {
		   fv->div_Gt[r] = 0.0;
		  for ( q=0; q<dim; q++)
		  {
			  fv->div_Gt[r] += 
			  fv->grad_G[q][r][q];
		  }
	  }
      
      if ( pd->CoordinateSystem != CARTESIAN )
	  {
          for ( s=0; s<VIM; s++)
		  {
              for ( r=0; r<VIM; r++)
			  {
                  for ( p=0; p<VIM; p++)
				  {
                      fv->div_Gt[s] +=
					  fv->G[s][p]*fv->grad_e[p][r][s];
				  }
			  }
		  }
		  
          for ( s=0; s<VIM; s++)
		  {
              for ( r=0; r<VIM; r++)
			  {
                  for ( q=0; q<VIM; q++)
				  {
                      fv->div_Gt[s] +=
					  fv->G[q][r]* fv->grad_e[q][r][s] ;
				  }
			  }
		  }
	  }
	  
  } else if ( zero_unused_grads &&  upd->vp[VELOCITY_GRADIENT11] == -1 ) {  
      for ( p=0; p<VIM; p++)
	  {
          fv->div_G[p] = 0.0;
          fv->div_Gt[p] = 0.0;
		  for ( q=0; q<VIM; q++)
		  {
			  for ( r=0; r<VIM; r++)
			  {
				  fv->grad_G[r][p][q]=0.;
			  }
		  }
	  }
  }
  
  
  /*
   * grad(vd)
   */
  if (pd->v[VORT_DIR1])
    {
      v = VORT_DIR1;;
      dofs = ei->dof[v];
		
      for (p = 0; p < DIM; p++)
	{
          for (q = 0; q < DIM; q++)
	    {
	      fv->grad_vd[p][q] = 0.0;
	      for (r = 0; r < DIM; r++)
		{
		  for (i = 0; i < dofs; i++)
		    {
		      fv->grad_vd[p][q] += 
			*esp->vd[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    } else if (zero_unused_grads &&  upd->vp[VORT_DIR1] == -1 ) {
      for (p=0; p<DIM; p++)
	{
          for (q=0; q<DIM; q++)
	    {
              fv->grad_vd[p][q] = 0.0;
	    }
	}
    }	   
  
  /*
   * div(vd)
   */
   
  if ( pd->v[VORT_DIR1] )
  {
	  fv->div_vd = 0.0;
      for ( p=0; p<VIM; p++)
	  {
		  fv->div_vd += fv->grad_vd[p][p];
	  }
  } else if ( zero_unused_grads && upd->vp[VORT_DIR1] == -1 ) 
	   fv->div_vd = 0.0;

/*
 * grad_sh_K
 * Gradient of curvature in structural shell
 */

  if ( pd->v[SHELL_CURVATURE] )
    {
      v = SHELL_CURVATURE;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
        {
          fv->grad_sh_K[p] = 0.0;
          for ( i=0; i<dofs; i++)
            {
              fv->grad_sh_K[p] += *esp->sh_K[i] * bf[v]->grad_phi[i] [p];
            }
        }
    } else if ( zero_unused_grads && upd->vp[SHELL_CURVATURE] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_K[p] = 0.0;
  }

  if ( pd->v[SHELL_CURVATURE2] )
    {
      v = SHELL_CURVATURE2;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
        {
          fv->grad_sh_K2[p] = 0.0;
          for ( i=0; i<dofs; i++)
            {
              fv->grad_sh_K2[p] += *esp->sh_K2[i] * bf[v]->grad_phi[i] [p];
            }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_CURVATURE2] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_K2[p] = 0.0;
  }  
  
  
  /*
   * grad(n_dot_curl_s_v)
   *   This is the normal gradient of a scalar field defined on a shell.
   */
  if (pd->v[N_DOT_CURL_V])
    {
      v = N_DOT_CURL_V;
      bfn = bf[v];
      dofs  = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_n_dot_curl_s_v[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_n_dot_curl_s_v[p] += *esp->n_dot_curl_s_v[i] * bfn->grad_phi[i] [p];
	    }
	}
    } else if (zero_unused_grads &&  upd->vp[N_DOT_CURL_V] == -1 ) {
      for (p = 0; p < VIM; p++) {
	fv->grad_n_dot_curl_s_v[p] = 0.0;
      }
    }

  /*
   * grad(div_s_v)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the surface divergence of the velocity
   */
  if (pd->v[SHELL_SURF_DIV_V])
    {
      v = SHELL_SURF_DIV_V;
      bfn = bf[v];
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_div_s_v[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_div_s_v[p] += *esp->div_s_v[i] * bfn->grad_phi[i][p];
	    }
	}
    } else if (zero_unused_grads &&  upd->vp[SHELL_SURF_DIV_V] == -1 ) {
      for (p = 0; p < VIM; p++) {
	fv->grad_div_s_v[p] = 0.0;
      }
    }

 /*
   * grad(curv)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the mean curvature of the surface
   */
  if (pd->v[SHELL_SURF_CURV])
    {
      v = SHELL_SURF_CURV;
      bfn = bf[v];
      dofs = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_curv[p] = 0.0;
	  for (i = 0; i < dofs; i++)
	    {
	      fv->grad_curv[p] += *esp->curv[i] * bfn->grad_phi[i][p];
	    }
	}
    } else if (zero_unused_grads &&  upd->vp[SHELL_SURF_CURV] == -1) {
      for (p = 0; p < VIM; p++) {
	fv->grad_curv[p] = 0.0;
      }
    }

  /*
   * grad(grad_s_v_dot_n)
   *   This is the normal gradient of a scalar field defined on a shell.
   *        The scalar field is the vector component of the grad_s_v_dot_n
   */
  if (pd->v[GRAD_S_V_DOT_N1])
    {
      for (r = 0; r < dim; r++) 
	{
	  v = GRAD_S_V_DOT_N1 + r;
	  bfn = bf[v];
	  dofs = ei->dof[v];
	  for (p = 0; p < VIM; p++)
	    {
	      fv->serialgrad_grad_s_v_dot_n[r][p] = 0.0;
	      for (i = 0; i < dofs; i++)
		{
		  fv->serialgrad_grad_s_v_dot_n[r][p] += *(esp->grad_v_dot_n[r][i]) * bfn->grad_phi[i][p];
		}
	    }
	}
    } 
  else  if (zero_unused_grads &&  upd->vp[GRAD_S_V_DOT_N1] == -1 ) 
    {
      for (r = 0; r < dim; r++) 
	{
	  for (p = 0; p < VIM; p++) {
	    fv->serialgrad_grad_s_v_dot_n[r][p] = 0.0;
	  }
	}
    }



  
  /*
   * grad(sh_J)
   */
  if (pd->v[SHELL_DIFF_FLUX])
    {
      v = SHELL_DIFF_FLUX;
      dofs  = ei->dof[v];
      for (p = 0; p < VIM; p++)
	{
	  fv->grad_sh_J[p] = 0.0;
          for (i = 0; i < dofs; i++)
	    {
	      // HKM -> I changed dphidxi to grad_phi below. There didn't seem
	      //        to be a case where the raw element derivative (dphidxi) which
	      //        doesn't even take into account of the size of the element should
	      //        be used at this level. 
	      //      fv->grad_sh_J[p] += *esp->sh_J[i] * bf[v]->dphidxi[i][p];
	      fv->grad_sh_J[p] += *esp->sh_J[i] * bf[v]->grad_phi[i][p];
	    }
	}
    }
  else if ( zero_unused_grads &&  upd->vp[SHELL_DIFF_FLUX] == -1 ){
    for (p=0; p<VIM; p++) fv->grad_sh_J[p] = 0.0;
  }
  
  /*
   * grad(APR)
   */
 
  if ( pd->v[ACOUS_PREAL] )
    {
      v = ACOUS_PREAL;
      dofs  = ei->dof[v];
      
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_apr[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_apr[p] += *esp->apr[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[ACOUS_PREAL] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_apr[p] = 0.0;
  }
  
  
  
  
  if ( pd->v[ACOUS_PIMAG] )
    {
      v = ACOUS_PIMAG;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_api[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_api[p] += *esp->api[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[ACOUS_PIMAG] == -1 ){
    for (p=0; p<VIM; p++) fv->grad_api[p] = 0.0;
  }
  
  
  
  
  if ( pd->v[ACOUS_REYN_STRESS] )
    {
      v = ACOUS_REYN_STRESS;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_ars[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_ars[p] += *esp->ars[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[ACOUS_REYN_STRESS] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_ars[p] = 0.0;
  }
  
  if ( pd->v[SHELL_BDYVELO] )
    {
      v = SHELL_BDYVELO;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_bv[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_bv[p] += *esp->sh_bv[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_BDYVELO] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_bv[p] = 0.0;
  }
  
  if ( pd->v[SHELL_LUBP] )
    {
      v = SHELL_LUBP;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_p[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_p[p] += *esp->sh_p[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_LUBP] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_p[p] = 0.0;
  }
  
  if ( pd->v[LUBP] )
    {
      v = LUBP;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_lubp[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_lubp[p] += *esp->lubp[i] * bf[v]->grad_phi[i] [p];
	      fv_old->grad_lubp[p] += *esp_old->lubp[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  else if ( zero_unused_grads && upd->vp[LUBP] == -1 ) 
    {
      for (p=0; p<VIM; p++) fv->grad_lubp[p] = 0.0;
    }

  if ( pd->v[LUBP_2])
    {
      v = LUBP_2;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_lubp_2[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_lubp_2[p] += *esp->lubp_2[i] * bf[v]->grad_phi[i] [p];
	      fv_old->grad_lubp_2[p] += *esp_old->lubp_2[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  else if ( zero_unused_grads && upd->vp[LUBP_2] == -1 ) 
    {
      for (p=0; p<VIM; p++) fv->grad_lubp_2[p] = 0.0;
    }
  
  if ( pd->v[SHELL_LUB_CURV] )
    {
      v = SHELL_LUB_CURV;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_l_curv[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_l_curv[p] += *esp->sh_l_curv[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_LUB_CURV] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_l_curv[p] = 0.0;
  }

   if ( pd->v[SHELL_LUB_CURV_2] )
    {
      v = SHELL_LUB_CURV_2;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_l_curv_2[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_l_curv_2[p] += *esp->sh_l_curv_2[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_LUB_CURV_2] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_l_curv_2[p] = 0.0;
  }

  if ( pd->v[SHELL_PRESS_OPEN] )
    {
      v = SHELL_PRESS_OPEN;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_p_open[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_p_open[p] += *esp->sh_p_open[i] * bf[v]->grad_phi[i] [p];
	      fv_old->grad_sh_p_open[p] += *esp_old->sh_p_open[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } 
  else if ( zero_unused_grads && upd->vp[SHELL_PRESS_OPEN] == -1 ) 
    {
      for (p=0; p<VIM; p++) fv->grad_sh_p_open[p] = 0.0;
    }

   if ( pd->v[SHELL_PRESS_OPEN_2] )
     {
       v = SHELL_PRESS_OPEN_2;
       dofs  = ei->dof[v];
       for ( p=0; p<VIM; p++)
	 {
	   fv->grad_sh_p_open_2[p] = 0.0;
	   for ( i=0; i<dofs; i++)
	     {
	       fv->grad_sh_p_open_2[p] += *esp->sh_p_open_2[i] * bf[v]->grad_phi[i] [p];
	       fv_old->grad_sh_p_open_2[p] += *esp_old->sh_p_open_2[i] * bf[v]->grad_phi[i] [p];
	     }
	 }
     } 
  else if ( zero_unused_grads && upd->vp[SHELL_PRESS_OPEN_2] == -1 ) 
    {
      for (p=0; p<VIM; p++) fv->grad_sh_p_open_2[p] = 0.0;
    }
  
  
  if ( pd->v[SHELL_FILMP] )
    {
      v = SHELL_FILMP;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_fp[p] = 0.0;
          fv_old->grad_sh_fp[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_fp[p] += *esp->sh_fp[i] * bf[v]->grad_phi[i] [p];
	      fv_old->grad_sh_fp[p] += *esp_old->sh_fp[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_FILMP] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_fp[p] = fv_old->grad_sh_fp[p] = 0.0;
  }
  
  if ( pd->v[SHELL_FILMH] )
    {
      v = SHELL_FILMH;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_fh[p] = 0.0;
	  fv_old->grad_sh_fh[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_fh[p] += *esp->sh_fh[i] * bf[v]->grad_phi[i] [p];
	      fv_old->grad_sh_fh[p] += *esp_old->sh_fh[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_FILMH] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_fh[p] = fv_old->grad_sh_fh[p] = 0.0;
  }
  
  if ( pd->v[SHELL_PARTC] )
    {
      v = SHELL_PARTC;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_pc[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_pc[p] += *esp->sh_pc[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_PARTC] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_pc[p] = 0.0;
  } 

  if ( pd->v[SHELL_TEMPERATURE] )
    {
      v = SHELL_TEMPERATURE;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_sh_t[p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_sh_t[p] += *esp->sh_t[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[SHELL_TEMPERATURE] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_sh_t[p] = 0.0;
  } 

  if ( pd->v[LIGHT_INTP] )
    {
      v = LIGHT_INTP;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_poynt[0][p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_poynt[0][p] += *esp->poynt[0][i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[LIGHT_INTP] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_poynt[0][p] = 0.0;
  } 

  if ( pd->v[LIGHT_INTM] )
    {
      v = LIGHT_INTM;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_poynt[1][p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_poynt[1][p] += *esp->poynt[1][i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[LIGHT_INTM] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_poynt[1][p] = 0.0;
  } 


  if ( pd->v[LIGHT_INTD] )
    {
      v = LIGHT_INTD;
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv->grad_poynt[2][p] = 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv->grad_poynt[2][p] += *esp->poynt[2][i] * bf[v]->grad_phi[i] [p];
	    }
	}
    } else if ( zero_unused_grads && upd->vp[LIGHT_INTD] == -1 ) {
    for (p=0; p<VIM; p++) fv->grad_poynt[2][p] = 0.0;
  } 

 /*
  * External 
  */
 if ( efv->ev )
   {
     v = EXTERNAL;
     for ( w=0; w<efv->Num_external_field; w++)
       {
	 dofs     = ei->dof_ext[w];
	 if( efv->i[w] != I_TABLE )
	   {
	     if (strcmp(efv->name[w], "F1") == 0 && !pd->v[PHASE1])
	       {
		 /* load up the gradient of the the phase function variables 
		  * Currently, this is the only field whose gradient is computed when
		  * is applied as externally
		  */
		 grad_scalar_fv_fill( evp->external_field[w], bfex[v]->grad_phi, dofs, fv->grad_pF[0] );
	       }
             else
               {
                for ( p=0; p<VIM; p++)
                   {
                    fv->grad_ext_field[w][p] = 0.0;
                    for ( i=0; i<dofs; i++)
                       {
                        fv->grad_ext_field[w][p] += *evp->external_field[w][i] * bfex[w]->grad_phi[i][p];
                       }
                   }
               }
	   }
       }
   }

 zero_unused_grads = FALSE; 
  return(status);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int 
load_fv_mesh_derivs(int okToZero)

/*******************************************************************************
 *
 * load_fv_mesh_derivs() -- ld mesh derivatives of field var grads at gauss pt
 *
 * input: (assume the appropriate parts of esp, bf, ei, pd are filled in)
 * ------
 *	
 *
 * output: ( this routine fills in the following parts of the fv structure)
 * -------
 *	fv->
 *  		d_grad_T_dmesh[p] [b][j] 	d( e_p . grad(T) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_P_dmesh[p] [b][j] 	d( e_p . grad(P) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_c_dmesh[p][w] [b][j] 	d( e_p . grad(c_w) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_p_liq_dmesh[p] [b][j] 	d( e_p . grad(p_liq) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_div_v_dmesh [b][j] 		d( grad.v ) 
 * 						-----------
 * 						d ( d_b,j )
 *
 *		d_grad_v_dmesh[p][q] [b][j] 	d( e_p e_q : grad(v) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_X_dmesh[p][q] [b][j] 	d( e_p e_q : grad(X) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_pv_dmesh[p][q] [b][j] 	d( e_p e_q : grad(pv) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_div_d_dmesh [b][j] 		d( grad.d ) 
 * 						-----------
 * 						d ( d_b,j )
 *
 *		d_grad_d_dmesh[p][q] [b][j] 	d( e_p e_q : grad(d) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_S_dmesh[p][q][r] [b][j] 	d( e_p e_q e_r : grad(S) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_G_dmesh[p][q][r] [b][j] 	d( e_p e_q e_r : grad(G) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_Gt_dmesh[p][q][r] [b][j] d( e_p e_q e_r : grad(G)^t ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 *		d_grad_vd_dmesh[p][q] [b][j] 	d( e_p e_q : grad(vd) ) 
 * 						------------------
 * 						    d ( d_b,j )
 *
 * NOTE: To accommodate shell elements, this function has been modified
 *       so that fv variables are not zeroed out when they are active
 *       on an element block other than the current one.
 *       The check done for variable v is then:
 *          if ( pd->v[v] || upd->vp[v] == -1 )
 *       In many cases below, this conditional zeroing is done in
 *       a separate small loop before the main one.
 *
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 *
 * Notes:
 * ------
 *		[1] There is an extra Langrangian contribution for the mesh
 *		    displacment unknowns, in addition to the parts that the
 *		    have to do with spatial gradient operators over a
 *		    deforming mesh that all variables have.
 *
 *
 * Created:	Fri Mar 18 10:18:03 MST 1994 pasacki@sandia.gov
 *
 * Modified:	Tue Feb 21 11:16 MST 1995 pasacki@sandia.gov
 ****************************************************************************/
{
  int v;			/* variable type indicator */
  int i, j;			/* index */
  int b;
  int p, q, r, s;		/* dimension index */
  int vdofs = 0;		/* degrees of freedom for var  in the elem */
  int mdofs;			/* degrees of freedom for mesh in the elem */
  int w;			/* concentration species */
  int status;
  int wim;                      /* Number of components in the velocity vector */
  int dimNonSym;                /* # of dimensions that don't have symmetry */
  int dim;                      /* # dimensions in the physical mesh */
  int siz;
  int mode;                     /* modal counter */
  int transient_run, discontinuous_porous_media;
 
#ifdef DEBUG_HKM
  int doExtra = 0;
#endif

  dbl T_i, v_ri, P_i, d_ri, F_i, d_dot_ri;

  struct Basis_Functions *bfm, *bfv;	/* For mesh variables. */
  
  static int is_initialized = FALSE;
  int VIMis3;

  status = 0;
  
  VIMis3 = (VIM == 3) ? TRUE : FALSE ;

#ifdef DEBUG_HKM
  if ((pd->v[SHELL_NORMAL1] == 0) && ei->dof[SHELL_NORMAL1] > 0) {
    doExtra = 1;
  }
#endif

  /*
   * If this is not a deforming mesh problem, then
   * leave now.
   */
  if (!ei->deforming_mesh)
    {
      return(status);
    }

  dim       =  pd->Num_Dim;
  dimNonSym =  pd->Num_Dim;
  if( (pd->CoordinateSystem == CARTESIAN) ||
      (pd->CoordinateSystem == CYLINDRICAL) )
    {
      wim = pd->Num_Dim;
    }
  else if ( (pd->CoordinateSystem == SWIRLING) ||
	    (pd->CoordinateSystem == PROJECTED_CARTESIAN) )
    {
      wim = 3;
    }
  else
    {
      wim = VIM;
    }


  mdofs = ei->dof[R_MESH1];

  bfm = bf[R_MESH1];
  
  transient_run = pd->TimeIntegration != STEADY ; 
  discontinuous_porous_media = mp->PorousMediaType !=CONTINUOUS ;

  /*
   * d(grad(T))/dmesh
   */

  if ( pd->v[TEMPERATURE] )
    {
      v = TEMPERATURE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_T_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->T[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for ( b=0; b < dim; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      fv->d_grad_T_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  T_i = *esp->T[0];

	  fv->d_grad_T_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_T_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_T_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_T_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_T_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_T_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_T_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_T_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_T_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
	      T_i = *esp->T[i];
			
	      fv->d_grad_T_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_T_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_T_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_T_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_T_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_T_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_T_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_T_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_T_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if (upd->vp[TEMPERATURE] != -1 && !is_initialized && okToZero) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_T_dmesh[0][0][0]),0, siz);
    }
  
  

  /*
   * d(grad(V))/dmesh
   */
  if ( pd->v[VOLTAGE] )
    {
      v = VOLTAGE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_V_dmesh[0][0][0]),0, siz);
      for (p = 0; p < dimNonSym; p++)
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for ( i=0; i<vdofs; i++)
		    {
		      fv->d_grad_V_dmesh[p] [b][j] +=
			*esp->V[i]  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_V_dmesh[0] [0][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_V_dmesh[1] [1][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_V_dmesh[1] [0][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_V_dmesh[0] [1][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_V_dmesh[2] [2][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_V_dmesh[2] [0][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_V_dmesh[2] [1][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_V_dmesh[0] [2][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_V_dmesh[1] [2][j] = *esp->V[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_V_dmesh[0] [0][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_V_dmesh[1] [1][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_V_dmesh[1] [0][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_V_dmesh[0] [1][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_V_dmesh[2] [2][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_V_dmesh[2] [0][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_V_dmesh[2] [1][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_V_dmesh[0] [2][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_V_dmesh[1] [2][j] += *esp->V[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[VOLTAGE] != -1  && !is_initialized  && okToZero) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_V_dmesh[0][0][0]),0, siz);
    }
  
  /*
   * d(grad(sh_K))/dmesh
   */

 if ( pd->v[SHELL_CURVATURE] )
    {
      v = SHELL_CURVATURE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_K_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
        {
          T_i = *esp->sh_K[i];
          for (p = 0; p < dimNonSym; p++)
            {
              for (b = 0; b < dim; b++)
                {
                  for (j = 0; j < mdofs; j++)
                    {
                      fv->d_grad_sh_K_dmesh[p] [b][j] +=
                        T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
                    }
                }
            }
	}
#else
      for ( j=0; j<mdofs; j++)
        {
          T_i = *esp->sh_K[0];

          fv->d_grad_sh_K_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
          fv->d_grad_sh_K_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
          fv->d_grad_sh_K_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
          fv->d_grad_sh_K_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];

          if (dimNonSym == 3)
            {
              fv->d_grad_sh_K_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
              fv->d_grad_sh_K_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
              fv->d_grad_sh_K_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
              fv->d_grad_sh_K_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
              fv->d_grad_sh_K_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
            }

          for (i = 1; i < vdofs; i++)
            {
              T_i = *esp->sh_K[i];

              fv->d_grad_sh_K_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
              fv->d_grad_sh_K_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
              fv->d_grad_sh_K_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
              fv->d_grad_sh_K_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];

              if (dimNonSym == 3)
                {
                  fv->d_grad_sh_K_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
                  fv->d_grad_sh_K_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
                  fv->d_grad_sh_K_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
                  fv->d_grad_sh_K_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
                  fv->d_grad_sh_K_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
                }
            }
	}
#endif
    } else   if ( upd->vp[SHELL_CURVATURE] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_K_dmesh[0][0][0]),0, siz);
    }

 if ( pd->v[SHELL_CURVATURE2] )
    {
      v = SHELL_CURVATURE2;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_K2_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
        {
          T_i = *esp->sh_K2[i];
          for (p = 0; p < dimNonSym; p++)
            {
              for (b = 0; b < dim; b++)
                {
                  for (j = 0; j < mdofs; j++)
                    {
                      fv->d_grad_sh_K2_dmesh[p] [b][j] +=
                        T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
                    }
                }
            }
	}
#else
      for ( j=0; j<mdofs; j++)
        {
          T_i = *esp->sh_K2[0];

          fv->d_grad_sh_K2_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
          fv->d_grad_sh_K2_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
          fv->d_grad_sh_K2_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
          fv->d_grad_sh_K2_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];

          if (dimNonSym == 3)
            {
              fv->d_grad_sh_K2_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
              fv->d_grad_sh_K2_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
              fv->d_grad_sh_K2_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
              fv->d_grad_sh_K2_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
              fv->d_grad_sh_K2_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
            }

          for (i = 1; i < vdofs; i++)
            {
              T_i = *esp->sh_K2[i];

              fv->d_grad_sh_K2_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
              fv->d_grad_sh_K2_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
              fv->d_grad_sh_K2_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
              fv->d_grad_sh_K2_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];

              if (dimNonSym == 3)
                {
                  fv->d_grad_sh_K2_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
                  fv->d_grad_sh_K2_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
                  fv->d_grad_sh_K2_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
                  fv->d_grad_sh_K2_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
                  fv->d_grad_sh_K2_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
                }
            }
	}
#endif
    } else   if ( upd->vp[SHELL_CURVATURE2] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_K2_dmesh[0][0][0]),0, siz);
    }


  /*
   * d(grad(apr))/dmesh
   */
  if ( pd->v[ACOUS_PREAL] )
    {
      v = ACOUS_PREAL;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_apr_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->apr[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      fv->d_grad_apr_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->apr[0];

	  fv->d_grad_apr_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_apr_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_apr_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_apr_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_apr_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_apr_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_apr_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_apr_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_apr_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->apr[i];
			
	      fv->d_grad_apr_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_apr_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_apr_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_apr_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_apr_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_apr_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_apr_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_apr_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_apr_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else if ( upd->vp[ACOUS_PREAL] != -1   && okToZero)  {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_apr_dmesh[0][0][0]),0, siz);
    }



  if ( pd->v[ACOUS_PIMAG] )
    {
      v = ACOUS_PIMAG;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_api_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->api[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      fv->d_grad_api_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->api[0];

	  fv->d_grad_api_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_api_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_api_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_api_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_api_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_api_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_api_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_api_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_api_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->api[i];
			
	      fv->d_grad_api_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_api_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_api_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_api_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_api_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_api_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_api_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_api_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_api_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[ACOUS_PIMAG] != -1  && okToZero) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_api_dmesh[0][0][0]),0, siz);
    }
  



  if ( pd->v[ACOUS_REYN_STRESS])
    {
      v = ACOUS_REYN_STRESS;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_ars_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->ars[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_ars_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->ars[0];

	  fv->d_grad_ars_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_ars_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_ars_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_ars_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_ars_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_ars_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_ars_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_ars_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_ars_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->ars[i];
			
	      fv->d_grad_ars_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_ars_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_ars_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_ars_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_ars_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_ars_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_ars_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_ars_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_ars_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[ACOUS_REYN_STRESS] != -1   && okToZero){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_ars_dmesh[0][0][0]),0, siz);
    }

  if (pd->v[SHELL_BDYVELO])
    {
      v = SHELL_BDYVELO;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_bv_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_bv[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_bv_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_bv[0];

	  fv->d_grad_sh_bv_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_bv_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_bv_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_bv_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_bv_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_bv_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_bv_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_bv_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_bv_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_bv[i];
			
	      fv->d_grad_sh_bv_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_bv_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_bv_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_bv_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_bv_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_bv_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_bv_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_bv_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_bv_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_BDYVELO] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_bv_dmesh[0][0][0]),0, siz);
    }


 if ( pd->v[LUBP] )
    {
      v = LUBP;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_lubp_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->lubp[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_lubp_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->lubp[0];

	  fv->d_grad_lubp_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_lubp_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_lubp_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_lubp_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_lubp_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_lubp_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_lubp_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_lubp_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_lubp_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->lubp[i];
			
	      fv->d_grad_lubp_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_lubp_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_lubp_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_lubp_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_lubp_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_lubp_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_lubp_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_lubp_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_lubp_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[LUBP] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_lubp_dmesh[0][0][0]),0, siz);
    }

if ( pd->v[LUBP_2] )
    {
      v = LUBP_2;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_lubp_2_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->lubp_2[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_lubp_2_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->lubp_2[0];

	  fv->d_grad_lubp_2_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_lubp_2_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_lubp_2_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_lubp_2_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_lubp_2_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_lubp_2_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_lubp_2_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_lubp_2_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_lubp_2_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->lubp_2[i];
			
	      fv->d_grad_lubp_2_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_lubp_2_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_lubp_2_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_lubp_2_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_lubp_2_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_lubp_2_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_lubp_2_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_lubp_2_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_lubp_2_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[LUBP_2] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_lubp_2_dmesh[0][0][0]),0, siz);
    }

 if ( pd->v[SHELL_PRESS_OPEN] )
    {
      v = SHELL_PRESS_OPEN;
      bfv = bf[v];
#ifdef DO_NOT_UNROLL
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_p_open_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_p_open[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_p_open_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_p_open[0];

	  fv->d_grad_sh_p_open_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_p_open_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_p_open_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_p_open_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_p_open_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_p_open_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_p_open_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_p_open_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_p_open_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_p_open[i];
			
	      fv->d_grad_sh_p_open_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_p_open_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_p_open_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_p_open_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_p_open_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_p_open_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_p_open_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_p_open_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_p_open_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_PRESS_OPEN] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_p_open_dmesh[0][0][0]),0, siz);
    }

 if ( pd->v[SHELL_PRESS_OPEN_2] )
    {
      v = SHELL_PRESS_OPEN_2;
      bfv = bf[v];
#ifdef DO_NOT_UNROLL
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_p_open_2_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_p_open_2[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_p_open_2_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_p_open_2[0];

	  fv->d_grad_sh_p_open_2_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_p_open_2_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_p_open_2_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_p_open_2_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_p_open_2_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_p_open_2_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_p_open_2_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_p_open_2_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_p_open_2_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_p_open_2[i];
			
	      fv->d_grad_sh_p_open_2_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_p_open_2_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_p_open_2_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_p_open_2_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_p_open_2_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_p_open_2_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_p_open_2_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_p_open_2_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_p_open_2_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_PRESS_OPEN_2] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_p_open_2_dmesh[0][0][0]),0, siz);
    }


if ( pd->v[SHELL_FILMP] )
    {
      v = SHELL_FILMP;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_fp_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_fp[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_fp_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_fp[0];

	  fv->d_grad_sh_fp_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_fp_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_fp_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_fp_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_fp_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_fp_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_fp_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_fp_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_fp_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_fp[i];
			
	      fv->d_grad_sh_fp_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_fp_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_fp_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_fp_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_fp_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_fp_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_fp_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_fp_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_fp_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_FILMP] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_fp_dmesh[0][0][0]),0, siz);
    }

 if ( pd->v[SHELL_FILMH] )
    {
      v = SHELL_FILMH;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_fh_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_fh[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_fh_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_fh[0];

	  fv->d_grad_sh_fh_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_fh_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_fh_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_fh_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_fh_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_fh_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_fh_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_fh_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_fh_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_fh[i];
			
	      fv->d_grad_sh_fh_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_fh_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_fh_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_fh_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_fh_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_fh_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_fh_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_fh_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_fh_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_FILMH] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_fh_dmesh[0][0][0]),0, siz);
    }

 
  if ( pd->v[SHELL_PARTC] )
    {
      v = SHELL_PARTC;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_pc_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_pc[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_pc_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_pc[0];

	  fv->d_grad_sh_pc_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_pc_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_pc_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_pc_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_pc_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_pc_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_pc_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_pc_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_pc_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_pc[i];
			
	      fv->d_grad_sh_pc_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_pc_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_pc_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_pc_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_pc_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_pc_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_pc_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_pc_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_pc_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_PARTC] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_pc_dmesh[0][0][0]),0, siz);
    } 

if ( pd->v[SHELL_TEMPERATURE] )
    {
      v = SHELL_TEMPERATURE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_t_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_t[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_t_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_t[0];

	  fv->d_grad_sh_t_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_t_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_t_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_t_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_t_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_t_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_t_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_t_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_t_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_t[i];
			
	      fv->d_grad_sh_t_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_t_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_t_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_t_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_t_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_t_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_t_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_t_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_t_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_TEMPERATURE] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_t_dmesh[0][0][0]),0, siz);
    }

if ( pd->v[LIGHT_INTP] )
    {
      v = LIGHT_INTP;
      bfv = bf[v];
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[0][0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->poynt[0][i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_poynt_dmesh[0][p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
    } else   if ( upd->vp[LIGHT_INTP] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[0][0][0][0]),0, siz);
    }

if ( pd->v[LIGHT_INTM] )
    {
      v = LIGHT_INTM;
      bfv = bf[v];
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[1][0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->poynt[1][i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_poynt_dmesh[1][p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
    } else   if ( upd->vp[LIGHT_INTM] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[1][0][0][0]),0, siz);
    }

if ( pd->v[LIGHT_INTD] )
    {
      v = LIGHT_INTD;
      bfv = bf[v];
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[2][0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->poynt[2][i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_poynt_dmesh[2][p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
    } else   if ( upd->vp[LIGHT_INTD] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_poynt_dmesh[2][0][0][0]),0, siz);
    }

if ( pd->v[SHELL_LUB_CURV] )
    {
      v = SHELL_LUB_CURV;
      bfv = bf[v];
#ifdef DO_NOT_UNROLL
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_l_curv_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_l_curv[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_l_curv_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_l_curv[0];

	  fv->d_grad_sh_l_curv_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_l_curv_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_l_curv_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_l_curv_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_l_curv_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_l_curv_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_l_curv_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_l_curv_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_l_curv_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_l_curv[i];
			
	      fv->d_grad_sh_l_curv_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_l_curv_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_l_curv_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_l_curv_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_l_curv_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_l_curv_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_l_curv_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_l_curv_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_l_curv_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_LUB_CURV] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_l_curv_dmesh[0][0][0]),0, siz);
    }

if ( pd->v[SHELL_LUB_CURV_2] )
    {
      v = SHELL_LUB_CURV_2;
      bfv = bf[v];
#ifdef DO_NOT_UNROLL
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_l_curv_2_dmesh[0][0][0]),0, siz);
      for ( i=0; i<vdofs; i++)
	{
	  T_i = *esp->sh_l_curv_2[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_sh_l_curv_2_dmesh[p] [b][j] +=
 			T_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
          T_i = *esp->sh_l_curv_2[0];

	  fv->d_grad_sh_l_curv_2_dmesh[0] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_l_curv_2_dmesh[1] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_l_curv_2_dmesh[1] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_l_curv_2_dmesh[0] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_l_curv_2_dmesh[2] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_l_curv_2_dmesh[2] [0][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_l_curv_2_dmesh[2] [1][j] = T_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_l_curv_2_dmesh[0] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_l_curv_2_dmesh[1] [2][j] = T_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
              T_i = *esp->sh_l_curv[i];
			
	      fv->d_grad_sh_l_curv_2_dmesh[0] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_l_curv_2_dmesh[1] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_l_curv_2_dmesh[1] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_l_curv_2_dmesh[0] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_l_curv_2_dmesh[2] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_l_curv_2_dmesh[2] [0][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_l_curv_2_dmesh[2] [1][j] += T_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_l_curv_2_dmesh[0] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_l_curv_2_dmesh[1] [2][j] += T_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if ( upd->vp[SHELL_LUB_CURV_2] != -1  ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_l_curv_2_dmesh[0][0][0]),0, siz);
    }


  /*
   * d(grad(qs))/dmesh
   */

  if (pd->v[SURF_CHARGE])
    {
      v = SURF_CHARGE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_qs_dmesh[0][0][0]),0, siz);
      for (p = 0; p < dimNonSym; p++)
	{
	  for ( b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for ( i=0; i<vdofs; i++)
		    {
		      fv->d_grad_qs_dmesh[p] [b][j] +=
			*esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	} 
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_qs_dmesh[0] [0][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_qs_dmesh[1] [1][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_qs_dmesh[1] [0][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_qs_dmesh[0] [1][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_qs_dmesh[2] [2][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_qs_dmesh[2] [0][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_qs_dmesh[2] [1][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_qs_dmesh[0] [2][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_qs_dmesh[1] [2][j] = *esp->qs[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_qs_dmesh[0] [0][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_qs_dmesh[1] [1][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_qs_dmesh[1] [0][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_qs_dmesh[0] [1][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_qs_dmesh[2] [2][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_qs_dmesh[2] [0][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_qs_dmesh[2] [1][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_qs_dmesh[0] [2][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_qs_dmesh[1] [2][j] += *esp->qs[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else  if ( upd->vp[SURF_CHARGE] != -1  && okToZero ) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_qs_dmesh[0][0][0]),0, siz);
    }

  /*
   * d(grad(sh_J))/dmesh
   */

  /* This is needed for bulk assembly of shell diffusion KBC */
  if (ei->dof[SHELL_DIFF_FLUX] > 0)
    {
      v = SHELL_DIFF_FLUX;
      bfv = bf[v];
#ifdef DO_NOT_UNROLL
      vdofs  = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_J_dmesh[0][0][0]),0, siz);
      for (p = 0; p < dimNonSym; p++)
	{
	  for (b = 0; b < dim; b++)
	    {
	      for ( j=0; j<mdofs; j++)
		{
		  for ( i=0; i<vdofs; i++)
		    {
		      fv->d_grad_sh_J_dmesh[p] [b][j] +=
			*esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	} 
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_sh_J_dmesh[0] [0][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_sh_J_dmesh[1] [1][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_sh_J_dmesh[1] [0][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_sh_J_dmesh[0] [1][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_sh_J_dmesh[2] [2][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_sh_J_dmesh[2] [0][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_sh_J_dmesh[2] [1][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_sh_J_dmesh[0] [2][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_sh_J_dmesh[1] [2][j] = *esp->sh_J[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_sh_J_dmesh[0] [0][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_sh_J_dmesh[1] [1][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_sh_J_dmesh[1] [0][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_sh_J_dmesh[0] [1][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_sh_J_dmesh[2] [2][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_sh_J_dmesh[2] [0][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_sh_J_dmesh[2] [1][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_sh_J_dmesh[0] [2][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_sh_J_dmesh[1] [2][j] += *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else  if ( upd->vp[SHELL_DIFF_FLUX] != -1  && okToZero ) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_J_dmesh[0][0][0]),0, siz);
    }

  /*
   * d(grad(SH))/dmesh
   */
  if (pd->v[SHEAR_RATE])
    {
      v = SHEAR_RATE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_SH_dmesh,0, siz);
      
      for ( p = 0; p < dimNonSym; p++)
	{
	  for (b = 0; b < dim; b++)
	    {
	      for ( j=0; j<mdofs; j++)
		{
		  for ( i=0; i<vdofs; i++)
		    {
		      fv->d_grad_SH_dmesh[p] [b][j] +=
			*esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_SH_dmesh[0] [0][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_SH_dmesh[1] [1][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_SH_dmesh[1] [0][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_SH_dmesh[0] [1][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_SH_dmesh[2] [2][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_SH_dmesh[2] [0][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_SH_dmesh[2] [1][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_SH_dmesh[0] [2][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_SH_dmesh[1] [2][j] = *esp->SH[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_SH_dmesh[0] [0][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_SH_dmesh[1] [1][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_SH_dmesh[1] [0][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_SH_dmesh[0] [1][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_api_dmesh[2] [2][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_api_dmesh[2] [0][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_api_dmesh[2] [1][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_api_dmesh[0] [2][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_api_dmesh[1] [2][j] += *esp->SH[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else   if (upd->vp[SHEAR_RATE] != -1  && okToZero ){
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_SH_dmesh,0, siz);
    }
  


  /*
   * d(grad(F))/dmesh
   */
  if (pd->v[FILL])
    {
      v = FILL;
      vdofs  = ei->dof[v];
#ifdef DO_NO_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_F_dmesh,0, siz);
      for (p = 0; p < dimNonSym; p++)
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_F_dmesh[p] [b][j] +=
			*esp->F[i]  *  bf[v]->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  F_i = *esp->F[0];

	  fv->d_grad_F_dmesh[0] [0][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_F_dmesh[1] [1][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_F_dmesh[1] [0][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_F_dmesh[0] [1][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_F_dmesh[2] [2][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_F_dmesh[2] [0][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_F_dmesh[2] [1][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_F_dmesh[0] [2][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_F_dmesh[1] [2][j] = F_i  *  bf[v]->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {
	      F_i = *esp->F[i];
			
	      fv->d_grad_F_dmesh[0] [0][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_F_dmesh[1] [1][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_F_dmesh[1] [0][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_F_dmesh[0] [1][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_F_dmesh[2] [2][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_F_dmesh[2] [0][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_F_dmesh[2] [1][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_F_dmesh[0] [2][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_F_dmesh[1] [2][j] += F_i  *  bf[v]->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
  
	  
    } else if (upd->vp[FILL] != -1) 
      {
	siz = sizeof(double)*DIM*DIM*MDE;
	memset(fv->d_grad_F_dmesh,0, siz);
      }
  
  /*
   * d(grad(P))/dmesh
   */

  if (pd->v[PRESSURE])
    {
      v = PRESSURE;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NO_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_P_dmesh,0, siz);
      for (i = 0; i < vdofs; i++)
	{
	  P_i = *esp->P[i];
	  for (p = 0; p < dimNonSym; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_P_dmesh[p] [b][j] +=
 			P_i  *  bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for (j = 0; j < mdofs; j++)
	{
	  P_i = *esp->P[0]; 
	  fv->d_grad_P_dmesh[0] [0][j] = P_i  *  bfv->d_grad_phi_dmesh[0][0] [0][j]; 
	  fv->d_grad_P_dmesh[1] [1][j] = P_i  *  bfv->d_grad_phi_dmesh[0][1] [1][j]; 
	  fv->d_grad_P_dmesh[1] [0][j] = P_i  *  bfv->d_grad_phi_dmesh[0][1] [0][j]; 
	  fv->d_grad_P_dmesh[0] [1][j] = P_i  *  bfv->d_grad_phi_dmesh[0][0] [1][j]; 
		  
	  if (dimNonSym == 3)
	    {
	      fv->d_grad_P_dmesh[2] [2][j] = P_i  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_P_dmesh[2] [0][j] = P_i  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_P_dmesh[2] [1][j] = P_i  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_P_dmesh[1] [2][j] = P_i  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	      fv->d_grad_P_dmesh[0] [2][j] = P_i  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	    }
		  
	  for (i = 1; i < vdofs; i++)
	    {
	      P_i = *esp->P[i];  
	      fv->d_grad_P_dmesh[0] [0][j] += P_i  *  bfv->d_grad_phi_dmesh[i][0] [0][j]; 
	      fv->d_grad_P_dmesh[1] [1][j] += P_i  *  bfv->d_grad_phi_dmesh[i][1] [1][j]; 
	      fv->d_grad_P_dmesh[1] [0][j] += P_i  *  bfv->d_grad_phi_dmesh[i][1] [0][j]; 
	      fv->d_grad_P_dmesh[0] [1][j] += P_i  *  bfv->d_grad_phi_dmesh[i][0] [1][j]; 
			  
	      if (dimNonSym == 3)
		{
		  fv->d_grad_P_dmesh[2] [2][j] += P_i  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_P_dmesh[2] [0][j] += P_i  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_P_dmesh[2] [1][j] += P_i  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_P_dmesh[1] [2][j] += P_i  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		  fv->d_grad_P_dmesh[0] [2][j] += P_i  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		}
	    } 
	}
#endif	  
    } else   if ( upd->vp[PRESSURE] != -1  && okToZero )
      {
	siz = sizeof(double)*DIM*DIM*MDE;
	memset(fv->d_grad_P_dmesh,0, siz);
      }

  /*
   * d(grad(nn))/dmesh
   */

  if ( pd->v[BOND_EVOLUTION] )
    {
      v = BOND_EVOLUTION;
      bfv = bf[v];
      vdofs  = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_nn_dmesh,0, siz);
      for (i = 0; i < vdofs; i++)
	{
	  for (p = 0; p < dim; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_nn_dmesh[p] [b][j] += *esp->nn[i] * bfv->d_grad_phi_dmesh[i][p] [b][j]; 
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_nn_dmesh[0] [0][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_nn_dmesh[1] [1][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_nn_dmesh[1] [0][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_nn_dmesh[0] [1][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_nn_dmesh[2] [2][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_nn_dmesh[2] [0][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_nn_dmesh[2] [1][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_nn_dmesh[0] [2][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_nn_dmesh[1] [2][j] = *esp->nn[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_nn_dmesh[0] [0][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_nn_dmesh[1] [1][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_nn_dmesh[1] [0][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_nn_dmesh[0] [1][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_nn_dmesh[2] [2][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_nn_dmesh[2] [0][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_nn_dmesh[2] [1][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_nn_dmesh[0] [2][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_nn_dmesh[1] [2][j] += *esp->nn[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    }

  /*
   * d(grad(c_w))/dmesh
   */
  if (pd->v[MASS_FRACTION])
    {	  
      v = MASS_FRACTION;
      bfv = bf[v];
      vdofs     = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MAX_CONC*MDE;
      memset(fv->d_grad_c_dmesh, 0, siz);
	
      for (w = 0; w < pd->Num_Species_Eqn; w++)
	{
	  for (p = 0; p < dim; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      for (i = 0; i < vdofs; i++)
			{
			  fv->d_grad_c_dmesh[p][w] [b][j] +=
			    *esp->c[w][i] 
			    * 
			    bfv->d_grad_phi_dmesh[i][p] [b][j];
			}
		    }
		}
	    }
	}
#else
      for (w = 0; w < pd->Num_Species_Eqn; w++)
        {
         for ( j=0; j<mdofs; j++)
	    {
	     fv->d_grad_c_dmesh[0][w] [0][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	     fv->d_grad_c_dmesh[1][w] [1][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	     fv->d_grad_c_dmesh[1][w] [0][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	     fv->d_grad_c_dmesh[0][w] [1][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	     if (dimNonSym == 3) 
	       {
	        fv->d_grad_c_dmesh[2][w] [2][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	        fv->d_grad_c_dmesh[2][w] [0][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	        fv->d_grad_c_dmesh[2][w] [1][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	        fv->d_grad_c_dmesh[0][w] [2][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	        fv->d_grad_c_dmesh[1][w] [2][j] = *esp->c[w][0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	       }
		
	     for (i = 1; i < vdofs; i++)
	       {			
	        fv->d_grad_c_dmesh[0][w] [0][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	        fv->d_grad_c_dmesh[1][w] [1][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	        fv->d_grad_c_dmesh[1][w] [0][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	        fv->d_grad_c_dmesh[0][w] [1][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	        if (dimNonSym == 3) 
		  {
		   fv->d_grad_c_dmesh[2][w] [2][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		   fv->d_grad_c_dmesh[2][w] [0][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		   fv->d_grad_c_dmesh[2][w] [1][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		   fv->d_grad_c_dmesh[0][w] [2][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		   fv->d_grad_c_dmesh[1][w] [2][j] += *esp->c[w][i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		  }
	       }
	    }
        }
#endif
    } else if (upd->vp[MASS_FRACTION] != -1  && okToZero )
      {
	siz = sizeof(double)*DIM*DIM*MAX_CONC*MDE;
	memset(fv->d_grad_c_dmesh,0, siz);
      }  

  /*
   * d(grad(porous_media_variables))/dmesh
   */
  if (pd->v[POR_LIQ_PRES])
    {
      v = POR_LIQ_PRES;
      bfv = bf[v];
      vdofs     = ei->dof[v];
#ifdef DO_NO_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_p_liq_dmesh,0, siz);
	  
      for (p = 0; p < dim; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_p_liq_dmesh[p] [b][j] +=
			*esp->p_liq[i]
			*
			bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for (j = 0; j < mdofs; j++)
	{
	  P_i = *esp->p_liq[0];
	  fv->d_grad_p_liq_dmesh[0] [0][j] =  P_i * bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_p_liq_dmesh[1] [1][j] =  P_i * bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_p_liq_dmesh[0] [1][j] =  P_i * bfv->d_grad_phi_dmesh[0][0] [1][j];
	  fv->d_grad_p_liq_dmesh[1] [0][j] =  P_i * bfv->d_grad_phi_dmesh[0][1] [0][j];
		  
	  if (dim == 3)
	    {
	      fv->d_grad_p_liq_dmesh[2] [2][j] =  P_i * bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_p_liq_dmesh[2] [0][j] =  P_i * bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_p_liq_dmesh[2] [1][j] =  P_i * bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_p_liq_dmesh[0] [2][j] =  P_i * bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_p_liq_dmesh[1] [2][j] =  P_i * bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
		  
	  for (i = 1; i < vdofs; i++)
	    {
	      P_i = *esp->p_liq[i];
	      /* fv->d_grad_p_liq_dmesh[p] [b][j] +=  P_i * bfv->d_grad_phi_dmesh[i][p] [b][j];*/
	      fv->d_grad_p_liq_dmesh[0] [0][j] +=  P_i * bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_p_liq_dmesh[1] [1][j] +=  P_i * bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_p_liq_dmesh[0] [1][j] +=  P_i * bfv->d_grad_phi_dmesh[i][0] [1][j];
	      fv->d_grad_p_liq_dmesh[1] [0][j] +=  P_i * bfv->d_grad_phi_dmesh[i][1] [0][j];
			  
	      if (dim == 3)
		{
		  fv->d_grad_p_liq_dmesh[2] [2][j] +=  P_i * bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_p_liq_dmesh[2] [0][j] +=  P_i * bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_p_liq_dmesh[2] [1][j] +=  P_i * bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_p_liq_dmesh[0] [2][j] +=  P_i * bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_p_liq_dmesh[1] [2][j] +=  P_i * bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
			  
	    }
	}
#endif

    } else if ( upd->vp[POR_LIQ_PRES] != -1  && okToZero ) 
      {
	siz = sizeof(double)*DIM*DIM*MDE;
	memset(fv->d_grad_p_liq_dmesh,0, siz);
      }



  if (pd->v[POR_GAS_PRES])
    {
      v = POR_GAS_PRES;
      bfv = bf[v];
      vdofs     = ei->dof[v];
#ifdef DO_NOT_UNROLL
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_p_gas_dmesh,0, siz);
	  
      for (p = 0; p < dim; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_p_gas_dmesh[p] [b][j] +=
			*esp->p_gas[i]
			*
			bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_p_gas_dmesh[0] [0][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_p_gas_dmesh[1] [1][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_p_gas_dmesh[1] [0][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_p_gas_dmesh[0] [1][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_p_gas_dmesh[2] [2][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_p_gas_dmesh[2] [0][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_p_gas_dmesh[2] [1][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_p_gas_dmesh[0] [2][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_p_gas_dmesh[1] [2][j] = *esp->p_gas[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_p_gas_dmesh[0] [0][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_p_gas_dmesh[1] [1][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_p_gas_dmesh[1] [0][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_p_gas_dmesh[0] [1][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_p_gas_dmesh[2] [2][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_p_gas_dmesh[2] [0][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_p_gas_dmesh[2] [1][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_p_gas_dmesh[0] [2][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_p_gas_dmesh[1] [2][j] += *esp->p_gas[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else if (upd->vp[POR_GAS_PRES] != -1  && okToZero )
      {
	siz = sizeof(double)*DIM*DIM*MDE;
	memset(fv->d_grad_p_gas_dmesh,0, siz);
      }


  if (pd->v[POR_POROSITY])
    {
      v = POR_POROSITY;
      bfv = bf[v];
      vdofs     = ei->dof[v];
#ifdef DO_NOT_UNROLL	  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_porosity_dmesh,0, siz);
	  
      for (p = 0; p < dim; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for ( j=0; j<mdofs; j++)
		{
		  for ( i=0; i<vdofs; i++)
		    {
		      fv->d_grad_porosity_dmesh[p] [b][j] +=
			*esp->porosity[i]
			*
			bfv->d_grad_phi_dmesh[i][p] [b][j];
		    }
		}
	    }
	}
#else
      for ( j=0; j<mdofs; j++)
	{
	  fv->d_grad_porosity_dmesh[0] [0][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][0] [0][j];
	  fv->d_grad_porosity_dmesh[1] [1][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][1] [1][j];
	  fv->d_grad_porosity_dmesh[1] [0][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][1] [0][j];
	  fv->d_grad_porosity_dmesh[0] [1][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][0] [1][j];
		
	  if (dimNonSym == 3) 
	    {
	      fv->d_grad_porosity_dmesh[2] [2][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][2] [2][j];
	      fv->d_grad_porosity_dmesh[2] [0][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][2] [0][j];
	      fv->d_grad_porosity_dmesh[2] [1][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][2] [1][j];
	      fv->d_grad_porosity_dmesh[0] [2][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][0] [2][j];
	      fv->d_grad_porosity_dmesh[1] [2][j] = *esp->porosity[0]  *  bfv->d_grad_phi_dmesh[0][1] [2][j];
	    }
		
	  for (i = 1; i < vdofs; i++)
	    {			
	      fv->d_grad_porosity_dmesh[0] [0][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][0] [0][j];
	      fv->d_grad_porosity_dmesh[1] [1][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][1] [1][j];
	      fv->d_grad_porosity_dmesh[1] [0][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][1] [0][j];
	      fv->d_grad_porosity_dmesh[0] [1][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][0] [1][j];
			
	      if (dimNonSym == 3) 
		{
		  fv->d_grad_porosity_dmesh[2] [2][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][2] [2][j];
		  fv->d_grad_porosity_dmesh[2] [0][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][2] [0][j];
		  fv->d_grad_porosity_dmesh[2] [1][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][2] [1][j];
		  fv->d_grad_porosity_dmesh[0] [2][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][0] [2][j];
		  fv->d_grad_porosity_dmesh[1] [2][j] += *esp->porosity[i]  *  bfv->d_grad_phi_dmesh[i][1] [2][j];
		}
	    }
	}
#endif
    } else if ( upd->vp[POR_POROSITY] != -1 && okToZero) {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_porosity_dmesh,0, siz);
    }
  
	      
  /*
   * d( grad(v))/dmesh
   */      
  if (pd->v[VELOCITY1])
    {
      v = VELOCITY1;
      bfv = bf[v];
      vdofs = ei->dof[v];
      
#ifdef DO_NO_UNROLL
      siz = sizeof(double)*DIM*VIM*DIM*MDE;
      memset(fv->d_grad_v_dmesh,0, siz);
	  
      for (r = 0; r < wim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {  
	      v_ri = *esp->v[r][i];  
	      for (p = 0; p < VIM; p++)
		{
		  for (q = 0; q < VIM; q++)
		    {	  
		      for (b = 0; b < dim; b++)
			{
			  for (j = 0; j < mdofs; j++)
			    {					  
			      fv->d_grad_v_dmesh[p][q] [b][j] +=
				v_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];
			    }
			}
		    }
		}
	    }
	}
#else
	  
      v_ri = *esp->v[0][0];	  
	  
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      fv->d_grad_v_dmesh[0][0] [b][j] = 0.0; 
	      fv->d_grad_v_dmesh[1][1] [b][j] = 0.0; 
	      fv->d_grad_v_dmesh[1][0] [b][j] = 0.0; 
	      fv->d_grad_v_dmesh[0][1] [b][j] = 0.0;
	      if (VIMis3)
		{
		  fv->d_grad_v_dmesh[2][2] [b][j] = 0.0;
		  fv->d_grad_v_dmesh[2][0] [b][j] = 0.0;
		  fv->d_grad_v_dmesh[2][1] [b][j] = 0.0;
		  fv->d_grad_v_dmesh[0][2] [b][j] = 0.0;
		  fv->d_grad_v_dmesh[1][2] [b][j] = 0.0;
		}
	    }
	}
      
      for (r = 0; r < wim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {
	      v_ri = *esp->v[r][i];
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_v_dmesh[0][0] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][0] [b][j]; 
		      fv->d_grad_v_dmesh[1][1] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][1] [b][j]; 
		      fv->d_grad_v_dmesh[1][0] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][0] [b][j]; 
		      fv->d_grad_v_dmesh[0][1] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][1] [b][j]; 
		      
		      if (VIMis3) 
			{
			  fv->d_grad_v_dmesh[2][2] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][2] [b][j];
			  fv->d_grad_v_dmesh[2][0] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][0] [b][j];
			  fv->d_grad_v_dmesh[2][1] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][1] [b][j];
			  fv->d_grad_v_dmesh[0][2] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][2] [b][j];
			  fv->d_grad_v_dmesh[1][2] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][2] [b][j];
			}
		    }
		}
			  
	    }
	}	 
#endif
      /*
       * d( div(v) )/dmesh
       */
      /*
       * There is currently no need for div(pv), so this isn't cloned.
       */
      siz = sizeof(double)*DIM*MDE;
      memset(fv->d_div_v_dmesh, 0, siz);

      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      for (p = 0; p < VIM; p++)
		{
		  fv->d_div_v_dmesh[b][j] += fv->d_grad_v_dmesh[p][p] [b][j];
		}
	    }
	}
    } else if (upd->vp[VELOCITY1] != -1 && okToZero) 
      {
	siz = sizeof(double)*DIM*DIM*DIM*MDE;
	memset(fv->d_grad_v_dmesh, 0, siz);
	siz = sizeof(double)*DIM*MDE;
	memset(fv->d_div_v_dmesh, 0, siz);
      }


  /*
   * d( grad(pv))/dmesh
   */      
  if ( pd->v[PVELOCITY1] )
    {
      v = PVELOCITY1;
      bfv = bf[v];
      vdofs     = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_pv_dmesh, 0, siz);	  
	  
      for (r = 0; r < wim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {
	      v_ri = *esp->pv[r][i];
	      for (p = 0; p < VIM; p++)
		{
		  for (q = 0; q < VIM; q++)
		    {
		      for (b = 0; b < dim; b++)
			{
			  for ( j=0; j<mdofs; j++)
			    {
			      fv->d_grad_pv_dmesh[p][q] [b][j] +=
  				v_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];  
			    }
			}
		    }
		}
	    }
	}
    } else   if ( upd->vp[PVELOCITY1] != -1  && okToZero ) {
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_pv_dmesh, 0, siz);
    }

  /*
   * d( grad(ext_v))/dmesh
   */
  if ( pd->v[EXT_VELOCITY] )
    {
      v = EXT_VELOCITY;
      bfv = bf[v];
      vdofs     = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_ext_v_dmesh, 0, siz);
	  
      for ( i=0; i<vdofs; i++)
	{
	  v_ri = *esp->ext_v[i]; 
	  for (p = 0; p < VIM; p++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      fv->d_grad_ext_v_dmesh[p][b][j] +=
			v_ri * bfv->d_grad_phi_dmesh[i][p] [b][j];  
		    }
		}
	    }
	}
    } else if ( upd->vp[EXT_VELOCITY] != -1 && okToZero ) {
      siz = sizeof(double) * DIM * DIM * MDE;
      memset(fv->d_grad_ext_v_dmesh, 0, siz);
    }
  
  /*
   * d(grad(n))/dmesh  -- ls normal vector wrt mesh
   *
   * This is the gradient of the surface normal with respect to the mesh displacements.
   * This is a tensor quantity.
   *     d_grad_n_dmesh[p][q] [b][j]  -
   *
   *           p is the first coordinate of the tensor (nominally the derivative)
   *           q is the second cordinate of the tensor (nominally the direction of n)
   *           b is the mesh displacement coordinate.
   *           j is the mesh displacement degree of freedom in the current element
   *             (note for SHELL_NORMAL_1 this element refers to the shell element dof)
   *            
   */
  if (pd->v[NORMAL1] || pd->v[SHELL_NORMAL1])
    {
      if (pd->v[NORMAL1]) v = NORMAL1;
      if (pd->v[SHELL_NORMAL1]) v = SHELL_NORMAL1;
#ifdef DEBUG_HKM
      if (ei->ielem == 165) {
	//printf("we are here 165\n");
      }
#endif	  
      bfv = bf[v];
      vdofs     = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_n_dmesh, 0, siz);
      siz = sizeof(double)*DIM*MDE;
      memset(fv->d_div_n_dmesh, 0, siz);
	  
      // HKM -> looks good according to d_grad_v_dmesh[][][][]
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      for (b = 0; b < dim; b++)
	        {
		  for (j = 0; j < mdofs; j++)
		    {
		      for (r = 0; r < dim; r++)
			{
			  for (i = 0; i < vdofs; i++)
			    {
			      v_ri = *esp->n[r][i];				  
			      fv->d_grad_n_dmesh[p][q] [b][j] += v_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];
			    }
			}
		    }
		}
	    }
	}
	  
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      for (p = 0; p < VIM; p++)
		{
		  fv->d_div_n_dmesh[b][j] += fv->d_grad_n_dmesh[p][p] [b][j];
		}
	    }
	}
    }
  else if (((upd->vp[NORMAL1] != -1) || (upd->vp[SHELL_NORMAL1] != -1)) && (okToZero))
    {
      siz = sizeof(double)*DIM*VIM*DIM*MDE;
      memset(fv->d_grad_n_dmesh,0, siz);
      siz = sizeof(double)*DIM*MDE;
      memset(fv->d_div_n_dmesh,0, siz);
    }


  /*
   * d( grad(E_field))/dmesh
   */

  if (pd->v[EFIELD1])
    {
      v = EFIELD1;
      bfv = bf[v];
      vdofs     = ei->dof[v];
	  
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_E_field_dmesh,0, siz);
	  
      for (r = 0; r < wim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {
	      v_ri = *esp->E_field[r][i];
	      for (p = 0; p < VIM; p++)
		{
		  for (q = 0; q < VIM; q++)
		    {
		      for (b = 0; b < dim; b++)
			{
			  for (j = 0; j < mdofs; j++)
			    {
			      fv->d_grad_E_field_dmesh[p][q] [b][j] +=
  				v_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];  
			    }
			}
		    }
		}
	    }
	}
    } else if (upd->vp[EFIELD1] != -1 && okToZero )
      {
	siz = sizeof(double)*DIM*VIM*DIM*MDE;
	memset(fv->d_grad_E_field_dmesh,0, siz);
      }
      
  /*
   * d( grad(d) )/dmesh
   * d( grad(d_dot)/dmesh
   */

  /* Yes, I know it should be each "+p", */
  /* but this way we can collect three */
  /* components of vectors and tensors for */
  /* less than 3 independent spatial dimensions*/



  if (pd->v[MESH_DISPLACEMENT1] == 1)
    {
      v = MESH_DISPLACEMENT1;
      bfv = bf[v];
      vdofs     = ei->dof[v];
	  
#ifdef DO_NO_UNROLL
      siz = sizeof(double)*DIM*MDE;
      memset(fv->d_div_d_dmesh,0,siz);
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_d_dmesh,0,siz);
      memset(fv->d_grad_d_dot_dmesh,0,siz);
      for (r = 0; r < dim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {
	      d_ri = *esp->d[r][i]; 
	      for (p = 0; p < VIM; p++)
		{
		  for (q = 0; q < VIM; q++)
		    {
		      for (b = 0; b < dim; b++)
			{
			  for (j = 0; j < mdofs; j++)
			    {
			      /*  fv->d_grad_d_dmesh[p][q] [b][j] = 0.; */
			      fv->d_grad_d_dmesh[p][q] [b][j] +=
				d_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j]; 
							  
			      fv->d_grad_d_dot_dmesh[p][q] [b][j] +=
				*esp_dot->d[r][i] * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j]; 
			    }
			}
		    }
		}
	    }
	}
#else

      d_ri = *esp->d[0][0]; 
      d_dot_ri = *esp_dot->d[0][0];
	  
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      /* The first element of the (r,i) matrix is split out so that it can be initialized by direct assignment
	       * This saves having to keep doing memsets on d_grad_d_mesh */
			   
			   
	      /*	  fv->d_grad_d_dmesh[p][q] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j]; */
			  
	      fv->d_grad_d_dmesh[0][0] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][0] [b][j]; 
	      fv->d_grad_d_dmesh[1][1] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][1] [b][j]; 
	      fv->d_grad_d_dmesh[1][0] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][0] [b][j]; 
	      fv->d_grad_d_dmesh[0][1] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][1] [b][j]; 
			  
			  
			  
	      /*  fv->d_grad_d_dot_dmesh[p][q] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];*/
			  
			  
	      fv->d_grad_d_dot_dmesh[0][0] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][0] [b][j]; 
	      fv->d_grad_d_dot_dmesh[1][1] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][1] [b][j]; 
	      fv->d_grad_d_dot_dmesh[1][0] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][0] [b][j]; 
	      fv->d_grad_d_dot_dmesh[0][1] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][1] [b][j]; 
			  
	      if( VIMis3 ) {
		fv->d_grad_d_dmesh[2][2] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][2] [b][j]; 
		fv->d_grad_d_dmesh[2][0] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][0] [b][j]; 
		fv->d_grad_d_dmesh[2][1] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][1] [b][j]; 
		fv->d_grad_d_dmesh[0][2] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][2] [b][j]; 
		fv->d_grad_d_dmesh[1][2] [b][j] = d_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][2] [b][j]; 
				  
				  
		fv->d_grad_d_dot_dmesh[2][2] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][2] [b][j];
		fv->d_grad_d_dot_dmesh[2][0] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][0] [b][j];
		fv->d_grad_d_dot_dmesh[2][1] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [2][1] [b][j];
		fv->d_grad_d_dot_dmesh[0][2] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [0][2] [b][j];
		fv->d_grad_d_dot_dmesh[1][2] [b][j] = d_dot_ri * bfv->d_grad_phi_e_dmesh[0][0] [1][2] [b][j];
	      } 
	    }
	}
	  
      for (r = 0; r < dim; r++)
	{
	  for (i = 0; i < vdofs; i++)
	    {
	      d_ri = *esp->d[r][i]; 
	      d_dot_ri = *esp_dot->d[r][i];
			  
	      for (b = 0; (r+i) && (b < dim); b++)   /* the (r+i) test is needed to make sure only the (r=0,i=0) element is excluded */
		{
		  for ( j=0; j<mdofs; j++)
		    {

		      /*	  fv->d_grad_d_dmesh[p][q] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j]; */
					  
		      fv->d_grad_d_dmesh[0][0] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][0] [b][j]; 
		      fv->d_grad_d_dmesh[1][1] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][1] [b][j]; 
		      fv->d_grad_d_dmesh[1][0] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][0] [b][j]; 
		      fv->d_grad_d_dmesh[0][1] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][1] [b][j]; 
					  
					  
					  
		      /*  fv->d_grad_d_dot_dmesh[p][q] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];*/
					  
					  
		      fv->d_grad_d_dot_dmesh[0][0] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][0] [b][j]; 
		      fv->d_grad_d_dot_dmesh[1][1] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][1] [b][j]; 
		      fv->d_grad_d_dot_dmesh[1][0] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][0] [b][j]; 
		      fv->d_grad_d_dot_dmesh[0][1] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][1] [b][j]; 
					  
		      if( VIMis3 ) {
			fv->d_grad_d_dmesh[2][2] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][2] [b][j]; 
			fv->d_grad_d_dmesh[2][0] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][0] [b][j]; 
			fv->d_grad_d_dmesh[2][1] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][1] [b][j]; 
			fv->d_grad_d_dmesh[0][2] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][2] [b][j]; 
			fv->d_grad_d_dmesh[1][2] [b][j] += d_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][2] [b][j]; 
						  
						  
			fv->d_grad_d_dot_dmesh[2][2] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][2] [b][j];
			fv->d_grad_d_dot_dmesh[2][0] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][0] [b][j];
			fv->d_grad_d_dot_dmesh[2][1] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [2][1] [b][j];
			fv->d_grad_d_dot_dmesh[0][2] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [0][2] [b][j];
			fv->d_grad_d_dot_dmesh[1][2] [b][j] += d_dot_ri * bfv->d_grad_phi_e_dmesh[i][r] [1][2] [b][j];
		      }
						  

		    }
		}
	    }
	}
#endif	  
	  
#ifdef DO_NO_UNROLL
      for (p = 0; p < VIM; p++)
	{
          for (q = 0; q < VIM; q++)
	    {
              for (b = 0; b < dim; b++)
		{
                  for ( j=0; j<mdofs; j++)
		    {
                      /*
                       * Special added Lagrangian piece for d grad_d / d mesh...
                       */
                      fv->d_grad_d_dmesh[p][q] [b][j] += bfv->grad_phi_e[j][b] [p][q];
					  
		      /*N.B. PRS: need to actually furbish this for Newmark Beta schemes */
		      if (transient_run && discontinuous_porous_media ) 
			{
			  fv->d_grad_d_dot_dmesh[p][q] [b][j] += (1.+2.*tran->theta)*
			    bfv->grad_phi_e[j][b] [p][q]/tran->delta_t;
			}
		    }
		}
	    }
	}
#else
	  
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      /*
	       * Special added Lagrangian piece for d grad_d / d mesh...
	       */
	      /*  fv->d_grad_d_dmesh[p][q] [b][j] += bfv->grad_phi_e[j][b] [p][q];*/
			  
	      fv->d_grad_d_dmesh[0][0] [b][j] += bfv->grad_phi_e[j][b] [0][0];
	      fv->d_grad_d_dmesh[1][1] [b][j] += bfv->grad_phi_e[j][b] [1][1];
	      fv->d_grad_d_dmesh[1][0] [b][j] += bfv->grad_phi_e[j][b] [1][0];
	      fv->d_grad_d_dmesh[0][1] [b][j] += bfv->grad_phi_e[j][b] [0][1];
			  
	      if ( VIMis3 ) {
		fv->d_grad_d_dmesh[2][2] [b][j] += bfv->grad_phi_e[j][b] [2][2];
		fv->d_grad_d_dmesh[2][0] [b][j] += bfv->grad_phi_e[j][b] [2][0];
		fv->d_grad_d_dmesh[2][1] [b][j] += bfv->grad_phi_e[j][b] [2][1];
		fv->d_grad_d_dmesh[0][2] [b][j] += bfv->grad_phi_e[j][b] [0][2];
		fv->d_grad_d_dmesh[1][2] [b][j] += bfv->grad_phi_e[j][b] [1][2];
	      }

	      /*N.B. PRS: need to actually furbish this for Newmark Beta schemes */
	      if (transient_run && discontinuous_porous_media ) 
		{
		  /* fv->d_grad_d_dot_dmesh[p][q] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [p][q]/tran->delta_t;*/
				 
		  fv->d_grad_d_dot_dmesh[0][0] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [0][0]/tran->delta_t;
		  fv->d_grad_d_dot_dmesh[1][1] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [1][1]/tran->delta_t;
		  fv->d_grad_d_dot_dmesh[1][0] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [1][0]/tran->delta_t;
		  fv->d_grad_d_dot_dmesh[0][1] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [0][1]/tran->delta_t;
				  
		  if ( VIMis3 ) {
				  
		    fv->d_grad_d_dot_dmesh[2][2] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [2][2]/tran->delta_t;
		    fv->d_grad_d_dot_dmesh[2][0] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [2][0]/tran->delta_t;
		    fv->d_grad_d_dot_dmesh[2][1] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [2][1]/tran->delta_t;
		    fv->d_grad_d_dot_dmesh[0][2] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [0][2]/tran->delta_t;
		    fv->d_grad_d_dot_dmesh[1][2] [b][j] += (1.+2.*tran->theta)*bfv->grad_phi_e[j][b] [1][2]/tran->delta_t;
		  }	  
		}
	    }
	}
	  
#endif
	  
      
      /*
       * d( div(d) )/dmesh [ just trace(grad(d)) ]
       */
#ifdef DO_NO_UNROLL
      for (b = 0; b < dim; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      for (p = 0; p < VIM; p++)
		{
		  fv->d_div_d_dmesh[b][j] += fv->d_grad_d_dmesh[p][p] [b][j];
		}
	    }
	}
#else
	  
	  
      fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[0][0] [0][0];
      fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[1][1] [1][0];
      fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[1][1] [0][0];
      fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[0][0] [1][0];
	  
      if( VIMis3 ) {
	fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[2][2] [2][0];
	fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[0][0] [2][0];
	fv->d_div_d_dmesh[2][0] = fv->d_grad_d_dmesh[1][1] [2][0];
	fv->d_div_d_dmesh[0][0] = fv->d_grad_d_dmesh[2][2] [0][0];
	fv->d_div_d_dmesh[1][0] = fv->d_grad_d_dmesh[2][2] [1][0];
      }			  
	  
      for (j = 1; j < mdofs; j++)
	{
		  
	  /* fv->d_div_d_dmesh[b][j] += fv->d_grad_d_dmesh[p][p] [b][j]; */

	  fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[0][0] [0][j];
	  fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[1][1] [1][j];
	  fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[1][1] [0][j];
	  fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[0][0] [1][j];
		  
	  if( VIMis3 ) {
	    fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[2][2] [2][j];
	    fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[0][0] [2][j];
	    fv->d_div_d_dmesh[2][j] += fv->d_grad_d_dmesh[1][1] [2][j];
	    fv->d_div_d_dmesh[0][j] += fv->d_grad_d_dmesh[2][2] [0][j];
	    fv->d_div_d_dmesh[1][j] += fv->d_grad_d_dmesh[2][2] [1][j];
	  }			  

		  
	}
	  
	  
#endif
	  
    } else   if ( upd->vp[MESH_DISPLACEMENT1] != -1 && okToZero )
      {
	siz = sizeof(double)*DIM*DIM*DIM*MDE;
	memset(fv->d_grad_d_dmesh,0,siz);
	memset(fv->d_grad_d_dot_dmesh,0,siz);
	siz = sizeof(double)*DIM*MDE;
	memset(fv->d_div_d_dmesh,0,siz);
      }

  /*
   * d( grad(d_rs) )/dmesh and d( grad(d_rs) ) /dd_rs
   */

  /* Yes, I know it should be each "+p", */
  /* but this way we can collect three */
  /* components of vectors and tensors for */
  /* less than 3 independent spatial dimensions*/



  if ( pd->v[SOLID_DISPLACEMENT1] )
    {
      v = SOLID_DISPLACEMENT1;
      bfv = bf[v];
      vdofs     = ei->dof[v];
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_d_rs_dmesh,0,siz);
      siz = sizeof(double)*DIM*MDE;
      memset(fv->d_div_d_rs_dmesh,0,siz);
	  
      for (p = 0; p < VIM; p++)
	{
	  for (q = 0; q < VIM; q++)
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      /*  fv->d_grad_d_rs_dmesh[p][q] [b][j] = 0.; */
		      for (r = 0; r < dim; r++)
			{
			  for (i = 0; i < vdofs; i++)
			    {
			      fv->d_grad_d_rs_dmesh[p][q] [b][j] +=
				*esp->d_rs[r][i] 
				*  bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];
			    }
			}
		    }
		}
	    }
	}
      
      /*
       * d( div(d_rs) )/dmesh [ just trace(grad(d_rs)) ]
       */
      for (b = 0; b < VIM; b++)
	{
	  for (j = 0; j < mdofs; j++)
	    {
	      for (p = 0; p < VIM; p++)
		{
		  fv->d_div_d_rs_dmesh[b][j] += fv->d_grad_d_rs_dmesh[p][p] [b][j];
		}
	    }
	}
    } else   if ( upd->vp[SOLID_DISPLACEMENT1] != -1  && okToZero)
      {
	siz = sizeof(double)*DIM*DIM*DIM*MDE;
	memset(fv->d_grad_d_rs_dmesh,0,siz);
	siz = sizeof(double)*DIM*MDE;
	memset(fv->d_div_d_rs_dmesh,0,siz);
      }

  /*
   * d(grad(S))/dmesh
   */


      
  if ( pd->v[POLYMER_STRESS11] )
    {
      v = POLYMER_STRESS11;
      vdofs  = ei->dof[v];
      bfv = bf[v];
      siz = sizeof(double)*MAX_MODES*DIM*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_S_dmesh,0,siz);
      siz = sizeof(double)*MAX_MODES*DIM*DIM*MDE;
      memset(fv->d_div_S_dmesh,0,siz);	  
	  
      for ( mode=0; mode<vn->modes; mode++)
	{
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  for ( r=0; r<VIM; r++)
		    {
		      for (b = 0; b < VIM; b++)
			{
			  for ( j=0; j<mdofs; j++)
			    {
			      for ( i=0; i<vdofs; i++)
				{
				  if( p <= q)
				    {
				      fv->d_grad_S_dmesh[mode][r][p][q] [b][j] +=
					*esp->S[mode][p][q][i]  *  bfv->d_grad_phi_dmesh[i][r] [b][j];
				    }
				  else
				    {
				      fv->d_grad_S_dmesh[mode][r][p][q] [b][j] +=
					*esp->S[mode][q][p][i]  *  bfv->d_grad_phi_dmesh[i][r] [b][j];
				    }
				}
			    }
			}
		    }
		}
	    }
		  
	  for ( r=0; r<VIM; r++)
	    { 
	      for ( q=0; q<VIM; q++)
		{
		  for ( b=0; b<VIM; b++)
		    {
		      for ( j=0; j<mdofs; j++)
			{
			  fv->d_div_S_dmesh[mode][r] [b][j] +=
			    fv->d_grad_S_dmesh[mode][q][q][r] [b][j];
			}
		    }
		}
	    }
		  
	  if ( pd->CoordinateSystem != CARTESIAN )
	    {
	      for ( s=0; s<VIM; s++)
		{
		  for ( r=0; r<VIM; r++)
		    {
		      for ( p=0; p<VIM; p++)
			{
			  for ( q=0; q<VIM; q++)
			    { 
			      for (b = 0; b < VIM; b++)
				{
				  for ( j=0; j<mdofs; j++)
				    {
				      fv->d_div_S_dmesh[mode][s] [b][j] += 
					fv->S[mode][p][q]*
					( fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s,q) +
					  fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r,p) );
				    }
				}
			    }
			}
		    }
		}
	    }
	} /* end of modal loop */
    } else if ( upd->vp[POLYMER_STRESS11] != -1  && okToZero) {
      siz = sizeof(double)*MAX_MODES*DIM*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_S_dmesh,0,siz);
      siz = sizeof(double)*MAX_MODES*DIM*DIM*MDE;
      memset(fv->d_div_S_dmesh,0,siz);
    }
  
  
  
  /*
   * d(grad(G))/dmesh
   */


      
  if ( pd->v[VELOCITY_GRADIENT11] )
    {
      v = VELOCITY_GRADIENT11;
      bfv = bf[v];
      vdofs  = ei->dof[v];
	  
      siz = sizeof(double)*DIM*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_G_dmesh,0,siz);
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_div_G_dmesh,0,siz);
      memset(fv->d_div_Gt_dmesh,0,siz);
	  
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      for ( r=0; r<VIM; r++)
		{
		  for ( b=0; b<VIM; b++)
		    {
		      for ( j=0; j<mdofs; j++)
			{
			  for ( i=0; i<vdofs; i++)
			    {
			      fv->d_grad_G_dmesh[r][p][q] [b][j] +=
				*esp->G[p][q][i]  *  bfv->d_grad_phi_dmesh[i][r] [b][j];
			    }
			}
		    }
		}
	    }
	}
      
      for ( r=0; r<VIM; r++)
	{ 
	  for ( q=0; q<VIM; q++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      fv->d_div_G_dmesh[r] [b][j] +=
			fv->d_grad_G_dmesh[q][q][r] [b][j];
		    }
		}
	    }
	}
      
      if ( pd->CoordinateSystem != CARTESIAN )
	{
	  for ( s=0; s<VIM; s++)
	    {
	      for ( r=0; r<VIM; r++)
		{
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{ 
			  for ( b=0; b<dim; b++)
			    {
			      for ( j=0; j<mdofs; j++)
				{
				  fv->d_div_G_dmesh[s] [b][j] += 
				    fv->G[p][q]*
				    ( fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s,q) +
				      fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r,p) );
				}
			    }
			}
		    }
		}
	    }
	}

      for ( r=0; r<VIM; r++)
	{ 
	  for ( q=0; q<VIM; q++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  for ( j=0; j<mdofs; j++)
		    {
		      fv->d_div_Gt_dmesh[r] [b][j] +=
			fv->d_grad_G_dmesh[q][r][q] [b][j];
		    }
		}
	    }
	}
      
      if ( pd->CoordinateSystem != CARTESIAN )
	{
	  for ( s=0; s<VIM; s++)
	    {
	      for ( r=0; r<VIM; r++)
		{
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{ 
			  for ( b=0; b<VIM; b++)
			    {
			      for ( j=0; j<mdofs; j++)
				{
				  fv->d_div_G_dmesh[s] [b][j] += 
				    fv->G[q][p]*
				    ( fv->d_grad_e_dq[p][r][q][b] * bfm->phi[j] * (double)delta(s,q) +
				      fv->d_grad_e_dq[q][p][s][b] * bfm->phi[j] * (double)delta(r,p) );
				}
			    }
			}
		    }
		}
	    }
	}

    } else   if ( upd->vp[VELOCITY_GRADIENT11] != -1  && okToZero)  {
      siz = sizeof(double)*DIM*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_G_dmesh,0,siz);
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_div_G_dmesh,0,siz);
      memset(fv->d_div_Gt_dmesh,0,siz);
    }

  /*
   * d(grad(vd))/dmesh
   */
  if ( pd->v[VORT_DIR1] )
    {
      v = VORT_DIR1;
      bfv = bf[v];
      vdofs = ei->dof[v];
	  
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_vd_dmesh,0, siz);
	  
      for (r = 0; r < DIM; r++)
	{
	  for(i = 0; i < vdofs; i++)
	    {
	      v_ri = *esp->vd[r][i];
	      for (p = 0; p < DIM; p++)
		{
		  for(q = 0; q < DIM; q++)
		    {
		      for(b = 0; b < VIM; b++)
			{
			  for(j = 0; j < mdofs; j++)
			    {
			      fv->d_grad_vd_dmesh[p][q] [b][j] +=
				v_ri * bfv->d_grad_phi_e_dmesh[i][r] [p][q] [b][j];  
			    }
			}
		    }
		}
	    }
	}
      
      /*
       * d( div(vd) )/dmesh
       */
	  
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      for ( p=0; p<DIM; p++)
		{
		  fv->d_div_vd_dmesh[b][j] += fv->d_grad_vd_dmesh[p][p] [b][j];
		}
	    }
	}
    } else if ( upd->vp[VORT_DIR1] != -1  && okToZero) {
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_grad_vd_dmesh,0, siz);
    }

  
  
  /*
   *   d(n_dot_curl_s_v)/dmesh
   *        This is carried out on a shell
   */
  if (pd->v[N_DOT_CURL_V])
    {
      v = N_DOT_CURL_V;
      bfv = bf[v];
      vdofs     = ei->dof[v];	  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_n_dot_curl_s_v_dmesh,0, siz);
      for (p = 0; p < VIM; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_n_dot_curl_s_v_dmesh[p] [b][j] += (*esp->n_dot_curl_s_v[i] * bfv->d_grad_phi_dmesh[i][p] [b][j]);
		    }
		}
	    }
	}
    } 
  else if (upd->vp[N_DOT_CURL_V] != -1 && okToZero) 
    {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_n_dot_curl_s_v_dmesh, 0, siz);
    }

  /*
   *   d(grad_div_s_v_dmesh[b][jvar][jShell]);
   *        This is carried out on a shell
   */
  if (pd->v[SHELL_SURF_DIV_V])
    {
      v = SHELL_SURF_DIV_V;
      bfv = bf[v];
      vdofs     = ei->dof[v];	  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_div_s_v_dmesh,0, siz);
      for (p = 0; p < VIM; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_div_s_v_dmesh[p] [b][j] += (*esp->div_s_v[i] * bfv->d_grad_phi_dmesh[i][p] [b][j]);
		    }
		}
	    }
	}
    }
  else if (upd->vp[SHELL_SURF_DIV_V] != -1 && okToZero) 
    {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_div_s_v_dmesh, 0, siz);
    }

  /*
   *   d_grad_curv_dmesh[b][jvar][jShell])
   *        This is carried out on a shell
   */
  if (pd->v[SHELL_SURF_CURV])
    {
      v = SHELL_SURF_CURV;
      bfv = bf[v];
      vdofs = ei->dof[v];	  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_curv_dmesh,0, siz);
      for (p = 0; p < VIM; p++ )
	{
	  for (b = 0; b < dim; b++)
	    {
	      for (j = 0; j < mdofs; j++)
		{
		  for (i = 0; i < vdofs; i++)
		    {
		      fv->d_grad_curv_dmesh[p] [b][j] += (*esp->curv[i] * bfv->d_grad_phi_dmesh[i][p] [b][j]);
		    }
		}
	    }
	}
    }
  else if (upd->vp[SHELL_SURF_CURV] != -1 && okToZero) 
    {
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(fv->d_grad_curv_dmesh, 0, siz);
    }

  /*   
   *   d_serialgrad_div_s_v_dmesh[b][jvar][jShell]);
   *        This is carried out on a shell
   */
  if (pd->v[GRAD_S_V_DOT_N1])
    {
      for (r = 0; r < dim; r++) 
	{
	  v = GRAD_S_V_DOT_N1 + r;
	  bfv = bf[v];
	  vdofs = ei->dof[v];	  
	  siz = sizeof(double)*DIM*DIM*DIM*MDE;
	  memset(fv->d_serialgrad_grad_s_v_dot_n_dmesh,0, siz);
	  for (p = 0; p < VIM; p++ )
	    {
	      for (b = 0; b < dim; b++)
		{
		  for (j = 0; j < mdofs; j++)
		    {
		      for (i = 0; i < vdofs; i++)
			{
			  fv->d_serialgrad_grad_s_v_dot_n_dmesh[r][p][b][j] += (*esp->grad_v_dot_n[r][i] * bfv->d_grad_phi_dmesh[i][p] [b][j]);
			}
		    }
		}
	    }
	}
    }
  else if (upd->vp[GRAD_S_V_DOT_N1] != -1 && okToZero) 
    {
      siz = sizeof(double)*DIM*DIM*DIM*MDE;
      memset(fv->d_serialgrad_grad_s_v_dot_n_dmesh, 0, siz);
    }


  /*
   * d(grad(sh_J))/dmesh
   */                                                                               
  if (pd->v[SHELL_DIFF_FLUX])
    {
      v = SHELL_DIFF_FLUX;
      bfv = bf[v];
      vdofs  = ei->dof[v];
	  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(&(fv->d_grad_sh_J_dmesh[0][0][0]),0, siz);
	  
      for (p = 0; p < VIM; p++)
        {
          for (b = 0; b < dim; b++)
            {
              for ( j=0; j<mdofs; j++)
                {
                  for ( i=0; i<vdofs; i++)
                    {
                      fv->d_grad_sh_J_dmesh[p] [b][j] +=
                        *esp->sh_J[i]  *  bfv->d_grad_phi_dmesh[i][p] [b][j];
                    }
                }
            }
        }
    } else if ( upd->vp[SHELL_DIFF_FLUX] == -1 && okToZero )
      {
	siz = sizeof(double)*DIM*DIM*MDE;
	memset(&(fv->d_grad_sh_J_dmesh[0][0][0]),0, siz);
      }

  /*
   * Now load up the sensitivity of the volume element to the
   * mesh variables...
   *
   * d ( h3 ) / d ( mesh_bj )
   */
  if (pd->v[MESH_DISPLACEMENT1] == 1)
    {
      for ( b=0; b<VIM; b++)
	{
	  for ( j=0; j<mdofs; j++)
	    {
	      fv->dh3dmesh[b][j] = fv->dh3dq[b] * bfm->phi[j];
	    }
	}
    }

  return(status);
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

double 
density(DENSITY_DEPENDENCE_STRUCT *d_rho, double time)
     
     /**************************************************************************
      *
      * density
      *
      *   Calculate the density and its derivatives at the local gauss point
      *
      * Output
      * -----
      *    dbl d_rho->T[j]    -> derivative of density wrt the jth
      *                          Temperature unknown in an element
      *    dbl d_rho->C[k][j] -> derivative of density wrt the jth
      *                          species unknown of ktype, k, in the element.
      *    dbl d_rho->F[j]    -> derivative of density wrt the jth
      *                          FILL unknown in an element
      *
      *  Return
      * --------
      *    Actual value of the density
      ***************************************************************************/
{
  int w, j, var, var_offset, matrl_species_var_type,
    dropped_last_species_eqn;
  int species, err;
  dbl vol=0, rho=0, rho_f, rho_s, pressureThermo, RGAS_CONST;
  dbl avgMolecWeight=0, tmp;
  double *phi_ptr;

  dbl sv[MAX_CONC];
  dbl sv_p, sum_sv=0;
  struct Level_Set_Data *ls_old;

  /*
   * Get the material's species variable type
   *    -> We will assume here that all species have the same var type.
   */
  matrl_species_var_type = mp->Species_Var_Type;
  
  /*
   * Decide whether the continuity equation for the last species
   * in the species list has been dropped in favor of an implicit
   * sum of mass fractions equals one constraint. If it has, then
   * that implicit constraint has to be reflected in the Jacobian terms
   * calculated from this routine. Basically, any change in the
   * mass fraction or concentration of one species is compensated by a
   * change in the mass fraction or concentration of the last species
   * in the mechanism.
   */
  dropped_last_species_eqn = mp->Dropped_Last_Species_Eqn;
  
  /* Zero out sensitivities */
  if (d_rho != NULL) {
    zeroStructures(d_rho, 1);
  }
  
  /*
   *   Branch according to the density model from the material properties
   *   structure
   */
  if (mp->DensityModel == CONSTANT) {
    rho   = mp->density;
  } else if (mp->DensityModel == USER)  {
    (void) usr_density(mp->u_density);      
    rho   = mp->density;

    if ( d_rho != NULL )
      {
	var = TEMPERATURE;
	for ( j=0; j<ei->dof[var]; j++)
	  {
	    d_rho->T[j]= mp->d_density[var]*bf[var]->phi[j];
	  }

	if (pd->v[MASS_FRACTION] )
	  {
	    for ( w=0; w<pd->Num_Species_Eqn; w++)
	      {
		var = MASS_FRACTION;
		var_offset = MAX_VARIABLE_TYPES + w;
		for ( j=0; j<ei->dof[var]; j++)
		  {
		    d_rho->C[w][j] =mp->d_density[var_offset]*bf[var]->phi[j];
		  }
	      }
	  }
      }
    
  } else if (mp->DensityModel == FILL ) {
      
    rho = mp->u_density[0];

    if(fv->F < 0.1)  {
      rho = mp->u_density[1];
    }

  } else if (mp->DensityModel == DENSITY_LEVEL_SET ) 
    {
      double rho0  = mp->u_density[0];
      double rho1  = mp->u_density[1];
      double width = mp->u_density[2];
      if ( d_rho == NULL )
        err = level_set_property(rho0, rho1, width, &rho, NULL);
      else
        err = level_set_property(rho0, rho1, width, &rho, d_rho->F);
      EH(err, "level_set_property() failed for density.");
    }

  else if (mp->DensityModel == DENSITY_CONST_PHASE_FUNC ) 
    {
      int num, a;
      double width;
      double rho1, rho2 = 0, tmp_rho;
      struct Level_Set_Data *ls_save = ls;
	
      num=pfd->num_phase_funcs;
      width = mp->u_density[num];
      //Major cludgey here.  this is fubar'd like viscosity    
		
      for( a = 0; a < num; a++)
	{
	  rho1 = mp->u_density[a];
	  rho1 = mp->u_density[a+2];
	  ls = pfd->ls[a];
			
	  if( d_rho != NULL)
	    err = level_set_property( rho1, rho2, width, &tmp_rho, d_rho->pf[a] );
	  else
	    err = level_set_property( rho1, rho2, width, &tmp_rho, NULL);
			
	  rho += tmp_rho;
	}
      if ( fabs(rho) < DBL_SMALL)  rho = mp->u_density[num + 1 ];
      ls = ls_save;

    }

  else if (mp->DensityModel == DENSITY_FOAM)
    {
      double x0, Rgas, MW, rho_epoxy, rho_fluor, T, vol, Press;
      species = (int) mp->u_density[0];  /* species number fluorinert */
      x0 = mp->u_density[1];         /* Initial fluorinert mass fraction */
      Rgas =  mp->u_density[2];      /* Gas constant in appropriate units */
      MW = mp->u_density[3];
      rho_epoxy = mp->u_density[4];  /* Density of epoxy resin */
      rho_fluor =  mp->u_density[5]; /* Density of liquid fluorinert */
      
      if(fv->c[species] > 0.)
	vol = fv->c[species];
      else
	vol = 0.;

      if(vol > x0) vol = x0;

      
      T = fv->T;
      Press = upd->Pressure_Datum;
      rho = 1./((x0-vol)*Rgas*T/(Press*MW)+(1.0-x0)/rho_epoxy+vol/rho_fluor);
      
      var = MASS_FRACTION;
      if (vol > 0. && d_rho != NULL )
	{
	  if (pd->v[var] )
	    {
	      if((vol > 0.)&& (vol < x0) )
		for ( j=0; j<ei->dof[var]; j++)
		  {
		    d_rho->C[species][j] = rho*rho*(Rgas*T/(Press*MW)-1./rho_fluor)*bf[var]->phi[j];
		  }
	    }
	}

      var = TEMPERATURE;
      if(d_rho != NULL )
	{
	  if (pd->v[var] )
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_rho->T[j] = -rho*rho*((x0-vol)*Rgas/(Press*MW))*bf[var]->phi[j];
		}
	    }
	}
      
    }
  else if (mp->DensityModel == DENSITY_FOAM_CONC)
    {
      double Rgas, MW_f, MW_a, rho_epoxy, rho_fluor, T, Press;
      dbl rho_v_inv,d_rho_v_inv_dT, rho_a_inv, d_rho_a_inv_dT, drho_T, drho_c_v,
	drho_c_a, drho_c_l;
      int species_l, species_v, species_a;

      species_l = (int) mp->u_density[0];  /* species number fluorinert liquid */
      species_v = (int) mp->u_density[1];  /* species number fluorinert vapor */
      species_a = (int) mp->u_density[2];  /* species number air vapor */
      Rgas =  mp->u_density[3];      /* Gas constant in appropriate units */
      MW_f = mp->u_density[4];         /* molecular weight fluorinert */
      MW_a = mp->u_density[5];         /* molecular weight air */
      rho_epoxy = mp->u_density[6];  /* Density of epoxy resin */
      rho_fluor =  mp->u_density[7]; /* Density of liquid fluorinert */
            
      T = fv->T;
      Press = upd->Pressure_Datum;
      rho_v_inv = Rgas*T/(Press*MW_f);
      d_rho_v_inv_dT = Rgas/(Press*MW_f);
      rho_a_inv = Rgas*T/(Press*MW_a);
      d_rho_a_inv_dT = Rgas/(Press*MW_a);
	
      drho_c_v =1.-rho_epoxy*rho_v_inv;
      drho_c_a =1.-rho_epoxy*rho_a_inv;
      drho_c_l =1.-rho_epoxy/rho_fluor;
 
      drho_T = -fv->c[species_v]*rho_epoxy*d_rho_v_inv_dT 
	- fv->c[species_a]*rho_epoxy*d_rho_a_inv_dT ;
      
      var = MASS_FRACTION;
      if(d_rho != NULL )
	{
	  if (pd->v[var] )
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_rho->C[species_v][j] = drho_c_v * bf[var]->phi[j];
		  d_rho->C[species_a][j] = drho_c_a * bf[var]->phi[j];
		  d_rho->C[species_l][j] = drho_c_l * bf[var]->phi[j];
		}
	    }

	  var = TEMPERATURE;
	  if(d_rho != NULL )
	    {
	      if (pd->v[var] )
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_rho->T[j] = drho_T*bf[var]->phi[j];
		    }
		}
	    }
	  
	}
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME)
    {
      double rho_init, rho_final, aexp, time_delay, realtime;
      rho_init = mp->u_density[0];  /* Initial density */
      rho_final = mp->u_density[1]; /* final density */
      aexp =  mp->u_density[2];      /* Arhennius constant for time exponent*/
      time_delay =  mp->u_density[3]; /* time delay before foaming starts */
      
      if (time > time_delay)
	{
	  realtime = time - time_delay;
	  rho = rho_final+ (rho_init-rho_final)*exp(-aexp*realtime);
	}
      else
	{
	  rho = rho_init;
	}        
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
    {
      double T;
      double rho_init  = mp->u_density[0];  /* Initial density (units of gm cm-3 */
      double rho_final = mp->u_density[1];  /* Final density (units of gm cm-3 */
      double cexp      = mp->u_density[2];  /* cexp in units of Kelvin sec */
      double coffset   = mp->u_density[3];  /* offset in units of seconds */
      double time_delay= mp->u_density[4];  /* time delay in units of seconds */
      double drhoDT;
      T = fv->T;
      if (time > time_delay)
	{
	  double realtime = time - time_delay;
	  double delRho = (rho_init - rho_final);
	  double cdenom = cexp - coffset * T;
	  double expT = exp(-realtime * T / cdenom);
	  rho = rho_final + delRho * expT;
	  drhoDT = delRho * expT *(- realtime * cexp / (cdenom * cdenom));
	}
      else
	{
	  rho = rho_init;
	  drhoDT = 0.0;
	}
      var = TEMPERATURE;
      if (d_rho != NULL)
	{
	  if (pd->v[var])
	    {
	      for (j = 0; j < ei->dof[var]; j++)
		{
		  d_rho->T[j] = drhoDT * bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->DensityModel == SUSPENSION) {
    species = (int) mp->u_density[0];
    rho_f   = mp->u_density[1];
    rho_s   = mp->u_density[2];

    vol = fv->c[species];
    if(vol < 0.)
      {
	vol = 0.;
      }
      
    rho = rho_f + (rho_s - rho_f)*vol;
    var = MASS_FRACTION;
    if(vol > 0. && d_rho != NULL )
      {
	if (pd->v[var] )
	  {
	    for ( j=0; j<ei->dof[var]; j++)
	      {
		d_rho->C[species][j] = (rho_s - rho_f)*bf[var]->phi[j];
	      }
	  }
      }
  } else if (mp->DensityModel == SUSPENSION_PM) {
    /* For the SUSPENSION_PM model, the density is "usually" just the
     * fluid density.  It is only in two cases that the density we want is
     * actually the particle phase density.  1) In the particle phase 
     * momentum equations. 2) In the particle phase advection equation
     * (like a chemical species).  In those two cases, this muse be "undone".
     * This seemed the best way to do it so there was a reasonable default
     * behavior, and you didn't need a bunch of if()'s bracketing every
     * density() call...
     */

    rho = mp->u_density[1];

    /* This doesn't account for other species present when calculating
     * density.  But they are currently just "collected" into one rho_f,
     * so all of the d_rho->C's should be zero!.
     */
    /*
      var = MASS_FRACTION;
      if(vol > 0.)
      {
      if (pd->v[var] )
      {
      for ( j=0; j<ei->dof[var]; j++)
      {
      d_rho->C[species][j] = (rho_s - rho_f)*bf[var]->phi[j];
      }
      }
      }
    */
  } else if (mp->DensityModel == THERMAL_BATTERY) {
    rho  = mp->u_density[0] - mp->u_density[1]*(2.0*fv->c[0]);
    /* Ref.: Pollard & Newman 1981, p.501 */
      
    if (pd->v[MASS_FRACTION] && d_rho != NULL)
      {
	var = MASS_FRACTION;
	for ( j=0; j<ei->dof[var]; j++)
	  {
	    d_rho->C[0][j] = -mp->u_density[1]*2.0*bf[var]->phi[j];
	  }
      }
    
  } else if (mp->DensityModel == DENSITY_IDEAL_GAS) {
    /*
     * Specify the thermodynamic pressure here
     *   -> For the moment we ignore any coupling between
     *      the dynamic pressure field and the thermodynamic
     *      pressure field.
     */
    pressureThermo = upd->Pressure_Datum;

    /*
     * Specify the gas constant according to CRC IUPAC conventions
     */
    RGAS_CONST = 8.314510E7;        /*    g cm^2/(sec^2 g-mole K)  */
  
    /*
     *  Base density calculation depends on species var type
     */
    switch (matrl_species_var_type) {
    case SPECIES_MASS_FRACTION:
      avgMolecWeight = wt_from_Yk(mp->Num_Species, fv->c, mp->molecular_weight);
      rho = pressureThermo * avgMolecWeight / ( RGAS_CONST * fv->T);
      break;
    case SPECIES_MOLE_FRACTION:
      avgMolecWeight = wt_from_Xk(mp->Num_Species, fv->c, mp->molecular_weight);
      rho = pressureThermo * avgMolecWeight / ( RGAS_CONST * fv->T);
      break;
    case SPECIES_CONCENTRATION:
      rho = 0.0;
      for (w = 0; w < mp->Num_Species; w++) {
	rho += fv->c[w] * mp->molecular_weight[w];
      }


      break;
    default:
      fprintf(stderr,"Density error: species var type not handled: %d\n",
	      matrl_species_var_type);
      exit(-1);
    }

    /*
     * Now do the Jacobian terms
     *          HKM -> No dependence on the pressure is calculated
     *
     * Temperature dependence is common to all species var types
     */
    if ( d_rho != NULL && pd->v[TEMPERATURE] ) {
      phi_ptr =  bf[MASS_FRACTION]->phi;
      tmp = -rho / fv->T;
      for (j = 0; j < ei->dof[TEMPERATURE]; j++) {
	d_rho->T[j] = tmp * phi_ptr[j];
      }
    }
    /*
     * Dependence on the species unknowns.
     */
    if ( d_rho != NULL && pd->v[MASS_FRACTION] ) {
      phi_ptr =  bf[MASS_FRACTION]->phi;
      switch (matrl_species_var_type) {
      case SPECIES_MASS_FRACTION:
	if (dropped_last_species_eqn) {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp = - rho * avgMolecWeight *
	      (1.0/(mp->molecular_weight[w]) -
	       1.0/(mp->molecular_weight[mp->Num_Species_Eqn]));
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	} else {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp = - rho * avgMolecWeight / (mp->molecular_weight[w]);
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	}
	break;
      case SPECIES_MOLE_FRACTION:
	if (dropped_last_species_eqn) {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp = rho / avgMolecWeight *
	      (mp->molecular_weight[w] - mp->molecular_weight[mp->Num_Species_Eqn]);
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	} else {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp = rho / avgMolecWeight * mp->molecular_weight[w];
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	}
	break;
      case SPECIES_CONCENTRATION:
	if (dropped_last_species_eqn) {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp =  mp->molecular_weight[w] -  mp->molecular_weight[mp->Num_Species_Eqn];
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	} else {
	  for (w = 0; w < mp->Num_Species_Eqn; w++) {
	    tmp =  mp->molecular_weight[w];
	    for (j = 0; j < ei->dof[MASS_FRACTION]; j++) { 
	      d_rho->C[w][j] = tmp * phi_ptr[j];
	    }
	  }
	}
	break;
      default:
	fprintf(stderr,"Density error: species var type not handled: %d\n",
		matrl_species_var_type);
	exit(-1);
      }
    }
  }
  else if (mp->DensityModel == REACTIVE_FOAM )
    {
      /* added ACSun 04/02 */
      if (mp->SpeciesSourceModel[0] == FOAM)
	{
	  foam_species_source(mp->u_species_source[0]);
	} else {
	  EH(-1, "Must specify FOAM species source in the material's file");
	}
    
      rho = 0.;
      for (w = 0; w < pd->Num_Species_Eqn; w++)
	{
	  sv[w] = mp->specific_volume[w];
	}
      sv_p = mp->specific_volume[pd->Num_Species_Eqn];

      for (w=0; w<pd->Num_Species_Eqn; w++)
	{  
	  sum_sv += (sv[w]-sv_p)*fv->c[w];
	}
      sum_sv += sv_p;
      rho = 1/sum_sv ;
    
      var = MASS_FRACTION;
      if (d_rho != NULL && pd->v[var]) {
	phi_ptr =  bf[var]->phi;
	for (j=0; j<ei->dof[var]; j++)
	  {
	    for (w = 0; w < pd->Num_Species_Eqn; w++)
	      {
		d_rho->C[w][j]= -rho*rho*(sv[w]-sv_p)
		  *phi_ptr[j] ; 
	      }
	  }
      }
      var = TEMPERATURE;
      if (d_rho != NULL && pd->v[var]) {
	phi_ptr =  bf[var]->phi;
	for (j=0; j<ei->dof[var]; j++)
	  {
	    d_rho->T[j] = 0.;
	    for (w = 0; w < pd->Num_Species_Eqn; w++)
	      {
		d_rho->T[j] -= rho*rho*(sv[w]-sv_p)
		  *mp->species_source[w]*3.;
	      }
	    d_rho->T[j] *= phi_ptr[j] ;
	  }
      }
    }
  else if (mp->DensityModel == SOLVENT_POLYMER)
    {
      double *param = mp->u_density;

      /* added ACSun 7/99 */
      rho = 0.;
      for (w = 0; w < pd->Num_Species_Eqn; w++)
	{
	  sv[w] = mp->specific_volume[w];
	}
      sv_p = param[0];
      var = MASS_FRACTION;

      switch (matrl_species_var_type) {
      case SPECIES_MASS_FRACTION:
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {  
	    sum_sv += (sv[w]-sv_p)*fv->c[w];
	  }
	sum_sv += sv_p;
	rho = 1/sum_sv ;

	if (d_rho != NULL && pd->v[MASS_FRACTION]) {
	  phi_ptr =  bf[MASS_FRACTION]->phi;
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      for (w = 0; w < pd->Num_Species_Eqn; w++)
		{
		  d_rho->C[w][j]= -rho*rho*(sv[w]-sv_p)
		    *phi_ptr[j] ; 
		}
	    }
	}
	break;
      case SPECIES_DENSITY:
	for (w = 0; w < pd->Num_Species_Eqn; w++)
	  {
	    rho += (1.- sv[w]/sv_p)*fv->c[w];
	  }
	rho += 1./sv_p;

	if (d_rho != NULL && pd->v[MASS_FRACTION]) {
	  phi_ptr =  bf[MASS_FRACTION]->phi;
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      for (w = 0; w < pd->Num_Species_Eqn; w++)
		{
		  d_rho->C[w][j] = (1.- sv[w]/sv_p)*phi_ptr[j];
		}
	    }
	}
	break;
      default:
	for (w=0; w<pd->Num_Species_Eqn; w++)
	  {  
	    sum_sv += (sv[w]-sv_p)*fv->c[w];
	  }
	sum_sv += sv_p;
	rho = 1/sum_sv ; 
      
	if (d_rho != NULL && pd->v[MASS_FRACTION]) {
	  phi_ptr =  bf[MASS_FRACTION]->phi;
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      for (w = 0; w < pd->Num_Species_Eqn; w++)
		{
		  d_rho->C[w][j]= -rho*rho*(sv[w]-sv_p)
		    *phi_ptr[j] ; 
		}
	    }
	}
	break;
      }
    }
  else if (mp->DensityModel == DENSITY_CONSTANT_LAST_CONC) {
    if (matrl_species_var_type != SPECIES_CONCENTRATION) {
      EH(-1,"unimplemented");
    }
    if (pd->Num_Species >  pd->Num_Species_Eqn) {
      w =  pd->Num_Species_Eqn;
      rho = mp->molecular_weight[w] * fv->c[w];
    } else {
      rho = 0.0;
    }
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      rho += mp->molecular_weight[w] * fv->c[w];
    }
    if (d_rho) {
      phi_ptr = bf[MASS_FRACTION]->phi;
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
	for (j = 0; j < ei->dof[MASS_FRACTION]; j++) {
	  d_rho->C[w][j] = mp->molecular_weight[w] * phi_ptr[j];
	}
      }
    }
  }
  else
    {
      EH(-1,"Unrecognized density model");
    }

  if( ls != NULL && mp->mp2nd != NULL &&
      ( mp->DensityModel != LEVEL_SET) &&
      ( mp->mp2nd->DensityModel == CONSTANT) )
    {
      double factor;

      if ( d_rho == NULL )
	{
	  /* kludge for solidification tracking with phase function 0 */
	  if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	    {
	      ls_old = ls;
	      ls = pfd->ls[0];
	      rho = ls_modulate_property( rho,
					  mp->mp2nd->density_phase[0],
					  ls->Length_Scale,
					  (double) mp->mp2nd->densitymask[0],
					  (double) mp->mp2nd->densitymask[1],
					  NULL,
					  &factor );
	      ls = ls_old;
	    }
          rho = ls_modulate_property( rho,
				      mp->mp2nd->density,
				      ls->Length_Scale,
				      (double) mp->mp2nd->densitymask[0],
				      (double) mp->mp2nd->densitymask[1],
				      NULL,
				      &factor );
	}
      else
	{
	  /* kludge for solidification tracking with phase function 0 */
	  if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	    {
	      ls_old = ls;
	      ls = pfd->ls[0];
	      rho = ls_modulate_property( rho,
					  mp->mp2nd->density_phase[0],
					  ls->Length_Scale,
					  (double) mp->mp2nd->densitymask[0],
					  (double) mp->mp2nd->densitymask[1],
					  d_rho->F,
					  &factor );
	      ls = ls_old;
	    }
          rho = ls_modulate_property( rho,
				      mp->mp2nd->density,
				      ls->Length_Scale,
				      (double) mp->mp2nd->densitymask[0],
				      (double) mp->mp2nd->densitymask[1],
				      d_rho->F,
				      &factor );
          if ( pd->v[MASS_FRACTION] )
	    {
	      for (w = 0; w < pd->Num_Species_Eqn; w++)
	        {
	          for (j = 0; j < ei->dof[MASS_FRACTION]; j++)
		    {
		      d_rho->C[w][j] *= factor;
		    }
	        }
	    }
		  
          if ( pd->v[TEMPERATURE] )
	    {
	      for (j = 0; j < ei->dof[TEMPERATURE]; j++)
	        {
	          d_rho->T[j] *= factor;
	        }
	    }
        }
    }
  
  return (rho);
}

double 
conductivity( CONDUCTIVITY_DEPENDENCE_STRUCT *d_k,
	      dbl time )

     
    /**************************************************************************
     *
     * conductivity
     *
     *   Calculate the thermal conductivity and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_k->T[j]    -> derivative of conductivity wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_k->C[k][j] -> derivative of conductivity wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_k->X[a][j] -> derivative of conductivity wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_k->F[j]    -> derivative of conductivity wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the conductivity
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double k, Tref, tmp;
  struct Level_Set_Data *ls_old;

  int i_therm_cond;

  k =0.;

  if ( d_k != NULL )
    {
      memset( d_k->T, 0, sizeof(double)*MDE);
      memset( d_k->X, 0, sizeof(double)*DIM*MDE);
      memset( d_k->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_k->F, 0, sizeof(double)*MDE);
    }

  if(mp->ConductivityModel == USER )
    {

      usr_thermal_conductivity(mp->u_thermal_conductivity,time);

      k   = mp->thermal_conductivity;

      var = TEMPERATURE;
      if( pd->v[TEMPERATURE] && d_k != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_k->T[j]= mp->d_thermal_conductivity[var]*bf[var]->phi[j];
	    }
	}

      if (pd->v[MESH_DISPLACEMENT1] && d_k != NULL )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_k->X[a][j] =mp->d_thermal_conductivity[var]*bf[var]->phi[j];
		}
	    }
	}

      if (pd->v[MASS_FRACTION] && d_k != NULL)
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_k->C[w][j] =mp->d_thermal_conductivity[var_offset]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->ConductivityModel == CONSTANT )
    {
      k    = mp->thermal_conductivity;
    }
  else if (mp->ConductivityModel == THERMAL_HEAT )
    {
      Tref = mp->u_thermal_conductivity[4];
      tmp = fv->T - Tref;
      mp->thermal_conductivity = mp->u_thermal_conductivity[0] +
		tmp*(mp->u_thermal_conductivity[1] + 
		tmp*(mp->u_thermal_conductivity[2] + 
		tmp*mp->u_thermal_conductivity[3]));
      k    = mp->thermal_conductivity;
      var = TEMPERATURE;
      if( pd->v[var] && d_k != NULL )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_k->T[j] = (mp->u_thermal_conductivity[1] + 
		tmp*(2.*mp->u_thermal_conductivity[2] + 
		tmp*3.*mp->u_thermal_conductivity[3]))*bf[var]->phi[j];
	    }
	}
    }
  else if (mp->ConductivityModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->thermal_conductivity_tableid];
      apply_table_mp( &mp->thermal_conductivity, table_local );
      k    = mp->thermal_conductivity;

      if ( d_k != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_k->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->ConductivityModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_thermal_conductivity[0], 
			     mp->u_thermal_conductivity[1], 
			     mp->u_thermal_conductivity[2], 
			     &mp->thermal_conductivity,
			     &mp->d_thermal_conductivity[FILL] );

      k = mp->thermal_conductivity;

      var = FILL;
      if( pd->v[var] && d_k != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_k->F[j] = mp->d_thermal_conductivity[var]*bf[var]->phi[j];
	    }
	}

    }
  else if (mp->ConductivityModel == EXTERNAL_FIELD)
    {
      i_therm_cond = mp->thermal_cond_external_field; 
      EH(i_therm_cond, "Thermal cond. external field not found!");
      k = fv->external_field[i_therm_cond]*mp->u_thermal_conductivity[0];
    }
  else
    {
      EH(-1,"Unrecognized thermal conductivity model");
    }

  if( ls != NULL && 
      mp->ConductivityModel != LEVEL_SET &&
      mp->ConductivityModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->ThermalConductivityModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
      /* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	  ls_old = ls;
	  ls = pfd->ls[0];
	  k = ls_modulate_thermalconductivity( k,
					       mp->mp2nd->thermalconductivity_phase[0],
					       ls->Length_Scale,
					       (double) mp->mp2nd->thermalconductivitymask[0],
					       (double) mp->mp2nd->thermalconductivitymask[1],
					       d_k );
	  ls = ls_old;
	}
      k = ls_modulate_thermalconductivity( k,
					   mp->mp2nd->thermalconductivity,
					   ls->Length_Scale,
					   (double) mp->mp2nd->thermalconductivitymask[0],
					   (double) mp->mp2nd->thermalconductivitymask[1],
					   d_k );

    }

  return(k);
}



double
heat_capacity( HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp,
	      dbl time )
    /**************************************************************************
     *
     * heat capacity
     *
     *   Calculate the thermal heat capacity and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_Cp->T[j]    -> derivative of heat capacity wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_Cp->V[j]    -> derivative of heat capacity wrt the jth
     *                          volage unknown in an element
     *    dbl d_Cp->v[a][j] -> derivative of heat capacity wrt the jth
     *                          velocity  in the ath direction.
     *    dbl d_Cp->C[k][j] -> derivative of heat capacity wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_Cp->X[a][j] -> derivative of heat capacity wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_Cp->F[j]    -> derivative of heat capacity wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the heat capacity
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double Cp, T_offset;
  struct Level_Set_Data *ls_old;



  Cp =0.;

  if ( d_Cp != NULL )
    {
      memset( d_Cp->T, 0, sizeof(double)*MDE);
      memset( d_Cp->V, 0, sizeof(double)*MDE);
      memset( d_Cp->X, 0, sizeof(double)*DIM*MDE);
      memset( d_Cp->v, 0, sizeof(double)*DIM*MDE);
      memset( d_Cp->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_Cp->F, 0, sizeof(double)*MDE);
    }

  if(mp->HeatCapacityModel == USER )
    {
      usr_heat_capacity(mp->u_heat_capacity,time);
      Cp = mp->heat_capacity;

      var = TEMPERATURE;
      if ( d_Cp != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_Cp->T[j]= mp->d_heat_capacity[var]*bf[var]->phi[j];
            }
        }

      if ( d_Cp != NULL && pd->v[MESH_DISPLACEMENT1] )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_Cp->X[a][j] =mp->d_heat_capacity[var]*bf[var]->phi[j];
		}
	    }
	}

      if ( d_Cp != NULL && pd->v[FILL] )
	{
	  var = FILL;

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_Cp->F[j] =mp->d_heat_capacity[var]*bf[var]->phi[j];
	    }
	}

      if ( d_Cp != NULL && pd->v[MASS_FRACTION] )
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_Cp->C[w][j] =mp->d_heat_capacity[var_offset]*bf[var]->phi[j];
		}
	    }
	}

    }
  else if (mp->HeatCapacityModel == CONSTANT )
    {
      Cp    = mp->heat_capacity;
    }
  else if (mp->HeatCapacityModel == THERMAL_HEAT )
    {
      T_offset = mp->u_heat_capacity[4];
      mp->heat_capacity = mp->u_heat_capacity[0] +
			mp->u_heat_capacity[1]*(fv->T+T_offset) +
			mp->u_heat_capacity[2]*SQUARE(fv->T+T_offset) +
			mp->u_heat_capacity[3]/SQUARE(fv->T+T_offset);
      Cp    = mp->heat_capacity;
      var = TEMPERATURE;
      if ( d_Cp != NULL && pd->v[var] )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_Cp->T[j] = ( mp->u_heat_capacity[1] +
			2.*(fv->T+T_offset)*mp->u_heat_capacity[2]
			- 2.*mp->u_heat_capacity[3]/(SQUARE(fv->T+T_offset)*(fv->T+T_offset)))
			*bf[var]->phi[j];
	    }
	}
    }
  else if (mp->HeatCapacityModel == ENTHALPY )
    {
      Cp = enthalpy_heat_capacity_model( d_Cp );
    }
  else if (mp->HeatCapacityModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_heat_capacity[0], 
			     mp->u_heat_capacity[1], 
			     mp->u_heat_capacity[2], 
			     &mp->heat_capacity,
			     &mp->d_heat_capacity[FILL] );

      Cp = mp->heat_capacity;

      var = FILL;
      if ( d_Cp != NULL && pd->v[var] )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_Cp->F[j] = mp->d_heat_capacity[var]*bf[var]->phi[j];
	    }
	}
      
    }
  else if (mp->HeatCapacityModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->heat_capacity_tableid];
      apply_table_mp(&mp->heat_capacity, table_local);
      Cp    = mp->heat_capacity;

      if ( d_Cp != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_Cp->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
            }
	}
    }
  else
    {
      EH(-1,"Unrecognized heat capacity model");
    }

  if( ls != NULL && 
      mp->HeatCapacityModel != LEVEL_SET &&
      mp->HeatCapacityModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->HeatCapacityModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        Cp = ls_modulate_heatcapacity( Cp,
				  mp->mp2nd->heatcapacity_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->heatcapacitymask[0],
				  (double) mp->mp2nd->heatcapacitymask[1],
				  d_Cp );
	ls = ls_old;
	}
      Cp = ls_modulate_heatcapacity( Cp,
				     mp->mp2nd->heatcapacity,
				     ls->Length_Scale,
				     (double) mp->mp2nd->heatcapacitymask[0],
				     (double) mp->mp2nd->heatcapacitymask[1],
				     d_Cp );

    }

  return (Cp);
}
double
ls_modulate_thermalconductivity ( double k1,
			double  k2,
			double width,
			double pm_minus,
			double pm_plus,
			CONDUCTIVITY_DEPENDENCE_STRUCT *d_k )
{
  double factor;
  int i,a, w, var;

  if ( d_k == NULL )
    {
      k1 = ls_modulate_property( k1, k2, width, pm_minus, pm_plus, NULL, &factor);
      return ( k1 );
    }

  k1 = ls_modulate_property( k1, k2, width, pm_minus, pm_plus, d_k->F, &factor);

  if ( pd->v[var=TEMPERATURE ] )
    {
      for(i=0; i<ei->dof[var]; i++)
	{
	  d_k->T[i] *= factor;
	}
    }
      
  if ( pd->v[var=MASS_FRACTION ] )
    {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  for( i=0; i<ei->dof[var]; i++) 
	    {
	      d_k->C[w][i] *= factor;
	    }
	}
    }

  if( pd->v[var=MESH_DISPLACEMENT1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
	{
	  for( i=0; i<ei->dof[var]; i++ )
	    {
	      d_k->X [a][i] *= factor;
	    }
	}
    }
  return ( k1 );
  
}
double
ls_modulate_heatcapacity ( double Cp1,
			double  Cp2,
			double width,
			double pm_minus,
			double pm_plus,
			HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp )
{
  double factor;
  int i,a, w, var;


  if ( d_Cp == NULL )
    {
      Cp1 = ls_modulate_property( Cp1, Cp2, width, pm_minus, pm_plus, NULL, &factor);
      return ( Cp1 );
    }

  Cp1 = ls_modulate_property( Cp1, Cp2, width, pm_minus, pm_plus, d_Cp->F, &factor);

  if ( pd->v[var=TEMPERATURE ] )
    {
      for(i=0; i<ei->dof[var]; i++)
	{
	  d_Cp->T[i] *= factor;
	}
    }
      
  if ( pd->v[var=MASS_FRACTION ] )
    {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	{
	  for( i=0; i<ei->dof[var]; i++) 
	    {
	      d_Cp->C[w][i] *= factor;
	    }
	}
    }

  if( pd->v[var=VELOCITY1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
	{
	  for( i=0; i<ei->dof[var]; i++ )
	    {
	      d_Cp->v[a][i] *= factor;
	    }
	}
    }

  if( pd->v[var=MESH_DISPLACEMENT1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
	{
	  for( i=0; i<ei->dof[var]; i++ )
	    {
	      d_Cp->X [a][i] *= factor;
	    }
	}
    }

  if( pd->v[var=VOLTAGE] )
    {
      for( i=0; i<ei->dof[var]; i++ )
	{
	  d_Cp->V[i] *= factor;
	}
    }

  return ( Cp1 );
  
}


double 
acoustic_impedance( CONDUCTIVITY_DEPENDENCE_STRUCT *d_R,
	      dbl time )

     
    /**************************************************************************
     *
     * acoustic impedance
     *
     *   Calculate the acoustic impedance and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the acoustic impedance
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double R;
  struct Level_Set_Data *ls_old;
  R =0.;

  if ( d_R != NULL )
    {
      memset( d_R->T, 0, sizeof(double)*MDE);
      memset( d_R->X, 0, sizeof(double)*DIM*MDE);
      memset( d_R->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_R->F, 0, sizeof(double)*MDE);
    }

  if(mp->Acoustic_ImpedanceModel == USER )
    {

	EH(-1,"user acoustic impedance code not ready yet");
	/*
	  usr_acoustic_impedance(mp->u_acoustic_impedance,time);
	*/

      R   = mp->acoustic_impedance;

      var = TEMPERATURE;
      if( pd->v[TEMPERATURE] && d_R != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_R->T[j]= mp->d_acoustic_impedance[var]*bf[var]->phi[j];
	    }
	}

      if (pd->v[MESH_DISPLACEMENT1] && d_R != NULL )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_R->X[a][j] =mp->d_acoustic_impedance[var]*bf[var]->phi[j];
		}
	    }
	}

      if (pd->v[MASS_FRACTION] && d_R != NULL)
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_R->C[w][j] =mp->d_acoustic_impedance[var_offset]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->Acoustic_ImpedanceModel == CONSTANT )
    {
      R    = mp->acoustic_impedance;
    }
  else if (mp->Acoustic_ImpedanceModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->acoustic_impedance_tableid];
      apply_table_mp( &mp->acoustic_impedance, table_local);
      R    = mp->acoustic_impedance;

      if ( d_R != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_R->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->Acoustic_ImpedanceModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_acoustic_impedance[0], 
			     mp->u_acoustic_impedance[1], 
			     mp->u_acoustic_impedance[2], 
			     &mp->acoustic_impedance,
			     &mp->d_acoustic_impedance[FILL] );

      R = mp->acoustic_impedance;

      var = FILL;
      if( pd->v[var] && d_R != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_R->F[j] = mp->d_acoustic_impedance[var]*bf[var]->phi[j];
	    }
	}

    }
  else
    {
      EH(-1,"Unrecognized acoustic impedance model");
    }

  if( ls != NULL && 
      mp->Acoustic_ImpedanceModel != LEVEL_SET &&
      mp->Acoustic_ImpedanceModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->AcousticImpedanceModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        R = ls_modulate_thermalconductivity( R,
				  mp->mp2nd->acousticimpedance_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->acousticimpedancemask[0],
				  (double) mp->mp2nd->acousticimpedancemask[1],
				  d_R );
	ls = ls_old;
	}
      R = ls_modulate_thermalconductivity( R,
				  mp->mp2nd->acousticimpedance,
				  ls->Length_Scale,
				  (double) mp->mp2nd->acousticimpedancemask[0],
				  (double) mp->mp2nd->acousticimpedancemask[1],
				  d_R );

    }

  return(R);
}


double 
wave_number( CONDUCTIVITY_DEPENDENCE_STRUCT *d_k,
	      dbl time )
{
  int var, j, i;
  double k;
  struct Level_Set_Data *ls_old;



  k =0.;

  if ( d_k != NULL )
    {
      memset( d_k->T, 0, sizeof(double)*MDE);
      memset( d_k->X, 0, sizeof(double)*DIM*MDE);
      memset( d_k->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_k->F, 0, sizeof(double)*MDE);
    }

  if (mp->wave_numberModel == CONSTANT )
    {
      k    = mp->wave_number;
    }
  else if (mp->wave_numberModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->wave_number_tableid];
      apply_table_mp( &mp->wave_number, table_local );
      k    = mp->wave_number;

      if ( d_k != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_k->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->wave_numberModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_wave_number[0], 
			     mp->u_wave_number[1], 
			     mp->u_wave_number[2], 
			     &mp->wave_number,
			     &mp->d_wave_number[FILL] );

      k = mp->wave_number;

      var = FILL;
      if( pd->v[var] && d_k != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_k->F[j] = mp->d_wave_number[var]*bf[var]->phi[j];
	    }
	}

    }
  else
    {
      EH(-1,"Unrecognized acoustic wave number model");
    }

  if( ls != NULL && 
      mp->wave_numberModel != LEVEL_SET &&
      mp->wave_numberModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->wavenumberModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        k = ls_modulate_thermalconductivity( k,
				  mp->mp2nd->wavenumber_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->wavenumbermask[0],
				  (double) mp->mp2nd->wavenumbermask[1],
				  d_k );
	ls = ls_old;
	}
      k = ls_modulate_thermalconductivity( k,
				  mp->mp2nd->wavenumber,
				  ls->Length_Scale,
				  (double) mp->mp2nd->wavenumbermask[0],
				  (double) mp->mp2nd->wavenumbermask[1],
				  d_k );
    }
  return(k);
}

/*********************************************************************************/
double 
acoustic_absorption( CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha,
	      dbl time )

     
    /**************************************************************************
     *
     * acoustic absorption
     *
     *   Calculate the acoustic absorption and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the acoustic absorption
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double alpha;
  struct Level_Set_Data *ls_old;
  alpha =0.;

  if ( d_alpha != NULL )
    {
      memset( d_alpha->T, 0, sizeof(double)*MDE);
      memset( d_alpha->X, 0, sizeof(double)*DIM*MDE);
      memset( d_alpha->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_alpha->F, 0, sizeof(double)*MDE);
    }

  if(mp->Acoustic_AbsorptionModel == USER )
    {

         mp->acoustic_absorption = mp->u_acoustic_absorption[0];
         mp->d_acoustic_absorption[TEMPERATURE] = 0;
      for ( a=0 ; a<DIM ; a++)
        {
         mp->d_acoustic_absorption[MESH_DISPLACEMENT1+a] = 0;
        }
      for ( a=0 ; a<pd->Num_Species_Eqn ; a++)
        {
         mp->d_acoustic_absorption[MAX_VARIABLE_TYPES+a] = mp->u_acoustic_absorption[a+1];
         mp->acoustic_absorption += mp->u_acoustic_absorption[a+1]*fv->c[a];
        }
      alpha   = mp->acoustic_absorption;

      var = TEMPERATURE;
      if( pd->v[TEMPERATURE] && d_alpha != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_alpha->T[j]= mp->d_acoustic_absorption[var]*bf[var]->phi[j];
	    }
	}

      if (pd->v[MESH_DISPLACEMENT1] && d_alpha != NULL )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_alpha->X[a][j] =mp->d_acoustic_absorption[var]*bf[var]->phi[j];
		}
	    }
	}

      if (pd->v[MASS_FRACTION] && d_alpha != NULL)
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_alpha->C[w][j] =mp->d_acoustic_absorption[var_offset]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->Acoustic_AbsorptionModel == CONSTANT )
    {
      alpha    = mp->acoustic_absorption;
    }
  else if (mp->Acoustic_AbsorptionModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->acoustic_absorption_tableid];
      apply_table_mp( &mp->acoustic_absorption, table_local );
      alpha    = mp->acoustic_absorption;

      if ( d_alpha != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_alpha->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->Acoustic_AbsorptionModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_acoustic_absorption[0], 
			     mp->u_acoustic_absorption[1], 
			     mp->u_acoustic_absorption[2], 
			     &mp->acoustic_absorption,
			     &mp->d_acoustic_absorption[FILL] );

      alpha = mp->acoustic_absorption;

      var = FILL;
      if( pd->v[var] && d_alpha != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_alpha->F[j] = mp->d_acoustic_absorption[var]*bf[var]->phi[j];
	    }
	}

    }
  else
    {
      EH(-1,"Unrecognized acoustic absorption model");
    }

  if( ls != NULL && 
      mp->Acoustic_AbsorptionModel != LEVEL_SET &&
      mp->Acoustic_AbsorptionModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->AcousticAbsorptionModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        alpha = ls_modulate_thermalconductivity( alpha,
				  mp->mp2nd->acousticabsorption_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->acousticabsorptionmask[0],
				  (double) mp->mp2nd->acousticabsorptionmask[1],
				  d_alpha );
	ls = ls_old;
	}
      alpha = ls_modulate_thermalconductivity( alpha,
				  mp->mp2nd->acousticabsorption,
				  ls->Length_Scale,
				  (double) mp->mp2nd->acousticabsorptionmask[0],
				  (double) mp->mp2nd->acousticabsorptionmask[1],
				  d_alpha );

    }

  return(alpha);
}

/*********************************************************************************/
double 
light_absorption( CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha,
	      dbl time )

     
    /**************************************************************************
     *
     * light absorption
     *
     *   Calculate the light absorption and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the light absorption
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double alpha;
  struct Level_Set_Data *ls_old;
  alpha =0.;

  if ( d_alpha != NULL )
    {
      memset( d_alpha->T, 0, sizeof(double)*MDE);
      memset( d_alpha->X, 0, sizeof(double)*DIM*MDE);
      memset( d_alpha->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_alpha->F, 0, sizeof(double)*MDE);
    }

  if(mp->Light_AbsorptionModel == USER )
    {

         mp->light_absorption = mp->u_light_absorption[0];
         mp->d_light_absorption[TEMPERATURE] = 0;
      for ( a=0 ; a<DIM ; a++)
        {
         mp->d_light_absorption[MESH_DISPLACEMENT1+a] = 0;
        }
      for ( a=0 ; a<pd->Num_Species_Eqn ; a++)
        {
         mp->d_light_absorption[MAX_VARIABLE_TYPES+a] = mp->u_light_absorption[a+1];
         mp->light_absorption += mp->u_light_absorption[a+1]*fv->c[a];
        }
      alpha   = mp->light_absorption;

      var = TEMPERATURE;
      if( pd->v[TEMPERATURE] && d_alpha != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_alpha->T[j]= mp->d_light_absorption[var]*bf[var]->phi[j];
	    }
	}

      if (pd->v[MESH_DISPLACEMENT1] && d_alpha != NULL )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_alpha->X[a][j] =mp->d_light_absorption[var]*bf[var]->phi[j];
		}
	    }
	}

      if (pd->v[MASS_FRACTION] && d_alpha != NULL)
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w +1;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_alpha->C[w][j] =mp->d_light_absorption[var_offset]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->Light_AbsorptionModel == CONSTANT )
    {
      alpha    = mp->light_absorption;
    }
  else if (mp->Light_AbsorptionModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->light_absorption_tableid];
      apply_table_mp( &mp->light_absorption, table_local );
      alpha    = mp->light_absorption;

      if ( d_alpha != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_alpha->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->Light_AbsorptionModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_light_absorption[0], 
			     mp->u_light_absorption[1], 
			     mp->u_light_absorption[2], 
			     &mp->light_absorption,
			     &mp->d_light_absorption[FILL] );

      alpha = mp->light_absorption;

      var = FILL;
      if( pd->v[var] && d_alpha != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_alpha->F[j] = mp->d_light_absorption[var]*bf[var]->phi[j];
	    }
	}

    }
  else
    {
      EH(-1,"Unrecognized light absorption model");
    }

  if( ls != NULL && 
      mp->Light_AbsorptionModel != LEVEL_SET &&
      mp->Light_AbsorptionModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->LightAbsorptionModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        alpha = ls_modulate_thermalconductivity( alpha,
				  mp->mp2nd->lightabsorption_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->lightabsorptionmask[0],
				  (double) mp->mp2nd->lightabsorptionmask[1],
				  d_alpha );
	ls = ls_old;
	}
      alpha = ls_modulate_thermalconductivity( alpha,
				  mp->mp2nd->lightabsorption,
				  ls->Length_Scale,
				  (double) mp->mp2nd->lightabsorptionmask[0],
				  (double) mp->mp2nd->lightabsorptionmask[1],
				  d_alpha );

    }

  return(alpha);
}


/*********************************************************************************/
double 
refractive_index( CONDUCTIVITY_DEPENDENCE_STRUCT *d_n,
	      dbl time )

     
    /**************************************************************************
     *
     * refractive index
     *
     *   Calculate the refractive index and its derivatives wrt to nodal dofs
     *   at the local gauss point
     *
     * Output
     * -----
     *    dbl d_R->T[j]    -> derivative of ref. index wrt the jth
     *                          Temperature unknown in an element
     *    dbl d_R->C[k][j] -> derivative of ref. index wrt the jth
     *                          species unknown of ktype, k, in the element.
     *    dbl d_R->X[a][j] -> derivative of ref. index wrt the jth
     *                          mesh displacement in the ath direction.
     *    dbl d_R->F[j]    -> derivative of ref. index wrt the jth
     *                          FILL unknown in an element
     *
     *  Return
     * --------
     *    Actual value of the refractive index
     ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double n;
  struct Level_Set_Data *ls_old;
  n =0.;

  if ( d_n != NULL )
    {
      memset( d_n->T, 0, sizeof(double)*MDE);
      memset( d_n->X, 0, sizeof(double)*DIM*MDE);
      memset( d_n->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_n->F, 0, sizeof(double)*MDE);
    }

  if(mp->Refractive_IndexModel == USER )
    {

         mp->refractive_index = mp->u_refractive_index[0];
         mp->d_refractive_index[TEMPERATURE] = 0;
      for ( a=0 ; a<DIM ; a++)
        {
         mp->d_refractive_index[MESH_DISPLACEMENT1+a] = 0;
        }
      for ( a=0 ; a<pd->Num_Species_Eqn ; a++)
        {
         mp->d_refractive_index[MAX_VARIABLE_TYPES+a] = mp->u_refractive_index[a+1];
         mp->refractive_index += mp->u_refractive_index[a+1]*fv->c[a];
        }
      n   = mp->refractive_index;

      var = TEMPERATURE;
      if( pd->v[TEMPERATURE] && d_n != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_n->T[j]= mp->d_refractive_index[var]*bf[var]->phi[j];
	    }
	}

      if (pd->v[MESH_DISPLACEMENT1] && d_n != NULL )
	{
	  for ( a=0; a<dim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_n->X[a][j] =mp->d_refractive_index[var]*bf[var]->phi[j];
		}
	    }
	}

      if (pd->v[MASS_FRACTION] && d_n != NULL)
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_n->C[w][j] =mp->d_refractive_index[var_offset]*bf[var]->phi[j];
		}
	    }
	}
    }
  else if (mp->Refractive_IndexModel == CONSTANT )
    {
      n    = mp->refractive_index;
    }
  else if (mp->Refractive_IndexModel == TABLE )
    {
      struct  Data_Table *table_local; 
      table_local = MP_Tables[mp->refractive_index_tableid];
      apply_table_mp( &mp->refractive_index, table_local );
      n    = mp->refractive_index;

      if ( d_n != NULL )
        {
          for(i=0;i<table_local->columns-1;i++)
	    {
	      var = table_local->t_index[i];
	      /* currently only set up to vary w.r.t. temperature */
	      switch (var)
	        {
	        case TEMPERATURE:
	          for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_n->T[j]= table_local->slope[i]*bf[var]->phi[j];
		    }
	          break;
	        default:
		      EH(-1, "Variable function not yet implemented in material property table");
	        }
	    }
        }
    }
  else if(mp->Refractive_IndexModel == LEVEL_SET )
    {
      ls_transport_property( mp->u_refractive_index[0], 
			     mp->u_refractive_index[1], 
			     mp->u_refractive_index[2], 
			     &mp->refractive_index,
			     &mp->d_refractive_index[FILL] );

      n = mp->refractive_index;

      var = FILL;
      if( pd->v[var] && d_n != NULL )
	{

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_n->F[j] = mp->d_refractive_index[var]*bf[var]->phi[j];
	    }
	}

    }
  else
    {
      EH(-1,"Unrecognized refractive index model");
    }

  if( ls != NULL && 
      mp->Refractive_IndexModel != LEVEL_SET &&
      mp->Refractive_IndexModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->RefractiveIndexModel == CONSTANT )  /* Only Newtonian constitutive equation allowed for 2nd phase */
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        n = ls_modulate_thermalconductivity( n,
				  mp->mp2nd->refractiveindex_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->refractiveindexmask[0],
				  (double) mp->mp2nd->refractiveindexmask[1],
				  d_n );
	ls = ls_old;
	}
      n = ls_modulate_thermalconductivity( n,
				  mp->mp2nd->refractiveindex,
				  ls->Length_Scale,
				  (double) mp->mp2nd->refractiveindexmask[0],
				  (double) mp->mp2nd->refractiveindexmask[1],
				  d_n );

    }

  return(n);
}
/*********************************************************************************/
/*********************************************************************************/
/******************************************************************************
 * momentum_source_term(): Computes the body force term for the momentum balance.
 *
 * Input
 * -----
 *   f    == Body force
 *   df->T == Derivative w.r.t. temperature
 *   df->X == Derivative w.r.t. mesh displacements
 *   df->C == Derivative w.r.t. concentration
 *   df->v == Derivative w.r.t. velocity
 *   df->F == Derivative w.r.t. FILL
 *   df->E == Derivative w.r.t. electric field
 *
 ******************************************************************************/
int
momentum_source_term(dbl f[DIM],                   /* Body force. */
                     MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df,
		     dbl time)
{
  int j, a, b, w;
  int eqn, var, var_offset;
  int err;
  const int dim = pd->Num_Dim;
  int siz;
  int status = 0;
  double *phi;
  struct Level_Set_Data *ls_old;
  
  /* initialize everything to zero */
  siz = sizeof(double)*DIM;
  memset(f,0,siz);

  if ( df != NULL )
    {
      siz = sizeof(double)*DIM*MDE;
      memset(df->T,0,siz);
      memset(df->F,0,siz);
      memset(df->ars,0,siz);
  
      siz = sizeof(double)*DIM*DIM*MDE;
      memset(df->X,0,siz);
      memset(df->v,0,siz);
      memset(df->E,0,siz);

      siz = sizeof(double)*DIM*MAX_CONC*MDE;
      memset(df->C,0,siz);
	  
    }
  
  /****Momentum Source Model******/
  if(mp->MomentumSourceModel == USER)
    {
      err = usr_momentum_source(mp->u_momentum_source);
      
      for ( a=0; a<dim; a++)
	{
	  eqn = R_MOMENTUM1 + a;
	  if ( pd->e[eqn] & T_SOURCE )
	    {
	      f[a] = mp->momentum_source[a];
              var  = TEMPERATURE;
              if ( df != NULL && pd->v[var] )
                {
	          phi = bf[var]->phi;
                  for (j=0; j<ei->dof[var]; j++)
		    {
		      df->T[a][j] = mp->d_momentum_source[a][var] * phi[j];
		    }
                }
#ifdef COUPLED_FILL
	      if ( df != NULL && pd->v[FILL] )
		{
		  var = FILL;
		  phi = bf[var]->phi;
		  for ( j=0; j < ei->dof[var]; j++ )
		    {
		      df->F[a][j] = mp->d_momentum_source[a][var] * phi[j];
		    }
		}
#endif /* COUPLED_FILL */
	       if ( df != NULL && pd->v[MESH_DISPLACEMENT1] )
		{
		  for ( b=0; b<dim; b++)
		    {
		      var = MESH_DISPLACEMENT1+b;
		      phi = bf[var]->phi;
		      for (j=0; j<ei->dof[var]; j++)
			{
			  df->X[a][b][j] = mp->d_momentum_source[a][var] * phi[j];
			}
		    }
		}
	      
	      if ( df != NULL && pd->v[MASS_FRACTION] )
		{
		  for ( w=0; w<pd->Num_Species_Eqn; w++)
		    {
		      var = MASS_FRACTION;
		      var_offset = MAX_VARIABLE_TYPES + w;
		      phi = bf[var]->phi;
		      for (j=0; j<ei->dof[var]; j++)
			{
			  df->C[a][w][j] = mp->d_momentum_source[a][var_offset] * phi[j];
			}
		    }
		}
	      
	      if ( df != NULL && pd->v[VELOCITY1] )
		{
		  for ( b=0; b<DIM; b++)
		    {
		      var = VELOCITY1 + b;;
		      phi = bf[var]->phi;
      		      for (j=0; j<ei->dof[var]; j++)
			{
			  df->v[a][b][j] = mp->d_momentum_source[a][var] * phi[j];
			}
		    }
		}
	    }
	}
    }
  else if (mp->MomentumSourceModel == CONSTANT )
    {
      for ( a=0; a<dim; a++)
	{
	  eqn   = R_MOMENTUM1+a;			
	  if ( pd->e[eqn] & T_SOURCE )
	    {
	      f[a] = mp->momentum_source[a];
 	    }
	}
    }
  else if (mp->MomentumSourceModel == VARIABLE_DENSITY)
    {
      if( mp->DensityModel == SOLVENT_POLYMER)
	{
	  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
	  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
	  double rho = density(d_rho, time);
	  for ( a=0; a<dim; a++)
	    {
	      eqn   = R_MOMENTUM1+a;			
	      if ( pd->e[eqn] & T_SOURCE )
		{
		  f[a] = rho*mp->momentum_source[a];
		}
	      if ( df != NULL && pd->v[MASS_FRACTION] )
		{
		  for ( w=0; w<pd->Num_Species_Eqn; w++)
		    {
		      var = MASS_FRACTION;
		      for (j=0; j<ei->dof[var]; j++)
			{
			  df->C[a][w][j] = d_rho->C[w][j]*mp->momentum_source[a];
			}
		    }
		}
	    }
	}
      else if ( mp->DensityModel == DENSITY_FOAM ||
		mp->DensityModel == DENSITY_FOAM_TIME ||
		mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
	{
	  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
	  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
	  double rho = density(d_rho, time);
	  for ( a=0; a<dim; a++)
	    {
	      eqn   = R_MOMENTUM1+a;			
	      if (pd->e[eqn] & T_SOURCE)
		{
		  f[a] = rho*mp->momentum_source[a];
		}
	      if (df != NULL && pd->v[MASS_FRACTION])
		{
		  for (w = 0; w < pd->Num_Species_Eqn; w++)
		    {
		      var = MASS_FRACTION;
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  df->C[a][w][j] = d_rho->C[w][j]*mp->momentum_source[a];
			}
		    }
		}
	      if ( df != NULL && pd->v[TEMPERATURE] )
		{
		  var = TEMPERATURE;
		  for (j=0; j<ei->dof[var]; j++)
		    {
		      df->T[a][j] = d_rho->T[j]*mp->momentum_source[a];
		    }
		}
	      if ( df != NULL && pd->v[FILL] )
		{
		  var = FILL;
		  phi = bf[var]->phi;
		  for (j=0; j<ei->dof[var]; j++)
		    {
		      df->F[a][j] = d_rho->F[j]*mp->momentum_source[a];
		    }
		}
	    }
	}

    }
  else if (mp->MomentumSourceModel == SUSPENSION_PM)
    {
      err = suspension_pm_fluid_momentum_source(f, df);
      EH(err,"Problems in suspension_pm_fluid_momentum_source");
    }
  else if (mp->MomentumSourceModel == SUSPEND || mp->MomentumSourceModel == SUSPENSION)
    {
      err = suspend_momentum_source(f, df);
      EH(err,"Problems in suspend_momentum_source");
    }
  else if (mp->MomentumSourceModel == BOUSS )
    {
      err = bouss_momentum_source(f, df, 0, TRUE);
      EH(err,"Problems in bouss_momentum_source");
    }
  else if (mp->MomentumSourceModel == BOUSS_JXB )
    {
      err = bouss_momentum_source(f, df, 1, FALSE);
      EH(err,"Problems in bouss_momentum_source");
    }
  else if (mp->MomentumSourceModel == BOUSSINESQ )
    {
      err = bouss_momentum_source(f, df, 0, FALSE);
      EH(err,"Problems in bouss_momentum_source");
    }
  else if (mp->MomentumSourceModel == EHD_POLARIZATION )
    {
      err = EHD_POLARIZATION_source(f, df);
      EH(err,"Problems in EHD_POLARIZATION force routine");
    }
  else if (mp->MomentumSourceModel == GRAV_VIBRATIONAL )
    {
      err = gravity_vibrational_source(f, df, time);
      EH(err,"Problems in GRAV_VIBRATIONAL force routine");
    }
  else if (mp->MomentumSourceModel == FILL_SRC)
    {
      err = fill_momentum_source(f);
    }
  else if (mp->MomentumSourceModel == LEVEL_SET  && pd->v[FILL]  )
    {
      DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
      double rho = density( d_rho , time);

      for ( a=0; a < pd->Num_Dim; a++ )
		{
	  f[a] = mp->momentum_source[a] * rho;
#ifdef COUPLED_FILL
          if ( df != NULL)
            {
	      for ( j=0; j < ei->dof[FILL]; j++)
	        df->F[a][j] = mp->momentum_source[a] * d_rho->F[j];
            }
#endif /* COUPLED_FILL */
		}
	}
	else if ( mp->MomentumSourceModel == LEVEL_SET  && pd->v[PHASE1]  )
	{
      DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
      double rho = density( d_rho, time );
	  
	  for( a=0; a < pd->Num_Dim; a++)
	  {
		f[a] = mp->momentum_source[a] * rho;
		
		if( df != NULL )
		{
			for( b=0 ; b < pfd->num_phase_funcs; b++)
			{
				for( j=0; j < ei->dof[PHASE1]; j++)
				{
					df->pf[a][b][j] =  mp->momentum_source[a] * d_rho->pf[b][j];
				}
			}
		}
	  }
	}
  else if (mp->MomentumSourceModel == ACOUSTIC )
    {
       	for ( a=0; a<dim; a++)
 		{
 	  	eqn   = R_MOMENTUM1+a;			
 	  	if ( pd->e[eqn] & T_SOURCE )
 	    		{
/*  Graviational piece	*/
 	      		f[a] = mp->momentum_source[a];
/*  Acoustic Reynolds Stress piece	*/
			f[a] += mp->u_momentum_source[0]*fv->grad_ars[a];
  	    		}
 	      	if ( df != NULL && pd->v[ACOUS_REYN_STRESS] )
 			{
 			var = ACOUS_REYN_STRESS;
    			phi = bf[var]->phi;
    			for (j=0; j<ei->dof[var]; j++)
 				{
 		  		df->ars[a][j] = mp->u_momentum_source[0]*bf[var]->grad_phi[j][a];
 				}
 			}
     		}
     }
  else
    {
      EH(-1,"No such Navier-Stokes Source Model");
    }

  if( ls != NULL &&
      mp->mp2nd != NULL &&
      mp->MomentumSourceModel != LEVEL_SET &&
      mp->mp2nd->MomentumSourceModel == CONSTANT &&
      (pd->e[R_MOMENTUM1] & T_SOURCE) )
    {
	/* kludge for solidification tracking with phase function 0 */
      if(pfd != NULL && pd->e[R_EXT_VELOCITY])
	{
	ls_old = ls;
	ls = pfd->ls[0];
        ls_modulate_momentumsource( f,
				  mp->mp2nd->momentumsource_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->momentumsourcemask[0],
				  (double) mp->mp2nd->momentumsourcemask[1],
				  df );
	ls = ls_old;
	}
      ls_modulate_momentumsource ( f,
				   mp->mp2nd->momentumsource,
				   ls->Length_Scale,
				   (double) mp->mp2nd->momentumsourcemask[0],
				   (double) mp->mp2nd->momentumsourcemask[1],
				   df );
    }

  return(status);
}

int
ls_modulate_momentumsource ( double f1[DIM],
			     double f2[DIM],
			     double width,
			     double pm_minus,
			     double pm_plus,
			     MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df )

{
  int i, a, b, var;
  int dim = pd->Num_Dim;
  double factor;

  
  for( a = 0 ; a< dim ; a++ )
    {

      if ( df == NULL )
	{
	  f1[a] = ls_modulate_property( f1[a],
					f2[a],
					width,
					pm_minus,
					pm_plus,
					NULL ,
					&factor );

	  return (0);
	}
      else
	{
	  f1[a] = ls_modulate_property( f1[a],
					f2[a],
					width,
					pm_minus,
					pm_plus,
					df->F[a] ,
					&factor );

	

	  if( pd->v[var=TEMPERATURE] )
	    {
	      for ( i=0; i<ei->dof[var]; i++)
		{
		  df->T[a][i] *=factor;
		}
	    }

	  if(  pd->v[var=MESH_DISPLACEMENT1] )
	    {
	      for( b=0; b<dim; b++ )
		{
		  for( i=0 ; i<ei->dof[var+b]; i++)
		    {
		      df->X[a][b][i] *= factor;
		    }
		}
	    }

	  if( pd->v[var=VELOCITY1] )
	    {
	      for( b=0; b<dim; b++ )
		{
		  for( i=0 ; i<ei->dof[var]; i++)
		    {
		      df->v[a][b][i] *= factor;
		    }
		}
	    }

	  if(  pd->v[var=MASS_FRACTION] )
	    {
	      for( b=0; b<pd->Num_Species; b++ )
		{
		  for( i=0 ; i<ei->dof[var]; i++)
		    {
		      df->C[a][b][i] *= factor;
		    }
		}
	    }
	
	  if( pd->v[var=EFIELD1] )
	    {
	      for( b=0; b<dim; b++ )
		{
		  for( i=0 ; i<ei->dof[var+b]; i++)
		    {
		      df->E[a][b][i] *= factor;
		    }
		}
	    }	  
	}
    }
  return (0 );
}


void
apply_table_mp( double *func, struct Data_Table *table)
{
  int i;
  double interp_val,var1[1],slope;

  for(i=0;i<table->columns-1;i++) 
  {
    if(strcmp(table->t_name[i], "TEMPERATURE")==0)
    {
      var1[i]=fv->T;
    }
    else if(strcmp(table->t_name[i], "MASS_FRACTION")==0)
    {
      var1[i]=fv->c[table->species_eq];
    }
    else if(strcmp(table->t_name[i], "FAUX_PLASTIC")==0)
    {
      EH(-1, "Oops, I shouldn't be using this call for FAUX_PLASTIC.");
    }
    else
    {
      EH(-1, "Material Table Model Error-Unknown Function Column ");
    }
  }

  interp_val = interpolate_table( table, var1, &slope, NULL);
  *func = interp_val;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

double
interpolate_table( struct Data_Table *table, 
			double x[], 
			double *sloper,
			double dfunc_dx[])
  /* 
   *      A general routine that uses data supplied in a Data_Table
	  structure to compute the ordinate of the table that corresponds
	  to the abscissa supplied in x.
	  
	  Author: Thomas A. Baer, Org 9111
	  Date  : July 16, 1998
	  Revised: raroach October 12, 1999

      Parameters:
          table = pointer to Data_Table structure
	  x     = array of abscissa(s) where ordinate is to be evaluated
	  slope = d(ordinate)/dx 

      returns:
          value of ordinate at x and the slope there too.
   */
{
  int i, N, iinter, istartx, istarty;
  double func=0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3],xleft,xright;
  int ngrid1,ngrid2,ngrid3;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  uu=0.0;

  switch ( table->interp_method )
  {
  case LINEAR:      /* This is knucklehead linear interpolation scheme */

/*  maybe someday you would prefer a more logical fem-based interpolation */

      if(0)
      {
	for( i=0; i < N ; i++)
	{
	  cee = (x[0]-t[i])/(t[i+1] - t[i]);
	  if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-1))
	  {
	    phi[0]=-cee+1.;
	    phi[1]=cee;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1];
	    break;
	  }
	}
	table->slope[2] = 0.0;
      }
      else
      {
	if ( x[0] < t[0] )
	{
	  table->slope[0] = ( f[1] - f[0] ) / ( t[1] - t[0] );
	  func = f[0] +  (table->slope[0]) * ( x[0] - t[0] );
	}
      
	for( i=0; x[0] >= t[i] &&  i < N ; i++)
	{
	  if ( x[0] >= t[i] && x[0] < t[i+1] )
	  {
	    table->slope[0] = ( f[i+1] - f[i] ) / ( t[i+1] - t[i] );
	    func  = f[i] + (table->slope[0]) * ( x[0] - t[i] );
	  }
	}
	if ( x[0]>= t[N] )
	{
	  table->slope[0] = ( f[N] - f[N-1] ) / ( t[N] - t[N-1] );
	  func = f[N] + (table->slope[0]) * ( x[0] - t[N] ) ;
	}
	table->slope[1]=0.0;
	*sloper=table->slope[0];
      }
      break;

  case QUADRATIC:      /* quadratic lagrangian interpolation scheme */

      if(table->columns == 3)
      {
	for( i=0; i < N ; i+=2)
        {
          cee = (x[0]-t[i])/(t[i+2] - t[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1]+f[N+3+i]*phi[2];
	    break;
	  }
        }
        table->slope[2] = 0.0;
      }
      else
      {
	for( i=0; i < N ; i+=2)
        {
          cee = (x[0]-t[i])/(t[i+2] - t[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
	    func = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    phi[0]=4.*cee-3.;
	    phi[1]=-8.*cee+4.;
	    phi[2]=4.*cee-1.;
	    table->slope[0] = (f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2])
		/(t[i+2]-t[i]);
	    break;
	  }
        }
	table->slope[1]=0.0;
	*sloper=table->slope[0];
      }
      break;

  case QUAD_GP:      /* quadratic lagrangian interpolation scheme */

      if(table->columns == 3)
      {
	for( i=0; i < N ; i+=3)
        {
          xleft=(5.+sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.-sqrt(15.))/6.*t[i+2];
          xright=(5.-sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.+sqrt(15.))/6.*t[i+2];
          cee = (x[0]-xleft)/(xright-xleft);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-3))
	  {
	    phi[0]=(20.*cee*cee-2.*(sqrt(15.)+10.)*cee+sqrt(15.)+5.)/6.;
	    phi[1]=(-10.*cee*cee+10.*cee-1.)*2./3.;
	    phi[2]=(20.*cee*cee+2.*(sqrt(15.)-10.)*cee-sqrt(15.)+5.)/6.;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1]+f[N+3+i]*phi[2];
	    break;
	  }
        }
        table->slope[2] = 0.0;
      }
      else
      {
	for( i=0; i < N ; i+=3)
        {
          xleft=(5.+sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.-sqrt(15.))/6.*t[i+2];
          xright=(5.-sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.+sqrt(15.))/6.*t[i+2];
          cee = (x[0]-xleft)/(xright-xleft);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-3))
	  {
	    phi[0]=(20.*cee*cee-2.*(sqrt(15.)+10.)*cee+sqrt(15.)+5.)/6.;
	    phi[1]=(-10.*cee*cee+10.*cee-1.)*2./3.;
	    phi[2]=(20.*cee*cee+2.*(sqrt(15.)-10.)*cee-sqrt(15.)+5.)/6.;
	    func = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    phi[0]=(20.*cee-sqrt(15.)+10.)/3.;
	    phi[1]=(-40.*cee+20.)/3.;
	    phi[2]=(20.*cee+sqrt(15.)-10.)/3.;
	    table->slope[0] = (f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2])
		/(xright-xleft);
	    break;
	  }
        }
	table->slope[1]=0.0;
	*sloper=table->slope[0];
      }
      break;

  case BIQUADRATIC:      /* biquadratic lagrangian interpolation scheme */

      ngrid1 = table->tablelength/table->ngrid;
      if(table->columns == 5)
      {
 	table->slope[0] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,f
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
 	table->slope[1] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,&f[N+1]
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
 	table->slope[2] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,&f[2*N+2]
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
      }
      else
      {
 	func = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,f
 				  ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
      }
      break;




  case BILINEAR:  /* BILINEAR Interpolation Scheme */
      /*Find Interval of Different Values of Abscissa #1*/
      iinter=0;
      for(i=0;i<N;i++)
      {
	if(table->t[i]!=table->t[i+1] && iinter==0)
	{ iinter=i+1;}
      }
      if(iinter == 1) 
      { 
	fprintf(stderr," MP Interpolate Error - Need more than 1 point per set");
	EH(-1, "Table interpolation not implemented");
      }

      istartx=iinter;

      for(i=iinter;t[i]<=x[0] && i<N-iinter+1;i=i+iinter)
      { istartx=i+iinter;}

      istarty=istartx;

      for(i=istartx+1;t2[i]<=x[1] && i<istartx+iinter-1;i++)
      { istarty=i;}

      y1=f[istarty];
      y2=f[istarty+1];
      y3=f[istarty+1-iinter];
      y4=f[istarty-iinter];

      tt=(x[1]-t2[istarty])/(t2[istarty+1]-t2[istarty]);
      uu=(x[0]-t[istarty])/(t[istarty+1-iinter]-t[istarty]);

      func=(1.-tt)*(1.-uu)*y1+tt*(1.-uu)*y2+tt*uu*y3+(1.-tt)*uu*y4;

      table->slope[1]= ( f[istarty+1] - f[istarty] ) / ( t2[istarty+1] - t2[istarty] );
      table->slope[0]= ( f[istarty+1] - f[istarty+1-iinter] ) / ( t[istarty+1] - t[istarty+1-iinter] );
      *sloper=table->slope[0];
	if(dfunc_dx != NULL)
	   {
		dfunc_dx[0] = table->slope[0];
		dfunc_dx[1] = table->slope[1];
	   }
      break;

   case TRILINEAR:      /* trilinear lagrangian interpolation scheme */
 
       ngrid1 = table->ngrid;
       ngrid2 = table->ngrid2/table->ngrid;
       ngrid3 = table->tablelength/table->ngrid2;
       func = quad_isomap_invert(x[0],x[1],x[2],t,t2,t3,f
 				  ,ngrid1,ngrid2,ngrid3,1,3,dfunc_dx);
       break;
 
   case TRIQUADRATIC:      /* triquadratic lagrangian interpolation scheme */
 
       ngrid1 = table->ngrid;
       ngrid2 = table->ngrid2/table->ngrid;
       ngrid3 = table->tablelength/table->ngrid2;
       func = quad_isomap_invert(x[0],x[1],x[2],t,t2,t3,f
 				  ,ngrid1,ngrid2,ngrid3,2,3,dfunc_dx);
       break;

  default:
      EH(-1, "Table interpolation order not implemented");
  }

  return(func);
}


/***************************************************************************/

double
table_distance_search( struct Data_Table *table, 
			double x[], 
			double *sloper,
			double dfunc_dx[])
  /* 
   *      A routine to find the shortest distance from a point
   *      to a table represented curve
	  
	  Author: Robert B. Secor
	  Date  : December 5, 2014
	  Revised: 

      Parameters:
          table = pointer to Data_Table structure
	  x     = array of point coordinates
	  slope = d(ordinate)/dx 

      returns:
          value of distance
   */
{
  int i, N, iinter, istartx, istarty,basis,ordinate,i_min=0,elem,iter;
  double func=0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3],xleft,xright,delta,epstol=0.0001,update,point;
  double phic[3]={0,0,0};
  int ngrid1,ngrid2,ngrid3;
  double dist, dist_min = BIG_PENALTY;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  if(table->t_index[0] == MESH_POSITION1)
     {basis = 0;ordinate =1;}
  else
  if(table->t_index[0] == MESH_POSITION2)
     {basis = 1;ordinate =0;}
  else
     {basis = 0;ordinate =1;}
  uu=0.0;

  switch ( table->interp_method )
  {
  case LINEAR:      /* This is knucklehead linear interpolation scheme */

/*  maybe someday you would prefer a more logical fem-based interpolation */

      if(0)
      {
	for( i=0; i < N ; i++)
	{
	  cee = (x[0]-t[i])/(t[i+1] - t[i]);
	  if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-1))
	  {
	    phi[0]=-cee+1.;
	    phi[1]=cee;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1];
	    break;
	  }
	}
	table->slope[2] = 0.0;
      }
      else
      {
	if ( x[0] < t[0] )
	{
	  table->slope[0] = ( f[1] - f[0] ) / ( t[1] - t[0] );
	  func = f[0] +  (table->slope[0]) * ( x[0] - t[0] );
	}
      
	for( i=0; x[0] >= t[i] &&  i < N ; i++)
	{
	  if ( x[0] >= t[i] && x[0] < t[i+1] )
	  {
	    table->slope[0] = ( f[i+1] - f[i] ) / ( t[i+1] - t[i] );
	    func  = f[i] + (table->slope[0]) * ( x[0] - t[i] );
	  }
	}
	if ( x[0]>= t[N] )
	{
	  table->slope[0] = ( f[N] - f[N-1] ) / ( t[N] - t[N-1] );
	  func = f[N] + (table->slope[0]) * ( x[0] - t[N] ) ;
	}
	table->slope[1]=0.0;
	*sloper=table->slope[0];
      }
      break;

  case QUADRATIC:      /* quadratic lagrangian interpolation scheme */

      if(table->columns == 3)
      {
	for( i=0; i < N ; i+=2)
        {
          cee = (x[0]-t[i])/(t[i+2] - t[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1]+f[N+3+i]*phi[2];
	    break;
	  }
        }
        table->slope[2] = 0.0;
      }
      else
      {
        /* search through points for minimum distance*/
        dist_min = dist = sqrt(SQUARE(x[basis]-t[0])+SQUARE(x[ordinate]-f[0]));
	for( i=1; i < N ; i++)
        {
          dist = sqrt(SQUARE(x[basis]-t[i])+SQUARE(x[ordinate]-f[i]));
          if(dist < dist_min)
              { i_min = i;  dist_min = dist; }
        }
        elem = i_min/2-1; cee = 0.;
        for(iter=0 ; iter<10 ; iter++)
          {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
	    func = f[2*elem]*phi[0]+f[2*elem+1]*phi[1]+f[2*elem+2]*phi[2];
	    phic[0]=4.*cee-3.;
	    phic[1]=-8.*cee+4.;
	    phic[2]=4.*cee-1.;
            delta = t[2*elem+2] - t[2*elem];
	    *sloper = (f[2*elem]*phic[0]+f[2*elem+1]*phic[1]+f[2*elem+2]*phic[2]);
            point = t[2*elem]+cee*delta;
            update = ((point-x[basis])*delta+(func-x[ordinate])*(*sloper))
                             /(SQUARE(delta)+SQUARE(*sloper));
            cee -= update;
            if(fabs(update) < epstol) break;
            if(cee < 0.0 && elem > 0) 
                   { cee += 1.0;    elem -= 1;}
            if(cee > 1.0 && elem <  ((N-1)/2-1)) 
                   { cee -= 1.0;    elem += 1;}
	  }
        dist_min = sqrt(SQUARE(x[basis]-point)+SQUARE(x[ordinate]-func));
        table->slope[0] = dfunc_dx[basis] = (point-x[basis])/dist_min;
        dfunc_dx[ordinate] = (func-x[ordinate])/dist_min;
	table->slope[1]=point;
	table->slope[2]=func;
	*sloper=table->slope[0];
      }
      break;

  case QUAD_GP:      /* quadratic lagrangian interpolation scheme */

      if(table->columns == 3)
      {
	for( i=0; i < N ; i+=3)
        {
          xleft=(5.+sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.-sqrt(15.))/6.*t[i+2];
          xright=(5.-sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.+sqrt(15.))/6.*t[i+2];
          cee = (x[0]-xleft)/(xright-xleft);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-3))
	  {
	    phi[0]=(20.*cee*cee-2.*(sqrt(15.)+10.)*cee+sqrt(15.)+5.)/6.;
	    phi[1]=(-10.*cee*cee+10.*cee-1.)*2./3.;
	    phi[2]=(20.*cee*cee+2.*(sqrt(15.)-10.)*cee-sqrt(15.)+5.)/6.;
	    table->slope[0] = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    table->slope[1] = f[N+1+i]*phi[0]+f[N+2+i]*phi[1]+f[N+3+i]*phi[2];
	    break;
	  }
        }
        table->slope[2] = 0.0;
      }
      else
      {
	for( i=0; i < N ; i+=3)
        {
          xleft=(5.+sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.-sqrt(15.))/6.*t[i+2];
          xright=(5.-sqrt(15.))/6.*t[i]-2./3.*t[i+1]+(5.+sqrt(15.))/6.*t[i+2];
          cee = (x[0]-xleft)/(xright-xleft);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == N-3))
	  {
	    phi[0]=(20.*cee*cee-2.*(sqrt(15.)+10.)*cee+sqrt(15.)+5.)/6.;
	    phi[1]=(-10.*cee*cee+10.*cee-1.)*2./3.;
	    phi[2]=(20.*cee*cee+2.*(sqrt(15.)-10.)*cee-sqrt(15.)+5.)/6.;
	    func = f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2];
	    phi[0]=(20.*cee-sqrt(15.)+10.)/3.;
	    phi[1]=(-40.*cee+20.)/3.;
	    phi[2]=(20.*cee+sqrt(15.)-10.)/3.;
	    table->slope[0] = (f[i]*phi[0]+f[i+1]*phi[1]+f[i+2]*phi[2])
		/(xright-xleft);
	    break;
	  }
        }
	table->slope[1]=0.0;
	*sloper=table->slope[0];
      }
      break;

  case BIQUADRATIC:      /* biquadratic lagrangian interpolation scheme */

      ngrid1 = table->tablelength/table->ngrid;
      if(table->columns == 5)
      {
 	table->slope[0] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,f
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
 	table->slope[1] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,&f[N+1]
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
 	table->slope[2] = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,&f[2*N+2]
				     ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
      }
      else
      {
 	func = quad_isomap_invert(x[0],x[1],0,t,t2,NULL,f
 				  ,ngrid1,table->ngrid,1,2,2,dfunc_dx);
      }
      break;




  case BILINEAR:  /* BILINEAR Interpolation Scheme */
      /*Find Interval of Different Values of Abscissa #1*/
      iinter=0;
      for(i=0;i<N;i++)
      {
	if(table->t[i]!=table->t[i+1] && iinter==0)
	{ iinter=i+1;}
      }
      if(iinter == 1) 
      { 
	fprintf(stderr," MP Interpolate Error - Need more than 1 point per set");
	EH(-1, "Table interpolation not implemented");
      }

      istartx=iinter;

      for(i=iinter;t[i]<=x[0] && i<N-iinter+1;i=i+iinter)
      { istartx=i+iinter;}

      istarty=istartx;

      for(i=istartx+1;t2[i]<=x[1] && i<istartx+iinter-1;i++)
      { istarty=i;}

      y1=f[istarty];
      y2=f[istarty+1];
      y3=f[istarty+1-iinter];
      y4=f[istarty-iinter];

      tt=(x[1]-t2[istarty])/(t2[istarty+1]-t2[istarty]);
      uu=(x[0]-t[istarty])/(t[istarty+1-iinter]-t[istarty]);

      func=(1.-tt)*(1.-uu)*y1+tt*(1.-uu)*y2+tt*uu*y3+(1.-tt)*uu*y4;

      table->slope[1]= ( f[istarty+1] - f[istarty] ) / ( t2[istarty+1] - t2[istarty] );
      table->slope[0]= ( f[istarty+1] - f[istarty+1-iinter] ) / ( t[istarty+1] - t[istarty+1-iinter] );
      *sloper=table->slope[0];
	if(dfunc_dx != NULL)
	   {
		dfunc_dx[0] = table->slope[0];
		dfunc_dx[1] = table->slope[1];
	   }
      break;

   case TRILINEAR:      /* trilinear lagrangian interpolation scheme */
 
       ngrid1 = table->ngrid;
       ngrid2 = table->ngrid2/table->ngrid;
       ngrid3 = table->tablelength/table->ngrid2;
       func = quad_isomap_invert(x[0],x[1],x[2],t,t2,t3,f
 				  ,ngrid1,ngrid2,ngrid3,1,3,dfunc_dx);
       break;
 
   case TRIQUADRATIC:      /* triquadratic lagrangian interpolation scheme */
 
       ngrid1 = table->ngrid;
       ngrid2 = table->ngrid2/table->ngrid;
       ngrid3 = table->tablelength/table->ngrid2;
       func = quad_isomap_invert(x[0],x[1],x[2],t,t2,t3,f
 				  ,ngrid1,ngrid2,ngrid3,2,3,dfunc_dx);
       break;

  default:
      EH(-1, "Table interpolation order not implemented");
  }

  return(dist_min);
}
/*******************************************************************************
 *
 * continuous_surface_tension(): This function computes a tensor (csf) which is
 *                               a continuous, volumetric representation of the
 *           capillary stress that arises at a fluid/fluid interface (viz., the
 *           traction condition).  This is only used in the level set
 *           formulation.
 *
 * Input
 * -----
 *   Nothing.
 *
 * Ouput
 * -----
 *   csf	  = The "Continuous Surface Force" contribution to the fluid
 *                  stress tensor.
 *   d_csf_dF	  = Jacobian derivative of csf w.r.t. FILL.  If this is NULL,
 *                  then no Jacobian terms are computed.
 *
 *   NB: It is assumed that csf and d_csf_dF are zeroed.
 *
 * Returns
 * -------
 *   0 = Success.
 *  -1 = Failure.
 *
 ******************************************************************************/
int 
continuous_surface_tension(double st, double csf[DIM][DIM],
			   double d_csf_dF[DIM][DIM][MDE],
				double d_csf_dX[DIM][DIM][DIM][MDE] )
{
  int    a, b, c, j;
  int    status = 0;
  int    do_deriv;
  
  int var;
  
  var =  ls->var;
  
/*  double st = mp->surface_tension; */
  
  /* See if we're supposed to compute derivatives. */
  do_deriv = d_csf_dF != NULL;

  /* Fetch the level set interface functions. */
  load_lsi(ls->Length_Scale);
    
  /* If we're not near the zero level set, then csf is a zero tensor. */
  if ( ! lsi->near ) return(status);

  /* Calculate the CSF tensor. */
  for ( a=0; a < VIM; a++)
    {
      for ( b=0; b < VIM; b++)
	{
 	  csf[a][b] = st * lsi->delta *
	    ( (double)delta(a,b) - lsi->normal[a] * lsi->normal[b] );
	}
    }

  /* Bail out if we're done. */
  if ( ! do_deriv ) return(status);

#ifdef COUPLED_FILL
  load_lsi_derivs();

  /* Calculate the derivatives. */
  for ( j=0 ; j < ei->dof[var] ; j++ )
  {
#ifdef DO_NO_UNROLL
      for ( a=0 ; a < VIM ; a++ )
	  {
		  for ( b=0; b < VIM; b++ )
		  {
			  d_csf_dF[a][b][j] = st * lsi->d_delta_dF[j] *
			  ( (double)delta(a,b) - lsi->normal[a] * lsi->normal[b] )
			  - st * lsi->delta * lsi->d_normal_dF[a][j] * lsi->normal[b]
			  - st * lsi->delta * lsi->normal[a] * lsi->d_normal_dF[b][j];
		  }
	  }
#else
	/*(0,0) */
		  d_csf_dF[0][0][j] = st * lsi->d_delta_dF[j] *
		  ( (double)delta(0,0) - lsi->normal[0] * lsi->normal[0] )
		  - st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[0]
		  - st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[0][j];
		  
	/* (1,1) */
		  d_csf_dF[1][1][j] = st * lsi->d_delta_dF[j] *
			  ( (double)delta(1,1) - lsi->normal[1] * lsi->normal[1] )
			  - st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[1]
			  - st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[1][j];
	/* (0,1) */
		  d_csf_dF[0][1][j] = st * lsi->d_delta_dF[j] *
			  ( (double)delta(0,1) - lsi->normal[0] * lsi->normal[1] )
			  - st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[1]
			  - st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[1][j];
	/* (1,0) */
		  d_csf_dF[1][0][j] = st * lsi->d_delta_dF[j] *
			  ( (double)delta(1,0) - lsi->normal[1] * lsi->normal[0] )
			  - st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[0]
			  - st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[0][j];
		  if( VIM == 3 )
		  {
	/* (2,2) */
			  d_csf_dF[2][2][j] = st * lsi->d_delta_dF[j] *
			  ( (double)delta(2,2) - lsi->normal[2] * lsi->normal[2] )
			  - st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[2]
			  - st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[2][j];
	/* (2,1) */		  
			  d_csf_dF[2][1][j] = st * lsi->d_delta_dF[j] *
				  ( (double)delta(2,1) - lsi->normal[2] * lsi->normal[1] )
				  - st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[1]
				  - st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[1][j];
	/* (2,0) */		  
			  d_csf_dF[2][0][j] = st * lsi->d_delta_dF[j] *
				  ( (double)delta(2,0) - lsi->normal[2] * lsi->normal[0] )
				  - st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[0]
				  - st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[0][j];
	/* (1,2) */
			  d_csf_dF[1][2][j] = st * lsi->d_delta_dF[j] *
				  ( (double)delta(1,2) - lsi->normal[1] * lsi->normal[2] )
				  - st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[2]
				  - st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[2][j];
	/* (0,2)  */
			  d_csf_dF[0][2][j] = st * lsi->d_delta_dF[j] *
				  ( (double)delta(0,2) - lsi->normal[0] * lsi->normal[2] )
				  - st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[2]
				  - st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[2][j];
		  }
	  
#endif
  }
  
#endif /* COUPLED_FILL */
  if( pd->v[MESH_DISPLACEMENT1] )
  {
	  for ( a=0; a<VIM; a++)
	  {
		  for( b=0; b<VIM; b++)
		  {
			  for( c=0; c<VIM; c++)
			  {
				  var = MESH_DISPLACEMENT1 + c;
				  
				  for( j=0; j<ei->dof[var]; j++)
				  {
					d_csf_dX[a][b][c][j] = st* lsi->delta * ( - lsi->d_normal_dmesh[a][c][j]*lsi->normal[b]
															   -lsi->normal[a]*lsi->d_normal_dmesh[b][c][j] );
				  }
			  }
		  }
	  }
  }
  return(status);	
}

/*      routine for inverting isoparametric map
        finds cee,eta coordinates given the x,y coords
        and interpolates function

	Points on a rectangular mapped grid are assumed.
	It should work for either linear or quadratic interpolation
 
 	Try to include 3D bricks also - 9/15/2004
*/
double
quad_isomap_invert( const double x1,
		    const double y1,
 		    const double z1,
		    const double x[],
		    const double y[],
 		    const double z[],
		    const double p[],
		    const int ngrid1,
		    const int ngrid2,
 		    const int ngrid3,
 		    const int elem_order,
 		    const int dim,
		    double dfunc_dx[])
{
        int i,j,k,l,iter;
        double pt[5] = {0.5,1.0,0.0,0.5,1.0};
        double phi[27],phic[27],phie[27],phig[27];
        double coord[3][27],xpt[3],pp[27];
        double jac[3][3],detjt,detjti = 0.0;
        static double xi[3] = {0.5,0.5,0.5};
        static int nell=0;
        int itp[27];
        int nell_xi[3],ne_xi[3];
        double dxi[3],eps,pvalue, pc, pe, pg;
        double epstol=1.0e-4;
  	double loc_tol=0.0001;
 	int max_its = 20;
	int elem_nodes;
        double bfi,bfj,dfi,dfj,bfl,dfl;

	/* Initialize dxi */
	dxi[0] = 0;
	dxi[1] = 0;
	dxi[2] = 0;

 	elem_nodes = pow(elem_order+1,dim);
 	switch (dim)
 		{
 		case 2:
         		ne_xi[0]=(ngrid1-1)/elem_order;
         		ne_xi[1]=(ngrid2-1)/elem_order;
         		nell_xi[0] = nell/ne_xi[1];
         		nell_xi[1] = nell%ne_xi[1];
 			break;
 		case 3:
         		ne_xi[0]=(ngrid1-1)/elem_order;
         		ne_xi[1]=(ngrid2-1)/elem_order;
         		ne_xi[2]=(ngrid3-1)/elem_order;
         		nell_xi[2] = nell/(ne_xi[0]*ne_xi[1]);
         		nell_xi[0] = nell%ne_xi[0];
         		nell_xi[1] = (nell-nell_xi[2]*ne_xi[0]*ne_xi[1])
 					/ne_xi[0];
 			break;
 		}

/*      iterate to find isoparametric coordinates  */

        iter=0;eps=10.*epstol;
        while ( iter < max_its && eps > epstol)
        {
        iter++;

/**     evaluate basis functions  **/

	switch (elem_order)
	{
	case 2:
 		switch(dim)
 		{
 		case 2:
         	for( i=0 ; i<3 ; i++ )
                 	{
                 	k=3*i;
                 	bfi = (xi[0]-pt[i])*(xi[0]-pt[i+1])/
                        	      ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 	dfi = (2.*xi[0]-pt[i]-pt[i+1])/
                        	      ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 	for (j=0 ; j<3 ; j++)
                         	{
                         	bfj = (xi[1]-pt[j])*(xi[1]-pt[j+1])/
                                     ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
                         	dfj = (2.*xi[1]-pt[j]-pt[j+1])/
                                     ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
                         	phi[k+j]=bfi*bfj;
                         	phic[k+j]=dfi*bfj;
                         	phie[k+j]=bfi*dfj;
                         	}
                 	}
 			break;
 		case 3:
         	for( i=0 ; i<3 ; i++ )
                 	{
                 	k=3*i;
                 	bfi = (xi[0]-pt[i])*(xi[0]-pt[i+1])/
                        	      ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 	dfi = (2.*xi[0]-pt[i]-pt[i+1])/
                        	      ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 	for (j=0 ; j<3 ; j++)
                         	{
                         	bfj = (xi[1]-pt[j])*(xi[1]-pt[j+1])/
                                     ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
                         	dfj = (2.*xi[1]-pt[j]-pt[j+1])/
                                     ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
 				for (l=0 ; l<3 ; l++)
 					{
                         		bfl = (xi[1]-pt[l])*(xi[1]-pt[l+1])/
                                     	   ((pt[l+2]-pt[l])*(pt[l+2]-pt[l+1]));
                         		dfl = (2.*xi[1]-pt[l]-pt[l+1])/
                                     	   ((pt[l+2]-pt[l])*(pt[l+2]-pt[l+1]));
                         		phi[k+j+9*l]=bfi*bfj*bfl;
                         		phic[k+j+9*l]=dfi*bfj*bfl;
                         		phie[k+j+9*l]=bfi*dfj*bfl;
                         		phig[k+j+9*l]=bfi*bfj*dfl;
 					}
                         	}
                 	}
 			break;
                }
	break;

	case 1:
 		switch(dim)
 		{
 		case 2:
 			phi[0]=(1.-xi[0])*(1.-xi[1]);
 			phi[1]=(1.-xi[0])*xi[1];
 			phi[2]=xi[0]*(1.-xi[1]);
 			phi[3]=xi[0]*xi[1];
 			phic[0]=xi[1]-1.;
 			phic[1]=-xi[1];
 			phic[2]=1.-xi[1];
 			phic[3]=xi[1];
 			phie[0]=xi[0]-1.;
 			phie[1]=1.-xi[0];
 			phie[2]=-xi[0];
 			phie[3]=xi[0];
 			break;
 		case 3:
 			phi[0]=(1.-xi[0])*(1.-xi[1])*(1.-xi[2]);
 			phi[1]=(1.-xi[0])*xi[1]*(1.-xi[2]);
 			phi[2]=xi[0]*(1.-xi[1])*(1.-xi[2]);
 			phi[3]=xi[0]*xi[1]*(1.-xi[2]);
 			phi[4]=(1.-xi[0])*(1.-xi[1])*xi[2];
 			phi[5]=(1.-xi[0])*xi[1]*xi[2];
 			phi[6]=xi[0]*(1.-xi[1])*xi[2];
 			phi[7]=xi[0]*xi[1]*xi[2];
 			phic[0]=-(1.-xi[1])*(1.-xi[2]);
 			phic[1]=-xi[1]*(1.-xi[2]);
 			phic[2]=(1.-xi[1])*(1.-xi[2]);
 			phic[3]=xi[1]*(1.-xi[2]);
 			phic[4]=(xi[1]-1.)*xi[2];
 			phic[5]=-xi[1]*xi[2];
 			phic[6]=(1.-xi[1])*xi[2];
 			phic[7]=xi[1]*xi[2];
 			phie[0]=-(1.-xi[0])*(1.-xi[2]);
 			phie[1]=(1.-xi[0])*(1.-xi[2]);
 			phie[2]=-xi[0]*(1.-xi[2]);
 			phie[3]=xi[0]*(1.-xi[2]);
 			phie[4]=-(1.-xi[0])*xi[2];
 			phie[5]=(1.-xi[0])*xi[2];
 			phie[6]=-xi[0]*xi[2];
 			phie[7]=xi[0]*xi[2];
 			phig[0]=-(1.-xi[0])*(1.-xi[1]);
 			phig[1]=-(1.-xi[0])*xi[1];
 			phig[2]=-xi[0]*(1.-xi[1]);
 			phig[3]=-xi[0]*xi[1];
 			phig[4]=(1.-xi[0])*(1.-xi[1]);
 			phig[5]=(1.-xi[0])*xi[1];
 			phig[6]=xi[0]*(1.-xi[1]);
 			phig[7]=xi[0]*xi[1];
 			break;
 		}
	break;

	default:
	fprintf(stderr,"\n Unknown element order - quad_isomap\n");
	EH(-1,"Fatal Error");
	}

/**  gather local node numbers  **/

 	switch(dim)
 	{
 	case 2:
         	itp[0] = elem_order*(elem_order*ne_xi[1]+1)*nell_xi[0] 
			+ elem_order*nell_xi[1];
         	itp[1] = itp[0] + 1;
 		if(elem_order == 2)itp[2] = itp[1] + 1;
         	for( i=1; i<=elem_order ; i++)
                	    {
                      for( j=0; j<=elem_order; j++)    
 		 	{
 			 itp[(elem_order+1)*i+j] = itp[j]+
 					(elem_order*ne_xi[1]+1)*i;
 			}
                      }
 		break;
 	case 3:
         	itp[0] = nell_xi[2]*elem_order*(elem_order*ne_xi[0]+1)
 					*(elem_order*ne_xi[1]+1)
 			+ elem_order*(elem_order*ne_xi[0]+1)*nell_xi[1] 
 			+ elem_order*nell_xi[0];
         	itp[1] = itp[0] + elem_order*(elem_order*ne_xi[0]+1);
 		if(elem_order == 2)
         		itp[2] = itp[1] + elem_order*(elem_order*ne_xi[0]+1);
         	for( i=0; i<=elem_order ; i++)
                	    {
                      for( j=0; j<=elem_order; j++)    
 		 	{
 			for( k=0 ; k<=elem_order ; k++)
 			   {
 			    itp[(elem_order+1)*i+j+SQUARE(elem_order+1)*k] = 
 				itp[j] + i + k*(elem_order*ne_xi[0]+1)
 					 *(elem_order*ne_xi[1]+1);
 			   }
 			}
                     }
 		break;
 	}
 
/**   gather global variables into local arrays   */
 
 	switch(dim)
 	{
 	case 3:
         	for( i=0; i<elem_nodes; i++)
                		{
                 	coord[2][i] = z[itp[i]];
                 	}
 	case 2:
         	for( i=0; i<elem_nodes; i++)
                		{
                 	coord[0][i] = x[itp[i]];
                 	coord[1][i] = y[itp[i]];
                 	}
 		break;
 	}
 
 	memset( xpt, 0, sizeof(dbl)*3);
 	memset( jac, 0, sizeof(dbl)*9);
 
 	switch(dim)
 	{
 	case 3:
         	for( i=0; i<elem_nodes ; i++)
                 	{
                 	jac[0][2] += coord[0][i]*phig[i];
                 	jac[1][2] += coord[1][i]*phig[i];
                 	jac[2][0] += coord[2][i]*phic[i];
                 	jac[2][1] += coord[2][i]*phie[i];
                 	jac[2][2] += coord[2][i]*phig[i];
                 	xpt[2] += coord[2][i]*phi[i];
                 	}
 	case 2:
         	for( i=0; i<elem_nodes ; i++)
                 	{
                 	jac[0][0] += coord[0][i]*phic[i];
                 	jac[0][1] += coord[0][i]*phie[i];
                 	jac[1][0] += coord[1][i]*phic[i];
                 	jac[1][1] += coord[1][i]*phie[i];
                 	xpt[0] += coord[0][i]*phi[i];
                 	xpt[1] += coord[1][i]*phi[i];
                 	}
 		break;
 	}
 
 
 	switch(dim)
 	{
 	case 2:
        		detjt=jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0];
        		detjti=1./detjt;
         	dxi[0]=(-jac[1][1]*(xpt[0]-x1)+jac[0][1]*(xpt[1]-y1))*detjti;
         	dxi[1]=(jac[1][0]*(xpt[0]-x1)-jac[0][0]*(xpt[1]-y1))*detjti;
 		break;
 	case 3:
         	detjt=jac[0][0]*(jac[1][1]*jac[2][2]-jac[2][1]*jac[1][2])
               		- jac[1][0]*(jac[0][1]*jac[2][2]-jac[2][1]*jac[0][2])
               		+ jac[2][0]*(jac[1][2]*jac[0][1]-jac[1][1]*jac[0][2]);
        		detjti=1./detjt;
         	dxi[0]=detjti*(
 			-(xpt[0]-x1)*(jac[2][2]*jac[1][1]-jac[2][1]*jac[1][2])
                 	-(xpt[1]-y1)*(jac[2][1]*jac[0][2]-jac[2][2]*jac[0][1])
                 	-(xpt[2]-z1)*(jac[1][2]*jac[0][1]-jac[1][1]*jac[0][2]));
         	dxi[1]=detjti*(
 			-(xpt[0]-x1)*(jac[2][0]*jac[1][2]-jac[2][2]*jac[1][0])
                		-(xpt[1]-y1)*(jac[2][2]*jac[0][0]-jac[2][0]*jac[0][2])
                 	-(xpt[2]-z1)*(jac[1][0]*jac[0][2]-jac[1][2]*jac[0][0]));
         	dxi[2]=detjti*(
 			-(xpt[0]-x1)*(jac[2][1]*jac[1][0]-jac[2][0]*jac[1][1])
                 	-(xpt[1]-y1)*(jac[2][0]*jac[0][1]-jac[2][1]*jac[0][0])
                 	-(xpt[2]-z1)*(jac[1][1]*jac[0][0]-jac[1][0]*jac[0][1]));
 		break;
 	}

/*newton updates  */

        eps=0.0;
        for(i=0;i<dim;i++)
                {
                xi[i] += dxi[i];
                eps += dxi[i]*dxi[i];
                }
 	eps=sqrt(eps)/((double)dim);

/*check to see if crossed element boundaries   */

        for (i=0; i<dim; i++)
                {
 		double fnum;
 		int dell;
                if(xi[i] < -loc_tol && nell_xi[i] > 0)
                        {
 			fnum = floor(xi[i]);
 			dell = MAX(0,nell_xi[i]+fnum-1)-nell_xi[i];
 			if(iter > max_its/2) dell = -1;
                        nell_xi[i]  += dell;
                        xi[i] -= ((double)dell);
 			if(dell != 0) eps=10.0*epstol;
                        }
                else if(xi[i] > 1.0+loc_tol && nell_xi[i] < ne_xi[i]-1)
                        {
 			fnum = floor(xi[i]);
 			dell = MIN(ne_xi[i]-1, nell_xi[i]+fnum)-nell_xi[i];
 			if(iter > max_its/2) dell = 1;
                        nell_xi[i]  += dell;
                        xi[i] -= ((double)dell);
 			if(dell != 0) eps=10.0*epstol;
                        }
                }
        }

/*check to see if iteration was successful      */

        if(eps > epstol)
                {
            fprintf(stderr," pquad iteration failed %d %g %g %g\n",iter,
 			xi[0],xi[1],xi[2]);
            fprintf(stderr," point %g %g %g\n",x1,y1,z1);
		EH(-1,"Fatal Error");
                }

/*evaluate function             */

	switch (elem_order)
	{
	case 2:
 		switch(dim)
 		{
 		case 2:
         		for( i=0 ; i<3 ; i++ )
                 		{
                 		k=3*i;
                 		bfi = (xi[0]-pt[i])*(xi[0]-pt[i+1])/
                         	    ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 		for (j=0 ; j<3 ; j++)
                         		{
                         		bfj = (xi[1]-pt[j])*(xi[1]-pt[j+1])/
                                 	    ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
                         		phi[k+j]=bfi*bfj;
                         		}
                 		}
 			break;
 		case 3:
         	for( i=0 ; i<3 ; i++ )
                 	{
                 	k=3*i;
                 	bfi = (xi[0]-pt[i])*(xi[0]-pt[i+1])/
                        	      ((pt[i+2]-pt[i])*(pt[i+2]-pt[i+1]));
                 	for (j=0 ; j<3 ; j++)
                         	{
                         	bfj = (xi[1]-pt[j])*(xi[1]-pt[j+1])/
                                     ((pt[j+2]-pt[j])*(pt[j+2]-pt[j+1]));
 				for (l=0 ; l<3 ; l++)
 					{
                         		bfl = (xi[1]-pt[l])*(xi[1]-pt[l+1])/
                                     	   ((pt[l+2]-pt[l])*(pt[l+2]-pt[l+1]));
                         		phi[k+j+9*l]=bfi*bfj*bfl;
 					}
                         	}
                 	}
 			break;
 		}
	break;

	case 1:
 		switch(dim)
 		{
 		case 2:
 			phi[0]=(1.-xi[0])*(1.-xi[1]);
 			phi[1]=(1.-xi[0])*xi[1];
 			phi[2]=xi[0]*(1.-xi[1]);
 			phi[3]=xi[0]*xi[1];
		    if(dfunc_dx != NULL)
			{
 			phic[0]=xi[1]-1.;
 			phic[1]=-xi[1];
 			phic[2]=1.-xi[1];
 			phic[3]=xi[1];
 			phie[0]=xi[0]-1.;
 			phie[1]=1.-xi[0];
 			phie[2]=-xi[0];
 			phie[3]=xi[0];
			}
 			break;
 		case 3:
 			phi[0]=(1.-xi[0])*(1.-xi[1])*(1.-xi[2]);
 			phi[1]=(1.-xi[0])*xi[1]*(1.-xi[2]);
 			phi[2]=xi[0]*(1.-xi[1])*(1.-xi[2]);
 			phi[3]=xi[0]*xi[1]*(1.-xi[2]);
 			phi[4]=(1.-xi[0])*(1.-xi[1])*xi[2];
 			phi[5]=(1.-xi[0])*xi[1]*xi[2];
 			phi[6]=xi[0]*(1.-xi[1])*xi[2];
 			phi[7]=xi[0]*xi[1]*xi[2];
		    if(dfunc_dx != NULL)
			{
 			phic[0]=-(1.-xi[1])*(1.-xi[2]);
 			phic[1]=-xi[1]*(1.-xi[2]);
 			phic[2]=(1.-xi[1])*(1.-xi[2]);
 			phic[3]=xi[1]*(1.-xi[2]);
 			phic[4]=(xi[1]-1.)*xi[2];
 			phic[5]=-xi[1]*xi[2];
 			phic[6]=(1.-xi[1])*xi[2];
 			phic[7]=xi[1]*xi[2];
 			phie[0]=-(1.-xi[0])*(1.-xi[2]);
 			phie[1]=(1.-xi[0])*(1.-xi[2]);
 			phie[2]=-xi[0]*(1.-xi[2]);
 			phie[3]=xi[0]*(1.-xi[2]);
 			phie[4]=-(1.-xi[0])*xi[2];
 			phie[5]=(1.-xi[0])*xi[2];
 			phie[6]=-xi[0]*xi[2];
 			phie[7]=xi[0]*xi[2];
 			phig[0]=-(1.-xi[0])*(1.-xi[1]);
 			phig[1]=-(1.-xi[0])*xi[1];
 			phig[2]=-xi[0]*(1.-xi[1]);
 			phig[3]=-xi[0]*xi[1];
 			phig[4]=(1.-xi[0])*(1.-xi[1]);
 			phig[5]=(1.-xi[0])*xi[1];
 			phig[6]=xi[0]*(1.-xi[1]);
 			phig[7]=xi[0]*xi[1];
			}
 			break;
 		}
	break;

	default:
	fprintf(stderr,"\n Unknown element order - quad_isomap\n");
	EH(-1,"Fatal Error");
	}

        for( i=0; i<elem_nodes; i++)
                {
                pp[i] = p[itp[i]];
                }
        pvalue=0.0; pc=0.; pe=0.;  pg=0.;
        for( i=0; i<elem_nodes ; i++)
                {
                pvalue += pp[i]*phi[i];
		    if(dfunc_dx != NULL)
			{
 				switch(dim)
 				{
 				case 3:
                			pg += pp[i]*phig[i];
 				case 2:
                			pc += pp[i]*phic[i];
                			pe += pp[i]*phie[i];
					break;
                		}
                     	}
		}
	if(dfunc_dx != NULL)
	{
 	switch(dim)
 		{
 		case 2:
         		dfunc_dx[0]=(jac[1][1]*pc-jac[0][1]*pe)*detjti;
         		dfunc_dx[1]=(-jac[1][0]*pc+jac[0][0]*pe)*detjti;
			break;
 		case 3:
         	dfunc_dx[0]=detjti*(
 			pc*(jac[2][2]*jac[1][1]-jac[2][1]*jac[1][2])
                 	+ pe*(jac[2][1]*jac[0][2]-jac[2][2]*jac[0][1])
                 	+ pg*(jac[1][2]*jac[0][1]-jac[1][1]*jac[0][2]));
         	dfunc_dx[1]=detjti*(
 			pc*(jac[2][0]*jac[1][2]-jac[2][2]*jac[1][0])
                	+ pe*(jac[2][2]*jac[0][0]-jac[2][0]*jac[0][2])
                 	+ pg*(jac[1][0]*jac[0][2]-jac[1][2]*jac[0][0]));
         	dfunc_dx[2]=detjti*(
 			pc*(jac[2][1]*jac[1][0]-jac[2][0]*jac[1][1])
                 	+ pe*(jac[2][0]*jac[0][1]-jac[2][1]*jac[0][0])
                 	+ pg*(jac[1][1]*jac[0][0]-jac[1][0]*jac[0][1]));
			break;
 		}
	}
        return(pvalue);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

static double
scalar_fv_fill_altmatrl(double **base, int lvdesc, int num_dofs,
			int var_type)

    /**************************************************************************
     *   scalar_fv_fill_altmatr:
     *
     *      Kernal operation for calculating the value of a
     *       scalar field variable at the quadrature point from
     *      the basis functions.
     *      Note for the optimal results, this function should
     *      be inlined during the compile.
     **************************************************************************/
{
  int ln, idof, lvdof, lvdof_active;
  double *phi_ptr = bf[var_type]->phi, sum = 0.0;
  for (idof = 0; idof < num_dofs; idof++) {
    /*
     *  Find the local variable dof for the degree of freedom
     *  that actually belongs to another material than the one
     *  we are currently in.
     */
    lvdof = ei->Lvdesc_to_lvdof[lvdesc][idof];
    /*
     *  Find the local node number of this degree of freedom
     */
    ln = ei->dof_list[var_type][lvdof];
    /*
     *  Find the local variable type degree of freedom 
     *  active in this element for this local node. This
     *  degree of freedom belongs to this material, and thus
     *  has a nonzero basis function associated with it.
     */
    lvdof_active = ei->ln_to_dof[var_type][ln];
    /*
     *  Add it to the sum to get the value of the variable.
     */
    sum += *(base[lvdof]) * phi_ptr[lvdof_active];
  }
  return sum;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double
scalar_fv_fill_adjmatrl(double **base, int lvdesc, int num_dofs,
			int var_type)

    /**************************************************************************
     *   scalar_fv_fill_adjmatr:
     *
     *      Identical to scalar_fv_fill_altmatrl but with phi_ptr array indexed
     *      by local node number.  This probably the way it should be done for interpolation
     *      orders that assign only one dof per node.  This mod was necessary to get it
     *      to work for look across BC's
     *      tabaer March 2004
     *     
     **************************************************************************/
{
  int ln, idof, lvdof;
  double *phi_ptr = bf[var_type]->phi, sum = 0.0;
  for (idof = 0; idof < num_dofs; idof++) {
    /*
     *  Find the local variable dof for the degree of freedom
     *  that actually belongs to another material than the one
     *  we are currently in.
     */
    lvdof = ei->Lvdesc_to_lvdof[lvdesc][idof];
    /*
     *  Find the local node number of this degree of freedom
     */
    ln = ei->dof_list[var_type][lvdof];

    /*
     *  Add it to the sum to get the value of the variable.
     */
    sum += *(base[lvdof]) * phi_ptr[ln];
  }
  return sum;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void
load_matrl_statevector(MATRL_PROP_STRUCT *mp_local)

/**************************************************************************
 * load_matrl_statevector():
 *
 *  load up values of all relevant state variables at the
 *  current surface quad pt for materials that are not associated with the
 *  current element.
 *
 * input:
 * ------
 *	We assume that the state vector for the material of the
 *      current element has already been filled in. We use its values as
 *      defaults.
 *
 * output: 
 * -------
 *       mp_local->StateVector[] is filled in as well, for 
 *                 pertinent entries that make up the specification of 2q
 *                 the state of the material.
 *
 * HKM -> NOTE , probably not complete in terms of covering all 
 *        needed state variable entries.
 *************************************************************************/
{
  int lvdesc, var_type, num_dofs, k, do_last_species = FALSE;
  double *sv = mp_local->StateVector;  
  double *sv_mp = mp->StateVector, sp_sum = 0.0, rho;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int matID = mp_local->MatID;

  sv[TEMPERATURE] = sv_mp[TEMPERATURE];
  sv[FILL] = sv_mp[FILL];
  sv[VOLTAGE] = sv_mp[VOLTAGE];
  sv[PRESSURE] = sv_mp[PRESSURE];
  for (k = 0; k < mp->Num_Species; k++) {
    sv[SPECIES_UNK_0 + k] = sv_mp[SPECIES_UNK_0 + k];
  }

  for (lvdesc = 0; lvdesc < ei->Num_Lvdesc; lvdesc++) {
    if (ei->Lvdesc_to_MatID[lvdesc] == matID) {

      var_type = ei->Lvdesc_to_Var_Type[lvdesc];
      num_dofs = ei->Lvdesc_Numdof[lvdesc];
     
      if (var_type == MASS_FRACTION) {
        vd = ei->Lvdesc_vd_ptr[lvdesc];
        k = vd->Subvar_Index;
        sv[SPECIES_UNK_0 + k] = scalar_fv_fill_altmatrl(esp->c[k], lvdesc, 
							num_dofs, var_type);
	sp_sum += sv[SPECIES_UNK_0 + k];
	mp_local->StateVector_speciesVT = upd->Species_Var_Type;
        if (mp_local->Num_Species > mp_local->Num_Species_Eqn) {
	  do_last_species = TRUE;
	}
      } else if (var_type == TEMPERATURE) {
	sv[var_type] = scalar_fv_fill_altmatrl(esp->T, lvdesc, num_dofs,
					       var_type);
      } else if (var_type >= VELOCITY1 && var_type <= VELOCITY3) {
        k = var_type - VELOCITY1;
	sv[var_type] = scalar_fv_fill_altmatrl(esp->v[k], lvdesc, num_dofs,
					       var_type);
      } else if (var_type == PRESSURE) {
	sv[var_type] = scalar_fv_fill_altmatrl(esp->P, lvdesc, num_dofs,
					       var_type);
      } else if (var_type == FILL) {
	sv[var_type] = scalar_fv_fill_altmatrl(esp->F, lvdesc, num_dofs,
					       var_type);
      } else if (var_type == LIGHT_INTP) {
	sv[var_type] = scalar_fv_fill_altmatrl(esp->poynt[0], lvdesc, num_dofs,
					       var_type);
      } else if (var_type == LIGHT_INTM) {
	sv[var_type] = scalar_fv_fill_altmatrl(esp->poynt[1], lvdesc, num_dofs,
					       var_type);
      } else {
        EH(-1,"Unimplemented");
      }

    }
  }
  if (do_last_species) {
    var_type = SPECIES_UNK_0 + mp_local->Num_Species - 1;
    switch (mp_local->Species_Var_Type) {
    case SPECIES_UNDEFINED_FORM:
    case SPECIES_MOLE_FRACTION:
    case SPECIES_MASS_FRACTION:
    case SPECIES_VOL_FRACTION:
	sv[var_type] = 1.0 - sp_sum; 
	break;
    case SPECIES_CONCENTRATION:
	switch(mp_local->DensityModel) {
	case DENSITY_CONSTANT_LAST_CONC:
	    sv[var_type] = mp_local->u_density[0];
	    break;
	default:	    
	    rho = calc_concentration(mp_local, FALSE, NULL);
	    sv[var_type] = rho - sp_sum;
	    break;
	}
	break;
    case SPECIES_DENSITY:
      /* Note, this won't work for time dependent densities - RRR */
	rho = calc_density(mp_local, FALSE, NULL, 0.0);
	sv[var_type] = rho - sp_sum;
	break;
    case SPECIES_CAP_PRESSURE:
	break;
    default:
	break;
    }
  }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
 * These density models locally permit a time and spatially varying
 *  density.  Consequently, neither the gradient nor the time derivative 
 *  of the density term in the continuity equation is zero. These terms
 *  must then be included here as a "source" term via the expresion:
 *
 *   0 = del dot V + FoamVolumeSource
 *
 *  Therefore, this implies that
 *
 *     FoamVolumeSource = 1/rho * drho/dt + 1/rho * (v dot del rho)
 *
 */
double FoamVolumeSource(double time,
			double dt,
			double tt,
			double dFVS_dv[DIM][MDE], 
			double dFVS_dT[MDE],
			double dFVS_dx[DIM][MDE],
			double dFVS_dC[MAX_CONC][MDE],
			double dFVS_dF[MDE])
{
  double DT_Dt, vol, rho, rho2, T, Press;
  double deriv_x, deriv_T;
  double source = 0.0, source_a, source_b, source_c;
  int var, j, a, dim, w, species;
  double  d_rho_dt, drho_T;

  /* density derivatives */
  /*   DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /\* density dependence *\/ */
  /*   DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct; */

  dim = pd->Num_Dim;

  memset( dFVS_dv, 0 , sizeof(double)*DIM*MDE );
  memset( dFVS_dT, 0 , sizeof(double)*MDE );
  memset( dFVS_dx, 0 , sizeof(double)*DIM*MDE );
  memset( dFVS_dF, 0 , sizeof(double)*MDE );
  memset( dFVS_dC, 0 , sizeof(double)*MDE*MAX_CONC );


  if (mp->DensityModel == DENSITY_FOAM)
    {
      dbl x0, MW, Rgas, rho_epoxy, rho_fluor, Dx_Dt, drho_x; 
      species = (int) mp->u_density[0];  /* species number fluorinert */
      x0 = mp->u_density[1];         /* Initial fluorinert mass fraction */
      Rgas =  mp->u_density[2];      /* Gas constant in appropriate units */
      MW = mp->u_density[3];
      rho_epoxy = mp->u_density[4];  /* Density of epoxy resin */
      rho_fluor =  mp->u_density[5]; /* Density of liquid fluorinert */

      if (fv->c[species] > 0.)
	vol = fv->c[species];
      else
	vol = 0.;
      if (vol > x0) vol = x0;
      
      T = fv->T;
      Press = upd->Pressure_Datum;
      
      rho = 1./((x0-vol)*Rgas*T/(Press*MW)+(1.0-x0)/rho_epoxy+vol/rho_fluor);
      rho2 = rho*rho;
      /*   rho = density(d_rho); */

      Dx_Dt = fv_dot->c[species];
      for( a=0; a<dim; a++)
	{
	  Dx_Dt += fv->v[a]*fv->grad_c[species][a];
	}
      
      DT_Dt = fv_dot->T;
      for( a=0; a<dim; a++)
	{
	  DT_Dt += fv->v[a]*fv->grad_T[a];
	}
      deriv_x = (-Rgas*T/(Press*MW)+1./rho_fluor);
      deriv_T =  (x0-vol)*Rgas/(Press*MW);
      drho_x =  -rho2*deriv_x;
      drho_T =  -rho2*deriv_T;
      if((vol > 0.)&& (vol < x0) )  
	source = -rho*(deriv_x*Dx_Dt +deriv_T*DT_Dt);
      else
	source = 0.;
      
      
      if (af->Assemble_Jacobian && ((vol > 0.)&& (vol < x0) ))
	{
	  
	  var = TEMPERATURE;
	  if (pd->v[var])
	    {
	      for( j=0; j<ei->dof[var]; j++)
		{
		  source_a = bf[var]->phi[j] * (1 + 2. * tt)/ dt;
		  for( a=0; a<dim; a++)
		    {
		      source_a += fv->v[a] * bf[var]->grad_phi[j][a];
		    }
		  source_a *= -rho*deriv_T;                                           /* derivative of convective terms */
		  source_b =  -drho_T*(deriv_x*Dx_Dt +deriv_T*DT_Dt)*bf[var]->phi[j]; /* derivative of density */
		  source_c =  -rho*Rgas/(Press*MW)*Dx_Dt*bf[var]->phi[j];             /* derivative of deriv_x wrt T */
		  dFVS_dT[j] = source_a + source_b + source_c;	
		}
	    }
	  
	  var = VELOCITY1;
	  if (pd->v[var])
	    {
	      for( a=0; a<dim; a++)
		{
		  var = VELOCITY1+a;
		  for( j=0; pd->v[var] && j<ei->dof[var]; j++)
		    {
		      dFVS_dv[a][j] = -rho*(deriv_x* bf[var]->phi[j]*fv->grad_c[species][a]
					    + deriv_T* bf[var]->phi[j]*fv->grad_T[a]);
		    }
		}
	    }
	  
	  var = MASS_FRACTION;
	  if (pd->v[var])
	    {
	      for( j=0; j<ei->dof[var]; j++)
		{
		  source_a = bf[var]->phi[j] * (1 + 2. * tt)/ dt;
		  for( a=0; a<dim; a++)
		    {
		      source_a += fv->v[a] * bf[var]->grad_phi[j][a];
		    }
		  source_a *= -rho*deriv_x;                                            /* derivative of convective terms */
		  source_b = -drho_x*(deriv_x*Dx_Dt +deriv_T*DT_Dt)*bf[var]->phi[j];  /* derivative of density */
		  source_c = -rho*Rgas/(Press*MW)*DT_Dt*bf[var]->phi[j];              /* derivative of deriv_T wrt x */
		  if((vol > 0.)&& (vol < x0) )
		    dFVS_dC[species][j] = source_a + source_b + source_c;
		}
	    }
	  
	}
    }
  else if (mp->DensityModel == DENSITY_FOAM_CONC)
    {
      double Rgas, MW_f, MW_a, rho_epoxy, rho_fluor, T;
      int species_l, species_v, species_a;
      double Dcl_Dt,  Dcv_Dt, Dca_Dt, DT_Dt;
      /* density derivatives */
      double drho_T, rho_v_inv, d_rho_v_inv_dT, rho_a_inv, d_rho_a_inv_dT;
      double drho_c_v, drho_c_a, drho_c_l, d2rho_T_c_v, d2rho_T_c_a; 
      double source_a_v, source_a_a, source_a_l, source_b_v, source_b_a, source_b_l, source_c_v, source_c_a;

      species_l = (int) mp->u_density[0];  /* species number fluorinert liquid */
      species_v = (int) mp->u_density[1];  /* species number fluorinert vapor */
      species_a = (int) mp->u_density[2];  /* species number air vapor */
      Rgas =  mp->u_density[3];      /* Gas constant in appropriate units */
      MW_f = mp->u_density[4];         /* molecular weight fluorinert */
      MW_a = mp->u_density[5];         /* molecular weight air */
      rho_epoxy = mp->u_density[6];  /* Density of epoxy resin */
      rho_fluor =  mp->u_density[7]; /* Density of liquid fluorinert */


      T = fv->T;
      Press = upd->Pressure_Datum;
      rho_v_inv = Rgas*T/(Press*MW_f);
      d_rho_v_inv_dT = Rgas/(Press*MW_f);
      rho_a_inv = Rgas*T/(Press*MW_a);
      d_rho_a_inv_dT = Rgas/(Press*MW_a);
	
      rho = rho_epoxy+ fv->c[species_v]*(1.-rho_epoxy*rho_v_inv) 
	             + fv->c[species_a]*(1.-rho_epoxy*rho_a_inv)
	             + fv->c[species_l]*(1.-rho_epoxy/rho_fluor);

      rho2 = rho * rho;
   
      drho_c_v =1.-rho_epoxy*rho_v_inv;
      drho_c_a =1.-rho_epoxy*rho_a_inv;
      drho_c_l =1.-rho_epoxy/rho_fluor;

      drho_T = -fv->c[species_v]*rho_epoxy*d_rho_v_inv_dT 
	- fv->c[species_a]*rho_epoxy*d_rho_a_inv_dT ;

      d2rho_T_c_v = -rho_epoxy*d_rho_v_inv_dT;
      d2rho_T_c_a = -rho_epoxy*d_rho_a_inv_dT;

      Dcl_Dt = fv_dot->c[species_l];
      Dcv_Dt = fv_dot->c[species_v];
      Dca_Dt = fv_dot->c[species_a];
      for( a=0; a<dim; a++)
	{
	  Dcl_Dt += fv->v[a]*fv->grad_c[species_l][a];
	  Dcv_Dt += fv->v[a]*fv->grad_c[species_v][a];
	  Dca_Dt += fv->v[a]*fv->grad_c[species_a][a];
	}
      
      DT_Dt = fv_dot->T;
      for( a=0; a<dim; a++)
	{
	  DT_Dt += fv->v[a]*fv->grad_T[a];
	}

      source = (drho_c_v*Dcv_Dt + drho_c_a*Dca_Dt + drho_c_l*Dcl_Dt + drho_T*DT_Dt)/rho;
      
      if (af->Assemble_Jacobian)
	{
	  var = TEMPERATURE;
	  if(pd->v[var] )
	    {
	      for( j=0; j<ei->dof[var]; j++)
		{
		  source_a = bf[var]->phi[j] * (1 + 2. * tt)/ dt;
		  for( a=0; a<dim; a++)
		    {
		      source_a += fv->v[a] * bf[var]->grad_phi[j][a];
		    }
		  source_a *= drho_T/rho;       /* derivative of convective terms */
		  source_b =  -drho_T*source*bf[var]->phi[j]/rho;      /* derivative of 1/density */
		  source_c =  + (d2rho_T_c_v*Dcv_Dt + d2rho_T_c_a*Dca_Dt)
		    *bf[var]->phi[j]/rho;        /* derivative of drho_c_v & drho_c_a wrt T */
		  dFVS_dT[j] = source_a + source_b + source_c;	
		}
	    }
	  
	  var = VELOCITY1;
	  if(pd->v[var] )
	    {
	      for( a=0; a<dim; a++)
		{
		  var = VELOCITY1+a;
		  for( j=0; pd->v[var] && j<ei->dof[var]; j++)
		    {
		      dFVS_dv[a][j] = + bf[var]->phi[j]*(    drho_c_v * fv->grad_c[species_v][a]
							   + drho_c_a * fv->grad_c[species_a][a]
							   + drho_c_l * fv->grad_c[species_l][a]
							   + drho_T * fv->grad_T[a])/rho;
		    }
		}
	    }
	  
	  var = MASS_FRACTION;
	  if(pd->v[var])
	    {
	      for( j=0; j<ei->dof[var]; j++)
		{
		  source_a = bf[var]->phi[j] * (1 + 2. * tt)/ dt;
		  for( a=0; a<dim; a++)
		    {
		      source_a += fv->v[a] * bf[var]->grad_phi[j][a];
		    }
		  source_a_v = source_a * drho_c_v/rho;    /* derivative of conc_v convective terms */
		  source_a_a = source_a * drho_c_a/rho;    /* derivative of conc_a convective terms */
		  source_a_l = source_a * drho_c_l/rho;    /* derivative of conc_l convective terms */

		  source_b_v = -drho_c_v*bf[var]->phi[j]/rho*source;  /* derivative of density */
		  source_b_a = -drho_c_a*bf[var]->phi[j]/rho*source;  /* derivative of density */
		  source_b_l = -drho_c_l*bf[var]->phi[j]/rho*source;  /* derivative of density */

		  source_c_v = d2rho_T_c_v *DT_Dt*bf[var]->phi[j]/rho;  /* cross-deriv: derivative of drho_T wrt x */
		  source_c_a = d2rho_T_c_a *DT_Dt*bf[var]->phi[j]/rho;  /* cross-deriv: derivative of drho_T wrt x */

		  dFVS_dC[species_v][j] = source_a_v + source_b_v + source_c_v;
		  dFVS_dC[species_a][j] = source_a_a + source_b_a + source_c_a;
		  dFVS_dC[species_l][j] = source_a_l + source_b_l;
		}
	    }
	  
	}
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME)
    {
      double rho_init, rho_final, aexp, time_delay, realtime;
      rho_init = mp->u_density[0];  /* Initial density */
      rho_final = mp->u_density[1]; /* final density */
      aexp =  mp->u_density[2];      /* Arhennius constant for time exponent*/
      time_delay =  mp->u_density[3]; /* time delay before foaming starts */
      
      if (time > time_delay)
	{
	  realtime = time - time_delay;
	  rho = rho_final + (rho_init - rho_final) * exp(-aexp*realtime);
	  d_rho_dt = -aexp * (rho_init - rho_final) * exp(-aexp*realtime);
	}
      else
	{
	  rho = rho_init;
	  d_rho_dt = 0.0;
	}
      
      source = d_rho_dt/rho;
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
    {
      double T;
      double rho_init  = mp->u_density[0];  /* Initial density (units of gm cm-3 */
      double rho_final = mp->u_density[1];  /* Final density (units of gm cm-3 */
      double cexp      = mp->u_density[2];  /* cexp in units of Kelvin sec */
      double coffset   = mp->u_density[3];  /* offset in units of seconds */
      double time_delay= mp->u_density[4];  /* time delay in units of seconds */
      double drhoDT;
      double d_source_dT = 0.0;
      double  d2_rho_dt_dT = 0.0;
      T = fv->T;
      if (time > time_delay)
	{
	  double realtime = time - time_delay;
	  double delRho = (rho_init - rho_final);
	  double cdenom = cexp - coffset * T;
	  double expT = exp(-realtime * T / cdenom);
	  rho = rho_final + delRho * expT;
	  drhoDT = delRho * expT *(- realtime * cexp / (cdenom * cdenom));
	  d_rho_dt =  - delRho * expT * T / cdenom;
	  d2_rho_dt_dT = - delRho * expT * ( - realtime * cexp * T / (cdenom * cdenom * cdenom) 
					    + cexp / (cdenom * cdenom));
	}
      else
	{
	  rho = rho_init;
	  drhoDT = 0.0;
	  d_rho_dt = 0.0;
	}
      source = d_rho_dt/rho;
      d_source_dT =  d2_rho_dt_dT/rho -  d_rho_dt / rho / rho * drhoDT;
      var = TEMPERATURE;
      if (pd->v[var])
	{
	  for (j = 0; j < ei->dof[var]; j++)
	    {
	      dFVS_dT[j] = d_source_dT * bf[var]->phi[j];
	    }
	}
	
    }

  
  if (ls != NULL && mp->mp2nd != NULL &&
      ( mp->DensityModel != LEVEL_SET) &&
      ( mp->mp2nd->DensityModel == CONSTANT) )
    {
      double factor;
      
      source = ls_modulate_property( source,
				     0.0,
				     ls->Length_Scale,
				     (double) mp->mp2nd->densitymask[0],
				     (double) mp->mp2nd->densitymask[1],
				     dFVS_dF,
				     &factor );
      
      var = TEMPERATURE;
      
      for( j=0; pd->v[var] && j<ei->dof[var]; j++)
	{
	  dFVS_dT[j] *= factor;
	}

      var = VELOCITY1;
      if(pd->v[var] )
	{
	  for( a=0; a<dim; a++)
	    {
	      var = VELOCITY1+a;
	      for( j=0; pd->v[var] && j<ei->dof[var]; j++)
		{
		  dFVS_dv[a][j] *= factor;
		}
	    }
	}
      
      var = MASS_FRACTION;
      if (pd->v[var])
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {

	      /* source = ls_modulate_property( source,
		 0.0,
		 ls->Length_Scale,
		 (double) mp->mp2nd->densitymask[0],
		 (double) mp->mp2nd->densitymask[1],
		 dFVS_dC[w],
		 &factor );*/
	      for( j=0;  j<ei->dof[var]; j++)
		{
 
		  dFVS_dC[w][j] *= factor;
		}
	    }
	}
    }

  return ( source );
}
/*******************************************************************************/
/*******************************************************************************/
double REFVolumeSource (double time,
			 double dt,
			 double tt,
			 double dFVS_dv[DIM][MDE], 
			 double dFVS_dT[MDE],
			 double dFVS_dx[DIM][MDE],
			 double dFVS_dC[MAX_CONC][MDE] )
{
  int a,b,j,w, w1, dim = pd->Num_Dim, var;
  double x_dot[DIM] = {0.,0.,0.};
  double source =0.;
  double phi_j;
  double rho = 0.;
  double sv[MAX_CONC],sv_p=0.;
  double dtmp[MAX_CONC], d2tmp[MAX_CONC][MAX_CONC];

  memset( dFVS_dv, 0 , sizeof(double)*DIM*MDE );
  memset( dFVS_dT, 0 , sizeof(double)*MDE );
  memset( dFVS_dx, 0 , sizeof(double)*DIM*MDE );
  memset( dtmp,  0 , sizeof(double)*MAX_CONC );
  memset( sv,    0 , sizeof(double)*MAX_CONC );
  memset( d2tmp, 0 , sizeof(double)*MAX_CONC*MAX_CONC );
  
   if(mp->SpeciesSourceModel[0] == FOAM)
     {
       foam_species_source(mp->u_species_source[0]);

     } else {
       EH(-1, "Must specify FOAM species source in the material's file");
     } 

  if (pd->TimeIntegration != STEADY ) 
    {
    if (pd->v[MESH_DISPLACEMENT1])
      {
	for (a = 0; a < dim; a++)
	  {
	    x_dot[a] = fv_dot->x[a];
	  }
      }
    }
  rho = density(NULL, time);

  /* must include ref reaction source */

  sv_p = mp->specific_volume[pd->Num_Species_Eqn];

  for (w=0; w<pd->Num_Species_Eqn; w++) 
    {
      sv[w] = mp->specific_volume[w];
      dtmp[w] = -rho*(sv[w]-sv_p) ;
      source += dtmp[w]*mp->species_source[w]; 

  /* must account for spatial variation of density */

  for( a=0; a<dim ; a++)
    {
      source += (fv->v[a] - x_dot[a])*dtmp[w]*fv->grad_c[w][a]; 
    } 

    }

  for(w=0; w<pd->Num_Species_Eqn; w++) {
    for(w1=0;w1<pd->Num_Species_Eqn; w1++)
	{
	  d2tmp[w][w1]=-rho*(sv[w]-sv_p)*dtmp[w1];
	}
  }

  if (  af->Assemble_Jacobian )
    {
      var = VELOCITY1;

      for( j=0; pd->v[var] && j<ei->dof[var]; j++)
	{
	  phi_j = bf[var]->phi[j];

	  for( a=0; a<dim; a++)
	    {
	      dFVS_dv[a][j] = 0.;
	      for(w=0; w<pd->Num_Species_Eqn; w++)
		{
		  dFVS_dv[a][j] += phi_j*dtmp[w]*fv->grad_c[w][a];
		}
	    }
	}

      var = TEMPERATURE;
  
      for( j=0; pd->v[var] && j<ei->dof[var]; j++)
	{
	  dFVS_dT[j] = 0.;
	  phi_j = bf[var]->phi[j];

	  for(w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      dFVS_dT[j] -= (sv[w]-sv_p)*(rho*mp->d_species_source[MAX_VARIABLE_TYPES+w]
				+ mp->species_source[w]*3.*rho*rho) ; 
	    }  
	  dFVS_dT[j] *= phi_j;
	}

      for( b=0; b < dim; b++)
	{
	  var = MESH_DISPLACEMENT1 + b;

	  for( j=0; pd->v[var] && j<ei->dof[var]; j++)
	    {
	      dFVS_dx[b][j] = 0.0;
	      phi_j = bf[var]->phi[j];

	      for(w=0 ; w<pd->Num_Species_Eqn; w++)
		{
		  for ( a=0; a<dim; a++)
		    {
		      dFVS_dx[b][j] += (fv->v[a] - x_dot[a])*dtmp[w]
			            *fv->d_grad_c_dmesh[a][w][b][j];
		    }
		      if( TimeIntegration != STEADY )
			      dFVS_dx[b][j] += -(1.0 + 2*tt)*phi_j*dtmp[w]*fv->grad_c[w][b]/dt; 
		}
	    }
	}

      var = MASS_FRACTION;

      for(w=0 ; w<pd->Num_Species_Eqn; w++)
	{
	  for( j=0; pd->v[var] && j <ei->dof[var]; j++ )
	    {
	      dFVS_dC[w][j] = 0.0;
	      phi_j = bf[var]->phi[j];

	      for(w1=0; w1<pd->Num_Species_Eqn; w1++)
		{
		  for( a=0 ; a<dim; a++ )
		    {
		      dFVS_dC[w][j] += (fv->v[a]-x_dot[a])
			*(d2tmp[w1][w]*fv->grad_c[w1][a] );

		      if(w1 == w) dFVS_dC[w][j] += (fv->v[a]-x_dot[a])
			 *dtmp[w]*bf[var]->grad_phi[j][a];

		    }

		  /* dFVS_dC[w][j] -= (sv[w1]-sv_p)*
		    (rho*mp->Jac_Species_Source[w1*pd->Num_Species+w] +
		    +dtmp[w]*mp->species_source[w1]); */

		  dFVS_dC[w][j] += (dtmp[w1]*mp->Jac_Species_Source[w1*pd->Num_Species+w] +
		    +d2tmp[w1][w]*mp->species_source[w1]); 
		}
	    }
	}
    }
  return ( source );
}

int
assemble_projection_stabilization(Exo_DB *exo)

     /* This routine applies the Dohrmann-Bochev Polynomial Projection Pressure Stabilization
      * to the continuity equation to help aid in iterative solutions the the Navier-Stokes equations
      * This routine also includes an expedient to ignore such terms near a level-set interface. 
      * Not sure how those work, though (PRS 2/27/2007).    this routine is basic, and doesn't necessarily
      * apply for graded meshes.  See assemble_PPPS_generalized routine for that. 
      */
{
  double phi_i, phi_j;
  int var, pvar, ip, i, j, a, b, ipass, num_passes;
  double xi[DIM];
   /* Variables for vicosity and derivative */
  double gamma[DIM][DIM];                  /* shrearrate tensor based on velocity */
  double mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  
  int ip_total = elem_info(NQUAD, ei->ielem_type);
  int eqn = R_PRESSURE;
  int peqn = upd->ep[eqn];
  double wt, det_J, h3, d_vol;
  


  if (ls == NULL)
    {
      num_passes = 1;
    }
  else
    {
      if ( ls->elem_overlap_state )
        num_passes = 2;
      else
        num_passes = 1;
    }
	
  for ( ipass=0; ipass < num_passes; ipass++ )
    {
      if ( ls != NULL )
        {
          if ( num_passes == 2 )
            ls->Elem_Sign = -1 + 2*ipass;
          else
            ls->Elem_Sign = 0;
        } 

      double P_avg = 0.;
      double phi_avg[MDE];
      double vol = 0.;
	  double ls_scale_factor = 1.0;  /* this is to turn off PSPP near the interface */
      
      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_avg[i] = 0.;
        }
        
      for (ip = 0; ip < ip_total; ip++)
        {
          find_stu(ip, ei->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
          wt = Gq_weight (ip, ei->ielem_type); /* find quadrature weights for */

          setup_shop_at_point( ei->ielem, xi, exo );
          fv->wt = wt;
          h3 = fv->h3;
          det_J = bf[eqn]->detJ;
          d_vol = wt*det_J*h3;
          
          P_avg += fv->P * d_vol;

          for ( i=0; i<ei->dof[eqn]; i++)
            {
              phi_i = bf[eqn]->phi[i];
              phi_avg[i] += phi_i * d_vol;
            }
          vol += d_vol;
        }
        
      double inv_vol = 1./vol;
      P_avg *= inv_vol;
      for ( i=0; i<ei->dof[eqn]; i++)
        {
          phi_avg[i] *= inv_vol;
        }
        
      for (ip = 0; ip < ip_total; ip++)
        {
          find_stu(ip, ei->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
          wt = Gq_weight (ip, ei->ielem_type); /* find quadrature weights for */

          setup_shop_at_point( ei->ielem, xi, exo );
          fv->wt = wt;
          h3 = fv->h3;
          det_J = bf[eqn]->detJ;
          d_vol = wt*det_J*h3;
		  
		  if ( ls != NULL && ls->PSPP_filter )
		  {
			load_lsi( ls->Length_Scale );
			
			ls_scale_factor = ls->on_sharp_surf ? 1.0 : 1.0 - lsi->delta/lsi->delta_max;  /*Punting in the case of subelementintegration */
			if ( af->Assemble_Jacobian ) load_lsi_derivs();
		  }
          
          /* load up shearrate tensor based on velocity */
          for ( a=0; a<VIM; a++)
            {
              for ( b=0; b<VIM; b++)
                {
                  gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
                }
            }
          mu = viscosity(gn, gamma, d_mu);
          
          if ( af->Assemble_Residual )
            {
              for ( i=0; i<ei->dof[eqn]; i++)
                {
                  phi_i = bf[eqn]->phi[i];
                  lec->R[peqn][i] += ls_scale_factor*PS_scaling*(fv->P - P_avg) * (phi_i - phi_avg[i]) / mu * d_vol;
                }
            }
          if ( af->Assemble_Jacobian )
            {

              for ( i=0; i<ei->dof[eqn]; i++)
                {
                  phi_i = bf[eqn]->phi[i];
                  
                  var = PRESSURE;
                  pvar = upd->vp[var];
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      lec->J[peqn][pvar][i][j] += ls_scale_factor*PS_scaling*(phi_j - phi_avg[j]) * (phi_i - phi_avg[i]) / mu * d_vol;
                    }
					

			if (pd->v[var = LS]) 
		       {
			 pvar = upd->vp[var];
			 for ( j=0; j<ei->dof[var]; j++)
			 {
			   phi_j = bf[var]->phi[j];

			   if( ls != NULL && ls->PSPP_filter ) {

		           lec->J[peqn][pvar][i][j] -=  ls_scale_factor*PS_scaling*(fv->P - P_avg) * (phi_i - phi_avg[i]) * d_mu->F[j] / (mu* mu) * d_vol ; 

			     lec->J[peqn][pvar][i][j] += -(lsi->d_delta_dF[j]/lsi->delta_max) * PS_scaling*(fv->P - P_avg) * 
						                              (phi_i - phi_avg[i]) / mu * d_vol;
			   }
			  }
		      }
				  
                }
            }
        }
    }
  return ( 1 );
}

int
assemble_PPPS_generalized(Exo_DB *exo)

     /* This routine applies the Dohrmann-Bochev Polynomial Projection Pressure Stabilization
      * to the continuity equation to help aid in iterative solutions the the Navier-Stokes equations
      * Unlike its couterpart assemble_projection_stabilization, this routine adds the expedients necessary to stabilized elements
      * that have significant aspect ratios.   Viz. they most solve an eigenproblem at the element level to do so.  
      */
{
  int var, pvar, ip, i, j, a, b, m, n, k;
  double xi[DIM];
   /* Variables for vicosity and derivative */
  double gamma[DIM][DIM];                  /* shrearrate tensor based on velocity */
  double mu, visc_e = 0.0, d_visc_e_dv[DIM][MDE], d_visc_e_dx[DIM][MDE];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  
  int ip_total = elem_info(NQUAD, ei->ielem_type);
  int eqn = R_PRESSURE;
  int peqn = upd->ep[eqn];
  int INFO, ipass, npass;
  double wt, det_J, h3, d_vol, sum;  
  
  double vol = 0.;
  double xi_dV = 0.;
  double wm_dV = 0.;
  double EeT_We_Ee = 0.;

  /* Declare and initialize some element-level matrices needed for this stabilization (cf. Dorhmann/Bochev notes/paper)
   * IJNMF 46:183-201. 2004.  
   */

  double x_bar[DIM], Vij[DIM][DIM], em[MDE], eim[DIM][MDE], eim_tilde[DIM][MDE];
  double Vij_tilde[DIM][DIM], Ee_tilde[DIM+1][MDE], Ee_hat[DIM+1][MDE];
  double Me[MDE][MDE], Ce[MDE][MDE], W[DIM], eval[DIM], evec[DIM][DIM];
  double De_hat[DIM+1][DIM+1], We[DIM+1][DIM+1], S[DIM+1][DIM+1], WeEe[DIM+1][MDE];
  double WORK[100];
  double *A;
  int LWORK = 100;
  double scaling_over_visc_e, my_multiplier=1.0;
  double R_new[MDE], R_old[MDE];  /*old and new resid vect pieces for the pressure equation */

  asdv(&A, ei->ielem_dim * ei->ielem_dim);

  npass = 1;

  if(pd->v[MESH_DISPLACEMENT1] && af->Assemble_Jacobian)
    {
      npass = 2;
      EH(-1,"need to work out numerical jacobian for PSPP and moving mesh");
    }

  for (ipass = 0; ipass < npass; ipass++)
    {

      memset( Me, 0, sizeof(double)*MDE*MDE); 
      memset( eim, 0, sizeof(double)*DIM*MDE); 
      memset( em, 0, sizeof(double)*MDE); 
      memset( x_bar, 0, sizeof(double)*DIM);
      memset( Vij, 0, sizeof(double)*DIM*DIM); 
      memset( d_visc_e_dv, 0, sizeof(double)*DIM*MDE);
      memset( d_visc_e_dx, 0, sizeof(double)*DIM*MDE);
      memset( R_new, 0, sizeof(double)*MDE);
      memset( R_old, 0, sizeof(double)*MDE);
      //phi_avg[MDE];
      vol = 0.;
      xi_dV = 0.;
      wm_dV = 0.;
      EeT_We_Ee = 0.;
      visc_e=0.;

      for (ip = 0; ip < ip_total; ip++)
	{

	  find_stu(ip, ei->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */

	  wt = Gq_weight (ip, ei->ielem_type); /* find quadrature weights for */

	  setup_shop_at_point( ei->ielem, xi, exo );

	  fv->wt = wt;
	  h3 = fv->h3;
	  det_J = bf[eqn]->detJ;
	  d_vol = wt*det_J*h3;
	  vol += d_vol;

	  /* load up shearrate tensor based on velocity */
	  for ( a=0; a<VIM; a++)
	    {
	      for ( b=0; b<VIM; b++)
		{
		  gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
		}
	    }
	  mu = viscosity(gn, gamma, d_mu);
	  if(gn->ConstitutiveEquation != NEWTONIAN ) EH(-1,"PSPP not available for non_Newtonian yet.   Sorry."); 

	  visc_e += mu*d_vol;
      


	  for (i=0; i < ei->ielem_dim; i++)
	    {
	      xi_dV = fv->x[i] * d_vol; 
	      x_bar[i] += xi_dV;
	      for(j=0; j < ei->ielem_dim; j++)
		{
		  Vij[i][j] += fv->x[j] * xi_dV;	      
		}
	      for(m=0; m < ei->dof[eqn]; m++)
		{
		  eim[i][m] += bf[eqn]->phi[m] * xi_dV;
		}
	    }
    

	  for(m=0; m < ei->dof[eqn]; m++)
	    {
	      wm_dV = bf[eqn]->phi[m] * d_vol;
	      em[m] += wm_dV;
	      for(n=0; n < ei->dof[eqn]; n++)
		{
		  Me[m][n] += bf[eqn]->phi[n] * wm_dV;
		}
	    }

	  /* Pick up some Jacobian pieces for later down in the routine */

	  if ( af->Assemble_Jacobian && ipass != 1)
	    {
	      var = VELOCITY1;
	      if(pd->v[var])
		{
		  for (a = 0; a < ei->ielem_dim; a++)
		    {
		      var = VELOCITY1 + a;
		      
		      for(j = 0; j < ei->dof[var]; j++)
			{
			  d_visc_e_dv[a][j] += d_vol*(*d_mu->v[a]);
			}
		    }
		}
	      
	      var = MESH_DISPLACEMENT1;
	      if(pd->v[var])
		{
		  /* We are doing all of mesh derivatives numerically */
		}
	    }
	}  /* End of integration loop */
      
      visc_e /= vol;  
      if ( af->Assemble_Jacobian)
	{
	  var = VELOCITY1;
	  if(pd->v[var])
	    {
	      for (a = 0; a < ei->ielem_dim; a++)
		{
		  var = VELOCITY1 + a;
		      
		  for(j = 0; j < ei->dof[var]; j++)
		    {
		      d_visc_e_dv[a][j] /= vol;
		    }
		}
	    }    
	} 
      
      for(i=0; i < ei->ielem_dim; i++)
	{
	  x_bar[i] /= vol;
	}
      
      for(i=0; i < ei->ielem_dim; i++)
	{
	  for(m=0; m < ei->dof[eqn]; m++)
	    {
	      eim_tilde[i][m] = eim[i][m] - x_bar[i] * em[m];
	    }
	  for(j = 0; j < ei->ielem_dim; j++)
	    {
	      Vij_tilde[i][j] = Vij[i][j] - vol * x_bar[i] * x_bar[j];
	    }
	}
      for(m=0; m < ei->dof[eqn]; m++)
	{
	  Ee_tilde[0][m] = em[m]; /* First row of Ee_tilde is special */
	  for (i=0; i < ei->ielem_dim; i++)
	    {
	      Ee_tilde[i+1][m] = eim_tilde[i][m];
	    }
	}
      
      /* Next we need to solve the eigen problem */
      
      for(i=0; i < ei->ielem_dim; i++)
	{
	  W[i] = 0;
	  for(j=0 ; j < ei->ielem_dim; j++)
	    {
	      A[j*ei->ielem_dim + i] = Vij_tilde[i][j];
	    }
	}
      
      dsyev_("V", "U", &(ei->ielem_dim), A, &(ei->ielem_dim), W, WORK, &LWORK, &INFO, 1, 1); 
      
      if (INFO > 0) EH(-1, "dsyev falied to converge");
      if (INFO < 0) EH(-1, "an argument of dsyev had an illegal value");
      
      for(j = 0; j < ei->ielem_dim; j++)
	{
	  eval[j] = W[j];
	  for( i=0; i<ei->ielem_dim; i++)
	    {
	      evec[j][i] = A[j*ei->ielem_dim + i];
	    }
	}
      
      /*user De_hat as scratch for calculating Ee_hat */
      
      memset(De_hat, 0, sizeof(double)*(DIM+1)*(DIM+1));
      De_hat[0][0] = 1.;
      for(i=0; i < ei->ielem_dim; i++)
	{
	  for(j=0; j < ei->ielem_dim; j++)
	    {
	      De_hat[i+1][j+1] = evec[i][j];
	    }
	}
      
      for(i=0; i < ei->ielem_dim+1; i++)
	{
	  for(j=0; j < ei->dof[eqn] ; j++)
	    {
	      sum = 0.;
	      for(k=0; k < ei->ielem_dim + 1; k++)
		{
		  sum += De_hat[i][k] * Ee_tilde[k][j];
		}
	      Ee_hat[i][j] = sum; 
	    }
	}
      
      
      /* Compute the real De_hat */
      memset(De_hat, 0, sizeof(double)*(DIM+1)*(DIM+1));
      De_hat[0][0] = vol;
      for(i=0; i<ei->ielem_dim; i++)
	{
	  De_hat[i+1][i+1] = eval[i];
	}
      
      memset(S     , 0, sizeof(double)*(DIM+1)*(DIM+1));
      S[0][0] = 1.;
      for(i=1; i<ei->ielem_dim; i++)
	{
	  S[i+1][i+1] = 1.0 - sqrt(eval[0]/eval[i]);
	}
      
      memset(We    , 0, sizeof(double)*(DIM+1)*(DIM+1));
      We[0][0] = 1./vol;
      for(i=1; i< ei->ielem_dim+1; i++)
	{
	  We[i][i] = S[i][i]*(2.-S[i][i])/De_hat[i][i];
	}
      
      memset(WeEe,0,sizeof(double) * (DIM+1) * MDE);
      for(i=0; i < ei->ielem_dim+1; i++)
	{
	  for(j=0; j < ei->dof[eqn] ; j++)
	    {
	      sum = 0.;
	      for(k=0; k < ei->ielem_dim + 1; k++)
		{
		  sum += We[i][k]* Ee_hat  [k][j];
		}
	      WeEe[i][j] = sum; 
	    }
	}
      
      
      scaling_over_visc_e = my_multiplier * PS_scaling / visc_e;
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  for(j=0; j < ei->dof[eqn]; j++)
	    {
	      EeT_We_Ee = 0.;
	      for(k=0; k < ei->ielem_dim +1; k++)
		{
		  EeT_We_Ee += Ee_hat[k][i] * WeEe[k][j];
		}
	      Ce[i][j] = scaling_over_visc_e * (Me[i][j] - EeT_We_Ee);
	    }
	}
      
      /* Now compute and add stabilizing term to residual equation for continuity*/
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  for(j=0; j < ei->dof[eqn]; j++)
	    {
	      lec->R[peqn][i] += Ce[i][j]* (*esp->P[j]);
	    }
	}

      if(pd->v[MESH_DISPLACEMENT1] && af->Assemble_Jacobian)
	{
	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      for(j=0; j < ei->dof[eqn]; j++)
		{
		  if(ipass == 0)
		    {
		      R_old[i] += Ce[i][j]* (*esp->P[j]);
		    }
		  else
		    {
		      R_new[i] +=  Ce[i][j]* (*esp->P[j]);
		    }
		}
	    }
	}

    } /* End ipass loop */
      
      /*Now compute and jacobian contributions */
  
  if ( af->Assemble_Jacobian )
    {
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  var = PRESSURE;
	  pvar = upd->vp[var];
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      lec->J[peqn][pvar][i][j] += Ce[i][j];
	    }
	}
	  
      var = MESH_DISPLACEMENT1;
	  
      if (pd->v[var]) 
	{
	  for (a = 0; a < ei->ielem_dim; a++)
	    {
	      var += a;
	      /* first the pieces for viscosity */
	      for ( i=0; i<ei->dof[eqn]; i++)
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      /* lec->J[peqn][pvar][i][j] += numerical piece; */
		      EH(-1,"need to finish off numerical jacobian for PPSP mesh displ");
		    }
		}
	    }
	}
      var = VELOCITY1;
      /* need this piece too for d(visc_e)/dv */
      
      for (a = 0; a < ei->ielem_dim; a++)
	{
	  var += a;
	  pvar = upd->vp[var];
	  for ( i=0; i<ei->dof[eqn]; i++)
	    {
		  
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  EeT_We_Ee = 0.;
		  for(k=0; k < ei->ielem_dim +1; k++)
		    {
		      EeT_We_Ee += Ee_hat[k][i] * WeEe[k][j];
		    }
		      
		  lec->J[peqn][pvar][i][j] -= my_multiplier * PS_scaling * 
		    d_visc_e_dv[a][j] * (Me[i][j] - EeT_We_Ee)/visc_e/visc_e;
		}
	    }
	}
		    
    }
	  
	  return ( 1 );
}

/*
 * apply_distributed_sources
 *  
 *    This function applies source terms to the conservation equation as narrow-banded
 *    source terms centered on the zero level set curve.  
 *
 */


int
apply_distributed_sources ( int elem, double width,
                            double x[], Exo_DB *exo,
                            double dt,      /* current time step size                       */
                            double theta,
                            double time,
			    const PG_DATA *pg_data,
                            int oAC,        /* Flag indicating calling function */
                            double *gAC,    /* Augmenting condition arrays */
                            double **bAC,
                            double **cAC )
{

  int ip, ip_total, ielem_type = ei->ielem_type, err;
  struct LS_Embedded_BC *bcref;
  double xi[DIM];
  int ipass, num_elem_passes = 1;
  double ls_F[MDE], ad_wt[10];    /* adaptive integration weights  */
  int i;

  if( pd->e[R_PHASE1] && !pd->e[R_FILL] ) 
	{
	  ls = pfd->ls[0];
	}
    
  if ( xfem != NULL )
    {
      num_elem_passes = 2;
    }
  else
    {
      ls->Elem_Sign = 0;
    }

  if( ls->Integration_Depth > 0 )
    {
      ip_total = Subgrid_Int.ip_total;      
    }
  else
    {
      ip_total = elem_info(NQUAD, ielem_type); 
    }
  if( ls->AdaptIntegration)
    {
	for(i=0 ; i<ei->num_local_nodes ; i++)	{ls_F[i]=*esp->F[i];}
	i = adaptive_weight( ad_wt, ip_total, ei->ielem_dim, ls_F, ls->Length_Scale,
				3, ielem_type);
	WH(i, "problem with adaptive weight routine");
    }

  for( ip=0; ip<ip_total ; ip++ )
    {

      if( ls->Integration_Depth > 0 )
        {
          xi[0] = Subgrid_Int.s[ip][0];
          xi[1] = Subgrid_Int.s[ip][1];
          xi[2] = Subgrid_Int.s[ip][2];
          fv->wt = Subgrid_Int.wt[ip];
        }
      else
        {
    	  find_stu( ip, ielem_type, &(xi[0]), &(xi[1]), &(xi[2]) );
	  if(ls->AdaptIntegration)
		{
		fv->wt = ad_wt[ip_total-1-ip];
		}
	  else
		{
    	  	fv->wt = Gq_weight( ip, ielem_type ) ;
		}
    	}

      for ( ipass = 0; ipass < num_elem_passes; ipass++ )
        {
          if ( num_elem_passes == 2 ) ls->Elem_Sign = -1 + 2*ipass;
               
          err = load_basis_functions( xi, bfd);
          EH( err, "problem from load_basis_functions");

          err = beer_belly();
          EH( err, "beer_belly");

          err = load_fv();
          EH( err, "load_fv");

          err = load_bf_grad();
          EH( err, "load_bf_grad");

          if ( pd->e[R_MESH1] )
            {
              err = load_bf_mesh_derivs();
              EH( err, "load_bf_mesh_derivs");
            }

          err = load_fv_grads();
          EH( err, "load_fv_grads");

          if ( pd->e[R_MESH1] )
            {
              err = load_fv_mesh_derivs(1);
              EH( err, "load_fv_mesh_derivs");
            }

          load_lsi( ls->Length_Scale );
          load_lsi_derivs();

          /* Now where are the material properties in all of this.  Well, here we propagate the
           * unfortunate practice of computing the material properties in the assembly routines
           * themselves.
           */

          bcref = ls->embedded_bc;

          while ( bcref != NULL )
            {
              BOUNDARY_CONDITION_STRUCT *bc;

              bc = BC_Types + bcref->bc_input_id ;

              if ( num_elem_passes == 2 ) 
			  {
				  if (  ( bc->BC_ID != 0 && ls->Elem_Sign != bc->BC_ID) || 
						( bc->BC_ID == 0 && ipass > 0 ) ) 
				  {
					  bcref = bcref->next;
					  continue;
				  }
			  }
                    
              switch (bc->BC_Name )
                {
                case LS_CAPILLARY_BC:
                  assemble_csf_tensor();
                  break;
                case LS_CAP_CURVE_BC:
		  if( pd->e[R_NORMAL1] ) assemble_curvature_with_normals_source () ;
                  else                   assemble_curvature_source ();
                  break;
		case LS_CAP_DIV_N_BC:
		  assemble_div_n_source ();
		  break;
		case LS_CAP_DIV_S_N_BC:
		  assemble_div_s_n_source ();
		  break;
                case LS_Q_BC:
                  assemble_q_source( bc->BC_Data_Float[0] );
                  break;
                case LS_QLASER_BC:
                  assemble_qlaser_source( bc->u_BC, time );
                  break;
                case LS_QVAPOR_BC:
                  assemble_qvapor_source( bc->u_BC);
                  break;
                case LS_QRAD_BC:
                  assemble_qrad_source( bc->BC_Data_Float[0],bc->BC_Data_Float[1],bc->BC_Data_Float[2],bc->BC_Data_Float[3]);
                  break;
                case LS_T_BC:
                  assemble_t_source( bc->BC_Data_Float[0], time );
                  break;
                case LS_LATENT_HEAT_BC:
                  assemble_ls_latent_heat_source ( bc->BC_Data_Float[0], 
			bc->BC_Data_Float[1], dt, theta, time, 
			bcref->bc_input_id, BC_Types );
                  break;
                case LS_YFLUX_BC:
                  assemble_ls_yflux_source ( bc->BC_Data_Int[0], bc->BC_Data_Float[0], bc->BC_Data_Float[1], dt, theta, time,
		                             bcref->bc_input_id, BC_Types );
                  break;
                case LS_CONT_T_BC:
                  assemble_cont_t_source( xi );
                  break;
                case LS_CONT_VEL_BC:
                  assemble_cont_vel_source( xi, exo );
                  break;
                case LS_FLOW_PRESSURE_BC:
                  assemble_p_source( bc->BC_Data_Float[0], bc->BC_Data_Int[0] );
                  break;
                case LS_ACOUSTIC_SOURCE_BC:
                  assemble_ars_source( bc->BC_Data_Float[0],bc->BC_Data_Float[1]);
                  break;
                case LS_RECOIL_PRESSURE_BC:
                  assemble_precoil_source( bc->BC_Data_Float );
                  break;
                case LS_U_BC:
                case LS_V_BC:
                case LS_W_BC:
                  assemble_uvw_source( bc->desc->equation, bc->BC_Data_Float[0] );
                  break;
                case BAAIJENS_FLUID_SOLID_BC:
                  assemble_LM_source( xi, oAC, gAC, bAC, cAC, x, exo );
                  break;
                case LS_EXTV_FLUID_SIC_BC:
                  assemble_interface_extension_velocity_sic( bc->BC_Data_Int[0] );
                  break;
                case LS_EXTV_KINEMATIC_BC:
                  assemble_extv_kinematic ( theta, dt, time, bcref->bc_input_id, BC_Types );
                  break;
                case LS_EIK_KINEMATIC_BC:
                  assemble_eik_kinematic ( theta, dt, time, bcref->bc_input_id, BC_Types );
                  break;
		case PF_CAPILLARY_BC:
		  assemble_pf_capillary ( &(bc->BC_Data_Float[0]));
				  
                  break;
                default:
                  break;
                }
              bcref = bcref->next;
            
			}
          /* equation path dependence terms */
#ifdef COUPLED_FILL
          if ( !ls->Ignore_F_deps )
            {
              if( pd->e[R_FILL] && ipass == 0)
	        {
	          if(  tran->Fill_Equation == FILL_EQN_EIKONAL )
	            {
	              assemble_fill_path_dependence();
	            }
	        }
              if( pd->e[R_MOMENTUM1] )
                {
                  err = assemble_momentum_path_dependence(time, theta, dt, pg_data);
                  EH( err, "assemble_momentum_path_dependence");
                }
              if( pd->e[R_PRESSURE] )
                {
#ifndef DARWIN_HACK
	           err = assemble_continuity_path_dependence(
                     time, theta, dt,pg_data);
                   EH( err, "assemble_continuity_path_dependence");
#endif
                }
	      if( pd->e[R_ENERGY] )
	        {
	          assemble_energy_path_dependence(
	            time, theta, dt, pg_data);
	        }
            }
#endif

        }
    }

  /* leave the place as tidy as it was before you came */
  ls->Elem_Sign = 0;
  
  return ( 2 );

}


int
assemble_pf_capillary (double *pf_surf_tens )
{
  int i,j,a,p,q, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;
  
  int pf;

  struct Basis_Functions *bfm;
  dbl (* grad_phi_i_e_a ) [DIM] = NULL;

  double wt, det_J, h3;

  double csf[DIM][DIM];
  double d_csf_dF[DIM][DIM][MDE];
  double d_csf_dX[DIM][DIM][DIM][MDE];

  
  double source;
  struct Level_Set_Data *ls_save = ls;
  
    eqn = R_MOMENTUM1;	
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;	
  h3 = fv->h3;	
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;



  for( pf=0; pf<pfd->num_phase_funcs; pf++)
  {
	memset( csf,          0, sizeof(double)*DIM*DIM);
	memset( d_csf_dF,     0, sizeof(double)*DIM*DIM*MDE);
	
	ls = pfd->ls[pf];

    continuous_surface_tension(mp->surface_tension*pf_surf_tens[pf], csf, d_csf_dF, d_csf_dX);

    
  /*
   * Wesiduals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  source = 0.;
		  
		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];
	  
		  for ( p=0; p<VIM; p++)
		    {
		      for ( q=0; q<VIM; q++)
			{
			  source += grad_phi_i_e_a[p][q] * csf[q][p]; 
			}
		    }
		      
		  source *= - det_J * wt;
		      
		  source *= h3;
		      
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];


		  grad_phi_i_e_a = bfm->grad_phi_e[i][a];


		  /* J_m_pF
		   */
		  var = PHASE1 + pf ;
		  if( pd->v[var])
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{
			  source = 0.;

			  for ( p=0; p<VIM; p++)
			    {
			      for ( q=0; q<VIM; q++)
				{
				  source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j]; 
				}
			    }

			  source *= -det_J*wt*h3;
			  source *= pd->etm[eqn][LOG2_SOURCE];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
		}
	    }
	}
    }
}

	ls = ls_save;
	
return (1);

}
		      


int
assemble_csf_tensor ( void )
{
	int i,j,a,b,p,q, ii, ledof;
	int eqn, peqn, var, pvar;
	int dim, wim;
	
	struct Basis_Functions *bfm;
	dbl (* grad_phi_i_e_a ) [DIM] = NULL;
	
	double wt, det_J, h3, phi_j, d_area;
	
	double csf[DIM][DIM];
	double d_csf_dF[DIM][DIM][MDE];
	double d_csf_dX[DIM][DIM][DIM][MDE];
	
	double source, source_etm;
	
	eqn = R_MOMENTUM1;	
	if ( ! pd->e[eqn] )
    {
		return(0);
    }
	
	wt = fv->wt;	
	h3 = fv->h3;	
	if ( ls->on_sharp_surf ) /* sharp interface */
    {
		det_J = fv->sdet;
    }
	else              /* diffuse interface */
    {
		det_J = bf[eqn]->detJ;
    }
	
	d_area = wt * det_J*h3;
	
	dim   = pd->Num_Dim;
	wim   = dim;
	if (pd->CoordinateSystem == SWIRLING ||
		pd->CoordinateSystem == PROJECTED_CARTESIAN)
		wim = wim+1;
	
	memset( csf,          0, sizeof(double)*DIM*DIM);
	memset( d_csf_dF,     0, sizeof(double)*DIM*DIM*MDE);
	memset( d_csf_dX,     0, sizeof(double)*DIM*DIM*DIM*MDE);
	
	continuous_surface_tension(mp->surface_tension, csf, d_csf_dF, d_csf_dX);
	
#ifdef COUPLED_FILL
	/* finite difference calculation of path dependencies for
		subelement integration
		*/
	if ( ls->CalcSurfDependencies )
    {
		for( a=0; a<wim; a++ )
		{
			eqn = R_MOMENTUM1 + a;
			peqn = upd->ep[eqn];
			bfm = bf[eqn];
			
			for (i = 0; i < ei->dof[eqn]; i++) 
			{
				
				ledof = ei->lvdof_to_ledof[eqn][i];
				
				if (ei->active_interp_ledof[ledof]) 
				{
					
					ii = ei->lvdof_to_row_lvdof[eqn][i];
					
					grad_phi_i_e_a = bfm->grad_phi_e[i][a];
					
					source = 0.;
					
					for ( p=0; p<VIM; p++)
					{
						for ( q=0; q<VIM; q++)
						{
							source += grad_phi_i_e_a[p][q] * csf[q][p]; 
						}
					}
					
					source *= -det_J * wt * h3;
					source *= pd->etm[eqn][LOG2_SOURCE];
					
					/* J_m_F
						*/
					var = LS ;
					pvar = upd->vp[var];
					
					for (j = 0; j < ei->dof[var]; j++)
					  {
					    phi_j = bf[var]->phi[j];
					    lec->J[peqn][pvar][ii][j] += source * phi_j;
					  }
                }
            }
        }
		return(0);
    }
#endif
	
	/*
	 * Residuals ____________________________________________________________________________
	 */
	
	if ( af->Assemble_Residual ) 
    {
		
		for( a=0; a<wim; a++ )
		{
			eqn = R_MOMENTUM1 + a;
			peqn = upd->ep[eqn];
			bfm = bf[eqn];
			
			for (i = 0; i < ei->dof[eqn]; i++) 
			{
				
				ledof = ei->lvdof_to_ledof[eqn][i];
				
				if (ei->active_interp_ledof[ledof])
				{
					
					ii = ei->lvdof_to_row_lvdof[eqn][i];
					
					grad_phi_i_e_a = bfm->grad_phi_e[i][a];
					
					source = 0.;
					
#ifdef DO_NO_UNROLL
					
					for ( p=0; p<VIM; p++)
					{
						for ( q=0; q<VIM; q++)
						{
							source += grad_phi_i_e_a[p][q] * csf[q][p]; 
						}
					}
#else
					source += grad_phi_i_e_a[0][0] * csf[0][0]; 
					source += grad_phi_i_e_a[1][1] * csf[1][1]; 
					source += grad_phi_i_e_a[0][1] * csf[1][0]; 
					source += grad_phi_i_e_a[1][0] * csf[0][1]; 
					if( VIM == 3)
					{
						source += grad_phi_i_e_a[2][2] * csf[2][2]; 
						source += grad_phi_i_e_a[2][1] * csf[1][2]; 
						source += grad_phi_i_e_a[2][0] * csf[0][2]; 
						source += grad_phi_i_e_a[1][2] * csf[2][1]; 
						source += grad_phi_i_e_a[0][2] * csf[2][0]; 
					}
#endif
					
					
					
					source *= -det_J * wt * h3;
					source *= pd->etm[eqn][LOG2_SOURCE];
					
					lec->R[peqn][ii] += source;
				}
			}
		}
    }
	
	
	/*
	 * Yacobian terms...
	 */
	
	if( af->Assemble_Jacobian )
    {
		for( a=0; a<wim; a++ )
		{
			eqn = R_MOMENTUM1 + a;
			peqn = upd->ep[eqn];
			bfm = bf[eqn];
			source_etm = pd->etm[eqn][LOG2_SOURCE];
			
			for (i = 0; i < ei->dof[eqn]; i++) 
			{
				
				ledof = ei->lvdof_to_ledof[eqn][i];
				
				if (ei->active_interp_ledof[ledof])
				{
					
					ii = ei->lvdof_to_row_lvdof[eqn][i];
					
					grad_phi_i_e_a = bfm->grad_phi_e[i][a];
					
					
					/* J_m_F
						*/
#ifdef COUPLED_FILL
					var = LS ;
					if( pd->v[var] )
					{
						pvar = upd->vp[var];
						
						for( j=0; j<ei->dof[var]; j++)
						{
							phi_j = bf[var]->phi[j];
							
							source = 0.;
							
#ifdef DO_NO_UNROLL 	
							
							for ( p=0; p<VIM; p++)
							{
								for ( q=0; q<VIM; q++)
								{
									source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j]; 
								}
							}
							
#else
							source += grad_phi_i_e_a[0][0] * d_csf_dF[0][0][j];
							source += grad_phi_i_e_a[0][1] * d_csf_dF[1][0][j];
							source += grad_phi_i_e_a[1][0] * d_csf_dF[0][1][j];
							source += grad_phi_i_e_a[1][1] * d_csf_dF[1][1][j];
							
							if( VIM == 3 )
							{
								source += grad_phi_i_e_a[2][0] * d_csf_dF[0][2][j];
								source += grad_phi_i_e_a[2][1] * d_csf_dF[1][2][j];
								source += grad_phi_i_e_a[2][2] * d_csf_dF[2][2][j];
								source += grad_phi_i_e_a[0][2] * d_csf_dF[2][0][j];
								source += grad_phi_i_e_a[1][2] * d_csf_dF[2][1][j];
							}
#endif							
	
							source *= -d_area;
							source *= source_etm;
							
							lec->J[peqn][pvar][ii][j] += source;
						}
					}
#endif
					var = MESH_DISPLACEMENT1;
					if(pd->v[var])
					{
						for( b=0; b<VIM; b++)
						{
							var = MESH_DISPLACEMENT1 + b;
							pvar = upd->vp[var];
							
							for( j=0; j<ei->dof[var]; j++)
							{
								phi_j = bf[var]->phi[j];
								
								source =0.;
#ifdef DO_NO_UNROLL
								for ( p=0; p<VIM; p++)
								{
									for ( q=0; q<VIM; q++)
									{
										source += grad_phi_i_e_a[p][q] * d_csf_dX[q][p][b][j]; 
									}
								} 
#else
								
								source += grad_phi_i_e_a[0][0] * d_csf_dX[0][0][b][j];
								source += grad_phi_i_e_a[0][1] * d_csf_dX[1][0][b][j];
								source += grad_phi_i_e_a[1][0] * d_csf_dX[0][1][b][j];
								source += grad_phi_i_e_a[1][1] * d_csf_dX[1][1][b][j];
								
								if( VIM == 3 )
								{
									source += grad_phi_i_e_a[2][0] * d_csf_dX[0][2][b][j];
									source += grad_phi_i_e_a[2][1] * d_csf_dX[1][2][b][j];
									source += grad_phi_i_e_a[2][2] * d_csf_dX[2][2][b][j];
									source += grad_phi_i_e_a[0][2] * d_csf_dX[2][0][b][j];
									source += grad_phi_i_e_a[1][2] * d_csf_dX[2][1][b][j];
								}
#endif								

								source *= -d_area;
								source *= source_etm;
								
								lec->J[peqn][pvar][ii][j] += source;
							}
						}
					}
								
				}
			}
		}
    }
	
	
	
	return ( 1 );
}


int
assemble_div_n_source ( )
{
  int i, j, a, b, p, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;
  double source, source1;


  eqn = R_MOMENTUM1;
  if (!pd->e[eqn])
    {
      return(-1);
    }

  wt = fv->wt;	
  h3 = fv->h3;	

  if (ls->on_sharp_surf) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ; 
    }
  
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;	

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies)
    {
      for (a = 0; a < wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = -mp->surface_tension*fv->div_n*fv->n[a]*lsi->delta;

		  source *= phi_i*wt*det_J*h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
	        }
            }
        }
      return(0);
    }
#endif

  /*
   * Residuals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = 0.0;

		  source = -mp->surface_tension*fv->div_n*fv->n[a]*lsi->delta;

		  source *= phi_i*wt*det_J*h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];
		  
		
		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }
		  
  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian)
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  /*
		   * J_m_n
		   */
 
		  for( b=0; b<dim; b++)
		    {
		      var = NORMAL1 + b;
		      if( pd->v[var])
			{

			  pvar = upd->vp[var];

			  for( j=0; j<ei->dof[var]; j++)
			    {
		      
			      source = fv->div_n * bf[var]->phi[j]* delta(a,b);
			  
			      for( p=0; p<VIM; p++)
				{
				  source += fv->n[a]*bf[var]->grad_phi_e[j][b][p][p];
				}
			      source *= -mp->surface_tension*lsi->delta;
			      source *= phi_i*wt*det_J*h3;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}
		    }

		  /*
		   * J_m_F
		   */
#ifdef COUPLED_FILL

		  var = LS;
		  if( pd->v[var])
		    {
		      pvar = upd->vp[var];

		      for( j=0; j<ei->dof[var]; j++)
			{
			  source = 0.0;

			  source = -mp->surface_tension*fv->div_n*fv->n[a];
			  source *= lsi->d_delta_dF[j];
			  source *= phi_i*wt*det_J*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
#endif

		  /*
		   * J_m_d
		   */

		  var = MESH_DISPLACEMENT1;

		  if( pd->v[var])
		    {
		      for( b=0; b<dim; b++)
			{
			  var = MESH_DISPLACEMENT1 + b;
			  pvar = upd->vp[var];
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {
			      source = -mp->surface_tension*fv->div_n*fv->n[a]*lsi->delta;
			      
			      source *=  bf[var]->d_det_J_dm[b][j]*h3 + fv->dh3dmesh[b][j]*det_J;
			      source *= wt*phi_i;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];

			      source1 = -mp->surface_tension*fv->d_div_n_dmesh[b][j]*fv->n[a]*lsi->delta;
			      source1 *= wt*phi_i*det_J*h3;
			      source1 *= pd->etm[eqn][(LOG2_SOURCE)];

			      lec->J[peqn][pvar][ii][j] += source + source1;

			    }
			}
		    }
		}
	    }
	}
    }
  return(0);
}



int
assemble_div_s_n_source ( )
{
  int i,j,a,b,p,q, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;
  double source, source1,source2,source3;


  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(-1);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;	

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = -mp->surface_tension*fv->div_s_n*fv->n[a]*lsi->delta;

		  source *= phi_i*wt*det_J*h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
	        }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = 0.0;

		  source = -mp->surface_tension*fv->div_s_n*fv->n[a]*lsi->delta;

		  source *= phi_i*wt*det_J*h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];
		  
		
		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }
		  
  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  /*
		   * J_m_n
		   */

		  for( b=0; b<dim; b++)
		    {
		      var = NORMAL1 + b;
		      if( pd->v[var])
			{

			  pvar = upd->vp[var];

			  for( j=0; j<ei->dof[var]; j++)
			    {
			      source1 = 0.0;
			      for(p=0; p<VIM;p++)
				{
				  source1 -= fv->n[p]*fv->grad_n[p][b] + fv->n[p]*fv->grad_n[b][p];
				}
			      source1 *= fv->n[a] * bf[var]->phi[j];

			      source2 = 0.0;
			      for(p=0; p<VIM;p++)
				{
				  for(q=0;q<VIM;q++)
				    {
				      source2 += (delta(p,q)-fv->n[p]*fv->n[q]) * bf[var]->grad_phi_e[j][b][p][q];
				    }
				}
                              source2 *= fv->n[a];

                              source3 = delta(a,b) * fv->div_s_n * bf[var]->phi[j];

			      source = source1+source2+source3;
			      source *= -mp->surface_tension*lsi->delta;
			      source *= phi_i*wt*det_J*h3;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}
		    }
#ifdef COUPLED_FILL
                  /*
                   * J_m_F
                   */
                  var = FILL;
                  if( pd->v[var])
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          source = 0.0;

                          source = -mp->surface_tension*fv->div_s_n*fv->n[a];
                          source *= lsi->d_delta_dF[j];
                          source *= phi_i*wt*det_J*h3;
                          source *= pd->etm[eqn][(LOG2_SOURCE)];

                         lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                }

	    }
	}
    }
  return(0);
}





int
assemble_curvature_with_normals_source ( )
{
  int i, j, a, b, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;

  /* Bail out with an error if CURVATURE equation not define 
   */

  if (!pd->e[R_CURVATURE]  )
    {
      EH(-1,"Error: Level set curvature equation needs to be activated to use LS_CAP_CURVE\n");
    }

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;	

  /*
   * Wesiduals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {
	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  
		  source =0.0;

		  source = mp->surface_tension*fv->H*fv->n[a]*lsi->delta;

		  source *= wt*phi_i*det_J*h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }


  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  /*
		   * J_m_H
		   */
		  var = CURVATURE;

		  if(pd->v[var])
		    {
		      pvar = upd->vp[var];


		      for( j=0;  j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];
			  
			  source = mp->surface_tension*phi_j*fv->n[a]*lsi->delta;
			  
			  source *= phi_i*wt*det_J*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }

		  /*
		   * J_m_n
		   */

		  var = NORMAL1 + a;

		  if( pd->v[var] )
		    {
		      pvar = upd->vp[var];

		      for( j=0; j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];

			  source = mp->surface_tension*fv->H*phi_j*lsi->delta;
			  source *= phi_i*wt*det_J*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
		  /* 
		   * J_m_F
		   */
		  
		  var = LS;

		  if( pd->v[var])
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{
			  source = mp->surface_tension*fv->H*fv->n[a]*lsi->d_delta_dF[j];
			  source *= phi_i*wt*det_J*h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }

		  /*
		   * J_m_d
		   */

		  var = MESH_DISPLACEMENT1;

		  if( pd->v[var])
		    {
		      for( b=0; b<dim; b++)
			{
			  var = MESH_DISPLACEMENT1 + b;
			  pvar = upd->vp[var];

			  for(j=0; j<ei->dof[var]; j++)
			    {
			      source = mp->surface_tension*fv->H*fv->n[a]*lsi->delta;
			      
			      source *= ( bf[var]->d_det_J_dm[b][j]*h3 + fv->dh3dmesh[b][j]*det_J);
			      source *= wt*phi_i;
			      source *= pd->etm[eqn][(LOG2_SOURCE)];
			    }
			}
		    }
		}
	    }
	}
    }
  return (1);
}



		  
int
assemble_curvature_source ( )
{
  int i, j, a, b, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double f[DIM], dfdX[DIM][DIM][MDE], dfdF[DIM][MDE], dfdH[DIM][MDE];

  /* Bail out with an error if CURVATURE equation not define 
   */

  if (! pd->e[R_CURVATURE] )
    {
      EH(-1,"Error: Level set curvature equation needs to be activated to use LS_CAP_CURVE\n");
    }	

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
    {det_J = fv->sdet;}
  else              /* diffuse interface */
    {det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  memset( dfdH, 0, DIM*MDE*sizeof(double) );
  memset( dfdX, 0, DIM*DIM*MDE*sizeof(double) );
  memset( dfdF, 0, DIM*MDE*sizeof(double) );
  memset(    f, 0, DIM*sizeof(double) );

  curvature_momentum_source( f, dfdX, dfdF, dfdH );

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies)
    {
      for (a = 0; a < wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  
		  source = phi_i*f[a];
		   		      
		  source *= det_J * wt;
		      
		  source *= h3;
		      
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
	        }
            }
        }
      return(0);
    }
#endif

  /*
   * Residuals ________________________________________________________________________________
   */

  if ( af->Assemble_Residual ) 
    {
      
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];
		  
		  source = phi_i*f[a];
		   		      
		  source *= det_J * wt;
		      
		  source *= h3;
		      
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }


  /*
   * Jacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++) 
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof]) 
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];


		  /* J_m_F
		   */
#ifdef COUPLED_FILL
		  var = LS ; 
		  if( pd->v[var])
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{			  
			  source = phi_i*dfdF[a][j];

			  source *= det_J*wt*h3;
			  source *= pd->etm[eqn][LOG2_SOURCE];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }
#endif
		  /*
		   * J_m_H
		   */
		  var = CURVATURE;

		  if ( pd->v[var] )
		    {
		      pvar = upd->vp[var];
		      
		      for( j=0; j<ei->dof[var]; j++)
			{
			  phi_j = bf[var]->phi[j];			  
			  source = phi_i*dfdH[a][j];

			  source *= det_J*wt*h3;
			  source *= pd->etm[eqn][LOG2_SOURCE];

			  lec->J[peqn][pvar][ii][j] += source;
			}
		    }

		  /* J_m_X
		   */

		  for( b=0; b<dim; b++ )
		    {
		      var = R_MESH1 + b;

		      if( pd->v[var] )
			{
			  pvar = upd->vp[var];
			  
			  for( j=0; j<ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];			  
			      source = phi_i*dfdX[a][b][j];
			  
			      source *= det_J*wt*h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}
		    }
		}
	    }
	}
    }
  return ( 1 );
}

int
assemble_q_source(double flux )
{
  int i,j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
        
  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{ det_J = bf[eqn]->detJ; }
	
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];

               source = phi_i * flux * lsi->delta;

               source *= det_J * wt * h3;
               source *= pd->etm[eqn][(LOG2_SOURCE)];

              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
            }
        }
      return(0);
    }
#endif
    
  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];

               source = phi_i * flux * lsi->delta;

               source *= det_J * wt * h3;
               source *= pd->etm[eqn][(LOG2_SOURCE)];

               lec->R[peqn][ii] += source;
            }
        }
    }
  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

#ifdef COUPLED_FILL
              /* diffuse interface version of path dependence integral */
              /*
               * J_e_F
               */
              var = FILL;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      source = phi_i * flux * lsi->d_delta_dF[j];

                      source *= det_J * wt * h3;
                      source *= pd->etm[eqn][LOG2_SOURCE];

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
#endif
            }
        }
    }

  return ( 1 );
}

#if 1
int
assemble_t_source ( double T, double time )
{
  int i,j, ii,ledof, p;
  int eqn, peqn, var, pvar;
  int xfem_active, extended_dof, base_interp, base_dof;
  int incomplete, other_side;
  int apply_NOBC[MDE], apply_SIC[MDE];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double sign = ls->Elem_Sign;
  

  double F_i;

  double q[DIM];
  HEAT_FLUX_DEPENDENCE_STRUCT d_q_struct;
  HEAT_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  double source;

  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if ( pd->i[eqn] != I_Q1_XV &&
       pd->i[eqn] != I_Q2_XV &&
       pd->i[eqn] != I_Q1_GP  &&
       pd->i[eqn] != I_Q2_GP  &&
       pd->i[eqn] != I_Q1_GN  &&
       pd->i[eqn] != I_Q2_GN  &&
       pd->i[eqn] != I_Q1_G  &&
       pd->i[eqn] != I_Q2_G )
    {
      EH(-1, "LS_T requires discontinuous temperature enrichment (Q1_XV, Q2_XV, Q1_G, Q2_G)");
    }

  wt = fv->wt;
  h3 = fv->h3;
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }

  /*
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  /* DRN: This was previously set up to add the penalty term everywhere.
     For the test cases I ran, it seemed to work great.  But when this was
     tried on the momentum eqns (uvw_source), it didn't work.  So the mixed
     nobc/SIC was invented.  So this now applied here too, even though it
     may be overkill.
   */

  heat_flux( q, d_q, time );

  for (i = 0; i < ei->dof[eqn]; i++)
    {
      ledof = ei->lvdof_to_ledof[eqn][i];

      if (ei->active_interp_ledof[ledof])
        {
	  xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                          &xfem_active, &extended_dof, &base_interp, &base_dof );
	  incomplete = dof_incomplete( base_dof, ei->ielem_type, base_interp, ei->ielem_shape);
	  F_i = dof_distance( eqn, i );
	  other_side = sign_change( F_i, (double) ls->Elem_Sign );
	      
	  if ( extended_dof ) {
	    if ( other_side ) {
	      if ( incomplete ) {
		/* extended_dof&&other_side&&incomplete -> NOBC + SIC */
		apply_NOBC[i] = TRUE;
		apply_SIC[i] = TRUE;
	      } else {
		/* extended_dof&&other_side&&!incomplete -> NOBC */
		apply_NOBC[i] = TRUE;
		apply_SIC[i] = FALSE;
	      }
	    } else {
	      /* extended_dof&&!other_side -> do nothing */
	      apply_NOBC[i] = FALSE;
	      apply_SIC[i] = FALSE;
	    }
	  } else {
	    if ( other_side ) {
	      /* !extended_dof&&other_side -> do nothing */
	      apply_NOBC[i] = FALSE;
	      apply_SIC[i] = FALSE;
	    } else {
	      /* !extended_dof&&!other_side -> NOBC */
	      apply_NOBC[i] = TRUE;
	      apply_SIC[i] = FALSE;
	    }
	  }
        }
    }
    
#if 0
{
   double nobc_flux = 0., sic_flux = 0.;
   double r = sqrt(fv->x[0]*fv->x[0] + fv->x[1]*fv->x[1]);
   double angle = atan2(fv->x[1], fv->x[0]);
   int side = 0;
   for ( p=0; p<VIM; p++) nobc_flux += q[p] * lsi->normal[p] * sign;
   sic_flux = (T - fv->T) * BIG_PENALTY;
   if ( r > 0.3 ) side = 1;
   if ( ls->Elem_Sign == -1 ) 
   fprintf(stderr,"FLUX: side, angle, T, NOBC, SIC: %d %g %g %g %g\n",side,angle,fv->T,nobc_flux,sic_flux);
}
#endif

    
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];
	      
	      if ( apply_NOBC[i] )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += q[p] * lsi->normal[p];
                  source *= phi_i * lsi->delta * sign;
                  source *= det_J * wt * h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
	            }
                }
	      if ( apply_SIC[i] )
	        {
                  source = phi_i * lsi->delta * (T - fv->T);
                  source *= det_J * wt * h3 * BIG_PENALTY;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
	            }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];
	      
	      if ( apply_NOBC[i] )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += q[p] * lsi->normal[p];
                  source *= phi_i * lsi->delta * sign;
                  source *= det_J * wt * h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
	      if ( apply_SIC[i] )
	        {
                  source = phi_i * lsi->delta * (T - fv->T);
                  source *= det_J * wt * h3 * BIG_PENALTY;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      eqn = R_ENERGY;
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              if ( apply_NOBC[i] )
                {
                   /*
                    * J_e_T
                    */
                  var = TEMPERATURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_q->T[p][j] * lsi->normal[p];
                          source *= phi_i * lsi->delta * sign;
                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][(LOG2_SOURCE)];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#ifdef COUPLED_FILL
                  /*
                   * J_e_F
                   */
                  var = FILL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_q->F[p][j] * lsi->normal[p] * lsi->delta +
                                                           q[p] * lsi->d_normal_dF[p][j] * lsi->delta +
                                                           q[p] * lsi->normal[p] * lsi->d_delta_dF[j];

                          source *= phi_i * sign;
                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                  /* NOTE: dependencies for X and C (and anything else that modifies q)
                     need to be added here */
                }
              if ( apply_SIC[i] )
	        {
                   /*
                    * J_T_T
                    */
                  var = TEMPERATURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = -phi_i * lsi->delta * phi_j;

                          source *= det_J * wt * h3 * BIG_PENALTY;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                }
	    }
	}
    }
  return ( 1 );
}
#endif

int
assemble_qlaser_source(const double p[], double time )
{
  int i,j, ii,ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux, qlaser;
  double absorp, absorp_base;
  double normal[3] = {0.,0.,0.};

  double source,d_laser_dx[MAX_PDIM];

  eqn = R_ENERGY;
  if (!pd->e[eqn])
    {
      return(0);
    }
  
  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  absorp_base = p[2];
  absorp = absorp_base;

  for (i=0; i<pd->Num_Dim; i++) normal[i] = lsi->normal[i] * ls->Elem_Sign;

  memset(d_laser_dx, 0, sizeof(double)*MAX_PDIM);

  /* Lets CALCULATE LASER FLUX  RARR */
  /* WARNING! sending p as x until we find what method can be used for use_pth and level sets */
  qlaser=calculate_laser_flux(p,time,p,normal,0,1,d_laser_dx);
    
  /* NO Derivatives for d_laser_dx used!!!
   * MUST FIX THIS LATER!!
   */
  
  flux =  absorp*qlaser;
  /* END CALCULATE LASER FLUX */

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
	    }
	}
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
	      lec->R[peqn][ii] += source;
	    }
	}
    }
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      /*
	       * J_T_T
	       */
	      var = TEMPERATURE;
	      
                if ( pd->v[var] )
                  {
                    pvar = upd->vp[var];
		    
                    for( j=0; j<ei->dof[var]; j++)
                      {
                        phi_j = bf[var]->phi[j];
			
                        source = phi_i * lsi->delta * (0.0) * phi_j;

                        source *= det_J * wt * h3;
                        source *= pd->etm[eqn][LOG2_SOURCE];
			
                        lec->J[peqn][pvar][ii][j] += source;
                      }
   		  }
#ifdef COUPLED_FILL
                /*
                 * J_e_F
                 */
                var = FILL;

                if ( pd->v[var] )
                  {
                    pvar = upd->vp[var];

                    for( j=0; j<ei->dof[var]; j++)
                      {
                        source = phi_i * flux * lsi->d_delta_dF[j];
			    
                        source *= det_J * wt * h3;
                        source *= pd->etm[eqn][LOG2_SOURCE];
			    
                        lec->J[peqn][pvar][ii][j] += source;
                      }
                  }
#endif
	    }
	}
    }    
  return ( 1 );
}

int
assemble_qvapor_source ( const double p[])
{
  int i,j, ii,ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j,flux;
  double qvaporloss,d_evap_loss=0.0;
  double source;

  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
    { det_J = fv->sdet; }
  else              /* diffuse interface */
    { det_J = bf[eqn]->detJ; }
   
  /***************************** EXECUTION BEGINS *******************************/
  /* Lets CALCULATE LASER VAPOR HEAT LOSS  RARR */
  qvaporloss=calculate_vapor_cool(p,d_evap_loss,0);
  
  flux = -qvaporloss;
  /* END CALCULATE Q VAPOR */
  
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
	    }
	}
      return(0);
    }
#endif
    
  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
	      lec->R[peqn][ii] += source;
	    }
	}
    }
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      /*
	       * J_T_T
	       */
	      var = TEMPERATURE;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      
		      source = phi_i * lsi->delta * (-d_evap_loss) * phi_j;
		      
		      source *= det_J * wt * h3;
		      source *= pd->etm[eqn][LOG2_SOURCE];
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#ifdef COUPLED_FILL
	      /*
	       * J_e_F
	       */
	      var = FILL;
		  
	      if ( pd->v[var] )
	        {
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
		    {
		      source = phi_i * flux * lsi->d_delta_dF[j];
			  
		      source *= det_J * wt * h3;
		      source *= pd->etm[eqn][LOG2_SOURCE];
			  
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#endif
	    }
	}
    }    
  return ( 1 );
}

int
assemble_qrad_source ( double htc, double Tref, double emiss, double sigma)
{
  int i,j, ii,ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j,flux;
  double d_qrad_dT;
  double source;

  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
    { det_J = fv->sdet; }
  else              /* diffuse interface */
    { det_J = bf[eqn]->detJ; }
   
  flux = htc * (Tref - fv->T) + 
           emiss*sigma * (pow(Tref,4.0) - pow(fv->T,4.0));
  d_qrad_dT = -htc-4.*emiss*sigma * pow(fv->T,3.0) ;

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
	    }
	}
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	      
	      lec->R[peqn][ii] += source;
	    }
	}
    }
  
  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      /*
	       * J_T_T
	       */
	      var = TEMPERATURE;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      
		      source = phi_i * lsi->delta * d_qrad_dT * phi_j;
		      
		      source *= det_J * wt * h3;
		      source *= pd->etm[eqn][LOG2_SOURCE];
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#ifdef COUPLED_FILL
	      /*
	       * J_e_F
	       */
	      var = FILL;
		  
	      if ( pd->v[var] )
	        {
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
		      source = phi_i * flux * lsi->d_delta_dF[j];
			  
		      source *= det_J * wt * h3;
		      source *= pd->etm[eqn][LOG2_SOURCE];
			  
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#endif
	    }
	}
    }    
  return ( 1 );
}

int
assemble_cont_t_source ( double *xi )
{
  int i,j, p, ii,ledof;
  int eqn, peqn, var, pvar;

  double phi_diff[MDE];
  double T_diff;
  double grad_T_diff[DIM];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  int xfem_active, extended_dof, base_interp, base_dof;

  double source;

  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  /* For LS matching conditions, just warn and exit if attempted on
     non-discontinuous interpolations */
  if ( pd->i[eqn] != I_Q1_XV &&
       pd->i[eqn] != I_Q2_XV &&
       pd->i[eqn] != I_Q1_G  &&
       pd->i[eqn] != I_Q2_G )
    {
      WH(-1,"Warning: attempting to apply LS_CONT_T without discontinuous enrichment.\n");
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  /*
   * non-xfem dofs do not need a bc and xfem
   * dofs will receive SIC
   */

  /* DRN: to get this to work for _G based enrichment
     will require making the "natural" dofs act like
     the continuous dofs in _XG type enrichment.
     I think that this could be accomplished via
     a collocated bc, but haven't implemented it yet.
   */
  eqn = R_ENERGY;
  if ( pd->i[eqn] == I_Q1_G || pd->i[eqn] == I_Q2_G )
    EH(-1,"LS_CONT_T_BC not yet implemented for _G type enrichment.");

  /* need difference in T between other side and this side */
  var = TEMPERATURE;
  xfem_var_diff( var, &T_diff, phi_diff, grad_T_diff );

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      /* this is currently handled below using the grad_T_diff */
      return(0);
    }
#endif
  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
               xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                               &xfem_active, &extended_dof, &base_interp, &base_dof );

              if ( extended_dof )
                {
                   ii = ei->lvdof_to_row_lvdof[eqn][i];
                   phi_i = bfm->phi[i];

                   source = phi_i * lsi->delta * T_diff;

                   source *= det_J * wt * h3 * BIG_PENALTY;
                   source *= pd->etm[eqn][(LOG2_SOURCE)];

                   lec->R[peqn][ii] += source;

                }
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      eqn = R_ENERGY;
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {

              xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );

              if ( extended_dof )
                {
                  ii = ei->lvdof_to_row_lvdof[eqn][i];
                  phi_i = bfm->phi[i];

                  /*
                   * J_T_T
                   */
                  var = TEMPERATURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * lsi->delta * phi_diff[j];

                          source *= det_J * wt * h3 * BIG_PENALTY;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
     		    }
                  /* path dependence integral */
                  if ( ls->on_sharp_surf && !ls->AdaptIntegration)
                    {
                      /*
                       * J_e_F
                       */
                      var = FILL;

                      if ( pd->v[var] )
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              phi_j = bf[var]->phi[j];

                              source = 0.;
                              for ( p=0; p < VIM; p++ )
                                {
                                  source -= phi_i * phi_j * lsi->normal[p] * grad_T_diff[p];
                                }

                              for ( p=0; p < VIM; p++ )
                                {
                                  source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * T_diff;
                                }

                              source /= lsi->gfmag;
                              source *= det_J * wt * h3 * BIG_PENALTY;
                              source *= pd->etm[eqn][(LOG2_SOURCE)];


                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }
                    }
                  else
                    {
                      /* diffuse interface version of path dependence integral */
                      /*
                       * J_e_F
                       */
                      var = FILL;
                      if( pd->v[var])
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              source = phi_i * T_diff;
                              source *= lsi->d_delta_dF[j];
                              source *= det_J * wt * h3 * BIG_PENALTY;
                              source *= pd->etm[eqn][(LOG2_SOURCE)];

                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }
                    }
    		}
	    }
	}
    }
  
  return ( 1 );
}

int
assemble_ls_yflux_source ( int wspec,	/* species number of this boundary condition */
		           double mass_tran_coeff, /* Mass transfer coefficient       */
		           double Y_c,	/* bath concentration 	                     */
		           double dt,	/* current value of the time step            */
		           double tt,
                           double time,
			   int bc_input_id,
	                   struct Boundary_Condition *BC_Types )
{
  int i,j, ii,ledof;
  int eqn, peqn, var, pvar, b, w;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j,flux;
  double source;
  
  double Y_w; /* local concentration of current species */
  double sign = 1.;
  double vnorm;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  struct Boundary_Condition *fluxbc;
  
  eqn = R_MASS;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
    { det_J = fv->sdet; }
  else              /* diffuse interface */
    { det_J = bf[eqn]->detJ; }
    
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;
   
  Y_w = fv->c[wspec];
  /*if ( Y_w > 1. ) Y_w = 1.;*/
  mass_flux_surf_mtc(mp->mass_flux, mp->d_mass_flux, fv->T,
		     fv->c, wspec, mass_tran_coeff, Y_c);     
		     
  fluxbc = BC_Types + bc_input_id;
  compute_leak_velocity(&vnorm, d_vnorm, tt, dt, NULL, fluxbc);

  flux = -mp->mass_flux[wspec];

  flux += Y_w * sign * vnorm;
  
  /*
  if ( fv->c[wspec] > 1. )
    fprintf(stderr,"flux=%g, Y_w=%g, mass_flux=%g, extv=%g, vnorm=%g\n",flux,fv->c[wspec],-mp->mass_flux[wspec],fv->ext_v,vnorm);
  */

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = MAX_PROB_VAR+wspec;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
	      
              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
	    }
	}
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = MAX_PROB_VAR+wspec;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
	      
	      lec->R[peqn][ii] += source;
	    }
	}
    }
  
  if ( af->Assemble_Jacobian )
    {
      peqn = MAX_PROB_VAR+wspec;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      /*
	       * J_W_T
	       */
	      var = TEMPERATURE;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      source = -mp->d_mass_flux[wspec][TEMPERATURE] * phi_j;
		      source += Y_w * d_vnorm->T[j] * sign;
		      
		      source *= phi_i * lsi->delta;
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
                
              /*
	       * J_W_V
	       */
	      var = VOLTAGE;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      source = -mp->d_mass_flux[wspec][VOLTAGE] * phi_j;
		      source += Y_w * d_vnorm->V[j] * sign;
		      
		      source *= phi_i * lsi->delta;
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}

              /*
	       * J_W_w
	       */
	      var = MASS_FRACTION;
	      
	      if ( pd->v[var] )
		{
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      for (w = 0; w < pd->Num_Species_Eqn; w++ )
                        {
                          pvar = MAX_PROB_VAR + w;
			  
			  
                          
                          source = -mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] * phi_j;
			  if ( w == wspec ) source += phi_j * vnorm * sign;
			  source += Y_w * d_vnorm->C[w][j] * sign;
		      
		          source *= phi_i * lsi->delta;
		      
		          source *= det_J * wt * h3;
		          /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		          lec->J[peqn][pvar][ii][j] += source;
                        }
		    }
		} 
                
              /*
	       * J_W_v
	       */
              for( b=0; b<wim; b++)
                {
                  var = VELOCITY1 + b;
	          if ( pd->v[var] )
		    {
                      pvar = upd->vp[var];
                      for( j=0; j<ei->dof[var]; j++)
		        {
		          phi_j = bf[var]->phi[j];
                          
                          source = Y_w * d_vnorm->v[b][j] * sign;
		      
		          source *= phi_i * lsi->delta;
		      
		          source *= det_J * wt * h3;
		          /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		          lec->J[peqn][pvar][ii][j] += source;
                        }
		    }
		} 
                
#ifdef COUPLED_FILL
              /*
	       * J_W_F
	       */
	      var = LS;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      source = Y_w * d_vnorm->F[j] * sign;
		      
		      source *= phi_i * lsi->delta;
                      
                      source += flux * phi_i * lsi->d_delta_dF[j];
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#endif
	    }
	}
    }    
  return ( 1 );
}

int
assemble_cont_vel_source ( double *xi, Exo_DB *exo )
{
  int a, i,j, ii,ledof;
  int dim, wim, p;
  int eqn, peqn, var, pvar;
  double source;

  double phi_diff[3][MDE];
  double v_diff[3];
  double grad_v_diff[3][3];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  int xfem_active, extended_dof, base_interp, base_dof;

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  /* For LS matching conditions, just warn and exit if attempted on
     non-discontinuous interpolations */
  if ( pd->i[eqn] != I_Q1_XV &&
       pd->i[eqn] != I_Q2_XV &&
       pd->i[eqn] != I_Q1_G  &&
       pd->i[eqn] != I_Q2_G )
    {
      WH(-1,"Warning: attempting to apply LS_CONT_VEL without discontinuous enrichment.\n");
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  /*
   * non-xfem dofs do not need a bc and xfem
   * dofs will receive SIC
   */

  /* DRN: to get this to work for _G based enrichment
     will require making the "natural" dofs act like
     the continuous dofs in _XG type enrichment.
     I think that this could be accomplished via
     a collocated bc, but haven't implemented it yet.
   */
  eqn = R_MOMENTUM1;
  if ( pd->i[eqn] == I_Q1_G || pd->i[eqn] == I_Q2_G )
    EH(-1,"LS_CONT_VEL_BC not yet implemented for _G type enrichment.");

   for ( a=0; a < wim; a++ )
    {
      var = VELOCITY1 + a;

      /* need difference in vel between other side and this side */
      /* note difference in convection for indices of grad_v_diff from grad_v */
      xfem_var_diff( var, &(v_diff[a]), &(phi_diff[a][0]), &(grad_v_diff[a][0]) );
    }

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      /* this is currently handled below using the grad_v_diff */
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  for ( a=0; a < wim; a++ )
    {
      eqn = R_MOMENTUM1 + a;

      if ( af->Assemble_Residual )
        {
          peqn = upd->ep[eqn];
          bfm = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              ledof = ei->lvdof_to_ledof[eqn][i];

              if (ei->active_interp_ledof[ledof])
                {
                  xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                                  &xfem_active, &extended_dof, &base_interp, &base_dof );

                  if ( extended_dof )
                    {
                      ii = ei->lvdof_to_row_lvdof[eqn][i];
                      phi_i = bfm->phi[i];

                      source = phi_i * lsi->delta * v_diff[a];
                      source *= det_J * wt * h3 * BIG_PENALTY;
                      source *= pd->etm[eqn][(LOG2_SOURCE)];

                      lec->R[peqn][ii] += source;
                    }
                }
            }
        }

      /*
       * Yacobian terms...
       */

      if( af->Assemble_Jacobian )
        {
          peqn = upd->ep[eqn];
          bfm = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              ledof = ei->lvdof_to_ledof[eqn][i];

              if (ei->active_interp_ledof[ledof])
                {
                  xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                                  &xfem_active, &extended_dof, &base_interp, &base_dof );

                  if ( extended_dof )
                    {
                      ii = ei->lvdof_to_row_lvdof[eqn][i];
                      phi_i = bfm->phi[i];

                      /*
                       * J_m_v
                       */
                      var = eqn;

                      if ( pd->v[var] )
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              source = phi_i * lsi->delta * phi_diff[a][j];

                              source *= det_J * wt * h3 * BIG_PENALTY;
                              source *= pd->etm[eqn][LOG2_SOURCE];

                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }

                      /* path dependence integral */
                      if ( ls->on_sharp_surf && !ls->AdaptIntegration)
                        {
                          /*
                           * J_m_F
                           */
                          var = FILL;

                          if ( pd->v[var] )
                            {
                              pvar = upd->vp[var];

                              for( j=0; j<ei->dof[var]; j++)
                                {
                                  phi_j = bf[var]->phi[j];

                                  source = 0.;
                                  for ( p=0; p < VIM; p++ )
                                    {
                                      source -= phi_i * phi_j * lsi->normal[p] * grad_v_diff[a][p];
                                    }

                                  for ( p=0; p < VIM; p++ )
                                    {
                                      source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * v_diff[a];
                                    }

                                  source *= lsi->gfmaginv;
                                  source *= det_J * wt * h3 * BIG_PENALTY;
                                  source *= pd->etm[eqn][(LOG2_SOURCE)];


                                  lec->J[peqn][pvar][ii][j] += source;
                                }
                            }
                        }
                      else
                        {
                          /* diffuse interface version of path dependence integral */
                          /*
                           * J_m_F
                           */
                          var = FILL;
                          if( pd->v[var])
                            {
                              pvar = upd->vp[var];

                              for( j=0; j<ei->dof[var]; j++)
                                {
                                  source = phi_i * v_diff[a];
                                  source *= lsi->d_delta_dF[j];
                                  source *= det_J * wt * h3 * BIG_PENALTY;
                                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                                  lec->J[peqn][pvar][ii][j] += source;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  return ( 1 );
}

int
assemble_extv_kinematic ( dbl tt,		/* parameter to vary time integration from 
		        			 * explicit (tt = 1) to implicit (tt = 0)    */
                          dbl dt,		/* current value of the time step            */
                          dbl time,
	                  int bc_input_id,
	                  struct Boundary_Condition *BC_Types )
{
  int a, i,j, ii,ledof, w ;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double vnorm = 0., coeff = 1.;
  double sign = 1.;
  double tau = 0.;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_coeff_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_coeff = &d_coeff_struct;
  
  struct Boundary_Condition *bc;
  struct Boundary_Condition *fluxbc;

  eqn = R_EXT_VELOCITY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }
    
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  /*  initialize derivatives upfront	*/
  	memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
  	memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
  	memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  	memset(d_vnorm->V, 0, sizeof(dbl)*MDE);
  	memset(d_vnorm->F, 0, sizeof(dbl)*MDE);
  	memset(d_vnorm->X, 0, sizeof(dbl)*MDE*DIM);
  	memset(d_coeff->v, 0, sizeof(dbl)*DIM*MDE);
  	memset(d_coeff->T, 0, sizeof(dbl)*MDE);
  	memset(d_coeff->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  	memset(d_coeff->V, 0, sizeof(dbl)*MDE);
  	memset(d_coeff->F, 0, sizeof(dbl)*MDE);
  	memset(d_coeff->X, 0, sizeof(dbl)*MDE*DIM);

  /* precompute vnorm and its dependencies */
  switch (BC_Types[bc_input_id].BC_Name) {
  case -1: 
#if 0
  case LS_EXTV_KINEMATIC_USER:
#endif
    {
      double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
      double theta = atan2( fv->x[1], fv->x[0] );
      double r0 = 3.;
      vnorm = ((r - r0) + 1.) * (2. + sin(4.*theta));
    }
    break;
   
  case LS_EXTV_KINEMATIC_BC:  /* fluid velocity */
    {
      /* sign is correct on this bc */
      sign = 1.;
      
      vnorm = 0.;
      for ( a=0; a < VIM; a++ ) vnorm += fv->v[a] * lsi->normal[a];

      for ( a=0; a < VIM; a++ )
        {
          var = VELOCITY1 + a;
          if ( pd->v[var] )
            {
              for( j=0; j<ei->dof[var]; j++)
                {
	          phi_j = bf[var]->phi[j];
		  d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
	        }
	    }
        }
    
      var = FILL;
      if ( pd->v[var] )
        {
          for( j=0; j<ei->dof[var]; j++)
            {
	      d_vnorm->F[j] = 0.;
              for ( a=0; a < wim; a++ ) 
                {
                  d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
                }
            }
	}
    }
    break;

  case LS_EXTV_KIN_LEAK_BC:
    {
      /* sign needs to be corrected for this bc */
      sign = 1.*ls->Elem_Sign;
      
      bc = BC_Types + bc_input_id;
      if ( bc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
      else fluxbc = BC_Types + bc->BC_Data_Int[1];
      compute_leak_velocity(&vnorm, d_vnorm, tt, dt, bc, fluxbc);
      
      /* DRN:  This isn't too pretty.  Sorry.
         Basically, we need to know if the fluid velocity is to be included
	 in the extension velocity or if it will be handled separately in 
	 the level set advection equation.
	 If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
	 the extension velocity to INCLUDE the fluid velocity.
	 Otherwise, we will EXCLUDE the fluid velocity and assume that
	 the level set advection equation will take this into account
      */
      if ( pd->v[VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V )
        {
          for ( a=0; a < VIM; a++ ) vnorm += sign * fv->v[a] * lsi->normal[a];
  
          for ( a=0; a < VIM; a++ )
            {
              var = VELOCITY1 + a;
              if ( pd->v[var] )
                {
                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->v[a][j] = sign * lsi->normal[a] * phi_j;
                    }
                }
            }
    
          var = FILL;
          if ( pd->v[var] )
            {
              for( j=0; j<ei->dof[var]; j++)
                {
                  d_vnorm->F[j] = 0.;
                  for ( a=0; a < wim; a++ ) 
                    {
                      d_vnorm->F[j] += sign * fv->v[a] * lsi->d_normal_dF[a][j];
                    }
                }
            }
        }
    }
    break;

  case LS_EXTV_LATENT_BC:
    {
      /* sign needs to be corrected for this bc */
      
      bc = BC_Types + bc_input_id;
      if ( bc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
      else fluxbc = BC_Types + bc->BC_Data_Int[1];

	tau = bc->BC_Data_Float[0];
	vnorm = bc->BC_Data_Float[1]*(fv->T - fluxbc->BC_Data_Float[0]);

	/* determine sign based on normal temperature gradient	*/
	coeff = 0.0;
        for ( a=0; a < dim; a++ ) coeff += fv->grad_T[a] * lsi->normal[a];
	coeff *= tran->delta_t_avg;

	sign = -1.;
	var = TEMPERATURE;
        if ( pd->v[var] )
        {
        for( j=0; j<ei->dof[var]; j++)
           {
            d_vnorm->T[j] += bc->BC_Data_Float[1]*bf[var]->phi[j];
            for ( a=0; a < dim; a++ ) 
		d_coeff->T[j] += bf[var]->grad_phi[j][a] * lsi->normal[a];
	    d_coeff->T[j] *= tran->delta_t_avg;
           }
        }

        var = ls->var;
        if ( pd->v[var] )
        {
        for( j=0; j<ei->dof[var]; j++)
           {
            for ( a=0; a < dim; a++ ) 
               {
                d_coeff->F[j] += fv->grad_T[a] * lsi->d_normal_dF[a][j];
               }
	    d_coeff->F[j] *= tran->delta_t_avg;
           }
        }

      
      /* DRN:  This isn't too pretty.  Sorry.
         Basically, we need to know if the fluid velocity is to be included
	 in the extension velocity or if it will be handled separately in 
	 the level set advection equation.
	 If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
	 the extension velocity to INCLUDE the fluid velocity.
	 Otherwise, we will EXCLUDE the fluid velocity and assume that
	 the level set advection equation will take this into account
      */
      if ( pd->v[VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V )
        {
          for ( a=0; a < VIM; a++ ) vnorm += fv->v[a] * lsi->normal[a];
  
          for ( a=0; a < VIM; a++ )
            {
              var = VELOCITY1 + a;
              if ( pd->v[var] )
                {
                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
                    }
                }
            }
    
          var = ls->var;
          if ( pd->v[var] )
            {
              for( j=0; j<ei->dof[var]; j++)
                {
                  d_vnorm->F[j] = 0.;
                  for ( a=0; a < wim; a++ ) 
                    {
                      d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
                    }
                }
            }
        }
    }
    break;
    
  default:
    sprintf(Err_Msg, "BC %s not found", BC_Types[bc_input_id].desc->name1);
    EH(-1, Err_Msg);
    break;
  }
#if 0
  DPRINTF(stderr,"kinematic extv at (%g,%g) vnorm=%g, fdot=%g\n",fv->x[0],fv->x[1],vnorm,fv_dot->F);
#endif

#define USE_GFMAG 1
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];
	  
          if (ei->active_interp_ledof[ledof])
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];

#if USE_GFMAG
	       source = 2. * lsi->delta * 
			(-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign) * lsi->gfmag;
#else
	       source = 2. * lsi->delta * 
			(-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign);
#endif

               source *= phi_i * det_J * wt * h3;
               /* no such term multiplier yet
                 source *= pd->etm[eqn][(LOG2_SOURCE)];
               */

	       /* J_m_F
		   */
              var = ls->var ;
              pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];
		  lec->J[peqn][pvar][ii][j] += source * phi_j;
		}
            }
        }
      return(0);  
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];
	  
          if (ei->active_interp_ledof[ledof])
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];
	       
#if USE_GFMAG
	       source = 2. * lsi->delta * 
			(-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign) * lsi->gfmag;
#else
	       source = 2. * lsi->delta * 
			(-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign);
#endif

               source *= phi_i * det_J * wt * h3;
               /* no such term multiplier yet
                 source *= pd->etm[eqn][(LOG2_SOURCE)];
               */
               lec->R[peqn][ii] += source;
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              /*
               * J_ext_v_ext_v
               */
/*              var = eqn;  */
              var = EXT_VELOCITY;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

#if USE_GFMAG
                      source = 2. * lsi->delta * lsi->gfmag * 
				(-coeff*(tau*(1.+2.*tt)*phi_j/dt + phi_j));
#else
                      source = 2. * lsi->delta * 
				(-coeff*(tau*(1.+2.*tt)*phi_j/dt + phi_j));
#endif

                      source *= phi_i * det_J * wt * h3;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }

              /*
               * J_ext_v_v
               */
              for ( a=0; a < VIM; a++ )
                {
                  var = VELOCITY1 + a;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

#if USE_GFMAG
                          source = 2. * lsi->delta * lsi->gfmag * 
				(-d_coeff->v[a][j]*(tau*fv_dot->ext_v + fv->ext_v) 
				+ sign*d_vnorm->v[a][j] );
#else
                          source = 2. * lsi->delta * 
				(-d_coeff->v[a][j]*(tau*fv_dot->ext_v + fv->ext_v) 
				+ sign*d_vnorm->v[a][j] );
#endif

                          source *= phi_i * det_J * wt * h3;
                          /* no such term multiplier yet
                            source *= pd->etm[eqn][(LOG2_SOURCE)];
                          */

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                }

              /*
               * J_ext_v_F
               */
              var = ls->var;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
		      
#if USE_GFMAG
		      source = 2. * lsi->delta * lsi->gfmag * 
			(-d_coeff->F[j]*(tau*fv_dot->ext_v +fv->ext_v) + d_vnorm->F[j]*sign) +
		               2. * (-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign) * 
		        (lsi->delta * lsi->d_gfmag_dF[j] + lsi->d_delta_dF[j] * lsi->gfmag);
#else
		      source = 2. * lsi->delta * 
			(-d_coeff->F[j]*(tau*fv_dot->ext_v +fv->ext_v) + d_vnorm->F[j]*sign) +
		               2. * (-coeff*(tau*fv_dot->ext_v + fv->ext_v) + vnorm*sign) * 
		        lsi->d_delta_dF[j];
#endif


                      source *= phi_i * det_J * wt * h3;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
                
              /*
               * J_ext_v_T
               */
              var = TEMPERATURE;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
		      
#if USE_GFMAG
	       	      source = 2. * lsi->delta * lsi->gfmag *
			(-d_coeff->T[j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->T[j]*sign);
#else
	       	      source = 2. * lsi->delta * 
			(-d_coeff->T[j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->T[j]*sign);
#endif
                      source *= phi_i * det_J * wt * h3;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }

              /*
               * J_ext_v_V
               */
              var = VOLTAGE;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
		      
#if USE_GFMAG
	       	      source = 2. * lsi->delta * lsi->gfmag *
			(-d_coeff->V[j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->V[j]*sign);
#else
	       	      source = 2. * lsi->delta * 
			(-d_coeff->V[j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->V[j]*sign);
#endif
                      source *= phi_i * det_J * wt * h3;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
                               
              /*
               * J_ext_v_C
               */
              var = MASS_FRACTION;
              if ( pd->v[var] )
                {
                  for( j=0; j<ei->dof[var]; j++)
                    {
                      for ( w=0; w<pd->Num_Species_Eqn; w++)
			{
                          phi_j = bf[var]->phi[j]; 
#if USE_GFMAG
	       	          source = 2. * lsi->delta * lsi->gfmag *
		(-d_coeff->C[w][j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->C[w][j]*sign);
#else
	       	          source = 2. * lsi->delta * 
		(-d_coeff->C[w][j]*(tau*fv_dot->ext_v + fv->ext_v) + d_vnorm->C[w][j]*sign);
#endif
                          source *= phi_i * det_J * wt * h3;
                          /* no such term multiplier yet
                             source *= pd->etm[eqn][(LOG2_SOURCE)];
                           */

                          lec->J[peqn][MAX_PROB_VAR + w][ii][j] += source;
                        }
                    }
                }
            }
        }
    }

  return ( 1 );
}

int
assemble_interface_extension_velocity_sic ( int ext_vel_sign )
{
  int a, i,j, ii,ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double F_i;

  eqn = R_EXT_VELOCITY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

	  F_i = dof_distance( eqn, i );
	  
          if (ei->active_interp_ledof[ledof] &&
              !sign_change( F_i, (double) ext_vel_sign ) )
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];

               source = phi_i * lsi->delta * fv->ext_v;

               for ( a=0; a < VIM; a++ )
                 {
                   source -= phi_i * lsi->delta * fv->v[a] * lsi->normal[a];
                 }

               source *= det_J * wt * h3 * BIG_PENALTY;
               /* no such term multiplier yet
                 source *= pd->etm[eqn][(LOG2_SOURCE)];
               */

               lec->R[peqn][ii] += source;
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          F_i = dof_distance( eqn, i );
          if (ei->active_interp_ledof[ledof] &&
              !sign_change( F_i, (double) ext_vel_sign ) )
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              /*
               * J_ext_v_ext_v
               */
              var = eqn;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

                      source = phi_i * lsi->delta * phi_j;

                      source *= det_J * wt * h3 * BIG_PENALTY;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }

              /*
               * J_ext_v_v
               */
              for ( a=0; a < VIM; a++ )
                {
                  var = VELOCITY1 + a;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = -phi_i * lsi->delta * phi_j * lsi->normal[a];

                          source *= det_J * wt * h3 * BIG_PENALTY;
                          /* no such term multiplier yet
                            source *= pd->etm[eqn][(LOG2_SOURCE)];
                          */

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                }

              /*
               * J_ext_v_F
               */
              var = FILL;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

                      source = phi_i * lsi->d_delta_dF[j] * fv->ext_v;

                      for ( a=0; a < VIM; a++ )
                        {
                          source -= phi_i * lsi->d_delta_dF[j] * fv->v[a] * lsi->normal[a];
                        }

                      for ( a=0; a < VIM; a++ )
                        {
                          source -= phi_i * lsi->delta * fv->v[a] * lsi->d_normal_dF[a][j];
                        }

                      source *= det_J * wt * h3 * BIG_PENALTY;
                      /* no such term multiplier yet
                        source *= pd->etm[eqn][(LOG2_SOURCE)];
                      */

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
            }
        }
    }

  return ( 1 );
}


int
assemble_eik_kinematic ( dbl tt,		/* parameter to vary time integration from 
		        			 * explicit (tt = 1) to implicit (tt = 0)    */
                         dbl dt,		/* current value of the time step            */
                         dbl time,
	                 int bc_input_id,
	                 struct Boundary_Condition *BC_Types )
{
  int a,  i, j, ii, ledof, w;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double vnorm;
  double sign = 1.;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  double mass = 0.;
  double advection = 0.;
  double source = 0.;
  double penalty = 1.;
  
  struct Boundary_Condition *bc;
  struct Boundary_Condition *fluxbc;


  eqn = R_FILL;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;
  h3 = fv->h3;
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }

  /* precompute vnorm and its dependencies */
  switch (BC_Types[bc_input_id].BC_Name) {
  case -1: 
#if 0
  case LS_EIK_KINEMATIC_USER:
#endif
    {
      double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
      double theta = atan2( fv->x[1], fv->x[0] );
      double r0 = 3.;
      vnorm = ((r - r0) + 1.) * (2. + sin(4.*theta));
      
      memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
      memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
      memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
      memset(d_vnorm->V, 0, sizeof(dbl)*MDE);
      memset(d_vnorm->F, 0, sizeof(dbl)*MDE);
    }
    break;
    
   case LS_EIK_KINEMATIC_BC:  /* fluid velocity */
    {
      /* sign is correct on this bc */
      sign = 1.;
      
      vnorm = 0.;
#if 0
      for ( a=0; a < VIM; a++ ) vnorm += fv->v[a] * lsi->normal[a];
#else
      for ( a=0; a < VIM; a++ ) vnorm += fv_old->v[a] * lsi->normal[a];
#endif

      memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
      memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
      memset(d_vnorm->V, 0, sizeof(dbl)*MDE);
  
#if 0
      for ( a=0; a < wim; a++ )
        {
          var = VELOCITY1 + a;
          if ( pd->v[var] )
            {
              for( j=0; j<ei->dof[var]; j++)
                {
	          phi_j = bf[var]->phi[j];
		  d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
	        }
	    }
        }
      var = FILL;
      if ( pd->v[var] )
        {
          for( j=0; j<ei->dof[var]; j++)
            {
	      d_vnorm->F[j] = 0.;
              for ( a=0; a < wim; a++ ) 
                {
                  d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
                }
            }
	}
#else
      memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
      
      var = FILL;
      if ( pd->v[var] )
        {
          for( j=0; j<ei->dof[var]; j++)
            {
	      d_vnorm->F[j] = 0.;
              for ( a=0; a < wim; a++ ) 
                {
                  d_vnorm->F[j] += fv_old->v[a] * lsi->d_normal_dF[a][j];
                }
            }
	}
#endif
    }
    break;

  case LS_EIK_KIN_LEAK_BC:
    {
      /* sign needs to be corrected for this bc */
      sign = 1.*ls->Elem_Sign;
      
      bc = BC_Types + bc_input_id;
      if ( bc->BC_Data_Int[1] == -1 ) fluxbc = NULL;
      else fluxbc = BC_Types + bc->BC_Data_Int[1];
      compute_leak_velocity(&vnorm, d_vnorm, tt, dt, bc, fluxbc);
      
      /* DRN:  This isn't too pretty.  Sorry.
         Basically, we need to know if the fluid velocity is to be included
	 in the extension velocity or if it will be handled separately in 
	 the level set advection equation.
	 If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
	 the extension velocity to INCLUDE the fluid velocity.
	 Otherwise, we will EXCLUDE the fluid velocity and assume that
	 the level set advection equation will take this into account
      */
      if ( pd->v[VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V )
        {
          for ( a=0; a < VIM; a++ ) vnorm += sign * fv->v[a] * lsi->normal[a];
  
          for ( a=0; a < VIM; a++ )
            {
              var = VELOCITY1 + a;
              if ( pd->v[var] )
                {
                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
                      d_vnorm->v[a][j] = sign * lsi->normal[a] * phi_j;
                    }
                }
            }
    
          var = FILL;
          if ( pd->v[var] )
            {
              for( j=0; j<ei->dof[var]; j++)
                {
                  d_vnorm->F[j] = 0.;
                  for ( a=0; a < wim; a++ ) 
                    {
                      d_vnorm->F[j] += sign * fv->v[a] * lsi->d_normal_dF[a][j];
                    }
                }
            }
        }
    }
    break;
    
  default:
    sprintf(Err_Msg, "BC %s not found", BC_Types[bc_input_id].desc->name1);
    EH(-1, Err_Msg);
    break;
  }
#if 0
  DPRINTF(stderr,"kinematic extv at (%g,%g) vnorm=%g, fdot=%g\n",fv->x[0],fv->x[1],vnorm,fv_dot->F);
#endif
    
  if ( tt != 0. ) EH(-1,"LS_EIK_KINEMATIC currently requires backward Euler");
  
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if ( ei->active_interp_ledof[ledof] )
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];
	       
#if 0
	       mass = lsi->delta * 2. * (-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 1
	       mass = lsi->delta * 2. * (fv->F-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 0
	       advection = lsi->delta * 2. * (dt * sign*vnorm);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
	       advection = lsi->delta * 2. * (dt * sign*vnorm * lsi->gfmag);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
	       
               /* J_m_F
                */
               var = LS ;
               pvar = upd->vp[var];
		      
               for( j=0; j<ei->dof[var]; j++)
                 {
                   phi_j = bf[var]->phi[j];
                   lec->J[peqn][pvar][ii][j] += (mass + advection + source) * phi_j;
                 }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if ( ei->active_interp_ledof[ledof] )
            {
               ii = ei->lvdof_to_row_lvdof[eqn][i];

               phi_i = bfm->phi[i];

#if 0
	       mass = lsi->delta * 2. * (-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 1
	       mass = lsi->delta * 2. * (fv->F-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 0
	       advection = lsi->delta * 2. * (dt * sign*vnorm);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
	       advection = lsi->delta * 2. * (dt * sign*vnorm * lsi->gfmag);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif

               lec->R[peqn][ii] += mass + advection + source;

            }
        }

#if 0
{
  fprintf(stderr,"checkkinematic %g %g %g %g %g %g\n",fv->x[0],fv->x[1],fv_dot->F + vnorm * lsi->gfmag,
          -(fv->F - fv_old->F)/dt, vnorm, lsi->gfmag);
}
#endif
#if 0
{
  double theta = atan2( fv->x[1], fv->x[0] );
  double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
  double rsoln = 2.+exp(time*(2. + sin(4.*theta)));
  rsoln = 2.-exp(time*(2. - sin(4.*theta))) + exp(time) + exp(3.*time);
  fprintf(stderr,"checkkinematic %g %g %g %g %g %g\n",theta,r,fv_dot->F + vnorm * lsi->gfmag, 
          -(fv->F - fv_old->F)/dt, vnorm, lsi->gfmag);
}
#endif
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if ( ei->active_interp_ledof[ledof] )
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              /*
               * J_ls_v
               */
              for ( a=0; a < VIM; a++ )
                {
                  var = VELOCITY1 + a;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

#if 0
                          advection = lsi->delta * 2. * dt * sign*d_vnorm->v[a][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
                          advection = lsi->delta * 2. * dt * lsi->gfmag * sign*d_vnorm->v[a][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif

                          lec->J[peqn][pvar][ii][j] += advection;
                        }
                    }
                }

              /*
               * J_ls_T
               */
              var = TEMPERATURE;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
		      
#if 0
                      advection = lsi->delta * 2. * dt * sign*d_vnorm->T[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
                      advection = lsi->delta * 2. * dt * lsi->gfmag * sign*d_vnorm->T[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif

                      lec->J[peqn][pvar][ii][j] += advection;
                    }
                }
		
              /*
               * J_ls_V
               */
              var = VOLTAGE;
              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];
		      
#if 0
                      advection = lsi->delta * 2. * dt * sign*d_vnorm->V[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
                      advection = lsi->delta * 2. * dt * lsi->gfmag * sign*d_vnorm->V[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif

                      lec->J[peqn][pvar][ii][j] += advection;
                    }
                }
		
              /*
               * J_ls_C
               */
              var = MASS_FRACTION;
              if ( pd->v[var] )
                {
                  for( j=0; j<ei->dof[var]; j++)
                    {
                      for ( w=0; w<pd->Num_Species_Eqn; w++)
                        {
                          phi_j = bf[var]->phi[j];

#if 0
                          advection = lsi->delta * 2. * dt * sign*d_vnorm->C[w][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
                          advection = lsi->delta * 2. * dt * lsi->gfmag * sign*d_vnorm->C[w][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
		      
                          lec->J[peqn][MAX_PROB_VAR + w][ii][j] += advection;
			}
                    }
                }
		
              /*
               * J_ls_F
               */
              var = FILL;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      phi_j = bf[var]->phi[j];

#if 0
		      mass = (lsi->d_delta_dF[j] * 2.) * (-fv_old->F);
		      mass *= phi_i * det_J * wt * h3 * penalty;
		      mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 1
		      mass = lsi->d_delta_dF[j] * 2. * (fv->F-fv_old->F) + lsi->delta * 2. * phi_j;
		      mass *= phi_i * det_J * wt * h3 * penalty;
		      mass *= pd->etm[eqn][(LOG2_MASS)];
#endif
#if 0
		      advection = 2. * dt *
		                  (lsi->d_delta_dF[j] * vnorm + 
		                   lsi->delta * sign*d_vnorm->F[j]);
		      advection *= phi_i * det_J * wt * h3 * penalty;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
#if 1
                      advection = 2. * dt *
		                  (vnorm * lsi->gfmag * lsi->d_delta_dF[j] +
				   lsi->delta * vnorm * lsi->d_gfmag_dF[j] +
				   lsi->gfmag * lsi->delta * sign*d_vnorm->F[j] );
		      advection *= phi_i * det_J * wt * h3 * penalty;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
#endif
		      
                      lec->J[peqn][pvar][ii][j] += mass + advection + source;
                    }
                }
            }
        }
    }

  return ( 1 );
}

int
assemble_p_source ( double pressure, const int bcflag )
{
  int i,j, a, b, p, q, ii,ledof, w;
  int eqn, peqn, var, pvar;
  int dim, wim;
  int err, elem_sign;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;

  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  double dTT_drs[MAX_PDIM][MAX_PDIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  double dTT_dp_liq[DIM][DIM][MDE];
  double dTT_dp_gas[DIM][DIM][MDE];
  double dTT_dporosity[DIM][DIM][MDE];
  double dTT_dsink_mass[DIM][DIM][MDE];
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[DIM][DIM][MDE];
  double dTT_dcur_strain[DIM][DIM][MDE];
  double elast_modulus;

  double Pi[DIM][DIM];		/** liquid stresses  **/
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;



  double force[MAX_PDIM], d_force[MAX_PDIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE];


  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;
         elem_sign=ls->Elem_Sign;
        }
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; 
         elem_sign=1.;
        }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  memset( force, 0, MAX_PDIM*sizeof(double));
  memset( d_force, 0, MAX_PDIM*(MAX_VARIABLE_TYPES+MAX_CONC)*MDE*sizeof(double));

  if ( bcflag == -1)		/**   add constant, isotropic stress  **/
	{
         for( a=0; a<wim; a++ )
		{ force[a] = -pressure*lsi->normal[a];}
  	 if ( af->Assemble_Jacobian )
    		{
                  var = ls->var;
                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];
                      for( j=0; j<ei->dof[var]; j++)
                        {
      			    for( a=0; a<wim; a++ )
				{
			  d_force[a][var][j] -= pressure*lsi->d_normal_dF[a][j];
				}
                        }
                    }
    		}

	}
  else if ( bcflag == 0 ) 	/**  subtract off viscous stress  **/
	{
	 if (!af->Assemble_Jacobian) d_Pi = NULL; 
  		/* compute stress tensor and its derivatives */
	  fluid_stress( Pi, d_Pi );

         for( a=0; a<wim; a++ )
	     {
              for( i=0; i<wim; i++ )
		{ force[a] = -pressure * Pi[a][i]*lsi->normal[i];}
	     }

  	  if (af->Assemble_Jacobian)
    		{
      var = TEMPERATURE;
      if (pd->v[var] )
        {
         for (p=0; p<pd->Num_Dim; p++)
      	    { 
              for (q=0; q<pd->Num_Dim; q++)
              	 { 
                    for (j=0; j<ei->dof[var]; j++)
                      { 
                      d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->T[p][q][j];
                      }
                  }
             }
        }
      var = BOND_EVOLUTION;
      if (pd->v[var] )
        {

          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                { 
                  for (j=0; j<ei->dof[var]; j++)
                    {
                      d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->nn[p][q][j];
                    }
                }
            }
        }

      var = ls->var;
      if (pd->v[var] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for (j=0; j<ei->dof[var]; j++)
                    {
                      d_force[p][var][j] -= pressure*(lsi->normal[q]*d_Pi->F[p][q][j]
				+  lsi->d_normal_dF[p][j]*Pi[p][q]);
                    }
                }
            }
        }
      var = MASS_FRACTION;
      if (pd->v[var] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for ( w=0; w<pd->Num_Species_Eqn; w++)
                    {
                      for (j=0; j<ei->dof[var]; j++)
                        {
                          d_force[p][MAX_VARIABLE_TYPES + w][j] -= pressure*lsi->normal[q]*d_Pi->C[p][q][w][j];
                        }
                    }
                }
            }
        }
          var = PRESSURE;
          if (pd->v[var] )
            {
              for (p=0; p<pd->Num_Dim; p++)
                {
                  for (q=0; q<pd->Num_Dim; q++)
                    {
                      for (j=0; j<ei->dof[var]; j++)
                        {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->P[p][q][j];
                        }
                    }
                }
            }
      if (pd->v[VELOCITY1] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for ( b=0; b<VIM; b++)
                    {
                      var = VELOCITY1+b;
                      for (j=0; j<ei->dof[var]; j++)
                        {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->v[p][q][b][j];
                        }
                    }
                }
            }
        }
     if (pd->v[VORT_DIR1] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for ( b=0; b<VIM; b++)
                    {
                      var = VORT_DIR1+b;
                      for (j=0; j<ei->dof[var]; j++)
                        {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->vd[p][q][b][j];
                        }
                    }
                }
            }
        }
      if (pd->v[MESH_DISPLACEMENT1] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for ( b=0; b<VIM; b++)
                    {
                      var = MESH_DISPLACEMENT1+b;
                      for (j=0; j<ei->dof[var]; j++)
                        {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*d_Pi->X[p][q][b][j];
                        }
                    }
                }
            }
        }
    		}	/*  if Jacobian  */
	}
  else if( bcflag == 1 && pd->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
     {
  /* initialize some arrays */
  memset( TT, 0, sizeof(double)*DIM*DIM);
  if (af->Assemble_Jacobian)
   {
    memset( dTT_dx,         0, sizeof(double)*DIM*DIM*DIM*MDE);
    memset( dTT_dp,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_drs,        0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dc,         0, sizeof(double)*DIM*DIM*MAX_CONC*MDE);
    memset( dTT_dp_liq,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dp_gas,     0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dporosity,  0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dsink_mass, 0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dT,         0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dmax_strain,0, sizeof(double)*DIM*DIM*MDE);
    memset( dTT_dcur_strain,0, sizeof(double)*DIM*DIM*MDE);
   }

      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if (err == 2) exit(-1);
                /*
                 * Total mesh stress tensor...
                 */
      err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc,
			       dTT_dp_liq, dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dmax_strain, dTT_dcur_strain,
			       dTT_dT,  elc->lame_mu, elc->lame_lambda,
			       tran->delta_t, ei->ielem, 0, 1);
       /* For LINEAR ELASTICITY */
         if (cr->MeshFluxModel == LINEAR)
            {
              if (dim == 2){
                            TT[2][2] = 1.;
                            TT[1][2] = 0.;
                            TT[0][2] = 0.;
                          }
             }

       /*  For Hookian Elasticity and shrinkage */
         else
             {
              if (dim == 2){
                     elast_modulus = elc->lame_mu;
                     if (cr->MeshMotion == ARBITRARY)
                       {
                        TT[2][2] = (1. - fv->volume_change) * elast_modulus;
                       }
                     else
                       {
                       if (cr->MeshFluxModel == NONLINEAR ||
                           cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
                           cr->MeshFluxModel == INCOMP_PSTRAIN )
       TT[2][2] = (1. - pow(fv->volume_change,2./3.)) * elast_modulus - fv->P;
       /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
                       else  TT[2][2] = 0.;
                        }
                          TT[1][2] = 0.;
                          TT[0][2] = 0.;
                         }
              }
         for( a=0; a<wim; a++ )
	     {
              for( i=0; i<wim; i++ )
		{ force[a] = -pressure*TT[a][i]*lsi->normal[i];}
	     }

  /*  Jacobian contributions            */
  if (af->Assemble_Jacobian)
        {
        for(b=0; b<dim; b++)
          {
            var = MESH_DISPLACEMENT1 + b;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dx[p][q][b][j];
                            }
                        }
                  }
              }
          }
/* For Lagrangian Mesh, add in the sensitivity to pressure  */
	      if (pd->MeshMotion == LAGRANGIAN ||
		  pd->MeshMotion == DYNAMIC_LAGRANGIAN)
		{
		  var = PRESSURE;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dp[p][q][j];
                            }
                        }
                  }
              }
		}
      var = MASS_FRACTION;
      if (pd->v[var] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for ( w=0; w<pd->Num_Species_Eqn; w++)
                    {
                      for (j=0; j<ei->dof[var]; j++)
                        {
      d_force[p][MAX_VARIABLE_TYPES + w][j] -= pressure*lsi->normal[q]*dTT_dc[p][q][w][j];
                        }
                    }
                }
            }
        }
		if (pd->v[POR_LIQ_PRES] &&  
		    (cr->MeshMotion == LAGRANGIAN ||
		     cr->MeshMotion == DYNAMIC_LAGRANGIAN) )
		{
		  var = POR_LIQ_PRES;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dp_liq[p][q][j];
                            }
                        }
                  }
              }
		}
	      if (pd->v[POR_GAS_PRES] &&
		      (cr->MeshMotion == LAGRANGIAN ||
			   cr->MeshMotion == DYNAMIC_LAGRANGIAN) )
		{
	      var = POR_GAS_PRES;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dp_gas[p][q][j];
                            }
                        }
                  }
              }
		}
	      var = POR_POROSITY;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dporosity[p][q][j];
                            }
                        }
                  }
              }

                  var = POR_SINK_MASS;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dsink_mass[p][q][j];
                            }
                        }
                  }
              }
 	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&  pd->e[eqn]  ) 
		{
		  var = TEMPERATURE;
            if (pd->v[var])
              {
                for (j = 0; j < ei->dof[var]; j++)
                  {
                    for(p=0; p<dim; p++)
                        {
                        for(q=0; q<dim; q++)
                            {
                          d_force[p][var][j] -= pressure*lsi->normal[q]*dTT_dT[p][q][j];
                            }
                        }
                  }
              }
                }

 	      if ((cr->MeshMotion == LAGRANGIAN ||
		   cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&  pd->e[eqn]  ) 
		{
		  var = MAX_STRAIN;
		  if (pd->v[var])
		    {
		      for (j = 0; j < ei->dof[var]; j++)
			{
			  for(p=0; p<dim; p++)
			    {
			      for(q=0; q<dim; q++)
				{
				  d_force[p][var][j] += pressure*lsi->normal[q]*dTT_dmax_strain[p][q][j];
				}
			    }
			}
		    }
                }
	      
      var = ls->var;
      if (pd->v[var] )
        {
          for (p=0; p<pd->Num_Dim; p++)
            {
              for (q=0; q<pd->Num_Dim; q++)
                {
                  for (j=0; j<ei->dof[var]; j++)
                    {
                      d_force[p][var][j] -= pressure*
				(lsi->d_normal_dF[p][j]*TT[p][q]);
                    }
                }
            }
        }
        }  /*  if Jacobian  */

     } /* end of STRESS_TENSOR */
    else
	{EH(-1, "Invalid LS_FLOW_PRESSURE boundary condition flag (-1|0|1)\n");}

#if 0
fprintf(stderr,"pf %g %g %g %d %g %g\n",fv->x[0],fv->x[1], lsi->delta, elem_sign, force[0],force[1]);
#endif

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = phi_i * force[a] * lsi->delta * elem_sign;

		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  /* J_m_F
		   */
                  var = ls->var ;
		  pvar = upd->vp[var];
		      
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if (af->Assemble_Residual )
    {
      for (a = 0; a < wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {
	      ledof = ei->lvdof_to_ledof[eqn][i];
	      if (ei->active_interp_ledof[ledof])
		{
		  ii = ei->lvdof_to_row_lvdof[eqn][i];
		  phi_i = bfm->phi[i];
		  source = phi_i * force[a] * lsi->delta * elem_sign;
		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];
		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }

  if (af->Assemble_Jacobian)
    {
      for (a = 0; a < wim; a++)
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              ledof = ei->lvdof_to_ledof[eqn][i];

              if (ei->active_interp_ledof[ledof])
                {
                  ii = ei->lvdof_to_row_lvdof[eqn][i];
                  phi_i = bfm->phi[i];

#ifdef COUPLED_FILL
		  /*
		   * J_m_F
		   */
                  var = ls->var;
                  if (pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for (j = 0; j < ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * elem_sign *
			    ( lsi->delta * d_force[a][var][j] +
			      lsi->d_delta_dF[j] * force[a] );

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
		  /* extra dependencies for bcflag = 0 or 1	*/
		  if( bcflag == 0 || bcflag ==1)
		    {
		      /*
		       * J_m_Pressure if bc_Int == 0
		       */

		      var = PRESSURE;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = TEMPERATURE;
		      if ( pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = BOND_EVOLUTION;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for(j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = MASS_FRACTION;
		      if (pd->v[var])
			{
			  for ( w=0; w<pd->Num_Species_Eqn; w++)
			    {
			      pvar = MAX_PROB_VAR + w;
			      for (j = 0; j < ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];

				  source = phi_i * d_force[a][MAX_VARIABLE_TYPES+w][j] * elem_sign *
				    lsi->delta;

				  source *= det_J * wt * h3;
				  source *= pd->etm[eqn][LOG2_SOURCE];

				  lec->J[peqn][pvar][ii][j] += source;
				}
			    }
			}

		      var = POR_LIQ_PRES;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = POR_GAS_PRES;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = POR_POROSITY;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      var = POR_SINK_MASS;
		      if (pd->v[var])
			{
			  pvar = upd->vp[var];

			  for (j = 0; j < ei->dof[var]; j++)
			    {
			      phi_j = bf[var]->phi[j];

			      source = phi_i * d_force[a][var][j] * elem_sign *
				lsi->delta;

			      source *= det_J * wt * h3;
			      source *= pd->etm[eqn][LOG2_SOURCE];

			      lec->J[peqn][pvar][ii][j] += source;
			    }
			}

		      for (b = 0; b < dim; b++)
			{
			  var = MESH_DISPLACEMENT1 + b;

			  if (pd->v[var])
			    {
			      pvar = upd->vp[var];

			      for(j = 0; j < ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];

				  source = phi_i * d_force[a][var][j] * elem_sign *
				    lsi->delta;

				  source *= det_J * wt * h3;
				  source *= pd->etm[eqn][LOG2_SOURCE];

				  lec->J[peqn][pvar][ii][j] += source;
				}
			    }
			}
		      for (b = 0; b < dim; b++)
			{
			  var = VELOCITY1 + b;

			  if (pd->v[var])
			    {
			      pvar = upd->vp[var];

			      for (j = 0; j < ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];

				  source = phi_i * d_force[a][var][j] * elem_sign *
				    lsi->delta;

				  source *= det_J * wt * h3;
				  source *= pd->etm[eqn][LOG2_SOURCE];

				  lec->J[peqn][pvar][ii][j] += source;
				}
			    }
			}

		      for (b = 0; b < dim; b++)
			{
			  var = VORT_DIR1 + b;

			  if (pd->v[var])
			    {
			      pvar = upd->vp[var];

			      for (j = 0; j < ei->dof[var]; j++)
				{
				  phi_j = bf[var]->phi[j];

				  source = phi_i * d_force[a][var][j] * elem_sign *
				    lsi->delta;

				  source *= det_J * wt * h3;
				  source *= pd->etm[eqn][LOG2_SOURCE];

				  lec->J[peqn][pvar][ii][j] += source;
				}
			    }
			}
		    }
                }
            }
        }
    }
  return ( 1 );
}

int
assemble_precoil_source (const double p[])
{
  int i, j, a, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double pressure;
  double source;

  double dprtemp, theta;
  double pabl_c0 = 0.0, pabl_c1 = 0.0, pabl_c2 = 0.0, pabl_c3 = 0.0;
  double T_boil, T_scale, P_scale;

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  /* Mike Kanouff's curve-fit for ablation pressure in Pascals               */
  /*                                              assume iron if T_boil>2000 */
  /*                                                      ice if T_boil<2000 */
  T_boil  = mp->melting_point_solidus;  /* T_boil stored as T_solidus in mat-file */
  P_scale = p[5];
  T_scale = p[6];

  if( T_boil > 2000.0*T_scale)
    {
      theta =  fv->T-T_boil;
      if(theta > 0.0 && theta <= 170.0*T_scale )
	{
	  pabl_c0 = 0.0;                                 /* iron coefficients */
	  pabl_c1 = 1.8272e-4*1.0133e5*(1.0/T_scale);
	  pabl_c2 =-1.9436e-6*1.0133e5*(1.0/T_scale)*(1.0/T_scale);
	  pabl_c3 = 1.5732e-8*1.0133e5*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);
	}
      if(theta > 170.0*T_scale)
	{
	  pabl_c0 = 0.0;                                 /* iron coefficients */
	  pabl_c1 =-5.7333e-4*1.0133e5*(1.0/T_scale);
	  pabl_c2 = 4.5500e-6*1.0133e5*(1.0/T_scale)*(1.0/T_scale);
	  pabl_c3 = 2.3022e-9*1.0133e5*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);
	}
      /*pabl_c0 =  0.0;          */               /* MKS coefficients for iron */
      /*pabl_c1 =  3.723086e+02*(1.0/T_scale);*/
      /*pabl_c2 =  -6.328050e-02*(1.0/T_scale)*(1.0/T_scale);*/
      /*pabl_c3 =  5.559470e-04*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);*/
    }
  else
    {
      pabl_c0 =  0.0;                          /* MKS coefficients for ice  */
      pabl_c1 =  3.294180e+03*(1.0/T_scale);
      pabl_c2 =  -7.726940e+00*(1.0/T_scale)*(1.0/T_scale);
      pabl_c3 =  5.480973e-01*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);
    }

 /* Calculate ablation pressure */
     theta =  fv->T-T_boil;

   if( theta < 0.)
     {
       theta = 0.;
       pressure = 0.;
       dprtemp = 0;
     }
   else
     {
       pressure = P_scale*(pabl_c0 + pabl_c1 * theta +
                             pabl_c2 * pow(theta,2.0)  +
                             pabl_c3 * pow(theta, 3.0));
       dprtemp =P_scale*( pabl_c1 + 2. * pabl_c2 * theta +
                                             3. * pabl_c3 * pow(theta,2.0));
     }
   /* END Kanouff Curve fir for pressure */

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = phi_i * lsi->normal[a] * pressure * lsi->delta * ls->Elem_Sign;

		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  /* J_m_F
		   */
                  var = LS ;
		  pvar = upd->vp[var];
		      
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = phi_i * lsi->normal[a] * pressure * lsi->delta * ls->Elem_Sign;

		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }

  if ( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              ledof = ei->lvdof_to_ledof[eqn][i];

              if (ei->active_interp_ledof[ledof])
                {
                  ii = ei->lvdof_to_row_lvdof[eqn][i];

                  phi_i = bfm->phi[i];

                   /*
                    * J_m_T
                    */
                  var = TEMPERATURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * lsi->normal[a] * lsi->delta * ls->Elem_Sign *
                                   dprtemp * phi_j;

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#ifdef COUPLED_FILL
                   /*
                    * J_m_F
                    */
                  var = FILL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * pressure * ls->Elem_Sign *
                                   ( lsi->delta * lsi->d_normal_dF[a][j] +
                                     lsi->d_delta_dF[j] * lsi->normal[a] );

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                }
            }
        }
    } 
  return ( 1 );
}

#if 0
int
assemble_uvw_source ( int eqn, double val )
{
  int i,j, a, ii,ledof;
  int peqn, var, pvar;
  int dim, wim, b, p;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;

  int xfem_active, extended_dof, base_interp, base_dof;
  double F_i;
  int sign = 1;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if ( pd->i[eqn] != I_Q1_XV &&
       pd->i[eqn] != I_Q2_XV &&
       pd->i[eqn] != I_Q1_GP  &&
       pd->i[eqn] != I_Q2_GP  &&
       pd->i[eqn] != I_Q1_GN  &&
       pd->i[eqn] != I_Q2_GN  &&
       pd->i[eqn] != I_Q1_G  &&
       pd->i[eqn] != I_Q2_G )
    {
      EH(-1, "LS_UVW requires discontinuous velocity enrichment (Q1_XV, Q2_XV, Q1_G, Q2_G)");
    }

  a = eqn - R_MOMENTUM1;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  /*
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  /* compute stress tensor and its derivatives */
  fluid_stress( Pi, d_Pi );

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      F_i = dof_distance( eqn, i );

	      if ( !extended_dof && !sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += Pi[p][a] * lsi->normal[p];
                  source *= phi_i * lsi->delta * ls->Elem_Sign;
                  source *= -det_J * wt * h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
	      if ( extended_dof && sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  source = phi_i * lsi->delta * (fv->v[a] - val);
                  source *= det_J * wt * h3 * BIG_PENALTY * sign;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      F_i = dof_distance( eqn, i );

	      if ( !extended_dof && !sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  /*
                   * J_m_m
                   */
                  for ( b=0; b<wim; b++)
                    {
                      var = VELOCITY1 + b;

                      if ( pd->v[var] )
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              phi_j = bf[var]->phi[j];

                              source = 0.;
                              for ( p=0; p<VIM; p++) source += d_Pi->v[b][p][a][j] * lsi->normal[p];
                              source *= phi_i * lsi->delta * ls->Elem_Sign;
                              source *= -det_J * wt * h3;
                              source *= pd->etm[eqn][LOG2_SOURCE];

                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }
                    }
                  /*
                   * J_m_P
                   */
                  var = PRESSURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->P[p][a][j] * lsi->normal[p];
                          source *= phi_i * lsi->delta * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
#ifdef COUPLED_FILL
                  var = FILL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->F[p][a][j] * lsi->normal[p] * lsi->delta +
                                                           Pi[p][a] * lsi->d_normal_dF[p][j] * lsi->delta +
                                                           Pi[p][a] * lsi->normal[p] * lsi->d_delta_dF[j];
                          
                          source *= phi_i * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                }

	      if ( extended_dof && sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  /*
                   * J_m_m
                   */
                  var = eqn;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * phi_j * lsi->delta;
                          source *= det_J * wt * h3 * BIG_PENALTY * sign;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
                  var = FILL;
                  if( pd->v[var])
                    {
                      pvar = upd->vp[var];
	      	      
	      	      /* path dependence integral */
                      if ( ls->on_sharp_surf && !ls->AdaptIntegration)
                	{
                	  for( j=0; j<ei->dof[var]; j++)
                	    {
                	      phi_j = bf[var]->phi[j];

                	      source = 0.;
/*
                	      for ( p=0; p < VIM; p++ )
                		{
                		  source -= phi_i * phi_j * lsi->normal[p] * fv->grad_v[p][a];
                		}
*/

                	      for ( p=0; p < VIM; p++ )
                		{
                		  source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * 
	      	 			   ( fv->v[a] - val );
                		}

                	      source *= lsi->gfmaginv;
                	      source *= det_J * wt * h3 * BIG_PENALTY * sign;
                	      source *= pd->etm[eqn][(LOG2_SOURCE)];

                	      lec->J[peqn][pvar][ii][j] += source;
                	    }
                	}
                      else
                	{
                	  /* diffuse interface version of path dependence integral */
                	  for( j=0; j<ei->dof[var]; j++)
                	    {
                	      source = phi_i * ( fv->v[a] - val );
                	      source *= lsi->d_delta_dF[j];
                	      source *= det_J * wt * h3 * BIG_PENALTY;
                	      source *= pd->etm[eqn][(LOG2_SOURCE)];

                	      lec->J[peqn][pvar][ii][j] += source;
                	    }
                	}
	      	    }
                }
            }
        }
    }
  return ( 1 );
}
#endif
#if 1
int
assemble_uvw_source ( int eqn, double val )
{
  int i,j, a, ii,ledof;
  int peqn, var, pvar;
  int dim, wim, b, p;
  int xfem_active, extended_dof, base_interp, base_dof;
  int incomplete, other_side;
  int apply_NOBC[MDE], apply_SIC[MDE];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double F_i;
  int sign = 1;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if ( pd->i[eqn] != I_Q1_XV &&
       pd->i[eqn] != I_Q2_XV &&
       pd->i[eqn] != I_Q1_XG &&
       pd->i[eqn] != I_Q2_XG &&
       pd->i[eqn] != I_Q1_GP  &&
       pd->i[eqn] != I_Q2_GP  &&
       pd->i[eqn] != I_Q1_GN  &&
       pd->i[eqn] != I_Q2_GN  &&
       pd->i[eqn] != I_Q1_G  &&
       pd->i[eqn] != I_Q2_G )
    {
      EH(-1, "LS_UVW requires discontinuous velocity enrichment");
    }

  a = eqn - R_MOMENTUM1;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  /*
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  if ( pd->i[eqn] == I_Q1_XG ||
       pd->i[eqn] == I_Q2_XG )
    {
      for (i = 0; i < ei->dof[eqn]; i++)
        {
	  apply_NOBC[i] = FALSE;
          apply_SIC[i] = TRUE;
	}
    }
  else
    {
      /* compute stress tensor and its derivatives */
      fluid_stress( Pi, d_Pi );
  
      /* determine if NOBC or SIC or both are to be applied */
      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      incomplete = dof_incomplete( base_dof, ei->ielem_type, base_interp, ei->ielem_shape);
	      F_i = dof_distance( eqn, i );
	      other_side = sign_change( F_i, (double) ls->Elem_Sign );
	      
	      if ( extended_dof ) {
	        if ( other_side ) {
	          if ( incomplete ) {
		    /* extended_dof&&other_side&&incomplete -> NOBC + SIC */
		    apply_NOBC[i] = TRUE;
		    apply_SIC[i] = TRUE;
	          } else {
		    /* extended_dof&&other_side&&!incomplete -> NOBC */
		    apply_NOBC[i] = TRUE;
		    apply_SIC[i] = FALSE;
	          }
	        } else {
	          /* extended_dof&&!other_side -> do nothing */
	          apply_NOBC[i] = FALSE;
	          apply_SIC[i] = FALSE;
	        }
	      } else {
	        if ( other_side ) {
	          /* !extended_dof&&other_side -> do nothing */
	          apply_NOBC[i] = FALSE;
	          apply_SIC[i] = FALSE;
	        } else {
	          /* !extended_dof&&!other_side -> NOBC */
	          apply_NOBC[i] = TRUE;
	          apply_SIC[i] = FALSE;
		}
	      }
	    }
        }
    }

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

	      if ( apply_NOBC[i] )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += Pi[p][a] * lsi->normal[p];
                  source *= phi_i * lsi->delta * ls->Elem_Sign;
                  source *= -det_J * wt * h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
	            }
                }
	      if ( apply_SIC[i] )
                {
                  source = phi_i * lsi->delta * (fv->v[a] - val);
                  source *= det_J * wt * h3 * BIG_PENALTY * sign;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  /* J_m_F
	           */
                  var = LS ;
	          pvar = upd->vp[var];
		      
	          for( j=0; j<ei->dof[var]; j++)
	            {
	              phi_j = bf[var]->phi[j];
	              lec->J[peqn][pvar][ii][j] += source * phi_j;
	            }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

	      if ( apply_NOBC[i] )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += Pi[p][a] * lsi->normal[p];
                  source *= phi_i * lsi->delta * ls->Elem_Sign;
                  source *= -det_J * wt * h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
	      if ( apply_SIC[i] )
                {
                  source = phi_i * lsi->delta * (fv->v[a] - val);
                  source *= det_J * wt * h3 * BIG_PENALTY * sign;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];

                  lec->R[peqn][ii] += source;
                }
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {
          ledof = ei->lvdof_to_ledof[eqn][i];

          if (ei->active_interp_ledof[ledof])
            {
              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      F_i = dof_distance( eqn, i );

	      if ( apply_NOBC[i] )
                {
                  /*
                   * J_m_v
                   */
                  for ( b=0; b<wim; b++)
                    {
                      var = VELOCITY1 + b;

                      if ( pd->v[var] )
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              phi_j = bf[var]->phi[j];

                              source = 0.;
                              for ( p=0; p<VIM; p++) source += d_Pi->v[p][a][b][j]* lsi->normal[p];
                              source *= phi_i * lsi->delta * ls->Elem_Sign;
                              source *= -det_J * wt * h3;
                              source *= pd->etm[eqn][LOG2_SOURCE];

                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }
                    }
                  /*
                   * J_m_vd
                   */
                  for ( b=0; b<DIM; b++)
                    {
                      var = VORT_DIR1 + b;

                      if ( pd->v[var] )
                        {
                          pvar = upd->vp[var];

                          for( j=0; j<ei->dof[var]; j++)
                            {
                              phi_j = bf[var]->phi[j];

                              source = 0.;
                              for ( p=0; p<VIM; p++) source += d_Pi->vd[p][a][b][j] * lsi->normal[p];
                              source *= phi_i * lsi->delta * ls->Elem_Sign;
                              source *= -det_J * wt * h3;
                              source *= pd->etm[eqn][LOG2_SOURCE];

                              lec->J[peqn][pvar][ii][j] += source;
                            }
                        }
                    }
                  /*
                   * J_m_P
                   */
                  var = PRESSURE;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->P[p][a][j] * lsi->normal[p];
                          source *= phi_i * lsi->delta * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
#ifdef COUPLED_FILL
                  var = FILL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->F[p][a][j] * lsi->normal[p] * lsi->delta +
                                                           Pi[p][a] * lsi->d_normal_dF[p][j] * lsi->delta +
                                                           Pi[p][a] * lsi->normal[p] * lsi->d_delta_dF[j];
                          
                          source *= phi_i * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                }

	      if ( apply_SIC[i] )
                {
                  /*
                   * J_m_m
                   */
                  var = eqn;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * phi_j * lsi->delta;
                          source *= det_J * wt * h3 * BIG_PENALTY * sign;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
                  var = FILL;
                  if( pd->v[var])
                    {
                      pvar = upd->vp[var];
	      	      
		      for( j=0; j<ei->dof[var]; j++)
		        {
		          source = phi_i * ( fv->v[a] - val );
		          source *= lsi->d_delta_dF[j];
		          source *= det_J * wt * h3 * BIG_PENALTY;
		          source *= pd->etm[eqn][(LOG2_SOURCE)];

		          lec->J[peqn][pvar][ii][j] += source;
                	}
	      	    }
                }
            }
        }
    }
  return ( 1 );
}
#endif

int
assemble_extension_velocity_path_dependence(void)
{
  int a, i, j, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i;

  double source = 0.;
  double resid;
  double *grad_F;

  eqn = R_EXT_VELOCITY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;
  h3 = fv->h3;
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }
  if(pfd == NULL)
        {
        grad_F = fv->grad_F;
        }
  else
        {
        grad_F = fv->grad_pF[ls->var-PHASE1];
        }

#define GRADF_GRADEXTV 1
#define NORMAL_GRADEXTV 0
#if GRADF_GRADEXTV
  resid = 0.;
  for ( a=0; a<dim; a++) resid += grad_F[a] * fv->grad_ext_v[a];
#endif
#if NORMAL_GRADEXTV
  resid = 0.;
  for ( a=0; a<VIM; a++) resid += lsi->normal[a] * fv->grad_ext_v[a];
#endif
  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if ( ei->active_interp_ledof[ledof] )
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              /* path dependence integral from volume equation */

              /*
               * J_ls_F
               */
              var = ls->var;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {

                      source = resid * (2. * lsi->d_H_dF[j]) * phi_i;
                      source *= det_J * wt * h3;

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
            }
        }
    }

  return ( 1 );
}

int
assemble_fill_path_dependence ( void )
{
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i;

  double source = 0.;


  eqn = R_FILL;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;
  h3 = fv->h3;
  if ( ls->on_sharp_surf ) /* sharp interface */
    {
      det_J = fv->sdet;
    }
  else              /* diffuse interface */
    {
      det_J = bf[eqn]->detJ;
    }
    
  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei->dof[eqn]; i++)
        {

          ledof = ei->lvdof_to_ledof[eqn][i];

          if ( ei->active_interp_ledof[ledof] )
            {

              ii = ei->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              /* path dependence integral from volume equation */

              /*
               * J_ls_F
               */
              var = FILL;

              if ( pd->v[var] )
                {
                  pvar = upd->vp[var];

                  for( j=0; j<ei->dof[var]; j++)
                    {
                      switch ( tran->Fill_Weight_Fcn ) {
                      case FILL_WEIGHT_G:
                      case FILL_WEIGHT_GLS:
		      
			source = (lsi->gfmag - 1.) * (2. * lsi->d_H_dF[j]) * phi_i;

                        break;
                      default:
                        EH(-1,"Unknown Fill_Weight_Fcn");
                      }

                      source *= det_J * wt * h3;
                      source *= pd->etm[eqn][(LOG2_SOURCE)];

                      lec->J[peqn][pvar][ii][j] += source;
                    }
                }
            }
        }
    }

  return ( 1 );
}

int
assemble_energy_path_dependence(
	        double time,	/* present time value */
		double tt,	/* parameter to vary time integration from
				   explicit (tt = 1) to implicit (tt = 0) */
		double dt,	/* current time step size */
		const PG_DATA *pg_data ) /* average element and velocity infor for SUPG and PSPG */
{
  int eqn, var, peqn, pvar, dim, p, i, j, status;

  dbl T_dot;				/* Temperature derivative wrt time. */

  dbl q[DIM];				/* Heat flux vector. */
  dbl rho;				/* Density */
  dbl Cp;				/* Heat capacity. */
  dbl h;				/* Heat source. */

  dbl mass;		         	/* For terms and their derivatives */
  dbl advection;			/* For terms and their derivatives */
  dbl diffusion;
  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

/* SUPG variables */
  dbl h_elem=0;
  dbl supg;
  const dbl *vcent, *hsquared;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3;			/* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  dbl *grad_T = fv->grad_T;

  double vconv[MAX_PDIM]; /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/

  int err;
  
  int sign = ls->Elem_Sign;
  double energy_residual;

  /*   static char yo[] = "assemble_energy";*/

  status = 0;
  
  eqn   = R_ENERGY;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] )
    {
      return(status);
    }
    
  dim   = pd->Num_Dim;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  supg = 0.0;

  if( mp->Ewt_funcModel == GALERKIN)
    {
      supg = 0.;
    }
  else if( mp->Ewt_funcModel == SUPG)
    {
      supg = mp->Ewt_func;
    }

  if (supg != 0.0)
    {
      h_elem = 0.;
      vcent = pg_data->v_avg;
      hsquared = pg_data->hsquared;

      for ( p=0; p<dim; p++)
	{
	  h_elem += vcent[p] * vcent[p] * hsquared[p];
	}
      h_elem = sqrt(h_elem)/2.;
    }
/* end Petrov-Galerkin addition */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are rho, Cp, k, and h.  For now we will
   * take rho as constant.   Cp, h, and k we will allow to vary with temperature,
   * spatial coordinates, and species concentration.
   */

  rho = density(NULL, time);

  /* CHECK FOR REMOVAL */
  conductivity( NULL, time );

  Cp = heat_capacity( NULL, time );

  h = heat_source( NULL, time, tt, dt );

  if (pd->TimeIntegration != STEADY) {
    T_dot = fv_dot->T;
  } else {
    T_dot = 0.0;
  }

  heat_flux( q, NULL, time );

/* get the convection velocity (it's different for arbitrary and
   lagrangian meshes) */
  if( cr->MeshMotion == ARBITRARY ||
      cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN)
    {
      err = get_convection_velocity(vconv, vconv_old, NULL, dt, tt);
      EH(err, "Error in calculating effective convection velocity");
    }
  else if ( cr->MeshMotion == TOTAL_ALE)
    {
      err = get_convection_velocity_rs(vconv, vconv_old, NULL, dt, tt);
      EH(err, "Error in calculating effective convection velocity_rs");
    }


  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Jacobian )
    {
      eqn = R_ENERGY;
      peqn = upd->ep[eqn];
      var = TEMPERATURE;
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  phi_i = bf[eqn]->phi[i];

	  mass = 0.;
	  if ( pd->TimeIntegration != STEADY )
	    {
	      if ( pd->e[eqn] & T_MASS )
		{
		  mass  = T_dot;
		  mass *= - phi_i * rho * Cp * det_J * wt;
		  mass *= h3;
		  mass *= pd->etm[eqn][(LOG2_MASS)];
		}
	    }

	  /* only use Petrov Galerkin on advective term - if required */
	  wt_func = bf[eqn]->phi[i];
	  /* add Petrov-Galerkin terms as necessary */
	  if(supg!=0.)
	    {
	      for(p=0; p<dim; p++)
		{
		  wt_func += supg * h_elem * vconv[p] * bf[eqn]->grad_phi[i][p];
		}
	    }

	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {

	      for ( p=0; p<VIM; p++)
		{
		  advection += vconv[p]     * grad_T[p];
		}

	      advection *= - wt_func * rho * Cp * det_J * wt;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }

	  diffusion = 0.;
	  if ( pd->e[eqn] & T_DIFFUSION )
	    {
	      for ( p=0; p<VIM; p++)
		{
		  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
		}

	      for ( p=0; p<VIM; p++)
		{
		  diffusion += grad_phi_i[p] * q[p];
		}
	      diffusion *= det_J * wt;
	      diffusion *= h3;
	      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  source = 0.;
	  if ( pd->e[eqn] & T_SOURCE )
	    {
	      source += phi_i * h * det_J * wt;
	      source *= h3;
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	    }

	  energy_residual = mass + advection +  diffusion + source;
	  
	  var = FILL;
          pvar = upd->vp[var];
          for ( j=0; j<ei->dof[var]; j++ )
            {
	      lec->J[peqn][pvar][i][j] += lsi->d_H_dF[j] * energy_residual * sign;
            }
	}
    }
  return(status);
}

int
assemble_momentum_path_dependence(dbl time,       /* currentt time step */
				  dbl tt,	  /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0) */
				  dbl dt,	  /* current time step size */
				  const PG_DATA *pg_data)
{
  int dim, wim;
  int i, j, a, p;
  int ledof, eqn, var, ii, peqn, pvar;
  int status;
  struct Basis_Functions *bfm;

  dbl zero[3] = {0.0, 0.0, 0.0};        /* A zero array, for convenience. */
  dbl *v_dot;		         	/* time derivative of velocity field. */
  dbl *x_dot;			        /* current position field derivative wrt time. */

  dbl h3;			        /* Volume element (scale factors). */
  dbl det_J;                            /* determinant of element Jacobian */

  /* field variables */
  dbl *grad_v[DIM];
  dbl *v = fv->v;

  dbl rho;				/* Density. */

  dbl f[DIM];				/* Body force. */

  dbl mass;			        /* For terms and their derivatives */
  dbl advection;
  dbl porous;
  dbl diffusion;
  dbl source;

  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl (* grad_phi_i_e_a ) [DIM] = NULL;
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl Pi[DIM][DIM];

  dbl wt;
  
  dbl d_area;

/* coefficient variables for the Brinkman Equation
   in flows through porous media: KSC on 5/10/95 */
  dbl por;         /* porosity of porous media */
  dbl por2;        /* square of porosity */
  dbl per = 0.0;         /* permeability of porous media */
  /* derivative of permeability wrt concentration */
  dbl d_per_dc[MAX_CONC][MDE];

  dbl vis;         /* flowing-liquid viscosity */
  dbl dvis_dT[MDE];/* sensitivities of flowing-liquid viscosity */
  dbl sc;          /* inertial coefficient */
  dbl speed;       /* magnitude of the velocity vector */

  int v_s[MAX_MODES][DIM][DIM];

  int sign;
  double momentum_residual;

  /* the residual is essentially calculated from
     R = R+ * H[F] + R- * (1-H[F]), so
     dR/dFj = delta[F] * Nj * (R+ - R-)

     So, here we are to assemble for one side of the interface,
     dR/Fj += delta[F] * Nj * sign * R[sign]
   */

  /* Variables used for the modified fluid momentum equations when looking
   * at the particle momentum model.  Note that pd->MomentumFluxModel ==
   * SUSPENSION_PM when the Buyevich particle momentum equations are active.
   */
  int particle_momentum_on;	/* Boolean for particle momentum eq.'s active */
  int species=-1;		/* Species number of particle phase */
  double p_vol_frac=0;		/* Particle volume fraction (phi) */
  double ompvf=1;		/* 1 - p_vol_frac "One Minus Particle */
				/*    Volume Fraction" */
  int mass_on;
  int advection_on;
  int diffusion_on;
  int source_on;
  int porous_brinkman_on;
  int transient_run = (pd->TimeIntegration != STEADY);
  
  int *pde = pd->e; 
  dbl mass_etm ;
  dbl advection_etm;
  dbl diffusion_etm ;
  dbl porous_brinkman_etm;
  dbl source_etm;

  dbl h_elem_avg;

  /*Continuity stabilization*/
  dbl cont_gls, continuity_stabilization;


  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_MOMENTUM1;
  var = FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if ( ! pd->e[eqn] || !pd->v[var] )
    {
      return(status);
    }
	
	
  sign = ls->Elem_Sign;
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }
  
  d_area = wt*h3*det_J;

  if( pd->v[POLYMER_STRESS11] )
  {
  (void)stress_eqn_pointer(v_s);
  }

  /* Set up variables for particle/fluid momentum coupling.
   */
  if(pd->e[R_PMOMENTUM1])
    {
      particle_momentum_on = 1;
      /* This is the species number of the particle phase. */
      species = (int) mp->u_density[0];
      p_vol_frac = fv->c[species];
      ompvf = 1.0 - p_vol_frac;
      /* Uncomment this to check for when the particle volume fraction
       * becomes non-physical.  Beware, however, that the intermediate
       * solutions may, indeed, become negative while converging to a
       * physical solution.
      if(p_vol_frac<0.0 || p_vol_frac>1.0)
	{
	  if(fabs(p_vol_frac)<1e-14)
	    p_vol_frac=0.0;
	  else
	    {
	      printf("assemble_momentum: p_vol_frac=%g, exiting\n",p_vol_frac);
	      exit(0);
	    }
	}
      */
    }
  else
    particle_momentum_on = 0;

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */
  h_elem_avg = pg_data->h_elem_avg;

  /*** Density ***/

  rho = density(NULL, time);

  if ( pd->e[eqn] & T_POROUS_BRINK )
    {
      if (mp->PorousMediaType != POROUS_BRINKMAN)
	WH(-1, "Set Porous term multiplier in continuous medium");
      if(mp->FlowingLiquidViscosityModel == CONSTANT)
	{
	  /* Do nothing */
	}
      else if (mp->FlowingLiquidViscosityModel == MOLTEN_GLASS)
	{
	  (void) molten_glass_viscosity(&(mp->FlowingLiquid_viscosity),
					dvis_dT, mp->u_FlowingLiquid_viscosity);
	}
      else if (mp->FlowingLiquidViscosityModel == USER)
	{
	  (void) usr_FlowingLiquidViscosity(mp->u_FlowingLiquid_viscosity);
	  var = TEMPERATURE;

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      dvis_dT[j]= mp->d_FlowingLiquid_viscosity[var]*bf[var]->phi[j];
	    }
	}
      else
	{
	  EH(-1,"Don't recognize your FlowingLiquidViscosity model");
	}

      if (mp->PermeabilityModel == SOLIDIFICATION)
	{
	  /* This is a permeability model that slows down
	   * the flow with increasing solid fraction.
	   * It should be useful for solidification problems.
	   */
	  per = solidification_permeability(h_elem_avg,d_per_dc);
	}
      else if (mp->PermeabilityModel == CONSTANT)
	{
	  per = mp->permeability;
	}
      else
	{
	  EH(-1,"Unrecognizable Permeability model");
	}

      /* Load up remaining parameters for the Brinkman Equation. */
      por = mp->porosity;
      vis = mp->FlowingLiquid_viscosity;
      /* Load variable FlowingLiquid_viscosity */
      sc  = mp->Inertia_coefficient;
    }
  else
    {
      por = 1.;
      per = 1.;
      vis = mp->viscosity;
      sc  = 0.;
    }

  eqn = R_MOMENTUM1;
  /*
   * Field variables...
   */

  if (  transient_run &&  pd->v[MESH_DISPLACEMENT1] )
    x_dot = fv_dot->x;
  else
    x_dot = zero;
  
  if (transient_run )
    v_dot = fv_dot->v;
  else
    v_dot = zero;


  /* for porous media stuff */
  speed = 0.0;
  for ( a=0; a<wim; a++)
    {
      speed += v[a]*v[a];
    }
  speed = sqrt(speed);

  for ( a=0; a<VIM; a++) grad_v[a] = fv->grad_v[a];

  /*
   * Stress tensor, but don't need dependencies
   */
  fluid_stress( Pi, NULL );

  (void) momentum_source_term(f, NULL, time);

  //Call continuity stabilization if desired
  if(Cont_GLS)
    {
      calc_cont_gls(&cont_gls, NULL, time, pg_data);
    }

  if ( af->Assemble_Jacobian )
  {
      /*
       * Assemble each component "a" of the momentum equation...
       */
      for ( a=0; a<wim; a++)
	  {
		  eqn  = R_MOMENTUM1 + a;
		  peqn = upd->ep[eqn];
		  bfm  = bf[eqn];
		  
		  mass_on = pde[eqn] & T_MASS;
		  advection_on =  pde[eqn] & T_ADVECTION ;
		  diffusion_on = pde[eqn] & T_DIFFUSION;
		  source_on = pde[eqn] & T_SOURCE;
		  porous_brinkman_on = pde[eqn] & T_POROUS_BRINK;
		  
		  mass_etm = pd->etm[eqn][(LOG2_MASS)];
		  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];
		  diffusion_etm = pd->etm[eqn][(LOG2_DIFFUSION)];
		  porous_brinkman_etm = pd->etm[eqn][(LOG2_POROUS_BRINK)];
		  source_etm = pd->etm[eqn][(LOG2_SOURCE)];
		  /*
		   * In the element, there will be contributions to this many equations
		   * based on the number of degrees of freedom...
		   */
		  
		  for (i = 0; i < ei->dof[eqn]; i++) {
			  ledof = ei->lvdof_to_ledof[eqn][i];
			  if (ei->active_interp_ledof[ledof]) {
				  /*
				   *  Here is where we figure out whether the row is to placed in
				   *  the normal spot (e.g., ii = i), or whether a boundary condition
				   *  require that the volumetric contribution be stuck in another
				   *  ldof pertaining to the same variable type.
				   */
				  ii = ei->lvdof_to_row_lvdof[eqn][i];
				  
				  phi_i = bfm->phi[i];
				  
				  mass = 0.;
				  
				  grad_phi_i_e_a = bfm->grad_phi_e[i][a];
				  
				  if ( transient_run  )
				  {
					  if ( mass_on )
					  {
						  mass = v_dot[a] * rho;
						  mass *= - phi_i * d_area;
						  mass *= mass_etm;
					  }
					  
					  if (porous_brinkman_on )
					  {
						  mass /= por;
					  }
					  
					  if(particle_momentum_on)
					  {
						  mass *= ompvf;
					  }
				  }
				  
				  advection = 0.;
				  if (advection_on  )
				  {
#ifdef DO_NO_UNROLL
					  for ( p=0; p<wim; p++)
					  {
						  advection += (v[p] - x_dot[p]) * grad_v[p][a];
					  }
					  
#else
					  advection += (v[0] - x_dot[0]) * grad_v[0][a];
					  advection += (v[1] - x_dot[1]) * grad_v[1][a];
					  if( wim == 3 ) advection += (v[2] - x_dot[2]) * grad_v[2][a];
					  
#endif
					  advection *= rho;
					  advection *= - phi_i * d_area ;
					  advection *= advection_etm;
					  
					  if ( porous_brinkman_on  )
					  {
						  por2 = por*por;
						  advection /= por2;
					  }
					  
					  if(particle_momentum_on)
					  {
						  advection *= ompvf;
					  }
				  }
				  
				  porous    = 0.;
				  if ( porous_brinkman_on  )
				  {
					  porous = v[a]*(rho*sc*speed/sqrt(per)+vis/per);
					  porous    *= - phi_i * d_area;
					  porous    *= porous_brinkman_etm;

				  }
				  
				  
				  diffusion = 0.;
				  if ( diffusion_on  )
				  {
#ifdef DO_NO_UNROLL
					  for ( p=0; p<VIM; p++)
					  {
						  for ( q=0; q<VIM; q++)
						  {
							  diffusion += grad_phi_i_e_a[p][q] * Pi[q][p];
						  }
					  }
#else
					  diffusion += grad_phi_i_e_a[0][0] * Pi[0][0];
					  diffusion += grad_phi_i_e_a[1][1] * Pi[1][1];
					  diffusion += grad_phi_i_e_a[0][1] * Pi[1][0];
					  diffusion += grad_phi_i_e_a[1][0] * Pi[0][1];
					  
					  if( VIM == 3 )
					  {
						  diffusion += grad_phi_i_e_a[2][2] * Pi[2][2];
						  diffusion += grad_phi_i_e_a[2][1] * Pi[1][2];
						  diffusion += grad_phi_i_e_a[2][0] * Pi[0][2];
						  diffusion += grad_phi_i_e_a[1][2] * Pi[2][1];
						  diffusion += grad_phi_i_e_a[0][2] * Pi[2][0];				  
					  }
					  
#endif
					  
					  diffusion *= - d_area;
					  diffusion *= diffusion_etm ;
				  }
				  
				  
				  source = 0.0;
				  if ( source_on )
				  {
					  source += f[a];
					  
					  source *= phi_i * d_area;
					  source *= source_etm;
				  }
				  
				  /* MMH For massful particles. */
				  if(Particle_Dynamics &&
					 (Particle_Model == SWIMMER_EXPLICIT || Particle_Model == SWIMMER_IMPLICIT))
				  {
					  if(a == pd->Num_Dim - 1)
						  /* These data values should hold the entire
						  * source term. */
						  source = element_particle_info[ei->ielem].source_term[i];
				  }

				  //Continuity residual
				  continuity_stabilization = 0.0;
				  if (Cont_GLS)
				    {
				      for (p=0; p<VIM; p++){
					continuity_stabilization += grad_phi_i_e_a[p][p];
				      }
				      continuity_stabilization *= cont_gls*d_area;
				    }
				  
				  
				  momentum_residual =
					  mass + advection + porous + diffusion + source + continuity_stabilization;
				  
				  var = FILL;
				  pvar = upd->vp[var];
				  for ( j=0; j<ei->dof[var]; j++ )
				  {
					  
					  lec->J[peqn][pvar][ii][j] += lsi->d_H_dF[j] * momentum_residual * sign;
				  }
				  
			  }  /*end if (active_dofs) */
		  } /* end of for (i=0,ei->dofs...) */
	  }
  }
return(status);

} /* end of function assemble_momentum_path_dependence                */
/**********************************************************************/

int
assemble_continuity_path_dependence (dbl time_value,
				     dbl tt,	/* parameter to vary time integration fromexplicit (tt = 1) to implicit (tt = 0)    */
				     dbl dt,	/* current time step size                    */
				     const PG_DATA *pg_data)
{
	int dim, wim;
	int p, a;
	int dofs;
	
	int eqn, var;
	int peqn, pvar;
	int w;
	
	int i;
	int j;
	
	int status, err;
	
	dbl time = 0.0; /*  RSL 6/6/02  */
	
	dbl *v = fv->v;			/* Velocity field. */
	dbl div_v = fv->div_v;		/* Divergence of v. */
	
	dbl epsilon=0, derivative, sum;		/*  RSL 7/24/00  */
	
	dbl advection;
	dbl source;
	dbl pressure_stabilization;
	
	dbl volsolvent=0;		/* volume fraction of solvent                */
	dbl initial_volsolvent=0;
	
	dbl det_J;
	dbl h3;
	dbl wt;
	
	
	
	/*
	 * Galerkin weighting functions...
	 */
	
	dbl phi_i;
	
	dbl ( *grad_phi )[DIM];               /* weight-function for PSPG term */
	
	
	/*
	 * Variables for Pressure Stabilization Petrov-Galerkin...
	 */
	int meqn;

	dbl pspg[DIM];
	
	dbl mass;
	
	
	dbl rho=0;
	DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
	DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
	
	struct Species_Conservation_Terms s_terms;
	dbl rhos=0, rhof=0;
	dbl h_flux=0;
	int w0=0;
	
	/* For particle momentum model.
		*/
	int species;			/* species number for particle phase,  */
	dbl ompvf;			/* 1 - partical volume fraction */
	int particle_momentum_on; 	/* boolean. */
	
	/* Foaming model TAB */
	double dFVS_dv[DIM][MDE];
	double dFVS_dT[MDE];
	double dFVS_dx[DIM][MDE];
	double dFVS_dC[MAX_CONC][MDE];
	double dFVS_dF[MDE];
	
	int sign;
	double continuity;


	status = 0;
	
	/*
	 * Unpack variables from structures for local convenience...
	 */
	
	eqn   = R_PRESSURE;
	peqn = upd->ep[eqn];
	
  /*
   * Bail out fast if there's nothing to do...
   */

	if ( ! pd->e[eqn] ||  !pd->v[FILL] )
    {
		return(status);
    }
	
	sign = ls->Elem_Sign;

	dim   = pd->Num_Dim;
	wim   = dim;
	if(pd->CoordinateSystem == SWIRLING ||
	   pd->CoordinateSystem == PROJECTED_CARTESIAN)
		wim = wim+1;
	
	wt = fv->wt;
	h3 = fv->h3;			/* Differential volume element (scales). */
	
	if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
	else
	{det_J = bf[eqn]->detJ; }
	
	grad_phi = bf[eqn]->grad_phi;
	
  /*
   * Get the deformation gradients and tensors if needed
   */

	if ((cr->MeshMotion == LAGRANGIAN ||
		 cr->MeshMotion == DYNAMIC_LAGRANGIAN) && pd->e[R_MESH1])
    {
		err = belly_flop(elc->lame_mu);
		EH(err, "error in belly flop");
		if (err == 2) return(err);
    }
	
	if ((cr->MeshMotion == TOTAL_ALE && !pd->v[VELOCITY1]) && pd->e[R_SOLID1])
    {
		err = belly_flop_rs(elc_rs->lame_mu);
		EH(err, "error in belly flop for real solid");
		if (err == 2) return(err);
    }
	
	if(pd->e[R_PMOMENTUM1])
    {
		particle_momentum_on = 1;
		species = (int) mp->u_density[0];
		ompvf = 1.0 - fv->c[species];
    }
	else
    {
		particle_momentum_on = 0;
		species = -1;
		ompvf = 1.0;
    }
	
	if(PSPG)
    {
		calc_pspg( pspg, NULL,
			   time_value, tt, dt,
			   pg_data);
    }
	

	if ((cr->MeshMotion == LAGRANGIAN ||
		 cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
		 (cr->MeshMotion == TOTAL_ALE && !pd->v[VELOCITY1]))
		&& (mp->PorousMediaType == CONTINUOUS))
    {
		initial_volsolvent = elc->Strss_fr_sol_vol_frac;
		volsolvent = 0.;
		for (w=0; w<pd->Num_Species_Eqn; w++) volsolvent += fv->c[w];
		if(particle_momentum_on)
			volsolvent -= fv->c[species];
    }
	
	if ( mp->MomentumSourceModel == SUSPENSION_PM ||
		 mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS )
	  {
	    err = get_continuous_species_terms(&s_terms, 0.0, tt, dt, pg_data->hsquared);
	    EH(err,"problem in getting the species terms");
	  }
	
	if ( mp->SpeciesSourceModel[0]  == ION_REACTIONS )
	{
		zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
		err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
		EH(err,"problem in getting the species terms");
	}
	
	if( (cr->MassFluxModel == HYDRODYNAMIC) && ( mp->DensityModel == SUSPENSION )
		&& ( mp->MomentumSourceModel == SUSPENSION ))
	{
		/*
		 * Compute hydrodynamic/sedimentation flux and sensitivities.
		 */
		
		
		w0 = (int) mp->u_density[0]; /* This is the species number that is transported HYDRODYNAMICally  */
		
		hydro_flux( &s_terms, w0, tt, dt, pg_data->hsquared);
		
		rhof = mp->u_density[1];
		rhos = mp->u_density[2];
		
	}
	
	if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
		mp->SpeciesSourceModel[0]  == ION_REACTIONS ) 
	{
	  rho  = density(d_rho, time_value);
	}
	
	if ( af->Assemble_Jacobian )
	  {
		for ( i=0; i<ei->dof[eqn]; i++)
        {
			phi_i      = bf[eqn]->phi[i];
			
			/* This terms shows up because the fluid phase is not assumed to be
			* incompressible.
			*/
			mass = 0.0;


			if( particle_momentum_on )
			{
				if( pd->TimeIntegration != STEADY )
				{
					mass = -s_terms.Y_dot[species];
					mass *= phi_i * det_J * h3 * wt;
				}
			}
			
			if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
				mp->SpeciesSourceModel[0]  == ION_REACTIONS )
			{
				mass = 0.0;
				if ( pd->TimeIntegration != STEADY )
				{
					if (mp->PorosityModel == CONSTANT)
					{
						epsilon = mp->porosity;
					}
					else if (mp->PorosityModel == THERMAL_BATTERY)
					{
						epsilon = mp->u_porosity[0];
					}
					else
					{
						EH(-1, "invalid porosity model");
					}
					var = MASS_FRACTION;
					for (w=0; w<pd->Num_Species-1; w++)
					{
						for ( j=0; j<ei->dof[var]; j++)
						{
							if ( bf[var]->phi[j] > 0.0 ) break;
						}
						derivative = d_rho->C[w][j]/bf[var]->phi[j];
						mass += derivative * s_terms.Y_dot[w];
					}
					mass *= epsilon/rho;
					mass *= phi_i * det_J * h3 * wt;
				}
			}
 
			advection = 0.;

			if ( pd->e[eqn] & T_ADVECTION )
			{
				if (pd->v[VELOCITY1]) /* then must be solving fluid mechanics in this material */
				{
					
					/*
					 * Standard incompressibility constraint means we have
					 * a solenoidal velocity field...if we don't, then
					 * we might be in serious trouble...
					 *
					 * HKM -> This is where a variable density implementation must go
					 */
					
					advection = div_v;
					
					/* We get a more complicated advection term because the
					* particle phase is not incompressible.
					*/
					if( particle_momentum_on ) advection *= ompvf;
					
					advection *= phi_i * h3 * det_J * wt;
					advection *= pd->etm[eqn][(LOG2_ADVECTION)];
				}
				else if  (cr->MeshMotion == LAGRANGIAN ||
						  cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
						  cr->MeshMotion == TOTAL_ALE)
					/* use divergence of displacement for linear elasticity */
				{
					advection = fv->volume_change;
					
					if( particle_momentum_on ) advection *= ompvf;
					
					advection *= phi_i * h3 * det_J * wt;
					advection *= pd->etm[eqn][(LOG2_ADVECTION)];
				}
				if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
					mp->SpeciesSourceModel[0]  == ION_REACTIONS)
					{
					advection = div_v;
					var = MASS_FRACTION;
					for (w=0; w<pd->Num_Species-1; w++)
					{
						for ( j=0; j<ei->dof[var]; j++)
						{
							if ( bf[var]->phi[j] > 0.0 ) break;
						}
						derivative = d_rho->C[w][j]/bf[var]->phi[j];
						sum = 0.;
						for (p=0; p<dim; p++)
						{
							sum += s_terms.conv_flux[w][p];
						}
						advection += derivative * sum / rho;
					}
					advection *= phi_i * h3 * det_J * wt;
					advection *= pd->etm[eqn][(LOG2_ADVECTION)];
					}
			}


			source = 0.;

			if ( pd->e[eqn] & T_SOURCE )
			{
				/*
				 * Maybe you want to penalize the equation to give a
				 * nonzero diagonal entry...
				 */
				
				if (pd->v[VELOCITY1])
				{
				  if (mp->DensityModel == CONSTANT ||
				      mp->DensityModel == DENSITY_FILL ||
				      mp->DensityModel == DENSITY_LEVEL_SET ||
				      mp->DensityModel == DENSITY_IDEAL_GAS )
				    {
				      /* These density models assume constant local density.  That is,
					 the density may be different in different parts of the flow
					 but in a sufficiently small region the density behavior is flat
				      */
				      /* DRN (07/13/05):
					 This was previously:
					 source     =  P;
					 But this messes with level set problems that have a density source
					 over part of the domain and constant in other regions.  If someone was
					 counting on this behavior as a form of a penalty method to give a non-zero
					 diagonal entry, we should implement a new density model that accomplishes this.
					 I really don't know what you want for DENSITY_IDEAL_GAS, though?!?
				      */
				      source = 0.;
				    }
				  else if (mp->DensityModel == DENSITY_FOAM || 
					   mp->DensityModel == DENSITY_FOAM_CONC || 
					   mp->DensityModel == DENSITY_FOAM_TIME ||
					   mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
				    {
				      /*  These density models locally permit a time and spatially varying
				       *  density.  Consequently, neither the gradient nor the time derivative 
				       *  of the density term in the continuity equation is zero. These terms
				       *  must then be included here as a "source" term
				       */
				      source = FoamVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx,
								dFVS_dC, dFVS_dF);
				    }
				  else if ( mp->DensityModel == REACTIVE_FOAM )
				    {
				      /* These density models locally permit a time and spatially varying
					 density.  Consequently, the Lagrangian derivative of the density
					 terms in the continuity equation are not zero and are
					 included here as a source term
				      */
				      source = REFVolumeSource( time_value,
								dt,
								tt,
								dFVS_dv,
								dFVS_dT,
								dFVS_dx,
								dFVS_dC );
				    }
				  else if
				    ( mp->DensityModel == SUSPENSION ||
				      mp->DensityModel == SUSPENSION_PM )
				    {
				      /* Although the suspension density models meet the definition of
					 a locally variable density model, the Lagrangian derivative
					 of their densities can be represented as a divergence of
					 mass flux.  This term must be integrated by parts and so is
					 included separately later on and is not include as part of the "source"
					 terms */
				      source = 0.0;
				    }
					
				  /* To include, or not to include, that is the question
				   * when considering the particle momentum coupled eq.'s...
				   */
					
				  /* We get this source term because the fluid phase is not
				   * incompressible.
				   */
				  if( particle_momentum_on )
				    {
				      source = 0.0;
				      for( a=0; a<wim; a++ )
					{
					  /* Cannot use s_terms.conv_flux[a] here because that
					   * is defined in terms of the particle phase
					   * velocities.
					   */
					  source -= fv->grad_c[species][a] * v[a];
					}
				    }
				  source     *= phi_i * h3 * det_J * wt;
				  source     *= pd->etm[eqn][(LOG2_SOURCE)];
				}
				
				if (( cr->MeshMotion == LAGRANGIAN ||
					  cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
					  (cr->MeshMotion == TOTAL_ALE && !pd->v[VELOCITY1]) ))
					/* add swelling as a source of volume */
				{
					if ( mp->PorousMediaType == CONTINUOUS )
					{
						source = - (1. - initial_volsolvent) / (1. - volsolvent);
						source *= phi_i * h3 * det_J * wt;
						source *= pd->etm[eqn][(LOG2_SOURCE)];
						
					}
				}
				if (mp->SpeciesSourceModel[0]  == ELECTRODE_KINETICS ||
					mp->SpeciesSourceModel[0]  == ION_REACTIONS)
				{
					source = 0.;
					for (j=0; j<pd->Num_Species; j++)
					{
						source -= s_terms.MassSource[j] * mp->molecular_weight[j];
					}
					source /= rho;
					source *= phi_i * h3 * det_J * wt;
					source *= pd->etm[eqn][(LOG2_SOURCE)];
				}
			}

			/* add Pressure-Stabilized Petrov-Galerkin term
				* if desired.
				*/
			
			pressure_stabilization = 0.;

			if(PSPG)
			{
				for ( a=0; a<wim; a++)
				{
					meqn = R_MOMENTUM1+a;
					if( pd->e[meqn])
					{
						pressure_stabilization += grad_phi[i][a] * pspg[a];
					}
				}
				pressure_stabilization *= h3 * det_J * wt;
			}

			
			h_flux = 0.0;

			
			if ( (cr->MassFluxModel == HYDRODYNAMIC) && ( mp->MomentumSourceModel == SUSPENSION ) )
			{
				/* add divergence of particle phase flux as source term */
				
				/* The particle flux terms has been integrated by parts.
				* No boundary integrals are included in this formulation
				* so it is tacitly assumed that the particle phase
				* relative mass flux over all boundaries is zero
				*/
				
				for( p=0; p<dim; p++)
				{
					h_flux += grad_phi[i][p]*s_terms.diff_flux[w0][p];
				}
				
				h_flux *= (rhos - rhof )/rhof;
				h_flux *= h3 * det_J * wt;
				h_flux *= pd->etm[eqn][(LOG2_ADVECTION)];
				/*  h_flux = 0.0; */
			}
			
			continuity = mass + advection + source + pressure_stabilization + h_flux;

			
			var = FILL;
			pvar = upd->vp[var];
			dofs = ei->dof[var];
			j=0;	
			while ( j < dofs ) 
			{
				
				lec->J[peqn][pvar][i][j] += lsi->d_H_dF[j] * continuity * sign;
			
				j++;
				
			}
		}
    }

  return(status);

} /* end of function assemble_continuity_path_dependence              */
/**********************************************************************/

int
assemble_LM_source ( double *xi,
                     int oAC,        /* Flag indicating calling function */
                     double *gAC,    /* Augmenting condition arrays */
                     double **bAC,
                     double **cAC,
                     double x[], Exo_DB *exo )
{
  int i, j, a, ii, ledof, v;
  int eqn, peqn;
  int dim, wim;
  int pass;
  int ac_lm = Do_Overlap;
  int id_side, nu, iAC=0, ioffset = 0;
  int dof_l[DIM];
  double phi_l[DIM][MDE];
  int ielem = ei->ielem;
  int iconn_ptr = exo->elem_ptr[ielem];
  double lagrange_mult[3] = {0.0, 0.0, 0.0};

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i;

  double source;

  load_lsi( ls->Length_Scale );
  if (lsi->delta == 0.) return(0);

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  pass = ( (Do_Overlap && oAC >= 0) ? 2 : 1);
  
  ielem = ei->ielem;
  if ( ( ls->Evolution == LS_EVOLVE_SLAVE ) && ( ls->init_surf_list->start->type == LS_SURF_SS ) )
    {
      struct LS_Surf *ss_surf;
      struct LS_Surf_Closest_Point *cp;

      ss_surf = closest_surf( ls->init_surf_list, x, exo, fv->x );
      cp = ss_surf->closest_point;

      /* use kitchen sink approach for now at solids location */
      if ( cp->elem == -1 ) EH(-1,"Invalid element at closest_point");

      /* Associate the correct AC number with this point */
      id_side = cp->elem_side;
      ioffset = first_overlap_ac(cp->elem, id_side);
      if (ioffset == -1) EH(-1,"Bad AC index");

      setup_shop_at_point( cp->elem, cp->xi, exo );

      for (a=0; a<ei->ielem_dim; a++)
        {
          if (ac_lm)
            {
              lagrange_mult[a] = augc[ioffset+a].lm_value;
            }
          else
            {
              lagrange_mult[a] = fv->lm[a];
            }

          /* Grab stuff for cross term evaluation */
          if (pass == 2)
            {
              v = LAGR_MULT1 + a;
              if (pd->v[v] && !ac_lm)
                {
                  dof_l[a] = ei->dof[v];
                  for (j = 0; j < dof_l[a]; j++)
                    {
                      phi_l[a][j] = bf[v]->phi[j];
                    }
                }
              else
                {
                  dof_l[a] = 1;
                  phi_l[a][0] = 1.0;
                }
            }
        }

      /* Yep, that's all I needed from closest_point.
       * Now go back to the original element
       */
      setup_shop_at_point( ielem, xi, exo );
    }
                      

  /*
   * Wesiduals ________________________________________________________________________________
   */
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

                  /* Fluid momentum residuals */
                  /*lec->R[peqn][i] += wt *(lagrange_mult[a]) *
                    det_J * 0.5 * length; */
                  source = -phi_i * (lagrange_mult[a]) * lsi->delta;
                  
                  source *= det_J * wt;
                  source *= h3;
                  source *= pd->etm[eqn][(LOG2_SOURCE)];
                                 
                  if (pass == 1) lec->R[peqn][ii] += source;

                  /* Sensitivities to solid Lagrange multiplier
                     (cross-mesh term) */
                  if (Do_Overlap && pass == 2)
                    {
                      v = LAGR_MULT1 + a;
                      iAC = ioffset + a;
                      nu = lookup_active_dof(eqn, i, iconn_ptr);

                      for (j = 0; j < dof_l[a]; j++)
                        {
                          if (nu >= 0)
                            {
                              source = -phi_i * phi_l[a][j] * lsi->delta;

                              source *= det_J * wt;
                              source *= h3;
                              source *= pd->etm[eqn][(LOG2_SOURCE)];
                              
                              bAC[iAC][nu] += source;
                            }
                        }
                    }
		}
	    }
	}

  return ( 1 );
}

/*
 * Calculate the total stress tensor for a fluid at a single gauss point
 *  This includes the diagonal pressure contribution
 *
 *  Pi = stress tensor
 *  d_Pi = dependence of the stress tensor on the independent variables
 */
void
fluid_stress( double Pi[DIM][DIM],
              STRESS_DEPENDENCE_STRUCT *d_Pi)
{

  /*
   * Variables for vicosity and derivative
   */
  dbl mu = 0.0;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /* polymer viscosity and derivatives */

  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;
  dbl mu_over_mu_num = 0.0;

  /* solvent viscosity and derivatives */

  dbl mus = 0.0;
  VISCOSITY_DEPENDENCE_STRUCT d_mus_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mus = &d_mus_struct;

  /* numerical "adaptive" viscosity and derivatives */

  dbl mu_num;
  dbl d_mun_dS[MAX_MODES][DIM][DIM][MDE];
  dbl d_mun_dG[DIM][DIM][MDE];

  /* Dilational viscosity */
  dbl kappa = 0.0;
  DILVISCOSITY_DEPENDENCE_STRUCT d_dilMu_struct;
  DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu = &d_dilMu_struct;
  int kappaWipesMu = 1;

  /* particle stress for suspension balance model*/
  dbl tau_p[DIM][DIM];
  dbl d_tau_p_dv[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dvd[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dy[DIM][DIM][MAX_CONC][MDE];
  dbl d_tau_p_dmesh[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dp[DIM][DIM][MDE];
  int w0; /* suspension species number */

  dbl *grad_v[DIM];		        /* Gradient of v. */
  dbl gamma[DIM][DIM];                  /* shrearrate tensor based on velocity */
  dbl s[DIM][DIM];                      /* polymer stress tensor */
  dbl gamma_cont[DIM][DIM];             /* shearrate tensor based on continuous gradient of velocity */
  dbl P;

  dbl evss_f;      /* flag to add in the continuous velocity
		      gradient tensor for Fortin's formulation */
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];
  int mode;       /* index for modal viscoelastic counter */


  //! Flag for doing dilational viscosity contributions.
  int do_dilational_visc = 0;
  int dim, wim;

  int a, b, p, q, j, w, c, var;

  dbl (* grad_phi_e ) [DIM][DIM][DIM] = NULL;

  int eqn = R_MOMENTUM1;

 
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  if( pd->v[POLYMER_STRESS11] )
    {
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
    }

  /*
   * Field variables...
   */

  P = fv->P;

  
  /* if d_Pi == NULL, then the dependencies aren't needed,
     so we won't need the viscosity dependencies, either
  */
  if ( d_Pi == NULL )
    {
      d_mu = NULL;
      d_mus = NULL;
      d_mup = NULL;
    }
      

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for ( a=0; a<VIM; a++ ) {
    grad_v[a] = fv->grad_v[a];
  }
  if ( pd->v[POLYMER_STRESS11] )
    {
      memset( s, 0, sizeof(dbl)*DIM*DIM);
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      for ( mode=0; mode<vn->modes; mode++)
		{
		  s[a][b] += fv->S[mode][a][b];
		}
	    }
	}
    }

  memset( tau_p,      0, sizeof(double) * DIM*DIM);
  memset( d_tau_p_dv, 0, sizeof(double) * DIM*DIM*DIM*MDE);
  memset( d_tau_p_dvd, 0, sizeof(double) * DIM*DIM*DIM*MDE);
  memset( d_tau_p_dy, 0, sizeof(double) * DIM*DIM*MAX_CONC*MDE);
  memset( d_tau_p_dmesh, 0, sizeof(double) * DIM*DIM*DIM*MDE);
  memset( d_tau_p_dp, 0, sizeof(double) * DIM*DIM*MDE);

/*   if( cr->MassFluxModel == DM_SUSPENSION_BALANCE|| mp->QTensorDiffType[w0] == CONSTANT) */
  if( cr->MassFluxModel == DM_SUSPENSION_BALANCE || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR)
    {
      w0 = gn->sus_species_no;
      particle_stress(tau_p, d_tau_p_dv, d_tau_p_dvd,d_tau_p_dy,d_tau_p_dmesh,d_tau_p_dp, w0);
    }

  if ( pd->v[POLYMER_STRESS11] && (vn->evssModel == EVSS_F) )
    {
      evss_f = 1.;
    }
  else
    {
      evss_f = 0.;
    }

  if ( evss_f )
    {
      for ( a=0; a<VIM; a++)
	{
	  for ( b=0; b<VIM; b++)
	    {
	      gamma_cont[a][b] = fv->G[a][b] + fv->G[b][a] ;
	    }
	}
    }
  else
    {
      memset(gamma_cont, 0, sizeof(double)*DIM*DIM);
    }

  /* load up shear rate tensor based on velocity */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }
  if (do_dilational_visc) {
    // use previously calculated div_v.
  }

  if ( pd->v[POLYMER_STRESS11] )
    {
      mus = viscosity(gn, gamma, d_mus);
    }
  else
    {
      mu = viscosity(gn, gamma, d_mu);
    }

  if ( pd->v[POLYMER_STRESS11] )
    {
      /* initialize the derivative wrt to stress and velocity gradient */

      memset( d_mun_dS, 0, sizeof(double) * MAX_MODES*DIM*DIM*MDE);
      memset( d_mun_dG, 0, sizeof(double) * DIM*DIM*MDE);

      /* This is the adaptive viscosity from Sun et al., 1999.
       * The term multiplies the continuous and discontinuous
       * shear-rate, so it should cancel out and not affect the
       * solution, other than increasing the stability of the
       * algorithm in areas of high shear and stress.
       */

      mu_num = numerical_viscosity(s, gamma_cont, d_mun_dS, d_mun_dG);

      mu_over_mu_num = mus;
      mu = mu_num * mus;

      /* first add the solvent viscosity to the total viscosity
       * including all the derivatives.
       */
      var = VELOCITY1;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( a=0; a<wim; a++)
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_mu->v[a][j] = mu_num * d_mus->v[a][j];
		}
	    }
	}

      var = MESH_DISPLACEMENT1;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( a=0; a<dim; a++)
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_mu->X[a][j] = mu_num * d_mus->X[a][j];
		}
	    }
	}

      var = TEMPERATURE;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_mu->T[j] = mu_num * d_mus->T[j];
	    }
	}

      var = BOND_EVOLUTION;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_mu->nn[j] = mu_num * d_mus->nn[j];
	    }
	}

#ifdef COUPLED_FILL
      var = FILL;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_mu->F[j] = mu_num * d_mus->F[j];
	    }
	}
#endif /* COUPLED_FILL */

      var = PRESSURE;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      d_mu->P[j] = mu_num * d_mus->P[j];
	    }
	}

      var = MASS_FRACTION;
      if ( d_Pi != NULL && pd->v[var] )
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_mu->C[w][j] = mu_num * d_mus->C[w][j];
		}
	    }
	}

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

      for ( mode=0; mode<vn->modes; mode++)
        {
          /* get polymer viscosity */
          mup = viscosity(ve[mode]->gn, gamma, d_mup);

          mu_over_mu_num += at * mup;
          mu += mu_num * at * mup;

          var = VELOCITY1;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( a=0; a<wim; a++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_mu->v[a][j] += mu_num * at * d_mup->v[a][j];
                    }
                }
            }

          var = MESH_DISPLACEMENT1;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( a=0; a<dim; a++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_mu->X[a][j] += mu_num * at * d_mup->X[a][j];
                    }
                }
            }

          var = TEMPERATURE;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_mu->T[j] += mu_num * (at * d_mup->T[j] + mup * d_at_dT[j]);
                }
            }

          var = BOND_EVOLUTION;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_mu->nn[j] += mu_num * at * d_mup->nn[j];
                }
            }

#ifdef COUPLED_FILL
          var = FILL;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_mu->F[j] += mu_num * at * d_mup->F[j];
                }
            }
#endif /* COUPLED_FILL */

          var = PRESSURE;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_mu->P[j] += mu_num * at * d_mup->P[j];
                }
            }

          var = MASS_FRACTION;
          if ( d_Pi != NULL && pd->v[var] )
            {
              for ( w=0; w<pd->Num_Species_Eqn; w++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_mu->C[w][j] += mu_num * at * d_mup->C[w][j];
                    }
                }
            }
        }
    }

  /*
   * Calculate the dilational viscosity, if necessary
   */
  if (mp->DilationalViscosityModel != DILVISCM_KAPPAWIPESMU) {
    kappa = dil_viscosity(gn, gamma, mu, d_mu, d_dilMu);
    kappa = 0.0;
    kappaWipesMu = 0;
  }

  /* xxxxx HKM xxxx  -> operational point */
  /*
   *  tau_p[a][b] is the particle stress contribution, so
   *              it's zero for non-particle cases
   *  gamma[a][b] is the usual gradV + gradV_t term
   *  delta(a,b) is a delta function
   */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
        {
          Pi[a][b]  = -P * (double)delta(a,b)
	    + mu * gamma[a][b] - tau_p[a][b];
        }

      // Add in the diagonal contribution
      if (!kappaWipesMu) {
	Pi[a][a] -= (mu / 3.0 - 0.5 * kappa) * gamma[a][a];
      }
    }

  if ( pd->v[POLYMER_STRESS11] )
    {
      for ( a=0; a<VIM; a++)
        {
          for ( b=0; b<VIM; b++)
            {
              Pi[a][b]  += - evss_f * (mu - mus) * gamma_cont[a][b]
		+ s[a][b];
            }
        }
    }

  /*
   * OK, FIND THE JACOBIAN
   */
  var = TEMPERATURE;
  if ( d_Pi != NULL && pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_Pi->T[p][q][j] = d_mu->T[j] * gamma[p][q];
                }
            }
        }
      if (!kappaWipesMu) {
	for (p = 0; p < VIM; p++) {
	  for (j = 0; j < ei->dof[var]; j++) {
	    d_Pi->T[p][p][j] -= (d_mu->T[j]/3.0 - 0.5 * d_dilMu->T[j]) * gamma[p][p];
	  }
	}
      }
      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_Pi->T[p][q][j] -= evss_f * ( d_mu->T[j] - d_mus->T[j] )* gamma_cont[p][q];
                    }
                }
	    }
        }
    }

  var = BOND_EVOLUTION;
  if ( d_Pi != NULL && pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_Pi->nn[p][q][j] = d_mu->nn[j] * gamma[p][q];
                }
            }
        }
      if (!kappaWipesMu) {
	for (p = 0; p < VIM; p++) {
	  for (j = 0; j < ei->dof[var]; j++) {
	    d_Pi->nn[p][p][j] -= (d_mu->nn[j]/3.0 - 0.5 * d_dilMu->nn[j]) * gamma[p][p];
	  }
	}
      }
      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_Pi->nn[p][q][j] -= evss_f * ( d_mu->nn[j] - d_mus->nn[j] )* gamma_cont[p][q];
                    }
                }
	    }
        }
    }

#ifdef COUPLED_FILL
  var = FILL;
  if ( d_Pi != NULL && pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_Pi->F[p][q][j] = d_mu->F[j] * gamma[p][q];
                }
            }
        }
      if (!kappaWipesMu) {
	for (p = 0; p < VIM; p++) {
	  for (j = 0; j < ei->dof[var]; j++) {
	    d_Pi->F[p][p][j] -= (d_mu->F[j] / 3.0 - 0.5 * d_dilMu->F[j]) * gamma[p][p];
	  }
	}
      }
      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_Pi->F[p][q][j] -= evss_f * ( d_mu->F[j] - d_mus->F[j] ) * gamma_cont[p][q];
                    }
                }
	    }
        }
    }
#endif /* COUPLED_FILL */

  
  if ( d_Pi != NULL && pd->v[PHASE1] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
	      for( a=0; a<pfd->num_phase_funcs; a++)
		{
		  var = PHASE1 + a;
		  for( j=0 ; j<ei->dof[var] ; j++)
		    {
		      d_Pi->pf[p][q][a][j] = d_mu->pf[a][j]*gamma[p][q];
		    }
		}
	    }
	}
      if (!kappaWipesMu) {
	for (p = 0; p < VIM; p++) {
	  for (a = 0; a < pfd->num_phase_funcs; a++) {
	    var = PHASE1 + a;
	    for (j = 0; j < ei->dof[var]; j++) {
	      d_Pi->pf[p][p][a][j] -= (d_mu->pf[a][j]/3.0 - 0.5 * d_dilMu->pf[a][j]) * gamma[p][p];
	    }
	  }
	}
      }
      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
		  for( a=0; a<pfd->num_phase_funcs; a++)
		    {
		      var = PHASE1 + a;
		      for ( j=0; j<ei->dof[var]; j++)
			{
			  d_Pi->pf[p][q][a][j] -= evss_f * ( d_mu->pf[a][j] - d_mus->pf[a][j] ) * gamma_cont[p][q];
			}
		    }
                }
	    }
        }

    }



  if ( d_Pi != NULL && pd->v[VELOCITY1] )
    {
      /* Damn... It is unfortunate that grad_phi_e is 1) assumed to be
       * the same for all velocity components (this is bad for 3d
       * stability of a 2d flow), and 2) there is no bf[] structure
       * for the theta-velocity in cylindrical coordinates (I believe
       * there is for swirling).  That's why this needed to be
       * split... */
      if(pd->CoordinateSystem != CYLINDRICAL) {
        for ( p=0; p<VIM; p++)
          {
            for ( q=0; q<VIM; q++)
              {
                for ( b=0; b<wim; b++)
                  {
                    for ( j=0; j<ei->dof[VELOCITY1]; j++)
                      {
                        /* grad_phi_e cannot be the same for all
                         * velocities for 3d stab of 2d flow!!
                         * Compare with the old way in the CYLINDRICAL
                         * chunk below... */
                        d_Pi->v[p][q][b][j] =
			  mu * (bf[VELOCITY1+q]->grad_phi_e[j][b][p][q] +
				bf[VELOCITY1+p]->grad_phi_e[j][b][q][p] )
                          + d_mu->v[b][j] * gamma[p][q]
                          - d_tau_p_dv[p][q][b][j];
                      }
                  }
              }
	    if (!kappaWipesMu) {
	      for (b = 0; b < wim; b++) {
		for (j = 0; j < ei->dof[VELOCITY1]; j++) {
		  d_Pi->v[p][p][b][j] -=
		    ((2.0 * mu / 3.0 - kappa) * (bf[VELOCITY1+p]->grad_phi_e[j][b][p][p]) +
		     (d_mu->v[b][j] / 3.0 - 0.5 * d_dilMu->v[b][j]) * gamma[p][p]);
		}
	      }
	    }
	  }
      } else {
	/* For CYLINDRICAL, we can assume that all of the velocity
	 * components share the same interpolating basis function.
	 * This was not allowable for 3d stability of 2d flow and
	 * the PROJECTED_CARTESIAN coordinate system.  In fact, the
	 * behavior above is "more" correct than this, so I
	 * defaulted everyrthing else to it.
	 */
	grad_phi_e = bf[eqn]->grad_phi_e;
	for ( p=0; p<VIM; p++) {
	  for ( q=0; q<VIM; q++)
	    {
	      for ( b=0; b<wim; b++)
		{
		  for ( j=0; j<ei->dof[VELOCITY1]; j++)
		    {
		      d_Pi->v[p][q][b][j] =
			mu * (grad_phi_e[j][b][p][q] + grad_phi_e[j][b][q][p])
			+ d_mu->v[b][j] * gamma[p][q]
			- d_tau_p_dv[p][q][b][j];
		    }
		}
	    }
	  if (!kappaWipesMu) {
	    for (b = 0; b < wim; b++) {
	      for (j = 0; j < ei->dof[VELOCITY1]; j++) {
		d_Pi->v[p][p][b][j] -=
		  ((2.0 * mu / 3.0 - kappa) * (grad_phi_e[j][b][p][p]) +
		   (d_mu->v[b][j] / 3.0 - 0.5 * d_dilMu->v[b][j]) * gamma[p][p]);
	      }
	    }
	  }
	}
      }

      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( b=0; b<wim; b++)
                    {
                      for ( j=0; j<ei->dof[var]; j++)
                        {
                          d_Pi->v[p][q][b][j] -= evss_f * ( d_mu->v[b][j] - d_mus->v[b][j] ) * gamma_cont[p][q];
                        }
                    }
                }
	    }
        }
    }

  // Vorticity direction dependence for qtensor
  if ( d_Pi != NULL && pd->v[VORT_DIR1] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( b=0; b<dim; b++)
                {
                  for ( j=0; j<ei->dof[VORT_DIR1]; j++)
                    {
                      d_Pi->vd[p][q][b][j] = - d_tau_p_dvd[p][q][b][j];
                    }
                }
            }
	}
    }

  // Mesh Dependence
  if ( d_Pi != NULL && pd->v[MESH_DISPLACEMENT1] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( b=0; b<dim; b++)
                {
                  for ( j=0; j<ei->dof[MESH_DISPLACEMENT1]; j++)
                    {
                      d_Pi->X[p][q][b][j] =
			mu * ( fv->d_grad_v_dmesh[p][q] [b][j]
                               + fv->d_grad_v_dmesh[q][p] [b][j] ) +
                        + d_mu->X [b][j] * gamma[p][q] - d_tau_p_dmesh[p][q][b][j];
                    }
                }
            }
	  if (!kappaWipesMu) {   
	    for (b = 0; b < dim; b++) {
	      for (j = 0; j<ei->dof[MESH_DISPLACEMENT1]; j++) {
		d_Pi->X[p][p][b][j] -=
		  ((2.0 * mu /3.0 - kappa) * (fv->d_grad_v_dmesh[p][p][b][j]) +
		   (d_mu->X[b][j] / 3.0 - 0.5 * d_dilMu->X[b][j]) * gamma[p][p]);
	      }
	    }    
	  }
	}
      
      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( b=0; b<dim; b++)
                    {
                      for ( j=0; j<ei->dof[MESH_DISPLACEMENT1]; j++)
                        {
                          d_Pi->X[p][q][b][j] -=
			    evss_f * (d_mu->X [b][j] - d_mus->X [b][j]) * gamma_cont[p][q];
                        }
                    }
                }
            }
        }
    }

  if ( d_Pi != NULL && pd->v[POLYMER_STRESS11] )
    {
      for ( mode=0; mode<vn->modes; mode++)
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( b=0; b<VIM; b++)
                    {
                      for ( c=0; c<VIM; c++)
                        {
                          var = v_s[mode][b][c];
                          for ( j=0; j<ei->dof[var]; j++)
                            {
                              d_Pi->S[p][q][mode][b][c][j] =
				((double)delta(b,p) * (double)delta(c,q)) * bf[var]->phi[j]
                                + mu_over_mu_num * d_mun_dS[mode][b][c][j] *
				( gamma[p][q] - evss_f * gamma_cont[p][q] );
                            }
                        }
                    }
                }
            }
        }
    }

  if ( d_Pi != NULL && pd->v[VELOCITY_GRADIENT11] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( b=0; b<VIM; b++)
                {
                  for ( c=0; c<VIM; c++)
                    {
                      var = v_g[b][c];
                      for ( j=0; j<ei->dof[var]; j++)
                        {
                          d_Pi->g[p][q][b][c][j] =
			    mu_over_mu_num * d_mun_dG[b][c][j] *
			    ( gamma[p][q] - evss_f * gamma_cont[p][q] )
                            + evss_f *
			    ( (mu - mus) * bf[var]->phi[j] *
			      ((double)delta(c,p) * (double)delta(b,q) +
			       (double)delta(b,p) * (double)delta(c,q)) );
                        }
                    }
                }
            }
        }
    }

  var = MASS_FRACTION;
  if ( d_mu != NULL && pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for (w = 0; w < pd->Num_Species_Eqn; w++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_Pi->C[p][q][w][j] = d_mu->C[w][j] * gamma[p][q]
			- d_tau_p_dy[p][q][w][j];
                    }
                }
            }
	  if (!kappaWipesMu) {   
	    for (w = 0; w < pd->Num_Species_Eqn; w++) {
	      for (j = 0; j < ei->dof[var]; j++) {
		d_Pi->C[p][p][w][j] -= (d_mu->C[w][j] /3.0 - 0.5 * d_dilMu->C[w][j]) * gamma[p][p];
	      }
	    }
	  }    
	}
        

      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( w=0; w<pd->Num_Species_Eqn; w++)
                    {
                      for ( j=0; j<ei->dof[var]; j++)
                        {
                          d_Pi->C[p][q][w][j] -=
                            evss_f * ( d_mu->C[w][j] - d_mus->C[w][j] ) * gamma_cont[p][q];
                        }
                    }
                }
            }
        }
    }

  var = PRESSURE;
  if ( d_Pi != NULL && pd->v[var] )
    {
      for ( p=0; p<VIM; p++)
        {
          for ( q=0; q<VIM; q++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_Pi->P[p][q][j] =  - (double)delta(p,q) * bf[var]->phi[j]
		    + d_mu->P[j] * gamma[p][q] - d_tau_p_dp[p][q][j];
                }
            }
	  if (!kappaWipesMu) {
	    for (j = 0; j < ei->dof[var]; j++) {
	      d_Pi->P[p][p][j] -= (d_mu->P[j] /3.0 - 0.5 * d_dilMu->P[j]) * gamma[p][p];
	    }
	  }  
        }

      if ( pd->v[POLYMER_STRESS11] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( q=0; q<VIM; q++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_Pi->P[p][q][j] -= evss_f * ( d_mu->P[j] - d_mus->P[j] )* gamma_cont[p][q];
                    }
                }
	    }
        }
    }
}

void
heat_flux( double q[DIM],
           HEAT_FLUX_DEPENDENCE_STRUCT *d_q,
           double time )
{
  dbl k;				/* Thermal conductivity. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl dq_gradT[DIM][DIM];               /*  Heat flux sensitivities  */
  dbl dq_dX[DIM][DIM];                  /*  Heat flux sensitivities wrt position  */
  dbl grad_T[DIM];			/* Temperature gradient. */

  int b, j = -1, p, a, w;
  int var;

  int dim = pd->Num_Dim;
  int wim;

  wim = dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN) wim = wim+1;

  if ( d_q == NULL ) d_k = NULL;
  
  k = conductivity( d_k, time );


  for ( p=0; p<VIM; p++)
    {
      grad_T[p] = fv->grad_T[p];
    }

  if ( cr->HeatFluxModel == CR_HF_FOURIER_0 )
    {
      for ( p=0; p<VIM; p++)
	{
	  q[p] = - k * grad_T[p];
	}

      var = TEMPERATURE;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->T[p][j] = - k * bf[var]->grad_phi[j][p]
		    - d_k->T[j] * grad_T[p];
                }
            }
        }

      var = MASS_FRACTION;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( w=0; w<pd->Num_Species_Eqn; w++)
            {
              for ( p=0; p<VIM; p++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_q->C[p][w][j] =  - d_k->C[w][j] * grad_T[p];
                    }
                }
            }
        }

      var = FILL;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->F[p][j] =  - d_k->F[j] * grad_T[p];
                }
            }
        }

      if (d_q != NULL && pd->v[MESH_DISPLACEMENT1])
        {
          for (p=0; p<VIM; p++)
	    {
	      for (b=0; b<wim; b++)
                {
                  var = MESH_DISPLACEMENT1+b;
                  for (j = 0; j < ei->dof[var]; j++)
                    {
                      d_q->X[p][b][j] =
                        -k * fv->d_grad_T_dmesh[p] [b][j]
                        -d_k->X[b][j]* grad_T[p];
                    }
                }
            }
        }
    }
  else if (cr->HeatFluxModel == CR_HF_USER)
    {
#ifdef SECOR_HEAT_FLUX
      double *hpar, h, dh_dX[DIM], Vb[DIM],Vt[DIM];
      double dq_dVb[DIM][DIM], dq_dVt[DIM][DIM];
#endif
#if 0
        int npadex,npadey;
        double dnum,dden,*pn,*pd;
        double hx,hy,dnumdx,ddendx;
#endif

/*     gap functions - product of pade functions */
#ifdef SECOR_HEAT_FLUX
        hpar = &mp->u_thermal_conductivity[0];
        h = hpar[0] + hpar[4]*fv->x[0]
                + (hpar[1]-hpar[5]*fv->x[0])*(hpar[3]-fv->x[1])
                        + 0.5*hpar[2]*SQUARE(hpar[3]-fv->x[1]);

        dh_dX[0] = hpar[4] - hpar[5]*(hpar[3]-fv->x[1]);
        dh_dX[1] = hpar[5]*fv->x[0]-hpar[1] - hpar[2]*(hpar[3]-fv->x[1]);
#endif
#if 0
        npadex = ((int)mp->u_thermal_conductivity[0]);
        pn = &mp->u_thermal_conductivity[1];
        pd = &mp->u_thermal_conductivity[npadex+1];
        dnum = pn[npadex-1];
        dden = pd[npadex-2];
        dnumdx = (npadex-1.)*pn[npadex-1];
        ddendx = (npadex-2.)*pd[npadex-2];
        for (i=npadex-2;i>=0;i--)
                { dnum = pn[i]+fv->x[0]*dnum; }
        for (i=npadex-3;i>=0;i--)
                { dden = pd[i]+fv->x[0]*dden;
                  dnumdx = (i+1.)*pn[i+1]+fv->x[0]*dnumdx;}
        for (i=npadex-4;i>=0;i--)
                { ddendx = (i+1.)*pd[i+1]+fv->x[0]*ddendx; }
        hx = dnum/dden;
        dhdx = (dden*dnumdx - dnum*ddendx)/(dden*dden);
        npadey = ((int)mp->u_thermal_conductivity[2*npadex]);
        pn = &mp->u_thermal_conductivity[2*npadex+1];
        pd = &mp->u_thermal_conductivity[2*npadex+1+npadey];
        dnum = pn[npadey-1];
        dden = pd[npadey-2];
        dnumdx = (npadey-1.)*pn[npadey-1];
        ddendx = (npadey-2.)*pd[npadey-2];
        for (i=npadey-2;i>=0;i--)
                { dnum = pn[i]+fv->x[1]*dnum; }
        for (i=npadey-3;i>=0;i--)
                { dden = pd[i]+fv->x[1]*dden;
                  dnumdx = (i+1.)*pn[i+1]+fv->x[1]*dnumdx;}
        for (i=npadey-4;i>=0;i--)
                { ddendx = (i+1.)*pd[i+1]+fv->x[1]*ddendx; }
        hy = dnum/dden;
        dhdy = hx*(dden*dnumdx - dnum*ddendx)/(dden*dden);
        dhdx *= hy;
        h = hx*hy;
#endif

#ifdef SECOR_HEAT_FLUX
/*     velocities of bottom and top surfaces   */
        Vb[0] = mp->u_heat_capacity[0];
        Vb[1] = mp->u_heat_capacity[1];
        Vt[0] = mp->u_heat_capacity[2];
        Vt[1] = mp->u_heat_capacity[3];

      usr_heat_flux(grad_T, q, dq_gradT, dq_dX, time, h, dh_dX, Vb, Vt,
                             dq_dVb, dq_dVt);
#else
      usr_heat_flux(grad_T, q, dq_gradT, dq_dX, time);
#endif

      var = TEMPERATURE;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->T[p][j] = 0.0;
                  for ( a=0; a<VIM; a++)
                    {
                      d_q->T[p][j] += dq_gradT[p][a] * bf[var]->grad_phi[j][a];
                    }
                }
            }
        }
      for ( b=0; b<wim; b++)
        {
          var = MESH_DISPLACEMENT1+b;
          if ( d_q != NULL && pd->v[var] )
            {
              for ( p=0; p<VIM; p++)
                {
                  d_q->X[p][b][j] = dq_dX[p][b] * bf[var]->phi[j];
                  for ( a=0 ; a<VIM ; a++ )
                    {
                      d_q->X[p][b][j] += dq_gradT[p][a] * fv->d_grad_T_dmesh[a][b][j];
                    }
                }
            }
        }
    }
  else
    {
      EH( -1, "Unimplemented thermal constitutive relation.");
    }
}

double heat_source( HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                    double time,/* current time */
                    double tt,	/* parameter to vary time integration from
				   explicit (tt = 1) to implicit (tt = 0) */
		    double dt )	/* current time step size */
{
  double h = 0.;
  double h_acous = 0.;
  int j, w, a, var_offset;
  int var, dim, wim;

  const double F = 96487.0;          /* Faraday's constant in units of C/euiv.; KSC: 2/17/99 */
  const double R = 8.314;            /* Universal gas constant in units of J/mole K */
  dbl TT;                            /* dummy variable */
  dbl ai0, ai0_anode, ai0_cathode;   /* exchange current density; KSC: 2/17/99 */
  int mn;

  struct Level_Set_Data *ls_old = ls;

  if (MAX_CONC < 3) {
    EH(-1, "heat_source expects MAX_CONC >= 3");
    return 0;
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim   = pd->Num_Dim;
  wim = dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN) wim = wim+1;

  /* initialize Heat Source sensitivities */
  if ( d_h != NULL )
    {
      memset( d_h->T, 0, sizeof(double)*MDE);
      memset( d_h->F, 0, sizeof(double)*MDE);
      memset( d_h->V, 0, sizeof(double)*MDE);
      memset( d_h->v, 0, sizeof(double)*DIM*MDE);
      memset( d_h->X, 0, sizeof(double)*DIM*MDE);
      memset( d_h->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_h->S, 0, sizeof(double)*MAX_MODES*DIM*DIM*MDE);
      memset( d_h->APR, 0, sizeof(double)*MDE);
      memset( d_h->API, 0, sizeof(double)*MDE);
      memset( d_h->INT, 0, sizeof(double)*MDE);
    }

  if(mp->HeatSourceModel == USER )
    {
      usr_heat_source(mp->u_heat_source,time);
      h = mp->heat_source;

      var = TEMPERATURE;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->T[j]= mp->d_heat_source[var]*bf[var]->phi[j];
            }
        }

      var = VOLTAGE;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->V[j]= mp->d_heat_source[var]*bf[var]->phi[j];
            }
        }

      if( d_h != NULL && pd->v[VELOCITY1])
	{
	  for ( a=0; a<wim; a++)
	    {
	      var = VELOCITY1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_h->v[a][j] =mp->d_heat_source[var]*bf[var]->phi[j];
		}
	    }
	}

      if( d_h != NULL && pd->v[MESH_DISPLACEMENT1] )
	{
	  for ( a=0; a<wim; a++)
	    {
	      var = MESH_DISPLACEMENT1 + a;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_h->X[a][j] =mp->d_heat_source[var]*bf[var]->phi[j];
		}
	    }
	}

      if ( d_h != NULL && pd->v[MASS_FRACTION] )
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_h->C[w][j] =mp->d_heat_source[var_offset]*bf[var]->phi[j];
		}
	    }
	}

    }
  else if(mp->HeatSourceModel == PHOTO_CURING )
    {
       double intensity;
       double k_prop, k_inh = 0, free_rad, d_free_rad_dI;
       double *param,dhdC[MAX_CONC],dhdT,dhdI,Conc[MAX_CONC];
       int model_bit, num_mon, O2_spec=-1, rad_spec=-1, init_spec = 0;
       double k_propX=1, k_propT=0, k_propX_num=0, k_propX_den=0;
       double intensity_cgs = 2.998e+10*8.85e-12/200.0;
       double dbl_small = 1.0e-15,Xconv_denom = 0.;
       double Xconv=0.0, Xconv_init=0.0, dXdC[MAX_CONC], sum_mon=0;

       param = mp->u_heat_source;
       model_bit = ((int)param[0]);
       h = 0;
       dhdT = 0;  dhdI = 0; d_free_rad_dI = 0;
       for(a=0; a<MAX_CONC; a++) dhdC[a]=0.;

       intensity = 0.;
       if(pd->e[R_LIGHT_INTP])
         {
         intensity += fv->poynt[0];
         if(pd->e[R_LIGHT_INTM])
          { intensity += fv->poynt[1];}
         if(pd->e[R_LIGHT_INTD])
          { intensity += fv->poynt[2];}
         intensity *= mp->u_species_source[init_spec][1];
         intensity = MAX(intensity,0.0);
         }
       else if(pd->e[R_ACOUS_PREAL])
         {
         intensity = mp->u_species_source[init_spec][1]*
                     intensity_cgs*
                     (SQUARE(fv->apr)+SQUARE(fv->api));
         }
       else
        { WH(-1,"No Intensity field found in PHOTO_CURING\n"); }

      /* insure concentrations are positive  */
	for ( j=0; j<pd->Num_Species_Eqn; j++)
	    {
	     Conc[j] = MAX(dbl_small,fv->c[j]);
	    }

/**  heat source from momomer heat of reaction     **/
        num_mon = model_bit>>2;

       if( (model_bit & 1)  && (model_bit & 2))
            {
            O2_spec = init_spec + num_mon +2;
            rad_spec = O2_spec + 1;
            free_rad = Conc[rad_spec];
            }
       else if( model_bit & 1)
            {
            O2_spec = init_spec + num_mon +2;
            k_inh = mp->u_species_source[O2_spec][1]*
                    exp(-mp->u_species_source[O2_spec][2]*
                    (1./fv->T - 1./mp->u_species_source[O2_spec][3]));
            free_rad = sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec])
                - k_inh*Conc[O2_spec]/2.;
            d_free_rad_dI += 0.5*Conc[init_spec]*
                              mp->u_species_source[init_spec+1][2]/
                           sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec]);
            }
       else if( model_bit & 2)
            {
            rad_spec = init_spec + num_mon +2;
            free_rad = Conc[rad_spec];
            }
       else
            {
            free_rad = sqrt(mp->u_species_source[init_spec+1][2]*
                       intensity*Conc[init_spec]);
            if (free_rad  > 0 )	{
            d_free_rad_dI += 0.5/free_rad
                     *mp->u_species_source[init_spec+1][2]*Conc[init_spec];
                    }
            }

       switch(mp->Species_Var_Type)   {
           case SPECIES_DENSITY:
                for ( w=init_spec+2; w<init_spec+2+num_mon; w++)
	            { 
                     Xconv += fv->external_field[w]/mp->molecular_weight[w]/mp->specific_volume[w];
                     sum_mon += Conc[w]/mp->molecular_weight[w];
                     dXdC[w] = -1.0/mp->molecular_weight[w];
                    }
                break;
           case SPECIES_CONCENTRATION:
                for ( w=init_spec+2; w<init_spec+2+num_mon; w++)
	            { 
                     Xconv += fv->external_field[w]/mp->specific_volume[w];
                     sum_mon += Conc[w];
                     dXdC[w] = -1.0;
                    }
                break;
           default:
                EH(-1,"invalid Species Type for PHOTO_CURING\n");
           }
       Xconv *= mp->specific_volume[pd->Num_Species_Eqn];
       Xconv_denom = Xconv + sum_mon;
       Xconv /= Xconv_denom;
       Xconv = MAX(dbl_small,Xconv);
       Xconv = MIN(1.0-dbl_small,Xconv);
       for ( w=init_spec+2; w<init_spec+2+num_mon; w++)
            { dXdC[w] *= Xconv/Xconv_denom; }
        if(Xconv <= dbl_small || Xconv >= (1.0-dbl_small) )
            { memset( dXdC, 0, sizeof(double) * MAX_CONC); }

      for(w=init_spec+2 ; w<init_spec+num_mon+2 ; w++)
         {
          k_prop = mp->u_species_source[w][1]*
                exp(-mp->u_species_source[w][2]*
                (1./fv->T - 1./mp->u_species_source[w][3]));
          k_propX_num = (1.0-mp->u_species_source[w][4])*(1.-Xconv)+mp->u_species_source[w][4]*(1.0-Xconv_init);
          k_propX_den = k_propX_num - (1.0-mp->u_species_source[w][4])*(1.-Xconv)*log((1.-Xconv)/(1.0-Xconv_init));
          k_propX = SQUARE(k_propX_num)/k_propX_den;
          k_propT = k_prop*k_propX;

          h += k_propT*Conc[w]*free_rad*param[w]*mp->molecular_weight[w];
          dhdC[w] = dXdC[w]*(mp->u_species_source[w][4]-1.0)*k_propX*(2./k_propX_num
                          +log((1.-Xconv)/(1.0-Xconv_init))/k_propX_den);
          dhdC[w] *= k_prop*Conc[w];  dhdC[w] += k_propT;
          dhdC[w] *= free_rad*param[w]*mp->molecular_weight[w];

          dhdT += k_propX*k_prop*mp->u_species_source[w][2]/SQUARE(fv->T)
                *Conc[w]*free_rad*param[w]*mp->molecular_weight[w];
          dhdI += k_propT*d_free_rad_dI*param[w]*mp->molecular_weight[w];

          if(model_bit & 2)
               { 
                dhdC[rad_spec] += k_propT*Conc[w]*param[w]*mp->molecular_weight[w];
               }
          else if(model_bit & 1)
               { 
                dhdC[O2_spec] += k_propT*Conc[w]*param[w]*
                      mp->molecular_weight[w]*(SQUARE(k_inh/2.)*Conc[O2_spec]/
                       sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
                mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec])
                       -k_inh/2.);
                dhdC[init_spec] += k_propT*Conc[w]*param[w]*
                       mp->molecular_weight[w]*
                       0.5*mp->u_species_source[init_spec+1][2]*intensity/
                       sqrt(SQUARE(k_inh*Conc[O2_spec])/4.+
           mp->u_species_source[init_spec+1][2]*intensity*Conc[init_spec]);

               }
          else 
               { 
                dhdC[init_spec] += k_propT*Conc[w]*param[w]*
                       mp->molecular_weight[w]*
                       0.5*sqrt(mp->u_species_source[init_spec+1][2]*intensity/
                       Conc[init_spec]);
               }
         }

/**  add heat generation from light absorption  **/

      h +=  param[1]*intensity*Conc[init_spec];
      dhdC[init_spec] += param[1]*intensity;
      dhdI =  param[1]*Conc[init_spec];

      var = TEMPERATURE;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->T[j]= dhdT*bf[var]->phi[j];
            }
        }
      var = LIGHT_INTP;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->INT[j]= dhdI*bf[var]->phi[j]
                           *mp->u_species_source[init_spec][1];
            }
        }
      var = ACOUS_PREAL;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->APR[j]= dhdI*bf[var]->phi[j]
              *mp->u_species_source[init_spec][1]*intensity_cgs*2.0*fv->apr;
            }
        }
      var = ACOUS_PIMAG;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->API[j]= dhdI*bf[var]->phi[j]
              *mp->u_species_source[init_spec][1]*intensity_cgs*2.0*fv->api;
            }
        }

      if ( d_h != NULL && pd->v[MASS_FRACTION] )
	{
	  for ( w=0; w<pd->Num_Species_Eqn; w++)
	    {
	      var = MASS_FRACTION;
	      var_offset = MAX_VARIABLE_TYPES + w;
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  d_h->C[w][j] =dhdC[w]*bf[var]->phi[j];
		}
	    }
	}

    }
  else if (mp->HeatSourceModel == CONSTANT )
    {
      h = mp->heat_source;
    }
  else if (mp->HeatSourceModel == ELECTRODE_KINETICS )  /* added by KSC/GHE: 10/20/98 */
    {
      electrolyte_temperature(time, dt, 0);  /* calculate electrolyte temperature at the present time */
      TT = mp->electrolyte_temperature;
      h = 0.0;
      for(w=0; w<pd->Num_Species_Eqn; w++)
        {
          electrode_species_source(w, time, dt);
          h -= F*mp->charge_number[w]*mp->species_source[w];
        }

      ai0_anode = mp->u_reaction_rate[0];
      ai0_cathode = mp->u_reaction_rate[2];
      mn = ei->mn;

      if(mn == 0)             /* KSC: 2/17/99 */
        {
	  ai0 = ai0_anode;     /* exchange current density for the anode (one rxn only) */
        }
      else if (mn == 2)
        {
	  ai0 = ai0_cathode;   /* exchange current density for the cathode (one rxn only) */
        }
      else
        {
	  ai0 = 0.0;           /* no electrochemical rxn for the separator */
        }

      var = TEMPERATURE;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->T[j]= -F*mp->charge_number[w]*(ai0/R/TT)*bf[var]->phi[j];
            }
        }

      var = VOLTAGE;
      if ( d_h != NULL && pd->e[var] )
        {
          for ( j=0; j<ei->dof[var]; j++)
            {
              d_h->V[j]= -F*mp->charge_number[w]*(-ai0/R/TT)*bf[var]->phi[j];
            }
        }
    }
  else if (mp->HeatSourceModel == BUTLER_VOLMER)      /* added by KSC: 4/28/06 */
    {
      h = butler_volmer_heat_source(d_h, mp->u_heat_source);
    }
  else if (mp->HeatSourceModel == JOULE )
    {
      h = joule_heat_source(d_h, tt);
    }
  else if (mp->HeatSourceModel == VISC_DISS )
    {
      h = visc_diss_heat_source(d_h, mp->u_heat_source);
    }
  else if (mp->HeatSourceModel == VISC_ACOUSTIC )
    {
      h = visc_diss_heat_source(d_h, mp->u_heat_source);
      h_acous = visc_diss_acoustic_source(d_h, mp->u_heat_source, mp->len_u_heat_source);
      h += h_acous;
    }
  else if (mp->HeatSourceModel == EPOXY )
    {
      h = epoxy_heat_source(d_h, tt, dt);
    }
  else if (mp->HeatSourceModel == VARIABLE_DENSITY )
    {
      h = vary_rho_heat_source(d_h, tt, dt);
    }
  else if (mp->HeatSourceModel == HS_FOAM )
  {
      h = foam_heat_source(d_h, tt, dt);
  }  
  else if (mp->HeatSourceModel == USER_GEN )
    {
      if ( d_h == NULL )
        {
          dbl dhdT[MDE];
          dbl dhdX[DIM][MDE];
          dbl dhdv[DIM][MDE];
          dbl dhdC[MAX_CONC][MDE];
          dbl dhdVolt[MDE];
          usr_heat_source_gen(&h, dhdT, dhdX, dhdv, dhdC, dhdVolt, mp->u_heat_source, time);
        }
      else
        {
          usr_heat_source_gen(&h, d_h->T, d_h->X, d_h->v, d_h->C, d_h->V, mp->u_heat_source, time);
        }
    }
  else
    {
      EH(-1,"Unrecognized heat source model");
    }
	
	if( ls != NULL &&
	  mp->mp2nd != NULL &&
	  mp->mp2nd->HeatSourceModel == CONSTANT &&
	  (pd->e[R_ENERGY] & T_SOURCE) )
		{
		/* kludge for solidification tracking with phase function 0 */
      		if(pfd != NULL && pd->e[R_EXT_VELOCITY])
		{
		ls_old = ls;
		ls = pfd->ls[0];
        	ls_modulate_heatsource( &h,
				  mp->mp2nd->heatsource_phase[0],
				  ls->Length_Scale,
				  (double) mp->mp2nd->heatsourcemask[0],
				  (double) mp->mp2nd->heatsourcemask[1],
				  d_h );
		ls = ls_old;
		}
	  ls_modulate_heatsource ( &h,
				mp->mp2nd->heatsource,
			   ls->Length_Scale,
			   (double) mp->mp2nd->heatsourcemask[0],
			   (double) mp->mp2nd->heatsourcemask[1],
			   d_h );
		}

  return ( h );
}

int
ls_modulate_heatsource(double *f,
		       double f2,
		       double width,
		       double pm_minus,
		       double pm_plus,
		       HEAT_SOURCE_DEPENDENCE_STRUCT *df )
{
  int i, b, var;
  int dim = pd->Num_Dim;
  double factor;
  double f1 = *f;
	
  if ( df == NULL )
    {
      *f = ls_modulate_property( f1,
				 f2,
				 width,
				 pm_minus,
				 pm_plus,
				 NULL ,
				 &factor );
		
      return (0);
    }
  else
    {
      *f = ls_modulate_property( f1,
				 f2,
				 width,
				 pm_minus,
				 pm_plus,
				 df->F ,
				 &factor );
		
		
		
      if( pd->v[var=TEMPERATURE] )
	{
	  for ( i=0; i<ei->dof[var]; i++)
	    {
	      df->T[i] *=factor;
	    }
	}
		
      if(  pd->v[var=MESH_DISPLACEMENT1] )
	{
	  for( b=0; b<dim; b++ )
	    {
	      for( i=0 ; i<ei->dof[var+b]; i++)
		{
		  df->X[b][i] *= factor;
		}
	    }
	}
		
      if( pd->v[var=VELOCITY1] )
	{
	  for( b=0; b<dim; b++ )
	    {
	      for( i=0 ; i<ei->dof[var]; i++)
		{
		  df->v[b][i] *= factor;
		}
	    }
	}
		
      if(  pd->v[var=MASS_FRACTION] )
	{
	  for( b=0; b<pd->Num_Species; b++ )
	    {
	      for( i=0 ; i<ei->dof[var]; i++)
		{
		  df->C[b][i] *= factor;
		}
	    }
	}
		
      if(  pd->v[var=VOLTAGE] )
	{
	  for ( i=0; i<ei->dof[var]; i++)
            {
	      df->V[i] *= factor;
            }
	} 
		
		
      if(  pd->v[var=POLYMER_STRESS11] )
	{
	  WH(-1,"LS modulation of heat source sensitivity wrt to polymer stress dofs not implemented.");
	}
    }
  return(0);
}

int
calc_pspg( dbl pspg[DIM],
	   PSPG_DEPENDENCE_STRUCT *d_pspg,
	   dbl time_value, /* current time */
	   dbl tt,	/* parameter to vary time integration from
					   explicit (tt = 1) to implicit (tt = 0)    */
	   dbl dt,	/* current time step size                    */
	   const PG_DATA *pg_data )
{
  const dbl h_elem_avg = pg_data->h_elem_avg;
  const dbl *hsquared = pg_data->hsquared;      	/* element size information for PSPG         */
  const dbl U_norm = pg_data->U_norm;	                /* global velocity norm for PSPG calcs       */
  const dbl mu_avg = pg_data->mu_avg;	                /* element viscosity for PSPG calculations   */
  const dbl rho_avg = pg_data->rho_avg;                 /* element density for PSPG calculations   */
  const dbl *v_avg = pg_data->v_avg;                    /* element velocity for PSPG calculations   */

  int dim, wim;
  int p, a, b, c;
  int var;
  int w, j;

  int pspg_global;
  int pspg_local;

  dbl *v = fv->v;				/* Velocity field. */
  dbl *grad_P = fv->grad_P;

  /*
   * Variables for vicosity and derivative
   */
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /*
   * density and sensitivity terms
   */
  dbl rho;
  dbl rho_t;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  dbl f[DIM];				/* Body force. */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct;  /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  /*
   * Interpolation functions...
   */

  dbl phi_j;

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn, meqn1;
  int r;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];

  dbl mass, advection, diffusion, source, porous;

  dbl momentum[DIM];           /* momentum residual for PSPG */
  dbl x_dot[DIM];
  dbl v_dot[DIM];
  dbl *grad_v[DIM];
  dbl div_s[DIM];
  dbl div_G[DIM];

  /* variables for Brinkman porous flow */
  dbl por=0, por2=0, per=0, vis=0, dvis_dT[MDE], sc=0, speed=0;

  dbl h_elem=0;
  dbl Re;
  dbl tau_pspg=0.;
  dbl tau_pspg1=0.;
  dbl d_tau_pspg_dX[DIM][MDE];
  dbl d_tau_pspg_dv[DIM][MDE];

  dbl hh_siz, vv_speed;

  dbl gamma[DIM][DIM];                  /* shrearrate tensor based on velocity */

  int w0=0;

  int mode;

  /* For particle momentum model.
   */
  int species;			/* species number for particle phase,  */
  dbl ompvf;			/* 1 - partical volume fraction */
  int particle_momentum_on; 	/* boolean. */
  /* particle stress for suspension balance model*/
  dbl tau_p[DIM][DIM];
  dbl d_tau_p_dv[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dvd[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dy[DIM][DIM][MAX_CONC][MDE];
  dbl d_tau_p_dmesh[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dp[DIM][DIM][MDE];
  
  static int is_initialized=FALSE;

  dim   = pd->Num_Dim;
  wim   = dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN) wim = wim+1;

  /* initialize */
  for ( a=0; a<DIM; a++ ) pspg[a] = 0.;

  if ( d_pspg == NULL )
    {
      /* I guess we won't be needing these! */
      d_mu = NULL;
      d_rho = NULL;
      df = NULL;
    }
  else if ( !is_initialized )
    {

      memset(d_pspg->v, 0, sizeof(double)*DIM*DIM*MDE);
      memset(d_pspg->X, 0, sizeof(double)*DIM*DIM*MDE);
      memset(d_pspg->T, 0, sizeof(double)*DIM*MDE);
      memset(d_pspg->P, 0, sizeof(double)*DIM*MDE);
      memset(d_pspg->C, 0, sizeof(double)*DIM*MAX_CONC*MDE);
      memset(d_pspg->S, 0, sizeof(double)*DIM*MAX_MODES*DIM*DIM*MDE);
      memset(d_pspg->g, 0, sizeof(double)*DIM*DIM*DIM*MDE);
    }

  /* This is the flag for the standard global PSPG */
  if(PSPG == 1)
    {
      pspg_global = TRUE;
      pspg_local = FALSE;
    }
  /* This is the flag for the standard local PSPG */
  else if(PSPG == 2)
    {
      pspg_global = FALSE;
      pspg_local = TRUE;
    }
  else
    {
      return 0;
    }

  if( pd->v[POLYMER_STRESS11] )
  {	
	stress_eqn_pointer(v_s);

	v_g[0][0] = VELOCITY_GRADIENT11;
	v_g[0][1] = VELOCITY_GRADIENT12;
	v_g[1][0] = VELOCITY_GRADIENT21;
	v_g[1][1] = VELOCITY_GRADIENT22;
	v_g[0][2] = VELOCITY_GRADIENT13;
	v_g[1][2] = VELOCITY_GRADIENT23;
	v_g[2][0] = VELOCITY_GRADIENT31;
	v_g[2][1] = VELOCITY_GRADIENT32;
	v_g[2][2] = VELOCITY_GRADIENT33;
  }

  /* initialize dependencies */
  memset( d_tau_pspg_dv, 0, sizeof(double) * DIM*MDE);
  memset( d_tau_pspg_dX, 0, sizeof(double) * DIM*MDE);
  

  if( cr->MassFluxModel == DM_SUSPENSION_BALANCE && PSPG)
    {
      w0 = gn->sus_species_no;
      /* This is the divergence of the particle stress  */
      /* divergence_particle_stress(div_tau_p, d_div_tau_p_dgd, d_div_tau_p_dy,
	 d_div_tau_p_dv, d_div_tau_p_dX, w0); */
    }


  if(pd->e[R_PMOMENTUM1])
    {
      particle_momentum_on = 1;
      species = (int) mp->u_density[0];
      ompvf = 1.0 - fv->c[species];
    }
  else
    {
      particle_momentum_on = 0;
      species = -1;
      ompvf = 1.0;
    }


  // Global average for pspg_global's element size
  h_elem = h_elem_avg;

  /*** Density ***/
  rho = density(d_rho, time_value);

  if(pspg_global)
    {
      
      /* Now calculate the element Reynolds number based on a global
       * norm of the velocity and determine tau_pspg discretely from Re 
       * The global version has no Jacobian dependencies
       */
      Re = rho * U_norm * h_elem / (2.0 * mu_avg);

      if (Re <= 3.0)
	{
	  tau_pspg = PS_scaling * h_elem * h_elem / (12.0 * mu_avg);
	}
      else if (Re > 3.0)
	{
	  tau_pspg = PS_scaling * h_elem / (2.0 * rho * U_norm);
	}
    }
  else if (pspg_local)
    {
      hh_siz = 0.;
      for ( p=0; p<dim; p++)
	{
	  hh_siz += hsquared[p];
	}
      // Average value of h**2 in the element
      hh_siz = hh_siz/ ((double )dim);

      // Average value of v**2 in the element 
      vv_speed = 0.0;
      for ( a=0; a<wim; a++)
	{
	  vv_speed += v_avg[a]*v_avg[a];
	}
      
      // Use vv_speed and hh_siz for tau_pspg, note it has a continuous dependence on Re
      tau_pspg1 = rho_avg*rho_avg*vv_speed/hh_siz + (9.0*mu_avg*mu_avg)/(hh_siz*hh_siz);
      if (  pd->TimeIntegration != STEADY)
	{
	  tau_pspg1 += 4.0/(dt*dt);
	}
      tau_pspg = PS_scaling/sqrt(tau_pspg1);

      // tau_pspg derivatives wrt v from vv_speed
      if ( d_pspg != NULL && pd->v[VELOCITY1] )
	{
	  for ( b=0; b<dim; b++)
	    {
	      var = VELOCITY1+b;
	      if ( pd->v[var] )
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      d_tau_pspg_dv[b][j] = -tau_pspg/tau_pspg1; 
		      d_tau_pspg_dv[b][j] *= rho_avg*rho_avg/hh_siz * v_avg[b]*pg_data->dv_dnode[b][j];
		    }
		}
	    }
	}

      // tau_pspg derivatives wrt mesh from hh_siz
      if ( d_pspg != NULL && pd->v[MESH_DISPLACEMENT1] )
	{
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {			
		      d_tau_pspg_dX[b][j] = tau_pspg/tau_pspg1;
		      d_tau_pspg_dX[b][j] *= (rho_avg*rho_avg*vv_speed + 18.0*mu_avg*mu_avg/hh_siz) / (hh_siz*hh_siz);
		      d_tau_pspg_dX[b][j] *= pg_data->hhv[b][b]*pg_data->dhv_dxnode[b][j]/((double)dim);
		      
		    }
		}
	    }
	}
    }

  for ( a=0; a<VIM; a++) grad_v[a] = fv->grad_v[a];
  
  /* load up shearrate tensor based on velocity */
  for ( a=0; a<VIM; a++)
    {
      for ( b=0; b<VIM; b++)
	{
	  gamma[a][b] = grad_v[a][b] + grad_v[b][a];
	}
    }

  /*
   * get viscosity for velocity second derivative/diffusion
   * term in PSPG stuff
   */
  mu = viscosity(gn, gamma, d_mu );


  /* get variables we will need for momentum residual */

  for ( a=0; a<wim; a++)
    {
      if (  pd->TimeIntegration != STEADY &&  pd->v[MESH_DISPLACEMENT1+a] )
	{
	  x_dot[a] = fv_dot->x[a];
	}
      else
	{
	  x_dot[a] = 0.;
	}

      if (  pd->TimeIntegration != STEADY )
	{
	  v_dot[a] = fv_dot->v[a];
	}
      else
	{
	  v_dot[a] = 0.;
	}
    }

  for ( p=0; p<wim; p++) div_s[p] = 0.;
  if ( pd->v[POLYMER_STRESS11] )
    {
      for ( p=0; p<wim; p++)
	{
	  for ( mode=0; mode<vn->modes; mode++)
	    {
	      div_s[p] += fv->div_S[mode][p];
	    }
	}
    }


	if( !is_initialized ) {
	memset( tau_p,      0, sizeof(double) * DIM*DIM);
	memset( d_tau_p_dv, 0, sizeof(double) * DIM*DIM*DIM*MDE);
	memset( d_tau_p_dvd, 0, sizeof(double) * DIM*DIM*DIM*MDE);
	memset( d_tau_p_dy, 0, sizeof(double) * DIM*DIM*MAX_CONC*MDE);
	memset( d_tau_p_dmesh, 0, sizeof(double) * DIM*DIM*DIM*MDE);
	memset( d_tau_p_dp, 0, sizeof(double) * DIM*DIM*MDE);
	}

  if( cr->MassFluxModel == DM_SUSPENSION_BALANCE
      || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD
      || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR)
    particle_stress(tau_p,d_tau_p_dv,d_tau_p_dvd,d_tau_p_dy,d_tau_p_dmesh,d_tau_p_dp, w0);

  if ( pd->v[VELOCITY_GRADIENT11] &&  pd->v[POLYMER_STRESS11])
    {
      for ( p=0; p<wim; p++)
	{
	  div_G[p] = fv->div_G[p];
	}
    }
  else
    {
      for ( p=0; p<wim; p++)
	{
	  div_G[p] = 0.;
	}
    }

  if (pd->e[R_MOMENTUM1] & T_POROUS_BRINK )
    {
      if (mp->PorousMediaType != POROUS_BRINKMAN)
	WH(-1, "Set Porous term multiplier in continuous medium");
      /* Short-hand notation for the four parameters in the Brinkman Equation. */
      por = mp->porosity;
      por2 = por * por;
      per = mp->permeability;
      if(mp->FlowingLiquidViscosityModel == CONSTANT)
	{
	  /* Do nothing */
	}
      else if (mp->FlowingLiquidViscosityModel == MOLTEN_GLASS)
	{
	  molten_glass_viscosity(&(mp->FlowingLiquid_viscosity), dvis_dT, mp->u_FlowingLiquid_viscosity);
	}
      else if (mp->FlowingLiquidViscosityModel == USER)
	{
	  usr_FlowingLiquidViscosity(mp->u_FlowingLiquid_viscosity);
	  var = TEMPERATURE;

	  for ( j=0; j<ei->dof[var]; j++)
	    {
	      dvis_dT[j]= mp->d_FlowingLiquid_viscosity[var]*bf[var]->phi[j];
	    }
	}
      else
	{
	  EH(-1,"Don't recognize your FlowingLiquidViscosity model");
	}
      vis = mp->FlowingLiquid_viscosity;
      sc  = mp->Inertia_coefficient;
    }
  else
    {
      por = 1.;
      por2 = 1.;
      per = 1.;
      vis = mp->viscosity;
      sc  = 0.;
    }

  /* for porous media stuff */
  speed = 0.0;
  for ( a=0; a<wim; a++)
    {
      speed += v[a]*v[a];
    }
  speed = sqrt(speed);

  /* get momentum source term */
  momentum_source_term(f, df, time_value);

  if(pd->e[R_PMOMENTUM1])
    {
      rho_t = ompvf * rho;
      meqn1 = R_PMOMENTUM1;
    }
  else
    {
      rho_t = rho;
      meqn1 = R_MOMENTUM1;
    }

  for ( a=0; a<wim; a++)
  {
      meqn = meqn1 + a;
	  
      mass = 0.;
      if ( ( pd->e[meqn] & T_MASS ) && ( pd->TimeIntegration != STEADY ) )
	  {
		  mass =  rho_t * v_dot[a] / por;
		  mass *= pd->etm[meqn][(LOG2_MASS)];
	  }
	  
      advection = 0.;
      if ( pd->e[meqn] & T_ADVECTION )
	  {
		  for ( p=0; p<wim; p++)
		  {
			  advection += rho_t * (v[p]- x_dot[p]) * grad_v[p][a] / por2;
		  }
		  advection *= pd->etm[meqn][(LOG2_ADVECTION)];
	  }
	  
      diffusion = 0.;
      if ( pd->e[meqn] & T_DIFFUSION )
	  {
		  diffusion = grad_P[a] - div_s[a]; 
		  /*diffusion  -= div_tau_p[a]  */
		   diffusion -= mu * div_G[a];
		  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
	  }
	  
      source = 0.;
      if ( pd->e[meqn] & LOG2_SOURCE )
	  {
		  source = -f[a];
		  source *= pd->etm[meqn][(LOG2_SOURCE)];
	  }
	  
      porous = 0.;
      if ( pd->e[meqn] & LOG2_POROUS_BRINK )
	  {
		  porous = v[a]*(rho_t*sc*speed/sqrt(per)+vis/per);
		  porous *= pd->etm[meqn][(LOG2_POROUS_BRINK)];
	  }
	  
      momentum[a] = mass + advection + diffusion + source + porous;
      pspg[a] = tau_pspg * momentum[a];
	  
	  
      if ( d_pspg != NULL && pd->v[VELOCITY1] )
	  {
		  for ( b=0; b<wim; b++ )
		  {
              var = VELOCITY1 + b;
			  for ( j=0; j<ei->dof[var]; j++ )
			  {
				  phi_j = bf[var]->phi[j];
				  
				  mass = 0.;
				  if ( ( pd->e[meqn] & T_MASS ) && ( pd->TimeIntegration != STEADY ) )
				  {
					  mass =  rho_t / por * (1.+2.*tt) * phi_j/dt * (double)delta(a,b);
					  mass *= pd->etm[meqn][(LOG2_MASS)];
				  }
				  
				  advection = 0.;
				  if ( pd->e[meqn] & T_ADVECTION )
				  {
					  advection = phi_j * grad_v[b][a];
					  for ( p=0; p<wim; p++)
					  {
						  advection += (v[p]- x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
					  }
					  advection *= rho_t / por2;
					  advection *= pd->etm[meqn][(LOG2_ADVECTION)];
				  }
				  
				  diffusion = 0.;
				  if ( pd->e[meqn] & T_DIFFUSION )
				  {
					  diffusion -= d_mu->v[b][j] * div_G[a];
					  /*diffusion -= d_div_tau_p_dv[a][b][j];*/
					  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
				  }
				  
				  source = 0.;
				  if ( pd->e[meqn] & T_SOURCE )
				  {
					  source -= df->v[a][b][j] * pd->etm[meqn][(LOG2_SOURCE)];
				  }
				  
				  porous = 0.;
				  if ( pd->e[meqn] & T_POROUS_BRINK )
				  {
					  porous = (rho_t*sc*speed/sqrt(per) + vis/per) * phi_j;
					  for ( p=0; p<wim; p++)
					  {
						  porous += rho_t*sc/sqrt(per) *  2.*v[p] * v[a];
					  }
					  porous *= pd->etm[meqn][(LOG2_POROUS_BRINK)];
				  }
				  d_pspg->v[a][b][j] = tau_pspg * ( mass + advection + diffusion + source + porous )
					  + d_tau_pspg_dv[b][j] * momentum[a];
				  
			  }
		  }
	  }
	  
      if ( d_pspg != NULL && pd->v[MESH_DISPLACEMENT1] )
	  {
		  for ( b=0; b<wim; b++ )
		  {
			  var = MESH_DISPLACEMENT1 + b;
			  for ( j=0; j<ei->dof[var]; j++ )
			  {
				  phi_j = bf[var]->phi[j];
				  
				  advection = 0.;
				  if  ( ( pd->e[meqn] & T_ADVECTION ) && ( pd->TimeIntegration != STEADY ) )
				  {
					  advection = -(1.+2.*tt) * phi_j/dt * grad_v[b][a];
					  for ( p=0; p<wim; p++)
					  {
						  advection += (v[p]- x_dot[p]) * fv->d_grad_v_dmesh[p][a] [b][j];
					  }
					  advection *= rho_t / por2;
					  advection *= pd->etm[meqn][(LOG2_ADVECTION)];
				  }
				  
				  diffusion = 0.;
				  if ( pd->e[meqn] & T_DIFFUSION )
				  {
					  diffusion = fv->d_grad_P_dmesh[a][b][j]
					  - d_mu->X[b][j] * div_G[a] 
					  - mu * fv->d_div_G_dmesh[a] [b][j] ;
					 /* diffusion -= d_div_tau_p_dX[a][b][j];*/
					  for ( mode=0; mode<vn->modes; mode++)
					  {
						  diffusion -= fv->d_div_S_dmesh[mode][a] [b][j];
					  }
					  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
				  }
				  
				  source = 0.;
				  if ( pd->e[meqn] & T_SOURCE )
				  {
					  source -= df->X[a][b][j] * pd->etm[meqn][(LOG2_SOURCE)];
				  }
				  
				  d_pspg->X[a][b][j] = tau_pspg * ( advection + diffusion + source ) +
					  + d_tau_pspg_dX[b][j] * momentum[a];
				  
			  }
		  }
	  }
	  
      var = TEMPERATURE;
      if ( d_pspg != NULL && pd->v[var] )
	  {
		  double d_rho_t_dT;
		  
		  for ( j=0; j<ei->dof[var]; j++ )
		  {
			  phi_j = bf[var]->phi[j];
			  
			  d_rho_t_dT = ompvf * d_rho->T[j];
			  
			  mass = 0.;
			  if ( ( pd->e[meqn] & T_MASS ) && ( pd->TimeIntegration != STEADY ) )
			  {
				  mass = d_rho_t_dT/por * v_dot[a];
				  mass *= pd->etm[meqn][(LOG2_MASS)];
			  }
			  
			  advection = 0.;
			  if ( pd->e[meqn] & T_ADVECTION )
			  {
				  advection = 0.;
				  for ( p=0; p<wim; p++)
				  {
					  advection += (v[p]- x_dot[p]) * grad_v[p][a];
				  }
				  advection *= d_rho_t_dT / por2;
				  advection *= pd->etm[meqn][(LOG2_ADVECTION)];
			  }
			  
			  diffusion = 0.;
			  if ( pd->e[meqn] & T_DIFFUSION )
			  {
				  diffusion -= d_mu->T[j] * div_G[a];
				  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
			  }
			  
			  source = 0.;
			  if ( pd->e[meqn] & T_SOURCE )
			  {
				  source -= df->T[a][j] * pd->etm[meqn][(LOG2_SOURCE)];
			  }
			  
			  porous = 0.;
			  if ( pd->e[meqn] & T_POROUS_BRINK )
			  {
				  porous = v[a] * (d_rho_t_dT*sc*speed/sqrt(per) + dvis_dT[j]/per);
				  porous *= pd->etm[meqn][(LOG2_POROUS_BRINK)];
			  }
			  d_pspg->T[a][j] = tau_pspg * ( mass + advection + diffusion + source + porous );
			  
			  
		  }
	  }
	  
      var = PRESSURE;
      if ( d_pspg != NULL && pd->v[var] )
	  {
		  for ( j=0; j<ei->dof[var]; j++ )
		  {
			  diffusion = 0.;
			  if ( pd->e[meqn] & T_DIFFUSION )
			  {
				  diffusion = bf[var]->grad_phi[j][a]
				  - d_mu->P[j] * div_G[a];
				  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
			  }
			  
			  d_pspg->P[a][j] = tau_pspg * ( diffusion );
			  
		  }
	  }
	  
      var = MASS_FRACTION;
      if ( d_pspg != NULL && pd->v[var] )
	  {
		  double d_rho_t_dC;
		  
		  for (w=0; w<pd->Num_Species_Eqn; w++)
		  {
			  for ( j=0; j<ei->dof[var]; j++ )
			  {
				  phi_j = bf[var]->phi[j];
				  
                  d_rho_t_dC = ompvf * d_rho->C[w][j];
				  if ( particle_momentum_on && w==species )
					  d_rho_t_dC -= phi_j * rho;
				  
				  mass = 0.;
				  if ( ( pd->e[meqn] & T_MASS ) && ( pd->TimeIntegration != STEADY ) )
				  {
					  mass = d_rho_t_dC/por * v_dot[a];
					  mass *= pd->etm[meqn][(LOG2_MASS)];
				  }
				  
				  advection = 0.;
				  if ( pd->e[meqn] & T_ADVECTION )
				  {
					  advection = 0.;
					  for ( p=0; p<wim; p++)
					  {
						  advection += (v[p]- x_dot[p]) * grad_v[p][a];
					  }
					  advection *= d_rho_t_dC / por2;
					  advection *= pd->etm[meqn][(LOG2_ADVECTION)];
				  }
				  
				  diffusion = 0.;
				  if ( pd->e[meqn] & T_DIFFUSION )
				  {
					  diffusion -= d_mu->C[w][j] * div_G[a] ;
					  /*diffusion  -= d_div_tau_p_dy[a][w][j]; */
					  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
				  }
				  
				  source = 0.;
				  if ( pd->e[meqn] & T_SOURCE )
				  {
					  source -= df->C[a][w][j] * pd->etm[meqn][(LOG2_SOURCE)];
				  }
				  
				  porous = 0.;
				  if ( pd->e[meqn] & T_POROUS_BRINK )
				  {
					  porous = v[a] * (d_rho_t_dC*sc*speed/sqrt(per));
					  porous *= pd->etm[meqn][(LOG2_POROUS_BRINK)];
				  }
				  d_pspg->C[a][w][j] = tau_pspg * ( mass + advection + diffusion + source + porous );
				  
			  }
		  }
	  }
	  
      var = POLYMER_STRESS11;
      if ( d_pspg != NULL && pd->v[var] )
	  {
		  for ( mode=0; mode<vn->modes; mode++)
		  {
			  for ( b=0; b<VIM; b++)
			  {
				  for ( c=0; c<VIM; c++)
				  {
					  var = v_s[mode][b][c];
					  for ( j=0; j<ei->dof[var]; j++ )
					  {
						  phi_j = bf[var]->phi[j];
						  
						  diffusion = 0.;
						  if ( pd->e[meqn] & T_DIFFUSION )
						  {
							  diffusion = -(double)delta(a,c) * bf[var]->grad_phi[j] [b];
							  
							  if ( pd->CoordinateSystem != CARTESIAN )
							  {
								  for ( r=0; r<VIM; r++)
								  {
									  diffusion -= (double)delta(a,c) * phi_j * fv->grad_e[b][r][c];
								  }
								  for ( r=0; r<wim; r++)
								  {
									  diffusion -= (double)delta(a,r) * phi_j * fv->grad_e[c][b][r];
								  }
							  }
							  
							  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
						  }
						  
						  d_pspg->S[a][mode][b][c][j] = tau_pspg * diffusion;
						  
					  }
				  }
			  }
		  }
	  }
	  
      var = VELOCITY_GRADIENT11;
      if ( d_pspg != NULL && pd->v[var] )
	  {
		  for ( b=0; b<VIM; b++)
		  {
			  for ( c=0; c<VIM; c++)
			  {
				  var = v_g[b][c];
				  for ( j=0; j<ei->dof[var]; j++ )
				  {
					  phi_j = bf[var]->phi[j];
					  
					  diffusion = 0.;
					  if ( pd->e[meqn] & T_DIFFUSION )
					  {
						  diffusion = -(double)delta(a,c) * bf[var]->grad_phi[j] [b];
						  
						  if ( pd->CoordinateSystem != CARTESIAN )
						  {
							  for ( r=0; r<VIM; r++)
							  {
								  diffusion -= (double)delta(a,c) * phi_j * fv->grad_e[b][r][c];
							  }
							  for ( r=0; r<wim; r++)
							  {
								  diffusion -= (double)delta(a,r) * phi_j * fv->grad_e[c][b][r];
							  }
						  }
						  diffusion *= mu;
						 /* diffusion -= d_div_tau_p_dgd0[a][b][c][j]; */
						  diffusion *= pd->etm[meqn][(LOG2_DIFFUSION)];
					  }
					  
					  d_pspg->g[a][b][c][j] = tau_pspg * diffusion;
					  
				  }
			  }
		  }
	  }
  }
is_initialized = TRUE; 
return 0;
}


int
calc_cont_gls( dbl *cont_gls,
                CONT_GLS_DEPENDENCE_STRUCT *d_cont_gls,
		dbl time_value,                          /*Current time value, needed for density */
 	        const PG_DATA *pg_data)                  /*Petrov-Galerkin data, needed for tau_cont   */
{
  dbl div_v = fv->div_v;
  dbl tau_cont, d_tau_dmesh[DIM][MDE], d_tau_dv[DIM][MDE];
  dbl advection_etm, advection, Re, div_phi_j_e_b, div_v_dmesh;
  int eqn; 
  int var, dim, wim, b, j, p;
  int advection_on=0;
  static int is_initialized = FALSE;

  //Density terms
  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  //Petrov-Galerkin values
  dbl h_elem = 0, U, hh_siz, vv;
  const dbl h_elem_avg = pg_data->h_elem_avg;
  const dbl *hsquared = pg_data->hsquared; 
  const dbl U_norm = pg_data->U_norm;
  const dbl *v_avg = pg_data->v_avg;
  const dbl mu_avg = pg_data->mu_avg;
  
  //Initialize and Define
  dim = pd->Num_Dim;
  wim=dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  tau_cont = 0;
  *cont_gls = 0.0;
  eqn = R_PRESSURE;
  advection_on = pd->e[eqn] & T_ADVECTION;
  advection_etm = pd->etm[eqn][(LOG2_ADVECTION)];

  if (d_cont_gls==NULL)
    {
      d_rho = NULL;
    }
  else if(!is_initialized)
    {
      memset(d_cont_gls->v, 0, sizeof(double)*DIM*MDE);
      memset(d_cont_gls->X, 0, sizeof(double)*DIM*MDE);
    }

  /* Calculate stabilization parameter tau_cont
   * From Wall:
   * tau_cont = h_elem*u_norm*max{Re,1}/2;
   * The Reynolds number, Re, is defined as in calc_pspg
   */
  if(Cont_GLS==1)
    {
      h_elem = h_elem_avg;
      U = U_norm;
    }
  else
    {
      hh_siz = 0.0;
      vv = 0.0;
      for(p=0; p<dim; p++)
	{
	  hh_siz += hsquared[p]/((double)dim);
	  vv += v_avg[p]*v_avg[p];
	}
      h_elem = sqrt(hh_siz);
      U = sqrt(vv);
    }
  
  tau_cont = U*h_elem /2.0;

  //We need density for the Reynolds number
  rho = density(d_rho, time_value);
  Re = rho*U*h_elem / (2.0*mu_avg);
  if (Re > 1.0)
    {
      tau_cont *= Re;
    }
  
  //d_tau terms for Jacobian
  if(d_cont_gls != NULL && pd->v[VELOCITY1])
    {
      for(b=0; b<dim; b++)
	{
	  var = VELOCITY1+b;
	  if(pd->v[var])
	    {
	      for(j=0; j<ei->dof[var]; j++)
		{
		  if(Cont_GLS==1 || U==0)
		    {
		      d_tau_dv[b][j] = 0.0;
		    }
		  else if(Re > 1.0)
		    {		      
		      d_tau_dv[b][j] = rho/(2.0*mu_avg) * h_elem*h_elem * v_avg[b]*pg_data->dv_dnode[b][j];			
		    }
		  else
		    {		     		
		      d_tau_dv[b][j] = h_elem/(2.0*U) * v_avg[b]*pg_data->dv_dnode[b][j];			
		    }
		}
	    }
	}
    }

  if(d_cont_gls != NULL && pd->v[MESH_DISPLACEMENT1])
    {
      for(b=0; b<dim; b++)
	{
	  var = MESH_DISPLACEMENT1+b;
	  if(pd->v[var])
	    {
	      for(j=0; j<ei->dof[var]; j++)
		{
		  if(Cont_GLS==1 || h_elem==0)
		    {
		      d_tau_dmesh[b][j] = 0.0;
		    }
		  else if(Re > 1.0)
		    {	      	
		      d_tau_dmesh[b][j] = rho/(2.0*mu_avg) * U*U * pg_data->hhv[b][b]*pg_data->dhv_dxnode[b][j]/((double)dim);		
		    }
		  else
		    {
		      d_tau_dmesh[b][j] = U/(2.0*h_elem) * pg_data->hhv[b][b]*pg_data->dhv_dxnode[b][j]/((double)dim);
		    }
		}
	    }
	}
    }

  
  /*
   * Calculate residual 
   * This term refers to the standard del dot v . 
   */
  advection = 0.0;
  if (advection_on)
    {
      if (pd->v[VELOCITY1])
	{
	  advection = div_v;
	  advection *= advection_etm;
	}
    }

  /*Multiply tau_cont by residual resulting in cont_gls
   *This eventually then gets combined with the stabilization functional
   *in assemble_momentum to form the GLS stabilization for continuity
   */
  *cont_gls = tau_cont*advection;


  //Determine Jacobian terms

  //J_v
  for (b=0; b<wim; b++)
    {
      var = VELOCITY1+b;
      if (pd->v[var] && d_cont_gls!=NULL)
	{	  
	for (j=0; j<ei->dof[var]; j++)
	  {				  
	    advection  = 0.;			  
	    if (advection_on)
	      {
		div_phi_j_e_b = 0.;
		for (p=0; p<VIM; p++)
		  {
		    div_phi_j_e_b += bf[var]->grad_phi_e[j][b] [p][p];
		  }
		advection = div_phi_j_e_b*advection_etm;
	      }
	    d_cont_gls->v[b][j] = tau_cont*advection + d_tau_dv[b][j]*div_v*advection_etm;
	  }
	}
    }


  //J_d 
  for (b=0; b<dim; b++)
    {
      var = MESH_DISPLACEMENT1+b;
      if (pd->v[var])
	{
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      advection = 0.0;
	      if (advection_on)
		{
		  if (pd->v[VELOCITY1] && d_cont_gls!=NULL)
		    {
		      div_v_dmesh = fv->d_div_v_dmesh[b][j];
		      advection += div_v_dmesh;
		    }
		}
	      advection *= advection_etm;
	      d_cont_gls->X[b][j] = tau_cont*advection + d_tau_dmesh[b][j]*div_v*advection_etm;
	    }
	}
    }


  is_initialized = TRUE;
  return 0;
}


int
assemble_ls_latent_heat_source ( double iso_therm,
				double latent_heat, 
		           double dt,	/* current value of the time step  */
		           double tt,
                           double time,
			   int bc_input_id,
	                   struct Boundary_Condition *BC_Types )
{
  int i,j, ii,ledof;
  int eqn, peqn, var, pvar, b, w;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j,flux;
  double source;
  
  double sign = 1.;
  double vnorm;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;

  /* struct Boundary_Condition *fluxbc; */
  
  eqn = R_ENERGY;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
    { det_J = fv->sdet; }
  else              /* diffuse interface */
    { det_J = bf[eqn]->detJ; }
    
  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;

  /*
  fluxbc = BC_Types + bc_input_id;
  compute_leak_velocity_heat(&vnorm, d_vnorm, tt, dt, NULL, fluxbc);  */

  vnorm = fv->ext_v;
  flux = latent_heat*vnorm;
  memset(d_vnorm->v, 0, sizeof(dbl)*DIM*MDE);
  memset(d_vnorm->T, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->C, 0, sizeof(dbl)*MAX_CONC*MDE);
  memset(d_vnorm->F, 0, sizeof(dbl)*MDE);
  memset(d_vnorm->X, 0, sizeof(dbl)*MDE*DIM);
  
  /*
  if ( fv->c[wspec] > 1. )
    fprintf(stderr,"flux=%g, Y_w=%g, mass_flux=%g, extv=%g, vnorm=%g\n",flux,fv->c[wspec],-mp->mass_flux[wspec],fv->ext_v,vnorm);
  */

#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      peqn = TEMPERATURE;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
	      
              /* J_m_F
	       */
              var = LS ;
	      pvar = upd->vp[var];
		      
	      for( j=0; j<ei->dof[var]; j++)
	        {
	          phi_j = bf[var]->phi[j];
	          lec->J[peqn][pvar][ii][j] += source * phi_j;
	        }
	    }
	}
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = TEMPERATURE;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      source = phi_i * flux * lsi->delta;
	      
	      source *= det_J * wt;
	      
	      source *= h3;
	      
	      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
	      
	      lec->R[peqn][ii] += source;
	    }
	}
    }
  
  if ( af->Assemble_Jacobian )
    {
      peqn = TEMPERATURE;
      bfm = bf[eqn];
      
      for (i = 0; i < ei->dof[eqn]; i++)
	{
	  ledof = ei->lvdof_to_ledof[eqn][i];
	  
	  if (ei->active_interp_ledof[ledof])
	    {
	      ii = ei->lvdof_to_row_lvdof[eqn][i];
	      
	      phi_i = bfm->phi[i];
	      
	      /*
	       * J_W_T
	       */
	      var = TEMPERATURE;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
		      source = latent_heat * d_vnorm->T[j] * sign;
		      
		      source *= phi_i * lsi->delta;
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
                
              /*
	       * J_W_vext
	       */
	      var = EXT_VELOCITY;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
		      source = latent_heat * phi_j * sign;
		      
		      source *= phi_i * lsi->delta;
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}

              /*
	       * J_W_w
	       */
	      var = MASS_FRACTION;
	      
	      if ( pd->v[var] )
		{
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      for (w = 0; w < pd->Num_Species_Eqn; w++ )
                        {
                          pvar = MAX_PROB_VAR + w;
			  
			  
                          
			  source = latent_heat * d_vnorm->C[w][j] * sign;
		      
		          source *= phi_i * lsi->delta;
		      
		          source *= det_J * wt * h3;
		          /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		          lec->J[peqn][pvar][ii][j] += source;
                        }
		    }
		} 
                
              /*
	       * J_W_v
	       */
              for( b=0; b<wim; b++)
                {
                  var = VELOCITY1 + b;
	          if ( pd->v[var] )
		    {
                      pvar = upd->vp[var];
                      for( j=0; j<ei->dof[var]; j++)
		        {
		          phi_j = bf[var]->phi[j];
                          
                          source = latent_heat * d_vnorm->v[b][j] * sign;
		      
		          source *= phi_i * lsi->delta;
		      
		          source *= det_J * wt * h3;
		          /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		          lec->J[peqn][pvar][ii][j] += source;
                        }
		    }
		} 
                
#ifdef COUPLED_FILL
              /*
	       * J_W_F
	       */
	      var = LS;
	      
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
                      
                      source = latent_heat * d_vnorm->F[j] * sign;
		      
		      source *= phi_i * lsi->delta;
                      
                      source += flux * phi_i * lsi->d_delta_dF[j];
		      
		      source *= det_J * wt * h3;
		      /*source *= pd->etm[eqn][(LOG2_SOURCE)];*/
		      
		      lec->J[peqn][pvar][ii][j] += source;
		    }
		}
#endif
	    }
	}
    }    
  return ( 1 );
}

int
assemble_max_strain ()   
{
  /*****************************************************************************
   * assemble_max_strain ()
   * 
   * Author:  Scott A Roberts, sarober@sandia.gov, 1514
   * 
   * Date:    April 12, 2012
   * 
   * Purpose: Calculates the maximum von Mises strain that the material has
   *          experienced over the length of the simulation.  This was primarily
   *          implemented for the FAUX_PLASTICITY model for the modulus in 
   *          solid mechanics.
   *****************************************************************************/

  // Define variables
  int eqn, peqn, var, pvar;
  int a, b, i, j, k;
  dbl phi_i, phi_j;                         
  dbl mass, source;
  dbl vmE, d_vmE_dstrain[DIM][DIM];

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_MAX_STRAIN;
  if (!pd->e[eqn]) return(status);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;                             // Gauss point weight
  dbl h3 = fv->h3;                             // Differential volume element
  dbl det_J = bf[eqn]->detJ;                   // Jacobian of transformation

  // FEM weights
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[eqn][(LOG2_MASS)];
  dbl etm_source = pd->etm[eqn][(LOG2_SOURCE)];

  // Calculate strain
  /* 
   * N.B.:  This strain is calculated in the same way the von Mises
   * stress is calculated, sqrt(3*II(T_dev)).  However, for this application,
   * we want to use a slightly different version, sqrt(4*II(E_dev)/3).
   * Therefore, we modify the returned values afterwards.
   */
  vmE = calc_tensor_invariant(fv->strain, d_vmE_dstrain, 4);
  vmE *= 2.0/3.0;
  for ( a = 0; a < DIM; a++) {
    for ( b = 0; b < DIM; b++) {
      d_vmE_dstrain[a][b] *= 2.0/3.0;
    }
  }

  // Prepare tests
  int use_new = vmE >= fv_old->max_strain;
  int mass_lump = 0;

  /* --- Assemble residuals --------------------------------------------------*/
  eqn = R_MAX_STRAIN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble mass term
      mass = 0.0;
      if ( T_MASS ) {
	if ( mass_lump ) {
	  mass -= *esp->max_strain[i] * phi_i;
	} else {
	  mass -= fv->max_strain * phi_i;
	}
      }
      mass *= dA * etm_mass;

      // Assemble source term
      source = 0.0;
      if ( T_SOURCE ) {
	if ( use_new ) {
	  source += vmE * phi_i;
	} else {
	  source += fv_old->max_strain * phi_i;
	}
      }
      source *= dA * etm_source;
      
      // Assemble full residual
      lec->R[peqn][i] += mass + source;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_MAX_STRAIN
    
  
  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_MAX_STRAIN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble sensitivities for MAX_STRAIN
      var = MAX_STRAIN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  if ( T_MASS ) { 
	    if ( mass_lump ) {
	      mass -= phi_i * delta(i,j);
	    } else {
	      mass -= phi_i * phi_j;
	    }
	  }
	  mass *= dA * etm_mass;
	  
	  // Assemble source term
	  source = 0.0;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + source;
	  
	} // End of loop over DOF (j)
	
      } // End of MAX_STRAIN sensitivities
      
      // Assemble sensitivities for MESH_DISPLACEMENT
      for ( k = 0; k < DIM; k++) {
	var = MESH_DISPLACEMENT1 + k;
	if (pd->v[var]) {
	  pvar = upd->vp[var];
	  
	  // Loop over DOF (j)
	  for ( j = 0; j < ei->dof[var]; j++) {
	    
	    // Load basis functions
	    phi_j = bf[eqn]->phi[j];
	    
	    // Assemble mass terms
	    mass = 0.0;
	    if ( T_MASS ) { 
	      if ( mass_lump ) {
		mass -= *esp->max_strain[i] * phi_i * fv->dh3dmesh[k][j] * det_J;
		mass -= *esp->max_strain[i] * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
	      } else {
		mass -= fv->max_strain * phi_i * fv->dh3dmesh[k][j] * det_J;
		mass -= fv->max_strain * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
	      }
	    }
	    mass *= wt * etm_mass;
	    
	    // Assemble source term
	    source = 0.0;
	    if ( T_SOURCE ) { 
	      if ( use_new ) {
		for ( a = 0; a < DIM; a++ ) {
		  for ( b = 0; b < DIM; b++ ) {
		    source += d_vmE_dstrain[a][b] * fv->d_strain_dx[a][b][k][j] * h3 * det_J;
		  }
		}
		source += vmE * fv->dh3dmesh[k][j] * det_J;
		source += vmE * h3 * bf[eqn]->d_det_J_dm[k][j];
	      } else {
		source += fv_old->max_strain * fv->dh3dmesh[k][j] * det_J;
		source += fv_old->max_strain * h3 * bf[eqn]->d_det_J_dm[k][j];
	      }
	    }
	    source *= phi_i * wt * etm_source;
	    
	    // Assemble full Jacobian
	    lec->J[peqn][pvar][i][j] += mass + source;
	    
	  } // End of loop over DOF (j)
	  
	} // End of MESH_DISPLACEMENTi sensitivities
	
      } // End loop over MESH_DISPLACEMENT vector components
      
    } // End of loop over DOF (i)
    
  } // End of Jacobian assembly
  
  return(status);
} // End of assemble_max_strain()

int
assemble_cur_strain ()   
{
  /*****************************************************************************
   * assemble_cur_strain ()
   * 
   * Author:  Scott A Roberts, sarober@sandia.gov, 1514
   * 
   * Date:    April 23, 2012
   * 
   * Purpose: Calculates the current von Mises strain.  This was primarily
   *          implemented for the FAUX_PLASTICITY model for the modulus in 
   *          solid mechanics.
   *****************************************************************************/

  // Define variables
  int eqn, peqn, var, pvar;
  int a, b, i, j, k;
  dbl phi_i, phi_j;                         
  dbl mass, source;
  dbl vmE, d_vmE_dstrain[DIM][DIM];

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_CUR_STRAIN;
  if (!pd->e[eqn]) return(status);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;                             // Gauss point weight
  dbl h3 = fv->h3;                             // Differential volume element
  dbl det_J = bf[eqn]->detJ;                   // Jacobian of transformation

  // FEM weights
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[eqn][(LOG2_MASS)];
  dbl etm_source = pd->etm[eqn][(LOG2_SOURCE)];

  // Calculate strain
  /* 
   * N.B.:  This strain is calculated in the same way the von Mises
   * stress is calculated, sqrt(3*II(T_dev)).  However, for this application,
   * we want to use a slightly different version, sqrt(4*II(E_dev)/3).
   * Therefore, we modify the returned values afterwards.
   */
  vmE = calc_tensor_invariant(fv->strain, d_vmE_dstrain, 4);
  vmE *= 2.0/3.0;
  for ( a = 0; a < DIM; a++) {
    for ( b = 0; b < DIM; b++) {
      d_vmE_dstrain[a][b] *= 2.0/3.0;
    }
  }

  // Prepare tests
  int mass_lump = 0;

  /* --- Assemble residuals --------------------------------------------------*/
  eqn = R_CUR_STRAIN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble mass term
      mass = 0.0;
      if ( T_MASS ) {
        if ( mass_lump ) {
	  mass -= *esp->cur_strain[i] * phi_i;
	} else {
	  mass -= fv->cur_strain * phi_i;
	}
      }
      mass *= dA * etm_mass;

      // Assemble source term
      source = 0.0;
      if ( T_SOURCE ) {
        source += vmE * phi_i;
      }
      source *= dA * etm_source;
      
      // Assemble full residual
      lec->R[peqn][i] += mass + source;
      
    }  // End of loop over DOF (i)

  } // End of residual assembly of R_CUR_STRAIN
    
  
  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_CUR_STRAIN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[eqn];
    
    // Loop over DOF (i)
    for ( i = 0; i < ei->dof[eqn]; i++) {
      
      // Load basis functions
      phi_i = bf[eqn]->phi[i];
      
      // Assemble sensitivities for CUR_STRAIN
      var = CUR_STRAIN;
      if (pd->v[var]) {
  	pvar = upd->vp[var];

  	// Loop over DOF (j)
  	for ( j = 0; j < ei->dof[var]; j++) {

	  // Load basis functions
	  phi_j = bf[eqn]->phi[j];
	  
	  // Assemble mass terms
	  mass = 0.0;
	  if ( T_MASS ) { 
	    if ( mass_lump ) {
	      mass -= phi_i * delta(i,j);
	    } else {
	      mass -= phi_i * phi_j;
	    }
	  }
	  mass *= dA * etm_mass;
	  
	  // Assemble source term
	  source = 0.0;
	  
	  // Assemble full Jacobian
	  lec->J[peqn][pvar][i][j] += mass + source;
	  
	} // End of loop over DOF (j)
	
      } // End of CUR_STRAIN sensitivities
      
      // Assemble sensitivities for MESH_DISPLACEMENT
      for ( k = 0; k < DIM; k++) {
	var = MESH_DISPLACEMENT1 + k;
	if (pd->v[var]) {
	  pvar = upd->vp[var];
	  
	  // Loop over DOF (j)
	  for ( j = 0; j < ei->dof[var]; j++) {
	    
	    // Load basis functions
	    phi_j = bf[eqn]->phi[j];
	    
	    // Assemble mass terms
	    mass = 0.0;
	    if ( T_MASS ) { 
	      if ( mass_lump ) {
		mass -= *esp->cur_strain[i] * phi_i * fv->dh3dmesh[k][j] * det_J;
		mass -= *esp->cur_strain[i] * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
	      } else {
		mass -= fv->cur_strain * phi_i * fv->dh3dmesh[k][j] * det_J;
		mass -= fv->cur_strain * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
	      }
	    }
	    mass *= wt * etm_mass;
	    
	    // Assemble source term
	    source = 0.0;
	    if ( T_SOURCE ) { 
              for ( a = 0; a < DIM; a++ ) {
                for ( b = 0; b < DIM; b++ ) {
                  source += d_vmE_dstrain[a][b] * fv->d_strain_dx[a][b][k][j] * h3 * det_J;
                }
              }
              source += vmE * fv->dh3dmesh[k][j] * det_J;
              source += vmE * h3 * bf[eqn]->d_det_J_dm[k][j];
	    }
	    source *= phi_i * wt * etm_source;
	    
	    // Assemble full Jacobian
	    lec->J[peqn][pvar][i][j] += mass + source;
	    
	  } // End of loop over DOF (j)
	  
	} // End of MESH_DISPLACEMENTi sensitivities
	
      } // End loop over MESH_DISPLACEMENT vector components
      
    } // End of loop over DOF (i)
    
  } // End of Jacobian assembly
  
  return(status);
} // End of assemble_cur_strain()


/*  _______________________________________________________________________  */

/* assemble_acoustic -- assemble terms (Residual &| Jacobian) for acoustic harmonic
 *				wave equations
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
 * Created:	Thu June 28, 2005 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_acoustic(double time,	/* present time value */
		  double tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
		  double dt,	/* current time step size */
		  const PG_DATA *pg_data,
		  const int ac_eqn,	/* acoustic eqn id and var id	*/
		  const int ac_var )
{
  int eqn, var, peqn, pvar, dim, p, b, w, i, j, status;
  int conj_var = 0;			/* identity of conjugate variable  */

  dbl P = 0, P_conj = 0, sign_conj = 0;		/* acoustic pressure	*/
  dbl q[DIM];				
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT d_q_struct; 
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  dbl R;				/* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl alpha;				/* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  dbl mass;		         	/* For terms and their derivatives */

  dbl advection;			/* For terms and their derivatives */

  dbl diffusion;
  dbl diff_a, diff_b, diff_c, diff_d;
  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM];		/* ditto.  */
  dbl wt;

  /* initialize grad_phi_i */
  for (i = 0; i < DIM; i++) {
    grad_phi_i[i] = 0;
  }

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim   = pd->Num_Dim;
  eqn   = ac_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary 
   * with temperature, spatial coordinates, and species concentration.
   */

  R = acoustic_impedance( d_R, time );

  k = wave_number( d_k, time );

  alpha = acoustic_absorption( d_alpha, time );

  acoustic_flux( q, d_q, time, ac_eqn, ac_var );

  if( ac_eqn == R_ACOUS_PREAL )
	{ P = fv->apr; P_conj = fv->api; conj_var = ACOUS_PIMAG; sign_conj=1.;}
  else if (eqn == R_ACOUS_PIMAG)
	{ P = fv->api; P_conj = -fv->apr; conj_var = ACOUS_PREAL; sign_conj=-1.;}
  else
	{ EH(-1,"Invalid Acoustic eqn");}
  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      eqn = ac_eqn;
      peqn = upd->ep[eqn];
      var = ac_var;
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  mass = 0.;

	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {
              advection += phi_i * (k/R)*2*alpha*P_conj * det_J*wt;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }

	  diffusion = 0.;
	  if ( pd->e[eqn] & T_DIFFUSION )
	    {
	      for ( p=0; p<VIM; p++)
		{
		  grad_phi_i[p] = bf[eqn]->grad_phi[i] [p];
		}

	      for ( p=0; p<VIM; p++)
		{
		  diffusion += grad_phi_i[p] * q[p];
		}
	      diffusion *= det_J * wt;
	      diffusion *= h3;
	      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  source = 0.;
	  if ( pd->e[eqn] & T_SOURCE )
	    {
	      source -= phi_i * (k/R)*P * det_J * wt;
	      source *= h3;
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	    }

	  lec->R[peqn][i] +=
	    mass + advection +  diffusion + source;


	}
    }


  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      eqn   = ac_eqn;
      peqn = upd->ep[eqn];
      for ( i=0; i<ei->dof[eqn]; i++)
	{
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  for ( p=0; p<VIM; p++)
	    {
	      grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
	    }

	  /*
	   * J_e_ap
	   */
	  var = ac_var;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  mass = 0.;

		  advection = 0.;

		  diffusion = 0.;
		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += d_q->P[p][j] * grad_phi_i[p];
			}
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source -= phi_i * (k/R)*phi_j*det_J*wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }
	/*
	 *  Conjugate pressure variable
	 */
	  var = conj_var;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
                      advection += phi_i * (k/R)*2*alpha*sign_conj*phi_j*det_J*wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }
		  lec->J[peqn][pvar][i][j] += advection;
		}
	    }
	  /*
	   * J_e_T
	   */
	  var = TEMPERATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  mass = 0.;

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
                      advection += phi_i * 2*P_conj*(R*(alpha*d_k->T[j]+k*d_alpha->T[j])
				-k*alpha*d_R->T[j])/(R*R)* det_J*wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  diffusion = 0.;
		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += d_q->T[p][j] * grad_phi_i[p];
			}
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
		      source -= phi_i * (R*d_k->T[j]-k*d_R->T[j])/(R*R)*P* det_J*wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		}
	    }

	  /*
	   * J_e_d
	   */
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

		      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

		      mass = 0.;

		      advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
                      advection += phi_i*2*P_conj*((k*alpha/R) * (d_det_J_dmeshbj*h3 +
						det_J*dh3dmesh_bj) +
					(R*(alpha*d_k->X[b][j]+k*d_alpha->X[b][j])
				-k*alpha*d_R->X[b][j])/(R*R)* det_J*h3) *wt;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

			  /*
			   * multiple parts:
			   * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
			   *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
			   *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
			   *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
			   */
		      diffusion = 0.;
                      if ( pd->e[eqn] & T_DIFFUSION )
		        {
                          diff_a = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      dgrad_phi_i_dmesh[p]
				= bf[eqn]->d_grad_phi_dmesh[i][p] [b][j];

			      diff_a += dgrad_phi_i_dmesh[p] * q[p];
			    }
			  diff_a *= det_J * h3 * wt;

                          diff_b = 0.;
                          for ( p=0; p<VIM; p++)
			    {
			      diff_b += d_q->X[p][b][j] * grad_phi_i[p];
			    }
		          diff_b *= det_J * h3 * wt;

			  diff_c = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      diff_c += grad_phi_i[p] * q[p];
			    }
			  diff_c *= d_det_J_dmeshbj * h3 * wt;

			  diff_d = 0.;
			  for ( p=0; p<dim; p++)
			    {
			      diff_d += grad_phi_i[p] * q[p];
			    }
			  diff_d *= det_J * dh3dmesh_bj * wt;

			  diffusion = diff_a + diff_b + diff_c + diff_d;

			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}


		      source = 0.;

		      if ( pd->e[eqn] & T_SOURCE )
			{
	      		  source -= phi_i * P*((k/R) * (d_det_J_dmeshbj*h3 +
						det_J*dh3dmesh_bj) +
					(R*d_k->X[b][j]-k*d_R->X[b][j])/(R*R)* 
					det_J *h3) * wt;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += mass + advection + diffusion + source;
		    }
		}
	    }

	  /*
	   * J_e_c
	   */
	  var = MASS_FRACTION;
	  if ( pd->e[eqn] && pd->v[var] )
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      mass = 0.;

		      advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
                      advection += phi_i*2*P_conj*(R*(alpha*d_k->C[w][j]+k*d_alpha->C[w][j])
				-k*alpha*d_R->C[w][j])/(R*R)* det_J*wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		      diffusion = 0.;
		      if ( pd->e[eqn] & T_DIFFUSION )
			{
			  for ( p=0; p<dim; p++)
			    {
			      diffusion += grad_phi_i[p] * d_q->C[p][w][j];
			    }
			  diffusion *= det_J * wt;
			  diffusion *= h3;
			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}

		      source = 0.;
		      if ( pd->e[eqn] & T_SOURCE )
			{
		          source -= phi_i * (R*d_k->C[w][j]-k*d_R->C[w][j])/(R*R)
					*P* det_J*wt;
			  source *= h3;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][MAX_PROB_VAR + w][i][j] += advection + mass + diffusion + source;
		    }
		}
	    }

	}
    }

  return(status);
} /* end of assemble_acoustic */

/* assemble_acoustic_reynolds_stress -- assemble terms for acoustic reynolds stress
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
 * Created:	Thu Aug. 18, 2005 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_acoustic_reynolds_stress(double time,	/* present time value */
		  double tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
		  double dt,	/* current time step size */
		  const PG_DATA *pg_data)
{
  int eqn, var, peqn, pvar, dim, p, b, w, i, j, status;

  dbl R;				/* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl omega, acous_pgrad=0;

  dbl mass;		         	/* For terms and their derivatives */

  dbl advection;			/* For terms and their derivatives */

  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM];

  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl wt;

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim   = pd->Num_Dim;
  eqn   = R_ACOUS_REYN_STRESS;
  /*
   * Bail out fast if there's nothing to do...
   */
  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary 
   * with temperature, spatial coordinates, and species concentration.
   */

  R = acoustic_impedance( d_R, time );

  k = wave_number( d_k, time );
  omega = upd->Acoustic_Frequency;

  for (p = 0; p < dim; p++)
      {
       acous_pgrad += fv->grad_api[p]*fv->grad_api[p];
       acous_pgrad += fv->grad_apr[p]*fv->grad_apr[p];
      }

  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[eqn];
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  mass = 0.;
	  if ( pd->e[eqn] & T_MASS )
	    {
	      mass += phi_i * fv->ars * det_J*wt;
	      mass *= h3;
	      mass *= pd->etm[eqn][(LOG2_MASS)];
	    }

	  advection = 0.;
	  if ( pd->e[eqn] & T_ADVECTION )
	    {
 	      advection -= phi_i * (1./(4.*k*R*omega))*(acous_pgrad) * det_J*wt;
	      advection *= h3;
	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
	    }

	  source = 0.;
	  if ( pd->e[eqn] & T_SOURCE )
	    {
 	      source -= phi_i * k/(4.*R*omega)*
			(fv->apr*fv->apr+fv->api*fv->api) * det_J * wt;
	      source *= h3;
	      source *= pd->etm[eqn][(LOG2_SOURCE)];
	    }

	  lec->R[peqn][i] += mass + advection + source;


	}
    }


  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      peqn = upd->ep[eqn];
      for ( i=0; i<ei->dof[eqn]; i++)
	{
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  /*
	   * J_e_ars
	   */
	  var = ACOUS_REYN_STRESS;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  mass = 0.;
		  if ( pd->e[eqn] & T_MASS )
		    {
	      		mass += phi_i * phi_j * det_J*wt;
	      		mass *= h3;
	      		mass *= pd->etm[eqn][(LOG2_MASS)];
		    }

		  advection = 0.;

		  source = 0.;

		  lec->J[peqn][pvar][i][j] += mass + advection + source;
		}
	    }
	  /*
	   * J_e_ap
	   */
	  var = ACOUS_PREAL;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  mass = 0.;

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  advection += 2.*fv->grad_apr[p]*grad_phi_j[p];
			}
 	      	      advection *= -phi_i * (1./(4.*k*R*omega))* det_J*wt;
	      	      advection *= h3;
	      	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
 	      		source -= phi_i * k/(4.*R*omega)*(2.*fv->apr*phi_j)
					 * det_J * wt;
	      		source *= h3;
		        source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + source;
		}
	    }
	/*
	 *  Conjugate pressure variable
	 */
	  var = ACOUS_PIMAG;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  advection += 2.*fv->grad_api[p]*grad_phi_j[p];
			}
 	      	      advection *= -phi_i * (1./(4.*k*R*omega))* det_J*wt;
	      	      advection *= h3;
	      	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
 	      		source -= phi_i * k/(4.*R*omega)*(2.*fv->api*phi_j)
					 * det_J * wt;
	      		source *= h3;
		        source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += advection + source;
		}
	    }
	  /*
	   * J_e_T
	   */
	  var = TEMPERATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  mass = 0.;

		  advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
 	      	      advection -= phi_i/(4.*omega)*
				(-1./(k*k*R*R))*(k*d_R->T[j]+R*d_k->T[j])
					*(acous_pgrad) * det_J*wt;
		      advection *= h3;
		      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		  source = 0.;
		  if ( pd->e[eqn] & T_SOURCE )
		    {
 	      	      source -= phi_i/(4.*omega)*
				(R*d_k->T[j]-k*d_R->T[j])/(R*R) *
			(fv->apr*fv->apr+fv->api*fv->api) * det_J * wt;
		      source *= h3;
		      source *= pd->etm[eqn][(LOG2_SOURCE)];
		    }

		  lec->J[peqn][pvar][i][j] += mass + advection + source;
		}
	    }

	  /*
	   * J_e_d
	   */
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

		      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

		      mass = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
	      	     mass += phi_i * fv->ars * 
			(d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj)*wt;
	      	     mass *= pd->etm[eqn][(LOG2_MASS)];
		    }
		      advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  advection += 2.*(fv->grad_apr[p]*fv->d_grad_apr_dmesh[p][b][j] + fv->grad_api[p]*fv->d_grad_api_dmesh[p][b][j]);
			}
 	      	      advection *= -phi_i * (1./(4.*k*R*omega))* det_J*h3*wt;
 	      	      advection -= phi_i * (1./(4.*k*R*omega))*(acous_pgrad) * 
				(d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj)*wt;
 	      	      advection -= phi_i/(4.*omega)*
				(-1./(k*R*k*R))*(k*d_R->X[b][j] + R*d_k->X[b][j])
				 *(acous_pgrad) * det_J * h3*wt;
	      	      advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		      source = 0.;

		      if ( pd->e[eqn] & T_SOURCE )
			{
 	      		  source -= phi_i/(4.*omega)*((k/R)*
					(fv->apr*fv->apr+fv->api*fv->api) * 
				(d_det_J_dmeshbj*h3 + det_J*dh3dmesh_bj)
				+ det_J*h3*(R*d_k->X[b][j]-k*d_R->X[b][j])/(R*R))
				 * wt;
			  source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][pvar][i][j] += mass + advection + source;
		    }
		}
	    }

	  /*
	   * J_e_c
	   */
	  var = MASS_FRACTION;
	  if ( pd->e[eqn] && pd->v[var] )
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      mass = 0.;

		      advection = 0.;
		  if ( pd->e[eqn] & T_ADVECTION )
		    {
 	      	     advection -= phi_i/(4.*omega)*(-1./(k*k*R*R))*
			(k*d_R->C[w][j] + R*d_k->C[w][j])*acous_pgrad*det_J*wt;
		     advection *= h3;
		     advection *= pd->etm[eqn][(LOG2_ADVECTION)];
		    }

		      source = 0.;
		      if ( pd->e[eqn] & T_SOURCE )
			{
 	      		 source -= phi_i/(4.*omega)*
				(R*d_k->C[w][j]-k*d_R->C[w][j])/(R*R)* 
				(fv->apr*fv->apr+fv->api*fv->api) * det_J * wt;
			 source *= h3;
			 source *= pd->etm[eqn][(LOG2_SOURCE)];
			}

		      lec->J[peqn][MAX_PROB_VAR + w][i][j] += advection + mass + source;
		    }
		}
	    }

	}
    }

  return(status);
} /* end of assemble_acoustic_reynolds_stress */

void
acoustic_flux( double q[DIM],
           ACOUSTIC_FLUX_DEPENDENCE_STRUCT *d_q,
           double time,
	   const int eqn,
	   const int ac_var )
{
  dbl R;				/* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; /* Acoustic impedance dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl kR_inv;				/* k*R product inverse	*/

  dbl grad_P[DIM];			/* Acoustic pressure. */

  int b, j, p, w;
  int var;

  int dim = pd->Num_Dim;
  int wim;

  wim = dim;
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN) wim = wim+1;

  if ( d_q == NULL ) {d_R = NULL;   d_k = NULL;}
  
  R = acoustic_impedance( d_R, time );
  k = wave_number( d_k, time );
  kR_inv = 1./(k*R);


  if( eqn == R_ACOUS_PREAL )
	{
  	for ( p=0; p<VIM; p++)
    		{
      		grad_P[p] = fv->grad_apr[p];
    		}
	}
  else if (eqn == R_ACOUS_PIMAG)
	{
  	for ( p=0; p<VIM; p++)
    		{
      		grad_P[p] = fv->grad_api[p];
    		}
	}
  else if (eqn == R_LIGHT_INTP)
	{
  	for ( p=0; p<VIM; p++)
    		{
      		grad_P[p] = fv->grad_poynt[0][p];
    		}
	}
  else if (eqn == R_LIGHT_INTM)
	{
  	for ( p=0; p<VIM; p++)
    		{
      		grad_P[p] = fv->grad_poynt[1][p];
    		}
	}
  else if (eqn == R_LIGHT_INTD)
	{
  	for ( p=0; p<VIM; p++)
    		{
      		grad_P[p] = fv->grad_poynt[2][p];
    		}
	}
  else
	{ EH(-1,"Invalid Acoustic eqn");}

  for ( p=0; p<VIM; p++)
	{
	  q[p] =  kR_inv * grad_P[p];
	}

  var = ac_var;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->P[p][j] = kR_inv * bf[var]->grad_phi[j][p];
                }
            }
        }

  var = TEMPERATURE;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->T[p][j] = -kR_inv*kR_inv*(k*d_R->T[j]+R*d_k->T[j]) * grad_P[p];
                }
            }
        }

 var = MASS_FRACTION;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( w=0; w<pd->Num_Species_Eqn; w++)
            {
              for ( p=0; p<VIM; p++)
                {
                  for ( j=0; j<ei->dof[var]; j++)
                    {
                      d_q->C[p][w][j] = -kR_inv*kR_inv*(k*d_R->C[w][j]+R*d_k->C[w][j]) * grad_P[p];
                    }
                }
            }
        }

 var = FILL;
      if ( d_q != NULL && pd->v[var] )
        {
          for ( p=0; p<VIM; p++)
            {
              for ( j=0; j<ei->dof[var]; j++)
                {
                  d_q->F[p][j] = -kR_inv*kR_inv*(k*d_R->F[j]+R*d_k->F[j]) * grad_P[p];
                }
            }
        }

 if ( d_q != NULL && pd->v[MESH_DISPLACEMENT1] )
        {
          for ( p=0; p<VIM; p++)
	    {
	      for ( b=0; b<wim; b++)
                {
                  var = MESH_DISPLACEMENT1+b;
                  for ( j=0; j<ei->dof[var]; j++)
                    {
			if(eqn == R_ACOUS_PREAL)
			{
                      	d_q->X[p][b][j] = kR_inv * fv->d_grad_apr_dmesh[p][b][j]
                                 -kR_inv*kR_inv*(k*d_R->X[b][j]+R*d_k->X[b][j])*grad_P[p];
			}
			else
			{
                      	d_q->X[p][b][j] = kR_inv * fv->d_grad_api_dmesh[p][b][j]
                                 -kR_inv*kR_inv*(k*d_R->X[b][j]+R*d_k->X[b][j])*grad_P[p];
			}
                    }
                }
            }
        }
return;
}

int
assemble_ars_source ( double ars_jump, double grad_jump )
{
  int i,j, a, p, ii,ledof;
  int eqn, peqn, var, pvar;
  int dim, wim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double omega, force, temp, acous_pgrad=0, force1, temp1;

  double R;				/* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  double k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  double time = tran->time_value;

  eqn = R_MOMENTUM1;
  if ( ! pd->e[eqn] )
    {
      return(0);
    }

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  dim   = pd->Num_Dim;
  wim   = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    wim = wim+1;
/* local force caused by Acoustic Reynolds Stress jump at interface	*/

  omega = upd->Acoustic_Frequency;
  temp = ars_jump/(4.*omega);
  force = (SQUARE(fv->apr)+SQUARE(fv->api))*temp;

  R = acoustic_impedance( d_R, time );

  k = wave_number( d_k, time );

  for (p = 0; p < dim; p++)
      {
       acous_pgrad += fv->grad_api[p]*fv->grad_api[p];
       acous_pgrad += fv->grad_apr[p]*fv->grad_apr[p];
      }
  temp1 = grad_jump/(4.*omega)/SQUARE(k*R);
  force1 = acous_pgrad*temp1;
  force += force1;
  
#ifdef COUPLED_FILL
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if ( ls->CalcSurfDependencies )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = phi_i * lsi->normal[a] * force * lsi->delta;

		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  /* J_m_F
		   */
                  var = LS ;
		  pvar = upd->vp[var];
		      
		  for( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];
		      lec->J[peqn][pvar][ii][j] += source * phi_j;
		    }
                }
            }
        }
      return(0);
    }
#endif

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

	  for (i = 0; i < ei->dof[eqn]; i++)
	    {

	      ledof = ei->lvdof_to_ledof[eqn][i];

	      if (ei->active_interp_ledof[ledof])
		{

		  ii = ei->lvdof_to_row_lvdof[eqn][i];

		  phi_i = bfm->phi[i];

		  source = phi_i * lsi->normal[a] * force * lsi->delta;

		  source *= det_J * wt * h3;
		  source *= pd->etm[eqn][(LOG2_SOURCE)];

		  lec->R[peqn][ii] += source;
		}
	    }
	}
    }

  if ( af->Assemble_Jacobian )
    {
      for( a=0; a<wim; a++ )
	{
	  eqn = R_MOMENTUM1 + a;
	  peqn = upd->ep[eqn];
	  bfm = bf[eqn];

          for (i = 0; i < ei->dof[eqn]; i++)
            {
              ledof = ei->lvdof_to_ledof[eqn][i];

              if (ei->active_interp_ledof[ledof])
                {
                  ii = ei->lvdof_to_row_lvdof[eqn][i];

                  phi_i = bfm->phi[i];

#ifdef COUPLED_FILL
                   /*
                    * J_m_F
                    */
                  var = FILL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * force *
                                   ( lsi->delta * lsi->d_normal_dF[a][j] +
                                     lsi->d_delta_dF[j] * lsi->normal[a] );

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
#endif
                   /*
                    * Acoustic Pressure - Real
                    */
                  var = ACOUS_PREAL;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];
                          source = 2*fv->apr*phi_j*temp;
		      	  for ( p=0; p<dim; p++)
			    {
		     source += 2.*fv->grad_apr[p]*bf[var]->grad_phi[j][p]*temp1;
			    }

                          source *= phi_i * lsi->delta * lsi->normal[a];

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                   /*
                    * Acoustic Pressure - Imag
                    */
                  var = ACOUS_PIMAG;

                  if ( pd->v[var] )
                    {
                      pvar = upd->vp[var];

                      for( j=0; j<ei->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];
                          source = 2*fv->api*phi_j*temp;
		      	  for ( p=0; p<dim; p++)
			    {
		     source += 2.*fv->grad_api[p]*bf[var]->grad_phi[j][p]*temp1;
			    }

                          source *= phi_i * lsi->delta * lsi->normal[a];

                          source *= det_J * wt * h3;
                          source *= pd->etm[eqn][LOG2_SOURCE];

                          lec->J[peqn][pvar][ii][j] += source;
                        }
                    }
                }
            }
        }
    }
  return ( 1 );
}
/*
 * Viscous Dissipation Heating Source due to Acoustic Waves Model 
 */

/*
 * double visc_diss_acoustic_source (dh, param)
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for filling up the following forces and sensitivities
 * at the current gauss point:
 *     intput:  
 *
 *     output:  h             - heat source
 *              d_h->T[j]    - derivative wrt temperature at node j.
 *              d_h->V[j]    - derivative wrt voltage at node j.
 *              d_h->C[i][j] - derivative wrt mass frac species i at node j
 *              d_h->v[i][j] - derivative wrt velocity component i at node j
 *              d_h->X[0][j] - derivative wrt mesh displacement components i at node j
 *              d_h->S[m][i][j][k] - derivative wrt stress mode m of component ij at node k
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
 * ---------------------------------------------------------------------------
 */

double
visc_diss_acoustic_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
		      dbl *param, int num_const) /* General multipliers */
{
  /* Local Variables */
  int var;

  int w, j;
  double omega, visc_first, visc_second, R_gas;
  double h, temp1, temp3, ap_square;
  double R;				/* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  double k;				/* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  double alpha;				/* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  double visc_cmb;				/* Combined viscosity term  */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_visc_cmb_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_visc_cmb = &d_visc_cmb_struct;

  dbl gamma[DIM][DIM];
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  double time = tran->time_value;

  /* Begin Execution */

  /**********************************************************/
  /* Source constant * viscosity* gammadot .. grad_v */

  /* get mu and grad_mu */
  omega = upd->Acoustic_Frequency;
  R = acoustic_impedance( d_R, time );
  k = wave_number( d_k, time );
  alpha = acoustic_absorption( d_alpha, time );

/* modulate viscosity term according to the LS field 
	(i.e. 4/3 shear viscosity plus bulk viscosity)
*/
  visc_first = param[1];
  visc_second = param[2];
  memset( gamma, 0, sizeof(dbl)*DIM*DIM);
/*  There are some options for the appropriate shear rate which
	to evaluate the viscosity - could choose zero-shear-rate,
	oscillation frequency or rate based on grad(v).
    some of the options require extra derivatives be calculated.
	I'll use oscillation frequency for now.
*/
#if 1
  gamma[0][0] = omega;
#endif
#if 0
  gamma[0][0] = k/R*(SQUARE(fv->apr)+SQUARE(fv->api))*
		sqrt(2.*(1+4.*SQUARE(alpha)/SQUARE(omega)));
#endif
  mu = viscosity(gn, gamma, d_mu);
  visc_cmb = 4.*mu/3. + visc_first;
      memset( d_visc_cmb->T, 0, sizeof(double)*MDE);
      memset( d_visc_cmb->X, 0, sizeof(double)*DIM*MDE);
      memset( d_visc_cmb->C, 0, sizeof(double)*MAX_CONC*MDE);
      memset( d_visc_cmb->F, 0, sizeof(double)*MDE);
  for (j=0; j<ei->dof[TEMPERATURE]; j++)
	    { d_visc_cmb->T[j] = d_mu->T[j]; }
  for (j=0; j<ei->dof[FILL]; j++)
	    { d_visc_cmb->F[j] = d_mu->F[j]; }
  for (j=0; j<ei->dof[MASS_FRACTION]; j++)
      {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
	    { d_visc_cmb->C[w][j] = d_mu->C[w][j]; }
      }
  if( ls != NULL)
  {
  visc_cmb = ls_modulate_thermalconductivity( visc_cmb,
                                  visc_second,
                                  ls->Length_Scale,
                                  (double) mp->mp2nd->wavenumbermask[0],
                                  (double) mp->mp2nd->wavenumbermask[1],
                                  d_visc_cmb );
/*  optionally set impedance to true value input on card	*/
  if(num_const == 4)
	{
	R_gas = param[3];
  	R = ls_modulate_thermalconductivity( mp->acoustic_impedance,
                                  R_gas,
                                  ls->Length_Scale,
                                  (double) mp->mp2nd->acousticimpedancemask[0],
                                  (double) mp->mp2nd->acousticimpedancemask[1],
                                  d_R );
	}
  }

#if 0
fprintf(stderr,"visc %g %g %g %d %g %g %g %d\n",mu,visc_cmb,gamma[0][0], num_const, param[3], fv->F, lsi->H, ls->SubElemIntegration);
#endif
  ap_square = SQUARE(fv->apr) + SQUARE(fv->api);
  temp1 = SQUARE(k)/SQUARE(R);
  temp3 = (1. + 4.*SQUARE(alpha)/SQUARE(omega));
  h = ap_square*temp1*visc_cmb*temp3;
  h *= 0.5*param[0];

  /* Now do sensitivies */
if(af->Assemble_Jacobian)
  {

  var = TEMPERATURE;
  if ( d_h != NULL && pd->v[var] )
    {
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      d_h->T[j] += 0.5*param[0]*ap_square*visc_cmb*
			(temp1*4.*alpha*d_alpha->T[j]/SQUARE(omega) 
			+ temp3*2.*k*(R*d_k->T[j]-k*d_R->T[j])/(R*R*R));
	    }
    }


  var = MASS_FRACTION;
  if ( d_h != NULL && pd->v[var] )
    {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
         {
	  for (j=0; j<ei->dof[var]; j++)
	    {
	     d_h->C[w][j] += 0.5*param[0]*ap_square*visc_cmb*
			(temp1*4.*alpha*d_alpha->C[w][j]/SQUARE(omega) 
			+ temp3*2.*k*(R*d_k->C[w][j]-k*d_R->C[w][j])/(R*R*R));
	    }
         }
    }

  var = FILL;
  if ( d_h != NULL && pd->v[var] )
    {
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      d_h->F[j] += 0.5*param[0]*ap_square*
			(temp1*visc_cmb*4.*alpha*d_alpha->F[j]/SQUARE(omega) 
			+ visc_cmb*temp3*2.*k*(R*d_k->F[j]-k*d_R->F[j])/(R*R*R)
			+ temp1*temp3*d_visc_cmb->F[j]);
	    }
    }

  var = ACOUS_PREAL;
  if ( d_h != NULL && pd->v[var] )
    {
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      d_h->APR[j] += 0.5*param[0]*temp1*visc_cmb*temp3*2.*fv->apr*bf[var]->phi[j]; 
	    }
    }
  var = ACOUS_PIMAG;
  if ( d_h != NULL && pd->v[var] )
    {
	  for (j=0; j<ei->dof[var]; j++)
	    {
	      d_h->API[j] += 0.5*param[0]*temp1*visc_cmb*temp3*2.*fv->api*bf[var]->phi[j]; 
	    }
    }
   
 
  }  /* end of if Assemble Jacobian  */

  return(h);
}
/*  _______________________________________________________________________  */

/* assemble_poynting -- assemble terms (Residual &| Jacobian) for Beer's law 
 *            type light intensity absorption equations
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
 * Created:	Thu February 16, 2011 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int
assemble_poynting(double time,	/* present time value */
		  double tt,	/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
		  double dt,	/* current time step size */
		  const PG_DATA *pg_data,
		  const int py_eqn,	/* eqn id and var id	*/
		  const int py_var )
{
  int eqn, var, peqn, pvar, dim, p, b, w, i, j, status, light_eqn = 0;

  dbl P;		                /* Light Intensity	*/
  dbl grad_P, Psign = 0;		/* grad intensity */

  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl alpha;				/* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  dbl diffusion = 0;
  //  dbl diff_a;
  dbl diff_b, diff_c, diff_d;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM];

  dbl h3;			/* Volume element (scale factors). */
  dbl dh3dmesh_bj;		/* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;			/* for specified (b,j) mesh dof */
  dbl wt;

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim   = pd->Num_Dim;
  eqn   = py_eqn;
  /*
 *    Radiative transfer equation variables - connect to input file someday
 */
  double svect[3]={0.,-1.,0.};
  double mucos=1.0;

  /*
   * Bail out fast if there's nothing to do...
   */
  if ( ! pd->e[eqn] )
    {
      return(status);
    }

  wt = fv->wt;				/* Gauss point weight. */
  h3 = fv->h3;			/* Differential volume element. */
  det_J = bf[eqn]->detJ;		/* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary 
   * with temperature, spatial coordinates, and species concentration.
   */

  /* CHECK FOR REMOVAL */
  acoustic_impedance( d_R, time );

  /* CHECK FOR REMOVAL */
  wave_number( d_k, time );

  alpha = light_absorption( d_alpha, time );

  /*acoustic_flux( q, d_q, time, py_eqn, py_var );*/

  switch(py_eqn)
   {
    case R_LIGHT_INTP:
         light_eqn = 0;
         Psign = 1.;
         break;
    case R_LIGHT_INTM:
         light_eqn = 1;
         Psign = -1.;
         break;
    case R_LIGHT_INTD:
         light_eqn = 2;
         Psign = 0.;
         break;
    default:
         EH(-1,"light intensity equation");
         break;
   }
  P = fv->poynt[light_eqn];
  grad_P = 0.;
  for(i=0 ; i<dim ; i++)
     {
      grad_P += svect[i]*fv->grad_poynt[light_eqn][i];
     }
  /*
   * Residuals___________________________________________________________
   */

  if ( af->Assemble_Residual )
    {
      eqn = py_eqn;
      peqn = upd->ep[eqn];
      var = py_var;
      for ( i=0; i<ei->dof[eqn]; i++)
	{
	  
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  diffusion = 0.;
	  if ( pd->e[eqn] & T_DIFFUSION )
	    {

	      diffusion += phi_i*(mucos*grad_P + Psign*alpha*P);
	      diffusion *= det_J * wt;
	      diffusion *= h3;
	      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
	    }

	  lec->R[peqn][i] += diffusion;
	}
    }


  /*
   * Jacobian terms...
   */

  if ( af->Assemble_Jacobian )
    {
      eqn   = py_eqn;
      peqn = upd->ep[eqn];
      for ( i=0; i<ei->dof[eqn]; i++)
	{
#if 1
          /* this is an optimization for xfem */
	  if ( xfem != NULL )
            {
	      int xfem_active, extended_dof, base_interp, base_dof;
	      xfem_dof_state( i, pd->i[eqn], ei->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      if ( extended_dof && !xfem_active ) continue;
            }
#endif
	  phi_i = bf[eqn]->phi[i];

	  /*
	   * Set up some preliminaries that are needed for the (a,i)
	   * equation for bunches of (b,j) column variables...
	   */

	  /*
	   * J_e_ap
	   */
	  var = py_var;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  diffusion = 0.;
		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
		      for ( p=0; p<VIM; p++)
			{
			  diffusion += grad_phi_j[p] * svect[p];
			}
		      diffusion *= phi_i*mucos;
		      diffusion += phi_i*Psign*alpha*phi_j;
		      diffusion *= h3*det_J*wt;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += diffusion;
		}
	    }
	  /*
	   * J_e_T
	   */
	  var = TEMPERATURE;
	  if ( pd->v[var] )
	    {
	      pvar = upd->vp[var];
	      for ( j=0; j<ei->dof[var]; j++)
		{
		  phi_j = bf[var]->phi[j];

		  for ( p=0; p<VIM; p++)
		    {
		      grad_phi_j[p] = bf[var]->grad_phi[j][p];
		    }

		  if ( pd->e[eqn] & T_DIFFUSION )
		    {
			  diffusion = phi_i*d_alpha->T[j]*P;
		      diffusion *= det_J * wt;
		      diffusion *= h3;
		      diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
		    }

		  lec->J[peqn][pvar][i][j] += diffusion;
		}
	    }

	  /*
	   * J_e_d
	   */
	  for ( b=0; b<dim; b++)
	    {
	      var = MESH_DISPLACEMENT1+b;
	      if ( pd->v[var] )
		{
		  pvar = upd->vp[var];
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

		      d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

			  /*
			   * multiple parts:
			   * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
			   *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
			   *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
			   *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
			   */
                      if ( pd->e[eqn] & T_DIFFUSION )
		        {
                          diff_b = 0.;
                          for ( p=0; p<VIM; p++)
			    {
			      diff_b += svect[p]*
                                 fv->d_grad_poynt_dmesh[light_eqn][p][b][j];
			    }
		          diff_b *= mucos;
                          diff_b += Psign*d_alpha->X[b][j]*P;
		          diff_b *= phi_i*det_J * h3 * wt;

			  diff_c = phi_i*(mucos*grad_P+Psign*alpha*P);
			  diff_c *= d_det_J_dmeshbj * h3 * wt;

			  diff_d = phi_i*(mucos*grad_P+Psign*alpha*P);
			  diff_d *= det_J * dh3dmesh_bj * wt;

			  diffusion = diff_b + diff_c + diff_d;

			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}


		      lec->J[peqn][pvar][i][j] += diffusion;
		    }
		}
	    }

	  /*
	   * J_e_c
	   */
	  var = MASS_FRACTION;
	  if ( pd->e[eqn] && pd->v[var] )
	    {
	      for ( w=0; w<pd->Num_Species_Eqn; w++)
		{
		  for ( j=0; j<ei->dof[var]; j++)
		    {
		      phi_j = bf[var]->phi[j];

		      if ( pd->e[eqn] & T_DIFFUSION )
			{
			  diffusion = phi_i*Psign*d_alpha->C[w][j]*P;
			  diffusion *= det_J * wt;
			  diffusion *= h3;
			  diffusion *= pd->etm[eqn][(LOG2_DIFFUSION)];
			}

		      lec->J[peqn][MAX_PROB_VAR + w][i][j] += diffusion;
		    }
		}
	    }

	}
    }

  return(status);
} /* end of assemble_poynting */
