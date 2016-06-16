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
 * $Id: mm_flux.c,v 5.25 2010-05-17 20:33:26 sarober Exp $
 */

/* Standard include files */
 
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
 
/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
 
#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h" 
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"
 
#include "mm_eh.h"
#include "mm_post_def.h"

#include "mm_std_models_shell.h"
#include "mm_std_models.h"

#define _MM_FLUX_C
#include "goma.h"

static int load_fv_sens
PROTO(( void ));

static int load_fv_grads_sens
PROTO(( void ));

/* evaluate_flux() -- sum up flux contributions along a sideset, print out
 *
 * Author:          P. R. Schunk
 * Date:            12 Jan 1996
 *
 * (o) Replaced EXODUS II I/O calls with references to data already read.
 *
 * (o) Considered distributed processing case where this SSID may appear
 *     on other processors and incompletely here.
 *
 * Modified: 1997/08/21 10:29 MDT pasacki@sandia.gov
 *   8/5/98 - Added normal, tangent, x, y, z, force calculations
 *          - volume, heat, species flux
 *	    - for right now, non-moving side sets (i.e. convective terms)
 *	things to do yet;
 *	(2) add diffusivity coeff models for terms for species flux
 *	(3) test heat flux, solid forces, species flux
 *	(4) write results to exodusII file as global variables
 *	(5) add viscoelastic stress terms
 *   10/99  - added sensitivity section, tabaer
 *   1/2002 - added sections on charged-species flux and electrical current, kschen 
 */

double
evaluate_flux(
			  const Exo_DB *exo, /* ptr to basic exodus ii mesh information */
			  const Dpi *dpi, /* distributed processing info */
			  const int side_set_id, /* on which SSID to evaluate flux */
			  const int quantity, /* to pick HEAT_FLUX, FORCE_NORMAL, etc. */
			  const char *qtity_str, /* quantity string */
			  const int blk_id,	/* material identification */
			  const int species_id, /* species identification */
			  const char *filenm, /* File name pointer */
			  const int profile_flag,  /*  flag for printing flux profiles  */
			  const double x[],	/* solution vector */
			  const double xdot[],	/* dx/dt vector */
			  double J_AC[], /* vector for augmenting condition sensitivities. 
			  May be NULL */
			  const double delta_t, /* time-step size */
			  const dbl time_value, /* current time */
			  const int print_flag)     /*  flag for printing results,1=print*/
			  {
  int j;			/* local index loop counter                 */
  int i;			/* Index for the local node number - row    */
  int ip = 0, a, b, c, p, w = -1;
  int mn;
  int var;
  int *n_dof=NULL;
  int dof_map[MDE];
  dbl H_lub; 
  dbl H_U, dH_U_dtime, H_L, dH_L_dtime;
  dbl dH_U_dX[DIM],dH_L_dX[DIM];
  dbl dH_U_dp, dH_U_ddh;
  dbl base_normal[DIM];

  double wt,weight;
  int err;			/* temp variable to hold diagnostic flags. */
  int ielem=-1;                 /* element number */
  int ip_total, ielem_type, gnn, ledof, matIndex;
  int num_local_nodes, iconnect_ptr, dim, ielem_dim, current_id;
  int nset_id, sset_id;
  double Tract[DIM], Torque[DIM], local_Torque[DIM];
  double Mag_real[DIM], Mag_imag[DIM], E_real[DIM], E_imag[DIM];

  /* 
   * Variables for vicosity and derivative 
   */
  dbl gamma[DIM][DIM];
  dbl mu = 0.0;
  dbl gamma_dot = 0.0;
  int elem_sign_org;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl q[DIM],dq_gradT[DIM][DIM];        /* heat flux vector, sensitvity */
  dbl dq_dX[DIM][DIM];
  dbl k=0;			        /* Thermal conductivity. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Thermal conductivity dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl Cp = 0.0;                               /* Heat capacity. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT d_Cp_struct; /* Heat capacity dependence. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp = &d_Cp_struct;

  dbl R_imped;                               /* Acoustic impedance */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl wnum, kR_inv;                               /* Acoustic wavenumber */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_wnum_struct; 
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_wnum = &d_wnum_struct;

  dbl dkdT[MDE];			/* Temperature derivative of electrical conductivity. */
  dbl dkdV[MDE];			/* Potential derivative of electrical conductivity. */
  dbl dkdC[MAX_CONC][MDE];		/* Concentration derivative of electrical conductivity. */
  dbl dkdX[DIM][MDE];   	       	/* Spatial derivatives of electrical conductivity. */


  double rho=0; /*  density variables */
  dbl e_theta[DIM] = {0.,0.,1.};  /* direction vector in theta direction for SWIRLING flow problems*/

  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  dbl dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  dbl dTT_dp_liq[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
				    to nodal porous liquid pressure*/
  dbl dTT_dp_gas[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
				    to nodal porous gas pressure*/
  dbl dTT_dporosity[DIM][DIM][MDE];/* Sensitivity of stress tensor... 
				    to nodal porosity*/
  dbl dTT_dsink_mass[DIM][DIM][MDE];
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];
  double elast_modulus;

  double I_T, II_T, III_T, coeff_a, coeff_b;
  double m_par = 0.0, theta1, evalue1 = 0.0, evalue2 = 0.0, evalue3 = 0.0;

  double es[MAX_PDIM][MAX_PDIM], efield[MAX_PDIM];	/* electric stress */
  double efield_sqr;				/* efield magnitude squared  */
  dbl perm = 0.0;		        /* electrical permittivity */
  dbl d_p_dT[MDE];		/* Temperature derivative of permittivity. */
  dbl d_p_dV[MDE];		/* Potential derivative of permittivity. */
  dbl d_p_dC[MAX_CONC][MDE];	/* Concentration derivative of permittivity. */
  dbl d_p_dX[DIM][MDE];         /* Spatial derivatives of permittivity. */

  double vs[MAX_PDIM][MAX_PDIM];	/* viscous stress */
  double dsigma_dx[DIM][MDE];		/* surface tension terms */
  double dsigmadT = 0, dsigmadC[MAX_CONC];

  double ves[MAX_PDIM][MAX_PDIM];	/* viscoelastic stress */
  int ve_mode;
  int v_s[MAX_MODES][DIM][DIM];

  double x_dot[MAX_PDIM];   /*    moving mesh stuff  */

  double s, t, u;             /* Gaussian-quadrature point locations          */
  double det;
  double local_flux = 0.0;    /* these are actually the integrated contact and*/
  double local_flux_conv = 0.0;       /*  convective fluxes  */
  double local_area= 0.0;
  double local_q = 0;               /* these are the actual fluxes  */
  double local_qconv = 0;

  double param[3] = {0.,0.,0.};

#ifdef PARALLEL
  double local_flux0 = 0.0;     /* old flux sum */
  double local_flux_conv0 = 0.0;/* old convective flux sum */
  double local_area0= 0.0;      /* old area sum */
  double proc_flux=0.0;         /* current flux sum on this proc */
  double proc_flux_conv=0.0;    /* current convective flux sum on this proc */
  double proc_area=0.0;         /* current area sum on this proc */
  double delta_flux=0.0;        /* increment of flux             */
  double delta_flux_conv=0.0;   /* increment of convective flux  */
  double delta_area=0.0;        /* increment of area */
  double global_flux=0.0;       /* flux sum over all procs */
  double global_flux_conv=0.0;  /* convective flux sum over all procs */
  double global_area=0.0;       /* area sum over all procs */
#endif

  double xi[DIM];             /* Local element coordinates of Gauss point. */

  int num_side_in_set;        /* returned number of sides(faces or edges) is  */
  int num_elem_in_set, num_node_in_set, num_nodes_on_side;
                              /* in the side set                              */
  int num_dist_fact_in_set;   /* returned number of distribution factors in   */
  int *elem_list;
  
  SGRID *element_search_grid=	NULL;	
  double surface_centroid[DIM];

  int id_side;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];

/* Local variables for the CHARGED_SPECIES_FLUX and CURRENT_FICKIAN cases            */ 
  double z[MAX_CONC];         /* species charge number                       */ 
  double kapta[MAX_CONC];     /* species electrical conductivity             */
  double d_kapta_dc[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];      /* dk/dc  */
  double d_kapta_dT[MAX_CONC];                                     /* dk/dT  */
  double d_kapta_dx[MAX_CONC][DIM];                                /* dk/dx  */
  double T = -1.0;                   /* electrolyte solution temperature            */
  const double R = 8.314;     /* Universal gas constant in units of J/mole K */
  const double FF = 96487.0;  /* Faraday's constant in units of C/equiv.     */
  double FRT;                 /* product of FF/R/T                           */

#ifndef PARALLEL
  static const char yo[] = "evaluate_flux";
#endif
 
/*  variables for adaptive quadrature weight calculation   */
 double ls_F[9],ad_wt[9], xf[2];
 int ierr = 0, wt_type;
 const int Jac_state = af->Assemble_Jacobian;

  for (a=0; a < DIM; a++) {
    base_normal[a] = 0;
  }
  memset( Torque, 0, sizeof(double)*DIM );

  /* load eqn and variable number in tensor form */
  err = stress_eqn_pointer(v_s);


  af->Assemble_Jacobian = TRUE;
  /* first right time stamp or run stamp to separate the sets */

        if (print_flag && ProcID == 0) {
        FILE  *jfp;
        if( (jfp=fopen(filenm,"a")) != NULL)
                {
      fprintf(jfp,"Time/iteration = %e \n", time_value);
      fprintf(jfp," %s  Side_Set  %d Block  %d Species %d\n", qtity_str,side_set_id,blk_id,species_id);
                if(profile_flag)
                   {
                    fprintf(jfp," x  y  z  flux  flux_conv");
                    if(profile_flag & 1)fprintf(jfp," shear_rate viscosity");
                    if(profile_flag & 2)fprintf(jfp," evalue1 evalue2 evalue3");
                    if(profile_flag & 4)fprintf(jfp," surf_tens");
                    if(profile_flag & 8)fprintf(jfp," temperature pressure");
                    if(profile_flag & 16)fprintf(jfp," x-normal y-normal z-normal");
                    fprintf(jfp," \n");
                   }
      fflush(jfp);
      fclose(jfp);
                }
        }

  current_id = in_list(side_set_id, 0, exo->num_side_sets, exo->ss_id);

  
  /*
   * Not finding this side set in the mesh is immediately serious if we're
   * running serial. Parallel execution is OK with this particular processor
   * not finding it, though, just so long as *some* processor(s) find this
   * SSID.
   */

  if ( current_id != -1 )
    {
      num_side_in_set      = exo->ss_num_sides[current_id];
      
      num_dist_fact_in_set = exo->ss_num_distfacts[current_id];
      
      num_elem_in_set      = num_side_in_set;

      num_node_in_set      = num_dist_fact_in_set;
      
      num_nodes_on_side    = num_node_in_set/num_elem_in_set; /* Well... */
      
      elem_list            = &exo->ss_elem_list[exo->ss_elem_index[current_id]];

      /* Now start element sweep */

      for (i=0,local_flux = local_flux_conv = local_area = 0.; 
	   i< num_elem_in_set; 
	   i++)
	{         
	  ei->ielem = elem_list[i];

	  mn = find_mat_number(ei->ielem, exo);


	/**    only do integration if the requested material number  **/

	if ( blk_id == -1 || 
             (mn == map_mat_index(blk_id) && dpi->elem_owner[ elem_list[i] ] == ProcID)) {
	  /*
	   * Yes, "x" gets recycled like garbage. Fortunately, this 
	   * routine should not write onto "x"...
	   */

	  err = load_elem_dofptr( elem_list[i], 
				  (Exo_DB*) exo,
				  (dbl *) x,
				  (dbl *) x,
				  (dbl *) xdot,
				  (dbl *) xdot,
				  (dbl *) x,
				  0);
	  EH(err, "load_elem_dofptr");

	  err = bf_mp_init(pd);
	  EH(err, "bf_mp_init");
  
	      iconnect_ptr        = ei->iconnect_ptr;
	      ielem_type          = ei->ielem_type;
	      ip_total            = elem_info(NQUAD_SURF, ielem_type);
	      num_local_nodes     = ei->num_local_nodes;
	      ielem_dim           = ei->ielem_dim;
	      
	 
	      id_side = find_id_side (ei->ielem, num_nodes_on_side,
				      &exo->ss_node_list[current_id]
				      [num_nodes_on_side*i],
				      id_local_elem_coord, exo);

	      /* Calculates the ID side correctly for tets */
	      if ( ielem_type == TRILINEAR_TET ) id_side = find_id_side_SS(ei->ielem, current_id, exo);

		  if( ls != NULL && ielem_dim == 2  &&
			(quantity == POS_LS_FLUX || quantity == NEG_LS_FLUX ||
				quantity == DELTA || quantity == LS_DCA ))
		  {
                    if ( ls->var != LS )
			WH(-1,"Level-set variable is not LS!");
			  switch (id_side)
			  {
				  case 1:
					  ls_F[0]=*esp_old->F[0];
					  ls_F[1]=*esp_old->F[1];
					  ls_F[2]=*esp_old->F[4];
					  break;
				  case 2:
					  ls_F[0]=*esp_old->F[1];
					  ls_F[1]=*esp_old->F[2];
					  ls_F[2]=*esp_old->F[5];
					  break;
				  case 3:
					  ls_F[0]=*esp_old->F[3];
					  ls_F[1]=*esp_old->F[2];
					  ls_F[2]=*esp_old->F[6];
					  break;
				  case 4:
					  ls_F[0]=*esp_old->F[0];
					  ls_F[1]=*esp_old->F[3];
					  ls_F[2]=*esp_old->F[7];
					  break;
			  }
			if(species_id > 0)
			{
			  var=species_id-1;
			  switch (id_side)
			  {
				  case 1:
					  ls_F[0]=*esp_old->pF[var][0];
					  ls_F[1]=*esp_old->pF[var][1];
					  ls_F[2]=*esp_old->pF[var][4];
					  break;
				  case 2:
					  ls_F[0]=*esp_old->pF[var][1];
					  ls_F[1]=*esp_old->pF[var][2];
					  ls_F[2]=*esp_old->pF[var][5];
					  break;
				  case 3:
					  ls_F[0]=*esp_old->pF[var][3];
					  ls_F[1]=*esp_old->pF[var][2];
					  ls_F[2]=*esp_old->pF[var][6];
					  break;
				  case 4:
					  ls_F[0]=*esp_old->pF[var][0];
					  ls_F[1]=*esp_old->pF[var][3];
					  ls_F[2]=*esp_old->pF[var][7];
					  break;
			  }
			}
			  
			  switch( quantity )
			  {
                                  case POS_LS_FLUX:
					  wt_type = 2;
					  break;
				  case NEG_LS_FLUX:
					  wt_type = 2;
					  ls_F[0] = -ls_F[0];
					  ls_F[1] = -ls_F[1];
					  ls_F[2] = -ls_F[2];
					  break;
                                  case DELTA:
                                  case LS_DCA:
					  wt_type = 3;
					  break;
                                  default:
                                          wt_type = 2;
                                          break;
			  }
			  ierr = adaptive_weight( ad_wt, ip_total, ielem_dim-1, ls_F, 0.0,
									  wt_type, ei->ielem_type);
			  if(ierr == -1)printf("adaptive wt problem %d %g %g %g\n",ierr,ls_F[0],ls_F[1],ls_F[2]);
		  }
		  
		  if( ls!=NULL && ls->elem_overlap_state && ls->Integration_Depth > 0 )
		  {
			  Subgrid_Int.active = TRUE;
			  
			  get_subgrid_integration_pts ( Subgrid_Tree, &element_search_grid,
											&Subgrid_Int.s, &Subgrid_Int.wt, ls->Length_Scale );
			  
			  find_surf_center_st ( ielem_type, id_side, ielem_dim, surface_centroid, &s, &t );
			  
			  ip_total = gather_surface_subgrid_integration_pts( element_search_grid, id_side, 
																 surface_centroid, Subgrid_Int.s, Subgrid_Int.wt, 0 );
			  
		  }
		  

 
	  /* Surface integration over element */
	  
	      for (ip = 0; ip < ip_total; ip++)
		  {
			  /*   zero local fluxes  */
			  local_q = local_qconv = 0.;
			  memset(local_Torque, 0, DIM*sizeof(double) );
			  
			  if( ls != NULL && ls->elem_overlap_state && ls->Integration_Depth > 0 )
			  {
				  ls->Elem_Sign = 0;
				  xi[0] = Subgrid_Int.s[ip][0];
				  xi[1] = Subgrid_Int.s[ip][1];
				  xi[2] = Subgrid_Int.s[ip][2];
				  
				  s = 1.e30; 
				  t = 1.e30;
				  wt = Subgrid_Int.wt[ip];
			  }
			  else
			  {
				  
				  /* find the quadrature point locations for current ip */
				  find_surf_st (ip, ielem_type, id_side, ielem_dim, xi, &s, &t, &u);
				  
				  /* find the quadrature weight for current ip */
				  wt = Gq_surf_weight (ip, ielem_type);    
			  }
			  
			  /* ****************************************/
			  err = load_basis_functions( xi, bfd);
			  EH( err, "problem from load_basis_functions");
			  
			  err = beer_belly();
			  EH( err, "beer_belly");
			  
			  /* precalculate variables at  current integration pt.*/
			  err = load_fv();
			  EH( err, "load_fv");
			  
			  err = load_bf_grad();
			  EH( err, "load_bf_grad");
			  
			  err = load_bf_mesh_derivs();
			  EH( err, "load_bf_mesh_derivs");
			  
			  surface_determinant_and_normal (ei->ielem, 
				iconnect_ptr, num_local_nodes, ielem_dim-1,  
				id_side, num_nodes_on_side,id_local_elem_coord);
			  
			  /*
			   * Load up physical space gradients of field variables at this
			   * Gauss point.
			   */
			  err = load_fv_grads();
			  EH( err, "load_fv_grads");
			  
			  err = load_fv_mesh_derivs(1);
			  EH( err, "load_fv_mesh_derivs");
			  
			  if (TimeIntegration != STEADY && pd->e[MESH_DISPLACEMENT1]) {
				  for (j = 0; j < VIM; j++) {
					  x_dot[j] = fv_dot->x[j];
				  } 
			  } else {
				  for (j = 0; j < VIM; j++) {
					  x_dot[j] = 0.;
				  }
			  }
			  
			  do_LSA_mods(LSA_SURFACE);

	    if (ielem_dim !=3)
	      {
		calc_surf_tangent (ei->ielem, iconnect_ptr, 
				   num_local_nodes, ielem_dim-1,
				   num_nodes_on_side,
				   id_local_elem_coord);
	      }

            weight = wt;
	    det = fv->sdet*fv->h3;

	    if(upd->CoordinateSystem == CYLINDRICAL ||
               upd->CoordinateSystem == SWIRLING   )
	      {
		weight = 2.*M_PIE*wt;
		det = fv->sdet;
	      }

/*	get density  */

	    rho = density( NULL, time_value );

/**   get solid stresses  **/

		  dim   = pd_glob[0]->Num_Dim;
		  
		  if(pd->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
		    {
		      err = belly_flop(elc->lame_mu);
		      EH(err, "error in belly flop");
		      if (err == 2) exit(-1);
		/*
		 * Total mesh stress tensor...
		 *
		 * Guess what the 4 new args will have to be and bloat accordingly.
		 */

		      err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, 
					       dTT_dp_liq, dTT_dp_gas, dTT_dporosity, dTT_dsink_mass,   
					       dTT_dT, dTT_dmax_strain,  dTT_dcur_strain, elc->lame_mu, elc->lame_lambda,
					       delta_t, ei->ielem, ip, ip_total);


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
			  if (dim == 2)
                          {
			    elast_modulus = elc->lame_mu;
				if (cr->MeshFluxModel == NONLINEAR ||
				    cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
				    cr->MeshFluxModel == INCOMP_PSTRAIN )
				  TT[2][2] = (1. - pow(fv->volume_change,2./3.)) * elast_modulus - fv->P;
			    /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
				else  TT[2][2] = 0.;
			    TT[1][2] = 0.;
			    TT[0][2] = 0.;
			  }
			}
   if(profile_flag & 2)
      {
      /*  try for Trig solution of cubic equation     */
      I_T = TT[0][0]+TT[1][1]+TT[2][2];
      II_T = TT[0][0]*TT[1][1]+TT[0][0]*TT[2][2]+TT[1][1]*TT[2][2]
              -(SQUARE(TT[0][1])+SQUARE(TT[0][2])+SQUARE(TT[1][2]));
      III_T = TT[0][0]*TT[1][1]*TT[2][2]+2.*(TT[0][1]*TT[1][2]*TT[0][2])
              -TT[0][0]*SQUARE(TT[1][2])-TT[1][1]*SQUARE(TT[0][2])
              -TT[2][2]*SQUARE(TT[0][1]);
      coeff_a = (3.*II_T - SQUARE(I_T))/3.;
      coeff_b = (2.*(-I_T)*SQUARE(I_T)-9.*(-I_T)*II_T+27.*(-III_T))/27.;
      if(coeff_a > 0)
              {fprintf(stderr,"trouble with principal stress - imaginary roots %g %g\n",coeff_a,coeff_b);}
      else
              {m_par = 2.*sqrt(-coeff_a/3.);}
      theta1 = acos(3.*coeff_b/(coeff_a*m_par))/3.;
      evalue1 = m_par*cos(theta1) + I_T/3.;
      evalue2 = m_par*cos(theta1+2.*M_PIE/3.) + I_T/3.;
      evalue3 = m_par*cos(theta1+4.*M_PIE/3.) + I_T/3.;
      theta1 = evalue2;
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue1=theta1;}
      theta1 = evalue3;
      if(fabs(theta1) > fabs(evalue2))
             { evalue3=evalue2; evalue2=theta1;}
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue2=theta1;}
      }
		    } /* end of STRESS_TENSOR */

/* calculate real-solid stress here !!*/
		  if(pd->e[R_SOLID1] && cr->MeshMotion != ARBITRARY)
		    {
		      mu = elc_rs->lame_mu;
		      err = belly_flop_rs(mu);
		      EH(err, "error in belly flop");
		      if (err == 2) return(err);
		  /*
		   * Total mesh stress tensor...
		   */
 		      err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp, 
						dTT_dc,  dTT_dp_liq, dTT_dp_gas, dTT_dporosity, 
						dTT_dT, dTT_dmax_strain, elc_rs->lame_mu, elc_rs->lame_lambda);
		      if (dim == 2)
			{
			  elast_modulus = elc_rs->lame_mu;
			  
			  if (cr->RealSolidFluxModel == NONLINEAR ||
			      cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
			      cr->RealSolidFluxModel == INCOMP_PSTRAIN )
			    TT[2][2] = (1. - pow(fv->volume_change,2./3.)) * elast_modulus - fv->P;
		      /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
			  else  TT[2][2] = 0.;
		      
			  TT[1][2] = 0.;
			  TT[0][2] = 0.;
			}
   if(profile_flag & 2)
      {
      /*  try for Trig solution of cubic equation     */
      I_T = TT[0][0]+TT[1][1]+TT[2][2];
      II_T = TT[0][0]*TT[1][1]+TT[0][0]*TT[2][2]+TT[1][1]*TT[2][2]
              -(SQUARE(TT[0][1])+SQUARE(TT[0][2])+SQUARE(TT[1][2]));
      III_T = TT[0][0]*TT[1][1]*TT[2][2]+2.*(TT[0][1]*TT[1][2]*TT[0][2])
              -TT[0][0]*SQUARE(TT[1][2])-TT[1][1]*SQUARE(TT[0][2])
              -TT[2][2]*SQUARE(TT[0][1]);
      coeff_a = (3.*II_T - SQUARE(I_T))/3.;
      coeff_b = (2.*(-I_T)*SQUARE(I_T)-9.*(-I_T)*II_T+27.*(-III_T))/27.;
      if(coeff_a > 0)
              {fprintf(stderr,"trouble with principal stress - imaginary roots %g %g\n",coeff_a,coeff_b);}
      else
              {m_par = 2.*sqrt(-coeff_a/3.);}
      theta1 = acos(3.*coeff_b/(coeff_a*m_par))/3.;
      evalue1 = m_par*cos(theta1) + I_T/3.;
      evalue2 = m_par*cos(theta1+2.*M_PIE/3.) + I_T/3.;
      evalue3 = m_par*cos(theta1+4.*M_PIE/3.) + I_T/3.;
      theta1 = evalue2;
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue1=theta1;}
      theta1 = evalue3;
      if(fabs(theta1) > fabs(evalue2))
             { evalue3=evalue2; evalue2=theta1;}
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue2=theta1;}
      }
		    } /* end of REAL_STRESS_TENSOR */

		/*
		 *  viscoelastic stress tensor
		 *  assume only EVSS_F formulation for now
		 */
  			memset( ves, 0, sizeof(dbl)*DIM*DIM);
  			if ( pd->v[POLYMER_STRESS11] )
    			  {
			    for ( a=0; a<VIM; a++)
        	              {
			        for ( b=0; b<VIM; b++)
            			   {
			             for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
                		       {
                  			ves[a][b] += fv->S[ve_mode][a][b];
                			}
            			   }
        		      }
    			   }

		/*
		 * OK, let's simplify things by computing the viscous
		 * stress tensor upfront.
		 */
  			memset( vs, 0, sizeof(dbl)*DIM*DIM);
  			memset( gamma, 0, sizeof(dbl)*DIM*DIM);
		      if(cr->MeshMotion == ARBITRARY)
    			 {
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
				  gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
				}
			    }
			  
			  mu = viscosity(gn, gamma, d_mu);

			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
				  vs[a][b] = mu*gamma[a][b]-fv->P*delta(a,b);
				}
			    }
                          if(profile_flag & 1)
                             {
   				gamma_dot = 0.;
   				for ( a=0; a<VIM; a++)
     				   {
       					for ( b=0; b<VIM; b++)
         				  {
           				   gamma_dot +=  gamma[a][b] * gamma[b][a];
         				  }
     				   }
   				gamma_dot  =  sqrt(gamma_dot/2.);
                             }
   if(profile_flag & 2)
      {
  			memset( TT, 0, sizeof(dbl)*DIM*DIM);
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
				  TT[a][b] = vs[a][b]+ves[a][b];
				}
			    }
      /*  Principal Stresses - try for Trig solution of cubic equation     */
      I_T = TT[0][0]+TT[1][1]+TT[2][2];
      II_T = TT[0][0]*TT[1][1]+TT[0][0]*TT[2][2]+TT[1][1]*TT[2][2]
              -(SQUARE(TT[0][1])+SQUARE(TT[0][2])+SQUARE(TT[1][2]));
      III_T = TT[0][0]*TT[1][1]*TT[2][2]+2.*(TT[0][1]*TT[1][2]*TT[0][2])
              -TT[0][0]*SQUARE(TT[1][2])-TT[1][1]*SQUARE(TT[0][2])
              -TT[2][2]*SQUARE(TT[0][1]);
      coeff_a = (3.*II_T - SQUARE(I_T))/3.;
      coeff_b = (2.*(-I_T)*SQUARE(I_T)-9.*(-I_T)*II_T+27.*(-III_T))/27.;
      if(coeff_a > 0)
              {fprintf(stderr,"trouble with principal stress - imaginary roots %g %g\n",coeff_a,coeff_b);}
      else
              {m_par = 2.*sqrt(-coeff_a/3.);}
      theta1 = acos(3.*coeff_b/(coeff_a*m_par))/3.;
      evalue1 = m_par*cos(theta1) + I_T/3.;
      evalue2 = m_par*cos(theta1+2.*M_PIE/3.) + I_T/3.;
      evalue3 = m_par*cos(theta1+4.*M_PIE/3.) + I_T/3.;
      theta1 = evalue2;
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue1=theta1;}
      theta1 = evalue3;
      if(fabs(theta1) > fabs(evalue2))
             { evalue3=evalue2; evalue2=theta1;}
      if(fabs(theta1) > fabs(evalue1))
             { evalue2=evalue1; evalue2=theta1;}
      }
    			 }
 		/*
 		 * computing the electric stress tensor upfront.
 		 */
  			memset( es, 0, sizeof(dbl)*DIM*DIM);
		      if(pd->e[R_POTENTIAL])
     			 {
 			  efield_sqr = 0.0;
 			  for ( a=0; a<VIM; a++)
 			    {
 				  efield[a] = - fv->grad_V[a];
 				  efield_sqr += efield[a]*efield[a];
 			    }
 
 			  for ( a=0; a<VIM; a++)
 			    {
 			      for ( b=0; b<VIM; b++)
 				{
 				  es[a][b] = efield[a]*efield[b] - 0.5*efield_sqr*delta(a,b);
 				}
 			    }
 		      perm   = mp->permittivity;
 	              if( J_AC != NULL)
 		        {
 			  memset( d_p_dV,0, sizeof(double)*MDE);
 			  memset( d_p_dT,0, sizeof(double)*MDE);
 			  memset( d_p_dX,0, sizeof(double)*MDE*DIM);
 			  memset( d_p_dC,0, sizeof(double)*MAX_CONC*MDE);
 			}
     			 }
			  
 		/*
 		 * load surface tension for variable models
 		 */
 		      if (mp->SurfaceTensionModel == CONSTANT) {
                            dsigmadT = 0.0;
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = 0.0;}
                            }
 		      if (mp->SurfaceTensionModel != CONSTANT) {
 		            load_surface_tension(dsigma_dx);
 		            }
 		      if (mp->SurfaceTensionModel == USER) {
                            dsigmadT = mp->d_surface_tension[TEMPERATURE];
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = mp->d_surface_tension[MAX_VARIABLE_TYPES+a];}
                            }

/*		  if (pd->e[R_POR_LIQ_PRES])
		    {
		       if (mp->PorousMediaType == POROUS_SATURATED) { 
 			err = get_porous_fully_sat_terms(&pm_terms, tt, dt); 
 			EH(err,"problem in getting the species terms"); 
 			if (neg_elem_volume) return(err);
 		      } else if (mp->PorousMediaType == POROUS_UNSATURATED || 
 				 mp->PorousMediaType == POROUS_TWO_PHASE) { 
 			err = get_porous_part_sat_terms(&pm_terms, tt, dt); 
 			EH(err,"problem in getting the species terms"); 
 			if (neg_elem_volume) return(err); 
 		      } 
		    }*/

		  switch( quantity )
		    {
		    case AREA:
		      
		      local_flux += weight*fv->sdet;
		      local_flux_conv = 0.0;

		      break;
		      
                    case VOL_REVOLUTION:

                      local_q +=  0.5 * fv->x[1] * fabs(fv->snormal[1]) ;
                      local_flux += weight*det* local_q ;
		      local_flux_conv = 0.0;
/*                      local_flux_conv += weight*det*mp->surface_tension ;*/

                      break;

		    case HEAT_FLUX:
	      
		  /* finally we can add up the outgoing flux and area */
		  /* but first evaluate properties if variable */
		      k = conductivity( d_k, time_value );
		      Cp = heat_capacity( d_Cp, time_value );

                      if ( cr->HeatFluxModel == CR_HF_FOURIER_0 )
                      {
                      for (a=0; a<VIM; a++)
                        {
                          local_q +=  -k * fv->snormal[a] * fv->grad_T[a] ;
                          local_qconv +=   rho * Cp*fv->T* fv->snormal[a] * (fv->v[a]-x_dot[a]) ;
                        }
                      }  else if ( cr->HeatFluxModel == CR_HF_USER )
                        {
                        double delP=0.;
#if defined SECOR_HEAT_FLUX
			double *hpar, h, dh_dX[DIM], Vb[DIM],Vt[DIM];
			double dq_dVb[DIM][DIM], dq_dVt[DIM][DIM];

			hpar = &mp->u_thermal_conductivity[0];
			h = hpar[0] + hpar[4]*fv->x[0]
			  + (hpar[1]-hpar[5]*fv->x[0])*(hpar[3]-fv->x[1])
			  + 0.5*hpar[2]*SQUARE(hpar[3]-fv->x[1]);

			dh_dX[0] = hpar[4] - hpar[5]*(hpar[3]-fv->x[1]);
			dh_dX[1] = hpar[5]*fv->x[0]-hpar[1] - hpar[2]*(hpar[3]-fv->x[1]);

			/*     velocities of bottom and top surfaces   */
			Vb[0] = mp->u_heat_capacity[0];
			Vb[1] = mp->u_heat_capacity[1];
			Vt[0] = mp->u_heat_capacity[2];
			Vt[1] = mp->u_heat_capacity[3];

			usr_heat_flux(fv->grad_T, q, dq_gradT, dq_dX, time_value, h, dh_dX, Vb, Vt
                             ,dq_dVb, dq_dVt);
#else
			usr_heat_flux(fv->grad_T, q, dq_gradT, dq_dX, time_value);
			printf("untested\n");
			exit(-1);
#endif

			for (a=0; a<VIM; a++)
			  {
                                delP += fv->snormal[a] * fv->grad_T[a];
				local_q +=  fv->snormal[a] * q[a] ;
                                for (b=0; b<VIM; b++)
				  {
				    local_qconv +=  fv->snormal[a] * dq_gradT[a][b]
				      * fv->snormal[b];
				  }
                          }
	       	/*                                local_qconv *= delP/fv->T;  */
                        }
		      local_flux +=  weight * det * local_q ;
		      local_flux_conv +=   weight * det * local_qconv ;
		      break;

		    case ACOUSTIC_FLUX_NORMAL:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                      for (a=0; a<dim; a++)
                        {
                          local_q +=  -kR_inv * fv->snormal[a] * fv->grad_api[a] ;
                          local_qconv +=  kR_inv * fv->snormal[a] * fv->grad_apr[a] ;
                        }
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case ACOUSTIC_FLUX_TANGENT1:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                      for (a=0; a<dim; a++)
                        {
                          local_q +=  -kR_inv * fv->stangent[0][a] * fv->grad_api[a] ;
                          local_qconv +=  kR_inv * fv->stangent[0][a] * fv->grad_apr[a] ;
                        }
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case ACOUSTIC_FLUX_TANGENT2:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                      for (a=0; a<dim; a++)
                        {
                          local_q +=  -kR_inv * fv->stangent[1][a] * fv->grad_api[a] ;
                          local_qconv +=  kR_inv * fv->stangent[1][a] * fv->grad_apr[a] ;
                        }
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case ACOUSTIC_FLUX_X:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                          local_q +=  -kR_inv * fv->grad_api[0] ;
                          local_qconv +=  kR_inv * fv->grad_apr[0] ;
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case ACOUSTIC_FLUX_Y:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                          local_q +=  -kR_inv * fv->grad_api[1] ;
                          local_qconv +=  kR_inv * fv->grad_apr[1] ;
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case ACOUSTIC_FLUX_Z:
	      
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);

                          local_q +=  -kR_inv * fv->grad_api[2] ;
                          local_qconv +=  kR_inv * fv->grad_apr[2] ;
                          local_flux +=  weight * det * local_q ;
                          local_flux_conv +=   weight * det * local_qconv ;
                      break;

		    case VOLUME_FLUX:
                      for(a=0; a<VIM; a++)
                        {
			  if (cr->MeshMotion == ARBITRARY)
			    local_q +=  fv->snormal[a]*( fv->v[a]-x_dot[a] ) ;
			  else if (pd->v[MESH_DISPLACEMENT1])
			    local_q +=  fv->snormal[a]*( fv->d[a] ) ;
			  else
			    EH(-1,"Inconsistency in volume-flux specification. Contact Developers");
                        }
                          local_flux += weight*det* local_q ;
		      break;

		    case SHELL_VOLUME_FLUX:

		      // shell_determinant_and_normal(ei->ielem, ei->iconnect_ptr, ei->num_local_nodes,
                      //         ei->ielem_dim, 1);

		      /* First save local normal to edge because lubrication_shell_init changes it */
		      for(a=0; a<VIM; a++) base_normal[a]=fv->snormal[a];

		      n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
		      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

		      H_lub = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time_value, 0); 
		      switch ( mp->FSIModel ) {
		      case FSI_MESH_CONTINUUM:
		      case FSI_MESH_UNDEF:
			for ( a = 0; a < dim; a++) {
			  H_lub -= fv->snormal[a] * fv->d[a];
			}
			break; 
		      case FSI_REALSOLID_CONTINUUM:
			for ( a = 0; a < dim; a++) {
			  H_lub -= fv->snormal[a] * fv->d_rs[a];
			}
			break;
		      }
		      /* Calculate the flow rate and its sensitivties */

		      calculate_lub_q_v(R_LUBP, time_value, 0, xi, exo);

                      for(a=0; a<VIM; a++)
                        {
			    local_q +=  base_normal[a]*LubAux->v_avg[a] * H_lub;
                        }
                          local_flux += weight*det* local_q ;
			  /* clean-up */
			  safe_free((void *) n_dof);
		      break;

     		    case NEG_LS_FLUX:
     		    case POS_LS_FLUX:
 
       			{
 			double H = 0.0;
 			load_lsi(ls->Length_Scale);
 
 			H =  ( quantity == POS_LS_FLUX ? lsi->H : ( 1.0 - lsi->H ) ) ;
 	
 			for( a=0; a<dim; a++)
 				{
 	  		local_q += H * (fv->v[a]-x_dot[a]) * fv->snormal[a] ;
  	  		local_qconv += (fv->v[a]-x_dot[a]) * fv->snormal[a] ;
 				}
 	  		local_flux += weight * det * local_q ;
  	  		local_flux_conv += ad_wt[ip_total-ip-1]*det*local_qconv;
 
       			}
       			break;
 
      
		    case PVELOCITY1:
		    case PVELOCITY2:
		    case PVELOCITY3:
                      for(a=0; a<VIM; a++)
                        {
                          local_q +=  fv->snormal[a]*fv->pv[a] ;
                        }
                          local_flux +=  weight*det*local_q ;
		      break;
		  
		    case EXT_VELOCITY: /* This is already a normal velocity */
		      local_q +=  fv->ext_v;
		      local_flux +=  weight*det*local_q ;
		      break;

		    case EFIELD1:
		    case EFIELD2:
		    case EFIELD3:
                      for(a=0; a<VIM; a++)
                        {
                          local_q +=  fv->snormal[a]*fv->E_field[a] ;
                        }
		      local_flux +=  weight*det*local_q ;
		      break;

		    case SPECIES_FLUX:

                      for (a = 0; a < VIM; a++) {
	                if ( cr->MassFluxModel == FICKIAN  ||
	                     cr->MassFluxModel == STEFAN_MAXWELL  ||
	                     cr->MassFluxModel == STEFAN_MAXWELL_CHARGED  ||
	                     cr->MassFluxModel == STEFAN_MAXWELL_VOLUME )
	                    {
	                     if ( Diffusivity() )  EH( -1, "Error in Diffusivity.");

	                     local_q += fv->snormal[a] * 
                               (-mp->diffusivity[species_id]*fv->grad_c[species_id][a]);
	                    }
	                else if ( cr->MassFluxModel == GENERALIZED_FICKIAN)
	                    {
	                     if ( Generalized_Diffusivity() )  EH( -1, "Error in Diffusivity.");
	                     for (w=0; w<pd->Num_Species_Eqn; w++)
	                       {
	                        local_q +=  fv->snormal[a] * 
                                           (-mp->diffusivity_gen_fick[species_id][w]
                                              *fv->grad_c[w][a]);
	                       }
	                    }
	                else  if ( cr->MassFluxModel == DARCY )
	                    { /* diffusion induced convection is zero */ }
	                else
	                    { EH( -1, "Unimplemented mass flux constitutive relation."); }
                          local_qconv += (fv->snormal[a]*(fv->v[a]-x_dot[a])
                                         *fv->c[species_id] );
                        }
                          local_qconv = 0;
                          local_flux +=  weight*det*local_q;
                          local_flux_conv += weight*det*local_qconv;
		      break;

		    case SPECIES_FLUX_REVOLUTION:

		      Diffusivity();

                      for(a=0; a<VIM; a++)
                        {
                          local_q += ( -mp->diffusivity[species_id]*
                                  fv->snormal[a]*fv->grad_c[species_id][a] );
                          local_qconv += ( fv->snormal[a]*(fv->v[a]-x_dot[a])*fv->c[species_id] );
                        }
                      local_q *=  0.5 * fv->x[1] * fabs(fv->snormal[1]) ;
                      local_qconv *=  0.5 * fv->x[1] * fabs(fv->snormal[1]) ;
                          local_flux +=  weight*det*local_q;
                          local_flux_conv += weight*det*local_qconv;
		      break;

		    case CHARGED_SPECIES_FLUX:

		      Diffusivity();

                      z[species_id] = mp->charge_number[species_id];
                      if (mp->SolutionTemperatureModel == CONSTANT)  
                        {
                          T = mp->solution_temperature;
                        } 
                      else 
                        {
                          EH(-1, "Solution-temperature model not yet implemented");
                        }
                      /* set solution temperature to 298 K if it is zero - safety feature */
                      if (T == 0.0)
                        {
                          T = 298.0;
                          fprintf(stderr, "Warning!: a default electrolyte temperature of 298 K is being used!");
                        }
 
                      FRT=FF/R/T;
                      kapta[species_id] = FRT*z[species_id]*mp->diffusivity[species_id]*fv->c[species_id]; 
		      for (w = 0; w < pd->Num_Species_Eqn; w++) 
                        {
                          d_kapta_dc[species_id][MAX_VARIABLE_TYPES + w] = FRT*z[species_id]*
				                                           mp->diffusivity[species_id];
                        }
                          d_kapta_dT[species_id] = (-FRT/T)*z[species_id]*mp->diffusivity[species_id]*
                                                   fv->c[species_id];
                       for(a=0; a<VIM; a++)
                        {
                          d_kapta_dx[species_id][a] = 0.0;
                        }

                      for(a=0; a<VIM; a++)
                        {
                          local_q += (-mp->diffusivity[species_id]*fv->snormal[a]*fv->grad_c[species_id][a] 
                                      -kapta[species_id]*fv->snormal[a]*fv->grad_V[a]);
                          local_qconv += ( fv->snormal[a]*(fv->v[a]-x_dot[a])*fv->c[species_id] );
                        }
                          local_flux +=  weight*det*local_q;
                          local_flux_conv += weight*det*local_qconv;
		      break;

		    case CURRENT_FICKIAN:

		      Diffusivity();

                      z[species_id] = mp->charge_number[species_id]; 

                      if (mp->SolutionTemperatureModel == CONSTANT)  
                        {
                          T = mp->solution_temperature;
                        } 
                      else 
                        {
                          EH(-1, "Solution-temperature model not yet implemented");
                        }
                      /* set solution temperature to 298 K if it is zero - safety feature */
                      if (T == 0.0)
                        {
                          T = 298.0;
                          fprintf(stderr, "Warning!: a default electrolyte temperature of 298 K is being used!");
                        }

                      FRT=FF/R/T;
                      kapta[species_id] = FRT*z[species_id]*mp->diffusivity[species_id]*fv->c[species_id]; 
		      for (w = 0; w < pd->Num_Species_Eqn; w++) 
                        {
                          d_kapta_dc[species_id][MAX_VARIABLE_TYPES + w] = FRT*z[species_id]*
				                                           mp->diffusivity[species_id];
                        }
                          d_kapta_dT[species_id] = (-FRT/T)*z[species_id]*mp->diffusivity[species_id]*
                                                   fv->c[species_id];
                       for(a=0; a<VIM; a++)
                        {
                          d_kapta_dx[species_id][a] = 0.0;
                        }

                      for(a=0; a<VIM; a++)
                        {
                          local_q += (-mp->diffusivity[species_id]*fv->snormal[a]*fv->grad_c[species_id][a] 
                                      -kapta[species_id]*fv->snormal[a]*fv->grad_V[a]);
                          local_qconv += ( fv->snormal[a]*(fv->v[a]-x_dot[a])*fv->c[species_id] );
                        }
                          local_flux +=  weight*det*local_q;
                          local_flux_conv += weight*det*local_qconv;
		      break;

		    case CURRENT:
		      k   = mp->electrical_conductivity;
	              if( J_AC != NULL)
		        {
			  memset( dkdV,0, sizeof(double)*MDE);
			  memset( dkdT,0, sizeof(double)*MDE);
			  memset( dkdX,0, sizeof(double)*MDE*DIM);
			  memset( dkdC,0, sizeof(double)*MAX_CONC*MDE);
			}
                      for(a=0; a<VIM; a++)
                        {
                          local_q += - k * fv->snormal[a]*fv->grad_V[a]; 
                        }
                          local_flux +=  weight*det*local_q;
                          local_flux_conv = 0.0;
		      break;

 		    case ELEC_FORCE_NORMAL:
 			  for ( a=0; a<VIM; a++)
 			    {
 			      for ( b=0; b<VIM; b++)
 				{
                                 local_q += (fv->snormal[a]*perm*es[a][b]*fv->snormal[b]);
 				}
 			    }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case ELEC_FORCE_TANGENT1:
 			  for ( a=0; a<VIM; a++)
 			    {
 			      for ( b=0; b<VIM; b++)
 				{
                                 local_q += (fv->stangent[0][a]*perm*es[a][b]*fv->snormal[b]);
 				}
 			    }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case ELEC_FORCE_TANGENT2:
 			if(pd->Num_Dim == 3)
 			  {
 			  for ( a=0; a<VIM; a++)
 			    {
 			     for ( b=0; b<VIM; b++)
 				{
                                 local_q += (fv->stangent[1][a]*perm*es[a][b]*fv->snormal[b]);
 				}
 			    }
 			   }
 			  else
 			   {
 			      EH(-1, "Illegal flux type");
 			   }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case ELEC_FORCE_X:
 			  for ( a=0; a<VIM; a++)
 			    {
                             local_q += (perm*es[0][a]*fv->snormal[a]);
 			    }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case ELEC_FORCE_Y:
 			  for ( a=0; a<VIM; a++)
 			    {
                                   local_q += (perm*es[1][a]*fv->snormal[a]);
 			    }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case ELEC_FORCE_Z:
 			if(pd->Num_Dim == 3)
 			  {
 			  for ( a=0; a<VIM; a++)
 			    {
                                   local_q += (perm*es[2][a]*fv->snormal[a]);
 			    }
 			  }
 			  else
 			   {
 			      EH(-1, "Illegal flux type");
 			   }
                                   local_flux += weight * det * local_q;
 		      break;
 		    case NET_SURF_CHARGE:
 			  for ( a=0; a<VIM; a++)
 			    {
                                   local_q += (-perm * fv->snormal[a] * efield[a]);
 			    }
                                   local_flux += weight * det * local_q;
 		      break;

		    case TORQUE:
		      if(pd->CoordinateSystem == PROJECTED_CARTESIAN)
			EH(-1, "TORQUE has not been updated for the PROJECTED_CARTESIAN coordinate system.");

		      if(pd->CoordinateSystem == SWIRLING || 
			 pd->CoordinateSystem == CYLINDRICAL)
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
			      /*
			       *  note that for CYLINDRICAL and SWIRLING coordinate systems
			       * the h3 factor has been incorporated already into sdet
			       * the moment arm is incorporated into sideset.
			       */

                                  local_q += ( fv->x[1] * e_theta[a]*
					(vs[a][b]+ves[a][b])*fv->snormal[b] );

				}
			    }
                                  local_flux += weight * det* local_q;
			}
		      else if( pd->CoordinateSystem == CARTESIAN ) 
			{
			  if( VIM == 2 ) 
			    {
			      fv->x[2] = 0.0;
			      fv->snormal[2] = 0.0;
			    }
			  for ( a=0; a<DIM; a++)
			    {
			      Tract[a] = 0.0;
			      for ( b=0; b<DIM; b++)
				{
				  Tract[a] += (vs[a][b]+ves[a][b]) * fv->snormal[b];
				}
			    }
			  for ( a=0; a<DIM; a++)
			    {
			      for ( b=0; b<DIM; b++)
				{
				  for ( c=0; c<DIM; c++)
				    {
				      local_Torque[a] += ( permute(b,c,a) *
						     fv->x[b]*Tract[c] );
				    }
				}
			    }
			  for ( a=0; a<DIM; a++)
			    { Torque[a] += weight * det * local_Torque[a];}
			  local_flux = Torque[2];
			}
		      else
			{
			  EH(-1,"Torque cannot be calculated in this case.");
			}
		      break;
		  
		    case PORE_LIQ_FLUX:
				
				err = load_porous_properties();
				EH( err, "load_porous_properties");

		      if(mp->PorousMediaType == POROUS_SATURATED)
			{

			  for ( a=0; a<VIM; a++)
			    {
			      local_flux += weight*det*mp->density * 
				pmv->liq_darcy_velocity[a]*fv->snormal[a];
			      local_flux_conv += 0.;
			    }
			}
		      else if (mp->PorousMediaType == POROUS_UNSATURATED ||
			       mp->PorousMediaType == POROUS_TWO_PHASE)
			{

			  for ( a=0; a<VIM; a++)
			    {
			      local_flux += weight*det* 
				pmv->rel_mass_flux[0][a]*fv->snormal[a];
			      local_flux_conv += 0.;
			    }
			}
		      else
			{
			  EH(-1,"unrecognized porous media type in mm_flux.c");
			}
		      break;

		    case FORCE_NORMAL:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
                                  local_q += (fv->snormal[a]*
                        		(vs[a][b]+ves[a][b] )*fv->snormal[b]);

                                  local_qconv += ( -rho*fv->snormal[a]*
                                (fv->v[a]-x_dot[a])*fv->v[b]*fv->snormal[b]);
				}
			    }
                                  local_flux += weight * det * local_q;
                                  local_flux_conv += weight*det*local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
                                  local_q += (fv->snormal[a]*TT[a][b]*fv->snormal[b]);
				}
			    }
                                  local_flux += weight * det * local_q;
			}
		      break;
		      
		    case FORCE_TANGENT1:
		      
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
                                  local_q += (fv->stangent[0][a]*
					(vs[a][b] + ves[a][b])*fv->snormal[b]);
                                  local_qconv += ( - rho*fv->stangent[0][a]
                                *(fv->v[a]-x_dot[a]) *fv->v[b]*fv->snormal[b]);
				}
			    }
                                  local_flux += weight * det * local_q;
                          local_flux_conv += weight * det*local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
                                  local_q += ( fv->stangent[0][a]*TT[a][b]
							*fv->snormal[b]);
				}
			    }
                                  local_flux += weight * det*local_q;
			}
		      break;
		      
		    case FORCE_TANGENT2:
		
		      if(cr->MeshMotion == ARBITRARY)
			{
			  if(pd->Num_Dim == 3)
			    {
			      for ( a=0; a<VIM; a++)
				{
				  for ( b=0; b<VIM; b++)
				    {
                                      local_q += (fv->stangent[1][a]
						*(vs[a][b] + ves[a][b])
							*fv->snormal[b]);
                                      local_qconv += ( - rho*fv->stangent[1][a]
	                        *(fv->v[a]-x_dot[a]) *fv->v[b]*fv->snormal[b]);
				    }
				}
                                  local_flux += weight * det * local_q;
                                  local_flux_conv += weight * det * local_qconv;
			    }
			  else
			    {
			      EH(-1, "Illegal flux type");
			    }
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
                                  local_q += ( fv->stangent[1][a]*TT[a][b]
							*fv->snormal[b]);
				}
			    }
                                  local_flux += weight * det * local_q;
			}
		      break;
		      
		    case FORCE_X:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
			      local_q += ((vs[0][a] + ves[0][a] ) *fv->snormal[a]) ;
			      local_qconv += ( -rho*(fv->v[0]-x_dot[0])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
			  local_flux += weight * det * local_q;
			  local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[0][a]*fv->snormal[a]);
			    }
                              local_flux += weight * det * local_q;
			}
		      break;
		      
		    case FORCE_X_POS:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ((vs[0][a] + ves[0][a] ) *fv->snormal[a]) ;
                              local_qconv += ( -rho*(fv->v[0]-x_dot[0])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
        	              if(local_q <0) {local_q=0;}
                              local_flux += weight * det * local_q;
                              local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[0][a]*fv->snormal[a]);
			    }
        	              if(local_q <0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}
		      break;
		      
		    case FORCE_X_NEG:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ((vs[0][a] + ves[0][a] ) *fv->snormal[a]) ;
                              local_qconv += ( -rho*(fv->v[0]-x_dot[0])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
        	              if(local_q >0) {local_q=0;}
                              local_flux += weight * det * local_q;
                              local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[0][a]*fv->snormal[a]);
			    }
        	              if(local_q >0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}
		      break;
		      
		    case FORCE_Y:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
                             local_q += ((vs[1][a] + ves[1][a]) *fv->snormal[a]);
                             local_qconv += ( -rho*(fv->v[1]-x_dot[1])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
                              local_flux += weight * det * local_q;
                              local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[1][a]*fv->snormal[a]);
			    }
                              local_flux += weight * det * local_q;
			}
		      break;

		    case FORCE_Y_POS:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
                             local_q += ((vs[1][a] + ves[1][a]) *fv->snormal[a]);
                             local_qconv += ( -rho*(fv->v[1]-x_dot[1])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
        	              if(local_q <0) {local_q=0;}
                              local_flux += weight * det * local_q;
                              local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[1][a]*fv->snormal[a]);
			    }
        	              if(local_q <0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}
		      break;

		    case FORCE_Y_NEG:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for ( a=0; a<VIM; a++)
			    {
                             local_q += ((vs[1][a] + ves[1][a]) *fv->snormal[a]);
                             local_qconv += ( -rho*(fv->v[1]-x_dot[1])
                                                   *fv->v[a]*fv->snormal[a]);
			    }
        	              if(local_q >0) {local_q=0;}
                              local_flux += weight * det * local_q;
                              local_flux_conv += weight * det * local_qconv;
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[1][a]*fv->snormal[a]);
			    }
        	              if(local_q >0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}
		      break;

		    case FORCE_Z:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  if(pd->Num_Dim == 3)
			    {
			      for ( a=0; a<VIM; a++)
				{
                                  local_q += ((vs[2][a] + ves[2][a]) *fv->snormal[a]);

                                  local_qconv += ( -rho*(fv->v[2]-x_dot[2])
                                                   *fv->v[a]*fv->snormal[a]);

				}
                                  local_flux += weight * det * local_q;
                                  local_flux_conv += weight *det*local_qconv;
			    }
			  else
			    {
			      EH(-1, "Illegal flux type");
			    }
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[2][a]*fv->snormal[a]);
			    }
                              local_flux += weight * det * local_q;
			}

		      break;

		    case FORCE_Z_POS:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  if(pd->Num_Dim == 3)
			    {
			      for ( a=0; a<VIM; a++)
				{
                                  local_q += ((vs[2][a] + ves[2][a]) *fv->snormal[a]);

                                  local_qconv += ( -rho*(fv->v[2]-x_dot[2])
                                                   *fv->v[a]*fv->snormal[a]);

				}
        	                  if(local_q <0) {local_q=0;}
                                  local_flux += weight * det * local_q;
                                  local_flux_conv += weight *det*local_qconv;
			    }
			  else
			    {
			      EH(-1, "Illegal flux type");
			    }
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[2][a]*fv->snormal[a]);
			    }
        	              if(local_q <0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}

		      break;

		    case FORCE_Z_NEG:
		      if(cr->MeshMotion == ARBITRARY)
			{
			  if(pd->Num_Dim == 3)
			    {
			      for ( a=0; a<VIM; a++)
				{
                                  local_q += ((vs[2][a] + ves[2][a]) *fv->snormal[a]);

                                  local_qconv += ( -rho*(fv->v[2]-x_dot[2])
                                                   *fv->v[a]*fv->snormal[a]);

				}
        	                  if(local_q >0) {local_q=0;}
                                  local_flux += weight * det * local_q;
                                  local_flux_conv += weight *det*local_qconv;
			    }
			  else
			    {
			      EH(-1, "Illegal flux type");
			    }
			}
		      else
			{
			  for ( a=0; a<VIM; a++)
			    {
                              local_q += ( TT[2][a]*fv->snormal[a]);
			    }
        	              if(local_q >0) {local_q=0;}
                              local_flux += weight * det * local_q;
			}

		      break;

		    case AVERAGE_CONC:

                      local_q += fv->c[species_id];
                      local_flux += weight*det*local_q;

		      break;

		    case REPULSIVE_FORCE:
                      for (a = 0; a < Num_BC; a++) {
                         if( BC_Types[a].BC_ID != side_set_id ||
                             BC_Types[a].BC_Name < CAP_REPULSE_BC ||
                             BC_Types[a].BC_Name > CAP_REPULSE_TABLE_BC)
                            {continue;}
                         else if(BC_Types[a].BC_Name == CAP_REPULSE_ROLL_BC )
                            {
		              double roll_rad, origin[3],dir_angle[3];
		              double hscale, repexp, P_rep, betainv;
                              double factor,t,axis_pt[3],R,kernel,tangent,dist;
                              double coord[3]={0,0,0};

	                      roll_rad = BC_Types[a].BC_Data_Float[1];
			      origin[0] =  BC_Types[a].BC_Data_Float[2];
			      origin[1] =  BC_Types[a].BC_Data_Float[3];
			      origin[2] =  BC_Types[a].BC_Data_Float[4];
			      dir_angle[0] =  BC_Types[a].BC_Data_Float[5];
			      dir_angle[1] =  BC_Types[a].BC_Data_Float[6];
			      dir_angle[2] =  BC_Types[a].BC_Data_Float[7];
	                      hscale = BC_Types[a].BC_Data_Float[8];
	                      repexp = BC_Types[a].BC_Data_Float[9];
	                      P_rep = BC_Types[a].BC_Data_Float[10];
	                      betainv = BC_Types[a].BC_Data_Float[11];
/*  initialize variables */

			  for( b=0 ; b<pd->Num_Dim ; b++)
			    {  coord[b] = fv->x[b];  }
/*  find intersection of axis with normal plane - i.e., locate point on
 *          axis that intersects plane normal to axis that contains local point. */
    factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
    t = (dir_angle[0]*(coord[0]-origin[0]) + dir_angle[1]*(coord[1]-origin[1])
        + dir_angle[2]*(coord[2]-origin[2]))/factor;
    axis_pt[0] = origin[0]+dir_angle[0]*t;
    axis_pt[1] = origin[1]+dir_angle[1]*t;
    axis_pt[2] = origin[2]+dir_angle[2]*t;

/*  compute radius and radial direction */

    R = sqrt( SQUARE(coord[0]-axis_pt[0]) + SQUARE(coord[1]-axis_pt[1]) +
                SQUARE(coord[2]-axis_pt[2]) );
    dist = R - roll_rad;

/*  repulsion function  */
                              kernel = P_rep/pow(dist/hscale, repexp); 
                              tangent = betainv/pow(dist/hscale, repexp); 
                              local_q = kernel;
                              local_qconv = tangent;
			      for( b=0 ; b<pd->Num_Dim ; b++)
                                 {local_Torque[b] = fv->x[b]*kernel;}
                            }
                         else if(BC_Types[a].BC_Name == CAP_REPULSE_TABLE_BC )
                            {
		              double hscale, repexp, P_rep, betainv, exp_scale;
                              double mod_factor,kernel,dist=0,tangent;
                              double dcl_dist, slope, d_tfcn[3];
                              int bc_table_id=-1,dcl_node,k,nsp,i1,i2;
                              double point[3]={0,0,0};
                              double coord[3]={0,0,0};

	                      hscale = BC_Types[a].BC_Data_Float[0];
	                      repexp = BC_Types[a].BC_Data_Float[1];
	                      P_rep = BC_Types[a].BC_Data_Float[2];
			      betainv =  BC_Types[a].BC_Data_Float[3];
			      exp_scale =  BC_Types[a].BC_Data_Float[4];
			      dcl_node =  BC_Types[a].BC_Data_Int[2];
/*  initialize variables */

			  for( b=0 ; b<pd->Num_Dim ; b++)
			    {  coord[b] = fv->x[b];  }

                          for (b = 0; b < Num_BC; b++) {
                              if( BC_Types[b].BC_Name == GD_TABLE_BC )
                                  { bc_table_id = b; }
                              }
                          if(bc_table_id != -1)
      dist = table_distance_search(BC_Types[bc_table_id].table, coord, 
                                                         &slope, d_tfcn);
                          if(dcl_node != -1)
                            {
                             nsp = match_nsid(dcl_node);
                             k = Proc_NS_List[Proc_NS_Pointers[nsp]];
                             for (b = 0; b < Proc_NS_Count[nsp]; b++)
                               {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+b];
               i1 = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1);
               EH(i1, "Could not resolve index_solution.");
                                for(i2=0 ; i2<pd->Num_Dim ; i2++)
                                  {
                                   point[i2] = Coor[i2][k] + x[i1+i2];
                                  }
                               }
                             dcl_dist = sqrt(SQUARE(coord[0]-point[0])
                                        +SQUARE(coord[1]-point[1])
                                        +SQUARE(coord[2]-point[2]));
                             }  else	{ dcl_dist = 0.; }
/*  modifying function for DCL  */
                           mod_factor = 1. - exp(-dcl_dist/exp_scale);
/*  repulsion function  */
                           kernel = P_rep*mod_factor/pow(dist/hscale, repexp); 
                           tangent = betainv*mod_factor/pow(dist/hscale,repexp); 
/*  repulsion function  */
                           local_q = kernel;
                           local_qconv = tangent;
			   for( b=0 ; b<pd->Num_Dim ; b++)
                                 {local_Torque[b] = fv->x[b]*kernel;}
                           if(pd->Num_Dim < DIM)
                                 {local_Torque[pd->Num_Dim] = dist*kernel;}
                            }
                         else if(BC_Types[a].BC_Name == CAP_REPULSE_BC )
                            {
		              double repexp, P_rep, factor, denom;
                              double kernel,dist=0;
                              double ap, bp, cp, dp;

	                      repexp = 4.;
	                      P_rep = BC_Types[a].BC_Data_Float[0];
			      ap =  BC_Types[a].BC_Data_Float[1];
			      bp =  BC_Types[a].BC_Data_Float[2];
			      cp =  BC_Types[a].BC_Data_Float[3];
			      dp =  BC_Types[a].BC_Data_Float[4];
                              denom = sqrt(ap*ap + bp*bp + cp*cp);
                              factor = ap*fv->x[0] + bp*fv->x[1] + cp*fv->x[2] +dp;
                              dist = fabs(factor)/denom;
/*  initialize variables */
/*  repulsion function  */
                           kernel = P_rep/pow(dist, repexp); 
/*  repulsion function  */
                           local_q = kernel;
			   for( b=0 ; b<pd->Num_Dim ; b++)
                                 {local_Torque[b] = fv->x[b]*kernel;}
                            }
                         else 
                            {EH(-1,"Repulsive force not found\n");} 
                         }

			      local_flux += weight * det * local_q;
                              local_flux_conv += weight *det*local_qconv;
			      for( b=0 ; b<DIM ; b++)
                                 {Torque[b] += local_Torque[b]*weight*det;}
		      break;

		    case SURF_DISSIP:
		      /* This is the energy dissipated at the surface due to surface tension 
		       *  See Batchelor, JFM, 1970 for details 
		       */
		      for( a=0; a<VIM ; a++)
			{
			  for( b=0 ; b<VIM ; b++)
			    {
			      local_q += mp->surface_tension * 
                                       ( fv->grad_v[a][b]*( delta(a,b)
					 -fv->snormal[a]*fv->snormal[b] ) );
			    }
			}

			      local_flux += weight * det * local_q;
                              local_flux_conv += weight *det*local_qconv;
		      break;

		    case POYNTING_X:
		    case POYNTING_Y:
		    case POYNTING_Z:
			/* For scalar e-field calculations, we will assume the e-vector 
                         * points out of the plane, i.e. normal to the plane of 
                         * incidence, Ez */
		      R_imped = acoustic_impedance( d_R, time_value );
		      wnum = wave_number( d_wnum, time_value );
		      kR_inv = 1./(wnum*R_imped);
                      memset( Mag_real,0, sizeof(double)*DIM);
                      memset( Mag_imag,0, sizeof(double)*DIM);
                      memset( E_real,0, sizeof(double)*DIM);
                      memset( E_imag,0, sizeof(double)*DIM);
		      if(pd->CoordinateSystem == PROJECTED_CARTESIAN)
			EH(-1, "POYNTING has not been updated for the PROJECTED_CARTESIAN coordinate system.");

		      if(pd->CoordinateSystem == SWIRLING || 
			 pd->CoordinateSystem == CYLINDRICAL)
			{
			EH(-1, "POYNTING has not been checked for CYLINDRICAL yet.");
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
			      /*
			       *  note that for CYLINDRICAL and SWIRLING coordinate systems
			       * the h3 factor has been incorporated already into sdet
			       * the moment arm is incorporated into sideset.
			       */

                                  local_q += ( fv->x[1] * e_theta[a]*
					(vs[a][b]+ves[a][b])*fv->snormal[b] );

				}
			    }
                                  local_flux += weight * det* local_q;
			}
		      else if( pd->CoordinateSystem == CARTESIAN ) 
			{
                          Mag_imag[0] = kR_inv*fv->grad_apr[1];
                          Mag_imag[1] = -kR_inv*fv->grad_apr[0];
                          Mag_real[0] = kR_inv*fv->grad_api[1];
                          Mag_real[1] = -kR_inv*fv->grad_api[0];
                          E_real[2] = fv->apr;  E_imag[2] = fv->api;
			  for ( a=0; a<DIM; a++)
			    {
			      for ( b=0; b<DIM; b++)
				{
				  for ( c=0; c<DIM; c++)
				    {
				      local_Torque[a] += 0.5*( permute(b,c,a) *
				(E_real[b]*Mag_real[c]-E_imag[b]*Mag_imag[c]) );
				    }
				}
			    }
			  for ( a=0; a<DIM; a++)
			    { Torque[a] += weight * det * local_Torque[a];}
			  local_flux = Torque[quantity-POYNTING_X];
			}
		      else
			{
			  EH(-1,"Torque cannot be calculated in this case.");
			}
		      break;
		  
		    case N_DOT_X:
		      /* 
		       * This is the position vector dotted into the local normal
		       * and has application in determining volumes of shapes 
		       * via surface integrals (Ask me how ! - tab )
		       */

		      for( a=0 ; a<dim; a++ )
			{
			  local_q += ( fv->x[a] - param[a] ) *fv->snormal[a];
			}

		      local_flux += weight*det*local_q/ (double) dim;

						
				break;
				
			case DELTA:
				
				load_lsi(ls->Length_Scale );
				local_q = lsi->delta;
	 			local_qconv = lsi->delta;
				local_flux += weight * det * lsi->delta;
  	  			local_flux_conv += ad_wt[ip_total-ip-1]*det*lsi->delta;
				break;

			case LS_DCA:
				load_lsi(ls->Length_Scale );
				ierr = interface_crossing_1DQ(ls_F, xf);
				if(ierr  && ip == 0)
					{
					 switch	(id_side)
			 			{
						 case 1:
							xi[0] = xf[0];
							xi[1] = -1;
							break;
						 case 2:
							xi[0] = 1;
							xi[1] = xf[0];
							break;
						 case 3:
					 		xi[0] = xf[0];
					 		xi[1] = 1;
							break;
						 case 4:
					 		xi[0] = -1;
					 		xi[1] = xf[0];
							break;
	  					 }
				 	xi[2] = 0;
       
			err = load_basis_functions( xi, bfd );
 			EH( err, "problem from load_basis_functions");
 
 			err = beer_belly();
 			EH( err, "beer_belly");
 	      
 			err = load_fv();
 			EH( err, "load_fv");
 
			err = load_bf_grad();
			EH( err, "load_bf_grad");
		  
			err = load_bf_mesh_derivs();
			EH( err, "load_bf_mesh_derivs");
			  
		surface_determinant_and_normal (ei->ielem, iconnect_ptr, 
		num_local_nodes, ielem_dim - 1,  id_side, num_nodes_on_side, 
		id_local_elem_coord );
			  
	 /*
	  * Load up physical space gradients of field variables at this
	  * Gauss point.
	  */
	 		err = load_fv_grads();
	 		EH( err, "load_fv_grads");
			  
	 		err = load_fv_mesh_derivs(1);
	 		EH( err, "load_fv_mesh_derivs");

			if (TimeIntegration != STEADY && 
					pd->e[MESH_DISPLACEMENT1]) {
				for (j = 0; j < VIM; j++) {
					x_dot[j] = fv_dot->x[j];
				  	} 
				} else {
				for (j = 0; j < VIM; j++) {
					x_dot[j] = 0.;
				  	}
				}
			  
	 		do_LSA_mods(LSA_SURFACE);

	 		if (ielem_dim !=3)
	      		{
			calc_surf_tangent (ei->ielem, iconnect_ptr, 
				   num_local_nodes, ielem_dim-1,
				   num_nodes_on_side,
				   id_local_elem_coord);
	      		}
			load_lsi(ls->Length_Scale );
                       if(profile_flag & 1)
                       {
			gamma_dot = 0.;
   			for ( a=0; a<VIM; a++)
     			   {
       				for ( b=0; b<VIM; b++)
       				  {
       				   gamma_dot +=  gamma[a][b] * gamma[b][a];
       				  }
     			   }
   			gamma_dot  =  sqrt(gamma_dot/2.);
        if(ls->SubElemIntegration)
        	{
        	elem_sign_org = ls->Elem_Sign;
        	switch (mp->mp2nd->viscositymask[1]-mp->mp2nd->viscositymask[0])
           		{
            		case 1:
                		ls->Elem_Sign = -1;
                		break;
            		case -1:
                		ls->Elem_Sign = 1;
                		break;
           		}
          	mu = viscosity(gn, gamma, d_mu);
        	ls->Elem_Sign = elem_sign_org;
        	}
        	else
        	{
          	mu = viscosity(gn, gamma, d_mu);
        	switch (mp->mp2nd->viscositymask[1]-mp->mp2nd->viscositymask[0])
           		{
            		case 1:
               		mu = (mu-mp->mp2nd->viscosity*lsi->H)/(1.-lsi->H);
                		break;
            		case -1:
               		mu = (mu-mp->mp2nd->viscosity*(1.-lsi->H))/lsi->H;
                		break;
           		}
        	}
                       }

         		for (a=0; a<VIM; a++)
             			{
				local_q -= fv->snormal[a]*lsi->normal[a];
				local_qconv += SQUARE(fv->v[a] - x_dot[a]);
	     			}
			local_q = acos(local_q)*180/M_PIE;
			local_qconv = sqrt(local_qconv);
			local_flux = local_q;
			local_flux_conv = local_qconv;
			}
				break;

		    default:
		      
		      EH(-1, "Illegal flux type");
		      break;
		    }  /*  end of switch */
	    
		  local_area += weight * det ;
#ifdef PARALLEL
        delta_flux = local_flux - local_flux0;
        delta_flux_conv = local_flux_conv - local_flux_conv0;
        delta_area = local_area - local_area0;

        if( Num_Proc > 1 &&  dpi->elem_owner[ elem_list[i] ] == ProcID )
         {
          proc_flux += delta_flux;
          proc_flux_conv += delta_flux_conv;
          proc_area += delta_area;
         }

        local_flux0 = local_flux;
        local_flux_conv0 = local_flux_conv;
        local_area0 = local_area;
#endif

        if (profile_flag && print_flag && 
		(quantity != LS_DCA || (ierr && ip == 0)) ) {
        FILE  *jfp;
        if( (jfp=fopen(filenm,"a")) != NULL)
            {
             fprintf(jfp," %g  %g  %g  %g  %g",
                        fv->x[0],fv->x[1],fv->x[2],local_q,local_qconv);
        if(profile_flag & 1) fprintf(jfp," %g  %g ", gamma_dot, mu);
        if(profile_flag & 4) fprintf(jfp," %g ", mp->surface_tension);
        if(profile_flag & 8)  fprintf(jfp," %g %g",fv->T, fv->P);
        if(profile_flag & 2) fprintf(jfp," %g %g %g",
                         evalue1,evalue2,evalue3);
        if(profile_flag & 16) fprintf(jfp," %g %g %g",
                         fv->snormal[0],fv->snormal[1],fv->snormal[2]);
              fprintf(jfp," \n");
              fflush(jfp);
            }
          fclose(jfp);
                }
  
       /* Compute sensitivities if requested */
	if ( J_AC != NULL)
		    {
		      int dir, sp = species_id;
		      double d_term=0, d_term1=0, d_term2=0, d_term3=0;
		      double (* d_diff)[MAX_VARIABLE_TYPES+MAX_CONC]  = mp->d_diffusivity;

		      switch (quantity)
			{
	 /* FORCE_X, FORCE_Y, FORCE_Z */
	 /* FORCE_X_POS, FORCE_Y_POS, FORCE_Z_POS */
	 /* FORCE_X_NEG, FORCE_Y_NEG, FORCE_Z_NEG */
			case FORCE_X:
			case FORCE_Y:
			case FORCE_Z:
			case FORCE_X_POS:
			case FORCE_Y_POS:
			case FORCE_Z_POS:
			case FORCE_X_NEG:
			case FORCE_Y_NEG:
			case FORCE_Z_NEG:
			   
			  dir = ( quantity == FORCE_X ? 0 :
				( quantity == FORCE_Y ? 1 : 
			        ( quantity == FORCE_Z ? 2 : 
			        ( quantity == FORCE_X_POS ? 0 :
				( quantity == FORCE_Y_POS ? 1 : 
			        ( quantity == FORCE_Z_POS ? 2 : 
			        ( quantity == FORCE_X_NEG ? 0 :
				( quantity == FORCE_Y_NEG ? 1 : 
			        ( quantity == FORCE_Z_NEG ? 2 : -1))))))))) ;  /* Can't do this with FORTRAN */

			  EH(dir, "PANIC in evaluate_flux sensitivity section");

		      if(cr->MeshMotion == ARBITRARY)
			{
			  for( b=0 ; b < dim ; b++)
			    {
			      var = VELOCITY1 + b ;

			      if( pd->v[var] )
				{
				  for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = 0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term += weight * det * fv->snormal[a] *
					    ( mu * ( bf[var]->grad_phi_e[j][b][dir][a] +
						     bf[var]->grad_phi_e[j][b][a][dir] ) + 
					      d_mu->v[b][j] * ( gamma[dir][a] )  +
					      -rho * ( bf[var]->phi[j]*delta(dir,b) * fv->v[a] +
						       bf[var]->phi[j]*delta(a,b) *( fv->v[dir] - x_dot[dir] ) ) );
					}
				      J_AC[ ei->gun_list[var][j] ] += d_term;
				    }
				}
			    }

			  var = PRESSURE ;
			  
			  if( pd->v[var] )
			    {
			      for(j=0 ; j < ei->dof[var]; j++)
				{				  
				  d_term = -weight * det * fv->snormal[dir] * bf[var]->phi[j] ;
				  J_AC[ ei->gun_list[var][j] ] += d_term;
				}
			    }


			  for( b=0 ; b < dim ; b++)
			    {
			      var = MESH_DISPLACEMENT1 + b ;

			      if( pd->v[var] )
				{
				  for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 = 0.0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term1 += weight*(vs[dir][a]+ves[dir][a])*fv->snormal[a] *
					             fv->dsurfdet_dx[b][j] ; 

				          d_term2 += weight*det*( vs[dir][a]+ves[a][b] ) * 
					             fv->dsnormal_dx[a][b][j];

					  d_term  += weight*det*fv->snormal[a] *
					             ( mu*( fv->d_grad_v_dmesh[dir][a][b][j] + fv->d_grad_v_dmesh[a][dir][b][j] ) +
						       d_mu->X[b][j] * gamma[dir][a] );
					}

				      for( a=0 ; a < dim ; a++ )
					{
					  d_term1 += -weight*rho*(fv->v[dir] - x_dot[dir])*fv->v[a] * fv->snormal[a] *
					              fv->dsurfdet_dx[b][j];

					  d_term2 += -weight*det*rho*(fv->v[dir] - x_dot[dir])*fv->v[a] *
					              fv->dsnormal_dx[a][b][j];
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1 + d_term2;

				    }
				}
			    }
			  
			  var = TEMPERATURE ;
			  
			  if( pd->v[var] )
			    {
			      for(j=0 ; j < ei->dof[var]; j++)
				{			
				  for( d_term=0., a=0; a < dim ; a++ )
				    {
				      d_term += weight *det * fv->snormal[a] * ( d_mu->T[j] * gamma[dir][a] ) ;
				    }
				  J_AC[ei->gun_list[var][j]] += d_term;
				}
			    }

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0; j < ei->dof[var]; j++) {
				/*
				 * Find the material index and gnn for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim ; a++ ) {
				  d_term += weight * det * fv->snormal[a] *
				      (d_mu->C[w][j] * gamma[dir][a]);
				}
				J_AC[c] += d_term;
			      }
			    }
			  }

  			if ( pd->v[POLYMER_STRESS11] )
    			  {
			    for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
                	       {
			         for ( b=0; b<VIM; b++)
            			   {
				      var = v_s[ve_mode][dir][b];
				      for(j=0 ; j < ei->dof[var]; j++)
				    	    {
				      d_term = weight *det * fv->snormal[b] * ( bf[var]->phi[j] ) ;
				      J_AC[ ei->gun_list[var][j] ] += d_term;
					    }
                		    }
        		      }
			  }
			}
			else	/*  Solid Stress terms	*/
			{
			  for( b=0 ; b < dim ; b++)
			    {
			      var = MESH_DISPLACEMENT1 + b ;

			      if( pd->v[var] )
				{
				  for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 = 0.0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term1 += weight*TT[dir][a]*fv->snormal[a] *
					             fv->dsurfdet_dx[b][j] ; 

				          d_term2 += weight*det*TT[dir][a]* 
					             fv->dsnormal_dx[a][b][j];

					  d_term  += weight*det*fv->snormal[a] *
					              dTT_dx[dir][a][b][j];
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1 + d_term2;

				    }
				}
			    }
			 var = PRESSURE;
			 if( pd->v[var] )
				{
				 for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term  += weight*det*fv->snormal[a] *
					              dTT_dp[dir][a][j];
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;

				    }
				}
			 var = TEMPERATURE;
			 if( pd->v[var] )
				{
				 for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term  += weight*det*fv->snormal[a] *
					              dTT_dT[dir][a][j];
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;

				    }
				}
			  for( b=0 ; b < dim ; b++)
			    {
			      var = SOLID_DISPLACEMENT1 + b ;

			      if( pd->v[var] )
				{
				  for(j=0 ; j < ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < dim ; a++) 
					{
					  d_term  += weight*det*fv->snormal[a] *
					              dTT_drs[dir][a][b][j];
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;

				    }
				}
			    }
			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0; j < ei->dof[var]; j++) {
				/*
				 * Find the material index and gnn for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim ; a++ ) {
				  d_term += weight * det * fv->snormal[a] *
				      (dTT_dc[dir][a][w][j] * gamma[dir][a]);
				}
				J_AC[c] += d_term;
			      }
			    }
			  }
			}

			  
			  break;

	  /* FORCE_NORMAL */
			case FORCE_NORMAL:
			  
		      if(cr->MeshMotion == ARBITRARY)
			{
			  for( p=0 ; p<dim; p++)
			    {
			      var = VELOCITY1 + p;
			      
			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for ( a=0; a<VIM; a++)
					{
					  for ( b=0; b<VIM; b++)
					    {
					      d_term  += weight*det *fv->snormal[a]*fv->snormal[b]*
						        ( mu*( bf[var]->grad_phi_e[j][p][a][b] + bf[var]->grad_phi_e[j][p][b][a] ) +
							  d_mu->v[p][j]*gamma[a][b] );
					      d_term1 += -weight*det* rho*fv->snormal[a]*fv->snormal[b]*
						          ( fv->v[b]*bf[var]->phi[j]*delta(p,a) +
							  ( fv->v[a] - x_dot[a])*bf[var]->phi[j]*delta(p,b) );

					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1;
				    }
				}
			    }

			  var = PRESSURE;
			      
			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;

				  for( a=0; a<dim; a++)
				    {
				      d_term += -weight*det*fv->snormal[a]*fv->snormal[a]*bf[var]->phi[j];
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term + d_term1;
				}
			    }

			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 = d_term3 = 0.0;

				      for( a=0; a < VIM; a++)
					{
					  for(b=0 ; b<VIM; b++)
					    {
					      d_term += weight * det *fv->snormal[a]*fv->snormal[b]*
						        (mu* ( fv->d_grad_v_dmesh[a][b][p][j] + fv->d_grad_v_dmesh[b][a][p][j]) + 
							 d_mu->X[p][j]* gamma[a][b] );
					      d_term1 += weight * fv->snormal[a]*(vs[a][b]+ves[a][b] )*fv->snormal[b] *
						         fv->dsurfdet_dx[p][j];
					      d_term2 += weight*det*(vs[a][b]+ves[a][b] )*
						         ( fv->snormal[a]*fv->dsnormal_dx[b][p][j] +
							   fv->dsnormal_dx[a][p][j]*fv->snormal[b] );
					      d_term3 += -weight* rho*(fv->v[a]-x_dot[a])*fv->v[b] *
						         ( fv->dsurfdet_dx[p][j]*fv->snormal[a]*fv->snormal[b] + 
							   det * fv->snormal[a]*fv->dsnormal_dx[b][p][j] +
							   det * fv->dsnormal_dx[a][p][j]*fv->snormal[b] );
 
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1 + d_term2 + d_term3;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;
				  
				  for( a=0; a < dim; a++)
				    {
				      for(b=0 ; b<dim; b++)
					{
					  d_term += weight*det*fv->snormal[a]*fv->snormal[b] * ( d_mu->T[j] * gamma[a][b] );
					}
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term;
				}
			    }
			

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0 ; j < ei->dof[var]; j++) {
				/*
				 * Find the material index for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim; a++) {
				  for (b = 0 ; b < dim; b++) {
				    d_term += weight * det * fv->snormal[a] *
					fv->snormal[b] * (d_mu->C[w][j] * gamma[a][b]);
				  }
				}
				J_AC[c] += d_term;
			      }
			    }
			  }

  			if ( pd->v[POLYMER_STRESS11] )
    			  {
			    for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
                	       {
			         for ( b=0; b<VIM; b++)
            			   {
				     for( a=0; a < dim; a++)
				       {
				      var = v_s[ve_mode][a][b];
				      for(j=0 ; j < ei->dof[var]; j++)
				    	    {
				      d_term = weight *det * fv->snormal[b] * ( bf[var]->phi[j] )* fv->snormal[a] ;
				      J_AC[ ei->gun_list[var][j] ] += d_term;
					    }
					}
                		    }
        		      }
			  }
			}
			else	/*  Solid Stress terms	*/
			{
			  var = PRESSURE;
			      
			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;

				  for( a=0; a<dim; a++)
				    {
				      d_term += -weight*det*fv->snormal[a]*fv->snormal[a]*bf[var]->phi[j];
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term + d_term1;
				}
			    }

			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 = 0.0;

				      for( a=0; a < VIM; a++)
					{
					  for(b=0 ; b<VIM; b++)
					    {
					      d_term += weight*det*fv->snormal[a]
						*fv->snormal[b]*dTT_dx[a][b][p][j];
					      d_term1 += weight * fv->snormal[a]
						*TT[a][b]*fv->snormal[b]*fv->dsurfdet_dx[p][j];
					      d_term2 += weight*det*TT[a][b]*
						         ( fv->snormal[a]*fv->dsnormal_dx[b][p][j] + fv->dsnormal_dx[a][p][j]*fv->snormal[b] );
 
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1 + d_term2;
				    }
				}
			    }
			  for( p=0; p<dim ; p++)
			    {
			      var = SOLID_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < VIM; a++)
					{
					  for(b=0 ; b<VIM; b++)
					    {
					      d_term += weight*det*fv->snormal[a]
					*fv->snormal[b]*dTT_drs[a][b][p][j];
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;
				  
				  for( a=0; a < dim; a++)
				    {
				      for(b=0 ; b<dim; b++)
					{
					  d_term += weight*det*fv->snormal[a]*fv->snormal[b] * dTT_dT[a][b][j];
					}
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term;
				}
			    }
			

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0 ; j < ei->dof[var]; j++) {
				/*
				 * Find the material index for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim; a++) {
				  for (b = 0 ; b < dim; b++) {
				    d_term += weight * det * fv->snormal[a] *
					fv->snormal[b] * dTT_dc[a][b][w][j];
				  }
				}
				J_AC[c] += d_term;
			      }
			    }
			  }

			}

			  break;


	  /* FORCE_TANGENT1, FORCE_TANGENT2  */

			case FORCE_TANGENT1:
			case FORCE_TANGENT2:

			  dir = quantity == FORCE_TANGENT1 ? 0 : 1;

		      if(cr->MeshMotion == ARBITRARY)
			{
			  for( p=0 ; p<dim; p++)
			    {
			      var = VELOCITY1 + p;
			      
			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for ( a=0; a<VIM; a++)
					{
					  for ( b=0; b<VIM; b++)
					    {
					      d_term  += weight*det *fv->stangent[dir][a]*fv->snormal[b]*
						        ( mu*( bf[var]->grad_phi_e[j][p][a][b] + bf[var]->grad_phi_e[j][p][b][a] ) +
							  d_mu->v[p][j]*gamma[a][b] );
					      d_term1 += -weight*det *rho*fv->stangent[dir][a]*fv->snormal[b]*
						         ( fv->v[b]*bf[var]->phi[j]*delta(p,a) +
							  (fv->v[a] - x_dot[a])*bf[var]->phi[j]*delta(p,b) );

					    }
					}
				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1;
				    }
				}
			    }

			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 = d_term3 = 0.0;

				      for( a=0; a < dim; a++)
					{
					  for(b=0 ; b<dim; b++)
					    {
					      d_term += weight * det *fv->stangent[dir][a]*fv->snormal[b]*
						        (mu* ( fv->d_grad_v_dmesh[a][b][p][j] + fv->d_grad_v_dmesh[b][a][p][j]) + 
							 d_mu->X[p][j]* gamma[a][b] );
					      d_term1 += weight * fv->stangent[dir][a]*(vs[a][b]+ves[a][b] )*fv->snormal[b] *
						         fv->dsurfdet_dx[p][j];
					      d_term2 += weight*det*(vs[a][b]+ves[a][b] )*
						         ( fv->stangent[dir][a]*fv->dsnormal_dx[b][p][j] +
							   fv->dstangent_dx[dir][a][p][j]*fv->snormal[b] );
					      d_term3 += -weight* rho * (fv->v[a]-x_dot[a])*fv->v[b] *
						         ( fv->dsurfdet_dx[p][j]*fv->stangent[dir][a]*fv->snormal[b] + 
							   det * fv->stangent[dir][a]*fv->dsnormal_dx[b][p][j] +
							   det * fv->dstangent_dx[dir][a][p][j]*fv->snormal[b] );
 
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term + d_term1 + d_term2 + d_term3;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;
				  
				  for( a=0; a < dim; a++)
				    {
				      for(b=0 ; b<dim; b++)
					{
					  d_term += weight*det*fv->stangent[dir][a]*fv->snormal[b] * ( d_mu->T[j] * gamma[a][b] );
					}
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term;
				}
			    }


			  var = MASS_FRACTION;	  
			  if (pd->v[var]) {
			      for(w = 0; w < pd->Num_Species_Eqn; w++) {
				for(j = 0 ; j < ei->dof[var]; j++) {
				  /*
				   * Find the material index for the current
				   * local variable degree of freedom
				   * (can't just query gun_list for MASS_FRACTION
				   *  unknowns -> have to do a lookup)
				   */		
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  for (d_term = 0.0, a = 0; a < dim; a++) {
				    for (b = 0; b < dim; b++) {
				      d_term += weight * det * 
					  fv->stangent[dir][a] *fv->snormal[b] * 
					  (d_mu->C[w][j] * gamma[a][b]);
				    }
				  }
				  J_AC[c] += d_term;
				}
			      }
			    }

  			if ( pd->v[POLYMER_STRESS11] )
    			  {
			    for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
                	       {
			         for ( b=0; b<VIM; b++)
            			   {
				     for( a=0; a < dim; a++)
				       {
				      var = v_s[ve_mode][a][b];
				      for(j=0 ; j < ei->dof[var]; j++)
				    	    {
				      d_term = weight *det * fv->snormal[b] * ( bf[var]->phi[j] )* fv->stangent[dir][a] ;
				      J_AC[ ei->gun_list[var][j] ] += d_term;
					    }
					}
                		    }
        		      }
			  }
			}
			else	/*  Solid Stress terms	*/
			{
			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < dim; a++)
					{
					  for(b=0 ; b<dim; b++)
					    {
					      d_term += weight*det*fv->stangent[dir][a]*fv->snormal[b]* dTT_dx[a][b][p][j];
 
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;
				    }
				}
			    }
			  for( p=0; p<dim ; p++)
			    {
			      var = SOLID_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = 0.0;

				      for( a=0; a < dim; a++)
					{
					  for(b=0 ; b<dim; b++)
					    {
					      d_term += weight*det*fv->stangent[dir][a]*fv->snormal[b]* dTT_drs[a][b][p][j];
 
					    }
					}

				      J_AC[ ei->gun_list[var][j] ] += d_term;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if(pd->v[var])
			    {
			      for( j=0 ; j<ei->dof[var]; j++)
				{
				  d_term = 0.0;
				  
				  for( a=0; a < dim; a++)
				    {
				      for(b=0 ; b<dim; b++)
					{
					  d_term += weight*det*fv->stangent[dir][a]*fv->snormal[b] *dTT_dT[a][b][j];
					}
				    }
				  J_AC[ ei->gun_list[var][j] ] += d_term;
				}
			    }


			  var = MASS_FRACTION;	  
			  if (pd->v[var]) {
			      for(w = 0; w < pd->Num_Species_Eqn; w++) {
				for(j = 0 ; j < ei->dof[var]; j++) {
				  /*
				   * Find the material index for the current
				   * local variable degree of freedom
				   * (can't just query gun_list for MASS_FRACTION
				   *  unknowns -> have to do a lookup)
				   */		
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  for (d_term = 0.0, a = 0; a < dim; a++) {
				    for (b = 0; b < dim; b++) {
				      d_term += weight * det * 
					  fv->stangent[dir][a] *fv->snormal[b] * 
					  dTT_dc[a][b][w][j];
				    }
				  }
				  J_AC[c] += d_term;
				}
			      }
			    }
			}

			  break;

	 /* VOLUME_FLUX */

			case VOLUME_FLUX:

			  for( p=0 ; p<dim; p++)
			    {
			      var = VELOCITY1 + p;
			      
			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for ( a=0; a<VIM; a++)
					{
					  d_term += weight*det*fv->snormal[a] * ( bf[var]->phi[j]*delta(p,a) ) ;
					}

				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    }
				}
			    }



			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for( a=0; a < dim; a++)
					{
					  d_term += weight * det * fv->dsnormal_dx[a][p][j]*( fv->v[a] - x_dot[a] ) ;
					  d_term1 += weight * fv->dsurfdet_dx[p][j]* fv->snormal[a]*( fv->v[a] - x_dot[a] ) ;
					}

				      J_AC[ ei->gun_list[var][j]  ] += d_term + d_term1;

				    }
				}
			    }
			      

			  break;

			case HEAT_FLUX:

			  var = TEMPERATURE;
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {
			      d_term=0;
			      
			      for( a=0; a<dim; a++)
				{
				  d_term += weight*det*fv->snormal[a]*
				           ( -d_k->T[j]*fv->grad_T[a] - k * bf[var]->grad_phi[j][a] +
					     rho*  ( d_Cp->T[j]*fv->T  +
						     Cp*bf[var]->phi[j] )*(fv->v[a] - x_dot[a] ) );
				}

			      J_AC[ ei->gun_list[var][j]  ] += d_term;
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = VELOCITY1 +p;
			      
			      if( pd->v[var])
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = 0.0;
				      
				      d_term += weight*det*rho*Cp*fv->T*fv->snormal[p]* bf[var]->phi[j];
				      
				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    } 
				  /* tabaer notes that d_Cp->v has not been included */
				}
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if( pd->v[var] )
				{
				  
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for( a=0; a<dim; a++)
					{
					  d_term += -weight * ( fv->dsurfdet_dx[p][j] * k * fv->snormal[a] * fv->grad_T[a] +
							       det * d_k->X[p][j] *  fv->snormal[a] * fv->grad_T[a] +
							       det * k * fv->dsnormal_dx[a][p][j] * fv->grad_T[a] + 
							       det * k * fv->snormal[a] * fv->d_grad_T_dmesh[a][p][j] );
					  d_term1 += weight * rho * fv->T * (fv->v[a]-x_dot[a]) *
					                   ( fv->dsurfdet_dx[p][j] * Cp * fv->snormal[a] +
							     det * d_Cp->X[p][j] * fv->snormal[a] +
							     det * Cp * fv->dsnormal_dx[a][p][j] );
					}
				      
				      J_AC[ ei->gun_list[var][j]  ] += d_term + d_term1;
				    }
				}
			    }

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0; j < ei->dof[var]; j++) {
				/*
				 * Find the material index for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim; a++) {
				  d_term += - weight * det * d_k->C[w][j] *
				      fv->snormal[a] * fv->grad_T[a] +
				      weight * det * rho * d_Cp->C[w][j] *
				      fv->T*fv->snormal[a]* (fv->v[a] - x_dot[a]);
				}
				J_AC[c] += d_term;
			      }
			    }
			  }
			  break;

			case SPECIES_FLUX:

			    var = MASS_FRACTION;
			    if (pd->v[var]) {
			      for (w = 0; w < pd->Num_Species_Eqn; w++) {
				for (j = 0; j < ei->dof[var]; j++) {
				  /*
				   * Find the material index for the current
				   * local variable degree of freedom
				   * (can't just query gun_list for MASS_FRACTION
				   *  unknowns -> have to do a lookup)
				   */		
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  d_term = d_term1 = 0.0;
				  for (a = 0; a < dim; a++) {
				    d_term += - weight * det * fv->snormal[a] *
					(d_diff[sp][MAX_VARIABLE_TYPES + w] * 
					 bf[var]->phi[j] * fv->grad_c[species_id][a] +
					 mp->diffusivity[species_id] *
					 bf[var]->grad_phi[j][a] * delta(w, species_id));

				    d_term1 += weight * det * fv->snormal[a] * 
					(fv->v[a]-x_dot[a]) *
					bf[var]->phi[j] * delta(w,species_id);
				  }
				  /*J_AC[c] += d_term + d_term1;*/
				  J_AC[c] += d_term;
				}
			      }
			    }

			  for(p=0; p<dim; p++)
			    {
			      var = VELOCITY1 + p;

			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = weight*det*fv->snormal[p]*fv->c[species_id]*bf[var]->phi[j];

				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    }
				}
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;
			      
			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 =0.0;

				      for( a=0; a<VIM; a++)
					{
					  d_term += -weight*(
						      fv->dsurfdet_dx[p][j]*mp->diffusivity[sp]*fv->snormal[a]*fv->grad_c[sp][a] + 
						      det* (d_diff[sp][var]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_c[sp][a] +
						            mp->diffusivity[sp]*(fv->dsnormal_dx[a][p][j]*fv->grad_c[species_id][a] +
										 fv->snormal[a]*fv->d_grad_c_dmesh[a][sp][p][j]  ) ) );
					  d_term1 += weight*(fv->v[a] - x_dot[a])*fv->c[sp] *
					                    (  fv->dsurfdet_dx[p][j] * fv->snormal[a] +
							       det * fv->dsnormal_dx[a][p][j] );
					}
				      J_AC[ ei->gun_list[var][j]  ] += d_term + d_term1;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if( pd->v[var])
			    {
			      for( j=0; j<ei->dof[var]; j++)
				{
				  d_term = 0;
				  
				  for( a=0; a<VIM; a++)
				    {
				      d_term += -weight*det*d_diff[sp][var]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_c[sp][a] ;
				    }

				  J_AC[ ei->gun_list[var][j]  ] += d_term;
				}
			    }				  

			  break;

			case CHARGED_SPECIES_FLUX:

			    var = MASS_FRACTION;
			    if (pd->v[var]) {
			      for (w = 0; w < pd->Num_Species_Eqn; w++) {
				for (j = 0; j < ei->dof[var]; j++) {
				  /*
				   * Find the material index for the current
				   * local variable degree of freedom
				   * (can't just query gun_list for MASS_FRACTION
				   *  unknowns -> have to do a lookup)
				   */		
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  d_term = d_term1 = 0.0;
				  for (a = 0; a < dim; a++) {
				    d_term += - weight * det * fv->snormal[a] *
					(d_diff[sp][MAX_VARIABLE_TYPES + w] * 
					 bf[var]->phi[j] * fv->grad_c[species_id][a] +
					 mp->diffusivity[species_id] *
					 bf[var]->grad_phi[j][a] * delta(w, species_id) 
                                         + d_kapta_dc[species_id][MAX_VARIABLE_TYPES + w]*
                                         bf[var]->phi[j]*fv->grad_V[a]);

				    d_term1 += weight * det * fv->snormal[a] * 
					(fv->v[a]-x_dot[a]) *
					bf[var]->phi[j] * delta(w,species_id);
				  }
				  J_AC[c] += d_term + d_term1;
				}
			      }
			    }

			  for(p=0; p<dim; p++)
			    {
			      var = VELOCITY1 + p;

			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = weight*det*fv->snormal[p]*fv->c[species_id]*bf[var]->phi[j];
				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    }
				}
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;
			      
			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 =0.0;

				      for( a=0; a<VIM; a++)
					{
					  d_term += -weight*(fv->dsurfdet_dx[p][j]*fv->snormal[a]*
				           (mp->diffusivity[sp]*fv->grad_c[sp][a] + kapta[species_id]*fv->grad_V[a]) +  
				           det*(d_diff[sp][var]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_c[sp][a] +
                                                d_kapta_dx[species_id][a]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_V[a] +
				                mp->diffusivity[sp]*fv->dsnormal_dx[a][p][j]*fv->grad_c[species_id][a] +
                                                kapta[species_id]*fv->dsnormal_dx[a][p][j]*fv->grad_V[a] +
				                mp->diffusivity[sp]*fv->snormal[a]*fv->d_grad_c_dmesh[a][sp][p][j] +
                                                kapta[species_id]*fv->snormal[a]*fv->d_grad_V_dmesh[a][p][j]));

					  d_term1 += weight*(fv->v[a] - x_dot[a])*fv->c[sp] *
					                    (fv->dsurfdet_dx[p][j] * fv->snormal[a] +
							     det * fv->dsnormal_dx[a][p][j] );
					}
				      J_AC[ ei->gun_list[var][j]  ] += d_term + d_term1;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if( pd->v[var])
			    {
			      for( j=0; j<ei->dof[var]; j++)
				{
				  d_term = 0;
				  
				  for( a=0; a<VIM; a++)
				    {
				      d_term += -weight*det*bf[var]->phi[j]*fv->snormal[a]*
                                                 (d_diff[sp][var]*fv->grad_c[sp][a]+
                                                  d_kapta_dT[species_id]*fv->grad_V[a]);
				    }
				  J_AC[ ei->gun_list[var][j]  ] += d_term;
				}
			    }				  

			  var = VOLTAGE;

			  if( pd->v[var])
			    {
			      for( j=0; j<ei->dof[var]; j++)
				{
				  d_term = 0;
				  
				  for( a=0; a<VIM; a++)
				    {
				      d_term += -weight*det*kapta[species_id]*bf[var]->grad_phi[j][a]*fv->snormal[a];
				    }
				  J_AC[ ei->gun_list[var][j]  ] += d_term;
				}
			    }				  

			  break;

			case CURRENT_FICKIAN:

			    var = MASS_FRACTION;
			    if (pd->v[var]) {
			      for (w = 0; w < pd->Num_Species_Eqn; w++) {
				for (j = 0; j < ei->dof[var]; j++) {
				  /*
				   * Find the material index for the current
				   * local variable degree of freedom
				   * (can't just query gun_list for MASS_FRACTION
				   *  unknowns -> have to do a lookup)
				   */		
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  d_term = d_term1 = 0.0;
				  for (a = 0; a < dim; a++) {
				    d_term += - weight * det * fv->snormal[a] *
					(d_diff[sp][MAX_VARIABLE_TYPES + w] * 
					 bf[var]->phi[j] * fv->grad_c[species_id][a] +
					 mp->diffusivity[species_id] *
					 bf[var]->grad_phi[j][a] * delta(w, species_id) 
                                         + d_kapta_dc[species_id][MAX_VARIABLE_TYPES + w]*
                                         bf[var]->phi[j]*fv->grad_V[a]);

				    d_term1 += weight * det * fv->snormal[a] * 
					(fv->v[a]-x_dot[a]) *
					bf[var]->phi[j] * delta(w,species_id);
				  }
				  J_AC[c] += d_term + d_term1;
				}
			      }
			    }

			  for(p=0; p<dim; p++)
			    {
			      var = VELOCITY1 + p;

			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = weight*det*fv->snormal[p]*fv->c[species_id]*bf[var]->phi[j];
				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    }
				}
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;
			      
			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = d_term2 =0.0;

				      for( a=0; a<VIM; a++)
					{
					  d_term += -weight*(fv->dsurfdet_dx[p][j]*fv->snormal[a]*
				           (mp->diffusivity[sp]*fv->grad_c[sp][a] + kapta[species_id]*fv->grad_V[a]) +  
				           det*(d_diff[sp][var]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_c[sp][a] +
                                                d_kapta_dx[species_id][a]*bf[var]->phi[j]*fv->snormal[a]*fv->grad_V[a] +
				                mp->diffusivity[sp]*fv->dsnormal_dx[a][p][j]*fv->grad_c[species_id][a] +
                                                kapta[species_id]*fv->dsnormal_dx[a][p][j]*fv->grad_V[a] +
				                mp->diffusivity[sp]*fv->snormal[a]*fv->d_grad_c_dmesh[a][sp][p][j] +
                                                kapta[species_id]*fv->snormal[a]*fv->d_grad_V_dmesh[a][p][j]));

					  d_term1 += weight*(fv->v[a] - x_dot[a])*fv->c[sp] *
					                    (fv->dsurfdet_dx[p][j] * fv->snormal[a] +
							     det * fv->dsnormal_dx[a][p][j] );
					}
				      J_AC[ ei->gun_list[var][j]  ] += d_term + d_term1;
				    }
				}
			    }

			  var = TEMPERATURE;

			  if( pd->v[var])
			    {
			      for( j=0; j<ei->dof[var]; j++)
				{
				  d_term = 0;
				  
				  for( a=0; a<VIM; a++)
				    {
				      d_term += -weight*det*bf[var]->phi[j]*fv->snormal[a]*
                                                 (d_diff[sp][var]*fv->grad_c[sp][a]+
                                                  d_kapta_dT[species_id]*fv->grad_V[a]);
				    }
				  J_AC[ ei->gun_list[var][j]  ] += d_term;
				}
			    }				  

			  var = VOLTAGE;

			  if( pd->v[var])
			    {
			      for( j=0; j<ei->dof[var]; j++)
				{
				  d_term = 0;
				  
				  for( a=0; a<VIM; a++)
				    {
				      d_term += -weight*det*kapta[species_id]*bf[var]->grad_phi[j][a]*fv->snormal[a];
				    }
				  J_AC[ ei->gun_list[var][j]  ] += d_term;
				}
			    }				  

			  break;

			case CURRENT:

			  var = VOLTAGE;
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {
			      d_term=0;
			      
			      for( a=0; a<dim; a++)
				{
				  d_term += weight*det*fv->snormal[a]*
				           ( -dkdV[j]*fv->grad_V[a] - k * bf[var]->grad_phi[j][a]);
				}

			      /* J_AC[ ei->gun_list[var][j]  ] += d_term; */
			      J_AC[ ei->gun_list[var][j]  ] = 1.0;
			    }

			  var = TEMPERATURE;
			  
			  for(j=0; j<ei->dof[var]; j++)
			    {
			      d_term=0;
			      
			      for( a=0; a<dim; a++)
				{
				  d_term += weight*det*fv->snormal[a]*
				           ( -dkdT[j]*fv->grad_V[a]);
				}

			      /* J_AC[ ei->gun_list[var][j]  ] += d_term; */
			      J_AC[ ei->gun_list[var][j]  ] = 1.0;
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = VELOCITY1 +p;
			      
			      if( pd->v[var])
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = 0.0;
				      /* J_AC[ ei->gun_list[var][j]  ] += d_term; */
				      J_AC[ ei->gun_list[var][j]  ] = 1.0;
				    } 
				}
			    }

			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if( pd->v[var] )
				{
				  
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = d_term1 = 0.0;

				      for( a=0; a<dim; a++)
					{
					  d_term += -weight * ( fv->dsurfdet_dx[p][j] * k * fv->snormal[a] * fv->grad_V[a] +
							       det * dkdX[p][j] *  fv->snormal[a] * fv->grad_V[a] +
							       det * k * fv->dsnormal_dx[a][p][j] * fv->grad_V[a] + 
							       det * k * fv->snormal[a] * fv->d_grad_V_dmesh[a][p][j] );
					}
				      
				      /* J_AC[ ei->gun_list[var][j]  ] += d_term; */
				      J_AC[ ei->gun_list[var][j]  ] = 1.0;
				    }
				}
			    }

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (w = 0; w < pd->Num_Species_Eqn; w++) {
			      for (j = 0; j < ei->dof[var]; j++) {
				/*
				 * Find the material index for the current
				 * local variable degree of freedom
				 * (can't just query gun_list for MASS_FRACTION
				 *  unknowns -> have to do a lookup)
				 */		
				gnn = ei->gnn_list[var][j];
				ledof = ei->lvdof_to_ledof[var][j];
				matIndex = ei->matID_ledof[ledof];
				c = Index_Solution(gnn, var, w, 0, matIndex);
				for (d_term = 0.0, a = 0; a < dim; a++) {
				  d_term += - weight * det * dkdC[w][j] * fv->snormal[a] * fv->grad_V[a];
				}
				/* J_AC[c] += d_term; */
				J_AC[c] = 1.0;
			      }
			    }
			  }
			  break;

			case AVERAGE_CONC:
			  
			  for( p=0; p<dim; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;
			      
			      if( pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      d_term = weight*fv->dsurfdet_dx[p][j]*fv->c[species_id];

				      J_AC[ ei->gun_list[var][j]  ] += d_term;
				    }
				}
			    }

			  var = MASS_FRACTION;
			  if (pd->v[var]) {
			    for (j = 0 ; j < ei->dof[var]; j++) {
			      /*
			       * Find the material index for the current
			       * local variable degree of freedom
			       * (can't just query gun_list for MASS_FRACTION
			       *  unknowns -> have to do a lookup)
			       */		
			      gnn = ei->gnn_list[var][j];
			      ledof = ei->lvdof_to_ledof[var][j];
			      matIndex = ei->matID_ledof[ledof];
			      c = Index_Solution(gnn, var, w, 0, matIndex);
			      d_term = weight * det * bf[var]->phi[j];
			      J_AC[c] += d_term;
			    }
			  }
			  break;

			case AREA:
			  
			  for(p=0; p<dim; p++ )
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if( pd->v[var])
				{
				  for (j=0; j<ei->dof[var]; j++)
				    {
				      J_AC[ ei->gun_list[var][j]  ] += weight*fv->dsurfdet_dx[p][j];
				    }
				}
			    }
			  break;

                        case VOL_REVOLUTION:

                          for(p=0; p<dim; p++ )
                            {
                              var = MESH_DISPLACEMENT1 + p;

                              if( pd->v[var])
                                {
                                  for (j=0; j<ei->dof[var]; j++)
                                    {
                                      J_AC[ ei->gun_list[var][j]  ] += 0.5*weight*
                     SGN(fv->snormal[1])*(fv->x[1]*(fv->snormal[1]*fv->dsurfdet_dx[p][j]
                                                   + det*fv->dsnormal_dx[1][p][j])
                         + det*fv->snormal[1]*bf[var]->phi[j]*delta(p,1));
                                    }
                                }
                            }
                          break;

			case N_DOT_X:
			  
			  for ( p=0; p<DIM; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if (pd->v[var] )
				{
				  for( j=0; j<ei->dof[var]; j++)
				    {
				      double term1=0.0, term2=0.0, term3 = 0.0;
				      
				      term1 = fv->snormal[p]*bf[var]->phi[j] / (double) dim;
				     
				      for( a=0 ; a < dim ; a++)
					{
					  term2 += ( fv->x[a] - param[a] ) * fv->dsnormal_dx[a][p][j]/(double) dim;
					  term3 += ( fv->x[a] - param[a] ) * fv->snormal[p]/(double) dim;
					}

				      J_AC[ ei->gun_list[var][j]  ] += weight * term1 * det +
					                               weight * term2 * det +
					                               weight * term3 * fv->dsurfdet_dx[p][j];
				    }
				}
			    }
				      

			  break;

 			case ELEC_FORCE_NORMAL:
 
 			  var = VOLTAGE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term1=0;
 			      for( a=0; a<dim; a++)
 				{
 				d_term1 += efield[a]*(-bf[var]->grad_phi[j][a]);
 				}
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 			        for( b=0; b<dim; b++)
 				  {
 				  d_term += weight*det*fv->snormal[a]*fv->snormal[b]*
 				      ( d_p_dV[j]* es[a][b]  + 
 					perm*((-bf[var]->grad_phi[j][a])*efield[b] +
 					efield[a]*(-bf[var]->grad_phi[j][b]) -
 					delta(a,b)*d_term1));
  				  }
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 			        for( b=0; b<dim; b++)
 				  {
 				  d_term += weight*det*fv->snormal[a]*fv->snormal[b]*
 				      ( d_p_dT[j]* es[a][b] );
 				  }
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  for( p=0; p<dim; p++)
 			    {
 			      var = MESH_DISPLACEMENT1 + p;
 
 			      if( pd->v[var] )
 				{
 				  
 				  for( j=0; j<ei->dof[var]; j++)
 				    {
 				      d_term = d_term1 = 0.0;
 			      		for( a=0; a<dim; a++)
 						{
 					d_term1 += efield[a]*(-fv->d_grad_V_dmesh[a][p][j]);
 						}
 				      for( a=0; a<dim; a++)
 					{
 			        	  for( b=0; b<dim; b++)
 				  	    {
 					  d_term += weight * ( fv->dsurfdet_dx[p][j] * 
 						fv->snormal[a]*perm*es[a][b]*fv->snormal[b] +
 						det * d_p_dX[p][j] *  fv->snormal[a] * 
 							es[a][b]*fv->snormal[b] +
 					det * perm * fv->dsnormal_dx[a][p][j] * es[a][b]*fv->snormal[b] + 
 					det * perm * fv->dsnormal_dx[b][p][j] * es[a][b]*fv->snormal[a] + 
 				       det * perm * fv->snormal[a] * fv->snormal[b] *
 		(efield[a]*(-fv->d_grad_V_dmesh[b][p][j]) + (-fv->d_grad_V_dmesh[a][p][j])*efield[b] - delta(a,b)*d_term1 ) );
 					    }
 					}
 				      
 				      J_AC[ ei->gun_list[var][j]  ] += d_term;
 				    }
 				}
 			    }
 
 			  break;
 
 			case ELEC_FORCE_TANGENT1:
 			case ELEC_FORCE_TANGENT2:
 
 			  dir = quantity == ELEC_FORCE_TANGENT1 ? 0 : 1;
 
 			  var = VOLTAGE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term1=0;
 			      for( a=0; a<dim; a++)
 				{
 				d_term1 += efield[a]*(-bf[var]->grad_phi[j][a]);
 				}
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 			        for( b=0; b<dim; b++)
 				  {
 				  d_term += weight*det*fv->stangent[dir][a]*-fv->snormal[b]*
 				      ( d_p_dV[j]* es[a][b]  + 
 					perm*((-bf[var]->grad_phi[j][a])*efield[b] +
 					efield[a]*(-bf[var]->grad_phi[j][b]) -
 					delta(a,b)*d_term1));
  				  }
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 			        for( b=0; b<dim; b++)
 				  {
  				  d_term += weight*det*fv->stangent[dir][a]*-fv->snormal[b]*
 				      ( d_p_dT[j]* es[a][b] );
 				  }
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  for( p=0; p<dim; p++)
 			    {
 			      var = MESH_DISPLACEMENT1 + p;
 
 			      if( pd->v[var] )
 				{
 				  
 				  for( j=0; j<ei->dof[var]; j++)
 				    {
 				      d_term = d_term1 = 0.0;
 			      		for( a=0; a<dim; a++)
 						{
 					d_term1 += efield[a]*(-fv->d_grad_V_dmesh[a][p][j]);
 						}
 				      for( a=0; a<dim; a++)
 					{
 			        	  for( b=0; b<dim; b++)
 				  	    {
 					  d_term += weight * ( fv->dsurfdet_dx[p][j] * 
  						fv->stangent[dir][a]*perm*es[a][b]*-fv->snormal[b] +
 						det * d_p_dX[p][j] *  fv->stangent[dir][a] * 
  							es[a][b]*-fv->snormal[b] +
 					det * perm * fv->dstangent_dx[dir][a][p][j] * es[a][b]*-fv->snormal[b] + 
 					det * perm * -fv->dsnormal_dx[b][p][j] * es[a][b]*fv->stangent[dir][a] + 
 				       det * perm * fv->stangent[dir][a] * -fv->snormal[b] *
 		(efield[a]*(-fv->d_grad_V_dmesh[b][p][j]) + (-fv->d_grad_V_dmesh[a][p][j])*efield[b] - delta(a,b)*d_term1 ) );
 					    }
 					}
 				      
 				      J_AC[ ei->gun_list[var][j]  ] += d_term;
 				    }
 				}
 			    }
 
 			  break;
 
 			case ELEC_FORCE_X:
 			case ELEC_FORCE_Y:
 			case ELEC_FORCE_Z:
 
 			  dir = ( quantity == ELEC_FORCE_X ? 0 :
 				( quantity == ELEC_FORCE_Y ? 1 : 
 			        ( quantity == ELEC_FORCE_Z ? 2 : -1) ) ) ;  /* Can't do this with FORTRAN */
 
 			  EH(dir, "PANIC in evaluate_flux sensitivity section");
 
 			  var = VOLTAGE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term1=0;
 			      for( a=0; a<dim; a++)
 				{
 				d_term1 += efield[a]*(-bf[var]->grad_phi[j][a]);
 				}
 			      d_term=0;
 			      
 			        for( b=0; b<dim; b++)
 				  {
 				  d_term += weight*det*-fv->snormal[b]*
 				      ( d_p_dV[j]* es[dir][b]  + 
 					perm*((-bf[var]->grad_phi[j][dir])*efield[b] +
 					efield[dir]*(-bf[var]->grad_phi[j][b]) -
 					delta(dir,b)*d_term1));
  				  }
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term=0;
 			      
 			        for( b=0; b<dim; b++)
 				  {
 				  d_term += weight*det*-fv->snormal[b]*
 				      ( d_p_dT[j]* es[dir][b] );
 				  }
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  for( p=0; p<dim; p++)
 			    {
 			      var = MESH_DISPLACEMENT1 + p;
 
 			      if( pd->v[var] )
 				{
 				  
 				  for( j=0; j<ei->dof[var]; j++)
 				    {
 				      d_term = d_term1 = 0.0;
 			      		for( a=0; a<dim; a++)
 						{
 					d_term1 += efield[a]*(-fv->d_grad_V_dmesh[a][p][j]);
 						}
 			        	  for( b=0; b<dim; b++)
 				  	    {
 					  d_term += weight * ( fv->dsurfdet_dx[p][j] * 
 						perm*es[dir][b]*-fv->snormal[b] +
 						det * d_p_dX[p][j] * es[dir][b]*-fv->snormal[b] +
 					det * perm * -fv->dsnormal_dx[b][p][j] * es[dir][b] + 
 				       det * perm * -fv->snormal[b] *
 		(efield[dir]*(-fv->d_grad_V_dmesh[b][p][j]) + (-fv->d_grad_V_dmesh[dir][p][j])*efield[b] - delta(dir,b)*d_term1 ) );
 					    }
 				      
 				      J_AC[ ei->gun_list[var][j]  ] += d_term;
 				    }
 				}
 			    }
 
 			  break;
 
 			case NET_SURF_CHARGE:
 
 			  var = VOLTAGE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 				  d_term += -weight*det*fv->snormal[a]*
 				      ( d_p_dV[j]*efield[a] + perm * (-bf[var]->grad_phi[j][a]));
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
 			      d_term=0;
 			      
 			      for( a=0; a<dim; a++)
 				{
 				  d_term += -weight*det*fv->snormal[a]*
 				           ( d_p_dT[j]*efield[a]);
 				}
 
 			      J_AC[ ei->gun_list[var][j]  ] += d_term; 
 			    }
 
 			  for( p=0; p<dim; p++)
 			    {
 			      var = MESH_DISPLACEMENT1 + p;
 
 			      if( pd->v[var] )
 				{
 				  
 				  for( j=0; j<ei->dof[var]; j++)
 				    {
 				      d_term = d_term1 = 0.0;
 
 				      for( a=0; a<dim; a++)
 					{
 					  d_term += -weight * ( fv->dsurfdet_dx[p][j] * fv->snormal[a]*perm*efield[a] +
 							       det * d_p_dX[p][j] *  fv->snormal[a] * efield[a] +
 							       det * perm * fv->dsnormal_dx[a][p][j] * efield[a] + 
 							       det * perm * fv->snormal[a] * (-fv->d_grad_V_dmesh[a][p][j]) );
 					}
 				      
 				      J_AC[ ei->gun_list[var][j]  ] += d_term;
 				    }
 				}
 			    }
 
 			  break;

			default:
			  EH(-1,"Constraint sensitivities haven't been implemented for that flux type.");
			  break;
			}
		    }

		}   /*  surface integration over element */

	      if( Subgrid_Int.active )
		{
		  safe_free ( (void *) Subgrid_Int.s  ); Subgrid_Int.s = NULL;
		  safe_free ( (void *) Subgrid_Int.wt ); Subgrid_Int.wt = NULL;
		  free_search_grid ( &element_search_grid );
		}


	    }   /*  material id conditional */
	}    /*   element loop */
    }      /*   sset id   */
  else
    {

/**  Apply end point conditions when the nset is not found   **/
  	nset_id = in_list(side_set_id, 0, exo->num_node_sets, exo->ns_id);
  	if ( nset_id != -1 )
    		{
		int corner_elem=-1,local_node, gnn, dir, kine_sset;
		double sign;
      		num_node_in_set      = exo->ns_num_nodes[nset_id];
      		num_dist_fact_in_set = exo->ns_num_distfacts[nset_id];
      		elem_list = &exo->ns_node_list[exo->ns_node_index[nset_id]];
		if (num_node_in_set != 1) 
		   EH(-1,"more than one node, this is for 2D only");
		gnn = elem_list[0];
		sign = (species_id < 0 ) ? -1. : 1.;
		kine_sset = abs(species_id);
	/*  Find the free surface element in the adjoining sset	*/
  		sset_id = in_list(kine_sset, 0, exo->num_side_sets, exo->ss_id);
		if( sset_id != -1)
		  {
      		  num_nodes_on_side    = exo->ss_num_distfacts[sset_id]/exo->ss_num_sides[sset_id]; /* Well... */
		  for(i=0 ; i<exo->ss_num_sides[sset_id] ; i++)
		     {
		      ielem = exo->ss_elem_list[exo->ss_elem_index[sset_id]+i];
		      if( in_list(gnn, exo->elem_node_pntr[ielem],
				exo->elem_node_pntr[ielem+1],
				exo->elem_node_list) != -1)
				{
				corner_elem = ielem;
     		local_node = in_list(gnn, exo->elem_node_pntr[corner_elem],
                                exo->elem_node_pntr[corner_elem+1],
                                exo->elem_node_list);
		local_node -= exo->elem_node_pntr[corner_elem];

	  	   mn = find_mat_number(corner_elem, exo);

	        if ( blk_id == -1 || 
                     (mn == map_mat_index(blk_id) && 
                      dpi->elem_owner[corner_elem] == ProcID)) {

	  	   err = load_elem_dofptr( corner_elem, 
				  (Exo_DB*) exo,
				  (dbl *) x,
				  (dbl *) x,
				  (dbl *) xdot,
				  (dbl *) xdot,
				  (dbl *) x,
				  0);
	  	   EH(err, "load_elem_dofptr");

	  	   err = bf_mp_init(pd);
	  	   EH(err, "bf_mp_init");
  
	           iconnect_ptr        = ei->iconnect_ptr;
	           ielem_type          = ei->ielem_type;
	      	   ip_total            = elem_info(NQUAD_SURF, ielem_type);
	      	   num_local_nodes     = ei->num_local_nodes;
	           ielem_dim           = ei->ielem_dim;
	      
	 
	           id_side = find_id_side (ei->ielem, num_nodes_on_side,
				      &exo->ss_node_list[sset_id]
				      [num_nodes_on_side*i],
				      id_local_elem_coord, exo);

            	   find_nodal_stu(local_node, ei->ielem_type, xi, xi+1, xi+2);

            	   err = load_basis_functions( xi, bfd);
            	   EH( err, "problem from load_basis_functions");

            	   err = beer_belly();
            	   EH( err, "beer_belly");

            	   err = load_fv();
            	   EH( err, "load_fv");

            	   err = load_bf_grad();
            	   EH( err, "load_bf_grad");

            	   err = load_bf_mesh_derivs();
            	   EH(err, "load_bf_mesh_derivs");

            	   surface_determinant_and_normal(corner_elem,
                                           exo->elem_node_pntr[corner_elem],
                                           ei->num_local_nodes,
                                           ei->ielem_dim-1,
                                           id_side,
                                           num_nodes_on_side,
                                           id_local_elem_coord);
	 	   err = load_fv_grads();
	 	   EH( err, "load_fv_grads");
			  
	 	   err = load_fv_mesh_derivs(1);
	 	   EH( err, "load_fv_mesh_derivs");

	    	   if (ielem_dim !=3)
	      		{
			 calc_surf_tangent (corner_elem, iconnect_ptr, 
				   num_local_nodes, ielem_dim-1,
				   num_nodes_on_side,
				   id_local_elem_coord);
	      		}
		dim = pd_glob[0]->Num_Dim;

 		/*
 		 * load surface tension for variable models
 		 */
 		      if (mp->SurfaceTensionModel == CONSTANT) {
                            dsigmadT = 0.0;
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = 0.0;}
                            }
 		      if (mp->SurfaceTensionModel != CONSTANT) {
 		            load_surface_tension(dsigma_dx);
 		            }
 		      if (mp->SurfaceTensionModel == USER) {
                            dsigmadT = mp->d_surface_tension[TEMPERATURE];
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = mp->d_surface_tension[MAX_VARIABLE_TYPES+a];}
                            }
		  switch( quantity )
		    {
		    case FORCE_X:
		    case FORCE_Y:
		    case FORCE_Z:
		    case FORCE_X_POS:
		    case FORCE_Y_POS:
		    case FORCE_Z_POS:
		    case FORCE_X_NEG:
		    case FORCE_Y_NEG:
		    case FORCE_Z_NEG:
			  dir = ( quantity == FORCE_X ? 0 :
				( quantity == FORCE_Y ? 1 : 
			        ( quantity == FORCE_Z ? 2 : 
			        ( quantity == FORCE_X_POS ? 0 :
				( quantity == FORCE_Y_POS ? 1 : 
			        ( quantity == FORCE_Z_POS ? 2 : 
			        ( quantity == FORCE_X_NEG ? 0 :
				( quantity == FORCE_Y_NEG ? 1 : 
			        ( quantity == FORCE_Z_NEG ? 2 : -1))))))))) ;  

		      if(cr->MeshMotion == ARBITRARY)
			{
		       if( pd->CoordinateSystem == CARTESIAN ) 
			  {
                         local_flux += sign*mp->surface_tension*fv->stangent[0][dir];
 	                if( J_AC != NULL)
 		           {
			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      J_AC[ ei->gun_list[var][j] ] += 
                         sign*mp->surface_tension*fv->dstangent_dx[0][dir][p][j];
				    }
				}
			    }
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
				      J_AC[ ei->gun_list[var][j] ] += 
                         	         sign*fv->stangent[0][dir]
				         *dsigmadT*bf[var]->phi[j];
 			    }
			    var = MASS_FRACTION;
			    if (pd->v[var]) {
			      for (w = 0; w < pd->Num_Species_Eqn; w++) {
				for (j = 0; j < ei->dof[var]; j++) {
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  J_AC[c] += sign*fv->stangent[0][dir]
                                    *dsigmadC[w]*bf[var]->phi[j];
				}
			      }
			    }
 		           }
			}
		      else if(pd->CoordinateSystem == SWIRLING || 
			 pd->CoordinateSystem == CYLINDRICAL)
			{
                         local_flux += 2*M_PIE*fv->x[1]*sign*mp->surface_tension*fv->stangent[0][dir];
 	                if( J_AC != NULL)
 		           {
			  for( p=0; p<dim ; p++)
			    {
			      var = MESH_DISPLACEMENT1 + p;

			      if(pd->v[var])
				{
				  for( j=0 ; j<ei->dof[var]; j++)
				    {
				      J_AC[ ei->gun_list[var][j] ] += 
                         2*M_PIE*sign*mp->surface_tension*
			(fv->x[1]*fv->dstangent_dx[0][dir][p][j]+
			delta(p,1)*bf[var]->phi[j]*fv->stangent[0][dir]);
				    }
				}
			    }
 			  var = TEMPERATURE;
 			  
 			  for(j=0; j<ei->dof[var]; j++)
 			    {
				      J_AC[ ei->gun_list[var][j] ] += 
                         2*M_PIE*sign*fv->x[1]*fv->stangent[0][dir]
				*dsigmadT*bf[var]->phi[j];
 			    }
			    var = MASS_FRACTION;
			    if (pd->v[var]) {
			      for (w = 0; w < pd->Num_Species_Eqn; w++) {
				for (j = 0; j < ei->dof[var]; j++) {
				  gnn = ei->gnn_list[var][j];
				  ledof = ei->lvdof_to_ledof[var][j];
				  matIndex = ei->matID_ledof[ledof];
				  c = Index_Solution(gnn, var, w, 0, matIndex);
				  J_AC[c] += 2*M_PIE*fv->x[1]*sign*fv->stangent[0][dir]
                                    *dsigmadC[w]*bf[var]->phi[j];
				}
			      }
			    }
 		           }
			}
			else
			    {
			EH(-1, "Point force has not been updated for the PROJECTED_CARTESIAN coordinate system.");
			    }
			}
		      break;
		    default:
		      EH(-1, "Illegal flux type");
		      break;
		    }  /*  end of switch */
#ifdef PARALLEL
        delta_flux = local_flux - local_flux0;
        delta_flux_conv = local_flux_conv - local_flux_conv0;
        delta_area = local_area - local_area0;

        if( Num_Proc > 1 &&  dpi->elem_owner[ elem_list[i] ] == ProcID )
         {
          proc_flux += delta_flux;
          proc_flux_conv += delta_flux_conv;
          proc_area += delta_area;
         }

        local_flux0 = local_flux;
        local_flux_conv0 = local_flux_conv;
        local_area0 = local_area;
#endif
        if (profile_flag && print_flag && 
		(quantity != LS_DCA || (ierr && ip == 0)) ) {
        FILE  *jfp;
        if( (jfp=fopen(filenm,"a")) != NULL)
            {
             fprintf(jfp," %g  %g  %g  %g  %g",
                        fv->x[0],fv->x[1],fv->x[2],local_q,local_qconv);
#if PRINT_VISCOSITY
             fprintf(jfp," %g  %g ",
                         gamma_dot, mu);
#endif
#if PRINT_SURFACE_TENSION
             fprintf(jfp," %g ",
                         fv->T);
/*                         mp->surface_tension); */
#endif
              fprintf(jfp," \n");
              fflush(jfp);
            }
          fclose(jfp);
                }
       /* Compute sensitivities if requested */
	if ( J_AC != NULL)
		    {
		      double d_term;

		      switch (quantity)
			{
	 /* FORCE_X, FORCE_Y, FORCE_Z */
	 /* FORCE_X_POS, FORCE_Y_POS, FORCE_Z_POS */
	 /* FORCE_X_NEG, FORCE_Y_NEG, FORCE_Z_NEG */
			case FORCE_X:
			case FORCE_Y:
			case FORCE_Z:
			case FORCE_X_POS:
			case FORCE_Y_POS:
			case FORCE_Z_POS:
			case FORCE_X_NEG:
			case FORCE_Y_NEG:
			case FORCE_Z_NEG:
			   
			  dir = ( quantity == FORCE_X ? 0 :
				( quantity == FORCE_Y ? 1 : 
				( quantity == FORCE_Z ? 2 : 
			        ( quantity == FORCE_X_POS ? 0 :
				( quantity == FORCE_Y_POS ? 1 : 
				( quantity == FORCE_Z_POS ? 2 : 
			        ( quantity == FORCE_X_NEG ? 0 :
				( quantity == FORCE_Y_NEG ? 1 : 
				( quantity == FORCE_Z_NEG ? 2 : -1)))))))));  

			  EH(dir, "PANIC in evaluate_flux sensitivity section");

		      if(cr->MeshMotion == ARBITRARY)
			{
			  for( b=0 ; b < dim ; b++)
			    {
			      var = MESH_DISPLACEMENT1 + b ;

			      if( pd->v[var] )
				{
				  for(j=0 ; j < ei->dof[var]; j++)
				    {
				     d_term = sign*mp->surface_tension * 
					      fv->dstangent_dx[0][dir][b][j];

				      J_AC[ ei->gun_list[var][j] ] += d_term;

				    }
				}
			     }
			}
			  break;
			default:
			  EH(-1,"Constraint sensitivities haven't been implemented for that flux type.");
			  break;
			}
		      }  /*  J_AC  */
					}	/*  blk_id  */
				}
		       }  /* ss_sides loop	*/
		   if(corner_elem == -1) EH(-1,"corner element not found");

		  }	/* if sset_id	*/
    		}	/* if nset_id	*/
  		else
    		{
#ifndef PARALLEL
      (void) sprintf(Err_Msg, "%s could not locate SSID %d.", yo, side_set_id);
      EH(-1, Err_Msg);
#endif      
    		}

    }
         

#ifdef PARALLEL
  if( Num_Proc > 1 )
   {
    MPI_Allreduce(&proc_flux, &global_flux, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD);
    MPI_Allreduce(&proc_flux_conv, &global_flux_conv, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD);

    MPI_Allreduce(&proc_area, &global_area, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD);
    local_flux = global_flux;
    local_flux_conv = global_flux_conv;
    local_area = global_area;
   }
#endif  

/**   skip this section if not printing - i.e. for AC conditions  **/

  if(print_flag && ProcID == 0)
    {
      FILE  *jfp;
      if( (jfp=fopen(filenm,"a")) != NULL)
	{ 
	  if(quantity==TORQUE)
	    {
              if(pd->CoordinateSystem == SWIRLING || 
                 pd->CoordinateSystem == CYLINDRICAL)
		{
                   fprintf(jfp," torque= %e   \n", local_flux);
		}
	      else   /* CARTESIAN */
		{
                   fprintf(jfp," torque= %e %e %e \n", Torque[0],Torque[1],Torque[2]);
		}
	      fprintf(jfp, "\n"); 
	      fflush(jfp);	
	    }
	  else
	    { 
              fprintf(jfp," flux=  %e %e  area= %e  ", local_flux,local_flux_conv,local_area);
	        if(quantity==REPULSIVE_FORCE)
	           {
		    for( b=0 ; b<DIM ; b++)
                       {Torque[b] /= local_flux;}
      fprintf(jfp," centers=  %e %e %e  ", Torque[0],Torque[1],Torque[2]);
                   }
	      fprintf(jfp, "\n\n"); 
	      fflush(jfp);	
	    }
	  fclose(jfp);
	}
    }
  af->Assemble_Jacobian = Jac_state;
  return( local_flux + local_flux_conv );		/* failsafe default? */
}


/* evalutate_volume_integral - integrate a volumetric quantity and print out
 *
 * Author: T.A. Baer
 * Date  : April , 2000
 *
 *
 *  Here is how to add another quantity
 *
 *  Step 1 :  Add an index name and an index value in mm_post_def.h.  For example I_VOLUME
 *  Step 2 :  Pick a representative string name for your quantity and add it with the new index
00 *            value in the structure array pp_vol_names structure in mm_post_def.h.  Don't forget
 *            to increase the dimension of this structure by 1.
 *  Step 3 :  Install a new switch branch in compute_volume_integrand in mm_flux.c and add code
 *            to evalute the integrand of your quantity.  Note that weight is the integration weight
 *            and det is the mapping jabobian determinant.  Reference your case with the index 
 *            name you defined earlier ( e.g. I_VOLUME ).  The code you define here will be 
 *            evaluate at each integration point in every element.
 *
 *            And that's it.  What could be easier ?  Use and enjoy your new volume integral.
 */

double
evaluate_volume_integral(const Exo_DB *exo, /* ptr to basic exodus ii mesh information */
			 const Dpi *dpi, /* distributed processing info */
			 const int quantity, /* to pick VOLUME, DISSIPATION, etc. */
			 const char *quantity_str, /* volume integral name */
			 const int blk_id,	/* material identification */
			 const int species_id, /* species identification */
			 const char *filenm, /* File name pointer */
			 const double *params,
                         const int num_params,
			 double *J_AC,   /* Pointer to AC sensitivity vector, may be NULL */
			 const double x[],	/* solution vector */
			 const double xdot[],	/* dx/dt vector */
			 const double delta_t, /* time-step size */
			 const double time_value, /* current time */
			 const int print_flag)     /*  flag for printing results,1=print*/
{
  int i,j;
  int eb,e_start,e_end, elem;

  int err=0;
  int mn, ip, ip_total;
  double sum = 0.0;
  extern int PRS_mat_ielem;
  extern int MMH_ip; 
#ifdef PARALLEL
  double sum0=0.0;
  double delta_sum=0.0;
  double proc_sum=0.0;
  double global_sum=0.0;
#endif
  NTREE *start_tree=NULL;
  SGRID *element_search_grid=NULL;


  int subgrid_integration_active = FALSE;
  int subelement_vol_integration_active = FALSE;
  int subelement_surf_integration_active = FALSE;
  int adaptive_integration_active = FALSE;

/*  variables for adaptive quadrature weight calculation   */
  double ls_F[MDE],ad_wt[MDE];
  int ierr = 0, wt_type;
  double xf3D[12][2],ecrd[12][MAX_PDIM];
  int side_id[12];
  double xf2D[6][2];
  int nint2D[6];


  adaptive_integration_active = ( ls != NULL && ls->AdaptIntegration &&
				 ( quantity == I_POS_FILL ||
				   quantity == I_NEG_FILL ||
				   quantity == I_POS_VOLPLANE ||
				   quantity == I_NEG_VOLPLANE ||
				   quantity == I_POS_CENTER_X ||
                                   quantity == I_POS_CENTER_Y ||
                                   quantity == I_POS_CENTER_Z ||
				   quantity == I_NEG_CENTER_X ||
                                   quantity == I_NEG_CENTER_Y ||
                                   quantity == I_NEG_CENTER_Z ||
				   quantity == I_POS_VX ||
				   quantity == I_POS_VY ||
				   quantity == I_POS_VZ ||
				   quantity == I_NEG_VX ||
				   quantity == I_NEG_VY ||
				   quantity == I_NEG_VZ ) );

  subgrid_integration_active = ( ls != NULL &&
				 ls->Integration_Depth > 0 &&
				 !adaptive_integration_active &&
				 ( quantity == I_POS_FILL ||
				   quantity == I_NEG_FILL ||
				   quantity == I_POS_VOLPLANE ||
				   quantity == I_NEG_VOLPLANE ||
                                   quantity == I_POS_CENTER_X ||
                                   quantity == I_POS_CENTER_Y ||
                                   quantity == I_POS_CENTER_Z ||
				   quantity == I_NEG_CENTER_X ||
                                   quantity == I_NEG_CENTER_Y ||
                                   quantity == I_NEG_CENTER_Z ||
				   quantity == I_LS_ARC_LENGTH ||
				   quantity == I_POS_VX ||
				   quantity == I_POS_VY ||
				   quantity == I_POS_VZ ||
				   quantity == I_NEG_VX || 
				   quantity == I_NEG_VY ||
				   quantity == I_NEG_VZ ) );
  subelement_vol_integration_active = ( ls != NULL &&
				        ls->SubElemIntegration &&
				 	!adaptive_integration_active &&
                                        ( params == NULL || params[0] == 0. ) &&
				        ( quantity == I_POS_FILL ||
				          quantity == I_NEG_FILL ||
				          quantity == I_POS_VOLPLANE ||
				          quantity == I_NEG_VOLPLANE ||
                                          quantity == I_POS_CENTER_X ||
                                          quantity == I_POS_CENTER_Y ||
                                          quantity == I_POS_CENTER_Z ||
				          quantity == I_NEG_CENTER_X ||
                                          quantity == I_NEG_CENTER_Y ||
                                          quantity == I_NEG_CENTER_Z ||
				          quantity == I_POS_VX ||
				          quantity == I_POS_VY ||
				          quantity == I_POS_VZ ||
				          quantity == I_NEG_VX ||
				          quantity == I_NEG_VY ||
				          quantity == I_NEG_VZ ) );
  subelement_surf_integration_active = ( ls != NULL &&
				         ls->SubElemIntegration &&
				 	 !adaptive_integration_active &&
                                         ( params == NULL || params[0] == 0. ) &&
				         ( quantity == I_LS_ARC_LENGTH ||
				           quantity == I_MAG_GRAD_FILL_ERROR ) );
                                   
  if( subgrid_integration_active )  start_tree = create_shape_fcn_tree( ls->Integration_Depth );


  /* first write time stamp or run stamp to separate the sets */

  if (print_flag && ProcID == 0) {
    FILE  *jfp;
    if( (jfp=fopen(filenm,"a")) != NULL) {
      if ( ppvi_type == PPVI_VERBOSE ) {
	fprintf(jfp,"Time/iteration = %e \n", time_value);
	fprintf(jfp,"\t  (%s) Volume Integral for block %d species %d\n", 
		quantity_str,blk_id, species_id);
      }
      if ( ppvi_type == PPVI_CSV ) {
	fprintf(jfp,"%e,", time_value);
      }
      fflush(jfp);
      fclose(jfp);
    }
  }

  mn = map_mat_index(blk_id);
  if( ( eb = in_list(blk_id, 0, exo->num_elem_blocks, exo->eb_id) ) != -1 )
    {

      e_start = exo->eb_ptr[eb];
      e_end   = exo->eb_ptr[eb+1];


      for( elem= e_start; elem < e_end; elem++)
        {
	  int Use_Subgrid_Integration = FALSE;
          int Use_Subelement_Integration = FALSE;

          double wt, xi[3];
	  double (*s)[DIM] = NULL, *weight = NULL;

          ei->ielem = elem;

	  /*needed for saturation hyst. func. */
	  PRS_mat_ielem = ei->ielem - exo->eb_ptr[find_elemblock_index(ei->ielem, exo)]; 

          err = load_elem_dofptr(elem, (Exo_DB*) exo,
			         (dbl *) x, (dbl *) x, (dbl *) xdot,
				 (dbl *) xdot, (dbl *) x, 0);
          EH(err, "load_elem_dofptr");

	  err = bf_mp_init(pd);
	  EH(err, "bf_mp_init");

          ip_total = elem_info(NQUAD, ei->ielem_type);

	  if ( subgrid_integration_active  )
	    {
	      double width = params == NULL ? ls->Length_Scale : 2.0*params[0];

	      if ((Use_Subgrid_Integration = current_elem_overlaps_interface(width)))
		{
		  ip_total = get_subgrid_integration_pts (  start_tree, &element_search_grid, &s, &weight, width );
		}
	    }
          else if ( subelement_vol_integration_active  )
	    {
              if ((Use_Subelement_Integration = current_elem_on_isosurface(FILL, 0.)))
		{
                  if ( quantity == I_NEG_FILL ||
		       quantity == I_NEG_VOLPLANE ||
		       quantity == I_NEG_VX ||
		       quantity == I_NEG_VY ||
		       quantity == I_NEG_VZ  ||
                       quantity == I_NEG_CENTER_X ||
                       quantity == I_NEG_CENTER_Y ||
                       quantity == I_NEG_CENTER_Z ) {
                    ip_total = get_subelement_integration_pts ( &s, &weight, NULL, 0., -2, -1 );
                  } else if ( quantity == I_POS_FILL ||
		              quantity == I_POS_VOLPLANE ||
		              quantity == I_POS_VX ||
		              quantity == I_POS_VY ||
		              quantity == I_POS_VZ  ||
                              quantity == I_POS_CENTER_X ||
                              quantity == I_POS_CENTER_Y ||
                              quantity == I_POS_CENTER_Z  ) {
                    ip_total = get_subelement_integration_pts ( &s, &weight, NULL, 0., -2, 1 );
                  } else {
                    ip_total = get_subelement_integration_pts ( &s, &weight, NULL, 0., -2, 0 );
                  }
                  ls->Elem_Sign = 0;
		}
	    }
          else if ( subelement_surf_integration_active  )
	    {
              ip_total = 0;
              if ((Use_Subelement_Integration = current_elem_on_isosurface(FILL, 0.)))
		{
                  ip_total = get_subelement_integration_pts ( &s, &weight, NULL, 0., -1, 0 );
                  ls->Elem_Sign = 0;
                  ls->on_sharp_surf = TRUE;
		}
	    }
 	  else if ( adaptive_integration_active )
 	    {
 		for(i=0;i<ei->num_local_nodes;i++)	{ls_F[i]=*esp_old->F[i];}
 		wt_type = 2;
 		if( quantity == I_NEG_FILL || 
				quantity == I_NEG_VOLPLANE ||
				quantity == I_NEG_CENTER_X ||
                                quantity == I_NEG_CENTER_Y ||
                                quantity == I_NEG_CENTER_Z ||
				quantity == I_NEG_VX ||
				quantity == I_NEG_VY ||
				quantity == I_NEG_VZ        )
 			{ 
			for(i=0;i<ei->num_local_nodes;i++)	{ls_F[i]=-ls_F[i];}
			 }
 		ierr = adaptive_weight( ad_wt, ip_total, ei->ielem_dim, ls_F, 0.0,
 				 wt_type, ei->ielem_type);
 		if(ierr == -1)printf("adaptive wt problem %d %g %g %g\n",ierr,ls_F[0],ls_F[1],ls_F[2]);
 	    }
/*  This section is a kludge for getting surface data for constant value surfaces out
*/
/*  element crossing interrogation */
	if( quantity == I_SURF_SPECIES ||
		quantity == I_SURF_TEMP )
	  {
 	   if(Num_Proc > 1)
 		EH(-1,"SURFACE_SPECIES not recommended in parallel\n");
	   if( quantity == I_SURF_SPECIES)
	     {
	   	for(i=0;i<ei->num_local_nodes;i++)	
		   {
 		    ls_F[i] = params[pd->Num_Species_Eqn] - params[pd->Num_Species_Eqn+1];
		    for(j=0;j<pd->Num_Species_Eqn;j++)
			{
			ls_F[i] += params[j]*(*esp_old->c[j][i]);
			}
 		    }
	     }
	   else if( quantity == I_SURF_TEMP)
	     {
		ierr=ei->dof[TEMPERATURE]/ei->num_local_nodes;
	   	for(i=0;i<ei->num_local_nodes;i++)	
		   {
 		    ls_F[i] = (*esp_old->T[ierr*i]) - params[0];
 		   }
	     }
	   else
	     EH(-1,"That SURF quantity not available.");

	  if(pd->Num_Dim == 3)
		{
		if( ei->ielem_type == TRILINEAR_HEX)
		   {
 		    ierr = interface_crossing_3DL( ls_F, xf3D, side_id, ecrd);
		   }
		else
	     	    EH(-1,"Only SURF 3D element is TRILINEAR_HEX.");
		}
	  else
		{
		if( ei->ielem_type == BIQUAD_QUAD)
		   {
 		    ierr = interface_crossing_2DQ( ls_F, xf2D, side_id, nint2D, ecrd);
		   }
		else
	     	    EH(-1,"Only SURF 2D element is BIQUAD_QUAD.");
		}
	if(ierr)
  	  {
 	  FILE *jfp;
	  sum += ((double)ierr);
 	  if( (jfp = fopen( filenm, "a")) != NULL )
 	     {
 	  for(i=0 ; i<ierr ; i++)
 		{
 		xi[0] = ecrd[i][0];
 		xi[1] = ecrd[i][1];
 		xi[2] = ecrd[i][2];
       
 	     	err = load_basis_functions( xi, bfd );
 	     	EH( err, "problem from load_basis_functions");
 
 	     	err = beer_belly();
 	     	EH( err, "beer_belly");
 	      
 	     	err = load_fv();
 	     	EH( err, "load_fv");
 
 		fprintf(jfp,"%d %d %g %g %g\n",elem,i,fv->x[0],fv->x[1],fv->x[2]);
 		}
   	   fclose(jfp);
 	      }
  	   }
	  }
	  
          for( ip=0; ip<ip_total; ip++)
            {

	      MMH_ip = ip;

	      if ( Use_Subgrid_Integration || Use_Subelement_Integration )
		{
/* 		  DPRINTF(stderr, "Subgrid integration on element : \t%d\n", elem ); */

		  xi[0] = s[ip][0];
		  xi[1] = s[ip][1];
		  xi[2] = s[ip][2];
		  fv->wt = wt = weight[ip];
		}
	      else
		{
		  /* find quadrature point locations for current ip */
		  find_stu (ip, ei->ielem_type, &xi[0], &xi[1], &xi[2]);

		  /* find quadrature weights for current ip */
 		 if(adaptive_integration_active)
 			{fv->wt = wt = ad_wt[ip_total-1-ip]; }
 		 else
 			{ fv->wt = wt = Gq_weight (ip, ei->ielem_type); }
		}
		  

	     /*
	      * Load up basis function information for ea variable...
	      */
      
	     err = load_basis_functions( xi, bfd );
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
	      
	     if (mp->PorousMediaType != CONTINUOUS){
	       err = load_porous_properties();
	       EH( err, "load_porous_properties");
	     }
	     do_LSA_mods(LSA_VOLUME);

             if ( subelement_surf_integration_active ) bf[pd->ShapeVar]->detJ = 1.;
             
             compute_volume_integrand( quantity, elem, species_id, 
                              params, num_params, &sum, J_AC, adaptive_integration_active,
				       time_value, delta_t, xi, exo);
#ifdef PARALLEL
             delta_sum = sum - sum0;

             if( Num_Proc > 1 && dpi->elem_owner[ ei->ielem ] == ProcID ) 
              {
               proc_sum += delta_sum;
              }

             sum0 = sum;
#endif
            }
	  if ( Use_Subgrid_Integration || Use_Subelement_Integration ) 
	    {
	      safe_free( s ) ;
	      safe_free( weight );
	      free_search_grid( &element_search_grid );
	    }

        }
    }

#ifdef PARALLEL
  if( Num_Proc > 1 ) {
    MPI_Allreduce( &proc_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 
		   MPI_COMM_WORLD);

    sum = global_sum;

  }
#endif

  if ( subgrid_integration_active ) 
  {
	  free_shape_fcn_tree( start_tree );
  }
  
  if ( ls != NULL ) ls->on_sharp_surf = FALSE;

  if (quantity == I_SPECIES_SOURCE)
      { 
       if(time_value <= tran->init_time+delta_t)
 	    { 
             Spec_source_inventory[mn][species_id] = sum;
	    }	else	{
            Spec_source_inventory[mn][species_id] += 0.5*sum*(delta_t+tran->delta_t);
	    }
      }
  if( print_flag && ProcID == 0 )
    {
      FILE *jfp;
      
      if( (jfp = fopen( filenm, "a")) != NULL )	{
	if (ppvi_type == PPVI_VERBOSE) {
           if(quantity == I_SPECIES_SOURCE)
	      {fprintf(jfp,"   volume= %10.7e \n", Spec_source_inventory[mn][species_id] );}
              else
	      {fprintf(jfp,"   volume= %10.7e \n", sum );}
	}
	if (ppvi_type == PPVI_CSV) {
	  fprintf(jfp,"%10.7e\n", sum );
	}
  	fclose(jfp);
      }
    }

  // Kind of a hack to keep track of the porous liquid inventory for a time-dependent BC
					      
  if(pd->e[R_SHELL_SAT_OPEN]) Porous_liq_inventory = sum; 
  
  return (sum); 
}
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

int
compute_volume_integrand(const int quantity, const int elem,
			 const int species_no, const double *params, 
                         const int num_params,
 			 double *sum, double *J_AC, const int adapt_int_flag,
			 const double time, const double delta_t, double xi[], const Exo_DB *exo)

    /******************************************************************************
     *
     * compute_volume_integrand():
     *
     *
     ******************************************************************************/
{
  double weight = fv->wt;
  double det = bf[pd->ShapeVar]->detJ * fv->h3;
  double det_J = bf[pd->ShapeVar]->detJ;
  double h3 = fv->h3;
  int a, b, p, var,j, dim=pd->Num_Dim, gnn, matIndex, ledof, c, WIM;
  int *n_dof=NULL;
  int dof_map[MDE];

  if(upd->CoordinateSystem == PROJECTED_CARTESIAN)
    WH(-1, "compute_volume_integrand has not been updated for the PROJECTED_CARTESIAN coordinate system.");

  if(upd->CoordinateSystem == CYLINDRICAL ||
     upd->CoordinateSystem == SWIRLING   )
    {
      weight *= 2.*M_PIE;
    }

  WIM = dim;
  if (pd->CoordinateSystem == SWIRLING ||
      pd->CoordinateSystem == PROJECTED_CARTESIAN)
    WIM = WIM+1;

  switch ( quantity )
    {
    case I_VOLUME:

      *sum += weight*det;

      if( J_AC != NULL )
	{
	  for( p=0; p<dim; p++)
	    {
	      var = MESH_DISPLACEMENT1 + p;

	      if( pd->v[var] )
		{

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      J_AC[ ei->gun_list[var][j] ] += weight * ( h3 * bf[pd->ShapeVar]->d_det_J_dm[p][j] +
								 fv->dh3dq[p]*bf[var]->phi[j] * det_J );
		    }
		}
	    }
	}
	       
      break;
    case I_VOLUME_PLANE:
      {
       double dist=0;
       if(num_params < 4)
         {
           WH(-1,"not enough plane parameters for VOL_PLANE\n");
         }
       dist = params[0]*fv->x[0]+params[1]*fv->x[1]
                                +params[2]*fv->x[2]+params[3];
       dist /= sqrt(SQUARE(params[0])+SQUARE(params[1])+SQUARE(params[2]));

       if(dist > 0)
       {
      *sum += weight*det;

      if( J_AC != NULL )
	{
	  for( p=0; p<dim; p++)
	    {
	      var = MESH_DISPLACEMENT1 + p;

	      if( pd->v[var] )
		{

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      J_AC[ ei->gun_list[var][j] ] += weight * ( h3 * bf[pd->ShapeVar]->d_det_J_dm[p][j] +
								 fv->dh3dq[p]*bf[var]->phi[j] * det_J );
		    }
		}
	    }
	}
      }
     }
	       
      break;
    case I_LUB_LOAD:
     {
      n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

      det = fv->sdet; //Different determinant since this is a shell 

      *sum += fv->lubp* weight * det;

      /* clean-up */
      safe_free((void *) n_dof);

      break;
     }
    case I_SHELL_VOLUME:
     {
      n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
      det = fv->sdet; //Different determinant since this is a shell

      dbl H, H_U, dH_U_dtime, H_L, dH_L_dtime;
      dbl dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
      H = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);

      *sum += H* weight * det;

      /* clean-up */
      safe_free((void *) n_dof);

      break;
     }
    case I_SPEED:
     {
      double vsq;
      vsq = 0;
      for( a=0; a<dim;a++) vsq += fv->v[a]*fv->v[a];

      *sum += vsq*weight*det;

      if( J_AC != NULL )
	{
	  for( p=0; p<dim; p++)
	    {
	      var = VELOCITY1 + p ;

	      if( pd->v[var] )
		{

		  for( j=0; j<ei->dof[var]; j++)
		    {
		      J_AC[ ei->gun_list[var][j] ] += weight * ( 2.0*fv->v[p]*bf[var]->phi[j] )*det;
		    }
		}
	    }
	}
      }
      break;
    case I_DISSIP:
      {
	int q;

	double gamma_dot[DIM][DIM];
        double mu;
        double vs[MAX_PDIM][MAX_PDIM];	/* viscous stress */
        double ves[MAX_PDIM][MAX_PDIM];	/* viscoelastic stress */
        int ve_mode;

  if (mp->HeatSourceModel == VISC_DISS )
    {
      *sum += weight*det*visc_diss_heat_source(NULL, mp->u_heat_source);
    }
  else if (mp->HeatSourceModel == VISC_ACOUSTIC )
    {
      *sum += weight*det*visc_diss_heat_source(NULL, mp->u_heat_source);
      *sum += weight*det*visc_diss_acoustic_source(NULL, mp->u_heat_source, mp->len_u_heat_source);
    }
  else
    {
	for( p=0; p<VIM; p++)
	  {
	    for( q=0; q<VIM; q++)
	      {
		gamma_dot[p][q] = fv->grad_v[p][q] + fv->grad_v[q][p];
	      }
	  }

	mu = viscosity( gn, gamma_dot, NULL);
	for( p=0; p<VIM; p++)
	  {
	    for( q=0; q<VIM; q++)
	      {
		vs[p][q] = mu*gamma_dot[p][q] - fv->P*delta(p,q);
	      }
	  }
	/*
	 *  viscoelastic stress tensor
	 *  assume only EVSS_F formulation for now
	 */
  	memset( ves, 0, sizeof(dbl)*DIM*DIM);
  	if ( pd->v[POLYMER_STRESS11] )
    	  {
	    for ( p=0; p<VIM; p++)
               {
	        for ( q=0; q<VIM; q++)
        	   {
	             for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
        	       {
               		ves[p][q] += fv->S[ve_mode][p][q];
               		}
            	   }
              }
    	}


	for(p=0; p<VIM; p++)
	  {
	    for ( q=0; q<WIM; q++)
	      {
		*sum += weight*( (vs[q][p]+ves[q][p])*fv->grad_v[p][q] )*det;
	      }
	  }
	if ( J_AC != NULL)
	  {
	    EH(-1,"Appropriate Jacobian entries for the DISSIP volume integral are not available.\n");
	  }

    }
      }
      break;

   case I_JOULE:
      {
	int q;

        double k;

	/* NOTE: MUST add volume_int sensitivities for use in augmenting conditions */

        if(mp->Elec_ConductivityModel == USER)
	  {
	    usr_electrical_conductivity(mp->u_electrical_conductivity, time);
	 
	    k   = mp->electrical_conductivity;
	  }
	else
	  {
	    k = mp->electrical_conductivity;
	  }

	for ( q=0; q<DIM; q++)
	  {
	    *sum += weight*(k * fv->grad_V[q] * fv->grad_V[q] )*det;
	  }

	if ( J_AC != NULL)
	  {
	    EH(-1,"Appropriate Jacobian entries for the I_JOULE volume integral are not available.\n");
	  }

      }
      break;

    case I_RATE_OF_DEF_II:

      /* Second invariant of the rate of deformation tensor */

      {
	int q;
	
	double gamma_dot[DIM][DIM], II=0.0 ;

	for( p=0; p<VIM; p++)
	  {
	    for( q=0; q<VIM; q++)
	      {
		gamma_dot[p][q] = fv->grad_v[p][q] + fv->grad_v[q][p];
	      }
	  }

	for( p=0; p<VIM; p++)
	  {
	    for( q=0; q<VIM; q++)
	      {
		II   += gamma_dot[p][q]*gamma_dot[q][p];
	      }
	  }


	II = pow ( 0.5*II, 0.5 );

	*sum += weight*det*II;

	if ( J_AC != NULL )
	  {
	    EH(-1,"Jacobian entries for II_GAMMA_DOT Volume integral not implemented.");
	  }
      }

    case I_SPECIES_MASS:
      {
 	double dens_factor, density_tot, untracked_spec=0;
        density_tot = calc_density(mp, FALSE, NULL, 0.0);
   		switch ( mp->Species_Var_Type )
     			{
 			case SPECIES_DENSITY:
 				{ dens_factor=1.; break;	}
 			case SPECIES_UNDEFINED_FORM:
 			case SPECIES_MASS_FRACTION:
 				{ dens_factor=density_tot; break;	}
 			default:
 				{ dens_factor=1.; break;	}
 			}
      if( species_no < pd->Num_Species_Eqn)
        { *sum += weight*det*( dens_factor*fv->c[species_no] ) ; }
      else  {
      switch(mp->Species_Var_Type)   {
      case SPECIES_CONCENTRATION:
        untracked_spec = density_tot;
        for (j=0 ; j < pd->Num_Species_Eqn ; j++)       {
                untracked_spec -= fv->c[j]*mp->molecular_weight[j];
                }
        untracked_spec /= mp->molecular_weight[pd->Num_Species_Eqn];
        break;
      case SPECIES_DENSITY:
        untracked_spec = density_tot;
        for (j=0 ; j < pd->Num_Species_Eqn ; j++)       {
                untracked_spec -= fv->c[j];
                }
        break;
      case SPECIES_MASS_FRACTION:
      case SPECIES_UNDEFINED_FORM:
        untracked_spec = 1.0;
        for (j=0 ; j < pd->Num_Species_Eqn ; j++)       {
                untracked_spec -= fv->c[j];
                }
        }
        *sum += weight*det*( dens_factor*untracked_spec ) ; 
      }

	if( J_AC != NULL )
	  {
	    for( a=0; a<dim ; a++)
	      {
		var = MESH_DISPLACEMENT1 + a;

		if( pd->v[var] )
		  {
		    for( j=0 ; j<ei->dof[var]; j++)
		      {
	    
			J_AC[ ei->gun_list[var][j] ] += weight*  ( h3 * bf[pd->ShapeVar]->d_det_J_dm[a][j] +
 					    fv->dh3dq[a]*bf[var]->phi[j] * det_J )*dens_factor*fv->c[species_no];
		      }
		  }
	      }

	    var = MASS_FRACTION;
	    if (pd->v[var]) {
	      for (j = 0; j < ei->dof[var]; j++) {
		/*
		 * Find the material index for the current
		 * local variable degree of freedom
		 * (can't just query gun_list for MASS_FRACTION
		 *  unknowns -> have to do a lookup)
		 */		
		gnn = ei->gnn_list[var][j];
		ledof = ei->lvdof_to_ledof[var][j];
		matIndex = ei->matID_ledof[ledof];
		c = Index_Solution(gnn, var, species_no, 0, matIndex);
 		J_AC[c] += weight * det * dens_factor * bf[var]->phi[j];
	      }
	    }
	  }
      }
      break;

    case I_SPECIES_SOURCE:
      {
        struct Species_Conservation_Terms s_terms; 
        int w1,i,ie,ldof,eqn=MASS_FRACTION;
        zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
        get_continuous_species_terms(&s_terms, time, tran->theta, delta_t, NULL);
        if(time > tran->init_time+delta_t)
 	   {
            *sum += weight*det*(-s_terms.MassSource[species_no])
                    *mp->specific_volume[species_no]/mp->specific_volume[pd->Num_Species_Eqn];
           }
  if (efv->ev) {
        if( efv->i[species_no] != I_TABLE)
          {       
           for (i = 0; i < ei->num_local_nodes; i++) {
                ie = Proc_Elem_Connect[ei->iconnect_ptr + i];
		ldof=ei->ln_to_dof[eqn][i];
                if(ldof >= 0 )
                   {
                    if(time <= tran->init_time+delta_t)
 	                 {
                          *sum += weight*det*bf[eqn]->phi[ldof]*
                          efv->ext_fld_ndl_val[species_no][ie];
                         }
                    efv->ext_fld_ndl_val[species_no][ie] += 
                              bf[eqn]->phi[ldof]*weight*det*
                              0.5*(delta_t+tran->delta_t)
                              *(-s_terms.MassSource[species_no])
                    *mp->specific_volume[species_no]/mp->specific_volume[pd->Num_Species_Eqn];
                    if(species_no == 0)	
                          {
                          Spec_source_lumped_mass[ie] += 
                               bf[eqn]->phi[ldof]*weight*det;
                          }
                   }
                } 
          }
  }

	if( J_AC != NULL )
	  {
	    for( a=0; a<dim ; a++)
	      {
		var = MESH_DISPLACEMENT1 + a;

		if( pd->v[var] )
		  {
		    for( j=0 ; j<ei->dof[var]; j++)
		      {
	    
			J_AC[ ei->gun_list[var][j] ] += weight*  ( h3 * bf[pd->ShapeVar]->d_det_J_dm[a][j] +
 			fv->dh3dq[a]*bf[var]->phi[j] * det_J )*s_terms.MassSource[species_no];
		      }
		  }
	      }

	    var = MASS_FRACTION;
	    if (pd->v[var]) {
	    for ( w1=0; w1<pd->Num_Species_Eqn; w1++)  {
	      for (j = 0; j < ei->dof[var]; j++) {
		/*
		 * Find the material index for the current
		 * local variable degree of freedom
		 * (can't just query gun_list for MASS_FRACTION
		 *  unknowns -> have to do a lookup)
		 */		
		gnn = ei->gnn_list[var][j];
		ledof = ei->lvdof_to_ledof[var][j];
		matIndex = ei->matID_ledof[ledof];
		c = Index_Solution(gnn, var, species_no, 0, matIndex);
 		J_AC[c] += weight * det * s_terms.d_MassSource_dc[species_no][w1][j];
	      }
	     }
	    }
	  }
      }
      break;
    case I_MOMX:
    case I_MOMY:
    case I_MOMZ:
      {
	int dir;
	double rho;

	dir = ( quantity == I_MOMX ? 0 :
	      ( quantity == I_MOMY ? 1 : 
	      ( quantity == I_MOMZ ? 2 : -1) ) ) ;

	EH(dir, "PANIC in compute_volume_integrand \n");

	rho = density( NULL, time );

	*sum += weight*det*rho*fv->v[dir];

	if ( J_AC != NULL)
	  {
	    EH(-1,"Appropriate Jacobian entries for the MOM volume integral are not available.\n");
	  }

      }
      break;
    case I_HEAT_ENERGY:
      {
      double rho, Cp;
      rho = density( NULL, time );
      Cp = heat_capacity( NULL, time );

      *sum += rho*Cp*fv->T*weight*det;

      if( J_AC != NULL )
      {
        for( p=0; p<dim; p++)
          {
            var = MESH_DISPLACEMENT1 + p;

            if( pd->v[var] )
              {

                for( j=0; j<ei->dof[var]; j++)
                  {
                    J_AC[ ei->gun_list[var][j] ] += rho*Cp*fv->T * weight *
                              ( h3 * bf[pd->ShapeVar]->d_det_J_dm[p][j] +
                              fv->dh3dq[p]*bf[var]->phi[j] * det_J );
                  }
              }
          }
            var = TEMPERATURE;

            if( pd->v[var] )
              {

                for( j=0; j<ei->dof[var]; j++)
                  {
                    J_AC[ ei->gun_list[var][j] ] +=
                              rho*Cp*bf[var]->phi[j]*weight*det;
                  }
              }
      }
/*      EH(-1,"This volumetric integral not yet implemented \n");  */
      }
      break;

    case I_KINETIC_ENERGY:
      {

      *sum += fv->T*weight*det;

      if( J_AC != NULL )
      {
        for( p=0; p<dim; p++)
          {
            var = MESH_DISPLACEMENT1 + p;

            if( pd->v[var] )
              {

                for( j=0; j<ei->dof[var]; j++)
                  {
                    J_AC[ ei->gun_list[var][j] ] += fv->T * weight *
                              ( h3 * bf[pd->ShapeVar]->d_det_J_dm[p][j] +
                              fv->dh3dq[p]*bf[var]->phi[j] * det_J );
                  }
              }
          }
            var = TEMPERATURE;

            if( pd->v[var] )
              {

                for( j=0; j<ei->dof[var]; j++)
                  {
                    J_AC[ ei->gun_list[var][j] ] +=
                              bf[var]->phi[j]*weight*det;
                  }
              }
      }
/*      EH(-1,"This volumetric integral not yet implemented \n");  */
      }
      break;

    case I_TRACE:
      {
	int q;

	double gamma_dot[DIM][DIM];
	double mu;
        VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
        VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;



	if( vn->ConstitutiveEquation != NOPOLYMER )
	  {
	    EH(-1,"Computation of stress integral for POLYMER models currently no supported.\n");
	  }
			       

	for( p=0; p<VIM; p++)
	  {
	    for( q=0; q<VIM; q++)
	      {
		gamma_dot[p][q] = fv->grad_v[p][q] + fv->grad_v[q][p];
	      }
	  }

	mu = viscosity( gn, gamma_dot, d_mu );

	*sum += weight*det*( -3.0*fv->P );

	if ( J_AC != NULL )
	  {
	    double dterm, dterm1;

	    var = PRESSURE;
	    
	    if( pd->v[var] )
	      {
		for ( j=0; j < ei->dof[var]; j++ )
		  {
		    J_AC[ ei->gun_list[var][j] ] += weight*det* ( -3.0* bf[var]->phi[j]);
		  }
	      }

	    
	    for( a=0; a<dim ; a++)
	      {
		var = VELOCITY1 + a;

		if( pd->v[var] )
		  {
		    for( j=0; j<ei->dof[var]; j++)
		      {
			for ( b=0; b<dim; b++)
			  {

			    J_AC [ ei->gun_list[var][j] ] += weight * det * ( d_mu->v[a][j] * 2.0*fv->grad_v[a][a] +
									     mu * 2.0*bf[var]->grad_phi_e[j][a][b][b] ) ;
			  }
		      }
		  }
	      }

	    
	    for( a=0; a<dim ; a++)
	      {
		var = MESH_DISPLACEMENT1 + a;

		if( pd->v[var] )
		  {
		    for( j=0 ; j<ei->dof[var]; j++)
		      {
			dterm = dterm1 = 0.0;
			for ( b=0; b<dim; b++ )
			  {
			    dterm += weight*  ( h3 * bf[pd->ShapeVar]->d_det_J_dm[a][j] +
					    fv->dh3dq[a]*bf[var]->phi[j] * det_J )*( -fv->P + mu*2.0*fv->grad_v[b][b] ) ;

			    dterm1 += weight * det * ( d_mu->X[a][j]*2.0*fv->grad_v[b][b] +
						   mu* 2.0* fv->d_grad_v_dmesh[b][b][a][j] );
			  }
			
			J_AC[ ei->gun_list[var][j] ] +=  dterm + dterm1;
		      }
		  }
	      }
	  }	  
      }
      
      break;

    case I_POS_FILL:
    case I_NEG_FILL:
     {
	double alpha, height = 1.0;
        double H_U, dH_U_dtime, H_L, dH_L_dtime;
        double dH_U_dX[DIM],dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
	if(pd->e[R_LUBP]) height = height_function_model
                          (&H_U, &dH_U_dtime, &H_L, &dH_L_dtime,
			   dH_U_dX, dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);
      if( adapt_int_flag)
      {
 	double dwt_dF=0;
 	*sum += height * weight * det;
 
 
 	if ( J_AC != NULL )
 	  {
 	    var = LS;
 	    for( j=0 ; j < ei->dof[var]; j++)
 	      {
 		J_AC[ ei->gun_list[var][j]  ] += dwt_dF*height*det;
 	      }
 	  } 
      }
      else
      {
	double H = 0.0, sign;

	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];

	load_lsi(alpha);

	H = quantity == I_POS_FILL ? lsi->H : ( 1.0 - lsi->H ) ;
	
	*sum += height * weight * det * H;

	if ( J_AC != NULL )
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();

	    var = LS;
	    for( j=0 ; j < ei->dof[var]; j++)
	      {
		sign = quantity == I_POS_FILL ? +1.0 : -1.0;
		J_AC[ ei->gun_list[var][j]  ] += height * weight * 
                                    det * ( sign * lsi->d_H_dF[j] );
	      }
	  } 
      }
     }
      break;
    case I_POS_VOLPLANE:
    case I_NEG_VOLPLANE:
     {
	double alpha, height = 1.0;
        double dist=0;
       if(num_params < 5)
         {
           WH(-1,"not enough plane parameters for POS-NEG_PLANE_FILL\n");
         }
       dist = params[1]*fv->x[0]+params[2]*fv->x[1]
                                +params[3]*fv->x[2]+params[4];
       dist /= sqrt(SQUARE(params[1])+SQUARE(params[2])+SQUARE(params[3]));
       if(dist > 0)
       {
      if( adapt_int_flag)
      {
 	double dwt_dF=0;
 	*sum += height * weight * det;
 
 
 	if ( J_AC != NULL )
 	  {
 	    var = LS;
 	    for( j=0 ; j < ei->dof[var]; j++)
 	      {
 		J_AC[ ei->gun_list[var][j]  ] += dwt_dF*height*det;
 	      }
 	  } 
      }
      else
      {
	double H = 0.0, sign;

	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];

	load_lsi(alpha);

	H = quantity == I_POS_VOLPLANE ? lsi->H : ( 1.0 - lsi->H ) ;
	
	*sum += height * weight * det * H;

	if ( J_AC != NULL )
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();

	    var = LS;
	    for( j=0 ; j < ei->dof[var]; j++)
	      {
		sign = quantity == I_POS_VOLPLANE ? +1.0 : -1.0;
		J_AC[ ei->gun_list[var][j]  ] += height * weight * 
                                    det * ( sign * lsi->d_H_dF[j] );
	      }
	  } 
      }
      }
     }
      break;

    case I_MAG_GRAD_FILL_ERROR:
      {
	double alpha;

	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];

	load_lsi(alpha);

	*sum += weight * det * lsi->delta * (lsi->gfmag - 1.) * (lsi->gfmag - 1.);

	if ( J_AC != NULL )
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();

	    var = ls->var;
	    for( j=0 ; j < ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * (lsi->gfmag - 1.) *
                                                 (lsi->delta * 2.*lsi->d_gfmag_dF[j] +
                                                  (lsi->gfmag - 1.) * lsi->d_delta_dF[j]);
	      }
	  }
      }
      break;
      
    case I_LS_ARC_LENGTH:
      {
	double alpha;

	alpha = (params == NULL ) ? ls->Length_Scale : 2.0*params[0];
	
	load_lsi( alpha );

	*sum += weight*det*lsi->delta;

	if ( J_AC != NULL )
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();

	    var = ls->var;
	    for( j=0 ; j < ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * lsi->d_delta_dF[j];
	      }
	  }
      }
      break;

    case I_POS_VX:
    case I_POS_VY:
    case I_POS_VZ:
    case I_NEG_VX:
    case I_NEG_VY:
    case I_NEG_VZ:

      /* 
	 These cases are for computing the integral of velocity in one
	 of the three direction over just the negative or positive
	 phase as defined by the level set function
      */

      {
	int dir;
	double H, sign, alpha;
	
	dir = ( quantity == I_NEG_VX ? 0 :
		( quantity == I_NEG_VY ? 1 :
		  ( quantity == I_NEG_VZ ? 2 : -1 ) ) );

	if ( dir == -1 )
	  {
	    dir = ( quantity == I_POS_VX ? 0 :
		    ( quantity == I_POS_VY ? 1 :
		      ( quantity == I_POS_VZ ? 2 : -1 ) ) );
	    sign = +1.0;
	  }
	else
	  {
	    /* OK, the negative phase was selected. */
	    sign = -1.0;
	  }

	if ( dir == -1 )
	  EH(-1, "compute_volume_integrand(): ACK! I'm confused!\n");
	  
	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];
	
	load_lsi(alpha);

	/* Positive or negative H? */
	H = sign > 0.0 ? lsi->H : 1.0 - lsi->H;
	
	*sum += weight * det * H * fv->v[dir];
	
	if( J_AC != NULL )        
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();
	    
	    var = LS ;
	    for( j=0; j<ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * (sign * lsi->d_H_dF[j]) * fv->v[dir];
	      }
	    
	    var = VELOCITY1 + dir;
	    for( j=0; j<ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * H * bf[var]->phi[j];
	      }
	  }
      }
      break;
    case I_POS_CENTER_X:
    case I_POS_CENTER_Y:
    case I_POS_CENTER_Z:
    case I_NEG_CENTER_X:
    case I_NEG_CENTER_Y:
    case I_NEG_CENTER_Z:

      /* 
	 These cases are for computing the integral of velocity in one
	 of the three direction over just the negative or positive
	 phase as defined by the level set function
      */

      {
	int dir;
	double H, sign, alpha;
	
	dir = ( quantity == I_NEG_CENTER_X ? 0 :
		( quantity == I_NEG_CENTER_Y ? 1 :
		  ( quantity == I_NEG_CENTER_Z ? 2 : -1 ) ) );

	if ( dir == -1 )
	  {
	    dir = ( quantity == I_POS_CENTER_X ? 0 :
		    ( quantity == I_POS_CENTER_Y ? 1 :
		      ( quantity == I_POS_CENTER_Z ? 2 : -1 ) ) );
	    sign = +1.0;
	  }
	else
	  {
	    /* OK, the negative phase was selected. */
	    sign = -1.0;
	  }

	if ( dir == -1 )
	  EH(-1, "compute_volume_integrand(): ACK! I'm confused!\n");
	  
	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];
	
	load_lsi(alpha);

	/* Positive or negative H? */
	H = sign > 0.0 ? lsi->H : 1.0 - lsi->H;
	
	*sum += weight * det * H * fv->x[dir];
	
	if( J_AC != NULL )        
	  {

	    /* Do this even for uncoupled fill problems. */
	    load_lsi_derivs();
	    
	    var = LS ;
	    for( j=0; j<ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * (sign * lsi->d_H_dF[j]) * fv->x[dir];
	      }
	    
            /* MOVING MESH NOT HANDLED YET */
            /*
	    var = VELOCITY1 + dir;
	    for( j=0; j<ei->dof[var]; j++)
	      {
		J_AC[ ei->gun_list[var][j]  ] += weight * det * H * bf[var]->phi[j];
	      }
            */
	  }
      }
      break;
    case I_POROUS_LIQUID_INV:
      {
	if(pd->e[R_SHELL_SAT_OPEN])
	  {
	    n_dof = (int *)array_alloc (1, MAX_VARIABLE_TYPES, sizeof(int));
	    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

	    det = fv->sdet; //Different determinant since this is a shell 

	    /* clean-up */
	    safe_free((void *) n_dof);
	  }

	*sum += weight*det*pmv->bulk_density[0];

	if ( J_AC != NULL)
	  {
	    EH(-1,"Appropriate Jacobian entries for the Porous Liquid INV volume integral are not available.\n");
	  }

      }
      break;
	  
	case I_ELOADX:
    case I_ELOADY:
    case I_ELOADZ:
      {
	int dir;
    double factor = 1.0;
	dir = ( quantity == I_ELOADX ? 0 :
	      ( quantity == I_ELOADY ? 1 :
          ( quantity == I_ELOADZ ? 2 : -1 ) ) );


	if( pd->e[R_MASS] ) factor = fv->c[species_no];

    if ( params != NULL ) factor *= params[0];

	if( pd->e[R_EFIELD1] )
		{	
           *sum += factor*weight*det*fv->E_field[dir];	
        }
	else if( pd->e[R_POTENTIAL] )
		{
           *sum += factor*weight*det*(-fv->grad_V[dir]);
		}		
	else
	  EH(-1,"You need to include the potential and/or electric field variables to compute this integrated quantity\n/");

      }

    default:
      break;
    }
  
  return ( 1 );
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/* evaluate_global_flux ()  -- This function evaluates the surface integral of the specified
 *                             quantity over the entire exterior surface of the specified
 *                             element block. 
 * Author: T.A. Baer
 * Date:   May 2002
 */


double
evaluate_global_flux (const Exo_DB *exo,
		      const Dpi *dpi,
		      const int quantity,
		      const int blk_id,
		      const int species_id,
		      const double *params,
		            double *J_AC, 
		      const double x[],
		      const dbl time_value,
		      const int print_flag )
{
  int j;
  int eb,e_start,e_end, elem;


  int err=0;
  int ip, ip_total;
  double sum = 0.0;

#ifdef PARALLEL
  double sum0=0.0;
  double delta_sum=0.0;
  double proc_sum=0.0;
  double global_sum=0.0;

#endif


  map_mat_index(blk_id);
  if( ( eb = in_list(blk_id, 0, exo->num_elem_blocks, exo->eb_id) ) != -1 )
    {

      e_start = exo->eb_ptr[eb];
      e_end   = exo->eb_ptr[eb+1];

      for( elem= e_start;   elem < e_end ;  elem++)
        {
          double xi[3], s,t,u;
	  int exterior_faces[8];
	  int num_exterior_faces;
	  int id_side;

	  int elem_type;
	  int id_local_elem_coord[MAX_NODES_PER_SIDE];
	  int nodes_per_side;
  

	  if (( num_exterior_faces = get_exterior_faces( elem, exterior_faces, exo, dpi ) ) > 0 ) /* and if it has exterior */
	    {

	      err = load_elem_dofptr(elem, (Exo_DB*) exo,
				     (dbl *) x, (dbl *) x, (dbl *) x,
				     (dbl *) x, (dbl *) x, 0);
	      EH(err, "load_elem_dofptr");
	      
	      elem_type = ei->ielem_type;
	      ip_total = elem_info( NQUAD_SURF, ei->ielem_type );

	      
	      for( j=0 ; j < num_exterior_faces; j++ )
		{
		  id_side = exterior_faces[j] + 1;	      

		  get_side_info ( elem_type, id_side, &nodes_per_side, id_local_elem_coord );

		  for( ip = 0 ; ip < ip_total; ip++ )
		    {

		      find_surf_st( ip, elem_type, id_side, ei->ielem_dim, xi, &s, &t, &u);

		      fv->wt = Gq_surf_weight (ip, elem_type);    
		  
		      err = load_basis_functions( xi, bfd);
		      EH( err, "problem from load_basis_functions");
	    
		      err = beer_belly();
		      EH( err, "beer_belly");
	    
		      /* precalculate variables at  current integration pt.*/
		      err = load_fv();
		      EH( err, "load_fv");
		  
		      err = load_bf_grad();
		      EH( err, "load_bf_grad");

		      surface_determinant_and_normal (elem, ei->iconnect_ptr, 
						      ei->num_local_nodes, 
						      ei->ielem_dim - 1,  
						      id_side,
						      nodes_per_side,
						      id_local_elem_coord );


		      /*
		       * Load up physical space gradients of field variables at this
		       * Gauss point.
		       */

		      err = load_fv_grads();
		      EH( err, "load_fv_grads");

		      compute_surface_integrand( quantity, elem, species_id, 
						params, &sum, J_AC );

#ifdef PARALLEL
		      delta_sum = sum - sum0;
		      
		      if( Num_Proc > 1 && dpi->elem_owner[ ei->ielem ] == ProcID ) 
			{
			  proc_sum += delta_sum;
			}

		      sum0 = sum;
#endif


		    }
		}
	    }
	}
    }

#ifdef PARALLEL
  if( Num_Proc > 1 ) 
    {

	MPI_Allreduce( &proc_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 
		       MPI_COMM_WORLD);

	sum = global_sum;
    }
#endif

return ( sum ) ;
}
		      

void
compute_surface_integrand (const int quantity,
			   int elem, 
			   const int species_id,
			   const double *params, 
			   double *sum,
			   double J_AC[] )
{
  double weight = fv->wt;
  double det;
  int a, dim=pd->Num_Dim;

  det = fv->sdet*fv->h3;

  if(upd->CoordinateSystem == CYLINDRICAL ||
     upd->CoordinateSystem == SWIRLING   )
    {
      weight *= 2.*M_PIE;
      det = fv->sdet;
    }

  switch (quantity) 
    {
 
    case VOLUME_FLUX:
      for( a = 0 ; a<dim; a++ )
	*sum += weight*fv->v[a]*fv->snormal[a]*det;
      
      if ( J_AC != NULL )
	{
	  EH(-1, "Global surface sentivities are on my to-do list.") ;
	}
      break;
    case NEG_LS_FLUX:
    case POS_LS_FLUX:

      {
	double H = 0.0;
	double alpha;

	alpha = params == NULL ? ls->Length_Scale : 2.0*params[0];

	load_lsi(alpha);

	H =  ( quantity == POS_LS_FLUX ? lsi->H : ( 1.0 - lsi->H ) ) ;
	
	for( a=0; a<dim; a++)
	  
	  *sum += weight * det * H*fv->v[a]*fv->snormal[a] ;

	if ( J_AC != NULL )
	  {
	    EH(-1, "Global surface sentivities are on my to-do list.") ;
	  }
      }
      break;

    default:
      break;
    }
}




/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/* evaluate_flux_sens() -- compute flux/force sensitivities along a sideset, print out
 *
 * Author:          P. R. Schunk
 * Date:            12 Jan 1996
 *
 * (o) Replaced EXODUS II I/O calls with references to data already read.
 *
 * (o) Considered distributed processing case where this SSID may appear
 *     on other processors and incompletely here.
 *
 *    Based on a clone of evaluate_flux on 6/25/99.  Field Variable sensitivities
 *	and gradient sensivities should be included but not sensitivities to
 *	mesh (i.e. sensitivity of normals, intervals of integration are
 *	ignored).  Also,  sensitivities of material properties are ignored and
 *	force sensitivities in a solid material are not done. - RBS
 */

double
evaluate_flux_sens(const Exo_DB *exo, /* ptr to basic exodus ii mesh information */
	      const Dpi *dpi, /* distributed processing info */
	      const int side_set_id, /* on which SSID to evaluate flux */
	      const int quantity, /* to pick HEAT_FLUX, FORCE_NORMAL, etc. */
              const char *qtity_str, /* quantity string */
	      const int mat_id,	/* material identification */
	      const int species_id, /* species identification */
	      const int sens_type,  /*  sensitivity type */
	      const int sens_id,  /* sensitivity id */
	      const int sens_flt,	/*  sensitivity float number */
	      const int sens_flt2,	/*  sensitivity float number for UM */
	      const int vector_id, /* sensitivity id */
	      const char *filenm, /* File name pointer */
              const int profile_flag,  /*  flux sens print flag  */
	      const double x[],	/* solution vector */
	      const double xdot[],	/* solution vector */
	      double **x_sens_p,	/* sensitivity vector */
	      const double delta_t, /* time-step size */
              const double time_value, /* current time */
              const int print_flag )   /*  printing control flag */
{
  int j;			/* local index loop counter                 */
  int i;			/* Index for the local node number - row    */
  int ip;
  int a,b;
  int mn;
  int var;
  double wt,weight;
  int err;                    /* temp variable to hold diagnostic flags. */
  int ielem=-1;			/* element number                            */
  int ip_total;       /* total number of gauss points */
  int ielem_type;          /* element type */
  int num_local_nodes;
  int iconnect_ptr;
  int ielem_dim;
  int current_id;
  int v, dofs;
  int nset_id, sset_id;

  /* 
   * Variables for vicosity and derivative 
   */
  dbl gamma[DIM][DIM];
  dbl gamma_sens[DIM][DIM];
  dbl mu = 0.0;

  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;  /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  dbl mu_sens;

  dbl q[DIM],dq_gradT[DIM][DIM],dq_dX[DIM][DIM];

  dbl Cp = 0.0;                               /* Heat capacity. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT d_Cp_struct; /* Heat capacity dependence. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp = &d_Cp_struct;

  double rho=0; /*  density variables */
  dbl e_theta[DIM] = {0.,0.,1.};  /* torque w.r.t theta only*/

  double vs_sens[MAX_PDIM][MAX_PDIM];	/* viscous stress */
  double vs[MAX_PDIM][MAX_PDIM];	

  double ves_sens[MAX_PDIM][MAX_PDIM];	/* viscoelastic stress */
  double ves[MAX_PDIM][MAX_PDIM];	
  int ve_mode;
  int v_s[MAX_MODES][DIM][DIM];

  double TT[MAX_PDIM][MAX_PDIM];   /**  solid stresses  **/
  dbl dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  dbl dTT_dp_liq[DIM][DIM][MDE];
  dbl dTT_dp_gas[DIM][DIM][MDE];
  dbl dTT_dporosity[DIM][DIM][MDE];
  dbl dTT_dsink_mass[DIM][DIM][MDE];
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];
  double TT_sens[MAX_PDIM][MAX_PDIM];
  double elast_modulus;
  int dim, velodim;

  double es[MAX_PDIM][MAX_PDIM], efield[MAX_PDIM];	/* electric stress */
  double efield_sqr;				/* efield magnitude squared  */
  double es_sens[MAX_PDIM][MAX_PDIM], efield_sens[MAX_PDIM];/* electric stress */
  double efield_sqr_sens;		/* efield magnitude squared  */

  double dsigma_dx[DIM][MDE];		/* surface tension terms */
  double dsigmadT, dsigmadC[MAX_CONC], sigma_sens;
  double x_dot[MAX_PDIM];   /*    moving mesh stuff  */
  double normal_sens[MAX_PDIM],tangent_sens[2][MAX_PDIM];

  double s, t, u;			/* Gaussian-quadrature point locations       */
  double det;
  double local_flux = 0.0;
  double local_flux_conv = 0.0;
  double local_q;
  double local_qconv;
  double local_area = 0.0;
#ifdef PARALLEL
  double local_flux0 = 0.0;     /* old flux sum */
  double local_flux_conv0 = 0.0;/* old convective flux sum */
  double local_area0= 0.0;      /* old area sum */
  double proc_flux=0.0;         /* current flux sum on this proc */
  double proc_flux_conv=0.0;    /* current convective flux sum on this proc */
  double proc_area=0.0;         /* current area sum on this proc */
  double delta_flux=0.0;        /* increment of flux             */
  double delta_flux_conv=0.0;   /* increment of convective flux  */
  double delta_area=0.0;        /* increment of area */
  double global_flux=0.0;       /* flux sum over all procs */
  double global_flux_conv=0.0;  /* convective flux sum over all procs */
  double global_area=0.0;       /* area sum over all procs */
#endif

  double xi[DIM];             /* Local element coordinates of Gauss point. */

  int num_side_in_set;        /* returned number of sides(faces or edges) is  */
  int num_elem_in_set, num_node_in_set, num_nodes_on_side;
                              /* in the side set                              */
  int num_dist_fact_in_set;   /* returned number of distribution factors in   */
  int *elem_list;


  int c;
  double Tract[DIM], Torque[DIM], local_Torque[DIM], Tract_sens[DIM];

  SGRID *element_search_grid=	NULL;	
  double surface_centroid[DIM];

  int id_side;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];

#ifndef PARALLEL
  static const char yo[] = "evaluate_flux_sens";
#endif

/*  variables for adaptive quadrature weight calculation   */
 double ls_F[9],ad_wt[9];
 int ierr = 0, wt_type;
 const int Jac_state = af->Assemble_Jacobian;

  memset( Torque, 0, sizeof(double)*DIM );
  memset(TT_sens, 0, sizeof(double) * MAX_PDIM * MAX_PDIM);

  /* load eqn and variable number in tensor form */
  err = stress_eqn_pointer(v_s);
  af->Assemble_Jacobian = TRUE;

  /* first right time stamp or run stamp to separate the sets */

        if (print_flag && ProcID == 0) {
        FILE  *jfp;
        if(  (filenm != NULL) && ( (jfp=fopen(filenm,"a")) != NULL))
                {
      fprintf(jfp,"Time/iteration = %e \n", time_value);
      fprintf(jfp," %s  Side_Set  %d Block  %d Species %d\n", qtity_str,side_set_id,mat_id,species_id);
            if(sens_type == 1)
              {
                fprintf(jfp," Sensitivity wrt  BC  %d  Float  %d\n",sens_id,sens_flt);
              }
            else if(sens_type == 2)
              {
                fprintf(jfp," Sensitivity wrt  MT  %d  Prop.  %d\n",sens_id,sens_flt);
              }
            else if(sens_type == 3)
              {
                fprintf(jfp," Sensitivity wrt  AC  %d  Float  %d\n",sens_id,sens_flt);
              }
            else if(sens_type == 4)
              {
                fprintf(jfp," Sensitivity wrt  UM  %d  Prop.  %d  %d\n",sens_id,sens_flt,sens_flt2);
              }
            else if(sens_type == 5)
              {
                fprintf(jfp," Sensitivity wrt  UF  %d  Float  %d\n",sens_id,sens_flt);
              }
            else if(sens_type == 6)
              {
                fprintf(jfp," Sensitivity wrt  AN  %d  Float  %d\n",sens_id,sens_flt);
              }
      fflush(jfp);
      fclose(jfp);
                }
        }

  current_id = in_list(side_set_id, 0, exo->num_side_sets, exo->ss_id);

  /*
   * Not finding this side set in the mesh is immediately serious if we're
   * running serial. Parallel execution is OK with this particular processor
   * not finding it, though, just so long as *some* processor(s) find this
   * SSID.
   */

  if ( current_id != -1 )
    {
      num_side_in_set      = exo->ss_num_sides[current_id];

      num_dist_fact_in_set = exo->ss_num_distfacts[current_id];

      num_elem_in_set      = num_side_in_set;

      num_node_in_set      = num_dist_fact_in_set;

      num_nodes_on_side    = num_node_in_set/num_elem_in_set; /* Well... */

      elem_list            = &exo->ss_elem_list[exo->ss_elem_index[current_id]];

      /* Now start element sweep */

      for (i=0; i< num_elem_in_set; i++)
	{         
	  ei->ielem = elem_list[i];

	  mn = find_mat_number(ei->ielem, exo);


	/**    only do integration if the requested material number  **/

	if ( mat_id == -1 || 
             (mn == map_mat_index(mat_id) && dpi->elem_owner[ elem_list[i] ] == ProcID))
	{
	  /*
	   * Yes, "x" gets recycled like garbage. Fortunately, this 
	   * routine should not write onto "x"...
	   */

	  err = load_elem_dofptr(elem_list[i], (Exo_DB*) exo,
				 (dbl *) x, (dbl *) &x_sens_p[vector_id][0],
				 (dbl *) xdot, (dbl *) xdot, (dbl *) x, 0);
	  EH(err, "load_elem_dofptr");

	  err = bf_mp_init(pd);
	  EH(err, "bf_mp_init");
  
	  iconnect_ptr        = ei->iconnect_ptr;
	  ielem_type          = ei->ielem_type;
	  ip_total            = elem_info(NQUAD_SURF, ielem_type);
	  num_local_nodes     = ei->num_local_nodes;
	  ielem_dim           = ei->ielem_dim;

	 
	  id_side = find_id_side (ei->ielem, num_nodes_on_side,
				  &exo->ss_node_list[current_id]
				                    [num_nodes_on_side*i],
				  id_local_elem_coord, exo);

	      /* Calculates the ID side correctly for tets */
	      if ( ielem_type == TRILINEAR_TET ) id_side = find_id_side_SS(ei->ielem,
 current_id, exo);

		  if( ls != NULL && ielem_dim == 2  &&
			(quantity == POS_LS_FLUX || quantity == NEG_LS_FLUX ||
				quantity == DELTA || quantity == LS_DCA ))
		  {
                    if ( ls->var != LS )
			WH(-1,"Level-set variable is not LS!");
			  switch (id_side)
			  {
				  case 1:
					  ls_F[0]=*esp->F[0];
					  ls_F[1]=*esp->F[1];
					  ls_F[2]=*esp->F[4];
					  break;
				  case 2:
					  ls_F[0]=*esp->F[1];
					  ls_F[1]=*esp->F[2];
					  ls_F[2]=*esp->F[5];
					  break;
				  case 3:
					  ls_F[0]=*esp->F[3];
					  ls_F[1]=*esp->F[2];
					  ls_F[2]=*esp->F[6];
					  break;
				  case 4:
					  ls_F[0]=*esp->F[0];
					  ls_F[1]=*esp->F[3];
					  ls_F[2]=*esp->F[7];
					  break;
			  }
			if(species_id > 0)
			{
			  var=species_id-1;
			  switch (id_side)
			  {
				  case 1:
					  ls_F[0]=*esp->pF[var][0];
					  ls_F[1]=*esp->pF[var][1];
					  ls_F[2]=*esp->pF[var][4];
					  break;
				  case 2:
					  ls_F[0]=*esp->pF[var][1];
					  ls_F[1]=*esp->pF[var][2];
					  ls_F[2]=*esp->pF[var][5];
					  break;
				  case 3:
					  ls_F[0]=*esp->pF[var][3];
					  ls_F[1]=*esp->pF[var][2];
					  ls_F[2]=*esp->pF[var][6];
					  break;
				  case 4:
					  ls_F[0]=*esp->pF[var][0];
					  ls_F[1]=*esp->pF[var][3];
					  ls_F[2]=*esp->pF[var][7];
					  break;
			  }
			}
			  
			  switch( quantity )
			  {
                                  case POS_LS_FLUX:
					  wt_type = 2;
					  break;
				  case NEG_LS_FLUX:
					  wt_type = 2;
					  ls_F[0] = -ls_F[0];
					  ls_F[1] = -ls_F[1];
					  ls_F[2] = -ls_F[2];
					  break;
                                  case DELTA:
                                  case LS_DCA:
					  wt_type = 3;
					  break;
                                  default:
                                          wt_type = 2;
                                          break;
			  }

			  ierr = adaptive_weight( ad_wt, ip_total, ielem_dim-1,
                                  ls_F, 0.0, wt_type, ei ->ielem_type);
			  if(ierr == -1)printf("adaptive wt problem %d %g %g %g\n",ierr,ls_F[0],ls_F[1],ls_F[2]);
		  }
		  
		  if( ls!=NULL && ls->elem_overlap_state && ls->Integration_Depth > 0)
		  {
			  Subgrid_Int.active = TRUE;
			  
			  get_subgrid_integration_pts ( Subgrid_Tree, &element_search_grid,
		&Subgrid_Int.s, &Subgrid_Int.wt, ls->Length_Scale );
			  
			  find_surf_center_st ( ielem_type, id_side, ielem_dim, surface_centroid, &s, &t );
			  ip_total = gather_surface_subgrid_integration_pts( element_search_grid, id_side, 
			surface_centroid, Subgrid_Int.s, Subgrid_Int.wt, 0 );
			  
		  }
		  



	  /* Surface integration over element */
	  
	  for (ip = 0; ip < ip_total; ip++)
	    {
                  /*   zero local fluxes  */
                  local_q = local_qconv = 0.;
		  memset(local_Torque, 0, DIM*sizeof(double) );
			  
		  if( ls != NULL && ls->elem_overlap_state && ls->Integration_Depth > 0 )
		  {
			  ls->Elem_Sign = 0;
			  xi[0] = Subgrid_Int.s[ip][0];
			  xi[1] = Subgrid_Int.s[ip][1];
			  xi[2] = Subgrid_Int.s[ip][2];
				  
			  s = 1.e30; 
			  t = 1.e30;
			  wt = Subgrid_Int.wt[ip];
		  }
		  else
		  {

	      /* find the quadrature point locations for current ip */
	      find_surf_st (ip, ielem_type, id_side, ielem_dim, xi, &s, &t, &u);
	    
	      /* find the quadrature weight for current ip */
	      wt = Gq_surf_weight (ip, ielem_type);    
		  }
			  
	      /* ****************************************/
	      err = load_basis_functions( xi, bfd);
	      EH( err, "problem from load_basis_functions");
	    
	      err = beer_belly();
	      EH( err, "beer_belly");
	    
	      /* precalculate variables at  current integration pt.*/
	      err = load_fv();
	      EH( err, "load_fv");

	      err = load_fv_sens();
	      EH( err, "load_fv_sens");
	    
	      err = load_bf_grad();
	      EH( err, "load_bf_grad");

	      surface_determinant_and_normal (ei->ielem, iconnect_ptr, num_local_nodes, 
					      ielem_dim - 1,  
					      id_side,
					      num_nodes_on_side,
					      id_local_elem_coord );
	    
	    /*
	     * Load up physical space gradients of field variables at this
	     * Gauss point.
	     */
	    err = load_fv_grads();
	    EH( err, "load_fv_grads");

	    err = load_fv_grads_sens();
	    EH( err, "load_fv_grads_sens");

	    err = load_fv_mesh_derivs(1);
	    EH( err, "load_fv_mesh_derivs");

            if(TimeIntegration != STEADY && pd->e[MESH_DISPLACEMENT1])
      		{
	        for(j=0; j<ielem_dim; j++ )
       		   {
            		x_dot[j] = fv_dot->x[j];
          	   }
      		}
	    else
      		{
        	for(j=0; j<ielem_dim; j++ )
          	     {
            		x_dot[j] = 0.;
          	     }
      		}

	    do_LSA_mods(LSA_SURFACE);

	    if (ielem_dim !=3) {
	      calc_surf_tangent(ei->ielem, iconnect_ptr, 
				num_local_nodes, ielem_dim-1,
				num_nodes_on_side,
				id_local_elem_coord);
	    }

            weight = wt;
	    det = fv->sdet*fv->h3;

	    if(pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN)
	      EH(-1, "evalute_flux_sens has not been updated for the PROJECTED_CARTESIAN coordinate system.");

	    if(pd_glob[0] ->CoordinateSystem == CYLINDRICAL ||
               pd_glob[0] ->CoordinateSystem == SWIRLING   )
	      {
		weight = 2.*M_PIE*wt;
		det = fv->sdet;
	      }

/*  normal and tangent vector sensitivities   */
        memset( normal_sens, 0, sizeof(dbl)*MAX_PDIM);
        memset( tangent_sens, 0, sizeof(dbl)*2*MAX_PDIM);
        for( a=0 ; a < ielem_dim ; a++)
           {
           for( b=0 ; b < ielem_dim ; b++)
              {
               var = MESH_DISPLACEMENT1 + b ;
               if( pd->v[var] )
                  {
                   for(j=0 ; j < ei->dof[var]; j++)
                      {
                       normal_sens[a] += fv->dsnormal_dx[a][b][j]*bf[var]->phi[j]
                                         *(*esp_old->d[b][j]);
                       tangent_sens[0][a] += fv->dstangent_dx[0][a][b][j]*
                                      bf[var]->phi[j]*(*esp_old->d[b][j]);
                             if (ielem_dim !=3) {
                       tangent_sens[1][a] += fv->dstangent_dx[1][a][b][j]*
                                      bf[var]->phi[j]*(*esp_old->d[b][j]);
                                }
                      }
                  }
               }
            }

/*	get density  */

	    rho = density( NULL, time_value );


/**   get solid stresses  **/

 dim   = pd_glob[0]->Num_Dim;

 if(pd->e[R_MESH1] && cr->MeshMotion != ARBITRARY)
    {
      err = belly_flop(elc->lame_mu);
      EH(err, "error in belly flop");
      if (err == 2) exit(-1);
      /*
       * Total mesh stress tensor...
       *
       * Guess what the 4 new args will have to be and bloat accordingly.
       */

      err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, 
			       dTT_dp_gas, dTT_dporosity, dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
                               elc->lame_mu, elc->lame_lambda,
			       delta_t, ielem, ip, ip_total);

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
                if (cr->MeshFluxModel == NONLINEAR ||
                    cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
                    cr->MeshFluxModel == INCOMP_PSTRAIN )
         TT[2][2] = (1. - pow(fv->volume_change,2./3.)) * elast_modulus - fv->P;
                /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
                else  TT[2][2] = 0.;
            TT[1][2] = 0.;
            TT[0][2] = 0.;
          }
        }
    } /* end of STRESS_TENSOR */

/* calculate real-solid stress here !!*/
  if(pd->e[R_SOLID1] && cr->MeshMotion != ARBITRARY)
    {
      mu = elc_rs->lame_mu;
      err = belly_flop_rs(mu);
      EH(err, "error in belly flop");
      if (err == 2) return(err);
      /*
       * Total mesh stress tensor...
       */
       err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp, dTT_dc,  
				 dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
				 dTT_dT, dTT_dmax_strain, elc_rs->lame_mu, elc_rs->lame_lambda);
      if (dim == 2){
        elast_modulus = elc_rs->lame_mu;

        if (cr->RealSolidFluxModel == NONLINEAR ||
            cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
            cr->RealSolidFluxModel == INCOMP_PSTRAIN )
          TT[2][2] = (1. - pow(fv->volume_change,2./3.)) * elast_modulus - fv->P
;
        /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
        else  TT[2][2] = 0.;

        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    } /* end of REAL_STRESS_TENSOR */

		/*
		 *  viscoelastic stress tensor
		 *  assume only EVSS_F formulation for now
		 */
  			memset( ves_sens, 0, sizeof(dbl)*DIM*DIM);
  			memset( ves, 0, sizeof(dbl)*DIM*DIM);
  			if ( pd->v[POLYMER_STRESS11] )
    			  {
			    for ( a=0; a<VIM; a++)
        	              {
			        for ( b=0; b<VIM; b++)
            			   {
			             for ( ve_mode=0; ve_mode<vn->modes; ve_mode++)
                		       {
                  			ves_sens[a][b] += fv_sens->S[ve_mode][a][b];
                  			ves[a][b] += fv->S[ve_mode][a][b];
                			}
            			   }
        		      }
    			   }


		/*
		 * OK, let's simplify things by computing the viscous
		 * stress tensor upfront.
		 */
		      if(cr->MeshMotion == ARBITRARY)
    			 {
			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
			gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
		  	gamma_sens[a][b] = fv_sens->grad_v[a][b] + fv_sens->grad_v[b][a];
				}
			    }
			  
			  mu = viscosity(gn, gamma, d_mu);

		/*  compute implicit viscosity sensitivities  */
			mu_sens = 0.;
  			v = TEMPERATURE;
  			if ( pd->v[v] )
    			   {
      			     dofs  = ei->dof[v];
      			     for ( j=0; j<dofs; j++)
				{
	  			 mu_sens += *esp_old->T[j] * d_mu->T[j];
				}
    			    }
		
		/*  velocity  */

  			velodim = dim;	
  			if(pd->CoordinateSystem == SWIRLING ||
			     pd->CoordinateSystem == PROJECTED_CARTESIAN)
    			velodim = dim + 1;

			for ( a=0; a<velodim; a++)
    			  {
      			   v = VELOCITY1 + a;
      			     if ( pd->v[v] )
				{
	  			 dofs     = ei->dof[v];
	  			   for ( j=0; j<dofs; j++)
	    				{
	      			mu_sens += *esp_old->v[a][j] * d_mu->v[a][j];
	    				}
				}
    			  }

		/* mesh */


  			for ( a=0; a<dim; a++)
    			   {
      			      v = MESH_DISPLACEMENT1 + a;
      				if ( pd->v[v] )
				   {
	  			     dofs     = ei->dof[v];
	  			     for ( j=0; j<dofs; j++)
	    				{
	      			mu_sens += *esp_old->d[a][j] * d_mu->X[a][j];
	    				}
				   }
    			   }

   		/* Pressure...  */

			v = PRESSURE;
  			if ( pd->v[v] )
    			  {
			      dofs  = ei->dof[v];
			      for ( j=0; j<dofs; j++)
				{
				  mu_sens += *esp_old->P[j] * d_mu->P[j];
				}
    			  }

		/*  Concentration  */

			v = MASS_FRACTION;
  			if (pd->v[v]) {
			    for (a = 0; a < pd->Num_Species_Eqn; a++) {
			      dofs = ei->dof[v];
			      for (j = 0; j < dofs; j++) {
				mu_sens += *esp_old->c[a][j] * d_mu->C[a][j];
      				}
    			     }
			  }
  		/* * Fill...  */

			  v = FILL;
			  if ( pd->v[v] )
			    {
			      dofs  = ei->dof[v];
			      for ( j=0; j<dofs; j++)
				{
				  mu_sens += *esp_old->F[j] * d_mu->F[j];
				}
			    }
/***   
 *
 *      end of material property sensitivity pre-calcs
 */

			  for ( a=0; a<VIM; a++)
			    {
			      for ( b=0; b<VIM; b++)
				{
				  vs_sens[a][b] = 
					mu*gamma_sens[a][b]-fv_sens->P*delta(a,b)
				+ gamma[a][b]*mu_sens;
				  vs[a][b] = mu*gamma[a][b]-fv->P*delta(a,b);
				}
			    }
    			 }

 		/*
   		 * computing the electric stress tensor upfront.
   		 */
  			memset( es, 0, sizeof(dbl)*DIM*DIM);
  			memset( es_sens, 0, sizeof(dbl)*DIM*DIM);
		      if(pd->e[R_POTENTIAL])
     			 {
 			  efield_sqr = 0.0;
 			  efield_sqr_sens = 0.0;
 			  for ( a=0; a<VIM; a++)
 			    {
 				  efield[a] = - fv->grad_V[a];
 				  efield_sens[a] = - fv_sens->grad_V[a];
 				  efield_sqr += efield[a]*efield[a];
 				  efield_sqr_sens += 2.*efield[a]*efield_sens[a];
 			    }
 
 			  for ( a=0; a<VIM; a++)
 			    {
 			      for ( b=0; b<VIM; b++)
 				{
 				  es[a][b] = efield[a]*efield[b] 
                                            - 0.5*efield_sqr*delta(a,b);
 				  es_sens[a][b] = efield[a]*efield_sens[b] 
                                                + efield_sens[a]*efield[b] 
                                                - 0.5*efield_sqr_sens*delta(a,b);
 				}
 			    }
     			 }
 		/*
   		 * load surface tension for variable models
   		 */
 		      if (mp->SurfaceTensionModel == CONSTANT) {
                            dsigmadT = 0.0;
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = 0.0;}
                            }
                      sigma_sens = 0.0;
 		      if (mp->SurfaceTensionModel != CONSTANT) {
 		            load_surface_tension(dsigma_dx);
 		            }
 		      if (mp->SurfaceTensionModel == USER) {
                            dsigmadT = mp->d_surface_tension[TEMPERATURE];
                            sigma_sens = dsigmadT*fv_sens->T;
 			    for ( a=0; a<MAX_CONC; a++)
                                {
                    dsigmadC[a] = mp->d_surface_tension[MAX_VARIABLE_TYPES+a];
                    sigma_sens += dsigmadC[a]*fv_sens->c[a];
                                }
                            }


	    switch( quantity )
	      {

	      case HEAT_FLUX:
	      
		/* finally we can add up the outgoing flux and area */
		/* but first evaluate properties if variable */
		if (mp->ConductivityModel == USER)
		  {
		    err = usr_thermal_conductivity(mp->u_thermal_conductivity, time_value);
		  }
		/*  heat capacity  */
		if(mp->HeatCapacityModel == USER )
		{
		      err = usr_heat_capacity(mp->u_heat_capacity, time_value);
		      Cp = mp->heat_capacity;
		}
		else if (mp->HeatCapacityModel == CONSTANT )
		{
		Cp    = mp->heat_capacity;
		}
		else if (mp->HeatCapacityModel == ENTHALPY )
		{
		  Cp = enthalpy_heat_capacity_model(d_Cp);
		}
		  else
		{
		EH(-1,"Unrecognized heat capacity model");
		}

                if ( cr->HeatFluxModel == CR_HF_FOURIER_0 )
                {
                for (j=0; j<VIM; j++)
                  {
                    local_q += (mp->thermal_conductivity *
                                    fv->snormal[j] * fv_sens->grad_T[j] );
                    local_qconv += ( rho * Cp*fv->T*
                                    fv->snormal[j] * (fv_sens->v[j]-x_dot[j]));
                  }
                }  else if ( cr->HeatFluxModel == CR_HF_USER )
                {

#if defined SECOR_HEAT_FLUX
        double *hpar, h, dh_dX[DIM], Vb[DIM],Vt[DIM];
        double dq_dVb[DIM][DIM], dq_dVt[DIM][DIM];

        hpar = &mp->u_thermal_conductivity[0];
        h = hpar[0] + hpar[4]*fv->x[0]
                + (hpar[1]-hpar[5]*fv->x[0])*(hpar[3]-fv->x[1])
                        + 0.5*hpar[2]*SQUARE(hpar[3]-fv->x[1]);

        dh_dX[0] = hpar[4] - hpar[5]*(hpar[3]-fv->x[1]);
        dh_dX[1] = hpar[5]*fv->x[0]-hpar[1] - hpar[2]*(hpar[3]-fv->x[1]);

/*     velocities of bottom and top surfaces   */
        Vb[0] = mp->u_heat_capacity[0];
        Vb[1] = mp->u_heat_capacity[1];
        Vt[0] = mp->u_heat_capacity[2];
        Vt[1] = mp->u_heat_capacity[3];

      usr_heat_flux(fv->grad_T, q, dq_gradT, dq_dX, time_value, h, dh_dX, Vb,Vt,
                             dq_dVb, dq_dVt);
#else
      usr_heat_flux(fv->grad_T, q, dq_gradT, dq_dX, time_value);
      printf("untested\n");
      exit(-1);
#endif
                for (j=0; j<VIM; j++)
                  {
                  for (a=0; a<VIM; a++)
                    {
            local_q += fv->snormal[j] * dq_gradT[j][a] * fv_sens->grad_T[a];
                    }
                  }
                }
                    local_flux += weight * fv->sdet * fv->h3 * local_q;
                    local_flux_conv += weight * fv->sdet * fv->h3 * local_qconv;
                break;

	      case VOLUME_FLUX:
                for(j=0; j<VIM; j++)
                  {
			  if (cr->MeshMotion == ARBITRARY)
                    local_q += ( fv->snormal[j]*(fv_sens->v[j]-x_dot[j])
                                 + normal_sens[j]*(fv->v[j]-x_dot[j]));
			  else if (pd->v[MESH_DISPLACEMENT1])
                    local_q += fv->snormal[j]*fv_sens->d[j]
                               + normal_sens[j]*fv->d[j];
			  else
			    EH(-1,"Inconsistency in volume-flux specification. Contact Developers");

                  }
                    local_flux += weight*det*local_q;
		break;

	      case PVELOCITY1:
	      case PVELOCITY2:
	      case PVELOCITY3:
                for(j=0; j<VIM; j++)
                  {
                    local_q += fv->snormal[j]*fv_sens->pv[j]
                               + normal_sens[j]*fv->pv[j];
                  }
		local_flux += weight*det*local_q;
		break;

	      case EXT_VELOCITY:
                local_q += fv_sens->ext_v;
                 
		local_flux += weight*det*local_q;
		break;

	      case EFIELD1:
	      case EFIELD2:
	      case EFIELD3:
                for(j=0; j<VIM; j++)
                  {
                    local_q += fv->snormal[j]*fv_sens->E_field[j]
                               +normal_sens[j]*fv->E_field[j];
                  }
		local_flux += weight*det*local_q;
		break;

	      case SPECIES_FLUX:
                for(j=0; j<VIM; j++)
                  {
                    local_q += ( mp->diffusivity[species_id]*(
                          fv->snormal[j]*fv_sens->grad_c[species_id][j]
                            + normal_sens[j]*fv->grad_c[species_id][j])  );

                    local_qconv += ( fv->snormal[j]*
				     ((fv_sens->v[j]-x_dot[j])*fv->c[species_id]
				      + (fv->v[j]-x_dot[j])*fv_sens->c[species_id] )
                             +normal_sens[j]*(fv->v[j]-x_dot[j])*fv->c[species_id]);
                  }
                    local_flux += weight*det*local_q;
                    local_flux_conv += weight*det*local_qconv;
                break;

	      case TORQUE:
		if(pd->CoordinateSystem == PROJECTED_CARTESIAN)
		  EH(-1, "TORQUE has not been updated for the PROJECTED_CARTESIAN coordinate system.");

		if(pd->CoordinateSystem == SWIRLING || 
		   pd->CoordinateSystem == CYLINDRICAL)
		  {
		    for ( a=0; a<VIM; a++)
		      {
			for ( b=0; b<VIM; b++)
			  {
			/*
			 *  note that for CYLINDRICAL and SWIRLING coordinate systems
			 * the h3 factor has been incorporated already into sdet
			 * the moment arm is incorporated into sideset.
			 */
                        local_q += e_theta[a]*((fv->x[1] * 
                        	(vs_sens[a][b]+ves_sens[a][b]) +
				fv_sens->x[1]*(vs[a][b]+ves[a][b]))
				*fv->snormal[b] +  
                                  normal_sens[b]*(fv->x[1]*(vs[a][b]+ves[a][b])));
			  }
		      }
                        local_flux += weight * det* local_q;
		  }
		else if( pd->CoordinateSystem == CARTESIAN ) 
		  {
                    if( VIM == 2 ) 
                       {
                         fv->x[2] = 0.0;
                         fv->snormal[2] = 0.0;
                         normal_sens[2] = 0.0;
                         vs[0][2] = 0.0;
                         vs[1][2] = 0.0;
                         vs[2][0] = 0.0;
                         vs[2][1] = 0.0;
                         vs[2][2] = -fv->P;
                       }
		    for ( a=0; a<DIM; a++)
		      {
                        Tract[a] = 0.0;
                        Tract_sens[a] = 0.0;
			for ( b=0; b<DIM; b++)
			  {
		            Tract[a] += (vs[a][b] + ves[a][b]) * fv->snormal[b];
			    Tract_sens[a] += (vs_sens[a][b] + ves_sens[a][b])
                                      * fv->snormal[b] + normal_sens[b]*
                                        (vs[a][b]+ves[a][b]);
			  }
		      }
		
		    for ( a=0; a<DIM; a++)
		      {
			for ( b=0; b<DIM; b++)
			  {
			    for ( c=0; c<DIM; c++)
			      {
                               local_Torque[a] += (permute(b,c,a) * 
                                             (  fv->x[b]*Tract_sens[c] 
                                               + fv_sens->x[b]*Tract[c] ) );
			      }
			  }
		      }
		    for ( a=0; a<DIM; a++)
		      { Torque[a] +=  weight * det * local_Torque[a]; }
                   local_flux = Torque[2];
		  }
		else
		  {
		    EH(-1,"Torque cannot be calculated in this case.");
		  }

		break;

      case FORCE_NORMAL:
		if(cr->MeshMotion == ARBITRARY)
		{
                for ( a=0; a<VIM; a++)
                  {
                    for ( b=0; b<VIM; b++)
                      {
                local_q += (fv->snormal[a]*
			(vs_sens[a][b] +ves_sens[a][b])
                                *fv->snormal[b]
                           +(vs[a][b]+ves[a][b])*(fv->snormal[a]*normal_sens[b]
                               +normal_sens[a]*fv->snormal[b]));
                local_qconv += rho*(fv->snormal[a]*((fv->v[a]-x_dot[a])*fv_sens->v[b]
                 + (fv_sens->v[a]-x_dot[a])*fv->v[b])*fv->snormal[b]
                    + (fv->v[a]-x_dot[a])*fv->v[b]*(fv->snormal[a]*normal_sens[b]
                         +normal_sens[a]*fv->snormal[b]));
                      }
                  }
                local_flux += weight * det * local_q;
                local_flux_conv += weight * det * local_qconv;
                }
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
                for ( a=0; a<VIM; a++)
                  {
                    for ( b=0; b<VIM; b++)
                      {
                        local_q += ( fv->snormal[a]*TT_sens[a][b]*fv->snormal[b]
                                + TT[a][b]*(fv->snormal[a]*normal_sens[b]+
                                       normal_sens[a]*fv->snormal[b]));
                      }
                  }
                        local_flux += weight * det * local_q;
                }
		break;

	      case FORCE_TANGENT1:

		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
		    for ( b=0; b<VIM; b++)
		      {
                        local_q += (fv->stangent[0][a]*(vs_sens[a][b]
				+ves_sens[a][b])*fv->snormal[b]);
                        local_qconv += ( rho*fv->stangent[0][a]*((fv->v[a]-x_dot[a])
        *fv_sens->v[b]+(fv_sens->v[a]-x_dot[a])*fv->v[b])*fv->snormal[b]);

		      }
		  }
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det* local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available"); 
		for ( a=0; a<VIM; a++)
		  {
		    for ( b=0; b<VIM; b++)
		      {
                local_q += ( fv->stangent[0][a]*TT_sens[a][b]*fv->snormal[b]);
		      }
		  }
                local_flux += weight * det* local_q;
		}
		break;

	      case FORCE_TANGENT2:

		if(cr->MeshMotion == ARBITRARY)
		{
		  if(pd->Num_Dim == 3)
		    {
                for ( a=0; a<VIM; a++)
                  {
                    for ( b=0; b<VIM; b++)
                      {
                        local_q += ( fv->stangent[1][a]*(vs_sens[a][b]
					+ves_sens[a][b])*fv->snormal[b]);
                        local_qconv += ( rho*fv->stangent[1][a]*((fv->v[a]-x_dot[a])
        *fv_sens->v[b]+(fv_sens->v[a]-x_dot[a])*fv->v[b])*fv->snormal[b]);
                      }
                  }
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det* local_qconv;
                }
			else
			{
			EH(-1, "Illegal flux type");
			}
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available"); 
                for ( a=0; a<VIM; a++)
                  {
                    for ( b=0; b<VIM; b++)
                      {
                local_q += ( fv->stangent[1][a]*TT_sens[a][b]*fv->snormal[b]);
                      }
                  }
                local_flux += weight * det* local_q;
                }
		break;

	      case FORCE_X:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[0][a] + ves_sens[0][a])
		   	    *fv->snormal[a] +normal_sens[a]*(vs[0][a]+ves[0][a]));

                        local_qconv += rho*(((fv->v[0]-x_dot[0])*fv_sens->v[a]
                            +(fv_sens->v[0]-x_dot[0])*fv->v[a])*fv->snormal[a]
                            +normal_sens[a]*(fv->v[0]-x_dot[0])*fv->v[a]);
		  }
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += (TT_sens[0][a]*fv->snormal[a]
                                    +TT[0][a]*normal_sens[a]);
		  }
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_X_POS:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[0][a] + ves_sens[0][a])
				*fv->snormal[a]) ;

                        local_qconv += (rho*((fv->v[0]-x_dot[0])
             *fv_sens->v[a]+(fv_sens->v[0]-x_dot[0])*fv->v[a])*fv->snormal[a]);
		  }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += (TT_sens[0][a]*fv->snormal[a]);
		  }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_X_NEG:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[0][a] + ves_sens[0][a])
				*fv->snormal[a]) ;

                        local_qconv += (rho*((fv->v[0]-x_dot[0])
             *fv_sens->v[a]+(fv_sens->v[0]-x_dot[0])*fv->v[a])*fv->snormal[a]);
		  }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += (TT_sens[0][a]*fv->snormal[a]);
		  }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Y:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[1][a] + ves_sens[1][a])
				*fv->snormal[a]);

                        local_qconv += (rho*((fv->v[1]-x_dot[1])
        *fv_sens->v[a]+(fv_sens->v[1]-x_dot[1])*fv->v[a])*fv->snormal[a]);
		  }
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ( TT_sens[1][a]*fv->snormal[a]);
		  }
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Y_POS:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[1][a] + ves_sens[1][a])
				*fv->snormal[a]);

                        local_qconv += (rho*((fv->v[1]-x_dot[1])
        *fv_sens->v[a]+(fv_sens->v[1]-x_dot[1])*fv->v[a])*fv->snormal[a]);
		  }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ( TT_sens[1][a]*fv->snormal[a]);
		  }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Y_NEG:
		if(cr->MeshMotion == ARBITRARY)
		{
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ((vs_sens[1][a] + ves_sens[1][a])
				*fv->snormal[a]);

                        local_qconv += (rho*((fv->v[1]-x_dot[1])
        *fv_sens->v[a]+(fv_sens->v[1]-x_dot[1])*fv->v[a])*fv->snormal[a]);
		  }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight * det * local_qconv;
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available");
		for ( a=0; a<VIM; a++)
		  {
                        local_q += ( TT_sens[1][a]*fv->snormal[a]);
		  }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Z:
		if(cr->MeshMotion == ARBITRARY)
		{
			if(pd->Num_Dim == 3)
			{
			for ( a=0; a<VIM; a++)
		  	   {
                        	local_q += ((vs_sens[2][a] + ves_sens[2][a])
					*fv->snormal[a]);

                       		local_qconv += (rho*((fv->v[2]-x_dot[2])
        *fv_sens->v[a]+(fv_sens->v[2]-x_dot[2])*fv->v[a])*fv->snormal[a]);
		  	   }
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight *det * local_qconv;
			}
			else
			{
			EH(-1, "Illegal flux type");
			}
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available"); 
		for ( a=0; a<VIM; a++)
		  {
                       local_q += (TT_sens[2][a]*fv->snormal[a]);
		  }
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Z_POS:
		if(cr->MeshMotion == ARBITRARY)
		{
			if(pd->Num_Dim == 3)
			{
			for ( a=0; a<VIM; a++)
		  	   {
                        	local_q += ((vs_sens[2][a] + ves_sens[2][a])
					*fv->snormal[a]);

                       		local_qconv += (rho*((fv->v[2]-x_dot[2])
        *fv_sens->v[a]+(fv_sens->v[2]-x_dot[2])*fv->v[a])*fv->snormal[a]);
		  	   }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight *det * local_qconv;
			}
			else
			{
			EH(-1, "Illegal flux type");
			}
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available"); 
		for ( a=0; a<VIM; a++)
		  {
                       local_q += (TT_sens[2][a]*fv->snormal[a]);
		  }
        	        if(local_q <0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      case FORCE_Z_NEG:
		if(cr->MeshMotion == ARBITRARY)
		{
			if(pd->Num_Dim == 3)
			{
			for ( a=0; a<VIM; a++)
		  	   {
                        	local_q += ((vs_sens[2][a] + ves_sens[2][a])
					*fv->snormal[a]);

                       		local_qconv += (rho*((fv->v[2]-x_dot[2])
        *fv_sens->v[a]+(fv_sens->v[2]-x_dot[2])*fv->v[a])*fv->snormal[a]);
		  	   }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
                        local_flux_conv += weight *det * local_qconv;
			}
			else
			{
			EH(-1, "Illegal flux type");
			}
		}
		else
		{
	EH(-1,"Force Sensitivity Calculation for Solids  not Available"); 
		for ( a=0; a<VIM; a++)
		  {
                       local_q += (TT_sens[2][a]*fv->snormal[a]);
		  }
        	        if(local_q >0) {local_q=0;}
                        local_flux += weight * det * local_q;
		}
		break;

	      default:

		EH(-1, "Illegal flux type");
		break;
	      }  /*  end of switch */
	    
	    local_area += weight * det ;
#ifdef PARALLEL
        delta_flux = local_flux - local_flux0;
        delta_flux_conv = local_flux_conv - local_flux_conv0;
        delta_area = local_area - local_area0;

        if( Num_Proc > 1 &&  dpi->elem_owner[ elem_list[i] ] == ProcID )
         {
          proc_flux += delta_flux;
          proc_flux_conv += delta_flux_conv;
          proc_area += delta_area;
         }

        local_flux0 = local_flux;
        local_flux_conv0 = local_flux_conv;
        local_area0 = local_area;
#endif

        if (profile_flag && print_flag ) {
        FILE  *jfp = 0;
        if(  (filenm != NULL) && ( (jfp=fopen(filenm,"a")) != NULL))
            {
             fprintf(jfp," %g  %g  %g  %g  %g\n",
                        fv->x[0],fv->x[1],fv->x[2],local_q,local_qconv);
              fflush(jfp);
            }
          fclose(jfp);
                }

	    }   /*  surface integration over element */
	      if( Subgrid_Int.active )
		{
		  safe_free ( (void *) Subgrid_Int.s  ); Subgrid_Int.s = NULL;
		  safe_free ( (void *) Subgrid_Int.wt ); Subgrid_Int.wt = NULL;
		  free_search_grid ( &element_search_grid );
		}

	  }   /*  material id conditional */
	}    /*   element loop */
    }      /*   sset id   */
  else
    {
/**  Apply end point conditions when the nset is not found   **/
  	nset_id = in_list(side_set_id, 0, exo->num_node_sets, exo->ns_id);
  	if ( nset_id != -1 )
    		{
		int corner_elem=-1,local_node, gnn, dir, kine_sset;
		double sign;
      		num_node_in_set      = exo->ns_num_nodes[nset_id];
      		num_dist_fact_in_set = exo->ns_num_distfacts[nset_id];
      		elem_list = &exo->ns_node_list[exo->ns_node_index[nset_id]];
		if (num_node_in_set != 1) 
		   EH(-1,"more than one node, this is for 2D only");
		gnn = elem_list[0];
		sign = (species_id < 0 ) ? -1. : 1.;
		kine_sset = abs(species_id);
	/*  Find the free surface element in the adjoining sset	*/
  		sset_id = in_list(kine_sset, 0, exo->num_side_sets, exo->ss_id);
		if( sset_id != -1)
		  {
      		  num_nodes_on_side    = exo->ss_num_distfacts[sset_id]/exo->ss_num_sides[sset_id]; /* Well... */
		  for(i=0 ; i<exo->ss_num_sides[sset_id] ; i++)
		     {
		      ielem = exo->ss_elem_list[exo->ss_elem_index[sset_id]+i];
		      if( in_list(gnn, exo->elem_node_pntr[ielem],
				exo->elem_node_pntr[ielem+1],
				exo->elem_node_list) != -1)
				{
				corner_elem = ielem;
     		local_node = in_list(gnn, exo->elem_node_pntr[corner_elem],
                                exo->elem_node_pntr[corner_elem+1],
                                exo->elem_node_list);
		local_node -= exo->elem_node_pntr[corner_elem];

	  	   mn = find_mat_number(corner_elem, exo);

	        if ( mat_id == -1 || 
                     (mn == map_mat_index(mat_id) && 
                      dpi->elem_owner[corner_elem] == ProcID)) {

	  	   err = load_elem_dofptr( corner_elem, 
				  (Exo_DB*) exo,
				  (dbl *) x,
				  (dbl *) &x_sens_p[vector_id][0],
				  (dbl *) xdot,
				  (dbl *) xdot,
				  (dbl *) x,
				  0);
	  	   EH(err, "load_elem_dofptr");

	  	   err = bf_mp_init(pd);
	  	   EH(err, "bf_mp_init");
  
	           iconnect_ptr        = ei->iconnect_ptr;
	           ielem_type          = ei->ielem_type;
	      	   ip_total            = elem_info(NQUAD_SURF, ielem_type);
	      	   num_local_nodes     = ei->num_local_nodes;
	           ielem_dim           = ei->ielem_dim;
	      
	 
	           id_side = find_id_side (ei->ielem, num_nodes_on_side,
				      &exo->ss_node_list[sset_id]
				      [num_nodes_on_side*i],
				      id_local_elem_coord, exo);

            	   find_nodal_stu(local_node, ei->ielem_type, xi, xi+1, xi+2);

            	   err = load_basis_functions( xi, bfd);
            	   EH( err, "problem from load_basis_functions");

            	   err = beer_belly();
            	   EH( err, "beer_belly");

            	   err = load_fv();
            	   EH( err, "load_fv");

                   err = load_fv_sens();
                   EH( err, "load_fv_sens");

            	   err = load_bf_grad();
            	   EH( err, "load_bf_grad");

            	   err = load_bf_mesh_derivs();
            	   EH(err, "load_bf_mesh_derivs");

            	   surface_determinant_and_normal(corner_elem,
                                           exo->elem_node_pntr[corner_elem],
                                           ei->num_local_nodes,
                                           ei->ielem_dim-1,
                                           id_side,
                                           num_nodes_on_side,
                                           id_local_elem_coord);
	 	   err = load_fv_grads();
	 	   EH( err, "load_fv_grads");
			  
                   err = load_fv_grads_sens();
                   EH( err, "load_fv_grads_sens");

	 	   err = load_fv_mesh_derivs(1);
	 	   EH( err, "load_fv_mesh_derivs");

	    	   if (ielem_dim !=3)
	      		{
			 calc_surf_tangent (corner_elem, iconnect_ptr, 
				   num_local_nodes, ielem_dim-1,
				   num_nodes_on_side,
				   id_local_elem_coord);
	      		}
		dim = pd_glob[0]->Num_Dim;

/*  normal and tangent vector sensitivities   */
        memset( normal_sens, 0, sizeof(dbl)*MAX_PDIM);
        memset( tangent_sens, 0, sizeof(dbl)*2*MAX_PDIM);
        for( a=0 ; a < ielem_dim ; a++)
           {
           for( b=0 ; b < ielem_dim ; b++)
              {
               var = MESH_DISPLACEMENT1 + b ;
               if( pd->v[var] )
                  {
                   for(j=0 ; j < ei->dof[var]; j++)
                      {
                       normal_sens[a] += fv->dsnormal_dx[a][b][j]*bf[var]->phi[j]
                                         *(*esp_old->d[b][j]);
                       tangent_sens[0][a] += fv->dstangent_dx[0][a][b][j]*
                                      bf[var]->phi[j]*(*esp_old->d[b][j]);
                             if (ielem_dim !=3) {
                       tangent_sens[1][a] += fv->dstangent_dx[1][a][b][j]*
                                      bf[var]->phi[j]*(*esp_old->d[b][j]);
                                }
                      }
                  }
               }
            }
 		/*
   		 * load surface tension for variable models
    		 */
 		      if (mp->SurfaceTensionModel == CONSTANT) {
                            dsigmadT = 0.0;
 			    for ( a=0; a<MAX_CONC; a++)
                                {dsigmadC[a] = 0.0;}
                            }
                      sigma_sens = 0.0;
 		      if (mp->SurfaceTensionModel != CONSTANT) {
 		            load_surface_tension(dsigma_dx);
 		            }
 		      if (mp->SurfaceTensionModel == USER) {
                            dsigmadT = mp->d_surface_tension[TEMPERATURE];
                            sigma_sens = dsigmadT*fv_sens->T;
 			    for ( a=0; a<MAX_CONC; a++)
                                {
                    dsigmadC[a] = mp->d_surface_tension[MAX_VARIABLE_TYPES+a];
                    sigma_sens += dsigmadC[a]*fv_sens->c[a];
                                }
                            }
		  switch( quantity )
		    {
		    case FORCE_X:
		    case FORCE_Y:
		    case FORCE_Z:
		    case FORCE_X_POS:
		    case FORCE_Y_POS:
		    case FORCE_Z_POS:
		    case FORCE_X_NEG:
		    case FORCE_Y_NEG:
		    case FORCE_Z_NEG:
			  dir = ( quantity == FORCE_X ? 0 :
				( quantity == FORCE_Y ? 1 : 
			        ( quantity == FORCE_Z ? 2 : 
			        ( quantity == FORCE_X_POS ? 0 :
				( quantity == FORCE_Y_POS ? 1 : 
			        ( quantity == FORCE_Z_POS ? 2 : 
			        ( quantity == FORCE_X_NEG ? 0 :
				( quantity == FORCE_Y_NEG ? 1 : 
			        ( quantity == FORCE_Z_NEG ? 2 : -1))))))))) ;  

		      if(cr->MeshMotion == ARBITRARY)
			{
		       if( pd->CoordinateSystem == CARTESIAN ) 
			  {
                         local_flux += sign*(mp->surface_tension*tangent_sens[0][dir] + sigma_sens*fv->stangent[0][dir]);
			}
		      else if(pd->CoordinateSystem == SWIRLING || 
			 pd->CoordinateSystem == CYLINDRICAL)
			{
                         local_flux += 2*M_PIE*sign*
                            (fv_sens->x[1]*mp->surface_tension*fv->stangent[0][dir] 
                            + fv->x[1]*sigma_sens*fv->stangent[0][dir] + 
                            fv->x[1]*mp->surface_tension*tangent_sens[0][dir]);
			}
			else
			    {
			EH(-1, "Point force has not been updated for the PROJECTED_CARTESIAN coordinate system.");
			    }
			}
		      break;
		    default:
		      EH(-1, "Illegal flux type");
		      break;
		    }  /*  end of switch */
#ifdef PARALLEL
        delta_flux = local_flux - local_flux0;
        delta_flux_conv = local_flux_conv - local_flux_conv0;
        delta_area = local_area - local_area0;

        if( Num_Proc > 1 &&  dpi->elem_owner[ elem_list[i] ] == ProcID )
         {
          proc_flux += delta_flux;
          proc_flux_conv += delta_flux_conv;
          proc_area += delta_area;
         }

        local_flux0 = local_flux;
        local_flux_conv0 = local_flux_conv;
        local_area0 = local_area;
#endif
                 }     /*  mat_id  */
              }        
           }           /*  ss_sides loop  */
	if(corner_elem == -1) EH(-1,"corner element not found");
         }             /*  if sset_id     */
      }                /*  if nset_id     */
  		else
    		{
#ifndef PARALLEL
      (void) sprintf(Err_Msg, "%s could not locate SSID %d.", yo, side_set_id);
      EH(-1, Err_Msg);
#endif      
    		}
    }
         

#ifdef PARALLEL
  MPI_Allreduce(&local_flux, &global_flux, 1, MPI_DOUBLE, MPI_SUM, 
		MPI_COMM_WORLD);
  MPI_Allreduce(&local_flux_conv, &global_flux_conv, 1, MPI_DOUBLE, MPI_SUM, 
		MPI_COMM_WORLD);

  MPI_Allreduce(&local_area, &global_area, 1, MPI_DOUBLE, MPI_SUM, 
		MPI_COMM_WORLD);

  local_flux = global_flux;
  local_flux_conv = global_flux_conv;
  local_area = global_area;
#endif  

/**   skip this section if not printing - i.e. for AC conditions  **/

if(print_flag && ProcID == 0)
{
	FILE  *jfp;
	if(  (filenm != NULL) && ( (jfp=fopen(filenm,"a")) != NULL)) 
	  { 
            if(quantity==TORQUE)
              {
                fprintf(jfp," torque_sens= %e  \n", local_flux);
                fprintf(jfp, "\n");
                fflush(jfp);
              }
            else
              {
                fprintf(jfp," flux_sens= %e %e  area= %e  \n", local_flux
					,local_flux_conv,local_area);
                fprintf(jfp, "\n");
                fflush(jfp);
              }
	    fclose(jfp);
	  }
}  /*  end of print_flag conditional */
  af->Assemble_Jacobian = Jac_state;
  return(local_flux + local_flux_conv);		/* failsafe default? */
}
/***  end of evaluate_flux_sens  **/



/* load_fv_sens() -- load up values of solution sensitivities that have 
 *	pointers loaded into the esp->R_structure
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
 *			P -- pressure		(scalar)
 *			S -- polymer stress	(tensor)
 *			G -- velocity gradient	(tensor)
 *                     pv -- particle velocity  (vector)
 *                     pG -- particle velocity gradient (tensor)
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 * Created:	Fri Mar 18 06:44:41 MST 1994 pasacki@sandia.gov
 *
 * Modified:	
 */

static int 
load_fv_sens(void)
{
  int v;			/* variable type indicator */
  int i;			/* index */
  int p, q;			/* dimension indeces */
  int dim;
  int velodim;			/* Someday...we might have more velocity */
				/* components than we have spatial dimensions */
				/* these are the 2.5 dimensional problems.*/
  int dofs;			/* degrees of freedom for a var in the elem */
  int w;			/* concentration species counter */
  int mode;

  int status;

  int v_s[MAX_MODES][DIM][DIM];

  const int v_g[DIM][DIM] = 
  { 
    { VELOCITY_GRADIENT11, VELOCITY_GRADIENT12, VELOCITY_GRADIENT13 }, 
    { VELOCITY_GRADIENT21, VELOCITY_GRADIENT22, VELOCITY_GRADIENT23 }, 
    { VELOCITY_GRADIENT31, VELOCITY_GRADIENT32, VELOCITY_GRADIENT33 }
  };

  status = 0;

  /* load eqn and variable number in tensor form */
  stress_eqn_pointer(v_s);

  dim = ei->ielem_dim;

  /*
   * Temperature...
   */

  v = TEMPERATURE;
  fv_sens->T = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->T += *esp_old->T[i] * bf[v]->phi[i];
	}
    }


  /*
   * Fill...
   */

  v = FILL;
  fv_sens->F = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->F += *esp_old->F[i] * bf[v]->phi[i];
	}
    }

  /*
   * Voltage...
   */

  v = VOLTAGE;
  fv_sens->V = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->V += *esp_old->V[i] * bf[v]->phi[i];
	}
    }

  /*
   *  Surface charge density ...
   */

  v = SURF_CHARGE;
  fv_sens->qs = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->qs += *esp_old->qs[i] * bf[v]->phi[i];
        }
    }

  /*
   *  Structural shell curvature
   */

  v = SHELL_CURVATURE;
  fv_sens->sh_K = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_K += *esp_old->sh_K[i] * bf[v]->phi[i];
        }
    }

  v = SHELL_CURVATURE2;
  fv_sens->sh_K2 = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_K2 += *esp_old->sh_K2[i] * bf[v]->phi[i];
        }
    }

  /*
   *  Structural shell tension
   */

  v = SHELL_TENSION;
  fv_sens->sh_tens = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_tens += *esp_old->sh_tens[i] * bf[v]->phi[i];
        }
    }

  /*
   *  Structural shell x coordinate
   */

  v = SHELL_X;
  fv_sens->sh_x = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_x += *esp_old->sh_x[i] * bf[v]->phi[i];
        }
    }

  /*
   *  Structural shell x coordinate
   */

  v = SHELL_Y;
  fv_sens->sh_y = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_y += *esp_old->sh_y[i] * bf[v]->phi[i];
        }
    }
 
  /*
   *  shell user
   */
 
  v = SHELL_USER;
  fv_sens->sh_u = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_u += *esp_old->sh_u[i] * bf[v]->phi[i];
        }
    }

  /*
   *  Shell orientation angles
   */

  velodim = dim;	
  for ( p=0; p<velodim-1; p++)
    {
 
      v = SHELL_ANGLE1 + p;
      fv_sens->sh_ang[p]     = 0.;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->sh_ang[p] += *esp_old->sh_ang[p][i] * bf[v]->phi[i];
	    }
	}
    }
    

  /*
   * Acoustic Pressure
   */

  v = ACOUS_PREAL;
  fv_sens->apr = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->apr += *esp_old->apr[i] * bf[v]->phi[i];
	}
    }
  v = ACOUS_PIMAG;
  fv_sens->api = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->api += *esp_old->api[i] * bf[v]->phi[i];
	}
    }

  v = ACOUS_REYN_STRESS;
  fv_sens->ars = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->ars += *esp_old->ars[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_BDYVELO;
  fv_sens->sh_bv = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_bv += *esp_old->sh_bv[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_LUBP;
  fv_sens->sh_p = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->sh_p += *esp_old->sh_p[i] * bf[v]->phi[i];
      }
    }

  v = LUBP;
  fv_sens->lubp = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->lubp += *esp_old->lubp[i] * bf[v]->phi[i];
      }
    }

  v = LUBP_2;
  fv_sens->lubp_2 = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->lubp_2 += *esp_old->lubp_2[i] * bf[v]->phi[i];
      }
    }

  v = SHELL_FILMP;
  fv_sens->sh_fp = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->sh_fp += *esp_old->sh_fp[i] * bf[v]->phi[i];
      }
    }

  v = SHELL_FILMH;
  fv_sens->sh_fh = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->sh_fh += *esp_old->sh_fh[i] * bf[v]->phi[i];
      }
    }

  v = SHELL_PARTC;
  fv_sens->sh_pc = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
      {
        fv_sens->sh_pc += *esp_old->sh_pc[i] * bf[v]->phi[i];
      }
    }


  v = SHELL_SAT_CLOSED;
  fv_sens->sh_sat_closed = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_sat_closed += *esp_old->sh_sat_closed[i] * bf[v]->phi[i];
	}
    }

  v = LIGHT_INTP;
  fv_sens->poynt[0] = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->poynt[0] += *esp_old->poynt[0][i] * bf[v]->phi[i];
	}
    }

  v = LIGHT_INTM;
  fv_sens->poynt[1] = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->poynt[1] += *esp_old->poynt[1][i] * bf[v]->phi[i];
	}
    }

  v = LIGHT_INTD;
  fv_sens->poynt[2] = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->poynt[2] += *esp_old->poynt[2][i] * bf[v]->phi[i];
	}
    }

  v = SHELL_PRESS_OPEN;
  fv_sens->sh_p_open = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_p_open += *esp_old->sh_p_open[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_PRESS_OPEN_2;
  fv_sens->sh_p_open_2 = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_p_open_2 += *esp_old->sh_p_open_2[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_TEMPERATURE;
  fv_sens->sh_t = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_t += *esp_old->sh_t[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_DELTAH;
  fv_sens->sh_dh = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_dh += *esp_old->sh_dh[i] * bf[v]->phi[i];
	}
    }
  
  v = SHELL_LUB_CURV;
  fv_sens->sh_l_curv = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_l_curv += *esp_old->sh_l_curv[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_LUB_CURV_2;
  fv_sens->sh_l_curv_2 = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_l_curv_2 += *esp_old->sh_l_curv_2[i] * bf[v]->phi[i];
	}
    }
  
  v = SHELL_SAT_GASN;
  fv_sens->sh_sat_gasn = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sh_sat_gasn += *esp_old->sh_sat_gasn[i] * bf[v]->phi[i];
	}
    }
  v = POR_SINK_MASS;
  fv_sens->sink_mass = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->sink_mass += *esp_old->sink_mass[i] * bf[v]->phi[i];
	}
    }

  v = SHELL_SHEAR_TOP;
  fv_sens->sh_shear_top = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_shear_top += *esp_old->sh_shear_top[i] * bf[v]->phi[i];
        }
    }

  v = SHELL_SHEAR_BOT;
  fv_sens->sh_shear_bot = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_shear_bot += *esp_old->sh_shear_bot[i] * bf[v]->phi[i];
        }
    }

  v = SHELL_CROSS_SHEAR;
  fv_sens->sh_cross_shear = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->sh_cross_shear += *esp_old->sh_cross_shear[i] * bf[v]->phi[i];
        }
    }

  v = MAX_STRAIN;
  fv_sens->max_strain = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->max_strain += *esp_old->max_strain[i] * bf[v]->phi[i];
        }
    }

  v = CUR_STRAIN;
  fv_sens->cur_strain = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->cur_strain += *esp_old->cur_strain[i] * bf[v]->phi[i];
        }
    }

  /*
   * Pressure...
   */

  v = PRESSURE;
  fv_sens->P = 0.;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->P += *esp_old->P[i] * bf[v]->phi[i];
	}
    }

	      
  /*
   * Mesh displacement (vector)...
   * and positions (vector)
   */
  
  /*
   * Insure default value of 3rd coordinate is zero for low
   * dimensional problems...DIM=3, while dim can be less...
   */

  for ( p=dim; p<DIM; p++)
    {
      fv_sens->x[p]     = 0.;
      fv_sens->d[p]     = 0.;
    }

  for ( p=0; p<dim; p++)
    {
      v = MESH_DISPLACEMENT1 + p;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  fv_sens->x[p] = 0.;
	  fv_sens->d[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->d[p] += *esp_old->d[p][i] * bf[v]->phi[i];
	      fv_sens->x[p] += *esp_old->d[p][i] * bf[v]->phi[i];
	    }
	}
    }

  /*
   * SOLID displacement (vector)...
   * and positions (vector)
   */
  
  for ( p=dim; p<DIM; p++)
    {
      fv_sens->d_rs[p]     = 0.;
    }

  for ( p=0; p<dim; p++)
    {
      v = SOLID_DISPLACEMENT1 + p;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  fv_sens->d_rs[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->d_rs[p] += *esp_old->d_rs[p][i] * bf[v]->phi[i];
	    }
	}
    }

  status = load_coordinate_scales(pd->CoordinateSystem, fv);
  EH(status, "load_coordinate_scales(fv)");

  /*
   * Velocity (vector)...
   */
  
  velodim = dim;		/* Later, this might include v_theta... */
  if(pd->CoordinateSystem == SWIRLING ||
     pd->CoordinateSystem == PROJECTED_CARTESIAN)
    velodim = dim + 1;		/* Later is Now!  Woo!!! */

  /*
   * Default: all velocities are zero...
   */
  for ( p=velodim; p<DIM; p++)
    {
      fv_sens->v[p]     = 0.;
    }

  for ( p=0; p<velodim; p++)
    {
      v = VELOCITY1 + p;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  fv_sens->v[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->v[p] += *esp_old->v[p][i] * bf[v]->phi[i];

	    }
	}
      else
	{
	  fv_sens->v[p] = 0.;
	}
    }

  /* MMH
   * Particle velocity (vector)...
   */
  
  for ( p=0; p<velodim; p++)
    {
 
      v = PVELOCITY1 + p;
      fv_sens->pv[p]     = 0.;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->pv[p] += *esp_old->pv[p][i] * bf[v]->phi[i];
	    }
	}
    }

  /* 
   * Extension velocity (vector)...
   */
  
  v = EXT_VELOCITY;
  fv_sens->ext_v     = 0.;
  if ( pd->v[v] )
    {
      dofs     = ei->dof[v];
      for ( i=0; i<dofs; i++)
	{
	  fv_sens->ext_v += *esp_old->ext_v[i] * bf[v]->phi[i];
	}
    }


  /* 
   *  Level set normal vecotr
   */
  
  for ( p=0; p<velodim; p++)
    {
 
      v = NORMAL1 + p;
      fv_sens->n[p]     = 0.;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->n[p] += *esp_old->n[p][i] * bf[v]->phi[i];
	    }
	}
    }


  /*
   * Shell surface diffusion flux ...
   */

  v = SHELL_DIFF_FLUX;
  fv_sens->ext_v     = 0.;
  if ( pd->v[v] )
    {
      dofs     = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->ext_v += *esp_old->sh_J[i] * bf[v]->phi[i];
        }
    }
 
 
  /*
   * Shell surface curvature with normal unknowns ...
   */
 
  v = SHELL_DIFF_CURVATURE;
  fv_sens->ext_v     = 0.;
  if ( pd->v[v] )
    {
      dofs     = ei->dof[v];
      for ( i=0; i<dofs; i++)
        {
          fv_sens->ext_v += *esp_old->sh_Kd[i] * bf[v]->phi[i];
        }
    }


  /*
   *  Shell surface normal vecotr
   */

  for ( p=0; p<pd->Num_Dim; p++)
    {

      v = SHELL_NORMAL1 + p;
      fv_sens->n[p]     = 0.;
      if ( pd->v[v] )
        {
          dofs     = ei->dof[v];
          for ( i=0; i<dofs; i++)
            {
              fv_sens->n[p] += *esp_old->n[p][i] * bf[v]->phi[i];
            }
        }
    }


  /* 
   * Electric Field (vector)...
   */
  
  for ( p=0; p<velodim; p++)
    {
 
      v = EFIELD1 + p;
      fv_sens->pv[p]     = 0.;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->E_field[p] += *esp_old->E_field[p][i] * bf[v]->phi[i];
	    }
	}
    }

  /*
   * Polymer Stress (tensor)...
   */
  for ( mode=0; mode<vn->modes; mode++)
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      if( p <= q)
		{
		  /* good default behavior */
		  fv_sens->S[mode][p][q] = 0.;
		  
		  v = v_s[mode][p][q];
		  if ( pd->v[v] )
		    {
		      dofs     = ei->dof[v];
		      for ( i=0; i<dofs; i++)
			{
			  fv_sens->S[mode][p][q] += *esp_old->S[mode][p][q][i] * bf[v]->phi[i];
			}
		    }
		  /* form the entire symmetric stress matrix for the momentum equation */
		  fv_sens->S[mode][q][p] = fv_sens->S[mode][p][q];
		}
	    }
	}
    }


  /*
   * Velocity Gradient (tensor)...
   */
  
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  /* good default behavior */
	  fv_sens->G[p][q] = 0.;

	  v = v_g[p][q];
	  if ( pd->v[v] )
	    {
	      dofs     = ei->dof[v];
	      for ( i=0; i<dofs; i++)
		{
		  fv_sens->G[p][q] += *esp_old->G[p][q][i] * bf[v]->phi[i];
		}
	    }
	}
    }

	      
  /*
   * Concentration...
   */

  v = MASS_FRACTION;
  if (pd->v[v]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      dofs = ei->dof[v];
      fv_sens->c[w] = 0.0;
      for (i = 0; i < dofs; i++) {
	fv_sens->c[w] += *esp_old->c[w][i] * bf[v]->phi[i];
      }
    }
  } else {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      fv_sens->c[w] = 0.0;
    }
  }
	
  /*
   * External...
   */

  v = EXTERNAL;
  if ( efv->ev )
    {
      for ( w=0; w<efv->Num_external_field; w++)
	{
	  dofs     = ei->dof_ext[w];
	  fv_sens->external_field[w] = 0.;
  	  if( efv->i[w] != I_TABLE )
  	  {
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->external_field[w] += *evp->external_field[w][i] * bfex[w]->phi[i];
	    }
  	  }

          /* 
           * If the variable name is velocity, and the momentum equations are not active,
           * load the external_fields into the velocity fv for use in Advection-diffusion analysis
           */
          if (strcmp(efv->name[w], "VX") == 0 && !pd->v[VELOCITY1])
	    {
	      fv_sens->v[0] = fv_sens->external_field[w];
	    }
          if (strcmp(efv->name[w], "VY") == 0 && !pd->v[VELOCITY2])
	    {
	      fv_sens->v[1] = fv_sens->external_field[w];
	    }
          if (strcmp(efv->name[w], "VZ") == 0 && !pd->v[VELOCITY3])
	    {
	      fv_sens->v[2] = fv_sens->external_field[w];
	    }
          /* 
           * If the variable name is mesh displacement, and the mesh equations are not active,
           * load the external_fields into the mesh fv
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
	}
    }
  else
    {
      for ( w=0; w<efv->Num_external_field; w++)
	{
	  fv_sens->external_field[w]     = 0.;
	}
    }

  /* initial displacements for TALE anneals. all for KINEMATIC DISPLACEMENT BC */
  if ( efv->TALE )
    {
      for ( w=0; w<dim; w++)
	{
	  v = MESH_DISPLACEMENT1 + w;
	  dofs     = ei->dof[v];
	  fv_sens->initial_displacements[w] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->initial_displacements[w] += *evp->initial_displacements[w][i] * bf[v]->phi[i];

	    }
	  v = SOLID_DISPLACEMENT1 + w;
	  dofs     = ei->dof[v];
	  fv_sens->initial_displacements[w+DIM] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->initial_displacements[w+DIM] += *evp->initial_displacements[w+DIM][i] * bf[v]->phi[i];

	    }
	}
    }
  else
    {
      for ( w=0; w<dim; w++)
	{
	  fv_sens->initial_displacements[w]     = 0.;
	  fv_sens->initial_displacements[w + DIM]     = 0.;
	}
    }
	      
  return(status);
}  /*   end of load_fv_sens  */

/* load_fv_grads_sens() -- load relevant field variable sensitivity
 *		gradients at this gauss pt
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
 *			grad_c -- species concentration gradient
 *			grad_F -- fill gradient
 *			grad_V -- voltage potential gradient
 *			div_v  -- divergence of velocity 
 *			grad_v -- velocity gradient tensor
 *                      curl_v -- curl of velocity, a.k.a. vorticity
 *			div_d  -- divergence of displacement ( dilatation )
 *			grad_d -- gradient of displacement ( "strain" )
 *			grad_S -- gradient of the polymer stress tensor
 *			grad_G -- gradient of the velocity gradient tensor
 *                     grad_pv -- gradient of particle velocity
 *
 * Return values:
 *		0 -- if things went OK
 *	       -1 -- if a problem occurred
 *
 * Created:	Fri Mar 18 07:36:14 MST 1994 pasacki@sandia.gov
 *
 * Modified:	Tue Feb 21 11:08 MST 1995 pasacki@sandia.gov
 */

static int 
load_fv_grads_sens(void)
{
  int v;			/* variable type indicator */
  int i, a;			/* index */
  int p, q, r, s;		/* dimension index */
  int dofs;			/* degrees of freedom for a var in the elem */
  int w;			/* concentration species counter */
  const int dim = pd->Num_Dim;
  int status;
  int mode;

  status = 0;

  /*
   * grad(T)
   */
  v = TEMPERATURE;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_T[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_T[p] += *esp_old->T[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
	      
  /*
   * grad(P)
   */
  v = PRESSURE;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_P[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_P[p] += *esp_old->P[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
	      
  /*
   * grad(F)
   */
  v = FILL;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_F[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_F[p] += *esp_old->F[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }

  /*
   * grad(V)
   */
  v = VOLTAGE;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_V[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_V[p] += *esp_old->V[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }

  /*
   * grad(SH)
  v = SHEAR_RATE;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_SH[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_SH[p] += *esp_old->SH[i] * bf[v]->grad_phi[i][p];
	    }
	}
    }
   */


  /*
   * grad(d)
   */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  fv_sens->grad_d[p][q] = 0.;
	  v = MESH_DISPLACEMENT1 + p;
	  if ( pd->v[v] )
	    {
	      dofs     = ei->dof[v];
	      for ( r=0; r<VIM; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_d[p][q] += 
			*esp_old->d[r][i] * bf[v]->grad_phi_e[i][r] [p][q];

		    }
		}
	    }
	}
    }

  /*
   * grad(d_rs)
   */
  for ( p=0; p<VIM; p++)
    {
      for ( q=0; q<VIM; q++)
	{
	  fv_sens->grad_d_rs[p][q] = 0.;
	  v = SOLID_DISPLACEMENT1 + p;
	  if ( pd->v[v] )
	    {
	      dofs     = ei->dof[v];
	      for ( r=0; r<VIM; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_d_rs[p][q] += 
			*esp_old->d_rs[r][i] * bf[v]->grad_phi_e[i][r] [p][q];

		    }
		}
	    }
	}
    }
	      
  /*
   * div(d)
   */
  fv_sens->div_d = 0.;
  for ( p=0; p<VIM; p++)
    {
      fv_sens->div_d += fv_sens->grad_d[p][p];
    }

  /*
   * div(d_rs)
   */
  fv_sens->div_d_rs = 0.;
  for ( p=0; p<VIM; p++)
    {
      fv_sens->div_d_rs += fv_sens->grad_d_rs[p][p];
    }

  /*
   * grad(v)
   */
  dofs     = ei->dof[VELOCITY1];
  v = VELOCITY1;
  if(pd->v[v])
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      fv_sens->grad_v[p][q] = 0.;
	      for ( r=0; r<VIM; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_v[p][q] += 
			*esp_old->v[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    }


	      
  /* MMH
   * grad(pv), particle velocity gradients.
   */
  dofs = ei->dof[PVELOCITY1];
  v = PVELOCITY1;
  if(pd->v[v])
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      fv_sens->grad_pv[p][q] = 0.0;
	      for ( r=0; r<VIM; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_pv[p][q] += 
			*esp_old->pv[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    }

 
  /* 
   * grad(ext_v), extension velocity gradients.
   */
  dofs = ei->dof[EXT_VELOCITY];
  v = EXT_VELOCITY;
  if(pd->v[v])
    {
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_ext_v[p]= 0.0;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_ext_v[p] += 
		*esp_old->ext_v[i] * bf[v]->grad_phi[i][p];
	    }
	}
    }


 /* 
   * grad(E_field), Electric Field gradients.
   */
  v = EFIELD1;
  dofs = ei->dof[v];
  if(pd->v[v])
    {
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      fv_sens->grad_E_field[p][q] = 0.0;
	      for ( r=0; r<VIM; r++)
		{
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_E_field[p][q] += 
			*esp_old->E_field[r][i] * bf[v]->grad_phi_e[i][r] [p][q];
		    }
		}
	    }
	}
    }

  /*
   * grad(apr,api)
   */
  v = ACOUS_PREAL;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_apr[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_apr[p] += *esp_old->apr[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  v = ACOUS_PIMAG;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_api[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_api[p] += *esp_old->api[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
  v = ACOUS_REYN_STRESS;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_ars[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_ars[p] += *esp_old->ars[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
	      
  v = SHELL_BDYVELO;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->grad_sh_bv[p] = 0.;
	  for ( i=0; i<dofs; i++)
	    {
	      fv_sens->grad_sh_bv[p] += *esp_old->sh_bv[i] * bf[v]->grad_phi[i] [p];
	    }
	}
    }
	      
  v = SHELL_LUBP;
  if ( pd->v[v] )
    {
      dofs  = ei->dof[v];
      for ( p=0; p<VIM; p++)
      {
        fv_sens->grad_sh_p[p] = 0.;
        for ( i=0; i<dofs; i++)
          {
            fv_sens->grad_sh_p[p] += *esp_old->sh_p[i] * bf[v]->grad_phi[i] [p];
          }
      }
    }

  
  /* MMH
   * curl(v)
   */
  if( CURL_V != -1 )
  {
  v = VELOCITY1;
  dofs = ei->dof[VELOCITY1];
  if(pd->v[v])
    {
      for(p = 0; p < VIM; p++)
	{
	  fv_sens->curl_v[p] = 0.0;
	  for(i = 0; i < dofs; i++)
	    for(a = 0; a < VIM; a++)
	     fv_sens->curl_v[p] += *esp_old->v[a][i] * bf[v]->curl_phi_e[i][a][p];
	}
    }
  }

	      
  /*
   * div(v)
   */
  fv_sens->div_v = 0.;
  v = VELOCITY1;
  if ( pd->v[v] )
    {
 
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->div_v += fv_sens->grad_v[p][p];
	}
    }


  /* MMH
   * div(pv)
   */
  fv_sens->div_pv = 0.;
  v = PVELOCITY1;
  if ( pd->v[v] )
    {
      for ( p=0; p<VIM; p++)
	{
	  fv_sens->div_pv += fv_sens->grad_pv[p][p];
	}
    }


  /*
   * Concentration...
   */

  for (w = 0; w < pd->Num_Species_Eqn; w++)
    {
      v = MASS_FRACTION;
      if ( pd->v[v] )
	{
	  dofs     = ei->dof[v];
	  for ( p=0; p<VIM; p++ )
	    {
	      fv_sens->grad_c[w][p] = 0.;
	      for ( i=0; i<dofs; i++)
		{
		  fv_sens->grad_c[w][p] += 
		    *esp_old->c[w][i] * bf[v]->grad_phi[i][p];
		}
	    }
	}
    }
	      
  /*
   * grad(S), mode 0
   */
  v = POLYMER_STRESS11;
  if ( pd->v[v] )
    {  
      dofs = ei->dof[v];
      for ( mode=0; mode<vn->modes; mode++)
	{
	  for ( p=0; p<VIM; p++)
	    {
	      for ( q=0; q<VIM; q++)
		{
		  for ( r=0; r<VIM; r++)
		    {
		      fv_sens->grad_S[mode][r][p][q]=0.;
		      for ( i=0; i<dofs; i++)
			{
			  if(p<=q)
			    {
			      fv_sens->grad_S[mode][r][p][q] += 
				*esp_old->S[mode][p][q][i] * bf[v]->grad_phi[i][r] ;
			    }
			  else
			    {
			      fv_sens->grad_S[mode][r][p][q] += 
				*esp_old->S[mode][q][p][i] * bf[v]->grad_phi[i][r] ;
			    }
			}
		      
		    }
		}
	    }

	  /*
	   * div(S) - this is a vector!
	   */
	  for ( r=0; r<VIM; r++)
	    {
	      fv_sens->div_S[mode][r]=0.;  
	    }
	  
	  for ( r=0; r<dim; r++)
	    {
	      for ( q=0; q<dim; q++)
		{
		  fv_sens->div_S[mode][r] += 
		    fv_sens->grad_S[mode][q][q][r];
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
			  fv_sens->div_S[mode][s] += 
			    fv_sens->S[mode][p][s]*fv_sens->grad_e[p][r][s];
			}
		    }
		}
	      
	      for ( s=0; s<VIM; s++)
		{
		  for ( r=0; r<VIM; r++)
		    {
		      for ( q=0; q<VIM; q++)
			{ 
			  fv_sens->div_S[mode][s] += 
			    fv_sens->S[mode][r][q]* fv_sens->grad_e[q][r][s] ;
			}
		    }
		}
	    }
	}
      
    }
  
  
  
  /*
   * grad(G)
   */
  
  v = VELOCITY_GRADIENT11;
  if ( pd->v[v] )
    {  
      dofs = ei->dof[v];
      for ( p=0; p<VIM; p++)
	{
	  for ( q=0; q<VIM; q++)
	    {
	      for ( r=0; r<VIM; r++)
		{
		  fv_sens->grad_G[r][p][q]=0.;
		  for ( i=0; i<dofs; i++)
		    {
		      fv_sens->grad_G[r][p][q] += 
			*esp_old->G[p][q][i] * bf[v]->grad_phi[i][r];
		    }
		}
	    }
	}
      
      /*
       * div(G) - this is a vector!
       */
      for ( r=0; r<VIM; r++)
	{
	  fv_sens->div_G[r]=0.;  
	}
      
      for ( r=0; r<dim; r++)
	{
	  for ( q=0; q<dim; q++)
	    {
	      fv_sens->div_G[r] += 
		fv_sens->grad_G[q][q][r];
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
		      fv_sens->div_G[s] += 
			fv_sens->G[p][s]*fv_sens->grad_e[p][r][s];
		    }
		}
	    }
	  
	  for ( s=0; s<VIM; s++)
	    {
	      for ( r=0; r<VIM; r++)
		{
		  for ( q=0; q<VIM; q++)
		    { 
		      fv_sens->div_G[s] += 
			fv_sens->G[r][q]* fv_sens->grad_e[q][r][s] ;
		    }
		}
	    }
	}
    }

  return(status);
}/*  end of load_fv_grads_sens  */

/*
	ADAPTIVE INTEGRATION WEIGHT ROUTINE
 */


 
int adaptive_weight (
			double w[],
 			const int ngp,
 			const int dim,
 			const double ls_F[],
 			const double alpha,
 			const int wt_type,
 			const int elem_type
 		     )
{
double bf_mom[9];
double a = 0.0, b = 0.0;
int i,j;
/*int wt_type=1;       	1 - unity
 			2 - sharp Heaviside
 			3 - sharp delta function
 			4 - diffuse Heaviside
			5 - diffuse delta function
 		*/
 
int F_type = 0;	/*	type of interface crossing  */
int interpolation=2,nint = 0;
int sharp_interface=1;
int same_sign;
double f0, f1;
double xf[2];
double xf2D[6][2];
int nint2D[6];
int is2D = 0, side_id[8];
int side_diff, side_ct;
int return_val = 1;
double ecrd[12][MAX_PDIM];

#ifndef NO_CHEBYSHEV_PLEASE
int chev_order=3;
#endif
double gauss_wt1D[3]={5/9.,8/9.,5/9.};
double gauss_wt2D[9]={25/81.,40/81.,25/81.,40/81.,64/81.,40/81.,25/81.,
 			40/81.,25/81.};
struct interfaces 
	{
	int F_type;
	double coeff[6];
	double sdet_c[6];
	double endpts[2];
	int sense;
	int integration_dir;
	double bf_mom[9];
	double hfactor;
	} m_int[4];

/* ANSI C note: function calls (e.g., sqrt) not allowed in intializer constant expressions */
/* The following lines are retained for identification purposes only:                */
/* double bfi_1D[3][3]={                                                             */
/*  	0.5*sqrt(5/3)+5/6.,-0.5*sqrt(5/3)+5/6.,0,                                    */
/*  	-2/3.,-2/3.,1,                                                               */
/*  	-0.5*sqrt(5./3.)+5/6.,0.5*sqrt(5./3.)+5/6.,0                                 */
/*  	};                                                                           */
/* double bfi_2D[9][9]={                                                             */
/*  	5/18.*(4+sqrt(15.)),5/18.,-5/18.*(-4+sqrt(15.)),5/18.,0,0,0,0,0,             */
/*  	(-5-sqrt(15.))/9.,(-5-sqrt(15.))/9.,(-5+sqrt(15.))/9.,(-5+sqrt(15.))/9.,(5+sqrt(15.))/6.,0,(5-sqrt(15.))/6.,0,0, */
/*  	5/18.,5/18.*(4+sqrt(15.)),5/18.,-5/18.*(-4+sqrt(15.)),0,0,0,0,0,             */
/*  	(-5-sqrt(15.))/9.,(-5+sqrt(15.))/9.,(-5+sqrt(15.))/9.,(-5-sqrt(15.))/9.,0,(5-sqrt(15.))/6.,0,(5+sqrt(15.))/6.,0, */
/*  	4/9.,4/9.,4/9.,4/9.,-2/3.,-2/3.,-2/3.,-2/3.,1,                               */
/*  	(-5+sqrt(15.))/9.,(-5-sqrt(15.))/9.,(-5-sqrt(15.))/9.,(-5+sqrt(15.))/9.,0,(5+sqrt(15.))/6.,0,(5-sqrt(15.))/6.,0, */
/*  	5/18.,-5/18.*(-4+sqrt(15.)),5/18.,5/18.*(4+sqrt(15.)),0,0,0,0,0,             */
/*  	(-5+sqrt(15.))/9.,(-5+sqrt(15.))/9.,(-5-sqrt(15.))/9.,(-5-sqrt(15.))/9.,(5-sqrt(15.))/6.,0,(5+sqrt(15.))/6.,0,0, */
/*  	-5/18.*(-4+sqrt(15.)),5/18.,5/18.*(4+sqrt(15.)),5/18.,0,0,0,0,0              */
/*  	};                                                                           */

double bfi_1D[3][3]={
 	{1.478830557701236147529878,0.1878361089654305191367891,0},
  	{-2/3.,-2/3.,1},
   	{0.1878361089654305191367891,1.478830557701236147529878,0}
  	};
double bfi_2D[9][9]={
	{2.1869398183909491347720182,5/18.,0.03528240383127308745020406,5/18.,0,0,0,0,0},
 	{-0.9858870384674907650199184,-0.9858870384674907650199184,
          -0.12522407264362034609119273,-0.1252240726436203460911927,
           1.4788305577012361475298776,0,0.1878361089654305191367891,0,0},
 	{5/18.,2.186939818390949134772018,5/18.,0.03528240383127308745020406,0,0,0,0,0},
 	{-0.9858870384674907650199184,-0.1252240726436203460911927,
          -0.1252240726436203460911927,-0.9858870384674907650199184,0,
          0.1878361089654305191367891,0,1.478830557701236147529878,0},
 	{4/9.,4/9.,4/9.,4/9.,-2/3.,-2/3.,-2/3.,-2/3.,1},
 	{-0.1252240726436203460911927,-0.9858870384674907650199184,
          -0.9858870384674907650199184,-0.1252240726436203460911927,0,
          1.478830557701236147529878,0,0.1878361089654305191367891,0},
 	{5/18.,0.03528240383127308745020406,5/18.,2.186939818390949134772018,0,0,0,0,0},
 	{-0.1252240726436203460911927,-0.1252240726436203460911927,
          -0.9858870384674907650199184,-0.9858870384674907650199184,
          0.1878361089654305191367891,0,1.478830557701236147529878,0,0},
 	{0.03528240383127308745020406,5/18.,2.186939818390949134772018,5/18.,0,0,0,0,0}
  	};  

int dupl_side = 0, dupl_id = 0, side1 = -1, side2 = -1;
double int_angle[8], xloc;

 
if( elem_type != BIQUAD_QUAD)
 	{EH(-1,"adaptive integration for 2D quads only!");}

#ifndef NO_CHEBYSHEV_PLEASE 
chev_order = ls->Adaptive_Order;
#endif
 
if( sharp_interface)
   {
 
/*	check to see if interface crosses the element  */
 
i=1;
do {
 	same_sign = SGN(ls_F[i]) == SGN(ls_F[i-1]) ? 1 : 0;
 	i++;
 	}	while (i < pow(3,dim) && same_sign);
if(same_sign)
 	{
 	if(SGN(ls_F[0]) == 1)
 		{
 		switch (dim)
 			{
 			case 1:
 				for(j=0;j<ngp;j++)	{w[j]=gauss_wt1D[j];}
 				break;
 			case 2:
 				for(j=0;j<ngp;j++)	{w[j]=gauss_wt2D[j];}
 				break;
 			}
 		}
 	else
 		{
 		for(j=0;j<ngp;j++)	{w[j]=0;}
 		}
 	return(2);
 	}
 
switch (dim)
{
case 1:
/*  determine interface location  */
 
switch (interpolation)
 	{
 	case 1:
 		f0 = 0.5*(ls_F[1]+ls_F[0]);
 		f1 = 0.5*(ls_F[1]-ls_F[0]);
 		xf[0] = -f0/f1;
 		break;
 	case 2:
		nint = interface_crossing_1DQ( ls_F, xf );
 		break;
 	default:
 		printf("unknown interpolation %d \n",interpolation);
 	}  /* end of interpolation switch */
 
/*  determine how the LS fcn intersects the element  */
 
if(nint == 1)
{
if( fabs(xf[0]) <= 1.0)	
 	{
 	if( (ls_F[1]-ls_F[0]) >= 0.0)
 		{F_type = 1;}
 		else
 		{F_type = 2;}
 	}
 	else
 	{
 	if( ls_F[2] >= 0.0)
 		{F_type = 4;}
 		else
 		{F_type = 3;}
 	}
}
else if (nint == 2)
{
if( ls_F[0] <= 0.0)	
 	{
 	F_type = 5;
 	}
else
 	{
 	F_type = 6;
 	}
}
 
switch (wt_type)
 	{
 	case 1:
 		bf_mom[0] = 1./3.;
 		bf_mom[1] = 1./3.;
 		bf_mom[2] = 4./3.;
 		break;
 	case 2:
 		switch (F_type)
 			{
 			case 1:
 				a = xf[0];  b = 1.;
 				break;
 			case 2:
 				a = -1.;  b = xf[0];
 				break;
 			case 3:
 				a = 0.;  b = 0.;
 				break;
 			case 4:
 				a = -1;  b = 1.;
 				break;
 			case 5:
 			case 6:
 				a = xf[0];  b = xf[1];
 				break;
 			default:
 			printf("unknown F_type %d \n",F_type);
 			}  /* end of F_type switch */
 		bf_mom[0] = (a*a*(3.-2.*a)+b*b*(2.*b-3.))/12.;
 		bf_mom[1] = (a*a*(-3.-2.*a)+b*b*(2.*b+3.))/12.;
 		bf_mom[2] = (a*a*a - b*b*b)/3. +b -a;
		if(F_type == 6 )
			{
			bf_mom[0] = 1./3. - bf_mom[0];
			bf_mom[1] = 1./3. - bf_mom[1];
			bf_mom[2] = 4./3. - bf_mom[2];
			}
 		break;
 	case 3:
 		switch (F_type)
 			{
 			case 1:
 			case 2:
 				bf_mom[0] = xf[0]*(xf[0]-1)*0.5;
 				bf_mom[1] = xf[0]*(xf[0]+1)*0.5;
 				bf_mom[2] = 1-xf[0]*xf[0];
 				break;
 			case 3:
 			case 4:
 				bf_mom[0] = 0;
 				bf_mom[1] = 0;
 				bf_mom[2] = 0;
 				break;
 			case 5:
 			case 6:
 				bf_mom[0] = xf[0]*(xf[0]-1)*0.5;
 				bf_mom[1] = xf[0]*(xf[0]+1)*0.5;
 				bf_mom[2] = 1-xf[0]*xf[0];
 				bf_mom[0] += xf[1]*(xf[1]-1)*0.5;
 				bf_mom[1] += xf[1]*(xf[1]+1)*0.5;
 				bf_mom[2] += 1-xf[1]*xf[1];
 				break;
 			default:
 			printf("unknown F_type %d \n",F_type);
 			}
 		break;
 	default:
 		printf("unknown weight fcn type %d \n",wt_type);
 	}  /* end of wt_type switch */
break;  /* end of dim switch, case 1 */

case 2:
/*  determine interface location for 2D elements */
 
switch (interpolation)
 	{
 	case 1:
 		f0 = 0.5*(ls_F[1]+ls_F[0]);
 		f1 = 0.5*(ls_F[1]-ls_F[0]);
 		xf[0] = -f0/f1;
 		break;
 	case 2:
 		is2D = interface_crossing_2DQ( ls_F, xf2D, side_id, nint2D, ecrd);
 		break;
 	default:
 		printf("unknown interpolation %d \n",interpolation);
 	}  /* end of interpolation switch */

/*  determine how the LS fcn intersects the element  */
 
/*	some initialization			*/
for(i=0 ; i<is2D/2 ; i++)
	{
	m_int[i].sense = FALSE;
	m_int[i].integration_dir = 0;
	m_int[i].hfactor = 1;
	m_int[i].F_type = 0;
	}

if(is2D == 2)
 	{
 	side_diff=abs(side_id[1]-side_id[0]);
 	switch (side_diff)
 		{
 		case 1:
 		case 3:
 			m_int[0].F_type = MAX(side_id[0],side_id[1])+side_diff;
			break;
		case 2:
			m_int[0].F_type = 10*MAX(side_id[0],side_id[1])+side_diff;
			break;
		case 0:
			m_int[0].F_type = 100+side_id[0];
			break;
		default:
			printf("shouldn't get here - side_diff switch\n");
			return_val = -1;
			break;
		}  /* end of side_diff */
	}
else if(is2D == 0)
	{
	m_int[0].sense = FALSE;
	m_int[0].integration_dir = 0;
	m_int[0].hfactor = 1;
	m_int[0].F_type = 0;
	if(ls_F[0] < 0.) m_int[0].sense=TRUE;
	}
else if (is2D == 4)
	{
	interface_inclination_2DQ( ls_F, is2D, int_angle, side_id, xf2D);

	fprintf(stderr,"two interface intersections %d\n",is2D);
	side_ct = 0;
	for(i=1;i<is2D;i++)	
		{
		if(side_id[i] == side_id[i-1])
			{
			dupl_side = side_id[i];
			dupl_id = i;
			}
			else
			{
			side_ct++;
			}
		}
	for(i=0;i<is2D;i++)	
		{
		j=side_id[i];
		xloc = xf2D[j][0];
		if( i>0 && j == side_id[i-1] )
			{	xloc = xf2D[j][1];	}
		}

	if( fabs(cos(M_PIE*(int_angle[dupl_id]-90*dupl_side)/180)) < 0.7 ||
		fabs(cos(M_PIE*(int_angle[dupl_id-1]-90*dupl_side)/180)) <  0.7 )
		fprintf(stderr,"\n  Adapt_int - duplicate side intersection(s) not near tangent. \n");
 	switch (side_ct)
 		{
		case 1:
			m_int[0].F_type = 100+side_id[0];
			m_int[1].F_type = 100+side_id[is2D-1];
			break;
		case 2:
			switch (dupl_id)
			    {
			     case 1:
				side1 = dupl_id+2;
				side2 = dupl_id+1;
				break;
			     case 2:
				side1 = dupl_id-2;
				side2 = dupl_id+1;
				break;
			     case 3:
				side1 = dupl_id-2;
				side2 = dupl_id-3;
				break;
				}
/* first intersection   */
			side_diff=abs(side_id[dupl_id-1]-side_id[side1]);
 				switch (side_diff)
 				   {
 					case 1:
 					case 3:
	m_int[0].F_type = MAX(side_id[dupl_id-1],side_id[side1]) +side_diff;
						break;
					case 2:
	m_int[0].F_type = 10*MAX(side_id[dupl_id-1],side_id[side1])+side_diff;
						break;
					default:
			printf("shouldn't get here - side_diff switch\n");
					return_val = -1;
						break;
				    }  /* end of side_diff */
/* second intersection   */
			side_diff=abs(side_id[dupl_id]-side_id[side2]);
 				switch (side_diff)
 				   {
 					case 1:
 					case 3:
	m_int[1].F_type = MAX(side_id[dupl_id],side_id[side2]) +side_diff;
						break;
					case 2:
	m_int[1].F_type = 10*MAX(side_id[dupl_id],side_id[side2])+side_diff;
						break;
					default:
			printf("shouldn't get here - side_diff switch\n");
					return_val = -1;
						break;
				    }  /* end of side_diff */
			break;
		case 3:
			if( (int_angle[0] < 90 && int_angle[0] >=0 ) ||
				(int_angle[0] <= -90 && int_angle[0] >= -180))
				{
				m_int[0].F_type = 2; m_int[1].F_type = 4;
				}
			else
				{
				m_int[0].F_type = 6; m_int[1].F_type = 3;
				}
			break;
		}
	}
else if (is2D == 6)
	{
	fprintf(stderr,"three interface intersections %d\n",is2D);
	interface_inclination_2DQ( ls_F, is2D, int_angle, side_id, xf2D);

	for(i=0;i<is2D;i++)	
		{
		j=side_id[i];
		xloc = xf2D[j][0];
		if( i>0 && j == side_id[i-1] )
			{	xloc = xf2D[j][1];	}
		fprintf(stderr,"%d %g %g %d\n",j,xloc,cos(M_PIE*(int_angle[i]-90*j)/180),nint2D[j]);
		}

 	for(j=0;j<ngp;j++)	{w[j]=gauss_wt2D[j];}
	return(2);
	}
else if (is2D == 8)
	{
	fprintf(stderr,"four interface intersections %d\n",is2D);
	interface_inclination_2DQ( ls_F, is2D, int_angle, side_id, xf2D);

	for(i=0;i<is2D;i++)	
		{
		j=side_id[i];
		xloc = xf2D[j][0];
		if( i>0 && j == side_id[i-1] )
			{	xloc = xf2D[j][1];	}
		fprintf(stderr,"%d %g %g %d\n",j,xloc,cos(M_PIE*(int_angle[i]-90*j)/180),nint2D[j]);
		}

 	for(j=0;j<ngp;j++)	{w[j]=gauss_wt2D[j];}
	return(2);
	}
else
	{
	fprintf(stderr,"odd number of intersections %d\n",is2D);
 	for(j=0;j<ngp;j++)	{w[j]=gauss_wt2D[j];}
	return(2);
	}

/* determine coefficients for interface representation	*/

for(i=0 ; i<is2D/2 ; i++)
	  {
#ifndef NO_CHEBYSHEV_PLEASE
	m_int[i].integration_dir = 
	chebyshev_coeff_2DQ( m_int[i].F_type, ls_F, xf2D, 
		chev_order , m_int[i].coeff , m_int[i].endpts,
		 &m_int[i].sense, nint2D);
#else	
	 EH(-1,"Turn off NO_CHEBYSHEV_PLEASE please.\n");
#endif
	  }

switch (wt_type)
	{
	case 1:
		bf_mom[0] = 1./9.;
		bf_mom[1] = 1./9.;
		bf_mom[2] = 1./9.;
		bf_mom[3] = 1./9.;
		bf_mom[4] = 4./9.;
		bf_mom[5] = 4./9.;
		bf_mom[6] = 4./9.;
		bf_mom[7] = 4./9.;
		bf_mom[8] = 16./9.;
		break;
	case 2:
		memset(bf_mom, 0, sizeof(double)*9);

/*  Try to insure correct addition for multiple heaviside integrals  */

		if( is2D == 4 )
			{
			if(m_int[0].sense && m_int[1].sense)
				{
			m_int[1].sense = FALSE;
			m_int[1].hfactor = -1;
				}
			else if( (m_int[0].sense && !m_int[1].sense) &&
			         (m_int[1].F_type > 10 && m_int[1].F_type <90)) 
				{
				m_int[1].sense = TRUE;
				m_int[1].hfactor = -1;
				}
			else if( (!m_int[0].sense && m_int[1].sense) &&
			         (m_int[0].F_type > 10 && m_int[0].F_type <90)) 
				{
				m_int[0].sense = TRUE;
				m_int[0].hfactor = -1;
				}
			}
		for(i=0 ; i<is2D/2 ; i++)
	  		{
#ifndef NO_CHEBYSHEV_PLEASE
			 heaviside_chev_moments_2DQ( m_int[i].F_type, 
				m_int[i].sense, m_int[i].bf_mom , 
				m_int[i].coeff, chev_order, m_int[i].endpts );
#else
			EH(-1,"Turn off NO_CHEBYSHEV_PLEASE please\n");
#endif				
			for(j=0 ; j<9 ; j++)	
				{bf_mom[j] += m_int[i].hfactor*m_int[i].bf_mom[j];}
			}
		break;
 	case 3:
		memset(bf_mom, 0, sizeof(double)*9);
		for(i=0 ; i<is2D/2 ; i++)
	  		{
#ifndef NO_CHEBYSHEV_PLEASE
			 surfdet_chev_coeff_2DQ( m_int[i].F_type, ls_F,
				m_int[i].coeff, m_int[i].sdet_c , chev_order,
				m_int[i].integration_dir, m_int[i].endpts);

			 delta_chev_moments_2DQ( m_int[i].F_type, 
				m_int[i].bf_mom , m_int[i].coeff, 
				m_int[i].sdet_c , chev_order, m_int[i].endpts);
#else
			EH(-1,"Turn off NO_CHEBYSHEV_PLEASE please at compile time.\n");
#endif				
			for(j=0 ; j<9 ; j++)	
				{bf_mom[j] += m_int[i].bf_mom[j];}
	  		}
 		break;
	default:
		printf("unknown weight fcn type %d \n",wt_type);
	}  /* end wt-type switch */

  break;

}    /* end of dim switch, case 2 */ /* dimension switch  */
  }	else	{
	printf("diffuse interface not done yet\n");
	return(-1);
  } 

/*  compute weights  */

switch(dim)
	{
	case 1:
		for (i=0 ; i<ngp ; i++)	{
			w[i] = 0.0;
			for(j=0; j<ngp ; j++)	{
				w[i] += bfi_1D[i][j]*bf_mom[j];
				}
			}
		break;
	case 2:
		for (i=0 ; i<ngp ; i++)	{
			w[i] = 0.0;
			for(j=0; j<ngp ; j++)	{
				w[i] += bfi_2D[i][j]*bf_mom[j];
				}
			}
		break;
	default:
		printf("3D not done yet!  %d  \n",dim);
	}

return(return_val);
} /* end of function adaptive_weight */
/**********************************************************************/

int solve_quadratic( 
			const double c,
			const double b,
			const double a,
			double root[]
			)
{
double disc, temp;
int nroot;
		disc = b*b-4.*a*c;
		temp = -0.5*(b+SGN(b)*sqrt(disc));
		if(disc < 0.0)
			{ nroot=0; }
		else if(disc > 0.0)
			{
			root[0] = c/temp;
			if( a != 0.0)
				{
				nroot=2;
				root[1] = temp/a;
				}
				else
				{
				nroot=1;
				}
			}
		else
			{
			nroot=1;
			root[0] = c/temp;
			}
	return(nroot);
}

#ifndef NO_CHEBYSHEV_PLEASE
/*   computation of delta fcn integration moments using
	Chebyshev polynomial representation of interface		  */

int 
chebyshev_coeff_2DQ( 
			const int interface_type,
			const double ls_F[],
			double xf2D[6][2],
			int n_chev,
			double coeff[],
			double endpts[],
			int *sense,
			const int nint2D[]
			)
{
int i, j;
double a=-1, b=1;
double x, y, sum, factor;
double fval[MAX_CHEV];
double ca=0, cb=0, cm;
int idir=0;
double f0, f1, f2, xint[2]={0,0};


switch (interface_type)
 	{
	case 0:
		for(i=0;i<n_chev;i++)	{coeff[i]=0;}
		return(idir);
		break;
 	case 6:
		b=xf2D[0][0];
 		ca = xf2D[3][0];	 cb = -1;
 		break;
 	case 2:
		a=xf2D[0][0];
		if(nint2D[0] == 2) a=xf2D[0][1];
 		ca = -1;	cb = xf2D[1][0];
 		break;
 	case 4:
		b=xf2D[2][0];
 		ca = xf2D[3][0];
		if(nint2D[3] == 2) ca=xf2D[3][1];
		cb = 1;
 		break;
 	case 3:
		a=xf2D[2][0];
		if(nint2D[2] == 2) a=xf2D[2][1];
 		ca = 1;	
		cb = xf2D[1][0];
		if(nint2D[1] == 2) cb=xf2D[1][1];
 		break;
	case 22:
		idir = 1;
 		ca = xf2D[0][0];	cb = xf2D[2][0];
		if(nint2D[0] == 2 && nint2D[3])	ca = xf2D[0][1]; 
		if(nint2D[2] == 2 && nint2D[3])	cb = xf2D[2][1]; 
		break;
	case 32:
 		ca = xf2D[3][0];	cb = xf2D[1][0];
		if(nint2D[3] == 2 && nint2D[0])	ca = xf2D[3][1]; 
		if(nint2D[1] == 2 && nint2D[0])	cb = xf2D[1][1]; 
		break;
	case 101:
		idir = 1;
		a=xf2D[1][0]; b=xf2D[1][1];
 		ca = 1;	cb = 1;
		break;
	case 103:
		idir = 1;
		a=xf2D[3][0]; b=xf2D[3][1];
 		ca = -1;	cb = -1;
		break;
	case 100:
		a=xf2D[0][0]; b=xf2D[0][1];
 		ca = -1;	cb = -1;
		break;
	case 102:
		a=xf2D[2][0]; b=xf2D[2][1];
 		ca = 1;	cb = 1;
		break;
 	default:
 		EH(-1,"shouldn't get here - F_type switch \n");
 		break;
 	} /* end of F_type switch */

/*  compute Chebyshev coefficients  */

for (i=0 ; i<n_chev ; i++)
	{
	y = chevpoly[n_chev-3].root[i];
	x = 0.5*(y*(b-a)+b+a);
		/*  evaluate interface location at x in [-1,1]	*/
	if( idir )
		{
		f0 = ls_F[4]*0.5*x*(x-1) + ls_F[8]*(1-x*x) + ls_F[6]*0.5*x*(1+x);
		f1 = 0.5*(ls_F[1]*0.5*x*(x-1) + ls_F[5]*(1-x*x) + 
			ls_F[2]*0.5*x*(1+x) - (ls_F[0]*0.5*x*(x-1) + 
			ls_F[7]*(1-x*x) + ls_F[3]*0.5*x*(1+x)) );
		f2 = 0.5*(ls_F[1]*0.5*x*(x-1) + ls_F[5]*(1-x*x) + 
			ls_F[2]*0.5*x*(1+x) + ls_F[0]*0.5*x*(x-1) + 
			ls_F[7]*(1-x*x) + ls_F[3]*0.5*x*(1+x))- f0;
		}
	else
		{
 		f0 = ls_F[7]*0.5*x*(x-1) + ls_F[8]*(1-x*x) + ls_F[5]*0.5*x*(1+x);
 		f1 = 0.5*(ls_F[3]*0.5*x*(x-1) + ls_F[6]*(1-x*x) + 
 			ls_F[2]*0.5*x*(1+x) - (ls_F[0]*0.5*x*(x-1) + 
 			ls_F[4]*(1-x*x) + ls_F[1]*0.5*x*(1+x)) );
		f2 = 0.5*(ls_F[3]*0.5*x*(x-1) + ls_F[6]*(1-x*x) + 
 			ls_F[2]*0.5*x*(1+x) + ls_F[0]*0.5*x*(x-1) + 
 			ls_F[4]*(1-x*x) + ls_F[1]*0.5*x*(1+x)) - f0;
		}
	cm = 0.5*(y*(cb-ca)+cb+ca);
 	j = solve_quadratic(f0,f1,f2,xint);
 	switch (j)
		{
		case 0:
			EH(-1," no interior intersections found\n"); break;
		case 1:
			fval[i]=xint[0]; break;
		case 2:
			if(fabs(xint[0]) <= 1.0 && fabs(xint[1]) <= 1.0)
				{ printf(" two interior intersections found\n");}
			if(fabs(xint[0]-cm) < fabs(xint[1]-cm))
				{ fval[i]=xint[0]; }
			else
				{ fval[i]=xint[1]; }
			break;
		}
	}
factor = 2./n_chev;
for (j=0 ; j<n_chev ; j++)
	{
	sum=0.;
	for(i=0 ; i<n_chev ; i++)
		{
		sum += fval[i]*chevpoly[n_chev-3].cosval[n_chev*j+i];
		}
	coeff[j] = sum*factor;
	}

switch (interface_type)
 	{
 	case 6:
		if(ls_F[0] < 0.0) *sense = TRUE;
		break;
 	case 2:
		if(ls_F[1] < 0.0) *sense = TRUE;
		break;
 	case 4:
		if(ls_F[3] < 0.0) *sense = TRUE;
		break;
 	case 3:
		if(ls_F[2] < 0.0) *sense = TRUE;
		break;
 	case 22:
		if(ls_F[7] > ls_F[5]) *sense = TRUE;
				break;
	case 32:
		if(ls_F[4] > ls_F[6]) *sense = TRUE;
				break;
	case 100:
		f1 = 0.5*(a+b);
		f0 = ls_F[0]*0.5*f1*(f1-1) + ls_F[4]*(1-f1*f1) 
			+ ls_F[1]*0.5*f1*(1+f1);
		if(f0 < 0) *sense = TRUE;
		break;
	case 101:
		f1 = 0.5*(a+b);
		f0 = ls_F[1]*0.5*f1*(f1-1) + ls_F[5]*(1-f1*f1) 
			+ ls_F[2]*0.5*f1*(1+f1);
		if(f0 < 0) *sense = TRUE;
		break;
	case 102:
		f1 = 0.5*(a+b);
		f0 = ls_F[3]*0.5*f1*(f1-1) + ls_F[6]*(1-f1*f1) 
			+ ls_F[2]*0.5*f1*(1+f1);
		if(f0 < 0) *sense = TRUE;
		break;
	case 103:
		f1 = 0.5*(a+b);
		f0 = ls_F[0]*0.5*f1*(f1-1) + ls_F[7]*(1-f1*f1) 
			+ ls_F[3]*0.5*f1*(1+f1);
		if(f0 < 0) *sense = TRUE;
		break;
	}
endpts[0] = a;	endpts[1] = b;
return(idir);
}
/*  end of chebyshev_coeff_2DQ*/

#endif

/*   determine how interface crosses the element	  */

int
interface_crossing_1DQ( 
			const double ls_F[],
			double xf[2]
			)
{
int nint, n_intersections;
double f0, f1, f2, xint[2];

f0 = ls_F[2];
f1 = 0.5*(ls_F[1]-ls_F[0]);
f2 = 0.5*(ls_F[0]+ls_F[1])-ls_F[2];
nint = solve_quadratic(f0,f1,f2,xint);
n_intersections = nint;
switch (nint)
	{
 	case 0:
 	/* return(0);*/
		n_intersections = 0;
 		break;
 	case 1:
 		if(fabs(xint[0]) > 1.0) n_intersections = 0;
 			xf[0] = xint[0];
 		break;
 	case 2:
 		if(fabs(xint[0]) <= 1.0 && fabs(xint[1]) <= 1.0)
 			{
			if(xint[0] < xint[1])
				{xf[0] = xint[0];	xf[1] = xint[1];}
			else
				{xf[0] = xint[1];	xf[1] = xint[0];}
 			}
		else if(fabs(xint[0]) <= 1.0)
			{
 			xf[0] = xint[0];
			n_intersections = 1;
			}
		else if(fabs(xint[1]) <= 1.0)
			{
 			xf[0] = xint[1];
			n_intersections = 1;
			}
		else
			{
			n_intersections = 0;
			}
 		break;
 	}  /* end of nint switch */
return(n_intersections);
}
/*  end of interface_crossing_1DQ	*/



int
interface_crossing_2DQ( 
			const double ls_F[],
			double xf2D[6][2],
			int side_id[],
			int nint2D[],
			double coord[12][MAX_PDIM]
			)
{
int i, j, iside, is2D=0;
int begnode[6]={0,1,3,0,7,4}, endnode[6]={1,2,2,3,5,6},midnode[6]={4,5,6,7,8,8};
double f0, f1, f2;
double xint[2]={0.,0.};
double dist,point, dist_tol=1.0E-06;
int dupl[8]={0,0,0,0,0,0,0,0};
int nint;


for(iside=0 ; iside<4 ; iside++)	{
f0 = ls_F[midnode[iside]];
f1 = 0.5*(ls_F[endnode[iside]]-ls_F[begnode[iside]]);
f2 = 0.5*(ls_F[begnode[iside]]+ls_F[endnode[iside]])-ls_F[midnode[iside]];
nint = solve_quadratic(f0,f1,f2,xint);
nint2D[iside] = nint;
switch (nint)
 	{
 	case 0:
 		break;
 	case 1:
 		if(fabs(xint[0]) <= 1.0)
 			{
 			xf2D[iside][0]=xint[0];
 			side_id[is2D]=iside;
 			is2D++;
 			}
 		break;
 	case 2:
 		if(fabs(xint[0]) <= 1.0 && fabs(xint[1]) <= 1.0)
 			{
/*	printf(" two interface intersections found on a side\n");*/
 			if(xint[0] < xint[1])
 				{
 				xf2D[iside][0]=xint[0];
 				xf2D[iside][1]=xint[1];
 				}
 			else
 				{
 				xf2D[iside][0]=xint[1];
 				xf2D[iside][1]=xint[0];
 				}
 			side_id[is2D]=iside;
 			is2D++;
 			side_id[is2D]=iside;
 			is2D++;
 			}
 		else if(fabs(xint[0]) <= 1.0)
 			{
 			xf2D[iside][0]=xint[0];
 			side_id[is2D]=iside;
 			is2D++;
			nint2D[iside] = 1;
 			}
 		else if(fabs(xint[1]) <= 1.0)
 			{
 			xf2D[iside][0]=xint[1];
 			side_id[is2D]=iside;
 			is2D++;
			nint2D[iside] = 1;
 			}
 		break;
 	}  /* end of nint2D[iside] switch */
} /* end of for iside loop  */


/**  check for interface duplication from corner nodes  **/
 for(i=0 ; i<is2D ; i++)
 	{
 	if( i > 0 && side_id[i] == side_id[i-1])
 		{ point = xf2D[side_id[i]][1]; }
 	else
 		{ point = xf2D[side_id[i]][0]; }
 	switch (side_id[i])
 		{
 		case 0:
 			coord[i][0] = point;
 			coord[i][1] = -1;
 			break;
 		case 1:
 			coord[i][1] = point;
 			coord[i][0] = 1;
 			break;
 		case 2:
 			coord[i][0] = point;
 			coord[i][1] = 1;
 			break;
 		case 3:
 			coord[i][1] = point;
 			coord[i][0] = -1;
 			break;
 		} /* end of side_id[i] switch */
 	coord[i][2] = 0;
 	}
 if(is2D > 2)
 {
 for(i=0 ; i<is2D ; i++)
 	{
 	for(j=i+1 ; j<is2D ; j++)
 		{
 		dist=pow(coord[i][0]-coord[j][0],2) +
 			pow(coord[i][1]-coord[j][1],2);
 		if(dist < dist_tol)
 			{
/* 			printf("duplicate point %d\n",j);  */
 			dupl[j]=1;
 			}
 		}
 	}
 for(i=is2D-1 ; i >= 0 ; i--)
 	{
 	if(dupl[i])
 		{
 		if( i < is2D-1 && side_id[i] == side_id[i+1])
 			{ xf2D[side_id[i]][0] = xf2D[side_id[i]][1]; }
 		for (j=i ; j < is2D ; j++)
 			{ side_id[j]=side_id[j+1]; }
 		is2D--;
		nint2D[side_id[i]]--;
 		}
 	}
 
 }
return(is2D);
}
/*  end of interface_crossing_2DQ*/

/*	computation of interface inclination	*/

void
interface_inclination_2DQ(
			const double ls_F[],
			const int num_inter,
			double int_angle[],
			const int side_id[],
			double xf2D[6][2]
			)
{ 
	double gradF[9][2], angle[9], c2;
	int i, iside;
	gradF[0][0]=0.5*(-3*ls_F[0]-ls_F[1]+4*ls_F[4]);
	gradF[1][0]=0.5*(ls_F[0]+3*ls_F[1]-4*ls_F[4]);
	gradF[2][0]=0.5*(3*ls_F[2]+ls_F[3]-4*ls_F[6]);
	gradF[3][0]=0.5*(-ls_F[2]-3*ls_F[3]+4*ls_F[6]);
	gradF[4][0]=0.5*(-ls_F[0]+ls_F[1]);
	gradF[5][0]=0.5*(3*ls_F[5]+ls_F[7]-4*ls_F[8]);
	gradF[6][0]=0.5*(ls_F[2]-ls_F[3]);
	gradF[7][0]=0.5*(-ls_F[5]-3*ls_F[7]+4*ls_F[8]);
	gradF[8][0]=0.5*(ls_F[5]-ls_F[7]);
	gradF[0][1]=0.5*(-3*ls_F[0]-ls_F[3]+4*ls_F[7]);
	gradF[1][1]=0.5*(-3*ls_F[1]-ls_F[2]+4*ls_F[5]);
	gradF[2][1]=0.5*(ls_F[1]+3*ls_F[2]-4*ls_F[5]);
	gradF[3][1]=0.5*(ls_F[0]+3*ls_F[3]-4*ls_F[7]);
	gradF[4][1]=0.5*(-3*ls_F[4]-ls_F[6]+4*ls_F[8]);
	gradF[5][1]=0.5*(-ls_F[1]+ls_F[2]);
	gradF[6][1]=0.5*(ls_F[4]+3*ls_F[6]-4*ls_F[8]);
	gradF[7][1]=0.5*(-ls_F[0]+ls_F[3]);
	gradF[8][1]=0.5*(-ls_F[4]+ls_F[6]);
	for(i=0 ; i<9 ; i++)
		{
		angle[i]=(180/M_PIE)*atan2(-gradF[i][0],gradF[i][1]);
		}

/* interpolate inclination at each interface intersection	*/

	for(i=0 ; i<num_inter ; i++)
		{
		iside=side_id[i];
		c2 = xf2D[iside][0];
		if( i>0  && iside == side_id[i-1] )
			{	c2 = xf2D[iside][1];	}
	
		switch (iside)
			{
			case 0:
 				int_angle[i] = angle[0]*0.5*c2*(c2-1) + 
					angle[4]*(1-c2*c2) + angle[1]*0.5*c2*(1+c2);
				break;
			case 1:
 				int_angle[i] = angle[1]*0.5*c2*(c2-1) + 
					angle[5]*(1-c2*c2) + angle[2]*0.5*c2*(1+c2);
				break;
			case 2:
 				int_angle[i] = angle[3]*0.5*c2*(c2-1) + 
					angle[6]*(1-c2*c2) + angle[2]*0.5*c2*(1+c2);
				break;
			case 3:
 				int_angle[i] = angle[0]*0.5*c2*(c2-1) + 
					angle[7]*(1-c2*c2) + angle[3]*0.5*c2*(1+c2);
				break;
			}
		}
return;
}
#ifndef NO_CHEBYSHEV_PLEASE
/*   computation of surface determinant coefficients using
	Chebyshev polynomial representation of interface		  */

void 
surfdet_chev_coeff_2DQ( 
			const int interface_type,
			const double ls_F[],
			const double coeff[],
			double sdet_c [],
			const int n_chev,
			const int idir,
			const double end_pts[]
			)
{
double a, b, y, factor;
double fval[MAX_CHEV], sum, cder[MAX_CHEV];
double xi[DIM];
int i, j, err;
int fill_interp=0;
double dxdc, dxde, dydc, dyde, dFdc, dFde, dedc, dcde;
int use_total_deriv = 0;


a = end_pts[0];
b = end_pts[1];

xi[2]=0;

/*  find basis function for level-set field - maybe it should be mesh? */

  for (j = 0; j < Num_Basis_Functions; j++) {
    if (pd_glob[ei->mn]->i[FILL] == bfd[j]->interpolation) {
      fill_interp = j;
    }
  }


/*  compute coefficients for derivative of polynomial   */

if( !use_total_deriv )
	{
	cder[n_chev-1] = 0;
	cder[n_chev-2] = 2*(n_chev-1)*coeff[n_chev-1];
	for( j=n_chev-3 ; j>-1 ; j-- )
		{
		cder[j] = cder[j+2] + 2*(j+1)*coeff[j+1];
		}
	factor = 2./(b-a);
	for( j=0 ; j<n_chev ; j++)
		{ cder[j] *= factor;	}
	}


for (i=0 ; i<n_chev ; i++)
	{
	y = chevpoly[n_chev-3].root[i];
	sum = -0.5*coeff[0];
	for(j=0 ; j<n_chev ; j++)	
		{ sum += coeff[j]*chevpoly[n_chev-3].cosval[n_chev*j+i]; }
	if(idir)
		{
		xi[0] = sum;
		xi[1] = 0.5*(y*(b-a) + a+b);
		}
	else
		{
		xi[0] = 0.5*(y*(b-a) + a+b);
		xi[1] = sum;
		}

	err = load_basis_functions(xi, bfd);
	EH( err, "problem from load_basis_functions");
      
	err = beer_belly();
	EH( err, "beer_belly");


	dxdc = bfd[fill_interp]->J[0][0];
	dxde = bfd[fill_interp]->J[1][0];
	dydc = bfd[fill_interp]->J[0][1];
	dyde = bfd[fill_interp]->J[1][1];

/*  Should have been able to get the interface slope from the level_set
	derivatives, but using the parabolic interface representation
	coefficients instead works better
*/


	if( use_total_deriv )
		{
		dFdc = dFde =0;
		for(j=0 ; j<9 ; j++)
			{
			dFdc += ls_F[j]*bfd[fill_interp]->dphidxi[j][0];
			dFde += ls_F[j]*bfd[fill_interp]->dphidxi[j][1];
			}
		}

	if(idir)
		{
		if( use_total_deriv )
			{ dcde = -dFde/dFdc;	}
		else
			{
			dcde = -0.5*cder[0];
			for(j=0 ; j<n_chev ; j++)	
			{ dcde += cder[j]*chevpoly[n_chev-3].cosval[n_chev*j+i]; }
			}

		fval[i] = sqrt( (SQUARE(dxdc)+SQUARE(dxde))*SQUARE(dcde)
			+2*(dxdc*dxde+dydc*dyde)*dcde
			+ SQUARE(dxde)+SQUARE(dyde) );
		}
	else
		{
		if( use_total_deriv )
			{ dedc = -dFdc/dFde;	}
		else
			{
			dedc = -0.5*cder[0];
			for(j=0 ; j<n_chev ; j++)	
			{ dedc += cder[j]*chevpoly[n_chev-3].cosval[n_chev*j+i]; }
			}

		fval[i] = sqrt( SQUARE(dxdc)+SQUARE(dxde)
			+2*(dxdc*dxde+dydc*dyde)*dedc
			+ (SQUARE(dxde)+SQUARE(dyde))*SQUARE(dedc));
		}
	}
factor = 2./n_chev;
for (j=0 ; j<n_chev ; j++)
	{
	sum=0.;
	for(i=0 ; i<n_chev ; i++)
		{
		sum += fval[i]*chevpoly[n_chev-3].cosval[n_chev*j+i];
		}
	sdet_c[j] = sum*factor;
	}
return;
}
/*  end of surfdet_chev_coeff_2DQ*/


/*   computation of heaviside integration moments using
	Chebyshev polynomial representation of interface		  */

void 
heaviside_chev_moments_2DQ( 
			const int interface_type,
			const int switch_sense,
			double bf_mom[],
			const double c[],
			const int chev_order,
			const double end_pts[]
			)
{
double gauss_mom[9]={1/9.,1/9.,1/9.,1/9.,4/9.,4/9.,4/9.,4/9.,16/9.};
int npts=9;
int i;
double a,b;

a = end_pts[0];
b = end_pts[1];
switch (chev_order)
{
case 3:
switch (interface_type)
	{
	case 0:
		for(i=0;i<npts;i++)	{bf_mom[i]=gauss_mom[i];}
		break;
	case 2:
		bf_mom[0] = ((a-1)*(a-1)*(21*(-5*(1 + 2*a)*(20 + 
	(-3 + c[0])*c[0]*c[0]) + 30*a*(-2 + c[0])*c[0]*c[1] - 12*(1 + 4*a)*(-1 + c[0])
	*c[1]*c[1] + 24*a*c[1]*c[1]*c[1]) + 18*(7*(3 + 2*a)*(-2 + c[0])*c[0] + 28*a*
	(-1 + c[0])*c[1] - 4*(-1 + 8*a)*c[1]*c[1])*c[2] - 36*((19 + 30*a)*(-1 + c[0]) 
	- 22*a*c[1])*c[2]*c[2] + 8*(47 + 34*a)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[1] = -((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*(20 + (-3 + c[0])*
	c[0]*c[0]) - 30*(-2 + a + a*a)*(-2 + c[0])*c[0]*c[1] + 12*(9 + a*(7 + 4*a))*
	(-1 + c[0])*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + a*
	(11 + 2*a))*(-2 + c[0])*c[0] + 28*(-2 + a + a*a)*(-1 + c[0])*c[1] - 4*(15 + 
	a*(5 + 8*a))*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*(-1 + c[0]) - 22*
	(-2 + a + a*a)*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[2] = -((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*(-1 + c[0])*
	pow(2 + c[0],2) - 30*(-2 + a + a*a)*c[0]*(2 + c[0])*c[1] + 12*(9 + a*(7 + 
	4*a))*(1 + c[0])*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + 
	a*(11 + 2*a))*c[0]*(2 + c[0]) + 28*(-2 + a + a*a)*(1 + c[0])*c[1] - 4*(15 + 
	a*(5 + 8*a))*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*(1 + c[0]) - 22*
	(-2 + a + a*a)*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[3] = ((a-1)*(a-1)*(21*(-5*(1 + 2*a)*(-1 + c[0])*
	pow(2 + c[0],2) + 30*a*c[0]*(2 + c[0])*c[1] - 12*(1 + 4*a)*(1 + c[0])
	*c[1]*c[1] + 24*a*c[1]*c[1]*c[1]) + 18*(7*(3 + 2*a)*c[0]*(2 + c[0]) + 28*a*
	(1 + c[0])*c[1] - 4*(-1 + 8*a)*c[1]*c[1])*c[2] - 36*((19 + 30*a)*(1 + c[0]) 
	- 22*a*c[1])*c[2]*c[2] + 8*(47 + 34*a)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[4] = ((a-1)*(a-1)*(21*(5*(2 + a)*(20 + (-3 + c[0])*c[0]*c[0]) 
	- 15*(1 + a)*(-2 + c[0])*c[0]*c[1] + 12*(3 + 2*a)*(-1 + c[0])*c[1]*c[1] - 12*
	(1 + a)*c[1]*c[1]*c[1]) - 18*(7*(4 + a)*(-2 + c[0])*c[0] + 14*(1 + a)*(-1 + 
	c[0])*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] + 36*((34 + 15*a)*(-1 + c[0]) - 11*
	(1 + a)*c[1])*c[2]*c[2] - 8*(64 + 17*a)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[5] = ((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*(-4 + c[0])*
	pow(2 + c[0],2) - 30*(-2 + a + a*a)*(-4 + c[0]*c[0])*c[1] + 12*(9 + a*
	(7 + 4*a))*c[0]*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + 
	a*(11 + 2*a))*(-4 + c[0]*c[0]) + 28*(-2 + a + a*a)*c[0]*c[1] - 4*(15 + a*
	(5 + 8*a))*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*c[0] - 22*(-2 + a + a*a)
	*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))*c[2]*c[2]*c[2]))/30240.;
		bf_mom[6] = ((a-1)*(a-1)*(21*(5*(2 + a)*(-1 + c[0])*pow(2 + c[0],2) 
	- 15*(1 + a)*c[0]*(2 + c[0])*c[1] + 12*(3 + 2*a)*(1 + c[0])*c[1]*c[1] - 12*
	(1 + a)*c[1]*c[1]*c[1]) - 18*(7*(4 + a)*c[0]*(2 + c[0]) + 14*(1 + a)*
	(1 + c[0])*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] + 36*((34 + 15*a)*(1 + c[0]) 
	- 11*(1 + a)*c[1])*c[2]*c[2] - 8*(64 + 17*a)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[7] = ((a-1)*(a-1)*(21*(5*(1 + 2*a)*(-4 + c[0])*pow(2 + c[0],2)
	 - 30*a*(-4 + c[0]*c[0])*c[1] + 12*(1 + 4*a)*c[0]*c[1]*c[1] - 24*a*c[1]*c[1]
	*c[1]) - 18*(7*(3 + 2*a)*(-4 + c[0]*c[0]) + 28*a*c[0]*c[1] - 4*(-1 + 8*a)*
	c[1]*c[1])*c[2] + 36*((19 + 30*a)*c[0] - 22*a*c[1])*c[2]*c[2] - 8*(47 + 34*a)
	*c[2]*c[2]*c[2]))/30240.;
		bf_mom[8] = ((a-1)*(a-1)*(21*(-5*(2 + a)*(-4 + c[0])*pow(2 + c[0],2) 
	+ 15*(1 + a)*(-4 + c[0]*c[0])*c[1] - 12*(3 + 2*a)*c[0]*c[1]*c[1] + 12*(1 + a)
	*c[1]*c[1]*c[1]) + 18*(7*(4 + a)*(-4 + c[0]*c[0]) + 14*(1 + a)*c[0]*c[1] - 4*
	(3 + 4*a)*c[1]*c[1])*c[2] - 36*((34 + 15*a)*c[0] - 11*(1 + a)*c[1])*c[2]*c[2] 
	+ 8*(64 + 17*a)*c[2]*c[2]*c[2]))/7560.;
		break;
	case 3:
		bf_mom[0] = ((a-1)*(a-1)*(21*(5*(1 + 2*a)*pow(-2 + c[0],2)*(1 + c[0])
	 - 30*a*(-2 + c[0])*c[0]*c[1] + 12*(1 + 4*a)*(-1 + c[0])*c[1]*c[1] - 24*a*c[1]
	*c[1]*c[1]) - 18*(7*(3 + 2*a)*(-2 + c[0])*c[0] + 28*a*(-1 + c[0])*c[1] - 4*
	(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*(-1 + c[0]) - 22*a*c[1])*c[2]
	*c[2] - 8*(47 + 34*a)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[1] = ((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*pow(-2 + c[0],2)*
	(1 + c[0]) - 30*(-2 + a + a*a)*(-2 + c[0])*c[0]*c[1] + 12*(9 + a*(7 + 4*a))*
	(-1 + c[0])*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + a*(11 
	+ 2*a))*(-2 + c[0])*c[0] + 28*(-2 + a + a*a)*(-1 + c[0])*c[1] - 4*(15 + a*
	(5 + 8*a))*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*(-1 + c[0]) - 22*
	(-2 + a + a*a)*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[2] = ((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*(-20 + c[0]*c[0]*
	(3 + c[0])) - 30*(-2 + a + a*a)*c[0]*(2 + c[0])*c[1] + 12*(9 + a*(7 + 4*a))*
	(1 + c[0])*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + a*
	(11 + 2*a))*c[0]*(2 + c[0]) + 28*(-2 + a + a*a)*(1 + c[0])*c[1] - 4*(15 + a*
	(5 + 8*a))*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*(1 + c[0]) - 22*
	(-2 + a + a*a)*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[3] = ((a-1)*(a-1)*(21*(5*(1 + 2*a)*(-20 + c[0]*c[0]*(3 + c[0]))
	 - 30*a*c[0]*(2 + c[0])*c[1] + 12*(1 + 4*a)*(1 + c[0])*c[1]*c[1] - 24*a*c[1]
	*c[1]*c[1]) - 18*(7*(3 + 2*a)*c[0]*(2 + c[0]) + 28*a*(1 + c[0])*c[1] - 4*
	(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*(1 + c[0]) - 22*a*c[1])*c[2]*c[2]
	 - 8*(47 + 34*a)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[4] = ((a-1)*(a-1)*(21*(-5*(2 + a)*pow(-2 + c[0],2)*(1 + c[0]) 
	+ 15*(1 + a)*(-2 + c[0])*c[0]*c[1] - 12*(3 + 2*a)*(-1 + c[0])*c[1]*c[1] + 
	12*(1 + a)*c[1]*c[1]*c[1]) + 18*(7*(4 + a)*(-2 + c[0])*c[0] + 14*(1 + a)*
	(-1 + c[0])*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] - 36*((34 + 15*a)*(-1 + c[0]) 
	- 11*(1 + a)*c[1])*c[2]*c[2] + 8*(64 + 17*a)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[5] = -((-1 + a)*(21*(5*(5 + a*(5 + 2*a))*pow(-2 + c[0],2)*
	(4 + c[0]) - 30*(-2 + a + a*a)*(-4 + c[0]*c[0])*c[1] + 12*(9 + a*(7 + 4*a))
	*c[0]*c[1]*c[1] - 24*(-2 + a + a*a)*c[1]*c[1]*c[1]) - 18*(7*(7 + a*(11 + 2*a))
	*(-4 + c[0]*c[0]) + 28*(-2 + a + a*a)*c[0]*c[1] - 4*(15 + a*(5 + 8*a))
	*c[1]*c[1])*c[2] + 36*((79 + 87*a + 30*a*a)*c[0] - 22*(-2 + a + a*a)
	*c[1])*c[2]*c[2] - 8*(115 + a*(175 + 34*a))*c[2]*c[2]*c[2]))/30240.;
		bf_mom[6] = ((a-1)*(a-1)*(21*(-5*(2 + a)*(-20 + c[0]*c[0]*(3 + c[0]))
	 + 15*(1 + a)*c[0]*(2 + c[0])*c[1] - 12*(3 + 2*a)*(1 + c[0])*c[1]*c[1] + 12*
	(1 + a)*c[1]*c[1]*c[1]) + 18*(7*(4 + a)*c[0]*(2 + c[0]) + 14*(1 + a)*
	(1 + c[0])*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] - 36*((34 + 15*a)*(1 + c[0]) 
	- 11*(1 + a)*c[1])*c[2]*c[2] + 8*(64 + 17*a)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[7] = ((a-1)*(a-1)*(21*(-5*(1 + 2*a)*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*a*(-4 + c[0]*c[0])*c[1] - 12*(1 + 4*a)*c[0]*c[1]*c[1] + 
	24*a*c[1]*c[1]*c[1]) + 18*(7*(3 + 2*a)*(-4 + c[0]*c[0]) + 28*a*c[0]*c[1] - 
	4*(-1 + 8*a)*c[1]*c[1])*c[2] - 36*((19 + 30*a)*c[0] - 22*a*c[1])*c[2]*c[2] + 
	8*(47 + 34*a)*c[2]*c[2]*c[2]))/30240.;
		bf_mom[8] = ((a-1)*(a-1)*(21*(5*(2 + a)*pow(-2 + c[0],2)*(4 + c[0]) 
	- 15*(1 + a)*(-4 + c[0]*c[0])*c[1] + 12*(3 + 2*a)*c[0]*c[1]*c[1] - 12*(1 + a)
	*c[1]*c[1]*c[1]) - 18*(7*(4 + a)*(-4 + c[0]*c[0]) + 14*(1 + a)*c[0]*c[1] - 
	4*(3 + 4*a)*c[1]*c[1])*c[2] + 36*((34 + 15*a)*c[0] - 11*(1 + a)*c[1])
	*c[2]*c[2] - 8*(64 + 17*a)*c[2]*c[2]*c[2]))/7560.;
		break;
	case 4:
		bf_mom[0] = -((1 + b)*(21*(5*(5 + b*(-5 + 2*b))*pow(-2 + c[0],2)*
	(1 + c[0]) + 30*(-2 + b)*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(9 + b*(-7 + 
	4*b))*(-1 + c[0])*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*
	(7 + b*(-11 + 2*b))*(-2 + c[0])*c[0] - 28*(-2 + b)*(1 + b)*(-1 + c[0])*c[1] 
	- 4*(15 + b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*(-1 + 
	c[0]) + 22*(-2 + b)*(1 + b)*c[1])*c[2]*c[2] - 8*(115 + b*(-175 + 34*b))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[1] = ((1+b)*(1+b)*(21*(-5*(-1 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) - 30*b*(-2 + c[0])*c[0]*c[1] - 12*(-1 + 4*b)*(-1 + c[0])
	*c[1]*c[1] - 24*b*c[1]*c[1]*c[1]) + 18*(7*(-3 + 2*b)*(-2 + c[0])*c[0] - 
	28*b*(-1 + c[0])*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] - 36*((-19 + 30*b)*
	(-1 + c[0]) + 22*b*c[1])*c[2]*c[2] + 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[2] = ((1+b)*(1+b)*(21*(-5*(-1 + 2*b)*(-20 + c[0]*c[0]*
	(3 + c[0])) - 30*b*c[0]*(2 + c[0])*c[1] - 12*(-1 + 4*b)*(1 + c[0])*c[1]*c[1]
	 - 24*b*c[1]*c[1]*c[1]) + 18*(7*(-3 + 2*b)*c[0]*(2 + c[0]) - 28*b*(1 + c[0])
	*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] - 36*((-19 + 30*b)*(1 + c[0]) + 
	22*b*c[1])*c[2]*c[2] + 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[3] = -((1 + b)*(21*(5*(5 + b*(-5 + 2*b))*(-20 + c[0]*c[0]*
	(3 + c[0])) + 30*(-2 + b)*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(9 + b*
	(-7 + 4*b))*(1 + c[0])*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 
	18*(7*(7 + b*(-11 + 2*b))*c[0]*(2 + c[0]) - 28*(-2 + b)*(1 + b)*(1 + c[0])
	*c[1] - 4*(15 + b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*
	(1 + c[0]) + 22*(-2 + b)*(1 + b)*c[1])*c[2]*c[2] - 8*(115 + b*(-175 + 34*b))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[4] = ((1+b)*(1+b)*(21*(5*(-2 + b)*pow(-2 + c[0],2)*(1 + c[0])
	 + 15*(-1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(-3 + 2*b)*(-1 + c[0])*c[1]*c[1] 
	+ 12*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-4 + b)*(-2 + c[0])*c[0] - 14*(-1 + b)
	*(-1 + c[0])*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] + 36*((-34 + 15*b)*
	(-1 + c[0]) + 11*(-1 + b)*c[1])*c[2]*c[2] - 8*(-64 + 17*b)
	*c[2]*c[2]*c[2]))/15120.;
		bf_mom[5] = ((1+b)*(1+b)*(21*(5*(-1 + 2*b)*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*b*(-4 + c[0]*c[0])*c[1] + 12*(-1 + 4*b)*c[0]*c[1]*c[1] + 
	24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-4 + c[0]*c[0]) - 28*b*c[0]*c[1] -
	 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*c[0] + 22*b*c[1])*c[2]*c[2] 
	- 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/30240.;
		bf_mom[6] = ((1+b)*(1+b)*(21*(5*(-2 + b)*(-20 + c[0]*c[0]*(3 + c[0]))
	 + 15*(-1 + b)*c[0]*(2 + c[0])*c[1] + 12*(-3 + 2*b)*(1 + c[0])*c[1]*c[1] + 
	12*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-4 + b)*c[0]*(2 + c[0]) - 14*(-1 + b)*
	(1 + c[0])*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] + 36*((-34 + 15*b)*(1 + c[0])
	 + 11*(-1 + b)*c[1])*c[2]*c[2] - 8*(-64 + 17*b)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[7] = ((1 + b)*(21*(5*(5 + b*(-5 + 2*b))*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*(-2 + b)*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(9 + b*
	(-7 + 4*b))*c[0]*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*
	(7 + b*(-11 + 2*b))*(-4 + c[0]*c[0]) - 28*(-2 + b)*(1 + b)*c[0]*c[1] - 4*
	(15 + b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*c[0] + 22*
	(-2 + b)*(1 + b)*c[1])* c[2]*c[2] - 8*(115 + b*(-175 + 34*b))
	*c[2]*c[2]*c[2]))/30240.;
		bf_mom[8] = ((1+b)*(1+b)*(21*(-5*(-2 + b)*pow(-2 + c[0],2)*
	(4 + c[0]) - 15*(-1 + b)*(-4 + c[0]*c[0])*c[1] - 12*(-3 + 2*b)*c[0]*c[1]
	*c[1] - 12*(-1 + b)*c[1]*c[1]*c[1]) + 18*(7*(-4 + b)*(-4 + c[0]*c[0]) - 
	14*(-1 + b)*c[0]*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] - 36*((-34 + 15*b)*c[0]
	 + 11*(-1 + b)*c[1])*c[2]*c[2] + 8*(-64 + 17*b)*c[2]*c[2]*c[2]))/7560.;
		break;
	case 6:
		bf_mom[0] = ((1 + b)*(-525*(-20 + 3*c[0]*c[0]) + 21*(5*(20*b*
	(-5 + 2*b) - 3*b*(-5 + 2*b)*c[0]*c[0] + (5 + b*(-5 + 2*b))*c[0]*c[0]*c[0])
	 + 30*(-2 + b)*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(9 + b*(-7 + 4*b))*
	(-1 + c[0])*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(7 + 
	b*(-11 + 2*b))*(-2 + c[0])*c[0] - 28*(-2 + b)*(1 + b)*(-1 + c[0])*c[1] - 
	4*(15 + b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*
	(-1 + c[0]) + 22*(-2 + b)*(1 + b)*c[1])*c[2]*c[2] - 8*(115 + b*
	(-175 + 34*b))*c[2]*c[2]*c[2]))/60480.;
		bf_mom[1] = (pow(1 + b,2)*(21*(5*(-1 + 2*b)*(20+(-3+c[0])*c[0]*c[0])
	 + 30*b*(-2 + c[0])*c[0]*c[1] + 12*(-1 + 4*b)*(-1 + c[0])*c[1]*c[1] + 
	24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-2 + c[0])*c[0] - 28*b*(-1 + c[0])
	*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*(-1 + c[0]) + 
	22*b*c[1])*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[2] = (pow(1 + b,2)*(21*(5*(-1 + 2*b)*(-1 + c[0])*
	pow(2 + c[0],2) + 30*b*c[0]*(2 + c[0])*c[1] + 12*(-1 + 4*b)*(1 + c[0])
	*c[1]*c[1] + 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*c[0]*(2 + c[0]) - 
	28*b*(1 + c[0])*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*
	(1 + c[0]) + 22*b*c[1])*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/60480.;
		bf_mom[3] = ((1 + b)*(525*(-4 + 3*c[0]*c[0]) + 21*(5*(5*c[0]*c[0]
	*c[0] - 5*b*(-1 + c[0])*pow(2 + c[0],2) + 2*b*b*(-1 + c[0])*pow(2 + c[0],2))
	 + 30*(-2 + b)*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(9 + b*(-7 + 4*b))*
	(1 + c[0])*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(7 + b*
	(-11 + 2*b))*c[0]*(2 + c[0]) - 28*(-2 + b)*(1 + b)*(1 + c[0])*c[1] - 4*(15 
	+ b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*(1 + c[0]) + 22*
	(-2 + b)*(1 + b)*c[1])*c[2]*c[2] - 8*(115 + b*(-175 + 34*b))
	*c[2]*c[2]*c[2]))/60480.;
		bf_mom[4] = (pow(1 + b,2)*(21*(-5*(-2 + b)*(20 + (-3 + c[0])*c[0]
	*c[0]) - 15*(-1 + b)*(-2 + c[0])*c[0]*c[1] - 12*(-3 + 2*b)*(-1 + c[0])*c[1]
	*c[1] - 12*(-1 + b)*c[1]*c[1]*c[1]) + 18*(7*(-4 + b)*(-2 + c[0])*c[0] - 14*
	(-1 + b)*(-1 + c[0])*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] - 36*((-34 + 15*b)
	*(-1 + c[0]) + 11*(-1 + b)*c[1])*c[2]*c[2] + 8*(-64 + 17*b)
	*c[2]*c[2]*c[2]))/15120.;
		bf_mom[5] = (pow(1 + b,2)*(21*(-5*(-1 + 2*b)*(-4 + c[0])*
	pow(2 + c[0],2) - 30*b*(-4 + c[0]*c[0])*c[1] - 12*(-1 + 4*b)*c[0]*c[1]*c[1]
	 - 24*b*c[1]*c[1]*c[1]) + 18*(7*(-3 + 2*b)*(-4 + c[0]*c[0]) - 28*b*c[0]*c[1]
	 - 4*(1 + 8*b)*c[1]*c[1])*c[2] - 36*((-19 + 30*b)*c[0] + 22*b*c[1])
	*c[2]*c[2] + 8*(-47 + 34*b)*c[2]*c[2]*c[2]))/30240.;
		bf_mom[6] = (pow(1 + b,2)*(21*(-5*(-2 + b)*(-1 + c[0])*
	pow(2 + c[0],2) - 15*(-1 + b)*c[0]*(2 + c[0])*c[1] - 12*(-3 + 2*b)*
	(1 + c[0])*c[1]*c[1] - 12*(-1 + b)*c[1]*c[1]*c[1]) + 18*(7*(-4 + b)*c[0]*
	(2 + c[0]) - 14*(-1 + b)*(1 + c[0])*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] - 
	36*((-34 + 15*b)*(1 + c[0]) + 11*(-1 + b)*c[1])*c[2]*c[2] + 8*(-64 + 17*b)
	*c[2]*c[2]*c[2]))/15120.;
		bf_mom[7] = -((1 + b)*(21*(5*(5 + b*(-5 + 2*b))*(-4 + c[0])*
	pow(2 + c[0],2) + 30*(-2 + b)*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(9 + b*
	(-7 + 4*b))*c[0]*c[1]*c[1] + 24*(-2 + b)*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*
	(7 + b*(-11 + 2*b))*(-4 + c[0]*c[0]) - 28*(-2 + b)*(1 + b)*c[0]*c[1] - 4*
	(15 + b*(-5 + 8*b))*c[1]*c[1])*c[2] + 36*((79 - 87*b + 30*b*b)*c[0] + 22*
	(-2 + b)*(1 + b)*c[1])* c[2]*c[2] - 8*(115 + b*(-175 + 34*b))*
	c[2]*c[2]*c[2]))/30240.;
		bf_mom[8] = (pow(1 + b,2)*(21*(5*(-2 + b)*(-4 + c[0])*pow(2 + c[0],2)
	 + 15*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-3 + 2*b)*c[0]*c[1]*c[1] + 12*
	(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-4 + b)*(-4 + c[0]*c[0]) - 14*(-1 + b)*
	c[0]*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] + 36*((-34 + 15*b)*c[0] + 11*
	(-1 + b)*c[1])*c[2]*c[2] - 8*(-64 + 17*b)*c[2]*c[2]*c[2]))/7560.;
		break;
	case 22:
		bf_mom[0] = (21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 
	30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) 
	- 18*(7*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
    	396*(-1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[1] = (21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 
	30*c[0]*(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 
    	18*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
    	396*(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[2] = (21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 
	30*c[0]*(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 
    	18*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
    	396*(1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[3] = (21*(-5*pow(-2 + c[0],2)*(1 + c[0]) 
	- 30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) 
	- 18*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
    	396*(-1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[4] = (21*(5*pow(-2 + c[0],2)*(4 + c[0]) 
	- 30*(-4 + c[0]*c[0])*c[1] + 36*c[0]*c[1]*c[1] - 24*pow(c[1],3)) + 
    	18*(-28 + 7*c[0]*c[0] - 28*c[0]*c[1] + 36*c[1]*c[1])*c[2] + 
    	396*(c[0] - 2*c[1])*c[2]*c[2] + 104*pow(c[2],3))/7560.;
		bf_mom[5] = (-105*(-20 + c[0]*c[0]*(3 + c[0])) - 
	252*(1 + c[0])*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 
    	684*(1 + c[0])*c[2]*c[2] + 376*pow(c[2],3))/3780.;
		bf_mom[6] = (105*pow(c[0],3) + 126*c[0]*c[0]*(5*c[1] + c[2]) + 
    	36*c[0]*(-35 + 21*c[1]*c[1] + 14*c[1]*c[2] + 11*c[2]*c[2]) + 
    	8*(210 + 63*pow(c[1],3) - 63*c[2] + 81*c[1]*c[1]*c[2] + 
       	13*pow(c[2],3) + 9*c[1]*(-35 + 11*c[2]*c[2])))/7560.;
		bf_mom[7] = (-105*pow(-2 + c[0],2)*(1 + c[0]) - 
	252*(-1 + c[0])*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] 
	- 684*(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3))/3780.;
		bf_mom[8] = (105*pow(c[0],3) - 378*c[0]*c[0]*c[2] + 
    	36*c[0]*(-35 + 7*c[1]*c[1] + 19*c[2]*c[2]) - 8*(-210 + 9*(-21 + 
	c[1]*c[1])*c[2] + 47*pow(c[2],3)))/1890.;
		break;
	case 32:
		bf_mom[0] = (21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 30*(-2 + c[0])
	*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 18*(7*(-2 + 
	c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 396*(-1 + c[0] 
	- 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[1] = (21*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 30*(-2 + c[0])
	*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*(7*(-2 + 
	c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 396*(-1 + c[0] 
	+ 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[2] = (21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 30*c[0]*(2 + 
	c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*(7*c[0]*
	(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 396*(1 + c[0] 
	+ 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[3] = (21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 396*
	(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3))/15120.;
		bf_mom[4] = (-105*pow(-2 + c[0],2)*(1 + c[0]) - 252*(-1 + c[0])
	*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] - 684*
	(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3))/3780.;
		bf_mom[5] = (105*pow(c[0],3) + 126*c[0]*c[0]*(5*c[1] + c[2]) 
	+ 36*c[0]*(-35 + 21*c[1]*c[1] + 14*c[1]*c[2] + 11*c[2]*c[2]) + 8*
	(210 + 63*pow(c[1],3) - 63*c[2] + 81*c[1]*c[1]*c[2] + 13*pow(c[2],3) 
	+ 9*c[1]*(-35 + 11*c[2]*c[2])))/7560.;
		bf_mom[6] = (-105*(-20 + c[0]*c[0]*(3 + c[0])) - 252*(1 + c[0])
	*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 684*(1 + c[0])
	*c[2]*c[2] + 376*pow(c[2],3))/3780.;
		bf_mom[7] = (21*(5*pow(-2 + c[0],2)*(4 + c[0]) - 30*(-4 + c[0]
	*c[0])*c[1] + 36*c[0]*c[1]*c[1] - 24*pow(c[1],3)) + 18*(-28 + 7*c[0]
	*c[0] - 28*c[0]*c[1] + 36*c[1]*c[1])*c[2] + 396*(c[0] - 2*c[1])
	*c[2]*c[2] + 104*pow(c[2],3))/7560.;
		bf_mom[8] = (105*pow(c[0],3) - 378*c[0]*c[0]*c[2] + 36*c[0]*
	(-35 + 7*c[1]*c[1] + 19*c[2]*c[2]) - 8*(-210 + 9*(-21 + c[1]*c[1])
	*c[2] + 47*pow(c[2],3)))/1890.;
		break;
	case 100:
		bf_mom[0] = -((a - b)*(b*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])
	*c[0]*c[0]) + 30*(-1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 4*b)*(-1 + c[0])
	*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-2 + c[0])*c[0]
	 - 28*(-1 + b)*(-1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])* c[2] + 36*((-49 +
	 30*b)*(-1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*
	c[2]) + a*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])*c[0]*c[0]) + 30*(-2 + c[0])*
	c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*
	(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 4*(7 + 2*b)*c[1]*c[1])*
	c[2] + 36*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)*
	c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*(5 + 5*c[1] + 
	2*c[2]) + 18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*c[2])
	 + 4*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*
	c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[1] = -((a - b)*(b*(21*(5*(3 + 2*b)*(20 + (-3 + c[0])*c[0]*c[0])
	 + 30*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(5 + 4*b)*(-1 + c[0])*c[1]*c[1] + 
	24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-2 + c[0])*c[0] - 28*(1 + b)*
	(-1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*(-1 + c[0]) 
	+ 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*(3 +
	 2*b)*(20 + (-3 + c[0])*c[0]*c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 2*b)
	*(-1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2 + c[0])*c[0]
	 + 28*(-1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*
	(-1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*
	(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*(5 + 5*c[1] + 2*c[2]) + 18*c[0]*(7*c[1]*
	(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*c[2]) + 4*(525 - 63*c[1]*c[1]*
	(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[2] = -((a - b)*(2*a*a*(21*(5*(-1 + c[0])*(2+c[0])*(2+c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 + 
	15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(3 + 2*b)*(-1 
	+ c[0])*(2+c[0])*(2+c[0]) + 30*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(5 + 4*b)*
	(1 + c[0])*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*c[0]*
	(2 + c[0]) - 28*(1 + b)*(1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*
	((49 + 30*b)*(1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*
	c[2]*c[2]) + a*(21*(5*(3 + 2*b)*(-1 + c[0])*(2+c[0])*(2+c[0]) - 30*c[0]*(2 
	+ c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*
	(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*
	c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*
	c[2]*c[2])))/60480.;
		bf_mom[3] = -((a - b)*(2*a*a*(21*(5*(-1 + c[0])*(2+c[0])*(2+c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 + 15*
	c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)*(-1 + 
	c[0])*(2+c[0])*(2+c[0]) + 30*(-1 + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 + 4*b)*
	(1 + c[0])*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*c[0]*
	(2 + c[0]) - 28*(-1 + b)*(1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])* c[2] + 
	36*((-49 + 30*b)*(1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*(-81 + 34*b)*
	c[2]*c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(-1 + c[0])*(2+c[0])*(2+c[0]) + 30*c[0]*
	(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 
	18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*(7 + 2*b)*c[1]*
	c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*c[2]*c[2] 
	- 8*(-81+94*b)* c[2]*c[2]*c[2])))/60480.;
		bf_mom[4] = (63*(a - b)*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 20*
	c[0]*c[2] + 28*c[2]*c[2]) + 9*(-35*(a*a*a - b*b*b)* (-4 + c[0]*c[0]) + 
	70*(a-b)*(a-b)*(a + b)*c[0]*c[1] - 28*(a - b)*(2*a*a + a*b + 2*b*b)*c[1]*c[1]
	 + 28*(a - b)*((a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*c[2] - 4*
	(a - b)*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2]) + (a - b)*(21*(5*(a*a + a*b + 
	b*b)*(8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]*c[0]*c[1] + 12*(2*a*a + 
	a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) - 18*(7*
	(a*a + 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)*c[0]*c[1] - 4*(4*a*a - 
	a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b + 15*b*b)*c[0] - 11*
	(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b +17*b*b)*c[2]*c[2]*c[2])
	 - 9*(a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*c[2] + 28*c[0]*(5*c[1]*c[1] +
	 7*c[2]*c[2]) + 8*(35 + 7*c[1]*c[1]*c[2] - 9*c[2]*c[2]*c[2])))/15120.;
		bf_mom[5] = ((a - b)*(b*(21*(5*(3 + 2*b)*(-4 + c[0])*(2+c[0])*(2+c[0])
	 + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*c[1]*c[1] + 24*
	(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) - 28*(1 + b)*
	c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0] + 22*(1 + b)*
	c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*(3 + 2*b)*(-4 + 
	c[0])*(2+c[0])*(2+c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*(5 + 2*b)*c[0]*
	c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*
	c[0]*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0] - 22*c[1])*
	c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 84*
	(20 + 3*c[1]*(-5 + c[1]*c[1])) + 72*(7 + 4*c[1]*c[1])*c[2] - 396*c[1]*c[2]*
	c[2] - 136*c[2]*c[2]*c[2] - 63*c[0]*c[0]*(5*c[1] + 2*c[2]) + 36*c[0]*
	(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 15*c[2]*c[2]))))/30240.;
		bf_mom[6] = ((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*(-1 + c[0])* 
	(2+c[0])*(2+c[0]) - 15*(a - b)*(a + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*a*a
	 + a*b + 2*b*b)*(1 + c[0])* c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) - 
	18*(7*(-5 + a*a + 3*a*b + b*b)*c[0]*(2 + c[0]) + 14*(a - b)*(a + b)*(1 + c[0])
	*c[1] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((-49 + 15*a*a + 
	19*a*b + 15*b*b)*(1 + c[0]) - 11*(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(-81 + 
	17*a*a + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[7] = ((a - b)*(b*(21*(5*(-3 + 2*b)*(-4 + c[0])*(2+c[0])*
	(2+c[0]) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]*c[1]*c[1] 
	+ 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0]) - 28*
	(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)*c[0] + 
	22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*
	(-3 + 2*b)*(-4 + c[0])*(2+c[0])*(2+c[0]) + 30*(-4 + c[0]*c[0])*c[1] + 12*
	(-5 + 2*b)*c[0]*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + c[0]*
	c[0]) - 28*c[0]*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*c[0] +
	 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*
	c[0] - 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 72*(7 + 4*c[1]*c[1])*c[2] - 
	396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] - 63*c[0]*c[0]*(5*c[1] + 2*c[2]) + 
	36*c[0]*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 15*c[2]*c[2]))))/30240.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + 
	((a - b)*(5*(a*a + a*b + b*b)*(2 + c[0]) - 5*(a - b)*(a + b)*c[1] - 2*
	(a*a + 3*a*b + b*b)*c[2]))/ 30. - ((a - b)*(21*(5*(a*a + a*b + b*b)* 
	(8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]*c[0]*c[1] + 12*(2*a*a + a*b
	 + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) - 18*(7*(a*a 
	+ 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)*c[0]*c[1] - 4*(4*a*a - a*b + 
	4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b + 15*b*b)*c[0] - 11*(a - b)*
	(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/7560.
	 + ((a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*c[2] + 28*c[0]*(5*c[1]*c[1] + 
	7*c[2]*c[2]) + 8*(35 + 7*c[1]*c[1]*c[2] - 9*c[2]*c[2]*c[2])))/840.;
		break;
	case 101:
		bf_mom[0] = ((a - b)*(b*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*(1+c[0])
	 + 30*(-1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 4*b)*(-1 + c[0])*c[1]*c[1] 
	+ 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-2 + c[0])*c[0] - 28*
	(-1 + b)*(-1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])* c[2] + 36*((-49+30*b)
	*(-1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2])
	 + a*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*(1 + c[0]) + 30*(-2 + c[0])*c[0]
	*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*
	(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 4*(7 + 2*b)*c[1]*c[1])
	*c[2] + 36*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)
	*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*(5 + 5*c[1] + 
	2*c[2]) + 18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*c[2])
	 + 4*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 
	11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[1] = ((a - b)*(2*a*a*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)
	*(-20 + c[0]*c[0]*(3 + c[0])) + 30*(-1 + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 
	+ 4*b)*(1 + c[0])*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5+2*b)
	*c[0]*(2 + c[0]) - 28*(-1 + b)*(1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])
	*c[2] + 36*((-49 + 30*b)*(1 + c[0]) + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 
	+ 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(-20 + c[0]*c[0]*(3 + c[0])) 
	+ 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] + 24*c[1]
	*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*
	(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*c[2]*c[2]
	 - 8*(-81 + 94*b)*c[2]*c[2]*c[2])))/60480.;
		bf_mom[2] = ((a - b)*(2*a*a*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b* (21*(5*(3 + 2*b)
	*(-20 + c[0]*c[0]*(3 + c[0])) + 30*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(5 + 
	4*b)*(1 + c[0])*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)
	*c[0]*(2 + c[0]) - 28*(1 + b)*(1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2]
	 + 36*((49 + 30*b)*(1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)
	*c[2]*c[2]*c[2]) + a*(21*(5*(3 + 2*b)*(-20 + c[0]*c[0]*(3 + c[0])) - 30*c[0]
	*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) 
	- 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]
	*c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81+94*b)
	*c[2]*c[2]*c[2])))/60480.;
		bf_mom[3] = ((a - b)*(b*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*(1 + c[0])
	 + 30*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(5 + 4*b)*(-1 + c[0])*c[1]*c[1] + 
	24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-2 + c[0])*c[0] - 28*(1 + b)*
	(-1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*(-1 + c[0])
	 + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*
	(3 + 2*b)*pow(-2 + c[0],2)*(1 + c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 
	2*b)*(-1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2+c[0])
	*c[0] + 28*(-1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49+38*b)
	*(-1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*
	(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*(5 + 5*c[1] + 2*c[2]) + 18*c[0]*(7*c[1]*
	(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*c[2]) + 4*(105 - 63*c[1]*c[1]*
	(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[4] = -((a - b)*(2*a*a*(21*(5*pow(-2 + c[0],2)*(4 + c[0]) - 
	15*(-4 + c[0]*c[0])*c[1] + 24*c[0]*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*c[0] + 14*c[0]*c[1] - 4*(7 + 4*c[1]*c[1]))* c[2] + 36*(15*c[0] - 
	11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)*
	pow(-2 + c[0],2)*(4 + c[0]) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 
	+ 4*b)*c[0]*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*
	(-4 + c[0]*c[0]) - 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 
	36*((-49 + 30*b)*c[0] + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]
	*c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*(4 + c[0]) + 30*(-4 + 
	c[0]*c[0])*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*
	(7*(-5 + 6*b)*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2]
	 + 36*((-49 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 
	8*(-81 + 94*b)*c[2]*c[2]*c[2])))/30240.;
		bf_mom[5] = (-6300*a + 6300*b + 63*(a - b)* (15*c[0]*c[0] + 20*c[1]
	*c[1] - 20*c[0]*c[2] + 28*c[2]*c[2]) + 9*(-35*(a*a*a - b*b*b)* (-4 + c[0]
	*c[0]) + 70*pow(a - b,2)*(a + b)*c[0]*c[1] - 28*(a - b)*(2*a*a + a*b+2*b*b)
	*c[1]*c[1] + 28*(a - b)*((a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])
	*c[2] - 4*(a - b)*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2]) - (a - b)*(21*(5*
	(a*a + a*b + b*b)*(-8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]*c[0]*c[1]
	 + 12*(2*a*a + a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]
	*c[1]) - 18*(7*(a*a + 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)*c[0]*c[1]
	 - 4*(4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b + 15*b*b)
	*c[0] - 11*(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b + 17*b*b)
	*c[2]*c[2]*c[2]) + 9*(a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*c[2] + 
	8*c[2]*(7*c[1]*c[1] - 9*c[2]*c[2]) + 28*c[0]*(5*c[1]*c[1] + 
	7*c[2]*c[2])))/15120.;
		bf_mom[6] = -((a - b)*(2*a*a*(21*(5*pow(-2 + c[0],2)*(4 + c[0]) - 
	15*(-4 + c[0]*c[0])*c[1] + 24*c[0]*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*c[0] + 14*c[0]*c[1] - 4*(7 + 4*c[1]*c[1]))* c[2] + 36*(15*c[0] - 
	11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(3 + 2*b)*
	pow(-2 + c[0],2)*(4 + c[0]) + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5+4*b)
	*c[0]*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]
	*c[0]) - 28*(1 + b)*c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49+30*b)
	*c[0] + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*
	(5*(3 + 2*b)*pow(-2 + c[0],2)*(4 + c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*
	(5 + 2*b)*c[0]*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]
	*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0]
	 - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2])))/30240.;
		bf_mom[7] = -((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*pow(-2+c[0],2)
	* (1 + c[0]) - 15*(a - b)*(a + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 2*a*a +
	 a*b + 2*b*b)*(-1 + c[0])* c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) 
	- 18*(7*(-5 + a*a + 3*a*b + b*b)*(-2 + c[0])*c[0] + 14*(a - b)*(a + b)*
	(-1 + c[0])*c[1] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*(
	(-49 + 15*a*a + 19*a*b + 15*b*b)*(-1 + c[0]) - 11*(a - b)*(a + b)*c[1])
	*c[2]*c[2] - 8*(-81 + 17*a*a + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/ 15120.;
		bf_mom[8] = ((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*pow(-2 + c[0],2)
	* (4 + c[0]) - 15*(a - b)*(a + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 2*a*a 
	+ a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) - 18*
	(7*(-5 + a*a + 3*a*b + b*b)*(-4 + c[0]*c[0]) + 14*(a - b)*(a + b)*c[0]*c[1]
	 - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((-49 + 15*a*a + 
	19*a*b + 15*b*b)*c[0] + 11*(-a*a + b*b)*c[1])*c[2]*c[2] - 8*(-81 + 17*a*a 
	+ 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/7560.;
		break;
	case 102:
		bf_mom[0] = ((a - b)*(b*(21*(5*(-3 + 2*b)*(c[0]-2)*(c[0]-2)*
	(1 + c[0]) + 30*(-1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 4*b)*(-1 + c[0])
	*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-2 + c[0])*c[0]
	 - 28*(-1 + b)*(-1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])* c[2] + 36*((-49 
	+ 30*b)*(-1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*(-81 + 34*b)*c[2]*
	c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(c[0]-2)*(c[0]-2)*(1 + c[0]) + 30*(-2 + 
	c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1] + 24*c[1]*c[1]*c[1]) 
	- 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 4*(7 + 2*b)*
	c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2]*c[2] - 8*
	(-81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*
	(5 + 5*c[1] + 2*c[2]) + 18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] 
	+ 30*c[2]*c[2]) + 4*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])
	*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[1] = ((a - b)*(b*(21*(5*(3 + 2*b)*(c[0]-2)*(c[0]-2)*
	(1 + c[0]) + 30*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(5 + 4*b)*(-1 + c[0])
	*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-2 + c[0])*c[0]
	 - 28*(1 + b)*(-1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 
	30*b)*(-1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*
	c[2]) + a*(21*(5*(3 + 2*b)*(c[0]-2)*(c[0]-2)*(1 + c[0]) - 30*(-2 + c[0])*
	c[0]*c[1] + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*
	(7*(5 + 6*b)*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*
	c[1])*c[2] + 36*((49 + 38*b)*(-1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*
	(81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*
	(5 + 5*c[1] + 2*c[2]) + 18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] 
	+ 30*c[2]*c[2]) + 4*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])
	*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[2] = ((a - b)*(2*a*a*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(3 + 2*b)*
	(-20 + c[0]*c[0]*(3 + c[0])) + 30*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(5 + 
	4*b)*(1 + c[0])*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)
	*c[0]*(2 + c[0]) - 28*(1 + b)*(1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2]
	 + 36*((49 + 30*b)*(1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)
	*c[2]*c[2]*c[2]) + a*(21*(5*(3 + 2*b)*(-20 + c[0]*c[0]*(3 + c[0])) - 30*
	c[0]*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*
	c[1]) - 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 + 2*b)
	*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*
	(81 + 94*b)*c[2]*c[2]*c[2])))/60480.;
		bf_mom[3] = ((a - b)*(2*a*a*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)
	*(-20 + c[0]*c[0]*(3 + c[0])) + 30*(-1 + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 
	+ 4*b)*(1 + c[0])*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)
	*c[0]*(2 + c[0]) - 28*(-1 + b)*(1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])
	*c[2] + 36*((-49 + 30*b)*(1 + c[0]) + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*
	(-81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(-20 + c[0]*c[0]*(3 + 
	c[0])) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] + 
	24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1]
	 + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*
	c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2])))/60480.;
		bf_mom[4] = -((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*(c[0]-2)*
	(c[0]-2)* (1 + c[0]) - 15*(a - b)*(a + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 
	+ 2*a*a + a*b + 2*b*b)*(-1 + c[0])* c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*
	c[1]*c[1]) - 18*(7*(-5 + a*a + 3*a*b + b*b)*(-2 + c[0])*c[0] + 14*(a - b)*
	(a + b)*(-1 + c[0])*c[1] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 
	36*((-49 + 15*a*a + 19*a*b + 15*b*b)*(-1 + c[0]) - 11*(a - b)*(a + b)*c[1])
	*c[2]*c[2] - 8*(-81 + 17*a*a + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/ 15120.;
		bf_mom[5] = -((a - b)*(2*a*a*(21*(5*(c[0]-2)*(c[0]-2)*(4 + c[0]) - 
	15*(-4 + c[0]*c[0])*c[1] + 24*c[0]*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*c[0] + 14*c[0]*c[1] - 4*(7 + 4*c[1]*c[1]))*c[2] + 36*(15*c[0] - 
	11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(3 + 2*b)*(c[0]-2)*
	(c[0]-2)*(4 + c[0]) + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]
	*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) 
	- 28*(1 + b)*c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0]
	 + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*
	(3 + 2*b)*(c[0]-2)*(c[0]-2)*(4 + c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*
	(5 + 2*b)*c[0]*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]
	*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0]
	 - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2])))/30240.;
		bf_mom[6] = (-6300*a + 6300*b + 63*(a - b)*(15*c[0]*c[0] + 20*c[1]
	*c[1] - 20*c[0]*c[2] + 28*c[2]*c[2]) + 9*(-35*(a*a*a - b*b*b)* (-4 + c[0]
	*c[0]) + 70*pow(a - b,2)*(a + b)*c[0]*c[1] - 28*(a - b)*(2*a*a + a*b + 
	2*b*b)*c[1]*c[1] + 28*(a - b)*((a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)
	*c[1])*c[2] - 4*(a - b)*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2]) - (a - b)*
	(21*(5*(a*a + a*b + b*b)*(-8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]
	*c[0]*c[1] + 12*(2*a*a + a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)
	*c[1]*c[1]*c[1]) - 18*(7*(a*a + 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)
	*c[0]*c[1] - 4*(4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b
	 + 15*b*b)*c[0] - 11*(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b
	 + 17*b*b)*c[2]*c[2]*c[2]) + 9*(a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]
	*c[2] + 8*c[2]*(7*c[1]*c[1] - 9*c[2]*c[2]) + 28*c[0]*(5*c[1]*c[1] 
	+ 7*c[2]*c[2])))/15120.;
		bf_mom[7] = -((a - b)*(2*a*a*(21*(5*(c[0]-2)*(c[0]-2)*(4 + c[0]) - 
	15*(-4 + c[0]*c[0])*c[1] + 24*c[0]*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*c[0] + 14*c[0]*c[1] - 4*(7 + 4*c[1]*c[1]))* c[2] + 36*(15*c[0] - 
	11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)*(c[0]-2)*
	(c[0]-2)*(4 + c[0]) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*
	c[0]*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]
	*c[0]) - 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 +
	 30*b)*c[0] + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2])
	 + a*(21*(5*(-3 + 2*b)*(c[0]-2)*(c[0]-2)*(4 + c[0]) + 30*(-4 + c[0]*c[0])
	*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5+6*b)
	*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49
	 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 8*(-81+94*b)*c[2]*c[2]*c[2])))/30240.;
		bf_mom[8] = ((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*(c[0]-2)*
	(c[0]-2)* (4 + c[0]) - 15*(a - b)*(a + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 
	+ 2*a*a + a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) 
	- 18*(7*(-5 + a*a + 3*a*b + b*b)*(-4 + c[0]*c[0]) + 14*(a - b)*(a + b)*c[0]
	*c[1] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((-49 + 15*a*a +
	 19*a*b + 15*b*b)*c[0] + 11*(-a*a + b*b)*c[1])*c[2]*c[2] - 8*(-81 + 17*a*a
	 + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/7560.;  
		break;
	case 103:
		bf_mom[0] = -((a - b)*(b*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])*c[0]*
	c[0]) + 30*(-1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 4*b)*(-1 + c[0])*c[1]*
	c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-2 + c[0])*c[0] - 
	28*(-1 + b)*(-1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])* c[2] + 36*((-49 + 
	30*b)*(-1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*(-81 + 34*b)*c[2]*
	c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])*c[0]*c[0]) + 30*(-2 + 
	c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1] + 24*c[1]*c[1]*c[1])
	 - 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 4*(7 + 2*b)*
	c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2]*c[2] - 8*
	(-81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*
	(5 + 5*c[1] + 2*c[2]) + 18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2]
	 + 30*c[2]*c[2]) + 4*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])
	*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[1] = -((a - b)*(2*a*a*(21*(5*(-1 + c[0])*pow(2 + c[0],2) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(-3 + 2*b)
	*(-1 + c[0])*pow(2 + c[0],2) + 30*(-1 + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 + 
	4*b)*(1 + c[0])*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)
	*c[0]*(2 + c[0]) - 28*(-1 + b)*(1 + c[0])*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*
	 c[2] + 36*((-49 + 30*b)*(1 + c[0]) + 22*(-1 + b)*c[1])* c[2]*c[2] - 8*
	(-81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*(-3 + 2*b)*(-1 + c[0])*
	pow(2 + c[0],2) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*
	c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1+c[0])
	*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0])+22*c[1])
	*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2])))/60480.;
		bf_mom[2] = -((a - b)*(2*a*a*(21*(5*(-1 + c[0])*pow(2 + c[0],2) - 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])* c[2] + 36*
	(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + b*(21*(5*(3+2*b)
	*(-1 + c[0])*pow(2 + c[0],2) + 30*(1 + b)*c[0]*(2 + c[0])*c[1] + 12*(5+4*b)
	*(1 + c[0])*c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*c[0]*
	(2 + c[0]) - 28*(1 + b)*(1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 
	36*((49 + 30*b)*(1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*
	c[2]*c[2]*c[2]) + a*(21*(5*(3 + 2*b)*(-1 + c[0])*pow(2 + c[0],2) - 30*c[0]*
	(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 
	18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*
	c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81+94*b)
	*c[2]*c[2]*c[2])))/60480.;
		bf_mom[3] = -((a - b)*(b*(21*(5*(3 + 2*b)*(20 + (-3 + c[0])*c[0]*
	c[0]) + 30*(1 + b)*(-2 + c[0])*c[0]*c[1] + 12*(5 + 4*b)*(-1 + c[0])*c[1]*
	c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-2 + c[0])*c[0] - 28*
	(1 + b)*(-1 + c[0])*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*
	(-1 + c[0]) + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) 
	+ a*(21*(5*(3 + 2*b)*(20 + (-3 + c[0])*c[0]*c[0]) - 30*(-2 + c[0])*c[0]*c[1]
	 + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5+6*b)
	*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 
	36*((49 + 38*b)*(-1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*
	c[2]) + 2*a*a*(105*c[0]*c[0]*c[0] - 63*c[0]*c[0]*(5 + 5*c[1] + 2*c[2]) + 
	18*c[0]*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*c[2]) + 4*
	(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15+11*c[1])
	*c[2]*c[2] - 34*c[2]*c[2]*c[2]))))/60480.;
		bf_mom[4] = ((a - b)*(b*(21*(5*(-3 + 2*b)*(-4 + c[0])*
	pow(2 + c[0],2) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]
	*c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0])
	 - 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)
	*c[0] + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) + a*
	(21*(5*(-3 + 2*b)*(-4 + c[0])*pow(2 + c[0],2) + 30*(-4 + c[0]*c[0])*c[1] + 
	12*(-5 + 2*b)*c[0]*c[1]*c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + 
	c[0]*c[0]) - 28*c[0]*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)
	*c[0] + 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*
	(105*c[0]*c[0]*c[0] - 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 72*(7 + 4*c[1]*
	c[1])*c[2] - 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] - 63*c[0]*c[0]*(5*c[1]
	 + 2*c[2]) + 36*c[0]*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 
	15*c[2]*c[2]))))/30240.;
		bf_mom[5] = ((a - b)*(21*(5*(-3 + a*a + a*b + b*b)*(-1 + c[0])* 
	pow(2 + c[0],2) - 15*(a - b)*(a + b)*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*a*a
	 + a*b + 2*b*b)*(1 + c[0])* c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) 
	- 18*(7*(-5 + a*a + 3*a*b + b*b)*c[0]*(2 + c[0]) + 14*(a - b)*(a + b)*
	(1 + c[0])*c[1] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*(
	(-49 + 15*a*a + 19*a*b + 15*b*b)*(1 + c[0]) - 11*(a - b)*(a + b)*c[1])*
	c[2]*c[2] - 8*(-81 + 17*a*a + 47*a*b + 17*b*b)*c[2]*c[2]*c[2]))/15120.;
		bf_mom[6] = ((a - b)*(b*(21*(5*(3 + 2*b)*(-4 + c[0])*pow(2+c[0],2)
	 + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*c[1]*c[1] + 24*
	(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) - 28*(1 + b)*
	c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0] + 22*
	(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) + a*(21*(5*
	(3 + 2*b)*(-4 + c[0])*pow(2 + c[0],2) - 30*(-4 + c[0]*c[0])*c[1] + 12*
	(5 + 2*b)*c[0]*c[1]*c[1] - 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + 
	c[0]*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)
	*c[0] - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 2*a*a*
	(105*c[0]*c[0]*c[0] - 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 72*(7 + 
	4*c[1]*c[1])*c[2] - 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] - 63*c[0]*
	c[0]*(5*c[1] + 2*c[2]) + 36*c[0]*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] 
	+ 15*c[2]*c[2]))))/30240.;
		bf_mom[7] = (63*(a - b)*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) + 9*(-35*(a*a*a - b*b*b)* (-4 + c[0]*c[0]) 
	+ 70*pow(a - b,2)*(a + b)*c[0]*c[1] - 28*(a - b)*(2*a*a + a*b + 2*b*b)*
	c[1]*c[1] + 28*(a - b)*((a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*
	c[2] - 4*(a - b)*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2]) + (a - b)*(21*(5*
	(a*a + a*b + b*b)*(8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]*c[0]*c[1] 
	+ 12*(2*a*a + a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*
	c[1]) - 18*(7*(a*a + 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)*c[0]*c[1] 
	- 4*(4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b + 15*b*b)
	*c[0] - 11*(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b + 17*b*b)
	*c[2]*c[2]*c[2]) - 9*(a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*c[2] + 
	28*c[0]*(5*c[1]*c[1] + 7*c[2]*c[2]) + 8*(35 + 7*c[1]*c[1]*c[2] - 
	9*c[2]*c[2]*c[2])))/15120.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + (
	(a - b)*(5*(a*a + a*b + b*b)*(2 + c[0]) - 5*(a - b)*(a + b)*c[1] - 2*
	(a*a + 3*a*b + b*b)*c[2]))/ 30. - ((a - b)*(21*(5*(a*a + a*b + b*b)* 
	(8 + c[0]*c[0]*c[0]) - 15*(a - b)*(a + b)*c[0]*c[0]*c[1] + 12*(2*a*a + 
	a*b + 2*b*b)*c[0]*c[1]*c[1] - 12*(a - b)*(a + b)*c[1]*c[1]*c[1]) - 18*
	(7*(a*a + 3*a*b + b*b)*c[0]*c[0] + 14*(a - b)*(a + b)*c[0]*c[1] - 4*
	(4*a*a - a*b + 4*b*b)*c[1]*c[1])*c[2] + 36*((15*a*a + 19*a*b + 15*b*b)*c[0]
	 - 11*(a - b)*(a + b)*c[1])*c[2]*c[2] - 8*(17*a*a + 47*a*b + 17*b*b)*
	c[2]*c[2]*c[2]))/7560. + ((a - b)*(35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*c[2] + 
	28*c[0]*(5*c[1]*c[1] + 7*c[2]*c[2]) + 8*(35 + 7*c[1]*c[1]*c[2] - 
	9*c[2]*c[2]*c[2])))/840.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
case 4:
switch (interface_type)
	{
	case 0:
		for(i=0;i<npts;i++)	{bf_mom[i]=gauss_mom[i];}
		break;
	case 2:
		bf_mom[0] = -((a-1)*(a-1)*(11*(21*(5*(1 + 2*a)*(20 + (-3 + c[0])
	*c[0]*c[0]) - 30*a*(-2 + c[0])*c[0]*c[1] + 12*(1 + 4*a)*(-1 + c[0])*c[1]*c[1] 
	- 24*a*c[1]*c[1]*c[1]) - 18*(7*(3 + 2*a)*(-2 + c[0])*c[0] + 28*a*(-1 + c[0])
	*c[1] - 4*(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*(-1 + c[0]) -22*a*c[1])
	*c[2]*c[2] - 8*(47 + 34*a)*c[2]*c[2]*c[2]) + 66*(9*(7*a*(-2 + c[0])*c[0] - 4*
	(3 + 4*a)*(-1 + c[0])*c[1] + 4*a*c[1]*c[1]) - 4*(45*a*(-1 + c[0]) - 2*(13 + 
	32*a)*c[1])*c[2] + 28*a*c[2]*c[2])*c[3] + 12*(11*((53 + 100*a)*(-1 + c[0]) - 
	94*a*c[1]) - 2*(321 + 328*a)*c[2])*c[3]*c[3]+4248*a*c[3]*c[3]*c[3]))/ 665280.;
		bf_mom[1] = -((-1 + a)*(5775*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(-25 + 
	20*c[1] - 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] + 
	6*c[1]*(-35 + 14*c[2] - 33*c[3]) + 3*c[2]*(49 + 60*c[3]) + c[3]*(126 + 
	253*c[3])) + 4*(2772*c[1]*c[1]*c[1] - 2530*c[2]*c[2]*c[2] + 297*c[1]*c[1]*
	(-21 + 10*c[2] - 4*c[3]) - 33*c[2]*c[2]*(237 + 28*c[3]) - 6*c[2]*c[3]*(990 + 
	977*c[3]) - 3*(-9625 + 2783*c[3]*c[3] + 708*c[3]*c[3]*c[3]) + 66*c[1]*
	(66*c[2]*c[2] + 14*c[2]*(-3 + 11*c[3]) + c[3]*(99 + 94*c[3]))) + 2*a*a*
	(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*
	(84*c[1]*c[1] + 90*c[2]*c[2] + c[2]*(42 - 90*c[3]) - 3*c[1]*(-35 + 14*c[2] + 
	24*c[3]) + c[3]*(-63 + 100*c[3])) - 4*(-5775 + 693*c[1]*c[1]*c[1] + 374*c[2]
	*c[2]*c[2] + 1650*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 99*c[1]*c[1]*(-14 + 8*c[2]
	 + 3*c[3]) - 33*c[2]*c[2]*(-45 + 7*c[3]) + 3*c[2]*c[3]*(-495 + 328*c[3]) + 
	33*c[1]*(33*c[2]*c[2] + c[3]*(-36 + 47*c[3]) - c[2]*(21 + 64*c[3])))) + a*
	(5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(25 + 10*c[1] + 22*c[2] - 6*c[3]) + 
	132*c[0]*(147*c[1]*c[1] + 261*c[2]*c[2] + c[2]*(231 - 90*c[3]) + 7*c[3]*
	(-9 + 37*c[3]) - 3*c[1]*(-35 + 14*c[2] + 78*c[3])) - 4*(1386*c[1]*c[1]*c[1] 
	+ 3850*c[2]*c[2]*c[2] + c[2]*c[2]*(8613 - 462*c[3]) - 99*c[1]*c[1]*(-49 + 
	10*c[2] + 6*c[3]) + 6*c[2]*c[3]*(-495 + 1291*c[3]) - 3*(9625 - 2849*c[3]*c[3]
	 + 354*c[3]*c[3]*c[3]) + 66*c[1]*(33*c[2]*c[2] + c[3]*(-117 + 47*c[3]) - c[2]*
	(21 + 142*c[3]))))))/665280.;
		bf_mom[2] = -((-1 + a)*(5775*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(25 + 
	20*c[1] - 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] + 
	6*c[1]*(35 + 14*c[2] - 33*c[3]) + 3*c[2]*(-49 + 60*c[3]) + c[3]*(-126 + 253*
	c[3])) + 4*(-5775 + 2772*c[1]*c[1]*c[1] - 2530*c[2]*c[2]*c[2] + c[2]*c[2]*
	(7821 - 924*c[3]) + 297*c[1]*c[1]*(21 + 10*c[2] - 4*c[3]) + 8349*c[3]*c[3] - 
	2124*c[3]*c[3]*c[3] - 6*c[2]*c[3]*(-990 + 977*c[3]) + 66*c[1]*(66*c[2]*c[2] +
	 14*c[2]*(3 + 11*c[3]) + c[3]*(-99 + 94*c[3]))) + 2*a*a*(1155*c[0]*c[0]*c[0] 
	- 693*c[0]*c[0]*(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 
	90*c[2]*c[2] - 6*c[2]*(7 + 15*c[3]) - 3*c[1]*(35 + 14*c[2] + 24*c[3]) + c[3]*
	(63 + 100*c[3])) - 4*(693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 99*c[1]*c[1]*
	(14 + 8*c[2] + 3*c[3]) - 33*c[2]*c[2]*(45 + 7*c[3]) + 3*c[2]*c[3]*(495 + 
	328*c[3]) - 3*(-385 + 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]
	*c[2] + c[2]*(21 - 64*c[3]) + c[3]*(36 + 47*c[3])))) + a*(5775*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]*(-25 + 10*c[1] + 22*c[2] - 6*c[3]) + 132*c[0]*(147*c[1]*c[1]
	 + 261*c[2]*c[2] - 3*c[2]*(77 + 30*c[3]) + 7*c[3]*(9 + 37*c[3]) - 3*c[1]*
	(35 + 14*c[2] + 78*c[3])) - 4*(5775 + 1386*c[1]*c[1]*c[1] + 3850*c[2]*c[2]
	*c[2] - 8547*c[3]*c[3] - 1062*c[3]*c[3]*c[3] - 99*c[1]*c[1]*(49 + 10*c[2] + 
	6*c[3]) - 33*c[2]*c[2]*(261 + 14*c[3]) + 6*c[2]*c[3]*(495 + 1291*c[3]) + 66*
	c[1]*(33*c[2]*c[2] + c[2]*(21 - 142*c[3]) + c[3]*(117 + 47*c[3]))))))/665280.;
		bf_mom[3] = -((a-1)*(a-1)*(11*(21*(5*(1 + 2*a)*(-1 + c[0])*
	pow(2 + c[0],2) - 30*a*c[0]*(2 + c[0])*c[1] + 12*(1 + 4*a)*(1 + c[0])*c[1]
	*c[1] - 24*a*c[1]*c[1]*c[1]) - 18*(7*(3 + 2*a)*c[0]*(2 + c[0]) + 28*a*(1 + 
	c[0])*c[1] - 4*(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*(1 + c[0]) - 
	22*a*c[1])*c[2]*c[2] - 8*(47 + 34*a)*c[2]*c[2]*c[2]) + 66*(9*(7*a*c[0]*(2 + 
	c[0]) - 4*(3 + 4*a)*(1 + c[0])*c[1] + 4*a*c[1]*c[1]) - 4*(45*a*(1 + c[0]) - 
	2*(13 + 32*a)*c[1])*c[2] + 28*a*c[2]*c[2])*c[3] + 12*(11*((53 + 100*a)*(1 + 
	c[0]) - 94*a*c[1]) - 2*(321 + 328*a)*c[2])*c[3]*c[3] + 
	4248*a*c[3]*c[3] *c[3]))/665280.;
		bf_mom[4] = ((a-1)*(a-1)*(2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(10 + 
	5*c[1] + 8*c[2] - 3*c[3]) + 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + c[2]*
	(168 - 90*c[3]) - 3*c[1]*(-35 + 14*c[2] + 60*c[3]) + c[3]*(-63 + 206*c[3])) - 
	4*(-11550 + 693*c[1]*c[1]*c[1] + 1408*c[2]*c[2]*c[2] + 3399*c[3]*c[3] - 531*
	c[3]*c[3]*c[3] - 297*c[1]*c[1]*(-7 + 2*c[2] + c[3]) - 33*c[2]*c[2]*(-102 + 
	7*c[3]) + 15*c[2]*c[3]*(-99 + 194*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(-90 
	+ 47*c[3]) - c[2]*(21 + 116*c[3]))) + a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	 (5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + c[2]*
	(42 - 90*c[3]) - 3*c[1]*(-35 + 14*c[2] + 24*c[3]) + c[3]*(-63 + 100*c[3])) - 
	4*(-5775 + 693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] + 1650*c[3]*c[3] - 531*
	c[3]*c[3]*c[3] - 99*c[1]*c[1]*(-14 + 8*c[2] + 3*c[3]) - 33*c[2]*c[2]*(-45 + 
	7*c[3]) + 3*c[2]*c[3]*(-495 + 328*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(-36 
	+ 47*c[3]) - c[2]*(21 + 64*c[3]))))))/166320.;
		bf_mom[5] = ((-1 + a)*(5775*c[0]*c[0]*c[0] + 1386*c[0]*c[0]*(10*c[1] 
	- 7*c[2] - 6*c[3]) + 132*c[0]*(-525 + 189*c[1]*c[1] + 237*c[2]*c[2] + 6*c[1]*
	(14*c[2] - 33*c[3]) + 180*c[2]*c[3] + 253*c[3]*c[3]) + 8*(1386*c[1]*c[1]*c[1]
	 - 1265*c[2]*c[2]*c[2] + 297*c[1]*c[1]*(5*c[2] - 2*c[3]) - 462*c[2]*c[2]*c[3]
	 + c[2]*(4851 - 2931*c[3]*c[3]) + 66*c[1]*(-105 + 33*c[2]*c[2] + 77*c[2]*c[3]
	 + 47*c[3]*c[3]) - 6*(1925 - 693*c[3] + 177*c[3]*c[3]*c[3])) + 2*a*a*(1155
	*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*
	c[1]*c[1] - 3*c[1]*(7*c[2] + 12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 
	10*c[3]*c[3])) - 4*(4620 + 693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] + 2079*c[3]
	 - 231*c[2]*c[2]*c[3] - 531*c[3]*c[3]*c[3] - 99*c[1]*c[1]*(8*c[2] + 3*c[3]) 
	+ 33*c[1]*(-105 + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3]) + 6*c[2]*(-231 
	+ 164*c[3]*c[3]))) + a*(5775*c[0]*c[0]*c[0] - 1386*c[0]*c[0]* (5*c[1] + 11*
	c[2] - 3*c[3]) + 132*c[0]*(-525 + 147*c[1]*c[1] + 261*c[2]*c[2] - 90*c[2]*c[3]
	 + 259*c[3]*c[3] - 6*c[1]*(7*c[2] + 39*c[3])) - 8*(11550 + 693*c[1]*c[1]*c[1]
	 + 1925*c[2]*c[2]*c[2] + 2079*c[3] - 231*c[2]*c[2]*c[3] - 531*c[3]*c[3]*c[3]
	 - 99*c[1]*c[1]*(5*c[2] + 3*c[3]) + 33*c[1]*(-105 + 33*c[2]*c[2] - 142*c[2]
	*c[3] + 47*c[3]*c[3]) + c[2]*(-7623 + 3873*c[3]*c[3])))))/ 332640.;
		bf_mom[6] = ((a-1)*(a-1)*(2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-10 + 
	5*c[1] + 8*c[2] - 3*c[3]) + 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] - 6*c[2]*
	(28 + 15*c[3]) - 3*c[1]*(35 + 14*c[2] + 60*c[3]) + c[3]*(63 + 206*c[3])) - 
	4*(2310 + 693*c[1]*c[1]*c[1] + 1408*c[2]*c[2]*c[2] - 3399*c[3]*c[3] - 531*
	c[3]*c[3]*c[3] - 297*c[1]*c[1]*(7 + 2*c[2] + c[3]) - 33*c[2]*c[2]*(102 + 
	7*c[3]) + 15*c[2]*c[3]*(99 + 194*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[2]*(21 - 
	116*c[3]) + c[3]*(90 + 47*c[3]))) + a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]* 
	(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] - 
	6*c[2]*(7 + 15*c[3]) - 3*c[1]*(35 + 14*c[2] + 24*c[3]) + c[3]*(63 + 100*c[3]))
	 - 4*(693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 99*c[1]*c[1]*(14 + 8*c[2] + 
	3*c[3]) - 33*c[2]*c[2]*(45 + 7*c[3]) + 3*c[2]*c[3]*(495 + 328*c[3]) - 3*
	(-385 + 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[2]*
	(21 - 64*c[3]) + c[3]*(36 + 47*c[3]))))))/166320.;
		bf_mom[7] = ((a-1)*(a-1)*(11*(21*(5*(1 + 2*a)*(-4 + c[0])*
	pow(2 + c[0],2) - 30*a*(-4 + c[0]*c[0])*c[1] + 12*(1 + 4*a)*c[0]*c[1]*c[1] 
	- 24*a*c[1]*c[1]*c[1]) - 18*(7*(3 + 2*a)*(-4 + c[0]*c[0]) + 28*a*c[0]*c[1] - 
	4*(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*c[0] - 22*a*c[1])*c[2]*c[2] - 
	8*(47 + 34*a)*c[2]*c[2]*c[2]) + 66*(4*c[1]*(-27*c[0] + 26*c[2]) + a*(-252 + 
	63*c[0]*c[0] - 36*c[0]*(4*c[1] + 5*c[2]) + 4*(9*c[1] + c[2])*(c[1] + 7*c[2])))
	*c[3] + 12*(11*(53 + 100*a)*c[0] - 2*(517*a*c[1] + (321 + 328*a)*c[2]))*c[3]
	*c[3] + 4248*a*c[3]*c[3]*c[3]))/332640.;
		bf_mom[8] = -((a-1)*(a-1)*(11*(21*(5*(2 + a)*(-4 + c[0])*
	pow(2 + c[0],2) - 15*(1 + a)*(-4 + c[0]*c[0])*c[1] + 12*(3 + 2*a)*c[0]*c[1]*
	c[1] - 12*(1 + a)*c[1]*c[1]*c[1]) - 18*(7*(4 + a)*(-4 + c[0]*c[0]) + 14*
	(1 + a)*c[0]*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] + 36*((34 + 15*a)*c[0] - 11*
	(1 + a)*c[1])*c[2]*c[2] - 8*(64 + 17*a)*c[2]*c[2]*c[2]) + 33*(9*(7*(1 + a)*
	(-4 + c[0]*c[0]) - 8*(5 + 2*a)*c[0]*c[1] + 4*(1 + a)*c[1]*c[1]) - 4*(45*
	(1 + a)*c[0] - 4*(29 + 16*a)*c[1])*c[2] + 28*(1 + a)*c[2]*c[2])*c[3] + 12*
	(11*(103 + 50*a)*c[0] - 517*(1 + a)*c[1] - 2*(485 + 164*a)*c[2])* c[3]*c[3] 
	+ 2124*(1 + a)*c[3]*c[3]*c[3]))/83160.;
		break;
	case 3:
		bf_mom[0] = ((a-1)*(a-1)*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5 + 6*c[2]) + 132*c[0]*(21*c[1]*c[1] + 63*c[2] + 57*c[2]*c[2] - 54*c[1]*c[4]
	 + 53*c[3]*c[3]) - 4*(1881*c[2]*c[2] + 1034*c[2]*c[2]*c[2] + 99*c[1]*c[1]*
	(7 + 2*c[2]) - 66*c[1]*(27 + 26*c[2])*c[4] + 1926*c[2]*c[3]*c[3] + 33*
	(-35 + 53*c[3]*c[3])) + 2*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1]
	 + 2*c[2] - 3*c[4]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + c[2]*(42 - 
	90*c[4]) - 3*c[1]*(-35 + 14*c[2] + 24*c[4]) + c[4]*(-63 + 100*c[4])) - 4*
	(693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 99*c[1]*c[1]*(-14 + 8*c[2]+3*c[4])
	 - 33*c[2]*c[2]*(-45 + 7*c[4]) + 3*c[2]*c[4]*(-495 + 328*c[4]) - 3*(385 - 
	550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[4]*(-36 + 
	47*c[4]) - c[2]*(21 + 64*c[4]))))))/665280.;
		bf_mom[1] = ((-1 + a)*(5775*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(-25 + 
	20*c[1] - 14*c[2] - 12*c[4]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] + 
	6*c[1]*(-35 + 14*c[2] - 33*c[4]) + 3*c[2]*(49 + 60*c[4]) + c[4]*(126 + 253*
	c[4])) + 4*(5775 + 2772*c[1]*c[1]*c[1] - 2530*c[2]*c[2]*c[2] + 297*c[1]*c[1]*
	(-21 + 10*c[2] - 4*c[4]) - 8349*c[3]*c[3] - 2124*c[3]*c[3]*c[3] - 33*c[2]*c[2]
	*(237 + 28*c[4]) - 6*c[2]*c[4]*(990 + 977*c[4]) + 66*c[1]*(66*c[2]*c[2] + 
	14*c[2]*(-3 + 11*c[4]) + c[4]*(99 + 94*c[4]))) + 2*a*a*(1155*c[0]*c[0]*c[0] -
	 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[4]) + 66*c[0]*(84*c[1]*c[1] + 90*
	c[2]*c[2] + c[2]*(42 - 90*c[4]) - 3*c[1]*(-35 + 14*c[2] + 24*c[4]) + c[4]*
	(-63 + 100*c[4])) - 4*(693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 99*c[1]*c[1]*
	(-14 + 8*c[2] + 3*c[4]) - 33*c[2]*c[2]*(-45 + 7*c[4]) + 3*c[2]*c[4]*(-495 + 
	328*c[4]) - 3*(385 - 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]
	*c[2] + c[4]*(-36 + 47*c[4]) - c[2]*(21 + 64*c[4])))) + a*(5775*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]* (25 + 10*c[1] + 22*c[2] - 6*c[4]) + 132*c[0]*(147*c[1]*c[1] 
	+ 261*c[2]*c[2] + c[2]*(231 - 90*c[4]) + 7*c[4]*(-9 + 37*c[4]) - 3*c[1]*(-35 
	+ 14*c[2] + 78*c[4])) - 4*(-5775 + 1386*c[1]*c[1]*c[1] + 3850*c[2]*c[2]*c[2] 
	+ c[2]*c[2]*(8613 - 462*c[4]) + 8547*c[3]*c[3] - 1062*c[3]*c[3]*c[3] - 99*c[1]
	*c[1]*(-49 + 10*c[2] + 6*c[4]) + 6*c[2]*c[4]*(-495 + 1291*c[4]) + 66*c[1]*
	(33*c[2]*c[2] + c[4]*(-117 + 47*c[4]) - c[2]*(21 + 142*c[4]))))))/665280.;
		bf_mom[2] = ((-1 + a)*(5775*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(25 + 20
	*c[1] - 14*c[2] - 12*c[4]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] + 6*c[1]*
	(35 + 14*c[2] - 33*c[4]) + 3*c[2]*(-49 + 60*c[4]) + c[4]*(-126 + 253*c[4])) + 
	4*(2772*c[1]*c[1]*c[1] - 2530*c[2]*c[2]*c[2] + c[2]*c[2]*(7821 - 924*c[4]) + 
	297*c[1]*c[1]*(21 + 10*c[2] - 4*c[4]) - 6*c[2]*c[4]*(-990 + 977*c[4]) - 3*
	(9625 - 2783*c[3]*c[3] + 708*c[3]*c[3]*c[3]) + 66*c[1]*(66*c[2]*c[2] + 14*c[2]
	*(3 + 11*c[4]) + c[4]*(-99 + 94*c[4]))) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*
	c[0]*c[0]*(-5 + 5*c[1] + 2*c[2] - 3*c[4]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]
	*c[2] - 6*c[2]*(7 + 15*c[4]) - 3*c[1]*(35 + 14*c[2] + 24*c[4]) + c[4]*
	(63 + 100*c[4])) - 4*(5775 + 693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 1650*
	c[3]*c[3] - 531*c[3]*c[3]*c[3] - 99*c[1]*c[1]*(14 + 8*c[2] + 3*c[4]) - 33*c[2]
	*c[2]*(45 + 7*c[4]) + 3*c[2]*c[4]*(495 + 328*c[4]) + 33*c[1]*(33*c[2]*c[2] + 
	c[2]*(21 - 64*c[4]) + c[4]*(36 + 47*c[4])))) + a*(5775*c[0]*c[0]*c[0] - 693*
	c[0]*c[0]* (-25 + 10*c[1] + 22*c[2] - 6*c[4]) + 132*c[0]*(147*c[1]*c[1] + 
	261*c[2]*c[2] - 3*c[2]*(77 + 30*c[4]) + 7*c[4]*(9 + 37*c[4]) - 3*c[1]*(35 + 
	14*c[2] + 78*c[4])) - 4*(1386*c[1]*c[1]*c[1] + 3850*c[2]*c[2]*c[2] - 99*c[1]
	*c[1]*(49 + 10*c[2] + 6*c[4]) - 33*c[2]*c[2]*(261 + 14*c[4]) + 6*c[2]*c[4]*
	(495 + 1291*c[4]) - 3*(-9625 + 2849*c[3]*c[3] + 354*c[3]*c[3]*c[3]) + 66*c[1]*
	(33*c[2]*c[2] + c[2]*(21 - 142*c[4]) + c[4]*(117 + 47*c[4]))))))/665280.;
		bf_mom[3] = ((a-1)*(a-1)*(11*(21*(5*(1 + 2*a)*(-20 + c[0]*c[0]*
	(3 + c[0])) - 30*a*c[0]*(2 + c[0])*c[1] + 12*(1 + 4*a)*(1 + c[0])*c[1]*c[1] 
	- 24*a*c[1]*c[1]*c[1]) - 18*(7*(3 + 2*a)*c[0]*(2 + c[0]) + 28*a*(1 + c[0])
	*c[1] - 4*(-1 + 8*a)*c[1]*c[1])*c[2] + 36*((19 + 30*a)*(1 + c[0]) - 22*a*c[1])
	*c[2]*c[2] - 8*(47 + 34*a)*c[2]*c[2]*c[2]) + 66*(9*(7*a*c[0]*(2 + c[0]) - 4*
	(3 + 4*a)*(1 + c[0])*c[1] + 4*a*c[1]*c[1]) - 4*(45*a*(1 + c[0]) - 2*(13+32*a)
	*c[1])*c[2] + 28*a*c[2]*c[2])*c[4] + 12*(11*((53 + 100*a)*(1 + c[0]) - 94*a*
	c[1]) - 2*(321 + 328*a)*c[2])* c[3]*c[3] + 4248*a*c[3]*c[3]*c[3]))/665280.;
		bf_mom[4] = -((a-1)*(a-1)*(2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(10 + 
	5*c[1] + 8*c[2] - 3*c[4]) + 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + c[2]*
	(168 - 90*c[4]) - 3*c[1]*(-35 + 14*c[2] + 60*c[4]) + c[4]*(-63 + 206*c[4])) - 
	4*(-2310 + 693*c[1]*c[1]*c[1] + 1408*c[2]*c[2]*c[2] + 3399*c[3]*c[3] - 531*
	c[3]*c[3]*c[3] - 297*c[1]*c[1]*(-7 + 2*c[2] + c[4]) - 33*c[2]*c[2]*(-102 + 
	7*c[4]) + 15*c[2]*c[4]*(-99 + 194*c[4]) + 33*c[1]*(33*c[2]*c[2] + c[4]*(-90 
	+ 47*c[4]) - c[2]*(21 + 116*c[4]))) + a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5 + 5*c[1] + 2*c[2] - 3*c[4]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + c[2]*
	(42 - 90*c[4]) - 3*c[1]*(-35 + 14*c[2] + 24*c[4]) + c[4]*(-63 + 100*c[4])) - 
	4*(693*c[1]*c[1]*c[1] + 374*c[2]*c[2]*c[2] - 99*c[1]*c[1]*(-14 + 8*c[2] + 
	3*c[4]) - 33*c[2]*c[2]*(-45 + 7*c[4]) + 3*c[2]*c[4]*(-495 + 328*c[4]) - 3*
	(385 - 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[4]*
	(-36 + 47*c[4]) - c[2]*(21 + 64*c[4]))))))/166320.;
		bf_mom[5] = -((-1 + a)*(5775*c[0]*c[0]*c[0] + 1386*c[0]*c[0]*(10*c[1] 
	- 7*c[2] - 6*c[4]) + 132*c[0]*(-525 + 189*c[1]*c[1] + 237*c[2]*c[2] + 6*c[1]*
	(14*c[2] - 33*c[4]) + 180*c[2]*c[4] + 253*c[3]*c[3]) + 8*(1386*c[1]*c[1]*c[1] 
	- 1265*c[2]*c[2]*c[2] + 297*c[1]*c[1]*(5*c[2] - 2*c[4]) - 462*c[2]*c[2]*c[4] 
	+ c[2]*(4851 - 2931*c[3]*c[3]) + 66*c[1]*(-105 + 33*c[2]*c[2] + 77*c[2]*c[4] 
	+ 47*c[3]*c[3]) + 6*(1925 + 693*c[4] - 177*c[3]*c[3]*c[3])) + a*(5775*c[0]*
	c[0]*c[0] - 1386*c[0]*c[0]*(5*c[1] + 11*c[2] - 3*c[4]) + 132*c[0]*(-525 + 
	147*c[1]*c[1] + 261*c[2]*c[2] - 90*c[2]*c[4] + 259*c[3]*c[3] - 6*c[1]*(7*c[2] 
	+ 39*c[4])) - 8*(-11550 + 693*c[1]*c[1]*c[1] - 7623*c[2] + 1925*c[2]*c[2]*c[2]
	 + 2079*c[4] - 231*c[2]*c[2]*c[4] + 3873*c[2]*c[3]*c[3] - 531*c[3]*c[3]*c[3] -
	 99*c[1]*c[1]*(5*c[2] + 3*c[4]) + 33*c[1]*(-105 + 33*c[2]*c[2] - 142*c[2]*c[4]
	 + 47*c[3]*c[3]))) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 
	2*c[2] - 3*c[4]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 12*c[4]) + 5*
	(-21 + 9*c[2]*c[2] - 9*c[2]*c[4] + 10*c[3]*c[3])) - 4*(-4620 + 693*c[1]*c[1]
	*c[1] + 374*c[2]*c[2]*c[2] + 2079*c[4] - 231*c[2]*c[2]*c[4] - 531*c[3]*c[3]
	*c[3] - 99*c[1]*c[1]*(8*c[2] + 3*c[4]) + 33*c[1]*(-105 + 33*c[2]*c[2] - 
	64*c[2]*c[4] + 47*c[3]*c[3]) + 6*c[2]*(-231 + 164*c[3]*c[3])))))/ 332640.;
		bf_mom[6] = -((-1 + a)*(-2310*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(-10 + 
	5*c[1] + 8*c[2] - 3*c[4]) - 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] - 6*c[2]*
	(28 + 15*c[4]) - 3*c[1]*(35 + 14*c[2] + 60*c[4]) + c[4]*(63 + 206*c[4])) + a*
	(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-5 + 6*c[2]) + 132*c[0]*(21*c[1]*c[1] - 
	63*c[2] + 57*c[2]*c[2] - 54*c[1]*c[4] + 53*c[3]*c[3]) - 4*(5775 - 1881*c[2]
	*c[2] + 1034*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(-7 + 2*c[2]) - 66*c[1]*(-27 + 
	26*c[2])*c[4] - 1749*c[3]*c[3] + 1926*c[2]*c[3]*c[3])) + 4*(11550 + 693*c[1]
	*c[1]*c[1] + 1408*c[2]*c[2]*c[2] - 3399*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 297*
	c[1]*c[1]*(7 + 2*c[2] + c[4]) - 33*c[2]*c[2]*(102 + 7*c[4]) + 15*c[2]*c[4]*
	(99 + 194*c[4]) + 33*c[1]*(33*c[2]*c[2] + c[2]*(21 - 116*c[4]) + c[4]*(90 + 
	47*c[4]))) + a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-5 + 5*c[1] + 2*c[2] - 
	3*c[4]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] - 6*c[2]*(7 + 15*c[4]) - 3*
	c[1]*(35 + 14*c[2] + 24*c[4]) + c[4]*(63 + 100*c[4])) - 4*(5775 + 693*c[1]*
	c[1]*c[1] + 374*c[2]*c[2]*c[2] - 1650*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 99*
	c[1]*c[1]*(14 + 8*c[2] + 3*c[4]) - 33*c[2]*c[2]*(45 + 7*c[4]) + 3*c[2]*c[4]*
	(495 + 328*c[4]) + 33*c[1]*(33*c[2]*c[2] + c[2]*(21 - 64*c[4]) + c[4]*(36 + 
	47*c[4]))))))/166320.;
		bf_mom[7] = -((a-1)*(a-1)*(11*(21*(5*(1 + 2*a)*pow(-2 + c[0],2)*(4 + 
	c[0]) - 30*a*(-4 + c[0]*c[0])*c[1] + 12*(1 + 4*a)*c[0]*c[1]*c[1] - 24*a*c[1]
	*c[1]*c[1]) - 18*(7*(3 + 2*a)*(-4 + c[0]*c[0]) + 28*a*c[0]*c[1] - 4*(-1 + 8*a)
	*c[1]*c[1])*c[2] + 36*((19 + 30*a)*c[0] - 22*a*c[1])*c[2]*c[2] - 8*(47 + 34*a)
	*c[2]*c[2]*c[2]) + 66*(4*c[1]*(-27*c[0] + 26*c[2]) + a*(-252 + 63*c[0]*c[0] - 
	36*c[0]*(4*c[1] + 5*c[2]) + 4*(9*c[1] + c[2])*(c[1] + 7*c[2])))*c[4] + 12*(11*
	(53 + 100*a)*c[0] - 2*(517*a*c[1] + (321 + 328*a)*c[2]))* c[3]*c[3] + 
	4248*a*c[3]*c[3]*c[3]))/332640.;
		bf_mom[8] = ((a-1)*(a-1)*(11*(21*(5*(2 + a)*pow(-2 + c[0],2)*
	(4 + c[0]) - 15*(1 + a)*(-4 + c[0]*c[0])*c[1] + 12*(3 + 2*a)*c[0]*c[1]*c[1] 
	- 12*(1 + a)*c[1]*c[1]*c[1]) - 18*(7*(4 + a)*(-4 + c[0]*c[0]) + 14*(1 + a)
	*c[0]*c[1] - 4*(3 + 4*a)*c[1]*c[1])*c[2] + 36*((34 + 15*a)*c[0] - 11*(1 + a)
	*c[1])*c[2]*c[2] - 8*(64 + 17*a)*c[2]*c[2]*c[2]) + 33*(9*(7*(1 + a)*(-4 + 
	c[0]*c[0]) - 8*(5 + 2*a)*c[0]*c[1] + 4*(1 + a)*c[1]*c[1]) - 4*(45*(1 + a)*c[0]
	 - 4*(29 + 16*a)*c[1])*c[2] + 28*(1 + a)*c[2]*c[2])*c[4] + 12*(11*(103 + 50*a)
	*c[0] - 517*(1 + a)*c[1] - 2*(485 + 164*a)*c[2])* c[3]*c[3] + 2124*(1 + a)
	*c[3]*c[3]*c[3]))/83160.;
		break;
	case 4:
		bf_mom[0] = -((1 + b)*(5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(25 + 
	20*c[1] + 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] - 
	6*c[1]*(-35 + 14*c[2] + 33*c[3]) - 3*c[2]*(-49 + 60*c[3]) + c[3]*(-126 + 
	253*c[3])) - 4*(-5775 + 2772*c[1]*c[1]*c[1] + 2530*c[2]*c[2]*c[2] + c[2]*c[2]*
	(7821 - 924*c[3]) + 8349*c[3]*c[3] - 2124*c[3]*c[3]*c[3] - 297*c[1]*c[1]*(-21 
	+ 10*c[2] + 4*c[3]) + 6*c[2]*c[3]*(-990 + 977*c[3]) + 66*c[1]*(66*c[2]*c[2] 
	- 14*c[2]*(3 + 11*c[3]) + c[3]*(-99 + 94*c[3]))) + 2*b*b*(1155*c[0]*c[0]
	*c[0] + 693*c[0]*c[0]*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1]
	 + 90*c[2]*c[2] + 3*c[1]*(-35 + 14*c[2] - 24*c[3]) + 6*c[2]*(7 + 15*c[3]) + 
	c[3]*(63 + 100*c[3])) + 4*(693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]
	*c[1]*(-14 + 8*c[2] - 3*c[3]) - 33*c[2]*c[2]*(45 + 7*c[3]) - 3*c[2]*c[3]*
	(495 + 328*c[3]) - 3*(-385 + 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*
	(33*c[2]*c[2] + c[3]*(36 + 47*c[3]) + c[2]*(-21 + 64*c[3])))) + b*(-5775*c[0]
	*c[0]*c[0] - 693*c[0]*c[0]*(-25 + 10*c[1] - 22*c[2] - 6*c[3]) - 132*c[0]*
	(147*c[1]*c[1] + 261*c[2]*c[2] + 3*c[1]*(-35 + 14*c[2] - 78*c[3]) + 3*c[2]*
	(77 + 30*c[3]) + 7*c[3]*(9 + 37*c[3])) - 4*(5775 + 1386*c[1]*c[1]*c[1] - 
	3850*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(-49 + 10*c[2] - 6*c[3]) - 8547*c[3]*c[3] 
	- 1062*c[3]*c[3]*c[3] - 33*c[2]*c[2]*(261 + 14*c[3]) - 6*c[2]*c[3]*(495 + 
	1291*c[3]) + 66*c[1]*(33*c[2]*c[2] + c[3]*(117 + 47*c[3]) + c[2]*(-21 + 
	142*c[3]))))))/665280.;
		bf_mom[1] = -((1+b)*(1+b)*(11*(21*(5*(-1 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) + 30*b*(-2 + c[0])*c[0]*c[1] + 12*(-1 + 4*b)*(-1 + c[0])*c[1]*c[1]
	 + 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-2 + c[0])*c[0] - 28*b*(-1 + c[0])
	*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*(-1 + c[0]) +22*b*c[1])
	*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*b*(-2 + c[0])*c[0] + 
	4*(-3 + 4*b)*(-1 + c[0])*c[1] + 4*b*c[1]*c[1]) - 4*(45*b*(-1 + c[0]) + 2*(-13
	 + 32*b)*c[1])*c[2] + 28*b*c[2]*c[2])*c[3] + 12*(11*((-53 + 100*b)*(-1+c[0]) 
	+ 94*b*c[1]) - 2*(-321 + 328*b)*c[2])*c[3]*c[3] - 4248*b*
	c[3]*c[3]*c[3]))/ 665280.;
		bf_mom[2] = -((1+b)*(1+b)*(11*(21*(5*(-1 + 2*b)*(-20 + c[0]*c[0]*
	(3 + c[0])) + 30*b*c[0]*(2 + c[0])*c[1] + 12*(-1 + 4*b)*(1 + c[0])*c[1]*c[1] 
	+ 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*c[0]*(2 + c[0]) - 28*b*(1 + c[0])
	*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*(1 + c[0]) + 22*b*c[1])
	*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*b*c[0]*(2 + c[0]) + 
	4*(-3 + 4*b)*(1 + c[0])*c[1] + 4*b*c[1]*c[1]) - 4*(45*b*(1 + c[0]) + 2*(-13 + 
	32*b)*c[1])*c[2] + 28*b*c[2]*c[2])*c[3] + 12*(11*((-53 + 100*b)*(1 + c[0]) + 
	94*b*c[1]) - 2*(-321 + 328*b)*c[2])*c[3]*c[3] - 4248*b*c[3]*c[3]*c[3]))
	/ 665280.;
		bf_mom[3] = -((1 + b)*(5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-25 + 
	20*c[1] + 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] - 
	6*c[1]*(35 + 14*c[2] + 33*c[3]) - 3*c[2]*(49 + 60*c[3]) + c[3]*(126+253*c[3]))
	 - 4*(2772*c[1]*c[1]*c[1] + 2530*c[2]*c[2]*c[2] - 297*c[1]*c[1]*(21 + 10*c[2] 
	+ 4*c[3]) - 33*c[2]*c[2]*(237 + 28*c[3]) + 6*c[2]*c[3]*(990 + 977*c[3]) - 3*
	(-9625 + 2783*c[3]*c[3] + 708*c[3]*c[3]*c[3]) + 66*c[1]*(66*c[2]*c[2] + c[2]*
	(42 - 154*c[3]) + c[3]*(99 + 94*c[3]))) + 2*b*b*(1155*c[0]*c[0]*c[0] + 693*
	c[0]*c[0]*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]
	*c[2] + 3*c[1]*(35 + 14*c[2] - 24*c[3]) + 6*c[2]*(-7 + 15*c[3]) + c[3]*(-63 + 
	100*c[3])) + 4*(-5775 + 693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]
	*c[1]*(14 + 8*c[2] - 3*c[3]) + 3*c[2]*(495 - 328*c[3])*c[3] + 1650*c[3]*c[3] 
	- 531*c[3]*c[3]*c[3] - 33*c[2]*c[2]*(-45 + 7*c[3]) + 33*c[1]*(33*c[2]*c[2] + 
	c[3]*(-36 + 47*c[3]) + c[2]*(21 + 64*c[3])))) - b*(5775*c[0]*c[0]*c[0] + 
	693*c[0]*c[0]*(25 + 10*c[1] - 22*c[2] - 6*c[3]) + 132*c[0]*(147*c[1]*c[1] + 
	261*c[2]*c[2] + 3*c[1]*(35 + 14*c[2] - 78*c[3]) + 3*c[2]*(-77 + 30*c[3]) + 
	7*c[3]*(-9 + 37*c[3])) + 4*(1386*c[1]*c[1]*c[1] - 3850*c[2]*c[2]*c[2] + c[2]
	*c[2]*(8613 - 462*c[3]) + 99*c[1]*c[1]*(49 + 10*c[2] - 6*c[3]) + 6*c[2]*
	(495 - 1291*c[3])*c[3] - 3*(9625 - 2849*c[3]*c[3] + 354*c[3]*c[3]*c[3]) + 
	66*c[1]*(33*c[2]*c[2] + c[3]*(-117 + 47*c[3]) + c[2]*
	(21 + 142*c[3]))))))/665280.;
		bf_mom[4] = ((1+b)*(1+b)*(-2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-10 
	+ 5*c[1] - 8*c[2] - 3*c[3]) - 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + 3*c[1]
	*(-35 + 14*c[2] - 60*c[3]) + 6*c[2]*(28 + 15*c[3]) + c[3]*(63 + 206*c[3])) - 
	4*(2310 + 693*c[1]*c[1]*c[1] - 1408*c[2]*c[2]*c[2] + 297*c[1]*c[1]*(-7 + 
	2*c[2] - c[3]) - 3399*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 33*c[2]*c[2]*(102 + 
	7*c[3]) - 15*c[2]*c[3]*(99 + 194*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(90 + 
	47*c[3]) + c[2]*(-21 + 116*c[3]))) + b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]* 
	(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + 
	3*c[1]*(-35 + 14*c[2] - 24*c[3]) + 6*c[2]*(7 + 15*c[3]) + c[3]*(63+100*c[3]))
	 + 4*(693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(-14 + 8*c[2] - 
	3*c[3]) - 33*c[2]*c[2]*(45 + 7*c[3]) - 3*c[2]*c[3]*(495 + 328*c[3]) - 3*(-385 
	+ 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(36 + 
	47*c[3]) + c[2]*(-21 + 64*c[3]))))))/166320.;
		bf_mom[5] = ((1+b)*(1+b)*(11*(21*(5*(-1 + 2*b)*pow(-2 + c[0],2)*(4 
	+ c[0]) + 30*b*(-4 + c[0]*c[0])*c[1] + 12*(-1 + 4*b)*c[0]*c[1]*c[1] + 24*b
	*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-4 + c[0]*c[0]) - 28*b*c[0]*c[1] - 4*
	(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*c[0] + 22*b*c[1])*c[2]*c[2] - 
	8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(4*c[1]*(-27*c[0] + 26*c[2]) + b*(-252 + 
	9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]
	*c[2]))*c[3] + 12*(11*(-53 + 100*b)*c[0] + 1034*b*c[1] + 642*c[2] - 656*b*
	c[2])* c[3]*c[3] - 4248*b*c[3]*c[3]*c[3]))/332640.;
		bf_mom[6] = ((1+b)*(1+b)*(-2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(10 + 
	5*c[1] - 8*c[2] - 3*c[3]) - 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + 3*c[1]*
	(35 + 14*c[2] - 60*c[3]) + 6*c[2]*(-28 + 15*c[3]) + c[3]*(-63 + 206*c[3])) - 
	4*(-11550 + 693*c[1]*c[1]*c[1] - 1408*c[2]*c[2]*c[2] + 297*c[1]*c[1]*(7 + 
	2*c[2] - c[3]) + 3399*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 33*c[2]*c[2]*(-102 + 
	7*c[3]) - 15*c[2]*c[3]*(-99 + 194*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(-90 
	+ 47*c[3]) + c[2]*(21 + 116*c[3]))) + b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*
	 (5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + 
	3*c[1]*(35 + 14*c[2] - 24*c[3]) + 6*c[2]*(-7 + 15*c[3]) + c[3]*(-63+100*c[3]))
	 + 4*(-5775 + 693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(14 + 
	8*c[2] - 3*c[3]) + 3*c[2]*(495 - 328*c[3])*c[3] + 1650*c[3]*c[3] - 531*c[3]*
	c[3]*c[3] - 33*c[2]*c[2]*(-45 + 7*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(-36 +
	 47*c[3]) + c[2]*(21 + 64*c[3]))))))/166320.;
		bf_mom[7] = ((1 + b)*(5775*c[0]*c[0]*c[0] - 1386*c[0]*c[0]*(10*c[1] + 
	7*c[2] - 6*c[3]) + 132*c[0]*(-525 + 189*c[1]*c[1] + 237*c[2]*c[2] - 180*c[2]
	*c[3] + 253*c[3]*c[3] - 6*c[1]*(14*c[2] + 33*c[3])) - 8*(1386*c[1]*c[1]*c[1] 
	+ 1265*c[2]*c[2]*c[2] - 462*c[2]*c[2]*c[3] - 297*c[1]*c[1]*(5*c[2] + 2*c[3]) 
	+ 66*c[1]*(-105 + 33*c[2]*c[2] - 77*c[2]*c[3] + 47*c[3]*c[3]) + 3*c[2]*(-1617 
	+ 977*c[3]*c[3]) - 6*(1925 - 693*c[3] + 177*c[3]*c[3]*c[3])) + 2*b*b*(1155*
	c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]
	*c[1] + 3*c[1]*(7*c[2] - 12*c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*
	c[3]*c[3])) + 4*(4620 + 693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]*c[1]
	*(8*c[2] - 3*c[3]) + 2079*c[3] - 231*c[2]*c[2]*c[3] - 531*c[3]*c[3]*c[3] + 
	c[2]*(1386 - 984*c[3]*c[3]) + 33*c[1]*(-105 + 33*c[2]*c[2] + 64*c[2]*c[3] + 
	47*c[3]*c[3]))) - b*(5775*c[0]*c[0]*c[0] + 1386*c[0]*c[0]* (5*c[1] - 11*c[2] 
	- 3*c[3]) + 132*c[0]*(-525 + 147*c[1]*c[1] + 261*c[2]*c[2] + 6*c[1]*(7*c[2] 
	- 39*c[3]) + 90*c[2]*c[3] + 259*c[3]*c[3]) + 8*(11550 + 693*c[1]*c[1]*c[1] - 
	1925*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(5*c[2] - 3*c[3]) + 2079*c[3] - 231*c[2]*
	c[2]*c[3] - 531*c[3]*c[3]*c[3] + c[2]*(7623 - 3873*c[3]*c[3]) + 33*c[1]*(-105 
	+ 33*c[2]*c[2] + 142*c[2]*c[3] + 47*c[3]*c[3])))))/332640.;
		bf_mom[8] = -((1+b)*(1+b)*(11*(21*(5*(-2 + b)*pow(-2 + c[0],2)*(4 + 
	c[0]) + 15*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-3 + 2*b)*c[0]*c[1]*c[1] + 
	12*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-4 + b)*(-4 + c[0]*c[0]) - 14*(-1 + b)
	*c[0]*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] + 36*((-34 + 15*b)*c[0] + 11*
	(-1 + b)*c[1])*c[2]*c[2] - 8*(-64 + 17*b)*c[2]*c[2]*c[2]) - 33*(9*(7*(-1 + b)
	*(-4 + c[0]*c[0]) + 8*(-5 + 2*b)*c[0]*c[1] + 4*(-1 + b)*c[1]*c[1]) - 4*(45*
	(-1 + b)*c[0] + 4*(-29 + 16*b)*c[1])*c[2] + 28*(-1 + b)*c[2]*c[2])*c[3] + 12*
	(11*(-103 + 50*b)*c[0] + 517*(-1 + b)*c[1] + 2*(485 - 164*b)*c[2])*c[3]*c[3] 
	- 2124*(-1 + b)*c[3]*c[3]*c[3]))/83160.;
		break;
	case 6:
		bf_mom[0] = ((1 + b)*(5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(25 + 
	20*c[1] + 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] - 
	6*c[1]*(-35 + 14*c[2] + 33*c[3]) - 3*c[2]*(-49 + 60*c[3]) + c[3]*(-126 + 
	253*c[3])) - 4*(2772*c[1]*c[1]*c[1] + 2530*c[2]*c[2]*c[2] + c[2]*c[2]*
	(7821 - 924*c[3]) - 297*c[1]*c[1]*(-21 + 10*c[2] + 4*c[3]) + 6*c[2]*c[3]*
	(-990 + 977*c[3]) - 3*(9625 - 2783*c[3]*c[3] + 708*c[3]*c[3]*c[3]) + 66*c[1]
	*(66*c[2]*c[2] - 14*c[2]*(3 + 11*c[3]) + c[3]*(-99 + 94*c[3]))) + 2*b*b*
	(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*
	c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + 3*c[1]*(-35 + 14*c[2] - 24*c[3]) + 
	6*c[2]*(7 + 15*c[3]) + c[3]*(63 + 100*c[3])) + 4*(5775 + 693*c[1]*c[1]*c[1] 
	- 374*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(-14 + 8*c[2] - 3*c[3]) - 1650*c[3]*c[3]
	 - 531*c[3]*c[3]*c[3] - 33*c[2]*c[2]*(45 + 7*c[3]) - 3*c[2]*c[3]*(495 + 
	328*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(36 + 47*c[3]) + c[2]*(-21 + 64*
	c[3])))) + b*(-5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-25 + 10*c[1] - 22*c[2] 
	- 6*c[3]) - 132*c[0]*(147*c[1]*c[1] + 261*c[2]*c[2] + 3*c[1]*(-35 + 14*c[2] 
	- 78*c[3]) + 3*c[2]*(77 + 30*c[3]) + 7*c[3]*(9 + 37*c[3])) - 4*(1386*c[1]*
	c[1]*c[1] - 3850*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(-49 + 10*c[2] - 6*c[3]) - 
	33*c[2]*c[2]*(261 + 14*c[3]) - 6*c[2]*c[3]*(495 + 1291*c[3]) - 3*(-9625 + 
	2849*c[3]*c[3] + 354*c[3]*c[3]*c[3]) + 66*c[1]*(33*c[2]*c[2] + c[3]*(117 + 
	47*c[3]) + c[2]*(-21 + 142*c[3]))))))/665280.;
		bf_mom[1] = (pow(1 + b,2)*(11*(21*(5*(-1 + 2*b)*(20 + (-3 + c[0])
	*c[0]*c[0]) + 30*b*(-2 + c[0])*c[0]*c[1] + 12*(-1 + 4*b)*(-1 + c[0])*c[1]
	*c[1] + 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-2 + c[0])*c[0] - 28*b*
	(-1 + c[0])*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*(-1 + c[0])
	 + 22*b*c[1])*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*b*
	(-2 + c[0])*c[0] + 4*(-3 + 4*b)*(-1 + c[0])*c[1] + 4*b*c[1]*c[1]) - 4*(45*b*
	(-1 + c[0]) + 2*(-13 + 32*b)*c[1])*c[2] + 28*b*c[2]*c[2])*c[3] + 12*(11*
	((-53 + 100*b)*(-1 + c[0]) + 94*b*c[1]) - 2*(-321 + 328*b)*c[2])*c[3]*c[3] 
	- 4248*b*c[3]*c[3]*c[3]))/ 665280.;
		bf_mom[2] = (pow(1 + b,2)*(11*(21*(5*(-1 + 2*b)*(-1 + c[0])*
	pow(2 + c[0],2) + 30*b*c[0]*(2 + c[0])*c[1] + 12*(-1 + 4*b)*(1 + c[0])
	*c[1]*c[1] + 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*c[0]*(2 + c[0]) - 
	28*b*(1 + c[0])*c[1] - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*
	(1 + c[0]) + 22*b*c[1])*c[2]*c[2] - 8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*
	(7*b*c[0]*(2 + c[0]) + 4*(-3 + 4*b)*(1 + c[0])*c[1] + 4*b*c[1]*c[1]) - 4*
	(45*b*(1 + c[0]) + 2*(-13 + 32*b)*c[1])*c[2] + 28*b*c[2]*c[2])*c[3] + 12*
	(11*((-53 + 100*b)*(1 + c[0]) + 94*b*c[1]) - 2*(-321 + 328*b)*c[2])
	*c[3]*c[3] - 4248*b*c[3]*c[3]*c[3]))/ 665280.;
		bf_mom[3] = ((1 + b)*(5775*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(-25 + 
	20*c[1] + 14*c[2] - 12*c[3]) + 132*c[0]*(189*c[1]*c[1] + 237*c[2]*c[2] - 
	6*c[1]*(35 + 14*c[2] + 33*c[3]) - 3*c[2]*(49 + 60*c[3]) + c[3]*(126 + 253*
	c[3])) - 4*(5775 + 2772*c[1]*c[1]*c[1] + 2530*c[2]*c[2]*c[2] - 8349*c[3]
	*c[3] - 2124*c[3]*c[3]*c[3] - 297*c[1]*c[1]*(21 + 10*c[2] + 4*c[3]) - 
	33*c[2]*c[2]*(237 + 28*c[3]) + 6*c[2]*c[3]*(990 + 977*c[3]) + 66*c[1]*
	(66*c[2]*c[2] + c[2]*(42 - 154*c[3]) + c[3]*(99 + 94*c[3]))) + 2*b*b*
	(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 
	66*c[0]*(84*c[1]*c[1] + 90*c[2]*c[2] + 3*c[1]*(35 + 14*c[2] - 24*c[3]) + 
	6*c[2]*(-7 + 15*c[3]) + c[3]*(-63 + 100*c[3])) + 4*(693*c[1]*c[1]*c[1] - 
	374*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(14 + 8*c[2] - 3*c[3]) + 3*c[2]*(495 - 
	328*c[3])*c[3] - 33*c[2]*c[2]*(-45 + 7*c[3]) - 3*(385 - 550*c[3]*c[3] + 
	177*c[3]*c[3]*c[3]) + 33*c[1]*(33*c[2]*c[2] + c[3]*(-36 + 47*c[3]) + c[2]*
	(21 + 64*c[3])))) - b*(5775*c[0]*c[0]*c[0] + 693*c[0]*c[0]* (25 + 10*c[1] 
	- 22*c[2] - 6*c[3]) + 132*c[0]*(147*c[1]*c[1] + 261*c[2]*c[2] + 3*c[1]*
	(35 + 14*c[2] - 78*c[3]) + 3*c[2]*(-77 + 30*c[3]) + 7*c[3]*(-9 + 37*c[3]))
	 + 4*(-5775 + 1386*c[1]*c[1]*c[1] - 3850*c[2]*c[2]*c[2] + c[2]*c[2]*(8613 
	- 462*c[3]) + 99*c[1]*c[1]*(49 + 10*c[2] - 6*c[3]) + 6*c[2]*(495-1291*c[3])
	*c[3] + 8547*c[3]*c[3] - 1062*c[3]*c[3]*c[3] + 66*c[1]*(33*c[2]*c[2] + 
	c[3]*(-117 + 47*c[3]) + c[2]*(21 + 142*c[3]))))))/665280.;
		bf_mom[4] = -(pow(1 + b,2)*(-2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(-10 + 5*c[1] - 8*c[2] - 3*c[3]) - 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + 
	3*c[1]*(-35 + 14*c[2] - 60*c[3]) + 6*c[2]*(28 + 15*c[3]) + c[3]*(63 + 206
	*c[3])) - 4*(11550 + 693*c[1]*c[1]*c[1] - 1408*c[2]*c[2]*c[2] + 297*c[1]
	*c[1]*(-7 + 2*c[2] - c[3]) - 3399*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 33*c[2]
	*c[2]*(102 + 7*c[3]) - 15*c[2]*c[3]*(99 + 194*c[3]) + 33*c[1]*(33*c[2]*c[2]
	 + c[3]*(90 + 47*c[3]) + c[2]*(-21 + 116*c[3]))) + b*(1155*c[0]*c[0]*c[0] + 
	693*c[0]*c[0]*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1] + 
	90*c[2]*c[2] + 3*c[1]*(-35 + 14*c[2] - 24*c[3]) + 6*c[2]*(7 + 15*c[3]) + 
	c[3]*(63 + 100*c[3])) + 4*(5775 + 693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 
	99*c[1]*c[1]*(-14 + 8*c[2] - 3*c[3]) - 1650*c[3]*c[3] - 531*c[3]*c[3]*c[3] 
	- 33*c[2]*c[2]*(45 + 7*c[3]) - 3*c[2]*c[3]*(495 + 328*c[3]) + 33*c[1]*
	(33*c[2]*c[2] + c[3]*(36 + 47*c[3]) + c[2]*(-21 + 64*c[3]))))))/166320.;
		bf_mom[5] = -(pow(1 + b,2)*(11*(21*(5*(-1 + 2*b)*(-4 + c[0])
	*pow(2 + c[0],2) + 30*b*(-4 + c[0]*c[0])*c[1] + 12*(-1 + 4*b)*c[0]*c[1]*c[1]
	 + 24*b*c[1]*c[1]*c[1]) - 18*(7*(-3 + 2*b)*(-4 + c[0]*c[0]) - 28*b*c[0]*c[1]
	 - 4*(1 + 8*b)*c[1]*c[1])*c[2] + 36*((-19 + 30*b)*c[0] +22*b*c[1])*c[2]*c[2]
	 - 8*(-47 + 34*b)*c[2]*c[2]*c[2]) - 66*(4*c[1]*(-27*c[0] + 26*c[2]) + b*
	(-252 + 9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 
	28*c[2]*c[2]))*c[3] + 12*(11*(-53 + 100*b)*c[0] + 1034*b*c[1] + 642*c[2] - 
	656*b*c[2])* c[3]*c[3] - 4248*b*c[3]*c[3]*c[3]))/332640.;
		bf_mom[6] = -(pow(1 + b,2)*(-2310*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(10 + 5*c[1] - 8*c[2] - 3*c[3]) - 66*c[0]*(126*c[1]*c[1] + 204*c[2]*c[2] + 
	3*c[1]*(35 + 14*c[2] - 60*c[3]) + 6*c[2]*(-28 + 15*c[3]) + c[3]*(-63 + 206
	*c[3])) - 4*(-2310 + 693*c[1]*c[1]*c[1] - 1408*c[2]*c[2]*c[2] + 297*c[1]
	*c[1]*(7 + 2*c[2] - c[3]) + 3399*c[3]*c[3] - 531*c[3]*c[3]*c[3] - 33*c[2]
	*c[2]*(-102 + 7*c[3]) - 15*c[2]*c[3]*(-99 + 194*c[3]) + 33*c[1]*(33*c[2]
	*c[2] + c[3]*(-90 + 47*c[3]) + c[2]*(21 + 116*c[3]))) + b*(1155*c[0]*c[0]
	*c[0] + 693*c[0]*c[0]*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 66*c[0]*(84*c[1]*c[1]
	 + 90*c[2]*c[2] + 3*c[1]*(35 + 14*c[2] - 24*c[3]) + 6*c[2]*(-7 + 15*c[3]) + 
	c[3]*(-63 + 100*c[3])) + 4*(693*c[1]*c[1]*c[1] - 374*c[2]*c[2]*c[2] + 99*c[1]
	*c[1]*(14 + 8*c[2] - 3*c[3]) + 3*c[2]*(495 - 328*c[3])*c[3] - 33*c[2]*c[2]*
	(-45 + 7*c[3]) - 3*(385 - 550*c[3]*c[3] + 177*c[3]*c[3]*c[3]) + 33*c[1]*
	(33*c[2]*c[2] + c[3]*(-36 + 47*c[3]) + c[2]*(21 + 64*c[3]))))))/166320.;
		bf_mom[7] = -((1 + b)*(5775*c[0]*c[0]*c[0] - 1386*c[0]*c[0]*(10*c[1]
	 + 7*c[2] - 6*c[3]) + 132*c[0]*(-525 + 189*c[1]*c[1] + 237*c[2]*c[2] - 
	180*c[2]*c[3] + 253*c[3]*c[3] - 6*c[1]*(14*c[2] + 33*c[3])) - 8*(1386*c[1]
	*c[1]*c[1] + 1265*c[2]*c[2]*c[2] - 462*c[2]*c[2]*c[3] - 297*c[1]*c[1]*
	(5*c[2] + 2*c[3]) + 66*c[1]*(-105 + 33*c[2]*c[2] - 77*c[2]*c[3] + 47*c[3]
	*c[3]) + 3*c[2]*(-1617 + 977*c[3]*c[3]) + 6*(1925 + 693*c[3] - 177*c[3]*c[3]
	*c[3])) + 2*b*b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 
	3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*c[3]) + 5*(-21 + 
	9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) + 4*(-4620 + 693*c[1]*c[1]*c[1] 
	- 374*c[2]*c[2]*c[2] + 99*c[1]*c[1]*(8*c[2] - 3*c[3]) + 2079*c[3] - 231*c[2]
	*c[2]*c[3] - 531*c[3]*c[3]*c[3] + c[2]*(1386 - 984*c[3]*c[3]) + 33*c[1]*
	(-105 + 33*c[2]*c[2] + 64*c[2]*c[3] + 47*c[3]*c[3]))) - b*(5775*c[0]*c[0]
	*c[0] + 1386*c[0]*c[0]*(5*c[1] - 11*c[2] - 3*c[3]) + 132*c[0]*(-525 + 
	147*c[1]*c[1] + 261*c[2]*c[2] + 6*c[1]*(7*c[2] - 39*c[3]) + 90*c[2]*c[3] 
	+ 259*c[3]*c[3]) + 8*(-11550 + 693*c[1]*c[1]*c[1] + 7623*c[2] - 1925*c[2]
	*c[2]*c[2] + 99*c[1]*c[1]*(5*c[2] - 3*c[3]) + 2079*c[3] - 231*c[2]*c[2]*c[3]
	 - 3873*c[2]*c[3]*c[3] - 531*c[3]*c[3]*c[3] + 33*c[1]*(-105 + 33*c[2]*c[2] 
	+ 142*c[2]*c[3] + 47*c[3]*c[3])))))/332640.;
		bf_mom[8] = (pow(1 + b,2)*(11*(21*(5*(-2 + b)*(-4 + c[0])*
	pow(2 + c[0],2) + 15*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-3 + 2*b)*c[0]
	*c[1]*c[1] + 12*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-4 + b)*(-4 + c[0]*c[0]) 
	- 14*(-1 + b)*c[0]*c[1] - 4*(-3 + 4*b)*c[1]*c[1])*c[2] + 36*((-34 + 15*b)
	*c[0] + 11*(-1 + b)*c[1])*c[2]*c[2] - 8*(-64 + 17*b)*c[2]*c[2]*c[2]) - 33*
	(9*(7*(-1 + b)*(-4 + c[0]*c[0]) + 8*(-5 + 2*b)*c[0]*c[1] + 4*(-1 + b)*c[1]
	*c[1]) - 4*(45*(-1 + b)*c[0] + 4*(-29 + 16*b)*c[1])*c[2] + 28*(-1 + b)*c[2]
	*c[2])*c[3] + 12*(11*(-103 + 50*b)*c[0] + 517*(-1 + b)*c[1] + 2*(485 - 164*b)
	*c[2])*c[3]*c[3] - 2124*(-1 + b)*c[3]*c[3]*c[3]))/ 83160.;
		break;
	case 22:
		bf_mom[0] = (11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 
	30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 
       	18*(7*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*
        c[2] - 396*(-1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 
	66*(63*c[0]*c[0] - 18*c[0]*(7 + 2*c[1] + 10*c[2]) + 4*(9*c[1]*
	(1 + c[1]) + (45 + 38*c[1])*c[2] + 7*c[2]*c[2]))*c[3] - 12*(517*
	(-1 + c[0] - 2*c[1]) - 14*c[2])*c[3]*c[3] - 4248*pow(c[3],3))/166320.;
		bf_mom[1] = (11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 
       	18*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])*
        c[2] - 396*(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 
    	66*(9*(7*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 
       	4*(45 + 45*c[0] - 38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*
	(1 + c[0] - 2*c[1]) - 14*c[2])*c[3]*c[3] - 4248*pow(c[3],3))/166320.;
		bf_mom[2] = (11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*c[0]*
	(2 + c[0]) + 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] + 2*c[1]) 
	- 14*c[2])*c[3]*c[3] + 4248*pow(c[3],3))/166320.;
		bf_mom[3] = (11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 30*
	(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 
	18*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 
	396*(-1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*
	(-2 + c[0])*c[0] + 4*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] 
	+ 38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(-1 + c[0] + 2*c[1]) - 
	14*c[2])*c[3]*c[3] + 4248*pow(c[3],3))/166320.;
		bf_mom[4] = (1155*pow(c[0],3) - 1386*c[0]*c[0]*(5*c[1] - c[2] 
	- 3*c[3]) + 132*c[0]*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 90*c[2]*c[3] 
	+ 47*c[3]*c[3] - 6*c[1]*(7*c[2] + 3*c[3])) - 8*(693*pow(c[1],3) - 
	143*pow(c[2],3) - 231*c[2]*c[2]*c[3] - 297*c[1]*c[1]*(3*c[2] + c[3]) + 
	21*c[2]*(33 + c[3]*c[3]) + 33*c[1]*(-105 + 33*c[2]*c[2] - 38*c[2]*c[3] 
	+ 47*c[3]*c[3]) - 3*(770 - 693*c[3] + 177*pow(c[3],3))))/83160.;
		bf_mom[5] = (11*(-105*(-20 + c[0]*c[0]*(3 + c[0])) - 252*
	(1 + c[0])*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 
       	684*(1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 264*c[1]*(27 + 27*c[0] - 
	26*c[2])*c[3] - 12*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3])/41580.;
		bf_mom[6] = (1155*pow(c[0],3) + 1386*c[0]*c[0]*(5*c[1] + c[2] - 
	3*c[3]) + 132*c[0]*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(7*c[2] 
	- 3*c[3]) + 90*c[2]*c[3] + 47*c[3]*c[3]) + 8*(2310 + 693*pow(c[1],3) + 
	143*pow(c[2],3) + 297*c[1]*c[1]*(3*c[2] - c[3]) + 2079*c[3] - 231*c[2]
	*c[2]*c[3] - 531*pow(c[3],3) - 21*c[2]*(33 + c[3]*c[3]) + 33*c[1]*
	(-105 + 33*c[2]*c[2] + 38*c[2]*c[3] + 47*c[3]*c[3])))/ 83160.;
		bf_mom[7] = (11*(-105*pow(-2 + c[0],2)*(1 + c[0]) - 252*
	(-1 + c[0])*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] - 
       	684*(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 264*c[1]*(-27 + 27*c[0] 
	- 26*c[2])*c[3] - 12*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3])/41580.;
		bf_mom[8] = (1155*pow(c[0],3) - 4158*c[0]*c[0]*c[2] + 
	132*c[0]*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 
	53*c[3]*c[3]) - 8*(-2310 + 517*pow(c[2],3) + 3*c[2]*(-693 + 
	33*c[1]*c[1] - 286*c[1]*c[3] + 321*c[3]*c[3])))/20790.;
		break;
	case 32:
		bf_mom[0] = (11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 30*
	(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 
       	18*(7*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*
        c[2] - 396*(-1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 
	66*(63*c[0]*c[0] - 18*c[0]*(7 + 2*c[1] + 10*c[2]) + 4*(9*c[1]*
	(1 + c[1]) + (45 + 38*c[1])*c[2] + 7*c[2]*c[2]))*c[3] - 12*(517*
	(-1 + c[0] - 2*c[1]) - 14*c[2])*c[3]*c[3] - 4248*pow(c[3],3))/166320.;
		bf_mom[1] = (11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 30*
	(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 
       18*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*
        c[2] - 396*(-1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*
	(9*(7*(-2 + c[0])*c[0] + 4*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*
	(-45 + 45*c[0] + 38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*
	(-1 + c[0] + 2*c[1]) - 14*c[2])*c[3]*c[3] + 4248*pow(c[3],3))/166320.;
		bf_mom[2] = (11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*c[0]*
	(2 + c[0]) + 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] + 2*c[1]) - 
	14*c[2])*c[3]*c[3] + 4248*pow(c[3],3))/166320.;
		bf_mom[3] = (11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 66*(9*(7*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] - 2*c[1]) - 
	14*c[2])*c[3]*c[3] - 4248*pow(c[3],3))/166320.;
		bf_mom[4] = (11*(-105*pow(-2 + c[0],2)*(1 + c[0]) - 252*
	(-1 + c[0])*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] - 
       	684*(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 264*c[1]*(-27 + 
	27*c[0] - 26*c[2])*c[3] - 12*(583*(-1 + c[0]) - 642*c[2])
	*c[3]*c[3])/41580.;
		bf_mom[5] = (1155*pow(c[0],3) + 1386*c[0]*c[0]*(5*c[1] + 
	c[2] - 3*c[3]) + 132*c[0]*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*
	(7*c[2] - 3*c[3]) + 90*c[2]*c[3] + 47*c[3]*c[3]) + 8*(2310 + 693*
	pow(c[1],3) + 143*pow(c[2],3) + 297*c[1]*c[1]*(3*c[2] - c[3]) + 
	2079*c[3] - 231*c[2]*c[2]*c[3] - 531*pow(c[3],3) - 21*c[2]*(33 + 
	c[3]*c[3]) + 33*c[1]*(-105 + 33*c[2]*c[2] + 38*c[2]*c[3] + 
	47*c[3]*c[3])))/ 83160.;
		bf_mom[6] = (11*(-105*(-20 + c[0]*c[0]*(3 + c[0])) - 252*
	(1 + c[0])*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 
	684*(1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 264*c[1]*(27 + 27*c[0] - 
	26*c[2])*c[3] - 12*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3])/41580.;
		bf_mom[7] = (1155*pow(c[0],3) - 1386*c[0]*c[0]*(5*c[1] - c[2] 
	- 3*c[3]) + 132*c[0]*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 90*c[2]*c[3] 
	+ 47*c[3]*c[3] - 6*c[1]*(7*c[2] + 3*c[3])) - 8*(693*pow(c[1],3) - 
	143*pow(c[2],3) - 231*c[2]*c[2]*c[3] - 297*c[1]*c[1]*(3*c[2] + c[3]) 
	+ 21*c[2]*(33 + c[3]*c[3]) + 33*c[1]*(-105 + 33*c[2]*c[2] - 38*c[2]
	*c[3] + 47*c[3]*c[3]) - 3*(770 - 693*c[3] + 177*pow(c[3],3))))/83160.;
		bf_mom[8] = (1155*pow(c[0],3) - 4158*c[0]*c[0]*c[2] + 132*c[0]*
	(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 
	8*(-2310 + 517*pow(c[2],3) + 3*c[2]*(-693 + 33*c[1]*c[1] - 286*c[1]*c[3]
	 + 321*c[3]*c[3])))/20790.;
		break;
	case 100:
		bf_mom[0] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])*
	c[0]*c[0]) + 30*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1] 
	+ 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1]
	 + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*
	c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(-2 + c[0])*c[0] + 4*
	(-7 + 6*b)*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + (-90 + 52*b)*
	c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*(-1 + c[0])+94*c[1])
	 - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*(99*(7*(-5*
	(20 + (-3 + c[0])*c[0]*c[0]) - 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 + c[0])*
	c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1]
	 - 4*c[1]*c[1])*c[2] - 4*(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*
	c[2]) + 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*
	(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*(-153 + 153*c[0]+94*c[1]
	 - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(100 + 5*(-3+c[0])
	*c[0]*c[0] + 15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*
	c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2]
	 + 36*(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*
	(-2 + c[0])*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 
	64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 
	328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3])) + 2*a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) 
	- 42*(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*
	c[3]*c[3]) + 4*(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] 
	- 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[1] = -((a - b)*(a*(11*(21*(5*(3 + 2*b)*(20 + (-3 + c[0])*c[0]*
	c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] - 24*
	c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*
	(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(-1 + c[0]) - 22*c[1])*c[2]*c[2]
	 - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] - 18*c[0]*(7 + 2*(7+6*b)*
	c[1] + 10*c[2]) + 4*(9*c[1]*(7 + 6*b + c[1]) + (45 + (90 + 52*b)*c[1])*c[2] 
	+ 7*c[2]*c[2]))*c[3] + 12*(11*((153 + 106*b)*(-1 + c[0]) - 94*c[1]) - 2*
	(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(99*(700 + 35*c[0]*
	c[0]*c[0] + 28*c[1]*c[1]*(-5 + 2*c[1]) + 35*c[0]*c[0]*(-3 + 2*c[1] - 2*c[2])
	 + 56*(-1 + c[1])*c[1]*c[2] + 4*(-49 + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]
	 + 28*c[0]*(5*(-1 + c[1])*c[1] + (5 + 2*c[1])*c[2] + 7*c[2]*c[2])) - 66*(9*
	(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 
	2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*(-153 + 153*c[0] + 94*c[1]-118*c[2])
	*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(100 + 5*(-3 + c[0])*c[0]*c[0]
	 + 15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 
	+ 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3]
	 - 2124*c[3]*c[3]*c[3])) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 
	5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*
	c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*
	(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 
	11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 
	64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] 
	+ 531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[2] = -((a - b)*(2*a*a*(11*(21* (5*(-1 + c[0])*(2+c[0])*(2+c[0])
	 - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 
	15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + c[0])
	 - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2] + 28*
	c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3] + 2124*
	c[3]*c[3]*c[3]) + a*(11*(21*(5*(3 + 2*b)*(-1 + c[0])*(2+c[0])*(2+c[0]) - 30*
	c[0]*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*c[1]*c[1])
	 - 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 + 2*b)*c[1]*
	c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*
	c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) - 4*(7 + 6*b)*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(45 + 45*c[0] - 2*(45 + 26*b)*c[1])*c[2] + 28*c[2]*c[2])*
	c[3] + 12*(11*((153 + 106*b)*(1 + c[0]) - 94*c[1]) - 2*(649 + 642*b)*c[2])*
	c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(99*(7*(5*(-1 + c[0])*(2+c[0])*(2+c[0])
	 + 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 
	14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])* c[2] + 4*(49 + 49*
	c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 66*(9*(7*c[0]*(2 + c[0]) + 
	28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*
	c[2])*c[3] + 132*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 4248*c[3]*
	c[3]*c[3] + 2*b*(11*(21*(5*(-1 + c[0])*(2+c[0])*(2+c[0]) + 15*c[0]*(2+c[0])*
	c[1] + 24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) 
	- 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*
	c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*
	c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 
	12*(550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 
	2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[3] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-1 + c[0])*(2+c[0])*
	(2+c[0]) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] + 
	24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 
	4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0])+22*c[1])*c[2]*c[2]
	 - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*c[0]*(2 + c[0]) + 4*(-7 + 6*b)*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + (-90 + 52*b)*c[1])*c[2] +
	 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*(1 + c[0]) + 94*c[1]) - 2*(-649
	 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + 2*a*a*(11*(21*(5*(-1+c[0])
	*(2+c[0])*(2+c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 
	12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1]-16*c[1]*c[1])
	*c[2] + 36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*
	(7*c[0]*(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 
	64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 
	328*c[2])*c[3]*c[3] + 2124*c[3]*c[3]*c[3]) + b*(99*(7*(20 - 5*c[0]*c[0]*
	(3 + c[0]) - 10*c[0]*(2 + c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]
	*c[1]) + 14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])* c[2] - 4*
	(49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 
	+ c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 
	28*c[2]*c[2])*c[3] - 132*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 
	4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(5*(-1 + c[0])*(2+c[0])*(2+c[0]) + 15*c[0]*
	(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*
	(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 
	11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]
	*c[2])*c[3] + 12*(550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*
	c[3]*c[3]*c[3]))))/665280.;
		bf_mom[4] = (99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) - 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)*
	(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*(-12*a*b*
	c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 10*c[2]))*
	c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 27720*(-a + b - ((a - b)*
	(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] - 27*c[2]*c[2]
	 + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 49*c[2]*c[2] - 
	42*c[1]*c[3] + 51*c[3]*c[3])))/840.) + (a - b)*(a*b*(9240 + 1155*c[0]*c[0]*
	c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*
	c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] 
	+ 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2]
	 - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*
	c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*
	c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])
	*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3]-531*c[3]*c[3]*c[3]))
	 + a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 
	132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*
	c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*
	c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 
	7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[5] = ((a - b)*(a*(11*(21*(5*(3 + 2*b)*(-4 + c[0])*(2+c[0])*
	(2+c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*(5 + 2*b)*c[0]*c[1]*c[1] - 24*c[1]*
	c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)*
	c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0] - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)*
	c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] + 36*(-7 + c[1]*c[1]) + 8*(45 + 26*b)*
	c[1]*c[2] + 28*c[2]*c[2] - 36*c[0]*((7 + 6*b)*c[1] + 5*c[2]))*c[3] + 12*
	(11*(153 + 106*b)*c[0] - 2*(517*c[1] + (649 + 642*b)*c[2]))* c[3]*c[3] + 
	4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(3 + 2*b)*(-4 + c[0])*(2+c[0])*(2+c[0])
	 + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*c[1]*c[1] + 24*(1+b)
	*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) - 28*(1 + b)*c[0]*c[1] - 
	4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0] + 22*(1+b)*c[1])*c[2]*c[2]
	 - 8*(81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(1 + b)*(-4 + c[0]*c[0]) + 4*
	(7 + 4*b)*c[0]*c[1] + 4*(1 + b)*c[1]*c[1]) - 4*(45*(1 + b)*c[0] + 2*(45+32*b)
	*c[1])*c[2] + 28*(1 + b)*c[2]*c[2])*c[3] + 12*(11*(153 + 100*b)*c[0] + 1034*
	(1 + b)*c[1] - 2*(649 + 328*b)*c[2])*c[3]*c[3] - 4248*(1 + b)*c[3]*c[3]*c[3])
	 + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 
	132*c[0]*(3*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 15*c[2]*c[2]) - 9*(4*c[1] + 
	5*c[2])*c[3] + 50*c[3]*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) - 
	18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])
	*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 332640.;
		bf_mom[6] = (-99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) + 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)*
	(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*(-12*a*b*
	c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 10*c[2]))
	*c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 27720*(-a + b - ((a - b)*
	(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1]-27*c[2]*c[2]
	 + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 49*c[2]*c[2] - 
	42*c[1]*c[3] + 51*c[3]*c[3])))/840.) + (a - b)*(a*b*(9240 + 1155*c[0]*c[0]*
	c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*
	c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3]
	 + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 
	2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 
	36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] 
	+ 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 
	7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]
	*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 
	3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*
	c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 
	72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + 
	c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[7] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-4 + c[0])*(2+c[0])*
	(2+c[0]) + 30*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 24*c[1]
	*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 4*(7 + 2*b)*
	c[1]*c[1])*c[2] + 36*((-49 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)
	*c[2]*c[2]*c[2]) - 66*(63*c[0]*c[0] + 36*c[1]*c[1] + 36*c[0]*((-7 + 6*b)*c[1]
	 - 5*c[2]) - 8*(-45 + 26*b)*c[1]*c[2] + 28*(-9 + c[2]*c[2]))*c[3] + 12*(11*
	(-153 + 106*b)*c[0] + 2*(517*c[1] + (649 - 642*b)*c[2]))* c[3]*c[3] - 4248*
	c[3]*c[3]*c[3]) + b*(11*(21*(5*(-3 + 2*b)*(-4 + c[0])*(2+c[0])*(2+c[0]) +30*
	(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]*c[1]*c[1] + 24*(-1 + b)*
	c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0]) - 28*(-1 + b)*c[0]*c[1] 
	- 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)*c[0] + 22*(-1 + b)*c[1])
	*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(-1 + b)*(-4 + c[0]
	*c[0]) + 4*(-7 + 4*b)*c[0]*c[1] + 4*(-1 + b)*c[1]*c[1]) - 4*(45*(-1 + b)*c[0]
	 + 2*(-45 + 32*b)*c[1])*c[2] + 28*(-1 + b)*c[2]*c[2])*c[3] + 12*(11*(-153 
	+ 100*b)*c[0] + 1034*(-1 + b)*c[1] + 2*(649 - 328*b)*c[2])*c[3]*c[3] - 4248*
	(-1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1]
	 + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) - 18*(7 + 
	4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(-63 + 
	(9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] +
	 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 12*c[3]) + 
	5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3])))))/ 332640.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + 
	((a - b)*(a*b*(10 + 5*c[0] - 6*c[2]) + b*b*(10 + 5*c[0] + 5*c[1] - 2*c[2] - 
	3*c[3]) + a*a*(10 + 5*c[0] - 5*c[1] - 2*c[2] + 3*c[3])))/30. + (a - b + 
	((a - b)*(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] - 
	27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 49*
	c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/840.)/3. - ((a - b)*(a*b*(9240 
	+ 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 57*
	c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 517*c[2]*
	c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] + 693*c[0]
	*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] +
	 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 
	63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2]-34*c[2]*c[2]*c[2])
	 - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3]
	 - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] +
	 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 
	36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210+63*c[1]*c[1]*c[1]
	 - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] 
	+ c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/ 83160.;
		break;
	case 101:
		bf_mom[0] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) + 30*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1]
	 + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])
	*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1 + c[0])+22*c[1])
	*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(-2 + c[0])*c[0] + 
	4*(-7 + 6*b)*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + (-90 + 
	52*b)*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*(-1 + c[0]) 
	+ 94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*
	(99*(7*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 10*(-2 + c[0])*c[0]*c[1] - 20*
	(-1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*(-2 + c[0])*c[0] - 4*(-1 
	+ c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] 
	+ 72*c[2]*c[2]*c[2]) + 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*
	(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 
	2*b*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*(-2 + c[0])*c[0]*c[1] + 24*
	(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7-8*c[1])
	*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] 
	- 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*
	(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3])) 
	+ 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3])
	 + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*
	(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*(105 - 63*c[1]*c[1]*
	(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 
	34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])*c[2] + 7*c[2]
	*c[2])* c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[1] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-20 + c[0]*c[0]*
	(3 + c[0])) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] 
	+ 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1]
	 + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*
	c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*c[0]*(2 + c[0]) + 4*
	(-7 + 6*b)*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + (-90 + 52*b)
	*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(11*((-153 + 106*b)*(1 + c[0]) + 
	94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + 2*a*a*
	(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*
	(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] - 11*c[1])*c[2]*
	c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + c[0]) - 16*(1 + c[0])*c[1]
	 + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 
	12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3] + 2124*c[3]*c[3]*c[3]) 
	+ b*(99*(7*(100 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + c[0])*c[1] - 20*
	(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) - 4*(1 + 
	c[0])*c[1] - 4*c[1]*c[1])* c[2] - 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 
	72*c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]
	*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*(153 + 
	153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 2*b* (11*
	(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1+c[0])
	*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1]
	 - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]
	*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) 
	- 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0]
	 + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[2] = ((a - b)*(2*a*a*(11*(21* (-100 + 5*c[0]*c[0]*(3 + c[0])
	 - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) 
	- 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + 
	c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3]
	 + 2124*c[3]*c[3]*c[3]) + a*(11*(21*(5*(3 + 2*b)*(-20 + c[0]*c[0]*(3+c[0]))
	 - 30*c[0]*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*
	c[1]*c[1]) - 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*
	(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2]
	 - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) - 4*(7 + 6*b)*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 2*(45 + 26*b)*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(11*((153 + 106*b)*(1 + c[0]) - 94*c[1]) - 2*
	(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(99*(7*(-100 + 
	5*c[0]*c[0]*(3 + c[0]) + 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] 
	+ 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*
	c[1])* c[2] + 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 
	66*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1+c[0]
	 + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*(153 + 153*c[0] + 94*c[1] - 
	118*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 2*b* (11*(21*(-100 + 5*c[0]*
	c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 
	12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*
	c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 
	33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 
	45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] + 517*c[1]
	 - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[3] = ((a - b)*(a*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] 
	- 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2 + c[0])*c[0] + 28*(-1 + c[0])
	*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(-1 + c[0])-22*c[1])
	*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] - 18*c[0]*
	(7 + 2*(7 + 6*b)*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + 6*b + c[1]) + (45 + 
	(90 + 52*b)*c[1])*c[2] + 7*c[2]*c[2]))*c[3] + 12*(11*((153 + 106*b)*
	(-1 + c[0]) - 94*c[1]) - 2*(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]
	*c[3]) + b*(99*(7*(5*pow(-2 + c[0],2)*(1 + c[0]) + 10*(-2 + c[0])*c[0]*c[1]
	 + 20*(-1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 14*(5*(-2 + c[0])*c[0] - 
	4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(-49 + 49*c[0] + 22*c[1])*c[2]
	*c[2] - 72*c[2]*c[2]*c[2]) - 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])
	*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] 
	+ 132*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]
	 + 2*b*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*(-2 + c[0])*c[0]*c[1] + 
	24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 
	8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*
	c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1+c[0])
	*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3]-2124*c[3]*c[3]
	*c[3])) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] 
	- 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 90*c[2]
	*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*(105 - 
	63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]
	*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])*c[2] 
	+ 7*c[2]*c[2])* c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]
	*c[3]*c[3]))))/665280.;
		bf_mom[4] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 
	24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 
	4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 
	8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(63*c[0]*c[0] + 36*c[1]*c[1] + 36*c[0]*
	((-7 + 6*b)*c[1] - 5*c[2]) - 8*(-45 + 26*b)*c[1]*c[2] + 28*(-9 + c[2]*c[2]))
	*c[3] + 12*(11*(-153 + 106*b)*c[0] + 2*(517*c[1] + (649 - 642*b)*c[2]))*
	 c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(-3 + 2*b)*pow(-2+c[0],2)
	*(4 + c[0]) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]*c[1]
	*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0]) - 
	28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)*c[0]
	 + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*
	(-1 + b)*(-4 + c[0]*c[0]) + 4*(-7 + 4*b)*c[0]*c[1] + 4*(-1 + b)*c[1]*c[1]) 
	- 4*(45*(-1 + b)*c[0] + 2*(-45 + 32*b)*c[1])*c[2] + 28*(-1 + b)*c[2]*c[2])
	*c[3] + 12*(11*(-153 + 100*b)*c[0] + 1034*(-1 + b)*c[1] + 2*(649 - 328*b)
	*c[2])*c[3]*c[3] - 4248*(-1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]
	*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20+3*c[1]*
	(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*
	c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*
	(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*
	c[1] - 3*c[1]*(7*c[2] + 12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 
	10*c[3]*c[3])))))/ 332640.;
		bf_mom[5] = (99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) - 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)
	*(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*
	(-12*a*b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 
	10*c[2]))*c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 83160*((-a+b)/3.
	 + ((a - b)* (105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1]
	 - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/2520.) - (a - b)*(a*b*(-9240 
	+ 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 
	57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] 
	+ 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1]
	 - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]
	*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*
	(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 
	34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] 
	+ 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[6] = -((a - b)*(a*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*
	(4 + c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*(5 + 2*b)*c[0]*c[1]*c[1] - 
	24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*c[0]*c[1] + 4*
	(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0] - 22*c[1])*c[2]*c[2] - 8*
	(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] + 36*(-7 + c[1]*c[1]) + 8*
	(45 + 26*b)*c[1]*c[2] + 28*c[2]*c[2] - 36*c[0]*((7 + 6*b)*c[1] + 5*c[2]))
	*c[3] + 12*(11*(153 + 106*b)*c[0] - 2*(517*c[1] + (649 + 642*b)*c[2]))*c[3]
	*c[3] + 4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*c[1]*c[1]
	 + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) - 28*(1+b)
	*c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0] + 22*(1+b)
	*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(1 + b)*(-4 + 
	c[0]*c[0]) + 4*(7 + 4*b)*c[0]*c[1] + 4*(1 + b)*c[1]*c[1]) - 4*(45*(1 + b)
	*c[0] + 2*(45 + 32*b)*c[1])*c[2] + 28*(1 + b)*c[2]*c[2])*c[3] + 12*(11*
	(153 + 100*b)*c[0] + 1034*(1 + b)*c[1] - 2*(649 + 328*b)*c[2])*c[3]*c[3] - 
	4248*(1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(3*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] 
	+ 15*c[2]*c[2]) - 9*(4*c[1] + 5*c[2])*c[3] + 50*c[3]*c[3]) + 4*(-11*(21*
	(-20 + 3*c[1]*(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]
	*c[2] + 34*c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3]
	 - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 332640.;
		bf_mom[7] = (a - b + ((a - b)*(21*(5*(a*a + a*b + b*b)* (-4 + c[0]
	*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) 
	- 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*c[2] + 12*(15*a*a + 
	19*a*b + 15*b*b)*c[2]*c[2] + 18*(-12*a*b*c[1] + a*a*(7*c[0] - 8*c[1] - 
	10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 10*c[2]))*c[3] + 4*(50*a*a + 53*a*b + 
	50*b*b)*c[3]*c[3]))/1260. - ((a - b)*(105*c[0]*c[0] - 140*c[0]*c[2] + 4*
	(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/420. + 2*(
	(-a + b)/3. + ((a - b)* (105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*
	(21*c[1]*c[1] - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*
	c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/2520.) - ((a - b)
	*(a*b*(-9240 + 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*
	c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*
	c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*
	c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*
	c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*
	c[3]) + 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*
	c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*
	(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*
	c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*
	c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*
	c[3]) + 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*
	c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*
	(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/41580.)/4.;
		bf_mom[8] = -a + b + ((a - b)*c[0])/2. - ((a - b)*c[2])/3. - 
	((a - b)*(a*b*(-10 + 5*c[0] - 6*c[2]) + b*b*(5*c[0] + 5*c[1] - 2*(5 + c[2])
	 - 3*c[3]) + a*a*(5*c[0] - 5*c[1] - 2*(5 + c[2]) + 3*c[3])))/30. + (a - b 
	- ((a - b)*(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] 
	- 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/840.)/3. + ((a - b)*(a*b*
	(-9240 + 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1]
	 + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0]
	 + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1]
	 - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*
	c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 
	4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 
	34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1]
	 + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 83160.;
		break;
	case 102:
		bf_mom[0] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) + 30*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*c[1]
	 + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*(-1 + c[0])
	*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1 + c[0])+22*c[1])
	*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(-2 + c[0])*c[0] + 
	4*(-7 + 6*b)*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + (-90 + 
	52*b)*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*(-1 + c[0]) 
	+ 94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*
	(99*(7*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 
	+ c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*(-2 + c[0])*c[0] - 4*(-1 + 
	c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 
	72*c[2]*c[2]*c[2]) + 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*
	(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 
	2*b*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*(-2 + c[0])*c[0]*c[1] + 24*
	(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7-8*c[1])
	*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 
	136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*
	(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3])) + 
	2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 
	66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 
	8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*(105 - 63*c[1]*c[1]*(2 + 
	c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*
	c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])*c[2] + 7*c[2]*c[2])*
	 c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[1] = ((a - b)*(a*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*
	(1 + c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] 
	- 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2 + c[0])*c[0] + 28*(-1 + c[0])
	*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(-1 + c[0])-22*c[1])
	*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] - 18*c[0]*
	(7 + 2*(7 + 6*b)*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + 6*b + c[1]) + (45 + (90 
	+ 52*b)*c[1])*c[2] + 7*c[2]*c[2]))*c[3] + 12*(11*((153 + 106*b)*(-1 + c[0])
	 - 94*c[1]) - 2*(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*
	(99*(7*(5*pow(-2 + c[0],2)*(1 + c[0]) + 10*(-2 + c[0])*c[0]*c[1] + 20*(-1 
	+ c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 14*(5*(-2 + c[0])*c[0] - 4*(-1 + 
	c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] - 
	72*c[2]*c[2]*c[2]) - 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*
	(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 
	2*b*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*(-2 + c[0])*c[0]*c[1] + 
	24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 
	8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*
	c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])
	*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])
	*c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]
	*c[3]*c[3])) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 
	2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 
	90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*(105
	 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])
	*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])
	*c[2] + 7*c[2]*c[2])* c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[2] = ((a - b)*(2*a*a*(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) 
	- 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 +
	 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + 
	c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3]
	 + 2124*c[3]*c[3]*c[3]) + a*(11*(21*(5*(3 + 2*b)*(-20 + c[0]*c[0]*(3+c[0]))
	 - 30*c[0]*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*
	c[1]*c[1]) - 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 
	+ 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 
	8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) - 4*(7 + 6*b)*(1 
	+ c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 2*(45 + 26*b)*c[1])*c[2] 
	+ 28*c[2]*c[2])*c[3] + 12*(11*((153 + 106*b)*(1 + c[0]) - 94*c[1]) - 2*
	(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(99*(7*(-100 + 
	5*c[0]*c[0]*(3 + c[0]) + 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1]
	 + 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]
	*c[1])* c[2] + 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 
	66*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1+c[0]
	 + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*(153 + 153*c[0] + 94*c[1] - 
	118*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(-100 + 5*c[0]*c[0]
	*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 12*c[1]
	*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*c[1])
	*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*
	(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0]
	 + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] + 517*c[1] - 
	328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[3] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-20 + c[0]*c[0]*
	(3 + c[0])) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]*c[1] 
	+ 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1]
	 + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*
	c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*c[0]*(2 + c[0]) + 4*
	(-7 + 6*b)*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + (-90 + 52*b)
	*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(11*((-153 + 106*b)*(1 + c[0]) + 
	94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + 2*a*a*
	(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*
	(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2]
	 - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + c[0]) - 16*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*
	(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3] + 2124*c[3]*c[3]*c[3]) + 
	b*(99*(7*(100 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + c[0])*c[1] - 20*(1 
	+ c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) - 4*(1+c[0])
	*c[1] - 4*c[1]*c[1])* c[2] - 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]
	*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*c[1])
	 - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*(153 + 153*c[0]
	 + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(-100
	 + 5*c[0]*c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]
	*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 
	16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]
	*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*
	(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] + 
	517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[4] = (-99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) + 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)
	*(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*
	(-12*a*b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 
	10*c[2]))*c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 83160*((-a+b)/3.
	 + ((a - b)* (105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1]
	 - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/2520.) - (a - b)*(a*b*(-9240 
	+ 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 
	57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 517*
	c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] + 
	693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*
	c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 
	4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 
	34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1]
	 - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]
	*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*
	(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 
	34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1]
	 + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[5] = -((a - b)*(a*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*(4 
	+ c[0]) - 30*(-4 + c[0]*c[0])*c[1] + 12*(5 + 2*b)*c[0]*c[1]*c[1] - 24*c[1]
	*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)
	*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0] - 22*c[1])*c[2]*c[2] - 8*(81 + 94*b)
	*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] + 36*(-7 + c[1]*c[1]) + 8*(45 + 26*b)
	*c[1]*c[2] + 28*c[2]*c[2] - 36*c[0]*((7 + 6*b)*c[1] + 5*c[2]))*c[3] + 12*
	(11*(153 + 106*b)*c[0] - 2*(517*c[1] + (649 + 642*b)*c[2]))* c[3]*c[3] + 
	4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(3 + 2*b)*pow(-2 + c[0],2)*(4 + c[0]) + 
	30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*c[1]*c[1] + 24*(1 + b)
	*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) - 28*(1 + b)*c[0]*c[1]
	 - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*c[0] + 22*(1 + b)*c[1])
	*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(1 + b)*(-4 + c[0]
	*c[0]) + 4*(7 + 4*b)*c[0]*c[1] + 4*(1 + b)*c[1]*c[1]) - 4*(45*(1 + b)*c[0]
	 + 2*(45 + 32*b)*c[1])*c[2] + 28*(1 + b)*c[2]*c[2])*c[3] + 12*(11*(153 + 
	100*b)*c[0] + 1034*(1 + b)*c[1] - 2*(649 + 328*b)*c[2])*c[3]*c[3] - 4248*
	(1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(3*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2]
	 + 15*c[2]*c[2]) - 9*(4*c[1] + 5*c[2])*c[3] + 50*c[3]*c[3]) + 4*(-11*(21*
	(-20 + 3*c[1]*(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]
	*c[2] + 34*c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3]
	 - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 332640.;
		bf_mom[6] = (99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) - 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)
	*(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*(-12*a*
	b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1]+10*c[2]))
	*c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 83160*((-a + b)/3. + (
	(a - b)* (105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] 
	- 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/2520.) - (a - b)*(a*b*(-9240
	 + 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 
	57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0]
	 + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3])
	 + 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2]
	 - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*
	c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3])
	 + 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2]
	 + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*
	c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[7] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)*
	(4 + c[0]) + 30*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 
	24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 
	4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 
	8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(63*c[0]*c[0] + 36*c[1]*c[1] + 36*c[0]*
	((-7 + 6*b)*c[1] - 5*c[2]) - 8*(-45 + 26*b)*c[1]*c[2] + 28*(-9+c[2]*c[2]))
	*c[3] + 12*(11*(-153 + 106*b)*c[0] + 2*(517*c[1] + (649 - 642*b)*c[2]))*
	 c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(-3 + 2*b)*pow(-2 + c[0],2)
	*(4 + c[0]) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]*c[1]
	*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0]) - 
	28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)*c[0]
	 + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*
	(7*(-1 + b)*(-4 + c[0]*c[0]) + 4*(-7 + 4*b)*c[0]*c[1] + 4*(-1+b)*c[1]*c[1])
	 - 4*(45*(-1 + b)*c[0] + 2*(-45 + 32*b)*c[1])*c[2] + 28*(-1 + b)*c[2]*c[2])
	*c[3] + 12*(11*(-153 + 100*b)*c[0] + 1034*(-1 + b)*c[1] + 2*(649 - 328*b)*
	c[2])*c[3]*c[3] - 4248*(-1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*
	(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]
	*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1]
	 + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 
	3*c[1]*(7*c[2] + 12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 
	10*c[3]*c[3])))))/ 332640.;
		bf_mom[8] = -a + b + ((a - b)*c[0])/2. - ((a - b)*c[2])/3. - (
	(a - b)*(a*b*(-10 + 5*c[0] - 6*c[2]) + b*b*(5*c[0] + 5*c[1] - 2*(5 + c[2]) 
	- 3*c[3]) + a*a*(5*c[0] - 5*c[1] - 2*(5 + c[2]) + 3*c[3])))/30. + (a - b 
	- ((a - b)*(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] 
	- 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/840.)/3. + ((a - b)*(a*b*
	(-9240 + 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1]
	 + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0]
	 + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*
	c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	+ 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*
	(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 83160.;
		break;
	case 103:
		bf_mom[0] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(20 + (-3 + c[0])
	*c[0]*c[0]) + 30*(-2 + c[0])*c[0]*c[1] + 12*(-5 + 2*b)*(-1 + c[0])*c[1]*
	c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-2 + c[0])*c[0] - 28*
	(-1 + c[0])*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(-1+c[0])
	 + 22*c[1])*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*(-2+c[0])
	*c[0] + 4*(-7 + 6*b)*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 
	(-90 + 52*b)*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*
	(-1 + c[0]) + 94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]
	*c[3]) + b*(99*(7*(-5*(20 + (-3 + c[0])*c[0]*c[0]) - 10*(-2 + c[0])*c[0]
	*c[1] - 20*(-1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*(-2 + c[0])
	*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49 + 49*c[0] +22*c[1])
	*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1+c[0])
	*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] 
	- 132*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]
	 + 2*b*(11*(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 15*(-2 + c[0])*c[0]*c[1] + 
	24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 
	8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*c[2]
	*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])
	*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])
	*c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]
	*c[3]*c[3])) + 2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 
	2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 
	90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*
	(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 +11*c[1])
	*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])
	*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[1] = -((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-1 + c[0])*
	pow(2 + c[0],2) + 30*c[0]*(2 + c[0])*c[1] + 12*(-5 + 2*b)*(1 + c[0])*c[1]
	*c[1] + 24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*c[0]*(2 + c[0]) - 28*(1+c[0])
	*c[1] + 4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*(1 + c[0])+22*c[1])
	*c[2]*c[2] - 8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(9*(7*c[0]*(2 + c[0]) + 
	4*(-7 + 6*b)*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + (-90+52*b)
	*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(11*((-153 + 106*b)*(1 + c[0]) + 
	94*c[1]) - 2*(-649 + 642*b)*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + 2*a*a*
	(11*(21*(5*(-1 + c[0])*pow(2 + c[0],2) - 15*c[0]*(2 + c[0])*c[1] + 24*
	(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2]
	 - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + c[0]) - 16*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*
	(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3] + 2124*c[3]*c[3]*c[3]) + 
	b*(99*(7*(20 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + c[0])*c[1] - 20*(1 + 
	c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])
	*c[1] - 4*c[1]*c[1])* c[2] - 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*
	c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*
	c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 132*(153 + 
	153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3] + 2*b*(11*
	(21*(5*(-1 + c[0])*pow(2 + c[0],2) + 15*c[0]*(2 + c[0])*c[1] + 24*(1+c[0])
	*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])
	*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*
	c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*
	c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 
	550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[2] = -((a - b)*(2*a*a*(11*(21*(5*(-1 + c[0])*pow(2 + c[0],2)
	 - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 
	15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2+c[0])
	 - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2] + 
	28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3] + 
	2124*c[3]*c[3]*c[3]) + a*(11*(21*(5*(3 + 2*b)*(-1 + c[0])*pow(2 + c[0],2) 
	- 30*c[0]*(2 + c[0])*c[1] + 12*(5 + 2*b)*(1 + c[0])*c[1]*c[1] - 24*c[1]*
	c[1]*c[1]) - 18*(7*(5 + 6*b)*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*(-7 
	+ 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2]*c[2] - 
	8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(9*(7*c[0]*(2 + c[0]) - 4*(7 + 6*b)*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 2*(45 + 26*b)*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(11*((153 + 106*b)*(1 + c[0]) - 94*c[1]) - 2*
	(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(99*(7*(5*(-1+c[0])
	*pow(2 + c[0],2) + 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] + 
	8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1]-4*c[1]*c[1])
	* c[2] + 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 66*(9*
	(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 
	2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*(153 + 153*c[0] + 94*c[1]-118*c[2])
	*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 2*b*(11*(21*(5*(-1 + c[0])*pow(2+c[0],2)
	 + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) 
	- 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*
	(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*
	(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] + 517*c[1] - 328*c[2])*
	c[3]*c[3] - 2124*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[3] = -((a - b)*(a*(11*(21*(5*(3 + 2*b)*(20 + (-3 + c[0])*
	c[0]*c[0]) - 30*(-2 + c[0])*c[0]*c[1] + 12*(5 + 2*b)*(-1 + c[0])*c[1]*c[1] 
	- 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-2 + c[0])*c[0] + 28*(-1 + c[0])
	*c[1] + 4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*(-1 + c[0])-22*c[1])
	*c[2]*c[2] - 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] - 18*c[0]*
	(7 + 2*(7 + 6*b)*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + 6*b + c[1]) + (45 + 
	(90 + 52*b)*c[1])*c[2] + 7*c[2]*c[2]))*c[3] + 12*(11*((153 + 106*b)*(-1 + 
	c[0]) - 94*c[1]) - 2*(649 + 642*b)*c[2])*c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + 
	b*(99*(700 + 35*c[0]*c[0]*c[0] + 28*c[1]*c[1]*(-5 + 2*c[1]) + 35*c[0]*c[0]*
	(-3 + 2*c[1] - 2*c[2]) + 56*(-1 + c[1])*c[1]*c[2] + 4*(-49 + 22*c[1])*c[2]*
	c[2] - 72*c[2]*c[2]*c[2] + 28*c[0]*(5*(-1 + c[1])*c[1] + (5 + 2*c[1])*c[2] 
	+ 7*c[2]*c[2])) - 66*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*
	c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 132*(-153 + 
	153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 4248*c[3]*c[3]*c[3] + 2*b*(11*
	(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 15*(-2 + c[0])*c[0]*c[1] + 24*(-1+c[0])
	*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 
	14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*
	c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]
	*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(-550 
	+ 550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3])) + 
	2*a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) 
	+ 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*
	(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 4*(11*(525 - 63*c[1]*c[1]*
	(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 11*c[1])*c[2]*c[2] - 34*
	c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 64*c[1])*c[2]+7*c[2]*c[2])
	*c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] +
	 531*c[3]*c[3]*c[3]))))/665280.;
		bf_mom[4] = ((a - b)*(a*(11*(21*(5*(-3 + 2*b)*(-4 + c[0])*
	pow(2+c[0],2) + 30*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 2*b)*c[0]*c[1]*c[1] + 
	24*c[1]*c[1]*c[1]) - 18*(7*(-5 + 6*b)*(-4 + c[0]*c[0]) - 28*c[0]*c[1] + 
	4*(7 + 2*b)*c[1]*c[1])*c[2] + 36*((-49 + 38*b)*c[0] + 22*c[1])*c[2]*c[2] - 
	8*(-81 + 94*b)*c[2]*c[2]*c[2]) - 66*(63*c[0]*c[0] + 36*c[1]*c[1] + 36*c[0]*
	((-7 + 6*b)*c[1] - 5*c[2]) - 8*(-45 + 26*b)*c[1]*c[2] + 28*(-9+c[2]*c[2]))
	*c[3] + 12*(11*(-153 + 106*b)*c[0] + 2*(517*c[1] + (649 - 642*b)*c[2]))*
	 c[3]*c[3] - 4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(-3 + 2*b)*(-4 + c[0])
	*pow(2 + c[0],2) + 30*(-1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(-5 + 4*b)*c[0]*
	c[1]*c[1] + 24*(-1 + b)*c[1]*c[1]*c[1]) - 18*(7*(-5 + 2*b)*(-4 + c[0]*c[0])
	 - 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1])*c[2] + 36*((-49 + 30*b)
	*c[0] + 22*(-1 + b)*c[1])*c[2]*c[2] - 8*(-81 + 34*b)*c[2]*c[2]*c[2]) - 66*
	(9*(7*(-1 + b)*(-4 + c[0]*c[0]) + 4*(-7 + 4*b)*c[0]*c[1] + 4*(-1 + b)*c[1]
	*c[1]) - 4*(45*(-1 + b)*c[0] + 2*(-45 + 32*b)*c[1])*c[2] + 28*(-1 + b)*c[2]
	*c[2])*c[3] + 12*(11*(-153 + 100*b)*c[0] + 1034*(-1 + b)*c[1] + 2*(649 - 
	328*b)*c[2])*c[3]*c[3] - 4248*(-1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*
	c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 
	3*c[1]*(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 
	34*c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*
	(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*
	c[1] - 3*c[1]*(7*c[2] + 12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 
	10*c[3]*c[3])))))/ 332640.;
		bf_mom[5] = (-99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) + 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)
	*(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*
	(-12*a*b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 
	10*c[2]))*c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 27720*(-a + b - 
	((a - b)*(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] - 
	27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/840.) + (a - b)*(a*b*(9240 
	+ 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 
	57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] 
	+ 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	- 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*
	c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0]
	 - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	+ 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*
	c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[6] = ((a - b)*(a*(11*(21*(5*(3 + 2*b)*(-4 + c[0])*
	pow(2 + c[0],2) - 30*(-4 + c[0]*c[0])*c[1] + 12*(5 + 2*b)*c[0]*c[1]*c[1] 
	- 24*c[1]*c[1]*c[1]) - 18*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*c[0]*c[1] + 
	4*(-7 + 2*b)*c[1]*c[1])*c[2] + 36*((49 + 38*b)*c[0] - 22*c[1])*c[2]*c[2] 
	- 8*(81 + 94*b)*c[2]*c[2]*c[2]) + 66*(63*c[0]*c[0] + 36*(-7 + c[1]*c[1]) + 
	8*(45 + 26*b)*c[1]*c[2] + 28*c[2]*c[2] - 36*c[0]*((7 + 6*b)*c[1] + 5*c[2]))
	*c[3] + 12*(11*(153 + 106*b)*c[0] - 2*(517*c[1] + (649 + 642*b)*c[2]))*
	 c[3]*c[3] + 4248*c[3]*c[3]*c[3]) + b*(11*(21*(5*(3 + 2*b)*(-4 + c[0])
	*pow(2 + c[0],2) + 30*(1 + b)*(-4 + c[0]*c[0])*c[1] + 12*(5 + 4*b)*c[0]*
	c[1]*c[1] + 24*(1 + b)*c[1]*c[1]*c[1]) - 18*(7*(5 + 2*b)*(-4 + c[0]*c[0]) 
	- 28*(1 + b)*c[0]*c[1] - 4*(7 + 8*b)*c[1]*c[1])*c[2] + 36*((49 + 30*b)*
	c[0] + 22*(1 + b)*c[1])*c[2]*c[2] - 8*(81 + 34*b)*c[2]*c[2]*c[2]) - 66*(9*
	(7*(1 + b)*(-4 + c[0]*c[0]) + 4*(7 + 4*b)*c[0]*c[1] + 4*(1 + b)*c[1]*c[1]) 
	- 4*(45*(1 + b)*c[0] + 2*(45 + 32*b)*c[1])*c[2] + 28*(1 + b)*c[2]*c[2])*c[3]
	 + 12*(11*(153 + 100*b)*c[0] + 1034*(1 + b)*c[1] - 2*(649 + 328*b)*c[2])*
	c[3]*c[3] - 4248*(1 + b)*c[3]*c[3]*c[3]) + 2*a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(3*(-35 + 14*c[1]*c[1] 
	- 7*c[1]*c[2] + 15*c[2]*c[2]) - 9*(4*c[1] + 5*c[2])*c[3] + 50*c[3]*c[3]) + 
	4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 
	99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 
	7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]))))/ 332640.;
		bf_mom[7] = (99*(a - b)*(7*(-60 + 15*c[0]*c[0] + 20*c[1]*c[1] - 
	20*c[0]*c[2] + 28*c[2]*c[2]) - 168*c[1]*c[3] + 204*c[3]*c[3]) - 33*(a - b)*
	 (21*(5*(a*a + a*b + b*b)*(-4 + c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 
	4*(2*a*a + a*b + 2*b*b)*c[1]*c[1]) - 84*((a*a + 3*a*b + b*b)*c[0] + (a - b)*
	(a + b)*c[1])* c[2] + 12*(15*a*a + 19*a*b + 15*b*b)*c[2]*c[2] + 18*(-12*a*
	b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1]+10*c[2]))*
	c[3] + 4*(50*a*a + 53*a*b + 50*b*b)*c[3]*c[3]) + 27720*(-a + b - ((a - b)*
	(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] - 27*c[2]*
	c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 49*c[2]*c[2] 
	- 42*c[1]*c[3] + 51*c[3]*c[3])))/840.) + (a - b)*(a*b*(9240 + 1155*c[0]*
	c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] + 57*c[2]*c[2] - 
	54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 
	858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*
	(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 
	45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 
	63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*
	c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])
	*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 
	45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 
	+ 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]
	*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])
	*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 166320.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + ((a 
	- b)*(a*b*(10 + 5*c[0] - 6*c[2]) + b*b*(10 + 5*c[0] + 5*c[1] - 2*c[2] - 
	3*c[3]) + a*a*(10 + 5*c[0] - 5*c[1] - 2*c[2] + 3*c[3])))/30. + (a - b + 
	((a - b)*(105*c[0]*c[0]*c[0] - 210*c[0]*c[0]*c[2] + 8*c[2]*(21*c[1]*c[1] 
	- 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) + 12*c[0]*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3])))/840.)/3. - ((a - b)*(a*b*
	(9240 + 1155*c[0]*c[0]*c[0] - 4158*c[0]*c[0]*c[2] + 132*c[0]*(21*c[1]*c[1] 
	+ 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 8*c[2]*(99*c[1]*c[1] + 
	517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) + b*b*(1155*c[0]*c[0]*c[0] 
	+ 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 
	34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] 
	- 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) + a*a*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 
	21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) 
	+ 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] 
	+ 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*
	c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3]))))/ 83160.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
case 5:
switch (interface_type)
	{
	case 0:
		for(i=0;i<npts;i++)	{bf_mom[i]=gauss_mom[i];}
		break;
	case 2:
		bf_mom[0] = -((-1 + a)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 
	+ 42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3]
	 + 583*c[3]*c[3] + c[2]*(693 - 682*c[4]) - 99*c[4] + 579*c[4]*c[4]) + 4*
	(-75075 + 13442*c[2]*c[2]*c[2] + 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] + 
	22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(21 + 6*c[2] + 22*c[4])
	 - 117*c[2]*c[2]*(-209 + 226*c[4]) - 78*c[1]*c[3]*(297 + 286*c[2] + 238*c[4])
	 + 6*c[2]*(4173*c[3]*c[3] + c[4]*(-4433 + 4465*c[4]))) + a*(-15015*c[0]*c[0]*
	c[0] + 1287*c[0]*c[0]* (35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4])-156*c[0]*
	(693*c[1]*c[1] + 363*c[2]*c[2] - 693*c[3] + 517*c[3]*c[3] - 33*c[1]*(-35 + 
	14*c[2] + 6*c[3] - 26*c[4]) - 11*c[2]*(21 + 90*c[3] - 10*c[4]) + 429*c[4] - 
	1078*c[3]*c[4] + 547*c[4]*c[4]) + 4*(-75075 + 18018*c[1]*c[1]*c[1] - 3718*
	c[2]*c[2]*c[2] + 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-63 
	+ 54*c[2] + 18*c[3] - 34*c[4]) - 42042*c[3]*c[4] + 12954*c[3]*c[3]*c[4] + 
	21333*c[4]*c[4] - 30822*c[3]*c[4]*c[4] + 8530*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(-363 + 154*c[3] + 118*c[4]) + 6*c[2]*(91*c[3]*c[3] + (715 - 987*c[4])*c[4] + 
	117*c[3]*(-55 + 14*c[4])) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*
	(21 + 38*c[3] - 10*c[4]) + c[4]*(429 + 547*c[4]) - c[3]*(99 + 602*c[4])))) + 
	2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(35 + 35*c[1] + 14*c[2] - 21*c[3]
	 + 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*
	c[3] - 33*c[1]*(-35 + 14*c[2] + 24*c[3] - 26*c[4]) + 330*c[4] - 1078*c[3]*c[4]
	 + 1126*c[4]*c[4] - 22*c[2]*(-21 + 45*c[3] + 26*c[4])) - 4*(-75075 + 9009*
	c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] 
	- 429*c[1]*c[1]*(-42 + 24*c[2] + 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] + 11292*
	c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4]
	 - 39*c[2]*c[2]*(-495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*
	c[3]*c[3] - 11*c[2]*(21 + 64*c[3] - 10*c[4]) - 12*c[3]*(33 + 70*c[4]) + c[4]*
	(429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 2*
	c[4]*(-1859 + 1739*c[4]))))))/8648640.;
		bf_mom[1] = -((-1 + a)*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-175 
	+ 140*c[1] - 98*c[2] - 84*c[3] - 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*
	c[2]*c[2] + 1386*c[3] + 2783*c[3]*c[3] + 11*c[2]*(147 + 180*c[3] - 166*c[4]) 
	+ 66*c[1]*(-35 + 14*c[2] - 33*c[3] - 26*c[4]) + 561*c[4] + 2156*c[3]*c[4] + 
	2831*c[4]*c[4]) + 4*(375375 + 36036*c[1]*c[1]*c[1] - 32890*c[2]*c[2]*c[2] - 
	108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(2607 + 308*c[3] - 
	2270*c[4]) + 429*c[1]*c[1]*(-189 + 90*c[2] - 36*c[3] - 134*c[4]) - 84084*c[3]
	*c[4] - 54798*c[3]*c[3]*c[4] - 110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] - 
	11258*c[4]*c[4]*c[4] + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + 22*c[2]*
	(-21 + 77*c[3] + 10*c[4]) + 2*c[4]*(429 + 547*c[4]) + c[3]*(1089 + 1918*c[4]))
	 - 6*c[2]*(12701*c[3]*c[3] - 234*c[3]*(-55 + 14*c[4]) + c[4]*(-11869 + 11421
	*c[4]))) + 2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]* (35 + 35*c[1] + 
	14*c[2] - 21*c[3] + 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*
	c[3] + 1100*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2] + 24*c[3] - 26*c[4]) + 330*c[4]
	 - 1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*(-21 + 45*c[3] + 26*c[4])) - 4*
	(-75075 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-42 + 24*c[2] + 9*c[3] - 28*c[4]) - 21021*
	c[3]*c[4] + 11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 
	3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-495 + 77*c[3] + 398*c[4]) + 39*c[1]*
	(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 64*c[3] - 10*c[4]) - 12*c[3]*
	(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*
	(-55 + 14*c[4]) + 2*c[4]*(-1859 + 1739*c[4])))) + a*(75075*c[0]*c[0]*c[0] - 
	1287*c[0]*c[0]*(175 + 70*c[1] + 154*c[2] - 42*c[3] + 2*c[4]) + 156*c[0]* 
	(1617*c[1]*c[1] + 2871*c[2]*c[2] - 693*c[3] + 2849*c[3]*c[3] - 33*c[1]*
	(-35 + 14*c[2] + 78*c[3] - 26*c[4]) + 33*c[4] - 1078*c[3]*c[4] + 2863*c[4]*
	c[4] - 11*c[2]*(-231 + 90*c[3] + 238*c[4])) - 4*(-375375 + 18018*c[1]*c[1]*
	c[1] + 50050*c[2]*c[2]*c[2] + 111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*
	c[1]*c[1]*(-147 + 30*c[2] + 18*c[3] - 122*c[4]) - 42042*c[3]*c[4] + 51474*
	c[3]*c[3]*c[4] + 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] + 794*c[4]*c[4]*c[4]
	 - 39*c[2]*c[2]*(-2871 + 154*c[3] + 2830*c[4]) + 78*c[1]*(363*c[2]*c[2]+517*
	c[3]*c[3] - 11*c[2]*(21 + 142*c[3] - 10*c[4]) - 3*c[3]*(429 + 518*c[4]) + 
	c[4]*(429 + 547*c[4])) + 6*c[2]*(16783*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) 
	+ c[4]*(-17017 + 16873*c[4]))))))/8648640.;
		bf_mom[2] = -((-1 + a)*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(175 +
	 140*c[1] - 98*c[2] - 84*c[3] - 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*
	c[2]*c[2] - 1386*c[3] + 2783*c[3]*c[3] + 11*c[2]*(-147 + 180*c[3] - 166*c[4])
	 + 66*c[1]*(35 + 14*c[2] - 33*c[3] - 26*c[4]) - 561*c[4] + 2156*c[3]*c[4] + 
	2831*c[4]*c[4]) + 4*(-75075 + 36036*c[1]*c[1]*c[1] - 32890*c[2]*c[2]*c[2] + 
	108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-2607 + 308*c[3] - 
	2270*c[4]) + 429*c[1]*c[1]*(189 + 90*c[2] - 36*c[3] - 134*c[4]) + 84084*c[3]*
	c[4] - 54798*c[3]*c[3]*c[4] + 110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] - 
	11258*c[4]*c[4]*c[4] + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + 22*c[2]*
	(21 + 77*c[3] + 10*c[4]) + 2*c[4]*(-429 + 547*c[4]) + c[3]*(-1089+1918*c[4]))
	 - 6*c[2]*(12701*c[3]*c[3] - 234*c[3]*(55 + 14*c[4]) + c[4]*(11869 + 11421*
	c[4]))) + 2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-35 + 35*c[1] + 14*
	c[2] - 21*c[3] + 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3]
	 + 1100*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 
	1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*(21 + 45*c[3] + 26*c[4])) - 4*
	(15015 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(42 + 24*c[2] + 9*c[3] - 28*c[4]) + 21021*
	c[3]*c[4] + 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 
	3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(495 + 77*c[3] + 398*c[4]) + 39*c[1]*
	(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(396 - 840*c[4]) + c[2]*(231 - 704*c[3]
	 + 110*c[4]) + c[4]*(-429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*
	(55 + 14*c[4]) + 2*c[4]*(1859 + 1739*c[4])))) + a*(75075*c[0]*c[0]*c[0] - 
	1287*c[0]*c[0]*(-175 + 70*c[1] + 154*c[2] - 42*c[3] + 2*c[4]) + 156*c[0]* 
	(1617*c[1]*c[1] + 2871*c[2]*c[2] + 693*c[3] + 2849*c[3]*c[3] - 33*c[1]*(35 + 
	14*c[2] + 78*c[3] - 26*c[4]) - 33*c[4] - 1078*c[3]*c[4] + 2863*c[4]*c[4] - 
	11*c[2]*(231 + 90*c[3] + 238*c[4])) - 4*(75075 + 18018*c[1]*c[1]*c[1] + 
	50050*c[2]*c[2]*c[2] - 111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*
	c[1]*(147 + 30*c[2] + 18*c[3] - 122*c[4]) + 42042*c[3]*c[4] + 51474*c[3]*
	c[3]*c[4] - 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] + 794*c[4]*c[4]*c[4] - 
	39*c[2]*c[2]*(2871 + 154*c[3] + 2830*c[4]) + 78*c[1]*(363*c[2]*c[2] + 517*
	c[3]*c[3] + c[2]*(231 - 1562*c[3] + 110*c[4]) - 3*c[3]*(-429 + 518*c[4]) + 
	c[4]*(-429 + 547*c[4])) + 6*c[2]*(16783*c[3]*c[3] + 117*c[3]*(55 + 14*c[4])
	 + c[4]*(17017 + 16873*c[4]))))))/8648640.;
		bf_mom[3] = -((-1 + a)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*
	(-35 + 42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*
	c[1]*c[3] + 583*c[3]*c[3] + 99*c[4] + 579*c[4]*c[4] - 11*c[2]*(63 + 62*c[4]))
	 + 4*(15015 + 13442*c[2]*c[2]*c[2] - 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] -
	 22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(-21 + 6*c[2]+22*c[4])
	 - 117*c[2]*c[2]*(209 + 226*c[4]) - 78*c[1]*c[3]*(-297 + 286*c[2] + 238*c[4])
	 + 6*c[2]*(4173*c[3]*c[3] + c[4]*(4433 + 4465*c[4])))-a*(15015*c[0]*c[0]*c[0]
	 - 1287*c[0]*c[0]*(-35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4]) + 156*c[0]* 
	(693*c[1]*c[1] + 363*c[2]*c[2] + 693*c[3] + 517*c[3]*c[3] - 33*c[1]*(35 + 
	14*c[2] + 6*c[3] - 26*c[4]) - 429*c[4] - 1078*c[3]*c[4] + 547*c[4]*c[4] + 
	c[2]*(231 - 990*c[3] + 110*c[4])) - 4*(15015 + 18018*c[1]*c[1]*c[1] - 3718*
	c[2]*c[2]*c[2] - 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*
	(63 + 54*c[2] + 18*c[3] - 34*c[4]) + 42042*c[3]*c[4] + 12954*c[3]*c[3]*c[4] - 
	21333*c[4]*c[4] - 30822*c[3]*c[4]*c[4] + 8530*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(363 + 154*c[3] + 118*c[4]) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*
	(99 - 602*c[4]) + c[2]*(231 - 418*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4])) 
	+ 6*c[2]*(91*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) - c[4]*(715 + 987*c[4])))) + 
	2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-35 + 35*c[1] + 14*c[2] - 
	21*c[3] + 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 
	1100*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 
	1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*(21 + 45*c[3] + 26*c[4])) - 4*
	(15015 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(42 + 24*c[2] + 9*c[3] - 28*c[4]) + 21021*
	c[3]*c[4] + 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 
	3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(495 + 77*c[3] + 398*c[4]) + 39*c[1]*
	(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(396 - 840*c[4]) + c[2]*(231 - 704*c[3]
	 + 110*c[4]) + c[4]*(-429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*
	(55 + 14*c[4]) + 2*c[4]*(1859 + 1739*c[4]))))))/8648640.;
		bf_mom[4] = ((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(70 +
	 35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]*
	c[2] - 693*c[3] + 2266*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2] + 60*c[3] - 26*c[4])
	 + 132*c[4] - 1078*c[3]*c[4] + 2284*c[4]*c[4] - 22*c[2]*(-84 + 45*c[3] + 88*
	c[4])) - 4*(-150150 + 9009*c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] + 44187*
	c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-63 + 18*c[2] + 9*c[3] - 
	50*c[4]) - 21021*c[3]*c[4] + 20922*c[3]*c[3]*c[4] + 44538*c[4]*c[4] - 15411*
	c[3]*c[4]*c[4] + 1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-1122 + 77*c[3] + 1076*
	c[4]) + 3*c[2]*(12610*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 88*c[4]*(-143 + 
	141*c[4])) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 116*c[3] 
	- 10*c[4]) + c[4]*(429 + 547*c[4]) - 2*c[3]*(495 + 658*c[4]))) + a*(15015*
	c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(35 + 35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) 
	+ 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*c[3] - 33*
	c[1]*(-35 + 14*c[2] + 24*c[3] - 26*c[4]) + 330*c[4] - 1078*c[3]*c[4] + 1126*
	c[4]*c[4] - 22*c[2]*(-21 + 45*c[3] + 26*c[4])) - 4*(-75075 + 9009*c[1]*c[1]*
	c[1] + 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*
	c[1]*c[1]*(-42 + 24*c[2] + 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] + 11292*c[3]*
	c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 
	39*c[2]*c[2]*(-495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*
	c[3] - 11*c[2]*(21 + 64*c[3] - 10*c[4]) - 12*c[3]*(33 + 70*c[4]) + c[4]*
	(429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 
	2*c[4]*(-1859 + 1739*c[4]))))))/2162160.;
		bf_mom[5] = ((-1 + a)*(75075*c[0]*c[0]*c[0] + 2574*c[0]*c[0]*(70*c[1]
	 - 49*c[2] - 42*c[3] - 17*c[4]) + 156*c[0]*(-5775 + 2079*c[1]*c[1] + 2607*
	c[2]*c[2] + 2783*c[3]*c[3] + 22*c[2]*(90*c[3] - 83*c[4]) + 66*c[1]*(14*c[2] 
	- 33*c[3] - 26*c[4]) + 2156*c[3]*c[4] + 2831*c[4]*c[4]) + 8*(-150150 + 18018*
	c[1]*c[1]*c[1] - 16445*c[2]*c[2]*c[2] + 54054*c[3] - 13806*c[3]*c[3]*c[3] + 
	429*c[1]*c[1]*(45*c[2] - 18*c[3] - 67*c[4]) + 21879*c[4] - 27399*c[3]*c[3]*
	c[4] - 30822*c[3]*c[4]*c[4] - 5629*c[4]*c[4]*c[4] + c[2]*c[2]*(-6006*c[3] + 
	44265*c[4]) + c[2]*(63063 - 38103*c[3]*c[3] + 9828*c[3]*c[4] - 34263*c[4]*
	c[4]) + 78*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 959*c[3]*c[4] + 
	547*c[4]*c[4] + 11*c[2]*(77*c[3] + 10*c[4]))) + 2*a*a*(15015*c[0]*c[0]*c[0] -
	 1287*c[0]*c[0]*(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 156*c[0]*(-1155 + 
	462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 12*c[3] - 
	13*c[4]) - 539*c[3]*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 26*c[4])) - 4*
	(60060 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*
	c[3]*c[3] - 429*c[1]*c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*c[4] + 11292*
	c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(77*c[3] + 398*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 
	22*c[2]*(32*c[3] - 5*c[4]) - 840*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-3003 
	+ 2132*c[3]*c[3] + 819*c[3]*c[4] + 1739*c[4]*c[4]))) + a*(75075*c[0]*c[0]*
	c[0] - 2574*c[0]*c[0]*(35*c[1] + 77*c[2] - 21*c[3] + c[4]) + 156*c[0]*(1617*
	c[1]*c[1] + 2871*c[2]*c[2] - 66*c[1]*(7*c[2] + 39*c[3] - 13*c[4]) - 22*c[2]*
	(45*c[3] + 119*c[4]) + 7*(-825 + 407*c[3]*c[3] - 154*c[3]*c[4] + 409*c[4]*
	c[4])) - 8*(150150 + 9009*c[1]*c[1]*c[1] + 25025*c[2]*c[2]*c[2] + 27027*c[3] 
	- 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(15*c[2] + 9*c[3] - 61*c[4]) - 1287*c[4]
	 + 25737*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 397*c[4]*c[4]*c[4] - 39*c[2]*
	c[2]*(77*c[3] + 1415*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 
	22*c[2]*(71*c[3] - 5*c[4]) - 1554*c[3]*c[4] + 547*c[4]*c[4]) + c[2]*(-99099 +
	 50349*c[3]*c[3] + 4914*c[3]*c[4] + 50619*c[4]*c[4])))))/4324320.;
		bf_mom[6] = ((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-70 
	+ 35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]
	*c[2] + 693*c[3] + 2266*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 60*c[3] - 26*c[4])
	 - 132*c[4] - 1078*c[3]*c[4] + 2284*c[4]*c[4] - 22*c[2]*(84+45*c[3]+88*c[4]))
	 - 4*(30030 + 9009*c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] - 44187*c[3]*c[3] -
	 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(63 + 18*c[2] + 9*c[3] - 50*c[4]) + 
	21021*c[3]*c[4] + 20922*c[3]*c[3]*c[4] - 44538*c[4]*c[4] - 15411*c[3]*c[4]
	*c[4] + 1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(1122 + 77*c[3] + 1076*c[4]) + 
	3*c[2]*(12610*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + 88*c[4]*(143 + 141*c[4]))
	 + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(990 - 1316*c[4]) + c[2]*
	(231 - 1276*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4]))) + a*(15015*c[0]*c[0]
	*c[0] - 1287*c[0]*c[0]*(-35 + 35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 78*
	c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 1100*c[3]*c[3] - 33*c[1]*
	(35 + 14*c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 1078*c[3]*c[4] + 1126*c[4]*
	c[4] - 22*c[2]*(21 + 45*c[3] + 26*c[4])) - 4*(15015 + 9009*c[1]*c[1]*c[1] + 
	4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*
	(42 + 24*c[2] + 9*c[3] - 28*c[4]) + 21021*c[3]*c[4] + 11292*c[3]*c[3]*c[4] - 
	21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*
	(396 - 840*c[4]) + c[2]*(231 - 704*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4]))
	 + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + 
	2*c[4]*(1859 + 1739*c[4]))))))/2162160.;
		bf_mom[7] = ((-1 + a)*(-15015*c[0]*c[0]*c[0] + 7722*c[0]*c[0]*
	(7*c[2] - c[4]) - 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*
	c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(30030 + 6721*c[2]*
	c[2]*c[2] - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*c[1]*c[3] + 
	4815*c[3]*c[3])*c[4] - 967*c[4]*c[4]*c[4] + 3*c[2]*(-9009 + 429*c[1]*c[1] - 
	3718*c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])) + a*(-15015*c[0]*c[0]*c[0]
	 + 2574*c[0]*c[0]*(35*c[1] - 7*c[2] - 21*c[3] + 13*c[4]) - 156*c[0]*(-1155 +
	 693*c[1]*c[1] + 363*c[2]*c[2] + 517*c[3]*c[3] - 66*c[1]*(7*c[2] + 3*c[3] -
	 13*c[4]) - 110*c[2]*(9*c[3] - c[4]) - 1078*c[3]*c[4] + 547*c[4]*c[4]) + 
	8*(30030 + 9009*c[1]*c[1]*c[1] - 1859*c[2]*c[2]*c[2] + 27027*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(27*c[2] + 9*c[3] - 17*c[4]) - 16731*c[4] + 
	6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 4265*c[4]*c[4]*c[4] - 39*c[2]*
	c[2]*(77*c[3] + 59*c[4]) + 21*c[2]*(429 + 13*c[3]*c[3] + 234*c[3]*c[4] - 
	141*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*
	(19*c[3] - 5*c[4]) - 602*c[3]*c[4] + 547*c[4]*c[4]))) + 2*a*a*(15015*c[0]*
	c[0]*c[0] - 1287*c[0]*c[0]*(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4])+156*c[0]*
	(-1155 + 462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 
	12*c[3] - 13*c[4]) - 539*c[3]*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 
	26*c[4])) - 4*(60060 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3]
	 - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*
	c[4] + 11292*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 
	39*c[2]*c[2]*(77*c[3] + 398*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*
	c[3]*c[3] - 22*c[2]*(32*c[3] - 5*c[4]) - 840*c[3]*c[4] + 547*c[4]*c[4]) + 
	6*c[2]*(-3003 + 2132*c[3]*c[3] + 819*c[3]*c[4] + 
	1739*c[4]*c[4])))))/4324320.;
		bf_mom[8] = -((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 156*c[0]*(-2310 + 693*c[1]*c[1] + 
	1122*c[2]*c[2] + 1133*c[3]*c[3] - 33*c[1]*(7*c[2] + 30*c[3] - 13*c[4]) - 539*
	c[3]*c[4] + 1142*c[4]*c[4] - 11*c[2]*(45*c[3] + 88*c[4])) - 4*(120120 + 9009*
	c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 
	429*c[1]*c[1]*(18*c[2] + 9*c[3] - 50*c[4]) - 5148*c[4] + 20922*c[3]*c[3]*c[4]
	 - 15411*c[3]*c[4]*c[4] + 1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 1076*
	c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(58*c[3] - 
	5*c[4]) - 1316*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-12012 + 6305*c[3]*c[3] + 
	819*c[3]*c[4] + 6204*c[4]*c[4])) + a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 156*c[0]*(-1155 + 462*c[1]*c[1] + 
	495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 12*c[3] - 13*c[4]) - 539*
	c[3]*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 26*c[4])) - 4*(60060 + 9009*
	c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 
	429*c[1]*c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*c[4] + 11292*c[3]*c[3]*c[4]
	 - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 
	398*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(32*c[3]
	 - 5*c[4]) - 840*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-3003 + 2132*c[3]*c[3]
	 + 819*c[3]*c[4] + 1739*c[4]*c[4])))))/1081080.;
		break;
	case 3:
		bf_mom[0] = ((-1 + a)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] + c[2]*(693 - 682*c[4]) - 99*c[4] + 579*c[4]*c[4]) + 4*(-15015 +
	 13442*c[2]*c[2]*c[2] + 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4]+22581*c[4]*c[4] 
	- 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(21 + 6*c[2] + 22*c[4]) - 117*c[2]*c[2]*
	(-209 + 226*c[4]) - 78*c[1]*c[3]*(297 + 286*c[2] + 238*c[4]) + 6*c[2]*(4173*
	c[3]*c[3] + c[4]*(-4433 + 4465*c[4]))) + a*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*
	c[0]*(35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4]) - 156*c[0]*(693*c[1]*c[1] +
	 363*c[2]*c[2] - 693*c[3] + 517*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2] + 6*c[3] - 
	26*c[4]) - 11*c[2]*(21 + 90*c[3] - 10*c[4]) + 429*c[4] - 1078*c[3]*c[4] + 
	547*c[4]*c[4]) + 4*(-15015 + 18018*c[1]*c[1]*c[1] - 3718*c[2]*c[2]*c[2] + 
	20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-63 + 54*c[2] + 
	18*c[3] - 34*c[4]) - 42042*c[3]*c[4] + 12954*c[3]*c[3]*c[4] + 21333*c[4]*c[4] 
	- 30822*c[3]*c[4]*c[4] + 8530*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-363 + 154*c[3] 
	+ 118*c[4]) + 6*c[2]*(91*c[3]*c[3] + (715 - 987*c[4])*c[4] + 117*c[3]*(-55 + 
	14*c[4])) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 38*c[3] - 
	10*c[4]) + c[4]*(429 + 547*c[4]) - c[3]*(99 + 602*c[4])))) + 2*a*a*(15015*c[0]
	*c[0]*c[0] - 1287*c[0]*c[0]*(35 + 35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 
	78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*c[3] - 33*c[1]*
	(-35 + 14*c[2] + 24*c[3] - 26*c[4]) + 330*c[4] - 1078*c[3]*c[4] + 1126*c[4]*
	c[4] - 22*c[2]*(-21 + 45*c[3] + 26*c[4])) - 4*(-15015 + 9009*c[1]*c[1]*c[1] + 
	4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*
	(-42 + 24*c[2] + 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] + 11292*c[3]*c[3]*c[4] + 
	21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(-495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]
	*(21 + 64*c[3] - 10*c[4]) - 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) + 
	3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*
	(-1859 + 1739*c[4]))))))/8648640.;
		bf_mom[1] = ((-1 + a)*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-175 + 
	140*c[1] - 98*c[2] - 84*c[3] - 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*c[2]
	*c[2] + 1386*c[3] + 2783*c[3]*c[3] + 11*c[2]*(147 + 180*c[3] - 166*c[4]) + 
	66*c[1]*(-35 + 14*c[2] - 33*c[3] - 26*c[4]) + 561*c[4] + 2156*c[3]*c[4] + 
	2831*c[4]*c[4]) + 4*(75075 + 36036*c[1]*c[1]*c[1] - 32890*c[2]*c[2]*c[2] - 
	108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(2607 + 308*c[3] - 2270
	*c[4]) + 429*c[1]*c[1]*(-189 + 90*c[2] - 36*c[3] - 134*c[4]) - 84084*c[3]*c[4]
	 - 54798*c[3]*c[3]*c[4] - 110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] - 11258*c[4]
	*c[4]*c[4] + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + 22*c[2]*(-21 + 77*c[3] 
	+ 10*c[4]) + 2*c[4]*(429 + 547*c[4]) + c[3]*(1089 + 1918*c[4])) - 6*c[2]*
	(12701*c[3]*c[3] - 234*c[3]*(-55 + 14*c[4]) + c[4]*(-11869 + 11421*c[4]))) + 
	2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(35 + 35*c[1] + 14*c[2] - 21*c[3]
	 + 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*
	c[3] - 33*c[1]*(-35 + 14*c[2] + 24*c[3] - 26*c[4]) + 330*c[4] - 1078*c[3]*c[4]
	 + 1126*c[4]*c[4] - 22*c[2]*(-21 + 45*c[3] + 26*c[4])) - 4*(-15015 + 9009*c[1]
	*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 
	429*c[1]*c[1]*(-42 + 24*c[2] + 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] + 11292*
	c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] 
	- 39*c[2]*c[2]*(-495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]
	*c[3] - 11*c[2]*(21 + 64*c[3] - 10*c[4]) - 12*c[3]*(33 + 70*c[4]) + c[4]*(429 
	+ 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*
	(-1859 + 1739*c[4])))) + a*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(175 + 70*
	c[1] + 154*c[2] - 42*c[3] + 2*c[4]) + 156*c[0]*(1617*c[1]*c[1]+2871*c[2]*c[2] 
	- 693*c[3] + 2849*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2] + 78*c[3] - 26*c[4]) + 
	33*c[4] - 1078*c[3]*c[4] + 2863*c[4]*c[4] - 11*c[2]*(-231 + 90*c[3] + 
	238*c[4])) - 4*(-75075 + 18018*c[1]*c[1]*c[1] + 50050*c[2]*c[2]*c[2] + 111111*
	c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-147 + 30*c[2] + 18*c[3] - 
	122*c[4]) - 42042*c[3]*c[4] + 51474*c[3]*c[3]*c[4] + 111657*c[4]*c[4] - 30822*
	c[3]*c[4]*c[4] + 794*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-2871 + 154*c[3] + 2830*
	c[4]) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 142*c[3] - 
	10*c[4]) - 3*c[3]*(429 + 518*c[4]) + c[4]*(429 + 547*c[4])) + 6*c[2]*(16783*
	c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + c[4]*
	(-17017 + 16873*c[4]))))))/8648640.;
		bf_mom[2] = ((-1 + a)*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(175 + 
	140*c[1] - 98*c[2] - 84*c[3] - 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*c[2]
	*c[2] - 1386*c[3] + 2783*c[3]*c[3] + 11*c[2]*(-147 + 180*c[3] - 166*c[4]) + 
	66*c[1]*(35 + 14*c[2] - 33*c[3] - 26*c[4]) - 561*c[4] + 2156*c[3]*c[4] + 
	2831*c[4]*c[4]) + 4*(-375375 + 36036*c[1]*c[1]*c[1] - 32890*c[2]*c[2]*c[2] + 
	108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-2607 + 308*c[3] - 
	2270*c[4]) + 429*c[1]*c[1]*(189 + 90*c[2] - 36*c[3] - 134*c[4]) + 84084*c[3]*
	c[4] - 54798*c[3]*c[3]*c[4] + 110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] - 
	11258*c[4]*c[4]*c[4] + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + 22*c[2]*
	(21 + 77*c[3] + 10*c[4]) + 2*c[4]*(-429 + 547*c[4]) + c[3]*(-1089+1918*c[4])) 
	- 6*c[2]*(12701*c[3]*c[3] - 234*c[3]*(55 + 14*c[4])+c[4]*(11869+11421*c[4])))
	 + 2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-35 + 35*c[1] + 14*c[2] - 
	21*c[3] + 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 
	1100*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 
	1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*(21 + 45*c[3] + 26*c[4])) - 4*
	(75075 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(42 + 24*c[2] + 9*c[3] - 28*c[4]) + 21021*c[3]*
	c[4] + 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*
	c[4]*c[4]*c[4] - 39*c[2]*c[2]*(495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*
	c[2] + 517*c[3]*c[3] + c[3]*(396 - 840*c[4]) + c[2]*(231 - 704*c[3]+110*c[4]) 
	+ c[4]*(-429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) 
	+ 2*c[4]*(1859 + 1739*c[4])))) + a*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(-175 + 70*c[1] + 154*c[2] - 42*c[3] + 2*c[4]) + 156*c[0]* (1617*c[1]*c[1] + 
	2871*c[2]*c[2] + 693*c[3] + 2849*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 78*c[3] 
	- 26*c[4]) - 33*c[4] - 1078*c[3]*c[4] + 2863*c[4]*c[4] - 11*c[2]*(231 +90*c[3]
	 + 238*c[4])) - 4*(375375 + 18018*c[1]*c[1]*c[1] + 50050*c[2]*c[2]*c[2] - 
	111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(147 + 30*c[2] + 18*
	c[3] - 122*c[4]) + 42042*c[3]*c[4] + 51474*c[3]*c[3]*c[4] - 111657*c[4]*c[4] 
	- 30822*c[3]*c[4]*c[4] + 794*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(2871 + 154*c[3] + 
	2830*c[4]) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[2]*(231 - 1562*c[3] + 
	110*c[4]) - 3*c[3]*(-429 + 518*c[4]) + c[4]*(-429 + 547*c[4])) + 6*c[2]*
	(16783*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + c[4]*
	(17017 + 16873*c[4]))))))/8648640.;
		bf_mom[3] = ((-1 + a)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] + 99*c[4] + 579*c[4]*c[4] - 11*c[2]*(63 + 62*c[4])) + 4*
	(75075 + 13442*c[2]*c[2]*c[2] - 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] - 22581*
	c[4]*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(-21 + 6*c[2] + 22*c[4]) - 
	117*c[2]*c[2]*(209 + 226*c[4]) - 78*c[1]*c[3]*(-297 + 286*c[2] + 238*c[4]) + 
	6*c[2]*(4173*c[3]*c[3] + c[4]*(4433 + 4465*c[4]))) - a*(15015*c[0]*c[0]*c[0] -
	 1287*c[0]*c[0]*(-35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4]) + 156*c[0]* 
	(693*c[1]*c[1] + 363*c[2]*c[2] + 693*c[3] + 517*c[3]*c[3] - 33*c[1]*(35 + 14*
	c[2] + 6*c[3] - 26*c[4]) - 429*c[4] - 1078*c[3]*c[4] + 547*c[4]*c[4] + c[2]*
	(231 - 990*c[3] + 110*c[4])) - 4*(75075 + 18018*c[1]*c[1]*c[1] - 3718*c[2]*
	c[2]*c[2] - 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(63 + 54*
	c[2] + 18*c[3] - 34*c[4]) + 42042*c[3]*c[4] + 12954*c[3]*c[3]*c[4] - 21333*
	c[4]*c[4] - 30822*c[3]*c[4]*c[4] + 8530*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(363 + 
	154*c[3] + 118*c[4]) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(99 - 
	602*c[4]) + c[2]*(231 - 418*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4])) + 6*
	c[2]*(91*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) - c[4]*(715 + 987*c[4])))) + 
	2*a*a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-35 + 35*c[1] + 14*c[2] - 
	21*c[3] + 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 
	1100*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 
	1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*(21 + 45*c[3] + 26*c[4])) - 4*
	(75075 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*
	c[3]*c[3]*c[3] - 429*c[1]*c[1]*(42 + 24*c[2] + 9*c[3] - 28*c[4]) + 21021*c[3]*
	c[4] + 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*
	c[4]*c[4]*c[4] - 39*c[2]*c[2]*(495 + 77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*
	c[2] + 517*c[3]*c[3] + c[3]*(396 - 840*c[4]) + c[2]*(231 - 704*c[3] +110*c[4])
	 + c[4]*(-429 + 547*c[4])) + 3*c[2]*(4264*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) 
	+ 2*c[4]*(1859 + 1739*c[4]))))))/8648640.;
		bf_mom[4] = -((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(70 
	+ 35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]*
	c[2] - 693*c[3] + 2266*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2] + 60*c[3] - 26*c[4])
	 + 132*c[4] - 1078*c[3]*c[4] + 2284*c[4]*c[4] - 22*c[2]*(-84 + 45*c[3] + 88*
	c[4])) - 4*(-30030 + 9009*c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] + 44187*c[3]*
	c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-63 + 18*c[2] + 9*c[3] - 50*c[4]) 
	- 21021*c[3]*c[4] + 20922*c[3]*c[3]*c[4] + 44538*c[4]*c[4] - 15411*c[3]*c[4]*
	c[4] + 1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-1122 + 77*c[3] + 1076*c[4]) + 3*
	c[2]*(12610*c[3]*c[3] + 117*c[3]*(-55 + 14*c[4]) + 88*c[4]*(-143 + 141*c[4])) 
	+ 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 116*c[3] - 10*c[4]) 
	+ c[4]*(429 + 547*c[4]) - 2*c[3]*(495 + 658*c[4]))) + a*(15015*c[0]*c[0]*c[0] 
	- 1287*c[0]*c[0]*(35 + 35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 78*c[0]*(924*
	c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*c[3] - 33*c[1]*(-35 + 14*c[2]
	 + 24*c[3] - 26*c[4]) + 330*c[4] - 1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*c[2]*
	(-21 + 45*c[3] + 26*c[4])) - 4*(-15015 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]
	*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-42 + 24*c[2] 
	+ 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] + 11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4]
	 - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-495 + 77*c[3] +
	 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] - 11*c[2]*(21 + 64*c[3] - 
	10*c[4]) - 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) + 3*c[2]*(4264*c[3]
	*c[3] + 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*(-1859 + 1739*c[4]))))))/2162160.;
		bf_mom[5] = -((-1 + a)*(75075*c[0]*c[0]*c[0] + 2574*c[0]*c[0]*(70*c[1]
	 - 49*c[2] - 42*c[3] - 17*c[4]) + 156*c[0]*(-5775 + 2079*c[1]*c[1] + 2607*c[2]
	*c[2] + 2783*c[3]*c[3] + 22*c[2]*(90*c[3] - 83*c[4]) + 66*c[1]*(14*c[2] - 33*
	c[3] - 26*c[4]) + 2156*c[3]*c[4] + 2831*c[4]*c[4]) + 8*(150150 + 18018*c[1]*
	c[1]*c[1] - 16445*c[2]*c[2]*c[2] + 54054*c[3] - 13806*c[3]*c[3]*c[3] + 429*
	c[1]*c[1]*(45*c[2] - 18*c[3] - 67*c[4]) + 21879*c[4] - 27399*c[3]*c[3]*c[4] - 
	30822*c[3]*c[4]*c[4] - 5629*c[4]*c[4]*c[4] + c[2]*c[2]*(-6006*c[3] + 44265*
	c[4]) + c[2]*(63063 - 38103*c[3]*c[3] + 9828*c[3]*c[4] - 34263*c[4]*c[4]) + 
	78*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 959*c[3]*c[4] + 547*c[4]*c[4]
	 + 11*c[2]*(77*c[3] + 10*c[4]))) + 2*a*a*(15015*c[0]*c[0]*c[0] -1287*c[0]*c[0]
	*(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 156*c[0]*(-1155 + 462*c[1]*c[1] + 
	495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 12*c[3] - 13*c[4]) -539*c[3]
	*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 26*c[4])) - 4*(-60060 + 9009*c[1]*
	c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*
	c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*c[4] + 11292*c[3]*c[3]*c[4] - 15411*
	c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 398*c[4]) + 39*
	c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(32*c[3] - 5*c[4]) - 
	840*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-3003 + 2132*c[3]*c[3] + 819*c[3]*
	c[4] + 1739*c[4]*c[4]))) + a*(75075*c[0]*c[0]*c[0] - 2574*c[0]*c[0]*(35*c[1] +
	 77*c[2] - 21*c[3] + c[4]) + 156*c[0]*(1617*c[1]*c[1] + 2871*c[2]*c[2] - 66*
	c[1]*(7*c[2] + 39*c[3] - 13*c[4]) - 22*c[2]*(45*c[3] + 119*c[4]) + 7*(-825 + 
	407*c[3]*c[3] - 154*c[3]*c[4] + 409*c[4]*c[4])) - 8*(-150150 + 9009*c[1]*c[1]*
	c[1] + 25025*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]
	*(15*c[2] + 9*c[3] - 61*c[4]) - 1287*c[4] + 25737*c[3]*c[3]*c[4] - 15411*c[3]*
	c[4]*c[4] + 397*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 1415*c[4]) + 39*c[1]*
	(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(71*c[3] - 5*c[4]) - 1554*
	c[3]*c[4] + 547*c[4]*c[4]) + c[2]*(-99099 + 50349*c[3]*c[3] + 4914*c[3]*c[4]
	 + 50619*c[4]*c[4])))))/4324320.;
		bf_mom[6] = -((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-70 
	+ 35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]*
	c[2] + 693*c[3] + 2266*c[3]*c[3] - 33*c[1]*(35 + 14*c[2] + 60*c[3] - 26*c[4]) 
	- 132*c[4] - 1078*c[3]*c[4] + 2284*c[4]*c[4] - 22*c[2]*(84 + 45*c[3]+88*c[4]))
	 - 4*(150150 + 9009*c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] - 44187*c[3]*c[3] - 
	6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(63 + 18*c[2] + 9*c[3] - 50*c[4]) + 21021*
	c[3]*c[4] + 20922*c[3]*c[3]*c[4] - 44538*c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 
	1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(1122 + 77*c[3] + 1076*c[4]) + 3*c[2]*
	(12610*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + 88*c[4]*(143 + 141*c[4])) + 39*
	c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(990 - 1316*c[4]) + c[2]*(231 - 
	1276*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4]))) + a*(15015*c[0]*c[0]*c[0] - 
	1287*c[0]*c[0]*(-35 + 35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 78*c[0]* 
	(924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 1100*c[3]*c[3] - 33*c[1]*(35 + 14*
	c[2] + 24*c[3] - 26*c[4]) - 330*c[4] - 1078*c[3]*c[4] + 1126*c[4]*c[4] - 22*
	c[2]*(21 + 45*c[3] + 26*c[4])) - 4*(75075 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*
	c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(42 + 24*
	c[2] + 9*c[3] - 28*c[4]) + 21021*c[3]*c[4] + 11292*c[3]*c[3]*c[4] - 21957*
	c[4]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(495 + 
	77*c[3] + 398*c[4]) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + c[3]*(396 - 
	840*c[4]) + c[2]*(231 - 704*c[3] + 110*c[4]) + c[4]*(-429 + 547*c[4])) + 3*
	c[2]*(4264*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + 2*c[4]*
	(1859 + 1739*c[4]))))))/2162160.;
		bf_mom[7] = -((-1 + a)*(-15015*c[0]*c[0]*c[0] + 7722*c[0]*c[0]*
	(7*c[2] - c[4]) - 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*
	c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-30030 + 6721*c[2]*
	c[2]*c[2] - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*c[1]*c[3] + 
	4815*c[3]*c[3])*c[4] - 967*c[4]*c[4]*c[4] + 3*c[2]*(-9009 + 429*c[1]*c[1] - 
	3718*c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])) - a*(15015*c[0]*c[0]*c[0] 
	- 2574*c[0]*c[0]*(35*c[1] - 7*c[2] - 21*c[3] + 13*c[4]) + 156*c[0]*(-1155 + 
	693*c[1]*c[1] + 363*c[2]*c[2] + 517*c[3]*c[3] - 66*c[1]*(7*c[2] + 3*c[3] - 
	13*c[4]) - 110*c[2]*(9*c[3] - c[4]) - 1078*c[3]*c[4] + 547*c[4]*c[4]) - 8*
	(-30030 + 9009*c[1]*c[1]*c[1] - 1859*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*
	c[3]*c[3] - 429*c[1]*c[1]*(27*c[2] + 9*c[3] - 17*c[4]) - 16731*c[4] + 6477*
	c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 4265*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(77*c[3] + 59*c[4]) + 21*c[2]*(429 + 13*c[3]*c[3] + 234*c[3]*c[4] - 141*c[4]*
	c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(19*c[3] - 
	5*c[4]) - 602*c[3]*c[4] + 547*c[4]*c[4]))) + 2*a*a*(15015*c[0]*c[0]*c[0] - 
	1287*c[0]*c[0]*(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 156*c[0]*(-1155 + 
	462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 12*c[3] - 
	13*c[4]) - 539*c[3]*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 26*c[4])) - 4*
	(-60060 + 9009*c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*
	c[3]*c[3] - 429*c[1]*c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*c[4] + 11292*
	c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(77*c[3] + 398*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*
	c[2]*(32*c[3] - 5*c[4]) - 840*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-3003 + 
	2132*c[3]*c[3] + 819*c[3]*c[4] + 1739*c[4]*c[4])))))/4324320.;
		bf_mom[8] = ((a-1)*(a-1)*(30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] + 56*c[2] - 21*c[3] + 4*c[4]) + 156*c[0]*(-2310 + 693*c[1]*c[1] + 
	1122*c[2]*c[2] + 1133*c[3]*c[3] - 33*c[1]*(7*c[2] + 30*c[3] - 13*c[4]) - 539*
	c[3]*c[4] + 1142*c[4]*c[4] - 11*c[2]*(45*c[3] + 88*c[4])) - 4*(-120120 + 9009*
	c[1]*c[1]*c[1] + 18304*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 
	429*c[1]*c[1]*(18*c[2] + 9*c[3] - 50*c[4]) - 5148*c[4] + 20922*c[3]*c[3]*c[4] 
	- 15411*c[3]*c[4]*c[4] + 1364*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 1076*
	c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(58*c[3] - 
	5*c[4]) - 1316*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-12012 + 6305*c[3]*c[3] + 
	819*c[3]*c[4] + 6204*c[4]*c[4])) + a*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] + 14*c[2] - 21*c[3] + 10*c[4]) + 156*c[0]*(-1155 + 462*c[1]*c[1] + 
	495*c[2]*c[2] + 550*c[3]*c[3] - 33*c[1]*(7*c[2] + 12*c[3] - 13*c[4]) - 539*
	c[3]*c[4] + 563*c[4]*c[4] - 11*c[2]*(45*c[3] + 26*c[4])) - 4*(-60060 + 9009*
	c[1]*c[1]*c[1] + 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 
	429*c[1]*c[1]*(24*c[2] + 9*c[3] - 28*c[4]) - 12870*c[4] + 11292*c[3]*c[3]*c[4]
	 - 15411*c[3]*c[4]*c[4] + 3298*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(77*c[3] + 398*
	c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(32*c[3] - 
	5*c[4]) - 840*c[3]*c[4] + 547*c[4]*c[4]) + 6*c[2]*(-3003 + 2132*c[3]*c[3] + 
	819*c[3]*c[4] + 1739*c[4]*c[4])))))/1081080.;
		break;
	case 4:
		bf_mom[0] = -((1 + b)*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(175 + 
	140*c[1] + 98*c[2] - 84*c[3] + 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*
	c[2]*c[2] - 1386*c[3] + 2783*c[3]*c[3] - 66*c[1]*(-35 + 14*c[2] + 33*c[3] - 
	26*c[4]) + 561*c[4] - 2156*c[3]*c[4] + 2831*c[4]*c[4] - 11*c[2]*(-147 + 180*
	c[3] + 166*c[4])) - 4*(-75075 + 36036*c[1]*c[1]*c[1] + 32890*c[2]*c[2]*c[2] + 
	108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(-189 + 90*c[2] + 
	36*c[3] - 134*c[4]) - 84084*c[3]*c[4] + 54798*c[3]*c[3]*c[4]+110409*c[4]*c[4]
	 - 61644*c[3]*c[4]*c[4] + 11258*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(-2607+308*c[3]
	 + 2270*c[4]) + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] - 22*c[2]*(21+77*c[3]
	 - 10*c[4]) + 2*c[4]*(429 + 547*c[4]) - c[3]*(1089 + 1918*c[4])) + 6*c[2]*
	(12701*c[3]*c[3] + 234*c[3]*(-55 + 14*c[4]) + c[4]*(-11869 + 11421*c[4]))) + 
	2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 35*c[1] - 14*c[2] - 
	21*c[3] - 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 
	1100*c[3]*c[3] + 33*c[1]*(-35 + 14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(21 + 
	45*c[3] - 26*c[4]) + 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(15015 + 
	9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]
	*c[3] - 39*c[2]*c[2]*(495 + 77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(-42 + 24*c[2]
	 - 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] - 11292*c[3]*c[3]*c[4]-21957*c[4]*c[4]
	 - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 
	517*c[3]*c[3] + 11*c[2]*(-21 + 64*c[3] + 10*c[4]) + 12*c[3]*(33 + 70*c[4]) + 
	c[4]*(429 + 547*c[4])) - 3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 
	2*c[4]*(-1859 + 1739*c[4])))) + b*(-75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(-175 + 70*c[1] - 154*c[2] - 42*c[3] - 2*c[4]) - 156*c[0]* (1617*c[1]*c[1] + 
	2871*c[2]*c[2] + 693*c[3] + 2849*c[3]*c[3] + 11*c[2]*(231 + 90*c[3]-238*c[4])
	 + 33*c[1]*(-35 + 14*c[2] - 78*c[3] - 26*c[4]) + 33*c[4] + 1078*c[3]*c[4] + 
	2863*c[4]*c[4]) - 4*(75075 + 18018*c[1]*c[1]*c[1] - 50050*c[2]*c[2]*c[2] - 
	111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(2871 + 154*c[3] - 
	2830*c[4]) + 429*c[1]*c[1]*(-147 + 30*c[2] - 18*c[3] - 122*c[4]) - 42042*c[3]
	*c[4] - 51474*c[3]*c[3]*c[4] - 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] - 794*
	c[4]*c[4]*c[4] + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(-21 + 142
	*c[3] + 10*c[4]) + 3*c[3]*(429 + 518*c[4]) + c[4]*(429 + 547*c[4])) - 6*c[2]*
	(16783*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + c[4]*
	(-17017 + 16873*c[4]))))))/8648640.;
		bf_mom[1] = -((1 + b)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] + c[2]*(693 - 682*c[4]) - 99*c[4] + 579*c[4]*c[4]) + 4*(-15015 
	+ 13442*c[2]*c[2]*c[2] + 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] + 22581*c[4]*
	c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(21 + 6*c[2] + 22*c[4]) - 117*c[2]
	*c[2]*(-209 + 226*c[4]) - 78*c[1]*c[3]*(297 + 286*c[2] + 238*c[4]) + 6*c[2]*
	(4173*c[3]*c[3] + c[4]*(-4433 + 4465*c[4]))) + b*(15015*c[0]*c[0]*c[0] + 
	1287*c[0]*c[0]*(-35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4]) + 156*c[0]* 
	(693*c[1]*c[1] + 363*c[2]*c[2] + 693*c[3] + 517*c[3]*c[3] + 33*c[1]*(-35 + 
	14*c[2] - 6*c[3] - 26*c[4]) + 429*c[4] + 1078*c[3]*c[4] + 547*c[4]*c[4] + 
	11*c[2]*(-21 + 90*c[3] + 10*c[4])) + 4*(15015 + 18018*c[1]*c[1]*c[1] + 3718
	*c[2]*c[2]*c[2] - 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*
	(363 + 154*c[3] - 118*c[4]) + 429*c[1]*c[1]*(-63 + 54*c[2] - 18*c[3]-34*c[4])
	 - 42042*c[3]*c[4] - 12954*c[3]*c[3]*c[4] - 21333*c[4]*c[4] - 30822*c[3]*c[4]
	*c[4] - 8530*c[4]*c[4]*c[4] - 6*c[2]*(91*c[3]*c[3] + (715 - 987*c[4])*c[4] - 
	117*c[3]*(-55 + 14*c[4])) + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*
	(-21 + 38*c[3] + 10*c[4]) + c[4]*(429 + 547*c[4]) + c[3]*(99 + 602*c[4])))) + 
	2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 35*c[1] - 14*c[2]-21*c[3]
	 - 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 1100*c[3]
	*c[3] + 33*c[1]*(-35 + 14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(21 + 45*c[3] - 
	26*c[4]) + 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(15015 + 9009*c[1]
	*c[1]*c[1] - 4862*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 
	39*c[2]*c[2]*(495 + 77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(-42 + 24*c[2] - 
	9*c[3] - 28*c[4]) - 21021*c[3]*c[4] - 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] 
	- 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517
	*c[3]*c[3] + 11*c[2]*(-21 + 64*c[3] + 10*c[4]) + 12*c[3]*(33 + 70*c[4]) + 
	c[4]*(429 + 547*c[4])) - 3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 
	2*c[4]*(-1859 + 1739*c[4]))))))/8648640.;
		bf_mom[2] = -((1 + b)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] + 99*c[4] + 579*c[4]*c[4] - 11*c[2]*(63 + 62*c[4])) + 4*(75075 
	+ 13442*c[2]*c[2]*c[2] - 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] - 22581*c[4]
	*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(-21 + 6*c[2] + 22*c[4]) - 117
	*c[2]*c[2]*(209 + 226*c[4]) - 78*c[1]*c[3]*(-297 + 286*c[2] + 238*c[4]) + 
	6*c[2]*(4173*c[3]*c[3] + c[4]*(4433 + 4465*c[4]))) + b*(15015*c[0]*c[0]*c[0] 
	+ 1287*c[0]*c[0]* (35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4]) + 156*c[0]*
	(693*c[1]*c[1] + 363*c[2]*c[2] - 693*c[3] + 517*c[3]*c[3] + 33*c[1]*(35 + 
	14*c[2] - 6*c[3] - 26*c[4]) - 429*c[4] + 1078*c[3]*c[4] + 547*c[4]*c[4] + 
	11*c[2]*(21 + 90*c[3] + 10*c[4])) + 4*(-75075 + 18018*c[1]*c[1]*c[1] + 3718
	*c[2]*c[2]*c[2] + 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*
	(-363 + 154*c[3] - 118*c[4]) + 429*c[1]*c[1]*(63 + 54*c[2] - 18*c[3]-34*c[4])
	 + 42042*c[3]*c[4] - 12954*c[3]*c[3]*c[4] + 21333*c[4]*c[4] - 30822*c[3]*c[4]
	*c[4] - 8530*c[4]*c[4]*c[4] + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]
	*(21 + 38*c[3] + 10*c[4]) + c[4]*(-429 + 547*c[4]) + c[3]*(-99 + 602*c[4])) + 
	6*c[2]*(-91*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + c[4]*(715 + 987*c[4])))) + 
	2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 35*c[1] - 14*c[2] - 21*c[3]
	 - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]
	*c[3] + 33*c[1]*(35 + 14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(-21 + 45*c[3] 
	- 26*c[4]) - 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(-75075 + 9009*
	c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] 
	- 39*c[2]*c[2]*(-495 + 77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(42 + 24*c[2] - 
	9*c[3] - 28*c[4]) + 21021*c[3]*c[4] - 11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4] 
	- 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*
	c[3]*c[3] + 11*c[2]*(21 + 64*c[3] + 10*c[4]) + 12*c[3]*(-33 + 70*c[4]) + c[4]
	*(-429 + 547*c[4])) - 3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + 
	2*c[4]*(1859 + 1739*c[4]))))))/8648640.;
		bf_mom[3] = -((1 + b)*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-175 + 
	140*c[1] + 98*c[2] - 84*c[3] + 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 2607*c[2]
	*c[2] + 1386*c[3] + 2783*c[3]*c[3] - 66*c[1]*(35 + 14*c[2] + 33*c[3]-26*c[4])
	 - 561*c[4] - 2156*c[3]*c[4] + 2831*c[4]*c[4] - 11*c[2]*(147 + 180*c[3] + 
	166*c[4])) - 4*(375375 + 36036*c[1]*c[1]*c[1] + 32890*c[2]*c[2]*c[2] - 108537
	*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(189 + 90*c[2] + 36*c[3] - 
	134*c[4]) + 84084*c[3]*c[4] + 54798*c[3]*c[3]*c[4] - 110409*c[4]*c[4] - 61644
	*c[3]*c[4]*c[4] + 11258*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(2607 + 308*c[3] + 
	2270*c[4]) + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + c[3]*(1089 - 1918*c[4])
	 - 22*c[2]*(-21 + 77*c[3] - 10*c[4]) + 2*c[4]*(-429 + 547*c[4])) + 6*c[2]*
	(12701*c[3]*c[3] + 234*c[3]*(55 + 14*c[4]) + c[4]*(11869 + 11421*c[4]))) + 
	2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]* (35 + 35*c[1] - 14*c[2] - 21*
	c[3] - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*
	c[3]*c[3] + 33*c[1]*(35 + 14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(-21+45*c[3]
	 - 26*c[4]) - 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(-75075 + 9009*
	c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] 
	- 39*c[2]*c[2]*(-495 + 77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(42 + 24*c[2] - 
	9*c[3] - 28*c[4]) + 21021*c[3]*c[4] - 11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4] 
	- 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*
	c[3]*c[3] + 11*c[2]*(21 + 64*c[3] + 10*c[4]) + 12*c[3]*(-33 + 70*c[4]) + c[4]*
	(-429 + 547*c[4])) - 3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + 2*c[4]
	*(1859 + 1739*c[4])))) - b*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(175 + 70*
	c[1] - 154*c[2] - 42*c[3] - 2*c[4]) + 156*c[0]* (1617*c[1]*c[1] + 2871*c[2]*
	c[2] - 693*c[3] + 2849*c[3]*c[3] + 11*c[2]*(-231 + 90*c[3] - 238*c[4]) + 33*
	c[1]*(35 + 14*c[2] - 78*c[3] - 26*c[4]) - 33*c[4] + 1078*c[3]*c[4] + 2863*
	c[4]*c[4]) + 4*(-375375 + 18018*c[1]*c[1]*c[1] - 50050*c[2]*c[2]*c[2] + 
	111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-2871 + 154*c[3] - 
	2830*c[4]) + 429*c[1]*c[1]*(147 + 30*c[2] - 18*c[3] - 122*c[4]) + 42042*c[3]*
	c[4] - 51474*c[3]*c[3]*c[4] + 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] - 794*
	c[4]*c[4]*c[4] + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 142*
	c[3] + 10*c[4]) + 3*c[3]*(-429 + 518*c[4]) + c[4]*(-429 + 547*c[4])) - 6*c[2]*
	(16783*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + c[4]*
	(17017 + 16873*c[4]))))))/8648640.;
		bf_mom[4] = ((1 + b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-70 + 
	35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]*
	c[2] + 693*c[3] + 2266*c[3]*c[3] + 22*c[2]*(84 + 45*c[3] - 88*c[4]) + 33*c[1]*
	(-35 + 14*c[2] - 60*c[3] - 26*c[4]) + 132*c[4] + 1078*c[3]*c[4] + 2284*c[4]*
	c[4]) - 4*(30030 + 9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*c[2] - 44187*c[3]*
	c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(1122 + 77*c[3] - 1076*c[4]) + 429*
	c[1]*c[1]*(-63 + 18*c[2] - 9*c[3] - 50*c[4]) - 21021*c[3]*c[4] - 20922*c[3]*
	c[3]*c[4] - 44538*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 1364*c[4]*c[4]*c[4] - 
	3*c[2]*(12610*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 88*c[4]*(-143+141*c[4]))
	 + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(-21 + 116*c[3] + 10*c[4])
	 + c[4]*(429 + 547*c[4]) + 2*c[3]*(495 + 658*c[4]))) + b*b*(15015*c[0]*c[0]*
	c[0] + 1287*c[0]*c[0]*(-35 + 35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]
	* (924*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(-35 + 
	14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(21 + 45*c[3] - 26*c[4]) + 330*c[4] +
	1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(15015 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*
	c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(495+77*c[3]
	 - 398*c[4]) + 429*c[1]*c[1]*(-42 + 24*c[2] - 9*c[3] - 28*c[4]) - 21021*c[3]*
	c[4] - 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*
	c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(-21 + 
	64*c[3] + 10*c[4]) + 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) - 3*c[2]*
	(4264*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*(-1859 + 1739*c[4])))) + 
	b*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 42*c[2] - 6*c[4]) - 156*c[0]*
	(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 583*c[3]*c[3] + c[2]*(693 - 
	682*c[4]) - 99*c[4] + 579*c[4]*c[4]) + 4*(-15015 + 13442*c[2]*c[2]*c[2] + 
	22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] + 22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4] 
	+ 429*c[1]*c[1]*(21 + 6*c[2] + 22*c[4]) - 117*c[2]*c[2]*(-209 + 226*c[4]) - 
	78*c[1]*c[3]*(297 + 286*c[2] + 238*c[4]) + 6*c[2]*(4173*c[3]*c[3] + c[4]*
	(-4433 + 4465*c[4]))))))/ 2162160.;
		bf_mom[5] = ((1 + b)*(-15015*c[0]*c[0]*c[0] + 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) - 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-30030 + 6721*c[2]*c[2]*
	c[2] - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*c[1]*c[3] + 4815*
	c[3]*c[3])*c[4] - 967*c[4]*c[4]*c[4] + 3*c[2]*(-9009 + 429*c[1]*c[1] - 3718*
	c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])) + b*(15015*c[0]*c[0]*c[0] + 
	2574*c[0]*c[0]*(35*c[1] + 7*c[2] - 21*c[3] - 13*c[4]) + 156*c[0]*(-1155 + 
	693*c[1]*c[1] + 363*c[2]*c[2] + 517*c[3]*c[3] + 66*c[1]*(7*c[2] - 3*c[3] - 
	13*c[4]) + 1078*c[3]*c[4] + 547*c[4]*c[4] + 110*c[2]*(9*c[3] + c[4])) + 8*
	(30030 + 9009*c[1]*c[1]*c[1] + 1859*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*
	c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 59*c[4]) + 429*c[1]*c[1]*(27*c[2] - 9*c[3]
	 - 17*c[4]) + 16731*c[4] - 6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 4265*
	c[4]*c[4]*c[4] - 21*c[2]*(429 + 13*c[3]*c[3] - 234*c[3]*c[4] - 141*c[4]*c[4])
	 + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 602*c[3]*c[4] + 547*c[4]
	*c[4] + 22*c[2]*(19*c[3] + 5*c[4])))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287*
	c[0]*c[0]*(35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*c[0]*(-1155 + 462*
	c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*c[3] - 26*c[4]) + 
	33*c[1]*(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 563*c[4]*c[4]) + 4*
	(60060 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*
	c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(24*c[2] - 
	9*c[3] - 28*c[4]) + 12870*c[4] - 11292*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] 
	- 3298*c[4]*c[4]*c[4] - 6*c[2]*(-3003 + 2132*c[3]*c[3] - 819*c[3]*c[4] + 
	1739*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 840*c[3]*
	c[4] + 547*c[4]*c[4] + 22*c[2]*(32*c[3] + 5*c[4]))))))/4324320.;
		bf_mom[6] = ((1 + b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(70 + 
	35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 78*c[0]*(1386*c[1]*c[1] + 2244*c[2]*
	c[2] - 693*c[3] + 2266*c[3]*c[3] + 22*c[2]*(-84 + 45*c[3] - 88*c[4]) + 33*c[1]
	*(35 + 14*c[2] - 60*c[3] - 26*c[4]) - 132*c[4] + 1078*c[3]*c[4] + 2284*c[4]*
	c[4]) - 4*(-150150 + 9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*c[2] + 44187*c[3]*
	c[3] - 6903*c[3]*c[3]*c[3] + 429*c[1]*c[1]*(63 + 18*c[2] - 9*c[3] - 50*c[4]) 
	+ 21021*c[3]*c[4] - 20922*c[3]*c[3]*c[4] + 44538*c[4]*c[4] - 15411*c[3]*c[4]*
	c[4] - 1364*c[4]*c[4]*c[4] - 3*c[2]*(12610*c[3]*c[3] - 117*c[3]*(55 + 14*c[4])
	 + 88*c[4]*(143 + 141*c[4])) - 39*c[2]*c[2]*(77*c[3] - 2*(561 + 538*c[4])) + 
	39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 116*c[3] + 10*c[4]) + 
	c[4]*(-429 + 547*c[4]) + 2*c[3]*(-495 + 658*c[4]))) + b*b*(15015*c[0]*c[0]*
	c[0] + 1287*c[0]*c[0]*(35 + 35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]*
	(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(35 + 
	14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(-21 + 45*c[3] - 26*c[4]) - 330*c[4] + 
	1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(-75075 + 9009*c[1]*c[1]*c[1] - 4862*c[2]
	*c[2]*c[2] + 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-495 + 
	77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(42 + 24*c[2] - 9*c[3] - 28*c[4]) + 21021*
	c[3]*c[4] - 11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 
	3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 
	64*c[3] + 10*c[4]) + 12*c[3]*(-33 + 70*c[4]) + c[4]*(-429 + 547*c[4])) - 
	3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + 2*c[4]* (1859 + 1739*c[4])
	))) - b*(15015*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-35 + 42*c[2] - 6*c[4]) + 156*
	c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 583*c[3]*c[3] + 99*c[4] 
	+ 579*c[4]*c[4] - 11*c[2]*(63 + 62*c[4])) - 4*(75075 + 13442*c[2]*c[2]*c[2] - 
	22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] - 22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4]
	 + 429*c[1]*c[1]*(-21 + 6*c[2] + 22*c[4]) - 117*c[2]*c[2]*(209 + 226*c[4]) 
	- 78*c[1]*c[3]*(-297 + 286*c[2] + 238*c[4]) + 6*c[2]*(4173*c[3]*c[3] + c[4]*
	(4433 + 4465*c[4]))))))/ 2162160.;
		bf_mom[7] = ((1 + b)*(75075*c[0]*c[0]*c[0] - 2574*c[0]*c[0]*(70*c[1] 
	+ 49*c[2] - 42*c[3] + 17*c[4]) + 156*c[0]*(-5775 + 2079*c[1]*c[1] + 2607*c[2]
	*c[2] + 2783*c[3]*c[3] - 66*c[1]*(14*c[2] + 33*c[3] - 26*c[4]) - 2156*c[3]*
	c[4] + 2831*c[4]*c[4] - 22*c[2]*(90*c[3] + 83*c[4])) - 8*(-150150 + 18018*c[1]
	*c[1]*c[1] + 16445*c[2]*c[2]*c[2] + 54054*c[3] - 13806*c[3]*c[3]*c[3] - 429*
	c[1]*c[1]*(45*c[2] + 18*c[3] - 67*c[4]) - 21879*c[4] + 27399*c[3]*c[3]*c[4] 
	- 30822*c[3]*c[4]*c[4] + 5629*c[4]*c[4]*c[4] - 39*c[2]*c[2]*(154*c[3] + 1135*
	c[4]) + c[2]*(-63063 + 38103*c[3]*c[3] + 9828*c[3]*c[4] + 34263*c[4]*c[4]) +
	 78*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 959*c[3]*c[4] + 547*c[4]*
	c[4] + c[2]*(-847*c[3] + 110*c[4]))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]
	*c[0]*(35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*c[0]*(-1155+462*c[1]*c[1]
	 + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*c[3] - 26*c[4]) + 33*c[1]*
	(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 563*c[4]*c[4]) + 4*(60060 + 
	9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] 
	- 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(24*c[2] - 9*c[3]-28*c[4])
	 + 12870*c[4] - 11292*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*
	c[4] - 6*c[2]*(-3003 + 2132*c[3]*c[3] - 819*c[3]*c[4] + 1739*c[4]*c[4]) + 39*
	c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 840*c[3]*c[4] + 547*c[4]*c[4] + 
	22*c[2]*(32*c[3] + 5*c[4])))) - b*(75075*c[0]*c[0]*c[0] + 2574*c[0]*c[0]*
	(35*c[1] - 77*c[2] - 21*c[3] - c[4]) + 156*c[0]*(1617*c[1]*c[1] + 2871*c[2]*
	c[2] + 22*c[2]*(45*c[3] - 119*c[4]) + 7*(-825 + 407*c[3]*c[3] + 154*c[3]*c[4]
	 + 409*c[4]*c[4]) + 66*c[1]*(7*c[2] - 13*(3*c[3] + c[4]))) + 8*(150150 + 9009
	*c[1]*c[1]*c[1] - 25025*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 
	39*c[2]*c[2]*(77*c[3] - 1415*c[4]) + 429*c[1]*c[1]*(15*c[2] - 9*c[3] - 
	61*c[4]) + 1287*c[4] - 25737*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 397*c[4]
	*c[4]*c[4] + c[2]*(99099 - 50349*c[3]*c[3] + 4914*c[3]*c[4] - 50619*c[4]*c[4])
	 + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 1554*c[3]*c[4] + 547*c[4]
	*c[4] + 22*c[2]*(71*c[3] + 5*c[4]))))))/4324320.;
		bf_mom[8] = -((1+b)*(1+b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 156*c[0]*(-2310 + 693*c[1]*c[1] + 
	1122*c[2]*c[2] + 1133*c[3]*c[3] + 11*c[2]*(45*c[3] - 88*c[4]) + 33*c[1]*
	(7*c[2] - 30*c[3] - 13*c[4]) + 539*c[3]*c[4] + 1142*c[4]*c[4]) - 4*(120120 + 
	9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] 
	- 39*c[2]*c[2]*(77*c[3] - 1076*c[4]) + 429*c[1]*c[1]*(18*c[2] - 9*c[3] - 
	50*c[4]) + 5148*c[4] - 20922*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 1364*c[4]
	*c[4]*c[4] - 6*c[2]*(-12012 + 6305*c[3]*c[3] - 819*c[3]*c[4] + 6204*c[4]*c[4])
	 + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 1316*c[3]*c[4] + 547*c[4]
	*c[4] + 22*c[2]*(58*c[3] + 5*c[4]))) + b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*
	c[0]*(35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*c[0]*(-1155 + 462*c[1]*c[1]
	 + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*c[3] - 26*c[4]) + 33*c[1]*
	(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 563*c[4]*c[4]) + 4*(60060 + 
	9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] 
	- 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(24*c[2] - 9*c[3] - 
	28*c[4]) + 12870*c[4] - 11292*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4]-3298*c[4]
	*c[4]*c[4] - 6*c[2]*(-3003 + 2132*c[3]*c[3] - 819*c[3]*c[4] + 1739*c[4]*c[4])
	 + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 840*c[3]*c[4] + 547*c[4]
	*c[4] + 22*c[2]*(32*c[3] + 5*c[4]))))))/1081080.;
		break;
	case 6:
		bf_mom[0] = ((1 + b)*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(175 + 
	140*c[1] + 98*c[2] - 84*c[3] + 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 
	2607*c[2]*c[2] - 1386*c[3] + 2783*c[3]*c[3] - 66*c[1]*(-35 + 14*c[2] + 
	33*c[3] - 26*c[4]) + 561*c[4] - 2156*c[3]*c[4] + 2831*c[4]*c[4] - 11*c[2]*
	(-147 + 180*c[3] + 166*c[4])) - 4*(-375375 + 36036*c[1]*c[1]*c[1] + 
	32890*c[2]*c[2]*c[2] + 108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 429*c[1]
	*c[1]*(-189 + 90*c[2] + 36*c[3] - 134*c[4]) - 84084*c[3]*c[4] + 54798*c[3]
	*c[3]*c[4] + 110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] + 11258*c[4]*c[4]*c[4]
	 - 39*c[2]*c[2]*(-2607 + 308*c[3] + 2270*c[4]) + 78*c[1]*(726*c[2]*c[2] + 
	1034*c[3]*c[3] - 22*c[2]*(21 + 77*c[3] - 10*c[4]) + 2*c[4]*(429 + 547*c[4])
	 - c[3]*(1089 + 1918*c[4])) + 6*c[2]*(12701*c[3]*c[3] + 234*c[3]*(-55 + 14
	*c[4]) + c[4]*(-11869 + 11421*c[4]))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287
	*c[0]*c[0]*(-35 + 35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]* (924
	*c[1]*c[1] + 990*c[2]*c[2] + 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(-35 + 14
	*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(21 + 45*c[3] - 26*c[4]) + 330*c[4] + 
	1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*(75075 + 9009*c[1]*c[1]*c[1] - 4862
	*c[2]*c[2]*c[2] - 21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*
	(495 + 77*c[3] - 398*c[4]) + 429*c[1]*c[1]*(-42 + 24*c[2] - 9*c[3] - 28*c[4])
	 - 21021*c[3]*c[4] - 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]
	*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 
	11*c[2]*(-21 + 64*c[3] + 10*c[4]) + 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 
	547*c[4])) - 3*c[2]*(4264*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*
	(-1859 + 1739*c[4])))) + b*(-75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-175 + 
	70*c[1] - 154*c[2] - 42*c[3] - 2*c[4]) - 156*c[0]* (1617*c[1]*c[1] + 2871
	*c[2]*c[2] + 693*c[3] + 2849*c[3]*c[3] + 11*c[2]*(231 + 90*c[3] - 238*c[4]) 
	+ 33*c[1]*(-35 + 14*c[2] - 78*c[3] - 26*c[4]) + 33*c[4] + 1078*c[3]*c[4] + 
	2863*c[4]*c[4]) - 4*(375375 + 18018*c[1]*c[1]*c[1] - 50050*c[2]*c[2]*c[2] -
	 111111*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(2871 + 154*c[3] - 
	2830*c[4]) + 429*c[1]*c[1]*(-147 + 30*c[2] - 18*c[3] - 122*c[4]) - 42042
	*c[3]*c[4] - 51474*c[3]*c[3]*c[4] - 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] 
	- 794*c[4]*c[4]*c[4] + 78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*
	(-21 + 142*c[3] + 10*c[4]) + 3*c[3]*(429 + 518*c[4]) + c[4]*(429 + 547*c[4]))
	 - 6*c[2]*(16783*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + c[4]*(-17017 + 
	16873*c[4]))))))/8648640.;
		bf_mom[1] = ((1 + b)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] 
	+ 583*c[3]*c[3] + c[2]*(693 - 682*c[4]) - 99*c[4] + 579*c[4]*c[4]) + 4* 
	(-75075 + 13442*c[2]*c[2]*c[2] + 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] + 
	22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(21 + 6*c[2] + 22*c[4])
	 - 117*c[2]*c[2]*(-209 + 226*c[4]) - 78*c[1]*c[3]*(297 + 286*c[2] + 238*c[4])
	 + 6*c[2]*(4173*c[3]*c[3] + c[4]*(-4433 + 4465*c[4]))) + b*(15015*c[0]*c[0]*
	c[0] + 1287*c[0]*c[0]*(-35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4]) + 156*
	c[0]* (693*c[1]*c[1] + 363*c[2]*c[2] + 693*c[3] + 517*c[3]*c[3] + 33*c[1]*
	(-35 + 14*c[2] - 6*c[3] - 26*c[4]) + 429*c[4] + 1078*c[3]*c[4] + 547*c[4]*
	c[4] + 11*c[2]*(-21 + 90*c[3] + 10*c[4])) + 4*(75075 + 18018*c[1]*c[1]*c[1] 
	+ 3718*c[2]*c[2]*c[2] - 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]
	*(363 + 154*c[3] - 118*c[4]) + 429*c[1]*c[1]*(-63 + 54*c[2] - 18*c[3] - 
	34*c[4]) - 42042*c[3]*c[4] - 12954*c[3]*c[3]*c[4] - 21333*c[4]*c[4] - 
	30822*c[3]*c[4]*c[4] - 8530*c[4]*c[4]*c[4] - 6*c[2]*(91*c[3]*c[3] + (715 - 
	987*c[4])*c[4] - 117*c[3]*(-55 + 14*c[4])) + 78*c[1]*(363*c[2]*c[2] + 517*
	c[3]*c[3] + 11*c[2]*(-21 + 38*c[3] + 10*c[4]) + c[4]*(429 + 547*c[4]) + 
	c[3]*(99 + 602*c[4])))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*
	(-35 + 35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 
	990*c[2]*c[2] + 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(-35 + 14*c[2] - 24*c[3]
	 - 26*c[4]) + 22*c[2]*(21 + 45*c[3] - 26*c[4]) + 330*c[4] + 1078*c[3]*c[4] +
	 1126*c[4]*c[4]) + 4*(75075 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] - 
	21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(495 + 77*c[3] - 398*
	c[4]) + 429*c[1]*c[1]*(-42 + 24*c[2] - 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] -
	 11292*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*
	c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(-21 + 64*c[3] 
	+ 10*c[4]) + 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) - 3*c[2]*
	(4264*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*
	(-1859 + 1739*c[4]))))))/8648640.;
		bf_mom[2] = ((1 + b)*(-15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 
	42*c[2] - 6*c[4]) - 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] 
	+ 583*c[3]*c[3] + 99*c[4] + 579*c[4]*c[4] - 11*c[2]*(63 + 62*c[4])) + 4*
	(15015 + 13442*c[2]*c[2]*c[2] - 22737*c[3]*c[3] + 9630*c[3]*c[3]*c[4] - 
	22581*c[4]*c[4] - 1934*c[4]*c[4]*c[4] + 429*c[1]*c[1]*(-21 + 6*c[2] + 
	22*c[4]) - 117*c[2]*c[2]*(209 + 226*c[4]) - 78*c[1]*c[3]*(-297 + 286*c[2] + 
	238*c[4]) + 6*c[2]*(4173*c[3]*c[3] + c[4]*(4433 + 4465*c[4]))) + b*(15015*
	c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4])
	 + 156*c[0]*(693*c[1]*c[1] + 363*c[2]*c[2] - 693*c[3] + 517*c[3]*c[3] + 33*
	c[1]*(35 + 14*c[2] - 6*c[3] - 26*c[4]) - 429*c[4] + 1078*c[3]*c[4] + 547*
	c[4]*c[4] + 11*c[2]*(21 + 90*c[3] + 10*c[4])) + 4*(-15015 + 18018*c[1]*c[1]
	*c[1] + 3718*c[2]*c[2]*c[2] + 20163*c[3]*c[3] - 13806*c[3]*c[3]*c[3] - 39*
	c[2]*c[2]*(-363 + 154*c[3] - 118*c[4]) + 429*c[1]*c[1]*(63 + 54*c[2] - 
	18*c[3] - 34*c[4]) + 42042*c[3]*c[4] - 12954*c[3]*c[3]*c[4] + 21333*c[4]*c[4]
	 - 30822*c[3]*c[4]*c[4] - 8530*c[4]*c[4]*c[4] + 78*c[1]*(363*c[2]*c[2] + 
	517*c[3]*c[3] + 11*c[2]*(21 + 38*c[3] + 10*c[4]) + c[4]*(-429 + 547*c[4]) + 
	c[3]*(-99 + 602*c[4])) + 6*c[2]*(-91*c[3]*c[3] + 117*c[3]*(55 + 14*c[4]) + 
	c[4]*(715 + 987*c[4])))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*
	(35 + 35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 
	990*c[2]*c[2] - 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(35 + 14*c[2] - 24*c[3] 
	- 26*c[4]) + 22*c[2]*(-21 + 45*c[3] - 26*c[4]) - 330*c[4] + 1078*c[3]*c[4] 
	+ 1126*c[4]*c[4]) + 4*(-15015 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 
	21450*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-495 + 77*c[3] - 398*
	c[4]) + 429*c[1]*c[1]*(42 + 24*c[2] - 9*c[3] - 28*c[4]) + 21021*c[3]*c[4] - 
	11292*c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]
	*c[4]*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 64*c[3] 
	+ 10*c[4]) + 12*c[3]*(-33 + 70*c[4]) + c[4]*(-429 + 547*c[4])) - 3*c[2]*
	(4264*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + 2*c[4]*
	(1859 + 1739*c[4]))))))/8648640.;
		bf_mom[3] = ((1 + b)*(75075*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*(-175 + 
	140*c[1] + 98*c[2] - 84*c[3] + 34*c[4]) + 156*c[0]*(2079*c[1]*c[1] + 
	2607*c[2]*c[2] + 1386*c[3] + 2783*c[3]*c[3] - 66*c[1]*(35 + 14*c[2] + 
	33*c[3] - 26*c[4]) - 561*c[4] - 2156*c[3]*c[4] + 2831*c[4]*c[4] - 11*c[2]*
	(147 + 180*c[3] + 166*c[4])) - 4*(75075 + 36036*c[1]*c[1]*c[1] + 32890*c[2]
	*c[2]*c[2] - 108537*c[3]*c[3] - 27612*c[3]*c[3]*c[3] - 429*c[1]*c[1]*(189 + 
	90*c[2] + 36*c[3] - 134*c[4]) + 84084*c[3]*c[4] + 54798*c[3]*c[3]*c[4] - 
	110409*c[4]*c[4] - 61644*c[3]*c[4]*c[4] + 11258*c[4]*c[4]*c[4] - 39*c[2]*c[2]
	*(2607 + 308*c[3] + 2270*c[4]) + 78*c[1]*(726*c[2]*c[2] + 1034*c[3]*c[3] + 
	c[3]*(1089 - 1918*c[4]) - 22*c[2]*(-21 + 77*c[3] - 10*c[4]) + 2*c[4]*(-429 
	+ 547*c[4])) + 6*c[2]*(12701*c[3]*c[3] + 234*c[3]*(55 + 14*c[4]) + c[4]*
	(11869 + 11421*c[4]))) + 2*b*b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 
	35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*
	c[2] - 693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(35 + 14*c[2] - 24*c[3] - 26*c[4])
	 + 22*c[2]*(-21 + 45*c[3] - 26*c[4]) - 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]
	*c[4]) + 4*(-15015 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 21450*c[3]*
	c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-495 + 77*c[3] - 398*c[4]) + 
	429*c[1]*c[1]*(42 + 24*c[2] - 9*c[3] - 28*c[4]) + 21021*c[3]*c[4] - 11292*
	c[3]*c[3]*c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*
	c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 64*c[3] + 
	10*c[4]) + 12*c[3]*(-33 + 70*c[4]) + c[4]*(-429 + 547*c[4])) - 3*c[2]*
	(4264*c[3]*c[3] - 117*c[3]*(55 + 14*c[4]) + 2*c[4]*(1859 + 1739*c[4])))) -
	 b*(75075*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(175 + 70*c[1] - 154*c[2] - 
	42*c[3] - 2*c[4]) + 156*c[0]*(1617*c[1]*c[1] + 2871*c[2]*c[2] - 693*c[3] + 
	2849*c[3]*c[3] + 11*c[2]*(-231 + 90*c[3] - 238*c[4]) + 33*c[1]*(35 + 14*c[2]
	 - 78*c[3] - 26*c[4]) - 33*c[4] + 1078*c[3]*c[4] + 2863*c[4]*c[4]) + 4*
	(-75075 + 18018*c[1]*c[1]*c[1] - 50050*c[2]*c[2]*c[2] + 111111*c[3]*c[3] - 
	13806*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-2871 + 154*c[3] - 2830*c[4]) + 429*
	c[1]*c[1]*(147 + 30*c[2] - 18*c[3] - 122*c[4]) + 42042*c[3]*c[4] - 51474*c[3]
	*c[3]*c[4] + 111657*c[4]*c[4] - 30822*c[3]*c[4]*c[4] - 794*c[4]*c[4]*c[4] + 
	78*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 142*c[3] + 10*c[4]) +
	 3*c[3]*(-429 + 518*c[4]) + c[4]*(-429 + 547*c[4])) - 6*c[2]*(16783*c[3]*c[3]
	 - 117*c[3]*(55 + 14*c[4]) + c[4]*(17017 + 16873*c[4]))))))/8648640.;
		bf_mom[4] = -((1+b)*(1+b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(-70 + 35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 78*c[0]*(1386*c[1]*c[1] + 
	2244*c[2]*c[2] + 693*c[3] + 2266*c[3]*c[3] + 22*c[2]*(84 + 45*c[3] - 88*c[4])
	 + 33*c[1]*(-35 + 14*c[2] - 60*c[3] - 26*c[4]) + 132*c[4] + 1078*c[3]*c[4] +
	 2284*c[4]*c[4]) - 4*(150150 + 9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*c[2] - 
	44187*c[3]*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(1122 + 77*c[3] - 
	1076*c[4]) + 429*c[1]*c[1]*(-63 + 18*c[2] - 9*c[3] - 50*c[4]) - 21021*c[3]
	*c[4] - 20922*c[3]*c[3]*c[4] - 44538*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 
	1364*c[4]*c[4]*c[4] - 3*c[2]*(12610*c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 
	88*c[4]*(-143 + 141*c[4])) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 
	11*c[2]*(-21 + 116*c[3] + 10*c[4]) + c[4]*(429 + 547*c[4]) + 2*c[3]*(495 + 
	658*c[4]))) + b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(-35 + 35*c[1] - 
	14*c[2] - 21*c[3] - 10*c[4]) + 78*c[0]* (924*c[1]*c[1] + 990*c[2]*c[2] + 
	693*c[3] + 1100*c[3]*c[3] + 33*c[1]*(-35 + 14*c[2] - 24*c[3] - 26*c[4]) + 
	22*c[2]*(21 + 45*c[3] - 26*c[4]) + 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]
	*c[4]) + 4*(75075 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] - 21450*c[3]
	*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(495 + 77*c[3] - 398*c[4]) + 
	429*c[1]*c[1]*(-42 + 24*c[2] - 9*c[3] - 28*c[4]) - 21021*c[3]*c[4] - 11292
	*c[3]*c[3]*c[4] - 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]
	*c[4] + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(-21 + 64*c[3] + 
	10*c[4]) + 12*c[3]*(33 + 70*c[4]) + c[4]*(429 + 547*c[4])) - 3*c[2]*(4264*
	c[3]*c[3] - 117*c[3]*(-55 + 14*c[4]) + 2*c[4]*
	(-1859 + 1739*c[4]))))))/2162160.;
		bf_mom[5] = -((1 + b)*(-15015*c[0]*c[0]*c[0] + 7722*c[0]*c[0]*
	(7*c[2] - c[4]) - 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]
	*c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(30030 + 6721*c[2]
	*c[2]*c[2] - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*c[1]*c[3] 
	+ 4815*c[3]*c[3])*c[4] - 967*c[4]*c[4]*c[4] + 3*c[2]*(-9009 + 429*c[1]*c[1] 
	- 3718*c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])) + b*(15015*c[0]*c[0]
	*c[0] + 2574*c[0]*c[0]*(35*c[1] + 7*c[2] - 21*c[3] - 13*c[4]) + 156*c[0]*
	(-1155 + 693*c[1]*c[1] + 363*c[2]*c[2] + 517*c[3]*c[3] + 66*c[1]*(7*c[2] - 
	3*c[3] - 13*c[4]) + 1078*c[3]*c[4] + 547*c[4]*c[4] + 110*c[2]*(9*c[3] + c[4])
	) + 8*(-30030 + 9009*c[1]*c[1]*c[1] + 1859*c[2]*c[2]*c[2] + 27027*c[3] - 
	6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 59*c[4]) + 429*c[1]*c[1]*(27*
	c[2] - 9*c[3] - 17*c[4]) + 16731*c[4] - 6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]
	*c[4] - 4265*c[4]*c[4]*c[4] - 21*c[2]*(429 + 13*c[3]*c[3] - 234*c[3]*c[4] - 
	141*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 602*c[3]*
	c[4] + 547*c[4]*c[4] + 22*c[2]*(19*c[3] + 5*c[4])))) + 2*b*b*(15015*c[0]*
	c[0]*c[0] + 1287*c[0]*c[0]*(35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*c[0]
	*(-1155 + 462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*c[3] 
	- 26*c[4]) + 33*c[1]*(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 563*c[4]
	*c[4]) + 4*(-60060 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 27027*c[3] 
	- 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 429*c[1]*c[1]*
	(24*c[2] - 9*c[3] - 28*c[4]) + 12870*c[4] - 11292*c[3]*c[3]*c[4] - 15411*
	c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] - 6*c[2]*(-3003 + 2132*c[3]*c[3] - 
	819*c[3]*c[4] + 1739*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]
	*c[3] + 840*c[3]*c[4] + 547*c[4]*c[4] + 22*c[2]*
	(32*c[3] + 5*c[4]))))))/4324320.;
		bf_mom[6] = -((1+b)*(1+b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(70 + 35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 78*c[0]*(1386*c[1]*c[1] + 
	2244*c[2]*c[2] - 693*c[3] + 2266*c[3]*c[3] + 22*c[2]*(-84 + 45*c[3] - 88*
	c[4]) + 33*c[1]*(35 + 14*c[2] - 60*c[3] - 26*c[4]) - 132*c[4] + 1078*c[3]*
	c[4] + 2284*c[4]*c[4]) - 4*(-30030 + 9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*
	c[2] + 44187*c[3]*c[3] - 6903*c[3]*c[3]*c[3] + 429*c[1]*c[1]*(63 + 18*c[2] 
	- 9*c[3] - 50*c[4]) + 21021*c[3]*c[4] - 20922*c[3]*c[3]*c[4] + 44538*c[4]*
	c[4] - 15411*c[3]*c[4]*c[4] - 1364*c[4]*c[4]*c[4] - 3*c[2]*(12610*c[3]*c[3] 
	- 117*c[3]*(55 + 14*c[4]) + 88*c[4]*(143 + 141*c[4])) - 39*c[2]*c[2]*(77*c[3]
	 - 2*(561 + 538*c[4])) + 39*c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*
	(21 + 116*c[3] + 10*c[4]) + c[4]*(-429 + 547*c[4]) + 2*c[3]*(-495 + 658*
	c[4]))) + b*(15015*c[0]*c[0]*c[0] + 1287*c[0]*c[0]*(35 + 35*c[1] - 14*c[2] 
	- 21*c[3] - 10*c[4]) + 78*c[0]*(924*c[1]*c[1] + 990*c[2]*c[2] - 693*c[3] + 
	1100*c[3]*c[3] + 33*c[1]*(35 + 14*c[2] - 24*c[3] - 26*c[4]) + 22*c[2]*(-21 
	+ 45*c[3] - 26*c[4]) - 330*c[4] + 1078*c[3]*c[4] + 1126*c[4]*c[4]) + 4*
	(-15015 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 21450*c[3]*c[3] - 
	6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(-495 + 77*c[3] - 398*c[4]) + 429*c[1]*
	c[1]*(42 + 24*c[2] - 9*c[3] - 28*c[4]) + 21021*c[3]*c[4] - 11292*c[3]*c[3]*
	c[4] + 21957*c[4]*c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] + 39*
	c[1]*(363*c[2]*c[2] + 517*c[3]*c[3] + 11*c[2]*(21 + 64*c[3] + 10*c[4]) + 
	12*c[3]*(-33 + 70*c[4]) + c[4]*(-429 + 547*c[4])) - 3*c[2]*(4264*c[3]*c[3] 
	- 117*c[3]*(55 + 14*c[4]) + 2*c[4]*(1859 + 1739*c[4]))))))/2162160.;
		bf_mom[7] = -((1 + b)*(75075*c[0]*c[0]*c[0] - 2574*c[0]*c[0]*
	(70*c[1] + 49*c[2] - 42*c[3] + 17*c[4]) + 156*c[0]*(-5775 + 2079*c[1]*c[1] 
	+ 2607*c[2]*c[2] + 2783*c[3]*c[3] - 66*c[1]*(14*c[2] + 33*c[3] - 26*c[4]) - 
	2156*c[3]*c[4] + 2831*c[4]*c[4] - 22*c[2]*(90*c[3] + 83*c[4])) - 8*(150150 
	+ 18018*c[1]*c[1]*c[1] + 16445*c[2]*c[2]*c[2] + 54054*c[3] - 13806*c[3]*c[3]
	*c[3] - 429*c[1]*c[1]*(45*c[2] + 18*c[3] - 67*c[4]) - 21879*c[4] + 27399*
	c[3]*c[3]*c[4] - 30822*c[3]*c[4]*c[4] + 5629*c[4]*c[4]*c[4] - 39*c[2]*c[2]*
	(154*c[3] + 1135*c[4]) + c[2]*(-63063 + 38103*c[3]*c[3] + 9828*c[3]*c[4] + 
	34263*c[4]*c[4]) + 78*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 959*c[3]
	*c[4] + 547*c[4]*c[4] + c[2]*(-847*c[3] + 110*c[4]))) + 2*b*b*(15015*c[0]
	*c[0]*c[0] + 1287*c[0]*c[0]*(35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*
	c[0]*(-1155 + 462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*
	c[3] - 26*c[4]) + 33*c[1]*(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 
	563*c[4]*c[4]) + 4*(-60060 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 
	27027*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 
	429*c[1]*c[1]*(24*c[2] - 9*c[3] - 28*c[4]) + 12870*c[4] - 11292*c[3]*c[3]*
	c[4] - 15411*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] - 6*c[2]*(-3003 + 2132
	*c[3]*c[3] - 819*c[3]*c[4] + 1739*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]
	*c[2] + 517*c[3]*c[3] + 840*c[3]*c[4] + 547*c[4]*c[4] + 22*c[2]*(32*c[3] + 
	5*c[4])))) - b*(75075*c[0]*c[0]*c[0] + 2574*c[0]*c[0]*(35*c[1] - 77*c[2] - 
	21*c[3] - c[4]) + 156*c[0]*(1617*c[1]*c[1] + 2871*c[2]*c[2] + 22*c[2]*
	(45*c[3] - 119*c[4]) + 7*(-825 + 407*c[3]*c[3] + 154*c[3]*c[4] + 409*c[4]*
	c[4]) + 66*c[1]*(7*c[2] - 13*(3*c[3] + c[4]))) + 8*(-150150 + 9009*c[1]*c[1]
	*c[1] - 25025*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]*c[3] - 39*c[2]
	*c[2]*(77*c[3] - 1415*c[4]) + 429*c[1]*c[1]*(15*c[2] - 9*c[3] - 61*c[4]) + 
	1287*c[4] - 25737*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 397*c[4]*c[4]*c[4]
	 + c[2]*(99099 - 50349*c[3]*c[3] + 4914*c[3]*c[4] - 50619*c[4]*c[4]) + 39*
	c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 1554*c[3]*c[4] + 547*c[4]*c[4]
	 + 22*c[2]*(71*c[3] + 5*c[4]))))))/4324320.;
		bf_mom[8] = ((1+b)*(1+b)*(-30030*c[0]*c[0]*c[0] - 1287*c[0]*c[0]*
	(35*c[1] - 56*c[2] - 21*c[3] - 4*c[4]) - 156*c[0]*(-2310 + 693*c[1]*c[1] + 
	1122*c[2]*c[2] + 1133*c[3]*c[3] + 11*c[2]*(45*c[3] - 88*c[4]) + 33*c[1]*
	(7*c[2] - 30*c[3] - 13*c[4]) + 539*c[3]*c[4] + 1142*c[4]*c[4]) - 4*(-120120 
	+ 9009*c[1]*c[1]*c[1] - 18304*c[2]*c[2]*c[2] + 27027*c[3] - 6903*c[3]*c[3]
	*c[3] - 39*c[2]*c[2]*(77*c[3] - 1076*c[4]) + 429*c[1]*c[1]*(18*c[2] - 9*c[3]
	 - 50*c[4]) + 5148*c[4] - 20922*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 
	1364*c[4]*c[4]*c[4] - 6*c[2]*(-12012 + 6305*c[3]*c[3] - 819*c[3]*c[4] + 
	6204*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 1316*c[3]
	*c[4] + 547*c[4]*c[4] + 22*c[2]*(58*c[3] + 5*c[4]))) + b*(15015*c[0]*c[0]
	*c[0] + 1287*c[0]*c[0]* (35*c[1] - 14*c[2] - 21*c[3] - 10*c[4]) + 156*c[0]*
	(-1155 + 462*c[1]*c[1] + 495*c[2]*c[2] + 550*c[3]*c[3] + 11*c[2]*(45*c[3] - 
	26*c[4]) + 33*c[1]*(7*c[2] - 12*c[3] - 13*c[4]) + 539*c[3]*c[4] + 563*c[4]
	*c[4]) + 4*(-60060 + 9009*c[1]*c[1]*c[1] - 4862*c[2]*c[2]*c[2] + 27027*c[3] 
	- 6903*c[3]*c[3]*c[3] - 39*c[2]*c[2]*(77*c[3] - 398*c[4]) + 429*c[1]*c[1]*
	(24*c[2] - 9*c[3] - 28*c[4]) + 12870*c[4] - 11292*c[3]*c[3]*c[4] - 15411
	*c[3]*c[4]*c[4] - 3298*c[4]*c[4]*c[4] - 6*c[2]*(-3003 + 2132*c[3]*c[3] - 
	819*c[3]*c[4] + 1739*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]
	*c[3] + 840*c[3]*c[4] + 547*c[4]*c[4] + 22*c[2]*
	(32*c[3] + 5*c[4]))))))/1081080.;
		break;
	case 22:
		bf_mom[0] = (13*(11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 
	30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) 
	- 18*(7*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
        396*(-1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 66*(63*c[0]*c[0]
	 - 18*c[0]*(7 + 2*c[1] + 10*c[2]) + 4*(9*c[1]*(1 + c[1]) + (45 + 
	38*c[1])*c[2] + 7*c[2]*c[2]))* c[3] - 12*(517*(-1 + c[0] - 2*c[1]) - 
	14*c[2])*c[3]*c[3] - 4248*pow(c[3],3)) + 6* (143*(39*(-2 + c[0])*c[0] 
	- 156*(-1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(-1 + c[0] - 2*c[1])*c[2] 
	- 3068*c[2]*c[2] + 364*(-77 + 77*c[0] - 86*c[1] + 18*c[2])* c[3] + 
	8636*c[3]*c[3])*c[4] - 12*(-7111 + 7111*c[0] - 14222*c[1] + 1974*c[2] 
	+ 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[1] = (13*(11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 66*(9*(7*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] - 2*c[1]) - 
	14*c[2])*c[3]*c[3] - 4248*pow(c[3],3)) + 6* (143*(39*c[0]*(2 + c[0]) 
	- 156*(1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(1 + c[0] - 2*c[1])*c[2] 
	- 3068*c[2]*c[2] + 364*(77 + 77*c[0] - 86*c[1] + 18*c[2])*c[3] + 
	8636*c[3]*c[3])* c[4] - 12*(7111 + 7111*c[0] - 14222*c[1] + 1974*c[2] 
	+ 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[2] = (13*(11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*c[0]*(2 + 
	c[0]) + 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 38*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] + 2*c[1]) - 14*c[2])
	*c[3]*c[3] + 4248*pow(c[3],3)) + 6* (143*(39*c[0]*(2 + c[0]) + 156*
	(1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(1 + c[0] + 2*c[1])*c[2] - 
	3068*c[2]*c[2] - 364*(77 + 77*c[0] + 86*c[1] + 18*c[2])*c[3] + 
	8636*c[3]*c[3])* c[4] - 12*(7111 + 7111*c[0] + 14222*c[1] + 1974*c[2] 
	- 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[3] = (13*(11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 
	30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) 
	- 18*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
 	396*(-1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*(-2 
	+ c[0])*c[0] + 4*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] 
	+ 38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(-1 + c[0] + 2*c[1]) 
	- 14*c[2])*c[3]*c[3] + 4248*pow(c[3],3)) + 6* (143*(39*(-2 + c[0])*c[0] 
	+ 156*(-1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(-1 + c[0] + 2*c[1])*c[2] 
	- 3068*c[2]*c[2] - 364*(-77 + 77*c[0] + 86*c[1] + 18*c[2])* c[3] + 
	8636*c[3]*c[3])*c[4] - 12*(-7111 + 7111*c[0] + 14222*c[1] + 1974*c[2] 
	- 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[4] = (15015*pow(c[0],3) - 2574*c[0]*c[0]* (35*c[1] - 
	7*c[2] - 21*c[3] + 13*c[4]) + 156*c[0]*(-1155 + 693*c[1]*c[1] + 
	363*c[2]*c[2] + 517*c[3]*c[3] - 66*c[1]*(7*c[2] + 3*c[3] - 13*c[4]) - 
       	110*c[2]*(9*c[3] - c[4]) - 1078*c[3]*c[4] + 547*c[4]*c[4]) - 8*(-30030 
	+ 9009*pow(c[1],3) - 1859*pow(c[2],3) + 27027*c[3] - 6903*pow(c[3],3) 
	- 429*c[1]*c[1]* (27*c[2] + 9*c[3] - 17*c[4]) - 16731*c[4] + 6477*c[3]
	*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 4265*pow(c[4],3) - 39*c[2]*c[2]*
	(77*c[3] + 59*c[4]) + 21*c[2]*(429 + 13*c[3]*c[3] + 234*c[3]*c[4] - 
         141*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] - 
         22*c[2]*(19*c[3] - 5*c[4]) - 602*c[3]*c[4] + 547*c[4]*c[4])))/
	  1081080.;
		bf_mom[5] = (143*(-105*(-20 + c[0]*c[0]*(3 + c[0])) - 252*
	(1 + c[0])*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 
       	684*(1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 3432*c[1]*(27 + 27*c[0] 
	- 26*c[2])*c[3] - 156*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3] - 6*(143*
	(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] + 17628*c[2]
	*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] - 12*(7527*(1 + c[0]) - 
	8930*c[2])*c[4]*c[4] - 7736*pow(c[4],3))/ 540540.;
		bf_mom[6] = (15015*pow(c[0],3) + 2574*c[0]*c[0]* (35*c[1] + 
	7*c[2] - 21*c[3] - 13*c[4]) + 156*c[0]*(-1155 + 693*c[1]*c[1] + 
	363*c[2]*c[2] + 517*c[3]*c[3] + 66*c[1]*(7*c[2] - 3*c[3] - 13*c[4]) + 
        1078*c[3]*c[4] + 547*c[4]*c[4] + 110*c[2]*(9*c[3] + c[4])) + 8*(30030 
	+ 9009*pow(c[1],3) + 1859*pow(c[2],3) + 27027*c[3] - 6903*pow(c[3],3) 
	- 39*c[2]*c[2]*(77*c[3] - 59*c[4]) + 429*c[1]*c[1]*(27*c[2] - 9*c[3] 
	- 17*c[4]) + 16731*c[4] - 6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 
      	 4265*pow(c[4],3) - 21*c[2]* (429 + 13*c[3]*c[3] - 234*c[3]*c[4] - 
	141*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 
        602*c[3]*c[4] + 547*c[4]*c[4] + 22*c[2]*(19*c[3] + 5*c[4]))))/
	  1081080.;
		bf_mom[7] = (143*(-105*pow(-2 + c[0],2)*(1 + c[0]) - 252*(-1 
	+ c[0])*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] - 684*
	(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 3432*c[1]*(-27 + 27*c[0] - 
	26*c[2])*c[3] - 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] - 6*(143*(9*
	(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*
	c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] - 12*(7527*(-1 + c[0]) - 
	8930*c[2])*c[4]*c[4] - 7736*pow(c[4],3))/ 540540.;
		bf_mom[8] = (15015*pow(c[0],3) - 7722*c[0]*c[0]*(7*c[2] - c[4]) 
	+ 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) - 8*(-30030 + 6721*
	pow(c[2],3) - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*
	c[1]*c[3] + 4815*c[3]*c[3])* c[4] - 967*pow(c[4],3) + 3*c[2]*(-9009 + 
	429*c[1]*c[1] - 3718*c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])))
	/270270.;
		break;
	case 32:
		bf_mom[0] = (13*(11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) + 
	30*(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 
        18*(7*(-2 + c[0])*c[0] - 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
        396*(-1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 66*
	(63*c[0]*c[0] - 18*c[0]*(7 + 2*c[1] + 10*c[2]) + 4*(9*c[1]*(1 + c[1]) 
	+ (45 + 38*c[1])*c[2] + 7*c[2]*c[2]))* c[3] - 12*(517*(-1 + c[0] - 
	2*c[1]) - 14*c[2])*c[3]*c[3] - 4248*pow(c[3],3)) + 6* (143*(39*
	(-2 + c[0])*c[0] - 156*(-1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*
	(-1 + c[0] - 2*c[1])*c[2] - 3068*c[2]*c[2] + 364*(-77 + 77*c[0] - 
	86*c[1] + 18*c[2])* c[3] + 8636*c[3]*c[3])*c[4] - 12*(-7111 + 7111*c[0] 
	- 14222*c[1] + 1974*c[2] + 10274*c[3])* c[4]*c[4] + 
	34120*pow(c[4],3))/2162160.;
		bf_mom[1] = (13*(11*(21*(-5*pow(-2 + c[0],2)*(1 + c[0]) - 30*
	(-2 + c[0])*c[0]*c[1] - 36*(-1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 
	18*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 36*c[1]*c[1])*c[2] - 
        396*(-1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*
	(-2 + c[0])*c[0] + 4*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] 
	+ 38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(-1 + c[0] + 2*c[1]) - 
	14*c[2])*c[3]*c[3] + 4248*pow(c[3],3)) + 6* (143*(39*(-2 + c[0])*c[0] + 
	156*(-1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(-1 + c[0] + 2*c[1])*c[2] - 
        3068*c[2]*c[2] - 364*(-77 + 77*c[0] + 86*c[1] + 18*c[2])* c[3] + 
	8636*c[3]*c[3])*c[4] - 12*(-7111 + 7111*c[0] + 14222*c[1] + 1974*c[2] 
	- 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[2] = (13*(11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) - 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] - 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] + 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) + 66*(9*(7*c[0]*
	(2 + c[0]) + 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] + 2*c[1]) - 
	14*c[2])*c[3]*c[3] + 4248*pow(c[3],3)) + 6* (143*(39*c[0]*(2 + c[0]) + 
	156*(1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(1 + c[0] + 2*c[1])*c[2] - 
	3068*c[2]*c[2] - 364*(77 + 77*c[0] + 86*c[1] + 18*c[2])*c[3] + 
	8636*c[3]*c[3])* c[4] - 12*(7111 + 7111*c[0] + 14222*c[1] + 1974*c[2] 
	- 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[3] = (13*(11*(21*(100 - 5*c[0]*c[0]*(3 + c[0]) + 30*c[0]*
	(2 + c[0])*c[1] - 36*(1 + c[0])*c[1]*c[1] + 24*pow(c[1],3)) - 18*
	(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 36*c[1]*c[1])* c[2] - 396*
	(1 + c[0] - 2*c[1])*c[2]*c[2] - 104*pow(c[2],3)) - 66*(9*(7*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 
	38*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 12*(517*(1 + c[0] - 2*c[1]) - 
	14*c[2])*c[3]*c[3] - 4248*pow(c[3],3)) + 6* (143*(39*c[0]*(2 + c[0]) 
	- 156*(1 + c[0])*c[1] + 68*c[1]*c[1]) - 2860*(1 + c[0] - 2*c[1])*c[2] 
	- 3068*c[2]*c[2] + 364*(77 + 77*c[0] - 86*c[1] + 18*c[2])*c[3] + 
	8636*c[3]*c[3])* c[4] - 12*(7111 + 7111*c[0] - 14222*c[1] + 1974*c[2] 
	+ 10274*c[3])* c[4]*c[4] + 34120*pow(c[4],3))/2162160.;
		bf_mom[4] = (143*(-105*pow(-2 + c[0],2)*(1 + c[0]) - 252*
	(-1 + c[0])*c[1]*c[1] + 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] - 
       	684*(-1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 3432*c[1]*(-27 + 27*c[0] 
	- 26*c[2])*c[3] - 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] - 6*(143*
	(9*(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*
	c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] - 12*(7527*
	(-1 + c[0]) - 8930*c[2])*c[4]*c[4] - 7736*pow(c[4],3))/ 540540.;
		bf_mom[5] = (15015*pow(c[0],3) + 2574*c[0]*c[0]* (35*c[1] + 
	7*c[2] - 21*c[3] - 13*c[4]) + 156*c[0]*(-1155 + 693*c[1]*c[1] + 
	363*c[2]*c[2] + 517*c[3]*c[3] + 66*c[1]*(7*c[2] - 3*c[3] - 13*c[4]) + 
	1078*c[3]*c[4] + 547*c[4]*c[4] + 110*c[2]*(9*c[3] + c[4])) + 8*(30030 + 
	9009*pow(c[1],3) + 1859*pow(c[2],3) + 27027*c[3] - 6903*pow(c[3],3) - 
	39*c[2]*c[2]*(77*c[3] - 59*c[4]) + 429*c[1]*c[1]*(27*c[2] - 9*c[3] - 
	17*c[4]) + 16731*c[4] - 6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] - 
       	4265*pow(c[4],3) - 21*c[2]* (429 + 13*c[3]*c[3] - 234*c[3]*c[4] - 
	141*c[4]*c[4]) + 39*c[1]*(-1155 + 363*c[2]*c[2] + 517*c[3]*c[3] + 
        602*c[3]*c[4] + 547*c[4]*c[4] + 22*c[2]*(19*c[3] + 5*c[4]))))
	/ 1081080.;
		bf_mom[6] = (143*(-105*(-20 + c[0]*c[0]*(3 + c[0])) - 252*
	(1 + c[0])*c[1]*c[1] + 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] - 
	684*(1 + c[0])*c[2]*c[2] + 376*pow(c[2],3)) + 3432*c[1]*(27 + 27*c[0] 
	- 26*c[2])*c[3] - 156*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3] - 6*(143*
	(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] + 17628*c[2]
	*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] - 12*(7527*(1 + c[0]) - 
	8930*c[2])*c[4]*c[4] - 7736*pow(c[4],3))/ 540540.;
		bf_mom[7] = (15015*pow(c[0],3) - 2574*c[0]*c[0]* (35*c[1] - 
	7*c[2] - 21*c[3] + 13*c[4]) + 156*c[0]*(-1155 + 693*c[1]*c[1] + 
	363*c[2]*c[2] + 517*c[3]*c[3] - 66*c[1]*(7*c[2] + 3*c[3] - 13*c[4]) - 
        110*c[2]*(9*c[3] - c[4]) - 1078*c[3]*c[4] + 547*c[4]*c[4]) - 
    	8*(-30030 + 9009*pow(c[1],3) - 1859*pow(c[2],3) + 27027*c[3] - 
       	6903*pow(c[3],3) - 429*c[1]*c[1]* (27*c[2] + 9*c[3] - 17*c[4]) - 
	16731*c[4] + 6477*c[3]*c[3]*c[4] - 15411*c[3]*c[4]*c[4] + 4265*
	pow(c[4],3) - 39*c[2]*c[2]*(77*c[3] + 59*c[4]) + 21*c[2]*(429 + 
	13*c[3]*c[3] + 234*c[3]*c[4] - 141*c[4]*c[4]) + 39*c[1]*(-1155 + 
	363*c[2]*c[2] + 517*c[3]*c[3] - 22*c[2]*(19*c[3] - 5*c[4]) - 602*c[3]*
	c[4] + 547*c[4]*c[4])))/ 1081080.;
		bf_mom[8] = (15015*pow(c[0],3) - 7722*c[0]*c[0]*(7*c[2] - c[4]) 
	+ 156*c[0]*(-1155 + 231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 
	583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) - 8*(-30030 + 6721*
	pow(c[2],3) - 13221*c[2]*c[2]*c[4] + (3861 + 4719*c[1]*c[1] - 9282*c[1]
	*c[3] + 4815*c[3]*c[3])* c[4] - 967*pow(c[4],3) + 3*c[2]*(-9009 + 
	429*c[1]*c[1] - 3718*c[1]*c[3] + 4173*c[3]*c[3] + 4465*c[4]*c[4])))
	/270270.;
		break;
	case 100:
		bf_mom[0] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*
	c[3]) + 4*(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*
	(-1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 
	796*c[2]*c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*
	c[3]*c[3])*c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(-5*(20 + (-3 + c[0])*
	c[0]*c[0]) - 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 + c[0])*c[1]*c[1] - 8*c[1]*
	c[1]*c[1]) + 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2]
	 - 4*(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*
	(-2 + c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 
	2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])
	*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + 
	c[0])*c[1] + 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*
	c[2] - 28*(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 
	12*(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] 
	+ 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 
	15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 
	+ 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*c[2]
	 + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])* c[3]*c[3]
	 - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] + 858*(-1 + c[0])*c[1]
	 + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]*c[2]) - 182*
	(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4])) + a*(-45045*c[0]*c[0]*c[0] + 9009*c[0]*c[0]*(15 + 10*c[1] 
	+ 10*c[2] - 6*c[3] + 2*c[4]) - 1716*c[0]*(21*(5*c[1]*(1 + c[1]) + 5*c[2] - 
	2*c[1]*c[2] + 7*c[2]*c[2]) - 9*(7 + 14*c[1] + 10*c[2])*c[3] + 153*c[3]*c[3] + 
	(21 + 78*c[1] - 114*c[2] - 98*c[3])*c[4] + 155*c[4]*c[4]) + 12*(13*(33*(-175
	 + 7*c[1]*c[1]*(5 + 2*c[1]) - 14*c[1]*(1 + c[1])*c[2] + (49 + 22*c[1])*c[2]*
	c[2] + 18*c[2]*c[2]*c[2]) - 22*(9*c[1]*(7 + c[1]) + 45*(1 + 2*c[1])*c[2] + 
	7*c[2]*c[2])*c[3] + 11*(153 + 94*c[1] + 118*c[2])*c[3]*c[3] - 354*c[3]*c[3]*
	c[3]) + 26*(11*(39*c[1]*(1 + c[1]) + (-57 + 10*c[1])*c[2] - 67*c[2]*c[2]) - 
	7*(77 + 154*c[1] - 18*c[2])*c[3] + 413*c[3]*c[3])*c[4] + (22165 + 14222*c[1] 
	+ 15886*c[2] - 10274*c[3])*c[4]*c[4] + 1554*c[4]*c[4]*c[4]) + 2*b*(143*(105*
	(20 + (-3 + c[0])*c[0]*c[0]) + 252*(-1 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])
	*c[0] + 4*c[1]*c[1])*c[2] + 684*(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 
	3432*c[1]*(-27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])
	*c[3]*c[3] + 6*(143*(9*(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])
	*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*
	(-1 + c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[1] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]
	*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*(-1 + 
	c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*c[3]) + 
	4*(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 9*(15 + 
	11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + (45 + 
	64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 + c[0])*c[1] + 
	616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(-77 
	+ 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(-7319 + 
	7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*
	c[4]) + b*(39*(33*(700 + 35*c[0]*c[0]*c[0] + 28*c[1]*c[1]*(-5 + 2*c[1]) + 
	35*c[0]*c[0]*(-3 + 2*c[1] - 2*c[2]) + 56*(-1 + c[1])*c[1]*c[2] + 4*(-49 + 
	22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2] + 28*c[0]*(5*(-1 + c[1])*c[1] + (5 + 
	2*c[1])*c[2] + 7*c[2]*c[2])) - 22*(9*(7*(-2 + c[0])*c[0] + 28*(-1 + c[0])
	*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 
	44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 
	78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(-57 + 
	57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(-1 + c[0] + 2*c[1]) + 
	18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(-22165 + 22165*c[0] + 14222*c[1] -
	 15886*c[2] - 10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*
	(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 15*(-2 + c[0])*c[0]*c[1] + 24*(-1+c[0])
	*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 
	14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]
	*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) 
	- 4*(-45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0]
	 + 517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + 
	c[0])*c[0] + 858*(-1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0]-5*c[1])
	*c[2] - 796*c[2]*c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*
	c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(700 + 35*c[0]*c[0]*c[0]
	 - 28*c[1]*c[1]*(5 + 2*c[1]) + 56*c[1]*(1 + c[1])*c[2] - 4*(49 + 22*c[1])*c[2]
	*c[2] - 72*c[2]*c[2]*c[2] - 35*c[0]*c[0]*(3 + 2*c[1] + 2*c[2]) + 28*c[0]*
	(5*c[1]*(1 + c[1]) + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2])) + 22*(63*c[0]*c[0] 
	- 18*c[0]*(7 + 14*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + c[1]) + 45*(1+2*c[1])*c[2] 
	+ 7*c[2]*c[2]))*c[3] + 44*(-153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]*c[3] + 
	1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*c[0] + 52*c[1]*(1 + c[1]) - 2*c[0]*(7 + 
	26*c[1])) + 44*(-57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] + 28*(77*(-1 
	+ c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(22165*c[0] - 
	13*(1705 + 1094*c[1] + 1222*c[2]) + 10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]
	*c[4] + 2*b*(143*(105*(20 + (-3 + c[0])*c[0]*c[0]) + 252*(-1 + c[0])*c[1]*c[1]
	 - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 684*(-1 + c[0])*c[2]*c[2] - 
	376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*
	(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*(-2 + c[0])*c[0] - 44*c[1]*c[1])
	 - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3]-6420*c[3]*c[3])
	*c[4] + 12*(7527*(-1 + c[0]) - 8930*c[2])*c[4]*c[4] + 
	7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[2] = -((a - b)*(2*a*a*(13*(11* (21*(5*(-1 + c[0])*(2+c[0])*
	(2+c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]
	*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*
	(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*
	c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*
	c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])*c[1]
	 + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 182*
	(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 + 
	7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*
	c[4]) + b*(39*(33*(7*(5*(-1 + c[0])*(2+c[0])*(2+c[0]) + 10*c[0]*(2 + c[0])*
	c[1] + 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) - 
	4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49 + 49*c[0] + 22*c[1])*c[2]*c[2]
	 - 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*(153 + 
	153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 78*(33*
	(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] - 
	10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 
	1652*c[3]*c[3])*c[4] + 12*(22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 
	10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*(-1+c[0])
	*(2+c[0])*(2+c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 
	12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*
	c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 
	33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0]
	 + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(550 + 550*c[0] + 517*c[1] - 
	328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) + 
	858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] - 5*c[1])*c[2] - 
	796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*
	c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(7*(5*(-1 + c[0])*(2+c[0])*(2+c[0])
	 - 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) - 
	14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49+49*c[0]
	 - 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) - 28*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] - 2*c[1])*c[2] + 28*c[2]*c[2])
	*c[3] + 44*(153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]
	*c[3]) - 78*(33*(7*c[0]*(2 + c[0]) - 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 
	44*(57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] + 28*(77*(1 +c[0]-2*c[1]) 
	+ 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(22165 + 22165*c[0] - 14222*c[1] 
	- 15886*c[2] + 10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*
	(-1 + c[0])*(2+c[0])*(2+c[0]) + 252*(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + 
	c[0]) + 4*c[1]*c[1])*c[2] + 684*(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 
	3432*c[1]*(27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])*
	c[3]*c[3] + 6*(143*(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] 
	+ 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(1 + 
	c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[3] = -((a - b)*(2*a*a*(13*(11* (21*(5*(-1 + c[0])*(2+c[0])*
	(2+c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]
	*c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*
	(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]
	*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])
	*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 
	182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4]) + b*(39*(33*(7*(20 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + 
	c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + 
	c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] + 22*c[1])
	*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])
	*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 
	44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 
	78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 
	57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 2*c[1]) + 
	18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] + 14222*c[1] - 
	15886*c[2] - 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*
	(5*(-1 + c[0])*(2+c[0])*(2+c[0]) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])
	*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] 
	- 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]
	*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*
	(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(550 + 550*c[0] + 
	517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*
	(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0]-5*c[1])
	*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528
	*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(7*(20 - 5*c[0]*c[0]*(3+c[0])
	 + 10*c[0]*(2 + c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) + 
	14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(49 + 
	49*c[0] - 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) 
	- 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] - 2*c[1])*c[2] + 28*c[2]
	*c[2])*c[3] - 44*(153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]
	*c[3]*c[3]) + 78*(33*(7*c[0]*(2 + c[0]) - 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) 
	+ 44*(57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] + 28*(77*(1 + c[0] - 
	2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] - 
	14222*c[1] - 15886*c[2] + 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 
	2*b*(143*(105*(-1 + c[0])*(2+c[0])*(2+c[0]) + 252*(1 + c[0])*c[1]*c[1] - 18*
	(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] + 684*(1 + c[0])*c[2]*c[2] - 376*c[2]
	*c[2]*c[2]) - 3432*c[1]*(27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(1 + c[0]) 
	- 642*c[2])*c[3]*c[3] + 6*(143*(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*
	(1 + c[0])*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 
	12*(7527*(1 + c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[4] = (429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*c[4] 
	- 84*c[0]*(5*c[2] + c[4])) - 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*c[0]*
	(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*c[3]) + 
	5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] + 39*c[1] 
	+ 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 5*c[0]*c[0] 
	- 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*
	(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 
	26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*c[0]*c[0] - 84*c[0]
	*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3])
	 + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4])) + 360360*(-a + b + ((a - b)*
	 (-15015*c[0]*c[0]*c[0] + 6006*c[0]*c[0]*(5*c[2] + c[4]) - 572*c[0]*(3*
	(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3]) - 114*c[2]*c[4] + 
	155*c[4]*c[4]) + 8*(143*c[2]*(-21*c[1]*c[1] + 27*c[2]*c[2] - 90*c[1]*c[3] + 
	59*c[3]*c[3]) + 13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]
	*c[3])*c[4] + 7943*c[2]*c[4]*c[4] + 777*c[4]*c[4]*c[4])))/120120.) + (a - b)*
	(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 
	132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*
	c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2]
	 + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])
	*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*
	(165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*
	(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] 
	+ 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 
	36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]*c[1]
	 - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + 
	c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]
	*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2]) 
	+ 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* 
	c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*
	c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3]
	 + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(15015 - 13*c[2]*
	(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*
	c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[5] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) - 18*
	(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(-63 + 
	(9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 12*c[3]) + 
	5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*(11*(-60 + 
	15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 5*c[1])*c[2] - 
	796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	 c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4]) + b*(39*(33*(35*c[0]*c[0]*c[0] + 70*c[0]*c[0]*(c[1] - 
	c[2]) + 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 7*c[2]*c[2]) + 8*(7*(-10 
	- 5*c[1] + c[1]*c[1]*c[1]) + 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2] - 
	9*c[2]*c[2]*c[2])) - 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] + 4*(-7 + c[1]*c[1])) 
	- 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*(153*c[0] + 94*c[1] - 
	118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 78*(231*c[0]*c[0] + 44*(-21 + 
	39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 44*c[0]*(39*c[1] + 57*c[2] - 
	49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(22165*c[0]
	 + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 
	2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*(-20 + 3*c[1]*(-5 + c[1]*c[1])) + 
	63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]
	*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 
	15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] 
	+ 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550*c[0] + 517*c[1] - 328*c[2])
	*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]
	*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*
	(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] 
	+ 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + 
	a*(3*(15015*c[0]*c[0]*c[0] - 6006*c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4])
	 + 572*c[0]*(3*(-105 + 35*c[1]*c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*
	c[3] - 14*c[1]*(c[2] + 3*c[3])) + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 
	155*c[4]*c[4]) + 8*(13*(-33*(7*(10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*
	c[1])*c[2] + 11*c[1]*c[2]*c[2] + 9*c[2]*c[2]*c[2]) + 11*(9*c[1]*c[1] + 
	90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))*c[3] - 11*(47*c[1] + 59*c[2])*c[3]*c[3] 
	+ 177*c[3]*c[3]*c[3]) - 13*(-231 + 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*
	(5*c[2] - 49*c[3]) + 126*c[2]*c[3] + 413*c[3]*c[3])*c[4] - (7111*c[1] + 
	7943*c[2] - 5137*c[3])*c[4]*c[4] - 777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*
	c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*
	c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*
	c[4]*c[4]) + 8*(-13*(2310 + c[2]* (-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 
	858*c[1]*c[3] + 963*c[3]*c[3])) - 3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 
	3094*c[1]*c[3] + 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*c[4] + 
	967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[6] = (-429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*c[4]
	 - 84*c[0]*(5*c[2] + c[4])) + 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*c[0]*
	(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*c[3]) + 
	5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] + 39*c[1] +
	 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 5*c[0]*c[0] 
	- 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*
	(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 
	26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*c[0]*c[0] - 84*c[0]
	*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3])
	 + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4])) + 360360*(-a + b + ((a - b)*
	 (-15015*c[0]*c[0]*c[0] + 6006*c[0]*c[0]*(5*c[2] + c[4]) - 572*c[0]*(3*(35*
	c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3]) - 114*c[2]*c[4] + 
	155*c[4]*c[4]) + 8*(143*c[2]*(-21*c[1]*c[1] + 27*c[2]*c[2] - 90*c[1]*c[3] + 
	59*c[3]*c[3]) + 13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*
	c[3]*c[3])*c[4] + 7943*c[2]*c[4]*c[4] + 777*c[4]*c[4]*c[4])))/120120.) + 
	(a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 
	3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3]
	 + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*c[1]
	*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*
	(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*
	c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]
	*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 
	45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 
	63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2])
	 + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3]
	 + 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*
	(3*c[1] - 2*c[2]) + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] 
	+ 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2]
	 + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] 
	- 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 
	594*c[1]*c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(15015 - 
	13*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*
	(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[7] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) - 
	18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*(11*
	(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 5*c[1])
	*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]
	*c[3])* c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4]
	 - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(-35*c[0]*c[0]*c[0] - 70*c[0]*c[0]*(c[1]
	 - c[2]) - 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 7*c[2]*c[2]) + 8*(70 + 
	35*c[1] - 7*c[1]*c[1]*c[1] - 7*(5 + c[1]*c[1])*c[2] - 11*c[1]*c[2]*c[2] + 
	9*c[2]*c[2]*c[2])) + 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] + 4*(-7 + c[1]*c[1])) 
	- 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*(153*c[0] + 94*c[1] - 
	118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*(231*c[0]*c[0] + 44*(-21 + 
	39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 44*c[0]*(39*c[1] + 57*c[2] - 
	49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165*c[0]
	 + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 
	2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*(-20 + 3*c[1]*(-5 + c[1]*c[1])) + 
	63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]
	*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 
	15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] 
	+ 64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550*c[0] + 517*c[1] - 328*c[2])
	*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]
	*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*
	(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] 
	+ 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + 
	a*(3*(-15015*c[0]*c[0]*c[0] + 6006*c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + 
	c[4]) - 572*c[0]*(3*(-105 + 35*c[1]*c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 
	51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3])) + 2*(39*c[1] - 57*c[2] - 49*c[3])
	*c[4] + 155*c[4]*c[4]) + 8*(13*(33*(7*(10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*
	(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2] + 9*c[2]*c[2]*c[2]) - 11*(9*c[1]
	*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))*c[3] + 11*(47*c[1] + 59*c[2])*
	c[3]*c[3] - 177*c[3]*c[3]*c[3]) + 13*(-231 + 429*c[1]*c[1] - 737*c[2]*c[2] + 
	22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3] + 413*c[3]*c[3])*c[4] + (7111*c[1]
	 + 7943*c[2] - 5137*c[3])*c[4]*c[4] + 777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]
	*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*(11*(-105 + 21*c[1]
	*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*
	c[4]*c[4]) + 8*(-13*(2310 + c[2]* (-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 
	858*c[1]*c[3] + 963*c[3]*c[3])) - 3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 
	3094*c[1]*c[3] + 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*c[4] + 
	967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + 
	((a - b)*c[4])/15. + ((a - b)*(b*b*(7*(10 + 5*c[0] + 5*c[1] - 2*c[2] - 3*c[3])
	 - 10*c[4]) + a*a* (7*(10 + 5*c[0] - 5*c[1] - 2*c[2] + 3*c[3]) - 10*c[4]) + 
	a*b*(70 + 35*c[0] - 42*c[2] + 6*c[4])))/210. + (a - b + ((a - b)*(15015*c[0]
	*c[0]*c[0] - 6006*c[0]*c[0]*(5*c[2] + c[4]) + 572*c[0]*(3*(35*c[1]*c[1] + 
	49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) 
	+ 8*(143*c[2]*(21*c[1]*c[1] - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) - 
	13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] - 
	7943*c[2]*c[4]*c[4] - 777*c[4]*c[4]*c[4])))/120120.)/3. - ((a - b)*(b*b*
	(13*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 
	132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]
	*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 
	99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])
	*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*
	(165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*
	(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] 
	+ 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 
	36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]
	*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 
	531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*
	(3*c[1] - 2*c[2]) + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] 
	+ 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2]
	 + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0]
	 - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 
	594*c[1]*c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(15015 - 
	13*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*
	(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/1081080.;
		break;
	case 101:
		bf_mom[0] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*
	c[3]) + 4*(11*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 +
	 c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*
	c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*
	c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(-5*pow(-2 + c[0],2)*(1 + c[0])
	 - 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) 
	+ 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49
	 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*(-2 + c[0])
	*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] 
	+ 28*c[2]*c[2])*c[3] - 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] 
	+ 1416*c[3]*c[3]*c[3]) + 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])*c[1] 
	+ 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*
	(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] + 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*
	(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*
	(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 
	+ c[0])*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 
	64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 
	328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] 
	+ 858*(-1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] 
	- 796*c[2]*c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*
	c[3]*c[3])*c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2]-5137*c[3])
	* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(-45045*c[0]*c[0]*c[0] + 9009*c[0]
	*c[0]*(15 + 10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]) - 1716*c[0]*(21*(5*c[1]*
	(1 + c[1]) + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2]) - 9*(7 + 14*c[1] + 10*c[2])
	*c[3] + 153*c[3]*c[3] + (21 + 78*c[1] - 114*c[2] - 98*c[3])*c[4] + 155*c[4]
	*c[4]) + 12*(13*(33*(-35 + 7*c[1]*c[1]*(5 + 2*c[1]) - 14*c[1]*(1 + c[1])
	*c[2] + (49 + 22*c[1])*c[2]*c[2] + 18*c[2]*c[2]*c[2]) - 22*(9*c[1]*
	(7 + c[1]) + 45*(1 + 2*c[1])*c[2] + 7*c[2]*c[2])*c[3] + 11*(153 + 94*c[1] 
	+ 118*c[2])*c[3]*c[3] - 354*c[3]*c[3]*c[3]) + 26*(11*(39*c[1]*(1 + c[1]) + 
	(-57 + 10*c[1])*c[2] - 67*c[2]*c[2]) - 7*(77 + 154*c[1] - 18*c[2])*c[3] + 
	413*c[3]*c[3])*c[4] + (22165 + 14222*c[1] + 15886*c[2] - 10274*c[3])*c[4]
	*c[4] + 1554*c[4]*c[4]*c[4]) + 2*b*(143*(105*pow(-2 + c[0],2)*(1 + c[0]) + 
	252*(-1 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 
	684*(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] 
	- 26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*
	(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2]
	 + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(-1 + c[0])-8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[1] = ((a - b)*(2*a*a*(13*(11* (21*(-100 + 5*c[0]*c[0]*(3 + 
	c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*
	c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]
	*(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0]-64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]
	*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])
	*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 
	182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4]) + b*(39*(33*(7*(100 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*
	(2 + c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] + 
	22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) + 28*
	(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2]+28*c[2]*c[2])
	*c[3] - 44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]
	*c[3]) + 78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 
	44*(57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 
	2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] 
	+ 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] 
	+ 2*b*(13*(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] 
	+ 24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 
	14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]
	*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1]
	 + 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 
	12*(550 + 550*c[0] + 517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3])
	 - 6*(13*(165*c[0]*(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*
	(13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] +120*c[1]
	 + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1]
	 - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*
	(7*(100 - 5*c[0]*c[0]*(3 + c[0]) + 10*c[0]*(2 + c[0])*c[1] - 20*(1 + c[0])
	*c[1]*c[1] + 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1]
	 - 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] - 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]
	*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 
	180*(1 + c[0] - 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*(153 + 153*c[0] - 
	94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) + 78*(33*(7*c[0]*(2 
	+ c[0]) - 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])
	*c[2] - 2948*c[2]*c[2] + 28*(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 
	1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] 
	+ 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-20 + 
	c[0]*c[0]*(3 + c[0])) + 252*(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0])
	 + 4*c[1]*c[1])*c[2] + 684*(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 
	3432*c[1]*(27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])
	*c[3]*c[3] + 6*(143*(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])
	*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*
	(7527*(1 + c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[2] = ((a - b)*(2*a*a*(13*(11* (21*(-100 + 5*c[0]*c[0]*(3 + 
	c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*
	c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]
	*(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0]-64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*
	c[3]*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]
	*c[2]) + 182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(-100 + 5*c[0]*c[0]*(3 + c[0])
	 + 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 
	14*(5*c[0]*(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49 + 
	49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]*(2+c[0])
	 + 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 
	28*c[2]*c[2])*c[3] + 44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 
	1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 
	52*c[1]*c[1]) + 44*(57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*
	(77*(1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] - 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 
	15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*
	(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*
	(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0]+64*c[1])
	*c[2] + 28*c[2]*c[2])* c[3] + 12*(550 + 550*c[0] + 517*c[1] - 328*c[2])*
	 c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) + 858*(1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]
	*c[2]) - 182*(77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(7*(-100 + 5*c[0]*c[0]*(3+c[0])
	 - 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) 
	- 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49 + 
	49*c[0] - 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2+c[0])
	 - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] - 2*c[1])*c[2] + 
	28*c[2]*c[2])*c[3] + 44*(153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]*c[3] + 
	1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*(2 + c[0]) - 52*(1 + c[0])*c[1] + 
	52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] + 28*
	(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] + 10274*c[3])* c[4]*c[4] - 
	18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-20 + c[0]*c[0]*(3 + c[0])) + 252*
	(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] + 684*
	(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(27 + 27*c[0] - 
	26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*c[0]
	*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] + 17628*c[2]*c[2] + 
	12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(1 + c[0]) - 8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[3] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]
	*c[3]) + 4*(11*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] 
	- 9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) 
	+ (45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])
	*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]
	*c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(5*pow(-2 + c[0],2)*(1 + c[0])
	 + 10*(-2 + c[0])*c[0]*c[1] + 20*(-1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) 
	- 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(-49
	 + 49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 22*(9*(7*(-2+c[0])
	*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3]
	 - 1416*c[3]*c[3]*c[3]) - 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])*c[1]
	 + 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*
	(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] -
	 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 
	15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*
	(-15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 
	+ c[0])*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 
	64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 
	328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] 
	+ 858*(-1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] 
	- 796*c[2]*c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*
	c[3]*c[3])*c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2]-5137*c[3])
	* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(140 + 35*c[0]*c[0]*c[0]
	 - 28*c[1]*c[1]*(5 + 2*c[1]) + 56*c[1]*(1 + c[1])*c[2] - 4*(49 + 22*c[1])
	*c[2]*c[2] - 72*c[2]*c[2]*c[2] - 35*c[0]*c[0]*(3 + 2*c[1] + 2*c[2]) + 
	28*c[0]*(5*c[1]*(1 + c[1]) + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2])) + 22*
	(63*c[0]*c[0] - 18*c[0]*(7 + 14*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + c[1]) + 
	45*(1 + 2*c[1])*c[2] + 7*c[2]*c[2]))*c[3] + 44*(-153 + 153*c[0] - 94*c[1] 
	- 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*c[0] + 52*
	c[1]*(1 + c[1]) - 2*c[0]*(7 + 26*c[1])) + 44*(-57 + 57*c[0] + 10*c[1])*c[2]
	 - 2948*c[2]*c[2] + 28*(77*(-1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] + 12*(22165*c[0] - 13*(1705 + 1094*c[1] + 1222*c[2]) + 
	10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*
	pow(-2 + c[0],2)*(1 + c[0]) + 252*(-1 + c[0])*c[1]*c[1] - 18*(21*(-2+c[0])
	*c[0] + 4*c[1]*c[1])*c[2] + 684*(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2])
	 - 3432*c[1]*(-27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*
	c[2])*c[3]*c[3] + 6*(143*(9*(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 
	+ c[0])*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 
	12*(7527*(-1 + c[0]) - 8930*c[2])*c[4]*c[4] + 
	7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[4] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*(-5 + c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*
	(11*(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 
	5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(-35*c[0]*c[0]*c[0] - 70*c[0]
	*c[0]*(c[1] - c[2]) - 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] +7*c[2]*c[2])
	 + 8*(-7*(10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2]-11*c[1]*
	c[2]*c[2] + 9*c[2]*c[2]*c[2])) + 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] + 4*
	(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*
	(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 
	1652*c[3]*c[3])*c[4] - 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*
	c[3]))* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] 
	+ 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*
	(7 + 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*
	(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(-15015*c[0]*c[0]*c[0] + 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) - 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(33*(7*
	(-10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2]
	 + 9*c[2]*c[2]*c[2]) - 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))
	*c[3] + 11*(47*c[1] + 59*c[2])*c[3]*c[3] - 177*c[3]*c[3]*c[3]) + 13*(-231 
	+ 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3]
	 + 413*c[3]*c[3])*c[4] + (7111*c[1] + 7943*c[2] - 5137*c[3])*c[4]*c[4] + 
	777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] 
	+ 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(-2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[5] = (429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*
	c[4] - 84*c[0]*(5*c[2] + c[4])) - 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] 
	+ 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 
	5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]
	*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*
	(15*c[0] - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*
	(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 
	54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4]+2316*c[4]*c[4])) 
	+ 1081080*((-a + b)/3. + ((a - b)* (15015*c[0]*c[0]*c[0] - 6006*c[0]*c[0]*
	(5*c[2] + c[4]) + 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] 
	+ 51*c[3]*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*(143*c[2]*(21*c[1]*c[1]
	 - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) - 13*(429*c[1]*c[1] - 737*
	c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] - 7943*c[2]*c[4]*c[4] - 
	777*c[4]*c[4]*c[4])))/360360.) - (a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 
	693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]
	*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*
	(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 
	34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] 
	- 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*
	c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 
	182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*
	c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])
	 + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3])
	 + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 
	45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*
	c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*
	(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]
	*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2])
	 + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3]
	 + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2]+5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*
	c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*
	c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-15015 -13*c[2]*
	 (99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*
	c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[6] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*(-5 + c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*
	(11*(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 
	5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(35*c[0]*c[0]*c[0] + 70*c[0]
	*c[0]*(c[1] - c[2]) + 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2]+7*c[2]*c[2])
	 + 8*(7*(10 - 5*c[1] + c[1]*c[1]*c[1]) + 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*
	c[2]*c[2] - 9*c[2]*c[2]*c[2])) - 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] + 4*
	(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*
	(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] + 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))*
	 c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*
	(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 
	+ 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*
	(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(15015*c[0]*c[0]*c[0] - 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) + 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(-33*(7*
	(-10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2]
	 + 9*c[2]*c[2]*c[2]) + 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))
	*c[3] - 11*(47*c[1] + 59*c[2])*c[3]*c[3] + 177*c[3]*c[3]*c[3]) - 13*(-231 
	+ 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3]
	 + 413*c[3]*c[3])*c[4] - (7111*c[1] + 7943*c[2] - 5137*c[3])*c[4]*c[4] - 
	777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] 
	+ 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(-2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[7] = (-429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*
	c[4] - 84*c[0]*(5*c[2] + c[4])) + 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*
	(15*c[0] + 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*
	(7*(-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])
	*c[2] + 60*c[2]*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]
	*c[3] - 44*(15*c[0] - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) 
	+ a*b*(33*(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]
	*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4] + 
	2316*c[4]*c[4])) + 1081080*((-a + b)/3. + ((a - b)* (15015*c[0]*c[0]*c[0] 
	- 6006*c[0]*c[0]*(5*c[2] + c[4]) + 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2]
	 - 42*c[1]*c[3] + 51*c[3]*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*
	(143*c[2]*(21*c[1]*c[1] - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) - 
	13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] - 
	7943*c[2]*c[4]*c[4] - 777*c[4]*c[4]*c[4])))/360360.) - (a - b)*(b*b*(13*
	(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]
	*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3]
	 + 50*c[3]*c[3]) + 4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 
	99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1]-c[2])
	*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*
	(165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*
	(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*
	c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] 
	- 13192*c[4]*c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*
	(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 
	45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(210 
	+ 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]
	*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])
	*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] 
	- 286*c[0]*(3*c[1] - 2*c[2]) + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*
	(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] 
	- 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + 
	a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*
	(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 583*c[3]*c[3] - 682*c[2]*
	c[4] + 579*c[4]*c[4]) + 8*(-15015 - 13*c[2]* (99*c[1]*c[1] + 517*c[2]*c[2] 
	- 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 
	3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*c[4] + 
	967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[8] = ((a - b)*(429*(-105*(16 + c[0]*c[0]*c[0] + 4*c[0]*
	(-3 + c[1]*c[1])) + 42*(5*c[0]*c[0] - 4*(5 + c[1]*c[1]))*c[2] - 588*c[0]*
	c[2]*c[2] + 216*c[2]*c[2]*c[2] + 72*c[1]*(7*c[0] - 10*c[2])*c[3] - 4*
	(153*c[0] - 118*c[2])*c[3]*c[3]) + 78*(11*(-84 + 21*c[0]*c[0] + 156*c[1]*
	c[1] + 228*c[0]*c[2] - 268*c[2]*c[2]) - 4312*c[1]*c[3] + 1652*c[3]*c[3])
	*c[4] - 156*(1705*c[0] - 1222*c[2])*c[4]*c[4] + 18648*c[4]*c[4]*c[4] + b*b*
	(13*(11*(105*c[0]*c[0]*c[0] + 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*
	c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 
	136*c[2]*c[2]*c[2] + 36*c[0]*(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*
	c[2])) - 33*(-252 + 9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 
	64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550*c[0] + 517*c[1] - 328*c[2])*
	c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*
	c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*
	(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])* c[4] + 12*(7319*c[0]
	 + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + 
	a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 
	132*c[0]*(3*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 15*c[2]*c[2]) - 9*(4*c[1] 
	+ 5*c[2])*c[3] + 50*c[3]*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*(-5 + c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 
	33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])
	*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(11*(-60 + 15*c[0]*c[0] - 78*c[0]
	*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 182*
	(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])* c[4] + 12*(7319*c[0]
	 - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + 
	a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*(11*
	(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 53*c[3]*c[3]) - 682*
	c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(-2310 + c[2]* (-2079 + 99*c[1]*c[1] 
	+ 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 3*(1287 + 1573*c[1]*
	c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])*c[4] - 13395*c[2]
	*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/1081080.;
		break;
	case 102:
		bf_mom[0] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]
	*c[3]) + 4*(11*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]
	*c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(-5*pow(-2 + c[0],2)*(1+c[0])
	 - 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1])
	 + 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49
	 + 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*(-2 + c[0])
	*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] - 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3]
	 + 1416*c[3]*c[3]*c[3]) + 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])*c[1]
	 + 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*
	(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] + 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 15*
	(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*
	(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 + 
	15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])
	*c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*
	 c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] + 858*(-1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]*
	c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4])) + a*(-45045*c[0]*c[0]*c[0] + 9009*c[0]*c[0]*
	(15 + 10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]) - 1716*c[0]*(21*(5*c[1]*(1+c[1])
	 + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2]) - 9*(7 + 14*c[1] + 10*c[2])*c[3] + 
	153*c[3]*c[3] + (21 + 78*c[1] - 114*c[2] - 98*c[3])*c[4] + 155*c[4]*c[4]) 
	+ 12*(13*(33*(-35 + 7*c[1]*c[1]*(5 + 2*c[1]) - 14*c[1]*(1 + c[1])*c[2] + 
	(49 + 22*c[1])*c[2]*c[2] + 18*c[2]*c[2]*c[2]) - 22*(9*c[1]*(7 + c[1]) + 
	45*(1 + 2*c[1])*c[2] + 7*c[2]*c[2])*c[3] + 11*(153 + 94*c[1] + 118*c[2])
	*c[3]*c[3] - 354*c[3]*c[3]*c[3]) + 26*(11*(39*c[1]*(1 + c[1]) + (-57 + 10*
	c[1])*c[2] - 67*c[2]*c[2]) - 7*(77 + 154*c[1] - 18*c[2])*c[3] + 413*c[3]
	*c[3])*c[4] + (22165 + 14222*c[1] + 15886*c[2] - 10274*c[3])*c[4]*c[4] + 
	1554*c[4]*c[4]*c[4]) + 2*b*(143*(105*pow(-2 + c[0],2)*(1 + c[0]) + 252*
	(-1 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 684*
	(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] - 
	26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*
	(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2] 
	+ 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(-1 + c[0])-8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[1] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]
	*c[3]) + 4*(11*(105 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 +
	 c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]
	*c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(5*pow(-2 + c[0],2)*(1 + c[0])
	 + 10*(-2 + c[0])*c[0]*c[1] + 20*(-1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1])
	 - 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*
	(-49 + 49*c[0] + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 22*(9*(7*(-2 + 
	c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]
	*c[3] - 1416*c[3]*c[3]*c[3]) - 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])
	*c[1] + 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] 
	- 28*(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] - 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*pow(-2 + c[0],2)*(1 + c[0]) + 
	15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 
	+ 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])
	*c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])* 
	c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] + 858*(-1 + 
	c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]
	*c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])
	*c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]
	*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(140 + 35*c[0]*c[0]*c[0] - 28*
	c[1]*c[1]*(5 + 2*c[1]) + 56*c[1]*(1 + c[1])*c[2] - 4*(49 + 22*c[1])*c[2]
	*c[2] - 72*c[2]*c[2]*c[2] - 35*c[0]*c[0]*(3 + 2*c[1] + 2*c[2]) + 28*c[0]*
	(5*c[1]*(1 + c[1]) + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2])) + 22*(63*c[0]*c[0]
	 - 18*c[0]*(7 + 14*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + c[1]) + 45*(1 + 2*c[1])
	*c[2] + 7*c[2]*c[2]))*c[3] + 44*(-153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]
	*c[3] + 1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*c[0] + 52*c[1]*(1 + c[1]) - 
	2*c[0]*(7 + 26*c[1])) + 44*(-57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] 
	+ 28*(77*(-1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(22165*c[0] - 13*(1705 + 1094*c[1] + 1222*c[2]) + 10274*c[3])* c[4]*c[4] 
	- 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*pow(-2 + c[0],2)*(1 + c[0]) + 252*
	(-1 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 684*
	(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] - 
	26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*(-2 
	+ c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2] + 
	12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(-1 + c[0]) - 8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[2] = ((a - b)*(2*a*a*(13*(11* (21*(-100 + 5*c[0]*c[0]*(3 + 
	c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*
	c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*
	(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]
	*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])
	*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 
	182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4]) + b*(39*(33*(7*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 10*c[0]*
	(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49 + 49*c[0] + 22*
	c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + 
	c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*
	c[3] + 44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*
	c[3]) - 78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 
	44*(57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 
	2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(22165 + 22165*c[0] + 
	14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 
	2*b*(13*(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 
	24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2]
	 - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*
	(550 + 550*c[0] + 517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(165*c[0]*(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(13 
	+ 13*c[0] - 5*c[1])*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1] + 
	18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] - 
	3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(7*
	(-100 + 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])
	*c[1]*c[1] - 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] 
	- 4*c[1]*c[1])*c[2] + 4*(49 + 49*c[0] - 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]
	*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 
	180*(1 + c[0] - 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*(153 + 153*c[0] - 
	94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*(2 + 
	c[0]) - 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])
	*c[2] - 2948*c[2]*c[2] + 28*(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 
	1652*c[3]*c[3])*c[4] + 12*(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] + 
	10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-20 + c[0]
	*c[0]*(3 + c[0])) + 252*(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0]) + 
	4*c[1]*c[1])*c[2] + 684*(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 
	3432*c[1]*(27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])
	*c[3]*c[3] + 6*(143*(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])
	*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*
	(7527*(1 + c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[3] = ((a - b)*(2*a*a*(13*(11* (21*(-100 + 5*c[0]*c[0]*(3 + 
	c[0]) - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*
	c[1]) - 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 
	36*(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]
	*(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]
	*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])
	*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 
	182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4]) + b*(39*(33*(7*(100 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 
	+ c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] + 22*
	c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + 
	c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*
	c[3] - 44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]
	*c[3]) + 78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 
	44*(57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 
	2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] + 
	14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 
	2*b*(13*(11*(21*(-100 + 5*c[0]*c[0]*(3 + c[0]) + 15*c[0]*(2 + c[0])*c[1] + 
	24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2]
	 - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*
	(550 + 550*c[0] + 517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(165*c[0]*(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(13 
	+ 13*c[0] - 5*c[1])*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1] +
	 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] -
	 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*
	(7*(100 - 5*c[0]*c[0]*(3 + c[0]) + 10*c[0]*(2 + c[0])*c[1] - 20*(1 + c[0])
	*c[1]*c[1] + 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] 
	- 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] - 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]
	*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 
	180*(1 + c[0] - 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*(153 + 153*c[0] - 
	94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) + 78*(33*(7*c[0]*(2 + 
	c[0]) - 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])
	*c[2] - 2948*c[2]*c[2] + 28*(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 
	1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] + 
	10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-20 + c[0]
	*c[0]*(3 + c[0])) + 252*(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0]) + 
	4*c[1]*c[1])*c[2] + 684*(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 
	3432*c[1]*(27 + 27*c[0] - 26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])
	*c[3]*c[3] + 6*(143*(9*c[0]*(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])
	*c[2] + 17628*c[2]*c[2] + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*
	(7527*(1 + c[0]) - 8930*c[2])*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[4] = (-429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]
	*c[4] - 84*c[0]*(5*c[2] + c[4])) + 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] 
	+ 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 
	5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]
	*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0]
	 - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*c[0]*
	c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*
	c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4])) + 
	1081080*((-a + b)/3. + ((a - b)* (15015*c[0]*c[0]*c[0] - 6006*c[0]*c[0]*
	(5*c[2] + c[4]) + 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] 
	+ 51*c[3]*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*(143*c[2]*(21*c[1]*c[1]
	 - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) - 13*(429*c[1]*c[1] - 737*
	c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] - 7943*c[2]*c[4]*c[4] - 
	777*c[4]*c[4]*c[4])))/360360.) - (a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 
	693*c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*
	c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 
	4*(11*(-210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 
	34*c[2]*c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] 
	- 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*
	c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 
	182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*
	c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])
	 + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) 
	+ 132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 
	45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*
	c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*
	(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*
	c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2])
	 + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3]
	 + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])
	* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]
	*c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]
	*c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-15015 - 13*
	c[2]* (99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*
	(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[5] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]
	*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*(-5+c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*
	(11*(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 
	5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(35*c[0]*c[0]*c[0] + 70*c[0]
	*c[0]*(c[1] - c[2]) + 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 7*c[2]*
	c[2]) + 8*(7*(10 - 5*c[1] + c[1]*c[1]*c[1]) + 7*(5 + c[1]*c[1])*c[2] + 
	11*c[1]*c[2]*c[2] - 9*c[2]*c[2]*c[2])) - 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] 
	+ 4*(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 
	44*(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] + 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))*
	 c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*
	(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 
	4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*
	(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(15015*c[0]*c[0]*c[0] - 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) + 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(-33*(7*
	(-10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2]
	 + 9*c[2]*c[2]*c[2]) + 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))
	*c[3] - 11*(47*c[1] + 59*c[2])*c[3]*c[3] + 177*c[3]*c[3]*c[3]) - 13*(-231 
	+ 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3]
	 + 413*c[3]*c[3])*c[4] - (7111*c[1] + 7943*c[2] - 5137*c[3])*c[4]*c[4] - 
	777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] 
	+ 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(-2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[6] = (429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*
	c[4] - 84*c[0]*(5*c[2] + c[4])) - 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] 
	+ 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 
	5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*
	c[2]*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*
	(15*c[0] - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*
	(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 
	54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4]))
	 + 1081080*((-a + b)/3. + ((a - b)* (15015*c[0]*c[0]*c[0] - 6006*c[0]*c[0]*
	(5*c[2] + c[4]) + 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] 
	+ 51*c[3]*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*(143*c[2]*(21*c[1]*c[1]
	 - 27*c[2]*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3]) - 13*(429*c[1]*c[1]-737*c[2]
	*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] - 7943*c[2]*c[4]*c[4] - 777*
	c[4]*c[4]*c[4])))/360360.) - (a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*
	c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2]
	 + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(-210
	 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]
	*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])
	*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] 
	- 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*
	c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 
	7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + 
	a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 
	132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*
	c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]
	*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] 
	+ 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3])) 
	- 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2]) + 
	220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] 
	+ 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])
	* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]
	*c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]
	*c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-15015 - 13*
	c[2]* (99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*
	(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 
	13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[7] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*(-5 + c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] 
	+ 12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*
	(11*(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 
	5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] 
	+ 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])
	* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(-35*c[0]*c[0]*c[0] - 
	70*c[0]*c[0]*(c[1] - c[2]) - 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 
	7*c[2]*c[2]) + 8*(-7*(10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2]
	 - 11*c[1]*c[2]*c[2] + 9*c[2]*c[2]*c[2])) + 22*(9*(7*c[0]*c[0] + 28*c[0]*
	c[1] + 4*(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] 
	- 44*(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] - 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))*
	 c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*
	(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 
	4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*(-35
	 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(-15015*c[0]*c[0]*c[0] + 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) - 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(33*(7*
	(-10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*
	c[2] + 9*c[2]*c[2]*c[2]) - 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*
	c[2]))*c[3] + 11*(47*c[1] + 59*c[2])*c[3]*c[3] - 177*c[3]*c[3]*c[3]) + 13*
	(-231 + 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 
	126*c[2]*c[3] + 413*c[3]*c[3])*c[4] + (7111*c[1] + 7943*c[2] - 5137*c[3])
	*c[4]*c[4] + 777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]
	*c[0]*(7*c[2] - c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 
	54*c[1]*c[3] + 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*
	(-2310 + c[2]* (-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 
	963*c[3]*c[3])) - 3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*
	c[3] + 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*c[4] + 
	967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[8] = ((a - b)*(429*(-105*(16 + c[0]*c[0]*c[0] + 4*c[0]*
	(-3 + c[1]*c[1])) + 42*(5*c[0]*c[0] - 4*(5 + c[1]*c[1]))*c[2] - 588*c[0]
	*c[2]*c[2] + 216*c[2]*c[2]*c[2] + 72*c[1]*(7*c[0] - 10*c[2])*c[3] - 4*
	(153*c[0] - 118*c[2])*c[3]*c[3]) + 78*(11*(-84 + 21*c[0]*c[0] + 156*c[1]
	*c[1] + 228*c[0]*c[2] - 268*c[2]*c[2]) - 4312*c[1]*c[3] + 1652*c[3]*c[3])
	*c[4] - 156*(1705*c[0] - 1222*c[2])*c[4]*c[4] + 18648*c[4]*c[4]*c[4] + b*b*
	(13*(11*(105*c[0]*c[0]*c[0] + 84*(20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]
	*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 
	136*c[2]*c[2]*c[2] + 36*c[0]*(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]
	*c[2])) - 33*(-252 + 9*(c[0] + 2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 
	64*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 12*(550*c[0] + 517*c[1] - 328*c[2])
	*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(-660 + 165*c[0]*c[0] + 616*
	c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 
	182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])* c[4] + 12*
	(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*
	c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2]
	 - 3*c[3]) + 132*c[0]*(3*(-35 + 14*c[1]*c[1] - 7*c[1]*c[2] + 15*c[2]*c[2])
	 - 9*(4*c[1] + 5*c[2])*c[3] + 50*c[3]*c[3]) + 4*(-11*(21*(-20 + 3*c[1]*
	(-5 + c[1]*c[1])) - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*
	c[2]*c[2]*c[2]) + 33*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*
	(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(11*(-60 + 
	15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 5*c[1])*c[2] - 
	796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	 c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2]
	 - c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3]
	 + 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(-2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/1081080.;
		break;
	case 103:
		bf_mom[0] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*
	c[3]) + 4*(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 + 
	c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*
	c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*
	c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(7*(-5*(20 + (-3+c[0])*c[0]*c[0]) 
	- 10*(-2 + c[0])*c[0]*c[1] - 20*(-1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) 
	+ 14*(5*(-2 + c[0])*c[0] - 4*(-1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(-49 
	+ 49*c[0] + 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*(-2 + c[0])
	*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])*c[2] 
	+ 28*c[2]*c[2])*c[3] - 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] 
	+ 1416*c[3]*c[3]*c[3]) + 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])*c[1] 
	+ 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*
	(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] + 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 
	15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 
	+ 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*
	c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*
	 c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] + 858*(-1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]*
	c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4]
	+ 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4])) + a*(-45045*c[0]*c[0]*c[0] + 9009*c[0]*c[0]*(15 + 
	10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]) - 1716*c[0]*(21*(5*c[1]*(1 + c[1]) + 
	5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2]) - 9*(7 + 14*c[1] + 10*c[2])*c[3] + 
	153*c[3]*c[3] + (21 + 78*c[1] - 114*c[2] - 98*c[3])*c[4] + 155*c[4]*c[4]) 
	+ 12*(13*(33*(-175 + 7*c[1]*c[1]*(5 + 2*c[1]) - 14*c[1]*(1 + c[1])*c[2] + 
	(49 + 22*c[1])*c[2]*c[2] + 18*c[2]*c[2]*c[2]) - 22*(9*c[1]*(7 + c[1]) + 
	45*(1 + 2*c[1])*c[2] + 7*c[2]*c[2])*c[3] + 11*(153 + 94*c[1] + 118*c[2])*
	c[3]*c[3] - 354*c[3]*c[3]*c[3]) + 26*(11*(39*c[1]*(1 + c[1]) + (-57 + 10*
	c[1])*c[2] - 67*c[2]*c[2]) - 7*(77 + 154*c[1] - 18*c[2])*c[3] + 413*c[3]*
	c[3])*c[4] + (22165 + 14222*c[1] + 15886*c[2] - 10274*c[3])*c[4]*c[4] + 
	1554*c[4]*c[4]*c[4]) + 2*b*(143*(105*(20 + (-3 + c[0])*c[0]*c[0]) + 252*(-1
	 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 684*(-1 
	+ c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] - 26*
	c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*(-2 + 
	c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2] + 
	12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(-1 + c[0]) - 8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[1] = -((a - b)*(2*a*a*(13*(11*(21*(5*(-1+c[0])*pow(2+c[0],2)
	 - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) 
	- 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 
	+ 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*(2 + 
	c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])*c[2]
	 + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]*c[3]
	 + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])*c[1] 
	+ 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 182*
	(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 
	+ 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*
	c[4]*c[4]) + b*(39*(33*(7*(20 - 5*c[0]*c[0]*(3 + c[0]) - 10*c[0]*(2 + c[0])
	*c[1] - 20*(1 + c[0])*c[1]*c[1] - 8*c[1]*c[1]*c[1]) + 14*(5*c[0]*(2 + c[0])
	 - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] - 4*(49 + 49*c[0] + 22*c[1])*c[2]*
	c[2] + 72*c[2]*c[2]*c[2]) + 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + c[0])*c[1] + 
	4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 44*(153 
	+ 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*(33*
	(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] 
	- 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 2*c[1]) + 18*c[2])*
	c[3] + 1652*c[3]*c[3])*c[4] - 12*(22165 + 22165*c[0] + 14222*c[1] - 15886*
	c[2] - 10274*c[3])* c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(5*
	(-1 + c[0])*pow(2 + c[0],2) + 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*
	c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*(1 + c[0])*c[1] - 
	16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*
	c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*
	(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 12*(550 + 550*c[0] + 
	517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*
	(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0]-5*c[1])
	*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137
	*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(7*(20 - 5*c[0]*c[0]
	*(3 + c[0]) + 10*c[0]*(2 + c[0])*c[1] - 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*
	c[1]*c[1]) + 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] 
	- 4*(49 + 49*c[0] - 22*c[1])*c[2]*c[2] + 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]
	*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] - 2*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] - 44*(153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]
	*c[3] - 1416*c[3]*c[3]*c[3]) + 78*(33*(7*c[0]*(2 + c[0]) - 52*(1 + c[0])
	*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] 
	+ 28*(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] - 12*
	(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] + 10274*c[3])* c[4]*c[4] + 
	18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-1 + c[0])*pow(2 + c[0],2) + 252*
	(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] + 684*
	(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(27 + 27*c[0] - 26*
	c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*c[0]*
	(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] + 17628*c[2]*c[2] + 
	12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(1 + c[0]) - 8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[2] = -((a - b)*(2*a*a*(13*(11*(21*(5*(-1+c[0])*pow(2+c[0],2)
	 - 15*c[0]*(2 + c[0])*c[1] + 24*(1 + c[0])*c[1]*c[1] - 12*c[1]*c[1]*c[1]) 
	- 18*(7*c[0]*(2 + c[0]) + 14*(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*
	(15 + 15*c[0] - 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) + 33*(9*(7*c[0]*
	(2 + c[0]) - 16*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(45 + 45*c[0] - 64*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 12*(550 + 550*c[0] - 517*c[1] - 328*c[2])*c[3]
	*c[3] + 2124*c[3]*c[3]*c[3]) - 6*(13*(165*c[0]*(2 + c[0]) - 858*(1 + c[0])
	*c[1] + 616*c[1]*c[1] + 44*(13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*c[2]) + 
	182*(77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*
	(7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*
	c[4]*c[4]*c[4]) + b*(39*(33*(7*(5*(-1 + c[0])*pow(2 + c[0],2) + 10*c[0]*
	(2 + c[0])*c[1] + 20*(1 + c[0])*c[1]*c[1] + 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*
	(2 + c[0]) - 4*(1 + c[0])*c[1] - 4*c[1]*c[1])*c[2] + 4*(49 + 49*c[0] + 22*
	c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2]) - 22*(9*(7*c[0]*(2 + c[0]) + 28*(1 + 
	c[0])*c[1] + 4*c[1]*c[1]) - 180*(1 + c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])
	*c[3] + 44*(153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*
	c[3]*c[3]) - 78*(33*(7*c[0]*(2 + c[0]) + 52*(1 + c[0])*c[1] + 52*c[1]*c[1])
	 + 44*(57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] - 28*(77*(1 + c[0] + 
	2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*(22165 + 22165*c[0] + 
	14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 
	2*b*(13*(11*(21*(5*(-1 + c[0])*pow(2 + c[0],2) + 15*c[0]*(2 + c[0])*c[1] + 
	24*(1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 18*(7*c[0]*(2 + c[0]) - 14*
	(1 + c[0])*c[1] - 16*c[1]*c[1])*c[2] + 36*(15 + 15*c[0] + 11*c[1])*c[2]*
	c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*c[0]*(2 + c[0]) + 16*(1 + c[0])*c[1] 
	+ 4*c[1]*c[1]) - 4*(45 + 45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])* c[3] + 
	12*(550 + 550*c[0] + 517*c[1] - 328*c[2])* c[3]*c[3] - 2124*c[3]*c[3]*c[3])
	 - 6*(13*(165*c[0]*(2 + c[0]) + 858*(1 + c[0])*c[1] + 616*c[1]*c[1] + 44*
	(13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]*c[2]) - 182*(77 + 77*c[0] + 120*c[1]
	 + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319 + 7319*c[0] + 7111*c[1] 
	- 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(39*(33*
	(7*(5*(-1 + c[0])*pow(2 + c[0],2) - 10*c[0]*(2 + c[0])*c[1] + 20*(1 + c[0])
	*c[1]*c[1] - 8*c[1]*c[1]*c[1]) - 14*(5*c[0]*(2 + c[0]) + 4*(1 + c[0])*c[1] 
	- 4*c[1]*c[1])*c[2] + 4*(49 + 49*c[0] - 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*
	c[2]) + 22*(9*(7*c[0]*(2 + c[0]) - 28*(1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*
	(1 + c[0] - 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*(153 + 153*c[0]-94*c[1] 
	- 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*(2 + c[0]) - 
	52*(1 + c[0])*c[1] + 52*c[1]*c[1]) + 44*(57 + 57*c[0] + 10*c[1])*c[2] - 
	2948*c[2]*c[2] + 28*(77*(1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*
	c[3])*c[4] + 12*(22165 + 22165*c[0] - 14222*c[1] - 15886*c[2] + 10274*c[3])
	*c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(-1 + c[0])*pow(2+c[0],2)
	 + 252*(1 + c[0])*c[1]*c[1] - 18*(21*c[0]*(2 + c[0]) + 4*c[1]*c[1])*c[2] + 
	684*(1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(27 + 27*c[0] - 
	26*c[2])*c[3] + 156*(583*(1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*c[0]*
	(2 + c[0]) - 44*c[1]*c[1]) - 17732*(1 + c[0])*c[2] + 17628*c[2]*c[2] + 
	12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(1 + c[0]) - 8930*c[2])*
	c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[3] = -((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + 66*c[0]*(21*c[1]*(5 + 4*c[1]) - 42*
	(-1 + c[1])*c[2] + 90*c[2]*c[2] - 9*(7 + 8*c[1] + 10*c[2])*c[3] + 100*c[3]*
	c[3]) + 4*(11*(525 - 63*c[1]*c[1]*(2 + c[1]) + 9*c[1]*(7 + 8*c[1])*c[2] - 
	9*(15 + 11*c[1])*c[2]*c[2] - 34*c[2]*c[2]*c[2]) + 33*(9*c[1]*(4 + c[1]) + 
	(45 + 64*c[1])*c[2] + 7*c[2]*c[2])*c[3] - 3*(550 + 517*c[1] + 328*c[2])*
	c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*(-2 + c[0])*c[0] - 858*(-1 +
	 c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] + 5*c[1])*c[2] - 796*c[2]*
	c[2]) + 182*(-77 + 77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	c[4] + 12*(-7319 + 7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*
	c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(700 + 35*c[0]*c[0]*c[0] + 28*
	c[1]*c[1]*(-5 + 2*c[1]) + 35*c[0]*c[0]*(-3 + 2*c[1] - 2*c[2]) + 56*(-1 + 
	c[1])*c[1]*c[2] + 4*(-49 + 22*c[1])*c[2]*c[2] - 72*c[2]*c[2]*c[2] + 28*c[0]
	*(5*(-1 + c[1])*c[1] + (5 + 2*c[1])*c[2] + 7*c[2]*c[2])) - 22*(9*(7*(-2 + 
	c[0])*c[0] + 28*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 180*(-1 + c[0] + 2*c[1])
	*c[2] + 28*c[2]*c[2])*c[3] + 44*(-153 + 153*c[0] + 94*c[1] - 118*c[2])*c[3]
	*c[3] - 1416*c[3]*c[3]*c[3]) - 78*(33*(7*(-2 + c[0])*c[0] + 52*(-1 + c[0])
	*c[1] + 52*c[1]*c[1]) + 44*(-57 + 57*c[0] - 10*c[1])*c[2] - 2948*c[2]*c[2] 
	- 28*(77*(-1 + c[0] + 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(-22165 + 22165*c[0] + 14222*c[1] - 15886*c[2] - 10274*c[3])* c[4]*c[4] - 
	18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(21*(100 + 5*(-3 + c[0])*c[0]*c[0] + 
	15*(-2 + c[0])*c[0]*c[1] + 24*(-1 + c[0])*c[1]*c[1] + 12*c[1]*c[1]*c[1]) - 
	18*(7*c[0]*c[0] + 2*(7 - 8*c[1])*c[1] - 14*c[0]*(1 + c[1]))*c[2] + 36*(-15 
	+ 15*c[0] + 11*c[1])*c[2]*c[2] - 136*c[2]*c[2]*c[2]) - 33*(9*(7*(-2 + c[0])
	*c[0] + 16*(-1 + c[0])*c[1] + 4*c[1]*c[1]) - 4*(-45 + 45*c[0] + 64*c[1])*
	c[2] + 28*c[2]*c[2])* c[3] + 12*(-550 + 550*c[0] + 517*c[1] - 328*c[2])*
	 c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 6*(13*(165*(-2 + c[0])*c[0] + 858*(-1 
	+ c[0])*c[1] + 616*c[1]*c[1] + 44*(-13 + 13*c[0] - 5*c[1])*c[2] - 796*c[2]
	*c[2]) - 182*(-77 + 77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*
	c[4] + 12*(-7319 + 7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4]
	 - 13192*c[4]*c[4]*c[4])) + a*(39*(33*(700 + 35*c[0]*c[0]*c[0] - 28*c[1]*
	c[1]*(5 + 2*c[1]) + 56*c[1]*(1 + c[1])*c[2] - 4*(49 + 22*c[1])*c[2]*c[2] - 
	72*c[2]*c[2]*c[2] - 35*c[0]*c[0]*(3 + 2*c[1] + 2*c[2]) + 28*c[0]*(5*c[1]*
	(1 + c[1]) + 5*c[2] - 2*c[1]*c[2] + 7*c[2]*c[2])) + 22*(63*c[0]*c[0] - 18*
	c[0]*(7 + 14*c[1] + 10*c[2]) + 4*(9*c[1]*(7 + c[1]) + 45*(1 + 2*c[1])*c[2] 
	+ 7*c[2]*c[2]))*c[3] + 44*(-153 + 153*c[0] - 94*c[1] - 118*c[2])*c[3]*c[3] 
	+ 1416*c[3]*c[3]*c[3]) - 78*(33*(7*c[0]*c[0] + 52*c[1]*(1 + c[1]) - 2*c[0]*
	(7 + 26*c[1])) + 44*(-57 + 57*c[0] + 10*c[1])*c[2] - 2948*c[2]*c[2] + 28*
	(77*(-1 + c[0] - 2*c[1]) + 18*c[2])*c[3] + 1652*c[3]*c[3])*c[4] + 12*
	(22165*c[0] - 13*(1705 + 1094*c[1] + 1222*c[2]) + 10274*c[3])* c[4]*c[4] - 
	18648*c[4]*c[4]*c[4] + 2*b*(143*(105*(20 + (-3 + c[0])*c[0]*c[0]) + 252*
	(-1 + c[0])*c[1]*c[1] - 18*(21*(-2 + c[0])*c[0] + 4*c[1]*c[1])*c[2] + 684*
	(-1 + c[0])*c[2]*c[2] - 376*c[2]*c[2]*c[2]) - 3432*c[1]*(-27 + 27*c[0] - 
	26*c[2])*c[3] + 156*(583*(-1 + c[0]) - 642*c[2])*c[3]*c[3] + 6*(143*(9*
	(-2 + c[0])*c[0] - 44*c[1]*c[1]) - 17732*(-1 + c[0])*c[2] + 17628*c[2]*c[2]
	 + 12376*c[1]*c[3] - 6420*c[3]*c[3])*c[4] + 12*(7527*(-1 + c[0])-8930*c[2])
	*c[4]*c[4] + 7736*c[4]*c[4]*c[4]))))/8648640.;
		bf_mom[4] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1])) 
	- 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 
	12*c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*
	(11*(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 
	5*c[1])*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])* c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(-35*c[0]*c[0]*c[0] - 70*
	c[0]*c[0]*(c[1] - c[2]) - 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 7*
	c[2]*c[2]) + 8*(70 + 35*c[1] - 7*c[1]*c[1]*c[1] - 7*(5 + c[1]*c[1])*c[2] 
	- 11*c[1]*c[2]*c[2] + 9*c[2]*c[2]*c[2])) + 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1]
	 + 4*(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] - 
	44*(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] + 1416*c[3]*c[3]*c[3]) + 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] - 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))*
	 c[4]*c[4] + 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*
	(-20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 
	+ 4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*
	(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(-15015*c[0]*c[0]*c[0] + 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) - 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(33*(7*(10 
	- 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2] + 
	9*c[2]*c[2]*c[2]) - 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))*
	c[3] + 11*(47*c[1] + 59*c[2])*c[3]*c[3] - 177*c[3]*c[3]*c[3]) + 13*(-231 + 
	429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3] 
	+ 413*c[3]*c[3])*c[4] + (7111*c[1] + 7943*c[2] - 5137*c[3])*c[4]*c[4] + 
	777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] 
	+ 53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[5] = (-429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*
	c[4] - 84*c[0]*(5*c[2] + c[4])) + 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] 
	+ 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 
	5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]
	*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0]
	 - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*c[0]*
	c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3]
	 + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4])) + 360360*
	(-a + b + ((a - b)* (-15015*c[0]*c[0]*c[0] + 6006*c[0]*c[0]*(5*c[2] + c[4]) 
	- 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3]) 
	- 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*(143*c[2]*(-21*c[1]*c[1] + 27*c[2]*c[2]
	 - 90*c[1]*c[3] + 59*c[3]*c[3]) + 13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*
	c[1]*c[3] + 413*c[3]*c[3])*c[4] + 7943*c[2]*c[4]*c[4] + 777*c[4]*c[4]*
	c[4])))/120120.) + (a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*
	(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*
	c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*
	c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) 
	- 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*
	c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 220*
	c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 
	120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] 
	- 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*a*(13*
	(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*
	(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 
	50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*
	c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*
	c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*
	(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2]) + 220*c[1]*c[2]
	 - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3]+7528*c[3]*c[3])
	*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 
	13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] 
	- c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 583*c[3]
	*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(15015 - 13*c[2]*(99*c[1]*c[1] 
	+ 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*c[1]*c[1] + 
	4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*
	c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[6] = ((a - b)*(2*a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*
	c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 4*(-11*(21*(20 + 3*c[1]*(-5 + c[1]*c[1]))
	 - 18*(7 + 4*c[1]*c[1])*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*
	(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2]))*c[3] - 3*(517*c[1] + 328*c[2])*c[3]
	*c[3] + 531*c[3]*c[3]*c[3]) + 132*c[0]*(42*c[1]*c[1] - 3*c[1]*(7*c[2] + 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] - 9*c[2]*c[3] + 10*c[3]*c[3]))) - 6*(13*(11*
	(-60 + 15*c[0]*c[0] - 78*c[0]*c[1] + 56*c[1]*c[1]) + 44*(13*c[0] + 5*c[1])
	*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]
	*c[3])* c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*
	c[4] - 13192*c[4]*c[4]*c[4]) + b*(39*(33*(35*c[0]*c[0]*c[0] + 70*c[0]*c[0]*
	(c[1] - c[2]) + 28*c[0]*(-15 + 5*c[1]*c[1] + 2*c[1]*c[2] + 7*c[2]*c[2]) + 
	8*(7*(-10 - 5*c[1] + c[1]*c[1]*c[1]) + 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*
	c[2]*c[2] - 9*c[2]*c[2]*c[2])) - 22*(9*(7*c[0]*c[0] + 28*c[0]*c[1] + 4*
	(-7 + c[1]*c[1])) - 180*(c[0] + 2*c[1])*c[2] + 28*c[2]*c[2])*c[3] + 44*
	(153*c[0] + 94*c[1] - 118*c[2])*c[3]*c[3] - 1416*c[3]*c[3]*c[3]) - 78*
	(231*c[0]*c[0] + 44*(-21 + 39*c[1]*c[1] - 10*c[1]*c[2] - 67*c[2]*c[2]) + 
	44*c[0]*(39*c[1] + 57*c[2] - 49*c[3]) - 56*(77*c[1] + 9*c[2])*c[3] + 1652*
	c[3]*c[3])*c[4] + 12*(22165*c[0] + 2*(7111*c[1] - 7943*c[2] - 5137*c[3]))*
	 c[4]*c[4] - 18648*c[4]*c[4]*c[4] + 2*b*(13*(11*(105*c[0]*c[0]*c[0] + 84*
	(-20 + 3*c[1]*(-5 + c[1]*c[1])) + 63*c[0]*c[0]*(5*c[1] - 2*c[2]) + 72*(7 + 
	4*c[1]*c[1])*c[2] + 396*c[1]*c[2]*c[2] - 136*c[2]*c[2]*c[2] + 36*c[0]*
	(-35 + 14*c[1]*c[1] + 7*c[1]*c[2] + 15*c[2]*c[2])) - 33*(-252 + 9*(c[0] + 
	2*c[1])*(7*c[0] + 2*c[1]) - 4*(45*c[0] + 64*c[1])*c[2] + 28*c[2]*c[2])*c[3]
	 + 12*(550*c[0] + 517*c[1] - 328*c[2])*c[3]*c[3] - 2124*c[3]*c[3]*c[3]) - 
	6*(13*(-660 + 165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*c[2]
	 + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 
	7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137*c[3])*
	 c[4]*c[4] - 13192*c[4]*c[4]*c[4])) + a*(3*(15015*c[0]*c[0]*c[0] - 6006*
	c[0]*c[0]*(5*c[1] + 5*c[2] - 3*c[3] + c[4]) + 572*c[0]*(3*(-105 + 35*c[1]*
	c[1] + 49*c[2]*c[2] - 30*c[2]*c[3] + 51*c[3]*c[3] - 14*c[1]*(c[2] + 3*c[3]))
	 + 2*(39*c[1] - 57*c[2] - 49*c[3])*c[4] + 155*c[4]*c[4]) + 8*(13*(-33*(7*
	(10 - 5*c[1] + c[1]*c[1]*c[1]) - 7*(5 + c[1]*c[1])*c[2] + 11*c[1]*c[2]*c[2]
	 + 9*c[2]*c[2]*c[2]) + 11*(9*c[1]*c[1] + 90*c[1]*c[2] + 7*(-9 + c[2]*c[2]))
	*c[3] - 11*(47*c[1] + 59*c[2])*c[3]*c[3] + 177*c[3]*c[3]*c[3]) - 13*(-231 
	+ 429*c[1]*c[1] - 737*c[2]*c[2] + 22*c[1]*(5*c[2] - 49*c[3]) + 126*c[2]*c[3]
	 + 413*c[3]*c[3])*c[4] - (7111*c[1] + 7943*c[2] - 5137*c[3])*c[4]*c[4] - 
	777*c[4]*c[4]*c[4])) + 2*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - 
	c[4]) + 156*c[0]*(11*(-105 + 21*c[1]*c[1] + 57*c[2]*c[2] - 54*c[1]*c[3] + 
	53*c[3]*c[3]) - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(-13*(2310 + c[2]* 
	(-2079 + 99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3])) - 
	3*(1287 + 1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3])
	*c[4] - 13395*c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4])))))/4324320.;
		bf_mom[7] = (429*(a - b)*(315*c[0]*c[0] + 420*(-3 + c[1]*c[1]) + 
	588*c[2]*c[2] - 504*c[1]*c[3] + 612*c[3]*c[3] - 456*c[2]*c[4] + 620*c[4]*
	c[4] - 84*c[0]*(5*c[2] + c[4])) - 39*(a - b)*(b*b*(1155*c[0]*c[0] + 462*
	c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*(7*c[2] - 12*
	c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 44*(15*c[0] 
	+ 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*a*(33*(7*(-20 + 
	5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]
	*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*
	c[0] - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*
	c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*
	c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4])) + 
	360360*(-a + b + ((a - b)* (-15015*c[0]*c[0]*c[0] + 6006*c[0]*c[0]*(5*c[2] 
	+ c[4]) - 572*c[0]*(3*(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]
	*c[3]) - 114*c[2]*c[4] + 155*c[4]*c[4]) + 8*(143*c[2]*(-21*c[1]*c[1] + 
	27*c[2]*c[2] - 90*c[1]*c[3] + 59*c[3]*c[3]) + 13*(429*c[1]*c[1] - 737*c[2]*
	c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3])*c[4] + 7943*c[2]*c[4]*c[4] + 777*
	c[4]*c[4]*c[4])))/120120.) + (a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*
	c[0]*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*
	c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*
	(210 + 63*c[1]*c[1]*c[1] + 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*
	c[2]*c[2]) - 33*(c[1] - 7*c[2])*(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*
	c[2])*c[3]*c[3] - 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*
	c[1] - 220*c[1]*c[2] - 796*c[2]*c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*
	(77*c[0] + 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] 
	+ 7111*c[1] - 3478*c[2] - 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + 
	a*a*(13*(1155*c[0]*c[0]*c[0] - 693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 
	132*c[0]*(42*c[1]*c[1] - 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*
	c[2]*c[3] + 50*c[3]*c[3]) + 4*(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*
	c[1]*c[2] + 99*c[1]*c[2]*c[2] + 34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*
	(c[1] + 7*c[2])*c[3] - 3*(517*c[1] + 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*
	c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2])
	 + 220*c[1]*c[2] - 796*c[2]*c[2]) + 182*(77*c[0] - 120*c[1] + 18*c[2])*c[3]
	 + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] - 7111*c[1] - 3478*c[2]+5137*c[3])
	*c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*
	c[0]*(7*c[2] - c[4]) + 156*c[0]*(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*
	c[3] + 583*c[3]*c[3] - 682*c[2]*c[4] + 579*c[4]*c[4]) + 8*(15015 - 13*c[2]*
	(99*c[1]*c[1] + 517*c[2]*c[2] - 858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*
	c[1]*c[1] + 4407*c[2]*c[2] + 3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 13395*
	c[2]*c[4]*c[4] + 967*c[4]*c[4]*c[4]))))/2162160.;
		bf_mom[8] = -a + b + ((-a + b)*c[0])/2. + ((a - b)*c[2])/3. + 
	((a - b)*c[4])/15. + ((a - b)*(b*b*(7*(10 + 5*c[0] + 5*c[1] - 2*c[2] - 3*
	c[3]) - 10*c[4]) + a*a* (7*(10 + 5*c[0] - 5*c[1] - 2*c[2] + 3*c[3]) - 
	10*c[4]) + a*b*(70 + 35*c[0] - 42*c[2] + 6*c[4])))/210. + (a - b + ((a - b)
	*(15015*c[0]*c[0]*c[0] - 6006*c[0]*c[0]*(5*c[2] + c[4]) + 572*c[0]*(3*
	(35*c[1]*c[1] + 49*c[2]*c[2] - 42*c[1]*c[3] + 51*c[3]*c[3]) - 114*c[2]*c[4]
	 + 155*c[4]*c[4]) + 8*(143*c[2]*(21*c[1]*c[1] - 27*c[2]*c[2] + 90*c[1]*c[3]
	 - 59*c[3]*c[3]) - 13*(429*c[1]*c[1] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 
	413*c[3]*c[3])*c[4] - 7943*c[2]*c[4]*c[4]-777*c[4]*c[4]*c[4])))/120120.)/3.
	 - ((a - b)*(b*b*(13*(1155*c[0]*c[0]*c[0] + 693*c[0]*c[0]*(5*c[1] - 2*c[2] 
	- 3*c[3]) + 132*c[0]*(42*c[1]*c[1] + 21*c[1]*c[2] + 45*c[2]*c[2] - 36*c[1]*
	c[3] + 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*(11*(210 + 63*c[1]*c[1]*c[1] + 72*
	c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] - 34*c[2]*c[2]*c[2]) - 33*(c[1]-7*c[2])*
	(9*c[1] - c[2])*c[3] + 3*(517*c[1] - 328*c[2])*c[3]*c[3] - 531*c[3]*c[3]*
	c[3])) - 6*(13*(165*c[0]*c[0] + 616*c[1]*c[1] - 220*c[1]*c[2] - 796*c[2]*
	c[2] + 286*c[0]*(3*c[1] + 2*c[2])) - 182*(77*c[0] + 120*c[1] + 18*c[2])*
	c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*c[0] + 7111*c[1] - 3478*c[2] - 5137
	*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4]) + a*a*(13*(1155*c[0]*c[0]*c[0] - 
	693*c[0]*c[0]*(5*c[1] + 2*c[2] - 3*c[3]) + 132*c[0]*(42*c[1]*c[1] - 21*c[1]
	*c[2] + 45*c[2]*c[2] - 36*c[1]*c[3] - 45*c[2]*c[3] + 50*c[3]*c[3]) + 4*
	(-11*(-210 + 63*c[1]*c[1]*c[1] - 72*c[1]*c[1]*c[2] + 99*c[1]*c[2]*c[2] + 
	34*c[2]*c[2]*c[2]) + 33*(9*c[1] + c[2])*(c[1] + 7*c[2])*c[3] - 3*(517*c[1] 
	+ 328*c[2])*c[3]*c[3] + 531*c[3]*c[3]*c[3])) - 6*(13*(165*c[0]*c[0] + 616*
	c[1]*c[1] - 286*c[0]*(3*c[1] - 2*c[2]) + 220*c[1]*c[2] - 796*c[2]*c[2]) + 
	182*(77*c[0] - 120*c[1] + 18*c[2])*c[3] + 7528*c[3]*c[3])*c[4] + 12*(7319*
	c[0] - 7111*c[1] - 3478*c[2] + 5137*c[3])* c[4]*c[4] - 13192*c[4]*c[4]*c[4])
	 + a*b*(15015*c[0]*c[0]*c[0] - 7722*c[0]*c[0]*(7*c[2] - c[4]) + 156*c[0]*
	(231*c[1]*c[1] + 627*c[2]*c[2] - 594*c[1]*c[3] + 583*c[3]*c[3] - 682*c[2]*
	c[4] + 579*c[4]*c[4]) + 8*(15015 - 13*c[2]*(99*c[1]*c[1] + 517*c[2]*c[2] - 
	858*c[1]*c[3] + 963*c[3]*c[3]) + 3*(-1573*c[1]*c[1] + 4407*c[2]*c[2] + 
	3094*c[1]*c[3] - 1605*c[3]*c[3])*c[4] - 13395*c[2]*c[4]*c[4] + 
	967*c[4]*c[4]*c[4]))))/1081080.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
default:
	EH(-1,"unknown Chebyshev polynomial order");
	break;
}

/*  if opposite side is needed  */
if(switch_sense)
	{
	for(i=0;i<npts;i++)	{bf_mom[i]=gauss_mom[i]-bf_mom[i];}
	}
return;
}

/*  end of heaviside_chev_moments_2DQ		*/

#endif
/*   determine how interface crosses 3D linear elements	  */

int
interface_crossing_3DL( 
			const double ls_F[],
			double xf3D[12][2],
			int side_id[],
			double coord[12][MAX_PDIM]
			)
{
int i, j, iside, is3D=0;
int begnode[12]={0,1,3,0,4,5,7,4,0,1,2,3};
int endnode[12]={1,2,2,3,5,6,6,7,4,5,6,7};
double f0, f1, xf;
double dist,point, dist_tol=1.0E-06;
int dupl[12]={0,0,0,0,0,0,0,0,0,0,0,0};


for(iside=0 ; iside<12 ; iside++)	
	{
 	f0 = 0.5*(ls_F[endnode[iside]] + ls_F[begnode[iside]]);
 	f1 = 0.5*(ls_F[endnode[iside]] - ls_F[begnode[iside]]);
 	xf = -f0/f1;

	if(fabs(xf) <= 1.0)
 		{
 		xf3D[iside][0]=xf;
 		side_id[is3D]=iside;
 		is3D++;
 		}
	} /* end of for iside loop  */

/*	compute intersection point coordinates  */

for(i=0 ; i<is3D ; i++)
	{
	point = xf3D[side_id[i]][0];
 	switch (side_id[i])
 		{
 		case 0:
 			coord[i][0] = point;
 			coord[i][1] = -1;
 			coord[i][2] = -1;
 			break;
 		case 1:
 			coord[i][0] = 1;
 			coord[i][1] = point;
 			coord[i][2] = -1;
 			break;
 		case 2:
 			coord[i][0] = point;
 			coord[i][1] = 1;
 			coord[i][2] = -1;
 			break;
 		case 3:
 			coord[i][0] = -1;
 			coord[i][1] = point;
 			coord[i][2] = -1;
 			break;
 		case 4:
 			coord[i][0] = point;
 			coord[i][1] = -1;
 			coord[i][2] = 1;
 			break;
 		case 5:
 			coord[i][0] = 1;
 			coord[i][1] = point;
 			coord[i][2] = 1;
 			break;
 		case 6:
 			coord[i][0] = point;
 			coord[i][1] = 1;
 			coord[i][2] = 1;
 			break;
 		case 7:
 			coord[i][0] = -1;
 			coord[i][1] = point;
 			coord[i][2] = 1;
 			break;
 		case 8:
 			coord[i][0] = -1;
 			coord[i][1] = -1;
 			coord[i][2] = point;
 			break;
 		case 9:
 			coord[i][0] = 1;
 			coord[i][1] = -1;
 			coord[i][2] = point;
 			break;
 		case 10:
 			coord[i][0] = 1;
 			coord[i][1] = 1;
 			coord[i][2] = point;
 			break;
 		case 11:
 			coord[i][0] = -1;
 			coord[i][1] = 1;
 			coord[i][2] = point;
 			break;
 		} /* end of side_id[i] switch */
 	}

/* check for duplicate points  */

for(i=0 ; i<is3D ; i++)
 	{
 	for(j=i+1 ; j<is3D ; j++)
 		{
 		dist=pow(coord[i][0]-coord[j][0],2) +
 			pow(coord[i][1]-coord[j][1],2) +
 			  pow(coord[i][2]-coord[j][2],2);
 		if(dist < dist_tol)
 			{
/* 			printf("duplicate point %d\n",j);  */
 			dupl[j]=i;
 			}
 		}
 	}

 for(i=is3D-1 ; i >= 0 ; i--)
 	{
 	if(dupl[i])
 		{
 		for (j=i ; j < is3D ; j++)
 			{ side_id[j]=side_id[j+1]; 
 			  coord[j][0] = coord[j+1][0];
 			  coord[j][1] = coord[j+1][1];
 			  coord[j][2] = coord[j+1][2];
			}
 		is3D--;
 		}
 	}
 
return(is3D);
}
/*  end of interface_crossing_3DL*/


#ifndef NO_CHEBYSHEV_PLEASE
/*   computation of delta fcn integration moments using
	Chebyshev polynomial representation of interface		  */

void 
delta_chev_moments_2DQ( 
			const int interface_type,
			double bf_mom[],
			const double c[],
			const double d[],
			const int chev_order,
			const double end_pts[]
			)
{
int npts=9;
double a,b;

a = end_pts[0];
b = end_pts[1];
switch (chev_order)
{
case 3:
switch (interface_type)
	{
	case 0:
 		memset(bf_mom, 0, sizeof(double)*npts);
		break;
	case 2:
	case 3:
		bf_mom[0] = -(pow(-1 + a,2)*(12*(7*c[1]*(c[1] + a*(5 + 4*c[1])) - 
	7*(-3 + 2*a*(-1 + c[1]))*c[2] + (19 + 30*a)*c[2]*c[2])*d[0] - 24*(21*a*
	c[1]*c[1] + c[1]*(7 + 4*a*(7 - 4*c[2]) + 2*c[2]) + a*c[2]*(-7 + 11*c[2]))
	*d[1] + 8*(3*c[1]*(7*a + (-1 + 8*a)*c[1]) - 3*(19 + a*(30 + 22*c[1]))*c[2]
	 - (47 + 34*a)*c[2]*c[2])*d[2] + 6*c[0]*(-7*(5 + 6*c[2] + 2*a*(5 + 5*c[1] 
	+ 2*c[2]))*d[0] + 14*(2*c[1] + a*(5 + 8*c[1] - 2*c[2]))*d[1] - 2*(7*(-3 + 
	2*a*(-1 + c[1])) - 2*(19 + 30*a)*c[2])*d[2]) + 21*c[0]*c[0]*(5*(1 + 2*a)*
	d[0] - 2*(5*a*d[1] + (3 + 2*a)*d[2]))))/ 20160.;
		bf_mom[1] = -((-1 + a)*(21*c[0]*c[0]*(5*(5 + a*(5 + 2*a))*d[0] - 
	10*(-2 + a + a*a)*d[1] - 2*(7 + a*(11 + 2*a))*d[2]) + 6*c[0]*(-7*(5*(5 + 
	5*a + 2*a*a + 2*(-2 + a + a*a)*c[1]) + 2*(7 + a*(11 + 2*a))*c[2])*d[0] + 
	14*(2*(-5 + 9*c[1] + 2*c[2]) + a*(5 + 5*a + 14*c[1] + 8*a*c[1] - 2*(1 + a)
	*c[2]))*d[1] + 2*(7*(7 + 4*c[1] + a*(11 + 2*a - 2*(1 + a)*c[1])) + 2*(79 
	+ 87*a + 30*a*a)*c[2])*d[2]) + 4*(3*c[1]*c[1]*(7*(9 + a*(7 + 4*a))*d[0] - 
	42*(-2 + a + a*a)*d[1] + 2*(15 + a*(5 + 8*a))*d[2]) + 3*c[1]*(-7*(-2 + a 
	+ a*a)*(-5 + 2*c[2])*d[0] + 2*(-7*(9 + a*(7 + 4*a)) + 2*(15 + a*(5 + 8*a))
	*c[2])*d[1] - 2*(-2 + a + a*a)*(-7 + 22*c[2])*d[2]) + c[2]*(3*(49 + 79*c[2]
	 + a*(77 + 14*a + 87*c[2] + 30*a*c[2]))* d[0] - 6*(-2 + a + a*a)*(-7 + 
	11*c[2])*d[1] - 2*(237 + 115*c[2] + a*(261 + 90*a + 175*c[2] + 34*a*c[2]))
	*d[2]))))/20160.;
		bf_mom[2] = -((-1 + a)*(21*c[0]*c[0]*(5*(5 + a*(5 + 2*a))*d[0] - 
	10*(-2 + a + a*a)*d[1] - 2*(7 + a*(11 + 2*a))*d[2]) + 6*c[0]*(-7*(-5*(5 + 
	a*(5 + 2*a)) + 10*(-2 + a + a*a)*c[1] + 2*(7 + a*(11 + 2*a))*c[2])* d[0] + 
	14*(2*(5 + 9*c[1] + 2*c[2]) + a*(-5 - 5*a + 14*c[1] + 8*a*c[1] - 2*(1 + a)
	*c[2]))*d[1] - 2*(7*(7 + 11*a + 2*a*a + 2*(-2 + a + a*a)*c[1]) - 2*(79 + 
	87*a + 30*a*a)*c[2])*d[2]) + 4*(3*c[1]*c[1]*(7*(9 + a*(7 + 4*a))*d[0] - 
	42*(-2 + a + a*a)*d[1] + 2*(15 + a*(5 + 8*a))*d[2]) + 3*c[1]*(-7*(-2 + a + 
	a*a)*(5 + 2*c[2])*d[0] + 2*(7*(9 + a*(7 + 4*a)) + 2*(15 + a*(5 + 8*a))
	*c[2])*d[1] - 2*(-2 + a + a*a)*(7 + 22*c[2])*d[2]) + c[2]*(3*(-7*(7 + a*
	(11 + 2*a)) + (79 + 87*a + 30*a*a)*c[2])*d[0] - 6*(-2 + a + a*a)*(7 + 11*
	c[2])*d[1] - 2*(-3*(79 + 87*a + 30*a*a) + (115 + a*(175 + 34*a))*c[2])
	*d[2]))))/20160.;
		bf_mom[3] = -(pow(-1 + a,2)*(12*(7*c[1]*(c[1] + a*(-5 + 4*c[1])) - 
	7*(3 + 2*a*(1 + c[1]))*c[2] + (19 + 30*a)*c[2]*c[2])*d[0] - 24*(21*a*c[1]*
	c[1] + a*c[2]*(7 + 11*c[2]) - c[1]*(7 - 2*c[2] + 4*a*(7 + 4*c[2])))*d[1] + 
	8*(3*c[1]*(-7*a + (-1 + 8*a)*c[1]) - 3*(-19 + a*(-30 + 22*c[1]))*c[2] - 
	(47 + 34*a)*c[2]*c[2])* d[2] + 6*c[0]*(-7*(-5 + 6*c[2] + 2*a*(-5 + 5*c[1] 
	+ 2*c[2]))* d[0] + 14*(2*c[1] + a*(-5 + 8*c[1] - 2*c[2]))*d[1] - 2*(7*(3 + 
	2*a*(1 + c[1])) - 2*(19 + 30*a)*c[2])*d[2]) + 21*c[0]*c[0]*(5*(1 + 2*a)*
	d[0] - 2*(5*a*d[1] + (3 + 2*a)*d[2]))))/ 20160.;
		bf_mom[4] = (pow(-1 + a,2)*(6*(14*(3 + 2*a)*c[1]*c[1] + 7*(1 + a)*
	c[1]*(5 - 2*c[2]) + 2*c[2]*(7*(4 + a) + (34 + 15*a)*c[2]))*d[0] - 12*(21*
	(1 + a)*c[1]*c[1] + (1 + a)*c[2]*(-7 + 11*c[2]) - 2*c[1]*(-21 - 14*a + 
	6*c[2] + 8*a*c[2]))*d[1] + 4*(6*(3 + 4*a)*c[1]*c[1] + 3*(1 + a)*c[1]*(7 - 
	22*c[2]) - 2*c[2]*(102 + 45*a + 64*c[2] + 17*a*c[2]))*d[2] + 21*c[0]*c[0]*
	(5*(2 + a)*d[0] - 5*(1 + a)*d[1] - 2*(4 + a)*d[2]) + 6*c[0]*(-7*(5*(2 + a 
	+ c[1] + a*c[1]) + 2*(4 + a)*c[2])*d[0] + 7*(5 + 5*a + 12*c[1] + 8*a*c[1] 
	- 2*(1 + a)*c[2])*d[1] - 2*(-7*(4 + a) + 7*(1 + a)*c[1] - 2*(34 + 15*a)
	*c[2])*d[2])))/5040.;
		bf_mom[5] = ((-1 + a)*(3*(7*(5*(5 + a*(5 + 2*a))*(-4 + c[0]*c[0]) - 
	20*(-2 + a + a*a)*c[0]*c[1] + 4*(9 + a*(7 + 4*a))*c[1]*c[1]) - 28*((7 + a*
	(11 + 2*a))*c[0] + 2*(-2 + a + a*a)*c[1])* c[2] + 4*(79 + 87*a + 30*a*a)
	*c[2]*c[2])*d[0] - 6*(35*(-2 + a + a*a)*(-4 + c[0]*c[0]) - 28*(9 + a*(7 + 
	4*a))*c[0]*c[1] + 84*(-2 + a + a*a)*c[1]*c[1] + 4*(7*(-2 + a + a*a)*c[0] 
	- 2*(15 + a*(5 + 8*a))*c[1])* c[2] + 44*(-2 + a + a*a)*c[2]*c[2])*d[1] - 
	2*(21*(7 + a*(11 + 2*a))*(-4 + c[0]*c[0]) + 84*(-2 + a + a*a)*c[0]*c[1] - 
	12*(15 + a*(5 + 8*a))*c[1]*c[1] - 12*((79 + 87*a + 30*a*a)*c[0] - 22*(-2 + 
	a + a*a)*c[1])*c[2] + 4*(115 + a*(175 + 34*a))*c[2]*c[2])*d[2]))/10080.;
		bf_mom[6] = (pow(-1 + a,2)*(6*(14*(3 + 2*a)*c[1]*c[1] - 7*(1 + a)*
	c[1]*(5 + 2*c[2]) + 2*c[2]*(-7*(4 + a) + (34 + 15*a)*c[2]))*d[0] - 12*(21*
	(1 + a)*c[1]*c[1] + (1 + a)*c[2]*(7 + 11*c[2]) - 2*c[1]*(21 + 14*a + 6*c[2]
	 + 8*a*c[2]))*d[1] + 4*(6*(3 + 4*a)*c[1]*c[1] - 3*(1 + a)*c[1]*(7 + 22*c[2])
	 + 2*c[2]*(102 + 45*a - 64*c[2] - 17*a*c[2]))*d[2] + 21*c[0]*c[0]*(5*(2+a)
	*d[0] - 5*(1 + a)*d[1] - 2*(4 + a)*d[2]) + 6*c[0]*(-7*(-5*(2 + a) + 5*(1+a)
	*c[1] + 2*(4 + a)*c[2])*d[0] + 7*(-5 - 5*a + 12*c[1] + 8*a*c[1] - 2*(1+a)
	*c[2])*d[1] - 2*(7*(4 + a + c[1] + a*c[1]) - 2*(34 + 15*a)*c[2])
	*d[2])))/5040.;
		bf_mom[7] = (pow(-1 + a,2)*(3*(7*(5*(1 + 2*a)*(-4 + c[0]*c[0]) - 
	20*a*c[0]*c[1] + 4*(1 + 4*a)*c[1]*c[1]) - 28*((3 + 2*a)*c[0] + 2*a*c[1])
	*c[2] + 4*(19 + 30*a)*c[2]*c[2])*d[0] - 6*(35*a*(-4 + c[0]*c[0]) - 28*(1 + 
	4*a)*c[0]*c[1] + 84*a*c[1]*c[1] + 4*(7*a*c[0] + 2*c[1] - 16*a*c[1])*c[2] + 
	44*a*c[2]*c[2])*d[1] - 2*(21*(3 + 2*a)*(-4 + c[0]*c[0]) + 84*a*c[0]*c[1] - 
	12*(-1 + 8*a)*c[1]*c[1] - 12*((19 + 30*a)*c[0] - 22*a*c[1])*c[2] + 4*(47 + 
	34*a)*c[2]*c[2])*d[2]))/10080.;
		bf_mom[8] = -(pow(-1 + a,2)*(3*(35*(2 + a)*(-4 + c[0]*c[0]) - 70*
	(1 + a)*c[0]*c[1] + 28*(3 + 2*a)*c[1]*c[1] - 28*((4 + a)*c[0] + (1 + a)
	*c[1])*c[2] + 4*(34 + 15*a)*c[2]*c[2])*d[0] - 3*(35*(1 + a)*(-4 + c[0]*c[0])
	 - 56*(3 + 2*a)*c[0]*c[1] + 84*(1 + a)*c[1]*c[1] + 4*(7*(1 + a)*c[0] - 4*
	(3 + 4*a)*c[1])*c[2] + 44*(1 + a)*c[2]*c[2])*d[1] - 2*(21*(4 + a)*(-4 + 
	c[0]*c[0]) + 42*(1 + a)*c[0]*c[1] - 12*(3 + 4*a)*c[1]*c[1] - 12*((34 + 15*a)
	*c[0] - 11*(1 + a)*c[1])*c[2] + 4*(64 + 17*a)*c[2]*c[2])*d[2]))/2520.;
		break;
	case 4:
	case 6:
		bf_mom[0] = ((1 + b)*(21*c[0]*c[0]*(5*(5 + b*(-5 + 2*b))*d[0] + 
	10*(-2 + b)*(1 + b)*d[1] - 2*(7 + b*(-11 + 2*b))*d[2]) + 6*c[0]*(7*(5*(-5 
	+ 5*b - 2*b*b + 2*(-2 + b)*(1 + b)*c[1]) - 2*(7 + b*(-11 + 2*b))*c[2])*d[0]
	 + 14*(2*(5 + 9*c[1] - 2*c[2]) + b*(5 - 5*b - 14*c[1] + 8*b*c[1] + 2*
	(-1 + b)*c[2]))*d[1] + 2*(7*(7 - 11*b + 2*b*b + 2*(-2 + b)*(1 + b)*c[1]) + 
	2*(79 - 87*b + 30*b*b)*c[2])*d[2]) + 4*(3*c[1]*c[1]*(7*(9 + b*(-7 + 4*b))
	*d[0] + 42*(-2 + b)*(1 + b)*d[1] + 2*(15 + b*(-5 + 8*b))*d[2]) + 3*c[1]*
	(7*(-2 + b)*(1 + b)*(-5 + 2*c[2])*d[0] + 2*(-7*(9 + b*(-7 + 4*b)) + 2*(15 
	+ b*(-5 + 8*b))*c[2])*d[1] + 2*(-2 + b)*(1 + b)*(-7 + 22*c[2])*d[2]) + 
	c[2]*(3*(49 + 79*c[2] + b*(-77 + 14*b - 87*c[2] + 30*b*c[2]))* d[0] + 6*
	(-2 + b)*(1 + b)*(-7 + 11*c[2])*d[1] - 2*(237 + 115*c[2] + b*(-261 + 90*b 
	- 175*c[2] + 34*b*c[2]))*d[2]))))/20160.;
		bf_mom[1] = (pow(1 + b,2)*(12*(7*c[1]*(-5*b + (-1 + 4*b)*c[1]) + 
	7*(-3 + 2*b*(1 + c[1]))*c[2] + (-19 + 30*b)*c[2]*c[2])*d[0] + 24*(21*b*c[1]
	*c[1] + b*c[2]*(-7 + 11*c[2]) + c[1]*(7 + 2*c[2] + 4*b*(-7 + 4*c[2])))*d[1]
	 + 8*(3*c[1]*(c[1] + b*(-7 + 8*c[1])) + 3*(19 + b*(-30 + 22*c[1]))*c[2] - 
	(-47 + 34*b)*c[2]*c[2])*d[2] + 21*c[0]*c[0]*(5*(-1 + 2*b)*d[0] + 10*b*d[1] 
	+ 6*d[2] - 4*b*d[2]) + 6*c[0]*(7*(5 + 2*b*(-5 + 5*c[1] - 2*c[2]) + 6*c[2])
	*d[0] + 14*(-2*c[1] + b*(-5 + 8*c[1] + 2*c[2]))*d[1] + 2*(-21 - 38*c[2] + 
	2*b*(7 + 7*c[1] + 30*c[2]))*d[2])))/20160.;
		bf_mom[2] = (pow(1 + b,2)*(12*(7*c[1]*(5*b + (-1 + 4*b)*c[1]) + 7*
	(3 + 2*b*(-1 + c[1]))*c[2] + (-19 + 30*b)*c[2]*c[2])*d[0] + 24*(21*b*c[1]*
	c[1] + b*c[2]*(7 + 11*c[2]) + c[1]*(-7 + 2*c[2] + 4*b*(7 + 4*c[2])))*d[1] + 
	8*(3*c[1]*(c[1] + b*(7 + 8*c[1])) + 3*(-19 + b*(30 + 22*c[1]))*c[2] - (-47 
	+ 34*b)*c[2]*c[2])*d[2] + 21*c[0]*c[0]*(5*(-1 + 2*b)*d[0] + 10*b*d[1] + 
	6*d[2] - 4*b*d[2]) + 6*c[0]*(7*(-5 + 2*b*(5 + 5*c[1] - 2*c[2]) + 6*c[2])
	*d[0] + 14*(-2*c[1] + b*(5 + 8*c[1] + 2*c[2]))*d[1] + 2*(21 - 38*c[2] + 
	2*b*(-7 + 7*c[1] + 30*c[2]))*d[2])))/20160.;
		bf_mom[3] = ((1 + b)*(21*c[0]*c[0]*(5*(5 + b*(-5 + 2*b))*d[0] + 10*
	(-2 + b)*(1 + b)*d[1] - 2*(7 + b*(-11 + 2*b))*d[2]) + 6*c[0]*(7*(5*(5 - 5*b
	 + 2*b*b + 2*(-2 + b)*(1 + b)*c[1]) - 2*(7 + b*(-11 + 2*b))*c[2])*d[0] + 
	14*(2*(-5 + 9*c[1] - 2*c[2]) + b*(-5 + 5*b - 14*c[1] + 8*b*c[1] + 2*(-1+b)
	*c[2]))*d[1] + 2*(7*(-7 + 11*b - 2*b*b + 2*(-2 + b)*(1 + b)*c[1]) + 2*(79 
	- 87*b + 30*b*b)*c[2])*d[2]) + 4*(3*c[1]*c[1]*(7*(9 + b*(-7 + 4*b))*d[0] + 
	42*(-2 + b)*(1 + b)*d[1] + 2*(15 + b*(-5 + 8*b))*d[2]) + 3*c[1]*(7*(-2 + b)
	*(1 + b)*(5 + 2*c[2])*d[0] + 2*(7*(9 + b*(-7 + 4*b)) + 2*(15 + b*(-5+8*b))
	*c[2])*d[1] + 2*(-2 + b)*(1 + b)*(7 + 22*c[2])*d[2]) + c[2]*(3*(-49 + 
	79*c[2] + b*(77 - 87*c[2] + 2*b*(-7 + 15*c[2])))* d[0] + 6*(-2 + b)*(1+b)
	*(7 + 11*c[2])*d[1] + 2*(237 - 115*c[2] + b*(-261 + 90*b + 175*c[2] - 
	34*b*c[2]))*d[2]))))/20160.;
		bf_mom[4] = -(pow(1 + b,2)*(84*c[1]*c[1]*((-3 + 2*b)*d[0] + 3*
	(-1 + b)*d[1]) + 12*c[2]*((7*(-4 + b) + (-34 + 15*b)*c[2])*d[0] + (-1 + b)
	*(-7 + 11*c[2])*d[1]) + 6*c[1]*(7*(-1 + b)*(-5 + 2*c[2])*d[0] + 4*(21 - 14*b
	 - 6*c[2] + 8*b*c[2])*d[1]) + 4*(6*(-3 + 4*b)*c[1]*c[1] + 3*(-1 + b)*c[1]*
	(-7 + 22*c[2]) + 2*c[2]*(102 - 45*b + 64*c[2] - 17*b*c[2]))*d[2] + 21*c[0]
	*c[0]*(5*(-2 + b)*d[0] + 5*(-1 + b)*d[1] - 2*(-4 + b)*d[2]) + 6*c[0]* (7*
	(10 - 5*c[1] + b*(-5 + 5*c[1] - 2*c[2]) + 8*c[2])*d[0] + 7*(5 - 5*b - 12*
	c[1] + 8*b*c[1] + 2*(-1 + b)*c[2])*d[1] + 2*(7*(-4 + b + (-1 + b)*c[1]) + 
	2*(-34 + 15*b)*c[2])*d[2])))/5040.;
		bf_mom[5] = -(pow(1 + b,2)*(3*(7*(5*(-1 + 2*b)*(-4 + c[0]*c[0]) + 
	20*b*c[0]*c[1] + 4*(-1 + 4*b)*c[1]*c[1]) - 28*((-3 + 2*b)*c[0] - 2*b*c[1])
	*c[2] + 4*(-19 + 30*b)*c[2]*c[2])*d[0] + 6*(-28*c[0]*c[1] + 7*b*(-20 + 
	(c[0] + 2*c[1])*(5*c[0] + 6*c[1])) + 4*(7*b*c[0] + 2*(1 + 8*b)*c[1])*c[2] 
	+ 44*b*c[2]*c[2])*d[1] - 2*(21*(-3 + 2*b)*(-4 + c[0]*c[0]) - 84*b*c[0]*c[1]
	 - 12*(1 + 8*b)*c[1]*c[1] - 12*((-19 + 30*b)*c[0] + 22*b*c[1])*c[2] + 4*
	(-47 + 34*b)*c[2]*c[2])*d[2]))/10080.;
		bf_mom[6] = -(pow(1 + b,2)*(84*c[1]*c[1]*((-3 + 2*b)*d[0] + 3*
	(-1 + b)*d[1]) + 12*c[2]*((28 - 7*b - 34*c[2] + 15*b*c[2])*d[0] + (-1 + b)
	*(7 + 11*c[2])*d[1]) + 6*c[1]*(7*(-1 + b)*(5 + 2*c[2])*d[0] + 4*(-21 + 
	14*b - 6*c[2] + 8*b*c[2])*d[1]) + 4*(6*(-3 + 4*b)*c[1]*c[1] + 3*(-1 + b)
	*c[1]*(7 + 22*c[2]) + 2*c[2]*(-102 + 45*b + 64*c[2] - 17*b*c[2]))*d[2] + 
	21*c[0]*c[0]*(5*(-2 + b)*d[0] + 5*(-1 + b)*d[1] - 2*(-4 + b)*d[2]) + 6*
	c[0]* (7*(5*(-2 + b + (-1 + b)*c[1]) - 2*(-4 + b)*c[2])*d[0] + 7*(-5 + 5*b
	 - 12*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2])*d[1] + 2*(28 - 7*c[1] - 68*c[2] 
	+ b*(-7 + 7*c[1] + 30*c[2]))*d[2])))/5040.;
		bf_mom[7] = -((1 + b)*(3*(7*(5*(5 + b*(-5 + 2*b))*(-4 + c[0]*c[0])
	 + 20*(-2 + b)*(1 + b)*c[0]*c[1] + 4*(9 + b*(-7 + 4*b))*c[1]*c[1]) - 28*(
	(7 + b*(-11 + 2*b))*c[0] - 2*(-2 + b)*(1 + b)*c[1])*c[2] + 4*(79 - 87*b + 
	30*b*b)*c[2]*c[2])*d[0] + 6*(35*(-2 + b)*(1 + b)*(-4 + c[0]*c[0]) + 28*
	(9 + b*(-7 + 4*b))*c[0]*c[1] + 84*(-2 + b)*(1 + b)*c[1]*c[1] + 4*(7*
	(-2 + b)*(1 + b)*c[0] + 2*(15 + b*(-5 + 8*b))*c[1])*c[2] + 44*(-2 + b)*
	(1 + b)*c[2]*c[2])*d[1] - 2*(21*(7 + b*(-11 + 2*b))*(-4 + c[0]*c[0]) - 84*
	(-2 + b)*(1 + b)*c[0]*c[1] - 12*(15 + b*(-5 + 8*b))*c[1]*c[1] - 12*((79 - 
	87*b + 30*b*b)*c[0] + 22*(-2 + b)*(1 + b)*c[1])*c[2] + 4*(115 + b*(-175 + 
	34*b))*c[2]*c[2])*d[2]))/10080.;
		bf_mom[8] = (pow(1 + b,2)*(3*(35*(-2 + b)*(-4 + c[0]*c[0]) + 70*
	(-1 + b)*c[0]*c[1] + 28*(-3 + 2*b)*c[1]*c[1] - 28*((-4 + b)*c[0] + c[1] - 
	b*c[1])*c[2] + 4*(-34 + 15*b)*c[2]*c[2])*d[0] + 3*(35*(-1 + b)*(-4 + c[0]*
	c[0]) + 56*(-3 + 2*b)*c[0]*c[1] + 84*(-1 + b)*c[1]*c[1] + 4*(7*(-1 + b)*
	c[0] + 4*(-3 + 4*b)*c[1])*c[2] + 44*(-1 + b)*c[2]*c[2])*d[1] + 2*(-21*(-4 
	+ b)*(-4 + c[0]*c[0]) + 42*(-1 + b)*c[0]*c[1] + 12*(-3 + 4*b)*c[1]*c[1] + 
	12*((-34 + 15*b)*c[0] + 11*(-1 + b)*c[1])*c[2] - 4*(-64 + 17*b)*c[2]
	*c[2])*d[2]))/2520.;
		break;
	case 22:
		bf_mom[0] = (3*(35*(-2 + c[0])*c[0] - 140*(-1 + c[0])*c[1] + 
	84*c[1]*c[1] + 28*(-1 + c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] - 6*
	(35*c[0]*c[0] + 84*c[1]*(1 + c[1]) - 14*c[0]*(5 + 6*c[1] - 2*c[2]) - 4*
	(7 + 18*c[1])*c[2] + 44*c[2]*c[2])*d[1] + 2*(21*(-2 + c[0])*c[0] - 84*(-1 
	+ c[0])*c[1] + 108*c[1]*c[1] + 132*(-1 + c[0] - 2*c[1])*c[2]+52*c[2]*c[2])
	*d[2])/5040.;
		bf_mom[1] = (3*(35*c[0]*(2 + c[0]) - 140*(1 + c[0])*c[1] + 84*c[1]
	*c[1] + 28*(1 + c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] - 6*(35*c[0]*(2 + 
	c[0]) - 84*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 18*c[1])*c[2] 
	+ 44*c[2]*c[2])*d[1] + 2*(21*c[0]*(2 + c[0]) - 84*(1 + c[0])*c[1] + 
	108*c[1]*c[1] + 132*(1 + c[0] - 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/5040.;
		bf_mom[2] = (3*(35*c[0]*(2 + c[0]) + 140*(1 + c[0])*c[1] + 84*c[1]
	*c[1] + 28*(1 + c[0] + 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(35*c[0]*(2 
	+ c[0]) + 84*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] + 18*c[1])*c[2]
	 + 44*c[2]*c[2])*d[1] + 2*(21*c[0]*(2 + c[0]) + 84*(1 + c[0])*c[1] + 
	108*c[1]*c[1] + 132*(1 + c[0] + 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/5040.;
		bf_mom[3] = (3*(35*(-2 + c[0])*c[0] + 140*(-1 + c[0])*c[1] + 84*
	c[1]*c[1] + 28*(-1 + c[0] + 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(35*(-2 
	+ c[0])*c[0] + 84*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0]+18*c[1])
	*c[2] + 44*c[2]*c[2])*d[1] + 2*(21*(-2 + c[0])*c[0] + 84*(-1 + c[0])*c[1] 
	+ 108*c[1]*c[1] + 132*(-1 + c[0] + 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/5040.;
		bf_mom[4] = (-3*(7*(-20 + 5*c[0]*c[0] - 20*c[0]*c[1] + 12*c[1]*c[1])
	 + 28*(c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(-140 + 35*c[0]*c[0] + 
	84*c[1]*c[1] - 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(-3*c[1] + c[2]))*d[1] 
	- 2*(3*(-28 + 7*c[0]*c[0] - 28*c[0]*c[1] + 36*c[1]*c[1]) + 132*(c[0] - 
	2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/2520.;
		bf_mom[5] = (3*(35*c[0]*(2 + c[0]) + 28*c[1]*c[1] - 84*(1 + c[0])
	*c[2] + 76*c[2]*c[2])*d[0] + 24*c[1]*(7 + 7*c[0] - 2*c[2])*d[1] - 2*(63*
	c[0]*(2 + c[0]) + 12*c[1]*c[1] - 228*(1 + c[0])*c[2] + 188*c[2]*c[2])*
	d[2])/1260.;
		bf_mom[6] = (-3*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] + 56*c[1]*c[2] 
	+ 44*c[2]*c[2] + 28*c[0]*(5*c[1] + c[2]))*d[0] - 6*(-140 + 35*c[0]*c[0] + 
	84*c[1]*c[1] + 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(3*c[1] + c[2]))*d[1] 
	- 2*(3*(-28 + 7*c[0]*c[0] + 28*c[0]*c[1] + 36*c[1]*c[1]) + 132*(c[0] + 
	2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/2520.;
		bf_mom[7] = (3*(35*(-2 + c[0])*c[0] + 28*c[1]*c[1] - 84*(-1+c[0])
	*c[2] + 76*c[2]*c[2])*d[0] + 24*c[1]*(-7 + 7*c[0] - 2*c[2])*d[1] - 2*(63*
	(-2 + c[0])*c[0] + 12*c[1]*c[1] - 228*(-1 + c[0])*c[2] + 188*c[2]*c[2])
	*d[2])/1260.;
		bf_mom[8] = (-3*(35*c[0]*c[0] + 28*(-5 + c[1]*c[1]) - 84*c[0]*c[2] 
	+ 76*c[2]*c[2])*d[0] + 24*c[1]*(-7*c[0] + 2*c[2])*d[1] + 2*(63*c[0]*c[0] 
	- 228*c[0]*c[2] + 4*(-63 + 3*c[1]*c[1] + 47*c[2]*c[2]))*d[2])/630.;
		break;
	case 32:
		bf_mom[0] = (3*(35*(-2 + c[0])*c[0] - 140*(-1 + c[0])*c[1] + 
	84*c[1]*c[1] + 28*(-1 + c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] - 6*
	(35*c[0]*c[0] + 84*c[1]*(1 + c[1]) - 14*c[0]*(5 + 6*c[1] - 2*c[2]) - 4*
	(7 + 18*c[1])*c[2] + 44*c[2]*c[2])*d[1] + 2*(21*(-2 + c[0])*c[0] - 84*
	(-1 + c[0])*c[1] + 108*c[1]*c[1] + 132*(-1 + c[0] - 2*c[1])*c[2] + 
	52*c[2]*c[2])*d[2])/5040.;
		bf_mom[1] = (3*(35*(-2 + c[0])*c[0] + 140*(-1 + c[0])*c[1] + 84*
	c[1]*c[1] + 28*(-1 + c[0] + 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(35*(-2 
	+ c[0])*c[0] + 84*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0] + 
	18*c[1])*c[2] + 44*c[2]*c[2])*d[1] + 2*(21*(-2 + c[0])*c[0] + 84*(-1+c[0])
	*c[1] + 108*c[1]*c[1] + 132*(-1 + c[0] + 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])
	/5040.;
		bf_mom[2] = (3*(35*c[0]*(2 + c[0]) + 140*(1 + c[0])*c[1] + 84*c[1]
	*c[1] + 28*(1 + c[0] + 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(35*c[0]*(2 + 
	c[0]) + 84*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] + 18*c[1])*c[2] 
	+ 44*c[2]*c[2])*d[1] + 2*(21*c[0]*(2 + c[0]) + 84*(1 + c[0])*c[1] + 108*
	c[1]*c[1] + 132*(1 + c[0] + 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/5040.;
		bf_mom[3] = (3*(35*c[0]*(2 + c[0]) - 140*(1 + c[0])*c[1] + 84*c[1]
	*c[1] + 28*(1 + c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] - 6*(35*c[0]*(2 + 
	c[0]) - 84*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 18*c[1])*c[2] 
	+ 44*c[2]*c[2])*d[1] + 2*(21*c[0]*(2 + c[0]) - 84*(1 + c[0])*c[1] + 108*
	c[1]*c[1] + 132*(1 + c[0] - 2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/5040.;
		bf_mom[4] = (3*(35*(-2 + c[0])*c[0] + 28*c[1]*c[1] - 84*(-1+c[0])
	*c[2] + 76*c[2]*c[2])*d[0] + 24*c[1]*(-7 + 7*c[0] - 2*c[2])*d[1] - 2*(63*
	(-2 + c[0])*c[0] + 12*c[1]*c[1] - 228*(-1 + c[0])*c[2]+188*c[2]*c[2])
	*d[2])/1260.;
		bf_mom[5] = (-3*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] + 56*c[1]*c[2] 
	+ 44*c[2]*c[2] + 28*c[0]*(5*c[1] + c[2]))*d[0] - 6*(-140 + 35*c[0]*c[0] + 
	84*c[1]*c[1] + 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(3*c[1] + c[2]))*d[1] 
	- 2*(3*(-28 + 7*c[0]*c[0] + 28*c[0]*c[1] + 36*c[1]*c[1]) + 132*(c[0] + 
	2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/2520.;
		bf_mom[6] = (3*(35*c[0]*(2 + c[0]) + 28*c[1]*c[1] - 84*(1 + c[0])
	*c[2] + 76*c[2]*c[2])*d[0] + 24*c[1]*(7 + 7*c[0] - 2*c[2])*d[1] - 2*(63*
	c[0]*(2 + c[0]) + 12*c[1]*c[1] - 228*(1 + c[0])*c[2] + 188*c[2]*c[2])
	*d[2])/1260.;
		bf_mom[7] = (-3*(7*(-20 + 5*c[0]*c[0] - 20*c[0]*c[1] + 12*c[1]*c[1])
	 + 28*(c[0] - 2*c[1])*c[2] + 44*c[2]*c[2])*d[0] + 6*(-140 + 35*c[0]*c[0] + 
	84*c[1]*c[1] - 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(-3*c[1] + c[2]))*d[1] 
	- 2*(3*(-28 + 7*c[0]*c[0] - 28*c[0]*c[1] + 36*c[1]*c[1]) + 132*(c[0] - 
	2*c[1])*c[2] + 52*c[2]*c[2])*d[2])/2520.;
		bf_mom[8] = (-3*(35*c[0]*c[0] + 28*(-5 + c[1]*c[1]) - 84*c[0]*c[2] 
	+ 76*c[2]*c[2])*d[0] + 24*c[1]*(-7*c[0] + 2*c[2])*d[1] + 2*(63*c[0]*c[0] - 
	228*c[0]*c[2] + 4*(-63 + 3*c[1]*c[1] + 47*c[2]*c[2]))*d[2])/630.;
		break;
	case 100:
	case 102:
		bf_mom[0] = -((a - b)*(2*a*a*(3*(35*c[0]*c[0] - 70*c[0]*(1 + c[1]) 
	+ 14*c[1]*(5 + 4*c[1]) - 28*(-1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 
	3*(35*(-2 + c[0])*c[0] - 112*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*
	c[0] - 16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*c[0] - 6*c[1]*(7 + 
	8*c[1]) + 6*c[0]*(-7 + 7*c[1] - 30*c[2]) + 12*(15 + 11*c[1])*c[2] + 68*c[2]
	*c[2])*d[2]) + a*(3*(7*(5*(-3 + 2*b)*(-2 + c[0])*c[0] + 20*(-1 + c[0])*c[1]
	 + 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5 + 6*b)*(-1 + c[0]) - 2*c[1])*c[2] + 4*
	(-49 + 38*b)*c[2]*c[2])*d[0] + 6*(35*(-2 + c[0])*c[0] + 28*(-5 + 2*b)*(-1 + 
	c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0] - 2*(7 + 2*b)*c[1])*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(21*(-5 + 6*b)*(-2 + c[0])*c[0] - 84*(-1 + c[0])*c[1] 
	+ 12*(7 + 2*b)*c[1]*c[1] - 12*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2] + 
	4*(-81 + 94*b)*c[2]*c[2])*d[2]) + b*(12*(7*(-5 + 4*b)*c[1]*c[1] + 7*(-1 + b)
	*c[1]*(-5 + 2*c[2]) + c[2]*(-7*(5 + 7*c[2]) + 2*b*(7 + 15*c[2])))*d[0] + 
	24*(21*(-1 + b)*c[1]*c[1] + (-1 + b)*c[2]*(-7 + 11*c[2]) + c[1]*(35 - 14*
	c[2] + 4*b*(-7 + 4*c[2])))*d[1] + 8*(3*(-7 + 8*b)*c[1]*c[1] + 3*(-1 + b)*
	c[1]*(-7 + 22*c[2]) + c[2]*(147 + 81*c[2] - 2*b*(45 + 17*c[2])))*d[2] + 21*
	c[0]*c[0]*(5*(-3 + 2*b)*d[0] + 10*(-1 + b)*d[1] + 2*(5 - 2*b)*d[2]) + 6*
	c[0]*(7*(5*(3 - 2*b + 2*(-1 + b)*c[1]) - 2*(-5 + 2*b)*c[2])* d[0] + 14*(5 
	- 5*b - 10*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2])* d[1] + 2*(7*(-5 + 2*b + 2*
	(-1 + b)*c[1]) + 2*(-49 + 30*b)*c[2])*d[2]))))/20160.;
		bf_mom[1] = -((a - b)*(2*a*a*(3*(35*c[0]*c[0] - 70*c[0]*(1 + c[1]) 
	+ 14*c[1]*(5 + 4*c[1]) - 28*(-1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 
	3*(35*(-2 + c[0])*c[0] - 112*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*
	c[0] - 16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*c[0] - 6*c[1]*(7 + 
	8*c[1]) + 6*c[0]*(-7 + 7*c[1] - 30*c[2]) + 12*(15 + 11*c[1])*c[2] + 68*c[2]
	*c[2])*d[2]) + a*(3*(7*(5*(3 + 2*b)*(-2 + c[0])*c[0] - 20*(-1 + c[0])*c[1] 
	+ 4*(5 + 2*b)*c[1]*c[1]) - 28*((5 + 6*b)*(-1 + c[0]) + 2*c[1])*c[2] + 4*
	(49 + 38*b)*c[2]*c[2])*d[0] - 6*(35*(-2 + c[0])*c[0] - 28*(5 + 2*b)*(-1 + 
	c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(21*(5 + 6*b)*(-2 + c[0])*c[0] + 84*(-1 + c[0])*c[1] 
	+ 12*(-7 + 2*b)*c[1]*c[1] - 12*((49 + 38*b)*(-1 + c[0]) - 22*c[1])*c[2] + 
	4*(81 + 94*b)*c[2]*c[2])*d[2]) + b*(12*(7*(5 + 4*b)*c[1]*c[1] + 7*(1 + b)*
	c[1]*(-5 + 2*c[2]) + c[2]*(35 + 14*b + 49*c[2] + 30*b*c[2]))*d[0] + 24*(21*
	(1 + b)*c[1]*c[1] + (1 + b)*c[2]*(-7 + 11*c[2]) + c[1]*(-7*(5 + 4*b) + 2*
	(7 + 8*b)*c[2]))*d[1] + 8*(3*(7 + 8*b)*c[1]*c[1] + 3*(1 + b)*c[1]*(-7 + 
	22*c[2]) - c[2]*(147 + 90*b + 81*c[2] + 34*b*c[2]))*d[2] + 21*c[0]*c[0]*(5*
	(3 + 2*b)*d[0] + 10*(1 + b)*d[1] - 2*(5 + 2*b)*d[2]) + 6*c[0]*(7*(5*(-3 - 
	2*b + 2*(1 + b)*c[1]) - 2*(5 + 2*b)*c[2])* d[0] + 14*(-5 - 5*b + 10*c[1] + 
	8*b*c[1] + 2*(1 + b)*c[2])* d[1] + 2*(7*(5 + 2*b + 2*(1 + b)*c[1]) + 2*
	(49 + 30*b)*c[2])* d[2]))))/20160.;
		bf_mom[2] = -((a - b)*(2*a*a*(3*(35*c[0]*(2 + c[0]) - 70*(1 + c[0])
	*c[1] + 56*c[1]*c[1] - 28*(1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*
	(35*c[0]*(2 + c[0]) - 112*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 
	16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*(2 + c[0]) + 42*(1 + c[0])
	*c[1] - 48*c[1]*c[1] - 12*(15 + 15*c[0] - 11*c[1])*c[2] + 68*c[2]*c[2])
	*d[2]) + a*(3*(7*(5*(3 + 2*b)*c[0]*(2 + c[0]) - 20*(1 + c[0])*c[1] + 4*(5 
	+ 2*b)*c[1]*c[1]) - 28*((5 + 6*b)*(1 + c[0]) + 2*c[1])*c[2] + 4*(49 + 38*b)
	*c[2]*c[2])*d[0] - 6*(35*c[0]*(2 + c[0]) - 28*(5 + 2*b)*(1 + c[0])*c[1] + 
	84*c[1]*c[1] + 4*(7 + 7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])*
	 d[1] - 2*(21*(5 + 6*b)*c[0]*(2 + c[0]) + 84*(1 + c[0])*c[1] + 12*(-7+2*b)
	*c[1]*c[1] - 12*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2] + 4*(81 + 94*b)*
	c[2]*c[2])*d[2]) + b*(12*(7*(5 + 4*b)*c[1]*c[1] + 7*(1 + b)*c[1]*(5+2*c[2])
	 + c[2]*(-7*(5 + 2*b) + (49 + 30*b)*c[2]))*d[0] + 24*(21*(1 + b)*c[1]*c[1] 
	+ (1 + b)*c[2]*(7 + 11*c[2]) + c[1]*(7*(5 + 4*b) + 2*(7 + 8*b)*c[2]))*d[1] 
	+ 8*(3*(7 + 8*b)*c[1]*c[1] + 3*(1 + b)*c[1]*(7 + 22*c[2]) + c[2]*(147 + 
	90*b - 81*c[2] - 34*b*c[2]))*d[2] + 21*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 10*(1 
	+ b)*d[1] - 2*(5 + 2*b)*d[2]) + 6*c[0]*(7*(5*(3 + 2*b + 2*(1 + b)*c[1]) - 
	2*(5 + 2*b)*c[2])* d[0] + 14*(5 + 5*b + 10*c[1] + 8*b*c[1] + 2*(1+b)*c[2])*
	 d[1] + 2*(7*(-5 - 2*b + 2*(1 + b)*c[1]) + 2*(49 + 30*b)*
	c[2])* d[2]))))/20160.;
		bf_mom[3] = -((a - b)*(2*a*a*(3*(35*c[0]*(2 + c[0]) - 70*(1 + c[0])
	*c[1] + 56*c[1]*c[1] - 28*(1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*
	(35*c[0]*(2 + c[0]) - 112*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 
	16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*(2 + c[0]) + 42*(1 + c[0])
	*c[1] - 48*c[1]*c[1] - 12*(15 + 15*c[0] - 11*c[1])*c[2]+68*c[2]*c[2])*d[2])
	 + a*(3*(7*(5*(-3 + 2*b)*c[0]*(2 + c[0]) + 20*(1 + c[0])*c[1] + 4*(-5+2*b)
	*c[1]*c[1]) - 28*((-5 + 6*b)*(1 + c[0]) - 2*c[1])*c[2] + 4*(-49 + 38*b)*
	c[2]*c[2])*d[0] + 6*(35*c[0]*(2 + c[0]) + 28*(-5 + 2*b)*(1 + c[0])*c[1] + 
	84*c[1]*c[1] + 4*(7 + 7*c[0] - 2*(7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])*d[1] 
	- 2*(21*(-5 + 6*b)*c[0]*(2 + c[0]) - 84*(1 + c[0])*c[1] + 12*(7 + 2*b)*c[1]
	*c[1] - 12*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*c[2] + 4*(-81 + 94*b)*c[2]
	*c[2])*d[2]) + b*(12*(7*(-5 + 4*b)*c[1]*c[1] + 7*(-1 + b)*c[1]*(5 + 2*c[2])
	 + c[2]*(35 - 49*c[2] + 2*b*(-7 + 15*c[2])))*d[0] + 24*(21*(-1 + b)*c[1]*
	c[1] + (-1 + b)*c[2]*(7 + 11*c[2]) + c[1]*(-7*(5 + 2*c[2]) +4*b*(7+4*c[2])))
	*d[1] + 8*(3*(-7 + 8*b)*c[1]*c[1] + 3*(-1 + b)*c[1]*(7 + 22*c[2]) + c[2]*
	(-147 + b*(90 - 34*c[2]) + 81*c[2]))*d[2] + 21*c[0]*c[0]*(5*(-3 + 2*b)*d[0]
	 + 10*(-1 + b)*d[1] + 2*(5 - 2*b)*d[2]) + 6*c[0]*(7*(5*(-3 + 2*b + 2*(-1+b)
	*c[1]) - 2*(-5 + 2*b)*c[2])* d[0] + 14*(-5 + 5*b - 10*c[1] + 8*b*c[1] + 2*
	(-1 + b)*c[2])* d[1] + 2*(7*(5 - 2*b + 2*(-1 + b)*c[1]) + 2*(-49 + 30*b)
	*c[2])*d[2]))))/20160.;
		bf_mom[4] = ((a - b)*(21*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)*d[0] 
	- 5*(a - b)*(a + b)*d[1] - 2*(-5 + a*a + 3*a*b + b*b)*d[2]) + 6*c[0]*(-7*
	(5*(-3 + a*a + a*b + b*b + (a - b)*(a + b)*c[1]) + 2*(-5 + a*a + 3*a*b+b*b)
	*c[2])*d[0] + 7*(-20*c[1] + 4*a*b*c[1] + a*a*(5 + 8*c[1] - 2*c[2]) + b*b*
	(-5 + 8*c[1] + 2*c[2]))*d[1] - 2*(35 + a*a*(-7 + 7*c[1] - 30*c[2]) + 98*c[2]
	 - b*b*(7 + 7*c[1] + 30*c[2]) - a*b*(21 + 38*c[2]))*d[2]) + 2*(-42*(5*c[1]
	*c[1] + c[2]*(5 + 7*c[2]))*d[0] - 84*c[1]*(-5 + 2*c[2])*d[1] + 6*a*b*
	(7*c[1]*c[1]*d[0] + c[2]*(21 + 19*c[2])*d[0] - 2*c[1]*(7 + 2*c[2])*d[1]) - 
	12*(7*c[1]*c[1] - c[2]*(49 + 27*c[2]))*d[2] - 4*a*b*(3*c[1]*c[1] + c[2]*(57
	 + 47*c[2]))*d[2] + a*a*(3*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 
	30*c[2]*c[2])*d[0] - 6*(7*c[1]*(4 + 3*c[1]) - (7 + 16*c[1])*c[2] + 11*c[2]
	*c[2])*d[1] + 2*(3*c[1]*(7 + 8*c[1]) - 6*(15 + 11*c[1])*c[2] - 34*c[2]*c[2])
	*d[2]) + b*b*(6*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2]) + 2*c[2]*(3*(7 + 
	15*c[2])*d[0] + 3*(-7 + 11*c[2])*d[1] - 2*(45 + 17*c[2])*d[2]) + 3*c[1]*(7*
	(-5 + 2*c[2])*d[0] + 8*(-7 + 4*c[2])*d[1] + 2*(-7 + 22*c[2])*d[2])))))/5040.;
		bf_mom[5] = ((a - b)*(2*a*a*(3*(7* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1]
	 + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*(-140 + 
	35*c[0]*c[0] + 84*c[1]*c[1] - 28*c[0]*(4*c[1] - c[2]) - 64*c[1]*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 
	30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2])* d[2]) + b*(3*(7*(5*(3 + 2*b)*(-4 
	+ c[0]*c[0]) + 20*(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*c[1]*c[1]) - 28*((5+2*b)
	*c[0] - 2*(1 + b)*c[1])*c[2] + 4*(49 + 30*b)*c[2]*c[2])*d[0] + 6*(35*(1+b)
	*(-4 + c[0]*c[0]) + 28*(5 + 4*b)*c[0]*c[1] + 84*(1 + b)*c[1]*c[1] + 4*(7*
	(1 + b)*c[0] + 2*(7 + 8*b)*c[1])*c[2] + 44*(1 + b)*c[2]*c[2])*d[1] + 2*
	(-21*(5 + 2*b)*(-4 + c[0]*c[0]) + 84*(1 + b)*c[0]*c[1] + 12*(7 + 8*b)*c[1]*
	c[1] + 12*((49 + 30*b)*c[0] + 22*(1 + b)*c[1])*c[2] - 4*(81 + 34*b)*c[2]*
	c[2])*d[2]) + a*(3*(7*(5*(3 + 2*b)*(-4 + c[0]*c[0]) - 20*c[0]*c[1] + 4*(5 
	+ 2*b)*c[1]*c[1]) - 28*((5 + 6*b)*c[0] + 2*c[1])*c[2] + 4*(49 + 38*b)*c[2]
	*c[2])*d[0] - 6*(35*(-4 + c[0]*c[0]) - 28*(5 + 2*b)*c[0]*c[1] + 84*c[1]*
	c[1] + 4*(7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*
	(5 + 6*b)*(-4 + c[0]*c[0]) + 84*c[0]*c[1] + 12*(-7 + 2*b)*c[1]*c[1] - 12*(
	(49 + 38*b)*c[0] - 22*c[1])*c[2] + 4*(81 + 94*b)*c[2]*c[2])*d[2])))/10080.;
		bf_mom[6] = ((a - b)*(21*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)*d[0] - 
	5*(a - b)*(a + b)*d[1] - 2*(-5 + a*a + 3*a*b + b*b)*d[2]) + 6*c[0]*(-7*(15 
	- 5*a*b + 5*a*a*(-1 + c[1]) - 5*b*b*(1 + c[1]) + 2*(-5 + a*a + 3*a*b + b*b)
	*c[2])*d[0] + 7*(-20*c[1] + 4*a*b*c[1] + a*a*(-5 + 8*c[1] - 2*c[2]) + b*b*
	(5 + 8*c[1] + 2*c[2]))*d[1] - 2*(-35 + a*b*(21 - 38*c[2]) + b*b*(7 - 7*c[1]
	 - 30*c[2]) + a*a*(7 + 7*c[1] - 30*c[2]) + 98*c[2])*d[2]) + 2*(-42*(5*c[1]
	*c[1] + c[2]*(-5 + 7*c[2]))*d[0] - 84*c[1]*(5 + 2*c[2])*d[1] + 6*a*b*(7*c[1]
	*c[1]*d[0] + c[2]*(-21 + 19*c[2])*d[0] - 2*c[1]*(-7 + 2*c[2])*d[1]) - 12*
	(7*c[1]*c[1] + (49 - 27*c[2])*c[2])*d[2] - 4*a*b*(3*c[1]*c[1] + c[2]*(-57 
	+ 47*c[2]))*d[2] + a*a*(6*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2]) + 2*c[2]*
	(3*(-7 + 15*c[2])*d[0] - 3*(7 + 11*c[2])*d[1] + 90*d[2] - 34*c[2]*d[2]) - 
	3*c[1]*(7*(5 + 2*c[2])*d[0] - 8*(7 + 4*c[2])*d[1] + 2*(7 + 22*c[2])*d[2])) 
	+ b*b*(6*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2]) + 2*c[2]*(3*(-7 + 15*c[2])
	*d[0] + 3*(7 + 11*c[2])*d[1] + 90*d[2] - 34*c[2]*d[2]) + 3*c[1]*(7*(5 + 
	2*c[2])*d[0] + 8*(7 + 4*c[2])*d[1] + 2*(7 + 22*c[2])*d[2])))))/5040.;
		bf_mom[7] = ((a - b)*(2*a*a*(3*(7* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1]
	 + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*(-140 + 
	35*c[0]*c[0] + 84*c[1]*c[1] - 28*c[0]*(4*c[1] - c[2]) - 64*c[1]*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 
	30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2])* d[2]) + b*(3*(7*(5*(-3 + 2*b)*(-4
	 + c[0]*c[0]) + 20*(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)*c[1]*c[1]) - 28*((-5 
	+ 2*b)*c[0] - 2*(-1 + b)*c[1])*c[2] + 4*(-49 + 30*b)*c[2]*c[2])*d[0] + 6*
	(35*(-1 + b)*(-4 + c[0]*c[0]) + 28*(-5 + 4*b)*c[0]*c[1] + 84*(-1 + b)*c[1]*
	c[1] + 4*(7*(-1 + b)*c[0] + 2*(-7 + 8*b)*c[1])*c[2] + 44*(-1+b)*c[2]*c[2])*
	d[1] - 2*(21*(-5 + 2*b)*(-4 + c[0]*c[0]) - 84*(-1 + b)*c[0]*c[1] - 12*(-7 
	+ 8*b)*c[1]*c[1] - 12*((-49 + 30*b)*c[0] + 22*(-1 + b)*c[1])*c[2] + 4*(-81 
	+ 34*b)*c[2]*c[2])*d[2]) + a*(3*(7*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*c[0]
	*c[1] + 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5 + 6*b)*c[0] - 2*c[1])*c[2] + 4*
	(-49 + 38*b)*c[2]*c[2])*d[0] + 6*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] - 8*
	(7 + 2*b)*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*((-5 + 2*b)*c[1] + c[2]))*d[1]
	 - 2*(21*(-5 + 6*b)*(-4 + c[0]*c[0]) - 84*c[0]*c[1] + 12*(7 + 2*b)*c[1]*c[1]
	 - 12*((-49 + 38*b)*c[0] + 22*c[1])*c[2] + 4*(-81 + 94*b)*c[2]*c[2])
	*d[2])))/10080.;
		bf_mom[8] = -((a - b)*(3*(35*(-3 + a*a + a*b + b*b)* (-4+c[0]*c[0])
	 - 70*(a - b)*(a + b)*c[0]*c[1] + 28*(-5 + 2*a*a + a*b + 2*b*b)*c[1]*c[1] - 
	28*((-5 + a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*c[2] + 4*(-49 + 
	15*a*a + 19*a*b + 15*b*b)*c[2]*c[2])* d[0] - 3*(35*(a - b)*(a + b)*(-4 + 
	c[0]*c[0]) - 56*(-5 + 2*a*a + a*b + 2*b*b)*c[0]*c[1] + 84*(a - b)*(a + b)
	*c[1]*c[1] + 4*(7*(a - b)*(a + b)*c[0] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1])
	*c[2] + 44*(a - b)*(a + b)*c[2]*c[2])*d[1] - 2*(21*(-5 + a*a + 3*a*b + b*b)
	*(-4 + c[0]*c[0]) + 42*(a - b)*(a + b)*c[0]*c[1] - 12*(-7 + 4*a*a - a*b + 
	4*b*b)*c[1]*c[1] - 12*((-49 + 15*a*a + 19*a*b + 15*b*b)*c[0] + 11*(-a*a + 
	b*b)*c[1])*c[2] + 4*(-81 + 17*a*a + 47*a*b + 17*b*b)*c[2]*c[2])*d[2]))/2520.;
		break;
	case 101:
	case 103:
		bf_mom[0] = -((a - b)*(2*a*a*(3*(35*c[0]*c[0] - 70*c[0]*(1 + c[1]) 
	+ 14*c[1]*(5 + 4*c[1]) - 28*(-1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 
	3*(35*(-2 + c[0])*c[0] - 112*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*
	c[0] - 16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*c[0] - 6*c[1]*(7 + 
	8*c[1]) + 6*c[0]*(-7 + 7*c[1] - 30*c[2]) + 12*(15 + 11*c[1])*c[2] + 68*c[2]*
	c[2])*d[2]) + a*(3*(7*(5*(-3 + 2*b)*(-2 + c[0])*c[0] + 20*(-1 + c[0])*c[1] 
	+ 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5 + 6*b)*(-1 + c[0]) - 2*c[1])*c[2] + 4*
	(-49 + 38*b)*c[2]*c[2])*d[0] + 6*(35*(-2 + c[0])*c[0] + 28*(-5 + 2*b)*(-1 +
	 c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0] - 2*(7 + 2*b)*c[1])*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(21*(-5 + 6*b)*(-2 + c[0])*c[0] - 84*(-1 + c[0])*c[1] 
	+ 12*(7 + 2*b)*c[1]*c[1] - 12*((-49 + 38*b)*(-1 + c[0]) + 22*c[1])*c[2] + 
	4*(-81 + 94*b)*c[2]*c[2])*d[2]) + b*(12*(7*(-5 + 4*b)*c[1]*c[1] + 7*(-1 + b)
	*c[1]*(-5 + 2*c[2]) + c[2]*(-7*(5 + 7*c[2]) + 2*b*(7 + 15*c[2])))*d[0] + 
	24*(21*(-1 + b)*c[1]*c[1] + (-1 + b)*c[2]*(-7 + 11*c[2]) + c[1]*(35 - 14*
	c[2] + 4*b*(-7 + 4*c[2])))*d[1] + 8*(3*(-7 + 8*b)*c[1]*c[1] + 3*(-1 + b)*
	c[1]*(-7 + 22*c[2]) + c[2]*(147 + 81*c[2] - 2*b*(45 + 17*c[2])))*d[2] + 21*
	c[0]*c[0]*(5*(-3 + 2*b)*d[0] + 10*(-1 + b)*d[1] + 2*(5 - 2*b)*d[2]) + 6*
	c[0]*(7*(5*(3 - 2*b + 2*(-1 + b)*c[1]) - 2*(-5 + 2*b)*c[2])* d[0] + 14*(5 - 
	5*b - 10*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2])* d[1] + 2*(7*(-5 + 2*b + 2*(-1 +
	 b)*c[1]) + 2*(-49 + 30*b)*c[2])*d[2]))))/20160.;
		bf_mom[1] = -((a - b)*(2*a*a*(3*(35*c[0]*(2 + c[0]) - 70*(1 + c[0])
	*c[1] + 56*c[1]*c[1] - 28*(1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*
	(35*c[0]*(2 + c[0]) - 112*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 
	16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*(2 + c[0]) + 42*(1 + c[0])*
	c[1] - 48*c[1]*c[1] - 12*(15 + 15*c[0] - 11*c[1])*c[2] + 68*c[2]*c[2])*d[2])
	 + a*(3*(7*(5*(-3 + 2*b)*c[0]*(2 + c[0]) + 20*(1 + c[0])*c[1] + 4*(-5 + 2*b)
	*c[1]*c[1]) - 28*((-5 + 6*b)*(1 + c[0]) - 2*c[1])*c[2] + 4*(-49 + 38*b)*c[2]
	*c[2])*d[0] + 6*(35*c[0]*(2 + c[0]) + 28*(-5 + 2*b)*(1 + c[0])*c[1] + 84*
	c[1]*c[1] + 4*(7 + 7*c[0] - 2*(7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*
	(21*(-5 + 6*b)*c[0]*(2 + c[0]) - 84*(1 + c[0])*c[1] + 12*(7 + 2*b)*c[1]*c[1]
	 - 12*((-49 + 38*b)*(1 + c[0]) + 22*c[1])*c[2] + 4*(-81 + 94*b)*c[2]*c[2])*
	d[2]) + b*(12*(7*(-5 + 4*b)*c[1]*c[1] + 7*(-1 + b)*c[1]*(5 + 2*c[2]) + c[2]*
	(35 - 49*c[2] + 2*b*(-7 + 15*c[2])))*d[0] + 24*(21*(-1 + b)*c[1]*c[1] + 
	(-1 + b)*c[2]*(7 + 11*c[2]) + c[1]*(-7*(5 + 2*c[2]) + 4*b*(7 + 4*c[2])))*
	d[1] + 8*(3*(-7 + 8*b)*c[1]*c[1] + 3*(-1 + b)*c[1]*(7 + 22*c[2]) + c[2]*
	(-147 + b*(90 - 34*c[2]) + 81*c[2]))*d[2] + 21*c[0]*c[0]*(5*(-3 + 2*b)*d[0] 
	+ 10*(-1 + b)*d[1] + 2*(5 - 2*b)*d[2]) + 6*c[0]*(7*(5*(-3 + 2*b + 2*(-1 + b)
	*c[1]) - 2*(-5 + 2*b)*c[2])* d[0] + 14*(-5 + 5*b - 10*c[1] + 8*b*c[1] + 2*
	(-1 + b)*c[2])* d[1] + 2*(7*(5 - 2*b + 2*(-1 + b)*c[1]) + 2*(-49 + 30*b)
	*c[2])*d[2]))))/20160.;
		bf_mom[2] = -((a - b)*(2*a*a*(3*(35*c[0]*(2 + c[0]) - 70*(1 + c[0])
	*c[1] + 56*c[1]*c[1] - 28*(1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*
	(35*c[0]*(2 + c[0]) - 112*(1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(7 + 7*c[0] - 
	16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*(2 + c[0]) + 42*(1 + c[0])*
	c[1] - 48*c[1]*c[1] - 12*(15 + 15*c[0] - 11*c[1])*c[2] + 68*c[2]*c[2])*d[2])
	 + a*(3*(7*(5*(3 + 2*b)*c[0]*(2 + c[0]) - 20*(1 + c[0])*c[1] + 4*(5 + 2*b)*
	c[1]*c[1]) - 28*((5 + 6*b)*(1 + c[0]) + 2*c[1])*c[2] + 4*(49 + 38*b)*c[2]*
	c[2])*d[0] - 6*(35*c[0]*(2 + c[0]) - 28*(5 + 2*b)*(1 + c[0])*c[1] + 84*c[1]*
	c[1] + 4*(7 + 7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])* d[1] - 2*
	(21*(5 + 6*b)*c[0]*(2 + c[0]) + 84*(1 + c[0])*c[1] + 12*(-7 + 2*b)*c[1]*c[1]
	 - 12*((49 + 38*b)*(1 + c[0]) - 22*c[1])*c[2] + 4*(81 + 94*b)*c[2]*c[2])*
	d[2]) + b*(12*(7*(5 + 4*b)*c[1]*c[1] + 7*(1 + b)*c[1]*(5 + 2*c[2]) + c[2]*
	(-7*(5 + 2*b) + (49 + 30*b)*c[2]))*d[0] + 24*(21*(1 + b)*c[1]*c[1] + (1 + b)
	*c[2]*(7 + 11*c[2]) + c[1]*(7*(5 + 4*b) + 2*(7 + 8*b)*c[2]))*d[1] + 8*(3*(7 
	+ 8*b)*c[1]*c[1] + 3*(1 + b)*c[1]*(7 + 22*c[2]) + c[2]*(147 + 90*b - 81*c[2]
	 - 34*b*c[2]))*d[2] + 21*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 10*(1 + b)*d[1] - 2*
	(5 + 2*b)*d[2]) + 6*c[0]*(7*(5*(3 + 2*b + 2*(1 + b)*c[1]) - 2*(5+2*b)*c[2])
	* d[0] + 14*(5 + 5*b + 10*c[1] + 8*b*c[1] + 2*(1 + b)*c[2])* d[1] + 2*(7*
	(-5 - 2*b + 2*(1 + b)*c[1]) + 2*(49 + 30*b)*c[2])* d[2]))))/20160.;
		bf_mom[3] = -((a - b)*(2*a*a*(3*(35*c[0]*c[0] - 70*c[0]*(1 + c[1]) 
	+ 14*c[1]*(5 + 4*c[1]) - 28*(-1 + c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 
	3*(35*(-2 + c[0])*c[0] - 112*(-1 + c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*
	c[0] - 16*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*c[0]*c[0] - 6*c[1]*(7 + 
	8*c[1]) + 6*c[0]*(-7 + 7*c[1] - 30*c[2]) + 12*(15 + 11*c[1])*c[2] + 68*c[2]*
	c[2])*d[2]) + a*(3*(7*(5*(3 + 2*b)*(-2 + c[0])*c[0] - 20*(-1 + c[0])*c[1] 
	+ 4*(5 + 2*b)*c[1]*c[1]) - 28*((5 + 6*b)*(-1 + c[0]) + 2*c[1])*c[2] + 4*
	(49 + 38*b)*c[2]*c[2])*d[0] - 6*(35*(-2 + c[0])*c[0] - 28*(5 + 2*b)*(-1 + 
	c[0])*c[1] + 84*c[1]*c[1] + 4*(-7 + 7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(21*(5 + 6*b)*(-2 + c[0])*c[0] + 84*(-1 + c[0])*c[1] 
	+ 12*(-7 + 2*b)*c[1]*c[1] - 12*((49 + 38*b)*(-1 + c[0]) - 22*c[1])*c[2] + 
	4*(81 + 94*b)*c[2]*c[2])*d[2]) + b*(12*(7*(5 + 4*b)*c[1]*c[1] + 7*(1 + b)*
	c[1]*(-5 + 2*c[2]) + c[2]*(35 + 14*b + 49*c[2] + 30*b*c[2]))*d[0] + 24*(21*
	(1 + b)*c[1]*c[1] + (1 + b)*c[2]*(-7 + 11*c[2]) + c[1]*(-7*(5 + 4*b) + 2*
	(7 + 8*b)*c[2]))*d[1] + 8*(3*(7 + 8*b)*c[1]*c[1] + 3*(1 + b)*c[1]*(-7 + 22*
	c[2]) - c[2]*(147 + 90*b + 81*c[2] + 34*b*c[2]))*d[2] + 21*c[0]*c[0]*(5*
	(3 + 2*b)*d[0] + 10*(1 + b)*d[1] - 2*(5 + 2*b)*d[2]) + 6*c[0]*(7*(5*(-3 - 
	2*b + 2*(1 + b)*c[1]) - 2*(5 + 2*b)*c[2])* d[0] + 14*(-5 - 5*b + 10*c[1] + 
	8*b*c[1] + 2*(1 + b)*c[2])* d[1] + 2*(7*(5 + 2*b + 2*(1 + b)*c[1]) + 2*(49 
	+ 30*b)*c[2])* d[2]))))/20160.;
		bf_mom[4] = ((a - b)*(2*a*a*(3*(7* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1]
	 + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*(-140 + 
	35*c[0]*c[0] + 84*c[1]*c[1] - 28*c[0]*(4*c[1] - c[2]) - 64*c[1]*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 
	30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2])* d[2]) + b*(3*(7*(5*(-3 + 2*b)*(-4 
	+ c[0]*c[0]) + 20*(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)*c[1]*c[1]) - 28*((-5 + 
	2*b)*c[0] - 2*(-1 + b)*c[1])*c[2] + 4*(-49 + 30*b)*c[2]*c[2])*d[0] + 6*(35*
	(-1 + b)*(-4 + c[0]*c[0]) + 28*(-5 + 4*b)*c[0]*c[1] + 84*(-1 + b)*c[1]*c[1]
	 + 4*(7*(-1 + b)*c[0] + 2*(-7 + 8*b)*c[1])*c[2] + 44*(-1 + b)*c[2]*c[2])*
	d[1] - 2*(21*(-5 + 2*b)*(-4 + c[0]*c[0]) - 84*(-1 + b)*c[0]*c[1] - 12*(-7 
	+ 8*b)*c[1]*c[1] - 12*((-49 + 30*b)*c[0] + 22*(-1 + b)*c[1])*c[2] + 4*(-81
	 + 34*b)*c[2]*c[2])*d[2]) + a*(3*(7*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*
	c[0]*c[1] + 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5 + 6*b)*c[0] - 2*c[1])*c[2] + 
	4*(-49 + 38*b)*c[2]*c[2])*d[0] + 6*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] - 8*
	(7 + 2*b)*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*((-5 + 2*b)*c[1] + c[2]))*d[1] 
	- 2*(21*(-5 + 6*b)*(-4 + c[0]*c[0]) - 84*c[0]*c[1] + 12*(7 + 2*b)*c[1]*c[1]
	 - 12*((-49 + 38*b)*c[0] + 22*c[1])*c[2] + 4*(-81 + 94*b)*c[2]*c[2])
	*d[2])))/10080.;
		bf_mom[5] = ((a - b)*(21*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)*d[0] - 
	5*(a - b)*(a + b)*d[1] - 2*(-5 + a*a + 3*a*b + b*b)*d[2]) + 6*c[0]*(-7*(15 
	- 5*a*b + 5*a*a*(-1 + c[1]) - 5*b*b*(1 + c[1]) + 2*(-5 + a*a + 3*a*b + b*b)
	*c[2])*d[0] + 7*(-20*c[1] + 4*a*b*c[1] + a*a*(-5 + 8*c[1] - 2*c[2]) + b*b*
	(5 + 8*c[1] + 2*c[2]))*d[1] - 2*(-35 + a*b*(21 - 38*c[2]) + b*b*(7 - 7*c[1]
	 - 30*c[2]) + a*a*(7 + 7*c[1] - 30*c[2]) + 98*c[2])*d[2]) + 2*(-42*(5*c[1]*
	c[1] + c[2]*(-5 + 7*c[2]))*d[0] - 84*c[1]*(5 + 2*c[2])*d[1] + 6*a*b*(7*c[1]*
	c[1]*d[0] + c[2]*(-21 + 19*c[2])*d[0] - 2*c[1]*(-7 + 2*c[2])*d[1]) - 12*(7*
	c[1]*c[1] + (49 - 27*c[2])*c[2])*d[2] - 4*a*b*(3*c[1]*c[1] + c[2]*(-57 + 47*
	c[2]))*d[2] + a*a*(6*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2]) + 2*c[2]*(3*
	(-7 + 15*c[2])*d[0] - 3*(7 + 11*c[2])*d[1] + 90*d[2] - 34*c[2]*d[2]) - 3*
	c[1]*(7*(5 + 2*c[2])*d[0] - 8*(7 + 4*c[2])*d[1] + 2*(7 + 22*c[2])*d[2])) + 
	b*b*(6*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2]) + 2*c[2]*(3*(-7 + 15*c[2])*
	d[0] + 3*(7 + 11*c[2])*d[1] + 90*d[2] - 34*c[2]*d[2]) + 3*c[1]*(7*(5 + 2*
	c[2])*d[0] + 8*(7 + 4*c[2])*d[1] + 2*(7 + 22*c[2])*d[2])))))/5040.;
		bf_mom[6] = ((a - b)*(2*a*a*(3*(7* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1]
	 + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2])*d[0] - 3*(-140 + 
	35*c[0]*c[0] + 84*c[1]*c[1] - 28*c[0]*(4*c[1] - c[2]) - 64*c[1]*c[2] + 44*
	c[2]*c[2])* d[1] - 2*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 
	30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2])* d[2]) + b*(3*(7*(5*(3 + 2*b)*(-4 
	+ c[0]*c[0]) + 20*(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*c[1]*c[1]) - 28*((5 + 2*b)
	*c[0] - 2*(1 + b)*c[1])*c[2] + 4*(49 + 30*b)*c[2]*c[2])*d[0] + 6*(35*(1 + b)
	*(-4 + c[0]*c[0]) + 28*(5 + 4*b)*c[0]*c[1] + 84*(1 + b)*c[1]*c[1] + 4*(7*
	(1 + b)*c[0] + 2*(7 + 8*b)*c[1])*c[2] + 44*(1 + b)*c[2]*c[2])*d[1] + 2*(-21*
	(5 + 2*b)*(-4 + c[0]*c[0]) + 84*(1 + b)*c[0]*c[1] + 12*(7 + 8*b)*c[1]*c[1] 
	+ 12*((49 + 30*b)*c[0] + 22*(1 + b)*c[1])*c[2] - 4*(81 + 34*b)*c[2]*c[2])
	*d[2]) + a*(3*(7*(5*(3 + 2*b)*(-4 + c[0]*c[0]) - 20*c[0]*c[1] + 4*(5 + 2*b)
	*c[1]*c[1]) - 28*((5 + 6*b)*c[0] + 2*c[1])*c[2] + 4*(49 + 38*b)*c[2]*c[2])
	*d[0] - 6*(35*(-4 + c[0]*c[0]) - 28*(5 + 2*b)*c[0]*c[1] + 84*c[1]*c[1] + 4*
	(7*c[0] + 2*(-7 + 2*b)*c[1])*c[2] + 44*c[2]*c[2])*d[1] - 2*(21*(5 + 6*b)*
	(-4 + c[0]*c[0]) + 84*c[0]*c[1] + 12*(-7 + 2*b)*c[1]*c[1] - 12*((49 + 38*b)
	*c[0] - 22*c[1])*c[2] + 4*(81 + 94*b)*c[2]*c[2])*d[2])))/10080.;
		bf_mom[7] = ((a - b)*(21*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)*d[0] - 
	5*(a - b)*(a + b)*d[1] - 2*(-5 + a*a + 3*a*b + b*b)*d[2]) + 6*c[0]*(-7*(5*
	(-3 + a*a + a*b + b*b + (a - b)*(a + b)*c[1]) + 2*(-5 + a*a + 3*a*b + b*b)*
	c[2])*d[0] + 7*(-20*c[1] + 4*a*b*c[1] + a*a*(5 + 8*c[1] - 2*c[2]) + b*b*
	(-5 + 8*c[1] + 2*c[2]))*d[1] - 2*(35 + a*a*(-7 + 7*c[1] - 30*c[2]) + 98*c[2]
	 - b*b*(7 + 7*c[1] + 30*c[2]) - a*b*(21 + 38*c[2]))*d[2]) + 2*(-42*(5*c[1]*
	c[1] + c[2]*(5 + 7*c[2]))*d[0] - 84*c[1]*(-5 + 2*c[2])*d[1] + 6*a*b*(7*c[1]*
	c[1]*d[0] + c[2]*(21 + 19*c[2])*d[0] - 2*c[1]*(7 + 2*c[2])*d[1]) - 12*(7*
	c[1]*c[1] - c[2]*(49 + 27*c[2]))*d[2] - 4*a*b*(3*c[1]*c[1] + c[2]*(57 + 47*
	c[2]))*d[2] + a*a*(3*(7*c[1]*(5 + 4*c[1]) - 14*(-1 + c[1])*c[2] + 30*c[2]*
	c[2])*d[0] - 6*(7*c[1]*(4 + 3*c[1]) - (7 + 16*c[1])*c[2] + 11*c[2]*c[2])*
	d[1] + 2*(3*c[1]*(7 + 8*c[1]) - 6*(15 + 11*c[1])*c[2] - 34*c[2]*c[2])*d[2])
	 + b*b*(6*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2]) + 2*c[2]*(3*(7 + 15*c[2])*
	d[0] + 3*(-7 + 11*c[2])*d[1] - 2*(45 + 17*c[2])*d[2]) + 3*c[1]*(7*(-5 + 2*
	c[2])*d[0] + 8*(-7 + 4*c[2])*d[1] + 2*(-7 + 22*c[2])*d[2])))))/5040.;
		bf_mom[8] = -((a - b)*(3*(35*(-3 + a*a + a*b + b*b)* (-4+c[0]*c[0])
	 - 70*(a - b)*(a + b)*c[0]*c[1] + 28*(-5 + 2*a*a + a*b + 2*b*b)*c[1]*c[1] - 
	28*((-5 + a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*c[2] + 4*(-49 + 
	15*a*a + 19*a*b + 15*b*b)*c[2]*c[2])* d[0] - 3*(35*(a - b)*(a + b)*(-4 + 
	c[0]*c[0]) - 56*(-5 + 2*a*a + a*b + 2*b*b)*c[0]*c[1] + 84*(a - b)*(a + b)*
	c[1]*c[1] + 4*(7*(a - b)*(a + b)*c[0] - 4*(-7 + 4*a*a - a*b + 4*b*b)*c[1])*
	c[2] + 44*(a - b)*(a + b)*c[2]*c[2])*d[1] - 2*(21*(-5 + a*a + 3*a*b + b*b)*
	(-4 + c[0]*c[0]) + 42*(a - b)*(a + b)*c[0]*c[1] - 12*(-7 + 4*a*a - a*b + 
	4*b*b)*c[1]*c[1] - 12*((-49 + 15*a*a + 19*a*b + 15*b*b)*c[0] + 11*(-a*a + 
	b*b)*c[1])*c[2] + 4*(-81 + 17*a*a + 47*a*b+17*b*b)*c[2]*c[2])*d[2]))/2520.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
case 4:
switch (interface_type)
	{
	case 0:
 		memset(bf_mom, 0, sizeof(double)*npts);
		break;
	case 2:
	case 3:
		bf_mom[0] = -(pow(-1 + a,2)*(231*c[0]*c[0]* (5*(1 + 2*a)*d[0] - 6*
	d[2] - 2*a*(5*d[1] + 2*d[2] - 3*d[3])) + 22*c[0]*(-21*(5 + 6*c[2] + 2*a*
	(5 + 5*c[1] + 2*c[2] - 3*c[3]))* d[0] + 6*(7*(2*c[1] + a*(5 + 8*c[1] - 2*
	c[2])) - 6*(3 + 4*a)*c[3])*d[1] + 6*(21 + 38*c[2] + 2*a*(7 - 7*c[1] + 30*
	c[2] - 15*c[3]))*d[2] - 2*(9*(6*c[1] + a*(7 + 8*c[1] + 10*c[2])) - 2*(53 + 
	100*a)*c[3])* d[3]) + 4*(693*c[2]*d[0] + 627*c[2]*c[2]*d[0] + 583*c[3]*c[3]
	*d[0] + 594*c[3]*d[1] + 572*c[2]*c[3]*d[1] + 33*c[1]*c[1]*(7*d[0] - 2*d[2])
	 - 1254*c[2]*d[2] - 1034*c[2]*c[2]*d[2] - 642*c[3]*c[3]*d[2] - 2*(583 + 
	642*c[2])*c[3]*d[3] - 22*c[1]*(3*(7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 26*
	d[2]) - (27 + 26*c[2])*d[3]) + a*(c[3]*(11*(-63 + 100*c[3])*d[0] - 2*c[3]*
	(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*d[1] + 45*d[2] - 100*d[3])) + 
	66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0] 
	- 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*(-33*(-7 + 15*c[3])*d[0] + 11*(21 + 
	64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2] + 495*d[3] - 656*c[3]*d[3]) + 11*
	c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0] + 2*(6*(-14 + 8*c[2] + 3*c[3])*d[1]
	 + (21 - 66*c[2] + 64*c[3])*d[2] + 2*(18 + 32*c[2] - 47*c[3])
	*d[3]))))))/221760.;
		bf_mom[1] = -((-1 + a)*(231*c[0]*c[0]* (5*(5 + a*(5 + 2*a))*d[0] - 
	10*(-2 + a + a*a)*d[1] - 2*(7 + a*(11 + 2*a))*d[2] + 6*(-2 + a + a*a)*d[3])
	 - 22*c[0]*(21*(25 - 20*c[1] + 14*c[2] + a*(25 + 10*c[1] + 22*c[2] + 2*a*
	(5 + 5*c[1] + 2*c[2] - 3*c[3]) - 6*c[3]) + 12*c[3])* d[0] - 6*(7*(2*(-5 + 
	9*c[1] + 2*c[2]) + a*(5 + 5*a + 14*c[1] + 8*a*c[1] - 2*(1 + a)*c[2])) - 
	6*(11 + a*(13 + 4*a))*c[3])*d[1] + 6*(-49 - 28*c[1] - 158*c[2] - 60*c[3] + 
	a*(-77 - 14*a + 14*c[1] + 14*a*c[1] - 174*c[2] - 60*a*c[2] + 30*(1 + a)*
	c[3]))*d[2] + 2*(9*(-14 + 22*c[1] - 20*c[2] + a*(7 + 7*a + 26*c[1] + 8*a*
	c[1] + 10*(1 + a)*c[2])) - 2*(253 + a*(259 + 100*a))*c[3])*d[3]) + 4*(11*
	c[2]*c[2]*(3*(79 + 87*a + 30*a*a)*d[0] - 66*(-2 + a + a*a)*d[1] - 2*(115 + 
	a*(175 + 34*a))*d[2] + 14*(-2 + a + a*a)*d[3]) + c[2]*(-33*(-7*(7 + a*(11 + 
	2*a)) + 30*(-2 + a + a*a)*c[3])*d[0] + 22*((21*(-2 + a + a*a) + 2*(77 + a*
	(71 + 32*a))*c[3])*d[1] + (-3*(79 + 87*a + 30*a*a) + 14*(-2 + a + a*a)*
	c[3])*d[2]) - 2*(-495*(-2 + a + a*a) + 2*(977 + a*(1291 + 328*a))*c[3])*
	d[3]) - 11*c[1]*(3*(7*(-2 + a + a*a)*(-5 + 2*c[2]) + 6*(11 + a*(13 + 4*a))
	*c[3])*d[0] - 6*(-63 + 30*c[2] - 12*c[3] + a*(-49 - 28*a + 10*c[2] + 16*a*
	c[2] + 6*(1 + a)*c[3]))* d[1] + 2*(3*(-2 + a + a*a)*(-7 + 22*c[2]) - 2*
	(77 + a*(71 + 32*a))*c[3])*d[2] - 2*(99 + 154*c[2] + 188*c[3] + a*(117 + 
	36*a + 142*c[2] + 64*a*c[2] - 94*(1 + a)*c[3]))* d[3]) + 33*c[1]*c[1]* (7*
	(9 + a*(7 + 4*a))*d[0] + 2*(-21*(-2 + a + a*a)*d[1] + (15 + a*(5 + 8*a))
	*d[2] + 3*(-2 + a + a*a)*d[3])) + c[3]*(11*(-63*(-2 + a + a*a) + (253 + a*
	(259 + 100*a))*c[3])*d[0] + 2*(-11*(-9*(11 + a*(13 + 4*a)) + 47*(-2 + a + 
	a*a)*c[3])*d[1] - (-495*(-2 + a + a*a) + (977 + a*(1291 + 328*a))*c[3])
	*d[2] + (-11*(253 + a*(259 + 100*a)) + 531*(-2 + a + a*a)*c[3])
	*d[3])))))/221760.;
		bf_mom[2] = -((-1 + a)*(231*c[0]*c[0]* (5*(5 + a*(5 + 2*a))*d[0] - 
	10*(-2 + a + a*a)*d[1] - 2*(7 + a*(11 + 2*a))*d[2] + 6*(-2 + a + a*a)*d[3])
	 - 22*c[0]*(21*(-25 - 20*c[1] + 14*c[2] + a*(-25 + 10*c[1] + 22*c[2] + 2*a*
	(-5 + 5*c[1] + 2*c[2] - 3*c[3]) - 6*c[3]) + 12*c[3])* d[0] - 6*(7*(2*(5 + 
	9*c[1] + 2*c[2]) + a*(-5 - 5*a + 14*c[1] + 8*a*c[1] - 2*(1 + a)*c[2])) - 
	6*(11 + a*(13 + 4*a))*c[3])*d[1] + 6*(49 - 28*c[1] - 158*c[2] - 60*c[3] + 
	a*(77 + 14*a + 14*c[1] + 14*a*c[1] - 174*c[2] - 60*a*c[2] + 30*(1+a)*c[3]))
	*d[2] + 2*(9*(14 + 22*c[1] - 20*c[2] + a*(-7 - 7*a + 26*c[1] + 8*a*c[1] + 
	10*(1 + a)*c[2])) - 2*(253 + a*(259 + 100*a))*c[3])*d[3]) + 4*(11*c[2]*c[2]
	*(3*(79 + 87*a + 30*a*a)*d[0] - 66*(-2 + a + a*a)*d[1] - 2*(115 + a*(175 + 
	34*a))*d[2] + 14*(-2 + a + a*a)*d[3]) + c[3]*(11*(63*(-2 + a + a*a) + (253 
	+ a*(259 + 100*a))*c[3])* d[0] - 22*(9*(11 + a*(13 + 4*a)) + 47*(-2 + a + 
	a*a)*c[3])*d[1] - 2*(495*(-2 + a + a*a) + (977 + a*(1291 + 328*a))*c[3])
	*d[2] + 2*(11*(253 + a*(259 + 100*a)) + 531*(-2 + a + a*a)*c[3])*d[3]) + 
	c[2]*(-33*(7*(7 + a*(11 + 2*a)) + 30*(-2 + a + a*a)*c[3])*d[0] + 22*((-21*
	(-2 + a + a*a) + 2*(77 + a*(71 + 32*a))*c[3])*d[1] + (3*(79 + 87*a + 30*a*a)
	 + 14*(-2 + a + a*a)*c[3])*d[2]) - 2*(495*(-2 + a + a*a) + 2*(977 + a*(1291
	 + 328*a))*c[3])*d[3]) - 11*c[1]*(3*(7*(-2 + a + a*a)*(5 + 2*c[2]) + 6*
	(11 + a*(13 + 4*a))*c[3])*d[0] - 6*(63 + 30*c[2] - 12*c[3] + a*(49 + 28*a 
	+ 10*c[2] + 16*a*c[2] + 6*(1 + a)*c[3]))*d[1] + 2*(3*(-2 + a + a*a)*(7 + 
	22*c[2]) - 2*(77 + a*(71 + 32*a))*c[3])* d[2] - 2*(-99 + 154*c[2] + 188*c[3]
	 + a*(-117 - 36*a + 142*c[2] + 64*a*c[2] - 94*(1 + a)*c[3]))* d[3]) + 33*
	c[1]*c[1]* (7*(9 + a*(7 + 4*a))*d[0] + 2*(-21*(-2 + a + a*a)*d[1] + (15 + 
	a*(5 + 8*a))*d[2] + 3*(-2 + a + a*a)*d[3])))))/221760.;
		bf_mom[3] = -(pow(-1 + a,2)*(231*c[0]*c[0]* (5*(1 + 2*a)*d[0] - 6*
	d[2] - 2*a*(5*d[1] + 2*d[2] - 3*d[3])) + 22*c[0]*(-21*(-5 + 6*c[2] + 2*a*
	(-5 + 5*c[1] + 2*c[2] - 3*c[3]))* d[0] + 6*(7*(2*c[1] + a*(-5 + 8*c[1] - 2*
	c[2])) - 6*(3 + 4*a)*c[3])*d[1] - 6*(21 - 38*c[2] + 2*a*(7 + 7*c[1] - 30*
	c[2] + 15*c[3]))*d[2] - 2*(9*(6*c[1] + a*(-7 + 8*c[1] + 10*c[2])) - 2*(53 
	+ 100*a)*c[3])* d[3]) + 4*(-693*c[2]*d[0] + 627*c[2]*c[2]*d[0] + 583*c[3]*
	c[3]*d[0] - 594*c[3]*d[1] + 572*c[2]*c[3]*d[1] + 33*c[1]*c[1]*(7*d[0] - 2*
	d[2]) + 1254*c[2]*d[2] - 1034*c[2]*c[2]*d[2] - 642*c[3]*c[3]*d[2] - 2*(-583
	 + 642*c[2])*c[3]*d[3] - 22*c[1]*(3*(-7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 
	26*d[2]) + (27 - 26*c[2])*d[3]) + a*(-22*c[2]*(3*(7 + 15*c[3])*d[0] + (21 
	- 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]*(495 + 656*c[3])*d[3] + 
	66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0] 
	- 33*d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 100*c[3])*d[0] - 22*(36 + 
	47*c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*(1100 + 531*c[3])*d[3]) - 11*
	c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*(14 + 8*c[2] + 3*c[3])*d[1] + 
	(-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 32*c[2] - 47*c[3])
	*d[3]))))))/221760.;
		bf_mom[4] = (pow(-1 + a,2)*(231*c[0]*c[0]* (5*(2 + a)*d[0] - 5*
	(1 + a)*d[1] - 2*(4 + a)*d[2] + 3*(1 + a)*d[3]) + 22*c[0]*(-21*(10 + 5*a + 
	5*c[1] + 5*a*c[1] + 8*c[2] + 2*a*c[2] - 3*(1 + a)*c[3])*d[0] + 3*(7*(5 + 
	5*a + 12*c[1] + 8*a*c[1] - 2*(1 + a)*c[2]) - 12*(5 + 2*a)*c[3])*d[1] - 6*
	(-28 - 7*a + 7*c[1] + 7*a*c[1] - 68*c[2] - 30*a*c[2] + 15*(1 + a)*c[3])*d[2]
	 - (9*(7 + 7*a + 20*c[1] + 8*a*c[1] + 10*(1 + a)*c[2]) - 4*(103+50*a)*c[3])
	*d[3]) + 2*(22*c[2]*((21*(4 + a) - 45*(1 + a)*c[3])*d[0] + (21*(1 + a) + 4*
	(29 + 16*a)*c[3])*d[1] + 2*(-3*(34 + 15*a) + 7*(1 + a)*c[3])*d[2]) - 2*c[2]*
	(-495*(1 + a) + 4*(485 + 164*a)*c[3])*d[3] + 66*c[1]*c[1]*(7*(3 + 2*a)*d[0]
	 - 21*(1 + a)*d[1] + 2*(3 + 4*a)*d[2] + 3*(1 + a)*d[3]) + 22*c[2]*c[2]*(3*
	(34 + 15*a)*d[0] - 33*(1 + a)*d[1] - 2*(64 + 17*a)*d[2] + 7*(1 + a)*d[3]) 
	+ 11*c[1]*(-3*(7*(1 + a)*(-5 + 2*c[2]) + 12*(5 + 2*a)*c[3])*d[0] + 2*(6*
	(-21 - 14*a + 6*c[2] + 8*a*c[2] + 3*(1 + a)*c[3])*d[1] + (-3*(1 + a)*(-7 + 
	22*c[2]) + 4*(29 + 16*a)*c[3])*d[2] + 2*(45 + 18*a + 58*c[2] + 32*a*c[2] 
	- 47*(1 + a)*c[3])*d[3])) + c[3]*(11*(-63*(1 + a) + 2*(103 + 50*a)*c[3])
	*d[0] + 2*(-11*(-18*(5 + 2*a) + 47*(1 + a)*c[3])*d[1] + (495*(1 + a) - 2*
	(485 + 164*a)*c[3])*d[2] + (-22*(103 + 50*a) + 531*(1 + a)*c[3])
	*d[3])))))/55440.;
		bf_mom[5] = ((-1 + a)*(11*(21*(5*(5 + a*(5 + 2*a))*(-4 + c[0]*c[0]) 
	- 20*(-2 + a + a*a)*c[0]*c[1] + 4*(9 + a*(7 + 4*a))*c[1]*c[1]) - 84*((7 + 
	a*(11 + 2*a))*c[0] + 2*(-2 + a + a*a)*c[1])* c[2] + 12*(79 + 87*a + 30*a*a)
	*c[2]*c[2] + 36*(7*(-2 + a + a*a)*c[0] - 2*(11 + a*(13 + 4*a))*c[1] - 10*
	(-2 + a + a*a)*c[2])* c[3] + 4*(253 + a*(259 + 100*a))*c[3]*c[3])*d[0] + 
	2*(-11*(3*(35*(-2 + a + a*a)*(-4 + c[0]*c[0]) - 28*(9 + a*(7 + 4*a))*c[0]*
	c[1] + 84*(-2 + a + a*a)*c[1]*c[1] + 4*(7*(-2 + a + a*a)*c[0] - 2*(15 + a*
	(5 + 8*a))*c[1])*c[2] + 44*(-2 + a + a*a)*c[2]*c[2]) + 4*(9*(11 + a*(13 + 
	4*a))*c[0] - 18*(-2 + a + a*a)*c[1] - 2*(77 + a*(71 + 32*a))*c[2])*c[3] + 
	188*(-2 + a + a*a)*c[3]*c[3])*d[1] + (33*(-7*(7 + a*(11 + 2*a))*(-4 + c[0]*
	c[0]) - 28*(-2 + a + a*a)*c[0]*c[1] + 4*(15 + a*(5 + 8*a))*c[1]*c[1]) + 
	132*((79 + 87*a + 30*a*a)*c[0] - 22*(-2 + a + a*a)*c[1])*c[2] - 44*(115 + 
	a*(175 + 34*a))*c[2]*c[2] - 44*(45*(-2 + a + a*a)*c[0] - 2*(77 + a*(71 + 
	32*a))*c[1] - 14*(-2 + a + a*a)*c[2])*c[3] - 4*(977 + a*(1291 + 328*a))*
	c[3]*c[3])*d[2] + (99*(7*(-2 + a + a*a)*(-4 + c[0]*c[0]) - 4*(11 + a*(13 + 
	4*a))*c[0]*c[1] + 4*(-2 + a + a*a)*c[1]*c[1]) - 44*(45*(-2 + a + a*a)*c[0] 
	- 2*(77 + a*(71 + 32*a))*c[1])*c[2] + 308*(-2 + a + a*a)*c[2]*c[2] + 4*(11*
	(253 + a*(259 + 100*a))*c[0] - 1034*(-2 + a + a*a)*c[1] - 2*(977 + a*(1291 
	+ 328*a))*c[2])*c[3] + 2124*(-2 + a + a*a)*c[3]*c[3])*d[3])))/110880.;
		bf_mom[6] = (pow(-1 + a,2)*(231*c[0]*c[0]* (5*(2 + a)*d[0] - 5*(1 
	+ a)*d[1] - 2*(4 + a)*d[2] + 3*(1 + a)*d[3]) + 22*c[0]*(-21*(-10 - 5*a + 
	5*c[1] + 5*a*c[1] + 8*c[2] + 2*a*c[2] - 3*(1 + a)*c[3])*d[0] + 3*(7*(-5 - 
	5*a + 12*c[1] + 8*a*c[1] - 2*(1 + a)*c[2]) - 12*(5 + 2*a)*c[3])*d[1] - 6*
	(28 + 7*a + 7*c[1] + 7*a*c[1] - 68*c[2] - 30*a*c[2] + 15*(1 + a)*c[3])*d[2] 
	- (9*(-7 - 7*a + 20*c[1] + 8*a*c[1] + 10*(1 + a)*c[2]) - 4*(103+50*a)*c[3])
	*d[3]) + 2*(22*c[2]*(-3*(7*(4 + a) + 15*(1 + a)*c[3])*d[0] + (-21*(1 + a) 
	+ 4*(29 + 16*a)*c[3])*d[1] + 2*(102 + 45*a + 7*c[3] + 7*a*c[3])*d[2]) - 
	2*c[2]*(495*(1 + a) + 4*(485 + 164*a)*c[3])*d[3] + 66*c[1]*c[1]*(7*(3+2*a)
	*d[0] - 21*(1 + a)*d[1] + 2*(3 + 4*a)*d[2] + 3*(1 + a)*d[3]) + 22*c[2]*c[2]*
	(3*(34 + 15*a)*d[0] - 33*(1 + a)*d[1] - 2*(64 + 17*a)*d[2] + 7*(1 + a)*d[3])
	 + c[3]*(11*(63*(1 + a) + 2*(103 + 50*a)*c[3])*d[0] - 22*(90 + 36*a + 47*
	(1 + a)*c[3])*d[1] - 2*(495*(1 + a) + 2*(485 + 164*a)*c[3])*d[2] + 2*(22*
	(103 + 50*a) + 531*(1 + a)*c[3])*d[3]) + 11*c[1]*(-3*(7*(1 + a)*(5+2*c[2])
	 + 12*(5 + 2*a)*c[3])*d[0] + 2*(6*(21 + 14*a + 6*c[2] + 8*a*c[2] + 3*(1+a)
	*c[3])*d[1] - (3*(1 + a)*(7 + 22*c[2]) - 4*(29 + 16*a)*c[3])*d[2] + 2*
	(-45 - 18*a + 58*c[2] + 32*a*c[2] - 47*(1 + a)*c[3])*d[3])))))/55440.;
		bf_mom[7] = (pow(-1 + a,2)*(11*(21*(5*(1 + 2*a)*(-4 + c[0]*c[0]) - 
	20*a*c[0]*c[1] + 4*(1 + 4*a)*c[1]*c[1]) - 84*((3 + 2*a)*c[0] + 2*a*c[1])
	*c[2] + 12*(19 + 30*a)*c[2]*c[2] - 36*(6*c[1] + a*(-7*c[0] + 8*c[1] + 
	10*c[2]))*c[3] + 4*(53 + 100*a)*c[3]*c[3])*d[0] + 2*(-693*c[0]*c[0]*d[2] + 
	a*(-11*(105*c[0]*c[0] + c[0]*(-336*c[1] + 84*c[2] + 144*c[3]) + 4*(-105 + 
	63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3] - 6*c[1]*(8*c[2]
	 + 3*c[3])))*d[1] - 2*(231*c[0]*c[0] - 44*(21 + 12*c[1]*c[1] - 33*c[1]*c[2]
	 - 17*c[2]*c[2]) - 44*(32*c[1] + 7*c[2])*c[3] + 656*c[3]*c[3] + 66*c[0]*
	(7*c[1] + 15*(-2*c[2] + c[3])))*d[2] + (693*c[0]*c[0] - 396*c[0]*(4*c[1] +
	 5*c[2]) + 44*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2])) + 8*(550*c[0] - 517*
	c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3])*d[3]) + 44*c[0]*(57*c[2]*d[2] + 
	3*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-27*d[1] + 53*d[3])) - 4*(-693*d[2] + 
	33*c[1]*c[1]*d[2] + 517*c[2]*c[2]*d[2] + 321*c[3]*c[3]*d[2] + c[2]*c[3]*
	(-286*d[1] + 642*d[3]) - 22*c[1]*(13*c[3]*d[2] + c[2]*
	(-3*d[1] + 13*d[3]))))))/110880.;
		bf_mom[8] = -(pow(-1 + a,2)*(11*(21*(5*(2 + a)*(-4 + c[0]*c[0]) - 
	10*(1 + a)*c[0]*c[1] + 4*(3 + 2*a)*c[1]*c[1]) - 84*((4 + a)*c[0] + (1 + a)
	*c[1])*c[2] + 12*(34 + 15*a)*c[2]*c[2] + 18*(7*(1 + a)*c[0] - 4*(5 + 2*a)
	*c[1] - 10*(1 + a)*c[2])*c[3] + 4*(103 + 50*a)*c[3]*c[3])*d[0] - 11*(3*(35*
	(1 + a)*(-4 + c[0]*c[0]) - 56*(3 + 2*a)*c[0]*c[1] + 84*(1 + a)*c[1]*c[1] + 
	4*(7*(1 + a)*c[0] - 4*(3 + 4*a)*c[1])*c[2] + 44*(1 + a)*c[2]*c[2]) + 8*(9*
	(5 + 2*a)*c[0] - 9*(1 + a)*c[1] - 2*(29 + 16*a)*c[2])* c[3] + 188*(1 + a)*
	c[3]*c[3])*d[1] - 2*(33*(7*(4 + a)*(-4 + c[0]*c[0]) + 14*(1 + a)*c[0]*c[1]
	 - 4*(3 + 4*a)*c[1]*c[1]) - 132*((34 + 15*a)*c[0] - 11*(1 + a)*c[1])*c[2] 
	+ 44*(64 + 17*a)*c[2]*c[2] + 22*(45*(1 + a)*c[0] - 4*(29 + 16*a)*c[1] - 
	14*(1 + a)*c[2])* c[3] + 4*(485 + 164*a)*c[3]*c[3])*d[2] + (99*(7*(1 + a)*
	(-4 + c[0]*c[0]) - 8*(5 + 2*a)*c[0]*c[1] + 4*(1 + a)*c[1]*c[1]) - 44*(45*
	(1 + a)*c[0] - 4*(29 + 16*a)*c[1])*c[2] + 308*(1 + a)*c[2]*c[2] + 8*(11*
	(103 + 50*a)*c[0] - 517*(1 + a)*c[1] - 2*(485 + 164*a)*c[2])*c[3] + 2124*
	(1 + a)*c[3]*c[3])*d[3]))/ 27720.;
		break;
	case 4:
	case 6:
		bf_mom[0] = ((1 + b)*(22*c[0]*(21*(-25 - 20*c[1] - 14*c[2] + 12*c[3]
	 + b*(25 - 10*c[1] + 22*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3])+6*c[3]))
	*d[0] + 6*(7*(2*(5 + 9*c[1] - 2*c[2]) + b*(5 - 5*b - 14*c[1] + 8*b*c[1] + 
	2*(-1 + b)*c[2])) - 6*(11 + b*(-13 + 4*b))*c[3])*d[1] + 6*(49 - 28*c[1] + 
	158*c[2] - 60*c[3] + b*(-77 + 14*b - 14*c[1] + 14*b*c[1] - 174*c[2] + 
	60*b*c[2] + 30*(-1 + b)*c[3]))*d[2] - 2*(9*(14 + 22*c[1] + 20*c[2] + b*
	(7 - 7*b - 26*c[1] + 8*b*c[1] - 10*(-1 + b)*c[2])) - 2*(253 + b*(-259 + 
	100*b))*c[3])*d[3]) + 231*c[0]*c[0]*(5*(5 + b*(-5 + 2*b))*d[0] + 2*(5*(-2 
	+ b)*(1 + b)*d[1] - (7 + b*(-11 + 2*b))*d[2] - 3*(-2 + b)*(1 + b)*d[3])) + 
	4*(c[2]*(33*(7*(7 + b*(-11 + 2*b)) + 30*(-2 + b)*(1 + b)*c[3])* d[0] + 22*
	(-21*(-2 + b)*(1 + b) + 2*(77 + b*(-71 + 32*b))*c[3])*d[1] - 22*(3*(79 - 
	87*b + 30*b*b) + 14*(-2 + b)*(1 + b)*c[3])* d[2] - 2*(495*(-2 + b)*(1 + b)
	 + 2*(977 + b*(-1291 + 328*b))*c[3])*d[3]) + 11*c[2]*c[2]*(3*(79 - 87*b + 
	30*b*b)*d[0] + 2*(33*(-2 + b)*(1 + b)*d[1] - (115 + b*(-175 + 34*b))*d[2]
	 - 7*(-2 + b)*(1 + b)*d[3])) + 33*c[1]*c[1]*(7*(9 + b*(-7 + 4*b))*d[0] + 
	2*(21*(-2 + b)*(1 + b)*d[1] + (15 + b*(-5 + 8*b))*d[2] - 3*(-2 + b)*(1 + b)
	*d[3])) + 11*c[1]*(3*(7*(-2 + b)*(1 + b)*(-5 + 2*c[2]) - 6*(11 + b*(-13 + 
	4*b))*c[3])*d[0] + 2*(3*(-63 + 30*c[2] + 12*c[3] + b*(49 - 28*b - 10*c[2] 
	+ 16*b*c[2] - 6*(-1 + b)*c[3]))* d[1] + (3*(-2 + b)*(1 + b)*(-7 + 22*c[2]) 
	+ 2*(77 + b*(-71 + 32*b))*c[3])*d[2] + (99 - 117*b + 36*b*b + 154*c[2] - 
	142*b*c[2] + 64*b*b*c[2] + 94*(-2 + b)*(1 + b)*c[3])*d[3])) + c[3]*(11*(63*
	(-2 + b)*(1 + b) + (253 + b*(-259 + 100*b))*c[3])* d[0] + 2*(11*(9*(11 + 
	b*(-13 + 4*b)) + 47*(-2 + b)*(1 + b)*c[3])*d[1] - (495*(-2 + b)*(1 + b) + 
	(977 + b*(-1291 + 328*b))*c[3])* d[2] - (11*(253 + b*(-259 + 100*b)) + 
	531*(-2 + b)*(1 + b)*c[3])*d[3])))))/221760.;
		bf_mom[1] = (pow(1 + b,2)*(231*c[0]*c[0]* (5*(-1 + 2*b)*d[0] + 
	6*d[2] + 2*b*(5*d[1] - 2*d[2] - 3*d[3])) + 22*c[0]*(21*(5 + 6*c[2] + 2*b*
	(-5 + 5*c[1] - 2*c[2] - 3*c[3]))* d[0] + 6*(7*(-2*c[1] + b*(-5 + 8*c[1] + 
	2*c[2])) - 6*(-3 + 4*b)*c[3])*d[1] + 6*(-21 - 38*c[2] + 2*b*(7 + 7*c[1] + 
	30*c[2] + 15*c[3]))*d[2] - 2*(-54*c[1] + b*(-63 + 72*c[1] - 90*c[2] - 
	200*c[3]) + 106*c[3])* d[3]) + 4*(-693*c[2]*d[0] - 627*c[2]*c[2]*d[0] - 
	583*c[3]*c[3]*d[0] - 594*c[3]*d[1] - 572*c[2]*c[3]*d[1] - 33*c[1]*c[1]*
	(7*d[0] - 2*d[2]) + 1254*c[2]*d[2] + 1034*c[2]*c[2]*d[2] + 642*c[3]*c[3]
	*d[2] + 2*(583 + 642*c[2])*c[3]*d[3] + 22*c[1]*(3*(7 + 2*c[2])*d[1] + c[3]*
	(27*d[0] - 26*d[2]) - (27 + 26*c[2])*d[3]) + b*(22*c[2]*(3*(7 + 15*c[3])*
	d[0] + (-21 + 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) + 22*c[2]*c[2]*(45*d[0]
	 + 33*d[1] - 34*d[2] - 7*d[3]) + 66*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2] 
	- 3*d[3]) - 2*c[2]*(495 + 656*c[3])*d[3] + 11*c[1]*(3*(-35 + 14*c[2] - 
	24*c[3])*d[0] + 2*(6*(-14 + 8*c[2] - 3*c[3])*d[1] + (-21 + 66*c[2] + 64*
	c[3])*d[2] + 2*(18 + 32*c[2] + 47*c[3])*d[3])) + c[3]*(11*(63 + 100*c[3])
	*d[0] + 2*c[3]*(517*d[1] - 328*d[2] - 531*d[3]) + 22*(36*d[1] - 5*(9*d[2]
	 + 20*d[3])))))))/221760.;
		bf_mom[2] = (pow(1 + b,2)*(231*c[0]*c[0]* (5*(-1 + 2*b)*d[0] + 
	6*d[2] + 2*b*(5*d[1] - 2*d[2] - 3*d[3])) + 22*c[0]*(21*(-5 + 6*c[2] + 2*b*
	(5 + 5*c[1] - 2*c[2] - 3*c[3]))* d[0] + 6*(7*(-2*c[1] + b*(5 + 8*c[1] + 
	2*c[2])) - 6*(-3 + 4*b)*c[3])*d[1] + 6*(21 - 38*c[2] + 2*b*(-7 + 7*c[1] + 
	30*c[2] + 15*c[3]))*d[2] - 2*(-54*c[1] + b*(63 + 72*c[1] - 90*c[2] - 200*
	c[3]) + 106*c[3])* d[3]) + 4*(693*c[2]*d[0] - 627*c[2]*c[2]*d[0] - 583*
	c[3]*c[3]*d[0] + 594*c[3]*d[1] - 572*c[2]*c[3]*d[1] - 33*c[1]*c[1]*(7*d[0]
	 - 2*d[2]) - 1254*c[2]*d[2] + 1034*c[2]*c[2]*d[2] + 642*c[3]*c[3]*d[2] + 
	2*(-583 + 642*c[2])*c[3]*d[3] + 22*c[1]*(3*(-7 + 2*c[2])*d[1] + c[3]*
	(27*d[0] - 26*d[2]) + (27 - 26*c[2])*d[3]) + b*(22*c[2]*(3*(-7 + 15*c[3])
	*d[0] + (21 + 64*c[3])*d[1] + 2*(45 - 7*c[3])*d[2]) + 22*c[2]*c[2]*(45*d[0]
	 + 33*d[1] - 34*d[2] - 7*d[3]) + 66*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2]
	 - 3*d[3]) - 2*c[2]*(-495 + 656*c[3])*d[3] + c[3]*(11*(-63 + 100*c[3])*d[0]
	 - 792*d[1] + 990*d[2] + 2*c[3]*(517*d[1] - 328*d[2] - 531*d[3])+2200*d[3])
	 + 11*c[1]*(3*(35 + 14*c[2] - 24*c[3])*d[0] + 2*(6*(14 + 8*c[2] - 3*c[3])
	*d[1] + (21 + 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 32*c[2] + 
	47*c[3])*d[3]))))))/221760.;
		bf_mom[3] = ((1 + b)*(22*c[0]*(21*(25 - 20*c[1] - 14*c[2] + 12*c[3]
	 + b*(-25 - 10*c[1] + 22*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3])+6*c[3]))
	*d[0] + 6*(7*(2*(-5 + 9*c[1] - 2*c[2]) + b*(-5 + 5*b - 14*c[1] + 8*b*c[1]
	 + 2*(-1 + b)*c[2])) - 6*(11 + b*(-13 + 4*b))*c[3])*d[1] + 6*(-49 - 28*c[1]
	 + 158*c[2] - 60*c[3] + b*(77 - 14*b - 14*c[1] + 14*b*c[1] - 174*c[2] + 
	60*b*c[2] + 30*(-1 + b)*c[3]))*d[2] - 2*(9*(-14 + 22*c[1] + 20*c[2] + b*
	(-7 + 7*b - 26*c[1] + 8*b*c[1] - 10*(-1 + b)*c[2])) - 2*(253 + b*(-259 + 
	100*b))*c[3])*d[3]) + 231*c[0]*c[0]*(5*(5 + b*(-5 + 2*b))*d[0] + 2*(5*
	(-2 + b)*(1 + b)*d[1] - (7 + b*(-11 + 2*b))*d[2] - 3*(-2 + b)*(1 + b)*d[3]))
	 + 4*(c[2]*(33*(-7*(7 + b*(-11 + 2*b)) + 30*(-2 + b)*(1 + b)*c[3])* d[0] 
	+ 22*(21*(-2 + b)*(1 + b) + 2*(77 + b*(-71 + 32*b))*c[3])*d[1] - 22*(-3*
	(79 - 87*b + 30*b*b) + 14*(-2 + b)*(1 + b)*c[3])*d[2] - 2*(-495*(-2 + b)*
	(1 + b) + 2*(977 + b*(-1291 + 328*b))*c[3])* d[3]) + 11*c[2]*c[2]* (3*
	(79 - 87*b + 30*b*b)*d[0] + 2*(33*(-2 + b)*(1 + b)*d[1] - (115 + b*(-175 
	+ 34*b))*d[2] - 7*(-2 + b)*(1 + b)*d[3])) + 33*c[1]*c[1]*(7*(9 + b*(-7 + 
	4*b))*d[0] + 2*(21*(-2 + b)*(1 + b)*d[1] + (15 + b*(-5 + 8*b))*d[2] - 3*
	(-2 + b)*(1 + b)*d[3])) + c[3]*(11*(-63*(-2 + b)*(1 + b) + (253 + b*(-259 
	+ 100*b))*c[3])* d[0] + 2*(11*(-9*(11 + b*(-13 + 4*b)) + 47*(-2 + b)*(1+b)
	*c[3])*d[1] - (-495*(-2 + b)*(1 + b) + (977 + b*(-1291 + 328*b))*c[3])*
	 d[2] + (11*(253 + b*(-259 + 100*b)) - 531*(-2 + b)*(1 + b)*c[3])*d[3])) + 
	11*c[1]*(3*(7*(-2 + b)*(1 + b)*(5 + 2*c[2]) - 6*(11 + b*(-13 + 4*b))*c[3])
	*d[0] + 2*(3*(63 + 30*c[2] + 12*c[3] + b*(-49 + 28*b - 10*c[2] + 16*b*c[2]
	 - 6*(-1 + b)*c[3]))* d[1] + (3*(-2 + b)*(1 + b)*(7 + 22*c[2]) + 2*(77 + 
	b*(-71 + 32*b))*c[3])*d[2] + (-99 + 117*b - 36*b*b + 154*c[2] - 142*b*c[2]
	 + 64*b*b*c[2] + 94*(-2 + b)*(1 + b)*c[3])*d[3])))))/ 221760.;
		bf_mom[4] = -(pow(1 + b,2)*(231*c[0]*c[0]* (5*(-2 + b)*d[0] + 5*
	(-1 + b)*d[1] - 2*(-4 + b)*d[2] - 3*(-1 + b)*d[3]) + 22*c[0]* (21*(10 - 
	5*c[1] + 8*c[2] + b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 3*c[3])*d[0] + 3* (7*
	(5 - 5*b - 12*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2]) - 12*(-5 + 2*b)*c[3])*d[1]
	 + 6*(-28 + 7*b - 7*c[1] + 7*b*c[1] - 68*c[2] + 30*b*c[2] + 15*(-1 + b)*
	c[3])*d[2] - (9*(7 - 7*b - 20*c[1] + 8*b*c[1] - 10*(-1 + b)*c[2]) - 4*
	(-103 + 50*b)*c[3])*d[3]) + 2*(22*c[2]*(3*(7*(-4 + b) + 15*(-1 + b)*c[3])
	*d[0] + (21 - 21*b - 116*c[3] + 64*b*c[3])*d[1] + 2*(102 - 45*b + 7*c[3] 
	- 7*b*c[3])*d[2]) - 2*c[2]*(495*(-1 + b) + 4*(-485 + 164*b)*c[3])*d[3] + 
	22*c[2]*c[2]*(3*(-34 + 15*b)*d[0] + 33*(-1 + b)*d[1] - 2*(-64 + 17*b)*d[2] 
	- 7*(-1 + b)*d[3]) + 66*c[1]*c[1]*(7*(-3 + 2*b)*d[0] + 21*(-1 + b)*d[1] + 
	2*(-3 + 4*b)*d[2] - 3*(-1 + b)*d[3]) + c[3]*(11*(63*(-1 + b) + 2*(-103 + 
	50*b)*c[3])*d[0] + 2*(11*(-90 + 36*b + 47*(-1 + b)*c[3])*d[1] + (495 - 
	495*b + 970*c[3] - 328*b*c[3])*d[2] + (2266 - 1100*b - 531*(-1 + b)*c[3])
	*d[3])) + 11*c[1]*(3*(7*(-1 + b)*(-5 + 2*c[2]) - 12*(-5 + 2*b)*c[3])*d[0] 
	+ 2*(6*(21 - 14*b - 6*c[2] + 8*b*c[2] - 3*(-1 + b)*c[3])*d[1] + (3*(-1 + b)
	*(-7 + 22*c[2]) + 4*(-29 + 16*b)*c[3])*d[2] + 2*(-45 + 18*b - 58*c[2] + 
	32*b*c[2] + 47*(-1 + b)*c[3])*d[3])))))/55440.;
		bf_mom[5] = -(pow(1 + b,2)*(11*(21*(5*(-1 + 2*b)*(-4 + c[0]*c[0]) 
	+ 20*b*c[0]*c[1] + 4*(-1 + 4*b)*c[1]*c[1]) - 84*((-3 + 2*b)*c[0] - 2*b*c[1])
	*c[2] + 12*(-19 + 30*b)*c[2]*c[2] - 36*(-6*c[1] + b*(7*c[0] + 8*c[1] - 
	10*c[2]))*c[3] + 4*(-53 + 100*b)*c[3]*c[3])*d[0] + 2*(693*c[0]*c[0]*d[2] + 
	b*(11*(3*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] + 64*c[1]*c[2] + 44*c[2]*c[2]
	 + 28*c[0]*(4*c[1] + c[2])) - 8*(18*c[0] + 9*c[1] - 32*c[2])*c[3] + 188*
	c[3]*c[3])* d[1] - 2*(231*c[0]*c[0] - 44*(21 + 12*c[1]*c[1] + 33*c[1]*c[2]
	 - 17*c[2]*c[2]) - 44*(32*c[1] - 7*c[2])*c[3] + 656*c[3]*c[3] - 66*c[0]*
	(7*c[1] + 30*c[2] + 15*c[3]))* d[2] - (99*(-28 + (c[0] + 2*c[1])*(7*c[0] + 
	2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2] + 308*c[2]*c[2] - 8*(550*c[0] + 
	517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3])*d[3]) - 44*c[0]*(57*c[2]*d[2]
	 + 3*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-27*d[1] + 53*d[3])) + 4*(-693*d[2] + 
	33*c[1]*c[1]*d[2] + 517*c[2]*c[2]*d[2] + 321*c[3]*c[3]*d[2] + c[2]*c[3]*
	(-286*d[1] + 642*d[3]) - 22*c[1]*(13*c[3]*d[2] + c[2]*
	(-3*d[1] + 13*d[3]))))))/110880.;
		bf_mom[6] = -(pow(1 + b,2)*(231*c[0]*c[0]* (5*(-2 + b)*d[0] + 5*
	(-1 + b)*d[1] - 2*(-4 + b)*d[2] - 3*(-1 + b)*d[3]) + 22*c[0]* (21*(-10 - 
	5*c[1] + 8*c[2] + b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 3*c[3])*d[0] + 3* (7*
	(-5 + 5*b - 12*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2]) - 12*(-5 + 2*b)*c[3])*d[1]
	 + 6*(28 - 7*b - 7*c[1] + 7*b*c[1] - 68*c[2] + 30*b*c[2] + 15*(-1 + b)*
	c[3])*d[2] - (9*(-7 + 7*b - 20*c[1] + 8*b*c[1] - 10*(-1 + b)*c[2]) - 4*
	(-103 + 50*b)*c[3])*d[3]) + 2*(22*c[2]*(3*(-7*(-4 + b) + 15*(-1 + b)*c[3])
	*d[0] + (21*(-1 + b) + 4*(-29 + 16*b)*c[3])*d[1] + 2*(-102 + 45*b + 7*c[3]
	 - 7*b*c[3])*d[2]) - 2*c[2]*(-495*(-1 + b) + 4*(-485 + 164*b)*c[3])*d[3] 
	+ 22*c[2]*c[2]*(3*(-34 + 15*b)*d[0] + 33*(-1 + b)*d[1] - 2*(-64 + 17*b)*d[2]
	 - 7*(-1 + b)*d[3]) + 66*c[1]*c[1]*(7*(-3 + 2*b)*d[0] + 21*(-1 + b)*d[1] + 
	2*(-3 + 4*b)*d[2] - 3*(-1 + b)*d[3]) + c[3]*(11*(63 - 206*c[3] + b*(-63 
	+ 100*c[3]))*d[0] + 2*(11*(90 - 47*c[3] + b*(-36 + 47*c[3]))*d[1] + (-495 
	+ 495*b + 970*c[3] - 328*b*c[3])*d[2] + (22*(-103 + 50*b) - 531*(-1 + b)
	*c[3])*d[3])) + 11*c[1]*(3*(7*(-1 + b)*(5 + 2*c[2]) - 12*(-5 + 2*b)*c[3])
	*d[0] + 2*(6*(-21 + 14*b - 6*c[2] + 8*b*c[2] - 3*(-1 + b)*c[3])*d[1] + 
	(3*(-1 + b)*(7 + 22*c[2]) + 4*(-29 + 16*b)*c[3])*d[2] + 2*(45 - 18*b - 
	58*c[2] + 32*b*c[2] + 47*(-1 + b)*c[3])*d[3])))))/55440.;
		bf_mom[7] = -((1 + b)*(11*(21*(5*(5 + b*(-5 + 2*b))*(-4 + c[0]*c[0])
	 + 20*(-2 + b)*(1 + b)*c[0]*c[1] + 4*(9 + b*(-7 + 4*b))*c[1]*c[1]) - 84*
	((7 + b*(-11 + 2*b))*c[0] - 2*(-2 + b)*(1 + b)*c[1])*c[2] + 12*(79 - 87*b 
	+ 30*b*b)*c[2]*c[2] - 36*(7*(-2 + b)*(1 + b)*c[0] + 2*(11 + b*(-13 + 4*b))
	*c[1] - 10*(-2 + b)*(1 + b)*c[2])*c[3] + 4*(253 + b*(-259 + 100*b))*c[3]
	*c[3])*d[0] + 2*(11*(3*(35*(-2 + b)*(1 + b)*(-4 + c[0]*c[0]) + 28*(9 + b*
	(-7 + 4*b))*c[0]*c[1] + 84*(-2 + b)*(1 + b)*c[1]*c[1] + 4*(7*(-2 + b)*
	(1 + b)*c[0] + 2*(15 + b*(-5 + 8*b))*c[1])*c[2] + 44*(-2 + b)*(1 + b)*
	c[2]*c[2]) - 4*(9*(11 + b*(-13 + 4*b))*c[0] + 18*(-2 + b)*(1 + b)*c[1] - 
	2*(77 + b*(-71 + 32*b))*c[2])* c[3] + 188*(-2 + b)*(1 + b)*c[3]*c[3])*d[1]
	 + (33*(-7*(7 + b*(-11 + 2*b))*(-4 + c[0]*c[0]) + 28*(-2 + b)*(1 + b)*c[0]
	*c[1] + 4*(15 + b*(-5 + 8*b))*c[1]*c[1]) + 132*((79 - 87*b + 30*b*b)*c[0] 
	+ 22*(-2 + b)*(1 + b)*c[1])*c[2] - 44*(115 + b*(-175 + 34*b))*c[2]*c[2] + 
	44*(45*(-2 + b)*(1 + b)*c[0] + 2*(77 + b*(-71 + 32*b))*c[1] - 14*(-2 + b)*
	(1 + b)*c[2])* c[3] - 4*(977 + b*(-1291 + 328*b))*c[3]*c[3])*d[2] + (99*
	(-7*(-2 + b)*(1 + b)*(-4 + c[0]*c[0]) - 4*(11 + b*(-13 + 4*b))*c[0]*c[1] 
	- 4*(-2 + b)*(1 + b)*c[1]*c[1]) + 44*(45*(-2 + b)*(1 + b)*c[0] + 2*(77 + 
	b*(-71 + 32*b))*c[1])* c[2] - 308*(-2 + b)*(1 + b)*c[2]*c[2] + 4*(11*
	(253 + b*(-259 + 100*b))*c[0] + 1034*(-2 + b)*(1 + b)*c[1] - 2*(977 + b*
	(-1291 + 328*b))*c[2])*c[3] - 2124*(-2 + b)*
	(1 + b)*c[3]*c[3])*d[3])))/110880.;
		bf_mom[8] = (pow(1 + b,2)*(11*(21*(5*(-2 + b)*(-4 + c[0]*c[0]) + 
	10*(-1 + b)*c[0]*c[1] + 4*(-3 + 2*b)*c[1]*c[1]) - 84*((-4 + b)*c[0] + c[1]
	 - b*c[1])*c[2] + 12*(-34 + 15*b)*c[2]*c[2] - 18*(7*(-1 + b)*c[0] + 4*
	(-5 + 2*b)*c[1] - 10*(-1 + b)*c[2])* c[3] + 4*(-103 + 50*b)*c[3]*c[3])*d[0]
	 + 11*(3*(35*(-1 + b)*(-4 + c[0]*c[0]) + 56*(-3 + 2*b)*c[0]*c[1] + 84*(-1 
	+ b)*c[1]*c[1] + 4*(7*(-1 + b)*c[0] + 4*(-3 + 4*b)*c[1])*c[2] + 44*(-1 + b)
	*c[2]*c[2]) - 8*(9*(-5 + 2*b)*c[0] + 9*(-1 + b)*c[1] + 58*c[2] - 32*b*c[2])
	* c[3] + 188*(-1 + b)*c[3]*c[3])*d[1] - 2*(33*(7*(-4 + b)*(-4 + c[0]*c[0]) 
	- 14*(-1 + b)*c[0]*c[1] - 4*(-3 + 4*b)*c[1]*c[1]) - 132*((-34 + 15*b)*c[0] 
	+ 11*(-1 + b)*c[1])*c[2] + 44*(-64 + 17*b)*c[2]*c[2] - 22*(45*(-1 + b)*c[0]
	 + 4*(-29 + 16*b)*c[1] - 14*(-1 + b)*c[2])* c[3] + 4*(-485 + 164*b)*c[3]
	*c[3])*d[2] - (99*(7*(-1 + b)*(-4 + c[0]*c[0]) + 8*(-5 + 2*b)*c[0]*c[1] + 
	4*(-1 + b)*c[1]*c[1]) - 44*(45*(-1 + b)*c[0] + 4*(-29 + 16*b)*c[1])*c[2] + 
	308*(-1 + b)*c[2]*c[2] - 8*(11*(-103 + 50*b)*c[0] + 517*(-1 + b)*c[1] + 
	970*c[2] - 328*b*c[2])*c[3] + 2124*(-1 + b)*c[3]*c[3])*d[3]))/27720.;
		break;
	case 22:
		bf_mom[0] = (22*c[0]*(-21*(5 + 10*c[1] - 2*c[2] - 6*c[3])*d[0] + 
	6*(35 + 42*c[1] - 14*c[2] - 6*c[3])*d[1] - 6*(7 + 14*c[1] - 22*c[2] + 30*
	c[3])*d[2] - 2*(63 + 18*c[1] + 90*c[2] - 94*c[3])*d[3]) + 231*c[0]*c[0]*
	(5*d[0] + 2*(-5*d[1] + d[2] + 3*d[3])) + 4*(c[3]*(11*(-63 + 47*c[3])*d[0] 
	- 2*c[3]*(517*d[1] + 7*d[2] - 531*d[3]) + 22*(9*d[1] + 45*d[2] - 47*d[3])) 
	+ 11*c[2]*c[2]*(33*d[0] - 66*d[1] + 26*d[2] + 14*d[3]) + 99*c[1]*c[1]*
	(7*d[0] + 2*(-7*d[1] + 3*d[2] + d[3])) + c[2]*(-33*(7 + 30*c[3])*d[0] + 
	4*c[3]*(209*d[1] + 77*d[2] - 7*d[3]) + 66*(7*d[1] - 11*d[2] + 15*d[3])) + 
	11*c[1]*(-3*(-35 + 14*c[2] + 6*c[3])*d[0] + 2*(9*(-7 + 6*c[2] + 2*c[3])*d[1]
	 + (21 - 66*c[2] + 38*c[3])*d[2] + (9 + 38*c[2] - 94*c[3])*d[3]))))/55440.;
		bf_mom[1] = (22*c[0]*(3*(7*(5 - 10*c[1] + 2*c[2] + 6*c[3])*d[0] + 
	2*(-35 + 42*c[1] - 14*c[2] - 6*c[3])*d[1] + 2*(7 - 14*c[1] + 22*c[2] - 
	30*c[3])*d[2]) - 2*(-63 + 18*c[1] + 90*c[2] - 94*c[3])*d[3]) + 231*c[0]*
	c[0]*(5*d[0] + 2*(-5*d[1] + d[2] + 3*d[3])) + 4*(11*c[2]*c[2]*(33*d[0] - 
	66*d[1] + 26*d[2] + 14*d[3]) + c[2]*(-33*(-7 + 30*c[3])*d[0] + 22*(-21 + 
	38*c[3])*d[1] + 22*(33 + 14*c[3])*d[2] - 2*(495 + 14*c[3])*d[3]) + c[3]*
	(11*(63 + 47*c[3])*d[0] - 22*(9 + 47*c[3])*d[1] - 2*(495 + 7*c[3])*d[2] + 
	2*(517 + 531*c[3])*d[3]) + 99*c[1]*c[1]*(7*d[0] + 2*(-7*d[1] + 3*d[2] + 
	d[3])) - 11*c[1]*(3*(35 + 14*c[2] + 6*c[3])*d[0] - 2*(9*(7 + 6*c[2] + 
	2*c[3])*d[1] + (-21 - 66*c[2] + 38*c[3])*d[2] + (-9 + 38*c[2] - 94*c[3])
	*d[3]))))/55440.;
		bf_mom[2] = (231*c[0]*c[0]*(5*d[0] + 2*(5*d[1] + d[2] - 3*d[3])) + 
	22*c[0]*(21*(5 + 10*c[1] + 2*c[2] - 6*c[3])*d[0] + 6*(35 + 42*c[1] + 14*c[2]
	 - 6*c[3])*d[1] + 6*(7 + 14*c[1] + 22*c[2] + 30*c[3])*d[2] - 2*(63 + 18*c[1]
	 - 90*c[2] - 94*c[3])*d[3]) + 4*(11*c[2]*c[2]*(33*d[0] + 66*d[1] + 26*d[2] 
	- 14*d[3]) + 99*c[1]*c[1]*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) + c[2]*(33*
	(7 + 30*c[3])*d[0] + 22*(21 + 38*c[3])*d[1] + 726*d[2] + 990*d[3] - 28*c[3]
	*(11*d[2] + d[3])) + c[3]*(11*(-63 + 47*c[3])*d[0] + 2*c[3]*(517*d[1] - 
	7*d[2] - 531*d[3]) + 22*(-9*d[1] + 45*d[2] + 47*d[3])) + 11*c[1]*(3*(35 + 
	14*c[2] - 6*c[3])*d[0] + 2*(9*(7 + 6*c[2] - 2*c[3])*d[1] + (21 + 66*c[2] + 
	38*c[3])*d[2] + (-9 + 38*c[2] + 94*c[3])*d[3]))))/55440.;
		bf_mom[3] = (231*c[0]*c[0]*(5*d[0] + 2*(5*d[1] + d[2] - 3*d[3])) + 
	22*c[0]*(21*(-5 + 10*c[1] + 2*c[2] - 6*c[3])*d[0] + 6*(-35 + 42*c[1] + 14*
	c[2] - 6*c[3])*d[1] + 6*(-7 + 14*c[1] + 22*c[2] + 30*c[3])*d[2] - 2*(-63 + 
	18*c[1] - 90*c[2] - 94*c[3])*d[3]) + 4*(11*c[2]*c[2]*(33*d[0] + 66*d[1] + 
	26*d[2] - 14*d[3]) + 99*c[1]*c[1]*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) + 
	11*c[1]*(3*(-35 + 14*c[2] - 6*c[3])*d[0] + 2*(9*(-7 + 6*c[2] - 2*c[3])*d[1]
	 + (-21 + 66*c[2] + 38*c[3])*d[2] + (9 + 38*c[2] + 94*c[3])*d[3])) + c[3]*
	(11*(63 + 47*c[3])*d[0] + 2*(11*(9 + 47*c[3])*d[1] - (495 + 7*c[3])*d[2] - 
	(517 + 531*c[3])*d[3])) + c[2]*(33*(-7 + 30*c[3])*d[0] - 66*(7*d[1] + 
	11*d[2] + 15*d[3]) + 4*c[3]*(209*d[1] - 7*(11*d[2] + d[3])))))/55440.;
		bf_mom[4] = (-11*(21*(-20 + 5*c[0]*c[0] - 20*c[0]*c[1] + 12*c[1]*
	c[1]) + 84*(c[0] - 2*c[1])*c[2] + 132*c[2]*c[2] + 36*(7*c[0] - 2*(c[1] + 
	5*c[2]))*c[3] + 188*c[3]*c[3])*d[0] - 1848*(5*d[1] - d[2] - 3*d[3]) + 2*
	(231*c[0]*c[0]*(5*d[1] - d[2] - 3*d[3]) + 44*c[0]*(c[3]*(9*d[1] + 45*d[2] 
	- 47*d[3]) + c[1]*(-63*d[1] + 21*d[2] + 9*d[3]) + 3*c[2]*(7*d[1] - 11*d[2] 
	+ 15*d[3])) + 4*(c[3]*c[3]*(517*d[1] + 7*d[2] - 531*d[3]) + 11*c[2]*c[2]*
	(33*d[1] - 13*d[2] - 7*d[3]) - 2*c[2]*c[3]*(209*d[1] + 77*d[2] - 7*d[3]) + 
	99*c[1]*c[1]*(7*d[1] - 3*d[2] - d[3]) - 22*c[1]*(c[3]*(9*d[1] + 19*d[2] - 
	47*d[3]) + c[2]*(27*d[1] - 33*d[2] + 19*d[3])))))/27720.;
		bf_mom[5] = (231*c[0]*c[0]*(5*d[0] - 6*d[2]) + 4*(11*c[2]*c[2]*
	(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] 
	- 594*d[1] - 642*c[3]*d[2] + 1166*d[3]) - 22*c[1]*(3*(-7 + 2*c[2])*d[1] + 
	c[3]*(27*d[0] - 26*d[2]) + (27 - 26*c[2])*d[3]) + c[2]*(-693*d[0] + 572*
	c[3]*d[1] + 1254*d[2] - 1284*c[3]*d[3])) + 22*c[0]*(-21*(-5 + 6*c[2])*d[0] 
	+ 2*(3*(-21 + 38*c[2])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1] 
	+ 53*d[3]))))/13860.;
		bf_mom[6] = (-11*(105*c[0]*c[0] + 84*c[0]*(5*c[1] + c[2] - 3*c[3]) 
	+ 4*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(7*c[2] - 3*c[3]) + 90*
	c[2]*c[3] + 47*c[3]*c[3]))*d[0] - 22*(3*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1]
	 + 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(3*c[1] + c[2])) - 4*(9*c[0] + 
	18*c[1] - 38*c[2])*c[3] + 188*c[3]*c[3])*d[1] - 2*(33*(-28 + 7*c[0]*c[0] + 
	28*c[0]*c[1] + 36*c[1]*c[1]) + 1452*(c[0] + 2*c[1])*c[2] + 572*c[2]*c[2] + 
	44*(45*c[0] + 38*c[1] - 14*c[2])*c[3] - 28*c[3]*c[3])*d[2] + 2*(693*c[0]*
	c[0] + 396*c[0]*(c[1] - 5*c[2]) + 44*(9*c[1]*c[1] - 38*c[1]*c[2] + 7*(-9 
	+ c[2]*c[2])) - 4*(517*(c[0] + 2*c[1]) - 14*c[2])*c[3] + 2124*c[3]*c[3])
	*d[3])/ 27720.;
		bf_mom[7] = (231*c[0]*c[0]*(5*d[0] - 6*d[2]) + 4*(11*c[2]*c[2]*
	(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] 
	+ 594*d[1] - 642*c[3]*d[2] - 1166*d[3]) - 22*c[1]*(3*(7 + 2*c[2])*d[1] + 
	c[3]*(27*d[0] - 26*d[2]) - (27 + 26*c[2])*d[3]) + c[2]*(693*d[0] + 572*
	c[3]*d[1] - 1254*d[2] - 1284*c[3]*d[3])) + 22*c[0]*(-21*(5 + 6*c[2])*d[0] 
	+ 2*(3*(21 + 38*c[2])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1]
	 + 53*d[3]))))/13860.;
		bf_mom[8] = (-11*(84*(-5 + c[1]*c[1]) + 3*(35*c[0]*c[0] - 84*c[0]
	*c[2] + 76*c[2]*c[2]) - 216*c[1]*c[3] + 212*c[3]*c[3])*d[0] + 2*(693*c[0]
	*c[0]*d[2] - 44*c[0]*(57*c[2]*d[2] + 3*c[1]*(7*d[1] - 9*d[3]) + c[3]*
	(-27*d[1] + 53*d[3])) + 4*(-693*d[2] + 33*c[1]*c[1]*d[2] + 517*c[2]*c[2]
	*d[2] + 321*c[3]*c[3]*d[2] + c[2]*c[3]*(-286*d[1] + 642*d[3]) - 22*c[1]*
	(13*c[3]*d[2] + c[2]*(-3*d[1] + 13*d[3])))))/6930.;
		break;
	case 32:
		bf_mom[0] = (22*c[0]*(-21*(5 + 10*c[1] - 2*c[2] - 6*c[3])*d[0] + 6*
	(35 + 42*c[1] - 14*c[2] - 6*c[3])*d[1] - 6*(7 + 14*c[1] - 22*c[2]+30*c[3])
	*d[2] - 2*(63 + 18*c[1] + 90*c[2] - 94*c[3])*d[3]) + 231*c[0]*c[0]*(5*d[0] 
	+ 2*(-5*d[1] + d[2] + 3*d[3])) + 4*(c[3]*(11*(-63 + 47*c[3])*d[0] - 2*c[3]*
	(517*d[1] + 7*d[2] - 531*d[3]) + 22*(9*d[1] + 45*d[2] - 47*d[3])) + 11*c[2]
	*c[2]*(33*d[0] - 66*d[1] + 26*d[2] + 14*d[3]) + 99*c[1]*c[1]*(7*d[0] + 2*
	(-7*d[1] + 3*d[2] + d[3])) + c[2]*(-33*(7 + 30*c[3])*d[0] + 4*c[3]*(209*d[1]
	 + 77*d[2] - 7*d[3]) + 66*(7*d[1] - 11*d[2] + 15*d[3])) + 11*c[1]*(-3*(-35 
	+ 14*c[2] + 6*c[3])*d[0] + 2*(9*(-7 + 6*c[2] + 2*c[3])*d[1] + (21 - 66*c[2]
	 + 38*c[3])*d[2] + (9 + 38*c[2] - 94*c[3])*d[3]))))/55440.;
		bf_mom[1] = (231*c[0]*c[0]*(5*d[0] + 2*(5*d[1] + d[2] - 3*d[3])) + 
	22*c[0]*(21*(-5 + 10*c[1] + 2*c[2] - 6*c[3])*d[0] + 6*(-35 + 42*c[1] + 
	14*c[2] - 6*c[3])*d[1] + 6*(-7 + 14*c[1] + 22*c[2] + 30*c[3])*d[2] - 2*(-63
	 + 18*c[1] - 90*c[2] - 94*c[3])*d[3]) + 4*(11*c[2]*c[2]*(33*d[0] + 66*d[1] 
	+ 26*d[2] - 14*d[3]) + 99*c[1]*c[1]*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) + 
	11*c[1]*(3*(-35 + 14*c[2] - 6*c[3])*d[0] + 2*(9*(-7 + 6*c[2] - 2*c[3])*d[1]
	 + (-21 + 66*c[2] + 38*c[3])*d[2] + (9 + 38*c[2] + 94*c[3])*d[3])) + c[3]*
	(11*(63 + 47*c[3])*d[0] + 2*(11*(9 + 47*c[3])*d[1] - (495 + 7*c[3])*d[2] - 
	(517 + 531*c[3])*d[3])) + c[2]*(33*(-7 + 30*c[3])*d[0] - 66*(7*d[1] + 
	11*d[2] + 15*d[3]) + 4*c[3]*(209*d[1] - 7*(11*d[2] + d[3])))))/55440.;
		bf_mom[2] = (231*c[0]*c[0]*(5*d[0] + 2*(5*d[1] + d[2] - 3*d[3])) + 
	22*c[0]*(21*(5 + 10*c[1] + 2*c[2] - 6*c[3])*d[0] + 6*(35 + 42*c[1] + 14*c[2]
	 - 6*c[3])*d[1] + 6*(7 + 14*c[1] + 22*c[2] + 30*c[3])*d[2] - 2*(63 +18*c[1] 
	- 90*c[2] - 94*c[3])*d[3]) + 4*(11*c[2]*c[2]*(33*d[0] + 66*d[1] + 26*d[2] 
	- 14*d[3]) + 99*c[1]*c[1]*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) + c[2]*(33*
	(7 + 30*c[3])*d[0] + 22*(21 + 38*c[3])*d[1] + 726*d[2] + 990*d[3]-28*c[3]*
	(11*d[2] + d[3])) + c[3]*(11*(-63 + 47*c[3])*d[0] + 2*c[3]*(517*d[1] - 
	7*d[2] - 531*d[3]) + 22*(-9*d[1] + 45*d[2] + 47*d[3])) + 11*c[1]*(3*(35 + 
	14*c[2] - 6*c[3])*d[0] + 2*(9*(7 + 6*c[2] - 2*c[3])*d[1] + (21 + 66*c[2] + 
	38*c[3])*d[2] + (-9 + 38*c[2] + 94*c[3])*d[3]))))/55440.;
		bf_mom[3] = (22*c[0]*(3*(7*(5 - 10*c[1] + 2*c[2] + 6*c[3])*d[0] + 
	2*(-35 + 42*c[1] - 14*c[2] - 6*c[3])*d[1] + 2*(7 - 14*c[1] + 22*c[2] - 30*
	c[3])*d[2]) - 2*(-63 + 18*c[1] + 90*c[2] - 94*c[3])*d[3]) + 231*c[0]*c[0]*
	(5*d[0] + 2*(-5*d[1] + d[2] + 3*d[3])) + 4*(11*c[2]*c[2]*(33*d[0] - 66*d[1]
	 + 26*d[2] + 14*d[3]) + c[2]*(-33*(-7 + 30*c[3])*d[0] + 22*(-21 + 38*c[3])
	*d[1] + 22*(33 + 14*c[3])*d[2] - 2*(495 + 14*c[3])*d[3]) + c[3]*(11*(63 + 
	47*c[3])*d[0] - 22*(9 + 47*c[3])*d[1] - 2*(495 + 7*c[3])*d[2] + 2*(517 + 
	531*c[3])*d[3]) + 99*c[1]*c[1]*(7*d[0] + 2*(-7*d[1] + 3*d[2] + d[3])) - 
	11*c[1]*(3*(35 + 14*c[2] + 6*c[3])*d[0] - 2*(9*(7 + 6*c[2] + 2*c[3])*d[1] 
	+ (-21 - 66*c[2] + 38*c[3])*d[2] + (-9 + 38*c[2] - 94*c[3])*d[3]))))/55440.;
		bf_mom[4] = (231*c[0]*c[0]*(5*d[0] - 6*d[2]) + 4*(11*c[2]*c[2]*
	(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] 
	+ 594*d[1] - 642*c[3]*d[2] - 1166*d[3]) - 22*c[1]*(3*(7 + 2*c[2])*d[1] + 
	c[3]*(27*d[0] - 26*d[2]) - (27 + 26*c[2])*d[3]) + c[2]*(693*d[0] + 572*c[3]
	*d[1] - 1254*d[2] - 1284*c[3]*d[3])) + 22*c[0]*(-21*(5 + 6*c[2])*d[0] + 2*
	(3*(21 + 38*c[2])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1] + 
	53*d[3]))))/13860.;
		bf_mom[5] = (-11*(105*c[0]*c[0] + 84*c[0]*(5*c[1] + c[2] - 3*c[3]) 
	+ 4*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(7*c[2] - 3*c[3]) + 90*
	c[2]*c[3] + 47*c[3]*c[3]))*d[0] - 22*(3*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1]
	 + 72*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(3*c[1] + c[2])) - 4*(9*c[0] + 18*
	c[1] - 38*c[2])*c[3] + 188*c[3]*c[3])*d[1] - 2*(33*(-28 + 7*c[0]*c[0] + 
	28*c[0]*c[1] + 36*c[1]*c[1]) + 1452*(c[0] + 2*c[1])*c[2] + 572*c[2]*c[2] + 
	44*(45*c[0] + 38*c[1] - 14*c[2])*c[3] - 28*c[3]*c[3])*d[2] + 2*(693*c[0]*
	c[0] + 396*c[0]*(c[1] - 5*c[2]) + 44*(9*c[1]*c[1] - 38*c[1]*c[2] + 7*(-9 
	+ c[2]*c[2])) - 4*(517*(c[0] + 2*c[1]) - 14*c[2])*c[3] + 2124*c[3]*c[3])
	*d[3])/ 27720.;
		bf_mom[6] = (231*c[0]*c[0]*(5*d[0] - 6*d[2]) + 4*(11*c[2]*c[2]*
	(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] 
	- 594*d[1] - 642*c[3]*d[2] + 1166*d[3]) - 22*c[1]*(3*(-7 + 2*c[2])*d[1] + 
	c[3]*(27*d[0] - 26*d[2]) + (27 - 26*c[2])*d[3]) + c[2]*(-693*d[0] + 572*
	c[3]*d[1] + 1254*d[2] - 1284*c[3]*d[3])) + 22*c[0]*(-21*(-5 + 6*c[2])*d[0] 
	+ 2*(3*(-21 + 38*c[2])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1] 
	+ 53*d[3]))))/13860.;
		bf_mom[7] = (-11*(21*(-20 + 5*c[0]*c[0] - 20*c[0]*c[1] + 12*c[1]
	*c[1]) + 84*(c[0] - 2*c[1])*c[2] + 132*c[2]*c[2] + 36*(7*c[0] - 2*(c[1] + 
	5*c[2]))*c[3] + 188*c[3]*c[3])*d[0] - 1848*(5*d[1] - d[2] - 3*d[3]) + 2*
	(231*c[0]*c[0]*(5*d[1] - d[2] - 3*d[3]) + 44*c[0]*(c[3]*(9*d[1] + 45*d[2] 
	- 47*d[3]) + c[1]*(-63*d[1] + 21*d[2] + 9*d[3]) + 3*c[2]*(7*d[1] - 11*d[2] 
	+ 15*d[3])) + 4*(c[3]*c[3]*(517*d[1] + 7*d[2] - 531*d[3]) + 11*c[2]*c[2]*
	(33*d[1] - 13*d[2] - 7*d[3]) - 2*c[2]*c[3]*(209*d[1] + 77*d[2] - 7*d[3]) + 
	99*c[1]*c[1]*(7*d[1] - 3*d[2] - d[3]) - 22*c[1]*(c[3]*(9*d[1] + 19*d[2] - 
	47*d[3]) + c[2]*(27*d[1] - 33*d[2] + 19*d[3])))))/27720.;
		bf_mom[8] = (-11*(84*(-5 + c[1]*c[1]) + 3*(35*c[0]*c[0] - 84*c[0]
	*c[2] + 76*c[2]*c[2]) - 216*c[1]*c[3] + 212*c[3]*c[3])*d[0] + 2*(693*c[0]
	*c[0]*d[2] - 44*c[0]*(57*c[2]*d[2] + 3*c[1]*(7*d[1] - 9*d[3]) + c[3]*
	(-27*d[1] + 53*d[3])) + 4*(-693*d[2] + 33*c[1]*c[1]*d[2] + 517*c[2]*c[2]
	*d[2] + 321*c[3]*c[3]*d[2] + c[2]*c[3]*(-286*d[1] + 642*d[3]) - 22*c[1]*
	(13*c[3]*d[2] + c[2]*(-3*d[1] + 13*d[3])))))/6930.;
		break;
	case 100:
	case 102:
		bf_mom[0] = -((a - b)*(a*(231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] - 12*b*
	d[2] + 10*(d[1] + d[2]) - 6*d[3]) + 22*c[0]*(-21*(-5*(3 + 2*c[1] + 2*c[2]) 
	+ 2*b*(5 + 6*c[2]) + 6*c[3])*d[0] + 6*(-35 + 14*(-5 + 2*b)*c[1] + 14*c[2] 
	+ 6*(7 - 6*b)*c[3])* d[1] + 6*(-35 + 14*c[1] - 98*c[2] + b*(42 + 76*c[2]) 
	+ 30*c[3])*d[2] - 2*(-63 + 18*(-7 + 6*b)*c[1] - 90*c[2]+(306 - 212*b)*c[3])*
	 d[3]) + 4*(11*c[2]*c[2]* (3*(-49 + 38*b)*d[0] + 2*(33*d[1] + (81-94*b)*d[2]
	 - 7*d[3])) + 11*c[1]*(3*(-35 + 14*c[2] + 6*(7 - 6*b)*c[3])*d[0] - 6*(-35 +
	 14*c[2] + 2*b*(7 + 2*c[2]) + 6*c[3])*d[1] + 2*(-21 + 66*c[2] + (-90 + 52*b)
	*c[3])*d[2] + 2*(-63 - 90*c[2] + b*(54 + 52*c[2]) + 94*c[3])*d[3]) + c[2]*
	(33*(-35 + 42*b + 30*c[3])*d[0] + 22*(-21 + (-90 + 52*b)*c[3])*d[1] - 22*
	(-147 + 114*b + 14*c[3])*d[2] - 2*(495 + 2*(-649 + 642*b)*c[3])*d[3]) + 33*
	c[1]*c[1]* (7*(-5 + 2*b)*d[0] - 2*(-21*d[1] + (7 + 2*b)*d[2] + 3*d[3])) + 
	c[3]*(11*(63 + (-153 + 106*b)*c[3])*d[0] + 2*(11*(-63 + 54*b + 47*c[3])*d[1]
	 + (-495 + (649 - 642*b)*c[3])*d[2] + (1683 - 1166*b - 531*c[3])*d[3])))) + 
	2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 2*d[2] + 3*d[3]) + 22*c[0]*(-21*(5
	 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*(35 + 56*c[1] - 14*c[2]-24*c[3])*d[1]
	 + 6*(7 - 7*c[1] + 30*c[2] - 15*c[3])*d[2] - (63 + 72*c[1] + 90*c[2] - 200*
	c[3])*d[3]) + 2*(c[3]*(11*(-63 + 100*c[3])*d[0] - 2*c[3]*(517*d[1]+328*d[2]
	 - 531*d[3]) + 22*(36*d[1] + 45*d[2] - 100*d[3])) + 66*c[1]*c[1]*(14*d[0] -
	 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 
	7*d[3]) + 2*c[2]*(-33*(-7 + 15*c[3])*d[0] + 11*(21 + 64*c[3])*d[1] + 22*
	(-45 + 7*c[3])*d[2] + 495*d[3] - 656*c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2]
	 + 24*c[3])*d[0] + 2*(6*(-14 + 8*c[2] + 3*c[3])*d[1] + (21 - 66*c[2] + 64*
	c[3])*d[2] + 2*(18 + 32*c[2] - 47*c[3])*d[3])))) + b*(22*c[0]*(21*(15 - 10*
	c[1] + 10*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3])*d[0] + 6*(7*
	(5 - 5*b - 10*c[1] + 8*b*c[1] + 2*(-1 + b)*c[2]) - 6*(-7 + 4*b)*c[3])*d[1] 
	+ 6*(-35 + 14*b - 14*c[1] + 14*b*c[1] - 98*c[2] + 60*b*c[2] + 30*(-1 + b)*
	c[3])*d[2] - 2*(b*(-63 + 72*c[1] - 90*c[2] - 200*c[3]) + 9*(7 - 14*c[1] + 
	10*c[2] + 34*c[3]))*d[3]) + 231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] + 2*(5*(-1+b)*
	d[1] + (5 - 2*b)*d[2] - 3*(-1 + b)*d[3])) + 4*(c[2]*(33*(7*(-5 + 2*b) + 
	30*(-1 + b)*c[3])*d[0] + 22*(21 - 21*b - 90*c[3] + 64*b*c[3])*d[1] - 22*
	(-147 + 90*b + 14*(-1 + b)*c[3])*d[2] - 2*(495*(-1 + b) + 2*(-649 + 328*b)
	*c[3])*d[3]) + 11*c[2]*c[2]* (3*(-49 + 30*b)*d[0] + 2*(33*(-1 + b)*d[1] + 
	(81 - 34*b)*d[2] - 7*(-1 + b)*d[3])) + 33*c[1]*c[1]*(7*(-5 + 4*b)*d[0] + 2*
	(21*(-1 + b)*d[1] + (-7 + 8*b)*d[2] - 3*(-1 + b)*d[3])) + 11*c[1]*(3*(7*(-1
	 + b)*(-5 + 2*c[2]) - 6*(-7 + 4*b)*c[3])* d[0] + 2*(3*(35 - 28*b - 14*c[2] 
	+ 16*b*c[2] - 6*(-1 + b)*c[3])*d[1] + (3*(-1 + b)*(-7 + 22*c[2]) + 2*(-45 +
	 32*b)*c[3])* d[2] + (-63 + 36*b - 90*c[2] + 64*b*c[2] + 94*(-1 + b)*c[3])
	*d[3])) + c[3]*(11*(63*(-1 + b) + (-153 + 100*b)*c[3])*d[0] + 2*(11*(9*(-7
	 + 4*b) + 47*(-1 + b)*c[3])*d[1] - (495*(-1 + b) + (-649 + 328*b)*c[3])*d[2]
	 - (-1683 + 1100*b + 531*(-1 + b)*c[3])*d[3]))))))/221760.;
		bf_mom[1] = -((a - b)*(2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 2*
	d[2] + 3*d[3]) + 22*c[0]*(-21*(5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*(35 
	+ 56*c[1] - 14*c[2] - 24*c[3])*d[1] + 6*(7 - 7*c[1] + 30*c[2] - 15*c[3])*
	d[2] - (63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(c[3]*(11*(-63 + 100*
	c[3])*d[0] - 2*c[3]*(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*d[1] + 45*d[2]
	 - 100*d[3])) + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]
	*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*(-33*(-7 + 15*c[3])
	*d[0] + 11*(21 + 64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2] + 495*d[3] - 656*
	c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0] + 2*(6*(-14 + 8*c[2]
	 + 3*c[3])*d[1] + (21 - 66*c[2] + 64*c[3])*d[2] + 2*(18 + 32*c[2] - 47*c[3])
	*d[3])))) + a*(231*c[0]*c[0]*(5*(3 + 2*b)*d[0] - 12*b*d[2] - 10*(d[1]+d[2])
	 + 6*d[3]) + 22*c[0]*(-21*(15 + 10*c[1] + 10*c[2] + 2*b*(5 + 6*c[2]) - 6*
	c[3])*d[0] + 6*(35 + 14*(5 + 2*b)*c[1] - 14*c[2] - 6*(7 + 6*b)*c[3])* d[1] 
	+ 6*(35 - 14*c[1] + 98*c[2] + b*(42 + 76*c[2]) - 30*c[3])*d[2] - 2*(63 + 18*
	(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 106*b)*c[3])* d[3]) + 4*(11*c[2]*c[2]* 
	(3*(49 + 38*b)*d[0] - 2*(33*d[1] + (81 + 94*b)*d[2] - 7*d[3])) + 33*c[1]*
	c[1]* (7*(5 + 2*b)*d[0] - 42*d[1] + 2*(7 - 2*b)*d[2] + 6*d[3]) + c[2]*(33*
	(35 + 42*b - 30*c[3])*d[0] + 22*(21 + (90 + 52*b)*c[3])*d[1] - 22*(147 + 
	114*b - 14*c[3])*d[2] + 990*d[3] - 4*(649 + 642*b)*c[3]*d[3]) + 11*c[1]*
	(-3*(-35 + 14*c[2] + 6*(7 + 6*b)*c[3])*d[0] + 2*(-3*(35 - 14*c[2] + 2*b*(7 
	+ 2*c[2]) - 6*c[3])*d[1] + (21 - 66*c[2] + (90 + 52*b)*c[3])*d[2] + (63 + 
	90*c[2] + b*(54 + 52*c[2]) - 94*c[3])*d[3])) + c[3]*(11*(-63 + (153 + 106*b)
	*c[3])*d[0] + 2*(11*(63 + 54*b - 47*c[3])*d[1] + (495 - (649 + 642*b)*c[3])
	*d[2] + (-1683 - 1166*b + 531*c[3])*d[3])))) + b*(22*c[0]*(21*(-15 + 10*c[1]
	 - 10*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3])*d[0] + 6*(7*(-5 
	- 5*b + 10*c[1] + 8*b*c[1] + 2*(1 + b)*c[2]) - 6*(7 + 4*b)*c[3])*d[1] + 6*
	(35 + 14*b + 14*c[1] + 14*b*c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3])*
	d[2] - 2*(9*(-7 - 7*b + 14*c[1] + 8*b*c[1] - 10*(1 + b)*c[2]) - 2*(153 + 
	100*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 2*(5*(1 + b)*d[1] - 
	(5 + 2*b)*d[2] - 3*(1 + b)*d[3])) + 4*(c[2]*(33*(7*(5 + 2*b) + 30*(1 + b)
	*c[3])*d[0] + 22*(-21*(1 + b) + 2*(45 + 32*b)*c[3])*d[1] - 22*(147 + 90*b 
	+ 14*(1 + b)*c[3])*d[2] - 2*(495*(1 + b) + 2*(649 + 328*b)*c[3])*d[3]) + 
	11*c[2]*c[2]* (3*(49 + 30*b)*d[0] + 2*(33*(1 + b)*d[1] - (81 + 34*b)*d[2] 
	- 7*(1 + b)*d[3])) + 33*c[1]*c[1]* (7*(5 + 4*b)*d[0] + 2*(21*(1 + b)*d[1] 
	+ (7 + 8*b)*d[2] - 3*(1 + b)*d[3])) + 11*c[1]*(3*(7*(1 + b)*(-5 + 2*c[2]) 
	- 6*(7 + 4*b)*c[3])* d[0] + 2*(3*(-35 - 28*b + 14*c[2] + 16*b*c[2] - 6*(1 
	+ b)*c[3])*d[1] + (3*(1 + b)*(-7 + 22*c[2]) + 2*(45 + 32*b)*c[3])*d[2] + 
	(63 + 36*b + 90*c[2] + 64*b*c[2] + 94*(1 + b)*c[3])*d[3])) + c[3]*(11*(63*
	(1 + b) + (153 + 100*b)*c[3])*d[0] + 2*(11*(9*(7 + 4*b) + 47*(1 + b)*c[3])
	*d[1] - (495*(1 + b) + (649 + 328*b)*c[3])*d[2] - (1683 + 1100*b + 531*(1 
	+ b)*c[3])*d[3]))))))/221760.;
		bf_mom[2] = -((a - b)*(a*(231*c[0]*c[0]* (5*(3 + 2*b)*d[0] - 12*b*
	d[2] - 10*(d[1] + d[2]) + 6*d[3]) + 22*c[0]*(-21*(-15 + 10*c[1] + 10*c[2] 
	+ 2*b*(-5 + 6*c[2]) - 6*c[3])*d[0] + 6*(14*(5 + 2*b)*c[1] - 7*(5 + 2*c[2]) 
	- 6*(7 + 6*b)*c[3])* d[1] - 6*(35 + 14*c[1] + b*(42 - 76*c[2]) - 98*c[2] + 
	30*c[3])*d[2] - 2*(-63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 106*b)*c[3])
	* d[3]) + 4*(11*c[2]*c[2]* (3*(49 + 38*b)*d[0] - 2*(33*d[1] + (81 + 94*b)*
	d[2] - 7*d[3])) + 33*c[1]*c[1]* (7*(5 + 2*b)*d[0] - 42*d[1] + 2*(7 - 2*b)*
	d[2] + 6*d[3]) + c[3]*(11*(63 + (153 + 106*b)*c[3])*d[0] - 22*(63 + 54*b 
	+ 47*c[3])*d[1] - 2*(495 + (649 + 642*b)*c[3])*d[2] + 2*(1683 + 1166*b + 
	531*c[3])*d[3]) + c[2]*(-33*(35 + 42*b + 30*c[3])*d[0] + 22*((-21 + (90 + 
	52*b)*c[3])*d[1] + (147 + 114*b + 14*c[3])*d[2]) - 2*(495 + 2*(649 + 642*b)
	*c[3])*d[3]) + 11*c[1]*(-3*(35 + 14*c[2] + 6*(7 + 6*b)*c[3])*d[0] + 2*(3*
	(35 + b*(14 - 4*c[2]) + 14*c[2] + 6*c[3])*d[1] + (-21 - 66*c[2] + (90+52*b)
	*c[3])*d[2] + (-63 + 90*c[2] + b*(-54 + 52*c[2]) - 94*c[3])*d[3])))) + 2*a*
	a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 2*d[2] + 3*d[3]) + 22*c[0]*(-21*(-5 + 
	5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*(-35 + 56*c[1] - 14*c[2] - 24*c[3])*d[1]
	 - 6*(7 + 7*c[1] - 30*c[2] + 15*c[3])*d[2] - (-63 + 72*c[1] + 90*c[2] - 200*
	c[3])*d[3]) + 2*(-22*c[2]*(3*(7 + 15*c[3])*d[0] + (21 - 64*c[3])*d[1] - 2*
	(45 + 7*c[3])*d[2]) - 2*c[2]*(495 + 656*c[3])*d[3] + 66*c[1]*c[1]*(14*d[0] 
	- 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 
	7*d[3]) + c[3]*(11*(63 + 100*c[3])*d[0] - 22*(36 + 47*c[3])*d[1] - 2*(495 
	+ 328*c[3])*d[2] + 2*(1100 + 531*c[3])*d[3]) - 11*c[1]*(3*(35 + 14*c[2] + 
	24*c[3])*d[0] - 2*(6*(14 + 8*c[2] + 3*c[3])*d[1] + (-21 - 66*c[2]+64*c[3])
	*d[2] + 2*(-18 + 32*c[2] - 47*c[3])*d[3])))) + b*(22*c[0]*(21*(15 + 10*c[1]
	 - 10*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3])*d[0] + 6*(7*(5 + 
	5*b + 10*c[1] + 8*b*c[1] + 2*(1 + b)*c[2]) - 6*(7 + 4*b)*c[3])*d[1] + 6*
	(-35 - 14*b + 14*c[1] + 14*b*c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3])
	*d[2] - 2*(9*(7 + 7*b + 14*c[1] + 8*b*c[1] - 10*(1 + b)*c[2]) - 2*(153 + 
	100*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 2*(5*(1 + b)*d[1] - 
	(5 + 2*b)*d[2] - 3*(1 + b)*d[3])) + 4*(c[2]*(33*(-7*(5 + 2*b) + 30*(1 + b)
	*c[3])*d[0] + 22*(21 + 21*b + 90*c[3] + 64*b*c[3])*d[1] - 22*(-3*(49 + 30*b)
	 + 14*(1 + b)*c[3])*d[2] - 2*(-495*(1 + b) + 2*(649 + 328*b)*c[3])*d[3]) + 
	11*c[2]*c[2]* (3*(49 + 30*b)*d[0] + 2*(33*(1 + b)*d[1] - (81 + 34*b)*d[2] 
	- 7*(1 + b)*d[3])) + 33*c[1]*c[1]* (7*(5 + 4*b)*d[0] + 2*(21*(1 + b)*d[1] 
	+ (7 + 8*b)*d[2] - 3*(1 + b)*d[3])) + c[3]*(11*(-63*(1 + b) + (153 + 100*b)
	*c[3])*d[0] + 2*(11*(-9*(7 + 4*b) + 47*(1 + b)*c[3])*d[1] + (495*(1 + b) 
	- (649 + 328*b)*c[3])*d[2] + (1683 + 1100*b - 531*(1 + b)*c[3])*d[3])) + 
	11*c[1]*(3*(7*(1 + b)*(5 + 2*c[2]) - 6*(7 + 4*b)*c[3])*d[0] + 2*(3*(35 + 
	28*b + 14*c[2] + 16*b*c[2] - 6*(1 + b)*c[3])* d[1] + (3*(1 + b)*(7+22*c[2])
	 + 2*(45 + 32*b)*c[3])* d[2] + (-63 - 36*b + 90*c[2] + 64*b*c[2] + 94*(1+b)
	*c[3])*d[3]))))))/221760.;
		bf_mom[3] = -((a - b)*(2*a*a*(231*c[0]*c[0]*(5*d[0] - 5*d[1] - 2*
	d[2] + 3*d[3]) + 22*c[0]*(-21*(-5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*
	(-35 + 56*c[1] - 14*c[2] - 24*c[3])*d[1] - 6*(7 + 7*c[1] - 30*c[2] + 15*
	c[3])*d[2] - (-63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(-22*c[2]*(3*
	(7 + 15*c[3])*d[0] + (21 - 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]*
	(495 + 656*c[3])*d[3] + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) 
	+ 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 100*
	c[3])*d[0] - 22*(36 + 47*c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*(1100 + 
	531*c[3])*d[3]) - 11*c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*(14 + 8*
	c[2] + 3*c[3])*d[1] + (-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 32*c[2] - 
	47*c[3])*d[3])))) + a*(231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] - 12*b*d[2] + 10*
	(d[1] + d[2]) - 6*d[3]) + 22*c[0]*(21*(-15 + 10*c[1] + 2*b*(5 - 6*c[2]) + 
	10*c[2] - 6*c[3])*d[0] + 6*(35 + 14*(-5 + 2*b)*c[1] + 14*c[2] + 6*(7 - 6*b)
	*c[3])* d[1] + 6*(35 + 14*c[1] - 98*c[2] + b*(-42 + 76*c[2]) + 30*c[3])*d[2]
	 - 2*(63 + 18*(-7 + 6*b)*c[1] - 90*c[2] + (306 - 212*b)*c[3])* d[3]) + 4*
	(11*c[2]*c[2]* (3*(-49 + 38*b)*d[0] + 2*(33*d[1] + (81 - 94*b)*d[2] - 7*
	d[3])) + c[2]*(-33*(42*b - 5*(7 + 6*c[3]))*d[0] + 22*(21 + (-90 + 52*b)*
	c[3])*d[1] + 22*(114*b - 7*(21 + 2*c[3]))*d[2] + 990*d[3] - 4*(-649 + 642*b)
	*c[3]*d[3]) + 33*c[1]*c[1]* (7*(-5 + 2*b)*d[0] - 2*(-21*d[1] + (7 + 2*b)
	*d[2] + 3*d[3])) + 11*c[1]*(3*(35 + 14*c[2] + 6*(7 - 6*b)*c[3])*d[0] + 2*
	(-3*(35 + 14*c[2] + 2*b*(-7 + 2*c[2]) + 6*c[3])*d[1] + (21 + 66*c[2] + (-90
	 + 52*b)*c[3])*d[2] + (63 - 90*c[2] + b*(-54 + 52*c[2]) + 94*c[3])*d[3])) 
	+ c[3]*(11*(-63 + (-153 + 106*b)*c[3])*d[0] - 2*(11*(-63 + 54*b - 47*c[3])
	*d[1] + (-495 + (-649 + 642*b)*c[3])*d[2] + (1683-1166*b+531*c[3])*d[3]))))
	 + b*(22*c[0]*(21*(-15 - 10*c[1] + 10*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*
	c[3]) + 6*c[3])*d[0] + 6*(7*(-5 + 5*b - 10*c[1] + 8*b*c[1] + 2*(-1+b)*c[2])
	 - 6*(-7 + 4*b)*c[3])*d[1] + 6*(35 - 14*b - 14*c[1] + 14*b*c[1] - 98*c[2] 
	+ 60*b*c[2] + 30*(-1 + b)*c[3])*d[2] - 2*(b*(63 + 72*c[1] - 90*c[2] - 200*
	c[3]) + 9*(-7 - 14*c[1] + 10*c[2] + 34*c[3]))*d[3]) + 231*c[0]*c[0]*(5*(-3 
	+ 2*b)*d[0] + 2*(5*(-1 + b)*d[1] + (5 - 2*b)*d[2] - 3*(-1 + b)*d[3])) + 4*
	(c[2]*(33*(35 - 14*b + 30*(-1 + b)*c[3])*d[0] + 22*(-21 + 21*b - 90*c[3] + 
	64*b*c[3])*d[1] - 22*(147 - 90*b - 14*c[3] + 14*b*c[3])*d[2] - 2*(-495*(-1 
	+ b) + 2*(-649 + 328*b)*c[3])*d[3]) + 11*c[2]*c[2]* (3*(-49 + 30*b)*d[0] + 
	2*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] - 7*(-1 + b)*d[3])) + 33*c[1]*c[1]*
	(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*d[1] + (-7 + 8*b)*d[2] - 3*(-1+b)*d[3]))
	 + c[3]*(11*(63 - 153*c[3] + b*(-63 + 100*c[3]))*d[0] + 2*(11*(63 - 36*b + 
	47*(-1 + b)*c[3])*d[1] + (-495 + 495*b + 649*c[3] - 328*b*c[3])*d[2] + 
	(-1683 + 1100*b - 531*(-1 + b)*c[3])*d[3])) + 11*c[1]*(3*(7*(-1 + b)*(5 + 
	2*c[2]) - 6*(-7 + 4*b)*c[3])* d[0] + 2*(3*(-35 + 28*b - 14*c[2] + 16*b*c[2]
	 - 6*(-1 + b)*c[3])*d[1] + (3*(-1 + b)*(7 + 22*c[2]) + 2*(-45 + 32*b)*c[3])
	*d[2] + (63 - 36*b - 90*c[2] + 64*b*c[2] + 94*(-1 + b)*c[3])
	*d[3]))))))/221760.;
		bf_mom[4] = ((a - b)*(-22*c[0]*(21*(-5*(3 + 2*c[2]) + a*b*(5 + 6*
	c[2]) + a*a*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + b*b*(5 - 5*c[1] + 2*c[2] + 
	3*c[3]))*d[0] - 3*(7*(-20*c[1] + 4*a*b*c[1] + a*a*(5 + 8*c[1] - 2*c[2]) + 
	b*b*(-5 + 8*c[1] + 2*c[2])) - 12*(-7 + 2*a*a + 3*a*b + 2*b*b)*c[3])*d[1] + 
	6*(35 + 98*c[2] - a*b*(21 + 38*c[2]) + a*a*(-7 + 7*c[1] - 30*c[2] + 15*c[3])
	 - b*b*(7 + 7*c[1] + 30*c[2] + 15*c[3]))*d[2] + (9*(-28*c[1] + 12*a*b*c[1] 
	+ b*b*(-7 + 8*c[1] - 10*c[2]) + a*a*(7 + 8*c[1] + 10*c[2])) - 4*(-153 + 50*
	a*a + 53*a*b + 50*b*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)
	*d[0] + 10*d[2] - 6*a*b*d[2] + b*b*(5*d[1] - 2*d[2] - 3*d[3]) + a*a*(-5*d[1]
	 - 2*d[2] + 3*d[3])) + 2*(22*(-3*(35*c[1]*c[1] + 7*c[2]*(5 + 7*c[2]) - 42*
	c[1]*c[3] + 51*c[3]*c[3])*d[0] - 6*(7*c[1]*(-5 + 2*c[2]) + 3*(7 + 10*c[2])*
	c[3])*d[1] - 2*(21*c[1]*c[1] - 3*c[2]*(49 + 27*c[2]) + 90*c[1]*c[3] - 59*
	c[3]*c[3])*d[2] - 2*(9*c[1]*(7 + 10*c[2]) - (153 + 118*c[2])*c[3])*d[3]) + 
	2*a*b*(11*c[2]*c[2]*(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + 
	c[3]*(583*c[3]*d[0] + 594*d[1] - 642*c[3]*d[2] - 1166*d[3]) - 22*c[1]*(3*
	(7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 26*d[2]) - (27 + 26*c[2])*d[3]) + c[2]*
	(693*d[0] + 572*c[3]*d[1] - 1254*d[2] - 1284*c[3]*d[3])) + a*a*(c[3]*(11*
	(-63 + 100*c[3])*d[0] - 2*c[3]*(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*
	d[1] + 45*d[2] - 100*d[3])) + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*
	d[3]) + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*(-33*
	(-7 + 15*c[3])*d[0] + 11*(21 + 64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2] + 
	495*d[3] - 656*c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0] + 2*
	(6*(-14 + 8*c[2] + 3*c[3])*d[1] + (21 - 66*c[2] + 64*c[3])*d[2] + 2*(18 + 
	32*c[2] - 47*c[3])*d[3]))) + b*b*(22*c[2]*(3*(7 + 15*c[3])*d[0] + (-21 + 
	64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) + 22*c[2]*c[2]*(45*d[0] + 33*d[1] - 
	34*d[2] - 7*d[3]) + 66*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2] - 3*d[3]) - 
	2*c[2]*(495 + 656*c[3])*d[3] + 11*c[1]*(3*(-35 + 14*c[2] - 24*c[3])*d[0] + 
	2*(6*(-14 + 8*c[2] - 3*c[3])*d[1] + (-21 + 66*c[2] + 64*c[3])*d[2] + 2*(18
	 + 32*c[2] + 47*c[3])*d[3])) + c[3]*(11*(63 + 100*c[3])*d[0] + 2*c[3]*(517*
	d[1] - 328*d[2] - 531*d[3]) + 22*(36*d[1]-5*(9*d[2]+20*d[3])))))))/55440.;
		bf_mom[5] = ((a - b)*(b*(11*(21*(5*(3 + 2*b)*(-4 + c[0]*c[0]) + 20*
	(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*c[1]*c[1]) - 84*((5 + 2*b)*c[0] - 2*(1 + b)
	*c[1])*c[2] + 12*(49 + 30*b)*c[2]*c[2] - 36*(7*(1 + b)*c[0] + 2*(7 + 4*b)
	*c[1] - 10*(1 + b)*c[2])* c[3] + 4*(153 + 100*b)*c[3]*c[3])*d[0] + 2*(11*
	(3*(35*(1 + b)*(-4 + c[0]*c[0]) + 28*(5 + 4*b)*c[0]*c[1] + 84*(1 + b)*c[1]
	*c[1] + 4*(7*(1 + b)*c[0] + 2*(7 + 8*b)*c[1])*c[2] + 44*(1 + b)*c[2]*c[2]) 
	- 4*(9*(7 + 4*b)*c[0] + 18*(1 + b)*c[1] - 2*(45 + 32*b)*c[2])*c[3] + 188*(1
	 + b)*c[3]*c[3])* d[1] + (33*(-7*(5 + 2*b)*(-4 + c[0]*c[0]) + 28*(1 + b)*
	c[0]*c[1] + 4*(7 + 8*b)*c[1]*c[1]) + 132*((49 + 30*b)*c[0] + 22*(1 + b)*
	c[1])*c[2] - 44*(81 + 34*b)*c[2]*c[2] + 44*(45*(1 + b)*c[0] + (90 + 64*b)*
	c[1] - 14*(1 + b)*c[2])*c[3] - 4*(649 + 328*b)*c[3]*c[3])* d[2] - (99*(7*
	(1 + b)*(-4 + c[0]*c[0]) + 4*(7 + 4*b)*c[0]*c[1] + 4*(1 + b)*c[1]*c[1]) - 
	44*(45*(1 + b)*c[0] + 2*(45 + 32*b)*c[1])*c[2] + 308*(1 + b)*c[2]*c[2] - 4*
	(11*(153 + 100*b)*c[0] + 1034*(1 + b)*c[1] - 2*(649 + 328*b)*c[2])*c[3] + 
	2124*(1 + b)*c[3]*c[3])* d[3])) + a*(11*(21* (5*(3 + 2*b)*(-4 + c[0]*c[0]) 
	- 20*c[0]*c[1] + 4*(5 + 2*b)*c[1]*c[1]) - 84*((5 + 6*b)*c[0] + 2*c[1])*c[2] 
	+ 12*(49 + 38*b)*c[2]*c[2] + 36*(7*c[0] - 2*((7 + 6*b)*c[1] + 5*c[2]))*c[3] 
	+ 4*(153 + 106*b)*c[3]*c[3])*d[0] + 2*(-11*(105*c[0]*c[0] + 12*c[0]*(-7*(5 
	+ 2*b)*c[1] + 7*c[2] + 3*(7 + 6*b)*c[3]) + 4*(-105 + 63*c[1]*c[1] + 6*(-7 
	+ 2*b)*c[1]*c[2] + 33*c[2]*c[2] - 2*(9*c[1] + (45 + 26*b)*c[2])*c[3] + 47*
	c[3]*c[3]))* d[1] - (33*(7*(5 + 6*b)*(-4 + c[0]*c[0]) + 28*c[0]*c[1] + 4*
	(-7 + 2*b)*c[1]*c[1]) - 132*((49 + 38*b)*c[0] - 22*c[1])*c[2] + 44*(81 + 
	94*b)*c[2]*c[2] + 44*(45*c[0] - 2*((45 + 26*b)*c[1] + 7*c[2]))*c[3] + 4*
	(649 + 642*b)*c[3]*c[3])*d[2] + (693*c[0]*c[0] + 44*(9*(-7 + c[1]*c[1]) + 
	2*(45 + 26*b)*c[1]*c[2] + 7*c[2]*c[2]) - 8*(517*c[1] + (649 + 642*b)*c[2])
	*c[3] + 2124*c[3]*c[3] - 44*c[0]*(9*(7 + 6*b)*c[1] + 45*c[2] - (153+106*b)
	*c[3]))* d[3])) + 2*a*a* (11*(21*(-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]
	*c[1]) - 84*(c[0] + c[1])*c[2] + 180*c[2]*c[2] + 18*(7*c[0] - 8*c[1] - 10*
	c[2])*c[3] + 200*c[3]*c[3])*d[0] + 924*(5*d[1] + 2*d[2] - 3*d[3]) - 231*
	c[0]*c[0]*(5*d[1] + 2*d[2] - 3*d[3]) - 44*c[0]*(c[3]*(36*d[1] + 45*d[2] - 
	100*d[3]) + 3*c[2]*(7*d[1] - 30*d[2] + 15*d[3]) + c[1]*(-84*d[1] + 21*d[2]
	 + 36*d[3])) + 4*(2*c[2]*c[3]*(352*d[1] + 77*d[2] - 328*d[3]) - 11*c[2]*
	c[2]*(33*d[1] + 34*d[2] - 7*d[3]) - 33*c[1]*c[1]*(21*d[1] - 8*d[2] - 3*d[3])
	 + c[3]*c[3]*(-517*d[1] - 328*d[2] + 531*d[3]) + 22*c[1]*(c[3]*(9*d[1] + 
	32*d[2] - 47*d[3]) + c[2]*(24*d[1] - 33*d[2] + 32*d[3]))))))/110880.;
		bf_mom[6] = ((a - b)*(-22*c[0]*(21*(15 - 10*c[2] + a*b*(-5 + 6*c[2])
	 + a*a*(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + b*b*(-5 - 5*c[1] + 2*c[2]+3*c[3]))
	*d[0] - 3*(7*(-20*c[1] + 4*a*b*c[1] + a*a*(-5 + 8*c[1] - 2*c[2]) + b*b*(5 
	+ 8*c[1] + 2*c[2])) - 12*(-7 + 2*a*a + 3*a*b + 2*b*b)*c[3])*d[1] + 6*(-35 
	+ a*b*(21 - 38*c[2]) + 98*c[2] + a*a*(7 + 7*c[1] - 30*c[2] + 15*c[3]) - b*b*
	(-7 + 7*c[1] + 30*c[2] + 15*c[3]))*d[2] + (9*(-28*c[1] + 12*a*b*c[1] + b*b*
	(7 + 8*c[1] - 10*c[2]) + a*a*(-7 + 8*c[1] + 10*c[2])) - 4*(-153 + 50*a*a + 
	53*a*b + 50*b*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(-3 + a*a + a*b + b*b)*d[0]
	 + 10*d[2] - 6*a*b*d[2] + b*b*(5*d[1] - 2*d[2] - 3*d[3]) + a*a*(-5*d[1] - 
	2*d[2] + 3*d[3])) + 2*(22*(-3*(35*c[1]*c[1] + 7*c[2]*(-5 + 7*c[2]) - 42*
	c[1]*c[3] + 51*c[3]*c[3])*d[0] - 6*(7*c[1]*(5 + 2*c[2]) + 3*(-7 + 10*c[2])*
	c[3])*d[1] - 2*(21*c[1]*c[1] + 3*(49 - 27*c[2])*c[2] + 90*c[1]*c[3] - 59*
	c[3]*c[3])*d[2] - 2*(9*c[1]*(-7 + 10*c[2]) + (153 - 118*c[2])*c[3])*d[3]) + 
	2*a*b*(11*c[2]*c[2]*(57*d[0] - 94*d[2]) + 33*c[1]*c[1]*(7*d[0] - 2*d[2]) + 
	c[3]*(583*c[3]*d[0] - 594*d[1] - 642*c[3]*d[2] + 1166*d[3]) - 22*c[1]*(3*
	(-7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 26*d[2]) + (27 - 26*c[2])*d[3]) + c[2]*
	(-693*d[0] + 572*c[3]*d[1] + 1254*d[2] - 1284*c[3]*d[3])) + a*a*(-22*c[2]*
	(3*(7 + 15*c[3])*d[0] + (21 - 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]
	*(495 + 656*c[3])*d[3] + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) 
	+ 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 
	100*c[3])*d[0] - 22*(36 + 47*c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*
	(1100 + 531*c[3])*d[3]) - 11*c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*
	(14 + 8*c[2] + 3*c[3])*d[1] + (-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 
	32*c[2] - 47*c[3])*d[3]))) + b*b*(22*c[2]*(3*(-7 + 15*c[3])*d[0] + (21 + 
	64*c[3])*d[1] + 2*(45 - 7*c[3])*d[2]) + 22*c[2]*c[2]*(45*d[0] + 33*d[1] - 
	34*d[2] - 7*d[3]) + 66*c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2] - 3*d[3]) - 
	2*c[2]*(-495 + 656*c[3])*d[3] + c[3]*(11*(-63 + 100*c[3])*d[0] - 792*d[1] 
	+ 990*d[2] + 2*c[3]*(517*d[1] - 328*d[2] - 531*d[3]) + 2200*d[3]) + 11*c[1]*
	(3*(35 + 14*c[2] - 24*c[3])*d[0] + 2*(6*(14 + 8*c[2] - 3*c[3])*d[1] + (21 
	+ 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 32*c[2] + 47*c[3])*d[3]))))))/55440.;
		bf_mom[7] = ((a - b)*(b*(11*(21*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*
	(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)*c[1]*c[1]) - 84*((-5 + 2*b)*c[0] - 2*(-1 
	+ b)*c[1])*c[2] + 12*(-49 + 30*b)*c[2]*c[2] - 36*(7*(-1 + b)*c[0] + 2*(-7 
	+ 4*b)*c[1] - 10*(-1 + b)*c[2])* c[3] + 4*(-153 + 100*b)*c[3]*c[3])*d[0] + 
	2*(11*(3*(35*(-1 + b)*(-4 + c[0]*c[0]) + 28*(-5 + 4*b)*c[0]*c[1] + 84*(-1 
	+ b)*c[1]*c[1] + 4*(7*(-1 + b)*c[0] + 2*(-7 + 8*b)*c[1])*c[2] + 44*(-1 + b)
	*c[2]*c[2]) - 4*(9*(-7 + 4*b)*c[0] + 18*(-1 + b)*c[1] + 2*(45 - 32*b)*c[2])
	*c[3] + 188*(-1 + b)*c[3]*c[3])* d[1] - (33*(7*(-5 + 2*b)*(-4 + c[0]*c[0]) 
	- 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1]) - 132*((-49 + 30*b)*c[0] 
	+ 22*(-1 + b)*c[1])*c[2] + 44*(-81 + 34*b)*c[2]*c[2] - 44*(45*(-1 + b)*c[0]
	 + (-90 + 64*b)*c[1] - 14*(-1 + b)*c[2])*c[3] + 4*(-649 + 328*b)*c[3]*c[3])
	*d[2] - (99*(7*(-1 + b)*(-4 + c[0]*c[0]) + 4*(-7 + 4*b)*c[0]*c[1] + 4*(-1 
	+ b)*c[1]*c[1]) - 44*(45*(-1 + b)*c[0] + 2*(-45 + 32*b)*c[1])*c[2] + 308*
	(-1 + b)*c[2]*c[2] - 4*(11*(-153 + 100*b)*c[0] + 1034*(-1 + b)*c[1] + 2*
	(649 - 328*b)*c[2])*c[3] + 2124*(-1 + b)*c[3]*c[3])*d[3])) + 2*a*a*(11*
	(21* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 84*(c[0] + c[1])*
	c[2] + 180*c[2]*c[2] + 18*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 200*c[3]*c[3])
	*d[0] + 924*(5*d[1] + 2*d[2] - 3*d[3]) - 231*c[0]*c[0]*(5*d[1] + 2*d[2] - 
	3*d[3]) - 44*c[0]*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + 3*c[2]*(7*d[1] - 
	30*d[2] + 15*d[3]) + c[1]*(-84*d[1] + 21*d[2] + 36*d[3])) + 4*(2*c[2]*c[3]*
	(352*d[1] + 77*d[2] - 328*d[3]) - 11*c[2]*c[2]*(33*d[1] + 34*d[2] - 7*d[3])
	 - 33*c[1]*c[1]*(21*d[1] - 8*d[2] - 3*d[3]) + c[3]*c[3]*(-517*d[1] - 328*
	d[2] + 531*d[3]) + 22*c[1]*(c[3]*(9*d[1] + 32*d[2] - 47*d[3]) + c[2]*(24*
	d[1] - 33*d[2] + 32*d[3])))) + a*(11*(21*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 
	20*c[0]*c[1] + 4*(-5 + 2*b)*c[1]*c[1]) - 84*((-5 + 6*b)*c[0] - 2*c[1])*c[2]
	 + 12*(-49 + 38*b)*c[2]*c[2] - 36*(7*c[0] + 2*(-7 + 6*b)*c[1] - 10*c[2])*
	c[3] + 4*(-153 + 106*b)*c[3]*c[3])*d[0] + 1848*(-5*(d[1] + d[2]) + 3*d[3])
	 + 2*(231*c[0]*c[0]*(5*d[1] + (5 - 6*b)*d[2] - 3*d[3]) + 44*c[0]*(3*c[2]*
	(7*d[1] + (-49 + 38*b)*d[2] + 15*d[3]) + 3*c[1]*(7*(-5 + 2*b)*d[1] + 7*d[2]
	 + 3*(7 - 6*b)*d[3]) + c[3]*(-9*(-7 + 6*b)*d[1] + 45*d[2] + (-153 + 106*b)
	*d[3])) + 4*(1386*b*d[2] + c[3]*c[3]* (517*d[1] + (649 - 642*b)*d[2] - 531*
	d[3]) + 11*c[2]*c[2]*(33*d[1] + (81 - 94*b)*d[2] - 7*d[3]) + 33*c[1]*c[1]*
	(21*d[1] - (7 + 2*b)*d[2] - 3*d[3]) + 2*c[2]*c[3]*(11*(-45 + 26*b)*d[1] - 
	77*d[2] + (649 - 642*b)*d[3]) - 22*c[1]*(c[3]*(9*d[1] + (45 - 26*b)*d[2] - 
	47*d[3]) + c[2]*(3*(7 + 2*b)*d[1] - 33*d[2] + 
	(45 - 26*b)*d[3])))))))/ 110880.;
		bf_mom[8] = -((a - b)*(11*(21*(5*(-3 + a*a + a*b + b*b)*(-4 + c[0]*
	c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 4*(-5 + 2*a*a + a*b + 2*b*b)*c[1]*
	c[1]) - 84*((-5 + a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*c[2] + 
	12*(-49 + 15*a*a + 19*a*b + 15*b*b)* c[2]*c[2] + 18* (28*c[1] - 12*a*b*c[1]
	 + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 10*c[2]))*c[3]
	 + 4*(-153 + 50*a*a + 53*a*b + 50*b*b)*c[3]*c[3])* d[0] - 4620*b*b*d[1] + 
	1155*b*b*c[0]*c[0]*d[1] - 9240*c[0]*c[1]*d[1] + 3696*b*b*c[0]*c[1]*d[1] + 
	2772*b*b*c[1]*c[1]*d[1] + 924*b*b*c[0]*c[2]*d[1] - 3696*c[1]*c[2]*d[1] + 
	2112*b*b*c[1]*c[2]*d[1] + 1452*b*b*c[2]*c[2]*d[1] + 5544*c[0]*c[3]*d[1] - 
	1584*b*b*c[0]*c[3]*d[1] - 792*b*b*c[1]*c[3]*d[1] - 7920*c[2]*c[3]*d[1] + 
	2816*b*b*c[2]*c[3]*d[1] + 2068*b*b*c[3]*c[3]*d[1] - 9240*d[2] + 1848*b*b*
	d[2] + 2310*c[0]*c[0]*d[2] - 462*b*b*c[0]*c[0]*d[2] + 924*b*b*c[0]*c[1]*
	d[2] - 1848*c[1]*c[1]*d[2] + 1056*b*b*c[1]*c[1]*d[2] - 12936*c[0]*c[2]*d[2]
	 + 3960*b*b*c[0]*c[2]*d[2] + 2904*b*b*c[1]*c[2]*d[2] + 7128*c[2]*c[2]*d[2] 
	- 1496*b*b*c[2]*c[2]*d[2] + 1980*b*b*c[0]*c[3]*d[2] - 7920*c[1]*c[3]*d[2] 
	+ 2816*b*b*c[1]*c[3]*d[2] - 616*b*b*c[2]*c[3]*d[2] + 5192*c[3]*c[3]*d[2] -
	 1312*b*b*c[3]*c[3]*d[2] - (b*b*(99*(-28 + (c[0] + 2*c[1])*(7*c[0]+2*c[1]))
	 - 44*(45*c[0] + 64*c[1])*c[2] + 308*c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 
	328*c[2])*c[3] + 2124*c[3]*c[3]) - 88*(9*c[0]*(7*c[1] - 17*c[3]) + 2*c[2]*
	(-45*c[1] + 59*c[3])))*d[3] + a*a*(-11*(105*c[0]*c[0] + c[0]*(-336*c[1] + 
	84*c[2] + 144*c[3]) + 4*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] 
	+ 47*c[3]*c[3] - 6*c[1]*(8*c[2] + 3*c[3])))*d[1] - 2*(231*c[0]*c[0] - 44*
	(21 + 12*c[1]*c[1] - 33*c[1]*c[2] - 17*c[2]*c[2]) - 44*(32*c[1] + 7*c[2])*
	c[3] + 656*c[3]*c[3] + 66*c[0]*(7*c[1] + 15*(-2*c[2] + c[3])))* d[2] + 
	(693*c[0]*c[0] - 396*c[0]*(4*c[1] + 5*c[2]) + 44*(-63 + (9*c[1] + c[2])*
	(c[1] + 7*c[2])) + 8*(550*c[0] - 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3])
	* d[3]) + 2*a*b*(-693*c[0]*c[0]*d[2] + 44*c[0]*(57*c[2]*d[2] + 3*c[1]*(7*
	d[1] - 9*d[3]) + c[3]*(-27*d[1] + 53*d[3])) - 4*(-693*d[2] + 33*c[1]*c[1]*
	d[2] + 517*c[2]*c[2]*d[2] + 321*c[3]*c[3]*d[2] + c[2]*c[3]*(-286*d[1] + 
	642*d[3]) - 22*c[1]*(13*c[3]*d[2] + c[2]*(-3*d[1] + 13*d[3]))))))/27720.;
		break;
	case 101:
	case 103:
		bf_mom[0] = -((a - b)*(a*(231*c[0]*c[0]* (5*(-3 + 2*b)*d[0] - 
	12*b*d[2] + 10*(d[1] + d[2]) - 6*d[3]) + 22*c[0]*(-21*(-5*(3 + 2*c[1] + 
	2*c[2]) + 2*b*(5 + 6*c[2]) + 6*c[3])*d[0] + 6*(-35 + 14*(-5 + 2*b)*c[1] 
	+ 14*c[2] + 6*(7 - 6*b)*c[3])* d[1] + 6*(-35 + 14*c[1] - 98*c[2] + b*(42 
	+ 76*c[2]) + 30*c[3])*d[2] - 2*(-63 + 18*(-7 + 6*b)*c[1] - 90*c[2] + (306 
	- 212*b)*c[3])* d[3]) + 4*(11*c[2]*c[2]* (3*(-49 + 38*b)*d[0] + 2*(33*d[1]
	 + (81 - 94*b)*d[2] - 7*d[3])) + 11*c[1]*(3*(-35 + 14*c[2] + 6*(7 - 6*b)
	*c[3])*d[0] - 6*(-35 + 14*c[2] + 2*b*(7 + 2*c[2]) + 6*c[3])*d[1] + 2*(-21 
	+ 66*c[2] + (-90 + 52*b)*c[3])*d[2] + 2*(-63 - 90*c[2] + b*(54 + 52*c[2]) 
	+ 94*c[3])*d[3]) + c[2]*(33*(-35 + 42*b + 30*c[3])*d[0] + 22*(-21 + (-90 +
	 52*b)*c[3])*d[1] - 22*(-147 + 114*b + 14*c[3])*d[2] - 2*(495 + 2*(-649 + 
	642*b)*c[3])*d[3]) + 33*c[1]*c[1]* (7*(-5 + 2*b)*d[0] - 2*(-21*d[1] + (7 +
	 2*b)*d[2] + 3*d[3])) + c[3]*(11*(63 + (-153 + 106*b)*c[3])*d[0] + 2*(11*
	(-63 + 54*b + 47*c[3])*d[1] + (-495 + (649 - 642*b)*c[3])*d[2] + (1683 - 
	1166*b - 531*c[3])*d[3])))) + 2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 2*
	d[2] + 3*d[3]) + 22*c[0]*(-21*(5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*(35 
	+ 56*c[1] - 14*c[2] - 24*c[3])*d[1] + 6*(7 - 7*c[1] + 30*c[2] - 15*c[3])*
	d[2] - (63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(c[3]*(11*(-63 + 100*
	c[3])*d[0] - 2*c[3]*(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*d[1] + 45*
	d[2] - 100*d[3])) + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 
	22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*(-33*(-7 + 
	15*c[3])*d[0] + 11*(21 + 64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2] + 495*d[3]
	 - 656*c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0] + 2*(6*
	(-14 + 8*c[2] + 3*c[3])*d[1] + (21 - 66*c[2] + 64*c[3])*d[2] + 2*(18 + 32*
	c[2] - 47*c[3])*d[3])))) + b*(22*c[0]*(21*(15 - 10*c[1] + 10*c[2] + 2*b*
	(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3])*d[0] + 6*(7*(5 - 5*b - 10*c[1] 
	+ 8*b*c[1] + 2*(-1 + b)*c[2]) - 6*(-7 + 4*b)*c[3])*d[1] + 6*(-35 + 14*b - 
	14*c[1] + 14*b*c[1] - 98*c[2] + 60*b*c[2] + 30*(-1 + b)*c[3])*d[2] - 2*(b*
	(-63 + 72*c[1] - 90*c[2] - 200*c[3]) + 9*(7 - 14*c[1] + 10*c[2]+34*c[3]))*
	d[3]) + 231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] + 2*(5*(-1 + b)*d[1] + (5 - 2*b)*
	d[2] - 3*(-1 + b)*d[3])) + 4*(c[2]*(33*(7*(-5 + 2*b) + 30*(-1 + b)*c[3])*
	d[0] + 22*(21 - 21*b - 90*c[3] + 64*b*c[3])*d[1] - 22*(-147 + 90*b + 14*
	(-1 + b)*c[3])*d[2] - 2*(495*(-1 + b) + 2*(-649 + 328*b)*c[3])*d[3]) + 11*
	c[2]*c[2]* (3*(-49 + 30*b)*d[0] + 2*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] 
	- 7*(-1 + b)*d[3])) + 33*c[1]*c[1]*(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*
	d[1] + (-7 + 8*b)*d[2] - 3*(-1 + b)*d[3])) + 11*c[1]*(3*(7*(-1 + b)*(-5 + 
	2*c[2]) - 6*(-7 + 4*b)*c[3])* d[0] + 2*(3*(35 - 28*b - 14*c[2] + 16*b*c[2]
	 - 6*(-1 + b)*c[3])*d[1] + (3*(-1 + b)*(-7 + 22*c[2]) + 2*(-45 + 32*b)*
	c[3])* d[2] + (-63 + 36*b - 90*c[2] + 64*b*c[2] + 94*(-1 + b)*c[3])*d[3]))
	 + c[3]*(11*(63*(-1 + b) + (-153 + 100*b)*c[3])*d[0] + 2*(11*(9*(-7 + 4*b)
	 + 47*(-1 + b)*c[3])*d[1] - (495*(-1 + b) + (-649 + 328*b)*c[3])*d[2] - 
	(-1683 + 1100*b + 531*(-1 + b)*c[3])*d[3]))))))/221760.;
		bf_mom[1] = -((a - b)*(2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 
	2*d[2] + 3*d[3]) + 22*c[0]*(-21*(-5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*
	(-35 + 56*c[1] - 14*c[2] - 24*c[3])*d[1] - 6*(7 + 7*c[1] - 30*c[2] + 15*
	c[3])*d[2] - (-63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(-22*c[2]*(3*
	(7 + 15*c[3])*d[0] + (21 - 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]*
	(495 + 656*c[3])*d[3] + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3])
	 + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 
	100*c[3])*d[0] - 22*(36 + 47*c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*
	(1100 + 531*c[3])*d[3]) - 11*c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*
	(14 + 8*c[2] + 3*c[3])*d[1] + (-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 
	32*c[2] - 47*c[3])*d[3])))) + a*(231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] - 12*b*
	d[2] + 10*(d[1] + d[2]) - 6*d[3]) + 22*c[0]*(21*(-15 + 10*c[1] + 2*b*(5 - 
	6*c[2]) + 10*c[2] - 6*c[3])*d[0] + 6*(35 + 14*(-5 + 2*b)*c[1] + 14*c[2] + 
	6*(7 - 6*b)*c[3])* d[1] + 6*(35 + 14*c[1] - 98*c[2] + b*(-42 + 76*c[2]) + 
	30*c[3])*d[2] - 2*(63 + 18*(-7 + 6*b)*c[1] - 90*c[2] + (306-212*b)*c[3])*
	 d[3]) + 4*(11*c[2]*c[2]* (3*(-49 + 38*b)*d[0] + 2*(33*d[1] + (81 - 94*b)
	*d[2] - 7*d[3])) + c[2]*(-33*(42*b - 5*(7 + 6*c[3]))*d[0] + 22*(21 + (-90
	 + 52*b)*c[3])*d[1] + 22*(114*b - 7*(21 + 2*c[3]))*d[2] + 990*d[3] - 4*
	(-649 + 642*b)*c[3]*d[3]) + 33*c[1]*c[1]* (7*(-5 + 2*b)*d[0] - 2*(-21*
	d[1] + (7 + 2*b)*d[2] + 3*d[3])) + 11*c[1]*(3*(35 + 14*c[2] + 6*(7 - 6*b)
	*c[3])*d[0] + 2*(-3*(35 + 14*c[2] + 2*b*(-7 + 2*c[2]) + 6*c[3])*d[1] + 
	(21 + 66*c[2] + (-90 + 52*b)*c[3])*d[2] + (63 - 90*c[2] + b*(-54 +52*c[2])
	 + 94*c[3])*d[3])) + c[3]*(11*(-63 + (-153 + 106*b)*c[3])*d[0] - 2*(11*
	(-63 + 54*b - 47*c[3])*d[1] + (-495 + (-649 + 642*b)*c[3])*d[2] + (1683 
	- 1166*b + 531*c[3])*d[3])))) + b*(22*c[0]*(21*(-15 - 10*c[1] + 10*c[2] + 
	2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3])*d[0] + 6*(7*(-5 + 5*b - 10*
	c[1] + 8*b*c[1] + 2*(-1 + b)*c[2]) - 6*(-7 + 4*b)*c[3])*d[1] + 6*(35 - 
	14*b - 14*c[1] + 14*b*c[1] - 98*c[2] + 60*b*c[2] + 30*(-1 + b)*c[3])*d[2] 
	- 2*(b*(63 + 72*c[1] - 90*c[2] - 200*c[3]) + 9*(-7 - 14*c[1] + 10*c[2] + 
	34*c[3]))*d[3]) + 231*c[0]*c[0]*(5*(-3 + 2*b)*d[0] + 2*(5*(-1 + b)*d[1] + 
	(5 - 2*b)*d[2] - 3*(-1 + b)*d[3])) + 4*(c[2]*(33*(35 - 14*b + 30*(-1 + b)*
	c[3])*d[0] + 22*(-21 + 21*b - 90*c[3] + 64*b*c[3])*d[1] - 22*(147 - 90*b 
	- 14*c[3] + 14*b*c[3])*d[2] - 2*(-495*(-1 + b) + 2*(-649 + 328*b)*c[3])*
	d[3]) + 11*c[2]*c[2]* (3*(-49 + 30*b)*d[0] + 2*(33*(-1 + b)*d[1] + (81 - 
	34*b)*d[2] - 7*(-1 + b)*d[3])) + 33*c[1]*c[1]*(7*(-5 + 4*b)*d[0] + 2*(21*
	(-1 + b)*d[1] + (-7 + 8*b)*d[2] - 3*(-1 + b)*d[3])) + c[3]*(11*(63 - 153*
	c[3] + b*(-63 + 100*c[3]))*d[0] + 2*(11*(63 - 36*b + 47*(-1 + b)*c[3])*
	d[1] + (-495 + 495*b + 649*c[3] - 328*b*c[3])*d[2] + (-1683 + 1100*b - 
	531*(-1 + b)*c[3])*d[3])) + 11*c[1]*(3*(7*(-1 + b)*(5 + 2*c[2]) - 6*(-7 
	+ 4*b)*c[3])* d[0] + 2*(3*(-35 + 28*b - 14*c[2] + 16*b*c[2] - 6*(-1 + b)*
	c[3])*d[1] + (3*(-1 + b)*(7 + 22*c[2]) + 2*(-45 + 32*b)*c[3])*d[2] + (63 
	- 36*b - 90*c[2] + 64*b*c[2] + 94*(-1 + b)*c[3])*d[3]))))))/221760.;
		bf_mom[2] = -((a - b)*(a*(231*c[0]*c[0]* (5*(3 + 2*b)*d[0] - 12*
	b*d[2] - 10*(d[1] + d[2]) + 6*d[3]) + 22*c[0]*(-21*(-15 + 10*c[1] + 10*
	c[2] + 2*b*(-5 + 6*c[2]) - 6*c[3])*d[0] + 6*(14*(5 + 2*b)*c[1] - 7*(5 + 
	2*c[2]) - 6*(7 + 6*b)*c[3])* d[1] - 6*(35 + 14*c[1] + b*(42 - 76*c[2]) - 
	98*c[2] + 30*c[3])*d[2] - 2*(-63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 
	+ 106*b)*c[3])* d[3]) + 4*(11*c[2]*c[2]* (3*(49 + 38*b)*d[0] - 2*(33*d[1]
	 + (81 + 94*b)*d[2] - 7*d[3])) + 33*c[1]*c[1]* (7*(5 + 2*b)*d[0] - 42*d[1]
	 + 2*(7 - 2*b)*d[2] + 6*d[3]) + c[3]*(11*(63 + (153 + 106*b)*c[3])*d[0] - 
	22*(63 + 54*b + 47*c[3])*d[1] - 2*(495 + (649 + 642*b)*c[3])*d[2] + 2*
	(1683 + 1166*b + 531*c[3])*d[3]) + c[2]*(-33*(35 + 42*b + 30*c[3])*d[0] + 
	22*((-21 + (90 + 52*b)*c[3])*d[1] + (147 + 114*b + 14*c[3])*d[2]) - 2*
	(495 + 2*(649 + 642*b)*c[3])*d[3]) + 11*c[1]*(-3*(35 + 14*c[2] + 6*(7 + 
	6*b)*c[3])*d[0] + 2*(3*(35 + b*(14 - 4*c[2]) + 14*c[2] + 6*c[3])*d[1] + 
	(-21 - 66*c[2] + (90 + 52*b)*c[3])*d[2] + (-63 + 90*c[2] + b*(-54 + 52*
	c[2]) - 94*c[3])*d[3])))) + 2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 2*
	d[2] + 3*d[3]) + 22*c[0]*(-21*(-5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*
	(-35 + 56*c[1] - 14*c[2] - 24*c[3])*d[1] - 6*(7 + 7*c[1] - 30*c[2] + 15*
	c[3])*d[2] - (-63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(-22*c[2]*(3*
	(7 + 15*c[3])*d[0] + (21 - 64*c[3])*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]*
	(495 + 656*c[3])*d[3] + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3])
	 + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 
	100*c[3])*d[0] - 22*(36 + 47*c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*
	(1100 + 531*c[3])*d[3]) - 11*c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*
	(14 + 8*c[2] + 3*c[3])*d[1] + (-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 
	32*c[2] - 47*c[3])*d[3])))) + b*(22*c[0]*(21*(15 + 10*c[1] - 10*c[2] + 
	2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3])*d[0] + 6*(7*(5 + 5*b + 10*
	c[1] + 8*b*c[1] + 2*(1 + b)*c[2]) - 6*(7 + 4*b)*c[3])*d[1] + 6*(-35 - 
	14*b + 14*c[1] + 14*b*c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3])*d[2] 
	- 2*(9*(7 + 7*b + 14*c[1] + 8*b*c[1] - 10*(1 + b)*c[2]) - 2*(153 + 100*b)
	*c[3])*d[3]) + 231*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 2*(5*(1 + b)*d[1] - (5 
	+ 2*b)*d[2] - 3*(1 + b)*d[3])) + 4*(c[2]*(33*(-7*(5 + 2*b) + 30*(1 + b)*
	c[3])*d[0] + 22*(21 + 21*b + 90*c[3] + 64*b*c[3])*d[1] - 22*(-3*(49+30*b)
	 + 14*(1 + b)*c[3])*d[2] - 2*(-495*(1 + b) + 2*(649 + 328*b)*c[3])*d[3]) 
	+ 11*c[2]*c[2]* (3*(49 + 30*b)*d[0] + 2*(33*(1 + b)*d[1] - (81 + 34*b)*
	d[2] - 7*(1 + b)*d[3])) + 33*c[1]*c[1]* (7*(5 + 4*b)*d[0] + 2*(21*(1+b)*
	d[1] + (7 + 8*b)*d[2] - 3*(1 + b)*d[3])) + c[3]*(11*(-63*(1 + b) + (153 
	+ 100*b)*c[3])*d[0] + 2*(11*(-9*(7 + 4*b) + 47*(1 + b)*c[3])*d[1] + (495*
	(1 + b) - (649 + 328*b)*c[3])*d[2] + (1683 + 1100*b - 531*(1 + b)*c[3])*
	d[3])) + 11*c[1]*(3*(7*(1 + b)*(5 + 2*c[2]) - 6*(7 + 4*b)*c[3])*d[0]+2*
	(3*(35 + 28*b + 14*c[2] + 16*b*c[2] - 6*(1 + b)*c[3])* d[1] + (3*(1 + b)*
	(7 + 22*c[2]) + 2*(45 + 32*b)*c[3])* d[2] + (-63 - 36*b + 90*c[2] + 64*b*
	c[2] + 94*(1 + b)*c[3])*d[3]))))))/221760.;
		bf_mom[3] = -((a - b)*(2*a*a*(231*c[0]*c[0]* (5*d[0] - 5*d[1] - 
	2*d[2] + 3*d[3]) + 22*c[0]*(-21*(5 + 5*c[1] + 2*c[2] - 3*c[3])*d[0] + 3*
	(35 + 56*c[1] - 14*c[2] - 24*c[3])*d[1] + 6*(7 - 7*c[1] + 30*c[2] - 15*
	c[3])*d[2] - (63 + 72*c[1] + 90*c[2] - 200*c[3])*d[3]) + 2*(c[3]*(11*
	(-63 + 100*c[3])*d[0] - 2*c[3]*(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*
	d[1] + 45*d[2] - 100*d[3])) + 66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 
	3*d[3]) + 22*c[2]*c[2]*(45*d[0] - 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*
	(-33*(-7 + 15*c[3])*d[0] + 11*(21 + 64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2]
	 + 495*d[3] - 656*c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0]
	 + 2*(6*(-14 + 8*c[2] + 3*c[3])*d[1] + (21 - 66*c[2] + 64*c[3])*d[2] + 2*
	(18 + 32*c[2] - 47*c[3])*d[3])))) + a*(231*c[0]*c[0]*(5*(3 + 2*b)*d[0] - 
	12*b*d[2] - 10*(d[1] + d[2]) + 6*d[3]) + 22*c[0]*(-21*(15 + 10*c[1] + 10*
	c[2] + 2*b*(5 + 6*c[2]) - 6*c[3])*d[0] + 6*(35 + 14*(5 + 2*b)*c[1] - 14*
	c[2] - 6*(7 + 6*b)*c[3])* d[1] + 6*(35 - 14*c[1] + 98*c[2] + b*(42 + 76*
	c[2]) - 30*c[3])*d[2] - 2*(63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 
	106*b)*c[3])* d[3]) + 4*(11*c[2]*c[2]* (3*(49 + 38*b)*d[0] - 2*(33*d[1] +
	 (81 + 94*b)*d[2] - 7*d[3])) + 33*c[1]*c[1]* (7*(5 + 2*b)*d[0] - 42*d[1] 
	+ 2*(7 - 2*b)*d[2] + 6*d[3]) + c[2]*(33*(35 + 42*b - 30*c[3])*d[0] + 22*
	(21 + (90 + 52*b)*c[3])*d[1] - 22*(147 + 114*b - 14*c[3])*d[2] + 990*d[3]
	 - 4*(649 + 642*b)*c[3]*d[3]) + 11*c[1]*(-3*(-35 + 14*c[2] + 6*(7 + 6*b)*
	c[3])*d[0] + 2*(-3*(35 - 14*c[2] + 2*b*(7 + 2*c[2]) - 6*c[3])*d[1] + (21 
	- 66*c[2] + (90 + 52*b)*c[3])*d[2] + (63 + 90*c[2] + b*(54 + 52*c[2]) - 
	94*c[3])*d[3])) + c[3]*(11*(-63 + (153 + 106*b)*c[3])*d[0] + 2*(11*(63 + 
	54*b - 47*c[3])*d[1] + (495 - (649 + 642*b)*c[3])*d[2] + (-1683 - 1166*b 
	+ 531*c[3])*d[3])))) + b*(22*c[0]*(21*(-15 + 10*c[1] - 10*c[2] + 2*b*(-5 
	+ 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3])*d[0] + 6*(7*(-5 - 5*b + 10*c[1] + 
	8*b*c[1] + 2*(1 + b)*c[2]) - 6*(7 + 4*b)*c[3])*d[1] + 6*(35 + 14*b + 14*
	c[1] + 14*b*c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3])*d[2] - 2*(9*(-7
	 - 7*b + 14*c[1] + 8*b*c[1] - 10*(1 + b)*c[2]) - 2*(153 + 100*b)*c[3])*
	d[3]) + 231*c[0]*c[0]*(5*(3 + 2*b)*d[0] + 2*(5*(1 + b)*d[1] - (5 + 2*b)*
	d[2] - 3*(1 + b)*d[3])) + 4*(c[2]*(33*(7*(5 + 2*b) + 30*(1 + b)*c[3])*d[0]
	 + 22*(-21*(1 + b) + 2*(45 + 32*b)*c[3])*d[1] - 22*(147 + 90*b + 14*(1+b)
	*c[3])*d[2] - 2*(495*(1 + b) + 2*(649 + 328*b)*c[3])*d[3]) + 11*c[2]*c[2]*
	 (3*(49 + 30*b)*d[0] + 2*(33*(1 + b)*d[1] - (81 + 34*b)*d[2] - 7*(1+b)*
	d[3])) + 33*c[1]*c[1]* (7*(5 + 4*b)*d[0] + 2*(21*(1 + b)*d[1] + (7+8*b)*
	d[2] - 3*(1 + b)*d[3])) + 11*c[1]*(3*(7*(1 + b)*(-5 + 2*c[2]) - 6*(7+4*b)
	*c[3])* d[0] + 2*(3*(-35 - 28*b + 14*c[2] + 16*b*c[2] - 6*(1 + b)*c[3])*
	d[1] + (3*(1 + b)*(-7 + 22*c[2]) + 2*(45 + 32*b)*c[3])*d[2] + (63 + 36*b 
	+ 90*c[2] + 64*b*c[2] + 94*(1 + b)*c[3])*d[3])) + c[3]*(11*(63*(1 + b) + 
	(153 + 100*b)*c[3])*d[0] + 2*(11*(9*(7 + 4*b) + 47*(1 + b)*c[3])*d[1] - 
	(495*(1 + b) + (649 + 328*b)*c[3])*d[2] - (1683 + 1100*b + 531*(1 + b)*
	c[3])*d[3]))))))/221760.;
		bf_mom[4] = ((a - b)*(b*(11*(21*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 
	20*(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)*c[1]*c[1]) - 84*((-5 + 2*b)*c[0] - 
	2*(-1 + b)*c[1])*c[2] + 12*(-49 + 30*b)*c[2]*c[2] - 36*(7*(-1 + b)*c[0] + 
	2*(-7 + 4*b)*c[1] - 10*(-1 + b)*c[2])* c[3] + 4*(-153 + 100*b)*c[3]*c[3])
	*d[0] + 2*(11*(3*(35*(-1 + b)*(-4 + c[0]*c[0]) + 28*(-5 + 4*b)*c[0]*c[1]
	 + 84*(-1 + b)*c[1]*c[1] + 4*(7*(-1 + b)*c[0] + 2*(-7 + 8*b)*c[1])*c[2] +
	 44*(-1 + b)*c[2]*c[2]) - 4*(9*(-7 + 4*b)*c[0] + 18*(-1 + b)*c[1] + 2*
	(45 - 32*b)*c[2])*c[3] + 188*(-1 + b)*c[3]*c[3])* d[1] - (33*(7*(-5 + 2*b)
	*(-4 + c[0]*c[0]) - 28*(-1 + b)*c[0]*c[1] - 4*(-7 + 8*b)*c[1]*c[1]) - 132*
	((-49 + 30*b)*c[0] + 22*(-1 + b)*c[1])*c[2] + 44*(-81 + 34*b)*c[2]*c[2] - 
	44*(45*(-1 + b)*c[0] + (-90 + 64*b)*c[1] - 14*(-1 + b)*c[2])*c[3] + 4*
	(-649 + 328*b)*c[3]*c[3])*d[2] - (99*(7*(-1 + b)*(-4 + c[0]*c[0]) + 4*(-7 
	+ 4*b)*c[0]*c[1] + 4*(-1 + b)*c[1]*c[1]) - 44*(45*(-1 + b)*c[0] + 2*(-45 
	+ 32*b)*c[1])*c[2] + 308*(-1 + b)*c[2]*c[2] - 4*(11*(-153 + 100*b)*c[0] + 
	1034*(-1 + b)*c[1] + 2*(649 - 328*b)*c[2])*c[3] + 2124*(-1 + b)*c[3]*c[3])
	*d[3])) + 2*a*a*(11*(21* (-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1])
	 - 84*(c[0] + c[1])*c[2] + 180*c[2]*c[2] + 18*(7*c[0] - 8*c[1] - 10*c[2])*
	c[3] + 200*c[3]*c[3])*d[0] + 924*(5*d[1] + 2*d[2] - 3*d[3]) - 231*c[0]*
	c[0]*(5*d[1] + 2*d[2] - 3*d[3]) - 44*c[0]*(c[3]*(36*d[1] + 45*d[2] - 100*
	d[3]) + 3*c[2]*(7*d[1] - 30*d[2] + 15*d[3]) + c[1]*(-84*d[1] + 21*d[2] + 
	36*d[3])) + 4*(2*c[2]*c[3]*(352*d[1] + 77*d[2] - 328*d[3]) - 11*c[2]*c[2]*
	(33*d[1] + 34*d[2] - 7*d[3]) - 33*c[1]*c[1]*(21*d[1] - 8*d[2] - 3*d[3]) + 
	c[3]*c[3]*(-517*d[1] - 328*d[2] + 531*d[3]) + 22*c[1]*(c[3]*(9*d[1] + 32*
	d[2] - 47*d[3]) + c[2]*(24*d[1] - 33*d[2] + 32*d[3])))) + a*(11*(21*(5*
	(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*c[0]*c[1] + 4*(-5 + 2*b)*c[1]*c[1]) - 84*
	((-5 + 6*b)*c[0] - 2*c[1])*c[2] + 12*(-49 + 38*b)*c[2]*c[2] - 36*(7*c[0] 
	+ 2*(-7 + 6*b)*c[1] - 10*c[2])*c[3] + 4*(-153 + 106*b)*c[3]*c[3])*d[0] + 
	1848*(-5*(d[1] + d[2]) + 3*d[3]) + 2*(231*c[0]*c[0]*(5*d[1] + (5 - 6*b)*
	d[2] - 3*d[3]) + 44*c[0]*(3*c[2]*(7*d[1] + (-49 + 38*b)*d[2] + 15*d[3]) + 
	3*c[1]*(7*(-5 + 2*b)*d[1] + 7*d[2] + 3*(7 - 6*b)*d[3]) + c[3]*(-9*(-7 + 
	6*b)*d[1] + 45*d[2] + (-153 + 106*b)*d[3])) + 4*(1386*b*d[2] + c[3]*c[3]*
	 (517*d[1] + (649 - 642*b)*d[2] - 531*d[3]) + 11*c[2]*c[2]*(33*d[1] + (81 
	- 94*b)*d[2] - 7*d[3]) + 33*c[1]*c[1]*(21*d[1] - (7 + 2*b)*d[2] - 3*d[3])
	 + 2*c[2]*c[3]*(11*(-45 + 26*b)*d[1] - 77*d[2] + (649 - 642*b)*d[3]) - 22*
	c[1]*(c[3]*(9*d[1] + (45 - 26*b)*d[2] - 47*d[3]) + c[2]*(3*(7 + 2*b)*d[1] 
	- 33*d[2] + (45 - 26*b)*d[3])))))))/ 110880.;
		bf_mom[5] = ((a - b)*(-22*c[0]*(21*(15 - 10*c[2] + a*b*(-5+6*c[2])
	 + a*a*(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + b*b*(-5 - 5*c[1] + 2*c[2] + 3*
	c[3]))*d[0] - 3*(7*(-20*c[1] + 4*a*b*c[1] + a*a*(-5 + 8*c[1] - 2*c[2]) + 
	b*b*(5 + 8*c[1] + 2*c[2])) - 12*(-7 + 2*a*a + 3*a*b + 2*b*b)*c[3])*d[1] + 
	6*(-35 + a*b*(21 - 38*c[2]) + 98*c[2] + a*a*(7 + 7*c[1] - 30*c[2] + 15*
	c[3]) - b*b*(-7 + 7*c[1] + 30*c[2] + 15*c[3]))*d[2] + (9*(-28*c[1] + 12*a*
	b*c[1] + b*b*(7 + 8*c[1] - 10*c[2]) + a*a*(-7 + 8*c[1] + 10*c[2])) - 4*
	(-153 + 50*a*a + 53*a*b + 50*b*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(-3 + 
	a*a + a*b + b*b)*d[0] + 10*d[2] - 6*a*b*d[2] + b*b*(5*d[1] - 2*d[2] - 
	3*d[3]) + a*a*(-5*d[1] - 2*d[2] + 3*d[3])) + 2*(22*(-3*(35*c[1]*c[1] + 
	7*c[2]*(-5 + 7*c[2]) - 42*c[1]*c[3] + 51*c[3]*c[3])*d[0] - 6*(7*c[1]*(5 +
	 2*c[2]) + 3*(-7 + 10*c[2])*c[3])*d[1] - 2*(21*c[1]*c[1] + 3*(49-27*c[2])
	*c[2] + 90*c[1]*c[3] - 59*c[3]*c[3])*d[2] - 2*(9*c[1]*(-7 + 10*c[2]) + 
	(153 - 118*c[2])*c[3])*d[3]) + 2*a*b*(11*c[2]*c[2]*(57*d[0] - 94*d[2]) + 
	33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] - 594*d[1] - 642*c[3]
	*d[2] + 1166*d[3]) - 22*c[1]*(3*(-7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 26*
	d[2]) + (27 - 26*c[2])*d[3]) + c[2]*(-693*d[0] + 572*c[3]*d[1] + 1254*d[2]
	 - 1284*c[3]*d[3])) + a*a*(-22*c[2]*(3*(7 + 15*c[3])*d[0] + (21 - 64*c[3])
	*d[1] - 2*(45 + 7*c[3])*d[2]) - 2*c[2]*(495 + 656*c[3])*d[3] + 66*c[1]*
	c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0] - 33*
	d[1] - 34*d[2] + 7*d[3]) + c[3]*(11*(63 + 100*c[3])*d[0] - 22*(36 + 47*
	c[3])*d[1] - 2*(495 + 328*c[3])*d[2] + 2*(1100 + 531*c[3])*d[3]) - 11*
	c[1]*(3*(35 + 14*c[2] + 24*c[3])*d[0] - 2*(6*(14 + 8*c[2] + 3*c[3])*d[1] 
	+ (-21 - 66*c[2] + 64*c[3])*d[2] + 2*(-18 + 32*c[2] - 47*c[3])*d[3]))) + 
	b*b*(22*c[2]*(3*(-7 + 15*c[3])*d[0] + (21 + 64*c[3])*d[1] + 2*(45-7*c[3])
	*d[2]) + 22*c[2]*c[2]*(45*d[0] + 33*d[1] - 34*d[2] - 7*d[3]) + 66*c[1]*
	c[1]*(14*d[0] + 21*d[1] + 8*d[2] - 3*d[3]) - 2*c[2]*(-495 + 656*c[3])*d[3]
	 + c[3]*(11*(-63 + 100*c[3])*d[0] - 792*d[1] + 990*d[2] + 2*c[3]*(517*d[1]
	 - 328*d[2] - 531*d[3]) + 2200*d[3]) + 11*c[1]*(3*(35 + 14*c[2]-24*c[3])*
	d[0] + 2*(6*(14 + 8*c[2] - 3*c[3])*d[1] + (21 + 66*c[2] + 64*c[3])*d[2] + 
	2*(-18 + 32*c[2] + 47*c[3])*d[3]))))))/55440.;
		bf_mom[6] = ((a - b)*(b*(11*(21*(5*(3 + 2*b)*(-4 + c[0]*c[0]) + 
	20*(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*c[1]*c[1]) - 84*((5 + 2*b)*c[0] - 2*
	(1 + b)*c[1])*c[2] + 12*(49 + 30*b)*c[2]*c[2] - 36*(7*(1 + b)*c[0] + 2*
	(7 + 4*b)*c[1] - 10*(1 + b)*c[2])* c[3] + 4*(153 + 100*b)*c[3]*c[3])*d[0]
	 + 2*(11*(3*(35*(1 + b)*(-4 + c[0]*c[0]) + 28*(5 + 4*b)*c[0]*c[1] + 84*(1
	 + b)*c[1]*c[1] + 4*(7*(1 + b)*c[0] + 2*(7 + 8*b)*c[1])*c[2] + 44*(1 + b)*
	c[2]*c[2]) - 4*(9*(7 + 4*b)*c[0] + 18*(1 + b)*c[1] - 2*(45 + 32*b)*c[2])*
	c[3] + 188*(1 + b)*c[3]*c[3])* d[1] + (33*(-7*(5 + 2*b)*(-4 + c[0]*c[0]) 
	+ 28*(1 + b)*c[0]*c[1] + 4*(7 + 8*b)*c[1]*c[1]) + 132*((49 + 30*b)*c[0] +
	 22*(1 + b)*c[1])*c[2] - 44*(81 + 34*b)*c[2]*c[2] + 44*(45*(1 + b)*c[0] +
	 (90 + 64*b)*c[1] - 14*(1 + b)*c[2])*c[3] - 4*(649 + 328*b)*c[3]*c[3])*
	 d[2] - (99*(7*(1 + b)*(-4 + c[0]*c[0]) + 4*(7 + 4*b)*c[0]*c[1] + 4*(1+b)*
	c[1]*c[1]) - 44*(45*(1 + b)*c[0] + 2*(45 + 32*b)*c[1])*c[2] + 308*(1 + b)*
	c[2]*c[2] - 4*(11*(153 + 100*b)*c[0] + 1034*(1 + b)*c[1] - 2*(649 + 328*b)
	*c[2])*c[3] + 2124*(1 + b)*c[3]*c[3])* d[3])) + a*(11*(21* (5*(3 + 2*b)*
	(-4 + c[0]*c[0]) - 20*c[0]*c[1] + 4*(5 + 2*b)*c[1]*c[1]) - 84*((5 + 6*b)*
	c[0] + 2*c[1])*c[2] + 12*(49 + 38*b)*c[2]*c[2] + 36*(7*c[0] - 2*((7 + 6*b)
	*c[1] + 5*c[2]))*c[3] + 4*(153 + 106*b)*c[3]*c[3])*d[0] + 2*(-11*(105*
	c[0]*c[0] + 12*c[0]*(-7*(5 + 2*b)*c[1] + 7*c[2] + 3*(7 + 6*b)*c[3]) + 4*
	(-105 + 63*c[1]*c[1] + 6*(-7 + 2*b)*c[1]*c[2] + 33*c[2]*c[2] - 2*(9*c[1] 
	+ (45 + 26*b)*c[2])*c[3] + 47*c[3]*c[3]))* d[1] - (33*(7*(5 + 6*b)*(-4 + 
	c[0]*c[0]) + 28*c[0]*c[1] + 4*(-7 + 2*b)*c[1]*c[1]) - 132*((49 + 38*b)*
	c[0] - 22*c[1])*c[2] + 44*(81 + 94*b)*c[2]*c[2] + 44*(45*c[0] - 2*((45 + 
	26*b)*c[1] + 7*c[2]))*c[3] + 4*(649 + 642*b)*c[3]*c[3])*d[2] + (693*c[0]*
	c[0] + 44*(9*(-7 + c[1]*c[1]) + 2*(45 + 26*b)*c[1]*c[2] + 7*c[2]*c[2]) - 
	8*(517*c[1] + (649 + 642*b)*c[2])*c[3] + 2124*c[3]*c[3] - 44*c[0]*(9*(7 + 
	6*b)*c[1] + 45*c[2] - (153 + 106*b)*c[3]))* d[3])) + 2*a*a* (11*(21*(-20 
	+ 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 84*(c[0] + c[1])*c[2] + 
	180*c[2]*c[2] + 18*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 200*c[3]*c[3])*d[0] 
	+ 924*(5*d[1] + 2*d[2] - 3*d[3]) - 231*c[0]*c[0]*(5*d[1] + 2*d[2]-3*d[3])
	 - 44*c[0]*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + 3*c[2]*(7*d[1]-30*d[2]
	 + 15*d[3]) + c[1]*(-84*d[1] + 21*d[2] + 36*d[3])) + 4*(2*c[2]*c[3]*(352*
	d[1] + 77*d[2] - 328*d[3]) - 11*c[2]*c[2]*(33*d[1] + 34*d[2] - 7*d[3]) - 
	33*c[1]*c[1]*(21*d[1] - 8*d[2] - 3*d[3]) + c[3]*c[3]*(-517*d[1] - 328*d[2]
	 + 531*d[3]) + 22*c[1]*(c[3]*(9*d[1] + 32*d[2] - 47*d[3]) + c[2]*(24*d[1]
	 - 33*d[2] + 32*d[3]))))))/110880.;
		bf_mom[7] = ((a - b)*(-22*c[0]*(21*(-5*(3 + 2*c[2]) + a*b*(5 + 
	6*c[2]) + a*a*(5 + 5*c[1] + 2*c[2] - 3*c[3]) + b*b*(5 - 5*c[1] + 2*c[2] +
	 3*c[3]))*d[0] - 3*(7*(-20*c[1] + 4*a*b*c[1] + a*a*(5 + 8*c[1] - 2*c[2]) 
	+ b*b*(-5 + 8*c[1] + 2*c[2])) - 12*(-7 + 2*a*a + 3*a*b + 2*b*b)*c[3])*d[1]
	 + 6*(35 + 98*c[2] - a*b*(21 + 38*c[2]) + a*a*(-7 + 7*c[1] - 30*c[2] + 15*
	c[3]) - b*b*(7 + 7*c[1] + 30*c[2] + 15*c[3]))*d[2] + (9*(-28*c[1] + 12*a*
	b*c[1] + b*b*(-7 + 8*c[1] - 10*c[2]) + a*a*(7 + 8*c[1] + 10*c[2])) - 4*
	(-153 + 50*a*a + 53*a*b + 50*b*b)*c[3])*d[3]) + 231*c[0]*c[0]*(5*(-3 + 
	a*a + a*b + b*b)*d[0] + 10*d[2] - 6*a*b*d[2] + b*b*(5*d[1] - 2*d[2] - 3*
	d[3]) + a*a*(-5*d[1] - 2*d[2] + 3*d[3])) + 2*(22*(-3*(35*c[1]*c[1] + 7*
	c[2]*(5 + 7*c[2]) - 42*c[1]*c[3] + 51*c[3]*c[3])*d[0] - 6*(7*c[1]*(-5 + 
	2*c[2]) + 3*(7 + 10*c[2])*c[3])*d[1] - 2*(21*c[1]*c[1] - 3*c[2]*(49 + 27*
	c[2]) + 90*c[1]*c[3] - 59*c[3]*c[3])*d[2] - 2*(9*c[1]*(7 + 10*c[2]) - 
	(153 + 118*c[2])*c[3])*d[3]) + 2*a*b*(11*c[2]*c[2]*(57*d[0] - 94*d[2]) + 
	33*c[1]*c[1]*(7*d[0] - 2*d[2]) + c[3]*(583*c[3]*d[0] + 594*d[1] - 642*
	c[3]*d[2] - 1166*d[3]) - 22*c[1]*(3*(7 + 2*c[2])*d[1] + c[3]*(27*d[0] - 
	26*d[2]) - (27 + 26*c[2])*d[3]) + c[2]*(693*d[0] + 572*c[3]*d[1] - 1254*
	d[2] - 1284*c[3]*d[3])) + a*a*(c[3]*(11*(-63 + 100*c[3])*d[0] - 2*c[3]*
	(517*d[1] + 328*d[2] - 531*d[3]) + 22*(36*d[1] + 45*d[2] - 100*d[3])) + 
	66*c[1]*c[1]*(14*d[0] - 21*d[1] + 8*d[2] + 3*d[3]) + 22*c[2]*c[2]*(45*d[0]
	 - 33*d[1] - 34*d[2] + 7*d[3]) + 2*c[2]*(-33*(-7 + 15*c[3])*d[0] + 11*(21
	 + 64*c[3])*d[1] + 22*(-45 + 7*c[3])*d[2] + 495*d[3] - 656*c[3]*d[3]) + 
	11*c[1]*(-3*(-35 + 14*c[2] + 24*c[3])*d[0] + 2*(6*(-14 + 8*c[2] + 3*c[3])
	*d[1] + (21 - 66*c[2] + 64*c[3])*d[2] + 2*(18 + 32*c[2] - 47*c[3])*d[3])))
	 + b*b*(22*c[2]*(3*(7 + 15*c[3])*d[0] + (-21 + 64*c[3])*d[1] - 2*(45 + 7*
	c[3])*d[2]) + 22*c[2]*c[2]*(45*d[0] + 33*d[1] - 34*d[2] - 7*d[3]) + 66*
	c[1]*c[1]*(14*d[0] + 21*d[1] + 8*d[2] - 3*d[3]) - 2*c[2]*(495 + 656*c[3])*
	d[3] + 11*c[1]*(3*(-35 + 14*c[2] - 24*c[3])*d[0] + 2*(6*(-14 + 8*c[2] - 
	3*c[3])*d[1] + (-21 + 66*c[2] + 64*c[3])*d[2] + 2*(18 + 32*c[2] + 47*c[3])
	*d[3])) + c[3]*(11*(63 + 100*c[3])*d[0] + 2*c[3]*(517*d[1] - 328*d[2] - 
	531*d[3]) + 22*(36*d[1] - 5*(9*d[2] + 20*d[3])))))))/55440.;
		bf_mom[8] = -((a - b)*(11*(21*(5*(-3 + a*a + a*b + b*b)* (-4 + 
	c[0]*c[0]) - 10*(a - b)*(a + b)*c[0]*c[1] + 4*(-5 + 2*a*a + a*b + 2*b*b)*
	c[1]*c[1]) - 84*((-5 + a*a + 3*a*b + b*b)*c[0] + (a - b)*(a + b)*c[1])*
	c[2] + 12*(-49 + 15*a*a + 19*a*b + 15*b*b)* c[2]*c[2] + 18* (28*c[1] - 
	12*a*b*c[1] + a*a*(7*c[0] - 8*c[1] - 10*c[2]) + b*b*(-7*c[0] - 8*c[1] + 
	10*c[2]))*c[3] + 4*(-153 + 50*a*a + 53*a*b + 50*b*b)*c[3]*c[3])* d[0] - 
	4620*b*b*d[1] + 1155*b*b*c[0]*c[0]*d[1] - 9240*c[0]*c[1]*d[1] + 3696*b*b*
	c[0]*c[1]*d[1] + 2772*b*b*c[1]*c[1]*d[1] + 924*b*b*c[0]*c[2]*d[1] - 3696*
	c[1]*c[2]*d[1] + 2112*b*b*c[1]*c[2]*d[1] + 1452*b*b*c[2]*c[2]*d[1] + 5544*
	c[0]*c[3]*d[1] - 1584*b*b*c[0]*c[3]*d[1] - 792*b*b*c[1]*c[3]*d[1] - 7920*
	c[2]*c[3]*d[1] + 2816*b*b*c[2]*c[3]*d[1] + 2068*b*b*c[3]*c[3]*d[1] - 9240*
	d[2] + 1848*b*b*d[2] + 2310*c[0]*c[0]*d[2] - 462*b*b*c[0]*c[0]*d[2] + 924*
	b*b*c[0]*c[1]*d[2] - 1848*c[1]*c[1]*d[2] + 1056*b*b*c[1]*c[1]*d[2] - 
	12936*c[0]*c[2]*d[2] + 3960*b*b*c[0]*c[2]*d[2] + 2904*b*b*c[1]*c[2]*d[2] 
	+ 7128*c[2]*c[2]*d[2] - 1496*b*b*c[2]*c[2]*d[2] + 1980*b*b*c[0]*c[3]*d[2] 
	- 7920*c[1]*c[3]*d[2] + 2816*b*b*c[1]*c[3]*d[2] - 616*b*b*c[2]*c[3]*d[2] 
	+ 5192*c[3]*c[3]*d[2] - 1312*b*b*c[3]*c[3]*d[2] - (b*b*(99*(-28 + (c[0] + 
	2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2] + 308*c[2]*c[2] 
	- 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3]) - 88*(9*c[0]*
	(7*c[1] - 17*c[3]) + 2*c[2]*(-45*c[1] + 59*c[3])))*d[3] + a*a*(-11*(105*
	c[0]*c[0] + c[0]*(-336*c[1] + 84*c[2] + 144*c[3]) + 4*(-105 + 63*c[1]*c[1]
	 + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3] - 6*c[1]*(8*c[2] + 3*c[3])))
	*d[1] - 2*(231*c[0]*c[0] - 44*(21 + 12*c[1]*c[1] - 33*c[1]*c[2] - 17*c[2]*
	c[2]) - 44*(32*c[1] + 7*c[2])*c[3] + 656*c[3]*c[3] + 66*c[0]*(7*c[1] + 15*
	(-2*c[2] + c[3])))* d[2] + (693*c[0]*c[0] - 396*c[0]*(4*c[1] + 5*c[2]) + 
	44*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2])) + 8*(550*c[0] - 517*c[1] - 
	328*c[2])*c[3] + 2124*c[3]*c[3])* d[3]) + 2*a*b*(-693*c[0]*c[0]*d[2] + 
	44*c[0]*(57*c[2]*d[2] + 3*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-27*d[1] + 53*
	d[3])) - 4*(-693*d[2] + 33*c[1]*c[1]*d[2] + 517*c[2]*c[2]*d[2] + 321*c[3]*
	c[3]*d[2] + c[2]*c[3]*(-286*d[1] + 642*d[3]) - 22*c[1]*(13*c[3]*d[2] + 
	c[2]*(-3*d[1] + 13*d[3]))))))/27720.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
case 5:
switch (interface_type)
	{
	case 0:
 		memset(bf_mom, 0, sizeof(double)*npts);
		break;
	case 2:
	case 3:
		bf_mom[0] = -(pow(-1 + a,2)*(-26*c[0]*(33* (7*(5 + 6*c[2] + 2*a*
	(5 + 5*c[1] + 2*c[2] - 3*c[3])) + 2*(-3 + 10*a)*c[4])*d[0] + 22*(-42*c[1]*
	d[1] + 54*c[3]*d[1] + (-63 - 114*c[2] + 62*c[4])*d[2] + 54*c[1]*d[3] - 
	106*c[3]*d[3] + a* (-3*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])*d[1] + 
	2*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])*d[2] + (63 + 72*c[1] + 90*
	c[2] - 200*c[3] + 98*c[4])*d[3])) - 2*(11*(-9 - 62*c[2] + a*(30 + 78*c[1] 
	- 52*c[2] - 98*c[3])) + 2*(579 + 1126*a)*c[4])*d[4]) + 429*c[0]*c[0]*(35*
	(1 + 2*a)*d[0] + 6*(-7*d[2] + d[4]) - 2*a*(35*d[1] + 14*d[2] - 21*d[3] + 
	10*d[4])) + 4*(9009*c[2]*d[0] + 8151*c[2]*c[2]*d[0] + 7579*c[3]*c[3]*d[0] 
	- 1287*c[4]*d[0] - 8866*c[2]*c[4]*d[0] + 7527*c[4]*c[4]*d[0] + 7722*c[3]*
	d[1] + 7436*c[2]*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 16302*c[2]*d[2] - 13442*
	c[2]*c[2]*d[2] - 8346*c[3]*c[3]*d[2] + 8866*c[4]*d[2] + 17628*c[2]*c[4]*d[2]
	 - 8930*c[4]*c[4]*d[2] - 15158*c[3]*d[3] - 16692*c[2]*c[3]*d[3] - 6420*c[3]
	*c[4]*d[3] - 26*c[1]*(11*(21 + 6*c[2] + 22*c[4])*d[1] - (297 + 286*c[2] + 
	238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) + 2*(4407*c[2]*c[2] - 1605*c[3]*c[3] + c[2]*
	(4433 - 8930*c[4]) + c[4]*(-7527 + 967*c[4]))*d[4] + a*(-9009*c[3]*d[0] + 
	14300*c[3]*c[3]*d[0] + 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]
	*c[4]*d[0] + 10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 
	21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]
	*c[3]*d[2] + 7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 
	28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]
	*d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] 
	- 11*(21 + 64*c[3] - 10*c[4])*d[1] - 2*(-495 + 77*c[3] + 398*c[4])*d[2] + 
	(-495 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 
	24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 3478*c[4])*d[4] - 
	2*(3764*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))
	*d[4] + 26*c[2]*c[2]* (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) 
	- 13*c[1]*(33*(-35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 24*
	c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2]
	 - 2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] - 840*
	c[3] + 1094*c[4])*d[4]))))))/ 2882880.;
		bf_mom[1] = -((-1 + a)*(429*c[0]*c[0]* (35*(5 + a*(5 + 2*a))*d[0] 
	- 70*(-2 + a + a*a)*d[1] - 14*(7 + a*(11 + 2*a))*d[2] + 42*(-2 + a + a*a)*
	d[3] - 2*(17 + a + 10*a*a)*d[4]) - 26*c[0]*(33*(7*(25 - 20*c[1] + 14*c[2] 
	+ a*(25 + 10*c[1] + 22*c[2] + 2*a*(5 + 5*c[1] + 2*c[2] - 3*c[3]) - 6*c[3]) 
	+ 12*c[3]) + 2*(17 + a + 10*a*a)*c[4])*d[0] + 2*(-33*(-70 + 126*c[1] + 
	28*c[2] - 66*c[3] - 52*c[4] + a*(35 + 98*c[1] - 14*c[2] - 78*c[3] + 26*c[4]
	 + a*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])))*d[1] + 11*(3*(-49 - 28*
	c[1] - 158*c[2] - 60*c[3] + a*(-77 - 14*a + 14*c[1] + 14*a*c[1] - 174*c[2] 
	- 60*a*c[2] + 30*(1 + a)*c[3])) + 2*(83 + a*(119 + 26*a))*c[4])*d[2] + 11*
	(2*(-63 + 99*c[1] - 90*c[2] - 253*c[3] - 98*c[4]) + a*(63 + 63*a + 234*c[1] 
	+ 72*a*c[1] + 90*c[2] + 90*a*c[2] - 518*c[3] - 200*a*c[3] + 98*(1+a)*c[4]))*
	 d[3] - (11*(51 - 156*c[1] - 166*c[2] + 196*c[3] + a*(3 + 30*a + 78*c[1] + 
	78*a*c[1] - 238*c[2] - 52*a*c[2] - 98*(1 + a)*c[3])) + 2*(2831 + a*(2863 + 
	1126*a))*c[4])*d[4])) + 4*(18018*c[3]*d[0] - 9009*a*c[3]*d[0] - 9009*a*a*
	c[3]*d[0] + 36179*c[3]*c[3]*d[0] + 37037*a*c[3]*c[3]*d[0] + 14300*a*a*c[3]*
	c[3]*d[0] + 7293*c[4]*d[0] + 429*a*c[4]*d[0] + 4290*a*a*c[4]*d[0] + 28028*
	c[3]*c[4]*d[0] - 14014*a*c[3]*c[4]*d[0] - 14014*a*a*c[3]*c[4]*d[0] + 36803*
	c[4]*c[4]*d[0] + 37219*a*c[4]*c[4]*d[0] + 14638*a*a*c[4]*c[4]*d[0] + 28314*
	c[3]*d[1] + 33462*a*c[3]*d[1] + 10296*a*a*c[3]*d[1] + 26884*c[3]*c[3]*d[1] 
	- 13442*a*c[3]*c[3]*d[1] - 13442*a*a*c[3]*c[3]*d[1] + 22308*c[4]*d[1] - 
	11154*a*c[4]*d[1] - 11154*a*a*c[4]*d[1] + 49868*c[3]*c[4]*d[1] + 40404*a*
	c[3]*c[4]*d[1] + 21840*a*a*c[3]*c[4]*d[1] + 28444*c[4]*c[4]*d[1] - 14222*a*
	c[4]*c[4]*d[1] - 14222*a*a*c[4]*c[4]*d[1] - 25740*c[3]*d[2] + 12870*a*c[3]*
	d[2] + 12870*a*a*c[3]*d[2] - 25402*c[3]*c[3]*d[2] - 33566*a*c[3]*c[3]*d[2] 
	- 8528*a*a*c[3]*c[3]*d[2] + 23738*c[4]*d[2] + 34034*a*c[4]*d[2] + 7436*a*a*
	c[4]*d[2] + 6552*c[3]*c[4]*d[2] - 3276*a*c[3]*c[4]*d[2] - 3276*a*a*c[3]*
	c[4]*d[2] - 22842*c[4]*c[4]*d[2] - 33746*a*c[4]*c[4]*d[2] - 6956*a*a*c[4]*
	c[4]*d[2] - 72358*c[3]*d[3] - 74074*a*c[3]*d[3] - 28600*a*a*c[3]*d[3] - 
	27612*c[3]*c[3]*d[3] + 13806*a*c[3]*c[3]*d[3] + 13806*a*a*c[3]*c[3]*d[3] - 
	28028*c[4]*d[3] + 14014*a*c[4]*d[3] + 14014*a*a*c[4]*d[3] - 36532*c[3]*c[4]
	*d[3] - 34316*a*c[3]*c[4]*d[3] - 15056*a*a*c[3]*c[4]*d[3] - 20548*c[4]*c[4]
	*d[3] + 10274*a*c[4]*c[4]*d[3] + 10274*a*a*c[4]*c[4]*d[3] - 2*((9133 + a*
	(8579 + 3764*a))*c[3]*c[3] - 11*(-2 + a + a*a)*c[3]*(637 + 934*c[4]) + c[4]*
	(13*(2831 + a*(2863 + 1126*a)) + (5629 + a*(397 + 3298*a))*c[4]))*d[4] + 
	143*c[1]*c[1]*(21*(9 + a*(7 + 4*a))*d[0] - 126*(-2 + a + a*a)*d[1] + 6*(15 
	+ a*(5 + 8*a))*d[2] + 18*(-2 + a + a*a)*d[3] - 2*(67 + a*(61 + 28*a))*d[4]) 
	+ 13*c[2]*c[2]*(33*(79 + 87*a + 30*a*a)*d[0] - 726*(-2 + a + a*a)*d[1] - 22*
	(115 + a*(175 + 34*a))*d[2] + 154*(-2 + a + a*a)*d[3] + 2*(1135 + a*(1415 
	+ 398*a))*d[4]) + c[2]*(-143*(-21*(7 + a*(11 + 2*a)) + 90*(-2 + a + a*a)*
	c[3] + 2*(83 + a*(119 + 26*a))*c[4])*d[0] + 286*(-42 + 154*c[3] + 20*c[4] 
	+ a*(21 + 21*a + 142*c[3] + 64*a*c[3] - 10*(1 + a)*c[4]))* d[1] + 26*(-33*
	(79 + 87*a + 30*a*a) + 154*(-2 + a + a*a)*c[3] + 2*(1135 + a*(1415 + 398*a))
	*c[4])*d[2] - 26*(990 + 1954*c[3] - 252*c[4] + a*(-495 - 495*a + 2582*c[3] 
	+ 656*a*c[3] + 126*(1 + a)*c[4]))*d[3] - 2*(-143*(83 + a*(119 + 26*a)) + 
	1638*(-2 + a + a*a)*c[3] + 94*(243 + a*(359 + 74*a))*c[4])*d[4]) - 13*c[1]*
	(33*(70 - 28*c[2] + 66*c[3] + 52*c[4] + a*(-35 - 35*a + 14*c[2] + 14*a*c[2]
	 + 78*c[3] + 24*a*c[3] - 26*(1 + a)*c[4]))*d[0] + 66*(63*d[1] + 14*d[2] - 
	33*d[3] - 26*d[4]) + 2*(a*(-11*(-147 + 30*c[2] + 18*c[3] - 122*c[4])*d[1] 
	+ 11*(-21 + 66*c[2] - 142*c[3] + 10*c[4])*d[2] - (1287 + 1562*c[2] - 1034*
	c[3] + 1554*c[4])*d[3] + (429 + 110*c[2] - 1554*c[3] + 1094*c[4])*d[4]) + 
	a*a*(-22*(-42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 
	64*c[3] + 10*c[4])*d[2] - 2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + 
	(429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4]) - 2*(11*c[2]*(45*d[1] + 66*
	d[2] + 77*d[3] + 10*d[4]) + c[3]*(-198*d[1] + 847*d[2] + 1034*d[3] + 959*
	d[4]) + c[4]*(-737*d[1] + 110*d[2] + 959*d[3] + 1094*d[4])))))))/ 2882880.;
		bf_mom[2] = -((-1 + a)*(429*c[0]*c[0]* (35*(5 + a*(5 + 2*a))*d[0] 
	- 70*(-2 + a + a*a)*d[1] - 14*(7 + a*(11 + 2*a))*d[2] + 42*(-2 + a + a*a)
	*d[3] - 2*(17 + a + 10*a*a)*d[4]) - 26*c[0]*(33*(7*(-25 - 20*c[1] + 14*c[2]
	 + a*(-25 + 10*c[1] + 22*c[2] + 2*a*(-5 + 5*c[1] + 2*c[2] - 3*c[3])-6*c[3])
	 + 12*c[3]) + 2*(17 + a + 10*a*a)*c[4])*d[0] + 2*(-33*(70 + 126*c[1] + 
	28*c[2] - 66*c[3] - 52*c[4] + a*(-35 + 98*c[1] - 14*c[2] - 78*c[3] + 26*c[4]
	 + a*(-35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])))*d[1] + 11*(3*(49 - 
	28*c[1] - 158*c[2] - 60*c[3] + a*(77 + 14*a + 14*c[1] + 14*a*c[1] - 174*c[2]
	 - 60*a*c[2] + 30*(1 + a)*c[3])) + 2*(83 + a*(119 + 26*a))*c[4])*d[2] + 
	11*(2*(63 + 99*c[1] - 90*c[2] - 253*c[3] - 98*c[4]) + a*(-63 - 63*a + 234*
	c[1] + 72*a*c[1] + 90*c[2] + 90*a*c[2] - 518*c[3] - 200*a*c[3] + 98*(1 + a)
	*c[4]))* d[3] - (11*(-51 - 156*c[1] - 166*c[2] + 196*c[3] + a*(-3 - 30*a + 
	78*c[1] + 78*a*c[1] - 238*c[2] - 52*a*c[2] - 98*(1 + a)*c[3])) + 2*(2831 + 
	a*(2863 + 1126*a))*c[4])*d[4])) + 4*(-18018*c[3]*d[0] + 9009*a*c[3]*d[0] + 
	9009*a*a*c[3]*d[0] + 36179*c[3]*c[3]*d[0] + 37037*a*c[3]*c[3]*d[0] + 14300*
	a*a*c[3]*c[3]*d[0] - 7293*c[4]*d[0] - 429*a*c[4]*d[0] - 4290*a*a*c[4]*d[0] 
	+ 28028*c[3]*c[4]*d[0] - 14014*a*c[3]*c[4]*d[0] - 14014*a*a*c[3]*c[4]*d[0] 
	+ 36803*c[4]*c[4]*d[0] + 37219*a*c[4]*c[4]*d[0] + 14638*a*a*c[4]*c[4]*d[0] 
	- 28314*c[3]*d[1] - 33462*a*c[3]*d[1] - 10296*a*a*c[3]*d[1] + 26884*c[3]*
	c[3]*d[1] - 13442*a*c[3]*c[3]*d[1] - 13442*a*a*c[3]*c[3]*d[1] - 22308*c[4]*
	d[1] + 11154*a*c[4]*d[1] + 11154*a*a*c[4]*d[1] + 49868*c[3]*c[4]*d[1] + 
	40404*a*c[3]*c[4]*d[1] + 21840*a*a*c[3]*c[4]*d[1] + 28444*c[4]*c[4]*d[1] - 
	14222*a*c[4]*c[4]*d[1] - 14222*a*a*c[4]*c[4]*d[1] + 25740*c[3]*d[2] - 12870
	*a*c[3]*d[2] - 12870*a*a*c[3]*d[2] - 25402*c[3]*c[3]*d[2] - 33566*a*c[3]*
	c[3]*d[2] - 8528*a*a*c[3]*c[3]*d[2] - 23738*c[4]*d[2] - 34034*a*c[4]*d[2] 
	- 7436*a*a*c[4]*d[2] + 6552*c[3]*c[4]*d[2] - 3276*a*c[3]*c[4]*d[2] - 3276*
	a*a*c[3]*c[4]*d[2] - 22842*c[4]*c[4]*d[2] - 33746*a*c[4]*c[4]*d[2] - 6956*a
	*a*c[4]*c[4]*d[2] + 72358*c[3]*d[3] + 74074*a*c[3]*d[3] + 28600*a*a*c[3]*
	d[3] - 27612*c[3]*c[3]*d[3] + 13806*a*c[3]*c[3]*d[3] + 13806*a*a*c[3]*c[3]*
	d[3] + 28028*c[4]*d[3] - 14014*a*c[4]*d[3] - 14014*a*a*c[4]*d[3] - 36532*
	c[3]*c[4]*d[3] - 34316*a*c[3]*c[4]*d[3] - 15056*a*a*c[3]*c[4]*d[3] - 20548*
	c[4]*c[4]*d[3] + 10274*a*c[4]*c[4]*d[3] + 10274*a*a*c[4]*c[4]*d[3] - 2*(
	(9133 + a*(8579 + 3764*a))*c[3]*c[3] + 11*(-2 + a + a*a)*c[3]*(637 - 934*
	c[4]) + c[4]*(-13*(2831 + a*(2863 + 1126*a)) + (5629 + a*(397 + 3298*a))*
	c[4]))*d[4] + 143*c[1]*c[1]*(21*(9 + a*(7 + 4*a))*d[0] - 126*(-2 + a + a*a)
	*d[1] + 6*(15 + a*(5 + 8*a))*d[2] + 18*(-2 + a + a*a)*d[3] - 2*(67 + a*(61 
	+ 28*a))*d[4]) + 13*c[2]*c[2]*(33*(79 + 87*a + 30*a*a)*d[0] - 726*(-2 + a 
	+ a*a)*d[1] - 22*(115 + a*(175 + 34*a))*d[2] + 154*(-2 + a + a*a)*d[3] + 2*
	(1135 + a*(1415 + 398*a))*d[4]) - c[2]*(143*(21*(7 + a*(11 + 2*a)) + 90*
	(-2 + a + a*a)*c[3] + 2*(83 + a*(119 + 26*a))*c[4])*d[0] - 286*(42 + 154*
	c[3] + 20*c[4] + a*(-21 - 21*a + 142*c[3] + 64*a*c[3] - 10*(1 + a)*c[4]))*
	 d[1] - 26*(33*(79 + 87*a + 30*a*a) + 154*(-2 + a + a*a)*c[3] + 2*(1135 + 
	a*(1415 + 398*a))*c[4])*d[2] + 26*(-990 + 1954*c[3] - 252*c[4] + a*(495 + 
	495*a + 2582*c[3] + 656*a*c[3] + 126*(1 + a)*c[4]))*d[3] + 2*(143*(83 + a*
	(119 + 26*a)) + 1638*(-2 + a + a*a)*c[3] + 94*(243 + a*(359 + 74*a))*c[4])
	*d[4]) - 13*c[1]*(33*(-70 - 28*c[2] + 66*c[3] + 52*c[4] + a*(35 + 35*a + 
	14*c[2] + 14*a*c[2] + 78*c[3] + 24*a*c[3] - 26*(1 + a)*c[4]))*d[0] - 66*
	(63*d[1] + 14*d[2] - 33*d[3] - 26*d[4]) + 2*(a*(-11*(147 + 30*c[2] + 18*
	c[3] - 122*c[4])*d[1] + 11*(21 + 66*c[2] - 142*c[3] + 10*c[4])*d[2] - 
	(-1287 + 1562*c[2] - 1034*c[3] + 1554*c[4])*d[3] + (-429 + 110*c[2] - 1554*
	c[3] + 1094*c[4])*d[4]) + a*a*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] 
	+ 11*(21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3]
	 + 420*c[4])*d[3] + (-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4]) - 2*
	(11*c[2]*(45*d[1] + 66*d[2] + 77*d[3] + 10*d[4]) + c[3]*(-198*d[1] + 847*
	d[2] + 1034*d[3] + 959*d[4]) + c[4]*(-737*d[1] + 110*d[2] + 959*d[3] + 
	1094*d[4])))))))/ 2882880.;
		bf_mom[3] = -(pow(-1 + a,2)*(-26*c[0]*(33* (7*(-5 + 6*c[2] + 2*a*
	(-5 + 5*c[1] + 2*c[2] - 3*c[3])) + 2*(-3 + 10*a)*c[4])*d[0] + 22*(-42*c[1]*
	d[1] + 54*c[3]*d[1] + (63 - 114*c[2] + 62*c[4])*d[2] + 54*c[1]*d[3] - 106*
	c[3]*d[3] + a* ((105 - 168*c[1] + 42*c[2] + 72*c[3] - 78*c[4])*d[1] + 2*
	(21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])*d[2] + (-63 + 72*c[1] + 90*
	c[2] - 200*c[3] + 98*c[4])*d[3])) - 2*(11*(9 - 62*c[2] + a*(-30 + 78*c[1] 
	- 52*c[2] - 98*c[3])) + 2*(579 + 1126*a)*c[4])*d[4]) + 429*c[0]*c[0]*(35*
	(1 + 2*a)*d[0] + 6*(-7*d[2] + d[4]) - 2*a*(35*d[1] + 14*d[2] - 21*d[3] + 
	10*d[4])) + 4*(-9009*c[2]*d[0] + 8151*c[2]*c[2]*d[0] + 7579*c[3]*c[3]*d[0] 
	+ 1287*c[4]*d[0] - 8866*c[2]*c[4]*d[0] + 7527*c[4]*c[4]*d[0] - 7722*c[3]*
	d[1] + 7436*c[2]*c[3]*d[1] + 6188*c[3]*c[4]*d[1] + 16302*c[2]*d[2] - 13442*
	c[2]*c[2]*d[2] - 8346*c[3]*c[3]*d[2] - 8866*c[4]*d[2] + 17628*c[2]*c[4]*d[2]
	 - 8930*c[4]*c[4]*d[2] + 15158*c[3]*d[3] - 16692*c[2]*c[3]*d[3] - 6420*c[3]*
	c[4]*d[3] - 26*c[1]*(11*(-21 + 6*c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 
	238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) + 2*(4407*c[2]*c[2] - 1605*c[3]*c[3] + c[4]*
	(7527 + 967*c[4]) - c[2]*(4433 + 8930*c[4]))*d[4] + a*(9009*c[3]*d[0] + 
	14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*
	c[4]*d[0] - 10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 
	21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*
	c[3]*d[2] - 7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 
	28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*
	d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(21 + 45*c[3] + 26*c[4])*d[0] + 
	(231 - 704*c[3] + 110*c[4])*d[1] - 2*(495 + 77*c[3] + 398*c[4])*d[2] + 
	(495 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*
	d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(1859 + 819*c[3] + 3478*c[4])*d[4] - 2*
	(c[3]*(7007 + 3764*c[3]) - 2*(7319 + 5137*c[3])*c[4] + 3298*c[4]*c[4])*d[4]
	 + 26*c[2]*c[2]* (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 
	13*c[1]*(33*(35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(42 + 24*c[2] 
	+ 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*
	(-198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] - 840*c[3] 
	+ 1094*c[4])*d[4]))))))/ 2882880.;
		bf_mom[4] = (pow(-1 + a,2)*(429*c[0]*c[0]* (35*(2 + a)*d[0] - 35*
	(1 + a)*d[1] - 14*(4 + a)*d[2] + 21*(1 + a)*d[3] - 2*(2 + 5*a)*d[4]) - 26*
	c[0]*(33*(7*(10 + 5*a + 5*c[1] + 5*a*c[1] + 8*c[2] + 2*a*c[2] - 3*(1 + a)*
	c[3]) + 2*(2 + 5*a)*c[4])*d[0] - 33*(35 + 35*a + 84*c[1] + 56*a*c[1] - 14*
	c[2] - 14*a*c[2] - 60*c[3] - 24*a*c[3] + 26*(1 + a)*c[4])*d[1] + 22*(3*
	(-28 - 7*a + 7*c[1] + 7*a*c[1] - 68*c[2] - 30*a*c[2] + 15*(1 + a)*c[3]) + 
	2*(44 + 13*a)*c[4])*d[2] + 11*(63 + 63*a + 180*c[1] + 72*a*c[1] + 90*c[2] 
	+ 90*a*c[2] - 412*c[3] - 200*a*c[3] + 98*(1 + a)*c[4])*d[3] - 2*(11*(6 + 
	15*a + 39*c[1] + 39*a*c[1] - 88*c[2] - 26*a*c[2] - 49*(1 + a)*c[3]) + 2*
	(1142 + 563*a)*c[4])*d[4]) + 2*(-9009*c[3]*d[0] - 9009*a*c[3]*d[0] + 29458*
	c[3]*c[3]*d[0] + 14300*a*c[3]*c[3]*d[0] + 1716*c[4]*d[0] + 4290*a*c[4]*d[0]
	 - 14014*c[3]*c[4]*d[0] - 14014*a*c[3]*c[4]*d[0] + 29692*c[4]*c[4]*d[0] + 
	14638*a*c[4]*c[4]*d[0] + 25740*c[3]*d[1] + 10296*a*c[3]*d[1] - 13442*c[3]*
	c[3]*d[1] - 13442*a*c[3]*c[3]*d[1] - 11154*c[4]*d[1] - 11154*a*c[4]*d[1] + 
	34216*c[3]*c[4]*d[1] + 21840*a*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 
	14222*a*c[4]*c[4]*d[1] + 12870*c[3]*d[2] + 12870*a*c[3]*d[2] - 25220*c[3]*
	c[3]*d[2] - 8528*a*c[3]*c[3]*d[2] + 25168*c[4]*d[2] + 7436*a*c[4]*d[2] - 
	3276*c[3]*c[4]*d[2] - 3276*a*c[3]*c[4]*d[2] - 24816*c[4]*c[4]*d[2] - 6956*a
	*c[4]*c[4]*d[2] - 58916*c[3]*d[3] - 28600*a*c[3]*d[3] + 13806*c[3]*c[3]*
	d[3] + 13806*a*c[3]*c[3]*d[3] + 14014*c[4]*d[3] + 14014*a*c[4]*d[3] - 
	27896*c[3]*c[4]*d[3] - 15056*a*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 
	10274*a*c[4]*c[4]*d[3] + 26*c[2]*(-11*(-84 + 45*c[3] + 88*c[4] + a*(-21 + 
	45*c[3] + 26*c[4]))*d[0] + 11*(21 + 21*a + 116*c[3] + 64*a*c[3] - 10*(1+a)
	*c[4])*d[1] + 2*(-33*(34 + 15*a) + 77*(1 + a)*c[3] + 2*(538 + 199*a)*c[4])*
	 d[2] - (-495 - 495*a + 1940*c[3] + 656*a*c[3] + 126*(1 + a)*c[4])*d[3]) - 
	4*c[2]*(-143*(44 + 13*a) + 819*(1 + a)*c[3] + 94*(132 + 37*a)*c[4])*d[4] - 
	2*((6974 + 3764*a)*c[3]*c[3] - 11*(1 + a)*c[3]*(637 + 934*c[4]) + 2*c[4]*
	(14846 + 682*c[4] + a*(7319 + 1649*c[4])))*d[4] + 286*c[1]*c[1]*(21*(3+2*a)
	*d[0] - 63*(1 + a)*d[1] + 18*d[2] + 24*a*d[2] + 9*d[3] + 9*a*d[3] - 50*d[4]
	 - 28*a*d[4]) + 26*c[2]*c[2]*(33*(34 + 15*a)*d[0] - 363*(1 + a)*d[1] - 22*
	(64 + 17*a)*d[2] + 77*(1 + a)*d[3] + 2*(538 + 199*a)*d[4]) - 13*c[1]*(33*
	(-35 - 35*a + 14*c[2] + 14*a*c[2] + 60*c[3] + 24*a*c[3] - 26*(1 + a)*c[4])
	*d[0] - 44*(9*(-7 + 2*c[2] + c[3]) + a*(-42 + 24*c[2] + 9*c[3] - 28*c[4]) 
	- 50*c[4])*d[1] + 22*(-21 - 21*a + 66*c[2] + 66*a*c[2] - 116*c[3] - 64*a*
	c[3] + 10*(1 + a)*c[4])*d[2] - 4*(11*(45 + 18*a + 58*c[2] + 32*a*c[2] - 47*
	(1 + a)*c[3]) + 14*(47 + 30*a)*c[4])*d[3] + 2*(429 + 429*a + 110*c[2] + 
	110*a*c[2] - 1316*c[3] - 840*a*c[3] + 1094*(1 + a)*c[4])*d[4]))))/720720.;
		bf_mom[5] = ((-1 + a)*(13*(5775*(-4 + c[0]*c[0]) + 44*(3*(7*c[1]*
	(10*c[0] + 9*c[1]) + 7*(-7*c[0] + 4*c[1])*c[2] + 79*c[2]*c[2]) - 18*(7*c[0]
	 + 11*c[1] - 10*c[2])*c[3] + 253*c[3]*c[3]) - 44*(51*c[0] + 2*(78*c[1] + 
	83*c[2] - 98*c[3]))*c[4] + 11324*c[4]*c[4] + 2*a*a*(33*(7* (-20 + 5*c[0]*
	c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) 
	+ 198*(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 
	39*c[1] + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4]) + a*(5775*c[0]*c[0] + 
	44*(-525 + 147*c[1]*c[1] + 261*c[2]*c[2] - 90*c[2]*c[3] + 259*c[3]*c[3] - 
	6*c[1]*(7*c[2] + 39*c[3])) + 88*(39*c[1] - 7*(17*c[2] + 7*c[3]))*c[4] + 
	11452*c[4]*c[4] - 132*c[0]*(35*c[1] + 77*c[2] - 21*c[3] + c[4])))*d[0] - 
	3432*(70*d[1] - 49*d[2] - 42*d[3] - 17*d[4]) + 2*(429*c[0]*c[0]*(70*d[1] - 
	49*d[2] - 42*d[3] - 17*d[4]) + 52*c[0]*(-11*(c[3]*(99*d[1] - 90*d[2] - 
	253*d[3]) + c[4]*(78*d[1] + 83*d[2] - 98*d[3])) + 11*c[2]*(42*d[1] + 237*
	d[2] + 90*d[3] - 83*d[4]) + 33*c[1]*(63*d[1] + 14*d[2] - 33*d[3] - 26*d[4])
	 + (1078*c[3] + 2831*c[4])*d[4]) - a*(13*(1155*c[0]*c[0] + 44*(-105 + 63*
	c[1]*c[1] + 33*c[2]*c[2] - 142*c[2]*c[3] + 47*c[3]*c[3] - 6*c[1]*(5*c[2] + 
	3*c[3])) + 8*(671*c[1] + 55*c[2] - 777*c[3])*c[4] + 2188*c[4]*c[4] - 132*
	c[0]*(49*c[1] - 7*c[2] - 39*c[3] + 13*c[4]))*d[1] + (33033*c[0]*c[0] + 52*
	(-33*(77 + 5*c[1]*c[1]) + 726*c[1]*c[2] + 1925*c[2]*c[2] - 22*(71*c[1] + 
	7*c[2])*c[3] + 1291*c[3]*c[3]) + 104*(55*c[1] - 1415*c[2] + 63*c[3])*c[4] 
	+ 67492*c[4]*c[4] + 572*c[0]*(21*c[1] - 261*c[2] + 45*c[3] + 119*c[4]))*d[2]
	 - (9009*c[0]*c[0] + 52*(99*c[1]*c[1] + 77*(-9 + c[2]*c[2]) + 22*c[1]*(71*
	c[2] - 47*c[3]) - 2582*c[2]*c[3] + 531*c[3]*c[3]) + 8*(10101*c[1] - 819*c[2]
	 - 8579*c[3])*c[4] + 20548*c[4]*c[4] - 572*c[0]*(117*c[1] + 45*c[2] - 259*
	c[3] + 49*c[4]))*d[3] + (429*c[0]*c[0] - 52*c[0]*(429*c[1] - 7*(187*c[2] +
	 77*c[3] - 409*c[4])) + 4*(-429 + 8723*c[1]*c[1] - 18395*c[2]*c[2] + 1638*
	c[2]*c[3] + 8579*c[3]*c[3] + 33746*c[2]*c[4] - 10274*c[3]*c[4] + 397*c[4]*
	c[4] + 26*c[1]*(55*c[2] - 777*c[3] + 547*c[4])))*d[4]) - a*a*(13*(1155*c[0]
	*c[0] + 44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3]
	 - 6*c[1]*(8*c[2] + 3*c[3])) + 8*(308*c[1] + 55*c[2] - 420*c[3])*c[4] + 
	2188*c[4]*c[4] - 132*c[0]*(28*c[1] - 7*c[2] - 12*c[3] + 13*c[4]))*d[1] + 
	2*(143*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 30*c[2]) + 
	132*c[1]*c[2] + 68*c[2]*c[2]) + 286*(45*c[0] - 64*c[1] - 14*c[2])*c[3] + 
	8528*c[3]*c[3] + 52*(143*c[0] + 55*c[1] - 398*c[2] + 63*c[3])*c[4] + 6956*
	c[4]*c[4])*d[2] - (9009*c[0]*c[0] + 52*(11*(-63 + (9*c[1] + c[2])*(c[1] + 
	7*c[2])) - 2*(517*c[1] + 328*c[2])*c[3] + 531*c[3]*c[3]) + 8*(5460*c[1] - 
	819*c[2] - 3764*c[3])*c[4] + 20548*c[4]*c[4] - 572*c[0]*(36*c[1] + 45*c[2] 
	- 100*c[3] + 49*c[4]))*d[3] + 2*(2145*c[0]*c[0] - 26*c[0]*(429*c[1] - 286*
	c[2] - 539*c[3] + 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] 
	+ 819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] - 5137*c[3]*c[4] + 1649*
	c[4]*c[4] + 13*c[1]*(55*c[2] - 420*c[3] + 547*c[4])))*d[4]) + 4*(13442*c[3]
	*c[3]*d[1] + 24934*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] - 12701*c[3]*c[3]
	*d[2] + 3276*c[3]*c[4]*d[2] - 11421*c[4]*c[4]*d[2] - 13806*c[3]*c[3]*d[3] 
	- 18266*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 143*c[1]*c[1]*(126*d[1] + 
	45*d[2] - 18*d[3] - 67*d[4]) - (9133*c[3]*c[3] + 20548*c[3]*c[4] + 5629*
	c[4]*c[4])* d[4] + 13*c[2]*c[2]* (726*d[1] - 1265*d[2] - 154*d[3] + 1135*
	d[4]) + 2*c[2]*(c[4]*(1430*d[1] + 14755*d[2] + 1638*d[3] - 11421*d[4]) + 
	13*c[3]*(847*d[1] - 154*d[2] - 977*d[3] + 126*d[4])) + 26*c[1]*(11*c[2]*
	(45*d[1] + 66*d[2] + 77*d[3] + 10*d[4]) + c[3]*(-198*d[1] + 847*d[2] + 
	1034*d[3] + 959*d[4]) + c[4]*(-737*d[1] + 110*d[2] + 959*d[3] + 1094*
	d[4]))))))/ 1441440.;
		bf_mom[6] = (pow(-1 + a,2)*(429*c[0]*c[0]* (35*(2 + a)*d[0] - 35*
	(1 + a)*d[1] - 14*(4 + a)*d[2] + 21*(1 + a)*d[3] - 2*(2 + 5*a)*d[4]) - 
	26*c[0]*(33*(7*(-10 - 5*a + 5*c[1] + 5*a*c[1] + 8*c[2] + 2*a*c[2] - 3*(1+a)
	*c[3]) + 2*(2 + 5*a)*c[4])*d[0] - 33*(-35 - 35*a + 84*c[1] + 56*a*c[1] - 
	14*c[2] - 14*a*c[2] - 60*c[3] - 24*a*c[3] + 26*(1 + a)*c[4])*d[1] + 22*(3*
	(28 + 7*a + 7*c[1] + 7*a*c[1] - 68*c[2] - 30*a*c[2] + 15*(1 + a)*c[3]) + 2*
	(44 + 13*a)*c[4])*d[2] + 11*(-63 - 63*a + 180*c[1] + 72*a*c[1] + 90*c[2] + 
	90*a*c[2] - 412*c[3] - 200*a*c[3] + 98*(1 + a)*c[4])*d[3] - 2*(11*(-6 - 
	15*a + 39*c[1] + 39*a*c[1] - 88*c[2] - 26*a*c[2] - 49*(1 + a)*c[3]) + 2*
	(1142 + 563*a)*c[4])*d[4]) + 2*(9009*c[3]*d[0] + 9009*a*c[3]*d[0] + 29458*
	c[3]*c[3]*d[0] + 14300*a*c[3]*c[3]*d[0] - 1716*c[4]*d[0] - 4290*a*c[4]*d[0]
	 - 14014*c[3]*c[4]*d[0] - 14014*a*c[3]*c[4]*d[0] + 29692*c[4]*c[4]*d[0] + 
	14638*a*c[4]*c[4]*d[0] - 25740*c[3]*d[1] - 10296*a*c[3]*d[1] - 13442*c[3]*
	c[3]*d[1] - 13442*a*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 11154*a*c[4]*d[1] + 
	34216*c[3]*c[4]*d[1] + 21840*a*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 
	14222*a*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 12870*a*c[3]*d[2] - 25220*c[3]*
	c[3]*d[2] - 8528*a*c[3]*c[3]*d[2] - 25168*c[4]*d[2] - 7436*a*c[4]*d[2] - 
	3276*c[3]*c[4]*d[2] - 3276*a*c[3]*c[4]*d[2] - 24816*c[4]*c[4]*d[2] - 6956*
	a*c[4]*c[4]*d[2] + 58916*c[3]*d[3] + 28600*a*c[3]*d[3] + 13806*c[3]*c[3]*
	d[3] + 13806*a*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 14014*a*c[4]*d[3] - 
	27896*c[3]*c[4]*d[3] - 15056*a*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 
	10274*a*c[4]*c[4]*d[3] + 26*c[2]*(-11*(84 + 45*c[3] + 88*c[4] + a*(21 + 
	45*c[3] + 26*c[4]))*d[0] + 11*(-21 - 21*a + 116*c[3] + 64*a*c[3] - 10*(1+a)
	*c[4])* d[1] + 2*(33*(34 + 15*a) + 77*(1 + a)*c[3] + 2*(538 + 199*a)*c[4])
	*d[2] - (495 + 495*a + 1940*c[3] + 656*a*c[3] + 126*(1 + a)*c[4])*d[3]) - 
	4*c[2]*(143*(44 + 13*a) + 819*(1 + a)*c[3] + 94*(132 + 37*a)*c[4])*d[4] - 
	2*((6974 + 3764*a)*c[3]*c[3] + 11*(1 + a)*c[3]*(637 - 934*c[4]) + 2*c[4]*
	(-14846 + 682*c[4] + a*(-7319 + 1649*c[4])))*d[4] + 286*c[1]*c[1]*(21*(3 + 
	2*a)*d[0] - 63*(1 + a)*d[1] + 18*d[2] + 24*a*d[2] + 9*d[3] + 9*a*d[3] - 50*
	d[4] - 28*a*d[4]) + 26*c[2]*c[2]*(33*(34 + 15*a)*d[0] - 363*(1 + a)*d[1] - 
	22*(64 + 17*a)*d[2] + 77*(1 + a)*d[3] + 2*(538 + 199*a)*d[4]) - 13*c[1]*
	(33*(35 + 35*a + 14*c[2] + 14*a*c[2] + 60*c[3] + 24*a*c[3] - 26*(1+a)*c[4])
	*d[0] + 2*(-22*(9*(7 + 2*c[2] + c[3]) + a*(42 + 24*c[2] + 9*c[3] - 28*c[4])
	 - 50*c[4])*d[1] + 11*(21 + 21*a + 66*c[2] + 66*a*c[2] - 116*c[3] - 64*a*
	c[3] + 10*(1 + a)*c[4])*d[2] - 2*(11*(-45 - 18*a + 58*c[2] + 32*a*c[2] - 
	47*(1 + a)*c[3]) + 14*(47 + 30*a)*c[4])*d[3] + (-429 - 429*a + 110*c[2] + 
	110*a*c[2] - 1316*c[3] - 840*a*c[3] + 1094*(1 + a)*c[4])*d[4])))))/720720.;
		bf_mom[7] = (pow(-1 + a,2)*(13*(33*(7*(5*(1 + 2*a)*(-4 + c[0]*c[0])
	 - 20*a*c[0]*c[1] + 4*(1 + 4*a)*c[1]*c[1]) - 28*((3 + 2*a)*c[0] + 2*a*c[1])
	*c[2] + 4*(19 + 30*a)*c[2]*c[2]) - 396*(6*c[1] + a*(-7*c[0] + 8*c[1] + 10*
	c[2]))*c[3] + 44*(53 + 100*a)*c[3]*c[3] - 44*((-9 + 30*a)*c[0] + 62*c[2] + 
	a*(-78*c[1] + 52*c[2] + 98*c[3]))*c[4] + 4*(579 + 1126*a)*c[4]*c[4])*d[0] + 
	2*(-104*(11*c[1]*(3*c[2] + 11*c[4]) - c[3]*(143*c[2] + 119*c[4]))* d[1] - 
	4*(429*c[1]*c[1] + 6721*c[2]*c[2] - 3718*c[1]*c[3] + 4173*c[3]*c[3] - 8814*
	c[2]*c[4] + 4465*c[4]*c[4])*d[2] + 8*(-321*c[3]*(13*c[2] + 5*c[4]) + 13*
	c[1]*(143*c[2] + 119*c[4]))* d[3] + 5148*(7*d[2] - d[4]) - 1287*c[0]*c[0]*
	(7*d[2] - d[4]) - 4*(1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 
	1605*c[3]*c[3] + 8930*c[2]*c[4] - 967*c[4]*c[4])*d[4] + 52*c[0]*(627*c[2]*
	d[2] - 341*c[4]*d[2] + 33*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-297*d[1] + 583*
	d[3]) - 341*c[2]*d[4] + 579*c[4]*d[4]) - a*(13*(1155*c[0]*c[0] + 44*(-105 
	+ 63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3] - 6*c[1]*(8*
	c[2] + 3*c[3])) + 8*(308*c[1] + 55*c[2] - 420*c[3])*c[4] + 2188*c[4]*c[4] 
	- 132*c[0]*(28*c[1] - 7*c[2] - 12*c[3] + 13*c[4]))*d[1] + 2*(143*(-84 + 21*
	c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 30*c[2]) + 132*c[1]*c[2] + 68*
	c[2]*c[2]) + 286*(45*c[0] - 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*
	(143*c[0] + 55*c[1] - 398*c[2] + 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - 
	(9009*c[0]*c[0] + 52*(11*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2])) - 2*(517*
	c[1] + 328*c[2])*c[3] + 531*c[3]*c[3]) + 8*(5460*c[1] - 819*c[2]-3764*c[3])
	*c[4] + 20548*c[4]*c[4] - 572*c[0]*(36*c[1] + 45*c[2] - 100*c[3]+49*c[4]))
	*d[3] + 2*(2145*c[0]*c[0] - 26*c[0]*(429*c[1] - 286*c[2] - 539*c[3] + 
	1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] + 819*c[2]*c[3] + 
	1882*c[3]*c[3] + 3478*c[2]*c[4] - 5137*c[3]*c[4] + 1649*c[4]*c[4] + 13*c[1]*
	(55*c[2] - 420*c[3] + 547*c[4])))*d[4]))))/ 1441440.;
		bf_mom[8] = -(pow(-1 + a,2)*(13*(33*(35*(2 + a)*(-4 + c[0]*c[0]) - 
	70*(1 + a)*c[0]*c[1] + 28*(3 + 2*a)*c[1]*c[1] - 28*((4 + a)*c[0] + (1 + a)*
	c[1])*c[2] + 4*(34 + 15*a)*c[2]*c[2]) + 198*(7*(1 + a)*c[0] - 4*(5 + 2*a)*
	c[1] - 10*(1 + a)*c[2])* c[3] + 44*(103 + 50*a)*c[3]*c[3] - 44*(3*(2 + 5*a)
	*c[0] - 39*(1 + a)*c[1] + 2*(44 + 13*a)*c[2] + 49*(1 + a)*c[3])*c[4] + 4*
	(1142 + 563*a)*c[4]*c[4])*d[0] - 15015*c[0]*c[0]*d[1] + 72072*c[0]*c[1]*
	d[1] - 36036*c[1]*c[1]*d[1] - 12012*c[0]*c[2]*d[1] + 20592*c[1]*c[2]*d[1] 
	- 18876*c[2]*c[2]*d[1] - 51480*c[0]*c[3]*d[1] + 10296*c[1]*c[3]*d[1] + 
	66352*c[2]*c[3]*d[1] - 26884*c[3]*c[3]*d[1] + 22308*c[0]*c[4]*d[1] - 57200*
	c[1]*c[4]*d[1] - 5720*c[2]*c[4]*d[1] + 68432*c[3]*c[4]*d[1] - 28444*c[4]*
	c[4]*d[1] - 24024*c[0]*c[0]*d[2] - 12012*c[0]*c[1]*d[2] + 10296*c[1]*c[1]*
	d[2] + 116688*c[0]*c[2]*d[2] - 37752*c[1]*c[2]*d[2] - 73216*c[2]*c[2]*d[2] 
	- 25740*c[0]*c[3]*d[2] + 66352*c[1]*c[3]*d[2] + 8008*c[2]*c[3]*d[2] - 50440*
	c[3]*c[3]*d[2] - 50336*c[0]*c[4]*d[2] - 5720*c[1]*c[4]*d[2] + 111904*c[2]*
	c[4]*d[2] - 6552*c[3]*c[4]*d[2] - 49632*c[4]*c[4]*d[2] + 9009*c[0]*c[0]*d[3]
	 - 51480*c[0]*c[1]*d[3] + 5148*c[1]*c[1]*d[3] - 25740*c[0]*c[2]*d[3] + 
	66352*c[1]*c[2]*d[3] + 4004*c[2]*c[2]*d[3] + 117832*c[0]*c[3]*d[3] - 53768*
	c[1]*c[3]*d[3] - 100880*c[2]*c[3]*d[3] + 27612*c[3]*c[3]*d[3] - 28028*c[0]*
	c[4]*d[3] + 68432*c[1]*c[4]*d[3] - 6552*c[2]*c[4]*d[3] - 55792*c[3]*c[4]*
	d[3] + 20548*c[4]*c[4]*d[3] - 4*(429*c[0]*c[0] - 13*c[0]*(429*c[1] - 968*
	c[2] - 539*c[3] + 2284*c[4]) + 2*(3575*c[1]*c[1] - 6994*c[2]*c[2] + 819*
	c[2]*c[3] + 3487*c[3]*c[3] + 12408*c[2]*c[4] - 5137*c[3]*c[4] + 682*c[4]*
	c[4] + 13*c[1]*(55*c[2] - 658*c[3] + 547*c[4])))* d[4] + 1716*(35*d[1] + 
	56*d[2] - 21*d[3] + 4*d[4]) - a*(13*(1155*c[0]*c[0] + 44*(-105 + 63*c[1]*
	c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] + 47*c[3]*c[3] - 6*c[1]*(8*c[2]+3*c[3]))
	 + 8*(308*c[1] + 55*c[2] - 420*c[3])*c[4] + 2188*c[4]*c[4] - 132*c[0]*(28*
	c[1] - 7*c[2] - 12*c[3] + 13*c[4]))*d[1] + 2*(143*(-84 + 21*c[0]*c[0] - 
	48*c[1]*c[1] + 6*c[0]*(7*c[1] - 30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2]) + 
	286*(45*c[0] - 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*(143*c[0] + 
	55*c[1] - 398*c[2] + 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (9009*c[0]*c[0]
	 + 52*(11*(-63 + (9*c[1] + c[2])*(c[1] + 7*c[2])) - 2*(517*c[1] + 328*c[2])
	*c[3] + 531*c[3]*c[3]) + 8*(5460*c[1] - 819*c[2] - 3764*c[3])*c[4] + 20548*
	c[4]*c[4] - 572*c[0]*(36*c[1] + 45*c[2] - 100*c[3] + 49*c[4]))*d[3] + 2*
	(2145*c[0]*c[0] - 26*c[0]*(429*c[1] - 286*c[2] - 539*c[3] + 1126*c[4]) + 
	4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] + 819*c[2]*c[3] + 1882*c[3]*c[3]
	 + 3478*c[2]*c[4] - 5137*c[3]*c[4] + 1649*c[4]*c[4] + 13*c[1]*(55*c[2] - 
	420*c[3] + 547*c[4])))*d[4])))/360360.;
		break;
	case 4:
	case 6:
		bf_mom[0] = ((1 + b)*(429*c[0]*c[0]*(35*(5 + b*(-5 + 2*b))*d[0] + 
	70*(-2 + b)*(1 + b)*d[1] - 14*(7 + b*(-11 + 2*b))*d[2] - 42*(-2 + b)*(1 + b)
	*d[3] - 2*(17 + b*(-1 + 10*b))*d[4]) + 26*c[0]*(33*(7*(-25 - 20*c[1] - 
	14*c[2] + 12*c[3] + b*(25 - 10*c[1] + 22*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] 
	- 3*c[3]) + 6*c[3])) - 2*(17 + b*(-1 + 10*b))*c[4])*d[0] + 2*(33*(70 + 
	126*c[1] - 28*c[2] - 66*c[3] + 52*c[4] + b*(35 - 98*c[1] - 14*c[2] + 78*c[3]
	 + b*(-35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4]))*d[1] + 11*
	(3*(49 - 28*c[1] + 158*c[2] - 60*c[3] + b*(-77 + 14*b - 14*c[1] + 14*b*c[1] 
	- 174*c[2] + 60*b*c[2] + 30*(-1 + b)*c[3])) - 2*(83 + b*(-119 + 26*b))*c[4])
	*d[2] - 11*(2*(63 + 99*c[1] + 90*c[2] - 253*c[3] + 98*c[4]) + b*(63 - 234*
	c[1] + 90*c[2] + 518*c[3] + b*(-63 + 72*c[1] - 90*c[2] - 200*c[3] - 98*c[4])
	 + 98*c[4]))*d[3] - (11*(-51 - 156*c[1] + 166*c[2] + 196*c[3] + b*(3 - 30*b
	 - 78*c[1] + 78*b*c[1] - 238*c[2] + 52*b*c[2] - 98*(-1 + b)*c[3])) - 2*
	(2831 + b*(-2863 + 1126*b))*c[4])*d[4])) + 4*(-18018*c[3]*d[0] - 9009*b*
	c[3]*d[0] + 9009*b*b*c[3]*d[0] + 36179*c[3]*c[3]*d[0] - 37037*b*c[3]*c[3]*
	d[0] + 14300*b*b*c[3]*c[3]*d[0] + 7293*c[4]*d[0] - 429*b*c[4]*d[0] + 4290*b*
	b*c[4]*d[0] - 28028*c[3]*c[4]*d[0] - 14014*b*c[3]*c[4]*d[0] + 14014*b*b*c[3]
	*c[4]*d[0] + 36803*c[4]*c[4]*d[0] - 37219*b*c[4]*c[4]*d[0] + 14638*b*b*c[4]
	*c[4]*d[0] + 28314*c[3]*d[1] - 33462*b*c[3]*d[1] + 10296*b*b*c[3]*d[1] - 
	26884*c[3]*c[3]*d[1] - 13442*b*c[3]*c[3]*d[1] + 13442*b*b*c[3]*c[3]*d[1] - 
	22308*c[4]*d[1] - 11154*b*c[4]*d[1] + 11154*b*b*c[4]*d[1] + 49868*c[3]*c[4]
	*d[1] - 40404*b*c[3]*c[4]*d[1] + 21840*b*b*c[3]*c[4]*d[1] - 28444*c[4]*c[4]
	*d[1] - 14222*b*c[4]*c[4]*d[1] + 14222*b*b*c[4]*c[4]*d[1] + 25740*c[3]*d[2]
	 + 12870*b*c[3]*d[2] - 12870*b*b*c[3]*d[2] - 25402*c[3]*c[3]*d[2] + 33566*
	b*c[3]*c[3]*d[2] - 8528*b*b*c[3]*c[3]*d[2] + 23738*c[4]*d[2] - 34034*b*c[4]
	*d[2] + 7436*b*b*c[4]*d[2] - 6552*c[3]*c[4]*d[2] - 3276*b*c[3]*c[4]*d[2] + 
	3276*b*b*c[3]*c[4]*d[2] - 22842*c[4]*c[4]*d[2] + 33746*b*c[4]*c[4]*d[2] - 
	6956*b*b*c[4]*c[4]*d[2] - 72358*c[3]*d[3] + 74074*b*c[3]*d[3] - 28600*b*b*
	c[3]*d[3] + 27612*c[3]*c[3]*d[3] + 13806*b*c[3]*c[3]*d[3] - 13806*b*b*c[3]*
	c[3]*d[3] + 28028*c[4]*d[3] + 14014*b*c[4]*d[3] - 14014*b*b*c[4]*d[3] - 
	36532*c[3]*c[4]*d[3] + 34316*b*c[3]*c[4]*d[3] - 15056*b*b*c[3]*c[4]*d[3] + 
	20548*c[4]*c[4]*d[3] + 10274*b*c[4]*c[4]*d[3] - 10274*b*b*c[4]*c[4]*d[3] - 
	2*((9133 + b*(-8579 + 3764*b))*c[3]*c[3] + 11*(-2 + b)*(1 + b)*c[3]*(637 + 
	934*c[4]) + c[4]*(13*(2831 + b*(-2863 + 1126*b)) + (5629 + b*(-397+3298*b))
	*c[4]))*d[4] + 143*c[1]*c[1]*(21*(9 + b*(-7 + 4*b))*d[0] + 126*(-2 + b)*
	(1 + b)*d[1] + 6*(15 + b*(-5 + 8*b))*d[2] - 18*(-2 + b)*(1 + b)*d[3] - 2*
	(67 + b*(-61 + 28*b))*d[4]) + 13*c[2]*c[2]*(33*(79 - 87*b + 30*b*b)*d[0] + 
	726*(-2 + b)*(1 + b)*d[1] - 22*(115 + b*(-175 + 34*b))*d[2] - 154*(-2 + b)*
	(1 + b)*d[3] + 2*(1135 + b*(-1415 + 398*b))*d[4]) + c[2]*(143*(21*(7 + b*
	(-11 + 2*b)) + 90*(-2 + b)*(1 + b)*c[3] - 2*(83 + b*(-119 + 26*b))*c[4])*
	d[0] + 26*(11*(42 + 154*c[3] - 20*c[4] + b*(21 - 21*b - 142*c[3] + 64*b*c[3]
	 + 10*(-1 + b)*c[4]))*d[1] + (-33*(79 - 87*b + 30*b*b) - 154*(-2 + b)*
	(1 + b)*c[3] + 2*(1135 + b*(-1415 + 398*b))*c[4])*d[2] + (990 + 495*b - 
	495*b*b - 1954*c[3] + 2582*b*c[3] - 656*b*b*c[3] + 126*(-2 + b)*(1 + b)*
	c[4])*d[3]) + 2*(143*(83 + b*(-119 + 26*b)) + 1638*(-2 + b)*(1 + b)*c[3] - 
	94*(243 + b*(-359 + 74*b))*c[4])*d[4]) + 13*c[1]*(33*(70 - 28*c[2] - 
	66*c[3] + 52*c[4] + b*(35 - 14*c[2] + 78*c[3] + b*(-35 + 14*c[2] - 24*c[3] 
	- 26*c[4]) + 26*c[4]))*d[0] - 66*(63*d[1] - 14*d[2] - 33*d[3] + 26*d[4]) + 
	2*(b*b*(22*(-42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] + 
	64*c[3] + 10*c[4])*d[2] + 2*(198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + 
	(429 + 110*c[2] + 840*c[3] + 1094*c[4])*d[4]) + b*(11*(147 - 30*c[2] + 
	18*c[3] + 122*c[4])*d[1] - 11*(-21 + 66*c[2] + 142*c[3] + 10*c[4])*d[2] - 
	2*(781*c[2] + 517*c[3] + 777*c[4])*d[3] - 2*(55*c[2] + 777*c[3] + 547*c[4])
	*d[4] - 429*(3*d[3] + d[4])) + 2*(11*c[2]*(45*d[1] - 66*d[2] + 77*d[3] - 
	10*d[4]) + c[3]*(198*d[1] + 847*d[2] - 1034*d[3] + 959*d[4]) - c[4]*
	(737*d[1] + 110*d[2] - 959*d[3] + 1094*d[4])))))))/ 2882880.;
		bf_mom[1] = (pow(1 + b,2)*(429*c[0]*c[0]* (7*(5*(-1 + 2*b)*d[0] + 
	6*d[2] + 2*b*(5*d[1] - 2*d[2] - 3*d[3])) - 2*(3 + 10*b)*d[4]) + 26*c[0]*
	 (11*(3*(7*(5 + 6*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3])) - 2*(3 + 10*b)
	*c[4])*d[0] + 2*(-42*c[1]*d[1] + 54*c[3]*d[1] + (-63 - 114*c[2] + 62*c[4])
	*d[2] + 54*c[1]*d[3] - 106*c[3]*d[3] + b*(3*(-35 + 56*c[1] + 14*c[2] - 
	24*c[3] - 26*c[4])*d[1] + 2*(21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])
	*d[2] + (63 - 72*c[1] + 90*c[2] + 200*c[3] + 98*c[4])*d[3]))) - 2*(-99 - 
	682*c[2] + 22*b*(-15 + 39*c[1] + 26*c[2] - 49*c[3]) + 1158*c[4] - 2252*b*
	c[4])*d[4]) + 4*(-9009*c[2]*d[0] - 8151*c[2]*c[2]*d[0] - 7579*c[3]*c[3]*
	d[0] + 1287*c[4]*d[0] + 8866*c[2]*c[4]*d[0] - 7527*c[4]*c[4]*d[0] - 
	7722*c[3]*d[1] - 7436*c[2]*c[3]*d[1] - 6188*c[3]*c[4]*d[1] + 16302*c[2]*
	d[2] + 13442*c[2]*c[2]*d[2] + 8346*c[3]*c[3]*d[2] - 8866*c[4]*d[2] - 
	17628*c[2]*c[4]*d[2] + 8930*c[4]*c[4]*d[2] + 15158*c[3]*d[3] + 16692*c[2]*
	c[3]*d[3] + 6420*c[3]*c[4]*d[3] + 26*c[1]*(11*(21 + 6*c[2] + 22*c[4])*d[1] 
	- (297 + 286*c[2] + 238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4]))
	 - 143*c[1]*c[1]*(21*d[0] - 6*d[2] - 22*d[4]) - 2*(4407*c[2]*c[2] - 
	1605*c[3]*c[3] + c[2]*(4433 - 8930*c[4]) + c[4]*(-7527 + 967*c[4]))*d[4] + 
	b*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] + 4290*c[4]*d[0] + 14014*c[3]*c[4]
	*d[0] + 14638*c[4]*c[4]*d[0] + 10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 
	11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] - 12870*c[3]*
	d[2] - 8528*c[3]*c[3]*d[2] + 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*
	c[4]*c[4]*d[2] - 28600*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 
	15056*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(21 + 45*c[3] - 
	26*c[4])*d[0] + 11*(-21 + 64*c[3] + 10*c[4])*d[1] - 2*(495 + 77*c[3] - 
	398*c[4])*d[2] - (495 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* 
	(42*d[0] + 63*d[1] + 24*d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(1859 + 819*c[3] 
	- 3478*c[4])*d[4] - 2*(3764*c[3]*c[3] + 11*c[3]*(637 + 934*c[4]) + 2*c[4]*
	(7319 + 1649*c[4]))*d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 
	77*d[3] + 398*d[4]) + 13*c[1]* (33*(-35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] 
	+ 2*(22*(-42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] + 64*
	c[3] + 10*c[4])*d[2] + 2*(198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + 
	(429 + 110*c[2] + 840*c[3] + 1094*c[4])*d[4]))))))/ 2882880.;
		bf_mom[2] = (pow(1 + b,2)*(429*c[0]*c[0]* (7*(5*(-1 + 2*b)*d[0] + 
	6*d[2] + 2*b*(5*d[1] - 2*d[2] - 3*d[3])) - 2*(3 + 10*b)*d[4]) + 26*c[0]* 
	(11*(3*(7*(-5 + 6*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3])) - 2*(3 + 10*b)
	*c[4])*d[0] + 2*(-42*c[1]*d[1] + 54*c[3]*d[1] + (63 - 114*c[2] + 62*c[4])
	*d[2] + 54*c[1]*d[3] - 106*c[3]*d[3] + b*(3*(35 + 56*c[1] + 14*c[2] - 24*
	c[3] - 26*c[4])*d[1] + 2*(-21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])*d[2]
	 + (-63 - 72*c[1] + 90*c[2] + 200*c[3] + 98*c[4])*d[3]))) - 2*(99 - 682*c[2]
	 + 22*b*(15 + 39*c[1] + 26*c[2] - 49*c[3]) + 1158*c[4] - 2252*b*c[4])*d[4])
	 + 4*(9009*c[2]*d[0] - 8151*c[2]*c[2]*d[0] - 7579*c[3]*c[3]*d[0] - 1287*c[4]
	*d[0] + 8866*c[2]*c[4]*d[0] - 7527*c[4]*c[4]*d[0] + 7722*c[3]*d[1] - 7436*
	c[2]*c[3]*d[1] - 6188*c[3]*c[4]*d[1] - 16302*c[2]*d[2] + 13442*c[2]*c[2]*
	d[2] + 8346*c[3]*c[3]*d[2] + 8866*c[4]*d[2] - 17628*c[2]*c[4]*d[2] + 8930*
	c[4]*c[4]*d[2] - 15158*c[3]*d[3] + 16692*c[2]*c[3]*d[3] + 6420*c[3]*c[4]*
	d[3] + 26*c[1]*(11*(-21 + 6*c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 238*
	c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) - 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) - 2*(4407*c[2]*c[2] - 1605*c[3]*c[3] + c[4]*
	(7527 + 967*c[4]) - c[2]*(4433 + 8930*c[4]))*d[4] + b*(-9009*c[3]*d[0] + 
	14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 14638*c[4]*
	c[4]*d[0] - 10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 
	21840*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*
	c[3]*d[2] - 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 
	28600*c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*
	d[3] - 10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(-21 + 45*c[3] - 26*c[4])*d[0] + 
	11*(21 + 64*c[3] + 10*c[4])*d[1] - 2*(-495 + 77*c[3] - 398*c[4])*d[2] - 
	(-495 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] + 63*d[1] + 
	24*d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(-1859 + 819*c[3] - 3478*c[4])*d[4] - 
	2*(3764*c[3]*c[3] + 11*c[3]*(-637 + 934*c[4]) + 2*c[4]*(-7319 + 1649*c[4]))
	*d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 77*d[3] + 398*d[4]) 
	+ 13*c[1]* (33*(35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] + 2*(22*(42 + 24*c[2]
	 - 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*
	(-198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] + 840*c[3]
	 + 1094*c[4])*d[4]))))))/ 2882880.;
		bf_mom[3] = ((1 + b)*(429*c[0]*c[0]*(35*(5 + b*(-5 + 2*b))*d[0] + 
	70*(-2 + b)*(1 + b)*d[1] - 14*(7 + b*(-11 + 2*b))*d[2] - 42*(-2 + b)*(1 + b)
	*d[3] - 2*(17 + b*(-1 + 10*b))*d[4]) + 26*c[0]*(33*(7*(25 - 20*c[1] - 14*
	c[2] + 12*c[3] + b*(-25 - 10*c[1] + 22*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 
	3*c[3]) + 6*c[3])) - 2*(17 + b*(-1 + 10*b))*c[4])*d[0] + 2*(33*(-70 + 126*
	c[1] - 28*c[2] - 66*c[3] + 52*c[4] + b*(-35 - 98*c[1] - 14*c[2] + 78*c[3] 
	+ b*(35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4]))*d[1] + 11*(3*
	(-49 - 28*c[1] + 158*c[2] - 60*c[3] + b*(77 - 14*b - 14*c[1] + 14*b*c[1] - 
	174*c[2] + 60*b*c[2] + 30*(-1 + b)*c[3])) - 2*(83 + b*(-119 + 26*b))*c[4])
	*d[2] - 11*(2*(-63 + 99*c[1] + 90*c[2] - 253*c[3] + 98*c[4]) + b*(-63 - 
	234*c[1] + 90*c[2] + 518*c[3] + b*(63 + 72*c[1] - 90*c[2] - 200*c[3] - 
	98*c[4]) + 98*c[4]))*d[3] - (11*(51 - 156*c[1] + 166*c[2] + 196*c[3] + b*
	(-3 + 30*b - 78*c[1] + 78*b*c[1] - 238*c[2] + 52*b*c[2] - 98*(-1 + b)*c[3]))
	 - 2*(2831 + b*(-2863 + 1126*b))*c[4])*d[4])) + 4*(18018*c[3]*d[0] + 9009*
	b*c[3]*d[0] - 9009*b*b*c[3]*d[0] + 36179*c[3]*c[3]*d[0] - 37037*b*c[3]*c[3]
	*d[0] + 14300*b*b*c[3]*c[3]*d[0] - 7293*c[4]*d[0] + 429*b*c[4]*d[0] - 4290*
	b*b*c[4]*d[0] - 28028*c[3]*c[4]*d[0] - 14014*b*c[3]*c[4]*d[0] + 14014*b*b*
	c[3]*c[4]*d[0] + 36803*c[4]*c[4]*d[0] - 37219*b*c[4]*c[4]*d[0] + 14638*b*
	b*c[4]*c[4]*d[0] - 28314*c[3]*d[1] + 33462*b*c[3]*d[1] - 10296*b*b*c[3]*
	d[1] - 26884*c[3]*c[3]*d[1] - 13442*b*c[3]*c[3]*d[1] + 13442*b*b*c[3]*c[3]*
	d[1] + 22308*c[4]*d[1] + 11154*b*c[4]*d[1] - 11154*b*b*c[4]*d[1] + 49868*
	c[3]*c[4]*d[1] - 40404*b*c[3]*c[4]*d[1] + 21840*b*b*c[3]*c[4]*d[1] - 
	28444*c[4]*c[4]*d[1] - 14222*b*c[4]*c[4]*d[1] + 14222*b*b*c[4]*c[4]*d[1] - 
	25740*c[3]*d[2] - 12870*b*c[3]*d[2] + 12870*b*b*c[3]*d[2] - 25402*c[3]*
	c[3]*d[2] + 33566*b*c[3]*c[3]*d[2] - 8528*b*b*c[3]*c[3]*d[2] - 23738*c[4]*
	d[2] + 34034*b*c[4]*d[2] - 7436*b*b*c[4]*d[2] - 6552*c[3]*c[4]*d[2] - 3276*
	b*c[3]*c[4]*d[2] + 3276*b*b*c[3]*c[4]*d[2] - 22842*c[4]*c[4]*d[2] + 33746*
	b*c[4]*c[4]*d[2] - 6956*b*b*c[4]*c[4]*d[2] + 72358*c[3]*d[3] - 74074*b*c[3]
	*d[3] + 28600*b*b*c[3]*d[3] + 27612*c[3]*c[3]*d[3] + 13806*b*c[3]*c[3]*d[3]
	 - 13806*b*b*c[3]*c[3]*d[3] - 28028*c[4]*d[3] - 14014*b*c[4]*d[3] + 14014*
	b*b*c[4]*d[3] - 36532*c[3]*c[4]*d[3] + 34316*b*c[3]*c[4]*d[3] - 15056*b*b*
	c[3]*c[4]*d[3] + 20548*c[4]*c[4]*d[3] + 10274*b*c[4]*c[4]*d[3] - 10274*b*
	b*c[4]*c[4]*d[3] - 2*((9133 + b*(-8579 + 3764*b))*c[3]*c[3] + 11*(-2 + b)*
	(1 + b)*c[3]*(-637 + 934*c[4]) + c[4]*(-13*(2831 + b*(-2863 + 1126*b)) + 
	(5629 + b*(-397 + 3298*b))*c[4]))*d[4] + 143*c[1]*c[1]*(21*(9 + b*(-7 + 
	4*b))*d[0] + 126*(-2 + b)*(1 + b)*d[1] + 6*(15 + b*(-5 + 8*b))*d[2] - 18*
	(-2 + b)*(1 + b)*d[3] - 2*(67 + b*(-61 + 28*b))*d[4]) + 13*c[2]*c[2]*(33*
	(79 - 87*b + 30*b*b)*d[0] + 726*(-2 + b)*(1 + b)*d[1] - 22*(115 + b*(-175 
	+ 34*b))*d[2] - 154*(-2 + b)*(1 + b)*d[3] + 2*(1135 + b*(-1415 + 398*b))
	*d[4]) + c[2]*(143*(-21*(7 + b*(-11 + 2*b)) + 90*(-2 + b)*(1 + b)*c[3] - 
	2*(83 + b*(-119 + 26*b))*c[4])*d[0] + 26*(11*(-42 + 154*c[3] - 20*c[4] + 
	b*(-21 + 21*b - 142*c[3] + 64*b*c[3] + 10*(-1 + b)*c[4]))*d[1] + (33*(79 - 
	87*b + 30*b*b) - 154*(-2 + b)*(1 + b)*c[3] + 2*(1135 + b*(-1415 + 398*b))*
	c[4])*d[2] - (990 + 1954*c[3] + 252*c[4] + b*(495 - 495*b - 2582*c[3] + 
	656*b*c[3] - 126*(-1 + b)*c[4]))*d[3]) + 2*(-143*(83 + b*(-119 + 26*b)) + 
	1638*(-2 + b)*(1 + b)*c[3] - 94*(243 + b*(-359 + 74*b))*c[4])*d[4]) + 13*
	c[1]*(33*(-70 - 28*c[2] - 66*c[3] + 52*c[4] + b*(-35 - 14*c[2] + 78*c[3] + 
	b*(35 + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4]))*d[0] + 66*(63*d[1] - 14*
	d[2] - 33*d[3] + 26*d[4]) + 2*(b*b*(22*(42 + 24*c[2] - 9*c[3] - 28*c[4])*
	d[1] + 11*(21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(-198 + 352*c[2] + 
	517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] + 840*c[3] + 1094*c[4])*d[4])
	 + b*(-11*(147 + 30*c[2] - 18*c[3] - 122*c[4])*d[1] - 11*(21 + 66*c[2] + 
	142*c[3] + 10*c[4])*d[2] - 2*(781*c[2] + 517*c[3] + 777*c[4])*d[3] - 2*
	(55*c[2] + 777*c[3] + 547*c[4])*d[4] + 429*(3*d[3] + d[4])) + 2*(11*c[2]*
	(45*d[1] - 66*d[2] + 77*d[3] - 10*d[4]) + c[3]*(198*d[1] + 847*d[2] - 
	1034*d[3] + 959*d[4]) - c[4]*(737*d[1] + 110*d[2] - 959*d[3] + 
	1094*d[4])))))))/ 2882880.;
		bf_mom[4] = -(pow(1 + b,2)*(429*c[0]*c[0]* (7*(5*(-2 + b)*d[0] + 
	5*(-1 + b)*d[1] - 2*(-4 + b)*d[2] - 3*(-1 + b)*d[3]) - 2*(-2 + 5*b)*d[4]) 
	+ 26*c[0]*(33*(7*(10 - 5*c[1] + 8*c[2] + b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) 
	+ 3*c[3]) - 2*(-2 + 5*b)*c[4])*d[0] + 11*(3*(35 - 84*c[1] - 14*c[2] + 60*
	c[3] + b*(-35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] + 
	2*(-3*(28 + 7*c[1] + 68*c[2] + 15*c[3]) + b*(21 + 21*c[1] + 90*c[2] + 
	45*c[3] - 26*c[4]) + 88*c[4])*d[2] + (-63 + 63*b + 180*c[1] - 72*b*c[1] - 
	90*c[2] + 90*b*c[2] - 412*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])*d[3]) - 
	2*(11*(6 - 39*c[1] - 88*c[2] + b*(-15 + 39*c[1] + 26*c[2] - 49*c[3]) + 
	49*c[3]) - 2*(-1142 + 563*b)*c[4])*d[4]) + 2*(-9009*c[3]*d[0] + 9009*b*c[3]
	*d[0] - 29458*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] - 1716*c[4]*d[0] + 
	4290*b*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 29692*
	c[4]*c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] - 25740*c[3]*d[1] + 10296*b*c[3]*
	d[1] - 13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 
	11154*b*c[4]*d[1] - 34216*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 
	14222*c[4]*c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 12870*b*
	c[3]*d[2] + 25220*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] - 25168*c[4]*d[2] 
	+ 7436*b*c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 24816*
	c[4]*c[4]*d[2] - 6956*b*c[4]*c[4]*d[2] + 58916*c[3]*d[3] - 28600*b*c[3]*
	d[3] + 13806*c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 
	14014*b*c[4]*d[3] + 27896*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 
	10274*c[4]*c[4]*d[3] - 10274*b*c[4]*c[4]*d[3] + 26*c[2]*(11*(-84 - 45*c[3] 
	+ b*(21 + 45*c[3] - 26*c[4]) + 88*c[4])*d[0] + 11*(21 - 21*b - 116*c[3] + 
	64*b*c[3] + 10*(-1 + b)*c[4])* d[1] - 2*(-1122 - 77*c[3] + b*(495 + 77*c[3]
	 - 398*c[4]) + 1076*c[4])*d[2] - (-495 + 495*b - 1940*c[3] + 656*b*c[3] - 
	126*(-1 + b)*c[4])* d[3]) + 4*c[2]*(143*(-44 + 13*b) + 819*(-1 + b)*c[3] - 
	94*(-132 + 37*b)*c[4])*d[4] - 2*((-6974 + 3764*b)*c[3]*c[3] + 11*(-1 + b)*
	c[3]*(637 + 934*c[4]) + 2*c[4]*(-14846 + 7319*b - 682*c[4] + 1649*b*c[4]))
	*d[4] + 286*c[1]*c[1]*(21*(-3 + 2*b)*d[0] + 63*(-1 + b)*d[1] - 18*d[2] + 
	24*b*d[2] + 9*d[3] - 9*b*d[3] + 50*d[4] - 28*b*d[4]) + 26*c[2]*c[2]*(11*
	(3*(-34 + 15*b)*d[0] + 33*(-1 + b)*d[1] - 2*(-64 + 17*b)*d[2] - 7*(-1 + b)
	*d[3]) + 2*(-538 + 199*b)*d[4]) + 13*c[1]*(33*(35 - 14*c[2] + 60*c[3] + 
	b*(-35 + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] + 66*(42*d[1] + 
	7*d[2] - 30*d[3] - 13*d[4]) + 2*(b*(22*(-42 + 24*c[2] - 9*c[3] - 28*c[4])
	*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(198 + 352*c[2] + 
	517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] + 840*c[3] + 1094*c[4])*d[4]) 
	- 2*(11*c[2]*(18*d[1] + 33*d[2] + 58*d[3] + 5*d[4]) + c[4]*(-550*d[1] + 
	55*d[2] + 658*d[3] + 547*d[4]) + c[3]*(-99*d[1] + 638*d[2] + 517*d[3] 
	+ 658*d[4])))))))/ 720720.;
		bf_mom[5] = -(pow(1 + b,2)*(13*(33*(7*(5*(-1 + 2*b)*(-4 + c[0]*c[0])
	 + 20*b*c[0]*c[1] + 4*(-1 + 4*b)*c[1]*c[1]) - 28*((-3 + 2*b)*c[0]-2*b*c[1])
	*c[2] + 4*(-19 + 30*b)*c[2]*c[2]) - 396*(-6*c[1] + b*(7*c[0] + 8*c[1] - 
	10*c[2]))*c[3] + 44*(-53 + 100*b)*c[3]*c[3] - 44*((9 + 30*b)*c[0] - 62*c[2]
	 + b*(78*c[1] + 52*c[2] - 98*c[3]))*c[4] + 4*(-579 + 1126*b)*c[4]*c[4])*d[0]
	 + 2*(-5148*(7*d[2] - d[4]) + 1287*c[0]*c[0]*(7*d[2] - d[4]) - 52*c[0]*
	(627*c[2]*d[2] - 341*c[4]*d[2] + 33*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-297*d[1]
	 + 583*d[3]) - 341*c[2]*d[4] + 579*c[4]*d[4]) + b*(13*(33*(-140 + 35*c[0]*
	c[0] + 84*c[1]*c[1] + 64*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2]))
	 - 88*(18*c[0] + 9*c[1] - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 
	616*c[1] - 110*c[2] - 840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 
	+ 21*c[0]*c[0] - 48*c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*
	(7*c[1] + 30*c[2])) - 286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*
	c[3] + 52*(143*c[0] - 55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])
	*d[2] - (13*(99*(-28 + (c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 
	64*c[1])*c[2] + 308*c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 
	2124*c[3]*c[3]) - 4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])*c[4] + 
	20548*c[4]*c[4])*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 
	539*c[3] - 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*
	c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*
	c[4] - 13*c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4]) + 4*(-3094*c[3]*c[4]
	*d[1] + 4173*c[3]*c[3]*d[2] + 4465*c[4]*c[4]*d[2] + 3210*c[3]*c[4]*d[3] + 
	13*c[2]*c[2]*(517*d[2] - 339*d[4]) + 1605*c[3]*c[3]*d[4] - 967*c[4]*c[4]*
	d[4] + 143*c[1]*c[1]*(3*d[2] + 11*d[4]) + 26*c[1]*(121*c[4]*d[1] - 143*c[3]
	*d[2] + 11*c[2]*(3*d[1] - 13*d[3]) - 119*c[4]*d[3] - 119*c[3]*d[4]) + c[2]
	*(-3718*c[3]*d[1] - 8814*c[4]*d[2] + 8346*c[3]*d[3] + 
	8930*c[4]*d[4])))))/1441440.;
		bf_mom[6] = -(pow(1 + b,2)*(429*c[0]*c[0]* (7*(5*(-2 + b)*d[0] + 
	5*(-1 + b)*d[1] - 2*(-4 + b)*d[2] - 3*(-1 + b)*d[3]) - 2*(-2 + 5*b)*d[4]) 
	+ 26*c[0]*(33*(7*(-10 - 5*c[1] + 8*c[2] + b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) 
	+ 3*c[3]) - 2*(-2 + 5*b)*c[4])*d[0] + 11*(3*(-35 - 84*c[1] - 14*c[2] + 
	60*c[3] + b*(35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] 
	+ 2*(84 - 21*c[1] - 204*c[2] - 45*c[3] + b*(-21 + 21*c[1] + 90*c[2] + 
	45*c[3] - 26*c[4]) + 88*c[4])*d[2] + (63 - 63*b + 180*c[1] - 72*b*c[1] - 
	90*c[2] + 90*b*c[2] - 412*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])*d[3]) - 
	2*(11*(-6 - 39*c[1] - 88*c[2] + b*(15 + 39*c[1] + 26*c[2] - 49*c[3]) + 
	49*c[3]) - 2*(-1142 + 563*b)*c[4])*d[4]) + 2*(9009*c[3]*d[0] - 9009*b*c[3]*
	d[0] - 29458*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] + 1716*c[4]*d[0] - 
	4290*b*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 29692*
	c[4]*c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] + 25740*c[3]*d[1] - 10296*b*c[3]*
	d[1] - 13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] + 11154*c[4]*d[1] - 
	11154*b*c[4]*d[1] - 34216*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 14222*
	c[4]*c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] - 12870*c[3]*d[2] + 12870*b*c[3]*
	d[2] + 25220*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] + 25168*c[4]*d[2] - 
	7436*b*c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 24816*
	c[4]*c[4]*d[2] - 6956*b*c[4]*c[4]*d[2] - 58916*c[3]*d[3] + 28600*b*c[3]*
	d[3] + 13806*c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] - 14014*c[4]*d[3] + 
	14014*b*c[4]*d[3] + 27896*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 
	10274*c[4]*c[4]*d[3] - 10274*b*c[4]*c[4]*d[3] + 26*c[2]*(11*(84 - 45*c[3] 
	+ b*(-21 + 45*c[3] - 26*c[4]) + 88*c[4])*d[0] + 11*(-21 + 21*b - 116*c[3] 
	+ 64*b*c[3] + 10*(-1 + b)*c[4])* d[1] - 2*(1122 - 77*c[3] + b*(-495 + 
	77*c[3] - 398*c[4]) + 1076*c[4])*d[2] - (495 - 495*b - 1940*c[3] + 656*b*
	c[3] - 126*(-1 + b)*c[4])* d[3]) + 4*c[2]*(6292 - 1859*b + 819*(-1 + b)*
	c[3] - 94*(-132 + 37*b)*c[4])*d[4] - 2*((-6974 + 3764*b)*c[3]*c[3] + 11*
	(-1 + b)*c[3]*(-637 + 934*c[4]) + 2*c[4]*(14846 - 682*c[4] + b*(-7319 + 
	1649*c[4])))*d[4] + 286*c[1]*c[1]*(21*(-3 + 2*b)*d[0] + 63*(-1 + b)*d[1] 
	- 18*d[2] + 24*b*d[2] + 9*d[3] - 9*b*d[3] + 50*d[4] - 28*b*d[4]) + 26*c[2]*
	c[2]*(11*(3*(-34 + 15*b)*d[0] + 33*(-1 + b)*d[1] - 2*(-64 + 17*b)*d[2] - 
	7*(-1 + b)*d[3]) + 2*(-538 + 199*b)*d[4]) + 13*c[1]*(33*(-35 - 14*c[2] + 
	60*c[3] + b*(35 + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] - 66*
	(42*d[1] + 7*d[2] - 30*d[3] - 13*d[4]) + 2*(b*(22*(42 + 24*c[2] - 9*c[3] 
	- 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(-198 + 
	352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] + 840*c[3] + 1094*
	c[4])*d[4]) - 2*(11*c[2]*(18*d[1] + 33*d[2] + 58*d[3] + 5*d[4]) + c[4]*
	(-550*d[1] + 55*d[2] + 658*d[3] + 547*d[4]) + c[3]*(-99*d[1] + 638*d[2] 
	+ 517*d[3] + 658*d[4])))))))/ 720720.;
		bf_mom[7] = -((1 + b)*(13*(5775*(-4 + c[0]*c[0]) + 44*(189*c[1]*c[1]
	 + 237*c[2]*c[2] - 21*c[0]*(10*c[1] + 7*c[2] - 6*c[3]) - 180*c[2]*c[3] + 
	253*c[3]*c[3] - 6*c[1]*(14*c[2] + 33*c[3])) - 44*(51*c[0] + 2*(-78*c[1] + 
	83*c[2] + 98*c[3]))*c[4] + 11324*c[4]*c[4] + b*(-5775*c[0]*c[0] - 924*c[0]*
	(5*c[1] - 11*c[2] - 3*c[3]) - 44*(-525 + 147*c[1]*c[1] + 261*c[2]*c[2] + 
	6*c[1]*(7*c[2] - 39*c[3]) + 90*c[2]*c[3] + 259*c[3]*c[3]) + 44*(3*c[0] + 
	78*c[1] + 238*c[2] - 98*c[3])*c[4] - 11452*c[4]*c[4]) + 2*b*b*(1155*c[0]*
	c[0] + 462*c[0]*(5*c[1] - 2*c[2] - 3*c[3]) + 44*(42*c[1]*c[1] + 3*c[1]*
	(7*c[2] - 12*c[3]) + 5*(-21 + 9*c[2]*c[2] + 9*c[2]*c[3] + 10*c[3]*c[3])) - 
	44*(15*c[0] + 39*c[1] + 26*c[2] - 49*c[3])*c[4] + 2252*c[4]*c[4]))*d[0] + 
	3432*(70*d[1] + 49*d[2] - 42*d[3] + 17*d[4]) + 2*(-429*c[0]*c[0]*(70*d[1] 
	+ 49*d[2] - 42*d[3] + 17*d[4]) + b*b*(13*(33* (-140 + 35*c[0]*c[0] + 84*
	c[1]*c[1] + 64*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2])) - 88*
	(18*c[0] + 9*c[1] - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 616*c[1]
	 - 110*c[2] - 840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 + 21*
	c[0]*c[0] - 48*c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*(7*c[1] 
	+ 30*c[2])) - 286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*
	(143*c[0] - 55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (13*
	(99*(-28 + (c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2]
	 + 308*c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3])
	 - 4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])*c[4] + 20548*c[4]*
	c[4])*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 539*c[3] 
	- 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*c[2]*c[3] 
	+ 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*c[4] - 
	13*c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4]) - b*(13*(1155*c[0]*c[0] + 
	44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(5*c[2] - 3*c[3]) + 142*
	c[2]*c[3] + 47*c[3]*c[3]) - 8*(671*c[1] - 55*c[2] - 777*c[3])*c[4] + 2188*
	c[4]*c[4] + 132*c[0]*(49*c[1] + 7*c[2] - 13*(3*c[3] + c[4])))*d[1] + 
	(-33033*c[0]*c[0] + 52*(11*(231 + 15*c[1]*c[1] + 66*c[1]*c[2] - 175*c[2]*
	c[2]) + 22*(71*c[1] - 7*c[2])*c[3] - 1291*c[3]*c[3]) + 572*c[0]*(21*c[1] 
	+ 261*c[2] + 45*c[3] - 119*c[4]) + 104*(55*c[1] + 1415*c[2] + 63*c[3])*c[4]
	 - 67492*c[4]*c[4])*d[2] - (9009*c[0]*c[0] + 52*(99*c[1]*c[1] + 77*(-9 + 
	c[2]*c[2]) + 2582*c[2]*c[3] + 531*c[3]*c[3] - 22*c[1]*(71*c[2] + 47*c[3]))
	 + 572*c[0]*(117*c[1] - 45*c[2] - 259*c[3] - 49*c[4]) - 8*(10101*c[1] + 
	819*c[2] - 8579*c[3])*c[4] + 20548*c[4]*c[4])*d[3] - (429*c[0]*c[0] + 52*
	c[0]*(429*c[1] + 7*(187*c[2] - 77*c[3] - 409*c[4])) + 4*(-429 + 8723*c[1]*
	c[1] - 18395*c[2]*c[2] - 1638*c[2]*c[3] + 8579*c[3]*c[3] + 33746*c[2]*c[4]
	 + 10274*c[3]*c[4] + 397*c[4]*c[4] - 26*c[1]*(55*c[2] + 777*c[3] + 547*
	c[4])))*d[4]) + 52*c[0]*(-11*(c[3]*(99*d[1] + 90*d[2] - 253*d[3]) + c[4]*
	(-78*d[1] + 83*d[2] + 98*d[3])) + (-1078*c[3] + 2831*c[4])*d[4] + 33*c[1]*
	(63*d[1] - 14*d[2] - 33*d[3] + 26*d[4]) - 11*c[2]*(42*d[1] - 237*d[2] + 
	90*d[3] + 83*d[4])) + 4*(-13442*c[3]*c[3]*d[1] + 24934*c[3]*c[4]*d[1] - 
	14222*c[4]*c[4]*d[1] - 12701*c[3]*c[3]*d[2] - 3276*c[3]*c[4]*d[2] - 11421*
	c[4]*c[4]*d[2] + 13806*c[3]*c[3]*d[3] - 18266*c[3]*c[4]*d[3] + 10274*c[4]*
	c[4]*d[3] - 13*c[2]*c[2]* (726*d[1] + 1265*d[2] - 154*d[3] - 1135*d[4]) - 
	(9133*c[3]*c[3] - 20548*c[3]*c[4] + 5629*c[4]*c[4])* d[4] - 143*c[1]*c[1]*
	 (126*d[1] - 45*d[2] - 18*d[3] + 67*d[4]) + 26*c[1]*(11*c[2]*(45*d[1] - 
	66*d[2] + 77*d[3] - 10*d[4]) + c[3]*(198*d[1] + 847*d[2] - 1034*d[3] + 
	959*d[4]) - c[4]*(737*d[1] + 110*d[2] - 959*d[3] + 1094*d[4])) + 2*c[2]*
	(13*c[3]*(847*d[1] + 154*d[2] - 977*d[3] - 126*d[4]) - c[4]*(1430*d[1] - 
	14755*d[2] + 1638*d[3] + 11421*d[4]))))))/ 1441440.;
		bf_mom[8] = (pow(1 + b,2)*(13*(33*(35*(-2 + b)*(-4 + c[0]*c[0]) + 
	70*(-1 + b)*c[0]*c[1] + 28*(-3 + 2*b)*c[1]*c[1] - 28*((-4 + b)*c[0] + c[1] 
	- b*c[1])*c[2] + 4*(-34 + 15*b)*c[2]*c[2]) - 198*(7*(-1 + b)*c[0] + 4*(-5 
	+ 2*b)*c[1] - 10*(-1 + b)*c[2])* c[3] + 44*(-103 + 50*b)*c[3]*c[3] - 44*(3*
	(-2 + 5*b)*c[0] + 39*(-1 + b)*c[1] + 2*(-44 + 13*b)*c[2] - 49*(-1 + b)*c[3])
	*c[4] + 4*(-1142 + 563*b)*c[4]*c[4])*d[0] - 15015*c[0]*c[0]*d[1] - 72072*
	c[0]*c[1]*d[1] - 36036*c[1]*c[1]*d[1] - 12012*c[0]*c[2]*d[1] - 20592*c[1]*
	c[2]*d[1] - 18876*c[2]*c[2]*d[1] + 51480*c[0]*c[3]*d[1] + 10296*c[1]*c[3]*
	d[1] - 66352*c[2]*c[3]*d[1] - 26884*c[3]*c[3]*d[1] + 22308*c[0]*c[4]*d[1] 
	+ 57200*c[1]*c[4]*d[1] - 5720*c[2]*c[4]*d[1] - 68432*c[3]*c[4]*d[1] - 
	28444*c[4]*c[4]*d[1] + 24024*c[0]*c[0]*d[2] - 12012*c[0]*c[1]*d[2] - 10296*
	c[1]*c[1]*d[2] - 116688*c[0]*c[2]*d[2] - 37752*c[1]*c[2]*d[2] + 73216*c[2]*
	c[2]*d[2] - 25740*c[0]*c[3]*d[2] - 66352*c[1]*c[3]*d[2] + 8008*c[2]*c[3]*
	d[2] + 50440*c[3]*c[3]*d[2] + 50336*c[0]*c[4]*d[2] - 5720*c[1]*c[4]*d[2] 
	- 111904*c[2]*c[4]*d[2] - 6552*c[3]*c[4]*d[2] + 49632*c[4]*c[4]*d[2] + 
	9009*c[0]*c[0]*d[3] + 51480*c[0]*c[1]*d[3] + 5148*c[1]*c[1]*d[3] - 25740*
	c[0]*c[2]*d[3] - 66352*c[1]*c[2]*d[3] + 4004*c[2]*c[2]*d[3] - 117832*c[0]*
	c[3]*d[3] - 53768*c[1]*c[3]*d[3] + 100880*c[2]*c[3]*d[3] + 27612*c[3]*c[3]
	*d[3] - 28028*c[0]*c[4]*d[3] - 68432*c[1]*c[4]*d[3] - 6552*c[2]*c[4]*d[3] 
	+ 55792*c[3]*c[4]*d[3] + 20548*c[4]*c[4]*d[3] + 1716*(35*d[1] - 56*d[2] - 
	21*d[3] - 4*d[4]) + 4*(429*c[0]*c[0] + 13*c[0]* (429*c[1] + 968*c[2] - 
	539*c[3] - 2284*c[4]) + 2*(3575*c[1]*c[1] - 6994*c[2]*c[2] - 819*c[2]*c[3] 
	+ 3487*c[3]*c[3] + 12408*c[2]*c[4] + 5137*c[3]*c[4] + 682*c[4]*c[4] - 
	13*c[1]*(55*c[2] + 658*c[3] + 547*c[4])))* d[4] + b*(13*(33*(-140 + 35*c[0]
	*c[0] + 84*c[1]*c[1] + 64*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1]+c[2]))
	- 88*(18*c[0] + 9*c[1] - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 
	616*c[1] - 110*c[2] - 840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 
	+ 21*c[0]*c[0] - 48*c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*
	(7*c[1] + 30*c[2])) - 286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]
	*c[3] + 52*(143*c[0] - 55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])
	*d[2] - (13*(99*(-28 + (c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 
	64*c[1])*c[2] + 308*c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 
	2124*c[3]*c[3]) - 4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])*c[4] + 
	20548*c[4]*c[4])*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] 
	- 539*c[3] - 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 
	819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*
	c[4]*c[4] - 13*c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4])))/360360.;
		break;
	case 22:
		bf_mom[0] = (429*c[0]*c[0]*(35*d[0] - 70*d[1] + 14*d[2] + 42*d[3] 
	- 26*d[4]) + 26*c[0]*(-33*(35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4])*d[0]
	 + 66*(35 + 42*c[1] - 14*c[2] - 6*c[3] + 26*c[4])*d[1] - 22*(21 + 42*c[1] 
	- 66*c[2] + 90*c[3] - 10*c[4])*d[2] - 22*(63 + 18*c[1] + 90*c[2] - 94*c[3] 
	+ 98*c[4])*d[3] + 2*(11*(39 + 78*c[1] + 10*c[2] - 98*c[3])+1094*c[4])*d[4])
	 + 4*(-9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] + 5577*c[4]*d[0] - 14014*c[3]*
	c[4]*d[0] + 7111*c[4]*c[4]*d[0] + 2574*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*
	d[2] - 182*c[3]*c[3]*d[2] - 1430*c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 1974*
	c[4]*c[4]*d[2] - 13442*c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] 
	- 8636*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 143*c[1]*c[1]*(9*(7*d[0] + 
	2*(-7*d[1] + 3*d[2] + d[3])) - 34*d[4]) - 2*(2159*c[3]*c[3] - 11*c[3]*
	(637 + 934*c[4]) + c[4]*(7111 + 4265*c[4]))*d[4] + 13*c[2]*c[2]*(363*d[0] 
	- 726*d[1] + 286*d[2] + 154*d[3] + 118*d[4]) - 13*c[1]* (33*(-35 + 14*c[2] 
	+ 6*c[3] - 26*c[4])*d[0] + 2*(-11*(-63 + 54*c[2] + 18*c[3] - 34*c[4])*d[1] 
	+ 11*(-21 + 66*c[2] - 38*c[3] + 10*c[4])*d[2] - (99 + 418*c[2] - 1034*c[3] 
	+ 602*c[4])*d[3] + (429 + 110*c[2] - 602*c[3] + 1094*c[4])*d[4])) + c[2]*
	(-143*(21 + 90*c[3] - 10*c[4])*d[0] - 4*c[4]*(715*d[1] - 767*d[2] + 819*
	d[3] - 987*d[4]) + 286*(21*d[1] - 33*d[2] + 45*d[3] - 5*d[4]) + 52*c[3]*
	(209*d[1] + 77*d[2] - 7*(d[3] + 9*d[4])))))/720720.;
		bf_mom[1] = (429*c[0]*c[0]*(35*d[0] - 70*d[1] + 14*d[2] + 42*d[3] 
	- 26*d[4]) + 26*c[0]*(11*(-3*(70*c[1] - 7*(5 + 2*c[2] + 6*c[3]) + 26*c[4])
	*d[0] + 2*(3*(-35 + 42*c[1] - 14*c[2] - 6*c[3] + 26*c[4])*d[1] + (21 - 42*
	c[1] + 66*c[2] - 90*c[3] + 10*c[4])*d[2] + (63 - 18*c[1] - 90*c[2] + 94*c[3]
	 - 98*c[4])*d[3])) + 2*(11*(-39 + 78*c[1] + 10*c[2] - 98*c[3]) + 1094*c[4])
	*d[4]) + 4*(9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] - 5577*c[4]*d[0] - 14014*
	c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] - 2574*c[3]*d[1] - 13442*c[3]*c[3]*d[1]
	 + 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 12870*
	c[3]*d[2] - 182*c[3]*c[3]*d[2] + 1430*c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 
	1974*c[4]*c[4]*d[2] + 13442*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 14014*c[4]*
	d[3] - 8636*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 143*c[1]*c[1]*(9*
	(7*d[0] + 2*(-7*d[1] + 3*d[2] + d[3])) - 34*d[4]) - 2*(2159*c[3]*c[3] - 
	11*c[3]*(-637 + 934*c[4]) + c[4]*(-7111 + 4265*c[4]))*d[4] + 13*c[2]*c[2]*
	(363*d[0] - 726*d[1] + 286*d[2] + 154*d[3] + 118*d[4]) - 13*c[1]* (33*(35 
	+ 14*c[2] + 6*c[3] - 26*c[4])*d[0] + 2*(-11*(63 + 54*c[2] + 18*c[3] - 34*
	c[4])*d[1] + 11*(21 + 66*c[2] - 38*c[3] + 10*c[4])*d[2] - (-99 + 418*c[2] 
	- 1034*c[3] + 602*c[4])*d[3] + (-429 + 110*c[2] - 602*c[3] + 1094*c[4])
	*d[4])) + c[2]*(-143*(-21 + 90*c[3] - 10*c[4])*d[0] - 4*c[4]*(715*d[1] 
	- 767*d[2] + 819*d[3] - 987*d[4]) - 286*(21*d[1] - 33*d[2] + 45*d[3] - 
	5*d[4]) + 52*c[3]*(209*d[1] + 77*d[2] - 7*(d[3] + 9*d[4])))))/720720.;
		bf_mom[2] = (429*c[0]*c[0]*(35*d[0] + 14*(5*d[1] + d[2] - 3*d[3]) 
	- 26*d[4]) + 26*c[0]*(11*(3*(35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4])
	*d[0] + 2*(3*(35 + 42*c[1] + 14*c[2] - 6*c[3] - 26*c[4])*d[1] + (21 + 42*
	c[1] + 66*c[2] + 90*c[3] + 10*c[4])*d[2] + (-63 - 18*c[1] + 90*c[2] + 
	94*c[3] + 98*c[4])*d[3])) - 2*(11*(39 + 78*c[1] - 10*c[2] - 98*c[3]) - 
	1094*c[4])*d[4]) + 4*(-9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] - 5577*c[4]*
	d[0] + 14014*c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] - 2574*c[3]*d[1] + 
	13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] + 14222*c[4]
	*c[4]*d[1] + 12870*c[3]*d[2] - 182*c[3]*c[3]*d[2] + 1430*c[4]*d[2] + 3276
	*c[3]*c[4]*d[2] + 1974*c[4]*c[4]*d[2] + 13442*c[3]*d[3] - 13806*c[3]*c[3]
	*d[3] + 14014*c[4]*d[3] - 8636*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 
	143*c[1]*c[1]*(9*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) - 34*d[4]) - 2*
	(2159*c[3]*c[3] + 11*c[3]*(-637 + 934*c[4]) + c[4]*(-7111 + 4265*c[4]))
	*d[4] + 13*c[2]*c[2]*(363*d[0] + 726*d[1] + 286*d[2] - 154*d[3] + 118*d[4])
	 + 13*c[1]* (33*(35 + 14*c[2] - 6*c[3] - 26*c[4])*d[0] + 2*(11*(63 + 54*c[2]
	 - 18*c[3] - 34*c[4])*d[1] + 11*(21 + 66*c[2] + 38*c[3] + 10*c[4])*d[2] + 
	(-99 + 418*c[2] + 1034*c[3] + 602*c[4])*d[3] + (-429 + 110*c[2] + 602*c[3] 
	+ 1094*c[4])*d[4])) + c[2]*(143*(21 + 90*c[3] + 10*c[4])*d[0] + 52*c[3]*
	(209*d[1] - 7*(11*d[2] + d[3] - 9*d[4])) + 4*c[4]*(715*d[1] + 767*d[2] + 
	819*d[3] + 987*d[4]) + 286*(21*d[1] + 33*d[2]+5*(9*d[3]+d[4])))))/720720.;
		bf_mom[3] = (429*c[0]*c[0]*(35*d[0] + 14*(5*d[1] + d[2] - 3*d[3]) 
	- 26*d[4]) + 26*c[0]*(11*(3*(-35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4])
	*d[0] + 2*(3*(-35 + 42*c[1] + 14*c[2] - 6*c[3] - 26*c[4])*d[1] + (-21 + 
	42*c[1] + 66*c[2] + 90*c[3] + 10*c[4])*d[2] + (63 - 18*c[1] + 90*c[2] + 
	94*c[3] + 98*c[4])*d[3])) - 2*(858*c[1] - 11*(39 + 10*c[2] + 98*c[3]) - 
	1094*c[4])*d[4]) + 4*(9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] + 5577*c[4]*d[0]
	 + 14014*c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] + 2574*c[3]*d[1] + 13442*c[3]
	*c[3]*d[1] + 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1]
	 - 12870*c[3]*d[2] - 182*c[3]*c[3]*d[2] - 1430*c[4]*d[2] + 3276*c[3]*c[4]
	*d[2] + 1974*c[4]*c[4]*d[2] - 13442*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 
	14014*c[4]*d[3] - 8636*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 143*c[1]*
	c[1]*(9*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) - 34*d[4]) - 2*(2159*c[3]*c[3]
	 + 11*c[3]*(637 + 934*c[4]) + c[4]*(7111 + 4265*c[4]))*d[4] + 13*c[2]*c[2]*
	(363*d[0] + 726*d[1] + 286*d[2] - 154*d[3] + 118*d[4]) + 13*c[1]* (33*(-35
	 + 14*c[2] - 6*c[3] - 26*c[4])*d[0] + 2*(11*(-63 + 54*c[2] - 18*c[3] - 
	34*c[4])*d[1] + 11*(-21 + 66*c[2] + 38*c[3] + 10*c[4])*d[2] + (99 + 418*c[2]
	 + 1034*c[3] + 602*c[4])*d[3] + (429 + 110*c[2] + 602*c[3] + 1094*c[4])
	*d[4])) + c[2]*(143*(-21 + 90*c[3] + 10*c[4])*d[0] + 52*c[3]*(209*d[1] - 
	7*(11*d[2] + d[3] - 9*d[4])) + 4*c[4]*(715*d[1] + 767*d[2] + 819*d[3] + 
	987*d[4]) - 286*(21*d[1] + 33*d[2] + 5*(9*d[3] + d[4])))))/720720.;
		bf_mom[4] = (-13*(1155*c[0]*c[0] + 44* (-105 + 63*c[1]*c[1] - 
	42*c[1]*c[2] + 33*c[2]*c[2] - 18*(c[1] + 5*c[2])*c[3] + 47*c[3]*c[3]) + 
	88*(39*c[1] + 5*c[2] - 49*c[3])*c[4] + 2188*c[4]*c[4] - 132*c[0]*(35*c[1] 
	- 7*(c[2] + 3*c[3]) + 13*c[4]))*d[0] + 2*(13*(1155*c[0]*c[0] + 44*(-105 + 
	63*c[1]*c[1] + 33*c[2]*c[2] - 38*c[2]*c[3] + 47*c[3]*c[3] - 18*c[1]*(3*c[2] 
	+ c[3])) + 8*(187*c[1] + 55*c[2] - 301*c[3])*c[4] + 2188*c[4]*c[4] - 
	132*c[0]*(21*c[1] - 7*c[2] - 3*c[3] + 13*c[4]))*d[1] - (3003*c[0]*c[0] + 
	52* (11*(-21 + 27*c[1]*c[1] - 66*c[1]*c[2] + 13*c[2]*c[2]) + 22*(19*c[1] 
	+ 7*c[2])*c[3] - 7*c[3]*c[3]) - 572*c[0]*(21*c[1] - 33*c[2] + 45*c[3] - 
	5*c[4]) - 104*(55*c[1] - 59*c[2] + 63*c[3])*c[4] + 3948*c[4]*c[4])*d[2] - 
	(143*(63*c[0]*c[0] - 36*c[0]*(c[1] + 5*c[2]) + 4*(9*c[1]*c[1] + 38*c[1]*c[2]
	 + 7*(-9 + c[2]*c[2]))) + 52*(517*(c[0] - 2*c[1]) - 14*c[2])*c[3] + 27612*
	c[3]*c[3] - 4*(91*(77*c[0] - 86*c[1] + 18*c[2]) + 4318*c[3])*c[4] + 20548*
	c[4]*c[4])*d[3] + (5577*c[0]*c[0] - 52*c[0]* (429*c[1] + 55*c[2] - 539*c[3]
	 + 547*c[4]) + 4*(-5577 + 2431*c[1]*c[1] - 767*c[2]*c[2] + 1638*c[2]*c[3] 
	+ 2159*c[3]*c[3] - 2*(987*c[2] + 5137*c[3])*c[4] + 4265*c[4]*c[4] + 
	26*c[1]*(55*c[2] - 301*c[3] + 547*c[4])))*d[4]))/360360.;
		bf_mom[5] = (26*c[0]*(-33*(-35 + 42*c[2] - 6*c[4])*d[0] + 22*(-63 
	+ 114*c[2] - 62*c[4])*d[2] + 44*(3*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-27*d[1] 
	+ 53*d[3])) + 2*(99 - 682*c[2] + 1158*c[4])*d[4]) + 429*c[0]*c[0]*(35*d[0] 
	+ 6*(-7*d[2] + d[4])) + 4*(7579*c[3]*c[3]*d[0] + 1287*c[4]*d[0] + 7527*c[4]
	*c[4]*d[0] - 7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] - 
	8866*c[4]*d[2] - 8930*c[4]*c[4]*d[2] + 15158*c[3]*d[3] - 6420*c[3]*c[4]*d[3]
	 - 26*c[1]* (11*(-21 + 6*c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 238*c[4])
	*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*(21*d[0] - 
	6*d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] - c[4]*(7527 + 967*c[4]))*d[4] + 
	13*c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(63 + 62*c[4])
	*d[0] + 78*(209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) - 2*(4433 
	+ 8930*c[4])*d[4])))/ 180180.;
		bf_mom[6] = (-13*(1155*c[0]*c[0] + 924*c[0]*(5*c[1] + c[2]-3*c[3]) 
	+ 44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(7*c[2] - 3*c[3]) + 90*
	c[2]*c[3] + 47*c[3]*c[3]) - 44*(39*c[0] + 78*c[1] - 10*c[2] - 98*c[3])*c[4]
	 + 2188*c[4]*c[4])*d[0] + 3432*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) - 2*
	(429*c[0]*c[0]*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) + 52*c[0]*(11*(3*c[2]*
	(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*
	(-39*d[1] + 5*d[2] + 49*d[3])) + 33*c[1]*(21*d[1] + 7*d[2] - 3*d[3] - 13*
	d[4]) + (55*c[2] + 539*c[3] + 547*c[4])*d[4]) + 4*(6721*c[3]*c[3]*d[1] + 
	7826*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] - 91*c[3]*c[3]*d[2] + 1638*c[3]*
	c[4]*d[2] + 987*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] - 4318*c[3]*c[4]*d[3] 
	- 5137*c[4]*c[4]*d[3] + 26*c[2]*(c[4]*(55*d[1] + 59*d[2] + 63*d[3]) + c[3]*
	(209*d[1] - 7*(11*d[2] + d[3]))) + 143*c[1]*c[1]*(63*d[1] + 27*d[2] - 9*d[3]
	 - 17*d[4]) + 42*c[2]*(39*c[3] + 47*c[4])*d[4] - (2159*c[3]*c[3] + 10274*
	c[3]*c[4] + 4265*c[4]*c[4])* d[4] + 13*c[2]*c[2]* (363*d[1] + 143*d[2] - 
	77*d[3] + 59*d[4]) + 26*c[1]*(11*c[2]*(27*d[1] + 33*d[2] + 19*d[3] + 5*d[4])
	 + c[3]*(-99*d[1] + 209*d[2] + 517*d[3] + 301*d[4]) + c[4]*(-187*d[1] + 
	55*d[2] + 301*d[3] + 547*d[4])))))/360360.;
		bf_mom[7] = (26*c[0]*(-33*(35 + 42*c[2] - 6*c[4])*d[0] + 22*((63 + 
	114*c[2] - 62*c[4])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1] + 
	53*d[3])) - 2*(99 + 682*c[2] - 1158*c[4])*d[4]) + 429*c[0]*c[0]*(35*d[0] + 
	6*(-7*d[2] + d[4])) + 4*(7579*c[3]*c[3]*d[0] - 1287*c[4]*d[0] + 7527*c[4]*
	c[4]*d[0] + 7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] + 
	8866*c[4]*d[2] - 8930*c[4]*c[4]*d[2] - 15158*c[3]*d[3] - 6420*c[3]*c[4]*d[3]
	 - 26*c[1]* (11*(21 + 6*c[2] + 22*c[4])*d[1] - (297 + 286*c[2] + 238*c[4])
	*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*(21*d[0] - 
	6*d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] + (7527 - 967*c[4])*c[4])*d[4] + 13*
	c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(-63 + 62*c[4])
	*d[0] + 78*(-209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) + 2*
	(4433 - 8930*c[4])*d[4])))/ 180180.;
		bf_mom[8] = (-13*(33*(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) 
	+ 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 
	62*c[2])*c[4] + 2316*c[4]*c[4])*d[0] + 2*(1287*c[0]*c[0]*(7*d[2] - d[4]) - 
	52*c[0]*(627*c[2]*d[2] - 341*c[4]*d[2] + 33*c[1]*(7*d[1] - 9*d[3]) + c[3]*
	(-297*d[1] + 583*d[3]) - 341*c[2]*d[4] + 579*c[4]*d[4]) + 4*(-3094*c[3]*
	c[4]*d[1] - 9009*d[2] + 4173*c[3]*c[3]*d[2] + 4465*c[4]*c[4]*d[2] + 3210*
	c[3]*c[4]*d[3] + 13*c[2]*c[2]*(517*d[2] - 339*d[4]) + (1287 + 1605*c[3]*
	c[3] - 967*c[4]*c[4])*d[4] + 143*c[1]*c[1]*(3*d[2] + 11*d[4]) + 26*c[1]*
	(121*c[4]*d[1] - 143*c[3]*d[2] + 11*c[2]*(3*d[1] - 13*d[3]) - 119*c[4]*d[3]
	 - 119*c[3]*d[4]) + c[2]*(-3718*c[3]*d[1] - 8814*c[4]*d[2] + 8346*c[3]*d[3]
	 + 8930*c[4]*d[4]))))/90090.;
		break;
	case 32:
		bf_mom[0] = (429*c[0]*c[0]*(35*d[0] - 70*d[1] + 14*d[2] + 42*d[3] 
	- 26*d[4]) + 26*c[0]*(-33*(35 + 70*c[1] - 14*c[2] - 42*c[3] + 26*c[4])*d[0]
	 + 66*(35 + 42*c[1] - 14*c[2] - 6*c[3] + 26*c[4])*d[1] - 22*(21 + 42*c[1] 
	- 66*c[2] + 90*c[3] - 10*c[4])*d[2] - 22*(63 + 18*c[1] + 90*c[2] - 94*c[3] 
	+ 98*c[4])*d[3] + 2*(11*(39 + 78*c[1] + 10*c[2] - 98*c[3])+1094*c[4])*d[4])
	 + 4*(-9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] + 5577*c[4]*d[0] - 14014*c[3]*
	c[4]*d[0] + 7111*c[4]*c[4]*d[0] + 2574*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]
	*d[2] - 182*c[3]*c[3]*d[2] - 1430*c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 1974*
	c[4]*c[4]*d[2] - 13442*c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] 
	- 8636*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 143*c[1]*c[1]*(9*(7*d[0] + 
	2*(-7*d[1] + 3*d[2] + d[3])) - 34*d[4]) - 2*(2159*c[3]*c[3] - 11*c[3]*
	(637 + 934*c[4]) + c[4]*(7111 + 4265*c[4]))*d[4] + 13*c[2]*c[2]*(363*d[0] 
	- 726*d[1] + 286*d[2] + 154*d[3] + 118*d[4]) - 13*c[1]* (33*(-35 + 14*c[2] 
	+ 6*c[3] - 26*c[4])*d[0] + 2*(-11*(-63 + 54*c[2] + 18*c[3] - 34*c[4])*d[1] 
	+ 11*(-21 + 66*c[2] - 38*c[3] + 10*c[4])*d[2] - (99 + 418*c[2] - 1034*c[3] 
	+ 602*c[4])*d[3] + (429 + 110*c[2] - 602*c[3] + 1094*c[4])*d[4])) + c[2]*
	(-143*(21 + 90*c[3] - 10*c[4])*d[0] - 4*c[4]*(715*d[1] - 767*d[2] + 819*d[3]
	 - 987*d[4]) + 286*(21*d[1] - 33*d[2] + 45*d[3] - 5*d[4]) + 52*c[3]*
	(209*d[1] + 77*d[2] - 7*(d[3] + 9*d[4])))))/720720.;
		bf_mom[1] = (429*c[0]*c[0]*(35*d[0] + 14*(5*d[1] + d[2] - 3*d[3]) 
	- 26*d[4]) + 26*c[0]*(11*(3*(-35 + 70*c[1] + 14*c[2] - 42*c[3] - 26*c[4])
	*d[0] + 2*(3*(-35 + 42*c[1] + 14*c[2] - 6*c[3] - 26*c[4])*d[1] + (-21 + 
	42*c[1] + 66*c[2] + 90*c[3] + 10*c[4])*d[2] + (63 - 18*c[1] + 90*c[2] + 
	94*c[3] + 98*c[4])*d[3])) - 2*(858*c[1] - 11*(39 + 10*c[2] + 98*c[3]) - 
	1094*c[4])*d[4]) + 4*(9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] + 5577*c[4]*d[0]
	 + 14014*c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] + 2574*c[3]*d[1] + 13442*c[3]
	*c[3]*d[1] + 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] 
	- 12870*c[3]*d[2] - 182*c[3]*c[3]*d[2] - 1430*c[4]*d[2] + 3276*c[3]*c[4]*
	d[2] + 1974*c[4]*c[4]*d[2] - 13442*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 
	14014*c[4]*d[3] - 8636*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 143*c[1]*
	c[1]*(9*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) - 34*d[4]) - 2*(2159*c[3]*c[3]
	 + 11*c[3]*(637 + 934*c[4]) + c[4]*(7111 + 4265*c[4]))*d[4] + 13*c[2]*c[2]*
	(363*d[0] + 726*d[1] + 286*d[2] - 154*d[3] + 118*d[4]) + 13*c[1]* (33*(-35 
	+ 14*c[2] - 6*c[3] - 26*c[4])*d[0] + 2*(11*(-63 + 54*c[2] - 18*c[3] - 
	34*c[4])*d[1] + 11*(-21 + 66*c[2] + 38*c[3] + 10*c[4])*d[2] + (99 + 418*c[2]
	 + 1034*c[3] + 602*c[4])*d[3] + (429 + 110*c[2] + 602*c[3] + 1094*c[4])
	*d[4])) + c[2]*(143*(-21 + 90*c[3] + 10*c[4])*d[0] + 52*c[3]*(209*d[1] - 
	7*(11*d[2] + d[3] - 9*d[4])) + 4*c[4]*(715*d[1] + 767*d[2] + 819*d[3] + 
	987*d[4]) - 286*(21*d[1] + 33*d[2] + 5*(9*d[3] + d[4])))))/720720.;
		bf_mom[2] = (429*c[0]*c[0]*(35*d[0] + 14*(5*d[1] + d[2] - 3*d[3]) 
	- 26*d[4]) + 26*c[0]*(11*(3*(35 + 70*c[1] + 14*c[2] - 42*c[3]-26*c[4])*d[0]
	 + 2*(3*(35 + 42*c[1] + 14*c[2] - 6*c[3] - 26*c[4])*d[1] + (21 + 42*c[1] + 
	66*c[2] + 90*c[3] + 10*c[4])*d[2] + (-63 - 18*c[1] + 90*c[2] + 94*c[3] + 
	98*c[4])*d[3])) - 2*(11*(39 + 78*c[1] - 10*c[2] - 98*c[3]) - 1094*c[4])
	*d[4]) + 4*(-9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] - 5577*c[4]*d[0] + 
	14014*c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] - 2574*c[3]*d[1] + 13442*c[3]
	*c[3]*d[1] - 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] 
	+ 12870*c[3]*d[2] - 182*c[3]*c[3]*d[2] + 1430*c[4]*d[2] + 3276*c[3]*c[4]
	*d[2] + 1974*c[4]*c[4]*d[2] + 13442*c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 
	14014*c[4]*d[3] - 8636*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 143*c[1]*
	c[1]*(9*(7*d[0] + 14*d[1] + 6*d[2] - 2*d[3]) - 34*d[4]) - 2*(2159*c[3]*c[3]
	 + 11*c[3]*(-637 + 934*c[4]) + c[4]*(-7111 + 4265*c[4]))*d[4] + 13*c[2]*
	c[2]*(363*d[0] + 726*d[1] + 286*d[2] - 154*d[3] + 118*d[4]) + 13*c[1]* 
	(33*(35 + 14*c[2] - 6*c[3] - 26*c[4])*d[0] + 2*(11*(63 + 54*c[2] - 18*c[3] 
	- 34*c[4])*d[1] + 11*(21 + 66*c[2] + 38*c[3] + 10*c[4])*d[2] + (-99 + 
	418*c[2] + 1034*c[3] + 602*c[4])*d[3] + (-429 + 110*c[2] + 602*c[3] + 
	1094*c[4])*d[4])) + c[2]*(143*(21 + 90*c[3] + 10*c[4])*d[0] + 52*c[3]*
	(209*d[1] - 7*(11*d[2] + d[3] - 9*d[4])) + 4*c[4]*(715*d[1] + 767*d[2] + 
	819*d[3] + 987*d[4]) + 286*(21*d[1] + 33*d[2]+5*(9*d[3]+d[4])))))/720720.;
		bf_mom[3] = (429*c[0]*c[0]*(35*d[0] - 70*d[1] + 14*d[2] + 42*d[3] 
	- 26*d[4]) + 26*c[0]*(11*(-3*(70*c[1] - 7*(5 + 2*c[2] + 6*c[3]) + 26*c[4])
	*d[0] + 2*(3*(-35 + 42*c[1] - 14*c[2] - 6*c[3] + 26*c[4])*d[1] + (21 - 
	42*c[1] + 66*c[2] - 90*c[3] + 10*c[4])*d[2] + (63 - 18*c[1] - 90*c[2] + 
	94*c[3] - 98*c[4])*d[3])) + 2*(11*(-39 + 78*c[1] + 10*c[2] - 98*c[3]) + 
	1094*c[4])*d[4]) + 4*(9009*c[3]*d[0] + 6721*c[3]*c[3]*d[0] - 5577*c[4]*d[0]
	 - 14014*c[3]*c[4]*d[0] + 7111*c[4]*c[4]*d[0] - 2574*c[3]*d[1] - 13442*c[3]
	*c[3]*d[1] + 11154*c[4]*d[1] + 15652*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1]
	 - 12870*c[3]*d[2] - 182*c[3]*c[3]*d[2] + 1430*c[4]*d[2] - 3276*c[3]*c[4]
	*d[2] + 1974*c[4]*c[4]*d[2] + 13442*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 
	14014*c[4]*d[3] - 8636*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] + 143*c[1]*
	c[1]*(9*(7*d[0] + 2*(-7*d[1] + 3*d[2] + d[3])) - 34*d[4]) - 2*(2159*c[3]*
	c[3] - 11*c[3]*(-637 + 934*c[4]) + c[4]*(-7111 + 4265*c[4]))*d[4] + 13*
	c[2]*c[2]*(363*d[0] - 726*d[1] + 286*d[2] + 154*d[3] + 118*d[4]) - 13*c[1]*
	 (33*(35 + 14*c[2] + 6*c[3] - 26*c[4])*d[0] + 2*(-11*(63 + 54*c[2] + 
	18*c[3] - 34*c[4])*d[1] + 11*(21 + 66*c[2] - 38*c[3] + 10*c[4])*d[2] - 
	(-99 + 418*c[2] - 1034*c[3] + 602*c[4])*d[3] + (-429 + 110*c[2] - 602*c[3] 
	+ 1094*c[4])*d[4])) + c[2]*(-143*(-21 + 90*c[3] - 10*c[4])*d[0] - 4*c[4]*
	(715*d[1] - 767*d[2] + 819*d[3] - 987*d[4]) - 286*(21*d[1] - 33*d[2] + 
	45*d[3] - 5*d[4]) + 52*c[3]*(209*d[1] + 77*d[2] - 
	7*(d[3] + 9*d[4])))))/720720.;
		bf_mom[4] = (26*c[0]*(-33*(35 + 42*c[2] - 6*c[4])*d[0] + 22*((63 
	+ 114*c[2] - 62*c[4])*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + 2*c[3]*(-27*d[1] 
	+ 53*d[3])) - 2*(99 + 682*c[2] - 1158*c[4])*d[4]) + 429*c[0]*c[0]*(35*d[0] 
	+ 6*(-7*d[2] + d[4])) + 4*(7579*c[3]*c[3]*d[0] - 1287*c[4]*d[0] + 7527*
	c[4]*c[4]*d[0] + 7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2]
	 + 8866*c[4]*d[2] - 8930*c[4]*c[4]*d[2] - 15158*c[3]*d[3] - 6420*c[3]*c[4]
	*d[3] - 26*c[1]* (11*(21 + 6*c[2] + 22*c[4])*d[1] - (297 + 286*c[2] + 
	238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] + (7527 - 967*c[4])*c[4])
	*d[4] + 13*c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(-63 
	+ 62*c[4])*d[0] + 78*(-209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3])
	 + 2*(4433 - 8930*c[4])*d[4])))/ 180180.;
		bf_mom[5] = (-13*(1155*c[0]*c[0] + 924*c[0]*(5*c[1] + c[2] - 3*c[3])
	 + 44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] + 6*c[1]*(7*c[2] - 3*c[3]) + 90*
	c[2]*c[3] + 47*c[3]*c[3]) - 44*(39*c[0] + 78*c[1] - 10*c[2] - 98*c[3])*c[4]
	 + 2188*c[4]*c[4])*d[0] + 3432*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) - 
	2*(429*c[0]*c[0]*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) + 52*c[0]*(11*
	(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) 
	+ c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) + 33*c[1]*(21*d[1] + 7*d[2] - 3*d[3] 
	- 13*d[4]) + (55*c[2] + 539*c[3] + 547*c[4])*d[4]) + 4*(6721*c[3]*c[3]*d[1]
	 + 7826*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] - 91*c[3]*c[3]*d[2] + 1638*
	c[3]*c[4]*d[2] + 987*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] - 4318*c[3]*c[4]
	*d[3] - 5137*c[4]*c[4]*d[3] + 26*c[2]*(c[4]*(55*d[1] + 59*d[2] + 63*d[3]) 
	+ c[3]*(209*d[1] - 7*(11*d[2] + d[3]))) + 143*c[1]*c[1]*(63*d[1] + 27*d[2] 
	- 9*d[3] - 17*d[4]) + 42*c[2]*(39*c[3] + 47*c[4])*d[4] - (2159*c[3]*c[3] + 
	10274*c[3]*c[4] + 4265*c[4]*c[4])* d[4] + 13*c[2]*c[2]* (363*d[1] + 143*d[2]
	 - 77*d[3] + 59*d[4]) + 26*c[1]*(11*c[2]*(27*d[1] + 33*d[2] + 19*d[3] + 
	5*d[4]) + c[3]*(-99*d[1] + 209*d[2] + 517*d[3] + 301*d[4]) + c[4]*
	(-187*d[1] + 55*d[2] + 301*d[3] + 547*d[4])))))/360360.;
		bf_mom[6] = (26*c[0]*(-33*(-35 + 42*c[2] - 6*c[4])*d[0] + 22*(-63 
	+ 114*c[2] - 62*c[4])*d[2] + 44*(3*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-27*d[1] 
	+ 53*d[3])) + 2*(99 - 682*c[2] + 1158*c[4])*d[4]) + 429*c[0]*c[0]*(35*d[0] 
	+ 6*(-7*d[2] + d[4])) + 4*(7579*c[3]*c[3]*d[0] + 1287*c[4]*d[0] + 7527*c[4]
	*c[4]*d[0] - 7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] - 
	8866*c[4]*d[2] - 8930*c[4]*c[4]*d[2] + 15158*c[3]*d[3] - 6420*c[3]*c[4]*
	d[3] - 26*c[1]* (11*(-21 + 6*c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 
	238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] - c[4]*(7527 + 967*c[4]))
	*d[4] + 13*c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(63 + 
	62*c[4])*d[0] + 78*(209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) - 
	2*(4433 + 8930*c[4])*d[4])))/ 180180.;
		bf_mom[7] = (-13*(1155*c[0]*c[0] + 44* (-105 + 63*c[1]*c[1] - 
	42*c[1]*c[2] + 33*c[2]*c[2] - 18*(c[1] + 5*c[2])*c[3] + 47*c[3]*c[3]) + 
	88*(39*c[1] + 5*c[2] - 49*c[3])*c[4] + 2188*c[4]*c[4] - 132*c[0]*(35*c[1] 
	- 7*(c[2] + 3*c[3]) + 13*c[4]))*d[0] + 2*(13*(1155*c[0]*c[0] + 44*(-105 + 
	63*c[1]*c[1] + 33*c[2]*c[2] - 38*c[2]*c[3] + 47*c[3]*c[3] - 18*c[1]*(3*c[2]
	 + c[3])) + 8*(187*c[1] + 55*c[2] - 301*c[3])*c[4] + 2188*c[4]*c[4] - 
	132*c[0]*(21*c[1] - 7*c[2] - 3*c[3] + 13*c[4]))*d[1] - (3003*c[0]*c[0] + 
	52* (11*(-21 + 27*c[1]*c[1] - 66*c[1]*c[2] + 13*c[2]*c[2]) + 22*(19*c[1] 
	+ 7*c[2])*c[3] - 7*c[3]*c[3]) - 572*c[0]*(21*c[1] - 33*c[2] + 45*c[3] - 
	5*c[4]) - 104*(55*c[1] - 59*c[2] + 63*c[3])*c[4] + 3948*c[4]*c[4])*d[2] - 
	(143*(63*c[0]*c[0] - 36*c[0]*(c[1] + 5*c[2]) + 4*(9*c[1]*c[1] + 38*c[1]*c[2]
	 + 7*(-9 + c[2]*c[2]))) + 52*(517*(c[0] - 2*c[1]) - 14*c[2])*c[3] + 27612
	*c[3]*c[3] - 4*(91*(77*c[0] - 86*c[1] + 18*c[2]) + 4318*c[3])*c[4] + 
	20548*c[4]*c[4])*d[3] + (5577*c[0]*c[0] - 52*c[0]* (429*c[1] + 55*c[2] - 
	539*c[3] + 547*c[4]) + 4*(-5577 + 2431*c[1]*c[1] - 767*c[2]*c[2] + 1638*
	c[2]*c[3] + 2159*c[3]*c[3] - 2*(987*c[2] + 5137*c[3])*c[4] + 4265*c[4]*c[4]
	 + 26*c[1]*(55*c[2] - 301*c[3] + 547*c[4])))*d[4]))/360360.;
		bf_mom[8] = (-13*(33*(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) 
	+ 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 
	62*c[2])*c[4] + 2316*c[4]*c[4])*d[0] + 2*(1287*c[0]*c[0]*(7*d[2] - d[4]) - 
	52*c[0]*(627*c[2]*d[2] - 341*c[4]*d[2] + 33*c[1]*(7*d[1] - 9*d[3]) + c[3]*
	(-297*d[1] + 583*d[3]) - 341*c[2]*d[4] + 579*c[4]*d[4]) + 4*(-3094*c[3]*
	c[4]*d[1] - 9009*d[2] + 4173*c[3]*c[3]*d[2] + 4465*c[4]*c[4]*d[2] + 3210*
	c[3]*c[4]*d[3] + 13*c[2]*c[2]*(517*d[2] - 339*d[4]) + (1287 + 1605*c[3]*
	c[3] - 967*c[4]*c[4])*d[4] + 143*c[1]*c[1]*(3*d[2] + 11*d[4]) + 26*c[1]*
	(121*c[4]*d[1] - 143*c[3]*d[2] + 11*c[2]*(3*d[1] - 13*d[3]) - 119*c[4]*d[3]
	 - 119*c[3]*d[4]) + c[2]*(-3718*c[3]*d[1] - 8814*c[4]*d[2] + 8346*c[3]*d[3]
	 + 8930*c[4]*d[4]))))/90090.;
		break;
	case 100:
	case 102:
		bf_mom[0] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] - 
	2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(5 + 5*c[1] + 2*c[2] - 3*c[3])
	 + 10*c[4])* d[0] + 33*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])* d[1] -
	 22*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(63 + 72*c[1] 
	+ 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(165 + 429*c[1] - 286*c[2] - 539*
	c[3] + 1126*c[4])* d[4]) + 2*(-9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] + 4290*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 10296*c[3]*d[1] - 
	13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 14222*c[4]*
	c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] + 7436*c[4]*d[2] - 3276*
	c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*d[3] + 13806*c[3]*c[3]*
	d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 26*
	c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] - 11*(21 + 64*c[3] - 10*c[4])*d[1] 
	- 2*(-495 + 77*c[3] + 398*c[4])*d[2] + (-495 + 656*c[3] + 126*c[4])*d[3]) 
	+ 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*
	(-1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(3764*c[3]*c[3] - 11*c[3]*(637 + 
	934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*d[4] + 26*c[2]*c[2]* (495*d[0] - 
	363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]*(33*(-35 + 14*c[2] + 
	24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 
	11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(198 + 352*c[2] - 517*c[3] 
	+ 420*c[4])*d[3] + (429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4])))) + a*
	(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*d[1] + (5 - 6*b)*d[2] - 3*d[3]) 
	+ 2*(7 + 6*b)*d[4]) - 26*c[0]*(33*(7*(-5*(3 + 2*c[1] + 2*c[2]) + 2*b*(5 + 
	6*c[2]) + 6*c[3]) - 2*(7 + 6*b)*c[4])*d[0] + 22*(-3*(-35 + 14*(-5 + 2*b)*
	c[1] + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*c[4])*d[1] - (2*b*(63 + 114*c[2] - 
	62*c[4]) + 3*(-35 + 14*c[1] - 98*c[2] + 30*c[3] + 38*c[4]))*d[2] + (-63 + 
	18*(-7 + 6*b)*c[1] - 90*c[2] + 306*c[3] - 212*b*c[3] - 98*c[4])* d[3]) + 2*
	(2*b*(99 + 682*c[2] - 1158*c[4]) + 11*(21 + 78*c[1] - 114*c[2] - 98*c[3] + 
	310*c[4]))*d[4]) + 4*(9009*c[3]*d[0] - 21879*c[3]*c[3]*d[0] + 15158*b*c[3]*
	c[3]*d[0] - 3003*c[4]*d[0] - 2574*b*c[4]*d[0] + 14014*c[3]*c[4]*d[0] - 
	22165*c[4]*c[4]*d[0] + 15054*b*c[4]*c[4]*d[0] - 18018*c[3]*d[1] + 15444*b*
	c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] - 28028*c[3]*c[4]*d[1] 
	+ 12376*b*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] + 16874*
	c[3]*c[3]*d[2] - 16692*b*c[3]*c[3]*d[2] - 16302*c[4]*d[2] + 17732*b*c[4]*
	d[2] + 3276*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] 
	+ 43758*c[3]*d[3] - 30316*b*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 14014*c[4]*
	d[3] + 21476*c[3]*c[4]*d[3] - 12840*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3]
	 - 2*((-5369 + 3210*b)*c[3]*c[3] + 11*c[3]*(637 + 934*c[4]) - c[4]*(22165 
	+ 2331*c[4] + 2*b*(-7527 + 967*c[4])))*d[4] + 13*c[2]*c[2]* (33*(-49 + 38*b)
	*d[0] + 22*(33*d[1] + (81 - 94*b)*d[2] - 7*d[3]) + 2*(-737 + 678*b)*d[4]) + 
	c[2]*(-143*(-3*(-35 + 42*b + 30*c[3]) + 2*(-57 + 62*b)*c[4])*d[0] + 52*
	(-627*b*d[2] + c[4]*(55*d[1] + (-737 + 678*b)*d[2] + 63*d[3]) + c[3]*(11*
	(-45 + 26*b)*d[1] - 77*d[2] + (649 - 642*b)*d[3])) + 4*(819*c[3] + b*(4433 
	- 8930*c[4]) + 7943*c[4])*d[4] - 858*(7*d[1] - 49*d[2] + 15*d[3] + 19*d[4]))
	 + 143*c[1]*c[1]* (21*(-5 + 2*b)*d[0] - 2*(-63*d[1] + 3*(7 + 2*b)*d[2] + 
	9*d[3] + (-39 + 22*b)*d[4])) + 13*c[1]*(33*(-35 + 14*c[2] + 6*(7 - 6*b)*
	c[3] - 26*c[4])* d[0] + 2*(-11* (3*(-35 + 14*c[2] + 6*c[3] - 26*c[4]) + 2*
	b*(21 + 6*c[2] + 22*c[4]))*d[1] + 11*(-21 + 66*c[2] + (-90 + 52*b)*c[3] + 
	10*c[4])*d[2] + (-11*(63 + 90*c[2] - 94*c[3] + 98*c[4]) + b*(594 + 572*c[2]
	 + 476*c[4]))*d[3] + (429 + 110*c[2] + 14*(-77 + 34*b)*c[3] + 1094*c[4])
	*d[4])))) + b*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*(-1 + b)*d[1] + 
	(5 - 2*b)*d[2] - 3*(-1 + b)*d[3]) - 2*(-7 + 10*b)*d[4]) + 26*c[0]*(33*(7*
	(15 - 10*c[1] + 10*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3]) - 
	2*(-7 + 10*b)*c[4])*d[0] + 22*(3*(35 - 70*c[1] - 14*c[2] + 42*c[3] + b*
	(-35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] + (-3*(35 + 
	14*c[1] + 98*c[2] + 30*c[3] - 38*c[4]) + 2*b*(21 + 21*c[1] + 90*c[2] + 45*
	c[3] - 26*c[4]))*d[2] + (-63 + 63*b + 126*c[1] - 72*b*c[1] - 90*c[2] + 90*
	b*c[2] - 306*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])*d[3]) - 2*(11*(21 - 78*
	c[1] - 114*c[2] + b*(-30 + 78*c[1] + 52*c[2] - 98*c[3]) + 98*c[3]) - 2*
	(-1705 + 1126*b)*c[4])*d[4]) + 4*(-9009*c[3]*d[0] + 9009*b*c[3]*d[0] - 
	21879*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] + 4290*b*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 22165*c[4]*c[4]
	*d[0] + 14638*b*c[4]*c[4]*d[0] - 18018*c[3]*d[1] + 10296*b*c[3]*d[1] - 
	13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 11154*b*
	c[4]*d[1] - 28028*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 14222*c[4]*c[4]
	*d[1] + 14222*b*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 12870*b*c[3]*d[2] + 
	16874*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] - 16302*c[4]*d[2] + 7436*b*
	c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*
	d[2] - 6956*b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] - 28600*b*c[3]*d[3] + 13806*
	c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 14014*b*c[4]*
	d[3] + 21476*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3]
	 - 10274*b*c[4]*c[4]*d[3] - 2*((-5369 + 3764*b)*c[3]*c[3] + 11*(-1 + b)*c[3]
	*(637 + 934*c[4]) + c[4]*(-22165 - 2331*c[4] + 2*b*(7319 + 1649*c[4])))*
	d[4] + 143*c[1]*c[1]* (3*(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*d[1] + (-7 + 
	8*b)*d[2] - 3*(-1 + b)*d[3])) - 2*(-39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*
	(-49 + 30*b)*d[0] + 22*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] - 7*(-1 + b)*
	d[3]) + 2*(-737 + 398*b)*d[4]) + c[2]*(143*(b*(42 + 90*c[3] - 52*c[4]) - 
	3*(35 + 30*c[3] - 38*c[4]))*d[0] + 286*(21 - 21*b - 90*c[3] + 64*b*c[3] + 
	10*(-1 + b)*c[4])* d[1] - 26*(-77*(21 + 2*c[3]) + 2*b*(495 + 77*c[3] - 398*
	c[4]) + 1474*c[4])*d[2] - 26*(-495 + 495*b - 1298*c[3] + 656*b*c[3] - 126*
	(-1 + b)*c[4])*d[3] + 2*(143*(-57 + 26*b) + 1638*(-1 + b)*c[3] - 94*(-169 
	+ 74*b)*c[4])*d[4]) + 13*c[1]*(33*(35 - 14*c[2] + 42*c[3] + b*(-35 + 14*
	c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] + 462*(5*d[1] + d[2] - 3*d[3]) - 
	858*d[4] + 2*(-22*(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] + 
	45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) - 2*(55*c[2] + 
	539*c[3] + 547*c[4])*d[4] + b*(22*(-42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] 
	+ 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(198 + 352*c[2] + 517*c[3]
	 + 420*c[4])*d[3] + (429 + 110*c[2] + 840*c[3] + 1094*c[4])
	*d[4])))))))/ 2882880.;
		bf_mom[1] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] 
	- 2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(5 + 5*c[1] + 2*c[2] - 3*
	c[3]) + 10*c[4])* d[0] + 33*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])*
	 d[1] - 22*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(63 + 
	72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(165 + 429*c[1] - 286*
	c[2] - 539*c[3] + 1126*c[4])* d[4]) + 2*(-9009*c[3]*d[0] + 14300*c[3]*c[3]*
	d[0] + 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 10296*
	c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] 
	- 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] + 7436*c[4]
	*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*d[3] + 13806
	*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]
	*d[3] - 26*c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] - 11*(21 + 64*c[3] - 
	10*c[4])*d[1] - 2*(-495 + 77*c[3] + 398*c[4])*d[2] + (-495 + 656*c[3] + 
	126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 
	28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(3764*c[3]*c[3] 
	- 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*d[4] + 26*c[2]*c[2]*
	 (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]*(33*(-35 
	+ 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 24*c[2] + 9*c[3] - 28*
	c[4])*d[1] + 11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(198 + 352*c[2]
	 - 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] - 840*c[3] + 1094*c[4])*
	d[4])))) + a*(429*c[0]*c[0]*(35*(3 + 2*b)*d[0] - 70*d[1] - 14*(5 + 6*b)*d[2]
	 + 42*d[3] + 2*(-7 + 6*b)*d[4]) - 26*c[0]*(33*(2*b*(35 + 42*c[2] - 6*c[4]) 
	+ 7*(15 + 10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]))*d[0] + 22*(-3*(35 + 14*(5 
	+ 2*b)*c[1] - 14*c[2] - 6*(7 + 6*b)*c[3] + 26*c[4])*d[1] + (-3*(35 - 14*c[1]
	 + 98*c[2] + b*(42 + 76*c[2]) - 30*c[3]) + 2*(57 + 62*b)*c[4])*d[2] + (63 
	+ 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 106*b)*c[3] + 98*c[4])*d[3]) + 2*
	(11*(-21 - 78*c[1] + 114*c[2] + 2*b*(9 + 62*c[2]) + 98*c[3]) - 2*(1705 + 
	1158*b)*c[4])*d[4]) + 4*(-9009*c[3]*d[0] + 21879*c[3]*c[3]*d[0] + 15158*b*
	c[3]*c[3]*d[0] + 3003*c[4]*d[0] - 2574*b*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 
	22165*c[4]*c[4]*d[0] + 15054*b*c[4]*c[4]*d[0] + 18018*c[3]*d[1] + 15444*b*
	c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 28028*c[3]*c[4]*d[1] 
	+ 12376*b*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 16874*
	c[3]*c[3]*d[2] - 16692*b*c[3]*c[3]*d[2] + 16302*c[4]*d[2] + 17732*b*c[4]*
	d[2] - 3276*c[3]*c[4]*d[2] - 15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] 
	- 43758*c[3]*d[3] - 30316*b*c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*
	d[3] - 21476*c[3]*c[4]*d[3] - 12840*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3]
	 - 2*((5369 + 3210*b)*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + c[4]*(22165 + 
	2*b*(7527 - 967*c[4]) + 2331*c[4]))*d[4] + 13*c[2]*c[2]* (33*(49 + 38*b)*
	d[0] - 22*(33*d[1] + (81 + 94*b)*d[2] - 7*d[3]) + 2*(737 + 678*b)*d[4]) + 
	c[2]*(-143*(-21*(5 + 6*b) + 90*c[3] + 2*(57 + 62*b)*c[4])* d[0] + 26*(11*
	(21 + (90 + 52*b)*c[3] - 10*c[4])*d[1] + (-33*(49 + 38*b) + 154*c[3] + 2*
	(737 + 678*b)*c[4])* d[2] + 495*d[3] - 2*((649 + 642*b)*c[3] + 63*c[4])*
	d[3]) - 2*(-13*(627 + 682*b - 126*c[3]) + 94*(169 + 190*b)*c[4])* d[4]) + 
	143*c[1]*c[1]* (21*(5 + 2*b)*d[0] - 2*(63*d[1] + 3*(-7 + 2*b)*d[2] - 9*d[3]
	 + (39 + 22*b)*d[4])) + 13*c[1]*(-33*(-35 + 14*c[2] + 6*(7 + 6*b)*c[3] - 
	26*c[4])* d[0] + 2*(-11* (105 - 42*c[2] - 18*c[3] + 78*c[4] + 2*b*(21 + 6*
	c[2] + 22*c[4]))*d[1] - 11*(-21 + 66*c[2] - 2*(45 + 26*b)*c[3] + 10*c[4])*
	 d[2] + (11* (63 + 54*b + 90*c[2] + 52*b*c[2] - 94*c[3]) + 14*(77 + 34*b)*
	c[4])*d[3] - (429 + 110*c[2] - 14*(77 + 34*b)*c[3] + 1094*c[4])*d[4])))) + 
	b*(429*c[0]*c[0]*(35*(3 + 2*b)*d[0] + 70*(1 + b)*d[1] - 14*(5 + 2*b)*d[2] - 
	42*(1 + b)*d[3] - 2*(7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(-15 + 10*c[1] - 10*
	c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3]) - 2*(7 + 10*b)*c[4])*
	d[0] + 22*(3*(-35 - 35*b + 70*c[1] + 56*b*c[1] + 14*c[2] + 14*b*c[2] - 42*
	c[3] - 24*b*c[3] - 26*(1 + b)*c[4])* d[1] + (3*(35 + 14*b + 14*c[1] + 14*b*
	c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3]) - 2*(57 + 26*b)*c[4])*d[2] + 
	(63 + 63*b - 126*c[1] - 72*b*c[1] + 90*c[2] + 90*b*c[2] + 306*c[3] + 200*b*
	c[3] + 98*(1 + b)*c[4])* d[3]) - 2*(11* (-21 - 30*b + 78*c[1] + 78*b*c[1] 
	+ 114*c[2] + 52*b*c[2] - 98*(1 + b)*c[3]) - 2*(1705 + 1126*b)*c[4])* d[4]) 
	+ 4*(9009*c[3]*d[0] + 9009*b*c[3]*d[0] + 21879*c[3]*c[3]*d[0] + 14300*b*
	c[3]*c[3]*d[0] + 3003*c[4]*d[0] + 4290*b*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 
	14014*b*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] + 
	18018*c[3]*d[1] + 10296*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 13442*b*c[3]*
	c[3]*d[1] + 11154*c[4]*d[1] + 11154*b*c[4]*d[1] + 28028*c[3]*c[4]*d[1] + 
	21840*b*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] - 
	12870*c[3]*d[2] - 12870*b*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 8528*b*c[3]*
	c[3]*d[2] + 16302*c[4]*d[2] + 7436*b*c[4]*d[2] + 3276*c[3]*c[4]*d[2] + 
	3276*b*c[3]*c[4]*d[2] - 15886*c[4]*c[4]*d[2] - 6956*b*c[4]*c[4]*d[2] - 
	43758*c[3]*d[3] - 28600*b*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 13806*b*c[3]*
	c[3]*d[3] - 14014*c[4]*d[3] - 14014*b*c[4]*d[3] - 21476*c[3]*c[4]*d[3] - 
	15056*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] - 10274*b*c[4]*c[4]*d[3] - 
	2*((5369 + 3764*b)*c[3]*c[3] + 11*(1 + b)*c[3]*(637 + 934*c[4]) + c[4]*
	(22165 + 2331*c[4] + 2*b*(7319 + 1649*c[4])))*d[4] + 143*c[1]*c[1]* (21*
	(5 + 4*b)*d[0] + 126*(1 + b)*d[1] + 6*(7 + 8*b)*d[2] - 18*(1 + b)*d[3] - 2*
	(39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*(49 + 30*b)*d[0] + 726*(1 + b)*d[1] - 
	22*(81 + 34*b)*d[2] - 154*(1 + b)*d[3] + 2*(737 + 398*b)*d[4]) + c[2]*(143*
	(b*(42 + 90*c[3] - 52*c[4]) + 3*(35 + 30*c[3] - 38*c[4]))*d[0] + 286*(-21 
	- 21*b + 90*c[3] + 64*b*c[3] + 10*(1 + b)*c[4])* d[1] - 26*(33*(49 + 30*b) 
	+ 154*(1 + b)*c[3] - 2*(737 + 398*b)*c[4])*d[2] - 26*(495 + 495*b + 1298*
	c[3] + 656*b*c[3] - 126*(1 + b)*c[4])*d[3] + 2*(143*(57 + 26*b) + 1638*(1 
	+ b)*c[3] - 94*(169 + 74*b)*c[4])*d[4]) + 13*c[1]*(33*(-35 - 35*b + 14*c[2]
	 + 14*b*c[2] - 42*c[3] - 24*b*c[3] - 26*(1 + b)*c[4])*d[0] - 462*(5*d[1] + 
	d[2] - 3*d[3]) + 858*d[4] + 2*(22*(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + 
	c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) + 
	2*(55*c[2] + 539*c[3] + 547*c[4])*d[4] + b*(22*(-42 + 24*c[2] - 9*c[3] - 
	28*c[4])*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(198 + 352*
	c[2] + 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] + 840*c[3] + 1094*c[4])
	*d[4])))))))/ 2882880.;
		bf_mom[2] = -((a - b)*(a*(429*c[0]*c[0]* (35*(3 + 2*b)*d[0] - 70*
	d[1] - 14*(5 + 6*b)*d[2] + 42*d[3] + 2*(-7 + 6*b)*d[4]) - 26*c[0]*(33*(2*b*
	(-35 + 42*c[2] - 6*c[4]) + 7*(-15 + 10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]))*
	d[0] + 22*(-3*(-35 + 14*(5 + 2*b)*c[1] - 14*c[2] - 6*(7 + 6*b)*c[3] + 26*
	c[4])*d[1] + (3*(35 + 14*c[1] + b*(42 - 76*c[2]) - 98*c[2] + 30*c[3]) + 2*
	(57 + 62*b)*c[4])*d[2] + (-63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153+106*b)
	*c[3] + 98*c[4])*d[3]) + 2*(11*(21 - 78*c[1] + 114*c[2] + 2*b*(-9 + 62*c[2])
	 + 98*c[3]) - 2*(1705 + 1158*b)*c[4])*d[4]) + 4*(9009*c[3]*d[0] + 21879*c[3]
	*c[3]*d[0] + 15158*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] + 2574*b*c[4]*d[0] - 
	14014*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 15054*b*c[4]*c[4]*d[0] - 
	18018*c[3]*d[1] - 15444*b*c[3]*d[1] - 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1]
	 + 28028*c[3]*c[4]*d[1] + 12376*b*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 
	12870*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 16692*b*c[3]*c[3]*d[2] - 16302*c[4]
	*d[2] - 17732*b*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 15886*c[4]*c[4]*d[2] - 
	17860*b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] + 30316*b*c[3]*d[3] + 13806*c[3]*
	c[3]*d[3] - 14014*c[4]*d[3] - 21476*c[3]*c[4]*d[3] - 12840*b*c[3]*c[4]*d[3]
	 + 10274*c[4]*c[4]*d[3] - 2*((5369 + 3210*b)*c[3]*c[3] - 11*c[3]*(-637 + 
	934*c[4]) - c[4]*(22165 - 2331*c[4] + 2*b*(7527 + 967*c[4])))*d[4] + 13*
	c[2]*c[2]* (33*(49 + 38*b)*d[0] - 22*(33*d[1] + (81 + 94*b)*d[2] - 7*d[3]) 
	+ 2*(737 + 678*b)*d[4]) + c[2]*(-143*(3*(35 + 42*b + 30*c[3]) + 2*(57+62*b)
	*c[4])* d[0] + 26*(11*(-21 + (90 + 52*b)*c[3] - 10*c[4])*d[1] + (11*(147 + 
	114*b + 14*c[3]) + 2*(737 + 678*b)*c[4])* d[2] - (2*(649 + 642*b)*c[3] + 9*
	(55 + 14*c[4]))*d[3]) - 2*(13*(627 + 682*b + 126*c[3]) + 94*(169 + 190*b)*
	c[4])*d[4]) + 143*c[1]*c[1]* (21*(5 + 2*b)*d[0] - 2*(63*d[1] + 3*(-7 + 2*b)
	*d[2] - 9*d[3] + (39 + 22*b)*d[4])) + 13*c[1]*(-33*(35 + 14*c[2] + 6*(7 + 
	6*b)*c[3] - 26*c[4])* d[0] + 2*(-11* (-3*(35 + 14*c[2] + 6*c[3] - 26*c[4]) 
	+ 2*b*(-21 + 6*c[2] + 22*c[4]))*d[1] - 11*(21 + 66*c[2] - 2*(45 + 26*b)*c[3]
	 + 10*c[4])*d[2] + (11*(-63 - 54*b + 90*c[2] + 52*b*c[2] - 94*c[3]) + 14*
	(77 + 34*b)*c[4])*d[3] + 429*d[4] - 2*(55*c[2] - 7*(77 + 34*b)*c[3] + 547*
	c[4])*d[4])))) + 2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] - 2*d[2]+3*d[3])
	 - 10*d[4]) + 26*c[0]*(-33*(7*(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + 10*c[4])* 
	d[0] + 33*(-35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])* d[1] - 22*(21 + 
	21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(-63 + 72*c[1] + 90*c[2] 
	- 200*c[3] + 98*c[4])* d[3] + 2*(429*c[1] - 11*(15 + 26*c[2] + 49*c[3]) + 
	1126*c[4])*d[4]) + 2*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0]
	 - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*d[1] - 13442*
	c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*
	d[1] - 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*d[2] - 3276*c[3]*
	c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 
	14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*
	(11*(21 + 45*c[3] + 26*c[4])*d[0] + (231 - 704*c[3] + 110*c[4])*d[1] - 2*
	(495 + 77*c[3] + 398*c[4])*d[2] + (495 + 656*c[3] + 126*c[4])*d[3]) + 286*
	c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(1859 
	+ 819*c[3] + 3478*c[4])*d[4] - 2*(c[3]*(7007 + 3764*c[3]) - 2*(7319 + 5137*
	c[3])*c[4] + 3298*c[4]*c[4])*d[4] + 26*c[2]*c[2]* (495*d[0] - 363*d[1] - 
	374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]*(33*(35 + 14*c[2] + 24*c[3] - 26*
	c[4])*d[0] + 2*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*
	c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3] + 420*c[4])
	*d[3] + (-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4])))) + b*(429*c[0]*
	c[0]*(35*(3 + 2*b)*d[0] + 70*(1 + b)*d[1] - 14*(5 + 2*b)*d[2] - 42*(1 + b)
	*d[3] - 2*(7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(15 + 10*c[1] - 10*c[2] + 2*b*
	(5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3]) - 2*(7 + 10*b)*c[4])*d[0] + 22*(3*
	(35 + 35*b + 70*c[1] + 56*b*c[1] + 14*c[2] + 14*b*c[2] - 42*c[3] - 24*b*c[3]
	 - 26*(1 + b)*c[4])* d[1] + (3*(-35 - 14*b + 14*c[1] + 14*b*c[1] + 98*c[2] 
	+ 60*b*c[2] + 30*(1 + b)*c[3]) - 2*(57 + 26*b)*c[4])*d[2] + (-63 - 63*b - 
	126*c[1] - 72*b*c[1] + 90*c[2] + 90*b*c[2] + 306*c[3] + 200*b*c[3] + 98*(1 
	+ b)*c[4])* d[3]) - 2*(11* (21 + 30*b + 78*c[1] + 78*b*c[1] + 114*c[2] + 
	52*b*c[2] - 98*(1 + b)*c[3]) - 2*(1705 + 1126*b)*c[4])* d[4]) + 4*(-9009*
	c[3]*d[0] - 9009*b*c[3]*d[0] + 21879*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0]
	 - 3003*c[4]*d[0] - 4290*b*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*
	c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] - 18018*c[3]*d[1]
	 - 10296*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] - 11154*b*c[4]*d[1] + 28028*c[3]*c[4]*d[1] + 21840*b*c[3]*
	c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] + 12870*c[3]*d[2]
	 + 12870*b*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] - 16302
	*c[4]*d[2] - 7436*b*c[4]*d[2] + 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2]
	 - 15886*c[4]*c[4]*d[2] - 6956*b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] + 28600*b*
	c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] + 14014*c[4]*d[3]
	 + 14014*b*c[4]*d[3] - 21476*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] - 
	10274*c[4]*c[4]*d[3] - 10274*b*c[4]*c[4]*d[3] - 2*((5369 + 3764*b)*c[3]*c[3]
	 + 11*(1 + b)*c[3]*(-637 + 934*c[4]) + c[4]*(-22165 + 2331*c[4] + 2*b*(-7319
	 + 1649*c[4])))*d[4] + 143*c[1]*c[1]* (21*(5 + 4*b)*d[0] + 126*(1 + b)*d[1]
	 + 6*(7 + 8*b)*d[2] - 18*(1 + b)*d[3] - 2*(39 + 28*b)*d[4]) + 13*c[2]*c[2]*
	 (33*(49 + 30*b)*d[0] + 726*(1 + b)*d[1] - 22*(81 + 34*b)*d[2] - 154*(1 + b)
	*d[3] + 2*(737 + 398*b)*d[4]) + c[2]*(143*(b*(-42 + 90*c[3] - 52*c[4]) + 3*
	(-35 + 30*c[3] - 38*c[4]))*d[0] + 286*(21 + 21*b + 90*c[3] + 64*b*c[3] + 
	10*(1 + b)*c[4])* d[1] - 26*(-33*(49 + 30*b) + 154*(1 + b)*c[3] - 2*(737 +
	 398*b)*c[4])*d[2] - 26*(-495 - 495*b + 1298*c[3] + 656*b*c[3] - 126*(1+b)*
	c[4])*d[3] + 2*(-143*(57 + 26*b) + 1638*(1 + b)*c[3] - 94*(169+74*b)*c[4])*
	d[4]) + 13*c[1]*(33*(35 + 35*b + 14*c[2] + 14*b*c[2] - 42*c[3] - 24*b*c[3] 
	- 26*(1 + b)*c[4])*d[0] + 462*(5*d[1] + d[2] - 3*d[3]) - 858*d[4] + 2*(22*
	(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) + 
	c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) + 2*(55*c[2] + 539*c[3] + 547*c[4])*d[4]
	 + b*(22*(42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3]
	 + 10*c[4])*d[2] + 2*(-198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 +
	 110*c[2] + 840*c[3] + 1094*c[4])*d[4])))))))/ 2882880.;
		bf_mom[3] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] - 
	2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(-5 + 5*c[1] + 2*c[2]-3*c[3])
	 + 10*c[4])* d[0] + 33*(-35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])* d[1] 
	- 22*(21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(-63 + 72*c[1]
	 + 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(429*c[1] - 11*(15 + 26*c[2] + 
	49*c[3]) + 1126*c[4])*d[4]) + 2*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 
	4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*
	d[1] - 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 
	14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*
	d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 28600*c[3]*d[3] + 13806*
	c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*
	d[3] - 26*c[2]*(11*(21 + 45*c[3] + 26*c[4])*d[0] + (231 - 704*c[3] + 110*
	c[4])*d[1] - 2*(495 + 77*c[3] + 398*c[4])*d[2] + (495 + 656*c[3] + 126*c[4])
	*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 
	4*c[2]*(1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(c[3]*(7007 + 3764*c[3]) - 2*
	(7319 + 5137*c[3])*c[4] + 3298*c[4]*c[4])*d[4] + 26*c[2]*c[2]* (495*d[0] - 
	363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]*(33*(35 + 14*c[2] + 24*
	c[3] - 26*c[4])*d[0] + 2*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*
	(21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3] + 
	420*c[4])*d[3] + (-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4])))) + a*
	(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*d[1] + (5 - 6*b)*d[2] - 3*d[3])
	 + 2*(7 + 6*b)*d[4]) - 26*c[0]*(11*(3*(7* (15 - 10*c[1] - 10*c[2] + 2*b*(-5
	 + 6*c[2]) + 6*c[3]) - 2*(7 + 6*b)*c[4])*d[0] - 6*(35 + 14*(-5 + 2*b)*c[1]
	 + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*c[4])*d[1] - 2*(2*b*(-63 + 114*c[2] - 
	62*c[4]) + 3*(35 + 14*c[1] - 98*c[2] + 30*c[3] + 38*c[4]))*d[2] + 2*(63 + 
	18*(-7 + 6*b)*c[1] - 90*c[2] + (306 - 212*b)*c[3] - 98*c[4])* d[3]) + 2*(2*
	b*(-99 + 682*c[2] - 1158*c[4]) + 11*(-21 + 78*c[1] - 114*c[2] - 98*c[3] + 
	310*c[4]))*d[4]) + 4*(-9009*c[3]*d[0] - 21879*c[3]*c[3]*d[0] + 15158*b*c[3]*
	c[3]*d[0] + 3003*c[4]*d[0] + 2574*b*c[4]*d[0] + 14014*c[3]*c[4]*d[0] - 
	22165*c[4]*c[4]*d[0] + 15054*b*c[4]*c[4]*d[0] + 18018*c[3]*d[1] - 15444*b*
	c[3]*d[1] + 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] - 28028*c[3]*c[4]*d[1] 
	+ 12376*b*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] + 16874*
	c[3]*c[3]*d[2] - 16692*b*c[3]*c[3]*d[2] + 16302*c[4]*d[2] - 17732*b*c[4]*
	d[2] + 3276*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] 
	- 43758*c[3]*d[3] + 30316*b*c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 14014*c[4]*
	d[3] + 21476*c[3]*c[4]*d[3] - 12840*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3]
	 - 2*((-5369 + 3210*b)*c[3]*c[3] + 11*c[3]*(-637 + 934*c[4]) - c[4]*(-22165
	 + 2331*c[4] + 2*b*(7527 + 967*c[4])))*d[4] + 13*c[2]*c[2]* (33*(-49 + 38*b)
	*d[0] + 22*(33*d[1] + (81 - 94*b)*d[2] - 7*d[3]) + 2*(-737 + 678*b)*d[4]) 
	+ c[2]*(143*(3*(35 + 30*c[3] + 38*c[4]) - 2*b*(63 + 62*c[4]))* d[0] + 26*
	(11*(21 + (-90 + 52*b)*c[3] + 10*c[4])*d[1] + (-11*(147 + 14*c[3]+134*c[4])
	 + 6*b*(209 + 226*c[4]))*d[2] + (2*(649 - 642*b)*c[3] + 9*(55 + 14*c[4]))
	*d[3]) + 2*(13*(627 - 682*b + 126*c[3]) - 94*(-169 + 190*b)*c[4])* d[4]) 
	+ 143*c[1]*c[1]* (21*(-5 + 2*b)*d[0] - 2*(-63*d[1] + 3*(7 + 2*b)*d[2] + 
	9*d[3] + (-39 + 22*b)*d[4])) + 13*c[1]*(33*(35 + 14*c[2] + 6*(7 - 6*b)*c[3]
	 - 26*c[4])* d[0] + 2*(-11* (3*(35 + 14*c[2] + 6*c[3] - 26*c[4]) + 2*b*
	(-21 + 6*c[2] + 22*c[4]))*d[1] + 11*(21 + 66*c[2] + (-90 + 52*b)*c[3] + 
	10*c[4])*d[2] + (11*(63 - 90*c[2] + b*(-54 + 52*c[2]) + 94*c[3]) + 14*(-77 
	+ 34*b)*c[4])*d[3] + (-429 + 110*c[2] + 14*(-77 + 34*b)*c[3] + 1094*c[4])*
	 d[4])))) + b*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*(-1 + b)*d[1] + 
	(5 - 2*b)*d[2] - 3*(-1 + b)*d[3]) - 2*(-7 + 10*b)*d[4]) + 26*c[0]*(33*(7*
	(-15 - 10*c[1] + 10*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3]) - 
	2*(-7 + 10*b)*c[4])*d[0] + 22*(3*(-7*(5 + 10*c[1] + 2*c[2] - 6*c[3]) + b*
	(35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] + (-3*(-35 + 
	14*c[1] + 98*c[2] + 30*c[3] - 38*c[4]) + 2*b*(-21 + 21*c[1] + 90*c[2] + 45*
	c[3] - 26*c[4]))* d[2] + (63 - 63*b + 126*c[1] - 72*b*c[1] - 90*c[2] + 90*
	b*c[2] - 306*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])* d[3]) - 2*(11* (-3*(7 
	+ 26*c[1] + 38*c[2]) + b*(30 + 78*c[1] + 52*c[2] - 98*c[3]) + 98*c[3]) - 
	2*(-1705 + 1126*b)*c[4])*d[4]) + 4*(9009*c[3]*d[0] - 9009*b*c[3]*d[0] - 
	21879*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] + 3003*c[4]*d[0] - 4290*b*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 22165*c[4]*c[4]
	*d[0] + 14638*b*c[4]*c[4]*d[0] + 18018*c[3]*d[1] - 10296*b*c[3]*d[1] - 
	13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] + 11154*c[4]*d[1] - 11154*b*
	c[4]*d[1] - 28028*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 14222*c[4]*
	c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] - 12870*c[3]*d[2] + 12870*b*c[3]*d[2] +
	 16874*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] + 16302*c[4]*d[2] - 7436*b*
	c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*
	d[2] - 6956*b*c[4]*c[4]*d[2] - 43758*c[3]*d[3] + 28600*b*c[3]*d[3] + 13806*
	c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] - 14014*c[4]*d[3] + 14014*b*c[4]*
	d[3] + 21476*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3]
	 - 10274*b*c[4]*c[4]*d[3] - 2*((-5369 + 3764*b)*c[3]*c[3] + 11*(-1 + b)*
	c[3]*(-637 + 934*c[4]) + c[4]*(22165 - 2331*c[4] + 2*b*(-7319 + 1649*c[4])))
	*d[4] + 143*c[1]*c[1]* (3*(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*d[1] + (-7 + 
	8*b)*d[2] - 3*(-1 + b)*d[3])) - 2*(-39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*
	(-49 + 30*b)*d[0] + 22*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] - 7*(-1 + b)*
	d[3]) + 2*(-737 + 398*b)*d[4]) + c[2]*(13*(11*(105 - 90*c[3] + b*(-42 + 
	90*c[3] - 52*c[4]) + 114*c[4])*d[0] + 2*(11*(-21 + 21*b - 90*c[3] + 64*b*
	c[3] + 10*(-1 + b)*c[4])*d[1] + (-1617 + 990*b + 154*c[3] - 154*b*c[3] - 
	1474*c[4] + 796*b*c[4])*d[2] + (-495 + 495*b + 1298*c[3] - 656*b*c[3] + 
	126*(-1 + b)*c[4])*d[3])) + 2*(8151 - 3718*b + 1638*(-1 + b)*c[3] - 94*
	(-169 + 74*b)*c[4])*d[4]) + 13*c[1]*(33*(-35 - 14*c[2] + 42*c[3] + b*(35 
	+ 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] - 462*(5*d[1] + d[2]-3*d[3])
	 + 858*d[4] + 2*(-22*(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] 
	+ 45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) - 2*(55*c[2] + 
	539*c[3] + 547*c[4])*d[4] + b*(22*(42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 
	11*(21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(-198 + 352*c[2] + 517*c[3]
	 + 420*c[4])*d[3] + (-429 + 110*c[2] + 840*c[3] + 1094*c[4])
	*d[4])))))))/ 2882880.;
		bf_mom[4] = ((a - b)*(429*c[0]*c[0]*(35*(-3 + a*a + a*b + b*b)*
	 d[0] - 7*(-10*d[2] + 6*a*b*d[2] + a*a*(5*d[1] + 2*d[2] - 3*d[3]) + b*b*
	(-5*d[1] + 2*d[2] + 3*d[3])) - 2*(-7 + 5*a*a - 3*a*b + 5*b*b)*d[4]) - 26*
	c[0]*(33*(7*(-5*(3 + 2*c[2]) + a*b*(5 + 6*c[2]) + a*a*(5 + 5*c[1] + 2*c[2] 
	- 3*c[3]) + b*b*(5 - 5*c[1] + 2*c[2] + 3*c[3])) + 2*(-7 + 5*a*a - 3*a*b + 
	5*b*b)*c[4])*d[0] + 11*(-(b*b*(3* (-35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*
	c[4])*d[1] + 2*(21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])*d[2] + (63 - 
	72*c[1] + 90*c[2] + 200*c[3] + 98*c[4])*d[3])) + 6*((35 + 98*c[2] - 38*c[4])
	*d[2] + 14*c[1]*(5*d[1] - 3*d[3]) + 6*c[3]*(-7*d[1] + 17*d[3]))) + 2*(231 
	- 1254*c[2] + b*b*(-165 + 429*c[1] + 286*c[2] - 539*c[3] - 1126*c[4]) + 
	3410*c[4])*d[4] + 2*a*b*(-11*((63 + 114*c[2] - 62*c[4])*d[2] + 6*c[1]*(7*
	d[1] - 9*d[3]) + c[3]*(-54*d[1] + 106*d[3])) + (99 + 682*c[2] - 1158*c[4])
	*d[4]) + a*a*(-33*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])* d[1] + 22*
	(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] + 11*(63 + 72*c[1] + 
	90*c[2] - 200*c[3] + 98*c[4])* d[3] - 2*(165 + 429*c[1] - 286*c[2] - 539*
	c[3] + 1126*c[4])* d[4])) + 2*(26*(-11* (3*(35*c[1]*c[1] + 7*c[2]*(5 + 7*
	c[2]) - 42*c[1]*c[3] + 51*c[3]*c[3]) - 3*(-7 + 38*c[2])*c[4] + 155*c[4]*
	c[4])*d[0] - 22*(3*c[1]*(-35 + 14*c[2] - 26*c[4]) + c[3]*(63 + 90*c[2] + 
	98*c[4]))*d[1] - 2*(231*c[1]*c[1] - 33*c[2]*(49 + 27*c[2]) + 990*c[1]*c[3] 
	- 649*c[3]*c[3] + 11*(57 + 134*c[2])*c[4] - 611*c[4]*c[4])*d[2] - 2*(11*
	c[1]*(63 + 90*c[2] + 98*c[4]) - c[3]*(1683 + 1298*c[2] + 826*c[4]))*d[3]) 
	+ 4*(13*(429*c[1]*c[1] - 11*c[2]*(57 + 67*c[2]) - 1078*c[1]*c[3] + 413*c[3]
	*c[3]) + 13*(1705 + 1222*c[2])*c[4] + 2331*c[4]*c[4])*d[4] + 2*a*b*(7579*
	c[3]*c[3]*d[0] - 1287*c[4]*d[0] + 7527*c[4]*c[4]*d[0] + 7722*c[3]*d[1] + 
	6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] + 8866*c[4]*d[2] - 8930*c[4]*c[4]
	*d[2] - 15158*c[3]*d[3] - 6420*c[3]*c[4]*d[3] - 26*c[1]*(11*(21 + 6*c[2] 
	+ 22*c[4])*d[1] - (297 + 286*c[2] + 238*c[4])*d[3] + c[3]*(297*d[0] - 286*
	d[2] - 238*d[4])) + 143*c[1]*c[1]*(21*d[0] - 6*d[2] - 22*d[4]) - 2*(1605*
	c[3]*c[3] + (7527 - 967*c[4])*c[4])*d[4] + 13*c[2]*c[2]*(627*d[0] - 1034*
	d[2] + 678*d[4]) + c[2]*(-143*(-63 + 62*c[4])*d[0] + 78*(-209 + 226*c[4])*
	d[2] + 52*c[3]*(143*d[1] - 321*d[3]) + 2*(4433 - 8930*c[4])*d[4])) + a*a*
	(-9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] + 4290*c[4]*d[0] - 14014*c[3]*c[4]*
	d[0] + 14638*c[4]*c[4]*d[0] + 10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*
	d[2] - 8528*c[3]*c[3]*d[2] + 7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*
	c[4]*c[4]*d[2] - 28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] 
	- 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(-21 + 45*c[3] 
	+ 26*c[4])*d[0] - 11*(21 + 64*c[3] - 10*c[4])*d[1] - 2*(-495 + 77*c[3] + 
	398*c[4])*d[2] + (-495 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*
	d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 
	3478*c[4])*d[4] - 2*(3764*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + 2*c[4]*
	(7319 + 1649*c[4]))*d[4] + 26*c[2]*c[2]*(495*d[0] - 363*d[1] - 374*d[2] + 
	77*d[3] + 398*d[4]) - 13*c[1]* (33*(-35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0]
	 + 2*(-22*(-42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 
	64*c[3] + 10*c[4])*d[2] - 2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + 
	(429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4]))) + b*b*(9009*c[3]*d[0] + 
	14300*c[3]*c[3]*d[0] + 4290*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 14638*c[4]*
	c[4]*d[0] + 10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 
	21840*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*
	c[3]*d[2] + 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 
	28600*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*
	d[3] - 10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(21 + 45*c[3] - 26*c[4])*d[0] + 
	11*(-21 + 64*c[3] + 10*c[4])*d[1] - 2*(495 + 77*c[3] - 398*c[4])*d[2] - 
	(495 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] + 63*d[1] + 
	24*d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(1859 + 819*c[3] - 3478*c[4])*d[4] - 
	2*(3764*c[3]*c[3] + 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*
	d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 77*d[3] + 398*d[4]) 
	+ 13*c[1]* (33*(-35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] + 2*(22*(-42 + 24*
	c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2]
	 + 2*(198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] + 840*
	c[3] + 1094*c[4])*d[4]))))))/720720.;
		bf_mom[5] = ((a - b)*(2*a*a*(13*(33* (7*(-20 + 5*c[0]*c[0] - 10*
	c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*(7*
	c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 
	26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4])*d[0] + 1716*(35*d[1] + 14*d[2] - 
	21*d[3] + 10*d[4]) - 429*c[0]*c[0]*(35*d[1] + 14*d[2] - 21*d[3] + 10*d[4]) 
	+ 52*c[0]*(-11*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + c[4]*(-39*d[1] + 26*
	d[2] + 49*d[3])) + (-539*c[3] + 1126*c[4])*d[4] + 33*c[1]*(28*d[1] - 7*d[2]
	 - 12*d[3] + 13*d[4]) - 11*c[2]*(21*d[1] - 90*d[2] + 45*d[3] + 26*d[4])) + 
	4*(-6721*c[3]*c[3]*d[1] + 10920*c[3]*c[4]*d[1] - 7111*c[4]*c[4]*d[1] - 
	4264*c[3]*c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 3478*c[4]*c[4]*d[2] + 26*c[2]*
	(c[3]*(352*d[1] + 77*d[2] - 328*d[3]) + c[4]*(-55*d[1] + 398*d[2] -63*d[3]))
	 + 6903*c[3]*c[3]*d[3] - 7528*c[3]*c[4]*d[3] + 5137*c[4]*c[4]*d[3] - 13*c[2]
	*c[2]*(363*d[1] + 374*d[2] - 77*d[3] - 398*d[4]) - 2*c[2]*(819*c[3] + 3478*
	c[4])*d[4] - 2*(1882*c[3]*c[3] - 5137*c[3]*c[4] + 1649*c[4]*c[4])* d[4] - 
	143*c[1]*c[1]* (63*d[1] - 24*d[2] - 9*d[3] + 28*d[4]) + 26*c[1]*(11*c[2]*
	(24*d[1] - 33*d[2] + 32*d[3] - 5*d[4]) + c[3]*(99*d[1] + 352*d[2] - 517*d[3]
	 + 420*d[4]) - c[4]*(308*d[1] + 55*d[2] - 420*d[3] + 547*d[4])))) + b*(13*
	(33*(7*(5*(3 + 2*b)*(-4 + c[0]*c[0]) + 20*(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*
	c[1]*c[1]) - 28*((5 + 2*b)*c[0] - 2*(1 + b)*c[1])*c[2] + 4*(49 + 30*b)*c[2]*
	c[2]) - 396*(7*(1 + b)*c[0] + 2*(7 + 4*b)*c[1] - 10*(1 + b)*c[2])* c[3] + 
	44*(153 + 100*b)*c[3]*c[3] - 44*(3*(7 + 10*b)*c[0] + 2*(39*(1 + b)*c[1] + 
	(57 + 26*b)*c[2] - 49*(1 + b)*c[3]))*c[4] + 4*(1705 + 1126*b)*c[4]*c[4])*
	d[0] - 24024*(5*d[1] - 5*d[2] - 3*d[3] - d[4]) + 2*(3003*c[0]*c[0]*(5*d[1] 
	- 5*d[2] - 3*d[3] - d[4]) + 572*c[0]*(-63*c[3]*d[1] - 39*c[4]*d[1] + 45*c[3]
	*d[2] - 57*c[4]*d[2] + 153*c[3]*d[3] + 49*c[4]*d[3] + 3*c[2]*(7*d[1] + 49*
	d[2] + 15*d[3] - 19*d[4]) + 3*c[1]*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) + 
	49*c[3]*d[4] + 155*c[4]*d[4]) + b*(13*(33*(-140 + 35*c[0]*c[0] + 84*c[1]*
	c[1] + 64*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2])) - 88*(18*c[0]
	 + 9*c[1] - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 616*c[1] - 110*
	c[2] - 840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 + 21*c[0]*c[0] 
	- 48*c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*(7*c[1] + 30*c[2])) 
	- 286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*(143*c[0] - 
	55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (13*(99*(-28 + 
	(c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2] + 308*c[2]
	*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3]) - 4*(91*
	(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])* c[4] + 20548*c[4]*c[4])*d[3] - 
	2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 539*c[3] - 1126*c[4]) + 
	4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*c[2]*c[3] + 1882*c[3]*c[3]
	 + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*c[4] - 13*c[1]*(55*c[2] + 
	420*c[3] + 547*c[4])))*d[4]) + 4*(6721*c[3]*c[3]*d[1] + 14014*c[3]*c[4]*d[1]
	 + 7111*c[4]*c[4]*d[1] - 8437*c[3]*c[3]*d[2] + 1638*c[3]*c[4]*d[2] - 7943*
	c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] - 10738*c[3]*c[4]*d[3] - 5137*c[4]*
	c[4]*d[3] + 429*c[1]*c[1]*(21*d[1] + 7*d[2] - 3*d[3] - 13*d[4]) - (5369*
	c[3]*c[3] + 10274*c[3]*c[4] + 2331*c[4]*c[4])*d[4] + 143*c[2]*c[2]*(33*d[1]
	 - 81*d[2] - 7*d[3] + 67*d[4]) + 26*c[2]*(c[4]*(55*d[1] + 737*d[2] + 63*d[3]
	 - 611*d[4]) + c[3]*(495*d[1] - 77*d[2] - 649*d[3] + 63*d[4])) + 26*c[1]*
	(11*c[2]*(21*d[1] + 33*d[2] + 45*d[3] + 5*d[4]) + 11*c[3]*(-9*d[1] + 45*d[2]
	 + 47*d[3] + 49*d[4]) + c[4]*(-429*d[1] + 55*d[2] + 539*d[3] + 547*d[4])))))
	 + a*(13*(33*(7*(5*(3 + 2*b)*(-4 + c[0]*c[0]) - 20*c[0]*c[1] + 4*(5 + 2*b)*
	c[1]*c[1]) - 28*((5 + 6*b)*c[0] + 2*c[1])*c[2] + 4*(49 + 38*b)*c[2]*c[2]) + 
	396*(7*c[0] - 2*((7 + 6*b)*c[1] + 5*c[2]))*c[3] + 44*(153 + 106*b)*c[3]*c[3]
	 + 44*(3*(-7 + 6*b)*c[0] - 2*(-39*c[1] + (57 + 62*b)*c[2] + 49*c[3]))*c[4] 
	+ 4*(1705 + 1158*b)*c[4]*c[4])*d[0] + 24024*(5*d[1] + 5*d[2] - 3*d[3]+d[4])
	 + 2*(-429*c[0]*c[0]* (35*d[1] + 7*(5 + 6*b)*d[2] - 21*d[3] + (7-6*b)*d[4])
	 + 52*c[0]*(-11*(3*c[2]*(7*d[1] - (49 + 38*b)*d[2] + 15*d[3]) + c[4]*(-39*
	d[1] + (57 + 62*b)*d[2] + 49*d[3]) + c[3]*(9*(7 + 6*b)*d[1] + 45*d[2] - 
	(153 + 106*b)*d[3])) - (11*(57 + 62*b)*c[2] + 539*c[3] - (1705 + 1158*b)*
	c[4])*d[4] + 33*c[1]*(7*(5 + 2*b)*d[1] - 7*d[2] - 3*(7 + 6*b)*d[3] + 13*
	d[4])) + 4* (-6721*c[3]*c[3]*d[1] + 14014*c[3]*c[4]*d[1] + 6188*b*c[3]*c[4]*
	d[1] - 7111*c[4]*c[4]*d[1] + 18018*b*d[2] - 8437*c[3]*c[3]*d[2] - 8346*b*
	c[3]*c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 7943*c[4]*c[4]*d[2] - 8930*b*c[4]*
	c[4]*d[2] + 6903*c[3]*c[3]*d[3] - 10738*c[3]*c[4]*d[3] - 6420*b*c[3]*c[4]*
	d[3] + 5137*c[4]*c[4]*d[3] + 26*c[2]*(c[4]*(-55*d[1] + (737 + 678*b)*d[2] 
	- 63*d[3]) + c[3]*(11*(45 + 26*b)*d[1] + 77*d[2] - (649 + 642*b)*d[3])) - 
	2*c[2]*(819*c[3] + 47*(169 + 190*b)*c[4])*d[4] - (5369*c[3]*c[3] + 6*b*(429 
	+ 535*c[3]*c[3]) - 10274*c[3]*c[4] - (-2331 + 1934*b)*c[4]*c[4])*d[4] - 
	143*c[1]*c[1]* (63*d[1] + 3*(-7 + 2*b)*d[2] - 9*d[3] + (39 + 22*b)*d[4]) - 
	13*c[2]*c[2]* (363*d[1] + 11*(81 + 94*b)*d[2] - 77*d[3] - (737 + 678*b)*
	d[4]) - 26*c[1]*(11*c[2]* (3*(-7 + 2*b)*d[1] + 33*d[2] - (45 + 26*b)*d[3] 
	+ 5*d[4]) + c[4]*(11*(39 + 22*b)*d[1] + 55*d[2] - 7*(77 + 34*b)*d[3] + 547*
	d[4]) + c[3]*(-99*d[1] - 11*(45 + 26*b)*d[2] + 517*d[3] - 
	7*(77 + 34*b)*d[4])))))))/1441440.;
		bf_mom[6] = ((a - b)*(429*c[0]*c[0]*(35*(-3 + a*a + a*b + b*b)* 
	d[0] - 7*(-10*d[2] + 6*a*b*d[2] + a*a*(5*d[1] + 2*d[2] - 3*d[3]) + b*b*(-5*
	d[1] + 2*d[2] + 3*d[3])) - 2*(-7 + 5*a*a - 3*a*b + 5*b*b)*d[4]) - 26*c[0]*
	(33*(7*(15 - 10*c[2] + a*b*(-5 + 6*c[2]) + a*a*(-5 + 5*c[1] + 2*c[2] - 3*
	c[3]) + b*b*(-5 - 5*c[1] + 2*c[2] + 3*c[3])) + 2*(-7 + 5*a*a - 3*a*b+5*b*b)
	*c[4])*d[0] + 11*(b*b*(-3*(35 + 56*c[1] + 14*c[2] - 24*c[3] - 26*c[4])*d[1]
	 - 2*(-21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])*d[2] + (63 + 72*c[1] - 
	90*c[2] - 200*c[3] - 98*c[4])*d[3]) + 6*((-35 + 98*c[2] - 38*c[4])*d[2] + 
	14*c[1]*(5*d[1] - 3*d[3]) + 6*c[3]*(-7*d[1] + 17*d[3]))) + 2*(-33*(7 + 38*
	c[2]) + b*b*(165 + 429*c[1] + 286*c[2] - 539*c[3] - 1126*c[4]) + 3410*c[4])
	*d[4] + 2*a*b*(11*((63 - 114*c[2] + 62*c[4])*d[2] + 2*c[3]*(27*d[1] - 53*
	d[3])) - 66*c[1]*(7*d[1] - 9*d[3]) + (-99 + 682*c[2] - 1158*c[4])*d[4]) + 
	a*a*(11*((105 - 168*c[1] + 42*c[2] + 72*c[3] - 78*c[4])* d[1] + 2*(21 + 21*
	c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] + (-63 + 72*c[1] + 90*c[2] - 200*
	c[3] + 98*c[4])*d[3]) - 2*(429*c[1] - 11*(15 + 26*c[2] + 49*c[3]) + 1126*
	c[4])*d[4])) + 2*(26*(-1683*c[3]*c[3]*d[0] + 231*c[4]*d[0] - 1705*c[4]*c[4]*
	d[0] + 1386*c[3]*d[1] - 2156*c[3]*c[4]*d[1] - 33*c[2]*c[2]*(49*d[0] - 54*
	d[2]) + 1298*c[3]*c[3]*d[2] + 1254*c[4]*d[2] + 1222*c[4]*c[4]*d[2] - 231*
	c[1]*c[1]*(5*d[0] + 2*d[2]) + 2*c[3]*(-1683 + 826*c[4])*d[3] + 22*c[1]*(-3*
	(35 + 14*c[2] - 26*c[4])*d[1] + 9*c[3]*(7*d[0] - 10*d[2]) + (63 - 90*c[2] 
	- 98*c[4])*d[3]) + 11*c[2]*(3*(35 + 38*c[4])*d[0] - 294*d[2] - 268*c[4]*
	d[2] + 4*c[3]*(-45*d[1] + 59*d[3]))) + 4*(13*(429*c[1]*c[1] + 11*(57 - 67*
	c[2])*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3]) + 13*(-1705 + 1222*c[2])*c[4] 
	+ 2331*c[4]*c[4])*d[4] + 2*a*b*(7579*c[3]*c[3]*d[0] + 1287*c[4]*d[0] + 
	7527*c[4]*c[4]*d[0] - 7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*
	d[2] - 8866*c[4]*d[2] - 8930*c[4]*c[4]*d[2] + 15158*c[3]*d[3] - 6420*c[3]*
	c[4]*d[3] - 26*c[1]*(11*(-21 + 6*c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 
	238*c[4])*d[3] + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*
	(21*d[0] - 6*d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] - c[4]*(7527 + 967*c[4]))*
	d[4] + 13*c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(63 + 
	62*c[4])*d[0] + 78*(209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) - 
	2*(4433 + 8930*c[4])*d[4])) + a*a*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 
	4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*
	d[1] - 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 
	14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*
	d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 28600*c[3]*d[3] + 13806*
	c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*
	d[3] - 26*c[2]*(11*(21 + 45*c[3] + 26*c[4])*d[0] + (231 - 704*c[3] + 110*
	c[4])*d[1] - 2*(495 + 77*c[3] + 398*c[4])*d[2] + (495 + 656*c[3] + 126*c[4])
	*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 
	4*c[2]*(1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(c[3]*(7007 + 3764*c[3]) - 2*
	(7319 + 5137*c[3])*c[4] + 3298*c[4]*c[4])*d[4] + 26*c[2]*c[2]*(495*d[0] - 
	363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]* (33*(35 + 14*c[2] + 
	24*c[3] - 26*c[4])*d[0] + 2*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 
	11*(21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3] 
	+ 420*c[4])*d[3] + (-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4]))) + b*b*
	(-9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] + 14014*c[3]*c[4]*
	d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 12870*c[3]*
	d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*
	c[4]*c[4]*d[2] + 28600*c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 
	15056*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(-21 + 45*c[3] - 
	26*c[4])*d[0] + 11*(21 + 64*c[3] + 10*c[4])*d[1] - 2*(-495 + 77*c[3] - 398*
	c[4])*d[2] - (-495 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] 
	+ 63*d[1] + 24*d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(-1859 + 819*c[3] - 3478*
	c[4])*d[4] - 2*(3764*c[3]*c[3] + 11*c[3]*(-637 + 934*c[4]) + 2*c[4]*(-7319 
	+ 1649*c[4]))*d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 77*d[3] 
	+ 398*d[4]) + 13*c[1]* (33*(35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] + 2*(22*
	(42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3] + 10*
	c[4])*d[2] + 2*(-198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 + 
	110*c[2] + 840*c[3] + 1094*c[4])*d[4]))))))/720720.;
		bf_mom[7] = ((a - b)*(2*a*a*(13*(33* (7*(-20 + 5*c[0]*c[0] - 10*c[0]
	*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*(7*c[0] 
	- 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 26*c[2]
	 + 49*c[3])*c[4] + 2252*c[4]*c[4])*d[0] + 1716*(35*d[1] + 14*d[2] - 21*d[3]
	 + 10*d[4]) - 429*c[0]*c[0]*(35*d[1] + 14*d[2] - 21*d[3] + 10*d[4]) + 52*
	c[0]*(-11*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + c[4]*(-39*d[1] + 26*d[2] 
	+ 49*d[3])) + (-539*c[3] + 1126*c[4])*d[4] + 33*c[1]*(28*d[1] - 7*d[2] - 
	12*d[3] + 13*d[4]) - 11*c[2]*(21*d[1] - 90*d[2] + 45*d[3] + 26*d[4])) + 4*
	(-6721*c[3]*c[3]*d[1] + 10920*c[3]*c[4]*d[1] - 7111*c[4]*c[4]*d[1] - 4264*
	c[3]*c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 3478*c[4]*c[4]*d[2] + 26*c[2]*(c[3]*
	(352*d[1] + 77*d[2] - 328*d[3]) + c[4]*(-55*d[1] + 398*d[2] - 63*d[3])) + 
	6903*c[3]*c[3]*d[3] - 7528*c[3]*c[4]*d[3] + 5137*c[4]*c[4]*d[3] - 13*c[2]*
	c[2]*(363*d[1] + 374*d[2] - 77*d[3] - 398*d[4]) - 2*c[2]*(819*c[3] + 3478*
	c[4])*d[4] - 2*(1882*c[3]*c[3] - 5137*c[3]*c[4] + 1649*c[4]*c[4])* d[4] - 
	143*c[1]*c[1]* (63*d[1] - 24*d[2] - 9*d[3] + 28*d[4]) + 26*c[1]*(11*c[2]*
	(24*d[1] - 33*d[2] + 32*d[3] - 5*d[4]) + c[3]*(99*d[1] + 352*d[2] - 517*d[3]
	 + 420*d[4]) - c[4]*(308*d[1] + 55*d[2] - 420*d[3] + 547*d[4])))) + b*(13*
	(33*(7*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)
	*c[1]*c[1]) - 28*((-5 + 2*b)*c[0] - 2*(-1 + b)*c[1])*c[2] + 4*(-49 + 30*b)*
	c[2]*c[2]) - 396*(7*(-1 + b)*c[0] + 2*(-7 + 4*b)*c[1] - 10*(-1 + b)*c[2])*
	c[3] + 44*(-153 + 100*b)*c[3]*c[3] - 44*(3*(-7 + 10*b)*c[0] + 78*(-1 + b)*
	c[1] + 2*(-57 + 26*b)*c[2] - 98*(-1 + b)*c[3])*c[4] + 4*(-1705 + 1126*b)*
	c[4]*c[4])*d[0] + 24024*(5*d[1] - 5*d[2] - 3*d[3] - d[4]) + 2*(-3003*c[0]*
	c[0]*(5*d[1] - 5*d[2] - 3*d[3] - d[4]) - 572*c[0]*(-63*c[3]*d[1] - 39*c[4]*
	d[1] + 45*c[3]*d[2] - 57*c[4]*d[2] + 153*c[3]*d[3] + 49*c[4]*d[3] + 3*c[2]*
	(7*d[1] + 49*d[2] + 15*d[3] - 19*d[4]) + 3*c[1]*(7*(5*d[1] + d[2] - 3*d[3])
	 - 13*d[4]) + 49*c[3]*d[4] + 155*c[4]*d[4]) + b*(13*(33*(-140 + 35*c[0]*c[0]
	 + 84*c[1]*c[1] + 64*c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2])) - 
	88*(18*c[0] + 9*c[1] - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 616*
	c[1] - 110*c[2] - 840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 + 21*
	c[0]*c[0] - 48*c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*(7*c[1] + 
	30*c[2])) - 286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*
	(143*c[0] - 55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - 
	(13*(99*(-28 + (c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])
	*c[2] + 308*c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]
	*c[3]) - 4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])* c[4] + 20548*
	c[4]*c[4])*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 539*
	c[3] - 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*c[2]*
	c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*c[4] - 
	13*c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4]) - 4*(6721*c[3]*c[3]*d[1] + 
	14014*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] - 8437*c[3]*c[3]*d[2] + 1638*c[3]
	*c[4]*d[2] - 7943*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3]-10738*c[3]*c[4]*d[3] 
	- 5137*c[4]*c[4]*d[3] + 429*c[1]*c[1]*(21*d[1] + 7*d[2] - 3*d[3] - 13*d[4])
	 - (5369*c[3]*c[3] + 10274*c[3]*c[4] + 2331*c[4]*c[4])*d[4] + 143*c[2]*c[2]*
	(33*d[1] - 81*d[2] - 7*d[3] + 67*d[4]) + 26*c[2]*(c[4]*(55*d[1] + 737*d[2] 
	+ 63*d[3] - 611*d[4]) + c[3]*(495*d[1] - 77*d[2] - 649*d[3] + 63*d[4])) + 
	26*c[1]*(11*c[2]*(21*d[1] + 33*d[2] + 45*d[3] + 5*d[4]) + 11*c[3]*(-9*d[1] 
	+ 45*d[2] + 47*d[3] + 49*d[4]) + c[4]*(-429*d[1] + 55*d[2] + 539*d[3] + 
	547*d[4]))))) + a*(13*(33*(7*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 20*c[0]*c[1] 
	+ 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5 + 6*b)*c[0] - 2*c[1])*c[2] + 4*(-49 + 
	38*b)*c[2]*c[2]) - 396*(7*c[0] + 2*(-7 + 6*b)*c[1] - 10*c[2])*c[3] + 44*
	(-153 + 106*b)*c[3]*c[3] + 44*(3*(7 + 6*b)*c[0] - 78*c[1] + 2*(57 - 62*b)*
	c[2] + 98*c[3])*c[4] + 4*(-1705 + 1158*b)*c[4]*c[4])*d[0] - 24024*(5*d[1] 
	+ 5*d[2] - 3*d[3] + d[4]) + 2*(429*c[0]*c[0]* (7*(5*d[1] + (5 - 6*b)*d[2] 
	- 3*d[3]) + (7 + 6*b)*d[4]) + 52*c[0]*(11*(63*c[3]*d[1] - 54*b*c[3]*d[1] - 
	39*c[4]*d[1] + 45*c[3]*d[2] + 57*c[4]*d[2] - 62*b*c[4]*d[2] + ((-153+106*b)
	*c[3] + 49*c[4])*d[3] + 3*c[2]*(7*d[1] + (-49 + 38*b)*d[2] + 15*d[3]) + 3*
	c[1]*(7*(-5 + 2*b)*d[1] + 7*d[2] + 3*(7 - 6*b)*d[3])) - (429*c[1] + 11*(-57
	 + 62*b)*c[2] - 539*c[3] + (1705 - 1158*b)*c[4])*d[4]) + 4*(6721*c[3]*c[3]*
	d[1] - 14014*c[3]*c[4]*d[1] + 6188*b*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] + 
	18018*b*d[2] + 8437*c[3]*c[3]*d[2] - 8346*b*c[3]*c[3]*d[2] + 1638*c[3]*c[4]
	*d[2] + 7943*c[4]*c[4]*d[2] - 8930*b*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] 
	+ 10738*c[3]*c[4]*d[3] - 6420*b*c[3]*c[4]*d[3] - 5137*c[4]*c[4]*d[3] - 
	(-5369*c[3]*c[3] + 6*b*(429 + 535*c[3]*c[3]) + 10274*c[3]*c[4] - (2331 + 
	1934*b)*c[4]*c[4])*d[4] + 143*c[1]*c[1]* (63*d[1] - 3*(7 + 2*b)*d[2] - 9*
	d[3] + (39 - 22*b)*d[4]) + 13*c[2]*c[2]* (11*(33*d[1] + (81 - 94*b)*d[2] - 
	7*d[3]) + (-737 + 678*b)*d[4]) + 2*c[2]*(13*c[3]* (11*(-45 + 26*b)*d[1] - 
	77*d[2] + (649 - 642*b)*d[3] + 63*d[4]) + c[4]*(715*d[1] + 13*(-737 + 678*b)
	*d[2] + 819*d[3] + 47*(169 - 190*b)*d[4])) + 26*c[1]*(-11*c[2]* (3*(7 + 2*b)
	*d[1] - 33*d[2] + (45 - 26*b)*d[3] - 5*d[4]) + c[4]*((429 - 242*b)*d[1] + 
	55*d[2] + 7*(-77 + 34*b)*d[3] + 547*d[4]) + c[3]*(-99*d[1] + 11*(-45 + 26*b)
	*d[2] + 517*d[3] + 7*(-77 + 34*b)*d[4])))))))/1441440.;
		bf_mom[8] = -((a - b)*(13*(33*(35*(-3 + b*b)*(-4 + c[0]*c[0]) + 70*
	b*b*c[0]*c[1] + 28*(-5 + 2*b*b)*c[1]*c[1] - 28*((-5 + b*b)*c[0] - b*b*c[1])
	*c[2] + 4*(-49 + 15*b*b)*c[2]*c[2]) - 198*(-28*c[1] + b*b*(7*c[0] + 8*c[1] 
	- 10*c[2]))*c[3] + 44*(-153 + 50*b*b)*c[3]*c[3] - 44*(3*(-7 + 5*b*b)*c[0] - 
	114*c[2] + b*b*(39*c[1] + 26*c[2] - 49*c[3]))*c[4] + 4*(-1705 + 563*b*b)*
	c[4]*c[4] + a*a*(33*(7*(-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]*c[1]) - 
	28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*(7*c[0] - 8*c[1] - 10*c[2])*
	c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 26*c[2] + 49*c[3])*c[4] + 
	2252*c[4]*c[4]) + a*b*(33*(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*c[2]) + 
	44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3]) + 44*(9*c[0] - 62*
	c[2])*c[4] + 2316*c[4]*c[4]))*d[0] - 60060*b*b*d[1] + 15015*b*b*c[0]*c[0]*
	d[1] - 120120*c[0]*c[1]*d[1] + 48048*b*b*c[0]*c[1]*d[1] + 36036*b*b*c[1]*
	c[1]*d[1] + 12012*b*b*c[0]*c[2]*d[1] - 48048*c[1]*c[2]*d[1] + 27456*b*b*
	c[1]*c[2]*d[1] + 18876*b*b*c[2]*c[2]*d[1] + 72072*c[0]*c[3]*d[1] - 20592*b*
	b*c[0]*c[3]*d[1] - 10296*b*b*c[1]*c[3]*d[1] - 102960*c[2]*c[3]*d[1] + 36608*
	b*b*c[2]*c[3]*d[1] + 26884*b*b*c[3]*c[3]*d[1] - 22308*b*b*c[0]*c[4]*d[1] + 
	89232*c[1]*c[4]*d[1] - 32032*b*b*c[1]*c[4]*d[1] + 5720*b*b*c[2]*c[4]*d[1] - 
	112112*c[3]*c[4]*d[1] + 43680*b*b*c[3]*c[4]*d[1] + 28444*b*b*c[4]*c[4]*d[1]
	 - 120120*d[2] + 24024*b*b*d[2] + 30030*c[0]*c[0]*d[2] - 6006*b*b*c[0]*c[0]*
	d[2] + 12012*b*b*c[0]*c[1]*d[2] - 24024*c[1]*c[1]*d[2] + 13728*b*b*c[1]*
	c[1]*d[2] - 168168*c[0]*c[2]*d[2] + 51480*b*b*c[0]*c[2]*d[2] + 37752*b*b*
	c[1]*c[2]*d[2] + 92664*c[2]*c[2]*d[2] - 19448*b*b*c[2]*c[2]*d[2] + 25740*b*
	b*c[0]*c[3]*d[2] - 102960*c[1]*c[3]*d[2] + 36608*b*b*c[1]*c[3]*d[2] - 8008*
	b*b*c[2]*c[3]*d[2] + 67496*c[3]*c[3]*d[2] - 17056*b*b*c[3]*c[3]*d[2] + 
	65208*c[0]*c[4]*d[2] - 14872*b*b*c[0]*c[4]*d[2] + 5720*b*b*c[1]*c[4]*d[2] 
	- 153296*c[2]*c[4]*d[2] + 41392*b*b*c[2]*c[4]*d[2] + 6552*b*b*c[3]*c[4]*d[2]
	 + 63544*c[4]*c[4]*d[2] - 13912*b*b*c[4]*c[4]*d[2] + 36036*b*b*d[3] - 9009*
	b*b*c[0]*c[0]*d[3] + 72072*c[0]*c[1]*d[3] - 20592*b*b*c[0]*c[1]*d[3] - 5148*
	b*b*c[1]*c[1]*d[3] + 25740*b*b*c[0]*c[2]*d[3] - 102960*c[1]*c[2]*d[3] + 
	36608*b*b*c[1]*c[2]*d[3] - 4004*b*b*c[2]*c[2]*d[3] - 175032*c[0]*c[3]*d[3] 
	+ 57200*b*b*c[0]*c[3]*d[3] + 53768*b*b*c[1]*c[3]*d[3] + 134992*c[2]*c[3]*
	d[3] - 34112*b*b*c[2]*c[3]*d[3] - 27612*b*b*c[3]*c[3]*d[3] + 28028*b*b*c[0]*
	c[4]*d[3] - 112112*c[1]*c[4]*d[3] + 43680*b*b*c[1]*c[4]*d[3] + 6552*b*b*
	c[2]*c[4]*d[3] + 85904*c[3]*c[4]*d[3] - 30112*b*b*c[3]*c[4]*d[3] - 20548*b*
	b*c[4]*c[4]*d[3] - 2*(-3003*(-4 + c[0]*c[0]) - 52*(429*c[1]*c[1] + 627*c[0]*
	c[2] - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3]) + 52*(1705*c[0] - 
	1222*c[2])*c[4] - 9324*c[4]*c[4] + b*b*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] 
	+ 286*c[2] - 539*c[3] - 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*
	c[2] - 819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 
	1649*c[4]*c[4] - 13*c[1]*(55*c[2] + 420*c[3] + 547*c[4]))))*d[4] - a*a*(13*
	(1155*c[0]*c[0] + 44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*c[3] + 
	47*c[3]*c[3] - 6*c[1]*(8*c[2] + 3*c[3])) + 8*(308*c[1] + 55*c[2] - 420*c[3])
	*c[4] + 2188*c[4]*c[4] - 132*c[0]*(28*c[1] - 7*c[2] - 12*c[3] + 13*c[4]))*
	d[1] + 2*(143*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*(7*c[1] - 30*c[2])
	 + 132*c[1]*c[2] + 68*c[2]*c[2]) + 286*(45*c[0] - 64*c[1] - 14*c[2])*c[3] + 
	8528*c[3]*c[3] + 52*(143*c[0] + 55*c[1] - 398*c[2] + 63*c[3])*c[4] + 6956*
	c[4]*c[4])*d[2] - (9009*c[0]*c[0] + 52*(11*(-63 + (9*c[1] + c[2])*(c[1] + 
	7*c[2])) - 2*(517*c[1] + 328*c[2])*c[3] + 531*c[3]*c[3]) + 8*(5460*c[1] - 
	819*c[2] - 3764*c[3])*c[4] + 20548*c[4]*c[4] - 572*c[0]*(36*c[1] + 45*c[2] 
	- 100*c[3] + 49*c[4]))*d[3] + 2*(2145*c[0]*c[0] - 26*c[0]*(429*c[1] - 286*
	c[2] - 539*c[3] + 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] + 
	819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] - 5137*c[3]*c[4] + 1649*
	c[4]*c[4] + 13*c[1]*(55*c[2] - 420*c[3] + 547*c[4])))*d[4]) + 2*a*b*(-104*
	(11*c[1]*(3*c[2] + 11*c[4]) - c[3]*(143*c[2] + 119*c[4]))*d[1] - 4*(429*
	c[1]*c[1] + 6721*c[2]*c[2] - 3718*c[1]*c[3] + 4173*c[3]*c[3] - 8814*c[2]*
	c[4] + 4465*c[4]*c[4])*d[2] + 8*(-321*c[3]*(13*c[2] + 5*c[4]) + 13*c[1]*
	(143*c[2] + 119*c[4]))*d[3] + 5148*(7*d[2] - d[4]) - 1287*c[0]*c[0]*(7*d[2]
	 - d[4]) - 4*(1573*c[1]*c[1] - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*
	c[3] + 8930*c[2]*c[4] - 967*c[4]*c[4])*d[4] + 52*c[0]*(627*c[2]*d[2] - 341*
	c[4]*d[2] + 33*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-297*d[1] + 583*d[3]) - 341*
	c[2]*d[4] + 579*c[4]*d[4]))))/ 360360.;
		break;
	case 101:
	case 103:
		bf_mom[0] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] 
	- 2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(5 + 5*c[1] + 2*c[2] - 3*
	c[3]) + 10*c[4])* d[0] + 33*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])*
	 d[1] - 22*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(63 +
	 72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(165 + 429*c[1] - 286*
	c[2] - 539*c[3] + 1126*c[4])* d[4]) + 2*(-9009*c[3]*d[0] + 14300*c[3]*c[3]
	*d[0] + 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 
	10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 21840*c[3]*c[4]
	*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] + 
	7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*
	d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 
	10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] - 11*
	(21 + 64*c[3] - 10*c[4])*d[1] - 2*(-495 + 77*c[3] + 398*c[4])*d[2] + (-495
	 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1]+24*d[2] 
	+ 9*d[3] - 28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(3764
	*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*d[4] + 
	26*c[2]*c[2]* (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*
	c[1]*(33*(-35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 24*c[2] 
	+ 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 
	2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] - 840*c[3]
	 + 1094*c[4])*d[4])))) + a*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*d[1]
	 + (5 - 6*b)*d[2] - 3*d[3]) + 2*(7 + 6*b)*d[4]) - 26*c[0]*(33*(7*(-5*(3 + 
	2*c[1] + 2*c[2]) + 2*b*(5 + 6*c[2]) + 6*c[3]) - 2*(7 + 6*b)*c[4])*d[0] + 
	22*(-3*(-35 + 14*(-5 + 2*b)*c[1] + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*c[4])*
	d[1] - (2*b*(63 + 114*c[2] - 62*c[4]) + 3*(-35 + 14*c[1] - 98*c[2] + 30*
	c[3] + 38*c[4]))*d[2] + (-63 + 18*(-7 + 6*b)*c[1] - 90*c[2] + 306*c[3] - 
	212*b*c[3] - 98*c[4])* d[3]) + 2*(2*b*(99 + 682*c[2] - 1158*c[4]) + 11*
	(21 + 78*c[1] - 114*c[2] - 98*c[3] + 310*c[4]))*d[4]) + 4*(9009*c[3]*d[0] 
	- 21879*c[3]*c[3]*d[0] + 15158*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] - 2574*b*
	c[4]*d[0] + 14014*c[3]*c[4]*d[0] - 22165*c[4]*c[4]*d[0] + 15054*b*c[4]*
	c[4]*d[0] - 18018*c[3]*d[1] + 15444*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 
	11154*c[4]*d[1] - 28028*c[3]*c[4]*d[1] + 12376*b*c[3]*c[4]*d[1] + 14222*
	c[4]*c[4]*d[1] - 12870*c[3]*d[2] + 16874*c[3]*c[3]*d[2] - 16692*b*c[3]*
	c[3]*d[2] - 16302*c[4]*d[2] + 17732*b*c[4]*d[2] + 3276*c[3]*c[4]*d[2] + 
	15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] - 30316*b*
	c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] + 21476*c[3]*c[4]*d[3] 
	- 12840*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] - 2*((-5369 + 3210*b)*c[3]
	*c[3] + 11*c[3]*(637 + 934*c[4]) - c[4]*(22165 + 2331*c[4] + 2*b*(-7527 + 
	967*c[4])))*d[4] + 13*c[2]*c[2]* (33*(-49 + 38*b)*d[0] + 22*(33*d[1] + 
	(81 - 94*b)*d[2] - 7*d[3]) + 2*(-737 + 678*b)*d[4]) + c[2]*(-143*(-3*(-35 
	+ 42*b + 30*c[3]) + 2*(-57 + 62*b)*c[4])*d[0] + 52*(-627*b*d[2] + c[4]*
	(55*d[1] + (-737 + 678*b)*d[2] + 63*d[3]) + c[3]*(11*(-45 + 26*b)*d[1] - 
	77*d[2] + (649 - 642*b)*d[3])) + 4*(819*c[3] + b*(4433 - 8930*c[4]) + 
	7943*c[4])*d[4] - 858*(7*d[1] - 49*d[2] + 15*d[3] + 19*d[4])) + 143*c[1]*
	c[1]* (21*(-5 + 2*b)*d[0] - 2*(-63*d[1] + 3*(7 + 2*b)*d[2] + 9*d[3] + 
	(-39 + 22*b)*d[4])) + 13*c[1]*(33*(-35 + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*
	c[4])* d[0] + 2*(-11* (3*(-35 + 14*c[2] + 6*c[3] - 26*c[4]) + 2*b*(21 + 
	6*c[2] + 22*c[4]))*d[1] + 11*(-21 + 66*c[2] + (-90 + 52*b)*c[3] + 10*c[4])
	*d[2] + (-11*(63 + 90*c[2] - 94*c[3] + 98*c[4]) + b*(594 + 572*c[2] + 476*
	c[4]))*d[3] + (429 + 110*c[2] + 14*(-77 + 34*b)*c[3] + 1094*c[4])*d[4]))))
	 + b*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*(-1 + b)*d[1] + (5 - 2*b)*
	d[2] - 3*(-1 + b)*d[3]) - 2*(-7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(15 - 10*
	c[1] + 10*c[2] + 2*b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3]) - 2*(-7 + 
	10*b)*c[4])*d[0] + 22*(3*(35 - 70*c[1] - 14*c[2] + 42*c[3] + b*(-35 + 56*
	c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] + (-3*(35 + 14*c[1] +
	 98*c[2] + 30*c[3] - 38*c[4]) + 2*b*(21 + 21*c[1] + 90*c[2] + 45*c[3] - 
	26*c[4]))*d[2] + (-63 + 63*b + 126*c[1] - 72*b*c[1] - 90*c[2] + 90*b*c[2] 
	- 306*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])*d[3]) - 2*(11*(21 - 78*c[1] - 
	114*c[2] + b*(-30 + 78*c[1] + 52*c[2] - 98*c[3]) + 98*c[3]) - 2*(-1705 + 
	1126*b)*c[4])*d[4]) + 4*(-9009*c[3]*d[0] + 9009*b*c[3]*d[0] - 21879*c[3]*
	c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] + 4290*b*c[4]*d[0] - 
	14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 22165*c[4]*c[4]*d[0] + 
	14638*b*c[4]*c[4]*d[0] - 18018*c[3]*d[1] + 10296*b*c[3]*d[1] - 13442*c[3]*
	c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 11154*b*c[4]*d[1] 
	- 28028*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] + 
	14222*b*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 12870*b*c[3]*d[2] + 16874*c[3]*
	c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] - 16302*c[4]*d[2] + 7436*b*c[4]*d[2] - 
	3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*d[2] - 6956*
	b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] - 28600*b*c[3]*d[3] + 13806*c[3]*c[3]*
	d[3] - 13806*b*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 14014*b*c[4]*d[3] + 
	21476*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 
	10274*b*c[4]*c[4]*d[3] - 2*((-5369 + 3764*b)*c[3]*c[3] + 11*(-1 + b)*c[3]*
	(637 + 934*c[4]) + c[4]*(-22165 - 2331*c[4] + 2*b*(7319 + 1649*c[4])))*
	d[4] + 143*c[1]*c[1]* (3*(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*d[1] + (-7 + 
	8*b)*d[2] - 3*(-1 + b)*d[3])) - 2*(-39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*
	(-49 + 30*b)*d[0] + 22*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] - 7*(-1 + b)*
	d[3]) + 2*(-737 + 398*b)*d[4]) + c[2]*(143*(b*(42 + 90*c[3] - 52*c[4]) - 
	3*(35 + 30*c[3] - 38*c[4]))*d[0] + 286*(21 - 21*b - 90*c[3] + 64*b*c[3] + 
	10*(-1 + b)*c[4])* d[1] - 26*(-77*(21 + 2*c[3]) + 2*b*(495 + 77*c[3] - 
	398*c[4]) + 1474*c[4])*d[2] - 26*(-495 + 495*b - 1298*c[3] + 656*b*c[3] - 
	126*(-1 + b)*c[4])*d[3] + 2*(143*(-57 + 26*b) + 1638*(-1 + b)*c[3] - 94*
	(-169 + 74*b)*c[4])*d[4]) + 13*c[1]*(33*(35 - 14*c[2] + 42*c[3] + b*(-35 
	+ 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] + 462*(5*d[1] + d[2] - 3*
	d[3]) - 858*d[4] + 2*(-22*(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*
	(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) - 2*
	(55*c[2] + 539*c[3] + 547*c[4])*d[4] + b*(22*(-42 + 24*c[2] - 9*c[3] - 28*
	c[4])*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(198 + 352*
	c[2] + 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] + 840*c[3] + 
	1094*c[4])*d[4])))))))/ 2882880.;
		bf_mom[1] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] -
	 2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(-5 + 5*c[1] + 2*c[2] - 3*
	c[3]) + 10*c[4])* d[0] + 33*(-35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])*
	 d[1] - 22*(21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(-63 + 
	72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(429*c[1] - 11*(15 + 26*
	c[2] + 49*c[3]) + 1126*c[4])*d[4]) + 2*(9009*c[3]*d[0] + 14300*c[3]*c[3]*
	d[0] - 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 
	10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 21840*c[3]*c[4]
	*d[1] - 14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 
	7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 28600*c[3]*
	d[3] + 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 
	10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(21 + 45*c[3] + 26*c[4])*d[0] + (231 - 
	704*c[3] + 110*c[4])*d[1] - 2*(495 + 77*c[3] + 398*c[4])*d[2] + (495 + 
	656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] 
	+ 9*d[3] - 28*d[4]) - 4*c[2]*(1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(c[3]*
	(7007 + 3764*c[3]) - 2*(7319 + 5137*c[3])*c[4] + 3298*c[4]*c[4])*d[4] + 
	26*c[2]*c[2]* (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*
	c[1]*(33*(35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(42 + 24*c[2] + 
	9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*
	(-198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] - 840*c[3]
	 + 1094*c[4])*d[4])))) + a*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*d[1]
	 + (5 - 6*b)*d[2] - 3*d[3]) + 2*(7 + 6*b)*d[4]) - 26*c[0]*(11*(3*(7* (15 
	- 10*c[1] - 10*c[2] + 2*b*(-5 + 6*c[2]) + 6*c[3]) - 2*(7 + 6*b)*c[4])*d[0]
	 - 6*(35 + 14*(-5 + 2*b)*c[1] + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*c[4])*d[1]
	 - 2*(2*b*(-63 + 114*c[2] - 62*c[4]) + 3*(35 + 14*c[1] - 98*c[2] + 30*c[3]
	 + 38*c[4]))*d[2] + 2*(63 + 18*(-7 + 6*b)*c[1] - 90*c[2] + (306 - 212*b)*
	c[3] - 98*c[4])* d[3]) + 2*(2*b*(-99 + 682*c[2] - 1158*c[4]) + 11*(-21 + 
	78*c[1] - 114*c[2] - 98*c[3] + 310*c[4]))*d[4]) + 4*(-9009*c[3]*d[0] - 
	21879*c[3]*c[3]*d[0] + 15158*b*c[3]*c[3]*d[0] + 3003*c[4]*d[0] + 2574*b*
	c[4]*d[0] + 14014*c[3]*c[4]*d[0] - 22165*c[4]*c[4]*d[0] + 15054*b*c[4]*
	c[4]*d[0] + 18018*c[3]*d[1] - 15444*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] - 28028*c[3]*c[4]*d[1] + 12376*b*c[3]*c[4]*d[1] + 14222*
	c[4]*c[4]*d[1] + 12870*c[3]*d[2] + 16874*c[3]*c[3]*d[2] - 16692*b*c[3]*
	c[3]*d[2] + 16302*c[4]*d[2] - 17732*b*c[4]*d[2] + 3276*c[3]*c[4]*d[2] + 
	15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] - 43758*c[3]*d[3] + 30316*b*
	c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] + 21476*c[3]*c[4]*d[3] 
	- 12840*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] - 2*((-5369 + 3210*b)*c[3]
	*c[3] + 11*c[3]*(-637 + 934*c[4]) - c[4]*(-22165 + 2331*c[4] + 2*b*(7527 
	+ 967*c[4])))*d[4] + 13*c[2]*c[2]* (33*(-49 + 38*b)*d[0] + 22*(33*d[1] + 
	(81 - 94*b)*d[2] - 7*d[3]) + 2*(-737 + 678*b)*d[4]) + c[2]*(143*(3*(35 + 
	30*c[3] + 38*c[4]) - 2*b*(63 + 62*c[4]))* d[0] + 26*(11*(21 + (-90 + 52*b)
	*c[3] + 10*c[4])*d[1] + (-11*(147 + 14*c[3] + 134*c[4]) + 6*b*(209 + 226*
	c[4]))*d[2] + (2*(649 - 642*b)*c[3] + 9*(55 + 14*c[4]))*d[3]) + 2*(13*
	(627 - 682*b + 126*c[3]) - 94*(-169 + 190*b)*c[4])* d[4]) + 143*c[1]*c[1]*
	 (21*(-5 + 2*b)*d[0] - 2*(-63*d[1] + 3*(7 + 2*b)*d[2] + 9*d[3] + (-39 + 
	22*b)*d[4])) + 13*c[1]*(33*(35 + 14*c[2] + 6*(7 - 6*b)*c[3] - 26*c[4])*
	 d[0] + 2*(-11* (3*(35 + 14*c[2] + 6*c[3] - 26*c[4]) + 2*b*(-21 + 6*c[2] 
	+ 22*c[4]))*d[1] + 11*(21 + 66*c[2] + (-90 + 52*b)*c[3] + 10*c[4])*d[2] + 
	(11*(63 - 90*c[2] + b*(-54 + 52*c[2]) + 94*c[3]) + 14*(-77 + 34*b)*c[4])*
	d[3] + (-429 + 110*c[2] + 14*(-77 + 34*b)*c[3] + 1094*c[4])* d[4])))) + 
	b*(429*c[0]*c[0]*(35*(-3 + 2*b)*d[0] + 14*(5*(-1 + b)*d[1] + (5 - 2*b)*
	d[2] - 3*(-1 + b)*d[3]) - 2*(-7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(-15 - 10*
	c[1] + 10*c[2] + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) + 6*c[3]) - 2*(-7 + 
	10*b)*c[4])*d[0] + 22*(3*(-7*(5 + 10*c[1] + 2*c[2] - 6*c[3]) + b*(35 + 56*
	c[1] + 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[1] + (-3*(-35 + 14*c[1] 
	+ 98*c[2] + 30*c[3] - 38*c[4]) + 2*b*(-21 + 21*c[1] + 90*c[2] + 45*c[3] - 
	26*c[4]))* d[2] + (63 - 63*b + 126*c[1] - 72*b*c[1] - 90*c[2] + 90*b*c[2] 
	- 306*c[3] + 200*b*c[3] + 98*(-1 + b)*c[4])* d[3]) - 2*(11* (-3*(7 + 26*
	c[1] + 38*c[2]) + b*(30 + 78*c[1] + 52*c[2] - 98*c[3]) + 98*c[3]) - 2*
	(-1705 + 1126*b)*c[4])*d[4]) + 4*(9009*c[3]*d[0] - 9009*b*c[3]*d[0] - 
	21879*c[3]*c[3]*d[0] + 14300*b*c[3]*c[3]*d[0] + 3003*c[4]*d[0] - 4290*b*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] - 22165*c[4]*
	c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] + 18018*c[3]*d[1] - 10296*b*c[3]*d[1] 
	- 13442*c[3]*c[3]*d[1] + 13442*b*c[3]*c[3]*d[1] + 11154*c[4]*d[1] - 11154*
	b*c[4]*d[1] - 28028*c[3]*c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] - 14222*c[4]*
	c[4]*d[1] + 14222*b*c[4]*c[4]*d[1] - 12870*c[3]*d[2] + 12870*b*c[3]*d[2] 
	+ 16874*c[3]*c[3]*d[2] - 8528*b*c[3]*c[3]*d[2] + 16302*c[4]*d[2] - 7436*b*
	c[4]*d[2] - 3276*c[3]*c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] + 15886*c[4]*c[4]*
	d[2] - 6956*b*c[4]*c[4]*d[2] - 43758*c[3]*d[3] + 28600*b*c[3]*d[3] + 13806
	*c[3]*c[3]*d[3] - 13806*b*c[3]*c[3]*d[3] - 14014*c[4]*d[3] + 14014*b*c[4]*
	d[3] + 21476*c[3]*c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*
	d[3] - 10274*b*c[4]*c[4]*d[3] - 2*((-5369 + 3764*b)*c[3]*c[3] + 11*(-1+b)
	*c[3]*(-637 + 934*c[4]) + c[4]*(22165 - 2331*c[4]+2*b*(-7319+1649*c[4])))
	*d[4] + 143*c[1]*c[1]* (3*(7*(-5 + 4*b)*d[0] + 2*(21*(-1 + b)*d[1] + (-7 
	+ 8*b)*d[2] - 3*(-1 + b)*d[3])) - 2*(-39 + 28*b)*d[4]) + 13*c[2]*c[2]* 
	(33*(-49 + 30*b)*d[0] + 22*(33*(-1 + b)*d[1] + (81 - 34*b)*d[2] - 7*(-1 + 
	b)*d[3]) + 2*(-737 + 398*b)*d[4]) + c[2]*(13*(11*(105 - 90*c[3] + b*(-42 
	+ 90*c[3] - 52*c[4]) + 114*c[4])*d[0] + 2*(11*(-21 + 21*b - 90*c[3] + 64*
	b*c[3] + 10*(-1 + b)*c[4])*d[1] + (-1617 + 990*b + 154*c[3] - 154*b*c[3] 
	- 1474*c[4] + 796*b*c[4])*d[2] + (-495 + 495*b + 1298*c[3] - 656*b*c[3] + 
	126*(-1 + b)*c[4])*d[3])) + 2*(8151 - 3718*b + 1638*(-1 + b)*c[3] - 94*
	(-169 + 74*b)*c[4])*d[4]) + 13*c[1]*(33*(-35 - 14*c[2] + 42*c[3] + b*(35 
	+ 14*c[2] - 24*c[3] - 26*c[4]) + 26*c[4])*d[0] - 462*(5*d[1] + d[2] - 3*
	d[3]) + 858*d[4] + 2*(-22*(3*c[2]*(7*d[1] + 11*d[2] + 15*d[3]) + c[3]*
	(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*d[2] + 49*d[3])) - 2*
	(55*c[2] + 539*c[3] + 547*c[4])*d[4] + b*(22*(42 + 24*c[2] - 9*c[3] - 28*
	c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3] + 10*c[4])*d[2] + 2*(-198 + 352*
	c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 + 110*c[2] + 840*c[3] + 
	1094*c[4])*d[4])))))))/ 2882880.;
		bf_mom[2] = -((a - b)*(a*(429*c[0]*c[0]* (35*(3 + 2*b)*d[0] - 70*
	d[1] - 14*(5 + 6*b)*d[2] + 42*d[3] + 2*(-7 + 6*b)*d[4]) - 26*c[0]*(33*(2*
	b*(-35 + 42*c[2] - 6*c[4]) + 7*(-15 + 10*c[1] + 10*c[2] - 6*c[3]+2*c[4]))*
	d[0] + 22*(-3*(-35 + 14*(5 + 2*b)*c[1] - 14*c[2] - 6*(7 + 6*b)*c[3] + 26*
	c[4])*d[1] + (3*(35 + 14*c[1] + b*(42 - 76*c[2]) - 98*c[2] + 30*c[3]) + 2*
	(57 + 62*b)*c[4])*d[2] + (-63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 
	106*b)*c[3] + 98*c[4])*d[3]) + 2*(11*(21 - 78*c[1] + 114*c[2] + 2*b*(-9 + 
	62*c[2]) + 98*c[3]) - 2*(1705 + 1158*b)*c[4])*d[4]) + 4*(9009*c[3]*d[0] +
	 21879*c[3]*c[3]*d[0] + 15158*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] + 2574*b*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 15054*b*c[4]*
	c[4]*d[0] - 18018*c[3]*d[1] - 15444*b*c[3]*d[1] - 13442*c[3]*c[3]*d[1] + 
	11154*c[4]*d[1] + 28028*c[3]*c[4]*d[1] + 12376*b*c[3]*c[4]*d[1] - 14222*
	c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 16692*b*c[3]*
	c[3]*d[2] - 16302*c[4]*d[2] - 17732*b*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 
	15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] + 43758*c[3]*d[3] + 30316*b*
	c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 21476*c[3]*c[4]*d[3]
	 - 12840*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 2*((5369 + 3210*b)*
	c[3]*c[3] - 11*c[3]*(-637 + 934*c[4]) - c[4]*(22165 - 2331*c[4] + 2*b*
	(7527 + 967*c[4])))*d[4] + 13*c[2]*c[2]* (33*(49 + 38*b)*d[0] - 22*(33*
	d[1] + (81 + 94*b)*d[2] - 7*d[3]) + 2*(737 + 678*b)*d[4]) + c[2]*(-143*(3*
	(35 + 42*b + 30*c[3]) + 2*(57 + 62*b)*c[4])* d[0] + 26*(11*(-21 + (90 + 
	52*b)*c[3] - 10*c[4])*d[1] + (11*(147 + 114*b + 14*c[3]) + 2*(737+678*b)*
	c[4])* d[2] - (2*(649 + 642*b)*c[3] + 9*(55 + 14*c[4]))*d[3]) - 2*(13*
	(627 + 682*b + 126*c[3]) + 94*(169 + 190*b)*c[4])*d[4]) + 143*c[1]*c[1]*
	 (21*(5 + 2*b)*d[0] - 2*(63*d[1] + 3*(-7 + 2*b)*d[2] - 9*d[3] + (39+22*b)
	*d[4])) + 13*c[1]*(-33*(35 + 14*c[2] + 6*(7 + 6*b)*c[3] - 26*c[4])* d[0] 
	+ 2*(-11* (-3*(35 + 14*c[2] + 6*c[3] - 26*c[4]) + 2*b*(-21 + 6*c[2] + 22*
	c[4]))*d[1] - 11*(21 + 66*c[2] - 2*(45 + 26*b)*c[3] + 10*c[4])*d[2] + (11*
	(-63 - 54*b + 90*c[2] + 52*b*c[2] - 94*c[3]) + 14*(77 + 34*b)*c[4])*d[3] 
	+ 429*d[4] - 2*(55*c[2] - 7*(77 + 34*b)*c[3] + 547*c[4])*d[4])))) + 2*a*a*
	(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] - 2*d[2] + 3*d[3]) - 10*d[4]) + 26*
	c[0]*(-33*(7*(-5 + 5*c[1] + 2*c[2] - 3*c[3]) + 10*c[4])* d[0] + 33*(-35 
	+ 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])* d[1] - 22*(21 + 21*c[1] - 90*
	c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(-63 + 72*c[1] + 90*c[2] - 200*c[3] 
	+ 98*c[4])* d[3] + 2*(429*c[1] - 11*(15 + 26*c[2] + 49*c[3]) + 1126*c[4])*
	d[4]) + 2*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] - 14014*
	c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*d[1] - 13442*c[3]*c[3]*
	d[1] + 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 
	12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*d[2] - 3276*c[3]*c[4]*
	d[2] - 6956*c[4]*c[4]*d[2] + 28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 
	14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*
	(11*(21 + 45*c[3] + 26*c[4])*d[0] + (231 - 704*c[3] + 110*c[4])*d[1] - 2*
	(495 + 77*c[3] + 398*c[4])*d[2] + (495 + 656*c[3] + 126*c[4])*d[3]) + 286*
	c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*
	(1859 + 819*c[3] + 3478*c[4])*d[4] - 2*(c[3]*(7007 + 3764*c[3]) - 2*(7319
	 + 5137*c[3])*c[4] + 3298*c[4]*c[4])*d[4] + 26*c[2]*c[2]* (495*d[0] - 363*
	d[1] - 374*d[2] + 77*d[3] + 398*d[4]) - 13*c[1]*(33*(35 + 14*c[2]+24*c[3] 
	- 26*c[4])*d[0] + 2*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(21 +
	 66*c[2] - 64*c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3] + 420*
	c[4])*d[3] + (-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4])))) + b*(429*
	c[0]*c[0]*(35*(3 + 2*b)*d[0] + 70*(1 + b)*d[1] - 14*(5 + 2*b)*d[2] - 42*
	(1 + b)*d[3] - 2*(7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(15 + 10*c[1] - 10*c[2]
	 + 2*b*(5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3]) - 2*(7 + 10*b)*c[4])*d[0]
	 + 22*(3*(35 + 35*b + 70*c[1] + 56*b*c[1] + 14*c[2] + 14*b*c[2] - 42*c[3] 
	- 24*b*c[3] - 26*(1 + b)*c[4])* d[1] + (3*(-35 - 14*b + 14*c[1] + 14*b*
	c[1] + 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3]) - 2*(57 + 26*b)*c[4])*d[2] 
	+ (-63 - 63*b - 126*c[1] - 72*b*c[1] + 90*c[2] + 90*b*c[2] + 306*c[3] + 
	200*b*c[3] + 98*(1 + b)*c[4])* d[3]) - 2*(11* (21 + 30*b + 78*c[1] + 78*b*
	c[1] + 114*c[2] + 52*b*c[2] - 98*(1 + b)*c[3]) - 2*(1705 + 1126*b)*c[4])*
	 d[4]) + 4*(-9009*c[3]*d[0] - 9009*b*c[3]*d[0] + 21879*c[3]*c[3]*d[0] + 
	14300*b*c[3]*c[3]*d[0] - 3003*c[4]*d[0] - 4290*b*c[4]*d[0] + 14014*c[3]*
	c[4]*d[0] + 14014*b*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 14638*b*c[4]*
	c[4]*d[0] - 18018*c[3]*d[1] - 10296*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 
	13442*b*c[3]*c[3]*d[1] - 11154*c[4]*d[1] - 11154*b*c[4]*d[1] + 28028*c[3]*
	c[4]*d[1] + 21840*b*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 14222*b*c[4]*
	c[4]*d[1] + 12870*c[3]*d[2] + 12870*b*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 
	8528*b*c[3]*c[3]*d[2] - 16302*c[4]*d[2] - 7436*b*c[4]*d[2] + 3276*c[3]*
	c[4]*d[2] + 3276*b*c[3]*c[4]*d[2] - 15886*c[4]*c[4]*d[2] - 6956*b*c[4]*
	c[4]*d[2] + 43758*c[3]*d[3] + 28600*b*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 
	13806*b*c[3]*c[3]*d[3] + 14014*c[4]*d[3] + 14014*b*c[4]*d[3] - 21476*c[3]*
	c[4]*d[3] - 15056*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] - 10274*b*c[4]*
	c[4]*d[3] - 2*((5369 + 3764*b)*c[3]*c[3] + 11*(1 + b)*c[3]*(-637 + 934*
	c[4]) + c[4]*(-22165 + 2331*c[4] + 2*b*(-7319 + 1649*c[4])))*d[4] + 143*
	c[1]*c[1]* (21*(5 + 4*b)*d[0] + 126*(1 + b)*d[1] + 6*(7 + 8*b)*d[2] - 18*
	(1 + b)*d[3] - 2*(39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*(49 + 30*b)*d[0] + 
	726*(1 + b)*d[1] - 22*(81 + 34*b)*d[2] - 154*(1 + b)*d[3] + 2*(737+398*b)*
	d[4]) + c[2]*(143*(b*(-42 + 90*c[3] - 52*c[4]) + 3*(-35 + 30*c[3] - 38*
	c[4]))*d[0] + 286*(21 + 21*b + 90*c[3] + 64*b*c[3] + 10*(1 + b)*c[4])* 
	d[1] - 26*(-33*(49 + 30*b) + 154*(1 + b)*c[3] - 2*(737 + 398*b)*c[4])*d[2]
	 - 26*(-495 - 495*b + 1298*c[3] + 656*b*c[3] - 126*(1 + b)*c[4])*d[3] + 2*
	(-143*(57 + 26*b) + 1638*(1 + b)*c[3] - 94*(169 + 74*b)*c[4])*d[4]) + 13*
	c[1]*(33*(35 + 35*b + 14*c[2] + 14*b*c[2] - 42*c[3] - 24*b*c[3] - 26*(1+b)
	*c[4])*d[0] + 462*(5*d[1] + d[2] - 3*d[3]) - 858*d[4] + 2*(22*(3*c[2]*(7*
	d[1] + 11*d[2] + 15*d[3]) + c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*
	(-39*d[1] + 5*d[2] + 49*d[3])) + 2*(55*c[2] + 539*c[3] + 547*c[4])*d[4] + 
	b*(22*(42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3]
	 + 10*c[4])*d[2] + 2*(-198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429
	 + 110*c[2] + 840*c[3] + 1094*c[4])*d[4])))))))/ 2882880.;
		bf_mom[3] = -((a - b)*(2*a*a*(429*c[0]*c[0]* (7*(5*d[0] - 5*d[1] 
	- 2*d[2] + 3*d[3]) - 10*d[4]) + 26*c[0]*(-33*(7*(5 + 5*c[1] + 2*c[2] - 3*
	c[3]) + 10*c[4])* d[0] + 33*(35 + 56*c[1] - 14*c[2] - 24*c[3] + 26*c[4])*
	 d[1] - 22*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] - 11*(63 + 
	72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])* d[3] + 2*(165 + 429*c[1] - 286*
	c[2] - 539*c[3] + 1126*c[4])* d[4]) + 2*(-9009*c[3]*d[0] + 14300*c[3]*c[3]
	*d[0] + 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 
	10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 21840*c[3]*c[4]
	*d[1] - 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] + 
	7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*
	d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 
	10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] - 11*
	(21 + 64*c[3] - 10*c[4])*d[1] - 2*(-495 + 77*c[3] + 398*c[4])*d[2] + 
	(-495 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 
	24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 3478*c[4])*d[4] 
	- 2*(3764*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319+1649*c[4]))
	*d[4] + 26*c[2]*c[2]* (495*d[0] - 363*d[1] - 374*d[2] + 77*d[3]+398*d[4]) 
	- 13*c[1]*(33*(-35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 
	24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])
	*d[2] - 2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] - 
	840*c[3] + 1094*c[4])*d[4])))) + a*(429*c[0]*c[0]*(35*(3 + 2*b)*d[0] - 70*
	d[1] - 14*(5 + 6*b)*d[2] + 42*d[3] + 2*(-7 + 6*b)*d[4]) - 26*c[0]*(33*(2*
	b*(35 + 42*c[2] - 6*c[4]) + 7*(15 + 10*c[1] + 10*c[2] - 6*c[3] + 2*c[4]))
	*d[0] + 22*(-3*(35 + 14*(5 + 2*b)*c[1] - 14*c[2] - 6*(7 + 6*b)*c[3] + 26*
	c[4])*d[1] + (-3*(35 - 14*c[1] + 98*c[2] + b*(42 + 76*c[2]) - 30*c[3]) + 
	2*(57 + 62*b)*c[4])*d[2] + (63 + 18*(7 + 6*b)*c[1] + 90*c[2] - 2*(153 + 
	106*b)*c[3] + 98*c[4])*d[3]) + 2*(11*(-21 - 78*c[1] + 114*c[2] + 2*b*(9 + 
	62*c[2]) + 98*c[3]) - 2*(1705 + 1158*b)*c[4])*d[4]) + 4*(-9009*c[3]*d[0] 
	+ 21879*c[3]*c[3]*d[0] + 15158*b*c[3]*c[3]*d[0] + 3003*c[4]*d[0] - 2574*b*
	c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 15054*b*c[4]*
	c[4]*d[0] + 18018*c[3]*d[1] + 15444*b*c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 
	11154*c[4]*d[1] + 28028*c[3]*c[4]*d[1] + 12376*b*c[3]*c[4]*d[1] - 14222*
	c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 16692*b*c[3]*
	c[3]*d[2] + 16302*c[4]*d[2] + 17732*b*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 
	15886*c[4]*c[4]*d[2] - 17860*b*c[4]*c[4]*d[2] - 43758*c[3]*d[3] - 30316*b*
	c[3]*d[3] + 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 21476*c[3]*c[4]*d[3] 
	- 12840*b*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 2*((5369 + 3210*b)*c[3]*
	c[3] - 11*c[3]*(637 + 934*c[4]) + c[4]*(22165 + 2*b*(7527 - 967*c[4]) + 
	2331*c[4]))*d[4] + 13*c[2]*c[2]* (33*(49 + 38*b)*d[0] - 22*(33*d[1] + 
	(81 + 94*b)*d[2] - 7*d[3]) + 2*(737 + 678*b)*d[4]) + c[2]*(-143*(-21*(5 + 
	6*b) + 90*c[3] + 2*(57 + 62*b)*c[4])* d[0] + 26*(11*(21 + (90 + 52*b)*
	c[3] - 10*c[4])*d[1] + (-33*(49 + 38*b) + 154*c[3] + 2*(737 + 678*b)*c[4])
	* d[2] + 495*d[3] - 2*((649 + 642*b)*c[3] + 63*c[4])*d[3]) - 2*(-13*(627 
	+ 682*b - 126*c[3]) + 94*(169 + 190*b)*c[4])* d[4]) + 143*c[1]*c[1]* (21*
	(5 + 2*b)*d[0] - 2*(63*d[1] + 3*(-7 + 2*b)*d[2] - 9*d[3] + (39 + 22*b)*
	d[4])) + 13*c[1]*(-33*(-35 + 14*c[2] + 6*(7 + 6*b)*c[3] - 26*c[4])* d[0] 
	+ 2*(-11* (105 - 42*c[2] - 18*c[3] + 78*c[4] + 2*b*(21 + 6*c[2]+22*c[4]))
	*d[1] - 11*(-21 + 66*c[2] - 2*(45 + 26*b)*c[3] + 10*c[4])* d[2] + (11* (63
	 + 54*b + 90*c[2] + 52*b*c[2] - 94*c[3]) + 14*(77 + 34*b)*c[4])*d[3] - 
	(429 + 110*c[2] - 14*(77 + 34*b)*c[3] + 1094*c[4])*d[4])))) + b*(429*c[0]*
	c[0]*(35*(3 + 2*b)*d[0] + 70*(1 + b)*d[1] - 14*(5 + 2*b)*d[2] - 42*(1 + b)
	*d[3] - 2*(7 + 10*b)*d[4]) + 26*c[0]*(33*(7*(-15 + 10*c[1] - 10*c[2] + 2*
	b*(-5 + 5*c[1] - 2*c[2] - 3*c[3]) - 6*c[3]) - 2*(7 + 10*b)*c[4])*d[0] + 
	22*(3*(-35 - 35*b + 70*c[1] + 56*b*c[1] + 14*c[2] + 14*b*c[2] - 42*c[3] - 
	24*b*c[3] - 26*(1 + b)*c[4])* d[1] + (3*(35 + 14*b + 14*c[1] + 14*b*c[1] 
	+ 98*c[2] + 60*b*c[2] + 30*(1 + b)*c[3]) - 2*(57 + 26*b)*c[4])*d[2] + (63 
	+ 63*b - 126*c[1] - 72*b*c[1] + 90*c[2] + 90*b*c[2] + 306*c[3] + 200*b*
	c[3] + 98*(1 + b)*c[4])* d[3]) - 2*(11* (-21 - 30*b + 78*c[1] + 78*b*c[1] 
	+ 114*c[2] + 52*b*c[2] - 98*(1 + b)*c[3]) - 2*(1705 + 1126*b)*c[4])* d[4])
	 + 4*(9009*c[3]*d[0] + 9009*b*c[3]*d[0] + 21879*c[3]*c[3]*d[0] + 14300*b*
	c[3]*c[3]*d[0] + 3003*c[4]*d[0] + 4290*b*c[4]*d[0] + 14014*c[3]*c[4]*d[0] 
	+ 14014*b*c[3]*c[4]*d[0] + 22165*c[4]*c[4]*d[0] + 14638*b*c[4]*c[4]*d[0] 
	+ 18018*c[3]*d[1] + 10296*b*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 13442*b*
	c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 11154*b*c[4]*d[1] + 28028*c[3]*c[4]*
	d[1] + 21840*b*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 14222*b*c[4]*c[4]*
	d[1] - 12870*c[3]*d[2] - 12870*b*c[3]*d[2] - 16874*c[3]*c[3]*d[2] - 8528*
	b*c[3]*c[3]*d[2] + 16302*c[4]*d[2] + 7436*b*c[4]*d[2] + 3276*c[3]*c[4]*
	d[2] + 3276*b*c[3]*c[4]*d[2] - 15886*c[4]*c[4]*d[2] - 6956*b*c[4]*c[4]*
	d[2] - 43758*c[3]*d[3] - 28600*b*c[3]*d[3] - 13806*c[3]*c[3]*d[3] - 13806*
	b*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 14014*b*c[4]*d[3] - 21476*c[3]*c[4]*
	d[3] - 15056*b*c[3]*c[4]*d[3] - 10274*c[4]*c[4]*d[3] - 10274*b*c[4]*c[4]*
	d[3] - 2*((5369 + 3764*b)*c[3]*c[3] + 11*(1 + b)*c[3]*(637 + 934*c[4]) + 
	c[4]*(22165 + 2331*c[4] + 2*b*(7319 + 1649*c[4])))*d[4] + 143*c[1]*c[1]* 
	(21*(5 + 4*b)*d[0] + 126*(1 + b)*d[1] + 6*(7 + 8*b)*d[2] - 18*(1 + b)*d[3]
	 - 2*(39 + 28*b)*d[4]) + 13*c[2]*c[2]* (33*(49 + 30*b)*d[0] + 726*(1 + b)
	*d[1] - 22*(81 + 34*b)*d[2] - 154*(1 + b)*d[3] + 2*(737 + 398*b)*d[4]) + 
	c[2]*(143*(b*(42 + 90*c[3] - 52*c[4]) + 3*(35 + 30*c[3] - 38*c[4]))*d[0] 
	+ 286*(-21 - 21*b + 90*c[3] + 64*b*c[3] + 10*(1 + b)*c[4])* d[1] - 26*(33*
	(49 + 30*b) + 154*(1 + b)*c[3] - 2*(737 + 398*b)*c[4])*d[2] - 26*(495 + 
	495*b + 1298*c[3] + 656*b*c[3] - 126*(1 + b)*c[4])*d[3] + 2*(143*(57 + 26*
	b) + 1638*(1 + b)*c[3] - 94*(169 + 74*b)*c[4])*d[4]) + 13*c[1]*(33*(-35 - 
	35*b + 14*c[2] + 14*b*c[2] - 42*c[3] - 24*b*c[3] - 26*(1 + b)*c[4])*d[0] 
	- 462*(5*d[1] + d[2] - 3*d[3]) + 858*d[4] + 2*(22*(3*c[2]*(7*d[1] + 11*
	d[2] + 15*d[3]) + c[3]*(-9*d[1] + 45*d[2] + 47*d[3]) + c[4]*(-39*d[1] + 5*
	d[2] + 49*d[3])) + 2*(55*c[2] + 539*c[3] + 547*c[4])*d[4] + b*(22*(-42 + 
	24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])*
	d[2] + 2*(198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] 
	+ 840*c[3] + 1094*c[4])*d[4])))))))/ 2882880.;
		bf_mom[4] = ((a - b)*(2*a*a*(13*(33* (7*(-20 + 5*c[0]*c[0] - 10*
	c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*
	(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1]
	 + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4])*d[0] + 1716*(35*d[1] + 14*
	d[2] - 21*d[3] + 10*d[4]) - 429*c[0]*c[0]*(35*d[1] + 14*d[2] - 21*d[3] + 
	10*d[4]) + 52*c[0]*(-11*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + c[4]*(-39*
	d[1] + 26*d[2] + 49*d[3])) + (-539*c[3] + 1126*c[4])*d[4] + 33*c[1]*(28*
	d[1] - 7*d[2] - 12*d[3] + 13*d[4]) - 11*c[2]*(21*d[1] - 90*d[2] + 45*d[3]
	 + 26*d[4])) + 4*(-6721*c[3]*c[3]*d[1] + 10920*c[3]*c[4]*d[1] - 7111*c[4]*
	c[4]*d[1] - 4264*c[3]*c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 3478*c[4]*c[4]*
	d[2] + 26*c[2]*(c[3]*(352*d[1] + 77*d[2] - 328*d[3]) + c[4]*(-55*d[1] + 
	398*d[2] - 63*d[3])) + 6903*c[3]*c[3]*d[3] - 7528*c[3]*c[4]*d[3] + 5137*
	c[4]*c[4]*d[3] - 13*c[2]*c[2]*(363*d[1] + 374*d[2] - 77*d[3] - 398*d[4]) 
	- 2*c[2]*(819*c[3] + 3478*c[4])*d[4] - 2*(1882*c[3]*c[3] - 5137*c[3]*c[4]
	 + 1649*c[4]*c[4])* d[4] - 143*c[1]*c[1]* (63*d[1] - 24*d[2] - 9*d[3] + 
	28*d[4]) + 26*c[1]*(11*c[2]*(24*d[1] - 33*d[2] + 32*d[3] - 5*d[4]) + c[3]*
	(99*d[1] + 352*d[2] - 517*d[3] + 420*d[4]) - c[4]*(308*d[1] + 55*d[2] - 
	420*d[3] + 547*d[4])))) + b*(13*(33*(7*(5*(-3 + 2*b)*(-4 + c[0]*c[0]) + 
	20*(-1 + b)*c[0]*c[1] + 4*(-5 + 4*b)*c[1]*c[1]) - 28*((-5 + 2*b)*c[0] - 
	2*(-1 + b)*c[1])*c[2] + 4*(-49 + 30*b)*c[2]*c[2]) - 396*(7*(-1 + b)*c[0] 
	+ 2*(-7 + 4*b)*c[1] - 10*(-1 + b)*c[2])*c[3] + 44*(-153 + 100*b)*c[3]*c[3]
	 - 44*(3*(-7 + 10*b)*c[0] + 78*(-1 + b)*c[1] + 2*(-57 + 26*b)*c[2] - 98*
	(-1 + b)*c[3])*c[4] + 4*(-1705 + 1126*b)*c[4]*c[4])*d[0] + 24024*(5*d[1] 
	- 5*d[2] - 3*d[3] - d[4]) + 2*(-3003*c[0]*c[0]*(5*d[1] - 5*d[2] - 3*d[3] 
	- d[4]) - 572*c[0]*(-63*c[3]*d[1] - 39*c[4]*d[1] + 45*c[3]*d[2] - 57*c[4]*
	d[2] + 153*c[3]*d[3] + 49*c[4]*d[3] + 3*c[2]*(7*d[1] + 49*d[2] + 15*d[3] 
	- 19*d[4]) + 3*c[1]*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) + 49*c[3]*d[4]
	 + 155*c[4]*d[4]) + b*(13*(33*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] + 64*
	c[1]*c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2])) - 88*(18*c[0] + 9*c[1]
	 - 32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 616*c[1] - 110*c[2] - 
	840*c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 + 21*c[0]*c[0] - 48*
	c[1]*c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*(7*c[1] + 30*c[2])) - 
	286*(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*(143*c[0] - 
	55*c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (13*(99*(-28 
	+ (c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2] + 308*
	c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3]) - 
	4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])* c[4] + 20548*c[4]*c[4])
	*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 539*c[3] - 
	1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*c[2]*c[3] + 
	1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*c[4] - 13*
	c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4]) - 4*(6721*c[3]*c[3]*d[1] + 
	14014*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] - 8437*c[3]*c[3]*d[2] + 1638*
	c[3]*c[4]*d[2] - 7943*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] - 10738*c[3]*
	c[4]*d[3] - 5137*c[4]*c[4]*d[3] + 429*c[1]*c[1]*(21*d[1] + 7*d[2] - 3*
	d[3] - 13*d[4]) - (5369*c[3]*c[3] + 10274*c[3]*c[4] + 2331*c[4]*c[4])*
	d[4] + 143*c[2]*c[2]*(33*d[1] - 81*d[2] - 7*d[3] + 67*d[4]) + 26*c[2]*
	(c[4]*(55*d[1] + 737*d[2] + 63*d[3] - 611*d[4]) + c[3]*(495*d[1] - 77*d[2]
	 - 649*d[3] + 63*d[4])) + 26*c[1]*(11*c[2]*(21*d[1] + 33*d[2] + 45*d[3] 
	+ 5*d[4]) + 11*c[3]*(-9*d[1] + 45*d[2] + 47*d[3] + 49*d[4]) + c[4]*(-429*
	d[1] + 55*d[2] + 539*d[3] + 547*d[4]))))) + a*(13*(33*(7*(5*(-3 + 2*b)*
	(-4 + c[0]*c[0]) + 20*c[0]*c[1] + 4*(-5 + 2*b)*c[1]*c[1]) - 28*((-5+6*b)*
	c[0] - 2*c[1])*c[2] + 4*(-49 + 38*b)*c[2]*c[2]) - 396*(7*c[0] + 2*(-7 + 
	6*b)*c[1] - 10*c[2])*c[3] + 44*(-153 + 106*b)*c[3]*c[3] + 44*(3*(7 + 6*b)
	*c[0] - 78*c[1] + 2*(57 - 62*b)*c[2] + 98*c[3])*c[4] + 4*(-1705 + 1158*b)
	*c[4]*c[4])*d[0] - 24024*(5*d[1] + 5*d[2] - 3*d[3] + d[4]) + 2*(429*c[0]*
	c[0]* (7*(5*d[1] + (5 - 6*b)*d[2] - 3*d[3]) + (7 + 6*b)*d[4]) + 52*c[0]*
	(11*(63*c[3]*d[1] - 54*b*c[3]*d[1] - 39*c[4]*d[1] + 45*c[3]*d[2] + 57*c[4]
	*d[2] - 62*b*c[4]*d[2] + ((-153 + 106*b)*c[3] + 49*c[4])*d[3] + 3*c[2]*
	(7*d[1] + (-49 + 38*b)*d[2] + 15*d[3]) + 3*c[1]*(7*(-5 + 2*b)*d[1] + 7*
	d[2] + 3*(7 - 6*b)*d[3])) - (429*c[1] + 11*(-57 + 62*b)*c[2] - 539*c[3] +
	 (1705 - 1158*b)*c[4])*d[4]) + 4*(6721*c[3]*c[3]*d[1] - 14014*c[3]*c[4]*
	d[1] + 6188*b*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] + 18018*b*d[2] + 8437*
	c[3]*c[3]*d[2] - 8346*b*c[3]*c[3]*d[2] + 1638*c[3]*c[4]*d[2] + 7943*c[4]*
	c[4]*d[2] - 8930*b*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] + 10738*c[3]*c[4]*
	d[3] - 6420*b*c[3]*c[4]*d[3] - 5137*c[4]*c[4]*d[3] - (-5369*c[3]*c[3] + 
	6*b*(429 + 535*c[3]*c[3]) + 10274*c[3]*c[4] - (2331 + 1934*b)*c[4]*c[4])*
	d[4] + 143*c[1]*c[1]* (63*d[1] - 3*(7 + 2*b)*d[2] - 9*d[3] + (39 - 22*b)*
	d[4]) + 13*c[2]*c[2]* (11*(33*d[1] + (81 - 94*b)*d[2] - 7*d[3]) + (-737 +
	 678*b)*d[4]) + 2*c[2]*(13*c[3]* (11*(-45 + 26*b)*d[1] - 77*d[2] + (649 
	- 642*b)*d[3] + 63*d[4]) + c[4]*(715*d[1] + 13*(-737 + 678*b)*d[2] + 819*
	d[3] + 47*(169 - 190*b)*d[4])) + 26*c[1]*(-11*c[2]* (3*(7 + 2*b)*d[1] - 
	33*d[2] + (45 - 26*b)*d[3] - 5*d[4]) + c[4]*((429 - 242*b)*d[1] + 55*d[2] 
	+ 7*(-77 + 34*b)*d[3] + 547*d[4]) + c[3]*(-99*d[1] + 11*(-45 + 26*b)*d[2]
	 + 517*d[3] + 7*(-77 + 34*b)*d[4])))))))/1441440.;
		bf_mom[5] = ((a - b)*(429*c[0]*c[0]*(35*(-3 + a*a + a*b + b*b)*
	 d[0] - 7*(-10*d[2] + 6*a*b*d[2] + a*a*(5*d[1] + 2*d[2] - 3*d[3]) + b*b*
	(-5*d[1] + 2*d[2] + 3*d[3])) - 2*(-7 + 5*a*a - 3*a*b + 5*b*b)*d[4]) - 
	26*c[0]*(33*(7*(15 - 10*c[2] + a*b*(-5 + 6*c[2]) + a*a*(-5 + 5*c[1] + 2*
	c[2] - 3*c[3]) + b*b*(-5 - 5*c[1] + 2*c[2] + 3*c[3])) + 2*(-7 + 5*a*a - 
	3*a*b + 5*b*b)*c[4])*d[0] + 11*(b*b*(-3*(35 + 56*c[1] + 14*c[2] - 24*c[3]
	 - 26*c[4])*d[1] - 2*(-21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])*d[2] + 
	(63 + 72*c[1] - 90*c[2] - 200*c[3] - 98*c[4])*d[3]) + 6*((-35 + 98*c[2] -
	 38*c[4])*d[2] + 14*c[1]*(5*d[1] - 3*d[3]) + 6*c[3]*(-7*d[1] + 17*d[3])))
	 + 2*(-33*(7 + 38*c[2]) + b*b*(165 + 429*c[1] + 286*c[2] - 539*c[3] - 
	1126*c[4]) + 3410*c[4])*d[4] + 2*a*b*(11*((63 - 114*c[2] + 62*c[4])*d[2] 
	+ 2*c[3]*(27*d[1] - 53*d[3])) - 66*c[1]*(7*d[1] - 9*d[3]) + (-99+682*c[2] 
	- 1158*c[4])*d[4]) + a*a*(11*((105 - 168*c[1] + 42*c[2] + 72*c[3] - 78*
	c[4])* d[1] + 2*(21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2] + (-63
	 + 72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])*d[3]) - 2*(429*c[1] - 11*(15 + 
	26*c[2] + 49*c[3]) + 1126*c[4])*d[4])) + 2*(26*(-1683*c[3]*c[3]*d[0] + 
	231*c[4]*d[0] - 1705*c[4]*c[4]*d[0] + 1386*c[3]*d[1] - 2156*c[3]*c[4]*d[1]
	 - 33*c[2]*c[2]*(49*d[0] - 54*d[2]) + 1298*c[3]*c[3]*d[2] + 1254*c[4]*d[2]
	 + 1222*c[4]*c[4]*d[2] - 231*c[1]*c[1]*(5*d[0] + 2*d[2]) + 2*c[3]*(-1683 +
	 826*c[4])*d[3] + 22*c[1]*(-3*(35 + 14*c[2] - 26*c[4])*d[1] + 9*c[3]*(7*
	d[0] - 10*d[2]) + (63 - 90*c[2] - 98*c[4])*d[3]) + 11*c[2]*(3*(35 + 38*
	c[4])*d[0] - 294*d[2] - 268*c[4]*d[2] + 4*c[3]*(-45*d[1] + 59*d[3]))) + 4*
	(13*(429*c[1]*c[1] + 11*(57 - 67*c[2])*c[2] - 1078*c[1]*c[3] + 413*c[3]*
	c[3]) + 13*(-1705 + 1222*c[2])*c[4] + 2331*c[4]*c[4])*d[4] + 2*a*b*(7579*
	c[3]*c[3]*d[0] + 1287*c[4]*d[0] + 7527*c[4]*c[4]*d[0] - 7722*c[3]*d[1] + 
	6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] - 8866*c[4]*d[2] - 8930*c[4]*
	c[4]*d[2] + 15158*c[3]*d[3] - 6420*c[3]*c[4]*d[3] - 26*c[1]*(11*(-21 + 6*
	c[2] + 22*c[4])*d[1] + (297 - 286*c[2] - 238*c[4])*d[3] + c[3]*(297*d[0] 
	- 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*(21*d[0] - 6*d[2] - 22*d[4]) - 2*
	(1605*c[3]*c[3] - c[4]*(7527 + 967*c[4]))*d[4] + 13*c[2]*c[2]*(627*d[0] -
	 1034*d[2] + 678*d[4]) + c[2]*(-143*(63 + 62*c[4])*d[0] + 78*(209 + 226*
	c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) - 2*(4433 + 8930*c[4])*d[4])) 
	+ a*a*(9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] - 14014*c[3]
	*c[4]*d[0] + 14638*c[4]*c[4]*d[0] - 10296*c[3]*d[1] - 13442*c[3]*c[3]*d[1]
	 + 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1] - 14222*c[4]*c[4]*d[1] - 12870*
	c[3]*d[2] - 8528*c[3]*c[3]*d[2] - 7436*c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 
	6956*c[4]*c[4]*d[2] + 28600*c[3]*d[3] + 13806*c[3]*c[3]*d[3] - 14014*c[4]*
	d[3] - 15056*c[3]*c[4]*d[3] + 10274*c[4]*c[4]*d[3] - 26*c[2]*(11*(21 + 
	45*c[3] + 26*c[4])*d[0] + (231 - 704*c[3] + 110*c[4])*d[1] - 2*(495 + 77*
	c[3] + 398*c[4])*d[2] + (495 + 656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]*
	 (42*d[0] - 63*d[1] + 24*d[2] + 9*d[3] - 28*d[4]) - 4*c[2]*(1859+819*c[3] 
	+ 3478*c[4])*d[4] - 2*(c[3]*(7007 + 3764*c[3]) - 2*(7319 + 5137*c[3])*c[4]
	 + 3298*c[4]*c[4])*d[4] + 26*c[2]*c[2]*(495*d[0] - 363*d[1] - 374*d[2] + 
	77*d[3] + 398*d[4]) - 13*c[1]* (33*(35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0]
	 + 2*(-22*(42 + 24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] - 64*
	c[3] + 10*c[4])*d[2] - 2*(-198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + 
	(-429 + 110*c[2] - 840*c[3] + 1094*c[4])*d[4]))) + b*b*(-9009*c[3]*d[0] + 
	14300*c[3]*c[3]*d[0] - 4290*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 14638*c[4]*
	c[4]*d[0] - 10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 
	21840*c[3]*c[4]*d[1] + 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*
	c[3]*d[2] - 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] + 
	28600*c[3]*d[3] - 13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*
	c[4]*d[3] - 10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(-21 + 45*c[3] - 26*c[4])
	*d[0] + 11*(21 + 64*c[3] + 10*c[4])*d[1] - 2*(-495 + 77*c[3] - 398*c[4])*
	d[2] - (-495 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] + 63*
	d[1] + 24*d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(-1859 + 819*c[3] - 3478*c[4])
	*d[4] - 2*(3764*c[3]*c[3] + 11*c[3]*(-637 + 934*c[4]) + 2*c[4]*(-7319 + 
	1649*c[4]))*d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 77*d[3]
	 + 398*d[4]) + 13*c[1]* (33*(35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] + 2*
	(22*(42 + 24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(21 + 66*c[2] + 64*c[3] + 
	10*c[4])*d[2] + 2*(-198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (-429 +
	 110*c[2] + 840*c[3] + 1094*c[4])*d[4]))))))/720720.;
		bf_mom[6] = ((a - b)*(2*a*a*(13*(33* (7*(-20 + 5*c[0]*c[0] - 10*
	c[0]*c[1] + 8*c[1]*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*
	(7*c[0] - 8*c[1] - 10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1]
	 + 26*c[2] + 49*c[3])*c[4] + 2252*c[4]*c[4])*d[0] + 1716*(35*d[1] + 14*
	d[2] - 21*d[3] + 10*d[4]) - 429*c[0]*c[0]*(35*d[1] + 14*d[2] - 21*d[3] +
	 10*d[4]) + 52*c[0]*(-11*(c[3]*(36*d[1] + 45*d[2] - 100*d[3]) + c[4]*
	(-39*d[1] + 26*d[2] + 49*d[3])) + (-539*c[3] + 1126*c[4])*d[4] + 33*c[1]*
	(28*d[1] - 7*d[2] - 12*d[3] + 13*d[4]) - 11*c[2]*(21*d[1] - 90*d[2] + 45*
	d[3] + 26*d[4])) + 4*(-6721*c[3]*c[3]*d[1] + 10920*c[3]*c[4]*d[1] - 7111*
	c[4]*c[4]*d[1] - 4264*c[3]*c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 3478*c[4]*
	c[4]*d[2] + 26*c[2]*(c[3]*(352*d[1] + 77*d[2] - 328*d[3]) + c[4]*(-55*d[1]
	 + 398*d[2] - 63*d[3])) + 6903*c[3]*c[3]*d[3] - 7528*c[3]*c[4]*d[3] + 
	5137*c[4]*c[4]*d[3] - 13*c[2]*c[2]*(363*d[1] + 374*d[2] - 77*d[3] - 398*
	d[4]) - 2*c[2]*(819*c[3] + 3478*c[4])*d[4] - 2*(1882*c[3]*c[3] - 5137*
	c[3]*c[4] + 1649*c[4]*c[4])* d[4] - 143*c[1]*c[1]* (63*d[1] - 24*d[2] - 
	9*d[3] + 28*d[4]) + 26*c[1]*(11*c[2]*(24*d[1] - 33*d[2] + 32*d[3]-5*d[4])
	 + c[3]*(99*d[1] + 352*d[2] - 517*d[3] + 420*d[4]) - c[4]*(308*d[1] + 55*
	d[2] - 420*d[3] + 547*d[4])))) + b*(13*(33*(7*(5*(3 + 2*b)*(-4+c[0]*c[0])
	 + 20*(1 + b)*c[0]*c[1] + 4*(5 + 4*b)*c[1]*c[1]) - 28*((5 + 2*b)*c[0] - 
	2*(1 + b)*c[1])*c[2] + 4*(49 + 30*b)*c[2]*c[2]) - 396*(7*(1 + b)*c[0] + 
	2*(7 + 4*b)*c[1] - 10*(1 + b)*c[2])* c[3] + 44*(153 + 100*b)*c[3]*c[3] - 
	44*(3*(7 + 10*b)*c[0] + 2*(39*(1 + b)*c[1] + (57 + 26*b)*c[2] - 49*(1+b)*
	c[3]))*c[4] + 4*(1705 + 1126*b)*c[4]*c[4])*d[0] - 24024*(5*d[1] - 5*d[2] 
	- 3*d[3] - d[4]) + 2*(3003*c[0]*c[0]*(5*d[1] - 5*d[2] - 3*d[3] - d[4]) + 
	572*c[0]*(-63*c[3]*d[1] - 39*c[4]*d[1] + 45*c[3]*d[2] - 57*c[4]*d[2] + 
	153*c[3]*d[3] + 49*c[4]*d[3] + 3*c[2]*(7*d[1] + 49*d[2] + 15*d[3] - 19*
	d[4]) + 3*c[1]*(7*(5*d[1] + d[2] - 3*d[3]) - 13*d[4]) + 49*c[3]*d[4] + 
	155*c[4]*d[4]) + b*(13*(33*(-140 + 35*c[0]*c[0] + 84*c[1]*c[1] + 64*c[1]*
	c[2] + 44*c[2]*c[2] + 28*c[0]*(4*c[1] + c[2])) - 88*(18*c[0] + 9*c[1] - 
	32*c[2])*c[3] + 2068*c[3]*c[3] - 4*(429*c[0] + 616*c[1] - 110*c[2] - 840*
	c[3])*c[4] + 2188*c[4]*c[4])*d[1] - 2*(143*(-84 + 21*c[0]*c[0] - 48*c[1]*
	c[1] - 132*c[1]*c[2] + 68*c[2]*c[2] - 6*c[0]*(7*c[1] + 30*c[2])) - 286*
	(45*c[0] + 64*c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*(143*c[0] - 55*
	c[1] - 398*c[2] - 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (13*(99*(-28 + 
	(c[0] + 2*c[1])*(7*c[0] + 2*c[1])) - 44*(45*c[0] + 64*c[1])*c[2] + 308*
	c[2]*c[2] - 8*(550*c[0] + 517*c[1] - 328*c[2])*c[3] + 2124*c[3]*c[3]) - 
	4*(91*(77*c[0] + 120*c[1] + 18*c[2]) - 7528*c[3])* c[4] + 20548*c[4]*c[4])
	*d[3] - 2*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 286*c[2] - 539*c[3] - 
	1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*c[2] - 819*c[2]*c[3] 
	+ 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] + 1649*c[4]*c[4] - 13*
	c[1]*(55*c[2] + 420*c[3] + 547*c[4])))*d[4]) + 4*(6721*c[3]*c[3]*d[1] + 
	14014*c[3]*c[4]*d[1] + 7111*c[4]*c[4]*d[1] - 8437*c[3]*c[3]*d[2] + 1638*
	c[3]*c[4]*d[2] - 7943*c[4]*c[4]*d[2] - 6903*c[3]*c[3]*d[3] - 10738*c[3]*
	c[4]*d[3] - 5137*c[4]*c[4]*d[3] + 429*c[1]*c[1]*(21*d[1] + 7*d[2] - 3*d[3]
	 - 13*d[4]) - (5369*c[3]*c[3] + 10274*c[3]*c[4] + 2331*c[4]*c[4])*d[4] + 
	143*c[2]*c[2]*(33*d[1] - 81*d[2] - 7*d[3] + 67*d[4]) + 26*c[2]*(c[4]*
	(55*d[1] + 737*d[2] + 63*d[3] - 611*d[4]) + c[3]*(495*d[1] - 77*d[2] - 
	649*d[3] + 63*d[4])) + 26*c[1]*(11*c[2]*(21*d[1] + 33*d[2] + 45*d[3] + 
	5*d[4]) + 11*c[3]*(-9*d[1] + 45*d[2] + 47*d[3] + 49*d[4]) + c[4]*(-429*
	d[1] + 55*d[2] + 539*d[3] + 547*d[4]))))) + a*(13*(33*(7*(5*(3 + 2*b)*(-4
	 + c[0]*c[0]) - 20*c[0]*c[1] + 4*(5 + 2*b)*c[1]*c[1]) - 28*((5 + 6*b)*c[0]
	 + 2*c[1])*c[2] + 4*(49 + 38*b)*c[2]*c[2]) + 396*(7*c[0] - 2*((7 + 6*b)*
	c[1] + 5*c[2]))*c[3] + 44*(153 + 106*b)*c[3]*c[3] + 44*(3*(-7 + 6*b)*c[0]
	 - 2*(-39*c[1] + (57 + 62*b)*c[2] + 49*c[3]))*c[4] + 4*(1705 + 1158*b)*
	c[4]*c[4])*d[0] + 24024*(5*d[1] + 5*d[2] - 3*d[3] + d[4]) + 2*(-429*c[0]*
	c[0]* (35*d[1] + 7*(5 + 6*b)*d[2] - 21*d[3] + (7 - 6*b)*d[4]) + 52*c[0]*
	(-11*(3*c[2]*(7*d[1] - (49 + 38*b)*d[2] + 15*d[3]) + c[4]*(-39*d[1] + (57
	 + 62*b)*d[2] + 49*d[3]) + c[3]*(9*(7 + 6*b)*d[1] + 45*d[2] - (153+106*b)*
	d[3])) - (11*(57 + 62*b)*c[2] + 539*c[3] - (1705 + 1158*b)*c[4])*d[4] + 
	33*c[1]*(7*(5 + 2*b)*d[1] - 7*d[2] - 3*(7 + 6*b)*d[3] + 13*d[4])) + 4* 
	(-6721*c[3]*c[3]*d[1] + 14014*c[3]*c[4]*d[1] + 6188*b*c[3]*c[4]*d[1] - 
	7111*c[4]*c[4]*d[1] + 18018*b*d[2] - 8437*c[3]*c[3]*d[2] - 8346*b*c[3]*
	c[3]*d[2] - 1638*c[3]*c[4]*d[2] - 7943*c[4]*c[4]*d[2] - 8930*b*c[4]*c[4]*
	d[2] + 6903*c[3]*c[3]*d[3] - 10738*c[3]*c[4]*d[3] - 6420*b*c[3]*c[4]*d[3]
	 + 5137*c[4]*c[4]*d[3] + 26*c[2]*(c[4]*(-55*d[1] + (737 + 678*b)*d[2] - 
	63*d[3]) + c[3]*(11*(45 + 26*b)*d[1] + 77*d[2] - (649 + 642*b)*d[3])) - 
	2*c[2]*(819*c[3] + 47*(169 + 190*b)*c[4])*d[4] - (5369*c[3]*c[3] + 6*b*
	(429 + 535*c[3]*c[3]) - 10274*c[3]*c[4] - (-2331 + 1934*b)*c[4]*c[4])*d[4]
	 - 143*c[1]*c[1]* (63*d[1] + 3*(-7 + 2*b)*d[2] - 9*d[3] + (39+22*b)*d[4])
	 - 13*c[2]*c[2]* (363*d[1] + 11*(81 + 94*b)*d[2] - 77*d[3] - (737+678*b)*
	d[4]) - 26*c[1]*(11*c[2]* (3*(-7 + 2*b)*d[1] + 33*d[2] - (45 + 26*b)*d[3]
	 + 5*d[4]) + c[4]*(11*(39 + 22*b)*d[1] + 55*d[2] - 7*(77 + 34*b)*d[3] +
	 547*d[4]) + c[3]*(-99*d[1] - 11*(45 + 26*b)*d[2] + 517*d[3] - 7*
	(77 + 34*b)*d[4])))))))/1441440.;
		bf_mom[7] = ((a - b)*(429*c[0]*c[0]*(35*(-3 + a*a + a*b + b*b)*
	 d[0] - 7*(-10*d[2] + 6*a*b*d[2] + a*a*(5*d[1] + 2*d[2] - 3*d[3]) + b*b*
	(-5*d[1] + 2*d[2] + 3*d[3])) - 2*(-7 + 5*a*a - 3*a*b + 5*b*b)*d[4]) - 
	26*c[0]*(33*(7*(-5*(3 + 2*c[2]) + a*b*(5 + 6*c[2]) + a*a*(5 + 5*c[1] + 
	2*c[2] - 3*c[3]) + b*b*(5 - 5*c[1] + 2*c[2] + 3*c[3])) + 2*(-7 + 5*a*a -
	 3*a*b + 5*b*b)*c[4])*d[0] + 11*(-(b*b*(3* (-35 + 56*c[1] + 14*c[2] - 
	24*c[3] - 26*c[4])*d[1] + 2*(21 + 21*c[1] + 90*c[2] + 45*c[3] - 26*c[4])
	*d[2] + (63 - 72*c[1] + 90*c[2] + 200*c[3] + 98*c[4])*d[3])) + 6*((35 + 
	98*c[2] - 38*c[4])*d[2] + 14*c[1]*(5*d[1] - 3*d[3]) + 6*c[3]*(-7*d[1] + 
	17*d[3]))) + 2*(231 - 1254*c[2] + b*b*(-165 + 429*c[1] + 286*c[2] - 539*
	c[3] - 1126*c[4]) + 3410*c[4])*d[4] + 2*a*b*(-11*((63 + 114*c[2]-62*c[4])
	*d[2] + 6*c[1]*(7*d[1] - 9*d[3]) + c[3]*(-54*d[1] + 106*d[3])) + (99 + 
	682*c[2] - 1158*c[4])*d[4]) + a*a*(-33*(35 + 56*c[1] - 14*c[2] - 24*c[3]
	 + 26*c[4])* d[1] + 22*(-21 + 21*c[1] - 90*c[2] + 45*c[3] + 26*c[4])* d[2]
	 + 11*(63 + 72*c[1] + 90*c[2] - 200*c[3] + 98*c[4])* d[3] - 2*(165 + 429*
	c[1] - 286*c[2] - 539*c[3] + 1126*c[4])* d[4])) + 2*(26*(-11* (3*(35*c[1]*
	c[1] + 7*c[2]*(5 + 7*c[2]) - 42*c[1]*c[3] + 51*c[3]*c[3]) - 3*(-7+38*c[2])
	*c[4] + 155*c[4]*c[4])*d[0] - 22*(3*c[1]*(-35 + 14*c[2] - 26*c[4]) + c[3]*
	(63 + 90*c[2] + 98*c[4]))*d[1] - 2*(231*c[1]*c[1] - 33*c[2]*(49 + 27*c[2])
	 + 990*c[1]*c[3] - 649*c[3]*c[3] + 11*(57 + 134*c[2])*c[4]-611*c[4]*c[4])*
	d[2] - 2*(11*c[1]*(63 + 90*c[2] + 98*c[4]) - c[3]*(1683 + 1298*c[2] + 826*
	c[4]))*d[3]) + 4*(13*(429*c[1]*c[1] - 11*c[2]*(57 + 67*c[2]) - 1078*c[1]*
	c[3] + 413*c[3]*c[3]) + 13*(1705 + 1222*c[2])*c[4] + 2331*c[4]*c[4])*d[4]
	 + 2*a*b*(7579*c[3]*c[3]*d[0] - 1287*c[4]*d[0] + 7527*c[4]*c[4]*d[0] + 
	7722*c[3]*d[1] + 6188*c[3]*c[4]*d[1] - 8346*c[3]*c[3]*d[2] + 8866*c[4]*
	d[2] - 8930*c[4]*c[4]*d[2] - 15158*c[3]*d[3] - 6420*c[3]*c[4]*d[3] - 26*
	c[1]*(11*(21 + 6*c[2] + 22*c[4])*d[1] - (297 + 286*c[2] + 238*c[4])*d[3]
	 + c[3]*(297*d[0] - 286*d[2] - 238*d[4])) + 143*c[1]*c[1]*(21*d[0] - 6*
	d[2] - 22*d[4]) - 2*(1605*c[3]*c[3] + (7527 - 967*c[4])*c[4])*d[4] + 13*
	c[2]*c[2]*(627*d[0] - 1034*d[2] + 678*d[4]) + c[2]*(-143*(-63 + 62*c[4])*
	d[0] + 78*(-209 + 226*c[4])*d[2] + 52*c[3]*(143*d[1] - 321*d[3]) + 2*
	(4433 - 8930*c[4])*d[4])) + a*a*(-9009*c[3]*d[0] + 14300*c[3]*c[3]*d[0] 
	+ 4290*c[4]*d[0] - 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 10296*
	c[3]*d[1] - 13442*c[3]*c[3]*d[1] - 11154*c[4]*d[1] + 21840*c[3]*c[4]*d[1]
	 - 14222*c[4]*c[4]*d[1] + 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] + 7436*
	c[4]*d[2] - 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*d[3] + 
	13806*c[3]*c[3]*d[3] + 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] + 10274*
	c[4]*c[4]*d[3] - 26*c[2]*(11*(-21 + 45*c[3] + 26*c[4])*d[0] - 11*(21 + 
	64*c[3] - 10*c[4])*d[1] - 2*(-495 + 77*c[3] + 398*c[4])*d[2] + (-495 + 
	656*c[3] + 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] - 63*d[1] + 24*d[2] 
	+ 9*d[3] - 28*d[4]) - 4*c[2]*(-1859 + 819*c[3] + 3478*c[4])*d[4] - 2*
	(3764*c[3]*c[3] - 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*
	d[4] + 26*c[2]*c[2]*(495*d[0] - 363*d[1] - 374*d[2] + 77*d[3] + 398*d[4])
	 - 13*c[1]* (33*(-35 + 14*c[2] + 24*c[3] - 26*c[4])*d[0] + 2*(-22*(-42 + 
	24*c[2] + 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] - 64*c[3] + 10*c[4])
	*d[2] - 2*(198 + 352*c[2] - 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2] 
	- 840*c[3] + 1094*c[4])*d[4]))) + b*b*(9009*c[3]*d[0] + 14300*c[3]*c[3]*
	d[0] + 4290*c[4]*d[0] + 14014*c[3]*c[4]*d[0] + 14638*c[4]*c[4]*d[0] + 
	10296*c[3]*d[1] + 13442*c[3]*c[3]*d[1] + 11154*c[4]*d[1] + 21840*c[3]*
	c[4]*d[1] + 14222*c[4]*c[4]*d[1] - 12870*c[3]*d[2] - 8528*c[3]*c[3]*d[2] 
	+ 7436*c[4]*d[2] + 3276*c[3]*c[4]*d[2] - 6956*c[4]*c[4]*d[2] - 28600*c[3]*
	d[3] - 13806*c[3]*c[3]*d[3] - 14014*c[4]*d[3] - 15056*c[3]*c[4]*d[3] - 
	10274*c[4]*c[4]*d[3] + 26*c[2]*(11*(21 + 45*c[3] - 26*c[4])*d[0] + 11*
	(-21 + 64*c[3] + 10*c[4])*d[1] - 2*(495 + 77*c[3] - 398*c[4])*d[2] - (495
	 + 656*c[3] - 126*c[4])*d[3]) + 286*c[1]*c[1]* (42*d[0] + 63*d[1] + 24*
	d[2] - 9*d[3] - 28*d[4]) + 4*c[2]*(1859 + 819*c[3] - 3478*c[4])*d[4] - 2*
	(3764*c[3]*c[3] + 11*c[3]*(637 + 934*c[4]) + 2*c[4]*(7319 + 1649*c[4]))*
	d[4] + 26*c[2]*c[2]*(495*d[0] + 363*d[1] - 374*d[2] - 77*d[3] + 398*d[4])
	 + 13*c[1]* (33*(-35 + 14*c[2] - 24*c[3] - 26*c[4])*d[0] + 2*(22*(-42 + 
	24*c[2] - 9*c[3] - 28*c[4])*d[1] + 11*(-21 + 66*c[2] + 64*c[3] + 10*c[4])
	*d[2] + 2*(198 + 352*c[2] + 517*c[3] + 420*c[4])*d[3] + (429 + 110*c[2]
	 + 840*c[3] + 1094*c[4])*d[4]))))))/720720.;
		bf_mom[8] = -((a - b)*(13*(33*(35*(-3 + b*b)*(-4 + c[0]*c[0]) + 
	70*b*b*c[0]*c[1] + 28*(-5 + 2*b*b)*c[1]*c[1] - 28*((-5 + b*b)*c[0] - b*b*
	c[1])*c[2] + 4*(-49 + 15*b*b)*c[2]*c[2]) - 198*(-28*c[1] + b*b*(7*c[0] + 
	8*c[1] - 10*c[2]))*c[3] + 44*(-153 + 50*b*b)*c[3]*c[3] - 44*(3*(-7+5*b*b)*
	c[0] - 114*c[2] + b*b*(39*c[1] + 26*c[2] - 49*c[3]))*c[4] + 4*(-1705 + 
	563*b*b)*c[4]*c[4] + a*a*(33*(7*(-20 + 5*c[0]*c[0] - 10*c[0]*c[1] + 8*c[1]
	*c[1]) - 28*(c[0] + c[1])*c[2] + 60*c[2]*c[2]) + 198*(7*c[0] - 8*c[1] - 
	10*c[2])*c[3] + 2200*c[3]*c[3] - 44*(15*c[0] - 39*c[1] + 26*c[2]+49*c[3])
	*c[4] + 2252*c[4]*c[4]) + a*b*(33*(35*c[0]*c[0] - 84*c[0]*c[2] + 76*c[2]*
	c[2]) + 44*(21*(-5 + c[1]*c[1]) - 54*c[1]*c[3] + 53*c[3]*c[3]) + 44*
	(9*c[0] - 62*c[2])*c[4] + 2316*c[4]*c[4]))*d[0] - 60060*b*b*d[1] + 15015*
	b*b*c[0]*c[0]*d[1] - 120120*c[0]*c[1]*d[1] + 48048*b*b*c[0]*c[1]*d[1] + 
	36036*b*b*c[1]*c[1]*d[1] + 12012*b*b*c[0]*c[2]*d[1] - 48048*c[1]*c[2]*
	d[1] + 27456*b*b*c[1]*c[2]*d[1] + 18876*b*b*c[2]*c[2]*d[1] + 72072*c[0]*
	c[3]*d[1] - 20592*b*b*c[0]*c[3]*d[1] - 10296*b*b*c[1]*c[3]*d[1] - 102960*
	c[2]*c[3]*d[1] + 36608*b*b*c[2]*c[3]*d[1] + 26884*b*b*c[3]*c[3]*d[1] - 
	22308*b*b*c[0]*c[4]*d[1] + 89232*c[1]*c[4]*d[1] - 32032*b*b*c[1]*c[4]*d[1]
	 + 5720*b*b*c[2]*c[4]*d[1] - 112112*c[3]*c[4]*d[1] + 43680*b*b*c[3]*c[4]*
	d[1] + 28444*b*b*c[4]*c[4]*d[1] - 120120*d[2] + 24024*b*b*d[2] + 30030*
	c[0]*c[0]*d[2] - 6006*b*b*c[0]*c[0]*d[2] + 12012*b*b*c[0]*c[1]*d[2] - 
	24024*c[1]*c[1]*d[2] + 13728*b*b*c[1]*c[1]*d[2] - 168168*c[0]*c[2]*d[2] 
	+ 51480*b*b*c[0]*c[2]*d[2] + 37752*b*b*c[1]*c[2]*d[2] + 92664*c[2]*c[2]*
	d[2] - 19448*b*b*c[2]*c[2]*d[2] + 25740*b*b*c[0]*c[3]*d[2] - 102960*c[1]*
	c[3]*d[2] + 36608*b*b*c[1]*c[3]*d[2] - 8008*b*b*c[2]*c[3]*d[2] + 67496*
	c[3]*c[3]*d[2] - 17056*b*b*c[3]*c[3]*d[2] + 65208*c[0]*c[4]*d[2] - 14872*
	b*b*c[0]*c[4]*d[2] + 5720*b*b*c[1]*c[4]*d[2] - 153296*c[2]*c[4]*d[2] + 
	41392*b*b*c[2]*c[4]*d[2] + 6552*b*b*c[3]*c[4]*d[2] + 63544*c[4]*c[4]*d[2]
	 - 13912*b*b*c[4]*c[4]*d[2] + 36036*b*b*d[3] - 9009*b*b*c[0]*c[0]*d[3] + 
	72072*c[0]*c[1]*d[3] - 20592*b*b*c[0]*c[1]*d[3] - 5148*b*b*c[1]*c[1]*d[3]
	 + 25740*b*b*c[0]*c[2]*d[3] - 102960*c[1]*c[2]*d[3] + 36608*b*b*c[1]*c[2]*
	d[3] - 4004*b*b*c[2]*c[2]*d[3] - 175032*c[0]*c[3]*d[3] + 57200*b*b*c[0]*
	c[3]*d[3] + 53768*b*b*c[1]*c[3]*d[3] + 134992*c[2]*c[3]*d[3] - 34112*b*b*
	c[2]*c[3]*d[3] - 27612*b*b*c[3]*c[3]*d[3] + 28028*b*b*c[0]*c[4]*d[3] - 
	112112*c[1]*c[4]*d[3] + 43680*b*b*c[1]*c[4]*d[3] + 6552*b*b*c[2]*c[4]*
	d[3] + 85904*c[3]*c[4]*d[3] - 30112*b*b*c[3]*c[4]*d[3] - 20548*b*b*c[4]*
	c[4]*d[3] - 2*(-3003*(-4 + c[0]*c[0]) - 52*(429*c[1]*c[1] + 627*c[0]*c[2]
	 - 737*c[2]*c[2] - 1078*c[1]*c[3] + 413*c[3]*c[3]) + 52*(1705*c[0] - 1222*
	c[2])*c[4] - 9324*c[4]*c[4] + b*b*(2145*c[0]*c[0] + 26*c[0]*(429*c[1] + 
	286*c[2] - 539*c[3] - 1126*c[4]) + 4*(-2145 + 2002*c[1]*c[1] - 2587*c[2]*
	c[2] - 819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*c[2]*c[4] + 5137*c[3]*c[4] 
	+ 1649*c[4]*c[4] - 13*c[1]*(55*c[2] + 420*c[3] + 547*c[4]))))*d[4] - a*a*
	(13*(1155*c[0]*c[0] + 44*(-105 + 63*c[1]*c[1] + 33*c[2]*c[2] - 64*c[2]*
	c[3] + 47*c[3]*c[3] - 6*c[1]*(8*c[2] + 3*c[3])) + 8*(308*c[1] + 55*c[2] 
	- 420*c[3])*c[4] + 2188*c[4]*c[4] - 132*c[0]*(28*c[1] - 7*c[2] - 12*c[3] 
	+ 13*c[4]))*d[1] + 2*(143*(-84 + 21*c[0]*c[0] - 48*c[1]*c[1] + 6*c[0]*
	(7*c[1] - 30*c[2]) + 132*c[1]*c[2] + 68*c[2]*c[2]) + 286*(45*c[0] - 64*
	c[1] - 14*c[2])*c[3] + 8528*c[3]*c[3] + 52*(143*c[0] + 55*c[1] - 398*c[2]
	 + 63*c[3])*c[4] + 6956*c[4]*c[4])*d[2] - (9009*c[0]*c[0] + 52*(11*(-63 +
	 (9*c[1] + c[2])*(c[1] + 7*c[2])) - 2*(517*c[1] + 328*c[2])*c[3] + 531*
	c[3]*c[3]) + 8*(5460*c[1] - 819*c[2] - 3764*c[3])*c[4] + 20548*c[4]*c[4] 
	- 572*c[0]*(36*c[1] + 45*c[2] - 100*c[3] + 49*c[4]))*d[3] + 2*(2145*c[0]*
	c[0] - 26*c[0]*(429*c[1] - 286*c[2] - 539*c[3] + 1126*c[4]) + 4*(-2145 + 
	2002*c[1]*c[1] - 2587*c[2]*c[2] + 819*c[2]*c[3] + 1882*c[3]*c[3] + 3478*
	c[2]*c[4] - 5137*c[3]*c[4] + 1649*c[4]*c[4] + 13*c[1]*(55*c[2] - 420*c[3]
	 + 547*c[4])))*d[4]) + 2*a*b*(-104*(11*c[1]*(3*c[2] + 11*c[4]) - c[3]*
	(143*c[2] + 119*c[4]))*d[1] - 4*(429*c[1]*c[1] + 6721*c[2]*c[2] - 3718*
	c[1]*c[3] + 4173*c[3]*c[3] - 8814*c[2]*c[4] + 4465*c[4]*c[4])*d[2] + 8*
	(-321*c[3]*(13*c[2] + 5*c[4]) + 13*c[1]*(143*c[2] + 119*c[4]))*d[3] + 
	5148*(7*d[2] - d[4]) - 1287*c[0]*c[0]*(7*d[2] - d[4]) - 4*(1573*c[1]*c[1]
	 - 4407*c[2]*c[2] - 3094*c[1]*c[3] + 1605*c[3]*c[3] + 8930*c[2]*c[4] - 
	967*c[4]*c[4])*d[4] + 52*c[0]*(627*c[2]*d[2] - 341*c[4]*d[2] + 33*c[1]*
	(7*d[1] - 9*d[3]) + c[3]*(-297*d[1] + 583*d[3]) - 341*c[2]*d[4] 
	+ 579*c[4]*d[4]))))/ 360360.;
			break;
	default:
			printf("unknown F_type %d\n",interface_type);
	}  /* end F-type switch */
	break;
default:
	EH(-1,"unknown Chebyshev polynomial order");
	break;
}

return;
}

/*  end of delta_chev_moments_2DQ		*/

#endif
