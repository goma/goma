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
 * Routines for calculating Boundary Conditions and adding them
 * to matrix_fill
 *
 *
 *$Id: bc_colloc.c,v 5.11 2010-07-21 16:39:26 hkmoffa Exp $
 */

/* Standard include files */
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "dpi.h"
 
#define _BC_COLLOC_C
#include "goma.h"



/*******************************************************************************/

int
apply_point_colloc_bc (
     double resid_vector[],   /* Residual vector for the current processor     */
     double delta_t,          /* current time step size                        */
     double theta,            /* parameter to vary time integration:
                               * explicit (theta = 1) -- implicit (theta = 0)  */
     int ielem,               /* element number                                */
     int ip_total,            /* total number of gauss points                  */
     int ielem_type,          /* element type                                  */
     int num_local_nodes,
     int ielem_dim,
     int iconnect_ptr,
     ELEM_SIDE_BC_STRUCT *elem_side_bc,
                              /* Pointer to an element side boundary
                               * condition structure                           */
     int num_total_nodes,
     int local_node_list_fs[],/* dimensioned [MDE]; list to keep track of
                               *  nodes at which solid contributions have been
                               *  transfered to liquid (fluid-solid
			       *  boundaries)                                  */
     double time_value,
     Exo_DB *exo)

/*******************************************************************************
  Function which fills the FEM stiffness matrices and the right-hand-side vector
  with contributions from boundary conditions.
 
  Authors:         Harry Moffat, Scott Hutchinson (1421) and others
                   Transferred to separate files by Rich Cairncross
  Date:            12 December 1994
******************************************************************************/
{
  int w, i, I, ibc, j, id, err, var, eqn, ldof_eqn, ldof_var, icount;
  int ieqn, pvar;
  int el1, el2, nf, jk;
  int status = 0;
  int index_eq, matID_apply;
  int bc_input_id;
  int var_flag;
  int dof_map[MDE];
  int n_dof[MAX_VARIABLE_TYPES];
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  int doMeshMapping = 0;
  double xi[DIM];             /* Gaussian-quadrature point locations         */
  double x_dot[MAX_PDIM];
  double dsigma_dx[DIM][MDE];
  double phi_j;
  double func;
  double f_time;
  double d_func[MAX_VARIABLE_TYPES + MAX_CONC];
  double kfunc[DIM];
  double d_kfunc[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double time_intermediate = time_value-theta*delta_t;
                   /* time at which bc's are evaluated */
  const double penalty = BIG_PENALTY;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int doFullJac = 0;
  double nwall[3];
  int contact_flag = FALSE;
  double xsurf[MAX_PDIM];


  /***************************************************************************/
  /*    START OF SURFACE LOOPS THAT DON'T REQUIRE INTEGRATION (STRONG SENSE) */
  /*      (POINTWISE COLLOcATION)                                            */
  /*              LOOP OVER THE SURFACES IN THE CURRENT ELEMENT THAT	     */
  /*         NEED THE APPLICATION OF A STONG SENSE CONDITION                 */
  /*     - INITIALIZATION THAT IS DEPENDENT ON THE IDENTITY OF THE SURFACE   */
  /***************************************************************************/
    
  /***************************************************************************/
  /*             LOOP OVER THE LOCAL ELEMENT NODE NUMBER                     */
  /*          FOR NODES THAT LIE ON THE CURRENT SURFACE		             */
  /*   - INITIALIZATION THAT IS DEPENDENT ON THE LOCAL ELEMENT NODE NUMBER   */
  /***************************************************************************/
#ifdef DEBUG_HKM
  if (ei->ielem == 1032) {
  //  printf("we are here - 1032\n");
  }
#endif
 
  for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) {
    
    /* Find the local element node number for the current node */
    id = (int) elem_side_bc->local_elem_node_id[i];
    
    /* Find the local node number given the local element node 
     * number,  'id' 
     * I is the global node number of this node !!  
     */
    I = Proc_Elem_Connect[iconnect_ptr + id];
    
    /*
     * Check to see if the current local node number is owned by the 
     *  processor - this is done by checking I < num_owned_nodes
     */
    if (I < DPI_ptr->num_owned_nodes) {
      
      /*
       * Load the field variable values at this node point
       */
      find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);

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
      
      /* calculate the shape functions and their gradients */

      /* calculate the determinant of the surface jacobian  and the normal to 
       * the surface all at one time */
      surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes, 
				     ielem_dim - 1,  
				     (int) elem_side_bc->id_side,
				     (int) elem_side_bc->num_nodes_on_side,
				     (elem_side_bc->local_elem_node_id) );

      if (ielem_dim !=3) {
	calc_surf_tangent(ielem, iconnect_ptr, num_local_nodes, ielem_dim-1,
			  (int) elem_side_bc->num_nodes_on_side,
			  (elem_side_bc->local_elem_node_id));      
      }

      /*
       * Load up porous media variables and properties, if needed 
       */
      if (mp->PorousMediaType == POROUS_UNSATURATED || 
	  mp->PorousMediaType == POROUS_TWO_PHASE   || 
	  mp->PorousMediaType == POROUS_SATURATED) {
	err = load_porous_properties(); 
	EH(err, "load_porous_properties"); 
      }

      if (TimeIntegration != STEADY) {
	if (mp->SurfaceTensionModel != CONSTANT)  {
	  load_surface_tension(dsigma_dx);
	  if (neg_elem_volume) return (status);
	}
      }
    
      if (TimeIntegration != STEADY) {
	for (icount = 0; icount < ielem_dim; icount++ ) {
	  x_dot[icount] = fv_dot->x[icount];
	}
      } else {
	for (icount = 0; icount < ielem_dim; icount++ ) {
	  x_dot[icount] = 0.0;
	}
      }
    
      do_LSA_mods(LSA_SURFACE);

      /**********************************************************************/
      /*                                                                    */
      /*         LOOP OVER THE INDIVIDUAL SURFACE CONDITION TERMS           */
      /*		Calculate the functions for each collocation        */
      /*                condition at this point                             */
      /*		                         		            */
      /**********************************************************************/
      
      for (ibc = 0; 
	   (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1 ;
	   ibc++) {
	/*
	 * This function only handles COLLOCATE_SURF boundary conditions.
	 * All other boundary conditions are weeded out here
	 */
	if (BC_Types[bc_input_id].desc->method == COLLOCATE_SURF) {
	  /* Initialize time function to unity */
	  f_time = 1.0;

	  /* 
	   * Check to see if this BC on this node is applicable
	   * (i.e. no other overriding Dirichlet conditions,
	   * And find the global unknown number for applying this condition
	   */
	  index_eq = bc_eqn_index(id, I, bc_input_id, ei->mn,
				  0, &eqn, &matID_apply, &vd);
	  if (index_eq >= 0) {

	    /* initialize the general function to zero */
	    func = 0.;
	    init_vec_value(d_func, 0.0, MAX_VARIABLE_TYPES + MAX_CONC);
	    doFullJac = 0;
	    doMeshMapping = 0;
/*
 * Here's a RECIPE for adding new boundary conditions so you don't have any
 * excuses not to add new ones.  The changes should be made in at least
 * four files (rf_bc_const.h, mm_names.h, mm_input.c, and bc_[method].c)
 * for some boundary conditions you may want to make additions elsewhere also.
 * One example of extra additions is in el_exoII_io.c where the ss_dup_list
 * is created - you may want to adapt the logic for rotation of new conditions.
 * (note that these lines are repeated at each place where you need to 
 * make changes):
 *  Boundary Condition  Step 1: add Macro Substitution for your new boundary
 *                              condition in rf_bc_const.h - this is how you
 *      rf_bc_const.h           will refer to the new boundary condition 
 *                              throughout the code.  Make sure the integer
 *                              you choose is unique from all the other BC
 *                              types.
 *  Boundary Condition  Step 2: add a description of your boundary condition
 *                              to the BC_Desc structure array in mm_names.h.
 *      mm_names.h              This structure includes information about the 
 *                              type of boundary condition, which equation it
 *                              applies to, what variables it is sensitive to,
 *                              whether to rotate mesh or momentum equations, 
 *                              etc.  It is very important that you fill out
 *                              this structure carefully, otherwise the code
 *                              won't know what to do.
 *  Boundary Condition  Step 3: add your BC case to the correct format listing 
 *                              for reading the necessary arguments from the
 *      mm_input_bc.c           input file in mm_input_bc.c.
 *
 *  Boundary Condition  Step 4: Add a function call (and a function) in the 
 *                              correct routine for evaluating your boundary
 *      bc_colloc.c             condition.  This will probably in bc_colloc.c
 *      bc_integ.c              for collocated conditions or bc_integ.c for 
 *                              strong or weak integrated conditions.
 *  Boundary Condition  Step 5: use and enjoy your new boundary condition
 *
 * Step 4 should be done below (or in bc_integ.c), and add your new function at 
 *     the end of this file (see fplane)
 */
		/* load in functional representation of this condition, and its sensitivity
		 * to all nodal unknowns */
	    switch (BC_Types[bc_input_id].BC_Name) {

	      /* Section for surface terms involving the mesh equations  */
	    case DXDISTNG_BC:
	    case DYDISTNG_BC:
	    case DZDISTNG_BC:
	    case DISTNG_BC:
		fTmelting (&func, d_func,
			   BC_Types[bc_input_id].BC_Data_Float[0]);
		break;
		    
	    case SPLINEX_BC:
	    case SPLINEY_BC:
	    case SPLINEZ_BC:
	    case SPLINE_BC:
		fspline(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
			time_intermediate);
		break;

 	    case SPLINEX_RS_BC:
 	    case SPLINEY_RS_BC:
 	    case SPLINEZ_RS_BC:
 	    case SPLINE_RS_BC:
 		fspline_rs(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
 			time_intermediate);
 		break;

            case UUSER_COLLOC_BC:
                uuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC,
                                  id, time_intermediate);
                break;
            case VUSER_COLLOC_BC:
               	vuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC,
                                  id, time_intermediate);
                break;
            case WUSER_COLLOC_BC:
               	wuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC,
                                  id, time_intermediate);
                break;

	    case T_USER_BC:
		tuser(&func, d_func, BC_Types[bc_input_id].u_BC,
		       time_intermediate);
		break;

	    case DX_USER_BC:
	      dx_user_surf(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	      break;

	    case DY_USER_BC:
	      dy_user_surf(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	      break;

	    case DZ_USER_BC:
	      dz_user_surf(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	      break;

	    case P_LIQ_USER_BC:
	      p_liq_user_surf(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	      break;

	    case SH_P_OPEN_USER_BC:
	      shell_p_open_user_surf(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	      break;

	    case POROUS_LIQ_PRESSURE_FILL_BC:
	      {
		MATRL_PROP_STRUCT *mp_2=NULL;

		porous_liq_fill(&func, d_func, BC_Types[bc_input_id].BC_Data_Int[0], 
				BC_Types[bc_input_id].BC_Data_Int[1],
				BC_Types[bc_input_id].BC_Data_Float[0], 
				BC_Types[bc_input_id].BC_Data_Float[1],
				BC_Types[bc_input_id].BC_Data_Float[2],
				mp_2);
	      }
		break;
		    
	    case PLANEX_BC:
	    case PLANEY_BC:
	    case PLANEZ_BC:
	    case PLANE_BC:
		fplane(ielem_dim, &func, d_func, 
		       BC_Types[bc_input_id].BC_Data_Float);
		break;


	    case FILLET_BC:
		f_fillet(ielem_dim, &func, d_func, 
		       BC_Types[bc_input_id].u_BC,BC_Types[bc_input_id].len_u_BC);
		break;

	    case DOUBLE_RAD_BC:
		f_double_rad(ielem_dim, &func, d_func, 
		       BC_Types[bc_input_id].u_BC,BC_Types[bc_input_id].len_u_BC);
		break;

	    case ROLL_FLUID_BC:
                icount = BC_Types[bc_input_id].BC_Data_Int[2];
xsurf[0] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt+1];
xsurf[1] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt+2];
xsurf[2] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt+3];
		f_roll_fluid(ielem_dim, &func, d_func, 
		       BC_Types[bc_input_id].u_BC,BC_Types[bc_input_id].len_u_BC, xsurf);
		break;
	    case MOVING_PLANE_BC:
	    {
	      double t = time_intermediate;

	      fplane (ielem_dim, &func, d_func,
		      BC_Types[bc_input_id].BC_Data_Float);

	      func +=  (BC_Types[bc_input_id].BC_Data_Float[4]*t +
			BC_Types[bc_input_id].BC_Data_Float[5]*t*t +
			BC_Types[bc_input_id].BC_Data_Float[6]*t*t*t);
	      if (af->Assemble_LSA_Mass_Matrix)
		  EH(-1, "LSA is not currently compatible with MOVING_PLANE_BC");
	    }
	    break;

	    case MOVING_PLANE_ETCH_BC:
	      fmesh_etch_bc (&func, d_func,
		      BC_Types[bc_input_id].BC_Data_Int[0], id, x_dot, theta, delta_t);
	    break;

	    case MESH_CONSTRAINT_BC:
	      fmesh_constraint(&func, d_func, bc_input_id);
	      break;

	    case UVARY_BC:
	    case VVARY_BC:
	    case WVARY_BC:
		var_flag = BC_Types[bc_input_id].desc->equation;
		fvelocity_profile(var_flag, ielem_dim, BC_Types[bc_input_id].BC_Name,
				  &func, d_func, BC_Types[bc_input_id].u_BC,
				  time_intermediate);
		break;
		    
	    case U_PARABOLA_BC:
	    case V_PARABOLA_BC:
	    case W_PARABOLA_BC:
		var_flag = BC_Types[bc_input_id].desc->equation;
		fvelocity_parabola(var_flag, ielem_dim, 
			BC_Types[bc_input_id].BC_Name,
			&func, d_func, BC_Types[bc_input_id].u_BC,
			time_intermediate,BC_Types[bc_input_id].len_u_BC);
		break;
		    
	    case GD_CONST_BC:
	    case GD_LINEAR_BC:
	    case GD_INVERSE_BC:
	    case GD_PARAB_BC:
	    case GD_PARAB_OFFSET_BC:
	    case GD_CIRC_BC:
	    case GD_POLYN_BC:
	    case GD_TABLE_BC:
		err = fgeneralized_dirichlet(&func, d_func, BC_Types[bc_input_id].BC_Name, 
					     bc_input_id, theta, delta_t);
		EH(err, "Illegal entry in Generalized Dirichlet Condition ");
		break;

	    case GD_TIME_BC:
		err = evaluate_time_func(time_intermediate, &f_time, bc_input_id);
		EH(err, "Problems in evaluating time function");
		if (af->Assemble_LSA_Mass_Matrix)
		    EH(-1, "LSA is not currently compatible with GD_TIME_BC");
		break;
                 
	    case POROUS_PRESSURE_BC:
		porous_pressure(&func, d_func, 
				(int) BC_Types[bc_input_id].BC_Data_Int[0],
				(int) BC_Types[bc_input_id].BC_Data_Int[1]);
		break;
	    case POROUS_PRESSURE_LUB_BC:
	      porous_pressure_lub(&func, d_func,
				   elem_side_bc->id_side, xi, exo, 
				   BC_Types[bc_input_id].BC_Data_Float[0]);
	      break;

	    case LUBP_SH_FP_FLUX_BC:
	      put_lub_flux_in_film(id, I, ielem_dim, resid_vector,
				   (int) BC_Types[bc_input_id].BC_Data_Int[0],
				   (int) BC_Types[bc_input_id].BC_Data_Int[1],
				   local_node_list_fs);
	      func = 0.;
	      break;

	    case FLUID_SOLID_BC:
	    case SOLID_FLUID_BC:
		put_liquid_stress_in_solid(id, I, 
					   ielem_dim, resid_vector,
					   (int) BC_Types[bc_input_id].BC_Data_Int[0],
					   (int) BC_Types[bc_input_id].BC_Data_Int[1],
					   local_node_list_fs,
					   BC_Types[bc_input_id].BC_Data_Float[0]);
		func = 0.; /* this boundary condition rearranges values already in res and jac,
				  and does not add anything into the residual */
		break;

	    case SOLID_FLUID_RS_BC:
		put_liquid_stress_in_solid_ALE(id, I, 
					       ielem_dim, resid_vector,
					       (int) BC_Types[bc_input_id].BC_Data_Int[0],
					       (int) BC_Types[bc_input_id].BC_Data_Int[1],
					       local_node_list_fs,
					       BC_Types[bc_input_id].BC_Data_Float[0]);
		func = 0.; /* this boundary condition rearranges values already in res and jac,
			    * and does not add anything into the residual */
		break;

	    case FLUID_SOLID_RS_BC:
		EH(-1,"FLUID_SOLID_RS bc not implemented yet");
		break;

	    case SH_FLUID_STRESS_BC:

	      /*Note that we send in i, with id, as this is the shell-element counterpart local num */
	      put_fluid_stress_on_shell(id, i , I, 
					ielem_dim, resid_vector,
					local_node_list_fs,
					BC_Types[bc_input_id].BC_Data_Float[0]);
	      func = 0.; /* this boundary condition rearranges values already in res and jac,
			    * and does not add anything into the residual */
	      break;

	    case KINEMATIC_COLLOC_BC:
	    case VELO_NORM_COLLOC_BC:
	      /* initialize the general function to zero may have more than 
	       * one entry for vector conditions like capillary */
	      memset(kfunc, 0, DIM*sizeof(double) );
	      memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	      fvelo_normal_bc(kfunc, d_kfunc, 
			      BC_Types[bc_input_id].BC_Data_Float[0], contact_flag = FALSE,
			      x_dot, theta, delta_t,
			      (int) BC_Types[bc_input_id].BC_Name,0,0, 135.0);
	      doFullJac = 1;
	      func = kfunc[0];
	      break;
	      
	    case VELO_NORMAL_LS_COLLOC_BC:
	      /* initialize the general function to zero may have more than 
	       * one entry for vector conditions like capillary */
	      memset(kfunc, 0, DIM*sizeof(double) );
	      memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
              contact_flag = (ls != NULL);
	      fvelo_normal_bc(kfunc, d_kfunc, 
			      BC_Types[bc_input_id].BC_Data_Float[0], contact_flag,
			      x_dot, theta, delta_t,
			      (int) BC_Types[bc_input_id].BC_Name,
                              BC_Types[bc_input_id].BC_Data_Float[1],
                              BC_Types[bc_input_id].BC_Data_Float[2],
                              BC_Types[bc_input_id].BC_Data_Float[3]);
	      doFullJac = 1;
	      func = kfunc[0];
	      break;
	      
	    case KIN_DISPLACEMENT_COLLOC_BC:
	      memset(kfunc, 0, DIM*sizeof(double) );
	      memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	      f_kinematic_displacement_bc(kfunc, d_kfunc, 
					  BC_Types[bc_input_id].BC_Data_Int[0],
					  BC_Types[bc_input_id].BC_ID,
                                          BC_Types[bc_input_id].u_BC, 
                                          BC_Types[bc_input_id].len_u_BC);
	      func = kfunc[0];
	      doFullJac = 1;
	      break;
	    case TABLE_BC:
	      apply_table_bc(&func, d_func, &BC_Types[bc_input_id], time_value );
	      if(neg_elem_volume) return(status);
	      break;

	    case SH_GAMMA1_DERIV_SYMM_BC:
	      memset(kfunc, 0, DIM*sizeof(double) );
	      memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	      fgamma1_deriv_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0]);
	      func = kfunc[0];
	      doFullJac = 1;
	      el1 = ei->ielem;
	      nf = num_elem_friends[el1];
	      if (nf == 0) {
		EH(-1, "no friends");
	      }; 
	      el2 = elem_friends[el1][0];
	      err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, -1, xi, exo); 
	      doMeshMapping = 1;
	      break;

	    case SH_GAMMA2_DERIV_SYMM_BC:
	      memset(kfunc, 0, DIM*sizeof(double) );
	      memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	      fgamma2_deriv_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0]);
	      func = kfunc[0];
	      doFullJac = 1;
	      el1 = ei->ielem;
	      nf = num_elem_friends[el1];
	      if (nf == 0) {
		EH(-1, "no friends");
	      }; 	 
	      el2 = elem_friends[el1][0];
	      err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, -1, xi, exo); 
	      doMeshMapping = 1;
	      
	      break;

            case DVZDR_ZERO_BC:
              memset(kfunc, 0, DIM*sizeof(double) );
              memset(d_kfunc,0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
          
              nwall[0] = BC_Types[bc_input_id].BC_Data_Float[1];
              nwall[1] = BC_Types[bc_input_id].BC_Data_Float[2];
              nwall[2] = BC_Types[bc_input_id].BC_Data_Float[3];
              dvzdr_zero_deriv_bc(kfunc, d_kfunc, nwall, BC_Types[bc_input_id].BC_Data_Float[0]);
              func = kfunc[0];
              doFullJac = 1;
              break;

	    } /* end of SWITCH statement */

	    /************************************************************************/
	    /*                                                                      */
	    /*         ADD the Boundary condition functions into the Residual       */
	    /*         vector and Jacobian Matrix                                   */
	    /*		                         			            */
	    /************************************************************************/
	    /*
	     * Collocated boundary conditions are always applied on the first 
	     * dof at a node. They are not discontinuous variables friendly
	     */
	    ldof_eqn = ei->ln_to_first_dof[eqn][id];


	    if (eqn == R_MASS) {
	      ieqn = MAX_PROB_EQN + BC_Types[bc_input_id].species_eq;
	    } else {
	      ieqn = upd->ep[eqn];
	    }

	    if (ldof_eqn != -1)   {
	      lec->R[ieqn][ldof_eqn] += penalty * func;
	      lec->R[ieqn][ldof_eqn] *= f_time;

	      /* 
	       * add sensitivities into matrix
	       *  - find index of sensitivity in matrix
	       *     (if variable is not defined at this node,
	       *      loop over all dofs in element)
	       *  - add into matrix
	       */
	      if (af->Assemble_Jacobian) {

	

		for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
		  pvar = upd->vp[var];
		  if (pvar != -1 && (BC_Types[bc_input_id].desc->sens[var] || 1)) {
		    /*
		     * Warning!!!!!!!!!!!!!!!!!!!!!!!
		     * This section of code will need to be revamped
		     * to work for discontinuous variable.
		     */
		    if (var != MASS_FRACTION) {
		      /*
		       * load sensitivity index at this node point 
		       * this routine determines the entry in the jacobian matrix which
		       * corresponds to this BC equation and this unknown - Jac_BC is a 
		       * pointer to this entry *
		       * if it exists, add sens into matrix
		       */
		      if (Dolphin[I][var] > 0) {				  
			if (! doFullJac) {
			  ldof_var = ei->ln_to_first_dof[var][id];
			  if (ldof_var != -1) {  
			    lec->J[ieqn][pvar][ldof_eqn][ldof_var] += penalty * d_func[var];
			    lec->J[ieqn][pvar][ldof_eqn][ldof_var] *= f_time;
			  }
			} else {
			  
			  if (doMeshMapping && 
			      (var == MESH_DISPLACEMENT1 || var == MESH_DISPLACEMENT2 || var == MESH_DISPLACEMENT3)) {
			    for (j = 0; j < ei->dof[var]; j++) 
			      {
				jk = dof_map[j];
				lec->J[ieqn][pvar][ldof_eqn][jk] += penalty * d_kfunc[0][var][j];
				lec->J[ieqn][pvar][ldof_eqn][jk] *= f_time;
			      }
			  } else {
			    for (j = 0; j < ei->dof[var]; j++) 
			      {
				lec->J[ieqn][pvar][ldof_eqn][j] += penalty * d_kfunc[0][var][j];
				lec->J[ieqn][pvar][ldof_eqn][j] *= f_time;
			      }
			  }
			}
		      } else {
			/*
			 *   if variable is not defined at this node, loop
			 *   over all dof for this variable in this element
			 */
			for (j = 0; j < ei->dof[var]; j++) {
			  phi_j = bf[var]->phi[j];
			  lec->J[ieqn][pvar] [ldof_eqn][j] += penalty * d_func[var] * phi_j;
			  lec->J[ieqn][pvar] [ldof_eqn][j] *= f_time;
			}
		      }
		    } else {
		      for (w = 0; w < pd->Num_Species_Eqn; w++) {
			pvar = MAX_PROB_VAR + w;
			if (Dolphin[I][var] > 0) {
			  ldof_var = ei->ln_to_first_dof[var][id];
			  if (ldof_var != -1) {
			    lec->J[ieqn][pvar] [ldof_eqn][ldof_var] += penalty * d_func[MAX_VARIABLE_TYPES + w];
			    lec->J[ieqn][pvar] [ldof_eqn][ldof_var] *= f_time;
			  }
			}
			/* if variable is not defined at this node,
			 * loop over all dof in this element */
			else {
			  for (j = 0; j < ei->dof[var]; j++) {
			    phi_j = bf[var]->phi[j];
			    lec->J[ieqn][pvar] [ldof_eqn][j] += penalty	* d_func[MAX_VARIABLE_TYPES + w] * phi_j;
			    lec->J[ieqn][pvar] [ldof_eqn][j] *= f_time;
			  }
			}
		      } /* end of loop over species */   
		    } /* end of if MASS_FRACTION */
		  } /* end of variable exists and condition is sensitive to it */
		} /* end of loop over variable types */
	      } /* end of NEWTON */
	    } /* if (ldof_eqn != -1) */
	  } /* END of if (Res_BC != NULL), i.e. (index_eqn != -1) */
	} /* END of if COLLOCATED BC */
/*****************************************************************************/
      } /* END for (ibc = 0; (int) elem_side_bc->BC_input_id[ibc] != ...*/
/*****************************************************************************/
    } /* END if (I < num_owned_nodes) 				      */
/*****************************************************************************/
  } /* END for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) */
/*****************************************************************************/
  return(status);
} /* end of routine apply_point_colloc_bc() */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
moving_plane( int ielem_dim,
	      double *func,
	      double d_func[],
	      dbl *aa,
	      double time )
{
  fplane ( ielem_dim, func, d_func, aa );
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
fplane (int ielem_dim,
        double *func,
        double d_func[],	/* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
        dbl *aa)		/*  function parameters from data card  */
{    
/**************************** EXECUTION BEGINS *******************************/
  int  i;

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  *func  = (double) aa[3];  
  for (i=0;i<ielem_dim;i++)
    {
      d_func[MESH_DISPLACEMENT1+i] = aa[i];
      *func +=  (double) aa[i]*fv->x[i];
    }
  
} /* END of routine fplane                                                   */
/*****************************************************************************/

void 
f_fillet (const int ielem_dim,
        double *func,
        double d_func[],	/* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
        const double *p,		/*  function parameters from data card  */
        const int num_const)           /* number of passed parameters   */
{    
/**************************** EXECUTION BEGINS *******************************/
  double xpt, ypt, theta1, theta2, rad, xcen , ycen, alpha;
  double theta;

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  if(num_const != 5)
       EH(-1,"Need 5 parameters for 2D fillet geometry bc!\n");

  xpt=p[0];
  ypt=p[1];
  theta1=p[2];
  theta2=p[3];
  rad=p[4];

  alpha = 0.5*(theta2-theta1);
  xcen = xpt + (rad/sin(alpha))*cos(theta1+alpha);
  ycen = ypt + (rad/sin(alpha))*sin(theta1+alpha);

  /**   compute angle of point on curve from arc center **/

  theta = atan2(fv->x[1]-ycen,fv->x[0]-xcen);
  theta = theta > theta2-1.5*M_PIE ? theta : theta + 2*M_PIE;

  /**  use different f depending on theta  **/

  if( (theta1-0.5*M_PIE) <= theta && theta <= (theta1+alpha))
     {
      *func = (fv->x[1]-ypt)*cos(theta1) - (fv->x[0]-xpt)*sin(theta1);
      d_func[MESH_DISPLACEMENT1] =  -sin(theta1);
      d_func[MESH_DISPLACEMENT2] =  cos(theta1);

     }
  else if ( (theta1+alpha) <= theta && (theta - 0.5*M_PIE) <= theta2)
     {
      *func = (fv->x[1]-ypt)*cos(theta2) - (fv->x[0]-xpt)*sin(theta2);
      d_func[MESH_DISPLACEMENT1] = -sin(theta2);
      d_func[MESH_DISPLACEMENT2] = cos(theta2);

     }
  else
     {
      *func = SQUARE(fv->x[0]-xcen)+SQUARE(fv->x[1]-ycen)-SQUARE(rad);
      d_func[MESH_DISPLACEMENT1] = 2.*(fv->x[0]-xcen);
      d_func[MESH_DISPLACEMENT2] = 2.*(fv->x[1]-ycen);

     }

  if(ielem_dim == 3)
      d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_fillet                                                   */
/*****************************************************************************/

void 
f_double_rad (const int ielem_dim,
        double *func,
        double d_func[],	/* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
        const double *p,		/*  function parameters from data card  */
        const int num_const)           /* number of passed parameters   */
{    
/**************************** EXECUTION BEGINS *******************************/
  double xpt1, ypt1, theta1, rad1, xcen1 , ycen1, alpha1;
  double xpt2, ypt2, theta2, rad2, xcen2 , ycen2, alpha2;
  double theta1m, theta2m, th1, th2, th2t, curv_mid, rad_curv;
  double beta=0, xcirc, ycirc, dist1, dist2, dist_mid;
  int is_curved = 0;

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  if(num_const < 8)
       EH(-1,"Need at least 8 parameters for Double Rad lip geometry bc!\n");

  xpt1=p[0];
  ypt1=p[1];
  theta1=p[2];
  rad1=p[3];
  xpt2=p[4];
  ypt2=p[5];
  theta2=p[6];
  rad2=p[7];
  if(num_const >= 8)
       { curv_mid=p[8];  }
  else
       { curv_mid=0.0;  }

  is_curved = DOUBLE_NONZERO(curv_mid);

  /*  slope of middle line                */

  theta1m = atan2(ypt2-ypt1,xpt2-xpt1);
  theta1m = theta1m > theta1 ? theta1m : theta1m + 2*M_PIE;
  theta2m = atan2(ypt1-ypt2,xpt1-xpt2);
  alpha1 = 0.5*(theta1m-theta1);
  alpha2 = 0.5*(theta2-theta2m);

  xcen1 = xpt1 + (rad1/sin(alpha1))*cos(theta1+alpha1);
  ycen1 = ypt1 + (rad1/sin(alpha1))*sin(theta1+alpha1);
  xcen2 = xpt2 + (rad2/sin(alpha2))*cos(theta2m+alpha2);
  ycen2 = ypt2 + (rad2/sin(alpha2))*sin(theta2m+alpha2);

  if(is_curved)
    {
     rad_curv = 1./curv_mid;
     dist_mid = sqrt(SQUARE(xpt1-xpt2)+SQUARE(ypt1-ypt2));
     beta = asin(0.5*dist_mid*curv_mid);
     xcirc = 0.5*(xpt1+xpt2)+rad_curv*cos(beta)*sin(theta1m);
     ycirc = 0.5*(ypt1+ypt2)-rad_curv*cos(beta)*cos(theta1m);
    /**   Shift fillet centers based on curvature  **/
    /*  Using approximate distance from 90 degree corner for simplicity */

     dist1 = rad1+sqrt((rad_curv-0.5*dist_mid)*(rad_curv+0.5*dist_mid-2*rad1))
                 -sqrt((rad_curv-0.5*dist_mid)*(rad_curv+0.5*dist_mid));
     dist2 = rad2+sqrt((rad_curv-0.5*dist_mid)*(rad_curv+0.5*dist_mid-2*rad2))
                 -sqrt((rad_curv-0.5*dist_mid)*(rad_curv+0.5*dist_mid));
#if 0
     dist1 = 0;  dist2 = 0; 
#endif
     xcen1 -= dist1*cos(theta1);  ycen1 -= dist1*sin(theta1);
     xcen2 -= dist2*cos(theta2);  ycen2 -= dist2*sin(theta2);
#if 0
fprintf(stderr,"arc distances %g %g \n",dist1,dist2);
fprintf(stderr,"thetas %g %g %g %g\n",theta1, theta2, alpha1, alpha2);
fprintf(stderr,"rads %g %g %g %g\n",rad1, rad2, rad_curv,beta);
fprintf(stderr,"circle %g %g %g %g\n",xcirc,ycirc,xcen2, ycen2);
#endif
    }

  /**   compute angle of point on curve from arc center **/

  th1 = atan2(fv->x[1]-ycen1,fv->x[0]-xcen1);
  th2 = atan2(fv->x[1]-ycen2,fv->x[0]-xcen2);
  th2t = th2 > 0.0 ? th2 : th2 + 2*M_PIE;

  /**  use different f depending on theta  **/

  if( (theta1-0.5*M_PIE) <= th1 && th1 <= (theta1+alpha1))
    {
     *func = (fv->x[1]-ypt1)*cos(theta1) - (fv->x[0]-xpt1)*sin(theta1);
      d_func[MESH_DISPLACEMENT1] =  -sin(theta1);
      d_func[MESH_DISPLACEMENT2] =  cos(theta1);
/*fprintf(stderr,"DR case 1 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
    }
  else if ( (theta2-alpha2) <= th2t && (th2t - 0.5*M_PIE) <= theta2)
    {
     *func = (fv->x[1]-ypt2)*cos(theta2) - (fv->x[0]-xpt2)*sin(theta2);
      d_func[MESH_DISPLACEMENT1] = -sin(theta2);
      d_func[MESH_DISPLACEMENT2] = cos(theta2);
/*fprintf(stderr,"DR case 2 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
    }
  else if ( theta2m <= (th1+0.5*M_PIE-beta) && th1 <= (theta1-0.5*M_PIE))
    {
     *func = SQUARE(fv->x[0]-xcen1)+SQUARE(fv->x[1]-ycen1)-SQUARE(rad1);
      d_func[MESH_DISPLACEMENT1] = 2.*(fv->x[0]-xcen1);
      d_func[MESH_DISPLACEMENT2] = 2.*(fv->x[1]-ycen1);
/*fprintf(stderr,"DR case 3 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
    }
  else if ( (theta2m-0.5*M_PIE-beta) <= th2 && 
                    ((th1+0.5*M_PIE-beta) <= theta2m || (th1+0.5*M_PIE-beta) <= (theta1m-M_PIE)))
    {
     if(is_curved)
       {
        *func = SQUARE(fv->x[0]-xcirc)+SQUARE(fv->x[1]-ycirc)-SQUARE(rad_curv);
        d_func[MESH_DISPLACEMENT1] = 2.*(fv->x[0]-xcirc);
        d_func[MESH_DISPLACEMENT2] = 2.*(fv->x[1]-ycirc);
       }  
     else
       {
        *func = (fv->x[1]-ypt1)*cos(theta1m) - (fv->x[0]-xpt1)*sin(theta1m);
        d_func[MESH_DISPLACEMENT1] =  -sin(theta1m);
        d_func[MESH_DISPLACEMENT2] =  cos(theta1m);
       }
/*fprintf(stderr,"DR case 4 %g %g %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t,theta1m,theta2m);*/
    }
  else if (  (th2t - 0.5*M_PIE) >= theta2 &&  (theta2m-0.5*M_PIE-beta) >= th2 )
    {
     *func = SQUARE(fv->x[0]-xcen2)+SQUARE(fv->x[1]-ycen2)-SQUARE(rad2);
      d_func[MESH_DISPLACEMENT1] = 2.*(fv->x[0]-xcen2);
      d_func[MESH_DISPLACEMENT2] = 2.*(fv->x[1]-ycen2);
/*fprintf(stderr,"DR case 5 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
    }
  else
    {
       fprintf(stderr,"Double Rad case not found... %g %g %g %g %g\n",fv->x[0],fv->x[1],th1,th2,th2t);
    }

  if(ielem_dim == 3)
      d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_double_rad                                             */
/*****************************************************************************/

void 
f_roll_fluid (int ielem_dim,
        double *func,
        double d_func[],	/* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
        const double *p,		/*  function parameters from data card  */
        const int num_const,           /* number of passed parameters   */
        double *xsurf)           /* number of passed parameters   */
{    
/**************************** EXECUTION BEGINS *******************************/
  double roll_rad; /* roll radius */
  double origin[3];        /* roll axis origin (x,y,z) */
  double dir_angle[3];     /* axis direction angles */
  double coord[3];     /* current coordinates */
  double axis_pt[3], rad_dir[3], d_dist[3], dist, R, factor, t;
  double omega,v_dir[3], v_roll[3];
  double velo_avg = 0.0,  pgrad=0.;
  double v_solid=0., res, jac, delta, flow, eps=1.0e-8, viscinv;
  double jacinv, thick;
  int Pflag = TRUE;
  double pg_factor=1.0, tang_sgn=1.0, v_mag=0.;;

  int j,var;
#if 0
  int jvar,k;
  double dthick_dV, dthick_dP;
#endif

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  if(num_const < 7)
       EH(-1,"Need at least 7 parameters for Roll geometry bc!\n");


  roll_rad=p[0];
  origin[0] = p[1];    origin[1] = p[2];  origin[2] = p[3];
  dir_angle[0] = p[4];   dir_angle[1] = p[5];   dir_angle[2] = p[6];

/* calculate distance from interface surface to solid surface for repulsion calculations */

      coord[0] = fv->x[0];
      coord[1] = fv->x[1];
      if( ielem_dim == 3)
        { coord[2] = fv->x[2];}
      else
        { coord[2] = 0.0;}

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
    rad_dir[0] = (coord[0]-axis_pt[0])/R;
    rad_dir[1] = (coord[1]-axis_pt[1])/R;
    rad_dir[2] = (coord[2]-axis_pt[2])/R;
    dist = R - roll_rad;
    d_dist[0] = rad_dir[0]*(1.-SQUARE(dir_angle[0])/factor)
          +rad_dir[1]*(-dir_angle[1]*dir_angle[0]/factor)
          +rad_dir[2]*(-dir_angle[2]*dir_angle[0]/factor);
    d_dist[1] = rad_dir[1]*(1.-SQUARE(dir_angle[1])/factor)
          +rad_dir[0]*(-dir_angle[0]*dir_angle[1]/factor)
          +rad_dir[2]*(-dir_angle[2]*dir_angle[1]/factor);
    d_dist[2] = rad_dir[2]*(1.-SQUARE(dir_angle[2])/factor)
          +rad_dir[0]*(-dir_angle[0]*dir_angle[2]/factor)
          +rad_dir[1]*(-dir_angle[1]*dir_angle[2]/factor);


    if(num_const < 10)
       WH(-1,"ROLL_FLUID: Less than 10 parameters - reverting to roll surface!\n");
  
    omega=p[7];
/* compute velocity direction as perpendicular to both axis and radial
 *         direction.  Positive direction is determined by right hand rule */

     v_dir[0] = dir_angle[1]*rad_dir[2]-dir_angle[2]*rad_dir[1];
     v_dir[1] = dir_angle[2]*rad_dir[0]-dir_angle[0]*rad_dir[2];
     v_dir[2] = dir_angle[0]*rad_dir[1]-dir_angle[1]*rad_dir[0];

     v_roll[0] =  omega*roll_rad*v_dir[0];
     v_roll[1] =  omega*roll_rad*v_dir[1];
     v_roll[2] =  omega*roll_rad*v_dir[2];

     if( TimeIntegration == TRANSIENT && pd->e[R_MESH1] )
          {
            /* Add the mesh motion to the substrate velocity */
            v_roll[0] += fv_dot->x[0];
            v_roll[1] += fv_dot->x[1];
            v_roll[2] += fv_dot->x[2];
          }

  /* quantities specific to FLUID bcs   */

  if(num_const > 8 && p[9] >= 0.0)
     {
      dist = 0.;
      for(var=0; var < pd->Num_Dim; var ++)
       {
        /* Uses undeformed node position */
         dist += SQUARE(fv->x0[var]-xsurf[var]);
       }
         dist /= SQUARE(p[10]);
/*if(dist < 10)fprintf(stderr,"roll_fl %g %g %g\n",fv->x0[0],xsurf[0],dist);
*/

     Pflag = (int)p[11];
     velo_avg = 0.0;  pgrad=0.;  v_mag = 0.;
     for (j = 0; j < pd->Num_Dim; j++)
        {
          velo_avg += fv->stangent[0][j]*(v_roll[j] + fv->v[j]);
          v_solid += fv->stangent[0][j]*v_roll[j];
          v_mag += SQUARE(v_roll[j]);
          if(Pflag)
              {
                pgrad += fv->stangent[0][j]*fv->grad_P[j];
              }
        }
     v_mag = sqrt(v_mag);
     tang_sgn = v_solid/v_mag;
     tang_sgn = (double)SGN(v_solid/v_mag);
     velo_avg *= 0.5;
   /* sometimes the tangent/normals flip causing havoc....*/
     if(v_solid < 0)
        {
         WH(-1,"fvelo_slip: normals and tangents have flipped! - try CONTACT_LINE model\n");
         velo_avg *= tang_sgn;
         v_solid *= tang_sgn;
         pgrad *= tang_sgn;
        }

     pg_factor = 1.0;
     if(dist < 10.0)
         {pg_factor = 1.0-exp(-dist);  }
     pgrad *= pg_factor;

     flow = MAX(0.,p[9]*v_solid);
     viscinv = 1./p[8];
     thick = flow/velo_avg;
     j=0;
     do {
         res = -CUBE(thick)*viscinv*pgrad/12. + thick*velo_avg - flow;
         jac = -0.25*SQUARE(thick)*viscinv*pgrad + velo_avg;
         jacinv = 1.0/jac;
         delta = -res*jacinv;
         thick += delta;
         j++;
        } while(fabs(delta) > eps && j<20);
#if 0
      dthick_dV = -0.5*jacinv;     /*  1/h*derivative  */
      dthick_dP = CUBE(thick)*viscinv/12.*jacinv;
#endif
#if 0
fprintf(stderr,"slip %d %g %g %g %g\n",Pflag,fv->x[0],thick,flow/v_solid,velo_avg);
fprintf(stderr,"more %g %g %g %g\n",res,jac, dthick_dV,dthick_dP);
#endif
	thick = 0.;
    *func = dist - thick;
    d_func[MESH_DISPLACEMENT1] =  d_dist[0];
    d_func[MESH_DISPLACEMENT2] =  d_dist[1];
    d_func[MESH_DISPLACEMENT3] =  d_dist[2];
#if 0
    for (jvar=0; jvar<pd->Num_Dim; jvar++)
      {
        var = VELOCITY1 + jvar;
        for (k=0; k<pd->Num_Dim; k++)
          {
           d_func[var] += -thick*dthick_dV*fv->stangent[0][k];
          }
       }
#endif
#if 0
/* Mesh motion Jacobian entries   */
        for (jvar=0; jvar<ei->ielem_dim; jvar++)
          {
            var = MESH_DISPLACEMENT1 + jvar;
            if (pd->v[var])
              {
                    for (k = 0; k < pd->Num_Dim; k++)
                      {
                        d_func[var] += -thick*dthick_dV*fv->v[k]
                                *fv->stangent[0][k];
                        if(Pflag)
                          {
                          d_func[var] += -dthick_dP*pg_factor*fv->grad_P[k]*fv->stangent[0][k];
                          }
                      }
              }
          }

#endif
#if 0
   var = PRESSURE;
    if (pd->v[var])
      {
        if(Pflag )
          {
               for (k = 0; k < pd->Num_Dim; k++)
                  {
                    d_func[var] += -dthick_dP*pg_factor*fv->stangent[0][k];
                  }
          }
      }
#endif
    }	else	{
    *func = dist;
    d_func[MESH_DISPLACEMENT1] =  d_dist[0];
    d_func[MESH_DISPLACEMENT2] =  d_dist[1];
    d_func[MESH_DISPLACEMENT3] =  d_dist[2];
    }



} /* END of routine f_roll_fluid                                                   */
/*****************************************************************************/


/*****************************************************************************/

void 
fvelocity_profile (int var_flag,
        int ielem_dim,
        int velo_condition,
        double *func,
        double d_func[],       /* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
        double p[],            /* p[] are the parameters passed in through the input deck*/
        double time)           /* time at which bc's are evaluated        */
{
  if(af->Assemble_LSA_Mass_Matrix)
    d_func[var_flag] = 0.0;
  else
    d_func[var_flag] = -1.0;
  if( pd->e[R_MESH1] )
  {
  d_func[MESH_DISPLACEMENT1] = 
    dvelo_vary_fnc_d1(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
  d_func[MESH_DISPLACEMENT2] = 
    dvelo_vary_fnc_d2(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
  if (ielem_dim == 3) d_func[MESH_DISPLACEMENT3] = 
    dvelo_vary_fnc_d3(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
  }
  *func = velo_vary_fnc(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);

  *func -= fv->v[var_flag-VELOCITY1];
  
} /* END of routine fvelocity_profile                                        */
/*****************************************************************************/

void 
fvelocity_parabola (const int var_flag,
        const int ielem_dim,
        const int velo_condition,
        double *func,
        double d_func[],       /* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
        const double p[],      /* parameters passed in from the input deck*/
        const double time,           /* time at which bc's are evaluated   */
        const int num_const)           /* number of passed parameters   */
{
/*    parabolic velocity profile
 *      p[0] = coordinate1
 *      p[1] = coordinate2
 *      p[2] = flow in positive coordinate direction
 */
double coord1, coord2, qflow, gap, pre_factor, temp, expon, time_factor;
double pl_index=1.0;
int i;
	coord1 = MIN(p[0],p[1]);
	coord2 = MAX(p[1],p[0]);
	qflow = p[2];
	gap = fabs(coord2-coord1);
	pre_factor = 6.*qflow/(gap*gap*gap);
        switch (pd->CoordinateSystem) {
          case CARTESIAN:
          case CARTESIAN_2pt5D:
	       pre_factor = 6.*qflow/(gap*gap*gap);
               break;
          case CYLINDRICAL:
          case SWIRLING:
               switch (velo_condition) {
                  case U_PARABOLA_BC:
                      if(coord1 <= DBL_SMALL)
                          { pre_factor = 2.*qflow/M_PIE/SQUARE(SQUARE(coord2));}
                      else
                          { 
                           pre_factor = 2.*qflow/M_PIE/
                                       (SQUARE(coord2)-SQUARE(coord1))/
                                       (SQUARE(coord2)+SQUARE(coord1)
                       -(SQUARE(coord2)-SQUARE(coord1))/log(coord2/coord1));
                          }
                      break;
                  case V_PARABOLA_BC:
	              pre_factor = 3.*qflow/M_PIE/(gap*gap*gap);
                      break;
                  }
               break;
          default:
              EH(-1,"Velo parabola not ready for that Coordinate System yet!\n");
          }

  if(ielem_dim > 2)
     {
      EH(-1,"Velo parabola not ready for 3D yet!\n");
      return;
     }
  for(i=0;i<ielem_dim;i++)
     {
      d_func[MESH_DISPLACEMENT1+i] = 0.0;
     }

  if(af->Assemble_LSA_Mass_Matrix)
    d_func[var_flag] = 0.0;
  else
    d_func[var_flag] = -1.0;


  if( gap > DBL_SMALL)
  {
    if(num_const == 3 || p[3] == 1.0)   /*  Newtonian solution   */
      {

       switch (pd->CoordinateSystem) {
          case CARTESIAN:
          case CARTESIAN_2pt5D:
               switch (velo_condition) {
                  case U_PARABOLA_BC:
                      *func = pre_factor*(fv->x[1]-coord1)*(coord2-fv->x[1]);
                      if( pd->e[R_MESH1] )
                         {
        d_func[MESH_DISPLACEMENT2] = pre_factor*(coord1+coord2-2.*fv->x[1]);
                         }
                      break;
                  case V_PARABOLA_BC:
	              *func = pre_factor*(fv->x[0]-coord1)*(coord2-fv->x[0]);
                      if( pd->e[R_MESH1] )
                         {
       d_func[MESH_DISPLACEMENT1] = pre_factor*(coord1+coord2-2.*fv->x[0]);
                         }
                      break;
                  case W_PARABOLA_BC:
	              *func = pre_factor*(fv->x[0]-coord1)*(coord2-fv->x[0]);
                      if( pd->e[R_MESH1] )
                         {
       d_func[MESH_DISPLACEMENT1] = pre_factor*(coord1+coord2-2.*fv->x[0]);
                         }
                      break;
                  default:
                      *func =0.; 
                  }
               break;
          case CYLINDRICAL:
          case SWIRLING:
               switch (velo_condition) {
                  case U_PARABOLA_BC:
                      if(coord1 <= DBL_SMALL)
                          {
                           *func = pre_factor*(SQUARE(coord2)-SQUARE(fv->x[1]));
                           if( pd->e[R_MESH1] )
                              { 
                      d_func[MESH_DISPLACEMENT2] = pre_factor*(-2.*fv->x[1]); 
                              }
                          }
                      else
                          {
                           *func = pre_factor*(SQUARE(coord1)-SQUARE(fv->x[1])
                                    +(SQUARE(coord2)-SQUARE(coord1))*
                                    (log(fv->x[1]/coord1)/log(coord2/coord1)));
                           if( pd->e[R_MESH1] )
                              { 
                      d_func[MESH_DISPLACEMENT2] = pre_factor*(-2.*fv->x[1]
                +(SQUARE(coord2)-SQUARE(coord1))/log(coord2/coord1)/fv->x[1]);
                              }
                          }
                      break;
                  case V_PARABOLA_BC:
	              *func = pre_factor/fv->x[1]*(fv->x[0]-coord1)*(coord2-fv->x[0]);
                      if( pd->e[R_MESH1] )
                         {
       d_func[MESH_DISPLACEMENT1] = pre_factor/fv->x[1]*(coord1+coord2-2.*fv->x[0]);
       d_func[MESH_DISPLACEMENT2] = -(*func)/fv->x[1];
                         }
                      break;
                  default:
                      *func =0.; 
                  }
               break;
          }
     }
    else if(num_const > 3 )   /*  Power-law  solution   */
      {
        if(p[3] < 0.0)
            {pl_index = gn->nexp;}
        else
            {pl_index = p[3];}
       switch (pd->CoordinateSystem) {
          case CARTESIAN:
          case CARTESIAN_2pt5D:
               expon = 1.+1./pl_index;
	       pre_factor = (2.*pl_index+1.)/(pl_index +1.)*qflow/pow(gap,expon+1.);
               switch (velo_condition) {
                  case U_PARABOLA_BC:
                      temp = 2*fv->x[1]-coord1-coord2;
                      *func = pre_factor*(pow(gap,expon) - pow(fabs(temp),expon));
                      if( pd->e[R_MESH1] )
                         {
    d_func[MESH_DISPLACEMENT2] = pre_factor*(-2.*SGN(temp)*expon*pow(fabs(temp),1./pl_index));
                         }
                      break;
                  case V_PARABOLA_BC:
                      temp = 2*fv->x[0]-coord1-coord2;
                      *func = pre_factor*(pow(gap,expon) - pow(fabs(temp),expon));
                      if( pd->e[R_MESH1] )
                         {
    d_func[MESH_DISPLACEMENT1] = pre_factor*(-2.*SGN(temp)*expon*pow(fabs(temp),1./pl_index));
                         }
                      break;
                  case W_PARABOLA_BC:
                      temp = 2*fv->x[0]-coord1-coord2;
                      *func = pre_factor*(pow(gap,expon) - pow(fabs(temp),expon));
                      if( pd->e[R_MESH1] )
                         {
    d_func[MESH_DISPLACEMENT1] = pre_factor*(-2.*SGN(temp)*expon*pow(fabs(temp),1./pl_index));
                         }
                      break;
                  default:
                      *func =0.; 
                  }
               break;
          case CYLINDRICAL:
          case SWIRLING:
               expon = 1.+1./pl_index;
               switch (velo_condition) {
                  case U_PARABOLA_BC:
                      if(coord1 <= DBL_SMALL)
                          {
	                   pre_factor = (3.*pl_index+1.)/(pl_index +1.)
                                         *qflow/M_PIE/pow(gap,expon+2.);
                           *func = pre_factor*(pow(gap,expon) - pow(fv->x[1],expon));
                           if( pd->e[R_MESH1] )
                               {
                                d_func[MESH_DISPLACEMENT2] = pre_factor*
                                          (-expon*pow(fv->x[1],expon-1.));
                               }
                          }  else  {
                              EH(-1,"Power-law annulus not done yet!\n");
                          }
                      break;
                  case V_PARABOLA_BC:
	              pre_factor = (2.*pl_index+1.)/(pl_index +1.)
                                         *qflow/M_PIE/pow(gap,expon+1.);
                      temp = 2*fv->x[0]-coord1-coord2;
                      *func = pre_factor/fv->x[1]*(pow(gap,expon) - pow(fabs(temp),expon));
                      if( pd->e[R_MESH1] )
                         {
                           d_func[MESH_DISPLACEMENT1] = pre_factor/fv->x[1]*
                                 (-2.*SGN(temp)*expon*pow(fabs(temp),expon-1.));
                           d_func[MESH_DISPLACEMENT2] = -(*func)/fv->x[1];
                         }
                      break;
                  default:
                      *func =0.; 
                  }
               break;
          }
     }
  }  else   {
       *func = 0.0;
  }
/*  Add sinusoidal time-varying pieces   */
  if(num_const > 3 && num_const % 3 == 1)
     {
      time_factor = 0.;
      for(i=4 ; i<num_const ; i=i+3)
         {
          time_factor += p[i]*sin(p[i+1]*time + p[i+2]);
         }
      *func *= (1.0 + time_factor);
      for(i=0;i<ielem_dim;i++)
         {
          d_func[MESH_DISPLACEMENT1+i] *= (1.0 + time_factor);
         }
     }

  *func -= fv->v[var_flag-VELOCITY1];
  
} /* END of routine fvelocity_parabola                                        */
/*****************************************************************************/

void
fspline (int ielem_dim,
         double *func, 
         double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
         double p[],          /* parameters to parameterize temperature eqn model*/
         double time)         /* time  at which bc's are evaluated     */
{
  if(af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[MESH_DISPLACEMENT1] =
    dfncd1(fv->x[0], fv->x[1], fv->x[2], p,  time);

  d_func[MESH_DISPLACEMENT2] =
    dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3) d_func[MESH_DISPLACEMENT3] =
			dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);
    
  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);
  
} /* END of routine fspline                                                  */
/*****************************************************************************/

void
fspline_rs (int ielem_dim,
         double *func, 
         double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
         double p[],        /* parameters to parameterize temperature eqn model*/
         double time)         /* time  at which bc's are evaluated     */
{
  if(af->Assemble_LSA_Mass_Matrix)
    return;
 
  d_func[SOLID_DISPLACEMENT1] =
    dfncd1(fv->x[0], fv->x[1], fv->x[2], p,  time);

  d_func[SOLID_DISPLACEMENT2] =
    dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3) d_func[SOLID_DISPLACEMENT3] =
			dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);
     
  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);
  
} /* END of routine fspline_rs                                               */
/*****************************************************************************/

void
fTmelting(double *func,
          double d_func[],     /* dimensioned MAX_VARIABLE_TYPES + MAX_CONC */
          double a1)             /*  function parameter from data card     */
{    
  if (af->Assemble_LSA_Mass_Matrix) return;
  d_func[TEMPERATURE] = 1.;
  *func = fv->T - a1;
} /* END of routine fTmelting                                                */
/*****************************************************************************/


/******************************************************************************

 fgeneralized_dirichlet() - simple pointwise collocation contribution for BCs

     Function which adds simple contributions to boundary conditions.
     This is used as a pointwise collocation boundary condition, and has 
     flexibility built in so that the condition can be applied to any equation
     and be sensitive with respect to any unknown

     Author:          Rich Cairncross (1511)
     Date:            22 December 1994
     Revised: 
******************************************************************************/

/*
 * There are inconsistencies in the prototype declaration and how this function
 * is invoked in bc_curve.c (circa line 1170).
 *
 * Based on other boundary conditions, I believe the declaration needs to be
 * more like:
 *
 *	double func[],
 *      double d_func[][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
 *      ...
 * Resolve this at some point. -PAS
 */

int

fgeneralized_dirichlet(double *func,
		       double d_func[],	/* MAX_VARIABLE_TYPES + MAX_CONC */
		       const int gd_condition, /* denoting which condition 
						* applied */
		       const int bc_input_id,
		       const double tt, /* parameter to vary time integration 
					 * from explicit (tt = 1) to 
					 * implicit (tt = 0) */
		       const double dt) /* current time step size          */
{
  int jvar, wspec, vector_sens, b;
  int index_var;                  /* Column index into the global stiffness matrix*/
  dbl x_var;                      /* value of variable at this node */
  dbl d_x_var;                    /* sensitivity of variable to nodal unknown */
  dbl d_vect_var[DIM];            /* sensitivity of vector variable to nodal unknown */
  dbl slope;                      /* slope of interpolated function in table */
  dbl x_var_mp[1];                /* dummy variable for table lookup subroutines */
  
  if(af->Assemble_LSA_Mass_Matrix)
    return 0;

/* ---- Find variable number and species number */

  jvar = BC_Types[bc_input_id].BC_Data_Int[2];
  
  wspec = BC_Types[bc_input_id].BC_Data_Int[3];
  
  /* put value of variable in GD Condition into x_var and put it's sensitivity in d_x_var */
  index_var = load_variable( &x_var, &d_x_var, jvar, wspec, tt, dt, d_vect_var);
  
  if(jvar == SPEED)
      { vector_sens = 1;}
  else
      { vector_sens = 0;}
  /* Now add in contributions to residual vector and jacobian matrix */
  
  switch(gd_condition)
    {
    case(GD_CONST_BC):  /* x - c0 */
      
      *func = (x_var - BC_Types[bc_input_id].BC_Data_Float[0] );
      
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	            d_func[index_var+b] = d_vect_var[b];
                    }
              }  else   {
	         d_func[index_var] = d_x_var;
              }
      }
      break;
      
    case(GD_LINEAR_BC):  /* C1 x + c0 */
      
      *func = (x_var* BC_Types[bc_input_id].BC_Data_Float[1] 
	       + BC_Types[bc_input_id].BC_Data_Float[0] );
      
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1] ;
                    }
              }  else   {
	d_func[index_var] = d_x_var * BC_Types[bc_input_id].BC_Data_Float[1] ;
              }
      }
      break;

    case(GD_INVERSE_BC):  /* C1/x + c0 */
      
      *func = (BC_Types[bc_input_id].BC_Data_Float[1] / x_var 
	       + BC_Types[bc_input_id].BC_Data_Float[0] );
      
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = -d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1]/(x_var*x_var) ;
                    }
              }  else   {
	d_func[index_var] = -d_x_var * BC_Types[bc_input_id].BC_Data_Float[1]/(x_var*x_var) ;
              }
      }
      break;
      
    case(GD_PARAB_BC):  /* C2 x^2 + C1 x + c0 */
      
      *func = (x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] 
	       + x_var * BC_Types[bc_input_id].BC_Data_Float[1] 
	       + BC_Types[bc_input_id].BC_Data_Float[0] );
      
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = d_vect_var[b] * ( BC_Types[bc_input_id].BC_Data_Float[1] 
					+ 2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] );
                    }
              }  else   {
	d_func[index_var] = d_x_var * ( BC_Types[bc_input_id].BC_Data_Float[1] 
					+ 2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] );
              }
      }
      break;
    case(GD_PARAB_OFFSET_BC):  /* C2 (x - C3)^2 + C1 (x - C3) + c0 */
      
      *func = ((x_var - BC_Types[bc_input_id].BC_Data_Float[3])
	       * (x_var - BC_Types[bc_input_id].BC_Data_Float[3])
	       * BC_Types[bc_input_id].BC_Data_Float[2]
	       + (x_var - BC_Types[bc_input_id].BC_Data_Float[3])
	       * BC_Types[bc_input_id].BC_Data_Float[1]
	       + BC_Types[bc_input_id].BC_Data_Float[0] );
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = d_vect_var[b] * ( BC_Types[bc_input_id].BC_Data_Float[1]
			                + 2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) 
					* BC_Types[bc_input_id].BC_Data_Float[2] );
                    }
              }  else   {
	d_func[index_var] = d_x_var * ( BC_Types[bc_input_id].BC_Data_Float[1]
			                + 2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) 
					* BC_Types[bc_input_id].BC_Data_Float[2] );
              }
      }
      break;
    case(GD_CIRC_BC):  /* C2 ( x - C1 )^2  - c0^2 */
      /* C2 represents ellipticity and C1 represents origin */
      /* C0 is radius and should enter only one of the BC's */
      
      *func = ( BC_Types[bc_input_id].BC_Data_Float[2] * 
	       ( x_var - BC_Types[bc_input_id].BC_Data_Float[1]) * 
	       ( x_var - BC_Types[bc_input_id].BC_Data_Float[1])
	       - BC_Types[bc_input_id].BC_Data_Float[0] * BC_Types[bc_input_id].BC_Data_Float[0]);
      
      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = d_vect_var[b] * ( 2. * BC_Types[bc_input_id].BC_Data_Float[2] * (
			x_var - BC_Types[bc_input_id].BC_Data_Float[1] ) );
                    }
              }  else   {
	d_func[index_var] = d_x_var * ( 2. * BC_Types[bc_input_id].BC_Data_Float[2] * (
			x_var - BC_Types[bc_input_id].BC_Data_Float[1] ) );
              }
      }
      break;

    case(GD_POLYN_BC):  /* up to 6th order polynomial */
                        /* C6 x^6 + C5 x^5 + C4 x^4 + C3 x^3 + C2 x^2 + C1 x + C0 */

      *func = (x_var * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6]
	       + x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5]
	       + x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4]
	       + x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3]
               + x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] 
	       + x_var * BC_Types[bc_input_id].BC_Data_Float[1] 
	       + BC_Types[bc_input_id].BC_Data_Float[0] );

      /* printf("POLYN fit X,F = %f %f\n", x_var, *func); */

      if (af->Assemble_Jacobian) {
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	d_func[index_var+b] = d_vect_var[b] * ( BC_Types[bc_input_id].BC_Data_Float[1] 
				       + 2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]
			       + 3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3]
		       + 4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4]
	       + 5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5]
       + 6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6] );
                    }
              }  else   {
	d_func[index_var] = d_x_var * ( BC_Types[bc_input_id].BC_Data_Float[1] 
				       + 2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]
			       + 3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3]
		       + 4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4]
	       + 5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5]
       + 6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6] );
              }
      }
      break;

    case(GD_TABLE_BC):

      *func = BC_Types[bc_input_id].BC_Data_Float[0];

      x_var_mp[0]=x_var;
      *func *= interpolate_table( BC_Types[bc_input_id].table, x_var_mp, &slope, NULL );

      if (af->Assemble_Jacobian) 
	{
          if (vector_sens)
              {
                for(b=0 ; b<DIM  ; b++)	{
	  d_func[index_var+b] = BC_Types[bc_input_id].BC_Data_Float[0]*slope*d_vect_var[b];
                    }
              }  else   {
	  d_func[index_var] = BC_Types[bc_input_id].BC_Data_Float[0]*slope*d_x_var;
              }
	}

      break;
                  
    default:
      return(-1);
    }

  return(0);
} /* END of routine fgeneralized_dirichlet                                   */
/*****************************************************************************/

void
fmesh_constraint(double *func,
		 double d_func[],
		 const int bc_input_id)
{
  EH(-1, "CGM not supported, MESH_CONSTRAINT_BC");
/*#endif  */

} /* END of routine fmesh_constraint                                   */


/*****************************************************************************/
/*****************************************************************************/

int
load_variable (double *x_var,        /* variable value */
               double *d_x_var,      /* sensitivities of variable value */
               int jvar,             /* variable number */
               int wspec,            /* species number */
               double tt,            /* parameter to vary time integration from 
                                        explicit (tt = 1) to implicit (tt = 0) */
               double dt,            /* current time step size */
               double d_vect_var[])	/* vector sensitivities  */

/******************************************************************************

   Function which calculates the value of a chosen variable at the current point
    (using the field variables); available variable names listed in mm_names.h

   Author:          Rich Cairncross (1511)
   Date:            22 MAY 1995
   Revised: 
******************************************************************************/

{
  int var=-1, b;
  *x_var = 0.;
  *d_x_var = 0.;

  memset(d_vect_var, 0, DIM*sizeof(double) );

  if (jvar >= D_VEL1_DT && jvar <= D_P_DT) {
    if ( pd->TimeIntegration == STEADY ) EH(-1, "Unsteady GD for Steady problem");
  }

  /* ---- Find variable value, sensitivities, and species number */

  switch(jvar)
    {
    case VELOCITY1:
      *x_var = fv->v[0];
      var = VELOCITY1;
      *d_x_var = 1.;
      break;
    case VELOCITY2:
      *x_var = fv->v[1];
      var = VELOCITY2;
      *d_x_var = 1.;
      break;
    case VELOCITY3:
      *x_var = fv->v[2];
      var = VELOCITY3;
      *d_x_var = 1.;
      break;
    case PVELOCITY1:
      *x_var = fv->pv[0];
      var = PVELOCITY1;
      *d_x_var = 1.;
      break;
    case PVELOCITY2:
      *x_var = fv->pv[1];
      var = PVELOCITY2;
      *d_x_var = 1.;
      break;
    case PVELOCITY3:
      *x_var = fv->pv[2];
      var = PVELOCITY3;
      *d_x_var = 1.;
      break;
    case TEMPERATURE:
      *x_var = fv->T;
      var = TEMPERATURE;
      *d_x_var = 1.;
      break;
    case VOLTAGE:
      *x_var = fv->V;
      var = VOLTAGE;
      *d_x_var = 1.;
      break;
    case SURF_CHARGE:
      *x_var = fv->qs;
      var = SURF_CHARGE;
      *d_x_var = 1.;
      break;
   case SHELL_CURVATURE:
      *x_var = fv->sh_K;
      var = SHELL_CURVATURE;
      *d_x_var = 1.;
      break;
   case SHELL_CURVATURE2:
      *x_var = fv->sh_K2;
      var = SHELL_CURVATURE2;
      *d_x_var = 1.;
      break;
    case SHELL_TENSION:
      *x_var = fv->sh_tens;
      var = SHELL_TENSION;
      *d_x_var = 1.;
      break;
    case SHELL_X:
      *x_var = fv->sh_x;
      var = SHELL_X;
      *d_x_var = 1.;
      break;
    case SHELL_Y:
      *x_var = fv->sh_y;
      var = SHELL_Y;
      *d_x_var = 1.;
      break;
    case SHELL_USER:
      *x_var = fv->sh_u;
      var = SHELL_USER;
      *d_x_var = 1.;
      break;
    case SHELL_ANGLE1:
      *x_var = fv->sh_ang[0];
      var = SHELL_ANGLE1;
      *d_x_var = 1.;
      break;
    case SHELL_ANGLE2:
      *x_var = fv->sh_ang[1];
      var = SHELL_ANGLE2;
      *d_x_var = 1.;
      break;
    case SHELL_SURF_DIV_V:
      *x_var = fv->div_s_v;
      var = SHELL_SURF_DIV_V;
      *d_x_var = 1.;
      break;
    case SHELL_SURF_CURV:
      *x_var = fv->curv;
      var = SHELL_SURF_CURV;
      *d_x_var = 1.;
      break;
    case N_DOT_CURL_V:
      *x_var = fv->n_dot_curl_s_v;
      var = N_DOT_CURL_V;
      *d_x_var = 1.;
      break;
    case GRAD_S_V_DOT_N1:
      *x_var = fv->grad_v_dot_n[0];
      var = GRAD_S_V_DOT_N1;
      *d_x_var = 1.;
      break;
    case GRAD_S_V_DOT_N2:
      *x_var = fv->grad_v_dot_n[1];
      var = GRAD_S_V_DOT_N2;
      *d_x_var = 1.;
      break;
    case GRAD_S_V_DOT_N3:
      *x_var = fv->grad_v_dot_n[2];
      var = GRAD_S_V_DOT_N3;
      *d_x_var = 1.;
      break;
    case SHELL_DIFF_FLUX:
      *x_var = fv->sh_J;
      var = SHELL_DIFF_FLUX;
      *d_x_var = 1.;
      break;
    case SHELL_DIFF_CURVATURE:
      *x_var = fv->sh_Kd;
      var = SHELL_DIFF_CURVATURE;
      *d_x_var = 1.;
      break;
    case SHELL_NORMAL1:
      *x_var = fv->n[0];
      var = SHELL_NORMAL1;
      *d_x_var = 1.;
      break;
    case SHELL_NORMAL2:
      *x_var = fv->n[1];
      var = SHELL_NORMAL2;
      *d_x_var = 1.;
      break;
    case SHELL_NORMAL3:
      *x_var = fv->n[2];
      var = SHELL_NORMAL3;
      *d_x_var = 1.;
      break;
    case ACOUS_PREAL:
      *x_var = fv->apr;
      var = ACOUS_PREAL;
      *d_x_var = 1.;
      break;
    case ACOUS_PIMAG:
      *x_var = fv->api;
      var = ACOUS_PIMAG;
      *d_x_var = 1.;
      break;
    case POR_SINK_MASS:
      *x_var = fv->sink_mass;
      var = POR_SINK_MASS;
      *d_x_var = 1.;
      break;
    case ACOUS_REYN_STRESS:
      *x_var = fv->ars;
      var = ACOUS_REYN_STRESS;
      *d_x_var = 1.;
      break;
    case SHELL_BDYVELO:
      *x_var = fv->sh_bv;
      var = SHELL_BDYVELO;
      *d_x_var = 1.;
      break;
    case SHELL_LUBP:
      *x_var = fv->sh_p;
      var = SHELL_LUBP;
      *d_x_var = 1.;
      break;
    case SHELL_FILMP:
      *x_var = fv->sh_fp;
      var = SHELL_FILMP;
      *d_x_var = 1.;
      break;
    case SHELL_FILMH:
      *x_var = fv->sh_fh;
      var = SHELL_FILMH;
      *d_x_var = 1.;
      break;
    case SHELL_PARTC:
      *x_var = fv->sh_pc;
      var = SHELL_PARTC;
      *d_x_var = 1.;
      break;
    case LUBP:
      *x_var = fv->lubp;
      var = LUBP;
      *d_x_var = 1.;
      break;
    case LUBP_2:
      *x_var = fv->lubp_2;
      var = LUBP_2;
      *d_x_var = 1.;
      break;
    case SHELL_SAT_CLOSED:
      *x_var = fv->sh_sat_closed;
      var = SHELL_SAT_CLOSED;
      *d_x_var = 1.;
      break;
    case SHELL_PRESS_OPEN:
      *x_var = fv->sh_p_open;
      var = SHELL_PRESS_OPEN;
      *d_x_var = 1.;
      break;
    case SHELL_PRESS_OPEN_2:
      *x_var = fv->sh_p_open_2;
      var = SHELL_PRESS_OPEN_2;
      *d_x_var = 1.;
      break;
    case SHELL_TEMPERATURE:
      *x_var = fv->sh_t;
      var = SHELL_TEMPERATURE;
      *d_x_var = 1.;
      break;
    case SHELL_DELTAH:
      *x_var = fv->sh_dh;
      var = SHELL_DELTAH;
      *d_x_var = 1.;
      break;
    case SHELL_LUB_CURV:
      *x_var = fv->sh_l_curv;
      var = SHELL_LUB_CURV;
      *d_x_var = 1.;
      break;
    case SHELL_LUB_CURV_2:
      *x_var = fv->sh_l_curv_2;
      var = SHELL_LUB_CURV_2;
      *d_x_var = 1.;
      break;
    case SHELL_SAT_GASN:
      *x_var = fv->sh_sat_gasn;
      var = SHELL_SAT_GASN;
      *d_x_var = 1.;
      break;
    case SHELL_SHEAR_TOP:
      *x_var = fv->sh_shear_top;
      var = SHELL_SHEAR_TOP;
      *d_x_var = 1.;
      break;
    case SHELL_SHEAR_BOT:
      *x_var = fv->sh_shear_bot;
      var = SHELL_SHEAR_BOT;
      *d_x_var = 1.;
      break;
    case SHELL_CROSS_SHEAR:
      *x_var = fv->sh_cross_shear;
      var = SHELL_CROSS_SHEAR;
      *d_x_var = 1.;
      break;
    case TFMP_PRES:
      *x_var = fv->tfmp_pres;
      var = TFMP_PRES;
      *d_x_var = 1;
      break;
    case TFMP_SAT:
      *x_var = fv->tfmp_sat;
      var = TFMP_SAT;
      *d_x_var = 1;
      break;
    case MAX_STRAIN:
      *x_var = fv->max_strain;
      var = MAX_STRAIN;
      *d_x_var = 1.;
      break;
    case CUR_STRAIN:
      *x_var = fv->cur_strain;
      var = CUR_STRAIN;
      *d_x_var = 1.;
      break;
    case LIGHT_INTP:
      *x_var = fv->poynt[0];
      var = LIGHT_INTP;
      *d_x_var = 1.;
      break;
    case LIGHT_INTM:
      *x_var = fv->poynt[1];
      var = LIGHT_INTM;
      *d_x_var = 1.;
      break;
    case LIGHT_INTD:
      *x_var = fv->poynt[2];
      var = LIGHT_INTD;
      *d_x_var = 1.;
      break;
    case RESTIME:
      *x_var = fv->restime;
      var = RESTIME;
      *d_x_var = 1.;
      break;  
    case MASS_FRACTION:
      *x_var = fv->c[wspec];
      var = MASS_FRACTION;
      *d_x_var = 1.;
      break;
    case MESH_DISPLACEMENT1:
      *x_var = fv->d[0];
      var = MESH_DISPLACEMENT1;
      *d_x_var = 1.;
      break;
    case MESH_DISPLACEMENT2:
      *x_var = fv->d[1];
      var = MESH_DISPLACEMENT2;
      *d_x_var = 1.;
      break;
    case MESH_DISPLACEMENT3:
      *x_var = fv->d[2];
      var = MESH_DISPLACEMENT3;
      *d_x_var = 1.;
      break;
    case SOLID_DISPLACEMENT1:
      *x_var = fv->d_rs[0];
      var = SOLID_DISPLACEMENT1;
      *d_x_var = 1.;
      break;
    case SOLID_DISPLACEMENT2:
      *x_var = fv->d_rs[1];
      var = SOLID_DISPLACEMENT2;
      *d_x_var = 1.;
      break;
    case SOLID_DISPLACEMENT3:
      *x_var = fv->d_rs[2];
      var = SOLID_DISPLACEMENT3;
      *d_x_var = 1.;
      break;
    case SURFACE:
      EH(-1,"SURFACE variables not defined yet");
      break;
    case PRESSURE:
      *x_var = fv->P;
      var = PRESSURE;
      *d_x_var = 1.;
      break;
    case SHEAR_RATE:
      *x_var = fv->SH;
      var = SHEAR_RATE;
      *d_x_var = 1.;
      break;

   case EXT_VELOCITY:
      *x_var = fv->ext_v;
      var = EXT_VELOCITY;
      *d_x_var = 1.;
      break;

   case EFIELD1:
      *x_var = fv->E_field[0];
      var = EFIELD1;
      *d_x_var = 1.;
      break;

   case EFIELD2:
      *x_var = fv->E_field[1];
      var = EFIELD2;
      *d_x_var = 1.;
      break;
    case EFIELD3:
      *x_var = fv->E_field[2];
      var = EFIELD3;
      *d_x_var = 1.;
      break;
      
    case ENORM:
      *x_var = fv->Enorm;
      var = ENORM;
      *d_x_var = 1.;
      break;

    case CURVATURE:
      *x_var = fv->H;
      var = CURVATURE;
      *d_x_var = 1.;
      break;

    case NORMAL1:
      *x_var = fv->n[0];
      var = NORMAL1;
      *d_x_var = 1;
      break;

    case NORMAL2:
      *x_var = fv->n[1];
      var = NORMAL2;
      *d_x_var = 1;
      break;

    case NORMAL3:
      *x_var = fv->n[2];
      var = NORMAL3;
      *d_x_var = 1;
      break;
      
    case POLYMER_STRESS11:
      *x_var = fv->S[0][0][0];
      var = POLYMER_STRESS11;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12:
      *x_var = fv->S[0][0][1];
      var = POLYMER_STRESS12;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22:
      *x_var = fv->S[0][1][1];
      var = POLYMER_STRESS22;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13:
      *x_var = fv->S[0][0][2];
      var = POLYMER_STRESS13;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23:
      *x_var = fv->S[0][1][2];
      var = POLYMER_STRESS23;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33:
      *x_var = fv->S[0][2][2];
      var = POLYMER_STRESS33;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_1:
      *x_var = fv->S[1][0][0];
      var = POLYMER_STRESS11_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_1:
      *x_var = fv->S[1][0][1];
      var = POLYMER_STRESS12_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_1:
      *x_var = fv->S[1][1][1];
      var = POLYMER_STRESS22_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_1:
      *x_var = fv->S[1][0][2];
      var = POLYMER_STRESS13_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_1:
      *x_var = fv->S[1][1][2];
      var = POLYMER_STRESS23_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_1:
      *x_var = fv->S[1][2][2];
      var = POLYMER_STRESS33_1;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_2:
      *x_var = fv->S[2][0][0];
      var = POLYMER_STRESS11_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_2:
      *x_var = fv->S[2][0][1];
      var = POLYMER_STRESS12_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_2:
      *x_var = fv->S[2][1][1];
      var = POLYMER_STRESS22_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_2:
      *x_var = fv->S[2][0][2];
      var = POLYMER_STRESS13_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_2:
      *x_var = fv->S[2][1][2];
      var = POLYMER_STRESS23_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_2:
      *x_var = fv->S[2][2][2];
      var = POLYMER_STRESS33_2;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_3:
      *x_var = fv->S[3][0][0];
      var = POLYMER_STRESS11_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_3:
      *x_var = fv->S[3][0][1];
      var = POLYMER_STRESS12_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_3:
      *x_var = fv->S[3][1][1];
      var = POLYMER_STRESS22_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_3:
      *x_var = fv->S[3][0][2];
      var = POLYMER_STRESS13_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_3:
      *x_var = fv->S[3][1][2];
      var = POLYMER_STRESS23_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_3:
      *x_var = fv->S[3][2][2];
      var = POLYMER_STRESS33_3;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_4:
      *x_var = fv->S[4][0][0];
      var = POLYMER_STRESS11_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_4:
      *x_var = fv->S[4][0][1];
      var = POLYMER_STRESS12_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_4:
      *x_var = fv->S[4][1][1];
      var = POLYMER_STRESS22_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_4:
      *x_var = fv->S[4][0][2];
      var = POLYMER_STRESS13_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_4:
      *x_var = fv->S[4][1][2];
      var = POLYMER_STRESS23_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_4:
      *x_var = fv->S[4][2][2];
      var = POLYMER_STRESS33_4;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_5:
      *x_var = fv->S[5][0][0];
      var = POLYMER_STRESS11_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_5:
      *x_var = fv->S[5][0][1];
      var = POLYMER_STRESS12_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_5:
      *x_var = fv->S[5][1][1];
      var = POLYMER_STRESS22_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_5:
      *x_var = fv->S[5][0][2];
      var = POLYMER_STRESS13_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_5:
      *x_var = fv->S[5][1][2];
      var = POLYMER_STRESS23_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_5:
      *x_var = fv->S[5][2][2];
      var = POLYMER_STRESS33_5;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_6:
      *x_var = fv->S[6][0][0];
      var = POLYMER_STRESS11_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_6:
      *x_var = fv->S[6][0][1];
      var = POLYMER_STRESS12_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_6:
      *x_var = fv->S[6][1][1];
      var = POLYMER_STRESS22_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_6:
      *x_var = fv->S[6][0][2];
      var = POLYMER_STRESS13_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_6:
      *x_var = fv->S[6][1][2];
      var = POLYMER_STRESS23_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_6:
      *x_var = fv->S[6][2][2];
      var = POLYMER_STRESS33_6;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS11_7:
      *x_var = fv->S[7][0][0];
      var = POLYMER_STRESS11_7;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS12_7:
      *x_var = fv->S[7][0][1];
      var = POLYMER_STRESS12_7;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS22_7:
      *x_var = fv->S[7][1][1];
      var = POLYMER_STRESS22_7;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS13_7:
      *x_var = fv->S[7][0][2];
      var = POLYMER_STRESS13_7;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS23_7:
      *x_var = fv->S[7][1][2];
      var = POLYMER_STRESS23_7;
      *d_x_var = 1.;
      break;

    case POLYMER_STRESS33_7:
      *x_var = fv->S[7][2][2];
      var = POLYMER_STRESS33_7;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT11:
      *x_var = fv->G[0][0];
      var = VELOCITY_GRADIENT11;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT12:
      *x_var = fv->G[0][1];
      var = VELOCITY_GRADIENT12;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT21:
      *x_var = fv->G[1][0];
      var = VELOCITY_GRADIENT21;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT22:
      *x_var = fv->G[1][1];
      var = VELOCITY_GRADIENT22;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT13:
      *x_var = fv->G[0][2];
      var = VELOCITY_GRADIENT13;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT23:
      *x_var = fv->G[1][2];
      var = VELOCITY_GRADIENT23;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT31:
      *x_var = fv->G[2][0];
      var = VELOCITY_GRADIENT31;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT32:
      *x_var = fv->G[2][1];
      var = VELOCITY_GRADIENT32;
      *d_x_var = 1.;
      break;

    case VELOCITY_GRADIENT33:
      *x_var = fv->G[2][2];
      var = VELOCITY_GRADIENT33;
      *d_x_var = 1.;
      break;

    case PHASE1:
    case PHASE2:
    case PHASE3:
    case PHASE4:
    case PHASE5:

      b = jvar - PHASE1;
      *x_var = fv->pF[b];
      var = PHASE1 + b;
      *d_x_var = 1.;
      break;

      /* if variable type is mesh position **not** mesh displacement*/
    case MESH_POSITION1:
      *x_var = fv->x[0];
      var = MESH_DISPLACEMENT1;
      *d_x_var = 1.;
      break;
    case MESH_POSITION2:
      *x_var = fv->x[1];
      var = MESH_DISPLACEMENT2;
      *d_x_var = 1.;
      break;
    case MESH_POSITION3:
      *x_var = fv->x[2];
      var = MESH_DISPLACEMENT3;
      *d_x_var = 1.;
      break;
    case SOLID_POSITION1:
      *x_var = fv->x[0];
      var = SOLID_DISPLACEMENT1;
      *d_x_var = 1.;
      break;
    case SOLID_POSITION2:
      *x_var = fv->x[1];
      var = SOLID_DISPLACEMENT2;
      *d_x_var = 1.;
      break;
    case SOLID_POSITION3:
      *x_var = fv->x[2];
      var = SOLID_DISPLACEMENT3;
      *d_x_var = 1.;
      break;
    case POR_LIQ_PRES:
      *x_var = fv->p_liq;
      var = POR_LIQ_PRES;
      *d_x_var = 1.;
      break;
    case POR_GAS_PRES:
      *x_var = fv->p_gas;
      var = POR_GAS_PRES;
      *d_x_var = 1.;
      break;
    case POR_POROSITY:
      *x_var = fv->porosity;
      var = POR_GAS_PRES;
      *d_x_var = 1.;
      break;
     /* adding velocity magnitude, i.e. SPEED  */
    case SPEED:
      for(b=0 ; b<pd->Num_Dim ; b++)	{
          *x_var += SQUARE(fv->v[b]);
          }
      if(pd->CoordinateSystem == SWIRLING || 
         pd->CoordinateSystem == PROJECTED_CARTESIAN ||
         pd->CoordinateSystem == CARTESIAN_2pt5D)
          { *x_var += SQUARE(fv->v[pd->Num_Dim]);  }
      *x_var = sqrt(*x_var);
      var = VELOCITY1;
      *d_x_var = 1./(*x_var);
      for(b=0 ; b<pd->Num_Dim ; b++)	{
          d_vect_var[b] += fv->v[b]*(*d_x_var);
          }
      if(pd->CoordinateSystem == SWIRLING || 
         pd->CoordinateSystem == PROJECTED_CARTESIAN ||
         pd->CoordinateSystem == CARTESIAN_2pt5D)
          { d_vect_var[pd->Num_Dim] += fv->v[pd->Num_Dim]*(*d_x_var);  }
      break;

      /* if variable type is a time derivative */
    case D_VEL1_DT:
      *x_var = fv_dot->v[0];
      var = VELOCITY1;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_VEL2_DT:
      *x_var = fv_dot->v[1];
      var = VELOCITY2;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_VEL3_DT:
      *x_var = fv_dot->v[2];
      var = VELOCITY3;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_T_DT:
      *x_var = fv_dot->T;
      var = TEMPERATURE;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_Y_DT:
      *x_var = fv_dot->c[wspec];
      var = MASS_FRACTION;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_X1_DT:
      *x_var = fv_dot->d[0];
      var = MESH_DISPLACEMENT1;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_X2_DT:
      *x_var = fv_dot->d[1];
      var = MESH_DISPLACEMENT2;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_X3_DT:
      *x_var = fv_dot->d[2];
      var = MESH_DISPLACEMENT3;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    case D_S_DT:
      EH(-1,"SURFACE variables not defined yet");
      break;
    case D_P_DT:
      *x_var = fv_dot->P;
      var = PRESSURE;
      *d_x_var = (1. + 2. * tt) / dt;
      break;
    default:
      EH(-1, "Illegal option in load_variable");
    } /* end of switch on jvar */

  if (wspec != 0 && var != MASS_FRACTION && var != POR_GAS_PRES) {
    EH(-1, "Non-zero species number for wrong variable");
  }
  if (var == MASS_FRACTION) var = MAX_VARIABLE_TYPES + wspec;

  return(var);
} /* END of routine load_variable()                                          */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
bc_eqn_index(int id,               /* local node number                 */
	     int I,                /* processor node number             */
	     int bc_input_id,      /* boundary condition number         */
	     int curr_matID,       /* Current material ID */
	     int kdir,       /* coordinate index for mesh or velocity
			      * - normally this is zero unless a
			      *   bc is a vector condition                 */
	     int *eqn,       /* eqn to which this condition is applied     */
	     int *matID_retn, /* material ID to apply this eqn on           */
	     VARIABLE_DESCRIPTION_STRUCT **vd_retn)
    
    /*************************************************************************
     *
     * bc_eqn_index():
     *
     * Check to see if this BC on this node is applicable
     *     (i.e. no other overriding Dirichlet conditions),
     * Find the global unknown number for applying this condition.
     *
     * NOTE: This must be called after the element basis functions have
     *       been set up.
     *
     * Input
     * -------
     *
     * id = local node number 
     * I  = processor node number
     * bc_input_id = boundary condition number
     * curr_matID = Current material index of the element that we
     *              calling from
     * BC_desc  = description structure for current boundary condition
     * 
     *
     * Output
     * -------
     *  *eqn = Variable type of the equation to which this boundary condition
     *        is applied.
     *  *matID_retn = Material number of the equation to which this
     *                boundary condition is applied.
     *
     * Return
     * --------
     *   This function returns the index in the solution vector of the
     *   equation to which this boundary condition will be applied.
     *   If the boundary condition is not to be applied, this function
     *   return the value of -1.
     *
     *  *matID_retn = returns the material id corresponding to the unknown
     *                on which this boundary condition is applied. Even if
     *                the underyling variable type has a value of -1, this
     *                variable will not be -1.
     *************************************************************************/
{
  int ieqn,  i_calc,  jeqn, bc_node, offset_j;
  int matID, index_eqn, num, node_offset;
  int  i, index;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd, *vd2 = NULL;
  nv = node->Nodal_Vars_Info;
    
  /*
   *  Find equation number and species number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  matID = curr_matID;
  
  /* need to change equation index for mesh or velocity */
  if (kdir != 0) {
    if      (ieqn == R_MESH1)     ieqn += kdir;
    else if (ieqn == R_MOMENTUM1) ieqn += kdir;
    else if (ieqn == R_SOLID1)    ieqn += kdir;
    else if (ieqn == R_LAGR_MULT1)ieqn += kdir;
    else EH(-1,"Can't have a rotated vector BC!");
  }

  /*
   *  Find out the offset under non-conflicting conditions
   *  If a boundary condition isn't to be applied under non-
   *  conflicting conditions, it won't be applied under
   *  conflicting conditions. So exit here if bc is not
   *  to be applied. 
   */
  node_offset = find_bc_unk_offset(bc, curr_matID, I, kdir,
				   &matID, &vd);
  if (node_offset < 0) {
    return -1;
  }

  /*
   * check the processor node number, I, against BC_dup_nodes[]
   * to see if this node is at an intersection of side-sets
   *
   * HKM -> the search below seems like a time-waste. In lieu of
   *        that, we could set up an added  bit field in
   *        Node_Info struct or elsewhere where we could store the
   *        results of this search (or just the need to do the
   *        search).
   */
  bc_node = in_list(I, 0, BC_dup_ptr+1, BC_dup_nodes);
  if (bc_node != -1) {
    /* 
     * current node is in the BC_duplication list.
     * The current boundary condition may have been deleted.
     * If the current boundary condition is a rotated boundary
     * condition, it may have been moved to another coordinate
     * direction.
     */
    i_calc = -1;
    if (ieqn <= LAST_REAL_EQ) {
      i_calc = search_bc_dup_list(bc_input_id, 
				  BC_dup_list[bc_node][node_offset]);

    } else if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
      /*
       * If it's a rotated boundary condition, the boundary
       * condition could have been moved to another coordinate
       * within check_for_bc_conflicts2D(). Therefore, we need
       * to check all coordinate directions for the presence of the
       * current boundary condition.
       */
      for (jeqn = R_MESH1; jeqn <= R_MESH3; jeqn++) {
	offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
	if (offset_j >= 0) {
	  i_calc = search_bc_dup_list(bc_input_id, 
				      BC_dup_list[bc_node][offset_j]);
	  if (i_calc != -1) {
	    node_offset = offset_j;
	    ieqn = jeqn;
	    vd = vd2;
	    break;
	  }
	}
      }

    } else if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
      for (jeqn = R_MOMENTUM1; jeqn <= R_MOMENTUM3; jeqn++) {
	offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
	if (offset_j >= 0) {
	  i_calc = search_bc_dup_list(bc_input_id, 
				      BC_dup_list[bc_node][offset_j]);
	  if (i_calc != -1) {
	    node_offset = offset_j;
	    ieqn = jeqn;
	    vd = vd2;
	    break;
	  }
	}
      }
    } else if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
      for (jeqn = R_SOLID1; jeqn <= R_SOLID3; jeqn++) {
	offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
	if (offset_j >= 0) {
	  i_calc = search_bc_dup_list(bc_input_id, 
				      BC_dup_list[bc_node][offset_j]);
	  if (i_calc != -1) {
	    node_offset = offset_j;
	    ieqn = jeqn;
	    vd = vd2;
	    break;
	  }
	}
      }
    } else {
      EH(-1,"Equation not in list ");
    }

    /*
     * The current boundary condition has been deleted from application
     * at the current node due to a conflict with another boundary
     * condition. Return a -1.
     */
    if (i_calc == -1) {
      return -1;
    }
  }
  /*
   * Not at an intersection:
   *    For boundary conditions which refer to rotated
   *    coordinates, assign an equation number to them
   *    base upon a (normal, tang1, tang2) = ( x, y, z)
   *    mapping.
   */

  if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
    ieqn = ieqn - R_MESH_NORMAL + R_MESH1;
  }
  if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
    ieqn = ieqn - R_SOLID_NORMAL + R_SOLID1;
  }
  if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
    ieqn = ieqn - R_MOM_NORMAL + R_MOMENTUM1;
  }

  /*
   * If the regular unknown equation is replaced by a Dirichlet 
   * condition at this node, return -1
   */
  if (node->DBC) {
    if (((int) node->DBC[node_offset]) != -1) {
      return -1;
    }
  } 

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown + node_offset;
  
  /*
   * HKM -> we can put this section into  
   *        a subroutine
   */
  if (vd_retn) {
    *vd_retn = NULL;
    num = nv->Num_Var_Desc_Per_Type[ieqn];
    for (i = 0; i < num; i++) {
      index = nv->Var_Type_Index[ieqn][i];
      vd = nv->Var_Desc_List[index];
      if (matID == vd->MatID ||
	  (! (*vd_retn) && vd->MatID == -1)) {
	*vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index                                             */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
bc_eqn_index_stress(int id,               /* local node number                 */
	            int I,                /* processor node number             */
	            int bc_input_id,      /* boundary condition number         */
	            int curr_matID,       /* Current material ID */
	            int kdir,             /* coordinate index for stress components */
                    int mode,       /* Stress mode number */
	            int *eqn,       /* eqn to which this condition is applied     */
	            int *matID_retn, /* material ID to apply this eqn on           */
	            VARIABLE_DESCRIPTION_STRUCT **vd_retn)

    /*************************************************************************
     *
     * bc_eqn_index_stress():
     *
     * Check to see if this BC on this node is applicable
     *     (i.e. no other overriding Dirichlet conditions),
     * Find the global unknown number for applying this condition. This is clone
     * of bc_eqn_index repurposed to handle stress equations
     *
     * NOTE: This must be called after the element basis functions have
     *       been set up.
     *
     * Input
     * -------
     *
     * id = local node number
     * I  = processor node number
     * bc_input_id = boundary condition number
     * curr_matID = Current material index of the element that we
     *              calling from
     * BC_desc  = description structure for current boundary condition
     *
     *
     * Output
     * -------
     *  *eqn = Variable type of the equation to which this boundary condition
     *        is applied.
     *  *matID_retn = Material number of the equation to which this
     *                boundary condition is applied.
     *
     * Return
     * --------
     *   This function returns the index in the solution vector of the
     *   equation to which this boundary condition will be applied.
     *   If the boundary condition is not to be applied, this function
     *   return the value of -1.
     *
     *  *matID_retn = returns the material id corresponding to the unknown
     *                on which this boundary condition is applied. Even if
     *                the underyling variable type has a value of -1, this
     *                variable will not be -1.
     *************************************************************************/
{
  int ieqn,  i_calc,  bc_node;
  int matID, index_eqn, num, node_offset;
  int  i, index;
  int ndofs_ieqn;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  nv = node->Nodal_Vars_Info;

  /*
   *  Find equation number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  /* Sanity check */
  if (ieqn != R_STRESS11) EH(-1,"You can't be here");

  matID = curr_matID;

  switch (mode) {

    case 0:
      ieqn = R_STRESS11 + kdir;
      break;

    case 1:
      ieqn = R_STRESS11_1 + kdir;
      break;

    case 2:
      ieqn = R_STRESS11_2 + kdir;
      break;

    case 3:
      ieqn = R_STRESS11_3 + kdir;
      break;


    case 4:
      ieqn = R_STRESS11_4 + kdir;
      break;

    case 5:
      ieqn = R_STRESS11_5 + kdir;
      break;

    case 6:
      ieqn = R_STRESS11_6 + kdir;
      break;

    case 7:
      ieqn = R_STRESS11_7 + kdir;
      break;

    default:
      EH(-1,"Maximum 8 modes allowed here!");
      break;
  }


  /*
   * If the equation has zero degrees of freedom at this node
   * return. Otherwise, get the number of variable description
   * structures with this variable type. Don't count more
   * than one species unknown per material at this point.
   */
  if ((ndofs_ieqn = get_nv_ndofs_modMF(nv, ieqn)) <= 0) {
    return -1;
  }


  /*
   *  Find out the offset under non-conflicting conditions
   *  If a boundary condition isn't to be applied under non-
   *  conflicting conditions, it won't be applied under
   *  conflicting conditions. So exit here if bc is not
   *  to be applied. 
   */
  node_offset = get_nodal_unknown_offset(nv, ieqn, matID, 0, &vd);
  if (node_offset < 0) {
    return -1;
  }

  /*
   * check the processor node number, I, against BC_dup_nodes[]
   * to see if this node is at an intersection of side-sets
   *
   * HKM -> the search below seems like a time-waste. In lieu of
   *        that, we could set up an added  bit field in
   *        Node_Info struct or elsewhere where we could store the
   *        results of this search (or just the need to do the
   *        search).
   */
  bc_node = in_list(I, 0, BC_dup_ptr+1, BC_dup_nodes);
  if (bc_node != -1) {
    /* 
     * current node is in the BC_duplication list.
     * The current boundary condition may have been deleted.
     * If the current boundary condition is a rotated boundary
     * condition, it may have been moved to another coordinate
     * direction.
     */
    i_calc = -1;
    if (ieqn <= LAST_REAL_EQ) {
      i_calc = search_bc_dup_list(bc_input_id,
				  BC_dup_list[bc_node][node_offset]);
    } else {
      EH(-1,"Equation not in list ");
    }

    /*
     * The current boundary condition has been deleted from application
     * at the current node due to a conflict with another boundary
     * condition. Return a -1.
     */
    if (i_calc == -1) {
      return -1;
    }
  }

  /*
   * If the regular unknown equation is replaced by a Dirichlet
   * condition at this node, return -1
   */
  if (node->DBC) {
    if (((int) node->DBC[node_offset]) != -1) {
      return -1;
    }
  } 

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown + node_offset;
  
  /*
   * HKM -> we can put this section into
   *        a subroutine
   */
  if (vd_retn) {
    *vd_retn = NULL;
    num = nv->Num_Var_Desc_Per_Type[ieqn];
    for (i = 0; i < num; i++) {
      index = nv->Var_Type_Index[ieqn][i];
      vd = nv->Var_Desc_List[index];
      if (matID == vd->MatID ||
	  (! (*vd_retn) && vd->MatID == -1)) {
	*vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index_stress                                      */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
evaluate_time_func(const double current_time,
		   double *f_time,      /* computed time function */
		   const int bc_input_id)
    
     /************************************************************************
      *
      * Function which multiplies a time function by previously
      * loaded GD conditions.
      ************************************************************************/
{
  double time = current_time;
  int time_function;
  /* Check if max time was specified and reset time if greater than max time */
  if (BC_Types[bc_input_id].BC_Data_Int[4] == GD_TIME_MAX) {
    if (time > BC_Types[bc_input_id].BC_Data_Float[2]) {
      time = BC_Types[bc_input_id].BC_Data_Float[2];
    }
  }

  /* ---- Find variable number and species number */
  time_function = BC_Types[bc_input_id].BC_Data_Int[2];
  
  /* Now compute time function */
  
  switch(time_function)
    {
    case(GD_TIME_LIN):  /* c0 +c1*t */
      
	*f_time = BC_Types[bc_input_id].BC_Data_Float[0]
	    + BC_Types[bc_input_id].BC_Data_Float[1] * time;
      
      break;
      
    case(GD_TIME_EXP):  /* exp(c0 +c1*t) */
      
	*f_time = exp(BC_Types[bc_input_id].BC_Data_Float[0]
		      + BC_Types[bc_input_id].BC_Data_Float[1] * time);
      
      break;
                  
    case(GD_TIME_SIN):  /* sin(c0 +c1*t) */
      
	*f_time = sin(BC_Types[bc_input_id].BC_Data_Float[0]
		      + BC_Types[bc_input_id].BC_Data_Float[1] * time);
      
      break;
    case (GD_TIME_TABLE):
      {
	double slope, _time[1];
	_time[0] = time;
	*f_time = interpolate_table(BC_Types[bc_input_id].table,
                                    _time, &slope, NULL);
      }
      break;
    default:
      return(-1);
    }
  
  return(0);
} /* END of routine evaluate_time_function                                   */
/*****************************************************************************/

void
apply_table_bc( double *func,
		    double d_func[MAX_VARIABLE_TYPES+MAX_CONC],
		    struct Boundary_Condition *BC_Type,
		    double time_value )
{
  /* 
        Compute the difference between the current value of the 
	field in question and its value obtained by interpolating
	a table of data points.  The abscissa of the table data is
	restricted to one of the three coordinates.  Return the
	the difference and its sensitivity.

	Author : Thomas A. Baer, Org 9111
	Date   : July 16, 1998

    Parameters:
        func = pointer to double that carries back the residual value
	d_func = array of sensitivities of residual wrt to all variables.
	BC_Type = pointer to Boundary Condition structure, i.e. &(BC_Types[bc_input_id]

   */

  int var, basis;
  double slope, interp_val,x_table[2];
  double dfunc_dx[3];

  if(af->Assemble_LSA_Mass_Matrix)
    return;

  basis = BC_Type->table->t_index[0];

  /*  Setup Dummy variable to pass array to interpolation function */
  if( basis != -1 )
    x_table[0] = fv->x[basis];
  else
    x_table[0] = time_value;

  if (BC_Type->table->interp_method == BIQUADRATIC
      || BC_Type->table->interp_method == BILINEAR)
      {
      x_table[1] = fv->x[BC_Type->table->t_index[1]];
      }

  interp_val = interpolate_table( BC_Type->table, x_table, &slope, dfunc_dx );
  interp_val *= BC_Type->BC_Data_Float[0];
  slope *= BC_Type->BC_Data_Float[0];
  
  var = BC_Type->table->f_index ;

  switch ( var ) 
    {
    case VELOCITY1:
      *func = fv->v[0] - interp_val;
      d_func[var] = 1.0;
      break;
    case VELOCITY2:
      *func = fv->v[1] - interp_val;
      d_func[var] = 1.0;
      break;
    case VELOCITY3:
      *func = fv->v[2] - interp_val;
      d_func[var] = 1.0;
      break;
    case TEMPERATURE:
      *func = fv->T - interp_val;
      d_func[var] = 1.0;
      break;
    case MESH_DISPLACEMENT1:
      *func = fv->d[0] - interp_val;
      d_func[var] = 1.0;
      break;
    case MESH_DISPLACEMENT2:
      *func = fv->d[1] - interp_val;
      d_func[var] = 1.0;
      break;
    case MESH_DISPLACEMENT3:
      *func = fv->d[2] - interp_val;
      d_func[var] = 1.0;
      break;
    case MESH_POSITION1:
      *func = fv->x[0] - interp_val;
      d_func[MESH_DISPLACEMENT1] = 1.0;
      break;
    case MESH_POSITION2:
      *func = fv->x[1] - interp_val;
      d_func[MESH_DISPLACEMENT2] = 1.0;
      break;
    case MESH_POSITION3:
      *func = fv->x[2] - interp_val;
      d_func[MESH_DISPLACEMENT3] = 1.0;
      break;
    case MASS_FRACTION:
      *func = fv->c[BC_Type->species_eq] - interp_val;
      d_func[MAX_VARIABLE_TYPES + BC_Type->species_eq] = 1.0;
      break;
    case PRESSURE:
      *func = fv->P - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11:
      *func = fv->S[0][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12:
      *func = fv->S[0][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22:
      *func = fv->S[0][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13:
      *func = fv->S[0][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23:
      *func = fv->S[0][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33:
      *func = fv->S[0][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_1:
      *func = fv->S[1][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_1:
      *func = fv->S[1][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_1:
      *func = fv->S[1][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_1:
      *func = fv->S[1][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_1:
      *func = fv->S[1][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_1:
      *func = fv->S[1][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_2:
      *func = fv->S[2][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_2:
      *func = fv->S[2][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_2:
      *func = fv->S[2][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_2:
      *func = fv->S[2][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_2:
      *func = fv->S[2][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_2:
      *func = fv->S[2][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_3:
      *func = fv->S[3][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_3:
      *func = fv->S[3][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_3:
      *func = fv->S[3][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_3:
      *func = fv->S[3][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_3:
      *func = fv->S[3][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_3:
      *func = fv->S[3][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_4:
      *func = fv->S[4][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_4:
      *func = fv->S[4][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_4:
      *func = fv->S[4][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_4:
      *func = fv->S[4][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_4:
      *func = fv->S[4][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_4:
      *func = fv->S[4][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_5:
      *func = fv->S[5][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_5:
      *func = fv->S[5][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_5:
      *func = fv->S[5][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_5:
      *func = fv->S[5][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_5:
      *func = fv->S[5][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_5:
      *func = fv->S[5][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_6:
      *func = fv->S[6][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_6:
      *func = fv->S[6][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_6:
      *func = fv->S[6][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_6:
      *func = fv->S[6][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_6:
      *func = fv->S[6][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_6:
      *func = fv->S[6][2][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS11_7:
      *func = fv->S[7][0][0] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS12_7:
      *func = fv->S[7][0][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS22_7:
      *func = fv->S[7][1][1] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS13_7:
      *func = fv->S[7][0][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS23_7:
      *func = fv->S[7][1][2] - interp_val;
      d_func[var] = 1.0;
      break;
    case POLYMER_STRESS33_7:
      *func = fv->S[7][2][2] - interp_val;
      d_func[var] = 1.0;
      break;

    default:
      EH(-1,"Variable not yet implemented in TABLE_BC");
      break;
      
    }

  /* And, at the last, we account for the dependence of the interpolation
   * on the basis coordinate
   */

  if (BC_Type->table->interp_method == BIQUADRATIC
              || BC_Type->table->interp_method == BILINEAR)
      {
                      d_func[R_MESH1 + BC_Type->table->t_index[0]]
                      -= dfunc_dx[0]*BC_Type->BC_Data_Float[0];
                      d_func[R_MESH1 + BC_Type->table->t_index[1]]
                      -= dfunc_dx[1]*BC_Type->BC_Data_Float[0];
      }
  else
      {
      if(  basis != -1 && pd->e[R_MESH1 + basis] )
         {
                      d_func[R_MESH1 + basis ] -= slope;
         }
      }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
