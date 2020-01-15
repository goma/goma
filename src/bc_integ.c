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
 
/* Routines for calculating Boundary Conditions and Adding them
 * to matrix_fill */

/*
 *$Id: bc_integ.c,v 5.24 2010-04-07 22:27:00 prschun Exp $
 */

/* Standard include files */
 
#include <stdio.h>
#include <string.h>
 
/* GOMA include files */
 
#include "std.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io.h"
#include "rf_bc_const.h"
#include "mm_elem_block_structs.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_fill_jac.h"
#include "mm_interface.h"
#include "el_elm_info.h"
#include "mm_fill_aux.h"
#include "mm_fill_species.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_fill_potential.h"
#include "mm_shell_bc.h"
#include "rotate_util.h"
#include "mm_eh.h"
#include "user_bc.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc_colloc.h"
#include "bc_integ.h"
#include "exo_struct.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_ns_bc.h"
#include "rd_mesh.h"
#include "rf_bc.h"
#include "rf_solver.h"


#define GOMA_BC_INTEG_C

int
apply_integrated_bc(double x[],           /* Solution vector for the current processor    */
		    double resid_vector[],/* Residual vector for the current processor    */
		    const double delta_t, /* current time step size                       */
		    const double theta,	/* parameter (0 to 1) to vary time integration
					 *  ( implicit - 0 to explicit - 1)             */
		    const PG_DATA *pg_data,
		    const int ielem,       /* element number */
		    const int ielem_type,  /* element type */
		    const int num_local_nodes,
		    const int ielem_dim,
		    const int iconnect_ptr,
		    ELEM_SIDE_BC_STRUCT *elem_side_bc, /* Pointer to an element side boundary condition
							* structure */
		    const int num_total_nodes,
		    const int bc_application, /* flag indicating whether to integrate
					       * strong or weak BC's */
		    const double time_value,
		    SGRID *grid,	
		    const Exo_DB *exo)

/****************************************************************************
 *
 * apply_integrated_bc():
 *
 *    Calculate the local element contributions to boundary conditions which
 *    involve integrations along element edges.
 *
 ****************************************************************************/
{
  int ip, w, i, I, ibc, k, j, id, icount, ss_index, is_ns, mn, lnn;
  int iapply, matID_apply, id_side, i_basis = -1, skip_other_side;
  int new_way = FALSE, ledof, mn_first;
  int eqn, ieqn, var, pvar, p, q, index_eq, ldof_eqn, lvdesc, jlv;
  int err, status = 0;
  int bc_input_id, ip_total; 
  int contact_flag = FALSE;
  int imode;
  int stress_bc = 0;
  double v_attach; 
  double phi_i, tmp;
  double *phi_ptr, *jac_ptr;
  double s, t, u;	      	/* Gaussian quadrature point locations  */
  double xi[DIM];             /* Local element coordinates of Gauss point. */
  double x_dot[MAX_PDIM];
  double x_rs_dot[MAX_PDIM];
  double wt, weight, pb;
  double xsurf[MAX_PDIM];
  double dsigma_dx[DIM][MDE];
  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double func_stress[MAX_MODES][6];
  double d_func_stress[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double cfunc[MDE][DIM];
  double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double time_intermediate = time_value-theta*delta_t; /* time at which bc's are
							  evaluated */
  static INTERFACE_SOURCE_STRUCT *is = NULL;
  static JACOBIAN_VAR_DESC_STRUCT jacCol;
  BOUNDARY_CONDITION_STRUCT *bc;
  MATRL_PROP_STRUCT *mp_2;
  struct BC_descriptions *bc_desc;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  double surface_centroid[DIM]; 
  int interface_id = -1;

  tran->time_value = time_intermediate;

  /***************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)        */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                      */
  /***************************************************************************/

  /* Find out the number of surface quadrature points 
   *  -this is assumed independent of the type of boundary condition
   *   applied at the surface 
   */

  if ( ls != NULL && ls->SubElemIntegration && ls->elem_overlap_state )
    {
      Subgrid_Int.ip_total = get_subelement_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, 0., elem_side_bc->id_side-1, 0 );
      ip_total = Subgrid_Int.ip_total;
    }
  else if ( ls != NULL && ls->Integration_Depth > 0 && ls->elem_overlap_state )
    {
#ifdef SUBELEMENT_FOR_SUBGRID
      /* DRN: you can also do subgrid integration with the following.
         it recursive divides element to create small subelements and then
         creates subgrid integration points */
      Subgrid_Int.ip_total = get_subelement_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, 0., elem_side_bc->id_side-1, 0 );
      ip_total = Subgrid_Int.ip_total;
#else
      /* Tom - We need subgrid integration along this element side here */
      /* for now just use regular gauss integration */
      find_surf_center_st ( ielem_type, elem_side_bc->id_side, ielem_dim, surface_centroid, &s, &t );
      
      ip_total = gather_surface_subgrid_integration_pts( grid, elem_side_bc->id_side, surface_centroid, Subgrid_Int.s, Subgrid_Int.wt, 0 );  
	  
      /* print_subgrid_surface_integration_pts( Subgrid_Int.s, Subgrid_Int.wt, ip_total);  */

#endif
    }
  else
    {
      ip_total = elem_info(NQUAD_SURF, ielem_type);
    }

  /*
   *  Loop over the quadrature points at the surface
   */
  for (ip = 0; ip < ip_total; ip++) {
    
    if (ls != NULL && ls->SubElemIntegration && ls->elem_overlap_state)
      {
        ls->Elem_Sign = Subgrid_Int.ip_sign[ip];
	xi[0] = Subgrid_Int.s[ip][0];
        xi[1] = Subgrid_Int.s[ip][1];
	xi[2] = Subgrid_Int.s[ip][2];
	/* DRN: Are these really needed? I'll put garbage in them so someone knows not to! */
	s = 1.e30; 
	t = 1.e30;
        wt = Subgrid_Int.wt[ip];
      }
    else if ( ls != NULL && ls->Integration_Depth > 0 && ls->elem_overlap_state )
      {
#ifdef SUBELEMENT_FOR_SUBGRID
        ls->Elem_Sign = 0;
	xi[0] = Subgrid_Int.s[ip][0];
        xi[1] = Subgrid_Int.s[ip][1];
	xi[2] = Subgrid_Int.s[ip][2];
	/* DRN: Are these really needed? I'll put garbage in them so someone knows not to! */
	s = 1.e30; 
	t = 1.e30;
        wt = Subgrid_Int.wt[ip];
#else
        /* Tom - We need subgrid integration along this element side here */
        /* for now just use regular gauss integration */
        if ( ls != NULL ) ls->Elem_Sign = 0;
	xi[0] = Subgrid_Int.s[ip][0];
        xi[1] = Subgrid_Int.s[ip][1];
	xi[2] = Subgrid_Int.s[ip][2];
	/* DRN: Are these really needed? I'll put garbage in them so someone knows not to! */
	s = 1.e30; 
	t = 1.e30;
        wt = Subgrid_Int.wt[ip];
	/* find the quadrature point locations (s, t) for current ip 
	   find_surf_st(ip, ielem_type, elem_side_bc->id_side, pd->Num_Dim, xi, &s, &t, &u);
	   find the quadrature weight for current surface ip 
	   wt = Gq_surf_weight(ip, ielem_type); */
#endif
      }
    else
      {
        if ( ls != NULL ) ls->Elem_Sign = 0;
	/* find the quadrature point locations (s, t) for current ip */
	find_surf_st(ip, ielem_type, elem_side_bc->id_side,
                     ei[pg->imtrx]->ielem_dim, xi, &s, &t, &u);
        /* find the quadrature weight for current surface ip */
        wt = Gq_surf_weight(ip, ielem_type);
      } 
    
    err = load_basis_functions(xi, bfd);
    EH(err, "problem from load_basis_functions");
    
    err = beer_belly();
    EH(err, "beer_belly");
    
    /*
     *  precalculate variables at current integration pt.
     *  for the current material comprising the current element
     */
    err = load_fv();
    EH(err, "load_fv");

    /* What's going on here */

    err = load_bf_grad();
    EH(err, "load_bf_grad");

    err = load_bf_mesh_derivs(); 
    EH(err, "load_bf_mesh_derivs");

    /* calculate the determinant of the surface jacobian and the normal to 
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
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    err = load_fv_grads();
    EH( err, "load_fv_grads");
    
    err = load_fv_mesh_derivs(1);
    EH(err, "load_fv_mesh_derivs");
    
    /*
     * Load up commonly used physical properties such as density at
     * the current quadrature point using the material state vector. 
     */
    load_properties(mp, time_value);
   
    /*
     * Determine the State Variable Vector for the material
     * on the "other" side of the interface, i.e., the material
     * located at the interface that is not the current material.
     * Then, calculate common properties using that state vector.
     * At the end of this loop, mp_2 should be pointing to the
     * material on the other side of the interface, and its
     * constituitive properties should be filled in.
     */
    for (i = 0; i < elem_side_bc->Num_MatID; i++) {
      mp_2 = mp_glob[elem_side_bc->MatID_List[i]];
      if (mp_2 != mp) {
        load_matrl_statevector(mp_2);
	load_properties(mp_2, time_value);
	break;
      }
    }    

    /*
     * Load up porous media variables and properties, if needed 
     */
    if (mp->PorousMediaType == POROUS_UNSATURATED ||
	mp->PorousMediaType == POROUS_SATURATED || 
	mp->PorousMediaType == POROUS_TWO_PHASE) {
      err = load_porous_properties(); 
      EH( err, "load_porous_properties"); 
    }

    if (mp->SurfaceTensionModel != CONSTANT) {
      load_surface_tension(dsigma_dx);
      if (neg_elem_volume) return(status);
    }
    
    if (TimeIntegration != STEADY && pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (icount = 0; icount < ielem_dim; icount++ ) {
	x_dot[icount] = fv_dot->x[icount];
	/* calculate surface position for wall repulsion/no penetration condition */
	xsurf[icount] = fv->x0[icount];
      }

      if (pd->e[pg->imtrx][SOLID_DISPLACEMENT1]) {
	for (icount = 0; icount < VIM; icount++)
	  x_rs_dot[icount] = 0.;
	for (icount = 0; icount < ielem_dim; icount++ )	{
	  /* Notice how I use here fv->d_rs instead of fv->x_rs
	   * (which doesn't exist) because
	   * x_rs_dot = d(X_rs)/dt = d(Coor[][] + d_rs)/dt = d(d_rs)/dt
	   */
	  x_rs_dot[icount] = fv_dot->d_rs[icount];
	}
      }
    }  else {
      for (icount = 0; icount < ielem_dim; icount++ ) {
	x_rs_dot[icount] = 0.;
	x_dot[icount] = 0.;
	xsurf[icount] = fv->x0[icount];
      }
    }
   
    /*
     *  Loop over all of the boundary conditions assigned to this side
     *  of the element
     */
    do_LSA_mods(LSA_SURFACE);

    for (ibc =0;
	 (bc_input_id = (int) elem_side_bc->BC_input_id[ibc]) != -1;
	 ibc++) {
      /*
       *  Create a couple of pointers to cut down on the
       *  amount of indirect addressing
       */
      bc = BC_Types + bc_input_id;
      bc_desc = bc->desc;

      ss_index = in_list(bc->BC_ID, 0, exo->num_side_sets,
                         &(ss_to_blks[0][0]));
      is_ns = strcmp(BC_Types[bc_input_id].Set_Type, "NS");
      if (ss_index == -1 && is_ns != 0) {
	sprintf(Err_Msg, "Could not find BC_ID %d in ss_to_blks",
	        BC_Types[bc_input_id].BC_ID);
	EH(-1, Err_Msg);
      }

      /*
       *  Set flag to indicate if we're in the right
       *  material (only one) to apply
       */
      if (bc->BC_Name == VL_EQUIL_PRXN_BC ||
	  bc->BC_Name == IS_EQUIL_PRXN_BC ||
	  bc->BC_Name == YFLUX_DISC_RXN_BC ||
	  bc->BC_Name == SDC_STEFANFLOW_BC ||
	  bc->BC_Name == SDC_KIN_SF_BC ) {
	  new_way = TRUE;
	  for (icount = 0; icount < Num_Interface_Srcs; icount++)  {
                if(bc->BC_ID == IntSrc_BCID[icount])    {
                     interface_id = icount;
                     }
                }
          if (is)
            {
	     for (icount = 0; icount < mp->Num_Species; icount++)  {
                    is[interface_id].Processed[icount] = FALSE;
                    }
	    }
      } else {
        new_way = FALSE;
      }
      iapply = 0;
      skip_other_side = FALSE;
      if (is_ns != 0) {
        if (ei[pg->imtrx]->elem_blk_id == ss_to_blks[1][ss_index]) {
          iapply = 1;
        }
      }

      /*
       *  However, override if the side set is an external one. In
       *  other words, if the side set is an external side set, we will
       *  apply the boundary condition no matter what.
       */
      if (SS_Internal_Boundary != NULL ) {
        if (SS_Internal_Boundary[ss_index] == -1) {
          iapply = 1;
        }
      }
	
      /*  check to see if this bc is an integrated bc and thus to be handled
       *  by this routine and not others. Also, this routine is called
       *  twice by matrix fill, once for weakly integrated bc's and once for
       *  strongly integrated bc's. Only go forward with the pertinent bc
       *  from here.
       */
      if (bc_desc->method == bc_application) {

	/*
	 * Check to see if this is a new wave boundary condition and thus
	 * to use Jacbobian_Var_Desc structure to hold jacobian terms
	 */
	if (bc->BC_Name > NEW_WAY_LOW_BC &&
	    bc->BC_Name < NEW_WAY_HIGH_BC) {
          new_way = TRUE;
	}

	/* Initialize the general function to zero may have more than one entry
	 * for vector conditions like capillary. 
	 * If not the new way, first do a
	 * really long and expense zeroing of a vast amount of basically unused
	 * array locations.
	 */
	if (!new_way) {
	  func[0] = 0.0; func[1] =0.0; func[2] = 0.0;
		
	  if ( bc_desc->vector == VECTOR )  /* save some initialization if not a vector cond */
	    memset(d_func, 0, DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	  else
	    memset(d_func,0, (MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));

          memset(func_stress, 0.0, MAX_MODES * 6 * sizeof(double));
          memset(d_func_stress, 0.0, MAX_MODES * 6 * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
	}
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
	 * Step 4 should be done below (or in bc_colloc.c), and add your new function at 
	 *     the end of this file (see fplane)
	 */

	switch (bc->BC_Name) {
	    
	case KINEMATIC_PETROV_BC:
	case KINEMATIC_BC:
	case VELO_NORMAL_BC:
	case VELO_NORMAL_LS_BC:
	case VELO_NORMAL_LS_PETROV_BC:
		
	  contact_flag = (ls != NULL);

	  /*  first all external boundaries with velocity
	      second - internal boundaries with an explicit block id
	      third  - internal boundaries with implicit iapply logic
	  */
	  if ( (SS_Internal_Boundary[ss_index] == -1 && pd->v[pg->imtrx][VELOCITY1])
	       || (SS_Internal_Boundary[ss_index] != -1 &&
		   bc->BC_Data_Int[0] == ei[pg->imtrx]->elem_blk_id)
	       || (SS_Internal_Boundary[ss_index] != -1 && 
		   bc->BC_Data_Int[0] == -1 && iapply && pd->v[pg->imtrx][VELOCITY1]))
	    {
	      fvelo_normal_bc(func, d_func, bc->BC_Data_Float[0], contact_flag,
			      x_dot, theta, delta_t, (int) bc->BC_Name,
                                bc->BC_Data_Float[1], bc->BC_Data_Float[2],
                                bc->BC_Data_Float[3]);
	    }
	  break;

	case VELO_TANGENT_LS_BC:
		
	  /*  first all external boundaries with velocity
	      second - internal boundaries with an explicit block id
	      third  - internal boundaries with implicit iapply logic
	  */
          if (ls != NULL)
          {
	  if ( (SS_Internal_Boundary[ss_index] == -1 && pd->v[pg->imtrx][VELOCITY1])
	       || (SS_Internal_Boundary[ss_index] != -1 &&
		   bc->BC_Data_Int[0] == ei[pg->imtrx]->elem_blk_id)
	       || (SS_Internal_Boundary[ss_index] != -1 && 
		   bc->BC_Data_Int[0] == -1 && iapply && pd->v[pg->imtrx][VELOCITY1]))
	    {
	      fvelo_tangential_ls_bc(func, 
                                     d_func, 
                                     bc->BC_Data_Float[0], 
			             x_dot, 
                                     theta, 
                                     delta_t, 
                                     bc->BC_Data_Float[1], 
                                     bc->BC_Data_Float[2],
                                     bc->BC_Data_Float[3]);
	    }
          }
	  break;

	case LS_ATTACH_BC:
	  fvelo_normal_bc(func, d_func, 0.0, contact_flag = FALSE, 
			  x_dot, theta, delta_t, (int) bc->BC_Name,0,0,135.0);
	  ls_attach_bc( func, d_func, bc->BC_Data_Float[0] );
	  break;
		
        case LS_WALL_ANGLE_BC:
          ls_wall_angle_bc(func, d_func, bc->BC_Data_Float[0]);
          break;

 	case KIN_DISPLACEMENT_PETROV_BC:
	case KIN_DISPLACEMENT_BC:
	  f_kinematic_displacement_bc(func, d_func, bc->BC_Data_Int[0], bc->BC_ID, bc->u_BC, bc->len_u_BC);
	  break;

	case KIN_DISPLACEMENT_RS_BC:
	  f_kinematic_displacement_rs_bc(func, d_func, bc->BC_Data_Int[0], 
					 elem_side_bc->id_side, xi, exo );
	  break;

	case KINEMATIC_DISC_BC:
	case VELO_NORMAL_DISC_BC:
	  fvelo_normal_disc_bc(func, d_func, bc->BC_Data_Float[0], 
			       x_dot, theta, delta_t);
	  break;

	case KINEMATIC_SPECIES_BC:
	  kinematic_species_bc(func, d_func, bc->BC_Data_Int[0], 
			       bc->BC_Data_Float[0], x_dot, theta, delta_t);
	  break;

	case CONT_TANG_VEL_BC:
	  continuous_tangent_velocity(func, d_func, ielem_dim);
	  break;

	case CONT_NORM_VEL_BC:
	  continuous_normal_velocity(func, d_func, ielem_dim);
	  break;

	case DISCONTINUOUS_VELO_BC:
	  discontinuous_velocity(func, d_func, x_dot, bc->BC_Data_Int[0],
				 bc->BC_Data_Int[1], bc->BC_Data_Int[2],
				 theta, delta_t);
	  break;

	case DARCY_CONTINUOUS_BC:
	case DARCY_LUB_BC:
	  
	  if(time_intermediate >= bc->BC_Data_Float[1])
	    {
	      v_attach = -0.0;  /*This will eventually be replaced by 
				 *a squeeze flow model */  

	      sat_darcy_continuous_bc(func, d_func, theta, delta_t, time_value,
				      bc->BC_Data_Int[0], bc->BC_Data_Int[1],
				      bc->BC_Data_Float[0], v_attach);
	      if (neg_elem_volume) return (status);
	    }

	  break;
	    
	case T_MELT_BC:
	  if (iapply) {
	    ftmelt_bc(func, d_func, bc->BC_Data_Float[0], 
		      x_dot, theta, delta_t);
	  }
	  break;
		
	case KIN_LEAK_BC:
	case KIN_LEAK_HEAT_BC:
	case VNORM_LEAK_BC:
        case KIN_CHEM_BC:
	  if(iapply) {
	      kin_bc_leak(func, d_func, x_dot, theta, delta_t,
			  bc_input_id, BC_Types);
	  }
	    break;

        case SHELL_SURFACE_CHARGE_BC:  /* Applies only to shell elements */
	  shell_surface_charge_bc(func, d_func, x_dot, theta, delta_t,
				  elem_side_bc->id_side, wt, xi, exo, 0);
	  break;
        case SHELL_SURFACE_CHARGE_SIC_BC:  /* Applies only to shell elements */
	  shell_surface_charge_bc(func, d_func, x_dot, theta, delta_t,
				  elem_side_bc->id_side, wt, xi, exo, 1);
	  break;

        case SHELL_DIFF_KINEMATIC_BC:  /* Applies only to shell elements */
	  shell_diff_kinematic_bc(func, d_func, x_dot, theta, delta_t,
				  elem_side_bc->id_side, wt, xi, exo, 0);
	  break;

        case SH_LUBP_SOLID_BC:
	case SH_LUBP_SOLID_RS_BC:
	  shell_lubr_solid_struct_bc(func, d_func, x_dot, theta, delta_t, 
				   elem_side_bc->id_side, wt, xi, exo, 
				   bc->BC_Data_Float[0]);

	  break;

	case LUBP_SH_FP_MATCH_BC:
	  match_lubrication_film_pressure(func, d_func,
					  bc->BC_Data_Int[0], 
					  bc->BC_Data_Int[1]);
	  break;

        case KIN_ELECTRODEPOSITION_BC:  /*  RSL 5/28/02  */
	  if (iapply) {
	    kin_bc_electrodeposition(func, d_func, time_value, theta, delta_t, bc_input_id, BC_Types);
	  }
	  break;

	case VNORM_ELECTRODEPOSITION_BC:
	  if (iapply) {
	    vnorm_bc_electrodeposition(func, d_func, theta, time_value, delta_t, bc_input_id, BC_Types);
	  }
	  break;

	case VELO_TANGENT_BC:
	case VELO_TANGENT_USER_BC:
 	case VELO_STREAMING_BC:
	  fvelo_tangential_bc(func, d_func, x,
			      bc->BC_Data_Int[0], 
			      bc->BC_Data_Float[0], 
			      bc->BC_Data_Float[1], 
			      bc->BC_Data_Float[2], 
			      xsurf, x_dot, theta, delta_t, (int) bc->BC_Name,
			      elem_side_bc->id_side, xi, exo,
			      time_intermediate, bc->u_BC, bc->len_u_BC);
	  if (neg_elem_volume) return (status);
	  break;

	case VELO_TANGENT_3D_BC:
	  fvelo_tangent_3d(func, d_func, x_dot, 
			   bc->BC_Data_Float[0], 
			   bc->BC_Data_Float[1], 
			   bc->BC_Data_Float[2], 
			   bc->BC_Data_Float[3], 
			   theta, delta_t);
	  if (neg_elem_volume) return (status);
	  break;
        case ZERO_VELO_TANGENT_3D_BC:
          fzero_velo_tangent_3d(func, d_func, elem_side_bc->id_side, 0);
          break;

	case VELO_TANGENT_SOLID_BC:
	case VELO_SLIP_SOLID_BC:
	  fvelo_tangential_solid_bc(func, d_func, x,  x_dot, x_rs_dot,
				    (int) bc->BC_Name,
				    bc->BC_Data_Int[0], 
				    bc->BC_Data_Int[1], 
				    bc->BC_Data_Float[0],
				    bc->BC_Data_Int[2],
				    bc->BC_Data_Float[1],
				    xsurf, theta, delta_t);
	  break;

	case VELO_NORMAL_SOLID_BC:

	  fvelo_normal_solid_bc(func, d_func, x,  x_dot, x_rs_dot,
				(int) bc->BC_Name,
				bc->BC_Data_Int[0], 
				bc->BC_Data_Int[1], 
				theta, delta_t);
	  break;

	case VELO_SLIP_BC:
	case VELO_SLIP_ROT_BC:
	case VELO_SLIP_FILL_BC:
 	case VELO_SLIP_ROT_FILL_BC:
	case VELO_SLIP_FLUID_BC:
	case VELO_SLIP_ROT_FLUID_BC:
	  fvelo_slip_bc(func, d_func, x, 
			(int) bc->BC_Name,
			(int) bc->max_DFlt,
			bc->BC_Data_Float,
			(int) bc->BC_Data_Int[0],
			xsurf, theta, delta_t);
	  break;
	case VELO_SLIP_POWER_CARD_BC:
	case VELO_SLIP_POWER_BC:
	  fvelo_slip_power_bc(func, d_func,
			      (int) bc->BC_Name,
			      (int) bc->max_DFlt,
			      bc->BC_Data_Float,
			      theta, delta_t);
	break;

	case AIR_FILM_BC:
	case AIR_FILM_ROT_BC:
	  fvelo_airfilm_bc(func, d_func, x, 
			(int) bc->BC_Name,
			bc->BC_Data_Float,
			(int) bc->BC_Data_Int[0],
			xsurf, theta, delta_t);
	  break;  

	case VELO_SLIP_EK_BC:
	  fvelo_slip_electrokinetic_bc(func, d_func, 
				       bc->BC_Data_Float[0],
				       bc->BC_Data_Float[1]);
	  break;  

	case VELO_EK_3D_BC:
	  fvelo_electrokinetic_3d(func, d_func, x_dot, 
				  bc->BC_Data_Float[0], 
				  bc->BC_Data_Float[1], 
				  bc->BC_Data_Float[2], 
				  bc->BC_Data_Float[3], 
				  bc->BC_Data_Float[4], 
				  theta, delta_t);
	  if (neg_elem_volume) return (status);
	  break;

	case VELO_SLIP_LEVEL_BC:
	case VELO_SLIP_LEVEL_SIC_BC:
	case VELO_SLIP_LS_ROT_BC:
	  fvelo_slip_level( func, d_func, 
		             (int) bc->BC_Name,
			     bc->BC_Data_Float[0], 
			     bc->BC_Data_Float[1], 
			     bc->BC_Data_Float[2], 
			     bc->BC_Data_Float[3], 
			     bc->BC_Data_Float[4],
			     bc->BC_Data_Float[5],
			     bc->BC_Data_Float[6],
			     bc->BC_Data_Float[7],
			     bc->BC_Data_Float[8],
                             theta, delta_t);
	  break;

	case VELO_SLIP_LS_HEAVISIDE_BC:
	  fvelo_slip_ls_heaviside( func, d_func,
				   bc->BC_Data_Float[0],
				   bc->BC_Data_Float[1],
				   bc->BC_Data_Float[2],
				   bc->BC_Data_Float[3],
				   bc->BC_Data_Float[4],
				   bc->BC_Data_Float[5],
				   theta, delta_t);
	  break;


	case Q_VELO_SLIP_BC:
	  q_velo_slip_bc(func, d_func, 
			 (int) BC_Types[bc_input_id].BC_Data_Int[0],
			 x, xsurf, theta, delta_t);
	  break;

	case NO_SLIP_BC:
	case NO_SLIP_RS_BC:    /*for backward compatibility even though
				 NO_SLIP_RS was removed on 12/17/98 */
	  no_slip(func, d_func, x_dot, x_rs_dot,
		  theta, delta_t, bc_input_id, BC_Types,
		  (int) bc->BC_Data_Int[0],
		  (int) bc->BC_Data_Int[1]);
	  break;

	case VN_POROUS_BC:
	  porous_normal_velocity_bc(func, d_func, x_dot, theta, delta_t, 
				    bc_input_id, BC_Types,
				    (int) bc->BC_Data_Int[0],
				    (int) bc->BC_Data_Int[1],
				    (int) bc->BC_Data_Int[2],
				    bc->BC_Data_Float[0]);
	  break;

	case UUSER_BC:
	  uuser_surf(func, d_func, bc->u_BC, time_intermediate);
	  break;

	case VUSER_BC:
	  vuser_surf(func, d_func, bc->u_BC, time_intermediate);
	  break;

	case WUSER_BC:
	  wuser_surf(func, d_func, bc->u_BC, time_intermediate);
	  break;

	case SLOPEX_BC:
	case SLOPEY_BC:
	case SLOPEZ_BC:
	case SLOPE_BC:
	  slope_n_dot_n0_bc(func, d_func,
			    bc->BC_Data_Float[0],
			    bc->BC_Data_Float[1],
			    bc->BC_Data_Float[2]);
	  break;
	      
	case FORCE_BC:
	case FORCE_RS_BC:
	  force_n_dot_f_bc(func, d_func,
			   bc->BC_Data_Float[0],
			   bc->BC_Data_Float[1],
			   bc->BC_Data_Float[2], 0,
			   delta_t, theta, ip, ip_total, time_intermediate);
	  break;
	case FORCE_SIC_BC:
	  force_n_dot_f_bc(func, d_func,
			   bc->BC_Data_Float[0],
			   bc->BC_Data_Float[1],
			   bc->BC_Data_Float[2], 1,
			   delta_t, theta, ip, ip_total, time_intermediate);
	  break;

	case FORCE_USER_BC:
	case FORCE_USER_RS_BC:
	  force_user_surf (func, d_func,
			   bc->u_BC,
			   time_intermediate);
	  break;

	case FORCE_USER_SIC_BC:
	  force_user_surf (func, d_func,
			   bc->u_BC,
			   time_intermediate);
	  force_n_dot_f_bc(func, d_func,
			   0, 0, 0, 2,
			   delta_t, theta, ip, ip_total, time_intermediate);
	  break;

	case REP_FORCE_BC:
	case REP_FORCE_RS_BC:
        case ATTR_FORCE_BC:
        case ATTR_FORCE_RS_BC:
	  rep_force_n_dot_f_bc(func, d_func,
			       bc->BC_Data_Float[0],
			       bc->BC_Data_Float[1],
			       bc->BC_Data_Float[2],
			       bc->BC_Data_Float[3],
			       bc->BC_Data_Float[4],
			       bc->BC_Data_Float[5],
			       bc->BC_Data_Float[6],
			       (int)bc->BC_Name );
	  break;
	case REP_FORCE_ROLL_BC:
	case REP_FORCE_ROLL_RS_BC:
	  rep_force_roll_n_dot_f_bc(func, d_func,
			       bc->BC_Data_Float[0],
			       &(bc->BC_Data_Float[1]),
			       &(bc->BC_Data_Float[4]),
			       bc->BC_Data_Float[7],
			       bc->BC_Data_Float[8],
			       bc->BC_Data_Float[9],
			       (int)bc->BC_Name );
	  break;

	case NORM_FORCE_BC:
	case NORM_FORCE_RS_BC:
	  norm_force_n_dot_f_bc(func, d_func,
				bc->BC_Data_Float[0],
				bc->BC_Data_Float[1],
				bc->BC_Data_Float[2]);
	  break;
	      
	case FRICTION_BC:
	case FRICTION_RS_BC:
	case FRICTION_ACOUSTIC_BC:
	case FRICTION_ACOUSTIC_RS_BC:
	  friction_n_dot_f_bc(func, d_func,
			      bc->BC_Data_Float[0],
			      bc->BC_Data_Int[0],
			      delta_t, theta, ip, ip_total, (int) bc->BC_Name,
			      time_intermediate, bc->u_BC, bc->len_u_BC);
	  break;

        case REP_FORCE_SHU_BC:  /* Applies only to shell elements */
	  rep_force_shell_n_dot_f_bc(func, d_func, x_dot, theta, delta_t,
				     ip, ip_total, elem_side_bc->id_side, wt, xi, exo, 0);
	  break;
        case REP_FORCE_SHU_SIC_BC:  /* Applies only to shell elements */
	  rep_force_shell_n_dot_f_bc(func, d_func, x_dot, theta, delta_t,
				     ip, ip_total, elem_side_bc->id_side, wt, xi, exo, 1);
	  break;


	case ELEC_TRACTION_BC:
	case ELEC_TRACTION_SOLID_BC:

	  /* There are no time derivatives in here. */
	  if (af->Assemble_LSA_Mass_Matrix) break;

	  memset(cfunc, 0, MDE*DIM*sizeof(double));
	  memset(d_cfunc, 0,
		 MDE*DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));

	  /*
	   * Electric Stress
	   *
	   * Add the electric stress if we're in the specified material
	   * (first integer parameter for the BC).
	   */

	  if (bc->BC_Data_Int[0] == ei[pg->imtrx]->elem_blk_id)
	    {
	      elec_surf_stress(cfunc,
			       d_cfunc,
			       elem_side_bc->id_side,
			       elem_side_bc,
			       iconnect_ptr,
  			       mp->permittivity,
			       bc->BC_Data_Float[0],
			       (int) bc->BC_Name);
	      
	    }
	  else
	    {
	      skip_other_side = TRUE;
	    }
	  break;
	  
	case SHEAR_TO_SHELL_BC:
	  memset(cfunc, 0, MDE*DIM*sizeof(double));
	  memset(d_cfunc, 0,MDE*DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
		
	  shear_to_shell ( cfunc, d_cfunc, elem_side_bc->id_side,	
			   bc->BC_Data_Float[0],
			   elem_side_bc,
			   iconnect_ptr,
			   xi,
			   exo);
	  break;
	  
	case TENSION_SHEET_BC:
	  memset(cfunc, 0, MDE*DIM*sizeof(double));
	  memset(d_cfunc, 0, MDE*DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));
	  sheet_tension( cfunc,d_cfunc,
			 elem_side_bc->id_side,
			 bc->BC_Data_Float[0],
			 elem_side_bc,
			 iconnect_ptr,
                         xi, exo);
	  break;

	case CAPILLARY_BC:
	case CAP_REPULSE_BC:
	case CAP_REPULSE_ROLL_BC:
	case CAP_REPULSE_USER_BC:
	case CAP_REPULSE_TABLE_BC:
	case CAP_RECOIL_PRESS_BC:
	case CAPILLARY_SHEAR_VISC_BC:
	case CAPILLARY_TABLE_BC:
	  memset(cfunc, 0, MDE*DIM*sizeof(double));
	  memset(d_cfunc, 0,
		 MDE*DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double));

	  /* MMH: I looked in the fncs below, and there were no time derivatives. */
	  if (af->Assemble_LSA_Mass_Matrix) break;

	  if ( bc->BC_Data_Int[0] == ei[pg->imtrx]->elem_blk_id ||
	       ( bc->BC_Data_Int[0] == -1 && iapply ) ) {

	    /* set up surface repulsion force and zero it if 
	     * bc is CAP_REPULSE or CAP_RECOIL_PRESS
	     * because then force is calculated
	     * later */
	    pb = BC_Types[bc_input_id].BC_Data_Float[1];
	    if (BC_Types[bc_input_id].BC_Name == CAPILLARY_TABLE_BC)
               {
	  apply_table_wic_bc(func, d_func, &BC_Types[bc_input_id], time_value);
              pb = func[0];
               }

	    fn_dot_T(cfunc, d_cfunc, elem_side_bc->id_side,
		     bc->BC_Data_Float[0], pb,
		     elem_side_bc, iconnect_ptr, dsigma_dx);

	    if (bc->BC_Name == CAP_REPULSE_BC) {
	      apply_repulsion(cfunc, d_cfunc, bc->BC_Data_Float[2],
			      bc->BC_Data_Float[3], bc->BC_Data_Float[4],
			      bc->BC_Data_Float[5], bc->BC_Data_Float[6],
			      elem_side_bc, iconnect_ptr);
	    }

	    if (bc->BC_Name == CAP_REPULSE_ROLL_BC) {
	      apply_repulsion_roll(cfunc, d_cfunc, x, 
                               bc->BC_Data_Float[2], 
			       &(bc->BC_Data_Float[3]),
			       &(bc->BC_Data_Float[6]),
                               bc->BC_Data_Float[9], 
                               bc->BC_Data_Float[10], 
                               bc->BC_Data_Float[11], 
                               bc->BC_Data_Float[12], 
                               bc->BC_Data_Float[13], 
                               bc->BC_Data_Float[14], 
                               bc->BC_Data_Int[2],
			      elem_side_bc, iconnect_ptr);
	    }
	    if (bc->BC_Name == CAP_REPULSE_USER_BC) {
	      apply_repulsion_user(cfunc, d_cfunc, 
                               bc->BC_Data_Float[2], 
			       &(bc->BC_Data_Float[3]),
			       &(bc->BC_Data_Float[6]),
                               bc->BC_Data_Float[9], 
                               bc->BC_Data_Float[10], 
                               bc->BC_Data_Float[11], 
                               bc->BC_Data_Float[12], 
                               bc->BC_Data_Float[13], 
			      elem_side_bc, iconnect_ptr);
	    }
	    if (bc->BC_Name == CAP_REPULSE_TABLE_BC) {
	      apply_repulsion_table(cfunc, d_cfunc, x, bc->BC_Data_Float[2], 
                               bc->BC_Data_Float[3], 
                               bc->BC_Data_Float[4], 
                               bc->BC_Data_Float[5], 
                               bc->BC_Data_Float[6], 
			       &(bc->BC_Data_Float[7]),
                               bc->BC_Data_Int[2],
			      elem_side_bc, iconnect_ptr);
	    }
	    if (BC_Types[bc_input_id].BC_Name == CAP_RECOIL_PRESS_BC)
	      {
		apply_vapor_recoil(cfunc, d_cfunc, 
				   BC_Types[bc_input_id].BC_Data_Float[2],
				   BC_Types[bc_input_id].BC_Data_Float[3],
				   BC_Types[bc_input_id].BC_Data_Float[4],
				   BC_Types[bc_input_id].BC_Data_Float[5],
				   BC_Types[bc_input_id].BC_Data_Float[6],
				   elem_side_bc, iconnect_ptr);
	      }
	    if (BC_Types[bc_input_id].BC_Name == CAPILLARY_SHEAR_VISC_BC)
	      {
		apply_surface_viscosity(cfunc, d_cfunc, bc->BC_Data_Float[0], bc->BC_Data_Float[1],
                                        bc->BC_Data_Float[2], bc->BC_Data_Float[3], time_value,
					elem_side_bc, wt, xi, exo, iconnect_ptr);
	      }
	  }
	  break;

	case FLOW_PRESS_USER_BC:
	case PRESSURE_USER_BC:
	  fn_dot_T_user(func, d_func, bc->u_BC,
			time_intermediate);
	  break;

	case FLOW_PRESSURE_BC:
	  flow_n_dot_T_hydro(func, d_func,0.,0.,0.,
			     bc->BC_Data_Float[0]);
	  break;

	case LGR_FLOWRATE_BC:
	  {
	    int iAC = bc->BC_Data_Int[0];  /* need this index so we can retrieve the current LM value from augc */

	    flow_n_dot_T_hydro(func, d_func, 0.,0.,0.,
			       -augc[iAC].DataFlt[0] );
	  }
	  break;


	case FLOW_HYDROSTATIC_BC:
	  flow_n_dot_T_hydro(func, d_func,
			     bc->BC_Data_Float[0],
			     bc->BC_Data_Float[1],
			     bc->BC_Data_Float[2],
			     bc->BC_Data_Float[3]);
	  break;

	case FLOW_PRESSURE_VAR_BC:
	  flow_n_dot_T_var_density(func, d_func,
			     bc->BC_Data_Float[0],
			     time_value);
	  break;

	case FLOW_STRESSNOBC_BC:
	  flow_n_dot_T_nobc(func, d_func,
			    bc->BC_Data_Float[0],
			    bc->BC_Data_Int[0]);
	  break;

	case FLOW_GRADV_BC:
	  flow_n_dot_T_gradv(func, d_func,
			     bc->BC_Data_Float[0],
			     bc->BC_Data_Int[0]);
	  break;


	case FLOW_GRADV_SIC_BC:
	  flow_n_dot_T_gradv_sic(func, d_func,
			     bc->BC_Data_Float[0],
			     bc->BC_Data_Int[0]);
	  break;
	  
        case STRESS_DEVELOPED_BC:
          if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV)
            {
              stress_no_v_dot_gradS_logc(func_stress, d_func_stress, delta_t, theta);
            }
          else
            {
              stress_no_v_dot_gradS(func_stress, d_func_stress, delta_t, theta);
            }
          break;


        case GRAD_LUB_PRESS_BC:
	  shell_n_dot_flow_bc_confined(func, d_func,
                                       bc->BC_Data_Float[0], 
                                       time_value, delta_t,
                                       xi, exo);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;

         case LUB_STATIC_BC:
         lub_static_pressure(func, d_func,
                             bc->BC_Data_Float[0],
                             time_value, delta_t,
                             xi, exo);
         surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                                       	ielem_dim - 1,
                                       	(int) elem_side_bc->id_side,
                                       	(int) elem_side_bc->num_nodes_on_side,
                                       	(elem_side_bc->local_elem_node_id) );
          break;


        case SHELL_GRAD_FP_BC:
	  shell_n_dot_flow_bc_film(func, d_func,
                                   bc->BC_Data_Float[0], 
                                   time_value, delta_t,
                                   xi, exo);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;


        case SHELL_FLOW_DEVELOPED_BC:
	  shell_n_dot_gradp_bc(func, d_func, 
                                time_value, delta_t,
                                xi, exo);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;


	case SHELL_GRAD_FP_NOBC_BC:
	  shell_n_dot_flow_bc_film(func, d_func, 0.0, 
                                   time_value, delta_t,
                                   xi, exo);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;



        case SHELL_GRAD_FH_BC:
	  shell_n_dot_gradh_bc(func, d_func,
                                  bc->BC_Data_Float[0]);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;

	case SHELL_GRAD_FH_NOBC_BC:
	  shell_n_dot_gradh_bc(func, d_func, 0.0);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;

        case SHELL_GRAD_PC_BC:
	  shell_n_dot_pflux_bc(func, d_func,
			       bc->BC_Data_Float[0], 
			       time_value, delta_t);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );          
	  break;

	case SHELL_GRAD_PC_NOBC_BC:
	    shell_n_dot_pflux_bc(func, d_func, 0.0, 
                                 time_value, delta_t);
            surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                                           ielem_dim - 1,
                                           (int) elem_side_bc->id_side,
                                           (int) elem_side_bc->num_nodes_on_side,
                                           (elem_side_bc->local_elem_node_id) );
            break;

	case SHELL_TFMP_FREE_LIQ_BC:
          if (pd->e[pg->imtrx][R_TFMP_MASS]) {
            shell_n_dot_liq_velo_bc_tfmp(func, d_func, 0.0,
                                         time_value, delta_t,
                                         xi, exo);

            surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                   ielem_dim - 1,
                   (int) elem_side_bc->id_side,
                   (int) elem_side_bc->num_nodes_on_side,
                 (elem_side_bc->local_elem_node_id) );
          }

	  break;
	case SHELL_TFMP_NUM_DIFF_BC:
	  shell_num_diff_bc_tfmp(func, d_func, 
                                 time_value, delta_t,
                                 xi, exo);
	  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );
	  break;
        case SHELL_LUBRICATION_OUTFLOW_BC:
            shell_lubrication_outflow(func, d_func,
                                      time_value, delta_t,
                                      xi, exo);
            surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                                           ielem_dim - 1,
                                           (int) elem_side_bc->id_side,
                                           (int) elem_side_bc->num_nodes_on_side,
                                           (elem_side_bc->local_elem_node_id) );
            break;

        case SH_S11_WEAK_BC:
        case SH_S22_WEAK_BC:
          apply_shell_traction_bc(func, d_func,
                                  bc->BC_Name,
                                  bc->BC_Data_Float[0],
                                  bc->BC_Data_Float[1],
                                  bc->BC_Data_Float[2]);

            break;
	case SHELL_TFMP_AVG_PLATE_VELO_BC:
          shell_tfmp_avg_plate_velo_liq(func, d_func,
                                        time_value,
                                        delta_t,
                                        xi,
                                        exo);
          surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                                         ielem_dim - 1,
                                         (int) elem_side_bc->id_side,
                                         (int) elem_side_bc->num_nodes_on_side,
                                         (elem_side_bc->local_elem_node_id) );
          break;

	case SHELL_TFMP_FREE_GAS_BC:
	  shell_n_dot_gas_velo_bc_tfmp(func, d_func, 0.0, 
				       time_value, delta_t,
				       xi, exo);
          surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
					 ielem_dim - 1,
					 (int) elem_side_bc->id_side,
					 (int) elem_side_bc->num_nodes_on_side,
					 (elem_side_bc->local_elem_node_id) );

	  break;
	case SHELL_TFMP_GRAD_S_BC:
          shell_tfmp_n_dot_grad_s(func, d_func,
                                  time_value,
                                  delta_t,
                                  xi,
                                                                                                                exo);
          surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes,
                                         ielem_dim - 1,
                                         (int) elem_side_bc->id_side,
                                         (int) elem_side_bc->num_nodes_on_side,
                                         (elem_side_bc->local_elem_node_id) );
          break;


	case HYDROSTATIC_SYMM_BC:
	  EH(-1, "HYDROSTATIC_SYMM is no longer supported.");
	  /* 	    hydrostatic_n_dot_T(func, d_func); */
	  break;

	case VL_EQUIL_BC:

	  raoults_law(func, d_func, bc->BC_Data_Int[0],
		      bc->BC_Data_Int[1], bc->BC_Data_Int[2],
		      bc->BC_Data_Float[0], bc->BC_Data_Float[1],
		      bc->BC_Data_Float[2], bc->BC_Data_Float[3],
		      bc->BC_Data_Float[4]);
	  break;

	case YFLUX_DISC_RXN_BC:

	  yflux_disc_rxn_bc(func, d_func, bc->BC_Data_Int[0],
			    bc->BC_Data_Int[1], bc->BC_Data_Int[2],
			    bc->BC_Data_Float[0], bc->BC_Data_Float[1], delta_t, theta) ;

	  break;

	case VL_EQUIL_PRXN_BC:
	  new_way = TRUE;
	  func[0] = raoults_law_prxn(&jacCol, bc, ip, elem_side_bc, 
		  x_dot, time_value, theta, delta_t,interface_id);
	  break;

	case IS_EQUIL_PRXN_BC:
	  new_way = TRUE;
	  func[0] = is_equil_prxn(&jacCol, bc, ip, elem_side_bc,
          	  x_dot, time_value, theta, delta_t,interface_id);
	  break;

        case SDC_STEFANFLOW_BC:
	  /*
	   * This is a one sided boundary condition. The side of the
	   * boundary condition is specified on the bc card as
	   * the first integer
	   */
	  if (ei[pg->imtrx]->elem_blk_id == bc->BC_Data_Int[0]) {
	    iapply = 1;
	  } else {
	    iapply = 0;
	  }
	  if (iapply) {
	    func[0] = sdc_stefan_flow(&jacCol, bc, ip, elem_side_bc,
		      x_dot, time_value, theta, delta_t,interface_id);
	  } else {
	    skip_other_side = TRUE;
	  }
	  new_way = TRUE;
	  break;
	case SDC_KIN_SF_BC:
	  /*
	   * This is a one sided boundary condition. The side of the
	   * boundary condition is specified on the bc card as
	   * the first integer
	   */
	  if (ei[pg->imtrx]->elem_blk_id == bc->BC_Data_Int[0]) {
	    iapply = 1;
	  } else {
	    iapply = 0;
	  }
	  if (iapply) {
	    func[0] = sdc_stefan_flow(&jacCol, bc, ip, elem_side_bc,
	             x_dot, time_value, theta, delta_t, interface_id);
	  } else {
	    skip_other_side = TRUE;
	  }
	  new_way = TRUE;
	  break;
	case SDC_KIN_SFV_BC:
	  /*
	   * This is a one sided boundary condition. The side of the
	   * boundary condition is specified on the bc card as
	   * the first integer
	   */
	  if (ei[pg->imtrx]->elem_blk_id == bc->BC_Data_Int[0]) {
	    iapply = 1;
	  } else {
	    iapply = 0;
	  }
	  if (iapply) {
	    func[0] = sdc_stefan_volume_flow(&jacCol, bc, ip, elem_side_bc,
					     x_dot, theta, delta_t);
	  } else {
	    skip_other_side = TRUE;
	  }
	  new_way = TRUE;
	  break;

	case VL_POLY_BC:
	  flory_huggins(func, d_func, bc->BC_Data_Int[0],
			bc->BC_Data_Int[1], bc->BC_Data_Int[2],
			bc->BC_Data_Int[3], bc->BC_Data_Float[0]);
	  break;

	case YTOTALFLUX_CONST_BC:  /*  RSL 6/7/00  */
	  if (iapply) {
	    const_mass_flux_surf_bc(func, d_func, bc->BC_Data_Int[0],
				    bc->BC_Data_Float[0], time_value, delta_t, theta);
	  }
	  break;

	case YFLUX_BC:
	  if (iapply) {
	    mass_flux_surf_bc(func, d_func, bc->BC_Data_Int[0],
			      bc->BC_Data_Float[0], bc->BC_Data_Float[1], 
			      delta_t, theta);
	  }
	  break;

	case YFLUX_ETCH_BC:
	  if (iapply) {
	    mass_flux_surf_etch(func, d_func, bc->BC_Data_Int[0],
			        bc->BC_Data_Int[1], time_value, delta_t, theta);
	  }
	  break;

	case YFLUX_BV_BC:
	  if (iapply) {
	    mass_flux_BV_surf_bc(func, d_func,
				 BC_Types[bc_input_id].BC_Data_Int[0],
				 BC_Types[bc_input_id].BC_Data_Float[0],
				 BC_Types[bc_input_id].BC_Data_Float[1], 
				 BC_Types[bc_input_id].BC_Data_Float[2], 
				 BC_Types[bc_input_id].BC_Data_Float[3], 
				 BC_Types[bc_input_id].BC_Data_Float[4], 
				 BC_Types[bc_input_id].BC_Data_Float[5], 
				 BC_Types[bc_input_id].BC_Data_Float[6], 
				 delta_t, theta);
	  }
	  break;

        case YFLUX_HOR_BC:
	  if (iapply) {
	    mass_flux_HOR_surf_bc(func, d_func,
				  BC_Types[bc_input_id].BC_Data_Int[0],
				  BC_Types[bc_input_id].BC_Data_Float[0],
				  BC_Types[bc_input_id].BC_Data_Float[1],
				  BC_Types[bc_input_id].BC_Data_Float[2],
				  BC_Types[bc_input_id].BC_Data_Float[3],
				  BC_Types[bc_input_id].BC_Data_Float[4],
				  BC_Types[bc_input_id].BC_Data_Float[5],
				  BC_Types[bc_input_id].BC_Data_Float[6],
				  BC_Types[bc_input_id].BC_Data_Float[7],
				  BC_Types[bc_input_id].BC_Data_Float[8],
				  BC_Types[bc_input_id].BC_Data_Float[9],
				  delta_t, theta);
	  }
	  break;

        case YFLUX_ORR_BC:
	  if (iapply) {
	    mass_flux_ORR_surf_bc(func, d_func,
				  BC_Types[bc_input_id].BC_Data_Int[0],
				  BC_Types[bc_input_id].BC_Data_Float[0],
				  BC_Types[bc_input_id].BC_Data_Float[1],
				  BC_Types[bc_input_id].BC_Data_Float[2],
				  BC_Types[bc_input_id].BC_Data_Float[3],
				  BC_Types[bc_input_id].BC_Data_Float[4],
				  BC_Types[bc_input_id].BC_Data_Float[5],
				  BC_Types[bc_input_id].BC_Data_Float[6],
				  BC_Types[bc_input_id].BC_Data_Float[7],
				  BC_Types[bc_input_id].BC_Data_Float[8],
				  delta_t, theta);
	  }
	  break;

        case YFLUX_H2O_ANODE_BC:
	  if (iapply) {
	    mass_flux_H2O_ANODE_surf_bc(func, d_func,
					BC_Types[bc_input_id].BC_Data_Int[0],
					BC_Types[bc_input_id].BC_Data_Float[0],
					BC_Types[bc_input_id].BC_Data_Float[1],
					BC_Types[bc_input_id].BC_Data_Float[2],
					BC_Types[bc_input_id].BC_Data_Float[3],
					BC_Types[bc_input_id].BC_Data_Float[4],
					BC_Types[bc_input_id].BC_Data_Float[5],
					BC_Types[bc_input_id].BC_Data_Float[6],
					BC_Types[bc_input_id].BC_Data_Float[7],
					delta_t, theta);
	  }
	  break;

        case YFLUX_H2O_CATHODE_BC:
	  if (iapply) {
	    mass_flux_H2O_CATHODE_surf_bc(func, d_func,
					  BC_Types[bc_input_id].BC_Data_Int[0],
					  BC_Types[bc_input_id].BC_Data_Float[0],
					  BC_Types[bc_input_id].BC_Data_Float[1],
					  BC_Types[bc_input_id].BC_Data_Float[2],
					  BC_Types[bc_input_id].BC_Data_Float[3],
					  BC_Types[bc_input_id].BC_Data_Float[4],
					  BC_Types[bc_input_id].BC_Data_Float[5],
					  BC_Types[bc_input_id].BC_Data_Float[6],
					  BC_Types[bc_input_id].BC_Data_Float[7],
					  delta_t, theta);
	  }
	  break;

	case YFLUX_SULFIDATION_BC: 
	  if (iapply) {
	    mass_flux_SULFIDATION_surf_bc(func, d_func, 
					  BC_Types[bc_input_id].BC_Data_Int[2],
					  BC_Types[bc_input_id].BC_Data_Int[0],
					  BC_Types[bc_input_id].BC_Data_Float[0],
					  BC_Types[bc_input_id].BC_Data_Float[1], 
					  BC_Types[bc_input_id].BC_Data_Float[2], 
					  BC_Types[bc_input_id].BC_Data_Float[3], 
					  BC_Types[bc_input_id].BC_Data_Float[4], 
					  BC_Types[bc_input_id].BC_Data_Float[5], 
					  BC_Types[bc_input_id].BC_Data_Float[6], 
					  BC_Types[bc_input_id].BC_Data_Float[7], 
					  BC_Types[bc_input_id].BC_Data_Float[8], 
					  BC_Types[bc_input_id].BC_Data_Float[9], 
					  delta_t, theta);
	  }
	  break;

        case YFLUX_BV2_BC:  /* RSL 1/15/01 */
	  if (iapply) {
	    mass_flux_BV2_surf_bc(func, d_func,
				  bc->BC_Data_Int[0],
				  bc->BC_Data_Float[0],
				  bc->BC_Data_Float[1],
				  bc->BC_Data_Float[2],
				  bc->BC_Data_Float[3],
				  bc->BC_Data_Float[4],
				  bc->BC_Data_Float[5],
				  time_value, delta_t, theta);
	  }
	  break;

        case CURRENT_BV2_BC:  /* RSL 1/22/01 */
	  if (iapply) {
	    current_BV2_surf_bc(func, d_func,
				bc->BC_Data_Int[0],
				bc->BC_Data_Float[0],
				bc->BC_Data_Float[1],
				bc->BC_Data_Float[2],
				bc->BC_Data_Float[3],
				bc->BC_Data_Float[4],
				bc->BC_Data_Float[5],
				time_value, delta_t, theta);
	  }
	  break;

        case YFLUX_NI_BC:  /* RSL 3/9/01 */
	  if (iapply) {
	    mass_flux_NI_surf_bc(func, d_func,
				 bc->BC_Data_Int[0],
				 bc->BC_Data_Float[0],
				 bc->BC_Data_Float[1],
				 time_value, delta_t, theta);
	  }
	  break;

        case CURRENT_NI_BC:  /* RSL 3/9/01 */
	  if (iapply) {
	    current_NI_surf_bc(func, d_func,
			       bc->BC_Data_Int[0],
			       bc->BC_Data_Float[0],
			       bc->BC_Data_Float[1],
			       time_value, delta_t, theta);
	  }
	  break;

	case YFLUX_SUS_BC:
	  sus_mass_flux_surf_bc (func, d_func,
				 bc->BC_Data_Int[0], 
				 time_value, delta_t, theta, pg_data->h);
	  break;

	case YFLUX_CONST_BC:
	  /* no need to call function - this is too easy */
	  *func = bc->BC_Data_Float[0];
	  break;

	case YFLUX_USER_BC:
	  mass_flux_surf_user_bc(func, d_func, bc->BC_Data_Int[0],
				 bc->u_BC, time_intermediate);
	  break;

	case YFLUX_EQUIL_BC:
	  get_equil_surf_bc(func, d_func, bc->BC_Data_Int[2],
			    bc->BC_Data_Int[0], bc->BC_Data_Float[0],
			    bc->BC_Data_Float[1], bc->BC_Data_Float[2],
			    bc->BC_Data_Float[3], bc->BC_Data_Float[4],
			    delta_t, theta);  
	  break;

	case YUSER_BC:
	  yuser_surf( func, d_func,
		      BC_Types[bc_input_id].BC_Data_Int[0], 
		      BC_Types[bc_input_id].u_BC,
		      time_intermediate);
	  break;
	    
	case YFLUX_ALLOY_BC:
	  mass_flux_alloy_surf( func, d_func,
				BC_Types[bc_input_id].BC_Data_Int[0], 
				BC_Types[bc_input_id].u_BC,
				time_intermediate);
	  break;

	case SURFACE_CHARGE_BC:   
	  surface_charge_surf(func, d_func, bc->BC_Data_Float[0]);
	  break;
	case FICK_CHRGD_SURF_GRAD_BC:
	  fickian_charged_gradient_bc(func, d_func, bc->BC_Data_Int[0],bc->BC_Data_Float[0]);
	  break; 

	case POROUS_FLUX_BC:
	  porous_mass_flux_surf_bc (func, d_func,
				    bc->BC_Data_Int[0],
				    bc->BC_Data_Float[0],
				    bc->BC_Data_Float[1], 
				    bc->BC_Data_Float[2],
				    bc->BC_Data_Float[3],
				    delta_t,theta);
	  if (neg_elem_volume) return (status);
	  break;
	    
	case POROUS_CONV_BC:
	  porous_convection_bc (func, d_func,
				bc->BC_Data_Int[0],
				delta_t,theta);
	  if( neg_elem_volume ) return(status);
	  break;
	    
	case VP_EQUIL_BC:
	  porous_vapor_equil_bc(func, d_func, x_dot, theta, delta_t,
				bc_input_id, BC_Types,
				(int) bc->BC_Data_Int[0],
				(int) bc->BC_Data_Int[1],
				(int) bc->BC_Data_Int[2],
				bc->BC_Data_Float[0]);
	  break;
	    
	case POROUS_GAS_BC:
	  put_gas_flux_in_pores (func, d_func, x_dot, theta, delta_t,
				 bc_input_id, BC_Types,
				 (int) bc->BC_Data_Int[0],
				 (int) bc->BC_Data_Int[1],
				 (int) bc->BC_Data_Int[2],
				 bc->BC_Data_Float[0],
				 bc->BC_Data_Float[1]);
	  break;
	case POR_LIQ_FLUX_FILL_BC:
	  por_liq_flux_fill ( func, d_func, theta, delta_t,
			      (int) bc->BC_Data_Int[0],
			      (int) bc->BC_Data_Int[1],
			      (double) bc->BC_Data_Float[0],
			      (double) bc->BC_Data_Float[1],
			      (double) bc->BC_Data_Float[2],
			      (double) bc->BC_Data_Float[3], 0 ); 

	  break;
        case POROUS_LIQ_FLUX_CONST_BC: /* no need to call function - this is too easy */
	  *func = BC_Types[bc_input_id].BC_Data_Float[0];
	  break;

        case POROUS_GAS_FLUX_CONST_BC: /* no need to call function - this is too easy */
	  *func = BC_Types[bc_input_id].BC_Data_Float[0];
	  break;

	case QRAD_BC:
	  qrad_surf (func, d_func,
		     BC_Types[bc_input_id].BC_Data_Float[0],
		     BC_Types[bc_input_id].BC_Data_Float[1],
		     BC_Types[bc_input_id].BC_Data_Float[2],
		     BC_Types[bc_input_id].BC_Data_Float[3]);
	  break;
            
	case QCONV_BC:
	  qrad_surf (func, d_func,
		     BC_Types[bc_input_id].BC_Data_Float[0],
		     BC_Types[bc_input_id].BC_Data_Float[1],
		     0.,0.);
	  break;

	case QRAD_REPULSE_ROLL_BC:
	  qrad_surf_repulse (func, d_func,
		     BC_Types[bc_input_id].BC_Data_Float[0],
		     BC_Types[bc_input_id].BC_Data_Float[1],
		     BC_Types[bc_input_id].BC_Data_Float[2],
		     BC_Types[bc_input_id].BC_Data_Float[3],
		     BC_Types[bc_input_id].BC_Data_Float[4],
		    &(BC_Types[bc_input_id].BC_Data_Float[5]),
		    &(BC_Types[bc_input_id].BC_Data_Float[8]),
		     BC_Types[bc_input_id].BC_Data_Float[11],
		     BC_Types[bc_input_id].BC_Data_Float[12],
		     BC_Types[bc_input_id].BC_Data_Float[13],
		     BC_Types[bc_input_id].BC_Data_Float[14]);
	  break;
            
	case QUSER_BC:
	  quser_surf(func, d_func, BC_Types[bc_input_id].u_BC,
		     time_intermediate);
	  break;
	case Q_LASER_WELD_BC:
	  qlaser_surf(func, d_func, BC_Types[bc_input_id].u_BC, x,
		      time_intermediate);
	  break;

	case Q_VAPOR_BC:
	  q_vapor(func, d_func, BC_Types[bc_input_id].u_BC, x,
		  time_intermediate);
	  break;
		  
	case QSIDE_BC: /* no need to call function - this is too easy */
	  *func = bc->BC_Data_Float[0];
	  break;

	case QSIDE_DIR_BC: 
	  qside_directional(func, d_func,
			    bc->BC_Data_Float[0],
			    bc->BC_Data_Float[1],
			    bc->BC_Data_Float[2]);
	  break;
	  
	case QSIDE_LS_BC: 
	   
	  qside_ls( func,d_func,
		    bc->BC_Data_Int[0],
		    bc->BC_Data_Int[1],
		    bc->BC_Data_Float[0],
		    bc->BC_Data_Float[1]);
		    
	  break;

	case QNOBC_BC:
	  qnobc_surf (func, d_func, time_intermediate);
	  break;

	case T_CONTACT_RESIS_BC:
	  qside_contact_resis(func, d_func,
			    bc->BC_Data_Int[0],
			    bc->BC_Data_Int[1],
			    bc->BC_Data_Float[0]);
	 break;
	case T_CONTACT_RESIS_2_BC:
	  qside_contact_resis(func, d_func,
			    bc->BC_Data_Int[0],
			    bc->BC_Data_Int[1],
			    bc->BC_Data_Float[0]);
	 break;

	case APR_PLANE_TRANS_BC:
	case API_PLANE_TRANS_BC:
	case APR_VELOCITY_BC:
	case API_VELOCITY_BC:
	    acoustic_plane_transmission (func, d_func, time_intermediate, 
                      (int)bc->BC_Name, bc->BC_Data_Float[0],
                      bc->BC_Data_Float[1], bc->BC_Data_Float[2],
                      bc->BC_Data_Float[3], bc->BC_Data_Int[0]);
	    break;

	case LIGHTP_TRANS_BC:
	case LIGHTM_TRANS_BC:
	case LIGHTD_TRANS_BC:
	    light_transmission (func, d_func, time_intermediate, 
                      (int)bc->BC_Name, bc->BC_Data_Float[0],
                      bc->BC_Data_Float[1], bc->BC_Data_Float[2],
                      bc->BC_Data_Int[0]);
	    break;
	case LIGHTP_JUMP_BC:
	case LIGHTM_JUMP_BC:
	  if (ei[pg->imtrx]->elem_blk_id == bc->BC_Data_Int[0]) {
	    iapply = 1;
	  } else {
	    iapply = 1;
	  }
	  if (iapply) {
	  qside_light_jump(func, d_func, time_intermediate,(int)bc->BC_Name,
			    bc->BC_Data_Int[0],
			    bc->BC_Data_Int[1]);
	  } else {
	    skip_other_side = FALSE;
	  }
	 break;
	case LIGHTP_JUMP_2_BC:
	case LIGHTM_JUMP_2_BC:
	  qside_light_jump(func, d_func, time_intermediate,(int)bc->BC_Name,
			    bc->BC_Data_Int[0],
			    bc->BC_Data_Int[1]);
	 break;
	case APR_NOBC_BC:
	case API_NOBC_BC:
	  acoustic_nobc_surf (func, d_func, time_intermediate, (int)bc->BC_Name);
	  break;

	case POTENTIAL_NOBC_BC:
	  potential_nobc_surf (func, d_func, time_intermediate);
	  break;
 
	case LATENT_HEAT_BC:
	  if (iapply) {
	    lat_heat_bc(func, d_func, x_dot, theta,
			delta_t, bc_input_id, BC_Types);
	    if (neg_elem_volume) return(status);
	  }
	  break;

	case LATENT_HEAT_INTERNAL_BC:
	  if (iapply) {
	    lat_heat_internal_bc(func, d_func, theta, delta_t, 
				 bc_input_id, BC_Types);
	  }
	  break;

	case HEAT_OF_RXN_BC:
	  /* NB. PRS (02/09/01). To implement this surface reaction, 
	   * nonisothermal effect you should start with a routine like
	   * lat_heat, in mm_fill_species.  In there you will have to make 
	   * sure your heats of reaction are deployed based on the flux 
	   * (and hence reaction) rates, either through YFLUX_BV, or 
	   * YFLUX_USER, etc. */
	  EH(-1, "HEAT_OF_RXN_BC: not yet implemented.");
	  break;
	      
	case PSPG_BC:
	  /* MMH: LSA LSA LSA Stopped here.  ALMOST done with bc_integ.c LSA LSA LSA */
	  if (!PSPG) {
	    EH(-1,
	       "You don't have PSPG turned on and you trying to apply a PSPG boundary condition");
	  }
	  PSPG_consistency_bc(func, d_func,  x_dot, time_value, delta_t,
			      theta, pg_data);
	  break;
                         
	case FILL_CA_BC:
	  apply_CA_FILL(func, d_func,
			BC_Types[bc_input_id].BC_Data_Float[0]);             
	  break;

	case WETTING_TENSION_BC :
	  apply_wetting_tension( func,d_func,
				 BC_Types[bc_input_id].BC_Data_Float[0]);
	  break;
		
	case WETTING_SPEED_LIN_BC:
	  apply_wetting_velocity( func, d_func,
				  elem_side_bc->id_side,
				  ielem_type,
				  BC_Types[bc_input_id].BC_Data_Float[3],
				  &(BC_Types[bc_input_id].BC_Data_Float[4]),
				  BC_Types[bc_input_id].BC_Data_Float[2],
				  BC_Types[bc_input_id].BC_Data_Float[0],
				  BC_Types[bc_input_id].BC_Data_Float[1]);
	  break;
	case LINEAR_WETTING_SIC_BC:
	  apply_linear_wetting_sic( func, d_func,delta_t, theta,
				    elem_side_bc->id_side,
				    ielem_type,
				    BC_Types[bc_input_id].BC_Data_Float[3],
				    &(BC_Types[bc_input_id].BC_Data_Float[4]),
				    BC_Types[bc_input_id].BC_Data_Float[2],
				    BC_Types[bc_input_id].BC_Data_Float[0],
				    BC_Types[bc_input_id].BC_Data_Float[1],
				    BC_Types[bc_input_id].BC_Data_Float[7]);
	  break;
	case WETTING_SPEED_BLAKE_BC:
	case WETTING_SPEED_HOFFMAN_BC:
	case WETTING_SPEED_COX_BC:
	case WETTING_SPEED_SHIK_BC:
	  apply_blake_wetting_velocity( func, d_func,
					(int) bc->BC_Name,
					elem_side_bc->id_side,
					ielem_type,
					BC_Types[bc_input_id].BC_Data_Float[4],
					BC_Types[bc_input_id].BC_Data_Float[3],
					BC_Types[bc_input_id].BC_Data_Float[0],
					BC_Types[bc_input_id].BC_Data_Float[1],
					BC_Types[bc_input_id].BC_Data_Float[2]);		
		
	  break;
	case BLAKE_DIRICH_ROLL_BC:
        case HOFFMAN_DIRICH_ROLL_BC:
        case COX_DIRICH_ROLL_BC:
        case BLAKE_DIRICHLET_BC:
        case HOFFMAN_DIRICHLET_BC:
        case COX_DIRICHLET_BC:
	  {	/* trying to find wall velocity, assume velo_normal=0 */
	    double wall_velo[DIM], theta_max = 180.0;
	    int i1;
	    wall_velo[0] = BC_Types[bc_input_id].BC_Data_Float[5];
	    wall_velo[1] = BC_Types[bc_input_id].BC_Data_Float[6];
	    wall_velo[2] = BC_Types[bc_input_id].BC_Data_Float[7];
	    for (i1 = 0; i1 < Num_BC; i1++) 
	      {
		if( BC_Types[i1].BC_ID == bc->BC_ID )
		  {
		    switch(BC_Types[i1].BC_Name)
		      {
		      case VELO_TANGENT_BC:
		      case VELO_STREAMING_BC:
			wall_velo[0] = fv->snormal[1]*BC_Types[i1].BC_Data_Float[0];
			wall_velo[1] = -fv->snormal[0]*BC_Types[i1].BC_Data_Float[0];
			break;
		      case VELO_SLIP_BC:
		      case VELO_SLIP_FLUID_BC:
		      case AIR_FILM_BC:
			wall_velo[0] = BC_Types[i1].BC_Data_Float[1];
			wall_velo[1] = BC_Types[i1].BC_Data_Float[2];
			break;
		      case VELO_SLIP_ROT_BC:
		      case VELO_SLIP_ROT_FLUID_BC:
		      case AIR_FILM_ROT_BC:
			wall_velo[0] = BC_Types[i1].BC_Data_Float[1]*(fv->x[1]-BC_Types[i1].BC_Data_Float[3]);
			wall_velo[1] =BC_Types[i1].BC_Data_Float[1]*(fv->x[0]-BC_Types[i1].BC_Data_Float[2]);
			break;
		      case VELO_TANGENT_USER_BC:
			wall_velo[0] = velo_vary_fnc(UVARY_BC, fv->x[0], fv->x[1], fv->x[2] , BC_Types[i1].u_BC, time_value);
			wall_velo[1] = velo_vary_fnc(VVARY_BC, fv->x[0], fv->x[1], fv->x[2] , BC_Types[i1].u_BC, time_value);
			wall_velo[2] = velo_vary_fnc(WVARY_BC, fv->x[0], fv->x[1], fv->x[2] , BC_Types[i1].u_BC, time_value);
			break;
		      default:
			if( Debug_Flag > 1 )WH(-1,"Wall velocity bc not found\n");
                      }  /* switch bc */
                    }       /*  if BC_Types  */
                  }               /*  Num_BC loop  */
		apply_blake_wetting_velocity_sic( func, d_func,delta_t,theta,
		             		  (int) bc->BC_Name,
					  elem_side_bc->id_side,
					  ielem_type,
					  BC_Types[bc_input_id].BC_Data_Float[3],
					  BC_Types[bc_input_id].BC_Data_Float[0],
					  BC_Types[bc_input_id].BC_Data_Float[1],
					  BC_Types[bc_input_id].BC_Data_Float[2],
                      BC_Types[bc_input_id].BC_Data_Float[4],
				  	  wall_velo, theta_max, 
                                          BC_Types[bc_input_id].BC_Data_Int[3],
       				  BC_Types[bc_input_id].BC_Data_Float[9],
				  BC_Types[bc_input_id].BC_Data_Float[10],
				  BC_Types[bc_input_id].BC_Data_Float[11]);		
		}
		break;
		
	case SHIK_DIRICH_ROLL_BC:
        case SHIK_DIRICHLET_BC:
		apply_blake_wetting_velocity_sic( func, d_func,delta_t,theta,
		             		  (int) bc->BC_Name,
					  elem_side_bc->id_side,
					  ielem_type,
					  BC_Types[bc_input_id].BC_Data_Float[3],
					  BC_Types[bc_input_id].BC_Data_Float[0],
					  BC_Types[bc_input_id].BC_Data_Float[1],
					  BC_Types[bc_input_id].BC_Data_Float[2],
					  BC_Types[bc_input_id].BC_Data_Float[4],
					&(BC_Types[bc_input_id].BC_Data_Float[5]),
					  BC_Types[bc_input_id].BC_Data_Float[8],
					  BC_Types[bc_input_id].BC_Data_Int[3],  
       				  BC_Types[bc_input_id].BC_Data_Float[9],
				  BC_Types[bc_input_id].BC_Data_Float[10],
				  BC_Types[bc_input_id].BC_Data_Float[11]);		
		break;
		
	case HYSTERESIS_WETTING_BC:
	  apply_hysteresis_wetting_sic ( func, d_func, delta_t, theta, &(BC_Types[bc_input_id].BC_Data_Float[0]) );
	  break;
		
	case STRONG_FILL_CA_BC :
	  /* nothing to be done here, this condition applied to FILL EQN */            
	  break;
		  
	case TABLE_WICV_BC:
	case TABLE_WICS_BC:
	  apply_table_wic_bc(func, d_func, &BC_Types[bc_input_id], time_value);
	  break;

	case H_FREE_BC :
	  boundary_curvature( func,d_func, NULL );
	  break;

	case LS_CA_H_BC :
	  boundary_curvature( func,d_func, &(BC_Types[bc_input_id].BC_Data_Float[0] ) );
	  break;

	case CURRENT_BC: /* no need to call function - this is too easy */
	  *func = bc->BC_Data_Float[0];
	  break;

 	case CURRENT_SIC_BC: /* strong integrated bc on potential eqn  */
	  apply_potential_grad_bc(func, d_func, &BC_Types[bc_input_id], time_intermediate, delta_t);
	  break;

	case CURRENT_USER_BC: 
	  current_user_surf(func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	  break;

	case CURRENT_USER_SIC_BC: 
	  current_user_surf(func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
	  apply_potential_grad_bc(func, d_func, &BC_Types[bc_input_id], time_intermediate, delta_t);
	  break;

        case VOLT_USER_BC:
	  volt_user_surf(func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);

	  break;

	case CURRENT_BV_BC: 
	  if (iapply) {
	    current_BV_surf(func, d_func, 
			    BC_Types[bc_input_id].BC_Data_Int[0],
			    BC_Types[bc_input_id].BC_Data_Float[0],
			    BC_Types[bc_input_id].BC_Data_Float[1], 
			    BC_Types[bc_input_id].BC_Data_Float[2], 
			    BC_Types[bc_input_id].BC_Data_Float[3], 
			    BC_Types[bc_input_id].BC_Data_Float[4], 
			    BC_Types[bc_input_id].BC_Data_Float[5], 
			    BC_Types[bc_input_id].BC_Data_Float[6], 
			    delta_t, theta);
	  }
	  break;

	case CURRENT_ORR_BC:
	  if (iapply) {
	    current_ORR_surf(func, d_func,
			     BC_Types[bc_input_id].BC_Data_Int[0],
			     BC_Types[bc_input_id].BC_Data_Float[0],
			     BC_Types[bc_input_id].BC_Data_Float[1],
			     BC_Types[bc_input_id].BC_Data_Float[2],
			     BC_Types[bc_input_id].BC_Data_Float[3],
			     BC_Types[bc_input_id].BC_Data_Float[4],
			     BC_Types[bc_input_id].BC_Data_Float[5],
			     BC_Types[bc_input_id].BC_Data_Float[6],
			     BC_Types[bc_input_id].BC_Data_Float[7],
			     delta_t, theta);
	  }
	  break;

	case CURRENT_HOR_BC:
	  if (iapply) {
	    current_HOR_surf(func, d_func,
			     BC_Types[bc_input_id].BC_Data_Int[0],
			     BC_Types[bc_input_id].BC_Data_Float[0],
			     BC_Types[bc_input_id].BC_Data_Float[1],
			     BC_Types[bc_input_id].BC_Data_Float[2],
			     BC_Types[bc_input_id].BC_Data_Float[3],
			     BC_Types[bc_input_id].BC_Data_Float[4],
			     BC_Types[bc_input_id].BC_Data_Float[5],
			     BC_Types[bc_input_id].BC_Data_Float[6],
			     BC_Types[bc_input_id].BC_Data_Float[7],
			     BC_Types[bc_input_id].BC_Data_Float[8],
			     delta_t, theta);
	  }
	  break;

        case SH_SDET_BC:
          apply_sdet(func, d_func, xi, exo);

          break;
        case SH_MESH2_WEAK_BC:
          apply_sh_weak(func, d_func, xi, exo,
                       BC_Types[bc_input_id].BC_Data_Float[0]);

          break;

	default:
	  sprintf(Err_Msg, "Integrated BC %s not found", bc_desc->name1);
	  EH(-1, Err_Msg);
	  break;


	} /* end of switch over bc type */

        /*
	 *  If we have determined that the current boundary condition need 
	 *  only be applied on one side of the interface and we are currently
	 *  on the other side of the interface, skip processing the rest
	 *  of this surface integral for this particular boundary condition.
	 *  Go directly to the next surface boundary condition on the current
	 *  side.
	 */
        if (skip_other_side) continue;
		
	/**********************************************************************/
	/*        LOOP OVER THE LOCAL ELEMENT NODE NUMBER                     */
	/*        FOR NODES THAT LIE ON THE CURRENT SURFACE		      */
	/* this is a loop designed to loop over equations that lie on         */
	/*  current surface                                                   */
	/*    ADD the Boundary condition functions into the Residual          */
	/*         vector and Jacobian Matrix                                 */
	/**********************************************************************/

	for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++)  {
	      
	  /* Find the local element node number for the current node */
	  id = (int) elem_side_bc->local_elem_node_id[i];
	      
	  /* Find the processor node number given the
	   * local element node number,  'i'
	   */
	  I = Proc_Elem_Connect[iconnect_ptr + id];

	  /*
	   * If BC is capillarity, convert the node-based cfunc to func
	   */

	  if ((BC_Types[bc_input_id].BC_Name == CAPILLARY_BC || 
	       BC_Types[bc_input_id].BC_Name == CAP_REPULSE_BC ||
	       BC_Types[bc_input_id].BC_Name == CAP_REPULSE_ROLL_BC ||
	       BC_Types[bc_input_id].BC_Name == CAP_REPULSE_USER_BC ||
	       BC_Types[bc_input_id].BC_Name == CAP_REPULSE_TABLE_BC ||
	       BC_Types[bc_input_id].BC_Name == CAP_RECOIL_PRESS_BC ||
	       BC_Types[bc_input_id].BC_Name == CAPILLARY_TABLE_BC ||
	       BC_Types[bc_input_id].BC_Name == ELEC_TRACTION_BC ||
	       BC_Types[bc_input_id].BC_Name == CAPILLARY_SHEAR_VISC_BC) &&
	      (Dolphin[pg->imtrx][I][VELOCITY1] > 0 )) { /* DRN: getting segfault for Q1P0 on Q2 mesh, does this fix it? */
	    for (p = 0; p < bc_desc->vector; p++) {
	      /* assuming momentum and species have same interpolation, for now */
	      func[p] = cfunc[ei[pg->imtrx]->ln_to_first_dof[R_MOMENTUM1][id]][p];
	      for (k = 0; k < MAX_VARIABLE_TYPES + MAX_CONC; k++) {
		for (j = 0; j < MDE; j++) {
		  d_func[p][k][j] = d_cfunc[ei[pg->imtrx]->ln_to_first_dof[R_MOMENTUM1][id]][p][k][j];
		}
	      }
	    }
	  }
	  
	  
	  if( (BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC ) && (Dolphin[pg->imtrx][I][MESH_DISPLACEMENT1] > 0 ) )
	    {
	      func[0] = cfunc[ei[pg->imtrx]->ln_to_first_dof[R_MESH1][id]][0];
	      for (k = 0; k < MAX_VARIABLE_TYPES + MAX_CONC; k++) {
		for (j = 0; j < MDE; j++) {
		  d_func[0][k][j] = d_cfunc[ ei[pg->imtrx]->ln_to_first_dof[R_MESH1][id]  ][0][k][j];
		}
	      }
	    }
	  
	  if( (BC_Types[bc_input_id].BC_Name == SHEAR_TO_SHELL_BC ) && (Dolphin[pg->imtrx][I][R_SHELL_TENSION] > 0 ) )
	    {
	      func[0] = cfunc[ei[pg->imtrx]->ln_to_first_dof[R_SHELL_TENSION][id]][0];
	      for (k = 0; k < MAX_VARIABLE_TYPES + MAX_CONC; k++) {
		for (j = 0; j < MDE; j++) {
		  d_func[0][k][j] = d_cfunc[ ei[pg->imtrx]->ln_to_first_dof[R_SHELL_TENSION][id]][0][k][j];
		}
	      }
	    }

          /* Here , we are going to determine whether it is a stress BCs or not */
          if (bc_desc->equation == R_STRESS11)
            {
             stress_bc = 1;
            }
          else
            {
             stress_bc = 0;
            }

          /* If it is not a stress BC go for the loop over vector components */
          if (stress_bc == 0) {
          /*
	   * Boundary condition may actually be a vector of
	   * bc's. Loop over that vector here.
	   */
	  for (p = 0; p < bc_desc->vector; p++) {
	    /* 
	     *   Check to see if this BC on this node is
	     *   applicable (i.e. no other overriding Dirichlet conditions),
	     *   And, find the global unknown number, index_eq, on which 
	     *   to applyi this additive boundary condition, eqn
	     */
            index_eq = bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn,
				    p, &eqn, &matID_apply, &vd);

	    if (index_eq >= 0) {
	      /*
	       * Obtain the first local variable degree of freedom
	       * at the current node, whether or not it actually an
	       * interpolating degree of freedom
	       */
	      ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

	      /*
	       *   for weakly integrated boundary conditions,
	       *   weight the function by wt
	       */
	      weight = wt;

	      /*
	       * Weight the function by the weighting function unless it is 
	       * integrated by parts twice
	       */
	      if (BC_Types[bc_input_id].BC_Name != CAPILLARY_BC  &&
		  BC_Types[bc_input_id].BC_Name != CAPILLARY_SHEAR_VISC_BC  &&
		  BC_Types[bc_input_id].BC_Name != ELEC_TRACTION_BC  &&
                  BC_Types[bc_input_id].BC_Name != YFLUX_USER_BC &&
		  BC_Types[bc_input_id].BC_Name != CAP_REPULSE_BC &&
		  BC_Types[bc_input_id].BC_Name != CAP_REPULSE_ROLL_BC &&
		  BC_Types[bc_input_id].BC_Name != CAP_REPULSE_USER_BC &&
		  BC_Types[bc_input_id].BC_Name != CAP_REPULSE_TABLE_BC &&
		  BC_Types[bc_input_id].BC_Name != CAP_RECOIL_PRESS_BC &&
		  BC_Types[bc_input_id].BC_Name != CAPILLARY_TABLE_BC &&
		  BC_Types[bc_input_id].BC_Name != TENSION_SHEET_BC &&
		  BC_Types[bc_input_id].BC_Name != SHEAR_TO_SHELL_BC ) {
		/*
		 * Do processing specific to CROSS_PHASE_DISCONTINUOUS
		 * boundary conditions.
		 */
		if (bc_desc->i_apply == CROSS_PHASE_DISCONTINUOUS) {
		  /*
		   * For these boundary conditions, we need to extract
		   * the correct nodal basis function for these
		   * discontinuous interpolations. Here we take
		   * all the "0" basis functions that Baby_Dolphin[pg->imtrx]
		   * refers to, as the "1" types have been zeroed.
		   *
		   *  In other words, the local variable dof for the
		   *  variable corresponding to the Elem_Blk 1 is always
		   *  the first degree of freedom at a local node. The
		   *  lvdof for the variable corresponding to Elem_Blk 2 is
		   *  always the second degree of freedom for that variable
		   *  type at the node. The basis function for dofs which
		   *  are not interpolating dofs are set to zero, previously
		   *  Therefore, in the loop below, we pick phi_i to
		   *  be the basis function which is not zeroed out in
		   *  the current element.
		   */
                  ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][ldof_eqn];
                  mn_first = ei[pg->imtrx]->matID_ledof[ledof];

                  if((BC_Types[bc_input_id].BC_Data_Int[1] && 
                        Current_EB_ptr->Elem_Blk_Id == 
                              BC_Types[bc_input_id].BC_Data_Int[1])  ||  
                      (!BC_Types[bc_input_id].BC_Data_Int[1] &&
                      ((Current_EB_ptr->Elem_Blk_Id  +1) % 2) == 0))   {

		    phi_i = bf[eqn]->phi[ldof_eqn];
		    weight *= phi_i;
		  } else {
		    phi_i = bf[eqn]->phi[ldof_eqn+1];
		    weight *= phi_i;
		  }

                  if (bc_desc->DV_Index_Default == DVI_MULTI_PHASE_VD) {
		    if (mn_first != vd->MatID) {
		      ldof_eqn++;
                    }
		  }

		  /*
		   *  TIE Boundary conditions:
		   *    These strongly integrated Dirichlet boundary
		   *    conditions get applied on the second degree
		   *    of freedom at the node.
		   */
		  if ((bf[eqn]->interpolation == I_Q2_D || 
		       bf[eqn]->interpolation == I_Q1_D || 
		       bf[eqn]->interpolation == I_Q2_D_LSA) &&
		      (bc->BC_Name == CONT_TANG_VEL_BC ||
		       bc->BC_Name == CONT_NORM_VEL_BC ||
		       bc->BC_Name == DISCONTINUOUS_VELO_BC ||
		       bc->BC_Name == VL_EQUIL_BC ||
		       bc->BC_Name == VL_POLY_BC ||
		       bc->BC_Name == SDC_STEFANFLOW_BC ||
		       bc->BC_Name == SDC_KIN_SF_BC ||
		       bc->BC_Name == LIGHTP_JUMP_BC ||
		       bc->BC_Name == LIGHTM_JUMP_BC ||  
		       bc->BC_Name == T_CONTACT_RESIS_2_BC)) {
		    ldof_eqn += 1;
		  }

		}
		/*
		 *  Handle the case of SINGLE_PHASE boundary conditions and
		 *  CROSS_PHASE boundary conditions for the phase in which
		 *  the eqn actually is solved for.
		 */
		else if (bc_desc->i_apply == SINGLE_PHASE || 
			 pd->e[pg->imtrx][eqn]) {
  		  if (bc->BC_Name == KINEMATIC_PETROV_BC ||
  		      bc->BC_Name == VELO_NORMAL_LS_PETROV_BC ||
  		      bc->BC_Name == KIN_DISPLACEMENT_PETROV_BC) {
		    if (pd->Num_Dim != 2) {
 		      EH(-1,"KINEMATIC_PETROV or KIN_DISPLACEMENT_PETROV not available in 3D yet");
		    }
		    id_side = elem_side_bc->id_side;
		    i_basis = 1 - id_side%2;
		    phi_i = bf[eqn]->dphidxi[ldof_eqn][i_basis];
		    weight *= phi_i;
		  } else {
		    phi_i = bf[eqn]->phi[ldof_eqn];
		    weight *= phi_i;
		  }
		}
		/*
		 *  Handle the case of CROSS_PHASE boundary conditions in the
		 *  phase in which the eqn isn't solved for.
		 */
		else if (bc_desc->i_apply == CROSS_PHASE) {
		  /*
		   *  We are here when the equation type doesn't exist
		   *  in this material. Find the intepolation from the
		   *  adjacent material
		   */

		  for (mn = 0; mn < upd->Num_Mat; mn++) {
                    if (pd_glob[mn]->e[pg->imtrx][eqn] &&
			(eb_in_matrl(BC_Types[bc_input_id].BC_Data_Int[0], mn) ||
			 eb_in_matrl(BC_Types[bc_input_id].BC_Data_Int[1], mn)))
		      {
			//type = pd_glob[mn]->w[eqn];
			//if (bfi[type] == NULL) EH(-1,"Illegal cross basis func");
			
			/* note that here, we don't have the ln_to_dof
			   array for the adjacent 
			   material - for now assume that ldof_eqn = id */
			/* This further means that if I_Q2_D or I_Q1_D
			   are used for velocity, then we
			   cannot apply these sorts of BCs
			   --ADD DIAGNOSTIC  */
			
			//phi_i = bfi[type]->phi[id];
			
			/* DSH 08/2016
			 * The above was the old way of loading up basis functions for 
			 * CROSS_PHASE boundary conditions when we are in the adjacent 
			 * material.  This approach breaks down when considering shells
			 * mixed with continuum elements.  In either case, bf[eqn] works,
			 * so I am not sure why the bfi[type] was used in the first place.
			 */
			
			phi_i = bf[eqn]->phi[id];
			weight *= phi_i;
		      }
		  }
		} else {
		  EH(-1,"Illegal bc phase definition");
		}
	      }

	      /*
	       * For strong conditions weight the function by BIG_PENALTY
	       */
	      if (bc_desc->method == STRONG_INT_SURF ) {
		weight *= BIG_PENALTY;
	      }
	      /*
	       *   Add in the multiplicative constant for corresponding to
	       *   all boundary conditions, except for certain special
	       *   cases
	       */
	      if (bc_desc->method == WEAK_INT_SURF
		  && bc->BC_Name != PSPG_BC
		  && bc->BC_Name != VELO_SLIP_SOLID_BC
		  && bc->BC_Name != ELEC_TRACTION_BC 
		  && bc->BC_Name != QSIDE_LS_BC
		  && bc->BC_Name != LS_ADC_BC
		  && bc->BC_Name != SHEAR_TO_SHELL_BC 
		  && bc->BC_Name != POR_LIQ_FLUX_FILL_BC 
		  && bc->BC_Name != DARCY_LUB_BC
		  && bc->BC_Name != SHELL_TFMP_FREE_LIQ_BC
		  && bc->BC_Name != SHELL_TFMP_NUM_DIFF_BC
                  && bc->BC_Name != SH_SDET_BC
                  && bc->BC_Name != SH_MESH2_WEAK_BC
                  && bc->BC_Name != SHELL_LUBRICATION_OUTFLOW_BC ) {
		weight *= pd->etm[pg->imtrx][eqn][(LOG2_BOUNDARY)];
	      }

	      /*
	       * Calculate the position in the local element residual
	       * vector to put the current contribution
	       * -> MASS_FRACTION unknowns get stuck at the end of
	       *    this vector.
	       */
	      if (eqn == R_MASS) {
		ieqn = MAX_PROB_EQN + bc->species_eq;
	      } else {
		ieqn = upd->ep[pg->imtrx][eqn];
                if (goma_automatic_rotations.automatic_rotations && (bc->desc->rotate != NO_ROT)) {
                  ieqn = equation_index_auto_rotate(elem_side_bc, I, eqn, p, ldof_eqn, bc);
                }
              }

	      /*
	       *  Add the current contribution to the local element
	       *  residual vector
	       */
	      if (ldof_eqn != -1) {
                lec->R[ieqn][ldof_eqn] += weight * fv->sdet * func[p];

		
#ifdef DEBUG_BC
		if (IFPD == NULL) IFPD = fopen("darcy.txt", "a");
		fprintf (IFPD,
			 "ielem = %d: BC_index = %d, lec->R[%d][%d] += weight"
			 "* fv->sdet * func[p]: weight = %g, fv->sdet = %g, func[%d] = %g\n",
			 ei[pg->imtrx]->ielem, bc_input_id, ieqn, ldof_eqn,
			 weight, fv->sdet, p, func[p]);
		fflush(IFPD);
#endif
	        /* 
	         *   Add sensitivities into matrix
	         *  - find index of sensitivity in matrix
	         *     (if variable is not defined at this node,
	         *      loop over all dofs in element)
	         *  - add into matrix
	         */
		      
		if (af->Assemble_Jacobian && ldof_eqn != -1) {
		
		  if (new_way) {
		    for (w = 0; w < jacCol.Num_lvdesc; w++) {
		      lvdesc = jacCol.Lvdesc_Index[w];
		      tmp = jacCol.JacCol[w] * weight * fv->sdet;
		      /*
		       *  Find the variable type that corresponds to 
		       *  the local variable description. This is the
		       *  variable type for the column.
		       */
		      if (lvdesc < 0 || lvdesc > MAX_LOCAL_VAR_DESC) {
			printf("we have an error\n");
		      }
		      var = ei[pg->imtrx]->Lvdesc_to_Var_Type[lvdesc];
		      if (var == MASS_FRACTION) {
			pvar = MAX_PROB_VAR + ei[pg->imtrx]->Lvdesc_to_MFSubvar[lvdesc];
		      } else {
			pvar = upd->vp[pg->imtrx][var];
		      }
		      /*
		       *  ieqn = upd->ep[pg->imtrx][eqn] (MF's put high)
		       *         Variable type of the row
		       *
		       *  HKM Note: This sum is over all of the basis
		       *            functions in an element. However,
		       *            we know that the only nonzero basis
		       *            functions will be ones corresponding
		       *            to nodes on the current side of the
		       *            element. We should make use of that
		       *            feature to cut down the amount of work.
		       */
		      jac_ptr = lec->J[ieqn][pvar][ldof_eqn];
		      phi_ptr = bf[var]->phi;
		      for (jlv = 0; jlv < ei[pg->imtrx]->Lvdesc_Numdof[lvdesc]; jlv++) {
			j = ei[pg->imtrx]->Lvdesc_to_lvdof[lvdesc][jlv];
			lnn = ei[pg->imtrx]->Lvdesc_to_Lnn[lvdesc][jlv];
			q = ei[pg->imtrx]->ln_to_dof[var][lnn];
			jac_ptr[j] += tmp * phi_ptr[q];
		      }
		    }
	
		    if (!af->Assemble_LSA_Mass_Matrix) {
		      tmp = weight * fv->sdet;
		      for (w = 0; w < jacCol.Num_lvdof; w++) {
			var = jacCol.Lvdof_var_type[w];
			pvar = upd->vp[pg->imtrx][var];
			j = jacCol.Lvdof_lvdof[w];
			lec->J[ieqn][pvar][ldof_eqn][j] += 
			  tmp * jacCol.Jac_lvdof[w];
		      }

		      for (q = 0; q < pd->Num_Dim; q++) {
			var = MESH_DISPLACEMENT1 + q;
			if (pd->v[pg->imtrx][var]) {
			  pvar = upd->vp[pg->imtrx][var];
			  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
			    lec->J[ieqn][pvar][ldof_eqn][j] +=
			      weight * func[p] * fv->dsurfdet_dx[q][j];
			  }
			}
		      }
		    }
		  } else {
		    /* OLD METHOD */

		    /* if mesh displacement is variable,
		     *  put in this sensitivity first
		     * ... unless we are computing the mass matrix
		     * for LSA.  In that case, we don't include
		     * this first term b/c it doesn't involve any
		     * primary time derivative variables.
		     */
		    if (!af->Assemble_LSA_Mass_Matrix) {
		      for (q = 0; q < pd->Num_Dim; q++) {
			var = MESH_DISPLACEMENT1 + q;
			pvar = upd->vp[pg->imtrx][var];
			if (pvar != -1) {
			  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
			    lec->J[ieqn][pvar][ldof_eqn][j] +=
			      weight * func[p] * fv->dsurfdet_dx[q][j];
			  }
			}
		      }
		    }
			
		    /* now add in sensitivity of BC function to
		     * variables
		     */
		    for (var=0; var < MAX_VARIABLE_TYPES; var++) {
		      pvar = upd->vp[pg->imtrx][var];
		      if (pvar != -1 &&
			  (BC_Types[bc_input_id].desc->sens[var] ||	1)) {
			if (var != MASS_FRACTION) {
			  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
			    lec->J[ieqn][pvar] [ldof_eqn][j] +=
			      weight * fv->sdet * d_func[p][var][j];
			  }
			} else {
			  /* variable type is MASS_FRACTION */
			  for (w = 0; w < pd->Num_Species_Eqn; w++) {
			    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
			      lec->J[ieqn][MAX_PROB_VAR + w][ldof_eqn][j] += 
				weight * fv->sdet *
				d_func[p][MAX_VARIABLE_TYPES + w][j];
			    }
			  } /* end of loop over species */
			} /* end of if MASS_FRACTION */
		      } /* end of variable exists and BC is sensitive to it */
		    } /* end of var loop over variable types */
		  } /* End of loop over new way */
		} /* end of NEWTON */
	      }
	    } /* end of if (Res_BC != NULL) - i.e. apply residual at this node */
	  } /* end of loop over equations that this condition applies to */
          } /* end of if it is not a stress BC */

          /* Stress BC is handled in different loop so that it is not too invasive to the
           * already overloaded loop */
          if (stress_bc == 1) {

             /* For a stress BC, we will loop over modes on top of loop over the stress components */
             for (imode = 0; imode < vn->modes; imode++) {
                for (p = 0; p < bc_desc->vector; p++) {
                   /*
                    *   Check to see if this BC on this node is
                    *   applicable (i.e. no other overriding Dirichlet conditions),
                    *   And, find the global unknown number, index_eq, on which
                    *   to applying this additive boundary condition, eqn
                    */
                    index_eq = bc_eqn_index_stress(id, I, bc_input_id, ei[pg->imtrx]->mn,
                                                   p, imode, &eqn, &matID_apply, &vd);

                    if (index_eq >= 0) {
                       /*
                        * Obtain the first local variable degree of freedom
                        * at the current node, whether or not it actually an
                        * interpolating degree of freedom
                        */
                        ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

                       /*
                        *   for weakly integrated boundary conditions,
                        *   weight the function by wt
                        */
                        weight = wt;

                       /*
                        *  Handle the case of SINGLE_PHASE boundary conditions
                        */

                        if (bc_desc->i_apply == SINGLE_PHASE) {
                           phi_i = bf[eqn]->phi[ldof_eqn];
                           weight *= phi_i;
                        }
                        else {
                           EH(-1,"Only SINGLE_PHASE is handled in stress BC implementation");
                        }

                       /*
                        * For strong conditions weight the function by BIG_PENALTY
                        */
                        if (bc_desc->method == STRONG_INT_SURF ) {
                           weight *= BIG_PENALTY;
                        }

                       /*
                        * Determine the position in the local element residual
                        * vector to put the current contribution
                        */
                        ieqn = upd->ep[pg->imtrx][eqn];

                       /*
                        *  Add the current contribution to the local element
                        *  residual vector
                        */

                        lec->R[ieqn][ldof_eqn] += weight * fv->sdet * func_stress[imode][p];

                       /*
                        *   Add sensitivities into matrix
                        *  - find index of sensitivity in matrix
                        *     (if variable is not defined at this node,
                        *      loop over all dofs in element)
                        *  - add into matrix
                        */

                        if (af->Assemble_Jacobian && ldof_eqn != -1) {

                          /* if mesh displacement is variable,
                           *  put in this sensitivity first
                           * ... unless we are computing the mass matrix
                           * for LSA.  In that case, we don't include
                           * this first term b/c it doesn't involve any
                           * primary time derivative variables.
                           */

                           if (!af->Assemble_LSA_Mass_Matrix) {
                              for (q = 0; q < pd->Num_Dim; q++) {
                                  var = MESH_DISPLACEMENT1 + q;
                                  pvar = upd->vp[pg->imtrx][var];
                                  if (pvar != -1) {
                                     for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                                         lec->J[ieqn][pvar][ldof_eqn][j] +=
                                         weight * func_stress[imode][p] * fv->dsurfdet_dx[q][j];
                                     }
                                 }
                             }
                           }

                          /* now add in sensitivity of BC function to
                           * variables
                           */

                           for (var=0; var < MAX_VARIABLE_TYPES; var++) {
                               pvar = upd->vp[pg->imtrx][var];
                               if (pvar != -1) {

                                  /* Case for variable type that is not MASS_FRACTION */
                                  if (var != MASS_FRACTION) {
                                     for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                                          lec->J[ieqn][pvar] [ldof_eqn][j] +=
                                          weight * fv->sdet * d_func_stress[imode][p][var][j];
                                     }
                                  }
                                 /* Case for variable type that is MASS_FRACTION */
                                  else {
                                     for (w = 0; w < pd->Num_Species_Eqn; w++) {
                                         for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                                             lec->J[ieqn][MAX_PROB_VAR + w][ldof_eqn][j] +=
                                             weight * fv->sdet * d_func_stress[imode][p][MAX_VARIABLE_TYPES + w][j];
                                         }
                                     }
                                  }
                               }
                           }
                        } /* End of if assemble Jacobian */
                    } /* end of if (Res_BC != NULL) - i.e. apply residual at this node */
                } /* End of loop over stress components */
             } /* End of loop over stress modes */
          }/* end of if it is a stress BC */

	}  /* end for (i=0; i< num_nodes_on_side; i++) */
      }  /*End (if INT) (CAPILLARY and KINEMATIC and VELO_NORMAL and VELO_TANGENT . . .) */
    } /*(end for ibc) */
  } /*End for ip = 1,...*/  
  return (status);
}
/* END of routine apply_integrated_bc */

int equation_index_auto_rotate(const ELEM_SIDE_BC_STRUCT *elem_side_bc,
                               int I,
                               int eqn,
                               int p,
                               int ldof_eqn,
                               const BOUNDARY_CONDITION_STRUCT *bc) {
  int ieqn;
  if (!goma_automatic_rotations.automatic_rotations) {
    EH(-1, "equation_index_auto_rotate requires 3D automatic rotations");
    return -1;
  }
  int node = ei[pg->imtrx]->gnn_list[eqn][ldof_eqn];
  int n_index = -1;
  for (int i = 0; i < goma_automatic_rotations.rotation_nodes[node].n_normals; i++) {
    if (goma_automatic_rotations.rotation_nodes[node].element[i] == ei[pg->imtrx]->ielem &&
        goma_automatic_rotations.rotation_nodes[node].face[i] == elem_side_bc->id_side) {
      n_index = i;
      break;
    }
  }
  EH(n_index, "Rotations incorrectly setup");
  int rot_dir = (int) goma_automatic_rotations.rotation_nodes[I].face_coordinate_association[n_index];

  int t1dir = 0;
  int t2dir = 0;
  for (int k = 0; k < DIM; k++) {
    if (k != rot_dir) {
      t1dir = k;
      break;
    }
  }
  for (int k = 0; k < DIM; k++) {
    if (k != rot_dir && k != t1dir) {
      t2dir = k;
    }
  }
  int eq_idx[DIM];
  eq_idx[0] = rot_dir;
  eq_idx[1] = t1dir;
  eq_idx[2] = t2dir;
  ieqn = upd->ep[pg->imtrx][eqn + eq_idx[p]];
  return ieqn;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void
apply_table_wic_bc( double func[],
		    double d_func[][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
		    struct Boundary_Condition *BC_Type,
		    double time_value)

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
  BC_Type = pointer to Boundary Condition structure,
  i.e. &(BC_Types[bc_input_id]
*/
{
  int  basis;
  double  slope, interp_val, x_table[2];

  if (af->Assemble_LSA_Mass_Matrix) return;

  basis = BC_Type->table->t_index[0];

  /*  Setup Dummy variable to pass array to interpolation function */

  if ( basis != -1 )
    x_table[0]=fv->x[basis];
  else
    x_table[0] = time_value;

  if (BC_Type->table->interp_method == BIQUADRATIC
      || BC_Type->table->interp_method == BILINEAR )
    {
      x_table[1] = fv->x[BC_Type->table->t_index[1]];
    }
  interp_val = interpolate_table(BC_Type->table, x_table, &slope, NULL);

  /*   the integrand for the weak integrated conditions is passed
   *     back through table->slope structure because it's
   *    convenient.
   */
  if (BC_Type->BC_Name == TABLE_WICV_BC) {
    func[0] = BC_Type->table->slope[0]*BC_Type->BC_Data_Float[0];
    func[1] = BC_Type->table->slope[1]*BC_Type->BC_Data_Float[0];
    func[2] = BC_Type->table->slope[2]*BC_Type->BC_Data_Float[0];
  } else {
    func[0] = interp_val*BC_Type->BC_Data_Float[0];
  }
}
/*******************************************************************************/
