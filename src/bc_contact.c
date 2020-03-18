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
 
/* Routines for calculating Boundary Conditions and Adding them to matrix_fill
 *
 * Modified by EDW to place cross-mesh sensitivities into
 * overlap augmenting condition terms (bAC and cAC)
 */

/*
 *$Id: bc_contact.c,v 5.13 2010-01-10 23:00:17 hkmoffa Exp $
 */

/* Standard include files */
 
#include <stdio.h>
#include <string.h>
#include <math.h>
 
/* GOMA include files */
 
#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_mp_const.h"
#include "mm_fill_common.h"
#include "mm_eh.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_as_alloc.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_species.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_ns_bc.h"
#include "mm_post_proc.h"
#include "mm_unknown_map.h"
#include "rf_bc.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "sl_util.h"

#define GOMA_BC_CONTACT_C


int
apply_contact_bc (
          double x[],           /* Solution vector for the current processor    */
          double resid_vector[],/* Residual vector for the current processor    */
          const double delta_t, /* current time step size                       */
          const double theta,	/* parameter (0 to 1) to vary time integration
                                 *  ( implicit - 0 to explicit - 1)             */
          const double h_elem_avg, /* global average element size */
          const double h[DIM],     /* average element size */
          const double mu_avg,     /* average element viscosity */
          const double U_norm,     /* global velocity norm */
          struct elem_side_bc_struct *first_elem_side_BC_array[],
	                         /* An array of pointers to the first surface
	                          * integral defined for each element.
                                  * It's length is equal to the total number of
                                  * elements defined on the current processor */
          const int ielem,       /* element number */
          const int ielem_type,  /* element type */
          const int num_local_nodes,
          const int ielem_dim,
          const int iconnect_ptr,
          struct elem_side_bc_struct *elem_side_bc,
	                        /* Pointer to an element side boundary condition
	                         * structure */
          const int num_total_nodes,
          const int bc_application, /* flag indicating whether to integrate
				     * strong or weak BC's */
          const int oAC,        /* Base augmenting condition number when
                                 * called for AC's; -1 when called for BC's */
          double *gAC,          /* Augmenting condition arrays */
          double **bAC,
          double **cAC,
          double **dAC,
          const double time_value,
	  const Exo_DB *exo)

    /****************************************************************************
     *
     * apply_contact_bc():
     *
     *    Calculate the local element contributions to boundary conditions which
     *    involve stresses and quantities coming from another disconnected mesh.
     *
     ****************************************************************************/
{ 
  int ip,  i, I, ibc, j, id, icount;
  int aa, bb, fluid_mesh_elem_id; /*new stuff for new contact algorithm */
  int ac_lm, pass, matID_apply, id_side;
  int apply_AC, ioffset=0, iAC, iconn_ptr,   nu;
  int eqn, ieqn, var, pvar, p, q, index_eq, ldof_eqn;
  int err, status = 0;
  int bc_input_id, ip_total; 
  int dof_l[DIM] = {0, 0, 0}, dof_q[DIM] = {0, 0, 0}, nunk[DIM][MDE];
 
			      
  /* PRS new stuff....clean up unused */
  double phi_l[DIM][MDE], phi_q[DIM][MDE];
  double xi_2[DIM];  /* PRS: new stuff */
  double lagrange_mult[3] = {0., 0.,  0.};
  double fluid_velocity[3] = {0., 0., 0.};
  double lm_fluid[3] = {0., 0., 0.};

  double  phi_i;
  double s, t, u;	      	/* Gaussian quadrature point locations  */
  double xi[DIM];             /* Local element coordinates of Gauss point. */
  double x_dot[MAX_PDIM];
  double wt, weight;
  double  res, jac;
  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
 

  
  VARIABLE_DESCRIPTION_STRUCT *vd;
  
  int num_dofs;

  /* Determine if Lagrange multiplier unknowns are in augmenting conditions */
  if ( Do_Overlap && augc[nAC-1].lm_eb == augc[nAC-1].solid_eb ) ac_lm = 1;
  else if ( Do_Overlap && augc[nAC-1].lm_eb == augc[nAC-1].fluid_eb ) ac_lm = 2;
  else ac_lm = 0;
  
  /* Determine if this is the augmenting condition pass */
  pass = ( (Do_Overlap && oAC >= 0) ? 2 : 1);

  /***************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)        */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                      */
  /***************************************************************************/

  /* Find out the number of surface quadrature points 
     -this is assumed independent of the surface */

  /****************************************/
  /* N.B. Reconcile all PRS queries below */
  /****************************************/

  /* Surface integration over element */
  
  id_side = elem_side_bc->id_side;

  if ( Do_Overlap && ls->CrossMeshQuadPoints > 0 ) ip_total = ls->CrossMeshQuadPoints;
  else ip_total = elem_info(NQUAD_SURF, ielem_type);

  for (ip = 0; ip < ip_total; ip++) {

    if ( Do_Overlap && ls->CrossMeshQuadPoints > 0 )
      {
        find_subsurf_st_wt(ip, ls->CrossMeshQuadPoints,
                           ielem_type, id_side, pd->Num_Dim, xi, &wt );
      }
    else
      {
        /* find the quadrature point locations (s, t) for current ip */
        find_surf_st(ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);

        /* find the quadrature weight for current surface ip */
        wt = Gq_surf_weight(ip, ielem_type);
      }

    err = load_basis_functions(xi, bfd);
    EH(err, "problem from load_basis_functions");
    
    err = beer_belly();
    EH( err, "beer_belly");
    
    /* precalculate variables at  current integration pt.*/
    err = load_fv();   /*PRS: not sure I need this for getting fluid stresses*/
    EH( err, "load_fv");

    err = load_bf_grad();  /*PRS: DITTO */
    EH( err, "load_bf_grad");

    err = load_bf_mesh_derivs(); 
    EH( err, "load_bf_mesh_derivs");
      
    /* calculate the determinant of the surface jacobian and the normal to 
     * the surface all at one time */
    surface_determinant_and_normal (ielem, iconnect_ptr, num_local_nodes, 
				    ielem_dim - 1,  
				    (int) elem_side_bc->id_side,
				    (int) elem_side_bc->num_nodes_on_side,
				    (elem_side_bc->local_elem_node_id) );
            
    if (ielem_dim !=3) {
      calc_surf_tangent (ielem, iconnect_ptr, num_local_nodes, ielem_dim-1,  
			 (int) elem_side_bc->num_nodes_on_side,
			 (elem_side_bc->local_elem_node_id));
    }
    
    /*
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    err = load_fv_grads();  /*PRS: DITTO, NOT SURE I NEED THESE HERE 
			      REMEMBER I  DON'T EVEN HAVE FLUID VELOCITIES
			      AT THIS NODE*/
    EH( err, "load_fv_grads");
    
    err = load_fv_mesh_derivs(1);
    EH( err, "load_fv_mesh_derivs");

    do_LSA_mods(LSA_SURFACE);
    
    /*
     * Load up porous media variables and properties, if needed 
     */
    if (mp->PorousMediaType == POROUS_UNSATURATED || 
	mp->PorousMediaType == POROUS_SATURATED || 
	mp->PorousMediaType == POROUS_TWO_PHASE){
      err = load_porous_properties(); 
      EH( err, "load_porous_properties"); 
    }

    if (TimeIntegration != STEADY && pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
      for(icount=0; icount<ielem_dim; icount++ ) {
	x_dot[icount] = (1 + 2. * theta)
	    * (fv->x[icount] - fv_old->x[icount])/delta_t
	    - 2. * theta  *  fv_dot->x[icount];
	/* calculate surface position for wall repulsion/no penetration condition */
      }

    }  else {
      for (icount=0; icount<ielem_dim; icount++ ) {
	x_dot[icount] = 0.;
      }
    }
    
    do_LSA_mods(LSA_SURFACE);

    for (ibc = 0;
	 (bc_input_id = (int) elem_side_bc-> BC_input_id[ibc]) != -1;
	 ibc++) {
	
	
      /* check to see if this bc is an integrated bc */

      if (BC_Types[bc_input_id].desc->method == bc_application) {
	    
	/* initialize the general function to zero may have more than one entry
	 * for vector conditions like capillary */

	memset(func, 0, DIM*sizeof(double));
	memset(d_func, 0,
	       DIM*(MAX_VARIABLE_TYPES + MAX_CONC)*MDE*sizeof(double) );

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

	switch (BC_Types[bc_input_id].BC_Name) {
	    
	case SOLID_LAGRANGE_MULT_BC:
          /*
           * If you apply this condition to determine the solid lagrange
           * multiplier, you must change the SOLID_FLUID_CONTACT_BC to use the
           * solid-version of the LM, viz. you don't grab it from the fluid. 
           */

	  /* We are at a gauss point on the surface of a structure 
	     and we need to pass the displaced coordinates of this
	     point in physical space to a search mechanism on the
	     background fluid block */

          /* Decide which terms to evaluate */
          apply_AC = FALSE;
          if (pass == 2) apply_AC = TRUE;

          /* New way: jump down to fluid */
          fluid_mesh_elem_id = jump_down_to_fluid(exo, bc_input_id, x, xi_2);
          iconn_ptr = Proc_Connect_Ptr[fluid_mesh_elem_id];

          /* Save the fluid domain variables */
          for (aa=0; aa < ielem_dim; aa++)
            {
              fluid_velocity[aa] = fv->v[aa];
              if (ac_lm)
                {
                  if (apply_AC) lm_fluid[aa] = augc[ioffset+aa].lm_value;
                  dof_l[aa] = 1;
                  phi_l[aa][0] = 1.0;
                }
              else if (pd->v[pg->imtrx][LAGR_MULT1])
                {
                  lm_fluid[aa] = fv->lm[aa];
                  dof_l[aa] = ei[pg->imtrx]->dof[LAGR_MULT1+aa];
                  for (j=0; j<dof_l[aa]; j++)
                    {
                      phi_l[aa][j] = bf[LAGR_MULT1+aa]->phi[j];
                    }
                }

              if (apply_AC)
                {
                  var = VELOCITY1 + aa;
                  dof_q[aa] = ei[pg->imtrx]->dof[var];
                  for (j=0; j<dof_q[aa]; j++)
                    {
                      phi_q[aa][j] = bf[var]->phi[j];
                      nu = lookup_active_dof(var, j, iconn_ptr);
                      nunk[aa][j] = nu;
                    }
                }
            }

	  /* Before we do the rest of this evaluation, restore world since
	     invert_isoparameteric_map changed this */
          setup_shop_at_point( ielem, xi, exo );
	  Lagrange_mult_equation(func, d_func, theta, delta_t, h_elem_avg,
                                 lm_fluid, fluid_velocity, apply_AC,
                                 dof_q, phi_q, x_dot);
	  break;   

	      
	case SOLID_FLUID_CONTACT_BC:

	  /* We are at a gauss point on the surface of a structure 
	     and we need to pass the displaced coordinates of this
	     point in physical space to a search mechanism on the
	     background fluid block */

          /* Decide which terms to evaluate */
          apply_AC = FALSE;
          if (pass == 2) apply_AC = TRUE;

          /* New way: jump down to fluid */
          fluid_mesh_elem_id = jump_down_to_fluid(exo, bc_input_id, x, xi_2);
          iconn_ptr = Proc_Connect_Ptr[fluid_mesh_elem_id];
          
          /* Save the fluid domain variables */
          for (aa=0; aa< ielem_dim; aa++)
            {
              if (ac_lm)
                {
                  lagrange_mult[aa] = augc[ioffset+aa].lm_value;
                }
              else if (pd->v[pg->imtrx][LAGR_MULT1])
                {
                  lagrange_mult[aa] = fv->lm[aa];
                }

              if (apply_AC)
                {
                  var = LAGR_MULT1 + aa;
                  if (pd->e[pg->imtrx][var] && !ac_lm)
                    {
                      dof_l[aa] = ei[pg->imtrx]->dof[var];
                      for (j=0; j<dof_l[aa]; j++)
                        {
                          phi_l[aa][j] = bf[var]->phi[j];
                        }
                    }
                  else
                    {
                      dof_l[aa] = 1;
                      phi_l[aa][0] = 1.0;
                    }
                }
            }

          /* there was a bunch of stuff here you can get back if you
	     need it to evaluate the fluid stress tensor. Check out
	     a version predating 4/28/02 if you revert to this non-LM
	     technique. */

	  /*before we do the rest of this evaluation, restore world since
	    invert_isoparameteric_map changed this */
          setup_shop_at_point(ielem, xi, exo);
          contact_fn_dot_T(func, d_func, delta_t,
                           lagrange_mult, apply_AC, dof_l, phi_l);
        
	  break;         
	case LAGRANGE_NO_SLIP_BC:

	  /* We are at a gauss point on the surface of a structure 
	     and we need to enforce the Lagrange constraint of
	     no slip with the fluid */
          if ( ac_lm == 2)
            {
               printf("For LAGRANGE_NO_SLIP_BC, AC Lagrange Multiplier must be on solid elements.\n");
               printf("Possibly, you want LS_NO_SLIP.\n");
               EH(-1,"LAGRANGE_NO_SLIP_BC error.");
            }
            
          if (pass == 1) ioffset = first_overlap_ac(ielem, elem_side_bc->id_side);
          else if (pass == 2) ioffset = oAC;
          
          /* Decide which terms to evaluate */
          apply_AC = FALSE;
          if (pass == 2) apply_AC = TRUE;
          
          /* New way: jump_down_to_fluid */
          fluid_mesh_elem_id = jump_down_to_fluid(exo, bc_input_id, x, xi_2);
          iconn_ptr = Proc_Connect_Ptr[fluid_mesh_elem_id];

          /* This call just evaluated fv->v for us. Let us store it
             so we don't lose it going back to the base state */
          for (aa = 0; aa < pd->Num_Dim; aa++)
            {
              fluid_velocity[aa] = fv->v[aa];

              if (apply_AC)
                {
                  var = VELOCITY1 + aa;
                  dof_q[aa] = ei[pg->imtrx]->dof[var];
                  for (j = 0; j < dof_q[aa]; j++)
                    {
                      phi_q[aa][j] = bf[var]->phi[j];
                      nu = lookup_active_dof(var, j, iconn_ptr);
                      nunk[aa][j] = nu;
                    }
                }
            }

	  /*before we do the rest of this evaluation, restore world since
	    invert_isoparameteric_map changed this */
          setup_shop_at_point(ielem, xi, exo);
	  solid_kinematic_bc(func, d_func, delta_t, theta, fluid_velocity,
                             apply_AC, dof_q, phi_q, x_dot);
	  break;

	case BAAIJENS_SOLID_FLUID_BC:
          if (ac_lm == 0)
            {
              setup_shop_at_point(ielem, xi, exo);

              for(aa = 0; aa < pd->Num_Dim; aa++)
                {
                  func[aa] += fv->lm[aa];
                }
              if (af->Assemble_Jacobian)
                {
                  for(bb = 0; bb < pd->Num_Dim; bb++)
                    {
                      var = LAGR_MULT1 + bb;

                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
                        {
                          d_func[bb][var][j] += bf[var]->phi[j];
                        }
                    }
                }
            }
          else
            {
              if ( ac_lm == 1 )
                {
                  /* We are at a gauss point on the surface of a structure
                     and we need to enforce the Lagrange constraint variation
                     on the mesh motion equation for the structure.  viz. we
                     can stay home on this one as we need nothing from the
                     fluid mesh */
                  if (pass == 1) ioffset = first_overlap_ac(ielem, elem_side_bc->id_side);
                  else if (pass == 2) ioffset = oAC;
                }
              else if ( ac_lm == 2 )
                {
                  /* We are at a gauss point on the surface of a structure
                     and we need to enforce the Lagrange constraint variation
                     on the mesh motion equation for the structure.  but the
                     lagrange multiplier lives on the fluid mesh */
                 
                  /* New way: jump_down_to_fluid */
                  fluid_mesh_elem_id = jump_down_to_fluid(exo, bc_input_id, x, xi_2);
                  iconn_ptr = Proc_Connect_Ptr[fluid_mesh_elem_id];
                  ioffset = first_overlap_ac(fluid_mesh_elem_id, -1);
                  
                  /* this could be dangerous.
                   * once in a while, a solid element can just slightly cut into a
                   * fluid element, but none of the fluid nodes change sign.
                   * One bad idea, perfect for trying here, is to just continue along
                   * without apply a force at all in this area
                   */
                  /*
                  if (ioffset == -1) EH(-1,"Bad AC index");
                   */
                  setup_shop_at_point(ielem, xi, exo);
                  if (ioffset == -1)
                    {
                      printf("Fluid element %d at point (%g,%g,%g) has no zero crossing\n",
                             fluid_mesh_elem_id,fv->x[0],fv->x[1],fv->x[2]);
                      continue;
                    }
                }
                
              for(aa = 0; aa < pd->Num_Dim; aa++)
                {
                  func[aa] += augc[ioffset+aa].lm_value;
                }
              if (af->Assemble_Jacobian)
                {
                  for(bb = 0; bb < pd->Num_Dim; bb++)
                    {
                      var = LAGR_MULT1 + bb;

                      d_func[bb][var][0] += 1.0;
                    }
                }
            }
	  break;

	} /* end of switch over bc type */
		
	/**********************************************************************/
	/*        LOOP OVER THE LOCAL ELEMENT NODE NUMBER                     */
	/*        FOR NODES THAT LIE ON THE CURRENT SURFACE		      */
	/* this is a loop designed to loop over equations that lie on         */
	/*  current surface                                                   */
	/*    ADD the Boundary condition functions into the Residual          */
	/*         vector and Jacobian Matrix                                 */
	/**********************************************************************/

        num_dofs = (int) elem_side_bc->num_nodes_on_side;
        if ( BC_Types[bc_input_id].BC_Name == LAGRANGE_NO_SLIP_BC ||
             BC_Types[bc_input_id].BC_Name == SOLID_LAGRANGE_MULT_BC ) num_dofs = 1;
                 
	for (i = 0; i < num_dofs; i++)  {
	      
	  /* Find the local element node number for the current node */
	  id = (int) elem_side_bc->local_elem_node_id[i];
	      
	  /* Find the processor node number given the
	   * local element node number,  'i'
	   */
	  I = Proc_Elem_Connect[iconnect_ptr + id];      

	  /*
	   * If BC is capillarity, convert the node-based cfunc to func 
	   PRS: Don't need
	   */
	
          /*
	   * Boundary condition may actually be a vector of
	   * bc's. Loop over that vector here.
	   */
	  for (p = 0; p < BC_Types[bc_input_id].desc->vector; p++) {
	    /* 
	     *   Check to see if this BC on this node is
	     *   applicable (i.e. no other overriding Dirichlet conditions),
	     *   And, find the global unknown number, index_eq, for
	     *   applying this condition, eqn
	     */
            index_eq = bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn,
				    p, &eqn, &matID_apply, &vd);	  
	    if (index_eq >= 0 || 
		(((BC_Types[bc_input_id].BC_Name == LAGRANGE_NO_SLIP_BC) ||
		 (BC_Types[bc_input_id].BC_Name == SOLID_LAGRANGE_MULT_BC)) &&
		 ((ac_lm == 1) || (pd->i[pg->imtrx][R_LAGR_MULT1] == I_P0))) ) {
	      /*
	       * Obtain the first local variable degree of freedom
	       * at the current node, whether or not it actually an
	       * interpolating degree of freedom
	       */
	      
	      if((BC_Types[bc_input_id].BC_Name == LAGRANGE_NO_SLIP_BC ||
		  BC_Types[bc_input_id].BC_Name == SOLID_LAGRANGE_MULT_BC) &&
		 (ac_lm == 1 || pd->i[pg->imtrx][R_LAGR_MULT1] == I_P0) )
		{
		  /* Sorry, but this is special compensation for boundary
		     integral conditions weigted by a piecewise constant
		     basis function.  bc_eqn_index above does not handle this */
		  ldof_eqn = 0;
		  eqn = R_LAGR_MULT1 + p;
		}
	      else
		{	      
		  ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
		}

	      /*
	       *   for weakly integrated boundary conditions,
	       *   weight the function by wt
	       */
	      weight = wt;

	
	      /*
	       * Do processing specific to CROSS_PHASE_DISCONTINUOUS
	       * boundary conditions.
	       */
	      if (BC_Types[bc_input_id].desc->method == CONTACT_SURF)
		{
		  /* Sorry, but this is special compensation for boundary
		     integral conditions weigted by a piecewise constant
		     basis function.  bc_eqn_index above does not handle this */
		  if((BC_Types[bc_input_id].BC_Name == LAGRANGE_NO_SLIP_BC ||
		     BC_Types[bc_input_id].BC_Name == SOLID_LAGRANGE_MULT_BC)&&
		    (ac_lm == 1 || pd->i[pg->imtrx][R_LAGR_MULT1] == I_P0) )
		    {
		      if(bf[R_LAGR_MULT1 + p ] != NULL) 
			{
			  phi_i = bf[R_LAGR_MULT1 + p]->phi[ldof_eqn];
			}
		      else
			{
                          phi_i = 1.;
			}
		      weight *= phi_i;
		      
		    }
		  else	 
		    {
		      /* normal way */
		      phi_i = bf[eqn]->phi[ldof_eqn];
		      weight *= phi_i;
		    }
		}
	      
	      else
		{
		  EH(-1,"Illegal bc method definition");
		}

	      /*
	       * For strong conditions weight the function by BIG_PENALTY
	       */
	      if (BC_Types[bc_input_id].desc->method == STRONG_INT_SURF) {
		weight *= BIG_PENALTY;
	      }
	      /*
	       *   Add in the multiplicative constant for corresponding to
	       *   all boundary conditions, except for certain special
	       *   cases
	       */
	      if (BC_Types[bc_input_id].desc->method == WEAK_INT_SURF
		  && BC_Types[bc_input_id].BC_Name != PSPG_BC
		  && BC_Types[bc_input_id].BC_Name != VELO_SLIP_SOLID_BC  
			) {
		weight *= pd->etm[pg->imtrx][eqn][(LOG2_BOUNDARY)];
	      }

              /*
               * At this point, term transfer will be done by one of two
               * methods. On pass 2, the computed terms will be handled
               * as augmenting conditions. These terms will be summed
               * into the gAC, bAC, and cAC arrays as needed.
               * On pass 1, everything else will be summed into lec to be
               * handled as boundary condition terms.
               *
               * EDW: Since fv_old->x values get reset between passes,
               *      the correct x_dot values cannot be obtained on
               *      pass 2. Therefore, the gAC terms involving x_dot
               *      are calculated and stored on pass 1, then are
               *      transferred to gAC on pass 2.
               */

              /* Handle AC transfer based on BC type (first pass only) */
              if (ac_lm == 1 && pass == 1 && p < 2)
                {
                  iAC = ioffset + p;

                  switch (BC_Types[bc_input_id].BC_Name) {

                    case SOLID_LAGRANGE_MULT_BC:
                    case LAGRANGE_NO_SLIP_BC:
                      eqn = R_LAGR_MULT1 + p;
                      res = func[p];
                      res *= weight * fv->sdet;
                      augc[iAC].lm_resid += res; 
                      break;
                    default:
                      break;
                  }
                }

              /* Handle AC transfer based on BC type (second pass only) */
              if (ac_lm && pass == 2 && p < 2)
                {
                  iAC = ioffset + p;

                  switch (BC_Types[bc_input_id].BC_Name) {

                    case SOLID_LAGRANGE_MULT_BC:
                    case LAGRANGE_NO_SLIP_BC:
                      eqn = R_LAGR_MULT1 + p;
                      var = VELOCITY1 + p;
                      res = func[p];

                      /* Solid LM residuals (assembled on first pass) */
                      gAC[iAC] = augc[iAC].lm_resid;

                      /* Sensitivities to fluid velocity (cross-mesh term) */
                      for (j = 0; j < dof_q[p]; j++)
                        {
                          nu = nunk[p][j];
                          /* DRN this was wrong, right? */
                          /*
                          jac = (d_func[p][var][j] * weight * fv->sdet
                                + res * weight * fv->dsurfdet_dx[p][j]);
                          */
                          jac = d_func[p][var][j] * weight * fv->sdet;
                          if ( nu >= 0 ) cAC[iAC][nu] += jac;
                        }

                      /* Sensitivities to solid displacement */
                      var = MESH_DISPLACEMENT1 + p;
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
                        {
                          nu = lookup_active_dof(var, j, iconnect_ptr);
                          jac = (d_func[p][var][j] * weight * fv->sdet
                                 + res * weight * fv->dsurfdet_dx[p][j]);
                          if ( nu >= 0 ) cAC[iAC][nu] += jac;
                        }
                      break;
                      
                    case SOLID_FLUID_CONTACT_BC:
                      eqn = R_MESH1 + p;
                      var = LAGR_MULT1 + p;
                      res = func[p];

                      /* Solid displacement residuals not affected */

                      /* Sensitivities to fluid Lagrange multiplier 
                         (cross-mesh term) */
                      jac = d_func[p][var][0] * weight * fv->sdet;
                      nu = index_eq;
                      if ( nu < 0 ) EH(-1,"Bad variable index");
                      bAC[iAC][nu] += jac;

                      /* Sensitivities to solid displacement not affected */

                      break;

                    case BAAIJENS_SOLID_FLUID_BC:
                      if (ac_lm)
                        {
                          eqn = R_MESH1 + p;
                          var = LAGR_MULT1 + p;

                          /* Sensitivities to solid Lagrangeo multiplier */
                          jac = d_func[p][var][0] * weight * fv->sdet;
                          nu = index_eq;
                          if ( nu < 0 ) EH(-1,"Bad variable index");
                          bAC[iAC][nu] += jac; 
                        }

                      break;
                  } /* End of BC Type switch block */
              }

              /*
               * Now proceed with loading lec->R and lec->J (first pass only).
	       *
	       * Calculate the position in the local element residual
	       * vector to put the current contribution
	       * -> MASS_FRACTION unknowns get stuck at the end of
	       *    this vector.
	       */
	      if (eqn == R_MASS) {
		ieqn = MAX_PROB_EQN + BC_Types[bc_input_id].species_eq;
	      } else {
		ieqn = upd->ep[pg->imtrx][eqn];
	      }

	      /*
	       *  Add the current contribution to the local element
	       *  residual vector (first pass only)
	       */
	      if (pass == 1 && ieqn != -1 && ldof_eqn != -1) {
		res = weight * fv->sdet * func[p];
		lec->R[ieqn][ldof_eqn] += res;
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
		      
		if (af->Assemble_Jacobian) {

		  /*regrab this from bc_integ when you need to do these terms, although
		    it is not going to be easy as these are nonlocal and so mush involve
		    direct injection into ija/a  */
                  
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
                            jac = weight * func[p] * fv->dsurfdet_dx[q][j];
                            lec->J[ieqn][pvar] [ldof_eqn][j] += jac;
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
			for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          jac = weight * fv->sdet * d_func[p][var][j];
                          lec->J[ieqn][pvar] [ldof_eqn][j] += jac;
			}
		      } /* end of variable exists and BC is sensitive to it */
		    } /* end of var loop over variable types */
                    

		} /* end of NEWTON */
	      }
	    } /* end of if (Res_BC != NULL) - i.e. apply residual at this node */
	  } /* end of loop over equations that this condition applies to */
	}  /* end for (i=0; i< num_nodes_on_side; i++) */
	    
      }  /*End (if INT) (CAPILLARY and KINEMATIC and VELO_NORMAL and VELO_TANGENT . . .) */
    } /*(end for ibc) */
  } /*End for ip = 1,...*/
  
  return(status);
} /* END of routine apply_contact_bc */
/******************************************************************************/
/******************************************************************************/

int
jump_down_to_fluid ( const Exo_DB *exo, /* Ptr to Exodus database */
                     int bcid,          /* Input BC index */
                     double x[],        /* Solution vector */
                     double xi_2[])     /* xi in fluid element */
/*
 * This function locates the current Gauss point from a solid boundary
 * in an overlapping fluid domain element. It evaluates the isoparametric
 * coordinates for that element, and fluid velocity and/or Lagrange multiplier
 * dofs as requested for that point. The fluid element number is returned.
 */
{
  double coordinates[DIM];
  int velo_interp, fluid_mesh_elem_id;
  int eb_index, e_start=0, e_end=0, ifound, j;
  int sign;


/*
 * We are at a gauss point on the surface of a structure and we need to
 * pass the displaced coordinates of this point in physical space to a
 * search mechanism on the background fluid block.
 */

/* First find the element block id you need to search */
  ifound = 0;
  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++)  
    {

/* Identify the materials index with the element block */
      if (!ifound && 
           exo->eb_id[eb_index] == BC_Types[bcid].BC_Data_Int[1])
        {
          ifound = 1;
          e_start = exo->eb_ptr[eb_index];
          e_end   = exo->eb_ptr[eb_index+1];
        }
    }
    
/* Note this last arg is the fluid block id, which is the
 * only one you want to search */
  fluid_mesh_elem_id = find_id_elem(fv->x[0], fv->x[1], fv->x[2],
                                    x, exo, e_start, e_end);
  load_ei(fluid_mesh_elem_id, exo, 0, pg->imtrx);
  
/* make sure this element has velocity defined */
  if ( pd_glob[ei[pg->imtrx]->mn]->i[pg->imtrx][VELOCITY1] <= 0 ) {
    EH( -1, "Element block specified in contact bc does not contain velocity!.");
  }
    
/* Find basis functions associated with velocity variables */
  for(j=0; j< Num_Basis_Functions; j++)
    {
      if (pd_glob[ei[pg->imtrx]->mn]->i[pg->imtrx][VELOCITY1] == bfd[j]->interpolation)
        {
          velo_interp = j;
        }
    }

/* Now determine the local isoparametric  coordinates in the new element */
/* First load up current coordinates */
  for(j=0; j<pd->Num_Dim; j++) coordinates[j]=fv->x[j];

/* Now find the local coordinates on the fluid mesh */
/* PRS fix: this routine core-dumps if the regions don't
            overlap.  Need a graceful exit here */
  for(j=0; j<DIM; j++) xi_2[j] = 0.0;
  invert_isoparametric_map(&fluid_mesh_elem_id, &coordinates[0],
                                  xi_2, exo, x, &velo_interp);

  /* noble hack here, this should be fixed eventually */
  /* currently, it is assumed that the fluid is on positive side of interface */
  sign = ls->Elem_Sign;
  ls->Elem_Sign = 1;
          
/* Now evaluate lm_fluid and/or fluid velocity with calls to load_fv, etc. */
/* For now use kitchen sink approach */
  setup_shop_at_point(fluid_mesh_elem_id, xi_2, exo);

  ls->Elem_Sign = sign;
/*
 * NOTE: Original element must be restored by calling function upon
 *       return from here (with a new call to setup_shop_at_point)!
 */
  /*
  printf("Original coords in solid elem %d = %g,%g, coords in fluid elem %d = %g,%g\n",
         ielem,coordinates[0],coordinates[1],fluid_mesh_elem_id,fv->x[0],fv->x[1]);
  */
  return fluid_mesh_elem_id;
}  /* END of function jump_down_to_fluid() */
/******************************************************************************/

void 
contact_fn_dot_T(double func[DIM],
		 double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                 double dt,
		 double lagrange_mult[3],
                 int cross,
                 int dof_l[DIM],
                 double phi_l[DIM][MDE])
  
/*******************************************************************************
*
*  Function which calculates the normal traction vector n.T, where
*  the stress tensor T is passed into this routine. 
*******************************************************************************/
     
{
/*    TAB certifies that this function conforms to the exo/patran side numbering convention 11/10/98. */
  int a, j;
  int  var, jvar;		/* Degree of freedom counter                 */

  int DeformingMesh;		/* Logical.                                  */

/***************************** EXECUTION BEGINS ******************************/
  /* Based on current element id_side, choose the correct curvature sign */

  DeformingMesh = pd->e[pg->imtrx][R_MESH1]; /* Use to catch bad references to moving */
				  /* mesh which isn't. */

  if(DeformingMesh) 
    {
      for (a = 0; a < pd->Num_Dim; a++)
	{

          /* Residuals */
	  func[a] += lagrange_mult[a]; /* remember to change back */
	}
      
          /* Fluid Lagrange multiplier sensitivities (cross-mesh terms) */
      if (cross)
        {
          for (jvar = 0; jvar < pd->Num_Dim; jvar++)
            {
              var = LAGR_MULT1 + jvar;
              for (j=0; j<dof_l[jvar]; j++)
                {
                  d_func[jvar][var][j] += phi_l[jvar][j];
                }
            }
        }

    } /* if velocity variable is defined.   */
  else 
    {
      EH(-1,"Must have a deforming mesh region for this bc");
    }
} /* END of routine contact_fn_dot_T */

/******************************************************************************/

void 
Lagrange_mult_equation(double func[DIM],
		       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		       double tt,
		       double dt,
		       double h_elem_avg,
		       double lm_fluid[3],
		       double fluid_velocity[DIM],
                       int cross,
                       int dof_u[DIM],
                       double phi_u[DIM][MDE],
		       double x_dot[DIM])  
/*******************************************************************************
*  Function which calculates the normal traction vector n.T,
*  where the stress tensor T is passed into this routine. 
*******************************************************************************/
     
{
/*    TAB certifies that this function conforms to the exo/patran side numbering convention 11/10/98. */
  int a, j;
  int var, jvar;		/* Degree of freedom counter                 */

  int DeformingMesh;		/* Logical.                                  */
  double A, phi_j;
  double penalty = 1.0e+6;
  double jac=0.0;

/***************************** EXECUTION BEGINS ******************************/
  /* Based on current element id_side, choose the correct curvature sign */

  DeformingMesh = pd->e[pg->imtrx][R_MESH1]; /* Use to catch bad references to moving */
				  /* mesh which isn't. */

  /* h_elem_avg = 0.03; */
  if(DeformingMesh) 
    {
      A = 0.;
      for (a = 0; a < pd->Num_Dim; a++)
	{
	  /* A  += dt*(x_dot[a] - fluid_velocity[a])*(x_dot[a] - fluid_velocity[a])/(h_elem_avg); */
	  A  += (x_dot[a] - fluid_velocity[a])*(x_dot[a] - fluid_velocity[a]);
	}

      /* Residuals */
      for (a = 0; a < pd->Num_Dim; a++)
	{
	  /*func[a] += (fv->lm[a] - lm_fluid[a]*(1.0 - A))*penalty; */
	  func[a] += (x_dot[a] - fluid_velocity[a]) * penalty; 
	}

      if(af->Assemble_Jacobian)
	{

          /* Solid displacement sensitivities */
	  for (jvar = 0; jvar < pd->Num_Dim; jvar++)
	    {
	      var = MESH_DISPLACEMENT1+jvar;
	      for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
		{
                  phi_j = bf[var]->phi[j];

                  jac = penalty * (1. - 2.*tt)*phi_j/dt;
                  d_func[jvar][var][j] += jac;
		}
	    }

          /* Fluid Lagrange multiplier sensitivities (disabled for now) */
	  for (jvar = 0; jvar < pd->Num_Dim; jvar++)
	    {
	      var = LAGR_MULT1 + jvar;
	      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)
		{
		  /*d_func[jvar][var][j] += bf[var]->phi[j] * 1.e6;*/
		}
	    }
	}
      
      /* Fluid velocity sensitivities (cross-mesh terms) */
      if (cross)
        {
          for (jvar = 0; jvar < pd->Num_Dim; jvar++)
            {
              var = VELOCITY1 + jvar;
              for (j = 0; j < dof_u[jvar]; j++)
                {
                  jac = phi_u[jvar][j] * penalty;
                  d_func[jvar][var][j] -= jac;
                }
            }
        }

    } /* if velocity variable is defined.   */
  else 
    {
      EH(-1,"Must have a deforming mesh region for this bc");
    }
} /* END of routine Lagrange_mult_equation */

/******************************************************************************/

void 
solid_kinematic_bc(double func[DIM],
		   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		   double dt,
		   double tt,
		   double fluid_velocity[DIM],
                   int cross,
                   int dof_u[DIM],
                   double phi_u[DIM][MDE],
		   double x_dot[DIM])
  
/*******************************************************************************
*  Function which calculates the lagrange condition n.(v_fluid - v_solid)
*  for the unknown lagrange multiplier, from the solid side. 
*******************************************************************************/
     
{
  int a, j;
  int var, jvar;		/* Degree of freedom counter                 */
  double phi_j;
  double penalty = 1.0e+0;
  double res=0.0, jac=0.0;

/***************************** EXECUTION BEGINS ******************************/
  /* Based on current element id_side, choose the correct curvature sign */

  /* Residuals */
  for (a = 0; a < pd->Num_Dim; a++)
    {
      res = x_dot[a] - fluid_velocity[a];
      res *= penalty;
      func[a] += res;
    }

  /* Solid displacement sensitivities */
  for (jvar = 0; jvar < pd->Num_Dim; jvar++)
    {
      var = MESH_DISPLACEMENT1+jvar;
      for (j=0; j<ei[pg->imtrx]->dof[var]; j++)
	{
	  phi_j = bf[var]->phi[j];
          jac = (1. + 2.*tt)*phi_j/dt;
          jac *= penalty;
          d_func[jvar][var][j] += jac;
		
	}
    }
      
  /* Fluid velocity sensitivities (cross-mesh terms) */
  if (cross)
    {
      for (jvar = 0; jvar < pd->Num_Dim; jvar++)
        {
          var = VELOCITY1 + jvar;
          for (j = 0; j < dof_u[jvar]; j++)
            {
              jac = phi_u[jvar][j];
              jac *= penalty;
              d_func[jvar][var][j] -= jac;
            }
        }
    }

} /* END of routine solid_kinematic_bc */
/*****************************************************************************/
int
apply_embedded_bc (
          int ielem,      /* element number */
	  double x[],     /* Solution vector for the current processor    */
          double dt,      /* current time step size                       */
          double theta,	  /* parameter (0 to 1) to vary time integration
                                   *  ( implicit - 0 to explicit - 1)     */
          double time_value,
	  const PG_DATA *pg_data,
          int oAC,        /* Flag indicating calling function */
          double *gAC,    /* Augmenting condition arrays */
          double **bAC,
          double **cAC,
	  Exo_DB *exo)

    /***************************************************************************
     *
     * apply_embedded_bc():
     *
     *    Calculate the local element contributions for "boundary" conditions
     *    that cut through the mesh (along level set surface).
     *
     **************************************************************************/
{
  double xi[DIM];
  int err, status = 0;
  int ip;
  double wt;
  int ipass, num_elem_passes = 1;
  int old_elem_sign;
  int  i, a;
  /* finite difference for path dependencies */
  double fplus, fminus, fmax, fmin;
  double h_elem;
  dbl ls_F[MDE], ad_wt[10];	/*  adaptive integration weights  */
  int dim = pd->Num_Dim;                  /* Number of dimensions. */

  if ( !ls->elem_overlap_state ) return(0);
    
  /* calc dphi */
  h_elem = 0.;
  for ( a=0; pg_data !=NULL && a<dim; a++)
    {
      h_elem += pg_data->hsquared[a];
    }
  /* This is the average value of h**2 in the element */
  h_elem = h_elem/ ((double )dim);

  /* This is the size of the element */
  h_elem = sqrt(h_elem);
  
  fmin = 0.;
  fmax = 0.;
  for ( i = 0; i< ei[pg->imtrx]->dof[LS]; i++ )
    {
      if ( *esp->F[i] < fmin ) fmin = *esp->F[i];
      if ( *esp->F[i] > fmax ) fmax = *esp->F[i];
    }
    
  fplus = 0.8 * fmax;
  if ( fplus > FD_FACTOR * h_elem ) fplus = FD_FACTOR * h_elem;
  fminus = 0.8 * fmin;
  if ( fminus < -FD_FACTOR * h_elem ) fminus = -FD_FACTOR * h_elem;
  
#if 0
  if ( fplus > -fminus )
    {
      fplus = -fminus;
      if ( fplus < 1.e-8 * h_elem ) fplus = 1.e-8 * h_elem;
    }
  else
    {
      fminus = -fplus;
      if ( fminus > -1.e-8 * h_elem ) fminus = -1.e-8 * h_elem;
    }
#endif

#if 0
  DPRINTF(stderr,"FINITEDIFF: ielem+1=%d, fmin=%g, fmax=%g, fminus=%g, fplus=%g\n",ielem+1,fmin,fmax,fminus,fplus);
#endif

  old_elem_sign = ls->Elem_Sign;
  if ( xfem != NULL )
    {
      num_elem_passes = 2;
    }
    
  /* on sharp interface */
  ls->on_sharp_surf = 1;

  Subgrid_Int.active = FALSE; /* this is for the volume */

#ifdef COUPLED_FILL
  /* finite difference method for forming path dependencies of source terms */
  if ( ls->SubElemIntegration && !ls->Ignore_F_deps )
    {
      ls->CalcSurfDependencies = TRUE;

      /* LS = fminus */
      Subgrid_Int.ip_total = get_subelement_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, fminus, -1, 0 );

      for ( ip = 0; ip < Subgrid_Int.ip_total; ip++ )
        {
          xi[0] = Subgrid_Int.s[ip][0];
          xi[1] = Subgrid_Int.s[ip][1];
          xi[2] = Subgrid_Int.s[ip][2];
          wt = Subgrid_Int.wt[ip];

          for ( ipass = 0; ipass < num_elem_passes; ipass++ )
            {
              if ( num_elem_passes == 2 ) ls->Elem_Sign = -1 + 2*ipass;

              setup_shop_at_point(ielem, xi, exo);

              if ( fv->h3 == 0. )
                {
                  fv->h3 = 1.;
                  fv->sdet = 1.;
                  fv->wt = wt/(fplus - fminus);
                }
              else
                {
                  fv->sdet = 1.;
                  fv->wt = wt/fv->h3/(fplus - fminus);
                }

              assemble_embedded_bc( ielem, x, dt, theta, time_value,
				    oAC, gAC, bAC, cAC, exo,
				    xi );

            }
        } /* loop over gauss points */

     /* LS = fplus */
     Subgrid_Int.ip_total = get_subelement_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, fplus, -1, 0 );


      for ( ip = 0; ip < Subgrid_Int.ip_total; ip++ )
        {
          xi[0] = Subgrid_Int.s[ip][0];
          xi[1] = Subgrid_Int.s[ip][1];
          xi[2] = Subgrid_Int.s[ip][2];
          wt = Subgrid_Int.wt[ip];

          for ( ipass = 0; ipass < num_elem_passes; ipass++ )
            {
              if ( num_elem_passes == 2 ) ls->Elem_Sign = -1 + 2*ipass;

              setup_shop_at_point(ielem, xi, exo);

              if ( fv->h3 == 0. )
                {
                  fv->h3 = 1.;
                  fv->sdet = 1.;
                  fv->wt = -wt/(fplus - fminus);
                }
              else
                {
                  fv->sdet = 1.;
                  fv->wt = -wt/fv->h3/(fplus - fminus);
                }

              assemble_embedded_bc( ielem, x, dt, theta, time_value,
                                    oAC, gAC, bAC, cAC, exo,
				    xi );

            }
        } /* loop over gauss points */
	
      ls->CalcSurfDependencies = FALSE;
    }
#endif

  /* now form the source terms themselves. Note that we want to do subelement for overset grid cases, which can only
     be using the phase-field struct */
  if ( ls->SubElemIntegration || ls->Evolution == LS_EVOLVE_SLAVE)
    {
      Subgrid_Int.ip_total = get_subelement_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, &Subgrid_Int.ip_sign, 0., -1, 0 );
    }
  else if (ls->AdaptIntegration )
    {
      int i;
      Subgrid_Int.ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
      for(i=0;i<ei[pg->imtrx]->num_local_nodes;i++)	{ls_F[i]=*esp->F[i];}
      i = adaptive_weight( ad_wt, Subgrid_Int.ip_total, ei[pg->imtrx]->ielem_dim, ls_F, 
                           ls->Length_Scale, 3, ei[pg->imtrx]->ielem_type);
      WH(i, "problem with adaptive weight routine");
    }
  else
    {
      Subgrid_Int.ip_total = get_facet_integration_pts ( &Subgrid_Int.s, &Subgrid_Int.wt, exo );
    }

  for ( ip = 0; ip < Subgrid_Int.ip_total; ip++ )
    {

      if( ls->AdaptIntegration )
	{
	  find_stu(ip, ei[pg->imtrx]->ielem_type, &(xi[0]), &(xi[1]), &(xi[2]));
	  wt = ad_wt[Subgrid_Int.ip_total-1-ip];
	}
      else
	{
	  xi[0] = Subgrid_Int.s[ip][0];
	  xi[1] = Subgrid_Int.s[ip][1];
	  xi[2] = Subgrid_Int.s[ip][2];
	  wt = Subgrid_Int.wt[ip];
	}
        
#define DEBUG_LS_INTEGRATION 0
#if DEBUG_LS_INTEGRATION
      ls->Surface_Area += wt * fv->h3;
#endif

      for ( ipass = 0; ipass < num_elem_passes; ipass++ )
        {
          if ( num_elem_passes == 2 ) ls->Elem_Sign = -1 + 2*ipass;

          setup_shop_at_point(ielem, xi, exo);

          if ( fv->h3 == 0. )
            {
              fv->h3 = 1.;
              fv->sdet = 1.;
              fv->wt = wt;
            }
          else
            {
              fv->sdet = 1.;
              fv->wt = wt/fv->h3;
            }
          assemble_embedded_bc( ielem, x, dt, theta, time_value,
                                oAC, gAC, bAC, cAC, exo,
				xi );

          /* equation path dependence terms, but ONLY if calling from volume assembly matrix_fill*/
#ifdef COUPLED_FILL
	  if (oAC != 0)
	    {
	      if( pd->e[pg->imtrx][R_FILL] && ipass == 0)
		{
		  if(  tran->Fill_Equation == FILL_EQN_EIKONAL )
		    {
	          assemble_fill_path_dependence();
		    }
		}
	      if ( !ls->Ignore_F_deps )
		{
		  if( pd->e[pg->imtrx][R_FILL] && ipass == 0)
		    {
		      if(  tran->Fill_Equation == FILL_EQN_EXT_V )
			{
			  err = assemble_extension_velocity_path_dependence();
			}
		    }
		  if( pd->e[pg->imtrx][R_MOMENTUM1] && !ls->AdaptIntegration )
		    {
		      err = assemble_momentum_path_dependence(time_value, theta, dt, pg_data);
		      EH( err, "assemble_momentum_path_dependence");
		    }
		  if( pd->e[pg->imtrx][R_PRESSURE] && !ls->AdaptIntegration )
		    {
		      err = assemble_continuity_path_dependence(
			                     time_value, theta, dt,
					     pg_data);
		      EH( err, "assemble_continuity_path_dependence");
		    }
		  if( pd->e[pg->imtrx][R_ENERGY] && !ls->AdaptIntegration )
		    {
		      assemble_energy_path_dependence(
				 time_value, theta, dt,
					        pg_data);
		    }
		  if( pd->e[pg->imtrx][R_MASS] && !ls->AdaptIntegration )
		    {
		      assemble_mass_transport_path_dependence(time_value, theta, dt, pg_data->hsquared,
							      pg_data->hhv, pg_data->dhv_dxnode,
							      pg_data->v_avg, pg_data->dv_dnode);
		    }
		}
	    } /* if(oAC != 0) */
#endif
        }

    } /* loop over gauss points */
    
  /* leave the place as tidy as it was before you came */
  ls->Elem_Sign = old_elem_sign;
  ls->on_sharp_surf = 0;

  return(status);
} /* END of routine apply_embedded_bc */
/*******************************************************************************/
int
assemble_embedded_bc (
  int ielem,      /* element number */
  double x[],     /* Solution vector for the current processor    */
  double dt,      /* current time step size                       */
  double theta,	  /* parameter (0 to 1) to vary time integration
                                     *  ( implicit - 0 to explicit - 1)     */
  double time_value,
  int oAC,        /* Flag indicating calling function */
  double *gAC,    /* Augmenting condition arrays */
  double **bAC,
  double **cAC,
  Exo_DB *exo,
  double xi[DIM] )

/***************************************************************************
 *
 * assemble_embedded_bc():
 *
 *    Calculate the local element contributions for "boundary" conditions
 *    that cut through the mesh (along level set surface).
 *
 **************************************************************************/
{
  double res, jac;
  struct LS_Embedded_BC *bcref;
  BOUNDARY_CONDITION_STRUCT *bc;
  /* New vars for overlap AC algorithm */
  int id_side, nu, iAC=0, ioffset = -1;
  int dof_l[DIM], dof_q[DIM], nunk[DIM][MDE];
  double phi_l[DIM][MDE], phi_q[DIM][MDE];
  int err, status = 0;

  int iconn_fptr = exo->elem_ptr[ielem];
  int iconn_sptr = exo->elem_ptr[ielem];
  int i, j, a, v;
  int  ldof;
  double phi_i, phi_j;
  double dsigma_dx[DIM][MDE];

  double lagrange_mult[3] = {0.0, 0.0, 0.0};
  int eqn, peqn, pvar;
  int ac_lm = FALSE;
  int pass;
  double wt = fv->wt;

  /***************************** EXECUTION BEGINS ******************************/

  /* perform any additional setup not done is setup_shop_at_point */
  
  load_lsi( ls->Length_Scale );
  load_lsi_derivs();
  
  /*
   * Load up commonly used physical properties such as density at
   * the current quadrature point using the material state vector. 
   */
  load_properties(mp, time_value);

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
  }

  /* Determine if Lagrange multiplier unknowns are in augmenting conditions */
  if ( Do_Overlap && augc[nAC-1].lm_eb == augc[nAC-1].solid_eb )
    ac_lm = 1;
  else if ( Do_Overlap && augc[nAC-1].lm_eb == augc[nAC-1].fluid_eb )
    ac_lm = 2;
  else
    ac_lm = 0;

  /* Determine if this is the augmenting condition pass */
  pass = ( (Do_Overlap && oAC >= 0) ? 2 : 1);

	if(ls->var == FILL)
		{ bcref = ls->embedded_bc;}
	else 
  		{bcref = pfd->ls[ls->var-PHASE1]->embedded_bc;}

  while (bcref)
    {
      bc = BC_Types + bcref->bc_input_id;

      if ( xfem != NULL && bc->BC_ID != 0 && ls->Elem_Sign != bc->BC_ID )
        {
          bcref = bcref->next;
          continue;
        }
      if ( (ls->var == LS && strcmp(bc->Set_Type, "LS")) ||
           (ls->var == PHASE1 && strcmp(bc->Set_Type, "PF")) )
        {
          bcref = bcref->next;
          continue;
        }

      switch (bc->BC_Name)
        {

        case LS_Q_BC:
          assemble_q_source( bc->BC_Data_Float[0] );
          break;
        case LS_QLASER_BC:
          assemble_qlaser_source( bc->u_BC, time_value );
          break;
        case LS_QVAPOR_BC:
          assemble_qvapor_source( bc->u_BC);
          break;
        case LS_QRAD_BC:
          assemble_qrad_source( bc->BC_Data_Float[0],bc->BC_Data_Float[1],bc->BC_Data_Float[2],bc->BC_Data_Float[3]);
          break;
        case LS_T_BC:
          assemble_t_source( bc->BC_Data_Float[0], time_value );
          break;
        case LS_LATENT_HEAT_BC:
          assemble_ls_latent_heat_source ( bc->BC_Data_Float[0], 
		bc->BC_Data_Float[1], dt, theta, time_value, 
		bcref->bc_input_id, BC_Types );
          break;
        case LS_YFLUX_BC:
          assemble_ls_yflux_source ( bc->BC_Data_Int[0], bc->BC_Data_Float[0], bc->BC_Data_Float[1], dt, theta, time_value,
                                 bcref->bc_input_id, BC_Types );
          break;
        case LS_CONT_T_BC:
          assemble_cont_t_source( xi );
          break;
        case LS_CONT_VEL_BC:
          assemble_cont_vel_source( xi, exo );
          break;
        case LS_CAPILLARY_BC:
          assemble_csf_tensor();
          break;
	case LS_CAP_HYSING_BC:
	  assemble_cap_hysing(dt, bc->BC_Data_Float[0]);
	  break;
	case LS_CAP_DENNER_DIFF_BC:
	  if (pd->gv[R_NORMAL1]) {
	    assemble_cap_denner_diffusion_n(dt, bc->BC_Data_Float[0]);
	  } else {
	    assemble_cap_denner_diffusion(dt, bc->BC_Data_Float[0]);
	  }
	  break;
        case LS_FLOW_PRESSURE_BC:
          assemble_p_source( bc->BC_Data_Float[0], bc->BC_Data_Int[0] );
          break;
        case LS_ACOUSTIC_SOURCE_BC:
          assemble_ars_source( bc->BC_Data_Float[0], bc->BC_Data_Float[1] );
          break;
        case LS_RECOIL_PRESSURE_BC:
          assemble_precoil_source( bc->BC_Data_Float );
          break;
        case LS_U_BC:
        case LS_V_BC:
        case LS_W_BC:
          assemble_uvw_source( bc->desc->equation, bc->BC_Data_Float[0] );
          break;
        case LS_EXTV_FLUID_SIC_BC:
          assemble_interface_extension_velocity_sic( bc->BC_Data_Int[0]  );
          break;
        case LS_EXTV_KINEMATIC_BC:
        case LS_EXTV_KIN_LEAK_BC:
        case LS_EXTV_LATENT_BC:
          assemble_extv_kinematic ( theta, dt, time_value, bcref->bc_input_id, BC_Types );
          break;
        case LS_EIK_KINEMATIC_BC:
        case LS_EIK_KIN_LEAK_BC:
          assemble_eik_kinematic ( theta, dt, time_value, bcref->bc_input_id, BC_Types );
          break;
        case LS_CAP_DIV_N_BC:
          assemble_div_n_source ();
          break;
        case LS_CAP_DIV_S_N_BC:
          assemble_div_s_n_source ();
          break;
        case LS_CAP_CURVE_BC:
          if( pd->gv[R_NORMAL1] )
            assemble_curvature_with_normals_source () ;
          else
            assemble_curvature_source ();
          break;

        case FLUID_SOLID_CONTACT_BC:
          {

            if (ac_lm)
              {
                if ( (ls->Evolution == LS_EVOLVE_SLAVE)
                     && (ls->init_surf_list->start->type == LS_SURF_SS) )
                  {
                    struct LS_Surf *ss_surf;
                    struct LS_Surf_Closest_Point *cp;
                    ss_surf = closest_surf( ls->init_surf_list, x, exo, fv->x );
                    cp = ss_surf->closest_point;
                    if ( cp->elem == -1 )
                      EH(-1,"Invalid element at closest_point");

                    /* Associate correct AC number with this point */
                    iconn_sptr = exo->elem_ptr[cp->elem];
                    id_side = cp->elem_side;
                    ioffset = first_overlap_ac(cp->elem, id_side);
                    if (ioffset == -1)
                      EH(-1,"Bad AC index");

                    for (a=0; a<pd->Num_Dim; a++)
                      {
                        lagrange_mult[a] = augc[ioffset+a].lm_value;
                      }
                  }
                else
                  {
                    EH(-1, "Level set must be slave SS for this BC!");
                  }
              }
            else
              {
                for (a=0; a<pd->Num_Dim; a++)
                  {
                    lagrange_mult[a] = fv->lm[a];
                  }
              }

            for ( a=0; a<ei[pg->imtrx]->ielem_dim; a++)
              {
                eqn = R_MOMENTUM1 + a;
                peqn = upd->ep[pg->imtrx][eqn];

                for(i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
                  {
                    phi_i = bf[eqn]->phi[i];

                    res = fv->wt * phi_i * lagrange_mult[a] * fv->sdet * fv->h3;
                    if (pass == 1)
                      lec->R[peqn][i] -= res;

                    if ( af->Assemble_Jacobian )
                      {
                        v = LAGR_MULT1 + a;
                        pvar = upd->vp[pg->imtrx][v];
                        ldof = ( (ac_lm) ? 1 : ei[pg->imtrx]->dof[v]);

                        for ( j=0; j<ldof; j++)
                          {
                            phi_j = ( (ac_lm) ? 1.0 : bf[v]->phi[j]);
                            jac = fv->wt * phi_i * phi_j * fv->sdet * fv->h3;

                            if (ac_lm && pass == 2)
                              {
                                iAC = ioffset + a;
                                nu = lookup_active_dof(eqn, i, iconn_fptr);
                                if ( nu >= 0 )
                                  bAC[iAC][nu] -= jac;
                              }
                            else if (!ac_lm && pass == 1)
                              {
                                lec->J[peqn][pvar][i][j] -= jac;
                              }
                          }
                      }
                  }
              }
          }
          break;
        case LS_NO_SLIP_BC:
          {
            double x_dot[3] = {0.,0.,0.};

            if(ls->CalcSurfDependencies) 
              { bcref = bcref->next;
                 continue;
              }

            if ( ac_lm == 1)
              {
                printf("For LS_LAGRANGE_NO_SLIP, AC Lagrange Multiplier must be on fluid elements.\n");
                printf("Possibly, you want LAGRANGE_NO_SLIP.\n");
                EH(-1,"LS_LAGRANGE_NO_SLIP error.");
              }
            /* get x_dot:
             * this may come from the motion of a mesh (slave to SS)
             * or from the rs motion (if rs is defined)
             * or zero (otherwise)
             */
            if ( ( ls->Evolution == LS_EVOLVE_SLAVE ) && ( ls->init_surf_list->start->type == LS_SURF_SS ) )
              {
                struct LS_Surf *ss_surf;
                struct LS_Surf_Closest_Point *cp;
                
                ss_surf = closest_surf( ls->init_surf_list, x, exo, fv->x );
                cp = ss_surf->closest_point;
                
                
                /* use kitchen sink approach for now at solids location */
                if ( cp->elem == -1 )
                  EH(-1,"Invalid element at closest_point");
                iconn_sptr = exo->elem_ptr[cp->elem];
                
                setup_shop_at_point(cp->elem, cp->xi, exo);
                
                for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
                  {
                    if (  pd->TimeIntegration != STEADY &&  pd->v[pg->imtrx][MESH_DISPLACEMENT1+a] )
                      {
                        x_dot[a] = (1.+2.*theta) * (fv->x[a] - fv_old->x[a])/dt -
                          2. * theta * fv_dot->x[a];
                      }
                    else
                      {
                        x_dot[a] = 0.;
                      }
                    
                    /* Grab some things for cross-mesh sensitivities */
                    if (pass == 2)
                      {
                        v = MESH_DISPLACEMENT1 + a;
                        dof_q[a] = ei[pg->imtrx]->dof[v];
                        for (j=0; j<dof_q[a]; j++)
                          {
                            phi_q[a][j] = bf[v]->phi[j];
                            nu = lookup_active_dof(v, j, iconn_sptr);
                            nunk[a][j] = nu;
                          }
                      }
                  }
                
                /* Yep, that's all I needed from closest_point.
                 * Now go back to the original element
                 */
                setup_shop_at_point(ielem, xi, exo);
              }
            else if (pd->e[pg->imtrx][R_SOLID1])
              {
                /* need to add x_dot calculation for rs */
              }
            
            if ( ac_lm )
              {
                if (pass == 1)
                  ioffset = first_overlap_ac(ielem, -1);
                else if (pass == 2)
                  ioffset = oAC;
                if (ioffset == -1)
                  EH(-1,"Bad AC index");
              }
            
            for ( a=0; a<ei[pg->imtrx]->ielem_dim; a++)
              {
                eqn = R_LAGR_MULT1 + a;
                
                /* forced to P0 for now */
                /*phi_i = bf[eqn]->phi[i];*/
                
                /* Fluid Lagrange multiplier residuals */
                res = wt * (x_dot[a] - fv->v[a]) * fv->sdet * fv->h3;
                
                if ( ac_lm )
                  {
                    iAC = ioffset + a;
                    if (pass == 1)
                      {
                        augc[iAC].lm_resid += res;
                      }
                    else if (pass == 2)
                      {
                        gAC[iAC] = augc[iAC].lm_resid;
                        /* Sensitivities to fluid velocity */
                        v = VELOCITY1 + a;
                        for ( j=0; j<ei[pg->imtrx]->dof[v]; j++)
                          {
                            phi_j = bf[v]->phi[j];
                            nu = lookup_active_dof(v, j, iconn_fptr);
                            jac = -wt * phi_j * fv->sdet * fv->h3;
                            if ( nu >= 0 )
                              cAC[iAC][nu] += jac;
                          }
                        
                        /* Sensitivities to solid displacement
                           (cross-mesh terms) */
                        v = MESH_DISPLACEMENT1 + a;
                        for (j=0; j<dof_q[a]; j++)
                          {
                            phi_j = phi_q[a][j] * (1.0 + 2.0 * theta) / dt;
                            jac = wt * phi_j * fv->sdet * fv->h3;
                            nu = nunk[a][j];
                            if ( nu >= 0 )
                              cAC[iAC][nu] += jac;
                          }
                      }
                  }
                else
                  {
                    peqn = upd->ep[pg->imtrx][eqn];
                    /* hard coded to P0 for now */
                    lec->R[peqn][0] += res;
                    
                    /* Sensitivities to fluid velocity */
                    v = VELOCITY1 + a;
                    pvar = upd->vp[pg->imtrx][v];
                    for ( j=0; j<ei[pg->imtrx]->dof[v]; j++)
                      {
                        phi_j = bf[v]->phi[j];
                        jac = -wt * phi_j * fv->sdet * fv->h3;
                        lec->J[peqn][pvar][0][j] += jac;
                      }
                  }
              }
          }
          break;
          
        case BAAIJENS_FLUID_SOLID_BC:
          {

            if(ls->CalcSurfDependencies) 
              { bcref = bcref->next;
                 continue;
              }

            if ( ac_lm == 0 )
              {
                /* local lagrange multiplier */
                for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
                  {
                    lagrange_mult[a] = fv->lm[a];
                  }
              }
            else
              {
                if ( ac_lm == 1 )
                  {
                    /* get solid-lagrange multiplier:
                     */
                    struct LS_Surf *ss_surf;
                    struct LS_Surf_Closest_Point *cp;
                    
                    ss_surf = closest_surf( ls->init_surf_list, x, exo, fv->x );
                    cp = ss_surf->closest_point;
                    
                    /* use kitchen sink approach for now at solids location */
                    if ( cp->elem == -1 )
                      EH(-1,"Invalid element at closest_point");
                    
                    /* Associate the correct AC number with this point */
                    id_side = cp->elem_side;
                    ioffset = first_overlap_ac(cp->elem, id_side);
                    if (ioffset == -1)
                      EH(-1,"Bad AC index");
                  }
                else if ( ac_lm == 2)
                  {
                    ioffset = first_overlap_ac(ielem, -1);
                    if (ioffset == -1)
                      EH(-1,"Bad AC index");
                  }
                
                for (a=0; a<ei[pg->imtrx]->ielem_dim; a++)
                  {
                    iAC = ioffset + a;
                    lagrange_mult[a] = augc[iAC].lm_value;
                    dof_l[a] = 1;
                    phi_l[a][0] = 1.0;
                  }
              }
            
            for ( a=0; a<ei[pg->imtrx]->ielem_dim; a++)
              {
                eqn = R_MOMENTUM1 + a;
                peqn = upd->ep[pg->imtrx][eqn];
                
                for(i=0; i<ei[pg->imtrx]->dof[eqn]; i++)
                  {
                    phi_i = bf[eqn]->phi[i];
                    
                    res = -wt * phi_i * lagrange_mult[a] * fv->sdet * fv->h3;
                    if (pass == 1)
                      lec->R[peqn][i] += res;
                    
                    if ( ac_lm == 0 && af->Assemble_Jacobian )
                      {
                        v = LAGR_MULT1 + a;
                        pvar = upd->vp[pg->imtrx][v];
                        
                        for (j = 0; j < ei[pg->imtrx]->dof[v]; j++)
                          {
                            phi_j = bf[v]->phi[j];
                            
                            jac = -wt * phi_i * phi_j * fv->sdet * fv->h3;
                            lec->J[peqn][pvar][i][j] += jac;
                          }
                      }
                    else if (pass == 2)
                      {
                        /* Sensitivities to solid Lagrange multiplier
                           (cross-mesh term) */
                        v = LAGR_MULT1 + a;
                        iAC = ioffset + a;
                        nu = lookup_active_dof(eqn, i, iconn_fptr);
                        
                        for (j = 0; j < dof_l[a]; j++)
                          {
                            jac = -wt * phi_i * phi_l[a][j] * fv->sdet * fv->h3;
                            if ( nu >= 0 )
                              bAC[iAC][nu] += jac;
                          }
                      }
                  }
              }
          }
          break;
      
        }
      bcref = bcref->next;
    } /*while bcref */
    

  return(status);
} /* END of routine apply_embedded_bc */

/*******************************************************************************/

int
apply_embedded_colloc_bc ( int ielem,      /* element number */
                           double x[],     /* Solution vector for the current processor    */
                           double dt,      /* current time step size                       */
                           double theta,   /* parameter (0 to 1) to vary time integration
                                            *  ( implicit - 0 to explicit - 1)     */
                           double time_value,
                           Exo_DB *exo,
                           Dpi *dpi )
{
  struct LS_Embedded_BC *bcref;
  BOUNDARY_CONDITION_STRUCT *bc;

  int a, i, j, idof_from, idof_to;
  int is_neg;
  int eqn, var, peqn, pvar;
  int status = 0;
  int ln;

/***************************** EXECUTION BEGINS ******************************/

  /* Bail out fast if there's nothing to do */

  if ( (ls == NULL) || (ls->embedded_bc == NULL) ||
       (xfem == NULL) || (xfem->elem_state == 0) ) return(0);


  bc = BC_Types + ls->embedded_bc->bc_input_id;

  /* Bail out for the special case of LS_CAPILLARY the only embedded_bc.  Don't need to process
   * anything here....
   */

  /*
  if ( ( bc->BC_Name == LS_CAPILLARY_BC || bc->BC_Name == LS_CAP_CURVE_BC )  &&
       (ls->embedded_bc->next == NULL ) ) return (0 ) ;
  */

  /* loop over embedded bc's */

  bcref = ls->embedded_bc;

  while (bcref)
    {
      bc = BC_Types + bcref->bc_input_id;
        
      switch (bc->BC_Name) {
        case LS_U_BC:
        case LS_V_BC:
        case LS_W_BC:
	case LS_T_BC:
          {
            eqn = bc->desc->equation;
            peqn = upd->ep[pg->imtrx][eqn];

            if ( pd->i[pg->imtrx][eqn] == I_Q1_G ||
	         pd->i[pg->imtrx][eqn] == I_Q2_G ||
		 pd->i[pg->imtrx][eqn] == I_Q1_GP ||
	         pd->i[pg->imtrx][eqn] == I_Q2_GP ||
		 pd->i[pg->imtrx][eqn] == I_Q1_GN ||
	         pd->i[pg->imtrx][eqn] == I_Q2_GN ||
		 pd->i[pg->imtrx][eqn] == I_Q1_XG ||
	         pd->i[pg->imtrx][eqn] == I_Q2_XG ||
		 pd->i[pg->imtrx][eqn] == I_Q1_XV ||
		 pd->i[pg->imtrx][eqn] == I_Q2_XV )
              {
                break; /* nothing needs to be done here */
              }
            else
              {
                EH(-1,"LS_UVW expects XFEM interpolation Q1_XV, Q2_XV, Q1_XG, Q2_XG, Q1_G, or Q2_G\n");
              }
            
            for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i+=2 )
              {
                ln = i/2;
                if ( TRUE )
                  {
                    is_neg = lnn_distance( ln ) < 0.;

                    if ( (is_neg && bc->BC_ID == 1) ||
                         (!is_neg && bc->BC_ID == -1) )
                      {
                        continue;
                      }
		    idof_to = i;
		    idof_from = i+1;

		    if (af->Assemble_Residual)
                      {
                        lec->R[peqn][idof_to] += bc->BC_ID  * lec->R[peqn][idof_from];
                      }

                    if (af->Assemble_Jacobian)
                      {
                        for (var=V_FIRST; var<V_LAST; var++)
                          {
                            pvar = upd->vp[pg->imtrx][var];
                            if (pvar != -1 && (Inter_Mask[pg->imtrx][eqn][var]))
                              {
                                for ( j = 0; j< ei[pg->imtrx]->dof[var]; j++ )
                                  {
                                    lec->J[peqn][pvar][idof_to][j] += bc->BC_ID  * lec->J[peqn][pvar][idof_from][j];
                                  }
                              }
                          }
                      }
                  }
              }
	  }
          break;
        case LS_CONT_FLUX_BC:
          {
            eqn = R_ENERGY;
            peqn = upd->ep[pg->imtrx][eqn];

            if ( pd->i[pg->imtrx][eqn] == I_P0_G ||
		 pd->i[pg->imtrx][eqn] == I_P1_G ||
		 pd->i[pg->imtrx][eqn] == I_Q1_G ||
		 pd->i[pg->imtrx][eqn] == I_Q2_G )
              {
		break;
              }
            else
              {
                EH(-1,"LS_CONT_FLUX expects XFEM interpolation P0_G, P1_G, Q1_G, or Q2_G\n");
              }
            
            for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i+=2 )
              {
                ln = i/2;
                if ( xfem->node_var_state[ln] == 1 )
                  {
                    is_neg = lnn_distance( ln ) < 0.;
                    if ( (is_neg && bc->BC_ID == -1) ||
                         (!is_neg && bc->BC_ID == 1) )
                      {
                        idof_from = i+1;
                        idof_to = i;
                      }
                    else
                      {
                        idof_from = i;
                        idof_to = i+1;
                      }
                    if (af->Assemble_Residual)
                      {
                        lec->R[peqn][idof_to] += lec->R[peqn][idof_from];
                      }

                    if (af->Assemble_Jacobian)
                      {
                        for (var=V_FIRST; var<V_LAST; var++)
                          {
                            pvar = upd->vp[pg->imtrx][var];
                            if (pvar != -1 && (Inter_Mask[pg->imtrx][eqn][var]))
                              {
                                for ( j = 0; j< ei[pg->imtrx]->dof[var]; j++ )
                                  {
                                    lec->J[peqn][pvar][idof_to][j] += lec->J[peqn][pvar][idof_from][j];
                                  }
                              }
                          }
                      }
                  }
              }
	  }
          break;
#if 1      
        case LS_CONT_TRACTION_BC:
          {
            for ( a=0; a<ei[pg->imtrx]->ielem_dim; a++ )
              {
                eqn = R_MOMENTUM1 + a;
                peqn = upd->ep[pg->imtrx][eqn];

                if ( pd->i[pg->imtrx][eqn] == I_P0_G ||
		     pd->i[pg->imtrx][eqn] == I_P1_G ||
		     pd->i[pg->imtrx][eqn] == I_Q1_G ||
		     pd->i[pg->imtrx][eqn] == I_Q2_G )
                  {
		    break;
                  }
                else
                  {
                    EH(-1,"LS_CONT_TRACTION_BC expects XFEM interpolation P0_G, P1_G, Q1_G, or Q2_G\n");
                  }

                for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i+=2 )
                  {
                    ln = i/2;
                    if ( xfem->node_var_state[ln] == 1 )
                      {
                        is_neg = lnn_distance( ln ) < 0.;
                        if ( (is_neg && bc->BC_ID == -1) ||
                             (!is_neg && bc->BC_ID == 1) )
                          {
                            idof_from = i+1;
                            idof_to = i;
                          }
                        else
                          {
                            idof_from = i;
                            idof_to = i+1;
                          }
                        if (af->Assemble_Residual)
                          {
                            lec->R[peqn][idof_to] += lec->R[peqn][idof_from];
                          }

                        if (af->Assemble_Jacobian)
                          {
                            for (var=V_FIRST; var<V_LAST; var++)
                              {
                                pvar = upd->vp[pg->imtrx][var];
                                if (pvar != -1 && (Inter_Mask[pg->imtrx][eqn][var]))
                                  {
                                    for ( j = 0; j< ei[pg->imtrx]->dof[var]; j++ )
                                      {
                                        lec->J[peqn][pvar][idof_to][j] += lec->J[peqn][pvar][idof_from][j];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
	  }
          break;
#endif
#if 0
        case LS_CONT_TRACTION_BC:
          {
            for ( a=0; a<VIM; a++ )
              {
                eqn = R_MOMENTUM1 + a;
                peqn = upd->ep[pg->imtrx][eqn];

                if ( pd->i[pg->imtrx][eqn] == I_P0_G ||
		     pd->i[pg->imtrx][eqn] == I_P1_G ||
		     pd->i[pg->imtrx][eqn] == I_Q1_G ||
		     pd->i[pg->imtrx][eqn] == I_Q2_G )
                  {
		    break;
                  }
                else
                  {
                    EH(-1,"LS_CONT_TRACTION_BC expects XFEM interpolation P0_G, P1_G, Q1_G, or Q2_G\n");
                  }

                for ( i = 0; i < ei[pg->imtrx]->dof[eqn]; i+=2 )
                  {
                    ln = i/2;
                    if ( xfem->active_node[ln] )
                      {
                        is_neg = lnn_distance( ln ) < 0.;
                        if ( (is_neg && bc->BC_ID == -1) ||
                             (!is_neg && bc->BC_ID == 1) )
                          {
                            idof_from = i+1;
                            idof_to = i;
                          }
                        else
                          {
                            idof_from = i;
                            idof_to = i+1;
                          }
                        if (af->Assemble_Residual)
                          {
                            lec->R[peqn][idof_to] += lec->R[peqn][idof_from];
                            lec->R[peqn][idof_from] = lec->R[peqn][idof_to];
                          }

                        if (af->Assemble_Jacobian)
                          {
                            for (var=V_FIRST; var<V_LAST; var++)
                              {
                                pvar = upd->vp[pg->imtrx][var];
                                if (pvar != -1 && (Inter_Mask[pg->imtrx][eqn][var]))
                                  {
                                    for ( j = 0; j< ei[pg->imtrx]->dof[var]; j++ )
                                      {
                                        lec->J[peqn][pvar][idof_to][j] += lec->J[peqn][pvar][idof_from][j];
                                        lec->J[peqn][pvar][idof_from][j] = lec->J[peqn][pvar][idof_to][j];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
	  }
          break;
#endif
      } /* switch over BC_NAME */

      bcref = bcref->next;
    } /* while (bcref) */

  return(status);
} /* END of routine apply_embedded_colloc_bc */

void
find_subsurf_st_wt( const int ip, const int ip_total,
                    const int ielem_type, const int id_side,
                    const int num_dim,
                    double *xi,
                    double *wt )
{
  double s, t, u;
  int shape = type2shape( ielem_type );
  int linear_type=0, quad_type=0;
  int linear_ip_total, quad_ip_total;

  switch( shape )
  {
  case QUADRILATERAL:
    linear_type = BILINEAR_QUAD;
    quad_type = BIQUAD_QUAD;
    break;

  case HEXAHEDRON:
    linear_type = TRILINEAR_HEX;
    quad_type = TRIQUAD_HEX;
    break;

  case TRIANGLE:
    linear_type = LINEAR_TRI;
    quad_type = QUAD_TRI;
    break;

  case SHELL:
  case TRISHELL:
    EH(-1,"subsurf integration for overset grids and ls subelement integration not supported for SHELLS");
    break;

  case LINE_SEGMENT:
    linear_type = LINEAR_BAR;
    quad_type = QUAD_BAR;
    break;

  default:
    EH(-1, "Unsupported element shape.");
    break;
  }

  linear_ip_total = elem_info(NQUAD_SURF, linear_type);
  quad_ip_total = elem_info(NQUAD_SURF, quad_type);

  if ( ip_total == linear_ip_total )
    {
      find_surf_st(ip, linear_type, id_side, num_dim, xi, &s, &t, &u);
      *wt = Gq_surf_weight(ip, linear_type);
    }
  else if ( ip_total == quad_ip_total )
    {
      find_surf_st(ip, quad_type, id_side, num_dim, xi, &s, &t, &u);
      *wt = Gq_surf_weight(ip, quad_type);
    }
  else
    {
      double xi0[3], xi1[3];
      int nodes_per_side;
      int nodes[3];

      /* this is currently hacked in for 2-D only */
      if ( num_dim != 2 ) EH(-1, "Cross Mesh sub-integration only supported in 2-D\n");

      switch( shape )
      {
      case QUADRILATERAL:
        *wt = 2./ip_total;
        break;

      case TRIANGLE:
        *wt = 1./ip_total;
        break;

      default:
        EH(-1, "Unsupported element shape.");
        break;
      }

      s = (2*ip-ip_total+1.)/ip_total;

      get_side_info( ielem_type, id_side, &nodes_per_side, nodes );
      find_nodal_stu( nodes[0], ielem_type, xi0, xi0+1, xi0+2 );
      find_nodal_stu( nodes[1], ielem_type, xi1, xi1+1, xi1+2 );

      xi[0] = xi0[0] + 0.5*(s + 1.) *(xi1[0]-xi0[0]);
      xi[1] = xi0[1] + 0.5*(s + 1.) *(xi1[1]-xi0[1]);
      xi[2] = 0.;
    }
}

void
find_segment_s_wt(const int iquad,           /* current GQ index */
                  const int nquad,           /* number of quadrature points */
                  double *s,                 /* Gaussian-quadrature point (s) */
                  double *wt )               /* Gaussian-quadrature wt */
{
  double t, u;
  double xi[DIM];
  if ( nquad == 2 )
    {
      find_stu( iquad, LINEAR_BAR, s, &t, &u );
      *wt = Gq_weight( iquad, LINEAR_BAR );
    }
  else if ( nquad == 3 )
    {
      find_stu( iquad, QUAD_BAR, s, &t, &u );
      *wt = Gq_weight( iquad, QUAD_BAR );
    }
  else if ( nquad == 5 )
    {
      find_surf_st(iquad, BIQUAD_QUAD_LS, 0, 2, xi, s, &t, &u);
      *wt = Gq_surf_weight(iquad, BIQUAD_QUAD_LS);
    }
  else
    {
      *s = (2*iquad-nquad+1.)/nquad;
      *wt = 2./nquad;
    }
}

double dof_distance ( int var_type,
                      int lvdof )
{
  int ln = ei[pd->mi[var_type]]->dof_list[var_type][lvdof];
  int F_dof = ei[pd->mi[ls->var]]->ln_to_dof[ls->var][ln];
  
  if ( F_dof != -1 ) {
    if (ls->var == LS)
      return *esp->F[F_dof];
    else
      return *(esp->pF[ls->var-PHASE1][F_dof]);
  }
  
  /* LS must be defined at the node except for discontinuous vars */
  /* This rules out Q1 LS and Q2 velocity !! */
  if ( pd->i[pg->imtrx][var_type] != I_P0 && pd->i[pg->imtrx][var_type] != I_P1 &&
       pd->i[pg->imtrx][var_type] != I_P0_GP && pd->i[pg->imtrx][var_type] != I_P1_GP &&
       pd->i[pg->imtrx][var_type] != I_P0_GN && pd->i[pg->imtrx][var_type] != I_P1_GN &&
       pd->i[pg->imtrx][var_type] != I_P0_XV && pd->i[pg->imtrx][var_type] != I_P1_XV &&
       pd->i[pg->imtrx][var_type] != I_P0_G && pd->i[pg->imtrx][var_type] != I_P1_G )
    {
      EH(-1, "This combination of LS and var interpolations is not supported!");
    }
      
  /* define element var distance to be distance to node 0 */
  F_dof = ei[pd->mi[ls->var]]->ln_to_dof[ls->var][0];
  if ( F_dof < 0 ) {
    EH(-1,"dof_distance expect LS var to be define at local node 0");
  }
  if (ls->var == LS)
    return *esp->F[F_dof];
  else
    return *(esp->pF[ls->var-PHASE1][F_dof]);
}

double lnn_distance ( int ln )
{
  int F_dof = ei[pd->mi[ls->var]]->ln_to_dof[ls->var][ln];
  
  if ( F_dof != -1 ) {
    if (ls->var == LS)
      return *esp->F[F_dof];
    else
      return *(esp->pF[ls->var-PHASE1][F_dof]);
  }
  /* define element var distance to be distance to node 0 */
  F_dof = ei[pd->mi[ls->var]]->ln_to_dof[ls->var][0];
  if ( F_dof < 0 ) {
    EH(-1,"dof_distance expect LS var to be define at local node 0");
  }
  if (ls->var == LS)
    return *esp->F[F_dof];
  else
    return *(esp->pF[ls->var-PHASE1][F_dof]);
}

/* determine distance at global node I */
void
gnn_distance( const int I,
              const double x[],
              const double xold[],
              const double delta_x[],
              double * F,
              double * Fold,
              double * Fprev )
{
  int ie = Index_Solution( I, ls->var, 0, 0, -1, ls->MatrixNum);
  
  if ( ie == -1 )
    {
      /* assume we are dealing with a discontinuous variable (like P0) */
      Exo_DB *exo = EXO_ptr;
      int elem = exo->centroid_list[I];
      
      if ( elem != -1 )
        {
	  /* use local node 0 from element */
	  int I0 = exo->elem_node_list[ exo->elem_node_pntr[elem] + 0 ];
          ie = Index_Solution( I0, ls->var, 0, 0, -1, ls->MatrixNum);
        }
      
      if ( elem == -1 || ie == -1 )
        {
	  EH(-1,"Combination of variable and LS interpolation types not supported.");
        }
    }

  *F = x[ie];
  if (Fold != NULL) *Fold = xold[ie];
  if (Fprev != NULL) *Fprev = x[ie] + damp_factor * var_damp[idv[pg->imtrx][ie][0]] * delta_x[ie];
  
}

int find_2d_side_id( double x0, double x1 )
/*
 * Quickly get a side ID based on xi[0] and xi[1]
 */
{
  double xa0, xa1;
  int id=0;

  xa0 = fabs(x0);
  xa1 = fabs(x1);

  if (xa1 > xa0)
    {
      if (x1 >  0.99) id = 3;
      if (x1 < -0.99) id = 1;
    }
  else
    {
      if (x0 >  0.99) id = 2;
      if (x0 < -0.99) id = 4;
    }

  return id;
}

int find_3d_side_id( double x0, double x1, double x2 )
/*
 * Quickly get a side ID based on xi[0], xi[1], and xi[2]
 */
{
  double xa0, xa1, xa2;
  int id=0;

  xa0 = fabs(x0);
  xa1 = fabs(x1);
  xa2 = fabs(x2);

  if (xa2 > xa0 && xa2 > xa1)
    {
      if (x2 >  0.99) id = 6;
      if (x2 < -0.99) id = 5;
    }
  else if (xa1 > xa0)
    {
      if (x1 >  0.99) id = 3;
      if (x1 < -0.99) id = 1;
    }
  else
    {
      if (x0 >  0.99) id = 2;
      if (x0 < -0.99) id = 4;
    }

  return id;
}

int
first_overlap_ac( int ielem,
                  int side )
/*
 * Search overlap AC's for this solid element number and side ID
 */
{
  int i = 0, iAC = -1;

  while (iAC == -1 && i < nAC)     
    {
      if (augc[i].Type == AC_OVERLAP && augc[i].lm_elem == ielem
                                     && (side == -1 || augc[i].lm_side == side)) iAC = i;
      i++;
    }

  return iAC;
}

int
lookup_active_dof(int var,
                  int j,
                  int eptr)
     /*
      * Find index into solution vector for DOF j of var on current element.
      * Must call load_ei or load_elem_dofptr for this element first.
      */
{
  int node, index, noffset, nu;
  NODE_INFO_STRUCT *nn;
  NODAL_VARS_STRUCT *nv;

  /* Get unknown index */
  node = ei[pg->imtrx]->dof_list[var][j];
  index = Proc_Elem_Connect[eptr+node];
  /* DRN: XFEM is duped by Index_Solution, use gun_list */
  /*
    nu = Index_Solution(index, var, 0, 0, -1, pg->imtrx);
  */
  nu = ei[pg->imtrx]->gun_list[var][j];

  EH(nu, "Bad unknown index!");

  /*
   * If there is an active Dirichlet on this var at this node,
   * set nu to -2 to signal this.
   */
  nn = Nodes[index];
  nv = nn->Nodal_Vars_Info[pg->imtrx];
  noffset = get_nodal_unknown_offset(nv, var, -2, 0, NULL);
  if (nn->DBC[pg->imtrx] != NULL)
    {
      /* dirichlet bc's only apply to natural dofs for xfem */
      if (pd->i[pg->imtrx][var] == I_P0_G ||
	  pd->i[pg->imtrx][var] == I_P1_G ||
	  pd->i[pg->imtrx][var] == I_P0_GP ||
	  pd->i[pg->imtrx][var] == I_P1_GP ||
	  pd->i[pg->imtrx][var] == I_P0_GN ||
	  pd->i[pg->imtrx][var] == I_P1_GN ||
	  pd->i[pg->imtrx][var] == I_P0_XV ||
	  pd->i[pg->imtrx][var] == I_P1_XV ||
	  pd->i[pg->imtrx][var] == I_P1_XG ||
	  pd->i[pg->imtrx][var] == I_Q1_G ||
	  pd->i[pg->imtrx][var] == I_Q2_G ||
	  pd->i[pg->imtrx][var] == I_Q1_GP ||
	  pd->i[pg->imtrx][var] == I_Q2_GP ||
	  pd->i[pg->imtrx][var] == I_Q1_GN ||
	  pd->i[pg->imtrx][var] == I_Q2_GN ||
	  pd->i[pg->imtrx][var] == I_Q1_XV ||
	  pd->i[pg->imtrx][var] == I_Q2_XV ||
	  pd->i[pg->imtrx][var] == I_Q1_XG ||
	  pd->i[pg->imtrx][var] == I_Q2_XG )
        {
          if (j%2 == 0 && nn->DBC[pg->imtrx][noffset] >= 0) nu = -2;
        }
      else
        {
          if (nn->DBC[pg->imtrx][noffset] >= 0) nu = -2;
        }
    }
  
  return nu;
}

void
setup_shop_at_point(int ielem, 
                    double *xi,
                    const Exo_DB * exo)
     /*
      * This routine moves to element, ielem, at parametric coords, xi[],
      * and populates all the usual structures at this point. Note, this
      * overwrites global data.
      *
      *  Output 
      * ------------
      *    ( I suspect, all of the external structures listed in mm_as.h)
      *
      *  struct Element_Indeces* ei
      *  struct Field_Variables *fv
      *  BASIS_FUNCTIONS_STRUCT    bf
      */
{
  int err;
  
  /* Hmmm, I think it is safe to try to save some time and
   * check if this is a new element.
   * If this is not a new element, don't redo dofptrs
   */
  if (ielem != ei[pg->imtrx]->ielem)
    {
      err = load_elem_dofptr(ielem, exo, x_static, x_old_static,
                             xdot_static, xdot_old_static, 0);
      EH(err, "load_elem_dofptr");
      
      err = bf_mp_init(pd);
      EH(err, "bf_mp_init");
    }

  /*
   * This may not be a gauss point, so minimize chance someone will use this
   * wrong by putting in ridiculous value
   */
  fv->wt = 1.e30;

  err = load_basis_functions(xi, bfd);
  EH( err, "problem from load_basis_functions");
      
  err = beer_belly();
  EH( err, "beer_belly");
  
  err = load_fv();
  EH( err, "load_fv");
  
  err = load_bf_grad();
  EH( err, "load_bf_grad");
  
  /*
   *  Just as in the main element assembly, we ensure that the current element
   *  actually has mesh equations associated with it before calculation
   *  d basis function d mesh equations.
   */
  if (ei[pg->imtrx]->deforming_mesh && (pd->e[pg->imtrx][R_MESH1] || pd->v[pg->imtrx][R_MESH1]))
    {
      if (!pd->e[pg->imtrx][R_MESH1]) {
     // printf(" We are here\n");
      }
      err = load_bf_mesh_derivs(); 
      EH( err, "load_bf_mesh_derivs");
    }
    //if (ei[pg->imtrx]->deforming_mesh && pd->e[pg->imtrx][R_MESH1])
    // {
    // err = load_bf_mesh_derivs(); 
    // EH( err, "load_bf_mesh_derivs");
    //}
  /*
   * HKM -> This has to be checked out for shell equations. Does it overwrite
   *        anything that shouldn't have been overwritten. 
   */
  err = load_fv_grads();
  EH( err, "load_fv_grads");
      
   /*
   *  Just as in the main element assembly, we ensure that the current element
   *  actually has mesh equations associated with it before calculation
   *  d value d mesh equations.
   */
  if (ei[pg->imtrx]->deforming_mesh && (pd->e[pg->imtrx][R_MESH1] || pd->v[pg->imtrx][R_MESH1]))
    {
      if (!pd->e[pg->imtrx][R_MESH1]) {
     //	printf(" We are here2\n");
      }
      err = load_fv_mesh_derivs(0);
      EH( err, "load_fv_mesh_derivs");
    }
    //if (ei[pg->imtrx]->deforming_mesh && pd->e[pg->imtrx][R_MESH1])
    // {
    //  err = load_fv_mesh_derivs(0);
    //  EH( err, "load_fv_mesh_derivs");
    // }

  if (mp->PorousMediaType != CONTINUOUS)
    {
      err = load_porous_properties();
      EH( err, "load_porous_properties");
    }

  do_LSA_mods(LSA_VOLUME);
  computeCommonMaterialProps_gp(tran->time_value);
}


double
fv_at_point(double *xi, int var)
{
  int i;
  double *esp, val, phi_i;
  
  for (i = 0, val=0.; i < ei[pg->imtrx]->dof[var]; i++)
    {
      phi_i = newshape( xi, ei[pg->imtrx]->ielem_type, PSI,
                        ei[pg->imtrx]->dof_list[var][i],
                        ei[pg->imtrx]->ielem_shape, pd->i[pg->imtrx][var], i );
      esp = x_static + ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[var][i]];
      val += phi_i * *esp;
    }
  return val;
}

double
fv_old_at_point( double *xi,
	     int var )
{
  int i;
  double *esp, val, phi_i;
  
  for (i = 0, val=0.; i < ei[pg->imtrx]->dof[var]; i++)
    {
      phi_i = newshape( xi, ei[pg->imtrx]->ielem_type, PSI,
                        ei[pg->imtrx]->dof_list[var][i],
                        ei[pg->imtrx]->ielem_shape, pd->i[pg->imtrx][var], i );
      esp = x_old_static + ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[var][i]];
      val += phi_i * *esp;
    }
  return val;
}

double
fv_and_phi_at_point( double *xi,
	             int var,
                     double *phi,
                     int sign )
{
  int i;
  double *esp, val;
  int old_sign;
  old_sign = ls->Elem_Sign;
  ls->Elem_Sign = sign;
  for (i = 0, val=0.; i < ei[pg->imtrx]->dof[var]; i++)
    {
      phi[i] = newshape( xi, ei[pg->imtrx]->ielem_type, PSI,
                        ei[pg->imtrx]->dof_list[var][i],
                        ei[pg->imtrx]->ielem_shape, pd->i[pg->imtrx][var], i );
      esp = x_static + ei[pg->imtrx]->ieqn_ledof[ei[pg->imtrx]->lvdof_to_ledof[var][i]];
      val += phi[i] * *esp;
    }
  ls->Elem_Sign = old_sign;
  return val;
}
