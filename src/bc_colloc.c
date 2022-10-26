/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 * Routines for calculating Boundary Conditions and adding them
 * to matrix_fill
 *
 *
 *$Id: bc_colloc.c,v 5.11 2010-07-21 16:39:26 hkmoffa Exp $
 */

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc/rotate_coordinates.h"
#include "bc_colloc.h"
#include "bc_integ.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "gds/gds_vector.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_porous.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_viscosity.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_node_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "std.h"
#include "user_bc.h"
#include "util/goma_normal.h"

#define GOMA_BC_COLLOC_C

/*******************************************************************************/

int apply_point_colloc_bc(double resid_vector[], /* Residual vector for the current processor     */
                          double delta_t,        /* current time step size                        */
                          double theta,          /* parameter to vary time integration:
                                                  * explicit (theta = 1) -- implicit (theta = 0)  */
                          int ielem,             /* element number                                */
                          int ip_total,          /* total number of gauss points                  */
                          int ielem_type,        /* element type                                  */
                          int num_local_nodes,
                          int ielem_dim,
                          int iconnect_ptr,
                          ELEM_SIDE_BC_STRUCT *elem_side_bc,
                          /* Pointer to an element side boundary
                           * condition structure                           */
                          int num_total_nodes,
                          int local_node_list_fs[], /* dimensioned [MDE]; list to keep track of
                                                     *  nodes at which solid contributions have been
                                                     *  transfered to liquid (fluid-solid
                                                     *  boundaries) */
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
  int ivar, jvar;
  int el1, el2, nf, jk;
  int inode;
  int status = 0;
  int index_eq, matID_apply;
  int bc_input_id;
  int var_flag;
  int dof_map[MDE];
  int n_dof[MAX_VARIABLE_TYPES];
  int n_dofptr[MAX_VARIABLE_TYPES][MDE];
  int doMeshMapping = 0;
  double xi[DIM];
  double x_dot[MAX_PDIM];
  double dsigma_dx[DIM][MDE];
  double phi_j;
  double func;
  double f_time;
  double d_func[MAX_VARIABLE_TYPES + MAX_CONC];
  double kfunc[DIM];
  double d_kfunc[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double time_intermediate = time_value - theta * delta_t;
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

  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {

    /* Find the local element node number for the current node */
    id = (int)elem_side_bc->local_elem_node_id[i];

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

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      err = load_fv();
      GOMA_EH(err, "load_fv");

      err = load_bf_grad();
      GOMA_EH(err, "load_bf_grad");

      err = load_bf_mesh_derivs();
      GOMA_EH(err, "load_bf_mesh_derivs");

      /* calculate the shape functions and their gradients */

      /* calculate the determinant of the surface jacobian  and the normal to
       * the surface all at one time */
      surface_determinant_and_normal(
          ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1, (int)elem_side_bc->id_side,
          (int)elem_side_bc->num_nodes_on_side, (elem_side_bc->local_elem_node_id));

      if (ielem_dim != 3 && ielem_dim == pd->Num_Dim) {
        calc_surf_tangent(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                          (int)elem_side_bc->num_nodes_on_side, (elem_side_bc->local_elem_node_id));
      }

      /*
       * Load up porous media variables and properties, if needed
       */
      if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_TWO_PHASE ||
          mp->PorousMediaType == POROUS_SATURATED) {
        err = load_porous_properties();
        GOMA_EH(err, "load_porous_properties");
      }

      if (TimeIntegration != STEADY) {
        if (mp->SurfaceTensionModel != CONSTANT) {
          load_surface_tension(dsigma_dx);
          if (neg_elem_volume)
            return (status);
        }
      }

      if (TimeIntegration != STEADY) {
        for (icount = 0; icount < ielem_dim; icount++) {
          x_dot[icount] = fv_dot->x[icount];
        }
      } else {
        for (icount = 0; icount < ielem_dim; icount++) {
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

      for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
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
          index_eq =
              bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn, 0, &eqn, &matID_apply, &vd);
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
              fTmelting(&func, d_func, BC_Types[bc_input_id].BC_Data_Float[0]);
              break;

            case SPLINEX_BC:
            case SPLINEY_BC:
            case SPLINEZ_BC:
            case SPLINE_BC:
              fspline(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
              break;

            case SPLINEX_RS_BC:
            case SPLINEY_RS_BC:
            case SPLINEZ_RS_BC:
            case SPLINE_RS_BC:
              fspline_rs(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
              break;

            case UUSER_COLLOC_BC:
              uuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC, id, time_intermediate);
              break;
            case VUSER_COLLOC_BC:
              vuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC, id, time_intermediate);
              break;
            case WUSER_COLLOC_BC:
              wuser_colloc_surf(&func, d_func, BC_Types[bc_input_id].u_BC, id, time_intermediate);
              break;

            case T_USER_BC:
              tuser(&func, d_func, BC_Types[bc_input_id].u_BC, time_intermediate);
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

            case POROUS_LIQ_PRESSURE_FILL_BC: {
              MATRL_PROP_STRUCT *mp_2 = NULL;

              porous_liq_fill(&func, d_func, BC_Types[bc_input_id].BC_Data_Int[0],
                              BC_Types[bc_input_id].BC_Data_Int[1],
                              BC_Types[bc_input_id].BC_Data_Float[0],
                              BC_Types[bc_input_id].BC_Data_Float[1],
                              BC_Types[bc_input_id].BC_Data_Float[2], mp_2);
            } break;

            case PLANEX_BC:
            case PLANEY_BC:
            case PLANEZ_BC:
            case PLANE_BC:
              fplane(ielem_dim, &func, d_func, BC_Types[bc_input_id].BC_Data_Float);
              break;

            case FILLET_BC:
              f_fillet(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
                       BC_Types[bc_input_id].len_u_BC);
              break;

            case DOUBLE_RAD_BC:
              f_double_rad(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
                           BC_Types[bc_input_id].len_u_BC);
              break;

            case FEATURE_ROLLON_BC:
#ifdef FEATURE_ROLLON_PLEASE
              f_feature_rollon(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
                               BC_Types[bc_input_id].len_u_BC, BC_Types[bc_input_id].BC_Data_Int[0],
                               time_intermediate);
#else
              GOMA_EH(-1, "FEATURE_ROLLON_PLEASE define needed and feature_rollon.h - talk to RBS");
#endif
              break;

            case ROLL_FLUID_BC:
              icount = BC_Types[bc_input_id].BC_Data_Int[2];
              xsurf[0] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt + 1];
              xsurf[1] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt + 2];
              xsurf[2] = BC_Types[icount].BC_Data_Float[BC_Types[icount].max_DFlt + 3];
              f_roll_fluid(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
                           BC_Types[bc_input_id].len_u_BC, xsurf);
              break;
            case MOVING_PLANE_BC: {
              double t = time_intermediate;

              fplane(ielem_dim, &func, d_func, BC_Types[bc_input_id].BC_Data_Float);

              func += (BC_Types[bc_input_id].BC_Data_Float[4] * t +
                       BC_Types[bc_input_id].BC_Data_Float[5] * t * t +
                       BC_Types[bc_input_id].BC_Data_Float[6] * t * t * t);
              if (af->Assemble_LSA_Mass_Matrix)
                GOMA_EH(GOMA_ERROR, "LSA is not currently compatible with MOVING_PLANE_BC");
            } break;

            case MOVING_PLANE_ETCH_BC:
              fmesh_etch_bc(&func, d_func, BC_Types[bc_input_id].BC_Data_Int[0], id, x_dot, theta,
                            delta_t);
              break;

            case MESH_CONSTRAINT_BC:
              fmesh_constraint(&func, d_func, bc_input_id);
              break;

            case UVARY_BC:
            case VVARY_BC:
            case WVARY_BC:
              var_flag = BC_Types[bc_input_id].desc->equation;
              fvelocity_profile(var_flag, ielem_dim, BC_Types[bc_input_id].BC_Name, &func, d_func,
                                BC_Types[bc_input_id].u_BC, time_intermediate);
              break;

            case U_PARABOLA_BC:
            case V_PARABOLA_BC:
            case W_PARABOLA_BC:
              var_flag = BC_Types[bc_input_id].desc->equation;
              fvelocity_parabola(var_flag, ielem_dim, BC_Types[bc_input_id].BC_Name, &func, d_func,
                                 BC_Types[bc_input_id].u_BC, time_intermediate,
                                 BC_Types[bc_input_id].len_u_BC);
              break;

            case U_VES11_PARABOLA_BC:
            case U_VES12_PARABOLA_BC:
            case U_VES22_PARABOLA_BC:
            case U_VES13_PARABOLA_BC:
            case U_VES23_PARABOLA_BC:
            case U_VES33_PARABOLA_BC:
            case U_VES11_1_PARABOLA_BC:
            case U_VES12_1_PARABOLA_BC:
            case U_VES22_1_PARABOLA_BC:
            case U_VES13_1_PARABOLA_BC:
            case U_VES23_1_PARABOLA_BC:
            case U_VES33_1_PARABOLA_BC:
            case U_VES11_2_PARABOLA_BC:
            case U_VES12_2_PARABOLA_BC:
            case U_VES22_2_PARABOLA_BC:
            case U_VES13_2_PARABOLA_BC:
            case U_VES23_2_PARABOLA_BC:
            case U_VES33_2_PARABOLA_BC:
            case U_VES11_3_PARABOLA_BC:
            case U_VES12_3_PARABOLA_BC:
            case U_VES22_3_PARABOLA_BC:
            case U_VES13_3_PARABOLA_BC:
            case U_VES23_3_PARABOLA_BC:
            case U_VES33_3_PARABOLA_BC:
            case U_VES11_4_PARABOLA_BC:
            case U_VES12_4_PARABOLA_BC:
            case U_VES22_4_PARABOLA_BC:
            case U_VES13_4_PARABOLA_BC:
            case U_VES23_4_PARABOLA_BC:
            case U_VES33_4_PARABOLA_BC:
            case U_VES11_5_PARABOLA_BC:
            case U_VES12_5_PARABOLA_BC:
            case U_VES22_5_PARABOLA_BC:
            case U_VES13_5_PARABOLA_BC:
            case U_VES23_5_PARABOLA_BC:
            case U_VES33_5_PARABOLA_BC:
            case U_VES11_6_PARABOLA_BC:
            case U_VES12_6_PARABOLA_BC:
            case U_VES22_6_PARABOLA_BC:
            case U_VES13_6_PARABOLA_BC:
            case U_VES23_6_PARABOLA_BC:
            case U_VES33_6_PARABOLA_BC:
            case U_VES11_7_PARABOLA_BC:
            case U_VES12_7_PARABOLA_BC:
            case U_VES22_7_PARABOLA_BC:
            case U_VES13_7_PARABOLA_BC:
            case U_VES23_7_PARABOLA_BC:
            case U_VES33_7_PARABOLA_BC:
              var_flag = BC_Types[bc_input_id].desc->equation;
              f_vestress_parabola(var_flag, ielem_dim, U_PARABOLA_BC, ei[pg->imtrx]->mn, &func,
                                  d_func, BC_Types[bc_input_id].u_BC, time_intermediate,
                                  BC_Types[bc_input_id].len_u_BC);
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
              GOMA_EH(err, "Illegal entry in Generalized Dirichlet Condition ");
              break;

            case GD_TIME_BC:
              err = evaluate_time_func(time_intermediate, &f_time, bc_input_id);
              GOMA_EH(err, "Problems in evaluating time function");
              if (af->Assemble_LSA_Mass_Matrix)
                GOMA_EH(GOMA_ERROR, "LSA is not currently compatible with GD_TIME_BC");
              break;

            case POROUS_PRESSURE_BC:
              porous_pressure(&func, d_func, (int)BC_Types[bc_input_id].BC_Data_Int[0],
                              (int)BC_Types[bc_input_id].BC_Data_Int[1]);
              break;
            case POROUS_PRESSURE_LUB_BC:
              porous_pressure_lub(&func, d_func, elem_side_bc->id_side, xi, exo,
                                  BC_Types[bc_input_id].BC_Data_Float[0]);
              break;

            case LUB_PRESS_HYDROSTATIC_BC:
              lub_press_hydro(&func, d_func, BC_Types[bc_input_id].BC_Data_Float[0],
                              BC_Types[bc_input_id].BC_Data_Float[1],
                              BC_Types[bc_input_id].BC_Data_Float[2],
                              BC_Types[bc_input_id].BC_Data_Float[3]);
              break;

            case LUBP_SH_FP_FLUX_BC:
              put_lub_flux_in_film(id, I, ielem_dim, resid_vector,
                                   (int)BC_Types[bc_input_id].BC_Data_Int[0],
                                   (int)BC_Types[bc_input_id].BC_Data_Int[1], local_node_list_fs);
              func = 0.;
              break;

            case FLUID_SOLID_BC:
            case SOLID_FLUID_BC:
              put_liquid_stress_in_solid(
                  id, I, ielem_dim, resid_vector, (int)BC_Types[bc_input_id].BC_Data_Int[0],
                  (int)BC_Types[bc_input_id].BC_Data_Int[1], local_node_list_fs,
                  BC_Types[bc_input_id].BC_Data_Float[0]);
              func = 0.; /* this boundary condition rearranges values already in res and jac,
                                and does not add anything into the residual */
              break;

            case SOLID_FLUID_RS_BC:
              put_liquid_stress_in_solid_ALE(
                  id, I, ielem_dim, resid_vector, (int)BC_Types[bc_input_id].BC_Data_Int[0],
                  (int)BC_Types[bc_input_id].BC_Data_Int[1], local_node_list_fs,
                  BC_Types[bc_input_id].BC_Data_Float[0]);
              func = 0.; /* this boundary condition rearranges values already in res and jac,
                          * and does not add anything into the residual */
              break;

            case FLUID_SOLID_RS_BC:
              GOMA_EH(GOMA_ERROR, "FLUID_SOLID_RS bc not implemented yet");
              break;

            case SH_FLUID_STRESS_BC: {
              int dof_map_curv[MDE] = {-1};
              int dof_map_tens[MDE] = {-1};
              /* Populate dof_map arrays */
              for (ivar = 0; ivar < ei[pg->imtrx]->dof[VELOCITY1]; ivar++) {
                inode = ei[pg->imtrx]->gnn_list[VELOCITY1][ivar];
                for (jvar = 0; jvar < ei[pg->imtrx]->dof[SHELL_CURVATURE]; jvar++) {
                  if (inode == ei[pg->imtrx]->gnn_list[SHELL_CURVATURE][jvar]) {
                    dof_map_curv[ivar] = jvar;
                  }
                }
                for (jvar = 0; jvar < ei[pg->imtrx]->dof[SHELL_TENSION]; jvar++) {
                  if (inode == ei[pg->imtrx]->gnn_list[SHELL_TENSION][jvar]) {
                    dof_map_tens[ivar] = jvar;
                  }
                }
              }

              /*Note that we send both local node numbers for bulk and shell elements */
              put_fluid_stress_on_shell(id, dof_map_curv[id], dof_map_tens[id], I, ielem_dim,
                                        resid_vector, local_node_list_fs,
                                        BC_Types[bc_input_id].BC_Data_Float[0]);

              func = 0.; /* this boundary condition rearranges values already in res and jac,
                          * and does not add anything into the residual */
            } break;

            case KINEMATIC_COLLOC_BC:
            case VELO_NORM_COLLOC_BC:
              /* initialize the general function to zero may have more than
               * one entry for vector conditions like capillary */
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
              //              if (goma_automatic_rotations.rotation_nodes == NULL) {
              fvelo_normal_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0],
                              contact_flag = FALSE, x_dot, theta, delta_t,
                              (int)BC_Types[bc_input_id].BC_Name, 0, 0, 135.0);
              //              } else {
              //              fvelo_normal_auto_bc(kfunc, d_kfunc,
              //                              BC_Types[bc_input_id].BC_Data_Float[0], contact_flag =
              //                              FALSE, x_dot, theta, delta_t, (int)
              //                              BC_Types[bc_input_id].BC_Name,0,0, 135.0,
              //                              elem_side_bc->id_side, I);
              //              }
              doFullJac = 1;
              func = kfunc[0];
              break;

            case VELO_TANG1_COLLOC_BC:
              GOMA_EH(GOMA_ERROR, "VELO_TANG1_COLLOC_BC not implemented");
              // fzero_velo_tangent_3d(kfunc, d_kfunc, elem_side_bc->id_side, I);
              doFullJac = 1;
              func = kfunc[1];
              for (int var = 0; var < MAX_VARIABLE_TYPES; var++) {
                for (int j = 0; j < MDE; j++) {
                  d_kfunc[0][var][j] = d_kfunc[1][var][j];
                }
              }
              break;
            case VELO_TANG2_COLLOC_BC:
              GOMA_EH(GOMA_ERROR, "VELO_TANG2_COLLOC_BC not implemented");
              // fzero_velo_tangent_3d(kfunc, d_kfunc, elem_side_bc->id_side, I);
              doFullJac = 1;
              func = kfunc[2];
              for (int var = 0; var < MAX_VARIABLE_TYPES; var++) {
                for (int j = 0; j < MDE; j++) {
                  d_kfunc[0][var][j] = d_kfunc[2][var][j];
                }
              }
              break;

            case VELO_NORMAL_LS_COLLOC_BC:
              /* initialize the general function to zero may have more than
               * one entry for vector conditions like capillary */
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
              contact_flag = (ls != NULL);
              fvelo_normal_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0], contact_flag,
                              x_dot, theta, delta_t, (int)BC_Types[bc_input_id].BC_Name,
                              BC_Types[bc_input_id].BC_Data_Float[1],
                              BC_Types[bc_input_id].BC_Data_Float[2],
                              BC_Types[bc_input_id].BC_Data_Float[3]);
              doFullJac = 1;
              func = kfunc[0];
              break;

            case KIN_DISPLACEMENT_COLLOC_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
              f_kinematic_displacement_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Int[0],
                                          BC_Types[bc_input_id].BC_ID, BC_Types[bc_input_id].u_BC,
                                          BC_Types[bc_input_id].len_u_BC);
              func = kfunc[0];
              doFullJac = 1;
              break;
            case TABLE_BC:
              apply_table_bc(&func, d_func, &BC_Types[bc_input_id], time_value);
              if (neg_elem_volume)
                return (status);
              break;

            case SH_GAMMA1_DERIV_SYMM_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
              fgamma1_deriv_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0]);
              func = kfunc[0];
              doFullJac = 1;
              el1 = ei[pg->imtrx]->ielem;
              nf = num_elem_friends[el1];
              if (nf == 0) {
                GOMA_EH(GOMA_ERROR, "no friends");
              };
              el2 = elem_friends[el1][0];
              err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, -1, xi, exo);
              doMeshMapping = 1;
              break;

            case SH_GAMMA2_DERIV_SYMM_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
              fgamma2_deriv_bc(kfunc, d_kfunc, BC_Types[bc_input_id].BC_Data_Float[0]);
              func = kfunc[0];
              doFullJac = 1;
              el1 = ei[pg->imtrx]->ielem;
              nf = num_elem_friends[el1];
              if (nf == 0) {
                GOMA_EH(GOMA_ERROR, "no friends");
              };
              el2 = elem_friends[el1][0];
              err = load_neighbor_var_data(el1, el2, n_dof, dof_map, n_dofptr, -1, xi, exo);
              doMeshMapping = 1;

              break;

            case DVZDR_ZERO_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              memset(d_kfunc, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

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
            ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];

            if (eqn == R_MASS) {
              ieqn = MAX_PROB_EQN + BC_Types[bc_input_id].species_eq;
            } else {
              ieqn = upd->ep[pg->imtrx][eqn];
              if (goma_automatic_rotations.automatic_rotations &&
                  (BC_Types[bc_input_id].desc->rotate != NO_ROT)) {
                int offset = offset_from_rotated_equation(BC_Types[bc_input_id].desc->equation);
                GOMA_EH(offset, "Error translating rotated equation to offset");
                ieqn =
                    equation_index_auto_rotate(elem_side_bc, I, BC_Types[bc_input_id].desc->rotate,
                                               offset, ldof_eqn, &(BC_Types[bc_input_id]));
                GOMA_EH(ieqn, "Could not find index from auto rotate eqn");
              }
            }

            if (ldof_eqn != -1) {
              lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] += penalty * func;
              lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] *= f_time;

              /*
               * add sensitivities into matrix
               *  - find index of sensitivity in matrix
               *     (if variable is not defined at this node,
               *      loop over all dofs in element)
               *  - add into matrix
               */
              if (af->Assemble_Jacobian) {

                for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                  pvar = upd->vp[pg->imtrx][var];
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
                      if (Dolphin[pg->imtrx][I][var] > 0) {
                        if (!doFullJac) {
                          ldof_var = ei[pg->imtrx]->ln_to_first_dof[var][id];
                          if (ldof_var != -1) {
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, ldof_var)] +=
                                penalty * d_func[var];
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, ldof_var)] *= f_time;
                          }
                        } else {

                          if (doMeshMapping &&
                              (var == MESH_DISPLACEMENT1 || var == MESH_DISPLACEMENT2 ||
                               var == MESH_DISPLACEMENT3)) {
                            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                              jk = dof_map[j];
                              lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, jk)] +=
                                  penalty * d_kfunc[0][var][j];
                              lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, jk)] *= f_time;
                            }
                          } else {
                            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                              lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                                  penalty * d_kfunc[0][var][j];
                              lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] *= f_time;
                            }
                          }
                        }
                      } else {
                        /*
                         *   if variable is not defined at this node, loop
                         *   over all dof for this variable in this element
                         */
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          phi_j = bf[var]->phi[j];
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                              penalty * d_func[var] * phi_j;
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] *= f_time;
                        }
                      }
                    } else {
                      for (w = 0; w < pd->Num_Species_Eqn; w++) {
                        pvar = MAX_PROB_VAR + w;
                        if (Dolphin[pg->imtrx][I][var] > 0) {
                          ldof_var = ei[pg->imtrx]->ln_to_first_dof[var][id];
                          if (ldof_var != -1) {
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, ldof_var)] +=
                                penalty * d_func[MAX_VARIABLE_TYPES + w];
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, ldof_var)] *= f_time;
                          }
                        }
                        /* if variable is not defined at this node,
                         * loop over all dof in this element */
                        else {
                          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                            phi_j = bf[var]->phi[j];
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                                penalty * d_func[MAX_VARIABLE_TYPES + w] * phi_j;
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] *= f_time;
                          }
                        }
                      } /* end of loop over species */
                    }   /* end of if MASS_FRACTION */
                  }     /* end of variable exists and condition is sensitive to it */
                }       /* end of loop over variable types */
              }         /* end of NEWTON */
            }           /* if (ldof_eqn != -1) */
          }             /* END of if (Res_BC != NULL), i.e. (index_eqn != -1) */
        }               /* END of if COLLOCATED BC */
        /*****************************************************************************/
      } /* END for (ibc = 0; (int) elem_side_bc->BC_input_id[ibc] != ...*/
        /*****************************************************************************/
    }   /* END if (I < num_owned_nodes) 				      */
        /*****************************************************************************/
  }     /* END for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) */
        /*****************************************************************************/
  return (status);
} /* end of routine apply_point_colloc_bc() */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void moving_plane(int ielem_dim, double *func, double d_func[], dbl *aa, double time) {
  fplane(ielem_dim, func, d_func, aa);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void fplane(int ielem_dim,
            double *func,
            double d_func[], /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
            dbl *aa)         /*  function parameters from data card  */
{
  /**************************** EXECUTION BEGINS *******************************/
  int i;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  *func = (double)aa[3];
  for (i = 0; i < ielem_dim; i++) {
    d_func[MESH_DISPLACEMENT1 + i] = aa[i];
    *func += (double)aa[i] * fv->x[i];
  }

} /* END of routine fplane                                                   */
/*****************************************************************************/

void f_fillet(const int ielem_dim,
              double *func,
              double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
              const double *p,     /*  function parameters from data card  */
              const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double pt[DIM], side_th[2], rad, center[DIM], alpha, theta_mid, theta_avg, tmp;
  double theta, siderad = -1., circ[DIM] = {0., 0., 0.}, dsign, beta = 0, theta_side[2];
  double cham_ang = -1., cham_rotn = 0., gamma = 0.;
  int iside = 0, chamfer = 0, i, dim = 2;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 5) {
    GOMA_EH(GOMA_ERROR, "Need at least 5 parameters for 2D fillet geometry bc!\n");
  }

  pt[0] = p[0];
  pt[1] = p[1];
  side_th[0] = p[2];
  side_th[1] = p[3];
  rad = p[4];
  if (num_const >= 7) {
    siderad = p[5];
    iside = ((int)p[6]);
  }
  if (num_const >= 8) {
    chamfer = ((int)p[7]);
  }
  if (num_const >= 9) {
    cham_ang = p[8];
  }
  if (ielem_dim > dim)
    GOMA_WH(-1, "FILLET_BC: Only z-invariant geometry available for now.\n");

  /**  find center of die face  **/

  alpha = 0.5 * (side_th[1] - side_th[0]);
  theta_avg = 0.5 * (side_th[1] + side_th[0]);
  if (iside) {
    dsign = ((double)(3 - 2 * iside));
    beta = side_th[2 - iside] - dsign * asin((rad - siderad * cos(2. * alpha)) / (rad + siderad));
    for (i = 0; i < dim; i++) {
      circ[i] = pt[i] + dsign * siderad * sin(side_th[iside - 1] - 0.5 * M_PIE * i);
      center[i] = circ[i] + (rad + siderad) * cos(beta - 0.5 * M_PIE * i);
    }
  } else {
    for (i = 0; i < dim; i++) {
      center[i] = pt[i] + (rad / sin(alpha)) * cos(theta_avg - 0.5 * M_PIE * i);
    }
  }

  /**   compute angle of point on curve from arc center **/

  theta = atan2(fv->x[1] - center[1], fv->x[0] - center[0]);
  theta = theta > side_th[1] - 1.5 * M_PIE ? theta : theta + 2 * M_PIE;
  theta_mid = atan2(center[1] - pt[1], center[0] - pt[0]);
  theta_side[0] = side_th[0] - 0.5 * M_PIE;
  theta_side[1] = side_th[1];
  if (iside == 1)
    theta_side[0] = beta - M_PIE;
  if (iside == 2)
    theta_side[1] = beta + 0.5 * M_PIE;
  if (cham_ang > 0) {
    cham_rotn = cham_ang - (theta_mid + 0.5 * M_PIE);
    gamma = atan(-2. * cos(alpha) * sin(cham_rotn) / cos(alpha + cham_rotn));
  }

  /**  use different f depending on theta  **/

  if ((theta_side[0] - gamma) <= theta && theta <= theta_avg) {
    if (iside == 1) {
      *func = SQUARE(fv->x[0] - circ[0]) + SQUARE(fv->x[1] - circ[1]) - SQUARE(siderad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - circ[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - circ[1]);
    } else {
      *func = (fv->x[1] - pt[1]) * cos(side_th[0]) - (fv->x[0] - pt[0]) * sin(side_th[0]);
      d_func[MESH_DISPLACEMENT1] = -sin(side_th[0]);
      d_func[MESH_DISPLACEMENT2] = cos(side_th[0]);
    }

  } else if (theta_avg <= theta && (theta - 0.5 * M_PIE) <= (theta_side[1] + 0 * gamma)) {
    if (iside == 2) {
      *func = SQUARE(fv->x[0] - circ[0]) + SQUARE(fv->x[1] - circ[1]) - SQUARE(siderad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - circ[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - circ[1]);
    } else {
      *func = (fv->x[1] - pt[1]) * cos(side_th[1]) - (fv->x[0] - pt[0]) * sin(side_th[1]);
      d_func[MESH_DISPLACEMENT1] = -sin(side_th[1]);
      d_func[MESH_DISPLACEMENT2] = cos(side_th[1]);
    }

  } else {
    if (chamfer) {
      tmp = theta_mid + cham_rotn;
      *func = (fv->x[1] - center[1]) * sin(tmp) + (fv->x[0] - center[0]) * cos(tmp) +
              rad * sin(alpha - cham_rotn);
      d_func[MESH_DISPLACEMENT1] = cos(tmp);
      d_func[MESH_DISPLACEMENT2] = sin(tmp);
    } else {
      *func = SQUARE(fv->x[0] - center[0]) + SQUARE(fv->x[1] - center[1]) - SQUARE(rad);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - center[0]);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - center[1]);
    }
  }

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_fillet                                                   */
/*****************************************************************************/

void f_double_rad(const int ielem_dim,
                  double *func,
                  double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,     /*  function parameters from data card  */
                  const int num_const) /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double xpt1, ypt1, theta1, rad1, xcen1, ycen1, alpha1;
  double xpt2, ypt2, theta2, rad2, xcen2, ycen2, alpha2;
  double theta1m, theta2m, th1, th2, th2t, curv_mid, rad_curv;
  double beta = 0, xcirc, ycirc, dist1, dist2, dist_mid;
  int is_curved = 0;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 8)
    GOMA_EH(GOMA_ERROR, "Need at least 8 parameters for Double Rad lip geometry bc!\n");

  xpt1 = p[0];
  ypt1 = p[1];
  theta1 = p[2];
  rad1 = p[3];
  xpt2 = p[4];
  ypt2 = p[5];
  theta2 = p[6];
  rad2 = p[7];
  if (num_const >= 8) {
    curv_mid = p[8];
  } else {
    curv_mid = 0.0;
  }

  is_curved = DOUBLE_NONZERO(curv_mid);

  /*  slope of middle line                */

  theta1m = atan2(ypt2 - ypt1, xpt2 - xpt1);
  theta1m = theta1m > theta1 ? theta1m : theta1m + 2 * M_PIE;
  theta2m = atan2(ypt1 - ypt2, xpt1 - xpt2);
  alpha1 = 0.5 * (theta1m - theta1);
  alpha2 = 0.5 * (theta2 - theta2m);

  xcen1 = xpt1 + (rad1 / sin(alpha1)) * cos(theta1 + alpha1);
  ycen1 = ypt1 + (rad1 / sin(alpha1)) * sin(theta1 + alpha1);
  xcen2 = xpt2 + (rad2 / sin(alpha2)) * cos(theta2m + alpha2);
  ycen2 = ypt2 + (rad2 / sin(alpha2)) * sin(theta2m + alpha2);

  if (is_curved) {
    rad_curv = 1. / curv_mid;
    dist_mid = sqrt(SQUARE(xpt1 - xpt2) + SQUARE(ypt1 - ypt2));
    beta = asin(0.5 * dist_mid * curv_mid);
    xcirc = 0.5 * (xpt1 + xpt2) + rad_curv * cos(beta) * sin(theta1m);
    ycirc = 0.5 * (ypt1 + ypt2) - rad_curv * cos(beta) * cos(theta1m);
    /**   Shift fillet centers based on curvature  **/
    /*  Using approximate distance from 90 degree corner for simplicity */

    dist1 = rad1 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad1)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
    dist2 = rad2 + sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid - 2 * rad2)) -
            sqrt((rad_curv - 0.5 * dist_mid) * (rad_curv + 0.5 * dist_mid));
#if 0
     dist1 = 0;  dist2 = 0;
#endif
    xcen1 -= dist1 * cos(theta1);
    ycen1 -= dist1 * sin(theta1);
    xcen2 -= dist2 * cos(theta2);
    ycen2 -= dist2 * sin(theta2);
#if 0
fprintf(stderr,"arc distances %g %g \n",dist1,dist2);
fprintf(stderr,"rads %g %g %g %g\n",rad1, rad2, rad_curv,beta);
fprintf(stderr,"circle %g %g %g %g\n",xcirc,ycirc,xcen2, ycen2);
#endif
  }

  /**   compute angle of point on curve from arc center **/

  th1 = atan2(fv->x[1] - ycen1, fv->x[0] - xcen1);
  th2 = atan2(fv->x[1] - ycen2, fv->x[0] - xcen2);
  th2t = th2 > 0.0 ? th2 : th2 + 2 * M_PIE;

  /**  use different f depending on theta  **/

  if ((theta1 - 0.5 * M_PIE) <= th1 && th1 <= (theta1 + alpha1)) {
    *func = (fv->x[1] - ypt1) * cos(theta1) - (fv->x[0] - xpt1) * sin(theta1);
    d_func[MESH_DISPLACEMENT1] = -sin(theta1);
    d_func[MESH_DISPLACEMENT2] = cos(theta1);
    /*fprintf(stderr,"DR case 1 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2 - alpha2) <= th2t && (th2t - 0.5 * M_PIE) <= theta2) {
    *func = (fv->x[1] - ypt2) * cos(theta2) - (fv->x[0] - xpt2) * sin(theta2);
    d_func[MESH_DISPLACEMENT1] = -sin(theta2);
    d_func[MESH_DISPLACEMENT2] = cos(theta2);
    /*fprintf(stderr,"DR case 2 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if (theta2m <= (th1 + 0.5 * M_PIE - beta) && th1 <= (theta1 - 0.5 * M_PIE)) {
    *func = SQUARE(fv->x[0] - xcen1) + SQUARE(fv->x[1] - ycen1) - SQUARE(rad1);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen1);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen1);
    /*fprintf(stderr,"DR case 3 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else if ((theta2m - 0.5 * M_PIE - beta) <= th2 &&
             ((th1 + 0.5 * M_PIE - beta) <= theta2m ||
              (th1 + 0.5 * M_PIE - beta) <= (theta1m - M_PIE))) {
    if (is_curved) {
      *func = SQUARE(fv->x[0] - xcirc) + SQUARE(fv->x[1] - ycirc) - SQUARE(rad_curv);
      d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcirc);
      d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycirc);
    } else {
      *func = (fv->x[1] - ypt1) * cos(theta1m) - (fv->x[0] - xpt1) * sin(theta1m);
      d_func[MESH_DISPLACEMENT1] = -sin(theta1m);
      d_func[MESH_DISPLACEMENT2] = cos(theta1m);
    }
    /*fprintf(stderr,"DR case 4 %g %g %g %g %g %g %g
     * %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t,theta1m,theta2m);*/
  } else if ((th2t - 0.5 * M_PIE) >= theta2 && (theta2m - 0.5 * M_PIE - beta) >= th2) {
    *func = SQUARE(fv->x[0] - xcen2) + SQUARE(fv->x[1] - ycen2) - SQUARE(rad2);
    d_func[MESH_DISPLACEMENT1] = 2. * (fv->x[0] - xcen2);
    d_func[MESH_DISPLACEMENT2] = 2. * (fv->x[1] - ycen2);
    /*fprintf(stderr,"DR case 5 %g %g %g %g %g %g\n",*func,fv->x[0],fv->x[1],th1,th2,th2t);*/
  } else {
    fprintf(stderr, "Double Rad case not found... %g %g %g %g %g\n", fv->x[0], fv->x[1], th1, th2,
            th2t);
  }

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = 0.0;

} /* END of routine f_double_rad                                             */
/*****************************************************************************/

#ifdef FEATURE_ROLLON_PLEASE
#include "feature_rollon.h"
#endif

void f_roll_fluid(int ielem_dim,
                  double *func,
                  double d_func[],     /* dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
                  const double *p,     /*  function parameters from data card  */
                  const int num_const, /* number of passed parameters   */
                  double *xsurf)       /* number of passed parameters   */
{
  /**************************** EXECUTION BEGINS *******************************/
  double roll_rad;     /* roll radius */
  double origin[3];    /* roll axis origin (x,y,z) */
  double dir_angle[3]; /* axis direction angles */
  double coord[3];     /* current coordinates */
  double axis_pt[3], rad_dir[3], d_dist[3], dist, R, factor, t;
  double omega, v_dir[3], v_roll[3];
  double velo_avg = 0.0, pgrad = 0.;
  double v_solid = 0., res, jac, delta, flow, eps = 1.0e-8, viscinv;
  double jacinv, thick;
  int Pflag = TRUE;
  double pg_factor = 1.0, tang_sgn = 1.0, v_mag = 0.;
  ;

  int j, var;
#if 0
  int jvar,k;
  double dthick_dV, dthick_dP;
#endif

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  if (num_const < 7)
    GOMA_EH(GOMA_ERROR, "Need at least 7 parameters for Roll geometry bc!\n");

  roll_rad = p[0];
  origin[0] = p[1];
  origin[1] = p[2];
  origin[2] = p[3];
  dir_angle[0] = p[4];
  dir_angle[1] = p[5];
  dir_angle[2] = p[6];

  /* calculate distance from interface surface to solid surface for repulsion calculations */

  coord[0] = fv->x[0];
  coord[1] = fv->x[1];
  if (ielem_dim == 3) {
    coord[2] = fv->x[2];
  } else {
    coord[2] = 0.0;
  }

  /*  find intersection of axis with normal plane - i.e., locate point on
   *          axis that intersects plane normal to axis that contains local point. */

  factor = SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]);
  t = (dir_angle[0] * (coord[0] - origin[0]) + dir_angle[1] * (coord[1] - origin[1]) +
       dir_angle[2] * (coord[2] - origin[2])) /
      factor;
  axis_pt[0] = origin[0] + dir_angle[0] * t;
  axis_pt[1] = origin[1] + dir_angle[1] * t;
  axis_pt[2] = origin[2] + dir_angle[2] * t;

  /*  compute radius and radial direction */

  R = sqrt(SQUARE(coord[0] - axis_pt[0]) + SQUARE(coord[1] - axis_pt[1]) +
           SQUARE(coord[2] - axis_pt[2]));
  rad_dir[0] = (coord[0] - axis_pt[0]) / R;
  rad_dir[1] = (coord[1] - axis_pt[1]) / R;
  rad_dir[2] = (coord[2] - axis_pt[2]) / R;
  dist = R - roll_rad;
  d_dist[0] = rad_dir[0] * (1. - SQUARE(dir_angle[0]) / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[0] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[0] / factor);
  d_dist[1] = rad_dir[1] * (1. - SQUARE(dir_angle[1]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[1] / factor) +
              rad_dir[2] * (-dir_angle[2] * dir_angle[1] / factor);
  d_dist[2] = rad_dir[2] * (1. - SQUARE(dir_angle[2]) / factor) +
              rad_dir[0] * (-dir_angle[0] * dir_angle[2] / factor) +
              rad_dir[1] * (-dir_angle[1] * dir_angle[2] / factor);

  if (num_const < 10)
    GOMA_WH(-1, "ROLL_FLUID: Less than 10 parameters - reverting to roll surface!\n");

  omega = p[7];
  /* compute velocity direction as perpendicular to both axis and radial
   *         direction.  Positive direction is determined by right hand rule */

  v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
  v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
  v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

  v_roll[0] = omega * roll_rad * v_dir[0];
  v_roll[1] = omega * roll_rad * v_dir[1];
  v_roll[2] = omega * roll_rad * v_dir[2];

  if (TimeIntegration == TRANSIENT && pd->gv[R_MESH1]) {
    /* Add the mesh motion to the substrate velocity */
    v_roll[0] += fv_dot->x[0];
    v_roll[1] += fv_dot->x[1];
    v_roll[2] += fv_dot->x[2];
  }

  /* quantities specific to FLUID bcs   */

  if (num_const > 8 && p[9] >= 0.0) {
    dist = 0.;
    for (var = 0; var < pd->Num_Dim; var++) {
      /* Uses undeformed node position */
      dist += SQUARE(fv->x0[var] - xsurf[var]);
    }
    dist /= SQUARE(p[10]);
    /*if(dist < 10)fprintf(stderr,"roll_fl %g %g %g\n",fv->x0[0],xsurf[0],dist);
     */

    Pflag = (int)p[11];
    velo_avg = 0.0;
    pgrad = 0.;
    v_mag = 0.;
    for (j = 0; j < pd->Num_Dim; j++) {
      velo_avg += fv->stangent[0][j] * (v_roll[j] + fv->v[j]);
      v_solid += fv->stangent[0][j] * v_roll[j];
      v_mag += SQUARE(v_roll[j]);
      if (Pflag) {
        pgrad += fv->stangent[0][j] * fv->grad_P[j];
      }
    }
    v_mag = sqrt(v_mag);
    tang_sgn = v_solid / v_mag;
    tang_sgn = (double)SGN(v_solid / v_mag);
    velo_avg *= 0.5;
    /* sometimes the tangent/normals flip causing havoc....*/
    if (v_solid < 0) {
      GOMA_WH(-1, "fvelo_slip: normals and tangents have flipped! - try CONTACT_LINE model\n");
      velo_avg *= tang_sgn;
      v_solid *= tang_sgn;
      pgrad *= tang_sgn;
    }

    pg_factor = 1.0;
    if (dist < 10.0) {
      pg_factor = 1.0 - exp(-dist);
    }
    pgrad *= pg_factor;

    flow = MAX(0., p[9] * v_solid);
    viscinv = 1. / p[8];
    thick = flow / velo_avg;
    j = 0;
    do {
      res = -CUBE(thick) * viscinv * pgrad / 12. + thick * velo_avg - flow;
      jac = -0.25 * SQUARE(thick) * viscinv * pgrad + velo_avg;
      jacinv = 1.0 / jac;
      delta = -res * jacinv;
      thick += delta;
      j++;
    } while (fabs(delta) > eps && j < 20);
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
    d_func[MESH_DISPLACEMENT1] = d_dist[0];
    d_func[MESH_DISPLACEMENT2] = d_dist[1];
    d_func[MESH_DISPLACEMENT3] = d_dist[2];
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
            if (pd->v[pg->imtrx][var])
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
    if (pd->v[pg->imtrx][var])
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
  } else {
    *func = dist;
    d_func[MESH_DISPLACEMENT1] = d_dist[0];
    d_func[MESH_DISPLACEMENT2] = d_dist[1];
    d_func[MESH_DISPLACEMENT3] = d_dist[2];
  }

} /* END of routine f_roll_fluid                                                   */
/*****************************************************************************/

/*****************************************************************************/

void fvelocity_profile(int var_flag,
                       int ielem_dim,
                       int velo_condition,
                       double *func,
                       double d_func[], /* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
                       double p[],      /* p[] are the parameters passed in through the input deck*/
                       double time)     /* time at which bc's are evaluated        */
{
  if (af->Assemble_LSA_Mass_Matrix)
    d_func[var_flag] = 0.0;
  else
    d_func[var_flag] = -1.0;
  if (pd->e[pg->imtrx][R_MESH1]) {
    d_func[MESH_DISPLACEMENT1] =
        dvelo_vary_fnc_d1(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
    d_func[MESH_DISPLACEMENT2] =
        dvelo_vary_fnc_d2(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
    if (ielem_dim == 3)
      d_func[MESH_DISPLACEMENT3] =
          dvelo_vary_fnc_d3(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);
  }
  *func = velo_vary_fnc(velo_condition, fv->x[0], fv->x[1], fv->x[2], p, time);

  *func -= fv->v[var_flag - VELOCITY1];

} /* END of routine fvelocity_profile                                        */
/*****************************************************************************/

void fvelocity_parabola(const int var_flag,
                        const int ielem_dim,
                        const int velo_condition,
                        double *func,
                        double d_func[],     /* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
                        const double p[],    /* parameters passed in from the input deck*/
                        const double time,   /* time at which bc's are evaluated   */
                        const int num_const) /* number of passed parameters   */
{
  /*    parabolic velocity profile
   *      p[0] = coordinate1
   *      p[1] = coordinate2
   *      p[2] = flow in positive coordinate direction
   */
  double coord1, coord2, qflow, gap, pre_factor, temp, expon, time_factor;
  double pl_index = 1.0;
  int i;
  coord1 = MIN(p[0], p[1]);
  coord2 = MAX(p[1], p[0]);
  qflow = p[2];
  gap = fabs(coord2 - coord1);
  pre_factor = 6. * qflow / (gap * gap * gap);
  switch (pd->CoordinateSystem) {
  case CARTESIAN:
  case CARTESIAN_2pt5D:
    pre_factor = 6. * qflow / (gap * gap * gap);
    break;
  case CYLINDRICAL:
  case SWIRLING:
    switch (velo_condition) {
    case U_PARABOLA_BC:
      if (coord1 <= DBL_SMALL) {
        pre_factor = 2. * qflow / M_PIE / SQUARE(SQUARE(coord2));
      } else {
        pre_factor = 2. * qflow / M_PIE / (SQUARE(coord2) - SQUARE(coord1)) /
                     (SQUARE(coord2) + SQUARE(coord1) -
                      (SQUARE(coord2) - SQUARE(coord1)) / log(coord2 / coord1));
      }
      break;
    case V_PARABOLA_BC:
      pre_factor = 3. * qflow / M_PIE / (gap * gap * gap);
      break;
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Velo parabola not ready for that Coordinate System yet!\n");
  }

  if (ielem_dim > 2) {
    GOMA_EH(GOMA_ERROR, "Velo parabola not ready for 3D yet!\n");
    return;
  }
  for (i = 0; i < ielem_dim; i++) {
    d_func[MESH_DISPLACEMENT1 + i] = 0.0;
  }

  if (af->Assemble_LSA_Mass_Matrix)
    d_func[var_flag] = 0.0;
  else
    d_func[var_flag] = -1.0;

  if (gap > DBL_SMALL) {
    if (num_const == 3 || p[3] == 1.0) /*  Newtonian solution   */
    {

      switch (pd->CoordinateSystem) {
      case CARTESIAN:
      case CARTESIAN_2pt5D:
        switch (velo_condition) {
        case U_PARABOLA_BC:
          *func = pre_factor * (fv->x[1] - coord1) * (coord2 - fv->x[1]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT2] = pre_factor * (coord1 + coord2 - 2. * fv->x[1]);
          }
          break;
        case V_PARABOLA_BC:
          *func = pre_factor * (fv->x[0] - coord1) * (coord2 - fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = pre_factor * (coord1 + coord2 - 2. * fv->x[0]);
          }
          break;
        case W_PARABOLA_BC:
          *func = pre_factor * (fv->x[0] - coord1) * (coord2 - fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = pre_factor * (coord1 + coord2 - 2. * fv->x[0]);
          }
          break;
        default:
          *func = 0.;
        }
        break;
      case CYLINDRICAL:
      case SWIRLING:
        switch (velo_condition) {
        case U_PARABOLA_BC:
          if (coord1 <= DBL_SMALL) {
            *func = pre_factor * (SQUARE(coord2) - SQUARE(fv->x[1]));
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] = pre_factor * (-2. * fv->x[1]);
            }
          } else {
            *func = pre_factor * (SQUARE(coord1) - SQUARE(fv->x[1]) +
                                  (SQUARE(coord2) - SQUARE(coord1)) *
                                      (log(fv->x[1] / coord1) / log(coord2 / coord1)));
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] =
                  pre_factor * (-2. * fv->x[1] + (SQUARE(coord2) - SQUARE(coord1)) /
                                                     log(coord2 / coord1) / fv->x[1]);
            }
          }
          break;
        case V_PARABOLA_BC:
          *func = pre_factor / fv->x[1] * (fv->x[0] - coord1) * (coord2 - fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = pre_factor / fv->x[1] * (coord1 + coord2 - 2. * fv->x[0]);
            d_func[MESH_DISPLACEMENT2] = -(*func) / fv->x[1];
          }
          break;
        default:
          *func = 0.;
        }
        break;
      }
    } else if (num_const > 3) /*  Power-law  solution   */
    {
      if (p[3] < 0.0) {
        pl_index = gn->nexp;
      } else {
        pl_index = p[3];
      }
      switch (pd->CoordinateSystem) {
      case CARTESIAN:
      case CARTESIAN_2pt5D:
        expon = 1. + 1. / pl_index;
        pre_factor = (2. * pl_index + 1.) / (pl_index + 1.) * qflow / pow(gap, expon + 1.);
        switch (velo_condition) {
        case U_PARABOLA_BC:
          temp = 2 * fv->x[1] - coord1 - coord2;
          *func = pre_factor * (pow(gap, expon) - pow(fabs(temp), expon));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT2] =
                pre_factor * (-2. * SGN(temp) * expon * pow(fabs(temp), 1. / pl_index));
          }
          break;
        case V_PARABOLA_BC:
          temp = 2 * fv->x[0] - coord1 - coord2;
          *func = pre_factor * (pow(gap, expon) - pow(fabs(temp), expon));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                pre_factor * (-2. * SGN(temp) * expon * pow(fabs(temp), 1. / pl_index));
          }
          break;
        case W_PARABOLA_BC:
          temp = 2 * fv->x[0] - coord1 - coord2;
          *func = pre_factor * (pow(gap, expon) - pow(fabs(temp), expon));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                pre_factor * (-2. * SGN(temp) * expon * pow(fabs(temp), 1. / pl_index));
          }
          break;
        default:
          *func = 0.;
        }
        break;
      case CYLINDRICAL:
      case SWIRLING:
        expon = 1. + 1. / pl_index;
        switch (velo_condition) {
        case U_PARABOLA_BC:
          if (coord1 <= DBL_SMALL) {
            pre_factor =
                (3. * pl_index + 1.) / (pl_index + 1.) * qflow / M_PIE / pow(gap, expon + 2.);
            *func = pre_factor * (pow(gap, expon) - pow(fv->x[1], expon));
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] = pre_factor * (-expon * pow(fv->x[1], expon - 1.));
            }
          } else {
            GOMA_EH(GOMA_ERROR, "Power-law annulus not done yet!\n");
          }
          break;
        case V_PARABOLA_BC:
          pre_factor =
              (2. * pl_index + 1.) / (pl_index + 1.) * qflow / M_PIE / pow(gap, expon + 1.);
          temp = 2 * fv->x[0] - coord1 - coord2;
          *func = pre_factor / fv->x[1] * (pow(gap, expon) - pow(fabs(temp), expon));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                pre_factor / fv->x[1] * (-2. * SGN(temp) * expon * pow(fabs(temp), expon - 1.));
            d_func[MESH_DISPLACEMENT2] = -(*func) / fv->x[1];
          }
          break;
        default:
          *func = 0.;
        }
        break;
      }
    }
  } else {
    *func = 0.0;
  }
  /*  Add sinusoidal time-varying pieces   */
  if (num_const > 3 && num_const % 3 == 1) {
    time_factor = 0.;
    for (i = 4; i < num_const; i = i + 3) {
      time_factor += p[i] * sin(p[i + 1] * time + p[i + 2]);
    }
    *func *= (1.0 + time_factor);
    for (i = 0; i < ielem_dim; i++) {
      d_func[MESH_DISPLACEMENT1 + i] *= (1.0 + time_factor);
    }
  }
  /*  Add initial time ramp   */
  else if (num_const == 5) {
    if (p[4] <= 0 || time >= p[4]) {
      time_factor = 1.;
    } else {
      time_factor = time / p[4];
    }
    *func *= time_factor;
    for (i = 0; i < ielem_dim; i++) {
      d_func[MESH_DISPLACEMENT1 + i] *= time_factor;
    }
  }

  *func -= fv->v[var_flag - VELOCITY1];

} /* END of routine fvelocity_parabola                                        */
/*****************************************************************************/
void f_vestress_parabola(const int var_flag,
                         const int ielem_dim,
                         const int velo_condition,
                         const int mn,
                         double *func,
                         double d_func[],     /* defined [MAX_VARIABLE_TYPES + MAX_CONC] */
                         const double p[],    /* parameters passed in from the input deck*/
                         const double time,   /* time at which bc's are evaluated   */
                         const int num_const) /* number of passed parameters   */
{
  /*    parabolic velocity profile
   *      p[0] = coordinate1
   *      p[1] = coordinate2
   *      p[2] = flow in positive coordinate direction
   */
  double coord1, coord2, qflow, gap, pre_factor, tmp, expon;
  double pl_index = 1.0, srate = 0.0, temp = 25.0, at = 1., d_at_dT = 0., wlf_denom;
  double alpha = -1., lambda = -1.;
  double Ws, A_alpha, Ksqr, f, mup;
  double gamma[DIM][DIM];
  int i, mode = 0, strs = 0;

  if (!pd->v[pg->imtrx][POLYMER_STRESS11]) {
    GOMA_EH(-1, "Polymer Stress needed for VE Stress PARABOLA BC.");
  }

  if (var_flag >= POLYMER_STRESS11_7) {
    mode = 7;
    strs = var_flag - POLYMER_STRESS11_7;
  } else if (var_flag >= POLYMER_STRESS11_6) {
    mode = 6;
    strs = var_flag - POLYMER_STRESS11_6;
  } else if (var_flag >= POLYMER_STRESS11_5) {
    mode = 5;
    strs = var_flag - POLYMER_STRESS11_5;
  } else if (var_flag >= POLYMER_STRESS11_4) {
    mode = 4;
    strs = var_flag - POLYMER_STRESS11_4;
  } else if (var_flag >= POLYMER_STRESS11_3) {
    mode = 3;
    strs = var_flag - POLYMER_STRESS11_3;
  } else if (var_flag >= POLYMER_STRESS11_2) {
    mode = 2;
    strs = var_flag - POLYMER_STRESS11_2;
  } else if (var_flag >= POLYMER_STRESS11_1) {
    mode = 1;
    strs = var_flag - POLYMER_STRESS11_1;
  } else if (var_flag >= POLYMER_STRESS11) {
    mode = 0;
    strs = var_flag - POLYMER_STRESS11;
  } else {
    GOMA_EH(-1, "Polymer Stress mode not found - VE Stress PARABOLA BC.");
  }

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }
  /*  shift factor  */
  if (vn->shiftModel == CONSTANT) {
    at = vn->shift[0];
    d_at_dT = 0.;
  } else if (vn->shiftModel == MODIFIED_WLF) {
    wlf_denom = vn->shift[1] + temp - mp->reference[TEMPERATURE];
    if (wlf_denom != 0.) {
      at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - temp) / wlf_denom);
      d_at_dT = -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom);
    }
  }

  coord1 = MIN(p[0], p[1]);
  coord2 = MAX(p[1], p[0]);
  qflow = p[2];
  gap = fabs(coord2 - coord1);
  pre_factor = 6. * qflow / (gap * gap * gap);
  switch (pd->CoordinateSystem) {
  case CARTESIAN:
  case CARTESIAN_2pt5D:
    pre_factor = 6. * qflow / (gap * gap * gap);
    break;
  case CYLINDRICAL:
  case SWIRLING:
    switch (velo_condition) {
    case U_PARABOLA_BC:
      if (coord1 <= DBL_SMALL) {
        pre_factor = 2. * qflow / M_PIE / SQUARE(SQUARE(coord2));
      } else {
        pre_factor = 2. * qflow / M_PIE / (SQUARE(coord2) - SQUARE(coord1)) /
                     (SQUARE(coord2) + SQUARE(coord1) -
                      (SQUARE(coord2) - SQUARE(coord1)) / log(coord2 / coord1));
      }
      break;
    case V_PARABOLA_BC:
      pre_factor = 3. * qflow / M_PIE / (gap * gap * gap);
      break;
    }
    break;
  default:
    GOMA_EH(-1, "Stress parabola not ready for that Coordinate System yet!\n");
  }

  if (ielem_dim > 2) {
    GOMA_EH(-1, "Stress parabola not ready for 3D yet!\n");
    return;
  }
  for (i = 0; i < ielem_dim; i++) {
    d_func[MESH_DISPLACEMENT1 + i] = 0.0;
  }

  if (af->Assemble_LSA_Mass_Matrix)
    d_func[var_flag] = 0.0;
  else
    d_func[var_flag] = -1.0;

  if (gap > DBL_SMALL) {
    if (num_const == 3 || p[3] == 1.0) /*  Newtonian solution   */
    {

      switch (pd->CoordinateSystem) {
      case CARTESIAN:
      case CARTESIAN_2pt5D:
        switch (velo_condition) {
        case U_PARABOLA_BC:
          srate = pre_factor * (coord1 + coord2 - 2. * fv->x[1]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT2] = -2. * pre_factor;
          }
          break;
        case V_PARABOLA_BC:
          srate = pre_factor * (coord1 + coord2 - 2. * fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = -2. * pre_factor;
          }
          break;
        case W_PARABOLA_BC:
          srate = pre_factor * (coord1 + coord2 - 2. * fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = -2. * pre_factor;
          }
          break;
        default:
          *func = 0.;
        }
        break;
      case CYLINDRICAL:
      case SWIRLING:
        switch (velo_condition) {
        case U_PARABOLA_BC:
          if (coord1 <= DBL_SMALL) {
            srate = pre_factor * (-2. * fv->x[1]);
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] = -2. * pre_factor;
            }
          } else {
            srate = pre_factor * (-2. * fv->x[1] + (SQUARE(coord2) - SQUARE(coord1)) /
                                                       log(coord2 / coord1) / fv->x[1]);
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] =
                  pre_factor * (-2. - (SQUARE(coord2) - SQUARE(coord1)) / log(coord2 / coord1) /
                                          SQUARE(fv->x[1]));
            }
          }
          break;
        case V_PARABOLA_BC:
          *func = pre_factor / fv->x[1] * (fv->x[0] - coord1) * (coord2 - fv->x[0]);
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] = pre_factor / fv->x[1] * (coord1 + coord2 - 2. * fv->x[0]);
            d_func[MESH_DISPLACEMENT2] = -(*func) / fv->x[1];
          }
          break;
        default:
          *func = 0.;
        }
        break;
      }
    } else if (num_const > 3) /*  Power-law  solution   */
    {
      if (p[3] < 0.0) {
        pl_index = gn->nexp;
      } else {
        pl_index = p[3];
      }
      switch (pd->CoordinateSystem) {
      case CARTESIAN:
      case CARTESIAN_2pt5D:
        expon = 1. + 1. / pl_index;
        pre_factor = (2. * pl_index + 1.) / (pl_index + 1.) * qflow / pow(gap, expon + 1.);
        switch (velo_condition) {
        case U_PARABOLA_BC:
          tmp = 2 * fv->x[1] - coord1 - coord2;
          srate = pre_factor * (-2. * SGN(tmp) * expon * pow(fabs(tmp), 1. / pl_index));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT2] =
                -4. * pre_factor * (expon * (expon - 1.) * pow(fabs(tmp), 1. / pl_index - 1.));
          }
          break;
        case V_PARABOLA_BC:
          tmp = 2 * fv->x[0] - coord1 - coord2;
          srate = pre_factor * (-2. * SGN(tmp) * expon * pow(fabs(tmp), 1. / pl_index));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                -4. * pre_factor * (expon * (expon - 1.) * pow(fabs(tmp), 1. / pl_index - 1.));
          }
          break;
        case W_PARABOLA_BC:
          tmp = 2 * fv->x[0] - coord1 - coord2;
          srate = pre_factor * (-2. * SGN(tmp) * expon * pow(fabs(tmp), 1. / pl_index));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                -4. * pre_factor * (expon * (expon - 1.) * pow(fabs(tmp), 1. / pl_index - 1.));
          }
          break;
        default:
          srate = 0.;
        }
        break;
      case CYLINDRICAL:
      case SWIRLING:
        expon = 1. + 1. / pl_index;
        switch (velo_condition) {
        case U_PARABOLA_BC:
          if (coord1 <= DBL_SMALL) {
            pre_factor =
                (3. * pl_index + 1.) / (pl_index + 1.) * qflow / M_PIE / pow(gap, expon + 2.);
            srate = pre_factor * (-expon * pow(fv->x[1], expon - 1.));
            if (pd->e[pg->imtrx][R_MESH1]) {
              d_func[MESH_DISPLACEMENT2] =
                  pre_factor * (-expon * (expon - 1.) * pow(fv->x[1], expon - 2.));
            }
          } else {
            GOMA_EH(-1, "Power-law annulus not done yet!\n");
          }
          break;
        case V_PARABOLA_BC:
          pre_factor =
              (2. * pl_index + 1.) / (pl_index + 1.) * qflow / M_PIE / pow(gap, expon + 1.);
          tmp = 2 * fv->x[0] - coord1 - coord2;
          srate = pre_factor / fv->x[1] * (-2. * SGN(tmp) * expon * pow(fabs(tmp), expon - 1.));
          if (pd->e[pg->imtrx][R_MESH1]) {
            d_func[MESH_DISPLACEMENT1] =
                pre_factor / fv->x[1] * (-4. * expon * (expon - 1.) * pow(fabs(tmp), expon - 2.));
            d_func[MESH_DISPLACEMENT2] = -(srate) / fv->x[1];
          }
          break;
        default:
          srate = 0.;
        }
        break;
      }
    }
  } else {
    srate = 0.0;
  }

  /* Compute stresses from shear rate */
  switch (vn->ConstitutiveEquation) {

  case GIESEKUS:
    alpha = ve_glob[mn][mode]->alpha;
    lambda = ve_glob[mn][mode]->time_const;
    Ws = at * lambda * fabs(srate);
    A_alpha = 8 * alpha * (1. - alpha);
    mup = viscosity(ve[mode]->gn, gamma, NULL);
    if (Ws >= 0.1) {
      Ksqr = (sqrt(1. + 2 * A_alpha * SQUARE(Ws)) - 1.) / (A_alpha * SQUARE(Ws));
    } else {
      Ksqr = 1. - 0.5 * A_alpha * SQUARE(Ws) + 0.5 * SQUARE(A_alpha) * pow(Ws, 4) -
             0.625 * pow(A_alpha, 3) * pow(Ws, 6) + 0.875 * pow(A_alpha, 4) * pow(Ws, 8);
    }
    f = (1. - sqrt(Ksqr)) / (1. + (1. - 2 * alpha) * sqrt(Ksqr));

    switch (strs) {
    case 0: /* S11 */
      *func = 2 * mup * lambda * f * (1. - alpha * f) / (SQUARE(at * lambda) * alpha * (1. - f));
      *func += -mup * lambda * f / SQUARE(at * lambda);
      break;
    case 1: /* S12 */
      *func = srate * at * mup * SQUARE(1. - f) / (1. + (1. - 2. * alpha) * f);
      d_func[TEMPERATURE] = d_at_dT * srate * mup;
      break;
    case 2: /* S22 */
      *func = -mup * lambda * f / SQUARE(at * lambda);
      break;
    case 3: /* S13 */
    case 4: /* S23 */
    case 5: /* S33 */
      *func = 0.0;
      break;
    }
    break;
  case OLDROYDB:
    lambda = ve_glob[mn][mode]->time_const;
    Ws = at * lambda * fabs(srate);
    mup = viscosity(ve[mode]->gn, gamma, NULL);
    switch (strs) {
    case 0: /* S11 */
      *func = 2 * mup * lambda * SQUARE(srate);
      break;
    case 1: /* S12 */
      *func = srate * at * mup;
      d_func[TEMPERATURE] = d_at_dT * srate * mup;
      break;
    case 2: /* S22 */
      *func = 0.0;
      break;
    case 3: /* S13 */
    case 4: /* S23 */
    case 5: /* S33 */
      *func = 0.0;
      break;
    }
    break;
  }

  switch (strs) {
  case 0: /* S11 */
    *func -= fv->S[mode][0][0];
    break;
  case 1: /* S12 */
    *func -= fv->S[mode][0][1];
    break;
  case 2: /* S22 */
    *func -= fv->S[mode][1][1];
    break;
  case 5: /* S33 */
    *func -= fv->S[mode][2][2];
    break;
  }

} /* END of routine f_vestress_parabola                                        */
/*****************************************************************************/

void fspline(int ielem_dim,
             double *func,
             double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
             double p[],      /* parameters to parameterize temperature eqn model*/
             double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[MESH_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[MESH_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[MESH_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline                                                  */
/*****************************************************************************/

void fspline_rs(int ielem_dim,
                double *func,
                double d_func[], /* dimensioned [MAX_VARIABLE_TYPES + MAX_CONC] */
                double p[],      /* parameters to parameterize temperature eqn model*/
                double time)     /* time  at which bc's are evaluated     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  d_func[SOLID_DISPLACEMENT1] = dfncd1(fv->x[0], fv->x[1], fv->x[2], p, time);

  d_func[SOLID_DISPLACEMENT2] = dfncd2(fv->x[0], fv->x[1], fv->x[2], p, time);

  if (ielem_dim == 3)
    d_func[SOLID_DISPLACEMENT3] = dfncd3(fv->x[0], fv->x[1], fv->x[2], p, time);

  *func = fnc(fv->x[0], fv->x[1], fv->x[2], p, time);

} /* END of routine fspline_rs                                               */
/*****************************************************************************/

void fTmelting(double *func,
               double d_func[], /* dimensioned MAX_VARIABLE_TYPES + MAX_CONC */
               double a1)       /*  function parameter from data card     */
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;
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
  int index_var;       /* Column index into the global stiffness matrix*/
  dbl x_var;           /* value of variable at this node */
  dbl d_x_var;         /* sensitivity of variable to nodal unknown */
  dbl d_vect_var[DIM]; /* sensitivity of vector variable to nodal unknown */
  dbl slope;           /* slope of interpolated function in table */
  dbl x_var_mp[1];     /* dummy variable for table lookup subroutines */

  if (af->Assemble_LSA_Mass_Matrix)
    return 0;

  /* ---- Find variable number and species number */

  jvar = BC_Types[bc_input_id].BC_Data_Int[2];

  wspec = BC_Types[bc_input_id].BC_Data_Int[3];

  /* put value of variable in GD Condition into x_var and put it's sensitivity in d_x_var */
  index_var = load_variable(&x_var, &d_x_var, jvar, wspec, tt, dt, d_vect_var);

  if (jvar == SPEED) {
    vector_sens = 1;
  } else {
    vector_sens = 0;
  }
  /* Now add in contributions to residual vector and jacobian matrix */

  switch (gd_condition) {
  case (GD_CONST_BC): /* x - c0 */

    *func = (x_var - BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = d_vect_var[b];
        }
      } else {
        d_func[index_var] = d_x_var;
      }
    }
    break;

  case (GD_LINEAR_BC): /* C1 x + c0 */

    *func =
        (x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1];
        }
      } else {
        d_func[index_var] = d_x_var * BC_Types[bc_input_id].BC_Data_Float[1];
      }
    }
    break;

  case (GD_INVERSE_BC): /* C1/x + c0 */

    *func =
        (BC_Types[bc_input_id].BC_Data_Float[1] / x_var + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              -d_vect_var[b] * BC_Types[bc_input_id].BC_Data_Float[1] / (x_var * x_var);
        }
      } else {
        d_func[index_var] = -d_x_var * BC_Types[bc_input_id].BC_Data_Float[1] / (x_var * x_var);
      }
    }
    break;

  case (GD_PARAB_BC): /* C2 x^2 + C1 x + c0 */

    *func =
        (x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
         x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (BC_Types[bc_input_id].BC_Data_Float[1] +
                               2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]);
        }
      } else {
        d_func[index_var] = d_x_var * (BC_Types[bc_input_id].BC_Data_Float[1] +
                                       2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2]);
      }
    }
    break;
  case (GD_PARAB_OFFSET_BC): /* C2 (x - C3)^2 + C1 (x - C3) + c0 */

    *func =
        ((x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
             (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
             BC_Types[bc_input_id].BC_Data_Float[2] +
         (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) * BC_Types[bc_input_id].BC_Data_Float[1] +
         BC_Types[bc_input_id].BC_Data_Float[0]);
    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (BC_Types[bc_input_id].BC_Data_Float[1] +
                               2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
                                   BC_Types[bc_input_id].BC_Data_Float[2]);
        }
      } else {
        d_func[index_var] = d_x_var * (BC_Types[bc_input_id].BC_Data_Float[1] +
                                       2. * (x_var - BC_Types[bc_input_id].BC_Data_Float[3]) *
                                           BC_Types[bc_input_id].BC_Data_Float[2]);
      }
    }
    break;
  case (GD_CIRC_BC): /* C2 ( x - C1 )^2  - c0^2 */
    /* C2 represents ellipticity and C1 represents origin */
    /* C0 is radius and should enter only one of the BC's */

    *func =
        (BC_Types[bc_input_id].BC_Data_Float[2] * (x_var - BC_Types[bc_input_id].BC_Data_Float[1]) *
             (x_var - BC_Types[bc_input_id].BC_Data_Float[1]) -
         BC_Types[bc_input_id].BC_Data_Float[0] * BC_Types[bc_input_id].BC_Data_Float[0]);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] * (2. * BC_Types[bc_input_id].BC_Data_Float[2] *
                               (x_var - BC_Types[bc_input_id].BC_Data_Float[1]));
        }
      } else {
        d_func[index_var] = d_x_var * (2. * BC_Types[bc_input_id].BC_Data_Float[2] *
                                       (x_var - BC_Types[bc_input_id].BC_Data_Float[1]));
      }
    }
    break;

  case (GD_POLYN_BC): /* up to 6th order polynomial */
                      /* C6 x^6 + C5 x^5 + C4 x^4 + C3 x^3 + C2 x^2 + C1 x + C0 */

    *func =
        (x_var * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6] +
         x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
         x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
         x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
         x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
         x_var * BC_Types[bc_input_id].BC_Data_Float[1] + BC_Types[bc_input_id].BC_Data_Float[0]);

    /* printf("POLYN fit X,F = %f %f\n", x_var, *func); */

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] =
              d_vect_var[b] *
              (BC_Types[bc_input_id].BC_Data_Float[1] +
               2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
               3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
               4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
               5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
               6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6]);
        }
      } else {
        d_func[index_var] =
            d_x_var *
            (BC_Types[bc_input_id].BC_Data_Float[1] +
             2. * x_var * BC_Types[bc_input_id].BC_Data_Float[2] +
             3. * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[3] +
             4. * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[4] +
             5. * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[5] +
             6. * x_var * x_var * x_var * x_var * x_var * BC_Types[bc_input_id].BC_Data_Float[6]);
      }
    }
    break;

  case (GD_TABLE_BC):

    *func = BC_Types[bc_input_id].BC_Data_Float[0];

    x_var_mp[0] = x_var;
    *func *= interpolate_table(BC_Types[bc_input_id].table, x_var_mp, &slope, NULL);

    if (af->Assemble_Jacobian) {
      if (vector_sens) {
        for (b = 0; b < DIM; b++) {
          d_func[index_var + b] = BC_Types[bc_input_id].BC_Data_Float[0] * slope * d_vect_var[b];
        }
      } else {
        d_func[index_var] = BC_Types[bc_input_id].BC_Data_Float[0] * slope * d_x_var;
      }
    }

    break;

  default:
    return (-1);
  }

  return (0);
} /* END of routine fgeneralized_dirichlet                                   */
/*****************************************************************************/

void fmesh_constraint(double *func, double d_func[], const int bc_input_id) {
  GOMA_EH(GOMA_ERROR, "CGM not supported, MESH_CONSTRAINT_BC");
  /*#endif  */

} /* END of routine fmesh_constraint                                   */

/*****************************************************************************/
/*****************************************************************************/

int load_variable(double *x_var,       /* variable value */
                  double *d_x_var,     /* sensitivities of variable value */
                  int jvar,            /* variable number */
                  int wspec,           /* species number */
                  double tt,           /* parameter to vary time integration from
                                          explicit (tt = 1) to implicit (tt = 0) */
                  double dt,           /* current time step size */
                  double d_vect_var[]) /* vector sensitivities  */

/******************************************************************************

   Function which calculates the value of a chosen variable at the current point
    (using the field variables); available variable names listed in mm_names.h

   Author:          Rich Cairncross (1511)
   Date:            22 MAY 1995
   Revised:
******************************************************************************/

{
  int var = -1, b;
  *x_var = 0.;
  *d_x_var = 0.;

  memset(d_vect_var, 0, DIM * sizeof(double));

  if (jvar >= D_VEL1_DT && jvar <= D_P_DT) {
    if (pd->TimeIntegration == STEADY)
      GOMA_EH(GOMA_ERROR, "Unsteady GD for Steady problem");
  }

  /* ---- Find variable value, sensitivities, and species number */

  switch (jvar) {
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
  case USTAR:
    *x_var = fv->v_star[0];
    var = USTAR;
    *d_x_var = 1.;
    break;
  case VSTAR:
    *x_var = fv->v_star[1];
    var = VSTAR;
    *d_x_var = 1.;
    break;
  case WSTAR:
    *x_var = fv->v_star[2];
    var = WSTAR;
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
  case EM_CONT_REAL:
    *x_var = fv->epr;
    var = EM_CONT_REAL;
    *d_x_var = 1.;
    break;
  case EM_CONT_IMAG:
    *x_var = fv->epi;
    var = EM_CONT_IMAG;
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
  case SHELL_SAT_1:
    *x_var = fv->sh_sat_1;
    var = SHELL_SAT_1;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_2:
    *x_var = fv->sh_sat_2;
    var = SHELL_SAT_2;
    *d_x_var = 1.;
    break;
  case SHELL_SAT_3:
    *x_var = fv->sh_sat_3;
    var = SHELL_SAT_3;
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
  case EM_E1_REAL:
    *x_var = fv->em_er[0];
    var = EM_E1_REAL;
    *d_x_var = 1.;
    break;
  case EM_E2_REAL:
    *x_var = fv->em_er[1];
    var = EM_E2_REAL;
    *d_x_var = 1.;
    break;
  case EM_E3_REAL:
    *x_var = fv->em_er[2];
    var = EM_E3_REAL;
    *d_x_var = 1.;
    break;
  case EM_E1_IMAG:
    *x_var = fv->em_ei[0];
    var = EM_E1_IMAG;
    *d_x_var = 1.;
    break;
  case EM_E2_IMAG:
    *x_var = fv->em_ei[1];
    var = EM_E2_IMAG;
    *d_x_var = 1.;
    break;
  case EM_E3_IMAG:
    *x_var = fv->em_ei[2];
    var = EM_E3_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H1_REAL:
    *x_var = fv->em_hr[0];
    var = EM_H1_REAL;
    *d_x_var = 1.;
    break;
  case EM_H2_REAL:
    *x_var = fv->em_hr[1];
    var = EM_H2_REAL;
    *d_x_var = 1.;
    break;
  case EM_H3_REAL:
    *x_var = fv->em_hr[2];
    var = EM_H3_REAL;
    *d_x_var = 1.;
    break;
  case EM_H1_IMAG:
    *x_var = fv->em_hi[0];
    var = EM_H1_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H2_IMAG:
    *x_var = fv->em_hi[1];
    var = EM_H2_IMAG;
    *d_x_var = 1.;
    break;
  case EM_H3_IMAG:
    *x_var = fv->em_hi[2];
    var = EM_H3_IMAG;
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
    GOMA_EH(GOMA_ERROR, "SURFACE variables not defined yet");
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
    for (b = 0; b < pd->Num_Dim; b++) {
      *x_var += SQUARE(fv->v[b]);
    }
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D) {
      *x_var += SQUARE(fv->v[pd->Num_Dim]);
    }
    *x_var = sqrt(*x_var);
    var = VELOCITY1;
    *d_x_var = 1. / (*x_var);
    for (b = 0; b < pd->Num_Dim; b++) {
      d_vect_var[b] += fv->v[b] * (*d_x_var);
    }
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D) {
      d_vect_var[pd->Num_Dim] += fv->v[pd->Num_Dim] * (*d_x_var);
    }
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
    GOMA_EH(GOMA_ERROR, "SURFACE variables not defined yet");
    break;
  case D_P_DT:
    *x_var = fv_dot->P;
    var = PRESSURE;
    *d_x_var = (1. + 2. * tt) / dt;
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Illegal option in load_variable");
  } /* end of switch on jvar */

  if (wspec != 0 && var != MASS_FRACTION && var != POR_GAS_PRES) {
    GOMA_EH(GOMA_ERROR, "Non-zero species number for wrong variable");
  }
  if (var == MASS_FRACTION)
    var = MAX_VARIABLE_TYPES + wspec;

  return (var);
} /* END of routine load_variable()                                          */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int bc_eqn_index(int id,          /* local node number                 */
                 int I,           /* processor node number             */
                 int bc_input_id, /* boundary condition number         */
                 int curr_matID,  /* Current material ID */
                 int kdir,        /* coordinate index for mesh or velocity
                                   * - normally this is zero unless a
                                   *   bc is a vector condition                 */
                 int *eqn,        /* eqn to which this condition is applied     */
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
  int ieqn, i_calc, jeqn, bc_node, offset_j;
  int matID, index_eqn, num, node_offset;
  int i, index;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd, *vd2 = NULL;
  nv = node->Nodal_Vars_Info[pg->imtrx];

  /*
   *  Find equation number and species number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  matID = curr_matID;

  /* need to change equation index for mesh or velocity */
  if (kdir != 0) {
    if (ieqn == R_MESH1)
      ieqn += kdir;
    else if (ieqn == R_MOMENTUM1)
      ieqn += kdir;
    else if (ieqn == R_SOLID1)
      ieqn += kdir;
    else if (ieqn == R_LAGR_MULT1)
      ieqn += kdir;
    else if (ieqn == R_EM_H1_REAL)
      ieqn += kdir; // AMC: These vector bcs are not rotated..
    else if (ieqn == R_EM_H1_IMAG)
      ieqn += kdir;
    else if (ieqn == R_EM_E1_REAL)
      ieqn += kdir;
    else if (ieqn == R_EM_E1_IMAG)
      ieqn += kdir;
  }

  /*
   *  Find out the offset under non-conflicting conditions
   *  If a boundary condition isn't to be applied under non-
   *  conflicting conditions, it won't be applied under
   *  conflicting conditions. So exit here if bc is not
   *  to be applied.
   */
  node_offset = find_bc_unk_offset(bc, curr_matID, I, kdir, &matID, &vd);
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
  if (BC_dup_ptr != NULL) {
    bc_node = in_list(I, 0, BC_dup_ptr[pg->imtrx] + 1, BC_dup_nodes[pg->imtrx]);
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
        i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][node_offset]);

      } else if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
        if (goma_automatic_rotations.automatic_rotations) {
          if (ieqn == R_MESH_NORMAL) {
            int best_dir = -1;
            double dot_max = 0;
            for (int dir = 0; dir < 3; dir++) {
              double dot = 0;
              for (int j = 0; j < 3; j++) {
                dot +=
                    goma_automatic_rotations.rotation_nodes[I].rotated_coord[dir]->normal->data[j] *
                    fv->snormal[j];
              }
              if (fabs(dot) > dot_max) {
                best_dir = dir;
                dot_max = fabs(dot);
              }
            }
            offset_j = get_nodal_unknown_offset(nv, R_MESH1 + best_dir, matID, 0, &vd2);
            if (offset_j >= 0) {
              node_offset = offset_j;
              ieqn = R_MESH1 + best_dir;
              vd = vd2;
              i_calc = 0;
            }
          } else {
            GOMA_EH(GOMA_ERROR, "Automatic rotations only for normal conditions");
          }
        } else {
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
              i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
              if (i_calc != -1) {
                node_offset = offset_j;
                ieqn = jeqn;
                vd = vd2;
                break;
              }
            }
          }
        }
      } else if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
        if (goma_automatic_rotations.automatic_rotations) {
          if (ieqn == R_MOM_NORMAL) {
            int best_dir = -1;
            double dot_max = 0;
            for (int dir = 0; dir < 3; dir++) {
              double dot = 0;
              for (int j = 0; j < 3; j++) {
                dot +=
                    goma_automatic_rotations.rotation_nodes[I].rotated_coord[dir]->normal->data[j] *
                    fv->snormal[j];
              }
              if (fabs(dot) > dot_max) {
                best_dir = dir;
                dot_max = fabs(dot);
              }
            }
            offset_j = get_nodal_unknown_offset(nv, R_MOMENTUM1 + best_dir, matID, 0, &vd2);
            if (offset_j >= 0) {
              node_offset = offset_j;
              ieqn = R_MOMENTUM1 + best_dir;
              vd = vd2;
              i_calc = 0;
            }
          } else {
            GOMA_EH(GOMA_ERROR, "Automatic rotations only for normal conditions");
          }
        } else {
          for (jeqn = R_MOMENTUM1; jeqn <= R_MOMENTUM3; jeqn++) {
            offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
            if (offset_j >= 0) {
              i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
              if (i_calc != -1) {
                node_offset = offset_j;
                ieqn = jeqn;
                vd = vd2;
                break;
              }
            }
          }
        }
      } else if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
        for (jeqn = R_SOLID1; jeqn <= R_SOLID3; jeqn++) {
          offset_j = get_nodal_unknown_offset(nv, jeqn, matID, 0, &vd2);
          if (offset_j >= 0) {
            i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][offset_j]);
            if (i_calc != -1) {
              node_offset = offset_j;
              ieqn = jeqn;
              vd = vd2;
              break;
            }
          }
        }
      } else {
        GOMA_EH(GOMA_ERROR, "Equation not in list ");
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
  if (node->DBC[pg->imtrx]) {
    if (((int)node->DBC[pg->imtrx][node_offset]) != -1) {
      return -1;
    }
  }

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown[pg->imtrx] + node_offset;

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
      if (matID == vd->MatID || (!(*vd_retn) && vd->MatID == -1)) {
        *vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index                                             */

int bc_eqn_index_stress(int id,          /* local node number                 */
                        int I,           /* processor node number             */
                        int bc_input_id, /* boundary condition number         */
                        int curr_matID,  /* Current material ID */
                        int kdir,        /* coordinate index for stress components */
                        int mode,        /* Stress mode number */
                        int *eqn,        /* eqn to which this condition is applied     */
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
  int ieqn, i_calc, bc_node;
  int matID, index_eqn, num, node_offset;
  int i, index;
  int ndofs_ieqn;
  NODAL_VARS_STRUCT *nv;
  struct Boundary_Condition *bc = BC_Types + bc_input_id;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node = Nodes[I];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  nv = node->Nodal_Vars_Info[pg->imtrx];

  /*
   *  Find equation number from the BC description
   *  structure
   */
  ieqn = bc_desc->equation;
  /* Sanity check */
  if (ieqn != R_STRESS11)
    GOMA_EH(GOMA_ERROR, "You can't be here");

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
    GOMA_EH(GOMA_ERROR, "Maximum 8 modes allowed here!");
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
  bc_node = in_list(I, 0, BC_dup_ptr[pg->imtrx] + 1, BC_dup_nodes[pg->imtrx]);
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
      i_calc = search_bc_dup_list(bc_input_id, BC_dup_list[pg->imtrx][bc_node][node_offset]);
    } else {
      GOMA_EH(GOMA_ERROR, "Equation not in list ");
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
  if (node->DBC[pg->imtrx]) {
    if (((int)node->DBC[pg->imtrx][node_offset]) != -1) {
      return -1;
    }
  }

  /*
   * Set up the return variables
   */
  *matID_retn = matID;
  *eqn = ieqn;
  index_eqn = node->First_Unknown[pg->imtrx] + node_offset;

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
      if (matID == vd->MatID || (!(*vd_retn) && vd->MatID == -1)) {
        *vd_retn = vd;
      }
    }
  }

  return index_eqn;
} /* END of routine bc_eqn_index_stress                                      */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int evaluate_time_func(const double current_time,
                       double *f_time, /* computed time function */
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

  /* Check if max time was specified  and reset time if so */
  if (BC_Types[bc_input_id].BC_Data_Int[3] == GD_TIME_MAX) {
    if (time > BC_Types[bc_input_id].BC_Data_Float[2]) {
      time = BC_Types[bc_input_id].BC_Data_Float[2];
    }
  }

  /* ---- Find variable number and species number */
  time_function = BC_Types[bc_input_id].BC_Data_Int[2];

  /* Now compute time function */

  switch (time_function) {
  case (GD_TIME_LIN): /* c0 +c1*t */

    *f_time =
        BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time;

    break;

  case (GD_TIME_EXP): /* exp(c0 +c1*t) */

    *f_time =
        exp(BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time);

    break;

  case (GD_TIME_SIN): /* sin(c0 +c1*t) */

    *f_time =
        sin(BC_Types[bc_input_id].BC_Data_Float[0] + BC_Types[bc_input_id].BC_Data_Float[1] * time);

    break;
  case (GD_TIME_TABLE): {
    double slope, _time[1];
    _time[0] = time;
    *f_time = interpolate_table(BC_Types[bc_input_id].table, _time, &slope, NULL);
  } break;
  default:
    return (-1);
  }

  return (0);
} /* END of routine evaluate_time_function                                   */
/*****************************************************************************/

void apply_table_bc(double *func,
                    double d_func[MAX_VARIABLE_TYPES + MAX_CONC],
                    struct Boundary_Condition *BC_Type,
                    double time_value) {
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
  double slope, interp_val, x_table[2];
  double dfunc_dx[3];

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  basis = BC_Type->table->t_index[0];

  /*  Setup Dummy variable to pass array to interpolation function */
  if (basis != -1)
    x_table[0] = fv->x[basis];
  else
    x_table[0] = time_value;

  if (BC_Type->table->interp_method == BIQUADRATIC || BC_Type->table->interp_method == BILINEAR) {
    x_table[1] = fv->x[BC_Type->table->t_index[1]];
  }

  interp_val = interpolate_table(BC_Type->table, x_table, &slope, dfunc_dx);
  interp_val *= BC_Type->table->yscale;
  slope *= BC_Type->table->yscale;

  var = BC_Type->table->f_index;

  switch (var) {
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
    GOMA_EH(GOMA_ERROR, "Variable not yet implemented in TABLE_BC");
    break;
  }

  /* And, at the last, we account for the dependence of the interpolation
   * on the basis coordinate
   */

  if (BC_Type->table->interp_method == BIQUADRATIC || BC_Type->table->interp_method == BILINEAR) {
    d_func[R_MESH1 + BC_Type->table->t_index[0]] -= dfunc_dx[0] * BC_Type->BC_Data_Float[0];
    d_func[R_MESH1 + BC_Type->table->t_index[1]] -= dfunc_dx[1] * BC_Type->BC_Data_Float[0];
  } else {
    if (basis != -1 && pd->e[pg->imtrx][R_MESH1 + basis]) {
      d_func[R_MESH1 + basis] -= slope;
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
