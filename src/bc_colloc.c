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
 */

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ad_turbulence.h"
#include "bc/eqn_index.h"
#include "bc/generalized_dirichlet.h"
#include "bc/geom.h"
#include "bc/rotate_coordinates.h"
#include "bc/user_geom.h"
#include "bc_colloc.h"
#include "bc_integ.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "gds/gds_vector.h"
#include "load_field_variables.h"
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
#include "table.h"
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
  const double penalty = upd->strong_penalty;
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

      if (upd->AutoDiff) {
#ifdef GOMA_ENABLE_SACADO
        fill_ad_field_variables();
#else
        GOMA_EH(GOMA_ERROR, "AutoDiff assembly enabled but Goma not compiled with Sacado support");
#endif
      }

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

            case SA_WALL_FUNC_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              ad_sa_wall_func(kfunc, d_kfunc);
              func = kfunc[0];
              d_func[EDDY_NU] = 1.0;
              break;
            case OMEGA_WALL_FUNC_BC:
              memset(kfunc, 0, DIM * sizeof(double));
              ad_omega_wall_func(kfunc, d_kfunc);
              func = kfunc[0];
              d_func[TURB_OMEGA] = 1.0;
              break;

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

            case DOUBLE_FILLET_BC:
              f_double_fillet(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
                              BC_Types[bc_input_id].len_u_BC);
              break;

            case DOUBLE_FILLET_GEOM_BASED_BC:
              f_double_fillet_geom_based(ielem_dim, &func, d_func, BC_Types[bc_input_id].u_BC,
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
                                               offset, &(BC_Types[bc_input_id]));
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
                          // I'm not sure this phi_j is ever going to be non-zero with the Dolphin
                          // check above check if we have a bf available before trying to use it
                          if (bf[var]) {
                            phi_j = bf[var]->phi[j];
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                                penalty * d_func[var] * phi_j;
                            lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] *= f_time;
                          }
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
  return (status);
} /* end of routine apply_point_colloc_bc() */

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
    lambda = ve_glob[mn][mode]->time_const_st->lambda0;
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
    lambda = ve_glob[mn][mode]->time_const_st->lambda0;
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