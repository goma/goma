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
 * Routines for calculating Boundary Conditions and Adding them to matrix_fill
 */

/*
 * $Id: bc_special.c,v 5.19 2010-07-21 16:39:26 hkmoffa Exp $
 */

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <string.h>

/* GOMA include files */
#define GOMA_BC_SPECIAL_C
#include "ac_stability_util.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_special.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_shell.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_shell_util.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "sl_util_structs.h"
#include "std.h"
#include "user_bc.h"

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */

/*
 * Function definitions.
 */
/*ARGSUSED*/
int apply_special_bc(struct GomaLinearSolverData *ams,
                     double x[],            /* Solution vector for the current processor      */
                     double resid_vector[], /* Residual vector for the current processor   */
                     double x_old[],        /* Solution vector at last time step (previous x) */
                     double x_older[],      /* Solution vector at 2nd to last time*/
                     double xdot[],         /* xdot of current solution                       */
                     double xdot_old[],     /* xdot at old solution                       */
                     double delta_t,        /* current time step size                         */
                     double theta,          /* parameter to vary time integration:
                                               explicit (theta = 1) -- implicit (theta = 0) */
                     struct elem_side_bc_struct *first_elem_side_BC_array[],
                     /*  An array of pointers to the first surface integral
                      *  defined for each element. It has a length equal to
                      *  the total number of elements defined on the current
                      *  processor
                      */
                     int ielem,      /* element number */
                     int ip_total,   /* total number of gauss points */
                     int ielem_type, /* element type */
                     int num_local_nodes,
                     int ielem_dim,
                     int iconnect_ptr,
                     struct elem_side_bc_struct *elem_side_bc, /* Pointer to an element side
                                                                  boundary condition structure   */
                     int num_total_nodes,
                     int bc_application,
                     int CA_id[],     /*  CA condition id array  */
                     int CA_fselem[], /* array of free surface elements for CA
                                         conditions, initialized to -1 */
                     int CA_sselem[], /* array of solid surface elements for CA
                                         conditions, initialized to -1 */
                     Exo_DB *exo,
                     double time_value)

/******************************************************************************
  Function which is used to coordinate the application of special conditions,
  like endpoint boundary conditions, etc.

  Author:          Randy Schunk et al.
  Date:            12 December 1994
******************************************************************************/

{
  double *a = ams->val;
  int *ija = ams->bindx;
  int w, i, I, ibc, k, j, id, icount, ss_index, i1, i2, i3; /* counters */
  /* HKM - worried that jflag shouldn't be initialized all the way up here */
  int jcnt, jflag = -1, local_node_id = -1, matID_apply;
  int GD_count;
  int Gibbs = 1;
  int iapply = 0;
  int ieqn, eqn, var, pvar, ldof_eqn, p, q, index_eq;
  int err; /* status variable for functions */
  int status = 0;
  int bc_input_id; /* bc # for side bc on this side */
  int j_bc_id;     /* bc # for point bc applied at current node */
  VARIABLE_DESCRIPTION_STRUCT *vd;

  double xi[DIM]; /* Local element coordinates of Gauss point. */
  double x_dot[MAX_PDIM];

  /* Normals for variable normal contact angle condition etc.       */
  static double fsnormal[MAX_CA][MAX_PDIM]; /* Free surface normal component
                                               vector   */
  static double dfsnormal_dx[MAX_CA][MAX_PDIM][MAX_PDIM][MDE];
  /* Free surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j */

  static double ssnormal[MAX_CA][MAX_PDIM]; /* Solid surface normal component
                                               vector  */
  static double dssnormal_dx[MAX_CA][MAX_PDIM][MAX_PDIM][MDE];
  /* Solid surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j */
  static double fsnrml[MAX_PDIM];
  static double dfsnrml_dx[MAX_PDIM][MAX_PDIM][MDE];
  static double ssnrml[MAX_PDIM];
  static double dssnrml_dx[MAX_PDIM][MAX_PDIM][MDE];

  int variable_wall_normal; /* flag to show whether to read a constant
                               solid surface normal from the input deck
                               (0) or calculate it as we go (1)           */
  double sum_normal;
  NODE_INFO_STRUCT *node;
  /***************************************************************************/
  int je, ja;
  double weight = 1e12;
  double rcoord = 0.;
  double dsigma_dx[DIM][MDE];
  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double d_func_ss[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double *f; /* shorthand pointer to BC_Data_Float */

  /* Before we start, initialize the "GD_count" flag for this side to zero */
  /* Recall that we enter this routine for each side on the current element.  */
  /*  This is really put in to handle cases in which more than one GD condition */
  /*  is used to describe geometry.  side_info is getting falsly incremented and */
  /*  fouls other CA apps*/

  GD_count = 0;

  for (ibc = 0; (bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1; ibc++) {
    // pass through for BC's that don't have side sets
    if (BC_Types[bc_input_id].BC_Name == SH_GAMMA1_DERIV_SYMM_BC ||
        BC_Types[bc_input_id].BC_Name == SH_GAMMA2_DERIV_SYMM_BC ||
        BC_Types[bc_input_id].BC_Name == SHELL_TFMP_GRAD_S_BC) {
      continue;
    }
    if ((ss_index = in_list(BC_Types[bc_input_id].BC_ID, 0, exo->num_side_sets, ss_to_blks[0])) ==
        -1) {
      GOMA_EH(GOMA_ERROR, "Cannot match side set id with that in ss_to_blks array");
    }

    /* Set flag to indicate if we're in the right material (only one) to apply*/
    iapply = 0;
    if (exo->eb_id[find_elemblock_index(ielem, exo)] == ss_to_blks[1][ss_index]) {
      iapply = 1;
    }

    /*
     *  However, override if the side set is an external one. In
     *  other words, if the side set is an external side set, we will
     *  apply the boundary condition no matter what.
     */
    if (SS_Internal_Boundary != NULL) {
      if (SS_Internal_Boundary[ss_index] == -1) {
        iapply = 1;
      }
    }
    /* We really need to put a clean indicator whether the bc type is
       for a user-prescribed geometry.  This would shorten these.
       Anyway, if the current condition is one of geometry, increment
       GD_count to indicate that we've got the geometry here */

    if (BC_Types[bc_input_id].BC_Name == PLANE_BC || BC_Types[bc_input_id].BC_Name == SPLINE_BC ||
        BC_Types[bc_input_id].BC_Name == FILLET_BC ||
        BC_Types[bc_input_id].BC_Name == DOUBLE_RAD_BC ||
        BC_Types[bc_input_id].BC_Name == FEATURE_ROLLON_BC ||
        BC_Types[bc_input_id].BC_Name == ROLL_FLUID_BC ||
        BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_BC ||
        BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_RS_BC ||
        BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_PETROV_BC ||
        BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_COLLOC_BC ||
        BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC ||
        BC_Types[bc_input_id].BC_Name == SOLID_FLUID_BC ||
        BC_Types[bc_input_id].BC_Name == SOLID_FLUID_RS_BC ||
        (BC_Types[bc_input_id].BC_Name == GD_CONST_BC ||
         BC_Types[bc_input_id].BC_Name == GD_LINEAR_BC ||
         BC_Types[bc_input_id].BC_Name == GD_INVERSE_BC ||
         BC_Types[bc_input_id].BC_Name == GD_PARAB_BC ||
         BC_Types[bc_input_id].BC_Name == GD_PARAB_OFFSET_BC ||
         BC_Types[bc_input_id].BC_Name == GD_CIRC_BC ||
         BC_Types[bc_input_id].BC_Name == GD_TIME_BC ||
         BC_Types[bc_input_id].BC_Name == GD_POLYN_BC ||
         BC_Types[bc_input_id].BC_Name == GD_TABLE_BC ||
         BC_Types[bc_input_id].BC_Name == MESH_CONSTRAINT_BC))
      GD_count++;

    if (((BC_Types[bc_input_id].BC_Name == KINEMATIC_BC ||
          BC_Types[bc_input_id].BC_Name == KINEMATIC_COLLOC_BC ||
          BC_Types[bc_input_id].BC_Name == KINEMATIC_PETROV_BC ||
          BC_Types[bc_input_id].BC_Name == KINEMATIC_DISC_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_LEAK_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_ELECTRODEPOSITION_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_RS_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_PETROV_BC ||
          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_COLLOC_BC ||
          BC_Types[bc_input_id].BC_Name == PLANE_BC || BC_Types[bc_input_id].BC_Name == SPLINE_BC ||
          BC_Types[bc_input_id].BC_Name == FILLET_BC ||
          BC_Types[bc_input_id].BC_Name == DOUBLE_RAD_BC ||
          BC_Types[bc_input_id].BC_Name == FEATURE_ROLLON_BC ||
          BC_Types[bc_input_id].BC_Name == ROLL_FLUID_BC ||
          BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC ||
          BC_Types[bc_input_id].BC_Name == SOLID_FLUID_BC ||
          BC_Types[bc_input_id].BC_Name == SOLID_FLUID_RS_BC ||
          (BC_Types[bc_input_id].BC_Name == GD_CONST_BC ||
           BC_Types[bc_input_id].BC_Name == GD_LINEAR_BC ||
           BC_Types[bc_input_id].BC_Name == GD_INVERSE_BC ||
           BC_Types[bc_input_id].BC_Name == GD_PARAB_BC ||
           BC_Types[bc_input_id].BC_Name == GD_PARAB_OFFSET_BC ||
           BC_Types[bc_input_id].BC_Name == GD_CIRC_BC ||
           BC_Types[bc_input_id].BC_Name == GD_TIME_BC ||
           BC_Types[bc_input_id].BC_Name == GD_POLYN_BC ||
           BC_Types[bc_input_id].BC_Name == GD_TABLE_BC ||
           BC_Types[bc_input_id].BC_Name == MESH_CONSTRAINT_BC)) &&
         (bc_application == SPECIAL)) ||

        ((BC_Types[bc_input_id].BC_Name == CAPILLARY_BC ||
          BC_Types[bc_input_id].BC_Name == ELEC_TRACTION_BC ||
          BC_Types[bc_input_id].BC_Name == CAP_REPULSE_BC ||
          BC_Types[bc_input_id].BC_Name == CAP_REPULSE_ROLL_BC ||
          BC_Types[bc_input_id].BC_Name == CAP_REPULSE_USER_BC ||
          BC_Types[bc_input_id].BC_Name == CAP_REPULSE_TABLE_BC ||
          BC_Types[bc_input_id].BC_Name == CAPILLARY_TABLE_BC ||
          BC_Types[bc_input_id].BC_Name == CAP_RECOIL_PRESS_BC) &&
         (bc_application == SPECIAL))) {
      /* Initialize jflag on a per-bc basis (untested)  */
      jflag = -1;

      /******************************************************************************/
      /*                LOOP OVER THE LOCAL ELEMENT NODE NUMBER                     */
      /*               FOR NODES THAT LIE ON THE CURRENT SURFACE		        */
      /* this is a loop designed to loop over equations that lie on current surface */
      /******************************************************************************/

      for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {

        /* Find the local element node number for the current node */
        id = (int)elem_side_bc->local_elem_node_id[i];

        /*  Find the processor node number given the local element
         *  node number, 'i' . Then store a pointer to the NODE_INFO
         *  structure pertaining to that node.
         */
        I = Proc_Elem_Connect[iconnect_ptr + id];
        node = Nodes[I];

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

        /* may want to load material properties here */

        /* calculate the determinant of the surface jacobian */
        surface_determinant_and_normal(
            ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1, (int)elem_side_bc->id_side,
            (int)elem_side_bc->num_nodes_on_side, (elem_side_bc->local_elem_node_id));

        if (ielem_dim != 3) {
          calc_surf_tangent(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                            (int)elem_side_bc->num_nodes_on_side,
                            (elem_side_bc->local_elem_node_id));
        }

        if (mp->SurfaceTensionModel != CONSTANT) {
          load_surface_tension(dsigma_dx);
          if (neg_elem_volume)
            return (status);
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

        /* initialize the general function to zero may have more than one entry
         * for vector conditions like capillary */

        memset(func, 0, MAX_PDIM * sizeof(double));
        memset(d_func, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
        memset(d_func_ss, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

        j_bc_id = -1;

        /* Here we will search for nodes which have contact angle conditions.
         * The idea is to check the flag Variable_Mask.DBCA, and then search
         * for the BC card to which this node corresponds.
         */

        /*
         *
         *  START OF Contact Angle CONDITION
         *
         *
         */
        if (node->DBCA == 1 &&
            (BC_Types[bc_input_id].BC_Name == KINEMATIC_BC ||
             BC_Types[bc_input_id].BC_Name == KINEMATIC_COLLOC_BC ||
             BC_Types[bc_input_id].BC_Name == KINEMATIC_PETROV_BC ||
             BC_Types[bc_input_id].BC_Name == KINEMATIC_DISC_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_LEAK_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_ELECTRODEPOSITION_BC || /*  RSL 5/28/02  */
             BC_Types[bc_input_id].BC_Name == PLANE_BC ||
             BC_Types[bc_input_id].BC_Name == SPLINE_BC ||
             BC_Types[bc_input_id].BC_Name == FILLET_BC ||
             BC_Types[bc_input_id].BC_Name == DOUBLE_RAD_BC ||
             BC_Types[bc_input_id].BC_Name == FEATURE_ROLLON_BC ||
             BC_Types[bc_input_id].BC_Name == ROLL_FLUID_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_RS_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_PETROV_BC ||
             BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_COLLOC_BC ||
             BC_Types[bc_input_id].BC_Name == MESH_CONSTRAINT_BC ||
             BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC ||
             BC_Types[bc_input_id].BC_Name == SOLID_FLUID_BC ||
             BC_Types[bc_input_id].BC_Name == SOLID_FLUID_RS_BC ||
             ((BC_Types[bc_input_id].BC_Name == GD_CONST_BC ||
               BC_Types[bc_input_id].BC_Name == GD_LINEAR_BC ||
               BC_Types[bc_input_id].BC_Name == GD_INVERSE_BC ||
               BC_Types[bc_input_id].BC_Name == GD_PARAB_BC ||
               BC_Types[bc_input_id].BC_Name == GD_PARAB_OFFSET_BC ||
               BC_Types[bc_input_id].BC_Name == GD_CIRC_BC ||
               BC_Types[bc_input_id].BC_Name == GD_TIME_BC ||
               BC_Types[bc_input_id].BC_Name == GD_POLYN_BC ||
               BC_Types[bc_input_id].BC_Name == GD_TABLE_BC) &&
              (BC_Types[bc_input_id].desc->equation == R_MESH_NORMAL ||
               BC_Types[bc_input_id].desc->equation ==
                   R_MESH2))) /* This is for the infamous GD_LINEAR trick */
                              /* && iapply   */
        ) {

          /* figure out which bc corresponds to this CA condition */
          j_bc_id = -1;
          for (j = 0; j < Num_BC; j++) {
            if (BC_Types[j].BC_Data_Int[0] == 1000000 * I + 1000000 &&
                ((BC_Types[j].BC_Data_Int[2] == -1 && iapply) ||
                 BC_Types[j].BC_Data_Int[2] == ei[pg->imtrx]->elem_blk_id)) {
              sum_normal = SQUARE(BC_Types[j].BC_Data_Float[1]) +
                           SQUARE(BC_Types[j].BC_Data_Float[2]) +
                           SQUARE(BC_Types[j].BC_Data_Float[3]);
              if (sum_normal) {
                /* if wall normals are input, use them for CA condition */
                variable_wall_normal = 0;
              } else {
                /* if no wall normals are input, calculate them for CA condition */
                variable_wall_normal = 1;
              }

              /*    determine which CA condition this is  */
              if (CA_id[0] == -1) {
                jcnt = 0;
                CA_id[0] = j;
              } else {
                jcnt = -1;
                for (p = 0; p < MAX_CA; p++) {
                  if (CA_id[p] == j) {
                    jcnt = p;
                  }
                }
                jflag = -1;
                for (p = 0; p < MAX_CA - 1; p++) {
                  if ((CA_id[p] != -1) && CA_id[p + 1] == -1) {
                    jflag = p + 1;
                  }
                }
                if (jcnt == -1) {
                  if (jflag == -1) {
                    GOMA_EH(GOMA_ERROR, "too many CA conditions - change MAX_CA\n");
                  } else {
                    jcnt = jflag;
                    CA_id[jcnt] = j;
                  }
                }
              }

              if (BC_Types[bc_input_id].BC_Name == KINEMATIC_BC ||
                  BC_Types[bc_input_id].BC_Name == KINEMATIC_COLLOC_BC ||
                  BC_Types[bc_input_id].BC_Name == KINEMATIC_PETROV_BC ||
                  BC_Types[bc_input_id].BC_Name == KINEMATIC_DISC_BC ||
                  BC_Types[bc_input_id].BC_Name == KIN_ELECTRODEPOSITION_BC || /*  RSL 5/27/02  */
                  BC_Types[bc_input_id].BC_Name == KIN_LEAK_BC) {

                CA_fselem[jcnt] = ielem;
                local_node_id = id;

                for (p = 0; p < ielem_dim; p++) {
                  fsnormal[jcnt][p] = fv->snormal[p];
                  for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                    for (q = 0; q < ielem_dim; q++) {
                      dfsnormal_dx[jcnt][p][q][k] = fv->dsnormal_dx[p][q][k];
                    }
                  }
                }

                /* if there is a constant normal read from the input deck */

                if (!(variable_wall_normal)) {
                  CA_sselem[jcnt] = ielem;
                  local_node_id = id;
                  for (p = 0; p < ielem_dim; p++) {
                    ssnormal[jcnt][p] = BC_Types[j].BC_Data_Float[p + 1];
                    for (q = 0; q < ielem_dim; q++) {
                      for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                        dssnormal_dx[jcnt][p][q][k] = 0.;
                      }
                    }
                  }
                }

              } else if ((variable_wall_normal) &&
                         (BC_Types[bc_input_id].BC_Name == PLANE_BC ||
                          BC_Types[bc_input_id].BC_Name == SPLINE_BC ||
                          BC_Types[bc_input_id].BC_Name == FILLET_BC ||
                          BC_Types[bc_input_id].BC_Name == DOUBLE_RAD_BC ||
                          BC_Types[bc_input_id].BC_Name == FEATURE_ROLLON_BC ||
                          BC_Types[bc_input_id].BC_Name == ROLL_FLUID_BC ||
                          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_BC ||
                          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_RS_BC ||
                          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_PETROV_BC ||
                          BC_Types[bc_input_id].BC_Name == KIN_DISPLACEMENT_COLLOC_BC ||
                          BC_Types[bc_input_id].BC_Name == MESH_CONSTRAINT_BC ||
                          BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC ||
                          BC_Types[bc_input_id].BC_Name == SOLID_FLUID_BC ||
                          BC_Types[bc_input_id].BC_Name == SOLID_FLUID_RS_BC ||
                          ((BC_Types[bc_input_id].BC_Name == GD_CONST_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_LINEAR_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_INVERSE_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_PARAB_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_PARAB_OFFSET_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_CIRC_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_TIME_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_POLYN_BC ||
                            BC_Types[bc_input_id].BC_Name == GD_TABLE_BC) &&
                           (BC_Types[bc_input_id].desc->equation == R_MESH_NORMAL ||
                            BC_Types[bc_input_id].desc->equation == R_MESH2)))) {

                /* However, we may have already gotten started, applied the
                   contact angle condition with say a KINEMATIC/GD_* combination,
                   and set side_info back to zero, only to increment here again
                   as we encounter another GD geometry condition on the same side
                   needed for a composite geometry.  This fouls up subsequent
                   entries in here for other contact angles. Hence only increment
                   if GD_count indicates zero or one geometry condition*/

                if (GD_count <= 1) {
                  CA_sselem[jcnt] = ielem;
                  local_node_id = id;
                }

                /* from calc_surf_normal we get the normal of the fluid instead of
                   the normal of the solid. We need to negate the normal and its
                   derivatives to get what we really want! */

                for (p = 0; p < ielem_dim; p++) {
                  ssnormal[jcnt][p] = -fv->snormal[p];
                  for (q = 0; q < ielem_dim; q++) {
                    for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                      dssnormal_dx[jcnt][p][q][k] = -fv->dsnormal_dx[p][q][k];
                    }
                  }
                }

              } /* endif BC_Types[bc_input_id].BC_Name == PLANE_BC || SPLINE_BC */
              j_bc_id = j;

              GOMA_EH(j_bc_id, "Contact angle condition not found at this node");
              /* First figure out where you are in the mesh */

              jflag = -1;
              for (p = 0; p < MAX_CA; p++) {
                if ((CA_fselem[p] == ielem || CA_sselem[p] == ielem) && CA_fselem[p] != -1 &&
                    CA_sselem[p] != -1 && CA_id[p] != -2)
                  jflag = p;
              }
              if (jflag != -1) {
                /* check for split element condition with frontal solver */
                if (Linear_Solver == FRONT && (CA_fselem[jflag] != CA_sselem[jflag])) {
                  GOMA_EH(GOMA_ERROR,
                          "Whoa!  frontal solver not equipped for split element CA conditions");
                }

                /*   make sure we are in the free surface element, etc.  */
                load_ei(CA_fselem[jflag], exo, 0, pg->imtrx);
                /*
                 * Load the field variable values at this node point
                 */
                find_nodal_stu(local_node_id, ielem_type, &xi[0], &xi[1], &xi[2]);

                err = load_basis_functions(xi, bfd);
                GOMA_EH(err, "problem from load_basis_functions");

                err = beer_belly();
                GOMA_EH(err, "beer_belly");

                err = load_fv();
                GOMA_EH(err, "load_fv");

                /* may want to load material properties here */

                /* calculate the determinant of the surface jacobian */
                surface_determinant_and_normal(
                    ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1, (int)elem_side_bc->id_side,
                    (int)elem_side_bc->num_nodes_on_side, (elem_side_bc->local_elem_node_id));

                if (ielem_dim != 3) {
                  calc_surf_tangent(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                    (int)elem_side_bc->num_nodes_on_side,
                                    (elem_side_bc->local_elem_node_id));
                }

                if (mp->SurfaceTensionModel != CONSTANT) {
                  load_surface_tension(dsigma_dx);
                  if (neg_elem_volume)
                    return (status);
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
                /* evaluate Gibbs inequality here with surface
                 * normal and solid surface normal */

                if (BC_Types[j_bc_id].BC_Name == CA_OR_FIX_BC) {
                  Gibbs = evaluate_gibbs_criterion(
                      &fsnormal[jflag][0], &ssnormal[jflag][0], &(BC_Types[j].BC_Data_Int[1]),
                      BC_Types[j].BC_Data_Float[0], BC_Types[j].BC_Data_Float[4],
                      BC_Types[j].BC_Data_Float[5], BC_Types[j].BC_Data_Float[6],
                      BC_Types[j].BC_Data_Float[7], BC_Types[j].BC_Data_Float[8],
                      BC_Types[j].BC_Data_Float[9]);

                  if (Gibbs) {
                    for (i1 = 0; i1 < MAX_PDIM; i1++) {
                      fsnrml[i1] = fsnormal[jflag][i1];
                      ssnrml[i1] = ssnormal[jflag][i1];
                      for (i2 = 0; i2 < MAX_PDIM; i2++) {
                        for (i3 = 0; i3 < MDE; i3++) {
                          dfsnrml_dx[i1][i2][i3] = dfsnormal_dx[jflag][i1][i2][i3];
                          dssnrml_dx[i1][i2][i3] = dssnormal_dx[jflag][i1][i2][i3];
                        }
                      }
                    }
                    fapply_CA(func, d_func, d_func_ss, fsnrml, dfsnrml_dx, ssnrml, dssnrml_dx,
                              BC_Types[j].BC_Data_Float[0]);
                  } else {
                    /* fix contact point where it sits right now */
                    for (p = 0; p < ielem_dim; p++) {
                      k = Index_Solution(I, MESH_DISPLACEMENT1 + p, 0, 0, -1, pg->imtrx);
                      if (!af->Assemble_LSA_Mass_Matrix) {
                        if (af->Assemble_Jacobian) {
                          d_func[p][MESH_DISPLACEMENT1 + p][id] = 1.;
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[MESH_DISPLACEMENT1 + p][id];
                          lec->R[LEC_R_INDEX(MESH_DISPLACEMENT1 + p, ldof_eqn)] = 0.;
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[MESH_DISPLACEMENT1 + p][id];
                          eqn = upd->ep[pg->imtrx][MESH_DISPLACEMENT1 + p];
                          zero_lec_row(lec->J, eqn, ldof_eqn);
                          lec->J[LEC_J_INDEX(eqn, eqn, ldof_eqn, ldof_eqn)] = DIRICHLET_PENALTY;
                        }
                      }
                      xdot[k] = 0.0;
                      xdot_old[k] = 0.0;
                      x_old[k] = x[k];
                      x_older[k] = x_old[k];
                    }
                  }
                }

                else if (BC_Types[j_bc_id].BC_Name == CA_BC ||
                         BC_Types[j_bc_id].BC_Name == CA_MOMENTUM_BC) {
                  /* Loop over the vector dimension i1 */
                  for (i1 = 0; i1 < MAX_PDIM; i1++) {
                    fsnrml[i1] = fsnormal[jflag][i1];
                    ssnrml[i1] = ssnormal[jflag][i1];
                    for (i2 = 0; i2 < MAX_PDIM; i2++) {
                      for (i3 = 0; i3 < MDE; i3++) {
                        dfsnrml_dx[i1][i2][i3] = dfsnormal_dx[jflag][i1][i2][i3];
                        dssnrml_dx[i1][i2][i3] = dssnormal_dx[jflag][i1][i2][i3];
                      }
                    }
                  }
                  fapply_CA(func, d_func, d_func_ss, fsnrml, dfsnrml_dx, ssnrml, dssnrml_dx,
                            BC_Types[j].BC_Data_Float[0]);
                } else if (BC_Types[j_bc_id].BC_Name == MOVING_CA_BC) {
                  for (i1 = 0; i1 < MAX_PDIM; i1++) {
                    fsnrml[i1] = fsnormal[jflag][i1];
                    ssnrml[i1] = ssnormal[jflag][i1];
                    for (i2 = 0; i2 < MAX_PDIM; i2++) {
                      for (i3 = 0; i3 < MDE; i3++) {
                        dfsnrml_dx[i1][i2][i3] = dfsnormal_dx[jflag][i1][i2][i3];
                        dssnrml_dx[i1][i2][i3] = dssnormal_dx[jflag][i1][i2][i3];
                      }
                    }
                  }
                  fapply_moving_CA(func, d_func, d_func_ss, fsnrml, dfsnrml_dx, ssnrml, dssnrml_dx,
                                   BC_Types[j].BC_Data_Float[0], BC_Types[j].BC_Data_Float[4],
                                   BC_Types[j].BC_Data_Float[5], BC_Types[j].BC_Data_Float[6],
                                   BC_Types[j].BC_Data_Float[7], BC_Types[j].BC_Data_Float[8],
                                   BC_Types[j].BC_Data_Float[9], x_dot, delta_t, theta);
                } else if (BC_Types[j_bc_id].BC_Name == VELO_THETA_TPL_BC ||
                           BC_Types[j_bc_id].BC_Name == VELO_THETA_HOFFMAN_BC ||
                           BC_Types[j_bc_id].BC_Name == VELO_THETA_COX_BC ||
                           BC_Types[j_bc_id].BC_Name == VELO_THETA_SHIK_BC) {
                  double wall_velocity = 0, velo[MAX_PDIM];
                  double dwall_velo_dx[DIM][MDE], dvelo_dx[DIM][DIM];
                  double lag_pars[2 * DIM + 1] = {0.0};
                  double t, axis_pt[DIM], R, rad_dir[DIM], v_dir[DIM];
                  int found_wall_velocity = 0, sfs_model = 0;
                  double theta_max = 180.0, dewet = 1.0, dcl_shearrate = -1.;

                  memset(dwall_velo_dx, 0, DIM * MDE * sizeof(double));

                  for (i1 = 0; i1 < MAX_PDIM; i1++) {
                    fsnrml[i1] = fsnormal[jflag][i1];
                    ssnrml[i1] = ssnormal[jflag][i1];
                    for (i2 = 0; i2 < MAX_PDIM; i2++) {
                      for (i3 = 0; i3 < MDE; i3++) {
                        dfsnrml_dx[i1][i2][i3] = dfsnormal_dx[jflag][i1][i2][i3];
                        dssnrml_dx[i1][i2][i3] = dssnormal_dx[jflag][i1][i2][i3];
                      }
                    }
                  }
                  /* try to find a VELO_TANGENT or similiar for wall velocity  */
                  for (i1 = 0; i1 < Num_BC; i1++) {

                    if ((BC_Types[i1].BC_Data_Int[0] == I &&
                         (BC_Types[i1].BC_Name == VELO_TANGENT_BC ||
                          BC_Types[i1].BC_Name == VELO_STREAMING_BC ||
                          BC_Types[i1].BC_Name == VELO_TANGENT_USER_BC ||
                          BC_Types[i1].BC_Name == VELO_SLIP_BC ||
                          BC_Types[i1].BC_Name == VELO_SLIP_ROT_BC ||
                          BC_Types[i1].BC_Name == VELO_SLIP_FLUID_BC ||
                          BC_Types[i1].BC_Name == VELO_SLIP_ROT_FLUID_BC ||
                          BC_Types[i1].BC_Name == AIR_FILM_BC ||
                          BC_Types[i1].BC_Name == AIR_FILM_ROT_BC)) ||
                        BC_Types[i1].BC_Name == VELO_TANGENT_SOLID_BC ||
                        (BC_Types[i1].BC_Data_Int[2] == I &&
                         BC_Types[i1].BC_Name == VELO_SLIP_SOLID_BC)) {
                      int mn;
                      switch (BC_Types[i1].BC_Name) {
                      case VELO_TANGENT_BC:
                      case VELO_STREAMING_BC:
                        wall_velocity = BC_Types[i1].BC_Data_Float[0] *
                                        (1.0 + BC_Types[i1].BC_Data_Float[3] * time_value);
                        found_wall_velocity = 1;
                        break;
                      case VELO_SLIP_BC:
                      case VELO_SLIP_FLUID_BC:
                      case AIR_FILM_BC:
                        velo[0] = BC_Types[i1].BC_Data_Float[1];
                        velo[1] = BC_Types[i1].BC_Data_Float[2];
                        wall_velocity = ssnrml[1] * velo[0] - ssnrml[0] * velo[1];
                        for (i2 = 0; i2 < ielem_dim; i2++) {
                          for (i3 = 0; i3 < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i3++) {
                            dwall_velo_dx[i2][i3] +=
                                dssnrml_dx[1][i2][i3] * velo[0] - dssnrml_dx[0][i2][i3] * velo[1];
                          }
                        }
                        found_wall_velocity = 1;
                        break;
                      case VELO_SLIP_ROT_BC:
                      case VELO_SLIP_ROT_FLUID_BC:
                      case AIR_FILM_ROT_BC:
                        wall_velocity = BC_Types[i1].BC_Data_Float[1] *
                                        (ssnrml[1] * (fv->x[1] - BC_Types[i1].BC_Data_Float[3]) -
                                         ssnrml[0] * (fv->x[0] - BC_Types[i1].BC_Data_Float[2]));
                        for (i2 = 0; i2 < ielem_dim; i2++) {
                          for (i3 = 0; i3 < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i3++) {
                            dwall_velo_dx[i2][i3] +=
                                BC_Types[i1].BC_Data_Float[1] *
                                (dssnrml_dx[1][i2][i3] *
                                     (fv->x[1] - BC_Types[i1].BC_Data_Float[3]) -
                                 dssnrml_dx[0][i2][i3] *
                                     (fv->x[0] - BC_Types[i1].BC_Data_Float[2]));
                          }
                        }
                        for (i3 = 0; i3 < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i3++) {
                          dwall_velo_dx[0][i3] += BC_Types[i1].BC_Data_Float[1] * (-ssnrml[0]) *
                                                  bf[MESH_DISPLACEMENT1]->phi[i3];
                          dwall_velo_dx[1][i3] += BC_Types[i1].BC_Data_Float[1] * (ssnrml[1]) *
                                                  bf[MESH_DISPLACEMENT1]->phi[i3];
                        }
                        found_wall_velocity = 1;
                        break;
                      case VELO_TANGENT_USER_BC:
                        velo[0] = velo_vary_fnc(UVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                BC_Types[i1].u_BC, time_value);
                        velo[1] = velo_vary_fnc(VVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                BC_Types[i1].u_BC, time_value);
                        wall_velocity = ssnrml[1] * velo[0] - ssnrml[0] * velo[1];
                        dvelo_dx[0][0] = dvelo_vary_fnc_d1(UVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                           BC_Types[i1].u_BC, time_value);
                        dvelo_dx[0][1] = dvelo_vary_fnc_d2(UVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                           BC_Types[i1].u_BC, time_value);
                        dvelo_dx[1][0] = dvelo_vary_fnc_d1(VVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                           BC_Types[i1].u_BC, time_value);
                        dvelo_dx[1][1] = dvelo_vary_fnc_d2(VVARY_BC, fv->x[0], fv->x[1], fv->x[2],
                                                           BC_Types[i1].u_BC, time_value);
                        for (i2 = 0; i2 < ielem_dim; i2++) {
                          for (i3 = 0; i3 < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i3++) {
                            dwall_velo_dx[i2][i3] +=
                                dssnrml_dx[1][i2][i3] * velo[0] - dssnrml_dx[0][i2][i3] * velo[1] +
                                (ssnrml[1] * dvelo_dx[0][i2] - ssnrml[0] * dvelo_dx[1][i2]) *
                                    bf[MESH_DISPLACEMENT1]->phi[i3];
                          }
                        }

                        found_wall_velocity = 1;
                        break;
                      case VELO_TANGENT_SOLID_BC:
                      case VELO_SLIP_SOLID_BC:
                        mn = map_mat_index(BC_Types[i1].BC_Data_Int[0]);
                        switch (pd_glob[mn]->MeshMotion) {
                        case LAGRANGIAN:
                          if (elc_glob[mn]->v_mesh_sfs_model == CONSTANT) {
                            sfs_model = CONSTANT;
                            for (i2 = 0; i2 < DIM; i2++) {
                              lag_pars[i2] = elc_glob[mn]->v_mesh_sfs[i2];
                            }
                          } else if (elc_glob[mn]->v_mesh_sfs_model == ROTATIONAL) {
                            sfs_model = ROTATIONAL;
                            for (i3 = 0; i3 <= DIM; i3++) {
                              lag_pars[i3] = elc_glob[mn]->u_v_mesh_sfs[i3];
                            }
                          } else if (elc_glob[mn]->v_mesh_sfs_model == ROTATIONAL_3D) {
                            sfs_model = ROTATIONAL_3D;
                            for (i3 = 0; i3 <= 2 * DIM; i3++) {
                              lag_pars[i3] = elc_glob[mn]->u_v_mesh_sfs[i3];
                            }
                          } else {
                            GOMA_EH(-1, "shouldn't be here\n");
                          }
                          break;
                        case TOTAL_ALE:
                          if (elc_rs_glob[mn]->v_mesh_sfs_model == CONSTANT) {
                            sfs_model = CONSTANT;
                            for (i2 = 0; i2 < DIM; i2++) {
                              lag_pars[i2] = elc_rs_glob[mn]->v_mesh_sfs[i2];
                            }
                          } else if (elc_rs_glob[mn]->v_mesh_sfs_model == ROTATIONAL) {
                            sfs_model = ROTATIONAL;
                            for (i3 = 0; i3 <= DIM; i3++) {
                              lag_pars[i3] = elc_rs_glob[mn]->u_v_mesh_sfs[i3];
                            }
                          } else if (elc_rs_glob[mn]->v_mesh_sfs_model == ROTATIONAL_3D) {
                            sfs_model = ROTATIONAL_3D;
                            for (i3 = 0; i3 <= 2 * DIM; i3++) {
                              lag_pars[i3] = elc_rs_glob[mn]->u_v_mesh_sfs[i3];
                            }
                          }
                          break;
                        }
                        switch (sfs_model) {
                        case CONSTANT:
                          for (i2 = 0; i2 < MAX_PDIM; i2++) {
                            wall_velocity += SQUARE(lag_pars[i2]);
                          }
                          wall_velocity = sqrt(wall_velocity);
                          break;
                        case ROTATIONAL:
                          for (i3 = 0; i3 < DIM; i3++) {
                            wall_velocity += SQUARE(fv->x[i3] - lag_pars[i3 + 1]);
                          }
                          wall_velocity = sqrt(wall_velocity) * lag_pars[0];

                          for (i3 = 0; i3 < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; i3++) {
                            dwall_velo_dx[0][i3] += SQUARE(lag_pars[0]) * (fv->x[0] - lag_pars[1]) *
                                                    bf[MESH_DISPLACEMENT1]->phi[i3] / wall_velocity;
                            dwall_velo_dx[1][i3] += SQUARE(lag_pars[0]) * (fv->x[1] - lag_pars[2]) *
                                                    bf[MESH_DISPLACEMENT1]->phi[i3] / wall_velocity;
                          }
                          break;
                        case ROTATIONAL_3D:
                          t = (lag_pars[4] * (fv->x0[0] - lag_pars[1]) +
                               lag_pars[5] * (fv->x0[1] - lag_pars[2]) +
                               lag_pars[6] * (fv->x0[2] - lag_pars[3])) /
                              (SQUARE(lag_pars[4]) + SQUARE(lag_pars[5]) + SQUARE(lag_pars[6]));
                          axis_pt[0] = lag_pars[1] + lag_pars[4] * t;
                          axis_pt[1] = lag_pars[2] + lag_pars[5] * t;
                          axis_pt[2] = lag_pars[3] + lag_pars[6] * t;
                          R = sqrt(SQUARE(fv->x0[0] - axis_pt[0]) + SQUARE(fv->x0[1] - axis_pt[1]) +
                                   SQUARE(fv->x0[2] - axis_pt[2]));
                          rad_dir[0] = (fv->x0[0] - axis_pt[0]) / R;
                          rad_dir[1] = (fv->x0[1] - axis_pt[1]) / R;
                          rad_dir[2] = (fv->x0[2] - axis_pt[2]) / R;
                          v_dir[0] = lag_pars[5] * rad_dir[2] - lag_pars[6] * rad_dir[1];
                          v_dir[1] = lag_pars[6] * rad_dir[0] - lag_pars[4] * rad_dir[2];
                          v_dir[2] = lag_pars[4] * rad_dir[1] - lag_pars[5] * rad_dir[0];
                          velo[0] = lag_pars[0] * R * v_dir[0];
                          velo[1] = lag_pars[0] * R * v_dir[1];
                          velo[2] = lag_pars[0] * R * v_dir[2];
                          for (i2 = 0; i2 < DIM; i2++) {
                            wall_velocity += SQUARE(velo[i2]);
                          }
                          wall_velocity = sqrt(wall_velocity);
                          break;
                        default:
                          wall_velocity = 0;
                        }
                        found_wall_velocity = 1;
                        break;
                      default:
                        GOMA_WH(-1, "Wall velocity bc not found\n");
                      } /* switch bc	*/
                    }   /*  if BC_Types  */
                  }     /*  Num_BC loop  */
                  if (!found_wall_velocity) {
                    wall_velocity = 0;
                    GOMA_WH(-1, "Wall velocity not found : setting to zero\n");
                  }
                  f = BC_Types[j].BC_Data_Float;
                  if (BC_Types[j].BC_Name == VELO_THETA_SHIK_BC) {
                    theta_max = f[8];
                  }
                  if (BC_Types[j].BC_Name == VELO_THETA_HOFFMAN_BC ||
                      BC_Types[j].BC_Name == VELO_THETA_COX_BC) {
                    theta_max = f[9];
                    dcl_shearrate = f[10];
                  }
                  if (BC_Types[j].BC_Name == VELO_THETA_HOFFMAN_BC ||
                      BC_Types[j].BC_Name == VELO_THETA_COX_BC ||
                      BC_Types[j].BC_Name == VELO_THETA_TPL_BC) {
                    dewet = f[8];
                  }
                  fapply_moving_CA_sinh(func, d_func, d_func_ss, fsnrml, dfsnrml_dx, ssnrml,
                                        dssnrml_dx, f[0], /* theta_0 */
                                        f[4],             /* v0 */
                                        f[5],             /* gamma/kT */
                                        f[6],             /* t_relax */
                                        f[7],             /* v_old */
                                        x_dot, delta_t, theta, time_value, wall_velocity, theta_max,
                                        dewet, dcl_shearrate, BC_Types[j_bc_id].BC_Name,
                                        dwall_velo_dx, local_node_id);
                } /* if VELO_THETA bc		*/
                else
                  GOMA_EH(GOMA_ERROR, "NO CA Condition applied ");
#if 0
		  load_ei(ielem, exo, 0, pg->imtrx);
            /*
	     * Load the field variable values at this node point
	     */
	    find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);

	    err = load_basis_functions( xi, bfd );
	    GOMA_EH( err, "problem from load_basis_functions");

	    err = beer_belly();
	    GOMA_EH( err, "beer_belly");
            
	    err = load_fv();
	    GOMA_EH( err, "load_fv");
            
	    /* may want to load material properties here */
	    
	    /* calculate the determinant of the surface jacobian */
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
	    
	    if (mp->SurfaceTensionModel != CONSTANT) 
              {
                load_surface_tension(dsigma_dx);
                if( neg_elem_volume ) return(status);
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
#endif
              }
            }
          }

        } /* end if (Variable_Mask[I].DBCA == 1 && ... */

        /*
         *
         *  START OF Surface Tangent CONDITIONS
         *
         *
         */
        if (node->DBST == 1 && iapply &&
            (BC_Types[bc_input_id].BC_Name == CAPILLARY_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_ROLL_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_USER_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_TABLE_BC ||
             BC_Types[bc_input_id].BC_Name == CAPILLARY_TABLE_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_RECOIL_PRESS_BC)) {
          j_bc_id = -1;
          for (j = 0; j < Num_BC; j++) {
            if (BC_Types[j].BC_Data_Int[0] == I &&
                (BC_Types[j].BC_Name == SURFTANG_BC || BC_Types[j].BC_Name == CAP_ENDFORCE_BC)) {
              j_bc_id = j;
            }
          }
          GOMA_EH(j_bc_id, "Surface Tangent condition not found at this node");
          jflag = 0; /* Sensitivies from this BC will be applied. */
          if (BC_Types[j_bc_id].BC_Name == SURFTANG_BC) {
            fapply_ST(func, d_func, BC_Types[j_bc_id].BC_Data_Float[0],
                      BC_Types[j_bc_id].BC_Data_Float[1], BC_Types[j_bc_id].BC_Data_Float[2],
                      BC_Types[j_bc_id].BC_Data_Float[3], rcoord, (int)elem_side_bc->id_side);
          } else {
            fapply_ST(func, d_func, BC_Types[j_bc_id].BC_Data_Float[0],
                      BC_Types[j_bc_id].BC_Data_Float[1], BC_Types[j_bc_id].BC_Data_Float[2],
                      BC_Types[j_bc_id].BC_Data_Float[3], rcoord, 0);
          }
        }

        if (node->DBSTS == 1 && iapply &&
            (BC_Types[bc_input_id].BC_Name == CAPILLARY_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_ROLL_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_USER_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_REPULSE_TABLE_BC ||
             BC_Types[bc_input_id].BC_Name == CAPILLARY_TABLE_BC ||
             BC_Types[bc_input_id].BC_Name == CAP_RECOIL_PRESS_BC)) {
          j_bc_id = -1;
          for (j = 0; j < Num_BC; j++) {
            if (BC_Types[j].BC_Data_Int[0] == I &&
                (BC_Types[j].BC_Name == SURFTANG_SCALAR_BC ||
                 BC_Types[j].BC_Name == CAP_ENDFORCE_SCALAR_BC)) {
              j_bc_id = j;
            }
          }
          GOMA_EH(j_bc_id, "Surface Tangent condition not found at this node");
          jflag = 0; /* Sensitivies from this BC will be applied. */
          if (BC_Types[j_bc_id].BC_Name == SURFTANG_SCALAR_BC) {
            apply_ST_scalar(func, d_func, (int)elem_side_bc->id_side,
                            BC_Types[j_bc_id].BC_Data_Float[0], BC_Types[j_bc_id].BC_Data_Int[1]);
          } else {
            apply_ST_scalar(func, d_func, 0, BC_Types[j_bc_id].BC_Data_Float[0],
                            BC_Types[j_bc_id].BC_Data_Int[1]);
          }
        }

        if (node->DBSES == 1 && iapply && (BC_Types[bc_input_id].BC_Name == TENSION_SHEET_BC)) {
          j_bc_id = -1;

          for (j = 0; j < Num_BC; j++) {
            if ((BC_Types[j].BC_Data_Int[0] == I) && (BC_Types[j].BC_Name == SHEET_ENDSLOPE_BC)) {
              j_bc_id = j;
            }
          }
          GOMA_EH(j_bc_id, "Surface Tangent condition not found at this node");

          apply_SES(func, d_func, elem_side_bc, BC_Types[bc_input_id].BC_Data_Float[0], I,
                    BC_Types[j_bc_id].BC_Data_Float[0], BC_Types[j_bc_id].BC_Data_Float[1], id,
                    iconnect_ptr, (int)elem_side_bc->id_side);

          jflag = 0; /* TAB:  this assignment is so the sensitivies from this BC will
                      * be applied.  Note it may be that SURFTANG_SCALAR, CAP_ENDFORCE might
                      * not be having their sensitivies added into the Jacobian.  If so
                      * that may be why they converge so poorly.
                      */
        }

        /*
         * !!!!!!!!!!!!!!! Add function into residual equations here !!!!!!!!!!!!!!!
         *        -> residual may be a vector of residuals
         */
        if (j_bc_id != -1) {
          for (p = 0; p < BC_Types[j_bc_id].desc->vector; p++) {
            /*
             * Check to see if this BC on this node is applicable
             *   (i.e. no other overriding Dirichlet conditions,
             * And find the global unknown number for applying this condition
             */
            if (jflag != -1 && (BC_Types[j_bc_id].BC_Name == CA_BC ||
                                BC_Types[j_bc_id].BC_Name == CA_MOMENTUM_BC ||
                                BC_Types[j_bc_id].BC_Name == CA_OR_FIX_BC ||
                                BC_Types[j_bc_id].BC_Name == VELO_THETA_TPL_BC ||
                                BC_Types[j_bc_id].BC_Name == VELO_THETA_HOFFMAN_BC ||
                                BC_Types[j_bc_id].BC_Name == VELO_THETA_COX_BC ||
                                BC_Types[j_bc_id].BC_Name == VELO_THETA_SHIK_BC ||
                                BC_Types[j_bc_id].BC_Name == MOVING_CA_BC)) {
              if (CA_fselem[jflag] != -1 && CA_sselem[jflag] != -1 && CA_id[jflag] == j_bc_id) {
                index_eq =
                    bc_eqn_index(id, I, j_bc_id, ei[pg->imtrx]->mn, p, &eqn, &matID_apply, &vd);
              } else {
                index_eq = -1;
              }
            } else {
              index_eq =
                  bc_eqn_index(id, I, j_bc_id, ei[pg->imtrx]->mn, p, &eqn, &matID_apply, &vd);
            }

            if (index_eq >= 0) {
              if (BC_Types[j_bc_id].BC_Name == SURFTANG_BC ||
                  BC_Types[j_bc_id].BC_Name == SURFTANG_SCALAR_BC ||
                  BC_Types[j_bc_id].BC_Name == CAP_ENDFORCE_BC ||
                  BC_Types[j_bc_id].BC_Name == CAP_ENDFORCE_SCALAR_BC)
                weight = 1.;
              if (BC_Types[j_bc_id].BC_Name == CA_BC ||
                  BC_Types[j_bc_id].BC_Name == CA_MOMENTUM_BC ||
                  BC_Types[j_bc_id].BC_Name == CA_OR_FIX_BC ||
                  BC_Types[j_bc_id].BC_Name == VELO_THETA_TPL_BC ||
                  BC_Types[j_bc_id].BC_Name == VELO_THETA_HOFFMAN_BC ||
                  BC_Types[j_bc_id].BC_Name == VELO_THETA_COX_BC ||
                  BC_Types[j_bc_id].BC_Name == VELO_THETA_SHIK_BC ||
                  BC_Types[j_bc_id].BC_Name == MOVING_CA_BC ||
                  BC_Types[j_bc_id].BC_Name == SHEET_ENDSLOPE_BC)
                weight = BIG_PENALTY;

              /*
               * Find the position in the local element residual to inject the
               * bc contribution into. Note, we always inject into the first
               * degree of freedom for a variable located at a node.
               */
              if (eqn == R_MASS) {
                ieqn = MAX_PROB_EQN + BC_Types[bc_input_id].species_eq;
              } else {
                ieqn = upd->ep[pg->imtrx][eqn];
              }
              ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
              lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] += weight * func[p];

              /*
               * add sensitivities into matrix
               *  - find index of sensitivity in matrix (if variable is not defined at this node,
               *      loop over all dofs in element)
               *  - add into matrix
               */

              if (af->Assemble_Jacobian) {

                /* now add in sensitivity of BC function */

                for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                  if (ei[pg->imtrx]->dof[var] && (BC_Types[j_bc_id].desc->sens[var] || 1) &&
                      jflag != -1)
                  /*			    if (pd->v[pg->imtrx][var] &&
                     (BC_Types[j_bc_id].desc->sens[var]
                     || 1 ) && jflag != -1) TAB 01-08-07 */
                  {
                    if (var != MASS_FRACTION) {
                      /* free surface normal derivatives :  if free surface element
                       * is not the same as current element, reload some element info
                       */
                      if (jflag != -1 && CA_fselem[jflag] != -1 && CA_fselem[jflag] != ielem &&
                          Linear_Solver != FRONT && CA_id[jflag] == j_bc_id && j_bc_id != -1) {
                        if (strcmp(Matrix_Format, "msr") != 0) {
                          GOMA_EH(GOMA_ERROR,
                                  "Unexpected matrix format in apply_bc_special, use msr");
                        }

                        load_ei(CA_fselem[jflag], exo, 0, pg->imtrx);

                        /* For nonlocal element information,
                         * we do a direct injection into a
                         * through a. TABAER BEWARE!!!
                         */
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          je = ei[pg->imtrx]->gun_list[var][j];
                          GOMA_EH(je, "Bad var index.");
                          ja = (index_eq == je)
                                   ? index_eq
                                   : in_list(je, ija[index_eq], ija[index_eq + 1], ija);
                          GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                          a[ja] += weight * d_func[p][var][j];
                        }
                      } else {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          pvar = upd->vp[pg->imtrx][var];
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                              weight * d_func[p][var][j];
                        }
                      }

                      /*    solid surface normal derivatives :  if solid surface element
                       *  is not the same as current element, reload some element info
                       */
                      if (jflag != -1 && CA_sselem[jflag] != -1 &&
                          CA_sselem[jflag] != ei[pg->imtrx]->ielem && Linear_Solver != FRONT &&
                          CA_id[jflag] == j_bc_id && j_bc_id != -1) {
                        if (strcmp(Matrix_Format, "msr") != 0) {
                          GOMA_EH(GOMA_ERROR,
                                  "Unexpected matrix format in apply_bc_special, use msr");
                        }

                        load_ei(CA_sselem[jflag], exo, 0, pg->imtrx);

                        /* For nonlocal element information, we do a
                         * direct injection into a through Jac_BC.
                         */
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          je = ei[pg->imtrx]->gun_list[var][j];
                          GOMA_EH(je, "Bad var index.");
                          ja = (index_eq == je)
                                   ? index_eq
                                   : in_list(je, ija[index_eq], ija[index_eq + 1], ija);
                          GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                          a[ja] += weight * d_func_ss[p][var][j];
                        }
                      } else {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          pvar = upd->vp[pg->imtrx][var];
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                              weight * d_func_ss[p][var][j];
                        }
                      }

                      /*    if changed elements, put element information back   */
                      if (ielem != ei[pg->imtrx]->ielem) {
                        load_ei(ielem, exo, 0, pg->imtrx);
                      }
                    } else /* variable type is MASS_FRACTION */
                    {
                      for (w = 0; w < pd->Num_Species_Eqn; w++) {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          // pvar = upd->vp[pg->imtrx][MAX_PROB_VAR + w];
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
                          lec->J[LEC_J_INDEX(ieqn, MAX_PROB_VAR + w, ldof_eqn, j)] +=
                              weight * d_func[p][MAX_PROB_VAR + w][j];
                        }
                      } /* end of loop over species */
                    }   /* end of if MASS_FRACTION */
                  }     /* end of variable exists and condition is sensitive to it */
                }       /* end of loop over variable types */
              }         /* end of NEWTON */

              if (jflag != -1 && CA_fselem[jflag] != -1 && CA_sselem[jflag] != -1 &&
                  CA_id[jflag] == j_bc_id && j_bc_id != -1) {
                CA_id[jflag] = -2;
              }

            } /* end of if Res_BC */
          }   /* end of loop over p for vector conditions */
        }     /* end of if j_bc_id != -1 */

      } /*End (if (CAPILLARY and KINEMATIC)) */
    }   /* end for (i=0; i< num_nodes_on_side; i++) */
  }     /*(end for ibc) */
  return status;
} /* END of apply_special_bc  */

int apply_shell_grad_bc(double x[],              /* Solution vector for the current processor */
                        double resid_vector[],   /* Residual vector for the current processor */
                        const double delta_t,    /* current time step size */
                        const double theta,      /* parameter (0 to 1) to vary time integration
                                                  *  ( implicit - 0 to explicit - 1)          */
                        const double h_elem_avg, /* global average element size */
                        const double h[DIM],     /* average element size */
                        const double mu_avg,     /* average element viscosity */
                        const double U_norm,     /* global velocity norm */
                        const int ielem,         /* element number */
                        const int ielem_type,    /* element type */
                        const int num_local_nodes,
                        const int ielem_dim,
                        const int iconnect_ptr,
                        ELEM_SIDE_BC_STRUCT *elem_side_bc,
                        /* Pointer to an element side boundary condition * structure */
                        const int num_total_nodes,
                        const int bc_application, /* flag indicating whether to integrate
                                                   * strong or weak BC's */
                        const double time_value,
                        const Exo_DB *exo)

/*************************************************************************
 *
 * apply_shell_grad_bc():
 *
 *    Calculate contributions to equations defined on shell element
 *    blocks which involve gradients of bulk field variables.
 *    In order to enable these sensitivities to be properly placed
 *    in the Jacobian, this part of the assembly must be done along
 *    with the bulk element(s).
 *    The BC assembly functions called from here load up local arrays
 *    local_r and local_j from the shell element Gauss points, after
 *    which they are summed into lec->R and lec->J. This function
 *    handles BC application types WEAK_SHELL_GRAD and STRONG_SHELL_GRAD
 *    only.
 *    Because these BC's are really parts of shell EQUATIONS, they
 *    are not handled the same way as other BC types. Some of the
 *    functions performed in other "apply_..._bc" functions are not
 *    present here, but could be invoked as necessary. These terms
 *    should normally be associated with the BOUNDARY etm for the
 *    equation.
 *
 *************************************************************************/
{
  double local_r[MAX_PROB_VAR + MAX_CONC][MAX_NODES_PER_SIDE];
  double local_j[MAX_PROB_VAR + MAX_CONC][MAX_PROB_VAR + MAX_CONC][MAX_NODES_PER_SIDE][MDE];
  double s, t, u, wt, xi[DIM], xi2[DIM];
  double coord[DIM];
  const int Nv = MAX_PROB_VAR + MAX_CONC;
  int ef, nf, elem1, eb1, ebid1, m1, elem2, eb2, ebid2, m2;
  int nbr_type, nbr_dim, ss_index, ibc, nbc;
  int eind, evec, pe, pv, ip, ip_total, i, ib, id, j, err;
  int id_side, nodes_per_side, bc_input_id;
  int idof, jdof, node, index, eqn, var;
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int shell_sv;
  int gnn_map[MAX_NODES_PER_SIDE];
  int *bulk_gnn_list = NULL, *n_dof = NULL;
  int dof_map[MDE], n_dofptr[MAX_VARIABLE_TYPES][MDE];
  int BCid;
  struct BC_descriptions *bc_desc;
  BOUNDARY_CONDITION_STRUCT *bc;
  ELEM_SIDE_BC_STRUCT *side = elem_side_bc;

  /* Basic data for bulk element */
  bc_input_id = 0;
  id_side = side->id_side;
  nbc = side->Num_BC;
  nodes_per_side = side->num_nodes_on_side;
  elem1 = ei[pg->imtrx]->ielem;
  eb1 = find_elemblock_index(ielem, exo);
  ebid1 = exo->eb_id[eb1];
  m1 = Matilda[eb1];

  /* Allocate the n_dof array and store DOF counts for bulk element */
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
  bulk_gnn_list = (int *)array_alloc(1, MAX_NODES_PER_SIDE, sizeof(int));
  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    n_dof[i] = ei[pg->imtrx]->dof[i];
  }

  /* Loop over friends of the current element */
  nf = num_elem_friends[ielem];
  for (ef = 0; ef < nf; ef++) {
    elem2 = elem_friends[ielem][ef];
    eb2 = find_elemblock_index(elem2, exo);
    ebid2 = exo->eb_id[eb2];
    m2 = Matilda[eb2];

    /* Determine if this element is bulk with shell friend(s) */
    nbr_type = Elem_Type(exo, elem2);
    nbr_dim = elem_info(NDIM, nbr_type);

    /* NOTE: this will not work for an edge of a 3D element yet */
    if ((ielem_dim - nbr_dim) == 1) {

      err = load_neighbor_var_data(elem1, elem2, n_dof, dof_map, n_dofptr, id_side, xi, exo);

      /*
       * DOF's may be indexed in a different order in the neighbor
       * element, so they will be matched by global node numbers.
       * Get ShapeVar for shell block, then determine LOCAL bulk
       * DOF nodes for it.
       */
      determine_ShapeVar(pd_glob[m2]);
      shell_sv = pd_glob[m2]->ShapeVar;

      /* Store these node numbers in bulk_gnn_list */
      for (i = 0; i < ei[pg->imtrx]->dof[shell_sv]; i++) {
        bulk_gnn_list[i] = ei[pg->imtrx]->gnn_list[shell_sv][i];
      }

      /* Assembly proceeds along SHELL element Gauss points */
      ip_total = elem_info(NQUAD, nbr_type);
      for (ip = 0; ip < ip_total; ip++) {

        /* Get basic data for the shell Gauss point */
        wt = Gq_weight(ip, nbr_type);
        find_stu(ip, nbr_type, &s, &t, &u);
        xi2[0] = s;
        xi2[1] = t;
        xi2[2] = u;

        /* First, setup shop at the bulk surface Gauss point */
        /* find_stu_on_bulk(id_side, ielem_dim, s, t, xi); */
        id = bulk_side_id_and_stu(elem1, elem2, xi2, xi, exo);
        if (id != id_side)
          GOMA_EH(GOMA_ERROR, "Bulk side ID mismatch!");
        setup_shop_at_point(elem1, xi, exo);

        /* Load surface determinant/normal data for bulk element side */
        err = get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
        surface_determinant_and_normal(elem1, ei[pg->imtrx]->iconnect_ptr, num_local_nodes,
                                       ielem_dim - 1, id_side, nodes_per_side, local_elem_node_id);

        if (ielem_dim != 3) {
          calc_surf_tangent(elem1, ei[pg->imtrx]->iconnect_ptr, num_local_nodes, ielem_dim - 1,
                            nodes_per_side, local_elem_node_id);
        }
        for (i = 0; i < ielem_dim; i++) {
          coord[i] = fv->x[i];
        }

        /* Then, setup shop on the shell element Gauss point */
        setup_shop_at_point(elem2, xi2, exo);

        /* Load mapping array gnn_map (shell DOF -> bulk DOF) */
        for (i = 0; i < ei[pg->imtrx]->dof[shell_sv]; i++) {
          node = ei[pg->imtrx]->gnn_list[shell_sv][i];
          index = in_list(node, 0, ei[pg->imtrx]->dof[shell_sv], bulk_gnn_list);
          GOMA_EH(index, "Mapping fault!");
          gnn_map[i] = index;
        }

        /* This is for error checking */
        err = 0;
        if (ei[pg->imtrx]->num_local_nodes != nodes_per_side)
          GOMA_EH(GOMA_ERROR, "Node number mismatch between elements on this side!");
        if (ei[pg->imtrx]->num_local_nodes != ei[pg->imtrx]->dof[pd->ShapeVar])
          GOMA_EH(GOMA_ERROR, "Cannot handle current shell ShapeVar!");

        do_LSA_mods(LSA_SURFACE);

        /* Loop over defined boundary conditions */
        for (ibc = 0; ibc < nbc; ibc++) {
          BCid = (int)side->BC_input_id[ibc];
          bc = BC_Types + BCid;
          bc_desc = bc->desc;
          if ((ss_index = in_list(bc->BC_ID, 0, exo->num_side_sets, &(ss_to_blks[0][0]))) == -1) {
            sprintf(Err_Msg, "Could not find BC_ID %d in ss_to_blks", BC_Types[bc_input_id].BC_ID);
            GOMA_EH(GOMA_ERROR, Err_Msg);
          }

          /*
           * Check prescribed application method and proceed
           * only for SHELL_GRAD type BC's, if the two element
           * block ID numbers match those from the input card.
           * BC_Data_Int[0] entry is bulk block ID
           * BC_Data_Int[1] entry is shell block ID
           * NOTE: When a shell equation requires contributions
           *       from two bulk blocks, separate BC cards
           *       are needed for each one!
           */
          if ((bc_desc->method == bc_application) && BCid != -1 && (bc->BC_Data_Int[0] == ebid1) &&
              (bc->BC_Data_Int[1] == ebid2)) {
            memset(local_r, 0, Nv * MAX_NODES_PER_SIDE * sizeof(double));
            memset(local_j, 0, Nv * Nv * MAX_NODES_PER_SIDE * MDE * sizeof(double));
            /* Default values for eind and evec */
            eind = 0;
            evec = MAX_VARIABLE_TYPES;

            /*
             * Call the appropriate "SHELL_GRAD" assembly function.
             * Also provide the following information to identify
             * the "real" shell equation:
             * eind: Equation (R_###) (first component if vector)
             * evec: Number of vector components (default = 1)
             * Otherwise, the defaults above will be used.
             */
            switch (bc->BC_Name) {

            case SURFACE_ELECTRIC_FIELD_BC:

              eind = R_SURF_CHARGE;
              evec = 1;
              surface_electric_field_bc(local_r, local_j, m1, n_dof, wt);
              break;

            case SURFACE_ACOUSTIC_VELOCITY_BC:

              eind = R_SHELL_BDYVELO;
              evec = 1;
              surface_acoustic_velocity_bc(local_r, local_j, m1, n_dof, wt, time_value);
              break;

            case SURFACE_USER_SHELL_BC:

              eind = R_SHELL_USER;
              evec = 1;
              surface_user_shell_bc(local_r, local_j, m1, n_dof, wt, theta, delta_t, coord);
              break;

            case SURFACE_LUBRICATION_BC:

              eind = R_SHELL_LUBP;
              evec = 1;
              surface_lubrication_shell_bc(local_r, local_j, m1, n_dof, wt, theta, delta_t, coord,
                                           dof_map, n_dofptr);
              break;

            default:
              /* Do nothing */
              break;
            }

            /* Now transfer these terms to lec->R and lec->J */
            /* This process takes several steps */

            /* Loop over applicable shell equations */
            for (eqn = eind; eqn < (eind + evec); eqn++) {
              pe = upd->ep[pg->imtrx][eqn];
              idof = ei[pg->imtrx]->dof[eqn];
              for (i = 0; i < idof; i++) {
                ib = gnn_map[i];

                /* Transfer residual term */
                lec->R[LEC_R_INDEX(pe, ib)] += local_r[pe][i];

                /* Loop through and transfer sensitivities */
                for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                  pv = upd->vp[pg->imtrx][var];
                  jdof = n_dof[var];
                  for (j = 0; j < jdof; j++) {
                    lec->J[LEC_J_INDEX(pe, pv, ib, j)] += local_j[pe][pv][i][j];
                  }
                } /* End of term transfer */
              }   /* End of ei[pg->imtrx]->dof[eqn] loop */
            }     /* End of eqn loop */
          }       /* End of applicable BC if block */
        }         /* End of loop over all BC's */
      }           /* End of loop over shell Gauss points */
    }             /* End of (ielem_dim - nbr_dim) == 1 if block */
  }               /* End of loop over friends of bulk element */

  /* Cleanup: return to bulk element */
  err = load_elem_dofptr(elem1, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
  safe_free((void *)n_dof);
  safe_free((void *)bulk_gnn_list);

  return err;
} /* END of function apply_shell_grad_bc() */

int apply_sharp_integrated_bc(double x[],            /* Solution vector for the current processor */
                              double resid_vector[], /* Residual vector for the current processor */
                              const double time,     /* current time  */
                              const double delta_t,  /* current time step size */
                              const double theta,    /* parameter (0 to 1) to vary time integration
                                                      *  ( implicit - 0 to explicit - 1)          */
                              const double hsquared[DIM],
                              const int ielem,      /* element number */
                              const int ielem_type, /* element type */
                              const int num_local_nodes,
                              const int ielem_dim,
                              const int iconnect_ptr,
                              ELEM_SIDE_BC_STRUCT *elem_side_bc,
                              /* Pointer to an element side boundary condition * structure */
                              const int bc_application,
                              const Exo_DB *exo)

/*
 * apply_sharp_integrated_bc
 *
 *   The task assigned to this routine is to apply conditions at a point where the
 *   level set surface intersects the boundary.  The path dependencies are computed
 *   via finite differences.
 */
{
  int ip, ip_total;
  double *wt, weight;
  double(*s)[DIM];
  double xi[DIM];

  double fmin, fmax, fplus, fminus;
  double F;
  double h_elem;

  int a, i, id;

  /* calc fplus and fminus for finite difference */
  h_elem = 0.;
  for (a = 0; a < ielem_dim; a++) {
    h_elem += hsquared[a];
  }
  /* This is the average value of h**2 in the element */
  h_elem = h_elem / ((double)ielem_dim);

  /* This is the size of the element */
  h_elem = sqrt(h_elem);

  fmin = 0.;
  fmax = 0.;
  for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {
    id = (int)elem_side_bc->local_elem_node_id[i];
    F = lnn_distance(id);
    if (F < fmin)
      fmin = F;
    if (F > fmax)
      fmax = F;
  }

  fplus = 0.8 * fmax;
  if (fplus > FD_FACTOR * h_elem)
    fplus = FD_FACTOR * h_elem;
  fminus = 0.8 * fmin;
  if (fminus < -FD_FACTOR * h_elem)
    fminus = -FD_FACTOR * h_elem;

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

  if (!ls->Ignore_F_deps) {
    /* finite difference evaluation at fminus */
    iso_contour_on_side(fminus, pd->Num_Dim, ielem_type, elem_side_bc->id_side, &ip_total, &s, &wt);

    for (ip = 0; ip < ip_total; ip++) {
      xi[0] = s[ip][0];
      xi[1] = s[ip][1];
      xi[2] = s[ip][2];
      weight = wt[ip] / (fplus - fminus);

      assemble_sharp_integrated_bc(x, resid_vector, time, delta_t, theta, ielem, ielem_type,
                                   num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                   bc_application, exo, xi, weight, TRUE);
    }

    /* finite difference evaluation at fplus */
    iso_contour_on_side(fplus, pd->Num_Dim, ielem_type, elem_side_bc->id_side, &ip_total, &s, &wt);

    for (ip = 0; ip < ip_total; ip++) {
      xi[0] = s[ip][0];
      xi[1] = s[ip][1];
      xi[2] = s[ip][2];
      weight = -wt[ip] / (fplus - fminus);

      assemble_sharp_integrated_bc(x, resid_vector, time, delta_t, theta, ielem, ielem_type,
                                   num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                   bc_application, exo, xi, weight, TRUE);
    }
  }

  /* Finally, evaluate bc at interface */
  iso_contour_on_side(0., pd->Num_Dim, ielem_type, elem_side_bc->id_side, &ip_total, &s, &wt);

  for (ip = 0; ip < ip_total; ip++) {
    xi[0] = s[ip][0];
    xi[1] = s[ip][1];
    xi[2] = s[ip][2];
    weight = wt[ip];

    assemble_sharp_integrated_bc(x, resid_vector, time, delta_t, theta, ielem, ielem_type,
                                 num_local_nodes, ielem_dim, iconnect_ptr, elem_side_bc,
                                 bc_application, exo, xi, weight, FALSE);
  }
  return (1);
}

void assemble_sharp_integrated_bc(
    double x[],            /* Solution vector for the current processor */
    double resid_vector[], /* Residual vector for the current processor */
    const double time,     /* current time  */
    const double delta_t,  /* current time step size */
    const double theta,    /* parameter (0 to 1) to vary time integration
                            *  ( implicit - 0 to explicit - 1)          */
    const int ielem,       /* element number */
    const int ielem_type,  /* element type */
    const int num_local_nodes,
    const int ielem_dim,
    const int iconnect_ptr,
    ELEM_SIDE_BC_STRUCT *elem_side_bc,
    /* Pointer to an element side boundary condition structure */
    const int bc_application,
    const Exo_DB *exo,
    double xi[DIM],
    double weight,
    int CalcPathDeps) {
  int err;
  int id, i, j, I, p;
  int index_eq, matID_apply, eqn, ldof_eqn;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  int ibc, bc_input_id;
  BOUNDARY_CONDITION_STRUCT *bc;
  struct BC_descriptions *bc_desc;

  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];

  double wt;

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

  /* calculate the determinant of the surface jacobian and the normal to
   * the surface all at one time */
  surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                 (int)elem_side_bc->id_side, (int)elem_side_bc->num_nodes_on_side,
                                 (elem_side_bc->local_elem_node_id));

  err = load_fv_grads();
  GOMA_EH(err, "load_fv_grads");

  err = load_fv_mesh_derivs(1);
  GOMA_EH(err, "load_fv_mesh_derivs");

  /*
   * Load up commonly used physical properties such as density at
   * the current quadrature point using the material state vector.
   */
  load_properties(mp, time);

  load_lsi(0.);
  load_lsi_derivs();

  ibc = 0;

  while ((bc_input_id = (int)elem_side_bc->BC_input_id[ibc]) != -1) {
    bc = BC_Types + bc_input_id;
    bc_desc = bc->desc;
    memset(func, 0, DIM * sizeof(double));
    memset(d_func, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

    switch (bc->BC_Name) {
    case SHARP_CA_2D_BC:
      apply_sharp_ca(func, d_func, BC_Types[bc_input_id].BC_Data_Float[0]);
      break;
    case SHARP_WETLIN_VELOCITY_BC:
      apply_sharp_wetting_velocity(
          func, d_func, bc->BC_Name, BC_Types[bc_input_id].BC_Data_Float[3],
          BC_Types[bc_input_id].BC_Data_Float[0], BC_Types[bc_input_id].BC_Data_Float[1], 0, 0, 0);
      break;
    case SHARP_BLAKE_VELOCITY_BC:
      apply_sharp_wetting_velocity(
          func, d_func, bc->BC_Name, BC_Types[bc_input_id].BC_Data_Float[3],
          BC_Types[bc_input_id].BC_Data_Float[0], BC_Types[bc_input_id].BC_Data_Float[2],
          BC_Types[bc_input_id].BC_Data_Float[1], BC_Types[bc_input_id].BC_Data_Float[4],
          BC_Types[bc_input_id].BC_Data_Float[5]);
      break;
    case SHARP_HOFFMAN_VELOCITY_BC:
      apply_sharp_wetting_velocity(
          func, d_func, bc->BC_Name, BC_Types[bc_input_id].BC_Data_Float[2],
          BC_Types[bc_input_id].BC_Data_Float[0], BC_Types[bc_input_id].BC_Data_Float[1], 0,
          BC_Types[bc_input_id].BC_Data_Float[3], BC_Types[bc_input_id].BC_Data_Float[4]);
      break;
    case SHARP_COX_VELOCITY_BC:
    case SHARP_SHIK_VELOCITY_BC:
      apply_sharp_wetting_velocity(
          func, d_func, bc->BC_Name, BC_Types[bc_input_id].BC_Data_Float[3],
          BC_Types[bc_input_id].BC_Data_Float[0], BC_Types[bc_input_id].BC_Data_Float[1],
          BC_Types[bc_input_id].BC_Data_Float[2], BC_Types[bc_input_id].BC_Data_Float[4],
          BC_Types[bc_input_id].BC_Data_Float[5]);
      break;

    default:
      break;
    }

    for (i = 0; i < (int)elem_side_bc->num_nodes_on_side; i++) {

      id = (int)elem_side_bc->local_elem_node_id[i];

      I = Proc_Elem_Connect[iconnect_ptr + id];

      for (p = 0; p < bc_desc->vector; p++) {

        index_eq = bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn, p, &eqn, &matID_apply, &vd);

        if (index_eq >= 0) {
          double phi_i;
          int ieqn;

          if ((ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id]) != -1) {

            phi_i = bf[eqn]->phi[ldof_eqn];
            wt = weight * phi_i * fv->sdet;
            /* DRN: I am a little confused, but I think that fv->sdet should not be used
               here though since we have transformed the 2D integral down to 1D.
            */
            /*wt = weight*phi_i*fv->sdet;*/
            wt = weight * phi_i;

            if (bc_desc->method == STRONG_SHARP_INT)
              wt *= BIG_PENALTY;

            if (eqn == R_MASS)
              ieqn = MAX_PROB_EQN + bc->species_eq;
            else
              ieqn = upd->ep[pg->imtrx][eqn];

            /* primary piece of path dependence from finite difference */
            if (CalcPathDeps) {
              if (af->Assemble_Jacobian) {
                int var = LS;
                int pvar = upd->vp[pg->imtrx][var];
                if (pvar != -1) {
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] += wt * func[p] * bf[var]->phi[j];
                  }
                }
              }
            } else {
              lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] += wt * func[p];

              if (af->Assemble_Jacobian) {
                int pvar, var, w;

                for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                  if ((pvar = upd->vp[pg->imtrx][var]) != -1) {
                    if (var != MASS_FRACTION) {
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] += wt * d_func[p][var][j];
                      }
                    } else {
                      for (w = 0; w < pd->Num_Species_Eqn; w++) {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          lec->J[LEC_J_INDEX(ieqn, MAX_PROB_VAR + w, ldof_eqn, j)] +=
                              wt * d_func[p][MAX_VARIABLE_TYPES + w][j];
                        }
                      }
                    }
                  }
                }
                /* last piece of path dependence */
                var = LS;
                pvar = upd->vp[pg->imtrx][var];
                if (pvar != -1 && !ls->Ignore_F_deps) {
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                        wt * func[p] * lsi->d_gfmag_dF[j] * lsi->gfmaginv;
                  }
                }
              } /* end of if( af->Assemble_Jacobian */
            }   /* end of if ( CalcSurfDeps */
          }     /* end of if ( (ldof_eqn = */
        }       /* end of if( index_eq >= */
      }         /* end of for ( p=0 ; ... */
    }           /* end of for( i=0 ... */
    ibc++;
  }
}

/*****************************************************************************/
