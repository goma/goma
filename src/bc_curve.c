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
 *$Id: bc_curve.c,v 5.5 2010-04-05 15:04:46 hkmoffa Exp $
 */

/* Standard include files */

#include <stdio.h>
#include <string.h>

/* GOMA include files */

#include "ac_stability.h"
#include "ac_stability_util.h"
#include "bc_colloc.h"
#include "bc_curve.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "load_field_variables.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_porous.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_vars_const.h"
#include "std.h"
#include "user_bc.h"

#define GOMA_BC_CURVE_C

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */

/*
 * Prototype declarations of static variables and functions for this file.
 */

/* End GOMA include files */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int apply_integrated_curve_bc(
    double x[],            /* Solution vector for the current processor */
    double resid_vector[], /* Residual vector for the current processor */
    const double delta_t,  /* current time step size                    */
    const double theta,    /* parameter to vary time integration from
                            * explicit (theta=1) to implicit (theta=0) */
    const int ielem,       /* element number */
    const int ielem_type,  /* element type */
    const int num_local_nodes,
    const int ielem_dim,
    const int iconnect_ptr,
    struct elem_edge_bc_struct *elem_edge_bc, /* Pointer to an element side
                                               * boundary condition structure */
    const int num_total_nodes,
    const int bc_application, /* flag indicating whether to integrate
                               * strong or weak BC's */
    const Exo_DB *exo)        /* ptr to FE database */

/************************************************************************
 *
 * apply_integrated_curve_bc():
 *
 *
 *  This function evaluates weakly integrated boundary conditions applied
 *  along curves in 3d (i.e., intersections of two side sets.
 ************************************************************************/
{
  int ip, w, i, I, k, j, id, icount, ss_index, type, mn, matID_apply;
  int param_dir;
  int eqn, ieqn, var, pvar, p, q, index_eq, ldof_eqn;
  int err; /* status variable for functions */
  int status = 0;
  int bc_input_id, BC_Name, ip_total;

  double s; /* Gaussian-quadrature point locations          */

  double phi_i;
  double xi[DIM]; /* Local element coordinates of Gauss point. */
  double x_dot[MAX_PDIM];
  /****************************************************************************/
  double wt; /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  double weight;
  double dsigma_dx[DIM][MDE];
  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double d_func_ss[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double theta_r; /* contact angle in radians */
  VARIABLE_DESCRIPTION_STRUCT *vd;
  /* Normals for variable normal contact angle condition etc.                     */
  static double fsnormal[DIM]; /* Free surface normal component vector   */
  static double dfsnormal_dx[DIM][DIM][MDE];
  /* Free surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j        */

  static double ssnormal[DIM]; /* Solid surface normal component vector  */
  static double dssnormal_dx[DIM][DIM][MDE];
  /* Solid surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j        */
  static double clnormal[MAX_PDIM]; /* contact line normal component vector  */
  static double dclnormal_dx[MAX_PDIM][MAX_PDIM][MDE];
  /* contact line normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j        */

  /***************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)        */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                      */
  /***************************************************************************/
  /* Find out the number of surface quadrature points
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_EDGE, ielem_type);

  /* Surface integration over element */

  for (ip = 0; ip < ip_total; ip++) {
    /* find the quadrature point locations for current ip */
    param_dir = find_edge_s(ip, ielem_type, elem_edge_bc->id_edge, pd->Num_Dim, xi, &s);

    /* find the quadrature weight for current ip */
    wt = Gq_edge_weight(ip, ielem_type);

    /* ****************************************/
    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "problem from load_basis_functions");

    err = beer_belly();
    GOMA_EH(err, "beer_belly");

    /* precalculate variables at  current integration pt.*/
    err = load_fv();
    GOMA_EH(err, "load_fv");

    err = load_bf_grad();
    GOMA_EH(err, "load_bf_grad");

    err = load_bf_mesh_derivs();
    GOMA_EH(err, "load_bf_mesh_derivs");

    /* use primary side to find edge vectors */
    edge_determinant_and_vectors(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                 (int)elem_edge_bc->elem_side_bc_1->id_side,
                                 (int)elem_edge_bc->elem_side_bc_1->num_nodes_on_side,
                                 elem_edge_bc->elem_side_bc_1->local_elem_node_id,
                                 /* id_side, num_nodes_on_side, local_elem_node_id,  */
                                 (int)elem_edge_bc->id_edge, (int)elem_edge_bc->num_nodes_on_edge,
                                 elem_edge_bc->edge_elem_node_id, param_dir);

    /*
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    err = load_fv_grads();
    GOMA_EH(err, "load_fv_grads");

    err = load_fv_mesh_derivs(1);
    GOMA_EH(err, "load_fv_mesh_derivs");

    /*
     * Load up porous media variables and properties, if needed
     */
    if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED ||
        mp->PorousMediaType == POROUS_TWO_PHASE) {
      err = load_porous_properties();
      GOMA_EH(err, "load_porous_properties");
    }

    if (mp->SurfaceTensionModel != CONSTANT) {
      load_surface_tension(dsigma_dx);
      if (neg_elem_volume)
        return (status);
    }

    if (TimeIntegration != STEADY && pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (icount = 0; icount < ielem_dim; icount++) {
        x_dot[icount] = fv_dot->x[icount];
        /* calculate surface position for wall repulsion/no penetration condition */
      }
    } else {
      for (icount = 0; icount < ielem_dim; icount++) {
        x_dot[icount] = 0.;
      }
    }

    do_LSA_mods(LSA_EDGE);

    if ((bc_input_id = (int)elem_edge_bc->BC_input_id[0]) != -1) {

      if ((ss_index = in_list(BC_Types[bc_input_id].BC_ID, 0, Proc_Num_Side_Sets, ss_to_blks[0])) ==
          -1) {
        GOMA_EH(GOMA_ERROR, "Cannot match side set id with that in ss_to_blks array");
      }

      /* check to see if this bc is an integrated bc */

      if (BC_Types[bc_input_id].desc->method == bc_application) {

        /* Get the initialization out of the way, for safty */
        memset(func, 0, DIM * sizeof(double));
        memset(d_func, 0, DIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

        memset(d_func_ss, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

        /* initialize the general function to zero may have more than one entry
         * for vector conditions like capillary */

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

        switch (BC_Name = BC_Types[bc_input_id].BC_Name) {

        case SURFTANG_EDGE_BC:
          /* calculate the tension exerted through the edge in the
           * direction of the binormal  -  binormal is stangent[0]
           */
          apply_ST_3D(func, d_func, BC_Types[bc_input_id].BC_Data_Float[0],
                      BC_Types[bc_input_id].BC_Data_Float[1],
                      BC_Types[bc_input_id].BC_Data_Float[2],
                      BC_Types[bc_input_id].BC_Data_Float[3]);
          break;

        case SURFTANG_SCALAR_EDGE_BC:
          /* calculate the tension exerted through the edge in the
           *  direction of the binormal
           *  - binormal is stangent[0]
           */
          apply_ST_scalar_3D(func, d_func, BC_Types[bc_input_id].BC_Data_Float[0]);
          break;

        case CA_EDGE_CURVE_INT_BC:

          if (elem_edge_bc->shared) {
            GOMA_EH(GOMA_ERROR, "CA_EDGE_CURVE_INT cannot be used with shared edges.");
          }
          break;

        case CA_EDGE_INT_BC:
          /* calculate contact angle boundary condition
           *     use first data entry as contact angle in degrees */
          theta_r = BC_Types[bc_input_id].BC_Data_Float[0] * M_PIE / 180.;

          /* collect the normals */
          /* Free surface normal */
          for (p = 0; p < ielem_dim; p++) {
            fsnormal[p] = fv->snormal[p];
            for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
              for (q = 0; q < ielem_dim; q++) {
                dfsnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
              }
            }
          }
          /* Wall surface normal */
          if (BC_Types[bc_input_id].BC_Name == CA_EDGE_INT_BC) {
            for (p = 0; p < ielem_dim; p++) {
              ssnormal[p] = BC_Types[bc_input_id].BC_Data_Float[p + 1];
              for (q = 0; q < ielem_dim; q++) {
                for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                  dssnormal_dx[q][p][k] = 0.;
                }
              }
            }
          } else if (BC_Types[bc_input_id].BC_Name == CA_EDGE_CURVE_INT_BC) {
            /* recompute surface normal based upon secondary sideset. */

            surface_determinant_and_normal(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                           (int)elem_edge_bc->elem_side_bc_2->id_side,
                                           (int)elem_edge_bc->elem_side_bc_2->num_nodes_on_side,
                                           elem_edge_bc->elem_side_bc_2->local_elem_node_id);

            for (p = 0; p < ielem_dim; p++) {
              ssnormal[p] = fv->snormal[p];
              for (q = 0; q < ielem_dim; q++) {
                for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                  dssnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
                }
              }
            }
          }

          fapply_CA(func, d_func, d_func_ss, fsnormal, dfsnormal_dx, ssnormal, dssnormal_dx,
                    theta_r);

          break;

        case VELO_NORMAL_EDGE_INT_BC:
        case KINEMATIC_EDGE_BC:

          fvelo_normal_edge_bc(func, d_func, BC_Types[bc_input_id].BC_Data_Float[0], x_dot, theta,
                               delta_t);
          break;

        case VELO_TANGENT_EDGE_INT_BC:

          fvelo_tangent_edge_bc(func, d_func, BC_Types[bc_input_id].BC_Data_Float, x_dot, theta,
                                delta_t);
          break;

        case VAR_CA_EDGE_BC:
        case VAR_CA_USER_BC:

          /* free surface normal and sensitivities.
           */

          for (p = 0; p < ielem_dim; p++) {
            fsnormal[p] = fv->snormal[p];
            for (q = 0; q < ielem_dim; q++) {
              for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                dfsnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
              }
            }
          }

          /* Wall surface normal *
           * for now assume that the contact angle condition is applied at a
           * fixed wall - where the wall surface normal is known */

          for (p = 0; p < ielem_dim; p++) {
            if (BC_Name == VAR_CA_EDGE_BC) {
              ssnormal[p] = BC_Types[bc_input_id].BC_Data_Float[p + 5];
            } else {
              ssnormal[p] = BC_Types[bc_input_id].u_BC[3 + p];
            }
          }

          memset(dssnormal_dx, 0, sizeof(double) * DIM * DIM * MDE);

          calc_CL_normal(ssnormal, dssnormal_dx, fsnormal, dfsnormal_dx, fv->stangent[1],
                         fv->dstangent_dx[1], ielem, elem_edge_bc->edge_elem_node_id, ielem_dim,
                         (int)elem_edge_bc->num_nodes_on_edge, clnormal, dclnormal_dx, exo);

          if (BC_Types[bc_input_id].BC_Name == VAR_CA_EDGE_BC) {
            fapply_var_CA(func, d_func, d_func_ss, fsnormal, dfsnormal_dx, ssnormal, dssnormal_dx,
                          clnormal, dclnormal_dx, BC_Types[bc_input_id].BC_Data_Float, x_dot, theta,
                          delta_t);
          } else if (BC_Types[bc_input_id].BC_Name == VAR_CA_USER_BC) {
            fapply_var_CA_user(func, d_func, d_func_ss, fsnormal, dfsnormal_dx, ssnormal,
                               dssnormal_dx, clnormal, dclnormal_dx, BC_Types[bc_input_id].len_u_BC,
                               BC_Types[bc_input_id].u_BC, x_dot, theta, delta_t);
          }
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Integrated BC not found");
          break;
        } /* end of switch over bc type */

        /******************************************************************************/
        /*                LOOP OVER THE LOCAL ELEMENT NODE NUMBER                     */
        /*               FOR NODES THAT LIE ON THE CURRENT SURFACE		      */
        /* this is a loop designed to loop over equations that lie on current surface */
        /*         ADD the Boundary condition functions into the Residual             */
        /*         vector and Jacobian Matrix                                         */
        /******************************************************************************/

        for (i = 0; i < (int)elem_edge_bc->num_nodes_on_edge; i++) {

          /* Find the local element node number for the current node */
          id = (int)elem_edge_bc->edge_elem_node_id[i];

          /* Find the local node number given the local element node number,  'i'     */
          I = Proc_Elem_Connect[iconnect_ptr + id];

          /* add function into residual equations here may a vector of functions */

          for (p = 0; p < BC_Types[bc_input_id].desc->vector; p++) {
            /*
             * Check to see if this BC on this node is applicable
             * (i.e. no other overriding Dirichlet conditions,
             * And find the global unknown number for applying this condition
             */
            index_eq =
                bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn, p, &eqn, &matID_apply, &vd);
            if (index_eq >= 0) {

              /* for weak conditions weight the function by temp */
              weight = wt;
              /* calculate the local variable dof number to apply
               * this boundary condition on. If the unknown doesn't
               * exist, skip calculation.
               */
              ldof_eqn = ei[pg->imtrx]->ln_to_dof[eqn][id];
              if (ldof_eqn >= 0) {

                if (BC_Types[bc_input_id].desc->i_apply == SINGLE_PHASE || pd->e[pg->imtrx][eqn]) {
                  phi_i = bf[eqn]->phi[ldof_eqn];
                  weight *= phi_i;
                } else if (BC_Types[bc_input_id].desc->i_apply == CROSS_PHASE) {
                  /* if the equation type doesn't exist in this material,
                   * then find the intepolation from the adjacent material
                   * using the local node number, id, as the dof number to
                   * look up the basis function.
                   */
                  for (mn = 0; mn < upd->Num_Mat; mn++) {
                    if (pd_glob[mn]->e[pg->imtrx][eqn] &&
                        (eb_in_matrl(BC_Types[bc_input_id].BC_Data_Int[0], mn) ||
                         eb_in_matrl(BC_Types[bc_input_id].BC_Data_Int[1], mn))) {
                      type = pd_glob[mn]->w[pg->imtrx][eqn];
                      if (bfi[type] == NULL)
                        GOMA_EH(GOMA_ERROR, "Illegal cross basis func");
                      /* note that here, we don't have the ln_to_dof array for the adjacent
                         material - for now assume that ldof_eqn = id */
                      phi_i = bfi[type]->phi[id];
                      weight *= phi_i;
                    }
                  }

                } else
                  GOMA_EH(GOMA_ERROR, "Illegal bc phase definition");

                /* for strong conditions weight the function by BIG_PENALTY */
                if (BC_Types[bc_input_id].desc->method == STRONG_INT_EDGE)
                  weight *= BIG_PENALTY;
                if (BC_Types[bc_input_id].desc->method == WEAK_INT_EDGE) {
                  weight *= pd->etm[pg->imtrx][eqn][(LOG2_BOUNDARY)];
                }

                if (BC_Types[bc_input_id].BC_Name == KINEMATIC_EDGE_BC)
                  weight *= BIG_PENALTY;

                /* if doing a weak BC - add into local element contribution also */

                ieqn = upd->ep[pg->imtrx][eqn];
                if (eqn == R_MASS)
                  ieqn = MAX_PROB_VAR + BC_Types[bc_input_id].species_eq;
                lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] += weight * fv->edge_det * func[p];

                /*
                 * add sensitivities into matrix
                 *  - find index of sensitivity in matrix (if variable is not
                 *    defined at this node, loop over all dofs in element)
                 *  - add into matrix
                 */

                if (af->Assemble_Jacobian) {
                  /* if mesh displacement is variable, put in this sensitivity first */
                  /* ... unless we are computing the mass matrix
                   * for LSA.  In that case, we don't include
                   * this first term b/c it doesn't involve any
                   * primary time derivative variables.
                   */
                  if (!af->Assemble_LSA_Mass_Matrix)
                    for (q = 0; q < pd->Num_Dim; q++) {
                      var = MESH_DISPLACEMENT1 + q;
                      if (pd->v[pg->imtrx][var]) {
                        pvar = upd->vp[pg->imtrx][var];
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                              weight * func[p] * fv->dedgedet_dx[q][j];
                        }
                      }
                    }

                  /* now add in sensitivity of BC function */

                  for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                    if (pd->v[pg->imtrx][var] && BC_Types[bc_input_id].desc->sens[var]) {
                      pvar = upd->vp[pg->imtrx][var];
                      if (var != MASS_FRACTION) {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                              weight * fv->edge_det * (d_func[p][var][j] + d_func_ss[p][var][j]);
                        }
                      } else /* variable type is MASS_FRACTION */
                      {
                        for (w = 0; w < pd->Num_Species_Eqn; w++) {
                          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                            lec->J[LEC_J_INDEX(ieqn, MAX_PROB_VAR + w, ldof_eqn, j)] +=
                                weight * fv->edge_det * d_func[p][MAX_VARIABLE_TYPES + w][j];
                          }
                        } /* end of loop over species */
                      }   /* end of if MASS_FRACTION */
                    }     /* end of variable exists and condition is sensitive to it */
                  }       /* end of loop over variable types */
                }         /* end of NEWTON */
              }           /* if (ldof_eqn >= 0) */
            }             /* end of if (Res_BC != NULL) - i.e. apply residual at this node */
          }               /* end of loop over equations that this condition applies to */
        }                 /* end for (i=0; i< num_nodes_on_side; i++) */

      } /*End (if INT) (CAPILLARY and KINEMATIC and VELO_NORMAL and VELO_TANGENT . . .) */
    }   /*(end for ibc) */
  }     /*End for ip = 1,...*/

  return (status);
} /* end of routine apply_integrated_curve_bc */
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

int apply_point_colloc_edge_bc(
    double x[],            /* Solution vector for the current processor    */
    double x_old[],        /* Solution vector at the previous timeon current proc */
    double x_older[],      /* Solution vector for the previous, provious time on current proc*/
    double xdot[],         /* xdot of current solution                     */
    double xdot_old[],     /* xdot of current solution at previous time         */
    double resid_vector[], /* Residual vector for the current processor         */
    const double delta_t,  /* current time step size                            */
    const double theta,    /* parameter to vary time integration:
                         explicit (theta = 1) -- implicit (theta = 0)    */
    const int ielem,       /* element number                                    */
    const int ip_total,    /* total number of gauss points                      */
    const int ielem_type,  /* element type                                      */
    const int num_local_nodes,
    const int ielem_dim,
    const int iconnect_ptr,
    struct elem_edge_bc_struct *elem_edge_bc, /* Pointer to an element side boundary condition
                                                 structure (writable by Gibbs ipin) */
    const int num_total_nodes,
    int local_node_list_fs[], /* dimensioned [MDE]; list to keep track of
                                  nodes at which solid contributions have been
                                  transfered to liquid (fluid-solid boundaries)  */
    const double time_value)

/*************************************************************************
 *
 *  apply_point_colloc_edge_bc():
 *
 *   This function evaluates collocated boundary conditions applied
 *    along curves in 2d and 3d (i.e., intersections of two side sets.)
 *************************************************************************/
{
  int Gibbs = 1;
  int w, i, I, k, j, id, err, var, eqn, icount, p, q, iquad; /* counters */
  int ieqn, pvar, ldof_eqn;
  int status = 0;
  int index_eq, param_dir, matID_apply;
  int bc_input_id;
  double xi[DIM], s, xiq[DIM]; /* Gaussian-quadrature point locations          */
  double x_dot[MAX_PDIM];
  double dsigma_dx[DIM][MDE];
  double func[DIM];
  double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double d_func_ss[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  double theta_r; /* contact angle in radians */
  VARIABLE_DESCRIPTION_STRUCT *vd;
  /* Normals for variable normal contact angle condition etc.                     */
  static double fsnormal[MAX_PDIM]; /* Free surface normal component vector   */
  static double dfsnormal_dx[MAX_PDIM][MAX_PDIM][MDE];
  /* Free surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j        */

  static double ssnormal[MAX_PDIM]; /* Solid surface normal component vector  */
  static double dssnormal_dx[MAX_PDIM][MAX_PDIM][MDE];
  /* Solid surface normal component vector
     derivatives ([i][j][k]) at node k
     (component i wrt displacement j        */

  /***************************************************************************/
  /*   START OF SURFACE LOOPS THAT DON'T REQUIRE INTEGRATION (STRONG SENSE)  */
  /*     (POINTWISE COLLOCATION)                                             */
  /*             LOOP OVER THE SURFACES IN THE CURRENT ELEMENT THAT	     */
  /*        NEED THE APPLICATION OF A STONG SENSE CONDITION                  */
  /*    - INITIALIZATION THAT IS DEPENDENT ON THE IDENTITY OF THE SURFACE    */
  /***************************************************************************/

  /***************************************************************************/
  /*                LOOP OVER THE LOCAL ELEMENT NODE NUMBER                  */
  /*               FOR NODES THAT LIE ON THE CURRENT SURFACE		     */
  /*   - INITIALIZATION THAT IS DEPENDENT ON THE LOCAL ELEMENT NODE NUMBER   */
  /***************************************************************************/

  for (i = 0; i < (int)elem_edge_bc->num_nodes_on_edge; i++) {

    /* Find the local element node number for the current node */
    id = (int)elem_edge_bc->edge_elem_node_id[i];

    /* Find the local node number given the local element node
     * number,  'id'
     * I is the global node number of this node !!
     */
    I = Proc_Elem_Connect[iconnect_ptr + id];

    /*
     * Check to see if the current local node number is owned by the
     * processor - this is done by checking I < num_total_nodes
     */
    if (I < DPI_ptr->num_owned_nodes) {

      /*
       * Load the field variable values at this node point
       */
      param_dir = find_edge_s(0, ielem_type, elem_edge_bc->id_edge, pd->Num_Dim, xi, &s);
      find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");
      err = load_fv();
      GOMA_EH(err, "load_fv");

      /*
       * Load up porous media variables and properties, if needed
       */
      if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_TWO_PHASE ||
          mp->PorousMediaType == POROUS_SATURATED) {
        err = load_porous_properties();
        GOMA_EH(err, "load_porous_properties");
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

      do_LSA_mods(LSA_EDGE);

      /******************************************************************************/
      /*                                                                            */
      /*         LOOP OVER THE INDIVIDUAL SURFACE CONDITION TERMS                   */
      /*		Calculate the functions for each collocation                */
      /*                condition at this point                                     */
      /*		                         			            */
      /******************************************************************************/

      if ((bc_input_id = (int)elem_edge_bc->BC_input_id[0]) != -1) {
        if (BC_Types[bc_input_id].desc->method == COLLOCATE_EDGE) {

          /*
           * Check to see if this BC on this node is applicable (i.e. no other overriding
           *  Dirichlet conditions,
           * And find the global unknown number for applying this condition
           */
          index_eq =
              bc_eqn_index(id, I, bc_input_id, ei[pg->imtrx]->mn, 0, &eqn, &matID_apply, &vd);
          if (index_eq >= 0) {
            memset(func, 0, MAX_PDIM * sizeof(double));
            memset(d_func, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));
            memset(d_func_ss, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

            ldof_eqn = ei[pg->imtrx]->ln_to_dof[eqn][id];

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

            case GD_CONST_BC:
            case GD_LINEAR_BC:
            case GD_PARAB_BC:
            case GD_PARAB_OFFSET_BC:
            case GD_CIRC_BC:
            case GD_TIME_BC:
            case GD_POLYN_BC:

              /*
               * The casting of the 1st and 2nd arguments bothers me.
               * I put them in to silence compiler warnings when
               * more rigorous prototype checking is done (cf, the
               * definition of fgeneralized_dirichlet() in bc_colloc.c
               *
               * This ought to be fixed to work properly.
               */

              err = fgeneralized_dirichlet((double *)(&func), (double *)(d_func),
                                           BC_Types[bc_input_id].BC_Name, bc_input_id, theta,
                                           delta_t);

              GOMA_EH(err, "Illegal entry in Generalized Dirichlet Condition ");
              break;

            case CA_EDGE_CURVE_BC:

              if (elem_edge_bc->shared) {
                GOMA_EH(GOMA_ERROR, "CA_EDGE_CURVE cannot be used with shared edges.");
              }
              /* fall through */
            case CA_EDGE_BC:
              /*
               * need surface vectors
               */
              err = beer_belly();
              GOMA_EH(err, "beer_belly");

              err = load_bf_grad();
              GOMA_EH(err, "load_bf_grad");

              err = load_bf_mesh_derivs();
              GOMA_EH(err, "load_bf_mesh_derivs");

              /* use primary side to find edge vectors */
              edge_determinant_and_vectors(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                           (int)elem_edge_bc->elem_side_bc_1->id_side,
                                           (int)elem_edge_bc->elem_side_bc_1->num_nodes_on_side,
                                           elem_edge_bc->elem_side_bc_1->local_elem_node_id,
                                           /* id_side, num_nodes_on_side, local_elem_node_id,  */
                                           (int)elem_edge_bc->id_edge,
                                           (int)elem_edge_bc->num_nodes_on_edge,
                                           elem_edge_bc->edge_elem_node_id, param_dir);

              /* calculate contact angle boundary condition
               *     use first data entry as contact angle in degrees */
              theta_r = BC_Types[bc_input_id].BC_Data_Float[0] * M_PIE / 180.;

              /* Free surface normal */
              for (p = 0; p < ielem_dim; p++) {
                fsnormal[p] = fv->snormal[p];
                for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                  for (q = 0; q < ielem_dim; q++) {
                    dfsnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
                  }
                }
              }
              /* Wall surface normal */

              if (BC_Types[bc_input_id].BC_Name == CA_EDGE_BC) {
                for (p = 0; p < ielem_dim; p++) {
                  ssnormal[p] = BC_Types[bc_input_id].BC_Data_Float[p + 1];
                  for (q = 0; q < ielem_dim; q++) {
                    for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                      dssnormal_dx[q][p][k] = 0.;
                    }
                  }
                }
              } else if (BC_Types[bc_input_id].BC_Name == CA_EDGE_CURVE_BC) {
                /* recompute surface normal based upon secondary sideset. */

                /* 			surface_determinant_and_normal(ielem,  */
                /* 						       iconnect_ptr,  */
                /* 						       num_local_nodes,  */
                /* 						       ielem_dim - 1,   */
                /* 						       (int)
                 * elem_edge_bc->elem_side_bc_2->id_side, */
                /* 						       (int)
                 * elem_edge_bc->elem_side_bc_2->num_nodes_on_side, */
                /* 						             elem_edge_bc->elem_side_bc_2->local_elem_node_id);
                 */
                edge_determinant_and_vectors(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                             (int)elem_edge_bc->elem_side_bc_2->id_side,
                                             (int)elem_edge_bc->elem_side_bc_2->num_nodes_on_side,
                                             elem_edge_bc->elem_side_bc_2->local_elem_node_id,
                                             /* id_side, num_nodes_on_side, local_elem_node_id,  */
                                             (int)elem_edge_bc->id_edge,
                                             (int)elem_edge_bc->num_nodes_on_edge,
                                             elem_edge_bc->edge_elem_node_id, param_dir);

                for (p = 0; p < ielem_dim; p++) {
                  ssnormal[p] = fv->snormal[p];
                  for (q = 0; q < ielem_dim; q++) {
                    for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                      dssnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
                    }
                  }
                }
              }

              fapply_CA(func, d_func, d_func_ss, fsnormal, dfsnormal_dx, ssnormal, dssnormal_dx,
                        theta_r);

              break;

            case CA_EDGE_OR_FIX_BC:
              /*
               * need surface vectors
               */
              err = beer_belly();
              GOMA_EH(err, "beer_belly");

              err = load_bf_grad();
              GOMA_EH(err, "load_bf_grad");

              err = load_bf_mesh_derivs();
              GOMA_EH(err, "load_bf_mesh_derivs");

              /* use primary side to find edge vectors */
              edge_determinant_and_vectors(ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                                           (int)elem_edge_bc->elem_side_bc_1->id_side,
                                           (int)elem_edge_bc->elem_side_bc_1->num_nodes_on_side,
                                           elem_edge_bc->elem_side_bc_1->local_elem_node_id,
                                           /* id_side, num_nodes_on_side, local_elem_node_id,  */
                                           (int)elem_edge_bc->id_edge,
                                           (int)elem_edge_bc->num_nodes_on_edge,
                                           elem_edge_bc->edge_elem_node_id, param_dir);

              /* calculate contact angle boundary condition
               *     use first data entry as contact angle in degrees */
              theta_r = (double)BC_Types[bc_input_id].u_BC[0] * M_PIE / 180.;

              /* Free surface normal */
              for (p = 0; p < ielem_dim; p++) {
                fsnormal[p] = fv->snormal[p];
                for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                  for (q = 0; q < ielem_dim; q++) {
                    dfsnormal_dx[q][p][k] = fv->dsnormal_dx[q][p][k];
                  }
                }
              }
              /* Wall surface normal */
              /* for now assume that the contact angle condition is applied at a
               * fixed wall - where the wall surface normal is known */
              for (p = 0; p < ielem_dim; p++) {
                ssnormal[p] = (double)BC_Types[bc_input_id].u_BC[p + 1];
                for (q = 0; q < ielem_dim; q++) {
                  for (k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; k++) {
                    dssnormal_dx[q][p][k] = 0.;
                  }
                }
              }
              Gibbs = user_gibbs_criterion(fsnormal, ssnormal, BC_Types[bc_input_id].BC_Data_Int[0],
                                           &(elem_edge_bc->ipin), BC_Types[bc_input_id].u_BC);

              if (Gibbs) {

                fapply_CA(func, d_func, d_func_ss, fsnormal, dfsnormal_dx, ssnormal, dssnormal_dx,
                          theta_r);
              } else {
                if (!af->Assemble_LSA_Mass_Matrix)
                  /* fix contact point where it sits right now */
                  for (p = 0; p < ielem_dim; p++) {
                    k = Index_Solution(I, MESH_DISPLACEMENT1 + p, 0, 0, -2, pg->imtrx);
                    ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[MESH_DISPLACEMENT1 + p][id];
                    lec->R[LEC_R_INDEX(MESH_DISPLACEMENT1 + p, ldof_eqn)] = 0.;
                    d_func[p][MESH_DISPLACEMENT1 + p][id] = 1.;

                    eqn = upd->ep[pg->imtrx][MESH_DISPLACEMENT1 + p];
                    zero_lec_row(lec->J, eqn, ldof_eqn);
                    lec->J[LEC_J_INDEX(eqn, eqn, ldof_eqn, ldof_eqn)] = DIRICHLET_PENALTY;

                    xdot[k] = 0.0;
                    xdot_old[k] = 0.0;
                    x_old[k] = x[k];
                    x_older[k] = x_old[k];
                  }
              }

              break;

            case VELO_NORMAL_EDGE_BC:
            case VELO_TANGENT_EDGE_BC:

              /* These conditions are applied to the Gauss points on
               *  the edge rather than at the nodes. The following
               *  code determines which guass point is nearest
               *  to the current node and computes the constraint
               *  at that gauss point.
               */

              if (elem_edge_bc->num_nodes_on_edge < 3) {
                iquad = (xi[param_dir] == -1.0 ? 1 : (xi[param_dir] == 1.0 ? 0 : -1));
              } else {
                iquad = (xi[param_dir] == -1.0
                             ? 2
                             : (xi[param_dir] == 0.0 ? 1 : (xi[param_dir] == 1.0 ? 0 : -1)));
              }

              GOMA_EH(iquad, "problem finding iquad in apply_point_collocated_edge_bc");

              find_edge_s(iquad, ielem_type, elem_edge_bc->id_edge, pd->Num_Dim, xiq, &s);

              err = load_basis_functions(xiq, bfd);
              GOMA_EH(err, "problem from load_basis_functions");

              err = load_fv();
              GOMA_EH(err, "load_fv");

              /*
               * need surface vectors
               */
              err = beer_belly();
              GOMA_EH(err, "beer_belly");

              err = load_bf_grad();
              GOMA_EH(err, "load_bf_grad");

              err = load_bf_mesh_derivs();
              GOMA_EH(err, "load_bf_mesh_derivs");

              /* use primary side to find edge vectors */
              edge_determinant_and_vectors(
                  ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                  (int)elem_edge_bc->elem_side_bc_1->id_side,
                  (int)elem_edge_bc->elem_side_bc_1->num_nodes_on_side,
                  elem_edge_bc->elem_side_bc_1->local_elem_node_id, (int)elem_edge_bc->id_edge,
                  (int)elem_edge_bc->num_nodes_on_edge, elem_edge_bc->edge_elem_node_id, param_dir);

              do_LSA_mods(LSA_EDGE);

              if (BC_Types[bc_input_id].BC_Name == VELO_NORMAL_EDGE_BC) {
                fvelo_normal_edge_bc(func, d_func, BC_Types[bc_input_id].BC_Data_Float[0], x_dot,
                                     theta, delta_t);
              } else {
                fvelo_tangent_edge_bc(func, d_func, BC_Types[bc_input_id].BC_Data_Float, x_dot,
                                      theta, delta_t);
              }
              break;

              /*
               * ADD more functions here for calculating collocation residuals
               */

            default:
              GOMA_EH(GOMA_ERROR, " Non-existant collocated edge condition");

            } /* end of SWITCH statement */

            /******************************************************************************/
            /*                                                                            */
            /*         ADD the Boundary condition functions into the Residual             */
            /*         vector and Jacobian Matrix                                         */
            /*		                         			            */
            /******************************************************************************/
            for (p = 0; p < BC_Types[bc_input_id].desc->vector; p++) {
              ieqn = upd->ep[pg->imtrx][eqn];
              if (eqn == R_MASS)
                ieqn = MAX_PROB_VAR + BC_Types[bc_input_id].species_eq;
              lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] += BIG_PENALTY * func[p];

              /*
               * add sensitivities into matrix
               *  - find index of sensitivity in matrix (if variable is not defined at this node,
               *      loop over all dofs in element)
               *  - add into matrix
               */

              if (af->Assemble_Jacobian) {
                for (var = 0; var < MAX_VARIABLE_TYPES; var++) {
                  if (pd->v[pg->imtrx][var] && BC_Types[bc_input_id].desc->sens[var]) {
                    if (var != MASS_FRACTION) {
                      /* load sensitivity index at this node point
                       * this routine determines the entry in the jacobian matrix which
                       * corresponds to this BC equation and this unknown
                       */
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        pvar = upd->vp[pg->imtrx][var];
                        ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
                        lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, j)] +=
                            BIG_PENALTY * (d_func[p][var][j] + d_func_ss[p][var][j]);
                      }
                    } else /* variable type is MASS_FRACTION */
                    {
                      for (w = 0; w < pd->Num_Species_Eqn; w++) {
                        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                          // pvar = upd->vp[pg->imtrx][MAX_PROB_VAR + w];
                          ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][id];
                          lec->J[LEC_J_INDEX(ieqn, MAX_PROB_VAR + w, ldof_eqn, j)] +=
                              BIG_PENALTY * d_func[p][MAX_VARIABLE_TYPES + w][j];
                        }
                      } /* end of loop over species */
                    }   /* end of if MASS_FRACTION */
                  }     /* end of variable exists and condition is sensitive to it */
                }       /* end of loop over variable types */
              }         /* end of NEWTON */
            }

          } /* END of if (Res_BC != NULL), i.e. (index_eqn != -1) */
        }   /* END of if COLLOCATE */
            /*****************************************************************************/
      }     /* END for (ibc = 0; (int) elem_side_bc->BC_input_id[ibc] != ...*/
            /*****************************************************************************/
    }       /* END if (I < num_total_nodes) 				      */
            /*****************************************************************************/
  }         /* END for (i = 0; i < (int) elem_side_bc->num_nodes_on_side; i++) */
            /*****************************************************************************/
  return (status);
} /* end of routine apply_collocated_edge_bc */
/*****************************************************************************/
