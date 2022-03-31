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

#ifndef GOMA_BC_INTEG_H
#define GOMA_BC_INTEG_H

#include "bc_dirich.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_INTEG_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_BC_INTEG_C
#define EXTERN extern
#endif

EXTERN int
apply_integrated_bc(double x[],            /* Solution vector for the current processor    */
                    double resid_vector[], /* Residual vector for the current processor    */
                    const double delta_t,  /* current time step size                       */
                    const double theta,    /* parameter (0 to 1) to vary time integration
                                            *  ( implicit - 0 to explicit - 1)             */
                    const PG_DATA *pg_data,

                    const int ielem,      /* element number */
                    const int ielem_type, /* element type */
                    const int num_local_nodes,
                    const int ielem_dim,
                    const int iconnect_ptr,
                    ELEM_SIDE_BC_STRUCT *elem_side_bc, /* Pointer to an element side boundary
                                                        * condition structure */
                    const int num_total_nodes,
                    const int bc_application, /* flag indicating whether to integrate
                                               * strong or weak BC's */
                    const double time_value,
                    SGRID *grid,
                    const Exo_DB *exo);

int apply_nedelec_bc(double x[],            /* Solution vector for the current processor    */
                     double resid_vector[], /* Residual vector for the current processor    */
                     const double delta_t,  /* current time step size                       */
                     const double theta,    /* parameter (0 to 1) to vary time integration
                                             *  ( implicit - 0 to explicit - 1)             */
                     const PG_DATA *pg_data,
                     const int ielem,      /* element number */
                     const int ielem_type, /* element type */
                     const int num_local_nodes,
                     const int ielem_dim,
                     const int iconnect_ptr,
                     ELEM_SIDE_BC_STRUCT *elem_side_bc, /* Pointer to an element side boundary
                                                         * condition structure */
                     const int num_total_nodes,
                     const int bc_application, /* flag indicating whether to integrate
                                                * strong or weak BC's */
                     const double time_value,
                     SGRID *grid,
                     const Exo_DB *exo);
EXTERN void apply_table_wic_bc(double[], /* func                                      */
                               double[][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
                               struct Boundary_Condition *, /* BC_Type                  */
                               double);                     /* time_value */

#ifdef STATIC

/*
 * Prototype declarations of static functions in bc_colloc.c...
 */

#endif

int equation_index_auto_rotate(const ELEM_SIDE_BC_STRUCT *elem_side_bc,
                               int I,
                               int eqn,
                               int p,
                               int ldof_eqn,
                               const BOUNDARY_CONDITION_STRUCT *bc);

#endif /* GOMA_BC_INTEG_H */
