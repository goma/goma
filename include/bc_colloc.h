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

#ifndef GOMA_BC_COLLOC_H
#define GOMA_BC_COLLOC_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_fem_const.h"
#include "rf_vars_const.h"
#include "std.h"

struct Boundary_Condition;
struct Data_Table;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_COLLOC_C
#define EXTERN
#endif

#ifndef GOMA_BC_COLLOC_C
#define EXTERN extern
#endif

EXTERN int
apply_point_colloc_bc(double[],                     /* resid_vector */
                      const double,                 /* delta_t - current time step size */
                      const double,                 /* theta - parameter to vary time integration:
                                                     * explicit (theta = 1) --
                                                     * implicit (theta = 0) */
                      const int,                    /* ielem - element number */
                      const int,                    /* ip_total - total number of gauss points */
                      const int,                    /* ielem_type - element type */
                      const int,                    /* num_local_nodes */
                      const int,                    /* ielem_dim */
                      const int,                    /* iconnect_ptr */
                      struct elem_side_bc_struct *, /* elem_side_bc - Pointer to an element
                                                     * side boundary condition structure */
                      const int,                    /* num_total_nodes */
                      int[],                        /* local_node_list_fs - dimensioned [MDE];
                                                     * list to keep track of nodes at which solid
                                                     * contributions have been transfered to
                                                     * liquid (fluid-solid boundaries) */
                      const double,                 /* time value */
                      Exo_DB *);

EXTERN void fvelocity_profile(const int,     /* var_flag */
                              const int,     /* ielem_dim */
                              const int,     /* velo_condition */
                              double *,      /* func */
                              double[],      /* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
                              double[],      /* p - parameters passed in thru input deck */
                              const double); /* time - time at which BC's are evaluated  */

EXTERN void fvelocity_parabola(const int,      /* var_flag */
                               const int,      /* ielem_dim */
                               const int,      /* velo_condition */
                               double *,       /* func */
                               double[],       /* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
                               const double[], /* p - parameters passed in thru input deck */
                               const double,   /* time - time at which BC's are evaluated  */
                               const int);     /* number of parameters */

EXTERN void f_vestress_parabola(const int,      /* var_flag */
                                const int,      /* ielem_dim */
                                const int,      /* velo_condition */
                                const int,      /* mn */
                                double *,       /* func */
                                double[],       /* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
                                const double[], /* p - parameters passed in thru input deck */
                                const double,   /* time - time at which BC's are evaluated  */
                                const int);     /* number of parameters */

EXTERN void fTmelting(double *, /* func */
                      double[], /* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
                      double);  /* a1 - function parameter from data card   */

#endif /* GOMA_BC_COLLOC_H */
