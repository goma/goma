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
 * mm_fill.h -- prototype declarations for mm_fill.c
 */

#ifndef GOMA_MM_FILL_H
#define GOMA_MM_FILL_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_bc.h"
#include "std.h"

struct GomaLinearSolverData;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_C
#define EXTERN extern
#endif

extern int matrix_fill_full(struct GomaLinearSolverData *,
                            double[], /* x - Solution vector                       */
                            double[], /* resid_vector - Residual vector            */
                            double[], /* x_old -  previous last time step          */
                            double[], /* x_older - previous prev time step         */
                            double[], /* xdot - xdot of current solution           */
                            double[], /* xdot_old - xdot_old of current soln       */
                            double[], /* x_update - last update vector             */
                            double *, /* delta_t - current time step size          */
                            double *, /* theta- parameter to vary time integration */
                            struct elem_side_bc_struct *[],
                            double *, /* time_value  */
                            Exo_DB *, /* exo - ptr to EXODUS II finite element db  */
                            Dpi *,    /* dpi - ptr to distributed processing info  */
                            int *,    /* num_total_nodes - Number of nodes that proc owns */
                            dbl *,    /* h_elem_avg - global average element size  for PSPG */
                            dbl *,    /* U_norm - global average velocity for PSPG */
                            dbl *);

EXTERN int matrix_fill(struct GomaLinearSolverData *,
                       double[], /* x - Solution vector                       */
                       double[], /* resid_vector - Residual vector            */
                       double[], /* x_old -  previous last time step          */
                       double[], /* x_older - previous prev time step         */
                       double[], /* xdot - xdot of current solution           */
                       double[], /* xdot_old - xdot_old of current soln       */
                       double[], /* x_update - last update vector             */
                       double *, /* delta_t - current time step size          */
                       double *, /* theta- parameter to vary time integration
                                  * from explicit (theta = 1) to implicit
                                  * (theta = 0)                               */
                       struct elem_side_bc_struct *[], /* first_elem_side_BC_array
                                                        * This is an array of pointers to the
                                                        * first surface integral defined for
                                                        * each element.
                                                        * It has a length equal to the total
                                                        * number of elements defined on
                                                        * current processor                  */
                       double *,                       /* time_value  */
                       Exo_DB *, /* exo - ptr to EXODUS II finite element db  */
                       Dpi *,    /* dpi - ptr to distributed processing info  */
                       int *,    /* ielem - element number                    */
                       int *,    /* num_total_nodes - Number of nodes that each
                                  * processor is responsible for              */
                       dbl *,    /* h_elem_avg - global average element size
                                  * for PSPG                                  */
                       dbl *,    /* U_norm - global average velocity for PSPG */
                       dbl *,    /* estifm - element stiffness Matrix for
                                  * frontal solver                            */
                       int);

EXTERN int matrix_fill_stress(struct GomaLinearSolverData *,
                              double[], /* x - Solution vector                       */
                              double[], /* resid_vector - Residual vector            */
                              double[], /* x_old -  previous last time step          */
                              double[], /* x_older - previous prev time step         */
                              double[], /* xdot - xdot of current solution           */
                              double[], /* xdot_old - xdot_old of current soln       */
                              double[], /* x_update - last update vector             */
                              double *, /* delta_t - current time step size          */
                              double *, /* theta- parameter to vary time integration
                                         * from explicit (theta = 1) to implicit
                                         * (theta = 0)                               */
                              struct elem_side_bc_struct *[], /* first_elem_side_BC_array
                                                               * This is an array of pointers to the
                                                               * first surface integral defined for
                                                               * each element.
                                                               * It has a length equal to the total
                                                               * number of elements defined on
                                                               * current processor */
                              double *,                       /* time_value  */
                              Exo_DB *, /* exo - ptr to EXODUS II finite element db  */
                              Dpi *,    /* dpi - ptr to distributed processing info  */
                              int *,    /* ielem - element number                    */
                              int *,    /* num_total_nodes - Number of nodes that each
                                         * processor is responsible for              */
                              dbl *,    /* h_elem_avg - global average element size
                                         * for PSPG                                  */
                              dbl *,    /* U_norm - global average velocity for PSPG */
                              dbl *,    /* estifm - element stiffness Matrix for
                                         * frontal solver                            */
                              int);     /* zeroCA */

EXTERN int checkfinite(const char *file,
                       const int line,       /* line                                      */
                       const char *message); /* message                                   */

EXTERN int load_pf_constraint(double pf_constraint, double d_pf_lm[][MDE], double d_lm_pf[][MDE]);

#if defined(CHECK_FINITE) || defined(DEBUG_NAN) || defined(DEBUG_INF)
#define CHECKFINITE(MESSAGE) checkfinite(__FILE__, __LINE__, MESSAGE)
#else
#define CHECKFINITE(MESSAGE) 0
#endif

#endif /* GOMA_MM_FILL_H */
