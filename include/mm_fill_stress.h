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

#ifndef GOMA_MM_FILL_STRESS_H
#define GOMA_MM_FILL_STRESS_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_stabilization.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "rf_allo.h"
#include "std.h"

struct GomaLinearSolverData;
struct Generalized_Newtonian;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_STRESS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_STRESS_C
#define EXTERN extern
#endif

EXTERN int assemble_stress /* mm_fill_stress.c                          */
    (dbl,                  /* tt - parm to vary time integration from
                            * explicit (tt = 1) to implicit (tt = 0)    */
     dbl,                  /* dt - current time step size               */
     dbl[DIM],             /* h - not scale factors methinks            */
     dbl[DIM][DIM],        /* hh                                        */
     dbl[DIM][MDE],        /* dh_dxnode                                 */
     dbl[DIM],             /* vcent - avg element velocity, which is the
                            * centroid velocity for Q2 and the average of
                            * the vertices for Q1. It comes from the
                            * routine "element_velocity."               */
     dbl[DIM][MDE]);       /* dvc_dnode                                 */

EXTERN int assemble_stress_fortin(dbl,        /* tt - parm to vary time integration from
                                               * explicit (tt = 1) to implicit (tt = 0)    */
                                  dbl,        /* dt - current time step size               */
                                  PG_DATA *); /* dvc_dnode                                 */

EXTERN int assemble_stress_log_conf(dbl tt, dbl dt, PG_DATA *pg_data);

EXTERN int assemble_stress_log_conf_transient(dbl tt, dbl dt, PG_DATA *pg_data);

EXTERN int assemble_stress_level_set(dbl,            /* tt - parm to vary time integration from
                                                      * explicit (tt = 1) to implicit (tt = 0)    */
                                     dbl,            /* dt - current time step size               */
                                     dbl[DIM],       /* h - not scale factors methinks            */
                                     dbl[DIM][DIM],  /* hh                                        */
                                     dbl[DIM][MDE],  /* dh_dxnode                                 */
                                     dbl[DIM],       /* vcent - avg element velocity, which is the
                                                      * centroid velocity for Q2 and the average of
                                                      * the vertices for Q1. It comes from the
                                                      * routine "element_velocity."               */
                                     dbl[DIM][MDE]); /* dvc_dnode                                 */

EXTERN int assemble_gradient /* mm_fill_stress.c                          */
    (dbl,                    /* tt - parm to vary time integration from
                              * explicit (tt = 1) to implicit (tt = 0)    */
     dbl);                   /* dt - current time step size               */

int assemble_rate_of_strain(dbl tt,  /* parameter to vary time integration from
                                      * explicit (tt = 1) to implicit (tt = 0) */
                            dbl dt); /* current time step size */

EXTERN int tensor_dot /* mm_fill_stress.c                          */
    (dbl[DIM][DIM],   /* t1                                        */
     dbl[DIM][DIM],   /* t2                                        */
     dbl[DIM][DIM],   /* t1_dot_t2                                 */
     const int);      /* dim                                       */

EXTERN dbl vec_dot /* mm_fill_stress.c                          */
    (const int,    /* n1                                        */
     dbl *,        /* v1                                        */
     dbl *);       /* v2                                        */

EXTERN void load_modal_pointers     /* mm_fill_stress.c                          */
    (int,                           /* ve_mode - mode number                     */
     dbl,                           /* tt                                        */
     dbl,                           /* dt                                        */
     dbl[DIM][DIM],                 /* s - stress tensor for mode ve_mode        */
     dbl[DIM][DIM],                 /* s_dot - d/dt stress tensor, mode ve_mode  */
     dbl[DIM][DIM][DIM],            /* grad_s - grad stress tensor mode ve_mode  */
     dbl[DIM][DIM][DIM][DIM][MDE]); /* d_grad_s_dm - mesh deriv of grad of
                                     *  stress tensor for mode ve_mode   */

EXTERN int modal_esp_alloc(void);

EXTERN int assemble_surface_stress(Exo_DB *, /* exo - ptr to basic exodus ii mesh info    */
                                   double[], /* x                                         */
                                   struct GomaLinearSolverData *,
                                   dbl[],  /* x_update - last update for x vector       */
                                   double, /* delta_t - current time step size          */
                                   double, /* t_ - parameter to vary time integration
                                            * from explicit (tt = 1) to
                                            * implicit (tt = 0)                         */
                                   int,    /* ielem_type - element type                 */
                                   int,    /* ielem_type_fill - elem type fill function */
                                   int,    /* id_side - id number of current side
                                            * according to EXODUS convention            */
                                   int,    /* neighbor - element neighboring this side  */
                                   int,    /* ielem - current element                   */
                                   int);   /* num_local_nodes - number of nodes per
                                            * element                                   */

EXTERN int neighbor_stress       /* mm_fill_stress.c                          */
    (Exo_DB *,                   /* exo - ptr to basic exodus ii mesh info    */
     dbl[],                      /* x                                         */
     dbl[],                      /* x_update                                  */
     int,                        /* current_elem                              */
     int,                        /* neighbor_elem                             */
     dbl[][MAX_MODES][DIM][DIM], /* stress_neighbor                       */
     dbl[][MAX_MODES][DIM][DIM], /* snv                                   */
     dbl[][MDE],                 /* phi_v                                     */
     int,                        /* num_local_nodes                           */
     int,                        /* nodes_per_side                            */
     int[],                      /* local_elem_node_id                        */
     int,                        /* ielem_type                                */
     int,                        /* ielem_type_fill                           */
     dbl **,                     /* x_n                                       */
     int[MAX_MODES][DIM][DIM]);  /* v_s                                   */

EXTERN int neighbor_stress_table /* mm_fill_stress.c                     */
    (Exo_DB *,                   /* exo - ptr to basic exodus ii mesh info    */
     dbl[],                      /* x                                         */
     dbl[],                      /* x_update                                  */
     int,                        /* current_elem                              */
     dbl[][MAX_MODES][DIM][DIM], /* stress_neighbor                       */
     dbl[][MAX_MODES][DIM][DIM], /* snv                                   */
     dbl[][MDE],                 /* phi_v                                     */
     int,                        /* num_local_nodes                           */
     int,                        /* nodes_per_side                            */
     int[],                      /* local_elem_node_id                        */
     int,                        /* ielem_type                                */
     int,                        /* ielem_type_fill                           */
     dbl **,                     /* x_n                                       */
     int[MAX_MODES][DIM][DIM],   /* v_s                                   */
     int[MAX_MODES][DIM][DIM]);  /* table_ibc                             */

EXTERN void load_neighbor_pointers         /* mm_fill_stress.c                       */
    (Exo_DB *,                             /* exo                                       */
     struct GomaLinearSolverData *,        /* pointer to matrix data */
     int,                                  /* ielem - neighbor element                  */
     int,                                  /* etype - element type                      */
     int,                                  /* mode - stress mode                        */
     int[MAX_MODES][DIM][DIM],             /* R_s - Equation number for mode ve_mode    */
     int[MAX_MODES][DIM][DIM],             /* v_s - Variable number for mode ve_mode    */
     dbl *[DIM][DIM][MDE][DIM][DIM][MDE]); /* J_S_S - Pointer array       */

EXTERN int segregate_stress_update /* mm_fill_stress.c                       */
    (double[]);                    /* x_update                                  */

EXTERN int stress_eqn_pointer(int[MAX_MODES][DIM][DIM]); /* v_s */

EXTERN
dbl numerical_viscosity(dbl[DIM][DIM],                 /* s - total stress */
                        dbl[DIM][DIM],                 /* gamma_cont - continuous shear rate */
                        dbl[MAX_MODES][DIM][DIM][MDE], /* d_mun_dS - derivative of mun wrt S*/
                        dbl[DIM][DIM][MDE]);           /* d_mun_dG - derivative of mun wrt G */

void compute_exp_s(double[DIM][DIM], double[DIM][DIM], double[DIM], double[DIM][DIM]);

void analytical_exp_s(double[DIM][DIM],
                      double[DIM][DIM],
                      double[DIM],
                      double[DIM][DIM],
                      double[DIM][DIM][DIM][DIM]); // d_exp_s_ds

void compute_d_exp_s_ds(dbl[DIM][DIM],            // s - stress
                        dbl[DIM][DIM],            // exp_s
                        dbl[DIM][DIM][DIM][DIM]); // d_exp_s_ds

void compute_saramito_model_terms(
    dbl *,                        // Saramito coefficient (S)
    SARAMITO_DEPENDENCE_STRUCT *, // struct for sCoeff sensitvities
    const dbl[DIM][DIM],          // stress
    const struct Generalized_Newtonian *,
    const int); // bounds S to [0,1] if TRUE. Only use this for postprocessing!

int assemble_stress_sqrt_conf(dbl tt, /* parameter to vary time integration from
                                       * explicit (tt = 1) to implicit (tt = 0) */
                              dbl dt, /* current time step size */
                              PG_DATA *pg_data);

int assemble_stress_conf(dbl tt, /* parameter to vary time integration from
                                  * explicit (tt = 1) to implicit (tt = 0) */
                         dbl dt, /* current time step size */
                         PG_DATA *pg_data);
int sqrt_conf_source(int mode,
                     dbl b[DIM][DIM],
                     dbl source_term[DIM][DIM],
                     dbl d_source_term_db[DIM][DIM][DIM][DIM]);
void compute_a_dot_b(dbl b[DIM][DIM],
                     dbl G[DIM][DIM],
                     dbl a_dot_b[DIM][DIM],
                     dbl d_a_dot_b_db[DIM][DIM][DIM][DIM],
                     dbl d_a_dot_b_dG[DIM][DIM][DIM][DIM]);
#endif /* GOMA_MM_FILL_STRESS_H */
