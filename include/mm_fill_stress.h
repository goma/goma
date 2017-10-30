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
 

#ifndef _MM_FILL_STRESS_H
#define _MM_FILL_STRESS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_STRESS_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_STRESS_C
#define EXTERN extern
#endif

EXTERN int assemble_stress	/* mm_fill_stress.c                          */
PROTO((dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current time step size               */
       dbl [DIM],		/* h - not scale factors methinks            */
       dbl [DIM][DIM],		/* hh                                        */
       dbl [DIM][MDE],		/* dh_dxnode                                 */
       dbl [DIM],		/* vcent - avg element velocity, which is the
				 * centroid velocity for Q2 and the average of
				 * the vertices for Q1. It comes from the 
				 * routine "element_velocity."               */
       dbl [DIM][MDE]));	/* dvc_dnode                                 */

EXTERN int assemble_stress_fortin
PROTO((dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current time step size               */
       dbl [DIM],		/* h - not scale factors methinks            */
       dbl [DIM][DIM],		/* hh                                        */
       dbl [DIM][MDE],		/* dh_dxnode                                 */
       dbl [DIM],		/* vcent - avg element velocity, which is the
				 * centroid velocity for Q2 and the average of
				 * the vertices for Q1. It comes from the 
				 * routine "element_velocity."               */
       dbl [DIM][MDE]));	/* dvc_dnode                                 */

EXTERN int assemble_stress_log_conf
PROTO((dbl ,                    /* tt - parm to vary time integration from 
                                 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,                    /* dt - current time step size               */
       dbl [DIM],               /* h - not scale factors methinks            */
       dbl [DIM][DIM],          /* hh                                        */
       dbl [DIM][MDE],          /* dh_dxnode                                 */
       dbl [DIM],               /* vcent - avg element velocity, which is the
                                 * centroid velocity for Q2 and the average of
                                 * the vertices for Q1. It comes from the 
                                 * routine "element_velocity."               */
       dbl [DIM][MDE]));        /* dvc_dnode    */


EXTERN int assemble_stress_level_set
PROTO((dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current time step size               */
       dbl [DIM],		/* h - not scale factors methinks            */
       dbl [DIM][DIM],		/* hh                                        */
       dbl [DIM][MDE],		/* dh_dxnode                                 */
       dbl [DIM],		/* vcent - avg element velocity, which is the
				 * centroid velocity for Q2 and the average of
				 * the vertices for Q1. It comes from the 
				 * routine "element_velocity."               */
       dbl [DIM][MDE]));	/* dvc_dnode                                 */


EXTERN int assemble_gradient	/* mm_fill_stress.c                          */
PROTO((dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ));			/* dt - current time step size               */

EXTERN int tensor_dot		/* mm_fill_stress.c                          */
PROTO((dbl [DIM][DIM],		/* t1                                        */
       dbl [DIM][DIM],		/* t2                                        */
       dbl [DIM][DIM],		/* t1_dot_t2                                 */
       const int ));		/* dim                                       */

EXTERN dbl vec_dot		/* mm_fill_stress.c                          */
PROTO((const int,		/* n1                                        */
       dbl *,			/* v1                                        */
       dbl *));			/* v2                                        */

EXTERN void load_modal_pointers	/* mm_fill_stress.c                          */
PROTO((int ,			/* ve_mode - mode number                     */
       dbl ,			/* tt                                        */
       dbl ,			/* dt                                        */
       dbl [DIM][DIM],		/* s - stress tensor for mode ve_mode        */
       dbl [DIM][DIM],		/* s_dot - d/dt stress tensor, mode ve_mode  */
       dbl [DIM][DIM][DIM],	/* grad_s - grad stress tensor mode ve_mode  */
       dbl [DIM][DIM][DIM][DIM][MDE])); /* d_grad_s_dm - mesh deriv of grad of
					 *  stress tensor for mode ve_mode   */

EXTERN int modal_esp_alloc
PROTO((void));

EXTERN int assemble_surface_stress
PROTO((Exo_DB *,		/* exo - ptr to basic exodus ii mesh info    */
       double [],		/* x                                         */
       struct Aztec_Linear_Solver_System *,
       dbl [],			/* x_update - last update for x vector       */
       double ,			/* delta_t - current time step size          */
       double ,			/* t_ - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       int ,			/* ielem_type - element type                 */
       int ,			/* ielem_type_fill - elem type fill function */
       int ,			/* id_side - id number of current side 
				 * according to EXODUS convention            */
       int ,			/* neighbor - element neighboring this side  */
       int ,			/* ielem - current element                   */
       int ));			/* num_local_nodes - number of nodes per 
				 * element                                   */

EXTERN int neighbor_stress	/* mm_fill_stress.c                          */
PROTO((Exo_DB *,		/* exo - ptr to basic exodus ii mesh info    */
       dbl [],			/* x                                         */
       dbl [],			/* x_update                                  */
       int ,			/* current_elem                              */
       int ,			/* neighbor_elem                             */
       dbl [][MAX_MODES][DIM][DIM], /* stress_neighbor                       */
       dbl [][MAX_MODES][DIM][DIM], /* snv                                   */
       dbl [][MDE],		/* phi_v                                     */
       int ,			/* num_local_nodes                           */
       int ,			/* nodes_per_side                            */
       int [],			/* local_elem_node_id                        */
       int ,			/* ielem_type                                */
       int ,			/* ielem_type_fill                           */
       dbl **,			/* x_n                                       */
       int [MAX_MODES][DIM][DIM])); /* v_s                                   */

EXTERN int neighbor_stress_table     /* mm_fill_stress.c                     */
PROTO((Exo_DB *,		/* exo - ptr to basic exodus ii mesh info    */
       dbl [],			/* x                                         */
       dbl [],			/* x_update                                  */
       int ,			/* current_elem                              */
       dbl [][MAX_MODES][DIM][DIM], /* stress_neighbor                       */
       dbl [][MAX_MODES][DIM][DIM], /* snv                                   */
       dbl [][MDE],		/* phi_v                                     */
       int ,			/* num_local_nodes                           */
       int ,			/* nodes_per_side                            */
       int [],			/* local_elem_node_id                        */
       int ,			/* ielem_type                                */
       int ,			/* ielem_type_fill                           */
       dbl **,			/* x_n                                       */
       int [MAX_MODES][DIM][DIM],   /* v_s                                   */
       int [MAX_MODES][DIM][DIM])); /* table_ibc                             */


EXTERN void load_neighbor_pointers /* mm_fill_stress.c                       */
PROTO((Exo_DB *,		/* exo                                       */
       struct Aztec_Linear_Solver_System *,  /* pointer to matrix data */
       int ,			/* ielem - neighbor element                  */
       int ,			/* etype - element type                      */
       int ,			/* mode - stress mode                        */
       int [MAX_MODES][DIM][DIM],		/* R_s - Equation number for mode ve_mode    */
       int [MAX_MODES][DIM][DIM],		/* v_s - Variable number for mode ve_mode    */
       dbl *[DIM][DIM][MDE][DIM][DIM][MDE])); /* J_S_S - Pointer array       */


EXTERN int segregate_stress_update /* mm_fill_stress.c                       */
PROTO(( double [] ));		/* x_update                                  */

EXTERN int stress_eqn_pointer
PROTO((int [MAX_MODES][DIM][DIM])); /* v_s */

EXTERN dbl numerical_viscosity
PROTO((dbl [DIM][DIM],		/* s - total stress */
       dbl [DIM][DIM],		/* gamma_cont - continuous shear rate */
       dbl [MAX_MODES][DIM][DIM][MDE], /* d_mun_dS - derivative of mun wrt S*/ 
       dbl [DIM][DIM][MDE]));	/* d_mun_dG - derivative of mun wrt G */

void
compute_exp_s(double [DIM][DIM],
	      double [DIM][DIM],
              double [DIM],
	      double [DIM][DIM]);

void
compute_d_exp_s_ds(dbl [DIM][DIM],                   //s - stress
                   dbl [DIM][DIM],                   // exp_s
                   dbl [DIM][DIM][DIM][DIM]);        // d_exp_s_ds

#endif /* _MM_FILL_STRESS_H */
