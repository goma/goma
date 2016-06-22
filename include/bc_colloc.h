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
 
#ifndef _BC_COLLOC_H
#define _BC_COLLOC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _BC_COLLOC_C
#define EXTERN
#
#endif

#ifndef _BC_COLLOC_C
#define EXTERN extern
#endif

#ifdef USE_CGM 
#include "gm_cgm_typedefs.h"
#endif

EXTERN int apply_point_colloc_bc
PROTO((double [],		/* resid_vector */
       const double ,		/* delta_t - current time step size */
       const double ,		/* theta - parameter to vary time integration:
                                 * explicit (theta = 1) -- 
				 * implicit (theta = 0) */
       const int ,		/* ielem - element number */
       const int ,		/* ip_total - total number of gauss points */
       const int ,		/* ielem_type - element type */
       const int ,		/* num_local_nodes */
       const int ,		/* ielem_dim */
       const int ,		/* iconnect_ptr */
       struct elem_side_bc_struct *, /* elem_side_bc - Pointer to an element 
				      * side boundary condition structure */
       const int ,		/* num_total_nodes */
       int [],			/* local_node_list_fs - dimensioned [MDE]; 
				 * list to keep track of nodes at which solid 
				 * contributions have been transfered to 
				 * liquid (fluid-solid boundaries) */
       const double,            /* time value */
       Exo_DB * ));		
     
EXTERN void moving_plane
PROTO((int ,			/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func */
       dbl *,			/* aa */
       double ));		/* time */

EXTERN void fmesh_constraint
PROTO((double *,		/* func */
       double [],		/* d_func */
       const int));		/* bc_id */


EXTERN void fplane
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       dbl *));		/* aa - function parameters from data card  */

EXTERN void f_fillet
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       const double *,		/* p - function parameters from data card  */
       const int ));		/* number of parameters from bc card  */

#ifdef USE_CGM
EXTERN void sm_fplane
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       PlaneHandle *));        /* pHdl - Handle to a VGI Plane object  */
#endif

EXTERN void fvelocity_profile
PROTO((const int ,		/* var_flag */
       const int ,		/* ielem_dim */
       const int ,		/* velo_condition */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameters passed in thru input deck */
       const double ));		/* time - time at which BC's are evaluated  */

EXTERN void fvelocity_parabola
PROTO((const int ,		/* var_flag */
       const int ,		/* ielem_dim */
       const int ,		/* velo_condition */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       const double [],		/* p - parameters passed in thru input deck */
       const double ,		/* time - time at which BC's are evaluated  */
       const int ));		/* number of parameters */

EXTERN void fspline
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameterize eqn model */
       const double ));		/* time - at which bc's are evaluated */

EXTERN void fspline_rs
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameterize eqn model */
       const double ));		/* time - at which bc's are evaluated */

EXTERN void fTmelting
PROTO((double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double ));		/* a1 - function parameter from data card   */

EXTERN int fgeneralized_dirichlet
PROTO((double *,		/* func */
       double [],		/* d_func - MAX_VARIABLE_TYPES + MAX_CONC */
       const int ,		/* gd_condition - denoting which condition 
				 * applied */
       const int ,		/* bc_input_id */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0) */
       const double ));		/* dt - current time step size          */

EXTERN int load_variable
PROTO((double *,		/* x_var - variable value */
       double *,		/* d_x_var - sensitivities of variable value */
       const int ,		/* jvar - variable number */
       const int ,		/* wspec - species number */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0) */
       const double ));		/* dt - current time step size */

extern int bc_eqn_index(int, int, int, int, int, int *, int *,
			VARIABLE_DESCRIPTION_STRUCT **);

EXTERN int evaluate_time_func
PROTO((const double ,		/* time                                      */
       double *,		/* f_time - computed time function           */
       const int ));		/* bc_input_id                               */

EXTERN void apply_table_bc
PROTO((double *,		/* func                                      */
       double [MAX_VARIABLE_TYPES+MAX_CONC], /* d_func                       */
       struct Boundary_Condition *, /* BC_Type                             */
       double  ));                  /* time _value */

EXTERN double interpolate_table
PROTO((struct Data_Table *,   /* table               */
       double [],                     /* x            */
       double *,              /* slope                 */
       double []));           /* gradient array         */

EXTERN double table_distance_search
PROTO((struct Data_Table *,   /* table               */
       double [],                     /* x            */
       double *,              /* slope                 */
       double []));           /* gradient array         */

EXTERN double interpolate_table_sat
PROTO((struct Data_Table *,	/* table                                     */
       double [3]));	      	/* slope                                     */

#endif /* _BC_COLLOC_H */
