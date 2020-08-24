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
 
#ifndef GOMA_BC_COLLOC_H
#define GOMA_BC_COLLOC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_COLLOC_C
#define EXTERN
#
#endif

#ifndef GOMA_BC_COLLOC_C
#define EXTERN extern
#endif


EXTERN int apply_point_colloc_bc
(double [],		/* resid_vector */
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
       Exo_DB * );		
     
EXTERN void moving_plane
(int ,			/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func */
       dbl *,			/* aa */
       double );		/* time */

EXTERN void fmesh_constraint
(double *,		/* func */
       double [],		/* d_func */
       const int);		/* bc_id */


EXTERN void fplane
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       dbl *);		/* aa - function parameters from data card  */

EXTERN void f_fillet
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       const double *,		/* p - function parameters from data card  */
       const int );		/* number of parameters from bc card  */

EXTERN void f_double_rad
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       const double *,		/* p - function parameters from data card  */
       const int );		/* number of parameters from bc card  */

#ifdef FEATURE_ROLLON_PLEASE
EXTERN void f_feature_rollon
PROTO((const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       const double *,		/* p - function parameters from data card  */
       const int ,		/* number of parameters from bc card  */
       const int ,		/* geometry model id  */
       const double ));		/* time - time at which BC's are evaluated  */
#endif

EXTERN void f_roll_fluid
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - dimensioned [MAX_VARIABLE_TYPES+MAX_CONC] */
       const double *,		/* p - function parameters from data card  */
       const int ,		/* number of parameters from bc card  */
       double * );		/* number of parameters from bc card  */


EXTERN void fvelocity_profile
(const int ,		/* var_flag */
       const int ,		/* ielem_dim */
       const int ,		/* velo_condition */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameters passed in thru input deck */
       const double );		/* time - time at which BC's are evaluated  */

EXTERN void fvelocity_parabola
(const int ,		/* var_flag */
       const int ,		/* ielem_dim */
       const int ,		/* velo_condition */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       const double [],		/* p - parameters passed in thru input deck */
       const double ,		/* time - time at which BC's are evaluated  */
       const int );		/* number of parameters */

EXTERN void f_vestress_parabola
(const int ,		/* var_flag */
       const int ,		/* ielem_dim */
       const int ,		/* velo_condition */
       const int ,		/* mn */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       const double [],		/* p - parameters passed in thru input deck */
       const double ,		/* time - time at which BC's are evaluated  */
       const int );		/* number of parameters */

EXTERN void fspline
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameterize eqn model */
       const double );		/* time - at which bc's are evaluated */

EXTERN void fspline_rs
(const int ,		/* ielem_dim */
       double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double [],		/* p - parameterize eqn model */
       const double );		/* time - at which bc's are evaluated */

EXTERN void fTmelting
(double *,		/* func */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC] */
       double );		/* a1 - function parameter from data card   */

EXTERN int fgeneralized_dirichlet
(double *,		/* func */
       double [],		/* d_func - MAX_VARIABLE_TYPES + MAX_CONC */
       const int ,		/* gd_condition - denoting which condition 
				 * applied */
       const int ,		/* bc_input_id */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0) */
       const double );		/* dt - current time step size          */

EXTERN int load_variable
(double *,		/* x_var - variable value */
       double *,		/* d_x_var - sensitivities of variable value */
       const int ,		/* jvar - variable number */
       const int ,		/* wspec - species number */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0) */
       const double ,		/* dt - current time step size */
       double [] );		/* vector sensitivity vector - SPEED */

extern int bc_eqn_index(int, int, int, int, int, int *, int *,
			VARIABLE_DESCRIPTION_STRUCT **);
EXTERN int
bc_eqn_index_stress(int id,               /* local node number                 */
	            int I,                /* processor node number             */
	            int bc_input_id,      /* boundary condition number         */
	            int curr_matID,       /* Current material ID */
	            int kdir,             /* coordinate index for stress components */
                    int mode,       /* Stress mode number */
	            int *eqn,       /* eqn to which this condition is applied     */
	            int *matID_retn, /* material ID to apply this eqn on           */
	            VARIABLE_DESCRIPTION_STRUCT **vd_retn);

EXTERN int evaluate_time_func
(const double ,		/* time                                      */
       double *,		/* f_time - computed time function           */
       const int );		/* bc_input_id                               */

EXTERN void apply_table_bc
(double *,		/* func                                      */
       double [MAX_VARIABLE_TYPES+MAX_CONC], /* d_func                       */
       struct Boundary_Condition *, /* BC_Type                             */
       double  );                  /* time _value */

EXTERN double interpolate_table
(struct Data_Table *,   /* table               */
       double [],                     /* x            */
       double *,              /* slope                 */
       double []);           /* gradient array         */

EXTERN double table_distance_search
(struct Data_Table *,   /* table               */
       double [],                     /* x            */
       double *,              /* slope                 */
       double []);           /* gradient array         */

EXTERN double interpolate_table_sat
(struct Data_Table *,	/* table                                     */
       double [3]);	      	/* slope                                     */

#endif /* GOMA_BC_COLLOC_H */
