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
 
#ifndef _BC_CONTACT_H
#define _BC_CONTACT_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _BC_CONTACT_C
#define EXTERN
#
#endif

#ifndef _BC_CONTACT_C
#define EXTERN extern
#endif

EXTERN int apply_contact_bc
PROTO((double [],               /* Solution vector for the current processor */
       double [],		/* resid_vector */
       const double ,		/* delta_t - current time step size */
       const double ,		/* theta - parameter to vary time integration:
                                 * explicit (theta = 1) --
				 * implicit (theta = 0) */
       const double,            /* global average element size */
       const double [],  	/* average element size */
       const double,	        /* average element viscosity */
       const double,	        /* global velocity norm */
       struct elem_side_bc_struct *[],
        		        /* An array of pointers to the first surface * integral defined for each element.
        		         * It's length is equal to the total number of
        		         * elements defined on the current processor */
       const int,               /* element number */
       const int,               /* element type */
       const int,               /* num_local_nodes */
       const int,               /* ielem_dim */
       const int,               /* iconnect_ptr */
       struct elem_side_bc_struct *,
        		        /* Pointer to an element side boundary condition
        		         * structure */
       const int,
       const int,               /* flag indicating whether to integrate
	     			 * strong or weak BC's */
       const int,               /* Base augmenting condition number */
       double *,                /* gAC array */
       double **,               /* bAC array */
       double **,               /* cAC array */
       double **,               /* dAC array */
       const double,            /* time_value */
       const Exo_DB * ));

EXTERN void find_segment_s_wt
PROTO(( const int,           /* current GQ index */
        const int,           /* number of quadrature points */
        double *,            /* Gaussian-quadrature point (s) (returned) */
        double * ));         /* Gaussian-quadrature wt (returned) */

EXTERN void find_subsurf_st_wt
PROTO(( const int,           /* current GQ index */
        const int,           /* number of quadrature points */
        const int,           /* element type */
        const int,           /* side number */
        const int,           /* number of dimensions */
        double *,            /* element coordinates (returned) */
        double *));           /* GQ weight (returned) */

EXTERN double dof_distance
PROTO(( int,		/* variable type    */
        int ));         /* dof number       */
	
EXTERN double lnn_distance
PROTO(( int ));         /* local node number       */


void gnn_distance
PROTO(( const int ,
        const double [],
        const double [],
        const double [],
        double * ,
        double *,
        double *  ));

EXTERN int apply_embedded_bc
PROTO((int,			/* ielem - element number */
       double [],		/* Solution vector for the current processor    */
       double ,		        /* delta_t - current time step size */
       double ,		        /* theta - parameter to vary time integration:
                                 * explicit (theta = 1) --
				 * implicit (theta = 0) */
       double,		        /* time_value */
       const PG_DATA *,
       int,                     /* Base augmenting condition number */
       double *,                /* gAC array */
       double **,               /* bAC array */
       double **,               /* cAC array */
       Exo_DB * ));

EXTERN void setup_shop_at_point
PROTO(( int,                    /* ielem - element number */
        double *,               /* xi */
        const Exo_DB * ));

EXTERN double fv_at_point
PROTO(( double *,
	int ));

EXTERN double fv_old_at_point
PROTO(( double *xi,
        int var ));


EXTERN double fv_and_phi_at_point
PROTO(( double *,
	int,
        double *,
        int ));

EXTERN int apply_embedded_colloc_bc
PROTO(( int,			/* ielem - element number */
        double [],		/* Solution vector for the current processor    */
        double ,		/* delta_t - current time step size */
        double ,		/* theta - parameter to vary time integration:
                                 * explicit (theta = 1) --
				 * implicit (theta = 0) */
        double,		        /* time_value */
        Exo_DB *,
        Dpi * ));
 
EXTERN void contact_fn_dot_T
PROTO((double [DIM],                                     /* func       */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       double,                                  /*  dt       */
       double [DIM],                            /*  lagrange_mult       */
       int,                                     /* cross */
       int [DIM],                               /* dof_l */
       double [DIM][MDE]));                     /* phi_l */

EXTERN void Lagrange_mult_equation
PROTO((double [DIM],                                     /* func       */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       double,                                  /*  tt       */
       double,                                  /*  dt       */
       double,                                  /*  h_elem_avg       */
       double [DIM],                            /* lm_fluid */
       double [DIM],                            /* x_dot */
       int,                                     /* cross */
       int [DIM],                               /* dof_u */
       double [DIM][MDE],                       /* phi_u */
       double [DIM]));                          /*  fluid_velocity       */

EXTERN void solid_kinematic_bc
PROTO((double [DIM],                                     /* func       */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       double,                                  /*  dt       */
       double,                                  /*  tt       */
       double [DIM],                            /*  fluid_velocity       */
       int,                                     /* cross */
       int [DIM],                               /* dof_u */
       double [DIM][MDE],                       /* phi_u */
       double [DIM]));                          /*  x_dot   (input)      */

EXTERN int jump_down_to_fluid
PROTO((const Exo_DB *,                          /* Ptr to Exodus database */
       int,                                     /* bc_input_id */
       double [],                               /* Solution vector */
       double []));                             /* xi in fluid element */

EXTERN int find_2d_side_id
PROTO((double,
       double));

EXTERN int find_3d_side_id
PROTO((double,
       double,
       double));

EXTERN int first_overlap_ac
PROTO((int,                                     /* solid element number */
       int));                                   /* side ID number */

EXTERN int lookup_active_dof
PROTO((int,					/* Current variable index */
       int,					/* Current DOF */
       int));					/* Proc_Elem_Connect */

EXTERN int assemble_qvapor_source
PROTO (( const double [] ));

EXTERN int assemble_qrad_source
PROTO ((double,
        double,
        double,
        double ));

EXTERN int assemble_embedded_bc 
PROTO(( int ,      /* element number */
	double [],     /* Solution vector for the current processor    */
	double ,      /* current time step size                       */
	double ,
	double ,
	int ,        /* Flag indicating calling function */
	double *,    /* Augmenting condition arrays */
	double **,
	double **,
	Exo_DB *,
	double [DIM] ));

#endif /* _BC_CONTACT_H */
