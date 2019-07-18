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
 

/*
 *$Id: mm_fill_shell.h,v 5.9 2010-04-07 22:27:00 prschun Exp $
 */

#ifndef _MM_FILL_SHELL_H
#define _MM_FILL_SHELL_H

extern int InShellElementWithParentElementCoverage;
extern int ShellElementParentElementCoverageForVariable[MAX_VARIABLE_TYPES];


EXTERN int assemble_surface_charge
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo,      /* ExodusII database struct pointer */
       const int eqn));      /* eqn applied to */

EXTERN int assemble_shell_structure
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_shell_tension
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_shell_coordinates
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));        /* ExodusII database struct pointer */

EXTERN int assemble_shell_diffusion
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));        /* ExodusII database struct pointer */

EXTERN int assemble_shell_geometry
PROTO((double time_value,        /* Time */
       double theta,             /* Time stepping parameter */
       double delta_t,           /* Time step size */
       const double wt,          /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));        /* ExodusII database struct pointer */

EXTERN void shell_surface_charge_bc
PROTO((double [DIM],             /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double [MAX_PDIM],  /* x_dot */
       const double,             /* theta or tt*/
       const double,             /* delta_t or dt */
       const int,                /* id_side for bulk BC's */
       const double,             /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo,        /* ExodusII database struct pointer */
       const int ));		 /* strong integrated bc flag  */

EXTERN void surface_electric_field_bc
PROTO((double [MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],  /* local_r */
       double [MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],  /* local_j */
       const int,                /* bulk element matrl index */
       const int *,              /* ei->dof for bulk element */
       const double));           /* Gauss surface point weight */

EXTERN void surface_acoustic_velocity_bc
PROTO((double [MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],  /* local_r */
       double [MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],  /* local_j */
       const int,                /* bulk element matrl index */
       const int *,              /* ei->dof for bulk element */
       const double,             /* Gauss surface point weight */
       const double));           /* time_value */

EXTERN void apply_surface_viscosity
PROTO((double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double ,           /* Surface shear viscosity */
       const double ,           /* surface dilational viscosity */
       const double ,           /* time_start */
       const double ,           /* time_full */
       const double ,           /* time_value */
       struct elem_side_bc_struct *, /* elem_side_bc                         */
       const double,            /* Gauss surface point weight */
       double xi[DIM],           /* local stu coords */
       const Exo_DB *exo,       /* ExodusII database struct pointer */
       const int ));		/* iconnect_ptr */

EXTERN int assemble_shell_angle
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* delta_t or dt */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_shell_surface_rheo_pieces
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* delta_t or dt */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_lubrication
PROTO((const int,                /* Equation type: R_LUBP or R_LUBP2 */
       double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */
EXTERN int assemble_shell_energy
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double xi[DIM],           /* Local stu coords */
       const PG_DATA *pg_data,   /* petrov galerkins stuff */
       const Exo_DB *exo));      /* ExodusII database struct pointer */
EXTERN int assemble_shell_deltah
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */
EXTERN int assemble_lubrication_curvature
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       const PG_DATA *pg_data,   /* petrov galerkins stuff */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */
EXTERN int assemble_lubrication_curvature_2
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       const PG_DATA *pg_data,   /* petrov galerkins stuff */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */
EXTERN int assemble_film
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double xi[DIM],           /* Local stu coordinates */
       const Exo_DB *exo));                


EXTERN int assemble_film_particles
PROTO(( double time,   /* present time value */
	double tt,     /* parameter to vary time integration from 
			 explicit (tt = 1) to implicit (tt = 0)    */
	double dt,     /* current time step size */  
        double xi[DIM],           /* Local stu coordinates */
        const PG_DATA *pg_data,   /* petrov galerkins stuff */
        const Exo_DB *exo));                


EXTERN int assemble_porous_shell_closed
PROTO((double,                   /* theta or tt */
       double,                   /* dt */
       double xi[DIM],
       const Exo_DB *exo
));

EXTERN int assemble_porous_shell_gasn
PROTO((
       double,
       double,
       double xi[DIM],
       const Exo_DB *exo
));
     
EXTERN int assemble_porous_shell_open
PROTO((double,                   /* theta or tt */
       double,                  /* dt */
       double xi[DIM],     
       const Exo_DB *exo
));
  
EXTERN int assemble_porous_shell_open_2
PROTO((double,                   /* theta or tt */
       double,                   /* dt */
       double xi[DIM],
       const Exo_DB *exo
));

EXTERN void shell_diff_kinematic_bc
PROTO((double [DIM],             /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double [MAX_PDIM],  /* x_dot */
       const double,             /* theta or tt*/
       const double,             /* delta_t or dt */
       const int,                /* id_side for bulk BC's */
       const double,             /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo,        /* ExodusII database struct pointer */
       const int ));		 /* strong integrated bc flag  */

EXTERN void shell_lubr_solid_struct_bc
PROTO((double [DIM],             /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double [MAX_PDIM],  /* x_dot */
       const double,             /* theta or tt*/
       const double,             /* delta_t */
       const int,                /* id_side for bulk BC's */
       const double,             /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo,        /* ExodusII database struct pointer */
       const double));           /* BC_data_float[0] */

EXTERN void rep_force_shell_n_dot_f_bc
PROTO((double [DIM],             /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double [MAX_PDIM],  /* x_dot */
       const double,             /* theta or tt*/
       const double,             /* delta_t or dt */
       const int ,                      /* ip         */
       const int ,                      /* ip_total         */
       const int,                /* id_side for bulk BC's */
       const double,             /* Gauss surface point weight */
       double xi[DIM],           /* Local stu coords */
       const Exo_DB *exo,        /* ExodusII database struct pointer */
       const int ));		 /* strong integrated bc flag  */

EXTERN void surface_user_shell_bc
PROTO((double [MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],  /* local_r */
       double [MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],  /* local_j */
       const int,                /* bulk element matrl index */
       const int *,              /* ei->dof for bulk element */
       const double,           /* Gauss surface point weight */
       const double,             /* theta or tt*/
       const double,           /* delta_t or dt */
       const double *));           /* coords */

EXTERN void surface_lubrication_shell_bc
PROTO((double [MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE],  /* local_r */
       double [MAX_PROB_VAR+MAX_CONC][MAX_PROB_VAR+MAX_CONC][MAX_NODES_PER_SIDE][MDE],  /* local_j */
       const int,                /* bulk element matrl index */
       const int *,              /* ei->dof for bulk element */
       const double,           /* Gauss surface point weight */
       const double,             /* theta or tt*/
       const double,           /* delta_t or dt */
       const double *,           /* coords */
       int [MDE],         /* dof_map */
       int [MAX_VARIABLE_TYPES][MDE]));   /* n_dofptr  */

EXTERN void dPdz_function
PROTO(( dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl ,
	dbl *,
	dbl *,
	dbl *,
	dbl *,
	dbl *,
	dbl *));

EXTERN void dPdz_calc
PROTO(( dbl ,
	dbl ,
	dbl *,
	dbl *,
	dbl *,
	dbl *,
	dbl *,
	dbl * ));


EXTERN int assemble_lubrication_power_law
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double [DIM],             /* Local stu coords */
       const Exo_DB * ));        /* ExodusII database struct pointer */  

EXTERN int assemble_shell_tfmp
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double [DIM],             /* Local stu coords */
       PG_DATA *,
       const Exo_DB * ));        /* ExodusII database struct pointer */  

EXTERN int load_lsi_shell_second
PROTO(( const double ));       /* width */


EXTERN int assemble_shell_normal
PROTO((double [DIM],             /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_shell_curvature
PROTO((double [DIM],             /* Local stu coords */
       const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_shell_mesh
PROTO((double,                   /* Time */
       double,                   /* theta or tt*/
       double,                   /* dt */
       double [DIM],             /* Local stu coords */
       const Exo_DB *exo ));     /* ExodusII database struct pointer */


#endif /* _MM_FILL_SHELL_H */
