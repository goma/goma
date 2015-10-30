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
 *$Id: mm_shell_util.h,v 5.2 2010-05-21 21:40:35 sarober Exp $
 */

#ifndef _MM_SHELL_UTIL_H
#define _MM_SHELL_UTIL_H


extern int **elem_friends;
extern int *num_elem_friends;
extern int num_shell_blocks;

EXTERN void init_shell_element_blocks
PROTO((const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int is_shell_element_type
PROTO((const int elem_type));   /* Element type index */

EXTERN int is_shell_element
PROTO((const int elem,          /* Element index */
       const Exo_DB *exo));     /* ExodusII database struct pointer */

EXTERN int is_shell_block
PROTO((const int block_id,      /* Element block ID */
       const Exo_DB *exo));     /* ExodusII database struct pointer */

EXTERN int solve_2x2            /* Solves A*x = b */
PROTO((double a[DIM][DIM],
       double b[DIM],
       double x[DIM]));

EXTERN int solve_3x3            /* Solves A*x = b */
PROTO((double a[DIM][DIM],
       double b[DIM],
       double x[DIM]));

EXTERN int find_stu_from_xyz    /* Inverts isoparametrix map to get xi */
PROTO((const int,               /* Element of interest */
       const double[DIM],       /* xyz coords in */
       double[DIM],             /* stu coords out */
       const Exo_DB *exo));	   /* Ptr to Exodus structure */

EXTERN int bulk_side_id_and_stu /* Determines bulk elem side shared by shell */
PROTO((const int,		/* Bulk element number */
       const int,		/* Shell element number */
       const double[DIM],	/* Shell stu coords in */
       double[DIM],		/* Bulk stu coords out */
       const Exo_DB *exo));	   /* Ptr to Exodus structure */

EXTERN int load_neighbor_var_data  /* Creates neighbor assembly structures */
PROTO((int el1,			   /* Local element number */
       int el2,			   /* Neighbor element number */
       int *ndofs,		   /* Neighbor variable DOF's */
       int *dof_map,               /* Shell-bulk local DOF conversion */
       int ndofptr[MAX_VARIABLE_TYPES][MDE],	/* neighbor dof pointer arrays	*/
       const int id_side,          /* id_side for bulk BC's */
       double xi[DIM],		   /* Local element stu coords */
       const Exo_DB *exo));	   /* Ptr to Exodus structure */

EXTERN int find_stu_on_shell    /* Converts bulk stu to shell stu */
PROTO((const int bulk_elem,     /* Bulk element number */
       const int id_side,       /* Bulk element side ID (Exo/Patran) */
       const int shell_elem,    /* Shell element number */
       const int bulk_dim,      /* 2D or 3D bulk element */
       const double s,          /* Bulk xi[0] */
       const double t,          /* Bulk xi[1] */
       const double u,          /* Bulk xi[2] */
       double xi2[DIM],         /* Shell stu coords */
       const Exo_DB *exo));       

EXTERN int find_stu_on_bulk     /* Converts shell stu to bulk stu */
PROTO((const int id_side,       /* Bulk element side ID (Exo/Patran) */
       const int bulk_dim,      /* 2D or 3D bulk element */
       const double s,          /* Shell xi[0] */
       const double t,          /* Shell xi[1] */
       double xi2[DIM]));       /* Bulk stu coords */

EXTERN int shell_normal_div_s   /* Loads shell normal surface divergence */
PROTO((double *div_s_nv,	/* Surface divergence of normal (scalar) */
       double d_div_s_nv_dnv[DIM][MDE],     /* Self-sensitivities */
       double d_div_s_nv_dmesh[DIM][MDE])); /* Mesh sensitivities */

EXTERN void shell_tangents
PROTO((
       double t0[DIM],
       double t1[DIM],
       double dt0_dx[DIM][DIM][MDE],
       double dt1_dx[DIM][DIM][MDE]
     ));

EXTERN void shell_tangents_seeded
PROTO((
       double t0[DIM],
       double t1[DIM],
       double dt0_dnormal[DIM][DIM][MDE],
       double dt1_dnormal[DIM][DIM][MDE]
     ));


EXTERN void shell_stress_tensor
PROTO((
       double TT[DIM][DIM],
       double dTT_dx[DIM][DIM][DIM][MDE],
       double dTT_dnormal[DIM][DIM][DIM][MDE]
     ));

EXTERN void shell_moment_tensor
PROTO((
       double M[DIM][DIM],
       double dM_dx[DIM][DIM][DIM][MDE],
       double dM_dnormal[DIM][DIM][DIM][MDE],
       double dM_curv1[DIM][DIM][MDE],
       double dM_curv2[DIM][DIM][MDE]
     ));

EXTERN void lubrication_shell_initialize
PROTO((
       int *n_dof,           // Degrees of freedom
       int *dof_map,         // Map of DOFs
       int id_side,          // Side ID
       double xi[DIM],       // Local STU coordinates
       const Exo_DB *exo,    // Exodus database
       int use_def           // Use deformed normal anyway
       ));

EXTERN void Inn
PROTO((
       double v[DIM],        // Input vector
       double w[DIM]         // Rotated output vector
       ));

EXTERN void ShellRotate
PROTO((
       double v[DIM],                 // Input vector
       double dv_dmx[DIM][DIM][MDE],  // Input vector mesh sensitivity
       double w[DIM],                 // Rotated output vector
       double dw_dmx[DIM][DIM][MDE],  // Rotated output vector mesh sensitivity
       int ndof                       // Number of DOFs for mesh equations
       ));

EXTERN void calculate_lub_q_v
PROTO((
       const int EQN, 
       double time,
       double dt,
       double xi[DIM],
       const Exo_DB *exo
     ));

EXTERN void calculate_lub_q_v_old
PROTO((
       const int EQN,
       double time_old,
       double dt_old,
       double xi[DIM],
       const Exo_DB *exo
     ));
EXTERN void ShellBF
PROTO((
       int ev,                     // Equation or variable to fetch basis functions
       int ii,                     // Integer for which DOF to fetch
       double *phi_i,
       double grad_phi_i[DIM],
       double gradII_phi_i[DIM],
       double d_gradII_phi_i_dx[DIM][DIM][MDE],
       int ndof,
       int *dof_map         // Map of DOFs
       ));

EXTERN double shell_saturation_pressure_curve
PROTO((
       double P,
       double *dSdP,
       double *dSdP_P
       ));


EXTERN void calculate_lub_q_v_nonnewtonian
PROTO((
       double time,
       double dt
     ));

EXTERN void calculate_lub_q_v_nonnewtonian_sens
PROTO((  double mu,
         double H,
         double veloU[DIM],
         double veloL[DIM],
         double grad_II_P[DIM],
         double phi_j,
         double grad_phi_j[DIM],
         double shear_top_plus,
         double shear_top_minus,
         double shear_bot_plus,
         double shear_bot_minus,
         double eps_top,
         double eps_bot  
     ));

#endif /* _MM_FILL_SHELL_H */
