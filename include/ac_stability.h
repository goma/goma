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
 
#ifndef GOMA_AC_STABILITY_H
#define GOMA_AC_STABILITY_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_STABILITY_C
#define EXTERN /* do nothing */
#else
#define EXTERN extern
#endif

#include "std.h"
#include "exo_struct.h"
#include "dpi.h"
#include "ac_update_parameter.h"

struct Aztec_Linear_Solver_System;
struct Results_Description;

/* Settings for Linear_Stability flag */
#define LSA_NONE          0
#define LSA_NORMAL        1
#define LSA_SAVE          2
#define LSA_3D_OF_2D      3
#define LSA_3D_OF_2D_SAVE 4

/* EDW: LSA integration type options */
#define LSA_VOLUME  111
#define LSA_SURFACE 222
#define LSA_EDGE    333

/* ARPACK method options */
#define LSA_DEFAULT 10
#define LSA_SI 20
#define LSA_CAYLEY 30

EXTERN int Linear_Stability;	/* what kind of linear stability, if any (see #defines above) */
EXTERN int LSA_3D_of_2D_pass;	/* pass # when constructing the jacobian/mass matrices */
EXTERN dbl LSA_3D_of_2D_wave_number; /* value of the wave number */
EXTERN dbl *LSA_wave_numbers;	/* the wave numbers (not necessarily integers!) to compute */
EXTERN int LSA_number_wave_numbers; /* length of LSA_wave_numbers */
EXTERN int LSA_current_wave_number; /* index of current LSA_wave_number */

EXTERN int solve_stability_problem /* ac_solve.c */
(struct Aztec_Linear_Solver_System *,           
       double [],               /* x - Value of the solution vector */
       double ,                 /* delta_t - time step size */
       double ,                 /* theta - parm to vary time integration
                                   from explicit (theta = 1) to
                                   implicit (theta = 0) */
       double [],               /* resid_vector */
       double [],               /* x_old - Value of the old solution vector */
       double [],               /* x_older */
       double [],               /* x_dot - dxdt predicted for new soln */
       double [],               /* x_dot_old - dxdt predicted for new soln */
       double [],               /* x_update - dxdt predicted for new soln */
       int *,                   /* converged - boolean for Newton converged */
       int *,                   /* nprint - counter for time step number */
       int ,                    /* tnv - total number nodal variables */
       int ,                    /* tnv_post - same for extra postproc vars */
       struct Results_Description  *, /* rd - for output */
       int *,                   /* gindex */
       int *,                   /* p_gsize */
       double *,                /* gvec */
       double,                  /* time_value */
       Exo_DB *,                /* ptr to finite element mesh database */
       Dpi *);                 /* distributed processing information */

EXTERN int solve_full_stability_problem /* ac_solve.c */
(struct Aztec_Linear_Solver_System *,           
       double [],               /* x - Value of the solution vector */
       double ,                 /* delta_t - time step size */
       double ,                 /* theta - parm to vary time integration
                                   from explicit (theta = 1) to
                                   implicit (theta = 0) */
       double [],               /* resid_vector */
       double [],               /* x_old - Value of the old solution vector */
       double [],               /* x_older */
       double [],               /* x_dot - dxdt predicted for new soln */
       double [],               /* x_dot_old - dxdt predicted for new soln */
       double [],               /* x_update - dxdt predicted for new soln */
       int *,                   /* converged - boolean for Newton converged */
       int *,                   /* nprint - counter for time step number */
       int ,                    /* tnv - total number nodal variables */
       int ,                    /* tnv_post - same for extra postproc vars */
       struct Results_Description  *, /* rd - for output */
       int *,                   /* gindex */
       int *,                   /* p_gsize */
       double *,                /* gvec */
       double,                  /* time_value */
       Exo_DB *,                /* ptr to finite element mesh database */
       Dpi *);                 /* distributed processing information */

EXTERN int solve_3D_of_2D_stability_problem /* ac_solve.c */
(struct Aztec_Linear_Solver_System *,           
       double [],               /* x - Value of the solution vector */
       double ,                 /* delta_t - time step size */
       double ,                 /* theta - parm to vary time integration
                                   from explicit (theta = 1) to
                                   implicit (theta = 0) */
       double [],               /* resid_vector */
       double [],               /* x_old - Value of the old solution vector */
       double [],               /* x_older */
       double [],               /* x_dot - dxdt predicted for new soln */
       double [],               /* x_dot_old - dxdt predicted for new soln */
       double [],               /* x_update - dxdt predicted for new soln */
       int *,                   /* converged - boolean for Newton converged */
       int *,                   /* nprint - counter for time step number */
       int ,                    /* tnv - total number nodal variables */
       int ,                    /* tnv_post - same for extra postproc vars */
       struct Results_Description  *, /* rd - for output */
       int *,                   /* gindex */
       int *,                   /* p_gsize */
       double *,                /* gvec */
       double,                  /* time_value */
       Exo_DB *,                /* ptr to finite element mesh database */
       Dpi *);                 /* distributed processing information */

EXTERN void output_stability_matrices
(double *,
       double *,
       int *,
       int,
       int,
       int);

EXTERN void compare_mass_matrices
(double *,
       double *,
       double *,
       int *,
       int);

extern void eggroll_init         /* replacement for sl_eggrollinit.c */
(int ,                    /* Number of equations */
       int ,                    /* Number of non-zeroes */
       int *,                   /* Info for eigenvalue extraction */
       dbl *);                 /* Info for eigenvalue extraction */

#endif /* GOMA_AC_STABILITY_H */
