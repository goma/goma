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
 

#ifndef GOMA_MM_SOL_NONLINEAR_H
#define GOMA_MM_SOL_NONLINEAR_H

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_prob_def.h"
#include "rf_io_structs.h"
#include "std.h"

struct Aztec_Linear_Solver_System;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_SOL_NONLINEAR_C
#define EXTERN
#
#endif

#ifndef GOMA_MM_SOL_NONLINEAR_C
#define EXTERN extern
#endif

EXTERN int solve_nonlinear_problem
(struct Aztec_Linear_Solver_System *, /* ams - ptrs to Aztec linear    *
					     * systems                       */
       double [],		/* x - soln vector on this proc              */
       double ,			/* delta_t - time step size                  */
       double ,			/* theta - parameter to vary time            * 
				 * integration from explicit (theta = 1) to  *
				 *   implicit (theta = 0)                    */
       double [],		/* x_old - soln vector @ previous time       */
       double [],		/* x_older - soln vector @ prev2 time        */
       double [],		/* xdot - dxdt predicted for new time        */
       double [],		/* xdot_old - dxdt for previous time         */
       double [],		/* resid_vector                              */
       double [],		/* x_update                                  */
       double [],		/* scale - Scale factor for modified newton  *
				 * resolves                                  */
       int *,			/* converged - whether the Newton iteration  *
				 * has converged (out)                       */
       int *,			/* nprint - counter for time step number     */
       int ,			/* tev - total number elem variables to      *
				 * output to EXODUS II file                  */
       int ,			/* tev_post - extra element post processing  *
				 * results                                   */
       double [],               /* global variable values                    */
       RESULTS_DESCRIPTION_STRUCT *, /* rd - details about post proc vars    */
       int *,			/* gindex                                    */
       int *,			/* gsize                                     */
       double *,		/* gvec                                      */
       double ***,		/* gvec_elem                                 */
       double ,			/* time_value                                */
       Exo_DB *,		/* exo                                       */
       Dpi *,			/* dpi                                       */
       Comm_Ex *,		/* cx                                        */
       int ,			/* nt                                        */
       int *,			/* time_step_reform                          */
       int ,			/* is_steady_state                           */
       double [],		/* x_AC - updating evp quantities            */
       double [],		/* x_AC_dot -                                */
       double ,			/* lambda                                    */
       double *,		/* resid_vector_sens                         */
       double *,		/* x_sens                                    */
       double **,  		/*  x_sens_p - solution sensitivities        */
     void *);                   /* con_ptr pointer                           */

EXTERN double L2_norm		/* mm_sol_nonlinear.c */
(double *,		/* vector */
       int );			/* nloc */

EXTERN double L2_norm_diff	/* mm_sol_nonlinear.c */
(double *,		/* vector 1 */
       double *,		/* vector 2 */
       int );			/* nloc */

EXTERN double L2_norm_r         /* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int );                  /* nloc */

EXTERN double L1_norm		/* mm_sol_nonlinear.c */
(double *,		/* vector */
       int );			/* nloc */

EXTERN double L1_norm_r         /* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int );                  /* nloc */

EXTERN double Loo_norm		/* mm_sol_nonlinear.c */
(double *,		/* vector */
       int ,			/* nloc */
       int *,  			/* num_unk - save the index! */
       char *);		/* dofname_x - dof name for num_unk  */

EXTERN double Loo_norm_r	/* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int ,                    /* nloc */
       int *,                   /* num_unk - save the index! */
       char *);		/* dofname_r - dof name for num_unk  */

EXTERN double L2_norm_1p        /* mm_sol_nonlinear.c */
(double *,                /* vector */
       int );                  /* nloc */

EXTERN double L2_norm_diff_1p   /* mm_sol_nonlinear.c */
(double *,                /* vector 1 */
       double *,                /* vector 2 */
       int );                  /* nloc */

EXTERN double L2_norm_r_1p      /* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int );                  /* nloc */

EXTERN double L1_norm_1p        /* mm_sol_nonlinear.c */
(double *,                /* vector */
       int );                  /* nloc */

EXTERN double L1_norm_r_1p      /* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int );                  /* nloc */

EXTERN double Loo_norm_1p       /* mm_sol_nonlinear.c */
(double *,                /* vector */
       int ,                    /* nloc */
       int *,                   /* num_unk - save the index! */
       char *);                /* dofname_x - dof name for num_unk  */

EXTERN double Loo_norm_r_1p     /* mm_sol_nonlinear.c */
(double *,                /* vector */
       double *,                /* vector scale */
       int ,                    /* nloc */
       int *,                   /* num_unk - save the index! */
       char *);                /* dofname_r - dof name for num_unk  */

EXTERN void print_array		/* mm_sol_nonlinear.c */
(const void *,		/* array - generic pointer */
       const int ,		/* length - of the array */
       const char *,		/* name - used to print name[23] = value */
       const Edt ,		/* datatype - std.h, type_int or type_double */
       const int);		/* procid  */

/*
 * EDW: Moved here from mm_sol_nonlinear.c ---
 * The interface to Bob Benner's frontal solver is here for now. If it
 * were fully road tested you'd expect a nice function prototype declaration
 * in the appropriate include file. To be fair, though, the interface
 * does involve passing a function pointer for a element stiffness matrix
 * assembly routine that could vary substantially depending on the application.
 */

#ifdef HAVE_FRONT
extern int mf_solve_lineqn
(int *,                   /* re_solve */
       double *,                /* rhs */
       int ,                    /* nrhs */
       int *,                   /* nsetbc */
       double *,                /* bcvalue */
       double *,                /* smallpiv */
       double *,                /* singpiv */
       int *,                   /* iautopiv */
       int *,                   /* iscale */
       void (*) ( struct Aztec_Linear_Solver_System *,
                  double [],
                  double [],
                  double [],
                  double [],
                  double [],
                  double [],
                  double [],
                  double *,
                  double *,
                  struct elem_side_bc_struct *[],
                  double *,
                  Exo_DB *,
                  Dpi *,
                  int *,
                  int *,
                  dbl *,
                  dbl *,
                  dbl *,
                  int),       /* matrix_fill(), anyone? see mm_fill.h */
       double *,                /* lhs */
       double *,                /* scaling_max */
       double *,                /* r_scale */
       struct Aztec_Linear_Solver_System *,
       double [],
       double [],
       double [],
       double [],
       double [],
       double [],
       double [],
       double *,
       double *,
       struct elem_side_bc_struct *[],
       double *,
       Exo_DB *,
       Dpi *,
       int *,               /* &num_total_nodes */
       dbl *,               /* &h_elem_avg */
       dbl *                /* &U_norm */
       );
#endif

#endif /* GOMA_MM_SOL_NONLINEAR_H */
