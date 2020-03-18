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
 
#ifndef GOMA_RF_SOLVE_H
#define GOMA_RF_SOLVE_H

#include "dp_types.h"
#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "rf_shape.h"
#include "std.h"

struct Results_Description;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RF_SOLVE_C
#define EXTERN
#endif

#ifndef GOMA_RF_SOLVE_C
#define EXTERN extern
#endif

/*
 * Declarations for top-level structures
 */
extern Comm_Ex **cx;

/*
 * rf_solve.c prototypes
 */

EXTERN void slow_square_dgemm         /* Matrix multiplication C = A X B */
( int ,                    /* Whether multiplying by transpose */
        int ,                    /* Dimensions of matrix */
        double [DIM][DIM] ,            /* A */
        double [DIM][DIM] ,            /* B */
        double [DIM][DIM]);           /* C */

EXTERN void initial_guess_stress_to_log_conf
(double *,                /* x array (solutions from initial guess) */
       int);                   /* num_total_nodes */

EXTERN void solve_problem
(Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *,  			/* dpi - ptr to distributed processing info */
       dbl *);			/* te_out - return actual end time */

EXTERN int anneal_mesh		/* rf_solve.c */
( double [],		/* x - solution vector */
	int ,			/* tev - number elem results */
	int ,			/* tev_post - xtra elem res */
	double *,               /* global variable values */
	struct Results_Description  *,	/* rd - info about results */
	double ,		/* time_value  */
	Exo_DB *,		/* exo - entire mesh desc. */
	Dpi *);		/* dpi - distr proc info */

EXTERN int anneal_mesh_with_external_field   /* rf_solve.c */
(const Exo_DB * );                    /* *exo  */

#ifdef LIBRARY_MODE
EXTERN int load_export_vars     /* rf_solve.c */
(const int,               /* num_nodes */
       dbl [],                  /* x - solution vector */
       dbl *);                /* x_pp - post processing varibale vector */
#endif

/*
 * rf_setup_problem.c prototypes
 */

EXTERN int setup_problem
(Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN int free_problem
(Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

void
set_matrix_index_and_global_v(void);

void
predict_solution(int N, double delta_t, double delta_t_old,
		 double delta_t_older, double theta_arg, double x[],
		 double x_old[], double x_older[], double x_oldest[],
		 double xdot[], double xdot_old[], double xdot_older[]);

extern int coordinate_discontinuous_variables(Exo_DB *,	Dpi *);
extern void determine_dvi_index(void);
extern void reconcile_bc_to_matrl(void);
extern void bc_set_index(Exo_DB *);
extern void bc_internal_boundary(Exo_DB *);
extern void bc_matrl_index(Exo_DB *);
extern int find_MaxMatrlPerNode(void);
extern void setup_external_nodal_matrls(Exo_DB *, Dpi *, Comm_Ex *);
extern int create_periodic_acs(Exo_DB *);

/*
 * main.c function prototypes
 */
extern void print_code_version(void);
extern void echo_command_line( int, char *[], char *);
#endif /* GOMA_RF_SOLVE_H */
