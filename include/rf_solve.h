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
 
#ifndef _RF_SOLVE_H
#define _RF_SOLVE_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RF_SOLVE_C
#define EXTERN
#endif

#ifndef _RF_SOLVE_C
#define EXTERN extern
#endif

/*
 * Declarations for top-level structures
 */
extern Comm_Ex *cx;

/*
 * rf_solve.c prototypes
 */

EXTERN void slow_square_dgemm         /* Matrix multiplication C = A X B */
PROTO(( int ,                    /* Whether multiplying by transpose */
        int ,                    /* Dimensions of matrix */
        double [DIM][DIM] ,            /* A */
        double [DIM][DIM] ,            /* B */
        double [DIM][DIM]));           /* C */

EXTERN void initial_guess_stress_to_log_conf
PROTO((double *,                /* x array (solutions from initial guess) */
       int));                   /* num_total_nodes */

EXTERN void solve_problem
PROTO((Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *,  			/* dpi - ptr to distributed processing info */
       dbl *));			/* te_out - return actual end time */

EXTERN int anneal_mesh		/* rf_solve.c */
PROTO(( double [],		/* x - solution vector */
	int ,			/* tev - number elem results */
	int ,			/* tev_post - xtra elem res */
	double *,               /* global variable values */
	struct Results_Description  *,	/* rd - info about results */
	double ,		/* time_value  */
	Exo_DB *,		/* exo - entire mesh desc. */
	Dpi *));		/* dpi - distr proc info */

EXTERN int anneal_mesh_with_external_field   /* rf_solve.c */
PROTO((const Exo_DB * ));                    /* *exo  */

#ifdef LIBRARY_MODE
EXTERN int load_export_vars     /* rf_solve.c */
PROTO((const int,               /* num_nodes */
       dbl [],                  /* x - solution vector */
       dbl *));                /* x_pp - post processing varibale vector */
#endif

/*
 * rf_setup_problem.c prototypes
 */

EXTERN int setup_problem
PROTO((Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *));			/* dpi - ptr to distributed processing info */

EXTERN int free_problem
PROTO((Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *));			/* dpi - ptr to distributed processing info */

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
#endif /* _RF_SOLVE_H */
