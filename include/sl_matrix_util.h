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
 * $Id: sl_matrix_util.h,v 5.1 2007-09-18 18:53:48 prschun Exp $ 
 */

#ifndef _SL_MATRIX_UTIL_H
#define _SL_MATRIX_UTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _SL_MATRIX_UTIL_C
#define EXTERN
#
#endif

#ifndef _SL_MATRIX_UTIL_C
#define EXTERN extern
#endif

EXTERN void canine_chaos
PROTO((int ,			/* local_order                          (in) */
       int ,			/* local_order_plus - order including 
				 * external rows                        (in) */
       int ,			/* local_nnz - number of nonzeroes      (in) */
       int ,			/* local_nnz_plus - including external 
				 * rows                                 (in) */
       int *,			/* global_order - the real order       (out) */
       int *,			/* global_order_plus - overcounting external 
				 * rows                                (out) */
       int *,			/* global_nnz - the strict count       (out) */
       int *));			/* global_nnz_plus - overcounting external 
				 * rows                                (out) */

EXTERN void print_msr_matrix	/* sl_matrix_util.c                          */
PROTO((int ,			/* n - order of matrix system                */
       int *,			/* ija - column pointer list                 */
       double *,     	        /* a - nonzero matrix entries                */
       double *));		/* x - solution vector                       */

EXTERN void print_vbr_matrix	/* sl_matrix_util.c                          */
PROTO(( struct Aztec_Linear_Solver_System *, /* matrix info                  */
        Exo_DB *,		/* ptr to the whole mesh                     */
        Dpi    *,		/* distributed processing info               */
        int [] ));		/* number of unknowns per node               */

EXTERN void row_sum_scaling_scale /* sl_matrix_util.c                        */
PROTO(( struct Aztec_Linear_Solver_System *,
	double [],		/* scaling matrix                            */
	double []));		/* b - rhs                                   */

EXTERN void row_sum_scaling_scale_AC
PROTO(( double **,              /* cAC */
        double **,              /* dAC */
        double *,               /* gAC */
        int  ));                /* nAC */
                          
EXTERN void row_sum_scale_MSR 
PROTO(( int ,
	double [],
	int   [],
	double  [],
	double  [] ));

EXTERN void row_sum_scale_VBR 
PROTO(( int ,
	double *,
	int    *,
	int    *,
	int    *,
	int    *,
	int    *,
	double *,
	double *));

void
row_sum_scale_epetra(struct Aztec_Linear_Solver_System *ams, double *b, double *scale);

EXTERN void matrix_scaling 
PROTO(( struct Aztec_Linear_Solver_System *,
        double *,
        double ,
        double *));

EXTERN void row_sum_scaling
PROTO(( struct Aztec_Linear_Solver_System *,
	double []));

EXTERN void vector_scaling
PROTO((const int ,		/* N */
       double [],		/* b */
       double []));		/* scale */

int
check_compatible_solver();

/*
 * Prototypes from sl_matrix_dump.c
 */
#ifdef MATRIX_DUMP
extern void matrix_dump_msr(struct Aztec_Linear_Solver_System *,
			    Exo_DB *, Dpi *, double *);
extern void matrix_dump_vbr(struct Aztec_Linear_Solver_System *,
			    Exo_DB *, Dpi *, double *);
#endif

#endif /* _SL_MATRIX_UTIL_H */
