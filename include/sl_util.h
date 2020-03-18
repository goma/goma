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

/* sl_util.h -- variables of interest during the Aztec soln of matrices
 *              and utility routines for manipulation of matrices 
 *		are declared  here for global usage, but are 
 *              defined and allocated in sl_util.c.
 *
 *
 * Notes:	[0] Look for more information about Aztec in SAND95-1559.
 *
 *
 *
 * Created:  1997/01/30 09:43 MST pasacki@sandia.gov
 *
 * Modified: 1997/02/20 07:36 MST pasacki@sandia.gov
 *
 * $Id: sl_util.h,v 5.1 2007-09-18 18:53:48 prschun Exp $
 */

#ifndef GOMA_SL_UTIL_H
#define GOMA_SL_UTIL_H

#include <mpi.h>

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "std.h"

struct Aztec_Linear_Solver_System;
/*
 *  Definitions of Constants for use in mm_fill_fill routines
 */

#define ADVECT                          0
#define CORRECT                         1
#define PROJECT                         2
#define EXO_READ                        4
#define HUYGENS                         5
#define SURFACES                        6
#define HUYGENS_C                       7
#define SM_OBJECT                       8
#define HUYGENS_MASS_ITER               9
#define SMOLIANSKI_ONLY                 10

#define LS_SURF_POINT                   0
#define LS_SURF_PLANE                   1
#define LS_SURF_CIRCLE                  2
#define LS_SURF_SPHERE                  3
#define LS_SURF_FACET                   4
#define LS_SURF_NS                      5
#define LS_SURF_SS                      6
#define LS_SURF_ISOSURFACE              7
#define LS_SURF_USER                    8
#define LS_SURF_ARC                     9

#define LS_EVOLVE_ADVECT_EXPLICIT       0 /* Subcycled advection equation */
#define LS_EVOLVE_ADVECT_COUPLED        1 /* Fully coupled advection equation */
#define LS_EVOLVE_SLAVE                 2 /* Isosurface slaved to other surface
                                             (current requires not using COUPLED_FILL  */
#define LS_EVOLVE_SEMILAGRANGIAN        3 /* Semi-lagrangian scheme for evolution
                                             (current requires not using COUPLED_FILL  */

#define MAX_NXN_RANK 6 /* To avoid many malloc()'s in solve_NxN_system(). */

#define GRID_SEARCH                     1  /* Search for points of zero level set via nested boxes */
#define SEGMENT_SEARCH                  0  /* Conventional search method on segments from connected nodes */

/*
 * Function prototypes used in sl_util.c
 */

extern void sl_init 
(unsigned int,		/* option flag */
	struct Aztec_Linear_Solver_System *[],
	Exo_DB *,		/* all the mesh information */
	Dpi *,			/* all the distributed processing information */
	Comm_Ex []);		/* after initialization, on a per proc basis */


extern void sl_free
(unsigned int ,		/* option_mask                               */
       struct Aztec_Linear_Solver_System *[]);	/* ams                       */

extern void free_ams
(struct Aztec_Linear_Solver_System *); /* ptr to ONE such system */

extern void set_aztec_options_params 
( int [],			/* options */
	double [] );		/* params */

extern void dump_aztec_status 
( double [] );		/* status - filled by AZ_solve() */

extern void hide_external
(int ,			/* n - order of the original system     (in) */
       int ,			/* m - order of the truncated system    (in) */
       int *,			/* ija - original column pointers       (in) */
       int *,			/* ijas - save ija area                      */
       double *);		/* a - original nonzero matrix values   (in) */

extern void show_external
(int ,			/* n - order of the original system     (in) */
       int ,			/* m - order of the truncated system    (in) */
       int *,			/* ija - original column pointers       (in) */
       int *,			/* ijas - save ija area                      */
       double *);		/* a - original nonzero matrix values   (in) */

extern int cmsr_ma28		/* sl_ma28.c                                 */
(const int ,		/* n - order of matrix system                */
       const int ,		/* nnz - nominal number of nonzeroes in a    */
       double [],		/* a - vector of nonzeroes in a matrix       */
       int [],			/* ija - column numbers nonzero entries      */
       double [],		/* x - space for solution vector             */
       double []);		/* b - right hand side vector                */

extern void row_scaling		/* sl_matrix_util.c                          */
(const int ,		/* N                                         */
       double [],		/* a                                         */
       int [],			/* ija                                       */
       double [],		/* b                                         */
       double []);		/* scale                                     */

/*****************************************************************************/

/*
 * Function prototypes used in sl_lu_fill.c
 */

extern void luf			/* sl_lu_fill.c */
(const int,		/* N - order of matrix */
       const int ,		/* NExt - number of external unknowns */
       const int ,		/* M - number of nonzero matrix entries in a */
       double [],		/* a - nonzero matrix values */
       int [],			/* ija - column ptrs, nonzero matrix values */
       double [],		/* x - RHS on input, solution vector on out */
       const int );		/* factor_flag - what to do */

extern void aztec_stringer      /* Called after AZ_solve(...) */
(int,                     /* status[AZ_why] */
       double,                  /* status[AZ_its] */
       char *);                /* status string  */

extern void solve_NxN_system
(dbl *,			/* A */ 
       dbl *,                   /* b */
       dbl *,                   /* x */
       const int,               /* rank (<= row_size) */
       const int);             /* row_size (A's 2nd dimension is row_size) */

#if defined(ENABLE_AMESOS) && defined(TRILINOS)
/* Use prototype in sl_amesos_interface.h */
#else
extern void amesos_solve_msr
(char *,
       struct Aztec_Linear_Solver_System *,
       double *,
       double *,
       int );
#endif

/*****************************************************************************/
#endif
