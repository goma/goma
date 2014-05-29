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
 

/* Header file for most of the functions in the suite of sl_eggroll*.c
 * source files.  These files constitute a generalized eigenvalue
 * problem solver with reverse communication.
 *
 * Originally written by Ian Gates, who adapted it from nsh01.c,
 * nsh02.c, nsh52.c, nsh59.c, and nsh71.c, also by Ian Gates.
 *
 * Modification history:
 *   - July 24, 1997, first version.  Altered to make subset of above codes.
 *   - July 31, 1997, completed reverse communication interface.
 *   - Jan 12, 2000, MMH rearranging and cleaning.
 */

#ifndef _SL_EGGROLL_H
#define _SL_EGGROLL_H

#include "rf_io_structs.h"	/* for struct Results_Description */
#include "exo_struct.h"		/* for typedef struct Exo_DB */
#include "dpi.h"		/* for typedef struct Dpi */
#include "sl_eggroll_def.h"	/* for typedef struct EV */
 
extern void gevp_solver_rc	/* sl_eggroll01.c */
PROTO((int ,
       int ,
       int ,
       dbl ,
       int ,
       int *,
       dbl *,
       dbl *,
       dbl *,
       dbl *,
       int *,
       int ,
       int ,
       int *,
       dbl **,
       dbl **,
       int ,
       dbl ,
       int,
       int *,
       int *,
       dbl *,
       dbl *,
       dbl * )); 

extern void gevp_arnoldi_rc	/* sl_eggroll02.c */
PROTO((int ,
       int ,
       int ,
       int ,
       int ,
       int ,
       int ,
       dbl *,
       dbl *,
       dbl ,
       int *,
       dbl *,
       dbl *,
       dbl *,
       dbl **,
       dbl **,
       int *,
       int *,
       dbl *,
       dbl *,
       dbl *));

extern void eigenvv		/* sl_eggroll04.c */
PROTO((int ,
       dbl **,
       EV *));

extern void balanc		/* sl_eggroll04.c */
PROTO((int ,
       dbl **,
       dbl *));

extern void elmhes		/* sl_eggroll04.c */
PROTO((int ,
       dbl **,
       int *));

extern void eltran		/* sl_eggroll04.c */
PROTO((int ,
       dbl **,
       int *,
       EV *));

extern void balbak		/* sl_eggroll04.c */
PROTO((int ,
       dbl *,
       EV *));

extern void eighqr		/* sl_eggroll04.c */
PROTO((int ,
       dbl **,
       EV *));

extern void _heapsort		/* sl_eggroll04.c */
PROTO((int ,
       dbl *,
       dbl *,
       dbl *,
       int *,
       int ));

extern void jeapsort		/* sl_eggroll04.c */
PROTO((int ,
       dbl *,
       dbl *,
       dbl *,
       int *,
       int ));

extern void cdiv		/* sl_eggrollutil.c */
PROTO((dbl ,
       dbl ,
       dbl ,
       dbl ,
       dbl *,
       dbl *));

extern void gevp_order		/* sl_eggroll05.c */
PROTO((int ,			/* nj */
       int ,			/* ev_n */
       dbl *,			/* ev_r */
       dbl *,			/* ev_i */
       dbl *,			/* ev_e */
       dbl *,			/* ev_x */
       dbl **,			/* evect */
       dbl **));		/* schur */

extern void gevp_transformation	/* sl_eggroll03.c */
PROTO((int ,			/* UMF_system_id */
       int ,			/* first */
       int ,			/* fflag */
       int ,			/* format */
       int ,			/* transformation */
       int ,			/* nj */
       int ,			/* nnz */
       int *,			/* ija */
       dbl *,			/* jac */
       dbl *,			/* mas */
       dbl *,			/* mat */
/*      int soln_tech,  */
       dbl *,			/* w */
       dbl *,			/* v */
       dbl ,			/* r_sigma */
       dbl ));			/* i_sigma */

extern void eggrollwrap		/* sl_eggrollwrap.c */
PROTO((int *,			/* Info for eigenvalue extraction */
       dbl *,			/* Info for eigenvalue extraction */
       int *,			/* Column pointer array */
       dbl *,			/* Nonzero array */
       dbl *,			/* Nonzero array - same structure 
				   as jac[] (ija[]) */
       dbl *,			/* Value of the solution vector */
       char *,			/* Name of exoII output file */
       int ,
       dbl ,			/* Time step size */
       dbl ,			/* Variable time integration parameter
				   explicit (theta = 1) to 
				   implicit (theta = 0) */
       dbl *,			/* Value of the old solution vector */
       dbl *,			/* Value of xdot predicted for new 
				   solution */
       dbl *,
       dbl *,
       int *,			/* Whether the Newton has converged */
       int *,			/* Counter for time step number */
       int ,			/* Number of nodal results */
       int ,			/* Number of post processing results */
       struct Results_Description *,
       int *,
       int *,
       dbl *,
       dbl ,
       Exo_DB *,		/* Ptr to finite element mesh db */
       int ,			/* Number of processors used */
       Dpi *));			/* Ptr to distributed processing info */
#endif
