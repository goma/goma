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
 *$Id: rf_solver_const.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */


#ifndef _RF_SOLVER_CONST_H
#define _RF_SOLVER_CONST_H

#ifndef solaris
#ifdef __sun
#ifdef __SVR4
#define solaris
#endif
#endif
#endif

#include "rf_io_structs.h"	/* To know about Results_Description */

#include "rf_bc_const.h"	/* To know about elem_side_bc_struct */

#include "sl_util_structs.h"

#include "exo_struct.h"		/* defn of Exo_DB */
#include "dpi.h"		/* defn of Dpi */
#include "dp_types.h"		/* defn of Comm_Ex */

/*
 * Kinds of solvers available...
 */

#define	SPARSE13a	0	/* Kundert's direct factorization package */
#define	AZTEC		1	/* Shadid's "distributed krysolve 2.0" */
#define MA28		2	/* old Harwell solver; MA32 is better costs $ 
				                       MA42 even better */
#define PIM20		3	/* Brazilian's iterative package */
#define PSP3		4	/* Saad's Minnesota package */
/*IGBRK*/
#define	UMFPACK2	5	/* Davis' multifrontal factorization package 
				   Nearly the same as MA42 */
#define	UMFPACK2F	6	/* Davis' multifrontal factorization package 
				   Force full analysis/factorization every time */
#define	FRONT   	7	/* Hood's Frontal Solver */
#define AMESOS          8       /* Heroux & Co. parallel direct solver package */
#define AZTECOO         9

/*
 * FORTRAN BLAS functions. Inside C, use "DCOPY" and the preprocessor to
 * make it look like the FORTRAN name for this routine.
 */

#define APPEND_UNDERSCORE	/* Default behavior... */

#if defined(_AIX) || defined(hpux)
#ifdef APPEND_UNDERSCORE
#undef APPEND_UNDERSCORE
#endif
#endif

#ifdef solaris
#ifndef APPEND_UNDERSCORE
#define APPEND_UNDERSCORE
#endif
#endif

#endif /* _RF_SOLVER_CONST_H */
