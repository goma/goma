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
	sl_auxutil.h

	ROUTINE PROTOTYPES

	BY IAN GATES

	AUGUST 8 1997

	FIRST VERSION FOR GOMA. 

	MODIFIED 

	IDG SEPT 28 1997
*/

#ifndef _SL_AUXUTIL_H
#define _SL_AUXUTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _SL_AUXUTIL_C
#define EXTERN
#
#endif

#ifndef _SL_AUXUTIL_C
#define EXTERN extern
#endif

#define USE_UPPER_CASE

#if defined (_AIX)
#define USE_LOWER_CASE
#undef USE_UPPER_CASE
#endif

#ifdef USE_LOWER_CASE
#define MV_CSR   mv_csr
#define MV_MSR   mv_msr
#endif

#ifdef USE_UPPER_CASE
#define MV_CSR   MV_CSR
#define MV_MSR   MV_MSR
#endif

EXTERN double **Dmatrix_birth
PROTO((const int ,              /* n */
       const int ));            /* m */

EXTERN int **Imatrix_birth
PROTO((const int ,              /* n */
       const int ));            /* m */

EXTERN double *Dvector_birth
PROTO(( const int  ));          /* n */

EXTERN int *Ivector_birth
PROTO(( const int  ));          /* n */

/*
 * ... the life of these is short, a tenuous span twixt birth & death ...
 */

EXTERN void Dmatrix_death
PROTO(( double **,              /* a */
        const int ,             /* n */
        const int ));           /* m */

EXTERN void Imatrix_death
PROTO(( int **,                 /* a */
        const int ,             /* n */
        const int ));           /* m */

EXTERN void Dvector_death
PROTO(( double *,               /* a */
        const int ));           /* n */

EXTERN void Ivector_death
PROTO(( int *,                  /* a */
        const int ));           /* n */

/*
 
        VECTOR MANIPULATION
 
*/

extern double nnorm PROTO((int, double*));

extern double dot_product PROTO((int, double*, double*));

extern void vcopy PROTO((int, double*, double, double*));

extern void vzero PROTO((int, double*));

extern void vinit PROTO((int, double*, double));

extern void v2sum PROTO((int, double*, double, double*, double, double*));

extern void v3sum 
PROTO((int, double*, double, double*, double, double*, double, double*));

extern void v1add PROTO((int, double*, double, double*));

extern void v2add PROTO((int, double*, double, double*, double, double*));

extern void v3add 
PROTO((int, double*, double, double*, double, double*, double, double*));

extern void vchange_sign PROTO((int, double*));

extern void v2product PROTO((int, double*, double, double*));

extern void vc_product PROTO((int, double*, double*, double*));

extern void vc_quotient PROTO((int, double*, double*, double*));

extern double vc_max PROTO((int, double*));

extern double vc_min PROTO((int, double*));

extern void vsproduct PROTO((int, double*, double));

extern void MV_CSR PROTO((int*, int*, int*, double*, double*, double*));

extern void MV_MSR PROTO((int*, int*, double*, double*, double*));

#endif
