/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
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

#ifndef GOMA_SL_AUXUTIL_H
#define GOMA_SL_AUXUTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_SL_AUXUTIL_C
#define EXTERN
#
#endif

#ifndef GOMA_SL_AUXUTIL_C
#define EXTERN extern
#endif

#define USE_UPPER_CASE

#if defined(_AIX)
#define USE_LOWER_CASE
#undef USE_UPPER_CASE
#endif

#ifdef USE_LOWER_CASE
#define MV_CSR mv_csr
#define MV_MSR mv_msr
#endif

#ifdef USE_UPPER_CASE
#define MV_CSR MV_CSR
#define MV_MSR MV_MSR
#endif

EXTERN double **Dmatrix_birth(const int,  /* n */
                              const int); /* m */

EXTERN int **Imatrix_birth(const int,  /* n */
                           const int); /* m */

EXTERN double *Dvector_birth(const int); /* n */

EXTERN int *Ivector_birth(const int); /* n */

/*
 * ... the life of these is short, a tenuous span twixt birth & death ...
 */

EXTERN void Dmatrix_death(double **,  /* a */
                          const int,  /* n */
                          const int); /* m */

EXTERN void Imatrix_death(int **,     /* a */
                          const int,  /* n */
                          const int); /* m */

EXTERN void Dvector_death(double *,   /* a */
                          const int); /* n */

EXTERN void Ivector_death(int *,      /* a */
                          const int); /* n */

/*

        VECTOR MANIPULATION

*/

extern double nnorm(int, double *);

extern double dot_product(int, double *, double *);

extern void vcopy(int, double *, double, double *);

extern void vzero(int, double *);

extern void vinit(int, double *, double);

extern void v2sum(int, double *, double, double *, double, double *);

extern void v3sum(int, double *, double, double *, double, double *, double, double *);

extern void v1add(int, double *, double, double *);

extern void v2add(int, double *, double, double *, double, double *);

extern void v3add(int, double *, double, double *, double, double *, double, double *);

extern void vchange_sign(int, double *);

extern void v2product(int, double *, double, double *);

extern void vc_product(int, double *, double *, double *);

extern void vc_quotient(int, double *, double *, double *);

extern double vc_max(int, double *);

extern double vc_min(int, double *);

extern void vsproduct(int, double *, double);

extern void MV_CSR(int *, int *, int *, double *, double *, double *);

extern void MV_MSR(int *, int *, double *, double *, double *);

#endif
