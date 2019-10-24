/************************************************************************ *
 * Goma - Multiphysics finite element software                             *
 * Sandia National Laboratories                                            *
 *                                                                         *
 * Copyright (c) 2019 GOMA                                                 *
 *                                                                         *
 * Authors: Robert Secor and Andrew Cochrane                               *
 *                                                                         *
 * This software is distributed under the GNU General Public License.      *
\************************************************************************/



#ifndef GOMA_MM_FILL_EM_H
#define GOMA_MM_FILL_EM_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_EM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_EM_C
#define EXTERN extern
#endif

#include <complex.h>
#undef I

EXTERN int assemble_emwave/* mm_fill_em.c                           */
      (double ,/* time - present time value         */
       double ,/* tt - parameter to vary time integration
                * from explicit (tt = 1) to
		* implicit (tt = 0)                   */
       double ,/* dt - current time step size               */
       const PG_DATA *,/* dvc_dnode                                 */
       const int ,/*  Light intensity eqn id and var id     */
       const int ,/*  Light intensity eqn id and var id     */
       const int );

EXTERN int apply_em_farfield_direct/* mm_fill_em.c                           */
(double [DIM],     // func
  double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE] , // d_func
  double [DIM] ,   // xi
  const int,
  double *bc_data);// bc_name

EXTERN void complex_cross_vectors(const complex [DIM],
                                  const complex [DIM],
                                  complex [DIM]);

#endif /* GOMA_MM_FILL_EM_H */
