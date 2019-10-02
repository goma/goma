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

EXTERN int assemble_emwave/* mm_fill_terms.c                           */
(double ,/* time - present time value         */
       double ,/* tt - parameter to vary time integration
		* from explicit (tt = 1) to 
		* implicit (tt = 0)                   */
       double ,/* dt - current time step size               */
       const PG_DATA *,/* dvc_dnode                                 */
       const int ,/*  Light intensity eqn id and var id     */
       const int ,/*  Light intensity eqn id and var id     */
       const int );

#endif /* GOMA_MM_FILL_EM_H */
