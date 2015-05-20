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

#ifndef _RF_SOLVE_SEGREGATED_H
#define _RF_SOLVE_SEGREGATED_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RF_SOLVE_SEGREGATED_C
#define EXTERN
#endif

#ifndef _RF_SOLVE_C
#define EXTERN extern
#endif

/*
 * rf_solve.c prototypes
 */

EXTERN void solve_problem_segregated(Exo_DB *, /* exo - ptr to finite element mesh database */
                                     Dpi *, /* dpi - ptr to distributed processing info */
                                     dbl *);
#endif /* _RF_SOLVE_H */
