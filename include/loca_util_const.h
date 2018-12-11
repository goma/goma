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
***************************************************************************/
 

/*
#ifndef lint
static char *cvs_utilconsth_id =
  "$Id";
#endif
*/
/*
 -----------------------------------------------------------------------------
   LOCA 1.0: Library of Continuation Algorithms
   Copyright (C) 2001, Sandia National Laboratories
 -----------------------------------------------------------------------------
*/

/*
 *  loca_util_const.h:
 *
 *      This include file contains function prototypes for functions in
 *      rf_util.c.
 */

#ifndef _LOCA_UTIL_CONST_H
#define _LOCA_UTIL_CONST_H

/*****************************************************************************/
/*     EXTERN STATEMENTS FOR GLOBAL FUNCTIONS IN rf_util.c                   */
/*****************************************************************************/

extern void    initialize_util_routines(int n_o, int n_t);
extern void    vec_init(double *u);
extern void    vec_copy(double *dx, double *dy);
extern double  dp(double *x, double *y);
extern double  ip(double *x, double *y);
extern double  ltransnorm(double *x, double *y);
extern double *alloc_vec(void);
extern void    free_vec(double **ptr);
extern double  null_vector_resid(double r_val, double i_val,
                                 double *r_vec, double *i_vec, int mm_flag);
extern void    sort_by_real(int nconv, int ncv, int ldv, double *d, double *v);

#endif
