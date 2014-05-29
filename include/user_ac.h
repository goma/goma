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
 * $Id: user_ac.h,v 5.2 2009-02-13 20:22:57 hkmoffa Exp $
 */

#ifndef _USER_AC_H
#define _USER_AC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _USER_AC_C
#define EXTERN
#
#endif

#ifndef _USER_AC_C
#define EXTERN extern
#endif

EXTERN void user_aug_cond_residuals
PROTO((int ,			/* nAC                                       */
       double *,		/* x                                         */
       double *,		/* xdot                                      */
       double ,			/* delta_t                                   */
       double ,			/* time_value                                */
       double **,		/* x_sens_p                                  */
       double *,		/* AC                                        */
       int *,			/* have_bAC                                  */
       int *,			/* have_cAC                                  */
       int *,			/* have_dAC                                  */
       double **,		/* bAC                                       */
       double **,		/* cAC                                       */
       double **,		/* dAC                                       */
       Exo_DB *,		/* exo                                       */
       Dpi *,			/* dpi                                       */
       Comm_Ex *));		/* cx                                        */

EXTERN void user_aug_cond_volume_residuals 
PROTO((
       const int iAC, 
       const double * const x, 
       const double * const xdot,
       const double delta_t,
       const double time_value,
       const double * const x_AC,
       double * const AC,
       double ** const cAC,
       double **const dAC,
       const int numProcUnknowns,
       const Exo_DB * const exo,
       const Dpi * const dpi,
       const Comm_Ex * const cx));


#endif /* _USER_AC_H */
