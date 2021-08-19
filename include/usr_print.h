/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2021 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* 
 *  user_print.h -- prototype declarations for user_print.c
 */

#ifndef GOMA_USER_PRINT_H
#define GOMA_USER_PRINT_H

extern int usr_print(double *, double, double *x, double **, int);
extern void usr_out_hkm(int, double, double, double *);

#endif /* GOMA_USER_PRINT_H */
