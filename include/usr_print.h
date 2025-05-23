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
 *  user_print.h -- prototype declarations for user_print.c
 */

#ifndef GOMA_USER_PRINT_H
#define GOMA_USER_PRINT_H

extern int usr_print(double *, double, double *x, double **, int);
extern void usr_out_hkm(int, double, double, double *);

#endif /* GOMA_USER_PRINT_H */
