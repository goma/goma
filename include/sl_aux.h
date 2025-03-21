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
 * Header for sl_aux.c of handy LU decompositon and back substitution routines.
 *
 * Randy Lober 9/98.
 */

#ifndef GOMA_SL_AUX_H
#define GOMA_SL_AUX_H

/*
 * Function prototypes used in sl_aux.c
 */

extern int
lu_decomp_backsub_driver(double **, /* coeff_matrix[0..sys_size-1][0..sys_size-1] */
                         double *,  /* rhs_vector  [0..sys_size-1] */
                         int *,     /* indx        [0..sys_size-1] (working space vector) */
                         int,       /* sys_size (the dimension size of coeff_matrix, rhs_vector) */
                         int);      /* lu_dcmp_flag (1 = perform
                                       lu decomposition, 0 simply back_sub) */

extern int lu_decomp(double **, /* coeff_matrix[0..sys_size-1][0..sys_size-1]*/
                     const int, /* sys_size (the dimension size of the       *
                                 * coeff_matrix, rhs_vector)                 */
                     int *,     /* indx[0..n-1] (vector recording row        *
                                 * permutation)                              */
                     double *); /* d[1] (row exchange counter)               */

extern void lu_backsub(double **, /* coeff_matrix[0..sys_size-1][0..sys_size-1] */
                       int,       /* sys_size (the dimension size of coeff_matrix, rhs_vector) */
                       int *,     /* indx[0..n-1] (vector recording row permutation) */
                       double *); /* rhs_vector  [0..sys_size-1] */

#endif
