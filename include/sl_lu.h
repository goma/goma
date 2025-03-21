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

#ifndef GOMA_SL_LU_H
#define GOMA_SL_LU_H

#include "sl_auxutil.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_SL_LU_C
#define EXTERN
#
#endif

#ifndef GOMA_SL_LU_C
#define EXTERN extern
#endif

EXTERN void lu(const int,  /* N                                         */
               const int,  /* NExt                                      */
               const int,  /* M                                         */
               double[],   /* a                                         */
               int[],      /* ija                                       */
               double[],   /* x                                         */
               const int); /* factor_flag                               */

#endif /* GOMA_SL_LU_H */
