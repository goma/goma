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

EXTERN void lu
(const int ,		/* N                                         */
       const int ,		/* NExt                                      */
       const int ,		/* M                                         */
       double [],		/* a                                         */
       int [],			/* ija                                       */
       double [],		/* x                                         */
       const int );		/* factor_flag                               */

#endif /* GOMA_SL_LU_H */
