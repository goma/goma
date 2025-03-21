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

#ifndef GOMA_USER_POST_H
#define GOMA_USER_POST_H

#include <stdio.h>

#include "std.h"
#include "user_mp_gen.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_USER_POST_C
#define EXTERN
#
#endif

#ifndef GOMA_USER_POST_C
#define EXTERN extern
#endif

EXTERN double user_post(dbl *); /* param - ptr to user-defined list          */

EXTERN int usr_ptracking(FILE *,         /*  jfp - filename for output                */
                         const int,      /*  part_id - particle id - starts at 1      */
                         const double[], /*  part_x - current particle coords         */
                         const double[], /*  part_v - current particle velocity       */
                         const double[], /*  part_xold - past coords                  */
                         const double[], /*  part_vold - past velocity                */
                         const int,      /*  heading - flag for writing headings      */
                         const double,   /*  time_value - particle time               */
                         const double);  /*  time_step - time step                    */

#endif /* GOMA_USER_POST_H */
