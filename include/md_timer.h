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

#ifndef GOMA_MD_TIMER_H
#define GOMA_MD_TIMER_H

#include "exo_conn.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MD_TIMER_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MD_TIMER_C
#define EXTERN extern
#endif

EXTERN dbl ut                 /* user time in seconds */
    (void);

EXTERN double ust             /* user + system time in seconds */
    (void);

EXTERN void get_date(char *); /* string - fill in with mm/dd/yy */

EXTERN void get_time(char *); /* string - fill in with hh:mm:ss */

#endif                        /* GOMA_MD_TIMER_H */
