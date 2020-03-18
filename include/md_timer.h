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

EXTERN dbl ut		/* user time in seconds */
(void);

EXTERN double ust		/* user + system time in seconds */
(void);

EXTERN void get_date
(char *);		/* string - fill in with mm/dd/yy */

EXTERN void get_time
(char *);		/* string - fill in with hh:mm:ss */

#endif /* GOMA_MD_TIMER_H */





