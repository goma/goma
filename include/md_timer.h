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
 
#ifndef _MD_TIMER_H
#define _MD_TIMER_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MD_TIMER_C
#define EXTERN /* do nothing */
#endif

#ifndef _MD_TIMER_C
#define EXTERN extern
#endif

EXTERN dbl ut		/* user time in seconds */
PROTO((void));

EXTERN double ust		/* user + system time in seconds */
PROTO((void));

EXTERN void get_date
PROTO((char *));		/* string - fill in with mm/dd/yy */

EXTERN void get_time
PROTO((char *));		/* string - fill in with hh:mm:ss */

#endif /* _MD_TIMER_H */





