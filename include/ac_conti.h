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
 
#ifndef _AC_CONTI_H
#define _AC_CONTI_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _AC_CONTI_C
#define EXTERN
#
#endif

#ifndef _AC_CONTI_C
#define EXTERN extern
#endif

EXTERN void continue_problem
(Comm_Ex *,		/* array of communications structures */
       Exo_DB  *,		/* ptr to the finite element mesh database */
       Dpi     *);		/* distributed processing information */

#endif /* _AC_CONTI_H */
