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
 

/*
 *$Id: user_pre.h,v 5.1 2007-09-18 18:53:49 prschun Exp $
 */

#ifndef _USER_PRE_H
#define _USER_PRE_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _USER_PRE_C
#define EXTERN
#
#endif

#ifndef _USER_PRE_C
#define EXTERN extern
#endif

EXTERN double user_surf_object
(int *,
       dbl *,    	/* param - ptr to user-defined list          */
       dbl *);		


#endif /* _USER_PRE_H */
