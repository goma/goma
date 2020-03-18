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

#ifndef GOMA_USER_PRE_H
#define GOMA_USER_PRE_H

#include "std.h"
#include "user_post.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_USER_PRE_C
#define EXTERN
#
#endif

#ifndef GOMA_USER_PRE_C
#define EXTERN extern
#endif

EXTERN double user_surf_object
(int *,
       dbl *,    	/* param - ptr to user-defined list          */
       dbl *);		


#endif /* GOMA_USER_PRE_H */
