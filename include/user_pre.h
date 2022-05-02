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

EXTERN double user_surf_object(int *,
                               dbl *, /* param - ptr to user-defined list          */
                               dbl *);

EXTERN double user_mat_init(const int,    /* variable                            */
                            const int,    /* node number                            */
                            const dbl,    /* Basic initial value                    */
                            const dbl[],  /* p                                      */
                            const dbl[],  /* nodal coordinates                      */
                            const int,    /* material ID                            */
                            const dbl[]); /* other variable values                 */

#endif /* GOMA_USER_PRE_H */
