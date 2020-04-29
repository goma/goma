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
 * $Id: user_continuation.h,v 5.1 2007-09-18 18:53:48 prschun Exp $
 */

#ifndef GOMA_USER_CONTINUATION_H
#define GOMA_USER_CONTINUATION_H

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "wr_side_data.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_USER_CONTINUATION_C
#define EXTERN
#
#endif

#ifndef GOMA_USER_CONTINUATION_C
#define EXTERN extern
#endif

EXTERN void update_user_parameter
(double,			/* PARAMETER VALUE */
       double*,                 /* UNKNOWN VECTOR */
       double*,                 /* UNKNOWN_DOT VECTOR */
       double*,                 /* x_AC VECTOR */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_user_TP_parameter
(double,			/* PARAMETER VALUE */
       double*,                 /* UNKNOWN VECTOR */
       double*,                 /* UNKNOWN_DOT VECTOR */
       double*,                 /* x_AC VECTOR */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN int do_user_update
(int,			/* Ptr to function number */
       int,			/* First update call (first param) */
       int,			/* First update call (second param) */
       int,			/* Type */
       int,			/* BCID */
       int,			/* DFID */
       int,			/* MTID */
       int,			/* MPID */
       int,			/* MDID */
       double,			/* Value */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

#endif /* GOMA_USER_CONTINUATION_H */
