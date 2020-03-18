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
 
#ifndef GOMA_AC_HUNT_H
#define GOMA_AC_HUNT_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_HUNT_C
#define EXTERN
#
#endif

#ifndef GOMA_AC_HUNT_C
#define EXTERN extern
#endif

#include "dp_types.h"
#include "exo_struct.h"
#include "dpi.h"
#include "ac_conti.h"

EXTERN void hunt_problem
(Comm_Ex *,		/* array of communications structures */
       Exo_DB *,		/* ptr to the finite element mesh database */
       Dpi *);			/* distributed processing information */

#endif /* GOMA_AC_HUNT_H */
