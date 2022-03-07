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

#include "ac_conti.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"

EXTERN void hunt_problem(Comm_Ex *, /* array of communications structures */
                         Exo_DB *,  /* ptr to the finite element mesh database */
                         Dpi *);    /* distributed processing information */

#endif /* GOMA_AC_HUNT_H */
