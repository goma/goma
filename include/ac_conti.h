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

#ifndef GOMA_AC_CONTI_H
#define GOMA_AC_CONTI_H

#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_shell_bc.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_CONTI_C
#define EXTERN
#
#endif

#ifndef GOMA_AC_CONTI_C
#define EXTERN extern
#endif

EXTERN void continue_problem(Comm_Ex *, /* array of communications structures */
                             Exo_DB *,  /* ptr to the finite element mesh database */
                             Dpi *);    /* distributed processing information */

#endif                                  /* GOMA_AC_CONTI_H */
