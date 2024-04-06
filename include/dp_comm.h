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

#ifndef DP_COMM_H
#define DP_COMM_H

#include "bc_surfacedomain.h"
#include "dp_types.h"
#include "dpi.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_DP_COMM_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_DP_COMM_C
#define EXTERN extern
#endif

EXTERN void exchange_dof(Comm_Ex *, /* cx - ptr to communications exchange info */
                         Dpi *,     /* dpi - distributed processing info */
                         double *,  /* x - local processor dof-based vector */
                         int);

EXTERN void exchange_dof_int(Comm_Ex *, /* cx - ptr to communications exchange info */
                             Dpi *,     /* dpi - distributed processing info */
                             int *,     /* x - local processor dof-based vector */
                             int);

EXTERN void exchange_dof_long_long(Comm_Ex *,   /* cx - ptr to communications exchange info */
                                   Dpi *,       /* dpi - distributed processing info */
                                   long long *, /* x - local processor dof-based vector */
                                   int);

EXTERN void exchange_node(Comm_Ex *cx, /* cx - ptr to communications exchange info */
                          Dpi *d,      /* dpi - distributed processing info */
                          double *a);  /* x - local processor node-based vector */

#endif /* GOMA_DP_COMM_H */
