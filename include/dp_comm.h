/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2026 Goma Developers, National Technology & Engineering   *
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

GOMA_EXTERN void exchange_dof(Comm_Ex *, /* cx - ptr to communications exchange info */
                         Dpi *,     /* dpi - distributed processing info */
                         double *,  /* x - local processor dof-based vector */
                         int);

GOMA_EXTERN void exchange_dof_int(Comm_Ex *, /* cx - ptr to communications exchange info */
                             Dpi *,     /* dpi - distributed processing info */
                             int *,     /* x - local processor dof-based vector */
                             int);

GOMA_EXTERN void exchange_dof_long_long(Comm_Ex *,   /* cx - ptr to communications exchange info */
                                   Dpi *,       /* dpi - distributed processing info */
                                   long long *, /* x - local processor dof-based vector */
                                   int);

GOMA_EXTERN void exchange_node(Comm_Ex *cx, /* cx - ptr to communications exchange info */
                          Dpi *d,      /* dpi - distributed processing info */
                          double *a);  /* x - local processor node-based vector */

#endif /* GOMA_DP_COMM_H */
