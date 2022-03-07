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

#ifndef GOMA_DP_MAP_COMM_VEC_H
#define GOMA_DP_MAP_COMM_VEC_H

#include "dp_comm.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_DP_MAP_COMM_VEC_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_DP_MAP_COMM_VEC_C
#define EXTERN extern
#endif

extern void exchange_neighbor_proc_info(int, COMM_NP_STRUCT *);

extern void setup_nodal_comm_map(Exo_DB *, Dpi *, Comm_Ex **);

EXTERN void setup_dof_comm_map(Exo_DB *,      /* exo - ptr to FE database */
                               Dpi *,         /* dpi - ptr to distrib proc db */
                               Comm_Ex **cx); /* array of structures, one for ea neighbor */

extern void output_comm_stats(Dpi *, Comm_Ex **);

#endif /* GOMA_DP_MAP_COMM_VEC_H */
