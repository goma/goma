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

#ifndef GOMA_DP_VIF_H
#define GOMA_DP_VIF_H

#include "dp_map_comm_vec.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_DP_VIF_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_DP_VIF_C
#define EXTERN extern
#endif

EXTERN void noahs_raven /* dp_vif.c */
    (void);

EXTERN void noahs_ark /* dp_vif.c */
    (void);

EXTERN void noahs_dove /* dp_vif.c */
    (void);

EXTERN void raven_landing /* dp_vif.c */
    (void);

EXTERN void ark_landing /* dp_vif.c */
    (void);

#endif /* GOMA_DP_VIF_H */
