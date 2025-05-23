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

#ifndef GOMA_MM_FILL_COMMON_H
#define GOMA_MM_FILL_COMMON_H

#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_COMMON_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_COMMON_C
#define EXTERN extern
#endif

EXTERN void computeCommonMaterialProps_gp(const dbl time // Current time (sec)
);

#endif
