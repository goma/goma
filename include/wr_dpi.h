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

#ifndef GOMA_WR_DPI_H
#define GOMA_WR_DPI_H

#include "dpi.h"
#include "user_pre.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_DPI_C
#define EXTERN
#
#endif

#ifndef GOMA_WR_DPI_C
#define EXTERN extern
#endif

EXTERN int wr_dpi(Dpi *d, char *filename); /* verbosity - how much to talk */

#endif                                     /* GOMA_WR_DPI_H */
