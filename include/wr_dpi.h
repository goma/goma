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

EXTERN int wr_dpi
(Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* filename */
       int );			/* verbosity - how much to talk */

#endif /* GOMA_WR_DPI_H */
