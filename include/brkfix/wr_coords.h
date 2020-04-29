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

/* wr_coords.h - prototype declarations for wr_coords.c
 */

#ifndef GOMA_WR_COORDS_H
#define GOMA_WR_COORDS_H

#include "exo_struct.h"
#include "wr_side_data.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_COORDS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_WR_COORDS_C
#define EXTERN extern
#endif

EXTERN void write_coords
(char *,			/* filename */
       Exo_DB *);		/* EXODUS II database */

#endif /* GOMA_WR_COORDS_H */
