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

/* wr_coords.h - prototype declarations for wr_coords.c
 */

#ifndef GOMA_WR_COORDS_H
#define GOMA_WR_COORDS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_COORDS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_WR_COORDS_C
#define EXTERN extern
#endif

EXTERN void write_coords(char *,    /* filename */
                         Exo_DB *); /* EXODUS II database */

#endif /* GOMA_WR_COORDS_H */
