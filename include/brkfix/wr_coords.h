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

#ifndef _WR_COORDS_H
#define _WR_COORDS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _WR_COORDS_C
#define EXTERN /* do nothing */
#endif

#ifndef _WR_COORDS_C
#define EXTERN extern
#endif

EXTERN void write_coords
PROTO((char *,			/* filename */
       Exo_DB *));		/* EXODUS II database */

#endif /* _WR_COORDS_H */
