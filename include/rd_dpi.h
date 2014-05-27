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
 
#ifndef _RD_DPI_H
#define _RD_DPI_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _RD_DPI_C
#define EXTERN
#
#endif

#ifndef _RD_DPI_C
#define EXTERN extern
#endif

EXTERN int rd_dpi		/* rd_dpi.c */
PROTO((Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* fn - filename */
       const int ));		/* verbosity - how much to talk */

EXTERN void getdid		/* rd_dpi.c */
PROTO((int ,			/* netcdf_unit */
       char *,			/* string_name */
       int,			/* boolean for hard error interpretation */
       int *));			/* dimension_identifier_address */


EXTERN void getvid		/* rd_dpi.c */
PROTO((int ,			/* netcdf_unit */
       char *,			/* string_name */
       int,			/* boolean for hard error interpretation */
       int *));			/* variable_identifier_address */

EXTERN void getdim		/* rd_dpi.c */
PROTO((int ,			/* netcdf_unit */
       int ,			/* dimension_id */
       int *));			/* where -- to put the dimension value */

EXTERN void uni_dpi		/* rd_dpi.c                                  */
PROTO((Dpi *,			/* dpi                                       */
       Exo_DB *));		/* exo                                       */

EXTERN void free_dpi		/* rd_dpi.c */
PROTO((Dpi *));			/* fantastic structure defd in "dpi.h" */

EXTERN void free_dpi_uni	/* rd_dpi.c */
PROTO((Dpi *));			/* fantastic structure defd in "dpi.h" */

EXTERN void init_dpi_struct	/* rd_dpi.c */
PROTO((Dpi *));			/* fantastic structure defd in "dpi.h" */

#endif /* _RD_DPI_H */
