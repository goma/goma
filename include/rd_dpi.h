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
 
#ifndef GOMA_RD_DPI_H
#define GOMA_RD_DPI_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RD_DPI_C
#define EXTERN
#
#endif

#ifndef GOMA_RD_DPI_C
#define EXTERN extern
#endif

EXTERN int rd_dpi		/* rd_dpi.c */
(Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* fn - filename */
       const int );		/* verbosity - how much to talk */

EXTERN void getdid		/* rd_dpi.c */
(int ,			/* netcdf_unit */
       char *,			/* string_name */
       int,			/* boolean for hard error interpretation */
       int *);			/* dimension_identifier_address */


EXTERN void getvid		/* rd_dpi.c */
(int ,			/* netcdf_unit */
       char *,			/* string_name */
       int,			/* boolean for hard error interpretation */
       int *);			/* variable_identifier_address */

EXTERN void getdim		/* rd_dpi.c */
(int ,			/* netcdf_unit */
       int ,			/* dimension_id */
       int *);			/* where -- to put the dimension value */

EXTERN void uni_dpi		/* rd_dpi.c                                  */
(Dpi *,			/* dpi                                       */
       Exo_DB *);		/* exo                                       */

EXTERN void free_dpi		/* rd_dpi.c */
(Dpi *);			/* fantastic structure defd in "dpi.h" */

EXTERN void free_dpi_uni	/* rd_dpi.c */
(Dpi *);			/* fantastic structure defd in "dpi.h" */

EXTERN void init_dpi_struct	/* rd_dpi.c */
(Dpi *);			/* fantastic structure defd in "dpi.h" */

void 
exo_dpi_clone(Exo_DB *exo, 
	      Dpi *dpi);

#endif /* GOMA_RD_DPI_H */
