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

/* bbb.h - prototype declarations for bbb.c
 */

#ifndef GOMA_BBB_H
#define GOMA_BBB_H

#include "dpi.h"
#include "exo_struct.h"
#include "wr_side_data.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BBB_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_BBB_C
#define EXTERN extern
#endif

EXTERN void build_big_bones	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_conn	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_attr	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_coords	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_ns	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_ss	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_res	/* bbb.c */
(Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *);		/* ptr to monolithic EXODUS II database */

EXTERN void mononame		/* bbb.c */
(char *,			/* in  - string to be mono-sized "a_b.c" */
       char *);		/* out - result string "a.c"  */

#endif /* GOMA_BBB_H */
