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

#ifndef _BBB_H
#define _BBB_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _BBB_C
#define EXTERN /* do nothing */
#endif

#ifndef _BBB_C
#define EXTERN extern
#endif

EXTERN void build_big_bones	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_conn	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_attr	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_coords	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_ns	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_ss	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void build_global_res	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

EXTERN void mononame		/* bbb.c */
PROTO((char *,			/* in  - string to be mono-sized "a_b.c" */
       char *));		/* out - result string "a.c"  */

#endif /* _BBB_H */
