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

/* bbb.h - prototype declarations for bbb.c
 */

#ifndef GOMA_BBB_H
#define GOMA_BBB_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BBB_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_BBB_C
#define EXTERN extern
#endif

#include "brkfix/fix.h"
EXTERN void build_big_bones     /* bbb.c */
    (Exo_DB *,                  /* EXODUS info from representative polylith */
     Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void build_global_conn   /* bbb.c */
    (Exo_DB *,                  /* EXODUS info from representative polylith */
     Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void build_global_attr   /* bbb.c */
    (Exo_DB *,                  /* EXODUS info from representative polylith */
     Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void build_global_coords /* bbb.c */
    (Exo_DB *,                  /* EXODUS info from representative polylith */
     Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *);                 /* ptr to monolithic EXODUS II database */

EXTERN void build_global_ns     /* bbb.c */
    (Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void build_global_ss     /* bbb.c */
    (Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void build_global_res    /* bbb.c */
    (Exo_DB *,                  /* EXODUS info from representative polylith */
     Dpi *,                     /* distributed processing info from polylith */
     Exo_DB *,                  /* ptr to monolithic EXODUS II database */
     struct fix_data *);

EXTERN void mononame            /* bbb.c */
    (char *,                    /* in  - string to be mono-sized "a_b.c" */
     char *);                   /* out - result string "a.c"  */

#endif                          /* GOMA_BBB_H */
