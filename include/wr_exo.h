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

#ifndef GOMA_WR_EXO_H
#define GOMA_WR_EXO_H

#include "exo_struct.h"
#include "wr_dpi.h"

struct Results_Description;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_EXO_C
#define EXTERN
#
#endif

#ifndef GOMA_WR_EXO_C
#define EXTERN extern
#endif

EXTERN int wr_mesh_exo /* wr_exo.c                                 */
    (Exo_DB *,         /* exo - ptr to full ripe EXODUS II fe db   */
     char *,           /* filename - where to write                */
     int);             /* verbosity - talk while writing           */

EXTERN void wr_result_prelim_exo   /* wr_exo.c                                 */
    (struct Results_Description *, /* rd - describe nodal variables        */
     Exo_DB *,                     /* exo - whole mesh                          */
     char *,                       /* filename - where to write                 */
     double ***);                  /* gvec_elem array, it gets sized & malloc'd *
                                    * in wr_exo.c                               */

EXTERN void wr_result_prelim_exo_segregated /* wr_exo.c                                 */
    (struct Results_Description **,         /* rd - describe nodal variables        */
     Exo_DB *,                              /* exo - whole mesh                          */
     char *,                                /* filename - where to write                 */
     double ****);                          /* gvec_elem array, it gets sized & malloc'd *
                                             * in wr_exo.c                               */

EXTERN void create_truth_table     /* wr_exo.c */
    (struct Results_Description *, /* rd - describe nodal variables        */
     Exo_DB *,                     /* filename - where to write            */
     double ***gvec_elem);         /* array holding elem values - final    *
                                    * dim gets malloc'd here               */

EXTERN void create_truth_table_segregated /* wr_exo.c */
    (struct Results_Description **,       /* rd - describe nodal variables        */
     Exo_DB *,                            /* filename - where to write            */
     double ****gvec_elem);               /* array holding elem values - final    *
                                           * dim gets malloc'd here               */

EXTERN void wr_nodal_result_exo /* wr_exo.c                                  */
    (Exo_DB *,                  /* exo - ptr to whole mesh                   */
     char *,                    /* filename - where to write this data       */
     double[],                  /* vector - of nodal values                  */
     int,                       /* nodal_variable_index - for exodus         *
                                 * (cf Results_Description                   */
     int,                       /* time_step                                 */
     double);                   /* time_value                                */

EXTERN void wr_elem_result_exo /* wr_exo.c                                  */
    (Exo_DB *,                 /* exo - ptr to whole mesh                   */
     const char *,             /* filename - where to write this data       */
     double ***,               /* vector [eb_indx][ev_indx][elem]           *
                                  - of element values                       */
     const int,                /* elem_variable_index - (0 based) for exodus*
                                * (cf Results_Description                   */
     const int,                /* time_step                                 */
     const double,             /* time_value                                */
     struct Results_Description *);

EXTERN void wr_global_result_exo /* wr_exo.c */
    (Exo_DB *,                   /* exo - ptr to mesh struct */
     const char *,               /* filename of exodus database file */
     const int,                  /* time_step */
     const int,                  /* number of globals to write */
     double[]);                  /* global value vector */

EXTERN void add_qa_stamp(Exo_DB *); /* exo                                       */

EXTERN void add_info_stamp(Exo_DB *); /* exo                                       */

void wr_result_exo(
    Exo_DB *exo, char *filename, int verbosity, int write_node_vars, int write_elem_vars);

void wr_resetup_exo(Exo_DB *exo, char *filename, int verbosity);

#endif /* GOMA_WR_EXO_H */
