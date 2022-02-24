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

/* wr_graph_file.h - prototype declarations for wr_graph_file.c
 */

#ifndef GOMA_WR_GRAPH_FILE_H
#define GOMA_WR_GRAPH_FILE_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_WR_GRAPH_FILE_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_WR_GRAPH_FILE_C
#define EXTERN extern
#endif

EXTERN void wr_graph_file /* wr_graph_file.c */
    (char *,              /* gfn - graph file name */
     char *,              /* ifn - input file name */
     char *,              /* efn - .exoII file name */
     int,                 /* number of elements */
     int,                 /* number of fe nodes */
     int,                 /* num dofs whole problem */
     int,                 /* total number of nonzeroes */
     int,                 /* total num assembled terms */
     int,                 /* boolean for graphfile fmt */
     int *,               /* ptr to node-node stuff */
     int[],               /* mins, maxes, and scales */
     int *,               /* to find neighbors */
     int *,               /* for vertex weights */
     int *);              /* for edge weights */

#endif /* GOMA_WR_GRAPH_FILE_H */
