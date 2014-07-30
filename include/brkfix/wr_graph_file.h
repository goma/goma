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

/* wr_graph_file.h - prototype declarations for wr_graph_file.c
 */

#ifndef _WR_GRAPH_FILE_H
#define _WR_GRAPH_FILE_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _WR_GRAPH_FILE_C
#define EXTERN /* do nothing */
#endif

#ifndef _WR_GRAPH_FILE_C
#define EXTERN extern
#endif

EXTERN void wr_graph_file			/* wr_graph_file.c */
PROTO((char *,					/* gfn - graph file name */
       char *,					/* ifn - input file name */
       char *,					/* efn - .exoII file name */
       int,					/* number of elements */
       int, 					/* number of fe nodes */
       int,					/* num dofs whole problem */
       int,					/* total number of nonzeroes */
       int,					/* total num assembled terms */
       int,					/* boolean for graphfile fmt */
       int *,					/* ptr to node-node stuff */
       int [],					/* mins, maxes, and scales */
       int *,					/* to find neighbors */
       int *,					/* for vertex weights */
       int *));					/* for edge weights */


#endif /* _WR_GRAPH_FILE_H */
