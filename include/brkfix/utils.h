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

/* utils.h - prototype declarations for utils.c
 */

#ifndef GOMA_UTILS_H
#define GOMA_UTILS_H

#include "brkfix/exo_utils.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_UTILS_C
#define EXTERN extern
#endif

EXTERN int gcf			/* utils.c */
(int ,			/* first integer */
       int );			/* second integer */


EXTERN int count_node_node_interactions	/* utils.c */
(int ,			/* num_nodes */
       int *,			/* node_ptr */
       int *,			/* elem_list */
       int *,			/* elem_ptr */
       int *);			/* node_list */

EXTERN int findex_mono		/* utils.c */
(int ,			/* val - what integer value to seek */
       int *,			/* start - where to begin looking */
       int );			/* length - how far to look from start */


EXTERN int get_node_index	/* utils.c */
(int ,			/* global_node - what we're searching for */
       int *,			/* node_list - 3 part list we're searching */
       int ,			/* num_internal_nodes - len part 1 of list */
       int ,			/* num_boundary_nodes - len part 2 of list */
       int );			/* num_external_nodes - len part 3 of list */

EXTERN int get_internal_boundary_index /* utils.c */
(int ,			/* global_node */
       int *,			/* node_list */
       int ,			/* num_internal_nodes */
       int );			/* num_boundary_nodes */

EXTERN void proc_sort		/* utils.c */
(int *,			/* node_list - like, eg, the external nodes */
       int ,			/* len - node_list[0...len-1] */
       int ,			/* len_assignment */
       int *);		        /* assigment - how we will sort nodes */

EXTERN void isort		/* utils.c                                   */
(const int ,		/* length                                    */
       int *);			/* array                                     */


EXTERN int get_filename_num_procs /* utils.c */
(const char *);		/* basename - of polylithic files */

EXTERN int get_max_val_index	/* utils.c */
(const int *,		/* array                                (in) */
       const int ,		/* length                               (in) */
       int *,			/* maximum_value                       (out) */
       int *);			/* maximum_index                       (out) */

EXTERN int get_min_val_index	/* utils.c */
(const int *,		/* array                                (in) */
       const int ,		/* length                               (in) */
       int *,			/* minimum_value                       (out) */
       int *);			/* minimum_index                       (out) */

EXTERN int is_shell_type 
(char *elem_type);       // String that is the element type from exo

#endif /* GOMA_UTILS_H */
