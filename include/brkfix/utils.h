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

#ifndef _UTILS_H
#define _UTILS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef _UTILS_C
#define EXTERN extern
#endif

EXTERN int gcf			/* utils.c */
PROTO((int ,			/* first integer */
       int ));			/* second integer */


EXTERN int count_node_node_interactions	/* utils.c */
PROTO((int ,			/* num_nodes */
       int *,			/* node_ptr */
       int *,			/* elem_list */
       int *,			/* elem_ptr */
       int *));			/* node_list */

EXTERN int in_list		/* utils.c */
PROTO((int ,			/* val    - what integer value to seek */
       int *,			/* start  - where to begin looking */
       int ));			/* length - how far to search from start */


EXTERN int fence_post		/* utils.c */
PROTO((int ,			/* val    - integer whose category we seek */
       int *,			/* array  - where to look  */
       int ));			/* length - how far to search in array */

EXTERN int findex_mono		/* utils.c */
PROTO((int ,			/* val - what integer value to seek */
       int *,			/* start - where to begin looking */
       int ));			/* length - how far to look from start */


EXTERN int get_node_index	/* utils.c */
PROTO((int ,			/* global_node - what we're searching for */
       int *,			/* node_list - 3 part list we're searching */
       int ,			/* num_internal_nodes - len part 1 of list */
       int ,			/* num_boundary_nodes - len part 2 of list */
       int ));			/* num_external_nodes - len part 3 of list */

EXTERN int get_internal_boundary_index /* utils.c */
PROTO((int ,			/* global_node */
       int *,			/* node_list */
       int ,			/* num_internal_nodes */
       int ));			/* num_boundary_nodes */

EXTERN void proc_sort		/* utils.c */
PROTO((int *,			/* node_list - like, eg, the external nodes */
       int ,			/* len - node_list[0...len-1] */
       int ,			/* len_assignment */
       int *));		        /* assigment - how we will sort nodes */

EXTERN void zero_base		/* utils.c */
PROTO((Exo_DB *));		/* E - pointer to an EXODUS II database */

EXTERN void one_base		/* utils.c */
PROTO((Exo_DB *));		/* x = pointer to an EXODUS II FE database */

EXTERN void isort		/* utils.c                                   */
PROTO((const int ,		/* length                                    */
       int *));			/* array                                     */


EXTERN int get_filename_num_procs /* utils.c */
PROTO((const char *));		/* basename - of polylithic files */

EXTERN int get_max_val_index	/* utils.c */
PROTO((const int *,		/* array                                (in) */
       const int ,		/* length                               (in) */
       int *,			/* maximum_value                       (out) */
       int *));			/* maximum_index                       (out) */

EXTERN int get_min_val_index	/* utils.c */
PROTO((const int *,		/* array                                (in) */
       const int ,		/* length                               (in) */
       int *,			/* minimum_value                       (out) */
       int *));			/* minimum_index                       (out) */

EXTERN int is_shell_type 
PROTO((char *elem_type));       // String that is the element type from exo

EXTERN int is_shell_element
PROTO((Exo_DB *,                // Exodus database
       int e));                 // Element number from database

#endif /* _UTILS_H */
