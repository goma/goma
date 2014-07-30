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

/* exo_utils.h - prototype declarations for exo_utils.c
 */

#ifndef _EXO_UTILS_H
#define _EXO_UTILS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _EXO_UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef _EXO_UTILS_C
#define EXTERN extern
#endif

EXTERN void zero_base
PROTO((Exo_DB *));		/* E - ptr to EXODUS II FE db                */

EXTERN void one_base
PROTO((Exo_DB *));		/* E - ptr to EXODUS II FE db                */

EXTERN Element_shape get_element_shape
PROTO((const Exo_DB *,		/* exo                                       */
       const int ));		/* element                                   */

#if 0
EXTERN Element_shape type2shape
PROTO((const Element_type ));
#endif

EXTERN int shape2sides
PROTO((const Element_shape ));

EXTERN int find_element_block
PROTO(( Exo_DB *,   // Exodus database
	int ));     // Element number

EXTERN int find_element_friends  // Outputs number of friends found
PROTO(( Exo_DB *,   // Exodus database
	int,        // Element of which I want to find friends
	int * ));   // (out) List of elements who are friends

EXTERN int find_local_node_number
PROTO(( Exo_DB *,   // Exodus database
	int,        // Element to inspect
	int ));     // Global node number

EXTERN void find_edge_connected_nodes
PROTO(( int,        // Local node number
	int * ));   // Three edge-connected nodes

EXTERN int is_node_in_element
PROTO(( Exo_DB *,   // Exodus database
	int,        // Node number
	int ));     // Element number

#endif /* _EXO_UTILS_H */
