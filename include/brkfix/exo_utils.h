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

#ifndef GOMA_EXO_UTILS_H
#define GOMA_EXO_UTILS_H

#include "exo_struct.h"

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_EXO_UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_EXO_UTILS_C
#define EXTERN extern
#endif

EXTERN int find_element_friends  // Outputs number of friends found
( Exo_DB *,   // Exodus database
	int,        // Element of which I want to find friends
	int * );   // (out) List of elements who are friends

EXTERN int find_local_node_number
( Exo_DB *,   // Exodus database
	int,        // Element to inspect
	int );     // Global node number

EXTERN void find_edge_connected_nodes
( int,        // Local node number
	int * );   // Three edge-connected nodes

EXTERN int is_node_in_element
( Exo_DB *,   // Exodus database
	int,        // Node number
	int );     // Element number

#endif /* GOMA_EXO_UTILS_H */
