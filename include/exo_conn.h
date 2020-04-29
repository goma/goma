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
 
#ifndef GOMA_EXO_CONN_H
#define GOMA_EXO_CONN_H

#include "dpi.h"
#include "el_quality.h"
#include "exo_struct.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_EXO_CONN_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_EXO_CONN_C
#define EXTERN extern
#endif

EXTERN void build_elem_node	/* exo_conn.c */
(Exo_DB *);		/* exo - ptr to EXODUS II database struct */

EXTERN void build_node_elem	/* exo_conn.c */
(Exo_DB *);		/* exo - ptr to EXODUS II database struct */

EXTERN void build_node_node	/* exo_conn.c */
(Exo_DB *);		/* exo - ptr to EXODUS II database struct */

EXTERN void build_elem_elem	/* exo_conn.c */
(Exo_DB *);		/* exo - ptr to EXODUS II database struct */

EXTERN int build_side_node_list
(int ,			/* elem - the element number */
       int ,			/* face - the face number */
       Exo_DB *,		/* exo - ptr to whole mesh structure FE db*/
       int *);			/* snl - node list for this side (out) */


EXTERN int get_exterior_faces
( int,
	int *,
	const Exo_DB *,
	const Dpi * );

void
brk_build_node_node(Exo_DB *exo);
#endif /* GOMA_EXO_CONN_H */
