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
 
#ifndef _EXO_CONN_H
#define _EXO_CONN_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _EXO_CONN_C
#define EXTERN /* do nothing */
#endif

#ifndef _EXO_CONN_C
#define EXTERN extern
#endif

EXTERN void build_elem_node	/* exo_conn.c */
PROTO((Exo_DB *));		/* exo - ptr to EXODUS II database struct */

EXTERN void build_node_elem	/* exo_conn.c */
PROTO((Exo_DB *));		/* exo - ptr to EXODUS II database struct */

EXTERN void build_node_node	/* exo_conn.c */
PROTO((Exo_DB *));		/* exo - ptr to EXODUS II database struct */

EXTERN void build_elem_elem	/* exo_conn.c */
PROTO((Exo_DB *));		/* exo - ptr to EXODUS II database struct */

EXTERN int get_exterior_faces
PROTO(( int,
	int *,
	const Exo_DB *,
	const Dpi * ));

#endif /* _EXO_CONN_H */
