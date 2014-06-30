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

/* exo_conn.h -- prototype declarations for exo_conn.c
 *
 */

#ifndef _EXO_CONN_H
#define _EXO_CONN_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _EXO_CONN_C
#define EXTERN
#
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


EXTERN int int_intersect
PROTO((int *,			/* a - 1st integer list			(in) */
       int *,			/* b - 2nd integer list			(in) */
       int ,			/* len_a - length of 1st integer list	(in) */
       int ,			/* len_b - length of 2nd integer list	(in) */
       int *,			/* ia - intersections indeces, 1st list (out)*/
       int *));			/* ib - intersections indeces, 2nd list (out)*/

EXTERN int sides2nodes
PROTO((int ,			/* face - 0,1,...,num_faces-1 */
       int ,			/* shape - LINE_SEGMENT, TRIANGLE, etc. */
       int *));			/* local_indeces - filled with right ones */

EXTERN int build_side_node_list
PROTO((const int ,		/* elem - the element number                 */
       const int ,		/* face - the face number                    */
       const Exo_DB *,		/* exo - ptr to whole mesh structure FE db*/
       int *));			/* snl - node list for this side (out) */

EXTERN int get_num_face_interactions
PROTO((char *elem_type));       // String for element type

#endif /* _EXO_CONN_H */
