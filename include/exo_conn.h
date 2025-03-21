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

EXTERN void build_elem_node               /* exo_conn.c */
    (Exo_DB *);                           /* exo - ptr to EXODUS II database struct */

EXTERN void build_node_elem               /* exo_conn.c */
    (Exo_DB *);                           /* exo - ptr to EXODUS II database struct */

EXTERN void build_node_node               /* exo_conn.c */
    (Exo_DB *);                           /* exo - ptr to EXODUS II database struct */

EXTERN void build_elem_elem               /* exo_conn.c */
    (Exo_DB *);                           /* exo - ptr to EXODUS II database struct */

EXTERN int build_side_node_list(int,      /* elem - the element number */
                                int,      /* face - the face number */
                                Exo_DB *, /* exo - ptr to whole mesh structure FE db*/
                                int *);   /* snl - node list for this side (out) */

EXTERN int get_exterior_faces(int, int *, const Exo_DB *, const Dpi *);

int sides2nodes(const int face,      /* assume face number 0,1,...,num_faces-1 */
                const int shape,     /* one of LINE_SEGMENT, etc. */
                int *local_indeces); /* get filled with right ones */

#endif                               /* GOMA_EXO_CONN_H */
