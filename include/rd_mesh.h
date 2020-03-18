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
 
#ifndef GOMA_RD_MESH_H
#define GOMA_RD_MESH_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "rd_exo.h"
#include "rf_bc_const.h"

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_RD_MESH_C
#define EXTERN
#
#endif

#ifndef GOMA_RD_MESH_C
#define EXTERN extern
#endif

EXTERN int * find_ss_internal_boundary(Exo_DB *e);

EXTERN int find_mat_number	/* rd_mesh.c */
(const int ,		/* ielem - element number (0 based) */
       const Exo_DB *);	/* exo - the whole mesh */


EXTERN int eb_in_matrl          /* rd_mesh.c */
(const int,               /* ebid */
       const int);             /* mn */

EXTERN int find_elemblock_index	/* rd_mesh.c */
(const int ,		/* element - number or index */
       const Exo_DB *);	/* exo - ptr to FE database*/

EXTERN int read_mesh_exoII	/* rd_mesh.c */
(Exo_DB *,		/* ptr to EXOII mesh datastructure */
       Dpi    *);		/* ptr to distributed processing info d.s. */

EXTERN void setup_old_exo	/* rd_mesh.c */
(Exo_DB *,		/* ptr to EXODUS II mesh database */
       Dpi *,
       int);

EXTERN void check_sidesets	/* rd_mesh.c */
(Exo_DB *,		/* EXODUS II FE db has all mesh info    (in) */
       struct Boundary_Condition [],	/* bct BC info                  (in) */
       int ,			/* number of boundary conditions        (in) */
       Dpi *);			/* distributed processing info          (in) */

EXTERN void check_nodesets	/* rd_mesh.c */
(Exo_DB *,		/* EXODUS II FE db has all mesh info    (in) */
       struct Boundary_Condition [], /* BC info                         (in) */
       int,			/* number of boundary conditions        (in) */
       Dpi *);			/* distributed processing info          (in) */


EXTERN void check_elemblocks	/* rd_mesh.c */
(Exo_DB *,		/* EXODUS II FE db has all mesh info (in) */
       int,			/* number of materials               (in) */
       struct Problem_Description *[], /* for material info          (in) */
       Dpi *);			/* distributed processing info       (in) */

EXTERN void setup_old_dpi	/* rd_mesh.c */
(Exo_DB *,		/* ptr to EXODUS II FE database */
       Dpi    *);		/* ptr to Distributed Processsing Info */

extern int ebID_to_ebIndex(const int);

EXTERN void setup_matilda	/* rd_mesh.c */
(Exo_DB *,		/* EXODUS II FE db has all mesh info (in) */
       int *);			/* distributed processing info       (out) */

EXTERN void sseb_conn		/* rd_mesh.c */
(Exo_DB *,		/* see exo_struct.h for full def         (in)*/
       int **,			/* side_set_pointers - ptrs to eb_list  (out)*/
       int **,			/* element_block_list - ebs for all ss's(out)*/
       int **,			/* element_block_pointers - ptrs ss_list(out)*/
       int **);		/* side_set_list - ss's for all eb's    (out)*/

EXTERN void build_list		/* rd_mesh.c */
(int,			/* prospective_member */
       int **,			/* list */
       int *,			/* current_size */
       int *);			/* current_max_size */

EXTERN void multiname		/* rd_mesh.c */
(char *,			/* in_name - generic global name "pref.suf" */
       int,			/* integer processor_name */
       int);			/* number_processors - total */

EXTERN void strip_suffix	/* rd_mesh.c */
(char *,			/* result - "a" */
       char *);		/* in -- input string "a.b" */

EXTERN int get_suffix		/* rd_mesh.c */
(char *,			/* result -- "b" */
       const char *);		/* in -- extract the tail of "a.b" -> ".b" */

EXTERN int get_prefix		/* rd_mesh.c */
(char *,			/* result -- "b" */
       const char *);		/* in -- extract the front of "a.b" -> "a" */

EXTERN int Elem_Type		/* rd_mesh.c */
(const Exo_DB *exo,	/* the mesh */
       const int );		/* element - the element number on this proc */

EXTERN void integer_sort		/* rd_mesh.c */
(int ,			/* length - of the integer vector */
       int *);			/* array - the vector to be sorted */

EXTERN int map_mat_index	/* rd_mesh.c                                 */
( const int );		/* ebid - element block ID                   */

#endif /* GOMA_RD_MESH_H */
