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
 
#ifndef GOMA_MM_BC_H
#define GOMA_MM_BC_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_augc_util.h"
#include "rf_vars_const.h"

struct BC_descriptions;
struct Boundary_Condition;
struct elem_edge_bc_struct;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_BC_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_BC_C
#define EXTERN extern
#endif

EXTERN void find_and_set_Dirichlet
(double [],		/* x - solution vector at this processor     */
       double [],		/* xdot - time derivative of solution vector */
       Exo_DB *,		/* exo - ptr to EXODUS II FE database        */
       Dpi *);			/* dpi - ptr to distrib proc FE database     */

EXTERN void alloc_First_Elem_BC
(struct elem_side_bc_struct ****, /* First_Elem_Side_BC_Array ptr       */
       struct elem_edge_bc_struct ****, /* First_Elem_Edge_BC_Array ptr       */
       const int );	       	       /* num_internal_elems                 */

EXTERN void set_up_Surf_BC
(struct elem_side_bc_struct **[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *);			/* dpi - ptr to dist proc info               */
	   
EXTERN void free_Surf_BC
        (struct elem_side_bc_struct **First_Elem_Side_BC_Array[], Exo_DB *exo);			/* dpi - ptr to dist proc info               */

	   
EXTERN void free_Edge_BC
( struct elem_edge_bc_struct **[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *);			/* dpi - ptr to dist proc info               */

EXTERN void setup_Point_BC
( struct elem_side_bc_struct **[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *);			/* dpi - ptr to dist proc info               */

EXTERN void set_up_Embedded_BC
(void );

EXTERN void set_up_Edge_BC
(struct elem_edge_bc_struct **[ ],	/* First_Elem_Edge_BC_Array          */
       Exo_DB *,		/* exo                                       */
       Dpi *);			/* dpi                                       */

EXTERN void check_for_bc_conflicts2D
(Exo_DB *,		/* exo                                       */
       Dpi *);			/* dpi                                       */

EXTERN void check_for_bc_conflicts3D
(Exo_DB *,		/* exo                                       */
       Dpi *);			/* dpi                                       */

extern void initialize_Boundary_Condition(struct Boundary_Condition *);

EXTERN int find_id_side
(const int ,		/* ielem - element index number              */
       const int ,		/* num_nodes_on_side - guess what?      (in) */
       const int [],		/* local_ss_node_list - nodes, ordered  (in) */
       int [],			/* id_local_elem_coord - (xi,eta,zeta)  (in) */
       const Exo_DB *);	/* exo - ptr to FE db                   (in) */

EXTERN int find_id_side_BC
(const int ,		/* ielem - element index number              */
       const int ,		/* num_nodes_on_side - guess what?      (in) */
       const int [],		/* local_ss_node_list - nodes, ordered  (in) */
       const int ,		/* ibc - BC index                       (in) */
       int [],			/* id_local_elem_coord - (xi,eta,zeta)  (in) */
       const Exo_DB *);	/* exo - ptr to FE db                   (in) */

EXTERN int find_id_side_SS
(const int,               /* ielem - element index number         (in) */
       const int,               /* iss - sideset index number           (in) */
       const Exo_DB *);        /* exo - ptr to FE db                   (in) */

EXTERN void print_setup_Surf_BC
(struct elem_side_bc_struct *[ ]); /* First_Elem_Side_BC_Array        */

EXTERN struct BC_descriptions *alloc_BC_description
(struct BC_descriptions *); /* old_ptr - serves as pattern for new    */

EXTERN char *rot_eq_type2string
(const int );		/* eq_type - change from into string    (in) */

EXTERN char *rotopology_type2string
(const int );		/* type - change from int to string     (in) */


extern int find_bc_unk_offset(struct Boundary_Condition *, int,
			      int, int, int *, 
			      VARIABLE_DESCRIPTION_STRUCT **);

extern int search_bc_dup_list(const int, int *);
extern void set_up_BC_connectivity(void);

int exchange_bc_info(void);
#endif /* GOMA_MM_BC_H */
