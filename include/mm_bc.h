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
 
#ifndef _MM_BC_H
#define _MM_BC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_BC_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_BC_C
#define EXTERN extern
#endif

EXTERN void find_and_set_Dirichlet
PROTO((double [],		/* x - solution vector at this processor     */
       double [],		/* xdot - time derivative of solution vector */
       Exo_DB *,		/* exo - ptr to EXODUS II FE database        */
       Dpi *));			/* dpi - ptr to distrib proc FE database     */

EXTERN void alloc_First_Elem_BC
PROTO((struct elem_side_bc_struct ***, /* First_Elem_Side_BC_Array ptr       */
       struct elem_edge_bc_struct ***, /* First_Elem_Edge_BC_Array ptr       */
       const int ));	       	       /* num_internal_elems                 */

EXTERN void set_up_Surf_BC
PROTO((struct elem_side_bc_struct *[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *));			/* dpi - ptr to dist proc info               */
	   
EXTERN void free_Surf_BC
PROTO(( struct elem_side_bc_struct *[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *));			/* dpi - ptr to dist proc info               */

	   
EXTERN void free_Edge_BC
PROTO(( struct elem_edge_bc_struct *[ ],	/* First_Elem_Side_BC_Array          */
       Exo_DB *,		/* exo - ptr to FE db                        */
       Dpi *));			/* dpi - ptr to dist proc info               */
       
EXTERN void set_up_Embedded_BC
PROTO((void ));

EXTERN void set_up_Edge_BC
PROTO((struct elem_edge_bc_struct *[ ],	/* First_Elem_Edge_BC_Array          */
       Exo_DB *,		/* exo                                       */
       Dpi *));			/* dpi                                       */

EXTERN void check_for_bc_conflicts2D
PROTO((Exo_DB *,		/* exo                                       */
       Dpi *));			/* dpi                                       */

EXTERN void check_for_bc_conflicts3D
PROTO((Exo_DB *,		/* exo                                       */
       Dpi *));			/* dpi                                       */

extern void initialize_Boundary_Condition(struct Boundary_Condition *);

EXTERN int find_id_side
PROTO((const int ,		/* ielem - element index number              */
       const int ,		/* num_nodes_on_side - guess what?      (in) */
       const int [],		/* local_ss_node_list - nodes, ordered  (in) */
       int [],			/* id_local_elem_coord - (xi,eta,zeta)  (in) */
       const Exo_DB *));	/* exo - ptr to FE db                   (in) */

EXTERN int find_id_side_BC
PROTO((const int ,		/* ielem - element index number              */
       const int ,		/* num_nodes_on_side - guess what?      (in) */
       const int [],		/* local_ss_node_list - nodes, ordered  (in) */
       const int ,		/* ibc - BC index                       (in) */
       int [],			/* id_local_elem_coord - (xi,eta,zeta)  (in) */
       const Exo_DB *));	/* exo - ptr to FE db                   (in) */

EXTERN int find_id_side_SS
PROTO((const int,               /* ielem - element index number         (in) */
       const int,               /* iss - sideset index number           (in) */
       const Exo_DB *));        /* exo - ptr to FE db                   (in) */

EXTERN void print_setup_Surf_BC
PROTO((struct elem_side_bc_struct *[ ])); /* First_Elem_Side_BC_Array        */

EXTERN struct BC_descriptions *alloc_BC_description
PROTO((struct BC_descriptions *)); /* old_ptr - serves as pattern for new    */

EXTERN char *rot_eq_type2string
PROTO((const int ));		/* eq_type - change from into string    (in) */

EXTERN char *rotopology_type2string
PROTO((const int ));		/* type - change from int to string     (in) */


extern int find_bc_unk_offset(struct Boundary_Condition *, int,
			      int, int, int *, 
			      VARIABLE_DESCRIPTION_STRUCT **);

extern int search_bc_dup_list(const int, int *);
extern void set_up_BC_connectivity(void);
#endif /* _MM_BC_H */
