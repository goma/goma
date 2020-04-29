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
 
#ifndef GOMA_BC_ROTATE_H
#define GOMA_BC_ROTATE_H

#include "bc_integ.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"

struct Aztec_Linear_Solver_System;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_ROTATE_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_BC_ROTATE_C
#define EXTERN extern
#endif

EXTERN int apply_rotated_bc
(double [],		/* resid_vector                              */
       struct elem_side_bc_struct *[], /* first_elem_side_BC_array
					* An array of pointers to the first 
					* surface integral defined for each 
					* element. Its length is equal to the 
					* total number of elements defined on 
					* the current processor.             */
       const int ,		/* ielem - element number                    */
       const int ,		/* ielem_type - element type                 */
       const int ,		/* num_local_nodes                           */
       const int ,		/* ielem_dim                                 */
       const int ,		/* iconnect_ptr                              */
       const int ,		/* num_total_nodes                           */
       const Exo_DB *);	/* exo - ptr to FE database                  */

EXTERN void rotate_mesh_eqn
(const int ,		/* id - Elemental stiffness matrix row index */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for the current element */
       const int ,		/* dim - physical dim of problem             */
       struct Aztec_Linear_Solver_System *);

EXTERN void rotate_res_jac_mesh
(int ,			/* irow_index - Elemental stiffness matrix   *
				 * row index                                 */
       int ,			/* Global node number                   */
       int ,			/* Pointer to beginning of connectivity list
				   for the current element	         */
       int ,			/* number of nodes in the element       */
       int ,			/* physical dim of elem surface - 0,1,2 */
       double [MAX_PDIM],	/* vector components of surface normal  */
       double [MAX_PDIM][MAX_PDIM][MDE], /*  dsnormal_dx */
				/* Array of derivatives ([i][j][k]) of
				   surface normal: component-i
				   wrt displacement-j at node-k         */
       double [2][MAX_PDIM],	/* stangent - vector components of 2 mutually*
				 * orthogonal surface tangent vectors        */
       double [2][MAX_PDIM][MAX_PDIM][MDE], /* dstangent_dx - Array of       *
					     * derivatives ([i][j][k]) of    *
					     * surface tangent vector 1 or 2 *
					     * component-i wrt 
					     * displacement-j at node-k      */
       double );		/* sign - of tangent vector                  */

EXTERN void rotate_res_jac_mom
(int ,			/* irow_index - Elemental stiffness matrix   *
				 * row index                                 */
       int ,			/* I - Global node number                    */
       int ,  		        /* iconnect_ptr - Pointer to beginning of    *
				 * connectivity list for current element     */
       int ,			/* nodes_per_elem - number of nodes in elem  */
       int ,			/* ielem_surf_dim - physical dim of elem     *
				 * surface - 0,1,2                           */
       double [MAX_PDIM],	/* snormal - vector components of surface    *
				 * normal                                    */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dsnormal_dx */
                                         /* Array of derivatives ([i][j][k]) of
                                            surface normal: component-i
                                            wrt displacement-j at node-k     */
       double [2][MAX_PDIM],	/* stangent - vector components of 2         *
				 * mutually orthogonal surface tangent       *
				 * vectors                                   */
       double [2][MAX_PDIM][MAX_PDIM][MDE], /* dstangent_dx  */
       /* Array of derivatives ([i][j][k]) of surface
	  tangent vector 1 or 2: component-i
	  wrt displacement-j at node-k         */
       double  );		/* sign - of tangent vector                  */

EXTERN void rotate_momentum_eqn
(const int ,		/* id - Elemental stiffness matrix row index */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for the current element */
       const int ,		/* dim - physical dim of problem             */
       struct Aztec_Linear_Solver_System *);	

EXTERN void calculate_all_rotation_vectors
(Exo_DB *,		/* exo - ptr to the whole FE mesh */
       double []);		/* x - Soln vector */

EXTERN void calculate_2D_rotation_vectors
(Exo_DB *,		/* exo - ptr to the whole FE mesh */
       double []);		/* x - Soln vector */

									
EXTERN void append_vectors
( int , 
		 int ,
		 int ,
		 int *,
		 ROTATION_VECTORS_STRUCT ** );

EXTERN void get_ss_element_list
(const int ,		/* ss_id1 - side set number for elementlist */
       const int ,		/* ss_id2 - side set number for elementlist */
       const int ,		/* ss_id3 - side set number for elementlist */
       int *,			/* num_elem - length of element list        */
       int **,			/* elems - ptr to element list              */
       int *,			/* ss_ptr - ss_index in EXODUS II db        */
       const int ,		/* eq - number of vector equation type      */
       const int ,		/* irc - counter of rotation conditions     */
       const Exo_DB *);	/* exo - the mesh */

EXTERN int set_pointers_to_vectors
(struct Rotation_Vectors **, /* normal                                 */
       struct Rotation_Vectors **, /* tangent1                               */
       struct Rotation_Vectors **, /* tangent2                               */
       struct Rotation_Vectors **, /* normal2                                */
       struct Rotation_Vectors **, /* normal3                                */
       struct Rotation_Vectors **, /* line_tangent                           */
       struct Rotation_Vectors **, /* binormal                               */
       struct Rotation_Vectors *, /* dum_vect                                */
       struct Rotation_Vectors *, /* vector                                  */
       const int ,		/* I - global node number                    */
       const int ,		/* eq - number of vector equation type       */
       const int ,		/* dim - dimensions                          */
       const int );		/* irc - counter for rotation Specs          */


EXTERN void rotate_eqns_at_node_2D
( int ,
		int ,
		int ,
		struct Aztec_Linear_Solver_System *);

#ifdef STATIC

/*
 * Prototype declarations of static functions in bc_rotate.c...
 */

static void rotate_res_jac_mesh
(const int ,		/* irow_index - Elemental stiffness row      */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for the current element */
       const int ,		/* nodes_per_elem - number of nodes in the 
				 * element                                   */
       const int ,		/* ielem_surf_dim - physical dim of element 
				 * surface - 0,1,2                           */
       double [DIM],		/* snormal - components of surface normal    */
       double [DIM][DIM][MDE],	/* dsnormal_dx - Array of derivatives 
				 * ([i][j][k]) of surface normal: component-i 
				 * wrt displacement-j at node-k              */
       double [2][DIM],		/* stangent - vector components of 2 mutually 
				 * orthogonal surface tangent vectors        */
       double [2][DIM][DIM][MDE], /* dstangent_dx - Array of derivatives 
				   * ([i][j][k]) of surface tangent vector 
				   * 1 or 2: component-i wrt displacement-j
				   * at node-k                               */
       const double sign);	/* sign of tangent vector                    */

static void rotate_res_jac_mom
(const int ,		/* irow_index - Elemental stiffness row      */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for the current element */
       const int ,		/* nodes_per_elem - number of nodes in the 
				 * element                                   */
       const int ,		/* ielem_surf_dim - physical dim of element 
				 * surface - 0,1,2                           */
       double [],		/* snormal - vector components of surface 
				 * normal                                    */
       double [][DIM][MDE],	/* dsnormal_dx - Array of derivatives 
				 * ([i][j][k]) of surface normal: 
				 * component-i wrt displacement-j at node-k  */
       double [][DIM],		/* stangent - vector components of 2 mutually 
				 * orthogonal surface tangent vectors        */
       double [][DIM][DIM][MDE], /* dstangent_dx - Array of derivatives 
				  * ([i][j][k]) of surface tangent vector 1 or
				  * 2: component-i wrt displacement-j at 
				  * node-k                                   */
       const double sign);	/* sign of tangent vector                    */


static void load_element_indices
(const Exo_DB *,		/* exo - ptr to FE db                        */
       const int ,		/* mn - material number                      */
       const int ,		/* ielem - element index number              */
       int *,			/* iconnect_ptr                              */
       int *,			/* ielem_type                                */
       int *,			/* ip_total                                  */
       int *,			/* num_local_nodes                           */
       int *,			/* ielem_dim                                 */
       int *);			/* dim                                       */

#endif

void
rotate_momentum_auto (
    int id,                             /* Elemental stiffness matrix row index */
    int I,                              /* Global node number                   */
    int dim,                            /* physical dim of problem              */
    struct Aztec_Linear_Solver_System *ams );

#endif /* GOMA_BC_ROTATE_H */
