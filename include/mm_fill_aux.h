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

/*                                                          
 * mm_fill_aux.h -- prototype declarations for mm_fill_aux.c
 */
 

#ifndef _MM_FILL_AUX_H
#define _MM_FILL_AUX_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_AUX_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_AUX_C
#define EXTERN extern
#endif


EXTERN int load_coordinate_scales
PROTO((const int ,		/* Coordinate system.                        */
       struct Field_Variables *)); /* f - Field_Variable pointer.            */

EXTERN int p0_vbl_dofptrs
PROTO((struct Element_Indices *, /* lei - copy of the element indeces        */
       struct Element_Stiffness_Pointers *)); /* lsp                         */

EXTERN double global_velocity_norm
PROTO((const dbl [],		/* x - solution vector                       */
       const Exo_DB *,		/* exo - ptr to FE db                        */
       const Dpi *));		/* dpi - ptr to dist. info                   */

EXTERN dbl element_viscosity
PROTO((void));

EXTERN void element_velocity
PROTO((dbl [DIM],		/* v_avg                                     */
       dbl [DIM][MDE],		/* dv_dnode                                  */
       const Exo_DB *));	/* exo - ptr to FE db                        */

EXTERN void h_elem_siz
PROTO((dbl [DIM],		/* h                                         */
       dbl [DIM][DIM],		/* hh                                        */
       dbl [DIM][MDE],		/* dh_dxnode                                 */
       const int DeformingMesh));

EXTERN void get_supg_stuff
PROTO((dbl *,                    /* supg_term                                 */
       dbl [DIM],               /* elemental centroid velocity               */
       dbl [DIM][MDE][DIM],     /* d_vcent_du                                */
       dbl [MDE][DIM],          /* d_supg_term_du                            */
       dbl [MDE][DIM],          /* d_supg_term_dx                            */
       const int ));            /* DeformingMesh                             */
      
EXTERN void const_h_elem_siz
PROTO((dbl [DIM],		/* h                                         */
       dbl [DIM][DIM],		/* hh                                        */
       dbl [DIM][MDE]));	/* dh_dxnode                                 */

EXTERN dbl global_h_elem_siz
PROTO((dbl [],			/* x                                         */
       dbl [],			/* x_old                                     */
       dbl [],			/* xdot                                      */
       dbl [],			/* resid_vector                              */
       Exo_DB *,		/* exo - ptr to EXODUS II FE database        */
       Dpi *));			/* dpi - distributed processing information  */

EXTERN void surface_determinant_and_normal
PROTO((const int ,		/* ielem - current element number            */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for current element     */
       const int ,		/* nodes_per_elem - number of nodes in elem  */
       const int ,		/* ielem_surf_dim - physical dimension of the 
				 * element surface (0, 1, 2)                 */
       const int ,		/* id_side - ID of element side              */
       const int ,		/* num_nodes_on_side - number of nodes on 
				 * side of element                           */
       const int  [] ));	/* local_elem_node_id - local element node 
				 * numbers on the side of the element        */

EXTERN void edge_determinant_and_vectors
PROTO((const int ,		/* ielem - current element number            */
       const int ,		/* iconnect_ptr - Pointer into the beginning 
				 * of the connectivity list for the current
				 * element			             */
       const int ,		/* nodes_per_elem - number of nodes elem     */
       const int ,		/* ielem_surf_dim - physical dimension of the
				 * surface of the element (0, 1, 2)          */
       const int ,		/* id_side - ID of the side of the element   */
       const int ,		/* num_nodes_on_side - number of nodes on the 
				 * side of the element			     */
       const int  [],		/* local_elem_node_id - local element node 
				 * numbers on the side of the element        */
       const int ,		/* id_edge - ID of the side of the element   */
       const int ,		/* num_nodes_on_edge - number of nodes on the 
				 * edge of the element	                     */
       const int [],		/* edge_elem_node_id - vector of the local 
				 * element node numbers on the side of the 
				 * element                                   */
       const int ));		/* param_dir - local coordinate which 
				 * parameterizes edge                        */

EXTERN void calc_CL_normal
PROTO(( double [DIM],
	double [DIM][DIM][MDE],
	double [DIM],
	double [DIM][DIM][MDE],
	double [DIM],
	double [DIM][DIM][MDE],
	int,
	int [],
	int,
	int,
	double [DIM],
	double [DIM][DIM][MDE],

	const Exo_DB *));

#endif /* _MM_FILL_AUX_H */
