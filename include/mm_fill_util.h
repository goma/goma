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
 

#ifndef _MM_FILL_UTIL_H
#define _MM_FILL_UTIL_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_UTIL_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_UTIL_C
#define EXTERN extern
#endif

EXTERN int beer_belly
PROTO((void));

EXTERN void calc_surf_tangent
PROTO((const int ,		/* ielem - current element number            */
       const int ,		/* iconnect_ptr - Ptr into the beginning of the
				 * connectivity list for the current element */
       const int ,		/* nodes_per_elem -                          */
       const int ,		/* ielem_surf_dim - physical dimension of the
				 * surface of the element (0, 1, 2)          */
       const int ,		/* num_nodes_on_side - number of nodes on this
				 * side of the element                       */
       const int []));		/* local_elem_node_id - vector of the local 
				 * element node numbers on this side of the 
				 * element                                   */

EXTERN void calc_tangent_from_seed
PROTO((struct Rotation_Vectors *, /* tangent                                 */
       struct Rotation_Vectors *, /* normal                                  */
       const double [DIM],	/* seed                                      */
       const int ));		/* dim                                       */

EXTERN void calc_tangent_along_basis
PROTO((struct Rotation_Vectors *, /* tangent                                 */
       struct Rotation_Vectors *, /* normal                                  */
       const int,		/* dim                                       */
       const int,		/* id_side                                   */
       const int,		/* num_nodes_on_side                         */
       const int [MAX_NODES_PER_SIDE])); /* local_elem_node_id               */

EXTERN void cross_vectors
PROTO((struct Rotation_Vectors *, /* tangent2                                */
       struct Rotation_Vectors *, /* tangent1                                */
       struct Rotation_Vectors *, /* normal                                  */
       const int ));		/* dim                                       */

EXTERN void calc_unseeded_edge_tangents
PROTO((struct Rotation_Vectors *, /* tangent                                 */
       const int ,		/* iconnect_ptr - Ptr to beginning of the
				 * elem-node connectivity list for the current
				 * element                                   */
       const int ,		/* dim - physical dimension of the surface 
				 * of the element i.e., (0, 1, 2)            */
       const int ,		/* id_side - ID of the side of the element   */
       const int ,		/* id_edge - ID of the edge of the element   */
       const int ,		/* num_nodes_on_edge - number of nodes on the 
				 * edge of the element                       */
       const int  [],		/* edge_elem_node_id - local element node 
				 * numbers on the edge of the element        */
       const int ));		/* param_dir - parametric direction of the 
				 * edge in local coordinates                 */

EXTERN void calc_unseeded_edge_tangents_TET
PROTO((struct Rotation_Vectors *, /* tangent                                 */
       const int ,		/* iconnect_ptr - Ptr to beginning of the
				 * elem-node connectivity list for the current
				 * element                                   */
       const int ,		/* dim - physical dimension of the surface 
				 * of the element i.e., (0, 1, 2)            */
       const int ,		/* id_side - ID of the side of the element   */
       const int ,		/* id_edge - ID of the edge of the element   */
       const int ,		/* num_nodes_on_edge - number of nodes on the 
				 * edge of the element                       */
       const int  [],		/* edge_elem_node_id - local element node 
				 * numbers on the edge of the element        */
       const int ));		/* param_dir - parametric direction of the 
				 * edge in local coordinates                 */
     
EXTERN void simple_normalize_vector
PROTO((struct Rotation_Vectors *, /* vector                                  */
       const int ));		/* dim                                       */

EXTERN int load_bf_grad
PROTO((void));

EXTERN int load_bf_mesh_derivs
PROTO((void));

EXTERN int load_basis_functions
PROTO((const double [],		/*  xi - local element coordinates [DIM]     */
       struct Basis_Functions **)); /* bfa - pointer to basis function       */

EXTERN void asdv
PROTO((double **,		/* v - vector to be allocated */
       const int ));		/* n - number of elements in vector */

EXTERN void dmemset
PROTO((double *,		/* v - vector to be allocated */
       const double,            /* value to inserted into memory */
       int ));	         	/* n - number of elements in vector */

EXTERN void alloc_extern_ija_buffer
PROTO((const int ,		/* n - Order full system                (in) */
       const int ,		/* m - Order trim system                (in) */
       const int *,		/* ija - Orig column pointers           (in) */
       int **));		/* ptr_ija_attic - where to hide       (out) */

EXTERN void alloc_sparse_arrays
PROTO((int **,			/* ptr_ija - column pointer array            */
       double **,		/* ptr_a - values of nonzero matrix entries  */
       double **,		/* ptr_a_old - backup array, matrix nonzeros */
       const int ,		/* Fill - flag to allocate spaces for either 
				 * the totally coupled problem or just
				 * explicit nodal fill eqn                   */
       int [],			/* node_to_fill - node index gives fill dof  */
       Exo_DB *,		/* exo - ptr to the whole mesh               */
       Dpi *));			/* dpi - distributed processing info         */

EXTERN void alloc_MSR_sparse_arrays		/* mm_fill_util.c */
PROTO((int    **,			/* ija - column pointer array */
       double **,			/* a - values nonzero matrix entries */
       double **,			/* a_old - backup for modified newton */
       int    ,				/* Fill - flag to allocate space for 
					 * either the totally coupled problem 
					 * or just  explicit nodal fill eqn */
       int [],				/* node_to_fill - node index gives 
					 * fill dof */
       Exo_DB *,			/* exo - ptr to the whole mesh */
       Dpi *  )); 			   /* dpi - distributed processing info */
                    
EXTERN void alloc_VBR_sparse_arrays     
PROTO((  struct Aztec_Linear_Solver_System *,
	 Exo_DB *,	                /* ptr to the whole mesh */
	 Dpi *));                       /* distributed processing info */

EXTERN void zero_lec_row	/* mm_fill_util.c                            */
PROTO((double [][MAX_PROB_VAR + MAX_CONC] [MDE][MDE], /* local_J             */
       int ,			/* eqn_type - Eqn Type of row to be zeroed   */
       int ));			/* ldof - Local dof of that equation type    */

EXTERN void zero_lec_column	/* mm_fill_util.c */
PROTO((double [][MAX_PROB_VAR + MAX_CONC][MDE][MDE], /* local_J */
       int ,			/* var_type - Variable type to be zeroed */
       int ));			/* ldof - Local dof of that variable */

EXTERN int find_VBR_index
PROTO((const int ,	/* Block row index */
       const int ,     /* Block column index */
       struct Aztec_Linear_Solver_System *) );

EXTERN double newshape
PROTO((const double [],		/* xi - local element coordinates            */
       const int ,		/* Ielem_type - element type                 */
       const int , 		/* Iquant - desired quantity (phi, phi_s, 
				 * etc.                                      */
       const int ,		/* Inode - current element node              */
       const int ,		/* eshape - element shape                    */
       const int ,		/* interpolation - interpolation             */
       const int ));		/* ledof - Typically, this is just the local
				 * node number, but for pressure basis 
				 * functions, this can be a particular dof at 
				 * the centroid node                         */
EXTERN double extended_shape
PROTO((const double [],		/* xi - local element coordinates            */
       const int ,		/* Ielem_type - element type                 */
       const int , 		/* Iquant - desired quantity (phi, phi_s,
				 * etc.                                      */
       const int ,		/* Inode - current element node              */
       const int ,		/* eshape - element shape                    */
       const int ,		/* interpolation - interpolation             */
       const int ));		/* ledof - Typically, this is just the local
				 * node number, but for pressure basis
				 * functions, this can be a particular dof at
				 * the centroid node                         */

EXTERN int calc_shearrate
PROTO((dbl *,			/* gammadot - strain rate scalar invariant   */
       dbl [DIM][DIM],		/* gamma_dot - strain rate tensor            */
       dbl [DIM][MDE],		/* d_gd_dv - deriv of gammadot wrt velocity  */
       dbl [DIM][MDE]));	/* d_gd_dmesh - deriv of gammadot wrt meshd  */

EXTERN void element_neighbor_list
PROTO((const int ,		/* num_sides - of the elements (homogeneous) */
       const int ,		/* num_total_nodes                           */
       int **));		/* neighbors                                 */

EXTERN void shell_determinant_and_normal
PROTO((const int ,		/* ielem - current element number            */
       const int ,		/* iconnect_ptr - Pointer to beginning of 
				 * connectivity list for current element     */
       const int ,		/* nodes_per_elem - number of nodes in elem  */
       const int ,		/* ielem_surf_dim - physical dimension of the 
				 * element surface (0, 1, 2)                 */
       const int ));		/* id_side - ID of element side              */

EXTERN double calc_tensor_invariant
PROTO(( dbl [DIM][DIM],        // Original tensor
	dbl [DIM][DIM],        // Sensitivities
	int ));                // Which invariant to calculate

extern void determine_ShapeVar(PROBLEM_DESCRIPTION_STRUCT *);
extern void determine_ProjectionVar(PROBLEM_DESCRIPTION_STRUCT *);
extern void set_solid_inertia(void);

extern int fill_variable_vector(int inode, int ivec_varType[], int ivec_matID[]);


#endif /* _MM_FILL_UTIL_H */
