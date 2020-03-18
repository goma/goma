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

#ifndef GOMA_MM_FILL_UTIL_H
#define GOMA_MM_FILL_UTIL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_UTIL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_UTIL_C
#define EXTERN extern
#endif

#include "el_elm.h"
#include "mm_as_structs.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_fill.h"
#include "std.h"

struct Aztec_Linear_Solver_System;

EXTERN int beer_belly(void);

EXTERN void
calc_surf_tangent(const int,    /* ielem - current element number            */
                  const int,    /* iconnect_ptr - Ptr into the beginning of the
                                 * connectivity list for the current element */
                  const int,    /* nodes_per_elem -                          */
                  const int,    /* ielem_surf_dim - physical dimension of the
                                 * surface of the element (0, 1, 2)          */
                  const int,    /* num_nodes_on_side - number of nodes on this
                                 * side of the element                       */
                  const int[]); /* local_elem_node_id - vector of the local
                                 * element node numbers on this side of the
                                 * element                                   */

EXTERN void calc_tangent_from_seed(struct Rotation_Vectors *, /* tangent */
                                   struct Rotation_Vectors *, /* normal */
                                   const double[DIM],         /* seed         */
                                   const int); /* dim                */

EXTERN void calc_tangent_along_basis(
    struct Rotation_Vectors *, /* tangent                                 */
    struct Rotation_Vectors *, /* normal                                  */
    const int,                 /* dim                                       */
    const int,                 /* id_side                                   */
    const int,                 /* num_nodes_on_side                         */
    const int[MAX_NODES_PER_SIDE]); /* local_elem_node_id               */

EXTERN void cross_vectors(struct Rotation_Vectors *, /* tangent2 */
                          struct Rotation_Vectors *, /* tangent1 */
                          struct Rotation_Vectors *, /* normal */
                          const int);                /* dim                */

EXTERN void calc_unseeded_edge_tangents(
    struct Rotation_Vectors *, /* tangent                                 */
    const int,                 /* iconnect_ptr - Ptr to beginning of the
                                * elem-node connectivity list for the current
                                * element                                   */
    const int,                 /* dim - physical dimension of the surface
                                * of the element i.e., (0, 1, 2)            */
    const int,                 /* id_side - ID of the side of the element   */
    const int,                 /* id_edge - ID of the edge of the element   */
    const int,                 /* num_nodes_on_edge - number of nodes on the
                                * edge of the element                       */
    const int[],               /* edge_elem_node_id - local element node
                                * numbers on the edge of the element        */
    const int);                /* param_dir - parametric direction of the
                                * edge in local coordinates                 */

EXTERN void calc_unseeded_edge_tangents_TET(
    struct Rotation_Vectors *, /* tangent                                 */
    const int,                 /* iconnect_ptr - Ptr to beginning of the
                                * elem-node connectivity list for the current
                                * element                                   */
    const int,                 /* dim - physical dimension of the surface
                                * of the element i.e., (0, 1, 2)            */
    const int,                 /* id_side - ID of the side of the element   */
    const int,                 /* id_edge - ID of the edge of the element   */
    const int,                 /* num_nodes_on_edge - number of nodes on the
                                * edge of the element                       */
    const int[],               /* edge_elem_node_id - local element node
                                * numbers on the edge of the element        */
    const int);                /* param_dir - parametric direction of the
                                * edge in local coordinates                 */

EXTERN void simple_normalize_vector(struct Rotation_Vectors *, /* vector */
                                    const int); /* dim                */

EXTERN int load_bf_grad(void);

EXTERN int load_bf_mesh_derivs(void);

EXTERN int load_basis_functions(
    const double[],             /*  xi - local element coordinates [DIM]     */
    struct Basis_Functions **); /* bfa - pointer to basis function       */

EXTERN void asdv(double **,  /* v - vector to be allocated */
                 const int); /* n - number of elements in vector */

EXTERN void dmemset(double *,     /* v - vector to be allocated */
                    const double, /* value to inserted into memory */
                    int);         /* n - number of elements in vector */

EXTERN void alloc_extern_ija_buffer(
    const int,   /* n - Order full system                (in) */
    const int,   /* m - Order trim system                (in) */
    const int *, /* ija - Orig column pointers           (in) */
    int **);     /* ptr_ija_attic - where to hide       (out) */

EXTERN void
alloc_sparse_arrays(int **,    /* ptr_ija - column pointer array            */
                    double **, /* ptr_a - values of nonzero matrix entries  */
                    double **, /* ptr_a_old - backup array, matrix nonzeros */
                    const int, /* Fill - flag to allocate spaces for either
                                * the totally coupled problem or just
                                * explicit nodal fill eqn                   */
                    int[],     /* node_to_fill - node index gives fill dof  */
                    Exo_DB *,  /* exo - ptr to the whole mesh               */
                    Dpi *);    /* dpi - distributed processing info         */

EXTERN void alloc_MSR_sparse_arrays /* mm_fill_util.c */
    (int **,                        /* ija - column pointer array */
     double **,                     /* a - values nonzero matrix entries */
     double **,                     /* a_old - backup for modified newton */
     int,                           /* Fill - flag to allocate space for
                                     * either the totally coupled problem
                                     * or just  explicit nodal fill eqn */
     int[],                         /* node_to_fill - node index gives
                                     * fill dof */
     Exo_DB *,                      /* exo - ptr to the whole mesh */
     Dpi *);                        /* dpi - distributed processing info */

EXTERN void alloc_VBR_sparse_arrays(struct Aztec_Linear_Solver_System *,
                                    Exo_DB *, /* ptr to the whole mesh */
                                    Dpi *);   /* distributed processing info */

EXTERN void zero_lec_row /* mm_fill_util.c                            */
    (double[][MAX_PROB_VAR + MAX_CONC][MDE][MDE], /* local_J             */
     int,  /* eqn_type - Eqn Type of row to be zeroed   */
     int); /* ldof - Local dof of that equation type    */

EXTERN void zero_lec_column                       /* mm_fill_util.c */
    (double[][MAX_PROB_VAR + MAX_CONC][MDE][MDE], /* local_J */
     int,  /* var_type - Variable type to be zeroed */
     int); /* ldof - Local dof of that variable */

EXTERN int find_VBR_index(const int, /* Block row index */
                          const int, /* Block column index */
                          struct Aztec_Linear_Solver_System *);

EXTERN double newshape(const double[], /* xi - local element coordinates */
                       const int,      /* Ielem_type - element type      */
                       const int,      /* Iquant - desired quantity (phi, phi_s,
                                        * etc.      */
                       const int,      /* Inode - current element node      */
                       const int,      /* eshape - element shape      */
                       const int,      /* interpolation - interpolation      */
                       const int); /* ledof - Typically, this is just the local
                                    * node number, but for pressure basis
                                    * functions, this can be a particular dof at
                                    * the centroid node */
EXTERN double
extended_shape(const double[], /* xi - local element coordinates            */
               const int,      /* Ielem_type - element type                 */
               const int,      /* Iquant - desired quantity (phi, phi_s,
                                * etc.                                      */
               const int,      /* Inode - current element node              */
               const int,      /* eshape - element shape                    */
               const int,      /* interpolation - interpolation             */
               const int);     /* ledof - Typically, this is just the local
                                * node number, but for pressure basis
                                * functions, this can be a particular dof at
                                * the centroid node                         */

EXTERN int
calc_shearrate(dbl *,          /* gammadot - strain rate scalar invariant   */
               dbl[DIM][DIM],  /* gamma_dot - strain rate tensor            */
               dbl[DIM][MDE],  /* d_gd_dv - deriv of gammadot wrt velocity  */
               dbl[DIM][MDE]); /* d_gd_dmesh - deriv of gammadot wrt meshd  */

EXTERN void
element_neighbor_list(const int, /* num_sides - of the elements (homogeneous) */
                      const int, /* num_total_nodes                           */
                      int **);   /* neighbors                                 */

EXTERN void shell_determinant_and_normal(
    const int,  /* ielem - current element number            */
    const int,  /* iconnect_ptr - Pointer to beginning of
                 * connectivity list for current element     */
    const int,  /* nodes_per_elem - number of nodes in elem  */
    const int,  /* ielem_surf_dim - physical dimension of the
                 * element surface (0, 1, 2)                 */
    const int); /* id_side - ID of element side              */

EXTERN double calc_tensor_invariant(dbl[DIM][DIM], // Original tensor
                                    dbl[DIM][DIM], // Sensitivities
                                    int); // Which invariant to calculate

extern void determine_ShapeVar(PROBLEM_DESCRIPTION_STRUCT *);
extern void determine_ProjectionVar(PROBLEM_DESCRIPTION_STRUCT *);
extern void set_solid_inertia(void);

extern int fill_variable_vector(int inode, int ivec_varType[],
                                int ivec_matID[]);

EXTERN void supg_tau_shakib(SUPG_terms *supg_terms, int dim, double dt,
                            int interp_eqn);

EXTERN void get_supg_tau(SUPG_terms *supg_terms,
                         int dim,
                         dbl diffusivity,
                         PG_DATA *pg_data);

EXTERN void supg_tau_gauss_point(SUPG_terms *supg_terms, int dim,
                                 dbl diffusivity, const PG_DATA *pg_data);

EXTERN void supg_tau(SUPG_terms *supg_terms, int dim, dbl diffusivity,
                     const PG_DATA *pg_data, double dt, int shakib, int interp_eqn);


EXTERN dbl yzbeta(dbl scale, int dim, dbl Y, dbl Z, dbl d_Z[MDE], dbl beta,
                       dbl u, dbl d_u[MDE], dbl grad_u[DIM],
                       dbl d_grad_u[MDE][DIM], dbl h_elem, int interp_eqn,
                       dbl deriv[MDE]);

EXTERN dbl yzbeta_model(int model, dbl scale, dbl beta, int dim,
                         dbl Y, dbl Z, dbl d_Z[MDE], dbl u,
                         dbl d_u[MDE], dbl grad_u[DIM],
                         dbl d_grad_u[MDE][DIM], dbl h_elem,
                         int interp_eqn, dbl deriv[MDE]);


#endif /* GOMA_MM_FILL_UTIL_H */
