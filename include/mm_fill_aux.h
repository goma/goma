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

/*
 * mm_fill_aux.h -- prototype declarations for mm_fill_aux.c
 */

#ifndef GOMA_MM_FILL_AUX_H
#define GOMA_MM_FILL_AUX_H

#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_fill.h"
#include "mm_fill_util.h"
#include "std.h"

struct Element_Indices;
struct Element_Stiffness_Pointers;
struct Field_Variables;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_AUX_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_AUX_C
#define EXTERN extern
#endif

EXTERN int load_coordinate_scales(const int, /* Coordinate system.                        */
                                  struct Field_Variables *); /* f - Field_Variable pointer. */

EXTERN int p0_vbl_dofptrs(struct Element_Indices *, /* lei - copy of the element indeces        */
                          struct Element_Stiffness_Pointers *); /* lsp                         */

EXTERN double global_velocity_norm(const dbl[],    /* x - solution vector                       */
                                   const Exo_DB *, /* exo - ptr to FE db                        */
                                   const Dpi *);   /* dpi - ptr to dist. info                   */

EXTERN dbl element_viscosity(void);

EXTERN void element_velocity(dbl[DIM],        /* v_avg                                     */
                             dbl[DIM][MDE],   /* dv_dnode                                  */
                             const Exo_DB *); /* exo - ptr to FE db                        */

EXTERN void h_elem_siz(dbl[DIM],      /* h                                         */
                       dbl[DIM][DIM], /* hh                                        */
                       dbl[DIM][MDE], /* dh_dxnode                                 */
                       const int DeformingMesh);

EXTERN void get_supg_stuff(dbl *,              /* supg_term                                 */
                           dbl[DIM],           /* elemental centroid velocity               */
                           dbl[DIM][MDE][DIM], /* d_vcent_du                                */
                           dbl[MDE][DIM],      /* d_supg_term_du                            */
                           dbl[MDE][DIM],      /* d_supg_term_dx                            */
                           const int);         /* DeformingMesh                             */

EXTERN void const_h_elem_siz(dbl[DIM],       /* h                                         */
                             dbl[DIM][DIM],  /* hh                                        */
                             dbl[DIM][MDE]); /* dh_dxnode                                 */

EXTERN dbl global_h_elem_siz(dbl[],    /* x                                         */
                             dbl[],    /* x_old                                     */
                             dbl[],    /* xdot                                      */
                             dbl[],    /* resid_vector                              */
                             Exo_DB *, /* exo - ptr to EXODUS II FE database        */
                             Dpi *);   /* dpi - distributed processing information  */

EXTERN void surface_determinant_and_normal(const int, /* ielem - current element number */
                                           const int, /* iconnect_ptr - Pointer to beginning of
                                                       * connectivity list for current element */
                                           const int, /* nodes_per_elem - number of nodes in elem */
                                           const int, /* ielem_surf_dim - physical dimension of the
                                                       * element surface (0, 1, 2) */
                                           const int, /* id_side - ID of element side */
                                           const int, /* num_nodes_on_side - number of nodes on
                                                       * side of element */
                                           const int[]); /* local_elem_node_id - local element node
                                                          * numbers on the side of the element */

EXTERN void edge_determinant_and_vectors(const int, /* ielem - current element number            */
                                         const int, /* iconnect_ptr - Pointer into the beginning
                                                     * of the connectivity list for the current
                                                     * element			             */
                                         const int, /* nodes_per_elem - number of nodes elem     */
                                         const int, /* ielem_surf_dim - physical dimension of the
                                                     * surface of the element (0, 1, 2)          */
                                         const int, /* id_side - ID of the side of the element   */
                                         const int, /* num_nodes_on_side - number of nodes on the
                                                     * side of the element */
                                         const int[], /* local_elem_node_id - local element node
                                                       * numbers on the side of the element */
                                         const int, /* id_edge - ID of the side of the element   */
                                         const int, /* num_nodes_on_edge - number of nodes on the
                                                     * edge of the element */
                                         const int[], /* edge_elem_node_id - vector of the local
                                                       * element node numbers on the side of the
                                                       * element */
                                         const int);  /* param_dir - local coordinate which
                                                       * parameterizes edge                        */

EXTERN void calc_CL_normal(double snormal[DIM],
                           double dsnormal_dx[DIM][DIM][MDE],
                           double fsnormal[DIM],
                           double dfsnormal_dx[DIM][DIM][MDE],
                           double tangent[DIM],
                           double dtangent_dx[DIM][DIM][MDE],
                           int elem,
                           int edge_elem_node_id[],
                           int elem_dim,
                           int num_edge_nodes,
                           double clnormal[DIM],
                           double dclnormal_dx[DIM][DIM][MDE],
                           const Exo_DB *exo);

#endif /* GOMA_MM_FILL_AUX_H */
