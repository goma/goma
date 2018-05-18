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


#ifndef _DG_UTILS_H
#define _DG_UTILS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _DG_UTILS_C
#define EXTERN
#endif

#ifndef _DG_UTILS_C
#define EXTERN extern
#endif

#include "mm_as_structs.h"
#include "rf_node_const.h"
#include "mm_mp_structs.h"

typedef struct {
  int *neighbor_send_elem_index;
  int *neighbor_recv_elem_index;

  int *local_elems_to_send;

  double *nodal_send_mass_fraction_values;
  double *nodal_send_stress_values;

  int *global_recv_elem_list;

  int mass_fraction_num_dof;
  double *nodal_recv_mass_fraction_values;

  int stress_num_dof;
  double *nodal_recv_stress_values;

  double **nodal_recv_coord;

} dg_neighbor_type;


EXTERN int has_discontinuous_interp(PROBLEM_DESCRIPTION_STRUCT *pd, int var, int imtrx);

EXTERN int parallel_discontinuous_galerkin_enabled(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs, UPD_STRUCT *upd);

EXTERN void count_neighbor_needed_elems(Exo_DB *exo, Dpi *dpi, int *neighbor_elem_count);

EXTERN int proc_to_neighbor(int proc, Dpi *dpi);

EXTERN int all_elements_of_same_type(Dpi *dpi);

// returns true if interp is unset as well
EXTERN int parallel_all_interp_same(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                 UPD_STRUCT *upd,
                 int var);

EXTERN int setup_dg_neighbor_data(PROBLEM_DESCRIPTION_STRUCT **pd_ptrs,
                                  UPD_STRUCT *upd,
                                  Exo_DB *exo,
                                  Dpi *dpi,
                                  dg_neighbor_type *dg_neighbor_data);

EXTERN void dg_alloc_and_fill_coord(Exo_DB *exo, Dpi *dpi, dg_neighbor_type *neighbor_data);

EXTERN
void dg_communicate_neighbor_data(Exo_DB *exo,
                                  Dpi *dpi,
                                  dg_neighbor_type *neighbor_data,
                                  UPD_STRUCT *upd,
                                  double *x,
                                  int imtrx,
                                  NODE_INFO_STRUCT **Nodes,
                                  struct Viscoelastic_Nonmodal **vn_ptrs);

#endif
