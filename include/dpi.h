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

#ifndef GOMA_DPI_H
#define GOMA_DPI_H

#include <stdbool.h>

struct Distributed_Processing_Information {
  int num_elem_blocks_global;
  int num_elems; /* New! */
  int num_neighbors;
  int num_node_sets_global;
  int num_side_sets_global;

  /*
   * variables (arrays, or scalars not used as dimensions).
   */

  int *eb_id_global;                 /* [neb_global] */
  int *eb_num_nodes_per_elem_global; /* New! - [neb_global] */

  int *elem_index_global; /* New! [num_elems_proc] */

  /*
   * Some new stuff to help in the assembly of element based methods, such
   * as discontinuous Galerkin.
   */

  int *elem_owner; /* [num_elems_proc] Which proc owns this e? */
  /* These are referenced using the pointers
   * from exo->elem_elem_pntr[elem] */

  int *elem_elem_list_global; /* Face-ordered names of elems facing this. */

  int *neighbor;          /* Array of neighboring processor id's to this
                             proc. Values range from 0 to Nproc - 1
                             length = [num_neighbors] */
  int *node_index_global; /* Array of Global node numbers.
                             The index is the processor node number
                             length - [num_nodes_proc] */
  int *ns_id_global;      /* [num_node_sets_global] */

  int *ss_internal_global;

  int num_dofs_global;    /* scalar
                           * Total number of unknowns in the global
                           * problem */
  int num_elems_global;   /* scalar
                           * Total number of elements in the global
                           * problem */
  int num_internal_nodes; /* scalar
                           * Number of nodes for which this processor
                           * has primary responsibility but which
                           * never need to be communicated to/from
                           * neighboring processors.  */
  int num_boundary_nodes; /* scalar - number of nodes for which this
                           * processor has primary responsibility but
                           * which are associated with unknowns that
                           * need to be sent out to other neighboring
                           * processors (external from their view).
                           */
  int num_external_nodes; /* scalar
                           * Number of nodes for which this processor
                           * has secondary responsibility. Such nodes
                           * are associated with values that need
                           * to be gathered from other processors and
                           * are used in calculations.
                           */
  int num_owned_nodes;    /* scalar - (calculated)
                           *  num_internal_nodes + num_boundary_nodes */
  int num_universe_nodes; /* This is equal to the number of internal plus
                           * boundary plus external nodes on
                           * the processor */
  int num_nodes_global;   /* The total number of nodes in the global
                           * problem
                           *   (scalar) */
  int *ss_id_global;      /* [num_side_sets_global] */

  int *ss_index_global; /* [num_side_sets] */

  // New Nemesis
  int *num_ns_global_node_counts;
  int *num_ns_global_df_counts;

  int *num_ss_global_side_counts;
  int *num_ss_global_df_counts;

  int *global_elem_block_ids;
  int *global_elem_block_counts;

  int num_proc;
  int num_proc_in_file;
  char ftype;

  int num_node_cmaps;
  int num_elem_cmaps;

  int num_internal_elems;
  int num_border_elems;

  int *proc_node_internal;
  int *proc_node_boundary;
  int *proc_node_external;

  int *proc_elem_internal;
  int *proc_elem_border;

  int *node_cmap_ids;
  int *elem_cmap_ids;
  int *node_cmap_node_counts;
  int *elem_cmap_elem_counts;

  int **node_map_node_ids;
  int **node_map_proc_ids;

  int **elem_cmap_elem_ids;
  int **elem_cmap_proc_ids;
  int **elem_cmap_side_ids;

  int *ss_block_index_global;
  int *ss_block_list_global;

  // old dpi
  int *node_owner;
  int *num_node_recv;
  int *num_node_send;

  // omega_h
  int *exodus_to_omega_h_node;

  // base mesh dpi info
  int base_internal_nodes;
  int base_boundary_nodes;
  int base_external_nodes;

  int base_internal_elems;
  int base_border_elems;

  // ns and ss consistency for Cubit
  bool goma_dpi_data;
  int global_ns_node_len;
  int *global_ns_nodes;
  int global_ss_elem_len;
  int *global_ss_elems;
  int *global_ss_sides;
};
typedef struct Distributed_Processing_Information Dpi;

extern Dpi *DPI_ptr;
#endif

// NETCDF definitions

#define GOMA_NC_DIM_LEN_NS_NODE_LIST "goma_len_ns_node_list"
#define GOMA_NC_DIM_LEN_SS_ELEM_LIST "goma_len_ss_elem_list"
#define GOMA_NC_VAR_NS_NODE_LIST     "goma_ns_node_list"
#define GOMA_NC_VAR_SS_ELEM_LIST     "goma_ss_elem_list"
#define GOMA_NC_VAR_SS_SIDE_LIST     "goma_ss_side_list"