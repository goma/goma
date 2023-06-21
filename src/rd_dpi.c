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

/* rd_dpi.c -- routines for reading distributed processing information
 * Uses Nemesis structure from SEACAS
 */

#define GOMA_RD_DPI_C

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <el_elm_info.h>
#include <string.h>

#include "dp_ghost.h"
#include "dpi.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "mm_eh.h"
#include "rd_dpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_mp.h"
#include "std.h"

// Helper for exodus return values
#define CHECK_EX_ERROR(err, format, ...)                              \
  do {                                                                \
    if (err < 0) {                                                    \
      goma_eh(GOMA_ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__); \
    }                                                                 \
  } while (0)

int rd_dpi(Exo_DB *exo, Dpi *d, char *fn, bool parallel_call) {
  init_dpi_struct(d);
  float version = -4.98; /* initialize. ex_open() changes this. */
  int comp_wordsize = sizeof(dbl);
  int io_wordsize = 0;
  int exoid = ex_open(fn, EX_READ, &comp_wordsize, &io_wordsize, &version);
  CHECK_EX_ERROR(exoid, "ex_open");
  int ex_error;

  ex_error = ex_get_init_global(exoid, &d->num_nodes_global, &d->num_elems_global,
                                &d->num_elem_blocks_global, &d->num_node_sets_global,
                                &d->num_side_sets_global);

  CHECK_EX_ERROR(ex_error, "ex_get_init_global");

  // Node Set Global
  d->ns_id_global = alloc_int_1(d->num_node_sets_global, 0);
  d->num_ns_global_node_counts = alloc_int_1(d->num_node_sets_global, 0);
  d->num_ns_global_df_counts = alloc_int_1(d->num_node_sets_global, 0);
  ex_error = ex_get_ns_param_global(exoid, d->ns_id_global, d->num_ns_global_node_counts,
                                    d->num_ns_global_df_counts);
  CHECK_EX_ERROR(ex_error, "ex_get_ns_param_global");
  // Side Set Global
  d->ss_id_global = alloc_int_1(d->num_side_sets_global, 0);
  d->num_ss_global_side_counts = alloc_int_1(d->num_side_sets_global, 0);
  d->num_ss_global_df_counts = alloc_int_1(d->num_side_sets_global, 0);
  ex_error = ex_get_ss_param_global(exoid, d->ss_id_global, d->num_ss_global_side_counts,
                                    d->num_ss_global_df_counts);

  CHECK_EX_ERROR(ex_error, "ex_get_ss_param_global");
  // Block Global
  d->global_elem_block_ids = alloc_int_1(d->num_elem_blocks_global, 0);
  d->global_elem_block_counts = alloc_int_1(d->num_elem_blocks_global, 0);
  ex_error = ex_get_eb_info_global(exoid, d->global_elem_block_ids, d->global_elem_block_counts);
  CHECK_EX_ERROR(ex_error, "ex_get_eb_info_global");

  // Nemesis Info
  d->rank = ProcID;
  ex_error = ex_get_init_info(exoid, &d->num_proc, &d->num_proc_in_file, &d->ftype);
  CHECK_EX_ERROR(ex_error, "ex_get_init_info");
  if (d->num_proc != Num_Proc) {
    GOMA_EH(GOMA_ERROR, "Nemesis mesh error num_proc != number of mpi processes");
  }
  if (d->num_proc_in_file != 1) {
    GOMA_EH(GOMA_ERROR, "Nemesis mesh error expected num_proc_in_file == 1");
  }
  if (d->num_proc != Num_Proc) {
    GOMA_EH(GOMA_ERROR, "Nemesis mesh error num_proc != number of mpi processes");
  }
  if (d->ftype != 'p') {
    GOMA_EH(GOMA_ERROR, "Nemesis mesh error ftype expected 'p'");
  }
  // global indices
  d->node_index_global = alloc_int_1(exo->num_nodes, 0);
  d->elem_index_global = alloc_int_1(exo->num_elems, 0);
  ex_error = ex_get_id_map(exoid, EX_NODE_MAP, d->node_index_global);
  CHECK_EX_ERROR(ex_error, "ex_get_id_map EX_NODE_MAP");
  ex_error = ex_get_id_map(exoid, EX_ELEM_MAP, d->elem_index_global);
  CHECK_EX_ERROR(ex_error, "ex_get_id_map EX_ELEM_MAP");

  // set base mesh  global indices
  exo->base_mesh->node_map = alloc_int_1(exo->num_nodes, 0);
  exo->base_mesh->elem_map = alloc_int_1(exo->num_elems, 0);
  memcpy(exo->base_mesh->node_map, d->node_index_global, sizeof(int) * exo->num_nodes);
  memcpy(exo->base_mesh->elem_map, d->elem_index_global, sizeof(int) * exo->num_elems);

  // Load Balance Information
  ex_error = ex_get_loadbal_param(
      exoid, &d->num_internal_nodes, &d->num_boundary_nodes, &d->num_external_nodes,
      &d->num_internal_elems, &d->num_border_elems, &d->num_node_cmaps, &d->num_elem_cmaps, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_loadbal_param");

  d->base_internal_nodes = d->num_internal_nodes;
  d->base_boundary_nodes = d->num_boundary_nodes;
  d->base_external_nodes = d->num_external_nodes;
  d->base_internal_elems = d->num_internal_elems;
  d->base_border_elems = d->num_border_elems;

  d->proc_node_internal = alloc_int_1(d->num_internal_nodes, 0);
  if (d->num_boundary_nodes > 0) {
    d->proc_node_boundary = alloc_int_1(d->num_boundary_nodes, 0);
  }
  if (d->num_external_nodes > 0) {
    d->proc_node_external = alloc_int_1(d->num_external_nodes, 0);
  }

  ex_error = ex_get_processor_node_maps(exoid, d->proc_node_internal, d->proc_node_boundary,
                                        d->proc_node_external, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_processor_node_maps");

  d->node_cmap_ids = alloc_int_1(d->num_node_cmaps, 0);
  d->node_cmap_node_counts = alloc_int_1(d->num_node_cmaps, 0);
  if (d->num_elem_cmaps > 0) {
    d->elem_cmap_ids = alloc_int_1(d->num_elem_cmaps, 0);
    d->elem_cmap_elem_counts = alloc_int_1(d->num_elem_cmaps, 0);
  }

  ex_error = ex_get_cmap_params(exoid, d->node_cmap_ids, d->node_cmap_node_counts, d->elem_cmap_ids,
                                d->elem_cmap_elem_counts, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_cmake_params");

  d->node_map_node_ids = calloc(sizeof(int *), d->num_node_cmaps);
  d->node_map_proc_ids = calloc(sizeof(int *), d->num_node_cmaps);

  for (int i = 0; i < d->num_node_cmaps; i++) {
    d->node_map_node_ids[i] = alloc_int_1(d->node_cmap_node_counts[i], 0);
    d->node_map_proc_ids[i] = alloc_int_1(d->node_cmap_node_counts[i], 0);
    ex_error = ex_get_node_cmap(exoid, d->node_cmap_ids[i], d->node_map_node_ids[i],
                                d->node_map_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_get_node_cmap %d", i);
  }

  if (d->num_elem_cmaps > 0) {
    d->elem_cmap_elem_ids = calloc(sizeof(int *), d->num_elem_cmaps);
    d->elem_cmap_side_ids = calloc(sizeof(int *), d->num_elem_cmaps);
    d->elem_cmap_proc_ids = calloc(sizeof(int *), d->num_elem_cmaps);
  }

  for (int i = 0; i < d->num_elem_cmaps; i++) {
    d->elem_cmap_elem_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    d->elem_cmap_side_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    d->elem_cmap_proc_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    ex_error = ex_get_elem_cmap(exoid, d->elem_cmap_ids[i], d->elem_cmap_elem_ids[i],
                                d->elem_cmap_side_ids[i], d->elem_cmap_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_get_elem_cmap %d", i);
  }

  // Setup old dpi information

  d->eb_id_global = calloc(d->num_elem_blocks_global, sizeof(int));
  int *eb_num_nodes_local = calloc(sizeof(int), d->num_elem_blocks_global);
  d->eb_num_nodes_per_elem_global = calloc(sizeof(int), d->num_elem_blocks_global);
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    for (int j = 0; j < d->num_elem_blocks_global; j++) {
      if (d->global_elem_block_ids[j] == exo->eb_id[i]) {
        eb_num_nodes_local[j] = exo->eb_num_nodes_per_elem[i];
        d->eb_id_global[i] = d->global_elem_block_ids[j];
      }
    }
  }

  d->ss_index_global = malloc(sizeof(int) * exo->num_side_sets);
  for (int i = 0; i < exo->num_side_sets; i++) {
    for (int j = 0; j < d->num_side_sets_global; j++) {
      if (exo->ss_id[i] == d->ss_id_global[j]) {
        d->ss_index_global[i] = j;
      }
    }
  }

  int *local_ss_internal = NULL;
  int *global_ss_internal = NULL;
  d->ss_internal_global = NULL;

  if (d->num_side_sets_global > 0) {
    if (parallel_call) {
      local_ss_internal = find_ss_internal_boundary(exo);
      global_ss_internal = calloc(d->num_side_sets_global, sizeof(int));
      d->ss_internal_global = calloc(d->num_side_sets_global, sizeof(int));

      for (int i = 0; i < exo->num_side_sets; i++) {
        int global_ss_index = d->ss_index_global[i];
        global_ss_internal[global_ss_index] = local_ss_internal[i];
      }

      d->ss_block_index_global = calloc(d->num_side_sets_global + 1, sizeof(int));
      d->ss_block_list_global =
          calloc(d->num_elem_blocks_global * d->num_side_sets_global, sizeof(int));
      int *ss_block_count_proc =
          calloc(d->num_elem_blocks_global * d->num_side_sets_global, sizeof(int));
      int *ss_block_count_global =
          calloc(d->num_elem_blocks_global * d->num_side_sets_global, sizeof(int));

      for (int ss_id = 0; ss_id < exo->num_side_sets; ss_id++) {
        for (int elem_index = exo->ss_elem_index[ss_id];
             elem_index < (exo->ss_elem_index[ss_id] + exo->ss_num_sides[ss_id]); elem_index++) {
          int elem = exo->ss_elem_list[elem_index];
          int block = find_elemblock_index(elem, exo);
          if (block == -1) {
            GOMA_EH(GOMA_ERROR, "Element block not found ss_block_list");
          }
          for (int j = 0; j < d->num_elem_blocks_global; j++) {
            if (d->eb_id_global[j] == exo->eb_id[block]) {
              ss_block_count_proc[ss_id * d->num_elem_blocks_global + j] = 1;
            }
          }
        }
      }

      MPI_Allreduce(ss_block_count_proc, ss_block_count_global,
                    d->num_elem_blocks_global * d->num_side_sets_global, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD);

      int ss_block_index = 0;
      for (int ss_id = 0; ss_id < d->num_side_sets_global; ss_id++) {
        int ss_start = ss_block_index;
        for (int j = 0; j < d->num_elem_blocks_global; j++) {
          int offset = ss_id * d->num_elem_blocks_global + j;
          if (ss_block_count_global[offset] > 0) {
            d->ss_block_list_global[ss_block_index] = j;
            ss_block_index++;
          }
        }
        d->ss_block_index_global[ss_id] = ss_start;
        d->ss_block_index_global[ss_id + 1] = ss_block_index;
      }

      free(ss_block_count_proc);
      free(ss_block_count_global);
    }
  }
  d->elem_owner = alloc_int_1(exo->num_elems, ProcID);

  if (parallel_call) {
    const int num_mpi_async = 2;
    MPI_Request request_array[2];

    MPI_Iallreduce(eb_num_nodes_local, d->eb_num_nodes_per_elem_global, d->num_elem_blocks_global,
                   MPI_INT, MPI_MAX, MPI_COMM_WORLD, &(request_array[0]));

    MPI_Iallreduce(global_ss_internal, d->ss_internal_global, d->num_side_sets_global, MPI_INT,
                   MPI_MAX, MPI_COMM_WORLD, &(request_array[1]));

    MPI_Waitall(num_mpi_async, request_array, MPI_STATUSES_IGNORE);

    free(local_ss_internal);
    free(global_ss_internal);

    int min_external;
    MPI_Allreduce(&d->num_external_nodes, &min_external, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (min_external > 0) {
      GOMA_EH(-1, "Found > 0 external nodes, use element decomposition");
    }
  }
  free(eb_num_nodes_local);
  ex_error = ex_close(exoid);
  CHECK_EX_ERROR(ex_error, "ex_close");

  d->num_universe_nodes = d->num_internal_nodes + d->num_boundary_nodes + d->num_external_nodes;

  // zero_base maps
  for (int i = 0; i < exo->base_mesh->num_nodes; i++) {
    exo->base_mesh->node_map[i]--;
  }
  for (int i = 0; i < exo->base_mesh->num_elems; i++) {
    exo->base_mesh->elem_map[i]--;
  }
  zero_dpi(d);

  // setup node owners
  if (parallel_call) {
    d->num_neighbors = d->num_node_cmaps;
    d->neighbor = alloc_int_1(d->num_neighbors, -1);
    d->node_owner = alloc_int_1(exo->num_nodes, ProcID);
    d->num_node_recv = alloc_int_1(d->num_neighbors, 0);
    d->num_node_send = alloc_int_1(d->num_neighbors, 0);
    for (int i = 0; i < d->num_node_cmaps; i++) {
      int neighbor = d->node_map_proc_ids[i][0];
      for (int j = 0; j < d->node_cmap_node_counts[i]; j++) {
        if (neighbor != d->node_map_proc_ids[i][j]) {
          GOMA_EH(GOMA_ERROR, "Unexpected proc id in node map proc ids");
        }
        if (neighbor < d->node_owner[d->node_map_node_ids[i][j]]) {
          d->node_owner[d->node_map_node_ids[i][j]] = neighbor;
          d->num_node_recv[i]++;
        }
      }
      d->neighbor[i] = neighbor;
    }

    int elem_offset = 0;
    for (int i = 0; i < exo->num_elem_blocks; i++) {
      for (int j = 0; j < exo->eb_num_elems[i]; j++) {
        int min_proc = ProcID;
        int nnode_per_elem = exo->eb_num_nodes_per_elem[i];
        for (int k = 0; k < nnode_per_elem; k++) {
          int local_node = exo->eb_conn[i][j * nnode_per_elem + k];
          int proc = d->node_owner[local_node];
          if (proc < min_proc) {
            min_proc = proc;
          }
        }
        d->elem_owner[elem_offset] = min_proc;
        elem_offset++;
      }
    }

    goma_error err = generate_ghost_elems(exo, d);
    GOMA_EH(err, "generate_ghost_elements");

    int *num_send_nodes = alloc_int_1(d->num_neighbors, 0);
    int **global_send_nodes = malloc(sizeof(int *) * d->num_neighbors);
    int *num_recv_nodes = alloc_int_1(d->num_neighbors, 0);
    int **global_recv_nodes = malloc(sizeof(int *) * d->num_neighbors);

    MPI_Request *requests = calloc(sizeof(MPI_Request), 2 * d->num_neighbors);

    for (int i = 0; i < d->num_neighbors; i++) {
      MPI_Irecv(&num_send_nodes[i], 1, MPI_INT, d->neighbor[i], 206, MPI_COMM_WORLD,
                &requests[d->num_neighbors + i]);
      // printf("Proc %d recv from Proc %d tag %d", d->neighbor[i], ProcID, 206);
    }
    for (int i = 0; i < d->num_neighbors; i++) {
      num_recv_nodes[i] = 0;
      for (int j = 0; j < d->num_external_nodes; j++) {
        if (d->node_owner[d->num_internal_nodes + d->num_boundary_nodes + j] == d->neighbor[i]) {
          num_recv_nodes[i] += 1;
        }
      }
      MPI_Isend(&num_recv_nodes[i], 1, MPI_INT, d->neighbor[i], 206, MPI_COMM_WORLD, &requests[i]);
      // printf("Proc %d sending to Proc %d tag %d", ProcID, d->neighbor[i], 206);
    }

    MPI_Waitall(d->num_neighbors * 2, requests, MPI_STATUSES_IGNORE);

    for (int i = 0; i < d->num_neighbors; i++) {
      global_send_nodes[i] = malloc(sizeof(int) * num_send_nodes[i]);
      MPI_Irecv(global_send_nodes[i], num_send_nodes[i], MPI_INT, d->neighbor[i], 207,
                MPI_COMM_WORLD, &requests[d->num_neighbors + i]);
      // printf("Proc %d recv from Proc %d tag %d", d->neighbor[i], ProcID, 207);
    }

    for (int i = 0; i < d->num_neighbors; i++) {
      global_recv_nodes[i] = malloc(sizeof(int) * num_recv_nodes[i]);
      int recv_index = 0;
      for (int j = 0; j < d->num_external_nodes; j++) {
        int local_node = d->num_internal_nodes + d->num_boundary_nodes + j;
        if (d->node_owner[local_node] == d->neighbor[i]) {
          global_recv_nodes[i][recv_index++] = d->node_index_global[local_node];
        }
      }
      MPI_Isend(global_recv_nodes[i], num_recv_nodes[i], MPI_INT, d->neighbor[i], 207,
                MPI_COMM_WORLD, &requests[i]);
      // printf("Proc %d sending to Proc %d tag %d", ProcID, d->neighbor[i], 207);
    }

    MPI_Waitall(d->num_neighbors * 2, requests, MPI_STATUSES_IGNORE);
    free(requests);

    // verify we have the send nodes
    for (int i = 0; i < d->num_neighbors; i++) {
      for (int j = 0; j < num_send_nodes[i]; j++) {
        int exists = in_list(global_send_nodes[i][j], d->num_internal_nodes,
                             d->num_internal_nodes + d->num_boundary_nodes, d->node_index_global);
        if (exists == -1) {
          GOMA_EH(GOMA_ERROR, "Required node to communicate doesn't exist in border nodes");
        }
      }
    }

#ifdef DEBUG_MPI

    int *global_owners_min = alloc_int_1(d->num_nodes_global, Num_Proc + 1);
    int *global_owners_max = alloc_int_1(d->num_nodes_global, -1);
    int *all_global_owners_min = alloc_int_1(d->num_nodes_global, Num_Proc + 1);
    int *all_global_owners_max = alloc_int_1(d->num_nodes_global, -1);
    for (int i = 0; i < exo->num_nodes; i++) {
      int gindex = d->node_index_global[i];
      int owner = d->node_owner[i];
      global_owners_max[gindex] = owner;
      global_owners_min[gindex] = owner;
    }

    MPI_Allreduce(global_owners_min, all_global_owners_min, d->num_nodes_global, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(global_owners_max, all_global_owners_max, d->num_nodes_global, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    for (int i = 0; i < d->num_nodes_global; i++) {
      if (all_global_owners_min[i] != all_global_owners_max[i]) {
        GOMA_EH(GOMA_ERROR, "Inconsistent node owners");
      }
    }
    free(global_owners_min);
    free(global_owners_max);
    free(all_global_owners_max);
    free(all_global_owners_min);

#endif

    // reorder external nodes
    int *new_external_node_order = alloc_int_1(d->num_external_nodes, 0);
    int *old_to_new_external_node_order = alloc_int_1(d->num_external_nodes, 0);
    int *old_node_owner = alloc_int_1(d->num_external_nodes, 0);
    int new_index = 0;
    int offset = d->num_internal_nodes + d->num_boundary_nodes;
    for (int j = 0; j < d->num_neighbors; j++) {
      int neighbor = d->neighbor[j];
      for (int i = 0; i < d->num_external_nodes; i++) {
        old_node_owner[i] = d->node_owner[i + offset];
        if (d->node_owner[d->num_internal_nodes + d->num_boundary_nodes + i] == neighbor) {
          new_external_node_order[new_index] = i;
          old_to_new_external_node_order[i] = new_index;
          new_index++;
        }
      }
    }
    if (new_index != d->num_external_nodes) {
      GOMA_EH(GOMA_ERROR, "incorrect new ordering %d != %d", new_index, d->num_external_nodes);
    }

    // node index global
    int *old_global_indices = alloc_int_1(d->num_external_nodes, 0);
    for (int i = 0; i < d->num_external_nodes; i++) {
      old_global_indices[i] =
          d->node_index_global[d->num_internal_nodes + d->num_boundary_nodes + i];
    }
    for (int i = 0; i < d->num_external_nodes; i++) {
      d->node_index_global[d->num_internal_nodes + d->num_boundary_nodes + i] =
          old_global_indices[new_external_node_order[i]];
    }
    free(old_global_indices);
    double *x_old = alloc_dbl_1(d->num_external_nodes, 0);
    double *y_old = alloc_dbl_1(d->num_external_nodes, 0);
    double *z_old = alloc_dbl_1(d->num_external_nodes, 0);
    for (int i = 0; i < d->num_external_nodes; i++) {
      x_old[i] = exo->x_coord[offset + i];
      if (exo->num_dim > 1) {
        y_old[i] = exo->y_coord[offset + i];
      }
      if (exo->num_dim > 2) {
        z_old[i] = exo->z_coord[offset + i];
      }
    }
    for (int i = 0; i < d->num_external_nodes; i++) {
      exo->x_coord[offset + i] = x_old[new_external_node_order[i]];
      if (exo->num_dim > 1) {
        exo->y_coord[offset + i] = y_old[new_external_node_order[i]];
      }
      if (exo->num_dim > 2) {
        exo->z_coord[offset + i] = z_old[new_external_node_order[i]];
      }
    }
    free(x_old);
    free(y_old);
    free(z_old);

    // conn
    for (int block = 0; block < exo->num_elem_blocks; block++) {
      int num_nodes_per_blk = exo->eb_num_nodes_per_elem[block] * exo->eb_num_elems[block];
      for (int i = 0; i < num_nodes_per_blk; i++) {
        if (exo->eb_conn[block][i] >= offset) {
          exo->eb_conn[block][i] =
              old_to_new_external_node_order[exo->eb_conn[block][i] - offset] + offset;
        }
      }
    }

    // ns
    for (int j = 0; j < exo->ns_node_len; j++) {
      if (exo->ns_node_list[j] >= offset) {
        exo->ns_node_list[j] =
            old_to_new_external_node_order[exo->ns_node_list[j] - offset] + offset;
      }
    }

    // ss
    for (int ins = 0; ins < exo->num_side_sets; ins++) {
      for (int side_index = 0; side_index < exo->ss_num_sides[ins]; side_index++) {
        for (int lni = exo->ss_node_side_index[ins][side_index];
             lni < exo->ss_node_side_index[ins][side_index + 1]; lni++) {
          int inode = exo->ss_node_list[ins][lni];
          if (inode >= offset) {
            exo->ss_node_list[ins][lni] = old_to_new_external_node_order[inode - offset] + offset;
          }
        }
      }
    }

    // update node owners
    for (int i = 0; i < d->num_external_nodes; i++) {
      d->node_owner[offset + old_to_new_external_node_order[i]] = old_node_owner[i];
    }

    // setup indexing from ghosted mesh to base mesh
    setup_ghost_to_base(exo, d);
    d->num_owned_nodes = d->num_internal_nodes + d->num_boundary_nodes;
    d->num_universe_nodes = d->num_internal_nodes + d->num_boundary_nodes + d->num_external_nodes;

    free(new_external_node_order);
    free(old_to_new_external_node_order);
    free(old_node_owner);
    free(num_send_nodes);
    free(num_recv_nodes);
    for (int i = 0; i < d->num_neighbors; i++) {
      free(global_send_nodes[i]);
      free(global_recv_nodes[i]);
    }
    free(global_send_nodes);
    free(global_recv_nodes);
  }

  d->goma_dpi_data = false;
  // read ns and ss consistency data
  int ncid;
  int err = nc_open(fn, NC_NOWRITE | NC_SHARE, &ncid);
  if (err)
    GOMA_EH(GOMA_ERROR, nc_strerror(err));

  int nc_ns_id;
  size_t nc_ns_len;
  bool goma_ns_found = true;
  err = nc_inq_dimid(ncid, GOMA_NC_DIM_LEN_NS_NODE_LIST, &nc_ns_id);
  if (err != NC_NOERR) {
    goma_ns_found = false;
  }

  int nc_ss_id;
  size_t nc_ss_len;
  bool goma_ss_found = true;
  err = nc_inq_dimid(ncid, GOMA_NC_DIM_LEN_SS_ELEM_LIST, &nc_ss_id);
  if (err != NC_NOERR) {
    goma_ss_found = false;
  }

  d->global_ns_node_len = 0;
  if (goma_ns_found) {
    err = nc_inq_dimlen(ncid, nc_ns_id, &nc_ns_len);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));

    int nc_node_list;
    err = nc_inq_varid(ncid, GOMA_NC_VAR_NS_NODE_LIST, &nc_node_list);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));

    d->goma_dpi_data = true;
    d->global_ns_node_len = nc_ns_len;
    d->global_ns_nodes = calloc(nc_ns_len, sizeof(int));
    err = nc_get_var(ncid, nc_node_list, d->global_ns_nodes);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));
  }

  d->global_ss_elem_len = 0;
  if (goma_ss_found) {
    err = nc_inq_dimlen(ncid, nc_ss_id, &nc_ss_len);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));

    int nc_elem_list;
    int nc_side_list;
    err = nc_inq_varid(ncid, GOMA_NC_VAR_SS_ELEM_LIST, &nc_elem_list);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));
    err = nc_inq_varid(ncid, GOMA_NC_VAR_SS_SIDE_LIST, &nc_side_list);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));

    d->goma_dpi_data = true;
    d->global_ss_elem_len = nc_ss_len;
    d->global_ss_elems = calloc(nc_ss_len, sizeof(int));
    d->global_ss_sides = calloc(nc_ss_len, sizeof(int));
    err = nc_get_var(ncid, nc_elem_list, d->global_ss_elems);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));
    err = nc_get_var(ncid, nc_side_list, d->global_ss_sides);
    if (err)
      GOMA_EH(GOMA_ERROR, nc_strerror(err));
  }

  err = nc_close(ncid);
  if (err)
    GOMA_EH(GOMA_ERROR, nc_strerror(err));

  return 0;
}
int zero_dpi(Dpi *d) {
  for (int i = 0; i < d->num_universe_nodes; i++) {
    d->node_index_global[i] -= 1;
  }

  for (int i = 0; i < (d->num_internal_elems + d->num_border_elems); i++) {
    d->elem_index_global[i] -= 1;
  }

  for (int cmap = 0; cmap < d->num_node_cmaps; cmap++) {
    for (int node = 0; node < d->node_cmap_node_counts[cmap]; node++) {
      d->node_map_node_ids[cmap][node] -= 1;
    }
  }

  return GOMA_SUCCESS;
}
int one_dpi(Dpi *d) {
  for (int i = 0; i < d->num_universe_nodes; i++) {
    d->node_index_global[i] += 1;
  }

  for (int i = 0; i < (d->num_internal_elems + d->num_border_elems); i++) {
    d->elem_index_global[i] += 1;
  }

  for (int cmap = 0; cmap < d->num_node_cmaps; cmap++) {
    for (int node = 0; node < d->node_cmap_node_counts[cmap]; node++) {
      d->node_map_node_ids[cmap][node] += 1;
    }
  }

  return GOMA_SUCCESS;
}
/* uni_dpi() -- setup distributed processing information for one processor
 *
 * When the problem is not partitioned, we need to set some acceptable
 * defaults.
 */
void uni_dpi(Dpi *dpi, Exo_DB *exo) {
  int i;
  int len;
  dpi->num_elems = exo->num_elems;

  if (exo->elem_elem_conn_exists) {
    dpi->elem_elem_list_global = exo->elem_elem_list;
  } else {
    dpi->elem_elem_list_global = NULL;
  }

  len = dpi->num_elems;
  dpi->elem_owner = alloc_int_1(len, ProcID);

  dpi->num_elems_global = exo->num_elems;

  dpi->num_elem_blocks_global = exo->num_elem_blocks;
  dpi->num_neighbors = 0;
  dpi->num_node_sets_global = exo->num_node_sets;
  dpi->num_side_sets_global = exo->num_side_sets;

  dpi->num_elems_global = exo->num_elems;
  dpi->num_internal_nodes = exo->num_nodes;
  dpi->num_boundary_nodes = 0;
  dpi->num_external_nodes = 0;
  dpi->num_owned_nodes = exo->num_nodes;
  dpi->num_universe_nodes = exo->num_nodes;
  dpi->num_nodes_global = exo->num_nodes;

  /*
   * Note! This aliasing of these pointers into the exo structure has
   * two advantages and one disadvantage.
   *
   *	(+) it's very easy to do
   *	(+) it makes more economic use of memory
   *	(-) it makes it too easy to free_dpi() and nuke your
   *        EXODUS information
   *
   * We'll just do it for now!
   */

  len = dpi->num_elems;
  dpi->elem_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < len; i++) {
    dpi->elem_index_global[i] = i;
  }

  len = dpi->num_universe_nodes;
  dpi->node_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < len; i++) {
    dpi->node_index_global[i] = i;
  }
  dpi->ss_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < exo->num_side_sets; i++) {
    dpi->ss_index_global[i] = i;
  }

  dpi->eb_id_global = exo->eb_id;

  dpi->eb_num_nodes_per_elem_global = exo->eb_num_nodes_per_elem;

  dpi->neighbor = alloc_int_1(1, ProcID);

  dpi->ns_id_global = exo->ns_id;

  dpi->ss_id_global = exo->ss_id;

  dpi->ss_internal_global = find_ss_internal_boundary(exo);

  // base mesh settings that can't be setup until dpi
  exo->base_mesh->node_map = malloc(sizeof(int) * exo->num_nodes);
  exo->base_mesh->elem_map = malloc(sizeof(int) * exo->num_elems);
  memcpy(exo->base_mesh->node_map, dpi->node_index_global, sizeof(int) * exo->num_nodes);
  memcpy(exo->base_mesh->elem_map, dpi->elem_index_global, sizeof(int) * exo->num_elems);

  exo->ghost_node_to_base = malloc(sizeof(int) * exo->num_nodes);
  for (int i = 0; i < exo->num_nodes; i++) {
    exo->ghost_node_to_base[i] = i;
  }
  exo->eb_ghost_elem_to_base = malloc(sizeof(int *) * exo->num_elem_blocks);
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    if (exo->eb_num_elems[i] > 0) {
      exo->eb_ghost_elem_to_base[i] = malloc(sizeof(int) * exo->eb_num_elems[i]);
      for (int j = 0; j < exo->eb_num_elems[i]; j++) {
        exo->eb_ghost_elem_to_base[i][j] = j;
      }
    } else {
      exo->eb_ghost_elem_to_base[i] = NULL;
    }
  }

  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
/* free_dpi() -- free internal dynamically allocated memory in a
 *               Dpi struct
 *
 * During rd_dpi(), various arrays are allocated. To cleanly
 * free up this  memory, this routine can be used. Typically,
 * if dpi is a (Dpi *), then
 *
 *	free_dpi(dpi);
 *	free(dpi);
 *
 * would accomplish the desired effect.
 *
 * Created: 1997/08/23 15:50 MDT pasacki@sandia.gov
 */

void free_dpi(Dpi *d) {
  free(d->ns_id_global);
  free(d->num_ns_global_node_counts);
  free(d->num_ns_global_df_counts);
  free(d->ss_id_global);
  free(d->num_ss_global_side_counts);
  free(d->num_ss_global_df_counts);
  free(d->global_elem_block_ids);
  free(d->global_elem_block_counts);
  free(d->node_index_global);
  free(d->elem_index_global);
  free(d->proc_elem_internal);
  if (d->num_border_elems > 0) {
    free(d->proc_elem_border);
  }
  free(d->proc_node_internal);
  if (d->num_boundary_nodes > 0) {
    free(d->proc_node_boundary);
  }
  if (d->num_external_nodes > 0) {
    free(d->proc_node_external);
  }
  free(d->node_cmap_ids);
  free(d->node_cmap_node_counts);
  if (d->num_elem_cmaps > 0) {
    free(d->elem_cmap_ids);
    free(d->elem_cmap_elem_counts);
  }

  for (int i = 0; i < d->num_node_cmaps; i++) {
    free(d->node_map_node_ids[i]);
    free(d->node_map_proc_ids[i]);
  }
  free(d->node_map_node_ids);
  free(d->node_map_proc_ids);

  for (int i = 0; i < d->num_elem_cmaps; i++) {
    free(d->elem_cmap_elem_ids[i]);
    free(d->elem_cmap_side_ids[i]);
    free(d->elem_cmap_proc_ids[i]);
  }

  if (d->num_elem_cmaps > 0) {
    free(d->elem_cmap_elem_ids);
    free(d->elem_cmap_side_ids);
    free(d->elem_cmap_proc_ids);
  }

  free(d->eb_id_global);
  free(d->eb_num_nodes_per_elem_global);
  free(d->ss_index_global);

  free(d->ss_internal_global);

  if (d->num_side_sets_global > 0) {
    free(d->ss_block_index_global);
    free(d->ss_block_list_global);
  }

  free(d->elem_owner);
  free(d->neighbor);
  free(d->node_owner);
  free(d->num_node_recv);
  free(d->num_node_send);
  free(d->exodus_to_omega_h_node);

  if (d->goma_dpi_data) {
    free(d->global_ns_nodes);
    free(d->global_ss_elems);
    free(d->global_ss_sides);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
/* free_dpi_uni() -- free internal dynamically allocated memory in
 *                   Dpi SERIAL!
 *
 * During rd_dpi(), various arrays are allocated. To cleanly free up this
 * memory, this routine can be used. Typically, if dpi is a (Dpi *), then
 *
 *	free_dpi(dpi);
 *	free(dpi);
 *
 * would accomplish the desired effect.
 *
 * !!!! IMPORTANT - DANGEROUS PROGRAMMING NOTE !!!!!
 * For serial processing, not all of the pieces of Dpi have been
 * allocated and some are merely aliases to parts of the exodus ii
 * database. To avoid munging that data, just free up what was
 * allocated in uni_dpi.
 *
 * Created: 1997/09/11 10:43 MDT pasacki@sandia.gov
 */

void free_dpi_uni(Dpi *d) {
  safer_free((void **)&(d->elem_owner));
  safer_free((void **)&(d->elem_index_global));
  safer_free((void **)&(d->node_index_global));
  safer_free((void **)&(d->neighbor));
  safer_free((void **)&(d->ss_index_global));
  free(d->ss_internal_global);

  return;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/* init_dpi_struct() -- initialize some defaults
 *
 * This is meant to be called right after allocation, to help setup some
 * reasonable defaults to describe an empty data structure. Call it a poor
 * man's constructor.
 *
 * Created: 1999/08/11 17:04 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void init_dpi_struct(Dpi *d) {
  if (d == NULL) {
    GOMA_EH(GOMA_ERROR, "Empty structure to initialize?");
  }
  memset((void *)d, 0, sizeof(Dpi));
  return;
}
/************************************************************************/
/************************************************************************/

/* exo_dpi_clone -- transfer needed global monolith data to child piece
 *
 *
 * Notes: Avert aliasing problems by allocating full-fledged arrays for dpi
 *        that won't disappear if the monolith does.
 *
 * Created: 1999/08/24 11:44 MDT pasacki@sandia.gov
 */

void exo_dpi_clone(Exo_DB *exo, Dpi *dpi) {
  int len;

  dpi->num_nodes_global = exo->num_nodes;
  dpi->num_elems_global = exo->num_elems;
  dpi->num_elem_blocks_global = exo->num_elem_blocks;
  dpi->num_node_sets_global = exo->num_node_sets;
  dpi->num_side_sets_global = exo->num_side_sets;

  /*
   * Allocate and fill arrays for element blocks...
   */

  len = dpi->num_elem_blocks_global * sizeof(int);

  dpi->eb_id_global = smalloc(len);
  memcpy(dpi->eb_id_global, exo->eb_id, len);

  /*
   * Allocate and fill arrays for node sets...
   */

  len = dpi->num_node_sets_global * sizeof(int);

  dpi->ns_id_global = smalloc(len);
  memcpy(dpi->ns_id_global, exo->ns_id, len);

  /*
   * Allocate and fill arrays for side sets...
   */
  len = dpi->num_side_sets_global * sizeof(int);

  dpi->ss_id_global = smalloc(len);
  memcpy(dpi->ss_id_global, exo->ss_id, len);

  return;
}
