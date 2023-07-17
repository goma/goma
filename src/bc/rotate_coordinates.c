#include "bc/rotate_coordinates.h"

#include <assert.h>
#include <mpi.h>
#include <rd_mesh.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bc/rotate_util.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "gds/gds_vector.h"
#include "load_field_variables.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#include "sl_util_structs.h"
#include "std.h"
#include "stdbool.h"

#ifndef GOMA_MAX_NORMALS_PER_NODE
#define GOMA_MAX_NORMALS_PER_NODE 50
#endif

#ifndef MAX_NUM_SIDESETS
#define MAX_NUM_SIDESETS 10000
#endif

goma_rotation_s goma_automatic_rotations = {false, NULL};

goma_error check_if_equation_is_rotation(int equation, bool *is_rotated) {
  switch (equation) {
  case R_MOM_NORMAL:
  case R_MOM_TANG1:
  case R_MOM_TANG2:
  case R_MESH_NORMAL:
  case R_MESH_TANG1:
  case R_MESH_TANG2:
    *is_rotated = true;
    break;
  default:
    *is_rotated = false;
  }
  if (!(*is_rotated) && (equation < V_FIRST || equation >= V_LAST)) {
    return GOMA_ERROR;
  }
  return GOMA_SUCCESS;
}

goma_error
setup_bc_is_rotated_list(struct Boundary_Condition *bc_types, int num_bc, bool **bc_rotate_list) {
  assert(num_bc > 0);
  *bc_rotate_list = calloc((size_t)num_bc, sizeof(bool));
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    (*bc_rotate_list)[bc_index] = false;
    if (!strcmp(bc_types[bc_index].Set_Type, "SS")) {
      bool bc_is_rotated = false;
      check_if_equation_is_rotation(BC_Desc[bc_types[bc_index].BC_Desc_index].equation,
                                    &bc_is_rotated);
      if (bc_is_rotated) {
        (*bc_rotate_list)[bc_index] = true;
      }
    }
  }
  return GOMA_SUCCESS;
}

goma_error allocate_rotations(Exo_DB *exo, goma_rotation_node_s **rotations) {
  if (goma_automatic_rotations.rotation_nodes != NULL) {
    for (int i = 0; i < exo->num_nodes; i++) {
      for (int j = 0; j < DIM; j++) {
        goma_normal_free((goma_automatic_rotations.rotation_nodes)[i].rotated_coord[j]);
      }
    }
    free(goma_automatic_rotations.rotation_nodes);
  }
  *rotations = calloc((size_t)exo->num_nodes, sizeof(goma_rotation_node_s));
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < DIM; j++) {
      (*rotations)[i].rotated_coord[j] = goma_normal_alloc(3);
    }
  }

  return GOMA_SUCCESS;
}
goma_error free_rotations(Exo_DB *exo, goma_rotation_node_s **rotations) {
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < DIM; j++) {
      free((*rotations)[i].rotated_coord[j]);
    }
  }
  free(*rotations);
  return GOMA_SUCCESS;
}

typedef struct {
  int ss;
  int num_sides;
  int total_nodes;
  int *node_per_side;
  int *global_node_ids;
} ss_edge_share;

typedef struct {
  int first;
  int second;
} int_pair;

static int int_pair_compare_first(const void *left, const void *right) {
  int_pair *left_val = ((int_pair *)left);
  int_pair *right_val = ((int_pair *)right);
  if (left_val->first > right_val->first) {
    return 1;
  } else if (left_val->first < right_val->first) {
    return -1;
  } else {
    return 0;
  }
}

static int int_pair_compare(const void *left, const void *right) {
  int_pair *left_val = ((int_pair *)left);
  int_pair *right_val = ((int_pair *)right);
  if (left_val->first > right_val->first) {
    return 1;
  } else if (left_val->first < right_val->first) {
    return -1;
  } else {
    if (left_val->second > right_val->second) {
      return 1;
    } else if (left_val->second < right_val->second) {
      return -1;
    } else {
      return 0;
    }
  }
}

goma_error exchange_neighbor_ss_edges(Exo_DB *exo, Dpi *dpi) {
  int *rotated_side_sets = calloc(exo->num_side_sets, sizeof(int));
  int num_rotated_side_sets = exo->num_side_sets;
  for (int i = 0; i < exo->num_side_sets; i++) {
    rotated_side_sets[i] = i;
  }

  ss_edge_share *ss_edge_info = calloc(dpi->num_side_sets_global, sizeof(ss_edge_share));
  int max_nodes_on_side = 0;
  for (int i = 0; i < num_rotated_side_sets; i++) {
    int ss_index = rotated_side_sets[i];
    ss_edge_info[i].ss = exo->ss_id[ss_index];
    ss_edge_info[i].num_sides = exo->ss_num_sides[ss_index];
    ss_edge_info[i].node_per_side = calloc(ss_edge_info[i].num_sides, sizeof(int));
    int total_nodes = 0;
    for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
      int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];
      int ielem_type = Elem_Type(exo, ielem);
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[ss_index] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
      if (num_nodes_on_side > max_nodes_on_side) {
        max_nodes_on_side = num_nodes_on_side;
      }
      ss_edge_info[i].node_per_side[e] = num_nodes_on_side;
      total_nodes += num_nodes_on_side;
    }
    ss_edge_info[i].global_node_ids = calloc(total_nodes, sizeof(int));
    ss_edge_info[i].total_nodes = total_nodes;
  }

  for (int i = 0; i < num_rotated_side_sets; i++) {
    int ss_index = rotated_side_sets[i];
    int offset = 0;
    for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
      int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

      int ielem_type = Elem_Type(exo, ielem);
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[ss_index] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
      int *elem_node_id = &(local_side_node_list[0]);
      /* use nodal points only!! */
      for (int k = 0; k < num_nodes_on_side; k++) {
        int id = elem_node_id[k];
        int I = Proc_Elem_Connect[Proc_Connect_Ptr[ielem] + id];
        int global_node = dpi->node_index_global[I];
        ss_edge_info[i].global_node_ids[offset++] = global_node;
      }
    }
    assert(offset == ss_edge_info[i].total_nodes);
  }

  // char *format = "gid%d.csv";
  // char fname[80];
  // snprintf(fname, 79, format, ProcID);
  // FILE *f = fopen(fname, "w");
  // fprintf(f,"x,y,z,gid,ss\n");
  // for (int k = 0; k < exo->num_side_sets; k++) {
  //   for (int i = 0; i < ss_edge_info[k].total_nodes; i++) {
  //     int ln = in_list(ss_edge_info[k].global_node_ids[i], 0, exo->num_nodes,
  //     dpi->node_index_global); fprintf(f, "%g,%g,%g,%d,%d\n",
  //         exo->x_coord[ln],
  //         exo->y_coord[ln],
  //         exo->z_coord[ln],
  //         ss_edge_info[k].global_node_ids[i],
  //         ss_edge_info[k].ss);
  //   }
  // }
  // fclose(f);
  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Finalize();
  // exit(0);

  // get num nodes

  MPI_Request *requests =
      malloc(sizeof(MPI_Request) * 2 * dpi->num_neighbors * dpi->num_side_sets_global);

  int *recv_nodes = calloc(dpi->num_neighbors * dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_neighbors; i++) {
    GOMA_CHECK_MPI_ERROR(MPI_Irecv(&recv_nodes[i * dpi->num_side_sets_global],
                                   dpi->num_side_sets_global, MPI_INT, dpi->neighbor[i],
                                   102 + dpi->neighbor[i], MPI_COMM_WORLD, &(requests[i])));
  }

  int *send_nodes = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int j = 0; j < dpi->num_side_sets_global; j++) {
    send_nodes[j] = ss_edge_info[j].total_nodes;
  }
  for (int i = 0; i < dpi->num_neighbors; i++) {
    GOMA_CHECK_MPI_ERROR(MPI_Send(send_nodes, dpi->num_side_sets_global, MPI_INT, dpi->neighbor[i],
                                  102 + ProcID, MPI_COMM_WORLD));
  }
  GOMA_CHECK_MPI_ERROR(MPI_Waitall(dpi->num_neighbors, requests, MPI_STATUSES_IGNORE));

  // get nodes

  int **global_nodes = calloc(dpi->num_neighbors * dpi->num_side_sets_global, sizeof(int *));
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    for (int j = 0; j < dpi->num_neighbors; j++) {
      int all_nodes = recv_nodes[j * dpi->num_side_sets_global + i];
      global_nodes[j * dpi->num_side_sets_global + i] = calloc(all_nodes, sizeof(int));
    }
  }

  // we have an odd warning with optimization silence with:
  GOMA_ASSERT_ALWAYS(dpi->num_side_sets_global < MAX_NUM_SIDESETS);
  int *offsets = calloc(dpi->num_side_sets_global, sizeof(int));
  int req_count = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      int count = recv_nodes[i * dpi->num_side_sets_global + j];
      if (count > 0) {
        GOMA_CHECK_MPI_ERROR(MPI_Irecv(global_nodes[i * dpi->num_side_sets_global + j], count,
                                       MPI_INT, dpi->neighbor[i], 103 + j + dpi->neighbor[i],
                                       MPI_COMM_WORLD, &(requests[req_count++])));
      }
      offsets[j] += count;
    }
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      if (ss_edge_info[j].total_nodes > 0) {
        GOMA_CHECK_MPI_ERROR(MPI_Send(ss_edge_info[j].global_node_ids, ss_edge_info[j].total_nodes,
                                      MPI_INT, dpi->neighbor[i], 103 + j + ProcID, MPI_COMM_WORLD));
      }
    }
  }

  GOMA_CHECK_MPI_ERROR(MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE));

  // get num sides
  int *recv_sides = calloc(dpi->num_neighbors * dpi->num_side_sets_global, sizeof(int));
  int *send_sides = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    send_sides[i] = ss_edge_info[i].num_sides;
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    GOMA_CHECK_MPI_ERROR(MPI_Irecv(&(recv_sides[i * dpi->num_side_sets_global]),
                                   dpi->num_side_sets_global, MPI_INT, dpi->neighbor[i],
                                   104 + dpi->neighbor[i], MPI_COMM_WORLD, &(requests[i])));
    // printf("-1 Proc %d IRECV %d fr %d size %d\n", ProcID, 104+ProcID, dpi->neighbor[i],
    // dpi->num_side_sets_global);
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    GOMA_CHECK_MPI_ERROR(MPI_Send(send_sides, dpi->num_side_sets_global, MPI_INT, dpi->neighbor[i],
                                  104 + ProcID, MPI_COMM_WORLD));
    // printf("-1 Proc %d ISEND %d to %d size %d\n", ProcID, 104+ProcID, dpi->neighbor[i],
    // dpi->num_side_sets_global);
  }
  GOMA_CHECK_MPI_ERROR(MPI_Waitall(dpi->num_neighbors, requests, MPI_STATUSES_IGNORE));

  MPI_Barrier(MPI_COMM_WORLD);

  int **recv_nodes_per_side = calloc(dpi->num_neighbors * dpi->num_side_sets_global, sizeof(int *));
  for (int i = 0; i < dpi->num_neighbors * dpi->num_side_sets_global; i++) {
    if (recv_sides[i] > 0) {
      recv_nodes_per_side[i] = calloc(recv_sides[i], sizeof(int));
    } else {
      recv_nodes_per_side[i] = NULL;
    }
  }

  req_count = 0;
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      int n_sides = recv_sides[i * dpi->num_side_sets_global + j];
      if (n_sides > 0) {
        GOMA_CHECK_MPI_ERROR(MPI_Irecv(recv_nodes_per_side[i * dpi->num_side_sets_global + j],
                                       n_sides, MPI_INT, dpi->neighbor[i],
                                       104 + j + dpi->neighbor[i], MPI_COMM_WORLD,
                                       &(requests[req_count++])));
        // printf("-2 Proc %d IRECV %d fr %d size %d\n", ProcID, 204 + j + dpi->neighbor[i],
        // dpi->neighbor[i], dpi->num_side_sets_global);
      }
    }
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      if (send_sides[j] > 0) {
        GOMA_CHECK_MPI_ERROR(MPI_Send(ss_edge_info[j].node_per_side, send_sides[j], MPI_INT,
                                      dpi->neighbor[i], 104 + j + ProcID, MPI_COMM_WORLD));
      }
      // printf("-2 Proc %d ISEND %d to %d size %d\n", ProcID, 204 + j + ProcID, dpi->neighbor[i],
      // dpi->num_side_sets_global);
    }
  }
  GOMA_CHECK_MPI_ERROR(MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE));

  int_pair *global_to_local = calloc(exo->num_nodes, sizeof(int_pair));
  for (int i = 0; i < exo->num_nodes; i++) {
    global_to_local[i].first = dpi->node_index_global[i];
    global_to_local[i].second = i;
  }
  qsort(global_to_local, exo->num_nodes, sizeof(int_pair), int_pair_compare);

  // convert all to local
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    for (int k = 0; k < dpi->num_neighbors; k++) {
      for (int j = 0; j < recv_nodes[i + k * dpi->num_side_sets_global]; j++) {
        int_pair global_node = {global_nodes[i + k * dpi->num_side_sets_global][j], -1};
        int_pair *match = (int_pair *)bsearch(&global_node, global_to_local, exo->num_nodes,
                                              sizeof(int_pair), int_pair_compare_first);
        if (match != NULL) {
          global_nodes[i + k * dpi->num_side_sets_global][j] = match->second;
        } else {
          global_nodes[i + k * dpi->num_side_sets_global][j] = -1;
        }
      }
    }
    for (int j = 0; j < ss_edge_info[i].total_nodes; j++) {
      int_pair global_node = {ss_edge_info[i].global_node_ids[j], -1};
      int_pair *match = (int_pair *)bsearch(&global_node, global_to_local, exo->num_nodes,
                                            sizeof(int_pair), int_pair_compare_first);
      if (match != NULL) {
        ss_edge_info[i].global_node_ids[j] = match->second;
      } else {
        GOMA_EH(GOMA_ERROR, "no mapping to local node for local ss node");
      }
    }
  }

  int_pair **ss_elem_sides_local = malloc(sizeof(int_pair *) * dpi->num_side_sets_global);
  int *ss_elem_sides_count_local = malloc(sizeof(int) * dpi->num_side_sets_global);

  // suppress warning I'm unsure about
  GOMA_ASSERT_ALWAYS((((size_t)dpi->num_side_sets_global) * sizeof(int)) < PTRDIFF_MAX);
  // end suppress
  int *ss_global_nodes = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      ss_global_nodes[j] += recv_nodes[i * dpi->num_side_sets_global + j];
    }
  }

  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    ss_elem_sides_count_local[i] = 0;
    if (ss_global_nodes[i] > 0) {
      ss_elem_sides_local[i] = calloc(sizeof(int_pair), ss_global_nodes[i]);
    } else {
      ss_elem_sides_local[i] = NULL;
    }
  }

  int_pair **known_ss_elem_sides_local = malloc(sizeof(int_pair *) * dpi->num_side_sets_global);
  int *known_ss_elem_sides_count = malloc(sizeof(int) * dpi->num_side_sets_global);
  for (int i = 0; i < exo->num_side_sets; i++) {
    known_ss_elem_sides_local[i] = calloc(exo->ss_num_sides[i], sizeof(int_pair));
    known_ss_elem_sides_count[i] = 0;
    for (int j = 0; j < exo->ss_num_sides[i]; j++) {
      int_pair pair = {exo->ss_elem_list[exo->ss_elem_index[i] + j],
                       exo->ss_side_list[exo->ss_elem_index[i] + j]};
      known_ss_elem_sides_local[i][j] = pair;
      known_ss_elem_sides_count[i] += 1;
    }
    qsort(known_ss_elem_sides_local[i], known_ss_elem_sides_count[i], sizeof(int_pair),
          int_pair_compare);
  }

  // loop over all sides and see if we find a matching side (start slow way)
  int other_side_nodes[MAX_NODES_PER_SIDE];
  for (int k = 0; k < dpi->num_neighbors; k++) {
    for (int i = 0; i < dpi->num_side_sets_global; i++) {
      int offset = 0;
      for (int side = 0; side < recv_sides[k * dpi->num_side_sets_global + i]; side++) {
        int n_nodes_other = recv_nodes_per_side[k * dpi->num_side_sets_global + i][side];
        // setup other side nodes
        bool skip = false;
        for (int n = 0; n < n_nodes_other; n++) {
          other_side_nodes[n] = global_nodes[k * dpi->num_side_sets_global + i][offset + n];
          if (other_side_nodes[n] == -1) {
            skip = true;
          }
        }
        offset += n_nodes_other;
        if (skip || n_nodes_other == 0) {
          continue;
        }

        // find matching elem side
        int node = other_side_nodes[0];
        for (int idx = exo->node_elem_pntr[node]; idx < exo->node_elem_pntr[node + 1]; idx++) {
          int elem = exo->node_elem_list[idx];
          int ielem_type = exo->eb_elem_itype[exo->elem_eb[elem]];
          int shape = type2shape(ielem_type);
          int n_sides = shape2sides(shape);
          for (int side = 0; side < n_sides; side++) {
            int local_side_node_list[MAX_NODES_PER_SIDE];
            int side_nodes[MAX_NODES_PER_SIDE];

            /* find SIDE info for primary side */
            int num_nodes_on_side = 0;
            get_side_info(ielem_type, side + 1, &num_nodes_on_side, local_side_node_list);

            /* use nodal points only!! */
            bool found = false;
            for (int k = 0; k < num_nodes_on_side; k++) {
              int id = local_side_node_list[k];
              int iconnect = Proc_Connect_Ptr[elem];
              int I = Proc_Elem_Connect[iconnect + id];
              side_nodes[k] = I;
              if (I == node) {
                found = true;
              }
            }

            if (found) {
              bool found_all = true;
              bool all_known = false;
              if (num_nodes_on_side != n_nodes_other) {
                found_all = false;
              }
              for (int k = 0; k < num_nodes_on_side; k++) {
                int sn = side_nodes[k];
                int found_ln = false;
                for (int n = 0; n < n_nodes_other; n++) {
                  if (sn == other_side_nodes[n]) {
                    found_ln = true;
                  }
                }
                if (!found_ln) {
                  found_all = false;
                }
              }
              int_pair elem_side = {elem, side};
              if (found_all) {
                int_pair *exists =
                    bsearch(&elem_side, known_ss_elem_sides_local[i], known_ss_elem_sides_count[i],
                            sizeof(int_pair), int_pair_compare);
                if (exists != NULL) {
                  all_known = true;
                }
              }

              if (found_all && !all_known) {
                if (ss_elem_sides_count_local[i] > 1) {
                  qsort(ss_elem_sides_local[i], ss_elem_sides_count_local[i], sizeof(int_pair),
                        int_pair_compare);
                }
                int_pair *exists =
                    bsearch(&elem_side, ss_elem_sides_local[i], ss_elem_sides_count_local[i],
                            sizeof(int_pair), int_pair_compare);
                if (exists == NULL) {
                  ss_elem_sides_local[i][ss_elem_sides_count_local[i]] = elem_side;
                  ss_elem_sides_count_local[i] += 1;
                }
              }
            }
          }
        }
      }
      assert(offset == recv_nodes[k * dpi->num_side_sets_global + i]);
    }
  }

  // fix exodus sidesets
  int *old_ss_elem_list = exo->ss_elem_list;
  int *old_ss_side_list = exo->ss_side_list;
  // find new sizes;
  int new_size = 0;
  for (int i = 0; i < exo->num_side_sets; i++) {
    int global_ss_index = -1;
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      if (dpi->ss_id_global[j] == exo->ss_id[i]) {
        global_ss_index = j;
        break;
      }
    }

    new_size += exo->ss_num_sides[i];
    if (ss_elem_sides_count_local[global_ss_index] > 0) {
      new_size += ss_elem_sides_count_local[global_ss_index];
    }
  }

  exo->ss_elem_list = calloc(new_size, sizeof(int));
  exo->ss_side_list = calloc(new_size, sizeof(int));
  // setup new lists
  int full_offset = 0;
  for (int i = 0; i < exo->num_side_sets; i++) {
    int global_ss_index = -1;
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      if (dpi->ss_id_global[j] == exo->ss_id[i]) {
        global_ss_index = j;
        break;
      }
    }

    for (int j = 0; j < exo->ss_num_sides[i]; j++) {
      exo->ss_elem_list[full_offset + j] = old_ss_elem_list[exo->ss_elem_index[i] + j];
      exo->ss_side_list[full_offset + j] = old_ss_side_list[exo->ss_elem_index[i] + j];
    }
    if (ss_elem_sides_count_local[global_ss_index] > 0) {
      int offset = exo->ss_num_sides[i];
      for (int j = 0; j < ss_elem_sides_count_local[global_ss_index]; j++) {
        exo->ss_elem_list[j + offset + full_offset] = ss_elem_sides_local[global_ss_index][j].first;
        exo->ss_side_list[j + offset + full_offset] =
            ss_elem_sides_local[global_ss_index][j].second + 1;
      }
    }
    exo->ss_elem_index[i] = full_offset;
    full_offset += exo->ss_num_sides[i] + ss_elem_sides_count_local[global_ss_index];
    exo->ss_num_sides[i] += ss_elem_sides_count_local[global_ss_index];
  }
  assert(full_offset == new_size);
  exo->ss_elem_len = full_offset;
  free(old_ss_elem_list);
  free(old_ss_side_list);

  for (int i = 0; i < exo->num_side_sets; i++) {
    free(exo->ss_node_cnt_list[i]);
    exo->ss_node_cnt_list[i] = (int *)malloc(exo->ss_num_sides[i] * sizeof(int));

    int ss_node_list_len = 0;
    for (int e = 0; e < exo->ss_num_sides[i]; e++) {
      int ielem = exo->ss_elem_list[exo->ss_elem_index[i] + e];
      int ielem_type = Elem_Type(exo, ielem);
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[i] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
      ss_node_list_len += num_nodes_on_side;
    }

    free(exo->ss_node_list[i]);
    exo->ss_node_list[i] = (int *)malloc(ss_node_list_len * sizeof(int));

    int offset = 0;
    for (int e = 0; e < exo->ss_num_sides[i]; e++) {
      int ielem = exo->ss_elem_list[exo->ss_elem_index[i] + e];
      int ielem_type = Elem_Type(exo, ielem);
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[i] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
      exo->ss_node_cnt_list[i][e] = num_nodes_on_side;
      for (int j = 0; j < num_nodes_on_side; j++) {
        int id = local_side_node_list[j];
        int iconnect = Proc_Connect_Ptr[ielem];
        int I = Proc_Elem_Connect[iconnect + id];
        exo->ss_node_list[i][j + offset] = I;
      }
      offset += num_nodes_on_side;
    }
    /*
     * Set up quick pointers for nodes on each given side of a sideset
     * that can be used later to find exactly where to go in the big
     * distribution factor list...
     */

    free(exo->ss_node_side_index[i]);
    exo->ss_node_side_index[i] = (int *)malloc((exo->ss_num_sides[i] + 1) * sizeof(int));

    exo->ss_node_side_index[i][0] = 0;

    for (int j = 0; j < exo->ss_num_sides[i]; j++) {
      exo->ss_node_side_index[i][j + 1] =
          (exo->ss_node_side_index[i][j] + exo->ss_node_cnt_list[i][j]);
    }
  }
  // char *format = "sgid%d.csv";
  // char fname[80];
  // snprintf(fname, 79, format, ProcID);
  // FILE *f = fopen(fname, "w");
  // fprintf(f,"x,y,z,gid,ss\n");
  ////for (int k = 0; k < exo->num_side_sets; k++) {
  ////  for (int i = 0; i < ss_edge_info[k].total_nodes; i++) {
  ////    int ln = in_list(ss_edge_info[k].global_node_ids[i], 0, exo->num_nodes,
  /// dpi->node_index_global); /  }
  ////}
  // for (int i = 0; i < exo->num_side_sets; i++) {
  //   int offset = 0;
  //   for (int e = 0; e < exo->ss_num_sides[i]; e++) {
  //     int ielem = exo->ss_elem_list[exo->ss_elem_index[i] + e];
  //     int ielem_type = Elem_Type(exo, ielem);
  //     int local_side_node_list[MAX_NODES_PER_SIDE];

  //    /* find SIDE info for primary side */
  //    int num_nodes_on_side = 0;

  //    int id_side = exo->ss_side_list[exo->ss_elem_index[i] + e];
  //    get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
  //    for (int j = 0; j < num_nodes_on_side; j++) {
  //      int ln = exo->ss_node_list[i][j + offset];
  //      fprintf(f, "%g,%g,%g,%d,%d\n",
  //          exo->x_coord[ln],
  //          exo->y_coord[ln],
  //          exo->z_coord[ln],
  //          dpi->node_index_global[ln],
  //          exo->ss_id[i]);
  //    }
  //    offset += num_nodes_on_side;
  //  }
  //}
  // fclose(f);
  // MPI_Finalize();
  // exit(0);

  free(requests);
  free(rotated_side_sets);
  for (int i = 0; i < num_rotated_side_sets; i++) {
    free(ss_edge_info[i].node_per_side);
    free(ss_edge_info[i].global_node_ids);
  }
  free(ss_edge_info);
  free(recv_nodes);
  free(send_nodes);

  for (int i = 0; i < dpi->num_neighbors * dpi->num_side_sets_global; i++) {
    if (recv_sides[i] > 0) {
      free(recv_nodes_per_side[i]);
    }
  }
  free(recv_nodes_per_side);
  free(recv_sides);
  free(send_sides);
  free(known_ss_elem_sides_count);
  for (int i = 0; i < exo->num_side_sets; i++) {
    free(known_ss_elem_sides_local[i]);
  }
  free(known_ss_elem_sides_local);
  for (int i = 0; i < dpi->num_neighbors * dpi->num_side_sets_global; i++) {
    free(global_nodes[i]);
  }
  free(global_nodes);
  free(offsets);
  free(global_to_local);
  free(ss_global_nodes);
  for (int i = 0; i < exo->num_side_sets; i++) {
    free(ss_elem_sides_local[i]);
  }
  free(ss_elem_sides_count_local);
  free(ss_elem_sides_local);

  return GOMA_SUCCESS;
}

goma_error setup_rotated_bc_nodes(
    Exo_DB *exo, Dpi *dpi, struct Boundary_Condition *bc_types, int num_bc, double *x) {

  if (num_bc == 0) {
    return GOMA_SUCCESS;
  }

  goma_rotation_node_s *rotations = NULL;
  bool *bc_is_rotated = NULL;
  goma_error error = setup_bc_is_rotated_list(bc_types, num_bc, &bc_is_rotated);
  GOMA_EH(error, "setup_bc_rotate_list");
  error = allocate_rotations(exo, &rotations);
  GOMA_EH(error, "allocate_rotations");

  bool *side_set_seen = malloc(sizeof(bool) * exo->num_side_sets);
  for (int i = 0; i < exo->num_side_sets; i++) {
    side_set_seen[i] = false;
  }

  struct node_normal {
    goma_normal **normals;
    int n_normals;
  };

  struct node_normal *node_normals = malloc(sizeof(struct node_normal) * exo->num_nodes);
  for (int i = 0; i < exo->num_nodes; i++) {
    node_normals[i].normals = NULL; // malloc(sizeof(goma_normal *) * GOMA_MAX_NORMALS_PER_NODE);
    node_normals[i].n_normals = 0;
    // for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
    //   node_normals[i].normals[j] = goma_normal_alloc(3);
    // }
  }

  int old_assemble_jacobian_setting = af->Assemble_Jacobian;
  af->Assemble_Jacobian = true;
  int err = 0;
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    if (bc_is_rotated[bc_index]) {
      int ss_index = bc_types[bc_index].Set_Index;
      if (ss_index == -1) { // ss isn't on this processor
        continue;
      }
      int vector_equation = vector_equation_from_equation(bc_types[bc_index].equation);

      for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
        int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

        err = load_elem_dofptr(ielem, exo, pg->matrices[pg->imtrx].x, pg->matrices[pg->imtrx].x_old,
                               pg->matrices[pg->imtrx].xdot, pg->matrices[pg->imtrx].xdot_old, 0);
        GOMA_EH(err, "load_elem_dofptr");
        err = bf_mp_init(pd);
        GOMA_EH(err, "bf_mp_init");

        int iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;
        int ielem_type = ei[pg->imtrx]->ielem_type;
        int num_local_nodes = ei[pg->imtrx]->num_local_nodes;
        int ielem_dim = ei[pg->imtrx]->ielem_dim;
        int local_side_node_list[MAX_NODES_PER_SIDE];

        /* find SIDE info for primary side */
        int num_nodes_on_side = 0;

        int id_side = exo->ss_side_list[exo->ss_elem_index[ss_index] + e];
        get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);

        /*
         * LOOP over NODES to which this condition applies
         */
        int num_ROT_nodes = 0;

        num_ROT_nodes = num_nodes_on_side;
        int *elem_node_id = &(local_side_node_list[0]);

        /* use nodal points only!! */
        for (int k = 0; k < num_ROT_nodes; k++) {
          int id = elem_node_id[k];
          int I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];
          if (!side_set_seen[ss_index]) {
            /* Find the local element node number for the current node */
            double xi[DIM];
            /* make sure we still need to calculate rotation vectors */
            find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);
            err = load_basis_functions(xi, bfd);
            GOMA_EH(err, "problem from load_basis_functions");
            err = beer_belly();
            GOMA_EH(err, "beer_belly");
            err = load_fv();
            GOMA_EH(err, "load_fv");
            err = load_bf_grad();
            GOMA_EH(err, "load_bf_grad");
            err = load_bf_mesh_derivs();
            GOMA_EH(err, "load_bf_mesh_derivs");

            /* put NORMAL vector into array */
            /* calculate the determinant of the surface jacobian  and the normal to
             * the surface all at one time */
            surface_determinant_and_normal(ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes,
                                           ielem_dim - 1, id_side, num_nodes_on_side,
                                           local_side_node_list);

            int n_index = node_normals[I].n_normals;
            node_normals[I].n_normals++;
            if (n_index == 0) {
              node_normals[I].normals = malloc(sizeof(goma_normal *) * GOMA_MAX_NORMALS_PER_NODE);
              for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
                node_normals[I].normals[j] = goma_normal_alloc(3);
              }
            }
            if (node_normals[I].n_normals > GOMA_MAX_NORMALS_PER_NODE) {
              GOMA_EH(GOMA_ERROR, "GOMA_MAX_NORMALS_PER_NODE too small, currently %d",
                      GOMA_MAX_NORMALS_PER_NODE);
            }
            gds_vector_set(node_normals[I].normals[n_index]->normal, 0, fv->snormal[0]);
            gds_vector_set(node_normals[I].normals[n_index]->normal, 1, fv->snormal[1]);
            gds_vector_set(node_normals[I].normals[n_index]->normal, 2, fv->snormal[2]);
            for (int i = 0; i < 3; i++) {
              for (int j = 0; j < MDE; j++) {
                gds_vector_set(node_normals[I].normals[n_index]->d_normal_dx[i][j], 0,
                               fv->dsnormal_dx[0][i][j]);
                gds_vector_set(node_normals[I].normals[n_index]->d_normal_dx[i][j], 1,
                               fv->dsnormal_dx[1][i][j]);
                gds_vector_set(node_normals[I].normals[n_index]->d_normal_dx[i][j], 2,
                               fv->dsnormal_dx[2][i][j]);
              }
            }
          }
          rotations[I].is_rotated = true;
          rotations[I].eqn_is_rotated[vector_equation] = true;
        }
      }

      side_set_seen[ss_index] = true;
    }
  }
  af->Assemble_Jacobian = old_assemble_jacobian_setting;

  free(side_set_seen);

  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotations[i].is_rotated) {
      goma_error err =
          goma_best_coordinate_system_3D(node_normals[i].normals, node_normals[i].n_normals,
                                         rotations[i].rotated_coord, &rotations[i].type);
      GOMA_EH(err, "find best coordinate error for node %d, %g %g %g", i, exo->x_coord[i],
              exo->y_coord[i], exo->z_coord[i]);
    }
  }

  // char fname[80] = "normals.csv";
  ////multiname(fname, ProcID, Num_Proc);
  // sprintf(fname, "normals-%d-%d.csv", ProcID, Num_Proc);

  // FILE *file = fopen(fname, "w");
  // fprintf(file, "x,y,z,nx,ny,nz,normal_id,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,ProcID,type\n");
  // for (int i = 0; i < exo->num_nodes; i++) {
  //   for (int j = 0; j < node_normals[i].n_normals; j++) {
  //   fprintf(file, "%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d\n",
  //       exo->x_coord[i],
  //       exo->y_coord[i],
  //       exo->z_coord[i],
  //       node_normals[i].normals[j]->normal->data[0],
  //       node_normals[i].normals[j]->normal->data[1],
  //       node_normals[i].normals[j]->normal->data[2],
  //       i,
  //       rotations[i].rotated_coord[0]->normal->data[0],
  //       rotations[i].rotated_coord[0]->normal->data[1],
  //       rotations[i].rotated_coord[0]->normal->data[2],
  //       rotations[i].rotated_coord[1]->normal->data[0],
  //       rotations[i].rotated_coord[1]->normal->data[1],
  //       rotations[i].rotated_coord[1]->normal->data[2],
  //       rotations[i].rotated_coord[2]->normal->data[0],
  //       rotations[i].rotated_coord[2]->normal->data[1],
  //       rotations[i].rotated_coord[2]->normal->data[2],
  //       ProcID, rotations[i].type);
  //   }
  // }
  // fprintf(file, "\n");
  // fclose(file);
  // MPI_Barrier(MPI_COMM_WORLD);
  // exit(0);

  goma_automatic_rotations.automatic_rotations = true;
  goma_automatic_rotations.rotation_nodes = rotations;
  for (int i = 0; i < exo->num_nodes; i++) {
    if (node_normals[i].normals != NULL) {
      for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
        goma_normal_free(node_normals[i].normals[j]);
      }
      free(node_normals[i].normals);
    }
  }
  free(node_normals);
  free(bc_is_rotated);
  return GOMA_SUCCESS;
}

int vector_equation_from_equation(int equation) {
  switch (equation) {
  case R_MESH1:
  case R_MESH2:
  case R_MESH3:
  case R_MESH_NORMAL:
  case R_MESH_TANG1:
  case R_MESH_TANG2:
    return VECT_EQ_MESH;
  case R_MOMENTUM1:
  case R_MOMENTUM2:
  case R_MOMENTUM3:
  case R_MOM_NORMAL:
  case R_MOM_TANG1:
  case R_MOM_TANG2:
    return VECT_EQ_MOM;
  default:
    return -1;
  }
}

int offset_from_rotated_equation(int eqn) {
  switch (eqn) {
  case R_MOM_NORMAL:
    return 0;
  case R_MOM_TANG1:
    return 1;
  case R_MOM_TANG2:
    return 2;
  case R_MESH_NORMAL:
    return 0;
  case R_MESH_TANG1:
    return 1;
  case R_MESH_TANG2:
    return 2;
  default:
    return -1;
  }
}

int first_equation_from_vector_equation(int eqn) {
  switch (eqn) {
  case R_MESH1:
  case R_MESH2:
  case R_MESH3:
    return R_MESH1;
  case R_MOMENTUM1:
  case R_MOMENTUM2:
  case R_MOMENTUM3:
    return R_MOMENTUM1;
  default:
    return -1;
  }
}
