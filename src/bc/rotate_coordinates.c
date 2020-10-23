#include "bc/rotate_coordinates.h"

#include <assert.h>
#include <math.h>
#include <rf_bc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "bc/rotate_util.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "gds/gds_vector.h"
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
#include "std.h"
#include "stdbool.h"

#ifndef GOMA_MAX_NORMALS_PER_NODE
#define GOMA_MAX_NORMALS_PER_NODE 20
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

static int
int_compare(const void *left, const void *right)
{
  int left_val = *((int *) left);
  int right_val = *((int *) right);
  if (left_val > right_val)
    return 1;
  else if (left_val < right_val)
    return -1;
  else
    return 0;
}

static int
int_pair_compare_first(const void *left, const void *right)
{
  int_pair * left_val = ((int_pair *) left);
  int_pair * right_val = ((int_pair *) right);
  if (left_val->first > right_val->first) {
    return 1;
  } else if (left_val->first < right_val->first) {
    return -1;
  } else {
    return 0;
  }
}


static int
int_pair_compare(const void *left, const void *right)
{
  int_pair * left_val = ((int_pair *) left);
  int_pair * right_val = ((int_pair *) right);
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

goma_error exchange_neighbor_ss_edges(
    Exo_DB *exo, Dpi *dpi, struct Boundary_Condition *bc_types, int num_bc, bool *bc_is_rotated, int_pair ***ss_elem_sides, int **ss_elem_sides_count) {
  int *rotated_side_sets = calloc(exo->num_side_sets, sizeof(int));
  int num_rotated_side_sets = 0;
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    if (bc_is_rotated[bc_index]) {
      int ss_index = bc_types[bc_index].Set_Index;
      if (ss_index == -1) { // ss isn't on this processor
        continue;
      }
      if ((in_list(ss_index, 0, num_rotated_side_sets, rotated_side_sets) == -1)) {
        rotated_side_sets[num_rotated_side_sets++] = ss_index;
      }
    }
  }

  ss_edge_share *ss_edge_info = calloc(dpi->num_side_sets_global, sizeof(ss_edge_share));
  for (int i = 0; i < num_rotated_side_sets; i++) {
    int ss_index = rotated_side_sets[i];
    ss_edge_info[i].ss = exo->ss_id[ss_index];
    ss_edge_info[i].num_sides = exo->ss_num_sides[ss_index];
    ss_edge_info[i].node_per_side = calloc(ss_edge_info[i].num_sides, sizeof(int));
    int total_nodes = 0;
    for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
      int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];
      load_ei(ielem, exo, 0, pg->imtrx);
      int ielem_type = ei[pg->imtrx]->ielem_type;
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[ss_index] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
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

      load_ei(ielem, exo, 0, pg->imtrx);

      int ielem_type = ei[pg->imtrx]->ielem_type;
      int local_side_node_list[MAX_NODES_PER_SIDE];

      /* find SIDE info for primary side */
      int num_nodes_on_side = 0;

      int id_side = exo->ss_side_list[exo->ss_elem_index[ss_index] + e];
      get_side_info(ielem_type, id_side, &num_nodes_on_side, local_side_node_list);
      int *elem_node_id = &(local_side_node_list[0]);

      /* use nodal points only!! */
      for (int k = 0; k < num_nodes_on_side; k++) {
        int id = elem_node_id[k];
        int I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];
        int global_node = dpi->node_index_global[I];
        ss_edge_info[i].global_node_ids[offset++] = global_node;
      }
    }
    assert(offset == ss_edge_info[i].total_nodes);
  }


  MPI_Request *requests = malloc(sizeof(MPI_Request) * 2 * dpi->num_neighbors);

  int *recv_nodes = calloc(dpi->num_neighbors * dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_neighbors; i++) {
    MPI_Irecv(&recv_nodes[i*dpi->num_side_sets_global], dpi->num_side_sets_global,
              MPI_INT, dpi->neighbor[i], 102+dpi->neighbor[i], MPI_COMM_WORLD, &(requests[i]));
  }

  int *send_nodes = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int j = 0; j < dpi->num_side_sets_global; j++) {
    send_nodes[j] = ss_edge_info[j].total_nodes;
  }
  for (int i = 0; i < dpi->num_neighbors; i++) {
    MPI_Isend(send_nodes, dpi->num_side_sets_global, MPI_INT, dpi->neighbor[i], 102+ProcID,
              MPI_COMM_WORLD, &(requests[dpi->num_neighbors + i]));
  }
  MPI_Waitall(2*dpi->num_neighbors, requests, MPI_STATUSES_IGNORE);

  int **global_nodes = calloc(dpi->num_side_sets_global, sizeof(int *));
  int *ss_global_nodes = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    int all_nodes = 0;
    for (int j = 0; j < dpi->num_neighbors; j++) {
      all_nodes += recv_nodes[j * dpi->num_side_sets_global + i];
    }
    global_nodes[i] = calloc(all_nodes, sizeof(int));
    ss_global_nodes[i] = all_nodes;
  }

  int *offsets = calloc(dpi->num_side_sets_global, sizeof(int));
  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      int count = recv_nodes[i * dpi->num_side_sets_global + j];
      printf("SS %d, %d recv %d from %d\n", j, ProcID, count, dpi->neighbor[i]);
      MPI_Irecv(&global_nodes[j][offsets[j]], count, MPI_INT, dpi->neighbor[i],
                103 + j + dpi->neighbor[i], MPI_COMM_WORLD, &(requests[i]));
      offsets[j] += count;
    }
  }

  for (int i = 0; i < dpi->num_neighbors; i++) {
    for (int j = 0; j < dpi->num_side_sets_global; j++) {
      printf("SS %d, %d send %d to %d\n", j, ProcID, ss_edge_info[j].total_nodes, dpi->neighbor[i]);
      MPI_Isend(ss_edge_info[j].global_node_ids, ss_edge_info[j].total_nodes, MPI_INT,
                dpi->neighbor[i], 103 + j + ProcID, MPI_COMM_WORLD,
                &(requests[dpi->num_neighbors + i]));
    }
  }

  MPI_Waitall(2*dpi->num_neighbors, requests, MPI_STATUSES_IGNORE);

  int_pair *global_to_local = calloc(exo->num_nodes, sizeof(int_pair));
  for (int i = 0; i < exo->num_nodes; i++) {
    global_to_local[i].first = dpi->node_index_global[i];
    global_to_local[i].second = i;
  }
  qsort(global_to_local, exo->num_nodes, sizeof(int_pair), int_pair_compare);

  // convert all to local
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    for (int j = 0; j < ss_global_nodes[i]; j++) {
      int_pair global_node = {global_nodes[i][j], -1};
      int_pair* match = (int_pair *) bsearch(&global_node, global_to_local, exo->num_nodes, sizeof(int_pair), int_pair_compare_first);
      if (match != NULL) {
        global_nodes[i][j] = match->second;
      } else {
        global_nodes[i][j] = -1;
      }
    }
    for (int j = 0; j < ss_edge_info[i].total_nodes; j++) {
      int_pair global_node = {ss_edge_info[i].global_node_ids[j], -1};
      int_pair* match = (int_pair *) bsearch(&global_node, global_to_local, exo->num_nodes, sizeof(int_pair), int_pair_compare_first);
      if (match != NULL) {
        ss_edge_info[i].global_node_ids[j] = match->second;
      } else {
        EH(GOMA_ERROR, "no mapping to local node for local ss node");
      }
    }
    qsort(global_nodes[i], ss_global_nodes[i], sizeof(int), int_compare);
    qsort(ss_edge_info[i].global_node_ids, ss_edge_info[i].total_nodes, sizeof(int), int_compare);
  }

  int_pair **ss_elem_sides_local = calloc(sizeof(int_pair *), dpi->num_side_sets_global);
  int *ss_elem_sides_count_local = calloc(sizeof(int), dpi->num_side_sets_global);
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    if (ss_global_nodes[i] > 0) {
      ss_elem_sides_local[i] = calloc(sizeof(int_pair), ss_global_nodes[i]);
    } else {
      ss_elem_sides_local[i] = NULL;
    }
  }

  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    for (int j = 0; j < ss_global_nodes[i]; j++) {
      int node = global_nodes[i][j];
      if (node != -1) {
        for (int idx = exo->node_elem_pntr[node]; idx < exo->node_elem_pntr[node+1]; idx++) {
          int elem = exo->node_elem_list[idx];
          int ielem_type = exo->eb_elem_itype[ exo->elem_eb[elem] ];
          int shape = type2shape(ielem_type);
          int n_sides = shape2sides(shape);
          for (int side = 0; side < n_sides; side++) {
            int local_side_node_list[MAX_NODES_PER_SIDE];
            int side_nodes[MAX_NODES_PER_SIDE];

            /* find SIDE info for primary side */
            int num_nodes_on_side = 0;
            get_side_info(ielem_type, side+1, &num_nodes_on_side, local_side_node_list);

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
              bool all_known = true;
              for (int k = 0; k < num_nodes_on_side; k++) {
                int sn = side_nodes[k];
                int *p = bsearch(&sn, global_nodes[i], ss_global_nodes[i], sizeof(int), int_compare);
                int *q = bsearch(&sn, ss_edge_info[i].global_node_ids, ss_edge_info[i].total_nodes, sizeof(int), int_compare);
                if (q == NULL) {
                  all_known = false;
                }
                if (p == NULL && q == NULL) {
                  found_all = false;
                }
              }

              if (found_all && !all_known) {
                if (ss_elem_sides_count_local[i] > 1) {
                  qsort(ss_elem_sides_local[i], ss_elem_sides_count_local[i], sizeof(int_pair), int_pair_compare);
                }
                int_pair elem_side = {elem, side};
                int_pair *exists = bsearch(&elem_side, ss_elem_sides_local[i], ss_elem_sides_count_local[i], sizeof(int_pair), int_pair_compare);
                if (exists == NULL) {
                  ss_elem_sides_local[i][ss_elem_sides_count_local[i]] = elem_side;
                  ss_elem_sides_count_local[i] += 1;
                }
              }
            }
          }

        }
      }
    }
  }

  *ss_elem_sides = ss_elem_sides_local;
  *ss_elem_sides_count = ss_elem_sides_count_local;
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
  EH(error, "setup_bc_rotate_list");
  error = allocate_rotations(exo, &rotations);
  EH(error, "allocate_rotations");

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
    node_normals[i].normals = malloc(sizeof(goma_normal *) * GOMA_MAX_NORMALS_PER_NODE);
    node_normals[i].n_normals = 0;
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      node_normals[i].normals[j] = goma_normal_alloc(3);
    }
  }

  int_pair **ss_elem_sides;
  int *ss_elem_sides_count;
  exchange_neighbor_ss_edges(exo, dpi, bc_types, num_bc, bc_is_rotated, &ss_elem_sides, &ss_elem_sides_count);

  int old_assemble_jacobian_setting = af->Assemble_Jacobian;
  af->Assemble_Jacobian = true;
  int err = 0;
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    if (bc_is_rotated[bc_index]) {
      int ss_index = bc_types[bc_index].Set_Index;
      int global_ss_index = -1;
      for (int i = 0; i < dpi->num_side_sets_global; i++) {
        if (dpi->ss_id_global[i] == exo->ss_id[ss_index]) {
          global_ss_index = i;
          break;
        }
      }
      if (ss_index == -1) { // ss isn't on this processor
        continue;
      }
      int vector_equation = vector_equation_from_equation(bc_types[bc_index].equation);

      for (int e = 0; e < ss_elem_sides_count[global_ss_index]; e++) {
        int ielem = ss_elem_sides[global_ss_index][e].first;

        err = load_elem_dofptr(ielem, exo, pg->matrices[pg->imtrx].x, pg->matrices[pg->imtrx].x_old,
                               pg->matrices[pg->imtrx].xdot, pg->matrices[pg->imtrx].xdot_old, 0);
        EH(err, "load_elem_dofptr");
        err = bf_mp_init(pd);
        EH(err, "bf_mp_init");

        int iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;
        int ielem_type = ei[pg->imtrx]->ielem_type;
        int num_local_nodes = ei[pg->imtrx]->num_local_nodes;
        int ielem_dim = ei[pg->imtrx]->ielem_dim;
        int local_side_node_list[MAX_NODES_PER_SIDE];

        /* find SIDE info for primary side */
        int num_nodes_on_side = 0;

        int id_side = ss_elem_sides[global_ss_index][e].second + 1;
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
            EH(err, "problem from load_basis_functions");
            err = beer_belly();
            EH(err, "beer_belly");
            err = load_fv();
            EH(err, "load_fv");
            err = load_bf_grad();
            EH(err, "load_bf_grad");
            err = load_bf_mesh_derivs();
            EH(err, "load_bf_mesh_derivs");

            /* put NORMAL vector into array */
            /* calculate the determinant of the surface jacobian  and the normal to
             * the surface all at one time */
            surface_determinant_and_normal(ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes,
                                           ielem_dim - 1, id_side, num_nodes_on_side,
                                           local_side_node_list);

            int n_index = node_normals[I].n_normals;
            node_normals[I].n_normals++;
            if (node_normals[I].n_normals > GOMA_MAX_NORMALS_PER_NODE) {
              EH(GOMA_ERROR, "GOMA_MAX_NORMALS_PER_NODE too small, currently %d",
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
      for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
        int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

        err = load_elem_dofptr(ielem, exo, pg->matrices[pg->imtrx].x, pg->matrices[pg->imtrx].x_old,
                               pg->matrices[pg->imtrx].xdot, pg->matrices[pg->imtrx].xdot_old, 0);
        EH(err, "load_elem_dofptr");
        err = bf_mp_init(pd);
        EH(err, "bf_mp_init");

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
            EH(err, "problem from load_basis_functions");
            err = beer_belly();
            EH(err, "beer_belly");
            err = load_fv();
            EH(err, "load_fv");
            err = load_bf_grad();
            EH(err, "load_bf_grad");
            err = load_bf_mesh_derivs();
            EH(err, "load_bf_mesh_derivs");

            /* put NORMAL vector into array */
            /* calculate the determinant of the surface jacobian  and the normal to
             * the surface all at one time */
            surface_determinant_and_normal(ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes,
                                           ielem_dim - 1, id_side, num_nodes_on_side,
                                           local_side_node_list);

            int n_index = node_normals[I].n_normals;
            node_normals[I].n_normals++;
            if (node_normals[I].n_normals > GOMA_MAX_NORMALS_PER_NODE) {
              EH(GOMA_ERROR, "GOMA_MAX_NORMALS_PER_NODE too small, currently %d",
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

  int num_rotated_nodes = 0;
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotations[i].is_rotated) {
      num_rotated_nodes++;
    }
  }

  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotations[i].is_rotated) {
      goma_best_coordinate_system_3D(node_normals[i].normals, node_normals[i].n_normals,
                                     rotations[i].rotated_coord);
    }
  }

  goma_automatic_rotations.automatic_rotations = true;
  goma_automatic_rotations.rotation_nodes = rotations;
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      goma_normal_free(node_normals[i].normals[j]);
    }
    free(node_normals[i].normals);
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
