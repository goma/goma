#include "bc/rotate_coordinates.h"

#include <assert.h>
#include <math.h>
#include <rf_bc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#define GOMA_MAX_NORMALS_PER_NODE 14
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

goma_error setup_bc_is_rotated_list(struct Boundary_Condition *bc_types,
                                    int num_bc,
                                    bool **bc_is_rotated_list) {
  assert(num_bc > 0);
  *bc_is_rotated_list = calloc((size_t)num_bc, sizeof(bool));
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    (*bc_is_rotated_list)[bc_index] = false;
    if (!strcmp(bc_types[bc_index].Set_Type, "SS")) {
      bool bc_is_rotated;
      check_if_equation_is_rotation(BC_Desc[bc_types[bc_index].BC_Desc_index].equation,
                                    &bc_is_rotated);
      if (bc_is_rotated) {
        (*bc_is_rotated_list)[bc_index] = true;
      }
    }
  }
  return GOMA_SUCCESS;
}

goma_error allocate_rotations(Exo_DB *exo, goma_rotation_node_s **rotations) {
  if (goma_automatic_rotations.rotation_nodes != NULL) {
    for (int i = 0; i < exo->num_nodes; i++) {
      for (int j = 0; j < DIM; j++) {
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].rotated_coord[j]);
      }
    }
    free(goma_automatic_rotations.rotation_nodes);
  }
  *rotations = calloc((size_t)exo->num_nodes, sizeof(goma_rotation_node_s));
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < DIM; j++) {
      (*rotations)[i].rotated_coord[j] = gds_vector_alloc(3);
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

goma_error
setup_rotated_bc_nodes(Exo_DB *exo, struct Boundary_Condition *bc_types, int num_bc, double *x) {

  if (num_bc == 0) {
    return GOMA_SUCCESS;
  }

  goma_rotation_node_s *rotations;
  bool *bc_is_rotated;
  goma_error error;
  error = setup_bc_is_rotated_list(bc_types, num_bc, &bc_is_rotated);
  EH(error, "setup_bc_rotate_list");
  error = allocate_rotations(exo, &rotations);
  EH(error, "allocate_rotations");

  bool *side_set_seen = malloc(sizeof(bool) * exo->num_side_sets);
  for (int i = 0; i < exo->num_side_sets; i++) {
    side_set_seen[i] = false;
  }

  struct node_normal {
    gds_vector **normals;
    int n_normals;
  };

  struct node_normal *node_normals = malloc(sizeof(struct node_normal) * exo->num_nodes);
  for (int i = 0; i < exo->num_nodes; i++) {
    node_normals[i].normals = malloc(sizeof(gds_vector *)*GOMA_MAX_NORMALS_PER_NODE);
    node_normals[i].n_normals = 0;
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      node_normals[i].normals[j] = gds_vector_alloc(3);
    }
  }


  int err;
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    if (bc_is_rotated[bc_index]) {
      int ss_index = bc_types[bc_index].Set_Index;
      if (ss_index == -1) { // ss isn't on this processor
        continue;
      }
      // only operate on a side set once
      if (side_set_seen[ss_index]) {
        continue;
      } else {
        side_set_seen[ss_index] = true;
      }
      for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
        int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

        err = load_elem_dofptr(ielem, exo, pg->matrices[pg->imtrx].x, pg->matrices[pg->imtrx].x_old,
                               pg->matrices[pg->imtrx].xdot, pg->matrices[pg->imtrx].xdot_old,
                               0);
        err = bf_mp_init(pd);

        int iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;
        int ielem_type = ei[pg->imtrx]->ielem_type;
        int num_local_nodes = ei[pg->imtrx]->num_local_nodes;
        int ielem_dim = ei[pg->imtrx]->ielem_dim;
        int local_side_node_list[MAX_NODES_PER_SIDE];

        /* find SIDE info for primary side */
        int num_nodes_on_side;

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
          /* Find the local element node number for the current node */
          int id = elem_node_id[k];
          int I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + id];
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
            EH(GOMA_ERROR, "GOMA_MAX_NORMALS_PER_NODE too small, currently %d", GOMA_MAX_NORMALS_PER_NODE);
          }
          gds_vector_set(node_normals[I].normals[n_index], 0, fv->snormal[0]);
          gds_vector_set(node_normals[I].normals[n_index], 1, fv->snormal[1]);
          gds_vector_set(node_normals[I].normals[n_index], 2, fv->snormal[2]);
          rotations[I].is_rotated = true;
        }
      }
    }
  }

  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotations[i].is_rotated) {
      goma_best_coordinate_system_3D(node_normals[i].normals, node_normals[i].n_normals, rotations[i].rotated_coord);
    }
  }

  goma_automatic_rotations.automatic_rotations = true;
  goma_automatic_rotations.rotation_nodes = rotations;
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      gds_vector_free(node_normals[i].normals[j]);
    }
    free(node_normals[i].normals);
  }
  free(node_normals);
  free(bc_is_rotated);
  return GOMA_SUCCESS;
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
