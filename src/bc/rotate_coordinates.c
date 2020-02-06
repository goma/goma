#include "bc/rotate_coordinates.h"

#include <assert.h>
#include <math.h>
#include <rf_bc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static const double critical_angle_radians = GOMA_ROTATION_CRITICAL_ANGLE * M_PI / 180;

#ifdef DEBUG_AUTO_ROTATE
static void
write_rotations_to_file(const char *filename, Exo_DB *exo, goma_rotation_node_s *rotations);
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
      for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].normals[j]);
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].average_normals[j]);
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].tangent1s[j]);
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].tangent2s[j]);
      }
      for (int j = 0; j < DIM; j++) {
        gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].rotated_coord[j]);
      }
    }
    free(goma_automatic_rotations.rotation_nodes);
  }
  *rotations = calloc((size_t)exo->num_nodes, sizeof(goma_rotation_node_s));
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      (*rotations)[i].normals[j] = gds_vector_alloc(3);
      (*rotations)[i].average_normals[j] = gds_vector_alloc(3);
      (*rotations)[i].tangent1s[j] = gds_vector_alloc(3);
      (*rotations)[i].tangent2s[j] = gds_vector_alloc(3);
    }
    for (int j = 0; j < DIM; j++) {
      (*rotations)[i].rotated_coord[j] = gds_vector_alloc(3);
    }
  }

  return GOMA_SUCCESS;
}
goma_error free_rotations(Exo_DB *exo, goma_rotation_node_s **rotations) {
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      free((*rotations)[i].normals[j]);
      free((*rotations)[i].average_normals[j]);
      free((*rotations)[i].tangent1s[j]);
      free((*rotations)[i].tangent2s[j]);
    }
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

        if (x == x_static) /* be the least disruptive possible */
        {
          err =
              load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static, 0);
        } else {
          err = load_elem_dofptr(ielem, exo, x, x, x, x, 0);
        }
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

          int n_index = rotations[I].n_normals;
          rotations[I].n_normals++;
          assert(rotations[I].n_normals < GOMA_MAX_NORMALS_PER_NODE);
          gds_vector_set(rotations[I].normals[n_index], 0, fv->snormal[0]);
          gds_vector_set(rotations[I].normals[n_index], 1, fv->snormal[1]);
          gds_vector_set(rotations[I].normals[n_index], 2, fv->snormal[2]);
          rotations[I].element[n_index] = ielem;
          rotations[I].face[n_index] = id_side;
        }
      }
    }
  }

  error = set_rotation_types(exo, rotations);
  EH(error, "set_rotation_types");
  error = set_average_normals_and_tangents(exo, rotations);
  EH(error, "set_average_normals_and_tangents");
  error = associate_directions(exo, rotations);
  EH(error, "associate_directions");
  error = set_rotated_coordinate_system(exo, rotations);
  EH(error, "set_rotated_coordinate_system");
  error = set_face_normal_association(exo, rotations);
  EH(error, "set_face_normal_association");

#ifdef DEBUG_AUTO_ROTATE
  write_rotations_to_file("normals.csv", exo, rotations);
#endif

  goma_automatic_rotations.automatic_rotations = true;
  goma_automatic_rotations.rotation_nodes = rotations;
  free(bc_is_rotated);
  return GOMA_SUCCESS;
}

#ifdef DEBUG_AUTO_ROTATE
static void
write_rotations_to_file(const char *filename, Exo_DB *exo, goma_rotation_node_s *rotations) {
  FILE *normalscsv;
  normalscsv = fopen(filename, "w");

  fprintf(normalscsv, "n_normals,type,node,x,y,z,nx,ny,nz,dir,ax,ay,az,t1x,t1y,t1z,t2x,t2y,t2z,c1x,"
                      "c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,element,face,face_association\n");
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < rotations[i].n_normals; j++) {
      fprintf(normalscsv,
              "%d,%d,%d,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,"
              "%d,%d,%d\n",
              rotations[i].n_normals, rotations[i].type, i, exo->x_coord[i], exo->y_coord[i],
              exo->z_coord[i], gds_vector_get(rotations[i].normals[j], 0),
              gds_vector_get(rotations[i].normals[j], 1),
              gds_vector_get(rotations[i].normals[j], 2), rotations[i].associate_direction[j],
              gds_vector_get(rotations[i].average_normals[j], 0),
              gds_vector_get(rotations[i].average_normals[j], 1),
              gds_vector_get(rotations[i].average_normals[j], 2),
              gds_vector_get(rotations[i].tangent1s[j], 0),
              gds_vector_get(rotations[i].tangent1s[j], 1),
              gds_vector_get(rotations[i].tangent1s[j], 2),
              gds_vector_get(rotations[i].tangent2s[j], 0),
              gds_vector_get(rotations[i].tangent2s[j], 1),
              gds_vector_get(rotations[i].tangent2s[j], 2),
              gds_vector_get(rotations[i].rotated_coord[0], 0),
              gds_vector_get(rotations[i].rotated_coord[0], 1),
              gds_vector_get(rotations[i].rotated_coord[0], 2),
              gds_vector_get(rotations[i].rotated_coord[1], 0),
              gds_vector_get(rotations[i].rotated_coord[1], 1),
              gds_vector_get(rotations[i].rotated_coord[1], 2),
              gds_vector_get(rotations[i].rotated_coord[2], 0),
              gds_vector_get(rotations[i].rotated_coord[2], 1),
              gds_vector_get(rotations[i].rotated_coord[2], 2), rotations[i].element[j],
              rotations[i].face[j], rotations[i].face_coordinate_association[j]);
    }
  }
  fclose(normalscsv);
}
#endif

/**
 * @brief goma_rotation_node_type_surface
 * @param node the node we are looking at rotating
 * @return rotation at this node is type surface
 */
bool goma_rotation_node_type_surface(const goma_rotation_node_s *node) {
  double min_dot = 1e300;
  for (int u_index = 0; u_index < node->n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < node->n_normals; v_index++) {
      double dot = gds_vector_dot(node->normals[u_index], node->normals[v_index]);
      if (dot < min_dot) {
        min_dot = dot;
      }
    }
  }

  if (min_dot > cos(critical_angle_radians)) {
    return true;
  }
  return false;
}

bool goma_rotation_node_type_corner(const goma_rotation_node_s *node) {
  if (node->n_normals < 3) {
    return false;
  }

  // If two normals are larger than critical angle then there should be a third normal larger
  // than a critical angle to both those to be a corner type
  for (int u_index = 0; u_index < node->n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < node->n_normals; v_index++) {
      double uv_dot = gds_vector_dot(node->normals[u_index], node->normals[v_index]);
      if (uv_dot < cos(critical_angle_radians)) {
        for (int w_index = 0; w_index < node->n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index)
            continue;
          double uw_dot = gds_vector_dot(node->normals[u_index], node->normals[w_index]);
          double vw_dot = gds_vector_dot(node->normals[v_index], node->normals[w_index]);
          if (uw_dot < cos(critical_angle_radians) && vw_dot < cos(critical_angle_radians)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool goma_rotation_node_type_edge(goma_rotation_node_s *node) {
  bool is_edge_type = false;
  if (node->n_normals < 2) {
    return false;
  }

  gds_vector *critical_normal[2];
  critical_normal[0] = gds_vector_alloc(3);
  critical_normal[1] = gds_vector_alloc(3);

  gds_vector_copy(critical_normal[0], node->normals[0]);
  // find another critical normal
  bool found = false;
  for (int u_index = 0; u_index < node->n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < node->n_normals; v_index++) {
      double dot = gds_vector_dot(node->normals[u_index], node->normals[v_index]);
      if (dot < cos(critical_angle_radians)) {
        gds_vector_copy(critical_normal[0], node->normals[u_index]);
        node->critical_normal_index[0] = u_index;
        gds_vector_copy(critical_normal[1], node->normals[v_index]);
        node->critical_normal_index[1] = v_index;
        found = true;
      }
    }
  }

  if (found) {
    // make sure we only have two critical normals
    is_edge_type = true;
    for (int u_index = 0; u_index < node->n_normals; u_index++) {
      if ((gds_vector_dot(critical_normal[0], node->normals[u_index]) <
           cos(critical_angle_radians)) &&
          (gds_vector_dot(critical_normal[1], node->normals[u_index]) <
           cos(critical_angle_radians))) {
        is_edge_type = false;
        break;
      }
    }

    if (is_edge_type) {
      // associate to a critical normal
      for (int u_index = 0; u_index < node->n_normals; u_index++) {
        if (gds_vector_dot(critical_normal[0], node->normals[u_index]) >=
            gds_vector_dot(critical_normal[1], node->normals[u_index])) {
          node->associate_critical_normal[u_index] = 0;
        } else {
          node->associate_critical_normal[u_index] = 1;
        }
      }
    }
  }

  gds_vector_free(critical_normal[0]);
  gds_vector_free(critical_normal[1]);
  return is_edge_type;
}

goma_error set_rotation_types(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotation[i].n_normals > 0) {
      if (rotation[i].n_normals == 1) {
        rotation[i].type = GOMA_ROTATION_SINGLE;
      } else if (goma_rotation_node_type_surface(&(rotation[i]))) {
        rotation[i].type = GOMA_ROTATION_SURFACE;
      } else if (goma_rotation_node_type_corner(&(rotation[i]))) {
        rotation[i].type = GOMA_ROTATION_CORNER;
      } else if (goma_rotation_node_type_edge(&(rotation[i]))) {
        rotation[i].type = GOMA_ROTATION_EDGE;
      } else {
        return GOMA_ERROR;
      }
    }
  }

  return GOMA_SUCCESS;
}

goma_error associate_directions(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotation[i].n_normals > 0) {
      bool dir_set[3] = {false, false, false};
      int associated_directions = 0;

      switch (rotation[i].type) {
      case GOMA_ROTATION_SINGLE:
      case GOMA_ROTATION_SURFACE:
      case GOMA_ROTATION_CORNER:
        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          unsigned int best_dir = 0;
          double max_dot = 0;
          for (unsigned int card = 0; card < DIM; card++) {
            if (dir_set[card]) {
              // check if we are within critical angle to already set
              // if not pick next best
              bool skip = false;
              for (int v_index = 0; v_index < u_index; v_index++) {
                if (rotation[i].associate_direction[v_index] == card) {
                  if (gds_vector_dot(rotation[i].normals[v_index], rotation[i].normals[u_index]) <
                      cos(critical_angle_radians)) {
                    skip = true;
                    break;
                  }
                }
              }
              if (skip)
                continue;
            }
            gds_vector *cvect = gds_vector_alloc(3);
            gds_vector_zero(cvect);
            gds_vector_set(cvect, card, 1.0);
            double dot = fabs(gds_vector_dot(cvect, rotation[i].average_normals[u_index]));
            if (dot > max_dot && fabs(max_dot - dot) > 1e-6) { // add tolerance
              best_dir = card;
              max_dot = dot;
            }
            gds_vector_free(cvect);
          }
          rotation[i].associate_direction[u_index] = best_dir;
          dir_set[best_dir] = true;
        }
        break;
      case GOMA_ROTATION_EDGE: {
        double ca_max[2] = {0.0, 0.0};
        unsigned int ca_coord[2] = {0, 0};
        gds_vector *ca1 = rotation[i].average_normals[rotation[i].critical_normal_index[0]];
        gds_vector *ca2 = rotation[i].average_normals[rotation[i].critical_normal_index[1]];
        for (unsigned int cord = 0; cord < DIM; cord++) {
          double dot1 = fabs(gds_vector_get(ca1, cord));
          double dot2 = fabs(gds_vector_get(ca2, cord));
          if (dot1 > ca_max[0] && dot1 >= dot2) {
            ca_max[0] = dot1;
            ca_coord[0] = cord;
          } else if (dot2 > ca_max[1]) {
            ca_max[1] = dot2;
            ca_coord[1] = cord;
          }
        }

        if (ca_max[0] < 1e-12 || ca_max[1] < 1e-12 || ca_coord[0] == ca_coord[1]) {
          return GOMA_ERROR;
        }

        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          if (rotation[i].associate_critical_normal[u_index] == 0) {
            rotation[i].associate_direction[u_index] = ca_coord[0];
            dir_set[ca_coord[0]] = true;
          } else {
            rotation[i].associate_direction[u_index] = ca_coord[1];
            dir_set[ca_coord[1]] = true;
          }
        }
      } break;
      default:
        return GOMA_ERROR;
      }

      for (int k = 0; k < 3; k++) {
        if (dir_set[k]) {
          associated_directions++;
        }
      }

      switch (rotation[i].type) {
      case GOMA_ROTATION_SURFACE:
        if (associated_directions != 1) {
          fprintf(stderr, "surface rotation isn't associated with 1 direction (%d directions)",
                  associated_directions);
          return GOMA_ERROR;
        }
        break;
      case GOMA_ROTATION_SINGLE:
        if (associated_directions != 1) {
          fprintf(stderr, "single rotation isn't associated with 1 direction (%d directions)",
                  associated_directions);
          return GOMA_ERROR;
        }
        break;
      case GOMA_ROTATION_EDGE:
        if (associated_directions != 2) {
          fprintf(stderr, "edge rotation isn't associated with 2 directions (%d directions)",
                  associated_directions);
          return GOMA_ERROR;
        }
        break;
      case GOMA_ROTATION_CORNER:
        if (associated_directions != 3) {
          fprintf(stderr, "corner rotation isn't associated with 3 directions (%d directions)",
                  associated_directions);
          return GOMA_ERROR;
        }
        break;
      }
    }
  }
  return GOMA_SUCCESS;
}

goma_error set_average_normals_and_tangents(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotation[i].n_normals > 0) {
      switch (rotation[i].type) {
      case GOMA_ROTATION_SURFACE:
      case GOMA_ROTATION_SINGLE: {
        gds_vector *average_normal = gds_vector_alloc(3);
        gds_vector_zero(average_normal);
        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          gds_vector_add(average_normal, rotation[i].normals[u_index]);
        }
        double invcount = 1.0 / rotation[i].n_normals;
        gds_vector_scale(average_normal, invcount);
        gds_vector_normalize(average_normal);
        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          gds_vector_copy(rotation[i].average_normals[u_index], average_normal);
        }
        gds_vector_free(average_normal);
      } break;
      case GOMA_ROTATION_EDGE: {
        for (unsigned int ca = 0; ca < 2; ca++) {
          gds_vector *average_normal = gds_vector_alloc(3);
          gds_vector_zero(average_normal);
          int count = 0;
          for (int v_index = 0; v_index < rotation[i].n_normals; v_index++) {
            if (rotation[i].associate_critical_normal[v_index] == ca) {
              gds_vector_add(average_normal, rotation[i].normals[v_index]);
              count++;
            }
          }
          gds_vector_scale(average_normal, 1.0 / (double)count);
          gds_vector_normalize(average_normal);
          for (int v_index = 0; v_index < rotation[i].n_normals; v_index++) {
            if (rotation[i].associate_critical_normal[v_index] == ca) {
              gds_vector_copy(rotation[i].average_normals[v_index], average_normal);
            }
          }
          gds_vector_free(average_normal);
        }
      } break;
      case GOMA_ROTATION_CORNER: {
        if (rotation[i].n_normals < 3) {
          EH(GOMA_ERROR, "Corner case incorrect number of normals");
          return GOMA_ERROR;
        } else {
          for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
            gds_vector_copy(rotation[i].average_normals[u_index], rotation[i].normals[u_index]);
          }
        }
      } break;
      }

      for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
        gds_vector_normalize(rotation[i].average_normals[u_index]);
        int tcross_card = -1;
        double min_card = 1e18;
        for (unsigned int card = 0; card < DIM; card++) {
          double dot = fabs(gds_vector_get(rotation[i].average_normals[u_index], card));
          if (dot < min_card) {
            min_card = dot;
            tcross_card = (int)card;
          }
        }

        assert(tcross_card != -1);
        gds_vector *tangent1 = rotation[i].tangent1s[u_index];
        gds_vector_zero(tangent1);
        gds_vector_set(tangent1, (unsigned int)tcross_card, 1.0);
        gds_vector *seed = gds_vector_alloc(3);
        gds_vector_copy(seed, tangent1);
        /* First tanget defined by seed vector gives
         *  t1 = tseed - n(tseed DOT n) normalized
         */
        for (unsigned int p = 0; p < DIM; p++) {
          for (unsigned int q = 0; q < DIM; q++) {
            double tangent_p = gds_vector_get(tangent1, p);
            tangent_p -= gds_vector_get(rotation[i].average_normals[u_index], p) *
                         gds_vector_get(rotation[i].average_normals[u_index], q) *
                         gds_vector_get(seed, q);
            gds_vector_set(tangent1, p, tangent_p);
          }
        }
        gds_vector_normalize(tangent1);
        gds_vector *tangent2 = rotation[i].tangent2s[u_index];
        gds_vector_cross(rotation[i].average_normals[u_index], tangent1, tangent2);
        gds_vector_free(seed);
      }
    }
  }

  return GOMA_SUCCESS;
}
goma_error set_rotated_coordinate_system(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotation[i].n_normals > 0) {
      switch (rotation[i].type) {
      case GOMA_ROTATION_SURFACE:
      case GOMA_ROTATION_SINGLE: {
        // lets just brute force all permutations see best association
        unsigned int best_n = 0, best_t1 = 0, best_t2 = 0;
        double maximum = 0;
        for (unsigned int n_index = 0; n_index < DIM; n_index++) {
          for (unsigned int t1_index = 0; t1_index < DIM; t1_index++) {
            if (n_index == t1_index)
              continue;
            for (unsigned int t2_index = 0; t2_index < DIM; t2_index++) {
              if (t2_index == n_index || t2_index == t1_index)
                continue;
              double idot = fabs(gds_vector_get(rotation[i].average_normals[0], n_index));
              double jdot = fabs(gds_vector_get(rotation[i].tangent1s[0], t1_index));
              double kdot = fabs(gds_vector_get(rotation[i].tangent2s[0], t2_index));
              double sum = idot + jdot + kdot;
              if (sum > maximum) {
                maximum = sum;
                best_n = n_index;
                best_t1 = t1_index;
                best_t2 = t2_index;
              }
            }
          }
        }

        assert(best_n != best_t1);
        assert(best_t1 != best_t2);
        gds_vector_copy(rotation[i].rotated_coord[best_n], rotation[i].average_normals[0]);
        gds_vector_copy(rotation[i].rotated_coord[best_t1], rotation[i].tangent1s[0]);
        gds_vector_copy(rotation[i].rotated_coord[best_t2], rotation[i].tangent2s[0]);
      } break;
      case GOMA_ROTATION_EDGE: {
        gds_vector *n1, *n2, *cross, *new_n1, *new_n2;
        n1 = gds_vector_alloc(3);
        n2 = gds_vector_alloc(3);
        cross = gds_vector_alloc(3);
        new_n1 = gds_vector_alloc(3);
        new_n2 = gds_vector_alloc(3);
        gds_vector_copy(n1, rotation[i].average_normals[0]);
        unsigned int n1_card = rotation[i].associate_direction[0];
        unsigned int n2_card = n1_card;
        for (int u_index = 1; u_index < rotation[i].n_normals; u_index++) {
          if (rotation[i].associate_direction[u_index] != n1_card) {
            gds_vector_copy(n2, rotation[i].average_normals[u_index]);
            n2_card = rotation[i].associate_direction[u_index];
            break;
          }
        }

        int o_card = -1;
        bool cards[DIM] = {false, false, false};
        cards[n1_card] = true;
        cards[n2_card] = true;
        for (int i = 0; i < DIM; i++) {
          if (!cards[i]) {
            o_card = i;
            break;
          }
        }
        assert((unsigned int)o_card != n1_card);
        assert((unsigned int)o_card != n2_card);
        assert(n1_card != n2_card);
        if (o_card == -1 || n1_card == n2_card) {
          return GOMA_ERROR;
        }

        gds_vector_cross(n1, n2, cross);
        gds_vector_normalize(cross);

        // find angle to offset rotation
        dbl dot = gds_vector_dot(n1, n2);
        const dbl ninety = 90 * M_PI / 180;

        dbl angle = acos(dot);
        dbl shift = 0.5 * (angle - ninety);
        gds_vector_rotate_around_vector(new_n1, n1, cross, shift);
        gds_vector_rotate_around_vector(new_n2, n2, cross, -shift);
        gds_vector_normalize(new_n1);
        gds_vector_normalize(new_n2);

        gds_vector_copy(rotation[i].rotated_coord[o_card], cross);
        gds_vector_copy(rotation[i].rotated_coord[n1_card], new_n1);
        gds_vector_copy(rotation[i].rotated_coord[n2_card], new_n2);
        gds_vector_free(n1);
        gds_vector_free(n2);
        gds_vector_free(cross);
        gds_vector_free(new_n1);
        gds_vector_free(new_n2);

      } break;
      case GOMA_ROTATION_CORNER: {
        const gds_vector *n1;
        const gds_vector *n2;
        const gds_vector *n3;
        gds_vector *cross = gds_vector_alloc(3);
        gds_vector *new_n1 = gds_vector_alloc(3);
        gds_vector *new_n2 = gds_vector_alloc(3);
        gds_vector *new_n3 = gds_vector_alloc(3);
        n1 = rotation[i].average_normals[0];
        n2 = rotation[i].average_normals[1];
        n3 = rotation[i].average_normals[2];
        bool found[DIM] = {true, false, false};
        unsigned int n1_card = rotation[i].associate_direction[0];
        unsigned int n2_card = n1_card;
        for (int j = 1; j < rotation[i].n_normals; j++) {
          if (rotation[i].associate_direction[j] != n1_card) {
            n2 = rotation[i].average_normals[j];
            n2_card = rotation[i].associate_direction[j];
            found[1] = true;
            break;
          }
        }
        unsigned int n3_card = n1_card;
        for (int j = 2; j < rotation[i].n_normals; j++) {
          if ((rotation[i].associate_direction[j] != n1_card) &&
              (rotation[i].associate_direction[j] != n2_card)) {
            n3 = rotation[i].average_normals[j];
            n3_card = rotation[i].associate_direction[j];
            found[2] = true;
            break;
          }
        }

        if (!found[0] || !found[1] || !found[2]) {
          EH(GOMA_ERROR, "Could not associate 3 directions for corner node");
          return GOMA_ERROR;
        }

        gds_vector_cross(n1, n2, cross);
        gds_vector_normalize(cross);

        if (gds_vector_dot(cross, n3) < 0) {
          gds_vector_scale(cross, -1.0);
        }

        gds_vector_copy(new_n3, cross);
        gds_vector_cross(n1, new_n3, cross);
        gds_vector_normalize(cross);

        if (gds_vector_dot(cross, n2) < 0) {
          gds_vector_scale(cross, -1.0);
        }
        gds_vector_copy(new_n2, cross);

        gds_vector_copy(new_n1, n1);

        gds_vector_copy(rotation[i].rotated_coord[n1_card], new_n1);
        gds_vector_copy(rotation[i].rotated_coord[n2_card], new_n2);
        gds_vector_copy(rotation[i].rotated_coord[n3_card], new_n3);
        gds_vector_free(cross);
        gds_vector_free(new_n1);
        gds_vector_free(new_n2);
        gds_vector_free(new_n3);

      } break;
      }
    }
  }
  return GOMA_SUCCESS;
}

goma_error set_face_normal_association(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    switch (rotation[i].type) {
    case GOMA_ROTATION_SINGLE:
    case GOMA_ROTATION_SURFACE:
    case GOMA_ROTATION_CORNER:
      for (int n_index = 0; n_index < rotation[i].n_normals; n_index++) {
        unsigned int best = 0;
        double max_dot = 0;
        for (unsigned int cord = 0; cord < DIM; cord++) {
          double dot =
              gds_vector_dot(rotation[i].normals[n_index], rotation[i].rotated_coord[cord]);
          if (dot > max_dot) {
            max_dot = dot;
            best = cord;
          }
        }
        rotation[i].face_coordinate_association[n_index] = best;
      }
      break;
    case GOMA_ROTATION_EDGE: {
      double ca_max[2] = {0.0, 0.0};
      unsigned int ca_coord[2] = {0, 0};
      gds_vector *ca1 = rotation[i].average_normals[rotation[i].critical_normal_index[0]];
      gds_vector *ca2 = rotation[i].average_normals[rotation[i].critical_normal_index[1]];
      for (unsigned int cord = 0; cord < DIM; cord++) {
        double dot1 = fabs(gds_vector_dot(ca1, rotation[i].rotated_coord[cord]));
        double dot2 = fabs(gds_vector_dot(ca2, rotation[i].rotated_coord[cord]));
        if (dot1 > ca_max[0] && dot1 >= dot2) {
          ca_max[0] = dot1;
          ca_coord[0] = cord;
        } else if (dot2 > ca_max[1]) {
          ca_max[1] = dot2;
          ca_coord[1] = cord;
        }
      }

      if (ca_max[0] < 1e-12 || ca_max[1] < 1e-12 || ca_coord[0] == ca_coord[1]) {
        return GOMA_ERROR;
      }

      for (int n_index = 0; n_index < rotation[i].n_normals; n_index++) {
        unsigned int ca = rotation[i].associate_critical_normal[n_index];
        rotation[i].face_coordinate_association[n_index] = ca_coord[ca];
      }
    } break;
    default:
      return GOMA_ERROR;
    }
  }
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

