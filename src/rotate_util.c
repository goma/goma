#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rotate_util.h"

#include "goma.h"

static const double critical_angle_radians = GOMA_ROTATION_CRITICAL_ANGLE * M_PI / 180;

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

goma_error setup_bc_is_rotated_list(struct Boundary_Condition *bc_types, int num_bc,
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
  *rotations = calloc((size_t)exo->num_nodes, sizeof(goma_rotation_node_s));
  return GOMA_SUCCESS;
}

goma_error setup_rotated_bc_nodes(Exo_DB *exo, struct Boundary_Condition *bc_types, int num_bc,
                                  double *x) {
  goma_rotation_node_s *rotations;
  bool *bc_is_rotated;
  goma_error error;
  error = setup_bc_is_rotated_list(bc_types, num_bc, &bc_is_rotated);
  EH(error, "setup_bc_rotate_list");
  error = allocate_rotations(exo, &rotations);
  EH(error, "allocate_rotations");

  //  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
  //    if (bc_is_rotated[bc_index]) {
  //      int ss_index = bc_types[bc_index].Set_Index;
  //      printf("SS %d has %d nodes and is rotated\n", exo->ss_id[ss_index], exo->ss_node_len);
  //      for (int eindex = exo->ss_elem_index[ss_index]; eindex < exo->ss_elem_index[ss_index + 1];
  //           eindex++) {
  //        printf("has elem %d\n", exo->ss_elem_list[eindex]);
  //      }

  //      for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
  //        int elem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

  //        printf("elem %d side %d\n", elem, exo->ss_side_list[e]);
  //        for (int k = exo->ss_node_side_index[ss_index][e];
  //             k < exo->ss_node_side_index[ss_index][e + 1]; k++) {
  //          int inode = exo->ss_node_list[ss_index][k];
  //          printf("node %d\n", inode);
  //        }
  //      }
  //    }
  //  }

  int err;
  for (int bc_index = 0; bc_index < num_bc; bc_index++) {
    if (bc_is_rotated[bc_index]) {
      int ss_index = bc_types[bc_index].Set_Index;
      printf("SS %d has %d nodes and is rotated\n", exo->ss_id[ss_index], exo->ss_node_len);
      for (int eindex = exo->ss_elem_index[ss_index]; eindex < exo->ss_elem_index[ss_index + 1];
           eindex++) {
        printf("has elem %d\n", exo->ss_elem_list[eindex]);
      }

      for (int e = 0; e < exo->ss_num_sides[ss_index]; e++) {
        int ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index] + e];

        if (x == x_static) /* be the least disruptive possible */
        {
          err = load_elem_dofptr(ielem, exo, x_static, x_old_static, xdot_static, xdot_old_static,
                                 x_static, 0);
        } else {
          err = load_elem_dofptr(ielem, exo, x, x, x, x, x, 0);
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
          assert(rotations[I].n_normals <= GOMA_MAX_NORMALS_PER_NODE);
          rotations[I].normals[n_index].e[0] = fv->snormal[0];
          rotations[I].normals[n_index].e[1] = fv->snormal[1];
          rotations[I].normals[n_index].e[2] = fv->snormal[2];
          rotations[I].element[n_index] = ielem;
          rotations[I].face[n_index] = id_side;
          printf("node %d, now has %d normals\n", I, rotations[I].n_normals);
        }
      }
    }
  }

  error = set_rotation_types(exo, rotations);
  EH(error, "set_rotation_types");
  error = associate_directions(exo, rotations);
  EH(error, "associate_directions");
  error = set_average_normals_and_tangents(exo, rotations);
  EH(error, "set_average_normals_and_tangents");
  error = set_rotated_coordinate_system(exo, rotations);
  EH(error, "set_rotated_coordinate_system");
  FILE *normalscsv;
  normalscsv = fopen("normals.csv", "w");

  fprintf(normalscsv, "n_normals,type,node,x,y,z,nx,ny,nz,dir,ax,ay,az,t1x,t1y,t1z,t2x,t2y,t2z,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z\n");
  for (int i = 0; i < exo->num_nodes; i++) {
    for (int j = 0; j < rotations[i].n_normals; j++) {
      fprintf(normalscsv, "%d,%d,%d,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
              rotations[i].n_normals, rotations[i].type, i, exo->x_coord[i], exo->y_coord[i],
              exo->z_coord[i], rotations[i].normals[j].e[0], rotations[i].normals[j].e[1],
              rotations[i].normals[j].e[2], rotations[i].associate_direction[j],
              rotations[i].average_normals[j].e[0], rotations[i].average_normals[j].e[1],
              rotations[i].average_normals[j].e[2], rotations[i].tangent1s[j].e[0],
              rotations[i].tangent1s[j].e[1], rotations[i].tangent1s[j].e[2],
              rotations[i].tangent2s[j].e[0], rotations[i].tangent2s[j].e[1],
              rotations[i].tangent2s[j].e[2],
              rotations[i].rotated_coord[0].e[0],
              rotations[i].rotated_coord[0].e[1],
              rotations[i].rotated_coord[0].e[2],
              rotations[i].rotated_coord[1].e[0],
              rotations[i].rotated_coord[1].e[1],
              rotations[i].rotated_coord[1].e[2],
              rotations[i].rotated_coord[2].e[0],
              rotations[i].rotated_coord[2].e[1],
              rotations[i].rotated_coord[2].e[2]
              );
    }
  }

  goma_automatic_rotations.automatic_rotations = true;
  goma_automatic_rotations.rotation_nodes = rotations;
  fclose(normalscsv);
  free(bc_is_rotated);
  return GOMA_SUCCESS;
}

double goma_rotation_vector_dot(const goma_rotation_vector_s *v, const goma_rotation_vector_s *u) {
  return v->e[0] * u->e[0] + v->e[1] * u->e[1] + v->e[2] * u->e[2];
}

/**
 * @brief goma_rotation_node_type_surface
 * @param node the node we are looking at rotating
 * @return rotation at this node is type surface
 */
bool goma_rotation_node_type_surface(const goma_rotation_node_s *node) {
  double min_dot = 1e300;
  for (int u_index = 0; u_index < node->n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < node->n_normals; v_index++) {
      double dot = goma_rotation_vector_dot(&(node->normals[u_index]), &(node->normals[v_index]));
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
      double uv_dot =
          goma_rotation_vector_dot(&(node->normals[u_index]), &(node->normals[v_index]));
      if (uv_dot < cos(critical_angle_radians)) {
        for (int w_index = 0; w_index < node->n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index)
            continue;
          double uw_dot =
              goma_rotation_vector_dot(&(node->normals[u_index]), &(node->normals[w_index]));
          double vw_dot =
              goma_rotation_vector_dot(&(node->normals[v_index]), &(node->normals[w_index]));
          if (uw_dot < cos(critical_angle_radians) && vw_dot < cos(critical_angle_radians)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

void goma_rotation_make_unit_vector(goma_rotation_vector_s *v) {
  double invmag = 1.0 / sqrt(goma_rotation_vector_dot(v, v));
  v->e[0] *= invmag;
  v->e[1] *= invmag;
  v->e[2] *= invmag;
}

goma_rotation_vector_s goma_rotation_unit_average_vectors(const goma_rotation_vector_s *v,
                                                          const goma_rotation_vector_s *u) {
  goma_rotation_vector_s average;
  average.e[0] = 0.5 * (v->e[0] + u->e[0]);
  average.e[1] = 0.5 * (v->e[1] + u->e[1]);
  average.e[2] = 0.5 * (v->e[2] + u->e[2]);
  goma_rotation_make_unit_vector(&average);
  return average;
}

bool goma_rotation_node_type_edge(const goma_rotation_node_s *node) {
  if (node->n_normals < 2) {
    return false;
  }

  goma_rotation_vector_s critical_normal[2];

  critical_normal[0] = node->normals[0];
  // find another critical normal
  bool found = false;
  for (int u_index = 1; u_index < node->n_normals; u_index++) {
    if (goma_rotation_vector_dot(&(critical_normal[0]), &(node->normals[u_index])) <
        cos(critical_angle_radians)) {
      found = true;
      critical_normal[1] = node->normals[u_index];
      break;
    }
  }

  if (!found) {
    return false;
  }

  // make sure we only have two critical normals
  for (int u_index = 0; u_index < node->n_normals; u_index++) {
    if ((goma_rotation_vector_dot(&(critical_normal[0]), &(node->normals[u_index])) <
         cos(critical_angle_radians)) &&
        (goma_rotation_vector_dot(&(critical_normal[1]), &(node->normals[u_index])) <
         cos(critical_angle_radians))) {
      return false;
    }
  }
  // if we made it here we should only have two angles
  return true;
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

      for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
        int best_dir = 0;
        double max_dot = 0;
        for (int card = 0; card < DIM; card++) {
          if (dir_set[card]) {
            // check if we are within critical angle to already set
            // if not pick next best
            bool skip = false;
            for (int v_index = 0; v_index < u_index; v_index++) {
              if (rotation[i].associate_direction[v_index] == card) {
                if (goma_rotation_vector_dot(&(rotation[i].normals[v_index]),
                                             &(rotation[i].normals[u_index])) <
                    cos(critical_angle_radians)) {
                  skip = true;
                  break;
                }
              }
            }
            if (skip)
              continue;
          }
          goma_rotation_vector_s cvect;
          cvect.e[0] = 0.0;
          cvect.e[1] = 0.0;
          cvect.e[2] = 0.0;
          cvect.e[card] = 1.0;
          double dot = fabs(goma_rotation_vector_dot(&cvect, &(rotation[i].normals[u_index])));
          if (dot > max_dot) {
            best_dir = card;
            max_dot = dot;
          }
        }
        rotation[i].associate_direction[u_index] = best_dir;
        for (int card = 0; card < DIM; card++) {
          if (best_dir != card) {
            rotation[i].tangent1_seeddir[u_index] = card;
          }
        }
        dir_set[best_dir] = true;
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
      for (int card = 0; card < DIM; card++) {
        goma_rotation_vector_s average_normal;
        average_normal.e[0] = 0;
        average_normal.e[1] = 0;
        average_normal.e[2] = 0;
        int n_card_normals = 0;
        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          if (rotation[i].associate_direction[u_index] == card) {
            average_normal.e[0] += rotation[i].normals[u_index].e[0];
            average_normal.e[1] += rotation[i].normals[u_index].e[1];
            average_normal.e[2] += rotation[i].normals[u_index].e[2];
            n_card_normals++;
          }
        }

        if (n_card_normals > 0) {
          average_normal.e[0] /= n_card_normals;
          average_normal.e[1] /= n_card_normals;
          average_normal.e[2] /= n_card_normals;
          goma_rotation_make_unit_vector(&average_normal);
          for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
            if (rotation[i].associate_direction[u_index] == card) {
              rotation[i].average_normals[u_index] = average_normal;
            }
          }
        }
      }
      for (int card = 0; card < DIM; card++) {
        for (int u_index = 0; u_index < rotation[i].n_normals; u_index++) {
          if (rotation[i].associate_direction[u_index] == card) {
            goma_rotation_vector_s tangent1;
            goma_rotation_vector_s seed;
            tangent1.e[0] = 0;
            tangent1.e[1] = 0;
            tangent1.e[2] = 0;
            tangent1.e[rotation[i].tangent1_seeddir[u_index]] = 1;
            seed = tangent1;
            for (int p = 0; p < DIM; p++) {
              for (int q = 0; q < DIM; q++) {
                tangent1.e[p] -= rotation[i].average_normals[u_index].e[p] *
                                 rotation[i].average_normals[u_index].e[q] * seed.e[q];
              }
            }
            goma_rotation_make_unit_vector(&tangent1);
            rotation[i].tangent1s[u_index] = tangent1;
            goma_rotation_vector_s tangent2;

            tangent2.e[0] = -(rotation[i].average_normals[u_index].e[1] * tangent1.e[2] -
                              rotation[i].average_normals[u_index].e[2] * tangent1.e[1]);
            tangent2.e[1] = -(rotation[i].average_normals[u_index].e[2] * tangent1.e[0] -
                              rotation[i].average_normals[u_index].e[0] * tangent1.e[2]);
            tangent2.e[2] = -(rotation[i].average_normals[u_index].e[0] * tangent1.e[1] -
                              rotation[i].average_normals[u_index].e[1] * tangent1.e[0]);
            rotation[i].tangent2s[u_index] = tangent2;
          }
        }
      }
    }
  }

  /* First tanget defined by seed vector gives
   *  t1 = tseed - n(tseed DOT n) normalized
   */

  return GOMA_SUCCESS;
}

void goma_rotation_cross_vectors(goma_rotation_vector_s *v, goma_rotation_vector_s *u, goma_rotation_vector_s *cross) {
  cross->e[0] = -(v->e[1] * u->e[2] - v->e[2] * u->e[1]);
  cross->e[1] = -(v->e[2] * u->e[0] - v->e[0] * u->e[2]);
  cross->e[2] = -(v->e[0] * u->e[1] - v->e[1] * u->e[0]);
}

void goma_rotation_scale_vector(goma_rotation_vector_s *v, dbl scale) {
v->e[0] *= scale;
v->e[1] *= scale;
v->e[2] *= scale;
}

void goma_rotation_sum_vectors(goma_rotation_vector_s *v, goma_rotation_vector_s *u, goma_rotation_vector_s *out) {
  out->e[0] = v->e[0] + u->e[0];
  out->e[1] = v->e[1] + u->e[1];
  out->e[2] = v->e[2] + u->e[2];
}


void goma_rotation_rotate_about_vector(goma_rotation_vector_s vector_to_rotate, goma_rotation_vector_s *vector_axis, goma_rotation_vector_s *rotated, dbl angle_radians) {
  goma_rotation_vector_s cross;
  goma_rotation_cross_vectors(vector_axis, &vector_to_rotate, &cross);

  goma_rotation_vector_s tmp1;
  tmp1 = vector_to_rotate;
  goma_rotation_scale_vector(&tmp1, cos(angle_radians));
  goma_rotation_vector_s tmp2 = cross;
  goma_rotation_scale_vector(&tmp2, sin(angle_radians));
  goma_rotation_vector_s tmp3 = *vector_axis;
  dbl dot = goma_rotation_vector_dot(&vector_to_rotate, vector_axis);
  goma_rotation_scale_vector(&tmp3, dot*(1-cos(angle_radians)));

  goma_rotation_sum_vectors(&tmp1, &tmp2, rotated);
  tmp1 = *rotated;
  goma_rotation_sum_vectors(&tmp1, &tmp3, rotated);
}


goma_error set_rotated_coordinate_system(Exo_DB *exo, goma_rotation_node_s *rotation) {
  for (int i = 0; i < exo->num_nodes; i++) {
    if (rotation[i].n_normals > 0) {
      switch (rotation[i].type) {
      case GOMA_ROTATION_SURFACE:
      case GOMA_ROTATION_SINGLE:
        rotation[i].rotated_coord[0] = rotation[i].average_normals[0];
        rotation[i].rotated_coord[1] = rotation[i].tangent1s[0];
        rotation[i].rotated_coord[2] = rotation[i].tangent2s[0];
        break;
      case GOMA_ROTATION_EDGE: {
        goma_rotation_vector_s n1, n2, cross, new_n1, new_n2;
        n1 = rotation[i].average_normals[0];
        int n1_card = rotation[i].associate_direction[0];
        int n2_card = n1_card;
        for (int u_index = 1; u_index < rotation[i].n_normals; u_index++) {
          if (rotation[i].associate_direction[u_index] != n1_card) {
            n2 = rotation[i].average_normals[u_index];
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
        assert(o_card != n1_card);
        assert(o_card != n2_card);
        assert(n1_card != n2_card);
        if (o_card == -1 || n1_card == n2_card) {
          return GOMA_ERROR;
        }

        cross.e[0] = -(n1.e[1] * n2.e[2] - n1.e[2] * n2.e[1]);
        cross.e[1] = -(n1.e[2] * n2.e[0] - n1.e[0] * n2.e[2]);
        cross.e[2] = -(n1.e[0] * n2.e[1] - n1.e[1] * n2.e[0]);

        goma_rotation_make_unit_vector(&cross);

        // find angle to offset rotation
        dbl dot = goma_rotation_vector_dot(&n1, &n2);
        const dbl ninety = 90*M_PI/180;

        dbl angle = acos(dot);
        dbl shift = 0.5*(angle - ninety);
        goma_rotation_rotate_about_vector(n1, &cross, &new_n1, shift);
        goma_rotation_rotate_about_vector(n2, &cross, &new_n2, -shift);

        goma_rotation_make_unit_vector(&cross);
        goma_rotation_make_unit_vector(&new_n1);
        goma_rotation_make_unit_vector(&new_n2);

        rotation[i].rotated_coord[o_card] = cross;
        rotation[i].rotated_coord[n1_card] = new_n1;
        rotation[i].rotated_coord[n2_card] = new_n2;

      } break;
      case GOMA_ROTATION_CORNER: {
        goma_rotation_vector_s n1, n2, n3, cross, new_n1, new_n2, new_n3;
        n1 = rotation[i].average_normals[0];
        int n1_card = rotation[i].associate_direction[0];
        n2 = rotation[i].average_normals[1];
        int n2_card = rotation[i].associate_direction[1];
        n3 = rotation[i].average_normals[2];
        int n3_card = rotation[i].associate_direction[2];

        cross.e[0] = -(n1.e[1] * n2.e[2] - n1.e[2] * n2.e[1]);
        cross.e[1] = -(n1.e[2] * n2.e[0] - n1.e[0] * n2.e[2]);
        cross.e[2] = -(n1.e[0] * n2.e[1] - n1.e[1] * n2.e[0]);

        goma_rotation_make_unit_vector(&cross);

        if (goma_rotation_vector_dot(&cross,&n3) < 0) {
          goma_rotation_scale_vector(&cross, -1.0);
        }

        new_n3 = cross;

        cross.e[0] = -(n1.e[1] * new_n3.e[2] - n1.e[2] * new_n3.e[1]);
        cross.e[1] = -(n1.e[2] * new_n3.e[0] - n1.e[0] * new_n3.e[2]);
        cross.e[2] = -(n1.e[0] * new_n3.e[1] - n1.e[1] * new_n3.e[0]);

        goma_rotation_make_unit_vector(&cross);

        if (goma_rotation_vector_dot(&cross,&n2) < 0) {
          goma_rotation_scale_vector(&cross, -1.0);
        }
        new_n2 = cross;

        new_n1 = n1;

        rotation[i].rotated_coord[n1_card] = new_n1;
        rotation[i].rotated_coord[n2_card] = new_n2;
        rotation[i].rotated_coord[n3_card] = new_n3;

      } break;
      }
    }
  }
  return GOMA_SUCCESS;
}
