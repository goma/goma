#include "bc/rotate_util.h"
#include "bc/rotate_coordinates.h"
#include "gds/gds_vector.h"
#include "mm_eh.h"
#include "std.h"
#include "util/goma_normal.h"
#include <assert.h>

static void goma_normal_assign_best_direction(goma_normal *normal,
    goma_normal *tangent,
    goma_normal *binormal,
    goma_normal *coord[3]) {
  double max_val = 0;
  int n_max = -1;
  int t_max = -1;
  int b_max = -1;
  double val = 0;
  if ((val = (fabs(gds_vector_get(normal->normal, 0)) + fabs(gds_vector_get(tangent->normal, 1)) +
              fabs(gds_vector_get(binormal->normal, 2)))) > max_val) {
    max_val = val;
    n_max = 0;
    t_max = 1;
    b_max = 2;
  } 

  if ((val = (fabs(gds_vector_get(normal->normal, 2)) +
                     fabs(gds_vector_get(tangent->normal, 1)) +
                     fabs(gds_vector_get(binormal->normal, 0)))) > max_val) {
    max_val = val;
    n_max = 2;
    t_max = 1;
    b_max = 0;
  } 

  if ((val = (fabs(gds_vector_get(normal->normal, 2)) +
                     fabs(gds_vector_get(tangent->normal, 0)) +
                     fabs(gds_vector_get(binormal->normal, 1)))) > max_val) {
    max_val = val;
    n_max = 2;
    t_max = 0;
    b_max = 1;
  } 

  if ((val = (fabs(gds_vector_get(normal->normal, 0)) +
                     fabs(gds_vector_get(tangent->normal, 2)) +
                     fabs(gds_vector_get(binormal->normal, 1)))) > max_val) {
    max_val = val;
    n_max = 0;
    t_max = 2;
    b_max = 1;
  } 

  if ((val = (fabs(gds_vector_get(normal->normal, 1)) +
                     fabs(gds_vector_get(tangent->normal, 0)) +
                     fabs(gds_vector_get(binormal->normal, 2)))) > max_val) {
    max_val = val;
    n_max = 1;
    t_max = 0;
    b_max = 2;
  } 

  if ((val = (fabs(gds_vector_get(normal->normal, 1)) +
                     fabs(gds_vector_get(tangent->normal, 2)) +
                     fabs(gds_vector_get(binormal->normal, 0)))) > max_val) {
    max_val = val;
    n_max = 1;
    t_max = 2;
    b_max = 0;
  }
  // Find best coordinate direction
  //for (int p = 0; p < 3; p++) {
  //  double val = fabs(gds_vector_get(normal->normal, (0 + p) % 3)) +
  //               fabs(gds_vector_get(tangent->normal, (1 + p) % 3)) +
  //               fabs(gds_vector_get(binormal->normal, (2 + p) % 3));
  //  if (val > max_val) {
  //    permute_max = p;
  //    max_val = val;
  //  }
  //}
  if (n_max == -1 || t_max == -1 || b_max == -1) {
    EH(GOMA_ERROR, "could not associate normals to coord");
  }

  // set best directions
  goma_normal_copy(coord[n_max], normal);
  goma_normal_copy(coord[t_max], tangent);
  goma_normal_copy(coord[b_max], binormal);
  
  // recompute binormal component to ensure righthanded
  // this lets us keep normal directions for surface and corners
  if (b_max == 0) {
    goma_normal_cross(coord[1], coord[2], coord[0]);
    goma_normal_normalize(coord[0]);
  } else if (b_max == 1) {
    goma_normal_cross(coord[2], coord[0], coord[1]);
    goma_normal_normalize(coord[1]);
  } else {
    goma_normal_cross(coord[0], coord[1], coord[2]);
    goma_normal_normalize(coord[2]);
  }
}

goma_error
goma_surface_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]) {

  goma_normal *normal = goma_get_average_normal(normals, n_normals);

  // find furthest coordinate
  int worst = -1;
  double dot_min = DBL_MAX;

  for (int i = 0; i < 3; i++) {
    goma_normal_val val = goma_normal_get(normal, i);
    goma_normal_val dot = fabs_goma_normal_val(&val);
    if (dot.val < dot_min) {
      worst = i;
      dot_min = dot.val;
    }
  }
  assert(worst > -1 && worst < 3);

  goma_normal *tangent = goma_normal_alloc(3);
  goma_normal *seed = goma_normal_alloc(3);
  goma_normal_zero(seed);
  // Use the worst as a seed for the tangent
  goma_normal_set_constant(seed, worst, 1.0);
  goma_normal_copy(tangent, seed);

  /* First tanget defined by seed vector gives
   *  tangent = seed - n(seed DOT n) normalized
   */

  goma_normal_val tmpval = goma_normal_dot(seed, normal);
  goma_normal *tmp_normal = goma_normal_alloc(3);
  goma_normal_copy(tmp_normal, normal);
  goma_normal_scale(tmp_normal, &tmpval);
  goma_normal_sub(tangent, tmp_normal);
  goma_normal_free(tmp_normal);

  // for (unsigned int p = 0; p < DIM; p++) {
  //  for (unsigned int q = 0; q < DIM; q++) {
  //    goma_normal_val tangent_p = goma_normal_get(tangent, p);
  //    tangent_p -= goma_normal_get(normal, p) * goma_normal_get(normal, q) * goma_normal_get(seed,
  //    q); goma_normal_set(tangent, p, tangent_p);
  //  }
  //}
  goma_normal_normalize(tangent);

  goma_normal *binormal = goma_normal_alloc(3);
  goma_normal_zero(binormal);
  // binormal is normal x tangent
  goma_normal_cross(normal, tangent, binormal);
  goma_normal_normalize(binormal);

  goma_normal_assign_best_direction(normal, tangent, binormal, coord);

  goma_normal_free(seed);
  goma_normal_free(normal);
  goma_normal_free(tangent);
  goma_normal_free(binormal);
  return GOMA_SUCCESS;
}

static goma_error get_average_edge_normals(goma_normal **normals,
                                           int n_normals,
                                           goma_normal **avg_1,
                                           goma_normal **avg_2) {
  bool found = false;
  int crit_normal[2];
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        crit_normal[0] = u_index;
        crit_normal[1] = v_index;
        found = true;
      }
    }
  }
  if (!found) {
    return GOMA_ERROR;
  }

  goma_normal **n1 = malloc(sizeof(goma_normal *) * n_normals);
  int n_n1 = 1;
  goma_normal **n2 = malloc(sizeof(goma_normal *) * n_normals);
  int n_n2 = 1;
  n1[0] = normals[crit_normal[0]];
  n2[0] = normals[crit_normal[1]];

  for (int i = 0; i < n_normals; i++) {
    if (i == crit_normal[0] || i == crit_normal[1]) {
      continue;
    }

    if (fabs(gds_vector_dot(normals[i]->normal, n1[0]->normal)) >=
        cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
      n1[n_n1++] = normals[i];
    } else if (fabs(gds_vector_dot(normals[i]->normal, n2[0]->normal)) >=
               cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
      n2[n_n2++] = normals[i];
    } else {
      free(n1);
      free(n2);
      return GOMA_ERROR;
    }
  }

  *avg_1 = goma_get_average_normal(n1, n_n1);
  *avg_2 = goma_get_average_normal(n2, n_n2);

  free(n1);
  free(n2);
  return GOMA_SUCCESS;
}

goma_error
goma_edge_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]) {

  goma_normal *n1 = NULL;
  goma_normal *n2 = NULL;
  goma_normal *cross = NULL;
  goma_normal *new_n1 = NULL;
  goma_normal *new_n2 = NULL;
  goma_error err = get_average_edge_normals(normals, n_normals, &n1, &n2);
  if (err != GOMA_SUCCESS) {
    return err;
  }

  cross = goma_normal_alloc(3);
  goma_normal_cross(n1, n2, cross);
  goma_normal_normalize(cross);

  // find angle to offset rotation
  goma_normal_val dot = goma_normal_dot(n1, n2);
  const double ninety = 90 * M_PI / 180;

  new_n1 = goma_normal_alloc(3);
  new_n2 = goma_normal_alloc(3);
  goma_normal_val angle = acos_goma_normal_val(&dot);
  goma_normal_val subval = sub_goma_normal_val(&angle, ninety);
  const double half = 0.5;
  goma_normal_val shift = scale_goma_normal_val(&subval, half);
  goma_normal_val minus_shift = scale_goma_normal_val(&subval, -half);
  goma_normal_rotate_around_vector(new_n1, n1, cross, shift);
  goma_normal_rotate_around_vector(new_n2, n2, cross, minus_shift);
  goma_normal_normalize(new_n1);
  goma_normal_normalize(new_n2);

  goma_normal_cross(new_n1, new_n2, cross);
  goma_normal_normalize(cross);

  goma_normal_assign_best_direction(new_n1, new_n2, cross, coord);

  goma_normal_free(new_n1);
  goma_normal_free(new_n2);
  goma_normal_free(cross);
  goma_normal_free(n1);
  goma_normal_free(n2);

  return GOMA_SUCCESS;
}

goma_error
goma_corner_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]) {
  int critical_angle[3];
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[w_index]->normal));
          double vw_dot = fabs(gds_vector_dot(normals[v_index]->normal, normals[w_index]->normal));
          if (uw_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE) &&
              vw_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
            critical_angle[0] = u_index;
            critical_angle[1] = v_index;
            critical_angle[2] = w_index;
            goto corner_coord_critical_found;
          }
        }
      }
    }
  }

  return GOMA_ERROR;

corner_coord_critical_found : {
  // find closest to x for first
  int first_crit = -1;
  double max_dot = 0;
  for (int i = 0; i < 3; i++) {
    double dot = fabs(gds_vector_get(normals[critical_angle[i]]->normal, 0));
    if (dot > max_dot) {
      first_crit = i;
      max_dot = dot;
    }
  }

  if (first_crit == -1) {
    return GOMA_ERROR;
  }

  // closest to y for second
  int second_crit = -1;
  max_dot = 0;
  for (int i = 0; i < 3; i++) {
    double dot = fabs(gds_vector_get(normals[critical_angle[i]]->normal, 1));
    if (dot > max_dot) {
      second_crit = i;
      max_dot = dot;
    }
  }

  if (second_crit == -1) {
    return GOMA_ERROR;
  }
  
  goma_normal *first = goma_normal_alloc(3);
  goma_normal *second = goma_normal_alloc(3);
  goma_normal *third = goma_normal_alloc(3);

  goma_normal_copy(first, normals[critical_angle[first_crit]]);
  goma_normal_normalize(first);
  goma_normal_cross(first, normals[critical_angle[second_crit]], second);
  goma_normal_normalize(second);
  goma_normal_cross(first, second, third);
  goma_normal_normalize(third);

  goma_normal_assign_best_direction(first, second, third, coord);

  goma_normal_free(first);
  goma_normal_free(second);
  goma_normal_free(third);
}

  return GOMA_SUCCESS;
}

bool goma_check_normals_within_critical_angle(goma_normal **normals, int n_normals) {
  double min_dot = DBL_MAX;
  for (int i = 0; i < n_normals; i++) {
    for (int j = 0; j < n_normals; j++) {
      double dot = fabs(gds_vector_dot(normals[i]->normal, normals[j]->normal));
      if (dot < min_dot) {
        min_dot = dot;
      }
    }
  }

  // True if all are within a critical angle
  return min_dot > cos(GOMA_ROTATION_CRITICAL_ANGLE);
}

bool goma_check_corner_rotation_case(goma_normal **normals, int n_normals) {
  if (n_normals < 3) {
    return false;
  }
  // If two normals are larger than critical angle then there should be a third normal larger
  // than a critical angle to both those to be a corner type
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[w_index]->normal));
          double vw_dot = fabs(gds_vector_dot(normals[v_index]->normal, normals[w_index]->normal));
          if (uw_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE) &&
              vw_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool goma_check_edge_rotation_case(goma_normal **normals, int n_normals) {
  if (n_normals < 2) {
    return false;
  }

  int crit_normal[2];
  bool found = false;
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        crit_normal[0] = u_index;
        crit_normal[1] = v_index;
        found = true;
      }
    }
  }

  if (!found) {
    return false;
  }

  // make sure all angles are near one of the critical angles
  for (int u_index = 0; u_index < n_normals; u_index++) {
    bool match_one =
        fabs(gds_vector_dot(normals[u_index]->normal, normals[crit_normal[0]]->normal)) >=
            cos(GOMA_ROTATION_CRITICAL_ANGLE) ||
        fabs(gds_vector_dot(normals[u_index]->normal, normals[crit_normal[1]]->normal)) >=
            cos(GOMA_ROTATION_CRITICAL_ANGLE);
    if (!match_one) {
      return false;
    }
  }

  return true;
}

goma_error
goma_best_coordinate_system_3D(goma_normal **normals, int n_normals, goma_normal *coord[3]) {

  // normalize just to be sure
  for (int i = 0; i < n_normals; i++) {
    goma_normal_normalize(normals[i]);
  }

  // check if all are within a critical angle
  // then we have a surface case
  if (goma_check_normals_within_critical_angle(normals, n_normals)) {
    return goma_surface_coordinate_system(normals, n_normals, coord);
  } else if (goma_check_corner_rotation_case(normals, n_normals)) {
    return goma_corner_coordinate_system(normals, n_normals, coord);
  } else if (goma_check_edge_rotation_case(normals, n_normals)) {
    return goma_edge_coordinate_system(normals, n_normals, coord);
  }

  return GOMA_ERROR;
}
