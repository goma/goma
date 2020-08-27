#include "bc/rotate_util.h"
#include "bc/rotate_coordinates.h"
#include "gds/gds_vector.h"
#include "mm_eh.h"
#include <assert.h>

#ifndef GOMA_ROTATIONS_CRITICAL_ANGLE
#define GOMA_ROTATIONS_CRITICAL_ANGLE (M_PI / 4.0) // 45 degrees
#endif

gds_vector *goma_get_average_normal(gds_vector **normals, int n_normals) {
  gds_vector *average = gds_vector_alloc(3);
  gds_vector_zero(average);
  for (int i = 0; i < n_normals; i++) {
    gds_vector_add(average, normals[i]);
  }
  gds_vector_normalize(average);
  return average;
}

goma_error
goma_surface_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]) {

  gds_vector *normal = goma_get_average_normal(normals, n_normals);

  // find furthest coordinate
  int worst = -1;
  double dot_min = DBL_MAX;

  for (int i = 0; i < 3; i++) {
    int dot = fabs(gds_vector_get(normal, i));
    if (dot < dot_min) {
      worst = i;
      dot = dot_min;
    }
  }
  assert(worst > -1 && worst < 3);

  gds_vector *tangent = gds_vector_alloc(3);
  gds_vector *seed = gds_vector_alloc(3);
  gds_vector_zero(seed);
  // Use the worst as a seed for the tangent
  gds_vector_set(seed, worst, 1.0);
  gds_vector_copy(tangent, seed);

  /* First tanget defined by seed vector gives
   *  tangent = seed - n(seed DOT n) normalized
   */
  for (unsigned int p = 0; p < DIM; p++) {
    for (unsigned int q = 0; q < DIM; q++) {
      double tangent_p = gds_vector_get(tangent, p);
      tangent_p -= gds_vector_get(normal, p) * gds_vector_get(normal, q) * gds_vector_get(seed, q);
      gds_vector_set(tangent, p, tangent_p);
    }
  }
  gds_vector_normalize(tangent);

  gds_vector *binormal = gds_vector_alloc(3);
  gds_vector_zero(binormal);
  // binormal is normal x tangent
  gds_vector_cross(normal, tangent, binormal);
  gds_vector_normalize(binormal);

  int permute_max = 0;
  double max_val = 0;
  // Find best coordinate direction
  for (int p = 0; p < 3; p++) {
    double val = gds_vector_get(normal, (0 + p) % 3) + gds_vector_get(tangent, (1 + p) % 3) +
                 gds_vector_get(binormal, (2 + p) % 3);
    if (val > max_val) {
      permute_max = p;
      max_val = val;
    }
  }

  // set best directions
  gds_vector_copy(coord[(0 + permute_max) % 3], normal);
  gds_vector_copy(coord[(1 + permute_max) % 3], tangent);
  gds_vector_copy(coord[(2 + permute_max) % 3], binormal);

  gds_vector_free(seed);
  gds_vector_free(normal);
  gds_vector_free(tangent);
  gds_vector_free(binormal);
  return GOMA_SUCCESS;
}

goma_error goma_edge_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]) {

  return GOMA_SUCCESS;
}

goma_error
goma_corner_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]) {
  int critical_angle[3];
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index], normals[v_index]));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = fabs(gds_vector_dot(normals[u_index], normals[w_index]));
          double vw_dot = fabs(gds_vector_dot(normals[v_index], normals[w_index]));
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
  gds_vector *first = gds_vector_alloc(3);
  gds_vector *second = gds_vector_alloc(3);
  gds_vector *third = gds_vector_alloc(3);

  gds_vector_copy(first, normals[critical_angle[0]]);
  gds_vector_normalize(first);
  gds_vector_cross(first, normals[critical_angle[1]], second);
  gds_vector_normalize(second);
  gds_vector_cross(first, second, third);
  gds_vector_normalize(third);

  int permute_max = 0;
  double max_val = 0;
  // Find best coordinate direction
  for (int p = 0; p < 3; p++) {
    double val = gds_vector_get(first, (0 + p) % 3) + gds_vector_get(second, (1 + p) % 3) +
                 gds_vector_get(third, (2 + p) % 3);
    if (val > max_val) {
      permute_max = p;
      max_val = val;
    }
  }

  // set best directions
  gds_vector_copy(coord[(0 + permute_max) % 3], first);
  gds_vector_copy(coord[(1 + permute_max) % 3], second);
  gds_vector_copy(coord[(2 + permute_max) % 3], third);

  gds_vector_free(first);
  gds_vector_free(second);
  gds_vector_free(third);
}

  return GOMA_SUCCESS;
}

bool goma_check_normals_within_critical_angle(gds_vector **normals, int n_normals) {
  double min_dot = DBL_MAX;
  for (int i = 0; i < n_normals; i++) {
    for (int j = 0; j < n_normals; j++) {
      double dot = fabs(gds_vector_dot(normals[i], normals[j]));
      if (dot < min_dot) {
        min_dot = dot;
      }
    }
  }

  // True if all are within a critical angle
  return min_dot > cos(GOMA_ROTATION_CRITICAL_ANGLE);
}

bool goma_check_corner_rotation_case(gds_vector **normals, int n_normals) {
  if (n_normals < 3) {
    return false;
  }
  // If two normals are larger than critical angle then there should be a third normal larger
  // than a critical angle to both those to be a corner type
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index], normals[v_index]));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = fabs(gds_vector_dot(normals[u_index], normals[w_index]));
          double vw_dot = fabs(gds_vector_dot(normals[v_index], normals[w_index]));
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

bool goma_check_edge_rotation_case(gds_vector **normals, int n_normals) {
  if (n_normals < 2) {
    return false;
  }

  int crit_normal[2];
  bool found = false;
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = fabs(gds_vector_dot(normals[u_index], normals[v_index]));
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
    bool match_one = fabs(gds_vector_dot(normals[u_index], normals[crit_normal[0]])) <
                         cos(GOMA_ROTATION_CRITICAL_ANGLE) ||
                     fabs(gds_vector_dot(normals[u_index], normals[crit_normal[1]])) <
                         cos(GOMA_ROTATION_CRITICAL_ANGLE);
    if (!match_one) {
      return false;
    }
  }

  return true;
}

goma_error
goma_best_coordinate_system_3D(gds_vector **normals, int n_normals, gds_vector *coord[3]) {

  // normalize just to be sure
  for (int i = 0; i < n_normals; i++) {
    gds_vector_normalize(normals[i]);
  }

  // check if all are within a critical angle
  // then we have a surface case
  if (goma_check_normals_within_critical_angle(normals, n_normals)) {
    return goma_surface_coordinate_system(normals, n_normals, coord);
  } else if (goma_check_corner_rotation_case(normals, n_normals)) {
    return goma_corner_coordinate_system(normals, n_normals, coord);
  } else if (goma_check_edge_rotation_case(normals, n_normals)) {
    return goma_edge_coordinate_system(normals, n_normals, coord);
  } else {
    return GOMA_ERROR;
  }

  return GOMA_SUCCESS;
}
