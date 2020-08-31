#include "bc/rotate_util.h"
#include "bc/rotate_coordinates.h"
#include "gds/gds_vector.h"
#include "mm_eh.h"
#include "std.h"
#include <assert.h>

#ifndef GOMA_ROTATIONS_CRITICAL_ANGLE
#define GOMA_ROTATIONS_CRITICAL_ANGLE (M_PI / 4.0) // 45 degrees
#endif

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

  int permute_max = 0;
  double max_val = 0;
  // Find best coordinate direction
  for (int p = 0; p < 3; p++) {
    double val = gds_vector_get(normal->normal, (0 + p) % 3) + gds_vector_get(tangent->normal, (1 + p) % 3) +
                 gds_vector_get(binormal->normal, (2 + p) % 3);
    if (val > max_val) {
      permute_max = p;
      max_val = val;
    }
  }

  // set best directions
  goma_normal_copy(coord[(0 + permute_max) % 3], normal);
  goma_normal_copy(coord[(1 + permute_max) % 3], tangent);
  goma_normal_copy(coord[(2 + permute_max) % 3], binormal);

  goma_normal_free(seed);
  goma_normal_free(normal);
  goma_normal_free(tangent);
  goma_normal_free(binormal);
  return GOMA_SUCCESS;
}

goma_error
goma_edge_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]) {

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
  goma_normal *first = goma_normal_alloc(3);
  goma_normal *second = goma_normal_alloc(3);
  goma_normal *third = goma_normal_alloc(3);

  goma_normal_copy(first, normals[critical_angle[0]]);
  goma_normal_normalize(first);
  goma_normal_cross(first, normals[critical_angle[1]], second);
  goma_normal_normalize(second);
  goma_normal_cross(first, second, third);
  goma_normal_normalize(third);

  int permute_max = 0;
  double max_val = 0;
  // Find best coordinate direction
  for (int p = 0; p < 3; p++) {
    double val = gds_vector_get(first->normal, (0 + p) % 3) + gds_vector_get(second->normal, (1 + p) % 3) +
                 gds_vector_get(third->normal, (2 + p) % 3);
    if (val > max_val) {
      permute_max = p;
      max_val = val;
    }
  }

  // set best directions
  goma_normal_copy(coord[(0 + permute_max) % 3], first);
  goma_normal_copy(coord[(1 + permute_max) % 3], second);
  goma_normal_copy(coord[(2 + permute_max) % 3], third);

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
    bool match_one = fabs(gds_vector_dot(normals[u_index]->normal, normals[crit_normal[0]]->normal)) <
                         cos(GOMA_ROTATION_CRITICAL_ANGLE) ||
                     fabs(gds_vector_dot(normals[u_index]->normal, normals[crit_normal[1]]->normal)) <
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
  } else {
    return GOMA_ERROR;
  }

  return GOMA_SUCCESS;
}
