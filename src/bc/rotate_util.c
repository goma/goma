#include "bc/rotate_util.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bc/rotate_coordinates.h"
#include "gds/gds_vector.h"
#include "mm_eh.h"
#include "std.h"
#include "util/goma_normal.h"

static bool within_critical_angle(double a_dot_b, double angle) {
  return a_dot_b <= (1 + 1e-14) && a_dot_b >= cos(angle);
}

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

  if ((val = (fabs(gds_vector_get(normal->normal, 2)) + fabs(gds_vector_get(tangent->normal, 1)) +
              fabs(gds_vector_get(binormal->normal, 0)))) > max_val) {
    max_val = val;
    n_max = 2;
    t_max = 1;
    b_max = 0;
  }

  if ((val = (fabs(gds_vector_get(normal->normal, 2)) + fabs(gds_vector_get(tangent->normal, 0)) +
              fabs(gds_vector_get(binormal->normal, 1)))) > max_val) {
    max_val = val;
    n_max = 2;
    t_max = 0;
    b_max = 1;
  }

  if ((val = (fabs(gds_vector_get(normal->normal, 0)) + fabs(gds_vector_get(tangent->normal, 2)) +
              fabs(gds_vector_get(binormal->normal, 1)))) > max_val) {
    max_val = val;
    n_max = 0;
    t_max = 2;
    b_max = 1;
  }

  if ((val = (fabs(gds_vector_get(normal->normal, 1)) + fabs(gds_vector_get(tangent->normal, 0)) +
              fabs(gds_vector_get(binormal->normal, 2)))) > max_val) {
    max_val = val;
    n_max = 1;
    t_max = 0;
    b_max = 2;
  }

  if ((val = (fabs(gds_vector_get(normal->normal, 1)) + fabs(gds_vector_get(tangent->normal, 2)) +
              fabs(gds_vector_get(binormal->normal, 0)))) > max_val) {
    max_val = val;
    n_max = 1;
    t_max = 2;
    b_max = 0;
  }
  // Find best coordinate direction
  // for (int p = 0; p < 3; p++) {
  //  double val = fabs(gds_vector_get(normal->normal, (0 + p) % 3)) +
  //               fabs(gds_vector_get(tangent->normal, (1 + p) % 3)) +
  //               fabs(gds_vector_get(binormal->normal, (2 + p) % 3));
  //  if (val > max_val) {
  //    permute_max = p;
  //    max_val = val;
  //  }
  //}
  if (n_max == -1 || t_max == -1 || b_max == -1) {
    GOMA_EH(GOMA_ERROR, "could not associate normals to coord [%g,%g,%g], [%g,%g,%g], [%g,%g,%g]",
            normal->normal->data[0], normal->normal->data[1], normal->normal->data[2],
            tangent->normal->data[0], tangent->normal->data[1], tangent->normal->data[2],
            binormal->normal->data[0], binormal->normal->data[1], binormal->normal->data[2]);
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
    goma_normal_val dot = goma_normal_get(normal, i);
    if (fabs(dot.val) < dot_min) {
      worst = i;
      dot_min = fabs(dot.val);
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
      double uv_dot = (gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
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

    if (within_critical_angle(gds_vector_dot(normals[i]->normal, n1[0]->normal),
                              GOMA_ROTATION_CRITICAL_ANGLE)) {
      n1[n_n1++] = normals[i];
    } else if (within_critical_angle(gds_vector_dot(normals[i]->normal, n2[0]->normal),
                                     GOMA_ROTATION_CRITICAL_ANGLE)) {
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
      double uv_dot = (gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = (gds_vector_dot(normals[u_index]->normal, normals[w_index]->normal));
          double vw_dot = (gds_vector_dot(normals[v_index]->normal, normals[w_index]->normal));
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

corner_coord_critical_found: {
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
    if (dot > max_dot && i != first_crit) {
      second_crit = i;
      max_dot = dot;
    }
  }

  if (second_crit == -1) {
    return GOMA_ERROR;
  }

  int third_crit = -1;
  for (int i = 0; i < 3; i++) {
    if (i != first_crit && i != second_crit) {
      third_crit = i;
    }
  }
  if (third_crit == -1) {
    return GOMA_ERROR;
  }

  goma_normal *first = goma_normal_alloc(3);
  goma_normal *second = goma_normal_alloc(3);
  goma_normal *third = goma_normal_alloc(3);

  goma_normal_copy(first, normals[critical_angle[first_crit]]);
  goma_normal_normalize(first);
  goma_normal_copy(second, normals[critical_angle[second_crit]]);
  goma_normal_normalize(second);
  goma_normal_copy(third, normals[critical_angle[third_crit]]);
  goma_normal_normalize(third);

  goma_normal_val dot = goma_normal_dot(first, second);
  const double ninety = 90 * M_PI / 180;

  goma_normal *new_n1 = goma_normal_alloc(3);
  goma_normal *new_n2 = goma_normal_alloc(3);
  goma_normal *cross = goma_normal_alloc(3);
  goma_normal *mid_n1_n2 = goma_normal_alloc(3);
  goma_normal_zero(cross);
  goma_normal_val angle = acos_goma_normal_val(&dot);
  angle = scale_goma_normal_val(&angle, 0);
  angle = add_goma_normal_val(&angle, ninety);
  const double half = 0.5;
  goma_normal_val shift = scale_goma_normal_val(&angle, half);
  goma_normal_val minus_shift = scale_goma_normal_val(&angle, -half);
  goma_normal_cross(first, second, cross);

  goma_normal_copy(mid_n1_n2, first);
  goma_normal_add(mid_n1_n2, second);
  goma_normal_normalize(mid_n1_n2);
  goma_normal_normalize(cross);
  goma_normal_rotate_around_vector(new_n1, mid_n1_n2, cross, shift);
  goma_normal_rotate_around_vector(new_n2, mid_n1_n2, cross, minus_shift);

#ifndef NDEBUG
  goma_normal_val dt1 = goma_normal_dot(new_n1, new_n2);
  GOMA_ASSERT(fabs(dt1.val) < 1e-14);
#endif

  goma_normal_normalize(new_n1);
  goma_normal_normalize(new_n2);

  goma_normal *final_n1 = goma_normal_alloc(3);
  goma_normal *final_n2 = goma_normal_alloc(3);
  goma_normal *final_n3 = goma_normal_alloc(3);
  goma_normal_copy(mid_n1_n2, new_n1);
  goma_normal_add(mid_n1_n2, new_n2);
  goma_normal_normalize(mid_n1_n2);
  goma_normal_zero(cross);

  goma_normal_cross(new_n1, new_n2, cross);
  goma_normal_normalize(cross);
  dot = goma_normal_dot(cross, third);
  if (dot.val < 0) {
    goma_normal_val val = {0};
    val.val = -1;
    goma_normal_scale(third, &val);
  }

  dot = goma_normal_dot(cross, third);
  // use mid_n1_n2 as cross X third
  goma_normal_cross(cross, third, mid_n1_n2);
  goma_normal_normalize(mid_n1_n2);
  angle = acos_goma_normal_val(&dot);
  shift = scale_goma_normal_val(&angle, half);
  goma_normal_rotate_around_vector(final_n1, new_n1, mid_n1_n2, shift);
  goma_normal_rotate_around_vector(final_n2, new_n2, mid_n1_n2, shift);
  goma_normal_rotate_around_vector(final_n3, cross, mid_n1_n2, shift);
  goma_normal_normalize(final_n1);
  goma_normal_normalize(final_n2);
  goma_normal_normalize(final_n3);

#ifndef NDEBUG
  goma_normal_val dt2 = goma_normal_dot(final_n1, final_n2);
  GOMA_ASSERT(fabs(dt2.val) < 1e-14);
  goma_normal_val dt3 = goma_normal_dot(final_n1, final_n3);
  GOMA_ASSERT(fabs(dt3.val) < 1e-14);
  goma_normal_val dt4 = goma_normal_dot(final_n2, final_n3);
  GOMA_ASSERT(fabs(dt4.val) < 1e-14);
#endif
  goma_normal_assign_best_direction(final_n1, final_n2, final_n3, coord);

  goma_normal_free(first);
  goma_normal_free(second);
  goma_normal_free(third);
  goma_normal_free(cross);
  goma_normal_free(new_n1);
  goma_normal_free(new_n2);
  goma_normal_free(mid_n1_n2);
  goma_normal_free(final_n1);
  goma_normal_free(final_n2);
  goma_normal_free(final_n3);
}

  return GOMA_SUCCESS;
}

bool goma_check_normals_within_critical_angle(goma_normal **normals, int n_normals) {
  double min_dot = DBL_MAX;
  for (int i = 0; i < n_normals; i++) {
    for (int j = 0; j < n_normals; j++) {
      double dot = (gds_vector_dot(normals[i]->normal, normals[j]->normal));
      if (dot < min_dot) {
        min_dot = dot;
      }
    }
  }

  // True if all are within a critical angle
  return within_critical_angle(min_dot, GOMA_ROTATION_CRITICAL_ANGLE);
}

bool goma_check_corner_rotation_case(goma_normal **normals, int n_normals) {
  if (n_normals < 3) {
    return false;
  }
  // If two normals are larger than critical angle then there should be a third normal larger
  // than a critical angle to both those to be a corner type
  for (int u_index = 0; u_index < n_normals; u_index++) {
    for (int v_index = u_index + 1; v_index < n_normals; v_index++) {
      double uv_dot = (gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
      if (uv_dot < cos(GOMA_ROTATION_CRITICAL_ANGLE)) {
        for (int w_index = 0; w_index < n_normals; w_index++) {
          if (w_index == u_index || w_index == v_index) {
            continue;
          }
          double uw_dot = (gds_vector_dot(normals[u_index]->normal, normals[w_index]->normal));
          double vw_dot = (gds_vector_dot(normals[v_index]->normal, normals[w_index]->normal));
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
      double uv_dot = (gds_vector_dot(normals[u_index]->normal, normals[v_index]->normal));
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
    bool match_one = within_critical_angle(
                         gds_vector_dot(normals[u_index]->normal, normals[crit_normal[0]]->normal),
                         GOMA_ROTATION_CRITICAL_ANGLE) ||
                     within_critical_angle(
                         gds_vector_dot(normals[u_index]->normal, normals[crit_normal[1]]->normal),
                         GOMA_ROTATION_CRITICAL_ANGLE);
    if (!match_one) {
      return false;
    }
  }

  return true;
}

goma_error goma_best_coordinate_system_3D(goma_normal **normals,
                                          int n_normals,
                                          goma_normal *coord[3],
                                          int *type) {

  // normalize just to be sure
  for (int i = 0; i < n_normals; i++) {
    goma_normal_normalize(normals[i]);
  }
  int err = GOMA_SUCCESS;
  // check if all are within a critical angle
  // then we have a surface case
  if (goma_check_normals_within_critical_angle(normals, n_normals)) {
    *type = 0;
    err = goma_surface_coordinate_system(normals, n_normals, coord);
    if (err == GOMA_ERROR) {
      fprintf(stderr, "best rot surface failure");
    }
  } else if (goma_check_corner_rotation_case(normals, n_normals)) {
    *type = 1;
    err = goma_corner_coordinate_system(normals, n_normals, coord);
    if (err == GOMA_ERROR) {
      fprintf(stderr, "best rot corner failure");
    }
  } else if (goma_check_edge_rotation_case(normals, n_normals)) {
    *type = 2;
    err = goma_edge_coordinate_system(normals, n_normals, coord);
    if (err == GOMA_ERROR) {
      fprintf(stderr, "best rot edge failure");
    }
  }
  return err;
}
