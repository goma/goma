#include "mm_eh.h"
#include "std.h"
#include <assert.h>
#include "util/goma_normal.h"

goma_normal *goma_normal_alloc(int size) {
  goma_normal *normal = malloc(sizeof(goma_normal));
  normal->normal = gds_vector_alloc(3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      normal->d_normal_dx[i][j] = gds_vector_alloc(3);
    }
  }
  return normal;
}

goma_normal * goma_get_average_normal(goma_normal **normals, int n_normals) {

  goma_normal *average = goma_normal_alloc(3);
  goma_normal_zero(average);

  for (int i = 0; i < n_normals; i++) {
    goma_normal_add(average, normals[i]);
  }

  goma_normal_normalize(average);

  return average;

}


void goma_normal_zero(goma_normal *normal) {

  gds_vector_zero(normal->normal);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_zero(normal->d_normal_dx[i][j]);
    }
  }
}

void goma_normal_add(goma_normal *y, goma_normal *x) {
  gds_vector_add(y->normal, x->normal);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_add(y->d_normal_dx[i][j], x->d_normal_dx[i][j]);
    }
  }
}

void goma_normal_sub(goma_normal *y, goma_normal *x) {
  gds_vector_sub(y->normal, x->normal);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_sub(y->d_normal_dx[i][j], x->d_normal_dx[i][j]);
    }
  }
}

goma_normal_val goma_normal_val_inverse(goma_normal_val *val) {

  goma_normal_val inv;
  inv.val = 1.0 / val->val;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      val->d_val[i][j] = - val->d_val[i][j] * inv.val * inv.val;
    }
  }

  return inv;
}

void goma_normal_normalize(goma_normal *normal) {
  goma_normal_val dot = goma_normal_dot(normal, normal);
  if (dot.val == 0) {
    return;
  }

  goma_normal_val inv = goma_normal_val_inverse(&dot);

  goma_normal_scale(normal, &inv);
}

goma_normal_val goma_normal_get(goma_normal *normal, int index) {
  goma_normal_val val;

  val.val = gds_vector_get(normal->normal, index);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      val.d_val[i][j] = gds_vector_get(normal->d_normal_dx[i][j], index);
    }
  }

  return val;
}

void goma_normal_set(goma_normal *normal, int index, const goma_normal_val *val) {

  gds_vector_set(normal->normal, index, val->val);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_set(normal->d_normal_dx[i][j], index, val->d_val[i][j]);
    }
  }
}

void goma_normal_set_constant(goma_normal *normal, int index, double val) {
  gds_vector_set(normal->normal, index, val);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_set(normal->d_normal_dx[i][j], index, 0.0);
    }
  }
}

void goma_normal_copy(goma_normal *dest, goma_normal *src) {
  gds_vector_copy(dest->normal, src->normal);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_copy(dest->d_normal_dx[i][j], src->d_normal_dx[i][j]);
    }
  }
}

void goma_normal_cross(goma_normal *u, goma_normal *v, goma_normal *cross) {
  gds_vector_cross(u->normal, v->normal, cross->normal);

  // d/dx (a X b) = d/dx a X b + a X d/dx b
  gds_vector *tmp = gds_vector_alloc(3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_cross(u->d_normal_dx[i][j], v->normal, cross->d_normal_dx[i][j]);
      gds_vector_cross(u->normal, v->d_normal_dx[i][j], tmp);
      gds_vector_add(cross->d_normal_dx[i][j], tmp);
    }
  }

  gds_vector_free(tmp);
}

void goma_normal_free(goma_normal *normal) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      free(normal->d_normal_dx[i][j]);
    }
  }
  free(normal);
}

goma_normal_val goma_normal_dot(goma_normal *u, goma_normal *v) {
  goma_normal_val val;
  val.val = gds_vector_dot(u->normal, v->normal);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      val.d_val[i][j] = gds_vector_dot(u->normal, v->d_normal_dx[i][j]) +
                        gds_vector_dot(u->d_normal_dx[i][j], v->normal);
    }
  }

  return val;
}

void goma_normal_scale(goma_normal *normal, const goma_normal_val *val) {
  gds_vector *tmp1 = gds_vector_alloc(3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      gds_vector_copy(tmp1, normal->normal);
      gds_vector_scale(normal->d_normal_dx[i][j], val->val);
      gds_vector_scale(tmp1, val->d_val[i][j]);
      gds_vector_add(normal->d_normal_dx[i][j], tmp1);
    }
  }
  gds_vector_scale(normal->normal, val->val);
  gds_vector_free(tmp1);
}

goma_normal_val fabs_goma_normal_val(const goma_normal_val *val) {

  goma_normal_val retval;

  retval.val = fabs(val->val);

  double sign = SGN(val->val);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < MDE; j++) {
      retval.d_val[i][j] = sign * val->d_val[i][j];
    }
  }
  return retval;
}

