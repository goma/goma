#ifndef UTIL_GOMA_NORMAL_H
#define UTIL_GOMA_NORMAL_H

#include "gds/gds_vector.h"

typedef struct {
  gds_vector *normal;
  gds_vector *d_normal_dx[3][MDE];
} goma_normal;

typedef struct {
  double val;
  double d_val[3][MDE];
} goma_normal_val;

goma_normal *goma_normal_alloc(int size);

goma_normal *goma_get_average_normal(goma_normal **normals, int n_normals);

void goma_normal_zero(goma_normal *normal);

void goma_normal_add(goma_normal *y, goma_normal *x);

void goma_normal_sub(goma_normal *y, goma_normal *x);

void goma_normal_normalize(goma_normal *normal);

goma_normal_val goma_normal_get(goma_normal *normal, int index);

void goma_normal_set(goma_normal *normal, int index, const goma_normal_val *val);

void goma_normal_set_constant(goma_normal *normal, int index, double val);

void goma_normal_copy(goma_normal *dest, goma_normal *src);

void goma_normal_cross(goma_normal *u, goma_normal *v, goma_normal *cross);

void goma_normal_free(goma_normal *normal);

goma_normal_val goma_normal_dot(goma_normal *u, goma_normal *v);

void goma_normal_scale(goma_normal *normal, const goma_normal_val *val);

void goma_normal_rotate_around_vector(goma_normal *rotated,
                                      goma_normal *vec_to_rotate,
                                      goma_normal *axis,
                                      goma_normal_val angle_radians);

goma_normal_val cos_goma_normal_val(const goma_normal_val *val);

goma_normal_val sin_goma_normal_val(const goma_normal_val *val);

goma_normal_val fabs_goma_normal_val(const goma_normal_val *val);

goma_normal_val goma_normal_val_inverse(goma_normal_val *val);

goma_normal_val acos_goma_normal_val(const goma_normal_val *val);

goma_normal_val add_goma_normal_val(const goma_normal_val *val, double value);

goma_normal_val sub_goma_normal_val(const goma_normal_val *val, double value);

goma_normal_val scale_goma_normal_val(const goma_normal_val *val, double value);

goma_normal_val mul_goma_normal_val(const goma_normal_val *left, const goma_normal_val *right);

goma_normal_val sqrt_goma_normal_val(const goma_normal_val *val);

#endif // UTIL_GOMA_NORMAL_H
