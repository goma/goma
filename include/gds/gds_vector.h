#ifndef GDS_VECTOR_H
#define GDS_VECTOR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double *data;
  size_t size;
} gds_vector;

gds_vector *gds_vector_alloc(size_t size);

double gds_vector_get(const gds_vector *v, const size_t index);

void gds_vector_set(gds_vector *v, const size_t index, const double value);

void gds_vector_set_all(gds_vector *v, const double value);

void gds_vector_zero(gds_vector *v);

void gds_vector_ones(gds_vector *v);

void gds_vector_copy(gds_vector *dest, const gds_vector *src);

void gds_vector_add(gds_vector *y, const gds_vector *x);

void gds_vector_sub(gds_vector *y, const gds_vector *x);

void gds_vector_mul(gds_vector *y, const gds_vector *x);

void gds_vector_div(gds_vector *y, const gds_vector *x);

void gds_vector_scale(gds_vector *y, const double b);

void gds_vector_add_constant(gds_vector *y, const double b);

void gds_vector_axpy(const double alpha, const gds_vector *x, const double beta, gds_vector *y);

void gds_vector_normalize(gds_vector *v);

void gds_vector_cross(const gds_vector *v, const gds_vector *u, gds_vector *cross);

void gds_vector_rotate_around_vector(gds_vector *rotated, const gds_vector *vec_to_rotate,
                                     const gds_vector *axis, double angle_radians);

double gds_vector_dot(const gds_vector *v, const gds_vector *u);

void gds_vector_free(gds_vector *v);

#ifdef __cplusplus
}
#endif

#endif
