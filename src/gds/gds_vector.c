#include "gds/gds_vector.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const size_t gds_vector_size = sizeof(gds_vector);
static const size_t double_size = sizeof(double);

gds_vector *gds_vector_alloc(size_t size) {
  gds_vector *v = malloc(gds_vector_size + double_size * size);
  return v;
}

double gds_vector_get(const gds_vector *v, const size_t index) {
  assert(index < v->size);
  return v->data[index];
}

void gds_vector_set(gds_vector *v, const size_t index, const double value) {
  assert(index < v->size);
  v->data[index] = value;
}

void gds_vector_set_all(gds_vector *v, const double value) {
  for (size_t i = 0; i < v->size; i++) {
    v->data[i] = value;
  }
}

void gds_vector_zero(gds_vector *v) { gds_vector_set_all(v, 0.0); }

void gds_vector_ones(gds_vector *v) { gds_vector_set_all(v, 1.0); }

void gds_vector_copy(gds_vector *dest, const gds_vector *src) {
  assert(dest->data != src->data);
  assert(src->size == dest->size);
  memcpy(dest->data, src->data, double_size * src->size);
}

void gds_vector_add(gds_vector *y, const gds_vector *x) {
  assert(x->size == y->size);
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] += x->data[i];
  }
}

void gds_vector_sub(gds_vector *y, const gds_vector *x) {
  assert(x->size == y->size);
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] -= x->data[i];
  }
}

void gds_vector_mul(gds_vector *y, const gds_vector *x) {
  assert(x->size == y->size);
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] *= x->data[i];
  }
}

void gds_vector_div(gds_vector *y, const gds_vector *x) {
  assert(x->size == y->size);
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] /= x->data[i];
  }
}

void gds_vector_scale(gds_vector *y, const double b) {
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] *= b;
  }
}

void gds_vector_add_constant(gds_vector *y, const double b) {
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] += b;
  }
}

void gds_vector_axpy(const double alpha, const gds_vector *x, const double beta, gds_vector *y) {
  assert(x->size == y->size);
  for (size_t i = 0; i < y->size; i++) {
    y->data[i] = y->data[i] * beta + alpha * x->data[i];
  }
}

void gds_vector_normalize(gds_vector *v) {
  double dot = gds_vector_dot(v, v);
  if (dot == 0) {
    return;
  }
  double inv_mag = 1.0 / sqrt(dot);
  gds_vector_scale(v, inv_mag);
}

void gds_vector_cross(const gds_vector *v, const gds_vector *u, gds_vector *cross) {
  assert(u->size == 3);
  assert(v->size == 3);
  assert(cross->size == 3);
  gds_vector * cross_tmp = gds_vector_alloc(3);
  cross_tmp->data[0] = (v->data[1] * u->data[2] - v->data[2] * u->data[1]);
  cross_tmp->data[1] = (v->data[2] * u->data[0] - v->data[0] * u->data[2]);
  cross_tmp->data[2] = (v->data[0] * u->data[1] - v->data[1] * u->data[0]);
  gds_vector_copy(cross, cross_tmp);
  gds_vector_free(cross_tmp);
}

double gds_vector_dot(const gds_vector *v, const gds_vector *u) {
  assert(v->size == u->size);
  double dot = 0;
  for (size_t i = 0; i < u->size; i++) {
    dot += v->data[i] * u->data[i];
  }
  return dot;
}

void gds_vector_rotate_around_vector(gds_vector *rotated, const gds_vector *vec_to_rotate,
                                     const gds_vector *axis, double angle_radians) {
  assert(rotated->size == 3);
  assert(rotated->size == vec_to_rotate->size);
  assert(rotated->size == axis->size);

  size_t size = axis->size;
  gds_vector *cross = gds_vector_alloc(size);
  gds_vector *temp1 = gds_vector_alloc(size);
  gds_vector *temp2 = gds_vector_alloc(size);
  gds_vector_cross(axis, vec_to_rotate, cross);

  gds_vector_copy(temp1, vec_to_rotate);
  gds_vector_scale(temp1, cos(angle_radians));

  gds_vector_scale(cross, sin(angle_radians));
  double dot = gds_vector_dot(vec_to_rotate, axis);

  gds_vector_copy(temp2, axis);
  gds_vector_scale(temp2, dot * (1 - cos(angle_radians)));

  gds_vector_zero(rotated);
  gds_vector_add(rotated, temp1);
  gds_vector_add(rotated, cross);
  gds_vector_add(rotated, temp2);

  gds_vector_free(cross);
  gds_vector_free(temp1);
  gds_vector_free(temp2);
}

void gds_vector_free(gds_vector *v) {
  free(v);
}
