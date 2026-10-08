#ifndef GDS_VECTOR_H
#define GDS_VECTOR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t size;
  size_t struct_size;
  // flexible array size (uses [1] so valid in c++, c99 would allow [])
  double data[1];
} gds_vector;

void gds_print_vector(const gds_vector *v);

/**
 * @brief gds_vector_alloc
 * @param size
 * @return  allocated vector with size elements
 */
gds_vector *gds_vector_alloc(size_t size);

/**
 * @brief gds_vector_get
 * @param v
 * @param index
 * @return  v[index]
 */
double gds_vector_get(const gds_vector *v, const size_t index);

/**
 * @brief gds_vector_set v[index] = value
 * @param v
 * @param index
 * @param value
 */
void gds_vector_set(gds_vector *v, const size_t index, const double value);

/**
 * @brief gds_vector_set_all set all values of v to value
 * @param v
 * @param value
 */
void gds_vector_set_all(gds_vector *v, const double value);

/**
 * @brief gds_vector_zero set all values of v to 0.0
 * @param v
 */
void gds_vector_zero(gds_vector *v);

/**
 * @brief gds_vector_ones set all values of v to 1.0
 * @param v
 */
void gds_vector_ones(gds_vector *v);

/**
 * @brief gds_vector_copy copy src to dest, must be of same size
 * @param dest dest[i] = src[i]
 * @param src
 */
void gds_vector_copy(gds_vector *dest, const gds_vector *src);

/**
 * @brief gds_vector_add
 * @param y y[i] += x[i]
 * @param x
 */
void gds_vector_add(gds_vector *y, const gds_vector *x);

/**
 * @brief gds_vector_sub
 * @param y y[i] -= x[i]
 * @param x
 */
void gds_vector_sub(gds_vector *y, const gds_vector *x);

/**
 * @brief gds_vector_mul
 * @param y y[i] *= x[i]
 * @param x
 */
void gds_vector_mul(gds_vector *y, const gds_vector *x);

/**
 * @brief gds_vector_div
 * @param y y[i] /= x[i]
 * @param x
 */
void gds_vector_div(gds_vector *y, const gds_vector *x);

/**
 * @brief gds_vector_scale
 * @param y y = y*b
 * @param b
 */
void gds_vector_scale(gds_vector *y, const double b);

/**
 * @brief gds_vector_add_constant
 * @param y y = y+b all elements of y have b added
 * @param b
 */
void gds_vector_add_constant(gds_vector *y, const double b);

/**
 * @brief gds_vector_axpy
 * @param alpha
 * @param x
 * @param beta
 * @param y y = alpha * x + beta * y
 */
void gds_vector_axpy(const double alpha, const gds_vector *x, const double beta, gds_vector *y);

/**
 * @brief gds_vector_normalize
 * @param v in/out normalizes v as v = v / mag(v)
 */
void gds_vector_normalize(gds_vector *v);

/**
 * @brief gds_vector_cross
 * @param v
 * @param u
 * @param cross return cross product of v and u
 */
void gds_vector_cross(const gds_vector *v, const gds_vector *u, gds_vector *cross);

/**
 * @brief gds_vector_rotate_around_vector
 * @param rotated rotated vector
 * @param vec_to_rotate vector to be rotated
 * @param axis axis to rotate around
 * @param angle_radians angle to rotate
 */
void gds_vector_rotate_around_vector(gds_vector *rotated,
                                     const gds_vector *vec_to_rotate,
                                     const gds_vector *axis,
                                     double angle_radians);

/**
 * @brief gds_vector_dot
 * @param v
 * @param u
 * @return u dot v
 */
double gds_vector_dot(const gds_vector *v, const gds_vector *u);

/**
 * @brief gds_vector_free free allocated vector v
 * @param v
 */
void gds_vector_free(gds_vector *v);

#ifdef __cplusplus
}
#endif

#endif
