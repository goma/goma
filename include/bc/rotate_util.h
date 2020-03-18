#ifndef GOMA_BC_ROTATE_UTIL_H
#define GOMA_BC_ROTATE_UTIL_H

#include "mm_eh.h"
#include "gds/gds_vector.h"

goma_error goma_rotated_coordinates_from_normals(gds_vector **normals,
                                                 unsigned int n_normals,
                                                 gds_vector *coord[3],
                                                 unsigned int normal_dirs[3]);
#endif // GOMA_BC_ROTATE_UTIL_H
