#ifndef GOMA_BC_ROTATE_UTIL_H
#define GOMA_BC_ROTATE_UTIL_H

#include "gds/gds_vector.h"
#include "mm_eh.h"

gds_vector *goma_get_average_normal(gds_vector **normals, int n_normals);

goma_error
goma_surface_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]);

goma_error goma_edge_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]);

goma_error goma_corner_coordinate_system(gds_vector **normals, int n_normals, gds_vector *coord[3]);

bool goma_check_normals_within_critical_angle(gds_vector **normals, int n_normals);

bool goma_check_corner_rotation_case(gds_vector **normals, int n_normals);

bool goma_check_edge_rotation_case(gds_vector **normals, int n_normals);

goma_error
goma_best_coordinate_system_3D(gds_vector **normals, int n_normals, gds_vector *coord[3]);

#endif // GOMA_BC_ROTATE_UTIL_H
