#ifndef GOMA_BC_ROTATE_UTIL_H
#define GOMA_BC_ROTATE_UTIL_H

#include "mm_eh.h"
#include "util/goma_normal.h"

goma_normal *goma_get_average_normal(goma_normal **normals, int n_normals);

goma_error
goma_surface_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]);

goma_error goma_edge_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]);

goma_error goma_corner_coordinate_system(goma_normal **normals, int n_normals, goma_normal *coord[3]);

bool goma_check_normals_within_critical_angle(goma_normal **normals, int n_normals);

bool goma_check_corner_rotation_case(goma_normal **normals, int n_normals);

bool goma_check_edge_rotation_case(goma_normal **normals, int n_normals);

goma_error
goma_best_coordinate_system_3D(goma_normal **normals, int n_normals, goma_normal *coord[3]);

#endif // GOMA_BC_ROTATE_UTIL_H
