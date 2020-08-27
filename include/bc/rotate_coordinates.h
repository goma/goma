#ifndef GOMA_ROTATE_UTIL
#define GOMA_ROTATE_UTIL
#include <stdbool.h>

#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_bc_const.h"
#include "el_elm.h"
#include "gds/gds_vector.h"

struct Boundary_Condition;

#ifndef GOMA_ROTATION_CRITICAL_ANGLE
#define GOMA_ROTATION_CRITICAL_ANGLE 45
#endif

typedef struct {
  gds_vector *rotated_coord[DIM];
  bool is_rotated;
} goma_rotation_node_s;

typedef struct {
  bool automatic_rotations;
  goma_rotation_node_s *rotation_nodes;
} goma_rotation_s;

extern goma_rotation_s goma_automatic_rotations;

goma_error check_if_equation_is_rotation(int equation, bool *is_rotated);

goma_error setup_bc_is_rotated_list(struct Boundary_Condition *bc_types, int num_bc,
                                    bool **bc_rotate_list);

goma_error setup_rotated_bc_nodes(Exo_DB *exo, struct Boundary_Condition *bc_types, int num_bc,
                                  double *x);

goma_error set_rotation_types(Exo_DB *exo, goma_rotation_node_s *rotation);

goma_error associate_directions(Exo_DB *exo, goma_rotation_node_s *rotation);
goma_error set_average_normals_and_tangents(Exo_DB *exo, goma_rotation_node_s *rotation);
goma_error set_rotated_coordinate_system(Exo_DB *exo, goma_rotation_node_s *rotation);
goma_error free_rotations(Exo_DB *exo, goma_rotation_node_s **rotations);
goma_error set_face_normal_association(Exo_DB *exo, goma_rotation_node_s *rotations);
int offset_from_rotated_equation(int eqn);
#endif
