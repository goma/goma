#ifndef GOMA_ROTATE_UTIL
#define GOMA_ROTATE_UTIL
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_bc_const.h"
#include "el_elm.h"

#include <stdbool.h>

#ifndef GOMA_MAX_NORMALS_PER_NODE
#define GOMA_MAX_NORMALS_PER_NODE 6
#endif

#ifndef GOMA_ROTATION_CRITICAL_ANGLE
#define GOMA_ROTATION_CRITICAL_ANGLE 45
#endif


typedef enum {
  GOMA_ROTATION_SURFACE, // all normals within GOMA_ROTATION_CRITICAL_ANGLE
  GOMA_ROTATION_EDGE,    // all normals are facing in 2 directions
  GOMA_ROTATION_CORNER,   // normals are facing in 3 directions
  GOMA_ROTATION_SINGLE   // single node is rotated, easy case
} goma_rotation_type_e;

typedef struct {
  double e[3];
} goma_rotation_vector_s;

typedef struct {
  goma_rotation_vector_s normals[GOMA_MAX_NORMALS_PER_NODE];
  goma_rotation_vector_s average_normals[GOMA_MAX_NORMALS_PER_NODE];
  goma_rotation_vector_s tangent1s[GOMA_MAX_NORMALS_PER_NODE];
  goma_rotation_vector_s tangent2s[GOMA_MAX_NORMALS_PER_NODE];
  goma_rotation_vector_s rotated_coord[DIM];
  int associate_direction[GOMA_MAX_NORMALS_PER_NODE];
  bool direction_is_associated[GOMA_MAX_NORMALS_PER_NODE];
  int tangent1_seeddir[GOMA_MAX_NORMALS_PER_NODE];
  int element[GOMA_MAX_NORMALS_PER_NODE];
  int face[GOMA_MAX_NORMALS_PER_NODE];
  int n_normals;
  goma_rotation_type_e type;
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
#endif
