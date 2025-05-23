

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <stddef.h>

#include "bc/rotate_util.h"
#include "gds/gds_vector.h"
#include "mm_eh.h"
#include "util/goma_normal.h"

#define INV_SQRT_2 M_SQRT1_2

static const auto zero_comp = Catch::Matchers::WithinAbsMatcher(0.0, 1e-15);

// TODO Unit Test Derivatives, need a good way to do this

static void helper_create_normals(goma_normal ***normals, int n_normals) {
  *normals = static_cast<goma_normal **>(malloc(sizeof(goma_normal *) * n_normals));
  for (int i = 0; i < n_normals; i++) {
    (*normals)[i] = goma_normal_alloc(3);
    goma_normal_zero((*normals)[i]);
  }
}

static void helper_free_normals(goma_normal ***normals, int n_normals) {
  for (int i = 0; i < n_normals; i++) {
    goma_normal_free((*normals)[i]);
  }
  free(*normals);
}

static void helper_set_normal_from_vector(goma_normal *normal, double vec[3]) {
  goma_normal_set_constant(normal, 0, vec[0]);
  goma_normal_set_constant(normal, 1, vec[1]);
  goma_normal_set_constant(normal, 2, vec[2]);
  goma_normal_normalize(normal);
}

static bool helper_normals_match(goma_normal *left, goma_normal *right) {
  for (int j = 0; j < 3; j++) {
    double cval = gds_vector_get(left->normal, j);
    double nval = gds_vector_get(right->normal, j);
    if (!zero_comp.match(fabs(cval - nval))) {
      return false;
    }
  }

  return true;
}

static bool helper_is_unit_normal(goma_normal *normal) {
  goma_normal *tmp = goma_normal_alloc(3);

  // check all normals are unit
  goma_normal_copy(tmp, normal);
  goma_normal_normalize(tmp);

  bool is_unit = true;
  if (!helper_normals_match(tmp, normal)) {
    is_unit = false;
  }

  goma_normal_free(tmp);
  return is_unit;
}

static bool helper_one_coordinate_matches_normal(goma_normal *coord[3], goma_normal *normal) {
  for (int i = 0; i < 3; i++) {
    if (helper_normals_match(coord[i], normal)) {
      return true;
    }
  }

  return false;
}

static bool helper_one_coordinate_near_normal(goma_normal *coord[3], goma_normal *normal) {
  for (int i = 0; i < 3; i++) {
    if (fabs(gds_vector_dot(normal->normal, coord[i]->normal)) <
        cos(GOMA_ROTATIONS_CRITICAL_ANGLE)) {
      return true;
    }
  }

  return false;
}

static bool helper_coordinate_is_coordinate_system(goma_normal *coord[3]) {
  // check all normals are unit
  for (int i = 0; i < 3; i++) {
    if (!helper_is_unit_normal(coord[i])) {
      return false;
    }
  }

  // coord[0] cross coord[1] should == coord[2]
  goma_normal *tmp = goma_normal_alloc(3);
  goma_normal_cross(coord[0], coord[1], tmp);
  bool is_coord_system = helper_normals_match(coord[2], tmp);

  goma_normal_free(tmp);
  return is_coord_system;
}

TEST_CASE("goma_check_normals_within_critical_angle", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 3;
  helper_create_normals(&normals, n_normals);

  REQUIRE(GOMA_ROTATIONS_CRITICAL_ANGLE == Catch::Approx(M_PI_4));

  double n1[3] = {1.0, 0.0, 0.0};
  double n2[3] = {0.0, 1.0, 0.0};
  double n3[3] = {0.0, 0.0, 1.0};

  // normals are original coordinate system
  helper_set_normal_from_vector(normals[0], n1);
  helper_set_normal_from_vector(normals[1], n2);
  helper_set_normal_from_vector(normals[2], n3);

  REQUIRE_FALSE(goma_check_normals_within_critical_angle(normals, n_normals));

  // all normals are the same
  helper_set_normal_from_vector(normals[0], n2);
  helper_set_normal_from_vector(normals[2], n2);
  REQUIRE(goma_check_normals_within_critical_angle(normals, n_normals));

  helper_set_normal_from_vector(normals[0], n1);
  helper_set_normal_from_vector(normals[1], n1);
  helper_set_normal_from_vector(normals[2], n1);
  REQUIRE(goma_check_normals_within_critical_angle(normals, n_normals));

  helper_set_normal_from_vector(normals[0], n3);
  helper_set_normal_from_vector(normals[1], n3);
  helper_set_normal_from_vector(normals[2], n3);
  REQUIRE(goma_check_normals_within_critical_angle(normals, n_normals));

  double sq2n1[3] = {INV_SQRT_2, INV_SQRT_2, 0.0};
  double sq2n2[3] = {INV_SQRT_2, 0.0, INV_SQRT_2};
  double sq2n3[3] = {0.0, INV_SQRT_2, INV_SQRT_2};

  // not within 45 degrees
  helper_set_normal_from_vector(normals[0], sq2n1);
  helper_set_normal_from_vector(normals[1], sq2n2);
  helper_set_normal_from_vector(normals[2], sq2n3);
  REQUIRE_FALSE(goma_check_normals_within_critical_angle(normals, n_normals));

  // not all within 45 degrees
  helper_set_normal_from_vector(normals[0], sq2n1);
  helper_set_normal_from_vector(normals[1], sq2n2);
  helper_set_normal_from_vector(normals[2], sq2n1);
  REQUIRE_FALSE(goma_check_normals_within_critical_angle(normals, n_normals));

  helper_set_normal_from_vector(normals[0], sq2n2);
  helper_set_normal_from_vector(normals[1], sq2n2);
  helper_set_normal_from_vector(normals[2], sq2n3);
  REQUIRE_FALSE(goma_check_normals_within_critical_angle(normals, n_normals));

  // bias normals slightly towards x, helper will renormalize
  sq2n1[0] = INV_SQRT_2 + 1e-1;
  sq2n2[0] = INV_SQRT_2 + 1e-1;
  helper_set_normal_from_vector(normals[0], sq2n1);
  helper_set_normal_from_vector(normals[1], sq2n2);
  helper_set_normal_from_vector(normals[2], sq2n3);
  REQUIRE_FALSE(goma_check_normals_within_critical_angle(normals, n_normals));

  // check that when they are near x direction they are a surface
  // n1 == n3 n2 is slightly off
  helper_set_normal_from_vector(normals[1], sq2n1);
  sq2n1[0] += 0.1;
  helper_set_normal_from_vector(normals[2], sq2n1);
  REQUIRE(goma_check_normals_within_critical_angle(normals, n_normals));

  // check that when they are near x direction they are a surface
  // n1 == n3 n2 is [1 0 0]
  sq2n3[0] = 1.0;
  sq2n3[1] = 0.0;
  sq2n3[2] = 0.0;
  helper_set_normal_from_vector(normals[2], sq2n3);
  REQUIRE(goma_check_normals_within_critical_angle(normals, n_normals));

  // A few checks with only 1 normal, should always pass
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], n1);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], n2);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], n3);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], sq2n1);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], sq2n2);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_set_normal_from_vector(normals[0], sq2n3);
  REQUIRE(goma_check_normals_within_critical_angle(normals, 1));

  helper_free_normals(&normals, n_normals);
}

TEST_CASE("goma_check_edge_rotation_case", "[bc][automatic_rotations]") {}

TEST_CASE("goma_check_corner_rotation_case", "[bc][automatic_rotations]") {}

TEST_CASE("goma_best_coordinate_system_3D 1 normal", "[bc][automatic_rotations]") {
  goma_error error;
  goma_normal **normals;
  int n_normals = 1;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  double n1[3] = {1.0, 0.0, 0.0};
  double n2[3] = {0.0, 1.0, 0.0};
  double n3[3] = {0.0, 0.0, 1.0};

  // Check that we keep our default coordinate system
  helper_set_normal_from_vector(normals[0], n1);
  int type;
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(gds_vector_get(rotated_coord[0]->normal, 0) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 1)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 0)));
  REQUIRE(gds_vector_get(rotated_coord[1]->normal, 1) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 0)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 1)));
  REQUIRE(gds_vector_get(rotated_coord[2]->normal, 2) == Catch::Approx(1.0));

  helper_set_normal_from_vector(normals[0], n2);
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(gds_vector_get(rotated_coord[0]->normal, 0) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 1)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 0)));
  REQUIRE(gds_vector_get(rotated_coord[1]->normal, 1) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 0)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 1)));
  REQUIRE(gds_vector_get(rotated_coord[2]->normal, 2) == Catch::Approx(1.0));

  helper_set_normal_from_vector(normals[0], n3);
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(gds_vector_get(rotated_coord[0]->normal, 0) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 1)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[0]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 0)));
  REQUIRE(gds_vector_get(rotated_coord[1]->normal, 1) == Catch::Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[1]->normal, 2)));

  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 0)));
  REQUIRE(zero_comp.match(gds_vector_get(rotated_coord[2]->normal, 1)));
  REQUIRE(gds_vector_get(rotated_coord[2]->normal, 2) == Catch::Approx(1.0));

  // Try some different normals now
  double mn1[3] = {-1.0, 0.0, 0.0};
  helper_set_normal_from_vector(normals[0], mn1);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double mn2[3] = {0.0, -1.0, 0.0};
  helper_set_normal_from_vector(normals[0], mn2);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double mn3[3] = {0.0, 0.0, -1.0};
  helper_set_normal_from_vector(normals[0], mn3);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff1[3] = {INV_SQRT_2, INV_SQRT_2, 0};
  helper_set_normal_from_vector(normals[0], noff1);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff2[3] = {INV_SQRT_2, 0, INV_SQRT_2};
  helper_set_normal_from_vector(normals[0], noff2);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff3[3] = {0.4242640687119285, 0.565685424949238, 0.7071067811865475};
  helper_set_normal_from_vector(normals[0], noff3);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff4[3] = {-0.565685424949238, 0.4242640687119285, 0.7071067811865475};
  helper_set_normal_from_vector(normals[0], noff4);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff5[3] = {-0.565685424949238, 0.4242640687119285, -0.7071067811865475};
  helper_set_normal_from_vector(normals[0], noff5);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double noff6[3] = {0.565685424949238, -0.4242640687119285, -0.7071067811865475};
  helper_set_normal_from_vector(normals[0], noff6);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));
  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}

TEST_CASE("goma_best_coordinate_system_3D 2 normals surface", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 2;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  // SURFACE
  double m1[3] = {1.0, 0.0, 0.0};
  helper_set_normal_from_vector(normals[0], m1);

  double m2[3] = {1.0, 0.0, 0.0};
  helper_set_normal_from_vector(normals[1], m2);

  goma_error error;
  int type;
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  double m3[3] = {0.98058068, 0.19611614, 0.0};
  helper_set_normal_from_vector(normals[1], m3);
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));

  double m4[3] = {0.8660254, -0.5, 0};
  helper_set_normal_from_vector(normals[1], m4);
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));

  double m5[3] = {0.8660254, 0.5, 0};
  helper_set_normal_from_vector(normals[0], m5);
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));

  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}

TEST_CASE("goma_best_coordinate_system_3D 2 normals edge", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 2;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  goma_normal *other_normal = goma_normal_alloc(3);
  double zn[3] = {0.0, 0.0, 1.0};
  helper_set_normal_from_vector(other_normal, zn);

  // EDGE
  //
  // COORD 1
  double m1[3] = {1.0, 0.0, 0.0};
  helper_set_normal_from_vector(normals[0], m1);

  double m2[3] = {0.0, 1.0, 0.0};
  helper_set_normal_from_vector(normals[1], m2);

  goma_error error;
  int type;
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[1]));
  // one normal should be z-vector
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 2
  double m3[3] = {INV_SQRT_2, INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[0], m3);

  double m4[3] = {-INV_SQRT_2, INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[1], m4);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[1]));
  // one normal should be z-vector
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 3
  double m5[3] = {INV_SQRT_2, INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[0], m5);

  double m6[3] = {INV_SQRT_2, -INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[1], m6);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[1]));
  // one normal should be minus z-vector
  zn[2] = -1.0;
  helper_set_normal_from_vector(other_normal, zn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 4
  double m7[3] = {INV_SQRT_2, INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[0], m7);

  double m8[3] = {INV_SQRT_2 + 0.2, -INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[1], m8);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  // not perfect coordinate normals, now need to use near check
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));
  // one normal should be z-vector
  zn[2] = 1.0;
  helper_set_normal_from_vector(other_normal, zn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 5
  double m9[3] = {INV_SQRT_2, INV_SQRT_2, 0.0};
  helper_set_normal_from_vector(normals[0], m9);

  double m10[3] = {INV_SQRT_2 + 0.2, -INV_SQRT_2 - 0.3, 0.0};
  helper_set_normal_from_vector(normals[1], m10);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  // not perfect coordinate normals, now need to use near check
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));
  // one normal should be minus z-vector
  zn[2] = -1.0;
  helper_set_normal_from_vector(other_normal, zn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // Test a couple that aren't rotating only x,y
  // COORD 6
  double m11[3] = {0.0, INV_SQRT_2, INV_SQRT_2};
  helper_set_normal_from_vector(normals[0], m11);

  double m12[3] = {0.0, INV_SQRT_2, -INV_SQRT_2};
  helper_set_normal_from_vector(normals[1], m12);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  // not perfect coordinate normals, now need to use near check
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[1]));
  // one normal should be x-vector
  double xn[3] = {1.0, 0.0, 0.0};
  helper_set_normal_from_vector(other_normal, xn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 7
  double m13[3] = {0.0, INV_SQRT_2, INV_SQRT_2 + 0.2};
  helper_set_normal_from_vector(normals[0], m13);

  double m14[3] = {0.0, -INV_SQRT_2 - 0.3, INV_SQRT_2};
  helper_set_normal_from_vector(normals[1], m14);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  // not perfect coordinate normals, now need to use near check
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));
  // one normal should be minus x-vector
  xn[0] = -1.0;
  helper_set_normal_from_vector(other_normal, xn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  // COORD 8
  double m15[3] = {INV_SQRT_2, 0.0, INV_SQRT_2 + 0.2};
  helper_set_normal_from_vector(normals[0], m15);

  double m16[3] = {-INV_SQRT_2 - 0.3, 0.0, INV_SQRT_2};
  helper_set_normal_from_vector(normals[1], m16);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);
  REQUIRE(error != GOMA_ERROR);

  // not perfect coordinate normals, now need to use near check
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[0]));
  REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[1]));
  // one normal should be minus y-vector
  double yn[3] = {0.0, -1.0, 0.0};
  helper_set_normal_from_vector(other_normal, yn);
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, other_normal));

  goma_normal_free(other_normal);
  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}

TEST_CASE("goma_best_coordinate_system_3D varying normals surface", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 8;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  // Generic surface cases
  // All same
  //
  double xn[3] = {1.0, 0.0, 0.0};
  double yn[3] = {0.0, 1.0, 0.0};
  double zn[3] = {0.0, 0.0, 1.0};

  // Set 1
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], xn);
  }

  goma_error error;
  int type;
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  // Set 2
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  // Set 3
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[0]));

  // Now some normals that are near but not the same
  // Set 4
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], xn);
  }

  double m1[3] = {0.90453403, 0.30151134, 0.30151134};
  double m2[3] = {0.90453403, -0.30151134, 0.30151134};
  double m3[3] = {0.90453403, -0.30151134, -0.30151134};
  double m4[3] = {0.90453403, 0.30151134, -0.30151134};
  helper_set_normal_from_vector(normals[0], m1);
  helper_set_normal_from_vector(normals[1], m2);
  helper_set_normal_from_vector(normals[2], m3);
  helper_set_normal_from_vector(normals[3], m4);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < n_normals; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 5
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }
  double m5[3] = {0.30151134, 0.90453403, 0.30151134};
  double m6[3] = {-0.30151134, 0.90453403, 0.30151134};
  double m7[3] = {-0.30151134, 0.90453403, -0.30151134};
  double m8[3] = {0.30151134, 0.90453403, -0.30151134};
  helper_set_normal_from_vector(normals[0], m5);
  helper_set_normal_from_vector(normals[1], m6);
  helper_set_normal_from_vector(normals[2], m7);
  helper_set_normal_from_vector(normals[3], m8);

  error = goma_best_coordinate_system_3D(normals, 5, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 5; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 6
  yn[1] = -1.0;
  for (int i = 0; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }
  double m9[3] = {0.30151134, -0.90453403, 0.30151134};
  double m10[3] = {0.34815531, -0.87038828, 0.34815531};
  double m11[3] = {-0.30151134, 0.90453403, -0.30151134};
  double m12[3] = {-0.34815531, -0.87038828, 0.34815531};
  helper_set_normal_from_vector(normals[0], m9);
  helper_set_normal_from_vector(normals[1], m10);
  helper_set_normal_from_vector(normals[2], m11);
  helper_set_normal_from_vector(normals[3], m12);

  error = goma_best_coordinate_system_3D(normals, 5, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);
  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));

  for (int i = 0; i < 5; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }
  // Set 7

  double m13[3] = {0.30151134, 0.30151134, -0.90453403};
  double m14[3] = {0.34815531, 0.34815531, -0.87038828};
  double m15[3] = {0, 0, -1.0};
  helper_set_normal_from_vector(normals[0], m13);
  helper_set_normal_from_vector(normals[1], m14);
  helper_set_normal_from_vector(normals[2], m15);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // cleanup

  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}

TEST_CASE("goma_best_coordinate_system_3D varying normals edge", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 8;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  // Generic edge cases
  // All same
  //
  double xn[3] = {1.0, 0.0, 0.0};
  double yn[3] = {0.0, 1.0, 0.0};
  double zn[3] = {0.0, 0.0, 1.0};

  // Set 1
  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], xn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  goma_error error;
  int type;
  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < n_normals; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 2
  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < n_normals; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 3
  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], xn);
  }

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < n_normals; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Varying edges, not all same direction
  // Set 3
  //

  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }
  double m1[3] = {0.0, INV_SQRT_2 + 0.3, INV_SQRT_2};
  double m2[3] = {0.0, INV_SQRT_2, INV_SQRT_2 + 0.3};
  double m3[3] = {0.0, INV_SQRT_2 + 0.2, INV_SQRT_2};
  double m4[3] = {0.0, INV_SQRT_2, INV_SQRT_2 + 0.2};
  helper_set_normal_from_vector(normals[0], m1);
  helper_set_normal_from_vector(normals[1], m2);
  helper_set_normal_from_vector(normals[4], m3);
  helper_set_normal_from_vector(normals[5], m4);

  error = goma_best_coordinate_system_3D(normals, n_normals, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < n_normals; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 4
  // same as 3 but only 1 is on another edge

  error = goma_best_coordinate_system_3D(normals, 5, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 5; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 5
  // 3 normals

  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }
  double m5[3] = {3, 1, 1};
  double m6[3] = {3.5, 1, 1};
  double m7[3] = {0, 3.5, 0.5};
  helper_set_normal_from_vector(normals[0], m5);
  helper_set_normal_from_vector(normals[1], m6);
  helper_set_normal_from_vector(normals[2], m7);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 5
  // 4 normals

  for (int i = 0; i < 4; i++) {
    helper_set_normal_from_vector(normals[i], yn);
  }

  for (int i = 4; i < n_normals; i++) {
    helper_set_normal_from_vector(normals[i], zn);
  }
  double m8[3] = {0, 3, 0.5};
  helper_set_normal_from_vector(normals[3], m8);

  error = goma_best_coordinate_system_3D(normals, 4, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 4; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // cleanup
  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}

TEST_CASE("goma_best_coordinate_system_3D varying normals corner", "[bc][automatic_rotations]") {
  goma_normal **normals;
  int n_normals = 8;
  helper_create_normals(&normals, n_normals);

  goma_normal *rotated_coord[3];
  for (int j = 0; j < 3; j++) {
    rotated_coord[j] = goma_normal_alloc(3);
  }

  // corner base cases
  double xn[3] = {1.0, 0.0, 0.0};
  double yn[3] = {0.0, 1.0, 0.0};
  double zn[3] = {0.0, 0.0, 1.0};

  // Set 1
  helper_set_normal_from_vector(normals[0], xn);
  helper_set_normal_from_vector(normals[1], yn);
  helper_set_normal_from_vector(normals[2], zn);

  goma_error error;
  int type;
  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 2
  // permuted set 1
  helper_set_normal_from_vector(normals[0], zn);
  helper_set_normal_from_vector(normals[1], yn);
  helper_set_normal_from_vector(normals[2], xn);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 3
  // permuted set 1
  helper_set_normal_from_vector(normals[0], yn);
  helper_set_normal_from_vector(normals[1], zn);
  helper_set_normal_from_vector(normals[2], xn);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 4
  // permuted set 1
  helper_set_normal_from_vector(normals[0], yn);
  helper_set_normal_from_vector(normals[1], xn);
  helper_set_normal_from_vector(normals[2], zn);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_matches_normal(rotated_coord, normals[i]));
  }

  // Set 5
  double m1[3] = {INV_SQRT_2 + 0.5, INV_SQRT_2, 0};
  double m2[3] = {INV_SQRT_2, -INV_SQRT_2 - 0.5, 0};
  double m3[3] = {0, 0.2, -INV_SQRT_2};

  helper_set_normal_from_vector(normals[0], m1);
  helper_set_normal_from_vector(normals[1], m2);
  helper_set_normal_from_vector(normals[2], m3);

  error = goma_best_coordinate_system_3D(normals, 3, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 3; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  // Set 6
  helper_set_normal_from_vector(normals[0], m1);
  helper_set_normal_from_vector(normals[1], m2);
  helper_set_normal_from_vector(normals[2], m3);
  helper_set_normal_from_vector(normals[3], zn);
  helper_set_normal_from_vector(normals[4], xn);

  error = goma_best_coordinate_system_3D(normals, 5, rotated_coord, &type);

  REQUIRE(error != GOMA_ERROR);

  REQUIRE(helper_coordinate_is_coordinate_system(rotated_coord));
  for (int i = 0; i < 5; i++) {
    REQUIRE(helper_one_coordinate_near_normal(rotated_coord, normals[i]));
  }

  helper_free_normals(&normals, n_normals);
  for (int j = 0; j < 3; j++) {
    goma_normal_free(rotated_coord[j]);
  }
}
