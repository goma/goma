#include <catch2/catch.hpp>

extern "C" {
#include "gds/gds_vector.h"
#include "bc/rotate_util.h"
}

static const auto zero_comp = Catch::Floating::WithinAbsMatcher(0.0, 1e-15);

TEST_CASE("rotate coordinates for single normal", "[bc][rotate]") {
  gds_vector *coord[3];
  gds_vector *normal;

  coord[0] = gds_vector_alloc(3);
  coord[1] = gds_vector_alloc(3);
  coord[2] = gds_vector_alloc(3);
  normal = gds_vector_alloc(3);

  gds_vector_zero(normal);
  gds_vector_set(normal, 0, 1);

  unsigned int normal_dirs[3];

  goma_error error = goma_rotated_coordinates_from_normals(&normal, 1, coord, normal_dirs);
  REQUIRE(error == GOMA_SUCCESS);
  REQUIRE(normal_dirs[0] == 0);

  REQUIRE(gds_vector_get(coord[0], 0) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 1)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 0)));
  REQUIRE(gds_vector_get(coord[1], 1) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 0)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 1)));
  REQUIRE(gds_vector_get(coord[2], 2) == Approx(1.0));

  gds_vector_zero(normal);
  gds_vector_set(normal, 2, 1);

  error = goma_rotated_coordinates_from_normals(&normal, 1, coord, normal_dirs);
  REQUIRE(error == GOMA_SUCCESS);
  REQUIRE(normal_dirs[0] == 2);

  REQUIRE(gds_vector_get(coord[0], 0) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 1)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 0)));
  REQUIRE(gds_vector_get(coord[1], 1) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 0)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 1)));
  REQUIRE(gds_vector_get(coord[2], 2) == Approx(1.0));


  double inv_root_2 = 1.0 / sqrt(2.0);
  gds_vector_zero(normal);
  gds_vector_set(normal, 0, inv_root_2);
  gds_vector_set(normal, 1, inv_root_2);

  error = goma_rotated_coordinates_from_normals(&normal, 1, coord, normal_dirs);
  error = goma_rotated_coordinates_from_normals(&normal, 1, coord, normal_dirs);
  REQUIRE(error == GOMA_SUCCESS);
  REQUIRE(normal_dirs[0] == 2);

  REQUIRE(gds_vector_get(coord[0], 0) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 1)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[0], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 0)));
  REQUIRE(gds_vector_get(coord[1], 1) == Approx(1.0));
  REQUIRE(zero_comp.match(gds_vector_get(coord[1], 2)));

  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 0)));
  REQUIRE(zero_comp.match(gds_vector_get(coord[2], 1)));
  REQUIRE(gds_vector_get(coord[2], 2) == Approx(1.0));
}
