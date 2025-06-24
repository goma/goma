#include "facet_helper.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

const double test_atol = 1e-12;
const double test_rtol = 1e-5;

#define REQUIRE_TOL(a, b)                                     \
  REQUIRE_THAT(a, Catch::Matchers::WithinRel(b, test_rtol) || \
                      Catch::Matchers::WithinAbs(b, test_atol))

TEST_CASE("point math", "[Point]") {
  using namespace goma::distance_tools;

  auto p1 = Point<2>({1, 2});
  auto p2 = Point<2>({3, 4});

  auto p3 = p1 + p2;

  REQUIRE_TOL(p3[0], 4.0);
  REQUIRE_TOL(p3[1], 6.0);

  auto p4 = p1 - p2;
  REQUIRE_TOL(p4[0], -2.0);
  REQUIRE_TOL(p4[1], -2.0);

  auto p5 = p1 * 2;
  REQUIRE_TOL(p5[0], 2.0);
  REQUIRE_TOL(p5[1], 4.0);

  auto p6 = dot_product(p1, p2);
  REQUIRE_TOL(p6, 11.0);
}

TEST_CASE("linear_interp", "[Point]") {
  using namespace goma::distance_tools;

  auto p0 = Point<2>({0, 0});
  auto p1 = Point<2>({0, 1});
  auto p2 = Point<2>({1, 0});

  auto p3 = linear_interp(p0, p1, 0, 1, 0.5);
  REQUIRE_TOL(p3[0], 0.0);
  REQUIRE_TOL(p3[1], 0.5);

  p3 = linear_interp(p0, p2, 0, 1, 0.5);
  REQUIRE_TOL(p3[0], 0.5);
  REQUIRE_TOL(p3[1], 0.0);

  p3 = linear_interp(p0, p2, 0, 1, 0.75);
  REQUIRE_TOL(p3[0], 0.75);
  REQUIRE_TOL(p3[1], 0.0);
}

TEST_CASE("line distance", "[Line]") {
  using namespace goma::distance_tools;

  auto p0 = Point<2>({0, 0});
  auto p1 = Point<2>({1, 1});
  auto p2 = Point<2>({2, 2});

  auto l1 = Line<2>(p0, p1);
  auto l2 = Line<2>(p0, p2);

  auto p3 = Point<2>({0.5, 0.5});
  auto p4 = Point<2>({1.5, 1.5});

  REQUIRE_TOL(l1.distance(p3), 0.0);
  REQUIRE_TOL(l1.distance(p4), distance(p4, p1));

  REQUIRE_TOL(l2.distance(p3), 0.0);
  REQUIRE_TOL(l2.distance(p4), 0.0);

  auto p5 = Point<2>({0.5, 0.0});
  auto p6 = Point<2>({1.5, 0.0});
}

TEST_CASE("line from triangle", "[Line]") {
  using namespace goma::distance_tools;

  auto p0 = Point<2>({1, 0});
  auto p1 = Point<2>({0, 0});
  auto p2 = Point<2>({0, 1});

  std::array<double, 3> values = {0, 0, 0};

  auto f1 = create_facet_from_triangle(p0, p1, p2, values, 0.5);

  REQUIRE(f1 == std::nullopt);

  values = {0, 0, 1};
  auto f2 = create_facet_from_triangle(p0, p1, p2, values, 0.5);

  std::cout << f2.value() << std::endl;
  REQUIRE(f2.has_value());
  auto l2 = f2.value();
  REQUIRE_TOL(l2.p0[0], 0.5);
  REQUIRE_TOL(l2.p0[1], 0.5);
  REQUIRE_TOL(l2.p1[0], 0.0);
  REQUIRE_TOL(l2.p1[1], 0.5);

  // shift towards p1
  values = {0, 0.75, 1};
  auto f3 = create_facet_from_triangle(p0, p1, p2, values, 0.75);

  std::cout << f3.value() << std::endl;
  REQUIRE(f3.has_value());
  auto l3 = f3.value();
  REQUIRE_TOL(l3.p0[0], 0.0);
  REQUIRE_TOL(l3.p0[1], 0.0);
  REQUIRE_TOL(l3.p1[0], 0.25);
  REQUIRE_TOL(l3.p1[1], 0.75);
}
