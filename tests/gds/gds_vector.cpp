#include <catch2/catch.hpp>
#include "gds/gds_vector.h"

static const auto zero_comp = Catch::Floating::WithinAbsMatcher(0.0, 1e-15);

TEST_CASE("gds vector Allocation", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(3);
  REQUIRE(v->size == 3);
  REQUIRE(v->data != NULL);
  gds_vector_free(v);
}

TEST_CASE("gds vector setting getting, setting, and copying", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(5);
  gds_vector *u = gds_vector_alloc(5);
  REQUIRE(v->size == 5);
  gds_vector_set(v, 0, 1);
  gds_vector_set(v, 1, -21131);
  gds_vector_set(v, 2, -0.89345);
  gds_vector_set(v, 3, 1e32);
  gds_vector_set(v, 4, -2.1e-32);
  REQUIRE(gds_vector_get(v,0) == 1);
  REQUIRE(gds_vector_get(v,1) == -21131);
  REQUIRE(gds_vector_get(v,2) == -0.89345);
  REQUIRE(gds_vector_get(v,3) == 1e32);
  REQUIRE(gds_vector_get(v,4) == -2.1e-32);
  gds_vector_copy(u, v);
  REQUIRE(gds_vector_get(u,0) == 1);
  REQUIRE(gds_vector_get(u,1) == -21131);
  REQUIRE(gds_vector_get(u,2) == -0.89345);
  REQUIRE(gds_vector_get(u,3) == 1e32);
  REQUIRE(gds_vector_get(u,4) == -2.1e-32);
  gds_vector_set_all(v, -3.0);
  gds_vector_set_all(u, 1e5);
  for (int i = 0; i < 5; i++) {
    REQUIRE(gds_vector_get(v,i) == -3.0);
    REQUIRE(gds_vector_get(u,i) == 1e5);
  }
  gds_vector_ones(v);
  for (int i = 0; i < 5; i++) {
    REQUIRE(gds_vector_get(v,i) == 1.0);
  }
  gds_vector_free(v);
  v = gds_vector_alloc(20);
  gds_vector_zero(v);
  for (int i = 0; i < 20; i++) {
    REQUIRE(gds_vector_get(v,i) == 0.0);
  }
  gds_vector_free(v);
  gds_vector_free(u);
}


TEST_CASE("gds_vector_add", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_add(u, v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 12);
    REQUIRE(gds_vector_get(v, i) == 2);
  }
  gds_vector_add(v, v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 4);
  }

  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -0.5);
  gds_vector_set(z, 2, 0.15);

  gds_vector_add(w, w);
  REQUIRE(gds_vector_get(w, 0) == 24);
  REQUIRE(gds_vector_get(w, 1) == -12);
  REQUIRE(gds_vector_get(w, 2) == -210);
  gds_vector_add(w, z);
  REQUIRE(gds_vector_get(w, 0) == 24.01);
  REQUIRE(gds_vector_get(w, 1) == -12.5);
  REQUIRE(gds_vector_get(w, 2) == -209.85);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_sub", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_sub(u,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 8);
    REQUIRE(gds_vector_get(v, i) == 2);
  }
  gds_vector_sub(v,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 0);
  }


  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -0.5);
  gds_vector_set(z, 2, 0.15);

  gds_vector_sub(w,z);
  REQUIRE(gds_vector_get(w,0) == 11.99);
  REQUIRE(gds_vector_get(w,1) == -5.5);
  REQUIRE(gds_vector_get(w,2) == -105.15);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_mul", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_mul(u,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 20);
    REQUIRE(gds_vector_get(v, i) == 2);
  }
  gds_vector_mul(v,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 4);
  }

  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -0.5);
  gds_vector_set(z, 2, 0.15);

  gds_vector_mul(w,z);
  REQUIRE(gds_vector_get(w,0) == 0.12);
  REQUIRE(gds_vector_get(w,1) == 3);
  REQUIRE(gds_vector_get(w,2) == -15.75);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_div", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_div(u,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 5);
    REQUIRE(gds_vector_get(v, i) == 2);
  }
  gds_vector_div(v,v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 1);
  }

  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -2.5);
  gds_vector_set(z, 2, 0.15);

  gds_vector_div(w,z);
  REQUIRE(gds_vector_get(w,0) == 1200);
  REQUIRE(gds_vector_get(w,1) == 2.4);
  REQUIRE(gds_vector_get(w,2) == -700);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}


TEST_CASE("gds_vector_scale", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_scale(u,12);
  gds_vector_scale(v,-8);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 120);
    REQUIRE(gds_vector_get(v, i) == -16);
  }
  gds_vector_scale(v,1);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == -16);
  }
  gds_vector_scale(v,0);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 0);
  }

  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_scale(w,0.5);
  REQUIRE(gds_vector_get(w,0) == 6);
  REQUIRE(gds_vector_get(w,1) == -3);
  REQUIRE(gds_vector_get(w,2) == -52.5);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
}

TEST_CASE("gds_vector_add_constant", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_add_constant(u,12);
  gds_vector_add_constant(v,-8);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(u, i) == 22);
    REQUIRE(gds_vector_get(v, i) == -6);
  }
  gds_vector_add_constant(v,1);
  gds_vector_add_constant(v,0);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == -5);
  }

  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_add_constant(w,0.5);
  REQUIRE(gds_vector_get(w,0) == 12.5);
  REQUIRE(gds_vector_get(w,1) == -5.5);
  REQUIRE(gds_vector_get(w,2) == -104.5);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
}

TEST_CASE("gds_vector_add_axpy", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(12);
  gds_vector *u = gds_vector_alloc(12);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_axpy(0, u, 2, v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 4);
  }
  gds_vector_axpy(0.5, u, 2, v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 13);
  }
  gds_vector_axpy(0.1, u, 0, v);
  for (int i = 0; i < 12; i++) {
    REQUIRE(gds_vector_get(v, i) == 1);
  }


  gds_vector_set(w, 0, 12);
  gds_vector_set(w, 1, -6);
  gds_vector_set(w, 2, -105);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -0.5);
  gds_vector_set(z, 2, 0.15);

  gds_vector_axpy(-0.1, w, 10, z);
  REQUIRE(gds_vector_get(z,0) == -1.1);
  REQUIRE(gds_vector_get(z,1) == -4.4);
  REQUIRE(gds_vector_get(z,2) == 12);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_add_normalize", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(4);
  gds_vector *u = gds_vector_alloc(4);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  gds_vector_normalize(u);
  gds_vector_normalize(v);
  for (int i = 0; i < 4; i++) {
    REQUIRE(gds_vector_get(v, i) == 0.5);
    REQUIRE(gds_vector_get(u, i) == 0.5);
  }

  gds_vector_set(w, 0, 0);
  gds_vector_set(w, 1, 0);
  gds_vector_set(w, 2, 1);

  gds_vector_set(z, 0, 3);
  gds_vector_set(z, 1, 5);
  gds_vector_set(z, 2, 4);

  gds_vector_normalize(w);
  gds_vector_normalize(z);

  REQUIRE(gds_vector_get(w,0) == 0);
  REQUIRE(gds_vector_get(w,1) == 0);
  REQUIRE(gds_vector_get(w,2) == 1);

  REQUIRE(gds_vector_get(z,0) == Approx(0.4242640687));
  REQUIRE(gds_vector_get(z,1) == Approx(0.7071067812));
  REQUIRE(gds_vector_get(z,2) == Approx(0.5656854249));

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_cross", "[gds][gds_vector]") {
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);
  gds_vector *cross = gds_vector_alloc(3);

  gds_vector_set(w, 0, 1);
  gds_vector_set(w, 1, 0);
  gds_vector_set(w, 2, 0);

  gds_vector_set(z, 0, 0);
  gds_vector_set(z, 1, 1);
  gds_vector_set(z, 2, 0);

  gds_vector_cross(w,z,cross);
  REQUIRE(gds_vector_get(cross,0) == 0);
  REQUIRE(gds_vector_get(cross,1) == 0);
  REQUIRE(gds_vector_get(cross,2) == Approx(1));

  gds_vector_cross(z,w,cross);
  REQUIRE(gds_vector_get(cross,0) == 0);
  REQUIRE(gds_vector_get(cross,1) == 0);
  REQUIRE(gds_vector_get(cross,2) == Approx(-1));

  gds_vector_cross(z,z,cross);
  REQUIRE(gds_vector_get(cross,0) == 0);
  REQUIRE(gds_vector_get(cross,1) == 0);
  REQUIRE(gds_vector_get(cross,2) == 0);

  gds_vector_set(w, 0, 3);
  gds_vector_set(w, 1, 4);
  gds_vector_set(w, 2, 5);

  gds_vector_set(z, 0, 5);
  gds_vector_set(z, 1, 1);
  gds_vector_set(z, 2, -8);

  gds_vector_cross(w,z,cross);
  REQUIRE(gds_vector_get(cross,0) == -37);
  REQUIRE(gds_vector_get(cross,1) == 49);
  REQUIRE(gds_vector_get(cross,2) == -17);

  gds_vector_free(cross);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_dot", "[gds][gds_vector]") {
  gds_vector *v = gds_vector_alloc(10);
  gds_vector *u = gds_vector_alloc(10);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set_all(u, 10);
  gds_vector_set_all(v, 2);

  REQUIRE(gds_vector_dot(u,v) == 200);
  REQUIRE(gds_vector_dot(v,u) == 200);
  REQUIRE(gds_vector_dot(u,u) == 1000);
  REQUIRE(gds_vector_dot(v,v) == 40);

  gds_vector_set(w, 0, 1);
  gds_vector_set(w, 1, -5);
  gds_vector_set(w, 2, -32);

  gds_vector_set(z, 0, 0.01);
  gds_vector_set(z, 1, -0.5);
  gds_vector_set(z, 2, 32.15);

  REQUIRE(gds_vector_dot(w,w) == 1050);
  REQUIRE(gds_vector_dot(z,z) == 1033.8726);
  REQUIRE(gds_vector_dot(w,z) == -1026.29);

  gds_vector_free(v);
  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}

TEST_CASE("gds_vector_rotate_around_vector", "[gds][gds_vector]") {
  gds_vector *u = gds_vector_alloc(3);
  gds_vector *w = gds_vector_alloc(3);
  gds_vector *z = gds_vector_alloc(3);

  gds_vector_set(w, 0, 1);
  gds_vector_set(w, 1, 0);
  gds_vector_set(w, 2, 0);

  gds_vector_set(z, 0, 0);
  gds_vector_set(z, 1, 1);
  gds_vector_set(z, 2, 0);

  gds_vector_rotate_around_vector(u, w, z, M_PI);
  REQUIRE(gds_vector_get(u,0) == Approx(-1));
  REQUIRE(zero_comp.match(gds_vector_get(u,1)));
  REQUIRE(zero_comp.match(gds_vector_get(u,2)));

  gds_vector_rotate_around_vector(u, z, w, M_PI);
  REQUIRE(zero_comp.match(gds_vector_get(u,0)));
  REQUIRE(gds_vector_get(u,1) == Approx(-1));
  REQUIRE(zero_comp.match(gds_vector_get(u,2)));

  gds_vector_rotate_around_vector(u, w, z, M_PI/2);
  REQUIRE(zero_comp.match(gds_vector_get(u,0)));
  REQUIRE(zero_comp.match(gds_vector_get(u,1)));
  REQUIRE(gds_vector_get(u,2) == Approx(-1));

  gds_vector_rotate_around_vector(u, w, z, -M_PI/2);
  REQUIRE(zero_comp.match(gds_vector_get(u,0)));
  REQUIRE(zero_comp.match(gds_vector_get(u,1)));
  REQUIRE(gds_vector_get(u,2) == Approx(1));

  gds_vector_free(u);
  gds_vector_free(w);
  gds_vector_free(z);
}
