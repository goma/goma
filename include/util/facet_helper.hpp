/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers                                      *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#ifndef GOMA_FACET_HELPER_HPP
#define GOMA_FACET_HELPER_HPP

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <ostream>
#include <vector>
#include <algorithm>

namespace goma {
namespace distance_tools {

template <int dim> class Point {
  std::array<double, dim> data;

public:
  Point() : data({}) {}
  Point(std::array<double, dim> data) : data(data) {}
  Point(double x, double y) : data({}) {
    data[0] = x;
    data[1] = y;
  }

  Point(const Point<dim> &other) : data(other.data) {}


  Point<dim>& operator=(const Point<dim>&) = default;
  Point<dim>& operator=(Point<dim>&&) = default;

  double operator[](int i) const { return data[i]; }
  double &operator[](int i) { return data[i]; }
};

template <int dim> constexpr Point<dim> operator-(const Point<dim> &a, const Point<dim> &b) {
  Point<dim> result;
  for (int i = 0; i < dim; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

template <int dim> constexpr Point<dim> operator+(const Point<dim> &a, const Point<dim> &b) {
  Point<dim> result;
  for (int i = 0; i < dim; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

template <int dim> constexpr Point<dim> operator*(const Point<dim> &a, const double &b) {
  Point<dim> result;
  for (int i = 0; i < dim; ++i) {
    result[i] = a[i] * b;
  }
  return result;
}

template <int dim> constexpr Point<dim> operator*(const double &a, const Point<dim> &b) {
  Point<dim> result;
  for (int i = 0; i < dim; ++i) {
    result[i] = b[i] * a;
  }
  return result;
}

template <int dim> constexpr double dot_product(const Point<dim> &a, const Point<dim> &b) {
  double result = 0;
  for (int i = 0; i < dim; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

template <int dim> constexpr Point<dim> cross(const Point<dim> &a, const Point<dim> &b) {
  static_assert(dim == 3, "Cross product is only defined for 3D vectors");
  return Point<dim>({a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]});
}

template <int dim> constexpr double squared_length(const Point<dim> &a) {
  return dot_product(a, a);
}

template <int dim> constexpr double distance(const Point<dim> &a, const Point<dim> &b) {
  return std::sqrt(squared_length(a - b));
}

template <int dim> class Facet {
  // virtual function for distance
public:
  virtual double distance(const Point<dim> &p) const = 0;
};

template <int dim> class Line : public Facet<dim> {
public:
  const Point<dim> p0, p1;

  Line(const Point<dim> &p0, const Point<dim> &p1)
      : p0(p0), p1(p1) {}

  Point<dim> closest_point(const Point<dim> &p) const {
    Point<dim> v = p1 - p0;
    Point<dim> u = p0 - p;

    double dot = dot_product(u, v);
    if (std::abs(dot) < 10 * std::numeric_limits<double>::epsilon()) { // degenerate case where p0
                                                                       // and p1 are the same point
      return Point<dim>(p0);
    }

    double length_sq = squared_length(v);

    Point<dim> closest_point;

    double projection_location = -dot / length_sq;

    if (projection_location < 0) {
      closest_point = p0;
    } else if (projection_location > 1) {
      closest_point = p1;
    } else {
      closest_point = (1 - projection_location) * p0 + projection_location * p1;
    }

    return closest_point;
  }

  double distance(const Point<dim> &p) const override {
    // calculate distance from point p to line segment p0p1

    auto closest = closest_point(p);

    return std::sqrt(squared_length(closest - p));
  }
};

template <int dim>
constexpr Point<dim> linear_interp(const Point<dim> &p1,
                                   const Point<dim> &p2,
                                   double value1,
                                   double value2,
                                   double isoval,
                                   double tolerance = 1e-15) {
  if (std::abs(isoval - value1) < tolerance) {
    return p1;
  }
  if (std::abs(isoval - value2) < tolerance) {
    return p2;
  }
  double scale = (isoval - value1) / (value2 - value1);

  Point<dim> result = p1 + scale*(p2 - p1);
  for (int i = 0; i < dim; ++i) {
    if (result[i] < -1 || result[i] > 1) {
      throw std::runtime_error("Interpolated point is outside of the unit cube");
    }
  }
  return result;
}

constexpr double sign(double x) {
  if (x < 0) return -1;
  return 1;
}

template <int dim> class Triangle : public Facet<dim> {
public:
  const Point<dim> p0, p1, p2;

  Triangle(const Point<dim> &p0, const Point<dim> &p1, const Point<dim> &p2)
      : p0(p0), p1(p1), p2(p2) {}

  Point<3> closest_point(const Point<dim> &p) const {
    Point<3> normal = cross(p1 - p0, p2 - p0);
    double length = std::sqrt(squared_length(normal));
    if (length < 1e-15) { // degenerate case where p0 is the same as p1 or p2 
      Line<dim> l(p1, p2);
      return l.closest_point(p);
    }
    normal = normal * (1.0 / std::sqrt(squared_length(normal)));

    double t = dot_product(normal, p0) - dot_product(normal, p);
    Point<3> intersection = p + t * normal;

    // check if intersection is inside the triangle
    Point<3> edge0 = p1 - p0;
    Point<3> edge1 = p2 - p1;
    Point<3> edge2 = p0 - p2;

    Point<3> c0 = intersection - p0;
    Point<3> c1 = intersection - p1;
    Point<3> c2 = intersection - p2;

    if (dot_product(cross(edge0, c0), normal) > 0 &&
        dot_product(cross(edge1, c1), normal) > 0 &&
        dot_product(cross(edge2, c2), normal) > 0) {
      return intersection;
    }

    // find closest point on each edge
    Line<dim> l0(p0, p1);
    Line<dim> l1(p1, p2);
    Line<dim> l2(p2, p0);

    Point<3> closest = l0.closest_point(p);
    double min_distance = squared_length(closest - p);

    Point<3> closest1 = l1.closest_point(p);
    double distance1 = squared_length(closest1 - p);
    if (distance1 < min_distance) {
      min_distance = distance1;
      closest = closest1;
    }

    Point<3> closest2 = l2.closest_point(p);
    double distance2 = squared_length(closest2 - p);
    if (distance2 < min_distance) {
      min_distance = distance2;
      closest = closest2;
    }

    return closest;
  }

  double distance(const Point<dim> &p) const override {
    auto closest = closest_point(p);
    return std::sqrt(squared_length(closest - p));

    // prepare data    
    // Point<3> v21 = p1 - p0; Point<3> v1 = p - p0;
    // Point<3> v32 = p2 - p1; Point<3> v2 = p - p1;
    // Point<3> v13 = p0 - p2; Point<3> v3 = p - p2;

    // double lv32 = squared_length(v32);
    // double lv13 = squared_length(v13);
    // double lv21 = squared_length(v21);
    // if (lv32 < 1e-15) {
    //   Line<dim> l(p0,p1);
    //   return l.distance(p);
    // }

    // if (lv13 < 1e-15) {
    //   Line<dim> l(p1,p2);
    //   return l.distance(p);
    // }

    // if (lv21 < 1e-15) {
    //   Line<dim> l(p2,p0);
    //   return l.distance(p);
    // }

    // Point<3> nor = cross( v21, v13 );

    // return sqrt( // inside/outside test    
    //              (sign(dot_product(cross(v21,nor),v1)) + 
    //               sign(dot_product(cross(v32,nor),v2)) + 
    //               sign(dot_product(cross(v13,nor),v3))<2.0) 
    //               ?
    //               // 3 edges    
    //               std::min( std::min( 
    //               squared_length(v21*std::clamp(dot_product(v21,v1)/squared_length(v21),0.0,1.0)-v1), 
    //               squared_length(v32*std::clamp(dot_product(v32,v2)/squared_length(v32),0.0,1.0)-v2) ), 
    //               squared_length(v13*std::clamp(dot_product(v13,v3)/squared_length(v13),0.0,1.0)-v3) )
    //               :
    //               // 1 face    
    //               dot_product(nor,p1)*dot_product(nor,p1)/squared_length(nor) );
  }

};

/*
 *
 *    3
 *     o
 *     * *
 *     *   *
 *     *     *
 *     *       *
 *     *         *
 *     o * * * * * o
 *    2               1
 */
template <int dim>
std::optional<Line<dim>> create_facet_from_triangle(const Point<dim> &p1,
                                                    const Point<dim> &p2,
                                                    const Point<dim> &p3,
                                                    std::array<double, 3> values,
                                                    double isoval) {

  // We have 8 cases for the triangle

  int cond = 0;
  for (int i = 0; i < 3; ++i) {
    if (values[i] >= isoval) {
      cond |= 1 << i;
    }
  }

  switch (cond) {
  case 0b000:
  case 0b111:
    return std::nullopt;
  // intersection on p1->p2 and p1->p3
  case 0b001:
  case 0b110:
    return Line<dim>(linear_interp(p1, p2, values[0], values[1], isoval),
                     linear_interp(p1, p3, values[0], values[2], isoval));
  // intersection on p1->p2 and p2->p3
  case 0b010:
  case 0b101:
    return Line<dim>(linear_interp(p1, p2, values[0], values[1], isoval),
                     linear_interp(p2, p3, values[1], values[2], isoval));
  // intersection on p2->p3 and p1->p3
  case 0b011:
  case 0b100:
    return Line<dim>(linear_interp(p1, p3, values[0], values[2], isoval),
                     linear_interp(p2, p3, values[1], values[2], isoval));
  default:
    return std::nullopt;
  }
}

/*
 *
 *
 *    4              3
 *     o * * * * * o
 *     *           *
 *     *           *
 *     *           *
 *     *           *
 *     o * * * * * o
 *    1              2
 */
template <int dim>
std::vector<Line<dim>> create_facet_from_quad(const Point<dim> &p0,
                                              const Point<dim> &p1,
                                              const Point<dim> &p2,
                                              const Point<dim> &p3,
                                              const std::array<double, 4> &values,
                                              double isoval) {

  // We have 16 cases for the square

  int cond = 0;
  for (int i = 0; i < 4; ++i) {
    if (values[i] > isoval) {
      cond |= 1 << i;
    }
  }

  std::vector<Line<dim>> result;

  // https://en.wikipedia.org/wiki/Marching_squares
  switch (cond) {
  case 0b0001: // 1
  case 0b1110: // 14
    result.push_back(Line<2>(linear_interp(p0, p1, values[0], values[1], isoval),
                             linear_interp(p0, p3, values[0], values[3], isoval)));
    break;
  case 0b0010: // 2
  case 0b1101: // 13
    result.push_back(Line<2>(linear_interp(p0, p1, values[0], values[1], isoval),
                             linear_interp(p1, p2, values[1], values[2], isoval)));
    break;
  case 0b0011: // 3
  case 0b1100: // 12
    result.push_back(Line<2>(linear_interp(p0, p3, values[0], values[3], isoval),
                             linear_interp(p1, p2, values[1], values[2], isoval)));
    break;
  case 0b0100: // 4
  case 0b1011: // 11
    result.push_back(Line<2>(linear_interp(p1, p2, values[1], values[2], isoval),
                             linear_interp(p2, p3, values[2], values[3], isoval)));
    break;
  case 0b0101: // 5
    result.push_back(Line<2>(linear_interp(p0, p1, values[0], values[1], isoval),
                             linear_interp(p1, p2, values[1], values[2], isoval)));
    result.push_back(Line<2>(linear_interp(p0, p3, values[0], values[3], isoval),
                             linear_interp(p2, p3, values[2], values[3], isoval)));
    break;
  case 0b1010: // 10
    result.push_back(Line<2>(linear_interp(p0, p1, values[0], values[1], isoval),
                             linear_interp(p0, p3, values[0], values[3], isoval)));
    result.push_back(Line<2>(linear_interp(p1, p2, values[1], values[2], isoval),
                             linear_interp(p2, p3, values[2], values[3], isoval)));
    break;
  case 0b0110: // 6
  case 0b1001: // 9
    result.push_back(Line<2>(linear_interp(p0, p1, values[0], values[1], isoval),
                             linear_interp(p2, p3, values[2], values[3], isoval)));
    break;
  case 0b0111: // 7
  case 0b1000: // 8
    result.push_back(Line<2>(linear_interp(p0, p3, values[0], values[3], isoval),
                             linear_interp(p2, p3, values[2], values[3], isoval)));
    break;
  case 0b0000: // 0
  case 0b1111: // 15
  default:
    break;
  }
  return result;
}

std::vector<Triangle<3>> create_facet_from_hex(const std::array<Point<3>, 8> &points,
                                               const std::array<double, 8> &values,
                                               double isoval);
} // namespace distance_tools
} // namespace goma

template <int dim>
std::ostream &operator<<(std::ostream &os, const goma::distance_tools::Point<dim> &p) {
  os << "(";
  for (int i = 0; i < dim; ++i) {
    os << p[i];
    if (i != dim - 1) {
      os << ", ";
    }
  }
  os << ")";
  return os;
}

template <int dim>
std::ostream &operator<<(std::ostream &os, const goma::distance_tools::Line<dim> &p) {
  os << "Line[" << p.p0 << " -> " << p.p1 << "]";
  return os;
}

#endif // GOMA_FACET_HELPER_HPP