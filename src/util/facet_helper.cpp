#include "util/facet_helper.hpp"
#include "util/marching_cubes.hpp"
#include <iostream>

namespace goma {
namespace distance_tools {
std::vector<Triangle<3>> create_facet_from_hex(const std::array<Point<3>, 8> &points,
                                               const std::array<double, 8> &values,
                                               double isoval) {

  using namespace goma::marching_cubes;

  // rotate the points for the lookup table
  std::array<int, 8> rotation = {3, 2, 1, 0, 7, 6,5, 4};
  std::array<Point<3>, 8> rotated_points;
  std::array<double, 8> rotated_values;

  for (int i = 0; i < 8; ++i) {
    rotated_points[i] = points[rotation[i]];
    rotated_values[i] = values[rotation[i]];
  }

  // We have 256 cases for the cube
  int cond = 0;
  if (rotated_values[0] < isoval) cond |= 1;
  if (rotated_values[1] < isoval) cond |= 2;
  if (rotated_values[2] < isoval) cond |= 4;
  if (rotated_values[3] < isoval) cond |= 8;
  if (rotated_values[4] < isoval) cond |= 16;
  if (rotated_values[5] < isoval) cond |= 32;
  if (rotated_values[6] < isoval) cond |= 64;
  if (rotated_values[7] < isoval) cond |= 128;


  std::array<Point<3>, 12> vertlist;

  // for each edge in the edge table add matching triangles
  if (edge_table[cond] & 1)
    vertlist[0] = linear_interp(rotated_points[0], rotated_points[1], rotated_values[0],
                                rotated_values[1], isoval);
  if (edge_table[cond] & 2)
    vertlist[1] = linear_interp(rotated_points[1], rotated_points[2], rotated_values[1],
                                rotated_values[2], isoval);
  if (edge_table[cond] & 4)
    vertlist[2] = linear_interp(rotated_points[2], rotated_points[3], rotated_values[2],
                                rotated_values[3], isoval);
  if (edge_table[cond] & 8)
    vertlist[3] = linear_interp(rotated_points[3], rotated_points[0], rotated_values[3],
                                rotated_values[0], isoval);
  if (edge_table[cond] & 16)
    vertlist[4] = linear_interp(rotated_points[4], rotated_points[5], rotated_values[4],
                                rotated_values[5], isoval);
  if (edge_table[cond] & 32)
    vertlist[5] = linear_interp(rotated_points[5], rotated_points[6], rotated_values[5],
                                rotated_values[6], isoval);
  if (edge_table[cond] & 64)
    vertlist[6] = linear_interp(rotated_points[6], rotated_points[7], rotated_values[6],
                                rotated_values[7], isoval);
  if (edge_table[cond] & 128)
    vertlist[7] = linear_interp(rotated_points[7], rotated_points[4], rotated_values[7],
                                rotated_values[4], isoval);
  if (edge_table[cond] & 256)
    vertlist[8] = linear_interp(rotated_points[0], rotated_points[4], rotated_values[0],
                                rotated_values[4], isoval);
  if (edge_table[cond] & 512)
    vertlist[9] = linear_interp(rotated_points[1], rotated_points[5], rotated_values[1],
                                rotated_values[5], isoval);
  if (edge_table[cond] & 1024)
    vertlist[10] = linear_interp(rotated_points[2], rotated_points[6], rotated_values[2],
                                 rotated_values[6], isoval);
  if (edge_table[cond] & 2048)
    vertlist[11] = linear_interp(rotated_points[3], rotated_points[7], rotated_values[3],
                                 rotated_values[7], isoval);

  std::vector<Triangle<3>> result;

   for (int i=0;tri_table[cond][i]!=-1;i+=3) {
      result.push_back(Triangle<3>(vertlist[tri_table[cond][i]], vertlist[tri_table[cond][i+1]], vertlist[tri_table[cond][i+2]]));

      // check vertlist values for values outside -1, 1
      // if so print error
      for (int j=0;j<3;j++) {
         for (int k=0;k<3;k++) {
            if (vertlist[tri_table[cond][i+j]][k]<-1 || vertlist[tri_table[cond][i+j]][k]>1) {
               std::cerr << "Error: vertlist[" << tri_table[cond][i+j] << "][" << k << "]=" << vertlist[tri_table[cond][i+j]][k] << std::endl;
            }
         }
      }
   }

   return result;
}
} // namespace distance_tools
} // namespace goma