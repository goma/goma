/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2025 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>
#include <mpi.h>
#include <nanoflann.hpp>
#include <numeric>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "ls/facet_based_reinit.h"
#include "util/facet_helper.hpp"

extern "C" {
#define DISABLE_CPP
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_fill_ptrs.h"
#include "mm_unknown_map.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#undef DISABLE_CPP
}

using namespace goma::distance_tools;

class BasicTimer {
  double start_;

public:
  BasicTimer() { start(); }
  void start() { start_ = MPI_Wtime(); }
  double elapsed() { return MPI_Wtime() - start_; }
  double elapsed_and_reset() {
    double elapsed = MPI_Wtime() - start_;
    start();
    return elapsed;
  }

  void print_elapsed(const char *msg) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcID == 0)
      std::cout << msg << " " << std::scientific << elapsed() << " sec" << std::endl;
  }

  void print_elapsed_and_reset(const char *msg) {
    print_elapsed(msg);
    start();
  }
};

void shape_to_triangle(int itype,
                       std::vector<std::tuple<Point<2>, Point<2>, Point<2>>> &triangles,
                       std::vector<std::tuple<int, int, int>> &order,
                       int dofs) {
  int shape = type2shape(itype);
  switch (shape) {
  case TRIANGLE:
    triangles.push_back(
        std::make_tuple(Point<2>(1.0, 0.0), Point<2>(0.0, 0.0), Point<2>(0.0, 1.0)));
    order.push_back({0, 1, 2});
    break;
  case QUADRILATERAL:
    if (dofs == 4) {
      triangles.push_back(
          std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(1.0, -1.0), Point<2>(1.0, 1.0)));
      order.push_back({0, 1, 2});
      triangles.push_back(
          std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), Point<2>(-1.0, 1.0)));
      order.push_back({0, 2, 3});
      triangles.push_back(
          std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(1.0, -1.0), Point<2>(-1.0, 1.0)));
      order.push_back({0, 1, 3});
      triangles.push_back(
          std::make_tuple(Point<2>(1.0, -1.0), Point<2>(1.0, 1.0), Point<2>(-1.0, 1.0)));
      order.push_back({1, 2, 3});
    } else if (dofs == 9) {
      triangles.push_back(
          std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(0.0, -1.0), Point<2>(0.0, 0.0)));
      order.push_back({0, 4, 8});
      triangles.push_back(
          std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(-1.0, 0.0), Point<2>(0.0, 0.0)));
      order.push_back({0, 7, 8});
      triangles.push_back(
          std::make_tuple(Point<2>(0.0, -1.0), Point<2>(1.0, -1.0), Point<2>(0.0, 0.0)));
      order.push_back({4, 1, 8});
      triangles.push_back(
          std::make_tuple(Point<2>(1.0, -1.0), Point<2>(1.0, 0.0), Point<2>(0.0, 0.0)));
      order.push_back({1, 5, 8});
      triangles.push_back(
          std::make_tuple(Point<2>(0.0, 0.0), Point<2>(1.0, 0.0), Point<2>(1.0, 1.0)));
      order.push_back({8, 5, 2});
      triangles.push_back(
          std::make_tuple(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), Point<2>(0.0, 1.0)));
      order.push_back({8, 2, 6});
      triangles.push_back(
          std::make_tuple(Point<2>(0.0, 0.0), Point<2>(0.0, 1.0), Point<2>(-1.0, 1.0)));
      order.push_back({8, 6, 3});
      triangles.push_back(
          std::make_tuple(Point<2>(0.0, 0.0), Point<2>(-1.0, 1.0), Point<2>(-1.0, 0.0)));
      order.push_back({8, 3, 7});
    } else {
      GOMA_EH(GOMA_ERROR, "Unsupported dofs %d for element type %d", dofs, itype);
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Unsupported element type %d", itype);
    break;
  }
}

void shape_to_quad(int itype,
                   std::vector<std::tuple<Point<2>, Point<2>, Point<2>, Point<2>>> &quads,
                   std::vector<std::tuple<int, int, int, int>> &order,
                   int dofs) {
  int shape = type2shape(itype);
  switch (shape) {
  case QUADRILATERAL:
    if (dofs == 4) {
      quads.push_back(std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(1.0, -1.0), Point<2>(1.0, 1.0),
                                      Point<2>(-1.0, 1.0)));
      order.push_back({0, 1, 2, 3});
    } else if (dofs == 9) {
      quads.push_back(std::make_tuple(Point<2>(-1.0, -1.0), Point<2>(0.0, -1.0), Point<2>(0.0, 0.0),
                                      Point<2>(-1.0, 0.0)));
      order.push_back({0, 4, 8, 7});
      quads.push_back(std::make_tuple(Point<2>(0.0, -1.0), Point<2>(1.0, -1.0), Point<2>(1.0, 0.0),
                                      Point<2>(0.0, 0.0)));
      order.push_back({4, 1, 5, 8});
      quads.push_back(std::make_tuple(Point<2>(0.0, 0.0), Point<2>(1.0, 0.0), Point<2>(1.0, 1.0),
                                      Point<2>(0.0, 1.0)));
      order.push_back({8, 5, 2, 6});
      quads.push_back(std::make_tuple(Point<2>(-1.0, 0.0), Point<2>(0.0, 0.0), Point<2>(0.0, 1.0),
                                      Point<2>(-1.0, 1.0)));
      order.push_back({7, 8, 6, 3});
    } else {
      GOMA_EH(GOMA_ERROR, "Unsupported dofs %d for element type %d", dofs, itype);
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Unsupported element type %d", itype);
    break;
  }
}

const std::array<std::array<double, 3>, 27> hex_points = {{{-1.0, -1.0, -1.0},
                                                           {1.0, -1.0, -1.0},
                                                           {
                                                               1.0,
                                                               1.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               1.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               -1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               1.0,
                                                               -1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               1.0,
                                                               1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               0.0,
                                                               -1.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               1.0,
                                                               0.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               0.0,
                                                               1.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               0.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               -1.0,
                                                               0.0,
                                                           },
                                                           {
                                                               1.0,
                                                               -1.0,
                                                               0.0,
                                                           },
                                                           {
                                                               1.0,
                                                               1.0,
                                                               0.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               1.0,
                                                               0.0,
                                                           },
                                                           {
                                                               0.0,
                                                               -1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               1.0,
                                                               0.0,
                                                               1.0,
                                                           },
                                                           {
                                                               0.0,
                                                               1.0,
                                                               1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               0.0,
                                                               1.0,
                                                           },
                                                           {
                                                               0.0,
                                                               0.0,
                                                               0.0,
                                                           },
                                                           {
                                                               0.0,
                                                               0.0,
                                                               -1.0,
                                                           },
                                                           {
                                                               0.0,
                                                               0.0,
                                                               1.0,
                                                           },
                                                           {
                                                               -1.0,
                                                               0.0,
                                                               0.0,
                                                           },
                                                           {
                                                               1.0,
                                                               0.0,
                                                               0.0,
                                                           },
                                                           {
                                                               0.0,
                                                               -1.0,
                                                               0.0,
                                                           },
                                                           {
                                                               0.0,
                                                               1.0,
                                                               0.0,
                                                           }}};

void shape_to_hex(int itype,
                  std::vector<std::array<Point<3>, 8>> &hexes,
                  std::vector<std::array<int, 8>> &order,
                  int dofs) {
  int shape = type2shape(itype);
  switch (shape) {
  case HEXAHEDRON:
    if (dofs == 8) {
      hexes.push_back({Point<3>(hex_points[0]), Point<3>(hex_points[1]), Point<3>(hex_points[2]),
                       Point<3>(hex_points[3]), Point<3>(hex_points[4]), Point<3>(hex_points[5]),
                       Point<3>(hex_points[6]), Point<3>(hex_points[7])});

      order.push_back({0, 1, 2, 3, 4, 5, 6, 7});
    } else if (dofs == 27) {
      hexes.push_back({Point<3>(hex_points[0]), Point<3>(hex_points[1]), Point<3>(hex_points[2]),
                       Point<3>(hex_points[3]), Point<3>(hex_points[4]), Point<3>(hex_points[5]),
                       Point<3>(hex_points[6]), Point<3>(hex_points[7])});

      order.push_back({0, 1, 2, 3, 4, 5, 6, 7});
    } else {
      GOMA_EH(GOMA_ERROR, "Unsupported dofs %d for element type %d", dofs, itype);
    }
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Unsupported element type %d", itype);
    break;
  }
}

void interp_quad_dofs(std::vector<double> &ls_values) {
  ls_values[4] = 0.5 * (ls_values[0] + ls_values[1]);
  ls_values[5] = 0.5 * (ls_values[1] + ls_values[2]);
  ls_values[6] = 0.5 * (ls_values[2] + ls_values[3]);
  ls_values[7] = 0.5 * (ls_values[3] + ls_values[0]);
  ls_values[8] = 0.25 * (ls_values[0] + ls_values[1] + ls_values[2] + ls_values[3]);
}

template <int pdim> struct PointCloud {
  std::vector<Point<pdim>> pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate
  // value, the
  //  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    if (dim == 0)
      return pts[idx][0];
    else if (dim == 1)
      return pts[idx][1];
    else {
      if (pdim == 3)
        return pts[idx][2];
      else
        return 0.0;
    }
  }

  // Optional bounding-box computation: return false to default to a standard
  // bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned
  //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
  //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX> bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
};

void generate_point_cloud(PointCloud<3> &pc, std::vector<Triangle<3>> &facets) {
  for (auto &facet : facets) {
    pc.pts.push_back(facet.centroid());
  }
}

void generate_point_cloud(PointCloud<2> &pc, std::vector<Line<2>> &facets) {
  for (auto &facet : facets) {
    pc.pts.push_back(facet.centroid());
  }
}

void facet_based_reinitialization_3D(
    double *x, Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, int num_total_nodes, double time) {
  using goma::distance_tools::Triangle;

  // We are going to loop over all elements and check if the level set interface exists on that
  // element.
  BasicTimer timer;
  timer.start();
  std::vector<Triangle<3>> facets;

  std::vector<std::array<Point<3>, 8>> hexes;
  std::vector<std::array<int, 8>> hex_nodes;

  std::unordered_set<int> level_set_nodes;
  for (int elem_block = 0; elem_block < exo->num_elem_blocks; elem_block++) {
    int mn = Matilda[elem_block];
    if (!(pd_glob[mn]->v[pg->imtrx][ls->var]))
      continue;
    int etype = exo->eb_elem_itype[elem_block];
    hexes.clear();
    hex_nodes.clear();
    int ielem_shape = type2shape(etype);
    int interpolation = pd_glob[mn]->i[pg->imtrx][ls->var];
    int dofs = getdofs(ielem_shape, interpolation);
    int shape = type2shape(etype);
    bool add_extra_dofs = false;
    switch (shape) {
    case HEXAHEDRON:
      shape_to_hex(etype, hexes, hex_nodes, dofs);
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Unsupported element type %d", etype);
      break;
    }

    for (int elem = exo->eb_ptr[elem_block]; elem < exo->eb_ptr[elem_block + 1]; elem++) {
      // If the level set interface exists on that element we will compute facets for that element.
      int iconnect = exo->elem_ptr[elem];
      int ls_dofs = add_extra_dofs ? 9 : dofs;
      std::vector<double> ls_values(ls_dofs);
      for (int ln = 0; ln < dofs; ln++) {
        int gnn = exo->node_list[iconnect + ln];
        ls_values[ln] = x[Index_Solution(gnn, ls->var, 0, 0, -2, pg->imtrx)];
        // if (std::abs(ls_values[ln]) < 1e-3) {
        //   level_set_nodes.insert(gnn);
        // }
      }

      double first_val = std::copysign(1.0, ls_values[0]);
      bool mark_elem = false;
      for (int ln = 1; ln < dofs; ln++) {
        if (std::abs(std::copysign(1.0, ls_values[ln]) + first_val) < 1e-15) {
          mark_elem = true;
          break;
        }
      }

      if (mark_elem) {
        std::vector<Triangle<3>> local_facets;
        for (size_t i = 0; i < hexes.size(); i++) {
          std::array<double, 8> values;
          for (int j = 0; j < 8; j++) {
            values[j] = ls_values[hex_nodes[i][j]];
          }
          auto facet_list = create_facet_from_hex(hexes[i], values, 0.0);

          for (auto &facet : facet_list) {
            local_facets.push_back(facet);
          }
        }

        // transform local facets to global facets
        if (!local_facets.empty()) {
          int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
          load_ei(elem, exo, NULL, pg->imtrx);
          if (DeformingMesh) {
            goma_error err = load_elem_dofptr(elem, exo, x, x, x, x, 0);
            GOMA_EH(err, "load_elem_dofptr");
          }
          int ShapeVar = pd->ShapeVar;
          int mdof = ei[pd->mi[ShapeVar]]->dof[ShapeVar];

          dbl x[DIM] = {0.0, 0.0, 0.0};
          dbl xi[DIM] = {0.0, 0.0, 0.0};
          for (size_t i = 0; i < local_facets.size(); i++) {
            Point<3> p0 = local_facets[i].p0;
            Point<3> p1 = local_facets[i].p1;
            Point<3> p2 = local_facets[i].p2;
            for (int p = 0; p < 3; p++) {
              if (p == 0) {
                xi[0] = p0[0];
                xi[1] = p0[1];
                xi[2] = p0[2];
              } else if (p == 1) {
                xi[0] = p1[0];
                xi[1] = p1[1];
                xi[2] = p1[2];
              } else {
                xi[0] = p2[0];
                xi[1] = p2[1];
                xi[2] = p2[2];
              }

              for (int a = 0; a < 3; a++) {
                x[a] = 0.0;

                if (!DeformingMesh) {
                  for (int j = 0; j < mdof; j++) {
                    int ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

                    int I = exo->node_list[iconnect + ln];

                    level_set_nodes.insert(I);

                    dbl phi_j =
                        newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln, ei[pg->imtrx]->ielem_shape,
                                 pd->i[pd->mi[ShapeVar]][ShapeVar], j);

                    x[a] += Coor[a][I] * phi_j;
                  }
                } else {
                  for (int j = 0; j < mdof; j++) {
                    int ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

                    int I = exo->node_list[iconnect + ln];

                    level_set_nodes.insert(I);

                    dbl phi_j =
                        newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln, ei[pg->imtrx]->ielem_shape,
                                 pd->i[pd->mi[ShapeVar]][ShapeVar], j);

                    x[a] += (Coor[a][I] + *esp->d[a][j]) * phi_j;
                  }
                }
              }

              if (p == 0) {
                p0[0] = x[0];
                p0[1] = x[1];
                p0[2] = x[2];
              } else if (p == 1) {
                p1[0] = x[0];
                p1[1] = x[1];
                p1[2] = x[2];
              } else {
                p2[0] = x[0];
                p2[1] = x[1];
                p2[2] = x[2];
              }
            }
            facets.push_back(Triangle<3>(p0, p1, p2));
          }
        }
      }
    }
  }

  timer.print_elapsed_and_reset("        Facet computation time:");

  // If the level set interface exists on that element we will compute facets for that element.

  // We will share our facets with other processors
  std::vector<double> facet_data(facets.size() * 9);
  for (size_t i = 0; i < facets.size(); i++) {
    facet_data[i * 9] = facets[i].p0[0];
    facet_data[i * 9 + 1] = facets[i].p0[1];
    facet_data[i * 9 + 2] = facets[i].p0[2];
    facet_data[i * 9 + 3] = facets[i].p1[0];
    facet_data[i * 9 + 4] = facets[i].p1[1];
    facet_data[i * 9 + 5] = facets[i].p1[2];
    facet_data[i * 9 + 6] = facets[i].p2[0];
    facet_data[i * 9 + 7] = facets[i].p2[1];
    facet_data[i * 9 + 8] = facets[i].p2[2];
  }

  std::vector<int> facet_sizes(Num_Proc, 0);
  int my_size = facet_data.size();
  MPI_Allgather(&my_size, 1, MPI_INT, facet_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<double> all_facets(std::accumulate(facet_sizes.begin(), facet_sizes.end(), 0));
  std::vector<int> offsets(Num_Proc, 0);
  int offset = 0;
  for (int i = 1; i < Num_Proc; i++) {
    offset += facet_sizes[i - 1];
    offsets[i] = offset;
  }

  MPI_Allgatherv(facet_data.data(), my_size, MPI_DOUBLE, all_facets.data(), facet_sizes.data(),
                 offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  facets.clear();
  facets.reserve(all_facets.size() / 9);
  for (size_t i = 0; i < all_facets.size() / 9; i++) {
    facets.push_back(Triangle<3>(
        Point<3>({all_facets[i * 9], all_facets[i * 9 + 1], all_facets[i * 9 + 2]}),
        Point<3>({all_facets[i * 9 + 3], all_facets[i * 9 + 4], all_facets[i * 9 + 5]}),
        Point<3>({all_facets[i * 9 + 6], all_facets[i * 9 + 7], all_facets[i * 9 + 8]})));
  }

  timer.print_elapsed_and_reset("       Facet communication time:");

  PointCloud<3> pc;
  generate_point_cloud(pc, facets);

  using my_kd_tree_t =
      nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<3>>,
                                          PointCloud<3>, 3 /* dim */
                                          >;

  my_kd_tree_t index(3, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));

  for (int node = 0; node < num_total_nodes; node++) {
    if (level_set_nodes.find(node) != level_set_nodes.end()) {
      continue;
    }
    int index_ls = Index_Solution(node, ls->var, 0, 0, -2, pg->imtrx);
    if (index_ls != -1) {
      Point<3> p({Coor[0][node], Coor[1][node], Coor[2][node]});
      if (pd->gv[MESH_DISPLACEMENT1]) {
        int index_dx = Index_Solution(node, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx);
        if (index_dx != -1) {
          p[0] += x[index_dx];
        }
        int index_dy = Index_Solution(node, MESH_DISPLACEMENT2, 0, 0, -2, pg->imtrx);
        if (index_dx != -1) {
          p[1] += x[index_dy];
        }
        int index_dz = Index_Solution(node, MESH_DISPLACEMENT3, 0, 0, -2, pg->imtrx);
        if (index_dz != -1) {
          p[2] += x[index_dz];
        }
      }

      const int num_results = 5;
      std::vector<size_t> ret_index(num_results);
      std::vector<double> out_dist_sqr(num_results);
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_index[0], &out_dist_sqr[0]);

      double query_pt[3] = {p[0], p[1], p[2]};
      index.findNeighbors(resultSet, query_pt, nanoflann::SearchParams());

      double min_distance = std::numeric_limits<double>::max();
      // for (size_t i = 0; i < facets.size(); i++) {
      for (size_t i = 0; i < resultSet.size(); i++) {
        int index = ret_index[i];
        double distance = facets[index].distance(p);
        min_distance = std::min(min_distance, distance);
      }

      // std::cout << "Node " << node << p[0] << " " << p[1] << " " << " distance: " << min_distance
      // << "\n";

      x[index_ls] = std::copysign(min_distance, x[index_ls]);
    }
  }

  timer.print_elapsed("        Reinitialization time:");
}

void facet_based_reinitialization_2D(
    double *x, Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, int num_total_nodes, double time) {
  using goma::distance_tools::Line;

  // We are going to loop over all elements and check if the level set interface exists on that
  // element.
  std::vector<Line<2>> facets;

  std::vector<std::tuple<Point<2>, Point<2>, Point<2>>> triangles;
  std::vector<std::tuple<int, int, int>> triangle_nodes;
  std::vector<std::tuple<Point<2>, Point<2>, Point<2>, Point<2>>> quads;
  std::vector<std::tuple<int, int, int, int>> quad_nodes;

  BasicTimer timer;
  timer.start();

  for (int elem_block = 0; elem_block < exo->num_elem_blocks; elem_block++) {
    int mn = Matilda[elem_block];
    if (!(pd_glob[mn]->v[pg->imtrx][ls->var]))
      continue;
    int etype = exo->eb_elem_itype[elem_block];
    triangles.clear();
    triangle_nodes.clear();
    quads.clear();
    quad_nodes.clear();
    int ielem_shape = type2shape(etype);
    int interpolation = pd_glob[mn]->i[pg->imtrx][ls->var];
    int dofs = getdofs(ielem_shape, interpolation);
    int shape = type2shape(etype);
    bool add_extra_dofs = false;
    switch (shape) {
    case TRIANGLE:
      shape_to_triangle(etype, triangles, triangle_nodes, dofs);
      break;
    case QUADRILATERAL:
      if (add_extra_dofs && dofs == 4) {
        shape_to_quad(etype, quads, quad_nodes, 9);
      } else {
        shape_to_quad(etype, quads, quad_nodes, dofs);
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Unsupported element type %d", etype);
      break;
    }

    for (int elem = exo->eb_ptr[elem_block]; elem < exo->eb_ptr[elem_block + 1]; elem++) {
      if (dpi->elem_owner[elem] != ProcID) {
        continue;
      }
      // If the level set interface exists on that element we will compute facets for that element.
      int iconnect = exo->elem_ptr[elem];
      int ls_dofs = add_extra_dofs ? 9 : dofs;
      std::vector<double> ls_values(ls_dofs);
      for (int ln = 0; ln < dofs; ln++) {
        int gnn = exo->node_list[iconnect + ln];
        ls_values[ln] = x[Index_Solution(gnn, ls->var, 0, 0, -2, pg->imtrx)];
      }
      if (add_extra_dofs && shape == QUADRILATERAL) {
        interp_quad_dofs(ls_values);
      }
      std::vector<Line<2>> local_facets;
      for (size_t i = 0; i < triangles.size(); i++) {
        auto [p1, p2, p3] = triangles[i];
        auto [v1, v2, v3] = triangle_nodes[i];
        std::array<double, 3> values = {ls_values[v1], ls_values[v2], ls_values[v3]};
        auto facet = create_facet_from_triangle<2>(p1, p2, p3, values, 0.0);

        if (facet.has_value()) {
          local_facets.push_back(facet.value());
        }
      }

      for (size_t i = 0; i < quads.size(); i++) {
        auto [p1, p2, p3, p4] = quads[i];
        auto [v1, v2, v3, v4] = quad_nodes[i];
        std::array<double, 4> values = {ls_values[v1], ls_values[v2], ls_values[v3], ls_values[v4]};
        auto facet_list = create_facet_from_quad<2>(p1, p2, p3, p4, values, 0.0);

        for (auto &facet : facet_list) {
          local_facets.push_back(facet);
        }
      }

      // transform local facets to global facets
      if (!local_facets.empty()) {
        int DeformingMesh = pd->e[pg->imtrx][R_MESH1];
        load_ei(elem, exo, NULL, pg->imtrx);
        if (DeformingMesh) {
          goma_error err = load_elem_dofptr(elem, exo, x, x, x, x, 0);
          GOMA_EH(err, "load_elem_dofptr");
        }
        int ShapeVar = pd->ShapeVar;
        int mdof = ei[pd->mi[ShapeVar]]->dof[ShapeVar];

        dbl x[DIM] = {0.0, 0.0, 0.0};
        dbl xi[DIM] = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < local_facets.size(); i++) {
          Point<2> p0 = local_facets[i].p0;
          Point<2> p1 = local_facets[i].p1;
          for (int p = 0; p < 2; p++) {
            if (p == 0) {
              xi[0] = p0[0];
              xi[1] = p0[1];
            } else {
              xi[0] = p1[0];
              xi[1] = p1[1];
            }

            for (int a = 0; a < 2; a++) {
              x[a] = 0.0;

              if (!DeformingMesh) {
                for (int j = 0; j < mdof; j++) {
                  int ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

                  int I = exo->node_list[iconnect + ln];

                  dbl phi_j =
                      newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln, ei[pg->imtrx]->ielem_shape,
                               pd->i[pd->mi[ShapeVar]][ShapeVar], j);

                  x[a] += Coor[a][I] * phi_j;
                }
              } else {
                for (int j = 0; j < mdof; j++) {
                  int ln = ei[pd->mi[ShapeVar]]->dof_list[ShapeVar][j];

                  int I = exo->node_list[iconnect + ln];

                  dbl phi_j =
                      newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ln, ei[pg->imtrx]->ielem_shape,
                               pd->i[pd->mi[ShapeVar]][ShapeVar], j);

                  x[a] += (Coor[a][I] + *esp->d[a][j]) * phi_j;
                }
              }
            }

            if (p == 0) {
              p0[0] = x[0];
              p0[1] = x[1];
            } else {
              p1[0] = x[0];
              p1[1] = x[1];
            }
          }
          facets.push_back(Line<2>(p0, p1));
        }
      }
    }
  }
  timer.print_elapsed_and_reset("Create facets");

  std::ofstream file("facets.txt");
  file << "x,y\n";
  for (size_t i = 0; i < facets.size(); i++) {
    file << facets[i].p0[0] << "," << facets[i].p0[1] << "\n"
         << facets[i].p1[0] << "," << facets[i].p1[1] << "\n";
  }
  file.close();

  // std::exit(0);

  // We will share our facets with other processors
  // pack facets into a vector
  std::vector<double> facet_data(facets.size() * 4);
  for (size_t i = 0; i < facets.size(); i++) {
    facet_data[i * 4] = facets[i].p0[0];
    facet_data[i * 4 + 1] = facets[i].p0[1];
    facet_data[i * 4 + 2] = facets[i].p1[0];
    facet_data[i * 4 + 3] = facets[i].p1[1];
  }

  std::vector<int> facet_sizes(Num_Proc, 0);
  int my_size = facets.size() * 4;
  MPI_Allgather(&my_size, 1, MPI_INT, facet_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<double> all_facets(std::accumulate(facet_sizes.begin(), facet_sizes.end(), 0));
  std::vector<int> offsets(Num_Proc, 0);
  int offset = 0;
  for (int i = 1; i < Num_Proc; i++) {
    offset += facet_sizes[i - 1];
    offsets[i] = offset;
  }

  MPI_Allgatherv(facet_data.data(), my_size, MPI_DOUBLE, all_facets.data(), facet_sizes.data(),
                 offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  facets.clear();
  facets.reserve(all_facets.size() / 4);
  for (size_t i = 0; i < all_facets.size() / 4; i++) {
    facets.push_back(Line<2>({all_facets[i * 4], all_facets[i * 4 + 1]},
                             {all_facets[i * 4 + 2], all_facets[i * 4 + 3]}));
  }

  timer.print_elapsed_and_reset("Facet MPI communication");
  PointCloud<2> pc;
  generate_point_cloud(pc, facets);

  using my_kd_tree_t =
      nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<2>>,
                                          PointCloud<2>, 2 /* dim */
                                          >;

  my_kd_tree_t index(2, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10));

  // For each node we will compute the distance to the closest facet
  for (int node = 0; node < num_total_nodes; node++) {
    int index_ls = Index_Solution(node, ls->var, 0, 0, -2, pg->imtrx);
    if (index_ls != -1) {
      Point<2> p = {Coor[0][node], Coor[1][node]};
      if (pd->gv[MESH_DISPLACEMENT1]) {
        int index_dx = Index_Solution(node, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx);
        if (index_dx != -1) {
          p[0] += x[index_dx];
        }
        int index_dy = Index_Solution(node, MESH_DISPLACEMENT2, 0, 0, -2, pg->imtrx);
        if (index_dx != -1) {
          p[1] += x[index_dy];
        }
      }

      // std::cout << "Node " << node << p[0] << " " << p[1] << " " << " distance: " << min_distance
      // << "\n";
      const int num_results = 5;
      std::vector<size_t> ret_index(num_results);
      std::vector<double> out_dist_sqr(num_results);
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_index[0], &out_dist_sqr[0]);

      double query_pt[2] = {p[0], p[1]};
      index.findNeighbors(resultSet, query_pt, nanoflann::SearchParams());

      double min_distance = std::numeric_limits<double>::max();
      // for (size_t i = 0; i < facets.size(); i++) {
      for (size_t i = 0; i < resultSet.size(); i++) {
        int index = ret_index[i];
        double distance = facets[index].distance(p);
        min_distance = std::min(min_distance, distance);
      }

      x[index_ls] = std::copysign(min_distance, x[index_ls]);
    }
  }

  timer.print_elapsed_and_reset("Reinitialization using facets");
}

extern "C" void facet_based_reinitialization(
    double *x, Exo_DB *exo, Comm_Ex *cx, Dpi *dpi, int num_total_nodes, double time) {
  if (ProcID == 0)
    std::cout << "Facet based reinitialization\n";
  if (exo->num_dim == 2) {
    facet_based_reinitialization_2D(x, exo, cx, dpi, num_total_nodes, time);
  } else if (exo->num_dim == 3) {
    facet_based_reinitialization_3D(x, exo, cx, dpi, num_total_nodes, time);
  } else {
    GOMA_EH(GOMA_ERROR, "Unsupported dimension %d", exo->num_dim);
  }
}