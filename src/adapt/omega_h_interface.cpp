#ifdef HAVE_OMEGA_H
#include "adapt/omega_h_interface.h"

#include <Omega_h_adapt.hpp>
#include <Omega_h_adj.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_config.h>
#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_macros.h>
#include <Omega_h_matrix.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_profile.hpp>
#include <Omega_h_remotes.hpp>
#include <Omega_h_vector.hpp>
#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <set>
#include <string.h>
#include <string>
#include <utility>
#include <vector>

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"

extern "C" {
#define DISABLE_CPP
#ifndef MAX_PDIM
#define MAX_PDIM 3 /* Maximum physical problem dimension    */
#endif
#include <mm_bc.h>

#include "dp_types.h"
#include "dpi.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_interface.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "std.h"
#include <exodusII.h>
#define IGNORE_CPP_DEFINE
#include "sl_util_structs.h"
#undef IGNORE_CPP_DEFINE
#include "adapt/resetup_problem.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_conn.h"
#include "mm_unknown_map.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "util/goma_normal.h"
#include "wr_dpi.h"
#include "wr_exo.h"

extern int ***Local_Offset;
extern int ***Dolphin;
extern int *NumUnknowns;    /* Number of unknown variables updated by this   */
extern int *NumExtUnknowns; /* Number of unknown variables updated by this   */
extern int ***idv;
extern int **local_ROT_list;
extern int rotation_allocated;
extern Comm_Ex **cx;
#include "bc/rotate_coordinates.h"

#undef DISABLE_CPP
}

//#define DEBUG_OMEGA_H

namespace Omega_h {

#define CALL(f)                                                                          \
  do {                                                                                   \
    auto f_err = (f);                                                                    \
    if (f_err != 0) {                                                                    \
      const char *errmsg;                                                                \
      const char *errfunc;                                                               \
      int errnum;                                                                        \
      ex_get_err(&errmsg, &errfunc, &errnum);                                            \
      Omega_h_fail("Exodus call %s failed (%d): %s: %s\n", #f, errnum, errfunc, errmsg); \
    }                                                                                    \
  } while (0)

namespace exodus {
#ifndef OMEGA_H_USE_SEACASEXODUS
enum ClassifyWith {
  NODE_SETS = 0x1,
  SIDE_SETS = 0x2,
};
#endif
static OMEGA_H_INLINE int side_osh2exo(int dim, int side) {
  switch (dim) {
  case 2:
    switch (side) {
    case 1:
      return 2;
    case 2:
      return 3;
    case 0:
      return 1;
    }
    return -1;
  case 3:
    switch (side) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    case 0:
      return 4;
    }
    return -1;
  }
  return -1;
}

static void setup_names(int nnames, std::vector<char> &storage, std::vector<char *> &ptrs) {
  constexpr auto max_name_length = MAX_STR_LENGTH + 1;
  storage = std::vector<char>(std::size_t(nnames * max_name_length), '\0');
  ptrs = std::vector<char *>(std::size_t(nnames), nullptr);
  for (int i = 0; i < nnames; ++i) {
    ptrs[std::size_t(i)] = storage.data() + max_name_length * i;
  }
}

static OMEGA_H_INLINE int side_exo2osh(Omega_h_Family family, int dim, int side) {
  switch (family) {
  case OMEGA_H_SIMPLEX:
    switch (dim) {
    case 2:
      return (side + 2) % 3;
    case 3:
      switch (side) {
      case 1:
        return 1;
      case 2:
        return 2;
      case 3:
        return 3;
      case 4:
        return 0;
      }
    }
    return -1;
  case OMEGA_H_HYPERCUBE:
    return -1;
  }
  return -1;
}

void convert_goma_to_omega_h(Exo_DB *exo, Dpi *dpi, double **x, Mesh *mesh, bool verbose) {
  std::vector<LO> local_to_global(exo->num_nodes);
  for (int i = 0; i < exo->num_nodes; i++) {
    local_to_global[i] = dpi->node_index_global[i];
  }

  std::vector<LO> owned_elems;
  for (int i = 0; i < exo->num_elems; i++) {
    if (dpi->elem_owner[i] == ProcID) {
      owned_elems.push_back(i);
    }
  }
  int nnodes_per_elem = exo->eb_num_nodes_per_elem[0];

  std::set<LO> local_verts;
  for (auto e : owned_elems) {
    for (int i = 0; i < nnodes_per_elem; i++) {
      local_verts.insert(exo->eb_conn[0][e * nnodes_per_elem + i]);
    }
  }

  std::vector<LO> old_locals(local_verts.begin(), local_verts.end());
  std::sort(old_locals.begin(), old_locals.end());
  std::vector<LO> new_local_verts(old_locals.size());
  std::iota(new_local_verts.begin(), new_local_verts.end(), 0);

  std::map<LO, LO> old_to_new;
  std::map<LO, LO> new_to_old;
  for (size_t i = 0; i < old_locals.size(); i++) {
    old_to_new.insert(std::pair<LO, LO>(old_locals[i], i));
    new_to_old.insert(std::pair<LO, LO>(i, old_locals[i]));
  }

  HostWrite<GO> vert_global;
  vert_global = decltype(vert_global)(GO(new_local_verts.size()), "global vertices");

  for (size_t i = 0; i < old_locals.size(); i++) {
    vert_global[i] = dpi->node_index_global[old_locals[i]];
  }

  HostWrite<LO> h_conn;
  h_conn = decltype(h_conn)(LO(owned_elems.size() * nnodes_per_elem), "host connectivity");

  for (size_t i = 0; i < owned_elems.size(); i++) {
    auto e = owned_elems[i];
    for (int j = 0; j < nnodes_per_elem; j++) {
      auto old_index = exo->eb_conn[0][e * nnodes_per_elem + j];
      auto new_index = old_to_new[old_index];
      h_conn[i * nnodes_per_elem + j] = new_index;
    }
  }

  auto dim = exo->num_dim;

  mesh->set_parting(OMEGA_H_ELEM_BASED);

  auto conn = LOs(h_conn.write());
  build_from_elems2verts(mesh, mesh->library()->world(), OMEGA_H_SIMPLEX, dim, conn,
                         GOs(vert_global));

  auto new_verts = mesh->globals(0);

  HostWrite<Real> h_coords(LO(new_verts.size() * dim));
  std::vector<LO> exo_from_omega;

  for (int i = 0; i < new_verts.size(); i++) {
    int idx = in_list(new_verts[i], 0, exo->num_nodes, dpi->node_index_global);
    assert(idx != -1);
    exo_from_omega.push_back(idx);

    double d_x = 0;
    double d_y = 0;
    double d_z = 0;

    // Get displacements to Anneal the mesh
    if (pd->gv[R_MESH1]) {
      int x_idx = Index_Solution(idx, R_MESH1, 0, 0, -1, pd->mi[R_MESH1]);
      int y_idx = Index_Solution(idx, R_MESH2, 0, 0, -1, pd->mi[R_MESH2]);
      int z_idx = -1;
      if (exo->num_dim == 3) {
        z_idx = Index_Solution(idx, R_MESH3, 0, 0, -1, pd->mi[R_MESH3]);
      }

      d_x = x[pd->mi[R_MESH1]][x_idx];
      d_y = x[pd->mi[R_MESH2]][y_idx];
      if (exo->num_dim == 3) {
        d_z = x[pd->mi[R_MESH3]][z_idx];
      }
    }

    h_coords[i * dim + 0] = exo->x_coord[idx] + d_x;
    h_coords[i * dim + 1] = exo->y_coord[idx] + d_y;
    if (exo->num_dim == 3) {
      h_coords[i * dim + 2] = exo->z_coord[idx] + d_z;
    }
  }

  auto coords = Reals(h_coords.write());
  mesh->add_coords(coords);
  Write<LO> elem_class_ids_w(LO(owned_elems.size()));
  for (size_t i = 0; i < owned_elems.size(); i++) {
    elem_class_ids_w[i] = 1;
  }
  std::map<LO, LO> exo_to_global;
  for (int node = 0; node < mesh->globals(0).size(); node++) {
    int exo_index = in_list(mesh->globals(0)[node], 0, exo->num_nodes, dpi->node_index_global);
    exo_to_global.insert(std::pair<LO, LO>(exo_index, node));
  }
  for (int j = V_FIRST; j < V_LAST; j++) {
    int imtrx = upd->matrix_index[j];

    if (imtrx >= 0) {
      if (j == MASS_FRACTION) {
        for (int mf = 0; mf < upd->Max_Num_Species; mf++) {
          auto var_values = Omega_h::Write<Omega_h::Real>(mesh->nverts());
          for (int i = 0; i < exo->num_nodes; i++) {
            if (exo_to_global.find(i) != exo_to_global.end()) {
              auto gnode = exo_to_global[i];
              int ja = Index_Solution(i, j, mf, 0, -2, imtrx);
              GOMA_EH(ja, "could not find solution");
              var_values[gnode] = x[imtrx][ja];
            }
          }
          std::string species_name = Exo_Var_Names[j].name2 + std::to_string(mf);
          mesh->add_tag(Omega_h::VERT, species_name, 1, Omega_h::Reals(var_values));
        }
      } else {
        auto var_values = Omega_h::Write<Omega_h::Real>(mesh->nverts());
        for (int i = 0; i < exo->num_nodes; i++) {
          if (exo_to_global.find(i) != exo_to_global.end()) {
            auto gnode = exo_to_global[i];
            int ja = Index_Solution(i, j, 0, 0, -2, imtrx);
            GOMA_EH(ja, "could not find solution");
            var_values[gnode] = x[imtrx][ja];
            if (tran->ale_adapt && (j >= R_MESH1 && j <= R_MESH3)) {
              var_values[gnode] = 0;
            }
          }
        }
        mesh->add_tag(Omega_h::VERT, Exo_Var_Names[j].name2, 1, Omega_h::Reals(var_values));

        if (tran->ale_adapt && j == R_MESH1 && ((ls == NULL) || (!ls->adapt))) {
          auto target_metrics =
              Omega_h::Write<Omega_h::Real>(mesh->nverts() * Omega_h::symm_ncomps(mesh->dim()));
          auto f0 = OMEGA_H_LAMBDA(Omega_h::LO index) {
            auto iso_size = tran->ale_adapt_iso_size;
            if (mesh->dim() == 2) {
              auto target_metric = Omega_h::compose_metric(Omega_h::identity_matrix<2, 2>(),
                                                           Omega_h::vector_2(iso_size, iso_size));
              Omega_h::set_vector(target_metrics, index, Omega_h::symm2vector(target_metric));
            } else {
              auto target_metric =
                  Omega_h::compose_metric(Omega_h::identity_matrix<3, 3>(),
                                          Omega_h::vector_3(iso_size, iso_size, iso_size));
              Omega_h::set_vector(target_metrics, index, Omega_h::symm2vector(target_metric));
            }
          };

          Omega_h::parallel_for(mesh->nverts(), f0, "set_iso_metric_values");
          mesh->add_tag(Omega_h::VERT, "iso_size_metric", Omega_h::symm_ncomps(mesh->dim()),
                        Omega_h::Reals(target_metrics));
        } else if (ls != NULL && ls->adapt && (j == FILL)) {
          auto target_metrics =
              Omega_h::Write<Omega_h::Real>(mesh->nverts() * Omega_h::symm_ncomps(mesh->dim()));
          auto f0 = OMEGA_H_LAMBDA(Omega_h::LO index) {
            auto F = var_values[index];
            auto iso_size = ls->adapt_outer_size;
            if (std::abs(F) < ls->adapt_width) {
              iso_size = ls->adapt_inner_size;
            }
            if (mesh->dim() == 2) {
              auto target_metric = Omega_h::compose_metric(Omega_h::identity_matrix<2, 2>(),
                                                           Omega_h::vector_2(iso_size, iso_size));
              Omega_h::set_vector(target_metrics, index, Omega_h::symm2vector(target_metric));
            } else {
              auto target_metric =
                  Omega_h::compose_metric(Omega_h::identity_matrix<3, 3>(),
                                          Omega_h::vector_3(iso_size, iso_size, iso_size));
              Omega_h::set_vector(target_metrics, index, Omega_h::symm2vector(target_metric));
            }
          };

          Omega_h::parallel_for(mesh->nverts(), f0, "set_iso_metric_values");
          mesh->add_tag(Omega_h::VERT, "iso_size_metric", Omega_h::symm_ncomps(mesh->dim()),
                        Omega_h::Reals(target_metrics));
        }
      }
    }
  }

  for (int w = 0; w < efv->Num_external_field; w++) {
    auto var_values = Omega_h::Write<Omega_h::Real>(mesh->nverts());
    for (int i = 0; i < exo->num_nodes; i++) {
      if (exo_to_global.find(i) != exo_to_global.end()) {
        auto gnode = exo_to_global[i];
        var_values[gnode] = efv->ext_fld_ndl_val[w][i];
      }
    }
    mesh->add_tag(Omega_h::VERT, efv->name[w], 1, Omega_h::Reals(var_values));
  }

  std::vector<int> side_set_ids(std::size_t(dpi->num_side_sets_global));
  for (int i = 0; i < dpi->num_side_sets_global; i++) {
    side_set_ids[i] = dpi->ss_id_global[i];
  }
  //  CALL(ex_get_ids(file, EX_SIDE_SET, side_set_ids.data()));
  Write<LO> side_class_ids_w(mesh->nents(dim - 1), -1);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  Write<I8> side_class_dims_w = deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_value_into(0, exposed_sides2side, side_class_ids_w);
#if 0
  if (dpi->num_side_sets_global) {
    int max_side_set_id = 0;
    if (side_set_ids.size()) {
      max_side_set_id = *std::max_element(side_set_ids.begin(), side_set_ids.end());
    }
    std::vector<int> node_set_ids(std::size_t(exo->num_node_sets));
    //    CALL(ex_get_ids(file, EX_NODE_SET, node_set_ids.data()));
    for (int i = 0; i < exo->num_node_sets; i++) {
      node_set_ids[i] = exo->ns_id[i];
    }
    std::vector<char> names_memory;
    std::vector<char *> name_ptrs;
    setup_names(int(dpi->num_node_sets_global + dpi->num_side_sets_global), names_memory, name_ptrs);
    //    CALL(ex_get_names(file, EX_NODE_SET, name_ptrs.data()));

    for (size_t i = 0; i < exo->num_node_sets; ++i) {
      int global_offset;
      for (global_offset = 0; global_offset < dpi->num_node_sets_global; global_offset++) {
        if (exo->ns_id[i] == dpi->ns_id_global[global_offset]) break;
      }
      assert(global_offset < dpi->num_node_sets_global);
      int nentries;
      //      CALL(ex_get_set_param(file, EX_NODE_SET, node_set_ids[i], &nentries, &ndist_factors));
      nentries = exo->ns_num_nodes[i];
      if (verbose) {
        std::cout << "node set " << node_set_ids[i] << " has " << nentries << " nodes\n";
      }

      int nnodes_owned = 0;
      for (int index = 0; index < nentries; index++) {
        int exoindex = exo->ns_node_list[exo->ns_node_index[i] + index];
        if (exo_to_global.find(exoindex) != exo_to_global.end()) {
          nnodes_owned++;
        }
      }

      HostWrite<LO> h_set_nodes2nodes(nnodes_owned);
      int index = 0;
      for (int idx = 0; idx < nentries; idx++) {
        int exoindex = exo->ns_node_list[exo->ns_node_index[i] + index];
        if (exo_to_global.find(exoindex) != exo_to_global.end()) {
          h_set_nodes2nodes[index] = exo_to_global.at(exoindex);
          index++;
        }
      }
      //      CALL(ex_get_set(file, EX_NODE_SET, node_set_ids[i], h_set_nodes2nodes.data(),
      //      nullptr));
      auto set_nodes2nodes = LOs(h_set_nodes2nodes.write());
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set = mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side = collect_marked(sides_are_in_set);
      auto surface_id = dpi->ns_id_global[global_offset] + dpi->num_node_sets_global;
      if (verbose) {
        std::cout << "node set #" << node_set_ids[i] << " \"" << name_ptrs[i]
                  << "\" will be surface " << surface_id << '\n';
      }
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
    }
    for (int offset = 0; offset < dpi->num_node_sets_global; offset++) {
      auto surface_id = dpi->ns_id_global[offset] + dpi->num_side_sets_global;
      mesh->class_sets[name_ptrs[offset]].push_back({I8(dim - 1), surface_id});
    }
  }
#endif
  if (1 && dpi->num_side_sets_global) {
    std::vector<char> names_memory;
    std::vector<char *> name_ptrs;
    setup_names(int(dpi->num_side_sets_global), names_memory, name_ptrs);
    //    CALL(ex_get_names(file, EX_SIDE_SET, name_ptrs.data()));
    for (int i = 0; i < dpi->num_side_sets_global; ++i) {
      assert(exo->ss_id[i] == dpi->ss_id_global[i]);
      int nsides;
      //      CALL(ex_get_set_param(file, EX_SIDE_SET, side_set_ids[i], &nsides, &ndist_factors));
      nsides = exo->ss_num_sides[i];
      if (verbose)
        std::cout << "Proc" << ProcID << " has " << nsides << "sides\n";
      int snl[MAX_NODES_PER_SIDE]; /* Side Node List - NOT Saturday Night Live! */

      std::set<int> ss_side_nodes;
      std::set<int> ss_side_nodes_global;
      for (int es = 0; es < exo->ss_num_sides[i]; es++) {
        for (int q = exo->ss_elem_index[i]; q < (exo->ss_elem_index[i] + exo->ss_num_sides[i]);
             q++) {
          int elem = exo->ss_elem_list[q];
          int side = exo->ss_side_list[q];
          int num_nodes = build_side_node_list(elem, side - 1, exo, snl);
          for (int j = 0; j < num_nodes; j++) {
            ss_side_nodes.insert(snl[j]);
            if (exo_to_global.find(snl[j]) != exo_to_global.end()) {
              ss_side_nodes_global.insert(mesh->globals(0)[exo_to_global[snl[j]]]);
              if (verbose)
                std::cout << "node found " << mesh->globals(0)[exo_to_global[snl[j]]] << " "
                          << dpi->node_index_global[snl[j]] << "\n";
            } else {
              if (verbose)
                std::cout << "node not found " << dpi->node_index_global[snl[j]] << "\n";
            }
          }
        }
      }

      if (verbose) {
        std::cout << "ss#" << dpi->ss_id_global[i] << " Proc " << ProcID << " ss_side_nodes ";

        for (auto n : ss_side_nodes) {
          std::cout << n << " ";
        }
        std::cout << "\n";
      }

      int nnodes = 0;
      for (int side = 0; side < nsides; side++) {
        nnodes += exo->ss_node_cnt_list[i][side];
      }
      int nnodes_owned = 0;
      for (int side = 0; side < nnodes; side++) {
        //        int exoindex = exo->ss_node_list[i][side];
        //        if (exo_to_global.find(exoindex) != exo_to_global.end()) {
        nnodes_owned++;
        //        }
      }

      int max_nnodes;
      MPI_Allreduce(&nnodes_owned, &max_nnodes, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      std::vector<int> all_owned(max_nnodes * Num_Proc);
      std::vector<int> my_owned(max_nnodes);
      std::fill(my_owned.begin(), my_owned.end(), -1);

      HostWrite<LO> h_set_nodes_q(nnodes_owned);
      int index = 0;
      for (int idx = 0; idx < nnodes; idx++) {
        int exoindex = exo->ss_node_list[i][idx];
        //        if (exo_to_global.find(exoindex) != exo_to_global.end()) {
        h_set_nodes_q[index] = dpi->node_index_global[exoindex];
        index++;
        //        }
      }
      for (int i = 0; i < index; i++) {
        my_owned[i] = h_set_nodes_q[i];
      }
      std::set<int> my_owned_set(my_owned.begin(), my_owned.end());

      MPI_Allgather(my_owned.data(), max_nnodes, MPI_INT, all_owned.data(), max_nnodes, MPI_INT,
                    MPI_COMM_WORLD);
      std::set<int> owned_set(all_owned.begin(), all_owned.end());
      std::vector<int> my_nodes_sorted(mesh->nverts());
      for (int i = 0; i < mesh->nverts(); i++) {
        my_nodes_sorted[i] = mesh->globals(0)[i];
      }
      std::sort(my_nodes_sorted.begin(), my_nodes_sorted.end());
      std::vector<int> all_owned_set(owned_set.begin(), owned_set.end());
      std::vector<int> mark_nodes;
      mark_nodes.reserve(all_owned_set.size());
      for (size_t j = 0; j < all_owned_set.size(); j++) {
        if (all_owned_set[j] != -1 &&
            std::binary_search(my_nodes_sorted.begin(), my_nodes_sorted.end(), all_owned_set[j])) {
          // TODO
          for (int k = 0; k < mesh->nverts(); k++) {
            if (mesh->globals(0)[k] == all_owned_set[j]) {
              mark_nodes.push_back(k);
              break;
            }
          }
        }
      }

      HostWrite<LO> h_set_nodes(mark_nodes.size());
      for (size_t j = 0; j < mark_nodes.size(); j++) {
        h_set_nodes[j] = mark_nodes[j];
      }

      auto set_nodes2nodes = LOs(h_set_nodes.write());
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set = mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side_tmp = collect_marked(sides_are_in_set);
      int max_sides;
      int my_size = set_sides2side_tmp.size();
      MPI_Allreduce(&my_size, &max_sides, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      std::vector<int> send_sides(max_sides);
      std::fill(send_sides.begin(), send_sides.end(), -1);
      for (int j = 0; j < my_size; j++) {
        send_sides[j] = mesh->globals(dim - 1)[set_sides2side_tmp[j]];
      }
      std::vector<int> all_sides(max_sides * Num_Proc);
      MPI_Allgather(send_sides.data(), max_sides, MPI_INT, all_sides.data(), max_sides, MPI_INT,
                    MPI_COMM_WORLD);

      std::vector<int> sides;

      for (size_t j = 0; j < all_sides.size(); j++) {
        for (int k = 0; k < mesh->globals(dim - 1).size(); k++) {
          if (all_sides[j] == mesh->globals(dim - 1)[k]) {
            sides.push_back(k);
          }
        }
      }
      std::sort(sides.begin(), sides.end());

      if (verbose)
        std::cout << "Proc" << ProcID << " has found " << sides.size() << "sides\n";
      Write<LO> set_sides2side(sides.size());

      for (size_t j = 0; j < sides.size(); j++) {
        set_sides2side[j] = sides[j];
      }

      auto surface_id = dpi->ss_id_global[i];
      if (verbose) {
        std::cout << "P" << ProcID << " side set #" << surface_id << " \"" << name_ptrs[i]
                  << "\" has " << nsides << " sides, will be surface " << surface_id << "\n";

        std::cout << "P" << ProcID << " side set #" << surface_id << " nodes = ";
        for (int i = 0; i < h_set_nodes.size(); i++) {
          std::cout << mesh->globals(0).data()[h_set_nodes.data()[i]] << " ";
        }
        std::cout << "P" << ProcID << " side set #" << surface_id << " sides = ";
        for (int i = 0; i < set_sides2side.size(); i++) {
          std::cout << mesh->globals(dim - 1).data()[set_sides2side.data()[i]] << " ";
        }
        std::cout << "\n\n\n";
      }
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
      mesh->class_sets[name_ptrs[i]].push_back({I8(dim - 1), surface_id});
    }
    for (int offset = 0; offset < dpi->num_side_sets_global; offset++) {
      auto surface_id = dpi->ss_id_global[offset];
      mesh->class_sets[name_ptrs[offset]].push_back({I8(dim - 1), surface_id});
    }
  }
  auto elem_class_ids = LOs(elem_class_ids_w);
  auto side_class_ids = LOs(side_class_ids_w);
  auto side_class_dims = Read<I8>(side_class_dims_w);
  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  /*
  classify_elements(mesh);
  auto elem_class_ids = LOs(elem_class_ids_w);
  //auto side_class_ids = LOs(side_class_ids_w);
  //auto side_class_dims = Read<I8>(side_class_dims_w);
//  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  //mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  //mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
//  mesh->balance();
//  mesh->set_parting(OMEGA_H_GHOSTED);
   */
  classify_elements(mesh);
  //  auto elem_class_ids = LOs(elem_class_ids_w);
  // auto side_class_ids = LOs(side_class_ids_w);
  // auto side_class_dims = Read<I8>(side_class_dims_w);
  //  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  // mesh->set_parting(OMEGA_H_GHOSTED);
  //  auto sides_are_exposed = mark_exposed_sides(mesh);
  //  classify_sides_by_exposure(mesh, sides_are_exposed);
  // mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  // mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
}

void convert_omega_h_to_goma(
    const char *path, Mesh *mesh, Exo_DB *exo, Dpi *dpi, bool verbose, int classify_with) {

  //  Omega_h::exodus::write(std::to_string(ProcID) + "tmp.e", mesh, true, classify_with);
  char out_par[MAX_FNL];
  strncpy(out_par, "tmp_oh.e", MAX_FNL - 1);
  multiname(out_par, ProcID, Num_Proc);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  //  Omega_h::exodus::write(out_par, mesh, true, classify_with);
  // mesh->set_parting(OMEGA_H_GHOSTED);
  strncpy(out_par, "tmp.e", MAX_FNL - 1);
  multiname(out_par, ProcID, Num_Proc);
  auto comp_ws = int(sizeof(Real));
  auto io_ws = comp_ws;
  auto exoid = ex_create(out_par, EX_CLOBBER, &comp_ws, &io_ws);
  std::set<LO> region_set;
  auto dim = mesh->dim();
  auto title = "Omega_h " OMEGA_H_SEMVER " Exodus Output";

  auto class_sets = mesh->class_sets;
  // TODO multiblock
  std::set<LO> surface_set;
  std::set<LO> node_set;
  for (auto &it : class_sets) {
    auto key = it.first;
    auto value = it.second;
    for (size_t i = 0; i < value.size(); i++) {
      if (static_cast<int>(i) < dpi->num_side_sets_global) {
        surface_set.insert(value[i].id);
      } else {
        node_set.insert(value[i].id);
      }
    }
  }
  auto elem_class_ids = mesh->get_array<ClassId>(dim, "class_id");
  auto h_elem_class_ids = HostRead<LO>(elem_class_ids);
  for (LO i = 0; i < h_elem_class_ids.size(); ++i) {
    region_set.insert(h_elem_class_ids[i]);
  }
  auto side_class_ids = mesh->get_array<ClassId>(dim - 1, "class_id");
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  auto h_side_class_ids = HostRead<LO>(side_class_ids);
  auto h_side_class_dims = HostRead<I8>(side_class_dims);
  auto nelem_blocks = int(region_set.size());
  auto nside_sets = (classify_with & exodus::SIDE_SETS) ? int(surface_set.size()) : 0;
  auto nnode_sets = (classify_with & exodus::NODE_SETS) ? int(node_set.size()) : 0;
  //  Omega_h::binary::write("tmp.osh", mesh);

  auto all_conn = mesh->ask_elem_verts();
  auto deg = element_degree(mesh->family(), dim, VERT);
  auto vert_owners = mesh->ask_owners(0).ranks;
  std::set<int> neighbors_set;
  for (int i = 0; i < vert_owners.size(); i++) {
    if (vert_owners[i] != ProcID) {
      neighbors_set.insert(vert_owners[i]);
    }
  }

  int my_neighbors = neighbors_set.size();
  int max_neighbors;
  MPI_Allreduce(&my_neighbors, &max_neighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  std::vector<int> neighbors(neighbors_set.begin(), neighbors_set.end());
  std::vector<int> neighbor_recv(max_neighbors);
  std::fill(neighbor_recv.begin(), neighbor_recv.end(), 0);
  std::sort(neighbors.begin(), neighbors.end());
  for (int i = neighbors.size(); i < max_neighbors; i++) {
    neighbors.push_back(-1);
  }
  for (int i = 0; i < vert_owners.size(); i++) {
    int proc = vert_owners[i];
    if (proc != ProcID) {
      for (size_t j = 0; j < neighbors.size(); j++) {
        if (neighbors[j] == proc) {
          neighbor_recv[j] += 1;
        }
      }
    }
  }

  std::vector<int> allneighbors(max_neighbors * Num_Proc);
  std::vector<int> allneighrecv(max_neighbors * Num_Proc);

  MPI_Allgather(neighbors.data(), max_neighbors, MPI_INT, allneighbors.data(), max_neighbors,
                MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(neighbor_recv.data(), max_neighbors, MPI_INT, allneighrecv.data(), max_neighbors,
                MPI_INT, MPI_COMM_WORLD);

  std::vector<MPI_Request> requests;

  std::vector<std::vector<int>> neighbor_needed(my_neighbors);
  for (int i = 0; i < my_neighbors; i++) {
    int proc = neighbors[i];
    for (int j = 0; j < vert_owners.size(); j++) {
      if (vert_owners[j] == proc) {
        neighbor_needed[i].push_back(mesh->globals(0)[j]);
      }
    }
    assert(neighbor_needed[i].size() == static_cast<size_t>(neighbor_recv[i]));
    requests.emplace_back();
    if (verbose) {
      std::cout << "Proc " << ProcID << " send proc " << proc << " count "
                << neighbor_needed[i].size() << " with tag " << ProcID * 100 + proc << "\n";
    }
    MPI_Isend(neighbor_needed[i].data(), neighbor_needed[i].size(), MPI_INT, proc,
              ProcID * 100 + proc, MPI_COMM_WORLD, &requests[i]);
  }
  std::vector<std::vector<int>> neighbor_send;
  std::vector<int> neighbor_send_id;
  int nindex = 0;
  int req_index = requests.size();
  for (size_t i = 0; i < allneighbors.size(); i++) {
    if (allneighbors[i] == ProcID) {
      int proc = i / max_neighbors;
      neighbor_send.emplace_back();
      neighbor_send[nindex].resize(allneighrecv[i]);
      neighbor_send_id.emplace_back(proc);
      requests.emplace_back();
      if (verbose) {
        std::cout << "Proc " << ProcID << " recv proc " << proc << " count "
                  << neighbor_send[nindex].size() << " with tag " << proc * 100 + ProcID << "\n";
      }
      MPI_Irecv(neighbor_send[nindex].data(), neighbor_send[nindex].size(), MPI_INT, proc,
                proc * 100 + ProcID, MPI_COMM_WORLD, &requests[req_index]);
      nindex++;
      req_index++;
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  std::set<int> boundary_nodes;
  for (auto &bound_nodes : neighbor_send) {
    for (auto gnode : bound_nodes) {
#ifndef NDEBUG
      bool found = false;
#endif
      for (int k = 0; k < mesh->globals(0).size(); k++) {
        if (mesh->globals(0)[k] == gnode) {
          boundary_nodes.insert(k);
#ifndef NDEBUG
          found = true;
#endif
          break;
        }
      }
      assert(found);
    }
  }
  std::vector<bool> keep_nodes(mesh->nverts());
  std::fill(keep_nodes.begin(), keep_nodes.end(), false);

  std::vector<int> mesh_owned_elements;
  for (int i = 0; i < mesh->nelems(); i++) {
    bool fully_other = true;
    for (int j = 0; j < deg; j++) {
      auto local = all_conn[i * deg + j];
      if (vert_owners[local] == ProcID) {
        fully_other = false;
      }
    }
    if (!fully_other) {
      mesh_owned_elements.push_back(i);
      for (int j = 0; j < deg; j++) {
        auto local = all_conn[i * deg + j];
        keep_nodes[local] = true;
      }
    }
  }
  std::vector<int> external_nodes;
  std::vector<int> proc_node_counts(my_neighbors);
  std::vector<int> proc_node_idx(my_neighbors + 1);
  std::fill(proc_node_counts.begin(), proc_node_counts.end(), 0);
  for (int i = 0; i < my_neighbors; i++) {
    int proc = neighbors[i];
    for (int j = 0; j < vert_owners.size(); j++) {
      if (vert_owners[j] == proc && keep_nodes[j]) {
        external_nodes.push_back(j);
        proc_node_counts[i] += 1;
      }
    }
  }
  proc_node_idx[0] = 0;
  std::vector<int> node_cmap_ids;
  std::vector<int> node_cmap_to_neighbor;
  for (int i = 0; i < my_neighbors; i++) {
    proc_node_idx[i + 1] = proc_node_idx[i] + proc_node_counts[i];
    if (proc_node_counts[i] > 0) {
      node_cmap_ids.push_back(neighbors[i]);
      node_cmap_to_neighbor.push_back(i);
    }
  }
  std::vector<int> external_nodes_sorted(external_nodes.begin(), external_nodes.end());
  std::vector<int> boundary_nodes_sorted(boundary_nodes.begin(), boundary_nodes.end());
  std::sort(boundary_nodes_sorted.begin(), boundary_nodes_sorted.end());
  std::vector<int> internal_nodes;
  for (int j = 0; j < vert_owners.size(); j++) {
    if (vert_owners[j] == ProcID) {
      if (!std::binary_search(boundary_nodes_sorted.begin(), boundary_nodes_sorted.end(), j)) {
        internal_nodes.push_back(j);
      }
    }
  }

  std::vector<int> reduced_conn;
  std::vector<int> global_elem;
  std::vector<int> old_to_new_node_map(mesh->nverts());
  std::vector<int> old_to_new_elem_map(mesh->nelems());
  std::fill(old_to_new_node_map.begin(), old_to_new_node_map.end(), -1);
  std::fill(old_to_new_elem_map.begin(), old_to_new_elem_map.end(), -1);
  std::set<int> new_nodes;
  for (auto offset : mesh_owned_elements) {
    for (int j = 0; j < deg; j++) {
      reduced_conn.push_back(all_conn[offset * deg + j]);
    }
    global_elem.push_back(mesh->globals(mesh->dim())[offset]);
  }
  std::vector<int> new_nodes_v;
  new_nodes_v.reserve(internal_nodes.size() + boundary_nodes_sorted.size() + external_nodes.size());
  for (auto inn : internal_nodes) {
    new_nodes_v.push_back(inn);
  }
  for (auto bnn : boundary_nodes_sorted) {
    new_nodes_v.push_back(bnn);
  }
  for (auto enn : external_nodes) {
    new_nodes_v.push_back(enn);
  }
  std::vector<int> global_node(new_nodes_v.size());
  for (size_t i = 0; i < new_nodes_v.size(); i++) {
    global_node[i] = mesh->globals(0)[new_nodes_v[i]];
    old_to_new_node_map[new_nodes_v[i]] = i;
  }
  for (size_t i = 0; i < mesh_owned_elements.size(); i++) {
    old_to_new_elem_map[mesh_owned_elements[i]] = i;
  }

  for (size_t i = 0; i < reduced_conn.size(); i++) {
    auto old = reduced_conn[i];
    auto new_node = old_to_new_node_map[old];
    assert(new_node != -1);
    reduced_conn[i] = new_node + 1;
  }

  ex_put_init(exoid, title, dim, global_node.size(), global_elem.size(), nelem_blocks, nnode_sets,
              nside_sets);

  Few<Write<Real>, 3> coord_blk;
  for (Int i = 0; i < dim; ++i)
    coord_blk[i] = Write<Real>(new_nodes_v.size());
  auto coords = mesh->coords();
  for (size_t i = 0; i < new_nodes_v.size(); i++) {
    auto local_node = new_nodes_v[i];
    for (Int j = 0; j < dim; ++j) {
      coord_blk[j][i] = coords[local_node * dim + j];
    }
  }
  HostRead<Real> h_coord_blk[3];
  for (Int i = 0; i < dim; ++i)
    h_coord_blk[i] = HostRead<Real>(coord_blk[i]);
  CALL(ex_put_coord(exoid, h_coord_blk[0].data(), h_coord_blk[1].data(), h_coord_blk[2].data()));
  auto elems2file_idx = Write<LO>(mesh->nelems());
  auto elem_file_offset = LO(0);

  // TODO
  assert(region_set.size() == 1);

  for (auto block_id : region_set) {
    auto type_name = (dim == 3) ? "tetra4" : "tri3";
    auto elems_in_block = each_eq_to(elem_class_ids, block_id);
    auto block_elems2elem = collect_marked(elems_in_block);
    auto nblock_elems = block_elems2elem.size();
    if (verbose) {
      std::cout << "element block " << block_id << " has " << nblock_elems << " of type "
                << type_name << '\n';
    }
    CALL(ex_put_block(exoid, EX_ELEM_BLOCK, block_id, type_name, global_elem.size(), deg, 0, 0, 0));
    auto block_conn = read(unmap(block_elems2elem, all_conn, deg));
    auto block_conn_ex = add_to_each(block_conn, 1);
    auto h_block_conn = HostRead<LO>(block_conn_ex);
    CALL(ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, reduced_conn.data(), nullptr, nullptr));
    auto f = OMEGA_H_LAMBDA(LO block_elem) {
      elems2file_idx[block_elems2elem[block_elem]] = elem_file_offset + block_elem;
    };
    parallel_for(nblock_elems, f);
    elem_file_offset += nblock_elems;
  }

  std::vector<int> nset_global_node;
  std::vector<int> sset_global_side;
  if (1 && classify_with) {
    for (auto set_id : surface_set) {
      auto sides_in_set =
          land_each(each_eq_to(side_class_ids, set_id), each_eq_to(side_class_dims, I8(dim - 1)));
      if (classify_with & exodus::SIDE_SETS) {
        std::vector<int> sides2elem_new;
        std::vector<int> sides2local_new;
        auto set_sides2side = collect_marked(sides_in_set);
        auto nset_sides = set_sides2side.size();
        if (verbose) {
          std::cout << "side set " << set_id << " has " << nset_sides << " sides\n";
        }
        auto sides2elems = mesh->ask_up(dim - 1, dim);
        Write<int> set_sides2elem(nset_sides);
        Write<int> set_sides2local(nset_sides);
        for (int set_side = 0; set_side < nset_sides; set_side++) {
          auto side = set_sides2side[set_side];
          auto side_elem = sides2elems.a2ab[side];
          auto elem = sides2elems.ab2b[side_elem];
          auto elem_in_file = elems2file_idx[elem];
          auto code = sides2elems.codes[side_elem];
          auto which_down = code_which_down(code);
          set_sides2elem[set_side] = elem_in_file + 1;
          set_sides2local[set_side] = side_osh2exo(dim, which_down);
          if (std::find(mesh_owned_elements.begin(), mesh_owned_elements.end(), elem_in_file) !=
              mesh_owned_elements.end()) {
            sides2elem_new.push_back(old_to_new_elem_map[elem_in_file] + 1);
            sides2local_new.push_back(side_osh2exo(dim, which_down));
          }
        }
        auto h_set_sides2elem = HostRead<int>(set_sides2elem);
        auto h_set_sides2local = HostRead<int>(set_sides2local);
        CALL(ex_put_set_param(exoid, EX_SIDE_SET, set_id, sides2elem_new.size(), 0));
        sset_global_side.push_back(sides2elem_new.size());
        if (sides2elem_new.size() > 0) {
          CALL(ex_put_set(exoid, EX_SIDE_SET, set_id, sides2elem_new.data(),
                          sides2local_new.data()));
        }
      }
      if (classify_with & exodus::NODE_SETS) {
        auto nodes_in_set = mark_down(mesh, dim - 1, VERT, sides_in_set);
        auto set_nodes2node = collect_marked(nodes_in_set);
        auto set_nodes2node_ex = add_to_each(set_nodes2node, 1);
        auto nset_nodes = set_nodes2node.size();
        std::vector<int> nset_nodes_new;
        for (int i = 0; i < set_nodes2node.size(); i++) {
          int local_node = old_to_new_node_map[set_nodes2node[i]];
          if (local_node != -1 && vert_owners[set_nodes2node[i]] == ProcID) {
            nset_nodes_new.push_back(local_node + 1);
          }
        }
        if (verbose) {
          std::cout << "node set " << set_id << " has " << nset_nodes << " nodes\n";
        }
        auto h_set_nodes2node = HostRead<LO>(set_nodes2node_ex);
        nset_global_node.push_back(nset_nodes_new.size());
        CALL(ex_put_set_param(exoid, EX_NODE_SET, set_id, nset_nodes_new.size(), 0));
        if (nset_nodes_new.size() > 0) {
          CALL(ex_put_set(exoid, EX_NODE_SET, set_id, nset_nodes_new.data(), nullptr));
        }
      }
      std::vector<std::string> set_names(surface_set.size());
      for (auto &pair : mesh->class_sets) {
        auto &name = pair.first;
        for (auto &cp : pair.second) {
          if (cp.dim != I8(dim - 1))
            continue;
          std::size_t index = 0;
          for (auto surface_id : surface_set) {
            if (surface_id == cp.id) {
              set_names[index] = name;
              if (verbose && (classify_with & exodus::NODE_SETS)) {
                std::cout << "node set " << surface_id << " will be called \"" << name << "\"\n";
              }
              if (verbose && (classify_with & exodus::SIDE_SETS)) {
                std::cout << "side set " << surface_id << " will be called \"" << name << "\"\n";
              }
            }
            ++index;
          }
        }
      }
      std::vector<char *> set_name_ptrs(surface_set.size(), nullptr);
      for (std::size_t i = 0; i < set_names.size(); ++i) {
        if (set_names[i].empty()) {
          std::stringstream ss;
          ss << "ss_" << i;
          set_names[i] = ss.str();
        }
        set_name_ptrs[i] = const_cast<char *>(set_names[i].c_str());
      }
      if (classify_with & exodus::SIDE_SETS) {
        CALL(ex_put_names(exoid, EX_SIDE_SET, set_name_ptrs.data()));
      }
      for (std::size_t i = 0; i < set_names.size(); ++i) {
        std::stringstream ss;
        ss << "ns_" << i;
        set_names[i] = ss.str();
        set_name_ptrs[i] = const_cast<char *>(set_names[i].c_str());
      }
      if (classify_with & exodus::NODE_SETS) {
        CALL(ex_put_names(exoid, EX_NODE_SET, set_name_ptrs.data()));
      }
    }
  }
  std::vector<int> node_map(global_node.begin(), global_node.end());
  auto add1 = [](int &v) { v++; };
  std::for_each(node_map.begin(), node_map.end(), add1);
  std::vector<int> elem_map(global_elem.begin(), global_elem.end());
  std::for_each(elem_map.begin(), elem_map.end(), add1);

  ex_put_id_map(exoid, EX_NODE_MAP, node_map.data());
  ex_put_id_map(exoid, EX_ELEM_MAP, elem_map.data());

  int max_sets = surface_set.size();
  int gmax_sets;
  MPI_Allreduce(&max_sets, &gmax_sets, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  std::vector<int> all_sets(gmax_sets * Num_Proc + 1);
  std::fill(all_sets.begin(), all_sets.end(), -1);
  std::vector<int> my_sets(surface_set.begin(), surface_set.end());
  my_sets.resize(gmax_sets);
  MPI_Allgather(my_sets.data(), my_sets.size(), MPI_INT, all_sets.data(), gmax_sets, MPI_INT,
                MPI_COMM_WORLD);
  std::set<int> global_set_uniq;
  for (auto c : all_sets) {
    if (c != -1) {
      global_set_uniq.insert(c);
    }
  }

  std::vector<int> global_sets(global_set_uniq.begin(), global_set_uniq.end());
  std::sort(global_sets.begin(), global_sets.end());
  std::vector<int> global_node_counts(global_sets.size());
  std::vector<int> global_side_counts(global_sets.size());
  std::vector<int> df_counts(global_sets.size());
  std::fill(df_counts.begin(), df_counts.end(), 0);

  for (size_t setidx = 0; setidx < global_sets.size(); setidx++) {
    auto set = global_sets[setidx];
    auto my_offset = std::find(my_sets.begin(), my_sets.end(), set);
    if (my_offset != my_sets.end()) {
      int index = std::distance(my_sets.begin(), my_offset);
      MPI_Allreduce(&nset_global_node[index], &global_node_counts[setidx], 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&sset_global_side[index], &global_side_counts[setidx], 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
    }
  }
  CALL(ex_put_init_global(exoid, mesh->nglobal_ents(0), mesh->nglobal_ents(mesh->dim()), 1,
                          global_sets.size(), global_sets.size()));

  CALL(ex_put_ns_param_global(exoid, global_sets.data(), global_node_counts.data(),
                              df_counts.data()));
  CALL(ex_put_ss_param_global(exoid, global_sets.data(), global_side_counts.data(),
                              df_counts.data()));

  std::vector<int> eb_id(1);
  eb_id[0] = 1;
  std::vector<int> eb_counts(1);
  eb_counts[0] = mesh->nglobal_ents(mesh->dim());

  CALL(ex_put_eb_info_global(exoid, eb_id.data(), eb_counts.data()));
  char ftype = 'p';
  CALL(ex_put_init_info(exoid, Num_Proc, 1, &ftype));

  CALL(ex_put_loadbal_param(exoid, internal_nodes.size(), boundary_nodes_sorted.size(),
                            external_nodes.size(), global_elem.size(), 0, node_cmap_ids.size(), 0,
                            ProcID));

  std::vector<int> elem_internal(global_elem.size());
  std::iota(elem_internal.begin(), elem_internal.end(), 1);

  CALL(ex_put_processor_elem_maps(exoid, elem_internal.data(), NULL, ProcID));
  std::vector<int> node_mesh(global_node.size());
  std::iota(node_mesh.begin(), node_mesh.end(), 1);

  CALL(ex_put_processor_node_maps(exoid, node_mesh.data(), &node_mesh[internal_nodes.size()],
                                  &node_mesh[internal_nodes.size() + boundary_nodes_sorted.size()],
                                  ProcID));

  std::vector<int> cmap_node_counts;
  cmap_node_counts.reserve(node_cmap_to_neighbor.size());
  for (auto idx : node_cmap_to_neighbor) {
    cmap_node_counts.push_back(proc_node_counts[idx]);
  }

  CALL(
      ex_put_cmap_params(exoid, node_cmap_ids.data(), cmap_node_counts.data(), NULL, NULL, ProcID));

  std::vector<int> node_map_node_ids(external_nodes.size());
  std::vector<int> node_map_proc_ids(external_nodes.size());
  for (auto nidx : node_cmap_to_neighbor) {
    std::fill(node_map_proc_ids.begin() + proc_node_idx[nidx],
              node_map_proc_ids.begin() + proc_node_idx[nidx + 1], neighbors[nidx]);
  }
  for (size_t i = 0; i < external_nodes.size(); i++) {
    node_map_node_ids[i] = old_to_new_node_map[external_nodes[i]] + 1;
  }
  for (auto nidx : node_cmap_to_neighbor) {
    CALL(ex_put_node_cmap(exoid, neighbors[nidx], &node_map_node_ids[proc_node_idx[nidx]],
                          &node_map_proc_ids[proc_node_idx[nidx]], ProcID));
  }

  CALL(ex_close(exoid));

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    free(idv[imtrx]);
  }
  free(idv);
  idv = NULL;

  free_Surf_BC(First_Elem_Side_BC_Array, exo);
  free_Edge_BC(First_Elem_Edge_BC_Array, exo, dpi);
  free_nodes();
  if (Num_Proc == 1) {
    free_dpi_uni(dpi);
  } else {
    free_dpi(dpi);
  }
  free_exo(exo);
  init_exo_struct(exo);
  init_dpi_struct(dpi);

  if (rotation_allocated) {
    for (int i = 0; i < exo->num_nodes; i++) {
      for (int j = 0; j < NUM_VECTOR_EQUATIONS; j++) {
        free(rotation[i][j]);
      }
      free(rotation[i]);
      free(local_ROT_list[i]);
    }
    free(rotation);
    rotation = NULL;
    free(local_ROT_list);
    local_ROT_list = NULL;
    rotation_allocated = FALSE;
  }

  if (goma_automatic_rotations.rotation_nodes != NULL) {
    for (int i = 0; i < exo->num_nodes; i++) {
      // for (int j = 0; j < GOMA_MAX_NORMALS_PER_NODE; j++) {
      //   gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].normals[j]);
      //   gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].average_normals[j]);
      //   gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].tangent1s[j]);
      //   gds_vector_free((goma_automatic_rotations.rotation_nodes)[i].tangent2s[j]);
      // }
      for (int j = 0; j < DIM; j++) {
        goma_normal_free((goma_automatic_rotations.rotation_nodes)[i].rotated_coord[j]);
      }
    }
    free(goma_automatic_rotations.rotation_nodes);
    goma_automatic_rotations.rotation_nodes = NULL;
  }

  //  const char * tmpfile = "tmp.e";
  strncpy(ExoFile, "tmp.e", 127);
  strncpy(ExoFileOutMono, path, 127);
  strncpy(ExoFileOut, path, 127);
  multiname(ExoFileOut, ProcID, Num_Proc);
  int num_total_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (int i = 0; i < num_total_nodes; i++) {
      free(Local_Offset[imtrx][i]);
      free(Dolphin[imtrx][i]);
    }
    free(Dolphin[imtrx]);
    free(Local_Offset[imtrx]);
  }
  safer_free((void **)&Local_Offset);
  safer_free((void **)&Dolphin);

  read_mesh_exoII(exo, dpi);
  one_base(exo);
  wr_mesh_exo(exo, ExoFileOut, 0);
  zero_base(exo);

  if (Num_Proc > 1) {
    wr_dpi(dpi, ExoFileOut);
  }
  dpi->exodus_to_omega_h_node = (int *)malloc(sizeof(int) * global_node.size());
  for (size_t i = 0; i < old_to_new_node_map.size(); i++) {
    if (old_to_new_node_map[i] != -1) {
      dpi->exodus_to_omega_h_node[old_to_new_node_map[i]] = i;
    }
  }
  if (dpi->num_neighbors > 0) {
    for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      free(cx[imtrx]);
      free(Request);
      free(Status);
      cx[imtrx] = alloc_struct_1(Comm_Ex, DPI_ptr->num_neighbors);
      Request = alloc_struct_1(MPI_Request, Num_Requests * DPI_ptr->num_neighbors);
      Status = alloc_struct_1(MPI_Status, Num_Requests * DPI_ptr->num_neighbors);
    }
  }

  end_code();
}
} // namespace exodus
} // namespace Omega_h

void adapt_mesh(Omega_h::Mesh &mesh) {
  Omega_h::MetricInput genopts;
  //  genopts.sources.push_back(
  //      Omega_h::MetricSource{OMEGA_H_VARIATION, 1e-3, "phi", OMEGA_H_ISO_SIZE});
  genopts.sources.push_back(Omega_h::MetricSource{OMEGA_H_GIVEN, 1.0, "iso_size_metric",
                                                  OMEGA_H_ISO_SIZE, OMEGA_H_ABSOLUTE});
  //      genopts.sources.push_back(Omega_h::MetricSource{OMEGA_H_GIVEN, 1.0,
  //    "initial_metric",OMEGA_H_ISO_SIZE,OMEGA_H_ABSOLUTE});
  genopts.should_limit_lengths = true;
  genopts.min_length = 1e-6;
  genopts.max_length = 0.6;
  genopts.should_limit_gradation = true;
  genopts.max_gradation_rate = 0.3;
  Omega_h::add_implied_isos_tag(&mesh);
  Omega_h::generate_target_metric_tag(&mesh, genopts);

  Omega_h::AdaptOpts opts(&mesh);
  for (int j = V_FIRST; j < V_LAST; j++) {
    int imtrx = upd->matrix_index[j];
    if (imtrx >= 0) {
      if (j == MASS_FRACTION) {
        for (int mf = 0; mf < upd->Max_Num_Species; mf++) {
          std::string species_name = Exo_Var_Names[j].name2 + std::to_string(mf);
          opts.xfer_opts.type_map[species_name] = OMEGA_H_LINEAR_INTERP;
        }
      } else {
        opts.xfer_opts.type_map[Exo_Var_Names[j].name2] = OMEGA_H_LINEAR_INTERP;
      }
    }
  }
  for (int w = 0; w < efv->Num_external_field; w++) {
    opts.xfer_opts.type_map[efv->name[w]] = OMEGA_H_LINEAR_INTERP;
  }
  opts.max_length_allowed = 5;
  opts.max_length_desired = 1.8;
  opts.should_coarsen_slivers = true;
  opts.should_refine = true;
  opts.min_quality_desired = 0.5;
  int count = 1000;

  for (int i = 0; i < count; i++) {
    if (Omega_h::approach_metric(&mesh, opts)) {
      std::cout << "Adapt count " << i << "\n";
      Omega_h::adapt(&mesh, opts);
    } else {
      break;
    }
  }
  auto imb = mesh.imbalance();
  if (ProcID == 0) {
    std::cout << "Mesh imbalance = " << imb << "\n";
  }
  if (imb > 2) {
    GOMA_EH(-1, "Mesh imbalance %g exiting", imb);
  }
  // mesh.balance();
  // imb = mesh.imbalance();
  // if (ProcID == 0) {
  //   std::cout << "Mesh imbalance after balance = " << imb << "\n";
  // }
}

extern "C" {

void copy_solution(Dpi *dpi, double **x, Omega_h::Mesh &mesh) {
  for (int j = V_FIRST; j < V_LAST; j++) {
    int imtrx = upd->matrix_index[j];
    if (imtrx >= 0) {
      if (j == MASS_FRACTION) {
        for (int mf = 0; mf < upd->Max_Num_Species; mf++) {
          std::string species_name = Exo_Var_Names[j].name2 + std::to_string(mf);
          auto var_values = mesh.get_array<Omega_h::Real>(Omega_h::VERT, species_name);

          for (int i = 0; i < dpi->num_universe_nodes; i++) {
            int ja = Index_Solution(i, j, mf, 0, -2, imtrx);
            GOMA_EH(ja, "could not find solution");
            int index = i;
            if (dpi->exodus_to_omega_h_node != NULL) {
              index = dpi->exodus_to_omega_h_node[i];
            }
            x[imtrx][ja] = var_values[index];
          }
        }
      } else {
        auto var_values = mesh.get_array<Omega_h::Real>(Omega_h::VERT, Exo_Var_Names[j].name2);

        for (int i = 0; i < dpi->num_universe_nodes; i++) {
          int ja = Index_Solution(i, j, 0, 0, -2, imtrx);
          GOMA_EH(ja, "could not find solution");
          int index = i;
          if (dpi->exodus_to_omega_h_node != NULL) {
            index = dpi->exodus_to_omega_h_node[i];
          }
          x[imtrx][ja] = var_values[index];
        }
      }
    }
  }
  for (int w = 0; w < efv->Num_external_field; w++) {
    auto var_values = mesh.get_array<Omega_h::Real>(Omega_h::VERT, efv->name[w]);
    for (int i = 0; i < dpi->num_universe_nodes; i++) {
      int index = i;
      if (dpi->exodus_to_omega_h_node != NULL) {
        index = dpi->exodus_to_omega_h_node[i];
      }
      efv->ext_fld_ndl_val[w][i] = var_values[index];
    }
  }
}

// start with just level set field
void adapt_mesh_omega_h(struct GomaLinearSolverData **ams,
                        Exo_DB *exo,
                        Dpi *dpi,
                        double **x,
                        double **x_old,
                        double **x_older,
                        double **xdot,
                        double **xdot_old,
                        double **x_oldest,
                        double **resid_vector,
                        double **x_update,
                        double **scale,
                        int step) {

  static std::string base_name;
  static bool first_call = true;
  int argc = 0;
  char argv[1][8];
  char **argvptr = (char **)argv;
  auto lib = Omega_h::Library(&argc, &argvptr);
  auto classify_with = Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS;
  auto verbose = false;
#ifdef DEBUG_OMEGA_H
  verbose = true;
#endif
  Omega_h::Mesh mesh(&lib);
  Omega_h::exodus::convert_goma_to_omega_h(exo, dpi, x, &mesh, verbose);
  adapt_mesh(mesh);

  if (first_call) {
    base_name = std::string(ExoFileOutMono);
    first_call = false;
  }

  std::stringstream ss2;

  if (step == 0) {
    ss2 << base_name;
  } else {
    ss2 << base_name << "-s." << step;
  }

  Omega_h::exodus::convert_omega_h_to_goma(ss2.str().c_str(), &mesh, exo, dpi, verbose,
                                           classify_with);

  resetup_problem(exo, dpi);

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    int numProcUnknowns = NumUnknowns[imtrx] + NumExtUnknowns[imtrx];
    realloc_dbl_1(&x[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&x_old[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&x_older[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&x_update[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&xdot[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&xdot_old[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&x_oldest[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&resid_vector[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&scale[imtrx], numProcUnknowns, 0);
    realloc_dbl_1(&x_update[imtrx], numProcUnknowns + numProcUnknowns, 0);
    pg->matrices[imtrx].ams = ams[imtrx];
    pg->matrices[imtrx].x = x[imtrx];
    pg->matrices[imtrx].x_old = x_old[imtrx];
    pg->matrices[imtrx].x_older = x_older[imtrx];
    pg->matrices[imtrx].xdot = xdot[imtrx];
    pg->matrices[imtrx].xdot_old = xdot_old[imtrx];
    pg->matrices[imtrx].x_update = x_update[imtrx];
    pg->matrices[imtrx].scale = scale[imtrx];
    pg->matrices[imtrx].resid_vector = resid_vector[imtrx];
  }

  for (int w = 0; w < efv->Num_external_field; w++) {
    realloc_dbl_1(&efv->ext_fld_ndl_val[w],
                  dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes, 0);
  }
  resetup_matrix(ams, exo, dpi);
  copy_solution(dpi, x, mesh);
  step++;
}

} // extern "C"
#endif

// vim: expandtab sw=2 ts=8
