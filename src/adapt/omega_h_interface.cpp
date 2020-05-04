#ifdef HAVE_OMEGA_H
#include "adapt/omega_h_interface.h"

#include <algorithm>
#include <iostream>
#include <sstream>

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include <Omega_h_adapt.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>

extern "C" {
#define DISABLE_CPP
#ifndef MAX_PDIM
#define MAX_PDIM         3   /* Maximum physical problem dimension    */
#endif
#include "std.h"

#include "el_geom.h"

#include "rf_allo.h"

#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "rf_io.h"

#include "rf_mp.h"

#include "rf_masks.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "rf_element_storage_struct.h"
#include "mm_elem_block_structs.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_species.h"

#include "mm_fill_jac.h"
#include "mm_interface.h"

#include "mm_post_def.h"

#include "mm_eh.h"

#include "exo_struct.h"
#include "dpi.h"
#include "dp_types.h"
#include <mm_bc.h>
#include "exo_struct.h"
#define IGNORE_CPP_DEFINE
#include "sl_util_structs.h"
#undef IGNORE_CPP_DEFINE
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_unknown_map.h"
#include "adapt/resetup_problem.h"
#include "rf_io.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rf_allo.h"
#include "rd_mesh.h"
#include "wr_exo.h"
#include "rf_bc.h"
#include "rf_node_const.h"
extern int ***Local_Offset;
extern int 	***Dolphin;
extern int *NumUnknowns; /* Number of unknown variables updated by this   */
extern int *NumExtUnknowns; /* Number of unknown variables updated by this   */
extern int ***idv;
#undef DISABLE_CPP
}

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
      // seeing files from CUBIT with triangle sides in {3,4,5}...
      // no clue what thats about, just modulo and move on
      return (side+2) % 3;
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
    return -1; // needs to be filled in!
  }
  return -1;
}
static void get_elem_type_info(std::string const &type, int *p_dim, Omega_h_Family *p_family) {
  if (type == "tri3") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TRI") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TRI3") {
    *p_dim = 2;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "tetra4") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TETRA") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else if (type == "TET4") {
    *p_dim = 3;
    *p_family = OMEGA_H_SIMPLEX;
  } else {
    Omega_h_fail("Unsupported Exodus element type \"%s\"\n", type.c_str());
  }
}

void convert_to_omega_h_mesh(Exo_DB *exo, Dpi *dpi, int file, Mesh *mesh, bool verbose, int classify_with) {
  begin_code("exodus::read_mesh");
  ex_init_params init_params;
  memcpy(init_params.title, exo->title, 81);
  init_params.num_dim = exo->num_dim;
  init_params.num_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  init_params.num_elem = exo->num_elems;
  init_params.num_elem_blk = exo->num_elem_blocks;
  init_params.num_node_sets = exo->num_node_sets;
  init_params.num_side_sets = exo->num_side_sets;
  if (verbose) {
    std::cout << "init params:\n";
    std::cout << " Exodus ID " << file << '\n';
    std::cout << " Title " << init_params.title << '\n';
    std::cout << " num_dim " << init_params.num_dim << '\n';
    std::cout << " num_nodes " << init_params.num_nodes << '\n';
    std::cout << " num_elem " << init_params.num_elem << '\n';
    std::cout << " num_elem_blk " << init_params.num_elem_blk << '\n';
    std::cout << " num_node_sets " << init_params.num_node_sets << '\n';
    std::cout << " num_side_sets " << init_params.num_side_sets << '\n';
  }
  std::vector<int> block_ids(std::size_t(init_params.num_elem_blk));
  for (int i = 0; i < exo->num_elem_blocks; i++) {
    block_ids[i] = exo->eb_id[i];
  }

  //  CALL(ex_get_ids(file, EX_ELEM_BLOCK, block_ids.data()));
  std::vector<char> block_names_memory;
  std::vector<char *> block_names;
  setup_names(int(init_params.num_elem_blk), block_names_memory, block_names);
  //  CALL(ex_get_names(file, EX_ELEM_BLOCK, block_names.data()));
  HostWrite<LO> h_conn;
  Write<LO> elem_class_ids_w(LO(init_params.num_elem));
  LO elem_start = 0;
  int family_int = -1;
  int dim = -1;
  for (size_t i = 0; i < block_ids.size(); ++i) {
    char elem_type[MAX_STR_LENGTH + 1];
    elem_type[MAX_STR_LENGTH] = '\0';
    int nentries;
    int nnodes_per_entry;
    int nedges_per_entry = 0;
    int nfaces_per_entry = 0;
    memcpy(elem_type, exo->eb_elem_type[i], MAX_STR_LENGTH);
    nentries = exo->eb_num_elems[i];
    nnodes_per_entry = exo->eb_num_nodes_per_elem[i];

    //    CALL(ex_get_block(file, EX_ELEM_BLOCK, block_ids[i], elem_type, &nentries,
    //    &nnodes_per_entry,
    //                      &nedges_per_entry, &nfaces_per_entry, &nattr_per_entry));
    if (verbose) {
      std::cout << "block " << block_ids[i] << " \"" << block_names[i] << "\""
                << " has " << nentries << " elements of type " << elem_type << '\n';
    }
    /* some pretty weird blocks from the CDFEM people... */
    if (std::string("NULL") == elem_type && nentries == 0)
      continue;
    int dim_from_type;
    Omega_h_Family family_from_type;
    get_elem_type_info(elem_type, &dim_from_type, &family_from_type);
    if (family_int == -1)
      family_int = family_from_type;
    OMEGA_H_CHECK(family_int == family_from_type);
    if (dim == -1)
      dim = dim_from_type;
    OMEGA_H_CHECK(dim == dim_from_type);
    auto deg = element_degree(Omega_h_Family(family_int), dim, VERT);
    OMEGA_H_CHECK(nnodes_per_entry == deg);
    if (!h_conn.exists())
      h_conn = decltype(h_conn)(LO(init_params.num_elem * deg), "host connectivity");
    if (nedges_per_entry < 0)
      nedges_per_entry = 0;
    if (nfaces_per_entry < 0)
      nfaces_per_entry = 0;
    //    CALL(ex_get_conn(file, EX_ELEM_BLOCK, block_ids[i],
    //                     h_conn.data() + elem_start * nnodes_per_entry, NULL,
    //                     NULL));
    for (int j = 0; j < nentries * nnodes_per_entry; j++) {
      h_conn.data()[j + elem_start * nnodes_per_entry] = exo->eb_conn[i][j];
      //      if (h_conn[j+elem_start * nnodes_per_entry]-1 != exo->eb_conn[i][j]) {
      //        std::cout << "Differing conn " << h_conn[j+elem_start*nnodes_per_entry] << " != " <<
      //        exo->eb_conn[i][j];
      //      }
    }
    auto region_id = block_ids[i];
    auto f0 = OMEGA_H_LAMBDA(LO entry) { elem_class_ids_w[elem_start + entry] = region_id; };
    parallel_for(nentries, f0, "set_elem_class_ids");
    mesh->class_sets[block_names[i]].push_back({I8(dim), region_id});
    elem_start += nentries;
  }
  OMEGA_H_CHECK(elem_start == init_params.num_elem);
  Omega_h_Family family = Omega_h_Family(family_int);
  auto conn = LOs(h_conn.write());
  HostWrite<Real> h_coords(LO(init_params.num_nodes * dim));
  for (LO i = 0; i < exo->num_nodes; ++i) {
    for (Int j = 0; j < exo->num_dim; ++j) {
      h_coords[i * dim + 0] = exo->x_coord[i];
      h_coords[i * dim + 1] = exo->y_coord[i];
      if (exo->num_dim == 3)
        h_coords[i * dim + 2] = exo->z_coord[i];
    }
  }
  HostWrite<GO> vert_global;
  auto coords = Reals(h_coords.write());
  vert_global = decltype(vert_global)(GO(exo->num_nodes), "global vertices");
  for (auto i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
    vert_global.data()[i] = dpi->node_index_global[i];
  }
  //  build_from_elems_and_coords(mesh, OMEGA_H_SIMPLEX, dim, conn, coords);
  build_from_elems2verts(mesh, mesh->library()->world(), OMEGA_H_SIMPLEX, dim, conn,
                         GOs(vert_global));
  mesh->add_coords(coords);
  classify_elements(mesh);
  std::vector<int> side_set_ids(std::size_t(init_params.num_side_sets));
  for (int i = 0; i < exo->num_side_sets; i++) {
    side_set_ids[i] = exo->ss_id[i];
  }
  //  CALL(ex_get_ids(file, EX_SIDE_SET, side_set_ids.data()));
  Write<LO> side_class_ids_w(mesh->nents(dim - 1), -1);
  auto sides_are_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, sides_are_exposed);
  Write<I8> side_class_dims_w = deep_copy(mesh->get_array<I8>(dim - 1, "class_dim"));
  auto exposed_sides2side = collect_marked(sides_are_exposed);
  map_value_into(0, exposed_sides2side, side_class_ids_w);
  if ((classify_with & NODE_SETS) && init_params.num_node_sets) {
    int max_side_set_id = 0;
    if ((classify_with & SIDE_SETS) && side_set_ids.size()) {
      max_side_set_id = *std::max_element(side_set_ids.begin(), side_set_ids.end());
    }
    std::vector<int> node_set_ids(std::size_t(init_params.num_node_sets));
    //    CALL(ex_get_ids(file, EX_NODE_SET, node_set_ids.data()));
    for (int i = 0; i < exo->num_node_sets; i++) {
      node_set_ids[i] = exo->ns_id[i];
    }
    std::vector<char> names_memory;
    std::vector<char *> name_ptrs;
    setup_names(int(init_params.num_node_sets), names_memory, name_ptrs);
    //    CALL(ex_get_names(file, EX_NODE_SET, name_ptrs.data()));

    for (size_t i = 0; i < node_set_ids.size(); ++i) {
      int nentries;
//      CALL(ex_get_set_param(file, EX_NODE_SET, node_set_ids[i], &nentries, &ndist_factors));
      nentries = exo->ns_num_nodes[i];
      if (verbose) {
        std::cout << "node set " << node_set_ids[i] << " has " << nentries << " nodes\n";
      }
      HostWrite<LO> h_set_nodes2nodes(nentries);
//      CALL(ex_get_set(file, EX_NODE_SET, node_set_ids[i], h_set_nodes2nodes.data(), nullptr));
      auto f0 = OMEGA_H_LAMBDA(LO index) { h_set_nodes2nodes[index] = exo->ns_node_list[exo->ns_node_index[i] + index]; };
      parallel_for(nentries, f0);
      auto set_nodes2nodes = LOs(h_set_nodes2nodes.write());
      auto nodes_are_in_set = mark_image(set_nodes2nodes, mesh->nverts());
      auto sides_are_in_set = mark_up_all(mesh, VERT, dim - 1, nodes_are_in_set);
      auto set_sides2side = collect_marked(sides_are_in_set);
      auto surface_id = node_set_ids[i] + max_side_set_id;
      if (verbose) {
        std::cout << "node set #" << node_set_ids[i] << " \"" << name_ptrs[i]
                  << "\" will be surface " << surface_id << '\n';
      }
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
      mesh->class_sets[name_ptrs[i]].push_back({I8(dim - 1), surface_id});
    }
  }
  if (classify_with & SIDE_SETS) {
    std::vector<char> names_memory;
    std::vector<char *> name_ptrs;
    setup_names(int(init_params.num_side_sets), names_memory, name_ptrs);
//    CALL(ex_get_names(file, EX_SIDE_SET, name_ptrs.data()));
    for (size_t i = 0; i < side_set_ids.size(); ++i) {
      int nentries;
//      CALL(ex_get_set_param(file, EX_SIDE_SET, side_set_ids[i], &nentries, &ndist_factors));
      nentries = exo->ss_num_sides[i];
      if (verbose) {
        std::cout << "side set #" << side_set_ids[i] << " \"" << name_ptrs[i] << "\" has "
                  << nentries << " sides, will be surface " << side_set_ids[i] << "\n";
      }
      HostWrite<LO> h_set_sides2elem(nentries);
      HostWrite<LO> h_set_sides2local(nentries);
//      CALL(ex_get_set(file, EX_SIDE_SET, side_set_ids[i], h_set_sides2elem.data(),
//                      h_set_sides2local.data()));
      auto f0 = OMEGA_H_LAMBDA(LO index) {
        int offset = exo->ss_elem_index[i];
        h_set_sides2elem[index] = exo->ss_elem_list[offset + index];
        h_set_sides2local[index] = exo->ss_side_list[offset+index];
      };
      parallel_for(nentries, f0);
      auto set_sides2elem = LOs(h_set_sides2elem.write());
      auto set_sides2local = LOs(h_set_sides2local.write());
      auto elems2sides = mesh->ask_down(dim, dim - 1).ab2b;
      auto nsides_per_elem = element_degree(family, dim, dim - 1);
      auto set_sides2side_w = Write<LO>(nentries);
      auto f2 = OMEGA_H_LAMBDA(LO set_side) {
        auto elem = set_sides2elem[set_side];
        auto side_of_element = side_exo2osh(family, dim, set_sides2local[set_side]);
        OMEGA_H_CHECK(side_of_element != -1);
        auto side = elems2sides[elem * nsides_per_elem + side_of_element];
        set_sides2side_w[set_side] = side;
      };
      parallel_for(nentries, f2, "set_sides2side");
      auto set_sides2side = LOs(set_sides2side_w);
      auto surface_id = side_set_ids[i];
      map_value_into(surface_id, set_sides2side, side_class_ids_w);
      map_value_into(I8(dim - 1), set_sides2side, side_class_dims_w);
      mesh->class_sets[name_ptrs[i]].push_back({I8(dim - 1), surface_id});
    }
  }
  auto elem_class_ids = LOs(elem_class_ids_w);
  auto side_class_ids = LOs(side_class_ids_w);
  auto side_class_dims = Read<I8>(side_class_dims_w);
  mesh->add_tag(dim, "class_id", 1, elem_class_ids);
  mesh->add_tag(dim - 1, "class_id", 1, side_class_ids);
  mesh->set_tag(dim - 1, "class_dim", side_class_dims);
  finalize_classification(mesh);
  end_code();
}

static OMEGA_H_INLINE int side_osh2exo(int dim, int side) {
  switch (dim) {
    case 2:
      switch (side) {
        case 1:
          return 1;
        case 2:
          return 2;
        case 0:
          return 3;
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
void convert_back_to_goma_exo(
    const char * path, Mesh* mesh, Exo_DB *exo, Dpi *dpi, bool verbose, int classify_with) {

  Omega_h::exodus::write("tmp.e", mesh, true, classify_with);



     for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
        {
          free(idv[imtrx]);
        }
        free(idv);
        idv = NULL;



  free_Surf_BC(First_Elem_Side_BC_Array, exo);
  free_Edge_BC(First_Elem_Edge_BC_Array, exo, dpi);
  free_nodes();
  free_dpi_uni(dpi);
  free_exo(exo);
  init_exo_struct(exo);
  init_dpi_struct(dpi);

//  const char * tmpfile = "tmp.e";
  strncpy(ExoFile, "tmp.e", 128);
  strncpy(ExoFileOutMono, path, 128);
  strncpy(ExoFileOut, path, 128);
  int num_total_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (int i = 0; i < num_total_nodes; i++) {
      free(Local_Offset[imtrx][i]);
      free(Dolphin[imtrx][i]);
    }
    free(Dolphin[imtrx]);
    free(Local_Offset[imtrx]);
  }
  safer_free((void **) &Local_Offset);
  safer_free((void **) &Dolphin);

  read_mesh_exoII(exo, dpi);
  one_base(exo);
  wr_mesh_exo(exo, ExoFileOutMono, 0);
  zero_base(exo);


  end_code();
}
// int open(filesystem::path const& path, bool verbose) {
//  auto comp_ws = int(sizeof(Real));
//  int io_ws = 0;
//  float version;
//  auto mode = EX_READ | EX_MAPS_INT64_API;
//  auto exodus_file = ex_open(path.c_str(), mode, &comp_ws, &io_ws, &version);
//  if (exodus_file < 0)
//    Omega_h_fail("can't open Exodus file %s\n", path.c_str());
//  if (verbose) {
//    std::cout << "ex_open(" << path << ")\n";
//    std::cout << "  comp_ws: " << comp_ws << '\n';
//    std::cout << "  io_ws: " << io_ws << '\n';
//    std::cout << "  version: " << version << '\n';
//  }
//  return exodus_file;
//}
} // namespace exodus
} // namespace Omega_h

#if 0
static int exo_side_to_osh(int side, int dim) {
  switch (dim) {
  case 2:
    return (side + 1) % 3;
  case 3:
    return (side + 1) % 4;
  default:
    EH(GOMA_ERROR, "Unknown dim exo_side_to_osh");
    return -1;
  }
}

static double smooth_H(double F, double eps) {
  if (F > eps) {
    return 1;
  } else if (F < -eps) {
    return 0;
  }
  return 0.5 * (1. + F / eps + sin(M_PI * F / eps) / M_PI);
}

static double indicator(double phi, double eps) { return (1.0 / eps) * phi * (1 - phi); }
#endif

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
        opts.xfer_opts.type_map[Exo_Var_Names[j].name2] = OMEGA_H_LINEAR_INTERP;
    }
  }
  for (int w=0; w<efv->Num_external_field; w++) {
    opts.xfer_opts.type_map[efv->name[w]] = OMEGA_H_LINEAR_INTERP;
  }
  opts.max_length_allowed = 5;
  opts.max_length_desired = 1.6;
  opts.should_coarsen_slivers = true;
  opts.should_refine = true;
  opts.min_quality_desired = 0.6;
  int count = 1000;

  for (int i = 0; i < count; i++) {
    if (Omega_h::approach_metric(&mesh, opts)) {
      Omega_h::adapt(&mesh, opts);
    } else {
      break;
    }
  }
}

extern "C" {

void copy_solution(Exo_DB *exo, Dpi *dpi, double **x, Omega_h::Mesh &mesh) {
  for (int j = V_FIRST; j < V_LAST; j++) {
    int imtrx = upd->matrix_index[j];
    if (imtrx >= 0) {
      auto var_values = mesh.get_array<Omega_h::Real>(Omega_h::VERT, Exo_Var_Names[j].name2);

      for (int i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
        int ja = Index_Solution(i, j, 0, 0, -2, imtrx);
        EH(ja, "could not find solution");
        x[imtrx][ja] = var_values[i];
      }
    }
  }
  for (int w=0; w<efv->Num_external_field; w++) {
    auto var_values = mesh.get_array<Omega_h::Real>(Omega_h::VERT, efv->name[w]);
      for (int i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
        efv->ext_fld_ndl_val[w][i] = var_values[i];
      }
  }
}

  // start with just level set field
void adapt_mesh_omega_h(struct Aztec_Linear_Solver_System **ams,
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
  static auto lib = Omega_h::Library();
  auto classify_with =  Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS;
  auto verbose = true;
  Omega_h::Mesh mesh(&lib);
  auto exodus_file = 0; // Omega_h::exodus::open("cup.g", verbose);
  Omega_h::exodus::convert_to_omega_h_mesh(exo, dpi, exodus_file, &mesh, verbose, classify_with);

  for (int j = V_FIRST; j < V_LAST; j++) {
    int imtrx = upd->matrix_index[j];
    if (imtrx >= 0) {
      auto var_values = Omega_h::Write<Omega_h::Real>(mesh.nverts());
      for (int i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
        int ja = Index_Solution(i, j, 0, 0, -2, imtrx);
        EH(ja, "could not find solution");
        var_values[i] = x[imtrx][ja];
      }
      mesh.add_tag(Omega_h::VERT, Exo_Var_Names[j].name2, 1, Omega_h::Reals(var_values));
      if (j == FILL) {
        //        auto H_values = Omega_h::Write<Omega_h::Real>(mesh.nverts());
        //        auto f0 = OMEGA_H_LAMBDA(Omega_h::LO index) {
        //          H_values[index] =
        //              indicator(smooth_H(var_values[index], ls->Length_Scale), ls->Length_Scale);
        //        };
        //        Omega_h::parallel_for(mesh.nverts(), f0, "set_indicator_values");
        //        mesh.add_tag(Omega_h::VERT, "indicator", 1, Omega_h::Reals(H_values));
        auto target_metrics =
            Omega_h::Write<Omega_h::Real>(mesh.nverts() * Omega_h::symm_ncomps(mesh.dim()));
        auto f0 = OMEGA_H_LAMBDA(Omega_h::LO index) {
          auto F = var_values[index];
          auto iso_size = 0.3;
          if (std::abs(F) < 0.5) {
            iso_size = 0.02;
          }
          auto target_metric = Omega_h::compose_metric(Omega_h::identity_matrix<2, 2>(),
                                                       Omega_h::vector_2(iso_size, iso_size));
          Omega_h::set_vector(target_metrics, index, Omega_h::symm2vector(target_metric));
        };

        Omega_h::parallel_for(mesh.nverts(), f0, "set_iso_metric_values");
        mesh.add_tag(Omega_h::VERT, "iso_size_metric", Omega_h::symm_ncomps(mesh.dim()),
                     Omega_h::Reals(target_metrics));
      }
    }
  }
  for (int w=0; w<efv->Num_external_field; w++) {
      auto var_values = Omega_h::Write<Omega_h::Real>(mesh.nverts());
      for (int i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
        var_values[i] = efv->ext_fld_ndl_val[w][i];
      }
      mesh.add_tag(Omega_h::VERT, efv->name[w], 1, Omega_h::Reals(var_values));
  }

  auto writer = Omega_h::vtk::Writer("transfer.vtk", &mesh);
  writer.write(step);
  adapt_mesh(mesh);

  std::string filename;
  std::stringstream ss;

  ss << "adapt." << step << ".vtk";

  auto writer_adapt = Omega_h::vtk::Writer(ss.str(), &mesh);
  writer_adapt.write(step);

  if (first_call) {
    base_name = std::string(ExoFileOutMono);
    first_call = false;
  }

  std::stringstream ss2;

  ss2 << base_name << "-s." << step;

  Omega_h::exodus::convert_back_to_goma_exo(ss2.str().c_str(), &mesh, exo, dpi, true, classify_with);

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
    realloc_dbl_1(&efv->ext_fld_ndl_val[w], dpi->num_internal_nodes + dpi->num_external_nodes, 0);
  }
  resetup_matrix(ams, exo, dpi);
  copy_solution(exo, dpi, x, mesh);
  step++;
}


} // extern "C"
#endif

// vim: expandtab sw=2 ts=8
