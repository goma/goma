/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "brkfix/setup_fix_data.h"

// include mpi before since it might get pulled in from Exodus
// an error might occur because of mpicxx headers being included as C
#include <mpi.h>

extern "C" {
#define DISABLE_CPP
#include "base_mesh.h"
#include "brkfix/bbb.h"
#include "brkfix/fix.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_mp_structs.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_element_storage_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "shell_tfmp_struct.h"
#include "shell_tfmp_util.h"
#include "std.h"
#include "wr_exo.h"
#undef DISABLE_CPP
}

#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>

extern "C" void
setup_fix_data(const char *mono_name, int num_procs, struct fix_data *fd, int *pmax) {

  std::vector<std::unordered_set<int>> ns_nodes;
  std::vector<std::vector<std::pair<int, int>>> ss_elem_sides;

  char monolith_file_name[FILENAME_MAX_ACK]; /* original mesh */
  char polylith_name[FILENAME_MAX_ACK];      /* "basename_1of2.exoII" */

  int num_node_var_max = 0;
  int ss_elem_count = 0;
  for (int p = 0; p < num_procs; p++) {
    strcpy(polylith_name, mono_name);

    strcpy(monolith_file_name, polylith_name);

    multiname(polylith_name, p, num_procs);

#ifdef DEBUG
    /*  err = sscanf(line, "%s", in_exodus_file_name); */
    fprintf(stderr, "monolith EXODUSII file = \"%s\"\n", monolith_file_name);
#endif

    /*
     * Make preliminary space for the monolithic database and the
     * polylithic pieces.
     *
     * The strategy will be to look closely at the first polylithic piece
     * in order to ascribe basic sizing information to the monolith.
     *
     * Then, the polyliths are traversed to slowly build the mesh and the
     * results of the monolith. Several passes may be necessary...
     */

    Exo_DB *poly = alloc_struct_1(Exo_DB, 1);
    Dpi *dpin = alloc_struct_1(Dpi, 1);

    /*
     * Open the first polylith to determine if there are many timeplanes
     * of data.
     *
     * Assume that the first polylith is representative, that the number
     * of nodal variables, the number of timeplanes, etc will be the same
     * for each and every polylith encountered hereafter.
     *
     * The first polylith is important, too, in that much of the sizing
     * information for the global problem will be derived from it.
     */

    init_exo_struct(poly);

    init_dpi_struct(dpin);

    /*
     * Some defaults are different for polylithic pieces - they have use these
     * maps to relate themselves to the monolith...
     *
     * Now, these maps are part of the Dpi information and no longer piggybacked
     * inside EXODUS II. Thus, don't attempt to read them if they aren't there.
     */

#ifdef DEBUG
    fprintf(stderr, "Fix: attempting to build a %d piece %s\n", num_procs, monolith_file_name);
#endif

    /*
     * Read everything in the 1st polylith's EXODUS information except for
     * the results data per se...
     *
     * Set to zero the variables that get used for allocating and writing...
     */

    rd_exo(poly, polylith_name, 0,
           (EXODB_ACTION_RD_INIT + EXODB_ACTION_RD_MESH + EXODB_ACTION_RD_RES0));

    if (poly->num_node_vars > num_node_var_max) {
      num_node_var_max = poly->num_node_vars;
      *pmax = p;
    }

    zero_base(poly);
    setup_base_mesh(dpin, poly, 1);

    rd_dpi(poly, dpin, polylith_name, false);

    if (p == 0) {
      ns_nodes.resize(dpin->num_node_sets_global);
      ss_elem_sides.resize(dpin->num_side_sets_global);

      fd->num_elem_blocks = dpin->num_elem_blocks_global;
      fd->ns_node_len_global = 0;
      for (int i = 0; i < dpin->num_node_sets_global; i++) {
        fd->ns_node_len_global += dpin->num_ns_global_node_counts[i];
        fd->ns_distfact_len_global += dpin->num_ns_global_df_counts[i];
      }
      for (int i = 0; i < dpin->num_side_sets_global; i++) {
        fd->ss_elem_len_global += dpin->num_ss_global_side_counts[i];
        fd->ss_distfact_len_global += dpin->num_ss_global_df_counts[i];
      }

      fd->eb_elem_type_global = (char **)malloc(sizeof(char *) * dpin->num_elem_blocks_global);
      fd->eb_index_global = (int *)malloc(sizeof(int) * dpin->num_elem_blocks_global);
      for (int i = 0; i < dpin->num_elem_blocks_global; i++) {
        fd->eb_index_global[i] = i;
        fd->eb_elem_type_global[i] = (char *)malloc(sizeof(char) * (MAX_STR_LENGTH + 1));
        strcpy(fd->eb_elem_type_global[i], "NULL");
      }

      fd->ns_node_index_global = (int *)calloc(sizeof(int), dpin->num_node_sets_global + 1);
      fd->ns_distfact_index_global = (int *)calloc(sizeof(int), dpin->num_node_sets_global + 1);
      fd->ns_node_index_global[0] = 0;
      fd->ns_distfact_index_global[0] = 0;
      for (int i = 0; i < dpin->num_node_sets_global; i++) {
        fd->ns_node_index_global[i + 1] =
            fd->ns_node_index_global[i] + dpin->num_ns_global_node_counts[i];
        fd->ns_distfact_index_global[i + 1] =
            fd->ns_distfact_index_global[i] + dpin->num_ns_global_df_counts[i];
      }
      fd->ss_elem_index_global = (int *)calloc(sizeof(int), dpin->num_side_sets_global + 1);
      fd->ss_distfact_index_global = (int *)calloc(sizeof(int), dpin->num_side_sets_global + 1);
      fd->ss_elem_index_global[0] = 0;
      fd->ss_distfact_index_global[0] = 0;
      for (int i = 0; i < dpin->num_side_sets_global; i++) {
        fd->ss_elem_index_global[i + 1] =
            fd->ss_elem_index_global[i] + dpin->num_ss_global_side_counts[i];
        fd->ss_distfact_index_global[i + 1] =
            fd->ss_distfact_index_global[i] + dpin->num_ss_global_df_counts[i];
      }

      fd->eb_num_nodes_per_elem = (int *)calloc(sizeof(int), dpin->num_elem_blocks_global);
      for (int i = 0; i < dpin->num_elem_blocks_global; i++) {
        fd->eb_num_nodes_per_elem[i] = 0;
      }

      fd->ns_node_list = NULL;
      fd->ns_distfact_list = NULL;
      if (dpin->num_node_sets_global > 0) {
        fd->ns_node_list = (int *)malloc(sizeof(int) * fd->ns_node_len_global);
        fd->ns_distfact_list = (double *)malloc(sizeof(double) * fd->ns_distfact_len_global);
        for (int i = 0; i < fd->ns_distfact_len_global; i++) {
          fd->ns_distfact_list[i] = 0.0;
        }
      }
      fd->ss_elem_list = NULL;
      fd->ss_side_list = NULL;
      fd->ss_distfact_list = NULL;
      if (dpin->num_side_sets_global > 0) {
        fd->ss_elem_list = (int *)malloc(sizeof(int) * fd->ss_elem_len_global);
        fd->ss_side_list = (int *)malloc(sizeof(int) * fd->ss_elem_len_global);
        fd->ss_distfact_list = (double *)malloc(sizeof(double) * fd->ss_distfact_len_global);
        // set distfact to 0.0;
        for (int i = 0; i < fd->ss_distfact_len_global; i++) {
          fd->ss_distfact_list[i] = 0.0;
        }
      }
    }

    for (int i = 0; i < poly->num_elem_blocks; i++) {
      if (strlen(poly->eb_elem_type[i]) > 0 && (strcmp(poly->eb_elem_type[i], "NULL") != 0)) {
        strcpy(fd->eb_elem_type_global[i], poly->eb_elem_type[i]);
      }
      fd->eb_num_nodes_per_elem[i] =
          MAX(fd->eb_num_nodes_per_elem[i], poly->eb_num_nodes_per_elem[i]);
    }

    for (int i = 0; i < poly->num_node_sets; i++) {
      for (int j = 0; j < poly->ns_num_nodes[i]; j++) {
        ns_nodes[i].insert(dpin->node_index_global[poly->ns_node_list[poly->ns_node_index[i] + j]]);
      }
    }
    for (int i = 0; i < poly->num_side_sets; i++) {
      for (int j = 0; j < poly->ss_num_sides[i]; j++) {
        int elem = dpin->elem_index_global[poly->ss_elem_list[poly->ss_elem_index[i] + j]];
        int side = poly->ss_side_list[poly->ss_elem_index[i] + j];
        ss_elem_sides[i].push_back(std::make_pair(elem, side));
        ss_elem_count++;
      }
    }

    free_dpi(dpin);
    free(dpin);

    free_element_blocks(poly);

    free_exo(poly);
    free(poly);
  }
  // lets check if our sides and nodes match
  std::vector<std::vector<int>> ns_nodes_vec;
  ns_nodes_vec.resize(ns_nodes.size());
  int ns_nodes_count = 0;
  for (unsigned int i = 0; i < ns_nodes.size(); i++) {
    ns_nodes_vec[i].assign(ns_nodes[i].begin(), ns_nodes[i].end());
    ns_nodes_count += ns_nodes_vec[i].size();
  }
  GOMA_ASSERT(ns_nodes_count == fd->ns_node_len_global);
  GOMA_ASSERT(ss_elem_count == fd->ss_elem_len_global);
  int offset = 0;
  for (unsigned int i = 0; i < ns_nodes.size(); i++) {
    for (unsigned int j = 0; j < ns_nodes[i].size(); j++) {
      fd->ns_node_list[offset] = ns_nodes_vec[i][j];
      offset++;
    }
  }

  offset = 0;
  for (unsigned int i = 0; i < ss_elem_sides.size(); i++) {
    for (unsigned int j = 0; j < ss_elem_sides[i].size(); j++) {
      int elem = ss_elem_sides[i][j].first;
      int side = ss_elem_sides[i][j].second;
      fd->ss_elem_list[offset] = elem;
      fd->ss_side_list[offset] = side;
      offset++;
    }
  }
}

extern "C" void free_fix_data(struct fix_data *fd) {

  free(fd->eb_num_nodes_per_elem);
  free(fd->ns_node_index_global);
  free(fd->ns_distfact_index_global);
  free(fd->ss_elem_index_global);
  free(fd->ss_distfact_index_global);
  for (int i = 0; i < fd->num_elem_blocks; i++) {
    free(fd->eb_elem_type_global[i]);
  }

  free(fd->eb_elem_type_global);
  free(fd->eb_index_global);
  free(fd->ns_node_list);
  free(fd->ns_distfact_list);
  free(fd->ss_elem_list);
  free(fd->ss_side_list);
  free(fd->ss_distfact_list);
}
