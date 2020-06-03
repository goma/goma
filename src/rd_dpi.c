/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* rd_dpi.c -- routines for reading distributed processing information
 * Uses Nemesis structure from SEACAS
 */

#define GOMA_RD_DPI_C

#include <stdio.h>
#include <stdlib.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

#include "dpi.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "mm_eh.h"
#include "rd_dpi.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_mp.h"
#include "std.h"

// Helper for exodus return values
#define CHECK_EX_ERROR(err, format, ...)                              \
  do {                                                                \
    if (err < 0) {                                                    \
      goma_eh(GOMA_ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__); \
    }                                                                 \
  } while (0)

int rd_dpi(Exo_DB *exo, Dpi *d, char *fn) {
  init_dpi_struct(d);
  float version = -4.98; /* initialize. ex_open() changes this. */
  int comp_wordsize = sizeof(dbl);
  int io_wordsize = 0;
  int exoid = ex_open(fn, EX_READ, &comp_wordsize, &io_wordsize, &version);

  int ex_error;

  ex_error = ex_get_init_global(exoid, &d->num_nodes_global, &d->num_elems_global,
                                &d->num_elem_blocks_global, &d->num_node_sets_global,
                                &d->num_side_sets_global);

  CHECK_EX_ERROR(ex_error, "ex_get_init_global");

  // Node Set Global
  d->ns_id_global = alloc_int_1(d->num_node_sets_global, 0);
  d->num_ns_global_node_counts = alloc_int_1(d->num_node_sets_global, 0);
  d->num_ns_global_df_counts = alloc_int_1(d->num_node_sets_global, 0);
  ex_error = ex_get_ns_param_global(exoid, d->ns_id_global, d->num_ns_global_node_counts,
                                    d->num_ns_global_df_counts);

  // Side Set Global
  d->ss_id_global = alloc_int_1(d->num_side_sets_global, 0);
  d->num_ss_global_side_counts = alloc_int_1(d->num_side_sets_global, 0);
  d->num_ss_global_df_counts = alloc_int_1(d->num_side_sets_global, 0);
  ex_error = ex_get_ns_param_global(exoid, d->ss_id_global, d->num_ss_global_side_counts,
                                    d->num_ss_global_df_counts);

  CHECK_EX_ERROR(ex_error, "ex_get_ns_param_global");
  // Block Global
  d->global_elem_block_ids = alloc_int_1(d->num_elem_blocks_global, 0);
  d->global_elem_block_counts = alloc_int_1(d->num_elem_blocks_global, 0);
  ex_error = ex_get_eb_info_global(exoid, d->global_elem_block_ids, d->global_elem_block_counts);
  CHECK_EX_ERROR(ex_error, "ex_get_eb_info_global");

  // Nemesis Info
  ex_error = ex_get_init_info(exoid, &d->num_proc, &d->num_proc_in_file, &d->ftype);
  CHECK_EX_ERROR(ex_error, "ex_get_init_info");
  if (d->num_proc != Num_Proc) {
    EH(GOMA_ERROR, "Nemesis mesh error num_proc != number of mpi processes");
  }
  if (d->num_proc_in_file != 1) {
    EH(GOMA_ERROR, "Nemesis mesh error expected num_proc_in_file == 1");
  }
  if (d->num_proc != Num_Proc) {
    EH(GOMA_ERROR, "Nemesis mesh error num_proc != number of mpi processes");
  }
  if (d->ftype != 'p') {
    EH(GOMA_ERROR, "Nemesis mesh error ftype expected 'p'");
  }
  // global indices
  d->node_index_global = alloc_int_1(exo->num_nodes, 0);
  d->elem_index_global = alloc_int_1(exo->num_elems, 0);
  ex_error = ex_get_id_map(exoid, EX_NODE_MAP, d->node_index_global);
  CHECK_EX_ERROR(ex_error, "ex_get_id_map EX_NODE_MAP");
  ex_error = ex_get_id_map(exoid, EX_ELEM_MAP, d->elem_index_global);
  CHECK_EX_ERROR(ex_error, "ex_get_id_map EX_ELEM_MAP");

  // Load Balance Information
  ex_error = ex_get_loadbal_param(
      exoid, &d->num_internal_nodes, &d->num_boundary_nodes, &d->num_external_nodes,
      &d->num_internal_elems, &d->num_border_elems, &d->num_node_cmaps, &d->num_elem_cmaps, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_loadbal_param");

  d->num_owned_nodes = d->num_internal_nodes + d->num_boundary_nodes;

  d->proc_elem_internal = alloc_int_1(d->num_internal_elems, 0);
  if (d->num_border_elems > 0) {
    d->proc_elem_border = alloc_int_1(d->num_border_elems, 0);
  }
  ex_error = ex_get_processor_elem_maps(exoid, d->proc_elem_internal, d->proc_elem_border, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_processor_elem_maps");

  d->proc_node_internal = alloc_int_1(d->num_internal_nodes, 0);
  if (d->num_boundary_nodes > 0) {
    d->proc_node_boundary = alloc_int_1(d->num_boundary_nodes, 0);
  }
  if (d->num_external_nodes > 0) {
    d->proc_node_external = alloc_int_1(d->num_external_nodes, 0);
  }

  ex_error = ex_get_processor_node_maps(exoid, d->proc_node_internal, d->proc_node_boundary,
                                        d->proc_node_external, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_processor_node_maps");

  d->node_cmap_ids = alloc_int_1(d->num_node_cmaps, 0);
  d->node_cmap_node_counts = alloc_int_1(d->num_node_cmaps, 0);
  if (d->num_elem_cmaps > 0) {
    d->elem_cmap_ids = alloc_int_1(d->num_elem_cmaps, 0);
    d->elem_cmap_elem_counts = alloc_int_1(d->num_elem_cmaps, 0);
  }

  ex_error = ex_get_cmap_params(exoid, d->node_cmap_ids, d->node_cmap_node_counts, d->elem_cmap_ids,
                                d->elem_cmap_elem_counts, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_get_cmake_params");

  d->node_map_node_ids = calloc(sizeof(int *), d->num_node_cmaps);
  d->node_map_proc_ids = calloc(sizeof(int *), d->num_node_cmaps);

  for (int i = 0; i < d->num_node_cmaps; i++) {
    d->node_map_node_ids[i] = alloc_int_1(d->node_cmap_node_counts[i], 0);
    d->node_map_proc_ids[i] = alloc_int_1(d->node_cmap_node_counts[i], 0);
    ex_error = ex_get_node_cmap(exoid, d->node_cmap_ids[i], d->node_map_node_ids[i],
                                d->node_map_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_get_node_cmap %d", i);
  }

  if (d->num_elem_cmaps > 0) {
    d->elem_cmap_elem_ids = calloc(sizeof(int *), d->num_elem_cmaps);
    d->elem_cmap_side_ids = calloc(sizeof(int *), d->num_elem_cmaps);
    d->elem_cmap_proc_ids = calloc(sizeof(int *), d->num_elem_cmaps);
  }

  for (int i = 0; i < d->num_elem_cmaps; i++) {
    d->elem_cmap_elem_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    d->elem_cmap_side_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    d->elem_cmap_proc_ids[i] = alloc_int_1(d->elem_cmap_elem_counts[i], 0);
    ex_error = ex_get_elem_cmap(exoid, d->elem_cmap_ids[i], d->elem_cmap_elem_ids[i],
                                d->elem_cmap_side_ids[i], d->elem_cmap_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_get_elem_cmap %d", i);
  }

  return 0;
}

/* uni_dpi() -- setup distributed processing information for one processor
 *
 * When the problem is not partitioned, we need to set some acceptable
 * defaults.
 */
void uni_dpi(Dpi *dpi, Exo_DB *exo) {
  int i;
  int len;
  dpi->num_elems = exo->num_elems;

  if (exo->elem_elem_conn_exists) {
    dpi->elem_elem_list_global = exo->elem_elem_list;
  } else {
    dpi->elem_elem_list_global = NULL;
  }

  len = dpi->num_elems;
  dpi->elem_owner = alloc_int_1(len, ProcID);

  dpi->num_elems_global = exo->num_elems;

  dpi->num_elem_blocks_global = exo->num_elem_blocks;
  dpi->num_neighbors = 0;
  dpi->num_node_sets_global = exo->num_node_sets;
  dpi->num_side_sets_global = exo->num_side_sets;

  dpi->num_elems_global = exo->num_elems;
  dpi->num_internal_nodes = exo->num_nodes;
  dpi->num_boundary_nodes = 0;
  dpi->num_external_nodes = 0;
  dpi->num_owned_nodes = exo->num_nodes;
  dpi->num_universe_nodes = exo->num_nodes;
  dpi->num_nodes_global = exo->num_nodes;

  /*
   * Note! This aliasing of these pointers into the exo structure has
   * two advantages and one disadvantage.
   *
   *	(+) it's very easy to do
   *	(+) it makes more economic use of memory
   *	(-) it makes it too easy to free_dpi() and nuke your
   *        EXODUS information
   *
   * We'll just do it for now!
   */

  len = dpi->num_elems;
  dpi->elem_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < len; i++) {
    dpi->elem_index_global[i] = i;
  }

  len = dpi->num_universe_nodes;
  dpi->node_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < len; i++) {
    dpi->node_index_global[i] = i;
  }

  dpi->eb_id_global = exo->eb_id;

  dpi->eb_num_nodes_per_elem_global = exo->eb_num_nodes_per_elem;

  dpi->neighbor = alloc_int_1(1, ProcID);

  dpi->ns_id_global = exo->ns_id;

  dpi->set_membership[0] = -1;
  dpi->ptr_set_membership[0] = 0;
  dpi->ptr_set_membership[1] = 1;

  dpi->ss_id_global = exo->ss_id;

  dpi->ss_internal_global = find_ss_internal_boundary(exo);

  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
/* free_dpi() -- free internal dynamically allocated memory in a
 *               Dpi struct
 *
 * During rd_dpi(), various arrays are allocated. To cleanly
 * free up this  memory, this routine can be used. Typically,
 * if dpi is a (Dpi *), then
 *
 *	free_dpi(dpi);
 *	free(dpi);
 *
 * would accomplish the desired effect.
 *
 * Created: 1997/08/23 15:50 MDT pasacki@sandia.gov
 */

void free_dpi(Dpi *d) {
  EH(GOMA_ERROR, "Not implemented free_dpi");
  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
/* free_dpi_uni() -- free internal dynamically allocated memory in
 *                   Dpi SERIAL!
 *
 * During rd_dpi(), various arrays are allocated. To cleanly free up this
 * memory, this routine can be used. Typically, if dpi is a (Dpi *), then
 *
 *	free_dpi(dpi);
 *	free(dpi);
 *
 * would accomplish the desired effect.
 *
 * !!!! IMPORTANT - DANGEROUS PROGRAMMING NOTE !!!!!
 * For serial processing, not all of the pieces of Dpi have been
 * allocated and some are merely aliases to parts of the exodus ii
 * database. To avoid munging that data, just free up what was
 * allocated in uni_dpi.
 *
 * Created: 1997/09/11 10:43 MDT pasacki@sandia.gov
 */

void free_dpi_uni(Dpi *d) {
  safer_free((void **)&(d->elem_owner));
  safer_free((void **)&(d->elem_index_global));
  safer_free((void **)&(d->node_index_global));
  safer_free((void **)&(d->neighbor));
  safer_free((void **)&(d->ptr_set_membership));
  safer_free((void **)&(d->set_membership));
  safer_free((void **)&(d->ss_index_global));
  free(d->ss_internal_global);

  return;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/* init_dpi_struct() -- initialize some defaults
 *
 * This is meant to be called right after allocation, to help setup some
 * reasonable defaults to describe an empty data structure. Call it a poor
 * man's constructor.
 *
 * Created: 1999/08/11 17:04 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void init_dpi_struct(Dpi *d) {
  if (d == NULL) {
    EH(GOMA_ERROR, "Empty structure to initialize?");
  }
  memset((void *)d, 0, sizeof(Dpi));
  return;
}
/************************************************************************/
/************************************************************************/

/* exo_dpi_clone -- transfer needed global monolith data to child piece
 *
 *
 * Notes: Avert aliasing problems by allocating full-fledged arrays for dpi
 *        that won't disappear if the monolith does.
 *
 * Created: 1999/08/24 11:44 MDT pasacki@sandia.gov
 */

void exo_dpi_clone(Exo_DB *exo, Dpi *dpi) {
  int len;

  dpi->num_nodes_global = exo->num_nodes;
  dpi->num_elems_global = exo->num_elems;
  dpi->num_elem_blocks_global = exo->num_elem_blocks;
  dpi->num_node_sets_global = exo->num_node_sets;
  dpi->num_side_sets_global = exo->num_side_sets;

  /*
   * Allocate and fill arrays for element blocks...
   */

  len = dpi->num_elem_blocks_global * sizeof(int);

  dpi->eb_id_global = smalloc(len);
  memcpy(dpi->eb_id_global, exo->eb_id, len);

  /*
   * Allocate and fill arrays for node sets...
   */

  len = dpi->num_node_sets_global * sizeof(int);

  dpi->ns_id_global = smalloc(len);
  memcpy(dpi->ns_id_global, exo->ns_id, len);

  /*
   * Allocate and fill arrays for side sets...
   */
  len = dpi->num_side_sets_global * sizeof(int);

  dpi->ss_id_global = smalloc(len);
  memcpy(dpi->ss_id_global, exo->ss_id, len);

  return;
}
