/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2021 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* wr_dpi() -- write distributed processing information
 *
 * Notes:
 *	    [1] Use netCDF to write out this great new data.
 *
 *	    [2] This should nicely augment EXODUS II finite element data.
 *
 *	    [3] Try to use names that are identical to the names of the
 *              structure elements defined in "dpi.h"
 *
 *	    [4] Write out arrays in one shot instead of an element at
 *		a time.
 *
 *
 * Created: 1997/05/16 14:31 MDT pasacki@sandia.gov
 *
 * Revised: 1997/05/18 12:54 MDT pasacki@sandia.gov
 */

#define GOMA_WR_DPI_C

#include "rd_dpi.h"
#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#ifndef lint
#endif

/* #define NO_NETCDF_2		 for pure netCDF 3 -- still impure since
 * EXODUS II v3.00 still requires backwards compatibility
 */

#include "netcdf.h"
#include "wr_dpi.h"

/*
 * Sigh, if you need to run netCDF 2 then here's some definitions to tide
 * you over until netCDF 3 is working for EXODUS II...
 */

/*
 * I like these symbols as more lucid indicators of which netCDF we are
 * using...
 */

#ifndef NC_MAX_VAR_DIMS
#define NC_MAX_VAR_DIMS MAX_VAR_DIMS
#endif

#ifndef NC_MAX_NAME
#define NC_MAX_NAME MAX_NC_NAME
#endif

#ifndef NC_INT
#define NC_INT NC_LONG
#endif

/*
 * Might need some NO_NETCDF_2 definitions here ...
 */

#include "dpi.h"
#include "mm_eh.h"
#include "std.h"

/*
 * Prototypes of functions defined here, but needed elsewhere.
 */
#define CHECK_EX_ERROR(err, format, ...)                              \
  do {                                                                \
    if (err < 0) {                                                    \
      goma_eh(GOMA_ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__); \
    }                                                                 \
  } while (0)

int wr_dpi(Dpi *d, char *filename) {
  float version = -4.98; /* initialize. ex_open() changes this. */
  int comp_wordsize = sizeof(dbl);
  int io_wordsize = 0;
  int exoid = ex_open(filename, EX_WRITE, &comp_wordsize, &io_wordsize, &version);
  CHECK_EX_ERROR(exoid, "ex_open");

  GOMA_EH(one_dpi(d), "one_dpi");
  int ex_error;

  ex_error =
      ex_put_init_global(exoid, d->num_nodes_global, d->num_elems_global, d->num_elem_blocks_global,
                         d->num_node_sets_global, d->num_side_sets_global);
  CHECK_EX_ERROR(ex_error, "ex_put_init_global");

  ex_error = ex_put_ns_param_global(exoid, d->ns_id_global, d->num_ns_global_node_counts,
                                    d->num_ns_global_df_counts);
  CHECK_EX_ERROR(ex_error, "ex_put_ns_param_global");

  ex_error = ex_put_ss_param_global(exoid, d->ss_id_global, d->num_ss_global_side_counts,
                                    d->num_ss_global_df_counts);
  CHECK_EX_ERROR(ex_error, "ex_put_ss_param_global");

  ex_error = ex_put_eb_info_global(exoid, d->global_elem_block_ids, d->global_elem_block_counts);
  CHECK_EX_ERROR(ex_error, "ex_put_eb_info_global");

  ex_error = ex_put_init_info(exoid, d->num_proc, d->num_proc_in_file, &d->ftype);
  CHECK_EX_ERROR(ex_error, "ex_put_init_info");

  ex_error = ex_put_id_map(exoid, EX_NODE_MAP, d->node_index_global);
  CHECK_EX_ERROR(ex_error, "ex_put_id_map EX_NODE_MAP");

  ex_error = ex_put_id_map(exoid, EX_ELEM_MAP, d->elem_index_global);
  CHECK_EX_ERROR(ex_error, "ex_put_id_map EX_ELEM_MAP");

  ex_error = ex_put_loadbal_param(exoid, d->num_internal_nodes, d->num_boundary_nodes,
                                  d->num_external_nodes, d->num_internal_elems, d->num_border_elems,
                                  d->num_node_cmaps, d->num_elem_cmaps, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_loadbal_param");

  ex_error = ex_put_processor_elem_maps(exoid, d->proc_elem_internal, d->proc_elem_border, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_processor_elem_maps");

  ex_error = ex_put_processor_node_maps(exoid, d->proc_node_internal, d->proc_node_boundary,
                                        d->proc_node_external, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_processor_node_maps");

  ex_error = ex_put_cmap_params(exoid, d->node_cmap_ids, d->node_cmap_node_counts, d->elem_cmap_ids,
                                d->elem_cmap_elem_counts, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_cmap_params");

  for (int i = 0; i < d->num_node_cmaps; i++) {
    ex_error = ex_put_node_cmap(exoid, d->node_cmap_ids[i], d->node_map_node_ids[i],
                                d->node_map_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_put_node_cmap cmap %d", i);
  }
  for (int i = 0; i < d->num_elem_cmaps; i++) {
    ex_error = ex_put_elem_cmap(exoid, d->elem_cmap_ids[i], d->elem_cmap_elem_ids[i],
                                d->elem_cmap_side_ids[i], d->elem_cmap_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_put_elem_cmap cmap %d", i);
  }

  ex_error = ex_close(exoid);
  CHECK_EX_ERROR(ex_error, "ex_close");

  GOMA_EH(zero_dpi(d), "one_dpi");
  return 0;
}

#ifdef YOU_NEED_IT
static char *string_type(nc_type t) {
  static char t_int[] = "INT";
  static char t_dbl[] = "DOUBLE";
  static char t_chr[] = "CHAR";
  static char t_flt[] = "FLOAT";
  static char t_byt[] = "BYTE";
  static char t_unk[] = "UNKNOWN";

  switch (t) {
  case NC_INT:
    return (t_int);
    break;

  case NC_DOUBLE:
    return (t_dbl);
    break;

  case NC_CHAR:
    return (t_chr);
    break;

  case NC_FLOAT:
    return (t_flt);
    break;

  case NC_BYTE:
    return (t_byt);
    break;

  default:
    return (t_unk);
    break;
  }

  return (t_unk);
}
#endif
