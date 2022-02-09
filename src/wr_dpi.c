/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Sandia Corporation.                                  *
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

#include <exodusII.h>

#include "rd_dpi.h"
#include "rf_mp.h"

#ifdef STDC_HEADERS
#include <stdlib.h>
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

  ex_error = ex_put_loadbal_param(exoid, d->base_internal_nodes, d->base_boundary_nodes,
                                  d->base_external_nodes, d->base_internal_elems,
                                  d->base_border_elems, d->num_node_cmaps, 0, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_loadbal_param");

  ex_error =
      ex_put_cmap_params(exoid, d->node_cmap_ids, d->node_cmap_node_counts, NULL, NULL, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_cmap_params");

  ex_error = ex_put_processor_node_maps(exoid, d->proc_node_internal, d->proc_node_boundary,
                                        d->proc_node_external, ProcID);
  CHECK_EX_ERROR(ex_error, "ex_put_processor_node_maps");

  for (int i = 0; i < d->num_node_cmaps; i++) {
    ex_error = ex_put_node_cmap(exoid, d->node_cmap_ids[i], d->node_map_node_ids[i],
                                d->node_map_proc_ids[i], ProcID);
    CHECK_EX_ERROR(ex_error, "ex_put_node_cmap cmap %d", i);
  }

  ex_error = ex_close(exoid);
  CHECK_EX_ERROR(ex_error, "ex_close");

  GOMA_EH(zero_dpi(d), "zero_dpi");
  return 0;
}