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
 *
 * Notes:
 *	    [1] The information is in netCDF format.
 *
 *	    [2] Typically, this augments EXODUS II finite element data.
 *		
 *	    [3] Try to use netCDF names identical to the names of the
 *              structure elements defined in "dpi.h"
 *
 *	    [4] Read in arrays in one shot instead of an element at
 *		a time.
 *
 *	    [5] Allocate space as needed to hold the information. Use
 *		free_dpi() to release memory that was dynamically allocated
 *              here.
 *
 *	    [6] Two routes through code - netCDF 3, as originally intended
 *              and netCDF 2 as demanded by exigencies of working with
 *		EXODUS II
 *
 *	    [7] Assume basic skeleton space for the Dpi structure has been
 *              allocated. However, individual arrays, etc, will have space
 *              allocated for them in this routine.
 *
 *	    [8] In the future, a flag might indicate that "globally similar"
 *		information is to be read or not read. Such a capability
 *		would permit processor 0 to read the globally similar
 *		information. MPI could be used to broadcast it to every other
 *		processor and the read step could be trimmed to read only
 *		that information that is unique to that processor. Note
 *		that some globally similar information resides in the
 *		EXODUS II portion of the data as well. My impression is that
 *		startup and latency issues are more important than message
 *		length issues. If substantial savings in I/O time justifies
 *		the expenditure in programming complexity, then a conditional
 *		switch could be inserted to enable the reading of purelylocal
 *		or globally-similar data.
 *		
 *	    [9] Updated per the new Dpi format used by brk/fix.
 *
 * Created: 1997/07/10 15:39 MDT pasacki@sandia.gov
 *
 * Revised: 1997/07/21 09:45 MDT pasacki@sandia.gov
 *
 * Revised: 1998/12/16 14:21 MST pasacki@sandia.gov
 */

#define GOMA_RD_DPI_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

#ifndef lint
#endif

/* for pure netCDF 3 - but still not yet. */

#include "netcdf.h"
#include "exodusII.h"
#include "rd_mesh.h"

/*
 * Sigh, if you need to run netCDF 2 then here's some definitions to tide
 * you over until netCDF 3 is working for EXODUS II...
 *
 * I had trouble linking EXODUS II with netcdf 3 as of June 1997, in that
 * execution started complaining about nonexistent global variables(?)
 *
 * Hence, the netCDF 3 implementation is usually off for backward compatibility
 * with EXODUS II.
 *
 * Perhaps later it will work better.
 *
 * 1998/07/23 10:24 MDT pasacki@sandia.gov - "Maybe now that EXODUS II v3.00
 * is out the netCDF 3 stuff will work well! Nope! It leans heavily on
 * backward compatability mode. Thus we're hosed and must rely on backward
 * compatibility mode indefinitely."
 *
 */

/*
 * I like these symbols as more lucid indicators of which netCDF we are
 * using...
 */



#define NETCDF_3


#ifndef NC_MAX_VAR_DIMS
#define NC_MAX_VAR_DIMS		MAX_VAR_DIMS    
#endif

#ifndef NC_MAX_NAME
#define NC_MAX_NAME		MAX_NC_NAME
#endif

/*
 * Might need some NO_NETCDF_2 definitions here ...
 */

#include "std.h"
#include "rf_allo.h"
#include "mm_eh.h"
#include "exo_struct.h"
#include "dpi.h"
#include "rf_mp.h"		/* to know ProcID */
#include "rd_dpi.h"

static int get_variable_call_count = 0;

/*
 * Prototypes of functions defined here and needed only here.
 */

static void get_variable
(const int netcdf_unit,
       const nc_type netcdf_type,
       const int num_dimensions,
       const int dimension_val_1,
       const int dimension_val_2,
       const int variable_identifier,
       void *variable_address);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

int rd_dpi(Dpi *d, char *fn) {
  int err, len, status = 0, u;
  struct Shadow_Identifiers si;
  zeroStructures(&si, 1);

#ifdef NETCDF_3			/* only since those code chunks are the only
				 * ones using more extensive error reporting
				 * for now. If you need them, unprotect `em. */
#endif

  /*
   * From the C interface guide, the basic calling sequence is given
   * for reading dimensions and variables from a netCDF dataset for the
   * case
   *
   *       when the names of the dimensions and variables are known
   *
   * NETCDF 3:
   *
   * nc_open();
   * nc_inq_dimid();
   * nc_inq_varid();
   * nc_get_att();
   * nc_get_var();
   * nc_close();
   *
   * NETCDF 2:
   *
   * ncopen();
   * ncdimid();
   * ncvarid();
   * ncattget();
   * ncvarget();
   * ncclose();
   *
   */

  /*
   * 1. Open the file.
   */

#ifdef NETCDF_3
  err = nc_open(fn, NC_NOWRITE, &u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_open() problem.");
    }
#endif
#ifdef NETCDF_2
  err = ncopen(fn, NC_NOWRITE);
  EH(err, "ncopen() problem.");
  u   = err;
#endif

  /*
   *  Alternative method: set ncopts to zero to stop fatal exits from
   *  netcdf
   */
  //ncopts = 0;

  /*
   * 2. Get dimension identifiers.
   *
   * These are determined from their names which are defined in dpi.h. 
   *
   * These integer dimension IDs, if read properly from the open file, 
   * are stuck into the Shadow Identifiers si structure for later use.
   * Briefly, the "TRUE" and "FALSE" arguments you see refer to whether
   * it is critical that this particular dimension exist in the database.
   * Often, a FALSE flag is desirable for objects of zero length that
   * do not exist in this case.
   */

  getdid(u, DIM_LEN_EB_NUM_PRIVATE_ELEMS,     FALSE,
	 &si.len_eb_num_private_elems);
  getdid(u, DIM_LEN_ELEM_VAR_TAB_GLOBAL,      FALSE,
	 &si.len_elem_var_tab_global);
  getdid(u, DIM_LEN_ELEM_ELEM_LIST,	      TRUE,
	 &si.len_elem_elem_list);
  getdid(u, DIM_LEN_NODE_DESCRIPTION,         TRUE,
	 &si.len_node_description);

  getdid(u, DIM_LEN_NS_NODE_LIST,             FALSE,
	 &si.len_ns_node_list);
  getdid(u, DIM_LEN_NS_DISTFACT_LIST,         FALSE,
	 &si.len_ns_distfact_list);
  getdid(u, DIM_LEN_SS_ELEM_LIST,             FALSE,
	 &si.len_ss_elem_list);
  getdid(u, DIM_LEN_SS_DISTFACT_LIST,         FALSE,
	 &si.len_ss_distfact_list);
  getdid(u, DIM_LEN_STRING,                   TRUE,
	 &si.len_string);

  getdid(u, DIM_LEN_PTR_SET_MEMBERSHIP,       TRUE,
	 &si.len_ptr_set_membership);
  getdid(u, DIM_LEN_SET_MEMBERSHIP,           TRUE,
	 &si.len_set_membership);
  getdid(u, DIM_NUM_ELEM_BLOCKS,              FALSE,
	 &si.num_elem_blocks);
  getdid(u, DIM_NUM_ELEM_BLOCKS_GLOBAL,       FALSE,
	 &si.num_elem_blocks_global);

  getdid(u, DIM_NUM_ELEMS,		      TRUE,
	 &si.num_elems);

  getdid(u, DIM_NUM_GLOBAL_NODE_DESCRIPTIONS, TRUE,
	 &si.num_global_node_descriptions);
  getdid(u, DIM_NUM_NEIGHBORS,                FALSE,
	 &si.num_neighbors);
  getdid(u, DIM_NUM_NODE_SETS,                FALSE,
	 &si.num_node_sets);
  getdid(u, DIM_NUM_NODE_SETS_GLOBAL,         FALSE,
	 &si.num_node_sets_global);

  getdid(u, DIM_NUM_NODES,                    TRUE,
	 &si.num_nodes);

  getdid(u, DIM_NUM_PROPS_EB,                FALSE,
	 &si.num_props_eb);
  getdid(u, DIM_NUM_PROPS_NS,                FALSE,
	 &si.num_props_ns);
  getdid(u, DIM_NUM_PROPS_SS,                FALSE,
	 &si.num_props_ss);

  getdid(u, DIM_NUM_SIDE_SETS,                FALSE,
	 &si.num_side_sets);
  getdid(u, DIM_NUM_SIDE_SETS_GLOBAL,         FALSE,
	 &si.num_side_sets_global);
  getdid(u, DIM_NUM_UNIVERSE_NODES,           TRUE,
	 &si.num_universe_nodes);

  getdid(u, DIM_LEN_SS_BLOCK_INDEX_GLOBAL, FALSE,
         &si.len_ss_block_list_global);

  getdid(u, DIM_LEN_SS_BLOCK_LIST_GLOBAL, FALSE,
         &si.len_ss_block_list_global);

  /*
   * 3. Using the dimension IDs, inquire of the dimension values. Load
   *    those values into the DP array.
   *
   *    These dimension values are important so that we know how much
   *    space to allocate to hold array variables below...
   *
   * getdim(netcdf_unit, dimension_id, address of the answer);
   */
  getdim(u, si.len_eb_num_private_elems,     &d->len_eb_num_private_elems);
  getdim(u, si.len_elem_var_tab_global,      &d->len_elem_var_tab_global);
  getdim(u, si.len_elem_elem_list,           &d->len_elem_elem_list);
  getdim(u, si.len_node_description,         &d->len_node_description);
  getdim(u, si.len_ns_node_list,             &d->len_ns_node_list);
  getdim(u, si.len_ns_distfact_list,         &d->len_ns_distfact_list);
  getdim(u, si.len_ss_elem_list,             &d->len_ss_elem_list);
  getdim(u, si.len_ss_distfact_list,         &d->len_ss_distfact_list);
  getdim(u, si.len_string,		     &d->len_string);

  getdim(u, si.len_ptr_set_membership,       &d->len_ptr_set_membership);
  getdim(u, si.len_set_membership,           &d->len_set_membership);
  getdim(u, si.num_elem_blocks,              &d->num_elem_blocks);
  getdim(u, si.num_elem_blocks_global,	     &d->num_elem_blocks_global);

  getdim(u, si.num_elems,		     &d->num_elems);

  getdim(u, si.num_global_node_descriptions, &d->num_global_node_descriptions);
  getdim(u, si.num_neighbors,                &d->num_neighbors);
  getdim(u, si.num_node_sets,                &d->num_node_sets);
  getdim(u, si.num_node_sets_global,         &d->num_node_sets_global);

  getdim(u, si.num_nodes,		     &d->num_nodes);

  getdim(u, si.num_props_eb,                 &d->num_props_eb);
  getdim(u, si.num_props_ns,                 &d->num_props_ns);
  getdim(u, si.num_props_ss,                 &d->num_props_ss);


  getdim(u, si.num_side_sets,                &d->num_side_sets);
  getdim(u, si.num_side_sets_global,         &d->num_side_sets_global);
  getdim(u, si.num_universe_nodes,           &d->num_universe_nodes);

  // I don't really want a variable in dpi struct for this but from what I see
  // netcdf expects dimensions to be defined in the netcdf file
  int len_ss_block_list_global = 0;
  getdim(u, si.len_ss_block_list_global, &len_ss_block_list_global);


  /*
   * 4. Get variable identifiers from netCDF.
   *
   *    The Booleans are set only roughly based on what "ought" to exist.
   *    If you deem looser or tighter criteria, then by all means change
   *    them to suit your needs.
   */

  getvid(u, VAR_EB_ELEM_TYPE_GLOBAL,       FALSE, /* New! */
	 &si.eb_elem_type_global);
  getvid(u, VAR_EB_ID_GLOBAL,              FALSE,
	 &si.eb_id_global);
  getvid(u, VAR_EB_INDEX_GLOBAL,           FALSE,
	 &si.eb_index_global);
  getvid(u, VAR_EB_NUM_ATTR_GLOBAL,        FALSE, /* New! */
	 &si.eb_num_attr_global);
  getvid(u, VAR_EB_NUM_ELEMS_GLOBAL,       FALSE,
	 &si.eb_num_elems_global);
  getvid(u, VAR_EB_NUM_NODES_PER_ELEM_GLOBAL, FALSE, /* New! */
	 &si.eb_num_nodes_per_elem_global);
  getvid(u, VAR_EB_NUM_PRIVATE_ELEMS,      FALSE,
	 &si.eb_num_private_elems);
  getvid(u, VAR_EB_PROP_GLOBAL,		   FALSE,
	 &si.eb_prop_global);

  getvid(u, VAR_ELEM_INDEX_GLOBAL,	   TRUE, /* New! */
	 &si.elem_index_global);

  getvid(u, VAR_ELEM_OWNER,		   TRUE, /* New! */
	 &si.elem_owner);
  getvid(u, VAR_ELEM_ELEM_LIST_GLOBAL,	   TRUE, /* New! */
	 &si.elem_elem_list_global);
  getvid(u, VAR_ELEM_ELEM_TWST_GLOBAL,	   TRUE, /* New! */
	 &si.elem_elem_twst_global);
  getvid(u, VAR_ELEM_ELEM_FACE_GLOBAL,	   TRUE, /* New! */
	 &si.elem_elem_face_global);
  getvid(u, VAR_ELEM_ELEM_PROC_GLOBAL,	   TRUE, /* New! */
	 &si.elem_elem_proc_global);

  getvid(u, VAR_ELEM_VAR_TAB_GLOBAL,	   FALSE, /* New! */
	 &si.elem_var_tab_global);
  getvid(u, VAR_GLOBAL_NODE_DESCRIPTION,   TRUE,
	 &si.global_node_description);
  getvid(u, VAR_MY_NAME,                   TRUE,
	 &si.my_name);
  getvid(u, VAR_NEIGHBOR,                  FALSE,
	 &si.neighbor);
  
  getvid(u, VAR_NODE_INDEX_GLOBAL,	   TRUE, /* New! */
	 &si.node_index_global);

  getvid(u, VAR_NS_DISTFACT_INDEX_GLOBAL,  FALSE, /* New! */
	 &si.ns_distfact_index_global);
  getvid(u, VAR_NS_DISTFACT_LEN_GLOBAL,  FALSE, /* New! */
	 &si.ns_distfact_len_global);
  getvid(u, VAR_NS_DISTFACT_LIST_INDEX_GLOBAL,  FALSE, /* New! */
	 &si.ns_distfact_list_index_global);

  getvid(u, VAR_NS_ID_GLOBAL,              FALSE,
	 &si.ns_id_global);
  getvid(u, VAR_NS_INDEX_GLOBAL,           FALSE,
	 &si.ns_index_global);

  getvid(u, VAR_NS_NODE_INDEX_GLOBAL,      FALSE, /* New! */
	 &si.ns_node_index_global);
  getvid(u, VAR_NS_NODE_LEN_GLOBAL,        FALSE, /* New! */
	 &si.ns_node_len_global);
  getvid(u, VAR_NS_NODE_LIST_INDEX_GLOBAL,  FALSE, /* New! */
	 &si.ns_node_list_index_global);


  getvid(u, VAR_NS_NUM_DISTFACTS_GLOBAL,   FALSE,
	 &si.ns_num_distfacts_global);
  getvid(u, VAR_NS_NUM_NODES_GLOBAL,       FALSE,
	 &si.ns_num_nodes_global);

  getvid(u, VAR_NS_PROP_GLOBAL,       FALSE,
	 &si.ns_prop_global);

  getvid(u, VAR_NUM_BOUNDARY_NODES,        TRUE,
	 &si.num_boundary_nodes);
  getvid(u, VAR_NUM_DOFS_GLOBAL,           TRUE,
	 &si.num_dofs_global);
  getvid(u, VAR_NUM_ELEMS_GLOBAL,          TRUE,
	 &si.num_elems_global);
  getvid(u, VAR_NUM_EXTERNAL_NODES,        TRUE,
	 &si.num_external_nodes);
  getvid(u, VAR_NUM_INTERNAL_NODES,        TRUE,
	 &si.num_internal_nodes);
  getvid(u, VAR_NUM_NODES_GLOBAL,          TRUE,
	 &si.num_nodes_global);
  getvid(u, VAR_PTR_SET_MEMBERSHIP,        TRUE,
	 &si.ptr_set_membership);
  getvid(u, VAR_SET_MEMBERSHIP,            TRUE,
	 &si.set_membership);

  getvid(u, VAR_SS_DISTFACT_INDEX_GLOBAL,    FALSE, /* New! */
	 &si.ss_distfact_index_global);

  getvid(u, VAR_SS_DISTFACT_LEN_GLOBAL,    FALSE, /* New! */
	 &si.ss_distfact_len_global);
  getvid(u, VAR_SS_DISTFACT_LIST_INDEX_GLOBAL,    FALSE, /* New! */
	 &si.ss_distfact_list_index_global);

  getvid(u, VAR_SS_ELEM_INDEX_GLOBAL,    FALSE, /* New! */
	 &si.ss_elem_index_global);

  getvid(u, VAR_SS_ELEM_LEN_GLOBAL,    FALSE, /* New! */
	 &si.ss_elem_len_global);
  getvid(u, VAR_SS_ELEM_LIST_INDEX_GLOBAL,    FALSE, /* New! */
	 &si.ss_elem_list_index_global);

  getvid(u, VAR_SS_NODE_LEN_GLOBAL,    FALSE, /* New! */
	 &si.ss_node_len_global);

  getvid(u, VAR_SS_ID_GLOBAL,              FALSE,
	 &si.ss_id_global);
  getvid(u, VAR_SS_INDEX_GLOBAL,           FALSE,
	 &si.ss_index_global);
  getvid(u, VAR_SS_NUM_DISTFACTS_GLOBAL,   FALSE,
	 &si.ss_num_distfacts_global);
  getvid(u, VAR_SS_NUM_SIDES_GLOBAL,       FALSE,
	 &si.ss_num_sides_global);

  getvid(u, VAR_SS_PROP_GLOBAL,       FALSE,
	 &si.ss_prop_global);

  getvid(u, VAR_UNDEFINED_BASIC_EQNVAR_ID, TRUE,
	 &si.undefined_basic_eqnvar_id);

  getvid(u, VAR_SS_INTERNAL_GLOBAL, FALSE,
         &si.ss_internal_global);

  getvid(u, VAR_SS_BLOCK_INDEX_GLOBAL, FALSE,
         &si.ss_block_index_global);

  getvid(u, VAR_SS_BLOCK_LIST_GLOBAL, FALSE,
         &si.ss_block_list_global);

  /*
   * 5. Allocate space to hold array variables.
   *
   *    a. Should verify that these dimensions are reasonable
   *       numbers 0 < dimval < humongous
   *
   *    b. Simple scalar variables are not allocated here - they got
   *       their space when *d became meaningful.
   */

  /*
   * HKM -> This is over the number of elements. Can we
   *        eliminate the need for these large arrays ?
   */
  len = d->num_elems;
  if ( len > 0 )
    {
      d->elem_index_global = alloc_int_1(len, INT_NOINIT);
      d->elem_owner        = alloc_int_1(len, INT_NOINIT);
    }

  len = d->len_elem_elem_list;
  if ( len > 0 )
    {
      d->elem_elem_list_global = alloc_int_1(len, INT_NOINIT);
      d->elem_elem_twst_global = alloc_int_1(len, INT_NOINIT);
      d->elem_elem_face_global = alloc_int_1(len, INT_NOINIT);
      d->elem_elem_proc_global = alloc_int_1(len, INT_NOINIT);
    }

  len = d->num_elem_blocks_global;
  if ( len > 0 )
    {
      d->eb_id_global                 = alloc_int_1(len, INT_NOINIT);
      d->eb_num_attr_global           = alloc_int_1(len, INT_NOINIT);
      d->eb_num_elems_global          = alloc_int_1(len, INT_NOINIT);
      d->eb_num_nodes_per_elem_global = alloc_int_1(len, INT_NOINIT);
      d->eb_elem_type_global =
	  alloc_VecFixedStrings(len, (MAX_STR_LENGTH+1));
    }

  len = d->num_props_eb;
  if (len > 0)		/* Props are special ID is implicit. */
    {
      d->eb_prop_global = alloc_int_2(len, d->num_elem_blocks_global, -1);
    }

  if (d->num_props_ns > 1)	/* Props are special, ID is implicit. */
    {
      d->ns_prop_global =
	  alloc_int_2(d->num_props_ns, d->num_node_sets_global, -1);
    }

  if (d->num_props_ss > 1)	/* Props are special, ID is implicit. */
    {
      d->ss_prop_global =
	  alloc_int_2(d->num_props_ss, d->num_side_sets_global, -1);
    }

  len = d->num_elem_blocks;
  if ( len > 0 )
    {
      d->eb_index_global      = alloc_int_1(len, INT_NOINIT);
      d->eb_num_private_elems = alloc_int_1(len, INT_NOINIT);
    }

  len = d->len_elem_var_tab_global;
  if ( len > 0 )
    {
      d->elem_var_tab_global = alloc_int_1(len, 0);
    }

  if (d->num_global_node_descriptions > 0)
    {
      d->global_node_description =
	  alloc_int_2(d->num_global_node_descriptions,
		      d->len_node_description, INT_NOINIT);
      //(int **) smalloc(d->num_global_node_descriptions*sizeof(int *));
    }

  len = d->num_nodes;
  if ( len > 0 )
    {
      d->node_index_global = alloc_int_1(len, INT_NOINIT);
    }

  if ( d->num_neighbors > 0 )
    {
      d->neighbor =  alloc_int_1(d->num_neighbors, INT_NOINIT);
    }


  len = d->num_node_sets_global;
  if ( len > 0 )
    {
      d->ns_id_global             = alloc_int_1(len, INT_NOINIT);
      d->ns_num_distfacts_global  = alloc_int_1(len, INT_NOINIT);
      d->ns_num_nodes_global      = alloc_int_1(len, INT_NOINIT);
      d->ns_node_index_global     = alloc_int_1(len, INT_NOINIT);
      d->ns_distfact_index_global = alloc_int_1(len, INT_NOINIT);
    }

  if (d->num_node_sets > 0)
    {
      d->ns_index_global =  alloc_int_1(d->num_node_sets, INT_NOINIT);
    }

  if ( d->len_ns_node_list > 0 )
    {
      d->ns_node_list_index_global = 
	  alloc_int_1(d->len_ns_node_list, INT_NOINIT);
    }

  if ( d->len_ns_distfact_list > 0 )
    {
      d->ns_distfact_list_index_global =
	  alloc_int_1(d->len_ns_distfact_list, INT_NOINIT);
    }

  if ( d->len_ptr_set_membership > 0 )
    {
      d->ptr_set_membership =
	  alloc_int_1(d->len_ptr_set_membership, INT_NOINIT);
    }

  if ( d->len_set_membership > 0 )
    {
      d->set_membership =
	  alloc_int_1(d->len_set_membership, INT_NOINIT);
    }

  len = d->num_side_sets_global;
  if ( len > 0 )
    {
      d->ss_distfact_index_global = alloc_int_1(len, INT_NOINIT);
      d->ss_elem_index_global     = alloc_int_1(len, INT_NOINIT);
      d->ss_id_global             = alloc_int_1(len, INT_NOINIT);
      d->ss_num_distfacts_global  = alloc_int_1(len, INT_NOINIT);
      d->ss_num_sides_global      = alloc_int_1(len, INT_NOINIT);
      d->ss_internal_global       = alloc_int_1(len, INT_NOINIT);
      d->ss_block_index_global    = alloc_int_1(len+1, INT_NOINIT);
    }

  if ( d->len_ss_elem_list > 0 )
    {
      d->ss_elem_list_index_global = 
	  alloc_int_1(d->len_ss_elem_list, INT_NOINIT);
    }

  if (d->len_ss_distfact_list > 0)
    {
      d->ss_distfact_list_index_global = 
	 alloc_int_1(d->len_ss_distfact_list, INT_NOINIT);
    }

  if (d->num_side_sets > 0)
    {
      d->ss_index_global = alloc_int_1(d->num_side_sets, INT_NOINIT);
    }

  if (len_ss_block_list_global > 0) {
    d->ss_block_list_global = alloc_int_1(len_ss_block_list_global, INT_NOINIT);
  }

  /*
   * 6. Get variables - this is messy. Netcdf 3 wants a different routine
   *    for each type. NetCDF 2 wants to know all the dimensions!
   */

  get_variable(u, NC_CHAR, 2,
               d->num_elem_blocks_global, d->len_string,
               si.eb_elem_type_global, &(d->eb_elem_type_global[0][0]));

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_id_global,			d->eb_id_global);

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks,		-1, 
	       si.eb_index_global,		d->eb_index_global);

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_attr_global,		d->eb_num_attr_global);

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_elems_global,		d->eb_num_elems_global);

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_nodes_per_elem_global,	
	       d->eb_num_nodes_per_elem_global);

  get_variable(u, NC_INT, 1, 
	       d->num_elem_blocks,		-1, 
	       si.eb_num_private_elems,		d->eb_num_private_elems);

  if ( d->num_props_eb > 1 )
    {
      get_variable(u, NC_INT, 2,
		   d->num_props_eb,		d->num_elem_blocks_global,
		   si.eb_prop_global,		&(d->eb_prop_global[0][0]));
    }

  if (d->len_elem_var_tab_global > 0)
    {
      get_variable(u, NC_INT, 1,
                   d->len_elem_var_tab_global, -1,
                   si.elem_var_tab_global, d->elem_var_tab_global);
    }

  if ( d->num_elems > 0 )
    {
      get_variable(u, NC_INT, 1, 
		   d->num_elems,	 -1, 
		   si.elem_index_global, d->elem_index_global);

      get_variable(u, NC_INT, 1, 
		   d->num_elems,	 -1, 
		   si.elem_owner, d->elem_owner);
    }

  get_variable(u, NC_INT, 2, 
	       d->num_global_node_descriptions,	d->len_node_description, 
	       si.global_node_description,&(d->global_node_description[0][0]));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.my_name,		&(d->my_name));

  get_variable(u, NC_INT, 1, 
	       d->num_neighbors,	-1, 
	       si.neighbor,		d->neighbor);

  if ( d->num_nodes > 0 )
    {
      get_variable(u, NC_INT, 1, 
		   d->num_nodes,	 -1, 
		   si.node_index_global, d->node_index_global);
    }

  get_variable(u, NC_INT, 0, 
	       -1, 	-1, 
	       si.ns_distfact_len_global, &(d->ns_distfact_len_global));

  get_variable(u, NC_INT, 0, 
	       -1, 	-1, 
	       si.ns_node_len_global, &(d->ns_node_len_global));

  get_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,	-1, 
	       si.ns_id_global,		d->ns_id_global);

  get_variable(u, NC_INT, 1, 
	       d->num_node_sets,	-1, 
	       si.ns_index_global,	d->ns_index_global);


  get_variable(u, NC_INT, 1, 
	       d->len_ns_distfact_list,	-1, 
	       si.ns_distfact_list_index_global, 
	       d->ns_distfact_list_index_global);

  get_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_distfact_index_global,	d->ns_distfact_index_global);

  get_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_node_index_global,	d->ns_node_index_global);

  get_variable(u, NC_INT, 1, 
	       d->len_ns_node_list,	-1, 
	       si.ns_node_list_index_global, 
	       d->ns_node_list_index_global);



  get_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_num_distfacts_global,	d->ns_num_distfacts_global);

  get_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,	-1, 
	       si.ns_num_nodes_global,	d->ns_num_nodes_global);

  if ( d->num_props_ns > 1 )
    {
      get_variable(u, NC_INT, 2,
		   d->num_props_ns,		d->num_node_sets_global,
		   si.ns_prop_global,		&(d->ns_prop_global[0][0]));
    }

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_boundary_nodes,	&(d->num_boundary_nodes));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_dofs_global,	&(d->num_dofs_global));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_elems_global,	&(d->num_elems_global));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_external_nodes,	&(d->num_external_nodes));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_internal_nodes,	&(d->num_internal_nodes));

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_nodes_global,	&(d->num_nodes_global));

  get_variable(u, NC_INT, 1, 
	       d->len_ptr_set_membership,	-1, 
	       si.ptr_set_membership,	d->ptr_set_membership);

  get_variable(u, NC_INT, 1, 
	       d->len_set_membership,	-1, 
	       si.set_membership,	d->set_membership);

  get_variable(u, NC_INT, 1, d->num_side_sets_global, -1, 
	       si.ss_distfact_index_global, d->ss_distfact_index_global);

  get_variable(u, NC_INT, 0, -1, -1, 
	       si.ss_distfact_len_global, &(d->ss_distfact_len_global));

  get_variable(u, NC_INT, 1, d->len_ss_distfact_list, -1, 
	       si.ss_distfact_list_index_global, 
	       d->ss_distfact_list_index_global);

  get_variable(u, NC_INT, 1, d->num_side_sets_global, -1, 
	       si.ss_elem_index_global, d->ss_elem_index_global);

  get_variable(u, NC_INT, 0, -1, -1, 
	       si.ss_elem_len_global, &(d->ss_elem_len_global));

  get_variable(u, NC_INT, 1, d->len_ss_elem_list, -1, 
	       si.ss_elem_list_index_global, d->ss_elem_list_index_global);

  get_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,	-1, 
	       si.ss_id_global,		d->ss_id_global);

  get_variable(u, NC_INT, 1, 
	       d->num_side_sets,	-1, 
	       si.ss_index_global,	d->ss_index_global);

  get_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,		-1, 
	       si.ss_num_distfacts_global,	d->ss_num_distfacts_global);

  get_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,		-1, 
	       si.ss_num_sides_global,		d->ss_num_sides_global);

  if ( d->len_elem_elem_list > 0 )
    {
      get_variable(u, NC_INT, 1, 
		   d->len_elem_elem_list,	-1, 
		   si.elem_elem_list_global,	d->elem_elem_list_global);
      
      get_variable(u, NC_INT, 1, 
		   d->len_elem_elem_list,	-1, 
		   si.elem_elem_twst_global,	d->elem_elem_twst_global);
      
      get_variable(u, NC_INT, 1, 
		   d->len_elem_elem_list,	-1, 
		   si.elem_elem_face_global,	d->elem_elem_face_global);
      
      get_variable(u, NC_INT, 1, 
		   d->len_elem_elem_list,	-1, 
		   si.elem_elem_proc_global,	d->elem_elem_proc_global);
    }

  if ( d->num_props_ss > 1 )
    {
      get_variable(u, NC_INT, 2,
		   d->num_props_ss,		d->num_side_sets_global,
		   si.ss_prop_global,		&(d->ss_prop_global[0][0]));
    }

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.undefined_basic_eqnvar_id, &(d->undefined_basic_eqnvar_id));

  if (d->num_side_sets_global > 0) {
      get_variable(u, NC_INT, 1,
                   d->num_side_sets_global, -1,
                   si.ss_internal_global, d->ss_internal_global);

      get_variable(u, NC_INT, 1,
                   d->num_side_sets_global + 1, -1,
                   si.ss_block_index_global, d->ss_block_index_global);

      get_variable(u, NC_INT, 1,
                   len_ss_block_list_global, -1,
                   si.ss_block_list_global, d->ss_block_list_global);
  }
  /*
   * 7. Close up.
   */

#ifdef NETCDF_3
  err = nc_close(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_close() problem.");
    }
#endif
#ifdef NETCDF_2
  err = ncclose(u);
  EH(err, "ncclose()");
#endif

  /*
   * 8. Calculate a few utility combinations
   */

  d->num_owned_nodes = d->num_internal_nodes + d->num_boundary_nodes;

  return(status);
}


/*
 * getdid() - get netCDF dimension identifiers
 *
 * Given the open unit and the string name of the dimension, load into the
 * input address the dimension identifier.
 *
 * Note: Intended to be netCDF 2 and netCDF 3 bi-compliant.
 *
 * Problem - some named dimensions might have zero length and not be recorded.
 * Need to account for this.
 *
 * 1997/08/23 09:24 MDT pasacki@sandia.gov
 * 1998/07/23 13:03 MDT pasacki@sandia.gov
 */

void
getdid(int netcdf_unit,			/* should already be open	(in) */
       char *string_name,		/* eg, "num_side_sets"          (in) */
       int hard_error_interpretation,	/* Boolean for ghosty dims      (in) */
       int *dimension_identifier_address) /*				(out)*/
{
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];  

#ifdef NETCDF_3  
  err  = nc_inq_dimid(netcdf_unit, string_name, dimension_identifier_address);

  /* Handle ghosty dims later */
  if (err != NC_NOERR)
    {
      *dimension_identifier_address = -1;
    }

  /*
   * Assume an error means this quanitity is not found here. That's OK for
   * some things, but not others.
   */
  if ( err != NC_NOERR && hard_error_interpretation )
    {
      sprintf(err_msg, "nc_inq_dimid() on %s id=%d", 
	      string_name, 
	      *dimension_identifier_address);
      EH(-1, err_msg);
    }
#endif

#ifdef NETCDF_2
  err  = ncdimid(netcdf_unit, string_name);
  if ( err == -1 && hard_error_interpretation )
    {

      sprintf(err_msg, "ncdimid() on %s rtn %d", string_name, err);
      EH(err, err_msg);
    }
  *dimension_identifier_address = err;
#endif  


  return;
}

/*
 * getvid() - get netCDF variable identifiers
 *
 * Given the open unit and the string name of the variable, load into the
 * input address the variable identifier.
 *
 * The Boolean tells whether this particular variable better be there or
 * if it is OK for the variable not to exist in the file. Nonexistence is
 * permissible under some circumstances. Some EXODUS II files might well
 * not have any sidesets or any nodesets, for example.
 *
 * Note: Intended to be netCDF 2 and netCDF 3 bi-compliant.
 *
 * 1997/08/23 09:25 MDT pasacki@sandia.gov
 *
 * Revised: 1998/07/27 10:44 MDT pasacki@sandia.gov
 */

void
getvid(int netcdf_unit,		         /* open netCDF unit identifier (in) */
       char *string_name,	         /* Name of the variable        (in) */
       int hard_error_interpretation,	 /* Boolean for ghosty vars     (in) */
       int *variable_identifier_address) /* integer id of variable     (out) */
{
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];  

#ifdef NETCDF_3  
  err  = nc_inq_varid(netcdf_unit, string_name, variable_identifier_address);
  if ( err != NC_NOERR && hard_error_interpretation )
    {
      sprintf(err_msg, "nc_inq_varid() on %s id=%d", 
		   string_name, 
		   *variable_identifier_address);
      EH(-1, err_msg);
    }
#endif

#ifdef NETCDF_2
  err  = ncvarid(netcdf_unit, string_name);
  if ( err == -1 && hard_error_interpretation )
    {
      sprintf(err_msg, "ncvarid() on %s rtn %d", string_name, err);
      EH(err, err_msg);
    }
  *variable_identifier_address = err;
#endif  


  return;
}

/*
 * getdim() - get netCDF dimension
 *
 * Given the open unit and the dimension ID of the dimension, load into the
 * input address the value of the dimension.
 *
 * Note: Intended to be netCDF 2 and netCDF 3 bi-compliant.
 *
 * 1997/08/23 09:26 MDT pasacki@sandia.gov
 */

void
getdim(int netcdf_unit,
       int dimension_id,
       int *where)
{
  int err;
#ifdef NETCDF_3
  size_t swhere;
#endif
#ifdef NETCDF_2
  long swhere;
  char junk[MAX_CHAR_ERR_MSG];
#endif
  char err_msg[MAX_CHAR_ERR_MSG];  
  
  /*
   * If the earlier attempt to find the dimension ID failed due to
   * nonexistence of a noncritical chunk, then we need a means of gracefully
   * getting the heck out of here and setting things most appropriately.
   * Here, let's guess that the value of a nonexistent dimension is zero.
   */

  if ( dimension_id == -1 )
    {
      *where = 0;
    }
  else
    {
#ifdef NETCDF_3  
      err  = nc_inq_dimlen(netcdf_unit, dimension_id, &swhere);
      if ( err != NC_NOERR )
	{
	  sprintf(err_msg, "nc_inq_dimlen() on did=%d", dimension_id);
	  EH(-1, err_msg);
	}
      *where = (int) swhere;
#endif

#ifdef NETCDF_2
      err  = ncdiminq(netcdf_unit, dimension_id, junk, &swhere);
      sprintf(err_msg, "ncdiminq() on did %d rtns %d", dimension_id, 
	      err);
      EH(err, err_msg);
      *where = (int) swhere;
#endif  
    }

  return;
}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/* uni_dpi() -- setup distributed processing information for one processor
 *
 * When the problem is not partitioned, we need to set some acceptable
 * defaults.
 *
 * Created: 1997/08/27 16:16 MDT pasacki@sandia.gov
 *
 * Revised: 1997/08/28 10:05 MDT pasacki@sandia.gov
 */

void
uni_dpi(Dpi *dpi, 
	Exo_DB *exo)
{
  int i;
  int len;

  dpi->num_elems                     = exo->num_elems;
  
  /*
   * These defaults for the element decomposition information might
   * be rather uninformative and not useful. However, to usefully reflect
   * information in exo->, we'd need to already have some of the elem_elem
   * connectivity information ready right now.
   */

  if (exo->elem_elem_conn_exists)
    {
      dpi->len_elem_elem_list            = exo->elem_elem_pntr[exo->num_elems];
      dpi->elem_elem_list_global         = exo->elem_elem_list;
      dpi->elem_elem_twst_global         = exo->elem_elem_twst;
      dpi->elem_elem_face_global         = exo->elem_elem_face;
    }
  else
    {
      dpi->len_elem_elem_list            = 0;
      dpi->elem_elem_list_global         = NULL;
      dpi->elem_elem_twst_global         = NULL;
      dpi->elem_elem_face_global         = NULL;
    }

  len = dpi->num_elems;
  dpi->elem_owner = alloc_int_1(len, ProcID);


  dpi->num_elems_global              = exo->num_elems;

  /*
   * Set default values for dimension identifiers.
   */

  dpi->len_eb_num_private_elems      = exo->num_elem_blocks;
  dpi->len_node_description          = 1;

  /*
   * You need an upper and lower bound to point to even 1 setmembership entry.
   */

  dpi->len_ptr_set_membership        = 2;
  dpi->len_set_membership            = 1;

  dpi->num_elem_blocks               = exo->num_elem_blocks;
  dpi->num_elem_blocks_global        = exo->num_elem_blocks;
  dpi->num_global_node_descriptions  = 1;
  dpi->num_neighbors                 = 0;
  dpi->num_node_sets                 = exo->num_node_sets;
  dpi->num_node_sets_global          = exo->num_node_sets;
  dpi->num_side_sets                 = exo->num_side_sets;
  dpi->num_side_sets_global          = exo->num_side_sets;

  dpi->num_elems_global              = exo->num_elems;
  dpi->num_internal_nodes	     = exo->num_nodes;
  dpi->num_boundary_nodes            = 0;
  dpi->num_external_nodes            = 0;
  dpi->num_owned_nodes               = exo->num_nodes;
  dpi->num_universe_nodes            = exo->num_nodes;
  dpi->num_nodes_global              = exo->num_nodes;

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
  for ( i=0; i<len; i++)
    {
      dpi->elem_index_global[i] = i;
    }

  len = dpi->num_universe_nodes;
  dpi->node_index_global = alloc_int_1(len, INT_NOINIT); 
  for ( i=0; i<len; i++)
    {
      dpi->node_index_global[i] = i;
    }

  dpi->eb_id_global = exo->eb_id;

  len = dpi->num_elem_blocks;
  dpi->eb_index_global = alloc_int_1(len, INT_NOINIT);
  for (i = 0; i < len; i++) {
    dpi->eb_index_global[i] = i;
  }


  dpi->eb_num_nodes_per_elem_global = exo->eb_num_nodes_per_elem;

  dpi->eb_num_elems_global  = exo->eb_num_elems;

  dpi->eb_num_private_elems = exo->eb_num_elems;

  len = dpi->num_global_node_descriptions;
  if (len > 0) {
    dpi->global_node_description = 
	alloc_int_2(len, dpi->len_node_description, INT_NOINIT);
  }

  dpi->neighbor                = alloc_int_1(1, ProcID);

  dpi->ns_id_global            = exo->ns_id;

  if (dpi->num_node_sets > 0) {
    dpi->ns_index_global = alloc_int_1(dpi->num_node_sets, INT_NOINIT);
    for (i = 0; i < dpi->num_node_sets; i++) {
      dpi->ns_index_global[i] = i;
    }
  }

  dpi->ns_num_distfacts_global = exo->ns_num_distfacts;

  dpi->ns_num_nodes_global     = exo->ns_num_nodes;

  if ( dpi->len_ptr_set_membership > 0 ) {
    dpi->ptr_set_membership = 
	alloc_int_1(dpi->len_ptr_set_membership, INT_NOINIT);
  }

  if ( dpi->len_set_membership > 0 ) {
    dpi->set_membership =
	alloc_int_1(dpi->len_set_membership, INT_NOINIT);
  }

  dpi->set_membership[0]       = -1;
  dpi->ptr_set_membership[0]   = 0;
  dpi->ptr_set_membership[1]   = 1;

  dpi->ss_id_global            = exo->ss_id;

  if (dpi->num_side_sets > 0) {
    dpi->ss_index_global = alloc_int_1(dpi->num_side_sets, INT_NOINIT);
    for (i = 0; i < dpi->num_side_sets; i++) {
      dpi->ss_index_global[i] = i;
    }
  }

  dpi->ss_num_distfacts_global = exo->ss_num_distfacts;

  dpi->ss_num_sides_global     = exo->ss_num_sides;

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

void
free_dpi(Dpi *d)
{
  safer_free((void **) &(d->elem_index_global));
  safer_free((void **) &(d->elem_owner));

  safer_free((void **) &(d->elem_elem_list_global));
  safer_free((void **) &(d->elem_elem_twst_global));
  safer_free((void **) &(d->elem_elem_face_global));
  safer_free((void **) &(d->elem_elem_proc_global));

  safer_free((void **) &(d->eb_id_global));
  safer_free((void **) &(d->eb_num_attr_global));
  safer_free((void **) &(d->eb_num_elems_global));
  safer_free((void **) &(d->eb_num_nodes_per_elem_global));
  safer_free((void **) &(d->eb_elem_type_global));

  safer_free((void **) &(d->eb_prop_global));
  safer_free((void **) &(d->ns_prop_global));
  safer_free((void **) &(d->ss_prop_global));
  safer_free((void **) &(d->eb_index_global)); 
  safer_free((void **) &(d->eb_num_private_elems));

  safer_free((void **) &(d->elem_var_tab_global));
  safer_free((void **) &(d->global_node_description));
  safer_free((void **) &(d->node_index_global));
  safer_free((void **) &(d->neighbor));

  safer_free((void **) &(d->ns_id_global));
  safer_free((void **) &(d->ns_num_distfacts_global));
  safer_free((void **) &(d->ns_num_nodes_global));
  safer_free((void **) &(d->ns_node_index_global));
  safer_free((void **) &(d->ns_distfact_index_global));

  safer_free((void **) &(d->ns_index_global));
  safer_free((void **) &(d->ns_node_list_index_global));
  safer_free((void **) &(d->ns_distfact_list_index_global));

  safer_free((void **) &(d->ptr_set_membership));
  safer_free((void **) &(d->set_membership));

  safer_free((void **) &(d->ss_distfact_index_global));
  safer_free((void **) &(d->ss_elem_index_global));
  safer_free((void **) &(d->ss_id_global));
  safer_free((void **) &(d->ss_num_distfacts_global));
  safer_free((void **) &(d->ss_num_sides_global));

  safer_free((void **) &(d->ss_elem_list_index_global));
  safer_free((void **) &(d->ss_distfact_list_index_global));
  safer_free((void **) &(d->ss_index_global));

  safer_free((void **) &(d->ss_block_index_global));
  safer_free((void **) &(d->ss_block_list_global));
  free(d->ss_internal_global);

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

void
free_dpi_uni(Dpi *d)
{
  safer_free((void **) &(d->elem_owner));
  safer_free((void **) &(d->elem_index_global));
  safer_free((void **) &(d->node_index_global));
  safer_free((void **) &(d->eb_index_global));
  safer_free((void **) &(d->global_node_description));
  safer_free((void **) &(d->neighbor));
  safer_free((void **) &(d->ns_index_global));
  safer_free((void **) &(d->ptr_set_membership));
  safer_free((void **) &(d->set_membership));
  safer_free((void **) &(d->ss_index_global));
  free(d->ss_internal_global);

  return;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
 * get_variable() -- provide interface to read netCDF vars via rev2 & rev3
 *
 * Created: 1998/08/21 09:59 MDT pasacki@sandia.gov
 *
 * Revised:
 */
static void
get_variable(const int netcdf_unit,
             const nc_type netcdf_type,
             const int num_dimensions,
             const int dimension_val_1,
             const int dimension_val_2,
             const int variable_identifier,
             void *variable_address)
{
  int err = 0;
  char err_msg[MAX_CHAR_ERR_MSG];

  get_variable_call_count++;

  /*
   * If a variable was really defined properly and doesn't have a valid
   * identifier, then don't even try to get it.
   */

  if ( variable_identifier < 0 )
    {
      return;
    }

  switch ( netcdf_type )
    {
    case NC_INT:
      err = nc_get_var_int(netcdf_unit, variable_identifier, variable_address);
      if ( err != NC_NOERR )
	{
	  sprintf(err_msg, "nc_get_var_int() varid=%d", 
	          variable_identifier);
          WH(-1, err_msg);
	}
      break;

    case NC_CHAR:
      err = nc_get_var_text(netcdf_unit, variable_identifier,
                            variable_address);
      if ( err != NC_NOERR )
        {
          sprintf(err_msg, "nc_get_var_text() varid=%d",
                  variable_identifier);
        }
      break;

    case NC_DOUBLE:
      err = nc_get_var_double(netcdf_unit, variable_identifier,
                              variable_address);
      if ( err != NC_NOERR )
        {
          sprintf(err_msg, "nc_get_var_double() varid=%d",
                  variable_identifier);
        }
        break;

    default:
      EH(-1, "Specified netCDF data type unrecognized or unimplemented.");
      break;
    }

  if (err != NC_NOERR)
    {
      EH(-1, err_msg);
    }


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

void
init_dpi_struct(Dpi *d)
{
  if ( d == NULL )
    {
      EH(-1, "Empty structure to initialize?");
    }
  /*
   * Initialize the entire structure to zeroes and NULL's
   */
  memset((void *)d, 0, sizeof(Dpi));

  /*
   * Initialize anything else that isn't zero.
   */

  d->undefined_basic_eqnvar_id    = UNDEFINED_EQNVARID;
  d->len_node_description         = LEN_NODE_DESCRIPTION;

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

void 
exo_dpi_clone(Exo_DB *exo, 
	      Dpi *dpi)
{
  int len;

  dpi->num_nodes_global         = exo->num_nodes;
  dpi->num_elems_global         = exo->num_elems;
  dpi->num_elem_blocks_global   = exo->num_elem_blocks;
  dpi->num_node_sets_global     = exo->num_node_sets;
  dpi->num_side_sets_global     = exo->num_side_sets;

  /*
   * Allocate and fill arrays for element blocks...
   */

  len = dpi->num_elem_blocks_global * sizeof(int);

  dpi->eb_id_global             = smalloc(len);
  memcpy(dpi->eb_id_global, exo->eb_id, len);

  dpi->eb_num_elems_global      = smalloc(len);
  memcpy(dpi->eb_num_elems_global, exo->eb_num_elems, len);

  /*
   * Allocate and fill arrays for node sets...
   */

  len = dpi->num_node_sets_global * sizeof(int);

  dpi->ns_id_global             = smalloc(len);
  memcpy(dpi->ns_id_global, exo->ns_id, len);

  dpi->ns_num_nodes_global      = smalloc(len);
  memcpy(dpi->ns_num_nodes_global, exo->ns_num_nodes, len);

  dpi->ns_num_distfacts_global  = smalloc(len);
  memcpy(dpi->ns_num_distfacts_global, exo->ns_num_distfacts, len);

  /*
   * Allocate and fill arrays for side sets...
   */
  len = dpi->num_side_sets_global * sizeof(int);
  
  dpi->ss_id_global             = smalloc(len);
  memcpy(dpi->ss_id_global, exo->ss_id, len);

  dpi->ss_num_sides_global      = smalloc(len);
  memcpy(dpi->ss_num_sides_global, exo->ss_num_sides, len);

  dpi->ss_num_distfacts_global  = smalloc(len);
  memcpy(dpi->ss_num_distfacts_global, exo->ss_num_distfacts, len);

  return;
}
