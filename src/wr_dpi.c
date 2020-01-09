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

#ifdef _HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#ifndef lint
#ifdef USE_RCSID
static char rcsid[] = "$Id: wr_dpi.c,v 5.1 2007-09-18 18:53:49 prschun Exp $";
#endif
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
#define NC_MAX_VAR_DIMS		MAX_VAR_DIMS    
#endif

#ifndef NC_MAX_NAME
#define NC_MAX_NAME		MAX_NC_NAME
#endif

#ifndef NC_INT
#define NC_INT			NC_LONG
#endif

/*
 * Might need some NO_NETCDF_2 definitions here ...
 */

#include "std.h"
#include "mm_eh.h"
#include "dpi.h"

/*
 * Prototypes of functions defined here, but needed elsewhere.
 */

static void define_dimension
(const int ,		/* unit */
       const char *,		/* string */
       const int ,		/* value */
       int *);			/* identifier */

static void define_variable
(const int ,		/* netcdf_unit */
       const char *,		/* name_string */
       const nc_type ,		/* netcdf_type */
       const int ,		/* num_dimensions (less than 3 for now) */
       const int ,		/* dimension_id_1 */
       const int ,		/* dimension_id_2 */
       const int ,		/* dimension_val_1 */
       const int ,		/* dimension_val_2 */
       int *);			/* identifier */

static void put_variable(const int netcdf_unit, const nc_type netcdf_type,
                         const int variable_identifier,
                         const void *variable_address);		/* variable_address */

int wr_dpi(Dpi *d, char *filename) {
  int err;
  int status;
  int u;			/* short hand for unit... */

  struct Shadow_Identifiers si;

  status = 0;

  /*
   * From the C interface guide the basic calling sequence is given
   * for the case of adding new dimensions, variables and attributes to
   * an existing netCDF dataset.
   *
   *  nc_open();
   *  nc_redef();
   *    nc_def_dim();
   *    nc_def_var();
   *    nc_put_att();
   *  nc_enddef();
   *  nc_put_var();
   *  nc_close();
   *
   */

  /*
   * Open the file.
   */

  err = nc_open(filename, NC_WRITE, &u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_open() problem.");
    }


  /*
   * Go into define mode.
   */

  err = nc_redef(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_redef() problem.");
    }

  /*
   * Define each of the netCDF dimensions that will be needed to describe
   * the extent of netCDF variables that are arrays.
   *
   * Dimensions that are "zero" are tricky. New way of doing it -- with
   * a function call to handle all the sticky details.
   */
 
  define_dimension(u, DIM_LEN_EB_NUM_PRIVATE_ELEMS, 
		   d->len_eb_num_private_elems,
		   &si.len_eb_num_private_elems);

  define_dimension(u, DIM_LEN_ELEM_VAR_TAB_GLOBAL,
		   d->len_elem_var_tab_global,
		   &si.len_elem_var_tab_global);

  define_dimension(u, DIM_LEN_ELEM_ELEM_LIST, /* new e-e */
		   d->len_elem_elem_list,
		   &si.len_elem_elem_list);

  define_dimension(u, DIM_LEN_NODE_DESCRIPTION,
		   d->len_node_description,
		   &si.len_node_description);

  define_dimension(u, DIM_LEN_NS_NODE_LIST,
		   d->len_ns_node_list,
		   &si.len_ns_node_list);

  define_dimension(u, DIM_LEN_NS_DISTFACT_LIST,
		   d->len_ns_distfact_list,
		   &si.len_ns_distfact_list);

  define_dimension(u, DIM_LEN_SS_ELEM_LIST,
		   d->len_ss_elem_list,
		   &si.len_ss_elem_list);

  define_dimension(u, DIM_LEN_SS_DISTFACT_LIST,
		   d->len_ss_distfact_list,
		   &si.len_ss_distfact_list);

  define_dimension(u, DIM_LEN_STRING,
		   d->len_string,
		   &si.len_string);

  define_dimension(u, DIM_LEN_PTR_SET_MEMBERSHIP,
		   d->len_ptr_set_membership,
		   &si.len_ptr_set_membership);

  define_dimension(u, DIM_LEN_SET_MEMBERSHIP,
		   d->len_set_membership,
		   &si.len_set_membership);

  define_dimension(u, DIM_NUM_ELEM_BLOCKS,
		   d->num_elem_blocks,
		   &si.num_elem_blocks);

  define_dimension(u, DIM_NUM_ELEM_BLOCKS_GLOBAL,
		   d->num_elem_blocks_global,
		   &si.num_elem_blocks_global);

  define_dimension(u, DIM_NUM_ELEMS,
		   d->num_elems,
		   &si.num_elems);

  define_dimension(u, DIM_NUM_GLOBAL_NODE_DESCRIPTIONS,
		   d->num_global_node_descriptions,
		   &si.num_global_node_descriptions);

  define_dimension(u, DIM_NUM_NEIGHBORS,
		   d->num_neighbors,
		   &si.num_neighbors);

  define_dimension(u, DIM_NUM_NODE_SETS,
		   d->num_node_sets,
		   &si.num_node_sets);

  define_dimension(u, DIM_NUM_NODE_SETS_GLOBAL,
		   d->num_node_sets_global,
		   &si.num_node_sets_global);

  define_dimension(u, DIM_NUM_NODES,
		   d->num_nodes,
		   &si.num_nodes);

  define_dimension(u, DIM_NUM_PROPS_EB, d->num_props_eb,
		   &si.num_props_eb);

      define_dimension(u, DIM_NUM_PROPS_NS,
		       d->num_props_ns,
		       &si.num_props_ns);

      define_dimension(u, DIM_NUM_PROPS_SS,
		       d->num_props_ss,
		       &si.num_props_ss);

  define_dimension(u, DIM_NUM_SIDE_SETS,
		   d->num_side_sets,
		   &si.num_side_sets);

  define_dimension(u, DIM_NUM_SIDE_SETS_GLOBAL,
		   d->num_side_sets_global,
		   &si.num_side_sets_global);

  define_dimension(u, DIM_NUM_UNIVERSE_NODES,
		   d->num_universe_nodes,
		   &si.num_universe_nodes);

  
  if (d->num_side_sets_global > 0) {
    define_dimension(u, DIM_LEN_SS_BLOCK_INDEX_GLOBAL,
                     d->num_side_sets_global + 1,
                     &si.len_ss_block_index_global);
    
    define_dimension(u, DIM_LEN_SS_BLOCK_LIST_GLOBAL,
                     d->ss_block_index_global[d->num_side_sets_global],
                     &si.len_ss_block_list_global);
  } else {
      define_dimension(u, DIM_LEN_SS_BLOCK_INDEX_GLOBAL,
                       0,
                       &si.len_ss_block_index_global);
      define_dimension(u, DIM_LEN_SS_BLOCK_LIST_GLOBAL,
                       0,
                       &si.len_ss_block_list_global);
  }

  /*
   * Define variables. Arrays only get defined if their respective dimensions
   * are greater than zero.
   *
   * Also, this handy routine uses two arguments for the possibility of
   * up to 2D arrays. Dummy arguments of "-1" are inserted for 1D arrays 
   * or for scalar variables ( zero dimensional arrays).
   */
  
  define_variable(u, VAR_EB_ELEM_TYPE_GLOBAL, NC_CHAR, 2, 
		  si.num_elem_blocks_global, si.len_string,
		  d->num_elem_blocks_global, d->len_string,
		  &si.eb_elem_type_global);

  define_variable(u, VAR_EB_ID_GLOBAL, NC_INT, 1, 
		  si.num_elem_blocks_global, -1,
		  d->num_elem_blocks_global, -1,
		  &si.eb_id_global);

  define_variable(u, VAR_EB_INDEX_GLOBAL, NC_INT, 1, 
		  si.num_elem_blocks, -1,
		  d->num_elem_blocks, -1,
		  &si.eb_index_global);

  define_variable(u, VAR_EB_NUM_ATTR_GLOBAL, NC_INT, 1, 
		  si.num_elem_blocks_global, -1,
		  d->num_elem_blocks_global, -1,
		  &si.eb_num_attr_global);

  define_variable(u, VAR_EB_NUM_ELEMS_GLOBAL, NC_INT, 1, 
		  si.num_elem_blocks_global, -1,
		  d->num_elem_blocks_global, -1,
		  &si.eb_num_elems_global);

  define_variable(u, VAR_EB_NUM_NODES_PER_ELEM_GLOBAL, NC_INT, 1, 
		  si.num_elem_blocks_global, -1,
		  d->num_elem_blocks_global, -1,
		  &si.eb_num_nodes_per_elem_global);

  define_variable(u, VAR_EB_NUM_PRIVATE_ELEMS, NC_INT, 1, 
		  si.num_elem_blocks, -1,
		  d->num_elem_blocks, -1,
		  &si.eb_num_private_elems);

  if ( d->num_props_eb > 1 )	/* Properties are weird, recall. */
    {
      define_variable(u, VAR_EB_PROP_GLOBAL, NC_INT, 2,
		      si.num_props_eb, si.num_elem_blocks_global,
		      d->num_props_eb, d->num_elem_blocks_global,
		      &si.eb_prop_global);
    }

  if ( d->num_elems > 0 )
    {
      define_variable(u, VAR_ELEM_INDEX_GLOBAL, NC_INT, 1,
		      si.num_elems, -1,
		      d->num_elems, -1,
		      &si.elem_index_global);
    }

  if ( d->len_elem_var_tab_global > 0 )
    {
      define_variable(u, VAR_ELEM_VAR_TAB_GLOBAL, NC_INT, 1,
		      si.len_elem_var_tab_global, -1,
		      d->len_elem_var_tab_global, -1,
		      &si.elem_var_tab_global);
    }

  if ( d->len_elem_elem_list > 0 )
    {
      define_variable(u, VAR_ELEM_OWNER, NC_INT, 1,
		      si.num_elems, -1,
		      d->num_elems, -1,
		      &si.elem_owner);

      define_variable(u, VAR_ELEM_ELEM_LIST_GLOBAL, NC_INT, 1,
		      si.len_elem_elem_list, -1,
		      d->len_elem_elem_list, -1,
		      &si.elem_elem_list_global);

      define_variable(u, VAR_ELEM_ELEM_TWST_GLOBAL, NC_INT, 1,
		      si.len_elem_elem_list, -1,
		      d->len_elem_elem_list, -1,
		      &si.elem_elem_twst_global);

      define_variable(u, VAR_ELEM_ELEM_FACE_GLOBAL, NC_INT, 1,
		      si.len_elem_elem_list, -1,
		      d->len_elem_elem_list, -1,
		      &si.elem_elem_face_global);

      define_variable(u, VAR_ELEM_ELEM_PROC_GLOBAL, NC_INT, 1,
		      si.len_elem_elem_list, -1,
		      d->len_elem_elem_list, -1,
		      &si.elem_elem_proc_global);
    }

  define_variable(u, VAR_GLOBAL_NODE_DESCRIPTION, NC_INT, 2, 
		  si.num_global_node_descriptions, si.len_node_description,
		  d->num_global_node_descriptions, d->len_node_description,
		  &si.global_node_description);

  define_variable(u, VAR_MY_NAME, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.my_name);

  define_variable(u, VAR_NEIGHBOR, NC_INT, 1, 
		  si.num_neighbors, -1,
		  d->num_neighbors, -1,
		  &si.neighbor);

  if ( d->num_nodes > 0 )
    {
      define_variable(u, VAR_NODE_INDEX_GLOBAL, NC_INT, 1,
		      si.num_nodes, -1,
		      d->num_nodes, -1,
		      &si.node_index_global);
    }

  define_variable(u, VAR_NS_DISTFACT_INDEX_GLOBAL, NC_INT, 1,
		  si.num_node_sets_global, -1,
		  d->num_node_sets_global, -1,
		  &si.ns_distfact_index_global);

  define_variable(u, VAR_NS_DISTFACT_LEN_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.ns_distfact_len_global);

  define_variable(u, VAR_NS_DISTFACT_LIST_INDEX_GLOBAL, NC_INT, 1,
		  si.len_ns_distfact_list, -1,
		  d->len_ns_distfact_list, -1,
		  &si.ns_distfact_list_index_global);

  define_variable(u, VAR_NS_ID_GLOBAL, NC_INT, 1, 
		  si.num_node_sets_global, -1,
		  d->num_node_sets_global, -1,
		  &si.ns_id_global);

  define_variable(u, VAR_NS_INDEX_GLOBAL, NC_INT, 1, 
		  si.num_node_sets, -1,
		  d->num_node_sets, -1,
		  &si.ns_index_global);

  define_variable(u, VAR_NS_NODE_INDEX_GLOBAL, NC_INT, 1,
		  si.num_node_sets_global, -1,
		  d->num_node_sets_global, -1,
		  &si.ns_node_index_global);


  define_variable(u, VAR_NS_NODE_LEN_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.ns_node_len_global);

  define_variable(u, VAR_NS_NODE_LIST_INDEX_GLOBAL, NC_INT, 1,
		  si.len_ns_node_list, -1,
		  d->len_ns_node_list, -1,
		  &si.ns_node_list_index_global);

  define_variable(u, VAR_NS_NUM_DISTFACTS_GLOBAL, NC_INT, 1, 
		  si.num_node_sets_global, -1,
		  d->num_node_sets_global, -1,
		  &si.ns_num_distfacts_global);

  define_variable(u, VAR_NS_NUM_NODES_GLOBAL, NC_INT, 1, 
		  si.num_node_sets_global, -1,
		  d->num_node_sets_global, -1,
		  &si.ns_num_nodes_global);

  if ( d->num_props_ns > 1 )
    {
      define_variable(u, VAR_NS_PROP_GLOBAL, NC_INT, 2, 
		      si.num_props_ns, si.num_node_sets_global,
		      d->num_props_ns, d->num_node_sets_global,
		      &si.ns_prop_global);
    }

  define_variable(u, VAR_NUM_BOUNDARY_NODES, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_boundary_nodes);

  define_variable(u, VAR_NUM_DOFS_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_dofs_global);

  define_variable(u, VAR_NUM_ELEMS_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_elems_global);

  define_variable(u, VAR_NUM_EXTERNAL_NODES, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_external_nodes);

  define_variable(u, VAR_NUM_INTERNAL_NODES, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_internal_nodes);

  define_variable(u, VAR_NUM_NODES_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.num_nodes_global);

  define_variable(u, VAR_PTR_SET_MEMBERSHIP, NC_INT, 1, 
		  si.len_ptr_set_membership, -1,
		  d->len_ptr_set_membership, -1,
		  &si.ptr_set_membership);

  define_variable(u, VAR_SET_MEMBERSHIP, NC_INT, 1, 
		  si.len_set_membership, -1,
		  d->len_set_membership, -1,
		  &si.set_membership);

  define_variable(u, VAR_SS_DISTFACT_INDEX_GLOBAL, NC_INT, 1, 
		  si.num_side_sets_global, -1,
		  d->num_side_sets_global, -1,
		  &si.ss_distfact_index_global);

  define_variable(u, VAR_SS_DISTFACT_LIST_INDEX_GLOBAL, NC_INT, 1, 
		  si.len_ss_distfact_list, -1,
		  d->len_ss_distfact_list, -1,
		  &si.ss_distfact_list_index_global);

  define_variable(u, VAR_SS_DISTFACT_LEN_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.ss_distfact_len_global);

  define_variable(u, VAR_SS_ELEM_INDEX_GLOBAL, NC_INT, 1, 
		  si.num_side_sets_global, -1,
		  d->num_side_sets_global, -1,
		  &si.ss_elem_index_global);

  define_variable(u, VAR_SS_ELEM_LEN_GLOBAL, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.ss_elem_len_global);

  define_variable(u, VAR_SS_ELEM_LIST_INDEX_GLOBAL, NC_INT, 1, 
		  si.len_ss_elem_list, -1,
		  d->len_ss_elem_list, -1,
		  &si.ss_elem_list_index_global);

  define_variable(u, VAR_SS_ID_GLOBAL, NC_INT, 1, 
		  si.num_side_sets_global, -1,
		  d->num_side_sets_global, -1,
		  &si.ss_id_global);

  define_variable(u, VAR_SS_INDEX_GLOBAL, NC_INT, 1, 
		  si.num_side_sets, -1,
		  d->num_side_sets, -1,
		  &si.ss_index_global);

  define_variable(u, VAR_SS_NUM_DISTFACTS_GLOBAL, NC_INT, 1, 
		  si.num_side_sets_global, -1,
		  d->num_side_sets_global, -1,
		  &si.ss_num_distfacts_global);

  define_variable(u, VAR_SS_NUM_SIDES_GLOBAL, NC_INT, 1, 
		  si.num_side_sets_global, -1,
		  d->num_side_sets_global, -1,
		  &si.ss_num_sides_global);

  if ( d->num_props_ss > 1 )
    {
      define_variable(u, VAR_SS_PROP_GLOBAL, NC_INT, 2, 
		      si.num_props_ss, si.num_side_sets_global,
		      d->num_props_ss, d->num_side_sets_global,
		      &si.ss_prop_global);
    }

  if ( d->num_side_sets_global > 0 )
    {
      define_variable(u, VAR_SS_INTERNAL_GLOBAL, NC_INT, 1,
                      si.num_side_sets_global, -1,
                      d->num_side_sets_global, -1,
                      &si.ss_internal_global);

      define_variable(u, VAR_SS_BLOCK_INDEX_GLOBAL, NC_INT, 1,
                      si.len_ss_block_index_global, -1,
                      d->num_side_sets_global + 1, -1,
                      &si.ss_block_index_global);

      define_variable(u, VAR_SS_BLOCK_LIST_GLOBAL, NC_INT, 1,
                      si.len_ss_block_list_global, -1,
                      d->ss_block_index_global[d->num_side_sets_global], -1,
                      &si.ss_block_list_global);
    }

  define_variable(u, VAR_UNDEFINED_BASIC_EQNVAR_ID, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.undefined_basic_eqnvar_id);

  /*
   * Leave define mode.
   */

  err = nc_enddef(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_enddef() problem.");
    }

  /*
   * Put variable values.
   *
   * This form is good for scalars, 1d arrays and 2d arrays. Any more and
   * you'll need to add another argument to the list for the backward
   * compatible to netCDF implementation to work properly. We'll assume
   * that start[] arrays that ncvarput() uses will be full of zeroes.
   * If not, then you'll need to do that case by hand.
   */

    put_variable(u, NC_CHAR, si.eb_elem_type_global, &(d->eb_elem_type_global[0][0]));

    put_variable(u, NC_INT, si.eb_id_global, d->eb_id_global);

    put_variable(u, NC_INT, si.eb_index_global, d->eb_index_global);

  put_variable(u, NC_INT, si.eb_num_attr_global, d->eb_num_attr_global);

  put_variable(u, NC_INT, si.eb_num_elems_global, d->eb_num_elems_global);

  put_variable(u, NC_INT, si.eb_num_nodes_per_elem_global, d->eb_num_nodes_per_elem_global);

  put_variable(u, NC_INT, si.eb_num_private_elems, d->eb_num_private_elems);

  if ( d->num_props_eb > 1 )
    {
    put_variable(u, NC_INT, si.eb_prop_global, &(d->eb_prop_global[0][0]));
    }

  if ( d->num_elems > 0 )
    {
    put_variable(u, NC_INT, si.elem_index_global, d->elem_index_global);

    put_variable(u, NC_INT, si.elem_owner, d->elem_owner);

    }

  if ( d->len_elem_var_tab_global > 0 )
    {

    put_variable(u, NC_INT, si.elem_var_tab_global, d->elem_var_tab_global);
    }

  if ( d->len_elem_elem_list > 0 )
    {
    put_variable(u, NC_INT, si.elem_elem_list_global, d->elem_elem_list_global);
    put_variable(u, NC_INT, si.elem_elem_face_global, d->elem_elem_face_global);
    put_variable(u, NC_INT, si.elem_elem_twst_global, d->elem_elem_twst_global);
      put_variable(u, NC_INT, si.elem_elem_proc_global, d->elem_elem_proc_global);
    }

    put_variable(u, NC_INT, si.global_node_description, &(d->global_node_description[0][0]));

    put_variable(u, NC_INT, si.my_name, &(d->my_name));

    put_variable(u, NC_INT, si.neighbor, d->neighbor);

  if ( d->num_nodes > 0 )
    {
    put_variable(u, NC_INT, si.node_index_global, d->node_index_global);
    }

    put_variable(u, NC_INT, si.ns_distfact_len_global, &(d->ns_distfact_len_global));

    put_variable(u, NC_INT, si.ns_node_len_global, &(d->ns_node_len_global));

    put_variable(u, NC_INT, si.ns_id_global, d->ns_id_global);

  put_variable(u, NC_INT, si.ns_index_global, d->ns_index_global);

  put_variable(u, NC_INT, si.ns_distfact_list_index_global, d->ns_distfact_list_index_global);

  put_variable(u, NC_INT, si.ns_distfact_index_global, d->ns_distfact_index_global);

  put_variable(u, NC_INT, si.ns_node_index_global, d->ns_node_index_global);

  put_variable(u, NC_INT, si.ns_node_list_index_global, d->ns_node_list_index_global);

  put_variable(u, NC_INT, si.ns_num_distfacts_global, d->ns_num_distfacts_global);

  put_variable(u, NC_INT, si.ns_num_nodes_global, d->ns_num_nodes_global);

  if ( d->num_props_ns > 1 )
    {
    put_variable(u, NC_INT, si.ns_prop_global, &(d->ns_prop_global[0][0]));
    }

    put_variable(u, NC_INT, si.num_boundary_nodes, &(d->num_boundary_nodes));

    put_variable(u, NC_INT, si.num_dofs_global, &(d->num_dofs_global));

    put_variable(u, NC_INT, si.num_elems_global, &(d->num_elems_global));

  put_variable(u, NC_INT, si.num_external_nodes, &(d->num_external_nodes));

  put_variable(u, NC_INT, si.num_internal_nodes, &(d->num_internal_nodes));

  put_variable(u, NC_INT, si.num_nodes_global, &(d->num_nodes_global));

  put_variable(u, NC_INT, si.ptr_set_membership, d->ptr_set_membership);

  put_variable(u, NC_INT, si.set_membership, d->set_membership);

  put_variable(u, NC_INT, si.ss_distfact_index_global, d->ss_distfact_index_global);

  put_variable(u, NC_INT, si.ss_distfact_len_global, &(d->ss_distfact_len_global));

  if ( d->len_ss_distfact_list > 0 )
    {
    put_variable(u, NC_INT, si.ss_distfact_list_index_global, d->ss_distfact_list_index_global);
    }

    put_variable(u, NC_INT, si.ss_elem_index_global, d->ss_elem_index_global);

    put_variable(u, NC_INT, si.ss_elem_len_global, &(d->ss_elem_len_global));

  if ( d->len_ss_elem_list > 0 )
    {
    put_variable(u, NC_INT, si.ss_elem_list_index_global, d->ss_elem_list_index_global);
    }

    put_variable(u, NC_INT, si.ss_id_global, d->ss_id_global);

  if ( d->num_side_sets > 0 )
    {
    put_variable(u, NC_INT, si.ss_index_global, d->ss_index_global);
    }

    put_variable(u, NC_INT, si.ss_num_distfacts_global, d->ss_num_distfacts_global);

    put_variable(u, NC_INT, si.ss_num_sides_global, d->ss_num_sides_global);

  if ( d->num_props_ss > 1 )
    {
    put_variable(u, NC_INT, si.ss_prop_global, &(d->ss_prop_global[0][0]));
    }

    put_variable(u, NC_INT, si.undefined_basic_eqnvar_id, &(d->undefined_basic_eqnvar_id));

  if (d->num_side_sets_global > 0) {

    put_variable(u, NC_INT, si.ss_internal_global, d->ss_internal_global);

    put_variable(u, NC_INT, si.ss_block_index_global, d->ss_block_index_global);

    put_variable(u, NC_INT, si.ss_block_list_global, d->ss_block_list_global);
  }

  /*
   * Close the file (flush buffers).
   */

  err = nc_close(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_close() problem.");
    }

  return(status);
}

/* define_dimension() - register a dimension name with a value assign ID
 *
 * NetCDF specifies that dimensions be defined with a name and value. When
 * that is done successfully, an integer identifier is assigned to that
 * dimension.
 *
 * We consider the present state of netCDF requiring backward compatability
 * with netCDF 2. Perhaps someday we'll be pure netCDF 3, but that must
 * wait until EXODUS II migrates successfully.
 *
 * Created: 1998/08/19 12:55 MDT pasacki@sandia.gov
 *
 * Revised:
 */

static void
define_dimension(const int unit,
		 const char *string,
		 const int value,
		 int *identifier)
{
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];

  /*
   * Dimensions with value zero are problematic, so filter them out.
   */

  if ( value <= 0 )
    {
      *identifier = -1;
      return;
    }
  err  = nc_def_dim(unit, string, value, identifier);
  if ( err != NC_NOERR )
    {
      sprintf(err_msg, "nc_def_dim() on %s [<%d] id=%d", string, 
		   value, *identifier);
      EH(-1, err_msg);
    }

  return;
}

/* define_variable() - define a netCDF variable type with dimensions
 *
 * If successful, the identifier will be filled with a registered ID for
 * this variable.
 *
 * This routine requires the values of the dimensions as well as their
 * netCDF identifiers so that we can easily determine if a trival zero
 * value is being used. In that case, no variable is created.
 *
 * Created: 1998/08/19 14:17 MDT pasacki@sandia.gov
 *
 * Revised:
 */

static void 
define_variable(const int netcdf_unit,
		const char *name_string,
		const nc_type netcdf_type,
		const int num_dimensions,
		const int dimension_id_1,
		const int dimension_id_2,
		const int dimension_val_1,
		const int dimension_val_2,
		int *identifier)
{			    
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];

  int dim[NC_MAX_VAR_DIMS];

  if ( num_dimensions < 0 )
    {
      EH(-1, "Bad dimension specified.");
    }

  if ( num_dimensions > 0 )
    {
      if ( dimension_id_1 < 0 )
	{
	  return;		/* Do not even create a variable. */
	  /*
	  sprintf(err_msg, "Bad 1st netCDF dimension ID %d for %s\n",
		       dimension_id_1, name_string);
	  EH(-1, err_msg);
	  */
	}
      if ( dimension_val_1 < 1 )
	{
	  return;
	}
    }

  if ( num_dimensions > 1 )
    {
      if ( dimension_id_2 < 0 )
	{
	  sprintf(err_msg, "Bad 2nd netCDF dimension ID %d for %s\n",
		       dimension_id_1, name_string);
	  EH(-1, err_msg);
	}
      if ( dimension_val_2 < 1 )
	{
	  return;
	}
    }

  dim[0] = MAX(dimension_id_1, 1);
  dim[1] = MAX(dimension_id_2, 1);

  err    = nc_def_var(netcdf_unit, name_string, netcdf_type, num_dimensions, 
		      dim, identifier);
  if ( err != NC_NOERR )
    {
      sprintf(err_msg, "nc_def_var on %s (%d-D) id=%d", name_string, 
		   num_dimensions, *identifier);
      EH(-1, err_msg);
    }

  return;
}

static void put_variable(const int netcdf_unit, const nc_type netcdf_type,
                         const int variable_identifier, const void *variable_address) {
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];

  /*
   * If a variable was really defined properly and doesn't have a valid
   * identifier, then don't even try to put it.
   */

  if ( variable_identifier < 0 )
    {
      return;
    }

  switch ( netcdf_type )
    {
    case NC_INT:
      err = nc_put_var_int(netcdf_unit, variable_identifier, variable_address);
      if ( err != NC_NOERR )
	{
	  sprintf(err_msg, "nc_put_var_int() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;
      
    case NC_CHAR:
      err = nc_put_var_text(netcdf_unit, variable_identifier, 
			    variable_address);
      if ( err != NC_NOERR )
	{
	  sprintf(err_msg, "nc_put_var_text() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;

    case NC_DOUBLE:
      err = nc_put_var_double(netcdf_unit, variable_identifier, 
			      variable_address);
      if ( err != NC_NOERR )
	{
	  sprintf(err_msg, "nc_put_var_double() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
	break;

    default:
      EH(-1, "Specified netCDF data type unrecognized or unimplemented.");
      break;
    }


  return;
}

#ifdef YOU_NEED_IT
static char *
string_type(nc_type t)
{
  static char t_int[] = "INT";
  static char t_dbl[] = "DOUBLE";
  static char t_chr[] = "CHAR";
  static char t_flt[] = "FLOAT";
  static char t_byt[] = "BYTE";
  static char t_unk[] = "UNKNOWN";

  switch (t)
    {
    case NC_INT:
      return(t_int);
      break;

    case NC_DOUBLE:
      return(t_dbl);
      break;

    case NC_CHAR:
      return(t_chr);
      break;

    case NC_FLOAT:
      return(t_flt);
      break;

    case NC_BYTE:
      return(t_byt);
      break;

    default:
      return(t_unk);
      break;
    }
      
  return(t_unk);
}
#endif
