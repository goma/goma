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

/* wr_dpi.c */

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

#define _WR_DPI_C

#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

/* #define NO_NETCDF_2		 for pure netCDF 3 -- still impure since
 * EXODUS II v3.00 still requires backwards compatibility 
 */
#define NETCDF_3

#include "netcdf.h"

/*
 * Sigh, if you need to run netCDF 2 then here's some definitions to tide
 * you over until netCDF 3 is working for EXODUS II...
 */

/*
 * I like these symbols as more lucid indicators of which netCDF we are
 * using...
 */



//#ifdef  NO_NETCDF_2
//#define NETCDF_3
//#endif

//#ifndef NO_NETCDF_2
//#define NETCDF_2
//#endif

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

#include "map_names.h"
#include "std.h"
#include "aalloc.h"
#include "eh.h"
#include "exo_struct.h"
#include "dpi.h"

/*
 * Prototypes of functions defined here, but needed here only.
 */

#if 0				/* Consign to attic... */
static char *string_type
PROTO((nc_type));		/* netcdf_type */
#endif

static void define_dimension
PROTO((const int ,		/* unit, netcdf */
       const char *,		/* string */
       const int ,		/* value */
       int *));			/* identifier */

static void define_variable
PROTO((const int ,		/* netcdf_unit */
       const char *,		/* name_string */
       nc_type ,		/* netcdf_type */
       const int ,		/* num_dimensions */
       const int ,		/* dimension_id_1 */
       const int ,		/* dimension_id_2 */
       const int ,		/* dimension_val_1 */
       const int ,		/* dimension_val_2 */
       int *));			/* identifier */

static void put_variable
PROTO((const int ,		/* netcdf_unit */
       nc_type ,		/* netcdf_type */
       const int ,		/* num_dimensions */
       const int ,		/* dimension_val_1 */
       const int ,		/* dimension_val_2 */
       const int ,		/* variable_identifier */
       const void *));		/* variable_address */

/*
 * Prototypes of functions defined elsewhere, but needed here.
 */

/*
 * Prototypes of functions defined here, but needed elsewhere.
 */

extern int wr_dpi
PROTO((Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* filename */
       int ));			/* verbosity - how much to talk */


int
wr_dpi(Dpi *d,
       char *filename,
       int verbosity)
{
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

#ifdef NETCDF_3
  err = nc_open(filename, NC_WRITE, &u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_open() problem.");
    }
#endif
#ifdef NETCDF_2
  err = ncopen(filename, NC_WRITE);
  EH(err, "ncopen() problem.");
  u   = err;
#endif

  /*
   * Go into define mode.
   */

#ifdef NETCDF_3
  err = nc_redef(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_redef() problem.");
    }
#endif
#ifdef NETCDF_2
  err = ncredef(u);
  EH(err, "ncredef() problem.");
#endif

  /*
   * Define each of the netCDF dimensions that will be needed to describe
   * the extent of netCDF variables that are arrays.
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

  /*
   * Define variables. Arrays only get defined if their respective dimensions
   * are greater than zero.
   *
   * Also, this handy routine uses two arguments for the possibility of
   * up to 2D arrays. Dummy arguments of "-1" are inserted for 1D arrays 
   * or for scalar variables ( zero dimensional arrays).
   */
  
  define_variable(u, VAR_DPI_VERSION_STRING, NC_CHAR, 1,
		  si.len_string, -1,
		  d->len_string, -1,
		  &si.dpi_version_string);

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

  define_variable(u, VAR_GLOBAL_NODE_DOF0, NC_INT, 1, 
		  si.num_universe_nodes, -1,
		  d->num_universe_nodes, -1,
		  &si.global_node_dof0);

  define_variable(u, VAR_GLOBAL_NODE_KIND, NC_INT, 1, 
		  si.num_universe_nodes, -1,
		  d->num_universe_nodes, -1,
		  &si.global_node_kind);

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

  define_variable(u, VAR_UNDEFINED_BASIC_EQNVAR_ID, NC_INT, 0, 
		  -1, -1,
		  -1, -1,
		  &si.undefined_basic_eqnvar_id);

  /*
   * Leave define mode.
   */

#ifdef NETCDF_3
  err = nc_enddef(u);
  if ( err != NC_NOERR )
    {
      EH(-1, "nc_enddef() problem.");
    }
#endif
#ifdef NETCDF_2
  err = ncendef(u);
  EH(err, "ncendef() problem.");
#endif

  /*
   * Put variable values.
   *
   * This form is good for scalars, 1d arrays and 2d arrays. Any more and
   * you'll need to add another argument to the list for the backward
   * compatible to netCDF implementation to work properly. We'll assume
   * that start[] arrays that ncvarput() uses will be full of zeroes.
   * If not, then you'll need to do that case by hand.
   */

  put_variable(u, NC_CHAR, 1, 
	       d->len_string,			-1, 
	       si.dpi_version_string,		d->dpi_version_string);

  put_variable(u, NC_CHAR, 2, 
	       d->num_elem_blocks_global,	d->len_string,
	       si.eb_elem_type_global, &(d->eb_elem_type_global[0][0]));

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_id_global,			d->eb_id_global);

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks,		-1, 
	       si.eb_index_global,		d->eb_index_global);

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_attr_global,		d->eb_num_attr_global);

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_elems_global,		d->eb_num_elems_global);

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks_global,	-1, 
	       si.eb_num_nodes_per_elem_global,	
	       d->eb_num_nodes_per_elem_global);

  put_variable(u, NC_INT, 1, 
	       d->num_elem_blocks,		-1, 
	       si.eb_num_private_elems,		d->eb_num_private_elems);

  if ( d->num_props_eb > 1 )
    {
      put_variable(u, NC_INT, 2,
		   d->num_props_eb, d->num_elem_blocks_global,
		   si.eb_prop_global, &(d->eb_prop_global[0][0]));
    }

  if ( d->num_elems > 0 )
    {
      put_variable(u, NC_INT, 1, 
		   d->num_elems, -1, 
		   si.elem_index_global,	d->elem_index_global);

      put_variable(u, NC_INT, 1, d->num_elems,	-1, si.elem_owner,
		   d->elem_owner);

    }

  if ( d->len_elem_var_tab_global > 0 )
    {
      put_variable(u, NC_INT, 1, 
		   d->len_elem_var_tab_global,	-1, 
		   si.elem_var_tab_global,	d->elem_var_tab_global);
    }

  if ( d->len_elem_elem_list > 0 )
    {
      put_variable(u, NC_INT, 1, d->len_elem_elem_list,	-1, 
		   si.elem_elem_list_global, d->elem_elem_list_global);
      put_variable(u, NC_INT, 1, d->len_elem_elem_list,	-1, 
		   si.elem_elem_face_global, d->elem_elem_face_global);
      put_variable(u, NC_INT, 1, d->len_elem_elem_list,	-1, 
		   si.elem_elem_twst_global, d->elem_elem_twst_global);
      put_variable(u, NC_INT, 1, d->len_elem_elem_list,	-1, 
		   si.elem_elem_proc_global, d->elem_elem_proc_global);
    }

  put_variable(u, NC_INT, 2, 
	       d->num_global_node_descriptions,	d->len_node_description, 
	       si.global_node_description,&(d->global_node_description[0][0]));

  put_variable(u, NC_INT, 1, 
	       d->num_universe_nodes,	-1, 
	       si.global_node_dof0,	d->global_node_dof0);

  put_variable(u, NC_INT, 1, 
	       d->num_universe_nodes,	-1, 
	       si.global_node_kind,	d->global_node_kind);

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.my_name,		&(d->my_name));

  put_variable(u, NC_INT, 1, 
	       d->num_neighbors,	-1, 
	       si.neighbor,		d->neighbor);

  if ( d->num_nodes > 0 )
    {
      put_variable(u, NC_INT, 1, 
		   d->num_nodes,	-1, 
		   si.node_index_global,	d->node_index_global);
    }

  put_variable(u, NC_INT, 0, 
	       -1, 	-1, 
	       si.ns_distfact_len_global, &(d->ns_distfact_len_global));

  put_variable(u, NC_INT, 0, 
	       -1, 	-1, 
	       si.ns_node_len_global, &(d->ns_node_len_global));

  put_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,	-1, 
	       si.ns_id_global,		d->ns_id_global);

  put_variable(u, NC_INT, 1, 
	       d->num_node_sets,	-1, 
	       si.ns_index_global,	d->ns_index_global);

  put_variable(u, NC_INT, 1, 
	       d->len_ns_distfact_list,	-1, 
	       si.ns_distfact_list_index_global, 
	       d->ns_distfact_list_index_global);

  put_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_distfact_index_global,	d->ns_distfact_index_global);

  put_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_node_index_global,	d->ns_node_index_global);

  put_variable(u, NC_INT, 1, 
	       d->len_ns_node_list,	-1, 
	       si.ns_node_list_index_global, 
	       d->ns_node_list_index_global);



  put_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,		-1, 
	       si.ns_num_distfacts_global,	d->ns_num_distfacts_global);

  put_variable(u, NC_INT, 1, 
	       d->num_node_sets_global,	-1, 
	       si.ns_num_nodes_global,	d->ns_num_nodes_global);

  if ( d->num_props_ns > 1 )
    {
      put_variable(u, NC_INT, 2,
		   d->num_props_ns,		d->num_node_sets_global,
		   si.ns_prop_global,	&(d->ns_prop_global[0][0]));
    }

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_boundary_nodes,	&(d->num_boundary_nodes));

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_dofs_global,	&(d->num_dofs_global));

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_elems_global,	&(d->num_elems_global));

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_external_nodes,	&(d->num_external_nodes));

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_internal_nodes,	&(d->num_internal_nodes));

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.num_nodes_global,	&(d->num_nodes_global));

  put_variable(u, NC_INT, 1, 
	       d->len_ptr_set_membership,	-1, 
	       si.ptr_set_membership,	d->ptr_set_membership);

  put_variable(u, NC_INT, 1, 
	       d->len_set_membership,	-1, 
	       si.set_membership,	d->set_membership);

  put_variable(u, NC_INT, 1, d->num_side_sets_global, -1, 
	       si.ss_distfact_index_global, d->ss_distfact_index_global);

  put_variable(u, NC_INT, 0, -1, -1, 
	       si.ss_distfact_len_global, &(d->ss_distfact_len_global));

  if ( d->len_ss_distfact_list > 0 )
    {
      put_variable(u, NC_INT, 1, d->len_ss_distfact_list, -1, 
		   si.ss_distfact_list_index_global, 
		   d->ss_distfact_list_index_global);
    }

  put_variable(u, NC_INT, 1, d->num_side_sets_global, -1, 
	       si.ss_elem_index_global, d->ss_elem_index_global);

  put_variable(u, NC_INT, 0, -1, -1, 
	       si.ss_elem_len_global, &(d->ss_elem_len_global));

  if ( d->len_ss_elem_list > 0 )
    {
      put_variable(u, NC_INT, 1, d->len_ss_elem_list, -1, 
		   si.ss_elem_list_index_global, d->ss_elem_list_index_global);
    }

  put_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,	-1, 
	       si.ss_id_global,		d->ss_id_global);

  if ( d->num_side_sets > 0 )
    {
      put_variable(u, NC_INT, 1, 
		   d->num_side_sets,	-1, 
		   si.ss_index_global,	d->ss_index_global);
    }

  put_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,		-1, 
	       si.ss_num_distfacts_global,	d->ss_num_distfacts_global);

  put_variable(u, NC_INT, 1, 
	       d->num_side_sets_global,		-1, 
	       si.ss_num_sides_global,		d->ss_num_sides_global);

  if ( d->num_props_ss > 1 )
    {
      put_variable(u, NC_INT, 2,
		   d->num_props_ss,		d->num_side_sets_global,
		   si.ss_prop_global,		&(d->ss_prop_global[0][0]));
    }

  put_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.undefined_basic_eqnvar_id, &(d->undefined_basic_eqnvar_id));

  /*
   * Close the file (flush buffers).
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
  Spfrtn sr=0;

  /*
   * Dimensions with value zero are problematic, so filter them out.
   */

  if ( value <= 0 )
    {
      *identifier = -1;
      return;
    }
#ifdef NETCDF_3  
  err  = nc_def_dim(unit, string, value, identifier);
  if ( err != NC_NOERR )
    {
      sr = sprintf(err_msg, "nc_def_dim() on %s [<%d] id=%d", string, 
		   value, *identifier);
      EH(-1, err_msg);
    }
#endif
#ifdef NETCDF_2
  err  = ncdimdef(unit, string, value);
  sr   = sprintf(err_msg, "ncdimdef() on %s [<%d] rtn %d", string, 
		 value, err);
  EH(err, err_msg);
  EH(sr, err_msg);
  *identifier = err;
#endif  

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
  Spfrtn sr=0;

#ifdef NETCDF_3
  int dim[NC_MAX_VAR_DIMS];
#endif
#ifdef NETCDF_2
  int  tim[NC_MAX_VAR_DIMS];			/* This is ludicrous... */
  long dim[NC_MAX_VAR_DIMS];
#endif

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
	  sr = sprintf(err_msg, "Bad 1st netCDF dimension ID %d for %s\n",
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
	  sr = sprintf(err_msg, "Bad 2nd netCDF dimension ID %d for %s\n",
		       dimension_id_1, name_string);
	  EH(-1, err_msg);
	  EH(sr, err_msg);
	}
      if ( dimension_val_2 < 1 )
	{
	  return;
	}
    }

  dim[0] = MAX(dimension_id_1, 1);
  dim[1] = MAX(dimension_id_2, 1);

#ifdef NETCDF_3
  err    = nc_def_var(netcdf_unit, name_string, netcdf_type, num_dimensions, 
		      dim, identifier);
  if ( err != NC_NOERR )
    {
      sr = sprintf(err_msg, "nc_def_var on %s (%d-D) id=%d", name_string, 
		   num_dimensions, *identifier);
      EH(-1, err_msg);
    }
#endif
#ifdef NETCDF_2
  for ( i=0; i<num_dimensions; i++ ) 
    {
      tim[i] = dim[i];
    }

  err    = ncvardef(netcdf_unit, name_string, netcdf_type, num_dimensions, 
		    tim);
  sr     = sprintf(err_msg, "ncvardef on %s (%d-D) id=%d", name_string, 
		   num_dimensions, err);
  EH(err, err_msg);
  *identifier   = err;
#endif

  return;
}

static void
put_variable(const int netcdf_unit,
	     const nc_type netcdf_type,
	     const int num_dimensions,
	     const int dimension_val_1,
	     const int dimension_val_2,
	     const int variable_identifier,
	     const void *variable_address)
{
  int err;
  char err_msg[MAX_CHAR_ERR_MSG];
  Spfrtn sr=0;

  /*
   * If a variable was really defined properly and doesn't have a valid
   * identifier, then don't even try to put it.
   */

  if ( variable_identifier < 0 )
    {
      return;
    }

#ifdef NETCDF_3		/* pure netCDF 3 calls... */
  switch ( netcdf_type )
    {
    case NC_INT:
      err = nc_put_var_int(netcdf_unit, variable_identifier, variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_put_var_int() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;
      
    case NC_CHAR:
      err = nc_put_var_text(netcdf_unit, variable_identifier, 
			    variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_put_var_text() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;

    case NC_DOUBLE:
      err = nc_put_var_double(netcdf_unit, variable_identifier, 
			      variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_put_var_double() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
	break;

    default:
      EH(-1, "Specified netCDF data type unrecognized or unimplemented.");
      break;
    }
#endif

#ifdef NETCDF_2		/* backward compatibility mode to netcdf2 */
  
  if ( num_dimensions < 0 || num_dimensions > 2 )
    {
      sr = sprintf(err_msg, "Bad or too large dimension value %d", 
		   num_dimensions);
      EH(-1, err_msg);
      EH(sr, err_msg);
    }

  for ( i=0; i<num_dimensions; i++)
    {
      count[i] = 1;
      start[i] = 0;
    }

  if ( num_dimensions > 0 )
    {
      if ( dimension_val_1 < 1 )
	{
	  /*
	   * This is a case where you didn't really want to put any
	   * variables here, right? I mean, really, a dimension of zero?
	   */
	  return;
	  /*
	  sr = sprintf(err_msg, 
	  "Data type %s (%d dimensional) has dimension values = %d %d",
		       string_type(netcdf_type), num_dimensions, 
		       dimension_val_1, dimension_val_2);
	  EH(-1, err_msg);
	  */
	}
      count[0] = dimension_val_1;
    }

  if ( num_dimensions > 1 ) 
    {
      if ( dimension_val_2 < 1 )
	{
	  /*
	   * Another case of where you were not serious. Here a 2D array
	   * has the second dimension of zero, so it is pretty darn flat.
	   * So flat, in fact, that we're going to call it done and leave
	   * right away.
	   */
	  return;
	  /*
	  sr = sprintf(err_msg, 
	  "Data type %s (%d dimensional) has dimension values = %d %d",
		       string_type(netcdf_type), num_dimensions, 
		       dimension_val_1, dimension_val_2);
	  EH(-1, err_msg);
	  */
	}
      count[1] = dimension_val_2;
    }

  /*
   * If we've made it here, then presumably the dimensions have some meat
   * to them.
   */

  err = ncvarput(netcdf_unit, variable_identifier, start, count, 
		 variable_address);
  EH(err, "ncvarput");
#endif

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
