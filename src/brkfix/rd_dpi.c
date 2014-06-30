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
 * Created: 1997/07/10 15:39 MDT pasacki@sandia.gov
 *
 * Revised: 1997/07/21 09:45 MDT pasacki@sandia.gov
 *
 * Revised: 1998/09/26 07:01 MDT pasacki@sandia.gov
 */

#define _RD_DPI_C

#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>
/*
*  Use NETCDF_3 with NO_NETCDF_2 not defined!
* #define NO_NETCDF_2		 for pure netCDF 3 - but still not yet.
*/
#define NETCDF_3


#include "netcdf.h"

/*
 * I like these symbols as more lucid indicators of which netCDF we are
 * using...
 *
 * Well, we need to maintain backwards compatibility...
 */

#define _RD_DPI_C

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
#include "rd_dpi.h"
#include "nodesc.h"

static int get_variable_call_count = 0;

#if 0
static int ProcID=-1;
#endif

/*
 * Prototypes of functions defined here and needed only here.
 */

static void get_variable
PROTO((const int netcdf_unit,
       const nc_type netcdf_type,
       const int num_dimensions,
       const int dimension_val_1,
       const int dimension_val_2,
       const int variable_identifier,
       void *variable_address));





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
   * Set all elements in the structure to NULL to ensure that all
   * pointers are set to NULL. Thus, we can check against NULL
   * to see if a deallocation is necessary.
   */
  (void) memset((void *)d, 0, sizeof(Dpi));
  
  /*
   * Initialize all of the dimensions to zero.
   */

  d->len_eb_num_private_elems     = 0;
  d->len_elem_var_tab_global      = 0;
  d->len_elem_elem_list           = 0;
  d->len_node_description         = 0;
  d->len_ns_node_list             = 0;
  d->len_ns_distfact_list         = 0;
  d->len_ss_elem_list             = 0;
  d->len_ss_distfact_list         = 0;

  d->len_string                   = MAX_STR_LENGTH;

  d->len_ptr_set_membership       = 0;
  d->len_set_membership           = 0;
  d->num_elem_blocks              = 0;  
  d->num_elem_blocks_global       = 0;
  d->num_elems                    = 0;
  d->num_global_node_descriptions = 0;
  d->num_neighbors                = 0;  
  d->num_node_sets                = 0;
  d->num_node_sets_global         = 0;
  d->num_nodes                    = 0;
  d->num_props_eb                 = 0;
  d->num_props_ns                 = 0;
  d->num_props_ss                 = 0;
  d->num_side_sets                = 0;
  d->num_side_sets_global         = 0;
  d->num_universe_nodes           = 0;

  d->undefined_basic_eqnvar_id    = UNDEFINED_EQNVARID;
  d->len_node_description         = LEN_NODE_DESCRIPTION;

  return;
}

/* init_dpi_version() -- make space and initialize the Dpi version string
 *
 * Take HKM's advice. Make space for a version string to manage backward
 * incompatibilities with more grace and less hacking than previously 
 * possible...
 *
 * Note: This only makes space! You still must put a version string in the
 *       variable that reflects how old or new is the data you're reading
 *       or writing.
 *
 * Created: 1999/11/15 13:28 MST pasacki@sandia.gov
 */

void 
init_dpi_version(Dpi *d)
{
  int i;
  int len;

  if ( d == NULL )
    {
      EH(-1, "Allocate b4 U initialize!");
    }

  len = d->len_string;

  if ( len < 1 )
    {
      EH(-1, "Not enough space to place a version string!");
    }

  d->dpi_version_string = (char *) smalloc((len+3)*sizeof(char));

  for ( i=0; i<len+3; i++)
    {
      d->dpi_version_string[i] = '\0';
    }

  return;
}


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


int
rd_dpi(Dpi *d,
       char *fn,
       int verbosity)
{
  int err;
  int i;
  int len;

#ifdef NETCDF_3
  int dim[NC_MAX_VAR_DIMS];
#endif
#ifdef NETCDF_2
  /* int  tim[NC_MAX_VAR_DIMS];	 This is ludicrous... */
  long dim[NC_MAX_VAR_DIMS];
#endif
  long count[NC_MAX_VAR_DIMS];
  long start[NC_MAX_VAR_DIMS];

  int status;
  int u;			/* short hand for unit... */

  struct Shadow_Identifiers si;

  status = 0;

  for ( i=0; i<NC_MAX_VAR_DIMS; i++)
    {
      start[i] = 0;
      dim[i]   = 1;
      count[i] = 1;
    }

  /* Gratuitous use to placate SGI cc warnings! */

  if ( start[0] != 0 || dim[0] != 1 || count[0] != 1 ) exit(2);

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

#ifdef DEBUG
  fprintf(stderr, "just ncopen-ed file \"%s\" on unit %d\n", fn, u);
#endif
  
  /*
   * 2. Get dimension identifiers.
   *
   * These are determined from their names which are defined in dpi.h. 
   *
   * These integer dimension IDs, if read properly from the open file, 
   * are stuck into the Shadow Identifiers si structure for later use.
   *
   * Briefly, the "TRUE" and "FALSE" arguments you see refer to whether
   * it is critical that this particular dimension exist in the database.
   * Often, a FALSE flag is desirable for objects of zero length that
   * do not exist in this case.
   */

  getdid(u, DIM_LEN_EB_NUM_PRIVATE_ELEMS,     FALSE,
	 &si.len_eb_num_private_elems);
  getdid(u, DIM_LEN_ELEM_VAR_TAB_GLOBAL,      FALSE,
	 &si.len_elem_var_tab_global);
  getdid(u, DIM_LEN_NODE_DESCRIPTION,         TRUE,
	 &si.len_node_description);

  /*
   * Found some cases of split up pieces that didn't have any lists,
   * i.e., there lengths are zero. In that case, netCDF likes it better
   * to not even ask.
   */

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

  /*
   * 3. Using the dimension IDs, inquire of the dimension values. Load
   *    those values into the Dpi structure.
   *
   *    These dimension values are important so that we know how much
   *    space to allocate to hold array variables below...
   *
   * getdim(netcdf_unit, dimension_id, address of the answer);
   */

  getdim(u, si.len_eb_num_private_elems,     &d->len_eb_num_private_elems);
  getdim(u, si.len_elem_var_tab_global,      &d->len_elem_var_tab_global);
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

#ifdef DEBUG
  fprintf(stderr, "Processor %d has d->num_side_sets = %d\n", ProcID,
	  d->num_side_sets);
  fprintf(stderr, "Processor %d has d->num_elem_blocks = %d\n", ProcID,
	  d->num_elem_blocks);
  fprintf(stderr, "Processor %d has d->num_side_sets_global = %d\n", ProcID,
	  d->num_side_sets_global);
#endif  

  /*
   * 4. Get variable identifiers from netCDF.
   *
   *    The Booleans are set only roughly based on what "ought" to exist.
   *    If you deem looser or tighter criteria, then by all means change
   *    them to suit your needs.
   */

  getvid(u, VAR_DPI_VERSION_STRING,        FALSE, 
	 &si.dpi_version_string);
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

  getvid(u, VAR_ELEM_VAR_TAB_GLOBAL,	   FALSE, /* New! */
	 &si.elem_var_tab_global);
  getvid(u, VAR_GLOBAL_NODE_DESCRIPTION,   TRUE,
	 &si.global_node_description);
  getvid(u, VAR_GLOBAL_NODE_DOF0,          TRUE,
	 &si.global_node_dof0);
  getvid(u, VAR_GLOBAL_NODE_KIND,          TRUE,
	 &si.global_node_kind);
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

  /*
   * 5. Allocate space to hold array variables.
   *
   *    a. Should verify that these dimensions are reasonable
   *       numbers 0 < dimval < humongous
   *
   *    b. Simple scalar variables are not allocated here - they got
   *       their space when *d became meaningful.
   */

  len = d->len_string;
  if ( len > 0 )
    {
      d->dpi_version_string = (char *) smalloc(len*sizeof(char));
      for ( i=0; i<len; i++)
	{
	  d->dpi_version_string[i] = '\0';
	}
    }

  len = d->num_elems;
  if ( len > 0 )
    {
      d->elem_index_global = (int *) smalloc(len*sizeof(int));
    }

  len = d->num_elem_blocks_global;
  if ( len > 0 )
    {
      d->eb_id_global                 = (int *) smalloc(len*sizeof(int));
      d->eb_num_attr_global           = (int *) smalloc(len*sizeof(int));
      d->eb_num_elems_global          = (int *) smalloc(len*sizeof(int));
      d->eb_num_nodes_per_elem_global = (int *) smalloc(len*sizeof(int));

      d->eb_elem_type_global          = (char **) smalloc(len*sizeof(char *));

      d->eb_elem_type_global[0] = (char *) smalloc(len*(MAX_STR_LENGTH+1)*
						   sizeof(char));
      for ( i=1; i<len; i++)
	{
	  d->eb_elem_type_global[i] = ( d->eb_elem_type_global[i-1] +
					(MAX_STR_LENGTH+1) );
	}
    }

  len = d->num_props_eb;
  if ( len > 1 )		/* Props are special ID is implicit. */
    {
      d->eb_prop_global = (int **) smalloc(len * sizeof(int *));
      d->eb_prop_global[0] = (int *) smalloc(len * (d->num_elem_blocks_global) * sizeof(int));
      for ( i=1; i<len; i++)
	{
	  d->eb_prop_global[i] = ( d->eb_prop_global[i-1] + d->num_elem_blocks_global );
	}
    }

  if ( d->num_props_ns > 1 )	/* Props are special, ID is implicit. */
    {
      d->ns_prop_global = (int **) smalloc(d->num_props_ns * sizeof(int *));
      d->ns_prop_global[0] = (int *) smalloc(d->num_props_ns*
					     d->num_node_sets_global*
					     sizeof(int));
      for ( i=1; i<d->num_props_ns; i++)
	{
	  d->ns_prop_global[i] = ( d->ns_prop_global[i-1] + 
				   d->num_node_sets_global );
	}
    }

  if ( d->num_props_ss > 1 )	/* Props are special, ID is implicit. */
    {
      d->ss_prop_global = (int **) smalloc(d->num_props_ss * sizeof(int *));

      d->ss_prop_global[0] = (int *) smalloc(d->num_props_ss*
					     d->num_side_sets_global*
					     sizeof(int));
      for ( i=1; i<d->num_props_ss; i++)
	{
	  d->ss_prop_global[i] = ( d->ss_prop_global[i-1] + 
				   d->num_side_sets_global );
	}
    }

  len = d->num_elem_blocks;
  if ( len > 0 )
    {
      d->eb_index_global      = (int *) smalloc(len * sizeof(int));
      d->eb_num_private_elems = (int *) smalloc(len * sizeof(int));
    }

  len = d->len_elem_var_tab_global;
  if ( len > 0 )
    {
      d->elem_var_tab_global = ( int *) smalloc(len*sizeof(int));
    }

  if ( d->num_global_node_descriptions > 0 )
    {
      d->global_node_description = 
	(int **) smalloc(d->num_global_node_descriptions * sizeof(int *));

      d->global_node_description[0] = 
	(int *) smalloc(d->num_global_node_descriptions*
			d->len_node_description*sizeof(int));

      for ( i=1; i<d->num_global_node_descriptions; i++)
	{
	  d->global_node_description[i] = ( d->global_node_description[i-1] +
					    d->len_node_description );
	}
    }

  len = d->num_nodes;
  if ( len > 0 )
    {
      d->node_index_global = (int *) smalloc(len*sizeof(int));
    }

  if ( d->num_universe_nodes > 0 ) 
    {
      d->global_node_dof0 = (int *) smalloc(d->num_universe_nodes*sizeof(int));
      d->global_node_kind = (int *) smalloc(d->num_universe_nodes*sizeof(int));
    }

  if ( d->num_neighbors > 0 )
    {
      d->neighbor = (int *) smalloc(d->num_neighbors*sizeof(int));
    }


  len = d->num_node_sets_global;
  if ( len > 0 )
    {
      d->ns_id_global             = (int *) smalloc(len*sizeof(int));
      d->ns_num_distfacts_global  = (int *) smalloc(len*sizeof(int));
      d->ns_num_nodes_global      = (int *) smalloc(len*sizeof(int));
      d->ns_node_index_global     = (int *) smalloc(len*sizeof(int));
      d->ns_distfact_index_global = (int *) smalloc(len*sizeof(int));
    }

  if ( d->num_node_sets > 0 )
    {
      d->ns_index_global = (int *) smalloc(d->num_node_sets*sizeof(int));
    }

  if ( d->len_ns_node_list > 0 )
    {
      d->ns_node_list_index_global = (int *) 
	smalloc(d->len_ns_node_list*sizeof(int));
    }

  if ( d->len_ns_distfact_list > 0 )
    {
      d->ns_distfact_list_index_global = (int *) 
	smalloc(d->len_ns_distfact_list*sizeof(int));
    }

  if ( d->len_ptr_set_membership > 0 )
    {
      d->ptr_set_membership = (int *) smalloc(d->len_ptr_set_membership*
					      sizeof(int));
    }

  if ( d->len_set_membership > 0 )
    {
      d->set_membership = (int *) smalloc(d->len_set_membership*sizeof(int));
    }

  len = d->num_side_sets_global;
  if ( len > 0 )
    {
      d->ss_distfact_index_global = (int *) smalloc(len*sizeof(int));
      d->ss_elem_index_global     = (int *) smalloc(len*sizeof(int));

      d->ss_id_global             = (int *) smalloc(len*sizeof(int));
      d->ss_num_distfacts_global  = (int *) smalloc(len*sizeof(int));
      d->ss_num_sides_global      = (int *) smalloc(len*sizeof(int));
    }

  if ( d->len_ss_elem_list > 0 )
    {
      d->ss_elem_list_index_global = (int *) 
	smalloc(d->len_ss_elem_list*sizeof(int));
    }

  if ( d->len_ss_distfact_list > 0 )
    {
      d->ss_distfact_list_index_global = (int *) 
	smalloc(d->len_ss_distfact_list*sizeof(int));
    }

  if ( d->num_side_sets > 0 )
    {
      d->ss_index_global         = (int *) smalloc(d->num_side_sets*
						   sizeof(int));
    }

  /*
   * 6. Get variables - this is messy. Netcdf 3 wants a different routine
   *    for each type. NetCDF 2 wants to know all the dimensions! Attempt
   *    to hide the messy details in get_variable()...
   */


  if ( si.dpi_version_string < 0 )
    {
      strcpy(d->dpi_version_string, "0.9");
    }
  else
    {
      get_variable(u, NC_CHAR, 1, 
	       d->len_string,			-1, 
	       si.dpi_version_string,		d->dpi_version_string);
    }

  get_variable(u, NC_CHAR, 2, 
	       d->num_elem_blocks_global,	d->len_string,
	       si.eb_elem_type_global,	&(d->eb_elem_type_global[0][0]));

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

  get_variable(u, NC_INT, 1, 
	       d->len_elem_var_tab_global,	-1, 
	       si.elem_var_tab_global,		d->elem_var_tab_global);

  if ( d->num_elems > 0 )
    {
      get_variable(u, NC_INT, 1, 
		   d->num_elems,	 -1, 
		   si.elem_index_global, d->elem_index_global);
    }

  get_variable(u, NC_INT, 2, 
	       d->num_global_node_descriptions,	d->len_node_description, 
	       si.global_node_description,&(d->global_node_description[0][0]));

  get_variable(u, NC_INT, 1, 
	       d->num_universe_nodes,	-1, 
	       si.global_node_dof0,	d->global_node_dof0);

  get_variable(u, NC_INT, 1, 
	       d->num_universe_nodes,	-1, 
	       si.global_node_kind,	d->global_node_kind);

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

  if ( d->num_props_ss > 1 )
    {
      get_variable(u, NC_INT, 2,
		   d->num_props_ss,		d->num_side_sets_global,
		   si.ss_prop_global,		&(d->ss_prop_global[0][0]));
    }

  get_variable(u, NC_INT, 0, 
	       -1,	-1, 
	       si.undefined_basic_eqnvar_id, &(d->undefined_basic_eqnvar_id));

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
  Spfrtn sr=0;

#ifdef DEBUG
  fprintf(stderr, "getdid(unit=%d, \"%s\", %d, 0x%x\n", netcdf_unit,
	  string_name, hard_error_interpretation, 
	  dimension_identifier_address);
#endif

#ifdef NETCDF_3  
  err  = nc_inq_dimid(netcdf_unit, string_name, dimension_identifier_address);
  /*
   * Assume an error means this quanitity is not found here. That's OK for
   * some things, but not others.
   */
  if ( err != NC_NOERR && hard_error_interpretation )
    {
      sr = sprintf(err_msg, "nc_inq_dimid() on %s id=%d", 
		   string_name, 
		   *dimension_identifier_address);
      EH(-1, err_msg);
    }
#endif

#ifdef NETCDF_2
  err  = ncdimid(netcdf_unit, string_name);
  if ( err == -1 && hard_error_interpretation )
    {
      sr   = sprintf(err_msg, "ncdimid() on %s rtn %d", string_name, err);
      EH(err, err_msg);
      EH(sr, err_msg);
    }
  *dimension_identifier_address = err;
#endif  

#ifdef DEBUG
  fprintf(stderr, "P_%d dim id is %d for %s\n", ProcID, 
	  *dimension_identifier_address, string_name);
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
  Spfrtn sr=0;


#ifdef NETCDF_3
  *variable_identifier_address = -1;
  err  = nc_inq_varid(netcdf_unit, string_name, variable_identifier_address);
  if ( err != NC_NOERR && hard_error_interpretation )
    {
      sr = sprintf(err_msg, "nc_inq_varid() on %s id=%d", 
		   string_name, 
		   *variable_identifier_address);
      EH(-1, err_msg);
    }
#endif

#ifdef NETCDF_2
  err  = ncvarid(netcdf_unit, string_name);
  if ( err == -1 && hard_error_interpretation )
    {
      sr   = sprintf(err_msg, "ncvarid() on %s rtn %d", string_name, err);
      EH(err, err_msg);
      EH(sr, err_msg);
    }
  *variable_identifier_address = err;
#endif  

#ifdef DEBUG
  fprintf(stderr, "P_%d var id is %d for %s\n", ProcID, 
	  *variable_identifier_address, string_name);
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
  char err_msg[MAX_CHAR_ERR_MSG];  
  Spfrtn sr=0;
#ifdef NETCDF_3
  size_t swhere;
#endif
#ifdef NETCDF_2
  long swhere;
#endif
  
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
	  sr = sprintf(err_msg, "nc_inq_dimlen() on did=%d", dimension_id);
	  EH(-1, err_msg);
	}
      *where = (int) swhere;
#endif

#ifdef NETCDF_2
      err  = ncdiminq(netcdf_unit, dimension_id, junk, &swhere);
      sr   = sprintf(err_msg, "ncdiminq() on did %d rtns %d", dimension_id, 
		     err);
      *where = (int) swhere;
      EH(err, err_msg);
      EH(sr, err_msg);
#endif  
    }
  return;
}

/*
 * These routines are not needed or wanted in brk/fix - define them only
 * for goma...
 */

#ifdef _RF_GOMA_H

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

  dpi->num_internal_nodes	     = exo->num_nodes;
  dpi->num_boundary_nodes            = 0;
  dpi->num_external_nodes            = 0;

  dpi->num_universe_nodes            = exo->num_nodes;
  
  /*
   * Note! This aliasing of these pointers into the exo structure has
   * two advantages and one disadvantage.
   *
   *	(+) it's very easy to do
   *	(+) it makes more economic use of memory
   *	(-) it makes it too easy to free_dpi() and nuke your EXODUS information
   *
   * We'll just do it for now!
   */

  dpi->eb_id_global = exo->eb_id;

  if ( dpi->num_elem_blocks > 0 )
    {
      dpi->eb_index_global = (int *) smalloc(dpi->num_elem_blocks*sizeof(int));
      for ( i=0; i<dpi->num_elem_blocks; i++)
	{
	  dpi->eb_index_global[i] = i;
	}
    }

  dpi->eb_num_elems_global  = exo->eb_num_elems;

  dpi->eb_num_private_elems = exo->eb_num_elems;

  if ( dpi->num_global_node_descriptions > 0 )
    {
      dpi->global_node_description = 
	(int **) smalloc(dpi->num_global_node_descriptions * sizeof(int *));

      dpi->global_node_description[0] = 
	(int *) smalloc(dpi->num_global_node_descriptions*
			dpi->len_node_description*
			sizeof(int));

      for ( i=1; i<dpi->num_global_node_descriptions; i++)
	{
	  dpi->global_node_description[i] = dpi->global_node_description[i-1] +
	    dpi->len_node_description;
	}
    }
  
  dpi->global_node_dof0        = First_Unknown;

  /*
   * Well, yes, this looks fishy. To avoid wasting memory, just point this
   * to some array with a length that is the number of nodes in the problem
   * and then hope that no one really depends on the global_node_kind for
   * serial problems. If and when you really make use of it, then go ahead
   * and allocate space to hold the proper integer labels.
   */

  dpi->global_node_kind        = First_Unknown;

  dpi->neighbor                = &ProcID;

  dpi->ns_id_global            = exo->ns_id;

  if ( dpi->num_node_sets > 0 )
    {
      dpi->ns_index_global = (int *) smalloc(dpi->num_node_sets*sizeof(int));
      for ( i=0; i<dpi->num_node_sets; i++)
	{
	  dpi->ns_index_global[i] = i;
	}
    }

  dpi->ns_num_distfacts_global = exo->ns_num_distfacts;

  dpi->ns_num_nodes_global     = exo->ns_num_nodes;

  if ( dpi->len_ptr_set_membership > 0 )
    {
      dpi->ptr_set_membership      = 
	(int *) smalloc(dpi->len_ptr_set_membership*sizeof(int));
    }

  if ( dpi->len_set_membership > 0 )
    {
      dpi->set_membership = (int *) smalloc(dpi->len_set_membership*
					    sizeof(int));
    }

  dpi->set_membership[0]       = -1;
  dpi->ptr_set_membership[0]   = 0;
  dpi->ptr_set_membership[1]   = 1;

  dpi->ss_id_global            = exo->ss_id;

  if ( dpi->num_side_sets > 0 )
    {
      dpi->ss_index_global = (int *) smalloc(dpi->num_side_sets*sizeof(int));
      for ( i=0; i<dpi->num_side_sets; i++)
	{
	  dpi->ss_index_global[i] = i;
	}
    }

  dpi->ss_num_distfacts_global = exo->ss_num_distfacts;

  dpi->ss_num_sides_global     = exo->ss_num_sides;

  return;
}

#endif

/* free_dpi() -- free internal dynamically allocated memory in a Dpi struct
 *
 * During rd_dpi(), various arrays are allocated. To cleanly free up this
 * memory, this routine can be used. Typically, if dpi is a (Dpi *), then
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
  safe_free(d->eb_elem_type_global[0]);
  safe_free(d->eb_elem_type_global);

  safe_free(d->eb_id_global);
  safe_free(d->eb_index_global);
  safe_free(d->eb_num_attr_global);
  safe_free(d->eb_num_elems_global);
  safe_free(d->eb_num_nodes_per_elem_global);
  safe_free(d->eb_num_private_elems);

  if ( d->num_props_eb > 1 )
    {
      safe_free(d->eb_prop_global[0]);
      safe_free(d->eb_prop_global);
    }
  
  if ( d->num_elems > 0 )
    {
      safe_free(d->elem_index_global);
    }

  safe_free(d->global_node_description[0]);
  safe_free(d->global_node_description);

  safe_free(d->global_node_dof0);
  safe_free(d->global_node_kind);
  safe_free(d->neighbor);

  if ( d->num_nodes > 0 )
    {
      safe_free(d->node_index_global);
    }

  if ( d->len_ns_node_list > 0 )
    {
      safe_free(d->ns_node_list_index_global);
    }

  if ( d->len_ns_distfact_list > 0 )
    {
      safe_free(d->ns_distfact_list_index_global);
    }

  safe_free(d->ns_distfact_index_global);
  safe_free(d->ns_id_global);

  if ( d->num_node_sets > 0 )
    {
      safe_free(d->ns_index_global);
    }

  safe_free(d->ns_node_index_global);
  safe_free(d->ns_num_distfacts_global);
  safe_free(d->ns_num_nodes_global);

  if ( d->num_props_ns > 1 )	
    {
      safe_free(d->ns_prop_global[0]);
      safe_free(d->ns_prop_global);
    }

  safe_free(d->ptr_set_membership);
  safe_free(d->set_membership);

  if ( d->len_ss_distfact_list > 0 )
    {
      safe_free(d->ss_distfact_list_index_global);
    }

  if ( d->len_ss_elem_list > 0 )
    {
      safe_free(d->ss_elem_list_index_global);
    }

  if ( d->num_side_sets_global > 0 )
    {
      safe_free(d->ss_distfact_index_global);
      safe_free(d->ss_elem_index_global);
      safe_free(d->ss_num_sides_global);
      safe_free(d->ss_num_distfacts_global);
      safe_free(d->ss_id_global);
    }

  if ( d->num_side_sets > 0 )
    {
      safe_free(d->ss_index_global);
    }

  if ( d->num_props_ss > 1 )	
    {
      safe_free(d->ss_prop_global[0]);
      safe_free(d->ss_prop_global);
    }

  if ( d->len_elem_elem_list > 0 )
    {
      safe_free(d->elem_owner);
      safe_free(d->elem_elem_list_global);
      safe_free(d->elem_elem_twst_global);
      safe_free(d->elem_elem_face_global);
      safe_free(d->elem_elem_proc_global);
    }

  return;
}

/* free_dpi_uni() -- free internal dynamically allocated memory in Dpi SERIAL!
 *
 * During rd_dpi(), various arrays are allocated. To cleanly free up this
 * memory, this routine can be used. Typically, if dpi is a (Dpi *), then
 *
 *	free_dpi(dpi);
 *	free(dpi);
 *
 * would accomplish the desired effect.
 *
 * For serial processing, not all of the pieces of Dpi have been allocated
 * and some are merely aliases to parts of the exodus ii database. To avoid
 * munging that data, just free up what was allocated in uni_dpi.
 *
 * Created: 1997/09/11 10:43 MDT pasacki@sandia.gov
 */

void
free_dpi_uni(Dpi *d)
{
  safe_free(d->eb_index_global);
  safe_free(d->global_node_description);
  safe_free(d->ns_index_global);
  safe_free(d->ptr_set_membership);
  safe_free(d->set_membership);
  safe_free(d->ss_index_global);

  return;
}

/*
 * get_variable() -- provide interface to read netCDF vars via rev2 & rev3
 *
 * NetCDF versions 2 and 3 provide different interfaces. We attempt to
 * accomodate the basic version 3 interface with backward compatible calls
 * surrounded by relevent preprocessor directives.
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
  int err;
  int i;
#ifndef NO_NETCDF_2
  long count[NC_MAX_VAR_DIMS];
  long start[NC_MAX_VAR_DIMS];
#endif
  char err_msg[MAX_CHAR_ERR_MSG];
  Spfrtn sr=0;

  get_variable_call_count++;

  /*
   * If a variable was really defined properly and doesn't have a valid
   * identifier, then don't even try to get it.
   */

  if ( variable_identifier < 0 )
    {
      return;
    }

#ifdef NO_NETCDF_2		/* pure netCDF 3 calls... */
  switch ( netcdf_type )
    {
    case NC_INT:
      err = nc_get_var_int(netcdf_unit, variable_identifier, variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_get_var_int() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;
      
    case NC_CHAR:
      err = nc_get_var_text(netcdf_unit, variable_identifier, 
			    variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_get_var_text() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
      break;

    case NC_DOUBLE:
      err = nc_get_var_double(netcdf_unit, variable_identifier, 
			      variable_address);
      if ( err != NC_NOERR )
	{
	  sr = sprintf(err_msg, "nc_get_var_double() varid=%d", 
		       variable_identifier);
	  EH(-1, err_msg);
	}
	break;

    default:
      EH(-1, "Specified netCDF data type unrecognized or unimplemented.");
      break;
    }
#endif

#ifndef NO_NETCDF_2		/* backward compatibility mode to netcdf2 */
  
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
	  EH(-1, "Bad dimension specification.");
	}
      count[0] = dimension_val_1;
    }

  if ( num_dimensions > 1 ) 
    {
      if ( dimension_val_2 < 1 )
	{
	  EH(-1, "Bad dimension specification.");
	}
      count[1] = dimension_val_2;
    }

  err = ncvarget(netcdf_unit, variable_identifier, start, count, 
		 variable_address);
  if ( err < 0 )
    {
      sr = sprintf(err_msg, 
		   "get_variable (%d call), varid %d (%d dim %d,%d)\n",
		   get_variable_call_count,
		   variable_identifier, num_dimensions, 
		   dimension_val_1, dimension_val_2);
      EH(err, err_msg);
    }

#endif

  return;
}
