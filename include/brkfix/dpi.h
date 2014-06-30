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

/* dpi.h -- definitions useful for distributed processing
 *
 * This defines useful data for distributed processing tasks used by GOMA
 * for mixed finite element solution of coupled conjugate problems.
 *
 * 
 * Notes:
 *	  [1]	Defines a data structure with information that should be
 *              useful in distributed computing contexts.
 *
 *        [2]   Defines some C preprocessor definitions useful for accessing
 *              same data through netCDF.
 *
 *	  [3]   Documention about each variable is still kind of sketchy.
 *		Fill it in more later...
 *
 *	  [4]   Attempt to handle both netCDF 3.3 and back down to netCDF 2.4.3
 *		that EXODUS II seems to rely upon.
 *
 *	  [5]   Note every piece of data in the DPI structure has a 
 *		corresponding entry in the Shadow Identifiers at the bottom.
 *		Those exists as a means of containing the integer identifiers
 *		that netCDF associates with each dimension or variable.
 *
 *	  [6]   Instead of usurping the elem_map[] and node_map[] from
 *		EXODUS II for our devious parallel purposes, make strict
 *		dpi versions of same and name accordingly.
 *
 * ! Where the names are different (eg., num_node_sets_proc ) from the
 *   structure components, it's because the name conflicts with one already
 *   in use by the EXODUS II database.
 *
 * Created: 1997/05/17 14:08 MDT pasacki@sandia.gov
 *
 * Revised: 1998/08/10 10:36 MDT pasacki@sandia.gov
 *
 * Revised: 1998/09/25 16:29 MDT pasacki@sandia.gov
 */

#ifndef _DPI_H
#define _DPI_H

/*
 * netCDF Dimensions.
 */

#define DIM_LEN_EB_NUM_PRIVATE_ELEMS		"len_eb_num_private_elems"
#define DIM_LEN_ELEM_ELEM_LIST			"len_elem_elem_list"
#define DIM_LEN_ELEM_VAR_TAB_GLOBAL		"len_elem_var_tab_global"
#define DIM_LEN_NODE_DESCRIPTION		"len_node_description"
#define DIM_LEN_NS_NODE_LIST			"len_ns_node_list"
#define DIM_LEN_NS_DISTFACT_LIST		"len_ns_distfact_list"
#define DIM_LEN_SS_ELEM_LIST			"len_ss_elem_list"
#define DIM_LEN_SS_DISTFACT_LIST		"len_ss_distfact_list"
#define DIM_LEN_STRING				"len_string_dpi"
#define DIM_LEN_PTR_SET_MEMBERSHIP		"len_ptr_set_membership"
#define DIM_LEN_SET_MEMBERSHIP			"len_set_membership"
#define DIM_NUM_ELEM_BLOCKS			"num_elem_blocks_proc"
#define DIM_NUM_ELEM_BLOCKS_GLOBAL		"num_elem_blocks_global"
#define DIM_NUM_ELEMS				"num_elems_proc"
#define DIM_NUM_GLOBAL_NODE_DESCRIPTIONS	"num_global_node_descriptions"
#define DIM_NUM_NEIGHBORS			"num_neighbors"
#define DIM_NUM_NODE_SETS			"num_node_sets_proc"
#define DIM_NUM_NODE_SETS_GLOBAL		"num_node_sets_global"
#define DIM_NUM_NODES				"num_nodes_proc"
#define DIM_NUM_PROPS_EB			"num_props_eb"
#define DIM_NUM_PROPS_NS			"num_props_ns"
#define DIM_NUM_PROPS_SS			"num_props_ss"
#define DIM_NUM_SIDE_SETS			"num_side_sets_proc"
#define DIM_NUM_SIDE_SETS_GLOBAL		"num_side_sets_global"
#define DIM_NUM_UNIVERSE_NODES			"num_universe_nodes"

/*
 * netCDF Variables.
 */

#define VAR_DPI_VERSION_STRING			"dpi_version_string"
#define VAR_EB_ELEM_TYPE_GLOBAL			"eb_elem_type_global"
#define VAR_EB_ID_GLOBAL			"eb_id_global"
#define VAR_EB_INDEX_GLOBAL			"eb_index_global"
#define VAR_EB_NUM_ATTR_GLOBAL			"eb_num_attr_global"
#define VAR_EB_NUM_ELEMS_GLOBAL			"eb_num_elems_global"
#define VAR_EB_NUM_NODES_PER_ELEM_GLOBAL	"eb_num_nodes_per_elem_global"
#define VAR_EB_NUM_PRIVATE_ELEMS		"eb_num_private_elems"
#define VAR_EB_PROP_GLOBAL			"eb_prop_global"
#define VAR_EB_PTR_GLOBAL			"eb_ptr_global"
#define VAR_ELEM_INDEX_GLOBAL			"elem_index_global"
#define VAR_ELEM_VAR_TAB_GLOBAL			"elem_var_tab_global"
#define VAR_ELEM_OWNER				"elem_owner"
#define VAR_ELEM_ELEM_LIST_GLOBAL		"elem_elem_list_global"
#define VAR_ELEM_ELEM_FACE_GLOBAL		"elem_elem_face_global"
#define VAR_ELEM_ELEM_TWST_GLOBAL		"elem_elem_twst_global"
#define VAR_ELEM_ELEM_PROC_GLOBAL		"elem_elem_proc_global"
#define VAR_GLOBAL_NODE_DESCRIPTION		"global_node_description"
#define VAR_GLOBAL_NODE_DOF0			"global_node_dof0"
#define VAR_GLOBAL_NODE_KIND			"global_node_kind"
#define VAR_MY_NAME				"my_name"
#define VAR_NEIGHBOR				"neighbor"
#define VAR_NODE_INDEX_GLOBAL			"node_index_global"
#define VAR_NS_DISTFACT_INDEX_GLOBAL		"ns_distfact_index_global"
#define VAR_NS_DISTFACT_LEN_GLOBAL		"ns_distfact_len_global"
#define VAR_NS_DISTFACT_LIST_INDEX_GLOBAL	"ns_distfact_list_index_global"
#define VAR_NS_ID_GLOBAL			"ns_id_global"
#define VAR_NS_INDEX_GLOBAL			"ns_index_global"
#define VAR_NS_NODE_INDEX_GLOBAL		"ns_node_index_global"
#define VAR_NS_NODE_LEN_GLOBAL			"ns_node_len_global"
#define VAR_NS_NODE_LIST_INDEX_GLOBAL		"ns_node_list_index_global"
#define VAR_NS_NUM_DISTFACTS_GLOBAL		"ns_num_distfacts_global"
#define VAR_NS_NUM_NODES_GLOBAL			"ns_num_nodes_global"
#define VAR_NS_PROP_GLOBAL			"ns_prop_global"
#define VAR_NUM_BOUNDARY_NODES			"num_boundary_nodes"
#define VAR_NUM_DOFS_GLOBAL			"num_dofs_global"
#define VAR_NUM_ELEMS_GLOBAL			"num_elems_global"
#define VAR_NUM_EXTERNAL_NODES			"num_external_nodes"
#define VAR_NUM_INTERNAL_NODES			"num_internal_nodes"
#define VAR_NUM_NODES_GLOBAL			"num_nodes_global"
#define VAR_PTR_SET_MEMBERSHIP			"ptr_set_membership"
#define VAR_SET_MEMBERSHIP			"set_membership"
#define VAR_SS_DISTFACT_INDEX_GLOBAL		"ss_distfact_index_global"
#define VAR_SS_DISTFACT_LEN_GLOBAL		"ss_distfact_len_global"
#define VAR_SS_DISTFACT_LIST_INDEX_GLOBAL	"ss_distfact_list_index_global"
#define VAR_SS_ELEM_INDEX_GLOBAL		"ss_elem_index_global"
#define VAR_SS_ELEM_LEN_GLOBAL			"ss_elem_len_global"
#define VAR_SS_ELEM_LIST_INDEX_GLOBAL		"ss_elem_list_index_global"
#define VAR_SS_ID_GLOBAL			"ss_id_global"
#define VAR_SS_INDEX_GLOBAL			"ss_index_global"
#define VAR_SS_NODE_LEN_GLOBAL			"ss_node_len_global"
#define VAR_SS_NUM_DISTFACTS_GLOBAL		"ss_num_distfacts_global"
#define VAR_SS_NUM_SIDES_GLOBAL			"ss_num_sides_global"
#define VAR_SS_PROP_GLOBAL			"ss_prop_global"
#define VAR_UNDEFINED_BASIC_EQNVAR_ID		"undefined_basic_eqnvar_id"


struct Distributed_Processing_Information
{
  /*
   * dimensions (array dimensions).
   */

  int len_eb_num_private_elems;
  int len_elem_var_tab_global;	/* New! */
  int len_elem_elem_list;	/* for connectivities */
  int len_node_description;
  int len_ns_node_list;		/* New! */
  int len_ns_distfact_list;	/* New! */
  int len_ss_elem_list;		/* New! */
  int len_ss_distfact_list;	/* New! */
  int len_string;		/* New! */

  int len_ptr_set_membership;
  int len_set_membership;
  int num_elem_blocks;
  int num_elem_blocks_global;	
  int num_elems;		/* New! */
  int num_global_node_descriptions;
  int num_neighbors;
  int num_node_sets;
  int num_node_sets_global;	
  int num_nodes;		/* New! */
  int num_props_eb;		/* New! */
  int num_props_ns;		/* New! */
  int num_props_ss;		/* New! */
  int num_side_sets;
  int num_side_sets_global;
  int num_universe_nodes;

  /*
   * variables (arrays, or scalars not used as dimensions).
   */
  
  char *dpi_version_string;	/* [len_dpi_all_purpose_string] */
  int *eb_id_global;		/* [neb_global] */
  char **eb_elem_type_global;	/* New! - [neb_global][MAX_STR_LENGTH] */
  int *eb_index_global;		/* [neb_proc] */
  int *eb_num_attr_global;	/* New! - [neb_global] */
  int *eb_num_elems_global;	/* [neb_global] */
  int *eb_num_nodes_per_elem_global; /* New! - [neb_global] */
  int *eb_num_private_elems;	/* [neb_proc] */
  int **eb_prop_global;		/* New! [num_props_eb][num_elem_blocks_glob] */

  int *elem_index_global;	/* New! [num_elems_proc] */

  int *elem_var_tab_global;	/* New! [num_elem_blocks_glob*num_elem_vars]
				 *    = [len_elem_var_tab_global]
				 * This flattened 2d->1d array has the index
				 * for the element variables cycling faster. */

  /*
   * Some new stuff to help in the assembly of element based methods, such
   * as discontinuous Galerkin.
   */



  int *elem_owner;		/* [num_elems_proc] Which proc owns this e? */
				/* These are referenced using the pointers
				 * from exo->elem_elem_pntr[elem] */

  int *elem_elem_list_global;	/* Face-ordered names of elems facing this. */
  int *elem_elem_twst_global;	/* Twist of facing elements. */
  int *elem_elem_face_global;	/* Face name of facing element. */
  int *elem_elem_proc_global;	/* Owning proc of facing element. */


  int **global_node_description;/* [num_global_node_descriptions]
				 *                  [len_node_description] */
  int *global_node_dof0;	/* [num_universe_nodes] */
  int *global_node_kind;	/* [num_universe_nodes] */
  int my_name;			/* scalar */
  int *neighbor;		/* [num_neighbors] */
  int *node_index_global;	/* New! - [num_nodes_proc] */
  int *ns_id_global;		/* [num_node_sets_global] */
  int *ns_index_global;		/* [num_node_sets] */

  int ns_distfact_len_global;	/* New! - Length of nodeset dist facts(glob) */
  int ns_node_len_global;	/* New! - Length of nodeset node lists(glob) */

  int *ns_num_distfacts_global;	/* [num_node_sets_global] */
  int *ns_num_nodes_global;	/* [num_node_sets_global] */


  int *ns_node_list_index_global; /* New! - where to map into ns_node_list[] 
				   * global - length is [ns_node_len] 
				   * Mucho convenient for concatenated arrays
				   * like this!
				   */

  int *ns_node_index_global;	/* New! [num_node_sets_global] */
  int *ns_distfact_index_global;/* New! [num_node_sets_global] */

  int *ns_distfact_list_index_global; /* New! - where to put distfacts from
				       * this processor in the global list
				       * length is [ns_distfact_len]
				       * Mucho convenient for concatenated 
				       * arrays like this!
				       */

  int **ns_prop_global;		/* New! [ns_num_props][num_node_sets_global] */

  int num_boundary_nodes;	/* scalar - number of nodes for which this 
				 * processor has primary responsibility but 
				 * which are associated with unknowns that 
				 * need to be sent out to other neighboring 
				 * processors (external from their view).
				 */

  int num_dofs_global;		/* scalar */

  int num_elems_global;		/* scalar */

  int num_external_nodes;	/* scalar
				 * Number of nodes for which this processor
				 * has secondary responsibility. Such nodes
				 * are associated with values that need
				 * to be gathered from other processors and
				 * are used in calculations. 
				 */

  int num_internal_nodes;	/* scalar
				 * Number of nodes for which this processor
				 * has primary responsibility but which
				 * never need to be communicated to/from
				 * neighboring processors. 
				 */

  int num_nodes_global;		/* scalar */

  int *ptr_set_membership;	/* [len_ptr_set_membership] */

  int *set_membership;		/* [len_set_membership] */

  int ss_distfact_len_global;	/* New! - Length of ss distfacts list (glob) */
  int ss_elem_len_global;	/* New! - Length of ss elem list (glob) */

  int *ss_id_global;		/* [num_side_sets_global] */
  int *ss_distfact_index_global; /* New! - [num_side_sets_global] */
  int *ss_elem_index_global;	/* New! - [num_side_sets_global] */

  int *ss_elem_list_index_global; /* New! - [ss_elem_len] maps to the right
				   * place in the global ss_elem_list[] array
				   * Mucho convenient for concatenated arrays
				   * like this!
				   */

  int *ss_distfact_list_index_global; /* New! - [ss_distfact_len] maps
				      * every local distfact into the right
				      * place in the global distfact list.
				      * I.e, these numbers have no local
				      * meaning but as a map to the global
				      * counterpart...
				      */

  int *ss_index_global;		/* [num_side_sets] */

  int *ss_num_distfacts_global;	/* [num_side_sets_global] */

  int *ss_num_sides_global;	/* [num_side_sets_global] */

  int **ss_prop_global;		/* New! - [ss_num_props][num_side_sets_glob] */

  int undefined_basic_eqnvar_id; /* scalar - pad ends of the rectangular shaped
				  * 2d arrays of node descriptions with these*/

};

typedef struct Distributed_Processing_Information Dpi;

/*
 * When you update the structure above, then you'll need to make room for
 * an idenfier below if you want to facilitate netCDF i/o later on.
 *
 * What are these things? Well, netCDF likes to assign integer names to
 * the dimensions and variables as well as the strings you like to use.
 *
 * These integer identifiers are stored in this structure.
 *
 * Yes, yes, I know, that a better alternative to this distributed
 * marbles on a sidewalk is possible. If it bothers you because you know
 * the RIGHT THING TO DO and have the time, please have at it...
 */

struct Shadow_Identifiers
{
  /*
   * dimensions (array dimensions).
   */
  
  int len_eb_num_private_elems;
  int len_elem_var_tab_global;
  int len_elem_elem_list;
  int len_node_description;

  int len_ns_node_list;
  int len_ns_distfact_list;
  int len_ss_elem_list;
  int len_ss_distfact_list;
  int len_string;

  int len_ptr_set_membership;
  int len_set_membership;
  int num_elem_blocks;
  int num_elem_blocks_global;	
  int num_elems;
  int num_global_node_descriptions;
  int num_neighbors;
  int num_node_sets;
  int num_node_sets_global;	
  int num_nodes;
  int num_props_eb;
  int num_props_ns;
  int num_props_ss;
  int num_side_sets;
  int num_side_sets_global;
  int num_universe_nodes;

  /*
   * variables (arrays).
   */
  
  int dpi_version_string;
  int eb_elem_type_global;
  int eb_id_global;
  int eb_index_global;
  int eb_num_attr_global;
  int eb_num_elems_global;
  int eb_num_nodes_per_elem_global;
  int eb_num_private_elems;
  int eb_prop_global;
  int elem_index_global;
  int elem_var_tab_global;
  int elem_owner;
  int elem_elem_list_global;
  int elem_elem_face_global;
  int elem_elem_twst_global;
  int elem_elem_proc_global;
  int global_node_description;
  int global_node_dof0;
  int global_node_kind;
  int my_name;
  int neighbor;
  int node_index_global;
  int ns_distfact_index_global;
  int ns_distfact_len_global;
  int ns_distfact_list_index_global;
  int ns_id_global;		/* List of global node set identifiers. */
  int ns_index_global;		/* Global nodeset INDEX for each local */
  int ns_node_index_global;
  int ns_node_len_global;
  int ns_node_list_index_global;
  int ns_num_distfacts_global;	/* Number of distribution factors in ea */
  int ns_num_nodes_global;	/* Number of nodes in ea global node set. */
  int ns_prop_global;
  int num_boundary_nodes;	/* Number of nodes for which this processor
				 * has primary responsibility but which are 
				 * associated with unknowns that need to be
				 * sent out to other neighboring processors. */
  int num_dofs_global;
  int num_elems_global;		  
  int num_external_nodes;	/* Number of nodes for which this processor
				 * has secondary responsibility. Such nodes
				 * are associated with values that need
				 * to be gathered from other processors and
				 * are used in calculations. */

  int num_internal_nodes;	/* Number of nodes for which this processor
				 * has primary responsibility but which
				 * never need to be communicated to/from
				 * neighboring processors. */

  int num_nodes_global;		

  int ptr_set_membership;

  int set_membership;	

  int ss_distfact_index_global;
  int ss_distfact_len_global;
  int ss_distfact_list_index_global;
  int ss_elem_index_global;
  int ss_elem_len_global;
  int ss_elem_list_index_global;

  int ss_id_global;		/* List of global side set identifiers. */
  int ss_index_global;		/* Global sideset INDEX for each local
				 * set/proc sideset index */
  int ss_node_len_global;
  int ss_num_distfacts_global;	/* Number of distribution factors in ea
				 * global side set. */
  int ss_num_sides_global;	/* Number of sides(&elems) in ea global ss. */
  int ss_prop_global;
  int undefined_basic_eqnvar_id; /* Pad ends of rect. 2d arrays w/ these. */
};

/*
 * Yes, some of this information will already be present as part of
 * the companion EXODUS II finite element data model. 
 *
 * It's replicated here for convenience, mostly in cases where the value 
 * in question is needed for an array dimension.
 *
 * Legend:
 *
 * Name				Description
 * ----				-----------
 *
 * len_eb_num_private_elems	Integer describes the length of the
 *				eb_num_private_elems[] array. Thus, the value
 *				of this variable is typically equal to the
 *				number of element blocks traversed by this
 *				set/proc.
 *
 * len_global_node_dof0		Integer describes the length of the
 *				global_node_dof0[] array. This should be the
 *				number of nodes traversed by this processor,
 *				internal, boundary and external.
 *
 * len_global_node_kind		Integer is the length of the global_node_kind[]
 *				array. Its value should be the number of nodes
 *				known to this processor: (int+bnd+ext).
 *
 * len_node_description		Integer is the length of a node description.
 *				For the current node description structure, the
 *				length of a node description is 4*MAXEQNVARS+1.
 *
 *  int num_universe_nodes;	 * The sum of all the previous three.
 *
 * Caution! The indeces appropriate for the arrays below are *global*
 *
 * On a given proc/set, there are two equivalent means for determining
 * the equivalent global identities:
 *
 *	eb_id = E->eb_id[local_index];
 *
 *    eb_id = D->eb_id_global[D->eb_index_global[local_index]];
 *
 * and likewise for nodesets and sidesets.
 */

#endif
