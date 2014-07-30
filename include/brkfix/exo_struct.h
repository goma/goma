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

/* exo_struct.h -- defines structure to hold an EXODUS II FE database
 *
 * Fri Apr 14 10:33 MDT 1995 pasacki@sandia.gov
 */

#ifndef _EXO_STRUCT_H
#define	_EXO_STRUCT_H 1

#include "std.h"		/* for "dbl" typedef, among other things.. */

#include "exodusII.h"

/*
 * For some reason this handy routine has been omitted from the
 * standard include file definitions of prototypes.
 *
 * Well, I think Larry put it back in lately...
 */

#if 0
extern int ex_get_side_set_node_list
PROTO((int,			/* exoid */
       int ,			/* side_set_id */
       int *,			/* side_set_node_cnt_list - length num sides */
       int *));			/* side_set_node_list - length num dfs */
#endif

#if 0
#ifndef _DBL_TYPEDEF
#define _DBL_TYPEDEF
typedef double dbl;
#endif
#endif

#ifndef _FLT_TYPEDEF
#define _FLT_TYPEDEF
typedef float flt;
#endif

typedef char  *(QA_Record[4]);
typedef char  *INFO_Record;

#ifndef MAX_INFO
#define MAX_INFO	(200)	/* Hopefully irrelevant... */
#endif

#ifndef MAX_QA
#define MAX_QA		(100)	/* Hopefully irrelevant... */
#endif

#ifndef MAX_SLENGTH
#define MAX_SLENGTH	(1024)	/* Do NOT use this to override EXODUSII's
				 * preference of MAX_STR_LENGTH or face
				 * severe pain. You were warned. */
#endif

#ifndef FILENAME_MAX_ACK      /* Avert name clash with HP's stdio.h ... */
#define FILENAME_MAX_ACK      (1024)
#endif

#define EXODUS_FILE_SUFFIX (".exoII")

/*
 * States of the Exo_DB database in memory - where's it at?
 * 
 * State must be unique. Some states are exclusive of other states, some
 * states do not care about other states...
 *
 * Someday this will all be truly object oriented...
 *
 */

#ifndef EXODB_STATE_GRND
#define EXODB_STATE_GRND	(1L<<0)	/* It just exists, nothing's done. */
#define EXODB_STATE_INIT	(1L<<1)	/* initialization info exists. */
#define EXODB_STATE_MESH	(1L<<2)	/* also the mesh info. */
#define EXODB_STATE_RES0	(1L<<3)	/* Preliminary result information. */
#define EXODB_STATE_NDVA	(1L<<4) /* Nodal results var space allocated.*/
#define EXODB_STATE_NDVR	(1L<<5) /* Nontrivial nodal results data read*/
#define EXODB_STATE_ELVA	(1L<<6) /* Elem results var space allocated.*/
#define EXODB_STATE_ELVR	(1L<<7) /* Nontrivial elem results data read*/
#define EXODB_STATE_GBVA	(1L<<8) /* Glob results var space allocated.*/
#define EXODB_STATE_GBVR	(1L<<9) /* Nontrivial glob results data read*/
#define EXODB_STATE_NDIA	(1L<<10) /* Indeces allocated for node vars. */
#define EXODB_STATE_ELIA	(1L<<11) /* Indeces allocated for elem vars. */
#define EXODB_STATE_GBIA	(1L<<12) /* Indeces allocated for glob vars. */
#endif

/*
 * Actions to be performed - "What is to be done?" - V. I. Lenin
 *
 * Ultimately, if some of the requested actions are compared with various
 * of the states of existence, unnecessary or inappropriate actions may
 * be flagged and handled accordingly.
 */

#ifndef EXODB_ACTION_RD_INIT
#define EXODB_ACTION_RD_INIT	(1L<<0)	/* Read the mesh meta data. */
#define EXODB_ACTION_RD_MESH	(1L<<1)	/* Read the mesh connectivity, etc. */
#define EXODB_ACTION_RD_RES0	(1L<<2)	/* Read results meta data. */
#define EXODB_ACTION_RD_RESN	(1L<<3)	/* Read results data big time (node) */
#define EXODB_ACTION_RD_RESE	(1L<<4)	/* Read results data big time (elem) */
#define EXODB_ACTION_RD_RESG	(1L<<5)	/* Read results data big time (glob) */
#define EXODB_ACTION_WR_INIT	(1L<<6)
#define EXODB_ACTION_WR_MESH	(1L<<7)
#define EXODB_ACTION_WR_RES0	(1L<<8)
#define EXODB_ACTION_WR_RESN	(1L<<9)
#define EXODB_ACTION_WR_RESE	(1L<<10)
#define EXODB_ACTION_WR_RESG	(1L<<11)
#endif

/*
 * For routine assignment from the src to the dst Exo_DB structure members
 * we have a macro that is used extensively in copy_exo().
 */

#define SRC_DST(member)		dst->member = src->member

/*
 * Note that some variable names are replicated no longer!
 */

struct Exodus_Database
{
  long  state;			/* of the database in memory */

  char	*path;			/* filename. Well, yes this could
				 * change, too, for reasons beyond our
				 * control.
				 */
  int	exoid;			/* netCDF open file identifier. Strictly,
				 * this only has meaning when the file is
				 * open and is not part of the database.
				 * In principle, this netCDF ID can be
				 * assigned to an arbitrary integer 
				 * that could change from call to call. 
				 */

  int	comp_wordsize;		/* native desired size for reals on machine */
  int	io_wordsize;		/* size of reals in database (maybe difft) */

  /*
   * General model information...
   */

  int   cmode;
  int	mode;
  int	num_dim;		/* Number of spatial dimensions */

  int	num_elems;		/* Number of elements, total. */

  int	num_elem_blocks;	/* Number of element blocks. */

  int	num_info;		/* Number of info[] records. */

  int	num_node_sets;		/* Number of node sets. */

  int	num_side_sets;		/* Number of side sets. */

  int	num_nodes;		/* Number of nodes, total. */

  int	num_qa_rec;		/* Number of quality assurance records. */

  int	num_times;		/* Number of time planes. */

  int   node_map_exists;
  int   elem_map_exists;
  int   elem_order_map_exists;
  int   ss_node_list_exists;

  int	*node_map;		/* returned by ex_get_node_num_map() */
  int	*elem_map;		/* returned by ex_get_elem_num_map() */
  int	*elem_order_map;	/* returned by ex_get_map() */

  dbl	*x_coord;
  dbl	*y_coord;
  dbl	*z_coord;

  flt	api_version;
  flt	db_version;
  flt	version;

  QA_Record *qa_record;

  /*   char	*qa_record[MAX_QA][4];  deadly poison! */

  INFO_Record	*info;

  /*  char	**info; */

  char	**coord_names;

  char	*title;			/* Title of the database. */

  /*
   * Element block data...
   */

  int	*eb_id;			/* Element block identifiers. */
  char	**eb_elem_type;		/* Type of element for a particular block. */
  int	*eb_num_elems;		/* Number of elements in ea block. */
  int	*eb_num_nodes_per_elem;	/* Number nodes/element in this block. */
  int	*eb_num_attr;		/* Number of attributes/elem in this block. */
  int	**eb_conn;		/* Connectivity in this eb. */
  dbl	**eb_attr;		/* Attributes in this eb. */

  /*
   * These are extras that are merely a convenience.
   */

  int   *eb_ptr;		/* [neb+1] - ptr sums up eb_num_elems[ieb] */
  int   *eb_elem_itype;		/* an integer instead of a character string */

  int   *elem_ptr;		/* [num_elems+1] - pt into node list */
  int   *node_list;		/* concatenated eb_conn for all elemblocks */

  /*
   * Additional element -> node connectivity information. Boolean indicates
   * arrays allocated and filled with meaningful data.
   *
   * These elem_node connectivities are just renamed versions of the elem_ptr
   * and node_list arrays above. The more descriptive name might help.
   */

  int    elem_node_conn_exists;
  int   *elem_node_pntr;
  int   *elem_node_list;

  /*
   * Additional node -> element connectivity information. Boolean indicates
   * arrays allocated and filled with meaningful data.
   */

  int    node_elem_conn_exists;
  int   *node_elem_pntr;
  int   *node_elem_list;

  /*
   * Additional node -> node connectivity information. 
   *
   * Boolean indicates arrays allocated and filled with meaningful data.
   *
   * These are nodes that are connected to other nodes by no more than one
   * element.
   *
   * This includes the self interaction - a node is connected to itself.
   *
   * The list of nodes to which a given node is connected is sorted in
   * ascending order. This makes it easy to implement a quick check
   * for whether an arbitrary given node is in the list - if its name is less
   * than the first name or greater than the last name you don't need to
   * check the other names.
   */
 
  int    node_node_conn_exists;
  int   *node_node_pntr;
  int   *node_node_list;

  /*
   * Additional element -> element connectivity information. Boolean indicates
   * arrays allocated and filled with meaningful data. These are elements
   * connected by faces of one less dimension. Thus, quadrilaterals connected
   * by line segments, or hexahedrons connected by quadrilaterals are valid
   * members of this connectivity, but quadrilaterals sharing but one point
   * or hexahedrons sharing but one line segment are invalid members of this
   * list.
   */

  int    elem_elem_conn_exists;
  int   *elem_elem_pntr;
  int   *elem_elem_list;	/* -2 means another proc - Look in dpi... */
  int   *elem_elem_twst;	/* How many twists? (Mainly for 3D elems.)*/
  int   *elem_elem_face;	/* Name of neighbor's face I am connected to.*/
  /*
   * Node set information...
   */

  int	ns_node_len;		/* Length of nodelist forall NS. */
  int	ns_distfact_len;	/* Length of df list forall NS. */
  int	*ns_id;			/* Node set IDs. */
  int   *ns_num_nodes;		/* Number of nodes in each ns. */
  int	*ns_num_distfacts;	/* Number of dfs in each ns.  */
  int	*ns_node_index;		/* Index in big list of nodes for ea ns.  */
  int	*ns_distfact_index;	/* Index in big list of dfs for ea ns. */
  int	*ns_node_list;		/* Big list of nds for all ns. */
  dbl	*ns_distfact_list;	/* Big list of dfs for all ns. */

  /*
   * Side set information...
   */

  int	ss_elem_len;		/* Length of ss element/side lists. */
  int	ss_distfact_len;	/* Length of ss df list */
  int	ss_node_len;		/* Length of ss node list. */

  int	*ss_id;			/* SS identifiers for ea side set. */
  int	*ss_num_sides;		/* SS num sides per set. */
  int	*ss_num_distfacts;	/* SS num of dfs per set. */
  int	*ss_elem_index;		/* SS index into element list. */
  int	*ss_distfact_index;	/* SS index into df list. */
  int	*ss_elem_list;		/* SS element list. */
  int	*ss_side_list;		/* SS side list. */

  dbl	*ss_distfact_list;	/* SS df list. */

  /*
   * These weird beasts help to reference node names and the correct
   * distribution factors in a sideset. This is problematic because the
   * distribution factors are associated with nodes in the sideset, not
   * sides, as you might expect.
   *
   * The arrays ss_node_cnt_list and ss_node_list are doubly indexed. First,
   * there is one for each side set. Second, they are indexed according to
   * the element/side counter within each sideset.
   */

  int   **ss_node_cnt_list;

  int   **ss_node_list;

  int   **ss_node_side_index;	/* For my own convenience... */


  /*
   * Properties...
   */
  
  int	ns_num_props;		/* Number of nodes set properites. */
  int	ss_num_props;		/* Number of side set properties. */
  int	eb_num_props;		/* Number of element block properties. */

  char	**ns_prop_name;		/* Names of node set properties. */
  char	**ss_prop_name;		/* Names of side set properties. */
  char	**eb_prop_name;		/* Names of element block properties. */

  int	**ns_prop;		/* Values of node set properties. */
  int	**ss_prop;		/* Values of side set properties. */
  int	**eb_prop;		/* Values of element block properties. */

  /*
   * Preliminary results data (meta-data) or about the results.
   */

  int	num_glob_vars;
  int	num_node_vars;
  int	num_elem_vars;

  char	**glob_var_names;
  char	**elem_var_names;
  char	**node_var_names;

  dbl	*time_vals;

  int	*elem_var_tab;		/* [neb*nev], [blk_index*nev+ev_index] */
  int   elem_var_tab_exists;
  int   *truth_table_existance_key; /* programurz cant rite or spel */

  /*
   * Nodal variables are chosen by index, by timeplane and by which node
   * we're at.
   */

  int	num_nv_time_indeces;
  int	*nv_time_indeces;
  int	num_nv_indeces;
  int	*nv_indeces;

  dbl	***nv;			/* nv[time_index][nv_index][node] */

  /*
   * Element variables are chosen by index, by timeplane and by element block
   * ID as well as which element we're in.
   */

  int	num_ev_time_indeces;
  int	*ev_time_indeces;

  /*
   * Use the element variable truth table to select which element variables
   * will be printed for which blocks.
   */
  
   int	num_ev_indeces;
   int   *ev_indeces;
   
   /*
   *  int   num_ev_blk_ids;
   *  int	*ev_blk_ids;
   */

  /*  dbl   ****ev;		 ev[time_index][blk_index][ev_index][elem] */

  dbl   ***ev;			/* ev[time_index][blk_index,ev_index][elem] */

  /*

   *ev == elem_var_vals; 

   */

  /*
   * Global variables are chosen by the timeplane - you get *all* of them -
   * no cherry picking of global variables.
   */

  int	num_gv_time_indeces;
  int	*gv_time_indeces;

  dbl	**gv;			/* gv[time_index][gv_index] */

  /*
   * Historical compatibility - these become arrays of pointers into
   * chunks for some given SINGLE timeplane. They are still quite useful
   * in cases where memory constraints limit the number of timeplanes of
   * results data that may be stored.
   *
   * Note the varying degrees of indirection...
   *
   * Note, too, the ordering convention for these indeces. The element
   * variables have the element block index appearing foremost, with
   * the element variable index appearing secondarily and the element
   * counter within the block appearing last. The basic idea is to allocate
   * the pointers pretty much unconditionally, with the final allocation
   * of the number of elements in a block being contingent on the
   * truth table. Thus, the 2d pointer array represented by elem_var_vals
   * can be "holey".
   */

  /*
   * These old variables are convenient for doing "just one timeplane" of
   * results.
   */

  dbl	**node_var_vals;	/* node_var_vals[nv_index][node] */
  dbl   *glob_var_vals;		/* glob_var_vals[gv_index] */
  dbl   **elem_var_vals;	/* elem_var_vals[block_ev_index][elem]*/
};

typedef struct Exodus_Database Exo_DB;

#endif /* _EXO_STRUCT_H */
