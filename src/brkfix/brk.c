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

/* brk/fix -- problem de/re-composition of EXODUS II database + goma probdesc
 *
 *
 *	brk -- examines a brk input file describing multiphysics interaction
 *             and a monolithic EXODUS II file. A problem graph is constructed,
 *             sent to Chaco to be properly fragmented and the resulting
 *             pieces are used to build the proper EXODUS II files, augmented
 *             with distributed processing information stored as additional
 *             netCDF variables.
 *
 *	fix -- examine a polylithic set of EXODUS II files and reconstruct
 *	       a monolithic EXODUS II file.
 *
 * Output:
 *		Create n little EXODUS IIv2 files. 
 *
 *		These little EXODUS IIv2 files have a new, local
 *		numbering of elements and nodes and also a map to the unique
 *		element and node names used for the global problem.
 *
 *		The little EXODUS IIv2 files will be AUGMENTED with
 *		communication variables stored using the netCDF format.
 *		These indicate which processors to send messages to, 
 *		to receive from, the global and local names of the dofs
 *		that are being sent and received.
 *		
 *		All of the dofs at any finite element node will belong to 
 *		a single processor. However, every processor will have some
 *		"external" nodes that are the primary responsibility of another
 *		processor.
 *		
 *		In order to properly assemble all contributions at a node that
 *		are needed, processors will need to loop over every element
 *		that contains each node for which it is responsible. That
 *		means that elements near the interprocessor boundary
 *		will appear in the element lists for more than one processor.
 *
 *		Lists of global degree of freedom names will be built on 
 *		each processor, including the corresponding local dof names.
 *		Each list will correspond to a recv or dest list and the
 *		name of the processor that will send or recv.
 *
 *		Assume that degrees of freedom either exist or not based solely
 *		upon the variable id. Consequences of this assumption are that
 *		if two element blocks adjoin, that the velocity has the same
 *		number of vector components and is interpolated so that the
 *		number of nodal degrees of freedom for each component is the
 *		same in each element block.
 *		
 *		Likewise, there are either zero or N concentrations in any
 *		element block.
 *
 * To Do:
 *		After bulk element block assignments of eqnvars and 
 *		interactivity between equations and variables, need to
 *		sweep over specified boundaries for lower dimensional
 *		specifications, like Randy's nodal dof jump for specified
 *		discontinous variables and like Rich's subparametric trick.
 *
 *		Also, the interactions of the boundary conditions themselves
 *		need to be expressed for the most accuracy in load balancing.
 *
 *		Provisions for dynamic balancing. How to perform an exchange
 *		and how does the application code do it? That's the crux of
 *		the matter - how does the application adapt to changes in
 *		the mesh ? Any initial mesh should do fine, but the code and
 *		the data structures should lend themselves to change.
 *             
 * 
 * Created:  1997/04/12 08:14 MDT pasacki@sandia.gov
 *
 * Revised:  1998/01/20 15:06 MST pasacki@sandia.gov
 */

#define _BRK_C

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _INCLUDE_POSIX_SOURCE	/* needed for HP-UX */

#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>		/* needed on HP-UX */

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <unistd.h>

#include <string.h>

#include "map_names.h"
#include "std.h"		/* useful general stuff */
#include "eh.h"			/* error handling */
#include "aalloc.h"		/* multi-dim array allocation */
#include "exo_struct.h"		/* some definitions for EXODUS II */
#include "dpi.h"		/* distributed processing information */
#include "nodesc.h"		/* node descriptions */
#include "brkfix_types.h"
#include "wr_graph_file.h"
#include "rd_in.h"
#include "wr_coords.h"
#include "ppi.h"
#include "mk_dm.h"
#include "sam_perea.h"
#include "utils.h"
#include "exo_utils.h"
#include "exo_conn.h"
#include "emuck.h"
#include "rd_dpi.h"
#include "rd_exo.h"

/*
 * The general dependency matrix exists for each element block.
 * It tells how the Basic constraint equations depend on corresponding
 * Basic variables.
 *
 *	evd[element_block][equation][variable]
 *
 *	   [NUM_ELEM_BLOCKS][EQNS_EB][EQNS_EB]
 *
 * For the interaction between equations and variables at nodes through
 * any element that is a member of element_block, this array tells the
 * strength of the interaction.
 *
 *  e.g., evd[0][TEMPERATURE][MESH_DISPLACEMENT] = 1
 *        evd[0][CONTINUITY][PRESSURE]           = 0
 *
 * While entries like "0" and "1" can denote the absence or existence
 * of a dependence, larger integer values may be used to ascribe greater
 * weight to the dependence. These larger values should translate
 * into greater weights for the edges and verteces according to the
 * magnitude of the dependence.
 */

int	***evd;

/*
 * Lucky -- the local node dof existence profile.
 * 
 * For the typical element, as we look at each node, different variables
 * can be active (-1 if no vbls active).
 *
 *	Lucky[element_block][local_node_number][eqnvar_index] = 0 or 1
 *
 * Note that these are PROVISIONAL! From the global perspective, 
 * nodes at the boundary between eb's will have different existence 
 * profiles. Expect that a careful search of nodes on inter-elementblock 
 * boundaries will reveal the need for different profiles than are in
 * either bounding block.
 *
 */

int ***Lucky;

static int show_help			= FALSE;
static int add_decomp_plot_vars		= FALSE;
static int be_quiet			= FALSE;
static int be_verbose			= FALSE;
static int preprocess_input_file	= TRUE;
static int rescale_edge_weights		= FALSE;
static int rescale_vertex_weights	= TRUE;
static int one_neighbor_per_line	= FALSE;
static int write_monolith_graph		= FALSE;
static int write_monolith_coords	= TRUE;

/*
 * Used to compute Sam's heuristic...
 */

static int total_internal_dofweight = 0;
static int total_boundary_dofweight = 0;
static int total_external_dofweight = 0;

char *program_name;		/* name this program was run with */

static char program_version[] = BRK_VERSION;

const char program_description[] = "GOMA distributed problem decomposition tool";

const char copyright[]="Copyright (c) 1999-2000 Sandia National Laboratories. All rights reserved.";

#define LAST_LEGAL_STRING "Last Legal String"

static char user_params_filename[] = "User_Params";

static char *chaco_user_params_file[] =
{
  "% Chaco 2.0 input file ",
  "%",
  "% This \"User_Params\" file was created by brk for decomposing a GOMA",
  "% finite element problem. If \"User_Params\" already exists, then brk",
  "% will not overwrite it.",
  "%",
  "% Chaco has MANY options. For a description of those options and",
  "% how they can be tailored to your problem, cf SAND95-2344, pp. 24--38",
  "%",
  "%",
  "check_input     = false",
  "echo            = 0",
  "output_metrics  = 1",
  "output_time     = 0",
  "output_assign   = true",
  "out_assign_inv  = true",
  "in_assign_inv   = false",
  "prompt          = false",
  "print_headers   = false",
  "architecture    = 1",
  LAST_LEGAL_STRING
};

/*
 * Make an array of pointers to structures like this. 
 * First index for the element block,
 * second index for the eqnvar we're dealing with. Note that the eqnvar
 * index is not the same as the eqnvar ID.
 */

Bevm ***mult; 

int *ep;			/* element pointers into node list */
int *np;			/* node pointers into element list */

int *nl;			/* node list */
int *el;			/* element list */

int *ebl;			/* element block list */

/*
 * Each node can be associated with one of a handful of prototypes.
 * Each prototype is extensively described by the number and names of
 * the active eqnvars at the node, as well as the cumulative weights
 * that apply to each eqnvar.
 */

Node_Description **pnd;	/* main list of prototype node descriptions */
Node_Description *tnd;	/* temporary pointer */
Node_Description *end;	/* for the eqnnode d'jour */
Node_Description *vnd;	/* for the varnode d'jour */

#define CHACO
#ifdef CHACO
extern int interface
PROTO((int,			/* nvtxs - number of vertices in full graph */
       int *,			/* start - start of edge list for ea vertex */
       int *,			/* adjacency - edge list data */
       int *,			/* vwgts - weights for all vertices */
       flt *,			/* ewgts - weights for all edges */
       flt *,			/* x coordinates for inertial method */
       flt *,			/* y coordinates for inertial method */
       flt *,			/* z coordinates for inertial method */
       char *,			/* outassignname - assignment out filename */
       char *,			/* outfilename - output file name */
       int *,			/* assignment - set num each vtx (length n) */
       int,			/* architecture - (0=hypercube, d=dim mesh) */
       int,			/* ndims_tot - number cube dims to divide */
       int [],			/* mesh_dims[3] - dims of mesh of processors */
       dbl *,			/* goal - desired set sizes for each set */
       int,			/* global_method - global partition alg */
       int,			/* local_method - local partition algorithm */
       int,                     /* rqi_flag - use RQI/Symmlq eigensolver? */
       int,			/* vmax - if so, coarsen down to ? vertices */
       int,			/* ndims - number of eigenvectors (2^d sets) */
       dbl,			/* eigtol - tolerance on eigenvectors */
       long));			/* seed - for random graph mutations */ 
#endif

extern int wr_mesh_exo		/* wr_exo.c */
PROTO((Exo_DB *,		/* exo - ptr to full ripe EXODUS II fe db */
       char *,			/* filename - where to write */
       int));			/* verbosity - talk while writing */

extern void wr_resetup_exo	/* wr_exo.c */
PROTO((Exo_DB *,		/* exo - ptr to full ripe EXODUS II fe db */
       char *,			/* filename - where to write */
       int ));			/* verbosity - 0 for quiet, more to talk */

extern void wr_result_exo	/* wr_exo.c */
PROTO((Exo_DB *,		/* exo */
       char *,			/* filename - where to write */
       int ));			/* verbosity - 0 for quiet, more to talk */

extern int wr_dpi		/* wr_dpi.c */
PROTO((Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* filename */
       int ));			/* verbosity - how much to talk */


static int integer_compare	/* used internally by qsort() brk.c */
PROTO((const void *, 
       const void *));


static void 
usage(const int status)
{
  if (status != 0)
    {
      fprintf(stderr, ("Try `%s -h' for help.\n"),
	       program_name);
    }
  else
    {
      fprintf(stdout, "Usage: %s [OPTIONS] -n p infile exofile\n", 
	      program_name);
      fprintf(stdout, "\n");
      fprintf(stdout, "Convert an EXODUS II finite element mesh and a GOMA multiphysics problem\n");
      fprintf(stdout, "description into an edge and vertex weighted graph, decompose the graph\n");
      fprintf(stdout, "using Chaco or Metis, and build little augmented EXODUS II files for\n");
      fprintf(stdout, "subsequent parallel processing simulations with GOMA.\n");
      fprintf(stdout, "\n");
      fprintf(stdout, "OPTIONS:\n");
      fprintf(stdout, "\t-1        format graph file 1 neighbor per line\n");
      fprintf(stdout, "\t-a file   created monolith w/ added decomp. plot vars\n");
      fprintf(stdout, "\t-c file   write nodal coordinates to specified file\n");
      /* fprintf(stdout, "\t-e file   break this external file, too\n"); */
      fprintf(stdout, "\t-g file   write graph to specified file\n");
      fprintf(stdout, "\t-h        summarize OPTIONS (i.e., this message)\n");
      fprintf(stdout, "\t-n P      break into P pieces\n");
      fprintf(stdout, "\t-p        do not run preprocessor on input file\n");
      /* fprintf(stdout, "\t-q        quiet mode\n"); */
      /* fprintf(stdout, "\t-s        re scale vertex and edge weights\n"); */
      fprintf(stdout, "\t-v        print version \n");
      fprintf(stdout, "\n");
    }
  exit(-1);
}

int 
main (int argc, char *argv[], char *envp[])
{
  int a;
  int at_least_one;

  int blk_index;
  int begin;
  int begin_df;
  int *boundary_nodes;
  int internal_dofweight, boundary_dofweight, external_dofweight;

  int c;			/* hold each option flag */
  int *ccs_contribute;		/* contribution to communication */
  int count;

  int *dragon;

  int numerator;
  int denominator;

  int e;
  int e_name;
  int eb;			/* element block counter */
  int eb_id;
  int eb_index;
  int *ebi;
  int edge_weight;
  int elem;
  int *element_procs;		/* list unique procs that assemble an elem */
  int *element_owner;
  int *element_owner_dist;
  int *element_bnd_stat;
  int enode;
  int *eqn_node_names;		/* for each node-node interaction */
  int err;
  int *external_nodes;

  int g_ns_index;
  int g_ss_index;
  int g_eb_index;
  int gbeg;
  int ge;
  int glob_eb_index;
  int global_node_index;
  int gnattr;
  int gndfpn;
  int gnnns;

  int I;
  int i;
  int ieb;
  int index;
  int index_boundary_node;
  int index_ccs;
  int index_dofs;
  int index_eoa;
  int index_eown;
  int index_external_node;
  int index_internal_node;
  int index_max_private_elems;
  int index_nat;
  int index_nknd;
  int index_nnz;
  int index_nown;
  int index_ns;
  int index_nscn;
  int index_priv;
  int index_proc_elems;
  int index_save;
  int index_shar;
  int index_ss;
  int inn;
  int ins;
  int *internal_nodes;
  int iss;

  int j;
  int jeb;

  int k;
  int keep_going;

  int l;
  int len;
  int len_sm;
  int length;
  int length_node_node;
  int lnattr;
  int local_elem;
  
  int m;
  int max_basic_eqnvars;
  int max_eweight;
  int max_private_elems;
  int max_vweight;
  int min_eweight;		/* of graph edges; rescale with this! */
  int min_vweight;		/* of graph verteces; rescale with this! */
  int my_set;
  
  int n;			/* node counter */
  /*   int nargs;		 number of nonoption arguments leftover */
  int *nat_contribute;		/* number of assembled terms for nn interact. */
  int ne;			/* number of elements */
  int neb;			/* number of element blocks */
  int nev;			/* number of element variables */
  int neighbor_set_name;
  int *new_proc_elems;
  int *new_proc_eb_ptr;
  int new_start;
  int *new_truth_table;		/* element variables get expanded for plotting
				 * of distribution */
  int new_weight;
  int nn;			/* number of nodes */
  int nns;			/* number of node sets */
  int nnv;			/* number of nodal variables */
  int *nnz_contribute;		/* to nonsparsity for a node-node interact. */
  int node;
  int *node_kind;
  int *node_dof0;
  int ns_touch;
  int nss;			/* number of side sets */
  int nsets;
  int *num_basic_eqnvars;	/* how many basic equation/variable pairs
				 * for each element block, indexed using the
				 * element block index, not the element block
				 * ID, which are not necessarily contiguous
				 * integers.
				 */

  int num_edges;		/* of the graph */
  int num_boundary_nodes;
  int num_element_procs;	/* count procs assembling an element */
  int num_external_nodes;
  int num_internal_nodes;
  int num_universe_nodes;
  int num_kinds_nodes;		/* which classic prototype */
  int num_pieces=1;		/* in which to divide the monolith */
  int num_sets=0;
  int num_recv_procs;
  int num_send_procs;
  int num_verteces;		/* of the graph */

  int old_start;
  int owner;			/* proc/set id of an element */

  int p;
  int *pnn;			/* pointer into node-node connectivity list */
  int private_elem;
  int *private_elem_count;
  int *proc_eb_id;
  int proc_eb_index;

  int *proc_eb_ptr;		/* frequency or offset */

  int *proc_elem_priv;
  int *proc_elem_shar;

  int *proc_elems;
  int proc_neb;
  int *proc_nodes;
  int *proc_ns_id;
  int *proc_ns_num_nodes;
  int *proc_ns_num_distfacts;
  int *proc_ns_node_index;
  int *proc_ns_distfact_index;
  int *proc_ns_node_list;
  int *proc_ns_node_list_index_global; /* New! */
  int proc_ne;
  int proc_num_node_sets;
  int proc_num_side_sets;
  int proc_ns_node_len;
  int proc_ns_distfact_len;
  int *proc_ns_distfact_list_index_global;
  int *proc_psm;		/* per proc version of ptr to set membership */
  int *proc_sm;			/* set membership list (of me & my neighbors) */
  int *proc_ss_id;
  int *proc_ss_num_sides;
  int *proc_ss_num_distfacts;
  int *proc_ss_elem_index;
  int *proc_ss_elem_list_index_global;
  int *proc_ss_distfact_index;
  int *proc_ss_distfact_list_index_global;
  int *proc_ss_elem_list;
  int proc_ss_elem_len;
  int proc_ss_distfact_len;
  int *proc_ss_side_list;
  int *psm;			/* ptrs to set memberships */

  int *recv_proc_names;		/* from the standpoint of a proc */

  int s;
  int scale_eweight;
  int scale_vweight;
  int *send_proc_names;
  int setalogue[MAX_ADJOINING_SETS];
  int set_name;
  int side;
  int size_proc_num_ns;
  int size_proc_ns_node_list;
  int size_proc_ns_distfact_list;
  int size_proc_num_ss;
  int size_proc_ss_elem_list;
  int size_proc_ss_distfact_list;
  int *sm;			/* set membership */
  int start_nd;			/* start of nodes in the ns_nodelist */
  int start_el;			/* start of elements in the ss_elemlist */
  int start_df;			/* start of distfacts in the ns_distfactlist */
  int status  = 0;		/* in case anything goes wrong */
  int sum_ccs;
  int sum_nat;
  int sum_nnz;
  
  int t;
  int tnnz;
  int tnat;
  int total_dofs;
  int touch;

  int *var_node_names;		/* for each node-node interaction */
  int vertex_weight;
  int vnode;

  int weight_e_v;
  int weight_v_e = 0;
  int where;

  int *assignment;

  dbl *proc_ns_distfact_list;
  dbl *proc_ss_distfact_list;

  Dpi *D;			/* distributed processing information */

  Exo_DB *mono;			/* monolithic EXODUS IIv2 database */

  Exo_DB *E;			/* for each polylithic piece */

  char **ptmp;

  struct stat stat_buf_struct;

  char  in_file_name[FILENAME_MAX_ACK];	/* to read in specifications */

  char  in_exodus_file_name[FILENAME_MAX_ACK]; /* input monolith EXODUS II file */

  char  out_coord_file_name[FILENAME_MAX_ACK]; /* coordinate file name */
  
  char  out_extra_file_name[FILENAME_MAX_ACK]; /* external field variables */

  char  out_graph_file_name[FILENAME_MAX_ACK]; /* graph file name */

  char  out_augplot_file_name[FILENAME_MAX_ACK]; /* decomposition plot file name */

  char  err_msg[MAX_CHAR_ERR_MSG];

  char *tmp;			/* char pointer junkyard of no interest */

  extern char *optarg;
  extern int   optind;

  Spfrtn sr=0;

#ifdef CHACO
  int nvtxs;
  int *start;
  int *adjacency;
  int *vwgts;
  flt *ewgts;
  flt *x = NULL;
  flt *y = NULL;
  flt *z = NULL;

  char *outassignname = NULL;
  char *outfilename = NULL;


  int architecture;
  int ndims_tot;
  int mesh_dims[3];
  
  dbl *goal = NULL;

  int global_method;
  int local_method;

  int rqi_flag;

  int user_params_file_exists;

  FILE *fs_up;			/* file stream for User_Params */

  int vmax;
  int ndims;

  dbl eigtol;
  long seed;

#endif /* CHACO */

  program_name = argv[0];

  /*
   * Defaults
   */

  tmp = strcpy(in_file_name,        "in");
  
  tmp = strcpy(in_exodus_file_name, "in.exoII");

  tmp = strcpy(out_coord_file_name, "coords");

  tmp = strcpy(out_graph_file_name, "graph");

  tmp = strcpy(out_augplot_file_name, "brk.exoII");

  /*
   * Process any command line options.
   */

  while ( ( c = getopt(argc, argv, "1a:c:e:g:hn:pqsv")) != EOF )
    {
      switch(c)
	{
	case '1':
	  one_neighbor_per_line = TRUE;
	  break;

	case 'a':
	  add_decomp_plot_vars = TRUE;
	  tmp = strcpy(out_augplot_file_name, optarg);
	  break;

	case 'c':
	  /*
	   * Name of coordinate file.
	   */
	  tmp = strcpy(out_coord_file_name, optarg);
	  break;

	case 'e':
	  /*
	   * Name of external data file - not yet implemented!
	   */
	  tmp = strcpy(out_extra_file_name, optarg);
	  break;

	case 'g':
	  /*
	   * Name of graph file.
	   */
	  write_monolith_graph = TRUE;
	  tmp = strcpy(out_graph_file_name, optarg);
	  break;

	case 'h':
	  /*
	   * Show synopsis of help.
	   */
	  show_help |= TRUE;
	  usage(status);
	  break;

	case 'n':
	  err = sscanf(optarg, "%d", &num_pieces);
#ifdef DEBUG
	  fprintf(stderr, "Processing with optind=%d, optarg=%s, argc=%d\n",
		  optind, optarg, argc);
#endif
	  if ( err != 1 || num_pieces < 1 )
	    {
	      EH(-1, 
		 "For number of pieces, specify one positive nonzero integer.");
	    }
	  break;

	case 'p':
	  preprocess_input_file = FALSE;
	  break;

	case 'q':
	  /*
	   * Quiet operation mode.
	   */
	  be_quiet |= TRUE;
	  if ( be_verbose )
	    {
	      show_help |= TRUE;
	      status = -1;
	    }
	  break;

	case 's':
	  rescale_edge_weights   = TRUE;
	  rescale_vertex_weights = TRUE;
	  break;

	case 'v':

	  /*
	   * Print code version and exit.
	   */

	  fprintf(stdout, "%s\n", BRK_VERSION);
	  exit(0);

#if 0
	  be_verbose |= TRUE;
	  if ( be_quiet )
	    {
	      show_help |= TRUE;
	      status = -1;
	    }
#endif
	  break;

	default:
	  show_help = TRUE;
	  status = -1;
	  break;
	}
    }

  if ( optind >= argc )
    {
      usage(-1);
    }

#ifdef DEBUG
  fprintf(stderr, "Opts done: optind=%d, optarg=%s, argc=%d, argv[%d]=%s\n",
	  optind, optarg, argc, optind, argv[optind]);
#endif

  tmp = strcpy(in_file_name, argv[optind]);
  optind++;

#ifdef DEBUG
  fprintf(stderr, "Opts done: optind=%d, optarg=%s, argc=%d, argv[%d]=%s\n",
	  optind, optarg, argc, optind, argv[optind]);
#endif
  tmp = strcpy(in_exodus_file_name, argv[optind]);
  optind++;

  /*
   * Save a default User_Params file for Chaco if the user has not already
   * created one.
   */

  user_params_file_exists = FALSE;
  err = stat(user_params_filename, &stat_buf_struct);

  /*
   * A normal return suggests the file already exists.
   */

  user_params_file_exists = ( err == 0 );

  /*
   * An abnormal return could be caused by many things. The file not existing
   * is what we're concerned with.
   */

  if ( err == -1 )
    {
      user_params_file_exists = ! ( errno == ENOENT );
    }

  if ( ! user_params_file_exists )
    {
      fs_up = fopen(user_params_filename, "w");
      if ( fs_up == NULL )
	{
	  sr = sprintf(err_msg, 
	  "Problem opening \"%s\" [Chaco input] for create/append.",
		       user_params_filename);
	  EH(-1, err_msg);
	}
      ptmp  = chaco_user_params_file;
      while ( strcmp(*ptmp, LAST_LEGAL_STRING) != 0 )
	{
	  fprintf(fs_up, "%s\n", *ptmp++);
	}
      fclose(fs_up);
    }

  if ( num_pieces < 2 )
    {
      sr = sprintf(err_msg, "Divide into more pieces than %d", num_pieces);
      EH(-1, err_msg);
    }

  if ( show_help )
    {
      usage(status);
    }
  
  if ( ! be_quiet )
    {
      fprintf(stdout, "%s %s - %s\n%s\n", program_name, program_version, 
	      program_description, copyright);
    }

  /*
   * Filter the input file into a temporary file whose name replaces the 
   * original.
   */

  if ( preprocess_input_file )
    {
      pre_process(in_file_name);
    }

#ifdef DEBUG
  /*  err = sscanf(line, "%s", in_exodus_file_name); */
  fprintf(stderr, "monolith EXODUSII file = \"%s\"\n", in_exodus_file_name);
#endif

  /*
   * Turn off annoying error reporting from within the EXODUS II API...
   */

  ex_opts(EX_VERBOSE); 

  /*
   * Make preliminary space for the monolithic database.
   */

  mono = (Exo_DB *) smalloc(sizeof(Exo_DB));
  init_exo_struct(mono);

  /*
   * Read everything except any results data. Do, however, read any
   * results metadata (that describes what it is and how much of it there is).
   *
   * WARNING! The RES0 means that nv_time_indeces[] etc will be allocated
   * and set.
   */

  status = rd_exo(mono, in_exodus_file_name, 0, (EXODB_ACTION_RD_INIT +
						 EXODB_ACTION_RD_MESH +
						 EXODB_ACTION_RD_RES0));
		  
  /*
   * Convenience variables...
   */

  ne  = mono->num_elems;
  nn  = mono->num_nodes;
  neb = mono->num_elem_blocks;
  nns = mono->num_node_sets;
  nss = mono->num_side_sets;

  /*
   * Push down the element names and node names by one so that C zero-based
   * arrays work more naturally. Fix up when outputting.
   */

  zero_base(mono);

  /*
   * The physical space coordinates of the finite element nodes corresponding
   * to each of the graph verteces are useful for inertial partitions that
   * are inexpensive, albeit of less quality than other partition algorithms.
   *
   * Write out those coordinates.
   */

  if ( write_monolith_coords )
    {
      write_coords(out_coord_file_name, mono);
    }

  /*
   * Read the input file for a better description of eqnvars etc than
   * is available in the EXODUS II file.
   */

  rd_input(in_file_name,	/* name of main input file           (in) */
	   mono,		/* EXODUS II fe database             (in) */
	   &mult,		/* multiplicity                     (out) */
	   &evd,		/* eqnvar dependencies              (out) */
	   &Lucky,		/* local nodedof existences         (out) */
	   &num_basic_eqnvars);	/* number of basic eqnvars in ea eb (out) */

  /*
   * Determine maximum number of basic eqns/vars in any element block.
   */

  max_basic_eqnvars = -999;

  for ( eb=0; eb<neb; eb++)
    {
      if ( num_basic_eqnvars[eb] > max_basic_eqnvars )
	{
	  max_basic_eqnvars = num_basic_eqnvars[eb];
	}
    }

  if ( max_basic_eqnvars > MAX_EQNVARS )
    {
      sr = sprintf(err_msg, "Try MAX_EQNVARS = %d\n", max_basic_eqnvars);
      EH(-1, err_msg);
    }

#ifdef DEBUG
  fprintf(stderr, "Maximum number of basic eqnvars in any block = %d\n",
	  max_basic_eqnvars);
#endif

  /*
   * Now these connectivities are created in the mono structure.
   */

  build_elem_node(mono);
  build_node_elem(mono);
  build_elem_elem(mono);
  build_node_node(mono);

  ep = mono->elem_node_pntr;
  nl = mono->elem_node_list;
  ebl = mono->eb_ptr;
  np = mono->node_elem_pntr;
  el = mono->node_elem_list;

  /*
   * Build the goma-like degree of freedom map for this problem.
   */

  node_kind = (int *) smalloc(nn*SZ_INT);
  INIT_IVEC(node_kind, -1, nn);

  node_dof0 = (int *) smalloc((nn+1)*SZ_INT);
  INIT_IVEC(node_dof0, -1, nn+1);

  pnd = (Node_Description **) smalloc(MAX_NODE_KINDS*
				      sizeof(Node_Description *));

  for ( i=0; i<MAX_NODE_KINDS; i++)
    {
      pnd[i] = (Node_Description *) smalloc(SZ_ND);
    }

  make_goma_dofmap(mono, mult, evd, Lucky, num_basic_eqnvars,	     /* (in) */
		   node_kind, node_dof0, pnd, &num_kinds_nodes);    /* (out) */


  total_dofs       = node_dof0[nn];	/* right after the last node */

  length_node_node = mono->node_node_pntr[mono->num_nodes];

  /*  length_node_node = count_node_node_interactions(nn, np, el, ep, nl);*/

#ifdef DEBUG
  fprintf(stderr, "Length node-node connectivity = %d\n", length_node_node);
#endif

  /*
   * Now allocate space for the 5 big arrays that will be used to assess
   * weights of graph verteces and graph edges for the problem. These are
   * filled in assess_weights().
   */

  eqn_node_names = (int *) smalloc(length_node_node*SZ_INT);
  var_node_names = (int *) smalloc(length_node_node*SZ_INT);
  nnz_contribute = (int *) smalloc(length_node_node*SZ_INT);
  nat_contribute = (int *) smalloc(length_node_node*SZ_INT);
  ccs_contribute = (int *) smalloc(length_node_node*SZ_INT);

  INIT_IVEC(eqn_node_names, UNDEFINED_EQNVARID, length_node_node);
  INIT_IVEC(var_node_names, UNDEFINED_EQNVARID, length_node_node);
  INIT_IVEC(nnz_contribute, 0,                  length_node_node);
  INIT_IVEC(nat_contribute, 0,                  length_node_node);
  INIT_IVEC(ccs_contribute, 0,                  length_node_node);
  
  /*
   * Number of nonzero matrix entries in the global matrix resulting
   * from this node-node interaction.
   *
   * Number of assembled terms resulting from this node-node interaction.
   *
   * Communication chunk size for this node-node interaction.
   */  

  assess_weights(mono, mult, evd, ebl, np, el, ep, nl, node_kind, /* (in) */
		 pnd, num_basic_eqnvars,			     /* (in) */
		 eqn_node_names, var_node_names, nnz_contribute,    /* (out) */
		 nat_contribute, ccs_contribute );                  /* (out) */

  /*
   * Create quick handy pointer into the node-node connectivity arrays.
   */

  pnn = (int *)smalloc((nn+1)*SZ_INT);

  index  = 0;
  for ( i=0; i<nn; i++)
    {
      pnn[i] = index;
      /*
       * Contortions to avoid any possible reference to eqn_node_names[dim+1].
       */
      keep_going = (index < length_node_node);
      if ( keep_going )
	{
	  keep_going &= ( eqn_node_names[index] == i );
	}
      while ( keep_going )
	{
	  index++;
	  keep_going = (index < length_node_node);
	  if ( keep_going )
	    {
	      keep_going &= ( eqn_node_names[index] == i );
	    }
	}
    }

  pnn[nn] = length_node_node;

  /*
   * The graph partitioner likes undirected graphs. Symmetrize the
   * communication cost vector to the maximum of each of the two cost estimates
   */

  for ( i=0; i<length_node_node; i++)
    {
      enode = eqn_node_names[i];
      vnode = var_node_names[i];
      weight_e_v = ccs_contribute[i];
#ifdef DEBUG      
      fprintf(stderr, "%d -> %d weighs %d\n", enode+1, vnode+1, weight_e_v);
#endif

      /*
       * Find the reciprocal interaction...
       */

      index_save = -1;
      for ( j=pnn[vnode]; j<pnn[vnode+1]; j++)
	{
	  if ( var_node_names[j] == enode )
	    {
	      weight_v_e = ccs_contribute[j];
	      index_save = j;
	    }
	}

#ifdef DEBUG      
      fprintf(stderr, "recipprocal: %d -> %d weighs %d\n", vnode+1, enode+1, weight_v_e);
#endif
      EH(index_save, "Could not find reciprocal edge..");

      new_weight = MAX(weight_e_v, weight_v_e);

      ccs_contribute[i]          = new_weight;
      ccs_contribute[index_save] = new_weight;
    }

#ifdef DEBUG
  fprintf(stderr, "Testing pnn:\n");
  for ( i=0; i<nn; i++)
    {
      fprintf(stderr, "Node (%d) -> ", i+1);
      for ( j=pnn[i]; j<pnn[i+1]; j++)
	{
	  fprintf(stderr, "%d ", var_node_names[j]+1);
	}
      fprintf(stderr, "\n");
    }
#endif

  tnnz = 0;
  tnat = 0;

#ifdef DEBUG
  fprintf(stderr, "Counting contributions to savings bonds and ECP...\n");
  fprintf(stderr, "tnnz starts at %d\n", tnnz);
  fprintf(stderr, "tnat starts at %d\n", tnat);
#endif

  for ( i=0; i<length_node_node; i++)
    {
      tnnz += nnz_contribute[i];
      tnat += nat_contribute[i];

#ifdef DEBUG
      fprintf(stderr, "%3d-%3d interaction contributes: %d %d -> %d <-\n",
	      eqn_node_names[i]+1, var_node_names[i]+1, nnz_contribute[i], 
	      nat_contribute[i], ccs_contribute[i]);
#endif
    }

#ifdef DEBUG
  err = fprintf(stdout, "nonzeroes = %d\n", tnnz);

  err = fprintf(stdout, "assembled terms = %d\n", tnat);
#endif


  /*
   * The graph itself is written out. The format is described on pp. 20-21
   * in SAND95-2344, "The Chaco User's Guide Version 2.0" by Hendrickson &
   * Leland. The format is very similar to the graph file formats accepted
   * by both Chaco and Metis.
   *
   * The number of verteces is likely the same as the number of nodes. 
   * However, sometimes we might have some lazy finite element nodes with
   * no connections to them.
   *
   * Lazy nodes with zero dofs are not handled well currently. While 
   * brk will zap out interactions where the edge weight is zero, if a node
   * has *no* eqns or dofs associated with it, then it can cause problems here.
   *
   *
   *    for ( i=0; i<length_node_node; i++)
   *   {
   *    index stuff...
   *   num_verteces += ( ) ? 1 : 0;
   *   }
   *
   *
   */

  num_verteces = nn;		/* i.e., the number of finite element nodes */

  /*
   * Count up the actual number of edges in the graph. Edges count
   * if they involve nontrival communication cost. Self interactions have
   * zero communication cost by assumption.
   */

  num_edges = 0;

  for ( i=0; i<length_node_node; i++)
    {
      num_edges += (ccs_contribute[i] != 0) ? 1 : 0;
    }

  /*  num_edges    = length_node_node - nn;  certainly fast. always right? */

  /*
   * Some heavily coupled problems have big vertex weights. Find the
   * lowest and highest vertex weights and, if desired, rescale weights 
   * so they are more reasonably sized.
   */

  max_vweight = -HUGE_INT;
  min_vweight = HUGE_INT;

  index = 0;
  for ( n=0; n<nn; n++)
    {
      vertex_weight = 0;
      for ( j=pnn[n]; j<pnn[n+1]; j++)
	{
	  vertex_weight += nat_contribute[j];
	}
      min_vweight = MIN(min_vweight, vertex_weight);
      max_vweight = MAX(max_vweight, vertex_weight);
    }

  if ( min_vweight < 1 )
    {
      fprintf(stderr, "Warning: Lazy nodes pulling no weight.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "         Consider remeshing with fewer nodes to\n");
      fprintf(stderr, "         avoid problems with effective noncontiguous\n");
      fprintf(stderr, "         vertex numbering choking graph partitioner.\n");
    }

  scale_vweight = gcf(MAX(1,min_vweight), max_vweight);

  if ( ! rescale_vertex_weights )
    {
      scale_vweight = 1;	/* despite all that work! */
    }


  /*
   * Likewise, some strongly coupled problems have large edge weights. Find
   * the lowest nontrival edge weight to represent the unit edge weight and 
   * rescale accordingly.
   */

  max_eweight = -HUGE_INT;
  min_eweight = HUGE_INT;

  for ( inn=0; inn<length_node_node; inn++)
    {
      edge_weight = ccs_contribute[inn];
      if ( edge_weight != 0 )
	{
	  min_eweight = MIN(min_eweight, edge_weight);
	}
      max_eweight = MAX(max_eweight, edge_weight);
    }

  /*
   * Use the greatest common factor of the min and max edge weights
   * as a scaling to make lowest meaningful integers. Intermediate
   * values could suffer some granularity due to integer division. Tough.
   */

  scale_eweight = gcf(min_eweight, max_eweight);

  if ( ! rescale_edge_weights )
    {
      scale_eweight = 1;	/* despite all that work! */
    }

  dragon = (int *) smalloc(LEN_DRAGON*SZ_INT);

  dragon[0] = num_verteces;
  dragon[1] = num_edges/2;

  dragon[2] = min_vweight;
  dragon[3] = scale_vweight;
  dragon[4] = max_vweight;

  dragon[5] = min_eweight;
  dragon[6] = scale_eweight;
  dragon[7] = max_eweight;

  if ( write_monolith_graph )
    {
      wr_graph_file(out_graph_file_name, in_file_name, in_exodus_file_name, 
		    ne, nn, total_dofs, tnnz, tnat, one_neighbor_per_line, pnn,
		    dragon, var_node_names, nat_contribute, ccs_contribute);
    }

#ifdef CHACO

  nvtxs        = num_verteces;

  start        = (int *) smalloc((nvtxs+1)*SZ_INT);

  adjacency    = (int *) smalloc(num_edges*SZ_INT);

  vwgts        = (int *) smalloc(nn*SZ_INT);

  ewgts        = (flt *) smalloc(num_edges*SZ_FLT);

  assignment   = (int *) smalloc(nvtxs*SZ_INT);

  /*
   * Some of these parameters will need to be read in from the inputfile.
   */

  architecture  = 1;		/* generic 1D mesh */

  ndims_tot     = 1;

  mesh_dims[0]  = num_pieces;
    
  global_method = 3;		/* tmp hardwire -- inertial */

  local_method  = 1;		/* tmp hardwire -- Kernighan-Lin */

  rqi_flag      = 1;

  vmax          = 0;

  ndims         = 1;

  eigtol        = 1e-3;

  seed          = 1234;

  /*
   * Fill in the coordinates of the graph verteces.
   */
  if ( mono->num_dim > 0 )
    {
      x = (flt *) smalloc(nn*SZ_FLT);
      for ( i=0; i<nn; i++)
	{
	  x[i] = mono->x_coord[i];
	}
    }
  if ( mono->num_dim > 1 )
    {
      y = (flt *) smalloc(nn*SZ_FLT);
      for ( i=0; i<nn; i++)
	{
	  y[i] = mono->y_coord[i];
	}
    }
  if ( mono->num_dim > 2 )
    {
      z = (flt *) smalloc(nn*SZ_FLT);
      for ( i=0; i<nn; i++)
	{
	  z[i] = mono->z_coord[i];
	}
    }


  /*
   * Fill in start pointer and vertex adjacency list.
   */

  index    = 0;
  start[0] = index;

  for ( n=0; n<nn; n++)
    {
      vwgts[n] = 0;

      for ( j=pnn[n]; j<pnn[n+1]; j++)
	{
	  vwgts[n] += nat_contribute[j];
	  if ( ccs_contribute[j] != 0 &&
	       var_node_names[j] != n ) 
	    {
	      adjacency[index] = var_node_names[j] + 1;
	      ewgts[index]     = (flt)ccs_contribute[j];
	      index++;
	      start[n+1] = index;
	    }
	}
    }

#ifdef DEBUG
  fprintf(stderr, "Verify start, nontrivial adjacency lists:\n");

  for ( i=0; i<nn; i++)
    {
      fprintf(stderr, "Node (%d) connects to: ", i+1);
      for ( j=start[i]; j<start[i+1]; j++)
	{
	  fprintf(stderr, "%d(%g) ", adjacency[j], ewgts[j]);
	}
      fprintf(stderr, "\n");
    }

#endif

  /*
   * This is the reference to the main routine for Chaco 2.0 usage.
   */

  if ( num_pieces > 1 )
    {
#ifdef DEBUG
      fprintf(stderr, "interface() called with:\n");
      fprintf(stderr, "nvtxs          = %d\n", nvtxs);
      fprintf(stderr, "start[0]       = %d\n", start[0]);
      fprintf(stderr, "architecture   = %d\n", architecture);
      fprintf(stderr, "ndims_tot      = %d\n", ndims_tot);
      fprintf(stderr, "mesh_dims[0]   = %d\n", mesh_dims[0]);
      fprintf(stderr, "ndims          = %d\n", ndims);
#endif
      err = interface(nvtxs, start, adjacency, vwgts, ewgts, x, y, z,
		      outassignname, outfilename, assignment, architecture,
		      ndims_tot, mesh_dims, goal, global_method, local_method,
		      rqi_flag, vmax, ndims, eigtol, seed);
    }
  else if ( num_pieces == 1 )
    {
      /*
       * The trivial partition...
       */
      for ( i=0; i<nvtxs; i++)
	{
	  assignment[i] = 1;
	}
    }
  else
    {
      sr = sprintf(err_msg, "? number of pieces = %d", num_pieces);
      EH(-1, err_msg);
    }

      
  if ( err != 0 )
    {
      EH(-1, "Problem return from Chaco interface().");
    }

#ifdef DEBUG
  fprintf(stderr, "Chaco interface() returns %d\n", err);
  fprintf(stderr, "assignments:\n");
  for ( i=0; i<nn; i++)
    {
      fprintf(stderr, "Node (%2d) in Set %2d\n", i+1, assignment[i]);
    }
#endif

#endif /* CHACO */



#ifdef POTATO_CHIP_DIET
  
  /* Comment:  Scott A Roberts, 2010-08-24
   * I really don't like doing this, but I'm going to muck around with
   * the nodal assignments for the special case where there are shells
   * with friends.  I want to ensure that the node that is "behind" 
   * a given shell node is owned by the same processor as the shell
   * node.  So if I find a disparity here, I will change the owner of the 
   * shell element to match that of the friend (as a shell node affects
   * many less other nodes than a friend does).  This is done to combat
   * the "potato chip" problem that we've seen.  Currently, it will only
   * work with SHELL4-HEX8 interactions, and only when there is only
   * one element friend.  If this screws up your problem, just undefine
   * "POTATO_CHIP_DIET" in the Makefile.
   */
  /* Loop through through all elements, investigating those shells with a friend */
  int fe, ei, ele, match, gnn;
  int shn, shnl;
  int num_friends;
  int *friend_list;
  friend_list = smalloc(30*sizeof(int));
  int num_unique_nodes, num_shared_nodes;
  int *unique_nodes, *shared_nodes;
  unique_nodes = smalloc(30*sizeof(int));
  shared_nodes = smalloc(30*sizeof(int));
  int *lcn;
  lcn = smalloc(3*sizeof(int));
  for ( e = 0; e < mono->num_elems; e++ ) {
    if (is_shell_element(mono,e)) {
      num_friends = find_element_friends(mono,e,friend_list);
      if (num_friends == 1) {
	fe = friend_list[0];
	
	/* Now that I have a friend element, figure out which nodes are not shared */
	num_unique_nodes = num_shared_nodes = 0;
	for ( n = mono->elem_node_pntr[fe]; n < mono->elem_node_pntr[fe+1]; n++ ) {
	  node = mono->elem_node_list[n];
	  match = 0;
	  for ( ei = mono->node_elem_pntr[node]; ei < mono->node_elem_pntr[node+1]; ei++) {
	    ele = mono->node_elem_list[ei];
	    if (ele == e) match = 1;
	  }
	  if (!match) BULL(node,unique_nodes,num_unique_nodes);
	  if ( match) BULL(node,shared_nodes,num_shared_nodes);
	}
	
	/* I now have the unique nodes, find which shared node is connected */
	for ( n = 0; n < num_shared_nodes; n++) {
	  shn = shared_nodes[n];
	  shnl = find_local_node_number( mono, fe, shn );
	  find_edge_connected_nodes(shnl,lcn);
	  for ( i = 0; i < 3; i++) {
	    gnn = mono->elem_node_list[mono->elem_node_pntr[fe]+lcn[i]];
	    if (!is_node_in_element(mono,gnn,e)) break;
	  }
	  
	  /* Let's set the shared assignment equal to the unique assignment */
	  if ( assignment[shn] != assignment[gnn] ) {
#ifdef DEBUG
	    printf("Moved assignment of node %d from proc %d to %d\n",shn,assignment[shn],assignment[gnn]);
#endif /* DEBUG */
	    assignment[shn] = assignment[gnn];
	  }
	  
	}
      } 
    }
  }


#endif /* POTATO_CHIP_DIET */



  /*
   * Just for convenience below, make an array that quickly tells
   * if a node is a boundary node or an internal node.
   *
   *       Internal nodes only connect to nodes in the same set.
   *
   *       Boundary nodes connect to nodes in more than one set.
   *
   *	   External nodes are in the "other" set and connect to more than 
   *       one set.
   *
   * (psm[node+1] -psm[node]) = 
   *    number of distinct sets to which me and my neighboring nodes belong
   *
   * psm[node] = first spot in sm where this nodes neighbor set 
   *             membership begins
   *
   * if ( psm[node+1] - psm[node] ) > 1, then this is more than an internal
   *                                     node, it's part of a boundary/external
   *                                     interaction.
   *
   * sm[i] - set membership list for all nodes in the global mesh
   */


  psm = (int *) smalloc((nn+1)*SZ_INT);

  psm[0] = 0;

  for ( i=0; i<nn; i++)
    {
      my_set = assignment[i];

      /*
       * For this node, clean out the list of distinct neighbor sets.
       */

      INIT_IVEC(setalogue, UNDEFINED_SET_NAME, MAX_ADJOINING_SETS)

      num_sets = 0;

      setalogue[num_sets] = my_set;

      num_sets++;

      /*
       * Look through all the neighbors. If they belong to a different
       * set, then add it to the list.
       */

#ifdef DEBUG
      fprintf(stderr, "At node (%d):", i+1);
#endif

      for ( j=pnn[i]; j<pnn[i+1]; j++ )
	{
#ifdef DEBUG
	  fprintf(stderr, " (%d)", var_node_names[j]+1);
#endif
	  neighbor_set_name = assignment[var_node_names[j]];

	  BULL(neighbor_set_name, setalogue, num_sets);

	}
#ifdef DEBUG
      fprintf(stderr, "\n");
#endif      
      psm[i+1] = num_sets;
    }


  /*
   * Convert occurrence count into a pointer list.
   */
  
  for ( i=0; i<nn; i++)
    {
      psm[i+1] += psm[i];
    }

  len_sm = psm[nn];

  sm = (int *) smalloc(len_sm*SZ_INT);

  INIT_IVEC(sm, UNDEFINED_SET_NAME, len_sm);

  for ( i=0; i<nn; i++)
    {
      my_set = assignment[i];

      /*
       * Clean out this node's list of distinct neighbor sets.
       */

      INIT_IVEC(setalogue, UNDEFINED_SET_NAME, MAX_ADJOINING_SETS)

      num_sets = 0;

      setalogue[num_sets] = my_set; /* It appears first! */

      num_sets++;

      /*
       * Look through all the neighbors. If they belong to a different
       * set, then add it to the list.
       */

      for ( j=pnn[i]; j<pnn[i+1]; j++ )
	{
	  neighbor_set_name = assignment[var_node_names[j]];

	  BULL(neighbor_set_name, setalogue, num_sets);
	}

      /*
       * Transcribe this node's setalogue into the global one.
       */

      for ( k=0; k<num_sets; k++)
	{
	  sm[psm[i]+k] = setalogue[k];
	}

    }

#ifdef DEBUG
  fprintf(stderr, "Set memberships...\n");

  for ( node=0; node<nn; node++)
    {
      fprintf(stderr, "Sets associated with node (%d): ", node+1);
      
      for ( j=psm[node]; j<psm[node+1]; j++)
	{
	  fprintf(stderr, "%d ", sm[j]);
	}
      fprintf(stderr, "\n");
    }
#endif

  /*
   * Element ownership determination.
   */

  assign_elem_ownership(mono, num_pieces, psm, sm, &element_owner, 
			&element_owner_dist, &element_bnd_stat);

#ifdef DEBUG
  fprintf(stderr, "Element ownership assigned.\n");
  for ( i=0; i<mono->num_elems; i++)
    {
      fprintf(stderr, "\tElement [%d] owned by set/proc %d (%s)\n", i+1, 
	      element_owner[i], (element_bnd_stat[i] ?"b":"i"));
    }

  fprintf(stderr, "Element owner distribution\n");
  for ( i=0; i<num_pieces; i++)
    {
      fprintf(stderr, "Set/proc [%d] owns %d elements.\n", i, 
	      element_owner_dist[i]);
    }

#endif

  /*
   * For the monolith, build the extra information about orientation and which
   * neighbor face is connected to this element.
   */

  build_elem_elem_xtra(mono);

  /*
   * Later, once dpi is under full construction, put together parts
   * to identify this set/proc's portion of the global info, plus
   * build the dpi->proc neighbor elem owner[] as well as, say, the
   * number and names of all procs to talk with about element data.
   */



  /*
   * Stalk through each set, counting the number of internal nodes
   * the number of boundary nodes and the external nodes.
   *
   * When we're done, make up little EXODUS II databases *and* extra
   * distributed information about:
   *
   *	(i)   who do I need to communicate with for send/recv?
   *		A: Make up a list of Sets to send to
   *                             and Sets to recv fr
   *
   *	(ii)  what nodes are sending/receiving from where?
   *
   *		A: Make up lists for ea sendset and recvset
   *		   listing the global node names I will send or recv.
   *
   *
   *	(iii) maps from local elem, node, dof to global elem, node, dof.
   *		A: For ea
   *
   *
   * Make little EXODUS II files with correct eb ids in the right places
   * correct ss's and ns's, as well as elem_maps, node_maps, and dof_maps
   * to the global problem.
   *
   * Can we do each set one by one? Yes, I think so...
   */

  E                        = (Exo_DB *) smalloc(sizeof(Exo_DB));
  
  D                        = (Dpi *) smalloc(sizeof(Dpi));

  send_proc_names = (int *) smalloc(MAX_SEND_PROCS*SZ_INT);
  recv_proc_names = (int *) smalloc(MAX_RECV_PROCS*SZ_INT);

  /*
   * There will be no more locally than appeared globally; probably much less.
   */

  proc_eb_id      = (int *) smalloc(neb*SZ_INT);

  proc_eb_ptr      = (int *) smalloc((neb+1)*SZ_INT);

  ebi              = (int *) smalloc(neb*SZ_INT);

  new_proc_eb_ptr  = (int *) smalloc((neb+1)*SZ_INT);

  for ( s=0; s<num_pieces; s++)
    {

      init_exo_struct(E);

      init_dpi_struct(D);

      init_dpi_version(D);

      /*
       * Just grab the BRK_VERSION string that was defined in config.h for the
       * brkfix distribution for now. When defining and writing Dpi we
       * call the shots.
       */

      strcpy(D->dpi_version_string, BRK_VERSION);

      D->num_global_node_descriptions = num_kinds_nodes;
  
      D->global_node_description      = (int **) smalloc(num_kinds_nodes*
							 SZPINT);

      D->global_node_description[0]   = (int *) smalloc(num_kinds_nodes*
							LEN_NODE_DESCRIPTION*
							SZ_INT);
      for ( i=1; i<num_kinds_nodes; i++)
	{
	  D->global_node_description[i] = ( D->global_node_description[i-1] +
					    LEN_NODE_DESCRIPTION );
	  
	}

      /*
       * Initialize...
       */

      for ( i=0; i<num_kinds_nodes; i++)
	{
	  for ( j=0; j<LEN_NODE_DESCRIPTION; j++)
	    {
	      D->global_node_description[i][j] = -7;
	    }
	}

      /*
       * Set to meaningful values...
       */
      
      for ( k=0; k<num_kinds_nodes; k++)
	{
	  D->global_node_description[k][0] = pnd[k]->num_basic_eqnvars;
	  for ( i=0; i<MAX_EQNVARS; i++)
	    {
	      D->global_node_description[k][1+i] = pnd[k]->eqnvar_ids[i];
	    }
	  for ( i=0; i<MAX_EQNVARS; i++)
	    {
	      D->global_node_description[k][1+MAX_EQNVARS+3*i]   = 
		pnd[k]->eqnvar_wts[i][0];
	      D->global_node_description[k][1+MAX_EQNVARS+3*i+1] = 
		pnd[k]->eqnvar_wts[i][1];
	      D->global_node_description[k][1+MAX_EQNVARS+3*i+2] = 
		pnd[k]->eqnvar_wts[i][2];
	    }
	}

      D->num_dofs_global          = node_dof0[nn];

      /*
       * Much of the monoliths information needs to be reproduced in the
       * child.
       */

      exo_dpi_clone(mono, D);


#ifdef DEBUG
      fprintf(stderr, "Set %d\n", s);
#endif
      /*
       * How many nodes from the monolith problem belong to Set s?
       *
       * There are:
       *
       *	internal nodes -- Set S has primary responsibility for
       *			  these nodes and equations, none of them
       *			  are needed by adjoining processors(sets)
       *
       *			  These nodes connect only to nodes that
       *                          are primarily in the same set.
       *
       *	boundary nodes -- Set S has primary responsibility for
       *			  these nodes and equations, each of which
       *			  is needed by one or more adjoining sets
       *
       *			  These nodes are primarily in this set
       *			  and secondarily in another set.
       *
       *
       *        external nodes -- Set S needs variables from these nodes
       *			  that are the primary responsibility of
       *			  other sets.
       *
       *			  These nodes are primarily in "other set"
       *			  and connect secondarily with "this set".
       *
       *
       * Count up each kind, then allocate space for lists of each kind.
       *
       * For the boundary nodes, make up a list of the other sets that will
       * need information from my boundary nodes. Each set will have a pointer
       * into a big list of boundary nodes. Note that some boundary nodes could
       * well be listed more than once, since more than one processor may need
       * that information.
       *
       * 
       */

      num_internal_nodes = 0;
      num_boundary_nodes = 0;
      num_external_nodes = 0;

      for ( n=0; n<nn; n++ )
	{
	  nsets = psm[n+1]-psm[n];

	  if ( nsets == 1 )
	    {
	      if ( assignment[n] == s )
		{
		  num_internal_nodes++;
		}
	    }
	  else if ( nsets > 1 )
	    {
	      if ( assignment[n] == s )
		{
		  num_boundary_nodes++;
		}
	      else
		{
		  /*
		   * Check - maybe this multi-set node whose primary assignment
		   * is a different set might have Set s as one of its 
		   * secondaries?
		   *
		   * No need to check the primary set again; we already 
		   * know that it is not s. Therefore the "+1" below.
		   */
		  
		  where = in_list(s, sm+psm[n], nsets);
		  
		  if (  where != -1 ) 
		    {
		      num_external_nodes++;
		    }
		}
	    }
	  else
	    {
	      EH(-1, "Every node needs to belong to at least one set.");
	    }
	}

#ifdef DEBUG
      fprintf(stderr, "Set %d:\n", s);
      fprintf(stderr, "        num_internal_nodes = %d\n", num_internal_nodes);
      fprintf(stderr, "        num_boundary_nodes = %d\n", num_boundary_nodes);
      fprintf(stderr, "        num_external_nodes = %d\n", num_external_nodes);
#endif      

      internal_nodes = (int *) smalloc(num_internal_nodes*SZ_INT);

      boundary_nodes = (int *) smalloc(num_boundary_nodes*SZ_INT);

      external_nodes = (int *) smalloc(num_external_nodes*SZ_INT);

      INIT_IVEC(internal_nodes, -1, num_internal_nodes);

      INIT_IVEC(boundary_nodes, -1, num_boundary_nodes);

      INIT_IVEC(external_nodes, -1, num_external_nodes);

      /*
       * The universe from this set/proc's perspective...
       */

      num_universe_nodes = ( num_internal_nodes + 
			     num_boundary_nodes + 
			     num_external_nodes );

      /*
       * Re-scan all nodes, this time collecting the 3 kinds for this set.
       */

      index_internal_node = 0;
      index_boundary_node = 0;
      index_external_node = 0;

      for ( n=0; n<nn; n++ )
	{
	  nsets = psm[n+1]-psm[n];

	  if ( nsets == 1 )
	    {
	      if ( assignment[n] == s )
		{
		  internal_nodes[index_internal_node++] = n;
		}
	    }
	  else if ( nsets > 1 )
	    {
	      if ( assignment[n] == s )
		{
		  boundary_nodes[index_boundary_node++] = n;
		}
	      else
		{
		  /*
		   * Check - maybe this multi-set node whose primary assignment
		   * is a different set might have Set s as one of its 
		   * secondaries?
		   *
		   * No need to check the primary set again; we already 
		   * know that it is not s. Therefore the "+1" below.
		   */
		  
		  where = in_list(s, sm+psm[n], nsets);
		  
		  if (  where != -1 ) 
		    {
		      external_nodes[index_external_node++] = n;
		    }
		}
	    }
	  else
	    {
	      EH(-1, "Every node needs to belong to at least one set.");
	    }
	}

      if ( index_internal_node != num_internal_nodes )
	{
	  fprintf(stderr, 
		  "Inconsistent count of internal nodes, set %d (%d != %d)\n",
		  s, index_internal_node, num_internal_nodes);
	  exit(0);
	}

      if ( index_boundary_node != num_boundary_nodes )
	{
	  fprintf(stderr, 
		  "Inconsistent count of boundary nodes, set %d (%d != %d)\n",
		  s, index_boundary_node, num_boundary_nodes);
	  exit(0);
	}

      if ( index_external_node != num_external_nodes )
	{
	  fprintf(stderr, 
		  "Inconsistent count of external nodes, set %d (%d != %d)\n",
		  s, index_external_node, num_external_nodes);
	  exit(0);
	}

#ifdef DEBUG
      printf("\nSET %d",s);
      printf("\n   Internal nodes:");
      for (i=0;i<num_internal_nodes;i++) printf("%4d",internal_nodes[i]+1);
      printf("\n   Boundary nodes:");
      for (i=0;i<num_boundary_nodes;i++) printf("%4d",boundary_nodes[i]+1);
      printf("\n   External nodes:");
      for (i=0;i<num_external_nodes;i++) printf("%4d",external_nodes[i]+1);
      printf("\n");
#endif

      /*
       * New - some up the dof weight for internal, boundary, and external
       * nodes for this set/proc.
       */

      internal_dofweight = 0;
      boundary_dofweight = 0;
      external_dofweight = 0;

      for ( i=0; i<num_internal_nodes; i++)
	{
	  global_node_index   = internal_nodes[i];
	  internal_dofweight += ( node_dof0[global_node_index+1] -
				  node_dof0[global_node_index] ); 
	}

      for ( i=0; i<num_boundary_nodes; i++)
	{
	  global_node_index   = boundary_nodes[i];
	  boundary_dofweight += ( node_dof0[global_node_index+1] -
				  node_dof0[global_node_index] ); 
	}

      for ( i=0; i<num_external_nodes; i++)
	{
	  global_node_index   = external_nodes[i];
	  external_dofweight += ( node_dof0[global_node_index+1] -
				  node_dof0[global_node_index] ); 
	}

      /*
       * Add to running sum accumulated over all processors.
       */

      total_internal_dofweight += internal_dofweight;
      total_boundary_dofweight += boundary_dofweight;
      total_external_dofweight += external_dofweight;

      /*
       * Sort external nodes so that, if more than one neighbor processor
       * owns external nodes, the nodes for each processor fall into a
       * contiguous numbering as far as this processor is concerned.
       *
       * Also, the clumps of external nodes are ordered so that lowest
       * processor ID's occur first.
       */

#ifdef DEBUG
      fprintf(stderr, "Before sorting external nodes for proc %d...\n", s);
      /*
      for ( i=0; i<num_external_nodes; i++)
	{
	  fprintf(stderr, "external_nodes[%d] = %d (owned by %d)\n",
		  i, external_nodes[i], assignment[external_nodes[i]]);
	}

      */
      fprintf(stderr, "num_external_nodes = %d\n", num_external_nodes);
#endif
      proc_sort(external_nodes, num_external_nodes, nvtxs, 
		assignment);
#ifdef DEBUG
      fprintf(stderr, "After sorting external nodes for proc %d...\n", s);
      /*
      for ( i=0; i<num_external_nodes; i++)
	{
	  fprintf(stderr, "external_nodes[%d] = %d (owned by %d)\n",
		  i, external_nodes[i], assignment[external_nodes[i]]);
	}
      */
#endif

#ifdef DEBUG
      /*
       * Verify that within each set of nodes (internal, boundary, external)
       * that monotonicity of the global node names is obeyed.
       */

      for ( i=1; i<num_internal_nodes; i++)
	{
	  if ( internal_nodes[i] < internal_nodes[i-1] )
	    {
	      sr = sprintf(err_msg, "internal_nodes[] nonmonotone in Set %d",
			   s);
	      EH(-1, err_msg);
	    }
	}

      for ( i=1; i<num_boundary_nodes; i++)
	{
	  if ( boundary_nodes[i] < boundary_nodes[i-1] )
	    {
	      sr = sprintf(err_msg, "boundary_nodes[] nonmonotone in Set %d",
			   s);
	      EH(-1, err_msg);
	    }
	}

      /*
       * The check below is now too rigorous. For ease of use with
       * Aztec, we'll reorder the external nodes so that each neighbor
       * processor's nodes appear contiguously.
       */

#if 0
      for ( i=1; i<num_external_nodes; i++)
	{
	  if ( external_nodes[i] < external_nodes[i-1] )
	    {
	      sr = sprintf(err_msg, "external_nodes[] nonmonotone in Set %d",
			   s);
	      EH(-1, err_msg);
	    }
	}
#endif /* 0 */

#endif

      /*
       * For convenience, make one big array of all nodes known to this
       * set/proc.
       */
      
      proc_nodes = (int *) smalloc(num_universe_nodes*SZ_INT);

      for ( i=0; i<num_internal_nodes; i++)
	{
	  proc_nodes[i] = internal_nodes[i];
	}

      for ( i=0; i<num_boundary_nodes; i++)
	{
	  proc_nodes[num_internal_nodes+i] = boundary_nodes[i];
	}

      for ( i=0; i<num_external_nodes; i++)
	{
	  proc_nodes[num_internal_nodes+num_boundary_nodes+i] = 
	    external_nodes[i];
	}

#ifdef DEBUG
      for ( i=0; i<num_universe_nodes; i++)
	{
	  fprintf(stderr, "proc_nodes[%d] = %d (%d)\n", i, 
		  proc_nodes[i], proc_nodes[i]+1);
	}
#endif

      /*
       * Figure out communication needs.
       *
       * Boundary nodes will need to be sent to other processors, in some
       * cases to more than one other processor.
       *
       * External nodes will need to be received from other processors, but
       * will only be originating from a distinct processor.
       *
       * First, count & catalog the names of the processors/sets that this set
       * will be sending to:
       */
      
      num_send_procs = 0;
      
      INIT_IVEC(send_proc_names, -1, MAX_SEND_PROCS);

      num_recv_procs = 0;
      
      INIT_IVEC(recv_proc_names, -1, MAX_RECV_PROCS);

      for ( i=0; i<num_boundary_nodes; i++)
	{
	  node = boundary_nodes[i];
	  for ( t=psm[node]+1; t<psm[node+1]; t++)
	    {
	      set_name = sm[t];

	      BULL(set_name, send_proc_names, num_send_procs);
	    }
	}

      /*
       * Since the external nodes and dofs are not only clumped by
       * processor, but ordered so lowest numbered processor appears
       * first, it is important to sequence the processor names properly
       * for GOMA -- lowest processors appear first.
       */

      qsort(send_proc_names, num_send_procs, sizeof(int), integer_compare);

#ifdef DEBUG
      fprintf(stderr, "Proc/Set %d has %d b nodes to send to:",
	      s, num_boundary_nodes);
      for ( i=0; i<num_send_procs; i++)
	{
	  fprintf(stderr, " %d", send_proc_names[i]);
	}
      fprintf(stderr, "\n");
#endif

      for ( i=0; i<num_external_nodes; i++)
	{
	  node = external_nodes[i];
	  
	  set_name = assignment[node];

	  BULL(set_name, recv_proc_names, num_recv_procs);
	}

#ifdef DEBUG
      fprintf(stderr, "Proc/Set %d has %d e nodes living on:",
	      s, num_external_nodes);
      for ( i=0; i<num_recv_procs; i++)
	{
	  fprintf(stderr, " %d", recv_proc_names[i]);
	}
      fprintf(stderr, "\n");
#endif

      /*
       * Gather together a list of the global element numbers that
       * will be traversed by this set.
       *
       * Any element that contains either an internal or boundary node
       * will be traversed.
       *
       * Count, allocate space for, and gather the list of elements for 
       * this set.
       */

      proc_ne = 0;

      for ( ge=0; ge<ne; ge++)
	{
	  /*
	   * Without contrary evidence, assume this element does not touch
	   * any of this set's nodes. Indeed, most elements do not.
	   */

	  touch = FALSE;

	  for ( n=ep[ge]; n<ep[ge+1]; n++ )
	    {
	      node   = nl[n];
	      touch |= IOL(node, internal_nodes, num_internal_nodes);
	      touch |= IOL(node, boundary_nodes, num_boundary_nodes);
	    }

	  if ( touch )
	    {
	      proc_ne++;
	    }
	}
      
#ifdef DEBUG      
      fprintf(stderr, "Set %d will traverse %d elems\n", s, proc_ne);
#endif

      if ( proc_ne < 1 || proc_ne > ne )
	{
	  sr = sprintf(err_msg, "Bad element count of %d in Set %d",
		       proc_ne, s);
	  EH(-1, err_msg);
	}

      proc_elems     = (int *) smalloc(proc_ne*SZ_INT);

      new_proc_elems = (int *) smalloc(proc_ne*SZ_INT);

      /*
       * Now collect those elements that will be traversed.
       */

      index_proc_elems = 0;
      for ( ge=0; ge<ne; ge++)
	{
	  touch = FALSE;

	  for ( n=ep[ge]; n<ep[ge+1]; n++ )
	    {
	      node   = nl[n];
	      touch |= IOL(node, internal_nodes, num_internal_nodes);
	      touch |= IOL(node, boundary_nodes, num_boundary_nodes);
	    }

	  if ( touch )
	    {
	      proc_elems[index_proc_elems] = ge;
	      index_proc_elems++;
	    }
	}
#ifdef DEBUG      
      fprintf(stderr, "Set %d has %d elems:\n", s, proc_ne);
      fprintf(stderr, "The global names, before ordering, are:\n");
      for ( i=0; i<proc_ne; i++)
	{
	  fprintf(stderr, "\tproc_elems[%d]=(%d)\n", i, proc_elems[i]+1);
	}
      fprintf(stderr, "\n");
#endif

      /*
       * Count up the number of distinct element blocks that this
       * set/proc will visit.
       */

      proc_neb = 0;

      /*
       * Initially, no valid eb ids have been found in this setprocs
       * element collection. They start nowhere
       */

      INIT_IVEC(proc_eb_id,    -1, neb);

      INIT_IVEC(proc_eb_ptr,     0, neb+1);

      INIT_IVEC(new_proc_eb_ptr, 0, neb+1);

      for ( e=0; e<proc_ne; e++)
	{
	  elem     = proc_elems[e];

	  eb_index = fence_post(elem, ebl, neb+1);

	  eb_id    = mono->eb_id[eb_index];

#ifdef DEBUG	  
	  fprintf(stderr, "Element (%d) belongs to EB[%d]=%d\n",
		  elem+1, eb_index, eb_id);
#endif

	  where    = in_list(eb_id, proc_eb_id, proc_neb);

	  if ( where == -1 )
	    {
	      /*
	       * A new eb id has been encountered for this proc/set.
	       */
	      proc_eb_id[proc_neb]    = eb_id;
	      where                    = proc_neb;
	      proc_neb++;
	    }
	  proc_eb_ptr[where]++;
	}

      /*
       * The beginning of each element blocks elements may be found by
       * summing up how many are in each part...
       */


#ifdef DEBUG
      /*
       * Consistency checking.
       */
      count = 0;
      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  count += proc_eb_ptr[ieb];
	  fprintf(stderr, "In Setproc %d, proc_eb_ptr[%d] for (ebid=%d) = %d\n",
		  s, ieb, proc_eb_id[ieb], proc_eb_ptr[ieb]);
	}

      if ( count != proc_ne )
	{
	  EH(-1, "Mismatch!");
	}
#endif      

      /*
       * Transform frequency of occurrence into pointers into proc_elem[].
       *
       * Shift each entry up by one to make space, then accumulate a 
       * running sum.
       */
      
      for ( i=proc_neb; i>0; i--)
	{
	  proc_eb_ptr[i] = proc_eb_ptr[i-1];
	}
      proc_eb_ptr[0] = 0;
      for ( i=0; i<proc_neb; i++)
	{
	  proc_eb_ptr[i+1] += proc_eb_ptr[i];
	}

      /*
       * Node numbering: internal, boundary, external.
       *
       * If there are no other constraints, then why not
       * order nodes so that the global node number increases
       * monotonically with local node number? This might
       * speed up searches as well as coordinating memory access
       * somewhat(?). Actually, given the way that the arrays
       * internal_nodes[], boundary_nodes[], and external_nodes[] 
       * were constructed, they should automatically be monotonically
       * ordered.
       *
       * Note, now that external_nodes are locally monotone on a per-processor
       * basis!
       */

      /*
       * Element numbering: renumber the elements and element block indeces
       * so that the element block with the greatest internal computational 
       * load appears first. Within that element block, the elements must,
       * of course, be numbered consecutively. However, do attempt to 
       * have low number elements require only internal information, then
       * proceeding to higher numbered elements that may contain external
       * nodes. 
       */

      /*
       * The true computational load count is complicated, requiring a
       * check against the active eqnvars at this node in each element
       * block. Use a quicker heuristic for the time being:
       *
       *     The element block with the greatest number of elements that
       *     strictly contain only internal and boundary nodes will become
       *     the element block with the lowest element block index.
       *
       *     Furthermore, the local element numbering scheme will be
       *     such that external nodes are avoided as long as possible.
       *
       * While we're counting, count up the total number of elements that
       * belongs to each element block this set/proc contains.
       */


      private_elem_count = (int *) smalloc(neb*SZ_INT);

      INIT_IVEC(private_elem_count, 0, neb);

      for ( e=0; e<proc_ne; e++)
	{
	  elem = proc_elems[e];	/* global element name */
		  
	  private_elem = TRUE;

	  for ( l=ep[elem]; l<ep[elem+1]; l++ )
	    {
	      node = nl[l];

	      /*
	       * If this global node appears in the list of external nodes,
	       * then clearly it's *NOT* a private element.
	       */ 

#ifdef DEBUG
	      if ( s == 0 )
		{
		  fprintf(stderr, "Looking for node %d in list of length %d\n",
			  node, num_external_nodes);
		  fprintf(stderr, "external_nodes[] = ");
		  for ( k=0; k<num_external_nodes; k++)
		    {
		      fprintf(stderr, "%d ", external_nodes[k]);
		    }
		  fprintf(stderr, "\n");
		}
#endif

	      if ( IUL(node, external_nodes, num_external_nodes) )
		{
		  private_elem = FALSE;
		}
	    }

	  glob_eb_index = fence_post(elem, ebl, neb+1);

	  eb_id         = mono->eb_id[glob_eb_index];

	  proc_eb_index = in_list(eb_id, proc_eb_id, proc_neb);

	  /* proc_eb_elem_count[proc_eb_index]++; */

#ifdef DEBUG
	  fprintf(stderr, "Set %d elem[%d](=(%d)) is %s in local EB[%d]=%d\n",
		  s, e, elem+1, (private_elem?("priv"):("shar")),
		  proc_eb_index, eb_id);
#endif

	  if ( private_elem )
	    {
	      private_elem_count[proc_eb_index]++;
	    }
	}

      /*
       * Within each element block, re-order the elements so that private
       * elements appear first, then shared elements.
       */

      /*
       * Build lists of private and shared elements for elements within
       * each element block index for this set/proc.
       */

      proc_elem_priv = (int *) smalloc(proc_ne*SZ_INT);
      proc_elem_shar = (int *) smalloc(proc_ne*SZ_INT);

      begin = 0;

      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  eb_id = proc_eb_id[ieb];
	  
	  INIT_IVEC(proc_elem_priv, -1, proc_ne);
	  
	  INIT_IVEC(proc_elem_shar, -1, proc_ne);
	  
	  index_priv = 0;
	  
	  index_shar = 0;

#ifdef DEBUG
	  fprintf(stderr, "proc_eb_id[%d] = %d\n", ieb, proc_eb_id[ieb]);
	  fprintf(stderr, "proc_eb_ptr[%d] = %d\n", ieb, proc_eb_ptr[ieb]);
#endif

	  for ( e=proc_eb_ptr[ieb]; e<proc_eb_ptr[ieb+1]; e++)
	    {
	      elem = proc_elems[e];

	      private_elem = TRUE;

	      for ( l=ep[elem]; l<ep[elem+1]; l++ )
		{
		  node = nl[l];

		  /*
		   * Private elements have no external nodes.
		   */

		  private_elem &= !(IUL(node, external_nodes, 
				       num_external_nodes));
		}

	      if ( private_elem )
		{
		  proc_elem_priv[index_priv] = elem;
		  index_priv++;
		}
	      else
		{
		  proc_elem_shar[index_shar] = elem;
		  index_shar++;
		}
#ifdef DEBUG
	      fprintf(stderr, "proc_elems[%d] = (%d) is %s.\n", e,
		      elem+1, (private_elem?"private":"shared"));
#endif
	    }

	  if ( index_priv != private_elem_count[ieb] )
	    {
	      sr = sprintf(err_msg, 
			   "index_priv=%d != private_elem_count[%d]=%d EBID %d",
			   index_priv, ieb, private_elem_count[ieb], eb_id);
	      EH(-1, err_msg);
	    }

	  if ( index_priv+index_shar != proc_eb_ptr[ieb+1]-proc_eb_ptr[ieb] )
	    {
	      sr = sprintf(err_msg, 
			   "#priv elems(%d)+#shar elems(%d) != total(%d)",
			   index_priv, index_shar, 
			   proc_eb_ptr[ieb+1]-proc_eb_ptr[ieb]);
	      EH(-1, err_msg);
	    }

	  /*
	   * Resequence the global element names for this element block
	   * in this set/proc.
	   */

	  begin = proc_eb_ptr[ieb];

	  for ( e=0; e<index_priv; e++)
	    {
	      proc_elems[begin+e] = proc_elem_priv[e];
	    }

	  begin = proc_eb_ptr[ieb] + index_priv;

	  for ( e=0; e<index_shar; e++)
	    {
	      proc_elems[begin+e] = proc_elem_shar[e];
	    }
	}

      /*
       * Re-order the element blocks for this set/proc so that the
       * element block with the greatest number of private elems appears
       * first, followed by the element blocks with successively fewer
       * quantities of private elements.
       */

#ifdef DEBUG
      fprintf(stderr, "BEFORE: Element block ordering for Set %d\n", s);
      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  fprintf(stderr, "EB[%d] = %d has %d/%d (private/total) elems.\n",
		  ieb, 
		  proc_eb_id[ieb], 
		  private_elem_count[ieb],
		  proc_eb_ptr[ieb+1]-proc_eb_ptr[ieb]);
	}
      fprintf(stderr, "\n");
#endif

      /*
       * Save the baby whale pointers!
       */

      for ( i=0; i<proc_neb; i++)
	{
	  ebi[i] = i;
	}

      
      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  index_max_private_elems = ieb;
	  max_private_elems = private_elem_count[ieb];

	  /*
	   * Search subsequent element blocks for more privacy.
	   */

	  for ( jeb=ieb; jeb<proc_neb; jeb++)
	    {
	      if ( private_elem_count[jeb] > max_private_elems )
		{
		  max_private_elems       = private_elem_count[jeb];
		  index_max_private_elems = jeb;
		}
	    }
	  
	  /*
	   * If the subsequent list had the max, then exchange each
	   * array element that depends on the local element block index.
	   */

	  if ( index_max_private_elems != ieb )
	    {
	      ISWAP(private_elem_count[ieb], 
		    private_elem_count[index_max_private_elems]);

	      ISWAP(proc_eb_id[ieb], 
		    proc_eb_id[index_max_private_elems]);

	      ISWAP(ebi[ieb],
		    ebi[index_max_private_elems]);
	    }
	}


#ifdef DEBUG      
      for ( i=0; i<=proc_neb; i++)
	{
	  fprintf(stderr, "\tproc_eb_ptr[%d]=%d\n", i, proc_eb_ptr[i]);
	}
#endif

      new_start = 0;

      for (ieb=0; ieb<proc_neb; ieb++)
	{
	  old_start            = proc_eb_ptr[ebi[ieb]];
	  
	  len                  = proc_eb_ptr[ebi[ieb]+1] - old_start;

	  new_proc_eb_ptr[ieb] = new_start;

	  for ( i=0; i<len ; i++)
	    {
	      new_proc_elems[new_start+i] = proc_elems[old_start+i];

#ifdef DEBUG
	      fprintf(stderr, "new_proc_elems[%d+%d] = %d (proc_elem[%d+%d])\n",
		      new_start,i, new_proc_elems[new_start+i],
		      old_start,i);
#endif	      
	    }

	  new_start += len;
	}

      new_proc_eb_ptr[proc_neb] = new_start;


      /*
       * Done. Now, overwrite the originals for convenience down below.
       */

      for ( i=0; i<proc_ne; i++)
	{
	  proc_elems[i] = new_proc_elems[i];
	}

      free(new_proc_elems);
      new_proc_elems = NULL;

      for ( i=0; i<=proc_neb; i++)
	{
	  proc_eb_ptr[i] = new_proc_eb_ptr[i];
	}

#ifdef DEBUG
      fprintf(stderr, "AFTER: Element block ordering for Set %d\n", s);
      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  fprintf(stderr, "EB[%d] = %d has %d/%d (private/total) elems.\n",
		  ieb, 
		  proc_eb_id[ieb], 
		  private_elem_count[ieb],
		  (proc_eb_ptr[ieb+1]-proc_eb_ptr[ieb]));
	}
#endif

      /*
       * This proc's element list should already be contiguous from the
       * global perspective. That is, the same element block ID will not
       * appear sandwiched by other element block ids.
       *
       * Now do a further refinement of the proc_elem[] array to order
       * the internal/private elements first, followed by those elements
       * that are not internal/private.
       */

#ifdef DEBUG
      for ( i=0; i<proc_neb; i++)
	{
	  fprintf(stderr, "proc_eb_start[%d]      = %d\n", 
		  i, proc_eb_ptr[i]);
	  fprintf(stderr, "private_elem_count[%d] = %d\n", 
		  i, private_elem_count[i]);
	  fprintf(stderr, "prob_eb_length[%d]     = %d\n",
		  i, proc_eb_ptr[i+1]-proc_eb_ptr[i]);
	}
      
      fprintf(stderr, "proc_elems after ordering...\n");
      
      for ( e=0; e<proc_ne; e++)
	{
	  fprintf(stderr, "\tproc_elems[%d] -> (%d)\n", e, proc_elems[e] + 1);
	}



      fprintf(stderr, "Starting nodesets...\n");
#endif

      /*
       * Node Sets:
       *
       * Look at every node in every global node set and check to see
       * if it is one of this set/proc's nodes.
       *
       * If it is, then increment the number of required to store this
       * proc's nodesets and add the global ns id to the unique list.
       * Count the number of global node sets and global nodes in nodesets
       * that appear in this set/proc's nodes...
       *
       * To avoid two passes, make a provisional list of those nodes and
       * of those nodes set ids.
       */


      proc_ns_id             = (int *) smalloc(INIT_PROC_NUM_NS*SZ_INT);

      proc_ns_num_nodes      = (int *) smalloc(INIT_PROC_NUM_NS*SZ_INT);

      proc_ns_num_distfacts  = (int *) smalloc(INIT_PROC_NUM_NS*SZ_INT);

      proc_ns_node_index     = (int *) smalloc(INIT_PROC_NUM_NS*SZ_INT);

      proc_ns_distfact_index = (int *) smalloc(INIT_PROC_NUM_NS*SZ_INT);

      size_proc_num_ns       = INIT_PROC_NUM_NS;

      INIT_IVEC(proc_ns_id,             -1, INIT_PROC_NUM_NS);

      INIT_IVEC(proc_ns_num_nodes,       0, INIT_PROC_NUM_NS);

      INIT_IVEC(proc_ns_num_distfacts,   0, INIT_PROC_NUM_NS);

      INIT_IVEC(proc_ns_node_index,     -1, INIT_PROC_NUM_NS);

      INIT_IVEC(proc_ns_distfact_index, -1, INIT_PROC_NUM_NS);


      proc_ns_node_list      = (int *) smalloc(INIT_PROC_NS_NODE_LIST_LENGTH*					 SZ_INT);

      proc_ns_node_list_index_global
	= (int *) smalloc(INIT_PROC_NS_NODE_LIST_LENGTH*
			  SZ_INT);

      size_proc_ns_node_list = INIT_PROC_NS_NODE_LIST_LENGTH;

      INIT_IVEC(proc_ns_node_list, -1, INIT_PROC_NS_NODE_LIST_LENGTH);


      proc_ns_distfact_list  = (dbl *) 
	smalloc(INIT_PROC_NS_DISTFACT_LIST_LENGTH*SZ_DBL);

      proc_ns_distfact_list_index_global = (int *)
	smalloc(INIT_PROC_NS_DISTFACT_LIST_LENGTH*SZ_INT);

      size_proc_ns_distfact_list = INIT_PROC_NS_DISTFACT_LIST_LENGTH;

      /*
       * Counters...
       */

      proc_num_node_sets   = 0;

      proc_ns_node_len     = 0;

      proc_ns_distfact_len = 0;

      for ( index_ns=0; index_ns<mono->num_node_sets; index_ns++)
	{
	  /*
	   * If *any* node in this nodeset belongs to this set/proc's nodes,
	   * then log the name of the nodeset and appropriate counters.
	   *
	   * Look over only the internal and boundary nodes. Don't worry
	   * about incomplete ideas of nodesets for nodes that some other
	   * processor is responsible for.
	   *
	   * Every node in this nodeset that belongs to this set/proc's nodes
	   * will be recorded in appropriate arrays and counters incremented.
	   *
	   * Now, too, the index location in the global concatenated list
	   * will be recorded, too, to make it easier to put it back later
	   * down the road.
	   */

	  begin    = mono->ns_node_index[index_ns];

	  begin_df = mono->ns_distfact_index[index_ns];

	  ns_touch = FALSE;

	  start_nd    = proc_ns_node_len;

	  start_df = proc_ns_distfact_len;

	  /*
	   * Number of nodes in this global nodeset.
	   */

	  gnnns    = mono->ns_num_nodes[index_ns];
	  if ( gnnns < 1 ) 
	    {
	      sr = sprintf(err_msg, "nodeset %d has %d nodes - punting!",
			   mono->ns_id[index_ns], 
			   mono->ns_num_nodes[index_ns]);
	      EH(-1, err_msg);
	    }

	  /*
	   * Number of distribution factors in this global nodeset PER NODE.
	   * Should be 0, 1, 2, ...
	   */

	  gndfpn    = mono->ns_num_distfacts[index_ns]/gnnns;

	  for ( j=0; j<gnnns; j++)
	    {
	      node = mono->ns_node_list[begin+j];
#ifdef DEBUG
	      fprintf(stderr, "NS[%d]=%d, Looking for gn[%d]=%d\n",
		      index_ns, mono->ns_id[index_ns], j, node+1);
#endif
	      n    = get_internal_boundary_index(node, proc_nodes,
						 num_internal_nodes,
						 num_boundary_nodes);
	      if ( n != -1 )
		{
		  proc_ns_node_list[proc_ns_node_len] = n;

		  proc_ns_node_list_index_global[proc_ns_node_len] = begin+j;

		  proc_ns_node_len++;



		  /*
		   * Segmentation fault next time unless we reallocate more
		   * space for the node list.
		   */

		  if ( proc_ns_node_len > size_proc_ns_node_list - 1 )
		    {
		      size_proc_ns_node_list += INIT_PROC_NS_NODE_LIST_LENGTH;
		      proc_ns_node_list = (int *) 
			realloc(proc_ns_node_list, 
				SZ_INT * size_proc_ns_node_list);
		      proc_ns_node_list_index_global = (int *) 
			realloc(proc_ns_node_list_index_global, 
				SZ_INT * size_proc_ns_node_list);
		    }

		  ns_touch = TRUE;

		  /*
		   * Are there any distribution factors to do here?
		   *
		   * All nodes in this nodeset assumed to have gndfpn
		   * number of distribution factors.
		   */

		  if ( gndfpn > 0 )
		    {
		      where = begin_df + j*gndfpn;
		      for ( m=0; m<gndfpn; m++)
			{
			  proc_ns_distfact_list[proc_ns_distfact_len] = 
			    mono->ns_distfact_list[where+m];

			  /*
			   * Mom: "Where did you get this, Johnny?"
			   *
			   * Johnny: "Aw shucks, ma, can't I keep `em?"
			   */

			  proc_ns_distfact_list_index_global
			    [proc_ns_distfact_len] = where + m;

			  proc_ns_distfact_len++;

			  if ( proc_ns_distfact_len > 
			       size_proc_ns_distfact_list - 1 )
			    {
			      size_proc_ns_distfact_list += 
				INIT_PROC_NS_DISTFACT_LIST_LENGTH;
			      proc_ns_distfact_list = (dbl *)
				realloc(proc_ns_distfact_list,
					SZ_DBL* size_proc_ns_distfact_list);
			      proc_ns_distfact_list_index_global = (int *)
				realloc(proc_ns_distfact_list_index_global,
					SZ_INT* size_proc_ns_distfact_list);
			    }
			}
		    }
		}
	    }

	  if ( ns_touch )
	    {
	      /*
	       * This "I" is just for convenience and brevity.
	       */

	      I                         = proc_num_node_sets;

	      proc_ns_id[I]             = mono->ns_id[index_ns];

	      proc_ns_num_nodes[I]      = proc_ns_node_len - start_nd;

	      proc_ns_num_distfacts[I]  = proc_ns_distfact_len - start_df;

	      proc_ns_node_index[I]     = start_nd;

	      proc_ns_distfact_index[I] = start_df;

	      proc_num_node_sets++;

	      if ( proc_num_node_sets > size_proc_num_ns - 1 )
		{
		  size_proc_num_ns      += INIT_PROC_NUM_NS;

		  proc_ns_id             = ( (int *) 
					     realloc(proc_ns_id, SZ_INT *
						     size_proc_num_ns));

		  proc_ns_num_nodes      = ( (int *) 
					     realloc(proc_ns_num_nodes, 
						     SZ_INT * 
						     size_proc_num_ns));

		  proc_ns_num_distfacts  = ( (int *) 
					     realloc(proc_ns_num_distfacts,
						     size_proc_num_ns *SZ_INT));

		  proc_ns_node_index     = ( (int *) 
					     realloc(proc_ns_node_index, SZ_INT*
						     size_proc_num_ns));

		  proc_ns_distfact_index = ( (int *) 
					     realloc(proc_ns_distfact_index,
						     SZ_INT *
						     size_proc_num_ns));
		}

	    }
	}

#ifdef DEBUG
      fprintf(stderr, "Assigning nodesets info to E...\n");
#endif

      E->num_node_sets     = proc_num_node_sets;

      E->ns_node_len       = proc_ns_node_len;

      E->ns_distfact_len   = proc_ns_distfact_len;

      if ( E->num_node_sets > 0 )
	{
	  E->ns_id             = (int *) smalloc(E->num_node_sets*SZ_INT);

	  E->ns_num_nodes      = (int *) smalloc(E->num_node_sets*SZ_INT);

	  E->ns_num_distfacts  = (int *) smalloc(E->num_node_sets*SZ_INT);

	  E->ns_node_index     = (int *) smalloc(E->num_node_sets*SZ_INT);

	  E->ns_distfact_index = (int *) smalloc(E->num_node_sets*SZ_INT);
	}

      if ( E->ns_node_len > 0 )
	{
	  E->ns_node_list      = (int *) smalloc(E->ns_node_len*SZ_INT);
	}

      if ( E->ns_distfact_len > 0 )
	{
	  E->ns_distfact_list  = (dbl *) smalloc(E->ns_distfact_len*SZ_DBL);
	}

      for ( i=0; i<E->num_node_sets; i++)
	{
	  E->ns_id[i]             = proc_ns_id[i];

	  E->ns_num_nodes[i]      = proc_ns_num_nodes[i];

	  E->ns_num_distfacts[i]  = proc_ns_num_distfacts[i];

	  E->ns_node_index[i]     = proc_ns_node_index[i];

	  E->ns_distfact_index[i] = proc_ns_distfact_index[i];
	}

      for ( i=0; i<E->ns_node_len; i++)
	{
	  E->ns_node_list[i] = proc_ns_node_list[i];
	}

      for ( i=0; i<E->ns_distfact_len; i++)
	{
	  E->ns_distfact_list[i] = proc_ns_distfact_list[i];
	}

#ifdef DEBUG
      fprintf(stderr, "Computing sidesets info to E...\n");
#endif

      /*
       * Side Sets:
       *
       * Look at every elem in every global side set and check to see
       * if it is one of this set/proc's elems (private or shared).
       *
       * If it is, and this element/side part of the sideset includes any
       * internal or boundary nodes, then add this side to the list of
       * sides for this set/proc. 
       *
       * Also, add the global ss id to the this proc's list of ssids.
       */

      proc_ss_id             = (int *) smalloc(INIT_PROC_NUM_SS*SZ_INT);
      proc_ss_num_sides      = (int *) smalloc(INIT_PROC_NUM_SS*SZ_INT);
      proc_ss_num_distfacts  = (int *) smalloc(INIT_PROC_NUM_SS*SZ_INT);
      proc_ss_elem_index     = (int *) smalloc(INIT_PROC_NUM_SS*SZ_INT);
      proc_ss_distfact_index = (int *) smalloc(INIT_PROC_NUM_SS*SZ_INT);

      size_proc_num_ss       = INIT_PROC_NUM_SS;


      proc_ss_elem_list      = (int *) smalloc(INIT_PROC_SS_SIDE_LIST_LENGTH*
					       SZ_INT);
      proc_ss_side_list      = (int *) smalloc(INIT_PROC_SS_SIDE_LIST_LENGTH*
					       SZ_INT);
      proc_ss_elem_list_index_global 
	                     = (int *) smalloc(INIT_PROC_SS_SIDE_LIST_LENGTH*
					       SZ_INT);


      size_proc_ss_elem_list = INIT_PROC_SS_SIDE_LIST_LENGTH;

      proc_ss_distfact_list      = (dbl *) 
	smalloc(INIT_PROC_SS_DISTFACT_LIST_LENGTH*SZ_DBL);

      proc_ss_distfact_list_index_global = (int *) 
	smalloc(INIT_PROC_SS_DISTFACT_LIST_LENGTH*SZ_INT);

      size_proc_ss_distfact_list = INIT_PROC_SS_DISTFACT_LIST_LENGTH;

      /*
       * Distribution factors: for each side set there are either
       *			(i)  0 distribution factors
       *			(ii) n distribution factors
       * where (n) is the number of nodes in the node list. Note that
       * nodes will often be replicated for elements along a side that
       * share nodes. You'll need to convert sides to nodes and count
       * all nodes up to the current element in the global side set
       * in order to find out which distribution factors apply to this
       * element.
       */



      /*
       * Loop through every elem/side in every global sideset. If that
       * element appears anywhere in this set/proc's list of elements, 
       * private or shared, then that sideset will become one that this
       * set/proc could see. The determination will finally be based on
       * whether the particular sideset elem/side combination includes
       * any internal or boundary nodes. If so, it will be included. If
       * the elem/side touches only external nodes, then it will *not*
       * be included in this set/proc's list.
       *
       * News flash!!! All elem/sides that have all nodes in this processors
       *               node list, including external nodes will be included.
       * 
       *
       * And, this particular element and side will
       * be added to this setproc's list for this sideset id. Also,
       * any associated distribution factors are added for this element.
       * (We need to find out where in the global scheme this elements
       *  distribution factors are located - just use the ss_cnt_nodes[].)
       */

      proc_num_side_sets   = 0;

      proc_ss_elem_len     = 0;	/* and the side length, too. */

      proc_ss_distfact_len = 0;

#ifdef DEBUG
      fprintf(stderr, "\n\nSet/proc = %d checking sidesets.......\n\n", s);
#endif      

      for ( index_ss=0; index_ss<mono->num_side_sets; index_ss++)
	{
	  begin              = mono->ss_elem_index[index_ss];

	  start_el           = proc_ss_elem_len;

	  begin_df           = mono->ss_distfact_index[index_ss];

	  start_df	     = proc_ss_distfact_len;


	  

	  at_least_one       = FALSE; /* assume this ssid is not here 
				       * unless contraindicated
				       */
#ifdef DEBUG
	  fprintf(stderr, 
		  "Checking %d elem/sides for global sideset SS[%d] = %d\n\n",
		  mono->ss_num_sides[index_ss], index_ss, 
		  mono->ss_id[index_ss]);
#endif
	  /*
	   * New quick paradigm:  if the element is in the processor's
	   * list of elements then we want it.
	   */

	  for ( j=0; j<mono->ss_num_sides[index_ss]; j++)
	    {
	      elem                 = mono->ss_elem_list[begin+j];
	      side                 = mono->ss_side_list[begin+j];
#ifdef DEBUG
	      fprintf(stderr, "Checking: elem=(%d) side=%d j=%d ",
		      elem+1, side, j);
#endif

	      local_elem = in_list(elem, proc_elems, proc_ne);

	      if ( local_elem > -1 )
		{
#ifdef DEBUG
		  fprintf(stderr, 
			  "Adding. global SS[%d]=%d, elem=(%d), sd=%d !\n",
			  index_ss, mono->ss_id[index_ss], elem+1,
			  mono->ss_side_list[begin+j]);
		  fprintf(stderr,
			  "Local element index = %d\n", local_elem);
#endif
		  at_least_one = TRUE;

		  /*
		   * Now, don't forget where you found this!
		   */
		  
		  proc_ss_elem_list_index_global[proc_ss_elem_len] = 
		    begin+j;

		  proc_ss_elem_list[proc_ss_elem_len] = local_elem;
		  proc_ss_side_list[proc_ss_elem_len] = side;

		  proc_ss_elem_len++;

		  if ( proc_ss_elem_len > size_proc_ss_elem_list - 1 )
		    {
		      size_proc_ss_elem_list += 
			INIT_PROC_SS_SIDE_LIST_LENGTH;
		      
		      proc_ss_elem_list = 
			(int *) realloc(proc_ss_elem_list, 
					SZ_INT*size_proc_ss_elem_list);
		      
		      proc_ss_elem_list_index_global = 
			(int *) realloc(proc_ss_elem_list_index_global,
					SZ_INT*size_proc_ss_elem_list);
		      
		      proc_ss_side_list = 
			(int *) realloc(proc_ss_side_list, 
					SZ_INT*size_proc_ss_elem_list);
		    }

		  /*
		   * If this global side set has distribution factors
		   * then tack their values into an accumulating local
		   * side set distribution factor list...
		   */
		  if ( mono->ss_num_distfacts[index_ss] > 0 )
		    {
		      for ( l=mono->ss_node_side_index[index_ss][j]; 
			    l<mono->ss_node_side_index[index_ss][j+1]; l++)
			{
			  proc_ss_distfact_list[proc_ss_distfact_len] =
			    mono->ss_distfact_list[ begin_df + l ];

			  /*
			   * And now record where we got it, too!
			   */
			  
			  proc_ss_distfact_list_index_global
			    [proc_ss_distfact_len] = begin_df + l;

			  proc_ss_distfact_len++;

			  if ( proc_ss_distfact_len > 
			       size_proc_ss_distfact_list - 1 )
			    {
			      size_proc_ss_distfact_list +=
				INIT_PROC_SS_DISTFACT_LIST_LENGTH;
			      
			      proc_ss_distfact_list = (dbl *)
				realloc(proc_ss_distfact_list, SZ_DBL * 
					size_proc_ss_distfact_list);
			      
			      proc_ss_distfact_list_index_global = (int *) 
				realloc(proc_ss_distfact_list_index_global,
					SZ_INT*size_proc_ss_distfact_list);
			    }
			}
		    }
		}
	    }

	  if ( at_least_one )	/* mark this ssid  */
	    {
	      /*
	       * This "I" is purely a brief convenience.
	       */
	      I                         = proc_num_side_sets;
	      
	      proc_ss_id[I]             = mono->ss_id[index_ss];

	      proc_ss_num_sides[I]      = proc_ss_elem_len     - start_el;

	      proc_ss_num_distfacts[I]  = proc_ss_distfact_len - start_df;

	      proc_ss_elem_index[I]     = start_el;

	      proc_ss_distfact_index[I] = start_df;

#ifdef DEBUG
	      fprintf(stderr, "At least one elem/side from SSID=%d was\n",
		      mono->ss_id[index_ss]);
	      fprintf(stderr, "found for this set/proc. Looks like:\n");
	      fprintf(stderr, "\tnum_sides = %d (vs %d globally)\n",
		      proc_ss_num_sides[I], mono->ss_num_sides[index_ss]);
	      fprintf(stderr, "\tnum_distfacts = %d (vs %d globally)\n",
		      proc_ss_num_distfacts[I], 
		      mono->ss_num_distfacts[index_ss]);
	      fprintf(stderr, "\telement, distfact indeces at %d, %d\n",
		      proc_ss_elem_index[I], proc_ss_distfact_index[I]);
#endif

	      proc_num_side_sets++;

	      if ( proc_num_side_sets > size_proc_num_ss - 1 )
		{
		  size_proc_num_ss      += INIT_PROC_NUM_SS;

		  proc_ss_id             = ( (int *)
					     realloc(proc_ss_id, SZ_INT*
						     size_proc_num_ss));

		  proc_ss_num_sides      = ( (int *)
					     realloc(proc_ss_num_sides, SZ_INT*
						     size_proc_num_ss));

		  proc_ss_num_distfacts  = ( (int *)
					     realloc(proc_ss_num_distfacts,
						     SZ_INT *
						     size_proc_num_ss));

		  proc_ss_elem_index     = ( (int *)
					     realloc(proc_ss_elem_index, SZ_INT*
						     size_proc_num_ss));

		  proc_ss_distfact_index = ( (int *)
					     realloc(proc_ss_distfact_index,
						     SZ_INT *
						     size_proc_num_ss));
		}
	    }

	}

#ifdef DEBUG
      fprintf(stderr, "Assigning this set/proc's sideset concepts to E...\n");
#endif


      E->num_side_sets     = proc_num_side_sets;

      E->ss_elem_len       = proc_ss_elem_len;

      E->ss_distfact_len   = proc_ss_distfact_len;

      if ( E->num_side_sets > 0 )
	{
	  E->ss_id             = (int *) smalloc(E->num_side_sets*SZ_INT);

	  E->ss_num_sides      = (int *) smalloc(E->num_side_sets*SZ_INT);

	  E->ss_num_distfacts  = (int *) smalloc(E->num_side_sets*SZ_INT);

	  E->ss_elem_index     = (int *) smalloc(E->num_side_sets*SZ_INT);

	  E->ss_distfact_index = (int *) smalloc(E->num_side_sets*SZ_INT);
	}

      if ( E->ss_elem_len > 0 )
	{
	  E->ss_elem_list      = (int *) smalloc(E->ss_elem_len*SZ_INT);

	  E->ss_side_list      = (int *) smalloc(E->ss_elem_len*SZ_INT);
	}

      if ( E->ss_distfact_len > 0 )
	{
	  E->ss_distfact_list = (dbl *) smalloc(E->ss_distfact_len*SZ_DBL);
	}

      for( i=0; i<E->num_side_sets; i++)
	{
	  E->ss_id[i]             = proc_ss_id[i];

	  E->ss_num_sides[i]      = proc_ss_num_sides[i];

	  E->ss_num_distfacts[i]  = proc_ss_num_distfacts[i];

	  E->ss_elem_index[i]     = proc_ss_elem_index[i];

	  E->ss_distfact_index[i] = proc_ss_distfact_index[i];
	}
      
      for ( i=0; i<E->ss_elem_len; i++)
	{
	  E->ss_elem_list[i] = proc_ss_elem_list[i];

	  E->ss_side_list[i] = proc_ss_side_list[i];
	}

      for ( i=0; i<E->ss_distfact_len; i++)
	{
	  E->ss_distfact_list[i] = proc_ss_distfact_list[i];
	}

      /*
       * Now, we can put together the elem_maps, node_maps, (dof_maps)
       *
       * proc_elem[] *is* the elem_map.
       *
       * The node map is internal,boundary,external.
       */

      /*
       * Element Blocks:
       */

      /*
       * New name of part: file.exoII becomes file.exoII.6.2...
       */
      

      E->path            = (char *) smalloc(FILENAME_MAX_ACK*SZ_CHR);

      sr                 = sprintf(E->path, "%s.%d.%d", in_exodus_file_name, num_pieces, s);

      E->title           = (char *) smalloc(MAX_LINE_LENGTH*SZ_CHR);

      tmp                = strcpy(E->title, mono->title);
      
      E->comp_wordsize   = SZ_DBL;

      E->io_wordsize     = SZ_DBL;

      E->num_dim         = mono->num_dim;

      E->num_nodes       = num_universe_nodes;

      E->num_elems       = proc_ne;

      E->num_elem_blocks = proc_neb;

      E->num_times       = mono->num_times;

      E->api_version     = mono->api_version;

      E->db_version      = mono->db_version;

      E->version         = mono->version;

      E->num_info        = mono->num_info;

      if ( E->num_info > 0 )
	{
	  E->info            = (INFO_Record *) smalloc(E->num_info*
						       sizeof(INFO_Record));

	  for ( i=0; i<E->num_info; i++)
	    {
	      E->info[i]     = (char *) smalloc((MAX_LINE_LENGTH+1)*SZ_CHR);
	      tmp            = strcpy(E->info[i], mono->info[i]);
	    }
	}

      E->num_qa_rec = mono->num_qa_rec;
      
      if ( E->num_qa_rec > 0 )
	{
	  E->qa_record = (QA_Record *)smalloc(E->num_qa_rec*sizeof(QA_Record));
	  for ( i=0; i<E->num_qa_rec; i++)
	    {
	      for ( j=0; j<4; j++)
		{
		  E->qa_record[i][j] = (char *)smalloc(MAX_STR_LENGTH*SZ_CHR);
		  tmp = strcpy(E->qa_record[i][j], mono->qa_record[i][j]);
		}
	    }
	}


      /*
       * Coordinate names.
       */

      E->coord_names = (char **)smalloc(E->num_dim*SZPCHR);

      for ( i=0; i<E->num_dim; i++ )
	{
	  E->coord_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	  tmp = strcpy(E->coord_names[i], mono->coord_names[i]);
	}

      /*
       * Node coordinates.
       */

      if ( E->num_dim > 0 )
	{
	  E->x_coord = (dbl *) smalloc(num_universe_nodes*SZ_DBL);
	  for ( i=0; i<num_universe_nodes; i++)
	    {
	      node          = proc_nodes[i];
	      E->x_coord[i] = mono->x_coord[node];
	    }
	}

      if ( E->num_dim > 1 )
	{
	  E->y_coord = (dbl *) smalloc(num_universe_nodes*SZ_DBL);
	  for ( i=0; i<num_universe_nodes; i++)
	    {
	      node          = proc_nodes[i];
	      E->y_coord[i] = mono->y_coord[node];
	    }
	}
      
      if ( E->num_dim > 2 )
	{
	  E->z_coord = (dbl *) smalloc(num_universe_nodes*SZ_DBL);
	  for ( i=0; i<num_universe_nodes; i++)
	    {
	      node          = proc_nodes[i];
	      E->z_coord[i] = mono->z_coord[node];
	    }
	}
      
      /*
       * Element Block data.
       */

      E->eb_id                 = (int *) smalloc(proc_neb*SZ_INT);

      E->eb_elem_type          = (char **) smalloc(proc_neb*SZPCHR);

      E->eb_num_elems          = (int *) smalloc(proc_neb*SZ_INT);

      E->eb_num_nodes_per_elem = (int *) smalloc(proc_neb*SZ_INT);

      E->eb_num_attr           = (int *) smalloc(proc_neb*SZ_INT);

      E->eb_conn               = (int **) smalloc(proc_neb*SZPINT);

      E->eb_attr               = (dbl **) smalloc(proc_neb*SZPDBL);

      E->eb_ptr                = (int *) smalloc((proc_neb+1)*SZ_INT);

      for ( i=0; i<proc_neb; i++)
	{
	  E->eb_elem_type[i]   = (char *) smalloc((MAX_STR_LENGTH+1)*SZ_CHR);
	}
      

      for ( i=0; i<proc_neb+1; i++)
	{
	  E->eb_ptr[i] = proc_eb_ptr[i];
	}

      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  E->eb_id[ieb]                 = proc_eb_id[ieb];

	  eb_index = in_list(E->eb_id[ieb], mono->eb_id, mono->num_elem_blocks);

	  strcpy(E->eb_elem_type[ieb], mono->eb_elem_type[eb_index]);

	  E->eb_num_nodes_per_elem[ieb] = mono->eb_num_nodes_per_elem[eb_index];

	  E->eb_num_attr[ieb]           = mono->eb_num_attr[eb_index];

	  E->eb_num_elems[ieb]          = proc_eb_ptr[ieb+1] - proc_eb_ptr[ieb];
	}

      for ( i=0; i<proc_neb; i++)
	{
	  if ( (E->eb_num_elems[i] * E->eb_num_nodes_per_elem[i]) > 0 )
	    {
	      E->eb_conn[i] = (int *) smalloc(E->eb_num_elems[i]*
					      E->eb_num_nodes_per_elem[i]*
					      SZ_INT);
	    }

	  if ( (E->eb_num_elems[i]*E->eb_num_attr[i]) > 0 )
	    {
	      E->eb_attr[i] = (dbl *) smalloc(E->eb_num_elems[i]* 
					      E->eb_num_attr[i]*
					      SZ_DBL);
	    }
	}


      /*
       * Construct connectivity and attributes.
       */

      for ( ieb=0; ieb<proc_neb; ieb++)
	{
	  eb_index = in_list(E->eb_id[ieb], mono->eb_id, mono->num_elem_blocks);

	  begin  = proc_eb_ptr[ieb];

	  length = proc_eb_ptr[ieb+1] - begin;

	  /*
	   * The number of element attributes should be the same from the
	   * global and local perspective...
	   */

	  gnattr = mono->eb_num_attr[eb_index];

	  lnattr = E->eb_num_attr[ieb];

	  count  = 0;

	  gbeg   = ebl[eb_index];

	  for ( e=0; e<length; e++)
	    {
	      elem   = proc_elems[begin+e]; /* global elem name */
	      e_name = elem - gbeg; /* from the global eb perspective */
	      for ( l=ep[elem]; l<ep[elem+1]; l++)
		{
		  node = nl[l];
		  /*
		   * Look for this node in the 3 parts of the proc_node list
		   * each of which is monotonic - the entirety is not 
		   * necessarily...
		   */

		  n = get_node_index(node, proc_nodes, num_internal_nodes,
				     num_boundary_nodes, num_external_nodes);
#ifdef DEBUG
		  fprintf(stderr, 
			  "gieb=%d, gelem=%d, (eb:%d), gnode=(%d)\n, n=[%d]\n",
			  eb_index, elem, e_name, node+1, n);
#endif
		  E->eb_conn[ieb][count] = n;		  
		  count++;
		}
	      
	      for ( a=0; a<lnattr; a++)
		{
		  E->eb_attr[ieb][e*lnattr+a] = 
		    mono->eb_attr[eb_index][e_name*gnattr+a];
		}
	    }
	}
      
      
      /*
       * Properties...
       *
       * Hmmm... Looks like the ID property might be hardwired into the
       * EXODUS II API and not need to be explicitly copied.
       *
       * The number of nodeset properties for the set/proc will be
       * equal to the number of nodeset properties in the global problem
       * for nodesets whose id's are in this set/procs' inventory.
       *
       * Similar principles govern the presence and ordering of
       * properties for sidesets and elementblocks.
       */

      E->ns_num_props  = mono->ns_num_props;

      if ( E->ns_num_props > 1 && 
	   E->num_node_sets > 0 &&
	   strcmp(mono->ns_prop_name[0], mono->ns_prop_name[1]) != 0 )
	{
	  E->ns_prop_name  = (char **) smalloc(E->ns_num_props*SZPCHR);

	  for ( i=0; i<E->ns_num_props; i++)
	    {
	      E->ns_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->ns_prop_name[i], mono->ns_prop_name[i]);
	    }

	  E->ns_prop = (int **) smalloc(E->ns_num_props*SZPINT);

	  for ( i=0; i<E->ns_num_props; i++)
	    {
	      E->ns_prop[i] = (int *)smalloc(E->num_node_sets*SZ_INT);
	    }

	  for ( ins=0; ins<E->num_node_sets; ins++)
	    {
	      g_ns_index = in_list(E->ns_id[ins], mono->ns_id, 
				   mono->num_node_sets);

	      for ( p=0; p<E->ns_num_props; p++)
		{
		  E->ns_prop[p][ins] = mono->ns_prop[p][g_ns_index];
		}
	    }
	}

      E->ss_num_props = mono->ss_num_props;

      if ( E->ss_num_props > 1 && 
	   E->num_side_sets > 0 &&
	   strcmp(mono->ss_prop_name[0], mono->ss_prop_name[1]) != 0 )
	{
	  E->ss_prop_name = (char **) smalloc(E->ss_num_props*SZPCHR);

	  for ( i=0; i<E->ss_num_props; i++)
	    {
	      E->ss_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->ss_prop_name[i], mono->ss_prop_name[i]);
	    }
	  
	  E->ss_prop = (int **) smalloc(E->ss_num_props*SZPINT);

	  for ( i=0; i<E->ss_num_props; i++)
	    {
	      E->ss_prop[i] = (int *)smalloc(E->num_side_sets*SZ_INT);
	    }
	  
	  for ( iss=0; iss<E->num_side_sets; iss++)
	    {
	      g_ss_index = in_list(E->ss_id[iss], mono->ss_id,
				   mono->num_side_sets);
	      for ( p=0; p<E->ss_num_props; p++)
		{
		  E->ss_prop[p][iss] = mono->ss_prop[p][g_ss_index];
		}
	    }
	}

      /*
       * The number of properties in for EXODUS II entities is very
       * slippery. The API will often fabricate one property anyway.
       * That's the "ID" property. Since it does this automatically
       * we have to be very careful not to replicate that one property.
       */


      E->eb_num_props = mono->eb_num_props;

      if ( E->eb_num_props > 1 && 
	   E->num_elem_blocks > 0 &&
	   strcmp(mono->eb_prop_name[0], mono->eb_prop_name[1]) != 0 )
	{
	  E->eb_prop_name = (char **) smalloc(E->eb_num_props*SZPCHR);

	  for ( i=0; i<E->eb_num_props; i++)
	    {
	      E->eb_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->eb_prop_name[i], mono->eb_prop_name[i]); 
	    }

	  E->eb_prop = (int **) smalloc(E->eb_num_props*SZPINT);
	  for ( i=0; i<E->eb_num_props; i++)
	    {
	      E->eb_prop[i] = (int *)smalloc(E->num_elem_blocks*SZ_INT);
	    }

	  for ( ieb=0; ieb<E->num_elem_blocks; ieb++)
	    {
	      g_eb_index = in_list(E->eb_id[ieb], mono->eb_id, 
				   mono->num_elem_blocks);

	      for ( p=0; p<E->eb_num_props; p++)
		{
		  E->eb_prop[p][ieb] = mono->eb_prop[p][g_eb_index];
		}
	      
	    }

	}

      /*
       * Build up auxiliary connectivity information for the little EXODUS II
       * piece that will be useful in construction of dpi->elem_elem_...
       */

        build_elem_node(E);
	build_node_elem(E);
	build_elem_elem(E);
	build_node_node(E);

      /*
       * Write the EXODUS II finite element MESH database information for this
       * set/proc piece of the global problem.
       */

      one_base(E);
      wr_mesh_exo(E, E->path, 0);
      zero_base(E);



      /*
       * Maybe we ought to wr_dpi() here, before the unlimited records start.
       */




      /*
       * Setup Distributed Processing Information...
       */

      /*
       * First, create a local version of the set membership stuff...
       */

      proc_psm = (int *) smalloc((num_universe_nodes+1)*SZ_INT);
      count    = 0;
      for ( n=0; n<num_universe_nodes; n++)
	{
	  proc_psm[n] = count;
	  node        = proc_nodes[n];
	  count      += (psm[node+1]-psm[node]);
	}
      proc_psm[num_universe_nodes] = count;


      proc_sm  = (int *) smalloc(count*SZ_INT);

      for ( n=0; n<num_universe_nodes; n++)
	{
	  node = proc_nodes[n];
	  len  = proc_psm[n+1]-proc_psm[n];
	  for ( i=0; i<len; i++)
	    {
	      proc_sm[proc_psm[n]+i] = sm[psm[node]+i];
	    }
	}

      /*
       * (Don't look at membership lists of the external nodes and their
       *  neighbors, except for the first one that indicates the primary
       *  owner of the node. Any remaining entries can indicate who owns
       *  nodes that this processor knows nothing about.
       */

      D->my_name               = s+1;

      if ( num_send_procs != num_recv_procs )
	{
	  sr = sprintf(err_msg, "! num_send_procs = %d, num_recv_procs = %d.",
		       num_send_procs, num_recv_procs);
	  EH(-1, err_msg);
	}

      D->num_neighbors                = num_send_procs;

      D->neighbor                     = (int *) 
	smalloc(D->num_neighbors*sizeof(int));
      memcpy(D->neighbor, send_proc_names, num_send_procs*sizeof(int));


      D->len_set_membership           = proc_psm[num_universe_nodes];

      D->set_membership               = (int *) 
	smalloc(D->len_set_membership*sizeof(int));

      for ( i=0; i<D->len_set_membership; i++)
	{
	  D->set_membership[i] = proc_sm[i];
	}

      D->len_ptr_set_membership       = num_universe_nodes+1;
      
      D->ptr_set_membership           = (int *) 
	smalloc(D->len_ptr_set_membership*sizeof(int));

      for ( i=0; i<num_universe_nodes+1; i++)
	{
	  D->ptr_set_membership[i]           = proc_psm[i];
	}
      

      D->global_node_dof0        = (int *) smalloc(num_universe_nodes*SZ_INT);

      for ( n=0; n<num_universe_nodes; n++)
	{
	  node = proc_nodes[n];
	  D->global_node_dof0[n]  = node_dof0[node];
	}

      D->global_node_kind         = (int *) smalloc(num_universe_nodes*SZ_INT);

      for ( n=0; n<num_universe_nodes; n++)
	{
	  node = proc_nodes[n];
	  D->global_node_kind[n]  = node_kind[node];
	}

      D->elem_index_global = (int *) smalloc(proc_ne*SZ_INT);
      for ( i=0; i<proc_ne; i++)
	{
	  D->elem_index_global[i] = proc_elems[i];
	}

      D->len_eb_num_private_elems = proc_neb;

      D->eb_num_private_elems     = (int *) smalloc(proc_neb*SZ_INT);

      for ( i=0; i<E->num_elem_blocks; i++)
	{
	  D->eb_num_private_elems[i] = private_elem_count[i];
	}

      D->num_internal_nodes       = num_internal_nodes;

      D->num_boundary_nodes       = num_boundary_nodes;

      D->num_external_nodes       = num_external_nodes;
      
      D->num_universe_nodes       = num_universe_nodes;

      /*
       * These are handy quick maps. NOTE these {eb,ns,ss}_index_global[]
       * arrays are dimensioned according to how many {eb,ns,ss} there are
       * on this set/processor. The other global arrays are dimensioned
       * according to the global number of {eb,ns,ss} entities that exist.
       *
       * Changes: While not a problem so much for element blocks, there
       * is definitely the possibility that a given processor will have
       * zero nodesets and/or zero sidesets. This has implications for
       * memory allocation, etc.
       */

      D->num_elem_blocks = proc_neb;

      if ( D->num_elem_blocks > 0 )
	{
	  D->eb_index_global = (int *) smalloc(D->num_elem_blocks*SZ_INT);
	  for ( i=0; i<D->num_elem_blocks; i++)
	    {
	      D->eb_index_global[i] = in_list(proc_eb_id[i], mono->eb_id, neb);
	    }
	}

      D->num_node_sets   = proc_num_node_sets;

      if ( D->num_node_sets > 0 )
	{
	  D->ns_index_global = (int *) smalloc(D->num_node_sets*SZ_INT);
	  for ( i=0; i<D->num_node_sets; i++)
	    {
	      D->ns_index_global[i] = in_list(proc_ns_id[i], mono->ns_id, nns);
	    }
	}
	  
      D->num_side_sets   = proc_num_side_sets;

      if ( D->num_side_sets > 0 )
	{
	  D->ss_index_global = (int *) smalloc(D->num_side_sets*SZ_INT);
	  for ( i=0; i<D->num_side_sets; i++)
	    {
	      D->ss_index_global[i] = in_list(proc_ss_id[i], mono->ss_id, nss);
	    }
	}

      /*
       * Setup some additional information helpful for fix, etc...
       */

      D->len_elem_var_tab_global = mono->num_elem_vars * mono->num_elem_blocks;

      D->len_ns_node_list        = E->ns_node_len;

      D->len_ns_distfact_list    = E->ns_distfact_len;

      D->len_ss_elem_list        = E->ss_elem_len;

      D->len_ss_distfact_list    = E->ss_distfact_len;

      D->len_string              = MAX_STR_LENGTH + 1; /* same overage that
							* EXODUS II likes to
							* claim! */
      D->num_elems               = proc_ne;

      D->num_nodes               = num_universe_nodes;

      D->num_props_eb            = mono->eb_num_props;

      D->num_props_ns            = mono->ns_num_props;

      D->num_props_ss            = mono->ss_num_props;

      /*
       * New stuff for discontinous Galerkin methods requiring element
       * to element connectivity information...
       */

      

      D->len_elem_elem_list      = E->elem_elem_pntr[E->num_elems];


      /*
       * xyzzy
       * Look here! for beginning of allocation...so that you may properly
       * free!
       *
       * char **eb_elem_type_global [num_elem_blocks_global] [MAX_STR_LENGTH]
       */

      D->eb_elem_type_global = (char **) smalloc(D->num_elem_blocks_global*
						 sizeof(char *));

      D->eb_num_attr_global  = (int *) smalloc(D->num_elem_blocks_global*
					       sizeof(char *));

      D->eb_num_nodes_per_elem_global = 
	(int *) smalloc(D->num_elem_blocks_global*sizeof(char *));


      /*
       * Ease the transport of 2D arrays down into the bowels of netCDF
       * by allocating a contiguous chunk. Then, set the pointers
       * accordingly...
       */

      D->eb_elem_type_global[0] = (char *) smalloc((D->num_elem_blocks_global)*
						   (D->len_string)*
						   sizeof(char));

      for ( i=0; i<(D->num_elem_blocks_global)*(D->len_string); i++)
	{
	  D->eb_elem_type_global[0][i] = '\0';
	}

      for ( i=1; i<D->num_elem_blocks_global; i++)
	{
	  D->eb_elem_type_global[i] = D->eb_elem_type_global[i-1] +
				      D->len_string;
	}

      for ( i=0; i<mono->num_elem_blocks; i++)
	{
	  /*
	   * D->eb_elem_type_global[i] = (char *) smalloc((MAX_STR_LENGTH+1)*
	   *				       sizeof(char));
	   *
	   *	  for ( j=0; j<MAX_STR_LENGTH+1; j++)
	   *	    {
	   *	      D->eb_elem_type_global[i][j] = '\0';
	   *	    }
	   */

	  strcpy(D->eb_elem_type_global[i], mono->eb_elem_type[i]);

#ifdef DEBUG
  fprintf(stderr, "elem type is now ...%s\n", mono->eb_elem_type[i]);
#endif

	  D->eb_num_attr_global[i]           = mono->eb_num_attr[i];

	  D->eb_num_nodes_per_elem_global[i] = mono->eb_num_nodes_per_elem[i];
	}

      
      if ( mono->eb_num_props > 1 )
	{
	  D->eb_prop_global 
	    = (int **) smalloc(mono->eb_num_props*sizeof(int *));

	  D->eb_prop_global[0] = (int *) smalloc(mono->eb_num_props*
						 mono->num_elem_blocks*
						 sizeof(int));

	  for ( i=1; i<mono->eb_num_props; i++)	  
	    {
	      D->eb_prop_global[i] = ( D->eb_prop_global[i-1] + 
				       mono->num_elem_blocks );
	    }

	  for ( i=0; i<mono->eb_num_props; i++)
	    {
	      for ( j=0; j<mono->num_elem_blocks; j++)
		{
		  D->eb_prop_global[i][j] = mono->eb_prop[i][j];
		}
	    }
	}

      if ( D->len_elem_var_tab_global > 0 )
	{
	  D->elem_var_tab_global = (int *) smalloc(D->len_elem_var_tab_global*
						   sizeof(int));
	}

      for ( i=0; i<D->len_elem_var_tab_global; i++)
	{
	  D->elem_var_tab_global[i] = mono->elem_var_tab[i];
	}

      /* Hey! It's flat anyway !
      for ( eb=0; eb<mono->num_elem_blocks; eb++)
	{
	  for ( ev=0; ev<mono->num_elem_vars; ev++)
	    {
	      D->elem_var_tab_global[eb*mono->num_elem_vars+ev] =
		mono->elem_var_tab[eb*mono->num_elem_vars+ev];
	    }
	}
	*/


      /*
       * Node number map.
       */

      if ( num_universe_nodes > 0 )
	{
	  D->node_index_global = (int *) smalloc(num_universe_nodes*SZ_INT);
	  for ( i=0; i<num_universe_nodes; i++)
	    {
	      D->node_index_global[i] = proc_nodes[i];
	    }
	}
      
      D->ns_distfact_len_global = mono->ns_distfact_len;

      D->ns_node_len_global     = mono->ns_node_len;

      D->ns_node_index_global   = (int *) smalloc(mono->num_node_sets*
						  sizeof(int));
      for ( i=0; i<mono->num_node_sets; i++)
	{
	  D->ns_node_index_global[i] = mono->ns_node_index[i];
	}

      D->ns_distfact_index_global   = (int *) smalloc(mono->num_node_sets*
						  sizeof(int));
      for ( i=0; i<mono->num_node_sets; i++)
	{
	  D->ns_distfact_index_global[i] = mono->ns_distfact_index[i];
	}

      
      if ( mono->ns_num_props > 1 )
	{
	  D->ns_prop_global = 
	    (int **) smalloc(mono->ns_num_props*sizeof(int *));

	  D->ns_prop_global[0] = (int *) smalloc(mono->ns_num_props*
						 mono->num_node_sets*
						 sizeof(int));
	  for ( i=1; i<mono->ns_num_props; i++)
	    {
	      D->ns_prop_global[i] = ( D->ns_prop_global[i-1] +
				       mono->num_node_sets );
	    }

	  for ( i=0; i<mono->ns_num_props; i++)
	    {
	      for ( j=0; j<mono->num_node_sets; j++)
		{
		  D->ns_prop_global[i][j] = mono->ns_prop[i][j];
		}
	    }
	}

      /*
       * OK, here's a tricky one. For each node in the concatenated node
       * list of nodesets for this little set/proc, we want to record the
       * index in the global concatenated node list.
       */

      D->ns_node_list_index_global = (int *) smalloc(E->ns_node_len*
						     sizeof(int));

      for ( i=0; i<E->ns_node_len; i++)
	{
	  D->ns_node_list_index_global[i] = proc_ns_node_list_index_global[i];
	}
      
      /*
       * Record the global indeces for the nodeset distribution factors that
       * were collected as well...
       */

      D->ns_distfact_list_index_global = (int *) smalloc(E->ns_distfact_len*
							 sizeof(int));
      for ( i=0; i<E->ns_distfact_len; i++)
	{
	  D->ns_distfact_list_index_global[i] = 
	    proc_ns_distfact_list_index_global[i];
	}


      /*
       * Side set information from the global problem that is indispensible
       * for later reconstruction...
       */

      D->ss_distfact_len_global = mono->ss_distfact_len;

      D->ss_elem_len_global     = mono->ss_elem_len;

      D->ss_distfact_index_global = (int *) smalloc(mono->num_side_sets*
						    sizeof(int));
      D->ss_elem_index_global = (int *) smalloc(mono->num_side_sets*
						    sizeof(int));
      for ( i=0; i<mono->num_side_sets; i++)
	{
	  D->ss_distfact_index_global[i] = mono->ss_distfact_index[i];
	  D->ss_elem_index_global[i]     = mono->ss_elem_index[i];
	}

      D->ss_elem_list_index_global = (int *) smalloc(E->ss_elem_len*
						     sizeof(int));
      D->ss_distfact_list_index_global = (int *) smalloc(E->ss_distfact_len*
							 sizeof(int));

      for ( i=0; i<E->ss_distfact_len; i++)
	{
	  D->ss_distfact_list_index_global[i] = 
	    proc_ss_distfact_list_index_global[i];
	}

      for ( i=0; i<E->ss_elem_len; i++)
	{
	  D->ss_elem_list_index_global[i] = 
	    proc_ss_elem_list_index_global[i];
	}

      if ( mono->ss_num_props > 1 )
	{
	  D->ss_prop_global = 
	    (int **) smalloc(mono->ss_num_props*sizeof(int *));

	  D->ss_prop_global[0] = (int *) smalloc(mono->ss_num_props*
						 mono->num_side_sets*
						 sizeof(int));
	  for ( i=1; i<mono->ss_num_props; i++)
	    {
	      D->ss_prop_global[i] = ( D->ss_prop_global[i-1] +
				       mono->num_side_sets );
	    }

	  for ( i=0; i<mono->ss_num_props; i++)
	    {
	      for ( j=0; j<mono->num_side_sets; j++)
		{
		  D->ss_prop_global[i][j] = mono->ss_prop[i][j];
		}
	    }
	}

#ifdef DEBUG
      for ( i=0; i<mono->eb_num_props; i++)
	{
	  for ( j=0; j<mono->num_elem_blocks; j++)
	    {
	      fprintf(stderr, "monolith eb props[%d][%d] = %d\n", i, j, 
		      mono->eb_prop[i][j]);
	    }
	}
      for ( i=0; i<mono->ns_num_props; i++)
	{
	  for ( j=0; j<mono->num_node_sets; j++)
	    {
	      fprintf(stderr, "monolith ns props[%d][%d] = %d\n", i, j, 
		      mono->ns_prop[i][j]);
	    }
	}
      for ( i=0; i<mono->ss_num_props; i++)
	{
	  for ( j=0; j<mono->num_side_sets; j++)
	    {
	      fprintf(stderr, "monolith ss props[%d][%d] = %d\n", i, j, 
		      mono->ss_prop[i][j]);
	    }
	}
#endif

      /*
       * Wait til now to build the dpi->elem_elem stuff since it depends
       * on having dpi->elem_index_global[] for example.
       */

      build_elem_elem_dpi(mono, element_owner, E, D);

      build_elem_elem_xtra(E);

      wr_dpi(D, E->path, 0);

      /*
       * Back to our regularly schedule transcription of E->results data
       * into the polyliths...
       */

      /*
       * Results data. This includes both the preliminary descriptions
       * such as nodal variable names as well as the vectors of data per se.
       *
       * Nodal variables are presumed to exist at every node by EXODUS II,
       * thus they will exist in a subset for each set/proc. Element variables
       * are more problematic, depending on the element truth table and which
       * element blocks are present for this processor, etc.
       *
       * Thus, for now, any element variables that exist in the monolithic
       * database are NOT transcribed into the polyliths.
       */

      E->num_glob_vars = mono->num_glob_vars;
      E->num_node_vars = mono->num_node_vars;
      E->num_elem_vars = mono->num_elem_vars;
      
      if ( E->num_glob_vars > 0 )
	{
	  E->glob_var_names = (char **) smalloc(E->num_glob_vars*SZPCHR);
	  for ( i=0; i<E->num_glob_vars; i++)
	    {
	      E->glob_var_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->glob_var_names[i], mono->glob_var_names[i]);
	    }

	}

      if ( E->num_node_vars > 0 )
	{
	  E->node_var_names = (char **) smalloc(E->num_node_vars*SZPCHR);
	  for ( i=0; i<E->num_node_vars; i++)
	    {
	      E->node_var_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->node_var_names[i], mono->node_var_names[i]);
	    }
	}

      if ( E->num_elem_vars > 0 )
	{
	  E->elem_var_names = (char **) smalloc(E->num_elem_vars*SZPCHR);
	  for ( i=0; i<E->num_elem_vars; i++)
	    {
	      E->elem_var_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	      tmp = strcpy(E->elem_var_names[i], mono->elem_var_names[i]);
	    }


	  E->elem_var_tab = (int *) smalloc(E->num_elem_vars*
					    E->num_elem_blocks*
					    sizeof(int));
	  /*
	   * Transcribe the element variable truth table from the monolith
	   * to the polylith for the element blocks that this polylith contains.
	   *
	   * Why even bother, we're not going to write it out!
	   */

	  nev = E->num_elem_vars; /* shorthand. Global monolith is same. */

	  for ( ieb=0; ieb<E->num_elem_blocks; ieb++)
	    {
	      g_eb_index = in_list(E->eb_id[ieb], mono->eb_id, 
				   mono->num_elem_blocks);
	      EH(g_eb_index, "Problem finding corresponding global eb index.");
	      for ( i=0; i<nev; i++)
		{
		  E->elem_var_tab[ieb*nev+i] = 
		    mono->elem_var_tab[g_eb_index*nev+i];
		}
	    }

	}

      /*
       * If there are multiple time planes of monolithic results, then
       * transcribe them into the polylith version, too.
       */

      if ( E->num_times > 0 )
	{
	  E->time_vals = (dbl *) smalloc(E->num_times*sizeof(dbl));

	  for ( i=0; i<E->num_times; i++)
	    {
	      E->time_vals[i] = mono->time_vals[i];
	    }
	}

      one_base(E);
      wr_resetup_exo(E, E->path, 0);
      zero_base(E);

      /*
       * If there are any nodal variables, then read their values from
       * each "timestep" and transcribe them into the right place
       * in this set/proc's repertoire.
       *
       * Now, we want to use "nv" instead, but still just do one step
       * at a time.
       *
       * Note! Element variables have not been done, only nodal results
       *       variables. To do element variables would be straightforward,
       *       but you do need to exercise caution when writing an augmented
       *       plot file with added element variables. Similar precautions
       *       were taken for nodal variables as you'll see below looking
       *       for how "nnv" and "len" are used...
       */

      if ( E->num_node_vars > 0 )
	{
	  E->num_nv_indeces      = E->num_node_vars;
	  E->nv_indeces          = (int *) smalloc(E->num_nv_indeces*
						   sizeof(int));
	  /*
	   * Even though there may be more than one time plane, we'll just
	   * walk through one at a time.
	   */

	  E->num_nv_time_indeces = 1;
	  E->nv_time_indeces     = (int *) smalloc(E->num_nv_time_indeces*
						   sizeof(int));
	  
	  E->state              |= EXODB_STATE_NDIA;

	  alloc_exo_nv(E);

	  /*
	   * This sweep through all the monolith's nodal results for each
	   * set/proc might be inefficient. Alternatives would require 
	   * having available the global node map for all the set/procs at the
	   * same time.
	   */

	  if ( mono->num_node_vars < 1 )
	    {
	      EH(-1, "Inconsistent nodal variable count?");
	    }

	  /* 
	   * These helper variables were set up for mono-> in the rd_exo
	   * module. The E-> memory representation was not read in, but
	   * synthesized from scratch, so it needs the rigmarole above...
	   */


	  if ( ! ( mono->state & EXODB_STATE_NDIA ) )
	    {
	      mono->num_nv_indeces      = mono->num_node_vars;
	      mono->nv_indeces          = (int *) smalloc(mono->num_nv_indeces*
						      sizeof(int));
	      mono->num_nv_time_indeces = 1;
	      mono->nv_time_indeces     = (int *) smalloc(mono->num_nv_time_indeces
						      * sizeof(int));
	      mono->state              |= EXODB_STATE_NDIA;

	    }

	  alloc_exo_nv(mono);

	  for ( i=0; i<mono->num_node_vars; i++)      
	    {
	      mono->nv_indeces[i]   = i+1;
	      E->nv_indeces[i]      = i+1;
	    }

	  for ( t=0; t<E->num_times; t++ )
	    {
	      /*
	       * Grab one monolith time plane for each nodal var...
	       */

	      mono->nv_time_indeces[0]  = t+1;
	      E->nv_time_indeces[0]     = t+1;

	      status = rd_exo(mono, in_exodus_file_name, 0, 
			      EXODB_ACTION_RD_RESN);

	      /*
	       * Pick out pieces for this set/proc...
	       */

	      for ( n=0; n<E->num_nodes; n++)
		{
		  node = proc_nodes[n];
		  if ( node < 0 || node > mono->num_nodes-1 )
		    {
		      EH(-1, "Bad map.");
		    }
		  for ( i=0; i<E->num_nv_indeces; i++)
		    {
		      E->nv[0][i][n] = mono->nv[0][i][node];
		    }
		}

	      /*
	       * Write out results.
	       */
	      
	      wr_result_exo(E, E->path, 0);
	    }

	  free_exo_nv(E);
	  free_exo_nv(mono);
	}
	
      if ( E->num_elem_vars > 0 )  /* Save me baby Jesus ! Imma try to clone the nv section and hope for the best */
	{
	  E->num_ev_indeces = E->num_elem_vars;
	  E->ev_indeces     = (int *) smalloc( E->num_ev_indeces * sizeof( int ) );
	  E->num_ev_time_indeces = 1;

	  /*
	   * Even though there may be more than one time plane, we'll just
	   * walk through one at a time.
	   */


	  E->state               |= EXODB_STATE_ELIA;

	  alloc_exo_ev (E,1);


	  if ( ! ( mono->state & EXODB_STATE_ELIA ) )
	    {
	      mono->num_ev_indeces      = mono->num_elem_vars;
	      mono->ev_indeces          = (int *) smalloc(mono->num_ev_indeces*sizeof(int));
	      mono->num_ev_time_indeces = 1;
	      mono->state              |= EXODB_STATE_ELIA;

	    }
      
	  alloc_exo_ev(mono,1);
      
	  if ( mono->num_elem_vars < 1 )
	    {
	      EH(-1, "Inconsistent element variable count?");
	    }


	  for ( i=0; i<mono->num_elem_vars; i++)      
	    {
	      mono->ev_indeces[i]   = i+1;
	      E->ev_indeces[i]      = i+1;
	    }


	  for ( t=0; t<E->num_times; t++ )
	    {
	      /*
	       * Grab one monolith time plane for each element var...
	       */

	      mono->ev_time_indeces[0]  = t+1;
	      E->ev_time_indeces[0]     = t+1;

	      status = rd_exo(mono, in_exodus_file_name, 0, EXODB_ACTION_RD_RESE);

	      /*
	       * Pick out pieces for this set/proc...
	       */

	      for ( ieb=0; ieb<E->num_elem_blocks; ieb++)
		{
		  eb_index = in_list(E->eb_id[ieb], mono->eb_id, mono->num_elem_blocks);
		  
		  gbeg = ebl[eb_index];

		  begin = E->eb_ptr[ieb];
		  
		  for( e=0; e<E->eb_num_elems[ieb]; e++)
		    {
		      elem = proc_elems[e + begin];
		  
		      if ( elem < 0 || elem > mono->num_elems-1 )
			{
			  EH(-1, "Bad map.");
			}

		      ge = elem - gbeg;

		      for ( i=0; i<E->num_elem_vars; i++)
			{
			  index = ieb*E->num_elem_vars + i;
			  g_eb_index = eb_index*mono->num_elem_vars + i;

			  if( E->elem_var_tab[index] != -1 )
			    {
			      E->ev[0][index][e] = mono->ev[0][g_eb_index][ge];
			    }
			}
		    }
		}

	      /*
	       * Write out results.
	       */
	      
	      wr_result_exo(E, E->path, 0);
	    }

	  free_exo_ev(E);
	  free_exo_ev(mono);
	}
	

      /*
       * Free up some arrays allocated during this set/proc loop...
       */

      free(internal_nodes);
      free(boundary_nodes);
      free(external_nodes);

      free(proc_nodes);
      free(proc_elems);

      free(proc_sm);
      free(proc_psm);

      free(private_elem_count);
      free(proc_elem_priv);
      free(proc_elem_shar);

      if ( D->num_node_sets > 0 )
	{
	  free(proc_ns_id);
	  free(proc_ns_num_nodes);
	  free(proc_ns_num_distfacts);
	  free(proc_ns_node_index);
	  free(proc_ns_distfact_index);
	  free(proc_ns_node_list);
	  free(proc_ns_node_list_index_global);
	  free(proc_ns_distfact_list);
	  free(proc_ns_distfact_list_index_global);
	}

      free(proc_ss_id);
      free(proc_ss_num_sides);
      free(proc_ss_num_distfacts);
      free(proc_ss_elem_index);
      free(proc_ss_distfact_index);
      free(proc_ss_elem_list);
      free(proc_ss_elem_list_index_global);
      free(proc_ss_side_list);
      free(proc_ss_distfact_list);
      free(proc_ss_distfact_list_index_global);

      /*
       * Build little versions of ep, nl, np, el, mult, bevm and output
       * node_descriptions with node_dof0 maps for this set/proc. 
       *
       * Then, on a node-by-node basis
       * compare the set/proc local node with the global perspective, which
       * might include additional unknowns at the node for external nodes. 
       *
       * Needed information for distributed processing:
       *
       *
       * num_proc_node_descriptions (int)
       *	proc_num_node_descriptions = number of distinct kinds of
       *				     nodes from the set/proc's
       *				     perspective.
       *
       *	proc_node_descriptions[]   = the descriptions themselves
       *
       *        proc_node_kind[pni=0,...,numnodes] = index into nodescriptions
       *					     for each node.
       *
       *	proc_node_dof0[proc_node_index=0...proc_num_nodes]
       *
       *		proc_node_dof0[local_node_number]
       *		global_node_name[local_node_number]
       *		global_node_dof0[local_node_number]
       *		global_node_kinds[]
       *		proc_node_kinds[]
       *
       * Obtain a map of the
       * global dof names corresponding to each local dof name. There can
       * be gaps at externalnodes. That is, the global perspective include
       * more dofs at the node than the setproc believes are appropriate.
       *
       * Also, write primary responsible set/proc for each global node.
       *
       */

      free_exo(E);
      free_dpi(D);

    } /* set loop */

  free(proc_eb_id);
  free(proc_eb_ptr);
  free(ebi);
  free(new_proc_eb_ptr);

  if ( mono->num_dim > 0 )
    {
      free(x);
    }
  if ( mono->num_dim > 1 )
    {
      free(y);
    }
  if ( mono->num_dim > 2 )
    {
      free(z);
    }

  free(D);

  /*
   * Replicate the monolithic EXODUS II file and augment the output
   * with information about the decomposition that was performed.
   */

  if ( add_decomp_plot_vars )
    {
      /*
       * Create a nodal variables that express for each node:
       *	(1) the number of nonzero matrix entries in the row
       *	(2) the number of assembled terms in the row
       *	(3) the communication cost with other nodes
       *	(4) the number of degrees of freedom or equations that
       *	    are associated with the node
       *	(5) the global node kind index
       *	(6) how many processors are interested in this node
       *	(7) the name of the processor that owns this node
       *
       * Allocate sufficient space for the new names, as well as any old
       * nodal variable names and set this up.
       */
      
      /*
       * These are the indeces for each of the new nodal post processing
       * variables.
       */

#ifdef DEBUG
      fprintf(stderr, "Adding decomposition plot vbls to %d already there.\n",
	      mono->num_node_vars);
#endif

      nnv        = mono->num_node_vars;

      index_nnz  = nnv + 0;
      index_nat  = nnv + 1;
      index_ccs  = nnv + 2;
      index_dofs = nnv + 3;
      index_nknd = nnv + 4;
      index_nscn = nnv + 5;
      index_nown = nnv + 6;

      len        = nnv + 7;

      /*
       * Allocate space for and read in old nodal results at the last timeplane
       * if they exist, then tack on our new results.
       */

      if ( mono->num_nv_time_indeces > 1 )
	{
	  EH(-1, 
	  "Reconcile new nodal results with multiple time previous results.");
	}

      /*
       * We want to boost the amount of allocated space for nodal variables
       * beyond what the meta data indicates from the original EXODUS II file.
       * After reading the existing nodal variables, we'll load in some new
       * results into the extra space we've allocated.
       */

      free_exo_nv(mono);

      mono->num_nv_time_indeces = 1;
      mono->num_nv_indeces      = nnv + 7;

      alloc_init_exo_nv_indeces(mono);

      alloc_exo_nv(mono);

      /*
       * Temporary retrograde values while we read in the existing nodal
       * results variables...
       */

      mono->num_nv_indeces = nnv; 

      status = rd_exo(mono, in_exodus_file_name, 0, EXODB_ACTION_RD_RESN);

      mono->num_nv_indeces = len; 
      mono->num_node_vars  = len;

      if ( nnv > 0 )
	{
	  mono->node_var_names = (char **) realloc(mono->node_var_names, 
						   len * sizeof(char *));
	}
      else
	{
	  mono->node_var_names = (char **) smalloc(len * sizeof(char *));
	}

      for ( i=nnv; i<len; i++)
	{
	  mono->node_var_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	}

      strcpy(mono->node_var_names[index_nnz],  "NNZ");
      strcpy(mono->node_var_names[index_nat],  "NAT");
      strcpy(mono->node_var_names[index_ccs],  "CCS");
      strcpy(mono->node_var_names[index_dofs], "DOFS");
      strcpy(mono->node_var_names[index_nknd], "GNKND");
      strcpy(mono->node_var_names[index_nscn], "NSCN");
      strcpy(mono->node_var_names[index_nown], "NOWN");

      /*
       * Now fill in the arrays of nodal values.
       */

      for ( i=0; i<nn; i++)
	{
	  sum_nnz = 0;
	  sum_nat = 0;
	  sum_ccs = 0;
	  for ( j=pnn[i]; j<pnn[i+1]; j++)
	    {
	      sum_nnz += nnz_contribute[j];
	      sum_nat += nat_contribute[j];
	      sum_ccs += ccs_contribute[j];
	    }
	  mono->nv[0][index_nnz][i]  = sum_nnz;
	  mono->nv[0][index_nat][i]  = sum_nat;
	  mono->nv[0][index_ccs][i]  = sum_ccs;
	  mono->nv[0][index_dofs][i] = node_dof0[i+1]-node_dof0[i];
	  mono->nv[0][index_nknd][i] = node_kind[i];
	  mono->nv[0][index_nscn][i] = psm[i+1]-psm[i];
	  mono->nv[0][index_nown][i] = sm[psm[i]];
	}
      
      /*
       * Create element variables that indicate:
       *	(1) element overassembly due to sharing of elements
       *	    across processors. index=1 means private, 2=shared
       *	    between 2 processors, 3...
       *
       *	(3) private elements get name of owning processor, shared
       *	    elements get a -1
       *
       * Note: With the new unambiguous element assignment needed for 
       *       discontinuous Galerkin methods, etc, we will assign the
       *       eown variable precisely to the correct processor.
       *
       * EXODUS II likes a truth table indicating which element variables
       * are active in which blocks. The element variables, too, are dispensed
       * on a per element block basis, so make provisions for doing so.
       *
       * elem_var_vals[ev_index][blk_index][elem_number]
       */

      nev       = mono->num_elem_vars;

      index_eoa = nev + 0;
      index_eown= nev + 1;

      len       = nev + 2;	

      if ( mono->num_ev_time_indeces > 1 )
	{
	  EH(-1, 
  "Reconcile new elemvar results with multiple time previous results.");
	}

      /*
       * Allocate or realloc, depending on whether the database already
       * had some element variables...
       */

      if ( nev > 0 )
	{
	  mono->elem_var_names = (char **) realloc(mono->elem_var_names, 
						   len * sizeof(char *));
	  mono->ev[0]          = (dbl **) realloc(mono->ev[0],
						  (len*neb) * sizeof(dbl *));
	}
      else
	{
	  mono->elem_var_names = (char **) smalloc(len * sizeof(char *));
	  mono->num_ev_time_indeces = 1;
	  mono->ev_time_indeces = (int *) smalloc(mono->num_ev_time_indeces*
						  sizeof(int));
	  mono->ev = (dbl ***) smalloc(mono->num_ev_time_indeces*
				       sizeof(dbl **));
	  for ( i=0; i<mono->num_ev_time_indeces; i++)
	    {
	      mono->ev[i] = (dbl **) smalloc(len*neb*sizeof(dbl *));
	      mono->ev_time_indeces[i] = i+1; /* Streetgang: "Fortran Rulz!" */
	    }

	  mono->state |= EXODB_STATE_ELIA; /* Keep accurate state info... */
	  mono->state |= EXODB_STATE_ELVA; /* Keep accurate state info... */
	}

      new_truth_table = (int *) smalloc(len*neb*sizeof(int));

      /*
       * Transcribe the old to the new.
       */

      if ( nev > 0 )
	{
	  for ( i=0; i<neb; i++)
	    {
	      for ( j=0; j<nev; j++)
		{
		  new_truth_table[i*len+j] = mono->elem_var_tab[i*nev+j];
		}
	    }
	}

      /*
       * Indicate the new element variables are to be active in every
       * element block.
       */

      for ( i=0; i<neb; i++)
	{
	  for ( j=nev; j<len; j++)
	    {
	      new_truth_table[i*len+j] = 1;
	    }
	}

      if ( nev > 0 )
	{
	  free(mono->elem_var_tab);
	}

      mono->elem_var_tab = new_truth_table;

      for ( i=nev; i<len; i++)
	{
	  mono->elem_var_names[i] = (char *) smalloc(MAX_STR_LENGTH*SZ_CHR);
	}

      for ( i=nev; i<len; i++)
	{
	  for ( j=0; j<neb; j++)
	    {
	      index = j * len + i;	      
	      mono->ev[0][index] = (dbl *) smalloc(mono->eb_num_elems[j]*
						       sizeof(dbl));
	    }
	}

      strcpy(mono->elem_var_names[index_eoa], "EOA");
      strcpy(mono->elem_var_names[index_eown], "EOWN");

      /*
       * Now fill in the arrays of element values.
       */

      element_procs = (int *) smalloc(MAX_ELEMENT_PROCS*SZ_INT);

      i=0;			/* overall element counter */
      for ( blk_index=0; blk_index<neb; blk_index++)
	{
	  for ( e=0; e<mono->eb_num_elems[blk_index]; e++, i++)
	    {

	      /*
	       * Look at every node this element touches and make a list
	       * of the processors that own nodes in this element.
	       * The length of that list is num_element_procs.
	       */
	      
	      num_element_procs = 0;
	      INIT_IVEC(element_procs, -1, MAX_ELEMENT_PROCS);	  

	      for ( l=ep[i]; l<ep[i+1]; l++)
		{
		  node = nl[l];
		  owner = sm[psm[node]];
		  BULL(owner, element_procs, num_element_procs);
		}

	      mono->ev[0][blk_index * len + index_eoa][e] = 
		(double) num_element_procs;

	      /*
	       * Ownership...
	       */
	      
	      mono->ev[0][blk_index*len+index_eown][e] = (double) 
		element_owner[i];

	    }
	}	  

      free(element_procs);

      mono->num_elem_vars = len;

      /*
       * Finally, let's pretend we have a time value. This is kind of
       * dirty and dangerous due to the potential for side effects on
       * memory management of the monolith and the polylithic time_vals[]...
       */

      if ( mono->num_times == 0 )
	{
	  mono->num_times    = 1;
	  mono->time_vals    = (dbl *) smalloc(mono->num_times*sizeof(dbl));
	  mono->time_vals[0] = 0.;
	}

      /*
       * Write out the modified monolith EXODUS II file.
       * REMEMBER that "mono" has been modified to include these
       * extra results variables!!!!
       */

      one_base(mono);
      wr_mesh_exo(mono, out_augplot_file_name, 0);
      wr_resetup_exo(mono, out_augplot_file_name, 0);
      wr_result_exo(mono, out_augplot_file_name, 0);
      zero_base(mono);

      free_exo_ev(mono);
    }

  free_exo(mono);
  free(mono);

  free(E);

  for ( i=0; i<MAX_NODE_KINDS; i++)
    {
      free(pnd[i]);
    }

  free(pnd);

  free(nnz_contribute);

  free(eqn_node_names);
  free(var_node_names);

  free(nat_contribute);
  free(ccs_contribute);

  free(pnn);

  free(node_dof0);
  free(node_kind);

  free(sm);
  free(psm);

  for ( i=0; i<neb; i++)
    {
      free(Lucky[i][0]);
      free(Lucky[i]);
      free(evd[i][0]);
      free(evd[i]);
    }
  free(Lucky);
  free(evd);

  free(dragon);
  free(assignment);

  for ( i=0; i<neb; i++)
    {
      for ( j=0; j<num_basic_eqnvars[i]; j++)
	{
	  free(mult[i][j]);
	}
      free(mult[i]);
    }
  free(mult);

  free(num_basic_eqnvars);

  free(send_proc_names);
  free(recv_proc_names);

  free(element_owner);
  free(element_owner_dist);
  free(element_bnd_stat);

  /*
   * Print Sam's heuristic...
   */

  numerator   = total_boundary_dofweight + total_external_dofweight;
  denominator = total_internal_dofweight + total_boundary_dofweight;

  if ( denominator != 0 )
    {
      fprintf(stdout, 
	      "Sam's dof weight figure of merit: (b+e)/(i+b) = %g\n",
	      ((double)(numerator))/((double)(denominator)));
    }

  fprintf(stdout, "-done.\n");

  if ( tmp == NULL || sr < 0 ) exit(2);

  return(0);
} /* end of main */


/* integer_compare() -- comparison function used by qsort. which is greater?
 *
 *
 * Created: 1998/02/26 07:27 MST pasacki@sandia.gov
 *
 * Revised:
 */

static int 
integer_compare(const void *arg1, 
                const void *arg2)
{
  int *a;
  int *b;
  
  a = (int *)arg1;
  b = (int *)arg2;

  /*
   * Primarily, sort according to the integer name of the owning processor.
   * Once that's done, then according to the global node number.
   */

  if ( *a < *b )
     {
        return(-1);
     }
  else if ( *a > *b )
     {
       return(1);
     }
  else
     {
       return(0);		/* should not normally happen! */
     }

}
