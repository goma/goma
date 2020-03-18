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
 
/* rd_exo.c --	open & read an EXODUS II finite element database. 
 *
 * rd_exo() --  input arguments
 *	x -- a pointer to a structure containing the fubdrmation from
 *	     an EXODUS IIv2 finite element database. A detailed description
 *	     of this structure may be found in "exo_struct.h"
 *
 *	fn -- points to a character string of the filename to be opened and
 *	      subsequently read
 *
 *	tasks -- there's a lot under the hood, dude. Things like read
 *		 initialization information, memory allocation, results
 *		 preliminaries, results data, for instance. These Boolean
 *		 flags can be added together to do multiple tasks on one
 *		 call, which is wanted at least as often as having only
 *		 selective tasks performed.
 *
 *	verbosity -- a integer that determines the level of output
 *		from rd_exo(). If zero, then no output is produced, larger
 *		values produce more detailed reporting of rd_exo()'s progress.
 *
 * Return values:
 *	integer -- A value of zero is returned upon normal completion of the
 *		   routine. A value of -1 is returned if abnormal conditions
 *		   were encountered.
 *
 * Notes:
 *
 *	1. rd_exo() allocates space as needed to accomodate the data from the
 *		    input file.
 *
 *	2. The structure gets filled in with as much data as it can find.
 *
 *	3. The emphasis is on reading mesh data and other preliminary
 *	   information. Little provision is made for reading results data
 *	   or history data from the file.
 *
 *      4. Implemented new variables for "state" of the database and for
 *	   "actions" to be performed.
 *
 * Modified: 1997/03/19 09:57 MST pasacki@sandia.gov
 *
 * Modified: 1998/12/18 09:57 MST pasacki@sandia.gov
 */

#define GOMA_RD_EXO_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef lint
#ifdef USE_RCSID
static char rcsid[] = "$Id: rd_exo.c,v 5.2 2008-05-02 19:07:57 hkmoffa Exp $";
#endif
#endif

#include "rd_exo.h"

#include "std.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_allo.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "mm_elem_block_structs.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "exodusII.h"
#include "mm_mp.h"
#include "rd_mesh.h"

struct Material_Properties;

struct Material_Properties;

/*
 * Variables used in several routines in this file.
 */

static const int sc  = sizeof(char);
static const int si  = sizeof(int);
static const int sd  = sizeof(dbl);
static const int spc = sizeof(char *);
static const int spi = sizeof(int *);
static const int spd = sizeof(dbl *);

Spfrtn sr;			/* sprintf() return type */
int 
rd_exo(Exo_DB *x,		/* def'd in exo_struct.h */
       const char *fn,
       const int verbosity,
       const int task)
{
  char err_msg[MAX_CHAR_ERR_MSG];
  int err;
  int i, j, ii, Iglobal = 0;
  int index;
  int k;
  int status=0;
  int len;

  char rc;			/* generic returned character (char)*/
  int ri;			/* generic returned integer (int)*/
  flt rf;			/* generic returned flt (float)*/
  /*  dbl rd;			 generic returned dbl (double)*/
  int nodal_var_index;
  int time_index;
  ELEM_BLK_STRUCT *eb_ptr = NULL;

  /*
   * Do not specify the ground state unless the action includes reading
   * initialization information...
   */

  if ( task & EXODB_ACTION_RD_INIT )
    {
      x->state = EXODB_STATE_GRND;
    }
  
  if ( verbosity > 0 )
    {
      fprintf(stderr, "rd_exo() begins.\n");
    }

  /*
   * Don't allocate memory for the file path unless this struct is in the
   * ground state. Thereafter, it won't be.
   */

  if ( x->state == EXODB_STATE_GRND )
    {
      len = strlen(fn) + 1;	/* Space for terminating null character. */
      if ( len > FILENAME_MAX_ACK )
	{
	  sr = sprintf(err_msg, "EXODUS II pathname \"%s\" too long (%d>%d).",
		       fn, len, FILENAME_MAX_ACK);
	  EH(-1, err_msg);
	}
      x->path = (char *) smalloc(len*sizeof(char));
    }
  strcpy(x->path, fn);

  /*
   * Setup some Boolean defaults. To begin, we have nothing.
   */

  x->node_map_exists       = FALSE;
  x->elem_map_exists       = FALSE;

  if( Linear_Solver == FRONT ) 
     {
       x->elem_order_map_exists = TRUE;
     }
   else
     x->elem_order_map_exists = FALSE;

  x->elem_var_tab_exists   = FALSE;

  /*
   * Christmastime 1998: Scrooge relents and doesn't force exodus to
   * carry information that parallel processing needs.
   */

  err = strlen(fn);

  if ( err > FILENAME_MAX_ACK )
    {
      sr = sprintf(err_msg, "EXODUS II pathname \"%s\" too long (%d>%d).",
		   fn, err, FILENAME_MAX_ACK);
      EH(-1, err_msg);
    }

  x->mode                  = EX_READ;
  x->version               = -4.98; /* initialize. ex_open() changes this. */

  /*
   * In GOMA, we intend to pass double values for our nodal variables, etc.
   */

  x->comp_wordsize = sizeof(dbl);

  /*
   * This is an existing file, so query it to find out what sized floating
   * point data exists in it.
   */

  x->io_wordsize   = 0;

  if ( verbosity > 1 )
    {
      fprintf(stderr, "ex_open() call...\n");

    }

  if ( verbosity > 2 )
    {
      fprintf(stderr, "\tx->path    = \"%s\"\n", x->path);
      fprintf(stderr, "\tx->mode    = %d\n", x->mode);
      fprintf(stderr, "\tx->comp_ws = %d\n", x->comp_wordsize);
      fprintf(stderr, "\tx->io_ws   = %d\n", x->io_wordsize);
      fprintf(stderr, "\tx->version = %g\n", x->version);
    }

  x->exoid = ex_open(x->path, x->mode, &(x->comp_wordsize), 
		     &(x->io_wordsize), &(x->version));
#ifdef PARALLEL

  /*
   *     Make sure the exodus files were all opened
   *    and if need be initiate a systematic shutdown 
   *       rd_exo -->  rd_mesh  -->  main
   *
   *  Error handling for serial is handled separately below
   */

  if( x->exoid >= 0 ) err = 0;
  parallel_err = err;
  if ( parallel_err )
    {
      fprintf(stderr,
              "\nProc %d: Exodus file read error, mesh files may not exist.\n", ProcID);
      return(-1);
    }
#endif

  if ( verbosity > 1 )
    {
      fprintf(stderr, "ex_open() rtn = %d\n", x->exoid);
    }
  if ( verbosity > 2 )
    {
      fprintf(stderr, "\tx->path    = \"%s\"\n", x->path);
      fprintf(stderr, "\tx->mode    = %d\n", x->mode);
      fprintf(stderr, "\tx->comp_ws = %d\n", x->comp_wordsize);
      fprintf(stderr, "\tx->io_ws   = %d\n", x->io_wordsize);
      fprintf(stderr, "\tx->version = %g\n", x->version);
    }

  if ( verbosity > 1 )
    {
      fprintf(stderr, "ex_open() rtn = %d\n", x->exoid);
    }

  EH(x->exoid, "ex_open");

  if ( task & EXODB_ACTION_RD_INIT )
    {

      /*
       * Check to see if we're trying to build on top of an existing
       * database? We really need to start from the ground zero state.
       */

      /*
	if ( x->state != EXODB_STATE_GRND )
	{
	EH(-1, "Attempt to rd init EXO db into previously used struct.");
	}
      */

      x->title = (char *) smalloc(MAX_SLENGTH*sc);

      for ( i=0; i<MAX_SLENGTH; i++)
	{
	  x->title[i] = '\0';
	}

      if ( verbosity > 1 )
	{
	  fprintf(stderr, "ex_get_init() call...\n");
	}
      status = ex_get_init(x->exoid, 
			   x->title, 
			   &x->num_dim, 
			   &x->num_nodes,
			   &x->num_elems, 
			   &x->num_elem_blocks, 
			   &x->num_node_sets,
			   &x->num_side_sets);
      EH(status, "ex_get_init");

      if ( verbosity > 0 )
	{
	  fprintf(stderr, "\tx->title           = \"%s\"\n", x->title);
	  fprintf(stderr, "\tx->num_nodes       = %d\n", x->num_nodes);
	  fprintf(stderr, "\tx->num_elems       = %d\n", x->num_elems);
	  fprintf(stderr, "\tx->num_elem_blocks = %d\n", x->num_elem_blocks);
	  fprintf(stderr, "\tx->num_node_sets   = %d\n", x->num_node_sets);
	  fprintf(stderr, "\tx->num_side_sets   = %d\n", x->num_side_sets);
	}

      status = ex_inquire(x->exoid, EX_INQ_QA, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_QA)");
      x->num_qa_rec = ri;

      if ( verbosity > 1 )
	{
	  fprintf(stderr, "\tx->num_qa_rec      = %d\n", x->num_qa_rec);
	}      

      if ( x->num_qa_rec > 0 )
	{

	  x->qa_record = (QA_Record *) smalloc(x->num_qa_rec*sizeof(QA_Record));

	  for ( i=0; i<x->num_qa_rec; i++)
	    {
	      for ( j=0; j<4; j++)
		{
		  x->qa_record[i][j] = (char *) smalloc(LEN_QA_RECORD*sc);

		  for ( k=0; k<LEN_QA_RECORD; k++)
		    {
		      x->qa_record[i][j][k] = '\0';
		    }
		}
	    }

	  status = ex_get_qa(x->exoid, x->qa_record);
	  EH(status, "ex_get_qa");
	}

      status = ex_inquire(x->exoid, EX_INQ_INFO,  &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_INFO)");
      x->num_info = ri;

      if ( x->num_info > 0 )
	{
	  /*
	    if ( x->num_info > MAX_INFO )
	    {
	    sr = sprintf(err_msg, 
	    "Number of info records in \"%s\" is %d > %d. Recompile w/ MAX_INFO larger.", 
	    fn, x->num_info, MAX_INFO);
	    EH(-1, err_msg);
	    }
	  */      

	  x->info = (INFO_Record *) smalloc(x->num_info*sizeof(INFO_Record));

	  for ( i=0; i<x->num_info; i++)
	    {
	      x->info[i] = (char *) smalloc((MAX_LINE_LENGTH+1)*sizeof(char));
	    }
	  status = ex_get_info(x->exoid, (x->info));
	  EH(status, "ex_get_info");
	}

      status = ex_inquire(x->exoid, EX_INQ_API_VERS, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_API_VERS)");
      x->api_version = rf;

      status = ex_inquire(x->exoid, EX_INQ_DB_VERS, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_DB_VERS)");
      x->db_version = rf;

      status = ex_inquire(x->exoid, EX_INQ_NS_NODE_LEN, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_NS_NODE_LEN)");
      x->ns_node_len = ri;

      status = ex_inquire(x->exoid, EX_INQ_NS_DF_LEN, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_NS_DF_LEN)");
      x->ns_distfact_len = ri;

      status = ex_inquire(x->exoid, EX_INQ_SS_ELEM_LEN, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_SS_ELEM_LEN)");
      x->ss_elem_len = ri;

      status = ex_inquire(x->exoid, EX_INQ_SS_DF_LEN, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_SS_DF_LEN)");
      x->ss_distfact_len = ri;

      status = ex_inquire(x->exoid, EX_INQ_SS_NODE_LEN, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_SS_NODE_LEN)");
      x->ss_node_len = ri;

      status = ex_inquire(x->exoid, EX_INQ_EB_PROP, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_EB_PROP)");
      x->eb_num_props = ri;

      status = ex_inquire(x->exoid, EX_INQ_NS_PROP, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_NS_PROP)");
      x->ns_num_props = ri;

      status = ex_inquire(x->exoid, EX_INQ_SS_PROP, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_SS_PROP)");
      x->ss_num_props = ri;

      status = ex_inquire(x->exoid, EX_INQ_TIME, &ri, &rf, &rc);
      EH(status, "ex_inquire()");
      x->num_times = ri;
      
      x->state |= EXODB_STATE_INIT;		/* Did it! */
    }

  if ( task & EXODB_ACTION_RD_MESH )
    {
      
      if ( ! ( x->state & EXODB_STATE_INIT ) )
	{
	  EH(-1, "Need to rd init info from EXO db before mesh.");
	}

      if ( x->state & EXODB_STATE_MESH )
	{
	  EH(-1, "Attempt to rd mesh EXO db into previously used struct.");
	}

      x->x_coord = NULL;
      x->y_coord = NULL;
      x->z_coord = NULL;

      if ( x->num_dim > 0 )
	{
	  x->x_coord = (dbl *) smalloc(x->num_nodes * sd);
	}

      if ( x->num_dim > 1 )
	{
	  x->y_coord = (dbl *) smalloc(x->num_nodes * sd);
	}
    
      if ( x->num_dim > 2 )
	{
	  x->z_coord = (dbl *) smalloc(x->num_nodes * sd);
	}

      status = ex_get_coord(x->exoid, x->x_coord, x->y_coord, x->z_coord);
      EH(status, "ex_get_coord");

#if DEBUG_LEVEL > 1
      for ( i=0; i<x->num_nodes; i++)
	{
	  log_msg("Node position x[%d] = %g", i, x->x_coord[i]);
	}
#endif

      if ( x->num_dim > 0 )
	{
	  x->coord_names = (char **)smalloc( x->num_dim * spc);
	  for ( i=0; i<x->num_dim; i++ )
	    {
	      x->coord_names[i] = (char *) calloc(MAX_STR_LENGTH+1, sc);
	    }
	  status = ex_get_coord_names(x->exoid, x->coord_names);
	  EH(status, "ex_get_coord_names");
	}

      if ( x->num_nodes > 0 )
	{
	  if ( x->node_map_exists )
	    {
	      x->node_map = (int *) smalloc(x->num_nodes * si);
	      status = ex_get_id_map(x->exoid, EX_NODE_MAP, x->node_map);

	      EH(status, "ex_get_node_num_map");
	    }
	}

      if ( x->num_elems > 0 )
	{
	  if ( x->elem_map_exists )
	    {
	      x->elem_map = (int *) smalloc(x->num_elems * si);
	      status      = ex_get_id_map(x->exoid, EX_ELEM_MAP, x->elem_map);
	      EH(status, "ex_get_elem_num_map");
	    }
	  if ( x->elem_order_map_exists )
	    {
	      x->elem_order_map = (int *) smalloc(x->num_elems * si);
	      status            = ex_get_map(x->exoid, x->elem_order_map);
	      EH(status, "ex_get_map");
	    }

	  x->elem_eb = (int *) smalloc(x->num_elems * si);
	  for (i = 0; i < x->num_elems; i++ ) x->elem_eb[i] = -1;
	}

      /*
       * ELEMENT BLOCKS...
       */

      /*
       * Allocate storage for element block structures and then zero
       * the storage space.
       */
      Element_Blocks = (ELEM_BLK_STRUCT *)
	  alloc_struct_1(ELEM_BLK_STRUCT, x->num_elem_blocks);

      if ( x->num_elem_blocks > 0 )
	{
	  x->eb_id        = (int *) smalloc(x->num_elem_blocks * si);
	  x->eb_elem_type = (char **) smalloc(x->num_elem_blocks * spc);
	  x->eb_num_elems = (int *) smalloc(x->num_elem_blocks * si);
	  x->eb_num_nodes_per_elem = (int *) smalloc(x->num_elem_blocks * si);
	  x->eb_num_attr  = (int *) smalloc(x->num_elem_blocks * si);
	  x->eb_conn = (int **) smalloc(x->num_elem_blocks * spi);
	  x->eb_attr = (dbl **) smalloc(x->num_elem_blocks * spd);

	  x->eb_ptr    = (int *) smalloc((x->num_elem_blocks+1) * si);
	  x->eb_ptr[0] = 0;

	  for ( i=0; i<x->num_elem_blocks; i++)
	    {
	      x->eb_elem_type[i] = (char *) smalloc(MAX_STR_LENGTH * sc);
	    }
          /*
	   *  Read the element block ID's from the exodus file
	   */
	  status = ex_get_ids(x->exoid, EX_ELEM_BLOCK, x->eb_id);
	  EH(status, "ex_get_ids elem_blk_ids");

	  /*
	   *  For each element block, specified by the id
	   *  Get the element type, number of elements in the block
	   *  the nodes per element, and the number of attributes
	   *  for a given element
	   */
	  for ( i=0; i<x->num_elem_blocks; i++)
	    {
	      status = ex_get_block(x->exoid,
				    EX_ELEM_BLOCK,
				    x->eb_id[i],
				    x->eb_elem_type[i],
				    &x->eb_num_elems[i],
				    &x->eb_num_nodes_per_elem[i],
				    0,
				    0,
				    &x->eb_num_attr[i]);
	      EH(status, "ex_get_block elem");

	      /*
	       *  Fill in the information in the current
	       *  element block structure
	       */
              eb_ptr = Element_Blocks + i;
	      eb_ptr->Elem_Blk_Num = i;
	      eb_ptr->Elem_Blk_Id = x->eb_id[i];
	      eb_ptr->Elem_Type = get_type(x->eb_elem_type[i], 
					   x->eb_num_nodes_per_elem[i],
					   x->eb_num_attr[i]);
	      eb_ptr->Num_Nodes_Per_Elem = x->eb_num_nodes_per_elem[i];
	      eb_ptr->Num_Attr_Per_Elem = x->eb_num_attr[i];
              int mindex = map_mat_index(eb_ptr->Elem_Blk_Id);
              if (mindex < 0) {
                eb_ptr->MatlProp_ptr = NULL;
              } else {
	        eb_ptr->MatlProp_ptr = mp_glob[mindex];
              }
	      eb_ptr->ElemStorage = NULL;
	      eb_ptr->Num_Elems_In_Block = x->eb_num_elems[i];
              eb_ptr->IP_total = elem_info(NQUAD, eb_ptr->Elem_Type);
	      
              /*
	       * Go on to read the information about the element block
	       * from the exodus file
	       */
	      if ( (x->eb_num_elems[i] * x->eb_num_nodes_per_elem[i]) > 0 )
		{
                  /*
		   *  Allocate storage for the element connectivity array.
		   *  Then, read the array from exodus data base
		   */
		  x->eb_conn[i] = 
		    (int *) smalloc( x->eb_num_elems[i] * 
				     x->eb_num_nodes_per_elem[i] * si);

		  status = ex_get_conn(x->exoid, EX_ELEM_BLOCK,
				       x->eb_id[i],
				       x->eb_conn[i], 0, 0);
		  EH(status, "ex_get_elem_conn");

		  /* Build the element - element block index map */
		  for (ii = 0; ii < x->eb_num_elems[i]; ii++ ) {	
		    x->elem_eb[Iglobal++] = i;
		  }
		}

	      if ( (x->eb_num_elems[i]*x->eb_num_attr[i]) > 0 )
		{
		  x->eb_attr[i] = (dbl *) 
		    smalloc(x->eb_num_elems[i] * x->eb_num_attr[i] * sd);
		  status = ex_get_attr(x->exoid, EX_ELEM_BLOCK, x->eb_id[i],
				       x->eb_attr[i]);
		  EH(status, "ex_get_attr elem");
		}
	  
	      x->eb_ptr[i+1] = x->eb_ptr[i] + x->eb_num_elems[i];
	    }
	}

      /*
       * NODE SETS...
       */

      if ( x->num_node_sets > 0 )
	{
	  x->ns_id             = (int *) smalloc(x->num_node_sets * si);
	  x->ns_num_nodes      = (int *) smalloc(x->num_node_sets * si);
	  x->ns_num_distfacts  = (int *) smalloc(x->num_node_sets * si);
	  x->ns_node_index     = (int *) smalloc(x->num_node_sets * si);
	  x->ns_distfact_index = (int *) smalloc(x->num_node_sets * si);


	  if ( x->ns_node_len > 0 )
	    {
	      x->ns_node_list = (int *) smalloc(x->ns_node_len * si);
	    }

	  if ( x->ns_distfact_len > 0 )
	    {
	      x->ns_distfact_list = (dbl *) smalloc(x->ns_distfact_len * sd);

	      ex_set_specs ns_specs;
	      ns_specs.sets_ids = x->ns_id;
	      ns_specs.num_entries_per_set = x->ns_num_nodes;
	      ns_specs.num_dist_per_set    = x->ns_num_distfacts;
	      ns_specs.sets_entry_index    = x->ns_node_index;
	      ns_specs.sets_dist_index     = x->ns_distfact_index;
	      ns_specs.sets_entry_list     = x->ns_node_list;
	      ns_specs.sets_extra_list     = NULL;
	      ns_specs.sets_dist_fact      = x->ns_distfact_list;

	      status = ex_get_concat_sets(x->exoid, EX_NODE_SET,
					  &ns_specs);
	      EH(status, "ex_get_concat_node_sets");
	    }
	}


      /*
       * SIDE SETS...
       */

      if ( x->num_side_sets > 0 ) 
	{
	  x->ss_id             = (int *) smalloc(x->num_side_sets * si);
	  x->ss_num_sides      = (int *) smalloc(x->num_side_sets * si);
	  x->ss_num_distfacts  = (int *) smalloc(x->num_side_sets * si);
	  x->ss_elem_index     = (int *) smalloc(x->num_side_sets * si);
	  x->ss_distfact_index = (int *) smalloc(x->num_side_sets * si);

	  x->ss_node_cnt_list  = (int **) smalloc(x->num_side_sets * spi);
	  x->ss_node_list      = (int **) smalloc(x->num_side_sets * spi);

	  x->ss_node_side_index  = (int **) smalloc(x->num_side_sets * spi);

	  if ( x->ss_elem_len > 0 )
	    {
	      x->ss_elem_list = (int *) smalloc(x->ss_elem_len * si);
	      x->ss_side_list = (int *) smalloc(x->ss_elem_len * si);
	    }

	  if ( x->ss_distfact_len > 0 )
	    {
	      x->ss_distfact_list = (dbl *) smalloc(x->ss_distfact_len * sd);
	    }	    

	  ex_set_specs ss_specs;
	  ss_specs.sets_ids = x->ss_id;
	  ss_specs.num_entries_per_set = x->ss_num_sides;
	  ss_specs.num_dist_per_set    = x->ss_num_distfacts;
	  ss_specs.sets_entry_index    = x->ss_elem_index;
	  ss_specs.sets_dist_index     = x->ss_distfact_index;
	  ss_specs.sets_entry_list     = x->ss_elem_list;
	  ss_specs.sets_extra_list     = x->ss_side_list;
	  ss_specs.sets_dist_fact      = x->ss_distfact_list;

	  status = ex_get_concat_sets(x->exoid, EX_SIDE_SET,
				      &ss_specs);

	  EH(status, "ex_get_concat_side_sets");

	  /*
	   * This information turns out to be useful in constructing more
	   * rapid indeces into the distribution factor array.
	   */

	  x->ss_node_list_exists = TRUE; /* does now! */
      

	  for ( i=0; i<x->num_side_sets; i++)
	    {
	      x->ss_node_cnt_list[i] = (int *) smalloc(x->ss_num_sides[i] * si);

	      x->ss_node_list[i] = (int *) smalloc(x->ss_num_distfacts[i] * si);

	      status = ex_get_side_set_node_list(x->exoid,
						 x->ss_id[i],
						 x->ss_node_cnt_list[i],
						 x->ss_node_list[i]);
#ifdef DEBUG
	      fprintf(stderr,"P_%d, SSID=%d has %d dfs/nds on -> %d <- sides.\n", 
                      ProcID,
		      x->ss_id[i],
		      x->ss_num_distfacts[i], x->ss_num_sides[i]);

	      fprintf(stderr, "Address x->ss_num_sides[%d] is %x\n", i, 
		      &(x->ss_num_sides[i]));
	  
	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "nodes for elem %d side %d = %d\n", 
			  x->ss_elem_list[x->ss_elem_index[i]+j],
			  x->ss_side_list[x->ss_elem_index[i]+j],
			  x->ss_node_cnt_list[i][j]);
		}

	      /*
	       * Did we read the distribution factors OK?
	       */

	      for ( j=0; j<x->ss_num_distfacts[i]; j++)
		{
		  fprintf(stderr, "(%d,%d) distfact = %g\n", i, j, 
			  x->ss_distfact_list[x->ss_distfact_index[i]+j]);
		}

#endif
	      /*
	       * Set up quick pointers for nodes on each given side of a sideset
	       * that can be used later to find exactly where to go in the big 
	       * distribution factor list...
	       */

	      x->ss_node_side_index[i] = (int *)smalloc((x->ss_num_sides[i]+1)*si);

	      x->ss_node_side_index[i][0] = 0;

	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  x->ss_node_side_index[i][j+1] = ( x->ss_node_side_index[i][j] +
						    x->ss_node_cnt_list[i][j] );
		}

#ifdef DEBUG
	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "P_%d SS[%d]=%d, nodes for elem %d, side %d:", i, 
                          ProcID,
			  x->ss_id[i], x->ss_elem_list[x->ss_elem_index[i]+j],
			  x->ss_side_list[x->ss_elem_index[i]+j]);


		  for ( k=x->ss_node_side_index[i][j]; 
			k<x->ss_node_side_index[i][j+1]; k++)
		    {
		      fprintf(stderr, " %d", x->ss_node_list[i][k]);
		    }
		  fprintf(stderr, "\n");
		}
#endif
	    }
	}


      /*
       * PROPERTIES...
       */

      /*
       * Node sets...
       */

      if ( x->ns_num_props > 0 ) 
	{
	  x->ns_prop_name = (char **) smalloc(x->ns_num_props * spc);
	  for ( i=0; i<x->ns_num_props; i++)
	    {
	      x->ns_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH * sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_NODE_SET, x->ns_prop_name);
	  EH(status, "ex_get_prop_names(EX_NODE_SET)");

	  x->ns_prop = (int **) smalloc(x->ns_num_props * spi);
	  for ( i=0; i<x->ns_num_props; i++)
	    {
	      x->ns_prop[i] = (int *)smalloc(x->num_node_sets * si);
	    }

	  for ( i=0; i<x->ns_num_props; i++)
	    {
	      status = ex_get_prop_array(x->exoid, EX_NODE_SET, 
					 x->ns_prop_name[i],
					 x->ns_prop[i]);
	      EH(status, "ex_get_prop_array(EX_NODE_SET)");
	    }
	}
      
      /*
       * Side sets...
       */

      if ( x->ss_num_props > 0 ) 
	{
	  x->ss_prop_name = (char **) smalloc(x->ss_num_props * spc);
	  for ( i=0; i<x->ss_num_props; i++)
	    {
	      x->ss_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH * sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_SIDE_SET, x->ss_prop_name);
	  EH(status, "ex_get_prop_names(EX_SIDE_SET)");

	  x->ss_prop = (int **) smalloc(x->ss_num_props * spi);
	  for ( i=0; i<x->ss_num_props; i++)
	    {
	      x->ss_prop[i] = (int *)smalloc(x->num_side_sets * si);
	    }

	  for ( i=0; i<x->ss_num_props; i++)
	    {
	      status = ex_get_prop_array(x->exoid, EX_SIDE_SET, 
					 x->ss_prop_name[i],
					 x->ss_prop[i]);
	      EH(status, "ex_get_prop_array(EX_SIDE_SET)");
	    }
	}
      
      /*
       * Element blocks...
       */

      if ( x->eb_num_props > 0 ) 
	{
	  x->eb_prop_name = (char **) smalloc(x->eb_num_props * spc);
	  for ( i=0; i<x->eb_num_props; i++)
	    {
	      x->eb_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH * sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_ELEM_BLOCK, x->eb_prop_name);
	  EH(status, "ex_get_prop_names(EX_ELEM_BLOCK)");

	  x->eb_prop = (int **) smalloc(x->eb_num_props * spi);
	  for ( i=0; i<x->eb_num_props; i++)
	    {
	      x->eb_prop[i] = (int *)smalloc(x->num_elem_blocks * si);
	    }

	  for ( i=0; i<x->eb_num_props; i++)
	    {
	      status = ex_get_prop_array(x->exoid, EX_ELEM_BLOCK, 
					 x->eb_prop_name[i],
					 x->eb_prop[i]);
	      EH(status, "ex_get_prop_array(EX_ELEM_BLOCK)");
	    }
	}
      
      x->state |= EXODB_STATE_MESH; /* Did it! */
    }

  /*
   * Read in result preliminaries. 
   *
   * In particular, read number of nodal, element, global variables 
   * as well as their names.
   *
   * Also, read the number of timeplanes and the element variable 
   * polygraph test results, Monica - Bill?
   */

  if ( task & EXODB_ACTION_RD_RES0 )
    {

      /*
       * Check for sanity - better not be doing this unless a fair
       * amount of information is already known about this database.
       */

      if ( ! ( x->state & EXODB_STATE_INIT ) )
	{
	  EH(-1, "Need to rd init info from EXO db before result metadata.");
	}

      /*
       * Hmmm... is it really just the element block IDs for element
       * variables that are required? Perhaps for nodal variables it suffices
       * to know merely the number of nodes...
       */

      if ( ! ( x->state & EXODB_STATE_MESH ) )
	{
	  EH(-1, "Need to rd mesh info from EXO db before result metadata.");
	}

      if ( x->state & EXODB_STATE_RES0 )
	{
	  EH(-1, "Attempt to rd mesh EXO db into previously used struct.");
	}

      /*
       * Initialize, then read what they really are.
       */

      x->num_glob_vars = 0;
      x->num_elem_vars = 0;
      x->num_node_vars = 0;

      /*
       * Sometimes EXODUS II v2.** interacts badly with netCDF 3.*
       * in the sense that "zero" answers are often equivalent to 
       * "undefined" answers. I think EXODUS II v3.00 fixes some
       * of these problems...
       */

      status = ex_get_variable_param(x->exoid, EX_GLOBAL, &x->num_glob_vars);
      EH(status, "ex_get_variable_param global");

      status = ex_get_variable_param(x->exoid, EX_ELEM_BLOCK, &x->num_elem_vars);
      EH(status, "ex_get_variable_param elem");

      status = ex_get_variable_param(x->exoid, EX_NODAL, &x->num_node_vars);
      EH(status, "ex_get_variable_param nodal");

      /*
       * Get the names of the results variables: global, element and nodal.
       */

      if ( x->num_glob_vars > 0 )
	{
	  x->glob_var_names = (char **) smalloc(x->num_glob_vars * spc);
	  for ( i=0; i<x->num_glob_vars; i++)
	    {
	      x->glob_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1) * sc);
	    }
	  status = ex_get_variable_names(x->exoid, EX_GLOBAL, x->num_glob_vars,
					 x->glob_var_names);
	  EH(status, "ex_get_variable_names global");
	}

      if ( x->num_elem_vars > 0 )
	{
	  x->elem_var_names = (char **) smalloc(x->num_elem_vars * spc);
	  for ( i=0; i<x->num_elem_vars; i++)
	    {
	      x->elem_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1) * sc);
	    }
	  status = ex_get_variable_names(x->exoid, EX_ELEM_BLOCK, x->num_elem_vars,
					 x->elem_var_names);
	  EH(status, "ex_get_variable_names elem_block");
	}

      if ( x->num_node_vars > 0 )
	{
	  x->node_var_names = (char **) smalloc(x->num_node_vars * spc);
	  for ( i=0; i<x->num_node_vars; i++)
	    {
	      x->node_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1) * sc);
	    }
	  status = ex_get_variable_names(x->exoid, EX_NODAL, x->num_node_vars,
				    x->node_var_names);
	  EH(status, "ex_get_variable_names nodal");
	}

      /*
       * Get time values...
       */

      if ( x->num_times > 0 )
	{
	  x->time_vals = (dbl *) smalloc(x->num_times * sd);
	  status = ex_get_all_times(x->exoid, x->time_vals);
	  EH(status, "ex_get_all_times");
	}

      /*
       * Any element variable truth table...
       */

      if ( x->num_elem_vars > 0 )
	{
	  x->elem_var_tab = (int *) smalloc(x->num_elem_vars * 
					    x->num_elem_blocks * si);

	  status = ex_get_truth_table(x->exoid, EX_ELEM_BLOCK, x->num_elem_blocks,
				      x->num_elem_vars, x->elem_var_tab);
	  EH(status, "ex_get_truth_table elem");
	}
      x->state |= EXODB_STATE_RES0;
    }


#if 0				/* old stuff relegated to scrap heap 12/21/98*/
  /*
   * Read nodal variable values at last time step...
   */

  if ( x->num_node_vars > 0 )
    {
      x->node_var_vals = (dbl **) smalloc(x->num_node_vars * spd);
      for ( i=0; i<x->num_node_vars; i++)
	{
	  x->node_var_vals[i] = (dbl *) smalloc(x->num_nodes * sd);
	}
	
      for ( i=0; i<x->num_node_vars; i++)      
	{
	  status = ex_get_nodal_var(x->exoid, 
				    x->num_times, 
				    i+1, 
				    x->num_nodes,
				    x->node_var_vals[i]);
	  EH(status, "ex_get_nodal_var");
	}
    }

#endif

  /*
   * New stuff from brkfix suite of {rd,wr}_exo...
   */
  if ( task & EXODB_ACTION_RD_RESN )
    {
      if ( !(x->state & EXODB_STATE_RES0) )
	{
	  EH(-1, "Need to rd result metadata from EXO db before result data.");
	}

      /*
       * Should be OK to read in another clump of data and to re-use
       * the allocated space if we want to save space by reading clumps
       * of data incrementally...
       *
       * The only hitch would be if the user changed x->num_nv_indeces
       * in the interim. Good opportunity for public and private data, eh?
       */
      /*
       * The idea is this. The variable x->num_node_vars will indicated what
       * is in the database file, while x->num_nv_indeces will be used to
       * indicate memory that has been allocated for the memory model.
       * Thus, if x->num_nv_indeces > 0, then we expect the arrays have
       * been allocated and places made so the request for data can be
       * satisfied.
       *
       * It is up to the user to say what, if any, nodal results data they
       * want.
       */

#if 0
      /*
       * It must be "OK" for use to retrieve one nodal variable at a time
       * from the database, so don't make this an error. It's more like
       * a weakly informative message...
       */

      if ( x->num_node_vars > x->num_nv_indeces )
	{
	  sr = sprintf(err_msg, "Database has %d node vars, you want %d.",
		       x->num_node_vars, x->num_nv_indeces);
	  EH(-1, err_msg);
	}
#endif

      /*
       * Don't do anything unless there's really something to do...
       */

      if ( x->num_nv_time_indeces > 0   &&
	   x->num_nv_indeces      > 0   &&
	   x->num_nodes           > 0	&&
	   x->num_node_vars       > 0 )
	{

	  /*
	   * Nodal variable results -- allocate space if needed.
	   */

	  if ( ! ( x->state & EXODB_STATE_NDVA ) )
	    {
	      alloc_exo_nv(x, 1, 1);
	    }

	  /*
	   * Nodal results -- get values.
	   */
	  
	  for ( i=0; i<x->num_nv_time_indeces; i++)
	    {
	      time_index = x->nv_time_indeces[i];
	      for ( j=0; j<x->num_nv_indeces; j++)
		{
		  nodal_var_index = x->nv_indeces[j];
		  
		  status = ex_get_var(x->exoid, time_index, EX_NODAL,
				      nodal_var_index, 1, x->num_nodes,
				      x->nv[i][j]);
		  EH(status, "ex_get_nodal_var");
		}
	    }
	}

      /*
       * Indicate that this memory representation of the database
       * has some meaty nodal results variables in it.
       */

      x->state |= EXODB_STATE_NDVR;
    }


  /*
   * Global variable results.
   */

  if ( task & EXODB_ACTION_RD_RESG )
    {

      /*
       * Should do some error checking here to make sure things are
       * setup and reasonable.
       */

      if ( x->num_gv_time_indeces > 0 &&
	   x->num_glob_vars > 0 )
	{

	  if ( ! ( x->state & EXODB_STATE_GBVA ) )
	    {
	      alloc_exo_gv(x, 1);
	    }

	  for ( i=0; i<x->num_gv_time_indeces; i++)
	    {
	      time_index = x->gv_time_indeces[i];
	      status = ex_get_var(x->exoid, time_index, EX_GLOBAL, 1, 1,
				  x->num_glob_vars,
				  x->gv[i]);
	      EH(status, "ex_get_glob_vars");
	    }

	  /*
	   * Indicate that this memory representation of the database
	   * has some meaty global results variables in it.
	   */

	  x->state |= EXODB_STATE_GBVR;
	}
    }

  /*
   * Element variable results.
   */

  if ( task & EXODB_ACTION_RD_RESE )
    {
   
      if ( x->num_ev_time_indeces > 0 &&
	   x->num_elem_vars       > 0 &&
	   x->num_elems	          > 0 )
	{

	  /*
	   * Element variable results -- allocate space.
	   */

	  if ( ! ( x->state & EXODB_STATE_ELVA ) )
	    {
	      alloc_exo_ev(x, 1);
	    }


	  /*
	   * Element results -- get values.
	   */
	  
	  for ( i=0; i<x->num_ev_time_indeces; i++)
	    {
	      time_index = x->ev_time_indeces[i];

	      for ( j=0; j<x->num_elem_blocks; j++)
		{
		  for ( k=0; k<x->num_elem_vars; k++)
		    {
		      index = j * x->num_elem_vars + k;

		      if ( x->elem_var_tab[index] != 0 )
			{
			  status = ex_get_var(x->exoid, time_index, EX_ELEM_BLOCK,
					      k+1,
					      x->eb_id[j],
					      x->eb_num_elems[j],
					      x->ev[i][index]);
			  if ( status < 0 )
			    {
			      sr = sprintf(err_msg,  "ex_get_elem_var() bad rtn: time %d, elemvar %d, EB ID %d",
					   time_index, k+1, x->eb_id[j]);
			      EH(-1, err_msg);
			    }
			}
		    }
		}
	    }

	  /*
	   * Leave a flag for anyone to see that this memory model of
	   * data really has some interesting element variable results
	   * in it now.
	   */

	  x->state |= EXODB_STATE_ELVR;
	}
    }




  status = ex_close(x->exoid);
  EH(status, "ex_close");

  return(status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  free_exo: Free memory previously malloced for storage of information
 *            from an exodus file
 *
 */
int
free_exo(Exo_DB *x)		/* pointer to EXODUS II FE db structure */
{
  int i;
  int j;

  free(x->path);

  free(x->title);

  if ( x->num_qa_rec > 0 )
    {
      for ( i=0; i<x->num_qa_rec; i++)
	{
	  for ( j=0; j<4; j++)
	    {
	      free(x->qa_record[i][j]);
	    }
	}
      free(x->qa_record);
    }

  if ( x->num_info > 0 )
    {
      for ( i=0; i<x->num_info; i++)
	{
	  free(x->info[i]);
	}
      free(x->info);
    }

  if ( x->elem_map_exists )
    {
      free(x->elem_map);
    }

  if ( x->node_map_exists )
    {
      free(x->node_map);
    }

  safer_free((void **) &(x->elem_order_map));

  if ( x->num_dim > 0 )
    {
      free(x->x_coord);
    }

  if ( x->num_dim > 1 )
    {
      free(x->y_coord);
    }
    
  if ( x->num_dim > 2 )
    {
      free(x->z_coord);
    }

  if ( x->num_dim > 0 )
    {
      for ( i=0; i<x->num_dim; i++ )
	{
	  free(x->coord_names[i]);
	}
      free(x->coord_names);
    }

  if ( x->num_elems > 0 ) {
    free( x->elem_eb );
  }

  if ( x->num_elem_blocks > 0 )
    {
      for ( i=0; i<x->num_elem_blocks; i++)
	{
	  free(x->eb_elem_type[i]);
	  free(x->eb_conn[i]);
	  if ( (x->eb_num_elems[i]*x->eb_num_attr[i]) > 0 )
	    {
	      free(x->eb_attr[i]);
	    }
	}

      free(x->eb_ptr);

      free(x->eb_id);
      free(x->eb_num_elems);
      free(x->eb_num_nodes_per_elem);
      free(x->eb_num_attr);

      free(x->eb_elem_type);
      free(x->eb_conn);
      free(x->eb_attr);
    }

  if ( x->num_node_sets > 0 )
    {
      free(x->ns_id);
      free(x->ns_num_nodes);
      free(x->ns_num_distfacts);
      free(x->ns_node_index);
      free(x->ns_distfact_index);

      if ( x->ns_node_len > 0 )
	{
	  free(x->ns_node_list);
	}

      if ( x->ns_distfact_len > 0 )
	{
	  free(x->ns_distfact_list);
	}
    }

  if ( x->num_side_sets > 0 ) 
    {
      
      if ( x->ss_node_list_exists )
	{
	  for ( i=0; i<x->num_side_sets; i++)
	    {
	      free(x->ss_node_cnt_list[i]);
	      free(x->ss_node_list[i]);
	      free(x->ss_node_side_index[i]);
	    }
	  free(x->ss_node_cnt_list);
	  free(x->ss_node_list);
	  free(x->ss_node_side_index);
	}

      free(x->ss_id);
      free(x->ss_num_sides);
      free(x->ss_num_distfacts);
      free(x->ss_elem_index);
      free(x->ss_distfact_index);

      if ( x->ss_elem_len > 0 )
	{
	  free(x->ss_elem_list);
	  free(x->ss_side_list);
	}

      if ( x->ss_distfact_len > 0 )
	{
	  free(x->ss_distfact_list);
	}
    }

  /*
   * PROPERTIES...
   */

  /*
   * Node sets...
   */

  if ( x->ns_num_props > 0 &&
       x->ns_prop_name != NULL && 
       x->ns_prop != NULL ) 
    {
      for ( i=0; i<x->ns_num_props; i++)
	{
	  free(x->ns_prop_name[i]);
	  free(x->ns_prop[i]);
	}
      free(x->ns_prop_name);
      free(x->ns_prop);
    }
      
  /*
   * Side sets...
   */

  if ( x->ss_num_props > 0 &&
       x->ss_prop_name != NULL &&
       x->ss_prop != NULL ) 
    {
      for ( i=0; i<x->ss_num_props; i++)
	{
	  free(x->ss_prop_name[i]);
	  free(x->ss_prop[i]);
	}
      free(x->ss_prop_name);
      free(x->ss_prop);
    }
      
  /*
   * Element blocks...
   */

  if ( x->eb_num_props > 0 &&
       x->eb_prop_name != NULL && 
       x->eb_prop != NULL ) 
    {
      for ( i=0; i<x->eb_num_props; i++)
	{
	  free(x->eb_prop_name[i]);
	  free(x->eb_prop[i]);
	}
      free(x->eb_prop_name);
      free(x->eb_prop);
    }
      
  /*
   * Results variables...
   */

  if ( x->num_glob_vars > 0 )
    {
      for ( i=0; i<x->num_glob_vars; i++)
	{
	  free(x->glob_var_names[i]);
	}
      free(x->glob_var_names);
    }

  if ( x->num_elem_vars > 0 )
    {
      for ( i=0; i<x->num_elem_vars; i++)
	{
	  /* Sanity check - we may have elem vars but when they are not
	     read in, these arrays were never malloc'd. The names are then
	     constructed in wr_result_prelim_exo */
	  free(x->elem_var_names[i]);
	}
      free(x->elem_var_names);
    }

  if ( x->num_node_vars > 0 )
    {
      for ( i=0; i<x->num_node_vars; i++)
	{
	  free(x->node_var_names[i]);
	}
      free(x->node_var_names);
    }

  if ( x->num_times > 0 )
    {
      free(x->time_vals);
    }

  if ( x->num_elem_vars > 0 )
    {
      free(x->elem_var_tab);
    }

#if 0
  /*
   * This is now done in free_exo_nv().
   */
  if ( x->num_node_vars > 0 )
    {
      for ( i=0; i<x->num_node_vars; i++)
	{
	  free(x->node_var_vals[i]);
	}
      free(x->node_var_vals);
    }
#endif

  /*
   * Free up any auxiliary connectivity arrays that might have been
   * created. Boolean variables tell us if that happened.
   */

  if ( x->elem_node_conn_exists )
    {
      free(x->elem_ptr);
      free(x->node_list);
    }

  if ( x->node_elem_conn_exists )
    {
      free(x->node_elem_pntr);
      free(x->node_elem_list);
    }

  if (x->elem_elem_conn_exists) {
    safer_free((void **) &(x->elem_elem_pntr));
    safer_free((void **) &(x->elem_elem_list));
    safer_free((void **) &(x->elem_elem_xadj));
    safer_free((void **) &(x->elem_elem_adjncy));
    safer_free((void **) &(x->elem_elem_twst));
    safer_free((void **) &(x->elem_elem_face));
  }

  if ( x->node_node_conn_exists )
    {
      free(x->node_node_pntr);
      free(x->node_node_list);
      free(x->centroid_list);
    }

  if ( x->elem_var_tab_exists ) {
    free(x->truth_table_existance_key);
    free(x->elem_var_tab);
  }

  return(0);
}	 

/*
 * zero_base() -- push down the element names and node names by one in an
 *                EXODUS II data base. This makes C language zero-based
 *                arrays work more naturally. Incoming data structure is
 *		  1-based FORTRAN like node names, etc. while the changed
 *		  data is 0-based.
 *
 * Created: 1997/08/26 09:17 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
zero_base(Exo_DB *E)
{
  int eb;
  int i,j;
  int l;
  int length_conn;

  /*
   * 1. Node numbers are named in the connectivity lists for each 
   *    element block. Decrement them.
   */

  for ( eb=0; eb<E->num_elem_blocks; eb++)
    {
      length_conn =  E->eb_num_elems[eb] * E->eb_num_nodes_per_elem[eb];
      for ( l=0; l<length_conn; l++)
	{
	  (E->eb_conn[eb][l])--;
	}
    }

  /*
   * 2. Node numbers are named in the node lists for each node set.
   *    Decrement them.
   */

  for ( l=0; l<E->ns_node_len; l++)
    {
      (E->ns_node_list[l])--;
    }

  /*
   * 3. Element numbers are named in the element list for the side sets.
   *    Decrement these.
   */

  for ( l=0; l<E->ss_elem_len; l++)
    {
      (E->ss_elem_list[l])--;
    }

  /*
   * 4. Node numbers in the special node list for the sidesets.
   *    Decrement these (if they exist).
   */

  if ( E->ss_node_list_exists )
    {
      for ( i=0; i<E->num_side_sets; i++)
	{
	  
	  /*
	   * Instead of doing a side at a time, do the whole list whose
	   * length has been counted up and stored in the last element of this
	   * index.
	   */
	  
	  for ( j=0; j<E->ss_node_side_index[i][ E->ss_num_sides[i] ]; j++)
	    {
	      (E->ss_node_list[i][j])--;
	    }
	}
    }

  /*
   * 5. If a node map exists, the names of the nodes might be too high.
   *    Decrement them. In some cases for very general node maps this
   *    step might be superfluous and even undesirable. Here, the idea is
   *    that the node map probably corresponds to an inherent contiguous
   *    node numbering scheme on a larger mesh. 
   */

  if ( E->node_map_exists )
    {
      for ( i=0; i<E->num_nodes; i++)
	{
	  (E->node_map[i])--;
	}
    }

  /*
   * 6. Likewise, if an element map exists, then it might well correspond to
   *    a contiguous integer sequence on a larger mesh. Those numbers could
   *    either be 1-base or 0-based arrays.
   */

  if ( E->elem_map_exists )
    {
      for ( i=0; i<E->num_elems; i++)
	{
	  (E->elem_map[i])--;
	}
    }

  return;
}

/*
 * one_base() -- push up the element names and node names by one in an
 *               EXODUS II data base. This makes C language zero-based
 *               arrays work more naturally.
 *
 *		 This undoes what the zero_base() routine does.
 *
 * Created: 1997/08/26 09:17 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
one_base(Exo_DB *E)
{
  int eb;
  int i,j,l;
  int length_conn;

  /*
   * 1. Node numbers are named in the connectivity lists for each 
   *    element block. Increment them.
   */

  for ( eb=0; eb<E->num_elem_blocks; eb++)
    {
      length_conn =  E->eb_num_elems[eb] * E->eb_num_nodes_per_elem[eb];
      for ( l=0; l<length_conn; l++)
	{
	  (E->eb_conn[eb][l])++;
	}
    }

  /*
   * 2. Node numbers are named in the node lists for each node set.
   *    Increment them.
   */

  for ( l=0; l<E->ns_node_len; l++)
    {
      (E->ns_node_list[l])++;
    }

  /*
   * 3. Element numbers are named in the element list for the side sets.
   *    Increment these.
   */

  for ( l=0; l<E->ss_elem_len; l++)
    {
      (E->ss_elem_list[l])++;
    }

  /*
   * 4. Node numbers in the special node list for the sidesets.
   *    Increment these (if they exist).
   */

  if ( E->ss_node_list_exists )
    {
      for ( i=0; i<E->num_side_sets; i++)
	{
	  
	  /*
	   * Instead of doing a side at a time, do the whole list, whose
	   * length has been counted up and stored in the last element of this
	   * index.
	   */
	  
	  for ( j=0; j<E->ss_node_side_index[i][ E->ss_num_sides[i] ]; j++)
	    {
	      (E->ss_node_list[i][j])++;
	    }
	}
    }



  /*
   * 5. If a node map exists, the names of the nodes might be too high.
   *    Increment them. In some cases for very general node maps this
   *    step might be superfluous and even undesirable. Here, the idea is
   *    that the node map probably corresponds to an inherent contiguous
   *    node numbering scheme on a larger mesh. 
   */

  if ( E->node_map_exists )
    {
      for ( i=0; i<E->num_nodes; i++)
	{
	  (E->node_map[i])++;
	}
    }

  /*
   * 6. Likewise, if an element map exists, then it might well correspond to
   *    a contiguous integer sequence on a larger mesh. Those numbers could
   *    either be 1-base or 0-based arrays.
   */

  if ( E->elem_map_exists )
    {
      for ( i=0; i<E->num_elems; i++)
	{
	  (E->elem_map[i])++;
	}
    }


}


/*
 * fence_post() -- return the index in an integer array where the 
 *                 provided integer value is bounded between fenceposts as
 *
 *			array[index] <= val < array[index+1]
 *
 *		   
 * Assumptions:
 *		[1] The array is monotonically increasing with index.
 *
 *		[2] If val < array[0], then -1 is returned.
 *
 *		[3] If val >= array[length-1], then -1 is returned.
 *
 * (This routine was created to quickly find that element block to which a
 *  given element belongs. 
 *
 * Created: 1997/04/03 07:28 MST pasacki@sandia.gov
 *
 * Revised: 1997/04/21 08:43 MDT pasacki@sandia.gov
 */

int
fence_post(const int val,		/* the integer we seek */
	   int *array,
	   const int length)
{
  int index;
  int first_val, last_val;
  int found;

  double frac;
#ifdef DEBUG
  int i;
#endif

  first_val = array[0];

  last_val = array[length-1];

  if ( first_val == last_val ) return(-1);

  if ( val < first_val )  return(-1);

  if ( val >= last_val )  return(-1);

  if ( val == first_val ) return(0);

  /*
   * Verify monotonicity. Turn off for efficiency later.
   */

#ifdef DEBUG

  for ( i=1; i<length; i++)
    {
      if ( array[i-1] > array[i] )
	{
	  sr = sprintf(err_msg, "Non monotone map where a[%d]=%d is > a[%d]=%d",
		       i-1, array[i-1], i, array[i]);
	  EH(-1, err_msg);
	}
    }

#endif

  /*
   * Linear approximation to first guess.
   */

  frac = ((double)(val - first_val))/((double)(last_val - first_val));

  index = (int)( (double)length * frac);

  index = MIN(index, length-2);

  /*
   * From this starting point, look up or down accordingly.
   */

  found = FALSE;

  while ( ! found )
    {
      if ( array[index+1] <= val ) 
	{
	  index++;
	}
      else if ( array[index] > val )
	{
	  index--;
	}
      else
	{
	  found = TRUE;
	}
    }

#ifdef DEBUG
  fprintf(stdout, "fencepost: val= %d in (a[%d]=%d, a[%d]=%d)\n",
	  val, index, array[index], index+1, array[index+1]);
#endif
  return(index);
}




/* init_exo_struct() -- initialize some defaults...
 * 
 * This is meant to be called right after allocation, to help setup some
 * reasonable defaults to describe an empty data structure. Call it a poor
 * man's constructor.
 *
 * Created: 1998/09/01 10:50 MDT pasacki@sandia.gov
 *
 * Revised:
 */
void
init_exo_struct(Exo_DB *x)
{
  if ( x == NULL )
    {
      EH(-1, "Empty structure to initialize?");
    }

  x->state                 = EXODB_STATE_GRND;
  
  x->num_dim               = 0;
  x->num_elems             = 0;
  x->num_elem_blocks       = 0;
  x->num_info              = 0;
  x->num_node_sets         = 0;
  x->num_side_sets         = 0;

  x->num_nodes             = 0;
  x->num_qa_rec            = 0;

  x->node_map              = NULL;
  x->elem_map              = NULL;
  x->elem_order_map        = NULL;
  x->x_coord               = NULL;
  x->y_coord               = NULL;
  x->z_coord               = NULL;
  x->title                 = NULL;
  x->eb_id                 = NULL;
  x->eb_num_elems          = NULL;
  x->eb_ptr                = NULL;
  x->elem_ptr              = NULL;
  x->node_list             = NULL;
  x->elem_node_pntr        = NULL;
  x->elem_node_list        = NULL;
  x->node_elem_pntr        = NULL;
  x->node_elem_list        = NULL;
  x->node_node_pntr        = NULL;
  x->node_node_list        = NULL;
  x->elem_elem_pntr        = NULL;
  x->elem_elem_list        = NULL;
  x->elem_elem_twst        = NULL;
  x->elem_elem_face        = NULL;
  x->ns_id                 = NULL;
  x->ss_id                 = NULL;

  x->elem_node_conn_exists = FALSE;
  x->node_elem_conn_exists = FALSE;
  x->node_node_conn_exists = FALSE;
  x->elem_elem_conn_exists = FALSE;

  x->node_map_exists       = FALSE;
  x->elem_map_exists       = FALSE;
  x->elem_order_map_exists = FALSE;
  x->ss_node_list_exists   = FALSE;

  x->num_elem_vars         = 0;
  x->num_node_vars         = 0;
  x->num_glob_vars         = 0;

  x->num_nv_time_indeces   = 0;
  x->num_nv_indeces        = 0;

  x->num_ev_time_indeces   = 0;

  x->num_gv_time_indeces   = 0;

  x->num_times             = 0;

  x->time_vals             = NULL;

  /*
   * Let this value indicate that the structure in memory is not currently
   * attached to any open netCDF file. Once open, this will become >-1.
   */

  x->exoid                 = -1;
  return;
}
  
void 
free_exo_ev(Exo_DB *x)
{
  int i;
  int index;
  int j;
  int k;

  if ( ! ( x->state & EXODB_STATE_ELVA ) )
    {
      return;
      /*      EH(-1, "Can't free what was never in chains.");*/
    }

  free(x->ev_time_indeces);

  for ( i=0; i<x->num_ev_time_indeces; i++)
    {


      for ( j=0; j<x->num_elem_blocks; j++)
	{
	  for ( k=0; k<x->num_elem_vars; k++)
	    {
	      index = j * x->num_elem_vars + k;
	      if ( x->elem_var_tab[index] != 0 )
		{
		  free(x->ev[i][index]);
		}
	    }
	}
      free(x->ev[i]);
    }
  free(x->ev);

  /*
   * Indicate that this memory structure no longer has space allocated for
   * element results variables.
   */

  x->state ^= EXODB_STATE_ELVA;

  return;
}

void 
free_exo_gv(Exo_DB *x)
{
  int i;

  if ( ! ( x->state & EXODB_STATE_GBVA ) )
    {
      return;
      /*      EH(-1, "Can't free what was never in chains.");*/
    }

  free(x->gv_time_indeces);

  for ( i=0; i<x->num_gv_time_indeces; i++)
    {
      free(x->gv[i]);
    }

  free(x->gv);

  /*
   * Indicate that this memory structure no longer has space allocated for
   * global results variables.
   */

  x->state ^= EXODB_STATE_GBVA;

  return;
}

void
free_exo_nv(Exo_DB *x)
{
  int i;
  int j;

  if ( ! ( x->state & EXODB_STATE_NDVA ) )
    {
      return;			/* This was a useless call... */
      /*       EH(-1, "Can't free what was never in chains.");*/
    }

  free(x->nv_indeces);

  free(x->nv_time_indeces);

  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  free(x->nv[i][j]);
	}
      free(x->nv[i]);
    }
  free(x->nv);

  /*
   * Indicate that this memory structure no longer has space allocated for
   * nodal results variables.
   */

  x->state ^= EXODB_STATE_NDVA;

  return;
}

/* alloc_exo_ev() -- allocate arrays in exodus structure for element variables
 *
 * 1998/12/21 13:43 MST pasacki@sandia.gov
 */

void
alloc_exo_ev(Exo_DB *x,
	     const int num_timeplanes)
{
  int i;
  int index;
  int j;
  int k;
  int len;

  if ( x->state & EXODB_STATE_ELVA )
    {
      EH(-1, "Please free before allocating...");
    }

  x->num_ev_time_indeces = num_timeplanes;

  x->ev_time_indeces = (int *) smalloc(x->num_ev_time_indeces*sizeof(int));

  x->ev = (dbl ***) smalloc(x->num_ev_time_indeces*sizeof(dbl **));

  /*
   * Unwritten hypothesis: elem_var_vals refers to the last single
   * timeplane that was read. This is a sop for backward usage.
   */
  
  x->elem_var_vals = x->ev[x->num_ev_time_indeces-1];

  for ( i=0; i<x->num_ev_time_indeces; i++)
    {
      x->ev_time_indeces[i] = i+1;

      len = x->num_elem_vars * x->num_elem_blocks;

      x->ev[i] = (dbl **) smalloc(len*sizeof(dbl *));

      for ( j=0; j<x->num_elem_blocks; j++)
	{
	  for ( k=0; k<x->num_elem_vars; k++)
	    {
	      index = j * x->num_elem_vars + k;
	      if ( x->elem_var_tab[index] != 0 )
		{
		  x->ev[i][index] = (dbl *) smalloc(x->eb_num_elems[j]*
						    sizeof(dbl));
		}
	    }
	}
    }

  x->state |= EXODB_STATE_ELVA;
  return;
}

/* alloc_exo_gv() -- allocate arrays in exodus structure for global variables
 *
 * 1998/12/21 13:44 MST pasacki@sandia.gov
 */

void
alloc_exo_gv(Exo_DB *x,
	     const int num_timeplanes)
{
  int i;

  if ( x->state & EXODB_STATE_GBVA )
    {
      EH(-1, "Please free before allocating...");
    }

  x->num_gv_time_indeces = num_timeplanes;

  x->gv_time_indeces     = (int *) smalloc(x->num_gv_time_indeces*sizeof(int));

  /*
   * Fill these indeces with reasonable defaults.
   */

  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      x->gv_time_indeces[i] = i+1;
    }

  /*
   * Allocate the mother of all arrays...
   */

  x->gv = (dbl **) smalloc(x->num_gv_time_indeces*spd);
  for ( i=0; i<x->num_gv_time_indeces; i++)
    {
      x->gv[i] = (dbl *) smalloc(x->num_glob_vars*sd);
    }
	  
  x->state |= EXODB_STATE_GBVA;
  return;
}

/* alloc_exo_nv() -- allocate arrays in exodus structure for nodal variables
 *
 * 1998/12/21 13:45 MST pasacki@sandia.gov
 */

void
alloc_exo_nv(Exo_DB *x,
	     const int num_timeplanes,
	     const int num_nodal_vars)
{
  int i;
  int j;

  if ( x->state & EXODB_STATE_NDVA )
    {
      EH(-1, "Please free before allocating...");
    }

  x->num_nv_indeces      = num_nodal_vars;
  x->num_nv_time_indeces = num_timeplanes;

  x->nv_indeces          = (int *) smalloc(x->num_nv_indeces*sizeof(int));
  x->nv_time_indeces     = (int *) smalloc(x->num_nv_time_indeces*sizeof(int));

  /*
   * Fill these indeces with reasonable defaults.
   */

  for ( i=0; i<x->num_nv_indeces; i++)
    {
      x->nv_indeces[i] = i+1;	/* FORTRAN rulz, dude! */
    }

  for ( j=0; j<x->num_nv_time_indeces; j++)
    {
      x->nv_time_indeces[j] = j+1;
    }

  /*
   * Allocate the mother of all arrays...
   */

  x->nv = (dbl ***) smalloc(x->num_nv_time_indeces*sizeof(dbl **));

  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      x->nv[i] = (dbl **) smalloc(x->num_nv_indeces*sizeof(dbl *));
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  x->nv[i][j] = (dbl *) smalloc(x->num_nodes*sizeof(dbl));
	}
    }

  x->state |= EXODB_STATE_NDVA;
  return;
}
