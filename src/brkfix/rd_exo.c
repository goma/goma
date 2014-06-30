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

/* rd_exo --	open & read an EXODUS II finite element database. 
 *
 * rd_exo() --  input arguments
 *	x -- a pointer to a structure containing the information from
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
 * Modified: 1998/08/04 12:45 MDT pasacki@sandia.gov
 */

#define _RD_EXO_C

#include <config.h>

#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

#include "map_names.h"
#include "std.h"
#include "exo_struct.h"
#include "eh.h"
#include "aalloc.h"
#include "rd_exo.h"

/*
 * Prototypes for functions needed here, defined elsewhere...
 */

extern int in_list              /* utils.c */
PROTO((int ,                    /* val    - what integer value to seek */
       int *,                   /* start  - where to begin looking */
       int ));                  /* length - how far to search from start */


/*
 * Variables used in several routines in this file.
 */

static int sc  = sizeof(char);
static int si  = sizeof(int);
static int sd  = sizeof(dbl);

static int spc = sizeof(char *);
static int spi = sizeof(int *);
static int spf = sizeof(flt *);
static int spd = sizeof(dbl *);

static Spfrtn sr=0;			/* sprintf() return type */
static char err_msg[MAX_CHAR_ERR_MSG];

int 
rd_exo(Exo_DB *x,		/* def'd in exo_struct.h */
       char *fn,
       int verbosity,
       int task)
{
  int i;
  int index;			/* translate elem var, block into 1D index */
  int j;
  int k;
  int len;
  int status=0;
  int nnn;

  char rc;			/* generic returned character (char)*/
  int ri;			/* generic returned integer (int)*/
  flt rf;			/* generic returned flt (float)*/
  /*  dbl rd;			 generic returned dbl (double)*/
  int nodal_var_index;
  int time_index;

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
      x->path = (char *) smalloc(len*sizeof(char));
    }

  strcpy(x->path, fn);

  x->mode          = EX_READ;
  x->comp_wordsize = sd;
  x->io_wordsize   = 0;		/* That is, you tell me. */

  if ( verbosity > 1 )
    {
      fprintf(stderr, "ex_open() call...\n");
    }
  x->exoid = ex_open(x->path, x->mode, &x->comp_wordsize, 
		     &x->io_wordsize, &x->version);
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

      x->title = (char *) smalloc(MAX_LINE_LENGTH*sc);

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

      /*
       * Consider Quality Assurance and Information to be vital initial
       * information. Hey, we're politically correct!
       */

      status = ex_inquire(x->exoid, EX_INQ_QA, &ri, &rf, &rc);
      EH(status, "ex_inquire(EX_INQ_QA)");
      x->num_qa_rec = ri;

      if ( verbosity > 1 )
	{
	  fprintf(stderr, "\tx->num_qa_rec      = %d\n", x->num_qa_rec);
	}      

      if ( x->num_qa_rec > 0 )
	{
	  x->qa_record = (QA_Record *)smalloc(x->num_qa_rec*sizeof(QA_Record));
	  for ( i=0; i<x->num_qa_rec; i++)
	    {
	      for ( j=0; j<4; j++)
		{
		  x->qa_record[i][j] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
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
      
	  x->info = (INFO_Record *) smalloc(x->num_info * sizeof(INFO_Record));

	  for ( i=0; i<x->num_info; i++)
	    {
	      x->info[i] = (char *) smalloc(MAX_LINE_LENGTH * sc);
	    }
	  status = ex_get_info(x->exoid, &(x->info[0]));
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

      if ( x->num_dim > 0 )
	{
	  x->coord_names = (char **) smalloc(x->num_dim * spc);
	  for ( i=0; i<x->num_dim; i++ )
	    {
	      x->coord_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
	    }
	  status = ex_get_coord_names(x->exoid, x->coord_names);
	  EH(status, "ex_get_coord_names");
	}

      if ( x->num_nodes > 0 )
	{
	  if ( x->node_map_exists )
	    {
	      x->node_map = (int *) smalloc(x->num_nodes * si);
	      
	      status = ex_get_node_num_map(x->exoid, x->node_map);
	      EH(status, "ex_get_node_num_map");
	    }
	}
      
      if ( x->num_elems > 0 )
	{
	  if ( x->elem_map_exists )
	    {
	      x->elem_map = (int *) smalloc(x->num_elems * si);
	      status = ex_get_elem_num_map(x->exoid, x->elem_map);
	      EH(status, "ex_get_elem_num_map");
	    }
	  
	  if ( x->elem_order_map_exists )
	    {
	      x->elem_order_map = (int *) smalloc(x->num_elems * si);

	      status = ex_get_map(x->exoid, x->elem_order_map);
	      EH(status, "ex_get_map");
	    }
	}

      /*
       * ELEMENT BLOCKS...
       */
      
      if ( x->num_elem_blocks > 0 )
	{
	  x->eb_id                 = (int *) smalloc(x->num_elem_blocks* si);
	  x->eb_elem_type          = (char **) smalloc(x->num_elem_blocks*spc);
	  x->eb_num_elems          = (int *) smalloc(x->num_elem_blocks*si);
	  x->eb_num_nodes_per_elem = (int *) smalloc(x->num_elem_blocks*si);
	  x->eb_num_attr           = (int *) smalloc(x->num_elem_blocks*si);
	  x->eb_conn               = (int **) smalloc(x->num_elem_blocks*spi);
	  x->eb_attr               = (dbl **) smalloc(x->num_elem_blocks*spf);
	  
	  x->eb_ptr                = (int *)smalloc((x->num_elem_blocks+1)*si);
	  x->eb_ptr[0] = 0;

	  for ( i=0; i<x->num_elem_blocks; i++)
	    {
	      x->eb_elem_type[i] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
	    }

	  status = ex_get_elem_blk_ids(x->exoid, x->eb_id);
	  EH(status, "ex_get_elem_blk_ids");
    
	  for ( i=0; i<x->num_elem_blocks; i++)
	    {
	      status = ex_get_elem_block(x->exoid, 
					 x->eb_id[i], 
					 x->eb_elem_type[i],
					 &x->eb_num_elems[i],
					 &x->eb_num_nodes_per_elem[i],
					 &x->eb_num_attr[i]);
	      EH(status, "ex_get_elem_blocks");

	      if ( (x->eb_num_elems[i] * x->eb_num_nodes_per_elem[i]) > 0 )
		{
		  x->eb_conn[i] = 
		    (int *) smalloc((x->eb_num_elems[i] * 
				     x->eb_num_nodes_per_elem[i])*
				    si);

		  status = ex_get_elem_conn(x->exoid, 
					    x->eb_id[i], 
					    x->eb_conn[i]);
		  EH(status, "ex_get_elem_conn");
		}
	      
	      if ( (x->eb_num_elems[i]*x->eb_num_attr[i]) > 0 )
		{
		  x->eb_attr[i] = (dbl *) 
		    smalloc( (x->eb_num_elems[i]*x->eb_num_attr[i])* sd);
		  status = ex_get_elem_attr(x->exoid, x->eb_id[i], 
					    x->eb_attr[i]);
		  EH(status, "ex_get_elem_attr");
		}
	      x->eb_ptr[i+1] = x->eb_ptr[i] + x->eb_num_elems[i];
	    }
	}

      /*
       * NODE SETS...
       */

      if ( x->num_node_sets > 0 )
	{
	  x->ns_id             = (int *) smalloc( x->num_node_sets* si);
	  x->ns_num_nodes      = (int *) smalloc(x->num_node_sets* si);
	  x->ns_num_distfacts  = (int *) smalloc(x->num_node_sets* si);
	  x->ns_node_index     = (int *) smalloc(x->num_node_sets* si);
	  x->ns_distfact_index = (int *) smalloc(x->num_node_sets* si);
	  
	  if ( x->ns_node_len > 0 )
	    {
	      x->ns_node_list = (int *) smalloc(x->ns_node_len* si);
	    }
	  
	  if ( x->ns_distfact_len > 0 )
	    {
	      x->ns_distfact_list = (dbl *) smalloc(x->ns_distfact_len* sd);
	      status = ex_get_concat_node_sets(x->exoid, x->ns_id, 
					       x->ns_num_nodes,
					       x->ns_num_distfacts,
					       x->ns_node_index,
					       x->ns_distfact_index,
					       x->ns_node_list,
					       x->ns_distfact_list);
	      EH(status, "ex_get_concat_node_sets");
	    }
	}

      /*
       * SIDE SETS...
       */

      if ( x->num_side_sets > 0 ) 
	{
	  x->ss_id             = (int *) smalloc(x->num_side_sets* si);
	  x->ss_num_sides      = (int *) smalloc(x->num_side_sets* si);
	  x->ss_num_distfacts  = (int *) smalloc(x->num_side_sets* si);
	  x->ss_elem_index     = (int *) smalloc(x->num_side_sets* si);
	  x->ss_distfact_index = (int *) smalloc(x->num_side_sets* si);
	  
	  x->ss_node_cnt_list  = (int **) smalloc(x->num_side_sets* spi);
	  x->ss_node_list      = (int **) smalloc(x->num_side_sets* spi);
	  
	  x->ss_node_side_index  = (int **) smalloc(x->num_side_sets* spi);
	  
	  if ( x->ss_elem_len > 0 )
	    {
	      x->ss_elem_list = (int *) smalloc(x->ss_elem_len* si);
	      x->ss_side_list = (int *) smalloc(x->ss_elem_len* si);
	    }

	  if ( x->ss_distfact_len > 0 )
	    {
	      x->ss_distfact_list = (dbl *) smalloc(x->ss_distfact_len* sd);
	    }	    

	  status = ex_get_concat_side_sets(x->exoid, x->ss_id, 
					   x->ss_num_sides,
					   x->ss_num_distfacts,
					   x->ss_elem_index,
					   x->ss_distfact_index,
					   x->ss_elem_list,
					   x->ss_side_list,
					   x->ss_distfact_list);
	  EH(status, "ex_get_concat_side_sets");

	  /*
	   * This information turns out to be useful in constructing more
	   * rapid indeces into the distribution factor array.
	   */

	  x->ss_node_list_exists = TRUE; /* does now! */
      
	  for ( i=0; i<x->num_side_sets; i++)
	    {
	      x->ss_node_cnt_list[i] = (int *) smalloc(x->ss_num_sides[i]* si);

	      x->ss_node_list[i] = (int *) smalloc(x->ss_num_distfacts[i]* si);

	      status = ex_get_side_set_node_list(x->exoid,
						 x->ss_id[i],
						 x->ss_node_cnt_list[i],
						 x->ss_node_list[i]);
#ifdef DEBUG
	      fprintf(stderr, "SSID=%d has %d dfs/nds on %d sides.\n", 
		      x->ss_id[i], x->ss_num_distfacts[i], x->ss_num_sides[i]);
	  
	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "nodes for elem %d side %d = %d\n", 
			  x->ss_elem_list[x->ss_elem_index[i]+j],
			  x->ss_side_list[x->ss_elem_index[i]+j],
			  x->ss_node_cnt_list[i][j]);
		}
#endif

	      /*
	       * Set up quick pointers for nodes on each given side of 
	       * a sideset that can be used later to find exactly where to go 
	       * in the big distribution factor list...
	       */

	      x->ss_node_side_index[i] = (int *)
		smalloc((x->ss_num_sides[i]+1)*si);

	      x->ss_node_side_index[i][0] = 0;

	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  x->ss_node_side_index[i][j+1] = (x->ss_node_side_index[i][j]+
						   x->ss_node_cnt_list[i][j]);
		}

#ifdef DEBUG
	      for ( j=0; j<x->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "SS[%d]=%d, nodes for elem %d, side %d:", i, 
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
       *
       * We try to not really count just "one" property. EXODUS II will
       * fake us out with the "ID" property that it will create anyway.
       * If we count it here, then we'll get double bookkeeping...
       */

      /*
       * Node sets...
       */

      if ( x->ns_num_props > 1 ) 
	{
	  x->ns_prop_name = (char **) smalloc(x->ns_num_props* spc);
	  for ( i=0; i<x->ns_num_props; i++)
	    {
	      x->ns_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_NODE_SET, x->ns_prop_name);
	  EH(status, "ex_get_prop_names(EX_NODE_SET)");

	  x->ns_prop = (int **) smalloc(x->ns_num_props* spi);
	  for ( i=0; i<x->ns_num_props; i++)
	    {
	      x->ns_prop[i] = (int *)smalloc(x->num_node_sets* si);
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

      if ( x->ss_num_props > 1 ) 
	{
	  x->ss_prop_name = (char **) smalloc(x->ss_num_props* spc);
	  for ( i=0; i<x->ss_num_props; i++)
	    {
	      x->ss_prop_name[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_SIDE_SET, x->ss_prop_name);
	  EH(status, "ex_get_prop_names(EX_SIDE_SET)");

	  x->ss_prop = (int **) smalloc(x->ss_num_props* spi);
	  for ( i=0; i<x->ss_num_props; i++)
	    {
	      x->ss_prop[i] = (int *)smalloc(x->num_side_sets* si);
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

      if ( x->eb_num_props > 1 )
	{
	  x->eb_prop_name = (char **) smalloc(x->eb_num_props* spc);
	  for ( i=0; i<x->eb_num_props; i++)
	    {
	      x->eb_prop_name[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	    }
	  status = ex_get_prop_names(x->exoid, EX_ELEM_BLOCK, x->eb_prop_name);
	  EH(status, "ex_get_prop_names(EX_ELEM_BLOCK)");

	  x->eb_prop = (int **) smalloc(x->eb_num_props* spi);
	  for ( i=0; i<x->eb_num_props; i++)
	    {
	      x->eb_prop[i] = (int *)smalloc(x->num_elem_blocks* si);
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
   *
   * Here might be a good place to allocate space for the little arrays
   * describing how much data to read?
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

      status = ex_get_var_param(x->exoid, "g", &x->num_glob_vars);
      EH(status, "ex_get_var_param(g)");

      status = ex_get_var_param(x->exoid, "e", &x->num_elem_vars);
      EH(status, "ex_get_var_param(e)");

      status = ex_get_var_param(x->exoid, "n", &x->num_node_vars);
      EH(status, "ex_get_var_param(n)");
      //printf("NUMBER OF NODAL VARIABLES  = %d\n", x->num_node_vars);

      /*
       * Get the names of the results variables: global, element and nodal.
       */

      if ( x->num_glob_vars > 0 )
	{
	  nnn = MAX(x->num_glob_vars, 1);
	  x->glob_var_names = (char **) smalloc(nnn* spc);
	  for ( i=0; i<x->num_glob_vars; i++)
	    {
	      x->glob_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	    }
	  status = ex_get_var_names(x->exoid, "g", x->num_glob_vars,
				    x->glob_var_names);
	  EH(status, "ex_get_var_names(g)");
	}

      if ( x->num_elem_vars > 0 )
	{
	  x->elem_var_names = (char **) smalloc(x->num_elem_vars* spc);
	  for ( i=0; i<x->num_elem_vars; i++)
	    {
	      x->elem_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	    }
	  status = ex_get_var_names(x->exoid, "e", x->num_elem_vars,
				    x->elem_var_names);
	  EH(status, "ex_get_var_names(e)");
	}

      if ( x->num_node_vars > 0 )
	{
	  x->node_var_names = (char **) smalloc(x->num_node_vars* spc);
	  for ( i=0; i<x->num_node_vars; i++)
	    {
	      x->node_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	    }
	  status = ex_get_var_names(x->exoid, "n", x->num_node_vars,
				    x->node_var_names);
	  //for ( i=0; i<x->num_node_vars; i++) {
	  //   printf("        name %d = %s\n", i, x->node_var_names[i]);
	  // }
	  EH(status, "ex_get_var_names(n)");
	}

	/*
	 * Get time values...
	 */

      if ( x->num_times > 0 )
	{
	  x->time_vals = (dbl *) smalloc(x->num_times* sizeof(dbl));
	  status = ex_get_all_times(x->exoid, x->time_vals);
	  EH(status, "ex_get_all_times");
	}

      /*
       * Get element variable truth table...
       */

      if ( x->num_elem_vars > 0 )
	{
	  x->elem_var_tab = (int *) smalloc(x->num_elem_vars *
					    x->num_elem_blocks*
					    si);

	  status = ex_get_elem_var_tab(x->exoid, x->num_elem_blocks,
				       x->num_elem_vars, x->elem_var_tab);
	  EH(status, "ex_get_elem_var_tab");
	}

      /*
       * There are auxiliary arrays that can be set up to direct the
       * slurping of multiple timeplanes of results in one shot.
       * 
       * Generally, though, it will be expedient to set up and allocate
       * for ONE timeplane and ALL nodal, elemental, and global variables.
       */

      /*
       * Some inactive alternatives for you...
       *
       * x->num_nv_time_indeces = x->num_times; 
       * x->num_nv_indeces      = 1;
       */

#ifdef DEBUG
      fprintf(stderr, 
"exo-> for %s has num_node_vars=%d, nv_time_indeces=%d, nv_indeces=%d\n",
	      x->path, x->num_node_vars, x->num_nv_time_indeces, 
	      x->num_nv_indeces);
#endif

      if ( x->state & EXODB_STATE_NDIA )
	{
	  EH(-1, "Attempt to make nd indeces twice.");
	}

      x->num_nv_time_indeces = 1;
      x->num_nv_indeces      = x->num_node_vars;

      if ( x->num_nv_indeces > 0 )
	{
	  x->nv_indeces          = (int *) smalloc(x->num_nv_indeces*
						   sizeof(int));
	}

      if ( x->num_nv_time_indeces > 0 )
	{
	  x->nv_time_indeces     = (int *) smalloc(x->num_nv_time_indeces*
						   sizeof(int));
	}

      x->state              |= EXODB_STATE_NDIA;

      /*
       * Fill these indeces with reasonable defaults.
       *
       * WARNING!!! On 2nd and subsequent calls, you will probably want
       * to increment x->nv_time_indeces[0] to look at the subsequent
       * timeplanes of results!
       */

      for ( i=0; i<x->num_nv_indeces; i++)
	{
	  x->nv_indeces[i] = i+1;	/* Names are 1,2,...,n -- FORTRAN rulz,
					 * dude! */
	}

      for ( j=0; j<x->num_nv_time_indeces; j++)
	{
	  x->nv_time_indeces[j] = j+1;
	}

      x->state |= EXODB_STATE_RES0;
    }

  /*
   * Get nodal results.
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
	      /*
	       * Read only one timeplane's worth of nodal variables at 
	       * a time, but read all of the nodal variables that are 
	       * available. To do so, allocate enough space. Note that
	       * allocation of x->nv[][][] will use info in
	       * x->num_nv_time_indeces, etc.
	       */

	      alloc_exo_nv(x);
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
		  
		  status = ex_get_nodal_var(x->exoid, time_index, 
					    nodal_var_index, x->num_nodes, 
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
	      status = ex_get_glob_vars(x->exoid, time_index, 
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
			  status = ex_get_elem_var(x->exoid, time_index, k+1,
						   x->eb_id[j], 
						   x->eb_num_elems[j],
						   x->ev[i][index]);
			  if ( status < 0 )
			    {
			      sr = sprintf(err_msg, 
   "ex_get_elem_var() bad rtn: time %d, elemvar %d, EB ID %d",
					   time_index, k+1, x->eb_id[j]);
			      EH(-1, err_msg);
			      EH(sr, err_msg);
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

/*
 * copy_exo() -- replicate via allocation and copying an EXODUS II db
 *
 * The intent is to fill in a second, destination EXODUS II db structure
 * with contents identical to the source. Then, you can muck with it
 * as you please.
 *
 * Be sure to free_exo() both the source and the destination databases
 * after you're done with them!
 *
 * This is *a lot* like rd_exo(), that also does allocation and sets values.
 * The main difference is that the source is another Exo_DB structure instead
 * of a disk file accessed through the EXODUS II API.
 *
 * The action directives are used here to determine what data gets transferred
 * from the source to the destination!
 *
 *
 * Created: 1998/08/05 10:44 MDT pasacki@sandia.gov
 *
 * Revised:
 */

int
copy_exo(Exo_DB *src, 
	 Exo_DB *dst,
	 int verbosity,
	 int task)
{
  int i;
  int index;
  int j;
  int k;
  int l;
  int elem_blk_index=0;
  int status=0;

  /*
   * First executable statement.
   */

  dst->state = EXODB_STATE_GRND;

  if ( verbosity > 0 )
    {
      fprintf(stderr, "copy_exo() begins.\n");
    }

  /*
   * Setup some Boolean defaults. To begin, we have nothing.
   *
   * Better do this at a higher level !!!!
   */

  if ( src->state == EXODB_STATE_GRND )
    {
      return(status);
    }
  
  if ( task & EXODB_ACTION_RD_INIT )
    {
      
      /*
       * The source should be reliable - 1st law of journalism.
       */

      if ( ! ( src->state & EXODB_STATE_INIT ) )
	{
	  EH(-1, "The source Exo DB is not init state; cannot make dst so.");
	}

      /*
       * Check to see if we're trying to build on top of an existing
       * database? We really need to start from the ground zero state.
       *
       * Well, maybe not. Sometimes you could want to overwrite, but
       * since we do allocation intermixed, we won't attempt to handle
       * the overwrite case here and now.
       */

      if ( dst->state != EXODB_STATE_GRND )
	{
	  EH(-1, "Attempt to copy into nonzero destination.");
	}

      /*
       * ex_open() -- preliminary information
       */

      dst->path = (char *) smalloc(FILENAME_MAX_ACK*sc);
      strcpy(dst->path, src->path);

      SRC_DST(mode);
      SRC_DST(comp_wordsize);
      SRC_DST(io_wordsize);
      SRC_DST(version);

      /*
       * ex_get_init() -- basic information
       */

      dst->title         = (char *) smalloc(MAX_LINE_LENGTH*sc);
      strcpy(dst->title, src->title);

      SRC_DST(num_dim);
      SRC_DST(num_nodes);
      SRC_DST(num_elems);
      SRC_DST(num_elem_blocks);
      SRC_DST(num_node_sets);
      SRC_DST(num_side_sets);

      if ( verbosity > 0 )
	{
	  fprintf(stderr, "\tdst->title           = \"%s\"\n", dst->title);
	  fprintf(stderr, "\tdst->num_nodes       = %d\n", dst->num_nodes);
	  fprintf(stderr, "\tdst->num_elems       = %d\n", dst->num_elems);
	  fprintf(stderr, "\tdst->num_elem_blocks = %d\n", 
		  dst->num_elem_blocks);
	  fprintf(stderr, "\tdst->num_node_sets   = %d\n", dst->num_node_sets);
	  fprintf(stderr, "\tdst->num_side_sets   = %d\n", dst->num_side_sets);
	}

      /*
       * ex_inquire() -- oddball quantities like number of QA records and
       *                 number of sideset distribution factors, etc.
       *
       * QA Records...
       */

      SRC_DST(num_qa_rec);
      if ( dst->num_qa_rec > 0 )
	{
	  dst->qa_record = (QA_Record *) smalloc(dst->num_qa_rec*
						 sizeof(QA_Record));
	  for ( i=0; i<dst->num_qa_rec; i++)
	    {
	      for ( j=0; j<4; j++)
		{
		  dst->qa_record[i][j] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
		  strcpy(dst->qa_record[i][j], src->qa_record[i][j]);
		}
	    }
	}

      /*
       * INFO Records...
       */

      SRC_DST(num_info);

      if ( dst->num_info > 0 )
	{
	  dst->info = (INFO_Record *) smalloc(dst->num_info*
					      sizeof(INFO_Record));
	  for ( i=0; i<dst->num_info; i++)
	    {
	      dst->info[i] = (char *) smalloc(MAX_LINE_LENGTH * sc);
	      strcpy(dst->info[i], src->info[i]);
	    }
	}

      /*
       * Miscellaneous...
       */

      SRC_DST(api_version);
      SRC_DST(db_version);
      SRC_DST(ns_node_len);
      SRC_DST(ns_distfact_len);
      SRC_DST(ss_elem_len);
      SRC_DST(ss_distfact_len);
      SRC_DST(ss_node_len);
      SRC_DST(eb_num_props);
      SRC_DST(ns_num_props);
      SRC_DST(ss_num_props);
      SRC_DST(num_times);

      dst->state |= EXODB_STATE_INIT;		/* Did it! */
    }

  if ( task & EXODB_ACTION_RD_MESH )
    {
      
      if ( ! ( dst->state & EXODB_STATE_INIT ) )
	{
	  EH(-1, "Need to rd init info from EXO db before mesh.");
	}

      if ( dst->state & EXODB_STATE_MESH )
	{
	  EH(-1, "Attempt to rd mesh EXO db into previously used struct.");
	}

      /*
       * Check for reliable source again...
       */

      if ( ! ( src->state & EXODB_STATE_MESH ) )
	{
	  EH(-1, "Attempt to copy mesh EXO db from unmeshed source.");
	}	  

      /*
       * Spatial coordinates of nodes...
       */

      if ( dst->num_dim > 0 )
	{
	  dst->x_coord = (dbl *) smalloc(dst->num_nodes * sd);
	  for ( i=0; i<dst->num_nodes; i++)
	    {
	      SRC_DST(x_coord[i]);
	    }
	}
      
      if ( dst->num_dim > 1 )
	{
	  dst->y_coord = (dbl *) smalloc(dst->num_nodes * sd);
	  for ( i=0; i<dst->num_nodes; i++)
	    {
	      SRC_DST(y_coord[i]);
	    }
	}
      
      if ( dst->num_dim > 2 )
	{
	  dst->z_coord = (dbl *) smalloc(dst->num_nodes * sd);
	  for ( i=0; i<dst->num_nodes; i++)
	    {
	      SRC_DST(z_coord[i]);
	    }
	}

      /*
       * Coordinate names.
       */

      if ( dst->num_dim > 0 )
	{
	  dst->coord_names = (char **) smalloc(dst->num_dim * spc);
	  for ( i=0; i<dst->num_dim; i++ )
	    {
	      dst->coord_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
	      strcpy(dst->coord_names[i], src->coord_names[i]);
	    }
	}

      /*
       * Node map.
       */

      SRC_DST(node_map_exists);

      if ( dst->num_nodes > 0 )
	{
	  if ( dst->node_map_exists )
	    {
	      dst->node_map = (int *) smalloc(dst->num_nodes * si);
	      for ( i=0; i<dst->num_nodes; i++)
		{
		  SRC_DST(node_map[i]);
		}
	    }
	}

      /*
       * Element map.
       */
      
      SRC_DST(elem_map_exists);

      if ( dst->num_elems > 0 )
	{
	  if ( dst->elem_map_exists )
	    {
	      dst->elem_map = (int *) smalloc(dst->num_elems * si);
	      for ( i=0; i<dst->num_elems; i++)
		{
		  SRC_DST(elem_map[i]);
		}
	    }
	}

      /*
       * Element order map.
       */
      SRC_DST(elem_order_map_exists);

      if ( dst->num_elems > 0 )
	{
	  if ( dst->elem_order_map_exists )
	    {
	      dst->elem_order_map = (int *) smalloc(dst->num_elems * si);
	      for ( i=0; i<dst->num_elems; i++)
		{
		  SRC_DST(elem_order_map[i]);
		}
	    }
	}

      /*
       * ELEMENT BLOCKS...
       */
      
      if ( dst->num_elem_blocks > 0 )
	{
	  dst->eb_id                 = (int *) smalloc(dst->num_elem_blocks*
						       si);
	  dst->eb_elem_type          = (char **) smalloc(dst->num_elem_blocks*
							 spc);
	  dst->eb_num_elems          = (int *) smalloc(dst->num_elem_blocks*
						       si);
	  dst->eb_num_nodes_per_elem = (int *) smalloc(dst->num_elem_blocks*
						       si);
	  dst->eb_num_attr           = (int *) smalloc(dst->num_elem_blocks*
						       si);
	  dst->eb_conn               = (int **) smalloc(dst->num_elem_blocks*
							spi);
	  dst->eb_attr               = (dbl **) smalloc(dst->num_elem_blocks*
							spf);
	  dst->eb_ptr                = (int *)smalloc((dst->num_elem_blocks+1)*
						      si);
	  dst->eb_ptr[0] = 0;

	  for ( i=0; i<dst->num_elem_blocks; i++)
	    {
	      dst->eb_elem_type[i] = (char *) smalloc((MAX_STR_LENGTH+1)*sc);
	    }

	  for ( i=0; i<dst->num_elem_blocks; i++)
	    {
	      SRC_DST(eb_id[i]);
	      strcpy(dst->eb_elem_type[i], src->eb_elem_type[i]);
	      SRC_DST(eb_num_elems[i]);
	      SRC_DST(eb_num_nodes_per_elem[i]);
	      SRC_DST(eb_num_attr[i]);
	    }

	  for ( i=0; i<dst->num_elem_blocks; i++)
	    {
	      if ( (dst->eb_num_elems[i] * dst->eb_num_nodes_per_elem[i]) > 0 )
		{
		  dst->eb_conn[i] = 
		    (int *) smalloc((dst->eb_num_elems[i] * 
				     dst->eb_num_nodes_per_elem[i])*si);
		  for ( j=0; 
			j<(dst->eb_num_elems[i]*dst->eb_num_nodes_per_elem[i]);
			j++ )
		    {
		      SRC_DST(eb_conn[i][j]);
		    }
		}
	      
	      if ( (dst->eb_num_elems[i]*dst->eb_num_attr[i]) > 0 )
		{
		  dst->eb_attr[i] = (dbl *) 
		    smalloc( (dst->eb_num_elems[i]*dst->eb_num_attr[i])* sd);

		  for ( j=0; j<(dst->eb_num_elems[i]*dst->eb_num_attr[i]); j++)
		    {
		      SRC_DST(eb_attr[i][j]);
		    }
		}
	      SRC_DST(eb_ptr[i+1]);
	    }
	}
    
      /*
       * Node sets...
       */

      if ( dst->num_node_sets > 0 )
	{
	  dst->ns_id             = (int *) smalloc(dst->num_node_sets*si);
	  dst->ns_num_nodes      = (int *) smalloc(dst->num_node_sets*si);
	  dst->ns_num_distfacts  = (int *) smalloc(dst->num_node_sets*si);
	  dst->ns_node_index     = (int *) smalloc(dst->num_node_sets*si);
	  dst->ns_distfact_index = (int *) smalloc(dst->num_node_sets*si);
	  
	  if ( dst->ns_node_len > 0 )
	    {
	      dst->ns_node_list = (int *) smalloc(dst->ns_node_len*si);
	    }
	  
	  if ( dst->ns_distfact_len > 0 )
	    {
	      dst->ns_distfact_list = (dbl *) smalloc(dst->ns_distfact_len*sd);
	    }

	  for ( i=0; i<dst->num_node_sets; i++)
	    {
	      SRC_DST(ns_id[i]);
	      SRC_DST(ns_num_nodes[i]);
	      SRC_DST(ns_num_distfacts[i]);
	      SRC_DST(ns_node_index[i]);
	      SRC_DST(ns_distfact_index[i]);
	    }
	      
	  for ( i=0; i<dst->ns_node_len; i++)
	    {
	      SRC_DST(ns_node_list[i]);
	    }

	  for ( i=0; i<dst->ns_distfact_len; i++)
	    {
	      SRC_DST(ns_distfact_list[i]);
	    }
	}

      /*
       * Side sets...
       */

      if ( dst->num_side_sets > 0 ) 
	{
	  dst->ss_id             = (int *) smalloc(dst->num_side_sets*si);
	  dst->ss_num_sides      = (int *) smalloc(dst->num_side_sets*si);
	  dst->ss_num_distfacts  = (int *) smalloc(dst->num_side_sets*si);
	  dst->ss_elem_index     = (int *) smalloc(dst->num_side_sets*si);
	  dst->ss_distfact_index = (int *) smalloc(dst->num_side_sets*si);
	  
	  dst->ss_node_cnt_list  = (int **) smalloc(dst->num_side_sets*spi);
	  dst->ss_node_list      = (int **) smalloc(dst->num_side_sets*spi);
	  
	  dst->ss_node_side_index  = (int **) smalloc(dst->num_side_sets*spi);
	  
	  if ( dst->ss_elem_len > 0 )
	    {
	      dst->ss_elem_list = (int *) smalloc(dst->ss_elem_len*si);
	      dst->ss_side_list = (int *) smalloc(dst->ss_elem_len*si);
	    }

	  if ( dst->ss_distfact_len > 0 )
	    {
	      dst->ss_distfact_list = (dbl *) smalloc(dst->ss_distfact_len*sd);
	    }	    

	  for ( i=0; i<dst->num_side_sets; i++)
	    {
	      SRC_DST(ss_id[i]);
	      SRC_DST(ss_num_sides[i]);
	      SRC_DST(ss_num_distfacts[i]);
	      SRC_DST(ss_elem_index[i]);
	      SRC_DST(ss_distfact_index[i]);
	    }

	  for ( i=0; i<dst->ss_elem_len; i++)
	    {
	      SRC_DST(ss_elem_list[i]);
	      SRC_DST(ss_side_list[i]);
	    }

	  for ( i=0; i<dst->ss_distfact_len; i++)
	    {
	      SRC_DST(ss_distfact_list[i]);
	    }

	  /*
	   * This information turns out to be useful in constructing more
	   * rapid indeces into the distribution factor array.
	   */

	  dst->ss_node_list_exists = TRUE;
      
	  for ( i=0; i<dst->num_side_sets; i++)
	    {
	      dst->ss_node_cnt_list[i] = 
		(int *) smalloc(dst->ss_num_sides[i]*si);

	      dst->ss_node_list[i] = 
		(int *) smalloc(dst->ss_num_distfacts[i]*si);

	      for ( j=0; j<dst->ss_num_sides[i]; j++)
		{
		  SRC_DST(ss_node_cnt_list[i][j]);
		}

	      for ( j=0; j<dst->ss_num_distfacts[i]; j++)
		{
		  SRC_DST(ss_node_list[i][j]);
		}

#ifdef DEBUG
	      fprintf(stderr, "SSID=%d has %d dfs/nds on %d sides.\n", 
		      dst->ss_id[i], dst->ss_num_distfacts[i], 
		      dst->ss_num_sides[i]);
	  
	      for ( j=0; j<dst->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "nodes for elem %d side %d = %d\n", 
			  dst->ss_elem_list[dst->ss_elem_index[i]+j],
			  dst->ss_side_list[dst->ss_elem_index[i]+j],
			  dst->ss_node_cnt_list[i][j]);
		}
#endif

	      /*
	       * Set up quick pointers for nodes on each given side of 
	       * a sideset that can be used later to find exactly where to go 
	       * in the big distribution factor list...
	       */

	      dst->ss_node_side_index[i] = (int *)
		smalloc((dst->ss_num_sides[i]+1)*si);

	      for ( j=0; j<(dst->ss_num_sides[i]+1); j++)
		{
		  SRC_DST(ss_node_side_index[i][j]);
		}
#ifdef DEBUG
	      for ( j=0; j<dst->ss_num_sides[i]; j++)
		{
		  fprintf(stderr, "SS[%d]=%d, nodes for elem %d, side %d:", i, 
			  dst->ss_id[i], 
			  dst->ss_elem_list[dst->ss_elem_index[i]+j],
			  dst->ss_side_list[dst->ss_elem_index[i]+j]);
		  for ( k=dst->ss_node_side_index[i][j]; 
			k<dst->ss_node_side_index[i][j+1]; k++)
		    {
		      fprintf(stderr, " %d", dst->ss_node_list[i][k]);
		    }
		  fprintf(stderr, "\n");
		}
#endif
	    }
	}

      /*
       * Properites of node sets...
       */

      if ( dst->ns_num_props > 0 ) 
	{
	  dst->ns_prop_name = (char **) smalloc(dst->ns_num_props* spc);
	  for ( i=0; i<dst->ns_num_props; i++)
	    {
	      dst->ns_prop_name[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->ns_prop_name[i], src->ns_prop_name[i]);
	    }

	  dst->ns_prop = (int **) smalloc(dst->ns_num_props* spi);
	  for ( i=0; i<dst->ns_num_props; i++)
	    {
	      dst->ns_prop[i] = (int *)smalloc(dst->num_node_sets* si);
	      for ( j=0; j<dst->num_node_sets; j++)
		{
		  SRC_DST(ns_prop[i][j]);
		}
	    }
	}
      
      /*
       * Properties of side sets...
       */

      if ( dst->ss_num_props > 0 ) 
	{
	  dst->ss_prop_name = (char **) smalloc(dst->ss_num_props* spc);
	  for ( i=0; i<dst->ss_num_props; i++)
	    {
	      dst->ss_prop_name[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->ss_prop_name[i], src->ss_prop_name[i]);
	    }

	  dst->ss_prop = (int **) smalloc(dst->ss_num_props* spi);
	  for ( i=0; i<dst->ss_num_props; i++)
	    {
	      dst->ss_prop[i] = (int *)smalloc(dst->num_side_sets* si);
	      for ( j=0; j<dst->num_side_sets; j++)
		{
		  SRC_DST(ss_prop[i][j]);
		}
	    }
	}
      
      /*
       * Properties of element blocks...
       */

      if ( dst->eb_num_props > 0 ) 
	{
	  dst->eb_prop_name = (char **) smalloc(dst->eb_num_props* spc);
	  for ( i=0; i<dst->eb_num_props; i++)
	    {
	      dst->eb_prop_name[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->eb_prop_name[i], src->eb_prop_name[i]);
	    }

	  dst->eb_prop = (int **) smalloc(dst->eb_num_props* spi);
	  for ( i=0; i<dst->eb_num_props; i++)
	    {
	      dst->eb_prop[i] = (int *)smalloc(dst->num_elem_blocks* si);
	      for ( j=0; j<dst->num_elem_blocks; j++)
		{
		  SRC_DST(eb_prop[i][j]);
		}
	    }
	}
      
      dst->state |= EXODB_STATE_MESH; /* Did it! */
    }

  /*
   * Result preliminaries... meta-data...
   *
   * In particular, read number of nodal, element, global variables 
   * as well as their names.
   *
   * Also, read the number of timeplanes and the element variable 
   * truth table.
   */

  if ( task & EXODB_ACTION_RD_RES0 )
    {

      /*
       * Check for sanity - better not be doing this unless a fair
       * amount of information is already known about this database.
       */

      if ( ! ( dst->state & EXODB_STATE_INIT ) )
	{
	  EH(-1, "Need to rd init info from EXO db before result metadata.");
	}

      /*
       * Hmmm... is it really just the element block IDs for element
       * variables that are required? Perhaps for nodal variables it suffices
       * to know merely the number of nodes...
       */

      if ( ! ( dst->state & EXODB_STATE_MESH ) )
	{
	  EH(-1, "Need to rd mesh info from EXO db before result metadata.");
	}

      if ( dst->state & EXODB_STATE_RES0 )
	{
	  EH(-1, "Attempt to rd mesh EXO db into previously used struct.");
	}

      SRC_DST(num_glob_vars);
      SRC_DST(num_elem_vars);
      SRC_DST(num_node_vars);

      /*
       * Names of the results variables: global, element and nodal.
       */

      if ( dst->num_glob_vars > 0 )
	{
	  dst->glob_var_names = (char **) smalloc(dst->num_glob_vars* spc);
	  for ( i=0; i<dst->num_glob_vars; i++)
	    {
	      dst->glob_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->glob_var_names[i], src->glob_var_names[i]);
	    }
	}

      if ( dst->num_elem_vars > 0 )
	{
	  dst->elem_var_names = (char **) smalloc(dst->num_elem_vars* spc);
	  for ( i=0; i<dst->num_elem_vars; i++)
	    {
	      dst->elem_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->elem_var_names[i], src->elem_var_names[i]);
	    }
	}

      if ( dst->num_node_vars > 0 )
	{
	  dst->node_var_names = (char **) smalloc(dst->num_node_vars* spc);
	  for ( i=0; i<dst->num_node_vars; i++)
	    {
	      dst->node_var_names[i] = (char *) smalloc((MAX_STR_LENGTH+1)* sc);
	      strcpy(dst->node_var_names[i], src->node_var_names[i]);
	    }
	}

      /*
       * Time values...
       */

      if ( dst->num_times > 0 )
	{
	  dst->time_vals = (dbl *) smalloc(dst->num_times* sizeof(dbl));
	  for ( i=0; i<dst->num_times; i++)
	    {
	      SRC_DST(time_vals[i]);
	    }
	}

      /*
       * Element variable truth table...
       */

      if ( dst->num_elem_vars > 0 )
	{
	  dst->elem_var_tab = (int *) smalloc(dst->num_elem_vars *
					      dst->num_elem_blocks * si);
	  for ( i=0; i<(dst->num_elem_vars*dst->num_elem_blocks); i++)
	    {
	      SRC_DST(elem_var_tab[i]);
	    }
	}

      dst->state |= EXODB_STATE_RES0;
    }

  /*
   * Results data: nodal, element and global according to what
   * has been indicated in the directions...
   */

  if ( ( task & EXODB_ACTION_RD_RESN ) ||
       ( task & EXODB_ACTION_RD_RESE ) ||
       ( task & EXODB_ACTION_RD_RESG ) )
    {

      if ( !(dst->state & EXODB_STATE_RES0) )
	{
	  EH(-1, "Need to rd result metadata from EXO db before result data.");
	}

    }

  /*
   * Nodal variable results -- kind of reading from src?
   */
  
  if ( task & EXODB_ACTION_RD_RESN )
    {

      if ( ! ( dst->state & EXODB_STATE_NDVA ) )
	{
	  dst->num_nv_time_indeces = dst->num_times;
	  dst->num_nv_indeces      = dst->num_node_vars;
	  alloc_exo_nv(dst);
	}

      if ( dst->num_node_vars > dst->num_nv_indeces )
	{
	  EH(-1, "Request for more nodal vars than exist in EXO db.");
	}

      /*
       * ??? Maybe need to copy over the specifications???
       */

      if ( dst->num_nv_time_indeces > 0 &&
	   dst->num_nv_indeces      > 0 &&
	   dst->num_nodes           > 0 &&
	   dst->num_node_vars       > 0 )
	{

	  for ( i=0; i<dst->num_nv_time_indeces; i++)
	    {
	      for ( j=0; j<dst->num_nv_indeces; j++)
		{
		  for ( k=0; k<dst->num_nodes; k++)
		    {
		      SRC_DST(nv[i][j][k]);
		    }
		}
	    }
	}
      dst->state |= EXODB_STATE_NDVR;
    }


  if ( task & EXODB_ACTION_RD_RESG )
    {

      if ( ! ( dst->state & EXODB_STATE_GBVA ) )
	{
	  alloc_exo_gv(dst, 1);
	}

      if ( dst->num_gv_time_indeces > 0 &&
	   dst->num_glob_vars > 0 )
	{

	  for ( i=0; i<dst->num_gv_time_indeces; i++)
	    {
	      for ( j=0; j<dst->num_glob_vars; j++)
		{
		  SRC_DST(gv[i][j]);
		}
	    }
	}
      dst->state |= EXODB_STATE_GBVR;
    }

  /*
   * Element variable results.
   */

  if ( task & EXODB_ACTION_RD_RESE )
    {

      if ( ! ( dst->state & EXODB_STATE_ELVA ) )
	{
	  alloc_exo_ev(dst, 1);
	}

      if ( dst->num_ev_time_indeces > 0 &&
	   dst->num_elem_vars       > 0 &&
	   dst->num_elems	    > 0 )
	{

	  for ( i=0; i<dst->num_ev_time_indeces; i++)
	    {
	      for ( j=0; j<dst->num_elem_blocks; j++)
		{
		  for ( k=0; k<dst->num_elem_vars; k++)
		    {
		      index = j * dst->num_elem_vars + k;
		      if ( dst->elem_var_tab[index] != 0 )
			{
			  for ( l=0; 
				l<dst->eb_num_elems[elem_blk_index]; l++)
			    {
			      SRC_DST(ev[i][index][l]);
			    }
			}
		    }
		}
	    }
	}
      dst->state |= EXODB_STATE_ELVR;
    }



  return(status);
} /* end copy_exo() */


int
free_exo(Exo_DB *x)		/* pointer to EXODUS II FE db structure */
{
  int i;
  int j;

  safe_free(x->title);

  safe_free(x->path);

  if ( x->num_qa_rec > 0 )
    {
      for ( i=0; i<x->num_qa_rec; i++)
	{
	  for ( j=0; j<4; j++)
	    {
	      safe_free(x->qa_record[i][j]);
	    }
	}
      safe_free(x->qa_record);
    }

  if ( x->num_info > 0 )
    {
      for ( i=0; i<x->num_info; i++)
	{
	  safe_free(x->info[i]);
	}
      safe_free(x->info);
    }

  if ( x->elem_map_exists )
    {
      safe_free(x->elem_map);
    }

  if ( x->node_map_exists )
    {
      safe_free(x->node_map);
    }

  if ( x->num_dim > 0 )
    {
      safe_free(x->x_coord);
    }

  if ( x->num_dim > 1 )
    {
      safe_free(x->y_coord);
    }
    
  if ( x->num_dim > 2 )
    {
      safe_free(x->z_coord);
    }

  if ( x->num_dim > 0 )
    {
      for ( i=0; i<x->num_dim; i++ )
	{
	  safe_free(x->coord_names[i]);
	}
      safe_free(x->coord_names);
    }

  if ( x->num_elem_blocks > 0 )
    {
      for ( i=0; i<x->num_elem_blocks; i++)
	{
	  safe_free(x->eb_elem_type[i]);
	  safe_free(x->eb_conn[i]);
	  if ( (x->eb_num_elems[i]*x->eb_num_attr[i]) > 0 )
	    {
	      safe_free(x->eb_attr[i]);
	    }
	}

      safe_free(x->eb_id);
      safe_free(x->eb_num_elems);
      safe_free(x->eb_num_nodes_per_elem);
      safe_free(x->eb_num_attr);

      safe_free(x->eb_ptr);
      safe_free(x->eb_elem_type);
      safe_free(x->eb_conn);
      safe_free(x->eb_attr);
    }

  if ( x->num_node_sets > 0 )
    {
      safe_free(x->ns_id);
      safe_free(x->ns_num_nodes);
      safe_free(x->ns_num_distfacts);
      safe_free(x->ns_node_index);
      safe_free(x->ns_distfact_index);

      if ( x->ns_node_len > 0 )
	{
	  safe_free(x->ns_node_list);
	}

      if ( x->ns_distfact_len > 0 )
	{
	  safe_free(x->ns_distfact_list);
	}
    }

  if ( x->num_side_sets > 0 ) 
    {
      
      if ( x->ss_node_list_exists )
	{
	  for ( i=0; i<x->num_side_sets; i++)
	    {
	      safe_free(x->ss_node_cnt_list[i]);
	      safe_free(x->ss_node_list[i]);
	      safe_free(x->ss_node_side_index[i]);
	    }
	  safe_free(x->ss_node_cnt_list);
	  safe_free(x->ss_node_list);
	  safe_free(x->ss_node_side_index);
	}

      safe_free(x->ss_id);
      safe_free(x->ss_num_sides);
      safe_free(x->ss_num_distfacts);
      safe_free(x->ss_elem_index);
      safe_free(x->ss_distfact_index);

      if ( x->ss_elem_len > 0 )
	{
	  safe_free(x->ss_elem_list);
	  safe_free(x->ss_side_list);
	}

      if ( x->ss_distfact_len > 0 )
	{
	  safe_free(x->ss_distfact_list);
	}
    }

  /*
   * PROPERTIES...
   */

  /*
   * Node sets...
   */

  if ( x->ns_num_props > 1 && x->num_node_sets > 0 ) 
    {
      for ( i=0; i<x->ns_num_props; i++)
	{
	  safe_free(x->ns_prop_name[i]);
	  safe_free(x->ns_prop[i]);
	}
      safe_free(x->ns_prop_name);
      safe_free(x->ns_prop);
    }
      
  /*
   * Side sets...
   */

  if ( x->ss_num_props > 1 && x->num_side_sets > 0 ) 
    {
      for ( i=0; i<x->ss_num_props; i++)
	{
	  safe_free(x->ss_prop_name[i]);
	  safe_free(x->ss_prop[i]);
	}
      safe_free(x->ss_prop_name);
      safe_free(x->ss_prop);
    }
      
  /*
   * Element blocks...
   */

  if ( x->eb_num_props > 1 && x->num_elem_blocks > 0 ) 
    {
      for ( i=0; i<x->eb_num_props; i++)
	{
	  safe_free(x->eb_prop_name[i]);
	  safe_free(x->eb_prop[i]);
	}
      safe_free(x->eb_prop_name);
      safe_free(x->eb_prop);
    }
      
  /*
   * Results variables...
   */

  if ( x->num_glob_vars > 0 )
    {
      for ( i=0; i<x->num_glob_vars; i++)
	{
	  safe_free(x->glob_var_names[i]);
	}
      safe_free(x->glob_var_names);
    }



  if ( x->num_elem_vars > 0 )
    {
      if ( x->num_ev_time_indeces > 0 )
	{
	  free_exo_ev(x);
	}
      for ( i=0; i<x->num_elem_vars; i++)
	{
	  safe_free(x->elem_var_names[i]);
	}
      safe_free(x->elem_var_names);
      safe_free(x->elem_var_tab);
    }

  if ( x->num_node_vars > 0 )
    {
      for ( i=0; i<x->num_node_vars; i++)
	{
	  safe_free(x->node_var_names[i]);
	}
      safe_free(x->node_var_names);
    }

  if ( x->num_times > 0 )
    {
      safe_free(x->time_vals);
    }

  free_exo_nv(x);
  free_exo_ev(x);
  free_exo_gv(x);
  
  if ( x->state & EXODB_STATE_NDVA )
    {
      EH(-1, "Nodal results values space still allocated!");
    }

  if ( x->state & EXODB_STATE_NDIA )
    {
      EH(-1, "Nodal results indeces space still allocated!");
    }

  if ( x->state & EXODB_STATE_ELIA )
    {
      EH(-1, "Element results indeces space still allocated!");
    }

  if ( x->state & EXODB_STATE_ELVA )
    {
      EH(-1, "Element results values space still allocated!");
    }

  if ( x->state & EXODB_STATE_GBVA )
    {
      EH(-1, "Global results values space still allocated!");
    }

  if ( x->state & EXODB_STATE_GBIA )
    {
      EH(-1, "Global results indeces space still allocated!");
    }

  /*
   * Free up any connectivity arrays...
   */

  if ( x->node_node_conn_exists )
    {
      free(x->node_node_pntr);
      free(x->node_node_list);
    }

  if ( x->elem_node_conn_exists )
    {
      free(x->elem_node_pntr);
      free(x->elem_node_list);
    }

  if ( x->elem_elem_conn_exists )
    {
      free(x->elem_elem_pntr);
      free(x->elem_elem_list);
      free(x->elem_elem_face);
      free(x->elem_elem_twst);
    }

  if ( x->node_elem_conn_exists )
    {
      free(x->node_elem_pntr);
      free(x->node_elem_list);
    }

  return(0);
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

  /*
   * Set all elements in the structure to NULL to ensure that all
   * pointers are set to NULL. Thus, we can check against NULL
   * to see if a deallocation is necessary.
   */
  (void) memset((void *)x, 0, sizeof(Exo_DB));
  
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

  if ( x->state &  EXODB_STATE_ELIA )
    {
      free(x->ev_time_indeces);
      x->state ^= EXODB_STATE_ELIA;
    }

  if ( x->state &  EXODB_STATE_ELVA )
    {
      for ( i=0; i<x->num_ev_time_indeces; i++)
	{
	  for ( j=0; j<x->num_elem_blocks; j++)
	    {
	      for ( k=0; k<x->num_elem_vars; k++)
		{
		  index = j * x->num_elem_vars + k;
		  if ( x->elem_var_tab[index] != 0 )
		    {
		      safe_free(x->ev[i][index]);
		    }
		}
	    }
	  safe_free(x->ev[i]);
	}
      safe_free(x->ev);

      /*
       * Indicate that this memory structure no longer has space allocated for
       * element results variables.
       */

      x->state ^= EXODB_STATE_ELVA;
    }

  return;
}

void free_exo_gv(Exo_DB *x)
{
  int i;

  if ( ! ( x->state & EXODB_STATE_GBVA ) )
    {
      return;
      /*      EH(-1, "Can't free what was never in chains.");*/
    }

  safe_free(x->gv_time_indeces);

  for ( i=0; i<x->num_gv_time_indeces; i++)
    {
      safe_free(x->gv[i]);
    }

  safe_free(x->gv);

  /*
   * Indicate that this memory structure no longer has space allocated for
   * global results variables.
   */

  x->state ^= EXODB_STATE_GBVA;

  return;
}

/*
 * Will free up exo->nv[][][]  and will also free up the auxiliary
 * arrays nv_indeces[] and nv_time_indeces[] if they have nonzero
 * lengths indicated by num_... variables.
 */

void
free_exo_nv(Exo_DB *x)
{
#if 0
  int i;
  int j;
#endif

  if ( x->state & EXODB_STATE_NDIA )
    {
#ifdef DEBUG
      fprintf(stderr, "in free_exo_nv with NDIA, nv_i=0x%x and nv_ti=0x%x\n",
	      (int)x->nv_indeces, (int)x->nv_time_indeces);
#endif
      if ( x->num_nv_indeces > 0 )
	{
	  safe_free(x->nv_indeces);
	}
      if ( x->num_nv_time_indeces > 0 )
	{
	  safe_free(x->nv_time_indeces);
	}
      x->state &= ~EXODB_STATE_NDIA; /* Flag that nv indeces are freed. */
    }

  if ( ! ( x->state & EXODB_STATE_NDVA ) )
    {
      return;
    }

  /*
   * New way to free this array since it's been allocated in large chunklets.
   */

  safe_free(x->nv[0][0]);
  safe_free(x->nv[0]);

#if 0
  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  safe_free(x->nv[i][j]);
	}
      safe_free(x->nv[i]);
    }
#endif

  safe_free(x->nv);

  /*
   * Indicate that this memory structure no longer has space allocated for
   * nodal results variables.
   */

  x->state ^= EXODB_STATE_NDVA;

  x->num_nv_indeces      = 0;
  x->num_nv_time_indeces = 0;

  return;
}

void
alloc_exo_ev(Exo_DB *x,
	     int num_timeplanes)
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

void
alloc_exo_gv(Exo_DB *x,
	     int num_timeplanes)
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


void 
alloc_init_exo_nv_indeces(Exo_DB *exo)
{
  int i;

  if ( exo->state & EXODB_STATE_NDIA )
    {
      EH(-1, "Do not allocate what has already been given!");
    }

  if ( exo->num_nv_time_indeces < 1 )
    {
      EH(-1, "Do not allocate nodal var indeces for no time steps!");
    }

  if ( exo->num_nv_indeces < 1 )
    {
      EH(-1, "Do not allocate nodal var indeces for no nodal vars!");
    }

  exo->nv_indeces = (int *) smalloc(exo->num_nv_indeces*sizeof(int));

  exo->nv_time_indeces = (int *) smalloc(exo->num_nv_time_indeces*sizeof(int));

  for ( i=0; i<exo->num_nv_time_indeces; i++) 
    {
      exo->nv_time_indeces[i] = i+1; /* fortran roolz */
    }   

  for ( i=0; i<exo->num_nv_indeces; i++) 
    {
      exo->nv_indeces[i] = i+1; /* fortran roolz */
    }   

  exo->state |= EXODB_STATE_NDIA;  

  return;
}



void
alloc_exo_nv(Exo_DB *x)
{
  int i;
  int j;

  if ( x->state & EXODB_STATE_NDVA )
    {
      EH(-1, "Please free before allocating...");
    }

  if ( x->num_nv_time_indeces < 1 )
    {
      fprintf(stderr, "x->num_nv_time_indeces = %d\n", x->num_nv_time_indeces);
      EH(-1, "Allocating space for EXODUS II nodal vars for <1 timeplane!");
    }

  if ( x->num_nv_indeces < 1 )
    {
      EH(-1, "Allocating space for EXODUS II nodal vars for <1 variable!");
    }

  x->nv       = (dbl ***) smalloc(x->num_nv_time_indeces*sizeof(dbl **));

#if 1
  /*
   * Allocate the mother of all arrays. For niceness, allocate in
   * big linear contiguous chunks so that strides are constant.
   *
   * The contiguousness helps performance a little and is needed
   * whenever netCDF is used to read multidimensional arrays (not often,
   * and not here, since nodal variables are read one shot at a time.)
   */

  x->nv[0]    = (dbl **)  smalloc(x->num_nv_time_indeces*x->num_nv_indeces*
			      sizeof(dbl *));

  x->nv[0][0] = (dbl *)   smalloc(x->num_nv_time_indeces * x->num_nv_indeces *
				  MAX(1,x->num_nodes) *sizeof(dbl));

  /*
   * The 1d pointers need to get set.
   */

  for ( i=1; i<x->num_nv_time_indeces; i++)
    {
      x->nv[i] = x->nv[0] + i * x->num_nv_indeces;
    }


  /*
   * A 2D array of pointers needs to be set.
   */
  
  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  x->nv[i][j] = x->nv[0][0] + j * x->num_nodes + 
	    i * x->num_nodes * x->num_nv_indeces; 
	}
    }

#endif
#if 0

  /*
   * The old easy way...
   */

  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      x->nv[i] = (dbl **) smalloc(x->num_nv_indeces*sizeof(dbl *));
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  x->nv[i][j] = (dbl *) smalloc(x->num_nodes*sizeof(dbl));
	}
    }


  for ( i=0; i<x->num_nv_time_indeces; i++)
    {
      for ( j=0; j<x->num_nv_indeces; j++)
	{
	  for ( k=0; k<x->num_nodes; k++)
	    {
	      x->nv[i][j][k] = 1e3*(double)i + (double)j + 1e-6*(double)k;
	    }
	}
    }

#endif

  x->state |= EXODB_STATE_NDVA;
  return;
}
