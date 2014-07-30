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

/* wr_exo --	open & write an EXODUS II finite element database. 
 *
 * wr_mesh_exo()
 *
 *	x -- a pointer to a structure containing the information from
 *	     an EXODUS IIv2 finite element database. A detailed description
 *	     of this structure may be found in "exo_struct.h"
 *
 *	filename -- where to write this mesh to
 *
 *	verbosity -- a integer that determines the level of chatty output
 *		     Larger values produce more output.
 *
 * Return values:
 *	integer -- A value of zero is returned upon normal completion of the
 *		   routine. A value of -1 is returned if abnormal conditions
 *		   were encountered.
 *
 * Notes:
 *
 *	1. wr_exo() does inquiries before attempting to write stuff out.
 *		    items need to exist before they can be written.
 *
 *	2. Assume that various memory allocation has already been done.
 *
 *	3. The emphasis is on reading mesh data and other preliminary
 *	   information. Little provision is made for reading results data
 *	   or history data from the file.
 *
 *
 * Created: 1997/05/09 14:12 MDT pasacki@sandia.gov
 *
 * Revised: 1998/01/26 13:34 MST pasacki@sandia.gov
 */

#define _WR_EXO_C

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


static Spfrtn sr=0;			/* sprintf() return type */
static char err_msg[MAX_CHAR_ERR_MSG];

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

int 
wr_mesh_exo(Exo_DB *x,		/* def'd in exo_struct.h */
	    char *filename,	/* where to write */
	    int verbosity)	/* how much to tell while writing */
{
  int i;
  int status=0;
  static char yo[] = "wr_mesh_exo";

  dbl dummy=0;

  if ( verbosity > 0 )
    {
      fprintf(stderr, "%s() begins.\n", yo);
    }

  /*
   * Mesh data is so fundamental that we'll create the file with clobber,
   * obliterating any existing file of the same name. That is, preserving
   * other data in an EXODUS II file while writing onto it new mesh information
   * is deemed too extraordinary. If mesh information is written, it causes
   * all information in the file to be superseded.
   */

  x->mode          = EX_WRITE;
  x->comp_wordsize = SZ_DBL;
  x->io_wordsize   = SZ_DBL;
  x->cmode         = EX_CLOBBER;

#ifdef DEBUG
  fprintf(stderr, "%s: ex_open with:\n", yo);
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", x->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", x->io_wordsize);
#endif

  x->exoid = ex_create(filename, x->cmode, &x->comp_wordsize, 
		       &x->io_wordsize);
  EH(x->exoid, "ex_create");
      
  if ( verbosity > 1 )
    {
      fprintf(stderr, "ex_open/create() rtn = %d\n", x->exoid);
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
      fprintf(stderr, "ex_put_init() call...\n");
    }

  status = ex_put_init(x->exoid,
		       x->title, 
		       x->num_dim, 
		       x->num_nodes,
		       x->num_elems, 
		       x->num_elem_blocks, 
		       x->num_node_sets,
		       x->num_side_sets);

  EH(status, "ex_put_init");

  if ( verbosity > 0 )
    {
      fprintf(stderr, "\tx->title           = \"%s\"\n", x->title);
      fprintf(stderr, "\tx->num_nodes       = %d\n", x->num_nodes);
      fprintf(stderr, "\tx->num_elems       = %d\n", x->num_elems);
      fprintf(stderr, "\tx->num_elem_blocks = %d\n", x->num_elem_blocks);
      fprintf(stderr, "\tx->num_node_sets   = %d\n", x->num_node_sets);
      fprintf(stderr, "\tx->num_side_sets   = %d\n", x->num_side_sets);
    }

  if ( verbosity > 1 )
    {
      fprintf(stderr, "\tx->num_qa_rec      = %d\n", x->num_qa_rec);
    }      

  if ( x->num_qa_rec > 0 )
    {
      status = ex_put_qa(x->exoid, x->num_qa_rec, x->qa_record);
      EH(status, "ex_put_qa");
    }

  if ( x->num_info > 0 )
    {
      status = ex_put_info(x->exoid, x->num_info, x->info);
      EH(status, "ex_put_info");
    }

  if ( verbosity > 0 )
    {
      fprintf(stderr, "ex_put_coord()...\n");
    }

  if ( x->num_dim < 3 )
    {
      x->z_coord = &dummy;
    }

  if ( x->num_dim < 2 )
    {
      x->y_coord = &dummy;
    }

  if ( x->num_dim < 1 )
    {
      x->x_coord = &dummy;
    }

  status = ex_put_coord(x->exoid, x->x_coord, x->y_coord, x->z_coord);
  EH(status, "ex_put_coord");

  status = ex_put_coord_names(x->exoid, x->coord_names);
  EH(status, "ex_get_coord_names");

  if ( x->num_nodes > 0 )
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_node_num_map()...\n");
	}
      if ( x->node_map_exists )
	{
	  status = ex_put_node_num_map(x->exoid, x->node_map);
	  EH(status, "ex_put_node_num_map");
	}
    }

  if ( x->num_elems > 0 )
    {
      
      if ( x->elem_map_exists )
	{
	  status = ex_put_elem_num_map(x->exoid, x->elem_map);	
	  EH(status, "ex_put_elem_num_map");
	}

      if ( x->elem_order_map_exists )
	{
	  status = ex_put_map(x->exoid, x->elem_order_map);
	  EH(status, "ex_put_map");
	}
    }

  /*
   * ELEMENT BLOCKS...
   */

  if ( x->num_elem_blocks > 0 )
    {
      for ( i=0; i<x->num_elem_blocks; i++)
	{
	  if ( verbosity > 0 )
	    {
	      fprintf(stderr, "ex_put_elem_block()...\n");
	    }
	  status = ex_put_elem_block(x->exoid, 
				     x->eb_id[i], 
				     x->eb_elem_type[i],
				     x->eb_num_elems[i],
				     x->eb_num_nodes_per_elem[i],
				     x->eb_num_attr[i]);
	  EH(status, "ex_put_elem_blocks");

	  if ( (x->eb_num_elems[i] * x->eb_num_nodes_per_elem[i]) > 0 )
	    {
	      status = ex_put_elem_conn(x->exoid, 
					x->eb_id[i], 
					x->eb_conn[i]);
	      EH(status, "ex_put_elem_conn");
	    }

	  if ( (x->eb_num_elems[i] * x->eb_num_attr[i]) > 0 )
	    {
	      status = ex_put_elem_attr(x->exoid, x->eb_id[i], x->eb_attr[i]);
	      EH(status, "ex_put_elem_attr");
	    }
	}
    }

  /*
   * NODE SETS...
   */

  if ( x->num_node_sets > 0 )
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_concat_node_sets()...\n");
	}
      status = ex_put_concat_node_sets(x->exoid, x->ns_id, 
				       x->ns_num_nodes,
				       x->ns_num_distfacts,
				       x->ns_node_index,
				       x->ns_distfact_index,
				       x->ns_node_list,
				       x->ns_distfact_list);
      EH(status, "ex_put_concat_node_sets");
    }


  /*
   * SIDE SETS...
   */

  if ( x->num_side_sets > 0 ) 
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_concat_side_sets()...\n");
	}
      status = ex_put_concat_side_sets(x->exoid, x->ss_id, 
				       x->ss_num_sides,
				       x->ss_num_distfacts,
				       x->ss_elem_index,
				       x->ss_distfact_index,
				       x->ss_elem_list,
				       x->ss_side_list,
				       x->ss_distfact_list);
      EH(status, "ex_put_concat_side_sets");
    }


  /*
   * PROPERTIES...
   */

  /*
   * Node sets...
   */

  if ( x->ns_num_props > 1 ) 
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_prop_names(nodesets) [%d]...\n", 
		  x->ns_num_props);
	}
      status = ex_put_prop_names(x->exoid, EX_NODE_SET, x->ns_num_props - 1,
				 &(x->ns_prop_name[1]) );
      EH(status, "ex_put_prop_names(EX_NODE_SET)");

      /* 
       * the following loop begins at 1 so as avoid writing
       * the first "ID" node set property table
       * This automatically added by ex_put_prop_array
       * as the first property table written to all exodus files
       * Consequently, if we were to write the "ID" table out
       * here it would continually be replicated as the file
       * is repeatedly rewritten
       */

      for ( i=1; i<x->ns_num_props; i++)
	{
	  status = ex_put_prop_array(x->exoid, EX_NODE_SET, 
				     x->ns_prop_name[i],
				     x->ns_prop[i]);
	  EH(status, "ex_put_prop_array(EX_NODE_SET)");

	}
    }
      
  /*
   * Side sets...
   */

  if ( x->ss_num_props > 1 ) 
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_prop_names(sidesets)...\n");
	}
      status = ex_put_prop_names(x->exoid, EX_SIDE_SET, x->ss_num_props - 1,
				 &(x->ss_prop_name[1]));
      EH(status, "ex_get_prop_names(EX_SIDE_SET)");

      /*
       * See node set prop table write for comment 
       */

      for ( i=1; i<x->ss_num_props; i++)
	{
	      
	  status = ex_put_prop_array(x->exoid, EX_SIDE_SET, 
				     x->ss_prop_name[i],
				     x->ss_prop[i]);
	  EH(status, "ex_put_prop_array(EX_SIDE_SET)");
	}
    }
      
  /*
   * Element blocks...
   */

  if ( x->eb_num_props > 1 )
    {
      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_prop_names(elemblocks)...\n");
	}
      status = ex_put_prop_names(x->exoid, EX_ELEM_BLOCK, x->eb_num_props - 1,
				 &(x->eb_prop_name[1]));
      EH(status, "ex_get_prop_names(EX_ELEM_BLOCK)");

      /*
       * See node set prop table write for comment 
       */

      for ( i=1; i<x->eb_num_props; i++)
	{
	  if ( strcmp( x->eb_prop_name[i], "ID") != 0 )
	    {

	      status = ex_put_prop_array(x->exoid, EX_ELEM_BLOCK, 
					 x->eb_prop_name[i],
					 x->eb_prop[i]);
	      EH(status, "ex_put_prop_array(EX_ELEM_BLOCK)");
	    }
	}
    }
      

  status = ex_close(x->exoid);
  EH(status, "ex_close()");
  
  return(status);
}


/* wr_resetup_exo() -- open/write/close EXODUS II db for results names
 *
 * Created: 1998/01/26 14:06 MST pasacki@sandia.gov
 *
 * Revised: 
 */

void 
wr_resetup_exo(Exo_DB *exo,
	       char *filename,
	       int verbosity)
{
  int error;
  int i;
  int status;

  /*
   * This file must already exist.
   */

  exo->cmode = EX_WRITE;

#ifdef DEBUG
  fprintf(stderr, "%s: begins\n", yo);
#endif

  exo->io_wordsize   = 0;	/* i.e., query */
  exo->comp_wordsize = sizeof(dbl);
  exo->exoid         = ex_open(filename, exo->cmode, &exo->comp_wordsize, 
			       &exo->io_wordsize, &exo->version);

#ifdef DEBUG
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  /*
   * Results setup...
   */

  if ( exo->num_glob_vars > 0 )
    {
      status = ex_put_var_param(exo->exoid, "g", exo->num_glob_vars);
      EH(status, "ex_put_var_param(g)");
      status = ex_put_var_names(exo->exoid, "g", exo->num_glob_vars,
				exo->glob_var_names);
      EH(status, "ex_put_var_names(g)");
    }

  if ( exo->num_elem_vars > 0 )
    {
      status = ex_put_var_param(exo->exoid, "e", exo->num_elem_vars);
      EH(status, "ex_put_var_param(e)");
      status = ex_put_var_names(exo->exoid, "e", 
				exo->num_elem_vars,
				exo->elem_var_names);
      EH(status, "ex_put_var_names(e)");
      status = ex_put_elem_var_tab(exo->exoid, 
				   exo->num_elem_blocks,
				   exo->num_elem_vars, 
				   exo->elem_var_tab);
      EH(status, "ex_put_elem_var_tab");
    }

  if ( exo->num_node_vars > 0 )
    {
      status = ex_put_var_param(exo->exoid, "n", exo->num_node_vars);
      EH(status, "ex_put_var_param(n)");
      status = ex_put_var_names(exo->exoid, "n", 
				exo->num_node_vars,
				exo->node_var_names);
      EH(status, "ex_put_var_names(n)");
    }

  if ( exo->num_times > 0 )
    {
      for ( i=0; i<exo->num_times; i++)
	{
	  status = ex_put_time(exo->exoid, i+1, &(exo->time_vals[i]));
	  EH(status, "ex_put_times");
	}
    }

  error      = ex_close(exo->exoid);
  if ( error != 0 ) exit(2);
  return;
}


void 
wr_result_exo(Exo_DB *exo,
	      char *filename,
	      int verbosity)
{
  int i;
  int index;
  int j;
  int k;
  int status;
  int time_index;

  /*
   * This file should already exist.
   */

  exo->cmode = EX_WRITE;

#ifdef DEBUG
  fprintf(stderr, "%s: begins\n", yo);
#endif

  exo->io_wordsize   = 0;	/* i.e., query */
  exo->comp_wordsize = sizeof(dbl);
  exo->exoid         = ex_open(filename, 
			       exo->cmode, 
			       &exo->comp_wordsize, 
			       &exo->io_wordsize, 
			       &exo->version);

#ifdef DEBUG
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  /*
   * Element variable truth table and values at ONE TIME ONLY.
   */

  if ( exo->num_elem_vars > 0 )
    {

#ifdef DEBUG
      fprintf(stderr, "\t\tneb         = %d\n", exo->num_elem_blocks);
      fprintf(stderr, "\t\tnev         = %d\n", exo->num_elem_vars);
      fprintf(stderr, "\t\tevt:        =   \n");
      for ( i=0; i<exo->num_elem_blocks; i++)
	{
	  for ( j=0; j<exo->num_elem_vars; j++)
	    {
	      fprintf(stderr, "block index %d, elem var index %d is %d\n",
		      i, j, exo->elem_var_tab[i*(exo->num_elem_vars)+j]);
	    }
	}
#endif

      /*
       * This has already been done.
       *

      status = ex_put_elem_var_tab(exo->exoid, 
				   exo->num_elem_blocks,
				   exo->num_elem_vars, 
				   exo->elem_var_tab);
      EH(status, "ex_put_elem_var_tab");
      */

      for ( i=0; i<exo->num_ev_time_indeces; i++)
	{
	  time_index = exo->ev_time_indeces[i];
	  
	  for ( j=0; j<exo->num_elem_blocks; j++)
	    {
	      for ( k=0; k<exo->num_elem_vars; k++)
		{
		  index = j * exo->num_elem_vars + k;
		  
		  if ( exo->elem_var_tab[index] != 0 )
		    {
		      status = ex_put_elem_var(exo->exoid, time_index, k+1,
					       exo->eb_id[j], 
					       exo->eb_num_elems[j],
					       &(exo->ev[i][index][0]));
		      if ( status < 0 )
			{
			  sr = sprintf(err_msg, 
				       "ex_put_elem_var() bad rtn: time %d, elemvar %d, EB ID %d",
				       time_index, k+1, exo->eb_id[j]);
			  EH(-1, err_msg);
			  EH(sr, err_msg);
			}
		    }
		}
	    }
	}
    }

  /*
   * Put nodal variable values at last time step...
   */

  if ( exo->num_node_vars > 0 )
    {
      for ( i=0; i<exo->num_nv_time_indeces; i++)
	{
	  time_index = exo->nv_time_indeces[i];

	  for ( j = 0; j < exo->num_nv_indeces; j++)
	    {
	      status = ex_put_nodal_var(exo->exoid, time_index, 
					exo->nv_indeces[j],
					exo->num_nodes, 
					&(exo->nv[i][j][0]));
	      EH(status, "ex_put_nodal_var");
	    }
	}
    }

  status = ex_close(exo->exoid);
  EH(status, "ex_close()");
  
  return;
}

