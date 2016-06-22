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
 
/* wr_exo -- open, write EXODUS II finite element mesh, close.
 *
 * Input:
 *	x -- a pointer to a structure containing the information from
 *	     an EXODUS IIv2 finite element database. A detailed description
 *	     of this structure may be found in "exo_struct.h"
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
 *	1. wr_exo() does inquiries before attempting to write stuff out.
 *		    items need to exist before they can be written.
 *
 *	2. Assume that various memory allocation has already been done.
 *
 *	3. Based on code for the more general wr_exo() routine, but
 *	   chopped down to do just the mesh.
 *
 * Created: 1997/08/04 10:02 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#include <sys/types.h>
#include <time.h>
#include <errno.h>
#include <sys/utsname.h>
#include <pwd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>		/* for getuid() */

#ifndef lint
#ifdef USE_RCSID
static char rcsid[] = "$Id: wr_exo.c,v 5.3 2008-05-02 19:07:58 hkmoffa Exp $";
#endif
#endif

static int has_been_called=0;

#include "std.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "rf_bc.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_allo.h"
#include "rf_mp.h"		/* are we serial or parallel? */
#include "rf_io_structs.h"	/* for Results_Description */

#define _WR_EXO_C
#include "goma.h"

extern char **Argv;		/* global shadow of argv, def'd in main.c */

static Spfrtn sr;		/* sprintf() return type, whatever it is. */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

int 
wr_mesh_exo(Exo_DB *x,		/* def'd in exo_struct.h */
	    char *filename,	/* where to write */
	    int verbosity)	/* how much to tell while writing */
{
#ifdef DEBUG
  char *yo = "wr_nodal_results_exo: ";
#endif
  int i;
  int status=0;

  /* 
   * This is a sad and pathetic little hack intended only for short term
   * use.  If its been longer than a year since 9/3/99 it should be replaced
   * with better code.  TABAER */

  dbl dummy=0;

  if ( verbosity > 0 )
    {
      fprintf(stderr, "wr_mesh_exo() begins.\n");
    }

  /*
   * Mesh data is so fundamental that we'll create the file with clobber,
   * obliterating any existing file of the same name. That is, preserving
   * other data in an EXODUS II file while writing onto it new mesh information
   * is deemed too extraordinary. If mesh information is written, it causes
   * all information in the file to be superseded.
   */

  x->io_wordsize = 8;

#ifdef DEBUG
  fprintf(stderr, "%s: ex_open with:\n", yo);
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", x->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", x->io_wordsize);
#endif

  x->cmode = EX_CLOBBER;
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
      fprintf(stderr, "\tx->num_info        = %d\n", x->num_info);
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

	  if ( (x->eb_num_elems[i]*x->eb_num_attr[i]) > 0 )
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
   * EXODUS II will write out the default table of one property called 
   * the "ID" for NS, SS, and EBs. Unless you actually have more properties
   * you want dumped, then we'll not write these out.
   */

  /*
   * Well, the damage s done. Very old EXODUS II data sets have spuriously
   * compounded ID properties already. 
   */

  /*
   * Node sets...
   */

  if ( x->ns_num_props > 1 ) 
    {

      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_prop_names(nodesets)...\n");
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
	  if( strcmp( x->ns_prop_name[i] , "ID" ) !=0 )
	    {
	      status = ex_put_prop_array(x->exoid, EX_NODE_SET, 
					 x->ns_prop_name[i],
					 x->ns_prop[i]);
	      EH(status, "ex_put_prop_array(EX_NODE_SET)");
	    }
	}
      
    }
      
  /*
   * Side sets...
   */

  if ( x->ss_num_props > 1 ) 
    {

      /*
       * Only write these out if the second property is not the same ole
       * "ID" like the first one...
       */

      if ( verbosity > 0 )
	{
	  fprintf(stderr, "ex_put_prop_names(sidesets)...\n");
	}
      status = ex_put_prop_names(x->exoid, EX_SIDE_SET, x->ss_num_props - 1,
				 &(x->ss_prop_name[1]));
      EH(status, "ex_get_prop_names(EX_SIDE_SET)");
	  
      for ( i=1; i<x->ss_num_props; i++)
	{
	  if( strcmp( x->ss_prop_name[i] , "ID" ) !=0 )
	    {

	      status = ex_put_prop_array(x->exoid, EX_SIDE_SET, 
					 x->ss_prop_name[i],
					 x->ss_prop[i]);
	      EH(status, "ex_put_prop_array(EX_SIDE_SET)");
	    }
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
				 &(x->eb_prop_name[1]) );
      EH(status, "ex_put_prop_names(EX_ELEM_BLOCK)");

      for ( i=1; i<x->eb_num_props; i++)
	{
	  if( strcmp( x->ss_prop_name[i] , "ID" ) !=0 )
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

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int file_is_readable(const char *filename)

    /***********************************************************************
     *
     * file_is_readable()
     *
     *  Attempts to determine if a file exists and is readable by this
     *  user
     *
     *  Input
     * -------
     *    filename = Relative pathname to the file.
     *
     *  Return
     * --------
     *    True (1) = File exists and is readable
     *    false(0) = File doesn't exist or is not readable.
     *
     *  (note: fopen will set the global variable, errno, as well. Therefore,
     *         there is more specific information if needed)
     ***********************************************************************/
{
  FILE *fp;
  if (filename == NULL) return 0;
  fp = fopen(filename, "r");
  if (fp == NULL) return 0;
  (void) fclose(fp);
  return 1;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* wr_result_prelim_exo() -- open, setup EXODUS II db for results, close
 *
 * One time initialization that writes out how many nodal variables to expect
 * in the results and what their names are. This information comes in via
 * the Results_Description.
 *
 * This routine is based on the setup_exoII_results() of long ago, but modified
 * to use the new exo_struct.
 *
 * The filename written to is "exo->path", so set it to what you need and
 * remember to set it back for the next person.
 *
 * in:
 * ---
 *	rd		pointer to a "Results_Description" structure that 
 * 			describes the results that are about to be dumped.
 *			for more information on this structure, see the
 *			file "rf_io.h"
 *
 *	exo		ptr to Exodus II finite element database that has
 *			been read in. Basically, expect that mesh data is
 *			known at this point, even if the results have not
 *			yet been computed.
 * out:
 * ----
 *	(none)
 *
 * Revised: 1997/08/26 14:06 MDT pasacki@sandia.gov
 */

void
wr_result_prelim_exo(struct Results_Description *rd, 
		     Exo_DB *exo,
		     char *filename,
		     double ***gvec_elem )
{
  int i, error;
  int filename_exists;		/* boolean */

  int	num_vars;		/* number of var_type variables to be written */
  char	*var_names[MAX_NNV];	/* NOTE: this array must be sized to be the 
				   greatest value of MAX_NNV, MAX_NEV, MAX_NGV,
				   etc. It is reused for each of these variable
				   types. Currently (8/13/98), MAX_NNV is the 
				   largest value @ 100 of them all - RRL   */

  char	*gvar_names[MAX_NGV];	
				/* array containing num_vars variable names */

#ifdef DEBUG  
  static char *yo="wr_result_prelim_exo: ";
#endif
  /*
   * We don't support history variables.
   */

  if ( rd->nhv > 0 )
    {
      EH(-1, "Not prepared to write history variables.");
    }

  if ( filename == NULL )
    {
      EH(-1, "No file specified to write EXODUS II info.");
    }

  /*
   *  Figure out whether the file exists and is readable by this
   *  user.
   */
  filename_exists = file_is_readable(filename);

  /*
   * If this is an entirely new file, then we'll take the opportunity to
   * write out 8 byte floating point data. If it's an old file, then
   * we'll query it by passing a zero and respect it's wishes.
   *
   * The I/O wordsize is really immaterial to the memory resident structure
   * and is often obsolete when reading and writing memory contents to and
   * from different EXODUS II files. Same applies for the filename.
   */
  if ( filename_exists )
    {
      exo->io_wordsize = 0;	/* query the file */
#ifdef DEBUG
      fprintf(stderr, "%s: \"%s\" already exists\n", yo, filename);
#endif
    }
  else
    {
      exo->io_wordsize = 8;
#ifdef DEBUG
      fprintf(stderr, "%s: \"%s\" does not exist, going for double\n", yo, 
	      filename);
#endif
    }

  exo->cmode = EX_WRITE;

#ifdef DEBUG
  fprintf(stderr, "%s: ex_open with:\n", yo);
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  if ( filename_exists )
    {
#ifdef DEBUG
      fprintf(stderr, "P_%d at barrier < ex_open in wr_result_prelim_exo()\n", 
	      ProcID);
#endif

      exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, 
			   &exo->io_wordsize, &exo->version);
      EH(exo->exoid, "ex_open");
    }
  else
    {
      exo->exoid = ex_create(filename, exo->cmode, &exo->comp_wordsize, 
			   &exo->io_wordsize);
      EH(exo->exoid, "ex_create");
    }


  /*
   * Analysis Results
   * -------- -------
   */

  /* --------------------- Global Variables -------------------------- */
  if ( rd->ngv > 0 )
    {
      num_vars = rd->ngv;
      error = ex_put_var_param(exo->exoid, "g", num_vars);
      EH(error, "ex_put_var_param(g)");
      for ( i=0; i<rd->ngv; i++)
	{
	 gvar_names[i] = rd->gvname[i];
	}
      error = ex_put_var_names(exo->exoid, "g", num_vars, gvar_names);
      EH(error, "ex_put_var_names(g)");
    }
      
  /* -------------------- Element Variables -------------------------- */
  if ( rd->nev > 0 )
    {
      num_vars = rd->nev;
      error = ex_put_var_param(exo->exoid, "e", num_vars);
      EH(error, "ex_put_var_param(e)");
      for ( i=0; i<rd->nev; i++)
	{
	  var_names[i] = rd->evname[i];
#ifdef DEBUG
	  printf("%s: elem var_name[%d] = %s\n", yo, i, var_names[i]);
#endif
	}
#ifdef DEBUG
      printf("%s: varnames loaded\n", yo);
#endif
      error = ex_put_var_names(exo->exoid, "e", num_vars, var_names);
      EH(error, "ex_put_var_names(e)");

      /* Create truth table at this time - saves mucho cycles later
         Also malloc the gvec_elem final dim. Easier to do right
         when the truth table is built. */

      create_truth_table(rd, exo, gvec_elem);
    }

  /* -------------------- Nodal Variables -------------------------- */
  if ( rd->nnv > 0 )
    {
      num_vars = rd->nnv;
      error = ex_put_var_param(exo->exoid, "n", num_vars);
      EH(error, "ex_put_var_param(n)");
      for ( i=0; i<rd->nnv; i++)
	{
	  var_names[i] = rd->nvname[i];
#ifdef DEBUG
	  printf("%s: nodal var_name[%d] = %s\n", yo, i, var_names[i]);
#endif
	}
#ifdef DEBUG
      printf("%s: varnames loaded\n", yo);
#endif
      error = ex_put_var_names(exo->exoid, "n", num_vars, var_names);
      EH(error, "ex_put_var_names(n)");
    }

  error = ex_close(exo->exoid);
  EH(error, "ex_close");

  return;
}
/* End of wr_result_prelim_exo() ------------------------------------------- */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void 
wr_nodal_result_exo(Exo_DB *exo, char *filename, double vector[],
		    int variable_index, int time_step, double time_value)

     /*****************************************************************
      * write_nodal_result_exo() 
      *     -- open/write/close EXODUS II db for 1 nodal var at one
      *        time step.
      *
      * The output EXODUS II database contains the original model
      * information with some minor QA and info additions, with new 
      * nodal value solution data written.
      * 
      ******************************************************************/
{
  char err_msg[MAX_CHAR_IN_INPUT];
  int error;
  exo->cmode = EX_WRITE;
  exo->io_wordsize = 0;		/* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, 
		       &exo->io_wordsize, &exo->version);
  if (exo->exoid < 0)
    {
      sr = sprintf(err_msg, 
		   "ex_open() = %d on \"%s\" failure @ step %d, time = %g",
		   exo->exoid, filename, time_step, time_value);
      EH(-1, err_msg);
    }
  error      = ex_put_time(exo->exoid, time_step, &time_value);
  EH(error, "ex_put_time");
  error      = ex_put_nodal_var(exo->exoid, time_step, variable_index, 
				exo->num_nodes, vector);
  EH(error, "ex_put_nodal_var");
  error      = ex_close(exo->exoid);
  return;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void 
wr_elem_result_exo(Exo_DB *exo, const char *filename, double ***vector, 
		   const int variable_index, const int time_step,
		   const double time_value, 
		   struct Results_Description *rd)
{
  int error, i;
  double local_time_value=time_value;
  /* static char *yo = "wr_elem_result_exo"; */

  /*
   * This file must already exist.
   */

  exo->cmode = EX_WRITE;

  exo->io_wordsize = 0;		/* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, 
		       &exo->io_wordsize, &exo->version);
  EH(exo->exoid, "ex_open");

#ifdef DEBUG
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  error = ex_put_time (exo->exoid, time_step, &local_time_value );
  EH(error, "ex_put_time");

  /* If the truth table has NOT been set up, this will be really slow... */

  for (i = 0; i < exo->num_elem_blocks; i++) {
    if (exo->elem_var_tab_exists == TRUE) {
      /* Only write out vals if this variable exists for the block */
      if (exo->elem_var_tab[i*rd->nev + variable_index] == 1) {
	error = ex_put_elem_var(exo->exoid, time_step, variable_index+1,
				exo->eb_id[i], exo->eb_num_elems[i],
				vector[i][variable_index]);
	EH(error, "ex_put_elem_var");
      }
    }
    else {
      /* write it anyway (not really recommended from a performance viewpoint) */
      error      = ex_put_elem_var ( exo->exoid,
				     time_step, 
				     variable_index+1, /* Convert to 1 based for 
							  exodus */
				     exo->eb_id[i],
				     exo->eb_num_elems[i],
				     vector[i][variable_index] );
      EH(error, "ex_put_elem_var");
    }
  }

  error      = ex_close ( exo->exoid );
  EH(error, "ex_close");

  return;
}


void
wr_global_result_exo( Exo_DB *exo, 
		      const char *filename, 
		      const int time_step, 
		      const int ngv, 
		      double u[] )
{
     /*****************************************************************
      * write_global_result_exo() 
      *     -- open/write/close EXODUS II db for all global values
      *
      * The output EXODUS II database contains the original model
      * information with some minor QA and info additions, with new 
      * global data written.
      * 
      ******************************************************************/
  int error;


  /* 
   * This capability is deactivated for parallel processing.
   * brkfix doesn't support global variables, when this
   * changes this restriction should be removed. TAB 3/2002 */

  if( u == NULL ) return ; /* Do nothing if this is NULL */

  exo->cmode = EX_WRITE;
  exo->io_wordsize = 0;		/* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, 
		       &exo->io_wordsize, &exo->version);
  if (exo->exoid < 0) {
    EH(-1,"wr_nodal_result_exo: could not open the output file");
  }

  error = ex_put_glob_vars ( exo->exoid, time_step, ngv, u );

  EH(error, "ex_put_glob_vars");

  error = ex_close( exo->exoid);


  return;
}
/* End of write_exoII_results ----------------------------------------- */



/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * create_truth_table() -- Test each proposed element against each element
 * block for existance and populate the truth table. This saves file
 * rewrites later (a netcdf phenom) when each variable is written sequentially
 * later.
 *
 * Created: 1998/08/19 08:17 MDT rrlober@sandia.gov
 */
void
create_truth_table(struct Results_Description *rd, Exo_DB *exo,
		   double ***gvec_elem )
{
  char err_msg[MAX_CHAR_IN_INPUT], if_ev;
  int   i, j, eb_indx, ev_indx, mat_num, error, check, iii;
  int	tev, found_match, ip_total;
  ELEM_BLK_STRUCT *eb_ptr;

  static const char yo[] = "create_truth_table";

  tev = 0;
  i = 0;
  exo->elem_var_tab = (int *) smalloc( (exo->num_elem_blocks*rd->nev)*sizeof(int));
  exo->truth_table_existance_key = (int *) smalloc( (V_LAST - V_FIRST)*sizeof(int));

  for ( i = 0; i < V_LAST - V_FIRST; i++ ) {
    exo->truth_table_existance_key[i] = 0;
  }

  /* This first cycle is a test to detect which potential elem variables 
     (treated as nodal in goma) exist in which block via tests below on
     whether the variables have been requested at all by the user and if
     requested, are of the appropriate integration order for conversion
     to an elem var.

     This is necessary since the array of the truth table cycles through
     the element var index fastest, and if a given element var is defined
     for one block but not another, the block in which it is undefined 
     will not know the difference between a defined variable that is 
     of the wrong interpolation order for this block, and a variable
     that is not defined at all for the problem. This first cycle scopes
     for these cases and sets up a temp array of all possible elem vars 
     model wide that must be treated (1 or 0) in the truth table. RRL */

  if_ev = FALSE;

  for ( eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++ ) {
    /* First test for all the potential elem vars from primary nodal vars
       for this block */
    mat_num = Matilda[eb_indx];
    if (mat_num < 0) {
      continue;
    }
    for ( j = V_FIRST; j < V_LAST; j++) {
      if ( pd_glob[mat_num]->v[j] != V_NOTHING ) {
	if ( pd_glob[mat_num]->i[j] == I_P0 ) {
	  if ( Num_Var_In_Type[j] > 1 ) {
	    fprintf(stderr,
		    "%s: Too many components in variable type for element variable %s (%s)\n",
		    yo,
		    Exo_Var_Names[j].name2,
		    Exo_Var_Names[j].name1 );
	    exit (-1);
	  }
	  if ( exo->truth_table_existance_key[j - V_FIRST] == 0 ) {
	    /* We just found a candidate for an element variable */
	    tev += Num_Var_In_Type[j];
	    exo->truth_table_existance_key[j - V_FIRST] = 1;
	  }
        }
      }
    }

    /* Now pick up all the post processing variables for this block 
       - yes, for now they must
       each be listed separately and painfully */
    if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[R_MOMENTUM1]) {
      tev++;
      if (ERROR_ZZ_VEL_ELSIZE  != -1) {
	tev++;	
      }
    }
    if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[R_ENERGY]) {
      tev++;
      if (ERROR_ZZ_Q_ELSIZE  != -1) {
	tev++;	
      }
    }
    check = 0;
    for ( i = 0; i < upd->Num_Mat; i++ ) {
      if( pd_glob[i]->MeshMotion == LAGRANGIAN ||
	  pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN) check = 1;
    }

    if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[R_MOMENTUM1] || check)) {
      tev++;
      if (ERROR_ZZ_P_ELSIZE  != -1) {
	tev++;	
      }
    }

    /* Finally pick up all of the element-level-storage continuation
     * variables, e.g. for saturation hysteresis function 
     */
    mat_num = Matilda[eb_indx];
    mp = mp_glob[mat_num];
    eb_ptr = Element_Blocks + eb_indx;
    ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
    if((mp->PorousMediaType == POROUS_UNSATURATED ||
	mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
	mp->PorousMediaType == POROUS_TWO_PHASE) &&
       mp->SaturationModel == TANH_HYST &&
       !if_ev)
      {
	for ( j = 0; j < ip_total; j++) {
	  if(SAT_CURVE_TYPE != -1) tev++; /*For Sat curve type */
	  if(CAP_PRESS_SWITCH != -1) tev++; /*For saturation switch */
	  if(SAT_QP_SWITCH != -1) tev++; /*for cap press switch point*/
	}
	if_ev = TRUE;
      }
  }

  /* Sanity check */
  if ( tev != rd->nev ) {
    sr = sprintf(err_msg, 
		 "%s: Elem var count mismatch: tev(%d)<>rd->nev(%d)!?",
		 yo, tev, rd->nev);
    EH(-1, err_msg);
    /*
    fprintf(stderr,
	    "%s: Disagreement over number of element variables\n",
	    yo );
    exit (-1);    
    */
  }
  
  /* Now do the real loop and populate the truth table */
  i = 0;
  for ( eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++ ) {
    /* First test for all the potential elem vars from primary nodal vars
       for this block */
    mat_num = Matilda[eb_indx];
    ev_indx = 0;
    for ( j = V_FIRST; j < V_LAST; j++) {
      found_match = FALSE;
      if ( pd_glob[mat_num]->v[j] != V_NOTHING ) {
	if ( pd_glob[mat_num]->i[j] == I_P0 ) {
	  if ( Num_Var_In_Type[j] > 1 ) {
	    fprintf(stderr,
		    "%s: Too many components in variable type for element variable %s (%s)\n",
		    yo,
		    Exo_Var_Names[j].name2,
		    Exo_Var_Names[j].name1 );
	    exit (-1);
	  }
	  /* We just found a candidate for an element variable */
	  exo->elem_var_tab[i++] = 1;
	  found_match = TRUE;
	  ev_indx++;
	  /* malloc the entry for this block by number of elems for this block 
	     but - only if the variable exists for this block! (by the truth table) */
	  if ( has_been_called == 0 ) {
	    /* NOTE: this final array dim is only to be malloc'd once; when a user 
	       is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	       and hence create_truth_table, which would realloc this dim of gvec_elem.
	       this test will prevent that. - RRL */
	    asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		   exo->eb_num_elems[eb_indx] );
	  }
	}
      }
      if ( found_match == FALSE && exo->truth_table_existance_key[j - V_FIRST] == 1 ) {
	exo->elem_var_tab[i++] = 0;
	ev_indx++;
      }
    }

    /* Now pick up all the post processing variables for this block 
       - yes, for now they must
       each be listed separately and painfully */

    if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[R_MOMENTUM1]) {
      exo->elem_var_tab[i++] = 1;   
      ev_indx++;
      /* malloc the entry for this block by number of elems for this block 
	 but - only if the variable exists for this block! (by the truth table) */
      if ( has_been_called == 0 ) {
	/* NOTE: this final array dim is only to be malloc'd once; when a user 
	   is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	   and hence create_truth_table, which would realloc this dim of gvec_elem.
	   this test will prevent that. - RRL */
	asdv ( &gvec_elem[eb_indx][ev_indx - 1],
	       exo->eb_num_elems[eb_indx] );
      }
      if (ERROR_ZZ_VEL_ELSIZE  != -1) {
	exo->elem_var_tab[i++] = 1;   
	ev_indx++;
	/* malloc the entry for this block by number of elems for this block 
	   but - only if the variable exists for this block! (by the truth table) */
	if ( has_been_called == 0 ) {
	  /* NOTE: this final array dim is only to be malloc'd once; when a user 
	     is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	     and hence create_truth_table, which would realloc this dim of gvec_elem.
	     this test will prevent that. - RRL */
	  asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		 exo->eb_num_elems[eb_indx] );
	}
      }
    }

    if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[R_ENERGY]) {
      exo->elem_var_tab[i++] = 1;   
      ev_indx++;
      /* malloc the entry for this block by number of elems for this block 
	 but - only if the variable exists for this block! (by the truth table) */
      if ( has_been_called == 0 ) {
	/* NOTE: this final array dim is only to be malloc'd once; when a user 
	   is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	   and hence create_truth_table, which would realloc this dim of gvec_elem.
	   this test will prevent that. - RRL */
	asdv ( &gvec_elem[eb_indx][ev_indx - 1],
	       exo->eb_num_elems[eb_indx] );
      }
      if (ERROR_ZZ_Q_ELSIZE  != -1) {
	exo->elem_var_tab[i++] = 1;   
	ev_indx++;
	/* malloc the entry for this block by number of elems for this block 
	   but - only if the variable exists for this block! (by the truth table) */
	if ( has_been_called == 0 ) {
	  /* NOTE: this final array dim is only to be malloc'd once; when a user 
	     is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	     and hence create_truth_table, which would realloc this dim of gvec_elem.
	     this test will prevent that. - RRL */
	  asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		 exo->eb_num_elems[eb_indx] );
	}
      }
    }

    check = 0;
    for ( iii = 0; iii < upd->Num_Mat; iii++ ) {
      if( pd_glob[iii]->MeshMotion == LAGRANGIAN ||
	  pd_glob[iii]->MeshMotion == DYNAMIC_LAGRANGIAN) check = 1;
    }

    if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[R_MOMENTUM1] || check)) {
      exo->elem_var_tab[i++] = 1;   
      ev_indx++;
      /* malloc the entry for this block by number of elems for this block 
	 but - only if the variable exists for this block! (by the truth table) */
      if ( has_been_called == 0 ) {
	/* NOTE: this final array dim is only to be malloc'd once; when a user 
	   is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	   and hence create_truth_table, which would realloc this dim of gvec_elem.
	   this test will prevent that. - RRL */
	asdv ( &gvec_elem[eb_indx][ev_indx - 1],
	       exo->eb_num_elems[eb_indx] );
      }
      if (ERROR_ZZ_P_ELSIZE  != -1) {
	exo->elem_var_tab[i++] = 1;   
	ev_indx++;
	/* malloc the entry for this block by number of elems for this block 
	   but - only if the variable exists for this block! (by the truth table) */
	if ( has_been_called == 0 ) {
	  /* NOTE: this final array dim is only to be malloc'd once; when a user 
	     is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
	     and hence create_truth_table, which would realloc this dim of gvec_elem.
	     this test will prevent that. - RRL */
	  asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		 exo->eb_num_elems[eb_indx] );
	  has_been_called++;
	}
      }
    }

    /*Now finally the saturation hysteresis variables */
    if(SAT_CURVE_TYPE != -1 || CAP_PRESS_SWITCH != -1 || SAT_QP_SWITCH != -1)
       {
	 eb_ptr = Element_Blocks + eb_indx;
	 ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
	 for(j=0; j < ip_total; j++)
	   {
	     /*Note that we will set these for all 3 var types because you
	      *will never see them individually. 
	      */
	     exo->elem_var_tab[i++] = 1;
	     ev_indx++;
	     if ( has_been_called == 0 ) 
	       {
		 asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		      exo->eb_num_elems[eb_indx] );
	       }
	     exo->elem_var_tab[i++] = 1;
	     ev_indx++;
	     if ( has_been_called == 0 ) 
	       {
		 asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		      exo->eb_num_elems[eb_indx] );
	       }
	     exo->elem_var_tab[i++] = 1;
	     ev_indx++;
	     if ( has_been_called == 0 ) 
	       {
		 asdv ( &gvec_elem[eb_indx][ev_indx - 1],
		      exo->eb_num_elems[eb_indx] );
	       }
	   }
       }
  }

  /* write out table */
  error      = ex_put_elem_var_tab ( exo->exoid,
				     exo->num_elem_blocks,
				     rd->nev,
				     exo->elem_var_tab );
  EH(error, "ex_put_elem_var_tab");

  /* Now set truth table exists flag */
  exo->elem_var_tab_exists = TRUE;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * add_qa_stamp() -- if space exists, add a goma QA record to the EXODUS II db.
 *
 * Notes:    [1] A new QA record is created that is slightly larger that the
 *		 old one. The old information is transcribed and a new record
 *		 is added that reflects this analysis.
 *		
 *	     [2] New records are not added if there are at least 4 QA records
 *		 already accumulated. The threshhold varies according to the
 *		 C preprocessor symbol MAX_QA.
 *		
 *		c.f. SAND83-0905
 *
 * Created: 1997/08/04 08:17 MDT pasacki@sandia.gov
 * 
 * Revised: 2000/02/07 12:40 MST pasacki@sandia.gov
 */

void
add_qa_stamp(Exo_DB *exo)
{
  int i;
  int j;
  int k;
  int n;

  /*  char *new[MAX_QA][4];*/

  QA_Record *Q;

  n = exo->num_qa_rec;

  if ( n > MAX_QA - 1 )
    {
      return;
    }

  exo->num_qa_rec++;

  Q = (QA_Record *) smalloc( (exo->num_qa_rec)*sizeof(QA_Record));

  for (i=0; i<exo->num_qa_rec; i++)
    {
      for (j=0; j<4; j++)
	{
	  Q[i][j] = (char *) smalloc(LEN_QA_RECORD*sizeof(char));

	  /*
	   * Initialize to null terminators...
	   */

	  for (k=0; k<LEN_QA_RECORD; k++)
	    {
	      Q[i][j][k] = '\0';
	    }
	}
    }

  /*
   * Copy any old stuff into the new area.
   */

  for ( i=0; i<exo->num_qa_rec-1; i++)
    {
      for ( j=0; j<4; j++)
	{
	  strcpy(Q[i][j], exo->qa_record[i][j]);
	}
    }

  /*
   * Preload new record with all terminating nulls.
   */

  for ( j=0; j<4; j++)
    {
      for ( k=0; k<MAX_STR_LENGTH; k++)
	{
	  Q[n][j][k] = '\0';
	}
    }

  /*
   * Concoct a new record.
   */

  strcpy(Q[n][0], "GOMA");
  strcpy(Q[n][1], VERSION); /* def'd in std.h for now */
  get_date(Q[n][2]);
  get_time(Q[n][3]);

  /*
   * Free the old record members.
   */

  for ( i=0; i<n; i++)
    {
      for ( j=0; j<4; j++)
	{
	  safer_free((void **)&(exo->qa_record[i][j]));
	}
    }

  if ( n > 0 )
    {
      safer_free((void **)&(exo->qa_record));
    }

  /*
   * Assign new pointer.
   */

  exo->qa_record = Q;

  return;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * add_info_stamp -- add an information record to EXODUS II db if space exists.
 *
 * Created: 1997/08/04 08:22 MDT pasacki@sandia.gov
 */

void
add_info_stamp(Exo_DB *exo)
{
  int i;
  int k;
  int n;
  char **a;
  char buf[MAX_LINE_LENGTH+1];
  time_t now, then;
#ifdef NO_LEAKY_GETPWUID
  struct passwd *pwe;
#endif
  struct utsname utsname;
  INFO_Record *I;

  n = exo->num_info;

  if ( n+10 > MAX_INFO )
    {
      return;
    }

  exo->num_info += 10;

  /*
  buf[0] = '\0';
  ni = (char **) calloc( exo->num_info, sizeof(char *));
  */

  /*
   * Allocate space for the new info record.
   */

  I = (INFO_Record *) smalloc(exo->num_info*sizeof(INFO_Record));

  for ( i=0; i<exo->num_info; i++)
    {
      I[i] = (char *) smalloc(MAX_LINE_LENGTH*sizeof(char));
    }

  /*
   * Transcribe any old records...
   */

  for ( i=0; i<n; i++)
    {
      strcpy(I[i], exo->info[i]);
    }

  /*
   * Initialize the new records to terminating nulls...
   */

  for ( i=n; i<exo->num_info; i++)
    {
      for ( k=0; k<MAX_LINE_LENGTH; k++)
	{
	  I[i][k] = '\0';
	}
    }

  /*
   * Fill in the new records with information about this run.
   */

  strcpy(I[n], "____");

  /*
   * -9 -- the command line issued for this simulation
   */

  /*
  for ( i=0; i<Argc; i++)
    {
      strcat(buf, Argv[i]);
      strcat(buf, " ");
    }
  */

  for ( i=0; i<MAX_LINE_LENGTH; i++)  
    {
      buf[i] = '\0';
    }

  a = Argv;

  if ( a != NULL )
    {
      i = 0;
      /* MMH: This needs a +2: +1 for the space, and +1 for the
       * null.
       */
      while (*a != NULL && 
	     ((strlen(buf) + strlen(*a) + 2) < MAX_LINE_LENGTH) )
	{
	  strcat(buf, *a);
	  strcat(buf, " ");
	  a++;
	}
    }

  strcpy(I[n+1], buf);

  /*
   * -8 -- the date and time of the simulation
   */

  now = time(&then);
  strftime(buf, MAX_LINE_LENGTH, "%C", localtime(&now));
  strcpy(I[n+2], buf);

  /*
   * -7 -- current working directory
   */

  if ( ProcID < 8 )		/* too much I/O overhead for many procs */
    {
      getcwd(buf, MAX_LINE_LENGTH+1);
    }
  else
    {
      strcpy(buf, ".");
    }

  strcpy(I[n+3], buf);

  /*
   * -6 -- the name of the user
   */

  sprintf(buf, "uid %d", (int)(getuid()));

#ifdef NO_LEAKY_GETPWUID
  pwe = getpwuid(getuid());
  strcpy(buf, pwe->pw_name);
#endif

  strcpy(I[n+4], buf);

  /*
   * -5 through -1 -- the POSIX system information
   */

  uname(&utsname);
  
  /* Lets be on the safe side put one char less than MAX_LINE_LENGTH into the info buffers */

  strncpy(I[n+5], utsname.sysname, MAX_LINE_LENGTH-1);
  strncpy(I[n+6], utsname.nodename, MAX_LINE_LENGTH-1);
  strncpy(I[n+7], utsname.release, MAX_LINE_LENGTH-1);
  strncpy(I[n+8], utsname.version, MAX_LINE_LENGTH-1);
  strncpy(I[n+9], utsname.machine, MAX_LINE_LENGTH-1);

  /*
   * Free the old beast and assign the new one.
   */
  for ( i=0; i<n; i++)
    {
      safer_free((void **)&(exo->info[i]));
    }
  if ( n > 0 )
    {
      safer_free((void **) &(exo->info));
    }
  exo->info = I;
  return;
}
/***********************************************************************/
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
  char err_msg[MAX_CHAR_IN_INPUT];
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
			  sprintf(err_msg, 
                                  "ex_put_elem_var() bad rtn: time %d, elemvar %d, EB ID %d",
                                  time_index, k+1, exo->eb_id[j]);
			  EH(-1, err_msg);
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
