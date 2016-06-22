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

/* fix -- GOMA monolith problem recomposition from polylithic EXODUS II data
 *
 * Copyright (c) 1999-2000 Sandia National Laboratories.  All rights reserved.
 *
 *	fix - using n augmented EXODUS II files with results, reconstruct a 
 *            monolithic EXODUS II file for the global problem that includes
 *            both mesh and results.
 *
 * Created:  1997/11/12 12:55 MST pasacki@sandia.gov
 *
 * Revised:  1998/08/29 07:22 MDT pasacki@sandia.gov
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "goma.h"

#include "brkfix/brkfix.h"		/* useful general stuff */
#include "mm_eh.h"			/* error handling */
#include "rf_allo.h"		/* multi-dim array allocation */
#include "exo_struct.h"		/* some definitions for EXODUS II */
#include "brkfix/bbb.h"
#include "dpi.h"		/* distributed processing information */
#include "brkfix/nodesc.h"		/* node descriptions */
#include "rd_dpi.h"

#include "brkfix/fix.h"

char *program_name;		/* name this program was run with */

static const char program_description[] = "GOMA distributed problem recomposition tool";

static const char copyright[] = 
"Copyright (c) 1999-2000 Sandia National Laboratories. All rights reserved.";

int *ep;			/* element pointers into node list */
int *np;			/* node pointers into element list */

int *nl;			/* node list */
int *el;			/* element list */

int *ebl;			/* element block list */

/*
 * Prototypes of functions defined in this file.
 */

static void setup_exo_res_desc	/* fix.c */
PROTO((Exo_DB *));		/* exo - ptr to database */


int
fix_exo_file(int num_procs, char* exo_mono_name)
{
  int i;
  int p, pmax = 0, num_node_var_max=0;
  int t;

  Exo_DB *mono;			/* monolith mesh */
  Exo_DB *poly;			/* polylith mesh+dpi+results in */

  Dpi *dpin;			/* polylith dpi in */

  char  monolith_file_name  [FILENAME_MAX_ACK]; /* original mesh */

  char  polylith_name       [FILENAME_MAX_ACK]; /* "basename_1of2.exoII" */

  char  err_msg[MAX_CHAR_ERR_MSG];
  char *tmp;			/* char pointer junkyard of no interest */

  Spfrtn sr=0;

  ELEM_BLK_STRUCT *element_blocks_save = Element_Blocks;
  /*
   * Defaults
   */

  for ( i=0; i<FILENAME_MAX_ACK; i++)
    {
      monolith_file_name[i] = '\0';
      polylith_name[i] = '\0';
    }

  tmp = strcpy(monolith_file_name, exo_mono_name);

  if ( num_procs < 1 ) {
    sr = sprintf(err_msg, "Bad number of processors specified: %d.",
                 num_procs);
    EH(-1, err_msg);
  } else if (num_procs == 1) {
    sprintf(err_msg, "fix_exo_file(), %s, called with %d processors", exo_mono_name, num_procs);
    WH(-1, err_msg);
    return -1;
  }


  /*
   * Turn off annoying error reporting from within the EXODUS II API...
   */

  ex_opts(EX_VERBOSE); 

  /* Here we will loop through the pieces and find that which represents the most
   * nodal variables to help with sizing of the monolith.  PRS-6/1/2010 
   */

  for ( p=0; p<num_procs; p++)
    { 
      strcpy(polylith_name, exo_mono_name);

      strcpy(monolith_file_name, polylith_name);

      multiname(polylith_name, p, num_procs);

#ifdef DEBUG
      /*  err = sscanf(line, "%s", in_exodus_file_name); */
      fprintf(stderr, "monolith EXODUSII file = \"%s\"\n", monolith_file_name);
#endif

      /*
       * Make preliminary space for the monolithic database and the
       * polylithic pieces.
       *
       * The strategy will be to look closely at the first polylithic piece
       * in order to ascribe basic sizing information to the monolith.
       *
       * Then, the polyliths are traversed to slowly build the mesh and the
       * results of the monolith. Several passes may be necessary...
       */

      poly = (Exo_DB *) smalloc(sizeof(Exo_DB));
      dpin = (Dpi *) smalloc(sizeof(Dpi));

      /*
       * Open the first polylith to determine if there are many timeplanes
       * of data.
       *
       * Assume that the first polylith is representative, that the number
       * of nodal variables, the number of timeplanes, etc will be the same
       * for each and every polylith encountered hereafter.
       *
       * The first polylith is important, too, in that much of the sizing
       * information for the global problem will be derived from it.
       */

      init_exo_struct(poly);

      init_dpi_struct(dpin);

      /*
       * Some defaults are different for polylithic pieces - they have use these
       * maps to relate themselves to the monolith...
       *
       * Now, these maps are part of the Dpi information and no longer piggybacked
       * inside EXODUS II. Thus, don't attempt to read them if they aren't there.
       */

#ifdef DEBUG  
      fprintf(stderr, "Fix: attempting to build a %d piece %s\n",
	      num_procs, monolith_file_name);
#endif

      /*
       * Read everything in the 1st polylith's EXODUS information except for
       * the results data per se...
       *
       * Set to zero the variables that get used for allocating and writing...
       */

      rd_exo(poly, polylith_name, 0, ( EXODB_ACTION_RD_INIT + 
				       EXODB_ACTION_RD_MESH + 
				       EXODB_ACTION_RD_RES0 ));

      if (poly->num_node_vars > num_node_var_max) 
	{
	  num_node_var_max = poly->num_node_vars;
	  pmax = p; 
	}

      zero_base(poly);

      rd_dpi(dpin, polylith_name, 0);
      free_dpi(dpin);
      free(dpin);

      free_element_blocks(poly);

      free_exo(poly);
      free(poly);
   
    }

  /*
   * Now that we know which piece to use, Build the monolithic skeleton...
   */
  strcpy(polylith_name, exo_mono_name);

  strcpy(monolith_file_name, polylith_name);

  multiname(polylith_name, pmax, num_procs);

  poly = (Exo_DB *) smalloc(sizeof(Exo_DB));
  dpin = (Dpi *) smalloc(sizeof(Dpi));

  init_exo_struct(poly);

  init_dpi_struct(dpin);

  rd_exo(poly, polylith_name, 0, ( EXODB_ACTION_RD_INIT + 
				       EXODB_ACTION_RD_MESH + 
				       EXODB_ACTION_RD_RES0 ));
  zero_base(poly);
  rd_dpi(dpin, polylith_name, 0);


  mono = (Exo_DB *) smalloc(sizeof(Exo_DB));
  memset(mono, 0, sizeof(Exo_DB));

  init_exo_struct(mono);

  /*
   * This fills in sketchy material like ex_get_init(), as well as
   * various array allocations and initializations in preparation for
   * the polylith sweep...
   */

  build_big_bones(poly, dpin, mono);

  free_dpi(dpin);
  free(dpin);

  free_element_blocks(poly);

  free_exo(poly);
  free(poly);

  for ( p=0; p<num_procs; p++)
    {
      poly = (Exo_DB *) smalloc(sizeof(Exo_DB));
      dpin = (Dpi *) smalloc(sizeof(Dpi));

      init_dpi_struct(dpin);
      init_exo_struct(poly);

      for ( i=0; i<FILENAME_MAX_ACK; i++)
	{
	  polylith_name[i] = '\0';
	}

      strcpy(polylith_name, exo_mono_name);
      multiname(polylith_name, p, num_procs);

      /*
       * Set actions...
       */

      rd_exo(poly, polylith_name, 0, ( EXODB_ACTION_RD_INIT + 
				       EXODB_ACTION_RD_MESH ) );
      zero_base(poly);

      rd_dpi(dpin, polylith_name, 0);

      build_global_coords(poly, dpin, mono);

      /*
       * Contribute to the element block data...
       */

      build_global_conn(poly, dpin, mono);
      build_global_attr(poly, dpin, mono);

      /*
       * Contribute to the node set node list and distribution factor list...
       */

      build_global_ns(poly, dpin, mono);

      /*
       * Contribute to the side set side, elem, and dist fact lists...
       */

      build_global_ss(poly, dpin, mono);

      free_element_blocks(poly);

      free_exo(poly);
      free(poly);

      free_dpi(dpin);
      free(dpin);
    }

  one_base(mono);
  wr_mesh_exo(mono, monolith_file_name, 0);
  wr_resetup_exo(mono, monolith_file_name, 0);
  zero_base(mono);

  /*
   * Now sweep through polyliths while there are timeplanes of results
   * and map those results into the monolith, then write them out one
   * at a time.
   */

  /* PRS Note (5/31/2010): this is where the memory for p->nv gets allocated
   * and it is based on the num_nod_vars set above */

  setup_exo_res_desc(mono);

#ifdef DEBUG
  fprintf(stderr, "mono->num_times = %d\n", mono->num_times);
#endif

  for ( t=0; t<mono->num_times; t++)
    {
      for ( p=0; p<num_procs; p++)
	{
	  poly = (Exo_DB *) smalloc(sizeof(Exo_DB));
	  dpin = (Dpi *) smalloc(sizeof(Dpi));

	  /*
	   * Initialize to semisanity...
	   */

	  init_dpi_struct(dpin);
	  init_exo_struct(poly);

	  for ( i=0; i<FILENAME_MAX_ACK; i++)
	    {
	      polylith_name[i] = '\0';
	    }
	  
	  strcpy(polylith_name, exo_mono_name);
	  multiname(polylith_name, p, num_procs);

	  /*
	   * Set actions - what to fetch from the polylith databases.
	   */
	  
	  rd_exo(poly, polylith_name, 0, ( EXODB_ACTION_RD_INIT +
					   EXODB_ACTION_RD_MESH +
					   EXODB_ACTION_RD_RES0 ));
	  zero_base(poly);
	  rd_dpi(dpin, polylith_name, 0);
	  
	  /*
	   * Now indicate what variables and time planes to read from the
	   * individual polyliths...
	   */

#ifdef DEBUG
	  fprintf(stderr, "\nBuilding results for proc=%d, time=%d\n", p, t);
#endif
	  setup_exo_res_desc(poly);

#ifdef DEBUG
	  fprintf(stderr, "A mono->nv_time_indeces[0] = %d\n", 
		  mono->nv_time_indeces[0]);
	  fprintf(stderr, "A poly->nv_time_indeces[0] = %d\n", 
		  poly->nv_time_indeces[0]);
#endif	  

	  /*
	   * Pick one timeplane to pick - this one!
	   */

	  if ( mono->num_glob_vars > 0 )
	    {
	      mono->gv_time_indeces[0] = t+1;
	    }

	  if ( mono->num_elem_vars > 0 )
	    {
	      mono->ev_time_indeces[0] = t+1;
	    }

	  if ( mono->num_node_vars > 0 )
	    {
	      mono->nv_time_indeces[0] = t+1;
	    }

	  if ( poly->num_glob_vars > 0 )
	    {
	      poly->gv_time_indeces[0] = t+1;
	    }

	  if ( poly->num_node_vars > 0 )
	    {
	      poly->nv_time_indeces[0] = t+1;
	    }

	  if ( poly->num_elem_vars > 0 )
	    {
	      poly->ev_time_indeces[0] = t+1;
	    }

#ifdef DEBUG
	  fprintf(stderr, "B mono->nv_time_indeces[0] = %d\n", 
		  mono->nv_time_indeces[0]);
	  fprintf(stderr, "B poly->nv_time_indeces[0] = %d\n", 
		  poly->nv_time_indeces[0]);
#endif
	  
	  /*
	   * This assignment will help rd_exo() figure out to allocate
	   * enough space to read in one timeplane with ALL the nodal
	   * variables that are in the database.
	   */

	  poly->num_nv_indeces = poly->num_node_vars;

	  rd_exo(poly, polylith_name, 0, ( EXODB_ACTION_RD_RESN +
					   EXODB_ACTION_RD_RESE +
					   EXODB_ACTION_RD_RESG ) );

#ifdef DEBUG
	  fprintf(stderr, "C mono->nv_time_indeces[0] = %d\n", 
		  mono->nv_time_indeces[0]);
	  fprintf(stderr, "C poly->nv_time_indeces[0] = %d\n", 
		  poly->nv_time_indeces[0]);
#endif	  

	  /*
	   * Map the polylith's results into the monolith...
	   */

	  build_global_res(poly, dpin, mono);

#ifdef DEBUG
	  fprintf(stderr, "D mono->nv_time_indeces[0] = %d\n", 
		  mono->nv_time_indeces[0]);
	  fprintf(stderr, "D poly->nv_time_indeces[0] = %d\n", 
		  poly->nv_time_indeces[0]);
#endif

	  free_dpi(dpin);
	  free(dpin);

          free_element_blocks(poly);

          free_exo_gv(poly);
          free_exo_nv(poly);
          free_exo_ev(poly);
	  free_exo(poly);
	  free(poly);
	}

      /*
       * Now, write out the global results at this particular timeplane.
       */

      one_base(mono);
#ifdef DEBUG
      fprintf(stderr, 
	      "About to wr_results at t=%d, mono->nv_time_indeces[0]=%d\n", 
	      t, mono->nv_time_indeces[0]);
#endif
      wr_result_exo(mono, monolith_file_name, 0);
      
      if ( mono->num_glob_vars > 0 ) {
	int status;
	
	mono->cmode = EX_WRITE;


	mono->io_wordsize   = 0;	/* i.e., query */
	mono->comp_wordsize = sizeof(dbl);
	mono->exoid         = ex_open(mono->path, 
				   mono->cmode, 
				   &mono->comp_wordsize, 
				   &mono->io_wordsize, 
				   &mono->version);

	/*status = ex_put_var_param(mono->exoid, "g", mono->num_glob_vars);
	EH(status, "ex_put_var_param(g)");
	*/
	status = ex_put_var_names(mono->exoid, "g", mono->num_glob_vars,
				  mono->glob_var_names);
	EH(status, "ex_put_var_names(g)");



	for (i = 0; i < mono->num_gv_time_indeces; i++) {
	  status = ex_put_glob_vars(mono->exoid, mono->gv_time_indeces[i],
			   mono->num_glob_vars,
			   mono->gv[i]);
	  EH(status, "ex_put_glob_vars");
	}

	status = ex_close(mono->exoid);
	EH(status, "ex_close()");

      }
      
      zero_base(mono);
    }

  free_exo_gv(mono);
  free_exo_nv(mono);
  free_exo_ev(mono);
  free_exo(mono);
  free(mono);

  if ( tmp == NULL ) exit(2);
  if ( sr < 0 ) exit(2);

  /* Restore Element_Blocks */

  Element_Blocks = element_blocks_save;

  return(0);
}

/* setup_exo_res_desc() -- allocate, set arrays to rd/wr all results, 1 time
 *
 * Allocate and setup arrays for reading and writing EXODUS II results
 * variables (nodal, global, element).
 *
 * Default values for some quantities are inserted after allocation. It is
 * the user's responsibility to fill them in later. I.e, things like
 * ev_time_indeces[], etc.
 *
 * These auxiliary arrays in the Exo_DB structure give instructions to
 * rd_exo and wr_exo as to what to read and write, without the need to
 * explicitly interact with the EXODUS II API.
 *
 *
 * Created: 1998/08/14 06:04 MDT pasacki@sandia.gov
 *
 * Revised: 1998/08/15 12:43 MDT pasacki@sandia.gov
 */

static void 
setup_exo_res_desc(Exo_DB *exo)
{
  if ( ! ( exo->state & EXODB_STATE_RES0 ) )
    {
      EH(-1, "Setup to read bulk results requires preliminary info.");
    }

  /*
   * Global variables...
   */

  if ( ! ( exo->state & EXODB_STATE_GBIA ) )
    {
      alloc_exo_gv(exo, 1);
    }

  /*
   * Nodal variables...
   */

  if ( ! ( exo->state & EXODB_STATE_NDVA ) )
    {
      /*
       * Read only one timeplane's worth of nodal variables at 
       * a time, but read all of the nodal variables that are 
       * available. To do so, allocate enough space. Note that
       * allocation of x->nv[][][] will use info in
       * x->num_nv_time_indeces, etc.
       */

      if ( exo->num_node_vars > 0 )
	{
	  if ( exo->state & EXODB_STATE_NDIA )
	    {
	      free(exo->nv_time_indeces);
	      free(exo->nv_indeces);
	      exo->state &= ~EXODB_STATE_NDIA;
	    }

	  exo->num_nv_time_indeces = 1;
	  exo->num_nv_indeces      = exo->num_node_vars;
	  /*
	   * Flag indicates we've allocated this memory...
	   */
	  
	  exo->state |= EXODB_STATE_NDIA;
      
	  if ( exo->num_nv_time_indeces > 0 &&
	       exo->num_nv_indeces      > 0 &&
	       exo->num_nodes           > 0 )
	    {
	      alloc_exo_nv(exo, exo->num_nv_time_indeces, exo->num_nv_indeces);
	    }
	}
    }

  /*
   * Element variables...
   */

  if ( ! ( exo->state & EXODB_STATE_ELVA ) )
    {
      alloc_exo_ev(exo, 1);
    }

  return;
}
