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

#define _FIX_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <unistd.h>

#include <string.h>

#include "std.h"		/* useful general stuff */
#include "eh.h"			/* error handling */
#include "aalloc.h"		/* multi-dim array allocation */
#include "exo_struct.h"		/* some definitions for EXODUS II */
#include "dpi.h"		/* distributed processing information */
#include "nodesc.h"		/* node descriptions */
#include "rd_dpi.h"

static int show_help			= FALSE;
static int be_quiet			= FALSE;
static int be_verbose			= FALSE;
static int num_procs_specified          = FALSE;

char *program_name;		/* name this program was run with */

static char program_version[] = BRK_VERSION;

const char program_description[] = "GOMA distributed problem recomposition tool";

const char copyright[] = 
"Copyright (c) 1999-2000 Sandia National Laboratories. All rights reserved.";

int *ep;			/* element pointers into node list */
int *np;			/* node pointers into element list */

int *nl;			/* node list */
int *el;			/* element list */

int *ebl;			/* element block list */

extern int rd_exo		/* def in "rd_exo.c" */
PROTO((Exo_DB *,		/* structure defined in exo_struct.h */
       char *,			/* filename of exo2 file */
       int,			/* verbosity flag */
       int));			/* task */

extern void init_exo_struct	/* "rd_exo.c" */
PROTO((Exo_DB *));		/* ptr to whole database in memory */

extern int free_exo		/* "rd_exo.c" */
PROTO((Exo_DB *));		/* pointer to EXODUS II FE db structure */

extern void alloc_exo_nv
PROTO((Exo_DB *));

extern void alloc_exo_gv
PROTO((Exo_DB *, int ));

extern void alloc_exo_ev
PROTO((Exo_DB *, int ));

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

extern void zero_base		/* utils.c */
PROTO((Exo_DB *));		/* E - pointer to an EXODUS II database */

extern void one_base		/* utils.c */
PROTO((Exo_DB *));		/* x = pointer to an EXODUS II FE database */

extern int get_filename_num_procs /* utils.c */
PROTO((const char *));		/* basename - of polylithic files */

extern int rd_dpi		/* rd_dpi.c */
PROTO((Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* fn - filename */
       int ));			/* verbosity - how much to talk */

extern int wr_dpi		/* wr_dpi.c */
PROTO((Dpi *,			/* fantastic structure defd in "dpi.h" */
       char *,			/* filename */
       int ));			/* verbosity - how much to talk */

extern void free_dpi		/* rd_dpi.c */
PROTO((Dpi *));			/* fantastic structure defd in "dpi.h" */

extern void free_dpi_uni	/* rd_dpi.c */
PROTO((Dpi *));			/* fantastic structure defd in "dpi.h" */


extern void multiname		/* rd_mesh.c */
PROTO((char *,			/* in_name - generic global name "pref.suf" */
       int,			/* integer processor_name */
       int));			/* number_processors - total */

extern void strip_suffix	/* rd_mesh.c */
PROTO((char *,			/* result - "a" */
       char *));		/* in -- input string "a.b" */

extern void get_suffix		/* rd_mesh.c */
PROTO((char *,			/* result -- "b" */
       char *));		/* in -- extract the tail of "a.b" -> "b" */

extern void build_big_bones
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_conn	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_attr	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_coords	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_ns	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_ss	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

extern void build_global_res	/* bbb.c */
PROTO((Exo_DB *,		/* EXODUS info from representative polylith */
       Dpi    *,		/* distributed processing info from polylith */
       Exo_DB *));		/* ptr to monolithic EXODUS II database */

/*
 * Prototypes of functions defined in this file.
 */

static void usage
PROTO((const int));		/* status */

static void setup_exo_res_desc	/* fix.c */
PROTO((Exo_DB *));		/* exo - ptr to database */

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
      fprintf(stdout, 
	      "Usage: %s [OPTIONS] basename\n", 
	      program_name);
      fprintf(stdout, "\n");
      fprintf(stdout, "Reconstruct a monolithic EXODUS II file from\n");
      fprintf(stdout, "multiple polylithic EXODUS II files of the form\n");
      fprintf(stdout, "basename_1ofn.exoII, ..., basename_nofn.exoII.\n");
      fprintf(stdout, "\n");
      fprintf(stdout, "OPTIONS:\n");
      fprintf(stdout, "\t-h        show options summary (this message)\n");
      fprintf(stdout, "\t-n P      explicitly assemble from P pieces\n");
      fprintf(stdout, "\t          (default deduces P from basename and files in cwd)\n");
      fprintf(stdout, "\n");
      fprintf(stdout, "\t-o file   name for reconstructed monolith\n");
      fprintf(stdout, "\t          (default is basename.exoII)\n");
      fprintf(stdout, "\t-v        print code version\n");
      fprintf(stdout, "\n");
    }
  exit(status);
}

int
main (int argc, char *argv[], char *envp[])
{
  int c;			/* hold each option flag */
  int err;
  int i;
  int num_procs=1;		/* from which to build the monolith */
  int p, pmax, num_node_var_max=0;
  int status  = 0;		/* in case anything goes wrong */
  int t;

  Exo_DB *mono;			/* monolith mesh */
  Exo_DB *poly;			/* polylith mesh+dpi+results in */

  Dpi *dpin;			/* polylith dpi in */

  char  monolith_file_name  [FILENAME_MAX_ACK]; /* original mesh */
  char  polylith_basename   [FILENAME_MAX_ACK]; /* broken results */

  char  polylith_name       [FILENAME_MAX_ACK]; /* "basename_1of2.exoII" */

  char  err_msg[MAX_CHAR_ERR_MSG];
  char *tmp;			/* char pointer junkyard of no interest */

  extern char *optarg;
  extern int   optind;

  Spfrtn sr=0;

  program_name = argv[0];

  /*
   * Defaults
   */

  for ( i=0; i<FILENAME_MAX_ACK; i++)
    {
      monolith_file_name[i] = '\0';
      polylith_basename[i] = '\0';
      polylith_name[i] = '\0';
    }

  tmp = strcpy(monolith_file_name, "in.exoII");

  tmp = strcpy(polylith_basename, "out");

  if ( argc < 2 )
    {
      usage(0);
    }

  /*
   * Process any command line options.
   */

  while ( ( c = getopt(argc, argv, "b:hn:o:v")) != EOF )
    {
#ifdef DEBUG
	  fprintf(stderr, "Processing with argc=%d\n", argc);
#endif
      switch(c)
	{

	case 'h':
	  /*
	   * Show synopsis of help.
	   */
	  show_help |= TRUE;
	  usage(status);
	  break;

	case 'n':
	  err = sscanf(optarg, "%d", &num_procs);
#ifdef DEBUG
	  fprintf(stderr, "Processing with optind=%d, optarg=%s, argc=%d\n",
		  optind, optarg, argc);
#endif
	  if ( err != 1 || num_procs < 1 )
	    {
	      EH(-1, 
		 "For number of pieces, specify one positive nonzero integer.");
	    }
	  num_procs_specified = TRUE;
	  break;

	case 'o':

	  tmp = strcpy(monolith_file_name, optarg);
#ifdef DEBUG
	  fprintf(stderr, "Processing with optind=%d, optarg=%s, argc=%d\n",
		  optind, optarg, argc);
#endif
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

  if ( optind > argc )
    {
      fprintf(stderr, "optind = %d, argc = %d\n", optind, argc);
      usage(-1);
    }

#ifdef DEBUG
  fprintf(stderr, "Opts done: optind=%d, argc=%d\n",
	  optind, argc);      /* , optind, argv[optind]);*/
#endif

  tmp = strcpy(polylith_basename, argv[optind]);
  optind++;

  /*
   * If the number of processors was not specified, then try to figure
   * it out...
   */

  if ( ! num_procs_specified )
    {
      num_procs = get_filename_num_procs(polylith_basename);
      fprintf(stderr, "Found %d pieces\n", num_procs);
    }

  if ( num_procs < 1 )
    {
      sr = sprintf(err_msg, "Bad number of processors specified: %d.",
		   num_procs);
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
   * Turn off annoying error reporting from within the EXODUS II API...
   */

  ex_opts(EX_VERBOSE); 

  /* Here we will loop through the pieces and find that which represents the most
   * nodal variables to help with sizing of the monolith.  PRS-6/1/2010 
   */

  for ( p=0; p<num_procs; p++)
    { 
      strcpy(polylith_name, polylith_basename);
      strcat(polylith_name, ".exoII");

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
      fprintf(stderr, "Fix: attempting to build a %d piece %s from %s pieces\n",
	      num_procs, monolith_file_name, polylith_basename);
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

      free_exo(poly);
      free(poly);
   
    }

  /*
   * Now that we know which piece to use, Build the monolithic skeleton...
   */
  strcpy(polylith_name, polylith_basename);
  strcat(polylith_name, ".exoII");

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

  init_exo_struct(mono);

  /*
   * This fills in sketchy material like ex_get_init(), as well as
   * various array allocations and initializations in preparation for
   * the polylith sweep...
   */

  build_big_bones(poly, dpin, mono);

  free_dpi(dpin);
  free(dpin);

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

      strcpy(polylith_name, polylith_basename);
      strcat(polylith_name, ".exoII");
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
	  
	  strcpy(polylith_name, polylith_basename);
	  strcat(polylith_name, ".exoII");
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
      zero_base(mono);
    }

  free_exo(mono);
  free(mono);

  if ( tmp == NULL ) exit(2);
  if ( sr < 0 ) exit(2);

  return(0);
}

/* multiname() -- translate filename string to distributed processing version
 *
 * 
 * Description:
 *
 * Many data file names will be unique to a given processor. Construct that
 * name for this processor. The names will be translate like this
 *
 * Old name: "basename.suffix"
 *
 * New name: "basename_13of437.suffix"
 *
 * Note: The input string is assumed to have sufficient space allocated to
 *       contain the revised name.
 *
 *	 Processors, while named beginning at zero, are incremented so that
 *       the string reads 1ofn, 2ofn, ...., nofn and NOT
 *	 0ofn, 1ofn, ..., n-1ofn.
 *
 *
 * Created: 1997/07/09 13:18 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
multiname(char *in_name, 
	  int processor_name, 
	  int number_processors)
{
  char err_msg[1024];

  Spfrtn sr=0;

  if ( number_processors == 1 ) return;

  if ( processor_name < 0 ) 
    {
      sr = sprintf(err_msg, "processor_name = %d ( less than zero).", 
		   processor_name);
      EH(-1, err_msg);
      EH(sr, err_msg);
    }
  else if ( processor_name > number_processors - 1 )
    {
      sr = sprintf(err_msg, "processor_name = %d ( too high ).", 
		   processor_name);
      EH(-1, err_msg);
    }

  if ( number_processors < 1 )
    {
      sr = sprintf(err_msg, "number_processors = %d ( less than one ).", 
		   number_processors);
      EH(-1, err_msg);
    }

  sprintf(in_name, "%s.%d.%d", in_name, number_processors, processor_name);
  return;
}

/*
 * strip_suffix() -- chop off trailing characters and the period in a string
 *
 * Description:
 *	The input string is examined from back to front until a period is
 * found or the beginning of the string is reached. A new string is allocated
 * and filled with everything but the suffix. Thus "a.b" -> "a"
 *
 * Created: 1997/07/09 16:10 MDT pasacki@sandia.gov
 */

void
strip_suffix(char *result, 
	     char *in)
{
  int i,j;
  int e;

  e = strlen(in);

  if ( e < 1 ) exit(-1);

  /*
   * Starting from the end of the string, look back until we find a
   * period "." character.
   */

  i = e;

  while ( i>0 && *(in+i) != '.' ) i--;

  for ( j=0; j<i; j++)
    {
      result[j] = *(in+j);
    }
  result[i] = '\0';

  /*
   * No suffix found? Then the whole string is the basename.
   */

  if ( i == 0 )
    {
      strcpy(result, in);
    }

  return;
}

/*
 * get_suffix() -- extract the tail substring of a string past the last period
 *
 * Created: 1997/07/09 16:14 MDT pasacki@sandia.gov
 */

void
get_suffix(char *result, 
	   char *in)
{
  int i;
  int b;
  int e;

  e = strlen(in);

  b = e;

  while ( b>0 && *(in+b) != '.' ) b--;

  if ( b == 0 )
    {
      *result = '\0';
    }
  else
    {
      for ( i=b; i<e; i++)
	{
	  result[i-b] = in[i+1];
	}
    }

  return;
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
  int i;

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
	  exo->nv_time_indeces     = (int *) 
	    smalloc(exo->num_nv_time_indeces*sizeof(int));

	  for ( i=0; i<exo->num_nv_time_indeces; i++)
	    {
	      exo->nv_time_indeces[i] = -1; /* Initialize to undefined. */
	    }
	  
	  exo->num_nv_indeces      = exo->num_node_vars;
	  exo->nv_indeces          = (int *) smalloc(exo->num_nv_indeces*
						     sizeof(int));

	  /*
	   * Flag indicates we've allocated this memory...
	   */
	  
	  exo->state |= EXODB_STATE_NDIA;
      
	  for ( i=0; i<exo->num_nv_indeces; i++)
	    {
	      exo->nv_indeces[i]   = i+1; /* 1-based naming sop to FORTRAN */
	    }

	  if ( exo->num_nv_time_indeces > 0 &&
	       exo->num_nv_indeces      > 0 &&
	       exo->num_nodes           > 0 )
	    {
	      alloc_exo_nv(exo);
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
