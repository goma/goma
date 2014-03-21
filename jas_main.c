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
 
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "std.h"


#ifdef HAVE_SUNMATH
#include <sunmath.h>
#endif

#ifdef PARALLEL
#include "az_aztec.h"
#endif
extern void handle_ieee(void );

#include "el_geom.h"
#include "rf_allo.h"
#include "rf_solver.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"

#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_io_defn.h"

#include "rf_bc_const.h"
#include "rf_element_storage_struct.h"
#include "mm_elem_block.h"
#include "rf_vars_const.h"

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h" 
#endif
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "rf_node.h"
#include "mm_as.h"

#include "mm_eh.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_chemkin.h"

/* #include "mm_names.h" -- use extern defs from rf_bc_const.h for these vars */

#include "exo_struct.h"
#include "dpi.h"

#include "dp_utils.h"
#include "dp_types.h"


#define _MAIN_C
#include "goma.h"

/*
 * Global variables defined here.
 */

Dpi *DPI_ptr = NULL;

Exo_DB *EXO_ptr = NULL;

#define LAST_LEGAL_STRING "Last Legal String"

/* 
 * Do not change the line below. It is automatically transformed
 * on its way in and out of the repository so that we have a variable
 * available that indicates when this version of goma came out of the
 * repository.
 */
/*
static int CVS_Extraction_Date_Stamp[] = { 2002, 9, 20, 21, 17, 15};
*/

static char *legal_notice[] = {
  "Copyright (c) 2013  Sandia Corporation.\n",
  "\n", 
  "Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive\n",
  "license for use of this work by or on behalf of the U.S. Government. Export\n",
  "of this program may require a license from the United States Government.\n",
  "\n", 
  "PRS, PAS, RRR, KSC, RAC, TAB, DRN, PLH, DAL, IDG, ACK, ACS, RBS, RRL, MMH, SRS\n",
  "HKM, RAR, EDW, PKN, SAR, EMB ...\n"
  "\n", 
  LAST_LEGAL_STRING
};

/*
 * Variables defined here but everywhere else are extern via "rf_mp.h".
 */

int parallel_err        = FALSE;
int parallel_err_global = FALSE;

int Unlimited_Output    = TRUE;          /* print limit flag */

int Num_Proc = 1;
int ProcID = 0;
int Dim = -1;

int unlink(const char *);
int unlerr ;

MPI_Request *Request = NULL;
MPI_Status *Status = NULL;
MPI_Aint type_extent;
int Num_Requests = 6;

/*
 * Data structures for transporting input data to other processors...
 */

DDD Noahs_Raven;		/* Stage 1 - preliminary sizing information */
DDD Noahs_Ark;			/* Stage 2 - big boatload - main body*/
DDD Noahs_Dove;			/* Stage 3 - doubly dynamic sized stuff */

Comm_Ex *cx = NULL;		/* communications info for ea neighbor proc */

double time_goma_started;	/* Save it here... */

/* Lets not be sloppy people.  Put in your prototypes */

void initialize_goma_export_vars__ 
PROTO(( int *,
	int *,
	int *,
	int *,
	int *,
	int *,
	double *,
	double *,
	double * ));


/*
 * Declare these copies so that other routines may access them as global
 * external variables instead of passing them through the argument lists.
 */

char **Argv;
int Argc;

int
goma_init_(dbl *time1, int *nnodes, int *nelems,
           int *nnv_in, int *nev_in, int *i_soln, int *i_post)
     
     /*
      * Initial main driver for GOMA. Derived from a (1/93) release of
      * the rf_salsa program by
      *        
      *        Original Authors: John  Shadid (1421)
      *		                 Scott Hutchinson (1421)
      *        		         Harry Moffat (1421)
      *       
      *        Date:		12/3/92
      * 
      *
      *        Updates and Changes by:
      *                           Randy Schunk (9111)
      *                           P. A. Sackinger (9111)
      *                           R. R. Rao       (9111)
      *                           R. A. Cairncross (Univ. of Delaware)
      *        Dates:           2/93 - 6/96
      *
      *       Modified for continuation
      *                           Ian Gates
      *       Dates:            2/98 - 10/98
      *       Dates:            7/99 - 8/99
      * 
      * Last modified: Wed  June 26 14:21:35 MST 1994 prschun@sandia.gov
      * Hello.
      * 
      * Note: Many modifications from an early 2/93 pre-release
      *	      version of rf_salsa were made by various persons 
      *       in order to test ideas about moving/deforming meshes...
      */ 
{
  /* Local Declarations */

  double time_start, total_time;   /* timing variables */
#ifndef PARALLEL
  struct tm *tm_ptr;               /* additional serial timing variables */
  time_t the_time;
#endif

  int error;
  int i;
  int j;
  static int first_goma_call=TRUE;

  char	**ptmp;
  static const char *yo="goma_init";

  struct Command_line_command **clc=NULL; /* point to command line structure */
  int           nclc = 0;		/* number of command line commands */

/********************** BEGIN EXECUTION ***************************************/
  
/* assume number of commands is less than or equal to the number of 
 * arguments in the command line minus 1 (1st is program name) */

  /*
  *  Get the name of the executable, yo
  */

#ifdef PARALLEL
if( first_goma_call ) {
	Argc = 1;
	Argv = (char **) smalloc( Argc*sizeof(char *) );
	Argv[0] = (char *) yo;
	MPI_Init(&Argc, &Argv);  /*PRS will have to fix this.  Too late TAB already did. */
  }
  time_start = MPI_Wtime();
#else /* PARALLEL */
  (void) time(&the_time);
  tm_ptr = gmtime(&the_time);
  time_start = (double)  ( tm_ptr->tm_sec
               + 60. * (   60. * ( tm_ptr->tm_yday * 24. + tm_ptr->tm_hour )
                                                         + tm_ptr->tm_min  )
                         );
#endif /* PARALLEL */
  *time1 = time_start;

/*   Argv = argv; */

/*   Argc = argc; */

  time_goma_started = time_start;

#ifdef PARALLEL
  /*
   * Determine the parallel processing status, if any. We need to know
   * pretty early if we're "one of many" or the only process.
   */

  error = MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);
  error = MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);

  /*
   * Setup a default Proc_config so we can use utility routines 
   * from Aztec
   */

  AZ_set_proc_config(Proc_Config, MPI_COMM_WORLD);

  /* set the output limit flag if need be */

  if( Num_Proc > DP_PROC_PRINT_LIMIT ) Unlimited_Output = FALSE;

#ifdef HAVE_MPE_H
  error = MPE_Init_log();
#endif /* HAVE_MPE_H */

  Dim = 0;			/* for any hypercube legacy code...  */

#endif /* PARALLEL */
  
#ifndef PARALLEL
  Dim        = 0;
  ProcID     = 0;
  Num_Proc   = 1;
#endif /* PARALLEL */


  /*
  *   HKM - Change the ieee exception handling based on the machine and
  *         the level of debugging/speed desired. This call currently causes
  *         core dumps for floating point exceptions.
  */

  handle_ieee();
  
  log_msg("--------------");
  log_msg("GOMA begins...");

#ifdef USE_CGM
  cgm_initialize();
#endif
  /*
   * Some initial stuff that only the master process does.
   */

/*PRS: Disable this command line stuff for the jas coupled version */
/*-----------------------------------------------------------------*/
/*   if ( ProcID == 0 ) */
/*     { */
/*       if (argc > 1) */
/* 	{ */
/* 	  log_msg("Preprocessing command line options."); */
/* 	  clc = (struct Command_line_command **)  */
/* 	    smalloc( argc * sizeof(struct Command_line_command *)); */
/* 	  for (i=0; i<argc; i++) */
/* 	    { */
/* 	      clc[i] = (struct Command_line_command *)  */
/* 		smalloc(sizeof(struct Command_line_command)); */
/* 	      clc[i]->type   = 0; /\* initialize command line structure *\/ */
/* 	      clc[i]->i_val  = 0; */
/* 	      clc[i]->r_val  = 0.; */
/* 	      clc[i]->string = (char *)  */
/* 		smalloc(MAX_COMMAND_LINE_LENGTH*sizeof(char)); */
/* 	      for ( j=0; j<MAX_COMMAND_LINE_LENGTH; j++) */
/* 		{ */
/* 		  clc[i]->string[j] = '\0'; */
/* 		} */
/* #ifdef DEBUG */
/* 	      fprintf(stderr, "clc[%d]->string is at 0x%x\n", i, clc[i]->string); */
/* 	      fprintf(stderr, "clc[%d]         is at 0x%x\n", i, clc[i]); */
/* #endif */
/* 	    } */
/* 	} */

/* PRS For the JAS version we will use the default input file name "input" */
      strcpy(Input_File, "input");

/* if (argc > 1) translate_command_line(argc, argv, clc, &nclc); */
      
/*       print_code_version(); */
/*       ptmp = legal_notice; */
/*       while ( strcmp(*ptmp, LAST_LEGAL_STRING) != 0 ) */
/* 	{ */
/* 	  fprintf(stderr, "%s", *ptmp++); */
/* 	} */
/* } */

  /*
   *  Allocate the uniform problem description structure and
   *  the problem description structures on all processors
   */
  error = pd_alloc();
  EH(error, "pd_alloc problem");

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier after pd_alloc\n", ProcID);
#ifdef PARALLEL
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  log_msg("Allocating mp, gn, ...");

  error = mp_alloc();
  EH(error, "mp_alloc problem");

  error = gn_alloc();
  EH(error, "gn_alloc problem");

  error = ve_alloc();
  EH(error, "ve_alloc problem");

  error = elc_alloc();
  EH(error, "elc_alloc problem");

  error = elc_rs_alloc();
  EH(error, "elc_alloc problem");

  error = cr_alloc();
  EH(error, "cr_alloc problem");

  error = evp_alloc();
  EH(error, "evp_alloc problem");

  error = tran_alloc();
  EH(error, "tran_alloc problem");

  error = libio_alloc();
  EH(error, "libio_alloc problem");

  error = eigen_alloc();
  EH(error, "eigen_alloc problem");

  error = cont_alloc();
  EH(error, "cont_alloc problem");

  error = loca_alloc();
  EH(error, "loca_alloc problem");

  error = efv_alloc();
  EH(error, "efv_alloc problem");

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier before read_input_file()\n", ProcID);
#ifdef PARALLEL
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

/*PRS AGAIN, NO COMMAND LINE OVERRIDES IN THIS JAS3D VERSION */
  /*
   * Read ASCII input file, data files, related exodusII FEM databases.
   */	
   if ( ProcID == 0 ) 
       { 
         log_msg("Reading input file ..."); 
         read_input_file(clc, nclc); 

       }

  /*
   * The user-defined material properties, etc. available to goma users
   * mean that some dynamically allocated data needs to be communicated.
   *
   * To handle this, sizing information from the input file scan is
   * broadcast in stages so that the other processors can allocate space
   * accordingly to hold the data.
   *
   * Note: instead of handpacking a data structure, use MPI derived datatypes
   * to gather and scatter. Pray this is done efficiently. Certainly it costs
   * less from a memory standpoint.
   */

#ifdef PARALLEL

  /*
   *  Make sure the input file was successully processed before moving on
   */
  check_parallel_error("Input file error");


  /*
   * This is some sizing information that helps fit a little bit more
   * onto the ark later on.
   */

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier before noahs_raven()\n", ProcID);
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif

  noahs_raven();

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier before MPI_Bcast of Noahs_Raven\n", ProcID);
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif

  MPI_Bcast(MPI_BOTTOM, 1, Noahs_Raven->new_type, 0, MPI_COMM_WORLD);

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier after Bcast/before raven_landing()\n", 
	  ProcID);
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif  
  /*
   * Get the other processors ready to handle ark data.
   */

  raven_landing();

#ifdef DEBUG
  fprintf(stderr, "P_%d at barrier before noahs_ark()\n", ProcID);
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  
  /*
   * This is the main body of communicated information, including some
   * whose sizes were determined because of advanced legwork by the raven.
   */

  noahs_ark();
  MPI_Bcast(MPI_BOTTOM, 1, Noahs_Ark->new_type, 0, MPI_COMM_WORLD);

  /*
   * Chemkin was initialized on processor zero during the input file
   * process. Now, distribute it to all processors
   */
#ifdef USE_CHEMKIN
  if (Chemkin_Needed) {
    chemkin_initialize_mp();
  }
#endif 

  /*
   * Once the ark has landed, there are additional things that will need to
   * be sent by dove. Example: BC_Types[]->u-BC arrays.
   *
   */

  ark_landing();

  noahs_dove();
  MPI_Bcast(MPI_BOTTOM, 1, Noahs_Dove->new_type, 0, MPI_COMM_WORLD);


#endif          /* End of ifdef PARALLEL */


  /*
   * We sent the packed line to all processors that contained geometry
   * creation commands.  Now we need to step through it and create
   * geometry as we go (including possibly reading an ACIS .sat file).
   *
   */
#ifdef USE_CGM
  create_cgm_geometry();
#endif

  /*
   * For parallel execution, assume the following variables will be changed
   * to reflect the multiple file aspect of the problem.
   *
   *	FEM file = file.exoII		--> file_3of15.exoII
   *
   *	Output EXODUS II file = out.exoII --> out_3of15.exoII
   *
   */


  /*
   * Allocate space for structures holding the EXODUS II finite element
   * database information and for the Distributed Processing information.
   *
   * These are mostly skeletons with pointers that get allocated in the
   * rd_exoII and rd_dpi routines. Remember to free up those arrays first
   * before freeing the major pointers.
   */

  EXO_ptr = alloc_struct_1(Exo_DB, 1);
  init_exo_struct(EXO_ptr);
  DPI_ptr = alloc_struct_1(Dpi, 1);
  init_dpi_struct(DPI_ptr);  

  log_msg("Reading mesh from EXODUS II file...");
  error = read_mesh_exoII(EXO_ptr, DPI_ptr);

  /*
   *   Missing files on any processor are detected at a lower level
   *   forcing a return to the higher level
   *         rd_exo -->  rd_mesh  -->  main
   *   Shutdown now, if any of the exodus files weren't found
   */
  if (error < 0) {
#ifdef PARALLEL
    MPI_Finalize();
#endif
    return(-1);
  }

  /*
   * All of the MPI_Type_commit() calls called behind the scenes that build
   * the dove, ark and raven really allocated memory. Let's free it up now that
   * the initial information has been communicated.
   */

#ifdef PARALLEL
  MPI_Type_free(&(Noahs_Raven->new_type));
  MPI_Type_free(&(Noahs_Ark->new_type));
  MPI_Type_free(&(Noahs_Dove->new_type));
#endif   

  /*
   * Setup the rest of the Problem Description structure that depends on
   * the mesh that was read in from the EXODUS II file...
   * 
   * Note that memory allocation and some setup has already been performed
   * in mm_input()...
   */

  error = setup_pd();
  EH( error, "Problem setting up Problem_Description.");
  /*
   * Let's check to see if we need the large elasto-plastic global tensors
   * and allocate them if so 
   */
  error = evp_tensor_alloc(EXO_ptr);
  EH( error, "Problems setting up evp tensors");
  
  /*
   * Now that we know about what kind of problem we're solving and the
   * mesh information, let's allocate space for elemental assembly structures
   *
   */
#ifdef DEBUG
  DPRINTF(stderr, "About to assembly_alloc()...\n");
#endif
  log_msg("Assembly allocation...");

  error = assembly_alloc(EXO_ptr);
  EH( error, "Problem from assembly_alloc");

  if (Debug_Flag)  {
    DPRINTF(stderr, "%s:  setting up EXODUS II output files...\n", yo);
  }

  /*
   * These are not critical - just niceties. Also, they should not overburden
   * your db with too much of this - they're capped verbiage compliant routines.
   */

  add_qa_stamp(EXO_ptr);

  add_info_stamp(EXO_ptr);

#ifdef DEBUG
  fprintf(stderr, "added qa and info stamps\n");
#endif

  /*
   * If the output EXODUS II database file is different from the input
   * file, then we'll need to replicate all the basic mesh information.
   * But, remember that if we're parallel, that the output file names must
   * be multiplexed first...
   */
  if ( Num_Proc > 1 )
    {
      multiname(ExoFileOut,     ProcID, Num_Proc);      
      multiname(Init_GuessFile, ProcID, Num_Proc);

      if ( strcmp( Soln_OutFile, "" ) != 0 )
	{
	  multiname(Soln_OutFile,   ProcID, Num_Proc);
	}

      if( strcmp( ExoAuxFile, "" ) != 0 )
        {
          multiname(ExoAuxFile,     ProcID, Num_Proc);
        }

      if( efv->Num_external_field != 0 )
        {
          for( i=0; i<efv->Num_external_field; i++ )
            {
              multiname(efv->file_nm[i], ProcID, Num_Proc);
            }
        }

    }



  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *   Preprocess the exodus mesh
   *        -> Allocate pointers to structures containing element
   *           side bc info, First_Elem_Side_BC_Array, and
   *           element edge info, First_Elem_Edge_BC_Array.
   *        -> Determine Unique_Element_Types[] array
   */
#ifdef DEBUG
  fprintf(stderr, "pre_process()...\n");
#endif
  log_msg("Pre processing of mesh...");
#ifdef PARALLEL
  error = MPI_Barrier(MPI_COMM_WORLD);
#endif
  pre_process(EXO_ptr);

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   * Load up a few key indeces in the bfd prototype basis function structures
   * and make sure that each active eqn/vbl has a bf[v] that points to the
   * right bfd[]...needs pre_process to find out the number of unique
   * element types in the problem.
   */

#ifdef DEBUG
  fprintf(stderr, "bf_init()...\n");
#endif
  log_msg("Basis function initialization...");
  error = bf_init(EXO_ptr);
  EH( error, "Problem from bf_init");

  /*
   * check for parallel errors before continuing
   */
  check_parallel_error("Error encountered in problem setup");

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/  
  /*
   * Allocate space for each communication exchange description.
   */
#ifdef PARALLEL
#ifdef DEBUG
  fprintf(stderr, "P_%d: Parallel cx allocation\n", ProcID);
#endif
  if (DPI_ptr->num_neighbors > 0) {
    cx = alloc_struct_1(Comm_Ex, DPI_ptr->num_neighbors);
    Request = alloc_struct_1(MPI_Request, 
			     Num_Requests * DPI_ptr->num_neighbors);
    Status = alloc_struct_1(MPI_Status, 
			    Num_Requests * DPI_ptr->num_neighbors);
  }
#endif

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *                           SET UP THE PROBLEM
   *
   * Setup node-based structures
   * Finalise how boundary conditions are to be handled
   * Determine what unknowns are at each owned node and then tell
   *  neighboring processors about your nodes
   * Set up communications pattern for fast unknown updates between
   *  processors.
   */
  (void) setup_problem(EXO_ptr, DPI_ptr);

  /*
   * check for parallel errors before continuing
   */
  check_parallel_error("Error encountered in problem setup");
  
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *                     WRITE OUT INITIAL INFO TO EXODUS FILE
   */

  /*
   *  Only have to initialize the exodus file if we are using different
   *  files for the output versus the input mesh
   */
  if (strcmp(ExoFile, ExoFileOut)) {
    /*
     * Temporarily we'll need to renumber the nodes and elements in the
     * mesh to be 1-based. After writing, return to the 0 based indexing
     * that is more convenient in C.
     */
#ifdef DEBUG
    fprintf(stderr, "1-base; wr_mesh; 0-base\n");
#endif
    one_base(EXO_ptr);
    wr_mesh_exo(EXO_ptr, ExoFileOut, 0);
    zero_base(EXO_ptr);

    /*
     * If running on a distributed computer, augment the plain finite
     * element information of EXODUS with the description of how this
     * piece fits into the global problem.
     */
    if (Num_Proc > 1) {
#ifdef PARALLEL
#ifdef DEBUG
      fprintf(stderr, "P_%d at barrier before wr_dpi()\n", ProcID);
      fprintf(stderr, "P_%d ExoFileOut = \"%s\"\n", ProcID, ExoFileOut);
      error = MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
      wr_dpi(DPI_ptr, ExoFileOut, 0);
    }
  }

  if (Num_Import_NV > 0 || Num_Import_EV > 0) printf
    (" Goma will import %d nodal and %d element variables.\n",
     Num_Import_NV, Num_Import_EV);
  if (Num_Export_XS > 0 || Num_Export_XP > 0) printf
    (" Goma will export %d solution and %d post-processing variables.\n",
     Num_Export_XS, Num_Export_XP);

  /* Return counts to calling program */
  *nnodes = EXO_ptr->num_nodes;
  *nelems = EXO_ptr->num_elems;
  *nnv_in = Num_Import_NV;
  *nev_in = Num_Import_EV;
  *i_soln = Num_Export_XS;
  *i_post = Num_Export_XP;

  return (0); /* Back to  animas*/
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/**************************GOMA SOLVE PROBLEM****************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

/* goma_solve_() */
int
goma_solve_(dbl *t_start_in,
            dbl *t_end_in,
            int *gfirst,
            int *ig,
            int *nprint,
            int *ngstep,
            dbl *previous_step,
            dbl  xsoln[],
            dbl  xpost[],
            dbl  xnv_in[],
            dbl  xev_in[])

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *                           SOLVE THE PROBLEM
   */
{

  /*declare some variables here */

  static const char *yo="goma_solve";
  double ts_in, te_in, old_dt;
  double te_out;
  double *p_te_out;
  int animas_step, print_flag, solve_steps, goma_first;
  int error;
int j;


/* Testing stuff....*/
printf("   goma_solve called for time %g to %g\n", *t_start_in, *t_end_in);
printf("   ANIMAS Step %3d:  Input times are %g to %g\n",
*ig, tran->init_time, tran->TimeMax);
printf("   Last ANIMAS step was %g\n", *previous_step);
printf("    Goma array check:\n");
printf(" Export var IDs:  %d  %d  %d\n",
Export_XS_ID[0], Export_XS_ID[1], Export_XS_ID[2]);
/*
printf("     NSOLN = %6d %6d %6d %6d %6d\n",
nsoln[0],nsoln[1],nsoln[2],nsoln[3],nsoln[4]);
printf("     NPOST = %6d %6d %6d %6d %6d\n",
npost[0],npost[1],npost[2],npost[3],npost[4]);
printf(" ..Nodes = %6d\n",EXO_ptr->num_nodes);
*/

  /*
   * This indicates that displacements are solved in another code
   * and imported to Goma for local mesh annealing.
   */
  efv->ev_porous_decouple = TRUE;
  animas_step = *ig;
  print_flag = *nprint;
  solve_steps = *ngstep;
  goma_first = *gfirst;
  p_te_out = &te_out;
  te_out = -1.0;

  /*
   * Save the start and end times passed in, then override the
   * input values in the "tran" structure.
   */
      ts_in = *t_start_in;
      te_in = *t_end_in;
      old_dt = *previous_step;
      tran->init_time = ts_in;
      tran->TimeMax = te_in;
printf("   Times passed in are  %g to %g\n", ts_in, te_in);

  /*
   * Load up libio structure with inputs from driver.
   */
  libio->xnv_in = xnv_in;
  libio->xev_in = xev_in;
  libio->xsoln = xsoln;
  libio->xpost = xpost;
  libio->animas_step = animas_step;
  libio->print_flag = print_flag;
  libio->goma_first = goma_first;
  libio->solve_steps = solve_steps;
  libio->t_start = ts_in;
  libio->t_end = te_in;
  libio->last_step = old_dt;

 /*
  * Set Exodus print intervals according to print_flag:
  * -1 = Never print
  *  0 = Print only at end
  *  1 = Print once per Goma call
  *  2 = Print once per Goma time step
  */
  if (goma_first == 1 && print_flag > 0) Write_Initial_Solution = TRUE;
  if (goma_first == 1 && print_flag == 0
                      && animas_step == -1) Write_Initial_Solution = TRUE;
  if (goma_first == 0) Write_Initial_Solution = FALSE;

  if (Debug_Flag) {
    switch (Continuation) {
    case ALC_ZEROTH:
        P0PRINTF("%s: continue_problem (zeroth order) ...\n", yo);
        break;
    case  ALC_FIRST:
        P0PRINTF("%s: continue_problem (first order) ...\n", yo);
        break;
    case HUN_ZEROTH:
        P0PRINTF("%s: hunt_problem (zeroth order) ...\n", yo);
        break;
    case  HUN_FIRST:
        P0PRINTF("%s: hunt_problem (first order) ...\n", yo);
        break;
    case LOCA:
        P0PRINTF("%s: do_loca ...\n", yo);
        break;
    default:
        P0PRINTF("%s: solve_problem...\n", yo);
        break;
    }
  }
#ifdef DEBUG
  switch (Continuation) {
  case ALC_ZEROTH:
      DPRINTF(stderr, "%s: continue_problem (zeroth order) ...\n", yo);
      break;
  case  ALC_FIRST:
      DPRINTF(stderr, "%s: continue_problem (first order) ...\n", yo);
      break;
  case HUN_ZEROTH:
      DPRINTF(stderr, "%s: hunt_problem (zeroth order) ...\n", yo);
      break;
  case  HUN_FIRST:
      DPRINTF(stderr, "%s: hunt_problem (first order) ...\n", yo);
      break;
  case LOCA:
      DPRINTF(stderr, "%s: do_loca ...\n", yo);
      break;
  default:
      DPRINTF(stderr, "%s: solve_problem...\n", yo);
      break;
  }
  fprintf(stderr, "P_%d: Before entering solve problem\n", ProcID);
#endif

  switch (Continuation) {
  case ALC_ZEROTH:
  case ALC_FIRST:
    log_msg("Solving continuation problem");
    continue_problem(cx, EXO_ptr, DPI_ptr);
    break;
  case HUN_ZEROTH:
  case HUN_FIRST:
    log_msg("Solving hunt problem");
    hunt_problem(cx, EXO_ptr, DPI_ptr);
    break;
  case LOCA:
    log_msg("Solving continuation problem with LOCA");
    error = do_loca(cx, EXO_ptr, DPI_ptr);
    break;
  default:
    log_msg("Solving problem");
    if (loca_in->Cont_Alg == LOCA_LSA_ONLY)
      {
        error = do_loca(cx, EXO_ptr, DPI_ptr);
      }
    else
      {
        solve_problem(EXO_ptr, DPI_ptr, p_te_out);
      }
    break;
  }

  /* Check actual vs. requested end time - warning if they differ */
  if ( (*t_end_in - te_out) > DBL_SMALL && solve_steps == 0)
    {
      P0PRINTF("Requested end time %g not reached - actual end time = %g!\n",
                *t_end_in, te_out);
    }
  *t_end_in = te_out;
  
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *  PRINT A MESSAGE TO STDOUT SAYING WE ARE DONE
   */
  P0PRINTF("\n-done\n\n");

  return(0);
}

/*************************************************************************/
/*****************************GOMA CLOSE**********************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

int goma_close_(dbl *time1)

{

  /*declare some variables here */
  /* Local Declarations */

     double time_start=*time1, total_time;   /* timing variables */
     static const char *yo="goma_solve";
#ifndef PARALLEL
     struct tm *tm_ptr;               /* additional serial timing variables */
     time_t the_time;
#endif
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *       FREE MEMORY ALLOCATED BY THE PROGRAM
   */
  /*
   * free the element block / element based structures
   */
  free_element_blocks();

  /*
   * free nodal based structures
   */
  free_nodes();

  /*
   * Free command line stuff
   */
   safer_free( ( void ** ) Argv );
/*   if ( ProcID == 0 ) */
/*     { */
/*       if ( argc > 1 )  */
/* 	{ */
/* 	  for (i=0; i<argc; i++) */
/* 	    { */
/* #ifdef DEBUG */
/* 	      fprintf(stderr, "clc[%d]->string &= 0x%x\n", i, clc[i]->string); */
/* 	      fprintf(stderr, "clc[%d]         &= 0x%x\n", i, clc[i]); */
/* #endif */
/* 	      safer_free((void **) &(clc[i]->string)); */
/* 	      safer_free((void **) (clc + i)); */
/* 	    } */
/* 	  safer_free((void **) &clc); */
/* 	} */
/*     } */

  /*
   * Free exodus database structures
   */
  free_exo(EXO_ptr);
  safer_free((void **) &EXO_ptr);

  if ( Num_Proc > 1 )
  {
    free_dpi(DPI_ptr);
  }
  else
  {
    free_dpi_uni(DPI_ptr);
  }

  safer_free((void **) &DPI_ptr);

  /*
   * Remove front scratch file [/tmp/lu.'pid'.0]
   */
  if (Linear_Solver == FRONT) 	
    {
  unlerr = unlink(front_scratch_directory);
  EH(unlerr, "Unlink problem with front scratch file");
    }


#ifdef PARALLEL
  total_time = ( MPI_Wtime() - time_start )/ 60. ;
  DPRINTF(stderr, "\nProc 0 runtime: %10.2f Minutes.\n\n",total_time);
  MPI_Finalize();
#else
  (void) time(&the_time);
  tm_ptr = gmtime(&the_time);
  total_time =  (double)  (   ( tm_ptr->tm_sec
                 + 60. * (   60. * ( tm_ptr->tm_yday * 24. + tm_ptr->tm_hour )
                                                           + tm_ptr->tm_min   )
                 - time_start ) / 60.
                          );
  fprintf(stderr, "\nProc 0 runtime: %10.2f Minutes.\n\n",total_time);
#endif  
  log_msg("GOMA ends normally.");
  return (0);
}

void
initialize_goma_export_vars_ ( int *nnodes,
				int *nelems,
				int *nnv_in,
				int *nev_in,
				int *nnod_r,
				int *i_soln,
				double *xnv_in,
				double *xev_in,
				double *xsoln )
{

  double *x;

  int ev;
  int I;
  int ie;
  int offset;

  asdv(&x, NumUnknowns);

  init_vec(x, cx, EXO_ptr, DPI_ptr, NULL, 0);

  /* 
   * Goma export variables.  may be read from material files or from preexisting exodus file
   */


  for( ev=0; ev< *i_soln; ev++)
    {
      offset = *nnod_r*ev;
      for( I=0; I< *nnod_r; I++ )
	{
	  ie = Index_Solution( I,Export_XS_ID[ev], 0, 0, -1);

	  if( ie != -1 ) xsoln[offset+I] = x[ie];
	}
    }

  safer_free( (void **) &x );
  return;
} 
	  



/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void
print_code_version(void)

    /*
     * Print the code version to standard out
     */
{
  printf("%s\n", VERSION);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
