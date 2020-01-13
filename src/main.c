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
#include <limits.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


#ifdef HAVE_SUNMATH
#include <sunmath.h>
#endif

#ifdef FP_EXCEPT
#define __USE_GNU
#include <fenv.h>
#endif


#ifdef PARALLEL
#include "az_aztec.h"
#endif
extern void handle_ieee(void );

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h" 
#endif
#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_particles.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_update_parameter.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_rotate.h"
#include "bc_special.h"
#include "bc_surfacedomain.h"
#include "brk_utils.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_elem_block.h"
#include "mm_elem_block_structs.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_solve_segregated.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "std.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"



/*
 * Global variables defined here.
 */

Dpi *DPI_ptr = NULL;

Exo_DB *EXO_ptr = NULL;

#define LAST_LEGAL_STRING "Last Legal String"

char	Input_File[MAX_FNL]="\0";  /* input EXODUS II database w/ problem defn */

char	ExoFile[MAX_FNL]="\0";	   /* input EXODUS II database w/ problem defn */

char	Exo_LB_File[MAX_FNL]="\0"; /* EXODUS II load balance info for mesh */

char	ExoFileOut[MAX_FNL]="\0";  /* output EXODUS II database w/ results */

char	ExoFileOutMono[MAX_FNL]="\0";  /* output EXODUS II database without per proc identifier */

char	ExoAuxFile[MAX_FNL]="\0";  /* auxiliary EXODUS II database for initguess */

int     ExoTimePlane = INT_MAX;    /* Time plane # or continuation # of soln to use as an initguess */

char	Echo_Input_File[MAX_FNL]="\0";	/* echo of problem def file  */

char    Brk_File[MAX_FNL]="\0";        /* input file for brk called as subroutine */

char    DomainMappingFile[MAX_FNL]="\0"; /* Domain Mapping file. Maps the materials
				       and names of material boundaries
				       specified in this file into this file
				       into chemkin domains and chemkin
				       surface and volumetric domains. */

int     CPU_word_size;
int     IO_word_size;		/*Precision variables for exodus II files*/

char	Init_GuessFile[MAX_FNL];/* ASCII file holding initial guess */

char	Soln_OutFile[MAX_FNL];	/* ASCII file holding solution vector, */
				/* same format as Init_GuessFile      */

int	Debug_Flag;		/* Flag to specify debug info is to be     */
				/* printed out. The value of this flag     */
				/* determines the level of diagnostic info */
				/* which is printed to stdout              */
				/* Debug_Flag == 0 	No output          */
			   	/*	 1	minimun output             */
				/*	 2	medium  output		   */
				/*	 3	maximum output             */
			   	/*	 -1	check jacobian             */
				/*	 -2	check jacobian with scaling*/

int	New_Parser_Flag;	/* New_Parser_Flag = 0	Parse with old parser */
				/* New_Parser_Flag = 1  Parse with new flex/bison parser */


#ifdef MATRIX_DUMP
int     Number_Jac_Dump = 0;    /* Number of jacobians to dump out
				 * If the value is negative, then the one
				 * jacobian, the -n'th jacobian, is dumped
				 * out */
#endif
int	Iout;			/* Flag to specify level of diagnostic output */
				/* which is to be printed out for the program */

int Write_Intermediate_Solutions = FALSE; /* Flag specifies whether to */
	      			          /* write out solution data at each */
                                          /* Newton iteration. */
int     Write_Initial_Solution = FALSE;
				/* Flag to indicate whether to write the
			         * initial solution to the ascii and exodus
				 * output files */
int     Num_Var_Init ;		/* number of variables to overwrite with
				 * global initialization */
int     Num_Var_LS_Init;        /* number of variables to overwirte with
				 * level set index initialization */
int     Num_Var_External;	/* number of total external variables (exoII or pixel)*/
int     Num_Var_External_pix;	/* number of external variables (pixel only)*/
int     Anneal_Mesh;            /* flag specifying creation of a special exodus
				 *  file with coordinates adjusted to the
				 * deformed coordinates (i.e. new displacements
				 * are set to zero and mesh is deemed stress-free */

double Porous_liq_inventory; /*global variable for finite-insult boundary condition*/
double **Spec_source_inventory; /*global variable for cumulative reacted source */
double *Spec_source_lumped_mass; /*global variable for species lumped mass */

const char anneal_file[] = ANNEAL_FILE_NAME;


/*
 * Benner's frontal solver wants to strcat() onto this directory name,
 * so leave enough space for dirname/lu.123456.0, for example.
 *
 * Here it is defined and initialized.
 */

char front_scratch_directory[MAX_FNL] = FRONT_SCRATCH_DIRECTORY;

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
  "/************************************************************************ * \n",
  "* Goma - Multiphysics finite element software.                            * \n",
  "* Sandia National Laboratories.                                           * \n",
  "*                                                                         * \n",
  "* Copyright (c) 2014  Sandia Corporation.                                 * \n",
  "*                                                                         * \n",
  "* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  * \n",
  "* the U.S. Government retains certain rights in this software.            * \n",
  "*                                                                         * \n",
  "* This software is distributed under the GNU General Public License.      * \n",
  "\\************************************************************************/  \n",
  "\n",
  "PRS, PAS, RRR, KSC, RAC, TAB, DRN, PLH, DAL, IDG, ACK, ACS, RBS, RRL, MMH, SRS\n",
  "HKM, RAR, EDW, PKN, SAR, EMB, KT, DSB, DSH ...\n"
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

Comm_Ex **cx = NULL;		/* communications info for ea neighbor proc */

double time_goma_started;	/* Save it here... */

/*
 * Declare these copies so that other routines may access them as global
 * external variables instead of passing them through the argument lists.
 */

char **Argv;
int Argc;

int
main(int argc, char **argv)
     
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
  /*  struct tm *tm_ptr;               additional serial timing variables */
  time_t now;
#endif

  int error;
  int i;
  int j;

  char	**ptmp;
  char *yo;

  struct Command_line_command **clc=NULL; /* point to command line structure */
  int           nclc = 0;		/* number of command line commands */

/********************** BEGIN EXECUTION ***************************************/

#ifdef FP_EXCEPT
  feenableexcept ((FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID));
#endif

/* assume number of commands is less than or equal to the number of 
 * arguments in the command line minus 1 (1st is program name) */

  /*
  *  Get the name of the executable, yo
  */
  yo = argv[0];

#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  time_start = MPI_Wtime();
#endif /* PARALLEL */
#ifndef PARALLEL
  (void)time(&now);
  time_start = (double)now;
#endif /* PARALLEL */

  time_goma_started = time_start;

  Argv = argv;

  Argc = argc;

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

  /*
   * Some initial stuff that only the master process does.
   */

  if ( ProcID == 0 )
    {
      if (argc > 1)
	{
	  log_msg("Preprocessing command line options.");
	  clc = (struct Command_line_command **) 
	    smalloc( argc * sizeof(struct Command_line_command *));
	  for (i=0; i<argc; i++)
	    {
	      clc[i] = (struct Command_line_command *) 
		smalloc(sizeof(struct Command_line_command));
	      clc[i]->type   = 0; /* initialize command line structure */
	      clc[i]->i_val  = 0;
	      clc[i]->r_val  = 0.;
	      clc[i]->string = (char *) 
		smalloc(MAX_COMMAND_LINE_LENGTH*sizeof(char));
	      for ( j=0; j<MAX_COMMAND_LINE_LENGTH; j++)
		{
		  clc[i]->string[j] = '\0';
		}
	    }
	}

      strcpy(Input_File, "input");
      strcpy(Echo_Input_File , "echo_input");

      if (argc > 1) translate_command_line(argc, argv, clc, &nclc);
	  	  
	  ECHO("OPEN", Echo_Input_File);
      
	  echo_command_line( argc, argv, Echo_Input_File );
      print_code_version();
      ptmp = legal_notice;
      while ( strcmp(*ptmp, LAST_LEGAL_STRING) != 0 )
	{
	  fprintf(stdout, "%s", *ptmp++);
	}
    }

  /*
   *  Allocate the uniform problem description structure and
   *  the problem description structures on all processors
   */
  error = pd_alloc();
  EH(error, "pd_alloc problem");


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

  error = eigen_alloc();
  EH(error, "eigen_alloc problem");

  error = cont_alloc();
  EH(error, "cont_alloc problem");

  error = loca_alloc();
  EH(error, "loca_alloc problem");

  error = efv_alloc();
  EH(error, "efv_alloc problem");


  /*
   * Read ASCII input file, data files, related exodusII FEM databases.
   */	
  if ( ProcID == 0 )
    {
      log_msg("Reading input file ...");
      read_input_file(clc, nclc); /* Read ascii input file get file names */

      /* update inputed data to account for command line arguments that
       * might override the input deck...
       */
      log_msg("Overriding any input file specs w/ any command line specs...");
      if (argc > 1) apply_command_line(clc, nclc);

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


  noahs_raven();


  MPI_Bcast(MPI_BOTTOM, 1, Noahs_Raven->new_type, 0, MPI_COMM_WORLD);

  /*
   * Get the other processors ready to handle ark data.
   */

  raven_landing();

  
  
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

  /* Check to see if BRK File option exists and if so check if file exits */
  if (Brk_Flag == 1) {
    check_for_brkfile(Brk_File);
  }
  check_parallel_error("Error encountered in check for brkfile");

  /* Now break the exodus files */
  if (Num_Proc > 1 && ProcID == 0 && Brk_Flag == 1) {
    call_brk();
  }
  check_parallel_error("Error in brking exodus files");
  MPI_Barrier(MPI_COMM_WORLD);

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
	  multiname(Soln_OutFile, ProcID, Num_Proc);
	}

      if( strcmp( ExoAuxFile, "" ) != 0 )
        {
          multiname(ExoAuxFile, ProcID, Num_Proc);
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
    cx = malloc(sizeof(Comm_Ex *) * upd->Total_Num_Matrices);
  if (DPI_ptr->num_neighbors > 0) {

    int imtrx;
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      cx[imtrx] = alloc_struct_1(Comm_Ex, DPI_ptr->num_neighbors);
      Request = alloc_struct_1(MPI_Request, 
                               Num_Requests * DPI_ptr->num_neighbors);
      Status = alloc_struct_1(MPI_Status, 
                              Num_Requests * DPI_ptr->num_neighbors);
    }
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
   *               CREATE BRK_FILE IF ONE DOES NOT EXIST
   *
   * If no Brk_File exists but the option was configured in the input or
   * optional command we create one now and exit from goma.
   */
  if ( Brk_Flag == 2 ) {
    write_brk_file(Brk_File, EXO_ptr);
    MPI_Finalize();
    exit(0);
  }
  
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
#endif
      wr_dpi(DPI_ptr, ExoFileOut);
    }
  }

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *                           SOLVE THE PROBLEM
   */

  if (upd->Total_Num_Matrices == 1){
     pg->imtrx = 0;

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

    
  if( TimeIntegration == TRANSIENT)
        {
        Continuation = ALC_NONE;
        if (Debug_Flag) {
          P0PRINTF("%s: solve_problem...TRANSIENT superceded Continuation...\n", yo);
          }
        solve_problem(EXO_ptr, DPI_ptr, NULL);
        }  

  switch (Continuation) {
  case ALC_ZEROTH:
  case ALC_FIRST:
    log_msg("Solving continuation problem");
    continue_problem(cx[0], EXO_ptr, DPI_ptr);
    break;
  case HUN_ZEROTH:
  case HUN_FIRST:
    log_msg("Solving hunt problem");
    hunt_problem(cx[0], EXO_ptr, DPI_ptr);
    break;
  case LOCA:
    log_msg("Solving continuation problem with LOCA");
    error = do_loca(cx[0], EXO_ptr, DPI_ptr);
    break;
  default:
    log_msg("Solving problem");
    if (loca_in->Cont_Alg == LOCA_LSA_ONLY)
      {
        error = do_loca(cx[0], EXO_ptr, DPI_ptr);
      }
    else if(TimeIntegration != TRANSIENT)
      {
        solve_problem(EXO_ptr, DPI_ptr, NULL);
      }
    break;
  }
  } else {/* End of if upd->Total_Num_Matrices === 1 */
        solve_problem_segregated(EXO_ptr, DPI_ptr, NULL);
  }

#ifdef PARALLEL
   MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (ProcID == 0 && Brk_Flag == 1 && Num_Proc > 1) {
    fix_output();
  }
  
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *  PRINT A MESSAGE TO STDOUT SAYING WE ARE DONE
   */
  P0PRINTF("\n-done\n\n");

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *       FREE MEMORY ALLOCATED BY THE PROGRAM
   */
  /*
   * free the element block / element based structures
   */
  free_element_blocks(EXO_ptr);

  /*
   * free nodal based structures
   */
  free_nodes();
#ifdef FREE_PROBLEM
  free_problem ( EXO_ptr, DPI_ptr );
#endif

  /*
   * Free command line stuff
   */
  if ( ProcID == 0 )
    {
      if ( argc > 1 ) 
	{
	  for (i=0; i<argc; i++)
	    {
	      safer_free((void **) &(clc[i]->string));
	      safer_free((void **) (clc + i));
	    }
	  safer_free((void **) &clc);
	}
    }

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
  WH(unlerr, "Unlink problem with front scratch file");
    }


#ifdef PARALLEL
  total_time = ( MPI_Wtime() - time_start )/ 60. ;
  DPRINTF(stdout, "\nProc 0 runtime: %10.2f Minutes.\n\n",total_time);
  MPI_Finalize();
#endif  
#ifndef PARALLEL
  (void)time(&now);
  total_time = (double)(now) - time_start;
  fprintf(stdout, "\nProc 0 runtime: %10.2f Minutes.\n\n",total_time/60);
#endif  
  fflush(stdout);
  fflush(stderr);
  log_msg("GOMA ends normally.");
  return (0);
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
  printf("%s\n", GOMA_VERSION);
}

void
echo_command_line( int argc, char *argv[], char *echo_file)
{
	int istr=0;
	char echo_string[MAX_CHAR_IN_INPUT]="\0";
	time_t start_time;
	char *time_string;
	
	start_time = time (NULL);
	time_string = asctime( localtime(&start_time) );
	
	SPF(echo_string,"Command line :");
	
	while( istr < argc )
	{
		SPF(endofstring(echo_string)," %s", argv[istr]);
		istr++;
	}
	
	ECHO(echo_string, echo_file);
	
	SPF(echo_string,"Run Time : %s",  time_string);
		
	SPF(endofstring(echo_string),"Version : %s", GOMA_VERSION);
	
	ECHO(echo_string,echo_file);
	
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
