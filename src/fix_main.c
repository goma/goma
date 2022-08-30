/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
extern void handle_ieee(void);

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h"
#endif
#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_stability_util.h"
#include "decomp_interface.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "exo_struct.h"
#ifdef GOMA_ENABLE_METIS
#include "metis_decomp.h"
#endif
#include "brkfix/fix.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_input.h"
#include "mm_prob_def.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_element_storage_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_solve.h"
#include "rf_solve_segregated.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "std.h"
#include "wr_dpi.h"
#include "wr_exo.h"

#ifdef GOMA_ENABLE_PETSC
#include <petscsys.h>
#endif

/*
 * Global variables defined here.
 */

Dpi *DPI_ptr = NULL;

Exo_DB *EXO_ptr = NULL;

#define LAST_LEGAL_STRING "Last Legal String"

char Input_File[MAX_FNL] = "\0"; /* input EXODUS II database w/ problem defn */

char ExoFile[MAX_FNL] = "\0"; /* input EXODUS II database w/ problem defn */

char Exo_LB_File[MAX_FNL] = "\0"; /* EXODUS II load balance info for mesh */

char ExoFileOut[MAX_FNL] = "\0"; /* output EXODUS II database w/ results */

char ExoFileOutMono[MAX_FNL] = "\0"; /* output EXODUS II database without per proc identifier */

char ExoAuxFile[MAX_FNL] = "\0"; /* auxiliary EXODUS II database for initguess */

int ExoTimePlane = INT_MAX; /* Time plane # or continuation # of soln to use as an initguess */

char Echo_Input_File[MAX_FNL] = "\0"; /* echo of problem def file  */

int Decompose_Flag = 1;
int Decompose_Type = 0;
int Skip_Fix = 0;

char *GomaPetscOptions = NULL;
int GomaPetscOptionsStrLen = 0;

char DomainMappingFile[MAX_FNL] = "\0"; /* Domain Mapping file. Maps the materials
                                      and names of material boundaries
                                      specified in this file into this file
                                      into chemkin domains and chemkin
                                      surface and volumetric domains. */

int CPU_word_size;
int IO_word_size; /*Precision variables for exodus II files*/

char Init_GuessFile[MAX_FNL]; /* ASCII file holding initial guess */

char Soln_OutFile[MAX_FNL]; /* ASCII file holding solution vector, */
                            /* same format as Init_GuessFile      */

int Debug_Flag; /* Flag to specify debug info is to be     */
                /* printed out. The value of this flag     */
                /* determines the level of diagnostic info */
                /* which is printed to stdout              */
                /* Debug_Flag == 0 	No output          */
                /*	 1	minimun output             */
                /*	 2	medium  output		   */
                /*	 3	maximum output             */
                /*	 -1	check jacobian             */
                /*	 -2	check jacobian with scaling*/

int New_Parser_Flag; /* New_Parser_Flag = 0	Parse with old parser */
                     /* New_Parser_Flag = 1  Parse with new flex/bison parser */

#ifdef MATRIX_DUMP
int Number_Jac_Dump = 0; /* Number of jacobians to dump out
                          * If the value is negative, then the one
                          * jacobian, the -n'th jacobian, is dumped
                          * out */
#endif
int Iout; /* Flag to specify level of diagnostic output */
          /* which is to be printed out for the program */

int Write_Intermediate_Solutions = FALSE; /* Flag specifies whether to */
                                          /* write out solution data at each */
                                          /* Newton iteration. */
int Write_Initial_Solution = FALSE;
/* Flag to indicate whether to write the
 * initial solution to the ascii and exodus
 * output files */
int Num_Var_Init;         /* number of variables to overwrite with
                           * global initialization */
int Num_Var_Bound;        /* number of variables to bound  */
int Num_Var_LS_Init;      /* number of variables to overwirte with
                           * level set index initialization */
int Num_Var_External;     /* number of total external variables (exoII or pixel)*/
int Num_Var_External_pix; /* number of external variables (pixel only)*/
int Anneal_Mesh;          /* flag specifying creation of a special exodus
                           *  file with coordinates adjusted to the
                           * deformed coordinates (i.e. new displacements
                           * are set to zero and mesh is deemed stress-free */

double Porous_liq_inventory;     /*global variable for finite-insult boundary condition*/
double **Spec_source_inventory;  /*global variable for cumulative reacted source */
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
/*
 * Variables defined here but everywhere else are extern via "rf_mp.h".
 */

int parallel_err = FALSE;
int parallel_err_global = FALSE;

int Unlimited_Output = TRUE; /* print limit flag */

int Num_Proc = 1;
int ProcID = 0;
int Dim = -1;

int unlerr;

MPI_Request *Request = NULL;
MPI_Status *Status = NULL;
MPI_Aint type_extent;
int Num_Requests = 6;

/*
 * Data structures for transporting input data to other processors...
 */

DDD Noahs_Raven; /* Stage 1 - preliminary sizing information */
DDD Noahs_Ark;   /* Stage 2 - big boatload - main body*/
DDD Noahs_Dove;  /* Stage 3 - doubly dynamic sized stuff */

Comm_Ex **cx = NULL; /* communications info for ea neighbor proc */

double time_goma_started; /* Save it here... */

/*
 * Declare these copies so that other routines may access them as global
 * external variables instead of passing them through the argument lists.
 */

char **Argv;

int Argc;

ELEM_BLK_STRUCT *Element_Blocks = NULL; /* Pointer to array of global
                                           element block information. */

ELEM_BLK_STRUCT *Current_EB_ptr = NULL; /* Pointer to the current element block
                                           structure. This is calculated in the
                                           fill */

static void print_usage(void) {
  DPRINTF(stdout, "usage: fix <exodus_output>.<numproc>.<single>\n\nexample: fix out.exoII.24.00");
  DPRINTF(stdout, "    example: fix out.exoII.24.00\n\n");
}

static void get_fix_info(char *filename, int *num_procs, char mono_name[MAX_FNL]) {

  char infile[MAX_FNL];
  strncpy(infile, filename, MAX_FNL - 1);
  char *tokens[100];

  char previous_token[MAX_FNL];
  char last_token[MAX_FNL];

  mono_name[0] = '\0';
  previous_token[0] = '\0';
  last_token[0] = '\0';

  char *parse_string = infile;

  char *token = strtok(parse_string, ".");
  int tk = 0;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncpy(last_token, token, MAX_FNL - 1);
  while (token) {
    GOMA_ASSERT_ALWAYS(tk < 100);
    tokens[tk++] = token;
    strncpy(previous_token, last_token, MAX_FNL - 1);
    strncpy(last_token, token, MAX_FNL - 1);
    token = strtok(NULL, ".");
  }
#pragma GCC diagnostic pop

  if (tk < 4) {
    GOMA_EH(GOMA_ERROR, "Expected to find a file with at least 4 parts, got %d from %s", tk,
            filename);
  }

  for (int i = 0; i < (tk - 2); i++) {
    strcat(mono_name, tokens[i]);
    if (i < (tk - 3)) {
      strcat(mono_name, ".");
    }
  }

  char *p_end;
  long num_proc_parse = strtol(previous_token, &p_end, 10);
  if (errno == ERANGE) {
    GOMA_EH(GOMA_ERROR, "Could not parse number of processors %s from %s", previous_token,
            filename);
  }

  *num_procs = num_proc_parse;
}

int main(int argc, char **argv) {
  dbl time_start, total_time;

#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  time_start = MPI_Wtime();
#endif /* PARALLEL */

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  ProcID = rank;
  Num_Proc = size;

  if (size > 1) {
    GOMA_EH(GOMA_ERROR, "Fix only supports running on one processor.");
  }

  if (rank == 0) {

    print_code_version();
    if (argc != 2) {
      print_usage();
      GOMA_EH(GOMA_ERROR, "Unknown number of arguments %d", argc - 1);
    }

    if (argc == 2) {
      if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage();
        goto fix_exit;
      }
    }

    int num_procs;
    char mono_name[MAX_FNL];

    get_fix_info(argv[1], &num_procs, mono_name);

    Num_Proc = num_procs;

    DPRINTF(stdout, "Fixing file %s with %d procs\n", mono_name, num_procs);
    fix_exo_file(num_procs, mono_name);
  }
  total_time = (MPI_Wtime() - time_start) / 60.;
  DPRINTF(stdout, "\nProc 0 runtime: %10.2f Minutes.\n\n", total_time);
fix_exit:
#ifdef PARALLEL
  MPI_Finalize();
#endif
  fflush(stdout);
  fflush(stderr);
  return (0);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void print_code_version(void)

/*
 * Print the code version to standard out
 */
{
  printf("fix (Goma %s)\n", GOMA_VERSION);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
