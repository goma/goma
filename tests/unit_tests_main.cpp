#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"

extern "C" {
#include <stdio.h>
#include <limits.h>

#include "dpi.h"
#include "dp_types.h"
#include "exo_struct.h"
#include "rf_io_const.h"
#include "mm_elem_block_structs.h"
}

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


ELEM_BLK_STRUCT *Element_Blocks = NULL;       /* Pointer to array of global
                                                 element block information. */

ELEM_BLK_STRUCT *Current_EB_ptr = NULL; /* Pointer to the current element block
                                           structure. This is calculated in the
                                           fill */
