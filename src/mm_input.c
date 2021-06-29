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

/*
 *$Id: mm_input.c,v 5.36 2010-06-29 22:23:42 prschun Exp $
 */

#ifdef USE_RCSID
static char rcsid[] =
"$Id: mm_input.c,v 5.36 2010-06-29 22:23:42 prschun Exp $";
#endif

#define _XOPEN_SOURCE /* POSIX WEXITSTATUS */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp and strncasecmp moved here for POSIX.1 */
#include <math.h>
#include <unistd.h>

#include <ctype.h>		/* for toupper(), isspace() */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_solver.h"
#include "rf_mp.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_bc_const.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"

#include "mm_eh.h"

#include "mm_post_def.h"

#include "sl_util_structs.h"
#include "mm_input.h"

#define _MM_INPUT_C
#include "goma.h"

extern void print_code_version	/* main.c -  */
PROTO((void ));

/*
 * Hey! This is the *one* place where these are defined. All other locations
 * have a mm_mp_structs and mm_mp.h to declare what these are.
 */

int Num_Var_Init_Mat[MAX_NUMBER_MATLS];	/* number of variables to overwrite  *
					 * with material-specific            *
					 * initialization                    */

struct Variable_Initialization	Var_init[MAX_VARIABLE_TYPES + MAX_CONC];

struct Variable_Initialization	Var_init_mat[MAX_NUMBER_MATLS]
						[MAX_VARIABLE_TYPES + MAX_CONC];

struct Boundary_Condition *BC_Types;

struct Rotation_Specs *ROT_Types;

extern struct AC_Information *augc;

struct HC_Information *hunt;

extern struct Eigensolver_Info *eigen;

struct Continuation_Conditions *cpcc;

struct Continuation_Conditions *tpcc;

struct User_Continuation_Info *cpuc;

struct User_Continuation_Info *tpuc;

struct Level_Set_Data *ls;
struct Level_Set_Interface *lsi;
struct Phase_Function_Data *pfd;

static char aprepro_command[1024];

static Spfrtn sr;


/*
 * How to blurt out what we found.
 */

static const char eoformat[MAX_CHAR_IN_INPUT] = "%s = %s";

/*
 * What to look for each time...
 */

static char search_string[MAX_CHAR_IN_INPUT];

static char default_string[MAX_CHAR_IN_INPUT] = "(default)";

//static char specify_string[MAX_CHAR_IN_INPUT] = "          ";


static Strcpy_rtn strcpy_rtn;		/* Data type def'd in std.h */

int	run_aprepro=0;

static char current_mat_file_name[MAX_FNL];

#define NO_USER  NULL
#define NO_INPUT 0
#define SCALAR_INPUT 1
#define VECTOR_INPUT 3

#define stringup(a) { char *p; for( p=a; *p != '\0'; *p=toupper(*p), p++); }

/*************** R O U T I N E S   I N   T H E   F I L E ***********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 * read_input_file			void		main (main.c)
 *    look_for				void		read_input_file
 *    look_for_optional			void		read_input_file
 *    strip				void		read_input_file,
 *							look_for
 *    count_parameters			void		rd_mp_props,
 *							
 *    read_string			void		read_input_file, 
 *							look_for
 *
 *    rd_file_spec()			void		read_input_file
 *
 *    rd_chemrx_specs()			void		read_input_file
 *
 *    rd_genl_specs()			void		read_input_file
 *
 *    rd_timeint_specs()		void		read_input_file
 *
 *    rd_track_specs()                  void            read_input_file
 * 
 *    rd_hunt_specs()                   void            read_input_file
 * 
 *    rd_ac_specs()                     void            read_input_file
 * 
 *    rd_solver_specs()			void		read_input_file	
 * 
 *    rd_eigen_specs()  		void		read_input_file	
 *
 *    rd_eq_specs()			void		read_input_file	
 *
 */
static double parse_press_datum_input PROTO((const char *));
static void read_MAT_line PROTO((FILE *, int, char *));

/*
 *	Read the input file for FEM reacting flow code
 *	according to the tumi/shadi format (not
 *	written down anywhere) and set algorithm and problem
 *	dependent parameters accordingly.
 *
 *	Moving mesh modification: split this big routine into several
 * 	smaller routines to handle each different section of the input file.
 *
 *	Add a new section for a new way to specify equations and terms in the
 *	equations that are to be activated for the problem.
 *    
 *	Unfinished:  ways of reading in the material property specifications
 *		     in a nice, extensible manner. The new EXODUS IIv2 has
 *		     advertised capabilities for specifying material properties,
 *		     but we shall see what the actual implementation in FASTQ
 *		     turns out to be...
 *
 *	Author:			John N. Shadid, 1421
 *	Revised:		Fri Oct 29 07:21:28 MDT 1993 pasacki@sandia.gov
 */

void
read_input_file(struct Command_line_command **clc,
		int nclc)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  char	input[MAX_CHAR_IN_INPUT];
  FILE	*ifp; 
  char *echo = Echo_Input_File;
#ifdef DEBUG
  static const char yo[] = "read_input_file";
#endif

  /* 
   * BEGIN EXECUTION
   */
  nAC = 0; augc = NULL;
  
  if( !New_Parser_Flag )	
  {			
    if ( (ifp=fopen(Input_File,"r")) != NULL)
    {

      rd_file_specs(ifp, input );

      rd_genl_specs(ifp, input );

      memset( tran->use_var_norm, 0, sizeof(int)*9 );
      rd_timeint_specs(ifp, input );
      
      rd_levelset_specs(ifp, input );
      
      rd_elem_quality_specs(ifp, input );

      rd_track_specs(ifp, input);

      rd_hunt_specs(ifp, input);

      rd_solver_specs(ifp, input );

      rd_ac_specs(ifp, input);

      rd_eigen_specs(ifp, input);

      if(look_for_optional(ifp, "Geometry Specifications", input, '=') == 1)
 	  EH(-1, "CGM not supported, there should be no Geometry Specifications section.");

      rd_particle_specs(ifp, input);

      rd_bc_specs(ifp, input);
      
      rd_matl_blk_specs(ifp, input);

      rd_post_process_specs(ifp, input);
    
      fclose(ifp);
	  echo_compiler_settings();
	  ECHO( "CLOSE", echo);	

    }
    else
    {
      sprintf(Err_Msg, "Could not open input file \"%s\"\n", Input_File);
      EH(-1, Err_Msg);
    }
  }  			
  else		
  {  			
#ifdef NEW_PARSER
    parse_input_file();   /* Call lex/yacc parser */
    exit(0);		  /* exit(0) is temporary */
#else	/* NEW_PARSER */
    sr = sprintf(err_msg, "This Goma executable was not built with the -newp option.");
    EH(-1, err_msg);      /* The -newp was used but this executable wasn'g built for it, exit */
    exit(0);
#endif	/* NEW_PARSER */    
  }    			    
  return;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

int 
count_parameters(const char string[])

/*
 *	This routine counts words in the string input.
 *      Allowed whitespace includes the normal definition plus
 *      the equals sign and a comma.
 *      Leading and trailing whitespace is ok
 *
 *	Author:			P. R. Schunk, 1511
 *	Date:			2/15/95
 *
 *	Parameter list:
 *
 *	string == On output 'string' contains the same characters as on
 *		  input 
 */
{
  int  i, np, state;
  char ch;

  i = 0;
  state = np = 0;

  /* '#' and '$' denote start of comment */
  while((ch = string[i]) != '\0'  && ch != '\n' && ch != '#' && ch != '$') {
    i++;
    if(ch == ' '  || ch == '\t' || ch == ',' || ch == '\f' ||
       ch == '\r' || ch == '\v' || ch == '=' ) {
      state = 0;
    } else if (state == 0) {
      state = 1;
      ++np;
    }
  } 
  return(np);
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void 
look_for(FILE  *ifp,
	 const char  *string,
	 char  input[],
	 const char ch_term)
/*
	Scan the input file (reading in strings according to 'read_string(ifp,)'
	specifications) until the character pattern in 'string' is matched.
	Leaves the file pointer positioned after the ch_term character on
	the line that contains the matched string.
	   An ERROR exit is created if the character string is not found.

	Author:			Ray S. Tuminaro Div 1422 
	Date:			8/8/90
	revised:		10/2/90 John N. Shadid

	Parameter list:

	ifp    == pointer to file "input"
	string == contains string pattern to be matched.
	input  == buffer array to hold characters that are read in.
	ch_term  == termination character. When scanning a line of input
		  is read until either a newline, the 'ch' termination
		  character is read, or the end-of-file is read.
*/
{
  char *yo = "look_for ERROR exit: ";
  if (read_string(ifp,input,ch_term) == -1) {
    fprintf(stderr,"%sEOF found in input file while searching for:\n",
	    yo);
    fprintf(stderr,"%s\n",string);
    exit(-1);
  }
  strip(input);
  while (strcmp(input,string) != 0 ) {
    if (read_string(ifp,input,ch_term) == -1) {
      fprintf(stderr,"%sEOF found in input file while searching for:\n",
	      yo);
      fprintf(stderr,"%s\n",string);
      exit(-1);
    }
    strip(input);
  }
}
/**************************************************************************/

int
look_for_either(FILE *ifp,
		const char *string1,
		const char *string2,
		char input[],
		const char ch_term)
{
/*
	Scan the input file (reading in strings according to 'read_string(ifp,)'
	specifications) until the character pattern in 'string1' or 'string2' is matched.

	Author:			Ray S. Tuminaro Div 1422 
	Date:			8/8/90
	revised:		5/9/95 Richard A. Cairncross

	Parameter list:

	ifp    == pointer to file "input"
	string1 == contains primary string pattern to be matched.
	string1 == contains secondary string pattern to be matched.
	input  == buffer array to hold characters that are read in.
	ch_term  == termination character. When scanning a line of input
		  is read until either a newline, the 'ch' termination
		  character is read, or the end-of-file is read.

        returns:  1  if string1 is matched
                  2  if string2 is matched
                 -1  if neither string is found before EOF
*/

	if (read_string(ifp,input,ch_term) == -1) {
	   fprintf(stderr,"EOF found in input file while searching for:\n");
	   fprintf(stderr,"%s or %s\n",string1, string2);
	   return -1;
	}

        strip(input);
 
        while ( (strcmp(input,string1) != 0) && (strcmp(input,string2) != 0) ) {
	   if (read_string(ifp,input,ch_term) == -1) {
	      fprintf(stderr,"EOF found in input file while searching for:\n");
	      fprintf(stderr,"%s or %s\n",string1, string2);
	      return -1;
	   }
           strip(input);
        }
	if (strcmp(input,string1) == 0) return 1;
	if (strcmp(input,string2) == 0) return 2;
	return -1;
}
/**************************************************************************/

int
count_list(FILE  *ifp,	/* file pointer to open goma input file */
	   const char  *string,
	   char  input[],
	   const char   ch_term,
	   const char  *stringend)
{
/*
	Scan the input file (reading in strings according to 'read_string(ifp,)'
	specifications) counting the occurences of 'string' before the 
        termination string 'stringend' is matched.  Then set the file 
        pointer back to the starting point, and return the number of 
        occurences of 'string' found.

	Author:			Richard A. Cairncross Div 1511 
	Date:			5/9/95

	Parameter list:

	ifp    == pointer to file "input"
	string == contains primary string pattern to be matched.
	input  == buffer array to hold characters that are read in.
	ch_term  == termination character. When scanning a line of input
		  is read until either a newline, the 'ch' termination
		  character is read, or the end-of-file is read.
	stringend == contains secondary string pattern which indicates end of list.

        returns:  number of occurences of string before stringend
*/

  int count = 0, status = 0;
  fpos_t file_position; /* position in file at start of search */
  
#ifndef tflop
  fgetpos(ifp, &file_position);
#else
  file_position = ftell( ifp );
#endif

  while( (status = look_for_either(ifp, string, stringend, input, ch_term)) == 1)
    {
      count++;
    }

#ifdef DEBUG
  if (status == 2) fprintf(stderr, "List of %s contains %d items\n", string, count);
#endif

  if (status == -1) fprintf(stderr, "### EOF found while counting list of %s ###\n"
	       "You should use termination string \"%s\" to indicate end of list\n"
			    "Will use %d as number of %s\n"
			    , string, stringend, count, string);

/* restore file position to position at start of search */
#ifndef tflop
  fsetpos(ifp, &file_position);
#else
  fseek(ifp, file_position, SEEK_SET);
#endif

  return count;
}
/**************************************************************************/
int 
look_for_optional(FILE *ifp,
		  const char *string,
		  char input[],
		  const char ch_term)
/*
	Scan the input file (reading in strings according to 'read_string(ifp,)'
	specifications) until the character pattern in 'string' is matched.
	If 'string' is not matched print a note to the screen and return
	a value of -1. If 'string' is matched, return 1.

	This search automatically starts at the beginning of the input
	file and searches to the end until it finds the input string.
	If this routine can't find the input string, it returns to the 
        beginning of the file

	Author:			Ray S. Tuminaro Div 1422 
	Date:			8/8/90
	revised:		2/2/95 Rich Cairncross

	Parameter list:

	ifp    == pointer to file "input"
	string == contains string pattern to be matched.
	input  == buffer array to hold characters that are read in.
	          The length of the buffer array is MAX_CHAR_IN_INPUT,
	          whose current value is 256.
	ch_term  == termination character. When scanning a line of input
		  is read until either a newline, the 'ch' termination
		  character is read, or the end-of-file is read.
*/
{
  rewind(ifp);
  if (read_string(ifp, input, ch_term) == -1) {
    /*
     * fprintf(stderr,"Didn't find \"%s\"; defaulting.\n", string);
     */
    rewind(ifp);
    return(-1);
  }

  strip(input);
 
  while (strcmp(input, string) != 0 ) {
    if (read_string(ifp, input, ch_term) == -1) {
      /*
       * fprintf(stderr,"Didn't find \"%s\"; defaulting.\n", string);
       */
      rewind(ifp);
      return(-1);
    }
    strip(input);
  }

  return 1;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int 
look_forward_optional(FILE *ifp, const char *string, char input[],
		      const char ch_term)
   /***********************************************************************
    *
    * look_forward_optional():
    *
    * Scan the input file (reading in strings according to 
    * 'read_string(ifp,)' specifications) 
    * until the character pattern in 'string' is matched.
    * If 'string' is not matched, return a value of -1.
    * If 'string' is matched, return 1.
    *
    * This search starts at the current position in the input
    * file and searches to the end until it finds the input string.
    * If this routine can't find the input string, it returns to the 
    * file position where the search was started. If it finds the
    * string pattern, it leaves the file pointer positioned after the
    * ch_term character.
    *
    * Author:			Ray S. Tuminaro Div 1422 
    * Date:			8/8/90
    * revised:		4/24/95 Rich Cairncross
    *
    * Parameter list:
    *  ifp      == pointer to file "input"
    *  string   == contains string pattern to be matched.
    *  input    == buffer array to hold characters that are read in.
    *              The length of the buffer array is MAX_CHAR_IN_INPUT,
    *              whose current value is 256.
    *  ch_term  == termination character. When scanning a line of input
    *	           is read until either a newline, the 'ch' termination
    *	           character is read, or the end-of-file is read.
    * Return Value:
    *    1 if string is matched
    *   -1 if no match
    *********************************************************************/
{
  int status = 1;
  fpos_t file_position;  /* position in file at start of search */
#ifndef tflop
  fgetpos(ifp, &file_position);
#else
  file_position = ftell(ifp);
#endif

  /*
   *  Read the first line to get the while look initialized
   */
  if (read_string(ifp, input, ch_term) == -1) {
#ifndef tflop
    fsetpos(ifp, &file_position);
#else
    fseek(ifp, file_position, SEEK_SET);
#endif
    return(-1);
  }
  strip(input);
  /*
   * Main loop -> process the input file until a match is found
   */
  while (strcmp(input, string) != 0 ) {
    if (read_string(ifp, input, ch_term) == -1) {
#ifndef tflop
      fsetpos(ifp, &file_position);
#else
      fseek(ifp, file_position, SEEK_SET);
#endif
      return(-1);
    }
    strip(input);
  }
  return(status);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*
 * rd_file_spec -- read problem specification section of input file
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Someday we need to comb through all these placeholder options
 *		to see if they are useful.
 *
 *
 * Revised:			Fri Oct 29 06:36:00 MDT 1993 pasacki@sandia.gov
 *
 * Revised: 1997/05/19 10:59 MDT pasacki@sandia.gov
 */

void 
rd_file_specs(FILE *ifp,
	      char *input )
{
#ifdef DEBUG
  static const char yo[] = "rd_file_specs";
#endif
  int foundMappingFile;
  int foundBrkFile;
  
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;
  
 
  look_for_optional(ifp,"FEM File Specifications",input,'=');

  ECHO("\n***FEM File Specifications***\n", echo_file);
  
  look_for(ifp,"FEM file",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input);
  (void) strcpy(ExoFile,input);

  SPF(echo_string, eoformat, "FEM file", ExoFile); ECHO(echo_string, echo_file);
  
  /*
   * Here's a new feature:  use a different file for the EXODUS II output
   * just in case we need to muck with anything that exists in the input
   * EXODUS II file...
   */
  
  look_for(ifp,"Output EXODUS II file",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input);
  strcpy(ExoFileOut,input);
  strcpy(ExoFileOutMono, input);
  SPF(echo_string,eoformat, "Output EXODUS II file", ExoFileOut); ECHO(echo_string, echo_file);
  
  look_for(ifp,"GUESS file",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input);
  if (strcasecmp(input, "NONE") && strcasecmp(input, "NO")) {
    strcpy(Init_GuessFile,input);
  } else {
    Init_GuessFile[0] = '\0';
  }

  SPF(echo_string, eoformat,"GUESS file" , Init_GuessFile); ECHO(echo_string, echo_file);
  
  look_for(ifp,"SOLN file",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input); 
  if (strcasecmp(input, "NONE") && strcasecmp(input, "NO")) {
    strcpy(Soln_OutFile,input);
  } else {
    Soln_OutFile[0] = '\0';
  }

  SPF(echo_string,eoformat, "SOLN file", Soln_OutFile); ECHO(echo_string, echo_file);

  /* Find BRK file */

  foundBrkFile = look_for_optional(ifp,"Brk file",input,'=');

  if ( foundBrkFile == 1 && Brk_Flag != 1 ) {
    Brk_Flag = 1;
    read_string(ifp,input,'\n');
    strip(input);
    strcpy(Brk_File,input);

    SPF(echo_string, eoformat, "Brk file", Brk_File); ECHO(echo_string, echo_file);
  }

  /*
   *   look_for Optional Domain mapping file, the usage of the default
   *   will be indicated by the null character string in the name.
   */
  DomainMappingFile[0] = '\0';
  foundMappingFile = look_for_optional_string(ifp, "Domain Mapping File", 
					      DomainMappingFile, MAX_FNL);
  if ( foundMappingFile == 1 )
    { SPF(echo_string,eoformat, "Domain Mapping File", DomainMappingFile ); ECHO(echo_string, echo_file); }
  
  /*
   * New feature: write out solution at each Newton iteration.
   *	Answer: "yes" or "no".
   */

  look_for(ifp, "Write intermediate results", input, '=');
  (void) read_string(ifp, input, '\n');
  strip(input);
  if (strcasecmp(input, "no") == 0 ) {
    Write_Intermediate_Solutions = FALSE;
    ECHO("Write intermediate results = no", echo_file);
  } else if (strcasecmp(input, "yes") == 0) {
    Write_Intermediate_Solutions = TRUE;
    ECHO("Write intermediate results = yes", echo_file);
  } else {
    EH( -1, "Bad specification for intermediate results");
  }
  
   /*
    * New feature: Write out an Initial guess or an
    *              initial solution at the start of a
    *              time dependent problem.
    *
    *	Answer: "yes" or "no".
    */ 
  if (look_forward_optional(ifp, "Write initial solution", input, '=')
      == 1) {
    (void) read_string(ifp, input,'\n');
    strip(input);
    if (! strcasecmp(input, "no")) {
      Write_Initial_Solution = FALSE;
      ECHO("Write initial solution = no", echo_file);
    } else if (!strcasecmp(input, "yes")) {
      Write_Initial_Solution = TRUE;
      ECHO("Write initial solution = yes", echo_file);
    } else {
      EH( -1, "Bad specification for intermediate results");
    }
  }
  

}
/* rd_file_spec -- read problem specification section of input file */

/*
 * rd_genl_specs -- read input file for general specifications
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Someday we need to comb through all these placeholder options
 *		to see if they are useful.
 *		
 *		This handles general purpose directions that do not fit nicely
 *		into other categories.
 *		
 *
 * Revised:			Fri Oct 29 07:48:27 MDT 1993 pasacki@sandia.gov
 */

void 
rd_genl_specs(FILE *ifp,
	      char *input )
{
  int nargs;
  char err_msg[MAX_CHAR_IN_INPUT];
  char first_string[MAX_CHAR_IN_INPUT];
  char second_string[MAX_CHAR_IN_INPUT];
  char third_string[MAX_CHAR_IN_INPUT];
  char *tmp; 
  char StringToSearch[]="Pixel"; /*used in strstr call below*/

  static const char yo[] = "rd_genl_specs";
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;
  char ftype[MAX_CHAR_IN_INPUT];

  int iread;
  int ex_id=0;
  
  /* Read in General Specifications section */
    
  iread = look_for_optional(ifp,"General Specifications", input, '=');

  ECHO("\n***General Specifications***\n", echo_file);
    
  iread = look_for_optional(ifp,"Output Level",input,'=');
  if (iread == 1) {
    if (fscanf(ifp,"%d",&Iout) != 1)
      {
	fprintf (stderr, "%s:\tError reading Output Level\n", yo);
	exit (-1);
      }
  } else {
    Iout = 0;
  }

  SPF(echo_string,"%s = %d", "Output Level",Iout ); ECHO(echo_string, echo_file);   
 
  iread = look_for_optional(ifp,"Debug",input,'=');
  if (iread == 1) {
    if (fscanf(ifp, "%d", &Debug_Flag) != 1)
      {
	fprintf(stderr, "%s:\tError reading Debug Level\n", yo);
	exit (-1);
      }
  } else {
    Debug_Flag = 0;
  }

  SPF(echo_string,"%s = %d","Debug", Debug_Flag); ECHO(echo_string, echo_file);
  
#ifdef MATRIX_DUMP
  (void) look_for_optional_int(ifp, "Number of Jacobian File Dumps",
			       &Number_Jac_Dump, 0
);

  SPF(echo_string,"%s = %d","Number of Jacobian File Dumps", Number_Jac_Dump); ECHO(echo_string, echo_file);
#endif

  iread = look_for_optional(ifp,"Initial Guess",input,'=');
  if(iread == 1)
    {
      (void) read_string(ifp,input,'\n');
      strip(input);
      nargs = sscanf(input, "%s %s %s", first_string, second_string, third_string);
      if ( nargs == 0 )
	{
	  EH(-1, "Found zero arguments for Initial Guess");
	}
      else if ( nargs == 1 )
	{
	  if (strcasecmp(first_string, "zero") == 0 )
	    {
	      Guess_Flag = 0;
	    }
	  else if (strcasecmp(first_string, "random") == 0 )
	    {
	      Guess_Flag = 1;
	    }
	  else if (strcasecmp(first_string,"one") == 0 )
	    {
	      Guess_Flag = 2;
	    }
	  else if (strcasecmp(first_string,"sines") == 0 )
	    {
	      Guess_Flag = 3;
	    }
	  else if (strcasecmp(first_string, "read") == 0 )
	    {
	      Guess_Flag = 4;
	    }
	  else if (strcasecmp(first_string, "read_exoII") == 0 )
	    {
	      Guess_Flag = 5;
	    }
	  else if (strcasecmp(first_string, "read_exoII_file") == 0 )
	    {
	      EH(-1, "Read from *what* exoII file?");
	    }

	  SPF(echo_string,eoformat,"Initial Guess", first_string); ECHO(echo_string, echo_file);

	}
      else if ( nargs == 2 )
	{
	  if (strcasecmp(first_string, "read_exoII_file") == 0 )
	    {
	      Guess_Flag = 6;
	      strcpy(ExoAuxFile, second_string);
	    }
	  else
	    {
	      EH(-1, "Undecipherable 2 options for Initial guess.");
	    }

	  SPF(echo_string,"%s = %s %s","Initial Guess", first_string, second_string); ECHO(echo_string, echo_file);

	}
      else if ( nargs == 3 )
	{
	  if (strcasecmp(first_string, "read_exoII_file") == 0 )
	    {
	      Guess_Flag = 6;
	      strcpy(ExoAuxFile, second_string);
	      if (sscanf(third_string, "%d", &ExoTimePlane) != 1) {
		EH(-1, "Time plane for read_exoII_file option is undecipherable");
	      }
	    }
	  else
	    {
	      EH(-1, "Undecipherable first 2 options for Initial guess.");
	    }

	  SPF(echo_string,"%s = %s %s %d","Initial Guess", first_string, second_string, ExoTimePlane);
	  ECHO(echo_string, echo_file);
	    
	}
      else
	{
	  fprintf(stderr,"%s:\tUnknown initial guess (%s)\n", yo, input);
	  exit(-1);
	}
    } else {
      Guess_Flag = 0;
      ECHO("Initial Guess card not read correctly", echo_file);
    }

  iread = look_for_optional(ifp,"Conformation Map",input,'=');
  if(iread == 1) 
    {
      (void) read_string(ifp,input,'\n');
      strip(input);
      nargs = sscanf(input, "%s %s %s", first_string, second_string, third_string);
      if ( nargs == 0 )
        {
          EH(-1, "Found zero arguments for Conformation Map");
        }
      else if ( nargs == 1 )
        {
          if (strcasecmp(first_string, "no") == 0 )
            {
              Conformation_Flag = 0;
            }
          else if (strcasecmp(first_string, "yes") == 0 )
            {
              Conformation_Flag = 1;
            }

          SPF(echo_string,eoformat,"Conformation Map", first_string); ECHO(echo_string, echo_file);

        }
      else  
        {
          fprintf(stderr,"%s:\tUnknown conformation map (%s)\n", yo, input);
          exit(-1);
        }
    } else {
      Conformation_Flag = 0;
    }

  /*
   *             Search for commands to initialize a specific variable
   */
  Num_Var_Init = 0;
  while ((iread = look_forward_optional(ifp,"Initialize",input,'=')) == 1)
  {
    /*
     *  Read the variable name to be fixed
     */
    if (fscanf(ifp, "%80s", input) != 1)
    {
      EH(-1, "Error reading variable for initialization");
    }
    /*
     *  Translate the string variable name to the internal integer value for
     *  that variable.
     */
    Var_init[Num_Var_Init].var =
	variable_string_to_int(input, "Initialize Keyword Error");

    if (fscanf(ifp, "%d %lf", &Var_init[Num_Var_Init].ktype,
	       &Var_init[Num_Var_Init].init_val)
	!= 2) {
      EH(-1,"Error reading initialization data");
    }

    SPF(echo_string,"%s = %s %d %f","Initialize",input, Var_init[Num_Var_Init].ktype, Var_init[Num_Var_Init].init_val); ECHO(echo_string, echo_file);

    Num_Var_Init++;
  }

  /* Search for commands to read in and hold fixed an external nodal field variable */
  Num_Var_External = 0;
  Num_Var_External_pix = 0;
  Num_Import_NV = 0;
  Num_Import_EV = 0;
  efv->ev = T_NOTHING;

  memset(efv->ipix, 0, sizeof(int)*MAX_EXTERNAL_FIELD);
  
  while ((iread = look_forward_optional(ifp,"External Field",input,'=')) == 1 || 
	 (iread = look_forward_optional(ifp,"External Pixel Field",input,'='))	 == 1 ||
	 (iread = look_forward_optional(ifp,"External Pixel Field_fast",input,'=')) == 1)
    {
       efv->ev = T_SOMETHING;

       tmp = strstr(input,StringToSearch);
       if(tmp != NULL)
	 {
	   strcpy(echo_string, "External Pixel Field = ");
	   Num_Var_External_pix++;
	   efv->ipix[Num_Var_External] = 1;
	   if (strstr(input,"Field_fast") != NULL) {
	     strcpy(echo_string, "External Pixel Field_fast = ");
	     Num_Var_External_pix++;
	     efv->ipix[Num_Var_External] = 2;
	   }
	 }
       else
	 {
	    strcpy(echo_string, "External Field = ");
	 }
      /*
       *  Read the variable name to be looked for on the exodus file
       */
      if (fscanf(ifp, "%20s", input) != 1)
	{
	  EH(-1, "Error reading variable for initialization");
	}

      (void) strcpy(efv->name[Num_Var_External] , input); 

      strcat(echo_string,input); strcat(echo_string,"  ");

      /*
       *  Read the interpolation of the variable
       */
      
      if ( fscanf(ifp, "%s", input) !=1 )
	{
	  fprintf(stderr, "%s: problem reading interpolation for external variable.\n",yo);
	  exit(-1);
	}
      /*
       * Check to see if this interpolation is recognizable and if it
       * is reasonable.
       */

      if ( !strcmp(input, "Q1") )
	{
	  efv->i[Num_Var_External] = I_Q1;
	}
      else if ( !strcmp(input, "Q2") )
	{
	  efv->i[Num_Var_External] = I_Q2;
	}
      else if ( !strcmp(input, "Q2_LSA") )
	{
	  efv->i[Num_Var_External] = I_Q2_LSA;
	}
      else if ( !strcmp(input, "Q1_D") )
	{
	  efv->i[Num_Var_External] = I_Q1_D;
	}
      else if ( !strcmp(input, "Q2_D") )
	{
	  efv->i[Num_Var_External] = I_Q2_D;
	}
      else if ( !strcmp(input, "Q2_D_LSA") )
	{
	  efv->i[Num_Var_External] = I_Q2_D_LSA;
	}
      else if ( !strcmp(input, "PQ1") )
	{
	  efv->i[Num_Var_External] = I_PQ1;
	}
      else if ( !strcmp(input, "PQ2") )
	{
	  efv->i[Num_Var_External] = I_PQ2;
	}
      else if ( !strcmp(input, "P0") )
	{
	  efv->i[Num_Var_External] = I_P0;
	}
      else if ( !strcmp(input, "P1") )
	{
	  efv->i[Num_Var_External] = I_P1;
	}
      else if ( !strcmp(input, "SP") )
	{
	  efv->i[Num_Var_External] = I_SP;
	}
      else if ( !strcmp(input, "TABLE") )
	{
	  efv->i[Num_Var_External] = I_TABLE;
	  num_ext_Tables++;
	}
      else
	{
	  sprintf(err_msg, 
	          "??? interpolation \"%s\" for external field var",
	          input);
          EH(-1, err_msg);
	}

      strcat(echo_string,input);strcat(echo_string,"  ");

      /* 
       * read in exodus II file name
       */

      if (fscanf(ifp, "%s", input) != 1)
	{
	  EH( -1, "error reading exodus II file Name or voxel file name for external field.");
	}
      strip(input);

      (void) strcpy (efv->file_nm[Num_Var_External], input);

      strcat( echo_string,input);

      /*
       * Count external (nodal and elemental) vars to be imported
       */
      if ( !strcmp(input, "IMPORT_EV") ) Num_Import_EV++;
      if ( !strcmp(input, "IMPORT") ) Num_Import_NV++;

      /* Option to initialize external field to some value (default = 0) if
       * if any part of the element block being mapped to lies outside of the
       * voxel field. DSB 7/26/13
       */
      //if ((iread = look_forward_optional(ifp, "t", input, 'E')) == 1)

      /*
       * read in type of external field variables
       * 
       * SMD 1/24/11
       * Currently the last time step within an external field is read
       * and that value held constant throughout the solution.
       * The ability to read an external exoII file at time values that
       * match the GOMA solution process has been added.
       *
       * To initiate time dependent reading of external fields add 
       * "time_dep" to the end of the External Field = line in the input deck.
       * Absense of that phrase causes the code to default to the standard fixed
       * external field reading that existed previously.
       */
      if ((iread = look_forward_optional(ifp, "time_dep", input, 'E')) == 1)
	{
	  (void) strcpy(ftype, "transient");
	}	
      else
	{
	  (void) strcpy(ftype, "steady");
	}

      strip(input);
     
      (void) strcpy (efv->field_type[Num_Var_External],ftype);
      
      strcat(echo_string,ftype);

      /*
       * read in block id for external pixel field.  Defaults to 1. 
       * 
       * PRS 7/27/2011
       */
      if(Num_Var_External+1 > MAX_EXTERNAL_FIELD)
  {
    SPF(err_msg,
        ">%d external field vars. Fix MAX_EXTERNAL_FIELD"
        " (rf_fem_const.h), recompile.",
        MAX_EXTERNAL_FIELD);
          EH(-1,err_msg);
  }
      if ( efv->ipix[Num_Var_External])
	{
    if (fscanf (ifp,"%d",&efv->ipix_matid[Num_Var_External]) != 1)
	    {
	      EH(-1, "Must specify a material ID for external pix field");
	    }


	  /* Optional floating point value to set for the field at any
	   * Gauss points that are outside of the pixel/voxel field.
	   * Default to 1.0 (0 can cause issues for e.g. conductivities).
	   * DSB 7/30/2013
	   */
	  if (fscanf (ifp,"%lf",&efv->empty_value[Num_Var_External]) != 1){
	    efv->empty_value[Num_Var_External] = 1.0;
	  }
	}

    
      Num_Var_External++;

      ECHO(echo_string, echo_file);
    }

  efv->Num_external_field = Num_Var_External;
  efv->Num_external_pixel_field = Num_Var_External_pix;


  /*
   * Read export variable cards
   */
  Num_Export_XS = 0;
  while ((iread = look_forward_optional(ifp,"Export Field",input,'='))
	 == 1)
    {
      
    /*
     *  Read the variable name to be exported from Goma
     */
    if (fscanf(ifp, "%d", &ex_id) != 1)
    {
      EH(-1, "Error reading variable for export!");
    }

    /*
     * Ensure the id corresponds to a valid Goma variable
     */
    if (ex_id < 0 || ex_id >= V_LAST)
    {
      fprintf(stderr, "Export Field %d is not a valid variable ID!", ex_id );
      exit(-1);
    }

    /*
     * Enter the id in the global array Export_XS_ID.
     */
    Export_XS_ID[Num_Export_XS] = ex_id;
    Num_Export_XS++;
    
    SPF(echo_string,"%s = %d","Export Field",ex_id); ECHO(echo_string, echo_file);

    }
#ifndef LIBRARY_MODE
    if (Num_Export_XS > 0) WH(-1, "Export Field only valid in LIBRARY_MODE!");
#endif


  /* XFEM turned off by default */
  upd->XFEM = FALSE;
  
  /*
   *  Pressure Datum - Optional
   *        Set the pressure Datum for use in thermodynamic
   *        equations of state calculations.
   *        This is an additive constant that get added onto the
   *        pressure field before calculation of all thermodynamic
   *        equations of state calculations. It is a constant
   *        over the entire domain. Therefore, it is not included
   *        in any one materials file. The default units for
   *        the quantity are cgs units, and the default value
   *        for the quantity is 1 atmosphere = gm cm-1 sec-2 (dyne cm-2).
   *        (conversion factor is to an exact standard atm)
   *
   */
  upd->Pressure_Datum =  1.0132500000E6;

  iread = look_for_optional(ifp, "Pressure Datum", input, '=');

  if (iread == 1) {

    (void) read_string(ifp, input,'\n'); strip(input);

    upd->Pressure_Datum = parse_press_datum_input(input);

    SPF(echo_string,"%s = %g","Pressure Datum", upd->Pressure_Datum); ECHO(echo_string, echo_file);

  }

  upd->Acoustic_Frequency =  2*M_PIE*20000.;
  iread = look_for_optional(ifp, "Acoustic Frequency", input, '=');
  if (iread == 1) {
      if ( fscanf(ifp,"%le",&upd->Acoustic_Frequency) != 1)
 	{
 	  EH( -1, "error reading Acoustic Frequency");
  	}
    SPF(echo_string,"%s = %g","Acoustic Frequency", upd->Acoustic_Frequency); ECHO(echo_string, echo_file);
   }
 
  upd->Process_Temperature =  25.0;
  iread = look_for_optional(ifp, "Process Temperature", input, '=');
  if (iread == 1) {
      if ( fscanf(ifp,"%le",&upd->Process_Temperature) != 1)
 	{
 	  EH( -1, "error reading Process Temperature");
  	}
    SPF(echo_string,"%s = %g","Process Temperature", upd->Process_Temperature); ECHO(echo_string, echo_file);
   }
 
  upd->Light_Cosmu =  1.0;
  iread = look_for_optional(ifp, "Light Cosmu", input, '=');
  if (iread == 1) {
      if ( fscanf(ifp,"%le",&upd->Light_Cosmu) != 1)
 	{
 	  EH( -1, "error reading Light Incident Cos(Angle)");
  	}
    SPF(echo_string,"%s = %g","Light Cosmu", upd->Light_Cosmu); ECHO(echo_string, echo_file);
   }
 
  /*
   *  Anneal Mesh on Output - Optional
   *        Go through a mesh anneal step at the end of the run
   */
  Anneal_Mesh = FALSE;
  iread = look_for_optional(ifp,"Anneal Mesh on Output",input,'=');
  if(iread == 1) {
    (void) read_string(ifp,input,'\n');
    strip(input);
    if ( strcmp(input,"no") == 0 )
      {
	Anneal_Mesh = FALSE;
      }
    else if ( strcmp(input,"yes") == 0 )
      {
	Anneal_Mesh = TRUE;
        fprintf(stderr, "Annealing mesh on output\n");
      }
    else
      {
      EH( -1, "Bad specification for annealing mesh");
      }

    SPF(echo_string,eoformat, "Anneal Mesh On Output", input ); ECHO(echo_string, echo_file);
  }
}
/* rd_genl_specs -- read input file for general specifications */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * rd_timeint_specs -- read input file for time integration specifications
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Someday we need to comb through all these placeholder options
 *		to see if they are useful.
 *		
 *
 * Revised:			Fri Oct 29 09:00:27 MDT 1993 pasacki@sandia.gov
 */

void 
rd_timeint_specs(FILE *ifp,
		 char *input )
{
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  int mn, iread, i;
/*  double init_time; this global variable appears in rf_fem.h */
/*  double print_delt; this global variable appears in rf_fem.h */
  /* Whoa, how did this work? yo = "rd_timeint_specs"; */
  nEQM = 0;

  /* Read in Time Integration Specifications */
    
  iread = look_for_optional(ifp,"Time Integration Specifications",input,'=');

  ECHO("\n***Time Integration Specifications*** \n", echo_file);
    
  look_for(ifp,"Time integration",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input);
  if ( strcmp(input,"steady") == 0 )
    {
      TimeIntegration = STEADY;
    }
  else if ( strcmp(input,"transient") == 0 )
    {
      TimeIntegration = TRANSIENT;
    }
  else
    {
      EH( -1, "unknown time integration option");
    }

  SPF(echo_string,eoformat, "Time integration", input); ECHO(echo_string, echo_file);

  for (mn=0; mn<MAX_NUMBER_MATLS; mn++) {
    pd_glob[mn]->TimeIntegration =  TimeIntegration;
  }

  /* An initialization needed to avoid umrs for linear
     stab problems.  Needed once PRS started using tran
     structure as a global variable for poroelastic probs */
  tran->theta = 0.0;

  /* set default frequency to 0 */
  tran->fix_freq = 0;

  if(pd_glob[0]->TimeIntegration != STEADY) {
    double delta_t0;
    double delta_t_min;
    double delta_t_max;
    int max_time_steps;
    double time_max;

    look_for(ifp,"delta_t",input,'=');
    if (fscanf (ifp,"%le",&delta_t0) != 1)
      {
	EH( -1, "error reading delta_t");
      }
    tran->Delta_t0 = delta_t0;

    SPF(echo_string,"%s = %.4g","delta_t",tran->Delta_t0); ECHO(echo_string, echo_file);
    
    look_for(ifp,"Maximum number of time steps",input,'=');
    if (fscanf (ifp,"%d",&max_time_steps) != 1)
      {
	EH( -1, "error reading max time steps");
      }
    tran->MaxTimeSteps = max_time_steps;

    SPF(echo_string,"%s = %d", "Maximum number of time steps", tran->MaxTimeSteps); ECHO(echo_string, echo_file);
  
    look_for(ifp,"Maximum time",input,'=');
    if ( fscanf(ifp,"%le",&time_max) != 1)
      {
	EH( -1, "error reading maximum time");
      }
    tran->TimeMax = time_max;

    SPF(echo_string,"%s = %.4g","Maximum time",tran->TimeMax); ECHO(echo_string, echo_file);

    look_for(ifp,"Minimum time step",input,'=');
    if ( fscanf(ifp,"%le",&delta_t_min) != 1)
      {
	EH( -1, "error reading minimum time step");
      }
    tran->Delta_t_min = delta_t_min;

    SPF(echo_string,"%s = %.4g","Minimum time step",tran->Delta_t_min); ECHO(echo_string, echo_file);

    delta_t_max = 1.e12;
    iread = look_for_optional(ifp,"Maximum time step",input,'=');
    if (iread == 1) {   
      if ( fscanf(ifp,"%le",&delta_t_max) != 1)
	{
	  EH( -1, "error reading Maximum time step");
	}
      
      SPF(echo_string,"%s = %.4g","Maximum time step", delta_t_max); ECHO(echo_string, echo_file);

    }
	
        tran->Delta_t_max = delta_t_max;

    look_for(ifp,"Time step parameter",input,'=');
    if ( fscanf(ifp,"%le",&(tran->theta) ) != 1)
      {
      EH( -1, "error reading time step parameter");
    }

    SPF(echo_string,"%s = %.4g", "Time step parameter",tran->theta); ECHO(echo_string, echo_file);

    look_for(ifp,"Time step error",input,'=');
    if(fscanf(ifp, "%le", &eps) != 1 )
      {
        EH( -1, "error reading Time step error, expected at least one float");
      }
    tran->eps = eps;

    /* initialize norm indicators */
    for (i = 0; i < 9; i++) tran->use_var_norm[i] = 0;
    tran->use_var_norm[9] = 1; /* for backwards Compatibility */
    read_line(ifp, input, FALSE);

    if ( sscanf(input,"%d %d %d %d %d %d %d %d %d %d", 
		&tran->use_var_norm[0], &tran->use_var_norm[1], 
		&tran->use_var_norm[2], &tran->use_var_norm[3], 
		&tran->use_var_norm[4], &tran->use_var_norm[5],        
		&tran->use_var_norm[6], &tran->use_var_norm[7], 
		&tran->use_var_norm[8], &tran->use_var_norm[9]) < 8)   
      {
	fprintf(stdout, "Warning: Time step error prefers 1 flt 10 ints\n");
	fprintf(stdout, 
		"%s   d=%1d, v=%1d, T=%1d, y=%1d, P=%1d, S=%1d, V=%1d, sd=%1d, ls=%1d, ac=%1d\n", 
		"Best guess:", 
		tran->use_var_norm[0], tran->use_var_norm[1], 
		tran->use_var_norm[2], tran->use_var_norm[3],
		tran->use_var_norm[4], tran->use_var_norm[5],
		tran->use_var_norm[6], tran->use_var_norm[7],
		tran->use_var_norm[8], tran->use_var_norm[9]);
      }

    SPF(echo_string,"%s = %.4g %d %d %d %d %d %d %d %d %d %d", "Time step error", tran->eps, 
	tran->use_var_norm[0], tran->use_var_norm[1],
	tran->use_var_norm[2], tran->use_var_norm[3],
	tran->use_var_norm[4], tran->use_var_norm[5],
	tran->use_var_norm[5], tran->use_var_norm[7],
	tran->use_var_norm[8], tran->use_var_norm[9]); ECHO(echo_string, echo_file);

    look_for(ifp,"Printing Frequency",input,'=');
    print_freq = read_int(ifp, "Printing Frequency");
    tran->print_freq = print_freq;
    if (print_freq == 0) 
      {
      if ( fscanf(ifp,"%le",&print_delt) != 1)
	{
	  EH( -1, "error reading Printing delta time");
	}
      tran->print_delt = print_delt;

      SPF(echo_string,"%s = %d %.4g","Printing Frequency",print_freq,print_delt ); ECHO(echo_string, echo_file);

      print_delt2 = -print_delt;
      print_delt2_time = time_max;
      iread = look_for_optional(ifp,"Second frequency time",input,'=');
      if (iread == 1) 
	  {   
		  if ( fscanf(ifp,"%le %le",&print_delt2_time, &print_delt2) != 2)
		  {
			  EH( -1, "error reading second frequency time");
		  }
		  
		  SPF(echo_string,"%s = %.4g %.4g","Second frequency time",print_delt2_time, print_delt2 );ECHO(echo_string, echo_file);
	  }
	  tran->print_delt2 = print_delt2;
	  tran->print_delt2_time = print_delt2_time;
	  	  
      }
    else
      {
		SPF(echo_string,"%s = %d","Printing Frequency",tran->print_freq); ECHO(echo_string, echo_file);
      }

    /* Look for fix frequency */

    /* only look for fix frequency in parallel */
    if (Num_Proc > 1) {
      iread = look_for_optional(ifp,"Fix Frequency",input,'=');
      if (iread == 1) {
        tran->fix_freq = read_int(ifp, "Fix Frequency");
        if (tran->fix_freq < 0) {
          EH(-1, "Expected Fix Frequency > 0");
        }
        SPF(echo_string, "%s = %d", "Fix Frequency", tran->fix_freq); ECHO(echo_string, echo_file);
      }
    }

    tran->resolved_delta_t_min = 0.;
    iread = look_for_optional(ifp,"Minimum Resolved Time Step",input,'=');
    if (iread == 1) {
      if ( fscanf(ifp,"%le",&tran->resolved_delta_t_min) != 1)
	{
	  EH( -1, "error reading Minimum Resolved Time Step");
	}
      SPF(echo_string,"%s = %.4g","Minimum Resolved Time Step",tran->resolved_delta_t_min); ECHO(echo_string, echo_file);
    }
    
    tran->Courant_Limit = 0.;
    iread = look_for_optional(ifp,"Courant Number Limit",input,'=');
    if (iread == 1) {
      if ( fscanf(ifp,"%le",&tran->Courant_Limit) != 1)
	{
	  EH( -1, "error reading Courant Number Limit");
	}
      SPF(echo_string,"%s = %.4g","Courant Number Limit",tran->Courant_Limit); ECHO(echo_string, echo_file);
    }
    
    tran->Restart_Time_Integ_After_Renorm = TRUE;
    iread = look_for_optional(ifp,"Restart Time Integration After Renormalization",input,'=');
    if (iread == 1) {
      if ( fscanf(ifp,"%s",input) != 1)
        {
          EH( -1, "error reading Restart Time Integration After Renormalization.");
        }
      strip(input); stringup(input);
      if (strcmp(input,"ON")   == 0 ||
          strcmp(input,"YES")  == 0 ||
          strcmp(input,"TRUE") == 0 )
        {
          tran->Restart_Time_Integ_After_Renorm = TRUE;
        }
      else if(strcmp(input,"OFF")   == 0 ||
              strcmp(input,"NO")    == 0 ||
              strcmp(input,"FALSE") == 0 )
        {
          tran->Restart_Time_Integ_After_Renorm = FALSE;
        }
      /* Sanity check while we're here. */
      else
        {
          EH(-1,"Invalid setting for Restart Time Integration After Renormalization.");
        }
      SPF(echo_string,"%s = %d","Restart Time Integration After Renormalization",tran->Restart_Time_Integ_After_Renorm); ECHO(echo_string, echo_file);
    }

    tran->init_time = 0.0;
    iread = look_for_optional(ifp,"Initial Time",input,'=');
    if (iread == 1) {   
      if ( fscanf(ifp,"%le",&(tran->init_time) ) != 1)
	{
	  EH( -1, "error reading Initial Time");
	}
      SPF(echo_string,"%s = %.4g","Initial Time",tran->init_time); ECHO(echo_string, echo_file);
      if( (time_max - tran->init_time) <= 0.0 ) {
      	EH( -1, "Your maximum time is less than or equal to your initial time!"); // a condition which may result in NAN's in the esp_dot struct
      }
    }

    tran->const_dt_after_failure = 0.;
    iread = look_for_optional(ifp,"Steps of constant delta_t after failure",input,'=');
    if (iread == 1) {
      if ( fscanf(ifp,"%d",&tran->const_dt_after_failure) != 1)
	{
	  EH( -1, "error reading Steps of constant delta_t after failure");
	}
      SPF(echo_string,"%s = %d","Steps of constant delta_t after failure", tran->const_dt_after_failure); ECHO(echo_string, echo_file);
    }
    
    /*
     * If a time step fails to converge, then the time step is multiplied
     * by this factor (to make it smaller) for the next attempt.
     */
    tran->time_step_decelerator = 0.5;
    iread = look_for_optional(ifp,"Time step decelerator",input,'=');
    if (iread == 1) {
      if ( fscanf(ifp,"%le",&tran->time_step_decelerator) != 1)
	{
	  EH( -1, "error reading Time step decelerator");
	}

      SPF(echo_string,"%s = %.4g", "Time step decelerator",tran->time_step_decelerator);ECHO(echo_string, echo_file);
    }

#ifndef COUPLED_FILL
    tran->exp_subcycle = 10;
    iread = look_for_optional(ifp,"Fill Subcycle",input,'=');
    if (iread == 1) {   
      if ( fscanf(ifp,"%d",&(tran->exp_subcycle) ) != 1)
	{
	  EH( -1, "error reading Fill Subcycle");
	}
      SPF(echo_string,"%s = %d","Fill Subcycle", tran->exp_subcycle); ECHO(echo_string, echo_file);
    }
#endif /* not COUPLED_FILL */

    tran->Fill_Weight_Fcn = FILL_WEIGHT_TG;
    
    iread = look_for_optional(ifp,"Fill Weight Function",input,'=');
    
    if ( iread == 1)
      {
	if( fscanf(ifp,"%s",input) != 1)
	  {
	    EH( -1, "error reading Fill Weight Function string");
	  }
	
	if ( strcmp( input, "Galerkin") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_G;
	  }
	else if ( strcmp( input, "Taylor-Galerkin") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_TG;
	  }
	else if ( strcmp( input, "SUPG") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_SUPG;
	  }
	else if ( strcmp( input, "GLS") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_GLS;
	  }
	else if ( strcmp( input, "SC") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_SC;
	    WH(-1, "This option only works for extension velocity.\n");
	  }
        else if ( strcmp( input, "Explicit") == 0 )
	  {
	    tran->Fill_Weight_Fcn = FILL_WEIGHT_EXPLICIT;
	  }
	else
	  {
	    EH(-1, "Fill Weight Function not known.\n");
	  }
	SPF(echo_string,eoformat,"Fill Weight Function",input);ECHO(echo_string, echo_file);
      }

    /* Default to advection equation for fill/level set */

    tran->Fill_Equation = FILL_EQN_ADVECT; 
    strcpy(input,"Advection");
    
    iread = look_for_optional(ifp,"Fill Equation",input,'=');
    
    if ( iread == 1)
      {
	if( fscanf(ifp,"%s",input) != 1)
	  {
	    EH( -1, "error reading Fill Equation string");
	  }
	
	if ( strcmp( input, "Advection") == 0 )
	  {
	    tran->Fill_Equation = FILL_EQN_ADVECT;
	  }
	else if ( strcmp( input, "Extension") == 0 )
	  {
	    tran->Fill_Equation = FILL_EQN_EXT_V;
	    WH(-1, "This option needs an extension velocity eqn.\n");
	  }
	else if ( strcmp( input, "Eikonal") == 0 )
	  {
	    tran->Fill_Equation = FILL_EQN_EIKONAL;
	  }
	else
	  {
	    EH(-1, "Fill Equation not known.\n");
	  }

	SPF(echo_string,eoformat,"Fill Equation",input);ECHO(echo_string, echo_file);
      }
  } /*   if(pd_glob[0]->TimeIntegration != STEADY) */

}    
/* rd_timeint_specs -- read input file for time integration specifications */

/*****************************************************************************/
/*
 * rd_levelset_specs -- read input file for level set specifications
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 */

void 
rd_levelset_specs(FILE *ifp,
		  char *input )
{
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  int iread, i;
  
  /* Level Set interface tracking parameters */

  iread = look_for_optional(ifp,"Level Set Interface Tracking",input,'=');
  ls  = NULL;
  lsi = NULL;
  if (iread == 1) 
    {   

      if (fscanf(ifp, "%s", input ) != 1 )
         {
           EH(-1, "Error reading Level Set Interface Tracking .");
         }

      strip(input); stringup(input);

      if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
        {
          ls =  alloc_struct_1(struct Level_Set_Data, 1);
          ls->var = FILL;
          ls->embedded_bc = NULL;
          ls->init_surf_list = NULL;
          ls->last_surf_list = NULL;
          ls->sm_object_name = NULL;
          ls->sm_object_type = NULL;
          lsi = alloc_struct_1(struct Level_Set_Interface, 1);
          ls->Use_Level_Set = Use_Level_Set = TRUE ;
          zero_lsi();
	  zero_lsi_derivs();
          ls->Elem_Sign = 0;
          ls->on_sharp_surf = 0;
          ls->Extension_Velocity = FALSE;
	      ls->CalcSurfDependencies = FALSE;
	      ls->Integration_Depth = 0;
		  ls->AdaptIntegration = FALSE;
	  	  ls->Sat_Hyst_Renorm_Lockout = 0;
        }

      ECHO("\n***Level Set Interface Tracking***\n", echo_file);
      SPF(echo_string,eoformat,"Level Set Interface Tracking",input); ECHO(echo_string,echo_file);

    }

  if (ls != NULL)
    {    
      /* for steady-state level set problems */
      if(pd_glob[0]->TimeIntegration == STEADY)
        {
          WH(-1, "Steady state level set problem.  Using Eikonal equation!\n");
          tran->Fill_Equation = FILL_EQN_EIKONAL;
        }
        

      /* Need to read in initialization data for fields that may be initially index wrt LS 
      *  Also need to run this first.   */

      ls->Num_Var_Init = 0;
      while ((iread = look_forward_optional(ifp,"Level Set Initialize",input,'=')) == 1)
        {         
	  int index = Num_Var_Init + ls->Num_Var_Init;
          if ( fscanf(ifp, "%80s", input ) != 1 )
            {
              EH(-1, "Error reading LS initialization field data \n");
            }

          /* OK.  Let us put the Variables initialized index by LS on the end of the original Var_init array. 
	  *       This is why we add Num_Var_Init to the ls->Num_Var_Init to define index 
	  */

          Var_init[index].var = 
            variable_string_to_int(input, "Initialize Keyword Error");

          if ( fscanf( ifp, "%d %lf %lf", 
                       &(Var_init[index].ktype), 
                       &(Var_init[index].init_val_minus),
                       &(Var_init[index].init_val_plus) ) != 3 )
            {
              EH(-1,"Need one int and two floats on LS initialization card \n");
            }

	  SPF(echo_string,"%s = %s %d %.4g %.4g","Level Set Initialize",input,
	      Var_init[index].ktype,Var_init[index].init_val_minus,Var_init[index].init_val_plus ); ECHO(echo_string,echo_file);

          ls->Num_Var_Init++;
        }

      Num_Var_LS_Init = ls->Num_Var_Init;
 
      ls->Length_Scale  = -1.0;
      iread = look_for_optional(ifp,"Level Set Length Scale",input,'=');
      if (iread == 1) 
        {   
          if ( fscanf(ifp,"%lf",&(ls->Length_Scale) ) != 1)
            {
             EH( -1, "error reading Level Set Length Scale");
            }
          if( fabs( ls->Length_Scale ) == 0.0 ) WH(-1,"Possible Syntax error.  Level set length scale is zero.");

	  SPF(echo_string,"%s = %.4g",input, ls->Length_Scale);
        }
      else
	{
	  SPF(echo_string," (%s = %f) %s","Level Set Length Scale", ls->Length_Scale, default_string); 
	}

      ECHO(echo_string,echo_file);

      ls->Init_Method  = -1;
      iread = look_for_optional(ifp,"Level Set Initialization Method",input,'=');
      if (iread == 1) 
        {   
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Level Set Initialization string");
            }
 
          if ( strcmp( input,"Projection") == 0 )
            {
              ls->Init_Method = PROJECT;  
	      ECHO("Level Set Initialization Method = Projection", echo_file);
            }
          else if  ( strcmp( input,"Exodus") == 0 )
            {
              ls->Init_Method =  EXO_READ;  
	      ECHO("Level Set Initialization Method = Exodus", echo_file);
            }
          else if  ( strcmp( input,"Nodeset") == 0 )
            {
              struct LS_Surf *surf;
              struct LS_Surf_NS_Data *s;
              ls->Init_Method = SURFACES;

	      SPF(echo_string,eoformat, "Level Set Initialization Method", input);

              ls->init_surf_list = create_surf_list();
              surf = create_surf( LS_SURF_NS );
              append_surf( ls->init_surf_list, surf );
              s = (struct LS_Surf_NS_Data *) surf->data;
              
              if ( fscanf(ifp, "%s %d %s %d", input, &(s->ns_id), input, &(s->PosEB_id)) != 4) 
                {
                  EH(-1,"Nodeset Initialization requires syntax:  Nodeset NS node_side_id EB e_block_id \n");
                }

	      SPF(endofstring(echo_string)," NS %d EB %d",s->ns_id, s->PosEB_id ); ECHO(echo_string,echo_file);

            }
          else if ( ( strcmp( input, "Surfaces") == 0 ) )
            {
              int num_surf;
              
              ls->Init_Method = SURFACES;

              if ( fscanf( ifp, "%d" , &num_surf ) != 1 )
                {
                  EH(-1,"Surfaces initialization method needs number of surfaces specified.");
                }

	      SPF(echo_string,"Level Set Initialization Method = Surfaces %d",num_surf); ECHO(echo_string,echo_file);

              ls->init_surf_list = create_surf_list();
              read_surface_objects ( ifp, input, ls->init_surf_list, num_surf); 
            }
          else if(!strcmp(input, "SM_object"))
	      {
	      EH(-1, "CGM not supported, SM_object");
            }
          else
            WH(-1, "Level Set Initialization method undefined");
        }
	
      ls->Periodic_Planes = FALSE;
      for ( i = 0; i < 6; i++ ) ls->Periodic_Plane_Loc[i] = 0.;
      iread = look_for_optional(ifp,"Level Set Periodic Planes",input,'=');
      if (iread == 1) 
        {   
	  if ( fscanf(ifp,"%lf %lf %lf %lf %lf %lf",
	              &(ls->Periodic_Plane_Loc[0]), &(ls->Periodic_Plane_Loc[1]),
	              &(ls->Periodic_Plane_Loc[2]), &(ls->Periodic_Plane_Loc[3]),
	              &(ls->Periodic_Plane_Loc[4]), &(ls->Periodic_Plane_Loc[5])) != 6 )
            {
              EH( -1, "error reading Level Set Periodic Planes");
            }
	  ls->Periodic_Planes = TRUE;

	  SPF(echo_string,"Level Set Periodic Planes = %.4g %.4g %.4g %.4g %.4g %.4g",
	      ls->Periodic_Plane_Loc[0], ls->Periodic_Plane_Loc[1], 
	      ls->Periodic_Plane_Loc[2], ls->Periodic_Plane_Loc[2], 
	      ls->Periodic_Plane_Loc[3], ls->Periodic_Plane_Loc[3]); ECHO(echo_string,echo_file); 
        }

      ls->Control_Width = 1.0;
  
      iread = look_for_optional(ifp,"Level Set Control Width",input,'=');
   
      if (iread == 1) 
        {   
          if ( fscanf(ifp,"%lf",&(ls->Control_Width) ) != 1)
            {
              EH( -1, "error reading Level Set Length Scale");
            }

          if ( ls->Control_Width > 10.0 ) WH(-1,"That's an awfully large Level Set Control Width.\n");
	  
	  SPF(echo_string,"%s = %.4g",input, ls->Control_Width); ECHO(echo_string,echo_file);

        }
      else
	SPF(echo_string," (%s = %f) %s", "Level Set Control Width",ls->Control_Width, default_string); 

      ECHO(echo_string,echo_file);
	  
      tran->use_var_norm[8] = 0;
	  
      iread = look_for_optional(ifp,"Level Set Timestep Control",input,'=');
	  
      if (iread == 1) 
	  {   
	    if ( fscanf(ifp,"%s",input ) != 1)
	      {
		EH( -1, "Expecting string for Level Set Timestep Control");
	      }
	    stringup(input);
	    if ( !strcmp(input,"YES") || !strcmp(input,"ON") ) tran->use_var_norm[8] = 1;

	    SPF(echo_string,eoformat,"Level Set Timestep Control",input); ECHO(echo_string,echo_file);
	  }
    
      ls->Renorm_Tolerance = 0.5;
  
      iread = look_for_optional(ifp,"Level Set Renormalization Tolerance",input,'=');
   
      if (iread == 1) 
        {   
          if ( fscanf(ifp,"%lf",&(ls->Renorm_Tolerance) ) != 1)
            {
              EH( -1, "error reading Level Set Renormalization tolerance");
            }
          if ( ls->Renorm_Tolerance > 1.0 ) WH(-1,"That's an awfully large Level Set Renormalization Tolerance.\n");

	  SPF(echo_string,"%s = %.4g","Level Set Renormalization Tolerance", ls->Renorm_Tolerance);
        }
      else
	SPF(echo_string," (%s = %f) ", "Level Set Renormalization Tolerance",ls->Renorm_Tolerance ); 

      ECHO(echo_string,echo_file);

      ls->Renorm_Method = FALSE;

      iread = look_for_optional(ifp,"Level Set Renormalization Method",input,'=');
  
      if ( iread == 1 )
        {
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Level Set Renormalization string");
            }

	  SPF(echo_string,"%s = ","Level Set Renormalization Method");

          if ( strcmp( input,"Correction") == 0 )
            {
              ls->Renorm_Method = CORRECT;

	      strcat(echo_string, "Correction");
            }
          else if  ( strcmp( input,"Huygens") == 0 )
            {
              ls->Renorm_Method = HUYGENS;
	      strcat(echo_string, "Huygens");
            }
          else if  ( strcmp( input,"Huygens_Constrained") == 0 )
            {

              ls->Renorm_Method = HUYGENS_C;

	      strcat(echo_string, "Huygens_Constrained");

              if( fscanf( ifp,"%lf", &(ls->Mass_Value)) == 1) 
                {
		  char *s = endofstring(echo_string);

                  ls->Mass_Sign = ls->Mass_Value <= 0 ? I_NEG_FILL : I_POS_FILL;
                  ls->Mass_Value = fabs( ls->Mass_Value);

		  SPF(s," %.4g",ls->Mass_Sign*ls->Mass_Value);
                }
              else
                {
                  ls->Mass_Value = 0.0;
                  ls->Mass_Sign  = I_NEG_FILL;
                }                 
            }
          else if  ( ( strcmp( input,"None") == 0 ) ||  (strcmp( input,"No") == 0) )
            {
              ls->Renorm_Method = FALSE;
	      strcat(echo_string, "None");
            }
          else
            {
              EH(-1,"Level Set Renorm Method not known.\n");
            }

	  ECHO(echo_string,echo_file);
        }
        
      ls->Ignore_F_deps = FALSE;
      iread = look_for_optional(ifp,"Ignore Level Set Dependencies",input,'=');
      if ( iread == 1 )
        {
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Ignore Level Set Dependencies");
            }
          strip(input); stringup(input);
          
          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
#ifdef COUPLED_FILL
              if ( tran->Fill_Equation == FILL_EQN_ADVECT &&
	           tran->Fill_Weight_Fcn != FILL_WEIGHT_EXPLICIT )
                {
                  EH(-1,"Maybe I should let you ignore LS dependencies even with implicit LS, but I won't!\n");
                }
#endif
              ls->Ignore_F_deps = TRUE;
            }
	  SPF(echo_string,eoformat, "Ignore Level Set Dependencies", input); ECHO(echo_string,echo_file);
        }

      ls->Interface_Output = FALSE;
      ls->output_file = NULL;

      iread = look_for_optional(ifp,"Level Set Output Filename",input,'=');

      if (iread == 1 )
        {
          if ( fscanf( ifp, "%s", input ) != 1 )
            {
              EH(-1,"Need string parameter for level set output filename .\n");
            }

          ls->output_file = smalloc( sizeof(char)* ( strlen(input) + 1 ) );

          strcpy (ls->output_file, input );

          ls->Interface_Output = TRUE;

	  SPF(echo_string,eoformat,"Level Set Output Filename", input); ECHO(echo_string,echo_file);
        }

      ls->Renorm_Freq = ls->Renorm_Countdown = -1;  /* Default is never to renormalize */

      iread = look_for_optional(ifp,"Level Set Renormalization Frequency",input,'=');
  
      if ( iread == 1 )
        {
          if( fscanf( ifp,"%d", &(ls->Renorm_Freq) ) != 1 )
            {
              EH(-1,"Error reading renormalization frequency.\n");
            }
          ls->Renorm_Countdown = ls->Renorm_Freq;

	  SPF(echo_string,"%s = %d","Level Set Renormalization Frequency", ls->Renorm_Freq); ECHO(echo_string,echo_file);
        }

      ls->Force_Initial_Renorm = FALSE;
      iread = look_for_optional(ifp,"Force Initial Level Set Renormalization",input,'=');
      if ( iread == 1 )
        {
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Force Initial Level Set Renormalization flag");
            }
          strip(input); stringup(input);
          
          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
              ls->Force_Initial_Renorm = TRUE;
            }

	  SPF(echo_string,eoformat,"Force Initial Level Set Renormalization", input); ECHO(echo_string,echo_file);
        }
	
      ls->Initial_LS_Displacement = 0.;
      iread = look_for_optional(ifp,"Initial Level Set Displacement",input,'=');
      if (iread == 1) 
        {   
          if ( fscanf(ifp,"%lf",&(ls->Initial_LS_Displacement) ) != 1)
            {
              EH( -1, "error reading Initial Level Set Displacement");
            }

	  SPF(echo_string,"%s = %.4g","Initial Level Set Displacement", ls->Initial_LS_Displacement); ECHO(echo_string,echo_file);
        }
        
      ls->Isosurface_Subsurf_Type = LS_SURF_POINT;
      iread = look_for_optional(ifp,"Level Set Reconstruction Method",input,'=');
      if ( iread == 1 )
        {
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Level Set Reconstruction Method");
            }
          strip(input); stringup(input);

          if ( strcmp( input,"POINTS") == 0 )
            {
              ls->Isosurface_Subsurf_Type = LS_SURF_POINT;
            }
          else if ( strcmp( input,"FACETS") == 0 )
            {
              ls->Isosurface_Subsurf_Type = LS_SURF_FACET;
            }
          else
            {
              EH(-1,"Level Set Reconstruction Method not known.\n");
            }
	  SPF(echo_string,eoformat, "Level Set Reconstruction Method", input); ECHO(echo_string,echo_file);
        }

      ls->Contact_Inflection = FALSE;
      iread = look_for_optional(ifp,"Level Set Contact Extension",input,'=');
      if ( iread == 1 )
        {
          if ( fscanf(ifp,"%s",input) != 1)
            {
              EH( -1, "error reading Level Set Contact Extension flag");
            }
          strip(input); stringup(input);
          
          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
              ls->Contact_Inflection = TRUE;
            }

	  SPF(echo_string,eoformat, "Level Set Contact Extension", input); ECHO(echo_string,echo_file);
        }

      /* set default evolution scheme based on compiler flag */
#ifdef COUPLED_FILL
      ls->Evolution = LS_EVOLVE_ADVECT_COUPLED;
#else
      ls->Evolution = LS_EVOLVE_ADVECT_EXPLICIT;
#endif

      iread = look_for_optional(ifp,"Level Set Slave Surface",input,'=');
      if (iread == 1)
        {
          if (fscanf(ifp, "%s", input ) != 1 )
             {
               EH(-1, "Error reading Level Set Slave Surface.");
             }
          strip(input); stringup(input);

          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
              ls->Evolution = LS_EVOLVE_SLAVE;
              ls->Renorm_Method = FALSE;
              if (ls->Init_Method != SURFACES)
                EH(-1,"Must define Level Set Initialization method to be SURFACES or NODESET for slave surfaces");
            }

	  SPF(echo_string,eoformat,"Level Set Slave Surface", input); ECHO(echo_string,echo_file);
        }

      ls->Contact_Tolerance = 0.5;

      /* Check if this is a fluid/solid interaction problem. */
      ls->Fluid_Solid = FALSE;
      iread = look_for_optional(ifp,"Level Set Fluid Solid",input,'=');
      if (iread == 1)
        {
          if (fscanf(ifp, "%s", input ) != 1 )
             {
               EH(-1, "Error reading Level Set Fluid Solid.");
             }
          strip(input); stringup(input);

          if (strcmp(input,"ON")   == 0 ||
              strcmp(input,"YES")  == 0 ||
              strcmp(input,"TRUE") == 0 )
            {
              ls->Fluid_Solid = TRUE;
            }
          /* Sanity check while we're here. */
          else if(strcmp(input,"OFF")   != 0 &&
                  strcmp(input,"NO")    != 0 &&
                  strcmp(input,"FALSE") != 0 )
            {
              EH(-1,"Invalid setting for Level Set Fluid Solid.");
            }

	  SPF(echo_string,eoformat,"Level Set Fluid Solid", input); ECHO(echo_string,echo_file);
        }

      /*
       * If this is a fluid/solid interaction problem, see which side of
       * the zero level set is the solid phase.  Note that Level Set Fluid
       * Solid must be specified for this to be used.
       */
      ls->Fluid_Sign = 0;
      ls->Solid_Sign = 0;
      if ( ls->Fluid_Solid )
        {
          look_for(ifp,"Level Set Solid Sign",input,'=');         
          ls->Solid_Sign = read_int(ifp, "Level Set Solid Sign");
          if ( ls->Solid_Sign > 0 )
              ls->Solid_Sign = 1;
          else if ( ls->Solid_Sign < 0 )
              ls->Solid_Sign = -1;
          else
            EH(-1,"Level Set Solid Sign requires a non-zero integer parameter.");

          ls->Fluid_Sign = -ls->Solid_Sign;

	  SPF(echo_string,"%s = %d","Level Set Solid Sign", ls->Solid_Sign); ECHO(echo_string,echo_file);
        }

      iread = look_for_optional(ifp,"Level Set Semi_Lagrange",input,'=');
      if (iread == 1)
        {
          if ( fscanf(ifp,"%s", input ) != 1)
            {
              EH(-1, "Need string parameter for Semi_Lagrange card.\n");
            }
          strip(input); stringup(input);

          if( strcmp(input, "YES" ) == 0 ||
              strcmp(input, "ON"  ) == 0 ||
              strcmp(input, "TRUE" ) == 0)
            {
#ifdef COUPLED_FILL
	      EH(-1, "Level set Semi_Lagrange not supported for COUPLED_FILL.");
#else
              ls->Evolution = LS_EVOLVE_SEMILAGRANGIAN;
              tran->exp_subcycle = 1.0;
#endif
            }
	  SPF(echo_string,eoformat,"Level Set Semi_Lagrange", input); ECHO(echo_string,echo_file);
        }

      ls->Search_Option = SEGMENT_SEARCH;
      ls->Grid_Search_Depth = 0;

      iread = look_for_optional(ifp,"Level Set Search Option",input,'=');
      if( iread == 1 ){
        if( fscanf(ifp, "%s", input) != 1)
          {
            EH(-1," Need string parameter for Search Option card.");
          }

      strip(input); stringup(input);

      SPF(echo_string,"%s = %s","Level Set Search Option", input);

      if( strcmp(input, "GRID_SEARCH") == 0 )
        {
	  int err;
	  EH(-1,"The Level Set Search Option : GRID_SEARCH is not functioning at this time.\n");

          ls->Search_Option = GRID_SEARCH;
          err = fscanf(ifp,"%d", &(ls->Grid_Search_Depth));
	  if (err != 1) {
	    EH(-1, "Expected to read one int for GRID_SEARCH");
	  }

	  SPF(endofstring(echo_string)," %d",ls->Grid_Search_Depth);
        }
      else if ( strcmp( input, "SEGMENT") == 0 )
        {
          ls->Search_Option = SEGMENT_SEARCH;
        }
      else 
        {
          EH(-1," Allowable choices for Search Option are: SEGMENT or GRID_SEARCH");
        }

      ECHO(echo_string,echo_file);
      }

      ls->Integration_Depth = 0;

      iread = look_for_optional( ifp, "Level Set Subgrid Integration Depth", input, '=' );

      if ( iread == 1 ) {
        if ( fscanf( ifp,"%d", &(ls->Integration_Depth) ) != 1 )
          {
            EH(-1, "Subgrid Integration Depth for Level Set needs a single positive integer parameter.");
          }
	SPF(echo_string,"%s = %d", "Level Set Subgrid Integration Depth", ls->Integration_Depth); ECHO(echo_string,echo_file);
      }

      ls->SubElemIntegration = FALSE;
      iread = look_for_optional(ifp,"Level Set Subelement Integration",input,'=');
      if (iread == 1)
        {
          if (fscanf(ifp, "%s", input ) != 1 )
             {
               EH(-1, "Error reading Level Set Subelement Integration Flag.");
             }
          strip(input); stringup(input);

          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
              ls->SubElemIntegration = TRUE;
              if (ls->Integration_Depth > 0)
                EH(-1,"Combination of Subgrid and Subelement integration is not supported.");
              if (ls->Length_Scale > 0.)
                EH(-1,"I don't think it makes sense to have subelement integration and non-zero length scale.");
            }

	  SPF(echo_string,eoformat,"Level Set Subelement Integration", input); ECHO(echo_string,echo_file);
        }

      ls->AdaptIntegration = FALSE;
      iread = look_for_optional(ifp,"Level Set Adaptive Integration",input,'=');
      if (iread == 1)
        {
          if (fscanf(ifp, "%s", input ) != 1 )
             {
               EH(-1, "Error reading Level Set Adaptive Integration Flag.");
             }
          strip(input); stringup(input);

          if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
            {
              ls->AdaptIntegration = TRUE;
              if (ls->Integration_Depth > 0)
                EH(-1,"Combination of Subgrid and Adaptive integration is not supported.");
              if (ls->Length_Scale > 0.)
                EH(-1,"Currently, Adaptive integration applies only to sharp interfaces (zero LS length scale).");
              if (ls->SubElemIntegration)
                EH(-1,"Combination of Subelement and Adaptive integration is not supported.");
            }
	  SPF(echo_string,eoformat,"Level Set Adaptive Integration", input); ECHO(echo_string,echo_file);
        }

      ls->Adaptive_Order = 3;
  
      iread = look_for_optional( ifp, "Level Set Adaptive Order", input, '=' );
  
      if ( iread == 1 ) {
        if ( fscanf( ifp,"%d", &(ls->Adaptive_Order) ) != 1 )
          {
            EH(-1, "Adaptive Integration Order for Level Set needs a single positive integer parameter.");
          }
	SPF(echo_string,"%s = %d", "Level Set Adaptive Order", ls->Adaptive_Order); ECHO(echo_string,echo_file);
      }

        
      ls->CrossMeshQuadPoints = 0;
      iread = look_for_optional(ifp,"Overlap Quadrature Points",input,'=');
      if ( iread == 1 ) {
        if ( fscanf( ifp,"%d", &(ls->CrossMeshQuadPoints) ) != 1 )
          {
            EH(-1, "Overlap Quadrature Points needs a single positive integer parameter.");
          }
	SPF(echo_string,"%s = %d","Overlap Quadrature Points" ,ls->CrossMeshQuadPoints ); ECHO(echo_string,echo_file);
      }
	  
	  /*
	   *  This flag is set so that the PSPP stabilization parameter will be weighted by 1.0 - delta * alpha at the interface
	   *  This will filter the PSPP stabilization at the interface region and improve mass conservatino
	   *  At the same time, convergence and matrix condition may be adversely affected
	   */
	  
	ls->PSPP_filter = FALSE;
	
	iread = look_for_optional(ifp,"Level Set PSPP filtering",input,'=');
	if (iread == 1)
	{
		if (fscanf(ifp, "%s", input ) != 1 )
		{
			EH(-1, "Error reading Level Set PSPP filtering flag.");
		}
		strip(input); stringup(input);
		
		if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
		{
			ls->PSPP_filter = TRUE;
		}
	}
	
	
    }  /* if ( ls != NULL ) */


  /*
   * Phase function input data 
   */
  
  iread = look_for_optional( ifp, "Number of phase functions", input, '=');
  
  pfd = NULL;
  if ( iread == 1 )
    {
      if( fscanf(ifp, "%d", &i) != 1 )
        {
          EH(-1,"Error scanning number of phase functions\n");
        }
      
      if( i > 0) 
        {
          pfd = alloc_struct_1(struct Phase_Function_Data, 1);
          pfd->num_phase_funcs = i;
	  pfd->Use_Phase_Field = Use_Phase_Field = i ;
	  
          pfd->ls = ( struct Level_Set_Data **) alloc_ptr_1(pfd->num_phase_funcs);
          for( i=0; i<pfd->num_phase_funcs;i++)
		  {
              pfd->ls[i] =  alloc_struct_1(struct Level_Set_Data, 1);
              pfd->ls[i]->var = PHASE1 + i;
              pfd->ls[i]->embedded_bc = NULL;
              pfd->ls[i]->Renorm_Method = FALSE;
              pfd->ls[i]->init_surf_list = NULL;
              pfd->ls[i]->last_surf_list = NULL;
              pfd->ls[i]->sm_object_name = NULL;
              pfd->ls[i]->sm_object_type = NULL;
	
	      if (ls) {
		pfd->ls[i]->Control_Width = ls->Control_Width;
		pfd->ls[i]->Renorm_Freq = ls->Renorm_Freq;
		pfd->ls[i]->Renorm_Countdown = ls->Renorm_Countdown;
		pfd->ls[i]->Renorm_Tolerance = ls->Renorm_Tolerance;
		pfd->ls[i]->Force_Initial_Renorm = ls->Force_Initial_Renorm;
	      } else {
		pfd->ls[i]->Control_Width = 1.0;
		pfd->ls[i]->Renorm_Freq = -1;
		pfd->ls[i]->Renorm_Countdown = -1;
		pfd->ls[i]->Renorm_Tolerance = 0.5;
		pfd->ls[i]->Force_Initial_Renorm = FALSE;
	      }
              pfd->ls[i]->Init_Method = -1;
              pfd->ls[i]->Elem_Sign = 0;
              pfd->ls[i]->on_sharp_surf = 0;
              pfd->ls[i]->Extension_Velocity = FALSE;

	      pfd->ls[i]->CalcSurfDependencies = FALSE;
	      pfd->ls[i]->Ghost_Integ_Active = FALSE;
	      pfd->ls[i]->Periodic_Planes = FALSE;
	      pfd->ls[i]->Num_Var_Init = 0;
	      if (ls) {
		pfd->ls[i]->Integration_Depth = ls->Integration_Depth;
		pfd->ls[i]->SubElemIntegration = ls->SubElemIntegration;
		pfd->ls[i]->AdaptIntegration = ls->AdaptIntegration;   	   
		pfd->ls[i]->Contact_Inflection = ls->Contact_Inflection;
		pfd->ls[i]->Ignore_F_deps = ls->Ignore_F_deps;
		pfd->ls[i]->Interface_Output = ls->Interface_Output;
		pfd->ls[i]->output_file = ls->output_file;
	      } else {
		pfd->ls[i]->Integration_Depth = 0;
		pfd->ls[i]->SubElemIntegration = FALSE;
		pfd->ls[i]->AdaptIntegration = FALSE;   	   
		pfd->ls[i]->Contact_Inflection = FALSE;
		pfd->ls[i]->Ignore_F_deps = FALSE;
		pfd->ls[i]->Interface_Output = FALSE;
		pfd->ls[i]->output_file = NULL;
	      }
	      pfd->ls[i]->Ghost_Integ = FALSE;

		if ( ls != NULL ) pfd->ls[i]->CrossMeshQuadPoints = ls->CrossMeshQuadPoints;
			  else pfd->ls[i]->CrossMeshQuadPoints = 0.;

#ifdef PHASE_COUPLED_FILL
	      pfd->ls[i]->Evolution = LS_EVOLVE_ADVECT_COUPLED;
#else
	      pfd->ls[i]->Evolution = LS_EVOLVE_ADVECT_EXPLICIT;
#endif
			  

            }
            
                    lsi = alloc_struct_1(struct Level_Set_Interface, 1);

		  ECHO("\n***Phase Function Section***\n", echo_file);
		  SPF(echo_string,"%s = %d","Number of phase functions", pfd->num_phase_funcs); ECHO(echo_string,echo_file);
        }
	  else
	  {
		  ECHO("\n***Phase Function Section***\n", echo_file);
		  SPF(echo_string,"%s = %d","Number of phase functions", i); ECHO(echo_string,echo_file);
	  }		  
	}

if( pfd != NULL ) {
    
	
	iread = look_for_optional(ifp,"Phase Function Slave Surface",input,'=');
	if (iread == 1)
    {
		if (fscanf(ifp, "%s", input ) != 1 )
		{
			EH(-1, "Error reading Phase Function Slave Surface.");
		}
		strip(input); stringup(input);
		
		if ((strcmp(input,"ON") == 0) || (strcmp(input,"YES") == 0 ))
		{
#ifdef PHASE_COUPLED_FILL
			EH(-1, "Phase function slave surfaces not currently supported for PHASE_COUPLED_FILL.");
#endif
			for( i=0; i<pfd->num_phase_funcs;i++)
			{
				pfd->ls[i]->Evolution = LS_EVOLVE_SLAVE;
				pfd->ls[i]->Renorm_Method = FALSE;
			}
		} 
		SPF(echo_string,eoformat,"Phase Function Slave Surface", input); ECHO(echo_string,echo_file);
    }
	
	if(pfd != NULL)
    {
		if(pfd->ls[0]->Evolution != LS_EVOLVE_SLAVE)
		{
#ifndef PHASE_COUPLED_FILL
			EH(-1,"Must have PHASE_COUPLED_FILL compile flag for phase-field functions that are not SLAVED");
#endif
		}
    }
	
    
	iread = look_for_optional(ifp,"Phase Function Initialization Method",input,'=');
	if (iread == 1) 
    {   
		for( i=0; i< pfd->num_phase_funcs; i++)
		{
			
			if ( fscanf(ifp,"%s",input) != 1)
			{
				EH( -1, "error reading Phase Function Initialization string");
			}
			
			
			if ( strcmp( input,"Projection") == 0 )
			{
				pfd->ls[i]->Init_Method = PROJECT;
				ECHO("Phase Function Initialization Method = Projection", echo_file);
			}
			else if  ( strcmp( input,"Exodus") == 0 )
			{
				pfd->ls[i]->Init_Method =  EXO_READ;
				ECHO("Phase Function Initialization Method = Exodus", echo_file);
			}
			else if  ( strcmp( input,"Nodeset") == 0 )
			{
				struct LS_Surf *surf;
				struct LS_Surf_NS_Data *s;
				pfd->ls[i]->Init_Method = SURFACES;
				
				SPF(echo_string, eoformat, "Phase Function Initialization Method", input);
				
				pfd->ls[i]->init_surf_list = create_surf_list();
				surf = create_surf( LS_SURF_NS );
				append_surf( pfd->ls[i]->init_surf_list, surf );
				s = (struct LS_Surf_NS_Data *) surf->data;
				
				if ( fscanf(ifp, "%s %d %s %d", input, &(s->ns_id), input, &(s->PosEB_id)) != 4) 
				{
					EH(-1,"Nodeset Initialization requires syntax:  Nodeset NS node_side_id EB e_block_id \n");
				}
				
				SPF(endofstring(echo_string),"NS %d EB %d", s->ns_id, s->PosEB_id); ECHO(echo_string,echo_file);
				
			}
			else if ( ( strcmp( input, "Surfaces") == 0 ) )
			{
				int num_surf;
				
				pfd->ls[i]->Init_Method = SURFACES;
				
				if ( fscanf( ifp, "%d" , &num_surf ) != 1 )
				{
					EH(-1,"Surfaces initialization method for phase function needs number of surfaces specified.");
				}
				
				SPF(echo_string,"%s = %s %d", "Phase Function Initialization Method", input, num_surf);ECHO(echo_string,echo_file);
				
				pfd->ls[i]->init_surf_list = create_surf_list();
				read_surface_objects ( ifp, input, pfd->ls[i]->init_surf_list, num_surf); 
			}
			
		}
    }
	
	iread = look_for_optional(ifp,"Phase Function Length Scale",input,'=');
	if (iread == 1) 
    {   
		double length_scale;
		
		if ( fscanf(ifp,"%lf",&length_scale) != 1)
		{
			EH( -1, "error reading Phase Function Initialization string");
		}
		
		for(i=0; i<pfd->num_phase_funcs; i++)
			pfd->ls[i]->Length_Scale = length_scale;
		
		SPF(echo_string,"%s = %.4g", "Phase Function Length Scale", length_scale); ECHO(echo_string,echo_file);
    }
	
	iread = look_for_optional(ifp,"Phase Function Renormalization Tolerance",input,'=');
	
	if (iread == 1) 
    {   
		double tolerance;
		if ( fscanf(ifp,"%lf",&tolerance) != 1)
		{
			EH( -1, "error reading Phase Function Renormalization tolerance");
		}
		if ( tolerance> 1.0 ) WH(-1,"That's an awfully large Level Set Renormalization Tolerance.\n");
		
		for(i=0; i<pfd->num_phase_funcs; i++) 
			pfd->ls[i]->Renorm_Tolerance = tolerance;
		
		SPF(echo_string,"%s = %.4g", "Phase Function Renormalization Tolerance", tolerance); ECHO(echo_string,echo_file);
    }
	
	iread = look_for_optional(ifp,"Phase Function Renormalization Method",input,'=');
	if (iread == 1) 
	{   
		int method=0, Mass_Sign = I_NEG_FILL;
		double Mass_Value = 0;
		
		if ( fscanf(ifp,"%s",input) != 1)
		{
			EH( -1, "error reading Phase Function Initialization string");
		}          
		
		if  ( strcmp( input,"Huygens") == 0 )
		{
			method = HUYGENS;
		}
		else if ( strcmp( input,"Huygen_Constrained") == 0 )
		{
			method = HUYGENS_C;
		}
		else if ( ( strcmp( input,"No") == 0) || ( strcmp( input,"None") == 0) )
		{
			method = FALSE;
			/* NOTE:  PRS tested these in the new era that all overset grid activities
			* will use RPhase0 and these don't work. The main reason is that there are
			* still direct references to "FILL" and "LS" and esp->F  in these routines
			* that need to be generalized.  3/24/05
			*/
		}
		else
		{
			EH(-1,"Unknown Phase Function Renorm method. n");
		}
		
		for(i=0; i<pfd->num_phase_funcs; i++)
		{
			pfd->ls[i]->Renorm_Method = method;
			pfd->ls[i]->Mass_Value = Mass_Value;
			pfd->ls[i]->Mass_Sign = Mass_Sign;
		}
		SPF(echo_string,eoformat,"Phase Function Renormalization Method", input); ECHO(echo_string,echo_file);
	}
}

#define LAGR_MULT 1

  iread = look_for_optional(ifp,"Phase Function Constraint Method",input,'=');
  if (iread == 1) 
  {   
    if ( fscanf(ifp,"%s",input) != 1)
      {
	EH( -1, "error reading Phase Function Initialization string");
      }          
    stringup(input);

    SPF(echo_string,eoformat, "Phase Function Constraint Method", input);
    
    if( ( strcmp(input,"NO" ) == 0 ) || ( strcmp(input,"NONE") == 0 ) )
      {
	pfd->Constraint_Method = FALSE;
      }
    else if( ( strcmp(input, "LM" ) == 0 ) )
      {
	pfd->Constraint_Method = LAGR_MULT;
      }
    else 
      {
	EH(-1,"Phase function Constraint Method not recognized.\n");
      }
	  
    if( pfd->Constraint_Method == LAGR_MULT )
      {
	
	pfd->Use_Constraint = TRUE;
	    
	if ( augc == NULL )
	  {
		  pfd->Constraint_Method = FALSE;
		  pfd->Use_Constraint = FALSE;
		  pfd->jac_info = NULL;
	  }
	else
	  {
	    EH(-1,"LM phase function constraint method must be first augmenting condition.\n");
	  }
	
	augc[nAC-1].Type = AC_PF_CONSTRAINT;
	
	pfd->jac_info = ( PF_JAC_INFO *) array_alloc( 1, 1, sizeof( PF_JAC_INFO ) );
	
	if( fscanf(ifp,"%d %lf", &(augc[nAC-1].MTID), &(augc[nAC-1].CONSTV) ) != 2 )
	  {
	    EH(-1,"LM constraint type for phase function needs integer and double arguments.\n");
	  }
	
	SPF(endofstring(echo_string)," %d %.4g", augc[nAC-1].MTID, augc[nAC-1].CONSTV); ECHO(echo_string,echo_file);
      }
  }

iread = look_for_optional(ifp,"Phase Function Initial Shift",input,'=');
if (iread == 1) 
  {   
	double *tmp;

	if( ( read_constants( ifp, &tmp, 0 ) != pfd->num_phase_funcs ) )
	{
		EH(-1,"Inconsistent number of phase function shift values \n.");
	}

	for(i=0; i<pfd->num_phase_funcs; i++ ) pfd->shift[i] = tmp[i]; 

	SPF(echo_string,"%s = ","Phase Function Initial Shift"); 

	for(i=0; i<pfd->num_phase_funcs; i++)
	  {
	    SPF(endofstring(echo_string),"%.4g ", pfd->shift[i]);
	  }

	ECHO(echo_string,echo_file);
	
	safe_free( (void *) tmp );
	
  }

}


/*
 * rd_elem_quality_specs -- read input file for element quality specifications
 *
 * Comments:    This code was lifted out of the ever-growing read_input_file()
 *              routine above. This makes things more modular, with division
 *              into "sections".
 */
void 
rd_elem_quality_specs(FILE *ifp,
                      char *input )
{
  const char yo[] = "rd_elem_quality_specs";
  
  int iread;
  nEQM = 0;
  
  /* Element quality specifications */
  
  eqm->do_jac = FALSE;
  eqm->wt_jac = 0.0;
  iread = look_for_optional
    (ifp,"Jacobian quality weight",input,'=');
  if (iread == 1)
    {   
      nEQM++;
      eqm->do_jac = TRUE;
      if ( fscanf(ifp,"%lf",&(eqm->wt_jac)) != 1)
        {
          EH( -1, "error reading Jacobian quality weight");
        }
    }
  
  eqm->do_vol = FALSE;
  eqm->wt_vol = 0.0;
  iread = look_for_optional
    (ifp,"Volume change quality weight",input,'=');
  if (iread == 1)
    {   
      nEQM++;
      eqm->do_vol = TRUE;
      if ( fscanf(ifp,"%lf",&(eqm->wt_vol)) != 1)
        {
          EH( -1, "error reading Volume change quality weight");
        }
    }
  
  eqm->do_ang = FALSE;
  eqm->wt_ang = 0.0;
  iread = look_for_optional
    (ifp,"Angle quality weight",input,'=');
  if (iread == 1)
    {   
      nEQM++;
      eqm->do_ang = TRUE;
      if ( fscanf(ifp,"%lf",&(eqm->wt_ang)) != 1)
        {
          EH( -1, "error reading Angle quality weight");
        }
    }
  
  eqm->do_tri = FALSE;
  eqm->wt_tri = 0.0;
  iread = look_for_optional
    (ifp,"Triangle quality weight",input,'=');
  if (iread == 1)
    {   
      nEQM++;
      eqm->do_tri = TRUE;
      if ( fscanf(ifp,"%lf",&(eqm->wt_tri)) != 1)
        {
          EH( -1, "error reading Triangle quality weight");
        }
    }
  
  eqm->eq_tol = 0.0;
  iread = look_for_optional
    (ifp,"Element quality tolerance",input,'=');
  if (iread == 1)
    {   
      eqm->do_jac = TRUE;
      if ( fscanf(ifp,"%lf",&(eqm->eq_tol)) != 1)
        {
          EH( -1, "error reading Element quality tolerance");
        }
    }
  
  iread=look_for_optional(ifp,"Element quality tolerance type",input,'=');
  
  if(iread == 1)
    {
      (void) read_string(ifp,input,'\n');
      strip(input);
      if ( strcmp(input,"jac") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tMinimum Jacobian quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_JAC;
	}
      else if ( strcmp(input,"vol") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tMinimum Volume change quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_VOL;
	}
      else if ( strcmp(input,"ang") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tMinimum Angle quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_ANG;
	}
      else if ( strcmp(input,"tri") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tMinimum Triangle quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_TRI;
	}
      else if ( strcmp(input,"min") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tMinimum overall quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_MIN;
	}
      else if ( strcmp(input,"avg") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tWeighted average quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_AVG;
	}
      else if ( strcmp(input,"wtmin") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tWeighted minimum quality tolerance\n", yo);
	    }
	  eqm->tol_type = EQM_WTMIN;
	}
      else
        {
          EH( -1, "error reading Element quality tolerance type");
        }
    }
  else
    {
      eqm->tol_type = EQM_MIN;
    }
}

/*
 * rd_track_specs -- read input file for continuation specifications
 *
 * Comments:    This code was lifted out of the ever-growing read_input_file()
 *              routine above. This makes things more modular, with division
 *              into "sections".
 *
 * Adapted by Ian Gates
 */

void 
rd_track_specs(FILE *ifp,
	       char *input)
{
  char *yo;
  
  int mn, iread;
  int TPContType=0, TPBdyCondID, TPDataFltID, TPMatID, TPMatPropID, TPMatPropSIID;
  int id1, id2, id3, iflag, iCC, iTC;
  double range;
  double beg_angle = 0.0, end_angle = 0.0;
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;
  
  yo = "rd_track_specs";
  

  /* Read in Continuation Specifications */

  iread = look_for_optional(ifp,"Continuation Specifications",input,'=');
  
  if( iread == 1 ) {
	  SPF(echo_string,"%s =",input); ECHO(echo_string,echo_file);
  }

  Continuation = ALC_NONE;
  loca_in->Cont_Alg = ALC_NONE;
  loca_in->Cont_Order = 0;
  nCC = 0;
  nTC = 0;
  cont->eps = 0.0;

  /*  initialize sensivity vector id for first order continuation  */

  cont->sensvec_id = -1;
  loca_in->NVSave = FALSE;

  iread=look_for_optional(ifp,"Continuation",input,'=');

  if(iread == 1)
    {
	  SPF(echo_string,"%s = ", input);
	  
      (void) read_string(ifp,input,'\n');
      strip(input);
      if ( strcmp(input,"none") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tno continuation\n", yo);
	    }
	  Continuation = ALC_NONE;
	}
      else
	if ( strcmp(input,"zero") == 0 )
	  {
	    if ( Debug_Flag && ProcID == 0 )
	      {
		printf("%s:\tzero order continuation\n", yo);
	      }
	    Continuation = ALC_ZEROTH;
	  }
	else
	  if ( strcmp(input,"hzero") == 0 )
	    {
	      if ( Debug_Flag && ProcID == 0 )
		{
		  printf("%s:\tzero order hunting continuation\n", yo);
		}
	      Continuation = HUN_ZEROTH;
	    }
	  else
	    if ( strcmp(input,"first") == 0 )
	      {
		if ( Debug_Flag && ProcID == 0 )
		  {
		    printf("%s:\tfirst order continuation\n", yo);
		  }
		Continuation = ALC_FIRST;
	      }
	    else
	      if ( strcmp(input,"hfirst") == 0 )
		{
		  if ( Debug_Flag && ProcID == 0 )
		    {
		      printf("%s:\tfirst order hunting continuation\n", yo);
		    }
		  Continuation = HUN_FIRST;
		}
	      else
		if ( strcmp(input,"second") == 0 )
		  {
		    if ( Debug_Flag && ProcID == 0 )
		      {
			printf("%s:\tsecond order continuation\n", yo);
		      }
		    Continuation = ALC_SECOND;
		  } else
                  if ( strcmp(input,"loca") == 0 )
                    {
                      if ( Debug_Flag && ProcID == 0 )
                        {
                          printf("%s:\tLOCA continuation library\n", yo);
                        }
                      Continuation = LOCA;
                    }
                  else
                    {
                      EH( -1, "unknown Continuation option");
                    }
	  SPF(endofstring(echo_string)," %s", input); ECHO(echo_string,echo_file);
      
      ContType = cont->upType = -1;
      look_for(ifp,"Continuation Type",input,'=');
	  
	  SPF(echo_string,"%s = ",input);
	  
      (void) read_string(ifp,input,'\n');
      strip(input);
      if ( strcmp(input,"none") == 0 )
	{
	  if ( Debug_Flag && ProcID == 0 )
	    {
	      printf("%s:\tno continuation\n", yo);
	    }
	  Continuation = ALC_NONE;
	}
      else
	if ( strcmp(input,"BC") == 0 )
	  {
	    if ( Debug_Flag && ProcID == 0 )
	      {
		printf("%s:\tBC continuation\n", yo);
	      }
	    ContType = 1;
	  }
	else
	  if ( strcmp(input,"MT") == 0 )
	    {
	      if ( Debug_Flag && ProcID == 0 )
		{
		  printf("%s:\tMATPROP continuation\n", yo);
		}
	      ContType = 2;
	    }
       else
	if ( strcmp(input,"AC") == 0 )
 	  {
 	    if ( Debug_Flag && ProcID == 0 )
 	      {
 		printf("%s:\tAC continuation\n", yo);
 	      }
 	    ContType = 3;
 	  }
	else
	  if ( strcmp(input,"UM") == 0 )
	    {
	      if ( Debug_Flag && ProcID == 0 )
		{
		  printf("%s:\tUser model continuation\n", yo);
		}
	      ContType = 4;
	    }
	else
	  if ( strcmp(input,"UF") == 0 )
	    {
	      if ( Debug_Flag && ProcID == 0 )
		{
		  printf("%s:\tUser-defined function continuation\n", yo);
		}
	      ContType = 5;
	    }
	else
	  if ( strcmp(input,"AN") == 0 )
	    {
	      if ( Debug_Flag && ProcID == 0 )
		{
		  printf("%s:\tAngular parameter continuation\n", yo);
		}
	      ContType = 6;
	    }
	  else
	    {
	      EH( -1, "Unknown ContType option (one of BC, MT, AC, UM, AN or UF)");
		}
	  
	  SPF(endofstring(echo_string)," %s", input); ECHO(echo_string,echo_file);
	  
      cont->upType = ContType;
      
      for (mn=0; mn<MAX_NUMBER_MATLS; mn++) {
	pd_glob[mn]->Continuation =  Continuation;
      }

/* Get number of user-defined continuation condition functions */
      nUC = 0;
      if (ContType == 5)
        {
          look_for(ifp, "Number of user continuation functions",input,'=');
          if (fscanf(ifp,"%d",&nUC) != 1)
            {
              EH(-1, "error reading Number of user continuation functions");
            }
		  
		  SPF(echo_string,"%s = %d",input, nUC); ECHO(echo_string,echo_file);

/* Allocate the structure cpuc for each function */
          cpuc = (struct User_Continuation_Info *)
                  array_alloc(1, nUC, sizeof(struct User_Continuation_Info));
        }
      
      cont->upBCID = -1;
      iread = look_for_optional(ifp,"Boundary condition ID",input,'=');
      if (iread == 1)
	  {   
		  if ( fscanf(ifp,"%d",&BdyCondID) != 1)
		  {
			  EH( -1, "error reading Boundary condition ID");
		  }
		  SPF(echo_string,"%s = %d", input,BdyCondID); ECHO(echo_string,echo_file);
	  }
      cont->upBCID = BdyCondID;
      
      look_for(ifp,"Boundary condition data float tag",input,'=');
      if (fscanf (ifp,"%d",&DataFltID) != 1)
	{
	  EH( -1, "error reading Boundary condition data float tag id");
	}
	  SPF(echo_string,"%s = %d", input, DataFltID); ECHO(echo_string,echo_file);
	  
      cont->upDFID = DataFltID;
      
      look_for(ifp,"Material id",input,'=');
      if (fscanf (ifp,"%d",&MatID) != 1)
	{
	  EH( -1, "error reading Material id");
	}

      cont->upMTID = MatID-1;
	  SPF(echo_string,"%s = %d", input, MatID); ECHO(echo_string,echo_file);
      
      look_for(ifp,"Material property tag",input,'=');
      if (fscanf (ifp,"%d",&MatPropID) != 1)
	{
	  EH( -1, "error reading Material property tag id");
	}
      cont->upMPID = MatPropID;
	  SPF(echo_string,"%s = %d", input, MatPropID); ECHO(echo_string,echo_file);
      
      look_for(ifp,"Material property tag subindex",input,'=');
      if (fscanf (ifp,"%d",&MatPropSIID) != 1)
	{
	  EH( -1, "error reading Material property tag subindex id");
	}
      cont->upMDID = MatPropSIID;
	  SPF(echo_string,"%s = %d", input, MatPropSIID); ECHO(echo_string,echo_file);
      
      BegParameterValue = 0.0;
      look_for(ifp,"Initial parameter value",input,'=');
      if ( fscanf(ifp,"%le",&BegParameterValue) != 1)
	{
	  EH( -1, "error reading BegParameterValue");
	}
      cont->BegParameterValue = BegParameterValue;
	  SPF(echo_string,"%s = %.4g", input, BegParameterValue); ECHO(echo_string,echo_file);
      
      EndParameterValue = 1.0;
      look_for(ifp,"Final parameter value",input,'=');
      if ( fscanf(ifp,"%le",&EndParameterValue) != 1)
	{
	  EH( -1, "error reading EndParameterValue");
	}
      cont->EndParameterValue = EndParameterValue;
	  SPF(echo_string,"%s = %.4g", input, EndParameterValue); ECHO(echo_string,echo_file);
      
      look_for(ifp,"delta_s",input,'=');
      if (fscanf (ifp,"%le",&Delta_s0) != 1)
	{
	  EH( -1, "error reading delta_s");
	}
      cont->Delta_s0 = Delta_s0;
	  SPF(echo_string,"%s = %.4g", input, Delta_s0); ECHO(echo_string,echo_file);
      
      look_for(ifp,"Maximum number of path steps",input,'=');
      if (fscanf (ifp,"%d",&MaxPathSteps) != 1)
	{
	  EH( -1, "error reading maximum number of path steps");
	}
      cont->MaxPathSteps = MaxPathSteps;
	  SPF(echo_string,"%s = %d", input, MaxPathSteps); ECHO(echo_string,echo_file);
      
      look_for(ifp,"Minimum path step",input,'=');
      if ( fscanf(ifp,"%le",&Delta_s_min) != 1)
	{
	  EH( -1, "error reading minimum path step");
	}
      cont->Delta_s_min = Delta_s_min;
	  SPF(echo_string,"%s = %.4g", input, Delta_s_min); ECHO(echo_string,echo_file);
      
      Delta_s_max = 1.e+12;
      look_for(ifp,"Maximum path step",input,'=');
      if ( fscanf(ifp,"%le",&Delta_s_max) != 1)
	{
	  EH( -1, "error reading Maximum path step");
	}
      cont->Delta_s_max = Delta_s_max;
	  SPF(echo_string,"%s = %.4g", input, Delta_s_max); ECHO(echo_string,echo_file);
      
      look_for(ifp,"Continuation Printing Frequency",input,'=');
      if ( fscanf(ifp,"%d",&print_freq) != 1)
	{
	  EH( -1, "error reading Continuation Printing Frequency");
	}
      cont->print_freq = print_freq;
	  SPF(echo_string,"%s = %d", input, print_freq); ECHO(echo_string,echo_file);

      cont->fix_freq = 0;
      if (Num_Proc > 1) {
        iread = look_for_optional(ifp,"Continuation Fix Frequency",input,'=');
        if (iread == 1) {
          cont->fix_freq = read_int(ifp, "Continuation Fix Frequency");
          if (cont->fix_freq < 0) {
            EH(-1, "Expected Fix Frequency > 0");
          }
          SPF(echo_string, "%s = %d", "Continuation Fix Frequency", cont->fix_freq); ECHO(echo_string, echo_file);
        }
      }


      if (Continuation == LOCA)
        {
        
          loca_in->Cont_Alg = CONTINUATION;
          iread=look_for_optional(ifp,"LOCA method",input,'=');
          if (iread == 1)
		  {
			  SPF(echo_string,"%s =",input);
			  
              (void) read_string(ifp,input,'\n');
              strip(input);
			  
              if ( strcmp(input,"ss") == 0 )
			  {
                  if ( Debug_Flag && ProcID == 0 )
				  {
                      printf("%s:\tsteady state continuation\n", yo);
				  }
				  
				  SPF(endofstring(echo_string)," %s",input);
				  
                  loca_in->Cont_Alg = CONTINUATION;
                  iread = look_for_optional(ifp,"Continuation order",input,'=');
                  if (iread == 1)
				  {
                      if ( fscanf(ifp,"%d",&(loca_in->Cont_Order)) != 1)
					  {
                          EH( -1, "error reading continuation order");
					  }
					  SPF(endofstring(echo_string),"/n%s = %d",input, loca_in->Cont_Order);
				  }
			  }
			  
              else
				  if ( strcmp(input,"zero") == 0 )
                  {
					  if ( Debug_Flag && ProcID == 0 )
                      {
						  printf("%s:\tzero order continuation\n", yo);
                      }
					  SPF(endofstring(echo_string)," %s",input);
                  }
			  
			  else
                  if ( strcmp(input,"first") == 0 )
				  {
                      if ( Debug_Flag && ProcID == 0 )
					  {
                          printf("%s:\tzero order continuation\n", yo);
					  }
                      loca_in->Cont_Order = 1;
					  SPF(endofstring(echo_string)," %s",input);
				  }
			  
			  else
				  if ( strcmp(input,"alc") == 0 )
				  {
					  if ( Debug_Flag && ProcID == 0 )
					  {
						  printf("%s:\tarc length continuation\n", yo);
					  }
					  loca_in->Cont_Order = 2;
					  SPF(endofstring(echo_string)," %s",input);
				  }
			  
			  else
				  if ( strcmp(input,"tp") == 0 )
				  {
					  if ( Debug_Flag && ProcID == 0 )
					  {
						  printf("%s:\tturning point continuation\n", yo);
					  }
					  loca_in->Cont_Alg = TP_CONTINUATION;
					  SPF(endofstring(echo_string)," %s",input);
				  }
			  
			  else
				  if ( strcmp(input,"pf") == 0 )
				  {
					  if ( Debug_Flag && ProcID == 0 )
					  {
						  printf("%s:\tpitchfork continuation\n", yo);
					  }
					  loca_in->Cont_Alg = PF_CONTINUATION;
					  SPF(endofstring(echo_string)," %s",input);
				  }
			  
			  else
				  if ( strcmp(input,"hp") == 0 )
				  {
					  if ( Debug_Flag && ProcID == 0 )
					  {
						  printf("%s:\tHopf continuation\n", yo);
					  }
					  loca_in->Cont_Alg = HP_CONTINUATION;
					  SPF(endofstring(echo_string)," %s",input);
				  }
			  
			  else
			  {
				  EH( -1, "Invalid LOCA method!");
			  }
		  }
		  
		  ECHO(echo_string,echo_file);

      
          loca_in->StepAggr = 0.0;
          iread = look_for_optional(ifp,"Step control aggressiveness",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->StepAggr)) != 1)
                {
                  EH( -1, "error reading step control aggressiveness");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->StepAggr); ECHO(echo_string,echo_file);
            }

 /* Use negative input for constant step size = range / (MaxSteps-1) */
          if (loca_in->StepAggr < -1.0e-6)
            {
              loca_in->StepAggr = 0.0;
              if (loca_in->Cont_Order != 2) cont->Delta_s0 = 
                (cont->EndParameterValue - cont->BegParameterValue)
                / ((double)(cont->MaxPathSteps - 1));
            }

          loca_in->perturb = 1.0e-9;
          iread = look_for_optional(ifp,"Perturbation magnitude",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->perturb)) != 1)
                {
                  EH( -1, "error reading Perturbation magnitude");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->perturb); ECHO(echo_string,echo_file);
            }

          loca_in->debug = 5;
          iread = look_for_optional (ifp,"LOCA print level",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%d",&(loca_in->debug)) != 1)
                {
                  EH( -1, "error reading LOCA print level");
                }
			  SPF(echo_string,"%s = %d", input, loca_in->debug); ECHO(echo_string,echo_file);
            }
      
          loca_in->DpDs2 = 0.5;
          iread = look_for_optional (ifp,"ALC Desired solution fraction",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->DpDs2)) != 1)
                {
                  EH( -1, "error reading ALC Desired solution fraction");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->DpDs2); ECHO(echo_string,echo_file);
            }
      
          loca_in->DpDsHi = 1.0;
          iread = look_for_optional (ifp,"ALC Max. parameter sensitivity",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->DpDsHi)) != 1)
                {
                  EH( -1, "error reading ALC Max. parameter sensitivity");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->DpDsHi); ECHO(echo_string,echo_file);
            }
      
          loca_in->Texp = 0.0;
          iread = look_for_optional (ifp,"ALC Tangent factor exponent",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->Texp)) != 1)
                {
                  EH( -1, "error reading ALC Tangent factor exponent");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->Texp); ECHO(echo_string,echo_file);
            }

          loca_in->MaxTS = 0.0;
          iread = look_for_optional(ifp,"ALC Tangent factor step limit",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->MaxTS)) != 1)
                {
                  EH( -1, "error reading ALC Tangent factor step limit");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->MaxTS); ECHO(echo_string,echo_file);			  
			}

  /* Bypass CC card section if LOCA is to use existing Hunting cards
     - this is for backward compatibility - EDW */
          nCC = 1;
          iread = look_for_optional(ifp,"Number of continuation conditions",input,'=');
          if (iread == 1)
            {
              if (fscanf(ifp,"%d",&nCC) != 1)
                {
                  EH( -1, "error reading Number of continuation conditions");
                }

              if (nCC == -1)
                {
                  nCC = count_list(ifp,"CC",input,'=',"END OF CC") + 1;
                }
			  

              if (nCC < 1 && nCC != -2)
                {
                  Continuation = ALC_NONE;
				  SPF(echo_string,"\t(%s)","No continuation conditions specified");
                  return;
                }
			  else
				  SPF(echo_string,"%s = %d", input, nCC);
			  
			   ECHO(echo_string,echo_file);
            }

  /*
   *  For angular continuation, ensure at least one angular CC is specified
   */
          if (cont->upType == 6 && nCC < 2)
            {
              EH(-1, "Angular continuation requires at least one CC card!");
            }

  /*
   *  Allocate space for the vector cpcc
   */
          if (nCC != -2)
            {

              cpcc = (struct Continuation_Conditions *)
                calloc(nCC, sizeof(struct Continuation_Conditions));

              if ( Debug_Flag && ProcID == 0 )
                {
                  printf("%s:\tallocated %d copies of the %lu sized\n", 
                         yo, nCC, (long unsigned int)sizeof(struct Continuation_Conditions));
                  printf("%s:\tCC structure.\n", yo);
                }

  /*
   *  Load up cpcc[0] with inputs from above
   *  For angular continuation, also get angle range in radians
   */
              cpcc[0].nCC = nCC;
              cpcc[0].ratio = 1.0;
              cpcc[0].Type = cont->upType;
              cpcc[0].BCID = cont->upBCID;
              cpcc[0].DFID = cont->upDFID;
              cpcc[0].MTID = cont->upMTID;
              cpcc[0].MPID = cont->upMPID;
              cpcc[0].MDID = cont->upMDID;
              cpcc[0].Beg_CC_Value = cont->BegParameterValue;
              cpcc[0].End_CC_Value = cont->EndParameterValue;
              range = cpcc[0].End_CC_Value - cpcc[0].Beg_CC_Value;
              if (cont->upType == 6)
                {
                  beg_angle = BegParameterValue * M_PIE / 180.0;
                  end_angle = EndParameterValue * M_PIE / 180.0;
                }
              cpcc[0].fn_flag = 0;
              cpcc[0].old_value = cpcc[0].Beg_CC_Value;
              cpcc[0].value = cpcc[0].Beg_CC_Value;

  /*
   *  Read in continuation parameter conditions and load up cpcc
   */

/* First check for UF/CC conflict (UF has priority) */
	      if (nCC > 1 && cont->upType == 5)
                {
			  SPF(echo_string,"\t(%s)","UF continuation overrides CC cards!");	ECHO(echo_string,echo_file);
                  nCC = 1;
                }

              if (nCC > 1)
                {
                  for (iCC=1; iCC<nCC; iCC++)
                    {

                  look_for(ifp, "CC", input, '=');
					  
					  SPF(echo_string,"%s = ", input);
					  
                  if (fscanf(ifp, "%80s", input) != 1)
                    {
                      fprintf (stderr, "%s:\tError reading CC[iCC].Set_Type ",yo);
                      fprintf (stderr, "(one of BC, MT, UM, or AC)\n");
                      exit (-1);
                    }
                  if (!strcmp(input,"BC"))
                    {
                      cpcc[iCC].Type = 1;
                    }
                  else if (!strcmp(input,"MT"))
                    {
                      cpcc[iCC].Type = 2;
                    }
                  else if (!strcmp(input,"AC"))
                    {
                      cpcc[iCC].Type = 3;
                    }
                  else if (!strcmp(input,"UM"))
                    {
                      cpcc[iCC].Type = 4;
                    }
                  else
                    {
                      fprintf(stderr, "%s:\tImproper CC set type - %s\n", yo, input);
                      if (!strcmp(input,"UF")) fprintf(stderr, "\tUM is not a valid option for CC cards!");
                      if (!strcmp(input,"AN")) fprintf(stderr, "\tAN is not a valid option for CC cards!");
                      exit (-1);
                    }
				  SPF(endofstring(echo_string)," %s",input);

                  cpcc[iCC].nCC = nCC;
  /*
   *  Read in three required integers (four if type is UM)
   */

                  if (cpcc[iCC].Type == 4)
                    {
                      if (fscanf(ifp, "%d %d %d %d", &id1, &id2, &id3, &iflag) != 4)
                        {
                          fprintf(stderr, "%s:\tError reading integer inputs, ",yo);
                          fprintf(stderr, "for CC[%d]\n",iCC+1);
                          fprintf(stderr, "\tFormat: CC = UM MTID MPID MDID FLAG ...\n");
                          exit(-1);
                        }
					  SPF(endofstring(echo_string)," %d %d %d %d", id1, id2, id3, iflag);

                    }
                  else 
				  {
					  if (fscanf(ifp, "%d %d %d", &id1, &id2, &iflag) != 3)
					  {
						  fprintf(stderr, "%s:\tError reading integer inputs ",yo);
						  fprintf(stderr, "for CC[%d]\n",iCC+1);
						  if (cpcc[iCC].Type == 1)
						  {
							  fprintf(stderr,"\tFormat: CC = BC BCID DFID FLAG ...\n");
						  }
						  if (cpcc[iCC].Type == 2)
						  {
							  fprintf(stderr,"\tFormat: CC = MT MTID MPID FLAG ...\n");
						  }
						  if (cpcc[iCC].Type == 3)
						  {
							  fprintf(stderr,"\tFormat: CC = AC ACID DFID FLAG ...\n");
						  }
						  exit (-1);
					  }
					  SPF(endofstring(echo_string)," %d %d %d", id1, id2, iflag);
                    }
				  /*
				   *  For continuation types other than "AN":
				   *  Third int input FLAG indicates meaning of second float input "vfloat":
				   *  0 - Value is same as continuation parameter, no floats needed
				   *  1 - vfloat is End_CC_Value[iCC]; use to find range ratio d(val[iCC])/d(val[0])
				   *  2 - vfloat is range ratio; use to find End_CC_Value[iCC]
				   *  SPECIAL CASE:
				   *  3 - Value is constructed from three floats:
				   *        value = coeff_0 + coeff_1 * lambda^coeff_2
				   */
				  
				  if (cont->upType != 6 && iflag == 0)
				  {
					  cpcc[iCC].ratio = 1.0;
					  cpcc[iCC].Beg_CC_Value = cpcc[0].Beg_CC_Value;
					  cpcc[iCC].End_CC_Value = cpcc[0].End_CC_Value;
				  }
				  else if (iflag == 3)  /* Three floats required */
				  {
					  if (fscanf(ifp, "%lf %lf %lf",
								 &cpcc[iCC].coeff_0, &cpcc[iCC].coeff_1, &cpcc[iCC].coeff_2) != 3)
					  {
						  fprintf(stderr, "%s:\tError reading CC[%d] floats\n",
								  yo, iCC+1);
						  fprintf(stderr, "\tThree floats are required!\n");
						  exit (-1);
					  }
					  SPF(endofstring(echo_string)," %.4g %.4g %.4g", cpcc[iCC].coeff_0, cpcc[iCC].coeff_1, cpcc[iCC].coeff_2);
				  }
				  else  /* Two floats are required */
				  {
					  if (fscanf(ifp, "%lf %lf",
								 &cpcc[iCC].coeff_0, &cpcc[iCC].coeff_1) != 2)
					  {
						  fprintf(stderr, "%s:\tError reading CC[%d] floats\n",
								  yo, iCC+1);
						  fprintf(stderr, "\tFormat: ...start end/range\n");
						  exit (-1);
					  }
					  SPF(endofstring(echo_string), " %.4g %.4g", cpcc[iCC].coeff_0, cpcc[iCC].coeff_1);
				  }
				  
				  ECHO(echo_string,echo_file);
				  
                  switch (cpcc[iCC].Type) {

                  case 1: /* BC */
                    cpcc[iCC].BCID = id1;
                    cpcc[iCC].DFID = id2;
                    SPF(echo_string, "\t( %3d. BC: BCID=%3d DFID=%5d)",iCC+1,
                            cpcc[iCC].BCID, cpcc[iCC].DFID);
                    break;

                  case 2: /* MT */
                    cpcc[iCC].MTID = id1;
                    cpcc[iCC].MPID = id2;
                    SPF(echo_string, "\t( %3d. MT: MTID=%3d MPID=%5d)",iCC+1,
                            cpcc[iCC].MTID, cpcc[iCC].MPID);
                    cpcc[iCC].MTID--;
                    break;

                  case 3: /* AC */
                    cpcc[iCC].BCID = id1;
                    cpcc[iCC].DFID = id2;
                    SPF(echo_string, "\t( %3d. AC: ACID=%3d DFID=%5d)",iCC+1,
                            cpcc[iCC].BCID, cpcc[iCC].DFID);
                    break;

                  case 4: /* UM */
                    cpcc[iCC].MTID = id1;
                    cpcc[iCC].MPID = id2;
                    cpcc[iCC].MDID = id3;
                    SPF(echo_string, "\t( %3d. UM: MTID=%3d MPID=%5d MDID=%3d)",
                            iCC+1,cpcc[iCC].MTID, cpcc[iCC].MPID, cpcc[iCC].MDID);
                    cpcc[iCC].MTID--;
                    break;
                  }
				  
				  ECHO(echo_string,echo_file);

				  
				  /*
   * The floats have different purposes in angular continuation,
   * so they will be stored in members coeff_0, coeff_1, and coeff_2.
   * Parameter value is then given by coeff_0 + coeff_1 * f(theta), where:
   *    fn_flag = 0:   f(theta) = sin(theta)
   *    fn_flag = 1:   f(theta) = cos(theta)
   *    fn_flag = 2:   f(theta) = tan(theta)
   * SPECIAL CASE:
   *    fn_flag = 3:   value = coeff_0 + coeff_1 * sin(theta) + coeff_2 * cos(theta)
   */
                      cpcc[iCC].fn_flag = iflag;

                      if (cont->upType == 6)
                        {
                          cpcc[iCC].ratio = -1.0;

                          switch (iflag) {
                          case 0:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(beg_angle);
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(end_angle);
                            break;
                          case 1:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * cos(beg_angle);
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * cos(end_angle);
                            break;
                          case 2:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * tan(beg_angle);
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * tan(end_angle);
                            break;
                          case 3:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(beg_angle)
                              + cpcc[iCC].coeff_2 * cos(beg_angle);
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1 * sin(end_angle)
                              + cpcc[iCC].coeff_2 * cos(end_angle);
                            break;
                          default:
                            fprintf(stderr, "%s:\tCC[%d] flag must be 0-3\n",
                                    yo, iCC+1);
                            exit (-1);
                            break;
                          }
                        }

  /* Other cases are handled here: */
                      else
                        {
                          switch (iflag) {
                          case 0:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0;
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_1;
                            break;
                          case 1:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0;
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_1;
                            cpcc[iCC].ratio =
                              (cpcc[iCC].coeff_1 - cpcc[iCC].coeff_0) / range;
                            break;
                          case 2:
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0;
                            cpcc[iCC].ratio = cpcc[iCC].coeff_1;
                            cpcc[iCC].End_CC_Value = cpcc[iCC].Beg_CC_Value
                                                     + cpcc[iCC].ratio * range;
                            break;
                          case 3:
                            cpcc[iCC].ratio = -1.0;
                            cpcc[iCC].Beg_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1
                              * pow(BegParameterValue, cpcc[iCC].coeff_2);
                            cpcc[iCC].End_CC_Value = cpcc[iCC].coeff_0
                              + cpcc[iCC].coeff_1
                              * pow(EndParameterValue, cpcc[iCC].coeff_2);
                            /* fall through */
                          default:
                            fprintf(stderr, "%s:\tCC[%d] flag must be 0, 1, or 2\n",
                                    yo, iCC+1);
                            exit (-1);
                            break;
                          }

                        }  /* End of else block (cont->upType == 6) */
                   
                      cpcc[iCC].old_value = cpcc[iCC].Beg_CC_Value;
                      cpcc[iCC].value = cpcc[iCC].Beg_CC_Value;
                    }  /* End of iCC loop */

                }  /* End of continuation condition if block */

            }  /* End of nCC!=2 if block (temporary) */

        }  /* End of continuation=LOCA if block */

      if ( Continuation == LOCA &&
           ( loca_in->Cont_Alg == TP_CONTINUATION ||
             loca_in->Cont_Alg == PF_CONTINUATION ||
             loca_in->Cont_Alg == HP_CONTINUATION ) )
        {

          look_for(ifp,"TP Continuation Type",input,'=');
		  
		  SPF(echo_string,"%s =", input); 
		  
          (void) read_string(ifp,input,'\n');
          strip(input);
          if ( strcmp(input,"none") == 0 )
		  {
              EH(-1, "TP parameter type is required!");
		  }
          else
			  if ( strcmp(input,"BC") == 0 )
              {
				  if ( Debug_Flag && ProcID == 0 )
                  {
					  printf("%s:\tTP BC continuation\n", yo);
                  }
				  TPContType = 1;
              }
		  else
              if ( strcmp(input,"MT") == 0 )
			  {
                  if ( Debug_Flag && ProcID == 0 )
				  {
                      printf("%s:\tTP MATPROP continuation\n", yo);
				  }
                  TPContType = 2;
			  }
		  else
			  if ( strcmp(input,"AC") == 0 )
			  {
				  if ( Debug_Flag && ProcID == 0 )
				  {
					  printf("%s:\tTP AC continuation\n", yo);
				  }
				  TPContType = 3;
			  }
		  else
			  if ( strcmp(input,"UM") == 0 )
			  {
				  if ( Debug_Flag && ProcID == 0 )
				  {
					  printf("%s:\tTP User model continuation\n", yo);
				  }
				  TPContType = 4;
			  }
		  else
			  if ( strcmp(input,"UF") == 0 )
			  {
				  if ( Debug_Flag && ProcID == 0 )
				  {
					  printf("%s:\tTP User-defined function continuation\n", yo);
				  }
				  TPContType = 5;
			  }
		  else
			  if ( strcmp(input,"AN") == 0 )
			  {
				  if ( Debug_Flag && ProcID == 0 )
				  {
					  printf("%s:\tAngular parameter continuation\n", yo);
				  }
				  TPContType = 6;
			  }
		  else
		  {
			  EH( -1, "unknown TP ContType option (one of BC, MT, AC, UM, AN, or UF)");
		  }
		  
          loca_in->TPupType = TPContType;

		  SPF(endofstring(echo_string)," %s",input); ECHO(echo_string,echo_file);
		  
/* Get number of user-defined TP continuation condition functions */
          nUTC = 0;
          if (TPContType == 5)
            {
              look_for(ifp, "Number of user TP continuation functions",input,'=');
              if (fscanf(ifp,"%d",&nUTC) != 1)
                {
                  EH(-1, "error reading Number of user TP continuation functions");
                }

/* Allocate the structure tpuc for each function */
              tpuc = (struct User_Continuation_Info *)
                      array_alloc(1, nUTC, sizeof(struct User_Continuation_Info));
			  
			  SPF(echo_string,"%s = %d", input,nUTC); ECHO(echo_string,echo_file);
            }
      
          loca_in->TPupBCID = -1;
          iread = look_for_optional(ifp,"TP Boundary condition ID",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%d",&TPBdyCondID) != 1)
                {
                  EH( -1, "error reading TP Boundary condition ID");
                }    
			  SPF(echo_string,"%s = %d", input,TPBdyCondID); ECHO(echo_string,echo_file);
            }
          loca_in->TPupBCID = TPBdyCondID;
      
          look_for(ifp,"TP BC data float tag",input,'=');
          if (fscanf (ifp,"%d",&TPDataFltID) != 1)
            {
              EH( -1, "error reading TP BC data float tag id");
            }
          loca_in->TPupDFID = TPDataFltID;
		  SPF(echo_string,"%s = %d", input,TPDataFltID); ECHO(echo_string,echo_file);
 
          look_for(ifp,"TP parameter material id",input,'=');
          if (fscanf (ifp,"%d",&TPMatID) != 1)
            {
              EH( -1, "error reading TP parameter material id");
            }
          loca_in->TPupMTID = TPMatID-1;
		  SPF(echo_string,"%s = %d", input,TPMatID); ECHO(echo_string,echo_file);
 
          look_for(ifp,"TP parameter material property tag",input,'=');
          if (fscanf (ifp,"%d",&TPMatPropID) != 1)
            {
              EH( -1, "error reading TP parameter material property tag id");
            }
          loca_in->TPupMPID = TPMatPropID;
		  SPF(echo_string,"%s = %d", input,TPMatPropID); ECHO(echo_string,echo_file);
 
          look_for(ifp,"TP Material property tag subindex",input,'=');
          if (fscanf (ifp,"%d",&TPMatPropSIID) != 1)
            {
              EH( -1, "error reading TP Material property tag subindex id");
            }
          loca_in->TPupMDID = TPMatPropSIID;
		  SPF(echo_string,"%s = %d", input,TPMatPropSIID); ECHO(echo_string,echo_file);
 
          loca_in->TPGuess = 0.0;
          iread = look_for_optional(ifp,"Initial guess of TP parameter",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->TPGuess)) != 1)
                {
                  EH( -1, "error reading initial guess of TP parameter");
                }
            }
		  SPF(echo_string,"%s = %.4g", input,loca_in->TPGuess); ECHO(echo_string,echo_file);
 
          loca_in->TPFinal = 1.0;
          iread = look_for_optional(ifp,"TP parameter final value",input,'=');
          if (iread == 1)
            {   
              if ( fscanf(ifp,"%lf",&(loca_in->TPFinal)) != 1)
                {
                  EH( -1, "error reading TP parameter final value");
                }
            }
		  SPF(echo_string,"%s = %.4g", input,loca_in->TPFinal); ECHO(echo_string,echo_file);

          loca_in->NVRestart = FALSE;
          iread = look_for_optional(ifp, "Null vector restart file", input, '=');
          if (iread == 1)
            {
			  SPF(echo_string,"%s =",input);
              (void) read_string(ifp, input, '\n');
              strip(input);
              (void) strcpy(loca_in->NV_exoII_infile, input);

              loca_in->NVRestart = TRUE;
			  SPF(endofstring(echo_string)," %s", loca_in->NV_exoII_infile); ECHO(echo_string,echo_file);
            }
          else
            {
              if (loca_in->Cont_Alg == PF_CONTINUATION)
                {
                  EH(-1, "Pitchfork algorithm requires null vector input!");
                }
              if (loca_in->Cont_Alg == HP_CONTINUATION)
                {
                  EH(-1, "Hopf algorithm requires null vector input!");
                }
            }

  /* These two cards pertain only to Hopf continuation */
          if (loca_in->Cont_Alg == HP_CONTINUATION)
            {
              look_for(ifp, "Null vector imaginary restart file", input, '='); SPF(echo_string,"%s =",input);
		
              (void) read_string(ifp, input, '\n');
              strip(input);
              (void) strcpy(loca_in->NV_imag_infile, input);
			  
			  SPF(endofstring(echo_string)," %s", loca_in->NV_imag_infile); ECHO(echo_string,echo_file); 

              look_for(ifp, "Initial Hopf frequency", input, '=');
              if (fscanf(ifp, "%le", &loca_in->omega) != 1)
                {
                  EH(-1, "error reading initial Hopf frequency");
                }
			  SPF(echo_string,"%s = %.4g", input, loca_in->omega);  ECHO(echo_string,echo_file); 
            }

          iread = look_for_optional(ifp, "Null vector save file", input, '=');
          if (iread == 1)
            {
			  SPF(echo_string,"%s =", input);
              (void) read_string(ifp, input, '\n');
              strip(input);
              (void) strcpy(loca_in->NV_exoII_outfile, input);
			  SPF(endofstring(echo_string)," %s", loca_in->NV_exoII_outfile); ECHO(echo_string,echo_file);
              loca_in->NVSave = TRUE;
            }

          if (loca_in->Cont_Alg == HP_CONTINUATION && loca_in->NVSave)
            {
              iread = look_for_optional(ifp, "Null vector imaginary save file", input, '=');
              if (iread == 1)
                {
				  SPF(echo_string,"%s =", input);
                  (void) read_string(ifp, input, '\n');
                  strip(input);
                  (void) strcpy(loca_in->NV_imag_outfile, input);
                  SPF(endofstring(echo_string)," %s", loca_in->NV_imag_outfile); ECHO(echo_string,echo_file);
                }
            }
          nTC = 1;
          iread = look_for_optional(ifp,"Number of TP continuation conditions",input,'=');
          if (iread == 1)
            {
              if (fscanf(ifp,"%d",&nTC) != 1)
                {
                  EH( -1, "error reading Number of TP continuation conditions");
                }

              if (nTC == -1)
                {
                  nTC = count_list(ifp,"TC",input,'=',"END OF TC") + 1;
                }
			  
              SPF(echo_string,"%s = %d", input, nTC); ECHO(echo_string,echo_file);

              if (nTC < 1)
                {
                  Continuation = ALC_NONE;
                  SPF(echo_string,"\t(%s)", "No TP continuation conditions specified"); ECHO(echo_string,echo_file);
                  return;
                }

            }

  /*
   *  For angular continuation, ensure at least one angular TC is specified
   */
          if (TPContType == 6 && nTC < 2)
            {
              EH(-1, "Angular TP continuation requires at least one TC card!");
            }

  /*
   *  Allocate space for the vector tpcc
   */

          tpcc = (struct Continuation_Conditions *)
                  array_alloc(1, nTC, sizeof(struct Continuation_Conditions));

          if ( Debug_Flag && ProcID == 0 )
            {
              printf("%s:\tallocated %d copies of the %lu sized\n", 
                     yo, nTC, (long unsigned int)sizeof(struct Continuation_Conditions));
              printf("%s:\tTC structure.\n", yo);
            }

  /*
   *  Load up tpcc[0] with inputs from above
   *  For angular continuation, also get angle range in radians
   */
          tpcc[0].nCC = nTC;
          tpcc[0].ratio = 1.0;
          tpcc[0].Type = loca_in->TPupType;
          tpcc[0].BCID = loca_in->TPupBCID;
          tpcc[0].DFID = loca_in->TPupDFID;
          tpcc[0].MTID = loca_in->TPupMTID;
          tpcc[0].MPID = loca_in->TPupMPID;
          tpcc[0].MDID = loca_in->TPupMDID;
          tpcc[0].Beg_CC_Value = loca_in->TPGuess;
          tpcc[0].End_CC_Value = loca_in->TPFinal;
          range = tpcc[0].End_CC_Value - tpcc[0].Beg_CC_Value;
          if (TPContType == 6)
            {
              beg_angle = tpcc[0].Beg_CC_Value * M_PIE / 180.0;
              end_angle = tpcc[0].End_CC_Value * M_PIE / 180.0;
            }
          tpcc[0].fn_flag = 0;
          tpcc[0].old_value = tpcc[0].Beg_CC_Value;
          tpcc[0].value = tpcc[0].Beg_CC_Value;

  /*
   *  Read in TP parameter conditions and load up rest of tpcc
   */

/* First check for UF/TC conflict (UF has priority) */
	      if (nTC > 1 && loca_in->TPupType == 5)
                {
			  SPF(echo_string,"\t(%s)", "UF continuation overrides TC cards!"); ECHO(echo_string,echo_file);
                  nTC = 1;
                }

          if (nTC > 1)
            {

              for (iTC=1; iTC<nTC; iTC++)
                {

                  tpcc[iTC].nCC = nTC;

                  look_for(ifp, "TC", input, '='); SPF(echo_string,"%s =", input);
                  if (fscanf(ifp, "%80s", input) != 1)
                    {
                      fprintf (stderr, "%s:\tError reading TC[iTC].Set_Type ",yo);
                      fprintf (stderr, "(one of BC, MT, UM, or AC)\n");
                      exit (-1);
                    }
                  if (!strcmp(input,"BC"))
                    {
                      tpcc[iTC].Type = 1;
                    }
                  else if (!strcmp(input,"MT"))
                    {
                      tpcc[iTC].Type = 2;
                    }
                  else if (!strcmp(input,"AC"))
                    {
                      tpcc[iTC].Type = 3;
                    }
                  else if (!strcmp(input,"UM"))
                    {
                      tpcc[iTC].Type = 4;
                    }
                  else
                    {
                      fprintf(stderr, "%s:\tImproper TC set type - %s\n", yo, input);
                      if (!strcmp(input,"UF")) fprintf(stderr, "\tUM is not a valid option for TC cards!");
                      if (!strcmp(input,"AN")) fprintf(stderr, "\tAN is not a valid option for TC cards!");
                      exit (-1);
                    }
				  

  /*
   *  Read in three required integers (four if type is UM)
   */

                  if (tpcc[iTC].Type == 4)
                    {
                      if (fscanf(ifp, "%d %d %d %d", &id1, &id2, &id3, &iflag) != 4)
                        {
                          fprintf(stderr, "%s:\tError reading integer inputs, ",yo);
                          fprintf(stderr, "for TC[%d]\n",iTC+1);
                          fprintf(stderr, "\tFormat: TC = UM MTID MPID MDID FLAG ...\n");
                          exit(-1);
                        }
					  SPF(endofstring(echo_string)," %s %d %d %d %d", input, id1, id2, id3, iflag);

                    }
                  else 
                    {
					  if (fscanf(ifp, "%d %d %d", &id1, &id2, &iflag) != 3)
					  {
						  fprintf(stderr, "%s:\tError reading integer inputs ",yo);
						  fprintf(stderr, "for TC[%d]\n",iTC+1);
						  if (tpcc[iTC].Type == 1)
						  {
							  fprintf(stderr,"\tFormat: TC = BC BCID DFID FLAG ...\n");
						  }
						  if (tpcc[iTC].Type == 2)
						  {
							  fprintf(stderr,"\tFormat: TC = MT MTID MPID FLAG ...\n");
						  }
						  if (tpcc[iTC].Type == 3)
						  {
							  fprintf(stderr,"\tFormat: TC = AC ACID DFID FLAG ...\n");
						  }

						  exit (-1);
					  }
					  SPF(endofstring(echo_string)," %s %d %d %d", input, id1, id2, iflag);
					}

  /*
   *  For TP continuation types other than "AN":
   *  Third int input FLAG indicates meaning of second float input "vfloat":
   *  0 - Value is same as turning point parameter, no floats needed
   *  1 - vfloat is End_CC_Value[iTC]; use to find range ratio d(val[iTC])/d(val[0])
   *  2 - vfloat is range ratio; use to find End_CC_Value[iTC]
   *  SPECIAL CASE:
   *  3 - Value is constructed from three floats:
   *        value = coeff_0 + coeff_1 * lambda^coeff_2
   */

                  if (TPContType != 6 && iflag == 0)
                    {
                      tpcc[iTC].ratio = 1.0;
                      tpcc[iTC].Beg_CC_Value = tpcc[0].Beg_CC_Value;
                      tpcc[iTC].End_CC_Value = tpcc[0].End_CC_Value;
                    }
                  else if (iflag == 3)  /* Three floats required */
                    {
                      if (fscanf(ifp, "%lf %lf %lf",
                          &tpcc[iTC].coeff_0, &tpcc[iTC].coeff_1, &tpcc[iTC].coeff_2) != 3)
                        {
                          fprintf(stderr, "%s:\tError reading TC[%d] floats\n",
                                  yo, iTC+1);
                          fprintf(stderr, "\tThree floats are required!\n");
                          exit (-1);
                        }
					  SPF(endofstring(echo_string)," %.4g %.4g %.4g",tpcc[iTC].coeff_0, tpcc[iTC].coeff_1, tpcc[iTC].coeff_2);
                    }
                  else  /* Two floats are required */
                    {
                      if (fscanf(ifp, "%lf %lf",
                          &tpcc[iTC].coeff_0, &tpcc[iTC].coeff_1) != 2)
                        {
                          fprintf(stderr, "%s:\tError reading TC[%d] floats\n",
                                  yo, iTC+1);
                          fprintf(stderr, "\tFormat: ...start end/range\n");
                          exit (-1);
                        }
					  SPF(endofstring(echo_string)," %.4g %.4g",tpcc[iTC].coeff_0, tpcc[iTC].coeff_1);
                    }
				  
				  ECHO(echo_string,echo_file);
				  
				  switch (tpcc[iTC].Type)
				  {
					  
					  case 1: 
						  /* BC */
						  tpcc[iTC].BCID = id1;
						  tpcc[iTC].DFID = id2;
						  SPF(echo_string, "\t( %3d. BC: BCID=%3d DFID=%5d)",iTC+1,tpcc[iTC].BCID, tpcc[iTC].DFID);
						  break;
						  
					  case 2: 
						  /* MT */
						  tpcc[iTC].MTID = id1;
						  tpcc[iTC].MPID = id2;
						  SPF(echo_string, "\t( %3d. MT: MTID=%3d MPID=%5d)",iTC+1,
							  tpcc[iTC].MTID, tpcc[iTC].MPID);
						  tpcc[iTC].MTID--;
						  break;
						  
					  case 3: 
						  /* AC */
						  tpcc[iTC].BCID = id1;
						  tpcc[iTC].DFID = id2;
						  SPF(echo_string, "\t( %3d. AC: ACID=%3d DFID=%5d)",iTC+1,
							  tpcc[iTC].BCID, tpcc[iTC].DFID);
						  break;
						  
					  case 4: 
						  /* UM */
						  tpcc[iTC].MTID = id1;
						  tpcc[iTC].MPID = id2;
						  tpcc[iTC].MDID = id3;
						  SPF(echo_string, "\t( %3d. MT: MTID=%3d MPID=%5d MDID=%3d)",
								  iTC+1,tpcc[iTC].MTID,tpcc[iTC].MPID,tpcc[iTC].MDID);
						  tpcc[iTC].MTID--;
						  break;
                  }
				  
				  ECHO(echo_string,echo_file);

  /*
   * The two floats have different purposes in angular continuation,
   * see comments above for CC cards.
   */
                  tpcc[iTC].fn_flag = iflag;

                  if (TPContType == 6)
                    {
                      tpcc[iTC].ratio = -1.0;
                                                                                
                      switch (iflag) {
                      case 0:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * sin(beg_angle);
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * sin(end_angle);
                        break;
                      case 1:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * cos(beg_angle);
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * cos(end_angle);
                        break;
                      case 2:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * tan(beg_angle);
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * tan(end_angle);
                        break;
                      case 3:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * sin(beg_angle)
                          + tpcc[iTC].coeff_2 * cos(beg_angle);
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1 * sin(end_angle)
                          + tpcc[iTC].coeff_2 * cos(beg_angle);
                        break;
                      default:
                        fprintf(stderr, "%s:\tCC[%d] flag must be 0, 1, or 2\n",
                                yo, iTC+1);
                        exit (-1);
                        break;
                      }
                    }

  /* Other cases are handled here: */
                  else
                    {
                      switch (iflag) {
                      case 0:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0;
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_1;
                        break;
                      case 1:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0;
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_1;
                        tpcc[iTC].ratio =
                          (tpcc[iTC].coeff_1 - tpcc[iTC].coeff_0) / range;
                        break;
                      case 2:
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0;
                        tpcc[iTC].ratio = tpcc[iTC].coeff_1;
                        tpcc[iTC].End_CC_Value = tpcc[iTC].Beg_CC_Value
                                               + tpcc[iTC].ratio * range;
                        break;
                      case 3:
                        tpcc[iTC].ratio = -1.0;
                        tpcc[iTC].Beg_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1
                          * pow(BegParameterValue, tpcc[iTC].coeff_2);
                        tpcc[iTC].End_CC_Value = tpcc[iTC].coeff_0
                          + tpcc[iTC].coeff_1
                          * pow(EndParameterValue, tpcc[iTC].coeff_2);
                        break;
                      default:
                        fprintf(stderr, "%s:\tTC[%d] flag must be 0-3\n",
                                yo, iTC+1);
                        exit (-1);
                        break;
                      }

                    }  /* End of else block (TPContType == 6) */ 
                   
                  tpcc[iTC].old_value = tpcc[iTC].Beg_CC_Value;
                  tpcc[iTC].value = tpcc[iTC].Beg_CC_Value;
                }  /* End of iTC loop */

            }  /* End of TP continuation condition if block */

        }  /* End of TP/PF/HP_CONTINUATION if block */

    }  /* End of continuation if block (iread = 1) */
  
}
/* rd_track_specs -- read input file for continuation specifications */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * rd_hunt_specs -- read input file for hunting specifications
 *
 * Comments:    This code was lifted out of the ever-growing read_input_file()
 *              routine above. This makes things more modular, with division
 *              into "sections".
 *
 * Adapted by Ian (my coding style sucks) Gates 
 *
 * FORMAT:   
 * 
 * BOUNDARY CONDITION DF:  HC = BC BCID DFID RAMP from_value to_value step_value step_min step_max
 *
 * MATERIAL PROPERTY:      HC = MT MTID MPID RAMP from_value to_value step_value step_min step_max
 */

void 
rd_hunt_specs(FILE *ifp,
	      char *input)
{
  int iHC;
  static const char yo[] = "rd_hunt_specs";
  double range, range_0;
  int iread;
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  if ((Continuation != HUN_ZEROTH) && (Continuation != HUN_FIRST)
    && (Continuation != LOCA)) {
    nHC = 0;
    return; }

  /* Read in Hunting Specifications */

  iread = look_for_optional(ifp,"Hunting Specifications",input,'=');

  if (iread != 1) {
    nHC = 0;
    return; }
  else {
	  SPF(echo_string,"%s =", input); ECHO(echo_string,echo_file);
  }

  look_for(ifp,"Number of hunting conditions",input,'=');
  if (fscanf (ifp,"%d",&nHC) != 1)
    {
      EH( -1, "error reading Number of hunting conditions");
    }
  
  /*
   * Count hunting conditions if nHC is set to -1
   */

  if (nHC == -1)
    {
      nHC = count_list(ifp, "HC", input, '=', "END OF HC");
    }
  
  SPF(echo_string,"%s = %d",  input, nHC); ECHO(echo_string,echo_file);

  /*
   *  Allocate space for the vector hunt
   */

  if (nHC == 0) 
    {
	  Continuation = -1;
	  SPF(echo_string,"\t(%s)", "No hunting - reverting to no continuation");ECHO(echo_string, echo_file);
      return; 
    }

  hunt = (struct HC_Information *)array_alloc(1, nHC, sizeof(struct HC_Information));

  if ( Debug_Flag && ProcID == 0 )
    {
      printf("%s:\tallocated %d copies of the %lu sized\n", 
	     yo, nHC, (long unsigned int)sizeof(struct HC_Information));
      printf("%s:\tHC structure.\n", yo);
    }

  hunt[0].nHC = nHC;

  for (iHC=0;iHC<nHC;iHC++)
  {
      
      hunt[iHC].nHC = nHC;
	  
      look_for(ifp, "HC", input, '=');
	  
	  SPF(echo_string,"%s =", input);
	  
      if (fscanf(ifp, "%80s", input) != 1)
	  {
		  fprintf (stderr, "%s:\tError reading HC[iHC].Set_Type (one of BC, MT, UM, or AC)\n", yo);
		  exit (-1);
	  }
	  
      if (!strcmp(input,"BC"))
	  { 
		  hunt[iHC].Type = 1;
	  }
      else
		  if (!strcmp(input,"MT"))
		  { 
			  hunt[iHC].Type = 2;
		  }
      else
		  if (!strcmp(input,"AC"))
		  { 
			  hunt[iHC].Type = 3;
		  }
      else
		  if (!strcmp(input,"UM"))
		  { 
			  hunt[iHC].Type = 4;
		  }
      else
	  {
		  fprintf(stderr, "%s:\tImproper set_type for hunting condition - %s\n", yo, input);
		  exit (-1);
	  }
	  
	  SPF(endofstring(echo_string)," %s", input);
	  
      switch (hunt[iHC].Type) {
		  case 3:
		  case 1:
			  /* AC */
			  /* BC */
			  if (fscanf(ifp, "%d %d %d %le %le %le %le %le", 
						 &hunt[iHC].BCID, &hunt[iHC].DFID, &hunt[iHC].ramp, 
						 &hunt[iHC].BegParameterValue, &hunt[iHC].EndParameterValue,
						 &hunt[iHC].Delta_s0, 
						 &hunt[iHC].Delta_s_min, &hunt[iHC].Delta_s_max	) != 8)
			  {
				  fprintf(stderr,"%s:\tError reading hunt[%d].BCID hunt[%d].DFID\n", yo, iHC, iHC);
				  fprintf(stderr,"%s:\tRecall Format: HC = BCID DFID RAMP from to step step_min step_max\n",yo);
				  exit (-1);
			  }
			  SPF(endofstring(echo_string)," %d %d %d %.4g %.4g %.4g %.4g %.4g", 
				  hunt[iHC].BCID, hunt[iHC].DFID, hunt[iHC].ramp, 
				  hunt[iHC].BegParameterValue, hunt[iHC].EndParameterValue,
				  hunt[iHC].Delta_s0, hunt[iHC].Delta_s_min, hunt[iHC].Delta_s_max );
			  SPF(endofstring(echo_string), "\n\t( %3d. BC: BCID=%3d DFID=%5d)", iHC+1, hunt[iHC].BCID, hunt[iHC].DFID);
			  break;
			  
		  case 2: 
			  /* MT */
			  if (fscanf(ifp, "%d %d %d %le %le %le %le %le", 
						 &hunt[iHC].MTID, &hunt[iHC].MPID, &hunt[iHC].ramp, 
						 &hunt[iHC].BegParameterValue, &hunt[iHC].EndParameterValue,
						 &hunt[iHC].Delta_s0, 
						 &hunt[iHC].Delta_s_min, &hunt[iHC].Delta_s_max) != 8)
			  {
				  fprintf(stderr,"%s:\tError reading hunt[%d].MTID hunt[%d].MPID\n", yo, iHC, iHC);
				  fprintf(stderr,"%s:\tRecall Format: HC = MTID MPID RAMP from to step step_min step_max\n",yo);
				  exit (-1);
			  }
			  
			  SPF(endofstring(echo_string)," %d %d %d %.4g %.4g %.4g %.4g %.4g",
				  hunt[iHC].MTID, hunt[iHC].MPID, hunt[iHC].ramp, 
				  hunt[iHC].BegParameterValue, hunt[iHC].EndParameterValue,
				  hunt[iHC].Delta_s0, hunt[iHC].Delta_s_min, hunt[iHC].Delta_s_max );
			  SPF(endofstring(echo_string), "\n\t( %3d. MT: MTID=%3d MPID=%5d)", iHC+1, hunt[iHC].MTID, hunt[iHC].MPID);
			  hunt[iHC].MTID--;
			  break;
			  
		  case 4:
			  /* UM */
			  if (fscanf(ifp, "%d %d %d %d %le %le %le %le %le", 
						 &hunt[iHC].MTID, &hunt[iHC].MPID, &hunt[iHC].MDID, &hunt[iHC].ramp, 
						 &hunt[iHC].BegParameterValue, &hunt[iHC].EndParameterValue,
						 &hunt[iHC].Delta_s0, 
						 &hunt[iHC].Delta_s_min, &hunt[iHC].Delta_s_max) != 9)
			  {
				  fprintf(stderr,"%s:\tError reading hunt[%d].MTID hunt[%d].MPID hunt[%d].MDID\n", yo, iHC, iHC, iHC);
				  fprintf(stderr,"%s:\tRecall Format: HC = MTID MPID MDID RAMP from to step step_min step_max\n",yo);
				  exit (-1);
			  }
			  SPF(endofstring(echo_string),  " %d %d %d %d %.4g %.4g %.4g %.4g %.4g",
				  hunt[iHC].MTID, hunt[iHC].MPID, hunt[iHC].MDID, hunt[iHC].ramp, 
				  hunt[iHC].BegParameterValue, hunt[iHC].EndParameterValue,
				  hunt[iHC].Delta_s0,hunt[iHC].Delta_s_min, hunt[iHC].Delta_s_max);
			  SPF(endofstring(echo_string), "\n\t( %3d. MT: MTID=%3d MPID=%5d MDID=%3d)", iHC+1, hunt[iHC].MTID, hunt[iHC].MPID, hunt[iHC].MDID);
			  hunt[iHC].MTID--;
			  break;
			  
      }
	  
	  ECHO(echo_string,echo_file);
          if(hunt[iHC].ramp == 2 && 
             ((hunt[iHC].BegParameterValue == hunt[iHC].EndParameterValue) ||
             (hunt[iHC].BegParameterValue<0 || hunt[iHC].EndParameterValue<0)))
             {
              hunt[iHC].ramp = 0;  
              fprintf(stderr, "%s:\tImproper Log ramp for hunting condition %d\n", yo, iHC);
             }
  }

/* This section is required for backward compatibility - EDW */
  if (Continuation == LOCA)
    {

/* Allocate cpcc for specified number of HC's here.*/
      nCC = nHC;
      cpcc = (struct Continuation_Conditions *)
                  array_alloc(1, nCC, sizeof(struct Continuation_Conditions));
      if ( Debug_Flag && ProcID == 0 )
        {
          printf("%s:\tallocated %d new copies of the %lu sized\n", 
                 yo, nCC, (long unsigned int)sizeof(struct Continuation_Conditions));
          printf("%s:\tCC structure.\n", yo);
        }

/* Transfer needed inputs from hunt to cpcc */
      range_0 = hunt[0].EndParameterValue - hunt[0].BegParameterValue;
      if (hunt[0].ramp == 1)
        {
          hunt[0].Delta_s0 = range_0
            / ( (double)(cont->MaxPathSteps - 1) );
          loca_in->StepAggr = 0.0;
        }
      for (iHC=0; iHC<nHC; iHC++)
        {
          cpcc[iHC].nCC = nHC;
          cpcc[iHC].Type = hunt[iHC].Type;
          switch(hunt[iHC].Type) {
          case 3:  /* AC */
          case 1:  /* BC */
            cpcc[iHC].BCID = hunt[iHC].BCID;
            cpcc[iHC].DFID = hunt[iHC].DFID;
            break;
          case 2:  /* MT */
            cpcc[iHC].MTID = hunt[iHC].MTID;
            cpcc[iHC].MPID = hunt[iHC].MPID;
            break;
          case 4:  /* UM */
            cpcc[iHC].MTID = hunt[iHC].MTID;
            cpcc[iHC].MPID = hunt[iHC].MPID;
            cpcc[iHC].MDID = hunt[iHC].MDID;
            break;
          }
          cpcc[iHC].Beg_CC_Value = hunt[iHC].BegParameterValue;
          cpcc[iHC].End_CC_Value = hunt[iHC].EndParameterValue;
          range = hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue;
          cpcc[iHC].ratio = range / range_0;
        }

/* Fill in cont-> entries from hunt[0]. entries */
      hunt[0].Delta_s0 *= SGN(range_0);
      cont->upType = hunt[0].Type;
      cont->upBCID = hunt[0].BCID;
      cont->upDFID = hunt[0].DFID;
      cont->upMTID = hunt[0].MTID;
      cont->upMPID = hunt[0].MPID;
      cont->upMDID = hunt[0].MDID;
      cont->BegParameterValue = hunt[0].BegParameterValue;
      cont->EndParameterValue = hunt[0].EndParameterValue;
      cont->Delta_s0 = hunt[0].Delta_s0;
      cont->Delta_s_min = hunt[0].Delta_s_min;
      cont->Delta_s_max = hunt[0].Delta_s_max;
    }

}
/* rd_hunt_specs -- read input file for hunting specifications */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * rd_ac_specs -- read input file for augmenting conditions specifications
 *
 * Comments:    This code was lifted out of the ever-growing read_input_file()
 *              routine above. This makes things more modular, with division
 *              into "sections".
 *
 * Adapted by Ian Gates
 *
 * FORMAT:   
 * 
 * BOUNDARY CONDITION DF:  AC = BC BCID DFID
 *
 * MATERIAL PROPERTY:      AC = MT MTID MPID
 */

void 
rd_ac_specs(FILE *ifp,
	    char *input)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  int iAC;
  char *yo; 
  int num_const = 0,i;
  int AC_rd_file;
  double augc_initial_value[10];

  int mn, iread;
  int nAC1;
  int do_alc = FALSE;
  char string[MAX_FNL];

  double z[MAX_CONC];        /* species charge number       */
  const double F = 96487.0;  /* Faraday's constant in C/mol */
  struct Data_Table *AC_table=NULL;
  char  line[255];
  char  *arguments[MAX_NUMBER_PARAMS];
  char filename[80],interpolation[80],table_var[80];

  yo = "rd_ac_specs";

  /* initialize volume and LSvelocity integral cards */
  for(mn = 0; mn < MAX_NUMBER_MATLS; mn++)
    {
      pd_glob[mn]->VolumeIntegral = -1;
      pd_glob[mn]->LSVelocityIntegral = -1;
    }

  /* Read in Augmenting Conditions Specifications */

/*   iread = look_for_optional(ifp,"Augmenting Conditions Specifications",input,'='); */
/*   if (iread != 1) { */
/*     nAC = 0; */
/*     augc = NULL; */
/*     return; } */


  iread = look_for_optional(ifp,"Number of augmenting conditions",input,'=');
  if (iread != 1) { return; }

  if (fscanf (ifp,"%d",&nAC) != 1)
    {
      EH( -1, "error reading Number of augmenting conditions");
    }
  ECHO("\n***Augmenting Conditions***\n", echo_file);


  /*
   * Count augmenting conditions if nAC is set to -1
   */

  if (nAC == -1)
    {
      nAC = count_list(ifp, "AC", input, '=', "END OF AC");
    }

  /*
   * Add one more AC if doing arc length continuation
   */

  nAC1 = nAC;
  if (nAC > 0 && loca_in->Cont_Order == 2)
    {
      do_alc = TRUE;
      nAC++;
    }

  SPF(echo_string,"%s = %d", "Number of augmenting conditions", nAC); ECHO(echo_string,echo_file);

  /*
   *  Allocate space for the vector augc
   */

  if (nAC == 0) 
    {
      return; 
    }

  AC_rd_file = 0;
  iread = look_for_optional(ifp,"Augmenting Conditions Initial Guess",input,'=');
  if(iread == 1)
    {
    if ( fscanf( ifp,"%s",string) != 1 )
      {
    fprintf(stderr,"%s:\tError reading Initial Guess string %s\n", yo, input);
      }
    if ( strcmp(string,"read") == 0 )
      {
        AC_rd_file = 1;
      }
    else if( strcmp(string,"initialize") == 0)
        {
        AC_rd_file = 2;
      for (iAC=0;iAC<nAC1;iAC++)
              {
                if( fscanf(ifp,"%lf",&augc_initial_value[iAC]) != 1 )
                      {
          fprintf(stderr,"%s:\tError reading augc_initial_value[%d]\n", yo, iAC);
          fprintf(stderr,"\tAdd AC value after initialize?\n");
                      }
              }
        }
     else
      {
          fprintf(stderr,"%s:\tread or initialize  - not %s\n", yo, string);
      }

    SPF(echo_string,eoformat, "Augmenting Conditions Initial Guess", input); ECHO(echo_string,echo_file);
    }


  augc = alloc_struct_1(struct AC_Information, nAC );

  if ( Debug_Flag && ProcID == 0 )
    {
      printf("%s:\tallocated %d copies of the %lu sized\n", 
	     yo, nAC, (long unsigned int)sizeof(struct AC_Information));
      printf("%s:\tAC structure.\n", yo);
    }

  augc[0].nAC = nAC;

  for (iAC=0;iAC<nAC1;iAC++)    /* Begin loop over ACs */
    {
      
      augc[iAC].nAC = nAC;
      augc[iAC].BCID = -99;
      augc[iAC].MTID = -99;
      augc[iAC].VOLID = -99;
      augc[iAC].MFID = -99;
      augc[iAC].iread=0;
      if(AC_rd_file == 1) {augc[iAC].iread=1;}
      if(AC_rd_file == 2) {augc[iAC].iread=2;}
      augc[iAC].tmp3 = augc_initial_value[iAC];

      look_for(ifp, "AC", input, '=');

      if (fscanf(ifp, "%80s", input) != 1)
	{
	  fprintf (stderr, "%s:\tError reading AC[iAC].Set_Type (one of BC, MT, VC, LSV, FC, LM, or OV)\n", yo);
	  exit (-1);
	}

      stringup( input );

      if (!strcmp(input,"BC") || !strcmp(input,"USERBC") )
        { 
	  augc[iAC].Type = AC_USERBC;
	}
      else
      if (!strcmp(input,"MT") || !strcmp(input,"USERMAT"))
        { 
	  augc[iAC].Type = AC_USERMAT;
	}
      else
      if (!strcmp(input,"VC") || !strcmp(input,"VOLUME") )
        { 
	  augc[iAC].Type = AC_VOLUME;
	}
      else
      if (!strcmp(input,"XY") || !strcmp(input,"POSITION") )
        { 
	  augc[iAC].Type = AC_POSITION;
	}
      else
      if (!strcmp(input,"LSV") || !strcmp(input,"LS_VELOCITY") )
        { 
	  augc[iAC].Type = AC_LS_VEL;

	}
      else
      if (!strcmp(input,"FC")|| !strcmp(input,"FLUX"))
        { 
	  augc[iAC].Type = AC_FLUX;
	}
      else
      if( !strcmp(input,"LM") || !strcmp(input, "LGR") )
	{
	  augc[iAC].Type = AC_LGRM;
	}
      else
      if( !strcmp(input,"OV") || !strcmp(input, "OVERLAP") )
        {
          augc[iAC].Type = AC_OVERLAP;
        }
      else
      if( !strcmp(input,"PBC") || !strcmp(input, "PERIODIC") )
        {
          augc[iAC].Type = AC_PERIODIC;
        }
      else
      if (!strcmp(input,"FCM")|| !strcmp(input,"FLUXMAT"))
        {
        augc[iAC].Type = AC_FLUX_MAT;
      }
      else
	{
	  fprintf(stderr, "%s:\tImproper set_type for augmenting condition - %s\n", yo, input);
	  exit (-1);
	}

      SPF(echo_string,eoformat, "AC", input);

      augc[iAC].len_AC = 0;
      augc[iAC].LewisNum = 0.0;

      switch (augc[iAC].Type) {

      case AC_USERBC: /* BC */
	if (fscanf(ifp, "%d %d", &augc[iAC].BCID, &augc[iAC].DFID) != 2)
	  {
	    fprintf(stderr,"%s:\tError reading augc[%d].BCID augc[%d].DFID\n", yo, iAC, iAC);
	    fprintf(stderr,"%s:\tRecall Format: AC = BCID DFID\n",yo);
	    exit (-1);
	  }
	SPF(endofstring(echo_string)," %d %d", augc[iAC].BCID, augc[iAC].DFID );
  
  	if( Num_Proc != 1 && augc[iAC].BCID == APREPRO_AC_BCID)
	  {
	    EH(-1,"Aprepro parameter AC not ready for parallel");
	  } 
	break;

      case AC_USERMAT: /* MT */
	if (fscanf(ifp, "%d %d", &augc[iAC].MTID, &augc[iAC].MPID) != 2)
	  {
	    fprintf(stderr,"%s:\tError reading augc[%d].MTID augc[%d].MPID\n", yo, iAC, iAC);
	    fprintf(stderr,"%s:\tRecall Format: AC = MTID MPID\n",yo);
	    exit (-1);
	  }
	SPF(endofstring(echo_string)," %d %d", augc[iAC].MTID, augc[iAC].MPID );

	break;

      case AC_VOLUME: /* VC */
	if (fscanf(ifp, "%d %d %d %d %d %lf", &augc[iAC].MTID, &augc[iAC].VOLID, &augc[iAC].BCID,
		   &augc[iAC].DFID, &augc[iAC].COMPID, &augc[iAC].CONSTV) != 6)
	  {
	    fprintf(stderr, "%s:\tError reading MTID, VOLID, BCID, DFID, COMPID, CONSTV\n", yo);
	    fprintf(stderr, "%s:\tRecall Format:AC=MTID VOLID BCID DFID COMPID CONSTV\n", yo);
	    exit (-1);
	  }

	SPF(endofstring(echo_string)," %d %d %d %d %d %.4g", augc[iAC].MTID, 
	    augc[iAC].VOLID, augc[iAC].BCID,
	    augc[iAC].DFID, augc[iAC].COMPID, augc[iAC].CONSTV);

	break;

      case AC_POSITION:
	/*
	 *  first  int - NSID Node set ID to be used
	 *  second int - CoordID direction to be used
	 *  third  int - BCID BC index to use for variable
	 *  4th    int - card float index to use for variable
	 *  5th    int - Form of the 1D position parameterization (0)
	 *  6th   flt  - Value of the coordinate position
	 */
        read_line(ifp, input, FALSE);
        augc[iAC].LewisNum = 0.0;
	if (sscanf(input, "%d %d %d %d %d %lf %lf", &augc[iAC].MTID, &augc[iAC].VOLID, &augc[iAC].BCID,
		   &augc[iAC].DFID, &augc[iAC].COMPID, &augc[iAC].CONSTV,
		   &augc[iAC].LewisNum) < 6)
	  {
	    fprintf(stderr,"%s:\tError reading NSID, CoordID, BCID, DFID, FORMID, CONSTV\n", yo);
	    fprintf(stderr,"%s:\tRecall Format:AC=NSID CoordID BCID DFID FORMID CONSTV\n",yo);
	    exit (-1);
	  }
	
	SPF(endofstring(echo_string)," %d %d %d %d %d %.4g", augc[iAC].MTID, 
	    augc[iAC].VOLID, augc[iAC].BCID,
	    augc[iAC].DFID, augc[iAC].COMPID, augc[iAC].CONSTV);
	
	break;  

      case AC_LS_VEL: /* LSV */
	if (fscanf(ifp, "%d %d %d %d %s", &augc[iAC].MTID, &augc[iAC].BCID,
		   &augc[iAC].DFID, &augc[iAC].LSPHASE, input) != 5)
	  {
	    fprintf(stderr,"%s:\tError reading MTID, BCID, DFID, PHASE, DIR \n", yo);
	    fprintf(stderr,"%s:\tRecall Format:AC=LSV MTID BCID DFID PHASE DIR \n",yo);
	    exit (-1);
	  }

	  if( strcmp( input, "VX" ) == 0 )
	    {
	      if( augc[iAC].LSPHASE > 0 ) { 
		augc[iAC].LSPHASE = I_POS_FILL;
		augc[iAC].DIR = I_POS_VX;
	      }
              else if( augc[iAC].LSPHASE < 0 ) {
		augc[iAC].LSPHASE=I_NEG_FILL;
		augc[iAC].DIR = I_NEG_VX;
	      }
	      else {
		EH(-1,"Error in reading LSV:PHASE; it can not equal zero.");
	      }

	      if( fscanf(ifp,"%lf", &augc[iAC].CONSTV) != 1 )
		{
		  EH(-1,"Error in reading LS_VEL - VX augmenting condition");
		}
	    }
	  else if ( strcmp( input, "VY" ) == 0 )
	    {
	      if( augc[iAC].LSPHASE > 0 ) { 
		augc[iAC].LSPHASE = I_POS_FILL;
		augc[iAC].DIR = I_POS_VY;
	      }
              else if( augc[iAC].LSPHASE < 0 ) {
		augc[iAC].LSPHASE = I_NEG_FILL;
		augc[iAC].DIR = I_NEG_VY;
	      }
	      else {
		EH(-1,"Error in reading LSV:PHASE; it can not equal zero.");
	      }

	      if( fscanf(ifp,"%lf", &augc[iAC].CONSTV) != 1 )
		{
		  EH(-1,"Error in reading LS_VEL - VY augmenting condition");
		}
	    }
	  else if ( strcmp( input, "VZ" ) == 0 )
	    {
	      if( augc[iAC].LSPHASE > 0 ) { 
		augc[iAC].LSPHASE = I_POS_FILL;
		augc[iAC].DIR = I_POS_VZ;
	      }
              else if( augc[iAC].LSPHASE < 0 ) {
		augc[iAC].LSPHASE = I_NEG_FILL;
		augc[iAC].DIR = I_NEG_VZ;
	      }
	      else {
		EH(-1,"Error in reading LSV:PHASE; it can not equal zero.");
	      }

	      if( fscanf(ifp,"%lf", &augc[iAC].CONSTV) != 1 )
		{
		  EH(-1,"Error in reading LS_VEL - VZ augmenting condition");
		}
	    }
 	  else
	    {
	      EH(-1, "Must specify the velocity component for LS_VEL augmenting condition as VX, VY or VZ.");
	    }


	  SPF(endofstring(echo_string)," %d %d %d %d %s %.4g", augc[iAC].MTID, augc[iAC].BCID,
	      augc[iAC].DFID, augc[iAC].LSPHASE, input, augc[iAC].CONSTV );

	break;  

      case AC_FLUX: /* FC */
      case AC_FLUX_MAT: /* FC */
	    /* BCID - Index of boundary condition controlling AC
	     * DFID - Index of floating point parameter to serve as additional unknown
	     * input - string that indicates the type of flux constraint.
	     * SSID - This is the sideset number where the flux is to be computed
	     * CONSTV - This is the value the flux should take over that side side
	     */

        if( augc[iAC].Type == AC_FLUX)
          {
	  if( fscanf( ifp, "%d %d %d %s", &augc[iAC].MTID, &augc[iAC].BCID, &augc[iAC].DFID, input) != 4 )
	    {
	      EH(-1, "Error in syntax for flux-type augmenting condition");
	    }

	  SPF(endofstring(echo_string)," %d %d %d %s", augc[iAC].MTID, augc[iAC].BCID, augc[iAC].DFID, input );

     	  if( Num_Proc != 1 && augc[iAC].BCID == APREPRO_AC_BCID)
  		{
  		EH(-1,"Aprepro parameter AC not ready for parallel");
  		} 
          }
        else
          {
        if( fscanf( ifp, "%d %d %d %s", &augc[iAC].MTID, &augc[iAC].MPID, &augc[iAC].MDID, input) != 4 )
           {
            fprintf(stderr,"%s:\tError reading augc[%d].MTID augc[%d].MPID augc[%d].MDID\n", yo, iAC, iAC, iAC);
            exit (-1);
           }
          SPF(endofstring(echo_string)," %d %d %d %s", augc[iAC].MTID, augc[iAC].MPID, augc[iAC].MDID, input );
          }

	  if( strcmp( input, "FORCE_X" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_X ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_X augmenting condition");
		}
	    }
          else if( strcmp( input, "FORCE_X_POS" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_X_POS ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_X_POS augmenting condition");
		}
	    }
	  else if( strcmp( input, "FORCE_X_NEG" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_X_NEG ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_X_NEG augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Y" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Y ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Y augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Y_POS" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Y_POS ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Y_POS augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Y_NEG" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Y_NEG ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Y_NEG augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Z" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Z ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Z augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Z_POS" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Z_POS ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Z_POS augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_Z_NEG" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_Z_NEG ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_Z_NEG augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_NORMAL" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_NORMAL ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_NORMAL augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_TANGENT1" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_TANGENT1 ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_TANGENT1 augmenting condition");
		}
	    }
	  else if ( strcmp( input, "FORCE_TANGENT2" ) == 0 )
	    {
	      augc[iAC].MFID = FORCE_TANGENT2 ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading FORCE_TANGENT2 augmenting condition");
		}
	    }
	  else if ( strcmp( input, "VOLUME_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = VOLUME_FLUX ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading VOLUME_FLUX augmenting condition");
		}
	    }
	  else if ( strcmp( input, "SHELL_VOLUME_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = SHELL_VOLUME_FLUX ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading SHELL_VOLUME_FLUX augmenting condition");
		}
  	     }
	  else if ( strcmp( input, "HEAT_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = HEAT_FLUX ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading HEAT_FLUX augmenting condition");
		}
	    }
	  else if ( strcmp( input, "SPECIES_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = SPECIES_FLUX ;

	      if( fscanf(ifp,"%d %d %lf", &augc[iAC].SSID, &augc[iAC].COMPID,&augc[iAC].CONSTV) != 3 )
		{
		  EH(-1,"Error in reading SPECIES_FLUX augmenting condition");
		}
	    }
	  else if ( strcmp( input, "CHARGED_SPECIES_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = CHARGED_SPECIES_FLUX ;

	      if( fscanf(ifp,"%d %d %lf", &augc[iAC].SSID, &augc[iAC].COMPID, &augc[iAC].CONSTV) != 3 )
		{
		  EH(-1,"Error in reading CHARGED_SPECIES_FLUX augmenting condition");
		}
	    }
	  else if ( strcmp( input, "CURRENT" ) == 0 )
	    {
	      augc[iAC].MFID = CURRENT ;

	      if( fscanf(ifp,"%d %d %lf", &augc[iAC].SSID, &augc[iAC].COMPID, &augc[iAC].CONSTV) != 3) 
		{
		  EH(-1,"Error in reading CURRENT augmenting condition");
		}
	    }
	  else if ( strcmp( input, "CURRENT_FICKIAN" ) == 0 )
	    {
	      augc[iAC].MFID = CURRENT_FICKIAN ;

	      if( fscanf(ifp,"%d %d %lf %lf",&augc[iAC].SSID, &augc[iAC].COMPID, &augc[iAC].CONSTV, 
                                             &z[augc[iAC].COMPID]) != 4 )
		{
		  EH(-1,"Error in reading CURRENT_FICKIAN augmenting condition");
		}
              augc[iAC].CONSTV /= z[augc[iAC].COMPID]*F; /* convert current density to molar flux */ 
	    }
	  else if ( strcmp( input, "PORE_LIQ_FLUX" ) == 0 )
	    {
	      augc[iAC].MFID = PORE_LIQ_FLUX ;

	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
		{
		  EH(-1,"Error in reading PORE_LIQ_FLUX augmenting condition");
		}
	    }
 	  else if ( strcmp( input, "TORQUE" ) == 0 )
 	    {
 	      augc[iAC].MFID = TORQUE ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading TORQUE augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "AVERAGE_CONC" ) == 0 )
 	    {
 	      augc[iAC].MFID = AVERAGE_CONC ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading AVERAGE_CONC augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "SURF_DISSIP" ) == 0 )
 	    {
 	      augc[iAC].MFID = SURF_DISSIP ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading SURF_DISSIP augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "AREA" ) == 0 )
 	    {
 	      augc[iAC].MFID = AREA ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading AREA augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "VOL_REVOLUTION" ) == 0 )
 	    {
 	      augc[iAC].MFID = VOL_REVOLUTION ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading VOL_REVOLUTION augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "NEG_LS_FLUX" ) == 0 )
 	    {
 	      augc[iAC].MFID = NEG_LS_FLUX ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading NEG_LS_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "POS_LS_FLUX" ) == 0 )
 	    {
 	      augc[iAC].MFID = POS_LS_FLUX ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading POS_LS_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "N_DOT_X" ) == 0 )
 	    {
 	      augc[iAC].MFID = N_DOT_X ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading N_DOT_X augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_NORMAL" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_NORMAL ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_NORMAL augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_TANGENT1" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_TANGENT1 ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_TANGENT1 augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_TANGENT2" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_TANGENT2 ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_TANGENT2 augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_X" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_X ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_X augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_Y" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_Y ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_Y augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ELEC_FORCE_Z" ) == 0 )
 	    {
 	      augc[iAC].MFID = ELEC_FORCE_Z ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ELEC_FORCE_Z augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "NET_CHARGE" ) == 0 )
 	    {
 	      augc[iAC].MFID = NET_SURF_CHARGE ;
 
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading NET_CHARGE augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_NORMAL" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_NORMAL ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_TANGENT1" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_TANGENT1 ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_TANGENT2" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_TANGENT2 ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_X" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_X ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_Y" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_Y ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "ACOUSTIC_FLUX_Z" ) == 0 )
 	    {
 	      augc[iAC].MFID = ACOUSTIC_FLUX_Z ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading ACOUSTIC_FLUX augmenting condition");
 		}
 	    }
 	  else if ( strcmp( input, "REPULSIVE_FORCE" ) == 0 )
 	    {
 	      augc[iAC].MFID = REPULSIVE_FORCE ;
 	      if( fscanf(ifp,"%d %lf",&augc[iAC].SSID, &augc[iAC].CONSTV) != 2 )
 		{
 		  EH(-1,"Error in reading REPULSIVE_FORCE augmenting condition");
 		}
 	    }
	  else if ( strcmp( input, "SPECIES_FLUX_REVOLUTION" ) == 0 )
	    {
	      augc[iAC].MFID = SPECIES_FLUX_REVOLUTION ;

	      if( fscanf(ifp,"%d %d %lf", &augc[iAC].SSID, &augc[iAC].COMPID,&augc[iAC].CONSTV) != 3 )
		{
		  EH(-1,"Error in reading SPECIES_FLUX_REVOLUTION augmenting condition");
		}
	    }
 	  else
	    {
	      EH(-1, "That flux-type augmenting condition has yet to be implemented");
	    }

	  switch( augc[iAC].MFID )
	    {
	    case CHARGED_SPECIES_FLUX:
	    case CURRENT:
	    case SPECIES_FLUX:
	    case SPECIES_FLUX_REVOLUTION:
	      SPF(endofstring(echo_string)," %d %d %.4g",  augc[iAC].SSID, augc[iAC].COMPID, augc[iAC].CONSTV );
	      break;
	    case CURRENT_FICKIAN:
	      SPF(endofstring(echo_string)," %d %d %.4g %.4g", augc[iAC].SSID, augc[iAC].COMPID, augc[iAC].CONSTV, z[augc[iAC].COMPID]);
	      break;
	    default:
	      SPF(endofstring(echo_string)," %d %.4g", augc[iAC].SSID, augc[iAC].CONSTV );
	      break;
	    }
	  break;

      case AC_LGRM:
	    /* MTID - Block identification 
	     * SSID - This is the sideset number over which the constraint is to be integrated
	     * CONSTV - This is the value the flux should take over that side side
	     */

	  if( fscanf( ifp, "%d %d %s %lf", &augc[iAC].MTID, &augc[iAC].SSID,  input, &augc[iAC].CONSTV ) != 4 )
	    {
	      EH(-1, "Error in syntax for Lagrange Multiplier augmenting condition");
	    }

	  stringup( input );

	  if ( !strcmp(input,"VOLUME_FLOW") || !strcmp(input,"VOLUME_FLUX") )
	    {
	      augc[iAC].MFID = VOLUME_FLUX;
	      augc[iAC].COMPID = 0.0;
	    }
 	  else
	    {
	      EH(-1,"Unknown Lagrange Multiplier constraint type.");
	    }

	  SPF(endofstring(echo_string)," %d %d %s %.4g",augc[iAC].MTID, augc[iAC].SSID,  input, augc[iAC].CONSTV );

	  break;

  /*
   * For overlapping grid problems, one OVERLAP AC card is provided.
   * The correct number of AC constraints will be created before
   * solving the problem, and augc[] will be reallocated accordingly.
   * This will take care of all embedded and contact BC types.
   */
      case AC_OVERLAP:
            /* SSID - Number of side set for moving mesh surface.
             * solid_eb - Solid material element block.
             * fluid_eb - Solid material element block.
             * lm_eb - Element block on which the Lagrange multiplier
             *         equation to be used for the constraints is defined
             *         (Note that there may be more than one LM equation).
             * --- THE ABOVE FOUR ENTRIES COME FROM THE EXODUSII FEM FILE ---
             *
             *  NOTE: "cross-mesh terms" are sensitivities of equations
             *        defined on one block to unknowns defined on the other.
             *        These terms have no space allocated in the Jacobian.
             */
        if (fscanf(ifp, "%d %d %d %d", &augc[iAC].SSID, &augc[iAC].solid_eb,
             &augc[iAC].fluid_eb, &augc[iAC].lm_eb) < 4)
           {
             EH(-1, "Error reading SSID, solid_eb, fluid_eb, lm_eb\n");
           }

  /* Do some error checking */
        if (augc[iAC].solid_eb == augc[iAC].fluid_eb)
          {
            EH(-1, "Fluid and solid element blocks cannot be the same!");
          }
        if (augc[iAC].solid_eb != augc[iAC].lm_eb &&
            augc[iAC].fluid_eb != augc[iAC].lm_eb)
          {
            EH(-1, "LM element block must be same as either fluid and solid!");
          }

  /* If these tests were passed, set global flag */
        Do_Overlap = TRUE;

	SPF(endofstring(echo_string)," %d %d %d %d", augc[iAC].SSID, augc[iAC].solid_eb, augc[iAC].fluid_eb, augc[iAC].lm_eb );

        break;
        
  /*
   * For periodic bc's, one PERIODIC AC card is provided per variable.
   * The correct number of AC constraints will be created before
   * solving the problem, and augc[] will be reallocated accordingly.
   */
      case AC_PERIODIC:
            /* SSID  - Number of side set of surface 1 to apply periodic bc on
	     * SSID2 - Number of side set of surface 2 to apply periodic bc on
             * --- THE ABOVE FOUR ENTRIES COME FROM THE EXODUSII FEM FILE ---
             *
             *  VAR - Variable to be made periodic
             */
	if( fscanf( ifp, "%s %d %d", input, &augc[iAC].SSID, &augc[iAC].SSID2 ) != 3 )
	    {
	      EH(-1, "Error in syntax for Periodic Boundary Condition augmenting condition");
	    }
	else
          {
	    int var, var_found = FALSE;
	    /* loop through equation names and compare to input */
	    stringup(input);
            for (var=0; var<MAX_VARIABLE_TYPES && !var_found; var++)
	      {
	        if ( !strcmp(input, Var_Name[var].name1) ||
	             !strcmp(input, Var_Name[var].name2)) var_found = TRUE;
	      } 
            if (!var_found)
	      {
	        EH(-1, "Could not identify variable type for Periodic Boundary Condition");
	      }
	    augc[iAC].VAR = var-1; /* this got incremented one more time than desired */
	  }

	SPF(endofstring(echo_string)," %s %d %d", input, augc[iAC].SSID, augc[iAC].SSID2);

        break;

	  
      default:
	EH(-1,"Unknown Augmenting Condition type.");

      }  /* End switch on AC_Type */

  /*
   * Fill in type for arc length AC
   */

  if (do_alc)
    {
      augc[nAC1].Type = AC_ARC_LENGTH;
      augc[nAC1].DFID = -1;
      augc[nAC1].CONSTV = cont->BegParameterValue;
    }

/*  read parameter filename for Aprepro parameter cases  */

  if( (augc[iAC].Type == AC_USERBC ||
       augc[iAC].Type == AC_FLUX ||
       augc[iAC].Type == AC_FLUX_MAT) &&
       augc[iAC].BCID == APREPRO_AC_BCID )
    {
      if (fscanf (ifp,"%s",string) != 1)
	{
	  EH( -1, "error reading Parameter File name");
	}
      strcpy(augc[iAC].Params_File, string);
      if (fscanf (ifp,"%s",string) != 1)
	{
	  EH( -1, "error reading Parameter name");
	}
      strcpy(augc[iAC].AP_param, string);


      SPF(endofstring(echo_string)," %s %s", augc[iAC].Params_File, augc[iAC].AP_param);
    }

/* add float list */

  if ( fgets(line, 255, ifp) != NULL)
    {
      strip(line);
      if ((num_const = count_parameters(line)) > 0) {
        iread = tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);	
  	augc[iAC].len_AC = 0;
	augc[iAC].DataFlt = alloc_dbl_1(num_const, 0.0);
	for(i = 0; i < num_const; i++) {
		if(!strncmp(arguments[i],"AC_TABLE",8))
		{
		if((i+3) < num_const)
		   {
      		strcpy(filename, arguments[i+1]);
      		strcpy(table_var, arguments[i+2]);
      		strcpy(interpolation, arguments[i+3]);
	  AC_Tables[0] = setup_table_AC( filename, AC_table, table_var, interpolation);
	  num_AC_Tables++;
		break;
		   }
		}	else	{
	  augc[iAC].DataFlt[i] = atof(arguments[i]);
  	  augc[iAC].len_AC++;
          SPF( endofstring(echo_string)," %.4g", augc[iAC].DataFlt[i]);
		}
	}
      }
    }
  if ( num_const < 0 )
    {
      sr = sprintf(err_msg, "? AC data floats\n");
      EH(-1, err_msg);
    }


  ECHO(echo_string,echo_file);
  
    }  /* End of iAC loop */

  ECHO("END OF AC", echo_file);
}
/* rd_ac_specs -- read input file for augmenting conditions specifications */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * rd_solver_specs -- read input file for equation solver specifications
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Krysolve has been ripped out by the roots. What remains
 *		is Aztec (krysolve 2.0) and a route to Kundert's sparse
 *		matrix direct factorization routines.
 *
 * Revised:			1996/07/09 07:33 MDT pasacki@sandia.gov
 */

void 
rd_solver_specs(FILE *ifp,
		char *input )
{
  char *c;
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  char def_form[MAX_CHAR_IN_INPUT]= " (%s = %s) %s";

  int iread, k, i;
  int is_Solver_Serial = TRUE;
  
  /*
   * Caution! This routine is gradually evolving to a dumb reader. Checking
   * for consistent and meaningful input is now done mostly in
   * set_aztec_options_params() in mm_sol_nonlinear.c
   *
   */

  /*
   * Default -- use Aztec, except for the linear solver.
   *            Insure that reasonable defaults are
   *		available in the "Matrix_" string variables. Otherwise,
   *		the parsing in set_aztec_options_params() in mm_sol_nonlinear.c
   *		will throw us out.
   *
   */

  strcpy(Matrix_Solver,			"lu"); /* Kundert's - not y12m ! */
  strcpy(Matrix_Format,			"msr"); 
  strcpy(Matrix_Scaling,		"none");
  strcpy(Matrix_Preconditioner,		"none");
  strcpy(Matrix_Subdomain_Solver,	"ilut"); /* Aztec 2 */
  strcpy(Matrix_Residual_Norm_Type,	"r0");
  strcpy(Matrix_Output_Type,		"none");
  strcpy(Matrix_Factorization_Reuse,	"recalc");
  strcpy(Matrix_Graph_Fillin,		"0"); /* Aztec 2 */
  strcpy(Matrix_Maximum_Iterations,	"500");
  strcpy(Matrix_Polynomial_Order,	"3");
  strcpy(Matrix_Factor_Overlap,		"none");
  strcpy(Matrix_Overlap_Type,		"standard");
  strcpy(Matrix_Krylov_Subspace,	"30");
  strcpy(Matrix_Orthogonalization,	"classic");
  strcpy(Matrix_Auxiliary_Vector,	"resid");
  strcpy(Matrix_Convergence_Tolerance,	"1e-6");
  strcpy(Matrix_Drop_Tolerance,		"0");
  strcpy(Matrix_Factorization_Save,	"0");
  strcpy(Matrix_ILUT_Fill_Factor,	"1");
  strcpy(Matrix_RILU_Relax_Factor,	"1");
  strcpy(Matrix_BILU_Threshold,         "0");
  strcpy(Matrix_Relative_Threshold,     "0");
  strcpy(Matrix_Absolute_Threshold,     "0");
  strcpy(Matrix_Reorder,                "none");
  strcpy(Amesos_Package,			    "KLU");

  /*  Read in Solver specifications */
  
  iread=look_for_optional(ifp,"Solver Specifications",input,'=');

  ECHO("\n***Solver Specifications***\n", echo_file);

  strcpy(search_string, "Solution Algorithm");

  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Solver, input); /* save string for aztec use */
      SPF(echo_string, eoformat, search_string, Matrix_Solver);
    }  
  else
    {    
      SPF(echo_string, eoformat, search_string, Matrix_Solver);
      strcat(echo_string, default_string);
    }

  ECHO(echo_string,echo_file);

  /*

   Which linear equation solver package will be used? 

   Most options for Matrix_Solver refer to Aztec names, 
   with 5 exceptions:

   "lu" means to use Kundert's Sparse matrix factorization package
   "ma28" means to use the old Harwell routines
   "umf" means that Davis' UMFPACK is used  
   "umff" means that Davis' UMFPACK is used 
               -- force full analysis/factorization every Newton step
   "front" means Hood's element-by-element frontal solver

   Note: to use the y12m direct solver in Aztec, set the solver to "y12m".

   */

  if (  strcmp(Matrix_Solver, "umf") == 0 )
    {
      Linear_Solver = UMFPACK2;
    }
  else 
  if (  strcmp(Matrix_Solver, "umff") == 0 )
    {
      Linear_Solver = UMFPACK2F;
   }
  else 
  if (  strcmp(Matrix_Solver, "lu") == 0 )
    {
      Linear_Solver = SPARSE13a;
    }
  else 
  if (  strcmp(Matrix_Solver, "ma28") == 0 )
    {
      Linear_Solver = MA28;
    }
  else
  if (  strcmp(Matrix_Solver, "front") == 0 )
    {
      i = sizeof(c);
      if (i == 4)
	{
	  Linear_Solver = FRONT;
	}
      else if (i == 8)
	{
	  EH(-1, "Benner's frontal solver not functional on 64-bit builds. We suggest using umf or at your own risk remove this error-handler statement");
	}
      else
	{
	  EH(-1,"frontal solver not available on this architecture. Please use umf");
	}
    }
  else
  if (strcmp(Matrix_Solver, "amesos") == 0) {
    Linear_Solver = AMESOS;
    is_Solver_Serial = FALSE;
  } else if (strcmp(Matrix_Solver, "aztecoo") == 0) {
    Linear_Solver = AZTECOO;
    is_Solver_Serial = FALSE;
  } else if (strcmp(Matrix_Solver, "stratimikos") == 0) {
    Linear_Solver = STRATIMIKOS;
    is_Solver_Serial = FALSE;
  } else {
    Linear_Solver = AZTEC;
    is_Solver_Serial = FALSE;
  }
  

  if ( strcmp(Matrix_Solver, "front") != 0 )
    {
      /*  Read in Matrix Format only if we are going to assemble one */
      
      strcpy(search_string, "Matrix storage format");
      iread=look_for_optional(ifp, search_string, input, '=');
      if ( iread == 1 ) 
	{
	  read_string(ifp,input,'\n');
	  strip(input);
	  strcpy(Matrix_Format, input); /* save string for aztec use */
	  SPF(echo_string,eoformat, search_string, Matrix_Format );
	}      
      else
	{
	  SPF(echo_string,def_form, search_string, Matrix_Format, default_string );
	}
      ECHO(echo_string,echo_file);
    }
  else
    {
      strcpy( Matrix_Format,"front");
    }


  /*

  OPTIONAL UMFPACK_IDIM SPEC

  */

  UMFPACK_IDIM = -1;
      
  strcpy(search_string, "UMF IDIM");
  iread = look_for_optional(ifp, search_string, input, '=');
  if (iread == 1) 
    { 
      if (fscanf(ifp, "%d", &UMFPACK_IDIM) != 1)
	{
	  iread = -1; 
	  SPF(echo_string," (%s = %d) %s", search_string, UMFPACK_IDIM, default_string); ECHO(echo_string,echo_file);
	}
      else
	SPF(echo_string,"%s = %d", search_string, UMFPACK_IDIM); ECHO(echo_string,echo_file);
    }



  /*

  OPTIONAL UMFPACK_XDIM SPEC

  */
  
  UMFPACK_XDIM = -1;

  iread = look_for_optional(ifp, "UMF XDIM", input, '=');
  if (iread == 1) 
    { 
      if (fscanf(ifp, "%d", &UMFPACK_XDIM) != 1)
	{ 
	  iread = -1; 
	  SPF(echo_string," (%s = %d) %s","UMF_XDIM",UMFPACK_XDIM, default_string ); ECHO(echo_string,echo_file);
	}
      else
	SPF(echo_string,"%s = %d","UMF_XDIM",UMFPACK_XDIM ); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "AztecOO Solver");
  iread = look_for_optional(ifp, search_string, input, '=');
  if (iread == 1) {
    read_string(ifp, input, '\n');
    strip(input);
    stringup(input);
    strcpy(AztecOO_Solver, input);
    SPF(echo_string, eoformat, search_string, input);
    ECHO(echo_string, echo_file);
  } else {
    // Set gmres as the default AztecOO Solver
    SPF(echo_string, def_form, search_string, "gmres", default_string);
    strcpy(AztecOO_Solver, "gmres");
    ECHO(echo_string, echo_file);
  }

  strcpy(search_string, "Stratimikos File");
  iread = look_for_optional(ifp, search_string, input, '=');
  if (iread == 1) {
    read_string(ifp, input, '\n');
    strip(input);
    strcpy(Stratimikos_File, input);
    SPF(echo_string, eoformat, search_string, input);
    ECHO(echo_string, echo_file);
  } else {
    // Set stratimikos.xml as defualt stratimikos file
    SPF(echo_string, def_form, search_string, "stratimikos.xml", default_string);
    strcpy(Stratimikos_File, "stratimikos.xml");
    ECHO(echo_string, echo_file);
  }

  strcpy(search_string, "Preconditioner");

  iread = look_for_optional(ifp, search_string, input, '=');

  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Preconditioner, input); /* save string for aztec use */

      SPF(echo_string,eoformat, search_string, Matrix_Preconditioner); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Preconditioner, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix subdomain solver");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Subdomain_Solver, input); /* save string for aztec use */
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Subdomain_Solver , default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix Scaling");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Scaling, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Scaling , default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix residual norm type");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Residual_Norm_Type, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Residual_Norm_Type, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix output type");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Output_Type, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Output_Type , default_string); ECHO(echo_string,echo_file);
     }

  strcpy(search_string, "Matrix factorization reuse");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Factorization_Reuse, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Factorization_Reuse , default_string); ECHO(echo_string,echo_file);
    }
  
  strcpy(search_string, "Matrix graph fillin");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Graph_Fillin, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Graph_Fillin, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix factorization overlap");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Factor_Overlap, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Factor_Overlap, default_string); ECHO(echo_string,echo_file);
    }
	
  strcpy(search_string, "Matrix overlap type");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Overlap_Type, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Overlap_Type, default_string); ECHO(echo_string,echo_file);
    }
	
  strcpy(search_string, "Matrix auxiliary vector");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Auxiliary_Vector, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Auxiliary_Vector, default_string); ECHO(echo_string,echo_file);
    }
	
  strcpy(search_string, "Matrix drop tolerance");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Drop_Tolerance, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Drop_Tolerance, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix polynomial order");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Polynomial_Order, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Polynomial_Order, default_string); ECHO(echo_string,echo_file);
    }

  /*
   * New Aztec 2.0 option for using RCM reordering, or not. While it seems
   * more problems benefit from reordering than are hurt by reordering, the
   * no reorder option will be our default, since that is backward compatible
   * with Aztec 1.1 behavior and "tradition is important".
   */

  strcpy(search_string, "Matrix reorder");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Reorder, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Reorder, default_string); ECHO(echo_string,echo_file);
    }

  /*
   * Matrix factorization save info...(Aztec 2)
   */

  strcpy(search_string, "Matrix factorization save");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Factorization_Save, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Factorization_Save, default_string); ECHO(echo_string,echo_file);
    }

  /*
   * Matrix ILUT factor fill ... (Aztec 2)
   */

  strcpy(search_string, "Matrix ILUT fill factor");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_ILUT_Fill_Factor, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_ILUT_Fill_Factor, default_string); ECHO(echo_string,echo_file);
    }

  /*
   * Matrix RILU relaxation factor ... (Aztec 2)
   */

  strcpy(search_string, "Matrix RILU relax factor");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_RILU_Relax_Factor, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_RILU_Relax_Factor, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix BILU Threshold");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 )
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_BILU_Threshold, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_BILU_Threshold, default_string); ECHO(echo_string,echo_file);
    }

#ifdef TRILINOS
  strcpy(search_string, "Matrix Relative Threshold");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 )
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Relative_Threshold, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Relative_Threshold, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Matrix Absolute Threshold");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 )
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Absolute_Threshold, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Absolute_Threshold, default_string); ECHO(echo_string,echo_file);
    }
#endif

  strcpy(search_string, "Size of Krylov subspace");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Krylov_Subspace, input);
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Krylov_Subspace, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Orthogonalization");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Orthogonalization, input); /* save for aztec */
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Orthogonalization, default_string); ECHO(echo_string,echo_file);
    }

  strcpy(search_string, "Maximum Linear Solve Iterations");
  iread = look_for_optional(ifp, search_string, input, '=');
  if ( iread == 1 ) 
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Maximum_Iterations, input); /* save for aztec */
      SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
    }
  else
    {
      SPF(echo_string, def_form,search_string, Matrix_Maximum_Iterations, default_string); ECHO(echo_string,echo_file);
    }
		 
	strcpy(search_string, "Amesos Solver Package");
	iread = look_for_optional(ifp, search_string, input, '=');
	if ( iread == 1 ) 
	{
		read_string(ifp,input,'\n');
		strip(input); stringup(input);
		strcpy( Amesos_Package, input);
		SPF(echo_string,eoformat, search_string, input); ECHO(echo_string,echo_file);
	}
	else
	{
		SPF(echo_string, def_form,search_string, "KLU", default_string); ECHO(echo_string,echo_file);
	}
		 
		 

  /* first initialize modified newton parameter to false */
  modified_newton = FALSE;

  look_for(ifp, "Number of Newton Iterations", input, '=');
  if (fscanf(ifp, "%d", &Max_Newton_Steps) != 1)
    {
      EH( -1, "error reading Number of Newton Iterations");
    }
  SPF(echo_string,"%s = %d", "Number of Newton Iterations", Max_Newton_Steps);

  if (fscanf(ifp, "%d", &Newt_Jacobian_Reformation_stride) == 1)
    {
      if(Newt_Jacobian_Reformation_stride > 1) modified_newton = TRUE;

      SPF(endofstring(echo_string), " %d", Newt_Jacobian_Reformation_stride );
    }
  else
    {
       Newt_Jacobian_Reformation_stride = 1;
    }

  ECHO(echo_string,echo_file);
    
  iread = look_for_optional(ifp, "Modified Newton Tolerance", input, '=');
  if (iread == 1) { 
    if (fscanf(ifp, "%le %le", &convergence_rate_tolerance, &modified_newt_norm_tol) != 2)
      {
	EH( -1, "ERROR reading Modified Newton Tolerance card");
      }
    else
      {
        if (convergence_rate_tolerance != 0. && modified_newt_norm_tol != 0.)
	  {
	    modified_newton=TRUE;
	  }
	SPF(echo_string,"%s = %.4g %.4g", "Modified Newton Tolerance",convergence_rate_tolerance , modified_newt_norm_tol ); 
	ECHO(echo_string,echo_file);
      }
  }
  else
    {
      convergence_rate_tolerance = 0.;
      modified_newt_norm_tol = 0.;
    }

  iread = look_for_optional(ifp, "Jacobian Reform Time Stride", input, '=');
  if (iread == 1) { 
    if (fscanf(ifp, "%d", &Time_Jacobian_Reformation_stride) != 1)
      {
	EH( -1, "ERROR reading Jacobian Reform Stride card");
      }
    else
      {
	if(Time_Jacobian_Reformation_stride > 1) modified_newton = TRUE;
	SPF(echo_string, "%s = %d", "Jacobian Reform Time Stride", Time_Jacobian_Reformation_stride);
	ECHO(echo_string,echo_file);
      }
  } 
  else
    {
      Time_Jacobian_Reformation_stride = 0;
    }



  look_for(ifp, "Newton correction factor", input, '=');
  if (fscanf(ifp, "%le", &damp_factor1) != 1 )
    {
      EH( -1, "error reading Newton correction factor, expected at least one float");
    }

  SPF(echo_string,"%s = %.4g", "Newton correction factor", damp_factor1 );

  if (fscanf(ifp, "%le %le %le %le %le", &custom_tol1, &damp_factor2, 
                  &custom_tol2, &damp_factor3, &custom_tol3) != 5 )
    {
      damp_factor2 = damp_factor3 =  -1.;
      custom_tol1 = custom_tol2 = custom_tol3 = -1;
      rewind (ifp);    /* Added to make ibm happy when single relaxation value input. dal/1-6-99 */
    }
  else
   {
      if((damp_factor1 <= 1. && damp_factor1 >= 0.) &&
         (damp_factor2 <= 1. && damp_factor2 >= 0.) &&
         (damp_factor3 <= 1. && damp_factor3 >= 0.))
        {
	  SPF(endofstring(echo_string)," %.4g %.4g %.4g %.4g %.4g", custom_tol1, damp_factor2,custom_tol2,damp_factor3,custom_tol3);
        }
      else
        {
          EH(-1,"All damping factors must be in the range 0 <= fact <=1");
        }
   }

  ECHO(echo_string,echo_file);


   /* initialize variable specific damping */
   for (k=0; k<MAX_VARIABLE_TYPES; k++)
      {
        var_damp[k] = 1.;
      }

  while ((iread = look_forward_optional(ifp,"Variable Type Newton correction factor",input,'=')) == 1)
    {
      int var;
      var = -1;

      if (fscanf(ifp, "%80s", input) != 1)
        {
          EH(-1, "Error reading variable type for Variable Type Newton correction factor");
        }

      /* loop through variable names and compare to input */
	
      for (k=0; k<MAX_VARIABLE_TYPES && var == -1; k++)
        {
          if ( !strcmp(input, Var_Name[k].name1) ||
    	       !strcmp(input, Var_Name[k].name2)) var = k;
        }

      if (var == -1)
        {
          EH(-1,"Could not identify variable type for Variable Type Newton Correction Factor!");
        }
	
      if (fscanf(ifp, "%le", &var_damp[var]) != 1)
        {
          EH( -1, "error reading value of variable type Newton correction factor");
        }
      SPF(echo_string,"%s = %s %.4g","Variable Type Newton correction factor" , Var_Name[k].name1, var_damp[var]);
      ECHO(echo_string,echo_file);
    }

  look_for(ifp, "Normalized Residual Tolerance", input, '=');
  if (fscanf(ifp, "%le", &Epsilon[0]) != 1)
    {
      EH( -1, "error reading Normalized (Newton) Residual Tolerance");
    }

  SPF(echo_string,"%s = %.4g", "Normalized Residual Tolerance", Epsilon[0]); ECHO(echo_string,echo_file);

  iread = look_for_optional(ifp, "Normalized Correction Tolerance", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le", &Epsilon[2]) != 1)
        {
	  EH( -1, "error reading Normalized (Newton) Correction Tolerance");
        }
      SPF(echo_string,"%s = %.4g", "Normalized Correction Tolerance", Epsilon[2] ); ECHO(echo_string,echo_file);
    }
    else
    {
      Epsilon[2]=1.0e+10;
    }

  iread = look_for_optional(ifp, "Residual Ratio Tolerance", input, '=');
  if ( iread == 1 )
    {
      read_string(ifp,input,'\n');
      strip(input);
      strcpy(Matrix_Convergence_Tolerance, input); /* save for aztec use */
      if ( sscanf(input, "%le", &Epsilon[1]) != 1 )
	{
	  EH( -1, "Bad residual ratio (matrix convergence) tolerance.");
	}
      SPF(echo_string,"%s = %.4g", "Residual Ratio Tolerance", Epsilon[1] ); ECHO(echo_string,echo_file);
    }

  iread = look_for_optional(ifp, "Pressure Stabilization", input, '=');
  if (iread == 1) { 
    (void) read_string(ifp, input, '\n');
    strip(input);
    if ( strcmp(input,"no") == 0 )
      {
	PSPG = 0;
        PSPP = 0;
      }
    else if ( strcmp(input,"yes") == 0 || strcmp(input,"global") == 0 )
      {
	PSPG = 1;
        PSPP = 0;
      }
    else if ( strcmp(input,"local") == 0 )
      {
	PSPG = 2;
        PSPP = 0;
      }
    else if ( strcmp(input,"pspp") == 0 )
      {
	PSPG = 0;
        PSPP = 1;
      }
    else if ( strcmp(input,"pspp_e") == 0 )
      {
	PSPG = 0;
        PSPP = 2;
      }
    else
      {
	EH( -1, "invalid choice: Pressure Stabilization yes, no, global or local, or pspp");
      }
    SPF(echo_string,eoformat, "Pressure Stabilization", input ); ECHO(echo_string,echo_file);
  }
  else
    {
      PSPG = 0;
      PSPP = 0;
      ECHO("(Pressure Stabilization = NONE) (default)", echo_file);
    }

  iread = look_for_optional(ifp, "Pressure Stabilization Scaling", input, '=');
  if (iread == 1) { 
    if (fscanf(ifp, "%le", &PS_scaling) != 1)
      {
	EH( -1, "error reading Pressure Stabilization Scaling");
      }
    SPF(echo_string,"%s = %.4g", "Pressure Stabilization Scaling", PS_scaling ); ECHO(echo_string,echo_file);
  }
  else
    {
      PS_scaling = 0.;
      ECHO("(Pressure Stabilization Scaling = 0.0) (default)", echo_file);
    }

  iread = look_for_optional(ifp, "Continuity Stabilization", input, '=');
  if (iread == 1)
    {
      (void) read_string(ifp, input, '\n');
      strip(input);
      if (strcmp(input,"no") == 0)
	{
	  Cont_GLS = 0;
	}
      else if(strcmp(input,"yes") == 0)
	{
	  Cont_GLS = 1;
	}

      else if(strcmp(input,"local")==0)
	{
          Cont_GLS = 2;
	}
       
      else
	{
	  EH( -1, "invalid choice: Continuity Stabilization yes, local, or no");
	}
      SPF(echo_string, eoformat, "Continuity Stabilization", input); ECHO(echo_string,echo_file);	  
	
    }
  else
    {
      Cont_GLS = 0;
      ECHO("(Continuity Stabilization = None) (default)", echo_file);
    }

  /*IGBRK*/
  iread = look_for_optional(ifp, "Linear Stability", input, '=');
  if (iread == 1)
    { 
      (void) read_string(ifp, input, '\n');
      strip(input);
      if(!strncmp(input, "no", 2))
	Linear_Stability = LSA_NONE;
      else if(!strncmp(input, "yes", 3) || !strncmp(input, "inline", 6))
	Linear_Stability = LSA_NORMAL;
      else if(!strncmp(input, "file", 4))
	Linear_Stability = LSA_SAVE;
      else if(!strncmp(input, "3D", 2))
	Linear_Stability = LSA_3D_OF_2D;
      else if(!strncmp(input, "3Dfile", 6))
	Linear_Stability = LSA_3D_OF_2D_SAVE;
      else
	EH( -1, "Invalid choice: Linear Stability must be yes, no, file, 3D, or 3Dfile");

      SPF(echo_string,eoformat, "Linear Stability", input); ECHO(echo_string,echo_file);

    }
  else
    Linear_Stability = LSA_NONE;

  iread = look_for_optional(ifp, "Filter Concentration", input, '=');
  if (iread == 1) 
    { 
      Filter_Species = TRUE;
      if (fscanf(ifp, "%d %le %le", &filter_species_material_number, &c_min, &c_max) != 3)
	{
	  EH( -1, "error reading Filter Concentration");
	}
      SPF(echo_string,"%s = %d %.4g %.4g","Filter Concentration", filter_species_material_number, c_min, c_max );
      ECHO(echo_string,echo_file); 
    }
  else
    {
      Filter_Species = FALSE;
    }

  iread = look_for_optional(ifp, "Disable Viscosity Sensitivities", input, '=');
  if (iread == 1)
    {
    (void) read_string(ifp, input, '\n');
    strip(input);
    if ( strncmp(input,"no",2) == 0 )
      {
        Include_Visc_Sens = TRUE;
      }
    else if ( strncmp(input,"yes",3) == 0 )
      {
	char temp_string[80];
        Include_Visc_Sens = FALSE;
        if (sscanf(input, "%s %d", temp_string, &Visc_Sens_Factor) != 2)
	    {  Visc_Sens_Factor = 2; }
      }
    else
      {
        EH( -1, "invalid choice for Disable Viscosity Sensitivities:must be yes or no");
      }
      SPF(echo_string,eoformat,"Disable Viscosity Sensitivities", input); ECHO(echo_string,echo_file);
    }
  else
    {
      Include_Visc_Sens = TRUE;
    }

  Visc_Sens_Copy = Include_Visc_Sens;
  
  /* look for optional flags specifying dependencies to ignore */
  {
    int err;
    int eq, var;
    
    for (eq=0; eq<MAX_VARIABLE_TYPES; eq++)
      {
	for (var=0; var<MAX_VARIABLE_TYPES; var++)
	  {
	    Ignore_Deps[eq][var] = FALSE;
	  }
      }
    
    while ((iread = look_forward_optional(ifp,"Ignore Dependency",input,'=')) == 1) {
      int eq_found = FALSE;
      int var_found = FALSE;
      int symmetric_flag = FALSE;

      if (fscanf(ifp, "%80s", input) != 1)
	{
	  EH(-1, "Error reading equation type for Ignore Dependency");
	}
      /* loop through equation names and compare to input */
      for (eq=0; eq<MAX_VARIABLE_TYPES && !eq_found; eq++)
	{
	  if ( !strcmp(input, EQ_Name[eq].name1) ||
	       !strcmp(input, EQ_Name[eq].name2))
	    {
	      eq_found = TRUE;
	    }
	}
      if (!eq_found)
	{
	  EH(-1, "Could not identify equation type for Ignore Dependency");
	}
      eq--; /* this got incremented one more time than desired */
      if ((fscanf(ifp, "%2s", input) != 1) || 
          (strcmp( input, "on") && strcmp( input, "ON")))
	{
	  EH(-1, "Syntax error for Ignore Dependency");
	}
      if (fscanf(ifp, "%80s", input) != 1)
	{
	  EH(-1, "Error reading variable type for Ignore Dependency");
	}
      /* look for optional flag to apply this card in a symmetric manner */
      err = fscanf(ifp, "%d",&symmetric_flag);

      errno = 0;
      if ((err != 0 || err != 1) && (err == EOF && errno != 0)) {
	EH(-1, "Error reading symmetric flag for Ignore Dependency");
      }
      
      if ( !strcmp(input, "all") || !strcmp(input, "ALL") )
        {
	  for (var=0; var<MAX_VARIABLE_TYPES; var++)
	    {
	      if ( var != eq ) /* keep this entry only */
	        {
	          if ( symmetric_flag )
	            {
	              Ignore_Deps[eq][var] = TRUE;
	              Ignore_Deps[var][eq] = TRUE;
	            }
                  else
	            {
	              Ignore_Deps[eq][var] = TRUE;
		    }
	        }
	    }
	}
      else
        {
	
          /* loop through equation names and compare to input */
          for (var=0; var<MAX_VARIABLE_TYPES && !var_found; var++)
	    {
	      if ( !strcmp(input, Var_Name[var].name1) ||
	           !strcmp(input, Var_Name[var].name2))
	        {
	          var_found = TRUE;
	        }
	    } 
          if (!var_found)
	    {
	      EH(-1, "Could not identify variable type for Ignore Dependency");
	    }
          var--; /* this got incremented one more time than desired */
      
          if ( symmetric_flag )
	    {
	      Ignore_Deps[eq][var] = TRUE;
	      Ignore_Deps[var][eq] = TRUE;
	    }
          else
	    {
	      Ignore_Deps[eq][var] = TRUE;
	    }
	}
    }
  }


#ifdef PARALLEL
  if( ( Num_Proc > 1 ) &&  ( is_Solver_Serial == TRUE ) )
    {
    EH( -1, "Your solver choice cannot be used for PARALLEL solutions");
    }
#endif

}
/* rd_solver_specs -- read input file for equation solver specifications */

/*
 * rd_eigen_specs -- read input file for eigensolver specifications
 *
 * Created:			1997/09/28 IDG
 * Adapted from rd_solver_specs()
 */

void 
rd_eigen_specs(FILE *ifp,
		 char *input)
{
  int i;
  int iread;
  char copy_of_input[MAX_CHAR_IN_INPUT];
  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;
  /*  */

/* READ IN EIGENSOLVER SPECIFICATIONS */
  
  if (Linear_Stability == LSA_NONE) return;
  iread = look_for_optional(ifp,"Eigensolver Specifications",input,'=');
  
  if( iread == 1) {
	  SPF(echo_string,"%s =", input);ECHO(echo_string,echo_file);
  }

/* SELECT ARPACK METHOD */
  eigen->Eigen_Algorithm = LSA_DEFAULT;
  iread = look_for_optional(ifp,"Eigen Algorithm",input,'=');
  if (iread == 1)
    {
	  SPF(echo_string,"%s =",input);
      (void) read_string(ifp, input, '\n');
      strip(input);

/* Handle based on selected continuation */
  switch (Continuation) {

/* LOCA continuation: default to shift and invert */
    case LOCA:
      if (strcmp(input, "cayley") == 0)
        {
          eigen->Eigen_Algorithm = LSA_CAYLEY;
        }
      else
        {
          eigen->Eigen_Algorithm = LSA_SI;
        }
      break;

/* No continuation: Use ARPACK if card present, otherwise eggroll */
    case ALC_NONE:
      if (strcmp(input, "cayley") == 0)
        {
          eigen->Eigen_Algorithm = LSA_CAYLEY;
          loca_in->Cont_Alg = LOCA_LSA_ONLY;
        }
      else if (strcmp(input, "si") == 0)
        {
          eigen->Eigen_Algorithm = LSA_SI;
          loca_in->Cont_Alg = LOCA_LSA_ONLY;
        }
      else eigen->Eigen_Algorithm = LSA_DEFAULT;

/* Initialize these to avoid UMR's: */
      if (loca_in->Cont_Alg == LOCA_LSA_ONLY)
        {
          cont->BegParameterValue = 0.0;
          cont->EndParameterValue = 0.0;
          cont->Delta_s0 = 0.0;
          cont->MaxPathSteps = 0;
          cont->print_freq = 1;
          cont->upType = LOCA_LSA_ONLY;
        }
      break;

    case HUN_ZEROTH:
    case HUN_FIRST:
fprintf(stderr,"HUN %d %d %d\n",nHC,hunt[0].Type,hunt[0].BCID);
      if (strcmp(input, "cayley") == 0)
        {
          eigen->Eigen_Algorithm = LSA_CAYLEY;
        }
      else if (strcmp(input, "si") == 0)
        {
          eigen->Eigen_Algorithm = LSA_SI;
        }
      else eigen->Eigen_Algorithm = LSA_DEFAULT;

/* Initialize these to avoid UMR's: */
      if (0 &&  loca_in->Cont_Alg == LOCA_LSA_ONLY)
        {
          cont->BegParameterValue = 0.0;
          cont->EndParameterValue = 0.0;
          cont->Delta_s0 = 0.0;
          cont->MaxPathSteps = 0;
          cont->print_freq = 1;
          cont->upType = LOCA_LSA_ONLY;
        }
      break;
/* For non-LOCA continuation, this card does not apply */
    default:
      fprintf(stdout, "LSA not available for continuation without LOCA!");
      Linear_Stability = LSA_NONE;  
      break;
    }
  }
  
  SPF(endofstring(echo_string)," %s",input); ECHO(echo_string,echo_file);

/* NUMBER OF WANTED MODES [eggroll or ARPACK] */

  iread = look_for_optional(ifp, "Eigen Number of modes", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%d", &eigen->Eigen_NEV_WANT) != 1) iread = -1;
    }
  if (iread != 1)
    {
      eigen->Eigen_NEV_WANT = 10;
      SPF(echo_string, "\t(%s = %d)", "Eigen Number of modes",eigen->Eigen_NEV_WANT);
    }
  else
  {
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_NEV_WANT);
  }
  
  ECHO(echo_string,echo_file);
 
/* NUMBER OF MODES TO BE RECORDED [eggroll or ARPACK] */
 
  iread = look_for_optional(ifp, "Eigen Record modes", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%d", &eigen->Eigen_Record_Modes) != 1) iread = -1;
    }
  if (iread != 1)
    {
      eigen->Eigen_Record_Modes = 0;
      SPF(echo_string, "\t(%s = %d)", "Eigen Record modes",eigen->Eigen_Record_Modes);
    }
  else 
  {
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Record_Modes);
  }
  
  ECHO(echo_string,echo_file);

/* KRYLOV SUBSPACE [eggroll or ARPACK] */

  iread = look_for_optional(ifp, "Eigen Size of Krylov subspace", input, '=');
  if (iread == 1)
    { 
      if (fscanf(ifp, "%d", &eigen->Eigen_Krylov_Subspace) != 1) iread = -1;
    }
  if (iread != 1)
    {
      eigen->Eigen_Krylov_Subspace = 30;
      SPF(echo_string, "\t(%s = %d)", "Eigen Size of Krylov subspace",eigen->Eigen_Krylov_Subspace);
    }
  else
  {
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Krylov_Subspace);	  
  }
  
  ECHO(echo_string,echo_file);

/* MAX ITERATIONS [eggroll only] */

  iread = look_for_optional(ifp, "Eigen Maximum Iterations", input, '=');
  if (iread == 1)
    { 
      if (fscanf(ifp, "%d", &eigen->Eigen_Maximum_Iterations) != 1) iread = -1;
    }
  if (iread != 1)
    {
      eigen->Eigen_Maximum_Iterations = 100;
      SPF(echo_string, "\t(%s = %d)", "Eigen Maximum Iterations",eigen->Eigen_Maximum_Iterations);
    }
  else
  {
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Maximum_Iterations)	;  	  
  }

  ECHO(echo_string,echo_file);

/* FILTER [eggroll only] */

  iread = look_for_optional(ifp, "Eigen Number of Filter Steps", input, '=');
  if (iread == 1)
    { 
      if (fscanf(ifp, "%d", &eigen->Eigen_Filter) != 1) iread = -1;
    }
  if (iread != 1)
    {
      eigen->Eigen_Filter = 2;
      SPF(echo_string, "\t(%s = %d)", "Eigen Number of Filter Steps",eigen->Eigen_Filter);
    }
  else
  {
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Filter)	;  	  	  
  }
  
  ECHO(echo_string,echo_file);

/* RECYCLE [eggroll only] */

  iread = look_for_optional(ifp, "Eigen Recycle", input, '=');
  if (iread == 1)
    { 
	  SPF(echo_string,"%s =", input);
      (void) read_string(ifp, input, '\n');
      strip(input);
      if ( strcmp(input,"no") == 0 )
        {
          eigen->Eigen_Recycle = 0;
        }
      else if ( strcmp(input,"yes") == 0 )
        {
          eigen->Eigen_Recycle = 1;
        }
      else
        {
          EH( -1, "invalid choice: Eigen Recycle must be yes or no");
        }
	  SPF(endofstring(echo_string)," %s",input);
    }
  else
    {
      eigen->Eigen_Recycle = 0;
      SPF(echo_string, "\t(%s = %s)", "Eigen Recycle","no");
    }

  ECHO(echo_string,echo_file);

  /* TOLERANCE */

  iread = look_for_optional(ifp, "Eigen Tolerance", input, '=');
  if (iread == 1)
    {
      eigen->Eigen_Tolerance = read_dbl(ifp, "Eigen Tolerance");
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_Tolerance)	;  	  	  
    }
  else
    {
      eigen->Eigen_Tolerance = 1.0e-06;
      SPF(echo_string, "\t(%s = %.4g)", "Eigen Tolerance",eigen->Eigen_Tolerance);
    }
  
  ECHO(echo_string,echo_file);

/* MATRIX OUTPUT FLAG [eggroll or ARPACK] */

  iread = look_for_optional(ifp, "Eigen Matrix Output", input, '=');
  if (iread == 1)
    {
	  SPF(echo_string,"%s =", input);
      (void) read_string(ifp, input, '\n');
      strip(input);
      if ( strcmp(input, "yes") == 0)
        {
          eigen->Eigen_Matrix_Output = 1;
        }
      else
        {
          eigen->Eigen_Matrix_Output = 0;
        }
	  SPF(endofstring(echo_string)," %s",input);
    }
  else
    {
      eigen->Eigen_Matrix_Output = 0;
      SPF(echo_string, "\t(%s = %s)","Eigen Matrix Output","no");
    }
  
  ECHO(echo_string,echo_file);

/* INITIAL VECTOR WEIGHT [eggroll only] */

  iread = look_for_optional(ifp, "Eigen Initial Vector Weight", input, '=');
  if (iread == 1)
    { 
      if (fscanf(ifp, "%le", &eigen->Eigen_IV_Wt) != 1)
        {
	  EH( -1, "error reading Eigen Initial Vector Weight");
        }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_IV_Wt)	;  	  	  
    }
  else
    {
      eigen->Eigen_IV_Wt = 0.5;
      SPF(echo_string, "\t(%s = %.4g)", "Eigen Tolerance",eigen->Eigen_IV_Wt);
    }
  
  ECHO(echo_string,echo_file);

/* INITIAL SHIFTS [eggroll only] */

  for (i=0; i<4; i++)
    {
      eigen->Eigen_Shifts[i] = -1.0;
    }
  iread = look_for_optional(ifp, "Eigen Initial Shifts", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le %le %le %le", 
                 &eigen->Eigen_Shifts[0],
                 &eigen->Eigen_Shifts[1],
                 &eigen->Eigen_Shifts[2],
                 &eigen->Eigen_Shifts[3]) != 4)
        {
          EH(-1,"Four float values required for Eigen Initial Shifts\n");
        }
	  SPF(echo_string,"%s =",input);
	  SPF_DBL_VEC(endofstring(echo_string), 4, eigen->Eigen_Shifts);
    }
  else
	 {
		SPF(echo_string,"\t(%s = ","Eigen Initial Shifts");
		SPF_DBL_VEC(endofstring(echo_string), 4, eigen->Eigen_Shifts);
		SPF(endofstring(echo_string),")");
	 }
  ECHO(echo_string,echo_file);


/* CAYLEY SIGMA AND MU [ARPACK only] */

  if (eigen->Eigen_Algorithm == LSA_CAYLEY)
    {
      eigen->Eigen_Cayley_Sigma = 100.0;
    }
  else
    {
      eigen->Eigen_Cayley_Sigma = 0.0;
    }
  iread = look_for_optional(ifp, "Eigen Cayley Sigma", input, '=');
  if (iread == 1)
  {
      if (fscanf(ifp, "%le", &eigen->Eigen_Cayley_Sigma) != 1)
	  {
          EH(-1, "Error reading Eigen Cayley Sigma");
	  }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_Cayley_Sigma)	;  	  	  
  }
  else
  {
	  SPF(echo_string, "\t(%s = %.4g)", "Eigen Cayley Sigma",eigen->Eigen_Cayley_Sigma);
  }
  
  ECHO(echo_string,echo_file);
  
  eigen->Eigen_Cayley_Mu = 1000.0;
  iread = look_for_optional(ifp, "Eigen Cayley Mu", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le", &eigen->Eigen_Cayley_Mu) != 1)
        {
          EH(-1, "Error reading Eigen Cayley Mu");
        }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_Cayley_Mu)	;  	  	  
    }
  else
  {
	  SPF(echo_string, "\t(%s = %.4g)", "Eigen Cayley Mu",eigen->Eigen_Cayley_Mu);
  }
  
  ECHO(echo_string,echo_file);
  
/* MAX OUTER ITERATIONS [ARPACK only] */
                                                                                
  iread = look_for_optional(ifp, "Eigen Maximum Outer Iterations", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%d", &eigen->Eigen_Maximum_Outer_Iterations) != 1) 
		  iread = -1;
	  else
		  SPF(echo_string,"%s = %d", input, eigen->Eigen_Maximum_Outer_Iterations)	;  	  	  
		  
    }
  
  if (iread != 1)
    {
      eigen->Eigen_Maximum_Outer_Iterations = 1;
	  SPF(echo_string, "\t(%s = %d)", "Eigen Maximum Outer Iterations",eigen->Eigen_Maximum_Outer_Iterations);
    }

  ECHO(echo_string,echo_file);

  
/* EIGENSOLVER SHIFT/INVERT TOLERANCE PARAMETER [ARPACK only] */
                                                                                
  eigen->Eigen_SI_Tol_Param = 2.0e-2;
  iread = look_for_optional(ifp, "Eigen SI tolerance parameter", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le", &eigen->Eigen_SI_Tol_Param) != 1)
        {
          EH(-1, "Error reading Eigen SI tolerance parameter");
        }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_SI_Tol_Param)	;  	  	  
    }
  else
	  SPF(echo_string, "\t(%s = %.4g)", "Eigen SI tolerance parameter",eigen->Eigen_SI_Tol_Param);

  ECHO(echo_string,echo_file);

/* EIGENSOLVER RELATIVE TOLERANCE [ARPACK only] */

  eigen->Eigen_Relative_Tol = 1.0e-6;
  iread = look_for_optional(ifp, "Eigen Relative tolerance", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le", &eigen->Eigen_Relative_Tol) != 1)
        {
          EH(-1, "Error reading Eigen Relative tolerance");
        }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_Relative_Tol)	;  	  	  
    }
  else
	  SPF(echo_string, "\t(%s = %.4g)", "Eigen Relative tolerance",eigen->Eigen_Relative_Tol);
  
  ECHO(echo_string,echo_file);
  

/* EIGEN LINEAR SOLVER TOLERANCE [ARPACK only] */

  eigen->Eigen_Linear_Tol = 1.0e-6;
  iread = look_for_optional(ifp, "Eigen Linear Solver tolerance", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%le", &eigen->Eigen_Linear_Tol) != 1)
        {
          EH(-1, "Error reading Eigen Linear Solver tolerance");
        }
	  SPF(echo_string,"%s = %.4g", input, eigen->Eigen_Linear_Tol)	;  	  	  
    }
  else
	  SPF(echo_string, "\t(%s = %.4g)", "Eigen Linear Solver tolerance",eigen->Eigen_Linear_Tol);
  
  ECHO(echo_string,echo_file);
  

/* EIGENVALUE FREQUENCY [ARPACK only] */
  eigen->Eigen_Solve_Freq = 1;
  iread = look_for_optional(ifp, "Eigenvalue output frequency", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%d", &eigen->Eigen_Solve_Freq) != 1)
        {
          EH(-1, "Error reading Eigenvalue output frequency");
        }
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Solve_Freq)	;  	  	  
    }
  else
	  SPF(echo_string, "\t(%s = %d)", "Eigenvalue output frequency",eigen->Eigen_Solve_Freq);
  
  ECHO(echo_string,echo_file);


/* EIGENVECTOR FREQUENCY [ARPACK only] */
  eigen->Eigen_Write_Freq = -1;
  iread = look_for_optional(ifp, "Eigenvector output frequency", input, '=');
  if (iread == 1)
    {
      if (fscanf(ifp, "%d", &eigen->Eigen_Write_Freq) != 1)
        {
          EH(-1, "Error reading Eigenvector output frequency");
        }
	  SPF(echo_string,"%s = %d", input, eigen->Eigen_Write_Freq)	;  	  	  
    }
  else
	  SPF(echo_string, "\t(%s = %d)", "Eigenvector output frequency",eigen->Eigen_Write_Freq);
  
  ECHO(echo_string,echo_file);
  

/* EIGENVECTOR OUTPUT EXODUS FILE [ARPACK only] */

  iread = look_for_optional(ifp, "Eigenvector output file", input, '=');
  if (iread == 1)
    {
	  SPF(echo_string,"%s =",input);
      (void) read_string(ifp, input, '\n');
      strip(input);
      strcpy(eigen->Eigen_Output_File, input);
      SPF(endofstring(echo_string)," %s",eigen->Eigen_Output_File);
    }
  else
    {
      strcpy(eigen->Eigen_Output_File, "LSA.exoII");	  
	  SPF(echo_string, "\t(%s = %s)", "Eigenvector output file",eigen->Eigen_Output_File);
    }
  
  ECHO(echo_string,echo_file);

  /* Which normal modes to compute?  Give the wave numbers as a
   * whitespace-separated list. */
  iread = look_for_optional(ifp, "Eigen Wave Numbers", input, '=');
  if(iread == 1)
  {
	  SPF(echo_string,"%s =",input);
      (void) read_string(ifp, input, '\n');
      strip(input);
      strncpy(copy_of_input, input, MAX_CHAR_IN_INPUT - 1); /* leave room for trailing \0 */
      if(strtok(copy_of_input, " \t\n"))
		  LSA_number_wave_numbers++;
      else
		  EH(-1, "Bad Eigen Wave Numbers card");
      while(strtok(NULL, " \t\n"))
		  LSA_number_wave_numbers++;
      LSA_wave_numbers = (dbl *)malloc(LSA_number_wave_numbers * sizeof(dbl));
      if(!LSA_wave_numbers)
		  EH(-1, "Could not allocate space for LSA_wave_numbers");
      LSA_wave_numbers[0] = atof(strtok(input, " \t\n"));
      if (LSA_wave_numbers[0] == 0.) EH(-1,"Cannot solve 3D of 2D for zero wave number");
      for(i = 1; i < LSA_number_wave_numbers; i++)
            {
		  LSA_wave_numbers[i] = atof(strtok(NULL, " \t\n"));
	          if (LSA_wave_numbers[i] == 0.) EH(-1, "cannot solve 3D of 2D for zero wave number");
            }
	  
	  SPF_DBL_VEC(endofstring(echo_string),LSA_number_wave_numbers,  LSA_wave_numbers);
  }
  else if(Linear_Stability == LSA_3D_OF_2D ||
	  Linear_Stability == LSA_3D_OF_2D_SAVE)
    EH(-1, "3D stability of 2D flow requested, but missing \"Eigen Wave Numbers\" card");
  
  ECHO(echo_string,echo_file);
}  /* End of rd_eigen_specs */



/*****************************************************************************/
/*  NOTE:  rd_bc_specs has been moved to file mm_input_bc.c!                 */
/*****************************************************************************/

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/*
 * rd_matl_blk_specs -- read input file to read and guide the problem
 *                      description input, equation input, and material
 *                      property input.  
 *
 *
 * Comments:	This routine calls rd_eq_specs and rd_mp_specs which used
 *              to be called at the highest level of mm_input.
 *
 * Created:			Thr Feb 16 15:24:26 MST 1995 prschun@sandia.gov 
 */

void 
rd_matl_blk_specs(FILE *ifp,
		  char *input)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  int	err, mn, i;
  FILE *imp;
  char MatFile[MAX_FNL];	/* Raw material database file. */
  char TmpMatFile[MAX_FNL];	/* Temporary copy of mat db after APREPRO. */

  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_input_file = Echo_Input_File;
  char echo_mat_file[MAX_FNL]="\0";

  static char MatFileSuffix[] = ".mat";
  static char System_Command[MAX_SYSTEM_COMMAND_LENGTH];

  /*
   * Identify section containing equation specification...
   */

  Use_DG = FALSE;
  look_for_optional(ifp, "Problem Description", input, '=');
  ECHO("\nProblem Description =\n", echo_input_file);

  look_for(ifp, "Number of Materials", input, '=');
  if ( fscanf(ifp,"%d",&(upd->Num_Mat)) != 1)
    {
      fprintf(stderr,
	      "Expected to read 1 int for \"Number of Materials\" - attempt to enumerate...");
      upd->Num_Mat = -1;
    }

    
  /*
   * Count Number of Materials if Num_Mat is set to -1
   */
  
  if (upd->Num_Mat == -1)
      upd->Num_Mat = count_list(ifp, "MAT", input, '=', "END OF MAT");
  /*
   * Check the number of materials against the default Maximum number of materials
   */
  if (upd->Num_Mat > MAX_NUMBER_MATLS) {
    sprintf(err_msg, 
	    "Found num materials = %d > %d - boost MAX_NUMBER_MATLS in mm_mp_const.h",
	    upd->Num_Mat, MAX_NUMBER_MATLS);
    EH( -1, err_msg);
  }

  SPF(echo_string,"%s = %d","Number of Materials",  upd->Num_Mat); ECHO(echo_string,echo_input_file);
       
  upd->Max_Num_Species_Eqn = -123456789; /* initialize */
  upd->Max_Num_Species     = -123456789; /* initialize */
  upd->Max_Num_Porous_Eqn  = -123456789; /* initialize */
  
  Num_Interpolations = 0; 	/* Initialize counter (cf rf_fem.h) */

  for ( i=0; i<MAX_INTERPOLATIONS; i++)
    {
      Unique_Interpolations[i] = I_NOTHING;
    }

  upd->Total_Num_EQ = 0;
  upd->Total_Num_Var = 0;


  for ( i=0; i<MAX_VARIABLE_TYPES+MAX_CONC; i++)
    {
      upd->ep[i] = -1;
      upd->vp[i] = -1;
    }

#ifndef COUPLED_FILL
  /* Initialized Explicit_Fill to zero. Only set it if
   * we have a fill equation in some material
   */  
  Explicit_Fill = 0;
#endif

  /*
   * ------------- Loop over the number of materials ----------
   */
  for (mn = 0; mn < upd->Num_Mat; mn++) 
    { /* mn */
      /* 
       *
       * Part 1 -- Look for Material Name (Arbitrary) and 
       *           corresponding ExodusII blk ID's associated
       *           with that material
       */
      read_MAT_line(ifp, mn, input);

      /*
       * Search for any volume ACs assigned to an element block
       * associated with this material. If a match is found, then
       * set the VolumeIntegral
       * flag in the material property structure to point to the
       * augmented condition index associated with the Volume Integral
       * constraint.
       */
      if (nAC > 0) {
	int iAC ;	
	for (iAC = 0; iAC < nAC; iAC++)
	  {
	    if (augc[iAC].Type == AC_VOLUME) 
	      {
		if (eb_in_matrl(augc[iAC].MTID, mn))
		  {
		    pd_glob[mn]->VolumeIntegral = iAC;
		  }
	      }
	    else if (augc[iAC].Type == AC_LS_VEL) 
	      {
		if (eb_in_matrl(augc[iAC].MTID, mn))
		  {
		    pd_glob[mn]->LSVelocityIntegral = iAC;
		  }
	      }
	  }
      }
      
      /*
       * Read Equation Specifications
       *
       */      
      rd_eq_specs(ifp, input, mn);

  /* run some further diagnostics to see if these variables are being used
   * for running decoupled code-to-code problem types.  The one type we are 
   * supporting right now are JAS-GOMA couples for poroelastic flow. 
   */
      efv->ev_porous_decouple = T_NOTHING;
      if(efv->ev)
	{

	  if ((strcmp(efv->name[0], "DMX") == 0 || strcmp(efv->name[0], "DISPLX") == 0) && 
	      (strcmp(efv->name[1], "DMY") == 0  || strcmp(efv->name[1], "DISPLY") == 0) &&
	      !pd_glob[0]->v[MESH_DISPLACEMENT1]       &&
	      pd_glob[0]->v[POR_LIQ_PRES]                )
	    {
	      efv->ev_porous_decouple = T_SOMETHING;
	    } 
	  else if (!pd_glob[0]->v[MESH_DISPLACEMENT1]  &&
		   pd_glob[0]->v[POR_LIQ_PRES]                )
	    {
	      EH(-1,"You have por_liq_press equation and efv's, but not the right order. Talk to PRS or RAR");
	    }
	}

      /*
       * Initialize mp structure to default unity properties and 
       * zero source terms and zero sensitivities
       */
      set_mp_to_unity(mn); 
      
      /*
       * Read Material Properties File.
       * If there is no file present, set all material properties to
       * unity and control weightings through multipliers in equation
       * specifications
       */       
      (void) strcpy (MatFile, pd_glob[mn]->MaterialName);
      (void) strcat (MatFile, MatFileSuffix);       

      strcpy( echo_mat_file, "echo_");
      strcat( echo_mat_file, MatFile);

      ECHO("OPEN", echo_mat_file);

      /*
       * If the "-a" or "-aprepro" flag was specified on the command line,
       * then assume that the MAT files will also be so preprocessed...
       *
       * Need to make up some temporary file names...
       *
       */

       if ( (imp=fopen(MatFile,"r")) != NULL )
       {
	 strcpy_rtn = strcpy(current_mat_file_name, MatFile);
	 if ( run_aprepro == 1 )
	   {
	     fclose(imp);
	     /* 
	      * Start fresh for each matl, build up command and temporary
	      * file name from available information...
	      */
	     System_Command[0] = '\0';
	     TmpMatFile[0] = '\0';
	     strcat(TmpMatFile, "tmp.");
	     strcat(TmpMatFile, MatFile);
	     sprintf(System_Command, "%s %s %s", aprepro_command,
		     MatFile, TmpMatFile);
#ifndef tflop
	     err = system(System_Command);
	     EH(err, "system() choked on mat file.");

	     if (WEXITSTATUS(err) == 127)
	       {
		 EH(-1, "System call failed, aprepro not found");
	       }

#else
             EH(-1, "aprepro the mat file prior to running goma.");
#endif
	     if ( Debug_Flag > 0 )
	       {
		 fprintf(stdout, "system: %s\n", System_Command);
	       }
	     imp = fopen(TmpMatFile, "r");
	     if ( imp == NULL )
	       {
		 EH(-1, "Problem opening temporary mat file.");	       }
	   }

         /*
	  * Now read the material properties from the correct file
	  */
	 rd_mp_specs(imp, input, mn, echo_mat_file);
	 fclose(imp);
       } else {
          EH(-1,"Not all Material Files found in current directory");
       }
       /*
	*  Print a description of the materials properties for the
	*  current material
	*/
	   if( Debug_Flag > 2 )
	   {
       (void) matrl_prop_print(mp_glob[mn], mn);
	   }
       ECHO("CLOSE",echo_mat_file);
    }  /* END OF for(mn=0; mn<upd->Num_Mat... */

  ECHO("\nEND OF MAT\n", echo_input_file);

  if (upd->Total_Num_EQ > MAX_PROB_VAR) {
    fprintf(stderr,
	    "%sTotal number of equations, %d, exceeds MAX_PROB_VAR, %d\n",
	    "input FATAL ERROR: ", upd->Total_Num_EQ, MAX_PROB_VAR);
    fprintf(stderr,
	    "\tIncrease MAX_PROB_VAR to %d and rerun\n",
	    upd->Total_Num_EQ);
      EH(-1, "MAX_PROB_VAR is too small for number of equations");
  }

}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

static int set_eqn(int eqnNum, PROBLEM_DESCRIPTION_STRUCT *pd_ptr)

     /******************************************************************
      *
      * set_eqn()
      *
      *   Utility function used repeatedly as a kernel algorithm
      *   in rd_eq_specs.
      ******************************************************************/
{
  if (pd_ptr->e[eqnNum]) {
    fprintf(stderr, 
	    "rd_eq_spec ERROR: eqnNum %d has already been activated.\n",
	    eqnNum);
    EH(-1,"rd_eq_spec ERROR fatal error in equation specification\n");
  } else {
    pd_ptr->e[eqnNum] = T_SOMETHING;
  }
  return eqnNum;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

static int set_var(int varType, PROBLEM_DESCRIPTION_STRUCT *pd_ptr)

     /******************************************************************
      *
      * set_var()
      *
      *   Utility function used repeatedly as a kernel algorithm
      *   in rd_eq_specs.
      ******************************************************************/
{
  if (pd_ptr->v[varType]) {
    fprintf(stderr, 
	    "rd_eq_spec ERROR: variable %d has already been activated.\n",
	    varType);
    EH(-1,"rd_eq_spec ERROR fatal error in variable specification\n");
  }  else {
    pd_ptr->v[varType] |= V_SOLNVECTOR;
  }
  return varType;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * rd_eq_specs -- read input file for equation & term specifications
 *
 * Comments:	This was written from scratch as an attempt to avoid
 *		permutation madness.
 *
 * Created:			Fri Oct 29 14:02:55 MDT 1993 pasacki@sandia.gov
 * Revised:			Wed Dec  1 15:24:26 MST 1993 pasacki@sandia.gov
 * Revised:			Thr Oct 20 15:24:26 MST 1994 racairn@sandia.gov 
 *					(adding species)
 * Revised:			Thr Feb 16 15:24:26 MST 1995 prschun@sandia.gov 
 *					(Multiple materials)
 */

void 
rd_eq_specs(FILE *ifp,
	    char *input,
	    const int mn)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  int numBulkSpecEqn;		/* temporary variable only!! */
  int numBulkSpec;              /* temporary variable only!! */
  int nonDilute = FALSE;        /* Flag to indicate if the species
				   equations use the dilute approximation or not */
  int	i, j, retn;
  int	ieqn;			/* Current equation linefrom input. */
  int   neqn;                   /* Current number of equations active within
				   the current material. This will be greater
				   or equal to the number of equation lines */
  int	ce;			/* Current identifiable equation. */
  int	cv = 0;			/* Current identifiable variable. */
  int cem, cvm, n, err = 0;
  int eqnMult;                  /* This variable refers to how many actual variable
				 * types the current line refers to. Normally, its
				 * one. However, for the SPECIES_UNK_# variables
				 * types, it can be greater than one
				 */
  int numEqnLines = 0;          /* Number of Equation lines to be read */
  int MeshMotion = -1;
  PROBLEM_DESCRIPTION_STRUCT *pd_ptr = pd_glob[mn];
  char	ts[MAX_BC_KEYWORD_LENGTH]   = "\0";

  char	tscs[MAX_CS_KEYWORD_LENGTH] = "\0";

  char echo_string[MAX_CHAR_ECHO_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  static char yo[] = "rd_eq_specs";

  /*
   * Identify section containing equation specification...
   */

  look_for(ifp, "Coordinate System", input, '=');

  if (fscanf(ifp, "%s", tscs) != 1)
    {
      EH( -1, "error reading Coordinate System");
    }
  if ( !strcasecmp(tscs, "CARTESIAN") )
    {
      CoordinateSystem = CARTESIAN;
    }
  else if ( !strcasecmp(tscs, "CYLINDRICAL") )
    {
      CoordinateSystem = CYLINDRICAL;
    }
  else if ( !strcasecmp(tscs, "SPHERICAL") )
    {
      CoordinateSystem = SPHERICAL;
    }
  else if ( !strcasecmp(tscs, "SWIRLING") )
    {
      CoordinateSystem = SWIRLING;
    }   
  else if ( !strcasecmp(tscs, "PROJECTED_CARTESIAN") )
    {
      CoordinateSystem = PROJECTED_CARTESIAN;
    }
  else if ( !strcasecmp(tscs, "CARTESIAN_2pt5D") )
    {
      CoordinateSystem = CARTESIAN_2pt5D;
    }
  else
    {
      EH(-1,"Unknown coordinate system.\n");
    }

  pd_ptr->CoordinateSystem = CoordinateSystem;
  upd->CoordinateSystem = CoordinateSystem;

  SPF(echo_string,"%s = %s","Coordinate System", tscs); ECHO(echo_string,echo_file);

  /*
   *  Element Mapping line
   */
  
  look_for(ifp,"Element Mapping",input,'=');
  (void) read_string(ifp,input,'\n');
  strip(input);
  if ( strcasecmp(input,"isoparametric") == 0 )
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  printf("%s:\t(T,c) = isoparametric\n", yo);
	}
	pd_ptr->IntegrationMap  = ISOPARAMETRIC;
    }
  else if ( strcasecmp(input,"subparametric") == 0 )
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  printf("%s:\t(T,c) = subparametric\n", yo);
	}
      fprintf (stderr, "%s:\tSubparametric mapping not active yet\n", yo);
      exit (-1);
	pd_ptr->IntegrationMap  = SUBPARAMETRIC;
    }
  else if (  strcasecmp(input,"Q1") == 0 )
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  printf("%s:\t(T,c) = Q1\n", yo);
	}
	pd_ptr->IntegrationMap  = I_Q1;
    }
  else if (  strcasecmp(input,"Q2") == 0 )
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  printf("%s:\t(T,c) = Q2\n", yo);
	}
	pd_ptr->IntegrationMap  = I_Q2;
    }      
  else if (  strcasecmp(input,"SP") == 0 )
    {
      if ( Debug_Flag && ProcID == 0 )
	{
	  printf("%s:\t(T,c) = SP\n", yo);
	}
	pd_ptr->IntegrationMap  = I_SP;
        ECHO("\n\t(CAUTION WITH SUBPARAMETRIC: DON'T USE WITH LAGRANGIAN OR TOTAL_ALE MESH MOTIONS.  Also not good for GD_TABLE conditions.)", echo_file);
    }      
  else
    {
      EH( -1, "invalid element mapping");
    }
  SPF(echo_string,"%s = %s", "Element Mapping", input); ECHO(echo_string, echo_file);

  /*
   *  Mesh Motion Line
   */
  
  look_for(ifp, "Mesh Motion", input, '=');
  if (fscanf(ifp, "%s", tscs) != 1)
  {
    EH( -1, "error reading Mesh Motion");
  }
  if ( !strcmp(tscs, "ARBITRARY") )
    {
      MeshMotion = ARBITRARY;
    }
  else if ( !strcmp(tscs, "LAGRANGIAN") )
    {
      MeshMotion = LAGRANGIAN;
    }   
  else if ( !strcmp(tscs, "DYNAMIC_LAGRANGIAN") )
    {
      MeshMotion = DYNAMIC_LAGRANGIAN;
    }   
  else if ( !strcmp(tscs, "TOTAL_ALE") )
    {
      MeshMotion = TOTAL_ALE;
    }
  else {
    fprintf(stderr, "%s: unrecognized mesh motion type: %s\n",
	    yo, tscs);
    exit(-1);
  }

  SPF(echo_string,"%s = %s","Mesh Motion", tscs); ECHO(echo_string,echo_file);

  pd_ptr->MeshMotion = MeshMotion;
  cr_glob[mn]->MeshMotion = pd_ptr->MeshMotion;


  /*
   *  Number of Bulk Species Line
   *
   * Get number of species for bulk convection-diffusion equation in
   * this material
   *      - this is a required parameter
   */
  look_for(ifp, "Number of bulk species", input, '=');
  numBulkSpec = read_int(ifp, "Number of bulk species");
  pd_ptr->Num_Species_Eqn = numBulkSpec;
  pd_ptr->Num_Species     = numBulkSpec;
  pd_ptr->Num_Porous_Eqn  = 0;

  SPF(echo_string,"%s = %d", "Number of bulk species", numBulkSpec); ECHO(echo_string, echo_file);

  /*
   * Backwards compatibility Line
   *
   */
  retn = look_forward_optional(ifp,
              "Really meant number of bulk species equations in line above",
			       input, '=');
  if (retn == 1) {
    retn = read_1_boolean(ifp,"Really meant number of bulk species equations in line above" );
    ECHO("Really meant number of bulk species equations in line above = TRUE", echo_file);
  } else {
    retn = FALSE;
  }
  if (retn) {
    numBulkSpecEqn = numBulkSpec;
    numBulkSpec = numBulkSpecEqn + 1;
    pd_ptr->Num_Species     = numBulkSpecEqn + 1;
  } else {
    numBulkSpecEqn = numBulkSpec;
  }

  /*
   *  The default is to assume we are solving a dilute
   *  system where the last species in the
   *  mechanism is in such excess that its continuity equation is taken
   *  care of via the main mass continuity equation (or by some other specified
   *  means). If the total number of species is one greater than the
   *  total number of species equations and the system is nondilute, then
   *  the total number of species equations is one less than the total
   *  number of species.
   *
   *  The non-dilute case will be indicated by setting Num_Species_Eqn
   *  to one less than Num_Species. In the input file there are two
   *  options to specify non-dilute: Either can be used. However, if
   *  they lead to contradictory input, the code will error exit.
   */
  
  retn = look_forward_optional(ifp, "Material is nondilute", input, '=');
  if (retn == 1) {
    nonDilute = read_1_boolean(ifp, "Material is nondilute"); 
    ECHO("Material is nondilute = TRUE", echo_file);
  }
  if (nonDilute) {
    pd_ptr->Num_Species_Eqn = numBulkSpec - 1;
    if (pd_ptr->Num_Species_Eqn < 0) pd_ptr->Num_Species_Eqn = 0;
  } 
  retn = look_forward_optional(ifp, "Number of bulk species equations", input, '=');
  if (retn == 1) {
    numBulkSpecEqn  = read_int(ifp, "Number of bulk species equations");
    if (nonDilute) {
      if (numBulkSpecEqn !=  pd_ptr->Num_Species_Eqn) {
	fprintf(stderr,"Confusing input: number of bulk species equations is illdefined\n");
	EH( -1, "\tNumber of bulk species equations is inconsistent with non-dilute assumption\n");
      }
      SPF(echo_string,"%s = %d","Number of bulk species equations", numBulkSpecEqn);
    }
    pd_ptr->Num_Species_Eqn =  numBulkSpecEqn;
  }
  if (pd_ptr->Num_Species > MAX_CONC) {
    sprintf(err_msg,
	    "Specified no. of species %d >= %d, must boost MAX_CONC to %d by using"
	    " either -DMAX_CONC option in Goma.mk or change it in rf_fem_const.h.",
	    pd_ptr->Num_Species, MAX_CONC, pd_ptr->Num_Species);
    EH(-1, err_msg);
    ABORTH(-1, "IMMEDIATE PARALLEL EXIT - can't recover gracefully");
  }

  /*
   * Now store the maximum values for the number of species equations and
   * number of species in the upd structure
   */
  upd->Max_Num_Species_Eqn = MAX(upd->Max_Num_Species_Eqn, pd_ptr->Num_Species_Eqn);
  upd->Max_Num_Species     = MAX(upd->Max_Num_Species    , pd_ptr->Num_Species);
 
  /*
   *  Write the resulting number of species equations out to the output file
   */
  if (pd_ptr->Num_Species_Eqn > 0)  {
    SPF(echo_string,"\t(Number of bulk C/D equations  = %d)", pd_ptr->Num_Species_Eqn); ECHO(echo_string,echo_file);
    if (pd_ptr->Num_Species_Eqn != pd_ptr->Num_Species) {
      SPF(echo_string, "\t(Nondilute approximation: # species is not equal to number of species equations)\n");
      SPF(endofstring(echo_string),"\t(Number of Species = %d)", pd_ptr->Num_Species); 
      ECHO(echo_string,echo_file);
   }
  }

  /* 
   * Currently, all materials must have either 0 or maximum number of species
   *
   * If this material wants more species than any previous material, update
   * the "Max_Num_Spec" variable for all materials up to and including this
   * one.
   */
  if (pd_ptr->Num_Species_Eqn > upd->Max_Num_Species_Eqn ) {
    upd->Max_Num_Species_Eqn = pd_ptr->Num_Species_Eqn;
    upd->Max_Num_Species     = pd_ptr->Num_Species;
  }

  /*
   * Specify the types of unknowns for the independent variables representing
   * the species unknowns. Currently, this is a placeholder that always
   * indicates SPECIES_UNDEFINED_FORM. The code may actually pay attention to
   * this value later.
   */
  retn = look_for_next_string(ifp, "Default Material Species Type", input, '=');
  if (retn == 1) {
    retn = read_string(ifp, input, '\n');
    strip(input);
    pd_ptr->Species_Var_Type = species_type_str_to_int(input);
    if (upd->Species_Var_Type == 0) upd->Species_Var_Type = pd_ptr->Species_Var_Type;
    SPF(echo_string,"%s = %s", "Default Material Species Type", input); ECHO(echo_string,echo_file);
  }
  
  /*
   *  Search for the number of viscoelastic modes
   */
  err = look_forward_optional(ifp, "Number of viscoelastic modes", input, '=');

  if (err == -1)
  {
	  vn_glob[mn]->modes = 0;
	  SPF(echo_string,"%s = %d", "Number of viscoelastic modes", vn_glob[mn]->modes); ECHO(echo_string,echo_file);
  }
  else
  {
    if ( fscanf(ifp,"%d",&vn_glob[mn]->modes) != 1)
    {
      EH( -1, "Expected to read 1 int for \"Number of viscoelastic modes\"");
    }

    SPF(echo_string,"%s = %d", "Number of viscoelastic modes", vn_glob[mn]->modes); ECHO(echo_string,echo_file);
  } 
  
  /*
   * User has to say how many equations will be specified; just like
   * boundary conditions...
   */

  look_for(ifp, "Number of EQ", input, '='); 

  if (fscanf(ifp, "%d", &numEqnLines) != 1) {
    fprintf(stderr, 
	    "Expected to read 1 int for \"Number of EQ\" - attempt to enumerate...");
    numEqnLines = -1;
  }


  /*
   * Count Number of Equations if Num_EQ is set to -1
   */
  if (numEqnLines == -1) {
    numEqnLines = count_list(ifp, "EQ", input, '=', "END OF EQ");
  }

  SPF(echo_string,"%s = %d", "Number of EQ", numEqnLines); ECHO(echo_string,echo_file);

  if ( (numEqnLines < 1 || numEqnLines > 25)) {
    EH( -1, "rd_eq_spec ERROR: Too many (>25) or too few eqn. lines");
  }

  /*
   * Here's the low down...
   *
   * [1] There are Num_EQ equations that are active. These will be a
   *     distinct list out the permissible equation types.
   * 
   * [2] There will be a "variable" unknown that is generally associated
   *     with each equation. For example the variable "T", or temperature
   *     is generally associated with the "energy" equation.
   * 
   * [3] Each variable will be interpolated using a particular kind of
   *     basis function.
   * 
   * [4] The coordinate system type, the number of spatial dimensions, and
   *     the time integration technique are loaded from the variables of the
   *     same name, many of which are read into the EXODUS II file.
   */

  /*
   * Initialize....
   */

  ce = 0;			/* Current equation. */

  for ( i=0; i<MAX_EQNS; i++ )
  {
    pd_ptr->e[i] = T_NOTHING;	/* No terms active in this equation. */
    pd_ptr->v[i] = V_NOTHING;	/* No active variables, either. */
    pd_ptr->w[i] = I_NOTHING;	/* No weighting function. */
    pd_ptr->i[i] = I_NOTHING;	/* Nothing to interpolate. */
    pd_ptr->m[i] = -1;		/* Map from i to R_EQNTYPE initially undefd. */
    for ( j=0; j<MAX_TERM_TYPES; j++)
    {
      pd_ptr->etm[i][j] = 0.;	/* Zero out all terms. */
    }
  }

  /*
   * Read in each equation and set up entries in pd...
   */
  neqn = 0;
  for (ieqn = 0; ieqn < numEqnLines; ieqn++) { /*ieqn*/
    /*
     * reset eqnMult
     *   This variable refers to how many actual variable
     * types the current line refers to. Normally, its one. However,
     * for the SPECIES_UNK_# variables types, it can be greater than one
     */
    eqnMult = 1;

    /*
     * Look for each EQ line and then parse the results
     */
     
    look_for(ifp, "EQ", input, '=');

    /* 
     *
     * Part 1 -- Look at the type of equation "energy", "mesh1", etc.
     * 
     */

    if (fscanf(ifp, "%s", ts) != 1)
    {
      fprintf (stderr, "%s:\tError reading EQ %d\n", yo, ieqn);
      exit (-1);
    }

    /*
     * Check to see if this equation type is recognizable...
     */

    if        (!strcasecmp(ts, "energy")) {
      ce = set_eqn(R_ENERGY, pd_ptr);
    } else if (!strcasecmp(ts, "momentum1")) {
      ce = set_eqn(R_MOMENTUM1, pd_ptr);
    } else if (!strcasecmp(ts, "momentum2")) {
      ce = set_eqn(R_MOMENTUM2, pd_ptr);
    } else if (!strcasecmp(ts, "momentum3")) {
      ce = set_eqn(R_MOMENTUM3, pd_ptr);
    } else if (!strcasecmp(ts, "pmomentum1")) {
      ce = set_eqn(R_PMOMENTUM1, pd_ptr);
    } else if (!strcasecmp(ts, "pmomentum2")) {
      ce = set_eqn(R_PMOMENTUM2, pd_ptr);
    } else if (!strcasecmp(ts, "pmomentum3")) {
      ce = set_eqn(R_PMOMENTUM3, pd_ptr);
    } else if (!strcasecmp(ts, "stress11")) {
      ce = set_eqn(R_STRESS11, pd_ptr);
      Use_DG = TRUE;
    } else if (!strcasecmp(ts, "stress12")) {
      ce = set_eqn(R_STRESS12, pd_ptr);
    } else if (!strcasecmp(ts, "stress13")) {
      ce = set_eqn(R_STRESS13, pd_ptr);
    } else if (!strcasecmp(ts, "stress22")) {
      ce = set_eqn(R_STRESS22, pd_ptr);
    } else if (!strcasecmp(ts, "stress23")) {
      ce = set_eqn(R_STRESS23, pd_ptr);
    } else if (!strcasecmp(ts, "stress33")) {
      ce = set_eqn(R_STRESS33, pd_ptr);
    } else if (!strcasecmp(ts, "gradient11")) {
      ce = set_eqn(R_GRADIENT11, pd_ptr);
    } else if (!strcasecmp(ts, "gradient12")) {
      ce = set_eqn(R_GRADIENT12, pd_ptr);
    } else if (!strcasecmp(ts, "gradient13")) {
      ce = set_eqn(R_GRADIENT13, pd_ptr);
    } else if (!strcasecmp(ts, "gradient21")) {
      ce = set_eqn(R_GRADIENT21, pd_ptr);
    } else if (!strcasecmp(ts, "gradient22")) {
      ce = set_eqn(R_GRADIENT22, pd_ptr);
    } else if (!strcasecmp(ts, "gradient23")) {
      ce = set_eqn(R_GRADIENT23, pd_ptr);
    } else if (!strcasecmp(ts, "gradient31")) {
      ce = set_eqn(R_GRADIENT31, pd_ptr);
    } else if (!strcasecmp(ts, "gradient32")) {
      ce = set_eqn(R_GRADIENT32, pd_ptr);
    } else if (!strcasecmp(ts, "gradient33")) {
      ce = set_eqn(R_GRADIENT33, pd_ptr);
    } else if (!strcasecmp(ts, "bond")) {
      ce = set_eqn(R_BOND_EVOLUTION, pd_ptr);
    } else if (!strcasecmp(ts, "species_bulk")) {
      ce = set_eqn(R_MASS, pd_ptr);
#ifdef DEBUG
      printf("Input Number of Species Equations =  %d\n",
	     pd_ptr->Num_Species_Eqn);
#endif
      if (pd_ptr->Num_Species_Eqn == 0) {
	EH(-1, "Cannot solve species transport without number of species");
      }
    } else if (!strcasecmp(ts, "mesh1")) {
      ce = set_eqn(R_MESH1, pd_ptr);
    } else if (!strcasecmp(ts, "mesh2")) {
      ce = set_eqn(R_MESH2, pd_ptr);
    } else if (!strcasecmp(ts, "mesh3")) {
      ce = set_eqn(R_MESH3, pd_ptr);
    } else if (!strcasecmp(ts, "mom_solid1")) {
      ce = set_eqn(R_SOLID1, pd_ptr);
    } else if (!strcasecmp(ts, "mom_solid2")) {
      ce = set_eqn(R_SOLID2, pd_ptr);
    } else if (!strcasecmp(ts, "mom_solid3")) {
      ce = set_eqn(R_SOLID3, pd_ptr);
    } else if (!strcasecmp(ts, "species_surf")) {
      ce = set_eqn(R_MASS_SURF, pd_ptr);
    } else if (!strcasecmp(ts, "continuity")) {
      ce = set_eqn(R_PRESSURE, pd_ptr);
    } else if (!strcasecmp(ts, "continuity_em_real")) {
      ce = set_eqn(R_EM_CONT_REAL, pd_ptr);
    } else if (!strcasecmp(ts, "continuity_em_imag")) {
      ce = set_eqn(R_EM_CONT_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "fill")) {
      ce = set_eqn(R_FILL, pd_ptr);
#ifndef COUPLED_FILL
      /* set explicit flag for fill equation */
      Explicit_Fill = 1;
#endif /* not COUPLED_FILL */
    } else if (!strcasecmp(ts, "level_set"))  {
      if ( ls == NULL ) EH(-1,"Level Set Interface tracking must be turned on.");
      ce = set_eqn(R_LEVEL_SET, pd_ptr);
#ifndef COUPLED_FILL
      /* set explicit flag for fill equation and turn on level set switch*/
      if (ls->Evolution == LS_EVOLVE_ADVECT_EXPLICIT  ||
          ls->Evolution == LS_EVOLVE_SEMILAGRANGIAN ) Explicit_Fill = 1;
#endif /* not COUPLED_FILL */
    } else if (!strcasecmp(ts, "curvature") ) {
      ce = set_eqn( R_CURVATURE, pd_ptr);
    } else if (!strcasecmp(ts,"normal1")) {
      if ( ls == NULL )  EH(-1," NORMAL1 equation requires Level Set Interface tracking on.");
      ce = set_eqn( R_NORMAL1, pd_ptr );
    } else if (!strcasecmp(ts,"normal2")) {
      if ( ls == NULL )  EH(-1," NORMAL2 equation requires Level Set Interface tracking on.");
      ce = set_eqn( R_NORMAL2, pd_ptr );
    } else if (!strcasecmp(ts,"normal3")) {
      if ( ls == NULL )  EH(-1," NORMAL3 equation requires Level Set Interface tracking on.");
      ce = set_eqn( R_NORMAL3, pd_ptr );
    } else if (!strcasecmp(ts, "voltage")) {
      ce = set_eqn(R_POTENTIAL, pd_ptr);
    } else if (!strcasecmp(ts, "surf_charge")) {
      ce = set_eqn(R_SURF_CHARGE, pd_ptr);
    } else if (!strcasecmp(ts, "shear_rate"))  {
      ce = set_eqn(R_SHEAR_RATE, pd_ptr);
    } else if (!strcasecmp(ts, "vort_dir1")) {
      ce = set_eqn(R_VORT_DIR1, pd_ptr);
    } else if (!strcasecmp(ts, "vort_dir2")) {
      ce = set_eqn(R_VORT_DIR2, pd_ptr);
    } else if (!strcasecmp(ts, "vort_dir3")) {
      ce = set_eqn(R_VORT_DIR3, pd_ptr);
    } else if (!strcasecmp(ts, "vort_lambda")) {
      ce = set_eqn(R_VORT_LAMBDA, pd_ptr);
    } else if (!strcasecmp(ts, "lagr_mult_1")) {
      ce = set_eqn(R_LAGR_MULT1, pd_ptr);
    } else if (!strcasecmp(ts, "lagr_mult_2")) {
      ce = set_eqn(R_LAGR_MULT2, pd_ptr);
    } else if (!strcasecmp(ts, "lagr_mult_3")) {
      ce = set_eqn(R_LAGR_MULT3, pd_ptr);
    } else if (!strcasecmp(ts, "shell_curvature")) {
      ce = set_eqn(R_SHELL_CURVATURE, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_curvature2")) {
      ce = set_eqn(R_SHELL_CURVATURE2, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_tension")) {
      ce = set_eqn(R_SHELL_TENSION, pd_ptr);
    } else if (!strcasecmp(ts, "shell_x")) {
      pd_ptr->Do_Surf_Geometry = 1;
      ce = set_eqn(R_SHELL_X, pd_ptr);
    } else if (!strcasecmp(ts, "shell_y")) {
      ce = set_eqn(R_SHELL_Y, pd_ptr);
    } else if (!strcasecmp(ts, "shell_user")) {
      ce = set_eqn(R_SHELL_USER, pd_ptr);
    } else if (!strcasecmp(ts, "shell_angle1")) {
      ce = set_eqn(R_SHELL_ANGLE1, pd_ptr);
    } else if (!strcasecmp(ts, "shell_angle2")) {
      ce = set_eqn(R_SHELL_ANGLE2, pd_ptr);
    } else if (!strcasecmp(ts, "shell_surf_div_v")) {
      ce = set_eqn(R_SHELL_SURF_DIV_V, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_surf_curv")) {
      ce = set_eqn(R_SHELL_SURF_CURV, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "n_dot_curl_v")) {
      ce = set_eqn(R_N_DOT_CURL_V, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "grad_v_dot_n1")) {
      ce = set_eqn(R_GRAD_S_V_DOT_N1, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "grad_v_dot_n2")) {
      ce = set_eqn(R_GRAD_S_V_DOT_N2, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "grad_v_dot_n3")) {
      ce = set_eqn(R_GRAD_S_V_DOT_N3, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "acous_preal")) {
      ce = set_eqn(R_ACOUS_PREAL, pd_ptr);
    } else if (!strcasecmp(ts, "acous_pimag")) {
      ce = set_eqn(R_ACOUS_PIMAG, pd_ptr);
    } else if (!strcasecmp(ts, "acous_reyn_stress")) {
      ce = set_eqn(R_ACOUS_REYN_STRESS, pd_ptr);
    } else if (!strcasecmp(ts, "shell_bdyvelo")) {
      ce = set_eqn(R_SHELL_BDYVELO, pd_ptr);
    } else if (!strcasecmp(ts, "shell_lubp")) {
      ce = set_eqn(R_SHELL_LUBP, pd_ptr);
    } else if (!strcasecmp(ts, "lubp")) {
      ce = set_eqn(R_LUBP, pd_ptr);
   } else if (!strcasecmp(ts, "lubp_2")) {
      ce = set_eqn(R_LUBP_2, pd_ptr);
    } else if (!strcasecmp(ts, "shell_filmp")) {
      ce = set_eqn(R_SHELL_FILMP, pd_ptr);
    } else if (!strcasecmp(ts, "shell_filmh")) {
      ce = set_eqn(R_SHELL_FILMH, pd_ptr);
    } else if (!strcasecmp(ts, "shell_partc")) {
      ce = set_eqn(R_SHELL_PARTC, pd_ptr);
    } else if (!strcasecmp(ts, "shell_sat_closed")) {
      ce = set_eqn(R_SHELL_SAT_CLOSED, pd_ptr);
    } else if (!strcasecmp(ts, "shell_sat_open")) {
      ce = set_eqn(R_SHELL_SAT_OPEN, pd_ptr);
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "shell_sat_open_2")) {
      ce = set_eqn(R_SHELL_SAT_OPEN_2, pd_ptr);
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "shell_energy")) {
      ce = set_eqn(R_SHELL_ENERGY, pd_ptr);
    } else if (!strcasecmp(ts, "shell_deltah")) {
      ce = set_eqn(R_SHELL_DELTAH, pd_ptr);
    } else if (!strcasecmp(ts, "shell_lub_curv")) {
      ce = set_eqn(R_SHELL_LUB_CURV, pd_ptr);
    } else if (!strcasecmp(ts, "shell_lub_curv_2")) {
      ce = set_eqn(R_SHELL_LUB_CURV_2, pd_ptr);
    } else if (!strcasecmp(ts, "shell_sat_gasn")) {
      ce = set_eqn(R_SHELL_SAT_GASN, pd_ptr);
    } else if (!strcasecmp(ts, "sink_mass")) {
      ce = set_eqn(R_POR_SINK_MASS, pd_ptr);
    } else if (!strcasecmp(ts, "shell_shear_top")) {
      ce = set_eqn(R_SHELL_SHEAR_TOP, pd_ptr);
    } else if (!strcasecmp(ts, "shell_shear_bot")) {
      ce = set_eqn(R_SHELL_SHEAR_BOT, pd_ptr);
    } else if (!strcasecmp(ts, "shell_cross_shear")) {
      ce = set_eqn(R_SHELL_CROSS_SHEAR, pd_ptr);
    } else if (!strcasecmp(ts, "max_strain")) {
      ce = set_eqn(R_MAX_STRAIN, pd_ptr);
    } else if (!strcasecmp(ts, "cur_strain")) {
      ce = set_eqn(R_CUR_STRAIN, pd_ptr);
    } else if (!strcasecmp(ts, "shell_diff_flux")) {
      ce = set_eqn(R_SHELL_DIFF_FLUX, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_diff_curv")) {
      ce = set_eqn(R_SHELL_DIFF_CURVATURE, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_normal1")) {
      ce = set_eqn(R_SHELL_NORMAL1, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_normal2")) {
      ce = set_eqn(R_SHELL_NORMAL2, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "shell_normal3")) {
      ce = set_eqn(R_SHELL_NORMAL3, pd_ptr);
      pd_ptr->Do_Surf_Geometry = 1;
    } else if (!strcasecmp(ts, "ext_v")) {
      ce = set_eqn(R_EXT_VELOCITY, pd_ptr);
      ls->Extension_Velocity = TRUE;
    } else if (!strcasecmp(ts, "efield1")) {
      ce = set_eqn(R_EFIELD1, pd_ptr);
    } else if (!strcasecmp(ts, "efield2")) {
      ce = set_eqn(R_EFIELD2, pd_ptr);
    } else if (!strcasecmp(ts, "efield3")) {
      ce = set_eqn(R_EFIELD3, pd_ptr);
    } else if (!strcasecmp(ts, "intp")) {
      ce = set_eqn(R_LIGHT_INTP, pd_ptr);
    } else if (!strcasecmp(ts, "intm")) {
      ce = set_eqn(R_LIGHT_INTM, pd_ptr);
    } else if (!strcasecmp(ts, "intd")) {
      ce = set_eqn(R_LIGHT_INTD, pd_ptr);
    } else if (!strcasecmp(ts, "phase1")) {
      ce = set_eqn(R_PHASE1, pd_ptr);
    } else if (!strcasecmp(ts, "phase2")) {
      ce = set_eqn(R_PHASE2, pd_ptr);
    } else if (!strcasecmp(ts, "phase3")) {
      ce = set_eqn(R_PHASE3, pd_ptr);
    } else if (!strcasecmp(ts, "phase4")) {
      ce = set_eqn(R_PHASE4, pd_ptr);
    } else if (!strcasecmp(ts, "phase5")) {
      ce = set_eqn(R_PHASE5, pd_ptr);
   } else if (!strcasecmp(ts, "Enorm"))  {
      ce = set_eqn(R_ENORM, pd_ptr);
   } else if (!strcasecmp(ts, "tfmp_mass"))  {
      ce = set_eqn(R_TFMP_MASS, pd_ptr);
   } else if (!strcasecmp(ts, "tfmp_bound"))  {
      ce = set_eqn(R_TFMP_BOUND, pd_ptr);
    } else if (!strcasecmp(ts, "restime")) {
      ce = set_eqn(R_RESTIME, pd_ptr);  
    } else if (!strcasecmp(ts, "em_e1_real")) {
      ce = set_eqn(R_EM_E1_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_e2_real")) {
      ce = set_eqn(R_EM_E2_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_e3_real")) {
      ce = set_eqn(R_EM_E3_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_e1_imag")) {
      ce = set_eqn(R_EM_E1_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "em_e2_imag")) {
      ce = set_eqn(R_EM_E2_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "em_e3_imag")) {
      ce = set_eqn(R_EM_E3_IMAG, pd_ptr);  
    } else if (!strcasecmp(ts, "em_h1_real")) {
      ce = set_eqn(R_EM_H1_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_h2_real")) {
      ce = set_eqn(R_EM_H2_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_h3_real")) {
      ce = set_eqn(R_EM_H3_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "em_h1_imag")) {
      ce = set_eqn(R_EM_H1_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "em_h2_imag")) {
      ce = set_eqn(R_EM_H2_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "em_h3_imag")) {
      ce = set_eqn(R_EM_H3_IMAG, pd_ptr);  

    } else if (!strcasecmp(ts, "porous_sat"))  {
      ce = set_eqn(R_POR_LIQ_PRES, pd_ptr);
      if( pd_ptr->e[R_POR_GAS_PRES] || pd_ptr->e[R_POR_SATURATION] 
          || pd_ptr->e[R_POR_ENERGY] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_unsat")) {
      ce = set_eqn(R_POR_LIQ_PRES, pd_ptr);
      if( pd_ptr->e[R_POR_GAS_PRES] || pd_ptr->e[R_POR_SATURATION] 
          || pd_ptr->e[R_POR_ENERGY] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_liq"))  {
      ce = set_eqn(R_POR_LIQ_PRES, pd_ptr);
      if( tran->theta != 0.0)EH(-1,"Please don't use Time step parameter of 0.5                                   for porous media problems.  Use 0.0");
      if( pd_ptr->e[R_POR_SATURATION] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_gas"))  {
      ce = set_eqn(R_POR_GAS_PRES, pd_ptr);
      if( pd_ptr->e[R_POR_SATURATION] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_energy"))  {
      ce = set_eqn(R_POR_ENERGY, pd_ptr);
      if( pd_ptr->e[R_POR_SATURATION] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_deform"))  {
      ce = set_eqn(R_POR_POROSITY, pd_ptr);
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strcasecmp(ts, "porous_saturation"))  {
      ce = set_eqn(R_POR_SATURATION, pd_ptr);
      if( pd_ptr->e[R_POR_GAS_PRES] || pd_ptr->e[R_POR_LIQ_PRES] 
          || pd_ptr->e[R_POR_ENERGY] )
        {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "More than one porous media equation type");
          EH(-1,"Too many porous media equations for this material block");
        }
      pd_ptr->Num_Porous_Eqn++;
    } else if (!strncasecmp(ts, "species_unk", 11)) {
      if (!strcasecmp(ts, "species_unk")) {
	ce = R_SPECIES_UNK_0;
#ifdef DEBUG
	printf("Input Number of Species Equations =  %d\n",
	       pd_ptr->Num_Species_Eqn);
#endif
	if (pd_ptr->Num_Species_Eqn == 0) {
	  fprintf(stderr, "%s%s\n",
		   "rd_eq_specs ERROR: ",
		   "Number of species unknowns is equal to zero");

          EH(-1,"Num_species is zero yet species equation has been selected");
	}
	j = R_SPECIES_UNK_0;
	eqnMult = pd_ptr->Num_Species_Eqn;
	for (i = 0; i < pd_ptr->Num_Species_Eqn; i++, j++) {
	  (void) set_eqn(j, pd_ptr);
	  pd_ptr->m[j] = R_SPECIES_UNK_0 + j;
	}
      } else {
	if (!interpret_int(ts + 12, &j)) {
	  fprintf(stderr,
		  "%s: Unable to decipher species eqn labeled: %s\n",
		  yo, ts);
	  exit(-1);
	}
	if (j < 0 || j > MAX_SPECIES_UNK_NUM) {
	  fprintf(stderr,
		  "%s: Species eqn %d is out of the permissible range\n",
		  yo, j);
	  exit(-1);
	}
	ce = set_eqn(R_SPECIES_UNK_0 + j, pd_ptr);
      }
    } else {
      fprintf(stderr, "%s:\tEQ %s not recognized.\n", yo, ts);
      exit(-1);
    }


    SPF(echo_string,"\t\t%s = %s","EQ", ts);

   upd->Max_Num_Porous_Eqn = MAX(upd->Max_Num_Porous_Eqn, pd_ptr->Num_Porous_Eqn);

    /*
     * This sets up a useful mapping over the active equations...
     */
    for (j = 0; j < eqnMult; j++) {
      pd_ptr->m[neqn] = ce + j;
      neqn++;
    }

    /*
     * Setup actual equation to problem equation index arrays
     * The ep[] array assigns a unique index to every active
     * equation in the problem. This index is the same for all
     * materials in the problem. 
     */
    for (j = 0; j < eqnMult; j++) {
      if (upd->ep[ce + j] == -1) {
	upd->ep[ce + j] = (upd->Total_Num_EQ)++;
      }
    }

    /*
     *
     * Part 2 -- Look at the Galerkin weight for this equation. This might
     *           be something like "Q1", "P0", etc.
     *
     */

    if (fscanf(ifp, "%s", ts) !=1) {
      fprintf(stderr, "%s: problem reading Galerkin weight for eqn.\n",yo);

    }
    /*
     * Check to see if this Galerkin weight is recognizable and if it
     * is reasonable.
     */

    if (!strcasecmp(ts, "Q1"))
    {
      pd_ptr->w[ce] = I_Q1;
    }
    else if (!strcasecmp(ts, "Q2"))
    {
      pd_ptr->w[ce] = I_Q2;
    }
    else if (!strcasecmp(ts, "Q1_G"))
    {
      pd_ptr->w[ce] = I_Q1_G;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_G"))
    {
      pd_ptr->w[ce] = I_Q2_G;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_GP"))
    {
      pd_ptr->w[ce] = I_Q1_GP;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_GP"))
    {
      pd_ptr->w[ce] = I_Q2_GP;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_GN"))
    {
      pd_ptr->w[ce] = I_Q1_GN;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_GN"))
    {
      pd_ptr->w[ce] = I_Q2_GN;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_XV"))
    {
      pd_ptr->w[ce] = I_Q1_XV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_XV"))
    {
      pd_ptr->w[ce] = I_Q2_XV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_XG"))
    {
      pd_ptr->w[ce] = I_Q1_XG;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_XG"))
    {
      pd_ptr->w[ce] = I_Q2_XG;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_HV"))
    {
      pd_ptr->w[ce] = I_Q1_HV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_HG"))
    {
      pd_ptr->w[ce] = I_Q1_HG;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q1_HVG"))
    {
      pd_ptr->w[ce] = I_Q1_HVG;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_HV"))
    {
      pd_ptr->w[ce] = I_Q2_HV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_HG"))
    {
      pd_ptr->w[ce] = I_Q2_HG;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "Q2_HVG"))
    {
      pd_ptr->w[ce] = I_Q2_HVG;
      upd->XFEM = TRUE;
    }
    else if ( !strcasecmp(ts, "Q2_LSA") )
    {
      pd_ptr->w[ce] = I_Q2_LSA;
    }
    else if ( !strcasecmp(ts, "Q1_D") )
    {
      pd_ptr->w[ce] = I_Q1_D;
    }
    else if ( !strcasecmp(ts, "PQ1") )
    {
      pd_ptr->w[ce] = I_PQ1;
    }
    else if ( !strcasecmp(ts, "Q2_D") )
    {
      pd_ptr->w[ce] = I_Q2_D;
    }
    else if ( !strcasecmp(ts, "Q2_D_LSA") )
    {
      pd_ptr->w[ce] = I_Q2_D_LSA;
    }
    else if ( !strcasecmp(ts, "PQ2") )
    {
      pd_ptr->w[ce] = I_PQ2;
    }
    else if ( !strcasecmp(ts, "P0") )
    {
      pd_ptr->w[ce] = I_P0;
    }
    else if ( !strcasecmp(ts, "P1") )
    {
      pd_ptr->w[ce] = I_P1;
    }
    else if (!strcasecmp(ts, "P0_G"))
    {
      pd_ptr->w[ce] = I_P0_G;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P1_G"))
    {
      pd_ptr->w[ce] = I_P1_G;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P0_GP"))
    {
      pd_ptr->w[ce] = I_P0_GP;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P1_GP"))
    {
      pd_ptr->w[ce] = I_P1_GP;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P0_GN"))
    {
      pd_ptr->w[ce] = I_P0_GN;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P1_GN"))
    {
      pd_ptr->w[ce] = I_P1_GN;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P0_XV"))
    {
      pd_ptr->w[ce] = I_P0_XV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P1_XV"))
    {
      pd_ptr->w[ce] = I_P1_XV;
      upd->XFEM = TRUE;
    }
    else if (!strcasecmp(ts, "P1_XG"))
    {
      pd_ptr->w[ce] = I_P1_XG;
      upd->XFEM = TRUE;
    }
    else if ( !strcasecmp(ts, "SP") )
    {
      pd_ptr->w[ce] = I_SP;
      WH(-1,
	 "CAUTION WITH SUBPARAMETRIC: DON'T USE WITH LAGRANGIAN OR TOTAL_ALE MESH MOTIONS");
    }
    else
    {
      fprintf(stderr, 
	      "%s:\tUnrecognized Galerkin weighting function.\n", yo);
      exit(-1);
    }

    SPF(endofstring(echo_string),"\t%4s",ts);
    /*
     * Now for eqnMult greater than one, copy the entry we just made to
     * the other variable type entries instantiated by this one EQ line
     */
    if (eqnMult > 1) {
      for (j = 1; j < eqnMult; j++) {
        pd_ptr->w[ce+j] = pd_ptr->w[ce];
      }
    }
    
    /*
     * If this is a new and different Galerkin weight function, then
     * we'd better add it to our list of unique basis functions so that
     * later we can allocate just the right amount of memory and compute
     * just the right kinds of basis function derivatives that we need...
     */

    if ( in_list( pd_ptr->w[ce], 
		  0, Num_Interpolations, Unique_Interpolations) == -1 )
    {
      Unique_Interpolations[Num_Interpolations] = pd_ptr->w[ce];
      if ( Debug_Flag && ProcID == 0 )
      {
	fprintf(stdout, "Unique_Interpolations[%d] = %d\n",
		Num_Interpolations, pd_ptr->w[ce]);
      }
      Num_Interpolations++;
    }

    /*
     *
     * Part 3 -- Look at the variable name for this equation. This might
     *           be something like "T", "u1", "d1", etc.
     *
     */

    if ( fscanf(ifp, "%s", ts) !=1 )
    {
      fprintf(stderr, "%s: problem reading variable name for eqn.\n",yo);
      exit(-1);
    }

    /*
     * Check to see if this variable name is recognizable and if it
     * is reasonable.
     */

    if        (!strcasecmp(ts, "T")) {
      cv = set_var(TEMPERATURE, pd_ptr);
    } else if (!strcasecmp(ts, "U1")) {
      cv = set_var(VELOCITY1, pd_ptr);
    } else if (!strcasecmp(ts, "U2")) {
      cv = set_var(VELOCITY2, pd_ptr);
    } else if (!strcasecmp(ts, "U3")) {
      cv = set_var(VELOCITY3, pd_ptr);
    } else if (!strcasecmp(ts, "PU1")) {
      cv = set_var(PVELOCITY1, pd_ptr);
    } else if (!strcasecmp(ts, "PU2")) {
      cv = set_var(PVELOCITY2, pd_ptr);
    } else if (!strcasecmp(ts, "PU3")) {
      cv = set_var(PVELOCITY3, pd_ptr);
    } else if (!strcasecmp(ts, "S11")) {
      cv = set_var(POLYMER_STRESS11, pd_ptr);
    } else if (!strcasecmp(ts, "S12")) {
      cv = set_var(POLYMER_STRESS12, pd_ptr);
    } else if (!strcasecmp(ts, "S13")) {
      cv = set_var(POLYMER_STRESS13, pd_ptr);
    } else if (!strcasecmp(ts, "S22")) {
      cv = set_var(POLYMER_STRESS22, pd_ptr);
    } else if (!strcasecmp(ts, "S23")) {
      cv = set_var(POLYMER_STRESS23, pd_ptr);
    } else if (!strcasecmp(ts, "S33")) {
      cv = set_var(POLYMER_STRESS33, pd_ptr);
    } else if (!strcasecmp(ts, "G11")) {
      cv = set_var(VELOCITY_GRADIENT11, pd_ptr);
    } else if (!strcasecmp(ts, "G12")) {
      cv = set_var(VELOCITY_GRADIENT12, pd_ptr);
    } else if (!strcasecmp(ts, "G13")) {
      cv = set_var(VELOCITY_GRADIENT13, pd_ptr);
    } else if (!strcasecmp(ts, "G21")) {
      cv = set_var(VELOCITY_GRADIENT21, pd_ptr);
    } else if (!strcasecmp(ts, "G22")) {
      cv = set_var(VELOCITY_GRADIENT22, pd_ptr);
    } else if (!strcasecmp(ts, "G23")) {
      cv = set_var(VELOCITY_GRADIENT23, pd_ptr);
    } else if (!strcasecmp(ts, "G31")) {
      cv = set_var(VELOCITY_GRADIENT31, pd_ptr);
    } else if (!strcasecmp(ts, "G32")) {
      cv = set_var(VELOCITY_GRADIENT32, pd_ptr);
    } else if (!strcasecmp(ts, "G33")) {
      cv = set_var(VELOCITY_GRADIENT33, pd_ptr);
    } else if (!strcasecmp(ts, "NN")) {
      cv = set_var(BOND_EVOLUTION, pd_ptr);
    } else if (!strcasecmp(ts, "Y")) {
      cv = set_var(MASS_FRACTION, pd_ptr);
    } else if (!strcasecmp(ts, "D1")) {
      cv = set_var(MESH_DISPLACEMENT1, pd_ptr);
    } else if (!strcasecmp(ts, "D2")) {
      cv = set_var(MESH_DISPLACEMENT2, pd_ptr);
    } else if (!strcasecmp(ts, "D3")) {
      cv = set_var(MESH_DISPLACEMENT3, pd_ptr);
    } else if (!strcasecmp(ts, "D1_RS")) {
      cv = set_var(SOLID_DISPLACEMENT1, pd_ptr);
    } else if (!strcasecmp(ts, "D2_RS")) {
      cv = set_var(SOLID_DISPLACEMENT2, pd_ptr);
    } else if (!strcasecmp(ts, "D3_RS")) {
      cv = set_var(SOLID_DISPLACEMENT3, pd_ptr);
    } else if (!strcasecmp(ts, "P")) {
      cv = set_var(PRESSURE, pd_ptr);
    } else if (!strcasecmp(ts, "F")) {
      cv = set_var(FILL, pd_ptr);
    } else if (!strcasecmp(ts, "V")) {
      cv = set_var(VOLTAGE, pd_ptr);
    } else if (!strcasecmp(ts, "QS")) {
      cv = set_var(SURF_CHARGE, pd_ptr);
    } else if (!strcasecmp(ts, "SH")) {
      cv = set_var(SHEAR_RATE, pd_ptr);
    } else if (!strcasecmp(ts, "H")) {
      cv = set_var(CURVATURE, pd_ptr);
    } else if (!strcasecmp(ts, "N1")) {
      cv = set_var(NORMAL1, pd_ptr);
    } else if (!strcasecmp(ts, "N2")) {
      cv = set_var(NORMAL2, pd_ptr);
    } else if (!strcasecmp(ts, "N3")) {
      cv = set_var(NORMAL3, pd_ptr);
    } else if (!strcasecmp(ts, "S")) {
      cv = set_var(SURFACE, pd_ptr);
    } else if (!strcasecmp(ts, "P_LIQ")) {
      cv = set_var(POR_LIQ_PRES, pd_ptr);
    } else if (!strcasecmp(ts, "P_GAS")) {
      cv = set_var(POR_GAS_PRES, pd_ptr);
    } else if (!strcasecmp(ts, "P_POR")) {
      cv = set_var(POR_POROSITY, pd_ptr);
    } else if (!strcasecmp(ts, "P_TEMP")) {
      cv = set_var(POR_TEMP, pd_ptr);
    } else if (!strcasecmp(ts, "P_SAT")) {
      cv = set_var(POR_SATURATION, pd_ptr);
    } else if (!strcasecmp(ts, "VD1")) {
      cv = set_var(VORT_DIR1, pd_ptr);
    } else if (!strcasecmp(ts, "VD2")) {
      cv = set_var(VORT_DIR2, pd_ptr);
    } else if (!strcasecmp(ts, "VD3")) {
      cv = set_var(VORT_DIR3, pd_ptr);
    } else if (!strcasecmp(ts, "VLAMBDA")) {
      cv = set_var(VORT_LAMBDA, pd_ptr);
    } else if (!strcasecmp(ts, "LM1")) {
      cv = set_var(LAGR_MULT1, pd_ptr);
    } else if (!strcasecmp(ts, "LM2")) {
      cv = set_var(LAGR_MULT2, pd_ptr);
    } else if (!strcasecmp(ts, "LM3")) {
      cv = set_var(LAGR_MULT3, pd_ptr);
    } else if (!strcasecmp(ts, "K")) {
      cv = set_var(SHELL_CURVATURE, pd_ptr);
    } else if (!strcasecmp(ts, "K2")) {
      cv = set_var(SHELL_CURVATURE2, pd_ptr);
    } else if (!strcasecmp(ts, "TENS")) {
      cv = set_var(SHELL_TENSION, pd_ptr);
    } else if (!strcasecmp(ts, "SH_X")) {
      cv = set_var(SHELL_X, pd_ptr);
    } else if (!strcasecmp(ts, "SH_Y")) {
      cv = set_var(SHELL_Y, pd_ptr);
    } else if (!strcasecmp(ts, "SH_U")) {
      cv = set_var(SHELL_USER, pd_ptr);
    } else if (!strcasecmp(ts, "SH_ANG1")) {
      cv = set_var(SHELL_ANGLE1, pd_ptr);
    } else if (!strcasecmp(ts, "SH_ANG2")) {
      cv = set_var(SHELL_ANGLE2, pd_ptr);
    } else if (!strcasecmp(ts, "gamma1")) {
      cv = set_var(SHELL_SURF_DIV_V, pd_ptr);
    } else if (!strcasecmp(ts, "gamma2")) {
      cv = set_var(SHELL_SURF_CURV, pd_ptr);
    } else if (!strcasecmp(ts, "gamma4")) {
      cv = set_var(N_DOT_CURL_V, pd_ptr);
    } else if (!strcasecmp(ts, "gamma3_1")) {
      cv = set_var(GRAD_S_V_DOT_N1, pd_ptr);
    } else if (!strcasecmp(ts, "gamma3_2")) {
      cv = set_var(GRAD_S_V_DOT_N2, pd_ptr);
    } else if (!strcasecmp(ts, "gamma3_3")) {
      cv = set_var(GRAD_S_V_DOT_N3, pd_ptr);  
    } else if (!strcasecmp(ts, "APR")) {
      cv = set_var(ACOUS_PREAL, pd_ptr);  
    } else if (!strcasecmp(ts, "API")) {
      cv = set_var(ACOUS_PIMAG, pd_ptr);  
    } else if (!strcasecmp(ts, "ARS")) {
      cv = set_var(ACOUS_REYN_STRESS, pd_ptr);  
    } else if (!strcasecmp(ts, "EPR")) {
      cv = set_var(EM_CONT_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EPI")) {
      cv = set_var(EM_CONT_IMAG, pd_ptr);  
    } else if (!strcasecmp(ts, "SH_BV")) {
      cv = set_var(SHELL_BDYVELO, pd_ptr);  
    } else if (!strcasecmp(ts, "SH_P")) {
      cv = set_var(SHELL_LUBP, pd_ptr);
    } else if (!strcasecmp(ts, "LUBP")) {
      cv = set_var(LUBP, pd_ptr);
    } else if (!strcasecmp(ts, "LUBP_2")) {
      cv = set_var(LUBP_2, pd_ptr);
    } else if (!strcasecmp(ts, "SHELL_FILMP")) {
      cv = set_var(SHELL_FILMP, pd_ptr);
    } else if (!strcasecmp(ts, "SHELL_FILMH")) {
      cv = set_var(SHELL_FILMH, pd_ptr);
    } else if (!strcasecmp(ts, "SHELL_PARTC")) {
      cv = set_var(SHELL_PARTC, pd_ptr);
    } else if (!strcasecmp(ts, "SH_SAT_CLOSED")) {
      cv = set_var(SHELL_SAT_CLOSED, pd_ptr);
    } else if (!strcasecmp(ts, "SH_P_OPEN")) {
      cv = set_var(SHELL_PRESS_OPEN, pd_ptr);
    } else if (!strcasecmp(ts, "SH_P_OPEN_2")) {
      cv = set_var(SHELL_PRESS_OPEN_2, pd_ptr);
    } else if (!strcasecmp(ts, "SH_T")) {
      cv = set_var(SHELL_TEMPERATURE, pd_ptr);
    } else if (!strcasecmp(ts, "SH_DH")) {
      cv = set_var(SHELL_DELTAH, pd_ptr);
    } else if (!strcasecmp(ts, "SH_L_CURV")) {
      cv = set_var(SHELL_LUB_CURV, pd_ptr);
    } else if (!strcasecmp(ts, "SH_L_CURV_2")) {
      cv = set_var(SHELL_LUB_CURV_2, pd_ptr);
    } else if (!strcasecmp(ts, "SH_SAT_GASN")) {
      cv = set_var(SHELL_SAT_GASN, pd_ptr);
    } else if (!strcasecmp(ts, "SINK_MASS")) {
      cv = set_var(POR_SINK_MASS, pd_ptr);
    } else if (!strcasecmp(ts, "SH_SHEAR_TOP")) {
      cv = set_var(SHELL_SHEAR_TOP, pd_ptr);
    } else if (!strcasecmp(ts, "SH_SHEAR_BOT")) {
      cv = set_var(SHELL_SHEAR_BOT, pd_ptr);
    } else if (!strcasecmp(ts, "SH_CROSS_SHEAR")) {
      cv = set_var(SHELL_CROSS_SHEAR, pd_ptr); 
    } else if (!strcasecmp(ts, "MAX_STRAIN")) {
      cv = set_var(MAX_STRAIN, pd_ptr);   
    } else if (!strcasecmp(ts, "CUR_STRAIN")) {
      cv = set_var(CUR_STRAIN, pd_ptr);       
    } else if (!strcasecmp(ts, "SH_J")) {
      cv = set_var(SHELL_DIFF_FLUX, pd_ptr);
    } else if (!strcasecmp(ts, "SH_KD")) {
      cv = set_var(SHELL_DIFF_CURVATURE, pd_ptr);
    } else if (!strcasecmp(ts, "SH_N1")) {
      cv = set_var(SHELL_NORMAL1, pd_ptr);
    } else if (!strcasecmp(ts, "SH_N2")) {
      cv = set_var(SHELL_NORMAL2, pd_ptr);
    } else if (!strcasecmp(ts, "SH_N3")) {
      cv = set_var(SHELL_NORMAL3, pd_ptr);
    } else if (!strcasecmp(ts, "EXT_V")) {
      cv = set_var(EXT_VELOCITY, pd_ptr);

    } else if (!strcasecmp(ts, "INTP")) {
      cv = set_var(LIGHT_INTP, pd_ptr);
    } else if (!strcasecmp(ts, "INTM")) {
      cv = set_var(LIGHT_INTM, pd_ptr);
    } else if (!strcasecmp(ts, "INTD")) {
      cv = set_var(LIGHT_INTD, pd_ptr);
    } else if (!strcasecmp(ts, "E1")) {
      cv = set_var(EFIELD1, pd_ptr);
    } else if (!strcasecmp(ts, "E2")) {
      cv = set_var(EFIELD2, pd_ptr);
    } else if (!strcasecmp(ts, "E3")) {
      cv = set_var(EFIELD3, pd_ptr);
    } else if (!strcasecmp(ts, "F1")) {
      cv = set_var(PHASE1, pd_ptr);
    } else if (!strcasecmp(ts, "F2")) {
      cv = set_var(PHASE2, pd_ptr);
    } else if (!strcasecmp(ts, "F3")) {
      cv = set_var(PHASE3, pd_ptr);
    } else if (!strcasecmp(ts, "F4")) {
      cv = set_var(PHASE4, pd_ptr);
    } else if (!strcasecmp(ts, "F5")) {
      cv = set_var(PHASE5, pd_ptr);

    } else if (!strcasecmp(ts, "ENORM")) {
      cv = set_var(ENORM, pd_ptr);

    } else if (!strcasecmp(ts, "TFMP_PRES")) {
      cv = set_var(TFMP_PRES, pd_ptr);
    } else if (!strcasecmp(ts, "TFMP_SAT")) {
      cv = set_var(TFMP_SAT, pd_ptr);
    } else if (!strcasecmp(ts, "RST")) {
      cv = set_var(RESTIME, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_E1_REAL")) {
      cv = set_var(EM_E1_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_E2_REAL")) {
      cv = set_var(EM_E2_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_E3_REAL")) {
      cv = set_var(EM_E3_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_E1_IMAG")) {
      cv = set_var(EM_E1_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "EM_E2_IMAG")) {
      cv = set_var(EM_E2_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "EM_E3_IMAG")) {
      cv = set_var(EM_E3_IMAG, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_H1_REAL")) {
      cv = set_var(EM_H1_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_H2_REAL")) {
      cv = set_var(EM_H2_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_H3_REAL")) {
      cv = set_var(EM_H3_REAL, pd_ptr);  
    } else if (!strcasecmp(ts, "EM_H1_IMAG")) {
      cv = set_var(EM_H1_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "EM_H2_IMAG")) {
      cv = set_var(EM_H2_IMAG, pd_ptr);
    } else if (!strcasecmp(ts, "EM_H3_IMAG")) {
      cv = set_var(EM_H3_IMAG, pd_ptr);  

    } else if (!strncasecmp(ts, "Sp", 2)) {
      if (!strcasecmp(ts, "Sp")) {
	cv = SPECIES_UNK_0;
	j =  SPECIES_UNK_0;
	for (i = 0; i < pd_ptr->Num_Species_Eqn; i++, j++) {
	  (void) set_var(j, pd_ptr);
	}
	if (eqnMult != pd_ptr->Num_Species_Eqn) {
          fprintf(stderr, "%s: Incompatible EQ Sp line: %d %d\n",
		  yo, eqnMult, pd_ptr->Num_Species_Eqn);
	  EH(-1, "rd_eq_specs error");
	}
      } else {
	if (!interpret_int(ts + 3, &j)) {
	  fprintf(stderr,
		  "%s: Unable to decipher species var labeled: %s\n",
		  yo, ts);
	  exit(-1);
	}
	if (j < 0 || j > MAX_SPECIES_UNK_NUM) {
	  fprintf(stderr,
		  "%s: Species var %d is out of the permissible range\n",
		  yo, j);
	  exit(-1);
	}
	if (eqnMult != 1) {
          fprintf(stderr, "%s: Incompatible EQ Sp_# line: %d %d\n",
		  yo, eqnMult, pd_ptr->Num_Species_Eqn);
	  EH(-1, "rd_eq_specs error");
	}
	cv = set_var(R_SPECIES_UNK_0 + j, pd_ptr);
      }
    } else {
      fprintf(stderr, "%s:\tUnrecognized variable.\n", yo);
      exit(-1);
    }

    SPF(endofstring(echo_string)," %7s",ts);


    /*
     * Setup actual variable to problem variable index arrays
     *  HKM -> I'm not sure these problem variable index arrays
     *         are valid any more for discontinuous variable
     *         problems. However, I will do my best to continue
     *         their implementation. 
     */
    for (j = 0; j < eqnMult; j++) {
      if (upd->vp[cv + j] == -1) {
	upd->vp[cv + j] = upd->Total_Num_Var;
	upd->Total_Num_Var++;
      }
    }

    /*
     *
     * Part 4 -- Obtain the interpolation function for the variable that
     * is associated with this equation. It might be something like
     * 	"P0" -- piecewise constant
     * 	"P1" -- piecewise linear
     * 	"Q1" -- Lagrangian { ,bi,tri}linear
     * 	"Q2" -- Lagrangian { ,bi,tri}quadratic
     * 	"SP" -- Subparametric { ,bi,tri}quadratic
     *
     */

    if ( fscanf(ifp, "%s", ts) !=1 )
    {
      fprintf(stderr, 
	      "%s: problem reading interpolation function for variable.\n",
	      yo);
      exit(-1);
    }
    
    /*
     *   Check to see if this interpolation function is recognizable
     *   and if it is reasonable.
     *
     * HKM -> For discontinuous interpolations that are flagged here,
     *        using Q1_D and Q2_D interpolations, we also flag the v[]
     *        field to denote the fact that they will be using
     *        discontinuous interpolations at material boundaries.
     *       NOTE -> P0 and P1 should also perhaps have v[] flagged,
     *         since eventually they will be written to a nodal based
     *         exodus file (thus necessatitating a discontinuous
     *         interpolation at the node). This will be addressed later.
     */
    if         (!strcasecmp(ts, "Q1")) {
      pd_ptr->i[cv] = I_Q1;
    } else if (!strcasecmp(ts, "Q2"))  {
      pd_ptr->i[cv] = I_Q2;
    } else if (!strcasecmp(ts, "Q1_G"))  {
      pd_ptr->i[cv] = I_Q1_G;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_G"))  {
      pd_ptr->i[cv] = I_Q2_G;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_GP"))  {
      pd_ptr->i[cv] = I_Q1_GP;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_GP"))  {
      pd_ptr->i[cv] = I_Q2_GP;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_GN"))  {
      pd_ptr->i[cv] = I_Q1_GN;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_GN"))  {
      pd_ptr->i[cv] = I_Q2_GN;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_XV"))  {
      pd_ptr->i[cv] = I_Q1_XV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_XV"))  {
      pd_ptr->i[cv] = I_Q2_XV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_XG"))  {
      pd_ptr->i[cv] = I_Q1_XG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_XG"))  {
      pd_ptr->i[cv] = I_Q2_XG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_HV"))  {
      pd_ptr->i[cv] = I_Q1_HV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_HG"))  {
      pd_ptr->i[cv] = I_Q1_HG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q1_HVG"))  {
      pd_ptr->i[cv] = I_Q1_HVG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_HV"))  {
      pd_ptr->i[cv] = I_Q2_HV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_HG"))  {
      pd_ptr->i[cv] = I_Q2_HG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_HVG"))  {
      pd_ptr->i[cv] = I_Q2_HVG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "Q2_LSA"))  {
      pd_ptr->i[cv] = I_Q2_LSA;
    } else if (!strcasecmp(ts, "Q1_D")) {
      pd_ptr->i[cv] = I_Q1_D;
      pd_ptr->v[cv] |= V_MATSPECIFIC;
    } else if (!strcasecmp(ts, "PQ1"))  {
      pd_ptr->i[cv] = I_PQ1;
    } else if (!strcasecmp(ts, "Q2_D")) {
      pd_ptr->i[cv] = I_Q2_D;
      pd_ptr->v[cv] |= V_MATSPECIFIC;
    } else if (!strcasecmp(ts, "Q2_D_LSA")) {
      pd_ptr->i[cv] = I_Q2_D_LSA;
      pd_ptr->v[cv] |= V_MATSPECIFIC;
    } else if (!strcasecmp(ts, "PQ2")) {
      pd_ptr->i[cv] = I_PQ2;
    } else if (!strcasecmp(ts, "P0")) {
      pd_ptr->i[cv] = I_P0;
    } else if (!strcasecmp(ts, "P1")) {
      pd_ptr->i[cv] = I_P1;
    } else if (!strcasecmp(ts, "P0_G"))  {
      pd_ptr->i[cv] = I_P0_G;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P1_G"))  {
      pd_ptr->i[cv] = I_P1_G;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P0_GP"))  {
      pd_ptr->i[cv] = I_P0_GP;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P1_GP"))  {
      pd_ptr->i[cv] = I_P1_GP;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P0_GN"))  {
      pd_ptr->i[cv] = I_P0_GN;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P1_GN"))  {
      pd_ptr->i[cv] = I_P1_GN;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P0_XV"))  {
      pd_ptr->i[cv] = I_P0_XV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P1_XV"))  {
      pd_ptr->i[cv] = I_P1_XV;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "P1_XG"))  {
      pd_ptr->i[cv] = I_P1_XG;
      upd->XFEM = TRUE;
    } else if (!strcasecmp(ts, "SP")) {
      pd_ptr->i[cv] = I_SP;
    } else {
      fprintf(stderr, 
	      "%s:\tUnrecognized interpolation function for variable.\n", 
	      yo);
      exit(-1);
    }

    SPF(endofstring(echo_string)," %4s",ts);

    /*
     * If this is a new and different interpolation function, then
     * we'd better add it to our list of unique basis functions so that
     * later we can allocate just the right amount of memory and compute
     * just the right kinds of basis function derivatives that we need...
     */

    if (in_list(pd_ptr->i[cv], 
		 0, Num_Interpolations, Unique_Interpolations) == -1)
    {
      Unique_Interpolations[Num_Interpolations] = pd_ptr->i[cv];
      Num_Interpolations++;
    }

    /*
     *
     * Part 5 -- These describe the constant multipliers in front of each
     *           kind of term of this equation. This may be used as a poor
     *           man's version of thermophysical property specification, but
     *           is actually intended to aid in debugging and to find more
     *           efficient paths through the assembly process. For example,
     *           if the mass matrix term of the energy equation is zero, then
     *           a bunch of code can be circumvented.
     * 
     *           Here's a general legend for the meaning of the different
     * 	   floating point terms...
     * 
     *				M_______MASS term (d/dt)
     *
     *        			A_______ADVECTIVE term (v.grad())
     *
     * 				B_______BOUNDARY term Boundary Integral
     *  					of n.flux
     * 
     *				D_______DIFFUSIVE term Volume Integral
     * 					of grad(phi).flux
     * 
     *				S_______SOURCE term (anything)
     * 
     *				V_______DIVERGENCE term (div())
     * 
     * 				P_______POROUS term (linear source)
     * 
     * 		energy equation:	M,A,B,D,S
     * 		bulk species equation:	M,A,B,D,S
     * 		momentum equation:	M,A,B,D,S,P
     * 		mesh equation:		M,A,B,D,S
     * 		continuity equation:	A,S
     * 		velocity gradient:	A,S
     *		surf species equation:	-
     * 		porous media equation:	M,A,B,D,S
     *          vorticity direction:    -
     *          vorticity lambda:       -
     *          Lagrange Mult Eqn:      -
     *
     */
    /*
     * Thus, depending on the type of equation, the number of term multipliers
     * will be different.
     */

    switch ( ce ) {
      /*
       * Zero terms ...
       */
    case R_LAGR_MULT1:
    case R_LAGR_MULT2:
    case R_LAGR_MULT3:
    case R_VORT_LAMBDA:
    case R_SHELL_ANGLE1:
    case R_SHELL_ANGLE2:
    case R_SHELL_SURF_DIV_V:
    case R_SHELL_SURF_CURV:
    case R_N_DOT_CURL_V:
    case R_GRAD_S_V_DOT_N1:
    case R_GRAD_S_V_DOT_N2:
    case R_GRAD_S_V_DOT_N3:
      /* In case this actually gets used for a boolean anywhere ... */
      pd_ptr->etm[ce][(LOG2_MASS)] = 1.0;
      break;

      /*
       * One term ...
       */
    case R_SHELL_CURVATURE:
    case R_SHELL_CURVATURE2:
    case R_SHELL_TENSION:
    case R_SHELL_X:
    case R_SHELL_Y:
    case R_SHELL_DIFF_FLUX:
    case R_SHELL_DIFF_CURVATURE:
    case R_SHELL_NORMAL1:
    case R_SHELL_NORMAL2:
    case R_SHELL_NORMAL3:

      /* add a little consistency check for any 3D or cylindrical problems */
      if ((pd_ptr->CoordinateSystem != CARTESIAN &&
	   pd_ptr->CoordinateSystem != PROJECTED_CARTESIAN &&
	   pd_ptr->CoordinateSystem != CARTESIAN_2pt5D))
	{
	  //  EH(-1,"Shell capability only avaible for 2D Cartesian systems");
          printf("WARNING: Shell Capability is experimental for non 2D Cartesian systems\n");
	}
      if (fscanf(ifp, "%lf", &(pd_ptr->etm[ce][(LOG2_DIFFUSION)])) != 1)
        {
          sr = sprintf(err_msg, 
                       "Provide 1 equation term multiplier (dif) on %s in %s\n",
                       EQ_Name[ce].name1, pd_ptr->MaterialName);
          EH(-1, err_msg);
        }
      if (pd_ptr->etm[ce][(LOG2_DIFFUSION)] == 0.0) 
	{
	  sr = sprintf(err_msg, 
                       "Equation term multiplier (dif) on %s in %s is equal to zero. Will get a singular jacobian\n",
                       EQ_Name[ce].name1, pd_ptr->MaterialName);
          EH(-1, err_msg);
	}
      SPF(endofstring(echo_string),"\t %.4g", pd_ptr->etm[ce][(LOG2_DIFFUSION)]);
      break;


      /* 
       * Two terms.... 
       */
    case R_PRESSURE:
    case R_EM_CONT_REAL:
    case R_EM_CONT_IMAG:

	if ( fscanf(ifp, "%lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 2 )
	{
            pd_ptr->etm[ce][(LOG2_ADVECTION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_SOURCE)] = 0.0;
	    sr = sprintf(err_msg, 
		       "Using default equation term multipliers (adv,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	    WH(-1, err_msg);
	  fprintf(stderr,"\t %s %.4g %.4g \n",EQ_Name[ce].name1, 
                    pd_ptr->etm[ce][(LOG2_ADVECTION)], pd_ptr->etm[ce][(LOG2_SOURCE)]);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g", pd_ptr->etm[ce][(LOG2_ADVECTION)],
	                                            pd_ptr->etm[ce][(LOG2_SOURCE)]);
      break;
    case R_GRADIENT11:
    case R_GRADIENT12:
    case R_GRADIENT13:
    case R_GRADIENT21:
    case R_GRADIENT22:
    case R_GRADIENT23:
    case R_GRADIENT31:
    case R_GRADIENT32:
    case R_GRADIENT33:
	if ( fscanf(ifp, "%lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 2 )
	{
            pd_ptr->etm[ce][(LOG2_ADVECTION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_SOURCE)] = 1.0;
	    sr = sprintf(err_msg, 
		       "Using default equation term multipliers (adv,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	    WH(-1, err_msg);
	  fprintf(stderr,"\t %s %.4g %.4g \n",EQ_Name[ce].name1, 
                    pd_ptr->etm[ce][(LOG2_ADVECTION)], pd_ptr->etm[ce][(LOG2_SOURCE)]);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g", pd_ptr->etm[ce][(LOG2_ADVECTION)],
	                                            pd_ptr->etm[ce][(LOG2_SOURCE)]);
      break;
    case R_EFIELD1:
    case R_EFIELD2:
    case R_EFIELD3:
    case R_ENORM:
    case R_NORMAL1:
    case R_NORMAL2:
    case R_NORMAL3:
    case R_EXT_VELOCITY:
    case R_VORT_DIR1:
    case R_VORT_DIR2:
    case R_VORT_DIR3:
    case R_SHELL_SHEAR_TOP:
    case R_SHELL_SHEAR_BOT:
    case R_SHELL_CROSS_SHEAR:

	if ( fscanf(ifp, "%lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 2 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 2 equation term multipliers (adv,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g", pd_ptr->etm[ce][(LOG2_ADVECTION)],
	                                            pd_ptr->etm[ce][(LOG2_SOURCE)]);


	break;
    case R_SHELL_SAT_CLOSED:
    case R_SHELL_DELTAH:
      if ( fscanf(ifp, "%lf %lf", 
		  &(pd_ptr->etm[ce][(LOG2_MASS)]),
		  &(pd_ptr->etm[ce][(LOG2_DIVERGENCE)]))
	   != 2 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 2 equation term multipliers (mass,div) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}
      
      SPF( endofstring(echo_string),"\t %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
	   pd_ptr->etm[ce][(LOG2_DIVERGENCE)]);    
      break;


    case R_SHELL_SAT_GASN:
    case R_MAX_STRAIN:
    case R_CUR_STRAIN:
      if ( fscanf(ifp, "%lf %lf", 
		  &(pd_ptr->etm[ce][(LOG2_MASS)]),
		  &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	   != 2 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 2 equation term multipliers (mass,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}
      
      SPF( endofstring(echo_string),"\t %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
	   pd_ptr->etm[ce][(LOG2_SOURCE)]);    
      break;

      /* 
       * Three terms.... 
       */
    case R_FILL:
    case R_PHASE1:
    case R_PHASE2:
    case R_PHASE3:
    case R_PHASE4:
    case R_PHASE5:
    case R_ACOUS_REYN_STRESS:
    case R_SHELL_LUBP:
    case R_POR_SINK_MASS:

	if ( fscanf(ifp, "%lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_MASS)]),
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 3 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 3 equation term multipliers (mas,adv,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
                                                         pd_ptr->etm[ce][(LOG2_ADVECTION)],
	                                                 pd_ptr->etm[ce][(LOG2_SOURCE)]);

	break;
    case R_SHEAR_RATE:
	if ( fscanf(ifp, "%lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 3 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 3 equation term multipliers (adv,diff,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_ADVECTION)],
                                                         pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                 pd_ptr->etm[ce][(LOG2_SOURCE)]);
	break;
    case R_SHELL_LUB_CURV:
    case R_SHELL_LUB_CURV_2:
	if ( fscanf(ifp, "%lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_MASS)]),
		    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		    &(pd_ptr->etm[ce][(LOG2_DIVERGENCE)]))
	     != 3 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 3 equation term multipliers (mas,dif,div) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}

	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
                                                         pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                 pd_ptr->etm[ce][(LOG2_DIVERGENCE)]);
	break;
    case R_CURVATURE:
    case R_LUBP:
    case R_LUBP_2:
        if ( fscanf(ifp, "%lf %lf %lf", 
                    &(pd_ptr->etm[ce][(LOG2_BOUNDARY)]),
                    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
                    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
             != 3 )
        {
          sr = sprintf(err_msg, 
                       "Provide 3 equation term multipliers (diff,src, bnd) on %s in %s",
                       EQ_Name[ce].name1, pd_ptr->MaterialName);
          EH(-1, err_msg);
        }
	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_BOUNDARY)],
                                                         pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                 pd_ptr->etm[ce][(LOG2_SOURCE)]);
        break;
    case R_SHELL_SAT_OPEN:
    case R_SHELL_SAT_OPEN_2:
        if ( fscanf(ifp, "%lf %lf %lf", 
                    &(pd_ptr->etm[ce][(LOG2_MASS)]),
                    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
                    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
             != 3 )
        {
          sr = sprintf(err_msg, 
                       "Provide 3 equation term multipliers (mass,diff,src) on %s in %s",
                       EQ_Name[ce].name1, pd_ptr->MaterialName);
          EH(-1, err_msg);
        }
	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
                                                         pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                 pd_ptr->etm[ce][(LOG2_SOURCE)]);
        break;

    case R_TFMP_MASS:
    	if ( fscanf(ifp, "%lf %lf %lf",
		  &(pd_ptr->etm[ce][(LOG2_MASS)]),
		  &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		  &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]))
	      != 3 )
    	{
    	  sr = sprintf(err_msg,
                       "Provide 3 equation term multipliers (mass,adv,dif) on %s in %s",
					   EQ_Name[ce].name1, pd_ptr->MaterialName);
    	  EH(-1, err_msg);
    	}
    	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
    		  	  	  	  	   	   	   	   	   	   	   	   pd_ptr->etm[ce][(LOG2_ADVECTION)],
														   pd_ptr->etm[ce][(LOG2_DIFFUSION)]);
      break;
    case R_TFMP_BOUND:
      if ( fscanf(ifp, "%lf %lf %lf",
		  &(pd_ptr->etm[ce][(LOG2_MASS)]),
		  &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		  &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	   != 3 ) {
	sr = sprintf(err_msg,
		     "Provide 3 equation term multipliers (mass,adv,dif) on %s in %s",
		     EQ_Name[ce].name1, pd_ptr->MaterialName);
	EH(-1, err_msg);
      }
      SPF( endofstring(echo_string),"\t %.4g %.4g %.4g",
	   pd_ptr->etm[ce][(LOG2_MASS)],
	   pd_ptr->etm[ce][(LOG2_ADVECTION)],
	   pd_ptr->etm[ce][(LOG2_SOURCE)]);
      break;
      /* 
       * Four terms.... 
       */
    case R_BOND_EVOLUTION:

      if ( fscanf(ifp, "%lf %lf %lf  %lf", 
		  &(pd_ptr->etm[ce][(LOG2_MASS)]),
		  &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		  &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		  &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	   != 4 )
	{
	  sr = sprintf(err_msg, 
		       "Provide 4 equation term multipliers (mas,adv,diff, src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	  EH(-1, err_msg);
	}

      SPF( endofstring(echo_string),"\t %.4g %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
	   pd_ptr->etm[ce][(LOG2_ADVECTION)],
	   pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	   pd_ptr->etm[ce][(LOG2_SOURCE)]);

	break;



	/* 
	 * Five terms....  mesh-like 
	 */
    case R_MESH1:
    case R_MESH2:
    case R_MESH3:
    case R_SOLID1:
    case R_SOLID2:
    case R_SOLID3:

	if ( fscanf(ifp, "%lf %lf %lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_MASS)]),
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_BOUNDARY)]),
		    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 5 )
	{
            pd_ptr->etm[ce][(LOG2_MASS)] = 0.0; 
            pd_ptr->etm[ce][(LOG2_ADVECTION)] = 0.0;
	    pd_ptr->etm[ce][(LOG2_BOUNDARY)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_DIFFUSION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_SOURCE)] = 0.0;
	    sr = sprintf(err_msg, 
		       "Using default equation term multipliers (mas,adv,bnd,dif,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	    WH(-1, err_msg);
	  fprintf(stderr,"\t %s %.4g %.4g %.4g %.4g %.4g \n", EQ_Name[ce].name1,
               pd_ptr->etm[ce][(LOG2_MASS)],
	       pd_ptr->etm[ce][(LOG2_ADVECTION)], pd_ptr->etm[ce][(LOG2_BOUNDARY)],
               pd_ptr->etm[ce][(LOG2_DIFFUSION)], pd_ptr->etm[ce][(LOG2_SOURCE)]);
	}
	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
                                                                   pd_ptr->etm[ce][(LOG2_ADVECTION)],
                                                                   pd_ptr->etm[ce][(LOG2_BOUNDARY)],
                                                                   pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                           pd_ptr->etm[ce][(LOG2_SOURCE)]);
	break;
	/* 
	 * Five terms....  other
	 */
    case R_ENERGY:
    case R_POTENTIAL:         /* KSC: 2/99 */ 
    case R_MASS:
    case R_POR_LIQ_PRES:
    case R_POR_GAS_PRES:
    case R_POR_POROSITY:
    case R_POR_ENERGY:
    case R_POR_SATURATION:
    case R_STRESS11:
    case R_STRESS12:
    case R_STRESS13:
    case R_STRESS22:
    case R_STRESS23:
    case R_STRESS33:
    case R_SURF_CHARGE:
    case R_SHELL_USER:
    case R_SHELL_BDYVELO:
    case R_ACOUS_PREAL:
    case R_ACOUS_PIMAG:
    case R_SHELL_FILMP:
    case R_SHELL_FILMH:
    case R_SHELL_PARTC:
    case R_SHELL_ENERGY:
    case R_LIGHT_INTP:
    case R_LIGHT_INTM:
    case R_LIGHT_INTD:
    case R_RESTIME:  
    case R_EM_E1_REAL:
    case R_EM_E2_REAL:
    case R_EM_E3_REAL:
    case R_EM_E1_IMAG:
    case R_EM_E2_IMAG:
    case R_EM_E3_IMAG:
    case R_EM_H1_REAL:
    case R_EM_H2_REAL:
    case R_EM_H3_REAL:
    case R_EM_H1_IMAG:
    case R_EM_H2_IMAG:
    case R_EM_H3_IMAG:

	if ( fscanf(ifp, "%lf %lf %lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_MASS)]),
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_BOUNDARY)]),
		    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]))
	     != 5 )
	{
            if(TimeIntegration == TRANSIENT)
                { pd_ptr->etm[ce][(LOG2_MASS)] = 1.0; }
            else
                { pd_ptr->etm[ce][(LOG2_MASS)] = .0; }
            pd_ptr->etm[ce][(LOG2_ADVECTION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_BOUNDARY)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_DIFFUSION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_SOURCE)] = 1.0;
	    sr = sprintf(err_msg, 
		       "Using default equation term multipliers (mas,adv,bnd,dif,src) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	    WH(-1, err_msg);
	  fprintf(stderr,"\t %s %.4g %.4g %.4g %.4g %.4g \n", EQ_Name[ce].name1,
               pd_ptr->etm[ce][(LOG2_MASS)],
	       pd_ptr->etm[ce][(LOG2_ADVECTION)], pd_ptr->etm[ce][(LOG2_BOUNDARY)],
               pd_ptr->etm[ce][(LOG2_DIFFUSION)], pd_ptr->etm[ce][(LOG2_SOURCE)]);
	}
	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
                                                                   pd_ptr->etm[ce][(LOG2_ADVECTION)],
                                                                   pd_ptr->etm[ce][(LOG2_BOUNDARY)],
                                                                   pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                           pd_ptr->etm[ce][(LOG2_SOURCE)]);
	break;

	/* 
	 * Six terms.... 
	 */
    case R_MOMENTUM1:
    case R_MOMENTUM2:
    case R_MOMENTUM3:
    case R_PMOMENTUM1:
    case R_PMOMENTUM2:
    case R_PMOMENTUM3:
	if ( fscanf(ifp, "%lf %lf %lf %lf %lf %lf", 
		    &(pd_ptr->etm[ce][(LOG2_MASS)]),
		    &(pd_ptr->etm[ce][(LOG2_ADVECTION)]),
		    &(pd_ptr->etm[ce][(LOG2_BOUNDARY)]),
		    &(pd_ptr->etm[ce][(LOG2_DIFFUSION)]),
		    &(pd_ptr->etm[ce][(LOG2_SOURCE)]),
		    &(pd_ptr->etm[ce][(LOG2_POROUS_BRINK)]))
	     != 6 )
	{
            if(TimeIntegration == TRANSIENT)
                { pd_ptr->etm[ce][(LOG2_MASS)] = 1.0; }
            else
                { pd_ptr->etm[ce][(LOG2_MASS)] = .0; }
            pd_ptr->etm[ce][(LOG2_ADVECTION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_BOUNDARY)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_DIFFUSION)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_SOURCE)] = 1.0;
	    pd_ptr->etm[ce][(LOG2_POROUS_BRINK)] = 0.0;
	    sr = sprintf(err_msg, 
		       "Using default equation term multipliers (mas,adv,bnd,dif,src,prs) on %s in %s",
		       EQ_Name[ce].name1, pd_ptr->MaterialName);
	    WH(-1, err_msg);
	  fprintf(stderr,"\t %s %.4g %.4g %.4g %.4g %.4g %.4g \n", EQ_Name[ce].name1,
               pd_ptr->etm[ce][(LOG2_MASS)], pd_ptr->etm[ce][(LOG2_ADVECTION)], 
               pd_ptr->etm[ce][(LOG2_BOUNDARY)], pd_ptr->etm[ce][(LOG2_DIFFUSION)],
               pd_ptr->etm[ce][(LOG2_SOURCE)], pd_ptr->etm[ce][(LOG2_POROUS_BRINK)]);
	}
	SPF( endofstring(echo_string),"\t %.4g %.4g %.4g %.4g %.4g %.4g", pd_ptr->etm[ce][(LOG2_MASS)],
	                                                                pd_ptr->etm[ce][(LOG2_ADVECTION)],
                                                                        pd_ptr->etm[ce][(LOG2_BOUNDARY)],
                                                                        pd_ptr->etm[ce][(LOG2_DIFFUSION)],
	                                                                pd_ptr->etm[ce][(LOG2_SOURCE)],
                                                                        pd_ptr->etm[ce][(LOG2_POROUS_BRINK)]);
	break;

    case R_MASS_SURF:
	/*
	 * No terms; not yet implemented...
	 */
	fprintf(stderr, "%s: No terms permitted here yet.\n",
		yo);
	exit(-1);
	break;
    default:
	fprintf(stderr, "%s: Unclassifiable EQ case. Exiting.\n", yo);
	exit(-1);
    }

    ECHO(echo_string,echo_file);
  }      
  

  /*
   *  Part 6 --
   * 
   *   Now figure out the mapping order, shape variables and
   *   the projection variable.
   *
   *   - pd_ptr->IntegrationMap variable is determined here
   *   - pd_ptr->ShapeVar       variable also
   *     pd_ptr->ProjectionVar
   */
  determine_ShapeVar(pd_ptr);
  determine_ProjectionVar(pd_ptr);
  
  /*
   *  Part 7 --
   * 
   *       In this section, we turn on the rest of the stress equations
   *       based on the number of modes specified and the base stress
   *       equation specified.
   */
   
  ce = R_STRESS11;
  cv = POLYMER_STRESS11;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS11_1;
    cvm = POLYMER_STRESS11_1;
    for(n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS11 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS11 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }

      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }

      /* move on to the 11 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS12;
  cv = POLYMER_STRESS12;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS12_1;
    cvm = POLYMER_STRESS12_1;
    for(n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS12 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS12 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 12 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS22;
  cv = POLYMER_STRESS22;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS22_1;
    cvm = POLYMER_STRESS22_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS22 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,
	 R_STRESS22 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 22 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS13;
  cv = POLYMER_STRESS13;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS13_1;
    cvm = POLYMER_STRESS13_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS13 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS13 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 13 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS23;
  cv = POLYMER_STRESS23;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS23_1;
    cvm = POLYMER_STRESS23_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS23 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS23 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 23 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS33;
  cv = POLYMER_STRESS33;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS33_1;
    cvm = POLYMER_STRESS33_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS33 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS33 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
 
      /* move on to the 33 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  /*
   *  Store the total number of equations that are active in this
   *   material in the Problem_Description structure
   */
  pd_ptr->Num_EQ = neqn;

  /*
   *  Check on tha bounds on the total number of equations permissible
   *  in a single material (NOTE: I don't know where the bounds comes
   *  from)
   */
  if ((neqn < 1 || neqn > 76)) {
    EH(-1,
       "rd_eq_specs: Too many (>76) or too few (<1) eqns. in this mat");
  }
  
  /* Consistency diagnostics for runs involving 3D stability of
   * 2D flows... */
  if(Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE)
    {
      if(!pd_glob[mn]->e[R_MOMENTUM3] || !pd_glob[mn]->v[VELOCITY3])
	{
	  fprintf(stderr, "\nR_MOMENTUM3/VELOCITY3 are required for a 3D stability of a 2D flow.\n");
	  fprintf(stderr, "They are missing in material %d (0-based)\n\n", mn);
	  EH(-1, "missing equation for 3D stability of 2D flow");
	}
    }

} /* END rd_eq_specs() -- read input file for equation & term specs *******/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*
 * look_for_mat_prop -- generic code to read a material property ( constant or 
 *                      user type variation)
 *
 * Comments:	This code was lifted out of the ever-growing rd_mp()
 *              which was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *
 *       Reads the string:
 *           
 *  "search_string" = "model_name"  [species_index]   value
 *  "search_string" = "model_name"  [species_index]   value1 value2 value3
 *
 *    It understands the following models:
 *        CONSTANT
 *        RATIO
 *        USER
 *        USER_GEN
 * 
 *  "search_string" =  CONSTANT  [species_index]   value
 *  "search_string" =  CONSTANT  [species_index]   value1 value2 value3
 *
 *  Return Value:
 *    -1 : This number is returned if the search string, "search_string"
 *         is not found in the material property file.
 *         This number is also returned if the model name is not one of the
 *         special cases it knows how to parse.
 *     1 : Successful completion of the routine.
 *
 * Revisions:
 *
 * 1997/06/23 13:21 MDT pasacki@sandia.gov
 *
 *   -- Added User_count argument to save the count of User_constants
 *	for variable number of user defined constants. Needed for parallel
 *      processing communication of the lengths of these lists.
 */

int 
look_for_mat_prop(FILE *imp,			/* ptr to input stream (in)*/
		  char *search_string,          /* search string (in) */
		  int  *MaterialModel,		/* int material model (out)*/
		  dbl  *Material_property,	/* double value matprop, 
						 * if constant (out) */
		  dbl  **User_constants,	/* ptr to vector of double 
						 * constants for user defined 
						 * material ppties (out) */
		  int  *User_count,		/* how many User_constants ? 
						 * (out) */
		  char *model_name,		/* ptr to model name (out) */
		  const int  num_values,        /* num. values for constant 
						 * models SCALAR=1, VECTOR=3 
						 * (in) */
		  int *read_species	,	/* flag whether to read species index
						 * or not 
						 * comes in here. If this is greater
						 * than zero, then a species ID
						 * is expected on the line.
						 * This number is
						 * treated as the Max Species
						 * number.
						 * On output this is equal to
						 * species ID, read from the
						 * file.*/ 
		  char *echo_string  )  /*           This char array will pass back the 
		                                 * reconstructed input card for output 
                                                 * to the echo file for this mat */
{
  char err_msg[MAX_CHAR_IN_INPUT];
#ifdef DEBUG
  static const char yo[] = "look_for_mat_prop";
#endif
  char	input[MAX_CHAR_IN_INPUT] = "zilch\0"; /* storage for input strings */
  int iread = -1; /* status flag  */
  double a0, a1, a2;			/* dummy for storing input properties */
  int DumModel=0;			/* dummy int for story input model */
  char  line[132];
  char  *arguments[MAX_NUMBER_PARAMS];

  int num_const, i, species_no, got_it;

  species_no = 0;		/* set species number to zero as default
				   so it works for non-species properties*/
  model_name[0] = '\0';        /* set the model name to empty */

  got_it = look_forward_optional(imp, search_string, input, '=');

#ifdef DEBUG
  fprintf(stderr, "%s: %d=look_forward_optional(imp, \"%s\", ... )\n", yo,
	  got_it, search_string);
#endif

  /*  look_for(imp, search_string, input, '='); */
  if (got_it == 1) 
    {
    /*
     * Read the Model name
     */
    if (fscanf(imp, "%s", model_name) != 1)
      {
	sr = sprintf(err_msg, 
	     "Error reading model name string, mat file \"%s\", property %s",
		     current_mat_file_name,
		     search_string);
	EH(-1, err_msg);
      }

    SPF(echo_string,"%s = %s",search_string, model_name );

    /*
     *   Read the species index, if requested to from the argument list
     */
    if (*read_species >= 0)
      {
	if (fscanf(imp, "%d", &species_no) != 1)
	  {
	    sr = sprintf(err_msg, 
                 "Error reading species number, mat file \"%s\", property %s",
			 current_mat_file_name,
			 search_string);
	    EH(-1, err_msg);
	  }
	if (species_no >= *read_species || species_no < 0)
	  {
	    sr = sprintf(err_msg,
		 "Illegal species number (%d), mat file \"%s\", property %s",
			 species_no,
			 current_mat_file_name,
			 search_string);
	    EH(-1, err_msg);
	  }
	*read_species = species_no;

	SPF(endofstring(echo_string)," %d",  *read_species);
      }
    
    /*
     * Branch depending upont the Model Name read on the line
     * 
     *    -> CONSTANT: Depending upon the type of input defined in the
     *                 argument list (num_values)
     *                 Material_Property[]
     *                 Material_Model
     */
    if ( !strcmp(model_name, "CONSTANT")  ||  !strcmp(model_name, "RATIO")  )
      {
      if(!strcmp(model_name, "CONSTANT"))DumModel = CONSTANT;
      if(!strcmp(model_name, "RATIO"))DumModel = RATIO;
	if (num_values == SCALAR_INPUT)
	  {
	    if ( fscanf(imp, "%lf ", 
			&a0)
		 != 1 )
	      {
		sr = sprintf(err_msg,
		     "Expected 1 flt for CONSTANT model %s, mat file \"%s\"",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    Material_property[species_no] = a0;
	    iread = 1;

	    SPF_DBL_VEC( endofstring(echo_string),1, Material_property);
	  }
	else if (num_values == VECTOR_INPUT)
	  {
	    if ( fscanf(imp, "%lf %lf %lf", 
			&a0, &a1, &a2)
		 != 3 )
	      {
		sr = sprintf(err_msg,
		     "Expected 3 flts for CONSTANT model %s, mat file \"%s\"",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    Material_property[0] = a0;
	    Material_property[1] = a1;
	    Material_property[2] = a2;
	    iread = 1;

	    SPF_DBL_VEC( endofstring(echo_string),3, Material_property);

	  }
	MaterialModel[species_no] = DumModel;
      }
    else if ( !strcmp(model_name, "USER") )
      {
	DumModel = USER;
	
	if ( fgets(line, 132, imp) != NULL)
	  {
	    strip(line);
	    num_const = count_parameters(line);
	    
	    /* allocate space */
	    if( User_constants == NULL) 
	      {
		sr = sprintf(err_msg, 
			     "USER model for %s in mat \"%s\" requires space!",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    User_constants[species_no] = 
	      (dbl *)array_alloc(1, num_const, sizeof(dbl));
	    
	    /*
	     * Save the count. Again, the non-species cases are "[0]" giving
	     * the effective *User_count behavior. For real species, this
	     * is one of an array of counts, for a particular species number. 
	     */

	    User_count[species_no] = num_const;

	    /* parse parameters into little strings */ 
  	    tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);
	    
	    for(i=0; i< num_const; i++)
	      {
		User_constants[species_no][i] = atof(arguments[i]);
	      }

	    SPF_DBL_VEC( endofstring(echo_string),num_const, User_constants[species_no] );

	  }
	else 
	  {
	    sr = sprintf(err_msg,
		 "Error reading USER constants for property %s, mat \"%s\"",
			 search_string,
			 current_mat_file_name);
	    EH(-1, err_msg);
	  }
	MaterialModel[species_no] = DumModel;
	iread = 1;
      }
    else if ( !strcmp(model_name, "USER_GEN") )
      {
	DumModel = USER_GEN;
	
	if ( fgets(line, 132, imp) != NULL)
	  {
	    strip(line);
	    num_const = count_parameters(line);
	    
	    /* allocate space */
	    if(User_constants == NULL) 
	      {
		sr = sprintf(err_msg, 
		"USER_GEN model for property %s in mat \"%s\" requires space!",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }

	    User_count[species_no] = num_const;

	    User_constants[species_no] = (dbl *)array_alloc(1, num_const, sizeof(dbl));
	    
	    /* parse parameters into little strings */ 
	    tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);
	    

	    for(i=0; i< num_const; i++)
	      {
		User_constants[species_no][i] = atof(arguments[i]);
	      }

	    SPF_DBL_VEC( endofstring(echo_string),num_const, User_constants[species_no] );
	  }
	else 
	  {
	    sr = sprintf(err_msg,
	    "Error reading USER_GEN constants for property %s, mat \"%s\"",
			 search_string,
			 current_mat_file_name);
	    EH(-1, err_msg);
	  }
	MaterialModel[species_no] = DumModel;
	iread = 1;
      }
    else 
      {
	iread = -1;
      }
    }
  else 
    {
      iread = -1;
      sr = sprintf(model_name," ");
      echo_string[0] = '\0';
    }
  return iread;
}
/* end of look_for_mat_prop */
/**************************************************************************/
/*
 * look_for_mat_proptable -- generic code to read a material property ( constant or 
 *                      user type variation)
 *
 * Comments:	This code is a variation of the look_for_mat_prop routine
 *              used to add the MP Table option.  This can eventually replace 
 *              look_for_mat_prop if all MP use the TABLE option.
 *
 *
 *       Reads the string:
 *           
 *  "search_string" = "model_name"  [species_index]   value
 *  "search_string" = "model_name"  [species_index]   value1 value2 value3
 *
 *    It understands the following models:
 *        CONSTANT
 *        USER
 *        USER_GEN
 *        TABLE
 * 
 *  "search_string" =  CONSTANT  [species_index]   value
 *  "search_string" =  CONSTANT  [species_index]   value1 value2 value3
 *
 *  Return Value:
 *    -1 : This number is returned if the search string, "search_string"
 *         is not found in the material property file.
 *         This number is also returned if the model name is not one of the
 *         special cases it knows how to parse.
 *     1 : Successful completion of the routine.
 *
 * Revisions:
 *
 * 1997/06/23 13:21 MDT pasacki@sandia.gov
 * raroach 11/16/99 - add TABLE lookup option
 *
 *   -- Added User_count argument to save the count of User_constants
 *	for variable number of user defined constants. Needed for parallel
 *      processing communication of the lengths of these lists.
 */

int 
look_for_mat_proptable(FILE *imp,			/* ptr to input stream (in)*/
		  char *search_string,          /* search string (in) */
		  int  *MaterialModel,		/* int material model (out)*/
		  dbl  *Material_property,	/* double value matprop, 
						 * if constant (out) */
		  dbl  **User_constants,	/* ptr to vector of double 
						 * constants for user defined 
						 * material ppties (out) */
		  int  *User_count,		/* how many User_constants ? 
						 * (out) */
		  int  *tableindex,             /* ptr to integer counter for #tables */
		  char *model_name,		/* ptr to model name (out) */
		  const int  num_values,        /* num. values for constant 
						 * models SCALAR=1, VECTOR=3 
						 * (in) */
		  int *read_species,		/* flag to read species or not 
						 * comes in as Max_Number 
						 * species output at species 
						 * number of property read 
						 * (in) */
		  char *echo_string )          /*   string to copy echoed data to */
{
  char err_msg[MAX_CHAR_IN_INPUT];
#ifdef DEBUG
  static const char yo[] = "look_for_mat_proptable";
#endif
  char	input[MAX_CHAR_IN_INPUT] = "zilch\0"; /* storage for input strings */
  int iread = -1; /* status flag  */
  double a0, a1, a2;			/* dummy for storing input properties */
  int DumModel;			/* dummy int for story input model */
  char  line[132];
  char  *arguments[MAX_NUMBER_PARAMS];
  int num_const=0, i, species_no, got_it;
  struct  Data_Table *table_local;

  species_no = 0;		/* set species number to zero as default
				   so it works for non-species properties*/

  got_it = look_forward_optional(imp, search_string, input, '=');

#ifdef DEBUG
  fprintf(stderr, "%s: %d=look_forward_optional(imp, \"%s\", ... )\n", yo, got_it,
	  search_string);
#endif

  if (got_it == 1) {
    if (fscanf(imp, "%s", model_name) != 1)
      {
	sr = sprintf(err_msg, 
	     "Error reading model name string, mat file \"%s\", property %s",
		     current_mat_file_name,
		     search_string);
	EH(-1, err_msg);
      }

    SPF(echo_string,"%s = %s", search_string,model_name);

    if (*read_species >= 0)
      {
	if (fscanf(imp, "%d", &species_no) != 1)
	  {
	    sr = sprintf(err_msg, 
                 "Error reading species number, mat file \"%s\", property %s",
			 current_mat_file_name,
			 search_string);
	    EH(-1, err_msg);
	  }
	if (species_no >= *read_species || species_no < 0)
	  {
	    sr = sprintf(err_msg,
		 "Illegal species number (%d), mat file \"%s\", property %s",
			 species_no,
			 current_mat_file_name,
			 search_string);
	    EH(-1, err_msg);
	  }
	*read_species = species_no;

	SPF(endofstring(echo_string)," %d", species_no);

      }
    
    if ( !strcmp(model_name, "CONSTANT") )
      {
	DumModel = CONSTANT;
	if (num_values == SCALAR_INPUT)
	  {
	    if ( fscanf(imp, "%lf ", 
			&a0)
		 != 1 )
	      {
		sr = sprintf(err_msg,
		     "Expected 1 flt for CONSTANT model %s, mat file \"%s\"",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    Material_property[species_no] = a0;
	    iread = 1;

	    SPF_DBL_VEC(endofstring(echo_string),1, Material_property);

	  }
	else if (num_values == VECTOR_INPUT)
	  {
	    if ( fscanf(imp, "%lf %lf %lf", 
			&a0, &a1, &a2)
		 != 3 )
	      {
		sr = sprintf(err_msg,
		     "Expected 3 flts for CONSTANT model %s, mat file \"%s\"",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    Material_property[0] = a0;
	    Material_property[1] = a1;
	    Material_property[2] = a2;
	    iread = 1;
	    SPF_DBL_VEC(endofstring(echo_string),3, Material_property);
	  }
	MaterialModel[species_no] = DumModel;
      }
    else if ( !strcmp(model_name, "USER") )
      {
	DumModel = USER;
	
	if ( fgets(line, 132, imp) != NULL)
	  {
	    strip(line);
	    num_const = count_parameters(line);
	    
	    /* allocate space */
	    if( User_constants == NULL) 
	      {
		sr = sprintf(err_msg, 
			     "USER model for %s in mat \"%s\" requires space!",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }
	    User_constants[species_no] = 
	      (dbl *)array_alloc(1, num_const, sizeof(dbl));
	    
	    /*
	     * Save the count. Again, the non-species cases are "[0]" giving
	     * the effective *User_count behavior. For real species, this
	     * is one of an array of counts, for a particular species number. 
	     */

	    User_count[species_no] = num_const;

	    /* parse parameters into little strings */ 
	    tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);	    
	    
	    for(i=0; i< num_const; i++)
	      {
		User_constants[species_no][i] = atof(arguments[i]);
	      }
	  }
	else 
	  {
	    sr = sprintf(err_msg,
		 "Error reading USER constants for property %s, mat \"%s\"",
			 search_string,
			 current_mat_file_name);
	    EH(-1, err_msg);
	  }
	MaterialModel[species_no] = DumModel;
	iread = 1;

	SPF_DBL_VEC(endofstring(echo_string),num_const, User_constants[species_no]);

      }
    else if ( !strcmp(model_name, "USER_GEN") )
      {
	DumModel = USER_GEN;
	
	if ( fgets(line, 132, imp) != NULL)
	  {
	    strip(line);
	    num_const = count_parameters(line);
	    
	    /* allocate space */
	    if(User_constants == NULL) 
	      {
		sr = sprintf(err_msg, 
		"USER_GEN model for property %s in mat \"%s\" requires space!",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }

	    User_count[species_no] = num_const;

	    User_constants[species_no] = (dbl *)array_alloc(1, num_const, sizeof(dbl));
	    
	    /* parse parameters into little strings */ 
	    tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);
	    
	    for(i=0; i< num_const; i++)
	      {
		User_constants[species_no][i] = atof(arguments[i]);
	      }

	  }
      	else 
	  {
	    sr = sprintf(err_msg,
	    "Error reading USER_GEN constants for property %s, mat \"%s\"",
			 search_string,
			 current_mat_file_name);
	    EH(-1, err_msg);
	  }
	MaterialModel[species_no] = DumModel;
	iread = 1;
	SPF_DBL_VEC(endofstring(echo_string),num_const, User_constants[species_no]);
      }
     else if ( !strcmp(model_name, "TABLE") )
      {
	DumModel = TABLE;

	 if(User_constants == NULL) 
          {
	    sr = sprintf(err_msg, 
		"TABLE model for property %s in mat \"%s\" requires space!",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	  }

	  /* 
	   * Fall through for all unimplemented cases 
	   */
	  if( num_MP_Tables == MAX_MP_TABLES )
	    {
	      EH(-1, "Maximum TABLE_MPs exceeded .");
	    }

          table_local = (struct Data_Table *) smalloc (sizeof( struct Data_Table));
          MP_Tables[num_MP_Tables] = setup_table_MP (imp,table_local, search_string);
	  *tableindex = num_MP_Tables++;
	  	  
  	  MaterialModel[species_no] = DumModel;
	  iread = 1;
      }
   else 
      {
	iread = -1;
      }
  }
  else {
    iread = -1;
    sprintf(model_name," ");
  }
  return iread;
}
/* end of look_for_mat_proptable */




/*
 * look_for_modal_prop -- generic code to read a material property for multimode
 *                        input data ( constant or user type variation)
 *                      
 *
 * Comments:	This code was handles the fact that we potentially have
 *              eight constants and eight time constants and eight
 *              mobilities to worry about.
 */

int 
look_for_modal_prop(FILE *imp,	/* ptr to input stream (in)*/
		    const char *search_string, /* search string (in) */
		    const int modes, /* number of viscoelastic modes (in) */
		    int  *MaterialModel, /* int material model (out)*/
		    dbl *modal_const,  /* modal data (out) */
		    char *echo_string)    /*character array to pass back echoed input */
{
  char err_msg[MAX_CHAR_IN_INPUT];
#ifdef DEBUG
  static const char yo[] = "look_for_modal_prop";
#endif
  char	input[MAX_CHAR_IN_INPUT];             /* dummy storage for input strings */
  int iread = -1; /* status flag  */
  int DumModel;			              /* dummy int for story input model */
  char  line[132];
  char  *arguments[MAX_NUMBER_PARAMS];
  int num_const, i, got_it;
  char	model_name[MAX_CS_KEYWORD_LENGTH];

  
  got_it = look_forward_optional(imp, search_string, input, '=');
  if (got_it == 1) 
    { 
    if (fscanf(imp, "%s", model_name) != 1)
      {
	sr = sprintf(err_msg, 
	     "Error reading model name string, mat file \"%s\", property %s",
		     current_mat_file_name,
		     search_string);
	EH(-1, err_msg);
      }

    SPF(echo_string, "%s = %s", search_string, model_name );

    if ( !strcmp(model_name, "CONSTANT") )
      {
	DumModel = CONSTANT;
	
	if ( fgets(line, 132, imp) != NULL)
	  {
	    strip(line);
	    num_const = count_parameters(line);
	    
	    if(modes > num_const)
	      {
		sr = sprintf(err_msg,
			     "Error: you must have data for each mode  %s, mat \"%s\"",
			     search_string,
			     current_mat_file_name);
		EH(-1, err_msg);
	      }

	    /* parse parameters into little strings */ 
	    tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);	    
	    
	    for(i=0; i< num_const; i++)
	      {
		modal_const[i] = atof(arguments[i]);
	      }
	    SPF_DBL_VEC( endofstring(echo_string), num_const, modal_const);
	  }
	else 
	  {
	    sr = sprintf(err_msg,
		 "Error reading CONSTANT mulitmodal constants for property %s, mat \"%s\"",
			 search_string,
			 current_mat_file_name);
	    EH(-1, err_msg);
	  }
	*MaterialModel= DumModel;
	iread = 1;
      }
    else if(!strcmp(model_name, "HERSCHEL_BULKLEY"))
    {
      *MaterialModel= HERSCHEL_BULKLEY;
      iread = 1;
      printf("HERSCHEL_BULKLEY model used for %s\n", search_string);
    }
    else 
      {
	iread = -2;
      }
    }
  else {
    iread = -1;
    sprintf(model_name," ");
    echo_string[0] = '\0';
  }
  return iread;
}
/* end of look_for_modal_prop */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 * read_constants -- generic code to read an arbitrary number of double constants
 *                   on a line. It reads the line, counts the number of numbers,
 *                   allocates space for storage of that many numbers, and then
 *                   reads the line, and stores the numbers in that storage. 
 *                   The file pointer is advanced to the end of the line.
 *
 *                   The number of constants read is returned in the return
 *                   variable.
 */

int 
read_constants(FILE *imp,	     /* pointer to file */
	       dbl **User_constants, /* array of double constants 
				      * for user mat props 
                                      * User_Constants[species_no] is the
                                      * pointer to the data to be filled in. */
	       const int species_no) /* species number (zero if no species) */

{
  static char yo[] = "read_constants";
  char  line[255];
  char  *arguments[MAX_NUMBER_PARAMS];
  int num_const, i;
  int wspec=-1;
  
  if (species_no >= 0) wspec = species_no;


  if (species_no < 0) wspec = 0;
  if (User_constants == NULL)
    {
      fprintf(stderr, 
	      "%s:\t User Defined not allowed yet\n",yo);
      exit(-1);
    }
  if ( fgets(line, 255, imp) != NULL)
    {
      strip(line);
      /*
       * If there are any more parameters, allocate space
       * and parse remainder of the line for additional parameters.
       */
      if ((num_const = count_parameters(line)) > 0) {
	User_constants[wspec] = alloc_dbl_1(num_const, 0.0);
	
        tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);	
	for(i = 0; i < num_const; i++) {
	  User_constants[wspec][i] = atof(arguments[i]);
	}
      }
    }
  else 
    {
      fprintf(stderr, 
	      "%s:\tError @ Constants: USER \n",yo);
      exit(-1);
    }
  return num_const;
}
/* end of read_constants */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*
 * set_mp_to_unity -- set material properties structure for current material 
 *		      to unity.  But do nothing with constitutive equation 
 *		      parameters.
 *                    Initialize mp structure to default unity transport properties, 
 *                    unity source terms, and zero sensitivities
 * Comments:
 *                    This function is used to initialize all material properties before
 *                    the material property file is parsed.
 *
 * Created:		       Valentine's Day 13:25:03 MST 1995 prschun@sandia.gov
 * Revised:
 *
 * 1997/06/23 13:52 MDT pasacki@sandia.gov
 *
 *	Set the lengths of the user-defined constant lists to zero. Later,
 *	if various numbers of user-defined constants lists are scanned in, then
 *	these are set to the appropriate number.
 */

void
set_mp_to_unity(const int mn)
{
#ifdef DEBUG
  static const char yo[] = "set_mp_to_unity";
#endif
  int	i;
  int p;
  int v;
  int w;
  struct Elastic_Constitutive *e;
  MATRL_PROP_STRUCT  *m = mp_glob[mn];

  /*
   * First load up the simplest constitutive relations for all fluxes...
   *
   * Constant scalar coefficients linearly dependent on gradient of
   * one kind of varible...
   */

  pd_glob[mn]->MomentumFluxModel = NEWTONIAN;


  gn_glob[mn]->ConstitutiveEquation = NEWTONIAN;


  pd_glob[mn]->MassFluxModel = FICKIAN;

  cr_glob[mn]->HeatFluxModel	= CR_HF_FOURIER_0; 

  cr_glob[mn]->MeshFluxModel	= LINEAR;

  cr_glob[mn]->MomentumFluxModel	= CR_MF_NEWTON_0;
  
  cr_glob[mn]->MassFluxModel	= pd_glob[mn]->MassFluxModel;

  cr_glob[mn]->MeshMotion = pd_glob[mn]->MeshMotion;


  mp_glob[mn]->DefaultDatabase = DB_GOMA_MAT;

  mp_glob[mn]->thermal_conductivity = 1.;
  mp_glob[mn]->ConductivityModel = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      mp_glob[mn]->d_thermal_conductivity[v] = 0.;
    }

  gn_glob[mn]->mu0 = 1.;
  m->viscosity = 1.;
  m->ViscosityModel = CONSTANT;
  gn_glob[mn]->len_u_tau_y = 0;
  m->surface_tension = 1.;
  m->SurfaceTensionModel = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      m->d_viscosity[v] = 0.;
      m->d2_viscosity[v] = 0.;
      m->d_surface_tension[v] = 0.;
    }

  m->heat_capacity = 1.;
  m->HeatCapacityModel = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      m->d_heat_capacity[v] = 0.;
    }
  
  m->density = 1.;
  m->DensityModel = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      m->d_density[v] = 0.;
    }

  m->porosity = 1.;
  m->PorosityModel = CONSTANT;
  m->permeability = 1.;
  m->PermeabilityModel = CONSTANT;
  m->rel_liq_perm = 1.;
  m->RelLiqPermModel = CONSTANT;
  m->rel_gas_perm = 1.;
  m->RelGasPermModel = CONSTANT;
  m->saturation = 1.;
  m->SaturationModel = CONSTANT;
  m->matrix_density = 1.0;
  m->specific_heat = 1.0;


  for ( v=0; v<MAX_PMV + MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      m->d_porosity[v] = 0.;
      m->d_permeability[v] = 0.;
      m->d_rel_liq_perm[v] = 0.;
      m->d_rel_gas_perm[v] = 0.;
      m->d_saturation[v] = 0.;
    }

  mp_glob[mn]->Porous_wt_func = 0.;
  mp_glob[mn]->Porous_wt_funcModel=GALERKIN;
  mp_glob[mn]->PorousGasConstantsModel=CONSTANT;
  
  for ( w=0; w<pd_glob[mn]->Num_Porous_Eqn; w++)
    {
      mp_glob[mn]->PorousDiffusivityModel[w]=CONSTANT;
      mp_glob[mn]->porous_diffusivity[w] = 1.;
      mp_glob[mn]->PorousLatentHeatVapModel[w] = CONSTANT;
      mp_glob[mn]->porous_latent_heat_vap[w] = 1.;
      mp_glob[mn]->PorousLatentHeatFusionModel[w] = CONSTANT;
      mp_glob[mn]->porous_latent_heat_fusion[w] = 1.;
      mp_glob[mn]->PorousVaporPressureModel[w] = CONSTANT;
      mp_glob[mn]->porous_vapor_pressure[w] = 1.;
    }

  for ( w=0; w < MAX_PMV; w++) {
    for ( v = 0; v < MAX_CONC + MAX_VARIABLE_TYPES; v++) {
      mp_glob[mn]->d_porous_diffusivity[w][v] = 0.;
      mp_glob[mn]->d_porous_vapor_pressure[w][v] = 0.;
    }
  }

  mp_glob[mn]->Spwt_func = 0.;
  mp_glob[mn]->Spwt_funcModel=GALERKIN;
  
  for ( w=0; w<pd_glob[mn]->Num_Species_Eqn; w++)
    {
      m->DiffusivityModel[w]=CONSTANT;
      m->diffusivity[w] = 1.;
      m->VaporPressureModel[w] = CONSTANT;
      m->vapor_pressure[w] = 1.;
      for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
	{
	  m->d_diffusivity[w][v] = 0.;
	  m->d_vapor_pressure[w][v] = 0.;
	  m->d_molecular_weight[w][v] = 0.;
	  m->d_charge_number[w][v] = 0.;
	}
    }

  for ( w=0; w<MAX_CONC; w++)
    {
      m->GamDiffType[w]=CONSTANT;
      m->MuDiffType[w]=CONSTANT;
      m->FickDiffType[w]=CONSTANT;
      m->CurvDiffType[w]=CONSTANT;
      m->GravDiffType[w]=CONSTANT;
      m->NSCoeffType[w]=CONSTANT;
      m->QTensorDiffType[w] = NO_MODEL;
      m->gam_diffusivity[w]=0.;
      m->mu_diffusivity[w]=0.;
      m->g_diffusivity[w]=0.;
      m->cur_diffusivity[w]=0.;
    }

  for ( w=0; w<MAX_CONC; w++)
    {
      m->MolecularWeightModel[w]=CONSTANT;
      m->MolarVolumeModel[w]=CONSTANT;
      m->SpecificVolumeModel[w]=CONSTANT;
      m->ChargeNumberModel[w]=CONSTANT;
      m->RefConcnModel[w]=CONSTANT;
      m->molecular_weight[w] = -1.;
      m->molar_volume[w] = -1.;
      m->specific_volume[w] = -1.;
      m->reference_concn[w] = 0.;
      m->charge_number[w] = 0.;
      for ( v=0; v<MAX_CONC; v++)
	{
	  m->flory_param[w][v] = -1.;
	  m->diffusivity_Stefan_Maxwell[w][v] = 0.;
	  m->u_diffusivity_Stefan_Maxwell[w][v][0] = 0.;   /* for Dij - KSC */
	  m->u_diffusivity_Stefan_Maxwell[w][v][1] = 0.;   /* for E   - KSC */
	  m->u_diffusivity_Stefan_Maxwell[w][v][2] = 298.; /* for T0  - KSC */
	}
    }
  for (w = 0; w < MAX_VARIABLE_TYPES; w++) {
   m->StateVector[w] = 0.0;
  }

  elc_glob[mn]->lame_mu = 1.;
  elc_glob[mn]->lame_lambda = 1.;
  elc_glob[mn]->lame_mu_model = CONSTANT;
  elc_glob[mn]->lame_lambda_model = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      elc_glob[mn]->d_lame_mu[v] = 0.;
      elc_glob[mn]->d_lame_lambda[v] = 0.;
    }
  elc_glob[mn]->thermal_expansion = 0.;
  elc_glob[mn]->solid_reference_temp = 25.;
  elc_glob[mn]->thermal_expansion_model = CONSTANT;
  elc_glob[mn]->solid_reference_temp = CONSTANT;

  elc_rs_glob[mn]->lame_mu = 1.;
  elc_rs_glob[mn]->lame_lambda = 1.;
  elc_rs_glob[mn]->lame_mu_model = CONSTANT;
  elc_rs_glob[mn]->lame_lambda_model = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      elc_rs_glob[mn]->d_lame_mu[v] = 0.;
      elc_rs_glob[mn]->d_lame_lambda[v] = 0.;
    }
  elc_rs_glob[mn]->thermal_expansion = 0.;
  elc_rs_glob[mn]->solid_reference_temp = 25.;
  elc_rs_glob[mn]->thermal_expansion_model = CONSTANT;
  elc_rs_glob[mn]->solid_reference_temp = CONSTANT;



  for ( v=0; v<MAX_VARIABLE_TYPES; v++)
    {
      m->ReferenceModel[v] = CONSTANT;
      m->reference[v] = 0.;
    }

/*
 *****Finally, default all source terms to unity
 */

  for ( p=0; p<DIM; p++)
    {
      m->momentum_source[p] = 1.;
      m->heat_source = 1.;
      m->mass_source = 1.;
      m->mesh_source[p] = 1.;

      for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
	{
	  m->d_momentum_source[p][v] = 0.;
	  m->d_heat_source[v] = 0.;
	  m->d_mass_source[v] = 0.;
	  m->d_mesh_source[p][v] = 0.;
	}
    }


      
  for ( w=0; w<MAX_CONC; w++)
    {
      m->species_source[w] = 1.;
    }
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      m->d_species_source[v] = 0.;
    }

  /*
   * Initialize lengths of user defined constant lists to zero for this matl.
   */

  m->len_u_thermal_conductivity	        = 0;
  m->len_u_electrical_conductivity	= 0;
  m->len_u_viscosity			= 0;
  m->len_u_surface_tension		= 0;
  m->len_u_heat_capacity		= 0;
  m->len_u_Volume_Expansion		= 0;
  m->len_u_density			= 0;
  m->len_u_permeability		        = 0;
  m->len_u_rel_gas_perm		        = 0;
  m->len_u_rel_liq_perm		        = 0;
  m->len_u_saturation			= 0;
  m->len_u_momentum_source		= 0;
  m->len_u_heat_source			= 0;
  m->len_u_mass_source			= 0;
  m->len_u_mesh_source			= 0;
  m->len_u_current_source		= 0;

  for ( i=0; i<MAX_CONC; i++)
    {
      m->len_u_diffusivity[i]		= 0;
      m->len_u_gadiffusivity[i]		= 0;
      m->len_u_mdiffusivity[i]		= 0;
      m->len_u_cdiffusivity[i]		= 0;
      m->len_u_fdiffusivity[i]		= 0;
      m->len_u_gdiffusivity[i]		= 0;
      m->len_u_species_source[i]	= 0;
      m->len_u_species_vol_expansion[i] = 0;
      m->len_u_vapor_pressure[i]	= 0;
      m->len_u_reference_concn[i]	= 0;
    }


  /*
   * Elastic Constitutive properties...
   */

  e = elc_glob[mn];

  e->len_u_mu                          = 0;
  e->len_u_lambda                      = 0;
  e->len_u_v_mesh_sfs                  = 0;
  e->len_u_thermal_expansion           = 0;

}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* Routines for translating and applying command line arguments */

/*
 * usage -- give informative brief synopsis of command line options
 */

void 
usage(const int exit_flag)
{
  fprintf(stdout, 
	  "Usage: goma [-options] [file]\n\n");
  fprintf(stdout, 
	  "Where: [file] is name of overriding input file. If an\n");
  fprintf(stdout, 
	  "overriding input file is not specified, then \"input\" must\n");
  fprintf(stdout, 
	  "be readable in the current working directory.\n\n");
  fprintf(stdout, 
	  "Options:\n");
  fprintf(stdout, 
	  "\t-a [aargs], -aprepro [aargs]    Input thru APREPRO [w/ aargs].\n");
  fprintf(stdout, 
	  "\t-brk FILE                       Read Brk file from FILE\n");
  fprintf(stdout, 
	  "\t-restart FILE, -rest FILE       Read initial guess from FILE.\n");
  fprintf(stdout, 
	  "\t-d INT,     -debug INT          Set debug flag to INT.\n");
  fprintf(stdout, 
	  "\t-h,         -help               Print this message.\n");
  fprintf(stdout, 
	  "\t-i FILE,    -input FILE         Input from FILE.\n");
  fprintf(stdout, 
	  "\t-ix FILE,   -inexoII FILE       Read FEM from FILE.\n");
  fprintf(stdout, 
	  "\t-n INT,     -newton INT         Change # of Newton it'ns to INT.\n");
  fprintf(stdout, 
	  "\t-ts FLT,                        Start time of simulation (initial time).\n");
  fprintf(stdout, 
	  "\t-te FLT,                        End time of simulation (maximum time).\n");
  fprintf(stdout, 
	  "\t-nd,        -nodisplay          Send stdout, stderr to null.\n");
  fprintf(stdout, 
	  "\t-ox FILE,   -outexoII FILE      Write FEM to FILE.\n");
  fprintf(stdout, 
	  "\t-r FLT,     -relax FLT          Set Newton relaxation to FLT.\n");
  fprintf(stdout, 
	  "\t-s FILE,    -soln FILE          Write solution to FILE.\n");
  fprintf(stdout, 
	  "\t-se FILE,   -stderr FILE        Redirect stderr to FILE\n");
  fprintf(stdout, 
	  "\t-so FILE,   -stdout FILE        Redirect stdout to FILE\n");
  fprintf(stdout, 
	  "\t-cb FLT                         Continuation: Start value\n");
  fprintf(stdout, 
	  "\t-ce FLT                         Continuation: Final value\n");
  fprintf(stdout, 
	  "\t-cd FLT                         Continuation: Path step, ds\n");
  fprintf(stdout, 
	  "\t-cn INT                         Continuation: Max number of path steps\n");
  fprintf(stdout, 
	  "\t-cmin FLT                       Continuation: Minimum path step\n");
  fprintf(stdout, 
	  "\t-cmax FLT                       Continuation: Maximum path step\n");
  fprintf(stdout, 
	  "\t-cm INT                         Continuation: Method\n");
  fprintf(stdout, 
	  "\t-ct INT                         Continuation: Type\n");
  fprintf(stdout, 
	  "\t-c_bc INT                       Continuation: Boundary condition ID\n");
  fprintf(stdout, 
	  "\t-c_df INT                       Continuation: BC Data Float ID\n");
  fprintf(stdout, 
	  "\t-c_mn INT                       Continuation: Material ID\n");
  fprintf(stdout, 
	  "\t-c_mp INT                       Continuation: Method property ID\n");
  fprintf(stdout, 
	  "\t-bc_list                        List BC tags for continuation\n");
  fprintf(stdout, 
	  "\t-wr_int                         Turn Write Intermediate Results On\n");
  fprintf(stdout, 
	  "\t-time_pl INT                    read_exoII_file time plane (default last)\n");
  fprintf(stdout, 
	  "\t-v          --version           Print code version and exit\n");

  exit(exit_flag);
}

/* translate_command_line */
void
translate_command_line( int argc,
                        char *argv[],
                        struct Command_line_command **clc,
                        int *nclc )
     /* 
      * Routine which takes the data from the command line and c
      * allowed a type, an integer value, a double value, and a string.
      * Most commands will involve redirecting the input or output of the routine,
      * and may override what is specified in the 'input' file.
      * See below for the specific formats and explanations of commands
      *
      *  Written by: Richard Cairncross
      *              4/3/95
      *  Appended to by Ian Gates Aug 1998
      *  
      */
{
  char err_msg[MAX_CHAR_IN_INPUT];
  int err;
  int istr;
  char command_line_ap[MAX_COMMAND_LINE_LENGTH];

  char temp_file_inp[MAX_FNL]=".tmp.input"; /* Hold APREPRO'd input file. */
  static char temp_file_out[]="tmp.std.output";
  static char temp_file_err[]="tmp.std.errput";

  static const char yo[] = "translate_command_line";

/* move through command line starting with second entry and separate into commands
 * and values */
  istr = 1;
  *nclc = -1;
  while(istr < argc)
    {
      if (*argv[istr] == '-')
	{
/* 
 * OPTION -nodisplay: redirect stdout and stderr to a temporary file 
 * that is destroyed at end of successful runs 
 */
	  if (strcmp(argv[istr], "-nodisplay") == 0 || 
	      strcmp(argv[istr], "-nd") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = NO_DISPLAY;
	      strcpy_rtn = strcpy(clc[*nclc]->string, temp_file_out);
	      (*nclc)++;
	      clc[*nclc]->type = NO_DISPLAY;
	      strcpy_rtn = strcpy(clc[*nclc]->string, temp_file_err);
	      /* create temporary file and redirect the output */
	      if (freopen(temp_file_out,"w",stdout) == NULL)
		EH(-1, "Problem redirecting stdout stream");
	      if (freopen(temp_file_err,"w",stderr) == NULL)
		EH(-1, "Problem redirecting stderr stream");
	    }
/* 
 * OPTION -h, -help:  print informative usage
 */
	  else if (strcmp(argv[istr], "-h") == 0 || 
		   strcmp(argv[istr], "-help") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = HELP_USAGE;
	      usage(0);
	    }
/* 
 * OPTION -v, --version:  print code version and exit
 */
	  else if (strcmp(argv[istr], "-v") == 0 || 
		   strcmp(argv[istr], "--version") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = PRINT_CODE_VERSION;
	      print_code_version(); /* don't come back */
	      log_msg("This is GOMA version %s", GOMA_VERSION);
	      log_msg("GOMA ends normally.");
	      echo_compiler_settings();
	      exit(0);
	    }
/*
 * OPTION -chdir: change working directory before starting
 */
	  else if(strcmp(argv[istr], "-chdir") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = WORKING_DIRECTORY;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	      fprintf(stderr,"Attempting to change working directory to %s\n",
               clc[*nclc]->string);
	      err = chdir( clc[*nclc]->string );\
          EH(err, "Could not change directory! Exiting..");
	    }
/* 
 * OPTION -stderr: redirect stderr to a file
 *         Print a warning message about expected behavior from MP
 *         programs.
 */
	  else if (strcmp(argv[istr], "-se") == 0 ||
		   strcmp(argv[istr], "-stderr") == 0 )
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = STDERR_OUT;
	      /* create temporary file and redirect the output */
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	      if (freopen(clc[*nclc]->string, "w", stderr) == NULL)
		{
		  sr = sprintf(err_msg, "Problem diverting stderr to \"%s\"",
			       clc[*nclc]->string);
		  EH(-1, err_msg);
		}
	      if (Num_Proc > 1 && Unlimited_Output ) {
                fprintf(stderr,
			"GOMA WARNING: Standard error was redirected to %s on Proc 0 only!\n",
			clc[*nclc]->string);
                fprintf(stderr,
			"\t\tThe other processors will still print to the terminal\n");
		fprintf(stderr,
			"\t\tA better method to redirect standard error would be to use shell redirects\n");
	      }
	    }
/* 
 * OPTION -stdout: redirect stdout to a file
 *         Print a warning message about expected behavior from MP
 *         programs. 
 */
	  else if (strcmp(argv[istr], "-so") == 0 ||
		   strcmp(argv[istr], "-stdout") == 0 )
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = STDOUT_OUT;
	      /* create temporary file and redirect the output */
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	      if (freopen(clc[*nclc]->string,"w",stdout) == NULL)
		{
		  sr = sprintf(err_msg, "Problem diverting stdout to \"%s\"",
			       clc[*nclc]->string);
		  EH(-1, err_msg);
		}
	      if (Num_Proc > 1) {
                fprintf(stderr,
			"GOMA WARNING: Standard output was redirected to %s on Proc 0 only!\n",
			clc[*nclc]->string);
                fprintf(stderr,
			"\t\tThe other processors will still print to the terminal\n");
		fprintf(stderr,
			"\t\tA better method to redirect standard output would be to use shell redirects\n");
	      }	      
	    }
/* 
 * OPTION -input: read input info from alternate file instead of 'input'
 */
	  else if(strcmp(argv[istr], "-input") == 0 || 
		  strcmp(argv[istr], "-i") == 0)
	    {
	/* read input info from alternate file instead of 'input' */
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = INPUT_FILE;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
#ifdef DEBUG
	      fprintf(stdout, "Change reading input from \"%s\" to \"%s\"\n",
		      Input_File, argv[istr]); 
#endif
	      strcpy_rtn = strcpy(Input_File, clc[*nclc]->string);
	      istr++;
	    }
/* 
 * OPTION -contin: read continuation (guessfile) from alternate file overriding file in 'input'
 */
	  else if(strcmp(argv[istr], "-contin") == 0 || strcmp(argv[istr], "-c") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONTIN_FILE;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -soln: write solution to alternate file overriding file in 'input'
 */
	  else if(strcmp(argv[istr], "-soln") == 0 || strcmp(argv[istr], "-s") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      clc[*nclc]->type = SOLN_FILE;
	      istr++;
	    }
/* 
 * OPTION -inexoII: read FEM file from alternate file overriding file in 'input'
 */
	  else if(strcmp(argv[istr], "-inexoII") == 0 || strcmp(argv[istr], "-ix") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = INEXOII_FILE;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -outexoII: read FEM file from alternate file overriding file in 'input'
 */
	  else if(strcmp(argv[istr], "-outexoII") == 0 || strcmp(argv[istr], "-ox") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = OUTEXOII_FILE;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -debug: override debug flag in 'input'
 */
	  else if(strcmp(argv[istr], "-debug") == 0 || strcmp(argv[istr], "-d") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = DEBUG_OPTION;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -newton: override number of newton iterations in 'input'
 */
	  else if(strcmp(argv[istr], "-newton") == 0 || strcmp(argv[istr], "-n") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = NEWTON_NUM;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -ts: override start time of the simulation (initial time)
 */
	  else if(strcmp(argv[istr], "-ts") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = TIME_START;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -te: override end time of the simulation (maximum time)
 */
	  else if(strcmp(argv[istr], "-te") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = TIME_END;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/*
 * OPTION -relax: override newton update factor in 'input'
 */
	  else if(strcmp(argv[istr], "-relax") == 0 || strcmp(argv[istr], "-r") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = RELAXATION;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -aprepro: convert input file with aprepro prior to reading it
 */
	  else if(strcmp(argv[istr], "-aprepro") == 0 || 
		  strcmp(argv[istr], "-a") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = APREPRO;
	      run_aprepro = 1;

	      /* check for aprepro options */
	      sr = sprintf(command_line_ap, "aprepro ");
	      /*
	       * assume all entries before next goma option (-[option_name]) 
	       * are arguments for APREPRO, and that --[option_name] are 
	       * options for aprepro 
	       */

	      while( istr < argc && (*argv[istr] != '-' || *(argv[istr] + 1) == '-'  ))
		{
		  if (*(argv[istr] + 1) == '-') 
		    {
		      strcat(command_line_ap, argv[istr++] + 1); /* remove first '-' */
		    }
		  else
		    {
		      strcat(command_line_ap, argv[istr++]);
		    }
		  strcat(command_line_ap," ");
		}
	      sprintf(aprepro_command, "%s", command_line_ap);
	      strcpy_rtn = strcpy(clc[*nclc]->string, command_line_ap);
	    } /*end of else if list */

/* 
 * OPTION -restart: read restart (guessfile) from alternate file overriding file in 'input'
 */
	  else if(strcmp(argv[istr], "-restart") == 0 || strcmp(argv[istr], "-rest") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONTIN_FILE;
	      strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -cb: CONTINUATION BEGIN VALUE
 */
	  else if(strcmp(argv[istr], "-cb") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_BEG_PVALUE;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -ce: CONTINUATION END VALUE
 */
	  else if(strcmp(argv[istr], "-ce") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_END_PVALUE;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -cd: CONTINUATION STEP
 */
	  else if(strcmp(argv[istr], "-cd") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_PATHSTEP;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -cn: CONTINUATION MAX PATH STEPS
 */
	  else if(strcmp(argv[istr], "-cn") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_PATH_STEPS;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -cmin: CONTINUATION MINIMUM STEP SIZE VALUE
 */
	  else if(strcmp(argv[istr], "-cmin") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_MIN_PVALUE;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -cmax: CONTINUATION MAXIMUM STEP SIZE VALUE
 */
	  else if(strcmp(argv[istr], "-cmax") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_MAX_PVALUE;
	      clc[*nclc]->r_val = strtod(argv[istr], NULL);
	      istr++;
	    }
/* 
 * OPTION -cm: CONTINUATION METHOD
 */
	  else if(strcmp(argv[istr], "-cm") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_METHOD;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -ct: CONTINUATION TYPE
 */
	  else if(strcmp(argv[istr], "-ct") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_TYPE;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -c_bc: CONTINUATION BOUNDARY CONDITION ID
 */
	  else if(strcmp(argv[istr], "-c_bc") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_BCID;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -c_df: CONTINUATION DATA FLOAT ID
 */
	  else if(strcmp(argv[istr], "-c_df") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_DFID;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -c_mn: CONTINUATION MATERIAL ID
 */
	  else if(strcmp(argv[istr], "-c_mn") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_MTID;
	      clc[*nclc]->i_val = atoi(argv[istr])-1;
	      istr++;
	    }
	    
/* 
 * OPTION -newp: USE THE FLEX/BISON PARSER 
 */
 	  else if(strcmp(argv[istr], "-newp") == 0)	
	    {	   					
	      (*nclc)++;				
	      istr++;					
	      clc[*nclc]->type = PARSER_OPTION;		
	      New_Parser_Flag = 1;			
	    } 					
	 	    	    
/* 
 * OPTION -c_mp: CONTINUATION MATERIAL PROPERTY ID
 */
	  else if(strcmp(argv[istr], "-c_mp") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_MPID;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/* 
 * OPTION -bc_list: LIST BOUNDARY CONDITIONS AND STOP
 */
	  else if(strcmp(argv[istr], "-bc_list") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = CONT_BC_LIST;
	    }
/* 
 * OPTION -wr_int: TURN INTERMEDIATE RESULTS ON
 */
	  else if(strcmp(argv[istr], "-wr_int") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = WRITE_INTERMEDIATE;
	    }
/* 
 * OPTION -time_pl: SPECIFY EXOII FILE STEP NUMBER TO READ
 */
	  else if(strcmp(argv[istr], "-time_pl") == 0)
	    {
	      (*nclc)++;
	      istr++;
	      clc[*nclc]->type = EXOII_TIME_PLANE;
	      clc[*nclc]->i_val = atoi(argv[istr]);
	      istr++;
	    }
/*
 * OPTION -ne:  Disable file echoing 
 */
 
		else if( strcmp( argv[istr],"-ne") == 0 )
		{
		  (*nclc)++;
	      istr++;
	      clc[*nclc]->type = NOECHO;
		  ECHO("NOECHO", NULL) ;
		}	
		else if( strcmp( argv[istr],"-brk") == 0 )
		{
                  Brk_Flag = 1;
		  (*nclc)++;
                  istr++;
                  clc[*nclc]->type = NOECHO;
		  strcpy_rtn = strcpy(clc[*nclc]->string, argv[istr]);
		  strcpy_rtn = strcpy( Brk_File, clc[*nclc]->string);
		  istr++;
		}
/*
 * Unknown '-' option: print an error and abort
 */
	  else {
	    sprintf(err_msg, "ERROR EXIT: unknown dash option: %s\n",
		    argv[istr]);
	    EH(-1, err_msg);
	    ABORTH(-1, "IMMEDIATE PARALLEL EXIT - can't recover gracefully");
	  }
	} /* end of if argv starts with '-' */
/* 
 * DEFAULT OPTION -input: read input info from alternate file instead of 'input'
 */
      else
	{
	  if (istr == 1 && istr == argc - 1) 
	    /* default assumes this string is the name of the input file */
	    {
	      clc[0]->type = INPUT_FILE;
	      strcpy_rtn = strcpy(clc[0]->string, argv[istr]);
	      fprintf(stdout, "Change reading input from \"%s\" to \"%s\"\n",
		      Input_File, argv[istr]); 

	      strcpy_rtn = strcpy(Input_File, clc[0]->string);
	      istr++;
	      *nclc = 0;
	    }
	  else 
	    {
	      usage(-1);
	    }
	}
    }
  (*nclc)++;
  if ( Debug_Flag > 0 )
    {
      fprintf(stdout, "Total of %d commands found\n", *nclc);
    }

  /* make Echo_Input_File consistent with Input_File */

  if(strcmp(Input_File,"input") != 0 )
    {
      strcpy(Echo_Input_File, "echo_");
      strcat(Echo_Input_File, Input_File);
    }
/* Perform system calls to aprepro, fastq, etc */

  if (run_aprepro == 1) 
    {
      strcpy(temp_file_inp, "tmp.");
      strcat(temp_file_inp, Input_File);

      strcat(command_line_ap, Input_File);
      strcat(command_line_ap, " ");
      strcat(command_line_ap, temp_file_inp);
      if ( Debug_Flag > 0 )
	{
	  fprintf(stdout, "system: %s\n", command_line_ap);
	}
#ifdef DEBUG
      fprintf(stderr, "system() = \"%s\"\n", command_line_ap);
#endif

#ifndef tflop
      err = system(command_line_ap);
      EH(err, "system() choked on input file.");

      if (WEXITSTATUS(err) == 127)
	{
	  EH(-1, "System call failed, aprepro not found");
	  return;
	}
#else
      EH(-1, "aprepro the input file prior to running goma.");
#endif

      strcpy(Input_File, temp_file_inp);
      if ( Debug_Flag > 0 )
	{
	  fprintf(stdout, "reading input from \"%s\"\n", Input_File);
	}
    }
  return;
}
/* END of routine translate_command_line */
/*****************************************************************************/

/* apply_command_line */
void
apply_command_line(struct Command_line_command **clc,
		   int nclc)
/*
 *  This routine overrides entries in input file with command-line options
 *
 *  Written by: Richard Cairncross
 *              4/4/95
 */
{
  static char a[] = "Applying command line override:  ";
  static char b[] = "                               ->";
  int i;

  for (i=0; i<nclc; i++)
    {
       if (clc[i]->type == CONTIN_FILE) {
	 fprintf(stdout, "%s%40s= %s\n%s%40s= %s\n",
		 a, "GUESS file", Init_GuessFile, 
		 b, "GUESS file", clc[i]->string); 
	 strcpy_rtn = strcpy(Init_GuessFile, clc[i]->string);
       }
       else if (clc[i]->type == HELP_USAGE) {
	 usage(0);
       }
       else if (clc[i]->type == SOLN_FILE) {
	 fprintf(stdout, "%s%40s= %s\n%s%40s= %s\n",
		 a, "SOLN file", Soln_OutFile,
		 b, "SOLN file", clc[i]->string); 
	 strcpy_rtn = strcpy(Soln_OutFile, clc[i]->string);
       }
       else if (clc[i]->type == INEXOII_FILE) {
	 fprintf(stdout, "%s%40s= %s\n%s%40s= %s\n",
		 a, "FEM file", ExoFile,
		 b, "FEM file", clc[i]->string); 
	 strcpy_rtn = strcpy(ExoFile, clc[i]->string);
       }
       else if (clc[i]->type == OUTEXOII_FILE) {
	 fprintf(stdout, "%s%40s= %s\n%s%40s= %s\n",
		 a, "Output EXODUS II file", ExoFileOut,
		 b, "Output EXODUS II file", clc[i]->string); 
	 strcpy_rtn = strcpy(ExoFileOut, clc[i]->string);
         strcpy_rtn = strcpy(ExoFileOutMono, clc[i]->string);
       }
       else if (clc[i]->type == DEBUG_OPTION) {
	 fprintf(stdout, "%s%40s= %d\n%s%40s= %d\n",
		 a, "Debug", Debug_Flag,
		 b, "Debug", clc[i]->i_val);
	 Debug_Flag = clc[i]->i_val;
       }
       else if (clc[i]->type == NEWTON_NUM) {
	 fprintf(stdout, "%s%40s= %d\n%s%40s= %d\n",
		 a, "Number of Newton Iterations", Max_Newton_Steps,
		 b, "Number of Newton Iterations", clc[i]->i_val);
	 Max_Newton_Steps = clc[i]->i_val;
       }
       else if (clc[i]->type == TIME_START) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Initial Time", tran->init_time,
		 b, "Initial Time", clc[i]->r_val);
	 tran->init_time = clc[i]->r_val;
       }
       else if (clc[i]->type == TIME_END) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Maximum Time", tran->TimeMax,
		 b, "Maximum TIme", clc[i]->r_val);
	 tran->TimeMax = clc[i]->r_val;
       }
       else if (clc[i]->type == RELAXATION) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Newton correction factor", damp_factor1,
		 b, "Newton correction factor", clc[i]->r_val);
	 damp_factor1 = clc[i]->r_val;
         damp_factor2 = -1.;
	 damp_factor3 = -1.;  /* To force the override default*/
       }
       else if (clc[i]->type == CONT_BEG_PVALUE) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Initial parameter value", cont->BegParameterValue,
		 b, "Initial parameter value", clc[i]->r_val);
	 cont->BegParameterValue = clc[i]->r_val;
       }
       else if (clc[i]->type == CONT_END_PVALUE) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Final parameter value", cont->EndParameterValue,
		 b, "Final parameter value", clc[i]->r_val);
	 cont->EndParameterValue = clc[i]->r_val;
       }
       else if (clc[i]->type == CONT_PATHSTEP) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "delta_s", cont->Delta_s0,
		 b, "delta_s", clc[i]->r_val);
	 cont->Delta_s0 = clc[i]->r_val;
       }
       else if (clc[i]->type == CONT_PATH_STEPS) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Maximum number of path steps", cont->MaxPathSteps,
		 b, "Maximum number of path steps", clc[i]->i_val);
	 cont->MaxPathSteps = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_MIN_PVALUE) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Minimum path step", cont->Delta_s_min,
		 b, "Minimum path step", clc[i]->r_val);
	 cont->Delta_s_min = clc[i]->r_val;
       }
       else if (clc[i]->type == CONT_MAX_PVALUE) {
	 fprintf(stdout, "%s%40s= %g ...\n%s%40s= %g ...\n",
		 a, "Maximum path step", cont->Delta_s_max,
		 b, "Maximum path step", clc[i]->r_val);
	 cont->Delta_s_max = clc[i]->r_val;
       }
       else if (clc[i]->type == CONT_METHOD) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Continuation", Continuation,
		 b, "Continuation", clc[i]->i_val);
	 Continuation = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_TYPE) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Continuation Type", cont->upType,
		 b, "Continuation Type", clc[i]->i_val);
	 cont->upType = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_BCID) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Boundary condition ID", cont->upBCID,
		 b, "Boundary condition ID", clc[i]->i_val);
	 cont->upBCID = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_DFID) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Boundary condition data float tag", cont->upDFID,
		 b, "Boundary condition data float tag", clc[i]->i_val);
	 cont->upDFID = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_MTID) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Material id", cont->upMTID,
		 b, "Material id", clc[i]->i_val);
	 cont->upMTID = clc[i]->i_val;
       }
       else if (clc[i]->type == CONT_MPID) {
	 fprintf(stdout, "%s%40s= %d ...\n%s%40s= %d ...\n",
		 a, "Material property tag", cont->upMPID,
		 b, "Material property tag", clc[i]->i_val);
	 cont->upMPID = clc[i]->i_val;
       }
       else if (clc[i]->type == PARSER_OPTION) {	
	 fprintf(stdout,"-newp request.");		
	 New_Parser_Flag = 1;				
       }       						     
       else if (clc[i]->type == CONT_BC_LIST) {
	 fprintf(stdout,"-bc_list request.\n\n\t goma done.\n\n"); 
	 exit(0);
       }
       else if (clc[i]->type == WRITE_INTERMEDIATE) {
	 fprintf(stdout,"Write Intermediate Solutions request.\n\n"); 
         Write_Intermediate_Solutions = TRUE;
       }
       else if (clc[i]->type == EXOII_TIME_PLANE) {
	 fprintf(stdout,"Exodus Time Plane = %d\n\n",clc[i]->i_val); 
	 ExoTimePlane = clc[i]->i_val;
       }
     }
  return;
}
/* end of apply_command_line */

struct Data_Table *
setup_table_BC(FILE *ifp,
			   char *input,
			   struct Boundary_Condition *BC_Type,
			   char *echo_string)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  /*           Function allocates memory for the BC table structure.
	       Parses the remainder of the TABLE BC card identifying
	       the abscissa and ordinate variables and setting the
	       members of the table structure and the BC_desc structure
	       accordingly.  Opens a separate file for reading the 
	       tabular data if indicated on the card and calls a function
	       to read the tabular data.

	       Author: Thomas A. Baer, Org. 9111
	       Date  : Jly 15, 1998
	       

	Parameters:
	       ifp     =  pointer to input deck FILE structure
	       input   =  character buffer array 
	       BC_Type =  pointer to TABLE BC structure that the table is being added to

        returns:
	       pointer to table structure created in the function
   */


  char *yo, *dataname = NULL;
  FILE *datafile =NULL;
  char *es;
  char * echo_file = Echo_Input_File;

  yo = "setup_table_BC";
  
  es = endofstring(echo_string);

  BC_Type->table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) ) ;

  if( BC_Type->BC_Name == TABLE_WICV_BC )
                {BC_Type->table->columns = 3;}
  else { BC_Type->table->columns = 2;}

  if ( fscanf(ifp, "%80s", input ) != 1 )
    {
      	sprintf (err_msg, "%s:\tError reading TABLE BC \n", yo);
	EH(-1,err_msg);
    }
  
  strip( input );
  stringup( input );

  /* 
   * "x-axis" of table
   */

  if ( strcmp( input, "TIME" ) == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"TIME");
      BC_Type->table->t_index[0] = -1;
      /*BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;      */
    }
  else if( strcmp( input, "X") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"X");
      BC_Type->table->t_index[0] = 0;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
    }
  else if ( strcmp( input, "Y") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Y");
      BC_Type->table->t_index[0] = 1;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
    }
  else if ( strcmp( input, "Z") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Z");
      BC_Type->table->t_index[0] = 2;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
    }
  else if ( strcmp( input, "XY") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"X");
      BC_Type->table->t_index[0] = 0;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
      strcpy( BC_Type->table->t_name[1],"Y");
      BC_Type->table->t_index[1] = 1;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
    }
  else if ( strcmp( input, "XZ") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"X");
      BC_Type->table->t_index[0] = 0;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
      strcpy( BC_Type->table->t_name[1],"Z");
      BC_Type->table->t_index[1] = 2;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
    }
  else if ( strcmp( input, "YZ") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Y");
      BC_Type->table->t_index[0] = 1;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
      strcpy( BC_Type->table->t_name[1],"Z");
      BC_Type->table->t_index[1] = 2;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
    }
  else if ( strcmp( input, "YX") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Y");
      BC_Type->table->t_index[0] = 1;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
      strcpy( BC_Type->table->t_name[1],"X");
      BC_Type->table->t_index[1] = 0;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
    }
  else if ( strcmp( input, "ZX") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Z");
      BC_Type->table->t_index[0] = 2;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
      strcpy( BC_Type->table->t_name[1],"X");
      BC_Type->table->t_index[1] = 0;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
    }
  else if ( strcmp( input, "ZY") == 0 )
    {
      strcpy( BC_Type->table->t_name[0],"Z");
      BC_Type->table->t_index[0] = 2;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
      strcpy( BC_Type->table->t_name[1],"Y");
      BC_Type->table->t_index[1] = 1;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
    }
  else
    {
      sprintf(err_msg,"\nInvalid choice for table abscissa.");
      EH(-1,err_msg);
    }

  SPF(es," %s", input );  

  /*
   * "y-axis" of table
   */

  if ( fscanf(ifp, "%80s", input ) != 1 )
    {
      	sprintf (err_msg, "%s:\tError reading TABLE BC \n", yo);
	EH(-1,err_msg);
    }
  
  strip( input );
  stringup( input );

  if ( ( strcmp( input, "VELOCITY1") == 0 ) || ( strcmp( input, "U") == 0 ) )
    {
      BC_Type->table->f_name = "VELOCITY1";
      BC_Type->table->f_index = VELOCITY1;
      BC_Type->desc->equation = R_MOMENTUM1;
      BC_Type->desc->sens[VELOCITY1] = 1;
    }
  else if ( ( strcmp( input, "VELOCITY2") == 0 ) || ( strcmp( input, "V") == 0 ) )
    {
      BC_Type->table->f_name = "VELOCITY2";
      BC_Type->table->f_index = VELOCITY2;
      BC_Type->desc->equation = R_MOMENTUM2;
      BC_Type->desc->sens[VELOCITY2] = 1;
    }
  else if ( ( strcmp( input, "VELOCITY3") == 0 ) || ( strcmp( input, "W") == 0 ) )
    {
      BC_Type->table->f_name = "VELOCITY3";
      BC_Type->table->f_index = VELOCITY3;
      BC_Type->desc->equation = R_MOMENTUM3;
      BC_Type->desc->sens[VELOCITY3] = 1;
    }
  else if ( ( strcmp( input, "TEMPERATURE") == 0 ) || ( strcmp( input, "T") == 0 ) )
    {
      BC_Type->table->f_name = "TEMPERATURE";
      BC_Type->table->f_index = TEMPERATURE;
      BC_Type->desc->equation = R_ENERGY;
      BC_Type->desc->sens[TEMPERATURE] = 1;
    }
  else if ( ( strcmp( input, "MASS_FRACTION") == 0 ) || ( strcmp( input, "Y") == 0 ) || 
	    ( strcmp( input, "SPECIES") == 0 ))
    {
      BC_Type->table->f_name = "MASS_FRACTION";
      BC_Type->table->f_index = MASS_FRACTION;
      BC_Type->desc->equation = R_MASS;
      BC_Type->desc->sens[MASS_FRACTION] = 1;
      if ( fscanf( ifp, "%d", &BC_Type->species_eq) != 1)
	{
	  sprintf (err_msg, "%s:\tError reading species number on TABLE BC \n", yo);
	  EH(-1,err_msg);
	}
    }

  else if ( ( strcmp( input, "MESH_DISPLACEMENT1") == 0 ) || ( strcmp( input, "DX") == 0 ) )
    {
      BC_Type->table->f_name = "MESH_DISPLACEMENT1";
      BC_Type->table->f_index = MESH_DISPLACEMENT1;
      BC_Type->desc->equation = R_MESH1;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
    }
  else if ( ( strcmp( input, "MESH_DISPLACEMENT2") == 0 ) || ( strcmp( input, "DY") == 0 ) )
    {
      BC_Type->table->f_name = "MESH_DISPLACEMENT2";
      BC_Type->table->f_index = MESH_DISPLACEMENT2;
      BC_Type->desc->equation = R_MESH2;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
    }
  else if ( ( strcmp( input, "MESH_DISPLACEMENT3") == 0 ) || ( strcmp( input, "DZ") == 0 ) )
    {
      BC_Type->table->f_name = "MESH_DISPLACEMENT3";
      BC_Type->table->f_index = MESH_DISPLACEMENT3;
      BC_Type->desc->equation = R_MESH3;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
    }
  else if ( ( strcmp( input, "MESH_POSITION1") == 0 )  )
    {
      BC_Type->table->f_name = "MESH_POSITION1";
      BC_Type->table->f_index = MESH_POSITION1;
      BC_Type->desc->equation = R_MESH1;
      BC_Type->desc->sens[MESH_DISPLACEMENT1] = 1;
    }
  else if ( ( strcmp( input, "MESH_POSITION2") == 0 ) )
    {
      BC_Type->table->f_name = "MESH_POSITION2";
      BC_Type->table->f_index = MESH_POSITION2;
      BC_Type->desc->equation = R_MESH2;
      BC_Type->desc->sens[MESH_DISPLACEMENT2] = 1;
    }
  else if ( ( strcmp( input, "MESH_POSITION3") == 0 ) )
    {
      BC_Type->table->f_name = "MESH_POSITION3";
      BC_Type->table->f_index = MESH_POSITION3;
      BC_Type->desc->equation = R_MESH3;
      BC_Type->desc->sens[MESH_DISPLACEMENT3] = 1;
    }
  else if ( ( strcmp( input, "SOLID_DISPLACEMENT1") == 0 ) || ( strcmp( input, "DX_RS") == 0 ) )
    {
      BC_Type->table->f_name = "SOLID_DISPLACEMENT1";
      BC_Type->table->f_index = SOLID_DISPLACEMENT1;
      BC_Type->desc->equation = R_SOLID1;
      BC_Type->desc->sens[SOLID_DISPLACEMENT1] = 1;
    }
  else if ( ( strcmp( input, "SOLID_DISPLACEMENT2") == 0 ) || ( strcmp( input, "DY_RS") == 0 ) )
    {
      BC_Type->table->f_name = "SOLID_DISPLACEMENT2";
      BC_Type->table->f_index = SOLID_DISPLACEMENT2;
      BC_Type->desc->equation = R_SOLID2;
      BC_Type->desc->sens[SOLID_DISPLACEMENT2] = 1;
    }
  else if ( ( strcmp( input, "SOLID_DISPLACEMENT3") == 0 ) || ( strcmp( input, "DZ_RS") == 0 ) )
    {
      BC_Type->table->f_name = "SOLID_DISPLACEMENT3";
      BC_Type->table->f_index = SOLID_DISPLACEMENT3;
      BC_Type->desc->equation = R_SOLID3;
      BC_Type->desc->sens[SOLID_DISPLACEMENT3] = 1;
    }
  else if ( ( strcmp( input, "PRESSURE") == 0 ) || ( strcmp( input, "P") == 0 ) )
    {
      BC_Type->table->f_name = "PRESSURE";
      BC_Type->table->f_index = PRESSURE;
      BC_Type->desc->equation = R_PRESSURE;
      BC_Type->desc->sens[PRESSURE] = 1;
    }
  else if ( ( strcmp( input, "SHEAR_RATE") == 0 ) || ( strcmp( input, "SH") == 0 ) )
    {
      BC_Type->table->f_name = "SHEAR_RATE";
      BC_Type->table->f_index = SHEAR_RATE;
      BC_Type->desc->equation = R_SHEAR_RATE;
      BC_Type->desc->sens[SHEAR_RATE] = 1;
    }
  else if ( ( strcmp( input, "EXT_VELOCITY") == 0 ) || ( strcmp( input, "EXT_V") == 0 ) )
    {
      BC_Type->table->f_name = "EXT_VELOCITY";
      BC_Type->table->f_index = EXT_VELOCITY;
      BC_Type->desc->equation = R_EXT_VELOCITY;
      BC_Type->desc->sens[EXT_VELOCITY] = 1;
    }
  else if ( ( strcmp( input, "EFIELD1") == 0 ) || ( strcmp( input, "E1") == 0 ) )
    {
      BC_Type->table->f_name = "EFIELD1";
      BC_Type->table->f_index = EFIELD1;
      BC_Type->desc->equation = R_EFIELD1;
      BC_Type->desc->sens[EFIELD1] = 1;
    }
  else if ( ( strcmp( input, "EFIELD2") == 0 ) || ( strcmp( input, "E2") == 0 ) )
    {
      BC_Type->table->f_name = "EFIELD2";
      BC_Type->table->f_index = EFIELD2;
      BC_Type->desc->equation = R_EFIELD2;
      BC_Type->desc->sens[EFIELD2] = 1;
    }
  else if ( ( strcmp( input, "EFIELD3") == 0 ) || ( strcmp( input, "E3") == 0 ) )
    {
      BC_Type->table->f_name = "EFIELD3";
      BC_Type->table->f_index = EFIELD3;
      BC_Type->desc->equation = R_EFIELD3;
      BC_Type->desc->sens[EFIELD3] = 1;
    }
  else if ( ( strcmp( input, "ENORM") == 0 ) )
    {
      BC_Type->table->f_name = "ENORM";
      BC_Type->table->f_index = ENORM;
      BC_Type->desc->equation = R_ENORM;
      BC_Type->desc->sens[ENORM] = 1;
    }
  else if ( ( strcmp( input, "VOLTAGE") == 0 ) || ( strcmp( input, "V") == 0 ) )
    {
      BC_Type->table->f_name = "VOLTAGE";
      BC_Type->table->f_index = VOLTAGE;
      BC_Type->desc->equation = R_POTENTIAL;
      BC_Type->desc->sens[VOLTAGE] = 1;
    }
  else if ( ( strcmp( input, "S11") == 0 ) )
    {
      BC_Type->table->f_name = "S11";
      BC_Type->table->f_index = POLYMER_STRESS11;
      BC_Type->desc->equation = R_STRESS11;
      BC_Type->desc->sens[POLYMER_STRESS11] = 1;
    }
  else if ( ( strcmp( input, "S12") == 0 ) )
    {
      BC_Type->table->f_name = "S12";
      BC_Type->table->f_index = POLYMER_STRESS12;
      BC_Type->desc->equation = R_STRESS12;
      BC_Type->desc->sens[POLYMER_STRESS12] = 1;
    }
  else if ( ( strcmp( input, "S22") == 0 ) )
    {
      BC_Type->table->f_name = "S22";
      BC_Type->table->f_index = POLYMER_STRESS22;
      BC_Type->desc->equation = R_STRESS22;
      BC_Type->desc->sens[POLYMER_STRESS22] = 1;
    }
  else if ( ( strcmp( input, "S13") == 0 ) )
    {
      BC_Type->table->f_name = "S13";
      BC_Type->table->f_index = POLYMER_STRESS13;
      BC_Type->desc->equation = R_STRESS13;
      BC_Type->desc->sens[POLYMER_STRESS13] = 1;
    }
  else if ( ( strcmp( input, "S23") == 0 ) )
    {
      BC_Type->table->f_name = "S23";
      BC_Type->table->f_index = POLYMER_STRESS23;
      BC_Type->desc->equation = R_STRESS23;
      BC_Type->desc->sens[POLYMER_STRESS23] = 1;
    }
  else if ( ( strcmp( input, "S33") == 0 ) )
    {
      BC_Type->table->f_name = "S33";
      BC_Type->table->f_index = POLYMER_STRESS33;
      BC_Type->desc->equation = R_STRESS33;
      BC_Type->desc->sens[POLYMER_STRESS33] = 1;
    }
  else if ( ( strcmp( input, "S11_1") == 0 ) )
    {
      BC_Type->table->f_name = "S11_1";
      BC_Type->table->f_index = POLYMER_STRESS11_1;
      BC_Type->desc->equation = R_STRESS11_1;
      BC_Type->desc->sens[POLYMER_STRESS11_1] = 1;
    }
  else if ( ( strcmp( input, "S12_1") == 0 ) )
    {
      BC_Type->table->f_name = "S12_1";
      BC_Type->table->f_index = POLYMER_STRESS12_1;
      BC_Type->desc->equation = R_STRESS12_1;
      BC_Type->desc->sens[POLYMER_STRESS12_1] = 1;
    }
  else if ( ( strcmp( input, "S22_1") == 0 ) )
    {
      BC_Type->table->f_name = "S22_1";
      BC_Type->table->f_index = POLYMER_STRESS22_1;
      BC_Type->desc->equation = R_STRESS22_1;
      BC_Type->desc->sens[POLYMER_STRESS22_1] = 1;
    }
  else if ( ( strcmp( input, "S13_1") == 0 ) )
    {
      BC_Type->table->f_name = "S13_1";
      BC_Type->table->f_index = POLYMER_STRESS13_1;
      BC_Type->desc->equation = R_STRESS13_1;
      BC_Type->desc->sens[POLYMER_STRESS13_1] = 1;
    }
  else if ( ( strcmp( input, "S23_1") == 0 ) )
    {
      BC_Type->table->f_name = "S23_1";
      BC_Type->table->f_index = POLYMER_STRESS23_1;
      BC_Type->desc->equation = R_STRESS23_1;
      BC_Type->desc->sens[POLYMER_STRESS23_1] = 1;
    }
  else if ( ( strcmp( input, "S33_1") == 0 ) )
    {
      BC_Type->table->f_name = "S33_1";
      BC_Type->table->f_index = POLYMER_STRESS33_1;
      BC_Type->desc->equation = R_STRESS33_1;
      BC_Type->desc->sens[POLYMER_STRESS33_1] = 1;
    }
  else if ( ( strcmp( input, "S11_2") == 0 ) )
    {
      BC_Type->table->f_name = "S11_2";
      BC_Type->table->f_index = POLYMER_STRESS11_2;
      BC_Type->desc->equation = R_STRESS11_2;
      BC_Type->desc->sens[POLYMER_STRESS11_2] = 1;
    }
  else if ( ( strcmp( input, "S12_2") == 0 ) )
    {
      BC_Type->table->f_name = "S12_2";
      BC_Type->table->f_index = POLYMER_STRESS12_2;
      BC_Type->desc->equation = R_STRESS12_2;
      BC_Type->desc->sens[POLYMER_STRESS12_2] = 1;
    }
  else if ( ( strcmp( input, "S22_2") == 0 ) )
    {
      BC_Type->table->f_name = "S22_2";
      BC_Type->table->f_index = POLYMER_STRESS22_2;
      BC_Type->desc->equation = R_STRESS22_2;
      BC_Type->desc->sens[POLYMER_STRESS22_2] = 1;
    }
  else if ( ( strcmp( input, "S13_2") == 0 ) )
    {
      BC_Type->table->f_name = "S13_2";
      BC_Type->table->f_index = POLYMER_STRESS13_2;
      BC_Type->desc->equation = R_STRESS13_2;
      BC_Type->desc->sens[POLYMER_STRESS13_2] = 1;
    }
  else if ( ( strcmp( input, "S23_2") == 0 ) )
    {
      BC_Type->table->f_name = "S23_2";
      BC_Type->table->f_index = POLYMER_STRESS23_2;
      BC_Type->desc->equation = R_STRESS23_2;
      BC_Type->desc->sens[POLYMER_STRESS23_2] = 1;
    }
  else if ( ( strcmp( input, "S33_2") == 0 ) )
    {
      BC_Type->table->f_name = "S33_2";
      BC_Type->table->f_index = POLYMER_STRESS33_2;
      BC_Type->desc->equation = R_STRESS33_2;
      BC_Type->desc->sens[POLYMER_STRESS33_2] = 1;
    }
  else if ( ( strcmp( input, "S11_3") == 0 ) )
    {
      BC_Type->table->f_name = "S11_3";
      BC_Type->table->f_index = POLYMER_STRESS11_3;
      BC_Type->desc->equation = R_STRESS11_3;
      BC_Type->desc->sens[POLYMER_STRESS11_3] = 1;
    }
  else if ( ( strcmp( input, "S12_3") == 0 ) )
    {
      BC_Type->table->f_name = "S12_3";
      BC_Type->table->f_index = POLYMER_STRESS12_3;
      BC_Type->desc->equation = R_STRESS12_3;
      BC_Type->desc->sens[POLYMER_STRESS12_3] = 1;
    }
  else if ( ( strcmp( input, "S22_3") == 0 ) )
    {
      BC_Type->table->f_name = "S22_3";
      BC_Type->table->f_index = POLYMER_STRESS22_3;
      BC_Type->desc->equation = R_STRESS22_3;
      BC_Type->desc->sens[POLYMER_STRESS22_3] = 1;
    }
  else if ( ( strcmp( input, "S13_3") == 0 ) )
    {
      BC_Type->table->f_name = "S13_3";
      BC_Type->table->f_index = POLYMER_STRESS13_3;
      BC_Type->desc->equation = R_STRESS13_3;
      BC_Type->desc->sens[POLYMER_STRESS13_3] = 1;
    }
  else if ( ( strcmp( input, "S23_3") == 0 ) )
    {
      BC_Type->table->f_name = "S23_3";
      BC_Type->table->f_index = POLYMER_STRESS23_3;
      BC_Type->desc->equation = R_STRESS23_3;
      BC_Type->desc->sens[POLYMER_STRESS23_3] = 1;
    }
  else if ( ( strcmp( input, "S33_3") == 0 ) )
    {
      BC_Type->table->f_name = "S33_3";
      BC_Type->table->f_index = POLYMER_STRESS33_3;
      BC_Type->desc->equation = R_STRESS33_3;
      BC_Type->desc->sens[POLYMER_STRESS33_3] = 1;
    }
  else if ( ( strcmp( input, "S11_4") == 0 ) )
    {
      BC_Type->table->f_name = "S11_4";
      BC_Type->table->f_index = POLYMER_STRESS11_4;
      BC_Type->desc->equation = R_STRESS11_4;
      BC_Type->desc->sens[POLYMER_STRESS11_4] = 1;
    }
  else if ( ( strcmp( input, "S12_4") == 0 ) )
    {
      BC_Type->table->f_name = "S12_4";
      BC_Type->table->f_index = POLYMER_STRESS12_4;
      BC_Type->desc->equation = R_STRESS12_4;
      BC_Type->desc->sens[POLYMER_STRESS12_4] = 1;
    }
  else if ( ( strcmp( input, "S22_4") == 0 ) )
    {
      BC_Type->table->f_name = "S22_4";
      BC_Type->table->f_index = POLYMER_STRESS22_4;
      BC_Type->desc->equation = R_STRESS22_4;
      BC_Type->desc->sens[POLYMER_STRESS22_4] = 1;
    }
  else if ( ( strcmp( input, "S13_4") == 0 ) )
    {
      BC_Type->table->f_name = "S13_4";
      BC_Type->table->f_index = POLYMER_STRESS13_4;
      BC_Type->desc->equation = R_STRESS13_4;
      BC_Type->desc->sens[POLYMER_STRESS13_4] = 1;
    }
  else if ( ( strcmp( input, "S23_4") == 0 ) )
    {
      BC_Type->table->f_name = "S23_4";
      BC_Type->table->f_index = POLYMER_STRESS23_4;
      BC_Type->desc->equation = R_STRESS23_4;
      BC_Type->desc->sens[POLYMER_STRESS23_4] = 1;
    }
  else if ( ( strcmp( input, "S33_4") == 0 ) )
    {
      BC_Type->table->f_name = "S33_4";
      BC_Type->table->f_index = POLYMER_STRESS33_4;
      BC_Type->desc->equation = R_STRESS33_4;
      BC_Type->desc->sens[POLYMER_STRESS33_4] = 1;
    }
  else if ( ( strcmp( input, "S11_5") == 0 ) )
    {
      BC_Type->table->f_name = "S11_5";
      BC_Type->table->f_index = POLYMER_STRESS11_5;
      BC_Type->desc->equation = R_STRESS11_5;
      BC_Type->desc->sens[POLYMER_STRESS11_5] = 1;
    }
  else if ( ( strcmp( input, "S12_5") == 0 ) )
    {
      BC_Type->table->f_name = "S12_5";
      BC_Type->table->f_index = POLYMER_STRESS12_5;
      BC_Type->desc->equation = R_STRESS12_5;
      BC_Type->desc->sens[POLYMER_STRESS12_5] = 1;
    }
  else if ( ( strcmp( input, "S22_5") == 0 ) )
    {
      BC_Type->table->f_name = "S22_5";
      BC_Type->table->f_index = POLYMER_STRESS22_5;
      BC_Type->desc->equation = R_STRESS22_5;
      BC_Type->desc->sens[POLYMER_STRESS22_5] = 1;
    }
  else if ( ( strcmp( input, "S13_5") == 0 ) )
    {
      BC_Type->table->f_name = "S13_5";
      BC_Type->table->f_index = POLYMER_STRESS13_5;
      BC_Type->desc->equation = R_STRESS13_5;
      BC_Type->desc->sens[POLYMER_STRESS13_5] = 1;
    }
  else if ( ( strcmp( input, "S23_5") == 0 ) )
    {
      BC_Type->table->f_name = "S23_5";
      BC_Type->table->f_index = POLYMER_STRESS23_5;
      BC_Type->desc->equation = R_STRESS23_5;
      BC_Type->desc->sens[POLYMER_STRESS23_5] = 1;
    }
  else if ( ( strcmp( input, "S33_5") == 0 ) )
    {
      BC_Type->table->f_name = "S33_5";
      BC_Type->table->f_index = POLYMER_STRESS33_5;
      BC_Type->desc->equation = R_STRESS33_5;
      BC_Type->desc->sens[POLYMER_STRESS33_5] = 1;
    }
  else if ( ( strcmp( input, "S11_6") == 0 ) )
    {
      BC_Type->table->f_name = "S11_6";
      BC_Type->table->f_index = POLYMER_STRESS11_6;
      BC_Type->desc->equation = R_STRESS11_6;
      BC_Type->desc->sens[POLYMER_STRESS11_6] = 1;
    }
  else if ( ( strcmp( input, "S12_6") == 0 ) )
    {
      BC_Type->table->f_name = "S12_6";
      BC_Type->table->f_index = POLYMER_STRESS12_6;
      BC_Type->desc->equation = R_STRESS12_6;
      BC_Type->desc->sens[POLYMER_STRESS12_6] = 1;
    }
  else if ( ( strcmp( input, "S22_6") == 0 ) )
    {
      BC_Type->table->f_name = "S22_6";
      BC_Type->table->f_index = POLYMER_STRESS22_6;
      BC_Type->desc->equation = R_STRESS22_6;
      BC_Type->desc->sens[POLYMER_STRESS22_6] = 1;
    }
  else if ( ( strcmp( input, "S13_6") == 0 ) )
    {
      BC_Type->table->f_name = "S13_6";
      BC_Type->table->f_index = POLYMER_STRESS13_6;
      BC_Type->desc->equation = R_STRESS13_6;
      BC_Type->desc->sens[POLYMER_STRESS13_6] = 1;
    }
  else if ( ( strcmp( input, "S23_6") == 0 ) )
    {
      BC_Type->table->f_name = "S23_6";
      BC_Type->table->f_index = POLYMER_STRESS23_6;
      BC_Type->desc->equation = R_STRESS23_6;
      BC_Type->desc->sens[POLYMER_STRESS23_6] = 1;
    }
  else if ( ( strcmp( input, "S33_6") == 0 ) )
    {
      BC_Type->table->f_name = "S33_6";
      BC_Type->table->f_index = POLYMER_STRESS33_6;
      BC_Type->desc->equation = R_STRESS33_6;
      BC_Type->desc->sens[POLYMER_STRESS33_6] = 1;
    }
  else if ( ( strcmp( input, "S11_7") == 0 ) )
    {
      BC_Type->table->f_name = "S11_7";
      BC_Type->table->f_index = POLYMER_STRESS11_7;
      BC_Type->desc->equation = R_STRESS11_7;
      BC_Type->desc->sens[POLYMER_STRESS11_7] = 1;
    }
  else if ( ( strcmp( input, "S12_7") == 0 ) )
    {
      BC_Type->table->f_name = "S12_7";
      BC_Type->table->f_index = POLYMER_STRESS12_7;
      BC_Type->desc->equation = R_STRESS12_7;
      BC_Type->desc->sens[POLYMER_STRESS12_7] = 1;
    }
  else if ( ( strcmp( input, "S22_7") == 0 ) )
    {
      BC_Type->table->f_name = "S22_7";
      BC_Type->table->f_index = POLYMER_STRESS22_7;
      BC_Type->desc->equation = R_STRESS22_7;
      BC_Type->desc->sens[POLYMER_STRESS22_7] = 1;
    }
  else if ( ( strcmp( input, "S13_7") == 0 ) )
    {
      BC_Type->table->f_name = "S13_7";
      BC_Type->table->f_index = POLYMER_STRESS13_7;
      BC_Type->desc->equation = R_STRESS13_7;
      BC_Type->desc->sens[POLYMER_STRESS13_7] = 1;
    }
  else if ( ( strcmp( input, "S23_7") == 0 ) )
    {
      BC_Type->table->f_name = "S23_7";
      BC_Type->table->f_index = POLYMER_STRESS23_7;
      BC_Type->desc->equation = R_STRESS23_7;
      BC_Type->desc->sens[POLYMER_STRESS23_7] = 1;
    }
  else if ( ( strcmp( input, "S33_7") == 0 ) )
    {
      BC_Type->table->f_name = "S33_7";
      BC_Type->table->f_index = POLYMER_STRESS33_7;
      BC_Type->desc->equation = R_STRESS33_7;
      BC_Type->desc->sens[POLYMER_STRESS33_7] = 1;
    }
  else
    {
      sprintf(err_msg, "%s:\tInvalid choice for table ordinate.",yo);
      EH(-1,err_msg);
    }

  SPF(endofstring(es)," %s", input); 
  
     /* Read scaling factor */
 
  if ( fscanf(ifp, "%lf", &BC_Type->table->yscale) != 1)
    {
      BC_Type->table->yscale = 1.0;
    }

  SPF(endofstring(es)," %.4g", BC_Type->table->yscale); 

  /* read interpolation order */

  if ( fscanf(ifp, "%80s", input ) != 1 )
    {
      	sprintf (err_msg, "%s:\tError reading interpolation order for table: %s\n", yo, BC_Type->table->f_name );
	EH(-1,err_msg);	
    }
  
  strip( input );
  stringup( input );

  if ( strcmp( input, "LINEAR") == 0 )
    {
      BC_Type->table->interp_method = LINEAR;
    }
  else if ( strcmp( input, "QUADRATIC") == 0 )
    {
      BC_Type->table->interp_method = QUADRATIC;
    }
  else if ( strcmp( input, "QUAD_GP") == 0 )
    {
      BC_Type->table->interp_method = QUAD_GP;
    }
  else if ( strcmp( input, "BIQUADRATIC") == 0 )
    {
      BC_Type->table->interp_method = BIQUADRATIC;
      BC_Type->table->columns = 3;
         if( BC_Type->BC_Name == TABLE_WICV_BC ) {BC_Type->table->columns = 5;}
    }
  else if ( strcmp( input, "BILINEAR") == 0 )
    {
      BC_Type->table->interp_method = BILINEAR;
      BC_Type->table->columns = 3;
         if( BC_Type->BC_Name == TABLE_WICV_BC ) {BC_Type->table->columns = 5;}
    }
  else
    {
      sprintf(err_msg, "\nUnknown table interpolation order for table: %s \n",BC_Type->table->f_name);
      EH(-1,err_msg);
    }

  SPF(endofstring(es)," %s", input); 
  
  if( look_for_next_string( ifp, "FILE", input, '=') )
  {
	  SPF(endofstring(es)," %s = ", input);
	  
      if ( fscanf( ifp, "%s", input ) != 1 )
	  {
		  sprintf(err_msg,"\n%s,\tError reading TABLE data filename for table: %s\n",yo,BC_Type->table->f_name);
		  EH(-1,err_msg);
	  }
	  
      strip(input);
	  
      if ( ( datafile = fopen_aprepro( input, "r") ) == NULL )
	  {
		  sprintf(err_msg, "\n%s,\tError opening TABLE data file %s\n",yo,input);
		  EH(-1,err_msg);
	  }
	  
	  SPF(endofstring(es)," %s", input);
		  
      if( look_for_next_string( ifp, "NAME", input, '=') )
	  {
		  SPF(endofstring(es)," %s = ", input);
		  
		  if ( fscanf( ifp, "%s", input ) != 1 )
		  {
			  sprintf(err_msg,"\n%s,\tError reading NAME for table: %s\n",yo,BC_Type->table->f_name);
			  EH(-1,err_msg);
		  }
		  strip(input);
		  
		  dataname = input;

		  SPF(endofstring(echo_string)," %s", input);
	  }  
  }
  
  
  ECHO(echo_string,echo_file);
  echo_string[0]='\0';

  if ( datafile != NULL )
    {
      if( dataname != NULL)
	{
	  dataname = smalloc( ( 1 + strlen(input) )*sizeof(char) ); /*OK, this is gratuitous */
	  strcpy( dataname,input);

	  /* Position file pointer to just after dataname in the datafile */

	  if( look_for_optional( datafile, dataname, input, ':') == -1)
	    {
	      sprintf(err_msg,"\n%s,\tError finding data name %s in file for table: %s\n",
		              yo, dataname, BC_Type->table->f_name);
	      free(dataname);
	      EH(-1,err_msg);
	    }

	  rd_table_data( datafile, input, BC_Type->table, "END TABLE" );

	  free(dataname);
	}
      else
	{
	  rd_table_data( datafile, input, BC_Type->table, NULL );
	}
      

      fclose( datafile );
    }
  else
    {
      rd_table_data( ifp, input, BC_Type->table, "END TABLE" );
    }

return( BC_Type->table );
  
}

int
count_datalines( FILE *ifp, 
		 char *input, 
		 const char *endlist)
{
/*
	Scan the input file and identify and count 
	lines that have at least one readable number
	as their first parameter.  Stop reading when the
	string *endlist is read.  If endlist is NULL, 
	stop reading when the EOF character is read.

	Author:			Thomas A. Baer, Org. 9111 
	Date:			7/9/98

	Parameter list:

	ifp    == pointer to file "input"
	input  == buffer array to hold characters that are read in.
	endlist == pointer to string that signals the end of the datalist 
	           has been reached
        returns:  number of occurences lines that start with at least one 
	          readable numerical value.
*/


  fpos_t file_position; 
  int datalines = 0;

#ifndef tflop
  fgetpos(ifp, &file_position); /* OK everybody, remember where we parked */
#else
  file_position = ftell( ifp ); /* OK everybody, remember where we parked */
#endif

  if ( endlist != NULL ) /* Table data is in input deck */
    {
      size_t len = strlen(endlist);
  
      if( read_string(ifp, input, '\n') == -1 )  /* read a line from the input file */
	{
	  EH(-1,"EOF found in file reading dataline for table \n");		  
	}
      
      strip(input);

      while ( strncmp( endlist, input, len ) != 0  )
	{
	  char *p,*q;

	  if( ( strtod( input, &p) != 0.0 || p != input) && 
	      ( strtod(p, &q) != 0.0 || q != p ) )       /* Failure of both these tests indicates
							 * that this line does not contain two valid doubles
							 */
	    {
	      datalines++;
	    }

	  if( read_string(ifp, input, '\n') == -1 )  /* read another line from the input file */
	    {
	      strip(input);
	      if( strncmp( input, endlist,len ) != 0 )
		EH(-1,"EOF found in input file reading dataline\n"); 
	    }

	  strip(input);

	}
    }
  else  /* reading table data from a separate file */
    {
      while ( read_string(ifp, input, '\n')  != -1 ) /* Read until you encounter EOF */
	{
	  char *p, *q;
	  
	  strip( input );

	  if( ( strtod( input, &p) != 0.0 || p != input) && 
	      ( strtod(p, &q) != 0.0 || q != p ) )       /* Failure of both these tests indicates
							 * that this line does not contain two valid doubles
							 */
	    {
	      datalines++;
	    }

	}
    }

#ifndef tflop
  fsetpos(ifp, &file_position);  /* Pick up the car */
#else
  fseek(ifp, file_position, SEEK_SET);  /* Pick up the car */
#endif

  return( datalines );

}


int
look_for_next_string ( FILE *ifp,
		       const char *string,
		       char *input,
		       const char ch_term)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  /*            Scan file ifp ignoring all white space characters
		until the first non-whitespace character is read.
		Read characters into input until the ch_term is
		encountered.  Compare input to the target string,
		string.  Return 1 if it matches 0 if it doesn't.
		Reset file pointer to original location if no match 
		is found.

		Author:             Thomas A. Baer, Org. 9111
		Date  :             July 14, 1998 
		
       Parameters:
	        ifp = FILE * to open file
		string = pointer to target string
		input  = character buffer to read file data into.
		ch_term = termination character to stop reading characters from file
		
       Returns:
                returns 1 if the next string in the input deck matches the target string
		returns 0 if the next string in the input deck does not match target string.
		
   */

  fpos_t file_position;
  int ch;
  /*  char ch;*/

#ifndef tflop
  fgetpos( ifp, &file_position );
#else
  file_position = ftell( ifp );
#endif

  while ( isspace ( ch = getc(ifp) ) );

  ungetc( ch, ifp );

  if ( read_string( ifp, input, ch_term ) == -1 )
    {
      sprintf( err_msg, "\nError detected while searching for : %s \n", string );
      EH(-1,err_msg);
    }

  strip( input );

  if ( strcmp( string, input ) == 0 )
    {
      return (1);
    }
  else
    {
#ifndef tflop
      fsetpos (ifp, &file_position );
#else
      fseek(ifp, file_position, SEEK_SET);
#endif
      return (0);
    }

}

FILE *
fopen_aprepro( const char *filename, const char *format )
{
  /*     Front end function for fopen.  If aprepro is not being
	 used, it simply calls fopen and returns.  If aprepro is
	 enabled, it assembles a temporary filename, preprocesses
	 the file specified in filename with APREPRO, and opens
	 the resulting file.  It returns a pointer to a file in 
	 either case.

	 Author: Thomas A. Baer, Org. 9111
	 Date  : July 14, 1998


     Parameters:

        filename = pointer to string containing name of file to be opened,
	format   = pointer to string identifying the format to use on opening,
	             i.e. "r", "rb", etc.  
     Returns:
                   pointer to FILE structure of the opened file.


   */
    
  int err;
  FILE *file;
  char Tmpfilename[MAX_FNL];
  static char System_Command[MAX_SYSTEM_COMMAND_LENGTH];

  if( run_aprepro == 1)
    {
      System_Command[0] = '\0';
      Tmpfilename[0] = '\0';
      (void) strcat(Tmpfilename, "tmp.");
      (void) strcat(Tmpfilename, filename);
      sprintf(System_Command, "%s %s %s", aprepro_command, filename, Tmpfilename);

#ifndef tflop
      err = system( System_Command);

      EH(err, "System call failed in fopen_aprepro.");

      if (WEXITSTATUS(err) == 127)
	{
	  EH(-1, "System call failed, aprepro not found");
	  return NULL;
	}
#else
      EH(-1, "aprepro the input file prior to running goma");
#endif

      if ( Debug_Flag > 0 )
	{
	  fprintf(stdout, "system: %s\n", System_Command);
	}

      file = fopen( Tmpfilename, format); 
    }
  else
    {
      file = fopen( filename, format );
    }

return (file);

}



struct Data_Table *
setup_gd_table_BC(
		  FILE *ifp,
		  char *input,
		  struct Boundary_Condition *BC_Type,
		  char * echo_string)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  /*           Function allocates memory for the BC table structure.
	       Parses the remainder of the GD_TABLE BC card.  Opens a separate file for reading the 
	       tabular data if indicated on the card and calls a function
	       to read the tabular data.

	       Author: Thomas A. Baer, Org. 9111
	       Date  : July 15, 1998
	       

	Parameters:
	       ifp     =  pointer to input deck FILE structure
	       input   =  character buffer array 
	       BC_Type =  pointer to GD_TABLE_BC structure that the table is being added to

        returns:
	       pointer to table structure created in the function
   */


  char *yo, *dataname=NULL;
  FILE *datafile =NULL;
  int k;
  char *es = endofstring(echo_string);
  char *echo_file = Echo_Input_File;

  yo = "setup_gd_table_BC";

  BC_Type->table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) ) ;
  BC_Type->table->columns = 2;

  BC_Type->table->t_index[0] = BC_Type->BC_Data_Int[2];
  BC_Type->table->f_index = BC_Type->BC_Data_Int[0];

  /* Find the names of the ordinate and the abscissa */

  for( k = 0; BC_Type->BC_Data_Int[0] != EQ_Name[k].Index && k < Num_EQ_Names; k++);

  BC_Type->table->f_name = EQ_Name[k].name1;

  if ( BC_Type->BC_Data_Int[2] == GD_TIME_TABLE )
    {
      strcpy( BC_Type->table->t_name[0], "time");
    }
  else
    {
      for( k = 0; BC_Type->BC_Data_Int[2] != Var_Name[k].Index && k < Num_Var_Names; k++);

      strcpy( BC_Type->table->t_name[0],Var_Name[k].name1 );
    }


  /* read interpolation method*/

  if ( fscanf(ifp, "%80s", input ) != 1 )
    {
      sr = sprintf(err_msg, "Cannot find interpolation method for %s on %sID=%d.\n",
		   BC_Type->desc->name1,
		   BC_Type->Set_Type,
		   BC_Type->BC_ID);
      EH(-1, err_msg);
    }
  
  strip( input );
  stringup( input );

  if ( strcmp( input, "LINEAR") == 0 )
    {
      BC_Type->table->interp_method = LINEAR;
    }
  else if ( strcmp( input, "QUADRATIC") == 0 )
    {
      BC_Type->table->interp_method = QUADRATIC;
    }
  else if ( strcmp( input, "QUAD_GP") == 0 )
    {
      BC_Type->table->interp_method = QUAD_GP;
    }
  else
    {
      sr = sprintf(err_msg, "Unknown interpolation method for %s on %sID=%d.\n",
		   BC_Type->desc->name1,
		   BC_Type->Set_Type,
		   BC_Type->BC_ID);
      EH(-1, err_msg);
    }
  
  SPF(es," %s", input);
  
  /*
   * Account for 2 different methods of specifying the table BC:
   *   1) Points in a 2D table specified in a file
   *   2) Using a solid model (CGM) Edge
   */
  /*if ( LINEAR == BC_Type->table->interp_method )
  {  */
      if( look_for_next_string( ifp, "FILE", input, '=') )
	  {
		  if ( fscanf( ifp, "%s", input ) != 1 )
		  {
			  sr = sprintf(err_msg, "Error  reading table file name for %s on %sID=%d.\n",
						   BC_Type->desc->name1,
						   BC_Type->Set_Type,
						   BC_Type->BC_ID);
			  EH(-1, err_msg);
		  }
		  
		  SPF(endofstring(es)," FILE = %s", input);
		  
		  strip(input);
		  
		  if ( ( datafile = fopen_aprepro( input, "r") ) == NULL )
		  {
			  sprintf(err_msg, "\n%s:\tError opening TABLE data file %s\n",yo,input);
			  EH(-1,err_msg);
		  }
		  
		  if( look_for_next_string( ifp, "NAME", input, '=') )
		  {
			  if ( fscanf( ifp, "%s", input ) != 1 )
			  {
				  sprintf(err_msg,"\n%s,\tError reading NAME for table: %s\n",yo,BC_Type->table->f_name);
				  EH(-1,err_msg);
			  }
			  strip(input);
			  
			  dataname = input;
			  
			  SPF(endofstring(es)," NAME = %s", input);
		  }
	  }
	  
	  ECHO(echo_string,echo_file);
	  
      if ( datafile != NULL )
	  {
		  if( dataname != NULL)
		  {
			  dataname = smalloc( ( 1 + strlen(input) )*sizeof(char) ); /*OK, this is gratuitous */
			  strcpy( dataname,input);
			  
			  /* Position file pointer to just after dataname in the datafile */
			  
			  if( look_for_optional( datafile, dataname, input, ':') == -1)
			  {
				  sprintf(err_msg,"\n%s,\tError finding data name %s in file for table: %s\n",
						  yo, dataname, BC_Type->table->f_name);
				  free(dataname);
				  EH(-1,err_msg);
			  }
			  
			  rd_table_data( datafile, input, BC_Type->table, "END TABLE" );
			  
			  free(dataname);
		  }
		  else
		  {
			  rd_table_data( datafile, input, BC_Type->table, NULL );
		  }
		  
		  fclose( datafile );
	  }
      else
	  {
		  rd_table_data( ifp, input, BC_Type->table, "END TABLE" );
	  }
 /* }  */

  return( BC_Type->table );
}



struct Data_Table *
setup_table_MP (FILE *imp, struct Data_Table * table, char *search_string)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  /*           Function allocates memory for the MP table structure.
	       Parses the remainder of the TABLE MP card identifying
	       the abscissa and ordinate variables and setting the
	       members of the table structure and the MP_desc structure
	       accordingly.  Opens a separate file for reading the 
	       tabular data if indicated on the card and calls a function
	       to read the tabular data.

	       Author: Thomas A. Baer, Org. 9111
	       Date  : July 15, 1998
	       revised:  raroach October 12, 1999
	       

	Parameters:
	       imp     =  pointer to input deck FILE structure
	       input   =  character buffer array 
	       BC_Type =  pointer to TABLE BC structure that the table is being added to

        returns:
	       pointer to table structure created in the function
   */
  char *yo, *dataname=NULL, line[132];
  int i,num_const,j;
  double scale;
  FILE *datafile = NULL;

  yo = "setup_table_MP";

  /* Set dependent material property */
  table->f_name = search_string;
  table->f_index = *search_string;
 
  /* read number of columns in table */
  if( fscanf(imp, "%80i", &num_const) != 1 )
     {
       sprintf( err_msg,
	"Error reading # columns for TABLE model for material property");
       EH(-1,err_msg);
     }
	
   if( (num_const > 3) && (strcmp( table->f_name, "Saturation") != 0))
     {
       sprintf( err_msg,
	"Error - Multi DOF table lookup limited to bilinear");
       EH(-1,err_msg);
     }
   if( (num_const > 4) && (strcmp( table->f_name, "Saturation") == 0))
     {
       sprintf( err_msg,
	"Error - Multi DOF table lookup limited to bilinear for Saturation");
       EH(-1,err_msg);
     }
	
   table->columns = num_const;

   /* read independent variable names */
     for(i=0;i<table->columns-1;i++)
       {
       if( fscanf(imp, "%80s", line) !=1 )
	 {
	   sprintf( err_msg,
	    "Error reading TABLE model for material property for Independent variables");
	   EH(-1,err_msg);
	 }

        strip(line);
	stringup(line);

	if( (strcmp( line, "TEMPERATURE") == 0) || (strcmp( line, "T") == 0))
	  {
	    strcpy( table->t_name[i],"TEMPERATURE");

	    table->t_index[i] = TEMPERATURE;
	    
	    for(j=0;j<i;j++)
	      {
		if(strcmp(table->t_name[j],table->t_name[i])==0)
		  {
		    sprintf (err_msg, "%s:Cannot set Temperature multiple times as Independent Variable \n", yo);		                EH(-1,err_msg);
		  }
	      }
	  }
	else if ( (strcmp( line, "MASS_FRACTION") == 0) || (strcmp( line, "SPECIES") == 0) ||
		  ( strcmp( line, "Y") == 0 ))
	  {
	    strcpy( table->t_name[i],"MASS_FRACTION");
	    table->t_index[i] = MASS_FRACTION;
	    if ( fscanf( imp, "%d", &table->species_eq) != 1)
	      {
		sprintf (err_msg, "%s:\tError reading species number on TABLE MP \n", yo);
		EH(-1,err_msg);
	      }
	  }
	else if (strcmp( line, "LINEAR_TIME") == 0) {
	  strcpy( table->t_name[i],"LINEAR_TIME");
	  table->t_index[i] = LINEAR_TIME;
	}
	else if( (strcmp( line, "CAP_PRES") == 0) )
	  {
	    if( (strcmp( table->f_name, "Saturation") != 0) ) 
	      {
		sprintf (err_msg, "%s:\tError: Variable CAP_PRES can only be used with Saturation \n", yo);
		EH(-1,err_msg);
	      }
	    strcpy( table->t_name[i],"CAP_PRES");
	    table->t_index[i] = CAP_PRES;
	    
	    for(j=0;j<i;j++)
	      {
		if(strcmp(table->t_name[j],table->t_name[i])==0)
		  {
		    sprintf (err_msg, "%s:Cannot set CAP_PRES multiple times as Independent Variable \n", yo);		                EH(-1,err_msg);
		  }
	      }
	  }
	else if (strcmp( line, "FAUX_PLASTIC") == 0)
	  {
	    strcpy( table->t_name[i],"FAUX_PLASTIC");
	    table->t_index[i] = FAUX_PLASTIC;
	    if ( fscanf( imp, "%d", &table->species_eq) != 1)
	      {
		sprintf (err_msg, "%s:\tError reading species number on TABLE MP \n", yo);
		EH(-1,err_msg);
	      }
	  }
	else if (strcmp( line, "LOWER_DISTANCE") == 0)
	  {
	    strcpy( table->t_name[i],"LOWER_DISTANCE");
	    table->t_index[i] = LOWER_DISTANCE;
	    if ( fscanf( imp, "%d", &table->species_eq) != 1)
	      {
		sprintf (err_msg, "%s:\tError reading species number on TABLE MP \n", yo);
		EH(-1,err_msg);
	      }
	  }
	else
	  {
	    sprintf (err_msg,"%s:\tInvalid choice for material property table independent variable",yo);
	    EH(-1, err_msg);
	  }
	if(i == 0) { strcpy( table->t_name[i+1],"NULL"); table->t_index[i+1]=-1;}
       }

     /* read interpolation scheme */
     if( fscanf(imp, "%80s", line) !=1 )
       {
	  sprintf( err_msg,
	    "Error reading TABLE %s for material property %s for interpolation", table->f_name, line);
	  EH(-1,err_msg);
       }

      strip(line);
      stringup(line);


      if( (strcmp( line, "LINEAR") == 0) )
	{
	  if( (table->columns == 2) || (table->columns == 3 && (strcmp( table->f_name, "Saturation") == 0)))
	    {
	      table->interp_method = LINEAR;
	    }
	  else 
	    {
	      sprintf( err_msg, " Incorrect number of columns for material property table lookup");
	      EH(-1, err_msg);
	    }
     	}
      else if( (strcmp( line, "BILINEAR") == 0) )
	{
	  if( table->columns >= 3)
	    {
	      table->interp_method = BILINEAR;
	    }
	  else 
	    {
	      sprintf( err_msg, " Incorrect number of columns for material property table lookup");
	      EH(-1, err_msg);
	    }
     	}
      else if( (strcmp( line, "QUADRATIC") == 0) )
	{
	  if( (table->columns == 2) )
	    {
	      table->interp_method = QUADRATIC;
	    }
	  else 
	    {
	      sprintf( err_msg, " Incorrect number of columns for material property table lookup");
	      EH(-1, err_msg);
	    }
     	}
      else if( (strcmp( line, "QUAD_GP") == 0) )
	{
	  if( table->columns == 2)
	    {
	      table->interp_method = QUAD_GP;
	    }
	  else 
	    {
	      sprintf( err_msg, " Incorrect number of columns for material property table lookup");
	      EH(-1, err_msg);
	    }
     	}
      else
	{
	  sprintf (err_msg,"%s:\tInvalid choice for material property table interpolation scheme",yo);
	  EH(-1, err_msg);
	}

     /* read file name (optional) */ 
    if( look_for_next_string( imp, "FILE", line, '=') )
    {
      if ( fscanf( imp, "%s", line ) != 1 )
	{
	  sprintf(err_msg,"\n%s,\tError reading TABLE data filename for material property \n",yo);
	  EH(-1,err_msg);
	}

      strip(line);

      if ( ( datafile = fopen_aprepro( line, "r") ) == NULL )
	{
	  sprintf(err_msg, "\n%s,\tError opening TABLE data file for material property \n",yo);
	  EH(-1,err_msg);
	}

      if( look_for_next_string( imp, "NAME", line, '=') )
	{
	  if ( fscanf( imp, "%s", line ) != 1 )
	    {
	      sprintf(err_msg,"\n%s,\tError reading NAME for table for material property \n",yo);
	      EH(-1,err_msg);
	    }
	  strip(line);

	  dataname = line;
	}  
      
      // Looks for optional scale factor for the y axis
      if ( look_for_next_string( imp, "YSCALE", line, '=') ) {
        if ( fscanf( imp, "%lf ", &scale ) != 1 ) {
          sprintf(err_msg,"\n%s,\tError reading YSCALE for table for material property \n",yo);
          EH(-1,err_msg);
        } else {
          printf("Table y-axis scaled by %f\n", scale);
        }
        strip(line);
      } else {
        scale = 1.0;
      }
      table->yscale = scale;
      
      // Look for optional EMin value for FAUX_PLASTICITY model of Lame MU
      if ( look_for_next_string( imp, "EMIN", line, '=') ) {
        if ( fscanf( imp, "%lf ", &scale ) != 1 ) {
          sprintf(err_msg,"\n%s,\tError reading EMIN for table for material property \n",yo);
          EH(-1,err_msg);
        } else {
          printf("Table Emin set to %e\n", scale);
        }
        strip(line);
      } else {
        scale = -1.0;
      }
      table->Emin = scale;
      
      
    }
  

  if ( datafile != NULL )
    {
      if( dataname != NULL)
	{
	  dataname = smalloc( ( 1 + strlen(line) )*sizeof(char) ); 
	  strcpy( dataname,line);

	  /* Position file pointer to just after dataname in the datafile */

	  if( look_for_optional( datafile, dataname, line, ':') == -1)
	    {
	      sprintf(err_msg,"\n%s,\tError finding data name %s in file for table for material property \n",
		              yo, dataname);
	      free(dataname);
	      EH(-1,err_msg);
	    }

	  rd_table_data( datafile, line, table, "END TABLE" );

	  free(dataname);
	}
      else
	{
	  rd_table_data( datafile, line, table, NULL );
	}
      

      fclose( datafile );
    }
  else
    {
      rd_table_data( imp, line, table, "END TABLE" );
    }

    return(table);
}


void
rd_table_data(FILE *ifp, char *input, struct Data_Table *table , char *endlist)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  /*       Scan input file for data lines, defined as any line in which the first
	   parameter encountered on the line can be converted to a valid double.
	   Count the number of such lines. Then read the data on each of this
	   lines and store it in the arrays in the table structure.  Reading will 
	   continue until either the endlist string is read or the EOF is encountered.

	   Author : Thomas A. Baer Org. 9111
	   Date   : July 14, 1998
	   Revised: rbsecor Sept. 15, 2004

      Parameters:

           ifp        = pointer to FILE structure of the input deck.
	   input      = character array that serves to buffer the input read
	   Data_Table = a pointer to a Data_Table structure.  The read data will
	                be stored there.

	   endlist    = a pointer to a string that markes the end of table data in the 
	                input deck.

   */
       
  double p,p2=0.0,p3,p4=0.0;
  int  i,j,k, Num_Pnts,ibegin,iend;
  int table_dim=0;
  char echo_string[MAX_CHAR_IN_INPUT]="\0";	
  char *echo_file = Echo_Input_File;

  

  /*
   * Count the number of PNTs in the TABLE
   */

  table->tablelength = Num_Pnts = count_datalines(ifp,input, endlist);
      
  if( table->tablelength == 0 )
    EH( -1, "Error reading tabular data . Can't find any points ");   
  
  if( table->tablelength == 1 )
    EH( -1, "Error reading tabular data . Need more than 1 point ");
  
  if( (table->interp_method == QUADRATIC || 
	table->interp_method == BIQUADRATIC ||
	table->interp_method == TRIQUADRATIC)
                 && Num_Pnts%2 != 1)
    EH( -1, "Need odd number of points for (bi,tri)quadratic interpolation ");
  if(table->interp_method == QUAD_GP && Num_Pnts%3 != 0)
    EH( -1, "Need 3N points for quadratic gauss point interpolation ");
  
  /* determine table dimension  */
 
  switch (table->interp_method)
	{
	case LINEAR:
	case QUADRATIC:
 	case QUAD_GP:
 		table_dim = 1;
 		break;
 	case BILINEAR:
 	case BIQUADRATIC:
 		table_dim = 2;
 		break;
 	case TRILINEAR:
 	case TRIQUADRATIC:
 		table_dim = 3;
 		break;
 	}
 
  /* 
   * Assign memory
   */

   switch (table_dim)
 	{
 	case 3:
   		table->t3 = (double *) smalloc( sizeof(double)*Num_Pnts );
   		/* fall through */
 	case 2:
 	case 1:
   		table->t2 = (double *) smalloc( sizeof(double)*Num_Pnts );
   		table->t = (double *) smalloc( sizeof(double)*Num_Pnts );
 		break;
 	}
 
   table->f = (double *) smalloc( sizeof(double)*Num_Pnts*
 					(table->columns-table_dim) );

  /* 
   * Now read all yer points
   */

  if ( endlist != NULL )  /* read data from input deck */
  {
	  
	  size_t len = strlen(endlist);
	  
	  /* read a line from the input file */
      if( read_string(ifp, input, '\n') == -1 )  
	  {
		  if( strncmp( input, endlist, strlen(endlist) ) != 0 )
			  EH(-1,"EOF found in input file reading dataline\n");		  
	  }
	  
	  strip(input);
	  
	  k = 0;
	  
	  while ( strncmp( endlist, input, len ) != 0 )  /* read until it finds the endlist pattern */
	  {
		  char *p;
		  
		  if( strtod( input, &p) != 0.0 || p != input ) /* Failure of both these tests indicates
			  * that strtod can't find a valid string 
			  * to convert to a double .i.e. this isn't 
			  * a dataline.
			  */
		  {
			  if( scan_table_columns( k, input, table, table_dim, err_msg, echo_string ) == 1 )
			  {
				  EH(-1, err_msg );
			  }
			  else
				  k++;
			  
			  ECHO(echo_string, echo_file);
			  
			/*  if( table->columns == 2)
			  {
				  if( sscanf( input, "%lf %lf", table->t+k, table->f+k ) != 2)
				  {
					  sprintf(err_msg,"\nCannot read two doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g)", table->t+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==3 && table_dim == 2 )
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->t2+k, table->f+k ) != 3)
				  {
					  sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  
				  SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==3 && 
					  (table_dim != 2 && strcmp(table->f_name, "Saturation") != 0) )
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+k, table->f+Num_Pnts+k ) != 3)
				  {
					  sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t+(k-1), table->f+(k-1), table->f+Num_Pnts+(k-1) ); 
			  }
			  else if(table->columns ==3 && table_dim == 1
					  && strcmp(table->f_name, "Saturation") == 0)
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+Num_Pnts+k, table->f+k ) != 3)
				  {
					  sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t+(k-1), table->f+Num_Pnts+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==4 && table_dim == 2 
					  && strcmp(table->f_name, "Saturation") == 0)
			  {
				  if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->f+Num_Pnts+k, table->f+k  ) != 4)
				  {
					  sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), 
					  table->f+Num_Pnts+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==4 && table_dim == 3 )
			  {
				  if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->t3+k, table->f+k  ) != 4)
				  {
					  sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), table->t3+(k-1), 
					  table->f+(k-1) ); 
			  }
			  else if(table->columns ==5 && table_dim == 2 )
			  {
				  if( sscanf( input, "%lf %lf %lf %lf %lf",
							  table->t+k, table->t2+k,
							  table->f+k, table->f+Num_Pnts+k, table->f+2*Num_Pnts+k) != 5)
				  {
					  sprintf(err_msg,"\nCannot read five doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), 
					  table->f+(k-1), table->f+Num_Pnts+(k-1), table->f+2*Num_Pnts+(k-1)); 
			  }
			  else
			  {
				  EH(-1,"invalid table column, dimension combination.");
			  }
			  */
		  }
		  /* read another line from the input file */
		  if( read_string(ifp, input, '\n') == -1 )  
		  {
			  /* The point of this next bit is that it is possible
			  * that when reading multiple functions from a single file by
			  * using the NAME = function, that the last END TABLE may not be followed
			  * by a \n character.  This will handle that possibility without bombing.
			  * Yes, I know it sounds like an esoteric bug, but the first time I 
			  * tried testing the NAME = stuff I did indeed leave off an endline 
			  * at the bottom of the table file. TAB
			  */
			  
			  if( strncmp( input, endlist, strlen(endlist) ) != 0 )
				  EH(-1,"EOF found in input file reading dataline\n");  
		  }
		  strip(input);
	  }
	  ECHO(endlist, echo_file);
  }
  else   /* read data from separate data file */
  {
      k = 0;
	  
	  /* read until you encounter the EOF */
      while ( read_string(ifp, input, '\n') != -1 ) 
	  {
		  char *p;
		  
		  strip(input);
		  
		  /* Failure of both these tests indicates
			  * that strtod can't find a valid string 
			  * to convert to a double .i.e. this isn't 
			  * a dataline.
			  */
		  
		  if( strtod( input, &p) != 0.0 || p != input ) 
			  
		  {
			  
			  if( scan_table_columns( k, input, table, table_dim, err_msg, echo_string ) == 1 )
			  {
				  EH(-1, err_msg );
			  }
			  else
				  k++;
			  
			  ECHO(echo_string, echo_file);
			  
			  /*
			  if( table->columns == 2)
			  {
				  if( sscanf( input, "%lf %lf", table->t+k, table->f+k ) != 2)
				  {
					  sprintf(err_msg,"\nCannot read two doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
				  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g)", table->t+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==3 && table_dim == 2 )
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->t2+k, table->f+k ) != 3)
				  {
                      sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
                      sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
                      WH(-1,err_msg);
				  }
                  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==3 && table_dim != 2 
					  && strcmp(table->f_name, "Saturation") != 0 )
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+k, table->f+Num_Pnts+k ) != 3)
				  {
                      sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
                      sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
                      WH(-1,err_msg);
				  }
                  else
					  k++;
			  }
			  else if(table->columns ==3 && table_dim == 1 
					  && strcmp(table->f_name, "Saturation") == 0)
			  {
				  if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+Num_Pnts+k, table->f+k ) != 3)
				  {
                      sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
                      sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
                      WH(-1,err_msg);
				  }
                  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t+(k-1), table->f+(k-1), table->f+Num_Pnts+(k-1) ); 
			  }
			  else if(table->columns ==4 && table_dim == 2 
					  && strcmp(table->f_name, "Saturation") == 0)
			  {
				  if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->f+Num_Pnts+k, table->f+k  ) != 4)
				  {
                      sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
                      sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
                      WH(-1,err_msg);
				  }
                  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), 
					  table->f+Num_Pnts+(k-1), table->f+(k-1) ); 
			  }
			  else if(table->columns ==4 && table_dim == 3 )
			  {
                  if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->t3+k, table->f+k  ) != 4)
				  {
					  sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
					  sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
					  WH(-1,err_msg);
				  }
                  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), table->t3+(k-1), 
					  table->f+(k-1) ); 
			  }
			  else if(table->columns ==5 && table_dim ==2 )
			  {
				  if( sscanf( input, "%lf %lf %lf %lf %lf",
							  table->t+k, table->t2+k,
							  table->f+k, table->f+Num_Pnts+k,
							  table->f+2*Num_Pnts+k) != 5)
				  {
                      sprintf(err_msg,"\nCannot read five doubles for Table: %s \n", table->f_name);
                      sprintf(err_msg + strlen(err_msg), "\nRead %d data values prior to error.\n", k);
                      WH(-1,err_msg);
				  }
                  else
					  k++;
				  SPF(echo_string,"\t(%.4g %.4g %.4g %.4g %.4g)", table->t+(k-1), table->t2+(k-1), 
					  table->f+(k-1), table->f+Num_Pnts+(k-1), table->f+2*Num_Pnts+(k-1)); 
			  }
			  else
			  {
				  EH(-1,"invalid table column, dimension combination.");
			  }
			   */
		  }
	  }
  }

  SPF(echo_string,"\t(Found %d data points in Table : %s)",k,table->f_name); ECHO(echo_string,echo_file);


  /* 
   * Now sort the points from lowest to highest abscissa (Because people forget)
   */

if((table_dim == 1 || table->interp_method == BILINEAR) && strcmp(table->f_name, "Saturation")!= 0 )  
   /* don't sort 2d bc tables or saturation mp tables*/
 
  /* Personally, I don't think we should be sorting bilinear tables either
	as this limits the type of tables we can read.  But for now, I'll 
	leave it in. - RBS		*/
    {
      for ( i = 1; i<Num_Pnts; i++)
	{
	  p=table->t[i];
 	  if(table_dim == 2) p2=table->t2[i];
	  p3=table->f[i];
          if((table->columns-table_dim) > 1)
                        {p4=table->f[i+Num_Pnts];}
	  j=i-1;
	  while (j>=0 && (p<table->t[j]) )
	    {
	      table->t[j+1]=table->t[j]; table->t[j]=p;
 	      if(table_dim == 2)
 		{
 		 table->t2[j+1]=table->t2[j]; table->t2[j]=p2;
 		}
	      table->f[j+1]=table->f[j]; table->f[j]=p3;
              if((table->columns-table_dim) > 1)
                        {
 			 table->f[j+1+Num_Pnts]=table->f[j+Num_Pnts];
                         table->f[j+Num_Pnts]=p4;
 			}
 	      if(( table->t[i] == table->t[j]) && table_dim == 1) 
		{
		  sprintf(err_msg, "\nMultivalued function detected in table :%s \n",table->f_name);
		  EH(-1,err_msg);
		}
	      j--;
	    }

  	}	      

      /* If 3 Columns then sort 2nd Column */
      if(table->columns == 3 && table->interp_method == BILINEAR)
	{
	  ibegin=1;
	  for(k=1;k<Num_Pnts;k++)
	    {
	      if((table->t[k-1] != table->t[k]) || (k==Num_Pnts-1) )
		{
		  if(k==Num_Pnts-1)
		    {
		      iend=Num_Pnts;
		    }
		    else
		    {
		      iend=k;
		    }
      		  for(i=ibegin;i<iend; i++)
		    {
		      p2=table->t2[i];
		      p3=table->f[i];
		      j=i-1;
		      while (j>=ibegin-1 && (p2<table->t2[j]) )
			{
			  table->t2[j+1]=table->t2[j]; table->t2[j]=p2;
			  table->f[j+1]=table->f[j]; table->f[j]=p3;
			  j--;
			}
		    }
		  ibegin=iend+1;
		}
	    }
	}
    }   /* if !BIQUAD */

	/*   How about if we do some error checking here
		i.e. check for duplicate abscissa values in 1D tables
	*/

 	if(table_dim == 1)
	{
      	     for( i=1; i < Num_Pnts ; i++)
        	{
                   if(table->t[i] == table->t[i-1])
                	{
                 fprintf(stderr, "\nMultivalued function detected in table \n");
                 EH(-1,"Fatal Error");
	                }
       		}
	}

 	/*  determine number of grid points for 2d & 3d tables */
 	if(table_dim > 1 && table->interp_method != BILINEAR)
	{
	double a, b, c, cosineC;
		for ( i=2 ; i < Num_Pnts ; i++ )
		   {
			a = SQUARE(table->t[0]-table->t[i-1]) +
				SQUARE(table->t2[0]-table->t2[i-1]);
			b = SQUARE(table->t[0]-table->t[i]) +
				SQUARE(table->t2[0]-table->t2[i]);
			c = SQUARE(table->t[i]-table->t[i-1]) +
				SQUARE(table->t2[i]-table->t2[i-1]);
 			if( table_dim == 3)
 				{
 				a += SQUARE(table->t3[0]-table->t3[i-1]);
 				b += SQUARE(table->t3[0]-table->t3[i]);
 				c += SQUARE(table->t3[i]-table->t3[i-1]);
 				}
			cosineC = (a+b-c)/(2.*sqrt(a*b));
			if(fabs(cosineC) < 0.7)
				{
				table->ngrid = i;
				break;
				}
		    }
 		if(Num_Pnts%table->ngrid != 0)
 		    {
 		     fprintf(stderr,"\n2D table is not a rectangular grid \n");
 		     fprintf(stderr,"%d x %d != %d points\n",
 			Num_Pnts/table->ngrid ,table->ngrid,Num_Pnts);
                     EH(-1,"Fatal Table Error");
 		     }
 
 	/*  for 3d tables, find grid points traversing first 2 directions */
 
 
 	     if( table_dim == 3)
		{
 		ibegin = table->ngrid - 1;
 		iend = 3*table->ngrid - 1;
 		for ( i=iend ; i < Num_Pnts ; i=i+table->ngrid )
 		   {
 			j = i - table->ngrid;
 			a = SQUARE(table->t[ibegin]-table->t[j]) +
 			    SQUARE(table->t2[ibegin]-table->t2[j]) +
 			    SQUARE(table->t3[ibegin]-table->t3[j]);
 			b = SQUARE(table->t[ibegin]-table->t[i]) +
 			    SQUARE(table->t2[ibegin]-table->t2[i]) +
 			    SQUARE(table->t3[ibegin]-table->t3[i]);
 			c = SQUARE(table->t[i]-table->t[j]) +
 			    SQUARE(table->t2[i]-table->t2[j]) +
 			    SQUARE(table->t3[i]-table->t3[j]);
 			cosineC = (a+b-c)/(2.*sqrt(a*b));
 			if(fabs(cosineC) < 0.7)
 				{
 				table->ngrid2 = i-table->ngrid+1;
 				break;
 				}
 		    }
 		if(Num_Pnts%table->ngrid2 != 0)
 		    {
 		     fprintf(stderr,"\n3D table is not a rectangular grid \n");
 		     fprintf(stderr,"%d x %d != %d points\n",table->ngrid ,
 				table->ngrid2 ,Num_Pnts);
                     EH(-1,"Fatal Table Error");
 		    }
		}
	}
}    


int
scan_table_columns( int k,
					char *input,
					struct Data_Table *table,
					int table_dim,
					char *err_msg,
					char *echo_string )
{
	int Num_Pnts = table->tablelength;
	int err_stat = 0;
	
	err_msg[0] = '\0';
	
	if( table->columns == 2)
	{
		if( sscanf( input, "%lf %lf", table->t+k, table->f+k ) != 2)
		{
			sprintf(err_msg,"\nCannot read two doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else	
			SPF(echo_string,"\t(%.4g %.4g)", table->t[k], table->f[k] ); 
	}
	else if(table->columns ==3 && table_dim == 2 )
	{
		if( sscanf( input, "%lf %lf %lf", table->t+k, table->t2+k, table->f+k ) != 3)
		{
			sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t[k], table->t2[k], table->f[k] ); 
	}
	else if(table->columns ==3 && 
			(table_dim != 2 && strcmp(table->f_name, "Saturation") != 0) )
	{
		if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+k, table->f+Num_Pnts+k ) != 3)
		{
			sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t[k], table->f[k], table->f[Num_Pnts+k] ); 
	}
	else if(table->columns ==3 && table_dim == 1
			&& strcmp(table->f_name, "Saturation") == 0)
	{
		if( sscanf( input, "%lf %lf %lf", table->t+k, table->f+Num_Pnts+k, table->f+k ) != 3)
		{
			sprintf(err_msg,"\nCannot read three doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g)", table->t[k], table->f[Num_Pnts+k], table->f[k] ); 
	}
	else if(table->columns ==4 && table_dim == 2 
			&& strcmp(table->f_name, "Saturation") == 0)
	{
		if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->f+Num_Pnts+k, table->f+k  ) != 4)
		{
			sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", *(table->t+(k-1)), *(table->t2+(k-1)), 
			    *(table->f+Num_Pnts+(k-1)), *(table->f+(k-1)) ); 
	}
	else if(table->columns ==4 && table_dim == 3 )
	{
		if( sscanf( input, "%lf %lf %lf %lf", table->t+k, table->t2+k, table->t3+k, table->f+k  ) != 4)
		{
			sprintf(err_msg,"\nCannot read four doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g %.4g)", table->t[k], table->t2[k], table->t3[k], 
			table->f[k] ); 
	}
	else if(table->columns ==5 && table_dim == 2 )
	{
		if( sscanf( input, "%lf %lf %lf %lf %lf",
					table->t+k, table->t2+k,
					table->f+k, table->f+Num_Pnts+k, table->f+2*Num_Pnts+k) != 5)
		{
			sprintf(err_msg,"\nCannot read five doubles for Table: %s \n", table->f_name);
			sprintf(endofstring(err_msg), "\nRead %d data values prior to error.\n", k);
			err_stat = 1;
		}
		else
			SPF(echo_string,"\t(%.4g %.4g %.4g %.4g %.4g)", *(table->t+(k-1)), *(table->t2+(k-1)), 
			    table->f[k], table->f[Num_Pnts+k], table->f[2*Num_Pnts+k]); 
	}
	else
	{
		EH(-1,"invalid table column, dimension combination.");
	}
	return(err_stat);
}	





/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 * Here's a small routine that might someday
 * be employed to perform as its name implies.
 * At the moment, it only prints out some information
 * on #define settings
 */

void 
echo_compiler_settings()
{
  FILE * echo_file;
  
  echo_file = ECHO("FILE_POINTER", Echo_Input_File );
  if( echo_file == NULL ) echo_file = stdout;

  ECHO("\n---Compiler defines---\n", Echo_Input_File );
	
  fprintf(echo_file, "%-30s= %d\n", "MDE", MDE);
  fprintf(echo_file, "%-30s= %d\n", "MAX_PROB_VAR", MAX_PROB_VAR);
  fprintf(echo_file, "%-30s= %d\n", "MAX_CONC", MAX_CONC);
  fprintf(echo_file, "%-30s= %d\n", "MAX_VARIABLE_TYPES", MAX_VARIABLE_TYPES);
#ifdef HAVE_CONFIG_H
  fprintf(echo_file, "%-30s= %s\n", "HAVE_CONFIG_H", "yes");
#else 
  fprintf(echo_file, "%-30s= %s\n", "HAVE_CONFIG_H", "no");
#endif 

#ifdef HAVE_MPI_H
       fprintf(echo_file, "%-30s= %d\n", "PARALLEL             Num_Proc", Num_Proc);
#else
       fprintf(echo_file, "%-30s= %d\n", "SERIAL               Num_Proc", Num_Proc);
#endif

#ifdef HAVE_MPI_H
       fprintf(echo_file, "%-30s= %s\n", "HAVE_MPI_H", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_MPI_H", "no");
#endif

#ifdef USE_CHEMKIN
       fprintf(echo_file, "%-30s= %s\n", "USE_CHEMKIN", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "USE_CHEMKIN", "no");
#endif

#ifdef HAVE_FRONT
       fprintf(echo_file, "%-30s= %s\n", "HAVE_FRONT", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_FRONT", "no");
#endif

#ifdef HAVE_UMFPACK
       fprintf(echo_file, "%-30s= %s\n", "HAVE_UMFPACK", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_UMFPACK", "no");
#endif

#ifdef HAVE_SPARSE
       fprintf(echo_file, "%-30s= %s\n", "HAVE_SPARSE", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_SPARSE", "no");
#endif

#ifdef GOMA_HAVE_BLAS
       fprintf(echo_file, "%-30s= %s\n", "HAVE_BLAS", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_BLAS", "no");
#endif

#ifdef GOMA_HAVE_LAPACK
       fprintf(echo_file, "%-30s= %s\n", "HAVE_LAPACK", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_LAPACK", "no");
#endif

#ifdef HAVE_Y12M
       fprintf(echo_file, "%-30s= %s\n", "HAVE_Y12M", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_Y12M", "no");
#endif

#ifdef HAVE_ARPACK
       fprintf(echo_file, "%-30s= %s\n", "HAVE_ARPACK", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_ARPACK", "no");
#endif

#ifdef HAVE_PARPACK
       fprintf(echo_file, "%-30s= %s\n", "HAVE_PARPACK", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_PARPACK", "no");
#endif

#ifdef HAVE_AZTEC
       fprintf(echo_file, "%-30s= %s\n", "HAVE_AZTEC", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "HAVE_AZTEC", "no");
#endif

#ifdef TRILINOS
       fprintf(echo_file, "%-30s= %s\n", "TRILINOS", "yes");
#else
       fprintf(echo_file, "%-30s= %s\n", "TRILINOS", "no");
#endif

#ifdef COUPLED_FILL
  fprintf(echo_file, "%-30s= %s\n", "COUPLED_FILL", "yes");
#else 
  fprintf(echo_file, "%-30s= %s\n", "COUPLED_FILL", "no");
#endif

#ifdef DEBUG
  fprintf(echo_file, "%-30s= %s\n", "DEBUG", "yes");
#else 
  fprintf(echo_file, "%-30s= %s\n", "DEBUG", "no");
#endif

  fprintf(echo_file, "%-30s= %s\n", "Pressure Stabilization (PSPG)", (PSPG > 0 ? "yes":"no") );
  fprintf(echo_file, "%-30s= %s\n", "Pressure Stabilization (PSPP)", (PSPP > 0 ? "yes":"no") );
  fprintf(echo_file, "%-30s= %f\n", "Stabilization Scaling ", PS_scaling );
  fprintf(echo_file, "%-30s= %s\n", "Linear Stability",
	  (Linear_Stability == LSA_NORMAL ? "yes":
	   (Linear_Stability == LSA_SAVE ? "file" :
	    (Linear_Stability == LSA_3D_OF_2D ? "3D" :
	     (Linear_Stability == LSA_3D_OF_2D_SAVE ? "3Dfile" :
	      "no")))));
  
  return;
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

static double
parse_press_datum_input(const char *input)
    
    /**********************************************************************
     *
     *  parse_press_datum_input:
     *
     *     Parse the pressure datum input field. One or two arguments are
     *  allowed. The first argument is the value expressed in atmospheres.
     *  The second optional value is a string representation of the units
     *  The default is atmospheres. Understood is "torr" and "cgs".
     *  Additions are welcome. On output the value is changed to cgs units.
     *  For no units conversion, specify cgs.
     **********************************************************************/
{
  double value;
  int numTok, units = 0;
  static const char yo[] = "parse_press_datum_input";

  TOKEN_STRUCT tok;
  numTok = fillTokStruct(&tok, input);
  if (numTok < 1) goto L_ERROR;
  if (numTok > 2) goto L_ERROR;
  if (! interpret_double(tok.tok_ptr[0], &value)) goto L_ERROR;
  if (numTok == 2) {
    if (!strncasecmp(tok.tok_ptr[1], "atm", 3)) {
      units = 0;
    } else if (!strncasecmp(tok.tok_ptr[1], "torr", 4)) {
      units = 1;
    } else if (!strncasecmp(tok.tok_ptr[1], "cgs", 3)) {
      units = 2;
    } else {
      goto L_ERROR;
    }
  }

  if (units == 0) {
    value *= 1.01325E6;
  } else if (units == 1) {
    value *= 1.01325E6 / 760;
  } else if (units == 2) {
    value *= 1.0;
  }
  return value;
  /*
   * Error 
   */
 L_ERROR:;
  (void) fprintf(stderr, 
		 "%s error: reading pressure datum, %s\n", yo, input);

  EH(-1, "Pressure Datum");
  return value;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_MAT_line(FILE *ifp, int mn, char *input)

    /************************************************************************
     *
     *  read_MAT_line()
     *
     *    This function reads and processes the MAT keyword line.
     *  Format of the line that it expects.
     *
     *          MAT = Material_Name Elem_Blk1 Elem_Blk2 ... ElemBlk_n
     *
     *  The number of element blocks is variable, but must be equal to or
     *  or great than 1.
     *
     *       Material_Name = Name of the material. There must be a file
     *                       named "Material_Name.mat" in the current
     *                       directory to read the results
     *       Elem_Blk[1-n] = ID's of the element blocks encompassing this
     *                       material. The ID's must match up with ID's from
     *                       the exodus file.
     *                       The total number of integers on the line
     *                       determines the total number of element blocks
     *                       encompassed by this material
     *
     *  Errors processing this line cause the program to error exit.
     *
     * Input
     * -------
     *   ifp -> file pointer. This routine will advance the file pointer to
     *          just past the "MAT = Material_Name Elem_Blks" line
     *   mn  == Current material index
     *   input == temporary character buffer of length MAX_CHAR_IN_LINE.
     *     
     ************************************************************************/
{
  int numTok, i;
  TOKEN_STRUCT tok;
  char *yo = "read_MAT_line";
  char echo_string[MAX_CHAR_IN_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  /*
   *  Search forward in the input deck for the next MAT keyword line
   */
  look_for(ifp, "MAT", input, '=');

  /*
   *  Store the remainder of the line in the character buffer, then tokenize
   *  it.
   */
  read_line(ifp, input, FALSE);
  numTok = fillTokStruct(&tok, input);

  /*
   * Check for a minimum number of parameters on the line
   */
  if (numTok < 2) {
    (void) fprintf(stderr,
		   "%s error: Need at least two arguments for MAT keyword\n%s\n",
		   yo, input);
    EH(-1, yo);
  }

  if( numTok == MAXTOKENS )
    {
      (void) fprintf(stderr, "%s error: Max tokens read parsing line\n %s \n", yo, input);
      (void) fprintf(stderr, "Recompile with larger value for MAXTOKENS.\n");
      EH(-1,yo);
    }

  /*
   * Store the material index into the material structure
   */
  mp_glob[mn]->MatID = mn;
  
  /*
   * Parse the material name
   */
  if (! interpret_string(tok.tok_ptr[0], pd_glob[mn]->MaterialName))
      goto L_ERROR;
  strcpy(mp_glob[mn]->Material_Name, pd_glob[mn]->MaterialName);

  SPF(echo_string,"%s = %s","MAT",pd_glob[mn]->MaterialName);

  /*
   * Store the number of element blocks having this material type
   *  -> Malloc space in the material structure based on this quantity
   */
  mp_glob[mn]->Num_Matrl_Elem_Blk = numTok - 1;
  mp_glob[mn]->Matrl_Elem_Blk_Ids =
      alloc_int_1(mp_glob[mn]->Num_Matrl_Elem_Blk, 0);
  
  /*
   * Store the element block id's of the element blocks in the list
   *  -> Note, we can't formulate  mp_glob[mn]->Matrl_Elem_Blks[] yet, since
   *     we haven't read the exodus file and thus haven't set up
   *     the Element Block structure list.
   */
  for (i = 0; i < mp_glob[mn]->Num_Matrl_Elem_Blk; i++) {
    if (! interpret_int(tok.tok_ptr[i+1],
			mp_glob[mn]->Matrl_Elem_Blk_Ids + i))
      goto L_ERROR;

    SPF_INT_VEC( endofstring(echo_string), mp_glob[mn]->Num_Matrl_Elem_Blk, mp_glob[mn]->Matrl_Elem_Blk_Ids);

    ECHO(echo_string, echo_file);
  }
  return;
  /*
   * Error exit
   */
 L_ERROR:;
  (void) fprintf(stderr, 
		 "%s error: reading MAT line, %s\n", yo, input);

  EH(-1, "MAT line");
}

void
read_surface_objects ( FILE* ifp,
		       char *input,
		       struct LS_Surf_List *list, 
		       int num_surf)

{
  int iread;
  char name[10];
  struct LS_Surf *surf;
  char echo_string[MAX_CHAR_IN_INPUT]="\0";
  char *echo_file = Echo_Input_File;

  while( num_surf > 0 )
    {
      iread = look_forward_optional( ifp, "SURF", input, '=');

      if ( iread == -1 )
	{
	  EH(-1,"Not enough Level Set initialization surfaces.");
	}

      if( fscanf(ifp, "%s", name ) != 1 )
	{
	  EH(-1,"Error reading Level Set initialization surface.");
	}

      stringup(name);

      if( strcmp( name, "PLANE") == 0 )
	{
	  struct LS_Surf_Plane_Data *s;
          
          surf = create_surf( LS_SURF_PLANE );
          append_surf( list, surf );
          s = (struct LS_Surf_Plane_Data *) surf->data;

	  if( fscanf( ifp, "%lf %lf %lf %lf", &(s->n[0]),
		                              &(s->n[1]),
		                              &(s->n[2]),
                                              &(s->d)) != 4 )
	    {
	      EH(-1,"PLANE Level Set initialization objects requires 4 float constants.");
	    }

	  SPF(echo_string,"\t\t%s = %s %.4g %.4g %.4g %.4g","SURF","PLANE",s->n[0],s->n[1],s->n[2],s->d);ECHO(echo_string,echo_file);

	}
      else if (  strcmp( name, "CIRCLE") == 0 )
	{
	  struct LS_Surf_Sphere_Data *s;
          
          surf = create_surf( LS_SURF_CIRCLE );
          append_surf( list, surf );
          s = (struct LS_Surf_Sphere_Data *) surf->data;
          
          s->center[2] = 0.;
          
	  if( fscanf( ifp, "%lf %lf %lf", &(s->center[0]),
		                          &(s->center[1]),
		                          &(s->r) ) != 3 )
	    {
	      EH(-1,"CIRCLE Level Set initialization objects requires 3 float constants.");
	    }
	  SPF(echo_string,"\t\t%s = %s %.4g %.4g %.4g","SURF","CIRCLE",s->center[0],s->center[1],s->r ); ECHO(echo_string,echo_file);
	}
      else if (  strcmp( name, "SPHERE") == 0 )
	{
	  struct LS_Surf_Sphere_Data *s;
          
          surf = create_surf( LS_SURF_SPHERE );
          append_surf( list, surf );
          s = (struct LS_Surf_Sphere_Data *) surf->data;
          
	  if( fscanf( ifp, "%lf %lf %lf %lf", &(s->center[0]),
		                              &(s->center[1]),
                                              &(s->center[2]),
		                              &(s->r) ) != 4 )
	    {
	      EH(-1,"SPHERE Level Set initialization objects requires 4 float constants.");
	    }

	  SPF(echo_string,"\t\t%s = %s %.4g %.4g %.4g %.4g","SURF","SPHERE",s->center[0],s->center[1],s->center[2],s->r); 
	  ECHO(echo_string,echo_file);
	}
      else if (  strcmp( name, "FACET_LIST") == 0 )
	{
	  int i, num, num_points, inflection;
          double r[3] = {0.,0.,0.};
          double xi[3] = {0.,0.,0.};
          struct LS_Surf *vertex[3], *facet;
          FILE *sfp;

	  (void) read_string(ifp,input,'\n');
          strip(input);
          
          if ( (sfp=fopen(input,"r")) == NULL) 
            EH( -1, "Could not open facet list file\n");

	  SPF(echo_string,"\t\t%s = %s %s","SURF","FACET_LIST",input);ECHO(echo_string,echo_file);

          while( look_forward_optional(sfp, "FACET", input, '=') == 1 )
            {
              if( fscanf( sfp, "%d %d", &num,
		                        &num_points) != 2 )
	        {
	          EH(-1,"Error reading FACET in facet file.");
	        }
                
              for (i=0; i<num_points; i++)
                {
                  if (look_forward_optional(sfp, "POINT", input, '=') == -1)
                    EH(-1,"Error finding POINT for FACET in facet file.");
                    
                  if( fscanf( sfp, "%d %lf %lf %lf %d", &num,
		                                        &(r[0]),
                                                        &(r[1]),
                                                        &(r[2]),
                                                        &inflection) != 5)
	            EH(-1,"Error reading POINT for FACET in facet file.");
                  
                  vertex[i] = create_surf_point ( r, -1, xi, inflection );
                }
                
              if (num_points == 2)
                {
                  facet = create_surf_facet_line( vertex[0], vertex[1], -1 , -1 );
                  append_surf( list, facet );
                }
              else if (num_points == 3)
                {
                  EH( -1, "Not yet implemented in 3-D\n");
                }
            }
          fclose( sfp );
	}
      else if (  strcmp( name, "SS") == 0 )
	{
	  struct LS_Surf_SS_Data *s;
          
          surf = create_surf( LS_SURF_SS );
          append_surf( list, surf );
          s = (struct LS_Surf_SS_Data *) surf->data;

	  if( fscanf( ifp, "%d", &(s->ss_id)) != 1 )
	    {
	      EH(-1,"SS Level Set initialization objects requires 1 int constant.");
	    }

	  SPF(echo_string,"\t\t%s = %s %d","SURF","SS",s->ss_id);ECHO(echo_string,echo_file);
	}
      else if (  strcmp( name, "ISOSURFACE") == 0 )
	{
	  struct LS_Surf_Iso_Data *s;
          int var_found, k;
          
          surf = create_surf( LS_SURF_ISOSURFACE );
          append_surf( list, surf );
          s = (struct LS_Surf_Iso_Data *) surf->data;

          if (fscanf(ifp, "%80s", input) != 1)
	    {
	      EH(-1, "Error reading variable type for Level Set initialization ISOSURFACE");
	    }

	    /* loop through variable names and compare to input */
	
	    var_found = FALSE;
	    for (k=0; k<MAX_VARIABLE_TYPES + 1; k++) /* only check variables (including porosity) */
	      {
	  	if ( !strcmp(input, Var_Name[k].name1) ||
	  	     !strcmp(input, Var_Name[k].name2))
	  	  {
	  	    s->isovar = Var_Name[k].Index;
	  	    var_found = TRUE;
	  	  }
	      }
	    if (!var_found)
	      {
	  	EH(-1, "Could not identify variable type for Level Set initialization ISOSURFACE");
	      }
              
	  if( fscanf( ifp, "%lf", &(s->isoval) ) != 1 )
	    {
	      EH(-1,"ISOSURFACE Level Set initialization objects requires varname and 1 float constant.");
	    }

	  SPF(echo_string,"\t\t%s = %s %s %.4g","SURF","ISOSURFACE",input,s->isoval);ECHO(echo_string,echo_file);

	}
      else if ( strcmp( name, "ARC" ) == 0 )
	{
	  struct LS_Surf_Arc_Data *s;
	  
	  surf = create_surf( LS_SURF_ARC );
	  append_surf( list, surf );
	  
	  s = (struct LS_Surf_Arc_Data *) surf->data;
		  
	  if( fscanf( ifp, "%lf %lf %lf %lf %lf %lf", &(s->center[0]), &(s->center[1]), &(s->r),
		      &(s->n[0]), &(s->n[1]), &(s->d) ) != 6) 
	    {
	      EH(-1," ARC initialization surface requires six floats,\n");
	    }  

	  SPF(echo_string,"\t\t%s = %s %.4g %.4g %.4g %.4g %.4g %.4g ","SURF","ARC",
	      s->center[0],s->center[1], s->r, s->n[0],s->n[1], s->d ); 

	 s->sign = 1.0;

	 read_string(ifp,input,'\n'); strip(input);

	 if( strlen(input) != 0 )
	    {
	      stringup(input);
	      
	      if( strcmp(input,"NEGATIVE") == 0 )
		{
		  s->sign = -1.0;
		}
	      else if ( strcmp(input,"POSITIVE") != 0 )
		{
		  EH(-1,"Sign designation string for ARC type of LS initialization must be POSITIVE or NEGATIVE");
		}
	      SPF(endofstring(echo_string)," %s", input);
	    }
	  
	  
	  ECHO(echo_string,echo_file);
	}
      else if (  strcmp( name, "USER") == 0 )
	{
	  int I;
	  double *a;
	  struct LS_Surf_User_Data *s;
	  char *s1;
          
          surf = create_surf( LS_SURF_USER );
          append_surf( list, surf );
          s = (struct LS_Surf_User_Data *) surf->data;

	  I = s->Int_Data[0] = read_constants( ifp, 
					       &a,
					       0 );

	  while ( I-- != 0 ) s->Real_Data[I] = a[I];

	  safe_free(a);
	  SPF(echo_string,"\t\t%s = %s","SURF","USER");

	  I = 0;

	  while( I < s->Int_Data[0])
	    {
	      s1 = endofstring(echo_string);
	      SPF(s1," %.4g",s->Real_Data[I]);
	      I++;
	    }
	  ECHO(echo_string,echo_file);
	}

      else
	{
	  EH(-1,"Cannot find Level Set initialization surface type.");
	}

      num_surf--;

    }

}

#if 0
void
echo_surface_objects( struct LS_Surf_List *list )
{
  struct LS_Surf *surf;

  surf = list->start;

  while( surf != NULL)
    {
      switch( surf->type )
	{

	case LS_SURF_PLANE:
	  {
	    struct LS_Surf_Plane_Data *s = (struct LS_Surf_Plane_Data *) surf->data;
	    SPF(echo_string,"\t\t%s  %f %f %f %f","PLANE", s->n[0],s->n[1],s->n[2], s->d ); ECHO(echo_sting);
	  }
	  break;
	case LS_SURF_CIRCLE:
	  {
	    struct LS_Surf_Sphere_Data *s = (struct LS_Surf_Sphere_Data *) surf->data;

	    SPF(echo_string,"\t\t%s  %f %f %f","CIRCLE",s->center[0], s->center[1], s->r ); ECHO(echo_string,echo_file)
	  }
	    break;
	case LS_SURF_SPHERE:
	  {
	    struct LS_Surf_Sphere_Data *s = (struct LS_Surf_Sphere_Data *) surf->data;
	    SPF(echo_string,"\t\t%s %f %f %f %f","SPHERE", s->center[0], s->center[1], s->center[2], s->r ); ECHO(echo_string,echo_file);
	  }
	  break;

	case LS_SURF_ARC:
	  {
	    struct LS_Surf_Arc_Data *s = (struct LS_Surf_Arc_Data *) surf->data;
	    SPF(echo_string,"\t\t%s %f %f %f %f %f %f", "ARC", s->center[0], s->center[1],  s->r, s->n[0], s->n[1], s->d ); ECHO(echo_string,echo_file);
	  }
	  break;

	case LS_SURF_POINT:
	  {
	    struct LS_Surf_Point_Data *s = (struct LS_Surf_Point_Data *) surf->data;

	    SPF(echo_string,"\t\t%s %f %f %f ","POINT", s->x[0], s->x[1], s->x[2] ); ECHO(echo_string,echo_file);
	  }

	  break;

	case LS_SURF_FACET:
	  ECHO("/t/tFACET", echo_file);
	  break;

	case LS_SURF_SS:
	  {
	    struct LS_Surf_SS_Data *s = (struct LS_Surf_SS_Data *) surf->data;
	    SPF(echo_string,"\t\t%s  %d","SS", s->ss_id ); ECHO(echo_string,echo_file);
	  }
	  break;

	case LS_SURF_ISOSURFACE:
	  {
	    int k=0;
	    struct LS_Surf_Iso_Data *s = (struct LS_Surf_Iso_Data *) surf->data;

	    while( s->isovar != Var_Name[k].Index && k<MAX_VARIABLE_TYPES ) k++;
	    
	    SPF(echo_string,"\t\t%s %s %.4g","ISOSURFACE", Var_Name[k].name1, s->isoval ); ECHO(echo_string,echo_file)  ;
	  }
	  break;
	case LS_SURF_NS:
	  {
	    struct LS_Surf_NS_Data *s = (struct LS_Surf_NS_Data *) surf->data;
	    DPRINTF(echo,
		    "Level Set Initialization Method -------> Surface:NODESET %s \n", 
		    s->ns_id );
	  }
	  break;

	}
      surf = surf->next;
    }
}
		      
#endif

struct Data_Table *
setup_table_external(
	       char *filename, 
	       struct Data_Table * table,
	       char *var_name)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  char *yo;
  FILE *datafile =NULL;
  yo = "setup_table_external";

  table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) ) ;

  table->columns = pd->Num_Dim+1;

  /*
   * "y-axis" of table
   */

  table->f_name = var_name;
  table->f_index = 0;

  /* read interpolation order */

  switch(pd->Num_Dim)
	{
	case 2:
  		table->interp_method = BILINEAR;
		break;
	case 3:
  		table->interp_method = TRILINEAR;
		break;
	}


  if ( ( datafile = fopen_aprepro( filename, "r") ) == NULL )
	{
	  sprintf(err_msg, "\n%s,\tError opening TABLE data file %s\n",yo,filename);
	  EH(-1,err_msg);
	}

  rd_table_data( datafile, filename, table, NULL );
  fclose( datafile );

return( table );
  
}
/*   end of setup_table_external	*/

struct Data_Table *
setup_table_AC(
	       char *filename, 
	       struct Data_Table * table,
	       char *var_name,
	       char *interpolation)
{
  char err_msg[MAX_CHAR_IN_INPUT];
  char *yo;
  FILE *datafile =NULL;
  yo = "setup_table_AC";

  table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) ) ;

  /* 
   * "x-axis" of table
   */

  if ( strcmp( var_name, "TIME" ) == 0 )
    {
      strcpy( table->t_name[0],"TIME");
      table->t_index[0] = -1;
      table->columns = 2;
    }
  else if( strcmp( var_name, "X") == 0 )
    {
      strcpy( table->t_name[0],"X");
      table->t_index[0] = 0;
      table->columns = 2;
    }
  else if ( strcmp( var_name, "Y") == 0 )
    {
      strcpy( table->t_name[0],"Y");
      table->t_index[0] = 1;
      table->columns = 2;
    }
  else if ( strcmp( var_name, "Z") == 0 )
    {
      strcpy( table->t_name[0],"Z");
      table->t_index[0] = 2;
      table->columns = 2;
    }
  else if ( strcmp( var_name, "XY") == 0 )
    {
      strcpy( table->t_name[0],"X");
      table->t_index[0] = 0;
      strcpy( table->t_name[1],"Y");
      table->t_index[1] = 1;
      table->columns = 3;
    }
  else if ( strcmp( var_name, "XZ") == 0 )
    {
      strcpy( table->t_name[0],"X");
      table->t_index[0] = 0;
      strcpy( table->t_name[1],"Z");
      table->t_index[1] = 2;
      table->columns = 3;
    }
  else if ( strcmp( var_name, "YZ") == 0 )
    {
      strcpy( table->t_name[0],"Y");
      table->t_index[0] = 1;
      strcpy( table->t_name[1],"Z");
      table->t_index[1] = 2;
      table->columns = 3;
    }
  else if ( strcmp( var_name, "YX") == 0 )
    {
      strcpy( table->t_name[0],"Y");
      table->t_index[0] = 1;
      strcpy( table->t_name[1],"X");
      table->t_index[1] = 0;
      table->columns = 3;
    }
  else if ( strcmp( var_name, "ZX") == 0 )
    {
      strcpy( table->t_name[0],"Z");
      table->t_index[0] = 2;
      strcpy( table->t_name[1],"X");
      table->t_index[1] = 0;
      table->columns = 3;
    }
  else if ( strcmp( var_name, "ZY") == 0 )
    {
      strcpy( table->t_name[0],"Z");
      table->t_index[0] = 2;
      strcpy( table->t_name[1],"Y");
      table->t_index[1] = 1;
      table->columns = 3;
    }
  else
    {
      sprintf(err_msg,"\nInvalid choice for table abscissa.");
      EH(-1,err_msg);
    }

  /*
   * "y-axis" of table
   */

  table->f_name = "table_variable";
  table->f_index = 0;

  /* read interpolation order */

  if ( strcmp( interpolation, "LINEAR") == 0 )
    {
      table->interp_method = LINEAR;
    }
  else if ( strcmp( interpolation, "QUADRATIC") == 0 )
    {
      table->interp_method = QUADRATIC;
    }
  else if ( strcmp( interpolation, "QUAD_GP") == 0 )
    {
      table->interp_method = QUAD_GP;
    }
  else if ( strcmp( interpolation, "BIQUADRATIC") == 0 )
    {
      table->interp_method = BIQUADRATIC;
    }
  else if ( strcmp( interpolation, "BILINEAR") == 0 )
    {
      table->interp_method = BILINEAR;
    }
  else
    {
      sprintf(err_msg, "\nUnknown table interpolation order for table: %s \n",interpolation);
      EH(-1,err_msg);
    }

  if ( ( datafile = fopen_aprepro( filename, "r") ) == NULL )
	{
	  sprintf(err_msg, "\n%s,\tError opening TABLE data file %s\n",yo,filename);
	  EH(-1,err_msg);
	}

  rd_table_data( datafile, filename, table, NULL );
  fclose( datafile );

return( table );
  
}
/*   end of setup_table_AC	*/


FILE *
handle_echo_string( char *echo_string, char *echo_file_name )
{
  static FILE *echo_file_ptr[FOPEN_MAX];
  static char *open_files[FOPEN_MAX];
  static int num_open_files=0;
  static int echo_off = FALSE;

  int i=0,len=0;
  if( (len=strlen(echo_string)) == 0 ) return NULL;
  if( echo_off == TRUE ) return NULL;
  if( len > MAX_CHAR_IN_INPUT)
    {
      DPRINTF(stderr,"Error: size of input card \n >>%s<< \n\t has exceeded MAX_CHAR_IN_INPUT",echo_string);
      exit(-1);
    }  

#ifdef DISABLE_ECHO
  return(NULL);
#else


  while ( (i<num_open_files) && (strcmp(open_files[i], echo_file_name) != 0) ) i++;

  if(strcmp(echo_string,"OPEN") == 0 )
    {
      if( i < num_open_files ) {
	fprintf(stderr,"Echo file %s already opened\n", echo_file_name);
	exit(-1);
      }
      if( (echo_file_ptr[i] = fopen(echo_file_name,"w")) == NULL ) {
	fprintf(stderr,"%s %s\n","Cannot open echo file for writing :", echo_file_name);
	exit(-1);
      }
      open_files[i] = alloc_copy_string(echo_file_name);
      num_open_files++;
      return( echo_file_ptr[i] );
    }
  else if ( strcmp( echo_string,"NOECHO") == 0 )
    {
      echo_off = TRUE;
      return( NULL );
    }
  else if ( i == num_open_files )
    {
      fprintf(stderr,"Echo file %s was never opened.\n",echo_file_name );
      return( NULL );
      /* exit(-1); */
    }
  else if ( strcmp(echo_string, "FILE_POINTER") == 0 )
    {
      return( echo_file_ptr[i] );
    }
  else if ( strcmp(echo_string,"CLOSE") == 0 )
    {
      fflush(echo_file_ptr[i]);
      fclose( echo_file_ptr[i] );
      safer_free((void **) &(open_files[i]) );

      while ( i + 1 < num_open_files )
	{
	  echo_file_ptr[i] = echo_file_ptr[i+1];
	  open_files[i] = open_files[i+1];
	}

      echo_file_ptr[i] = NULL;
      num_open_files--;
    }
  else
    {
      DPRINTF( echo_file_ptr[i],"%s\n", echo_string);
    }

  return(NULL);
#endif
}



int
sprintf_int_vec( char *string, int n, int *v)
{
  int i=0;
  char *s1 = string;

  while(i<n)
    {
      sprintf(s1," %d",v[i]);
      s1 = strchr(string,'\0');
      i++;
    }
  return(0);
}



int
sprintf_flt_vec( char *string, int n, float *v)
{
  int i=0;
  char *s1 = string;

  while(i<n)
    {
      sprintf(s1," %.4g",v[i]);
      s1 = strchr(string,'\0');
      i++;
    }
  return(0);
}

int
sprintf_dbl_vec( char *string, int n, double *v)
{
  int i=0;
  char *s1 = string;

  while(i<n)
    {
      sprintf(s1," %.4g",v[i]);
      s1 = strchr(string,'\0');
      i++;
    }
  return(0);
}


/*****************************************************************************/
/*****************************************************************************/
