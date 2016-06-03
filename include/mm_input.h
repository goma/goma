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
 *$Id: mm_input.h,v 5.3 2008-12-19 22:54:26 rbsecor Exp $
 */

#ifndef _MM_INPUT_H
#define _MM_INPUT_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_INPUT_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_INPUT_C
#define EXTERN extern
#endif

#include <stdio.h>

#ifndef MAX_INPUT_LINE_LENGTH
#define MAX_INPUT_LINE_LENGTH 2048
#endif

/*
 *  Definition of a Token structure
 */

#ifndef MAXTOKENS
#  define MAXTOKENS 30
#endif
struct Token {
  char  orig_str[MAX_CHAR_IN_INPUT + 1];
  char  tok_str[MAX_CHAR_IN_INPUT + 1];
  char *tok_ptr[MAXTOKENS];
  int   ntokes;
};
typedef struct Token TOKEN_STRUCT;

#ifndef SPF
#define SPF sprintf
#endif

#ifndef ECHO
#define ECHO handle_echo_string
#endif

extern void read_input_file	/* mm_input.c                                */
PROTO((struct Command_line_command **, /*clc                                */
       int ));			/* nclc                                      */

extern int read_line	        /* mm_input_util.c                           */
PROTO((FILE *,			/* ifp                                       */
       char [],			/* string                                    */
       const int ));		/* print_flag                                */

extern int read_string		/* mm_input_util.c                           */
PROTO((FILE *,			/* ifp                                       */
       char [],			/* string                                    */
       const char ));		/* ch                                        */

extern int strip		/* mm_input_util.c                           */
PROTO((char []));		/* string                                    */

extern int count_parameters	/* mm_input.c                                */
PROTO((const char []));		/* string                                    */

extern void look_for		/* mm_input.c                                */
PROTO((FILE  *,			/* ifp                                       */
       const char *,		/* string                                    */
       char [],			/* input                                     */
       const char ));		/* ch_term                                   */

extern int look_for_either	/* mm_input.c                                */
PROTO((FILE  *,			/* ifp                                       */
       const char *,		/* string1                                   */
       const char *,		/* string2                                   */
       char [],			/* input                                     */
       const char ));		/* ch_term                                   */

extern int count_list		/* mm_input.c                                */
PROTO((FILE  *,			/* ifp - file pointer, open goma input file  */
       const char *,		/* string -                                  */
       char [],			/* input -                                   */
       const char ,		/* ch_term                                   */
       const char *));		/* stringend                                 */

extern int look_for_n_doubles(FILE *ifp, int n, double *buf);

extern int look_for_optional	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       const char *,		/* string                                    */
       char [],			/* input                                     */
       const char ));		/* ch_term                                   */

extern int look_forward_optional /* mm_input.c                               */
PROTO((FILE *,			/* ifp                                       */
       const char *,		/* string                                    */
       char [],			/* input                                     */
       const char ));		/* ch_term                                   */

extern void rd_file_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char * ));  	        /* input                                     */

extern void rd_genl_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char * ));  		/* input                                     */

extern void rd_timeint_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char * ));  		/* input                                     */

extern void rd_levelset_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char * ));  		/* input                                     */

extern void rd_elem_quality_specs	/* mm_input.c                        */
PROTO((FILE *,			/* ifp                                       */
       char * ));   		/* input                                     */

extern void rd_track_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_hunt_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_ac_specs		/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_solver_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));   		/* input                                     */

extern void rd_eigen_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_geometry_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_bc_specs		/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_matl_block_specs	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *));		/* input                                     */

extern void rd_eq_specs		/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *,			/* input                                     */
       const int));		/* mn                                        */

extern void rd_mp_specs		/* mm_input.c                                */
PROTO((FILE *,			/* imp - input mat file ptr (open)           */
       char *,			/* input                                     */
       const int ,		/* mn                                        */
      char * ));      /*          echo file name for this mat               */
      
extern int look_for_mat_prop	/* mm_input.c                                */
PROTO((FILE *,			/* imp - ptr to input stream            (in) */
       char *,		        /* search_string -                      (in) */
       int *,			/* MaterialModel - int material model  (out) */
       dbl *,			/* Material_property - double value matprop, *
				 * if constant                         (out) */
       dbl **,			/* User_constants - ptr to vector of double  *
				 * constants for user defined material       *
				 * properties                          (out) */
       int *,			/* User_count - how many of the above        *
				 * User_constants ?                    (out) */
       char *,			/* model_name - ptr to model name      (out) */
       const int ,		/* num_values - num. values for constant     *
				 * models SCALAR=1, VECTOR=3            (in) */
		int * ,			/* read_species - flag to read species or    *
				 * not comes in as Max_Number species but    *
				 * output at species number of property that *
				 * is read.                             (in) */
       char *));            /*     pointer to echo string character array    */

EXTERN int look_for_mat_proptable	/* mm_input.c                                */
PROTO((FILE *,			/* imp - ptr to input stream            (in) */
       char *,		        /* search_string -                      (in) */
       int *,			/* MaterialModel - int material model  (out) */
       dbl *,			/* Material_property - double value matprop, *
				 * if constant                         (out) */
       dbl **,			/* User_constants - ptr to vector of double  *
				 * constants for user defined material       *
				 * properties                          (out) */
       int *,			/* User_count - how many of the above        *
				 * User_constants ?                    (out) */
       int *,                   /* table counter */
       char *,			/* model_name - ptr to model name      (out) */
       const int ,		/* num_values - num. values for constant     *
				 * models SCALAR=1, VECTOR=3            (in) */
		int *,			/* read_species - flag to read species or    *
				 * not comes in as Max_Number species but    *
				 * output at species number of property that *
				 * is read.                             (in) */
       char *));          /*      pointer to echo string                    */

extern int look_for_modal_prop	/* mm_input.c                                */
PROTO((FILE *,			/* imp - ptr to open input stream       (in) */
       const char *,		/* search_string -                      (in) */
       const int ,		/* modes - number of viscoelastic modes (in) */
       int *,			/* MaterialModel - int material model  (out) */
	   dbl * ,			/* modal_const - modal data            (out) */
       char * ));        /*      echo string char array              (out) */

extern int read_constants	/* mm_input.c                                */
PROTO((FILE *,			/* imp - pointer to (open) file stream       */
       dbl **,			/* User_constants - array of user mat props  */
       const int ));		/* species_no - species number (zero if no   * 
				 * species)                                  */

extern void set_mp_to_unity	/* mm_input.c                                */
PROTO((const int ));		/* mn                                        */

extern void usage		/* mm_input.c                                */
PROTO((const int ));		/* exit_flag                                 */

extern void translate_command_line /* mm_input.c                             */
PROTO((int ,			/* argc                                      */
       char *[],		/* argv                                      */
       struct Command_line_command **, /* clc                                */
       int *));			/* nclc                                      */

extern void apply_command_line	/* mm_input.c                                */
PROTO((struct Command_line_command **, /* clc                                */
       int ));			/* nclc                                      */

extern struct Data_Table *setup_table_BC /* mm_input.c                       */
PROTO((FILE *,			/* ifp                                       */
       char *,			/* input                                     */
       struct Boundary_Condition *, /* BC_Type                             */
		char * ));      /*   echo string */ 

EXTERN struct Data_Table *setup_table_MP  /*mm_input.c                    */
PROTO((FILE *,			 /*ifp                                       */
       struct Data_Table *,      /*mp                                        */
       char *));		 /*search string                             */

EXTERN struct Data_Table *setup_table_external /* mm_input.c              */
PROTO((char *,			/* filename                          */
       struct Data_Table *,	/* table structure			*/
       char * ));		/* table variable name			*/

EXTERN struct Data_Table *setup_table_AC /* mm_input.c              */
PROTO((char *,			/* filename                          */
       struct Data_Table *,	/* table structure			*/
       char *,			/* table variable name			*/
       char * ));		/* table interpolation			*/

extern void rd_matl_blk_specs   /* mm_input.c                                */
PROTO((FILE *,                  /* ifp                                       */
       char *));                /* input                                     */

extern void rd_table_data	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *,			/* input                                     */
       struct Data_Table * ,	/* table                                     */
       char *));		/* endlist                                   */

extern int scan_table_columns
PROTO(( int,                  /* table row index */
		char *,               /* input character array to be scanned */
		struct Data_Table *,  /* data table pointer */
		int ,	              /* table dimension    */
		char *,               /* error string */
		char *));             /* echo string */

extern int count_datalines	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       char *,			/* input                                     */
       const char *));		/* endlist                                   */

extern int look_for_next_string	/* mm_input.c                                */
PROTO((FILE *,			/* ifp                                       */
       const char *,		/* string                                    */
       char *,			/* input                                     */
       const char ));		/* ch_term                                   */

extern FILE *fopen_aprepro	/* mm_input.c                                */
PROTO((const char *,		/* filename                                  */
       const char *));		/* format                                    */

extern struct Data_Table *setup_gd_table_BC /* mm_input.c                    */
PROTO((FILE *,			/* ifp                                       */
       char *,			/* input                                     */
       struct Boundary_Condition *,  /* BC_Type                             */
       char * ));                    /* echo string */

extern void echo_compiler_settings		/* mm_input.c                                */
PROTO(( void ));

extern int look_for_optional_string	/* mm_input_mp.c                     */
PROTO((FILE *,			        /* imp - ptr to open input stream       (in) */
       const char *,	                /* search_string -                      (in) */
       char *,   	                /* retn_string - output string         (out) */
       int ));                          /* Length or retn_string                (in) */


extern void read_surface_objects 
PROTO (( FILE *,
	 char *, 
	 struct LS_Surf_List *,
	 int ));

extern void echo_surface_objects
PROTO (( FILE *,
	 struct LS_Surf_List * ));

extern FILE *handle_echo_string
PROTO (( char *,    /* echo_string */
	 char * )); /* base filename */


extern int sprintf_int_vec 
PROTO (( char *,           /* character array to recieve input */
	 int,              /* number of elements in vector */
	 int * ));         /* vector */

#define SPF_INT_VEC  sprintf_int_vec

extern int sprintf_flt_vec 
PROTO (( char *,           /* character array to recieve input */
	 int,              /* number of elements in vector */
	 float * ));       /* vector */

#define SPF_FLT_VEC  sprintf_flt_vec

extern int sprintf_dbl_vec 
PROTO (( char *,           /* character array to recieve input */
	 int,              /* number of elements in vector */
	 double * ));      /* vector */

#define SPF_DBL_VEC  sprintf_dbl_vec


/****************************************************************************
 *  Prototypes for functions in mm_input_util.c
 ****************************************************************************/

extern int look_for_optional_int
PROTO((FILE *, const char *, int *, const int));

extern int read_int          
PROTO((FILE *,                  /* ifp                                       */
       const char *));          /* intName                                   */

extern double read_dbl      
PROTO((FILE *,                  /* ifp                                       */
       const char *));          /* intName                                   */

extern int read_1_boolean      
PROTO((FILE *,                  /* ifp                                       */
       const char *));          /* boolName                                  */

extern int variable_string_to_int 
PROTO((const char *,            /* input                                      */
       const char *));          /* err_string                                 */

extern int species_type_str_to_int
PROTO((char *));                /* input                                      */

extern void species_type_int_to_str
PROTO((char *,                  /* str                                        */
       const int));             /* var                                        */

extern int interpret_string
PROTO((const char *,            /* string                                     */
       char *));                /* return string                              */

extern int interpret_int
PROTO((const char *,            /* string                                     */
       int *));                 /* return_value                               */

extern int interpret_double
PROTO((const char *,            /* string                                     */
       double *));              /* return_value                               */

extern int indentify_species_ID_string
PROTO((const char *,            /* input                                      */
       const char *,            /* list                                       */
       const int));             /* numList                                    */

extern int look_for_species_prop
PROTO((FILE *,                 /* imp                                         */
       const char *,           /* search_string                               */
       MATRL_PROP_STRUCT *,    /* mat_ptr                                     */
       int  *,                 /* materialModel                               */
       dbl  *,                 /* material_property                           */	
       dbl  **,                /* User_constants                              */
       int  *,                 /* User_count                                  */		
       char *,                 /* model_name                                  */	
       const int,              /* num_values_expected                         */   
       int * ,                  /* species_index                               */
       char *));               /* echo string                                 */

extern int look_for_porous_prop
PROTO((FILE *,                 /* imp                                         */
       const char *,           /* search_string                               */
       MATRL_PROP_STRUCT *,    /* mat_ptr                                     */
       int  *,                 /* materialModel                               */
       dbl  *,                 /* material_property                           */	
       dbl  **,                /* User_constants                              */
       int  *,                 /* User_count                                  */		
       char *,                 /* model_name                                  */	
       const int,              /* num_values_expected                         */   
       int * ,                  /* porous_phase_index                          */
       char *));             /*     echo string                                 */

extern int look_for_species_proptable
PROTO((FILE *,                 /* imp                                         */
       char *,                 /* search_string                               */
       MATRL_PROP_STRUCT *,    /* mat_ptr                                     */
       int  *,                 /* materialModel                               */
       dbl  *,                 /* material_property                           */	
       dbl  **,                /* User_constants                              */
       int  *,                 /* User_count                                  */		
       int *,                   /* table counter */
       char *,                 /* model_name                                  */	
       const int,              /* num_values_expected                         */   
       int * ,                  /* species_index                               */
       char * ));           /*    echo string                                 */

extern int stokenize
PROTO((char *,                  /* string                                     */
       const char *,            /* delimiters                                 */
       char **,                 /* token_pts                                  */
       const int));             /* token_pts                                  */

extern int tokenize_by_whsp
PROTO((char *,                  /* string                                     */
       char **,                 /* token_pts                                  */
       const int));             /* max_toks                                   */

extern int fillTokStruct
PROTO((TOKEN_STRUCT *,          /* tp                                         */
       const char *));          /* string                                     */

extern int in_char_list
PROTO((const char *,            /* target                                     */
       const char *,            /* list                                       */
       const int));             /* num_list                                   */

    
#endif /* _MM_INPUT_H */
