/************************************************************************\
 * Copyright (c) 2002 Sandia Corporation.                               *
 *					         			*
 * Under the terms of Contract DE-AC04-94AL85000, there is a            *
 * non-exclusive license for use of this work by or on behalf of the    *
 * U.S. Government. Export of this program may require a license from   *
 * the United States Government.                                        *
 *					         			*
 * This software is the property of Sandia Corporation and discloses    *
 * material protectable under copyright laws of the United States.      *
 * Use, Duplication, or Disclosure is prohibited, unless allowed        *
 * subject to the terms of a separate license agreement.                *
 *					         			*
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES          *
 * DEPARTMENT OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR       *
 * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY    *
 * LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,    *
 * OR USEFULNESS OF ANY INFORMATION, APPARATUS OR PROCESS DISCLOSED,    *
 * OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.*
 *					         			*
\************************************************************************/

%{



#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>	
#include <time.h>
	
#include "rf_io_const.h"	
#include "std.h"
#include "mm_input.h"	
#include "rf_vars_const.h"
#include "rf_io_const.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_solver.h"
#include "rf_mp.h"
#include "rf_io_structs.h"
/*#include "rf_io.h"*/
#include "rf_bc_const.h"
#include "rf_allo.h"
#include "rf_bc.h"

#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_as_alloc.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_mp_const.h"

#include "mm_eh.h"

#include "mm_post_def.h"

#include "sl_util_structs.h"
#include "mm_input.h"

/*#include "goma.h"*/
/*#include "rf_io_defn.h"*/
#include "mm_species.h"
#include "rd_mesh.h"
#include "mm_parser.h"
 
extern char yytext[];
extern struct cards;

#undef twrap()

/***********************/
/* Parser Declarations */
/***********************/

/* declare variables between these lines */

float  bc_float_array[MAX_NUMBER_PARAMS];
int    bc_float_array_index = 0;
char * bcstring;
FILE * bc_file;
extern int accept_table_data_cards;
extern float	table_array[100][5];
extern int	number_of_table_cards_read;
extern int	local_BC_index;
int    print_strings_to_log = FALSE;
extern FILE * parser_log;
extern int line_number;
char   tag[50];
int    columns_expected;
char   msg[100];
int reading_from_table_file = 0;
extern int error_found_in_last_line;
extern int error_number;
extern float floating_point_constant_list_array[MAX_NUMBER_PARAMS];
extern int  floating_point_constant_list_index;
int 	k;
int	tag_needed;

/* end the declarations ehre */

#define YYERROR_VERBOSE 1 	

%}

%union	{
	char 	string[100];
	int 	integer;
	char	character;
	int	integer_array[50];
	float	floating;
	}
	
/* Parser File Mode Tokens: */


%token <string> INTEGER_
%token <string> BC_ 
%token <string> END_
%token <string> CR_
%token <string> FLOAT_
%token <string> STRING_
%token <string> COLON_
%token <string> TABLE_


%type <string>	end_bc_line
%type <floating> table_data_line
%type <string> 	lines
%type <string> 	tag_line
%type <string> 	table_file
%type <floating> numerical_constant
%type <floating> numerical_const_list
%type <string> blank_line
%%

/************************************* BC Table File Rulse: **************************************/

table_file:	lines {} 
		| table_file lines {}
		| error {}
		;

lines:
	/*empty */ {}
	| table_data_line {}	
	| end_bc_line {}
	| tag_line {}
	| blank_line {}
	| error {}
	;
	
tag_line:
	 STRING_ COLON_ CR_ 
	{ 
	  if(!strcmp(tag,$1))
	  {
	    floating_point_constant_list_index = 0;
	    accept_table_data_cards = TRUE;
	    tag_needed = 0;
	    if (ProcID == 0) fprintf(parser_log,"%s:\n   FILE->  ",$1);
	  }	    
	}
	| error CR_ {}	
	;

table_data_line:
	numerical_const_list CR_
	{ 
          if (accept_table_data_cards )
          { 
            if (columns_expected == floating_point_constant_list_index)
            {		
              /*line_read(cards[GD_TABLE_DATA_CARD], 0);*/
              /*table_array[number_of_table_cards_read][0] = atof($1);*/
              /*line_read(cards[GD_TABLE_DATA_CARD], 0);*/
              for(k=0;k<floating_point_constant_list_index;k++)
              {
                table_array[number_of_table_cards_read][k] = floating_point_constant_list_array[k];
              }
              floating_point_constant_list_index = 0;
              number_of_table_cards_read++;
            }
            else
            {	
              sprintf(msg,"   ^^^<<< ERROR %i: %i columns expected, %i columns found. >>>\n   FILE->  ",error_number++,columns_expected, floating_point_constant_list_index);
              fprintf(parser_log,"%s", msg);
              floating_point_constant_list_index = 0;
            }
          }
	}
	| error {declare_error("Invalid GD_TABLE data card."); floating_point_constant_list_index = 0;}
	;

end_bc_line:
	END_ TABLE_ CR_
	{ 	   
	    accept_table_data_cards = FALSE;
	}
	| error CR_ {}
	;
	

numerical_const_list:
		numerical_constant {floating_point_constant_list_array[floating_point_constant_list_index] = $1; floating_point_constant_list_index++;}
		| numerical_const_list numerical_constant {floating_point_constant_list_array[floating_point_constant_list_index] = $2; floating_point_constant_list_index++;}
		| error {}
		;
		
numerical_constant:	FLOAT_ {$$ = atof($1);}
		| INTEGER_ {$$ = (float)(atoi($1));}
		| error {}
		;

blank_line: CR_ {}
	;

%%

int parse_bc_table_file
(
char * tablefile,
char * tabletag,
int    columns
)
{ 
  reading_from_table_file = 1;
  columns_expected = columns;
  strcpy(tag,tabletag);
  if(!strcmp(tag,"NA"))
  {
    accept_table_data_cards = TRUE;
  } 
  else
  {
    accept_table_data_cards = FALSE;
    tag_needed = 1;
  }
  if (ProcID == 0) 
  {
    bcstring=(char *)malloc(10000*sizeof(char));
    if(!bcstring)
    {
      declare_error("Could not allocate space for bc data string");
      reading_from_table_file = 0;
      return(-1);
    }
    if((bc_file = fopen_aprepro( tablefile, "r")) == NULL)
    {
      sprintf(msg,"<<< ERROR %i: Open of %s file failed.>>>\n",error_number++,tablefile);
      fprintf( parser_log,"%s",msg);
      reading_from_table_file = 0;
      return(1);
    } else {
      bc_file_to_string(bc_file, bcstring, 10000);
      fprintf(parser_log,"\n\n   FILE->  BEGIN READING %s\n   FILE->  ",tablefile);
      tparse();
      accept_table_data_cards = FALSE;
      if (tag_needed)
      {
        sprintf(msg,"   <<< ERROR %i: %s name tag expected but not found. >>>\n   FILE->  ",error_number++,tag);
        fprintf(parser_log,"%s", msg);
      }
      fprintf(parser_log,"END READING %s\n\n",tablefile);      
      /*fprintf(parser_log,"\n%i ", ++line_number);*/
      reading_from_table_file = 0;
      return(0);
    } 
    reading_from_table_file = 0;
    return(0);
  }  
  else 
  {
    /* this is not proc 0.  Do whatever it takes to get the string */
    /* then parse the string */
    tparse();
    return(0);  
  }
}

int twrap(const char *msg)
{ 
  return 1; 
}

int terror()
{ 
  yyclearin;
  if (ProcID == 0 && accept_table_data_cards) fprintf(parser_log, " <---<<< ERROR %d>>> ", error_number++);
  return 0;  
}

int bc_file_to_string 
  (
  FILE *file, char string[], 
  int max_number_chars_in
  )
{	 
  int i = 0;
  int new_ch;
  trestart();
  bzero(bcstring, 10000*sizeof(char));
  while ( (i < max_number_chars_in) && (new_ch = getc(file))
          && (new_ch != EOF) ) 
    {
    /*if (new_ch != '\n')*/ 
  	string[i++] = new_ch;
    /*else
  	string[i++] = ' ';*/
    }
  if (i == max_number_chars_in - 1)
    { 
    fprintf(stderr, "file_to_string:%d The maximum number of characters to be read is exceeded.\n", 
    max_number_chars_in);
    return(-1);
    } 
  string[i] = '\0';
  if(ProcID == 0) { /* Code to mpi broadcast the string goes here.  To be added later by geniuses; will work first time. */ }
  return (i+1);
}

