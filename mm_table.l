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

%{


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#include "t.tab.h"
#include "std.h"
#include "rf_io_const.h"
#include "mm_input.h"
#include "mm_parser.h"

#undef YY_INPUT
#define YY_INPUT(bcbuf,result,max) (result = bcinput(bcbuf,max))


/*jjj#undef yyinput*/
/*jjj#undef yyunput*/

/*#define input() (*mystring++)jjj*/
/*#define unput(c) (*--mystring = c)jjj*/

/*jjj#define yyinput(void) (return input(); )*/
/*jjj#define yyunput(int c) (unput(c);)jjj*/

extern char *bcstring;
extern parser_file_mode;
/*extern int line_number  first_tab_flag*/ 
extern FILE *parser_log;
extern int ProcID;
extern int print_strings_to_log;
extern int count_of_characters_in_this_line;
extern int error_found_in_last_line;
extern int error_number;
extern int line_number;
extern int first_tab_flag;
extern int accept_table_data_cards;
int first_cr = 1;





void bc_print_to_log();

%}
%e 100000
%p 200000
%a 100000
%n 40000
%k 2000

%%




"BC"	 				{bc_print_to_log(); strcpy(yylval.string, yytext);   return(BC_);}  
"END" 					{bc_print_to_log(); strcpy(yylval.string, yytext);   return(END_);}  
"TABLE"					{bc_print_to_log(); strcpy(yylval.string, yytext);   return(TABLE_);}
":"					{bc_print_to_log(); strcpy(yylval.string, yytext);   return(COLON_);}


[+-]?[0-9]+				{bc_print_to_log();  strcpy(yylval.string, yytext); return (INTEGER_);}

[+-]?(([0-9]+)|([0-9]*\.[0-9]*)([eE][-+]?[0-9]+)?)	{bc_print_to_log();  strcpy(yylval.string, yytext);   return (FLOAT_);}



^\n					{
					bc_print_to_log();		/* blank line */ 
					}


\s\n					{
					bc_print_to_log();		/* blank line */ 
					}



^\n*\n					{
					bc_print_to_log();		/* blank line */ 
					}
					

^(" "|\t)(" "|\t)*\n			{
					bc_print_to_log();		/* another blank line */
					}




				

[-|*|$|%|#|/].*\n       		{
					bc_print_to_log(); 
					if(ProcID == 0 && accept_table_data_cards) fprintf(parser_log,"   FILE->  ");
					}	/* commnet line */
			

^(" "|\t)*[-|*|$|%|#|/].*\n		{
					bc_print_to_log(); 
					if(ProcID == 0 && accept_table_data_cards) fprintf(parser_log,"   FILE->  ");
					}	/* commnet line */
					
[a-zA-Z_][0-9a-zA-Z_@\.]*		{bc_print_to_log(); strcpy(yylval.string, yytext); return (STRING_);}
					

\n					{
					/*bc_print_to_log();*/
					strcpy(yylval.string, yytext);
					if(ProcID == 0 && accept_table_data_cards) fprintf(parser_log,"\n   FILE->  ");
					return(CR_);
					
					}
					
									
\t					{
					if (ProcID == 0) fprintf(parser_log, " ");
					count_of_characters_in_this_line++;
					while ( (count_of_characters_in_this_line % 8) > 0)
						{ 
						if (ProcID == 0) fprintf(parser_log," ");
						count_of_characters_in_this_line++;
						}
					}
   
.					{bc_print_to_log();}

%%

/*
<COMMENT>([a-zA-Z_]|.)?([0-9a-zA-Z_@\.]|.)?\n				{printf("%s", yytext);   BEGIN 0;}
*/

bcinput(char *bcbuf, int max)
{  
 
  if(bcstring && strlen(bcstring))
    {

      *bcbuf = *bcstring++;

      return 1;
	/*
      if (len < strlen(bcstring))
	{	
	  memcpy(bcbuf, bcstring++, sizeof(char));
	  len++;
	  return 1;
	}
      else 
	return 0;
	*/
    }
  else
    {
      fflush(stdout);
     /* fprintf(stdout,"Uh oh! I was passed a null string!\n");*/
      fflush(stdout);
      return 0;
    }
}

void bc_print_to_log()
{
   stall(10000000);
        if (ProcID == 0 && accept_table_data_cards) 
 	{ fflush(parser_log);
  	  fprintf(parser_log, "%s", yytext);
 	}


}

int stall (int count)
{
  int s;
  for(s=0; s<count; s++) {/* stall */ s++;}
  return(0);
}


