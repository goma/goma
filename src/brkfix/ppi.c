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

/* pre_process() -- filter the input file into a tmp file; return tmp filename
 *
 *
 * Allow for nice comments, etc. in the input file. Pre-process
 * as follows:
 *		(i)   remove comments (everything between # and newline)
 *		(ii)  delete blank lines
 *		(iii) remove trailing whitespace
 *		(iv)  remove leading whitespace
 *
 *
 * Created: 1997/05/08 08:49 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#define _PPI_C

#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

#include "map_names.h"
#include "std.h"
#include "eh.h"
#include "ppi.h"

void
pre_process(char *fn)
{

  FILE *infile;
  FILE *temp_file;

  char temp[] = TEMP_PREFIX;
  char blank_string[] = " ";
  char * nfn;

  char buffer[80];
  char in_buffer[80];

  int fn_len;
  int prefix_len;

  int i;
  int j;
  int k;
  char * str1;
  char * str2;

  int err1=0;
  int err2=0;

  char * tmp;

  if ( ( infile = fopen( fn, "r" ) ) == NULL )
    {
      sprintf( buffer, "Error opening input file: %s\n", fn );
      fprintf(stdout, buffer );
      fflush( stdout );
      exit(1);
    }

  prefix_len = strlen( temp );
  fn_len = strlen( fn );
  nfn = malloc( ( fn_len + prefix_len +3 ) * sizeof( char ) );
  nfn[0] = '\0';

  /*
   * Temporary file name will be "tmp.foo"...
   */
  strcat(nfn, temp);
  strcat(nfn, fn);
  if ((temp_file = fopen( nfn, "w" ) ) == NULL)
    {
      sprintf( buffer,
               "Error opening temporary input file: %s\n", nfn );
      fprintf(stdout,buffer);
      fflush( stdout );
      exit(1);
    }

  str1 = malloc( 2 * sizeof( char ) );
  str1[1] = '\0';
  str2 = malloc( 2 * sizeof( char ) );
  str2[1] = '\0';

  for( k=0; k< 80; ++k ) buffer[k] = '\0';

  fgets( in_buffer, 80, infile );
  while( feof( infile ) == 0 )
    {
      str1[0] = in_buffer[0];
      while( strcmp( str1, blank_string ) == 0 ||
             strcmp( str1, "\t" ) == 0 ) {  /* strip leading blanks and tabs */
         for( i=0; i < ( strlen( in_buffer ) - 1 ) ; i++ ) {
             in_buffer[i] = in_buffer[i+1];
             str1[0] = in_buffer[i];
             str2[0] = in_buffer[i+1];
             if( strcmp( str1, "\n" ) == 0 && strcmp( str2, "\n") == 0 )
                 in_buffer[i+1] = '\0';   /* guard against double newlines */
         }
         str1[0] = in_buffer[0];
      }
    /* process the line if it is not a blank line and not a comment line */
      if( strcmp( str1, "\n" ) != 0 && strcmp ( str1, "#" ) != 0 ) {
        for( j=0; j < strlen( in_buffer ); j++ ) {
          str1[0] = in_buffer[j];
          if( strcmp( str1, "#" ) != 0 ) {    /* not a comment line */
            if( strcmp( str1, "\n" ) != 0 ) { /* not a newline */
              if( strcmp( str1, "\t" ) != 0 ) { /* not a tab */
                buffer[j] = in_buffer[j];
              }
              else                           /* replace tab with space */
                buffer[j] = ' ';
            } else {                         /* terminate with newline */
            buffer[j] = '\n';
            buffer[j+1] = '\0';
            }
          } else {                           /* comment - end the line */
            buffer[j] = '\n';
            buffer[j+1] = '\0';
            break;
          }
        }
        i = strlen( buffer ) - 2 ;  /* extra offset for the newline */
        str1[0] = buffer[ i ];
        while( strcmp( str1, blank_string ) == 0 ) { /* strip trailing blanks */
           buffer[i] = '\n';
           buffer[i+1] = '\0';
           str1[0] = buffer[--i];
        }
        str1[0] = buffer[0];
        if( strcmp( str1,"\n" ) != 0 ) {     /* save the line to temp file */
          fprintf( temp_file, buffer );
          for( k=0; k< 80; ++k ) buffer[k] = '\0';
        }
      }   /* end blank line check  and  comment line check */
      fgets( in_buffer, 80, infile );
    }   /* end of eof while */

  err1 = fclose( infile );

  fflush( temp_file );
  err2 = fclose( temp_file );

  if( err1 != 0 || err2 != 0 ) {
    fprintf(stdout," error closing files\n");
    fflush( stdout );
  }

  tmp = strcpy(fn, nfn);

  if ( tmp == NULL ) {
    fprintf( stdout," error creating temporary input file name\n");
    fflush( stdout );
    exit(2);
  }

  free( str1 );
  free( str2 );
  free( nfn );

  return;
}
