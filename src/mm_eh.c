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
 *$Id: mm_eh.c,v 5.3 2008-01-11 00:47:14 hkmoffa Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: mm_eh.c,v 5.3 2008-01-11 00:47:14 hkmoffa Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "std.h"
#include "rf_mp.h"
#include "mpi.h"
#include "rf_io_const.h"

#define GOMA_MM_EH_C
#include "mm_eh.h"

#ifdef PRINT_STACK_TRACE_ON_EH
#include <execinfo.h>

/* Obtain a backtrace and print it to stdout. */
/* https://www.gnu.org/software/libc/manual/html_node/Backtraces.html */
/* compile with -rdynamic */
static void
print_stacktrace(void)
{
  void *bt_array[20];
  size_t size;
  char **stack_strings;

  size = backtrace(bt_array, 20);
  stack_strings = backtrace_symbols(bt_array, size);

  printf("Obtained %zd stack frames from Proc %d.\n", size, ProcID);

  for (size_t i = 0; i < size; i++) {
     printf("%s\n", stack_strings[i]);
  }

  free (stack_strings);
}
#endif

/*
 * These variables are extern and included everywhere else, defined and
 * initialized here. They are important to proper switching during error
 * handling.
 */

char current_routine[MAX_FNL] = "?";
char current_file[MAX_FNL]    = "?";
int  current_line             = -1; 
int  current_severity         = 0;

/*
 *  This variable is included as a utility to writing informative
 *  strings before error exiting. Including it as a global variable
 *  means that you don't need to include a similar variable in
 *  each file or routine that needs error handling.
 */
char Err_Msg[MAX_CHAR_ERR_MSG];

/*
 * log - stuff to print out routinely that is generally of interest.
 *       sequence and time is important, so a /var/adm/log format
 *       is used.
 */
#ifdef ENABLE_LOGGING
static FILE *log_strm=NULL;
static char log_filename[MAX_FNL] = DEFAULT_GOMA_LOG_FILENAME;
#endif


/*
 * Global variable helps to put the brakes on a parallel train wreck, but
 * only if you take time to broadcast it...
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
void 
eh(const int error_flag, const char *file,
   const int line, const char *message)
{
  static char yo[] = "eh";
  if (error_flag == -1) { 
     log_msg("GOMA ends with an error condition.");
#ifdef PRINT_STACK_TRACE_ON_EH
    print_stacktrace();
#endif
#ifndef PARALLEL
    fprintf(stderr,"ERROR EXIT: %s:%d: %s\n", file, line, message); 
    exit(-1);
#else
    fprintf(stderr, "P_%d ERROR EXIT:%s:%d: %s\n", ProcID, file, line, message);
    if (Num_Proc == 1) {
      MPI_Finalize();
      exit(-1);
    }
    parallel_err = TRUE;
#endif
  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
wh(const int error_flag, const char * const file, const int line,
   const char * const message, int * const iw)
{
  if ( error_flag == -1 )
    {
      DPRINTF(stderr, "WARNING:  %s:%d: %s\n", file, line, message);
      *iw = 1;
      return;
    }
  else
    {
      return;
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void
aborth(const int error_flag, const char *  const file, const int line,
       const char * const message)

    /************************************************************************
     *  Wrapper around the command MPI_Abort: This is used when a
     *  graceful exit from MPI isn't possible -> for example from the
     *  input routines.
     ************************************************************************/
{
  fprintf(stderr, "P_%d ABORT:%s:%d: %s\n", ProcID, file, line, message);
  fflush(stderr);
#ifdef PARALLEL
  MPI_Abort(MPI_COMM_WORLD, error_flag);
#endif
  exit(error_flag);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* save_place() -- record the current filename, linenumber, routinename
 *
 * Notes: This saves local information about where a problem is occurring
 *        so that it may accessed from error handler routines by inspecting
 *        the corresponding global variables.
 *
 *	  If no local static routinename has been defined then references
 *        revert to a global "unassigned".
 *
 * Created: 1999/09/15 13:01 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void 
save_place(const int severity,
	   const char * const routine_name, 
           const char * const file_name, 
           const int line_number)
{
  strcpy(current_routine, routine_name);
  strcpy(current_file, file_name);
  current_line     = line_number;
  current_severity = severity;
  return;
}

void 
logprintf(const char *format, ... )
{
#ifndef ENABLE_LOGGING
  if ( current_severity < 0 )
    {
      va_list va;
      va_start(va, format);
      vfprintf(stderr, format, va);
      va_end(va);
      DPRINTF(stderr, "\nAbnormal termination in logprintf\n");
      exit(current_severity);
    }
  return;
#else
  time_t now;
  static char time_format[] = "%b %d %T";
  static char time_result[TIME_STRING_SIZE];
  int n;

  static char new_format[1024];

  static char old_buffer[1024];	/* For legacy help... */

  va_list va;

  /* Proceed only if there are no print limits set */

  if( !Unlimited_Output ) return;

  /*
   * For now direct all log entries to one file. In parallel you may want
   * to multiplex into separate files for each processor. Finally, you might
   * want to separate heavy debugging messages.
   */


  /*
   * Check for an open file output stream, open one if there isn't one yet...
   */

  if ( log_strm == NULL )
    {
      if ( Num_Proc > 1 )
	{
	  multiname(log_filename, ProcID, Num_Proc);
	}

      log_strm = fopen(log_filename, "a");
      
      if ( log_strm == NULL )
	{
	  current_severity=-1;
	  exit(current_severity);
	}
    }

  time(&now);

  (void) strftime(time_result, TIME_STRING_SIZE, 
                  time_format, localtime(&now));

  /*
   * Some niceties: if the user appended a carriage return, fine. If not,
   * then we'll add one.
   */

  n = strlen(format);

  if ( n > 1024 )
    {
      exit(0);
    }

  strcpy(new_format, format);

  if ( format[n] != '\n' )
    {
      strcat(new_format, "\n");
    }

  fprintf(log_strm, "%s ", time_result);

  if ( current_severity < 0 )
    {
      if ( Num_Proc > 1 )
	{
	  fprintf(log_strm, "P_%d ", ProcID);
	}
      fprintf(log_strm, "error at %s:%d %s() \n", current_file, current_line, 
	      current_routine);

      fprintf(log_strm, "%s ", time_result);
    }

  if ( Num_Proc > 1 )
    {
      if ( current_severity == 0 )
	{
	  fprintf(log_strm, "P_%d ", ProcID);
	}
    }

  if ( current_severity > 0 )
    {
      fprintf(log_strm, "%s:%d %s() ", current_file, current_line, 
	      current_routine);
      if ( Num_Proc > 1 )
	{
	  fprintf(log_strm, "P_%d ", ProcID);
	}
    }

  va_start(va, format);
  vfprintf(log_strm, new_format, va);
  va_end(va); 

  /*
   * This flushing is great for precise error determination, but could
   * cause big time performance loss on many processors, for example.
   */

  if ( Num_Proc < 16 )
    {
      fflush(log_strm);
    }

  /*
   * Need to pipe out to stderr for legacy usage...delete this once
   * user trust the log file.
   */

  va_start(va, format);
  vsprintf(old_buffer, new_format, va);
  va_end(va);



  if ( current_severity < 0 )
    {
      fprintf(log_strm, "%s ", time_result);
      fprintf(log_strm, "GOMA terminates abnormally.\n");


#ifdef PARALLEL
      DPRINTF(stderr, "\nAbnormal termination -- see log file for details.\n");
#endif
      exit(current_severity);
    }

  return;
#endif
}


void
smooth_stop_with_msg( const char *msg )
{

#ifndef PARALLEL
	fprintf(stderr, "\n\t -- Goma stops smoothly-- \n");
	fprintf(stderr, msg );
	exit(-1);
#else
	if( Num_Proc == 1 )
	{
		fprintf(stderr, "\n\t -- Goma stops smoothly--\n");
		fputs(msg, stderr);
	}
	else
	{
		fprintf(stderr,"\n\n Proc %d -- Goma stops smoothly --\n", ProcID);
		DFPUTS(msg, stderr);
	}	
	MPI_Finalize();
	exit(-1);
#endif
}
	
