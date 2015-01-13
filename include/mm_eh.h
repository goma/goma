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
 *$Id: mm_eh.h,v 5.2 2007-09-18 18:53:42 prschun Exp $
 */

#ifndef _MM_EH_H
#define _MM_EH_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_EH_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_EH_C
#define EXTERN extern
#endif

#include <stdarg.h>		/* for var args... */

/* Needed to use mm_eh without std.h */
#ifndef MAX_CHAR_ERR_MSG
  #define MAX_CHAR_ERR_MSG 1024
#endif

/*
 *  This variable is included as a utility to writing informative
 *  strings before error exiting.
 */
extern char Err_Msg[MAX_CHAR_ERR_MSG];

EXTERN void eh
PROTO((const int ,	       	/* error_flag                                */
       const char * const,      /* file                                      */
       const int ,	       	/* line                                      */
       const char * const));	/* message                                   */
      
EXTERN void wh
PROTO((const int ,	       	/* error_flag                                */
       const char * const,      /* file                                      */
       const int ,	       	/* line                                      */
       const char * const,      /* message                                   */
       int * const));		/* iw                                        */
      
extern void aborth(const int, const char * const, const int,
		   const char * const);

EXTERN void save_place		/* mm_eh.c                                   */
PROTO((const int ,		/* severity                                  */
       const char * const,	/* routine_name                              */
       const char * const,	/* file_name                                 */
       const int ));		/* line_number                               */

EXTERN void logprintf		/* mm_eh.c                                   */
PROTO((const char *,		/* format                                    */
       ... ));			/* var args */
	   
	   
EXTERN void smooth_stop_with_msg
PROTO((const char *  msg));



/* This macro expands to a function call to the error handler eh() with
 * four arguments:
 *		    IERR -- an int, if IERR=-1 there is a problem
 *		    MESSAGE -- a string, which gets printed if the
 *			       error dump is activated by IERR=-1
 *
 *	eh() finally exits with -1 if an error occurred.
 *
 * Intent:
 * -------
 *		EH(return_code, "I am informative.");
 */

#define     EH(IERR, MESSAGE)	eh(IERR, __FILE__, __LINE__, MESSAGE)
#define	    WH(IERR, MESSAGE)	{static int iw=0; if (iw == 0) wh(IERR, __FILE__, __LINE__, MESSAGE, &iw);}
#define     ABORTH(IERR, MESSAGE) aborth(IERR, __FILE__, __LINE__, MESSAGE)

#define GOMA_ERR	(-1)	/* definitely print; definitely exit */
#define GOMA_MSG	(0)	/* usually print; continue */
#define GOMA_DBG	(1)	/* maybe print; continue */

#define TIME_STRING_SIZE	(256)

extern char current_routine[];	/* name of current routine. */
extern char current_file[];	/* name of current file */
extern int current_line;	/* line number in file */
extern int current_severity;	/* global error signal (-1=die,0=prnt,1=dbg) */

/*
 * Yes, this is an accident waiting to happen to the unwary that don't
 * use curly braces for one liners like
 *
 *	if ( condition ) log_err("blah");
 *
 * that will break unless they are written as
 *
 *      if ( condition ) { log_err("blah"); }
 *
 * If you can figure out a way to do it that permits varargs and throws in
 * the extra args so the user doesn't have to do it by hand, then have at it.
 *
 * Requires use of a local static const char yo[] ="Routine_name".
 */

#define log_err	save_place(GOMA_ERR, yo, __FILE__, __LINE__); logprintf
#define log_msg save_place(GOMA_MSG, yo, __FILE__, __LINE__); logprintf
#define log_dbg save_place(GOMA_DBG, yo, __FILE__, __LINE__); logprintf

#ifndef DEFAULT_GOMA_LOG_FILENAME
#define DEFAULT_GOMA_LOG_FILENAME ".log"
#endif

#endif
