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

#ifndef GOMA_MM_EH_H
#define GOMA_MM_EH_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_EH_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_EH_C
#define EXTERN extern
#endif

typedef int goma_error;
#define GOMA_ERROR   -1
#define GOMA_SUCCESS 0
#define GOMA_DEBUG   1

#include <stdarg.h> /* for var args... */

#include "std.h"

/* Needed to use mm_eh without std.h */
#ifndef MAX_CHAR_ERR_MSG
#define MAX_CHAR_ERR_MSG 1024
#endif

/*
 *  This variable is included as a utility to writing informative
 *  strings before error exiting.
 */
extern char Err_Msg[MAX_CHAR_ERR_MSG];

EXTERN void goma_eh(const int error_flag, const char *file, const int line, const char *message);

EXTERN void
goma_wh(const int error_flag, const char *const file, const int line, const char *const message);

EXTERN void save_place  /* mm_eh.c                                   */
    (const int,         /* severity                                  */
     const char *const, /* routine_name                              */
     const char *const, /* file_name                                 */
     const int);        /* line_number                               */

EXTERN void logprintf /* mm_eh.c                                   */
    (const char *,    /* format                                    */
     ...);            /* var args */

EXTERN void smooth_stop_with_msg(const char *msg);

/* This macro expands to a function call to the error handler eh() with
 * four arguments:
 *		    IERR -- an int, if IERR=-1 there is a problem
 *		    MESSAGE -- a string, which gets printed if the
 *			       error dump is activated by IERR=-1
 *
 *	goma_eh() finally exits with -1 if an error occurred.
 *
 * Intent:
 * -------
 *		EH(return_code, "I am informative.");
 */

#define EH(IERR, MESSAGE) goma_eh(IERR, __FILE__, __LINE__, MESSAGE)
#define WH(IERR, MESSAGE) goma_wh(IERR, __FILE__, __LINE__, MESSAGE)

#define TIME_STRING_SIZE (256)

extern char current_routine[]; /* name of current routine. */
extern char current_file[];    /* name of current file */
extern int current_line;       /* line number in file */
extern int current_severity;   /* global error signal (-1=die,0=prnt,1=dbg) */

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
#define log_err(...)                                \
  do {                                              \
    save_place(GOMA_ERROR, yo, __FILE__, __LINE__); \
    logprintf(__VA_ARGS__);                         \
  } while (0)

#define log_msg(...)                                  \
  do {                                                \
    save_place(GOMA_SUCCESS, yo, __FILE__, __LINE__); \
    logprintf(__VA_ARGS__);                           \
  } while (0)

#define log_dbg(...)                                \
  do {                                              \
    save_place(GOMA_DEBUG, yo, __FILE__, __LINE__); \
    logprintf(__VA_ARGS__);                         \
  } while (0)

#ifndef DEFAULT_GOMA_LOG_FILENAME
#define DEFAULT_GOMA_LOG_FILENAME ".log"
#endif

#endif
