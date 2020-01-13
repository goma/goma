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

#ifdef hpux
/*
 * This file handles ieee floating point exceptions on the HP.  Files for
 * other computers should be similar to this one.
 */

#ifndef _HPUX_SOURCE
#define _HPUX_SOURCE
#endif

#endif

/*
 * To explicitly control ieee_handler() provided by Sun...
 */

#ifdef HAVE_SUNMATH_H
#include <sunmath.h>		/* Need "-lsunmath -lm" libs too! */
#endif

void handle_ieee (void );

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef hpux
void handle_ieee (void)

/***********************************************************************
 *    This function changes the default ieee floating point exceptions
 *    so that common ieee errors will cause an immediate exit
 *    fpsetdefaults() changes the default environment on Series 700
 *    workstations, which is
 *
 *           Round to nearest (FP_RN)
 *           All exception flags cleared (FP_X_CLEAR)
 *           All exception traps disabled
 *           Fast underflow mode disabled
 *
 *      fpsetdefaults() changes these defaults to more useful values.
 *      Specifically, it enables traps for the invalid operation, divide-by-
 *      zero, and overflow exceptions, while leaving the underflow and
 *      inexact-result exception traps disabled.  It sets the environment as
 *      follows:
 *
 *            Round to nearest (FP_RN)
 *            All exception flags cleared (FP_X_CLEAR)
 *                    (FP_X_INV+FP_X_DZ+FP_X_OFL)
 *            All exception traps enabled except underflow and inexact result
 *                    (FP_X_INV+FP_X_DZ+FP_X_OFL)
 *            Fast underflow mode enabled (if the system supports it)
 *
 * NOTE -> HKM 2/29/00:
 *         hpux' cc compiler has a major hassle associated with it.
 *         If floating point exception traps are enabled and if the
 *         default optimization is turned on, then spurious SIGFPE 8
 *         signals may occur, causing your program to prematurely
 *         terminate.
 *         There are two solutions to this (note already called HP on 
 *         this and was told that this is a undocumented "feature" not
 *         (gasp) a bug):
 *           1) Change the compile line options to include +Onomoveflops
 *              optimization option
 *           2) Turn off exception traps when optimization is used.
 *         The code below implements choice #2, keying off of the DEBUG
 *         compiler define.
 ***********************************************************************/
{
#ifdef DEBUG
   (int) fpsetfastmode(0);  
   (void) printf("ieee: fp exceptions are ignored\n");
   (void) printf("ieee: Fast underflow mode is not enabled. IEEE-754 arithmetic\n");
#else
   (void) fpsetfastmode(1);
   (void) printf("ieee: fp exceptions are ignored\n");
   (void) printf("ieee: Fast underflow mode is enabled\n");
#endif
}
/*****************************************************************************/

#endif

#ifdef solaris
void 
handle_ieee(void )
{
  static int err;
  static const char yo[] = "handle_ieee";

  err = 0;

  /*
   * Clear any previous existing settings...
   */

  err = ieee_handler ("clear", "all", SIGFPE_DEFAULT);

  /*
   * Call these harmless exceptions...
   */

  err = ieee_handler ("set", "inexact",   SIGFPE_IGNORE);
  err = ieee_handler ("set", "underflow", SIGFPE_IGNORE);

  /*
   * Call these exceptions worthy of halting right away...
   */

  err = ieee_handler ("set", "division",  SIGFPE_ABORT);
  err = ieee_handler ("set", "overflow",  SIGFPE_ABORT);
  err = ieee_handler ("set", "invalid",   SIGFPE_ABORT);

  /*
   * Insure none of these set attempts was thwarted...
   */

  if ( err != 0 )
    {
      log_msg("Trouble initializing ANSI/IEEE Std 754-1985");
      log_err("arithmetic handler on solaris");
    }

  return;
}

#if 0
/*
 * Just in case, a deprecated SVID3 implementation of matherr()
 */

int 
matherr(struct exception *e)
{
  if ( e->type == OVERFLOW )
    {
      log_err("Overflow from %s with args %g %g", e->name, e->arg1, e->arg2);
      abort();
    }
  else if ( e->type == UNDERFLOW )
    {
      log_msg("Underflow from %s with args %g %g", e->name, e->arg1, e->arg2);
      e->retval = 0;
    }
  return(0);
}
#endif

#endif

#if !defined(hpux) && !defined(solaris)

/*
 * Default for this: do nothing.
 */

void 
handle_ieee (void)
{
#ifdef DEBUG
  (void) printf("ieee: fp exceptions are not changed: generic block\n");
#endif
}
#endif
