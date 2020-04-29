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
 *$Id: md_timer.c,v 5.1 2009-02-26 23:28:24 hkmoffa Exp $
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef tflop
#ifndef hpux
#ifndef _AIX
#ifndef sun4
#include <sys/times.h>
#endif
#endif
#endif
#endif

#ifdef hpux
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
#endif

#ifdef _AIX
#include <sys/times.h>
#include <time.h>
#endif

#ifdef OSF1
#include <sys/limits.h>
#include <sys/time.h>
#endif

#ifdef dec_osf1
#include <time.h>
#endif

#ifdef sun4
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
#endif

#ifdef linux
#include <time.h>
#endif

#ifdef tflop
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
#endif

#define GOMA_MD_TIMER_C
#include "std.h"
#include "md_timer.h"

/*
 * ut -- return user time in seconds (double).
 */
dbl
ut(void)
{
  clock_t err;
  dbl u;
  static struct tms *tm = NULL;
  static long  tickPerSecond;
  static dbl secondsPerTick;
  if (!tm)
    {
      tm = (struct tms *) calloc(1, sizeof(struct tms));
      tickPerSecond = sysconf(_SC_CLK_TCK);
      secondsPerTick = 1.0 / ((double) tickPerSecond);
    }
  if ((err = times(tm)) == -1)
    {
      fprintf(stderr, "error from times().\n");
      exit(-1);
    }
  u  = ((dbl) (tm->tms_utime)) * secondsPerTick; 
  return(u);
} /* END of routine ut */
/*****************************************************************************/

/*
 * ust -- return sum of user and system time in seconds (double).
 */
dbl
ust(void)
{
  int err;
  dbl u;
  static long  tickPerSecond;
  static dbl secondsPerTick;
  static struct tms *tm = NULL;
  if (!tm)
    {
      tm = (struct tms *) calloc(1, sizeof(struct tms));
      tickPerSecond = sysconf(_SC_CLK_TCK);
      secondsPerTick = 1.0 / ((double) tickPerSecond);
    }

  if ((err = times(tm)) == -1)
    {
      fprintf(stderr, "error from times().\n");
      exit(-1);
    }
  u  = ((double)(tm->tms_utime + tm->tms_stime)) * secondsPerTick;
  return(u);
} /* END of routine ust */
/*****************************************************************************/

/*
 * get_date() -- fill a pre-allocated string with the date as "mm/dd/yy" 
 */
void 
get_date(char *string)
{
  char date_format[80];
  int i;
  char buffer[80];
  time_t now;

  int month_number;
  int day_of_month;
  int year_number;

  for ( i=0; i<80; i++ ) buffer[i] = '\0';

  strcpy(date_format, "%m/%d/%y");
  time(&now);

  /*
   * Yes, I know that %y only returns two digits for the year and that's not
   * in the best interests of Y2K compliance, but EXODUS II species date fields
   * that are only this format. Since we don't do any date arithmetic based
   * on years, this should not cause problems.
   *
   * I've become annoyed at the warnings from this and have implemented
   * a certified warning-free Y2K-complacent version.
   */

  strftime(buffer, 80, "%m", localtime(&now) );
  month_number = atoi(buffer);

  strftime(buffer, 80, "%d", localtime(&now) );
  day_of_month = atoi(buffer);

  strftime(buffer, 80, "%Y", localtime(&now) );
  year_number  = atoi(buffer);

  sprintf(string, "%2d/%2d/%2d", month_number, day_of_month, year_number%100);

  return;
} /* END of get_date() */

/******************************************************************************/

/*
 * get_time() -- fill a pre-allocated string with the time as "hh:mm:ss"
 */

void 
get_time(char *string)
{
  char time_format[80];
  size_t bufsize;
  time_t now;
  
  strcpy(time_format, "%H:%M:%S");
  bufsize = strlen(time_format) + 1;
  time(&now);
  strftime(string, bufsize, "%H:%M:%S", localtime(&now) );
  return;
} /* end of get_time() */

/*****************************************************************************/

/* END of file md_timer.c */
/*****************************************************************************/
