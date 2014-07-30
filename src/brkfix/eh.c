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

/* eh -- error handling routine that puts out filename and linenumber*/

#include <stdio.h>

#define _EH_C

#include "std.h"
#include "eh.h"
#include <stdlib.h>

void eh(const int error_flag,
	const char *file,
	const int line,
	const char *message)
{
  if ( error_flag == -1 )
    {
      fprintf(stderr, "%s:%d: %s\n", file, line, message);
      exit(-1);
    }
  else
    {
      return;
    }
}
