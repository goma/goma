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

/* ppi.h - prototype declarations for ppi.c
 */

#ifndef _PPI_H
#define _PPI_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _PPI_C
#define EXTERN /* do nothing */
#endif

#ifndef _PPI_C
#define EXTERN extern
#endif

#ifdef _PPI_C
const char filter[]="sed -e 's/#.*$//' -e '/^[ 	]*$/d' -e 's/[ 	]*$//'";
#define TEMP_PREFIX	"tmp."
#endif

EXTERN void pre_process
PROTO((char *));		/* input filename */

#endif /* _PPI_H */


