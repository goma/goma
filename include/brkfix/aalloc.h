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


/* aalloc.h - function declarations for array and genl memory allocation 
 */

#ifndef _AALLOC_H
#define _AALLOC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _AALLOC_C
#define EXTERN /* do nothing */
#endif

#ifndef _AALLOC_C
#define EXTERN extern
#endif

#if 0
typedef unsigned long ULONG;

EXTERN void *aalloc 
PROTO((int numdim,
	...));			/* variable number of arguments */

#endif

#define smalloc(arg)	safe_malloc((arg), __FILE__, __LINE__)

EXTERN void *safe_malloc 
PROTO((const int ,		/* numbytes                                  */
       const char *,		/* filename                                  */
       const int ));		/* line                                      */

EXTERN void safe_free
PROTO((void *));		/* ptr to previously allocated memory */

#endif /* _AALLOC_H */
