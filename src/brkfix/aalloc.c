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

/* aalloc -- multi-dimensional array allocator 
 *
 * The multidimensional array allocator is from Ray Tuminaro. Occassionally,
 * I find it's behavior suspicious, so I generally use smalloc().
 *
 */

#define _AALLOC_C
#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "map_names.h"
#include "std.h"
#include "aalloc.h"
#include "eh.h"

#if 0


typedef struct 
{
  size_t index;	/* Number of elements in the dimension	*/
  size_t total;	/* Total number of elements 		*/
  size_t size;	/* Size of a single element in bytes	*/
  size_t off;	/* offset from beginning of array	*/
} Adim;

/*
 * Function definitions
 */

#ifdef STDC_HEADERS
void *
aalloc( int n , ... )
#else
void *
aalloc(va_alist)
va_dcl
#endif
{
  int	i;
  int	j;
  Adim *dim;
#ifndef STDC_HEADERS
  int n;		/* Number of dimensions			*/
#endif
  size_t total;		/* Total size of the array		*/
  void *dfield;	/* ptr to avoid lint complaints		*/
  void *field;		/* The multi-dimensional array		*/
  void **ptr;		/* Pointer offset			*/
  void *data;		/* Data offset				*/
  va_list va;		/* Current pointer in the argument list	*/
  
#ifdef STDC_HEADERS
  va_start(va, n);
#else
  va_start(va);
  n = va_arg(va, int);
#endif
  
  if ( n < 1 )
    {
      EH(-1, "Multidimensional array with less than 1 dimension.");
    }

  dim = (Adim *) safe_malloc((int) (n * sizeof(Adim)),  __FILE__, __LINE__);
  
  dim[0].index = va_arg(va, int);

  if ( dim[0].index < 1 )
    {
      EH(-1, "First array dimension < 1");
    }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  
  for (i=1; i<n; i++)
    {
      dim[i].index = va_arg(va, int);
      if ( dim[i].index < 1 )
	{
	  fprintf (stderr, "array dimension %d < 1, %d\n", i+1, dim[i].index);
	}
      dim[i].total = dim[i-1].total * dim[i].index;
      dim[i].size  = sizeof(void *);
      dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
    }
  
   dim[n-1].size = va_arg(va, size_t);
   va_end(va);
   
   total  = dim[n-1].off + dim[n-1].total * dim[n-1].size;
   dfield = (void *) safe_malloc((int) total, __FILE__, __LINE__);

   for ( i=0; i<n-1; i++)
     {
       ptr  = (void *)dfield + dim[i].off;
       data = (void *)dfield + dim[i+1].off;
       for (j=0; j<dim[i].total; j++)
	 {
	   ptr[j] = (void *)data + j * dim[i+1].size * dim[i+1].index;
	 }
     }

   safe_free((void *) dim);

   return(dfield);
}

#endif

/* Safe version of malloc.  Does not initialize memory .*/

/* Modified by Scott Hutchinson (1421) 20 January 1993 */

/*
 * safe_malloc() - minor protective wrapper around malloc()
 *
 * Normally, the C preprocessor macro smalloc() is used, which automatically
 * appends the filename and line number arguments so the user does not need
 * to do this. Then, if DEBUG is defined below, you can watch how memory
 * gets hogged.
 *
 * Revised: 1998/08/05 08:27 MDT pasacki@sandia.gov
 */

void *
safe_malloc(const int n,	/* numbytes requested */
	    const char *filename, /* of caller */
	    const int line)	/* of caller */
{
  void *pntr;

  if ( n < 0 )
    {
      fprintf(stderr, "Negative memory specification (%d) at %s:%d\n",
	      n, filename, line);
      exit(-1);
    }
  else if ( n == 0 )
    {
#ifdef AR
      fprintf(stderr, "Warning - zero memory requested at %s:%d\n",
	      filename, line);
#endif
      pntr = NULL;
    }
  else
    {
      pntr = malloc((size_t) n);
    }

  if ( pntr == NULL && n != 0)
    {
      fprintf(stderr, "Memory allocation failure for %d bytes at %s:%d\n", 
	      n, filename, line);
      exit(-1);
    }

#ifdef DEBUG
  fprintf(stderr, "%s:%d: allocating %d bytes @ %x\n", filename, line, n,
	  pntr);

  /*
   * Verify we can write something to this piece...
   */
  
  if ( n > 0 && pntr != NULL )
    {
      *pntr = 0;
    }

#endif

  return(pntr);
} /* END of routine safe_malloc() */

/* safe_free() -- don't free any NULL pointers
 *
 *
 * Revised: 1998/08/05 08:32 MDT pasacki@sandia.gov
 */

void
safe_free (void *ptr)
{
  if ( NULL == ptr ) return;
#ifdef DEBUG
  fprintf(stderr, "Attempt to free a NULL ptr.\n");
#endif
  free(ptr);
  /*
  if ( ptr != NULL )
    {
      free(ptr);
    }
  return;
  */
}

