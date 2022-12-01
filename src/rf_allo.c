/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 * Some low level support utilities for memory allocation for multidim arrays.
 */

/*
 *$Id: rf_allo.c,v 5.5 2008-03-22 00:55:50 hkmoffa Exp $
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
*/
#include <stdarg.h>

/*
 * Default behavior on out of memory is to print an error message and return
 */
static int ALLO_errorOption = 1;
#include "rf_allo.h"
#include "rf_io.h"
#include "std.h"

extern int ProcID;

#ifndef ALLIGNMENT_BOUNDARY
#define ALLIGNMENT_BOUNDARY 8
#endif

#define ALLOC_INTERFACE_ERROR -23456

#ifdef DEBUG_MEMORY
#define ALLOC_PROBLEM_ADDRESS 0x3662448
#endif
/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef	struct
 *       {	int	bus1;
 *              int	bus2;
 *              int	dest;
 *      }       POINT;
 *
 *      POINT    **points,corner;
 *
 *      points = (POINT **)array(2,x,y,sizeof(POINT));
 *                               ^ ^ ^
 *                               | | |
 *         number of dimensions--+ | |
 *                                 | |
 *          first dimension max----+ |
 *                                   |
 *         second dimension max------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free(points); (frees the entire structure in one fell swoop)
 *****************************************************************************/
/******************************************************************************
        The following section is a commented section containing
        an example main code:
*******************************************************************************
double *array_alloc();
main()
{
   int ***temp;
   int *temp2;
   int i, j, k;
   int il, jl, kl;

   malloc_debug(2);
   il = 2;
   jl = 3;
   kl = 3;
   temp = (int ***) array_alloc(3,il,jl,kl,sizeof(int));
   for (i=0; i<il; i++) {
      for (j=0; j<jl; j++) {
         for (k=0; k<kl; k++) temp[i][j][k] = 1;
      }
   }

   temp2 = (int *) malloc(10*sizeof(int));
   for (i=0; i<10; i++) temp2[i] = 0;

   if (Unlimited_Output) {
     for (i=0; i<il; i++) {
        for (j=0; j<jl; j++) {
           for (k=0; k<kl; k++) (void) fprintf(stderr," %d\n", temp[i][j][k]);
        }
     }
   }
   malloc_verify();
}
******************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef __STDC__

double *array_alloc(int numdim, ...)
#else

double *array_alloc(va_alist) va_dcl

#endif
/*****************************************************************************/

{
  int i;
  int j;
  struct dim {
    int index; /* Number of elements in the dimension	*/
    int total; /* Total number of elements 		*/
    int size;  /* Size of a single element in bytes	*/
    int off;   /* offset from beginning of array	*/
  } * dim;     /* Info about each dimension 		*/
#ifndef __STDC__
  int numdim; /* Number of dimensions			*/
#endif
  int total;      /* Total size of the array		*/
  double *dfield; /* ptr to avoid lint complaints		*/
  char *field;    /* The multi-dimensional array		*/
  char **ptr;     /* Pointer offset			*/
  char *data;     /* Data offset				*/
  va_list va;     /* Current pointer in the argument list	*/

  /*
   * Pull off the first argument.
   */

#ifdef __STDC__
  va_start(va, numdim);
#else
  va_start(va);
  numdim = va_arg(va, int);
#endif

  if (numdim <= 0) {
    (void)fprintf(stderr, "array_alloc ERROR: number of dimensions, %d, is <=0\n", numdim);
    exit(-1);
  }

  dim = (struct dim *)smalloc((int)(numdim * sizeof(struct dim)));

  dim[0].index = va_arg(va, int);
  /* Sanity checking */
  if (dim[0].index <= 0) {
    (void)fprintf(stderr, "WARNING: array_alloc called with first dimension <= 0, %d\n",
                  dim[0].index);
    (void)fprintf(stderr, "\tProc_%d: will return the nil pointer\n", ProcID);
  }
  dim[0].total = dim[0].index;
  dim[0].size = sizeof(void *);
  dim[0].off = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, int);
    /* Sanity checking */
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: array_alloc called with dimension %d <= 0, %d\n", i + 1,
              dim[i].index);
      fprintf(stderr, "\twill return the nil pointer\n");
    }
    dim[i].total = dim[i - 1].total * dim[i].index;
    dim[i].size = sizeof(void *);
    dim[i].off = dim[i - 1].off + dim[i - 1].total * dim[i - 1].size;
  }
  dim[numdim - 1].size = (int)va_arg(va, size_t);
  va_end(va);

  /*
   * Make sure that the data to be accessed is aligned appropriately on the
   * desired alignment boundary. Without this statement, array_alloc may
   * cause bus errors on some systems for multidimensioned arrays of doubles,
   * whose beginning dimensions are all odd. This is because the double
   * will be alligned on an address which is not a multiple of 8.
   */
  if (dim[numdim - 1].off > 0) {
    dim[numdim - 1].off =
        (((dim[numdim - 1].off - 1) / ALLIGNMENT_BOUNDARY + 1) * ALLIGNMENT_BOUNDARY);
  }

  total = dim[numdim - 1].off + dim[numdim - 1].total * dim[numdim - 1].size;
  dfield = (double *)smalloc(total);
  field = (char *)dfield;

  for (i = 0; i < numdim - 1; i++) {
    ptr = (char **)(field + dim[i].off);
    data = (char *)(field + dim[i + 1].off);
    for (j = 0; j < dim[i].total; j++) {
      ptr[j] = data + j * dim[i + 1].size * dim[i + 1].index;
    }
  }
  safe_free((void *)dim);
  return (dfield);
} /* END of routine array_alloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void fprint_location(const char *filename, const int line)

/************************************************************************
 *
 * fprint_location:
 *
 *    Central location to print file location string. Protect against
 *    NULL pointer
 *************************************************************************/
{
  if (filename) {
    (void)fprintf(stderr, " from %s: line %d\n", filename, line);
  } else {
    (void)fprintf(stderr, " from unknown location\n");
  }
  fflush(stderr);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void *safe_malloc(const int n, const char *filename, const int line)

/************************************************************************
 *
 * safe_malloc():
 *
 * Safe version of malloc().  Does not initialize memory. If given
 * a size of zero bytes, this program WILL NOT INITIALIZE MEMORY.
 * If an out of memory error is received by this routine, the
 * ENTIRE PROGRAM will exit. (Note: Goma isn't really set up to
 * recover from an out of memory error anyway).
 *
 * Renamed so smalloc() is a C preprocessor defn that passes down
 * file and line number information so the messages from here can
 * be more helpful.
 *
 * Return
 * --------
 *  Will return the pointer to the initialized malloced memory
 *
 *  Will write to the log file if the DEBUG_LEVEL preprocessor
 *   macro is turned on.
 *
 *  If the DEBUG_MEMORY define is turned on, this routine will write
 *  a message stating that it has malloced at a particular memory
 *  location hard-coded at compile time (sometimes useful when
 *  tracking down specific memory allocation errors such as freeing
 *  something twice).
 ***********************************************************************/
{
  void *pntr;
  size_t nn;
  if (n < 0) {
    fprintf(stderr, "smalloc ERROR P_%d: Negative memory specification (%d)", ProcID, n);
    fprint_location(filename, line);
    exit(-1);
  } else if (n == 0) {
    pntr = NULL;
  } else {
    nn = ((n - 1) / ALLIGNMENT_BOUNDARY + 1) * ALLIGNMENT_BOUNDARY;
    pntr = (void *)malloc(nn);
    if (pntr == NULL) {
      fprintf(stderr, "smalloc ERROR P_%d:  Memory allocation failure for %ld bytes", ProcID,
              (long int)nn);
      fprint_location(filename, line);
      exit(-1);
    }
#ifdef DEBUG_MEMORY
    if ((int)pntr == ALLOC_PROBLEM_ADDRESS) {
      (void)fprintf(stderr, "smalloc: FOUND address %x malloced at location", pntr);
      fprint_location(filename, line);
    }
#endif
  }
#if DEBUG_LEVEL > 0
  log_msg("%s:%d: allocating %d bytes @ %x\n", filename, line, n, pntr);
#endif
  return (pntr);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void safe_free(void *ptr)

/************************************************************************
 *
 * safe_free():
 *
 *    This version of free calls the system's free function with
 *    maximum error checking. It doesn't call free if ptr is the
 *    NULL pointer.
 ************************************************************************/
{
#ifdef DEBUG_MEMORY
  if ((int)ptr == ALLOC_PROBLEM_ADDRESS) {
    fprintf(stderr, "safe_free: FOUND IT!\n");
  }
#endif
  if (ptr != NULL) {
    free(ptr);
  }
#if DEBUG_LEVEL > 1
  else {
    log_err("safe_free: Attempt to free NULL ptr.");
  }
#endif
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void alloc_eh(char *rname, int bytes)

/**************************************************************************
 *
 * alloc_eh:
 *
 *      Error Handling for out of memory conditions:
 *
 *  Currently, the handling is determined by a global variable,
 *  ALLO_errorOption, whose scope is limited to this file.
 *  Depending upon its values, this function does various things:
 *
 *         3 create a divide by zero for stack trace analysis.
 *         2 or above -> print a message and exit the executing code.
 *         1 print a message and return with the NULL pointer
 *         0 Keep completely silent about the matter.
 *
 *  Currently, the variable is set to 1 as a default.
 ***************************************************************************/
{
  static const char yo[] = "alloc_eh";
  double cd = 0.0;
  if (ALLO_errorOption) {
    if (bytes == ALLOC_INTERFACE_ERROR) {
      fprintf(stderr, "alloc Interface ERROR P_%d: %s\n", ProcID, rname);
    } else {
      fprintf(stderr, "%s ERROR P_%d: out of memory while mallocing %s %u bytes\n", yo, ProcID,
              rname, bytes);
    }
  }
  fflush(stderr);
  if (ALLO_errorOption == 3)
    cd = 1.0 / cd;
  if (ALLO_errorOption > 1)
    exit(-1);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void safer_free(void **ptr)

/***************************************************************************
 *
 * safer_free():
 *
 *    Wrapper around the free() function: This takes the handle to the
 *    pointer. So that if the pointer value is changed, then the
 *    change can be passed back to the calling program. In this way, we can
 *    ensure that pointers which don't point to valid memory locations
 *    always have the value of NULL.
 ***************************************************************************/
{
  if (ptr == NULL) {
    alloc_eh("NULL handle passed to safer_free()", ALLOC_INTERFACE_ERROR);
    return;
  }
  if (*ptr != NULL) {
#ifdef DEBUG_MEMORY
    if ((int)*ptr == ALLOC_PROBLEM_ADDRESS) {
      fprintf(stderr, "safer_free: FOUND IT!\n");
    }
#endif
    free(*ptr);
    *ptr = NULL;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void zero_structure(void *ptr, const size_t struct_size, const int num_structs)

/****************************************************************************
 *
 * zero_structure:
 *
 *    Intializes a contiguous portion of memory, which constitues an
 *    array of structures to zeroes. This has
 *    the effect of setting all null pointers to zero, all doubles
 *    to the value of zero, and all int numbers to zero within the
 *    structure.
 *
 *  Input
 * ---------
 *  ptr         Pointer to the first structure
 *  struct_size Size of the structure
 *  num_structs Number of structures continguous in memory to be zeroed.
 ***************************************************************************/
{
  int bytes;
  if (ptr != NULL) {
    bytes = struct_size * num_structs;
    if (bytes > 0)
      (void)memset(ptr, 0, bytes);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int *alloc_int_1_FL(const int nvalues, const int val, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_int_1_FL:
 *  alloc_int_1(const int nvalues, const int val):
 *
 *    Allocate and initialize a one dimensional array of ints.
 *    This routine always allocates space for at least one int.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = intialization value
 *    Return
 *    ------
 *        Pointer to the intialized array of ints.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int i, bytenum, nval1 = nvalues;
  int *array;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_int_1 P_%d: called with n = %d ", ProcID, nvalues);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(int);
  array = (int *)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (val == 0) {
      (void)memset((void *)array, 0, bytenum);
    } else if (val != INT_NOINIT) {
      for (i = 0; i < nval1; i++)
        array[i] = val;
    }
  } else {
    alloc_eh("alloc_int_1", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

short int *alloc_short_1_FL(const int nvalues, const int val, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_short_1_FL:
 *  alloc_short_1(const int nvalues, const int val):
 *
 *    Allocate and initialize a one dimensional array of short ints.
 *    This routine always allocates space for at least one short int.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = intialization value
 *    Return
 *    ------
 *        Pointer to the intialized array of short ints.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int i, bytenum, nval1 = nvalues;
  short int *array, val_short;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_short_1 P_%d: called with n = %d ", ProcID, nvalues);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(short int);
  array = (short int *)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (val == 0) {
      (void)memset((void *)array, 0, bytenum);
    } else if (val != INT_NOINIT) {
      if (val < SHRT_MAX && val > SHRT_MIN) {
        val_short = (short int)val;
        for (i = 0; i < nval1; i++)
          array[i] = val_short;
      } else {
        fprintf(stderr, "Error: alloc_short_1 P_%d: called with val out of bounds, %d", ProcID,
                val);
        fprint_location(filename, line);
        alloc_eh("alloc_short_1 outofbounds", bytenum);
      }
    }
  } else {
    alloc_eh("alloc_short_1", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int **alloc_int_2_FL(
    const int ndim1, const int ndim2, const int val, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_int_2_FL:
 *  alloc_int_2(const int ndim1, const int ndim2, const int val):
 *
 *    Allocate and initialize a two dimensional array of ints.
 *    This routine always allocates space for at least one int.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        ndim1 = Length of first dimension of the array (outer)
 *        ndim2 = Length of the second dimension of the array (inner)
 *        val   = intialization value
 *    Return
 *    ------
 *        Pointer to the intialized array of ints.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int i, bytenum, nval1 = ndim1, nval2 = ndim2;
  int **array, *iloc;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_int_2 P_%d: called with ndim1 = %d ", ProcID, ndim1);
      fprint_location(filename, line);
    }
  }
  if (nval2 <= 0) {
    nval2 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_int_2 P_%d: called with ndim2 = %d ", ProcID, ndim2);
      fprint_location(filename, line);
    }
  }
  array = (int **)array_alloc(2, nval1, nval2, sizeof(int));
  if (array != NULL) {
    iloc = array[0];
    if (val == 0) {
      bytenum = sizeof(int) * nval1 * nval2;
      (void)memset((void *)iloc, 0, bytenum);
    } else if (val != INT_NOINIT) {
      bytenum = nval1 * nval2;
      for (i = 0; i < bytenum; i++)
        iloc[i] = val;
    }
  } else {
    bytenum = sizeof(int) * nval1 * nval2;
    alloc_eh("alloc_int_2", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double *alloc_dbl_1_FL(const int nvalues, const double val, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_dbl_1_FL:
 *  alloc_dbl_1(const int nvalues, const double val):
 *
 *    Allocate and initialize a one dimensional array of doubles.
 *    This routine always allocates space for at least one double.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = intialization value (default = DBL_NOINIT)
 *    Return
 *    ------
 *        Pointer to the intialized array of doubles
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int i, bytenum, nval1 = nvalues;
  double *array;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      (void)fprintf(stderr, "Warning: alloc_dbl_1 P_%d: called with n = %d ", ProcID, nvalues);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(double);
  array = (double *)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (val == 0.0) {
      (void)memset((void *)array, 0, bytenum);
    } else if (val != DBL_NOINIT) {
      for (i = 0; i < nval1; i++)
        array[i] = val;
    }
  } else {
    alloc_eh("alloc_dbl_1", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double **alloc_dbl_2_FL(
    const int ndim1, const int ndim2, const double val, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_dbl_2_FL:
 *  alloc_dbl_2(const int ndim1, const int ndim2, const double val):
 *
 *    Allocate and initialize a two dimensional array of doubles.
 *    This routine always allocates space for at least one double.
 *    Thus, it will ensure that ndim1 and ndim2 are greater than or
 *    equal to 1.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        ndim1 = Length of first dimension of the array (outer)
 *        ndim2 = Length of the second dimension of the array (inner)
 *        val   = intialization value
 *    Return
 *    ------
 *        Pointer to the intialized array of doubles.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int i, bytenum = 0, nval1 = ndim1, nval2 = ndim2;
  double **array, *iloc;
  if (nval1 <= 0) {
    nval1 = 1;
    fprintf(stderr, "Warning: alloc_dbl_2 P_%d: called with ndim1 = %d", ProcID, ndim1);
    fprint_location(filename, line);
  }
  if (nval2 <= 0) {
    nval2 = 1;
    fprintf(stderr, "Warning: alloc_dbl_2 P_%d: called with ndim2 = %d", ProcID, ndim2);
    fprint_location(filename, line);
  }
  array = (double **)array_alloc(2, nval1, nval2, sizeof(double));
  if (array != NULL) {
    iloc = array[0];
    if (val == 0.0) {
      bytenum = sizeof(double) * nval1 * nval2;
      (void)memset((void *)iloc, 0, bytenum);
    } else if (val != DBL_NOINIT) {
      bytenum = nval1 * nval2;
      for (i = 0; i < bytenum; i++)
        iloc[i] = val;
    }
  } else {
    alloc_eh("alloc_dbl_2", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void *alloc_void_struct_1_FL(const size_t structsize,
                             const int nvalues,
                             const char *filename,
                             const int line)

/**************************************************************************
 *
 *  alloc_void_struct_1_FL:
 *  alloc_void_struct_1(const size_t structsize, const int nvalues)
 *  alloc_struct_1(structure_name, nvalues)
 *
 *    Allocate and initialize a one dimensional array of
 *    structures. Intialization to zero is always done.
 *    This routine always allocates space for at least one structure.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location.
 *
 *    Input
 *    -------
 *        structsize = Size of each of the structures
 *        nvalues = number of the structures
 *
 *    Return
 *    ------
 *        Pointer to the intialized array of structures
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  size_t bytenum = 0;
  void *array = NULL;
  if (structsize == 0) {
    if (Unlimited_Output) {
      fprintf(stderr, "alloc_void_struct_1 WARNING size of structure is zero\n");
    }
    alloc_eh("alloc_void_struct_1", bytenum);
    return NULL;
  }
  if (nvalues <= 0) {
    bytenum = structsize;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_void_struct_1 P_%d: called with n = %d ", ProcID, nvalues);
      fprint_location(filename, line);
    }
  } else {
    bytenum = nvalues * structsize;
  }
  array = safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    (void)memset(array, 0, bytenum);
  } else {
    alloc_eh("alloc_void_struct_1", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void **alloc_ptr_1_FL(const int num_ptrs, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_ptr_1_FL:
 *  alloc_ptr_1(const int num_ptrs):
 *
 *    Allocate and initialize a one dimensional array of pointers to void.
 *    This routine always allocates space for at least one pointer.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location. This routine also always initializes
 *    all pointers to the NULL value.
 *
 *    Input
 *    -------
 *        num_ptrs = Length of the array
 *    Return
 *    ------
 *        Pointer to the intialized array of pointers.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int nval1 = num_ptrs;
  size_t bytenum;
  void **array = NULL;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_ptr_1 P_%d: called with n = %d ", ProcID, num_ptrs);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(void *);
  array = (void **)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    (void)memset(array, 0, bytenum);
  } else {
    alloc_eh("alloc_ptr_1", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void realloc_int_1_FL(int **list_hndl,
                      const int new_length,
                      const int old_length,
                      const char *filename,
                      const int line)

/**************************************************************************
 *
 *  realloc_int_1_(list_hndl, new_num_ptrs, old_num_ptrs);
 *  realloc_int_1_FL():
 *
 *    Reallocates a one dimensional array of ints.
 *    This routine always allocates space for at least one int.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        list_undl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of integers to be reallocated
 *        new_length   = Length of the array
 *        old_length   = Length of the old array
 ***************************************************************************/
{
  int nval1 = new_length;
  size_t bytenum;
  int *array = NULL;
  if (nval1 == old_length)
    return;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: realloc_int_1 P_%d: called with n = %d ", ProcID, new_length);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(int);
  array = (int *)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (old_length > 0) {
      bytenum = sizeof(int) * old_length;
    } else {
      bytenum = 0;
    }
    if (new_length < old_length) {
      bytenum = sizeof(int) * new_length;
    }
    if (bytenum > 0) {
      (void)memcpy((void *)array, (void *)*list_hndl, bytenum);
    }
    safe_free(*list_hndl);
    *list_hndl = array;
  } else {
    alloc_eh("realloc_int_1", bytenum);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void realloc_dbl_1_FL(double **dvec_hndl,
                      const int new_length,
                      const int old_length,
                      const char *filename,
                      const int line)

/**************************************************************************
 *
 *  realloc_dbl_1_(dvec_hndl, new_num_ptrs, old_num_ptrs);
 *  realloc_dbl_1_FL():
 *
 *    Reallocates a one dimensional array of doubles.
 *    This routine always allocates space for at least one double.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        dvec_hndl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of doubles to be reallocated
 *        new_length   = Length of the array
 *        old_length   = Length of the old array
 ***************************************************************************/
{
  int nval1 = new_length;
  size_t bytenum;
  double *array = NULL;
  if (nval1 == old_length)
    return;
  if (nval1 <= 0) {
    nval1 = 1;
    fprintf(stderr, "Warning: realloc_dbl_1 P_%d: called with n = %d ", ProcID, new_length);
    fprint_location(filename, line);
  }
  bytenum = nval1 * sizeof(double);
  array = (double *)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (old_length > 0) {
      bytenum = sizeof(double) * old_length;
    } else {
      bytenum = 0;
    }
    if (new_length < old_length) {
      bytenum = sizeof(double) * new_length;
    }
    if (bytenum > 0) {
      (void)memcpy((void *)array, (void *)*dvec_hndl, bytenum);
    }
    safe_free(*dvec_hndl);
    *dvec_hndl = array;
  } else {
    alloc_eh("realloc_dbl_1", bytenum);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void realloc_ptr_1_FL(void ***list_hndl,
                      const int new_num_ptrs,
                      const int old_num_ptrs,
                      const char *filename,
                      const int line)

/**************************************************************************
 *
 *  realloc_ptr_1_(list_hndl, new_num_ptrs, old_num_ptrs);
 *  realloc_ptr_1_FL():
 *
 *    Reallocates a one dimensional array of pointers to void.
 *    This routine always allocates space for at least one pointer.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        list_undl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of pointers to be realloced
 *        new_num_ptrs = Length of the array
 *        old_num_ptrs = Length of the old array
 ***************************************************************************/
{
  int nval1 = new_num_ptrs;
  size_t bytenum;
  void **array = NULL;
  if (nval1 == old_num_ptrs)
    return;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: realloc_ptr_1 P_%d: called with n = %d ", ProcID, new_num_ptrs);
      fprint_location(filename, line);
    }
  }
  bytenum = nval1 * sizeof(void *);
  array = (void **)safe_malloc(bytenum, filename, line);
  if (array != NULL) {
    if (old_num_ptrs > 0) {
      bytenum = sizeof(void *) * old_num_ptrs;
    } else {
      bytenum = 0;
    }
    if (new_num_ptrs < old_num_ptrs) {
      bytenum = sizeof(void *) * new_num_ptrs;
    }
    if (old_num_ptrs > 0) {
      (void)memcpy(array, *list_hndl, bytenum);
    }
    if (new_num_ptrs > old_num_ptrs) {
      bytenum = (new_num_ptrs - old_num_ptrs) * sizeof(void *);
      (void)memset(array + old_num_ptrs, 0, bytenum);
    }
    safe_free(*list_hndl);
    *list_hndl = array;
  } else {
    alloc_eh("realloc_ptr_1", bytenum);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void *realloc_void_struct_1_FL(void *vec_ptr,
                               const size_t structsize,
                               const int new_length,
                               const int old_length,
                               const char *filename,
                               const int line)

/**************************************************************************
 *
 *  realloc_void_struct_1_FL(vec_ptr, size, new_num_ptrs, old_num_ptrs);
 *  realloc_struct_1(vec_ptr, struct name, new_length, old_length):
 *
 *    Reallocates a one dimensional array of structures with size structsize.
 *    This routine always allocates space for at least one struct.
 *    Calls the safe_malloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        vec_ptr      = Pointer to the first element of the structure
 *                       array to be reallocated
 *        structsize   = number of bytes in each structure array member
 *        new_length   = Length of the array
 *        old_length   = Length of the old array
 *    Output
 *    -------
 *        returns pointer to new structure array as void *.
 *                Calling function responsible for appropriate cast.
 ***************************************************************************/
{
  int nval1 = new_length;
  size_t bytenum;
  void *array = NULL;
  if (nval1 == old_length)
    return (vec_ptr);
  if (nval1 <= 0) {
    nval1 = 1;
    fprintf(stderr, "Warning: realloc_struct_1 P_%d: called with n = %d ", ProcID, new_length);
    fprint_location(filename, line);
  }
  bytenum = nval1 * structsize;
  array = safe_malloc(bytenum, filename, line);
  memset(array, 0, bytenum);
  if (array != NULL) {
    if (old_length > 0) {
      bytenum = structsize * old_length;
    } else {
      bytenum = 0;
    }
    if (new_length < old_length) {
      bytenum = sizeof(double) * new_length;
    }
    (void)memcpy(array, vec_ptr, bytenum);
    safe_free(vec_ptr);
  } else {
    alloc_eh("realloc_dbl_1", bytenum);
  }
  return (array);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ***alloc_ptr_2_FL(const int numDim1, const int numDim2, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_ptr_2_FL:
 *  alloc_ptr_2(const int numDim1, const int numDim2):
 *
 *    Allocate and initialize a two dimensional array of pointers to void.
 *    This routine always allocates space for at least one pointer in each
 *    dimension.
 *    Calls the array_alloc() routine which calls the safe malloc routine
 *    to ensure that all malloc
 *    calls go through one location. This routine also always initializes
 *    all pointers to the NULL value.
 *
 *    Input
 *    -------
 *        numDim1 = number of pointers in the first dimension
 *        numDim2 = number of pointers in the second dimension
 *    Return
 *    ------
 *        Pointer to the intialized array of pointers.
 *        Failures are first handled by calling an error handler, alloc_eh().
 *        The error handler may or may not print a message, and may or
 *        or may not terminate the program. If the error handler doesn't
 *        terminate the program, then control is returned to the calling
 *        program, with the error indicated by returning the NULL pointer.
 ***************************************************************************/
{
  int nval1 = numDim1;
  int nval2 = numDim2;
  int bytenum;
  void ***array;
  if (nval1 <= 0) {
    nval1 = 1;
    if (Unlimited_Output) {
      (void)fprintf(stderr, "Warning: alloc_ptr_2 P_%d: called with numDim1 = %d ", ProcID,
                    numDim1);
      fprint_location(filename, line);
    }
  }
  if (nval2 <= 0) {
    nval2 = 1;
    if (Unlimited_Output) {
      fprintf(stderr, "Warning: alloc_ptr_2 P_%d: called with numDim2 = %d ", ProcID, numDim2);
      fprint_location(filename, line);
    }
  }
  array = (void ***)array_alloc(2, nval1, nval2, sizeof(void *));
  bytenum = nval1 * nval2 * sizeof(void *);
  if (array != NULL) {
    (void)memset(array[0], 0, bytenum);
  } else {
    alloc_eh("alloc_ptr_2", bytenum);
  }
  return array;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char **alloc_VecFixedStrings(const int numStringsC, const int lenStringC)

/**************************************************************************
 *
 *  alloc_VecFixedStrings:
 *
 *    Allocate and initialize a vector of fixed-length
 *    strings. Each string is initialized to the NULL string.
 *    Something is always malloced by this routine, since arguments
 *    are increased to 1, if they are zero or less.
 *
 *    Input
 *    -------
 *        numStringsC = Number of strings
 *        lenStringC  = Length of each string including the trailing null
 *                      character
 *    Return
 *    ------
 *        This value is initialized to the correct address of the array.
 *        A NULL value in the position indicates an error.
 ***************************************************************************/
{
  int i, numStrings, lenString;
  char **array;
  if (numStringsC <= 0)
    numStrings = 1;
  else
    numStrings = numStringsC;
  if (lenStringC <= 0)
    lenString = 1;
  else
    lenString = lenStringC;
  array = (char **)array_alloc(2, numStrings, lenString, sizeof(char));
  if (array != NULL) {
    for (i = 0; i < numStrings; i++)
      array[i][0] = '\0';
  } else {
    alloc_eh("alloc_VecFixedStrings",
             sizeof(char) * numStrings * lenString + numStrings * sizeof(void *));
  }
  return array;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char *alloc_copy_string_FL(const char *copyFrom, const char *filename, const int line)

/**************************************************************************
 *
 *  alloc_copy_string():
 *  alloc_copy_string_FL():
 *
 *    Allocates space for and then copies a string
 *
 *    Input
 *    -------
 *        copyFrom = null terminated string. If NULL is supplied, then
 *                   nothing is malloced and a NULL value is returned.
 *    Return
 *    ------
 *        This value is initialized to the correct address of the new
 *        string.
 *        A NULL value in the position either indicates an error, or
 *        that the original pointer to the string was NULL.
 ***************************************************************************/
{
  char *cptr;
  int bytenum;
  if (copyFrom == NULL)
    return NULL;
  bytenum = (strlen(copyFrom) + 1) * sizeof(char);
  cptr = (char *)safe_malloc(bytenum, filename, line);
  if (cptr == NULL) {
    alloc_eh("alloc_copy_string", bytenum);
  }
  return strcpy(cptr, copyFrom);
}

void checkFinite(double tmp) {
  if (!isfinite(tmp)) {
    if (isnan(tmp)) {
      printf("ERROR: we have encountered a nan!\n");
    } else if (isinf(tmp) == 1) {
      printf("ERROR: we have encountered a pos inf!\n");
    } else {
      printf("ERROR: we have encountered a neg inf!\n");
    }
    exit(-1);
  }
}

/*****************************************************************************/
/*                       END of rf_allo.c                                    */
/*****************************************************************************/
