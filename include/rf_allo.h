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
 *$Id: rf_allo.h,v 5.4 2008-10-02 15:36:29 hkmoffa Exp $
 */

#ifndef GOMA_RF_ALLO_H
#define GOMA_RF_ALLO_H

#ifdef EXTERN
#undef EXTERN
#endif

#include <stdlib.h>

#include "mm_eh.h"

#ifndef GOMA_RF_ALLO_C
#define EXTERN extern
#else
#define EXTERN
#endif

/*
 *  These are a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers which can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define INT_NOINIT                      -68361
#define DBL_NOINIT                      -34.39383E11

/*
 * Wrapper helps identify where allocations take place and for how much.
 * Someday need to wrap varargs array_alloc() too.
 *  -> just put the two extra arguments at the beginning of the argument
 *     list.
 */
#define smalloc(arg)                    safe_malloc((arg), __FILE__, __LINE__)
#define alloc_ptr_1(arg)                alloc_ptr_1_FL((arg), __FILE__, __LINE__)
#define alloc_ptr_2(arg1, arg2)         alloc_ptr_2_FL((arg1), (arg2), __FILE__, __LINE__)
#define alloc_int_1(nvalues, value)     alloc_int_1_FL((nvalues), (value), __FILE__, __LINE__)
#define alloc_short_1(arg1, arg2)       alloc_short_1_FL((arg1), (arg2), __FILE__, __LINE__)
#define alloc_dbl_1(arg1, arg2)         alloc_dbl_1_FL((arg1), (arg2), __FILE__, __LINE__)
#define alloc_void_struct_1(arg1, arg2) alloc_void_struct_1_FL((arg1), (arg2), __FILE__, __LINE__)
#define alloc_struct_1(x, num)          (x *)alloc_void_struct_1_FL(sizeof(x), (num), __FILE__, __LINE__)

#define zeroStructures(xptr, num)       zero_structure((void *)(xptr), sizeof((*xptr)), (num))
#define realloc_ptr_1(arg1, arg2, arg3) realloc_ptr_1_FL((arg1), (arg2), (arg3), __FILE__, __LINE__)
#define realloc_int_1(arg1, arg2, arg3) realloc_int_1_FL((arg1), (arg2), (arg3), __FILE__, __LINE__)
#define realloc_dbl_1(arg1, arg2, arg3) realloc_dbl_1_FL((arg1), (arg2), (arg3), __FILE__, __LINE__)
#define realloc_struct_1(arg1, x, arg3, arg4) \
  (x *)realloc_void_struct_1_FL((arg1), sizeof(x), (arg3), (arg4), __FILE__, __LINE__)
#define alloc_int_2(arg1, arg2, arg3) alloc_int_2_FL((arg1), (arg2), (arg3), __FILE__, __LINE__)
#define alloc_dbl_2(arg1, arg2, arg3) alloc_dbl_2_FL((arg1), (arg2), (arg3), __FILE__, __LINE__)
#define alloc_copy_string(arg1)       alloc_copy_string_FL((arg1), __FILE__, __LINE__)

#define copy_dbl_1(copyTo, copyFrom, len) \
  (dbl *)memcpy((void *)(copyTo), (const void *)(copyFrom), (sizeof(dbl) * len))

#define copy_int_1(copyTo, copyFrom, len) \
  (int *)memcpy((void *)(copyTo), (const void *)(copyFrom), (sizeof(int) * len))

#define zero_dbl_1(copyTo, len) memset((void *)(copyTo), 0, (sizeof(dbl) * len))
#define zero_int_1(copyTo, len) memset((void *)(copyTo), 0, (sizeof(int) * len))
/*
 * Prototypes for functions in rf_allo.c
 */
EXTERN double *array_alloc(int,        /* numdim - number of dimensions */
                           ...);       /* all remaing varargs, last 2 are file:line */

EXTERN void *safe_malloc(const int,    /* numbytes */
                         const char *, /* filename */
                         const int);   /* line */

EXTERN void safe_free(void *);         /* ptr to block being freed */

extern void safer_free(void **);
extern int *alloc_int_1_FL(const int, const int, const char *, const int);
extern short int *alloc_short_1_FL(const int, const int, const char *, const int);
extern double *alloc_dbl_1_FL(const int, const double, const char *, const int);
extern void *alloc_void_struct_1_FL(const size_t, const int, const char *, const int);
extern void zero_structure(void *, const size_t, const int);
extern void **alloc_ptr_1_FL(const int, const char *, const int);
extern void realloc_ptr_1_FL(void ***, const int, const int, const char *, const int);
extern void ***alloc_ptr_2_FL(const int, const int, const char *, const int);
extern int **alloc_int_2_FL(const int, const int, const int, const char *, const int);
extern double **alloc_dbl_2_FL(const int, const int, const double, const char *, const int);
extern char **alloc_VecFixedStrings(const int, const int);
extern char *alloc_copy_string_FL(const char *, const char *, const int);

extern void realloc_int_1_FL(int **, const int, const int, const char *, const int);
extern void realloc_dbl_1_FL(double **, const int, const int, const char *, const int);
extern void *
realloc_void_struct_1_FL(void *, const size_t, const int, const int, const char *, const int);
extern void checkFinite(double tmp);

#endif
