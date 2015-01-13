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
 *$Id: std.h,v 5.3 2010-07-01 17:29:33 ebenner Exp $
 */

#ifndef _RF_GOMA_H
#define _RF_GOMA_H

#include <ctype.h>

/* If we're using autoconf, then include the config.h file. */
#ifdef HAVE_CONFIG_H
#include "config.h"

/*
 * If not, make sure some of these HAVE_* are def'ed for the test suite.
 * The test suite uses these -- via 'goma -v' -- to check for various
 * capabilities which may or may not be compiled in.
 */
#else /* HAVE_CONFIG_H */

#define HAVE_SPARSE 1
#define HAVE_UMFPACK 1
#define HAVE_AZTEC 1
/* Trilinos now seems to define these */
#ifndef HAVE_BLAS
  #define HAVE_BLAS 1
#endif
#ifndef HAVE_LAPACK
  #define HAVE_LAPACK 1
#endif
#define HAVE_Y12M 1

#ifdef EIGEN_SERIAL
#define HAVE_ARPACK 1
#endif /* EIGEN_SERIAL */

#ifdef EIGEN_PARALLEL
#define HAVE_ARPACK 1
#define HAVE_PARPACK 1
#endif /* EIGEN_PARALLEL */

#endif /* HAVE_CONFIG_H */

#ifndef VERSION
#define VERSION "6.0.0"
#endif

/*****************************************************************************/
/*
 *                  PROTOTYPES SECTION
 *
 *        Define the availability of prototypes. 
 */

#ifndef HAVE_PROTOTYPES
# if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#     define	HAVE_PROTOTYPES
# endif
#endif

#undef PROTO

#ifdef HAVE_PROTOTYPES 
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

/***********************************************************************************/
/*
 *                   MACHINE NAME SECTION
 *
 *   Section to determine standard names for machines
 * The idea is to localize  machine name determinations in this one section of the
 * code.
 *
 *  Possible values:
 *
 *        solaris
 *        ASCI_RED
 *        hpux
 *        tflop
 *        aix
 *        BSD
 */

#ifndef solaris
#if defined(__sun) && defined(__SVR4)
#define solaris 1
#endif
#endif

#if defined(__PUMAGON__)
#define ASCI_RED 1
#endif

#ifndef hpux
#ifdef __hpux
#define hpux
#endif
#endif

#if 0
#if defined(__PARAGON__)
#define tflop
#endif
#endif

#if defined (hpux) || defined (ASCI_RED)
#ifdef solaris
#undef solaris
#endif
#endif

/*****************************************************************************/
/*
 *   HAVE_SUNMATH_H -- have the sunmath.h include file in the default search
 *                     path for include files.
 */

#ifdef solaris
#define HAVE_SUNMATH_H
#endif

/*****************************************************************************/
/*
 *      SECTION TO DEFINE SOME BASIC MATH MACROS AND MATH SYMBOLS FOR CONSTANTS
 */

/*
 * Including this for hpux removes two compiler warnings for
 * redeclarations of MAX and MIN below
 */
#ifdef solaris
#include <ieeefp.h>
#endif

#ifdef hpux
#include <sys/param.h>
#endif

#ifndef MAX
#define MAX(x,y)     (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y)     (( (x) < (y) ) ? (x) : (y))
#endif

#ifndef SGN
#define SGN(x)               (( (x) < 0 ) ? -1 : 1 )
#endif

#ifndef ABS
#define ABS(x)               (( (x) < 0 ) ? -(x) : (x) )
#endif

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif
/*
 * This one lets sign(0) = -1...
 */

#define sign_of(x)	((x) <=  0. ?  -1 : 1 )

#ifndef SQUARE
#define SQUARE(x)	( (x) * (x) )
#endif
#ifndef CUBE
#define CUBE(x)	( (x) * (x) * (x) )
#endif

/*
 *  cbrt() is a common math function. however it is not an iso c function.
 *  Therefore, we provide a workaround here using the pow() function.
 */
#define HAVE_CBRT

#ifndef HAVE_CBRT
#define cbrt(arg)	( arg == 0 ? 0 : ( arg < 0. ? -pow(-arg,1./3.) : pow(arg, 1./3.) ))
#endif

/*
 * Some handy large integers...
 */

#define KILO_2                          1024
#define MEGA_2                          1048576
#define GIGA_2                          1073741824
#define TERA_2                          1099511627776
#define PETA_2                          1125899906842624
#define EXA_2                           1152921504606846976

#ifndef M_PIE
#define M_PIE 3.14159265358979323846
#endif

#define delta(m,n)	((m) == (n) ? 1 : 0 ) /* Kroenecker delta */
#define permute(i,j,k) ( ((i)-(j))*((j)-(k))*((k)-(i))/2 ) /* Permutation symbol (epsilon) */
#define stringup(a) { char *p; for( p=a; *p != '\0'; *p=toupper(*p), p++); }
#define stringlow(a) { char *p; for( p=a; *p != '\0'; *p=tolower(*p), p++); }
#define endofstring(a) strchr(a,'\0')
/*
 * Comparisons against zero should be done using this double value
 * in order to avoid machine dependent issues involving how
 * underflows are handled.
 */
#ifndef DBL_SMALL
#define DBL_SMALL 1.0E-298
#endif

#ifdef sgi
#ifdef jkc_pg
#include <limits.h>
#else
#include <internal/limits_core.h>
#endif
#endif
#ifndef DBL_MAX
#define DBL_MAX 1.0E300
#endif

/*
 * This definition has a result of true if the double precision argument is
 * nonzero, non-zero underflow, and false otherwise
 */
#ifndef DOUBLE_NONZERO
#define DOUBLE_NONZERO(value)  ((fabs((value)) > DBL_SMALL))
#endif
#ifndef DOUBLE_ZERO
#define DOUBLE_ZERO(value)  ((fabs((value)) < DBL_SMALL))
#endif

/*************************************************************************/
/*           DEFAULT STRING AND MALLOC SIZE SECTION
 *
 *     This section is devoted to specifying the default lengths
 *     of strings and malloced objects
 */

/*
 * How long are the strings that sprintf fills with more detailed error message
 * information?
 */

#ifndef MAX_CHAR_ERR_MSG
#  define MAX_CHAR_ERR_MSG	1024
#endif

/*
 * Allocate and reallocate memory in chunks this large.
 */

#ifndef LIST_CHUNK_SIZE
#define LIST_CHUNK_SIZE		100
#endif

/*************************************************************************/
/*           SOLVER AVAILABILITY SECTION
 *
 *     This section is devoted to specifying the default availability
 *     of optional solvers.
 */

/*
 * By default assume the Aztec linear solver package is available.
 * If not, then comment out this line.
 */
#ifndef HAVE_AZTEC
#define HAVE_AZTEC
#endif


/*************************************************************************/
/*            MPI IMPLEMENTATION SECTION
 *
 *
 *     If PARALLEL is defined, then we assume that the mpi.h include file
 *     is available. Also we turn on the parallel implementation of the
 *     cpc package.
 */

#ifdef PARALLEL
#ifndef MPI
#define MPI
#endif /* MPI */
#define CPC_PARALLEL
#define HAVE_MPI_H 1
#endif /* PARALLEL */

/*
 * HAVE_MPE_H
 *     Some implementations of mpi don't have the mpe package. We
 *     differentiate these by settting and unsetting the HAVE_MPE_H
 *     symbol
 */

#if defined(_AIX) || defined(tflop)
#undef HAVE_MPE_H
#endif

#ifdef HAVE_MPI_H
#include <mpi.h>
#endif
#ifdef HAVE_MPE_H
#include <mpe.h>
#endif

/*
 *    If MPI is not used, then we need to define some MPI data types
 *    used liberally in GOMA outside of the MPI package
 */

#ifndef HAVE_MPI_H
typedef int MPI_Datatype;
typedef int MPI_Request;
typedef int MPI_Status;
typedef int MPI_Aint;
#define MPI_CHAR 1
#define MPI_INT 2
#define MPI_FLOAT 3
#define MPI_DOUBLE 4
#endif

/*************************************************************************/
/*
 *               DEBUGGING SECTION
 *
 *
 *     Use this section to specify debugging parameters.
 */
 /*
 * Use this to conditionally compile code. The idea is that you can compile
 * with more V&V that might slow and bloat your code, but can be toggled on
 * and off with, ultimately, parameters like Debug_Level as the first
 * argument to log_dbg(). In essence, then, log_dbg() can be used to print
 * out as much as you specify, provided you have compiled in that level of
 * checking.
 */

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif


/*
 * For parallel processing, restrict most diagnostic output to a single
 * processor.
 */

#ifdef PARALLEL
#  define DPRINTF if ( ProcID == 0 ) fprintf
#  define P0PRINTF if (ProcID == 0) printf
#else
#  define DPRINTF fprintf
#  define P0PRINTF printf
#endif


/****************************************************************************/
/*
 *                   SYSTEM FUNCTION SECTION
 */
/*
 *            Use this section to define differences in system
 *            functions between various platforms.
 *
 */

#ifdef BSD
typedef char *Spfrtn;
#else
typedef int Spfrtn;		/* most System V systems... */
#endif


/*
 * HKM -> I'm not sure why the following is in there
 *        iso c standard has set strcpy() to return a char pointer.
 */
#define STRCPY_RTN_IS_STRING
#ifdef STRCPY_RTN_IS_STRING
typedef char *Strcpy_rtn;
#endif

#ifdef STRCPY_RTN_IS_INT
typedef int  Strcpy_rtn;
#endif

/****************************************************************************/
/*                   	TYPEDEF DECLARATIONS SECTION
 *
 *
 *         Use this section to define typedef's needed throughout goma.
 */
/*
 *    We use typedefs for several unsigned variable types:
 * The following names are needed by Goma, and can usually be found in
 * <sys/types.h>
 *
 *      u_short
 *      u_int
 *      u_long
 */
#include <sys/types.h>

/*
 *  There are several machine/compiler cases where u_int and u_long are not defined
 *  in  <sys/types.h>. Take care of these special cases here.
 */

#ifndef tflop
#ifndef dec_osf1
#if !defined(_H_TYPES) && !defined(_SYS_TYPES_INCLUDED) && !defined(_SYS_TYPES_H) && !defined(__sys_types_h) && !defined(_SYS_TYPES_H_)
typedef unsigned       int		u_int;
typedef unsigned long  int		u_long;
#endif
#endif
#endif

/*
 * This might ultimately be handy for changing the precision of all the
 * floating point values used in the calculations. However, until ALL the
 * variable declarations make use of this type, such usefulness must be
 * a dream.
 */

typedef	double	dbl;

/*
 * A boolean type named bool
 */
#ifndef __cplusplus
typedef unsigned int bool;
#endif

enum datatype 
{
  type_void, 
  type_char, 
  type_unsigned_char,
  type_short_int, 
  type_int, 
  type_unsigned_int,
  type_long_int, 
  type_float, 
  type_double, 
  type_long_long, 
  type_long_double
};

typedef enum datatype Edt;

/*
 *  typedef shortcuts for structures defined later in other include files.
 */

typedef struct Node_Struct     NODE_STRUCT;
typedef struct Node_Var_Struct NODE_VAR_STRUCT;
typedef struct Var_Struct      VAR_STRUCT;
typedef struct Material_Properties MATRL_PROP_STRUCT;

/***************************************************************************/
/*                  PROTOTYPES SECTION
 *
 *     Use this section to define prototyes used in goma.
 *     Variables declared here are defined elsewhere.
 *
 */
/*
 * Global error signal flags. These are used mostly for sending error signals
 * when running in parallel so that all processors can be shut down
 * gracefully if any one of them goes south.
 */

/*
 * This is a generic error that is anything that a call to EH(-1... would
 * trigger. It gets checked several places...
 */

extern int parallel_err;
extern int parallel_err_global;

/*
 * This particular error is triggered by moving mesh problems that get
 * too enthusiastic and require spending the night in the drunk tank -
 * let the whole town know.
 */

extern int neg_elem_volume;
extern int neg_elem_volume_global;


/*
 * This particular error is triggered by moving mesh problems in shell that get
 * too enthusiastic and require spending the night in the drunk tank -
 * let the whole town know.
 */

extern int neg_lub_height;
extern int neg_lub_height_global;
extern int zero_detJ;
extern int zero_detJ_global;


/***************************************************************************/
/*                       std.h end                                         */
/***************************************************************************/
#endif

#if ( defined (linux) && defined (COMPILER_64BIT) )
/* This ifndef was extracted from _RF_GOMA_H ifndef because of goma failures
 * when compiled with 64-bit gnu compilers 4.1+. (It may be a trilinos 9
 * versus 6 issue.)  This is not needed if goma is compiled with autoconf.
 * (Eric Benner 8/6/09)
 */ 
#ifndef VERSION
#define VERSION "6.0.0"
#endif

#endif
