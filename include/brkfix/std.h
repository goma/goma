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

/* std.h -- include file with typical stuff needed for brk/fix routines
 *
 * Created: 1997/04/12 10:19 MDT pasacki@sandia.gov
 */

#ifndef _STD_H
#define _STD_H

#ifndef FILENAME_MAX_ACK
#define FILENAME_MAX_ACK		1024
#endif

#define LINE_BUFFER_LENGTH		1024

#define MAX_CHAR_ERR_MSG		1024
#define MAX_NEIGHBOR_NODES		1000
#define MAX_SYSTEM_COMMAND_LENGTH	1024
#define MAX_ADJOINING_SETS		10


/*
 * Needed for node-node comparisons of element faces in exo_conn.c
 */

#ifndef MAX_EPN
#define MAX_EPN			(50)
#endif


#ifndef MAX_NODES_PER_SIDE
#define MAX_NODES_PER_SIDE     (9)
#endif

/*
 * Element types for GOMA. These mongrels include some indications of
 * interpolation as well as basic element shape.
 */

typedef enum
{
  BILINEAR_QUAD,
  C_BILINEAR_QUAD, 
  S_BIQUAD_QUAD, 
  BIQUAD_QUAD, 
  P1_QUAD, 
  P0_QUAD, 
  TRILINEAR_HEX, 
  C_TRILINEAR_HEX, 
  S_TRIQUAD_HEX, 
  TRIQUAD_HEX,
  P1_HEX,
  P0_HEX,
  BILINEAR_SHELL,
  BQUAD_SHELL,
  UNDEFINED_ELEMENT_TYPE
} Element_type;

typedef enum
{
  LINE_SEGMENT,
  TRIANGLE,
  QUADRILATERAL,
  SHELL,
  TETRAHEDRON,
  PRISM,
  HEXAHEDRON,
  PYRAMID,
  UNDEFINED_ELEMENT_SHAPE
} Element_shape;

/*
 * Not what you might think! This is the maximum number of processors
 * that a given processor might communicate with. The overall number of
 * processors can be much larger.
 */

#define MAX_SEND_PROCS			32
#define MAX_RECV_PROCS			32

#define MAX_ELEMENT_PROCS		16 /* how many difft procs could
					    * assemble an element? Typically,
					    * 1, sometimes 2, rarely 3 or more
					    */
/* 
 * hold various ints for wr_graph_file.c
 */

#define LEN_DRAGON			8 

#define KILO				1024
#define MEGA				1048576
#define GIGA				1073741824
#define TERA				1099511627776
#define PETA				1125899906842624
#define EXA				1152921504606846976

#define HUGE_INT			1000000000

/*
 * Expect people to pick integer names like 0, 1, 2, .. and even -1, but
 * do not expect people to pick integer names like this one...
 * 
 * Machines aren't supposed to pick -1 for the set/proc name coming out of
 * the graph partitioner, so it will suffice to mean "unassigned" for us.
 */

#define UNDEFINED_SET_NAME -1
#define UNDEFINED_EQNVARID -55555

#define UNASSIGNED	-55555

/*
 * Initial sizes for provisional storage arrays collecting nodeset information
 * for each set/processor. By using realloc(), any hard limit to needed sizes
 * is not reached. Same story for sidesets.
 */

#define INIT_PROC_NS_DISTFACT_LIST_LENGTH	2048
#define INIT_PROC_NS_NODE_LIST_LENGTH		2048
#define INIT_PROC_NUM_NS			32

#define INIT_PROC_SS_DISTFACT_LIST_LENGTH	2048
#define INIT_PROC_SS_SIDE_LIST_LENGTH		2048
#define INIT_PROC_NUM_SS			32


#ifndef TRUE
#define TRUE    1
#define FALSE   0
#endif

#ifndef HAVE_PROTOTYPES
#if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#define  HAVE_PROTOTYPES
#endif
#endif

#undef PROTO
#ifdef HAVE_PROTOTYPES
#define PROTO(x) x
#else
#define PROTO(x) ()
#endif     

/*
 * sprintf() returns different types on different machines...
 */

#ifdef BSD
#undef SPRINTF_RETURNS_INT
#else
#define SPRINTF_RETURNS_INT
#endif

/* #if defined(_AIX) || defined (__hpux) || defined ()
 * #define SPRINTF_RETURNS_INT
 * #endif
 */

#ifdef SPRINTF_RETURNS_INT
typedef int	Spfrtn;
#endif

#ifndef SPRINTF_RETURNS_INT	/* Probably a ptr to char, then, eh? */
typedef char   *Spfrtn;
#endif

/*
 * strcpy() returns different types on different machines...
 */

#define STRCPY_RTN_IS_STRING

#ifdef STRCPY_RTN_IS_STRING
typedef char *Strcpy_rtn;
#endif
#ifdef STRCPY_RTN_IS_INT
typedef int  Strcpy_rtn;
#endif

#ifndef _DBL_TYPEDEF
#define _DBL_TYPEDEF
typedef double dbl;
#endif

#ifndef _FLT_TYPEDEF
#define _FLT_TYPEDEF
typedef float flt;
#endif

#define SZ_INT	(sizeof(int))
#define SZ_FLT	(sizeof(flt))
#define SZ_SHT	(sizeof(short))
#define SZ_CHR	(sizeof(char))
#define SZ_DBL	(sizeof(dbl))
#define SZ_LNG	(sizeof(long))
#define SZ_LLG	(sizeof(long long))
#define SZ_LDB	(sizeof(long double))
#define SZ_DDB	(sizeof(double double))

#define SZPINT	(sizeof(int *))
#define SZPFLT	(sizeof(flt *))
#define SZPSHT	(sizeof(short *))
#define SZPCHR	(sizeof(char *))
#define SZPDBL	(sizeof(dbl *))
#define SZPLNG	(sizeof(long *))
#define SZPLLG	(sizeof(long long *))
#define SZPLDB	(sizeof(long double *))
#define SZPDDB	(sizeof(double double *))

/*
 * Initialize Integer Vector -- fill an integer array with values
 */

#define INIT_IVEC(vec, val, len) \
{int i_; int *p_; for(i_=0,p_=(vec);i_<(len);i_++,p_++) *p_=(val);}

/*
 * Integer Swap -- exchange two integers
 */

#define ISWAP(i,j)	{int _itmp; _itmp=(i); (i)=(j); (j)=_itmp;}

/*
 * Build Unique List - if value is not in list of length listlen, then add it
 *                     to the end of the list and increment the listlen.
 */

#define BULL(value, list, listlen) \
{ int _tmp; _tmp = in_list((value), (list), (listlen)); if ( _tmp == -1 ) { (list)[(listlen)] = (value); (listlen)++;}}

/*
 * Boolean conveniences -- is or is not in ordered list.
 */

/* 
 * In Unordered List -- return TRUE or FALSE if integer value appears anywhere
 *                      in an integer array of length listlen.
 */

#define IUL(val, list, len) (in_list(val, list, len) != -1)

/* 
 * In Ordered List   -- return TRUE or FALSE if integer value *appears*
 *                      in a strictly monotone integer array of length len.
 */


#define IOL(val, list, len) (findex_mono(val, list, len) != -1)

/*
 * In Fenced List    -- return TRUE or FALSE if integer value is *bounded*
 *	                in a monotone integer array of length len.
 */

#define IFL(val, list, len) (fence_post(val, list, len) != -1)

/*
 * Comparison -- return the maximum or minimum of two values
 */

#ifndef MAX
#define MAX(x,y)	(( (x) > (y) ) ? (x) : (y)) 
#endif

#ifndef MIN
#define MIN(a,b)        (( (a) < (b) ) ? (a) : (b)) 
#endif

#endif
