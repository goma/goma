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
 * Created: 1997/06/17 07:07 MDT pasacki@sandia.gov
 *
 * Revised:1999/04/27 14:10 MDT pasacki@sandia.gov
 */

#ifndef GOMA_DP_UTILS_H
#define GOMA_DP_UTILS_H

#include <stddef.h>

#include "az_aztec.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "mpi.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_DP_UTILS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_DP_UTILS_C
#define EXTERN extern
#endif

EXTERN DDD ddd_alloc
(void);

EXTERN void ddd_add_member
(DDD,			/* which collection? */
       void *,			/* ptr new datachunk */
       int,			/* how many of type? */
       MPI_Datatype);		/* what type? */

EXTERN void ddd_add_member2
(void *,			/* address */
       int ,			/* blockcount */
       size_t );		/* byte_size*/

EXTERN void ddd_set_commit2
(void);		

EXTERN void ddd_set_commit 
(DDD);

EXTERN void ddd_free 
(DDD);

EXTERN char *type2string
(MPI_Datatype );		/* type - MPI data type MPI_INT, etc */

extern int ProcWithMaxInt(const int, int *);
#define check_parallel_error(arg1) \
	check_parallel_error_FL((arg1), __FILE__, __LINE__)
extern void check_parallel_error_FL(char *, char *, int);
extern void ReduceBcast_BOR(int *, int);
extern int  gmaxloc_int(const int, const int, int *);
extern int  gminloc_int(const int, const int, int *);
extern int  gmax_int(int);
extern int  gmin_int(int);
extern int  gsum_Int(int);
extern double gavg_double(double);
extern void print_sync_start(int);
extern void print_sync_end(int);
extern void sync_processors(void);

#ifdef PARALLEL
extern int Proc_Config[AZ_PROC_SIZE];
#else
#ifdef TRILINOS
extern int Proc_Config[AZ_PROC_SIZE];
#endif
#endif

#endif
