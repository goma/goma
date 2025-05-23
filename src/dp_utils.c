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

/* dp_utils.c -- distributed processing utility routines
 *
 * These routines are meant as a convenience for using MPI with goma. We
 * make extensive use of MPI derived datatypes in goma. To do this easily,
 * the following data structures and routines have been created.
 *
 * Created: 1997/06/17 06:51 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dp_types.h"
#include "dp_utils.h"
#include "mm_eh.h"
#include "mpi.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "std.h"

#define GOMA_DP_UTILS_C

static char mpistringbuffer[80];

static Spfrtn sr;

#include "az_aztec.h"

int Proc_Config[AZ_PROC_SIZE];

/************************************************************************/
/************************************************************************/
/************************************************************************/

/************************************************************************/
/************************************************************************/
/************************************************************************/

DDD ddd_alloc(void) {
  int n = SZ_DDD_INIT_MAX_MEMBERS;
  DDD p;
  p = (DDD)calloc(1, sizeof(struct Derived_Datatype_Description));
  p->num_members = 0;
  p->max_members = n;
  p->block_count = (int *)calloc(n, sizeof(int));
  p->data_type = (MPI_Datatype *)calloc(n, sizeof(MPI_Datatype));
  p->address = (MPI_Aint *)calloc(n, sizeof(MPI_Aint));
  return (p);
}

void ddd_free(DDD p) {
  MPI_Type_free(&p->new_type);
  free(p->block_count);
  free(p->data_type);
  free(p->address);
  free(p);
  return;
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#ifdef PARALLEL
static int ddd_internal_count = 0;
#endif

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void ddd_add_member2(void *address, int blockcount, size_t byte_size)

/*****************************************************************
 *
 * ddd_add_member2:
 *
 *   Buffered broadcast routine. This routine will force an
 *   immediate exchange when the buffer gets larger than a fixed
 *   size (currently 8000 bytes). At this point, its worthwhile
 *   to exchange immediately, because further buffering won't yield
 *   speed improvements. Also, memory requirements for reallocs()
 *   in AZ_broadcast() are kept to a minimum.
 *
 *   This routine uses two static variables (listed below)
 *   and must be used in combination with the routine
 *   ddd_set_commit2() listed below. ddd_set_commit2() does
 *   the final broadcast, if needed.
 *
 *   Inbetween the calls to ddd_add_commit2() and
 *   ddd_set_commit2(), no other mp communication should be
 *   carried out that involves syncing processors.
 *   This is due to the fact that the send is carried out in
 *   Proc 0 in the ddd_set_commit2() routine, while the
 *   read is carried out on all other processors in the
 *   ddd_add_commit2() routine.
 *
 *   Note: Errors in the broadcast lead to immediate error exits.
 *
 *  Input
 * -------
 *   address    : base address of the item to be broadcast
 *   blockcount : Number of items in the object to be sent
 *   byte_size  : Number of bytes per item in the object to be
 *                sent.
 *
 *  Static Variables
 * ------------------
 *   dd_internal_count : Internal count of the buffer length
 *   Proc_Config[3] : MP configuration information, in a form
 *                    needed by Aztec.
 *****************************************************************/
{
#ifdef PARALLEL
  int length;

  if (byte_size <= 0) {
    fprintf(stderr, " ddd_add_member2 ERROR: byte_size = %ld\n", (long int)byte_size);
    GOMA_EH(GOMA_ERROR, "ddd_add_member2 parameter error");
  }
  length = blockcount * byte_size;
  if (length > 0) {
    AZ_broadcast(address, length, Proc_Config, AZ_PACK);
    ddd_internal_count += length;
    if (ddd_internal_count > 8000) {
      AZ_broadcast(NULL, 0, Proc_Config, AZ_SEND);
      ddd_internal_count = 0;
    }
  }
#endif
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void ddd_set_commit2(void)

/******************************************************************
 *
 * ddd_set_member2:
 *
 *   See description for ddd_add_member2 above().
 ******************************************************************/
{
#ifdef PARALLEL
  if (ddd_internal_count > 0) {
    AZ_broadcast(NULL, 0, Proc_Config, AZ_SEND);
    ddd_internal_count = 0;
  }
#endif
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void ddd_add_member(DDD p, void *address, int blockcount, MPI_Datatype type) {
  char err_msg[MAX_CHAR_IN_INPUT];
  int i;
  int n;

  i = p->num_members;

  if (blockcount == 0) {
    sr = sprintf(err_msg, "Attempt to add member %d type %s w/ 0 length!", i, type2string(type));
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  if (address == NULL) {
    sprintf(err_msg, "attempt to add member %d type %s blockcount %d with a NULL address!", i,
            type2string(type), blockcount);
    GOMA_EH(GOMA_ERROR, err_msg);
  }

  if (i > p->max_members - 1) {
    p->max_members += SZ_DDD_INIT_MAX_MEMBERS;
    n = p->max_members;

    p->block_count = (int *)realloc(p->block_count, n * sizeof(int));
    p->data_type = (MPI_Datatype *)realloc(p->data_type, n * sizeof(MPI_Datatype));
    p->address = (MPI_Aint *)realloc(p->address, n * sizeof(MPI_Aint));
  }

  p->block_count[i] = blockcount;
  p->data_type[i] = type;
#ifdef PARALLEL
  /*
   *  This is just the identity operator on most systems
   */
  MPI_Get_address(address, &p->address[i]);

#endif
  p->num_members++;

  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void ddd_set_commit(DDD p) {
#ifdef PARALLEL
  MPI_Type_create_struct(p->num_members, p->block_count, p->address, p->data_type, &p->new_type);
  MPI_Type_commit(&p->new_type);
  MPI_Type_get_extent(p->new_type, &p->lb, &p->extent);
  MPI_Type_size(p->new_type, &p->size);
  /*  rtn = MPI_Type_count(p->new_type, &p->count); */
#endif

  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

char *type2string(MPI_Datatype type) {
  if (type == MPI_INT) {
    strcpy(mpistringbuffer, "MPI_INT");
  } else if (type == MPI_CHAR) {
    strcpy(mpistringbuffer, "MPI_CHAR");
  } else if (type == MPI_FLOAT) {
    strcpy(mpistringbuffer, "MPI_FLOAT");
  } else if (type == MPI_DOUBLE) {
    strcpy(mpistringbuffer, "MPI_DOUBLE");
  } else {
    strcpy(mpistringbuffer, "MPI_WHAT?");
  }

  return (mpistringbuffer);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int ProcWithMaxInt(const int value, int *maxval)

/********************************************************************
 *
 * ProcWithMaxInt():
 *
 *   This routine will find the maximum integer value input from
 *   all processors. In case of ties, it will return the lowest
 *   numbered processor ID with the max value.
 *
 *  Output
 *  -------
 *  *maxval = Maximum int found on all processors
 *
 *  Return
 *  -------
 *  Routine returns the processor ID that had the max value. In case
 *  of ties, it will return the lowest numbered Processer ID.
 ********************************************************************/
{
  int in_buf[2], out_buf[2];
#ifdef PARALLEL
  in_buf[0] = value;
  in_buf[1] = ProcID;
  MPI_Allreduce((void *)in_buf, (void *)out_buf, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
  *maxval = out_buf[0];
  return out_buf[1];
#else
  *maxval = value;
  return 0;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void check_parallel_error_FL(char *errstring, char *file_name, int lineno)

/********************************************************************
 *
 * check_parallel_error_FL
 * check_parallel_error(string):
 *
 *   Checks the global variable parallel_error to see if a
 *   parallel error condition has been posted. If it has then
 *   this routine posts an error message and terminates GOMA.
 ********************************************************************/
{
#ifdef PARALLEL

  MPI_Allreduce(&parallel_err, &parallel_err_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  if (parallel_err_global) {
    DPRINTF(stderr, "PARALLEL ERROR: ");
    if (errstring)
      DPRINTF(stderr, "%s", errstring);
    DPRINTF(stderr, " check previous output for cause of error: %s line %d", file_name, lineno);
    fflush(stderr);
    if (parallel_err) {
      fprintf(stderr, "\tP_%d had a flagged error\n", ProcID);
    }
    (void)MPI_Finalize();
    exit(-1);
  }
#else
  if (parallel_err || parallel_err_global) {
    DPRINTF(stderr, "SERIAL PROG ERROR: ");
    DPRINTF(stderr, " check previous output for cause of error: %s line %d", file_name, lineno);
    fflush(stderr);
    if (parallel_err) {
      fprintf(stderr, "\tP_%d had a flagged error\n", ProcID);
    }
    exit(-1);
  }
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void ReduceBcast_BOR(int *ivec, int length)

/********************************************************************
 *
 * ReduceBcast_BOR()
 *
 *  Does a logical bitwize OR operation on a vector of ints
 *  distibuted across different processors.
 ********************************************************************/
{
#ifdef PARALLEL
  int k, err, *ivec_recv;
  if (length <= 0) {
    printf(" ReduceBcast_BOR Warning, length = %d\n", length);
    return;
  }
  if (ivec == NULL) {
    GOMA_EH(GOMA_ERROR, " ReduceBcast_BOR fatal Interface error");
  }
  ivec_recv = alloc_int_1(length, INT_NOINIT);
  err = MPI_Reduce(ivec, ivec_recv, length, MPI_INT, MPI_BOR, 0, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, " ReduceBcast_BOR fatal MPI Error");
  }
  err = MPI_Bcast(ivec_recv, V_LAST, MPI_INT, 0, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, " ReduceBcast_BOR fatal MPI Error");
  }
  for (k = 0; k < V_LAST; k++) {
    ivec[k] = ivec_recv[k];
  }
  safer_free((void **)&ivec_recv);
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int gmaxloc_int(const int value, const int local_loc, int *global_maxloc)

/********************************************************************
 *
 * gmaxloc_int():
 *
 *   This routine will find the maximum integer value input from
 *   all processors. It will also return the location from which
 *   the maximum occurred. In case of ties, it will return the
 *   location from the lowest numbered processor ID with the
 *   max value.
 *
 *  Output
 *  -------
 *  *global_maxloc = Returns the value of local_loc from the
 *                   processor containing the max value.
 *
 *  Return
 *  -------
 *  Routine returns the maximum value of "value" input from
 *  any of the processors.
 ********************************************************************/
{
#ifdef PARALLEL
  int in_buf[2], out_buf[2], err;
  in_buf[0] = value;
  in_buf[1] = local_loc;
  err = MPI_Allreduce((void *)in_buf, (void *)out_buf, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gmaxloc_int: MPI_Allreduce returned an error");
  }
  *global_maxloc = out_buf[1];
  return out_buf[0];
#else
  *global_maxloc = 0;
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int gminloc_int(const int value, const int local_loc, int *global_minloc)

/********************************************************************
 *
 * gminloc_int():
 *
 *   This routine will find the minimum integer value input from
 *   all processors. It will also return the location from which
 *   the minimum occurred. In case of ties, it will return the
 *   location from the lowest numbered processor ID with the
 *   min value.
 *
 *  Output
 *  -------
 *  *global_minloc = Returns the value of local_loc from the
 *                   processor containing the min value.
 *
 *  Return
 *  -------
 *  Routine returns the minimum value of "value" input from
 *  any of the processors.
 ********************************************************************/
{
#ifdef PARALLEL
  int in_buf[2], out_buf[2], err;
  in_buf[0] = value;
  in_buf[1] = local_loc;
  err = MPI_Allreduce((void *)in_buf, (void *)out_buf, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gminloc_int: MPI_Allreduce returned an error");
  }
  *global_minloc = out_buf[1];
  return out_buf[0];
#else
  *global_minloc = 0;
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int gmin_int(int value)

/********************************************************************
 *
 * gmin_int
 *
 *   This routine will return the minimum of a single integer
 *   distributed across the processors.
 *
 *
 *  Return
 *  -------
 *  Routine returns the minimum value.
 ********************************************************************/
{
#ifdef PARALLEL
  int out_buf, err;
  err = MPI_Allreduce((void *)&value, (void *)&out_buf, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gmin_int: MPI_Allreduce returned an error");
  }
  return out_buf;
#else
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int gmax_int(int value)

/********************************************************************
 *
 * gmax_int
 *
 *   This routine will return the maximum of a single integer
 *   distributed across the processors.
 *
 *
 *  Return
 *  -------
 *  Routine returns the maximum value.
 ********************************************************************/
{
#ifdef PARALLEL
  int out_buf, err;
  err = MPI_Allreduce((void *)&value, (void *)&out_buf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gmax_int: MPI_Allreduce returned an error");
  }
  return out_buf;
#else
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int gsum_Int(int value)

/********************************************************************
 *
 * gsum_Int
 *
 *   This routine will return the sum of a single integer
 *   distributed across the processors.
 *
 *
 *  Return
 *  -------
 *  Routine returns the sum.
 ********************************************************************/
{
#ifdef PARALLEL
  int out_buf, err;
  err = MPI_Allreduce((void *)&value, (void *)&out_buf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gsum_Int: MPI_Allreduce returned an error");
  }
  return out_buf;
#else
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

double gavg_double(double value)

/********************************************************************
 *
 * gavg_double
 *
 *   This routine will find the average value of an input
 *   across all of the processors
 *
 *
 *  Return
 *  -------
 *  Routine returns the average value
 ********************************************************************/
{
#ifdef PARALLEL
  int err;
  double out_buf;
  err = MPI_Allreduce((void *)&value, (void *)&out_buf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    GOMA_EH(GOMA_ERROR, "gavg_double: MPI_Allreduce returned an error");
  }
  return out_buf / Num_Proc;
#else
  return value;
#endif
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void print_sync_start(int do_print)

/*********************************************************************
 *
 * print_sync_start:
 *
 *  Wrapper around AZ_print_sync_start()
 *********************************************************************/
{
  if (Num_Proc > 1) {
    AZ_print_sync_start(ProcID, do_print, Proc_Config);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void print_sync_end(int do_print)

/*********************************************************************
 *
 * print_sync_end:
 *
 *  Wrapper around AZ_print_sync_end()
 *********************************************************************/
{
  if (Num_Proc > 1) {
    AZ_print_sync_end(Proc_Config, do_print);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void sync_processors(void)

/*********************************************************************
 *
 * sync_processors:
 *
 * This function will synchronize all processors
 * it is a wrapper around AZ_sync(). Also, since sync is primarily
 * used for printing applications, we flush both the standard out and
 * standard error here as well.
 *********************************************************************/
{
  fflush(stdout);
  fflush(stderr);
  if (Num_Proc > 1) {
    AZ_sync(Proc_Config);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
