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
 *
 *
 * This conveniently encapsulates all the data necessary to describe a
 * matrix system to the Aztec solver package. Since Aztec is equipped to
 * handle matrices in the DMSR and DVBR formats, there are a number of
 * index arrays to keep track of. Here, just stick these into a data structure.
 *
 * Created: 1997/02/20 07:40 MST pasacki@sandia.gov
 * Modified: 1998/12/24 (why am I not home spending time with my family? GEEK!)  prschun@sandia.gov
 *
 * Added structure to handle transition to a frontal solver.  Memory is allocated
 * only if the frontal solver is requested.
 *
 * Modified:
 */

#ifndef GOMA_SL_UTIL_STRUCTS_H
#define GOMA_SL_UTIL_STRUCTS_H

#include <stdio.h>

#ifdef PARALLEL
#ifndef MPI
#define MPI /* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "az_aztec.h"

/*
 * NUM_ALSS - Number of Aztec Linear Solver Systems. Here we have one.
 * Namely, JAC=0, is for the fully coupled Jacobian system.
 */

#define JAC      0
#define NUM_ALSS 1

struct GomaLinearSolverData {
  int proc_config[AZ_PROC_SIZE];

  int options[AZ_OPTIONS_SIZE];

  double params[AZ_PARAMS_SIZE];

  int *data_org;

  double status[AZ_STATUS_SIZE];

  int N;
  int N_update;
  int N_external;

  /*
   * These 4 mapping arrays are only used if you have chosen
   * a non canonical ordering for your unknowns on a given
   * processor.
   */

  int *update;
  int *update_index;
  int *external;
  int *extern_index;

  int *indx;
  int *rpntr;
  int *cpntr;
  int *bpntr;

  int mat_type;

  int *bindx; /* "ija" or bindx */

  int *belfry; /* for storing ija[] extern pieces to hide
                * them from the solver that only wants
                * rows owned by this processor
                * use of this thing betrays an
                * inherent inefficiency... */

  double *val;     /* "a" */
  double *val_old; /* "a_old" */
  int npn;         /* number of processor nodes, excluding external
                      nodes */
  int npn_plus;    /* number of processor nodes, including external
                      nodes */
  int npu;         /* number of processor dofs, excluding external
                      dofs */
  int npu_plus;    /* number of processor dofs, including external
                      dofs */
  int nnz;         /* length of "a" vector, excluding external dofs */
  int nnz_plus;    /* length of "a" vector, including external dofs */
#ifdef MATRIX_DUMP
  int Number_Jac_Dump;
#endif

  int solveSetup;

  void *PetscMatrixData;
  void *GomaMatrixData;
  void *SolverData;
  void (*DestroySolverData)(struct GomaLinearSolverData *ams);
};

struct Matrix_Data {
  struct GomaLinearSolverData *ams;
  double *x;       /* Solution vector */
  double *x_prev;  /* Solution vector */
  double *x_old;   /* Solution vector , previous last time step */
  double *x_older; /* Solution vector , previous prev time step */
  double *x_oldest;
  double *xdot;     /* xdot of current solution                  */
  double *xdot_old; /* xdot_old of current solution              */
  double *xdot_older;
  double *x_update;     /* last update vector */
  double *resid_vector; /* Residual vector */
  double *scale;
};
#endif
