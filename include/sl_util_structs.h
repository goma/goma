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
#define MPI			/* otherwise az_aztec.h trounces MPI_Request */
#endif
#endif

#include "sl_epetra_interface.h"
#include "az_aztec.h"

#ifdef COUPLED_FILL
/*
 * NUM_ALSS - Number of Aztec Linear Solver Systems. Here we have one.
 * Namely, JAC=0, is for the fully coupled Jacobian system.
 */

#define	JAC				0
#define NUM_ALSS			1
#else /* COUPLED_FILL */
/*
 * NUM_ALSS - Number of Aztec Linear Solver Systems. Here we have two.
 * First, JAC=0, is for the fully coupled Jacobian system. Second, FIL=1,
 * is for the smaller matrix resulting from a segregated VOF solution.
 */

#define	JAC				0
#define FIL				1
#define NUM_ALSS			2
#endif /* COUPLED_FILL */

struct Matrix_Data {
  struct Aztec_Linear_Solver_System *ams;
  double *x;                 /* Solution vector */
  double *x_old;             /* Solution vector , previous last time step */
  double *x_older;           /* Solution vector , previous prev time step */
  double *x_oldest;
  double *xdot;                      /* xdot of current solution                  */
  double *xdot_old;          /* xdot_old of current solution              */
  double *xdot_older;
  double *x_update;             /* last update vector */
  double *resid_vector;              /* Residual vector */
  double *scale;
};

struct Aztec_Linear_Solver_System
{
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

  int *bindx;		   /* "ija" or bindx */

  int *belfry;		   /* for storing ija[] extern pieces to hide
			    * them from the solver that only wants
			    * rows owned by this processor 
			    * use of this thing betrays an
			    * inherent inefficiency... */

  double *val;		   /* "a" */
  double *val_old;	   /* "a_old" */
  int     npn;             /* number of processor nodes, excluding external
			      nodes */
  int     npn_plus;        /* number of processor nodes, including external
			      nodes */
  int     npu;             /* number of processor dofs, excluding external
			      dofs */
  int     npu_plus;        /* number of processor dofs, including external
		              dofs */
  int     nnz;             /* length of "a" vector, excluding external dofs */
  int     nnz_plus;        /* length of "a" vector, including external dofs */
#ifdef MATRIX_DUMP  
  int     Number_Jac_Dump;
#endif

  C_Epetra_RowMatrix_t *RowMatrix; /* This is a Epetra_RowMatrix object */
  int *GlobalIDs;                  /* Pointer to global ids of DOFs (only available with epetra) */
};

/* See paper by Hood (1976, J. Numer. Meth Engr.) for details of this
 * struct
 */
struct Frontal_Solver_System  
{
  double *bc;              /* for each nodal dof this is assigned zero
		            * unless a boundary condition is applied, in
			    * which case it gets the bc value. */
  int    *ncn;             /* The total number of dofs in each element */
  int    *ncn_base;        /* The total number of dofs in each element */
  int    *ncod;            /* The entries are coded zero for a nodal dof with
			    * no applied bc, and unity for an applied bc */
  int   *nop;              /* The traditional nodal connectivity array,
		            * NE X NBN */
  int   *nopdof;           /* The unknown connectivity array, NE X NCN */
  int   *nopp;             /* Coded value of first dof at each node */
  int   *mdf;              /* Number of dofs at each node */
  int   ntra;              /* 1 on first entry into front, zero on
			    * subsequent entries */
  int   *el_proc_assign;   /* processor assignment of element i.  */
  int   *constraint;       /* true of false array for constraint equations */
  int   *level;            /* dissection level for elimination (see
			    * mpfront doc  */
  int   *perm;             /* rcm aztec reordering list */
  int   *iperm;            /* rcm aztec reordering list, inverse */
  int   *mask;             /* something needed for rcm-aztec */
};

extern struct Frontal_Solver_System *fss;



#endif
