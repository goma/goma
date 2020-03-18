/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.	  *
\************************************************************************/
#if defined(PARALLEL) && !defined(EPETRA_MPI)
#define EPETRA_MPI
#endif

#ifndef __cplusplus
#define __cplusplus
#endif

#include "mpi.h"
#include "AztecOO.h"
#include "Epetra_DataAccess.h"
#include "Epetra_RowMatrix.h"
#include "az_aztec.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else 
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "sl_util_structs.h"
#include "sl_aztecoo_interface.h"

extern "C" {

/**
 * Solve using the AztecOO solver with given aztec options
 *
 * A x = b
 *
 * @param ams ams with RowMatrix to be solved
 * @param x_ solution vector (initial guess), solution placed in this vector
 * @param b_ residual vector
 */
void
aztecoo_solve_epetra(struct Aztec_Linear_Solver_System *ams,
		  double *x_, 
		  double *b_) {

  /* Initialize MPI communications */
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  /* Define internal variables */
  Epetra_RowMatrix *A = ams->RowMatrix;
  Epetra_LinearProblem Problem;
    
  Epetra_Map map = A->RowMatrixRowMap();
  Epetra_Vector x(Copy, map, x_);
  Epetra_Vector b(Copy, map, b_);
  
  /* Assemble linear problem */
  Problem.SetOperator(A);
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);
  Problem.CheckInput();

  /* Create Amesos base package */
  AztecOO solver(Problem);

  solver.SetAllAztecOptions(ams->options);
  solver.SetAllAztecParams(ams->params);

  double tolerance = ams->params[AZ_tol];
  int max_iterations = ams->options[AZ_max_iter];
  /* Solve problem */
  solver.Iterate(max_iterations,tolerance);
  solver.GetAllAztecStatus(ams->status);

  /* Convert solution vector */
  int NumMyRows = map.NumMyElements();
  for(int i=0; i< NumMyRows; i++) {
    x_[i] = x[i]; 
  }
}

} /* extern "C" */
