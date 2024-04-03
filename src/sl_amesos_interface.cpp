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

#include <Amesos_config.h>
#include <stdlib.h>
#include <string>

#include "Amesos_BaseSolver.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_DataAccess.h"
#include "Epetra_RowMatrix.h"
#include "az_aztec.h"
#include "rf_fem_const.h"
#ifndef GOMA_SL_AMESOS_INTERFACE_CC
#define GOMA_SL_AMESOS_INTERFACE_CC
#endif

/* This removes the entire file if Amesos & Trilinos are not defined */
#if defined(GOMA_ENABLE_AMESOS) && defined(TRILINOS)

#if defined(PARALLEL) && !defined(EPETRA_MPI)
#define EPETRA_MPI
#endif

#ifndef __cplusplus
#define __cplusplus
#endif

#include <iostream>

#include "mpi.h"
#include "sl_amesos_interface.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#if 1
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorOut.h"
#endif

#include "Amesos.h"
#include "Amesos_ConfigDefs.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Trilinos_Util.h"
#include "linalg/sparse_matrix.h"
#include "linalg/sparse_matrix_epetra.h"
#include "sl_util_structs.h"

static void GomaMsr2EpetraCsr(struct GomaLinearSolverData *ams, Epetra_CrsMatrix *A, int newmatrix);
void amesos_solve(char *choice,
                  struct GomaLinearSolverData *ams,
                  double *x_,
                  double *b_,
                  int NewMatrix,
                  int imtrx) {

  /* Initialize MPI communications */
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  /* Define internal variables */
  static std::string Pkg_Name;
  static Epetra_CrsMatrix *A[MAX_NUM_MATRICES]{nullptr};
  static Epetra_LinearProblem Problem[MAX_NUM_MATRICES];
  static Amesos_BaseSolver *A_Base[MAX_NUM_MATRICES] = {nullptr};
  Amesos A_Factory;

  /* Convert to Epetra format */
  if (ams->GomaMatrixData == NULL) {
    if (!ams->solveSetup) {
      if (A[imtrx] != nullptr)
        delete A[imtrx];
      A[imtrx] = (Epetra_CrsMatrix *)construct_Epetra_CrsMatrix(ams);
    }
    GomaMsr2EpetraCsr(ams, A[imtrx], !ams->solveSetup);
  } else {
    GomaSparseMatrix matrix = (GomaSparseMatrix)ams->GomaMatrixData;
    EpetraSparseMatrix *epetra_matrix = static_cast<EpetraSparseMatrix *>(matrix->data);
    A[imtrx] = dynamic_cast<Epetra_CrsMatrix *>(epetra_matrix->matrix.get());
  }
  const Epetra_Map &map = A[imtrx]->RowMatrixRowMap();
  Epetra_Vector x(Copy, map, x_);
  Epetra_Vector b(Copy, map, b_);

#if 0
  EpetraExt::RowMatrixToMatrixMarketFile("Jep.mm", *A[imtrx]);
  EpetraExt::VectorToMatrixMarketFile("xep.mm", x);
  EpetraExt::VectorToMatrixMarketFile("bep.mm", b);
  MPI_Finalize();
  exit(0);
#endif
  /* Choose correct solver */
  /* Assemble linear problem */
  Problem[imtrx].SetOperator(A[imtrx]);
  Problem[imtrx].SetLHS(&x);
  Problem[imtrx].SetRHS(&b);
  Problem[imtrx].CheckInput();
  if (!ams->solveSetup) {
    std::string Pkg_Choice = choice;
    if (Pkg_Choice == "KLU")
      Pkg_Name = "Amesos_Klu";
#ifdef HAVE_AMESOS_UMFPACK
    else if (Pkg_Choice == "UMF")
      Pkg_Name = "Amesos_Umfpack";
#endif
    else if (Pkg_Choice == "LAPACK")
      Pkg_Name = "Amesos_Lapack";
#if defined(HAVE_AMESOS_SUPERLUDIST)
    else if (Pkg_Choice == "SUPERLU")
      Pkg_Name = "Amesos_Superludist";
    else if (Pkg_Choice == "SUPERLU_PARALLEL")
      Pkg_Name = "Amesos_Superludist";
#elif defined(HAVE_AMESOS_SUPERLU)
    else if (Pkg_Choice == "SUPERLU")
      Pkg_Name = "Amesos_Superlu";
#endif
#ifdef HAVE_AMESOS_SCALAPACK
    else if (Pkg_Choice == "SCALAPACK")
      Pkg_Name = "Amesos_Scalapack";
#endif
#ifdef HAVE_AMESOS_MUMPS
    else if (Pkg_Choice == "MUMPS")
      Pkg_Name = "Amesos_Mumps";
#endif
    else {
      std::cerr << "Error: Unsupported Amesos solver package" << std::endl;
      exit(-1);
    }

    /* Create Amesos base package */
    if (A_Base[imtrx] != nullptr)
      delete A_Base[imtrx];

    A_Base[imtrx] = A_Factory.Create(Pkg_Name.c_str(), Problem[imtrx]);
    if (A_Base[imtrx] == 0) {
      std::cout << "Error in amesos_solve_msr" << std::endl;
      std::cout << "It is likely that the solver package: " << Pkg_Name
                << " has not been linked into the amesos library " << std::endl;
      exit(-1);
    }
  }

  /* Solve problem */
  if (!ams->solveSetup) {
    A_Base[imtrx]->SymbolicFactorization();
  }
  A_Base[imtrx]->NumericFactorization();
  A_Base[imtrx]->Solve();

  /* Convert solution vector */
  int NumMyRows = map.NumMyElements();
  for (int i = 0; i < NumMyRows; i++) {
    x_[i] = x[i];
  }

  /* Cleanup problem */
  ams->solveSetup = 1;
}

static void GomaMsr2EpetraCsr(struct GomaLinearSolverData *ams, Epetra_CrsMatrix *A, int newmatrix)

{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int *bindx = ams->bindx;
  double *val = ams->val;

  int NumMyRows = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border];
  int NumExternal = ams->data_org[AZ_N_external];
  int NumMyCols = NumMyRows + NumExternal;

  (*A).PutScalar(0.0);

  const Epetra_Map &RowMap = (*A).RowMatrixRowMap();

  int *MyGlobalElements = RowMap.MyGlobalElements();

  double *dblColGIDs = new double[NumMyCols];
  int *ColGIDs = new int[NumMyCols];

  for (int i = 0; i < NumMyRows; i++)
    dblColGIDs[i] = (double)MyGlobalElements[i];

  AZ_exchange_bdry(dblColGIDs, ams->data_org, ams->proc_config);

  {
    for (int j = 0; j < NumMyCols; j++)
      ColGIDs[j] = (int)dblColGIDs[j];
  }

  int MaxNNZ = 0;
  for (int i = 0; i < NumMyRows; i++) {
    int NumNz = bindx[i + 1] - bindx[i] + 1;

    if (NumNz > MaxNNZ)
      MaxNNZ = NumNz;
  }

  int *Indices = new int[MaxNNZ];
  double *Values;

  for (int i = 0; i < NumMyRows; i++) {
    int NumNz = bindx[i + 1] - bindx[i];

    int *k = bindx + bindx[i];
    for (int j = 0; j < NumNz; j++) {
      Indices[j] = ColGIDs[*k];
      k++;
    }

    Values = val + bindx[i];

    if (newmatrix) {
      (*A).InsertGlobalValues(MyGlobalElements[i], NumNz, Values, Indices);
      (*A).InsertGlobalValues(MyGlobalElements[i], 1, &(val[i]), MyGlobalElements + i);
    } else {
      (*A).SumIntoGlobalValues(MyGlobalElements[i], NumNz, Values, Indices);
      (*A).SumIntoGlobalValues(MyGlobalElements[i], 1, &(val[i]), MyGlobalElements + i);
    }
  }

  (*A).FillComplete();

  delete[] dblColGIDs;
  delete[] ColGIDs;
  delete[] Indices;

  return;
}

void *construct_Epetra_CrsMatrix(struct GomaLinearSolverData *ams) {
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // if ( ams->EpetraMat != NULL ) return ( ams->EpetraMat );

  int *bindx = ams->bindx;

  int NumMyRows = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border];

  Epetra_Map RowMap(-1, NumMyRows, 0, comm);

  int *NumNz = new int[NumMyRows];

  for (int i = 0; i < NumMyRows; i++) {
    NumNz[i] = bindx[i + 1] - bindx[i] + 1;
  }

  Epetra_CrsMatrix *A = new Epetra_CrsMatrix(Copy, RowMap, NumNz);

  delete[] NumNz;

  return (A);
}

#endif // if defined(GOMA_ENABLE_AMESOS) && defined(TRILINOS)

#if 0

#include <iostream.h>
#include <stdio.h>
#include <string.h>

#include "sl_amesos_interface.h"
#include "sl_util_structs.h"


void
amesos_solve_msr(char *choice,  
		 struct Aztec_Linear_Solver_System *ams,
		 double *x_, 
		 double *b_,
		 int NewMatrix )
{
	
	
  std::string err_msg("Error: Need to compile with GOMA_ENABLE_AMESOS flag before using AMESOS solver packages.");
  std::string err_msg2("   Also make sure appropriate libraries are linked for solver packages.");
	
  cout << err_msg << endl;
  cout << err_msg2 << endl;
	
  exit(-1);
	
  return;
}

/* End of if statement for Amesos and Trilinos */
#endif
