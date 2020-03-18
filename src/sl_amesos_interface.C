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

#include <stdlib.h>
#include <string>

#include "Amesos_BaseSolver.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_DataAccess.h"
#include "Epetra_RowMatrix.h"
#include "az_aztec.h"
#ifndef GOMA_SL_AMESOS_INTERFACE_CC
#define GOMA_SL_AMESOS_INTERFACE_CC
#endif

/* This removes the entire file if Amesos & Trilinos are not defined */
#if defined(ENABLE_AMESOS) && defined(TRILINOS)

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

#include "Trilinos_Util.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Amesos.h"
#include "Amesos_ConfigDefs.h"
#include "sl_util_structs.h"

static void GomaMsr2EpetraCsr ( struct Aztec_Linear_Solver_System *,
				Epetra_CrsMatrix *);

void
amesos_solve_msr( char *choice,
		  struct Aztec_Linear_Solver_System *ams,
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

  static int prev_matrix = 0;

  /* Define internal variables */
  static int FirstRun = 1;
  static std::string Pkg_Name;
  static Epetra_CrsMatrix *A;
  Epetra_LinearProblem Problem;
  static Amesos_BaseSolver *A_Base;
  static Amesos A_Factory;

  if (prev_matrix != imtrx) {
    if (!FirstRun) {
      delete A;
    }
    FirstRun = 1;
    prev_matrix = imtrx;
  }
    
  /* Convert to Epetra format */
  if (NewMatrix) {
    if (!FirstRun) delete A;
    A = (Epetra_CrsMatrix * ) construct_Epetra_CrsMatrix ( ams ) ;  
    GomaMsr2EpetraCsr( ams, A);
  }
  const Epetra_Map &map = (*A).RowMatrixRowMap();
  Epetra_Vector x(Copy, map, x_);
  Epetra_Vector b(Copy, map, b_);
  
  /* Choose correct solver */
  if (FirstRun) {
    std::string Pkg_Choice  = choice;
    if( Pkg_Choice == "UMF") 
      Pkg_Name = "Amesos_Umfpack";
    else if ( Pkg_Choice == "KLU") 
      Pkg_Name = "Amesos_Klu";
    else if ( Pkg_Choice == "LAPACK") 
      Pkg_Name = "Amesos_Lapack";
    else if ( Pkg_Choice == "SUPERLU") 
      Pkg_Name = "Amesos_Superludist";
    else if ( Pkg_Choice == "SUPERLU_PARALLEL") 
      Pkg_Name = "Amesos_Superludist";
    else if ( Pkg_Choice == "SCALAPACK") 
      Pkg_Name = "Amesos_Scalapack";
    else if ( Pkg_Choice == "MUMPS") 
      Pkg_Name = "Amesos_Mumps";
    else {
      std::cout << "Error: Unsupport Amesos solver package"<<std::endl ;
      exit(-1);
    }
  }

  /* Assemble linear problem */
  Problem.SetOperator(A);
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);
  Problem.CheckInput();

  /* Create Amesos base package */

  A_Base = A_Factory.Create( Pkg_Name.c_str(), Problem );
  if( A_Base == 0 ) {
    std::cout << "Error in amesos_solve_msr" <<std::endl;
    std::cout << "It is likely that the solver package: " << Pkg_Name << " has not been linked into the amesos library " << std::endl;
    exit (-1 ) ;
  }

  /* Solve problem */
  A_Base->SymbolicFactorization();
  A_Base->NumericFactorization();
  A_Base->Solve();

  /* Convert solution vector */
  int NumMyRows = map.NumMyElements();
  for(int i=0; i< NumMyRows; i++) {
    x_[i] = x[i]; 
  }
  
  /* Cleanup problem */
  FirstRun = 0;
  delete A_Base;
}

/**
 * Solve linear system using amesos using epetra  (C interface)
 * A x = b
 * @param choice Amesos solver choice
 * @param ams ams containing Row Matrix for A
 * @param x_ array of values for solution to be put in (containing initial guess)
 * @param resid_vector array representing residual vector (b)
 * @return 0 on success
 */
int amesos_solve_epetra( char *choice,
			 struct Aztec_Linear_Solver_System *ams,
			 double *x_,
			 double *resid_vector,
			 int imtrx ) {

  /* Initialize MPI communications */
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  static int prev_matrix = 0;
  Epetra_RowMatrix *A = ams->RowMatrix;
  static Epetra_LinearProblem Problem;
  static Amesos_BaseSolver *Solver;
  static Amesos A_Factory;
  std::string Pkg_Name;
  static bool firstSolve = true;

  if (prev_matrix != imtrx) {
    if (!firstSolve) {
      delete Solver;
    }
    firstSolve = true;
  }

  const Epetra_Map &map = (*A).RowMatrixRowMap();

  Epetra_Vector x(Copy, map, x_);
  Epetra_Vector b(Copy, map, resid_vector);

  std::string Pkg_Choice  = choice;

  if (Pkg_Choice == "UMF")
    Pkg_Name = "Amesos_Umfpack";
  else if (Pkg_Choice == "KLU")
    Pkg_Name = "Amesos_Klu";
  else if (Pkg_Choice == "LAPACK")
    Pkg_Name = "Amesos_Lapack";
  else if (Pkg_Choice == "SUPERLU")
    Pkg_Name = "Amesos_Superludist";
  else if (Pkg_Choice == "SUPERLU_PARALLEL")
    Pkg_Name = "Amesos_Superludist";
  else if (Pkg_Choice == "SCALAPACK")
    Pkg_Name = "Amesos_Scalapack";
  else if (Pkg_Choice == "MUMPS")
    Pkg_Name = "Amesos_Mumps";
  else {
    std::cout << "Error: Unsupport Amesos solver package" << std::endl;
    exit(-1);
  }

  /* Assemble linear problem */
  Problem.SetOperator(A);
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);

  AMESOS_CHK_ERR(Problem.CheckInput())

  /* Create Amesos base package */
  if (firstSolve) {

    Solver = A_Factory.Create(Pkg_Name.c_str(), Problem);
    if (Solver == 0) {
      std::cout << "Error in amesos_solve_epetra" << std::endl;
      std::cout << "It is likely that the solver package: " << Pkg_Name
          << " has not been linked into the amesos library " << std::endl;
      exit(-1);
    }

    /* Solve problem */
    Solver->SymbolicFactorization();
  }
  Solver->NumericFactorization();
  Solver->Solve();

  /* Convert solution vector */
  int NumMyRows = map.NumMyElements();
  for(int i=0; i< NumMyRows; i++) {
    x_[i] = x[i];
  }

  /* Success! */
  firstSolve = false;
  return 0;
}


static void GomaMsr2EpetraCsr ( struct Aztec_Linear_Solver_System *ams,
				Epetra_CrsMatrix *A )

{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  static int newmatrix = 1;

  int *bindx = ams->bindx;
  double *val = ams->val;

  int NumMyRows = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border]; 
  int NumExternal = ams->data_org[AZ_N_external];
  int NumMyCols = NumMyRows + NumExternal;

  (*A).PutScalar(0.0);

  const Epetra_Map & RowMap = (*A).RowMatrixRowMap();

  int *MyGlobalElements = RowMap.MyGlobalElements ();

  double *dblColGIDs = new double[NumMyCols];
  int *ColGIDs = new int[NumMyCols];

  for( int i=0; i<NumMyRows; i++) dblColGIDs[i] = (double) MyGlobalElements[i];

  AZ_exchange_bdry( dblColGIDs, ams->data_org, ams->proc_config );

  {for ( int j=0; j<NumMyCols ; j++ ) ColGIDs[j] = (int) dblColGIDs[j]; }

  int MaxNNZ = 0;
  for( int i=0; i<NumMyRows; i++)
    {
      int NumNz =  bindx[i+1] - bindx[i]+1;
      
      if( NumNz > MaxNNZ ) MaxNNZ = NumNz;
    }

  int *Indices = new int[MaxNNZ];
  double *Values;

  for( int i=0; i<NumMyRows; i++)
    {
      int NumNz = bindx[i+1] - bindx[i];

      int *k= bindx + bindx[i];
      for( int j=0; j<NumNz; j++ )
	{
	  Indices[j] = ColGIDs[ *k]; k++;
	}
	 
      Values = val + bindx[i];

      if( newmatrix ) 
	{
	  (*A).InsertGlobalValues( MyGlobalElements[i], NumNz, Values, Indices);
	  (*A).InsertGlobalValues( MyGlobalElements[i], 1, &(val[i]), MyGlobalElements+i );
	}
      else
	{
	  (*A).SumIntoGlobalValues( MyGlobalElements[i], NumNz, Values, Indices);
	  (*A).SumIntoGlobalValues( MyGlobalElements[i], 1, &(val[i]), MyGlobalElements+i );
	}      
    }

  (*A).FillComplete();

  delete [] dblColGIDs;
  delete [] ColGIDs;
  delete [] Indices;

  return;

}

void * construct_Epetra_CrsMatrix(  struct Aztec_Linear_Solver_System *ams )
{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  //if ( ams->EpetraMat != NULL ) return ( ams->EpetraMat );

  int *bindx = ams->bindx;
 
  int NumMyRows = ams->data_org[AZ_N_internal] + ams->data_org[AZ_N_border]; 

  Epetra_Map RowMap( -1, NumMyRows, 0, comm );

  int *NumNz = new int[NumMyRows];

  for( int i=0; i<NumMyRows; i++ ) 
    {
      NumNz[i] = bindx[i+1] - bindx[i] + 1;
    }

  Epetra_CrsMatrix *A = new Epetra_CrsMatrix(Copy, RowMap, NumNz);

  delete [] NumNz;

  return( A);
}

#endif  // if defined(ENABLE_AMESOS) && defined(TRILINOS)

#if 0

#include <stdio.h>
#include <iostream.h>
#include <string.h>

#include "sl_util_structs.h"
#include "sl_amesos_interface.h"


void
amesos_solve_msr(char *choice,  
		 struct Aztec_Linear_Solver_System *ams,
		 double *x_, 
		 double *b_,
		 int NewMatrix )
{
	
	
  std::string err_msg("Error: Need to compile with ENABLE_AMESOS flag before using AMESOS solver packages.");
  std::string err_msg2("   Also make sure appropriate libraries are linked for solver packages.");
	
  cout << err_msg << endl;
  cout << err_msg2 << endl;
	
  exit(-1);
	
  return;
}

/* End of if statement for Amesos and Trilinos */
#endif
