#ifndef __cplusplus
#define __cplusplus
#endif

#if defined(PARALLEL) && !defined(EPETRA_MPI)
#define EPETRA_MPI
#endif

#include <stdio.h>

#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"

class Epetra_RowMatrix;

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "sl_epetra_interface.h"

/**
 * Create a row matrix with a row map of size NumberProcRows
 * @param NumberProcRows the number of rows on this processor
 * @return C compatible pointer to row matrix
 */
C_Epetra_RowMatrix_t* EpetraCreateRowMatrix(int NumberProcRows) {
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  Epetra_RowMatrix *AMatrix;
  Epetra_Map RowMap(-1, NumberProcRows, 0, comm);

  AMatrix = new Epetra_CrsMatrix(Copy, RowMap, 0);

  return AMatrix;
}

/**
 * Insert values into an epetra row matrix (C interface)
 *
 * Insert if not filled
 *
 * Replace if filled
 *
 * @param AMatrix Matrix to insert values
 * @param GlobalRow Global Row
 * @param NumEntries Number of Entries
 * @param Values Values to insert
 * @param Indices Indices of values
 */
void EpetraInsertGlobalRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int GlobalRow,
    int NumEntries, const double *Values, const int *Indices) {
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(AMatrix);
  if (CrsMatrix->Filled()) {
    CrsMatrix->ReplaceGlobalValues(GlobalRow, NumEntries, Values, Indices);
  } else {
    CrsMatrix->InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
  }
}

/**
 * Sum into global values for epetra row matrix (C interface)
 * @param AMatrix Matrix to sum into values
 * @param GlobalRow Global Row
 * @param NumEntries Number of Entries
 * @param Values Values to sum
 * @param Indices Indices for values
 */
void EpetraSumIntoGlobalRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int GlobalRow,
    int NumEntries, const double* Values, const int* Indices) {
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(AMatrix);
  int ierr = CrsMatrix->SumIntoGlobalValues(GlobalRow, NumEntries, Values,
      Indices);
  if (ierr)
    printf("Error in sum into epetra %d\n", ierr);
}

/**
 * Set an Epetra Row Matrix to a specified scalar value for all NZ entries (C interface)
 * @param AMatrix Row matrix to set
 * @param scalar scalar value to set
 */
void EpetraPutScalarRowMatrix(C_Epetra_RowMatrix_t *AMatrix, double scalar) {
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(AMatrix);
  CrsMatrix->PutScalar(scalar);
}

/**
 * Call fill complete on an Epetra Row Matrix (C interface)
 * @param AMatrix Matrix to fill completely
 */
void EpetraFillCompleteRowMatrix(C_Epetra_RowMatrix_t *AMatrix) {
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(AMatrix);
  CrsMatrix->FillComplete();
}

int EpetraExtractRowValuesRowMatrix(C_Epetra_RowMatrix_t *AMatrix,
    int GlobalRow, int size, double *values, int *indices) {
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(AMatrix);
  int NumEntries;
  CrsMatrix->ExtractGlobalRowCopy(GlobalRow, size, NumEntries, values, indices);
  return NumEntries;
}

/**
 * Call deconstructor for row matrix (C interface)
 * @param AMatrix Row Matrix to delete
 */
void EpetraDeleteRowMatrix(C_Epetra_RowMatrix_t *AMatrix) {
  delete AMatrix;
}

