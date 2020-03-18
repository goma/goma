/*
 * sl_epetra_interface.h
 *
 *  Created on: Oct 8, 2014
 *      Author: wortiz
 */

#ifndef INCLUDE_SL_EPETRA_INTERFACE_H_
#define INCLUDE_SL_EPETRA_INTERFACE_H_

#ifdef __cplusplus
#include "Epetra_RowMatrix.h"
typedef Epetra_RowMatrix C_Epetra_RowMatrix_t;
#else
typedef struct Epetra_RowMatrix C_Epetra_RowMatrix_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

C_Epetra_RowMatrix_t* EpetraCreateRowMatrix(int NumberProcRows);

void EpetraInsertGlobalRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int GlobalRow,
    int NumEntries, const double *Values, const int *Indices);

void EpetraSumIntoGlobalRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int GlobalRow,
    int NumEntries, const double* Values, const int* Indices);

void EpetraPutScalarRowMatrix(C_Epetra_RowMatrix_t *AMatrix, double scalar);

void EpetraFillCompleteRowMatrix(C_Epetra_RowMatrix_t *AMatrix);

int EpetraExtractRowValuesRowMatrix(C_Epetra_RowMatrix_t *AMatrix,
    int GlobalRow, int size, double *values, int *indices);

void EpetraDeleteRowMatrix(C_Epetra_RowMatrix_t *AMatrix);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* INCLUDE_SL_EPETRA_INTERFACE_H_ */
