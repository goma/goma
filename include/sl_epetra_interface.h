/*
 * sl_epetra_interface.h
 *
 *  Created on: Oct 8, 2014
 *      Author: wortiz
 */

#ifndef INCLUDE_SL_EPETRA_INTERFACE_H_
#define INCLUDE_SL_EPETRA_INTERFACE_H_

#ifdef __cplusplus
typedef Epetra_RowMatrix C_Epetra_RowMatrix_t;
#else
typedef struct Epetra_RowMatrix C_Epetra_RowMatrix_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

C_Epetra_RowMatrix_t* createEpetraRowMatrix(int NumberProcRows);

void insertEpetraRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int LocalRow,
    int NumEntries, const double *Values, const int *Indices);

void sumIntoEpetraRowMatrix(C_Epetra_RowMatrix_t *AMatrix, int LocalRow,
    int NumEntries, const double* Values, const int* Indices);

void putScalarEpetraRowMatrix(C_Epetra_RowMatrix_t *AMatrix, double scalar);

void fillCompleteEpetraRowMatrix(C_Epetra_RowMatrix_t *AMatrix);

void deleteEpetraRowMatrix(C_Epetra_RowMatrix_t *AMatrix);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* INCLUDE_SL_EPETRA_INTERFACE_H_ */
