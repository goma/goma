/*
 * sl_epetra_util.h
 *
 */

#ifndef INCLUDE_SL_EPETRA_UTIL_H_
#define INCLUDE_SL_EPETRA_UTIL_H_

#include "dpi.h"
#include "exo_struct.h"
#ifdef __cplusplus
extern "C" {
#endif

void EpetraCreateGomaProblemGraph(struct GomaLinearSolverData *ams, Exo_DB *exo, Dpi *dpi);

void EpetraLoadLec(int ielem, struct GomaLinearSolverData *ams, double resid_vector[]);

void EpetraRowSumScale(struct GomaLinearSolverData *ams, double *b, double *scale);

void EpetraSetDiagonalOnly(struct GomaLinearSolverData *ams, int GlobalRow);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* INCLUDE_SL_EPETRA_UTIL_H_ */
