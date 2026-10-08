#ifndef INCLUDE_SL_STRATIMIKOS_INTERFACE_H_
#define INCLUDE_SL_STRATIMIKOS_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_as_structs.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"

int stratimikos_solve_tpetra(struct GomaLinearSolverData *ams,
                             double *x_,
                             double *b_,
                             int *iterations,
                             char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                             int imtrx);

int stratimikos_solve(struct GomaLinearSolverData *ams,
                      double *x_,
                      double *b_,
                      int *iterations,
                      char stratimikos_file[MAX_NUM_MATRICES][MAX_CHAR_IN_INPUT],
                      int imtrx);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* INCLUDE_SL_STRATIMIKOS_INTERFACE_H_ */
