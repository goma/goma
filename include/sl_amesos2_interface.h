#ifndef GOMA_SL_AMESOS2_INTERFACE_H_
#define GOMA_SL_AMESOS2_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "rf_fem_const.h"
#include "rf_io_const.h"

int amesos2_solve(struct GomaLinearSolverData *ams,
                  double *x_,
                  double *b_,
                  char *amesos2_solver,
                  char *amesos2_file);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* GOMA_SL_AMESOS2_INTERFACE_H_ */
