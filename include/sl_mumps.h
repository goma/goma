#ifndef GOMA_SL_MUMPS_H
#define GOMA_SL_MUMPS_H

#include "std.h"
#include "sl_util_structs.h"

void mumps_solve(struct GomaLinearSolverData *data, dbl *x, dbl *rhs);

#endif // GOMA_SL_MUMPS_H