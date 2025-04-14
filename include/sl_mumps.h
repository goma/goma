#ifndef GOMA_SL_MUMPS_H
#define GOMA_SL_MUMPS_H

#include "mm_eh.h"
#include "sl_util_structs.h"
#include "std.h"

goma_error mumps_solve(struct GomaLinearSolverData *data, dbl *x, dbl *rhs);

#endif // GOMA_SL_MUMPS_H