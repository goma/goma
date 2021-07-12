#ifndef GOMA_ADAPT_RESETUP_PROBLEM_H
#define GOMA_ADAPT_RESETUP_PROBLEM_H

#include "exo_struct.h"
#include "dpi.h"
#include "mm_as_structs.h"

int resetup_problem(Exo_DB *exo, Dpi *dpi);	/* ptr to the finite element mesh database */
int resetup_matrix(struct GomaLinearSolverData **ams, Exo_DB *exo, Dpi* dpi);

#endif // GOMA_ADAPT_RESETUP_PROBLEM_H

// vim: expandtab sw=2 ts=8
