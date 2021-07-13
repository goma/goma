#ifndef GOMA_SL_PETSC_H
#define GOMA_SL_PETSC_H

#include "dpi.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_eh.h"

goma_error goma_setup_petsc_matrix(struct GomaLinearSolverData *ams,
                                   Exo_DB *exo,
                                   Dpi *dpi,
                                   int internal_dof,
                                   int boundary_dof,
                                   int external_dof,
                                   int imtrx);

void petsc_load_lec(int ielem, struct GomaLinearSolverData *ams,
                   double resid_vector[]);
void
petsc_solve(struct GomaLinearSolverData *ams,
		  double *x_, 
		  double *b_,
                  int *its);

goma_error petsc_scale_matrix(struct GomaLinearSolverData *ams,
		  double *b_,
                  double *scale);
#endif // GOMA_SL_PETSC_H
