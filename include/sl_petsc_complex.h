#ifndef GOMA_SL_PETSC_COMPLEX
#define GOMA_SL_PETSC_COMPLEX
#ifdef GOMA_ENABLE_PETSC
#include <petscconf.h>
#if PETSC_USE_COMPLEX
#include "dpi.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "rf_io_structs.h"
#include <petscksp.h>
#include <petscsys.h>
#ifdef I
#undef I
#endif

typedef struct PetscMatrixData {
  PetscOptions options;
  Mat mat;
  Vec residual;
  Vec update;
  KSP ksp;
  PetscInt local_dof;
  PetscInt global_dof;
  PetscInt *real_local_to_goma_local;
  PetscInt *imag_local_to_goma_local;
  PetscInt *local_to_global_complex;
  PetscInt *local_to_global;
  PetscBool *is_imag;
  PetscBool matrix_setup;
  int *o_nnz;
  int *d_nnz;
} PetscMatrixData;

PetscErrorCode petsc_reset_ksp_mat(struct GomaLinearSolverData *ams);
goma_error goma_setup_petsc_matrix_complex(struct GomaLinearSolverData *ams,
                                           Exo_DB *exo,
                                           Dpi *dpi,
                                           dbl *x,
                                           dbl *x_old,
                                           dbl *xdot,
                                           dbl *xdot_old,
                                           int internal_dof,
                                           int boundary_dof,
                                           int external_dof,
                                           int imtrx);

void petsc_load_lec_complex(int ielem, struct GomaLinearSolverData *ams, double resid_vector[]);

int petsc_solve_complex(struct GomaLinearSolverData *ams, double *x_, double *b_, int *its);
goma_error goma_petsc_free_matrix(struct GomaLinearSolverData *ams);

#endif
#endif
#endif