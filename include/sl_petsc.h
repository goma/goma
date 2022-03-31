#ifndef GOMA_SL_PETSC_H
#define GOMA_SL_PETSC_H

#ifdef GOMA_ENABLE_PETSC
#include <petscsystypes.h>
#ifdef I
#undef I
#endif
#if !(PETSC_USE_COMPLEX)
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

typedef struct PetscPCDData {
  Mat SchurS;
  Mat Ap;
  Mat Mp;
  Mat Mp_mu;
  Mat Fp;
  KSP ksp_Ap;
  KSP ksp_Mp;
  KSP ksp_Mp_mu;
  PetscInt *schur_s_local_to_global;
  PetscBool pcd_inverse_diag;
} PetscPCDData;

typedef struct PetscMatrixData {
  PetscOptions options;
  Mat mat;
  Mat SchurS;
  Vec residual;
  Vec update;
  KSP ksp;
  PetscInt *local_to_global;
  PetscInt *schur_s_local_to_global;
  PetscBool matrix_setup;
  PetscBool user_schur;
  PetscBool user_pcd;
  PetscBool user_pcd_inverse_diag;
  PetscPCDData *pcd_data;
  PetscInt pcd_ss_remove_n;
  PetscInt *pcd_ss_remove;
  PetscInt pcd_ns_remove_n;
  PetscInt *pcd_ns_remove;
  PetscInt *real_to_complex;
} PetscMatrixData;

PetscErrorCode petsc_PCD_setup(PC ppc,
                               PetscMatrixData *matrix_data,
                               Exo_DB *exo,
                               Dpi *dpi,
                               dbl *x,
                               dbl *x_old,
                               dbl *xdot,
                               dbl *xdot_old);

goma_error count_pressure_nodes(PetscMatrixData *matrix_data,
                                Exo_DB *exo,
                                Dpi *dpi,
                                dbl *x,
                                dbl *x_old,
                                dbl *xdot,
                                dbl *xdot_old,
                                PetscInt *local_nodes,
                                PetscInt *global_nodes);

goma_error goma_setup_petsc_matrix(struct GomaLinearSolverData *ams,
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

goma_error goma_setup_petsc_post_proc_matrix(
    Exo_DB *exo, Dpi *dpi, dbl *x, dbl *x_old, dbl *xdot, dbl *xdot_old);
void petsc_solve_post_proc(double **post_proc_vect, RESULTS_DESCRIPTION_STRUCT *rd, Dpi *dpi);

void petsc_load_lec(int ielem, struct GomaLinearSolverData *ams, double resid_vector[]);
int petsc_solve(struct GomaLinearSolverData *ams, double *x_, double *b_, int *its);

goma_error petsc_scale_matrix(struct GomaLinearSolverData *ams, double *b_, double *scale);

goma_error goma_petsc_free_matrix(struct GomaLinearSolverData *ams);
#endif
#endif
#endif // GOMA_SL_PETSC_H
