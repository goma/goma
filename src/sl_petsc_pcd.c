#ifdef GOMA_ENABLE_PETSC
#include <petscsystypes.h>
#if !(PETSC_USE_COMPLEX)
#include "sl_petsc.h"
#include <petscksp.h>

#include "dp_comm.h"
#include "dpi.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "rf_vars_const.h"
#include "sl_petsc.h"

static goma_error initialize_PCD_matrices(PetscMatrixData *matrix_data,
                                          Exo_DB *exo,
                                          Dpi *dpi,
                                          dbl *x,
                                          dbl *x_old,
                                          dbl *xdot,
                                          dbl *xdot_old,
                                          PetscInt local_nodes,
                                          PetscInt global_nodes) {
  PetscInt global_offset = 0;
  if (Num_Proc > 1) {
    MPI_Scan(&local_nodes, &global_offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
    global_offset -= local_nodes;
  }
  PetscInt num_rows = dpi->num_internal_nodes + dpi->num_boundary_nodes;

  PetscInt *d_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt *o_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));

  for (int eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++) {
    int mn = Matilda[eb_index];

    pd = pd_glob[mn];
    cr = cr_glob[mn];
    elc = elc_glob[mn];
    elc_rs = elc_rs_glob[mn];
    gn = gn_glob[mn];
    mp = mp_glob[mn];
    vn = vn_glob[mn];
    evpl = evpl_glob[mn];

    for (int mode = 0; mode < vn->modes; mode++) {
      ve[mode] = ve_glob[mn][mode];
    }

    int e_start = exo->eb_ptr[eb_index];
    int e_end = exo->eb_ptr[eb_index + 1];

    for (int iel = e_start; iel < e_end; iel++) {
      int ielem = iel;

      int err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);
      GOMA_EH(err, "load_elem_dofptr");
      err = bf_mp_init(pd);
      GOMA_EH(err, "bf_mp_init");

      int eqn = PRESSURE;
      for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
        int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
        int ldof_i = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][i];
        for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
          int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
          int ldof_j = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][j];
          if (ldof_i >= 0 && ldof_j >= 0 &&
              ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
            if (gnn_i < num_rows) {
              PetscInt global_i = matrix_data->schur_s_local_to_global[gnn_i];
              if (gnn_j >= num_rows) {
                int index = global_i - global_offset;
                o_nnz[index] += 1;
              } else {
                int index = global_i - global_offset;
                d_nnz[index] += 1;
              }
            }
          }
        }
      }
    } /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }   /* END for (ieb loop) */

  if (Num_Proc == 1) {
    MatSeqAIJSetPreallocation(matrix_data->pcd_data->Mp, 0, d_nnz);
    MatSeqAIJSetPreallocation(matrix_data->pcd_data->Mp_mu, 0, d_nnz);
    MatSeqAIJSetPreallocation(matrix_data->pcd_data->Fp, 0, d_nnz);
    MatSeqAIJSetPreallocation(matrix_data->pcd_data->Ap, 0, d_nnz);
  } else {
    MatMPIAIJSetPreallocation(matrix_data->pcd_data->Mp, 0, d_nnz, 0, o_nnz);
    MatMPIAIJSetPreallocation(matrix_data->pcd_data->Mp_mu, 0, d_nnz, 0, o_nnz);
    MatMPIAIJSetPreallocation(matrix_data->pcd_data->Fp, 0, d_nnz, 0, o_nnz);
    MatMPIAIJSetPreallocation(matrix_data->pcd_data->Ap, 0, d_nnz, 0, o_nnz);
  }

  free(d_nnz);
  free(o_nnz);
  return GOMA_SUCCESS;
}

PetscErrorCode PCDShellPCApply(PC pc, Vec x, Vec y) {
  PetscPCDData *pcd_data;
  PCShellGetContext(pc, (void **)&pcd_data);
  Vec a, b, x_tmp, my;
  VecDuplicate(x, &a);
  VecDuplicate(x, &b);
  VecDuplicate(x, &x_tmp);
  VecCopy(x, x_tmp);
  VecDuplicate(x, &my);
  if (pcd_data->pcd_inverse_diag) {
    Vec diag;
    VecDuplicate(x, &diag);
    MatGetDiagonal(pcd_data->Mp, diag);
    VecReciprocal(diag);
    VecPointwiseMult(a, diag, x);
    VecDestroy(&diag);
  } else {
    KSPSolve(pcd_data->ksp_Mp, x, a);
  }
  MatMult(pcd_data->Fp, a, b);
  KSPSolve(pcd_data->ksp_Ap, b, y);
  KSPSolve(pcd_data->ksp_Mp_mu, x_tmp, my);
  VecAXPY(y, 1.0, my);
  VecDestroy(&a);
  VecDestroy(&b);
  VecDestroy(&x_tmp);
  VecDestroy(&my);
  return 0;
}

PetscErrorCode petsc_PCD_setup(PC pc,
                               PetscMatrixData *matrix_data,
                               Exo_DB *exo,
                               Dpi *dpi,
                               dbl *x,
                               dbl *x_old,
                               dbl *xdot,
                               dbl *xdot_old) {
  PetscPCDData *data;
  PetscNew(&data);
  matrix_data->pcd_data = data;

  data->pcd_inverse_diag = matrix_data->user_pcd_inverse_diag;

  PetscInt global_nodes;
  PetscInt local_nodes;
  count_pressure_nodes(matrix_data, exo, dpi, x, x_old, xdot, xdot_old, &local_nodes,
                       &global_nodes);

  PetscErrorCode err = MatCreate(PETSC_COMM_WORLD, &data->Fp);
  CHKERRQ(err);
  err = MatCreate(PETSC_COMM_WORLD, &data->Mp);
  CHKERRQ(err);
  err = MatCreate(PETSC_COMM_WORLD, &data->Mp_mu);
  CHKERRQ(err);
  err = MatCreate(PETSC_COMM_WORLD, &data->Ap);
  CHKERRQ(err);
  err =
      MatSetSizes(matrix_data->pcd_data->Mp, local_nodes, local_nodes, global_nodes, global_nodes);
  CHKERRQ(err);
  err = MatSetSizes(matrix_data->pcd_data->Mp_mu, local_nodes, local_nodes, global_nodes,
                    global_nodes);
  CHKERRQ(err);
  err =
      MatSetSizes(matrix_data->pcd_data->Ap, local_nodes, local_nodes, global_nodes, global_nodes);
  CHKERRQ(err);
  err =
      MatSetSizes(matrix_data->pcd_data->Fp, local_nodes, local_nodes, global_nodes, global_nodes);
  CHKERRQ(err);

  err = KSPCreate(MPI_COMM_WORLD, &data->ksp_Mp);
  CHKERRQ(err);
  err = KSPCreate(MPI_COMM_WORLD, &data->ksp_Mp_mu);
  CHKERRQ(err);
  err = KSPCreate(MPI_COMM_WORLD, &data->ksp_Ap);
  CHKERRQ(err);
  err = MatSetOptionsPrefix(data->Mp, "pcdMp_");
  CHKERRQ(err);
  err = MatSetOptionsPrefix(data->Ap, "pcdAp_");
  CHKERRQ(err);
  err = MatSetOptionsPrefix(data->Fp, "pcdFp_");
  CHKERRQ(err);
  err = MatSetOptionsPrefix(data->Fp, "pcdMp_mu_");
  CHKERRQ(err);
  err = KSPSetOptionsPrefix(data->ksp_Ap, "innerAp_");
  CHKERRQ(err);
  err = KSPSetOptionsPrefix(data->ksp_Mp, "innerMp_");
  CHKERRQ(err);
  err = KSPSetOptionsPrefix(data->ksp_Mp_mu, "innerMp_mu_");
  CHKERRQ(err);
  err = MatSetFromOptions(data->Mp);
  CHKERRQ(err);
  err = MatSetFromOptions(data->Mp_mu);
  CHKERRQ(err);
  err = MatSetFromOptions(data->Fp);
  CHKERRQ(err);
  err = MatSetFromOptions(data->Ap);
  CHKERRQ(err);
  err = PCSetType(pc, PCFIELDSPLIT);
  CHKERRQ(err);

  initialize_PCD_matrices(matrix_data, exo, dpi, x, x_old, xdot, xdot_old, local_nodes,
                          global_nodes);

  err = KSPSetOperators(data->ksp_Ap, data->Ap, data->Ap);
  CHKERRQ(err);
  err = KSPSetOperators(data->ksp_Mp, data->Mp, data->Mp);
  CHKERRQ(err);
  err = KSPSetOperators(data->ksp_Mp_mu, data->Mp_mu, data->Mp_mu);
  CHKERRQ(err);
  err = KSPSetFromOptions(data->ksp_Mp_mu);
  CHKERRQ(err);
  err = KSPSetFromOptions(data->ksp_Mp);
  CHKERRQ(err);
  err = KSPSetFromOptions(data->ksp_Ap);
  CHKERRQ(err);
  err = PCSetUp(pc);
  CHKERRQ(err);
  KSP *subksp;
  PetscInt n_splits;
  err = PCFieldSplitSchurGetSubKSP(pc, &n_splits, &subksp);
  CHKERRQ(err);
  PC spc;
  GOMA_ASSERT_ALWAYS(n_splits > 1);
  err = KSPGetPC(subksp[1], &spc);
  CHKERRQ(err);
  PCSetType(spc, PCSHELL);
  PCShellSetApply(spc, PCDShellPCApply);
  PCShellSetName(spc, "GomaPCD");
  PCShellSetContext(spc, (void *)data);

  return 0;
}
#endif
#endif