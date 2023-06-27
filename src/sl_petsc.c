#include "load_field_variables.h"
#include <mpi.h>
#include <stddef.h>
#ifdef GOMA_ENABLE_PETSC
#include <petscsystypes.h>
#if !(PETSC_USE_COMPLEX)
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>
#ifdef I
#undef I
#endif

#include "density.h"
#include "dp_comm.h"
#include "dpi.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_fill_aux.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rf_bc.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "rf_vars_const.h"
#include "sl_petsc.h"

static bool GomaPetscOptionsInserted = false;
static bool GomaPetscOptionsPrinted = false;

static int pint_compare(const void *a, const void *b) {
  const PetscInt *ia = (const PetscInt *)a;
  const PetscInt *ib = (const PetscInt *)b;
  return *ia - *ib;
}

static goma_error get_pressure_velocity_is(PetscMatrixData *matrix_data,
                                           Exo_DB *exo,
                                           Dpi *dpi,
                                           dbl *x,
                                           dbl *x_old,
                                           dbl *xdot,
                                           dbl *xdot_old,
                                           PetscInt *vel_n,
                                           PetscInt **vel_is,
                                           PetscInt *pres_n,
                                           PetscInt **pres_is) {

  PetscInt num_cols = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

  PetscInt *nzp = (PetscInt *)calloc(num_cols, sizeof(PetscInt));
  PetscInt *nzv[3];
  nzv[0] = (PetscInt *)calloc(num_cols, sizeof(PetscInt));
  nzv[1] = (PetscInt *)calloc(num_cols, sizeof(PetscInt));
  nzv[2] = (PetscInt *)calloc(num_cols, sizeof(PetscInt));
  for (int i = 0; i < num_cols; i++) {
    nzp[i] = -1;
    nzv[0][i] = -1;
    nzv[1][i] = -1;
    nzv[2][i] = -1;
  }

  int index_pres = 0;
  int index_vel = 0;
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
        if (ldof_i >= 0 &&
            ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
          if (nzp[gnn_i] == -1) {
            nzp[gnn_i] = index_pres;
            index_pres++;
          }
        }
      }
      eqn = VELOCITY1;
      for (int var = eqn; var <= VELOCITY3; var++) {
        if (upd->vp[pg->imtrx][var] >= 0) {
          for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
            int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
            int ldof_i = ei[upd->matrix_index[var]]->ln_to_dof[var][i];
            if (ldof_i >= 0 &&
                ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
              if (nzv[var - VELOCITY1][gnn_i] == -1) {
                nzv[var - VELOCITY1][gnn_i] = index_vel;
                index_vel++;
              }
            }
          }
        }
      }
    } /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }   /* END for (ieb loop) */

  *vel_n = index_vel;
  *pres_n = index_pres;
  PetscMalloc1(index_vel, vel_is);
  PetscMalloc1(index_pres, pres_is);
  for (int i = 0; i < num_cols; i++) {
    nzp[i] = -1;
    nzv[0][i] = -1;
    nzv[1][i] = -1;
    nzv[2][i] = -1;
  }
  index_pres = 0;
  index_vel = 0;
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

      int v = PRESSURE;
      for (int j = 0; j < ei[pg->imtrx]->dof[v]; j++) {
        int ledof = ei[pg->imtrx]->lvdof_to_ledof[v][j];
        int kv = 0;
        int row_index =
            Index_Solution(ei[pg->imtrx]->gnn_list[v][j], v, kv, ei[pg->imtrx]->Baby_Dolphin[v][j],
                           ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
        GOMA_EH(row_index, "Bad var index.");
        PetscInt global_row = matrix_data->local_to_global[row_index];
        int gnn_i = ei[pg->imtrx]->gnn_list[v][j];
        if ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1) {
          if (nzp[gnn_i] == -1) {
            nzp[gnn_i] = 0;
            (*pres_is)[index_pres] = global_row;
            index_pres++;
          }
        }
      }
      for (int v = VELOCITY1; v <= VELOCITY3; v++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[v]; j++) {
          int ledof = ei[pg->imtrx]->lvdof_to_ledof[v][j];
          int kv = 0;
          int row_index = Index_Solution(ei[pg->imtrx]->gnn_list[v][j], v, kv,
                                         ei[pg->imtrx]->Baby_Dolphin[v][j],
                                         ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
          GOMA_EH(row_index, "Bad var index.");
          PetscInt global_row = matrix_data->local_to_global[row_index];
          int gnn_i = ei[pg->imtrx]->gnn_list[v][j];
          if ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1) {
            if (nzv[v - VELOCITY1][gnn_i] == -1) {
              nzv[v - VELOCITY1][gnn_i] = 0;
              (*vel_is)[index_vel] = global_row;
              index_vel++;
            }
          }
        }
      }
    } /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }   /* END for (ieb loop) */
  GOMA_ASSERT(index_vel == *vel_n);
  GOMA_ASSERT(index_pres == *pres_n);
  qsort(*vel_is, index_vel, sizeof(PetscInt), pint_compare);
  qsort(*pres_is, index_pres, sizeof(PetscInt), pint_compare);
  free(nzp);
  free(nzv[0]);
  free(nzv[1]);
  free(nzv[2]);
  return GOMA_SUCCESS;
}
goma_error count_pressure_nodes(PetscMatrixData *matrix_data,
                                Exo_DB *exo,
                                Dpi *dpi,
                                dbl *x,
                                dbl *x_old,
                                dbl *xdot,
                                dbl *xdot_old,
                                PetscInt *local_nodes,
                                PetscInt *global_nodes) {

  PetscInt num_cols = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

  PetscInt *nz = (PetscInt *)calloc(num_cols, sizeof(PetscInt));
  for (int i = 0; i < num_cols; i++) {
    nz[i] = -1;
  }

  int index = 0;
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
        if (ldof_i >= 0 &&
            ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
          if (nz[gnn_i] == -1) {
            nz[gnn_i] = index;
            index++;
          }
        }
      }
    } /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }   /* END for (ieb loop) */

  *local_nodes = index;

  MPI_Allreduce(local_nodes, global_nodes, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
  PetscInt global_offset = 0;
  if (Num_Proc > 1) {
    MPI_Scan(local_nodes, &global_offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
    global_offset -= *local_nodes;
  }
  for (int i = 0; i < (dpi->num_internal_nodes + dpi->num_boundary_nodes); i++) {
    if (nz[i] != -1) {
      nz[i] += global_offset;
    }
  }
  dbl *dblColGIDs = (dbl *)calloc(num_cols, sizeof(dbl));
  for (int i = 0; i < num_cols; i++) {
    dblColGIDs[i] = (double)nz[i];
  }
  matrix_data->schur_s_local_to_global = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);
  exchange_node(cx[0], dpi, dblColGIDs);
  for (int i = 0; i < num_cols; i++) {
    matrix_data->schur_s_local_to_global[i] = (PetscInt)dblColGIDs[i];
  }
  free(dblColGIDs);
  free(nz);
  return GOMA_SUCCESS;
}

static goma_error initialize_petsc_pressure_matrix(PetscMatrixData *matrix_data,
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
    MatSeqAIJSetPreallocation(matrix_data->SchurS, 0, d_nnz);
  } else {
    MatMPIAIJSetPreallocation(matrix_data->SchurS, 0, d_nnz, 0, o_nnz);
  }

  free(d_nnz);
  free(o_nnz);
  return GOMA_SUCCESS;
}

static goma_error set_pcd_matrices(PetscMatrixData *matrix_data,
                                   Exo_DB *exo,
                                   Dpi *dpi,
                                   dbl *x,
                                   dbl *x_old,
                                   dbl *xdot,
                                   dbl *xdot_old) {
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
      int ielem_type = ei[pg->imtrx]->ielem_type;
      int ip_total = elem_info(NQUAD, ielem_type); /* number of
                                                    * quadrature pts */

      for (int ip = 0; ip < ip_total; ip++) {
        dbl xi[3];
        dbl s, t, u;

        find_stu(ip, ielem_type, &s, &t, &u);
        xi[0] = s;
        xi[1] = t;
        xi[2] = u;

        /*
         * find quadrature weights for current ip
         */
        dbl wt = Gq_weight(ip, ielem_type);
        fv->wt = wt;

        /*
         * Load up basis function information for ea variable...
         * Old usage: fill_shape
         */
        err = load_basis_functions(xi, bfd);
        GOMA_EH(err, "problem from load_basis_functions");

        /*
         * This has elemental Jacobian transformation and some
         * basic mesh derivatives...
         * Old usage: calc_Jac, jelly_belly
         */
        err = beer_belly();
        GOMA_EH(err, "beer_belly");

        err = load_fv();
        GOMA_EH(err, "load_fv");

        err = load_bf_grad();
        GOMA_EH(err, "load_bf_grad");

        err = load_fv_grads();
        GOMA_EH(err, "load_fv_grads");

        double gamma[DIM][DIM];
        for (int a = 0; a < VIM; a++) {
          for (int b = 0; b < VIM; b++) {
            gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
          }
        }
        double mu = viscosity(gn, gamma, NULL);
        double invmu = (1.0 / mu);

        dbl rho = density(NULL, tran->time_value);

        int eqn = PRESSURE;
        for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
          int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          int ldof_i = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][i];
          for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            int ldof_j = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][j];
            if (ldof_i >= 0 && ldof_j >= 0 &&
                ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {

              PetscInt global_row = matrix_data->schur_s_local_to_global[gnn_i];
              PetscInt global_col = matrix_data->schur_s_local_to_global[gnn_j];
              GOMA_EH(global_row, "Row not found");
              GOMA_EH(global_col, "Col not found");
              // Mp
              dbl mp_contrib = bf[eqn]->phi[ldof_i] * bf[eqn]->phi[ldof_j] * fv->wt * bf[eqn]->detJ;
              PetscErrorCode err = MatSetValue(matrix_data->pcd_data->Mp, global_row, global_col,
                                               mp_contrib, ADD_VALUES);
              CHKERRQ(err);

              // Mp_mu
              dbl mu_contrib = 0.5 * invmu * bf[eqn]->phi[ldof_i] * bf[eqn]->phi[ldof_j] * fv->wt *
                               bf[eqn]->detJ;
              err = MatSetValue(matrix_data->pcd_data->Mp_mu, global_row, global_col, mu_contrib,
                                ADD_VALUES);
              CHKERRQ(err);

              // Ap
              dbl ap_contrib = 0.0;
              for (int a = 0; a < pd->Num_Dim; a++) {
                ap_contrib += mu * bf[eqn]->grad_phi[ldof_i][a] * bf[eqn]->grad_phi[ldof_j][a] *
                              fv->wt * bf[eqn]->detJ;
              }
              err = MatSetValue(matrix_data->pcd_data->Ap, global_row, global_col, ap_contrib,
                                ADD_VALUES);
              CHKERRQ(err);

              // Fp
              dbl fp_contrib = ap_contrib;
              for (int a = 0; a < pd->Num_Dim; a++) {
                fp_contrib += rho * fv->v[a] * bf[eqn]->grad_phi[ldof_j][a] * bf[eqn]->phi[ldof_i] *
                              fv->wt * bf[eqn]->detJ;
              }
              if (pd->TimeIntegration != STEADY) {
                fp_contrib += (pd->etm[pg->imtrx][VELOCITY1][LOG2_MASS] / tran->delta_t) * rho *
                              bf[eqn]->phi[ldof_i] * bf[eqn]->phi[ldof_j] * fv->wt * bf[eqn]->detJ;
              }
              err = MatSetValue(matrix_data->pcd_data->Fp, global_row, global_col, fp_contrib,
                                ADD_VALUES);
              CHKERRQ(err);
            }
          }
        }
      } /* END  for (ip = 0; ip < ip_total; ip++)                      */
    }   /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }     /* END for (ieb loop) */

  PetscInt n_ss_rows = 0;
  for (int ss_idx = 0; ss_idx < matrix_data->pcd_ss_remove_n; ss_idx++) {
    int set_index = -1;
    for (int i = 0; i < exo->num_side_sets; i++) {
      if (exo->ss_id[i] == matrix_data->pcd_ss_remove[ss_idx]) {
        set_index = i;
        break;
      }
    }
    if (set_index == -1) {
      GOMA_WH(GOMA_ERROR, "SS not found for PCD removal %d", matrix_data->pcd_ss_remove[ss_idx]);
      continue;
    }
    int elem_start = exo->ss_elem_index[set_index];
    for (int q = 0; q < exo->ss_num_sides[set_index]; q++) {
      int ielem = exo->ss_elem_list[elem_start + q];
      int iside = exo->ss_side_list[elem_start + q];
      int snl[MAX_NODES_PER_SIDE];
      int num_nodes = build_side_node_list(ielem, iside - 1, exo, snl);
      // Dirichlet nodes on side
      for (int n = 0; n < num_nodes; n++) {
        int gnn = snl[n];
        if (gnn < dpi->num_internal_nodes + dpi->num_boundary_nodes) {
          PetscInt global_row = matrix_data->schur_s_local_to_global[gnn];
          if (global_row != -1) {
            n_ss_rows++;
          }
        }
      }
    }
  }
  PetscInt *ss_rows;
  ss_rows = malloc(sizeof(PetscInt) * n_ss_rows);
  PetscInt index = 0;
  for (int ss_idx = 0; ss_idx < matrix_data->pcd_ss_remove_n; ss_idx++) {
    int set_index = -1;
    for (int i = 0; i < exo->num_side_sets; i++) {
      if (exo->ss_id[i] == matrix_data->pcd_ss_remove[ss_idx]) {
        set_index = i;
        break;
      }
    }
    if (set_index == -1) {
      continue;
    }
    int elem_start = exo->ss_elem_index[set_index];
    for (int q = 0; q < exo->ss_num_sides[set_index]; q++) {
      int ielem = exo->ss_elem_list[elem_start + q];
      int iside = exo->ss_side_list[elem_start + q];
      int snl[MAX_NODES_PER_SIDE];
      int num_nodes = build_side_node_list(ielem, iside - 1, exo, snl);
      // Dirichlet nodes on side
      for (int n = 0; n < num_nodes; n++) {
        int gnn = snl[n];
        if (gnn < dpi->num_internal_nodes + dpi->num_boundary_nodes) {
          PetscInt global_row = matrix_data->schur_s_local_to_global[gnn];
          if (global_row != -1) {
            ss_rows[index] = global_row;
            index++;
          }
        }
      }
    }
  }
  PetscInt n_ns_rows = 0;
  for (int ns_idx = 0; ns_idx < matrix_data->pcd_ns_remove_n; ns_idx++) {
    int set_index = -1;
    for (int i = 0; i < exo->num_node_sets; i++) {
      if (exo->ns_id[i] == matrix_data->pcd_ns_remove[ns_idx]) {
        set_index = i;
        break;
      }
    }
    if (set_index == -1) {
      GOMA_WH(GOMA_ERROR, "NS not found for PCD removal %d", matrix_data->pcd_ns_remove[ns_idx]);
      continue;
    }
    for (int i = 0; i < exo->ns_num_nodes[set_index]; i++) {
      int gnn = exo->ns_node_list[exo->ns_node_index[set_index] + i];
      PetscInt global_row = matrix_data->schur_s_local_to_global[gnn];
      if (global_row != -1) {
        n_ns_rows++;
      }
    }
  }

  PetscInt *ns_rows = malloc(sizeof(PetscInt) * n_ns_rows);
  index = 0;
  for (int ns_idx = 0; ns_idx < matrix_data->pcd_ns_remove_n; ns_idx++) {
    int set_index = -1;
    for (int i = 0; i < exo->num_node_sets; i++) {
      if (exo->ns_id[i] == matrix_data->pcd_ns_remove[ns_idx]) {
        set_index = i;
        break;
      }
    }
    if (set_index == -1) {
      continue;
    }

    for (int i = 0; i < exo->ns_num_nodes[set_index]; i++) {
      int gnn = exo->ns_node_list[exo->ns_node_index[set_index] + i];
      PetscInt global_row = matrix_data->schur_s_local_to_global[gnn];
      if (global_row != -1) {
        ns_rows[index] = global_row;
        index++;
      }
    }
  }
  MatAssemblyBegin(matrix_data->pcd_data->Ap, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->pcd_data->Ap, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(matrix_data->pcd_data->Mp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->pcd_data->Mp, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(matrix_data->pcd_data->Fp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->pcd_data->Fp, MAT_FINAL_ASSEMBLY);
  // MatZeroRowsColumns(matrix_data->pcd_data->Fp, n_ss_rows, ss_rows, 1.0, NULL, NULL);
  if (matrix_data->pcd_ss_remove_n > 0) {
    MatZeroRowsColumns(matrix_data->pcd_data->Ap, n_ss_rows, ss_rows, 1.0, NULL, NULL);
  }
  if (matrix_data->pcd_ns_remove_n > 0) {
    MatZeroRowsColumns(matrix_data->pcd_data->Ap, n_ns_rows, ns_rows, 1.0, NULL, NULL);
  }
  // MatZeroRowsColumns(matrix_data->pcd_data->Mp, n_ss_rows, ss_rows, 1e-4, NULL, NULL);
  free(ss_rows);
  free(ns_rows);
  return GOMA_SUCCESS;
}
static goma_error set_petsc_pressure_matrix(PetscMatrixData *matrix_data,
                                            Exo_DB *exo,
                                            Dpi *dpi,
                                            dbl *x,
                                            dbl *x_old,
                                            dbl *xdot,
                                            dbl *xdot_old) {
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
      int ielem_type = ei[pg->imtrx]->ielem_type;
      int ip_total = elem_info(NQUAD, ielem_type); /* number of
                                                    * quadrature pts */

      for (int ip = 0; ip < ip_total; ip++) {
        dbl xi[3];
        dbl s, t, u;

        find_stu(ip, ielem_type, &s, &t, &u);
        xi[0] = s;
        xi[1] = t;
        xi[2] = u;

        /*
         * find quadrature weights for current ip
         */
        dbl wt = Gq_weight(ip, ielem_type);
        fv->wt = wt;

        /*
         * Load up basis function information for ea variable...
         * Old usage: fill_shape
         */
        err = load_basis_functions(xi, bfd);
        GOMA_EH(err, "problem from load_basis_functions");

        /*
         * This has elemental Jacobian transformation and some
         * basic mesh derivatives...
         * Old usage: calc_Jac, jelly_belly
         */
        err = beer_belly();
        GOMA_EH(err, "beer_belly");

        err = load_fv();
        GOMA_EH(err, "load_fv");

        err = load_bf_grad();
        GOMA_EH(err, "load_bf_grad");

        err = load_fv_grads();
        GOMA_EH(err, "load_fv_grads");

        double gamma[DIM][DIM];
        for (int a = 0; a < VIM; a++) {
          for (int b = 0; b < VIM; b++) {
            gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
          }
        }
        double mu = viscosity(gn, gamma, NULL);
        double invmu = (1.0 / mu);
        int eqn = PRESSURE;
        for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
          int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          int ldof_i = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][i];
          for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            int ldof_j = ei[upd->matrix_index[eqn]]->ln_to_dof[eqn][j];
            if (ldof_i >= 0 && ldof_j >= 0 &&
                ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
              dbl mm_contrib = 0.5 * invmu * bf[eqn]->phi[ldof_i] * bf[eqn]->phi[ldof_j] * fv->wt *
                               bf[eqn]->detJ;
              PetscInt global_row = matrix_data->schur_s_local_to_global[gnn_i];
              PetscInt global_col = matrix_data->schur_s_local_to_global[gnn_j];
              GOMA_EH(global_row, "Row not found");
              GOMA_EH(global_col, "Col not found");
              PetscErrorCode err =
                  MatSetValue(matrix_data->SchurS, global_row, global_col, mm_contrib, ADD_VALUES);
              CHKERRQ(err);
            }
          }
        }
      } /* END  for (ip = 0; ip < ip_total; ip++)                      */
    }   /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }     /* END for (ieb loop) */

  return GOMA_SUCCESS;
}

static goma_error initialize_petsc_post_proc_matrix(Exo_DB *exo,
                                                    Dpi *dpi,
                                                    dbl *x,
                                                    dbl *x_old,
                                                    dbl *xdot,
                                                    dbl *xdot_old,
                                                    PetscInt local_nodes,
                                                    PetscInt global_nodes) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)upd->petsc_post_proc_data;
  PetscInt global_offset;
  MPI_Scan(&local_nodes, &global_offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
  global_offset -= local_nodes;
  PetscInt num_rows = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  PetscInt num_cols = dpi->num_internal_nodes + dpi->num_boundary_nodes + dpi->num_external_nodes;

  PetscInt *d_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt *o_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));

  dbl *dblColGIDs = (dbl *)calloc(num_cols, sizeof(dbl));
  for (int i = 0; i < num_cols; i++) {
    dblColGIDs[i] = (double)global_offset + i;
  }
  matrix_data->local_to_global = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);

  exchange_node(cx[0], dpi, dblColGIDs);
  // convert back to Int
  for (int i = 0; i < num_cols; i++) {
    matrix_data->local_to_global[i] = (PetscInt)dblColGIDs[i];
  }
  free(dblColGIDs);
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

      int eqn = pd->ProjectionVar;
      if (bf[eqn]->interpolation == I_N1) {
        for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
          int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            if (gnn_i < num_rows) {
              if (gnn_j >= num_rows) {
                o_nnz[gnn_i] += 1;
              } else {
                d_nnz[gnn_i] += 1;
              }
            }
          }
        }

      } else {
        for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
          int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
          int ldof_i = ei[upd->matrix_index[pd->ProjectionVar]]->ln_to_dof[eqn][i];
          for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            int ldof_j = ei[upd->matrix_index[pd->ProjectionVar]]->ln_to_dof[eqn][j];
            if (ldof_i >= 0 && ldof_j >= 0 &&
                ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
              if (gnn_i < num_rows) {
                if (gnn_j >= num_rows) {
                  o_nnz[gnn_i] += 1;
                } else {
                  d_nnz[gnn_i] += 1;
                }
              }
            }
          }
        }
      }
    } /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }   /* END for (ieb loop) */

  if (Num_Proc == 1) {
    MatSeqAIJSetPreallocation(matrix_data->mat, 0, d_nnz);
  } else {
    MatMPIAIJSetPreallocation(matrix_data->mat, 0, d_nnz, 0, o_nnz);
  }

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
      int ielem_type = ei[pg->imtrx]->ielem_type;
      int ip_total = elem_info(NQUAD, ielem_type); /* number of
                                                    * quadrature pts */

      for (int ip = 0; ip < ip_total; ip++) {
        dbl xi[3];
        dbl s, t, u;

        find_stu(ip, ielem_type, &s, &t, &u);
        xi[0] = s;
        xi[1] = t;
        xi[2] = u;

        /*
         * find quadrature weights for current ip
         */
        dbl wt = Gq_weight(ip, ielem_type);
        fv->wt = wt;

        /*
         * Load up basis function information for ea variable...
         * Old usage: fill_shape
         */
        err = load_basis_functions(xi, bfd);
        GOMA_EH(err, "problem from load_basis_functions");

        /*
         * This has elemental Jacobian transformation and some
         * basic mesh derivatives...
         * Old usage: calc_Jac, jelly_belly
         */
        err = beer_belly();
        GOMA_EH(err, "beer_belly");

        err = load_coordinate_scales(pd->CoordinateSystem, fv);
        GOMA_EH(err, "load_coordinate_scales(fv)");

        err = load_bf_grad();
        GOMA_EH(err, "load_bf_grad");

        int eqn = pd->ProjectionVar;
        if (bf[eqn]->interpolation == I_N1) {
          for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
            int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
            for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
              int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
              dbl mm_contrib = bf[eqn]->phi[i] * bf[eqn]->phi[j] * fv->wt * bf[eqn]->detJ;
              PetscInt global_row = matrix_data->local_to_global[gnn_i];
              PetscInt global_col = matrix_data->local_to_global[gnn_j];
              PetscErrorCode err =
                  MatSetValue(matrix_data->mat, global_row, global_col, mm_contrib, ADD_VALUES);
              CHKERRQ(err);
            }
          }
        } else {
          for (int i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
            int gnn_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
            int ldof_i = ei[upd->matrix_index[pd->ProjectionVar]]->ln_to_dof[eqn][i];
            for (int j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
              int gnn_j = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
              int ldof_j = ei[upd->matrix_index[pd->ProjectionVar]]->ln_to_dof[eqn][j];
              if (ldof_i >= 0 && ldof_j >= 0 &&
                  ((gnn_i < dpi->num_internal_nodes + dpi->num_boundary_nodes) || Num_Proc == 1)) {
                dbl mm_contrib =
                    bf[eqn]->phi[ldof_i] * bf[eqn]->phi[ldof_j] * fv->wt * bf[eqn]->detJ;
                PetscInt global_row = matrix_data->local_to_global[gnn_i];
                PetscInt global_col = matrix_data->local_to_global[gnn_j];
                PetscErrorCode err =
                    MatSetValue(matrix_data->mat, global_row, global_col, mm_contrib, ADD_VALUES);
                CHKERRQ(err);
              }
            }
          }
        }
      } /* END  for (ip = 0; ip < ip_total; ip++)                      */
    }   /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }     /* END for (ieb loop) */

  free(d_nnz);
  free(o_nnz);
  return GOMA_SUCCESS;
}

goma_error goma_setup_petsc_post_proc_matrix(
    Exo_DB *exo, Dpi *dpi, dbl *x, dbl *x_old, dbl *xdot, dbl *xdot_old) {
  PetscBool petsc_initialized = PETSC_FALSE;
  PetscErrorCode err;
  PetscInitialized(&petsc_initialized);
  if (!petsc_initialized) {
    err = PetscInitializeNoArguments();
    CHKERRQ(err);
  }

  PetscMatrixData *matrix_data = malloc(sizeof(struct PetscMatrixData));

  upd->petsc_post_proc_data = (void *)matrix_data;

  if (GomaPetscOptions != NULL && !GomaPetscOptionsInserted) {
    err = PetscOptionsInsertString(NULL, GomaPetscOptions);
    CHKERRQ(err);
    GomaPetscOptionsInserted = true;
  }

  if (pg->imtrx == 0 && !GomaPetscOptionsPrinted) {
    CHKERRQ(PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD));
    GomaPetscOptionsPrinted = true;
  }

  err = MatCreate(MPI_COMM_WORLD, &matrix_data->mat);
  CHKERRQ(err);
  err = VecCreate(MPI_COMM_WORLD, &matrix_data->residual);
  CHKERRQ(err);
  err = VecCreate(MPI_COMM_WORLD, &matrix_data->update);
  CHKERRQ(err);
  err = KSPCreate(MPI_COMM_WORLD, &matrix_data->ksp);
  CHKERRQ(err);

  PC pc;
  err = KSPGetPC(matrix_data->ksp, &pc);
  CHKERRQ(err);
  err = PCSetType(pc, PCBJACOBI);
  CHKERRQ(err);
  err = PCSetReusePreconditioner(pc, PETSC_TRUE);
  CHKERRQ(err);
  err = KSPSetType(matrix_data->ksp, KSPGMRES);
  CHKERRQ(err);
  err = KSPSetTolerances(matrix_data->ksp, 1e-14, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(err);

  char ksp_prefix[25];
  snprintf(ksp_prefix, 24, "pp_mm_");
  err = MatSetOptionsPrefix(matrix_data->mat, ksp_prefix);
  CHKERRQ(err);
  snprintf(ksp_prefix, 24, "pp_");
  err = KSPSetOptionsPrefix(matrix_data->ksp, ksp_prefix);
  CHKERRQ(err);
  snprintf(ksp_prefix, 24, "pp_f_");
  err = VecSetOptionsPrefix(matrix_data->residual, ksp_prefix);
  CHKERRQ(err);
  snprintf(ksp_prefix, 24, "pp_sol_");
  err = VecSetOptionsPrefix(matrix_data->update, ksp_prefix);
  CHKERRQ(err);

  // find universe nodes
  PetscInt global_n_dof;
  PetscInt local_dof = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  MPI_Allreduce(&local_dof, &global_n_dof, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);

  err = VecSetSizes(matrix_data->update, local_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetSizes(matrix_data->residual, local_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->residual);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->update);
  CHKERRQ(err);
  err = MatSetSizes(matrix_data->mat, local_dof, local_dof, global_n_dof, global_n_dof);
  CHKERRQ(err);
  err = MatSetFromOptions(matrix_data->mat);
  CHKERRQ(err);
  initialize_petsc_post_proc_matrix(exo, dpi, x, x_old, xdot, xdot_old, local_dof, global_n_dof);
  err = KSPSetOperators(matrix_data->ksp, matrix_data->mat, matrix_data->mat);
  CHKERRQ(err);
  err = KSPSetFromOptions(matrix_data->ksp);
  CHKERRQ(err);

  matrix_data->matrix_setup = PETSC_FALSE;
  return GOMA_SUCCESS;
}

static goma_error initialize_petsc_matrix(struct GomaLinearSolverData *ams,
                                          Exo_DB *exo,
                                          Dpi *dpi,
                                          int internal_dof,
                                          int boundary_dof,
                                          int external_dof,
                                          int imtrx) {
  PetscInt num_rows = internal_dof + boundary_dof;
  PetscInt num_cols = internal_dof + boundary_dof + external_dof;
  double *dblColGIDs = malloc(sizeof(double) * num_cols);
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  matrix_data->local_to_global = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);
  PetscInt *d_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt *o_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt global_offset;
  PetscInt local_dof = internal_dof + boundary_dof;
  PetscInt nnz = 0;
  MPI_Scan(&local_dof, &global_offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
  global_offset -= local_dof;
  // copy global id's and convert to double for boundary exchange
  for (int i = 0; i < num_cols; i++) {
    dblColGIDs[i] = (double)global_offset + i;
  }
  exchange_dof(cx[pg->imtrx], dpi, dblColGIDs, pg->imtrx);
  // convert back to Int
  for (int i = 0; i < num_cols; i++) {
    matrix_data->local_to_global[i] = (PetscInt)dblColGIDs[i];
  }
  free(dblColGIDs);

  int total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;

  PetscInt irow_index = 0;

  /*
   * loop over all of the nodes on this processor
   */
  for (int inode = 0; inode < total_nodes; inode++) {
    NODAL_VARS_STRUCT *nv = Nodes[inode]->Nodal_Vars_Info[pg->imtrx];
    /*
     * Fill the vector list which points to the unknowns defined at this
     * node...
     */
    int inode_varType[MaxVarPerNode], inode_matID[MaxVarPerNode];
    int row_num_unknowns = fill_variable_vector(inode, inode_varType, inode_matID);
    /*
     * Do a check against the number of unknowns at this
     * node stored in the global array
     */
    if (row_num_unknowns != nv->Num_Unknowns) {
      GOMA_EH(GOMA_ERROR, "Inconsistency counting unknowns.");
    }

    /*
     * Loop over the unknowns defined at this row node
     */
    for (int iunknown = 0; iunknown < row_num_unknowns; iunknown++) {
      PetscInt row_nnz = 0;
      PetscInt off_nnz = 0;
      /*
       * Retrieve the var type of the current unknown
       */
      int rowVarType = inode_varType[iunknown];

      /*
       * Loop over the nodes which are determined to have an interaction
       * with the current row node
       */
      for (int j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1]; j++) {
        int inter_node = exo->node_node_list[j];
        NODE_INFO_STRUCT *nodeCol = Nodes[inter_node];
        NODAL_VARS_STRUCT *nvCol = nodeCol->Nodal_Vars_Info[pg->imtrx];

        /*
         * fill the vector list which points to the unknowns
         * defined at this interaction node
         */
        int inter_node_varType[MaxVarPerNode], inter_node_matID[MaxVarPerNode];
        int col_num_unknowns =
            fill_variable_vector(inter_node, inter_node_varType, inter_node_matID);
        if (col_num_unknowns != nvCol->Num_Unknowns) {
          GOMA_EH(GOMA_ERROR, "Inconsistency counting unknowns.");
        }

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be an interaction.
         */

        for (int inter_unknown = 0; inter_unknown < col_num_unknowns; inter_unknown++) {

          /*
           * HKM ->
           *  The entire process below is designed to find out
           *  what unknown we are processing. This coding is very
           *  convoluted and fraught with pitfalls. It should
           *  be rewritten using the structured approach in MPSalsa
           *  to variable indentification.
           */
          int colVarType = inter_node_varType[inter_unknown];

          /*
           * Query the Interaction mask to determine if a jacobian entry
           * should be created
           */
          int add_var = Inter_Mask[pg->imtrx][rowVarType][colVarType];

          /* The following code should be activated when solving DG viscoelastic problems
           * with full Jacobian treatment of upwind element stress terms
           */
          if (exo->centroid_list[inode] != -1 && inode != inter_node &&
              exo->centroid_list[inter_node] != -1) {
            int eb1 = exo->elem_eb[exo->centroid_list[inode]];

            if (vn_glob[Matilda[eb1]]->dg_J_model == FULL_DG) {
              int i1 = pd_glob[Matilda[eb1]]->i[pg->imtrx][rowVarType];
              int i2 = pd_glob[Matilda[eb1]]->i[pg->imtrx][colVarType];

              if ((rowVarType == colVarType) &&
                  (i1 == I_P0 || i1 == I_P1 || i1 == I_PQ1 || i1 == I_PQ2) &&
                  (i2 == I_P0 || i2 == I_P1 || i2 == I_PQ1 || i2 == I_PQ2) &&
                  (rowVarType != PRESSURE) &&
                  (rowVarType > VELOCITY_GRADIENT33 || rowVarType < VELOCITY_GRADIENT11)) {
                add_var = Inter_Mask[pg->imtrx][rowVarType][colVarType];
              } else {
                add_var = 0;
              }
            }
          }

          if (Debug_Flag < 0)
            add_var = TRUE; /* add all vars for checking jacobian */

          if (add_var) {
            /*
             * Determine the equation number of the current unknown
             */
            int icol_index = nodeCol->First_Unknown[pg->imtrx] + inter_unknown;
            if (icol_index < num_rows) {
              row_nnz++;
            } else {
              off_nnz++;
            }
          }
        }
      }
      // only do if we actually own the node for petsc
      if (irow_index < num_rows) {
        d_nnz[irow_index] += row_nnz;
        o_nnz[irow_index] += off_nnz;
      }
      irow_index++;
      nnz += row_nnz;
    }
  }

  if (Num_Proc == 1) {
    MatSeqAIJSetPreallocation(matrix_data->mat, 0, d_nnz);
  } else {
    MatMPIAIJSetPreallocation(matrix_data->mat, 0, d_nnz, 0, o_nnz);
  }

  /*
   * Add ams values that are needed elsewhere
   */
  ams->nnz = nnz;
  ams->nnz_plus = nnz;

  ams->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  ;
  ams->npn_plus = dpi->num_universe_nodes;

  ams->npu = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
  ams->npu_plus = num_universe_dofs[pg->imtrx];

  size_t local_dofs = num_universe_dofs[pg->imtrx];
  size_t global_dofs = num_universe_dofs[pg->imtrx];
  size_t local_nnz = nnz;
  size_t global_nnz = nnz;

  MPI_Allreduce(&local_dofs, &global_dofs, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_nnz, &global_nnz, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

  DPRINTF(stdout, "\n%-30s= %ld\n", "Number of unknowns", global_dofs);
  DPRINTF(stdout, "\n%-30s= %ld\n", "Number of matrix nonzeroes", global_nnz);

  free(d_nnz);
  free(o_nnz);

  return GOMA_SUCCESS;
}

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
                                   int imtrx) {
  PetscBool petsc_initialized = PETSC_FALSE;
  PetscErrorCode err;
  PetscInitialized(&petsc_initialized);
  if (!petsc_initialized) {
    err = PetscInitializeNoArguments();
    CHKERRQ(err);
  }

  PetscMatrixData *matrix_data = malloc(sizeof(struct PetscMatrixData));

  ams->PetscMatrixData = (void *)matrix_data;

  if (GomaPetscOptions != NULL && !GomaPetscOptionsInserted) {
    err = PetscOptionsInsertString(NULL, GomaPetscOptions);
    CHKERRQ(err);
    GomaPetscOptionsInserted = true;
  }

  PetscBool enable_log;
  PetscBool enable_log_set;
  PetscOptionsGetBool(NULL, NULL, "-enable_log", &enable_log, &enable_log_set);
  if (enable_log_set && enable_log) {
    PetscLogDefaultBegin();
  }
  PetscInt stokes_matrix;
  PetscBool stokes_matrix_set;
  PetscOptionsGetInt(NULL, NULL, "-stokes_matrix", &stokes_matrix, &stokes_matrix_set);

  PetscInt pcd_ns_remove_n = 0;
  PetscBool pcd_ns_remove_n_set;
  PetscOptionsGetInt(NULL, NULL, "-pcd_ns_remove_n", &pcd_ns_remove_n, &pcd_ns_remove_n_set);

  PetscInt pcd_ss_remove_n = 0;
  PetscBool pcd_ss_remove_n_set;
  PetscOptionsGetInt(NULL, NULL, "-pcd_ss_remove_n", &pcd_ss_remove_n, &pcd_ss_remove_n_set);

  PetscBool pcd_ns_remove_set;
  PetscInt *pcd_ns_remove = NULL;
  PetscInt pcd_ns_found = pcd_ns_remove_n;
  if (pcd_ns_remove_n_set) {
    PetscMalloc1(pcd_ns_remove_n, &pcd_ns_remove);
    PetscOptionsGetIntArray(NULL, NULL, "-pcd_ns_remove", pcd_ns_remove, &pcd_ns_found,
                            &pcd_ns_remove_set);
  }

  PetscBool pcd_ss_remove_set;
  PetscInt *pcd_ss_remove = NULL;
  PetscInt pcd_ss_found = pcd_ss_remove_n;
  if (pcd_ss_remove_n_set) {
    PetscMalloc1(pcd_ss_remove_n, &pcd_ss_remove);
    PetscOptionsGetIntArray(NULL, NULL, "-pcd_ss_remove", pcd_ss_remove, &pcd_ss_found,
                            &pcd_ss_remove_set);
  }

  matrix_data->pcd_ss_remove_n = pcd_ss_remove_n;
  matrix_data->pcd_ss_remove = pcd_ss_remove;
  matrix_data->pcd_ns_remove_n = pcd_ns_remove_n;
  matrix_data->pcd_ns_remove = pcd_ns_remove;
  PetscBool user_schur;
  PetscBool user_schur_set;
  PetscOptionsGetBool(NULL, NULL, "-user_schur", &user_schur, &user_schur_set);

  PetscBool user_pcd;
  PetscBool user_pcd_set;
  PetscOptionsGetBool(NULL, NULL, "-user_pcd", &user_pcd, &user_pcd_set);

  PetscBool user_pcd_inverse_diag;
  PetscBool user_pcd_inverse_diag_set;
  PetscOptionsGetBool(NULL, NULL, "-user_pcd_inverse_diag", &user_pcd_inverse_diag,
                      &user_pcd_inverse_diag_set);

  if (imtrx == 0 && !GomaPetscOptionsPrinted) {
    CHKERRQ(PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD));
    GomaPetscOptionsPrinted = true;
  }

  err = MatCreate(MPI_COMM_WORLD, &matrix_data->mat);
  CHKERRQ(err);
  err = VecCreate(MPI_COMM_WORLD, &matrix_data->residual);
  CHKERRQ(err);
  err = VecCreate(MPI_COMM_WORLD, &matrix_data->update);
  CHKERRQ(err);
  err = KSPCreate(MPI_COMM_WORLD, &matrix_data->ksp);
  CHKERRQ(err);

  char ksp_prefix[25];
  if (upd->Total_Num_Matrices > 1) {
    snprintf(ksp_prefix, 24, "mat%d_", imtrx);
    err = MatSetOptionsPrefix(matrix_data->mat, ksp_prefix);
    CHKERRQ(err);
    snprintf(ksp_prefix, 24, "sys%d_", imtrx);
    err = KSPSetOptionsPrefix(matrix_data->ksp, ksp_prefix);
    CHKERRQ(err);
    snprintf(ksp_prefix, 24, "res%d_", imtrx);
    err = VecSetOptionsPrefix(matrix_data->residual, ksp_prefix);
    CHKERRQ(err);
    snprintf(ksp_prefix, 24, "upd%d_", pg->imtrx);
    err = VecSetOptionsPrefix(matrix_data->update, ksp_prefix);
    CHKERRQ(err);
  } else {
    snprintf(ksp_prefix, 24, "res_");
    err = VecSetOptionsPrefix(matrix_data->residual, ksp_prefix);
    CHKERRQ(err);
    snprintf(ksp_prefix, 24, "upd_");
    err = VecSetOptionsPrefix(matrix_data->update, ksp_prefix);
    CHKERRQ(err);
  }

  // find universe dof
  PetscInt global_n_dof;
  PetscInt local_dof = internal_dof + boundary_dof;
  MPI_Allreduce(&local_dof, &global_n_dof, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);

  err = VecSetSizes(matrix_data->update, internal_dof + boundary_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetSizes(matrix_data->residual, internal_dof + boundary_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->residual);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->update);
  CHKERRQ(err);
  err = MatSetSizes(matrix_data->mat, internal_dof + boundary_dof, internal_dof + boundary_dof,
                    global_n_dof, global_n_dof);
  CHKERRQ(err);
  err = MatSetFromOptions(matrix_data->mat);
  CHKERRQ(err);
  initialize_petsc_matrix(ams, exo, dpi, internal_dof, boundary_dof, external_dof, imtrx);
  err = KSPSetOperators(matrix_data->ksp, matrix_data->mat, matrix_data->mat);
  CHKERRQ(err);
  err = KSPSetFromOptions(matrix_data->ksp);
  CHKERRQ(err);

  if (stokes_matrix_set && stokes_matrix == imtrx) {
    GOMA_WH(GOMA_ERROR, "stokes matrix set to %ld", stokes_matrix, pd_glob[0]->Num_Dim);
    IS ufields;
    IS pfields;

    PetscInt vel_n, pres_n;
    PetscInt *vel_is, *pres_is;
    get_pressure_velocity_is(matrix_data, exo, dpi, x, x_old, xdot, xdot_old, &vel_n, &vel_is,
                             &pres_n, &pres_is);

    ISCreateGeneral(PETSC_COMM_WORLD, vel_n, vel_is, PETSC_OWN_POINTER, &ufields);
    ISCreateGeneral(PETSC_COMM_WORLD, pres_n, pres_is, PETSC_OWN_POINTER, &pfields);
    PC pc;
    KSPGetPC(matrix_data->ksp, &pc);
    PCFieldSplitSetIS(pc, "u", ufields);
    PCFieldSplitSetIS(pc, "p", pfields);

    // if (pd_glob[0]->Num_Dim == 3) {
    //   PC pc;
    //   const PetscInt ufields[] = {0, 1, 2}, pfields[] = {3};
    //   KSPGetPC(matrix_data->ksp, &pc);
    //   PCFieldSplitSetBlockSize(pc, 4);
    //   PCFieldSplitSetFields(pc, "u", 3, ufields, ufields);
    //   PCFieldSplitSetFields(pc, "p", 1, pfields, pfields);
    // } else if (pd_glob[0]->Num_Dim == 2) {
    //   PC pc;
    //   const PetscInt ufields[] = {0, 1}, pfields[] = {3};
    //   KSPGetPC(matrix_data->ksp, &pc);
    //   PCFieldSplitSetBlockSize(pc, 3);
    //   PCFieldSplitSetFields(pc, "u", 2, ufields, ufields);
    //   PCFieldSplitSetFields(pc, "p", 1, pfields, pfields);
    // }

    PetscInt local_pressure_nodes = 0, global_pressure_nodes = 0;
    if ((user_schur_set && user_schur) || (user_pcd_set && user_pcd)) {
      count_pressure_nodes(matrix_data, exo, dpi, x, x_old, xdot, xdot_old, &local_pressure_nodes,
                           &global_pressure_nodes);
    }
    matrix_data->user_schur = PETSC_FALSE;
    if (user_schur_set && user_schur) {
      matrix_data->user_schur = user_schur;
      MatCreate(MPI_COMM_WORLD, &matrix_data->SchurS);
      MatSetOptionsPrefix(matrix_data->SchurS, "schur_s_");
      MatSetFromOptions(matrix_data->SchurS);
      err = MatSetSizes(matrix_data->SchurS, local_pressure_nodes, local_pressure_nodes,
                        global_pressure_nodes, global_pressure_nodes);
      CHKERRQ(err);
      PC pc;
      KSPGetPC(matrix_data->ksp, &pc);
      initialize_petsc_pressure_matrix(matrix_data, exo, dpi, x, x_old, xdot, xdot_old,
                                       local_pressure_nodes, global_pressure_nodes);
      PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, matrix_data->SchurS);
    }

    matrix_data->user_pcd = PETSC_FALSE;
    if (user_pcd_set && user_pcd) {
      matrix_data->user_pcd = user_pcd;
      matrix_data->user_pcd_inverse_diag = PETSC_FALSE;
      if (user_pcd_inverse_diag_set) {
        matrix_data->user_pcd_inverse_diag = user_pcd_inverse_diag;
      }
      // PC pc;
      // err = KSPGetPC(matrix_data->ksp, &pc);
      // CHKERRQ(err);
      // petsc_PCD_setup(pc, matrix_data, exo, dpi, x, x_old, xdot, xdot_old);
    }

    matrix_data->matrix_setup = PETSC_FALSE;
  } else {
    matrix_data->user_schur = PETSC_FALSE;
    matrix_data->user_pcd = PETSC_FALSE;
  }
  return GOMA_SUCCESS;
}

void petsc_load_lec(int ielem, struct GomaLinearSolverData *ams, double resid_vector[]) {
  int e, v, i, j, pe, pv;
  int dofs;
  int gnn, row_index, ke, kv, nvdof;
  int col_index, ledof;
  int je_new;
  struct Element_Indices *ei_ptr;

  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  for (e = V_FIRST; e < V_LAST; e++) {
    pe = upd->ep[pg->imtrx][e];
    if (pe != -1) {
      if (e == R_MASS) {
        for (ke = 0; ke < upd->Max_Num_Species_Eqn; ke++) {
          pe = MAX_PROB_VAR + ke;
          dofs = ei[pg->imtrx]->dof[e];
          for (i = 0; i < dofs; i++) {
            /*
             * Check to see whether this dof is on a
             * node owned by this processor. If it isn't
             * there is no need to fill in the processor
             * residual and Jacobian row for this variable.
             */
            ledof = ei[pg->imtrx]->lvdof_to_ledof[e][i];
            if (ei[pg->imtrx]->owned_ledof[ledof]) {
              gnn = ei[pg->imtrx]->gnn_list[e][i];
              nvdof = ei[pg->imtrx]->Baby_Dolphin[e][i];
              je_new = ei[pg->imtrx]->ieqn_ledof[ledof] + ke;
              row_index =
                  Index_Solution(gnn, e, ke, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
              resid_vector[row_index] += lec->R[LEC_R_INDEX(MAX_PROB_VAR + ke, i)];

              if (af->Assemble_Jacobian) {
                for (v = V_FIRST; v < V_LAST; v++) {
                  pv = upd->vp[pg->imtrx][v];
                  if (pv != -1 && (Inter_Mask[pg->imtrx][e][v])) {
                    ei_ptr = ei[pg->imtrx];
                    if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                      if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                        ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                        if (ei_ptr == 0) {
                          printf("ei_ptr == 0\n");
                          exit(-1);
                        }
                      }
                    }
                    if (v == MASS_FRACTION) {
                      for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                        pv = MAX_PROB_VAR + kv;
                        for (j = 0; j < ei_ptr->dof[v]; j++) {
                          ledof = ei_ptr->lvdof_to_ledof[v][j];
                          je_new = ei_ptr->ieqn_ledof[ledof] + kv;
                          col_index = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                     ei_ptr->Baby_Dolphin[v][j],
                                                     ei_ptr->matID_ledof[ledof], pg->imtrx);
                          GOMA_EH(col_index, "Bad var index.");

                          PetscInt global_row = matrix_data->local_to_global[row_index];
                          PetscInt global_col = matrix_data->local_to_global[col_index];
                          MatSetValue(matrix_data->mat, global_row, global_col,
                                      lec->J[LEC_J_INDEX(pe, pv, i, j)], ADD_VALUES);
                        }
                      }
                    } else {
                      pv = upd->vp[pg->imtrx][v];
                      kv = 0;
                      for (j = 0; j < ei_ptr->dof[v]; j++) {
                        ledof = ei_ptr->lvdof_to_ledof[v][j];
                        je_new = ei_ptr->ieqn_ledof[ledof];
                        col_index = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                   ei_ptr->Baby_Dolphin[v][j],
                                                   ei_ptr->matID_ledof[ledof], pg->imtrx);
                        if (col_index != je_new) {
                          fprintf(stderr, "Oh fiddlesticks: je = %d, je_new = %d\n", col_index,
                                  je_new);
                          GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                        }
                        GOMA_EH(col_index, "Bad var index.");
                        PetscInt global_row = matrix_data->local_to_global[row_index];
                        PetscInt global_col = matrix_data->local_to_global[col_index];
                        MatSetValue(matrix_data->mat, global_row, global_col,
                                    lec->J[LEC_J_INDEX(pe, pv, i, j)], ADD_VALUES);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        pe = upd->ep[pg->imtrx][e];
        dofs = ei[pg->imtrx]->dof[e];
        for (i = 0; i < dofs; i++) {
          /*
           * Check to see whether this dof is on a
           * node owned by this processor. If it isn't
           * there is no need to fill in the processor
           * residual and Jacobian row for this variable.
           */
          ledof = ei[pg->imtrx]->lvdof_to_ledof[e][i];
          if (ei[pg->imtrx]->owned_ledof[ledof]) {
            row_index = ei[pg->imtrx]->gun_list[e][i];
            resid_vector[row_index] += lec->R[LEC_R_INDEX(pe, i)];

            if (af->Assemble_Jacobian) {
              for (v = V_FIRST; v < V_LAST; v++) {
                pv = upd->vp[pg->imtrx][v];
                if (pv != -1 && (Inter_Mask[pg->imtrx][e][v])) {
                  if (v == MASS_FRACTION) {

                    ei_ptr = ei[pg->imtrx];
                    if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                      if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                        ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                        if (ei_ptr == 0) {
                          GOMA_EH(GOMA_ERROR, "ei_ptr == 0\n");
                          exit(-1);
                        }
                      }
                    }
                    for (kv = 0; kv < upd->Max_Num_Species_Eqn; kv++) {
                      pv = MAX_PROB_VAR + kv;
                      for (j = 0; j < ei_ptr->dof[v]; j++) {
                        ledof = ei_ptr->lvdof_to_ledof[v][j];
                        je_new = ei_ptr->ieqn_ledof[ledof] + kv;
                        col_index = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                   ei_ptr->Baby_Dolphin[v][j],
                                                   ei_ptr->matID_ledof[ledof], pg->imtrx);
                        if (col_index != je_new) {
                          /*
                           * HKM -> another special case. Until we delineate
                           *        the subspecies as individual local
                           *        dofs, je_new will be wrong here. je
                           *        is correct.
                           */
                          if (Nodes[ei[pg->imtrx]->gnn_list[v][j]]->Mat_List.Length < 2) {
                          }
                        }
                        GOMA_EH(col_index, "Bad var index.");
                        PetscInt global_row = matrix_data->local_to_global[row_index];
                        PetscInt global_col = matrix_data->local_to_global[col_index];
                        MatSetValue(matrix_data->mat, global_row, global_col,
                                    lec->J[LEC_J_INDEX(pe, pv, i, j)], ADD_VALUES);
                      }
                    }
                  } else {
                    kv = 0;
                    // NEW IMPLEMENTATION
                    ei_ptr = ei[pg->imtrx];
                    if (ei[pg->imtrx]->owningElementForColVar[v] != ielem) {
                      if (ei[pg->imtrx]->owningElementForColVar[v] != -1) {
                        ei_ptr = ei[pg->imtrx]->owningElement_ei_ptr[v];
                        if (ei_ptr == 0) {
                          GOMA_EH(GOMA_ERROR, "ei slave pointer is null");
                          exit(-1);
                        }
                      }
                    }
                    for (j = 0; j < ei_ptr->dof[v]; j++) {
                      ledof = ei_ptr->lvdof_to_ledof[v][j];
                      je_new = ei_ptr->ieqn_ledof[ledof];
                      col_index =
                          Index_Solution(ei_ptr->gnn_list[v][j], v, kv, ei_ptr->Baby_Dolphin[v][j],
                                         ei_ptr->matID_ledof[ledof], pg->imtrx);
                      if (col_index != je_new) {
                        fprintf(stderr, "Oh fiddlesticks: je = %d, je_new = %d\n", col_index,
                                je_new);
                        GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                      }
                      GOMA_EH(col_index, "Bad var index.");
                      PetscInt global_row = matrix_data->local_to_global[row_index];
                      PetscInt global_col = matrix_data->local_to_global[col_index];
                      MatSetValue(matrix_data->mat, global_row, global_col,
                                  lec->J[LEC_J_INDEX(pe, pv, i, j)], ADD_VALUES);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

goma_error petsc_scale_matrix(struct GomaLinearSolverData *ams, double *b_, double *scale) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  Vec row_sums;

  MatAssemblyBegin(matrix_data->mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->mat, MAT_FINAL_ASSEMBLY);

  CHKERRQ(MatCreateVecs(matrix_data->mat, &row_sums, NULL));
  VecZeroEntries(row_sums);
  for (int i = 0; i < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]; i++) {
    PetscInt ncols;
    const PetscScalar *vals;
    MatGetRow(matrix_data->mat, matrix_data->local_to_global[i], &ncols, NULL, &vals);
    dbl sum = 0.0;
    for (PetscInt c = 0; c < ncols; c++) {
      sum += PetscAbsScalar(vals[c]);
    }
    MatRestoreRow(matrix_data->mat, matrix_data->local_to_global[i], &ncols, NULL, &vals);
    scale[i] = 1.0 / sum;
  }
  VecSetValues(row_sums, num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx],
               matrix_data->local_to_global, scale, INSERT_VALUES);
  MatDiagonalScale(matrix_data->mat, row_sums, NULL);
  for (int i = 0; i < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]; i++) {
    b_[i] *= scale[i];
    scale[i] = 1.0 / scale[i];
  }
  exchange_dof(cx[pg->imtrx], DPI_ptr, b_, pg->imtrx);
  exchange_dof(cx[pg->imtrx], DPI_ptr, scale, pg->imtrx);
  CHKERRQ(VecDestroy(&row_sums));
  return GOMA_SUCCESS;
}

// vim: expandtab sw=2 ts=8
int petsc_solve(struct GomaLinearSolverData *ams, double *x_, double *b_, int *its) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  MatAssemblyBegin(matrix_data->mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->mat, MAT_FINAL_ASSEMBLY);
  if (matrix_data->user_schur) {
    MatZeroEntries(matrix_data->SchurS);
    set_petsc_pressure_matrix(matrix_data, EXO_ptr, DPI_ptr, pg->matrices[pg->imtrx].x,
                              pg->matrices[pg->imtrx].x_old, pg->matrices[pg->imtrx].xdot,
                              pg->matrices[pg->imtrx].xdot_old);

    MatAssemblyBegin(matrix_data->SchurS, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->SchurS, MAT_FINAL_ASSEMBLY);
  } else if (matrix_data->user_pcd) {
    if (!matrix_data->matrix_setup) {
      matrix_data->matrix_setup = PETSC_TRUE;
      PC pc;
      int err = KSPGetPC(matrix_data->ksp, &pc);
      CHKERRQ(err);
      petsc_PCD_setup(pc, matrix_data, EXO_ptr, DPI_ptr, pg->matrices[pg->imtrx].x,
                      pg->matrices[pg->imtrx].x_old, pg->matrices[pg->imtrx].xdot,
                      pg->matrices[pg->imtrx].xdot_old);
    }
    MatZeroEntries(matrix_data->pcd_data->Fp);
    MatZeroEntries(matrix_data->pcd_data->Mp);
    MatZeroEntries(matrix_data->pcd_data->Mp_mu);
    MatZeroEntries(matrix_data->pcd_data->Ap);
    set_pcd_matrices(matrix_data, EXO_ptr, DPI_ptr, pg->matrices[pg->imtrx].x,
                     pg->matrices[pg->imtrx].x_old, pg->matrices[pg->imtrx].xdot,
                     pg->matrices[pg->imtrx].xdot_old);
    MatAssemblyBegin(matrix_data->pcd_data->Fp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->pcd_data->Fp, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(matrix_data->pcd_data->Mp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->pcd_data->Mp, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(matrix_data->pcd_data->Mp_mu, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->pcd_data->Mp_mu, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(matrix_data->pcd_data->Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->pcd_data->Ap, MAT_FINAL_ASSEMBLY);
  }
  for (PetscInt i = 0; i < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]; i++) {
    VecSetValue(matrix_data->residual, matrix_data->local_to_global[i], b_[i], INSERT_VALUES);
    VecSetValue(matrix_data->update, matrix_data->local_to_global[i], x_[i], INSERT_VALUES);
  }
  VecAssemblyBegin(matrix_data->residual);
  VecAssemblyEnd(matrix_data->residual);
  VecAssemblyBegin(matrix_data->update);
  VecAssemblyEnd(matrix_data->update);

  KSPSolve(matrix_data->ksp, matrix_data->residual, matrix_data->update);
  PetscInt pits;
  KSPGetIterationNumber(matrix_data->ksp, &pits);
  *its = pits;
  VecGetValues(matrix_data->update, num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx],
               matrix_data->local_to_global, x_);
  MatZeroEntries(matrix_data->mat);
  return 0;
}

void petsc_solve_post_proc(double **post_proc_vect, RESULTS_DESCRIPTION_STRUCT *rd, Dpi *dpi) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)upd->petsc_post_proc_data;
  if (!matrix_data->matrix_setup) {
    MatAssemblyBegin(matrix_data->mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_data->mat, MAT_FINAL_ASSEMBLY);
    matrix_data->matrix_setup = PETSC_TRUE;
  }
  P0PRINTF("\nPETSc Post Processing Projection:\n");
  double start = MPI_Wtime();
  for (int pp = 0; pp < rd->TotalNVPostOutput; pp++) {
    double pp_start = MPI_Wtime();
    exchange_node(cx[0], dpi, post_proc_vect[pp]);
    for (PetscInt i = 0; i < dpi->num_internal_nodes + dpi->num_boundary_nodes; i++) {
      VecSetValue(matrix_data->residual, matrix_data->local_to_global[i], post_proc_vect[pp][i],
                  INSERT_VALUES);
      VecSetValue(matrix_data->update, matrix_data->local_to_global[i], 0.0, INSERT_VALUES);
    }
    VecAssemblyBegin(matrix_data->residual);
    VecAssemblyEnd(matrix_data->residual);
    VecAssemblyBegin(matrix_data->update);
    VecAssemblyEnd(matrix_data->update);

    KSPSolve(matrix_data->ksp, matrix_data->residual, matrix_data->update);
    PetscInt pits;
    KSPGetIterationNumber(matrix_data->ksp, &pits);
    VecGetValues(matrix_data->update, dpi->num_internal_nodes + dpi->num_boundary_nodes,
                 matrix_data->local_to_global, post_proc_vect[pp]);
    double pp_end = MPI_Wtime();
    P0PRINTF("PP %*.*s %ld iterations %4.4e seconds\n", 10, 10,
             rd->nvname[rd->TotalNVSolnOutput + pp], (long)pits, pp_end - pp_start);
  }
  double end = MPI_Wtime();
  P0PRINTF("Post Proc Time = %4.4e s\n\n", end - start);
}

goma_error goma_petsc_free_matrix(struct GomaLinearSolverData *ams) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  PetscErrorCode err;

  err = MatDestroy(&matrix_data->mat);
  CHKERRQ(err);
  err = VecDestroy(&matrix_data->residual);
  CHKERRQ(err);
  err = VecDestroy(&matrix_data->update);
  CHKERRQ(err);
  err = KSPDestroy(&matrix_data->ksp);
  CHKERRQ(err);

  return GOMA_SUCCESS;
}
#endif
#endif
