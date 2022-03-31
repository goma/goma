#include <mpi.h>
#include <petscmath.h>
#include <stddef.h>
#ifdef GOMA_ENABLE_PETSC
#include <petscsystypes.h>
#if (PETSC_USE_COMPLEX)
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>
#ifdef I
#undef I
#endif

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
#include "sl_petsc_complex.h"

static bool GomaPetscOptionsInserted = false;
static bool GomaPetscOptionsPrinted = false;

static bool goma_is_imag_eqn(int eqn) {
  switch (eqn) {
  case EM_E1_IMAG:
  case EM_E2_IMAG:
  case EM_E3_IMAG:
  case EM_H1_IMAG:
  case EM_H2_IMAG:
  case EM_H3_IMAG:
    return true;
  default:
    return false;
  }
}
static int goma_real_from_imag(int eqn) {
  switch (eqn) {
  case EM_E1_IMAG:
    return EM_E1_REAL;
  case EM_E2_IMAG:
    return EM_E2_REAL;
  case EM_E3_IMAG:
    return EM_E3_REAL;
  case EM_H1_IMAG:
    return EM_H1_REAL;
  case EM_H2_IMAG:
    return EM_H2_REAL;
  case EM_H3_IMAG:
    return EM_H3_REAL;
  default:
    GOMA_EH(GOMA_ERROR, "goma_real_from_imag given non-imaginary equation");
    return -1;
  }
}

static void create_complex_mapping(struct GomaLinearSolverData *ams,
                                   Exo_DB *exo,
                                   Dpi *dpi,
                                   int internal_dof,
                                   int boundary_dof,
                                   int external_dof,
                                   int imtrx) {
  PetscInt num_cols = internal_dof + boundary_dof + external_dof;
  double *dblColGIDs = malloc(sizeof(double) * num_cols);
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  matrix_data->local_to_global = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);
  matrix_data->imag_local_to_goma_local = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);
  matrix_data->real_local_to_goma_local = (PetscInt *)malloc(sizeof(PetscInt) * num_cols);
  PetscInt global_offset;
  PetscInt local_dof = internal_dof + boundary_dof;

  PetscBool *col_is_real = (PetscBool *)calloc(num_cols, sizeof(PetscBool));

  PetscInt irow_index = 0;
  PetscInt rrow_index = 0;

  /*
   * loop over all of the nodes on this processor
   */
  for (int inode = 0; inode < (Num_Internal_Nodes + Num_Border_Nodes); inode++) {
    NODAL_VARS_STRUCT *nv = Nodes[inode]->Nodal_Vars_Info[pg->imtrx];
    /*
     * Fill the vector list which points to the unknowns defined at this
     * node...
     */
    int inode_varType[MaxVarPerNode], inode_matID[MaxVarPerNode];
    int inode_varIndex[MaxVarPerNode];
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
      /*
       * Retrieve the var type of the current unknown
       */
      int rowVarType = inode_varType[iunknown];

      if (goma_is_imag_eqn(rowVarType)) {
        // create a mapping from the real_eqn

        int eqn = goma_real_from_imag(rowVarType);
        int offset = -1;
        for (int i = 0; i < iunknown; i++) {
          if (inode_varType[i] == eqn) {
            offset = i;
            break;
          }
        }
        GOMA_EH(offset, "Could not find matching real equation for imaginary equation");
        inode_varIndex[iunknown] = inode_varIndex[offset];
        matrix_data->imag_local_to_goma_local[inode_varIndex[offset]] = rrow_index;
      } else {
        inode_varIndex[iunknown] = irow_index;
        matrix_data->real_local_to_goma_local[irow_index] = rrow_index;
        col_is_real[rrow_index] = true;
        irow_index++;
      }
      rrow_index++;
    }
  }

  local_dof = irow_index;
  MPI_Scan(&local_dof, &global_offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
  global_offset -= local_dof;
  // copy global id's and convert to double for boundary exchange
  int local_offset = 0;
  for (int i = 0; i < num_cols; i++) {
    if (col_is_real[i]) {
      dblColGIDs[i] = (double)global_offset + local_offset;
      local_offset++;
    } else {
      dblColGIDs[i] = -1;
    }
  }
  exchange_dof(cx[pg->imtrx], dpi, dblColGIDs, pg->imtrx);
  matrix_data->local_to_global_complex = malloc(sizeof(PetscInt) * num_universe_dofs[pg->imtrx]);
  // convert back to Int
  int offset = 0;
  for (int i = 0; i < num_cols; i++) {
    PetscInt id = (PetscInt)dblColGIDs[i];
    matrix_data->local_to_global[i] = id;
    if (id > -1) {
      matrix_data->local_to_global_complex[offset] = id;
      offset++;
    }
  }
  matrix_data->local_dof = local_dof;
  PetscInt global_dof = 0;
  MPI_Allreduce(&local_dof, &global_dof, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
  matrix_data->global_dof = global_dof;
  matrix_data->is_imag = col_is_real;
  free(dblColGIDs);
}

static goma_error initialize_petsc_matrix_complex(struct GomaLinearSolverData *ams,
                                                  Exo_DB *exo,
                                                  Dpi *dpi,
                                                  int internal_dof,
                                                  int boundary_dof,
                                                  int external_dof,
                                                  int imtrx) {
  PetscInt num_rows = internal_dof + boundary_dof;
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  PetscInt *d_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt *o_nnz = (PetscInt *)calloc(num_rows, sizeof(PetscInt));
  PetscInt nnz = 0;

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

      if (!goma_is_imag_eqn(rowVarType)) {

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

            if (!goma_is_imag_eqn(colVarType)) {

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
  }

  matrix_data->d_nnz = d_nnz;
  matrix_data->o_nnz = o_nnz;

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

  // free(d_nnz);
  // free(o_nnz);

  return GOMA_SUCCESS;
}

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
                                           int imtrx) {
  PetscBool petsc_initialized = PETSC_FALSE;
  PetscErrorCode err;
  PetscInitialized(&petsc_initialized);
  if (!petsc_initialized) {
    err = PetscInitializeNoArguments();
    CHKERRQ(err);
  }

  PetscMatrixData *matrix_data = calloc(1, sizeof(struct PetscMatrixData));

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
  create_complex_mapping(ams, exo, dpi, internal_dof, boundary_dof, external_dof, imtrx);

  err = VecSetSizes(matrix_data->update, matrix_data->local_dof, matrix_data->global_dof);
  CHKERRQ(err);
  err = VecSetSizes(matrix_data->residual, matrix_data->local_dof, matrix_data->global_dof);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->residual);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->update);
  CHKERRQ(err);
  err = MatSetSizes(matrix_data->mat, matrix_data->local_dof, matrix_data->local_dof,
                    matrix_data->global_dof, matrix_data->global_dof);
  CHKERRQ(err);
  err = MatSetFromOptions(matrix_data->mat);
  CHKERRQ(err);
  initialize_petsc_matrix_complex(ams, exo, dpi, internal_dof, boundary_dof, external_dof, imtrx);
  err = KSPSetOperators(matrix_data->ksp, matrix_data->mat, matrix_data->mat);
  CHKERRQ(err);
  err = KSPSetFromOptions(matrix_data->ksp);
  CHKERRQ(err);

  return GOMA_SUCCESS;
}

void petsc_load_lec_complex(int ielem, struct GomaLinearSolverData *ams, double resid_vector[]) {
  int e, v, i, j, pe, pv;
  int dofs;
  int gnn, row_index, ke, kv, nvdof;
  int col_index, ledof;
  int je_new;
  struct Element_Indices *ei_ptr;
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  if (ielem == 0) {

    MatDestroy(&matrix_data->mat);
    MatCreate(MPI_COMM_WORLD, &matrix_data->mat);
    MatSetSizes(matrix_data->mat, matrix_data->local_dof, matrix_data->local_dof,
                matrix_data->global_dof, matrix_data->global_dof);
    MatSetFromOptions(matrix_data->mat);
    if (Num_Proc == 1) {
      MatSeqAIJSetPreallocation(matrix_data->mat, 0, matrix_data->d_nnz);
    } else {
      MatMPIAIJSetPreallocation(matrix_data->mat, 0, matrix_data->d_nnz, 0, matrix_data->o_nnz);
    }
  }

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
            ei_ptr = ei[pg->imtrx];
            row_index = ei[pg->imtrx]->gun_list[e][i];
            resid_vector[row_index] += lec->R[LEC_R_INDEX(pe, i)];
            if (goma_is_imag_eqn(e)) {
              row_index =
                  Index_Solution(ei_ptr->gnn_list[e][i], goma_real_from_imag(e), 0,
                                 ei_ptr->Baby_Dolphin[e][i], ei_ptr->matID_ledof[ledof], pg->imtrx);
            }

            if (af->Assemble_Jacobian) {
              for (v = V_FIRST; v < V_LAST; v++) {
                pv = upd->vp[pg->imtrx][v];
                // We don't add these as they are already in the matrix from imag(e)
                if (goma_is_imag_eqn(v))
                  continue;
                if (pv != -1 && (Inter_Mask[pg->imtrx][e][v])) {
                  if (v == MASS_FRACTION) {

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
                      PetscComplex va;
                      PetscInt global_row;
                      PetscInt global_col;
                      ledof = ei_ptr->lvdof_to_ledof[v][j];
                      je_new = ei_ptr->ieqn_ledof[ledof];
                      global_row = matrix_data->local_to_global[row_index];
                      if (goma_is_imag_eqn(v)) {
                        col_index = Index_Solution(ei_ptr->gnn_list[goma_real_from_imag(v)][j],
                                                   goma_real_from_imag(v), kv,
                                                   ei_ptr->Baby_Dolphin[goma_real_from_imag(v)][j],
                                                   ei_ptr->matID_ledof[ledof], pg->imtrx);
                        GOMA_EH(col_index, "Bad var index.");
                        global_col = matrix_data->local_to_global[col_index];
                      } else {
                        col_index = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                                                   ei_ptr->Baby_Dolphin[v][j],
                                                   ei_ptr->matID_ledof[ledof], pg->imtrx);
                        if (col_index != je_new) {
                          fprintf(stderr, "Oh fiddlesticks: je = %d, je_new = %d\n", col_index,
                                  je_new);
                          GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                        }
                        GOMA_EH(col_index, "Bad var index.");
                        global_col = matrix_data->local_to_global[col_index];
                      }
                      bool add_value = false;
                      if (goma_is_imag_eqn(e) && !goma_is_imag_eqn(v)) {
                        va = lec->J[LEC_J_INDEX(pe, pv, i, j)] * PETSC_i;
                        add_value = true;
                      } else if (!goma_is_imag_eqn(e) && !goma_is_imag_eqn(v)) {
                        va = lec->J[LEC_J_INDEX(pe, pv, i, j)];
                        add_value = true;
                      }
                      if (add_value)
                        MatSetValue(matrix_data->mat, global_row, global_col, va, ADD_VALUES);
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

int petsc_solve_complex(struct GomaLinearSolverData *ams, double *x_, double *b_, int *its) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  KSPDestroy(&matrix_data->ksp);
  KSPCreate(MPI_COMM_WORLD, &matrix_data->ksp);
  KSPSetFromOptions(matrix_data->ksp);
  KSPSetOperators(matrix_data->ksp, matrix_data->mat, matrix_data->mat);
  MatAssemblyBegin(matrix_data->mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->mat, MAT_FINAL_ASSEMBLY);
  VecSet(matrix_data->residual, 0.0);
  VecSet(matrix_data->update, 0.0);
  for (int i = 0; i < matrix_data->local_dof; i++) {
    int imag_offset = matrix_data->imag_local_to_goma_local[i];
    int real_offset = matrix_data->real_local_to_goma_local[i];
    if (imag_offset > -1) {
      PetscComplex rval = b_[real_offset] + PETSC_i * b_[imag_offset];
      PetscComplex xval = x_[real_offset] + PETSC_i * x_[imag_offset];
      VecSetValue(matrix_data->residual, matrix_data->local_to_global_complex[i], rval,
                  INSERT_VALUES);
      VecSetValue(matrix_data->update, matrix_data->local_to_global_complex[i], xval,
                  INSERT_VALUES);
    } else {
      VecSetValue(matrix_data->residual, matrix_data->local_to_global_complex[i], b_[real_offset],
                  INSERT_VALUES);
      VecSetValue(matrix_data->update, matrix_data->local_to_global_complex[i], x_[real_offset],
                  INSERT_VALUES);
    }
  }
  VecAssemblyBegin(matrix_data->residual);
  VecAssemblyEnd(matrix_data->residual);
  VecAssemblyBegin(matrix_data->update);
  VecAssemblyEnd(matrix_data->update);

  KSPSolve(matrix_data->ksp, matrix_data->residual, matrix_data->update);
  PetscInt pits;
  KSPGetIterationNumber(matrix_data->ksp, &pits);
  *its = pits;

  PetscComplex *x_complex = malloc(sizeof(PetscComplex) * matrix_data->local_dof);
  VecGetValues(matrix_data->update, matrix_data->local_dof, matrix_data->local_to_global_complex,
               x_complex);

  for (int i = 0; i < matrix_data->local_dof; i++) {
    int imag_index = matrix_data->imag_local_to_goma_local[i];
    int real_index = matrix_data->real_local_to_goma_local[i];
    if (imag_index > -1) {
      x_[imag_index] = cimag(x_complex[i]);
    }
    x_[real_index] = creal(x_complex[i]);
  }

  free(x_complex);
  MatZeroEntries(matrix_data->mat);
  return 0;
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

  free(matrix_data->imag_local_to_goma_local);
  free(matrix_data->real_local_to_goma_local);
  free(matrix_data->local_to_global_complex);
  free(matrix_data->local_to_global);
  free(matrix_data->is_imag);

  return GOMA_SUCCESS;
}
#endif
#endif
