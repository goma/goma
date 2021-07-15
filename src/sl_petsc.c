#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include "dpi.h"
#include "mm_as.h"
#include "dp_comm.h"
#include "el_geom.h"
#include "mm_as_const.h"
#include "rf_fem.h"
#include "rf_io.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "rf_masks.h"
#include "rf_vars_const.h"
#include "sl_petsc.h"
#include "mm_fill_util.h"
#include "mm_unknown_map.h"
#include "mm_mp.h"

typedef struct PetscMatrixData {
  PetscOptions options;
  Mat mat;
  Vec residual;
  Vec update;
  KSP ksp;
  PetscInt *local_to_global;
  PetscBool matrix_setup;
} PetscMatrixData;

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
      for (int j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1];
          j++) {
        int inter_node = exo->node_node_list[j];
        NODE_INFO_STRUCT *nodeCol = Nodes[inter_node];
        NODAL_VARS_STRUCT *nvCol = nodeCol->Nodal_Vars_Info[pg->imtrx];

        /*
         * fill the vector list which points to the unknowns
         * defined at this interaction node
         */
        int inter_node_varType[MaxVarPerNode], inter_node_matID[MaxVarPerNode];
        int col_num_unknowns = fill_variable_vector(inter_node, inter_node_varType,
            inter_node_matID);
        if (col_num_unknowns != nvCol->Num_Unknowns) {
          GOMA_EH(GOMA_ERROR, "Inconsistency counting unknowns.");
        }

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be an interaction.
         */

        for (int inter_unknown = 0; inter_unknown < col_num_unknowns;
            inter_unknown++) {

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
          if (exo->centroid_list[inode] != -1 && inode != inter_node
              && exo->centroid_list[inter_node] != -1) {
            int eb1 = exo->elem_eb[exo->centroid_list[inode]];

            if (vn_glob[Matilda[eb1]]->dg_J_model == FULL_DG) {
              int i1 = pd_glob[Matilda[eb1]]->i[pg->imtrx][rowVarType];
              int i2 = pd_glob[Matilda[eb1]]->i[pg->imtrx][colVarType];

              if ((rowVarType == colVarType)
                  && (i1 == I_P0 || i1 == I_P1 || i1 == I_PQ1 || i1 == I_PQ2)
                  && (i2 == I_P0 || i2 == I_P1 || i2 == I_PQ1 || i2 == I_PQ2)
                  && (rowVarType != PRESSURE)
                  && (rowVarType > VELOCITY_GRADIENT33
                      || rowVarType < VELOCITY_GRADIENT11)) {
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

  ams->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;;
  ams->npn_plus = dpi->num_universe_nodes;

  ams->npu = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
  ams->npu_plus = num_universe_dofs[pg->imtrx];

  DPRINTF(stdout, "\n%-30s= %d\n", "Number of unknowns", num_universe_dofs[pg->imtrx]);
  DPRINTF(stdout, "\n%-30s= %ld\n", "Number of matrix nonzeroes", nnz);

  free(d_nnz);
  free(o_nnz);

  return GOMA_SUCCESS;
}

goma_error goma_setup_petsc_matrix(struct GomaLinearSolverData *ams,
                                   Exo_DB *exo,
                                   Dpi *dpi,
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

  if (GomaPetscOptions != NULL) {
    err = PetscOptionsInsertString(NULL, GomaPetscOptions);
    CHKERRQ(err);
  }

  PetscInt stokes_matrix;
  PetscBool stokes_matrix_set;
  PetscOptionsGetInt(NULL,NULL,"-stokes_matrix",&stokes_matrix, &stokes_matrix_set);

  if (imtrx == 0) {
    CHKERRQ(PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD));
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
    snprintf(ksp_prefix, 24, "sys%d_", imtrx);
    err = MatSetOptionsPrefix(matrix_data->mat, ksp_prefix);
    CHKERRQ(err);
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
  MPI_Allreduce(&local_dof, &global_n_dof, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  err = VecSetSizes(matrix_data->update, internal_dof + boundary_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetSizes(matrix_data->residual, internal_dof + boundary_dof, global_n_dof);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->residual);
  CHKERRQ(err);
  err = VecSetFromOptions(matrix_data->update);
  CHKERRQ(err);
  err = MatSetSizes(matrix_data->mat, internal_dof + boundary_dof,
              internal_dof + boundary_dof, global_n_dof, global_n_dof);
  CHKERRQ(err);
  err = MatSetFromOptions(matrix_data->mat);
  CHKERRQ(err);
  initialize_petsc_matrix(ams, exo, dpi, internal_dof, boundary_dof, external_dof, imtrx);
  err = KSPSetOperators(matrix_data->ksp, matrix_data->mat, matrix_data->mat);
  CHKERRQ(err);
  err = KSPSetFromOptions(matrix_data->ksp);
  CHKERRQ(err);

  if (stokes_matrix_set && stokes_matrix == imtrx) {
    GOMA_WH(GOMA_ERROR, "stokes matrix set to %ld, assuming equal order interpolation and dim %d", stokes_matrix, pd_glob[0]->Num_Dim);
    if (pd_glob[0]->Num_Dim == 3) {
      PC             pc;
      const PetscInt ufields[] = {0,1,2},pfields[] = {3};
      KSPGetPC(matrix_data->ksp, &pc);
      PCFieldSplitSetBlockSize(pc,4);
      PCFieldSplitSetFields(pc,"u",3,ufields,ufields);
      PCFieldSplitSetFields(pc,"p",1,pfields,pfields);
    } else if (pd_glob[0]->Num_Dim == 2) {
      PC             pc;
      const PetscInt ufields[] = {0,1},pfields[] = {3};
      KSPGetPC(matrix_data->ksp, &pc);
      PCFieldSplitSetBlockSize(pc,3);
      PCFieldSplitSetFields(pc,"u",2,ufields,ufields);
      PCFieldSplitSetFields(pc,"p",1,pfields,pfields);
    }
  }

  matrix_data->matrix_setup = PETSC_FALSE;
  return GOMA_SUCCESS;
}

void petsc_load_lec(int ielem, struct GomaLinearSolverData *ams,
                   double resid_vector[])
{
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
              row_index = Index_Solution(gnn, e, ke, nvdof,
					 ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
              resid_vector[row_index] += lec->R[LEC_R_INDEX(MAX_PROB_VAR + ke,i)];

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
                          col_index = Index_Solution(ei_ptr->gnn_list[v][j], v,
                              kv, ei_ptr->Baby_Dolphin[v][j],
						     ei_ptr->matID_ledof[ledof],pg->imtrx);
                          GOMA_EH(col_index, "Bad var index.");

                          PetscInt global_row = matrix_data->local_to_global[row_index];
                          PetscInt global_col = matrix_data->local_to_global[col_index];
                          MatSetValue(matrix_data->mat, global_row, global_col,lec->J[LEC_J_INDEX(pe,pv,i,j)], ADD_VALUES);
                        }
                      }
                    } else {
                      pv = upd->vp[pg->imtrx][v];
                      kv = 0;
                      for (j = 0; j < ei_ptr->dof[v]; j++) {
                        ledof = ei_ptr->lvdof_to_ledof[v][j];
                        je_new = ei_ptr->ieqn_ledof[ledof];
                        col_index = Index_Solution(ei_ptr->gnn_list[v][j], v,
                            kv, ei_ptr->Baby_Dolphin[v][j],
						   ei_ptr->matID_ledof[ledof],pg->imtrx);
                        if (col_index != je_new) {
                          fprintf(stderr,
                              "Oh fiddlesticks: je = %d, je_new = %d\n",
                              col_index, je_new);
                          GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                        }
                        GOMA_EH(col_index, "Bad var index.");
                        PetscInt global_row = matrix_data->local_to_global[row_index];
                        PetscInt global_col = matrix_data->local_to_global[col_index];
                        MatSetValue(matrix_data->mat, global_row, global_col,lec->J[LEC_J_INDEX(pe,pv,i,j)], ADD_VALUES);
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
            resid_vector[row_index] += lec->R[LEC_R_INDEX(pe,i)];

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
                        col_index = Index_Solution(ei_ptr->gnn_list[v][j], v,
                            kv, ei_ptr->Baby_Dolphin[v][j],
						   ei_ptr->matID_ledof[ledof],pg->imtrx);
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
                        MatSetValue(matrix_data->mat, global_row, global_col,lec->J[LEC_J_INDEX(pe,pv,i,j)], ADD_VALUES);
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
                      col_index = Index_Solution(ei_ptr->gnn_list[v][j], v, kv,
                          ei_ptr->Baby_Dolphin[v][j],
						 ei_ptr->matID_ledof[ledof],pg->imtrx);
                      if (col_index != je_new) {
                        fprintf(stderr,
                            "Oh fiddlesticks: je = %d, je_new = %d\n",
                            col_index, je_new);
                        GOMA_EH(GOMA_ERROR, "LEC Indexing error");
                      }
                      GOMA_EH(col_index, "Bad var index.");
                      PetscInt global_row = matrix_data->local_to_global[row_index];
                      PetscInt global_col = matrix_data->local_to_global[col_index];
                      MatSetValue(matrix_data->mat, global_row, global_col,lec->J[LEC_J_INDEX(pe,pv,i,j)], ADD_VALUES);
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

goma_error petsc_scale_matrix(struct GomaLinearSolverData *ams,
		  double *b_,
                  double *scale)  {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  Vec row_sums;

  MatAssemblyBegin(matrix_data->mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->mat,MAT_FINAL_ASSEMBLY);

  CHKERRQ(MatCreateVecs(matrix_data->mat, &row_sums, NULL)); 
  VecZeroEntries(row_sums);
  for (int i = 0; i < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]; i++) {
    PetscInt ncols;
    const PetscScalar *vals;
    MatGetRow(matrix_data->mat, matrix_data->local_to_global[i], &ncols, NULL, &vals);
    dbl sum = 0.0;
    for (PetscInt c = 0; c < ncols; c++)  {
      sum += PetscAbsScalar(vals[c]);
    }
    MatRestoreRow(matrix_data->mat, matrix_data->local_to_global[i], &ncols, NULL, &vals);
    scale[i] = 1.0 / sum;
  }
  VecSetValues(row_sums, num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx], matrix_data->local_to_global, scale, INSERT_VALUES);
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

goma_error goma_petsc_free_matrix(struct GomaLinearSolverData *ams) {
  PetscMatrixData *matrix_data = (PetscMatrixData *) ams->PetscMatrixData;
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

// vim: expandtab sw=2 ts=8
void
petsc_solve(struct GomaLinearSolverData *ams,
		  double *x_, 
		  double *b_,
                  int *its) {
  PetscMatrixData *matrix_data = (PetscMatrixData *)ams->PetscMatrixData;
  matrix_data->matrix_setup = PETSC_TRUE;
  MatAssemblyBegin(matrix_data->mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix_data->mat,MAT_FINAL_ASSEMBLY);
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
  VecGetValues(matrix_data->update, num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx], matrix_data->local_to_global, x_);
  MatZeroEntries(matrix_data->mat);
}
