#include <cstdlib>
#include <vector>

#include "linalg/sparse_matrix.h"
#ifdef GOMA_ENABLE_TPETRA
#include "linalg/sparse_matrix_tpetra.h"
#endif
#ifdef GOMA_ENABLE_EPETRA
#include "linalg/sparse_matrix_epetra.h"
#endif

extern "C" {
#define DISABLE_CPP
#include "dp_comm.h"
#include "dp_types.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_eh.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_unknown_map.h"
#include "rf_masks.h"
#include "rf_node_const.h"
#include "sl_util_structs.h"
#undef DISABLE_CPP
}

extern "C" goma_error GomaSparseMatrix_CreateFromFormat(GomaSparseMatrix *matrix,
                                                        char *matrix_format) {
  if (strcmp(matrix_format, "tpetra") == 0) {
    return GomaSparseMatrix_Create(matrix, GOMA_SPARSE_MATRIX_TYPE_TPETRA);
  } else if (strcmp(matrix_format, "epetra") == 0) {
    return GomaSparseMatrix_Create(matrix, GOMA_SPARSE_MATRIX_TYPE_EPETRA);
  }
  return GOMA_ERROR;
}

extern "C" goma_error GomaSparseMatrix_Create(GomaSparseMatrix *matrix,
                                              enum GomaSparseMatrixType type) {
  *matrix = (GomaSparseMatrix)malloc(sizeof(struct g_GomaSparseMatrix));
  switch (type) {
#ifdef GOMA_ENABLE_TPETRA
  case GOMA_SPARSE_MATRIX_TYPE_TPETRA:
    return GomaSparseMatrix_Tpetra_Create(matrix);
    break;
#endif
#ifdef GOMA_ENABLE_EPETRA
  case GOMA_SPARSE_MATRIX_TYPE_EPETRA:
    return GomaSparseMatrix_Epetra_Create(matrix);
    break;
#endif
  default:
    GOMA_EH(GOMA_ERROR, "Unknown matrix type, GomaSparseMatrix_Create");
    return GOMA_ERROR;
    break;
  }
  return GOMA_SUCCESS;
}

extern "C" goma_error GomaSparseMatrix_SetProblemGraph(
    GomaSparseMatrix matrix,
    int num_internal_dofs,
    int num_boundary_dofs,
    int num_external_dofs,
    int local_nodes,
    NODE_INFO_STRUCT **Nodes,
    int MaxVarPerNode,
    int *Matilda,
    int Inter_Mask[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES],
    Exo_DB *exo,
    Dpi *dpi,
    Comm_Ex *cx,
    int imtrx,
    int Debug_Flag,
    struct GomaLinearSolverData *ams) {
  int j, inode, i1, i2, eb1;
  int iunknown, inter_unknown, inter_node, row_num_unknowns, col_num_unknowns;
  int irow_index = 0;
  int icol_index, rowVarType, colVarType;
  int add_var = 0;
  NODE_INFO_STRUCT *nodeCol;
  NODAL_VARS_STRUCT *nv, *nvCol;
  int inode_varType[MaxVarPerNode], inode_matID[MaxVarPerNode];
  int inter_node_varType[MaxVarPerNode], inter_node_matID[MaxVarPerNode];
  int nnz = 0;

  GomaGlobalOrdinal NumMyRows = num_internal_dofs + num_boundary_dofs;
  GomaGlobalOrdinal NumExternal = num_external_dofs;
  GomaGlobalOrdinal NumMyCols = NumMyRows + NumExternal;
  matrix->n_rows = NumMyRows;
  matrix->n_cols = NumMyCols;

  GomaGlobalOrdinal RowOffset;
  MPI_Scan(&NumMyRows, &RowOffset, 1, MPI_GOMA_ORDINAL, MPI_SUM, MPI_COMM_WORLD);
  RowOffset -= NumMyRows;

  std::vector<GomaGlobalOrdinal> GlobalIDs(NumMyCols);

  for (int i = 0; i < NumMyRows; i++) {
    GlobalIDs[i] = i + RowOffset;
  }
  matrix->global_ids = (GomaGlobalOrdinal *)malloc(sizeof(GomaGlobalOrdinal) * NumMyCols);

#ifdef GOMA_MATRIX_GO_LONG_LONG
  exchange_dof_long_long(cx, dpi, GlobalIDs.data(), imtrx);
#else
  exchange_dof_int(cx, dpi, GlobalIDs.data(), imtrx);
#endif

  for (size_t i = 0; i < GlobalIDs.size(); i++) {
    matrix->global_ids[i] = GlobalIDs[i];
  }

  std::vector<GomaGlobalOrdinal> rows(GlobalIDs.begin(), GlobalIDs.begin() + NumMyRows);
  std::vector<GomaGlobalOrdinal> cols(GlobalIDs.begin(), GlobalIDs.end());
  std::vector<GomaGlobalOrdinal> coo_rows;
  std::vector<GomaGlobalOrdinal> coo_cols;

  int max_nz_per_row = 0;
  int row_nz;

  /*
   * loop over all of the nodes on this processor
   */
  for (inode = 0; inode < local_nodes; inode++) {
    nv = Nodes[inode]->Nodal_Vars_Info[pg->imtrx];
    /*
     * Fill the vector list which points to the unknowns defined at this
     * node...
     */
    row_num_unknowns = fill_variable_vector(inode, inode_varType, inode_matID);
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
    for (iunknown = 0; iunknown < row_num_unknowns; iunknown++) {
      row_nz = 0;
      /*
       * Retrieve the var type of the current unknown
       */
      rowVarType = inode_varType[iunknown];

      /*
       * Loop over the nodes which are determined to have an interaction
       * with the current row node
       */
      for (j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1]; j++) {
        inter_node = exo->node_node_list[j];
        nodeCol = Nodes[inter_node];
        nvCol = nodeCol->Nodal_Vars_Info[pg->imtrx];

        /*
         * fill the vector list which points to the unknowns
         * defined at this interaction node
         */
        col_num_unknowns = fill_variable_vector(inter_node, inter_node_varType, inter_node_matID);
        if (col_num_unknowns != nvCol->Num_Unknowns) {
          GOMA_EH(GOMA_ERROR, "Inconsistency counting unknowns.");
        }

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be an interaction.
         */

        for (inter_unknown = 0; inter_unknown < col_num_unknowns; inter_unknown++) {

          /*
           * HKM ->
           *  The entire process below is designed to find out
           *  what unknown we are processing. This coding is very
           *  convoluted and fraught with pitfalls. It should
           *  be rewritten using the structured approach in MPSalsa
           *  to variable indentification.
           */
          colVarType = inter_node_varType[inter_unknown];

          /*
           * Query the Interaction mask to determine if a jacobian entry
           * should be created
           */
          add_var = Inter_Mask[pg->imtrx][rowVarType][colVarType];

          /* The following code should be activated when solving DG viscoelastic problems
           * with full Jacobian treatment of upwind element stress terms
           */
          if (exo->centroid_list[inode] != -1 && inode != inter_node &&
              exo->centroid_list[inter_node] != -1) {
            eb1 = exo->elem_eb[exo->centroid_list[inode]];

            if (vn_glob[Matilda[eb1]]->dg_J_model == FULL_DG) {
              i1 = pd_glob[Matilda[eb1]]->i[pg->imtrx][rowVarType];
              i2 = pd_glob[Matilda[eb1]]->i[pg->imtrx][colVarType];

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
            icol_index = nodeCol->First_Unknown[pg->imtrx] + inter_unknown;
            coo_rows.push_back(GlobalIDs[irow_index]);
            coo_cols.push_back(GlobalIDs[icol_index]);
            nnz++;
            row_nz++;
          }
        }
      }
      if (row_nz > max_nz_per_row) {
        max_nz_per_row = row_nz;
      }
      irow_index++;
    }
  }

  matrix->create_graph(matrix, NumMyRows, rows.data(), NumMyCols, cols.data(), nnz, max_nz_per_row,
                       coo_rows.data(), coo_cols.data());

  if (matrix->complete_graph != NULL) {
    //  matrix->complete_graph(matrix);
  }

  /*
   * Add ams values that are needed elsewhere
   */
  ams->nnz = nnz;
  ams->nnz_plus = nnz;

  ams->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  ;
  ams->npn_plus = dpi->num_universe_nodes;

  ams->npu = num_internal_dofs + num_boundary_dofs;
  ams->npu_plus = num_internal_dofs + num_boundary_dofs + num_external_dofs;

  int64_t num_unknowns;
  int64_t my_unknowns = ams->npu_plus;

  int64_t num_nzz_global;
  int64_t my_nnz = nnz;

  MPI_Allreduce(&my_unknowns, &num_unknowns, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&my_nnz, &num_nzz_global, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  DPRINTF(stdout, "\n%-30s= %ld\n", "Number of unknowns", num_unknowns);
  DPRINTF(stdout, "\n%-30s= %ld\n", "Number of matrix nonzeroes", num_nzz_global);
  return GOMA_SUCCESS;
}

extern "C" goma_error GomaSparseMatrix_LoadLec(GomaSparseMatrix matrix,
                                               int ielem,
                                               struct Local_Element_Contributions *lec,
                                               double resid_vector[]) {
  int e, v, i, j, pe, pv;
  int dofs;
  int gnn, row_index, ke, kv, nvdof;
  int col_index, ledof;
  int je_new;
  struct Element_Indices *ei_ptr;
  std::vector<GomaGlobalOrdinal> Indices;
  std::vector<double> Values;

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
                          Indices.push_back(matrix->global_ids[col_index]);
                          Values.push_back(lec->J[LEC_J_INDEX(pe, pv, i, j)]);
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
                        Indices.push_back(matrix->global_ids[col_index]);
                        Values.push_back(lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                      }
                    }
                  }
                }
                matrix->sum_into_row_values(matrix, matrix->global_ids[row_index], Indices.size(),
                                            &Values[0], &Indices[0]);
                Indices.clear();
                Values.clear();
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
                        Indices.push_back(matrix->global_ids[col_index]);
                        Values.push_back(lec->J[LEC_J_INDEX(pe, pv, i, j)]);
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
                      Indices.push_back(matrix->global_ids[col_index]);
                      Values.push_back(lec->J[LEC_J_INDEX(pe, pv, i, j)]);
                    }
                  }
                }
              }
              matrix->sum_into_row_values(matrix, matrix->global_ids[row_index], Indices.size(),
                                          &Values[0], &Indices[0]);
              Indices.clear();
              Values.clear();
            }
          }
        }
      }
    }
  }
  return GOMA_SUCCESS;
}

extern "C" goma_error GomaSparseMatrix_Destroy(GomaSparseMatrix *matrix) {
  if (*matrix == NULL) {
    return GOMA_SUCCESS;
  }
  if ((*matrix)->destroy != NULL) {
    (*matrix)->destroy(*matrix);
  }
  free((*matrix)->global_ids);
  free(*matrix);
  *matrix = NULL;
  return GOMA_SUCCESS;
}