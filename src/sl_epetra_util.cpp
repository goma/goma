#ifndef __cplusplus
#define __cplusplus
#endif

#if defined(PARALLEL) && !defined(EPETRA_MPI)
#define EPETRA_MPI
#endif

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_ConfigDefs.h"
#include "dpi.h"
#include "exo_struct.h"

#ifdef EPETRA_MPI
#else
#include "Epetra_SerialComm.h"
#endif

extern "C" {
#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io.h"
#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "rf_node_const.h"
#include "mm_eh.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "dp_comm.h"
#include "dp_types.h"
#include "mm_unknown_map.h"
#include "rf_solve.h"
#include "sl_util_structs.h"
#include "mm_fill_util.h"
}

#include "sl_epetra_util.h"

#include "sl_epetra_interface.h"

extern "C" {

/**
 * Create the goma problem graph in the epetra matrix ams->RowMatrix
 * @param ams ams structure containing appropriate RowMatrix
 * @param exo exodus file for this processor
 */
void EpetraCreateGomaProblemGraph(struct Aztec_Linear_Solver_System *ams, Exo_DB *exo, Dpi *dpi) {
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
  int total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;
  std::vector<int> Indices;
  std::vector<double> Values;

  int NumMyRows = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
  int NumExternal = num_external_dofs[pg->imtrx];
  int NumMyCols = NumMyRows + NumExternal;

  /*
   * This is kind of hacky the way it is being done,
   * modified from the amesos interface for gomamsr to epetra
   *
   * Creates an array to be communicated of ints cast to double to use the exchange_dof
   * communicator, then casts doubles back to int and sets those as global ids
   * for the epetra array
   *
   * TODO: replace with non-double conversion routine
   */

  // get the row map from the row matrix
  Epetra_Map RowMap = ams->RowMatrix->RowMatrixRowMap();

  // get the global elements for this processor
  int *MyGlobalElements = RowMap.MyGlobalElements ();

  double *dblColGIDs = new double[NumMyCols];
  ams->GlobalIDs = (int *) malloc(sizeof(int)*NumMyCols);

  // copy global id's and convert to double for boundary exchange
  for( int i=0; i<NumMyRows; i++) dblColGIDs[i] = (double) MyGlobalElements[i];

  exchange_dof(cx[pg->imtrx], dpi, dblColGIDs, pg->imtrx);

  // convert back to int with known global id's from all processors
  for (int j = 0; j < NumMyCols; j++) {
    ams->GlobalIDs[j] = (int) dblColGIDs[j];
  }

  /*
   * loop over all of the nodes on this processor
   */
  for (inode = 0; inode < total_nodes; inode++) {
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
      EH(GOMA_ERROR, "Inconsistency counting unknowns.");
    }

    /*
     * Loop over the unknowns defined at this row node
     */
    for (iunknown = 0; iunknown < row_num_unknowns; iunknown++) {
      /*
       * Retrieve the var type of the current unknown
       */
      rowVarType = inode_varType[iunknown];

      Indices.clear();
      Values.clear();

      /*
       * Loop over the nodes which are determined to have an interaction
       * with the current row node
       */
      for (j = exo->node_node_pntr[inode]; j < exo->node_node_pntr[inode + 1];
          j++) {
        inter_node = exo->node_node_list[j];
        nodeCol = Nodes[inter_node];
        nvCol = nodeCol->Nodal_Vars_Info[pg->imtrx];

        /*
         * fill the vector list which points to the unknowns
         * defined at this interaction node
         */
        col_num_unknowns = fill_variable_vector(inter_node, inter_node_varType,
            inter_node_matID);
        if (col_num_unknowns != nvCol->Num_Unknowns) {
          EH(GOMA_ERROR, "Inconsistency counting unknowns.");
        }

        /*
         * Loop over the unknowns associated with the
         * interacting node and see if there should be an interaction.
         */

        for (inter_unknown = 0; inter_unknown < col_num_unknowns;
            inter_unknown++) {

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
          if (exo->centroid_list[inode] != -1 && inode != inter_node
              && exo->centroid_list[inter_node] != -1) {
            eb1 = exo->elem_eb[exo->centroid_list[inode]];

            if (vn_glob[Matilda[eb1]]->dg_J_model == FULL_DG) {
              i1 = pd_glob[Matilda[eb1]]->i[pg->imtrx][rowVarType];
              i2 = pd_glob[Matilda[eb1]]->i[pg->imtrx][colVarType];

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
            icol_index = nodeCol->First_Unknown[pg->imtrx] + inter_unknown;
            Indices.push_back(ams->GlobalIDs[icol_index]);
            Values.push_back(0);
          }
        }
      }
      EpetraInsertGlobalRowMatrix(ams->RowMatrix, ams->GlobalIDs[irow_index], Indices.size(), &Values[0],
          &Indices[0]);
      nnz += Indices.size();
      irow_index++;
    }
  }

  EpetraFillCompleteRowMatrix(ams->RowMatrix);
  EpetraPutScalarRowMatrix(ams->RowMatrix, 0);

  /*
   * Add ams values that are needed elsewhere
   */
  ams->nnz = nnz;
  ams->nnz_plus = nnz;

  ams->npn = dpi->num_internal_nodes + dpi->num_boundary_nodes;;
  ams->npn_plus = dpi->num_universe_nodes;

  ams->npu = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
  ams->npu_plus = num_universe_dofs[pg->imtrx];

  delete[] dblColGIDs;

  DPRINTF(stdout, "\n%-30s= %d\n", "Number of unknowns", num_universe_dofs[pg->imtrx]);
  DPRINTF(stdout, "\n%-30s= %d\n", "Number of matrix nonzeroes", nnz);
}

/**
 * Load local element contributions into the local matrix on this processor
 *
 * Modified from MSR version in mm_fill load_lec
 *
 * @param exo ptr to EXODUS II finite element mesh db
 * @param ielem Element number we are working on
 * @param ams Matrix contianer
 * @param x Solution vector
 * @param resid_vector residual vector
 */
void EpetraLoadLec(int ielem, struct Aztec_Linear_Solver_System *ams,
                   double resid_vector[])
{
  int e, v, i, j, pe, pv;
  int dofs;
  int gnn, row_index, ke, kv, nvdof;
  int col_index, ledof;
  int je_new;
  struct Element_Indices *ei_ptr;
  std::vector<int> Indices;
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
              row_index = Index_Solution(gnn, e, ke, nvdof,
					 ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
              resid_vector[row_index] += lec->R[MAX_PROB_VAR + ke][i];

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
                          EH(col_index, "Bad var index.");
                          Indices.push_back(ams->GlobalIDs[col_index]);
                          Values.push_back(lec->J[pe][pv][i][j]);
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
                          EH(GOMA_ERROR, "LEC Indexing error");
                        }
                        EH(col_index, "Bad var index.");
                        Indices.push_back(ams->GlobalIDs[col_index]);
                        Values.push_back(lec->J[pe][pv][i][j]);
                      }
                    }
                  }
                }
                EpetraSumIntoGlobalRowMatrix(ams->RowMatrix, ams->GlobalIDs[row_index],
                    Indices.size(), &Values[0], &Indices[0]);
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
            resid_vector[row_index] += lec->R[pe][i];

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
                          EH(GOMA_ERROR, "ei_ptr == 0\n");
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
                        EH(col_index, "Bad var index.");
                        Indices.push_back(ams->GlobalIDs[col_index]);
                        Values.push_back(lec->J[pe][pv][i][j]);
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
                          EH(GOMA_ERROR, "ei slave pointer is null");
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
                        EH(GOMA_ERROR, "LEC Indexing error");
                      }
                      EH(col_index, "Bad var index.");
                      Indices.push_back(ams->GlobalIDs[col_index]);
                      Values.push_back(lec->J[pe][pv][i][j]);
                    }
                  }
                }
              }
              EpetraSumIntoGlobalRowMatrix(ams->RowMatrix, ams->GlobalIDs[row_index], Indices.size(),
                  &Values[0], &Indices[0]);
              Indices.clear();
              Values.clear();
            }
          }
        }
      }
    }
  }
}

/**
 * Inverse row sum scale, used by goma for scaling matrix and residual
 * @param ams Aztec structure containing RowMatrix
 * @param b b from Ax = b
 * @param scale array for scale values to be placed for further usage (b[i] will equal old b[i] / scale[i])
 */
void EpetraRowSumScale(struct Aztec_Linear_Solver_System *ams, double *b, double *scale)
{
  Epetra_Vector vector_scale(ams->RowMatrix->RowMatrixRowMap(), false);
  ams->RowMatrix->InvRowSums(vector_scale);
  ams->RowMatrix->LeftScale(vector_scale);
  int local;
  for (int i = 0; i < ams->RowMatrix->NumMyRows(); i++) {
    local = ams->RowMatrix->RowMatrixRowMap().LID(ams->GlobalIDs[i]);
    b[i] *= vector_scale[local];
    scale[i] = 1 / vector_scale[local];
  }
}

/**
 * Set the GlobalRow to zero and set the diagonal column to 1 in that row
 *
 * @param ams Aztec_Linear_Solver_System matrix struct containing epetra matrix
 * @param GlobalRow global row to set to diagonal only
 */
void EpetraSetDiagonalOnly(struct Aztec_Linear_Solver_System *ams, int GlobalRow)
{
  Epetra_CrsMatrix* CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(ams->RowMatrix);
  int size = CrsMatrix->NumGlobalEntries(GlobalRow);
  double *values = new double[size];
  int *indices = new int[size];
  int NumEntries;
  CrsMatrix->ExtractGlobalRowCopy(GlobalRow, size, NumEntries, values, indices);
  for (int i = 0; i < NumEntries; i++) {
    if (indices[i] == GlobalRow) {
      values[i] = 1;
    } else {
      values[i] = 0;
    }
  }
  CrsMatrix->ReplaceGlobalValues(GlobalRow, NumEntries, values, indices);
  delete [] indices;
  delete [] values;
}

}
/* End extern "C" */
