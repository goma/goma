#ifdef GOMA_ENABLE_EPETRA
#include "std.h"
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCPDecl.hpp>
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <mpi.h>
#include <unordered_set>

extern "C" {
#include "mm_eh.h"
#include "rf_io.h"
}
#include "linalg/sparse_matrix.h"
#include "linalg/sparse_matrix_epetra.h"

using Teuchos::RCP;

extern "C" goma_error GomaSparseMatrix_Epetra_Create(GomaSparseMatrix *matrix) {
  EpetraSparseMatrix *tmp = new EpetraSparseMatrix();
  (*matrix)->type = GOMA_SPARSE_MATRIX_TYPE_EPETRA;
  (*matrix)->data = reinterpret_cast<void *>(tmp);
  (*matrix)->create_graph = g_epetra_create_graph;
  (*matrix)->complete_graph = g_epetra_complete_graph;
  (*matrix)->insert_row_values = g_epetra_insert_row_values;
  (*matrix)->sum_into_row_values = g_epetra_sum_into_row_values;
  (*matrix)->put_scalar = g_epetra_put_scalar;
  (*matrix)->row_sum_scaling = g_epetra_row_sum_scaling;
  (*matrix)->zero_global_row = g_epetra_zero_row;
  (*matrix)->zero_global_row_set_diag = g_epetra_zero_row_set_diag;
  (*matrix)->destroy = g_epetra_destroy;
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_create_graph(GomaSparseMatrix matrix,
                                            GomaGlobalOrdinal n_rows,
                                            GomaGlobalOrdinal *row_list,
                                            GomaGlobalOrdinal n_cols,
                                            GomaGlobalOrdinal *col_list,
                                            GomaGlobalOrdinal local_nnz,
                                            GomaGlobalOrdinal max_per_row,
                                            GomaGlobalOrdinal *coo_rows,
                                            GomaGlobalOrdinal *coo_cols) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  GomaGlobalOrdinal global_n_rows;
  GomaGlobalOrdinal global_n_cols;
  MPI_Allreduce(&n_rows, &global_n_rows, 1, MPI_GOMA_ORDINAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&n_cols, &global_n_cols, 1, MPI_GOMA_ORDINAL, MPI_SUM, MPI_COMM_WORLD);

  tmp->row_map = Teuchos::rcp(new Epetra_Map(global_n_rows, n_rows, row_list, 0, comm));
  tmp->crs_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *(tmp->row_map), max_per_row));

  GomaGlobalOrdinal current_row = row_list[0];
  std::vector<GomaGlobalOrdinal> indices(max_per_row);
  std::vector<double> values(max_per_row);
  GomaGlobalOrdinal row_count = 0;
  for (GomaGlobalOrdinal i = 0; i < local_nnz; i++) {
    if (coo_rows[i] != current_row) {
      tmp->crs_graph->InsertGlobalIndices(current_row, row_count, indices.data());
      current_row = coo_rows[i];
      row_count = 0;
    }
    indices[row_count] = coo_cols[i];
    row_count++;
  }
  tmp->crs_graph->InsertGlobalIndices(current_row, row_count, indices.data());

  tmp->crs_graph->FillComplete();

  tmp->matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(tmp->crs_graph)));
  tmp->matrix->FillComplete();
  tmp->matrix->PutScalar(0.0);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_complete_graph(GomaSparseMatrix matrix) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  tmp->matrix->FillComplete();
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_insert_row_values(GomaSparseMatrix matrix,
                                                 GomaGlobalOrdinal global_row,
                                                 GomaGlobalOrdinal num_entries,
                                                 double *values,
                                                 GomaGlobalOrdinal *indices) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  GO global_row_t = static_cast<GO>(global_row);
  tmp->matrix->InsertGlobalValues(global_row_t, num_entries, values, indices);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_sum_into_row_values(GomaSparseMatrix matrix,
                                                   GomaGlobalOrdinal global_row,
                                                   GomaGlobalOrdinal num_entries,
                                                   double *values,
                                                   GomaGlobalOrdinal *indices) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  GO global_row_t = static_cast<GO>(global_row);
  tmp->matrix->SumIntoGlobalValues(global_row_t, num_entries, values, indices);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_put_scalar(GomaSparseMatrix matrix, double scalar) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  tmp->matrix->PutScalar(scalar);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_row_sum_scaling(GomaSparseMatrix matrix, double *b, double *scale) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  tmp->matrix->RowMatrixRowMap();
  Epetra_Vector vector_scale(tmp->matrix->RowMatrixRowMap(), false);
  tmp->matrix->InvRowSums(vector_scale);
  tmp->matrix->LeftScale(vector_scale);
  int local;
  for (int i = 0; i < tmp->matrix->NumMyRows(); i++) {
    local = tmp->matrix->RowMatrixRowMap().LID(matrix->global_ids[i]);
    b[i] *= vector_scale[local];
    scale[i] = 1 / vector_scale[local];
  }
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_zero_row(GomaSparseMatrix matrix, GomaGlobalOrdinal global_row) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);

  int size = tmp->matrix->NumGlobalEntries(global_row);
  std::vector<double> values(size);
  std::vector<int> indices(size);
  int NumEntries;
  tmp->matrix->ExtractGlobalRowCopy(global_row, size, NumEntries, values.data(), indices.data());
  for (int i = 0; i < NumEntries; i++) {
    values[i] = 0;
  }
  tmp->matrix->ReplaceGlobalValues(global_row, NumEntries, values.data(), indices.data());

  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_zero_row_set_diag(GomaSparseMatrix matrix,
                                                 GomaGlobalOrdinal global_row) {
  auto *tmp = static_cast<EpetraSparseMatrix *>(matrix->data);
  int size = tmp->matrix->NumGlobalEntries(global_row);
  std::vector<double> values(size);
  std::vector<int> indices(size);
  int NumEntries;
  tmp->matrix->ExtractGlobalRowCopy(global_row, size, NumEntries, values.data(), indices.data());
  for (int i = 0; i < NumEntries; i++) {
    if (indices[i] == global_row) {
      values[i] = 1;
    } else {
      values[i] = 0;
    }
  }
  tmp->matrix->ReplaceGlobalValues(global_row, NumEntries, values.data(), indices.data());
  return GOMA_SUCCESS;
}

extern "C" goma_error g_epetra_destroy(GomaSparseMatrix matrix) {
  delete static_cast<EpetraSparseMatrix *>(matrix->data);
  return GOMA_SUCCESS;
}
#endif // GOMA_ENABLE_EPETRA