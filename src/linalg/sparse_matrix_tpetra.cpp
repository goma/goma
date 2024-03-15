#include "Tpetra_computeRowAndColumnOneNorms.hpp"
#include "std.h"
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_computeRowAndColumnOneNorms_decl.hpp>
#include <cstdint>
#include <mpi.h>
#include <unordered_set>
#include <iostream>

extern "C" {
#include "rf_io.h"
#include "mm_eh.h"
}
#include "linalg/sparse_matrix.h"
#include "linalg/sparse_matrix_tpetra.h"

using Teuchos::RCP;


extern "C" goma_error GomaSparseMatrix_Tpetra_Create(GomaSparseMatrix *matrix) {
  TpetraSparseMatrix *tmp = new TpetraSparseMatrix();
  (*matrix)->type = GOMA_SPARSE_MATRIX_TYPE_TPETRA;
  (*matrix)->data = reinterpret_cast<void *>(tmp);
  (*matrix)->create_graph = g_tpetra_create_graph;
  (*matrix)->complete_graph = complete_graph;
  (*matrix)->insert_row_values = g_tpetra_insert_row_values;
  (*matrix)->sum_into_row_values = g_tpetra_sum_into_row_values;
  (*matrix)->put_scalar = g_tpetra_put_scalar;
  (*matrix)->row_sum_scaling = g_tpetra_row_sum_scaling;
  (*matrix)->destroy = g_tpetra_destroy;
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_create_graph(GomaSparseMatrix matrix,
                                      GomaGlobalOrdinal n_rows,
                                      GomaGlobalOrdinal *row_list,
                                      GomaGlobalOrdinal n_cols,
                                      GomaGlobalOrdinal *col_list,
                                      GomaGlobalOrdinal local_nnz,
                                      GomaGlobalOrdinal max_per_row,
                                      GomaGlobalOrdinal *coo_rows,
                                      GomaGlobalOrdinal *coo_cols) {
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  RCP<const Teuchos::MpiComm<int>> comm(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  GomaGlobalOrdinal global_n_rows;
  GomaGlobalOrdinal global_n_cols;
  MPI_Allreduce(&n_rows, &global_n_rows, 1, MPI_GOMA_ORDINAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&n_cols, &global_n_cols, 1, MPI_GOMA_ORDINAL, MPI_SUM, MPI_COMM_WORLD);

  Teuchos::ArrayView<GomaGlobalOrdinal> row_list_view(row_list, n_rows);
  Teuchos::ArrayView<GomaGlobalOrdinal> col_list_view(col_list, n_cols);

  tmp->row_map = Teuchos::rcp(new Tpetra::Map<LO, GO>(global_n_rows, row_list_view, 0, comm));
  tmp->col_map = Teuchos::rcp(new Tpetra::Map<LO, GO>(global_n_rows, row_list_view, 0, comm));
  tmp->crs_graph = Teuchos::rcp(new Tpetra::FECrsGraph<LO, GO>(tmp->row_map, tmp->col_map, max_per_row));


  GomaGlobalOrdinal current_row = row_list[0];
  std::vector<GomaGlobalOrdinal> indices(max_per_row);
  std::vector<double> values(max_per_row);
  GomaGlobalOrdinal row_count = 0;
  for (GomaGlobalOrdinal i = 0; i < local_nnz; i++) {
    if (coo_rows[i] != current_row) {
      tmp->crs_graph->insertGlobalIndices(current_row, Teuchos::ArrayView<GomaGlobalOrdinal>(indices.data(), row_count));
      current_row = coo_rows[i];
      row_count = 0;
    }
    indices[row_count] = coo_cols[i];
    values[row_count] = coo_cols[i];
    row_count++;
  }
   tmp->crs_graph->insertGlobalIndices(current_row, Teuchos::ArrayView<GomaGlobalOrdinal>(indices.data(), row_count));

  if (Debug_Flag > 2) {
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    tmp->crs_graph->describe(*fos, Teuchos::VERB_EXTREME);
  }
  tmp->crs_graph->fillComplete();

  tmp->matrix = Teuchos::rcp(new Tpetra::FECrsMatrix<double, LO, GO>((tmp->crs_graph)));
  tmp->matrix->beginAssembly();
  return GOMA_SUCCESS;
}

extern "C" goma_error complete_graph(GomaSparseMatrix matrix) {
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  tmp->matrix->fillComplete();
  tmp->matrix->resumeFill();
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_insert_row_values(GomaSparseMatrix matrix,
                                           GomaGlobalOrdinal global_row,
                                           GomaGlobalOrdinal num_entries,
                                           double *values,
                                           GomaGlobalOrdinal *indices) {
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  GO global_row_t = static_cast<GO>(global_row);
  tmp->matrix->insertGlobalValues(global_row_t, Teuchos::ArrayView<GomaGlobalOrdinal>(indices, num_entries),
                                  Teuchos::ArrayView<double>(values, num_entries));
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_sum_into_row_values(GomaSparseMatrix matrix,
                                             GomaGlobalOrdinal global_row,
                                             GomaGlobalOrdinal num_entries,
                                             double *values,
                                             GomaGlobalOrdinal *indices) {
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  GO global_row_t = static_cast<GO>(global_row);
  tmp->matrix->sumIntoGlobalValues(global_row_t,
                                   num_entries,
                                   values,
                                   indices);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_put_scalar(GomaSparseMatrix matrix, double scalar) {
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  tmp->matrix->setAllToScalar(scalar);
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_row_sum_scaling(GomaSparseMatrix matrix, double *b, double *scale) {
  // noop
  return GOMA_SUCCESS;
  auto *tmp = static_cast<TpetraSparseMatrix *>(matrix->data);
  RCP<Tpetra::Vector<double, LO, GO>> b_vec =
      Teuchos::rcp(new Tpetra::Vector<double, LO, GO>(tmp->row_map));
  auto n_rows = tmp->matrix->getLocalNumRows();
  using crs_t = Tpetra::CrsMatrix<double, LO, GO>;
  typename crs_t::local_inds_host_view_type indices;
  typename crs_t::values_host_view_type values;
  for (size_t i = 0; i < n_rows; i++) {
    tmp->matrix->getLocalRowView(i, indices, values);
    double row_sum = 0;
    for (size_t j = 0; j < values.size(); j++) {
      row_sum += std::abs(values[j]);
    }
    b_vec->replaceLocalValue(i, 1.0 / row_sum);
  }
  tmp->matrix->leftScale(*b_vec);
  auto data = b_vec->getData();
  for (size_t i = 0; i < n_rows; i++) {
    scale[i] = data[i];
  }
  return GOMA_SUCCESS;
}

extern "C" goma_error g_tpetra_destroy(GomaSparseMatrix matrix) {
  delete static_cast<TpetraSparseMatrix *>(matrix->data);
  return GOMA_SUCCESS;
}