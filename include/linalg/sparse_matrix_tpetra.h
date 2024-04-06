#ifndef GOMA_SPARSE_MATRIX_TPETRA
#define GOMA_SPARSE_MATRIX_TPETRA
#ifdef GOMA_ENABLE_TPETRA
#ifdef __cplusplus
#include "Teuchos_RCP.hpp"
#include <Tpetra_FECrsMatrix.hpp>

#include "linalg/sparse_matrix.h"

using LO = int;
using GO = GomaGlobalOrdinal;

struct TpetraSparseMatrix {
  Teuchos::RCP<Tpetra::FECrsMatrix<double, LO, GO>> matrix;
  Teuchos::RCP<Tpetra::Map<LO, GO>> row_map;
  Teuchos::RCP<Tpetra::Map<LO, GO>> col_map;
  Teuchos::RCP<Tpetra::FECrsGraph<LO, GO>> crs_graph;
  TpetraSparseMatrix() = default;
};

extern "C" {
#endif

goma_error GomaSparseMatrix_Tpetra_Create(GomaSparseMatrix *matrix);

goma_error g_tpetra_create_graph(GomaSparseMatrix matrix,
                                 GomaGlobalOrdinal n_rows,
                                 GomaGlobalOrdinal *row_list,
                                 GomaGlobalOrdinal n_cols,
                                 GomaGlobalOrdinal *col_list,
                                 GomaGlobalOrdinal local_nnz,
                                 GomaGlobalOrdinal max_per_row,
                                 GomaGlobalOrdinal *coo_rows,
                                 GomaGlobalOrdinal *coo_cols);

goma_error g_tpetra_complete_graph(GomaSparseMatrix matrix);

goma_error g_tpetra_insert_row_values(GomaSparseMatrix matrix,
                                      GomaGlobalOrdinal global_row,
                                      GomaGlobalOrdinal num_entries,
                                      double *values,
                                      GomaGlobalOrdinal *indices);

goma_error g_tpetra_sum_into_row_values(GomaSparseMatrix matrix,
                                        GomaGlobalOrdinal global_row,
                                        GomaGlobalOrdinal num_entries,
                                        double *values,
                                        GomaGlobalOrdinal *indices);

goma_error g_tpetra_put_scalar(GomaSparseMatrix matrix, double scalar);

goma_error g_tpetra_row_sum_scaling(GomaSparseMatrix matrix, double *b, double *scale);

goma_error g_tpetra_zero_row(GomaSparseMatrix matrix, GomaGlobalOrdinal global_row);

goma_error g_tpetra_zero_row_set_diag(GomaSparseMatrix matrix, GomaGlobalOrdinal global_row);

goma_error g_tpetra_destroy(GomaSparseMatrix matrix);

#ifdef __cplusplus
}
#endif
#endif
#endif // GOMA_SPARSE_MATRIX_TPETRA