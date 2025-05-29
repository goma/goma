#ifndef GOMA_SPARSE_MATRIX_EPETRA
#define GOMA_SPARSE_MATRIX_EPETRA
#ifdef GOMA_ENABLE_EPETRA
#ifdef __cplusplus
#include "Teuchos_RCP.hpp"
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include "linalg/sparse_matrix.h"

using LO = int;
using GO = GomaGlobalOrdinal;

struct EpetraSparseMatrix {
  Teuchos::RCP<Epetra_CrsMatrix> matrix;
  Teuchos::RCP<Epetra_Map> row_map;
  Teuchos::RCP<Epetra_CrsGraph> crs_graph;
  EpetraSparseMatrix() = default;
};

extern "C" {
#endif

goma_error GomaSparseMatrix_Epetra_Create(GomaSparseMatrix *matrix);

goma_error g_epetra_create_graph(GomaSparseMatrix matrix,
                                 GomaGlobalOrdinal n_rows,
                                 GomaGlobalOrdinal *row_list,
                                 GomaGlobalOrdinal n_cols,
                                 GomaGlobalOrdinal *col_list,
                                 GomaGlobalOrdinal local_nnz,
                                 GomaGlobalOrdinal max_per_row,
                                 GomaGlobalOrdinal *coo_rows,
                                 GomaGlobalOrdinal *coo_cols);

goma_error g_epetra_complete_graph(GomaSparseMatrix matrix);

goma_error g_epetra_insert_row_values(GomaSparseMatrix matrix,
                                      GomaGlobalOrdinal global_row,
                                      GomaGlobalOrdinal num_entries,
                                      double *values,
                                      GomaGlobalOrdinal *indices);

goma_error g_epetra_sum_into_row_values(GomaSparseMatrix matrix,
                                        GomaGlobalOrdinal global_row,
                                        GomaGlobalOrdinal num_entries,
                                        double *values,
                                        GomaGlobalOrdinal *indices);

goma_error g_epetra_put_scalar(GomaSparseMatrix matrix, double scalar);

goma_error g_epetra_row_sum_scaling(GomaSparseMatrix matrix, double *b, double *scale);

goma_error g_epetra_zero_row(GomaSparseMatrix matrix, GomaGlobalOrdinal global_row);

goma_error g_epetra_zero_row_set_diag(GomaSparseMatrix matrix, GomaGlobalOrdinal global_row);

goma_error g_epetra_destroy(GomaSparseMatrix matrix);

#ifdef __cplusplus
}
#endif
#endif
#endif // GOMA_SPARSE_MATRIX_EPETRA