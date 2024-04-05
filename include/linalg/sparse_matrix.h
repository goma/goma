#ifndef GOMA_SPARSE_MATRIX
#define GOMA_SPARSE_MATRIX

#ifdef __cplusplus
extern "C" {
#endif
#ifdef GOMA_MATRIX_GO_LONG_LONG
typedef long long GomaGlobalOrdinal;
#define MPI_GOMA_ORDINAL MPI_LONG_LONG
#else
typedef int GomaGlobalOrdinal;
#define MPI_GOMA_ORDINAL MPI_INT
#endif

#ifdef DISABLE_CPP
#define CPP_ALREADY_DISABLED
#endif
#define DISABLE_CPP
#include "dp_types.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_node_const.h"
#include "sl_util_structs.h"
#ifndef CPP_ALREADY_DISABLED
#undef DISABLE_CPP
#endif

enum GomaSparseMatrixType {
  GOMA_SPARSE_MATRIX_TYPE_EPETRA,
  GOMA_SPARSE_MATRIX_TYPE_AZTEC_MRS,
  GOMA_SPARSE_MATRIX_TYPE_TPETRA,
  GOMA_SPARSE_MATRIX_TYPE_PETSC
};

struct g_GomaSparseMatrix {
  enum GomaSparseMatrixType type;
  void *data;
  GomaGlobalOrdinal *global_ids;
  GomaGlobalOrdinal n_rows;
  GomaGlobalOrdinal n_cols;
  GomaGlobalOrdinal nnz;
  // Create matrix with given rows and columns
  // coo_rows and coo_cols are the COO representation of the matrix ordered by row
  // local_nnz is the number of non-zero entries in the local partition
  // max_nz_per_row is the maximum number of non-zero entries in any row
  goma_error (*create_graph)(struct g_GomaSparseMatrix *matrix,
                             GomaGlobalOrdinal n_rows,
                             GomaGlobalOrdinal *row_list,
                             GomaGlobalOrdinal n_cols,
                             GomaGlobalOrdinal *col_list,
                             GomaGlobalOrdinal local_nnz,
                             GomaGlobalOrdinal max_nz_per_row,
                             GomaGlobalOrdinal *coo_rows,
                             GomaGlobalOrdinal *coo_cols);
  // optional, run after create graph to finalize structure
  goma_error (*complete_graph)(struct g_GomaSparseMatrix *matrix);
  // Insert values into matrix row replacing existing values
  goma_error (*insert_row_values)(struct g_GomaSparseMatrix *matrix,
                                  GomaGlobalOrdinal global_row,
                                  GomaGlobalOrdinal num_entries,
                                  double *values,
                                  GomaGlobalOrdinal *indices);
  // Add values into a matrix row, sums existing values
  goma_error (*sum_into_row_values)(struct g_GomaSparseMatrix *matrix,
                                    GomaGlobalOrdinal global_row,
                                    GomaGlobalOrdinal num_entries,
                                    double *values,
                                    GomaGlobalOrdinal *indices);
  // set matrix non-zeros to specified scalar value (commonly re-zero for next assembly)
  goma_error (*put_scalar)(struct g_GomaSparseMatrix *matrix, double scalar);
  // row sum scaling, compute row sum scale, scale matrix and b, and return scaling vector
  goma_error (*row_sum_scaling)(struct g_GomaSparseMatrix *matrix, double *b, double *scale);
  // Zeros a global row
  goma_error (*zero_global_row)(struct g_GomaSparseMatrix *matrix, GomaGlobalOrdinal global_row);
  // Zeros a global row and sets diagonal to 1.0
  goma_error (*zero_global_row_set_diag)(struct g_GomaSparseMatrix *matrix,
                                         GomaGlobalOrdinal global_row);
  // delete the allocated matrix;
  goma_error (*destroy)(struct g_GomaSparseMatrix *matrix);
};
typedef struct g_GomaSparseMatrix *GomaSparseMatrix;

goma_error GomaSparseMatrix_CreateFromFormat(GomaSparseMatrix *matrix, char *matrix_format);
goma_error GomaSparseMatrix_Create(GomaSparseMatrix *matrix, enum GomaSparseMatrixType type);
// populates:
//   - global_ids
//   - n_rows
//   - n_cols
goma_error GomaSparseMatrix_SetProblemGraph(
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
    struct GomaLinearSolverData *ams);

goma_error GomaSparseMatrix_LoadLec(GomaSparseMatrix matrix,
                                    int ielem,
                                    struct Local_Element_Contributions *lec,
                                    double resid_vector[]);

goma_error GomaSparseMatrix_Destroy(GomaSparseMatrix *matrix);

#ifdef __cplusplus
}
#endif

#endif // GOMA_SPARSE_MATRIX