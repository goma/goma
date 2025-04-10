#include "rf_mp.h"
#ifdef GOMA_ENABLE_MUMPS

#include <stdlib.h>

#include "dmumps_c.h"
#include "mm_eh.h"
#include "sl_mumps.h"
#include "sl_util_structs.h"
#include "std.h"

#define JOB_END          -2
#define JOB_INIT         -1
#define JOB_NULL         0
#define JOB_FACTSYMBOLIC 1
#define JOB_FACTNUMERIC  2
#define JOB_SOLVE        3
#define USE_COMM_WORLD   -987654

// C to Fortran macros for MUMPS
#define ICNTL(I)  icntl[(I) - 1]
#define CNTL(I)   cntl[(I) - 1]
#define INFOG(I)  infog[(I) - 1]
#define INFO(I)   info[(I) - 1]
#define RINFOG(I) rinfog[(I) - 1]
#define RINFO(I)  rinfo[(I) - 1]

struct MUMPS_data {
  DMUMPS_STRUC_C mumps;
  int N;
  int nnz;
  int *irn;
  int *jcn;
  double *val;
  double *rhs;
};

void msr_to_triplet(struct GomaLinearSolverData *data, dbl *rhs) {
  struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;

  int N = data->N + data->N_update;
  int local_nnz = data->nnz;

  int RowOffset;
  MPI_Scan(&N, &RowOffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  RowOffset -= N;
  int *colgids = (int *)malloc((N + data->N_external) * sizeof(int));
  for (int i = 0; i < N; i++) {
    colgids[i] = RowOffset + i + 1;
  }

  mumps_data->irn = (int *)malloc(local_nnz * sizeof(int));
  mumps_data->jcn = (int *)malloc(local_nnz * sizeof(int));
  mumps_data->val = (double *)malloc(local_nnz * sizeof(double));
  int nnz_check = 0;
  int offset = 0;
  for (int i = 0; i < N; i++) {
    int row_nnz = data->bindx[i + 1] - data->bindx[i] + 1;

    // diagonal
    mumps_data->irn[offset] = colgids[i];
    mumps_data->jcn[offset] = colgids[i];
    mumps_data->val[offset] = data->val[i];
    offset++;

    // off diagonal
    for (int j = data->bindx[i]; j < data->bindx[i + 1]; j++) {
      int col = data->bindx[j];
      mumps_data->irn[offset] = colgids[i];
      mumps_data->jcn[offset] = colgids[col];
      mumps_data->val[offset] = data->val[col];
      offset++;
    }
    nnz_check += row_nnz;
  }

  GOMA_ASSERT_ALWAYS(nnz_check == local_nnz);
  mumps_data->nnz = local_nnz;
  mumps_data->N = N;
}

void mumps_solve(struct GomaLinearSolverData *data, dbl *x, dbl *rhs) {
  if (data->SolverData == NULL) {

    struct MUMPS_data *mumps_data = malloc(sizeof(struct MUMPS_data));
    DMUMPS_STRUC_C *mumps = &mumps_data->mumps;
    data->SolverData = (void *)mumps_data;

    mumps->comm_fortran = USE_COMM_WORLD;
    mumps->job = JOB_INIT;
    mumps->par = 1;
    mumps->sym = 0;

    if (Num_Proc == 1) {
      mumps->ICNTL(18) = 0;
      mumps->ICNTL(7) = 7;
    } else {
      mumps->ICNTL(18) = 3;
      mumps->ICNTL(21) = 1;
    }
    // mumps->CNTL(1) = -1;
    // mumps->ICNTL(14) = 20;
    // mumps->ICNTL(23) = 0;

    mumps->job = JOB_INIT;
    dmumps_c(mumps);
  }
  struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;
  DMUMPS_STRUC_C *mumps = &mumps_data->mumps;
  msr_to_triplet(data, rhs);
  if (Num_Proc == 1) {
    mumps->n = mumps_data->N;
    mumps->nnz = mumps_data->nnz;
    mumps->irn = mumps_data->irn;
    mumps->jcn = mumps_data->jcn;
    mumps->a = mumps_data->val;
    mumps->rhs = rhs;
  } else {
    mumps->nnz_loc = mumps_data->nnz;
    mumps->irn_loc = mumps_data->irn;
    mumps->jcn_loc = mumps_data->jcn;
    mumps->a_loc = mumps_data->val;
  }

  mumps->job = JOB_FACTSYMBOLIC;
  dmumps_c(mumps);

  mumps->job = JOB_FACTNUMERIC;
  dmumps_c(mumps);

  mumps->job = JOB_SOLVE;
  dmumps_c(mumps);
}

#endif // GOMA_ENABLE_MUMPS