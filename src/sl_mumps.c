#include "ac_conti.h"
#include "dp_comm.h"
#include "mm_as.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_solve.h"
#include <mpi.h>
#ifdef GOMA_ENABLE_MUMPS

#include <stdlib.h>

#include "dmumps_c.h"
#include "mm_eh.h"
#include "sl_mumps.h"
#include "sl_util_structs.h"
#include "std.h"

#define JOB_END                     -2
#define JOB_INIT                    -1
#define JOB_NULL                    0
#define JOB_ANALYSIS                1
#define JOB_FACTORIZATION           2
#define JOB_SOLVE                   3
#define JOB_FACTORIZATION_AND_SOLVE 5
#define USE_COMM_WORLD              -987654

// C to Fortran macros for MUMPS
#define ICNTL(I)  icntl[(I) - 1]
#define CNTL(I)   cntl[(I) - 1]
#define INFOG(I)  infog[(I) - 1]
#define INFO(I)   info[(I) - 1]
#define RINFOG(I) rinfog[(I) - 1]
#define RINFO(I)  rinfo[(I) - 1]

struct MUMPS_data {
  DMUMPS_STRUC_C mumps;
  int N_global;
  int N;
  int nloc_rhs;
  int *irhs_loc;
  int nnz;
  int *irn;
  int *jcn;
  double *val;
  double *rhs;
  int *colgids;
  int *row_offsets;
  int *row_counts;
};

void free_solver_data(struct GomaLinearSolverData *data) {
  if (data->SolverData != NULL) {
    struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;
    if (mumps_data != NULL) {
      free(mumps_data->irhs_loc);
      free(mumps_data->irn);
      free(mumps_data->jcn);
      free(mumps_data->val);
      free(mumps_data->rhs);
      free(mumps_data->colgids);
      free(mumps_data->row_offsets);
      free(mumps_data->row_counts);
    }
    free(data->SolverData);
    data->SolverData = NULL;
  }
}

void msr_to_triplet(struct GomaLinearSolverData *data, dbl *rhs) {
  struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;

  int N = num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx];
  int local_nnz = data->nnz;

  int RowOffset;
  MPI_Scan(&N, &RowOffset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  RowOffset -= N;
  int *colgids = (int *)malloc((N + num_external_dofs[pg->imtrx]) * sizeof(int));
  for (int i = 0; i < N; i++) {
    colgids[i] = RowOffset + i + 1;
  }

  if (ProcID == 0) {
    mumps_data->row_offsets = (int *)malloc((Num_Proc + 1) * sizeof(int));
    MPI_Gather(&RowOffset, 1, MPI_INT, mumps_data->row_offsets, 1, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gather(&RowOffset, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
  }

  MPI_Allreduce(&N, &(mumps_data->N_global), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (ProcID == 0) {
    mumps_data->row_counts = (int *)malloc(Num_Proc * sizeof(int));
    mumps_data->row_offsets[Num_Proc] = mumps_data->N_global;
    for (int i = 0; i < Num_Proc - 1; i++) {
      mumps_data->row_counts[i] = mumps_data->row_offsets[i + 1] - mumps_data->row_offsets[i];
    }
    mumps_data->row_counts[Num_Proc - 1] =
        mumps_data->N_global - mumps_data->row_offsets[Num_Proc - 1];
  }

  exchange_dof_int(cx[pg->imtrx], DPI_ptr, colgids, pg->imtrx);
  mumps_data->colgids = colgids;

  mumps_data->irhs_loc = (int *)malloc(N * sizeof(int));
  mumps_data->irn = (int *)malloc(local_nnz * sizeof(int));
  mumps_data->jcn = (int *)malloc(local_nnz * sizeof(int));
  mumps_data->val = (double *)malloc(local_nnz * sizeof(double));
  mumps_data->rhs = (double *)malloc(N * sizeof(double));
  int nnz_check = 0;
  int offset = 0;
  for (int i = 0; i < N; i++) {
    int row_nnz = data->bindx[i + 1] - data->bindx[i] + 1;
    mumps_data->rhs[i] = rhs[i];
    mumps_data->irhs_loc[i] = colgids[i];

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
      mumps_data->val[offset] = data->val[j];
      offset++;
    }
    nnz_check += row_nnz;
  }

  if (nnz_check > local_nnz) {
    GOMA_EH(GOMA_ERROR, "MUMPS: nnz check failed, expected %d, got %d", local_nnz, nnz_check);
  }
  mumps_data->nnz = nnz_check;
  mumps_data->N = N;
}

void copy_rhs_and_jac(struct GomaLinearSolverData *data, dbl *rhs) {
  struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;

  int N = data->N + data->N_update;

  int offset = 0;
  for (int i = 0; i < N; i++) {
    mumps_data->rhs[i] = rhs[i];

    // diagonal
    mumps_data->val[offset] = data->val[i];
    offset++;

    // off diagonal
    for (int j = data->bindx[i]; j < data->bindx[i + 1]; j++) {
      mumps_data->val[offset] = data->val[j];
      offset++;
    }
  }
}

static goma_error check_mumps_error(DMUMPS_STRUC_C *mumps) {
  if (mumps->infog[0] < 0) {
    if (mumps->icntl[4] == 0) {
      if (ProcID == 0) {
        fprintf(stderr, "\nMUMPS Errors:\n");
        fprintf(stderr, "INFOG[1]: %d\n", mumps->infog[0]);
        fprintf(stderr, "INFOG[2]: %d\n", mumps->infog[1]);
      }
      fflush(stderr);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mumps->info[0] < 0) {
        fprintf(stderr, "P%d INFO[1]: %d\n", ProcID, mumps->info[0]);
        fprintf(stderr, "P%d INFO[2]: %d\n", ProcID, mumps->info[1]);
      }
    }
    return GOMA_ERROR;
  }
  return GOMA_SUCCESS;
}

static void set_icntl_and_cntl(DMUMPS_STRUC_C *mumps) {
  // Goma defaults
  mumps->ICNTL(4) = 0;
  if (Num_Proc == 1) {
    mumps->ICNTL(7) = 7;
    mumps->ICNTL(18) = 0;
  } else {
    mumps->ICNTL(18) = 3;
    mumps->ICNTL(20) = 11;
    mumps->ICNTL(21) = 0;
  }

  // User defined options
  for (int i = 0; i < 60; i++) {
    if (upd->solver_info->icntl_user_set[i]) {
      mumps->icntl[i] = upd->solver_info->icntl[i];
    }
  }
  for (int i = 0; i < 15; i++) {
    if (upd->solver_info->cntl_user_set[i]) {
      mumps->cntl[i] = upd->solver_info->cntl[i];
    }
  }
}

goma_error mumps_solve(struct GomaLinearSolverData *data, dbl *x, dbl *b) {
  goma_error err;
  if (data->SolverData == NULL) {

    struct MUMPS_data *mumps_data = calloc(1, sizeof(struct MUMPS_data));
    DMUMPS_STRUC_C *mumps = &mumps_data->mumps;
    data->SolverData = (void *)mumps_data;
    data->DestroySolverData = free_solver_data;

    mumps->comm_fortran = USE_COMM_WORLD;
    mumps->job = JOB_INIT;
    mumps->par = 1;
    mumps->sym = 0;

    set_icntl_and_cntl(mumps);

    mumps->job = JOB_INIT;
    dmumps_c(mumps);
    err = check_mumps_error(mumps);
    if (err != GOMA_SUCCESS) {
      free_solver_data(data);
      return err;
    }

    msr_to_triplet(data, b);
    if (Num_Proc == 1) {
      mumps->n = mumps_data->N;
      mumps->nnz = mumps_data->nnz;
      mumps->irn = mumps_data->irn;
      mumps->jcn = mumps_data->jcn;
      mumps->a = mumps_data->val;
      mumps->rhs = mumps_data->rhs;
    } else {
      mumps->n = mumps_data->N_global;
      mumps->nnz_loc = mumps_data->nnz;
      MPI_Allreduce(&(mumps->nnz_loc), &(mumps->nnz), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      mumps->irn_loc = mumps_data->irn;
      mumps->jcn_loc = mumps_data->jcn;
      mumps->a_loc = mumps_data->val;
    }

    set_icntl_and_cntl(mumps);
    mumps->job = JOB_ANALYSIS;
    dmumps_c(mumps);
    err = check_mumps_error(mumps);
    if (err != GOMA_SUCCESS) {
      free_solver_data(data);
      return err;
    }
  }
  struct MUMPS_data *mumps_data = (struct MUMPS_data *)data->SolverData;
  DMUMPS_STRUC_C *mumps = &mumps_data->mumps;

  copy_rhs_and_jac(data, b);

  if (Num_Proc == 1) {
    mumps->nrhs = 1;
    mumps->rhs_loc = mumps_data->rhs;
    mumps->nloc_rhs = mumps_data->N;
    mumps->irhs_loc = mumps_data->irhs_loc;
    mumps->lrhs_loc = mumps_data->N;
  }

  if (Num_Proc != 1) {
    mumps->nrhs = 1;
    mumps->rhs_loc = mumps_data->rhs;
    mumps->nloc_rhs = mumps_data->N;
    mumps->irhs_loc = mumps_data->irhs_loc;
    mumps->lrhs_loc = mumps_data->N;

    if (ProcID == 0)
      mumps->rhs = malloc(sizeof(double) * mumps_data->N_global);
  }
  mumps->job = JOB_FACTORIZATION_AND_SOLVE;
  set_icntl_and_cntl(mumps);
  dmumps_c(mumps);
  err = check_mumps_error(mumps);
  if (err != GOMA_SUCCESS) {
    return err;
  }

  // copy solution back to rhs
  double *rhs;
  if (Num_Proc == 1) {
    for (int i = 0; i < mumps_data->N; i++) {
      x[i] = mumps->rhs[i];
    }
  } else {
    if (ProcID == 0) {
      rhs = mumps->rhs;
    } else {
      rhs = malloc(sizeof(double) * mumps_data->N);
    }
    MPI_Scatterv(rhs, mumps_data->row_counts, mumps_data->row_offsets, MPI_DOUBLE, rhs,
                 mumps_data->N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < mumps_data->N; i++) {
      x[i] = rhs[i];
    }
  }
  if (Num_Proc != 1) {
    if (ProcID == 0) {
      free(mumps->rhs);
    } else {
      free(rhs);
    }
  }
  return GOMA_SUCCESS;
}

#endif // GOMA_ENABLE_MUMPS