/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: sl_matrix_util.c,v 5.2 2007-12-07 17:14:37 hkmoffa Exp $
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "dpi.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_unknown_map.h"
#include "mpi.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "sl_epetra_util.h"
#include "sl_matrix_util.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "std.h"
#ifdef GOMA_ENABLE_PETSC
#include "sl_petsc.h"
#endif

#define GOMA_SL_MATRIX_UTIL_C

/* canine_chaos() - return useful information about the matrix problem
 *
 * Most of this information is trivial for the serial case and a quick
 * return may be expected.
 *
 * For the parallel case, some communication is required to determine the
 * various quantities.
 *
 *
 * Created: 1997/11/04 06:26 MST pasacki@sandia.gov
 *
 * Revised:
 */

void canine_chaos(int local_order,        /* (in) */
                  int local_order_plus,   /* order including external rows (in) */
                  int local_nnz,          /* number of nonzeroes (in) */
                  int local_nnz_plus,     /* including external rows (in) */
                  int *global_order,      /* the real order (out) */
                  int *global_order_plus, /* overcounting external rows (out) */
                  int *global_nnz,        /* the strict count (out) */
                  int *global_nnz_plus)   /* overcounting external rows (out) */
{

  /*
   * Defaults work for serial mode.
   */

  *global_order = local_order;
  *global_order_plus = local_order_plus;
  *global_nnz = local_nnz;
  *global_nnz_plus = local_nnz_plus;

  /*
   * If we're parallel, then communicate and sum up each of the quantities
   * over all of the processors.
   */

#ifdef PARALLEL
  MPI_Allreduce(&local_order, global_order, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_order_plus, global_order_plus, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_nnz, global_nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_nnz_plus, global_nnz_plus, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_msr_matrix(int n, int *ija, double *a, double *x)

/*************************************************************************
 *
 *  We print out the matrix and the solution vector here
 *
 *************************************************************************/
{
  static int num_call = 0;
  int i, row, col;
  char filename[80];
  FILE *of;
#ifdef PARALLEL
#if 0
  char col_kind, row_kind;
#endif
  int global_row, global_col;
#endif

  if (n < 1) {
    GOMA_EH(GOMA_ERROR, "Bad matrix order.");
  }
  sprintf(filename, "A%d_of_%d.%d", ProcID + 1, Num_Proc, num_call);
  of = fopen(filename, "w");

  /*   fprintf(of, "# row col value \n");*/
  /* a[local_row, local_col] = A[global_row, global_col]\n"); */

  /*
   * First do the diagonal entries in a[], then do the off diagonal entries.
   */

  for (row = 0; row < n; row++) {
#ifndef PARALLEL
    fprintf(of, "%d %d %23.16e\n", row, row, a[row]);
#endif
#ifdef PARALLEL
    global_row = row;

    fprintf(of, "%d %d %23.16e %23.16e\n", global_row, global_row, a[row], x[row]);
#if 0
    if (row < num_internal_dofs[pg->imtrx] ) {
      row_kind = 'I';
    } else if (row < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]) {
      row_kind = 'B';
    } else {
      row_kind = 'E';
    }
    fprintf(of, "a[%d,%d] = %23.16e\t(A[%d,%d]) (%c,%c)\n", row, row, 
	    a[row],
	    global_row, global_row, row_kind, row_kind);
#endif
#endif
  }

  for (row = 0; row < n; row++) {
    for (i = ija[row]; i < ija[row + 1]; i++) {
      col = ija[i];
#ifndef PARALLEL
      fprintf(of, "%d %d %23.16e\n", row, col, a[i]);
#endif
#ifdef PARALLEL
      global_row = row;
      global_col = col;

      fprintf(of, "%d %d %23.16e\n", global_row, global_col, a[i]);
#if 0
      if (row < num_internal_dofs[pg->imtrx] ) {
	row_kind = 'I';
      } else if ( row < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx] ) {
	row_kind = 'B';
      } else {
	row_kind = 'E';
      }
      if (col < num_internal_dofs[pg->imtrx]) {
	col_kind = 'I';
      } else if (col < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx] ) {
	col_kind = 'B';
      } else {
	col_kind = 'E';
      }
      fprintf(of, "a[%d,%d] = %23.16e (A[%d,%d]) (%c,%c)\n", row, col, a[i],
	      global_row, global_col, row_kind, col_kind);
#endif
#endif
    }
  }
  num_call++;
  fclose(of);
  return;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void print_vbr_matrix(struct GomaLinearSolverData *ams, /* matrix info */
                      Exo_DB *exo,                      /* ptr to the whole mesh */
                      Dpi *dpi,                         /* distributed processing info */
                      int unknowns_per_node[]) {
  int i, j, k, m;
  int ii;
  int col_dof;
  int row;
  int col;
  int row_nodes;
  int col_nodes;
  int col_node = 0;
  int block;
  int block_rows, block_cols;
  int indx_kount;
  int offset;
  double a_val;

  char filename[80];

  FILE *of;

  char col_kind;

  row_nodes = dpi->num_internal_nodes + dpi->num_boundary_nodes;
  col_nodes = dpi->num_universe_nodes;

  if (ams->nnz < 1) {
    GOMA_EH(GOMA_ERROR, "Bad matrix order.");
  }

  sprintf(filename, "A%d_of_%d", ProcID + 1, Num_Proc);

  of = fopen(filename, "w");

  fprintf(of, "%d\n", ams->nnz);

  indx_kount = 0;

  for (i = 0; i < row_nodes; i++) { /* loop over row blocks */
    for (j = ams->bpntr[i]; j < ams->bpntr[i + 1]; j++) {

      block = ams->bindx[j];

      if (Num_Proc > 1) {
        /* account for possibility of different  */
        /*   number of unknowns per element      */
        col_dof = 0;
        for (ii = 0; ii < col_nodes; ii++) {
          col_node = ii;
          if (col_dof == ams->cpntr[block])
            break;
          col_dof += unknowns_per_node[ii];
        }

        row = Nodes[i]->First_Unknown[pg->imtrx];
        col = Nodes[col_node]->First_Unknown[pg->imtrx];
      } else {
        row = ams->rpntr[i];
        col = ams->cpntr[block];
      }

      block_rows = ams->rpntr[i + 1] - ams->rpntr[i];
      block_cols = ams->cpntr[block + 1] - ams->cpntr[block];

      offset = 0;

      for (k = 0; k < block_cols; k++) {
        for (m = 0; m < block_rows; m++) {
          a_val = ams->val[ams->indx[indx_kount] + offset];
          offset++;

          if (Num_Proc > 1) {
            if (col < num_internal_dofs[pg->imtrx]) {
              col_kind = 'I';
            } else if (col < num_internal_dofs[pg->imtrx] + num_boundary_dofs[pg->imtrx]) {
              col_kind = 'B';
            } else {
              col_kind = 'E';
            }
          } else {
            col_kind = 'I';
          }
          fprintf(of, "%d %2d %2d %9f %c\n", ProcID, row + m, col, a_val, col_kind);
        }
        col++;
      }
      indx_kount++;
    }
  }
  return;
}

/* consign this chunk of code to reading reference material only and
 * do not compile it in to bulk up the code.
 */

#if 0


/*      This function returns matrix and problem statistics:
 *
 *  input
 *  -----
 *	a[]	-	Vector of non-zero entries in the coefficient matrix
 *			for the current processor.
 *	ija[]	-	Vector of column pointer information, containing the 
 *			positional information for the coefficient matrix,
 *			specific to the current processor
 *	N	-	Number of unknowns on the current processor
 *
 *  output
 *  ------
 *	*nzeroes	Number of nonzeroes in a matrix (on this proc).
 *
 *	*gnzeros	Total number of non-zero entries in the global matrix
 *
 *	*gN		Total number of global unknowns in the problem.
 *
 * Notes:
 *	This needs to be modified to reflect the local overage issue.
 */


void 
matrix_stats (double a[],
	      int ija[], 
	      int N, 
	      int *local_num_nonzeroes, 
	      int *global_num_nonzeroes,
	      int *global_num_unknowns)
{

  *local_num_nonzeroes = ija[N]-1;
  
  /*
   * Determine the order of the global system. Easy if we're serial, not
   * too hard if we're parallel...
   */

#ifdef PARALLEL
  MPI_Allreduce(&N, global_num_unknowns, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
#ifndef PARALLEL
  *global_num_unknowns = N;
#endif


  /*
   * How many entries in this sparse matrix? Sum up the total for all
   * processors if applicable...
   */

#ifdef PARALLEL
  MPI_Allreduce(global_num_nonzeroes, local_num_nonzeroes, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
#endif
#ifndef PARALLEL
  *global_num_nonzeroes = *local_num_nonzeroes;
#endif

  return;
}

/*
	this routine returns maximum matrix norm for 
	C-MSR matrix a: this is a parallel version
	for a distributed matrix.


        Author:         John N. Shadid Div 1421 SNL
        Date:           10/4/1990
        revised:        10/4/1990

        Paramter list:

        N   ==  order of linear system
        a   ==  matrix A in sparse format (C-MSR)
        ija ==  pointers to nonzeros of A (C-MSR)

*/
double
pmax_matrix_norm ( int N,
                   double a[],
                   int ija[]  )
{
/* LOCAL VARIABLES */
  int 		 k;
  register int 	 j;
  int 		 j_last, ija_row, irow;
  double         row_sum, row_max;

/* EXTERNAL FUNCTIONS and PROTOTYPES */
extern double
gmax_double
   ( double,              /* var  */
            int,                 /* me   */
            int  );             /* dim  */


  row_max = 0.0;
  for( irow = 0; irow < N; irow++){

    /* compute diagonal contribution */

    row_sum  = fabs(a[irow]);

    /* nonzero off diagonal contibution */

    j_last = ija[irow+1] - ija[irow];
    ija_row = ija[irow];

    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      row_sum += fabs(a[k]);
    }
    row_max = MAX( row_sum, row_max);
  }
  row_max = gmax_double (row_max, ProcID, Dim);
  return(row_max);
} /* END of routine pmax_matrix_norm */
#endif

/******************************************************************************/

/*
 *
 *      routine to row sum scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire
 *      matrix problem Ax = b and return scaling vector
 *
 *
 *              John N. Shadid Div 1421 SNL
 *              date:    1/20/90
 *              revised: 1/6/95 Richard Cairncross
 *
 *
 *      N   ==  order of linear system
 *      a   ==  matrix A in sparse format (C-PCGPAK)
 *      ija ==  pointers to nonzeros of A (C-PCGPAK)
 *      b   ==  right hand side(rhs) of matrix problem
 *      scale == vector of scaling
 *
 *
 */

void row_sum_scaling_scale(struct GomaLinearSolverData *ams, double b[], double scale[]) {
  if (strcmp(Matrix_Format, "msr") == 0) {
    row_sum_scale_MSR(ams->npu, ams->val, ams->bindx, b, scale);
  } else if (strcmp(Matrix_Format, "vbr") == 0) {
    row_sum_scale_VBR(ams->npn, ams->val, ams->bpntr, ams->bindx, ams->indx, ams->rpntr, ams->cpntr,
                      b, scale);
  } else if (strcmp(Matrix_Format, "epetra") == 0) {
    row_sum_scale_epetra(ams, b, scale);
#ifdef GOMA_ENABLE_PETSC
#if PETSC_USE_COMPLEX
  } else if (strcmp(Matrix_Format, "petsc_complex") == 0) {
    // Skip
#else
  } else if (strcmp(Matrix_Format, "petsc") == 0) {
    petsc_scale_matrix(ams, b, scale);
#endif
#endif
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown sparse matrix format");
  }
}

/* end of row_sum_scaling */

void row_sum_scaling_scale_AC(double **cAC, double **dAC, double *gAC, int nAC) {
  int iAC, jAC, i;
  int numProcUnknowns = NumUnknowns[pg->imtrx] + NumExtUnknowns[pg->imtrx];
  double *row_sums, *local_sums, row_sum_inv;

  row_sums = alloc_dbl_1(nAC, 0.0);
  for (iAC = 0; iAC < nAC; iAC++) {
    row_sums[iAC] = 0.;
    for (jAC = 0; jAC < nAC; jAC++)
      row_sums[iAC] += fabs(dAC[iAC][jAC]);
    for (i = 0; i < numProcUnknowns; i++)
      row_sums[iAC] += fabs(cAC[iAC][i]);
  }

#ifdef PARALLEL
  if (Num_Proc > 1) {
    local_sums = alloc_dbl_1(nAC, 0.0);
    for (iAC = 0; iAC < nAC; iAC++)
      local_sums[iAC] = row_sums[iAC];
    MPI_Allreduce(local_sums, row_sums, nAC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    safe_free(local_sums);
  }
#endif /* PARALLEL */

  for (iAC = 0; iAC < nAC; iAC++) {
    if (fabs(row_sums[iAC]) < DBL_SMALL) {
      GOMA_EH(GOMA_ERROR, "Zero row sum scale for AC!\n");
    }
    row_sum_inv = 1. / row_sums[iAC];

    gAC[iAC] *= row_sum_inv;
    for (jAC = 0; jAC < nAC; jAC++)
      dAC[iAC][jAC] *= row_sum_inv;
    for (i = 0; i < numProcUnknowns; i++)
      cAC[iAC][i] *= row_sum_inv;
  }

  safe_free(row_sums);
}

void row_sum_scale_MSR(int N, double a[], int ija[], double b[], double scale[]) {
  int j, k, irow;
  int j_last, ija_row;
  double row_sum = 0.0;
#ifdef DEBUG_ZEROROW
  int etag = FALSE, inode, i_Var_Desc, i_offset, idof;
  VARIABLE_DESCRIPTION_STRUCT *vd = NULL;
#endif

  /* index through rows of matrix */

  for (irow = 0; irow < N; irow++) {
    /* scale nonzero off diagonal elements */
    j_last = ija[irow + 1] - ija[irow];
    ija_row = ija[irow];
    row_sum = fabs(a[irow]);
    for (j = 0; j < j_last; j++) {
      k = ija_row + j;
      row_sum += fabs(a[k]);
    }

#ifdef DEBUG_ZEROROW
    if (fabs(row_sum) == 0.0) {
      if (dofname) {
        printf("row_sum_scaling_scale ERROR: Row %d is zero, dofname = %s\n", irow,
               dofname[pg->imtrx][irow]);
      } else {
        printf("row_sum_scaling_scale ERROR: Row %d is zero, dofname = unknown\n", irow);
        vd = Index_Solution_Inv(irow, &inode, &i_Var_Desc, &i_offset, &idof);
        printf("\t var_type = %d, matid = %d, Node = %d\n", vd->Variable_Type, vd->MatID, inode);
      }
      fflush(stdout);
      row_sum = 1.0;
      etag = TRUE;
    }
#endif
    /*
     *  If a is nonzero, we make the diagonal component
     *  of the matrix positive. However, if the diagonal
     *  is basically zero, let's not change the sign of
     *  the other values. We do this so as to not have the
     *  matrix change values arbitrarily due to round off error.
     *  HKM -> Found that this change enhances the ability
     *         of the jacobian checker to compare matrices.
     */
    if (fabs(a[irow]) > 1.0E-200) {
      row_sum = row_sum * SGN(a[irow]);
    }

    scale[irow] = row_sum;

    /* scale elements & rhs*/
    if (row_sum == 0.0) {
#define KEEP_GOING_ON_ZERO_ROW_SUM 1
#define WARNING_ON_ZERO_ROW_SUM    1
#if KEEP_GOING_ON_ZERO_ROW_SUM
#if WARNING_ON_ZERO_ROW_SUM
      int i_Var_Desc, i_offset, idof;
      double x[3] = {0., 0., 0.};
      int inode = 0;
      VARIABLE_DESCRIPTION_STRUCT *vd;
      GOMA_WH(-1, "row_sum_scale_MSR ERROR: row_sum = 0.0,");
      vd = Index_Solution_Inv(irow, &inode, &i_Var_Desc, &i_offset, &idof, pg->imtrx);
      x[0] = Coor[0][inode];
      x[1] = Coor[1][inode];
      if (pd_glob[0]->Num_Dim == 3)
        x[2] = Coor[2][inode];
      if (dofname) {
        fprintf(stderr, "row_sum_scaling_scale ERROR: Row %d is zero, dofname = %s, x=(%g,%g,%g)\n",
                irow, dofname[pg->imtrx][irow], x[0], x[1], x[2]);
      } else {
        fprintf(stderr, "row_sum_scaling_scale ERROR: Row %d is zero, dofname = unknown\n", irow);
        fprintf(stderr, "\t var_type = %d, matid = %d, Node = %d, x=(%g,%g,%g)\n",
                vd->Variable_Type, vd->MatID, inode, x[0], x[1], x[2]);
      }
#else
      a[irow] = 1.;
      b[irow] = 0.;
      row_sum = 1.;
#endif
#else
      GOMA_EH(GOMA_ERROR, "row_sum_scale_MSR ERROR: row_sum = 0.0,");
#endif
    }

    b[irow] = b[irow] / row_sum;
    a[irow] = a[irow] / row_sum;
    for (j = 0; j < j_last; j++) {
      k = ija_row + j;
      a[k] /= row_sum;
    }
  }
#ifdef DEBUG_ZEROROW
  if (etag)
    exit(-1);
#endif
} /* END of routine row_sum_scaling_scale */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void row_sum_scale_VBR(int N,
                       double *a,
                       int *bpntr,
                       int *bindx,
                       int *indx,
                       int *rpntr,
                       int *cpntr,
                       double *b,
                       double *scale) {

  int index;
  int blk_rows;
  int I, J, K;
  int i, j;
  int ib, jb;
  double sign;

  for (I = 0; I < N; I++) {
    blk_rows = rpntr[I + 1] - rpntr[I];

    memset(&(scale[rpntr[I]]), 0, blk_rows * sizeof(double));

    for (i = rpntr[I], ib = 0; i < rpntr[I + 1]; i++, ib++) {
      sign = 0.0;

      for (K = bpntr[I]; K < bpntr[I + 1]; K++) {
        J = bindx[K];

        index = indx[K] + ib;

        for (j = cpntr[J], jb = 0; j < cpntr[J + 1]; j++, jb++) {
          scale[i] += fabs(a[index + jb * blk_rows]);

          if (j == i)
            sign = SGN(a[index + jb * blk_rows]);
        }
      }

      if (sign != 1.0 && sign != -1.0)
        GOMA_EH(GOMA_ERROR, "Can't find diagonal sign in row_scale_VBR");

      scale[i] *= sign;

      b[i] /= scale[i];

      for (K = bpntr[I]; K < bpntr[I + 1]; K++) {
        J = bindx[K];

        index = indx[K] + ib;

        for (j = cpntr[J], jb = 0; j < cpntr[J + 1]; j++, jb++) {
          a[index + jb * blk_rows] /= scale[i];
        }
      }
    }
  }
} /* END of routine row_sum_scale_VBR */

void row_sum_scale_epetra(struct GomaLinearSolverData *ams, double *b, double *scale) {
  EpetraRowSumScale(ams, b, scale);
}

/******************************************************************************/
/******************************************************************************/ /******************************************************************************/
/*
 * This routine scales a matrix using a previously-calculated
 * scale vector and an optional constant factor. It is useful
 * for scaling a mass matrix after the corresponding Jacobian
 * has been scaled, and allows the sign to be changed.
 */

void matrix_scaling(struct GomaLinearSolverData *ams, double *a, double factor, double *scale) {

  int index;
  int blk_rows;
  int I, J, K;
  int i, j;
  int ib, jb;
  int *bpntr = ams->bpntr;
  int *bindx = ams->bindx;
  int *indx = ams->indx;
  int *rpntr = ams->rpntr;
  int *cpntr = ams->cpntr;

  if (strcmp(Matrix_Format, "msr") == 0) {
    for (i = 0; i < ams->npu; i++) {
      a[i] *= factor / scale[i];
      for (j = bindx[i]; j < bindx[i + 1]; j++) {
        a[j] *= factor / scale[i];
      }
    }
  } else if (strcmp(Matrix_Format, "vbr") == 0) {
    for (I = 0; I < ams->npn; I++) {
      blk_rows = rpntr[I + 1] - rpntr[I];
      for (i = rpntr[I], ib = 0; i < rpntr[I + 1]; i++, ib++) {
        for (K = bpntr[I]; K < bpntr[I + 1]; K++) {
          J = bindx[K];
          index = indx[K] + ib;
          for (j = cpntr[J], jb = 0; j < cpntr[J + 1]; j++, jb++) {
            a[index + jb * blk_rows] *= factor / scale[i];
          }
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Matrix format must be MSR or VBR!");
  }
} /* END of routine matrix_scaling */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 *
 *      routine to scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire

 *      matrix problem Ax = b using vector scale
 *
 *
 *              John N. Shadid Div 1421 SNL
 *              date:    1/20/90
 *              revised: 1/6/95 Richard Cairncross
 *
 *
 *      N   ==  order of linear system
 *      a   ==  matrix A in sparse format (C-PCGPAK)
 *      ija ==  pointers to nonzeros of A (C-PCGPAK)
 *      b   ==  right hand side(rhs) of matrix problem
 *
 *
*/
void row_scaling(const int N, double a[], int ija[], double b[], double scale[])

{
  /* LOCAL VARIABLES */
  register int j, k, irow;
  int j_last, ija_row;
  double row_sum;

  /* index through rows of matrix */

  for (irow = 0; irow < N; irow++) {

    /* scale nonzero off diagonal elements */

    j_last = ija[irow + 1] - ija[irow];
    ija_row = ija[irow];

    /* scale elements & rhs*/
    row_sum = scale[irow];

    b[irow] = b[irow] / row_sum;
    a[irow] = a[irow] / row_sum;
    for (j = 0; j < j_last; j++) {
      k = ija_row + j;
      a[k] /= row_sum;
    }
  }
} /* END of routine row_scaling */
/******************************************************************************/

/*
 *
 */

void row_sum_scaling(struct GomaLinearSolverData *ams, double b[]) {
  double *scale;

  scale = (double *)smalloc(ams->npu * sizeof(double));

  row_sum_scaling_scale(ams, b, scale);

  safe_free((void *)scale);

} /* END of routine row_sum_scaling */

/*
 *
 *      routine to scale vector with previously determined scale
 *      coming from, for instance, row_sum_scaling_scale above.
 *
 *              P. R. Schunk Div 9111 SNL
 *              date:    5/30/97
 *
 *      N     ==  order of linear system
 *      a     ==  matrix A in sparse format (C-PCGPAK)
 *      ija   ==  pointers to nonzeros of A (C-PCGPAK)
 *      b     ==  right hand side(rhs) of matrix problem
 *      scale == scale vector (each component for each row) determined
 *               elsewhere
 *
 *
 */
void vector_scaling(const int N, double b[], double scale[]) {
  register int irow;
  double row_sum;

  /* index through rows of matrix */

  for (irow = 0; irow < N; irow++) {

    /* scale elements & rhs*/
    row_sum = scale[irow];

    b[irow] = b[irow] / row_sum;
  }
}

/**
 * Check to see if the linear solver is compatible with
 * the matrix format
 *
 * Currently only checks epetra matrix format
 *
 * @return 0 if compatible, -1 if not
 */
int check_compatible_solver(void) {
  if (strcmp(Matrix_Format, "epetra") == 0) {
    switch (Linear_Solver) {
    case AZTECOO:
    case AMESOS:
    case STRATIMIKOS:
      return GOMA_SUCCESS;
    default:
      return GOMA_ERROR;
    }
  }

#ifdef GOMA_ENABLE_PETSC
  if (strcmp(Matrix_Format, "petsc") == 0) {
    switch (Linear_Solver) {
    case PETSC_SOLVER:
      return GOMA_SUCCESS;
    default:
      return GOMA_ERROR;
    }
  }
  if (strcmp(Matrix_Format, "petsc_complex") == 0) {
    switch (Linear_Solver) {
    case PETSC_COMPLEX_SOLVER:
      return GOMA_SUCCESS;
    default:
      return GOMA_ERROR;
    }
  }
#endif

  return GOMA_ERROR;
}

/* unused routine, but keep for reference reading... */
#if 0

/*
 *
 *      routine to symmetricaly diagonaly scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire
 *      matrix problem Ax = b, the routine sym_rescale must be used
 *	to transform solution back to recover soution to original problem
 *
 *
 *              John N. Shadid Div 1421 SNL
 *              date:    1/20/90
 *              revised: 1/20/90
 *
 *
 *      N       ==  order of linear system
 *      a       ==  matrix A in sparse format (C-PCGPAK)
 *      ija     ==  pointers to nonzeros of A (C-PCGPAK)
 *      b       ==  right hand side(rhs) of matrix problem
 *	d_mhalf == vector containing one over the square root of
 *		   the diagonals of the original problem
 *
 *
*/
void
sym_diagonal_scaling ( int N,
                       double a[],
                       int ija[],
                       double b[],
                       double d_mhalf[] )

{
/* LOCAL VARIABLES */
  register int j,k,irow;
  int j_last,ija_row;
  double sign;

/* EXTERNAL FUNCTIONS and PROTOTYPES */

  /* do left diagonal scaling of matrix and rhs */

  /* index through rows of matrix */
  
  for( irow = 0; irow < N; irow++){

    /* scale nonzero off diagonal elements */

    j_last = ija[irow+1] - ija[irow];
    ija_row = ija[irow];
	
    sign = SGN(a[irow]);
    d_mhalf[irow]= 1./sqrt(fabs(a[irow]));

    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      a[k] *= sign*d_mhalf[irow]; 
    }

    /* scale diagonal elements */
	
    b[irow] *= sign*d_mhalf[irow];
    a[irow] *= sign*d_mhalf[irow];
  }

  /* do right diagonal scaling */

  /* index through rows of matrix */

  for( irow = 0; irow < N; irow++){
    /* scale diagonal elements */
	
    a[irow] *= d_mhalf[irow];

    /* scale nonzero off diagonal elements */

    j_last = ija[irow+1] - ija[irow];
    ija_row = ija[irow];
    
    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      a[k] *= d_mhalf[ija[k]]; 
    }
  }
} /* END of routine sym_diagonal_scaling */
/******************************************************************************/

/*
 * 
 *      routine to symmetricaly row sum scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire
 *      matrix problem Ax = b, the routine sym_rescale must be used  
 *      to transform solution back to recover soution to original problem 
 * 
 * 
 *              John N. Shadid Div 1421 SNL 
 *              date:    1/20/90
 *              revised: 1/20/90
 *
 * 
 *      N       ==  order of linear system 
 *      a       ==  matrix A in sparse format (C-PCGPAK) 
 *      ija     ==  pointers to nonzeros of A (C-PCGPAK) 
 *      b       ==  right hand side(rhs) of matrix problem 
 *      d_mhalf == vector containing one over the square root of 
 *                 the row sum of the original problem 
 * 
 * 
*/

void
kr_sym_row_sum_scaling ( int N,
                         double a[],
                         int ija[],
                         double b[],
                         double d_mhalf[] )

{
/* LOCAL VARIABLES */
  register int 	j, k, irow;
  int 		j_last, ija_row;
  double 	sign, row_sum;

/* EXTERNAL FUNCTIONS and PROTOTYPES */


  /* do left diagonal scaling of matrix and rhs */

  /* index through rows of matrix */

  for( irow = 0; irow < N; irow++){

    j_last = ija[irow+1] - ija[irow];
    ija_row = ija[irow];

    /* obtain row sum factor */

    row_sum = fabs(a[irow]);
    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      row_sum += fabs(a[k]);
    }

    sign = SGN(a[irow]);
    d_mhalf[irow]= 1./sqrt(row_sum);

    /* scale nonzero off diagonal elements */

    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      a[k] *= sign*d_mhalf[irow]; 
    }

    /* scale diagonal elements */
	
    b[irow] *= sign*d_mhalf[irow];
    a[irow] *= sign*d_mhalf[irow];
  }


  /* do right diagonal scaling */
  
  /* index through rows of matrix */

  for( irow = 0; irow < N; irow++){

    /* scale diagonal elements */
	
    a[irow] *= d_mhalf[irow];

    /* scale nonzero off diagonal elements */

    j_last = ija[irow+1] - ija[irow];
    ija_row = ija[irow];
    
    for( j = 0; j < j_last; j++){
      k = ija_row + j;
      a[k] *= d_mhalf[ija[k]]; 
    }
  }
} /* END of routine kr_sym_row_sum_scaling */
/******************************************************************************/

/*
 * 
 *      routine to symmetricaly diagonaly scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire
 *      matrix problem Ax = b, the routine sym_rescale must be used  
 *      to transform solution back to recover soution to original problem 
 * 
 * 
 *              John N. Shadid Div 1421 SNL 
 *              date:    1/20/90
 *              revised: 1/20/90
 *
 * 
 *      N       ==  order of linear system 
 *      x       ==  on input solution of transformed problem 
 *		    on output the solution to the original problem
 *      d_mhalf == vector containing one over the square root of 
 *                 the scaling factors of the original problem 
 *		   (see above symmetric scalings)
 * 
 * 
*/

void
x_scale ( int N,
          double x[],
          double d_mhalf[] )
{
/* LOCAL VARIABLE */
  int i;

  for(i = 0; i < N; i++){
    x[i] = x[i]/d_mhalf[i];
  }
} /* END of routine x_scale */
/******************************************************************************/

/*
 * 
 *      routine to symmetricaly diagonaly scale sparse matrix problem in
 *      sparse C-PCGPAK notation; Note: this scales the entire
 *      matrix problem Ax = b, the routine sym_rescale must be used  
 *      to transform solution back to recover soution to original problem 
 * 
 * 
 *              John N. Shadid Div 1421 SNL 
 *              date:    1/20/90
 *              revised: 1/20/90
 *
 * 
 *      N       ==  order of linear system 
 *      x       ==  on input solution of transformed problem 
 *		    on output the solution to the original problem
 *      d_mhalf == vector containing one over the square root of 
 *                 the scaling factors of the original problem 
 *		   (see above symmetric scalings)
 * 
 * 
 */
void
sym_rescale ( int N,
              double x[],
              double d_mhalf[] )
{
/* LOCAL VARIABLE */
  int i;

  for(i = 0; i < N; i++){
    x[i] = x[i]*d_mhalf[i];
  }
} /* END of routine sym_rescale */

#endif

/******************************************************************************/
/* END of file sl_matrix_util.c */
/******************************************************************************/
