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
 *$Id: sl_lu.c,v 5.2 2007-12-07 17:14:37 hkmoffa Exp $
 */

#include <stdlib.h>

#include "mm_as.h"
#include "mm_as_structs.h"
#include "rf_fem_const.h"
#include "sl_lu.h"

#define GOMA_SL_LU_C

#ifdef GOMA_ENABLE_SPARSE
#include "spMatrix.h"
/*
 * Uncomment the following line to write out the RHS, a, and ija vectors
 * to an ASCII file...
 */

/* #define DUMP_CMSR_MATRIX */

int first_time = TRUE;

static int call = 0;

void lu(const int N,
        const int NExt,
        const int M,
        double a[],
        int ija[],
        double x[],
        const int factor_flag) {
  static char *matrix;
  static int prev_matrix = 0;
  int error, type;
  int j, i, n, k, nzeros, ija_col;
  static spREAL **element, *b;
#ifdef MATRIX_STATISTICS
  spREAL norm;
#endif

  if (pg->imtrx != prev_matrix) {
    if (first_time == FALSE) {
      spDestroy(matrix);
    }
    first_time = TRUE;
    prev_matrix = pg->imtrx;
  }

  call++;

  /* allocate rhs and copy */

  b = (spREAL *)malloc((N + NExt + 1) * sizeof(spREAL));
  b[0] = 0.0;
  for (i = 0; i < N + NExt; i++)
    b[i + 1] = x[i];

  /* factor matrix if values have changed since last call */
  if (factor_flag <= 2) {

    if (first_time == TRUE && factor_flag == 1) {
      /* N = order of real matrix */
      type = 0;
      matrix = spCreate(N + NExt, type, &error);

      element = (spREAL **)malloc((M + NExt + 2) * sizeof(spREAL *));
    }

    /* n = 1 produce graph of matrix, n = 2 fill matrix */
    for (n = 1; n <= 2; n++) { /* define and fill */

      /* copy MSR matrix into linked list format */
      k = 1;
      for (i = 0; i < N; i++) { /* lop through rows */
        nzeros = ija[i + 1] - ija[i];
        ija_col = ija[i];

        if (n == 1) { /* create nonzero pattern */

          /* diagonal element */
          element[k++] = spGetElement(matrix, i + 1, i + 1);

          /* nonzero off diagonal elements */
          for (j = 0; j < nzeros; j++) {
            element[k++] = spGetElement(matrix, i + 1, ija[ija_col++] + 1);
          }
        } else { /* fill in numerical values */

          /* diagonal element */
          spADD_REAL_ELEMENT(element[k++], a[i]);

          /* nonzero off diagonal elements */
          for (j = 0; j < nzeros; j++) {
            spADD_REAL_ELEMENT(element[k++], a[ija_col++]);
          }
        }
      } /* end row loop */

      /* fill in Dirchlet equations for the interface unknowns */
      for (i = N; i < N + NExt; i++) {
        if (n == 1) { /* create nonzero pattern */
          element[k++] = spGetElement(matrix, i + 1, i + 1);
        } else { /* fill in numerical values */
          spADD_REAL_ELEMENT(element[k++], 1.0);
        }
      }

      if (n == 1)
        spClear(matrix); /* zero entries in graph */

    } /* end of define and fill */

#ifdef MATRIX_STATISTICS
    norm = spNorm(matrix);
#endif
    if (first_time == TRUE && factor_flag == 1) {

      /*         spFileMatrix(matrix,"matrix_file","channel",0,1,1); */
      error = spOrderAndFactor(matrix, b, (spREAL)-1, (spREAL)0, (int)1);
      first_time = FALSE;

    } else {

      spFactor(matrix);
    }

  } /* alocate matrix and factor matrix */

  /* solve lu factored system */

  spSolve(matrix, b, x - 1);

#ifdef MATRIX_STATISTICS
  if (call == 4) {
    lustat(N, M, a, ija, x, norm, matrix);
  }
#endif

  /*
    for(i = 0; i < N; i++) printf("i: %d\t%f\t%f\n",i,x[i],b[i+1]);
  */
  free(b);
  /*  spDestroy(matrix);*/

} /* END of routine lu */
#else // GOMA_ENABLE_SPARSE
#include "mm_eh.h"

void lu(const int N,
        const int NExt,
        const int M,
        double a[],
        int ija[],
        double x[],
        const int factor_flag) {
  GOMA_EH(GOMA_ERROR, "Goma not configured with sparse support");
}
#endif // GOMA_ENABLE_SPARSE
/*****************************************************************************/
/* END of file sl_lu.c */
/*****************************************************************************/
