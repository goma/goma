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
        UMFPACK CALLING ROUTINE

        ARGUMENTS:

        first  = +1   : FIRST CALL
        first !=  0   : NOT FIRST CALL
        first  = -1   : LAST CALL

        fact_optn =-2 : LOAD MATRIX AND ANALYSIS/DECOMPOSITION ONLY
                        NO BACK SUBSTITUTION
        fact_optn =-1 : LOAD MATRIX AND DECOMPOSITION USING PAST ANALYSIS
                        NO BACK SUBSTITUTION
        fact_optn = 0 : LOAD MATRIX AND ANALYSIS/DECOMPOSITION
                        WITH BACK SUBSTITUTION
        fact_optn = 1 : LOAD MATRIX AND DECOMPOSITION USING PAST ANALYSIS
                        WITH BACK SUBSTITUTION
        fact_optn > 2 : BACK SUBSTITUTION ONLY

        matr_form = 0 : INPUT MATRIX IN COORDINATE FORMAT
        matr_form = 1 : INPUT MATRIX IN MSR FORMAT
        matr_form = 2 : INPUT MATRIX IN CSR FORMAT

        n             : DIMENSION OF SYSTEM
        nnz           : NUMBER OF NON-ZEROES IN MATRIX (LENGTH OF a[])
        row           : ROW DATA (IF matr_form = 0,2)
        col           : COL DATA (OR ija IF matr_form = 1)
        a             : VALUES IN MATRIX
        b             : RHS
        x             : SOLUTION VECTOR

        BY IAN GATES AUG 8 1997

        MODIFIED:

        IDG SEPT 28 1997
        PUT IN FACTORIZE ONLY FLAG

        IDG OCT 18 1997
        LOAD MATRIX IN FORTRAN WRAPPER INTO IMEM/XMEM RATHER THAN HERE
        NEEDED FOR FAILED FACTORIZATION BASED ON OLD ANALYSIS

        IDG OCT 19 1997
        RE-ALLOCATE MEMORY AS NEEDED WITH CALL BACK TO FACTORIZER

        DRN AUG 26 2003
        NO MORE STATIC MEMORY ALLOCATION, WOOHOO!
        UPGRADE TO UMFPACKv4.1

*/

#include <stdio.h>

#ifdef GOMA_ENABLE_UMFPACK
#include <umfpack.h>
#endif

#define GOMA_SL_UMF_C
#include "mm_eh.h"
#include "sl_auxutil.h"
#include "sl_umf.h"

/* how many different linear systems might UMF be used for? */
#ifndef UMF_MAX_SYSTEMS
#define UMF_MAX_SYSTEMS 20
#endif

#ifdef GOMA_ENABLE_UMFPACK

#include "umfpack.h"

/* This returns an integer identifier that should be unique to your
 * system.  There were problems with UMF mixing up systems becuase it
 * would identify unique systems just by its size.
 *
 * This unique identifier is passed in as system_id.  If you're
 * creating the matrix for the first time, then you should pass in a
 * -1, otherwise you should pass in the returned value from SL_UMF
 * when you created your system.
 *
 * Note that we don't do this very intelligently.  We simply use
 * indices sequentially.  There is no mechanism to allow re-use.
 */
int SL_UMF(int system_id,
           int *first,
           int *fact_optn,
           int *matr_form,
           int *nj,
           int *nnz_j,
           int *row,
           int *col,
           double *a,
           double *b,
           double *x) {
  /* Static struct holds all linear systems also keep track of number
   * of systems we have set up */
  static struct UMF_Linear_Solver_System ums_a[UMF_MAX_SYSTEMS];
  static int number_systems = 0;

  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];

  struct UMF_Linear_Solver_System *ums = 0; /* pointer to current system */

  int i, j, k, umf_option = 0;
  int hit_diag, err;

  for (i = 0; i < UMFPACK_CONTROL; i++) {
    Control[i] = 0;
  }

  for (i = 0; i < UMFPACK_INFO; i++) {
    Info[i] = 0;
  }

#ifdef DEBUG_SL_UMF
  fprintf(stderr, "SL_UMF: system_id = %d, *first = %d, *fact_optn = %d\n", system_id, *first,
          *fact_optn);
#endif

  /* MEMORY */
  switch (*first) {
  case 1:
    /* If *first == 1, then we're creating a new matrix. */

    /* If system_id isn't -1, then we're probably making some sort of mistake... */
    if (system_id != -1)
      GOMA_EH(GOMA_ERROR, "Entered SL_UMF with *first == 1, but system_id != -1");
    /* If we've already gone through all of our slots, get out. */
    if (number_systems == UMF_MAX_SYSTEMS)
      GOMA_EH(GOMA_ERROR, "Already created UMF_MAX_SYSTEMS systems");

    system_id = number_systems;
    ums = &ums_a[number_systems++];
    ums->n = *nj;
    ums->nnz = *nnz_j;

    /* MATRIX VECTORS */
    ums->ap = Ivector_birth(ums->n + 1);
    ums->ai = Ivector_birth(ums->nnz);
    ums->ax = Dvector_birth(ums->nnz);

    /* MSR needs extra allocation for A-transpose */
    ums->atp = NULL;
    ums->ati = NULL;
    ums->atx = NULL;
    if (*matr_form == 1) {
      ums->atp = Ivector_birth(ums->n + 1);
      ums->ati = Ivector_birth(ums->nnz);
      ums->atx = Dvector_birth(ums->nnz);
    }

    break;

  case 0:
    /* If *first == 0, then we want to just reuse a previously created
     * system. */

    /* system_id should have the appropriate identifier. */
    if (system_id == -1)
      GOMA_EH(GOMA_ERROR, "Conflicting orders: system_id == -1 and *first != 1");
    if (system_id < 0 || system_id >= UMF_MAX_SYSTEMS)
      GOMA_EH(GOMA_ERROR, "Index out of range: system_id");

    /* Grab the hopeful system. */
    ums = &ums_a[system_id];

    /* Run through some sanity checks to help ensure we're dealing
     * with the correct system. */
    if (ums->n != *nj || ums->nnz != *nnz_j)
      GOMA_EH(GOMA_ERROR, "Tried to access a bad system");
    break;

  case -1:
    /* If *first == -1, then we want to free space. */

    /* system_id should have the appropriate identifier. */
    if (system_id == -1)
      GOMA_EH(GOMA_ERROR, "Conflicting orders: system_id == -1 and *first != 1");
    if (system_id < 0 || system_id >= UMF_MAX_SYSTEMS)
      GOMA_EH(GOMA_ERROR, "Index out of range: system_id");

    ums = &ums_a[system_id];
    /* Run through some sanity checks to help ensure we're dealing
     * with the correct system. */
    if (ums->n != *nj || ums->nnz != *nnz_j)
      GOMA_EH(GOMA_ERROR, "Tried to free a bad system");

    umfpack_di_free_symbolic(&ums->symbolic);
    ums->symbolic = NULL;
    umfpack_di_free_numeric(&ums->numeric);
    ums->numeric = NULL;
    Ivector_death(ums->ap, ums->n + 1);
    Ivector_death(ums->ai, ums->nnz);
    Dvector_death(ums->ax, ums->nnz);

    if (ums->atp != NULL) {
      Ivector_death(ums->atp, ums->n + 1);
      Ivector_death(ums->ati, ums->nnz);
      Dvector_death(ums->atx, ums->nnz);
    }

    /* MMH: The fix that changed the world... */
    ums->n = 0;
    ums->nnz = 0;

    /* So things break later in case we actually use the return value
     * after deallocating space. */
    system_id = -1;

    break;
  }

  /* CONVERT MSR FORMAT TO MATLAB FORMAT IF NEEDED */
  if (abs(*fact_optn) < 3) {
    switch (*matr_form) {
    case 0: /* COORDINATE FORMAT */
      umfpack_di_triplet_to_col(ums->n, ums->n, ums->nnz, row, col, a, ums->ap, ums->ai, ums->ax,
                                NULL);
      break;
    case 1: /* MSR FORMAT */
      /* Note: MSR is row-oriented and UMF wants column-oriented data.
         So, assemble A-transpose in UMF format, and use umf utility
         to get back A in UMF format.
         Note also that UMF can operate directly on A-transpose.  This
         can save having to make another copy of the matrix, but it limited
         experiments, I found it to be slower. -DRN

         To form A-transpose in UMF format, merge the diagonal entries
         back into the rows.
      */
      k = 0;
      for (i = 0; i < ums->n; i++) { /* loop over rows */
        ums->atp[i] = k;
        hit_diag = FALSE;
        for (j = col[i]; j < col[i + 1]; j++) { /* loop over colums within row */
          /* if we get to the spot where the diagonal term belongs, merge it in */
          if (!hit_diag && col[j] > i) {
            ums->ati[k] = i;
            ums->atx[k] = a[i];
            k++;
            hit_diag = TRUE;
          }
          ums->ati[k] = col[j];
          ums->atx[k] = a[j];
          k++;
        }
        /* if we never got to the diagonal, merge it in now */
        if (!hit_diag) {
          ums->ati[k] = i;
          ums->atx[k] = a[i];
          k++;
          hit_diag = TRUE;
        }
      }
      ums->atp[ums->n] = ums->nnz;

      if (ums->nnz != k) {
        DPRINTF(stderr, "E: NNZ=%12d CT=%12d\n", ums->nnz, k);
        exit(0);
      }

      /* transpose matrix */
      err = umfpack_di_transpose(ums->n, ums->n, ums->atp, ums->ati, ums->atx, (int *)NULL,
                                 (int *)NULL, ums->ap, ums->ai, ums->ax);
      if (err != UMFPACK_OK) {
        fprintf(stderr, "UMFPACK error = %d\n", err);
        GOMA_EH(GOMA_ERROR, "Error computing matrix transpose using umfpack_di_transpose\n");
      }

      break;
    case 2: /* CSR FORMAT - NOT DONE YET */
      GOMA_EH(GOMA_ERROR, "Sorry, cannot convert CSR systems");
      break;
    }

    /* SET OPTIONS */
    switch (*fact_optn) {
    case -2: /* FULL ANALYSIS AND FACTORIZATION */
      umf_option = 1;
      break;
    case -1: /* FACTORIZATION WITH PAST ANALYSIS */
      umf_option = 0;
      break;
    case 0: /* FULL ANALYSIS AND FACTORIZATION */
      umf_option = 1;
      break;
    case 1: /* FACTORIZATION WITH PAST ANALYSIS */
      umf_option = 0;
      break;
    case 3:
      umf_option = 0;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Bad *fact_optn");
    }

    /* load default control parameters for UMF */
    umfpack_di_defaults(Control);
    /* optionally can ask for feedback from routines by uncommenting below */
    /*Control[UMFPACK_PRL] = 2.;*/
    /* optionally force solution strategy */
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;

    if (umf_option == 1) {
      /* analysis */
      if (ums->symbolic != NULL) {
        umfpack_di_free_symbolic(&ums->symbolic);
        ums->symbolic = NULL;
      }
      err = umfpack_di_symbolic(ums->n, ums->n, ums->ap, ums->ai, ums->ax, &ums->symbolic, Control,
                                Info);
      umfpack_di_report_status(Control, err);
      umfpack_di_report_info(Control, Info);
    }

    /* factorization */
    if (ums->numeric != NULL) {
      umfpack_di_free_numeric(&ums->numeric);
      ums->numeric = NULL;
    }
    err =
        umfpack_di_numeric(ums->ap, ums->ai, ums->ax, ums->symbolic, &ums->numeric, Control, Info);
    umfpack_di_report_status(Control, err);
    umfpack_di_report_info(Control, Info);
  }

  /* solve */
  if (*fact_optn >= 0) {
    err = umfpack_di_solve(UMFPACK_A, ums->ap, ums->ai, ums->ax, x, b, ums->numeric, Control, Info);
    umfpack_di_report_status(Control, err);
    umfpack_di_report_info(Control, Info);
  }

  return system_id;

} /* END of routine SL_UMF */
/*****************************************************************************/
/* END of file sl_umf.c */
/*****************************************************************************/

#else

int SL_UMF(int system_id,
           int *first,
           int *fact_optn,
           int *matr_form,
           int *nj,
           int *nnz_j,
           int *row,
           int *col,
           double *a,
           double *b,
           double *x) {
  GOMA_EH(GOMA_ERROR, "Goma not compiled with UMFPACK support");
  return -1;
}

#endif
