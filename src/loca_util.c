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
#ifndef lint
static char *cvs_util_id =
  "$Id: loca_util.c,v 5.1 2007-09-18 18:53:42 prschun Exp $";
#endif
*/
/*
 -----------------------------------------------------------------------------
   LOCA 1.0: Library of Continuation Algorithms
   Copyright (C) 2001, Sandia National Laboratories
 -----------------------------------------------------------------------------
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "loca_const.h"
#include "loca_util_const.h"
#include "mm_eh.h"

int N_o = -1;    /* Length of vector that is acted on (owned unknowns)  */
int N_t = -1;    /* Length of vector that is allcoated (total unknowns) */
double N_g = -1; /* Total number of unknowns on all procs, cast to a
                    double, which is the sum of N_o's */

/********** R O U T I N E S   I N   T H I S   F I L E  ***********************

       NAME                     TYPE                    CALL BY
----------------------------------------------------------------------
        intiialize_util_routines        void
        vec_init                        void
        vec_copy                        void
        alloc_vec                       double *
        dp                              double
        ip                              double
******************************************************************************/

/******************* PROTOTYPES FOR STATIC FUNCTIONS *************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void initialize_util_routines(int n_o, int n_t)
/*
 *  This routine sets up information that is used by the rouitnes in this file
 * and must be called before any of the routines within it are called.
 * The current input are "int n_o" the number of owned unknowns for this
 * processors, which is the length of a vectorr that is acted upton, and
 * "int n_t" which is the length that the vector must be allocated
 */
{
  N_o = n_o;
  N_t = n_t;
  N_g = gsum_double_conwrap((double)N_o);
} /* END of routine initialize_util_routines *********************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void vec_init(double *u)

/* Initialize the vector to zero.  */

{
  register int i;

  for (i = 0; i < N_t; i++)
    u[i] = 0.0;

} /* END of routine vec_init *************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double *alloc_vec(void)

/* Allocate a vector of length N_t */

{
  double *x;
  x = (double *)malloc(N_t * sizeof(double));
  vec_init(x);
  return x;

} /* END of routine alloc_vec *******************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void free_vec(double **ptr) {
  /*
   *  This version of free calls the system's free function
   *  with maximum error checking. It also doesn't call free if ptr is
   *  the NULL pointer.
   */

  if (*ptr != NULL) {

    free((void *)*ptr);

    /*
     *  Set the value of ptr to NULL, so that further references
     *  to it will be flagged.
     */

    *ptr = NULL;
  }
} /* free_vec */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void vec_copy(double *dx, double *dy)

/*
 * This function copies the value of the vector of type double, dx, into the
 * vector of type double, dy.  No error checking is done. If N is negative or
 * zero, dy is not changed.  A stride of 1 is assumed, unlike the linpack
 * routine called DCOPY (which is why the name is dcopy1).
 */

{
  register int i;
  for (i = 0; i < N_t; i++)
    dy[i] = dx[i];

} /* END of routine vec_copy *************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double dp(double *x, double *y)
/* simple dot product */
{
  int i;
  double sum = 0.0;

  for (i = 0; i < N_o; i++)
    sum += x[i] * y[i];
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double ip(double *x, double *y)
/* inner product for pitchfork tracking. This inner product must
   have the same symmetry as the pitchfork, so a simple dot product
   only works if the mesh is symmetric. Otherwise, this product
   must be weighted by the lumped mass at each node */
{
  return dp(x, y);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double ltransnorm(double *x, double *scale_vec) {
  /* rescale by factor to make averages array element size one */
  return (dp(x, scale_vec) / N_g);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double null_vector_resid(double r_val, double i_val, double *r_vec, double *i_vec, int mm_flag) {
  /* Calculate bifurcation equation residual norm
     (as Rayleigh quotient of null vector). */
  double *vr, *vi, *v1;
  double *Jy, *Jz, *By, *Bz;
  double yJy, yJz, zJy, zJz, yBy, yBz, zBy, zBz;
  double num_r, num_i, denom, norm_r, norm_i, norm;

  /* Proceed according to bifurcation method */
  if (i_val == 0.0) {

    /*
     * For these cases, the null vector is real:
     *   When a mass matrix is available:
     *       R.Q. = (n.J.n) / (n.B.n)
     *
     *   When a mass matrix is NOT available:
     *       R.Q. = (n.J.n) / (n.n)
     */

    /* r_vec needs to be copied to vector with space for externals */
    vr = alloc_vec();
    vec_copy(r_vec, vr);

    v1 = alloc_vec();
    matvec_mult_conwrap(vr, v1);
    num_r = dp(vr, v1);

    /* v1 is either n or B.n, depending on mm_flag */

    if (mm_flag)
      mass_matvec_mult_conwrap(vr, v1);
    else
      vec_copy(vr, v1);

    denom = dp(vr, v1);
    if (fabs(denom) > 1.0e-20) {
      norm_r = r_val - num_r / denom;
      norm = fabs(norm_r);
    } else {
      norm = -1.0;
    }
    free_vec(&vr);
    free_vec(&v1);
  } else {
    /*
     * At a Hopf bifurcation, the eigenvalue is imaginary and the null vector is
     * complex, so the following formula is used for the complex Rayleigh quotient:
     *
     *                           h          h
     *       R.Q. = (omega)i - (n .J.n) / (n .B.n)
     *
     *        h
     * where n  is the conjugate transpose of null vector n, J is the LSA
     * Jacobian matrix, and B is the mass matrix.
     *
     * Note: In this routine, y and z refer to the real and imaginary parts
     * of n, respectively.
     */
    /* r_vec needs to be copied to vector with space for externals */
    vr = alloc_vec();
    vec_copy(r_vec, vr);
    vi = alloc_vec();
    vec_copy(i_vec, vi);

    /* Allocate work vectors */
    Jy = alloc_vec();
    Jz = alloc_vec();
    By = alloc_vec();
    Bz = alloc_vec();

    /* Use work vectors for matrix-vector products */
    matvec_mult_conwrap(vr, Jy);
    matvec_mult_conwrap(vi, Jz);
    mass_matvec_mult_conwrap(vr, By);
    mass_matvec_mult_conwrap(vi, Bz);
    /* DEBUG: */

    /* Intermediate terms are then obtained by dot products */
    yJy = dp(vr, Jy);
    yJz = dp(vr, Jz);
    zJy = dp(vi, Jy);
    zJz = dp(vi, Jz);
    yBy = dp(vr, By);
    yBz = dp(vr, Bz);
    zBy = dp(vi, By);
    zBz = dp(vi, Bz);

    /* First, get denominator and check is it is too small to divide by */
    denom = (yBy + zBz) * (yBy + zBz) + (yBz - zBy) * (yBz - zBy);
    if (fabs(denom) == 0.0) {
      norm = -1.0;
    } else {

      /* Get real and imaginary parts of numerator */
      num_r = (yJy + zJz) * (yBy + zBz) + (yJz - zJy) * (yBz - zBy);
      num_i = (yJz - zJy) * (yBy + zBz) - (yJy + zJz) * (yBz - zBy);

      /* Construct the complex norm: (norm_r) + (norm_i)i */
      norm_r = r_val - num_r / denom;
      norm_i = i_val - num_i / denom;

      norm = sqrt(norm_r * norm_r + norm_i * norm_i);
    }

    /* Clean up and return */
    free_vec(&vr);
    free_vec(&vi);
    free_vec(&Jy);
    free_vec(&Jz);
    free_vec(&By);
    free_vec(&Bz);
  }

  /*
   * Flag value of -1.0 means denominator of Rayleigh Quotient is zero,
   * which is a zero mass matrix or zero eigenvectors
   */

  return norm;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void sort_by_real(int nconv, int ncv, int ldv, double *d, double *v) {
  /*
   * This routine takes the converged eigenvalue array d and assigns
   * index array j_index for sorting from highest to lowest real part.
   * Complex eigenvalue pairs will be identified with equal index values.
   */
  int i, j, k, l, m, skip, skip_next;
  int ntd, ntv;
  int *count;
  double *td, *tv;

  /* First check if nconv is smaller than two. */
  if (nconv < 2)
    return;

  /* Array count[i] will record how many eigenvalues are more than the i'th */
  ntd = 2 * ncv;
  ntv = (1 + nconv) * ldv;
  count = (int *)calloc(ncv + 1, sizeof(int));
  for (i = 0; i < ncv; i++)
    count[i] = 0;

  /* Mark final element of count array with -1 */
  count[ncv] = -1;
  // remove faulty warning GCC 12
  GOMA_ASSERT_ALWAYS(((size_t)(ntd * sizeof(double))) < PTRDIFF_MAX);
  /* Arrays td and tv will hold the sorted eigenvalues and eigenvectors */
  td = (double *)calloc(ntd, sizeof(double));
  for (i = 0; i < ntd; i++)
    td[i] = 0.0;
  tv = (double *)calloc(ntv, sizeof(double));
  for (i = 0; i < ntv; i++)
    tv[i] = 0.0;

  /* Compare real parts for all i,j eigenvalue pairs, where j has the larger
   * initial index. Increment count for the smaller real part.
   * Here, when eigenvalue j is the first of a complex eigenpair, "skip_next"
   * will be TRUE; if it is the second of the pair, "skip" will be TRUE. */
  skip = FALSE;
  skip_next = FALSE;
  for (i = 0; i < nconv; i++) {

    /* Determine if this is the first of a complex eigenpair */
    /* If this is the second of a complex eigenpair, reset the skip_next flag */
    if (skip)
      skip_next = FALSE;
    else
      skip_next = ((d[i + ncv] == 0.0) ? FALSE : TRUE);
    for (j = 0; j < i; j++) {

      /* Do not compare complex conjugates - this ensures
       * that both will have the same value of count */
      if (!skip || i != j - 1) {

        /* If d values are different, increment count for the smaller one. */
        if (d[j] < d[i]) {
          count[j]++;
        } else if (d[j] > d[i]) {
          count[i]++;
        }

        /* If d values are the same but not a complex eigenpair,
         * increment count for the larger index (always j) */
        else {
          count[i]++;
        }
      }
    } /* End of inner loop */

    /* set skip for next pass of j */
    skip = skip_next;
  } /* End of outer loop */

  /* Now copy eigenvalues and eigenvectors into temporary arrays td and tv
   * in their sorted positions as determined by the count array */
  skip = FALSE;
  skip_next = FALSE;
  for (j = 0; j < nconv; j++) {

    /* Initial position: j --> Sorted position: i = count[j] */
    i = count[j];

    /* Determine if this is the second of a complex eigenpair;
     * if so, copying was done on previous pass */
    if (!skip) {

      /* Determine if this is the first of a complex eigenpair, set skip_next */
      skip_next = ((count[j + 1] == i) ? TRUE : FALSE);

      /* Copy eigenvalue into td */
      td[i] = d[j];
      td[i + ncv] = d[j + ncv];

      /* Copy complex conjugate eigenvalue into next td position if applicable */
      if (skip_next) {
        td[i + 1] = d[j];
        td[i + ncv + 1] = d[j + ncv + 1];
      }

      /* Copy eigenvector into tv */
      for (k = 0; k < ldv; k++) {

        /* Assign indices into v and tv */
        l = i * ldv + k;
        m = j * ldv + k;
        tv[l] = v[m];

        /* Copy corresponding element of complex conjugate eigenvector
         * into next tv position if applicable */
        if (skip_next) {
          tv[l + ldv] = v[m + ldv];
        }
      } /* End of eigenvector loop */
    }   /* End of if (!skip) */

    /* If this is the second of a complex eigenpair, just reset skip_next */
    else {
      skip_next = FALSE;
    }

    /* set skip_next for next pass */
    skip = skip_next;
  } /* End of eigenvalue loop */

  /* Now rewrite d and v arrays in sorted order */
  for (i = 0; i < nconv; i++) {
    d[i] = td[i];
    d[i + ncv] = td[i + ncv];
    for (j = 0; j < ldv; j++) {
      k = i * ldv + j;
      v[k] = tv[k];
    }
  }

  /* Free temporary arrays. */
  free((void *)count);
  free((void *)td);
  free((void *)tv);

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
