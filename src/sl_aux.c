/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* 
 *  sl_aux.c -- Handy routines, LU decomposition and back sub 
 *
 *  Randy Lober 9/98.
 *
 */

#ifndef lint
#endif

#include <stdio.h>
#include <math.h>

#include "rf_allo.h"
#include "sl_aux.h"

#define TINY 1.0e-20

/* The lu_decomp_backsub_driver function manages the access to the 
   lu_decomp and lu_backsub routines for the solution of a system of linear 
   equations, Ax = b.

   Using the lu_dcmp_flag value (1 = perform lu_decomp & back_sub, 0 = just back_sub.),
   this driver can be called over and over to solve new rhs vectors to the same
   coeff matrix A. Note that when back_sub is only being employed, the input
   coeff_matrix must be the resulting lu_decomp form produced by a previous
   lu_decomp_backsub_driver call.

   If lu_dcmp_flag == 1: (perform LU decompositon and backsub)

   Inputs: coeff_matrix[0..sys_size-1][0..sys_size-1] (initial matrix A of the system)
           rhs_vector  [0..sys_size-1]                (rhs vector of the system)
	   indx        [0..sys_size-1] (empty working space vector)
	   sys_size (the dimension size of coeff_matrix, rhs_vector, & indx)
	   lu_dcmp_flag (flag to control lu_decomp production - 1 = perform
                         lu decomposition, 0 simply back_sub)

   Output: coeff_matrix[0..sys_size-1][0..sys_size-1] (LU decompositon of coeff_matrix)
           rhs_vector[0..n-1] (containing the x solution vector)

   If lu_dcmp_flag == 0: (Just backsub the input rhs - LU already exists in coeff_matrix)

   Inputs: coeff_matrix[0..sys_size-1][0..sys_size-1] (LU decompositon of coeff_matrix)
           rhs_vector  [0..sys_size-1]                (rhs vector of the system)
	   indx        [0..sys_size-1] (filled working space vector)
	   sys_size (the dimension size of coeff_matrix, rhs_vector, & indx)
	   lu_dcmp_flag (flag to control lu_decomp production - 1 = perform
                         lu decomposition, 0 simply back_sub)

   Output: rhs_vector[0..n-1] (containing the x solution vector)

   lu_decomp_backsub_driver returns 0 if successful completion, -1 otherwise.

   Randy R. Lober 9/98.
*/

int
lu_decomp_backsub_driver ( double **coeff_matrix, 
			   double *rhs_vector,
			   int *indx,
			   int sys_size,
			   int lu_dcmp_flag ) 
{
  static double *d=NULL;
  static int first_call=1;

  if( first_call ) 
     {
       d = (double *) smalloc (1 * sizeof(double));
       first_call = 0;
     }

  if ( lu_dcmp_flag == 1 ) {
    if ( lu_decomp ( coeff_matrix,
		     sys_size,
		     indx,
		     d ) == -1 ) {
      fprintf ( stdout,
		" Error occurred in lu_decomp_backsub_driver\n");

      return (-1);
    }
  }

  lu_backsub ( coeff_matrix,
	       sys_size,
	       indx,
	       rhs_vector );
  return (0);
}

/* The lu_decomp function contributes to the solution of a system of linear 
   equations, Ax = b by breaking the A matrix into two successive linear 
   systems, the lower & upper triangular matrixes (L & U) that solve the same 
   linear system such that Ax = LUx = b. Used together with lu_backsub to solve
   multiple linear systems.

   Taken directly from Numerical Recipes in C, Press, W. H.,
   Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P., 2nd Ed.,
   Cambridges University Press, ISBN 0 521 43108 5, pp. 46 - 47. 

   Minor changes include the usage of double instead of float, etc., and converting
   from 1 based arrays to 0 based (This is supposed to a book on algorithms
   in C...? )

   Inputs: a[0..n-1][0..n-1] (initial matrix A of the system)
           n (system dimension)
	   indx[0..n-1] (working vector)
	   d[1] (row exchange counter)

   Output: a[0..n-1][0..n-1] (LU decompositon of A)
           indx[0..n-1] (vector recording row permutation effected by partial pivoting
	   d (+-1 depending on whether number of row exchanges was even or odd)

   lu_decomp returns 0 if successful completion, -1 otherwise.

   Randy R. Lober 9/98.
*/

int
lu_decomp ( double **a, 
	    const int n, 
	    int *indx, 
	    double *d ) 
{
  int i, imax=-1, j, k;
  double big, dum, sum, temp;
  /* vv stores the implicit scaling of each row */
  static double *vv=NULL;  
  static int n_max=0;

  if (n > n_max) {
    safer_free((void **) &vv);
    vv = alloc_dbl_1(n, DBL_NOINIT); 
    n_max = n;    
  }

  *d = 1.0; /* No row interchanges yet */
  for ( i = 0; i < n; i++ ) { /* loop over rows to get the implicit scaling info */
    big = 0.0;
    for ( j = 0; j < n; j++ ) {
      if ( ( temp = fabs( a[i][j] ) ) > big) {
	big = temp; 
      }
    }

    if ( big == 0.0 ) {
      printf("Singular matrix in routine lu_decomp - aborting\n");
      return (-1);
    }
    vv[i] = 1.0/big; /* save the implicit scaling info */
  }

  for ( j = 0; j < n; j++ ) {    /* Looping over columns (Crout's method) */
    for ( i = 0; i < j; i++ ) {  /* This is Eq. 2.3.12 except for i == j */
      sum = a[i][j];
      for ( k = 0; k < i; k++ ) {
	sum -= a[i][k]*a[k][j];
      }
      a[i][j] = sum;
    }

    big = 0.0;                  /* Initialize the search for the largest pivot element */
    for ( i = j; i < n; i++ ) { /* This is i == j of Eq. 2.3.12 & i = j+1..N for Eq. 2.3.13 */
      sum = a[i][j];
      for ( k = 0; k < j; k++ ) {
	sum -= a[i][k]*a[k][j];
      }
      a[i][j] = sum;
      if( ( dum = vv[i]*fabs(sum) ) >= big ) { /* Is the figure of merit for the pivot */
	big = dum;                             /* better than the best so far? */
	imax = i;
      }
    }

    if( j != imax ) {              /* Do we need to interchange rows? */     
      for ( k = 0; k < n; k++ ) {  /* Interchange them */
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = dum;
      }
      (*d) = -(*d);             /* ...and change the parity of d */
      vv[imax] = vv[j];         /* also change the scale factor */
    }

    indx[j] = imax;
    if( a[j][j] == 0.0 ) {
      a[j][j] = TINY;
    }                         /* If the pivot element is 0, the matrix is singular
			         (at least to the precision of the algorithm). For
			         some applications on singular matrices, it is desirable 
			         to substitute TINY for 0. */

    if( j != (n - 1) ) {      /* Divide by the pivot element */
      dum = 1.0/(a[j][j]);
      for( i = j + 1; i < n; i++ ) {
	a[i][j] *= dum;
      }
    }
  }       /* Get the next column */
  return (0);
}  

/* The lu_backsub function contributes to the solution of a system of linear 
   equations, Ax = b by [erforming the forward and backward substitution
   segment of the sequence, generating the solution vector x which
   replaces the input vector for the RHS, b.  lu_backsub is used together
   with lu_decomp to solve multiple linear systems.

   Taken directly from Numerical Recipes in C, Press, W. H.,
   Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P., 2nd Ed.,
   Cambridges University Press, ISBN 0 521 43108 5, pp. 47 - 48. 

   Minor changes include the usage of double instead of float, etc., and converting
   the input arrays so the lame 1 based code can be used (This is supposed to a 
   book on algorithms in C...? )

   Inputs: a[0..n-1][0..n-1] (LU decompositon of A)
           n (system dimension)
           indx[0..n-1] (vector recording row permutation effected by partial pivoting
	   b[0..n-1] (RHS vector of the system).

   Output:
           b[0..n-1] (containing the x solution vector)

   Randy R. Lober 9/98.
*/

void 
lu_backsub ( double **a,
	     int n, 
	     int *indx, 
	     double *b )  
{
  int i, ii, ip, j;
  double sum;

  ii = -1;
  
  for ( i = 0; i < n; i++ ) {       /* When ii is set to a value >= 0, it will */
                                    /* become the index of the first nonvanishing
				       element of b. We now do the forward 
				       substitution, Eq. 2.3.6. The only new
				       wrinkle is to unscramble the permutation
				       as we go. */
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if ( ii != -1 ) {
      for( j = ii; j <= i - 1; j++ ) { 
	sum -= a[i][j]*b[j];
      }
    } else if( sum ) {            /* A nonzero element was encountered, so from
				     now on, we have to do the sums within the
				     above loop. */
      ii = i;
    }
    b[i] = sum;
  }

  for ( i = n - 1; i >= 0; i-- ) { /* Now do the backsubstitution (Eq. 2.3.7) */
    sum = b[i];
    for ( j = i + 1; j <= n - 1; j++ ) {
      sum -= a[i][j]*b[j];
    }
    b[i] = sum/a[i][i];            /* Store a component of the solution
				      vector x */
  }
}

/******************************************************************************/
/* END of file sl_aux.c */
/******************************************************************************/

