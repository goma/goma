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
 * sl_squash -- auxiliary routine squashes out zeroes in the sparse matrix...
 *		the a[] and ija[] vectors are trimmed down to size
 */


/* static int first_time; */

/*
 * Return values:
 *
 *		new estimate for number of real nonzeroes
 *		-1 if something went wrong.
 *
 *		the a[] matrix in MSR format will be smaller,
 *		as will the ija[] pointer matrix.
 *
 *		This will NOT get rid of rows; assume the order of
 *		system will remain the same.
 *
 *		If there are zeroes on the diagonal, then we assume
 *		there are some nonzeroes hanging around on the off-diagonals
 *		somewhere...
 *
 *		The intended use is to further pack a[] in case the
 *		initial allotment of space for nonzero entries was overly
 *		generous.
 */
#if 0
int
squash ( int N,
         int M, 
         double a[],
         int ija[],
         double tol  )
{
/* LOCAL VARIABLES */
  int     i,  j,  r,  e;
  int     *ija_new;
  double  *a_new;


  if ( tol < 0 )
    {
      EH(-1, "Zero threshhold must be >= 0.");
    }

  a_new   = (dbl *)calloc(M, sizeof(dbl));
  ija_new = (int *)calloc(M, sizeof(int));

  if ( a_new == NULL || ija_new == NULL )
    {
      EH(-1, "Problem allocating atmp.");
    }

  /*
   * Diagonal terms stay the same; even zeroes are saved.
   */

  for ( i=0; i<N; i++)
    {
      a_new[i] = a[i];
    }

  /*
   * Loop through each row of the matrix, examining each of the
   * purported nonzero entries whose column numbers are recorded
   * in the appropriate segment of ija[]. Values are in a[].
   */

  j = N+1;
  for ( r=0; r<N; r++)
    {
      ija_new[r] = j;
      for ( e=ija[r]; e<ija[r+1]; e++)
	{
	  /*
	   * For really nonzero entries...
	   */
	  if ( abs(a[e]) > tol )
	    {
	      a_new[j]   = a[e]; /* Copy the value. */
	      ija_new[j] = ija[e]; /* Copy the column index. */
	      j++;		/* Increment the counter for really nonzeroes.*/
	    }
	}
    }

  ija_new[N] = j;

  /*
   * Check to insure there are no all zero rows...
   */

  for ( r=0; r<N; r++)
    {
      if ( a_new[r] == 0 )
	{
	  if ( ija_new[r] == ija_new[r+1] )
	    {	  
	      fprintf(stderr, "row = %d\n", r );
	      EH(-1, "Zero diagonal has no nonzero off-diagonals");
	    }
	}
    }

  fprintf(stdout, "Compression: %d -> %d (%.2f %%)\n", M, j, 
	  1e2 * (1. - (dbl)j / (dbl)M)   );
  
  /*
   * Copy new packed vector into old space...
   */

  for ( i=0; i<j; i++)
    {
      a[i]   = a_new[i];
      ija[i] = ija_new[i];
    }

  free(a_new);

  free(ija_new);

  return(j);
} /* END of routine squash */
#endif

/*****************************************************************************/
/* END of file sl_squash.c */
/*****************************************************************************/
