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
 * Created:  Wed Nov 23 08:58 MST 1994 pasacki@sandia.gov
 */

/*
 * Default: do not attempt to use Harwell MA28 linear solver. Kundert's is
 *          more robust and Harwell has a better successor to MA28 that you
 *          can buy with money.
 *
 * #define HARWELL
 */

#include "sl_util.h"
#ifdef USE_RCSID
static char rcsid[] = "$Id: sl_ma28.c,v 5.1 2007-09-18 18:53:48 prschun Exp $";
#endif




#ifdef HARWELL
static int call = 0;

static int n_previous;		/* To assist in error checking; check for */
				/* consistency in the order of the system */
				/* and the number of nonzero entries */
				/* between different solves.*/
static int nnz_previous;
#endif

/*
 * These definitions will depend on the C/Fortran calling convention
 * for the system. 
 *
 * These should be keyed off a centralized determination and not hidden
 * in the intestines like this.
 */

#ifdef _AIX
#define MA28AD	ma28ad
#define MA28BD  ma28bd
#define MA28CD	ma28cd

#define MANGLE	mangle
#endif

#ifdef solaris
#define MA28AD	ma28ad_
#define MA28BD  ma28bd_
#define MA28CD	ma28cd_

#define MANGLE	mangle_
#endif

#ifdef sun
#ifndef solaris
#define MA28AD	ma28ad_
#define MA28BD  ma28bd_
#define MA28CD	ma28cd_

#define MANGLE	mangle_
#endif
#endif

#ifdef hpux
#define MA28AD	ma28ad
#define MA28BD  ma28bd
#define MA28CD	ma28cd

#define MANGLE	mangle
#endif

/*
 * return values:  0 -- if everything went OK
 *		  -1 -- if something went wrong
 */


int
cmsr_ma28 ( const int n,		/* order of matrix system */
            const int nnz,		/* nominal number of nonzeroes in a */
            double a[],		/* vector of nonzeroes in a matrix */
            int ija[],		/* pointers for finding position of entries */
            double x[],		/* space for solution vector */
            double b[]  )	/* right hand side vector */

{
/* LOCAL VARIABLES */
  int status;

#ifdef HARWELL
  int err;
  int iaction;

  static int licn;		/* length of A[] matrix */
  static int *irn;		/* holds row indeces on input */
  static int lirn;		/* length of irn[] array */
  static int *icn;		/* holds column indeces on input */
  static int *ikeep;		/* workspace; preserve between calls */
  static int *iw;

  static int *ivect, *jvect;	/* used for 2nd and subsequent factorizations */

  int analyse, factor, solve;

  int iflag;

  int index;
  int mtype;
  int begin, end;
  int e, r, i;


  dbl u;			/* control pivoting strategy to numeric */
				/* stability*/
  static dbl *w;		/* work array for MA28 */
  static dbl *A;		/* holds nonzero entries in Harwell's format */
#endif

  status = 0;


#ifdef HARWELL
  /*
   * Begin execution.
   */
  
  call++;

  analyse = ( call < 3 );
  factor  = TRUE;
  solve   = TRUE;

  u = 1.0;			/* 0 -- emphasize sparseness */
				/* 1 -- numerical stability */
  if ( call == 1 )
    {
      n_previous   = n;
      nnz_previous = nnz;

      /*
       * Need to access the true count for how much storage is needed for
       * the factored LU matrices...
       *
       * The first numbers are too low; the second figures are too high...
       */

      iaction = 3;
      MANGLE (&iaction);	/* set the tolerance nice and low */

      lirn = nnz;
      licn = nnz;

      licn = n*n;
      lirn = n*n;


      irn  = (int *)calloc(lirn, sizeof(int));
      if ( irn == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for irn\n",
		  lirn, sizeof(int));
	  exit(-1);
	}
      
      ivect  = (int *)calloc(nnz, sizeof(int));
      if ( ivect == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for ivect\n",
		  nnz, sizeof(int));
	  exit(-1);
	}
      
      icn = (int *)calloc(licn, sizeof(int));
      if ( icn == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for icn\n",
		  licn, sizeof(int));
	  exit(-1);
	}
      
      jvect = (int *)calloc(nnz, sizeof(int));
      if ( jvect == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for jvect\n",
		  nnz, sizeof(int));
	  exit(-1);
	}
      
      A = (dbl *)calloc(licn, sizeof(dbl));
      if ( A == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for A\n",
		  licn, sizeof(dbl));
	  exit(-1);
	}
      
      ikeep = (int *)calloc(5*n, sizeof(int));
      if ( ikeep == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for ikeep\n",
		  5*n, sizeof(int));
	  exit(-1);
	}
      
      iw = (int *)calloc(8*n, sizeof(int));
      if ( iw == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for iw\n",
		  8*n, sizeof(int));
	  exit(-1);
	}
      
      w = (dbl *)calloc(n, sizeof(dbl));
      if ( iw == NULL )
	{
	  fprintf(stderr, 
		  "Problem allocating %d elements of size %d for w\n",
		  n, sizeof(dbl));
	  exit(-1);
	}
      
      iflag = 0;
    }

  if ( n != n_previous )
    {
      fprintf(stderr, "sl_ma28: n: %d != %d\n", n, n_previous);
      return(-1);
    }

  if ( nnz != nnz_previous )
    {
      fprintf(stderr, "sl_ma28: nnz: %d != %d\n", nnz, nnz_previous);
      return(-1);
    }

  /*
   * Fill new A[] matrix with entries stored in a[] and ija[]...
   */

  /*
   * Diagonal entries...
   */

  index = 0;

  for ( i=0; i<n; i++ )
    {
      A[i]   = a[i];
      if ( analyse )
	{
	  irn[i] = i+1;
	  icn[i] = i+1;
	}
      else
	{
	  ivect[index] = i+1;
	  jvect[index] = i+1;
	}      
      index++;
    }

  /*
   * Row by row pick out the off-diagonal nonzero entries...
   */

  for (r=0; r<n; r++)
    {
      for ( e=ija[r]; e<ija[r+1]; e++ )
	{
	  if ( analyse )
	    {
	      icn[index] = ija[e]+1;
	      irn[index] = r+1;
	    }
	  else
	    {
	      ivect[index] = r+1;
	      jvect[index] = ija[e]+1;
	    }
	  A[index]   = a[e];
	  index++;
	}
    }	  
  
  if ( analyse )
    {
      /*
       * Analyze and factor...
       */

      MA28AD (&n, &nnz, A, &licn, irn, &lirn, icn, &u, ikeep, iw, w, &iflag);

      if ( iflag != 0 )
	{
	  fprintf(stderr, "sl_ma28: MA28AB: iflag = %d\n", iflag);
	  status = -1;
	  iaction = 0;
	  MANGLE (&iaction);
	}
    }
  else if ( factor )
    {
      /*
       * Factor matrix with same sparsity pattern but different
       * numerical values for the matrix entries...
       */

      MA28BD (&n, &nnz, A, &licn, ivect, jvect, icn, ikeep, iw, w, &iflag);

      if ( iflag != 0 )
	{
	  fprintf(stderr, "sl_ma28: MA28BD: iflag = %d\n", iflag);
	  status = -1;
	  iaction = 0;
	  MANGLE (&iaction);
	}
    }

  /*
   * Solve step; back-substitution requires loading the RHS vector...
   */

  mtype = 1;			/* Solve direct equation; not transpose */

  if ( solve )
    {
      MA28CD (&n, A, &licn, icn, ikeep, b, w, &mtype);
      for ( r=0; r<n; r++)
	{
	  x[r] = b[r];
	}
    }

#endif
  return(status);

} /* END of routine cmsr_ma28 */
/******************************************************************************/
/* END of file sl_ma28.c */
/******************************************************************************/

