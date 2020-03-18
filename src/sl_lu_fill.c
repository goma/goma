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
 * $Id: sl_lu_fill.c,v 5.1 2007-09-18 18:53:48 prschun Exp $
 */

#include <stdlib.h>

#include "spMatrix.h"
#include "std.h"
#include "sl_util.h"

/*
 * HKM -> I wonder if first_time_fill needs to have global scope?
 *        It is only used in the function below.
 *        Shouldn't it be declared a static variable inside luf()?
 *        -> will leave it as is right now.
 */
int first_time_fill = TRUE;
#ifdef HAVE_SPARSE

void
luf (const int N,
     const int NExt,
     const int  M,
     double  a[],
     int  ija[],
     double x[],
     const int factor_flag )

/*
  	
  	C driver for direct lu factor of sparse matrix.
  	
  		This routine takes a C-MSR format
  		sparse matrix and determines the lu
  		decomposition. The package sparse1.3
  		is used to factor and solve the matrix.
  		This routine is somewhat inefficient in its
  		use of memory.
  		
  		J.N. Shadid SNL Div 1421
  	
  	   date: 2/27/93		
  	
  
         Parameter list:
  
 
       N    ==         Number of unknowns updated by this processor
       NExt ==         Number of unknowns updated by neighboring procs
                       for which copies are kept on this processor
       M    ==         upper bound on the number of nonzeros in matrix A
                        in sparse format (C-MSR).
       a[]   ==        matrix A in sparse format (C-MSR).
       ija[] ==        pointers to nonzeros of A (C-MSR).

       x   ==          On input 'x' contains the right hand side. On output
                        'x' contains the solution to our linear system.

  
                factor_flag: 1, 2, 3
  
                        = 1  Dynamically allocate memory for tmp vectors.
                             Calculate lu factors, solve system with LU.
  
                        = 2  Calculate lu factors, solve system with LU.
  
                        = 3  use previus lu factors to solve system with LU.
  
                x ==    on input is the residual(rhs) of the set of 
                        equations, on output is the result.
  
  
    --------------------------------------------------------------------
          routines called: (sparse1.3 package is used)
  		
    --------------------------------------------------------------------
  
*/
{
  static char *matrix;
  int  error, type;
  int j, i, n, k, nzeros, ija_col;
  static spREAL **element, *b;
#ifdef MATRIX_STATISTICS
  spREAL norm;
#endif
  static int call = 0;

  call++;

   /* allocate rhs and copy */

   b = (spREAL *)malloc((N+NExt+1)*sizeof(spREAL));
   for(i = 0; i < N+NExt; i++) b[i+1] = x[i];


   /* factor matrix if values have changed since last call */
   if( factor_flag <= 2 ){

     if( first_time_fill == TRUE && factor_flag == 1){
        /* N = order of real matrix */
        type = 0;
        matrix = spCreate(N+NExt,type,&error);

	element = (spREAL **)malloc((M+NExt+2)*sizeof(spREAL *));
     }

     /* n = 1 produce graph of matrix, n = 2 fill matrix */
     for(n = 1; n <= 2; n++){ /* define and fill */


	/* copy MSR matrix into linked list format */
        k = 1;
        for(i = 0; i <  N; i++){ /* lop through rows */
      	   nzeros  = ija[i+1] - ija[i];
	   ija_col = ija[i];

           if( n == 1 ){           /* create nonzero pattern */

	     /* diagonal element */
             element[k++] = spGetElement(matrix,i+1,i+1);

	     /* nonzero off diagonal elements */
	     for(j = 0; j < nzeros; j++){
                element[k++] = spGetElement(matrix,i+1,ija[ija_col++]+1);
	     }	
           }
           else {                  /* fill in numerical values */

	     /* diagonal element */
             spADD_REAL_ELEMENT(element[k++],a[i]);

	     /* nonzero off diagonal elements */
	     for(j = 0; j < nzeros; j++){
                spADD_REAL_ELEMENT(element[k++],a[ija_col++]);
	     }
           }
        } /* end row loop */

	/* fill in Dirchlet equations for the interface unknowns */
	for(i = N; i <  N+NExt; i++){
	   if( n == 1 ){           /* create nonzero pattern */
	     element[k++] = spGetElement(matrix,i+1,i+1);
           }
           else {                  /* fill in numerical values */
	     spADD_REAL_ELEMENT(element[k++],1.0);
	   }
        }

        if (n == 1) spClear(matrix);/* zero entries in graph */

     } /* end of define and fill */
#ifdef MATRIX_STATISTICS
     norm = spNorm(matrix);
#endif

      if( first_time_fill == TRUE && factor_flag == 1){

/*         spFileMatrix(matrix,"matrix_file","channel",0,1,1); */
         spOrderAndFactor(matrix,b,-1.0,0.0,1);
	 first_time_fill = FALSE;

       } else {

         spFactor(matrix); 

       }


  } /* alocate matrix and factor matrix */

  /* solve lu factored system */

  spSolve(matrix,b,x-1);

#ifdef MATRIX_STATISTICS
  if ( call == 4 )
    {
      lustat(N, M, a, ija, x, norm, matrix);
    }
#endif

/*
  for(i = 0; i < N; i++) printf("i: %d\t%f\t%f\n",i,x[i],b[i+1]);
*/
  free(b);
/*  spDestroy(matrix);*/

} /* END of routine luf */

#else //HAVE_SPARSE
void
luf (const int N,
     const int NExt,
     const int  M,
     double  a[],
     int  ija[],
     double x[],
     const int factor_flag )
{
  EH(-1, "Goma not configured with sparse solver support");
}
#endif // HAVE_SPARSE
/******************************************************************************/
/* END of file sl_lu_fill.c */
/******************************************************************************/


