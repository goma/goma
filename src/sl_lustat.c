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

/*   provide a variety of interesting statistics and output 
 *   about the matrix that sparse gets...
 */

/*
 * $Id: sl_lustat.c,v 5.1 2007-09-18 18:53:48 prschun Exp $
 */

/*
 * $Log: not supported by cvs2svn $
 * Revision 5.0  2006/07/06 16:18:57  edwilke
 *
 * New Goma version: 'Farewell CRMPR, hello $3 a gallon gas!'
 *
 * Revision 4.1  2003/11/25 23:16:03  dalabre
 * The copyright statement has been updated for Goma and the version ratcheted
 * to 4.4.0.
 *
 * The makefile Goma.mk has been updated for Linux and Sun so that these
 * versions can be built easily (until the configure-make production is
 * complete); the Linux version is default.
 *
 * Added a prototype for function assemble_interface_extension_velocity_sic
 * in mm_fill_terms.h.
 *
 * Revision 4.0  2001/12/21 06:01:58  dalabre
 * Up Goma source code repository to V4.0. This identification 'coincides'
 * with our documentation upgrade. It will be tagged "Tora_Bora" in
 * recognition of the deep, unknown recesses that still remain in Goma.
 *
 * Revision 3.3  1999/12/22 16:25:18  pasacki
 * Last stages of cleanup. Now VP 6 works OK. GOMA now compiles cleanly
 * using gcc -Wall -O with AZTEC_2, HAVE_FRONT, PARALLEL, USE_CHEMKIN
 * on sparc-sun-solaris2.6.
 *
 * Revision 3.2  1999/12/20 17:36:43  pasacki
 * Kode Kleen nearing completion. Removed unused variables, redundant
 * declarations, initialized variables.
 *
 * Revision 3.1  1999/11/16 15:50:30  hkmoffa
 * Chemkin-Goma merger checkin.
 * This passed the FULL test suite last night, except for the single problem
 * that Duane pointed out had NaN's on. It also passed my test suite which
 * actually checks the goma answer against a blessed solution file.
 * That test suite has two mp problems in it, as well as 2 chemkin
 * problems.
 *
 * Revision 1.1.1.1  1999/03/11 22:01:20  hkmoffa
 * Import of Goma code from engsci 031199
 *
 * Revision 3.0  1999/01/07 06:14:49  dalabre
 * The IMPEACHMENT Version.
 *
 * Revision 2.4  1998/11/10 16:28:14  pasacki
 * Cleanup - unused variables eliminated; more prototypes.
 *
 * Revision 2.3  1998/03/05 21:57:32  pasacki
 * o Fixed up bugs in distributed processing modules.
 * o Minor cosmetic fix of ex_opts in mm_post_proc.c eliminates warnings of
 *   nonexistent element order map.
 * o New logic in mm_bc.c diverts duplicate BC listing to a file if
 *   threshhold criteria are exceeded (too many nodesets, sidesets or BCs).
 *
 * Revision 2.2  1997/09/10 19:16:55  dalabre
 * This check-in corrects some problems with umfpack from my last check-in and completes the prototyping for the sl_*.* files.
 *
 * Revision 2.1  1997/03/12 17:06:36  pasacki
 * Bugfixes, mainly to cure memory cancer arising from sloppy Aztec initialization.
 * Also, stepping forward in the direction of function prototypes. Eliminated
 * need to explicitly "-Dsolaris" and "-Daix".
 *
 * Revision 2.0  1996/09/30 22:03:13  prschun
 * FY97-REVISION 2.0
 *
 * Revision 1.2  1996/09/27 22:35:45  prschun
 * Hold on, I must commit all of this files to get rid of salsa.h
 *
 * Revision 1.1  1995/04/25  16:24:23  pasacki
 * o "Parts is parts." -- These are some new ones.
 *
 */






#ifdef MATRIX_STATISTICS


FILE *mfp;			/* for matrix file statistics; each iter */
FILE *mhp;			/* for matrix histogram plots; 2nd iter */
FILE *mpp;			/* for matrix profile plots; 2nd iter */

static int mhp_open = FALSE;
static int mpp_open = FALSE;

static int mpp_written = FALSE;

static void plot_a
  ( int,                     /* n  */
	   int,                     /* nnz  */
	   double [],               /* a[]  */
	   int []  );              /* ija[] */

static void histogram
  ( int,                     /* n  */
	   int,                     /* nnz  */
	   double [],               /* a[]  */
	   int [],                  /* ija[] */
	   int,                     /* lo  */
	   int,                     /* hi  */
	   int *  );               /* *d   */

/* static int call=0; */

void
lustat ( int n,
         int nnz,
         double a[],
         int ija[],
         double x[],
         spREAL norm,
         char *matrix  )
{
/* LOCAL VARIABLES */
  /*
  int error;
  int i;
  int ne;
  int nf;
  int first_time = TRUE;
  int decades;
  int lo, hi;			 For decade histogram... 
  int distribution[100];
  */

/*  char  *fsf; */
/*  char  *label; */

  spREAL cn;

  int    det_exp;
  double det_man;


  call++;

  fsf   = "lu1";
  label = "from goma";

  if ( ! mfp_open )
    {
      mfp = fopen("lu2","w");
      mfp_open = TRUE;
    }

  ne = spElementCount(matrix); 
  nf = spFillinCount(matrix);
  cn = spCondition(matrix, norm, &error); 

  spDeterminant(matrix, &det_exp, &det_man);
  
#ifndef _AIX
  spFileStats(matrix, fsf, label);
#endif

  lo = -25;
  hi =  8;
  decades = hi-lo+2;

  histogram(n, nnz, a, ija, lo, hi, distribution);

  fprintf(mfp, "# Statistics about direct matrix solution\n");
  fprintf(mfp, "\tOrder of system, n = %d\n", n);
  fprintf(mfp, "\tNum elements, spmatrix (factored) = %d\n", ne);
  fprintf(mfp, "\tNum nonzeroes, a (nominal) = %d\n", nnz);
  fprintf(mfp, "\tNum nonzeroes, a (actual) = %d\n", 
	  nnz-distribution[0]);
  fprintf(mfp, "\tFillin count = %d\n", nf);      
  fprintf(mfp, "\tNorm of matrix = %g\n", norm);       
  fprintf(mfp, "\tCondition number estimate = %g\n", cn);      
  fprintf(mfp, "\tDeterminant of matrix = %g x 10^%d\n", 
	  det_man, det_exp);

  plot_a(n, nnz, a, ija);

} /* END of routine lustat */
/*****************************************************************************/

static void
plot_a ( int n,
         int nnz,
         double a[],
         int ija[]  )
{
  int rt, ct, et;
  int ij_is_zero;
  int ji_is_zero;
  int row, col;
  int r, e;
  int sym;
  int i;

  if ( ! mpp_open )
    {
      mpp = fopen("ij.xy", "w");
      mpp_open = TRUE;
    }

  if ( ! mpp_written )
    {
      /*
       * Diagonal entries...
       */

      for ( i=0; i<n; i++)
	{
	  row = i+1;
	  col = i+1;
	  if ( a[i] != 0 )
	    {
	      fprintf(mpp, "1e32 150.\n%d. %d.\n", col, -row);
	    }
	  if ( a[i] == 0 )
	    {
	      fprintf(mpp, "1e32 1190.\n%d. %d.\n", col, -row);
	    }
      	}

      /*
       * Off-diagonal entries...
       */

      for (r=0; r<n; r++)
	{
	  for ( e=ija[r]; e<ija[r+1]; e++ )
	    {
	      row = r+1;
	      col = ija[e]+1;

	      ij_is_zero = TRUE;

	      if ( a[e] != 0 )
		{
		  ij_is_zero = FALSE;
		}
	      
	      /*
	       * Find out what a_ji is doing ...
	       */

	      ji_is_zero = TRUE;

	      rt = ija[e];	/* row of transpose */
	      ct = r;		/* col of transpose */
	      et = in_list(ct, ija[rt], ija[rt+1], ija);

	      if ( et != -1 )
		{
		  if ( a[et] != 0 )
		    {
		      ji_is_zero = FALSE;
		    }
		}
	      
	      /*
	       * Color asymmetric sparseness blue...
	       */
	      if ( ! ij_is_zero && ji_is_zero )
		{
		  fprintf(mpp, "1e32 21150.\n%d. %d.\n", col, -row);
		}
	      /*
	       * Color symmetric sparseness green...
	       */
	      if ( ! ij_is_zero && ! ji_is_zero )
		{
		  fprintf(mpp, "1e32 11150.\n%d. %d.\n", col, -row);
		}
	      
	    }
	}

      /*
       * Now, write out the inverse dof map info so that we can find out
       * who is what...do all 4 sides of the matrix...
       */
      for ( i=0; i<n; i++)
	{
	  if ( idv[pg->imtrx][i][0] == PRESSURE )
	    {
	      sym = 26101;
	    }
	  else
	    {
	      sym = 76101;
	    }
	  fprintf(mpp, "1e32 %d.\n%d. %d.\n%d. %d.\n", 
		  sym, -15, -i, -5, -i); /* left */
	  fprintf(mpp, "1e32 %d.\n%d. %d.\n%d. %d.\n", 
		  sym, n+5, -i, n+15, -i); /* right */
	  fprintf(mpp, "1e32 %d.\n%d. %d.\n%d. %d.\n", 
		  sym, i, 5, i, 15); /* top */
	  fprintf(mpp, "1e32 %d.\n%d. %d.\n%d. %d.\n", 
		  sym, i, -(n+5), i, -(n+15)); /* bot */
	}

      mpp_written = TRUE;
    }
  fclose(mpp);

} /* END of routine plot_a */
/*****************************************************************************/

static void
histogram ( int n,
            int nnz,
            double a[],
            int ija[],
            int lo,
            int hi,
            int *d )
{
  int i;
  int l;
  int index;
  double val;


  if ( ! mhp_open )
    {
      mhp = fopen("mh.d", "w");
      mhp_open = TRUE;
    }

  l = hi-lo+2;

  /*
   * Initialize before counting elements...
   */

  for ( i=0; i<l; i++)
    {
      d[i] = 0;
    }

  for ( i=0; i<nnz; i++)
    {
      if ( a[i] == 0 )
	{
	  d[0]++;
	}
      else
	{
	  val = log10(ABS(a[i]));
	  index = (int)(val-(dbl)lo) + 1;
	  if ( index < 1 )
	    {
	      index = 1;
	    }
	  if ( index > l )
	    {
	      index = hi;
	    }
	  d[index]++;
	}
    }

  /*
   * Histogram: [0] -- has number of elements that are zero
   *	    [1] -- has number of elements that are between
   *		   1 x 10**lo and 10 x 10**lo
   */
  
  for ( i=0; i<l; i++)
    {
      if ( i == 0 )
	{
	  fprintf(mhp, "%d %d\n", lo-2, d[0]);
	}
      else
	{
	  fprintf(mhp, "%d %d\n", lo+i-1, d[i]);
	}
    }

  /* fclose(mhp);*/

} /* END of routine histogram */
#endif

/******************************************************************************/
/* END of file sl_lustat.c */
/******************************************************************************/
