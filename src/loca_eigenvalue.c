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
#ifndef lint
static char *cvs_solnonlin_id =
  "$Id: loca_eigenvalue.c,v 5.1 2007-09-18 18:53:42 prschun Exp $";
#endif
*/
/*
 -----------------------------------------------------------------------------
   LOCA 1.0: Library of Continuation Algorithms
   Copyright (C) 2001, Sandia National Laboratories
 -----------------------------------------------------------------------------
*/

/*
 * The following automates switches between ARPACK and PARPACK calls.
 * Only compiler directives EIGEN_PARALLEL or EIGEN_SERIAL apply.
 */

#ifdef EIGEN_PARALLEL
# include <mpi.h>
#else
# if !defined(MPI) && !defined(DLB)
#   define MPI_Comm int
#   define MPI_COMM_WORLD 0
# else
#include <mpi.h>
# endif
#endif

/* Function prototypes go here! */

extern void poleze_(int* m, double* sigma, double* mu, double* realpt,
                    double* imagpt, double* ritzes, double* trealpt,
                    double* timagpt, double* tritzes, int* tm, int* nconv,
                    double* tol, int* info, double* workreal, double* workimag,
                    double* writzes);
extern void polez2_(int* m, double* sigma, double* mu, double* delta,
                    double* realpt, double* imagpt, double* ritzes,
                    double* trealpt, double* timagpt, double* tritzes, int* tm,
                    int* nconv, double* tol, int* info);
extern void polez3_(int* m, double* sigma, double* mu, double* delta,
                    double* zeta,double* realpt, double* imagpt, double* ritzes,
                    double* trealpt, double* timagpt, double* tritzes, int* tm,
                    int* nconv, double* tol, int* info);
extern void stslct_(int* m, int* nconv, double* tol, double* sigma, double* mu,
                    double* cutoff, double* ritzr, double* ritzi, double* errbds,
                    int* select);
extern void stslc2_(int* m, int* nconv, double* tol, double* sigma, double* mu,
                    double* delta, double* cutoff, double* ritzr, double* ritzi,
                    double* errbds, int* select);
extern void stslc3_(int* m, int* nconv, double* tol, double* sigma, double* mu,
                    double* delta, double* zeta, double* cutoff, double* ritzr,
                    double* ritzi, double* errbds, int* select);
extern void cpdnaupc_(MPI_Comm* comm, int* ido, char* bmat, int* nloc,
                       char* which, int* nev, double* tol, double* resid,
                       int* ncv, double *v, int* ldv, int* iparam, int* ipntr,
                       double* workd, double* workl, int* lworkl, int* info);
extern void cpdneupc_(MPI_Comm* comm, int* ivec, char* howmny, int* celect,
                       double *d, double* v, int* ldv, double *sigma,
                       double *mu, double* delta, double* workev, char* bmat,
                       int* n, int* n2, char* which, int* nev, double* tol,
                       double* resid, int* ncv, int* iparam, int* ipntr,
                       double* workd, double* workl, int* lworkl, int* ierr,
                       int* select);
extern void cpdmout_  (MPI_Comm* comm, int* lout, int* m, int* n, double *A,
                      int* lda, int* idigit);
extern double second_();


#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>

#include "loca_const.h"
#include "loca_util_const.h"

#define SHIFTS   0
#define MAX_ITRS 2
#define MODE     6

#define ORD_SYM  0  /* meaningless for this driver */
#define GEN_SYM  1  /* meaningless for this driver */
#define SHIFTI   3  /* Shift and Invert */
#define CAYLEY   4  /* mode 4 in dnaupc, used for Cayley transform */

/*********** G L O B A L   F U N C T I O N S   I N   T H I S   F I L E ********
*
*     function                type                       called_by
*    -----------           -----------                 ---------------
*
*   calc_eigenvalues            void                 solve_nonlinear_problem()
*
******************************************************************************/
/*****************************************************************************/
/******* P R O T O T Y P T E  S   F O R   S T A T I C   F U N C T I O N S ****/
/*****************************************************************************/

#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
static int  eig_driver(char which[], char bmat[], int iparam[], int mode,
                       double sigma, double mu, double delta, double zeta,
                       int nev, int ncv, int info, double tol, double eta,
                       int printproc, int numOwnedUnks, int numUnks,
                       int con_step_num, int jmax, int sort, MPI_Comm comm);
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void calc_eigenvalues_loca(struct con_struct *con)

{
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)

  /* shorthand for long con sub-structure */
  struct general_info_struct *cgi = &(con->general_info);
  struct eigen_info_struct   *cei = &(con->eigen_info);

  double *x = cgi->x, *rhs;

  /* Eigenvalue definitions */
  int      iparam[11];
  int      nev, ncv, info, mode, az_fail_cnt;
  int      jmax;
  double   tol, sigma, eta, mu, delta, zeta;
  char     bmat[2], which[3];
  MPI_Comm comm;           /* MPI communicator                   */

/* --------- Execution begins  ---------------- */

  if (cgi->printproc > 1) {
    printf("\n");
    printf("\tStarting Eigenvalue Calculation");
    printf(" (2 matrix fills and an ARPACK call):\n");
  }

  /* Initialize some communications stuff for eig_driver */

  comm = MPI_COMM_WORLD;

  /* space needed for residual vector even for matrix-only fills */
  rhs = alloc_vec();

  /* Calculate Jacobian and mass matrices */

  matrix_residual_fill_conwrap(x, rhs, MATRIX_ONLY);
  mass_matrix_fill_conwrap(x, rhs);
  create_shifted_matrix_conwrap();

  free_vec (&rhs);

  /* Set parameters for ARPACK: (1)first those that come from the input file */

  iparam[MAX_ITRS] = cei->Max_Iter;
  nev              = cei->Num_Eigenvalues;
  jmax             = cei->Num_Eigenvectors;
  sigma            = cei->Shift_Point[0];
  mu               = cei->Shift_Point[1];
  delta            = cei->Shift_Point[2];
  zeta             = cei->Shift_Point[3];
  if (delta != 1.0) {
    if (cgi->printproc > 7)
       printf("In eigensolver, delta set to one from %g\n",delta);
    delta = 1.0;
  }
  ncv              = cei->Arnoldi;
  tol              = cei->Residual_Tol[0];
  eta              = cei->Residual_Tol[1];

  /*
   * shift and invert implemented on top of Cayley, by setting flag
   * value of mu = sigma. This capability was added after the Cayley
   * version was mature, so it just piggybacks off of the Cayley
   * version. This may explain some inefficiencies, such as passing
   * the mu parameter to polez3_ where it is never used.
   * AGS 1/29/02
   */

  if (cgi->printproc > 4) {
    if (sigma==mu)
      printf("\tSHIFT-n-INVERT: sigma,ncv = %g %d\n",sigma, ncv);
    else
      printf("\tCAYLEY: sigma, mu, ncv = %g %g %d\n",sigma, mu, ncv);
  }


  /* Set parameters for ARPACK: (2)then those that are fixed for MPSalsa */

  if (sigma==mu) {
    mode = SHIFTI;
    mu = 0;
  }
  else
    mode = CAYLEY;

  which[0] = 'L';
  which[1] = 'R';
  /* which[1] = 'M'; */
  which[2] = '\0';

  bmat[0] = 'I';
  bmat[1] = '\0';

  iparam[3] = 1;
  iparam[4] = 0;
  iparam[5] = 1; /* We will check for convergence ourselves */
  iparam[7] = 0;
  iparam[8] = 0;
  iparam[9] = 0;
  iparam[10] = 0;

  iparam[SHIFTS] = 0;
  iparam[MODE] = mode;

  info = 0;

  /* Have Aztec/ARPACK Calculate Eigenvlaues */

  az_fail_cnt = eig_driver(which, bmat, iparam, mode, sigma, mu, delta, zeta,
                           nev, ncv, info, tol, eta, cgi->printproc,
                           cgi->numOwnedUnks, cgi->numUnks,
                           con->private_info.step_num, jmax, cei->sort, comm);

  if( info != 0){
    if (cgi->printproc > 1) printf("  Error %d in eigensolver\n",info);
    exit(-1);
  }

  /* Before leaving back to MPSalsa, turn off time dependent terms */

  destroy_shifted_matrix_conwrap();

  if (cgi->printproc > 1) {
    if (az_fail_cnt) {
      printf("\tWARNING: Aztec failed to reach it convergence criterion\n");
      printf("\t\t %d times during eigenvalue computation!\n",az_fail_cnt);
    }

  }
}  /* end calc_eigenvalues  */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int eig_driver(char which[], char bmat[], int iparam[], int mode,
   double sigma, double mu, double delta, double zeta, int nev, int ncv,
   int info, double tol, double eta, int printproc, int numOwnedUnks,
   int numUnks, int con_step_num, int jmax, int sort, MPI_Comm comm)

/* matShifted is temp matrix, nnz is # nonzeros in matrix */
{
  int      j, kk, ldv, lworkl;
  int      nloc, nloc_max, nloc2, ido, flag;
  int      count, nconv=0, ierr;
  int      ipntr[14];
  int      dummy1, dummy2, dummy3, dummy4;
  double   *rhs_orig;
  double   norm_M;
  char     string[4];
  int      az_fail_cnt=0;
  double   *v, *workl, *workd, *workev, *d, *resid, *vecx, *vecy, *rhs, *mxx;
  int      *select, rvec, *work;
  extern   void sort2_double(int, double *);
   /* variables for ido=4 loop */
  double *trans=NULL; /* space for eigenvalues transformed to real system */
  double *workpol=NULL; /*space for eigenvales transformed in poleze*/
  double new_sigma=0.0, new_mu = 0.0, solve_tol;
  int    temp_ncv=0, temp_nconv=0, info_p=0, nrows=nev;

   /******************************************************
    * A standard eigenvalue problem is solved (BMAT = 'I').
    * NEV is the number of eigenvalues to be approximated.
    * NCV is the number of Krylov vectors kept. WHICH
    * determine what part of the spectrum is obtained.
    * The following conditions must be satisfied:
    *                  N <= MAXN,
    *                NEV <= MAXNEV,
    *             NEV + 2 <= NCV <= MAXNCV
    *
    ******************************************************/

  nloc = numOwnedUnks;
  ldv  = nloc;

   /******************************************************
    * The work array WORKL is used in P??AUPD as
    * workspace.  Its dimension LWORKL is set as
    * illustrated below.  The parameter TOL determines
    * the stopping criterion.  If TOL<=0, machine
    * precision is used.  The variable IDO is used for
    * reverse communication and is initially set to 0.
    * Setting INFO=0 indicates that a random vector is
    * generated in PNAUPD to start the Arnoldi iteration.
    ******************************************************/

  lworkl = 3*ncv*(ncv+2);
  ido    = 0;

   /******************************************************
    * Use exact shifts with respect to the current Hessenberg
    * matrix (iparam[SHIFTS] = 1) where IPARAM[MAX_ITRS] is
    * the maximum number of Arnoldi iterations allowed.
    * Mode 1 of P??AUPD is used (IPARAM[MODE] = 1). For
    * details, see the documentation in P??AUPD.
    ******************************************************/

  nloc2  = numUnks;
  nloc_max = gmax_int_conwrap(nloc2);

  select = (int    *) malloc(ncv*sizeof(int));
  work   = (int    *) malloc(ncv*sizeof(int));
  vecx   = (double *) malloc(nloc2*sizeof(double));
  vecy   = (double *) malloc(nloc2*sizeof(double));
  rhs    = (double *) malloc(nloc2*sizeof(double));
  rhs_orig = rhs;
  mxx    = (double *) malloc(nloc2*sizeof(double));
  d      = (double *) malloc(3*ncv*sizeof(double) );
  resid  = (double *) malloc(nloc_max*sizeof(double) );
  workd  = (double *) malloc(3*nloc2*sizeof(double));
  workev = (double *) malloc(3*ncv*sizeof(double));
  workl  = (double *) calloc(lworkl,sizeof(double));
  v      = (double *) malloc(ncv*ldv*sizeof(double)) ;
  if (v == NULL) {
    fprintf(stderr,"Not enough space to allocate workl\n");
    exit(1);
  }

   /* Generate smart initial guess in resid vector by applying */
   /* resid = Inv(J)Mx  for random x (called vecx) */

  info = 1;

  if (printproc > 4) printf("\n\t    Eigensolver Initial Guess Generation\n");
  random_vector_conwrap(vecx, nloc);

   /* Mx  = rhs */
  for (kk = 0 ; kk < nloc2 ; kk++)  rhs[kk] = 0.0;
  mass_matvec_mult_conwrap(vecx,rhs);

   /* Inv(J)Mx = resid */

  norm_M = sqrt(dp(rhs,rhs));
  solve_tol = eta * norm_M / 100.0;
  if (printproc > 4)
    printf("\t    Requiring 2 extra orders resid reduction: %g\n", solve_tol);

  /* Use shifted matrix space and solver to solve usual Jacobian so that
   * the Jacobian matrix doesn't get scaled
   */

  shifted_matrix_fill_conwrap(0.0);
  shifted_linear_solver_conwrap(rhs, vecx, NEW_JACOBIAN, solve_tol);

  for (kk = 0 ; kk < nloc2 ; kk++) resid[kk] = vecx[kk];

   /* main loop */

  for (j = 0; j < 3*ncv ; j++ ) d[j] = 0.0; /* zero d[:] */
  flag = 1;
  count = 0;
  while ( flag == 1 ) {

      /*****************************************************
       * Repeatedly call P??AUPD and take actions indicated
       * by parameter IDO until either convergence is indicated
       * or the maximum iterations have been exceeded.
       *****************************************************/

     cpdnaupc_( &comm,
               &ido, bmat, &nloc, which, &nev, &tol, resid,
               &ncv, v, &ldv, iparam, ipntr, workd, workl,
               &lworkl, &info );

/*  CGS vs. MGS tests used this
 *  if (printproc > 4) printf("sigma,mu = 10,50, stopping after 3\n");
 *  sigma = 10.0;
 *  mu = 50.0;
 *  if (count==3) exit(-1);
 */

    if ( (ido == -1) || (ido == 1) ) {
      count++;
      if (printproc > 7) printf("\n\t    Eigensolver iteration: %d",count);

         /***********************************************
          * matrix vector multiplication (using inverse)
          *   workd[ipntr[1]-1] <-- OP * workd[ipntr[0]-1]
          ***********************************************/

      if (mode==CAYLEY) {
         /* OP = inv(J-sM)(J-mM) */
         /* We do the solve below, at this line rhs := (J-mM)*vecx as above */

         /*    for (kk = 0 ; kk < nloc ; kk++) vecx[kk] = workd[ipntr[0]+kk-1];
          */

               /* rhs = Mx, norm_M = ||rhs|| */
          mass_matvec_mult_conwrap(&workd[ipntr[0]-1], rhs);
          norm_M = sqrt(dp(rhs,rhs));

               /* rhs = -mMx */
          for (kk = 0 ; kk < nloc ; kk++) rhs[kk] *= -mu;

               /* vecx = Jx */
          matvec_mult_conwrap(&workd[ipntr[0]-1], vecx);

               /* rhs = (J-mM)x */
          for (kk = 0 ; kk < nloc ; kk++) rhs[kk] += vecx[kk];
      }
      else if (mode==SHIFTI) {
         /* OP = inv(J-sM)M */

               /* rhs = Mx, norm_M = ||rhs|| */
          mass_matvec_mult_conwrap(&workd[ipntr[0]-1], rhs);
          norm_M = sqrt(dp(rhs,rhs));
      }
      else {
        if (printproc > 1)
          printf("eig_driver ERROR: bad value of mode! %d\n", mode);
	exit(-1);
      }

         /*  (J-sM) vecx = rhs  */

      for (kk = 0 ; kk < nloc ; kk++)  vecx[kk] = 0.0;

      if (mode == CAYLEY || mode == SHIFTI) {

        /* set linear solver tolerance based on norm_M */
        solve_tol = eta * norm_M;
        if (printproc > 7) printf("\tLinear Solve Tol = %g\n",solve_tol);

        /* After first iter, reuse preconditioner */

        if (count%ncv == 1) {
          shifted_matrix_fill_conwrap(sigma);
          shifted_linear_solver_conwrap(rhs, vecx, NEW_JACOBIAN, solve_tol);
        }
        else {
          shifted_linear_solver_conwrap(rhs, vecx, OLD_JACOBIAN, solve_tol);
        }
      }

      for (kk = 0 ; kk < nloc ; kk++) workd[ipntr[1]+kk-1] = vecx[kk];

    }
    else if ( ido == 2) {

      /* Need to fix this if this routine is called */

      if (mode != -1){
              if (printproc > 1)
                fprintf(stderr,"\nError:dnaupd:ido=2 & mode=ord\n");
              exit(-1);
      }
         /* rhs := workd */
      for (kk = 0 ; kk < nloc ; kk++) vecx[kk] = workd[ipntr[0]+kk-1];

      mass_matvec_mult_conwrap(vecx, rhs);

      vec_copy(rhs, &(workd[ipntr[1]-1]));


    }
    else if ( ido == 3) {

        /* set shifts */

      for (kk = 0 ; kk < iparam[7] ; kk++) {
        workl[ipntr[13] -1 + kk] = 0.0;
        workl[ipntr[13] -1 + iparam[7] + kk] = 0.0;
      }
      workl[ipntr[13] -1 + iparam[7] - 1] = 1.0;

      if (iparam[7] == ncv)  {
        sigma = new_sigma;
        mu    = new_mu   ;
      }
      else {
        if (printproc > 1)
          printf("\n\tDeflation in ARPACK occured, so sigma and mu reused"
                " (%g %g)\n", sigma, mu);
      }
    }
    else if ( ido == 4) {

      if (printproc > 4) {
        printf("\n\tBefore poleze: sigma and mu = %g  %g\n", sigma, mu);
      }
      temp_nconv = nev; /* How many we want converged */
      new_sigma = sigma;
      new_mu    = mu   ;
        /* set aside storage space */
      if (trans == NULL)
        trans = (double *) malloc(3*ncv*sizeof(double));

      if (workpol == NULL)
        workpol = (double *) malloc(3*ncv*sizeof(double));
 /*
  * if sigma is less than mu, do the cayley invert and shift. if
  * sigma is greater than mu, do the cayley transform
  */

      if (mode==SHIFTI) {
        polez3_(&ncv, &new_sigma, &new_mu, &delta, &zeta, &workl[ipntr[5] -1],
              &workl[ipntr[6] -1], &workl[ipntr[7] -1], trans, trans+ncv,
              trans+2*ncv, &temp_ncv, &temp_nconv, &tol, &info_p);
      }
      /* Cayley */
      else if (sigma < mu){
        poleze_(&ncv, &new_sigma, &new_mu, &workl[ipntr[5] -1],
              &workl[ipntr[6] -1], &workl[ipntr[7] -1], trans, trans+ncv,
              trans+2*ncv, &temp_ncv, &temp_nconv, &tol, &info_p,
              workpol, workpol+ncv, workpol+2*ncv);
      }
      else if (sigma > mu){
        polez2_(&ncv, &new_sigma, &new_mu, &delta, &workl[ipntr[5] -1],
              &workl[ipntr[6] -1], &workl[ipntr[7] -1], trans, trans+ncv,
              trans+2*ncv, &temp_ncv, &temp_nconv, &tol, &info_p);
      }

      iparam[4] = temp_nconv;

      if (printproc > 4) {
        printf("\tAfter  poleze: sigma and mu = %g  %g \n\n",
                 new_sigma, new_mu);
        printf("\t%d converged of %d candidate eigenvalues found\n",
                temp_nconv, temp_ncv);
      }
      if (printproc > 7) {
        printf("\t Eigenvalues and error estimates in the lambda plane\n");
        for (kk = 0 ; kk < temp_ncv ; kk++) {
          printf("\t %2d. %g  %g  %g\n",kk, trans[kk],
                  trans[kk+ncv],trans[kk+2*ncv]);
        }

        if (info_p && printproc > 1) {
          printf("ERROR - poleze returned info = %d\n", info_p);
          if (info_p==-1) {
            printf("\tEigenvalue suspected to the right of sigma\n");
            printf("\tTry increasing sigma and/or improving accuracy\n");
          }
          else if (info_p==-2) {
            printf("\tIncrease Arnoldi space size or move shift closer\n");
          }
        }
      }
    }
    else flag = 0;
  }

   /* Either convergence or an error */

  if ( info < 0 ) {
    if ( printproc > 1 ) {
      fprintf(stderr,"\nError with _naupd, info = %d\n",info);
      fprintf(stderr,"Check documentation in _naupd.\n\n");
      nconv = 0;
      if ( info == -9999 ) {
        fprintf(stderr,"Size of Arnoldi factorization:%d\n",iparam[4]);
        fprintf(stderr,"Decrease ncv=%d to %d & rerun\n",ncv, iparam[4]);
      }
    }
  }
  else {
      /***********************************************
       * No fatal errors occurred.
       * Post-Process using PSNEUPD.
       *
       * Extract computed eigenvalues.  Eigenvectors
       * may also be computed by setting rvec = 1.
       **********************************************/

      /* Form select array, telling which eigenvalues are converged */

    if (mode==SHIFTI) {
      stslc3_( &ncv, &iparam[4], &tol, &sigma, &mu, &delta, &zeta,
               &trans[temp_ncv-iparam[4]],
             &workl[lworkl-3*ncv], &workl[lworkl-2*ncv], &workl[lworkl-ncv],
             select );
    }
    else if (sigma < mu){
      stslct_( &ncv, &iparam[4], &tol, &sigma, &mu, &trans[temp_ncv-iparam[4]],
             &workl[lworkl-3*ncv], &workl[lworkl-2*ncv], &workl[lworkl-ncv],
             select );
    }
    else{
      stslc2_( &ncv, &iparam[4], &tol, &sigma, &mu,&delta,
               &trans[temp_ncv-iparam[4]],
             &workl[lworkl-3*ncv], &workl[lworkl-2*ncv], &workl[lworkl-ncv],
             select );
    }

    rvec = 1;
    ierr = 0;
    sprintf(string,"A");
    cpdneupc_  (&comm, &rvec, string, select, d, v, &ldv, &sigma,
               &mu, &delta,  workev, bmat, &nloc, &nloc2, which, &nev,
               &tol, resid, &ncv, iparam, ipntr, workd, workl,
               &lworkl, &ierr, work);

      /*----------------------------------------------
      | The real part of the eigenvalue is returned   |
      | in the first column of the two dimensional    |
      | array D, and the imaginary part is returned   |
      | in the second column of D.  The corresponding |
      | eigenvectors are returned in the first NEV    |
      | columns of the two dimensional array V if     |
      | requested.  Otherwise, an orthogonal basis    |
      | for the invariant subspace corresponding to   |
      | the eigenvalues in D is returned in V.        |
       -----------------------------------------------*/

    if ( ierr !=  0) {
         /*-----------------------------------
         | Error condition:                   |
         | Check the documentation of PDNEUPC.|
          ------------------------------------*/
      if ( printproc > 1) {
           fprintf(stderr,"\nError with _neupc, info = %d", ierr);
           fprintf(stderr,"\nCheck the documentation of _neupc.\n\n");
/*           exit(1); */
      }
    }
    /* if poleze also worked, go ahead and print output */
    else if (!info_p) {
      nconv =  iparam[4];

    /* Do not print more rows than requested number of eigenvalues */
      if (nconv > 0 && nconv < nev) nrows = nconv;

    /* Call sort_by_real if requested, otherwise use existing order */
      if(sort) {
        sort_by_real(nconv, ncv, ldv, d, v);
      }

      for (j = 0; j < nconv ; j++ ) {

            /*--------------------------
            | Compute the residual norm |
            |                           |
            |   ||  J*x - lambda*x ||   |
            |                           |
            | for the NCONV accurately  |
            | computed eigenvalues and  |
            | eigenvectors.  (iparam(5) |
            | indicates how many are    |
            | accurate to the requested |
            | tolerance)                |
            ---------------------------*/

        if (d[j+ncv] == 0.0){
               /*-------------------
               | Ritz value is real |
               --------------------*/
       /*
        *  d[j]     : j-th eigenvalue
        *
        *  v[j*ldv] : j-th eigenvector
        */

          /* Print out eigenvectors, if requested */

          if (j < jmax) {
            if (printproc > 1)
              printf("Printing real eigenvector with value  %g\n", d[j]);
            eigenvector_output_conwrap(j, 1, &v[(j)*ldv], d[j],
                                       NULL, 0, con_step_num);
          }

          /* now calculate Rayleigh quotient */
          d[j+2*ncv] = null_vector_resid(d[j], 0.0, &v[j*ldv], NULL, TRUE);
        }
        else{
               /*----------------------
               | Ritz value is complex |
               -----------------------*/
         /*
          *  d[j]     : real part j-th eigenvalue
          *  d[j+ncv] : imag part j-th eigenvalue
          *
          *  v[j*ldv]     : real part j-th eigenvector
          *  v[(j+1)*ldv] : imag part j-th eigenvector
          */

          /* If requested, print out eigenvectors */

          if (j < jmax) {
             if (printproc > 1) printf
               ("Printing eigenvector pair for complex eigenvalues:  %g +- %g i\n",
                 d[j], fabs(d[j+ncv]));

            eigenvector_output_conwrap(j, 2, &v[(j)*ldv], d[j], &v[(j+1)*ldv],
                                       fabs(d[j+ncv]), con_step_num);
          }

          d[j+2*ncv] = null_vector_resid
                       (d[j], d[j+ncv], &v[(j)*ldv], &v[(j+1)*ldv], TRUE);
          d[j+2*ncv+1] = d[j+2*ncv];

         /* end of Rayleigh quotient stuff */

               /*-----------------------
               | Ritz value is complex. |
               | Residual of one Ritz   |
               | value of the conjugate |
               | pair is computed.      |
               ------------------------*/
          if( j+1 < nconv ){
            d[j+1]       =  d[j];
            d[j+1+ncv]   = -d[j+ncv];
            j = j + 1;
          }
        }
      }
      /*   Display computed residuals   */

      dummy1 = 6; dummy2 = 3; dummy3 = ncv; dummy4 = -6;
      cpdmout_(&comm, &dummy1, &nrows, &dummy2, d, &dummy3, &dummy4);
    }

    /*  Print additional convergence information */

    if (printproc > 4) {
      if ( info == 1 ){
        printf("\nMaximum number of iterations reached.\n");
      }
      else if ( info == 3 ){
        printf("\nNo shifts could be applied during implicit\n");
        printf("Arnoldi update, try increasing NCV.\n\n");
      }

      printf("\nEigenvalue Calculation Summary\n\n");
      printf("The number of Ritz values requested is %d\n", nev);
      printf("The number of Arnoldi vectors generated (NCV) is %d\n",ncv);
      printf("What portion of the spectrum: %s\n", which);
      printf("The number of converged Ritz values is %d\n",nconv);
      printf("Number of Implicit Arnoldi update iterations taken is %d\n",
              iparam[2]-1);
      printf("The number of OP*x is %d\n",iparam[8]);
      printf("The convergence criterion is %e\n", tol);
    }
  }

  free_vec ((double **) &select);
  free_vec ((double **) &work);
  free_vec (&vecx);
  free_vec (&vecy);
  rhs = rhs_orig;
  free_vec (&rhs);
  free_vec (&mxx);
  free_vec (&d);
  free_vec (&resid);
  free_vec (&workd);
  free_vec (&workev);
  free_vec (&workl);
  free_vec (&v);
  if (trans != NULL) free_vec (&trans);
  if (workpol != NULL) free_vec (&workpol);

  return (az_fail_cnt);
}  /* end eig_driver */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
#else
}
#endif
