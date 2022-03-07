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

#ifndef ARSCND_NAME
#define ARSCND_NAME second
#endif

C loca_eigen_cayley.F,v

C====================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C-----------------------------------------------------------------------------
C  LOCA 1.0: Library of Continuation Algorithms
C  Copyright (C) 2001, Sandia National Laboratories
C-----------------------------------------------------------------------------
C
c------------------------------------------------------------------
c This file contains 8 concatenated fortran files that will someday
c be part of Rich Lehoucq's eigen analysis library. These files add
c the capability to ARPACK to do Cayley transformations.
c
c This file was created by:
c (1) cat dnaup2c.F dnaupc.F dneupc.F pdnaup2c.F 
c         pdnaupc.F pdneupc.F poleze.F stslct.F > rf_eigen_cayley.F
c (2) all 6 instances of including stat.h and debug.h were commented out
c     and the contents of those 2 files were inserted (6 times).
c
c These 8 files are in /Net/hondo/raid1/rf_crew/people/agsalin/CAYLEY
c AGS 9/98
c--------------------------------------------------------------------
c
c\BeginDoc
c
c\Name: dnaup2c
c
c\Description: 
c  Intermediate level interface called by dnaupd.
c
c  Modified by R.B.Lehoucq to apply up to NCV implicit shifts
c  when ISHIFT=0 and to exit via reverse communication so that
c  the user may check for convergence.
c
c\Usage:
c  call dnaup2c
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, ICHK,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user has passed
c          in via reverse comunication. 0 < NP <= NCV.
c
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  ICHK    Integer.  (INPUT)
c          If IDO is not equal to 4:
c          ICHK = 1, the user will check the number of converged Ritz
c          values.
c          ICHK = 0, dnaup2c checks the number of converged Ritz values.
c          If IDO = 4 then ICHK is the number of "converged" Ritz values.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine. 
c     dnaitr  ARPACK Arnoldi factorization routine.
c     dnapps  ARPACK application of implicit shifts routine.
c     dnconv  ARPACK convergence of Ritz values routine.
c     dneigh  ARPACK compute Ritz values and error bounds routine.
c     dngets  ARPACK reorder Ritz values and error bounds routine.
c     dsortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     ARSCND_NAME  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c 
c\SCCS Information: @(#) 
c FILE: naup2.F   SID: 2.4   DATE OF SID: 7/30/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaup2c
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, ichk, 
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
#ifdef EIGEN_SERIAL
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
      real       t0, t1, t2, t3
      save       t0, t1, t2, t3
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, ichk, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift, checkc, lastis
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0, 
     &           np0, nptemp, numcnv
      Double precision
     &           rnorm, temp, eps23
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
      save       eps23, nev0, np0, kplusp, nconv, iter, getv0, update,
     &           ushift, msglvl, cnorm, checkc, initv, numcnv
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dgetv0, dnaitr, dnconv, dneigh, dngets, dnapps,
     &           dvout, ivout, ARSCND_NAME
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlapy2, dlamch
      external   ddot, dnrm2, dlapy2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      lastis = .false.
      if (ido .eq. 0) then
c 
         call ARSCND_NAME (t0)
c 
         msglvl = mnaup2
c 
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0 / 3.0D+0)
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 1
c 
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
         checkc   = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c
c     Back from checking the number of converged Ritz values.
c
      if (checkc) go to 25
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call dnaitr (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, [iter], ndigit, 
     &           '_naup2: **** Start of major iteration number ****')
         end if
c 
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine dnapps.                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, [nev], ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, [np], ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v, ldv,
     &                h, ldh, ipntr, workd, info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call dvout (logfil, 1, [rnorm], ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call dneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,
     &                q, ldq, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from dneigh.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        | NOTE: The last two arguments of dngets are no     |
c        | longer used as of version 2.1.                    |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &                bounds, workl, workl(np+1))
         if ((nev .eq. nev0+1) .and. (ichk .eq. 0)) numcnv = nev0+1
c 
c        %-------------------%
c        | Convergence test. | 
c        %-------------------%
c
         if ( ichk .eq. 0 ) then
            call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
            call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), 
     &          tol, nconv)
         else if ( ichk .eq. 1 ) then
            checkc = .true.
            ido = 4
            go to 9000
         end if        
c 
c        Back from reverse communication when the user wants to 
c        check the number of converged Ritz values.
c
 25      continue
c
         if (ido .eq. 4) then
            checkc = .false.
            nconv = ichk
            ichk = 1
         end if
c
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit, 
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call dvout (logfil, kplusp, [ritzr], ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call dvout (logfil, kplusp, [ritzi], ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call dvout (logfil, kplusp, [bounds], ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. numcnv) .or. 
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
            if (msglvl .gt. 4) then
               call dvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz eistmates computed by _neigh:')
            end if
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%

            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc(wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0
                temp = max(eps23,dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc(wprime, .true., nev0, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                temp = max(eps23, dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc(which, .true., nconv, ritzr, ritzi, bounds)
c
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &              bounds, workl, workl(np+1))
c
         end if              
c     
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, [nconv], ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call dvout (logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call dvout (logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: reverse comminucation to       |
c           | compute the shifts. They are returned in the first    |
c           | 2*NP locations of WORKL.                              |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            np = kplusp
            do 55 j=1, kplusp
               if (bounds(j) .eq. zero) then
                  np = np - 1
               end if
 55            continue
            go to 9000
         end if
c 
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:2*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZR, RITZI to free up WORKL    |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call dcopy (np, workl,       1, ritzr, 1)
             call dcopy (np, workl(np+1), 1, ritzi, 1)
             nev = kplusp - np
             if ( nev.eq.0 ) then 
                nev = 1
                np = np - 1
                lastis = .true.
             else
                lastis = .false.
             end if
         end if
c
         if (msglvl .gt. 2) then 
            call ivout (logfil, 1, [np], ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call dvout (logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call dvout (logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 ) 
     &          call dvout (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call dnapps (n, nev, np, ritzr, ritzi, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
         if ( lastis ) then
            nev = 0
            np = np + 1
            call daxpy ( n, h(1,1)-ritzr(np), v(1,1), 1, resid, 1 )
         end if
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call ARSCND_NAME (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call ARSCND_NAME (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
c
         if (msglvl .gt. 2) then
            call dvout (logfil, 1, [rnorm], ndigit, 
     &      '_naup2: B-norm of residual for compressed factorization')
            call dmout (logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
c
         go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = numcnv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call ARSCND_NAME (t1)
      tnaup2 = t1 - t0
c     
 9000 continue
c
c     %---------------%
c     | End of dnaup2c |
c     %---------------%
c
#endif
      return
      end
c\BeginDoc
c
c\Name: dnaupc
c
c\Description: 
c
c  This version has been modified by R.B. Lehoucq to allow:
c  * a Cayley transformation to be used 
c  * the user to pass in up to NCV implicit shifts.
c  * the user to determine the number of "converged" Ritz values.
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs 
c  of a linear operator "OP" with respect to a semi-inner product defined 
c  by a symmetric positive semi-definite real matrix B. B may be the 
c  identity matrix. NOTE: If the linear operator "OP" is real and symmetric 
c  with respect to the real positive semi-definite symmetric matrix B, 
c  i.e. B*OP = (OP')*B, then subroutine dsaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dnaupc is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = inv[A - sigma*M]*M   and  B = M or I. 
c           ===> shift-and-invert mode with a real shift.
c           If OP*x = amu*x, then 
c           amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = inv[A - sigma*M]*[A - mu*M]  and  B = M or I. 
c           ===> Cayley transformation mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           lambda = (sigma*amu-mu)/(amu-1)
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dnaupc
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dnaupc.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dnaupc with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: The IPARAM(8) real and imaginary parts 
c                    of the implicit shifts where INPTR(14) 
c                    is the pointer into WORKL for placing the 
c                    shifts. See Remark 5 below.
c          IDO =  4: The user returns the number of "converged"
c                    Ritz values in IPARAM(5).  See Remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
c          values are kept together. (See remark 4 below)
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL 
c                      array corresponding to RITZR and RITZI. See remark 
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = 0 if dnaupc is to check the number of converged 
c                        Ritz values.
c                    otherwise ido = 4 and the user will detemine the
c                    number of converged Ritz values. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = ICHK.
c          ICHK = 1; the user wants to check for the number of converged 
c                    Ritz values.
c          ICHK = 0; dnaupc will check for convergence. 
c          See remark 6 below.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of dnaupc for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3, the user provides shifts through reverse
c          communication (IPARAM(1)=0). dnaupc returns NP, the maximum
c          number of shifts the user is to provide. 0 < NP <= NCV. The
c          user must make sure and possibly re-set IPARAM(8) with the 
c          number of shifts. See Remark 5 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array 
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the 
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine dneupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: IDO=3 and IPARAM(8) is returned larger than the
c                 value provided by dnaupc.f. See remark 5.
c          = -14: IDO=3 and IPARAM(8) is returned set equal to NCV
c                 and the imaginary portion of the last pair of shifts
c                 passed in is a complex conjugate pair. See remark 5.
c          = -15: IDO=4 and IPARAM(5) is less than 0 or greater than NCV.
c          = -16: IPARAM(6) must be equal to 0 or 1 when IDO=0.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine dneupd.
c
c  2. If a basis for the invariant subspace corresponding to the "converged" 
c     Ritz values is needed, the user must call dneupd immediately following 
c     completion of dnaupc. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL'
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
c     linear systems should be solved with L and L' rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L'z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that 
c     NCV >= NEV + 2. However, it is recommended that NCV .ge. 2*NEV+1.  
c     If many problems of the same type are to be solved,  one should 
c     experiment with increasing NCV while keeping NEV fixed for a given 
c     test problem.  This will usually decrease the required number of OP*x 
c     operations but it also increases the work and storage required to 
c     maintain the orthogonal basis vectors.  The optimal "cross-over" with 
c     respect to CPU time is problem dependent and must be determined 
c     empirically. See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user provides up to 
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs 
c     must be placed in consecutive locations except in the NP location
c     when NP=NCV. The real part of the 
c     eigenvalues of the current upper Hessenberg matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located 
c     in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c     The user must also set IPARAM(8) to be equal to the number of
c     shifts passed in (i.e. degree of the polynomial filter to be
c     implicitly applied.)
c     Note that NP may be returned by dnaupc less than NCV because 
c     a leading block of the current upper Hessenberg matrix has split 
c     off and hence contains "converged" Ritz values.
c
c  6. If IPARAM(6) = ICHK = 1 and IDO=4, the user needs to determine the 
c     number of "converged" Ritz values. The user must also set the set
c     IPARAM(5)=NCONV equal to the number that are "converged". The real 
c     part of the eigenvalues of the current upper Hessenberg matrix are 
c     located in WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the 
c     imaginary part in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). 
c     They are ordered according to the order defined by WHICH. The complex 
c     conjugate pairs are kept together and the associated Ritz estimates 
c     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , 
c     WORKL(IPNTR(8)+NCV-1).
c     
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     dnaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     ivout   ARPACK utility routine that prints integers.
c     ARSCND_NAME  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#) 
c FILE: naupc.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaupc
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, info )
#ifdef EIGEN_SERIAL
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
      real       t0, t1
      save       t0, t1
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, ichk, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, ichk, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dnaup2, dvout, ivout, ARSCND_NAME, dstatn
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call dstatn
         call ARSCND_NAME (t0)
         msglvl = mnaupd
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
         levec  = iparam(2)
         mxiter = iparam(3)
         nb     = iparam(4)
         ichk   = iparam(6)
         mode   = iparam(7)
c
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (ichk .lt. 0 .or. ichk .gt. 1) then
                                                ierr = -16
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = dlamch('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine dneigh called   |
c        | by dnaup2. Subroutine dneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw 
c
      end if
c
c     Do some error checking when ido=3.
c
      if (ido .eq. 3 .and. iparam(8) .gt. np ) then
         info = -13
         ido  = 99
         go to 9000
      else if (ido.eq.3 .and. iparam(8).eq.ncv .and. 
     &         workl(ipntr(14)+2*ncv-1).eq.-workl(ipntr(14)+2*ncv-2)
     &         .and. workl(ipntr(14)+ncv-1).eq.workl(ipntr(14)+ncv-2) )
     &         then
         info = -14
         ido  = 99
         go to 9000
      else if (ido .eq. 3) then
         np = iparam(8)
      end if
c
c     Do some error checking when ido=4.
c
      if (ido .eq. 4 .and. iparam(5) .ge. 0 
     &    .and. iparam(5) .le. ncv ) then
         ichk = iparam(5)
      else if (ido .eq. 4) then
         info = -15
         ido  = 99
         go to 9000
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call dnaup2c 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, ichk,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), 
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      iparam(5) = np
      iparam(3) = mxiter      
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
      if (ido .ne. 99) go to 9000
c 
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dnaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [mxiter], ndigit,
     &               '_naupc: Number of update iterations taken')
         call ivout (logfil, 1, [np], ndigit,
     &               '_naupc: Number of wanted "converged" Ritz values')
         call dvout (logfil, np, workl(ritzr), ndigit, 
     &               '_naupc: Real part of the final Ritz values')
         call dvout (logfil, np, workl(ritzi), ndigit, 
     &               '_naupc: Imaginary part of the final Ritz values')
         call dvout (logfil, np, workl(bounds), ndigit, 
     &               '_naupc: Associated Ritz estimates')
      end if
c
      call ARSCND_NAME (t1)
      tnaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number: ', ' 2.4', 21x, ' =',/
     &      5x, '= Version Date:   ', ' 07/31/96', 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in naup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
      end if
c
 9000 continue
c
#endif
      return
c
c     %---------------%
c     | End of dnaupc |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dneupc
c
c\Description: 
c
c  This version has been modified by R.B. Lehoucq to allow a Cayley
c  transformation to be used and for the user to determine which
c  Ritz values are "converged".
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors. An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to DNAUPC.  DNAUPC must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to 
c  as a Schur basis.
c
c  See documentation in the header of the subroutine DNAUPC for definition
c  of OP as well as other terms and the relation of computed Ritz values
c  and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of IPARAM(7),
c  MODE and WHICH in the documentation of DNAUPC.
c
c\Usage:
c  call dneupc 
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMA, MU, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c 
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c          If IPARAM(6)=ICHK=1, then SELECT(j) must be set to .true. if
c          WORKL(IPNTR(6)+j-1) and WORKL(IPNTR(7)+j-1) corresponds to a
c          Ritz value that the user has decided is "converged". Otherwise
c          SELECT(J) must be set to .false..
c
c  DR      Double precision array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
c          the real part of the Ritz  approximations to the eigenvalues of 
c          A*z = lambda*B*z. 
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by 
c          DNAUPC. A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPC to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Double precision N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. 
c          (OUTPUT) On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns
c          of Z represent approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem 
c          A*z = lambda*B*z. 
c 
c          The complex Ritz vector associated with the Ritz value 
c          with positive imaginary part is stored in two consecutive 
c          columns.  The first column holds the real part of the Ritz 
c          vector and the second column holds the imaginary part.  The 
c          Ritz vector associated with the Ritz value with negative 
c          imaginary part is simply the complex conjugate of the Ritz vector 
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPC.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMA  Double precision  (INPUT)
c         Represents the pole of the spectral transformation used if
c         IPARAM(7) = 3 or 4.
c
c  MU     Double precision  (INPUT)
c         If IPARAM(7) = 4, represents the zero of the Cayley transformation
c         imaginary part of the shift. 
c         Not referenced if IPARAM(7) = 1,2 or 3. 
c
c  WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPC that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPC following the last call
c         to DNAUPC.  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPC and the call to DNEUPC.
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPC.
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dNAUPC.  They are not changed by dneupc.
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupc.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupc if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen.
c                Re-enter subroutine dneupc with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least dimension NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPC did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: IPARAM(7)=4 and sigma = mu. Not allowed.
c          = -16: IPARAM(6) must be equal to 0 or 1.
c          = -17: IPARAM(6)=1 but the number of locations of SELECT
c                 is not equal to IPARAM(5). IPARAM(5) needs to contain
c                 the number of "converged" Ritz values as decided by 
c                 the user.
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dgeqr2  LAPACK routine that computes the QR factorization of 
c             a matrix.
c     dlacpy  LAPACK matrix copy routine.
c     dlahqr  LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset  LAPACK matrix initialization routine.
c     dorm2r  LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen  LAPACK routine that re-orders the Schur form.
c     dtrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     dger    Level 2 BLAS rank one update to a matrix.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let X' denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the real 
c     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics          
c     Rice University           
c     Houston, Texas            
c 
c\SCCS Information: @(#) 
c FILE: neupc.F   SID: 2.5   DATE OF SID: 7/31/96   RELEASE: 2 
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dneupc (rvec, howmny, select, dr, di, z, ldz, sigma, 
     &             mu, delta, workev, bmat, n, which, nev, tol, 
     &                   resid, ncv, v, ldv, iparam, ipntr, workd, 
     &                   workl, lworkl, info)
#ifdef EIGEN_SERIAL
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     real       t0, t1, t2, t3, t4, t5
c     save       t0, t1, t2, t3, t4, t5
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision     
     &           sigma, mu, delta, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1), di(nev+1), resid(n), v(ldv,ncv), z(ldz,*), 
     &           workd(3*n), workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      complex*16 c1, c2, c3, c4, c5
      integer    bounds, ierr, ih, ihbds, iheigr, iheigi, iconj, nconv, 
     &           invsub, iuptri, iwev, iwork(1), j, k, ktrord, ichk,
     &           ldh, ldq, mode, msglvl, outncv, ritzr, ritzi, wri, wrr,
     &           irr, iri, ibd
      logical    reord
      Double precision
     &           conds, rnorm, sep, temp, thres, vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dger, dgeqr2, dlacpy, dlahqr, dlaset, dmout, 
     &           dorm2r, dtrevc, dtrmm, dtrsen, dscal, dvout, ivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dnrm2, dlamch, ddot
      external   dlapy2, dnrm2, dlamch, ddot
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt, dble, dcmplx, dimag
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      ichk = iparam(6)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 .and. mu .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'CAYLEY'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G') then
         ierr = -11
      else if (mode .eq. 4 .and. (sigma .eq. mu)) then
         ierr = -15
      else if (ichk .lt. 0 .or. ichk .gt. 1) then
                                              ierr = -16
      end if
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPC.                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c     
      if (rvec) then
c
c        If the user checks for the number of "converged" Ritz
c        values, then SELECT has been set.
c
         if (ichk.eq.1) then
            ktrord = 0
            reord = .false.
            do 15 j=1, ncv
               if (select(j)) then
                  ktrord = ktrord+1
                  if (j .gt. nconv ) reord = select(j)
               end if
 15         continue
            if (ktrord .ne. nconv) then
       print*,-17,ktrord,nconv
               info = -17
               go to 9000
            end if
            go to 25
         end if 
c     
c        %-------------------------------------------%
c        | Get converged Ritz value on the boundary. |
c        | Note: converged Ritz values have been     |
c        | placed in the first NCONV locations in    |
c        | workl(ritzr) and workl(ritzi).  They have |
c        | been sorted (in _naup2) according to the  |
c        | WHICH selection criterion.                |
c        %-------------------------------------------%
c
         if (which .eq. 'LM' .or. which .eq. 'SM') then
            thres = dlapy2( workl(ritzr), workl(ritzi) )
         else if (which .eq. 'LR' .or. which .eq. 'SR') then
            thres = workl( ritzr )
         else if (which .eq. 'LI' .or. which .eq. 'SI') then
            thres = abs( workl(ritzi) )
         end if
c
         if (msglvl .gt. 2) then
            call dvout(logfil, 1, [thres], ndigit,
     &           '_neupc: Threshold eigenvalue used for re-ordering')
         end if
c
c        %----------------------------------------------------------%
c        | Check to see if all converged Ritz values appear at the  |
c        | top of the upper quasi-triangular matrix computed by     |
c        | _neigh in _naup2.  This is done in the following way:    |
c        |                                                          |
c        | 1) For each Ritz value obtained from _neigh, compare it  |
c        |    with the threshold Ritz value computed above to       |
c        |    determine whether it is a wanted one.                 |
c        |                                                          | 
c        | 2) If it is wanted, then check the corresponding Ritz    |
c        |    estimate to see if it has converged.  If it has, set  |
c        |    correponding entry in the logical array SELECT to     |
c        |    .TRUE..                                               |
c        |                                                          |
c        | If SELECT(j) = .TRUE. and j > NCONV, then there is a     |
c        | converged Ritz value that does not appear at the top of  |
c        | the upper quasi-triangular matrix computed by _neigh in  |
c        | _naup2.  Reordering is needed.                           |
c        %----------------------------------------------------------%
c
         reord = .false.
         ktrord = 0
         do 10 j = 0, ncv-1
            select(j+1) = .false.
            if (which .eq. 'LM') then
               if (dlapy2(workl(irr+j), workl(iri+j)) .ge. thres) then
                  temp1 = max( eps23, 
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SM') then
               if (dlapy2(workl(irr+j), workl(iri+j))
     &            .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'LR') then
               if (workl(irr+j) .ge. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SR') then
               if (workl(irr+j) .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'LI') then
               if (abs(workl(iri+j)) .ge. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SI') then
               if (abs(workl(iri+j)) .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            end if
            if (j+1 .gt. nconv ) reord = ( select(j+1) .or. reord )
            if (select(j+1)) ktrord = ktrord + 1
 10      continue 
c
c        Branch here directly if the user has set select.
c
 25      continue
c
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, [ktrord], ndigit,
     &            '_neupc: Number of specified eigenvalues')
             call ivout(logfil, 1, [nconv], ndigit,
     &            '_neupc: Number of "converged" eigenvalues')
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPC.        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c     
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv, zero, one, workl(invsub), ldq)
         call dlahqr (.true., .true., ncv, 1, ncv, workl(iuptri), ldh,
     &        workl(iheigr), workl(iheigi), 1, ncv, 
     &        workl(invsub), ldq, ierr)
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c     
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c     
         if (msglvl .gt. 1) then
            call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupc: Real part of the eigenvalues of H')
            call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupc: Imaginary part of the Eigenvalues of H')
            call dvout (logfil, ncv, workl(ihbds), ndigit,
     &           '_neupc: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call dmout (logfil, ncv, ncv, workl(iuptri), ldh, ndigit,
     &              '_neupc: The upper quasi-triangular matrix ')
            end if
         end if 
c
         if (reord) then
c     
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. | 
c           %-----------------------------------------------------%
c     
            call dtrsen ('None', 'V', select, ncv, workl(iuptri), ldh, 
     &           workl(invsub), ldq, workl(iheigr), workl(iheigi), 
     &           nconv, conds, sep, workl(ihbds), ncv, iwork, 1, ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
            if (msglvl .gt. 2) then
                call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupc: Real part of the eigenvalues of H--reordered')
                call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupc: Imag part of the eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call dmout (logfil, ncv, ncv, workl(iuptri), ldq, 
     &                  ndigit,
     &              '_neupc: Quasi-triangular matrix after re-ordering')
                end if
            end if
c     
         end if
         call dcopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then 
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
c     
c        %---------------------------------------------------%
c        | Compute the QR factorization of the Schur matrix. |
c        %---------------------------------------------------%
c     
         call dgeqr2 (ncv, ncv, workl(invsub), ldq, workev, 
     &        workev(ncv+1), ierr)
         call dscal (ncv-1, 0.0d0, workl(ihbds), 1)
         workl(ihbds+ncv-1) = one
         call dorm2r( 'Left', 'Transpose', ncv, 1, ncv, 
     &        workl(invsub), ldq, workev, workl(ihbds), ncv,
     &        workev(ncv+1), ierr)
         if (ierr .ne. 0) then
             info = -9
             go to 9000
         end if
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r.                     |   
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n, ncv, ncv,
     &        workl(invsub), ldq, workev, v, ldv, workd(n+1), ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, ncv
c     
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c     
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (ncv, -one, workl(iuptri+j-1), ldq)
               call dscal (ncv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c     
 20      continue
c     
         if (howmny .eq. 'A') then
c     
c           %-----------------------------------%
c           | Compute the NCV eigenvectors of T | 
c           | located in workl(iuptri,ldq).     |
c           %-----------------------------------%
c     
            call dtrevc ('Right', 'All', select, ncv, workl(iuptri), 
     &           ldq, vl, 1, workl(invsub), ldq, ncv, outncv, workev,
     &           ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c     
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c     
            iconj = 0
            do 40 j=1, ncv
c
               if ( workl(iheigi+j-1) .eq. zero ) then
c     
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c     
                  temp = dnrm2( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal ( ncv, one / temp, 
     &                 workl(invsub+(j-1)*ldq), 1 )
                  workev(j) = ddot( j,  workl(ihbds), 1, 
     &                              workl(invsub+(j-1)*ldq), 1 )
c
               else
c     
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2( dnrm2( ncv, workl(invsub+(j-1)*ldq),
     &                      1 ), dnrm2( ncv, workl(invsub+j*ldq),  1) )
                     call dscal ( ncv, one / temp, 
     &                      workl(invsub+(j-1)*ldq), 1 )
                     call dscal ( ncv, one / temp, 
     &                      workl(invsub+j*ldq), 1 )
                     workev(j) = ddot( j+1,  workl(ihbds), 1, 
     &                                 workl(invsub+(j-1)*ldq), 1 )
                     workev(j+1) = ddot( j+1,  workl(ihbds), 1, 
     &                                   workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
c           %-----------------------------------------------------%
c           | Copy last row of the eigenvectors into workl(ihbds) |
c           %-----------------------------------------------------%
c
            call dcopy(ncv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c     
            call dgeqr2 (ncv, nconv, workl(invsub), ldq, workev, 
     &                   workev(ncv+1), ierr)
c     
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |   
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  | 
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c     
            call dorm2r ('Right', 'Notranspose', n, ncv, nconv,
     &           workl(invsub), ldq, workev, z, ldz, workd(n+1), ierr)
c     
            call dtrmm ('Right', 'Upper', 'No transpose', 'Non-unit',
     &                  n, nconv, one, workl(invsub), ldq, z, ldz)
c     
         end if
c     
      else 
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPC into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
c 
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call dscal (ncv, rnorm, workl(ihbds), 1)     
c     
      else 
c     
c        %-------------------------------------%
c        | A spectral transformation was used. |
c        | Determine the Ritz estimates of the |
c        | Ritz values in the original system. |
c        %-------------------------------------%
c     
         if (rvec) then
c     
c           %------------------------------------------%
c           | Save the eigenvalues and last components |
c           | of the eigenevctors of H.                |
c           %------------------------------------------%
c     
            call dcopy (ncv, workl(ihbds), 1, workev, 1) 
            call dcopy (ncv, workl(iheigr), 1, workev(ncv+1), 1) 
            call dcopy (ncv, workl(iheigi), 1, workev(2*ncv+1), 1) 
c
            call dscal (ncv, rnorm, workl(ihbds), 1)
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workl(ihbds+j-1) = dlapy2( workl(ihbds+j-1), 
     &                                          workl(ihbds+j) )
                     workl(ihbds+j) = workl(ihbds+j-1)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               else
                  workl(ihbds+j-1) = abs( workl(ihbds+j-1) )
               end if
 45         continue
c
         end if
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workev, ndigit,
     &              '_neupc: Ritz estimates')
            call dvout (logfil, ncv, workev(ncv+1), ndigit,
     &              '_neupc: real part of the Ritz values')
            call dvout (logfil, ncv, workev(2*ncv+1), ndigit,
     &              '_neupc: imaginary part of the Ritz values')
         end if
c
         if (type .eq. 'SHIFTI') then
c
            do 50 k=1, ncv
               temp = dlapy2( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = workl(ihbds+k-1) / temp
 50         continue
c
         else if (type .eq. 'CAYLEY' ) then
c
            do 60 k=1, ncv
               c1 = dcmplx( workl(iheigr+k-1)-one, workl(iheigi+k-1) )
               c3 =  c1 * c1 
               temp = dlapy2( dble( c3 ), dimag( c3 ) )
               workl(ihbds+k-1) = workl(ihbds+k-1) / temp
 60         continue
c
         end if
c
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.| 
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c     
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               temp = dlapy2( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1) / temp / temp 
     &                           + sigma
               workl(iheigi+k-1) = -workl(iheigi+k-1) / temp / temp
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         else if (type .eq. 'CAYLEY') then
c
            c1 = dcmplx( sigma, 0.0D+0 )
            c2 = dcmplx( mu, 0.0D+0 )
            do 90 k=1, ncv
               c3 = dcmplx( workl(iheigr+k-1), workl(iheigi+k-1) )
               c4 = dcmplx( delta*workl(iheigr+k-1)-one, 
     &              delta*workl(iheigi+k-1) )
               c5 = (c1*c3 - c2)/c4
               workl(iheigr+k-1) = dble( c5 )
               workl(iheigi+k-1) = dimag( c5 )
 90         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         end if
c
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupc: Untransformed real part of the Ritz values.')
         call dvout (logfil, nconv, di, ndigit,
     &   '_neupc: Untransformed imag part of the Ritz values.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupc: Ritz estimates of untransformed Ritz values.')
      else if (type .eq. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupc: Real parts of converged Ritz values.')
         call dvout (logfil, nconv, di, ndigit,
     &   '_neupc: Imag parts of converged Ritz values.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupc: Associated Ritz estimates.')
      end if
c 
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------%
c
      if ( rvec .and. howmny .eq. 'A' .and. 
     &     type .eq. 'SHIFTI' ) then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j = 1, nconv
            if ( workev(2*ncv+j) .eq. zero ) then
               workev(j) = workev(j) / workev(ncv+j)
            else if (iconj .eq. 0 ) then
               iconj = 1
               c1 = dcmplx( workev(ncv+j), workev(2*ncv+j) )
               c2 = dcmplx( workev(j), workev(j+1) )
               c3 = c2 / c1
               workev(j) = dble( c3 )
               workev(j+1) = dimag( c3 )
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      else if ( rvec .and. howmny .eq. 'A' .and.
     &          type .eq. 'CAYLEY' ) then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / (theta - 1)     |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 120 j = 1, nconv
            if ( workev(2*ncv+j) .eq. zero ) then
               workev(j) = workev(j) / (workev(ncv+j)-one)
            else if (iconj .eq. 0 ) then
               iconj = 1
               c1 = dcmplx( workev(ncv+j)-one, workev(2*ncv+j) )
               c2 = dcmplx( workev(j), workev(j+1) )
               c3 = c2 / c1
               workev(j) = dble( c3 )
               workev(j+1) = dimag( c3 )
            else
               iconj = 0
            end if
 120     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
#endif
      return
c     
c     %---------------%
c     | End of DNEUPC |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: pdnaup2c
c
c Message Passing Layer: MPI
c
c\Description: 
c  Intermediate level interface called by pdnaupd.
c
c  Modified by R.B.Lehoucq to apply up to NCV implicit shifts
c  when ISHIFT=0 and to exit via reverse communication so that
c  the user may check for convergence.
c
c\Usage:
c  call pdnaup2c
c     ( COMM, IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, ICHK,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  COMM, IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in pdnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in pdnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user has passed
c          in via reverse comunication. 0 < NP <= NCV.
c
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  ICHK    Integer.  (INPUT)
c          If IDO is not equal to 4:
c          ICHK = 1, the user will check the number of converged Ritz
c          values.
c          ICHK = 0, pdnaup2c checks the number of converged Ritz values.
c          If IDO = 4 then ICHK is the number of "converged" Ritz values.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     pdgetv0  Parallel ARPACK initial vector generation routine. 
c     pdnaitr  Parallel ARPACK Arnoldi factorization routine.
c     pdnapps  Parallel ARPACK application of implicit shifts routine.
c     dnconv   ARPACK convergence of Ritz values routine.
c     pdneigh  Parallel ARPACK compute Ritz values and error bounds routine.
c     pdngets  Parallel ARPACK reorder Ritz values and error bounds routine.
c     dsortc   ARPACK sorting routine.
c     pivout   Parallel ARPACK utility routine that prints integers.
c     ARSCND_NAME   ARPACK utility routine for timing.
c     pdmout   Parallel ARPACK utility routine that prints matrices
c     pdvout   ARPACK utility routine that prints vectors.
c     pdlamch  ScaLAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     ddot     Level 1 BLAS that computes the scalar product of two vectors. 
c     pdnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
c     dswap    Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              Cray Research, Inc. &
c     Dept. of Computational &     CRPC / Rice University
c     Applied Mathematics          Houston, Texas
c     Rice University           
c     Houston, Texas    
c 
c\Revision history:
c     Starting Point: Serial Code FILE: naup2.F   SID: 2.2
c
c\SCCS Information: 
c FILE: naup2.F   SID: 1.3   DATE OF SID: 8/13/96   
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine pdnaup2c
     &   ( comm, ido, bmat, n, which, nev, np, tol, resid, mode, ichk, 
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
#ifdef EIGEN_PARALLEL
c
      include   'mpif.h'
c
c     %---------------%
c     | MPI Variables |
c     %---------------%
c
      integer    comm
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
      real       t0, t1, t2, t3
      save       t0, t1, t2, t3
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, ichk, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift, checkc, lastis
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0, 
     &           np0, nptemp, numcnv
      Double precision
     &           rnorm, temp, eps23
c
      Double precision
     &           rnorm_buf
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
      save       eps23, nev0, np0, kplusp, nconv, iter, getv0, update,
     &           ushift, msglvl, cnorm, checkc, initv, numcnv
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, pdgetv0, pdnaitr, dnconv, 
     &           pdneigh, pdngets, pdnapps,
     &           pdvout, pivout, ARSCND_NAME
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, pdnorm2, dlapy2, pdlamch
      external   ddot, pdnorm2, dlapy2, pdlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
         call ARSCND_NAME (t0)
c 
         msglvl = mnaup2
c 
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = pdlamch(comm, 'Epsilon-Machine')
         eps23 = eps23**(2.0 / 3.0)
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 1
c 
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
         checkc   = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call pdgetv0 (comm, ido, bmat, 1, initv, n, 1, v, ldv, 
     &                resid, rnorm, ipntr, workd, workl, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c
c     Back from checking the number of converged Ritz values.
c
      if (checkc) go to 25
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call pdnaitr (comm, ido, bmat, n, 0, nev, mode, 
     &             resid, rnorm, v, ldv, h, ldh, ipntr, 
     &             workd, workl, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call pivout (comm, logfil, 1, iter, ndigit, 
     &           '_naup2: **** Start of major iteration number ****')
         end if
c 
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine pdnapps.                 |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call pivout (comm, logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call pivout (comm, logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call pdnaitr (comm, ido, bmat, n, nev, np, mode, 
     &                resid, rnorm, v, ldv,
     &                h, ldh, ipntr, workd, workl, info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call pdvout (comm, logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call pdneigh ( comm, rnorm, kplusp, h, ldh, ritzr, ritzi, 
     &                  bounds, q, ldq, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from pdneigh.                      |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        | NOTE: The last two arguments of pdngets are no    |
c        | longer used as of version 2.1.                    |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call pdngets ( comm, ishift, which, nev, np, ritzr, ritzi, 
     &                  bounds, workl, workl(np+1))
         if ((nev .eq. nev0+1) .and. (ichk .eq. 0)) numcnv = nev0+1
c 
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         if ( ichk .eq. 0 ) then
            call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
            call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), 
     &          tol, nconv)
         else if ( ichk .eq. 1 ) then
            checkc = .true.
            ido = 4
            go to 9000
         end if        
c 
c        Back from reverse communication when the user wants to 
c        check the number of converged Ritz values.
c
 25      continue
c
         if (ido .eq. 4) then
            checkc = .false.
            nconv = ichk
            ichk = 1
         end if
c 
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call pivout (comm, logfil, 4, kp, ndigit, 
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call pdvout (comm, logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call pdvout (comm, logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call pdvout (comm, logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. numcnv) .or. 
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c     
            if (msglvl .gt. 4) then
               call dvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz estimates computed by _neigh:')
            end if
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%
 
            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc(wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0
                temp = max(eps23,dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc(wprime, .true., nev0, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                temp = max(eps23, dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc(which, .true., nconv, ritzr, ritzi, bounds)
c
c
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call pdngets(comm, ishift, which, nev, np, ritzr, ritzi, 
     &                      bounds, workl, workl(np+1))
c
         end if              
c     
         if (msglvl .gt. 0) then
            call pivout (comm, logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call pivout (comm, logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call pdvout (comm, logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call pdvout (comm, logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call pdvout (comm, logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: reverse comminucation to       |
c           | compute the shifts. They are returned in the first    |
c           | 2*NP locations of WORKL.                              |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            np = kplusp
            do 55 j=1, kplusp
               if (bounds(j) .eq. zero) then
                  np = np - 1
               end if
 55            continue
            go to 9000
         end if
c 
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:2*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZR, RITZI to free up WORKL    |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call dcopy (np, workl,       1, ritzr, 1)
             call dcopy (np, workl(np+1), 1, ritzi, 1)
             nev = kplusp - np
             if ( nev.eq.0 ) then 
                nev = 1
                np = np - 1
                lastis = .true.
             else
                lastis = .false.
             end if
         end if
c
         if (msglvl .gt. 2) then 
            call pivout (comm, logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call pdvout (comm, logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call pdvout (comm, logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 ) 
     &          call pdvout (comm, logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call pdnapps (comm, n, nev, np, ritzr, ritzi, v, ldv, 
     &                 h, ldh, resid, q, ldq, workl, workd)
c
         if ( lastis ) then
            nev = 0
            np = np + 1
            call daxpy ( n, h(1,1)-ritzr(np), v(1,1), 1, resid, 1 )
         end if
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to pdnaitr. |
c        %---------------------------------------------%
c
         cnorm = .true.
         call ARSCND_NAME (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call ARSCND_NAME (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif
c
         if (bmat .eq. 'G') then
            rnorm_buf = ddot (n, resid, 1, workd, 1)
            call MPI_ALLREDUCE( rnorm_buf, rnorm, 1,
     &                MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = pdnorm2( comm, n, resid, 1 )
         end if
         cnorm = .false.
c
         if (msglvl .gt. 2) then
            call pdvout (comm, logfil, 1, rnorm, ndigit, 
     &      '_naup2: B-norm of residual for compressed factorization')
            call pdmout (comm, logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
c 
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = numcnv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call ARSCND_NAME (t1)
      tnaup2 = t1 - t0
c     
 9000 continue
c
c     %----------------%
c     | End of pdnaup2c |
c     %----------------%
c
#endif
      return
      end
c\BeginDoc
c
c\Name: pdnaupc
c
c\Description: 
c
c  This version has been modified by R.B. Lehoucq to allow:
c  * a Cayley transformation to be used 
c  * the user to pass in up to NCV implicit shifts.
c  * the user to determine the number of "converged" Ritz values.
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs 
c  of a linear operator "OP" with respect to a semi-inner product defined by 
c  a symmetric positive semi-definite real matrix B. B may be the identity 
c  matrix. NOTE: If the linear operator "OP" is real and symmetric 
c  with respect to the real positive semi-definite symmetric matrix B, 
c  i.e. B*OP = (OP')*B, then subroutine ssaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  pdnaupc is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = inv[A - sigma*M]*M   and  B = M or I. 
c           ===> shift-and-invert mode with a real shift.
c           If OP*x = amu*x, then 
c           amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = inv[A - sigma*M]*[A - mu*M]  and  B = M or I. 
c           ===> Cayley transformation mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           lambda = (sigma*amu-mu)/(amu-1)
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call pdnaupc
c     ( COMM, IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid.  (INPUT)
c
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to pdnaupc.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          pdnaupc with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: The IPARAM(8) real and imaginary parts 
c                    of the implicit shifts where INPTR(14) 
c                    is the pointer into WORKL for placing the 
c                    shifts. See Remark 5 below.
c          IDO =  4: The user returns the number of "converged"
c                    Ritz values in IPARAM(5).  See Remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
c          values are kept together. (See remark 4 below)
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL 
c                      array corresponding to RITZR and RITZI. See remark 
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = 0 if dnaupc is to check the number of converged 
c                        Ritz values.
c                    otherwise ido = 4 and the user will detemine the
c                    number of converged Ritz values. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = ICHK.
c          ICHK = 1; the user wants to check for the number of converged 
c                    Ritz values.
c          ICHK = 0; dnaupc will check for convergence. 
c          See remark 6 below.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of pdnaupc for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3, the user provides shifts through reverse
c          communication (IPARAM(1)=0). dnaupc returns NP, the maximum
c          number of shifts the user is to provide. 0 < NP <= NCV. The
c          user must make sure and possibly re-set IPARAM(8) with the 
c          number of shifts. See Remark 5 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array 
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by pdneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the 
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     pdneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine pdneupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: IDO=3 and IPARAM(8) is returned larger than the
c                 value provided by dnaupc.f. See remark 5.
c          = -14: IDO=3 and IPARAM(8) is returned set equal to NCV
c                 and the imaginary portion of the last pair of shifts
c                 passed in is a complex conjugate pair. See remark 5.
c          = -15: IDO=4 and IPARAM(5) is less than 0 or greater than NCV.
c          = -16: IPARAM(6) must be equal to 0 or 1 when IDO=0.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine pdneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call pdneupd immediately following 
c     completion of pdnaupc. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL'
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
c     linear systems should be solved with L and L' rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L'z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
c     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user provides up to 
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs 
c     must be placed in consecutive locations except in the NP location
c     when NP=NCV. The real part of the 
c     eigenvalues of the current upper Hessenberg matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located 
c     in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c     The user must also set IPARAM(8) to be equal to the number of
c     shifts passed in (i.e. degree of the polynomial filter to be
c     implicitly applied.)
c     Note that NP may be returned by dnaupc less than NCV because 
c     a leading block of the current upper Hessenberg matrix has split 
c     off and hence contains "converged" Ritz values.
c
c  6. If IPARAM(6) = ICHK = 1 and IDO=4, the user needs to determine the 
c     number of "converged" Ritz values. The user must also set the set
c     IPARAM(5)=NCONV equal to the number that are "converged". The real 
c     part of the eigenvalues of the current upper Hessenberg matrix are 
c     located in WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the 
c     imaginary part in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). 
c     They are ordered according to the order defined by WHICH. The complex 
c     conjugate pairs are kept together and the associated Ritz estimates 
c     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , 
c     WORKL(IPNTR(8)+NCV-1).
c
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     pdnaup2  Parallel ARPACK routine that implements the implicitly 
c              restarted Arnoldi Iteration.
c     pivout   Parallel ARPACK utility routine that prints integers.
c     ARSCND_NAME   ARPACK utility routine for timing.
c     pdvout   Parallel ARPACK utility routine that prints vectors.
c     pdlamch  ScaLAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#) 
c FILE: naupc.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine pdnaupc
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &     iparam, ipntr, workd, workl, lworkl, info )
#ifdef EIGEN_PARALLEL
c
      include  'mpif.h'
c
c     %------------------%
c     | MPI Variables    |
c     %------------------%
c
      integer    comm, myid 
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
      real       t0, t1
      save       t0, t1
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, ichk, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, ichk, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   pdnaup2, pdvout, pivout, ARSCND_NAME, dstatn
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           pdlamch
      external   pdlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call dstatn
         call ARSCND_NAME (t0)
         msglvl = mnaupd
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
         levec  = iparam(2)
         mxiter = iparam(3)
         nb     = iparam(4)
         ichk   = iparam(6)
         mode   = iparam(7)
c
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (ichk .lt. 0 .or. ichk .gt. 1) then
                                                ierr = -16
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = pdlamch(comm,'EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine dneigh called   |
c        | by pdnaup2. Subroutine dneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw 
c
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
c     if ido = 3, then set np = iparam(8) and set nev0 accordingly.
c
c
c     Do some error checking when ido=3.
c
      if (ido .eq. 3 .and. iparam(8) .gt. np ) then
         info = -13
         ido  = 99
         go to 9000
      else if (ido.eq.3 .and. iparam(8).eq.ncv .and. 
     &         workl(ipntr(14)+2*ncv-1).eq.-workl(ipntr(14)+2*ncv-2)
     &         .and. workl(ipntr(14)+ncv-1).eq.workl(ipntr(14)+ncv-2) )
     &         then
         info = -14
         ido  = 99
         go to 9000
      else if (ido .eq. 3) then
         np = iparam(8)
      end if
c
c     Do some error checking when ido=4.
c
      if (ido .eq. 4 .and. iparam(5) .ge. 0 
     &    .and. iparam(5) .le. ncv ) then
         ichk = iparam(5)
      else if (ido .eq. 4) then
         info = -15
         ido  = 99
         go to 9000
      end if
c
      call pdnaup2c
     &   ( comm, ido, bmat, n, which, nev0, np, tol, resid, mode, ichk,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), 
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
      if (ido .ne. 99) go to 9000
c 
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within pdnaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
      if (msglvl .gt. 0) then
         call pivout (comm, logfil, 1, mxiter, ndigit,
     &               '_naupd: Number of update iterations taken')
         call pivout (comm, logfil, 1, np, ndigit,
     &               '_naupd: Number of wanted "converged" Ritz values')
         call pdvout (comm, logfil, np, workl(ritzr), ndigit, 
     &               '_naupd: Real part of the final Ritz values')
         call pdvout (comm, logfil, np, workl(ritzi), ndigit, 
     &               '_naupd: Imaginary part of the final Ritz values')
         call pdvout (comm, logfil, np, workl(bounds), ndigit, 
     &               '_naupd: Associated Ritz estimates')
      end if
c
      call ARSCND_NAME (t1)
      tnaupd = t1 - t0
c
      if (msglvl .gt. 0) then
         call MPI_COMM_RANK( comm, myid, ierr )
         if ( myid .eq. 0 ) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number: ', ' 2.4', 21x, ' =',/
     &      5x, '= Version Date:   ', ' 07/31/96', 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in p_naup2 routine              = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
         end if
      end if
c
 9000 continue
c
#endif
      return
c
c     %---------------%
c     | End of pdnaupc |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: pdneupc
c
c Message Passing Layer: MPI
c
c\Description: 
c
c  This version has been modified by R.B. Lehoucq to allow a Cayley
c  transformation to be used and for the user to determine which
c  Ritz values are "converged".
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to PDNAUPC.  PDNAUPC must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine PDNAUPC for 
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of 
c  IPARAM(7), MODE and WHICH in the documentation of PDNAUPC.
c
c\Usage:
c  call pdneupc 
c     ( COMM, RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMA, MU, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
c       LWORKL, INFO )
c
c\Arguments:
c  COMM    MPI Communicator for the processor grid.  (INPUT)
c
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c 
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
c          the real part of the Ritz  approximations to the eigenvalues of 
c          A*z = lambda*B*z. 
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by 
c          PDNAUPC. A further computation must be performed by the user
c          to transform the Ritz values computed for OP by PDNAUPC to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Double precision N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. 
c          (OUTPUT) On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns
c          of Z represent approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem 
c          A*z = lambda*B*z. 
c 
c          The complex Ritz vector associated with the Ritz value 
c          with positive imaginary part is stored in two consecutive 
c          columns.  The first column holds the real part of the Ritz 
c          vector and the second column holds the imaginary part.  The 
c          Ritz vector associated with the Ritz value with negative 
c          imaginary part is simply the complex conjugate of the Ritz vector 
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by PDNAUPC.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMA  Double precision  (INPUT)
c         Represents the pole of the spectral transformation used if
c         IPARAM(7) = 3 or 4.
c
c  MU     Double precision  (INPUT)
c         If IPARAM(7) = 4, represents the zero of the Cayley transformation
c         imaginary part of the shift. 
c         Not referenced if IPARAM(7) = 1,2 or 3. 
c
c  WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to PDNAUPC that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to PDNEUPC following the last call
c         to PDNAUPC.  These arguments MUST NOT BE MODIFIED between
c         the the last call to PDNAUPC and the call to PDNEUPC.
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by PDNAUPC .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          pdnaupc.  They are not changed by pdneupc.
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by pdneupc.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     pdneupc if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen.
c                Re-enter subroutine pdneupc with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least dimension NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: PDNAUPC did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: IPARAM(7)=4 and sigma = mu. Not allowed.
c          = -16: IPARAM(6) must be equal to 0 or 1.
c          = -17: IPARAM(6)=1 but the number of locations of SELECT
c                 is not equal to IPARAM(5). IPARAM(5) needs to contain
c                 the number of "converged" Ritz values as decided by 
c                 the user.
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     pivout  Parallel ARPACK utility routine that prints integers.
c     pdmout  Parallel ARPACK utility routine that prints matrices
c     pdvout  Parallel ARPACK utility routine that prints vectors.
c     dgeqr2  LAPACK routine that computes the QR factorization of 
c             a matrix.
c     dlacpy  LAPACK matrix copy routine.
c     dlahqr  LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     pdlamch ScaLAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset  LAPACK matrix initialization routine.
c     dorm2r  LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen  LAPACK routine that re-orders the Schur form.
c     dtrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     dger    Level 2 BLAS rank one update to a matrix.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let X' denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the real 
c     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics          
c     Rice University           
c     Houston, Texas            
c 
c\SCCS Information: @(#) 
c FILE: neupc.F   SID: 2.5   DATE OF SID: 7/31/96   RELEASE: 2 
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine pdneupc (comm, rvec, howmny, select, dr, di, z, ldz, 
     &               sigma, mu, delta, workev, bmat, n, which, nev, tol, 
     &                    resid, ncv, v, ldv, iparam, ipntr, workd, 
     &                    workl, lworkl, info)
#ifdef EIGEN_PARALLEL
c
c     %--------------------%
c     | MPI Communicator |
c     %--------------------%
c
      integer   comm
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cc      include   'debug.h'
cc      include   'stat.h'
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     real       t0, t1, t2, t3, t4, t5
c     save       t0, t1, t2, t3, t4, t5
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
c
c\SCCS Information: @(#) 
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision     
     &           sigma, mu, delta, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1), di(nev+1), resid(n), v(ldv,ncv), z(ldz,*), 
     &           workd(3*n), workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      complex*16 c1, c2, c3, c4, c5
      integer    bounds, ierr, ih, ihbds, iheigr, iheigi, iconj, nconv, 
     &           invsub, iuptri, iwev, iwork(1), j, k, ktrord, ichk,
     &           ldh, ldq, mode, msglvl, outncv, ritzr, ritzi, wri, wrr,
     &           irr, iri, ibd
      logical    reord
      Double precision
     &           conds, rnorm, sep, temp, thres, vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dger, dgeqr2, dlacpy, dlahqr, dlaset, pdmout, 
     &           dorm2r, dtrevc, dtrmm, dtrsen, dscal, pdvout, pivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dnrm2, pdlamch, ddot
      external   dlapy2, dnrm2, pdlamch, ddot
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt, dble, dcmplx, dimag
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      ichk = iparam(6)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = pdlamch(comm, 'Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 .and. mu .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'CAYLEY'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G') then
         ierr = -11
      else if (mode .eq. 4 .and. (sigma .eq. mu)) then
         ierr = -15
      else if (ichk .lt. 0 .or. ichk .gt. 1) then
                                              ierr = -16
      end if
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by PDNEUPC.                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c     
      if (rvec) then
c
c        If the user checks for the number of "converged" Ritz
c        values, then SELECT has been set.
c
         if (ichk.eq.1) then
            ktrord = 0
            reord = .false.
            do 15 j=1, ncv
               if (select(j)) then
                  ktrord = ktrord+1
                  if (j .gt. nconv ) reord = select(j)
               end if
 15         continue
            if (ktrord .ne. nconv) then
       print*,-17,ktrord,nconv
               info = -17
               go to 9000
            end if
            go to 25
         end if 
c     
c        %-------------------------------------------%
c        | Get converged Ritz value on the boundary. |
c        | Note: converged Ritz values have been     |
c        | placed in the first NCONV locations in    |
c        | workl(ritzr) and workl(ritzi).  They have |
c        | been sorted (in _naup2) according to the  |
c        | WHICH selection criterion.                |
c        %-------------------------------------------%
c
         if (which .eq. 'LM' .or. which .eq. 'SM') then
            thres = dlapy2( workl(ritzr), workl(ritzi) )
         else if (which .eq. 'LR' .or. which .eq. 'SR') then
            thres = workl(ritzr)
         else if (which .eq. 'LI' .or. which .eq. 'SI') then
            thres = abs( workl(ritzi) )
         end if
c
         if (msglvl .gt. 2) then
            call pdvout(comm, logfil, 1, thres, ndigit,
     &           '_neupc: Threshold eigenvalue used for re-ordering')
         end if
c
c        %----------------------------------------------------------%
c        | Check to see if all converged Ritz values appear at the  |
c        | top of the upper quasi-triangular matrix computed by     |
c        | _neigh in _naup2.  This is done in the following way:    |
c        |                                                          |
c        | 1) For each Ritz value obtained from _neigh, compare it  |
c        |    with the threshold Ritz value computed above to       |
c        |    determine whether it is a wanted one.                 |
c        |                                                          | 
c        | 2) If it is wanted, then check the corresponding Ritz    |
c        |    estimate to see if it has converged.  If it has, set  |
c        |    correponding entry in the logical array SELECT to     |
c        |    .TRUE..                                               |
c        |                                                          |
c        | If SELECT(j) = .TRUE. and j > NCONV, then there is a     |
c        | converged Ritz value that does not appear at the top of  |
c        | the upper quasi-triangular matrix computed by _neigh in  |
c        | _naup2.  Reordering is needed.                           |
c        %----------------------------------------------------------%
c
         reord = .false.
         ktrord = 0
         do 10 j = 0, ncv-1
            select(j+1) = .false.
            if (which .eq. 'LM') then
               if (dlapy2(workl(irr+j), workl(iri+j)) .ge. thres) then
                  temp1 = max( eps23, 
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SM') then
               if (dlapy2(workl(irr+j), workl(iri+j)) .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'LR') then
               if (workl(irr+j) .ge. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SR') then
               if (workl(irr+j) .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'LI') then
               if (abs(workl(iri+j)) .ge. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            else if (which .eq. 'SI') then
               if (abs(workl(iri+j)) .le. thres) then
                  temp1 = max( eps23,
     &                         dlapy2( workl(irr+j), workl(iri+j) ) )
                  if (workl(ibd+j) .le. tol*temp1)
     &               select(j+1) = .true.
               end if
            end if
            if (j+1 .gt. nconv ) reord = ( select(j+1) .or. reord )
            if (select(j+1)) ktrord = ktrord + 1
 10      continue 
c
c        Branch here directly if the user has set select.
c
 25      continue
c
         if (msglvl .gt. 2) then
             call pivout(comm, logfil, 1, ktrord, ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call pivout(comm, logfil, 1, nconv, ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by PDNAUPC.        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c     
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv, zero, one, workl(invsub), ldq)
         call dlahqr (.true., .true., ncv, 1, ncv, workl(iuptri), ldh,
     &        workl(iheigr), workl(iheigi), 1, ncv, 
     &        workl(invsub), ldq, ierr)
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c     
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c     
         if (msglvl .gt. 1) then
            call pdvout (comm, logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H')
            call pdvout (comm, logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imaginary part of the Eigenvalues of H')
            call pdvout (comm, logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call pdmout (comm, logfil, ncv, ncv, workl(iuptri), ldh, 
     &              ndigit,'_neupd: The upper quasi-triangular matrix ')
            end if
         end if 
c
         if (reord) then
c     
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. | 
c           %-----------------------------------------------------%
c     
            call dtrsen ('None', 'V', select, ncv, workl(iuptri), ldh, 
     &           workl(invsub), ldq, workl(iheigr), workl(iheigi), 
     &           nconv, conds, sep, workl(ihbds), ncv, iwork, 1, ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
            if (msglvl .gt. 2) then
                call pdvout (comm, logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H--reordered')
                call pdvout (comm, logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imag part of the eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call pdmout (comm, logfil, ncv, ncv, workl(iuptri), 
     &                  ldq, ndigit,
     &              '_neupd: Quasi-triangular matrix after re-ordering')
                end if
            end if
c     
         end if
         call dcopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then 
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
c     
c        %---------------------------------------------------%
c        | Compute the QR factorization of the Schur matrix. |
c        %---------------------------------------------------%
c     
         call dgeqr2 (ncv, ncv, workl(invsub), ldq, workev, 
     &        workev(ncv+1), ierr)
         call dscal (ncv-1, 0.0d0, workl(ihbds), 1)
         workl(ihbds+ncv-1) = one
         call dorm2r( 'Left', 'Transpose', ncv, 1, ncv, 
     &        workl(invsub), ldq, workev, workl(ihbds), ncv,
     &        workev(ncv+1), ierr)
         if (ierr .ne. 0) then
             info = -9
             go to 9000
         end if
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r.                     |   
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n, ncv, ncv,
     &        workl(invsub), ldq, workev, v, ldv, workd(n+1), ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, ncv
c     
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c     
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (ncv, -one, workl(iuptri+j-1), ldq)
               call dscal (ncv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c     
 20      continue
c     
         if (howmny .eq. 'A') then
c     
c           %-----------------------------------%
c           | Compute the NCV eigenvectors of T | 
c           | located in workl(iuptri,ldq).     |
c           %-----------------------------------%
c     
            call dtrevc ('Right', 'All', select, ncv, workl(iuptri), 
     &           ldq, vl, 1, workl(invsub), ldq, ncv, outncv, workev,
     &           ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c     
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c     
            iconj = 0
            do 40 j=1, ncv
c
               if ( workl(iheigi+j-1) .eq. zero ) then
c     
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c     
                  temp = dnrm2( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal ( ncv, one / temp, 
     &                 workl(invsub+(j-1)*ldq), 1 )
                  workev(j) = ddot( j,  workl(ihbds), 1, 
     &                              workl(invsub+(j-1)*ldq), 1 )
c
               else
c     
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2( dnrm2( ncv, workl(invsub+(j-1)*ldq),
     &                      1 ), dnrm2( ncv, workl(invsub+j*ldq),  1) )
                     call dscal ( ncv, one / temp, 
     &                      workl(invsub+(j-1)*ldq), 1 )
                     call dscal ( ncv, one / temp, 
     &                      workl(invsub+j*ldq), 1 )
                     workev(j) = ddot( j+1,  workl(ihbds), 1, 
     &                                 workl(invsub+(j-1)*ldq), 1 )
                     workev(j+1) = ddot( j+1,  workl(ihbds), 1, 
     &                                   workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
c           %-----------------------------------------------------%
c           | Copy last row of the eigenvectors into workl(ihbds) |
c           %-----------------------------------------------------%
c
            call dcopy(ncv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c     
            call dgeqr2 (ncv, nconv, workl(invsub), ldq, workev, 
     &                   workev(ncv+1), ierr)
c     
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |   
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  | 
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c     
            call dorm2r ('Right', 'Notranspose', n, ncv, nconv,
     &           workl(invsub), ldq, workev, z, ldz, workd(n+1), ierr)
c     
            call dtrmm ('Right', 'Upper', 'No transpose', 'Non-unit',
     &                  n, nconv, one, workl(invsub), ldq, z, ldz)
c     
         end if
c     
      else 
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed PDNAUPC into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
c 
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call dscal (ncv, rnorm, workl(ihbds), 1)     
c     
      else 
c     
c        %-------------------------------------%
c        | A spectral transformation was used. |
c        | Determine the Ritz estimates of the |
c        | Ritz values in the original system. |
c        %-------------------------------------%
c     
         if (rvec) then
c     
c           %------------------------------------------%
c           | Save the eigenvalues and last components |
c           | of the eigenevctors of H.                |
c           %------------------------------------------%
c     
            call dcopy (ncv, workl(ihbds), 1, workev, 1) 
            call dcopy (ncv, workl(iheigr), 1, workev(ncv+1), 1) 
            call dcopy (ncv, workl(iheigi), 1, workev(2*ncv+1), 1) 
            call dscal (ncv, rnorm, workl(ihbds), 1)
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workl(ihbds+j-1) = dlapy2( workl(ihbds+j-1), 
     &                                          workl(ihbds+j) )
                     workl(ihbds+j) = workl(ihbds+j-1)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               else
                  workl(ihbds+j-1) = abs( workl(ihbds+j-1) )
               end if
 45         continue
c
         end if
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workev, ndigit,
     &              '_neupc: Ritz estimates')
            call dvout (logfil, ncv, workev(ncv+1), ndigit,
     &              '_neupc: real part of the Ritz values')
            call dvout (logfil, ncv, workev(2*ncv+1), ndigit,
     &              '_neupc: imaginary part of the Ritz values')
         end if
c
         if (type .eq. 'SHIFTI') then
c
            do 50 k=1, ncv
               temp = dlapy2( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = workl(ihbds+k-1) / temp
 50         continue
c
         else if (type .eq. 'CAYLEY' ) then
c
            do 60 k=1, ncv
               c1 = dcmplx( workl(iheigr+k-1)-one, workl(iheigi+k-1) )
               c3 =  c1 * c1 
               temp = dlapy2( dble( c3 ), dimag( c3 ) )
               workl(ihbds+k-1) = workl(ihbds+k-1) / temp
 60         continue
c
         end if
c
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.| 
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c     
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               temp = dlapy2( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1) / temp / temp 
     &                           + sigma
               workl(iheigi+k-1) = -workl(iheigi+k-1) / temp / temp
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         else if (type .eq. 'CAYLEY') then
c
            c1 = dcmplx( sigma, 0.0D+0 )
            c2 = dcmplx( mu, 0.0D+0 )
            do 90 k=1, ncv
               c3 = dcmplx( workl(iheigr+k-1), workl(iheigi+k-1) )
               c4 = dcmplx( delta* workl(iheigr+k-1)-one, 
     &              delta*workl(iheigi+k-1) )
               c5 = (c1*c3 - c2)/c4
               workl(iheigr+k-1) = dble( c5 )
               workl(iheigi+k-1) = dimag( c5 )
 90         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         end if
c
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call pdvout (comm, logfil, nconv, dr, ndigit,
     &   '_neupd: Untransformed real part of the Ritz valuess.')
         call pdvout (comm, logfil, nconv, di, ndigit,
     &   '_neupd: Untransformed imag part of the Ritz valuess.')
         call pdvout (comm, logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Ritz estimates of untransformed Ritz values.')
      else if (type .eq. 'REGULR' .and. msglvl .gt. 1) then
         call pdvout (comm, logfil, nconv, dr, ndigit,
     &   '_neupd: Real parts of converged Ritz values.')
         call pdvout (comm, logfil, nconv, di, ndigit,
     &   '_neupd: Imag parts of converged Ritz values.')
         call pdvout (comm, logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Associated Ritz estimates.')
      end if
c 
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------%
c
      if ( rvec .and. howmny .eq. 'A' .and. 
     &     type .eq. 'SHIFTI' ) then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j = 1, nconv
            if ( workev(2*ncv+j) .eq. zero ) then
               workev(j) = workev(j) / workev(ncv+j)
            else if (iconj .eq. 0 ) then
               iconj = 1
               c1 = dcmplx( workev(ncv+j), workev(2*ncv+j) )
               c2 = dcmplx( workev(j), workev(j+1) )
               c3 = c2 / c1
               workev(j) = dble( c3 )
               workev(j+1) = dimag( c3 )
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      else if ( rvec .and. howmny .eq. 'A' .and.
     &          type .eq. 'CAYLEY' ) then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / (theta - 1)     |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 120 j = 1, nconv
            if ( workev(2*ncv+j) .eq. zero ) then
               workev(j) = workev(j) / (workev(ncv+j)-one)
            else if (iconj .eq. 0 ) then
               iconj = 1
               c1 = dcmplx( workev(ncv+j)-one, workev(2*ncv+j) )
               c2 = dcmplx( workev(j), workev(j+1) )
               c3 = c2 / c1
               workev(j) = dble( c3 )
               workev(j+1) = dimag( c3 )
            else
               iconj = 0
            end if
 120     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
#endif
      return
c     
c     %---------------%
c     | End of PDNEUPC |
c     %---------------%
c
      end
c
c   Name: poleze.F (m, sigma, mu, realpt, imagpt, ritzes, trealpt,
c                timagpt, tritzes, tm, nconv, tol, info, workreal,
c                workimag, writzes )
c   Description:
c   This routine takes eigenvalues nu=(lambda-mu)/(lambda-sigma) and
c   transforms them to the lambda plane. It determines how many
c   of the lambdas satisfy the approximation criterion. It picks
c   new values for mu and sigma to be used in the next iteration of pdnaupc
c
c
c   SUBROUTINE ARGUMENTS
c
c    m (integer):        ncv (INPUT)
c    sigma (double):     pole of the Cayley transform (nu) (INPUT/OUTPUT)
c    mu (double):        root of the Cayley transform (nu) (INPUT/OUTPUT)
c    realpt (array(tm)): the real part of the numerically
c                        computed eigenvalue (nu plane) (INPUT/OUTPUT)
c    imagpt (array(tm)): the imaginary part of the numerically
c                        computed eigenvalue (nu)  (INPUT/OUTPUT)
c    ritzes (array(tm)): the error in the computed eigenvalue
c                        (nu) (INPUT/OUTPUT)
c    trealpt (array(tm)):the real part of the transformed eigenvalue in the
c                        lambda plane, sorted by largest real part (OUTPUT)
c    timagpt (array(tm)):the imag part of the transformed eigenvalue
c                       n the lambda plane, sorted by largest real part (OUTPUT)
c    tritzes (array(tm)): the reordered ritz values (OUTPUT)
c    tm (integer):    number of "candidate evals" in lambda plane (OUTPUT)
c    nconv (integer): on input, number of eigenvalues requested,
c                     on output, number of converged lambdas (INPUT/OUTPUT)
c                     Iparam(5) in f77; Iparam(4) in C
c    tol (double) :   convergence criterion for eigenvalues (INPUT)
c    info (integer):  provides error information (OUTPUT)
c         info = 1:  could not achieve criterion for distance between x1, x2
c         info = 0 :  routine has performed normally
c         info = -1 :  the choice of sigma was not to the
c                     right of the largest eigenvalue in the lambda plane
c         info = -2 :  the number of candidate eigenvalues is less than
c                     the number of converged eigenvalues
c         info = -3 : could not find which transformed eigenvalue in the
c                     nu plane corresponds to the eigenvalue of largest real
c                     part in the lambda plane
c         info = -4 : all candidate eigenvalues lambda have the same real part
c         info = -5 : nconv on exit must be <= ncv/2 
c
c    workreal (array(tm)): work array used to hold real part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    workimag (array(tm)): work array used to hold imaginary part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    writzes (array(tm)):  work array to hold ritz values
c
c
c    VARIABLES
c
c    zsigma (complex):    sigma for complex calculations
c    zmu (complex):       mu for complex calculations
c    c3,c4,c5 (complex):  work scalars
c    mx (double):         used in computing a new shift
c    mn (double):         used in computing a new shift
c    x1 (double):         real part of rightmost eval in lambda plane
c    x2 (double):         real part of the [max(tm, inncnv+1)]th
c                         largest magnitude eigenvalue in the nu plane
c    inncnv (integer):    the number of eigenvalues (lambda) we wish to
c                         have converged 
c    jj (integer):        a counter
c
c    CONSTANTS
c
c    one: 1.0
c    zero: 0.0
c    stripw: .15
c
c Note: the above constants are used in transforming the Cayley eigenvalues
c from the nu to the lambda plane. "stripw" is the strip width within which
c we do not consider Cayley transformed eigenvalues (nu), for they
c correspond to infinite eigenvalues (lambda)
c
c  Routines called
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dsortc  ARPACK sorting routine.
c             uses the following:
c             'LM' -> want the NEV eigenvalues of largest magnitude.
c             'SM' -> want the NEV eigenvalues of smallest magnitude.
c             'LR' -> want the NEV eigenvalues of largest real part.
c             'SR' -> want the NEV eigenvalues of smallest real part.
c             'LI' -> want the NEV eigenvalues of largest imaginary part.
c             'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c-----------------------------------------------------------------------
      subroutine poleze ( m, sigma, mu, realpt, imagpt, ritzes, trealpt,
     &           timagpt, tritzes, tm, nconv, tol, info,
     &           workreal, workimag, writzes )
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)

      integer m, tm, nconv, info
      double precision realpt(*), imagpt(*), ritzes(*), trealpt(*),
     &     timagpt(*), tritzes(*),  sigma, mu, tol,
     &     workreal(*), workimag(*), writzes(*)
      complex*16 zsigma, zmu, c3, c4, c5
      double precision temp, mx, mn, x1, x2, y1
      integer j, inncnv, jj
c
      double precision one, zero, stripw
      parameter (one = 1.0D+0, zero = 0.0D+0, stripw=.15)
c
      double precision dlapy2
      external   dlapy2
c
      intrinsic  dcmplx, dimag, dble, max, min
c
c     Start executable statements
c
      info = 0
      zsigma = dcmplx( sigma, zero )
      zmu = dcmplx( mu, zero )
      tm = 0
      inncnv = nconv
c
c     Eigenvalues nu=(lambda - mu)/(lambda-sigma) are transformed
c     to the lambda plane, and the values of the real and imaginary
c     parts of lambdas are held in trealpt, timagpt.
c     The loop runs through m times, which is the number of eigenvalues
c     requested. It only stores tm values of lambda, requiring
c     that the real part of nu is greater than one.  
c
      do 10 j=1, m
       if (  realpt(j) .gt. one ) then
            c3 = dcmplx( realpt(j), imagpt(j) )
            c4 = dcmplx( realpt(j)-one, imagpt(j) )
            c5 = (zsigma*c3 - zmu)/c4
            tm = tm + 1
            trealpt(tm) = dble( c5 )
            timagpt(tm) = dimag( c5 )
            tritzes(tm) = ritzes(j)
c
c     if the magnitude of the eigenvalue nu is outside the unit circle
c     with real part less than one,
c     we have hit the wrong branch of the curve
c
       elseif (dlapy2(realpt(j), imagpt(j)).gt. one+stripw) then
         info=-1
       endif
 10   continue
c
c     Sort the lambda eigenvalues and the error bounds, store in trealpt,
c     timagpt, tritzes
c     'LR' gives the eigenvalues of largest real part
c
      call dsortc ('LR', .true., tm, trealpt, timagpt, tritzes)
c
c    error if the number of candidate eigenvalues is
c    less than or equal to the number of eigenvalues we want converged
c
      if ( tm .le. nconv ) then
         info = -2
         return
      endif
c
c     Determine how many consecutive lambda eigenvalues
c     satisfy the error tolerance and store the result in nconv.
c     This is regardless of how many eigenvalues were requested
c     (inncnv), but will not return more than m/2+1. 
c
      jj = 0
      flag = 0  
      j=tm
 15   if(flag.eq.1) goto 20
        if (j.le.tm .and. j.ge.1) then
         c3 = dcmplx( trealpt(j)-mu, timagpt(j) )
         c4 = dcmplx( trealpt(j)-sigma, timagpt(j) )
         c5 = c3 / c4
         temp = dlapy2( dble(c5), dimag(c5) )
         if  (tritzes(j) .le. tol*temp) then
           jj = jj + 1
         else
           flag=1
         endif
c
c     ensures loop exit
c
         if (j.eq.1) flag=1
         j= j-1
        endif
        goto 15
 20   continue
c
c     jj is the counter; decides how many evals are within
c     the tolerance and assigns the value to nconv
c
c    nconv now means how many consecutive evals are converged, not
c    how many we wanted. inncnv contains how many evals we wanted.
c
        nconv = jj
c
c     x1 is the real part of the rightmost eigenvalue in the lambda plane.
c     store y1, corresponding imaginary part, in order to track this eval later
c
      x1 = trealpt(tm)
      y1 = timagpt(tm)
c
c     Copy the tm transformed lambda eigenvalues from trealpt, timagpt
c     into workreal, workimag in order to sort
c
      do 40 j=1, tm
         workreal(j)=trealpt(j)
         workimag(j)=timagpt(j)
         writzes(j)=tritzes(j)
 40   continue
c
c     Transform the tm lambdas from workreal & workimag back to the
c     nu plane and store the result in workreal & workimag
c
      do 45 j=1, tm
         c3 = dcmplx( workreal(j), workimag(j) ) - zsigma
         c4 = dcmplx( workreal(j), workimag(j) ) - zmu
         c5 = c4 / c3
         workreal(j) = dble( c5 )
         workimag(j) = dimag( c5 )
 45   continue
c
c     Sort the nu eigenvalues from largest magnitude to smallest
c     Note that SM gives the eigenvalues 
c     sorted from largest magnitude to smallest
c
      call dsortc ('SM', .true., tm, workreal, workimag, writzes)
c
c     Re-Transform the tm Cayley transformed eigenvalues back to the
c     lambda plane just for output and store results in workreal, workimag
c
      do 50 j=1, tm
            c3 = dcmplx( workreal(j), workimag(j) )
            c4 = dcmplx( workreal(j)-one, workimag(j) )
            c5 = (zsigma*c3 - zmu) / c4
            workreal(j) = dble( c5 )
            workimag(j) = dimag( c5 )
 50   continue
c
c     x2 is either the real part of lambda corresponding to the
c     transformed eigenvalue nu of inncnv+1 th largest magnitude, or
c     if that is identical to the value x1, instead we pick the
c     real part of lambda that corresponds to the nu which has magnitude
c     next larger than the nu corresponding to x1.
c
c     j is a counter for this loop
c
      j=tm
c
c     we want to find which position in the nu plane that x1 (from
c     the lambda plane) got mapped to. When this first if-statement
c     is satisfied, the value of j corresponds to that desired position
c     Unless there is an error, this first statement will eventually
c     be satisfied and we will pick the desired value of x2
c
c     A remark in terms of ordering: workreal & workimag hold
c     the values of lambda sorted by magnitude of nu, where the
c     largest  magnitude is in the top position and the smallest
c     magnitude is in the last position. treal & timag hold
c     the values of lambda sorted by real part of lambda, where the
c     leftmost (smallest) is in the top position and the
c     rightmost (largest) is in the last position
c
 60   if (abs(x1 - workreal(j)) .le. (1.0d-10)*(dlapy2(x1,y1)+1.0))
     &  then
         jj = max(j+1, inncnv+1)
 70      x2 = workreal(min(tm, jj))
c
c     Now a separate embedded if-loop to guard against error cases
c
c     first check that x2 satisfies |x1-x2| > dlapy(x1,y1)*10^-7. if not
c     increment jj and repick x2
c
       if ((abs(x1-x2) .le. dlapy2(x1,x2)*1.0d-7) .and. jj .ne. 1) then
          jj=jj-1
          goto 70
c
c     this is the case that x1 equals x2, info = -4 and end subroutine
c
        elseif ((abs(x1 - x2) .le. (1.0d-10)*dlapy2(x1,y1))
     &   .and. (jj .eq. 1)) then
            info = -4
            return
c
c     this is the case that none of the x1, x2 are far enough apart
c     not a fatal error, but flag the user with info = 1
c
        elseif ((abs(x1-x2) .le. dlapy2(x1,x2)*1.0d-7) .and.
     &           jj .eq. 1) then
           info = 1
        endif
c
c     this continues the loop looking for the correct position (j)
c
      elseif ((abs(x1 - workreal(j)) .gt. (1.0d-10)*dlapy2(x1,y1))
     & .and. (j .ne. 1)) then
         j=j-1
         goto 60
c
c     this is the case that we can't find where x1 got mapped. info = -3
c      and end subroutine
c
      else
         info = -3
         return
      endif
c
c     Pick a new zero and pole for the Cayley transformation. x2 is selected
c     so | nu - 1 | .ge. mn for the inncnv eigenvalues of interest
c     and the rightmost eigenvalue in the lambda plane. The rightmost
c     eigenvalue may not satisfy this if it is a Hopf bifurcation
c     with large imaginary part and there are other eigenvalues closer
c     to the pole sigma.
c
c     pick sigma and mu to satisfy the system:
c     (x1-mu)/(x1-sigma) = mx
c     (x2-mu)/(x1-sigma)= mn +1
c    (this choice of mn + 1 is partly motivated by an error bound)
c
      mx = 10.0
      mn = 2.0
      sigma  = ( (mx - one) * x1 - mn * x2 ) / (mx - mn - one)
      mu     = ( mn + one ) * sigma - mn * x2
#endif
      return
      end

c-----------------------------------------------------------------------

      subroutine stslct ( m, nconv, tol, sigma, mu, cutoff, ritzr, 
     &           ritzi, errbds, select )
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
      integer m, nconv
      double precision cutoff, tol, sigma, mu
      double precision ritzr(m), ritzi(m), errbds(m)
      integer select(m)
c
      complex*16 zsigma, zmu, c3, c4, c5
      integer j
c     
      double precision one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0 )
c     
      double precision dlapy2
      external   dlapy2
c     
      intrinsic  dcmplx, dble, abs
c
c     Start executable statements
c
      zsigma = dcmplx( sigma, zero )
      zmu = dcmplx( mu, zero )
      do 10 j=1, m
         c3 = dcmplx( ritzr(j), ritzi(j) )
         c4 = dcmplx( ritzr(j)-one, ritzi(j) )
         c5 = (zsigma*c3 - zmu) / c4
         if ( (dble( c5 ) - cutoff .ge. -1.0d-10*abs(cutoff))
     & .and. (dble(c5) .lt.sigma)) then
c
c .and.  &        (errbds(j) .le. tol*dlapy2( ritzr(j), ritzi(j) )) ) then
c
            select(j)=1
         else
            select(j)=0
         end if
 10   continue
c
#endif
      return
      end

c
c   Name: polez2 (m, sigma, mu, delta, realpt, imagpt, ritzes, trealpt,
c                timagpt, tritzes, tm, nconv, tol, info, workreal,
c                workimag, writzes )
c   Description:
c   This routine takes eigenvalues nu=(lambda-mu)/(lambda-sigma) and
c   transforms them to the lambda plane. It determines how many
c   of the lambdas satisfy the approximation criterion. It picks
c   new values for mu and sigma to be used in the next iteration of pdnaupc
c
c
c   SUBROUTINE ARGUMENTS
c
c    m (integer):        ncv (INPUT)
c    sigma (double):     pole of the Cayley transform (nu) (INPUT/OUTPUT)
c    mu (double):        root of the Cayley transform (nu) (INPUT/OUTPUT)
c    realpt (array(tm)): the real part of the numerically
c                        computed eigenvalue (nu plane) (INPUT/OUTPUT)
c    imagpt (array(tm)): the imaginary part of the numerically
c                        computed eigenvalue (nu)  (INPUT/OUTPUT)
c    ritzes (array(tm)): the error in the computed eigenvalue
c                        (nu) (INPUT/OUTPUT)
c    trealpt (array(tm)):the real part of the transformed eigenvalue in the
c                        lambda plane, sorted by largest real part (OUTPUT)
c    timagpt (array(tm)):the imag part of the transformed eigenvalue
c                       n the lambda plane, sorted by largest real part (OUTPUT)
c    tritzes (array(tm)): the reordered ritz values (OUTPUT)
c    tm (integer):    number of "candidate evals" in lambda plane (OUTPUT)
c    nconv (integer): on input, number of eigenvalues requested,
c                     on output, number of converged lambdas (INPUT/OUTPUT)
c                     Iparam(5) in f77; Iparam(4) in C
c    tol (double) :   convergence criterion for eigenvalues (INPUT)
c    info (integer):  provides error information (OUTPUT)
c         info = 1:  could not achieve criterion for distance between x1, x2
c         info = 0 :  routine has performed normally
c         info = -1 :  the choice of sigma was not to the
c                     right of the largest eigenvalue in the lambda plane
c         info = -2 :  the number of candidate eigenvalues is less than
c                     the number of converged eigenvalues
c         info = -3 : could not find which transformed eigenvalue in the
c                     nu plane corresponds to the eigenvalue of largest real
c                     part in the lambda plane
c         info = -4 : all candidate eigenvalues lambda have the same real part
c         info = -5 : nconv on exit must be <= ncv/2 
c
c    workreal (array(tm)): work array used to hold real part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    workimag (array(tm)): work array used to hold imaginary part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    writzes (array(tm)):  work array to hold ritz values
c
c
c    VARIABLES
c
c    zsigma (complex):    sigma for complex calculations
c    zmu (complex):       mu for complex calculations
c    c3,c4,c5 (complex):  work scalars
c    mx (double):         used in computing a new shift
c    mn (double):         used in computing a new shift
c    x1 (double):         real part of rightmost eval in lambda plane
c    x2 (double):         real part of the [max(tm, inncnv+1)]th
c                         largest magnitude eigenvalue in the nu plane
c    inncnv (integer):    the number of eigenvalues (lambda) we wish to
c                         have converged 
c    jj (integer):        a counter
c
c    CONSTANTS
c
c    one: 1.0
c    zero: 0.0
c    stripw: .15
c
c Note: the above constants are used in transforming the Cayley eigenvalues
c from the nu to the lambda plane. "stripw" is the strip width within which
c we do not consider Cayley transformed eigenvalues (nu), for they
c correspond to infinite eigenvalues (lambda)
c
c  Routines called
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dsortc  ARPACK sorting routine.
c             uses the following:
c             'LM' -> want the NEV eigenvalues of largest magnitude.
c             'SM' -> want the NEV eigenvalues of smallest magnitude.
c             'LR' -> want the NEV eigenvalues of largest real part.
c             'SR' -> want the NEV eigenvalues of smallest real part.
c             'LI' -> want the NEV eigenvalues of largest imaginary part.
c             'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c-----------------------------------------------------------------------
      subroutine polez2 ( m, sigma, mu, delta, realpt, imagpt, 
     &           ritzes, trealpt,
     &           timagpt, tritzes, tm, nconv, tol, info )
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
      integer m, tm, nconv, info
      double precision delta
      double precision realpt(*), imagpt(*), ritzes(*), trealpt(*),
     &     timagpt(*), tritzes(*),  sigma, mu, tol
      complex*16 zsigma, zmu, c3, c4, c5
c     double precision temp, mx, mn, x1, x2, y1
      integer j, inncnv, jj
c
      double precision one, zero, stripw
      parameter (one = 1.0D+0, zero = 0.0D+0, stripw=.15)
c
      double precision dlapy2
      external   dlapy2
c
      intrinsic  dcmplx, dimag, dble, max, min
c
c     Start executable statements
c
      info = 0
      zsigma = dcmplx( sigma, zero )
      zmu = dcmplx( mu, zero )
      tm = 0
      inncnv = nconv
c
c     Eigenvalues nu=(lambda - mu)/(lambda-sigma) are transformed
c     to the lambda plane, and the values of the real and imaginary
c     parts of lambdas are held in trealpt, timagpt.
c     The loop runs through m times, which is the number of eigenvalues
c     requested. It only stores tm values of lambda, requiring
c     that the real part of nu is greater than one.  
c
      do 10 j=1, m
         if (  realpt(j) .lt. 1.d0/delta ) then
            c3 = dcmplx( realpt(j), imagpt(j) )
            c4 = dcmplx( delta*realpt(j)-one, delta*imagpt(j) )
            c5 = (zsigma*c3 - zmu)/c4
            tm = tm + 1
            trealpt(tm) = dble( c5 )
            timagpt(tm) = dimag( c5 )
            tritzes(tm) = ritzes(j)
         endif
 10   continue
c
c     Sort the lambda eigenvalues and the error bounds, store in trealpt,
c     timagpt, tritzes
c     'LR' gives the eigenvalues of largest real part
c
      call dsortc ('LR', .true., tm, trealpt, timagpt, tritzes)
c
c    error if the number of candidate eigenvalues is
c    less than or equal to the number of eigenvalues we want converged
c
      if ( tm .le. nconv ) then
         info = -2
         return
      endif
c
c     Determine how many consecutive lambda eigenvalues
c     satisfy the error tolerance and store the result in nconv.
c     This is regardless of how many eigenvalues were requested
c     (inncnv), but will not return more than m/2+1. 
c
      jj = 0
      flag = 0  
      j=tm
 15   if(flag.eq.1) goto 20
        if (j.le.tm .and. j.ge.1) then
         c3 = dcmplx( trealpt(j)-mu, timagpt(j) )
         c4 = dcmplx( delta*trealpt(j)-sigma, delta*timagpt(j) )
         c5 = c3 / c4
         temp = dlapy2( dble(c5), dimag(c5) )
         if  (tritzes(j) .le. tol*temp) then
           jj = jj + 1
         else
           flag=1
         endif
c
c     ensures loop exit
c
         if (j.eq.1) flag=1
         j= j-1
        endif
        goto 15
 20   continue
c
c     jj is the counter; decides how many evals are within
c     the tolerance and assigns the value to nconv
c
c    nconv now means how many consecutive evals are converged, not
c    how many we wanted. inncnv contains how many evals we wanted.
c
        nconv = jj
c
      sigma  = sigma 
      mu     = mu
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine stslc2 ( m, nconv, tol, sigma, mu, delta,
     &           cutoff, ritzr, 
     &           ritzi, errbds, select )
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
      integer m, nconv
      double precision cutoff, tol, sigma, mu, delta
      double precision ritzr(m), ritzi(m), errbds(m)
      integer select(m)
c
      complex*16 zsigma, zmu, c3, c4, c5
      integer j
c     
      double precision one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0 )
c     
      double precision dlapy2
      external   dlapy2
c     
      intrinsic  dcmplx, dble, abs
c
c     Start executable statements
c
      zsigma = dcmplx( sigma, zero )
      zmu = dcmplx( mu, zero )
      do 10 j=1, m
         c3 = dcmplx( ritzr(j), ritzi(j) )
         c4 = dcmplx( delta*ritzr(j)-one, delta*ritzi(j) )
         c5 = (zsigma*c3 - zmu) / c4
         if ( (dble( c5 ) - cutoff .ge. -1.0d-10*abs(cutoff))
     & .and. ( ritzr(j) .lt. 1.d0/delta)) then
c
c .and.  &        (errbds(j) .le. tol*dlapy2( ritzr(j), ritzi(j) )) ) then
c
            select(j)=1
         else
            select(j)=0
         end if
 10   continue
c
#endif
      return
      end
c
c   Name: polez3 (m, sigma, mu, delta, zeta, realpt, imagpt, ritzes, trealpt,
c                timagpt, tritzes, tm, nconv, tol, info, workreal,
c                workimag, writzes )
c   Description:
c   This routine takes eigenvalues nu=1/(lambda-sigma) and
c   transforms them to the lambda plane. It determines how many
c   of the lambdas satisfy the approximation criterion. It picks
c   new values for mu and sigma to be used in the next iteration of pdnaupc
c   This is the shift-and-invert version of polez2.
c
c
c   SUBROUTINE ARGUMENTS
c
c    m (integer):        ncv (INPUT)
c    sigma (double):     pole of the Cayley transform (nu) (INPUT/OUTPUT)
c    mu (double):        root of the Cayley transform (nu) (INPUT/OUTPUT)
c    realpt (array(tm)): the real part of the numerically
c                        computed eigenvalue (nu plane) (INPUT/OUTPUT)
c    imagpt (array(tm)): the imaginary part of the numerically
c                        computed eigenvalue (nu)  (INPUT/OUTPUT)
c    ritzes (array(tm)): the error in the computed eigenvalue
c                        (nu) (INPUT/OUTPUT)
c    trealpt (array(tm)):the real part of the transformed eigenvalue in the
c                        lambda plane, sorted by largest real part (OUTPUT)
c    timagpt (array(tm)):the imag part of the transformed eigenvalue
c                       n the lambda plane, sorted by largest real part (OUTPUT)
c    tritzes (array(tm)): the reordered ritz values (OUTPUT)
c    tm (integer):    number of "candidate evals" in lambda plane (OUTPUT)
c    nconv (integer): on input, number of eigenvalues requested,
c                     on output, number of converged lambdas (INPUT/OUTPUT)
c                     Iparam(5) in f77; Iparam(4) in C
c    tol (double) :   convergence criterion for eigenvalues (INPUT)
c    info (integer):  provides error information (OUTPUT)
c         info = 1:  could not achieve criterion for distance between x1, x2
c         info = 0 :  routine has performed normally
c         info = -1 :  the choice of sigma was not to the
c                     right of the largest eigenvalue in the lambda plane
c         info = -2 :  the number of candidate eigenvalues is less than
c                     the number of converged eigenvalues
c         info = -3 : could not find which transformed eigenvalue in the
c                     nu plane corresponds to the eigenvalue of largest real
c                     part in the lambda plane
c         info = -4 : all candidate eigenvalues lambda have the same real part
c         info = -5 : nconv on exit must be <= ncv/2 
c
c    workreal (array(tm)): work array used to hold real part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    workimag (array(tm)): work array used to hold imaginary part of
c                      eigenvalues as they are
c                      transformed between lambda and nu plane
c    writzes (array(tm)):  work array to hold ritz values
c
c
c    VARIABLES
c
c    zsigma (complex):    sigma for complex calculations
c    zmu (complex):       mu for complex calculations
c    c3,c4,c5 (complex):  work scalars
c    mx (double):         used in computing a new shift
c    mn (double):         used in computing a new shift
c    x1 (double):         real part of rightmost eval in lambda plane
c    x2 (double):         real part of the [max(tm, inncnv+1)]th
c                         largest magnitude eigenvalue in the nu plane
c    inncnv (integer):    the number of eigenvalues (lambda) we wish to
c                         have converged 
c    jj (integer):        a counter
c
c    CONSTANTS
c
c    one: 1.0
c    zero: 0.0
c    stripw: .15
c
c Note: the above constants are used in transforming the Cayley eigenvalues
c from the nu to the lambda plane. "stripw" is the strip width within which
c we do not consider Cayley transformed eigenvalues (nu), for they
c correspond to infinite eigenvalues (lambda)
c
c  Routines called
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dsortc  ARPACK sorting routine.
c             uses the following:
c             'LM' -> want the NEV eigenvalues of largest magnitude.
c             'SM' -> want the NEV eigenvalues of smallest magnitude.
c             'LR' -> want the NEV eigenvalues of largest real part.
c             'SR' -> want the NEV eigenvalues of smallest real part.
c             'LI' -> want the NEV eigenvalues of largest imaginary part.
c             'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c-----------------------------------------------------------------------
      subroutine polez3 ( m, sigma, mu, delta, zeta, realpt, imagpt, 
     &           ritzes, trealpt,
     &           timagpt, tritzes, tm, nconv, tol, info )
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
      integer m, tm, nconv, info
      double precision delta, zeta
      double precision realpt(*), imagpt(*), ritzes(*), trealpt(*),
     &     timagpt(*), tritzes(*),  sigma, mu, tol
      complex*16 zsigma, c3, c4, c5
      double precision temp
      integer j, inncnv, jj, itf
c
      double precision one, zero, stripw
      parameter (one = 1.0D+0, zero = 0.0D+0, stripw=.15)
c
      double precision dlapy2
      external   dlapy2
c
      intrinsic  dcmplx, dimag, dble, max, min
c
c     Start executable statements
c
      info = 0
      zsigma = dcmplx( sigma, zero )
      tm = 0
      inncnv = nconv
c
c     Eigenvalues nu=1.0/(lambda-sigma) are transformed
c     to the lambda plane, and the values of the real and imaginary
c     parts of lambdas are held in trealpt, timagpt.
c     The loop runs through m times, which is the number of eigenvalues
c     requested. It only stores tm values of lambda, requiring
c     that the magnitude of nu is bigger then some threshhold.
c     The threshhold is in subroutine sitest()
c
      do 10 j=1, m
         call sitest(realpt(j), imagpt(j), zeta, m, itf)
         if ( itf.eq.1 ) then
            c3 = dcmplx( realpt(j), imagpt(j) )
            c5 = zsigma + 1.0/c3
            tm = tm + 1
            trealpt(tm) = dble( c5 )
            timagpt(tm) = dimag( c5 )
            tritzes(tm) = ritzes(j)
         endif
 10   continue
c
c     Sort the lambda eigenvalues and the error bounds, store in trealpt,
c     timagpt, tritzes
c     'LR' gives the eigenvalues of largest real part
c
      call dsortc ('LR', .true., tm, trealpt, timagpt, tritzes)
c
c    error if the number of candidate eigenvalues is
c    less than or equal to the number of eigenvalues we want converged
c
      if ( tm .le. nconv ) then
         info = -2
         return
      endif
c
c     Determine how many consecutive lambda eigenvalues
c     satisfy the error tolerance and store the result in nconv.
c     This is regardless of how many eigenvalues were requested
c     (inncnv), but will not return more than m/2+1. 
c
      jj = 0
      flag = 0  
      j=tm
 15   if(flag.eq.1) goto 20
        if (j.le.tm .and. j.ge.1) then
         c3 =  1.0
         c4 = dcmplx( delta*trealpt(j)-sigma, delta*timagpt(j) )
         c5 = c3 / c4
         temp = dlapy2( dble(c5), dimag(c5) )
         if  (tritzes(j) .le. tol*temp) then
           jj = jj + 1
         else
           flag=1
         endif
c
c     ensures loop exit
c
         if (j.eq.1) flag=1
         j= j-1
        endif
        goto 15
 20   continue
c
c     jj is the counter; decides how many evals are within
c     the tolerance and assigns the value to nconv
c
c    nconv now means how many consecutive evals are converged, not
c    how many we wanted. inncnv contains how many evals we wanted.
c
        nconv = jj
c
      sigma  = sigma 
      mu     = mu
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine stslc3 ( m, nconv, tol, sigma, mu, delta, zeta,
     &           cutoff, ritzr, 
     &           ritzi, errbds, select )
c
c    This is the shift-and-invert version of stslc2
c
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
      integer m, nconv
      double precision cutoff, tol, sigma, mu, delta, zeta
      double precision ritzr(m), ritzi(m), errbds(m)
      integer select(m)
c
      complex*16 zsigma, c3, c5
      integer j, itf
c     
      double precision one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0 )
c     
      double precision dlapy2
      external   dlapy2
c     
      intrinsic  dcmplx, dble, abs
c
c     Start executable statements
c
      zsigma = dcmplx( sigma, zero )
      do 10 j=1, m
         c3 = dcmplx( ritzr(j), ritzi(j) )
         c5 = zsigma + 1.0 / c3
c select rightmost except for extraneous values
         call sitest(ritzr(j), ritzi(j), zeta, m, itf)
         if ( (dble( c5 ) - cutoff .ge. -1.0d-10*abs(cutoff))
     &          .and. itf.eq.1 ) then
c
            select(j)=1
         else
            select(j)=0
         end if
 10   continue
c
#endif
      return
      end
c-----------------------------------------------------------------------
#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)
c
c This function is a heuristic for weeding out the unphysical
c eigenvalue candidates from the shift and invert runs. If
c an eigenvalue in the transformed plane is very far from the
c shift point, and therefore has a very small magnitude in the
c transformed system, then we throw it out. The tolerance is a
c function of the Arnoldi space size, m, because we can trust
c small values better when we have a large arnoldi space.
c This is a first cut heuristic from running one problem, so
c feel free to modify it.
c
      subroutine sitest(tr, ti, zeta, m, itf)
      double precision tr, ti, zeta
      integer m, itf

      if ((tr*tr+ti*ti) .gt. (zeta / (100.0 + m*m))) then
         itf=1
      else 
         itf = 0
      endif

      return
      end
#endif
