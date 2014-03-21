/************************************************************************ *
 * Copyright (c) 2013 Sandia Corporation. All rights reserved.            *
 *                                                                        *
 * Under the terms of Contract DE-AC04-94AL85000, there is a              *
 * non-exclusive license for use of this work by or on behalf of the      *
 * U.S. Government.                                                       *
 *                                                                        *
 *                                NOTICE:				  * 
 * For five (5) years from 10/03/2013, the United States Government is 	  *
 * granted for itself and others acting on its behalf a paid-up,    	  *
 * nonexclusive, irrevocable worldwide license in this data to reproduce, *
 * prepare derivative works, and perform publicly and display publicly,   *
 * by or on behalf of the Government. There is provision for the possible *
 * extension of the term of this license. Subsequent to that period or    *
 * any extension granted, the United States Government is granted for     *
 * itself and others acting on its behalf a paid-up, nonexclusive,        *
 * irrevocable worldwide license in this data to reproduce, prepare       *
 * derivative works, distribute copies to the public, perform publicly    *
 * and display publicly, and to permit others to do so. The specific      *
 * term of the license can be identified by inquiry made to               *
 * Sandia Corporation or DOE.                                             *
 *	                                                                  *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT *
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES   *
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY  *
 * FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION,      *
 * APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE   *
 * WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.                             *
 *
 * Any licensee of this software has the obligation and responsibility to *
 * abide by the applicable export control laws, regulations, and general  *
 * prohibitions relating to the export of technical data. Failure to 	  *
 * obtain an export control license or other authority from the Government*
 *  may result in criminal liability under U.S. laws.  	    	          *
 *	 	                                                          *
 *                             (End of Notice)                            *
C  *					         			*
C \************************************************************************/
C  loca_eigen_c2f.F,v
C
C====================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C-----------------------------------------------------------------------------
C  LOCA 1.0: Library of Continuation Algorithms
C  Copyright (C) 2001, Sandia National Laboratories
C-----------------------------------------------------------------------------
C
C
C 
C  The following automates switches between ARPACK and PARPACK calls.
C  Only compiler directives EIGEN_PARALLEL or EIGEN_SERIAL apply.
C  

#if defined (EIGEN_SERIAL) || defined (EIGEN_PARALLEL)

      subroutine cpdnaupc(comm, ido, bmat, nloc, which, nev, tol, 
     &resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )

c
c     Purpose:
c     -------
c
      integer    comm, ido, nloc, info, ldv, lworkl, ncv, nev
      character  bmat, which*2
c      integer    info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
      integer    iparam(7), ipntr(11)
      Double precision
     &           resid(*), v(ldv,ncv), 
     &           workd(*), workl(lworkl)

#ifdef EIGEN_PARALLEL
      call pdnaupc( comm,
#else
      call dnaupc(
#endif
     &ido, bmat, nloc, which, nev, tol, resid,
     &ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
      return
      end

cs
      SUBROUTINE CPDMOUT( COMM, LOUT, M, N, A, LDA, IDIGIT)
      INTEGER            COMM, IDIGIT, LDA, LOUT, M, N
      DOUBLE PRECISION   A( LDA, * )

#ifdef EIGEN_PARALLEL
      call PDMOUT( comm,
#else
      call DMOUT(
#endif
     &             LOUT, M, N, A, LDA, IDIGIT, 
     &            'Ritz values (Real,Imag) and direct residuals')

      return
      end
cs
      subroutine cpdneupc( comm, ivec, howmny, celect, d, v, ldv,
     &   sigma, mu, delta, workev, bmat, n, n2, which, nev, tol, resid,
     &   ncv, iparam, ipntr, workd, workl, lworkl, ierr, select)
c
c                         C to Fortran conversion
c
      logical    select(ncv), rvec 
      character  howmny, bmat, which*2
      integer    comm, ivec, celect(ncv), n, n2, ncv, nev, ldv,
     &           iparam(11), ipntr(14), lworkl, ierr, i
      Double precision
     &           sigma, mu, tol, d(3*ncv), resid(n), v(ldv,ncv),
     &           workd(3*n2), workl(lworkl), workev(3*ncv)
c     ..
c     .. Executable statements ..
      if (ivec .eq. 0) then
         rvec = .FALSE.
      else
         rvec = .TRUE.
      endif
      if (howmny .eq. 'S' .or. iparam(6).eq.1) then
         do 100 i = 1,ncv
            if (celect(i) .eq. 0) then
               select(i) = .FALSE.
            else
               select(i) = .TRUE.
            endif
 100     continue
      elseif (howmny .ne. 'A') then
         print *,'Error in c_pdneupd'
         print *,'unknown value of howmny ',howmny
         stop
      endif

#ifdef EIGEN_PARALLEL
      call pdneupc ( comm,
#else
      call dneupc (
#endif
     &        rvec, 'A', select, d, d(ncv+1), v, ldv,
     &        sigma, mu, delta, workev, bmat, n, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
      return
      end
#else
      subroutine empty()
        i=0
      end
#endif
