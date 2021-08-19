/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2021 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.	  *
\************************************************************************/

C user_senkin.F,v

C====================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C========================================================================
C
#ifdef SENKIN_OUTPUT
*
      integer function senkin_open()
*
*       Fortran function to open save.dat
* 
      open (10, file='goma_senkin_save.dat', status = 'UNKNOWN',
     $      form = 'UNFORMATTED', err = 20)

      senkin_open = 10
      return
20    continue
      senkin_open = 0
      end
C
C========================================================================
C
      SUBROUTINE write_header (LOUT, LSAVE, LSENS, NSYS, KK, II, ICASE,
     $                         NEQ, TSTOP, TLIM, TOLS, ICKWRK, RCKWRK, 
     $                         CCKWRK)
*
* Parameters:
* 
      INTEGER                LOUT, LSAVE, NSYS, KK, II, ICASE, NEQ,
     $                       ICKWRK(*)
      LOGICAL                LSENS
      CHARACTER              CCKWRK(*)*(*)
*
C*****precision > double
      DOUBLE PRECISION       TSTOP, TLIM, TOLS(4), RCKWRK(*)
C*****END precision > double
*
C*****precision > single
C      REAL                  TSTOP, TLIM, TOLS(4), RCKWRK(*)
C*****END precision > single
* ----------------------------------------------------------------------
C
C write_header: This subroutine writes the header information onto the
C               binary output file
C               bin_ver is the version number of the binary output file.
C               It should be incremented whenever the format of the
C               output file has changed.
C
C Description of output:
C
C     NSYS  = Number of equations to be solved (depends on ICASE)
C     KK    = Number of gas-phase species
C     II    = Number of gas-phas reactions
C     ICASE = 1 to 5 - type of problem being solved:
C                   1 == CONSTANT PRESSURE
C                   2 == CONSTANT VOLUME
C                   3 == VOLUME GIVEN AS A FUNCTION OF TIME
C                   4 == CONSTANT TEMPERATURE, PRESSURE
C                   5 == TEMPERATURE GIVEN AS A FUNCTION OF TIME
C     NEQ   = Number of equations solved by DASAC.
C     TSTOP = Integration out to this value
C     TLIM  = Definition of the temperature at which the ignition 
C             delay will be based on
C     TOLS(1) = RTOL = Absolute tolerance for the solution components
C     TOLS(2) = ATOL = Relative tolerance for the solution components
C     TOLS(3) = ATLS = Absolute tolerance for the sensitivity coefficients
C     TOLS(4) = RTLS = Relative tolerance for the sensitivity coefficients
C
C   -----------------------------------------------------------------------

      INTEGER            bin_ver
      PARAMETER         (bin_ver = 1)
      CHARACTER          SECTION_HEAD*16
      EXTERNAL           CKSAVE

      IF (LSAVE .EQ. 0) THEN
        WRITE (LOUT,*)'LSAVE = 0 - something is wrong'
        STOP 2
      END IF

      SECTION_HEAD = '                '
      SECTION_HEAD = 'SENKIN bin file'
      WRITE (LSAVE) SECTION_HEAD
      WRITE (LSAVE) bin_ver
      SECTION_HEAD = 'SENKIN ProbDef'
      WRITE (LSAVE) SECTION_HEAD
      WRITE (LSAVE) LSENS
      WRITE (LSAVE) NSYS, KK, II, ICASE, NEQ
      WRITE (LSAVE) TSTOP, TLIM, TOLS(1), TOLS(2), TOLS(3), TOLS(4)
      SECTION_HEAD = 'CKLINK'
      WRITE (LSAVE) SECTION_HEAD
*
*     Call the chemkin routine that write out the common blocks and
*     work arrays - This routine creates 4 records of information
*
      CALL CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
*
* Write the start of the data section
*
      SECTION_HEAD = 'SENKIN SOLUTION'
      WRITE (LSAVE) SECTION_HEAD
*
      RETURN
      END
C
C========================================================================
C
      SUBROUTINE write_record (LSAVE, ICASE, TIM, Y, NSYS, LSENS, II,
     $                         ICKWRK, RCKWRK, T, P, VOL)
*
      INTEGER                LSAVE, ICASE, NSYS, II, ICKWRK(*)
      LOGICAL                LSENS
C*****precision > double
      DOUBLE PRECISION       TIM, Y(NSYS), RCKWRK(*),
     $                       T, P, VOL
C*****END precision > double
*
C*****precision > single
C      REAL                  TIM, Y(NSYS), RCKWRK(*),
C     $                       T, P, VOL
C*****END precision > single
*
* write_record: This subroutine writes a single output record out to
*               the binary output file.
*
*   Variables:
*
*        TIM  = Current time of the record
*        P    = Current value of the system pressure (dynes cm-2)
* -------------------------------------------------------------------
*
* Local Variables:
*
      INTEGER                  I
C*****precision > double
      DOUBLE PRECISION         RHO
C*****END precision > double
*
C*****precision > single
C      REAL                    RHO
C*****END precision > single
*
* --------------------------------------------------------------------
*
* Find the correct density:
*
      CALL CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
*
* Write out the current value of the solution =
*  - Note the record will have the same length and form
*    no matter what ICASE is. Temperature will be the 4th 
*    entry, and the mass fractions will be after.
*
 
      WRITE (LSAVE) TIM, P, VOL, RHO, T, (Y(I), I = 1, NSYS-1)
*
* Write out the sensitivity coefficients -> option not implemented
*
      IF (LSENS) THEN
        STOP 1
*       WRITE (LSAVE) (( Z(I,J), I=1,NSYS), J = 2, II+1)
      END IF
*
      RETURN
      END
C
#else
*
      integer function senkin_open()
*
*       dummy Fortran function
* 
      senkin_open = 0
      end
C======================================================================
#endif
C
