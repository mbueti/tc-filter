  SUBROUTINE XERABT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERABT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Aborts program execution and prints error message.
!***DESCRIPTION
!     Abstract
!        ***Note*** machine dependent routine
!        XERABT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG and NMESSG are as in XERROR, except that NMESSG may
!        be zero, in which case no message is being supplied.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERABT
    CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERABT
      ! STOP
  END SUBROUTINE XERABT