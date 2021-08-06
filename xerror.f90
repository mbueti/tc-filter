  SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
!***BEGIN PROLOGUE  XERROR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes an error (diagnostic) message.
!***DESCRIPTION
!     Abstract
!        XERROR processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed, containing
!                no more than 72 characters.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!
!     Examples
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
!                    43,2,1)
!        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
!    1ULLY COLLAPSED.',65,3,0)
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  XERRWV
!***END PROLOGUE  XERROR
    CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERROR
    CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
    RETURN
  END SUBROUTINE XERROR