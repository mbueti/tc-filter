
  SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
!***BEGIN PROLOGUE  XERCTL
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Allows user control over handling of individual errors.
!***DESCRIPTION
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCTL.
!        If the user has provided his own version of XERCTL, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        MESSG1 - the first word (only) of the error message.
!        NMESSG - same as in the call to XERROR or XERRWV.
!        NERR   - same as in the call to XERROR or XERRWV.
!        LEVEL  - same as in the call to XERROR or XERRWV.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERCTL
    IMPLICIT NONE
    CHARACTER*20, INTENT(IN) :: MESSG1
    INTEGER, INTENT(IN) :: NMESSG, NERR, LEVEL
    INTEGER, INTENT(INOUT) :: KONTRL
!***FIRST EXECUTABLE STATEMENT  XERCTL
    RETURN
  END SUBROUTINE XERCTL