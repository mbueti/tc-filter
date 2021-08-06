  SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
!***BEGIN PROLOGUE  XERRWV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes error message allowing 2 integer and two real
!            values to be included in the message.
!***DESCRIPTION
!     Abstract
!        XERRWV processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!        In addition, up to two integer values and two real
!        values may be printed along with the message.
!
!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed.
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
!        NI    - number of integer values to be printed. (0 to 2)
!        I1    - first integer value.
!        I2    - second integer value.
!        NR    - number of real values to be printed. (0 to 2)
!        R1    - first real value.
!        R2    - second real value.
!
!     Examples
!        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
!    1   1,NUM,0,0,0.,0.)
!        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
!    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
!
!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
!                    XGETUA
!***END PROLOGUE  XERRWV
    CHARACTER*(*) MESSG
    CHARACTER*20 LFIRST
    CHARACTER*37 FORM
    DIMENSION LUN(5)
!   GET FLAGS
!***FIRST EXECUTABLE STATEMENT  XERRWV
    LKNTRL = J4SAVE(2,0,.FALSE.)
    MAXMES = J4SAVE(4,0,.FALSE.)
!   CHECK FOR VALID INPUT
    IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND. &
        (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
    IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
    CALL XERPRT('XERROR -- INVALID INPUT',23)
    IF (LKNTRL.GT.0) CALL FDUMP
    IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',29)
    IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
    CALL XERABT('XERROR -- INVALID INPUT',23)
    RETURN
10  CONTINUE
!   RECORD MESSAGE
    JUNK = J4SAVE(1,NERR,.TRUE.)
    CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
!   LET USER OVERRIDE
    LFIRST = MESSG
    LMESSG = NMESSG
    LERR = NERR
    LLEVEL = LEVEL
    CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
!   RESET TO ORIGINAL VALUES
    LMESSG = NMESSG
    LERR = NERR
    LLEVEL = LEVEL
    LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
    MKNTRL = IABS(LKNTRL)
!   DECIDE WHETHER TO PRINT MESSAGE
    IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
    IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES))) &
        .OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES)) &
        .OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1)) &
        .OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
    IF (LKNTRL.LE.0) GO TO 20
    CALL XERPRT(' ',1)
!   INTRODUCTION
    IF (LLEVEL.EQ.(-1)) CALL XERPRT('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
    IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
    IF (LLEVEL.EQ.1) CALL XERPRT('RECOVERABLE ERROR IN...',23)
    IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
20  CONTINUE
!   MESSAGE
    CALL XERPRT(MESSG,LMESSG)
    CALL XGETUA(LUN,NUNIT)
    ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
    ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
    DO 50 KUNIT=1,NUNIT
      IUNIT = LUN(KUNIT)
      IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
      DO 22 I=1,MIN(NI,2)
        WRITE (FORM,21) I,ISIZEI
21      FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
        IF (I.EQ.1) WRITE (IUNIT,FORM) I1
        IF (I.EQ.2) WRITE (IUNIT,FORM) I2
22    CONTINUE
      DO 24 I=1,MIN(NR,2)
        WRITE (FORM,23) I,ISIZEF+10,ISIZEF
23      FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',I2,'.',I2,')')
        IF (I.EQ.1) WRITE (IUNIT,FORM) R1
        IF (I.EQ.2) WRITE (IUNIT,FORM) R2
24    CONTINUE
      IF (LKNTRL.LE.0) GO TO 40
!     ERROR NUMBER
      WRITE (IUNIT,30) LERR
30    FORMAT (15H ERROR NUMBER =,I10)
40    CONTINUE
50  CONTINUE
!   TRACE-BACK
    IF (LKNTRL.GT.0) CALL FDUMP
100 CONTINUE
    IFATAL = 0
    IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2))) IFATAL = 1
!   QUIT HERE IF MESSAGE IS NOT FATAL
    IF (IFATAL.LE.0) RETURN
    IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
!   PRINT REASON FOR ABORT
    IF (LLEVEL.EQ.1) CALL XERPRT('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
    IF (LLEVEL.EQ.2) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',29)
!   PRINT ERROR SUMMARY
    CALL XERSAV(' ',-1,0,0,KDUMMY)
120 CONTINUE
!   ABORT
    IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
    CALL XERABT(MESSG,LMESSG)
    RETURN
  END SUBROUTINE XERRWV