SUBROUTINE WNNLS(W,MDW,ME,MA,N,L,PRGOPT,X,RNORM,MODE,IWORK,WORK)
!***BEGIN PROLOGUE  WNNLS
!***DATE WRITTEN   790701   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  K1A2A
!***KEYWORDS  CONSTRAINED LEAST SQUARES,CURVE FITTING,DATA FITTING,
!             EQUALITY CONSTRAINTS,INEQUALITY CONSTRAINTS,
!             NONNEGATIVITY CONSTRAINTS,QUADRATIC PROGRAMMING
!***AUTHOR  HANSON, R. J., (SNLA)
!           HASKELL, K. H., (SNLA)
!***PURPOSE  Solve a linearly constrained least squares problem with
!            equality constraints and nonnegativity constraints on
!            selected variables.
!***DESCRIPTION
!
!     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),IWORK(M+N),WORK(M+5*N)
!
!     Written by Karen H. Haskell, Sandia Laboratories,
!     and R.J. Hanson, Sandia Laboratories.
!
!     Abstract
!
!     This subfoo solves a linearly constrained least squares
!     problem.  Suppose there are given matrices E and A of
!     respective dimensions ME by N and MA by N, and vectors F
!     and B of respective lengths ME and MA.  This subroutine
!     solves the problem
!
!               EX = F, (equations to be exactly satisfied)
!
!               AX = B, (equations to be approximately satisfied,
!                        in the least squares sense)
!
!               subject to components L+1,...,N nonnegative
!
!     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
!
!     The problem is reposed as problem WNNLS
!
!               (WT*E)X = (WT*F)
!               (   A)    (   B), (least squares)
!               subject to components L+1,...,N nonnegative.
!
!     The subfoo chooses the heavy weight (or penalty parameter) WT.
!
!     The parameters for WNNLS are
!
!     INPUT..
!
!     W(*,*),MDW,  The array W(*,*) is double subscripted with first
!     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
!                  discussion let us call M = ME + MA.  Then MDW
!                  must satisfy MDW.GE.M.  The condition MDW.LT.M
!                  is an error.
!
!                  The array W(*,*) contains the matrices and vectors
!
!                       (E  F)
!                       (A  B)
!
!                  in rows and columns 1,...,M and 1,...,N+1
!                  respectively.  Columns 1,...,L correspond to
!                  unconstrained variables X(1),...,X(L).  The
!                  remaining variables are constrained to be
!                  nonnegative. The condition L.LT.0 or L.GT.N is
!                  an error.
!
!     PRGOPT(*)    This real-valued array is the option vector.
!                  If the user is satisfied with the nominal
!                  subfoo features set
!
!                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!                  Otherwise PRGOPT(*) is a linked list consisting of
!                  groups of data of the following form
!
!                  LINK
!                  KEY
!                  DATA SET
!
!                  The parameters LINK and KEY are each one word.
!                  The DATA SET can be comprised of several words.
!                  The number of items depends on the value of KEY.
!                  The value of LINK points to the first
!                  entry of the next group of data within
!                  PRGOPT(*).  The exception is when there are
!                  no more options to change.  In that
!                  case LINK=1 and the values KEY and DATA SET
!                  are not referenced. The general layout of
!                  PRGOPT(*) is as follows.
!
!               ...PRGOPT(1)=LINK1 (link to first entry of next group)
!               .  PRGOPT(2)=KEY1 (key to the option change)
!               .  PRGOPT(3)=DATA VALUE (data value for this change)
!               .       .
!               .       .
!               .       .
!               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
!               .                       next group)
!               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
!               .  PRGOPT(LINK1+2)=DATA VALUE
!               ...     .
!               .       .
!               .       .
!               ...PRGOPT(LINK)=1 (no more options to change)
!
!                  Values of LINK that are nonpositive are errors.
!                  A value of LINK.GT.NLINK=100000 is also an error.
!                  This helps prevent using invalid but positive
!                  values of LINK that will probably extend
!                  beyond the foo limits of PRGOPT(*).
!                  Unrecognized values of KEY are ignored.  The
!                  order of the options is arbitrary and any number
!                  of options can be changed with the following
!                  restriction.  To prevent cycling in the
!                  processing of the option array a count of the
!                  number of options changed is maintained.
!                  Whenever this count exceeds NOPT=1000 an error
!                  message is printed and the subfoo returns.
!
!                  OPTIONS..
!
!                  KEY=6
!                         Scale the nonzero columns of the
!                  entire data matrix
!                  (E)
!                  (A)
!                  to have length one. The DATA SET for
!                  this option is a single value.  It must
!                  be nonzero if unit length column scaling is
!                  desired.
!
!                  KEY=7
!                         Scale columns of the entire data matrix
!                  (E)
!                  (A)
!                  with a user-provided diagonal matrix.
!                  The DATA SET for this option consists
!                  of the N diagonal scaling factors, one for
!                  each matrix column.
!
!                  KEY=8
!                         Change the rank determination tolerance from
!                  the nominal value of SQRT(SRELPR).  This quantity
!                  can be no smaller than SRELPR, The arithmetic-
!                  storage precision.  The quantity used
!                  here is internally restricted to be at
!                  least SRELPR.  The DATA SET for this option
!                  is the new tolerance.
!
!                  KEY=9
!                         Change the blow-up parameter from the
!                  nominal value of SQRT(SRELPR).  The reciprocal of
!                  this parameter is used in rejecting solution
!                  components as too large when a variable is
!                  first brought into the active set.  Too large
!                  means that the proposed component times the
!                  reciprocal of the parameter is not less than
!                  the ratio of the norms of the right-side
!                  vector and the data matrix.
!                  This parameter can be no smaller than SRELPR,
!                  the arithmetic-storage precision.
!
!                  For example, suppose we want to provide
!                  a diagonal matrix to scale the problem
!                  matrix and change the tolerance used for
!                  determining linear dependence of dropped col
!                  vectors.  For these options the dimensions of
!                  PRGOPT(*) must be at least N+6.  The FORTRAN
!                  statements defining these options would
!                  be as follows.
!
!                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
!                  PRGOPT(2)=7 (user-provided scaling key)
!
!                  CALL SCOPY(N,D,1,PRGOPT(3),1) (copy the N
!                  scaling factors from a user array called D(*)
!                  into PRGOPT(3)-PRGOPT(N+2))
!
!                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
!                  PRGOPT(N+4)=8 (linear dependence tolerance key)
!                  PRGOPT(N+5)=... (new value of the tolerance)
!
!                  PRGOPT(N+6)=1 (no more options to change)
!
!
!     IWORK(1),    The amounts of working storage actually allocated
!     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
!                  respectively.  These quantities are compared with
!                  the actual amounts of storage needed for WNNLS( ).
!                  Insufficient storage allocated for either WORK(*)
!                  or IWORK(*) is considered an error.  This feature
!                  was included in WNNLS( ) because miscalculating
!                  the storage formulas for WORK(*) and IWORK(*)
!                  might very well lead to subtle and hard-to-find
!                  execution errors.
!
!                  The length of WORK(*) must be at least
!
!                  LW = ME+MA+5*N
!                  This test will not be made if IWORK(1).LE.0.
!
!                  The length of IWORK(*) must be at least
!
!                  LIW = ME+MA+N
!                  This test will not be made if IWORK(2).LE.0.
!
!     OUTPUT..
!
!     X(*)         An array dimensioned at least N, which will
!                  contain the N components of the solution vector
!                  on output.
!
!     RNORM        The residual norm of the solution.  The value of
!                  RNORM contains the residual vector length of the
!                  equality constraints and least squares equations.
!
!     MODE         The value of MODE indicates the success or failure
!                  of the subfoo.
!
!                  MODE = 0  Subfoo completed successfully.
!
!                       = 1  Max. number of iterations (equal to
!                            3*(N-L)) exceeded. Nearly all problems
!                            should complete in fewer than this
!                            number of iterations. An approximate
!                            solution and its corresponding residual
!                            vector length are in X(*) and RNORM.
!
!                       = 2  Usage error occurred.  The offending
!                            condition is noted with the error
!                            processing subfoo, XERROR( ).
!
!     User-designated
!     Working arrays..
!
!     WORK(*)      A real-valued working array of length at least
!                  M + 5*N.
!
!     IWORK(*)     An integer-valued working array of length at least
!                  M+N.
!***REFERENCES  K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
!                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
!                 NONNEGATIVITY CONSTRAINTS*, SAND77-0552, JUNE 1978.
!               K.H. HASKELL AND R.J. HANSON, *SELECTED ALGORITHMS FOR
!                 THE LINEARLY CONSTRAINED LEAST SQUARES PROBLEM--
!                 A USERS GUIDE*, SAND78-1290, AUGUST 1979.
!               K.H. HASKELL AND R.H. HANSON, *AN ALGORITHM FOR
!                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
!                 NONNEGATIVITY CONSTRAINTS*, MATH. PROG. 21 (1981),
!                 PP. 98-118.
!               R.J. HANSON AND K.H. HASKELL, *TWO ALGORITHMS FOR THE
!                 LINEARLY CONSTRAINED LEAST SQUARES PROBLEM*, ACM
!                 TRANS. ON MATH. SOFTWARE, SEPT. 1982.
!***ROUTINES CALLED  WNLSM,XERROR,XERRWV
!***END PROLOGUE  WNNLS
!
!     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
!     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
!     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
!     (START AT LINE WITH C++ IN COLS. 1-3.)
!     /REAL (12 BLANKS)/DOUBLE PRECISION/,/, DUMMY/,SNGL(DUMMY)/
!
!     WRITTEN BY KAREN H. HASKELL, SANDIA LABORATORIES,
!     AND R.J. HANSON, SANDIA LABORATORIES.
!     REVISED FEB.25, 1982.
!
!     SUBROUTINES CALLED BY WNNLS( )
!
!++
!     WNLSM         COMPANION SUBROUTINE TO WNNLS( ), WHERE
!                   MOST OF THE COMPUTATION TAKES PLACE.
!
!     XERROR,XERRWV FROM SLATEC ERROR PROCESSING PACKAGE.
!                   THIS IS DOCUMENTED IN SANDIA TECH. REPT.,
!                   SAND78-1189.
!
!     REFERENCES
!
!     1. SOLVING LEAST SQUARES PROBLEMS, BY C.L. LAWSON
!        AND R.J. HANSON.  PRENTICE-HALL, INC. (1974).
!
!     2. BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE, BY
!        C.L. LAWSON, R.J. HANSON, D.R. KINCAID, AND F.T. KROGH.
!        TOMS, V. 5, NO. 3, P. 308.  ALSO AVAILABLE AS
!        SANDIA TECHNICAL REPORT NO. SAND77-0898.
!
!     3. AN ALGORITHM FOR LINEAR LEAST SQUARES WITH EQUALITY
!        AND NONNEGATIVITY CONSTRAINTS, BY K.H. HASKELL AND
!        R.J. HANSON.  AVAILABLE AS SANDIA TECHNICAL REPORT NO.
!        SAND77-0552, AND MATH. PROGRAMMING, VOL. 21, (1981), P. 98-118.
!
!     4. SLATEC COMMON MATH. LIBRARY ERROR HANDLING
!        PACKAGE.  BY R. E. JONES.  AVAILABLE AS SANDIA
!        TECHNICAL REPORT SAND78-1189.
!
  implicit NONE

  integer, intent(in) :: mdw, me, ma, n, l
  integer, intent(out) :: mode
  REAL, dimension(me+ma+5*n), intent(out) :: WORK
  INTEGER, dimension(N+ME+MA), intent(out) :: iwork
  real, dimension(mdw,1), intent(in) :: W
  real, dimension(1), intent(in) :: PRGOPT, x
  real, intent(in) :: rnorm
  integer :: iopt, l1, l2, l3, l4, l5, liw, lw, nerr
  REAL :: DUMMY
!
!
!***FIRST EXECUTABLE STATEMENT  WNNLS
  MODE = 0
  iwork(1)=mdw*6
  iwork(2)=mdw*2
  IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
  IF (.NOT.(IWORK(1).GT.0)) GO TO 20
  LW = ME + MA + 5*N
  IF (.NOT.(IWORK(1).LT.LW)) GO TO 10
  NERR = 2
  IOPT = 1
  print*,'work array',iwork(1),lw
  CALL XERRWV('WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR WORK(*), NEED LW=I1 BELOW', &
              70, NERR, IOPT, 1, LW, 0, 0, DUMMY, DUMMY)
  MODE = 2
  RETURN
10 CONTINUE
20 IF (.NOT.(IWORK(2).GT.0)) GO TO 40
  LIW = ME + MA + N
  IF (.NOT.(IWORK(2).LT.LIW)) GO TO 30
  NERR = 2
  IOPT = 1
  CALL XERRWV('WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR IWORK(*), NEED LIW=I1 BELOW', &
              72, NERR, IOPT, 1, LIW, 0, 0, DUMMY, DUMMY)
  MODE = 2
  RETURN
30 CONTINUE
40 IF (.NOT.(MDW.LT.ME+MA)) GO TO 50
  NERR = 1
  IOPT = 1
  CALL XERROR( 'WNNLS( ), THE VALUE MDW.LT.ME+MA IS AN ERROR', 44, NERR, IOPT)
  MODE = 2
  RETURN
50 IF (0.LE.L .AND. L.LE.N) GO TO 60
  NERR = 2
  IOPT = 1
  CALL XERROR( 'WNNLS( ), L.LE.0.AND.L.LE.N IS REQUIRED', 39, NERR, IOPT)
  MODE = 2
  RETURN
!
!     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
!     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
!     REQUIRED BY THE MAIN SUBROUTINE WNLSM( ).
!
60 L1 = N + 1
  L2 = L1 + N
  L3 = L2 + ME + MA
  L4 = L3 + N
  L5 = L4 + N
!
  CALL WNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK, IWORK(L1), WORK(1), &
             WORK(L1), WORK(L2), WORK(L3), WORK(L4),WORK(L5))
  RETURN
END subroutine wnnls
