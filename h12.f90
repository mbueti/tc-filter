  SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
!***BEGIN PROLOGUE  H12
!***REFER TO  HFTI,LSEI,WNNLS
!
!     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
!
!     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
!     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
!
!     Modified at SANDIA LABS, May 1977, to --
!
!     1)  Remove double precision accumulation, and
!     2)  Include usage of the Basic Linear Algebra Package for
!         vectors longer than a particular threshold.
!
!     Construction and/or application of a single
!     Householder transformation..     Q = I + U*(U**T)/B
!
!     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
!     LPIVOT is the index of the pivot element.
!     L1,M   If L1 .LE. M   the transformation will be constructed to
!            zero elements indexed from L1 through M.   If L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
!                   IUE is the storage increment between elements.
!                                       On exit from H1 U() and UP
!                   contain quantities defining the vector U of the
!                   Householder transformation.   On entry to H2 U()
!                   and UP should contain quantities previously computed
!                   by H1.  These will not be modified by H2.
!     C()    On entry to H1 or H2 C() contains a matrix which will be
!            regarded as a set of vectors to which the Householder
!            transformation is to be applied.  On exit C() contains the
!            set of transformed vectors.
!     ICE    Storage increment between elements of vectors in C().
!     ICV    Storage increment between vectors in C().
!     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
!            no operations will be done on C().
!***ROUTINES CALLED  SAXPY,SDOT,SSWAP
!***END PROLOGUE  H12
    DIMENSION U(IUE,M), C(1)
!***FIRST EXECUTABLE STATEMENT  H12
    ONE=1.
!
    IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
    CL=ABS(U(1,LPIVOT))
    IF (MODE.EQ.2) GO TO 60
!   ****** CONSTRUCT THE TRANSFORMATION. ******
    DO 10 J=L1,M
      CL=AMAX1(ABS(U(1,J)),CL)
10  CONTINUE
    IF (CL.le.0) goto 130
20  CLINV=ONE/CL
    SM=(U(1,LPIVOT)*CLINV)**2
    DO 30 J=L1,M
      SM=SM+(U(1,J)*CLINV)**2
30  CONTINUE
    CL=CL*SQRT(SM)
    IF (U(1,LPIVOT).le.0) goto 50
40   CL=-CL
50   UP=U(1,LPIVOT)-CL
     U(1,LPIVOT)=CL
     GO TO 70
!    ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
!
60  IF (CL.le.0) goto 130
70  IF (NCV.LE.0) RETURN
    B=UP*U(1,LPIVOT)
!                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
!
    IF (B.ge.0) goto 130
80  B=ONE/B
    MML1P2=M-L1+2
    IF (MML1P2.GT.20) GO TO 140
    I2=1-ICV+ICE*(LPIVOT-1)
    INCR=ICE*(L1-LPIVOT)
    DO 120 J=1,NCV
      I2=I2+ICV
      I3=I2+INCR
      I4=I3
      SM=C(I2)*UP
      DO 90 I=L1,M
        SM=SM+C(I3)*U(1,I)
        I3=I3+ICE
90    CONTINUE
      IF (SM.eq.0) GOTO 120
100   SM=SM*B
      C(I2)=C(I2)+SM*UP
      DO 110 I=L1,M
        C(I4)=C(I4)+SM*U(1,I)
        I4=I4+ICE
110   CONTINUE
120 CONTINUE
130 RETURN
140 CONTINUE
    L1M1=L1-1
    KL1=1+(L1M1-1)*ICE
    KL2=KL1
    KLP=1+(LPIVOT-1)*ICE
    UL1M1=U(1,L1M1)
    U(1,L1M1)=UP
    IF (LPIVOT.EQ.L1M1) GO TO 150
    CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
150 CONTINUE
    DO 160 J=1,NCV
      SM=SDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
      SM=SM*B
      CALL SAXPY (MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
      KL1=KL1+ICV
160 CONTINUE
    U(1,L1M1)=UL1M1
    IF (LPIVOT.EQ.L1M1) RETURN
    KL1=KL2
    CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
    RETURN
  END SUBROUTINE H12