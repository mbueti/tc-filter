  FUNCTION MAXMAG(A,LENGTH)
!
!*********************************************************************
!                                                                    *
!   THIS FUNCTION OBTAINS THE INDEX (OR LOCATION)  OF THE MAXIMUM    *
!   OR MINIMUM VALUE OF AN ARRAY                                     *
!                                                                    *
!*********************************************************************
!                                                                    *
!
!       FUNCTION RETURNS ZERO ORIGIN INDEX!!!!!
    DIMENSION A(LENGTH)
    HTEMP=-1.E30
    DO 10 I=1,LENGTH
    IF(A(I).GT.HTEMP) IB=I
    IF(A(I).GT.HTEMP) HTEMP=A(I)
10  CONTINUE
    MAXMAG=IB-1
    RETURN
  ENTRY MINMAG(A,LENGTH)
    HTEMP=1.E30
    DO 11 I=1,LENGTH
      IF(A(I).LT.HTEMP) IB=I
      IF(A(I).LT.HTEMP) HTEMP=A(I)
11  CONTINUE
    MINMAG=IB-1
    RETURN
  ENTRY MINVAL(A,LENGTH)
    HTEMP=1.E30
    DO 20 I=1,LENGTH
      IF(A(I).LT.HTEMP) IB=I
      IF(A(I).LT.HTEMP) HTEMP=A(I)
20  CONTINUE
    MINVAL=IB-1
    RETURN
  ENTRY MAXVAL(A,LENGTH)
    HTEMP=-1.E30
    DO 21 I=1,LENGTH
      IF(A(I).GT.HTEMP) IB=I
      IF(A(I).GT.HTEMP) HTEMP=A(I)
21  CONTINUE
    MAXVAL=IB-1
    RETURN
    END FUNCTION MAXMAG