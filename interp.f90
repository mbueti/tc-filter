  SUBROUTINE INTERP(UO,VO,SIGP,ILEV,UN,VN)
!
    PARAMETER ( KMAX = 18, ILT=20 )
!
    DIMENSION  UO(ILT),VO(ILT)
    DIMENSION  SIGP(ILT),Q(KMAX)
    DIMENSION  UN(KMAX),VN(KMAX)
    DIMENSION  UF(KMAX),VF(KMAX)
!
    DATA Q/ .9949968,.9814907,.9604809,.9204018,.8563145,.7772229, &
            .6881255,.5935378,.4974484,.4248250,.3748014,.3247711, &
            .2747291,.2246687,.1745733,.1244004,.0739862,.0207469/ 
!
    TMASK = 1.0e20
!
    ILEVM = ILEV -1
    IF(ILEV.LT.2)THEN
      UN = TMASK
      VN = TMASK
      RETURN
    ENDIF
    DO 10 KK = 1 , ILEVM
!
      SGO1 = SIGP(KK)
      SG02 = SIGP(KK+1)
!
      DO 20 K  = 1 ,  KMAX
!
        SGN = Q(K)
        IF(SGN.GT.SIGP(1))THEN
          UF(K) = TMASK
          VF(K) = TMASK
        ENDIF
        IF(SGN.LT.SIGP(ILEV))THEN
          UF(K) = TMASK
          VF(K) = TMASK
        ENDIF
!
        IF(SGN.LE.SGO1.AND.SGN.GE.SG02)THEN
          DX  = SGO1  - SG02
          DX1 = SGO1  - SGN
          DX2 = SGN  - SGO2
          X1 = DX1/DX
          X2 = DX2/DX
          UF(K) = (1.-X1)*UO(KK) + X1*UO(KK+1)
          VF(K) = (1.-X1)*VO(KK) + X1*VO(KK+1)
        ENDIF
!
20    CONTINUE
10  CONTINUE
    DO 30 K = 1 , KMAX
      KK = KMAX + 1 - K
      UN(K) = UF(KK)
      VN(K) = VF(KK)
30  CONTINUE
    RETURN
  END SUBROUTINE interp