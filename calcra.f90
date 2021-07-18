  SUBROUTINE CALCRa(RO,RTAN,iang)
    ! PARAMETER ( NMX=64)
    ! PARAMETER (IMX=640 , JMX=320)
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
!
    COMMON  /TOTAL/ DDEL,DTHA
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
!
    ! PI = 4.*ATAN(1.0)
    ! PI180 = 4.*ATAN(1.0)/180.
    FACT =  COS(YOLD)
!
    DX=DDEL/PI180
    DY=DTHA/PI180
    XC = (XOLD-XCORN)*DX
    YC = (YOLD-YCORN)*DY

!      
    THETA= 2.*PI*FLOAT(iang-1)/FLOAT(NMX)
      !   X=RO*fact*COS(THETA)+XC
    X=RO*COS(THETA)+XC
    Y=RO*SIN(THETA)+YC
    IX=FLOOR(X/DX)+1
    IY=FLOOR(Y/DY)+1
      !   print *, 'RO=', RO
      !   print *, 'theta=', THETA
      !   print *, 'XC=', XC
      !   print *, 'YC=', YC
      !   print *, 'Y=',Y
    IX1=IX+1
    IY1=IY+1
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
      !   print *, 'RO=',RO
      !   print *, 'XC=',XC
      !   print *, 'X=',X
      !   print *, 'DX=',DX
      !   print *, "IX=",IX
      !   print *, 'YC=',YC
      !   print *, "Y=",Y
      !   print *, "DY=",DY
      !   print *, "IY=",IY
      !   print *, 'mody=',modulo(IY+1,IY)
      !   print *, XF(IX+1,IY)
      !   print *, XF(IX+1,IY+1)
!        rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,MODULO(IY+1, IY))
!      1      +  (1.-Q)*P*XF(max(MODULO(IX+1, IX),1),IY) +
!      1      P*Q*XF(max(MODULO(IX+1,IX),1),max(MODULO(IY+1,IY),1))
      !   print *, XF(max(MODULO(IX+1,IX),1),max(MODULO(IY+1,IY),1))
    rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY+1) +  (1.-Q)*P*XF(IX+1,IY)+P*Q*XF(IX+1,IY+1)
! 10     CONTINUE
!
!
    RETURN
  END SUBROUTINE CALCRa