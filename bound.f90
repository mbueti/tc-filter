  SUBROUTINE BOUND(NMX,XR,ro)
! 
    ! PARAMETER (IMX=640 , JMX=320)
! 
    DIMENSION XR(NMX),ro(nmx)
!       COMMON  /IFACT/NNN,RO(nmx),RB,IENV
    COMMON  /XXX/  XF(IMX,JMX),XC,YC,DX,DY
    COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
!       COMMON  /COOR/ XV,YV,XOLD,YOLD
!       COMMON  /POSIT/ XOLD,YOLD
    ! PI = 4.*ATAN(1.0)
    fact=cos(yold*pi/180.)
    DO 10 I=1,NMX
      THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
      X=RO(i)/fact*COS(THETA)+XC +1.
      Y=RO(i)*SIN(THETA)+YC +1.
      IX=INT(X/DX)
      IY=INT(Y/DY)
      ! print *,'IX=',IX
      ! print *,'IY=',IY
      if(IX.le.0) IX=IX+IMX
      if(IX.gt.IMX) IX=MOD(IX,IMX)+1
      if(IY.le.0) IY=IY+JMX
      if(IY.gt.JMX) IY=MOD(IY,JMX)+1
      ! print *,'IX1=',IX1
      ! print *,'IY1=',IY1

      IX1=IX+1
      IY1=IY+1
      
      if(IX1.gt.IMX) IX1=MOD(IX1,IMX)+1
      if(IY1.gt.JMX) IY1=MOD(IY1,JMX)+1
      !print *,'IX=',IX
      !print *,'IY=',IY
      !print *,'IX1=',IX1
      ! print *,'IY1=',IY1
      P=X/DX-FLOAT(IX)
      Q=Y/DY-FLOAT(IY)
      XR(I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY1)+(1.-Q)*P*XF(IX1,IY) + P*Q*XF(IX1,IY1)
10  CONTINUE
    RETURN
  END SUBROUTINE BOUND