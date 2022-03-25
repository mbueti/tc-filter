  SUBROUTINE calcr(RO,RTAN,xc,yc,yold,u,v)
    ! PARAMETER ( nmx=64)
    ! PARAMETER (IMX=640, JMX=320)
  ! REAL, DIMENSION(NMX) :: XR
    DIMENSION XR(NMX),u(imx,jmx),v(imx,jmx)
!   
    COMMON  /TOTAL/ DDEL,dtha
!        COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
         !
    ! PI = 4.*ATAN(1.0)
    ! PI180 = PI/180.
!        FACT =  COS(YOLD*PI180)
!        yo=yold*pi180
    fact=cos(yold)
!   
    DX=DDEL/PI180
    DY=DTHA/PI180
        !  XC = (XOLD-XCORN)*DX
        !  YC = (YOLD-YCORN)*DY
    !DO I=1,NMX
    !  THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    !  X=RO*COS(THETA)+XC
    !  Y=RO*SIN(THETA)+YC
    !  if (X.lt.0) then
    !    print *, 'RO=',RO
    !    print *, 'theta=',THETA
    !    print *, 'sin=',sin(THETA)
    !    print *, 'XC=',YC
    !    print *, 'X=', X
    !  end if
    !  if (Y.lt.0) then
    !    print *, 'RO=',RO
    !    print *, 'theta=',THETA
    !    print *, 'cos=',cos(THETA)
    !    print *, 'YC=',YC
    !    print *, 'Y=', Y
    !  end if
    ! end do
    ! stop
    ! print *, 'yc=',yc
!    DO 10 I=1,NMX
!      THETA= 2.*PI*FLOAT(I)/FLOAT(NMX)
!      !   X=RO*fact*COS(THETA)+XC
!      X=RO*COS(THETA)+XC
!      Y=RO*SIN(THETA)+YC
!      IX=CEILING(X/DX)+1
!      IY=CEILING(Y/DY)+1
!      print *, 'THETA=',THETA
!      print *, 'IX=',IX
!      print *, 'IY=',IY
!10  continue
    DO 10 I=1,NMX
      THETA= 2.*PI*FLOAT(I)/FLOAT(NMX)
      !   X=RO*fact*COS(THETA)+XC
      X=RO*COS(THETA)+XC
      Y=RO*SIN(THETA)+YC
      IX=CEILING(X/DX)+1
      IY=CEILING(Y/DY)+1
      IX1=IX+1
      IY1=IY+1
      P=X/DX-FLOAT(IX)
      Q=Y/DY-FLOAT(IY)
!        XR(I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY+1)
!      1      +  (1.-Q)*P*XF(IX+1,IY) + P*Q*XF(IX+1,IY+1)
      ! print *, 'RO=',RO
      ! print *, 'XC=',XC
      ! print *, 'X=',X
      ! print *, 'DX=',DX
      ! print *, "IX=",IX
      ! print *, 'YC=',YC
      ! print *, "Y=",Y
      ! print *, "DY=",DY
      ! print *, "IY=",IY
      !if (IX.lt.1) then
      !  !print *, 'RO=',RO
      !  print *, 'theta=',THETA
      !  print *, 'cos=',cos(THETA)
      !  print *, 'XC=',YC
      !  print *, 'DX=', DX
      !  print *, 'X=', X
      !  print *, 'IX=', IX
      !end if
      !if (IY.lt.1) then
      !  print *, 'RO=',RO
      !  print *, 'theta=',THETA
      !  print *, 'sin=',sin(THETA)
      !  print *, 'YC=',YC
      !  print *, 'DY=', DY
      !  print *, 'Y=', Y
      !  print *, 'IY=', IY
      !end if
      if(IX.le.0) IX=IX+IMX
      if(IX.gt.IMX) IX=MOD(IX,IMX)+1
      if(IY.le.0) IY=IY+JMX
      if(IY.gt.JMX) IY=MOD(IY,JMX)+1

      IX1=IX+1
      IY1=IY+1
      
      if(IX1.gt.IMX) IX1=MOD(IX1,IMX)+1
      if(IY1.gt.JMX) IY1=MOD(IY1,JMX)+1
      XR(I) = -sin(THETA)*((1.-P)*(1.-Q)*u(IX,IY)+(1.-P)*Q*u(IX,IY1)+(1.-Q)*P*u(IX1,IY)+P*Q*u(IX1,IY1))+ &
               cos(THETA)*((1.-P)*(1.-Q)*v(IX,IY)+(1.-P)*Q*v(IX,IY1)+(1.-Q)*P*v(IX1,IY)+P*Q*v(IX1,IY1))
10  CONTINUE
    RTAN = 0.0
!
! calculate azimuthally averaged tangential wind at radius ro
!
    DO 40 I=1,NMX
      RTAN = RTAN + XR(I)
40  CONTINUE
    RTAN = RTAN/FLOAT(NMX)
    RETURN
  END SUBROUTINE calcr