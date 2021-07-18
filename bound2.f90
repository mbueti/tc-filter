  subroutine bound2(u,v,tanuv,r0,xc,yc,yyo)
    DIMENSION u(imx,jmx),v(imx,jmx),tani(nmx)
!   COMMON  /POSIT/ XOLD,YOLD,XCORN,YCORN,Rxx,XV,YV
    COMMON  /TOTAL/ DDEL,dtha
!   COMMON  /XXX/    XC,YC,DX,DY
    ! PI = 4.*ATAN(1.0)
    ! pi180=pi/180.
    dx=ddel/pi180
    dy=dtha/pi180
    fact=cos(yyo)
    arad=6371.
     r0=r0*111.19
    DO 10 I=1,NMX
    THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    X=(R0*COS(THETA))/(arad*fact*pi180)+XC
    Y=(R0*SIN(THETA))/(arad*pi180)+YC
    IX=INT(X/DX)
    IY=INT(Y/DY)
    IX1=IX+1
    IY1=IY+1
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
!      ui=(1.-P)*(1.-Q)*u(IX,IY) +(1.-P)*Q*u(IX,IY+1)
!    1      +  (1.-Q)*P*u(IX+1,IY) + P*Q*u(IX+1,IY+1)
!      vi=(1.-P)*(1.-Q)*v(IX,IY) +(1.-P)*Q*v(IX,IY+1)
!    1      +  (1.-Q)*P*v(IX+1,IY) + P*Q*v(IX+1,IY+1)
!      tani(i)=-sin(theta)*ui +cos(theta)*vi
    tani(i)=-sin(theta)*((1.-P)*(1.-Q)*u(IX,IY) +(1.-P)*Q*u(IX,IY+1)+  (1.-Q)*P*u(IX+1,IY) + P*Q*u(IX+1,IY+1)) &
            +cos(theta)*((1.-P)*(1.-Q)*v(IX,IY) +(1.-P)*Q*v(IX,IY+1)+  (1.-Q)*P*v(IX+1,IY) + P*Q*v(IX+1,IY+1))
10  CONTINUE
!
    tanuv=0.
    do 20 i=1,nmx
      tanuv=tanuv +tani(i)/float(nmx)
20  continue
    RETURN
  END SUBROUTINE bound2