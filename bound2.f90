subroutine bound2(u,v,tanuv,r0,xc,yc,yyo)
  IMPLICIT NONE

  integer :: ix, iy, ix1, iy1, i, j
  real :: dx, dy, theta, p, q
  real ::  fact, x, y, arad, ddel, dtha, r
  real, intent(out) :: tanuv
  real, intent(in) :: r0, xc, yc, yyo
  REAL, DIMENSION(imx, jmx), intent(in) :: u,v
  REAL, DIMENSION(nmx) :: tani
  COMMON  /TOTAL/ DDEL,dtha

  dx=ddel/pi180
  dy=dtha/pi180
  fact=cos(yyo)
  arad=6371.
  r=r0*111.19
  DO I=1,NMX
    THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    X=(R*COS(THETA))/(arad*fact*pi180)+XC
    Y=(R*SIN(THETA))/(arad*pi180)+YC
    IX=NINT(X/DX)
    IY=NINT(Y/DY)
    IX1=IX+1
    IY1=IY+1
    IX = MODULO(IX-1,imx) +1
    IX1 = MODULO(IX1-1,imx) +1
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
    tani(i)=-sin(theta)*((1.-P)*(1.-Q)*u(IX,IY)+(1.-P)*Q*u(IX,IY1)+(1.-Q)*P*u(IX1,IY)+P*Q*u(IX1,IY1)) &
            +cos(theta)*((1.-P)*(1.-Q)*v(IX,IY)+(1.-P)*Q*v(IX,IY1)+(1.-Q)*P*v(IX1,IY)+P*Q*v(IX1,IY1))
  end do

  tanuv=0.
  do i=1,nmx
    tanuv=tanuv +tani(i)/float(nmx)
  end do
  RETURN
END subroutine bound2
