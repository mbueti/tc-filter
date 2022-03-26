SUBROUTINE calcr(RO,RTAN,xc,yc,yold,u,v)
  implicit NONE

  integer, PARAMETER :: nmx=64
  integer :: ix, iy, ix1, iy1, i
  real :: ddel, dtha, pi, pi180, fact, dx, dy, theta, x, y, x1, y1, p, q
  real, DIMENSION(nmx) :: XR(NMX)
  real, dimension(imx,jmx), intent(in) :: u,v
  real, intent(in) :: ro, xc, yc, yold
  real, intent(out) :: rtan
  
  COMMON  /TOTAL/ DDEL,dtha
  
  fact=cos(yold)
  DX=DDEL/PI180
  DY=DTHA/PI180
  DO I=1,NMX
    THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    X=RO*COS(THETA)/fact +XC +1.
    if (NINT(X).gt.360) X = X - 360
    Y=RO*SIN(THETA)+YC +1.
    X1=X+DX
    if (NINT(X1).gt.360) X1 = X1 - 360
    Y1=Y+DY
    IX=NINT(X/DX)
    IY=NINT(Y/DY)
    IX1=NINT(X1/DX)
    IY1=NINT(Y1/DY)
    IX = MODULO(IX-1,imx) +1
    IX1 = MODULO(IX1-1,imx) +1
    IY = MAX(IY,1)
    IY1 = MAX(IY1,1)
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
    xr(i)=-sin(theta)*((1.-P)*(1.-Q)*u(IX,IY)+(1.-P)*Q*u(IX,IY1)+(1.-Q)*P*u(IX1,IY)+P*Q*u(IX1,IY1)) &
          +cos(theta)*((1.-P)*(1.-Q)*v(IX,IY)+(1.-P)*Q*v(IX,IY1)+(1.-Q)*P*v(IX1,IY)+P*Q*v(IX1,IY1))
  end do
  RTAN = 0.0
  
! calculate azimuthally averaged tangential wind at radius ro

  DO I = 1 , NMX
    RTAN = RTAN + XR(I)
  end do
  RTAN = RTAN/FLOAT(NMX)
  RETURN
END SUBROUTINE calcr
