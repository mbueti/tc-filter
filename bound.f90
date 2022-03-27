SUBROUTINE BOUND(NMX,XR,ro)
  IMPLICIT NONE

  integer, intent(in) :: nmx
  real, dimension(nmx), intent(in) :: ro
  real, dimension(nmx), intent(out) :: xr
  real :: xv, yv, xc, yc, dx, dy, xold, yold, xcorn,ycorn, factr, fact, theta, x, y, p, q
  integer :: ix, iy, ix1, iy1, i
  real, dimension(imx, jmx) :: xf
  
  COMMON  /XXX/  XF,XC,YC,DX,DY
  COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY

  fact=cos(yold*pi/180.)
  DO I=1,NMX
    THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    X=RO(i)/fact*COS(THETA)+XC +1.
    Y=RO(i)*SIN(THETA)+YC +1.
    IX=NINT(X/DX)
    IY=NINT(Y/DY)
    IX1=IX+1
    IY1=IY+1
    IX = MODULO(IX-1,imx) +1
    IX1 = MODULO(IX1-1,imx) +1
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
    print *, 'pi=',pi
    print *, 'fact=',fact
    print *, 'theta=',theta
    print *, 'ro=',RO(i)
    print *, 'Y=',Y
    print *, 'YC=',YC
    print *, 'YOLD=',YOLD
    print *, 'IY=', IY
    print *, 'IY1=', IY1
    print *, 'sin=', sin(theta)
    XR(I)=(1.-P)*(1.-Q)*XF(IX,IY)+(1.-P)*Q*XF(IX,IY1)+(1.-Q)*P*XF(IX1,IY)+P*Q*XF(IX1,IY1)
  END DO
  RETURN
END SUBROUTINE BOUND
