SUBROUTINE CALCRa(RO,RTAN,iang,dist)
  implicit none

  real, intent(in) :: ro, dist
  integer, intent(in) :: iang
  real, intent(out) :: rtan

  integer :: id1, id2, ix, ix1, iy, iy1
  real :: xv, yv, ddel, dtha, pi, pi180, fact, dx, dy, xc, yc, theta, x, dmmm, tang, del, tha, xf, ds, xold, yold, &
         xcorn, ycorn, factr, p, q, x1, y, y1
       COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX), DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
!
       COMMON  /TOTAL/ DDEL,DTHA
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
!
          PI = 4.*ATAN(1.0)
       PI180 = 4.*ATAN(1.0)/180.
       FACT =  COS(YOLD)
!
       DX=DDEL/PI180
       DY=DTHA/PI180
       XC = (XOLD-XCORN)
       YC = (YOLD-YCORN)

!
        THETA= 2.*PI*FLOAT(iang-1)/FLOAT(NMX)
        X=RO*COS(THETA)+XC
!        print *, cos(theta)
!        print *, 'x=', x, ' ro=', ro, ' xc=', xc, ' theta=', theta
        x = max(mod(x, 361.0),1.0)
        if (NINT(X).gt.360) X = X - 360
        if (NINT(X).le.0) X = X + 360
        Y=RO*SIN(THETA)+YC
        X1=X+DX
        x1 = max(mod(x1, 361.0),1.0)
        if (NINT(X1).gt.360) X1 = X1 - 360
        if (NINT(X1).le.0) X1 = X1 + 360
        Y1=Y+DY
        IX=nint(X/DX)
        IY=nint(Y/DY)
        IX1=NINT(X1/dx)
        IY1=nint(Y1/dy)
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
!        print *, 'ix=', ix, ' iy=', iy
       rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY1) + (1.-Q)*P*XF(IX1,IY) + P*Q*XF(IX1,IY1)
10     CONTINUE
!
!
         RETURN
         END
