SUBROUTINE SEPAR(XD)
!
!  SEPERATES A FIELD INTO HURRICANE COMPONENT AND REMAINDER
!
  implicit NONE
  integer, PARAMETER :: nmx1=nmx+1,nmx2=nmx*2,nmx6=nmx*6
  real, dimension(imx, jmx), intent(inout) :: xd
  integer :: ix, iy, i, ie, ienv, im, ip, is, j, je, js, jy, jp,n1, n2, nnn, md
  real :: xv, yv, xold, yold, xcorn, ycorn, capd2, ddel, dtha,delth, delx, dely, dr, dpij, ro, romax, temp, x, y,dx, dy, fact, &
          factr, theta, xc, yc, xro,rb, rnm
  real, dimension(nmx) :: xr
  real, dimension(nmx,nmx) :: a
  real, dimension(nmx,nmx1) :: ab
  real, dimension(nmx) :: b, w, xvect, yvect, rovect
  integer, dimension(nmx) :: ipvt
  real, dimension(121) :: wrk
  integer, dimension(25) :: iwrk
  real, dimension(imx, jmx) :: del, tha, tanp, ds, tang, xf
  real, dimension(imx, jmx, 2) :: dmmm
  
  COMMON /WINDS/ DMMM,TANG,DEL,THA,TANP,DS
  COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
  COMMON  /IFACT/NNN,rovect,RB,IENV
  COMMON /XXX/  XF,XC,YC,DX,DY
  COMMON /TOTAL/ DDEL,DTHA
  
! new arrays
  common /matrix/ a,capd2
  common /vect/xvect,yvect
  
  DATA XR/64*0./
!
!  XC,YC ARE HURRICANE COORDINATES
!  RO  IS RADIUS AT WHICH HURRICANE COMPONENT OF FIELD GOES TO ZERO
!  XR ARRAY CONTAINS THE FIELD VALUES OF 12 EQUALLY SPACED POINTS
!     ON CIRCLE OF RADIUS RO CENTERED AT XC,YC
!
!  set ro to be max value of rovect
!
!
  ro=0.
  DO i=1,nmx
    ro=amax1(ro,rovect(i))
  END DO
  ro=amin1(ro, 30.0)
!        print*,'rovect=',rovect
!        print*,'ro=',ro,capd2,a(1,1),a(2,1)
!        ro = ro*1.5
  FACT =  COS(YOLD*PI180)
!
!   XC IS THE I POSITION OF THE CENTER OF THE OLD VORTEX
!   YC IS THE J POSITION OF THE CENTER OF THE OLD VORTEX
!   DDEL IS THE LONG. IN RADIANS OF THE OUTER NEST
!   DTHA IS THE LAT.  IN RADIANS OF THE OUTER NEST
!
! no fact here
!      DX=FACT*DDEL/PI180
!
  dx=ddel/pi180
  DY=DTHA/PI180
!
!
  XC = (XOLD-XCORN)*DX
  YC = (YOLD-YCORN)*DY
  IS=NINT((XC-RO/fact)/DX) +1.
  IE=NINT((XC+RO/fact)/DX + 1.)
  JS=NINT((YC-RO)/DY) +1.
  JE=NINT((YC+RO)/DY + 1.)
!
  DO J = 1 , JMX
    DO I = 1 , IMX
      XF(I,J)  = XD(I,J)
    END DO
  END DO
!
!  SUBROUTINE BOUND COMPUTES FIELD VALUES OF ARRAY XR USING
!         BILINEAR INTERPOLATION
!
!
  Print*, 'calling BOUND from SEPAR '
!        print*,'here is xr ',xr
  CALL BOUND(NMX,XR,rovect)
!        print*,'here is rovect ',rovect
!
!  xrop(nmx) are the interpolated values of the disturbance
!   field at the rovect pts
!
! romax is the maximum value in rovect(nmx). Within the loop a local
! ro is computed for use in the separation. At the start of the loop
! ro is again set to romax to define the domain.
!
!
!
  w=0.
  romax=ro
!
  IF (ie.gt.imx.or.je.gt.jmx.or.is.lt.1.or.js.lt.1) THEN
    print *, 'ro=',ro
    print *, 'fact=',fact
    print *, 'xold=', xold, 'yold=', yold
    print *, 'xc=', xc, 'yc=', yc
    print *, 'js=', JS, ' je=', JE
  END IF
  DO 10 IX=IS,IE
    DO 11 JY=JS,JE
      ro=romax
      im = modulo(ix-1, imx)+1
      x=del(im,jy)/pi180 -xcorn
      y=tha(im,jy)/pi180 -ycorn
      delx=(x-xc)*fact
      dely=(y-yc)
      DR=SQRT((delx)**2 +(dely)**2)
      IF(DR.GT.RO)GOTO11
      IF(delx.ne.0.) THETA=ATAN((dely)/(delx))
      if(delx.eq.0..and.dely.lt.0.)theta=270.*pi180
      if(delx.eq.0..and.dely.gt.0.)theta=90. *pi180
      IF(delx.LT.0.)THETA=THETA+PI
      IF(THETA.LT.0.)THETA=2.*PI+THETA
      N1=INT(THETA*NMX/(2.*PI))
      IF(N1.GT.nmx)PRINT*,N1,THETA*57.296
      IF(N1.LT.0)PRINT*,N1,THETA*57.296
      N2=N1+2
      IF(N2.GT.NMX)N2=N2-NMX
      DELTH=THETA- 2.*PI*FLOAT(N1)/FLOAT(NMX)
!
      ro=delth*float(nmx)/(2.*pi)*(rovect(n2)-rovect(n1+1))+rovect(n1+1)
      IF(DR.GT.ro)GOTO11
      XRO=DELTH*FLOAT(NMX)/(2.*PI)*(XR(N2)-XR(N1+1)) +XR(N1+1)
!
! Now add new code to compute distance from each gridpt. to rovect pts
!
      DO ip=1,nmx
        dpij= (fact*(x-xvect(ip)))**2 +(y-yvect(ip))**2
        b(ip)=exp(-dpij/capd2)
      END DO

      DO ip=1,nmx
        DO jp=1,nmx
          ab(ip,jp)=a(ip,jp)
        END DO
        ab(ip,nmx1)=b(ip)
      END DO
!
!  solve system using constrained least squares method
!
      call wnnls(ab,nmx,0,nmx,nmx,0,[1.],w,rnm,md,iwrk,wrk)
            
      temp=0.
      DO ip=1,nmx
        temp=temp +w(ip)*xr(ip)
      END DO
      xd(ix,jy)=temp
11  CONTINUE
10 CONTINUE
  RETURN
END SUBROUTINE SEPAR
