  SUBROUTINE SEPAR(XD)
!
!  SEPERATES A FIELD INTO HURRICANE COMPONENT AND REMAINDER
!
    PARAMETER(nmx1=nmx+1,nmx2=nmx*2,nmx6=nmx*6)
    DIMENSION XR(NMX),XD(IMX,JMX)
!
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX), &
                   DEL(IMX,JMX),THA(IMX,JMX),TANP(IMX,JMX),DS(IMX,JMX)
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
    COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
    COMMON /XXX/  XF(IMX,JMX),XC,YC,DX,DY
    COMMON /TOTAL/ DDEL,DTHA
!
! new arrays
    dimension b(nmx),w(nmx),ab(nmx,nmx1),ipvt(nmx),wrk(nmx6),iwrk(nmx2)
    common /matrix/ a(nmx,nmx),capd2
    common /vect/xvect(nmx),yvect(nmx)
!
    DATA XR/nmx*0./
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
    do 22 i=1,nmx
      ro=amax1(ro,rovect(i))
22  continue
    print*,'ro=',ro,capd2,a(1,1),a(2,1)
    FACT = COS(YOLD*PI180)
!C
!C   XC IS THE I POSITION OF THE CENTER OF THE OLD VORTEX
!C   YC IS THE J POSITION OF THE CENTER OF THE OLD VORTEX
!C   DDEL IS THE LONG. IN RADIANS OF THE OUTER NEST
!C   DTHA IS THE LAT.  IN RADIANS OF THE OUTER NEST
!C
! no fact here
!      DX=FACT*DDEL/PI180
!
    dx=ddel/pi180
    DY=DTHA/PI180
!C
!c
    XC = (XOLD-XCORN)*DX
    YC = (YOLD-YCORN)*DY
    IS=INT((XC-RO/fact)/DX) +1.
    IE=INT((XC+RO/fact)/DX + 1.)
    JS=INT((YC-RO)/DY) +1.
    JE=INT((YC+RO)/DY + 1.)
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
    CALL BOUND(NMX,XR,rovect)
    print*,'here is rovect ',rovect
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
    DO 10 IX=IS,IE
      DO 11 JY=JS,JE
        ro=romax
!       X=XC-RO +DX*(IX-IS)
!       Y=YC-RO +DY*(JY-JS)
!       X= DX*float(IX)
!       Y= DY*float(JY)
        x=del(ix,jy)/pi180 -xcorn
        y=tha(ix,jy)/pi180 -ycorn
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
        print *, 'XRO=', XRO
!
! Now add new code to compute distance from each gridpt. to rovect pts
!
        do 12 ip=1,nmx
          dpij= (fact*(x-xvect(ip)))**2 +(y-yvect(ip))**2
          print *, 'ip=', ip
          print *, 'dpij=', dpij
          print *, 'capd2=', capd2
          print *, -dpij/capd2
          if (dpij/capd2 > 100) then
            b(ip) = 0
          else
            b(ip) = exp(-dpij/capd2)
          endif
          print *, 'b(ip)=', b(ip)
12      continue
        print *, 'b=', b
!
!
        do 44 ip=1,nmx
          do 43 jp=1,nmx
            ab(ip,jp)=a(ip,jp)
43        continue
          ab(ip,nmx1)=b(ip)
44      continue
        print *, 'ab=', ab
!
!  solve system using constrained least squares method
!
        print *, 'calling wnnls'
        call wnnls(ab,nmx,0,nmx,nmx,0,[1.],w,rnm,md,iwrk,wrk,nmx6,nmx2)
!
        temp=0.
        do 20 ip=1,nmx
          temp=temp +w(ip)*xr(ip)
20      continue
!       xh(ix,jy)=xf(ix,jy)-temp
        xd(ix,jy)=temp
11    CONTINUE
10  CONTINUE
    RETURN
  END SUBROUTINE SEPAR