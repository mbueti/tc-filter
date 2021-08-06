  SUBROUTINE rodist
      !   parameter(nmx=64)
    common /vect/xvect(nmx),yvect(nmx)
    COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY

      !   pi=4.0*atan(1.0)
      !   PI180 = 4.*ATAN(1.0)/180.
    yo=yold*pi180
    fact=cos(yo)
    xc=xold-xcorn
    yc=yold-ycorn

    do 10 ip=1,nmx
      theta=float(ip-1)/float(nmx)*2.*pi
      r=rovect(ip)
      xvect(ip)=r*cos(theta)/fact +xc
      yvect(ip)=r*sin(theta) +yc
10  continue
    return
  END SUBROUTINE rodist