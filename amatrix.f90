  SUBROUTINE amatrix
    common /matrix/ a(nmx,nmx),capd2
    common /vect/xvect(nmx),yvect(nmx)
    COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
    
    yo=yold*pi180
    fact=cos(yo)
    capd2=(2.25)*(2.25)
    do ip=1,nmx
      do jp=ip,nmx
        dpij=(fact*(xvect(ip)-xvect(jp)))**2 +(yvect(ip)-yvect(jp))**2
        a(ip,jp)= exp(-dpij/capd2)
        a(jp,ip)= a(ip,jp)
      end do
    end do
    return
  END SUBROUTINE amatrix