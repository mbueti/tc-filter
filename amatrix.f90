subroutine amatrix
  implicit none
  real :: xv, yv, capd2, fact, dpij, factr, rb, xcorn, ycorn, xold, yold, yo
  integer :: ienv, ip, jp, ix, iy, nnn
  real, dimension(nmx,nmx) :: a
  real, dimension(nmx) :: xvect, yvect, rovect
  common /matrix/ a,capd2
  common /vect/xvect,yvect
  COMMON /IFACT/NNN,rovect,RB,IENV
  COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
       
  yo=yold*pi180
  fact=cos(yo)
  capd2=(2.25)*(2.25)
  do 10 ip=1,nmx
  do 10 jp=ip,nmx
    dpij=(fact*(xvect(ip)-xvect(jp)))**2 +(yvect(ip)-yvect(jp))**2
    a(ip,jp)= exp(-dpij/capd2)
    a(jp,ip)= a(ip,jp)
10 continue
100 format(5f8.4)
  return
end subroutine amatrix
