subroutine rodist
  implicit none

  real :: xv, yv, rb, xold, yold, xcorn, ycorn, factr, yo, fact, xc, yc, r, theta
  integer :: nnn, ienv, ix, iy, ip
  real, dimension(nmx) :: xvect, yvect, rovect
  common /vect/ xvect, yvect
  COMMON  /IFACT/ NNN,rovect,RB,IENV
  COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
  
  yo=yold*pi180
  fact=cos(yo)
  xc=xold-xcorn
  yc=yold-ycorn
  
  do 10 ip=1,nmx
  theta=float(ip-1)/float(nmx)*2.*pi
  r=rovect(ip)
  xvect(ip)=r*cos(theta)/fact +xc
  yvect(ip)=r*sin(theta) +yc
10 continue
  return
end subroutine rodist
