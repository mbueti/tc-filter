SUBROUTINE calct(deltar,xc,yc,yo,xf,rmxlim)
  implicit none
  
!  calculates the radial profile for eight azimuthal angles

  real, intent(in) :: deltar, xc, yc, yo
  real, intent(out) :: rmxlim
  real :: xv, yv, arad, bfact, ddel, dtha,ro, theta, x, x1, y, y1, p, q, dx, dy, xcorn, ycorn, xold, yold,rmxavg, rfind, rdistl, &
          fact, factr
  integer :: ix, ix1, iy, iy1, i, iaa, iravg, irend, ir, irvgu, id1, id2, mdw, m, n, l
  real, DIMENSION(imx,jmx),intent(in) :: xf
  integer, dimension(nmx) :: idst
  real, dimension(nmx) :: hmax,rmax,dist
  real, dimension(iimx,nmx) :: xr
  
  common /scale/rmxavg,rfind
  common /pass/xr,dist
  COMMON  /TOTAL/ DDEL,dtha
  COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
  
  arad=6371.
  print*,'calct',deltar,xc,yc,yo,yold
  fact=cos(yo)
!
!      set the factor to determine lower limit of dist search
!
  bfact = .5
  rdistl =  bfact*rmxavg
!
  print*
  print*,'b factor to determine rdistl: ',bfact
!
!
! assume the maximum wind is within rfavg of center (but <10.8 also)
!
  irend =int(rmxlim/deltar)
!
!
  iravg =int(rdistl/deltar)
  print *, 'irend=', irend, ' iravg=', iravg
  if (irend.eq.0.or.iravg.eq.0) THEN
    return
  end if
!
  print*,'lower limit and radius of dist search: ',rdistl,iravg
  print*,'upper limit and radius of dist search: ',rmxlim,irend
!
  DX=DDEL*(1./PI180)
  DY=DTHA/PI180
!
!  angle loop
!
  DO I=1,NMX
    THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
    do ir=1,iimx
      ro=float(ir)*deltar
      X=(RO*COS(THETA)) +XC
      Y=(RO*SIN(THETA)) +YC + 91.0
      X1=X+DX
      Y1=Y+DY
      IX=NINT(X/DX)
      IY=NINT(Y/DY)
      IX1=NINT(X1/DX)
      IY1=NINT(Y1/DY)
      IX = MODULO(IX-1,IMX) + 1
      IX1 = MODULO(IX1-1,IMX) + 1
      P=X/DX-FLOAT(IX)
      Q=Y/DY-FLOAT(IY)
      if (IY.gt.jmx.or.iy1.gt.jmx) THEN
        CYCLE
      end if
       XR(ir,I)=(1.-P)*(1.-Q)*XF(IX,IY)+(1.-P)*Q*XF(IX,IY1)+(1.-Q)*P*XF(IX1,IY)+P*Q*XF(IX1,IY1)
    end do
!
! find relative max after which ro check begins
!
    idst = 0
!
    hmax=-10.e10
    do ir=iravg,irend
      hmax(i)=amax1(hmax(i),xr(ir,i))
      if(hmax(i).eq.xr(ir,i))dist(i)=float(ir)*deltar
      if(hmax(i).eq.xr(ir,i))idst(i)=ir
    end do
!
!  if the max. value is also the endpt it maynot be a relative max
!  check backwards from irend for the last relative max
!
!
    irvgu = iravg+1
    if(irvgu.le.2) then
      irvgu = 3
    end if
!
    if(irend.eq.idst(i).or.iravg.eq.idst(i))then
      do ir=irend,irvgu,-1
        if(xr(ir-1,i).lt.0.) exit
        if(xr(ir-1,i).gt.xr(ir,i).and.xr(ir-1,i).ge.xr(ir-2,i)) then
          dist(i)=float(ir-1)*deltar
          print*,'readjusting dist',dist(i),hmax(i),xr(ir-1,i),rmxlim
          exit
        end if
      end do
    end if
    
    if(idst(i).lt.iravg)then
      print*,'lower limit check, dist changed to rmxlim, i is: ',i
      dist(i) = rmxlim
    end if
  end do

  do iaa = 1 , nmx
    WRITE(94,*)iaa,dist(iaa)
  end do
  
  write(94,*)rdistl,rmxlim
  
  dist=dist*1.1
  print*,'relative max found '
  print 4400,dist
4400 format(25f4.1)
  write(4,400) (float(ir)*deltar,(xr(ir,i)/100.,i=1,1),ir=1,iimx)
400 format(25f5.1)
  RETURN
END SUBROUTINE calct
