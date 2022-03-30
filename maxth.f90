SUBROUTINE maxth(dumu,dumv,dxc,dyc,rmxlim,tw)
  IMPLICIT NONE

  integer, parameter :: lgth=60
  real, intent(inout) :: dxc, dyc, rmxlim
  real, dimension(imx,jmx), intent(in) :: dumu, dumv
  real, dimension(imx,jmx), intent(inout) :: tw
  real :: xv, yv, alim, ddel, dtha, deltar, dx, dy, fact, factr, hmax, rbd, rfavg, rfind, rmxavg, tanp, rxx, theta, &
          xcorn, ycorn, xold, yold, rmxpos, rtan, xcen,xcn, xctest, ycen, ycn, yctest, yyo
  integer :: i, iend, if, iflag, im, ipos, ir, ist, ix, iy, ixc,iyc, ixc1, iyc1, j, jend, jst, npts
  real, dimension(imx,jmx) :: del, tha, xf, ds, tang, th, tanmx
  real, dimension(imx,jmx,2) :: dmmm
  integer, dimension(7,7) :: itpos
  real, dimension(7,7) :: tmax
  real, dimension(7,7,lgth) :: tprof
  real, dimension(iimx) :: tanavg
  logical :: iexit
  COMMON  /TOTAL/ DDEL,dtha
  COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
  COMMON /WINDS/ DMMM,TANG,DEL,THA,XF,DS
  COMMON /scale/rmxavg,rfind
  
  fact=cos(yold)
  deltar=0.1
  dxc=xold/pi180-xcorn
  dyc=yold/pi180-ycorn
  ixc=nint(dxc)+1
  iyc=nint(dyc)+1
  ist=ixc-3
  jst=iyc-3
  iend=ixc+3
  jend=iyc+3
  npts=7
  print *, 'xold, yold', xold, yold
  print *, 'xcorn, ycorn', xcorn, ycorn
  print *, 'ixc, iyc', ixc, iyc
  print*,'ist,iend',ist,jst,iend,jend
!
!  compute radial profile of azimuthal avg. tang. wind at each pt
!
  do i=ist,iend
    do j=jst,jend
      IM = MODULO(I-1,IMX)+1
      xcen=(del(im,j)-del(1,1))/pi180 +1.
      ycen=(tha(im,j)-tha(1,1))/pi180 +1.
      yyo=tha(im,j)
      do ir=1,lgth
        rbd=float(ir)*0.2
        call bound2(dumu,dumv,tanp,rbd,xcen,ycen,yyo)
        tprof(i-ist+1,j-jst+1,ir)=tanp
      END DO
    END DO
  END DO
!
!  find the first relative maximum along each azimuthal direction
!  find the position of the largest relative maximum
!
  hmax=0.
  do i=1,npts
    iexit = .FALSE.
    do j=1,npts
      do ir=2,lgth-1
        if(tprof(i,j,ir).gt.tprof(i,j,ir-1).and.tprof(i,j,ir).gt.tprof(i,j,ir+1)) then
          tmax(i,j)=tprof(i,j,ir)
          itpos(i,j)=jmx*(ist+i-1)+j+jst-1
          if (tmax(i,j).gt.hmax) THEN
            hmax = tmax(i,j)
            ipos = itpos(i,j)
            rmxpos = float(if)*0.2
          end if
          IEXIT = .TRUE.
          exit
        end if
      END DO
      if (.not.iexit) then
        tmax(i,j)=tprof(i,j,1)
        itpos(i,j)=101
        hmax=amax1(hmax,tmax(i,j))
        if(hmax.eq.tmax(i,j)) ipos=itpos(i,j)
      end if
    end do
  end do
  print *, 'hmax, rmxpos, 5x rmxpos are'
  print*,hmax,rmxpos,rmxpos/0.2
  print *, 'tmax is'
  print *,((tmax(i,j),i=1,npts),j=1,npts)
  print *, 'itpos is'
  print *,((itpos(i,j),i=1,npts),j=1,npts)
  print *, 'ipos is ', ipos, mod(ipos,100)
!
!
!
!  use position of the largest relative maximum as the adjusted
!  center location
!
  ycn=float(mod(ipos,jmx))-1.
  xcn=float(ipos/jmx)-1.
  ixc=int(xcn)+1
  iyc=int(ycn)+1
  ixc = modulo(ixc-1, imx)+1
  ixc1 = modulo(ixc, imx)+1
  print *, 'ipos=', ipos
  print *, 'xcn=', xcn
  print *, 'ycn=', ycn
  xctest=(xcn+xcorn)*pi180
  yctest=(ycn+ycorn)*pi180
  print *, 'xctest=', xctest/pi180
  print *, 'yctest=', yctest/pi180
!
!  recompute the tangential wind component based on new center
!
  fact = cos(tha(1,iyc))
  print *,'in maxth',ycn,xcn,fact,xctest/pi180,yctest/pi180
  do j=1,jmx
    do i=1,imx
      dx=(del(i,j)-xctest)*fact
      dy=(tha(i,j)-yctest)
      if(dx.ne.0.)theta =atan2(dy,dx)
      if(dy.gt.0..and.dx.eq.0.)theta =90.*pi180
      if(dy.lt.0..and.dx.eq.0.)theta =270.*pi180
      tw(i,j)=-dumu(i,j)*sin(theta) +dumv(i,j)*cos(theta)
      if(i.eq.ixc.and.j.eq.iyc)print*,i,j,dumu(i,j),dumv(i,j),theta,tw(i,j),dx,dy,'check everything'
    END do
  END DO 
  write(77,*)ixc,iyc
 
  iflag=0
  hmax=0.
  rmxavg=0.
  do ir=3,iimx
    rxx=float(ir)*deltar
    call calcr( rxx,rtan,xcn,ycn,yctest,dumu,dumv )
    tanavg(ir)=rtan
    if(tanavg(ir-2).lt.tanavg(ir-1).and.tanavg(ir).lt.tanavg(ir-1).and.iflag.eq.0)then
      hmax=tanavg(ir-1)
      rmxavg=rxx-deltar
      iflag=1
    end if
  end do
  print*,'found rmxavg ',rmxavg,hmax
  dxc=xcn
  dyc=ycn
      
700 format(10f6.1)
  print 700,tanavg
  call findra( dxc,dyc,yctest,rmxavg,rfavg,tanavg)
!
  alim = .75
  print*
  print*,'a factor to determine rmxlim: ',alim
!
  rmxlim = alim*rmxavg + (1.-alim)*rfavg
  print*,'found rfavg ',rfavg,rmxlim,dxc,dyc
  return
END SUBROUTINE MAXTH
