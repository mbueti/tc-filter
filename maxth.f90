  SUBROUTINE maxth(dumu,dumv,dxc,dyc,rmxlim,tw)
    parameter(lgth=30,iimx=110)
    ! parameter(nmx=64,imx=640,jmx=320,lgth=30,iimx=110)
    dimension dumu(imx,jmx),dumv(imx,jmx),tw(imx,jmx)
    dimension tprof(7,7,lgth),itpos(7,7),tmax(7,7),tanavg(iimx)
    COMMON  /TOTAL/ DDEL,dtha
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
!
    common /scale/rmxavg,rfind
!
    ! pi=4.*atan(1.0)
    ! PI180 = pi/180.
    fact=cos(yold)
    deltar=0.1
    dxc=xold/pi180-xcorn
    dyc=yold/pi180-ycorn
    ixc=int(dxc)+1
    iyc=int(dyc)+1
    print *,'ixc=',ixc
    print *,'iyc=',iyc
    print *,'dxc=',dxc
    print *,'dyc=',dyc
    print *,'xold=',xold
    print *,'xcorn=',xcorn
    print *,'yold=',yold
    print *,'ycorn=',ycorn
    ist=ixc-3
    jst=iyc-3
    iend=ixc+3
    jend=iyc+3
    npts=7
    print*,'ist,iend',ist,jst,iend,jend
!
!  compute radial profile of azimuthal avg. tang. wind at each pt
!
    do i=ist,iend
      do j=jst,jend
        xcen=(del(i,j)-del(1,1))/pi180 +1.
        ycen=(tha(i,j)-tha(1,1))/pi180 +1.
        yyo=tha(i,j)
        do ir=1,lgth
          rbd=float(ir)*0.2
          call bound2(dumu,dumv,tanp,rbd,xcen,ycen,yyo)
          tprof(i-ist+1,j-jst+1,ir)=tanp
        end do                     
      end do
    end do
    print *, 'let us see the value of tprof'
    print 333,((tprof(i,1,k),i=4,7),(tprof(i,2,k),i=4,7),(tprof(i,3,k),i=4,7),(tprof(i,4,k),i=4,7),k=1,lgth)
333 format(16f7.1)
!
!  find the first relative maximum along each azimuthal direction
!  find the position of the largest relative maximum
!
    hmax=0.
    do 53 i=1,npts
      do j=1,npts
        do ir=2,lgth-1
          if(tprof(i,j,ir).gt.tprof(i,j,ir-1).and.tprof(i,j,ir).gt.tprof(i,j,ir+1))then
            tmax(i,j)=tprof(i,j,ir)
!          itpos(i,j)=100*(ist+i)+j+jst
!
!       fixed the bug found april 22, 1994.......>>>
!
            itpos(i,j)=100*(ist+i-1)+j+jst-1
            hmax=amax1(tmax(i,j),hmax)
            if(hmax.eq.tmax(i,j))ipos=itpos(i,j)
            if(hmax.eq.tmax(i,j)) rmxpos=float(ir)*0.2
            goto 53
          endif
        end do
        tmax(i,j)=tprof(i,j,1)
        itpos(i,j)=101
        hmax=amax1(hmax,tmax(i,j))
        if(hmax.eq.tmax(i,j)) ipos=itpos(i,j)
      end do
53  continue
    print *, 'hmax, rmxpos, 5x rmxpos are'
    print*,hmax,rmxpos,rmxpos/0.2
    print *, 'tmax is'
    print *,((tmax(i,j),i=1,npts),j=1,npts)
    print *, 'itpos is'
    print *,((itpos(i,j),i=1,npts),j=1,npts)
!
!
!
!  use position of the largest relative maximum as the adjusted
!  center location
!
    print *, 'ipos=',ipos
    ycn=float(mod(ipos,imx))-1.
    xcn=float(ipos/imx)-1.
    ixc=int(xcn)+1
    iyc=int(ycn)+1
    xctest=(xcn+xcorn)*pi180
    yctest=(ycn+ycorn)*pi180
!
!  recompute the tangential wind component based on new center
!
    fact = cos(tha(1,iyc))
    !  print *,'in maxth',ycn,xcn,fact,xctest/pi180,yctest/pi180
    ! print*, 'ixc=',ixc
    ! print*, 'iyc=',iyc
    ! print *, 'ist=', ist
    ! print *, 'jst=', jst
    ! print *, 'iend=', iend
    ! print *, 'jend=', jend
    ! print*, shape(del)
    ! print*, 'del=', del(ixc+1,iyc+1)
    print*,ixc,iyc,del(ixc+1,iyc+1)/pi180,tha(ixc+1,iyc+1)/pi180
    do j=1,jmx
      do i=1,imx
        dx=(del(i,j)-xctest)*fact
        dy=(tha(i,j)-yctest)
        if(dx.ne.0.)theta =atan2(dy,dx)
        if(dy.gt.0..and.dx.eq.0.)theta =90.*pi180
        if(dy.lt.0..and.dx.eq.0.)theta =270.*pi180
        tw(i,j)=-dumu(i,j)*sin(theta) +dumv(i,j)*cos(theta)
        if(i.eq.ixc.and.j.eq.iyc) print*,i,j,dumu(i,j),dumv(i,j),theta,tw(i,j),dx,dy,'check everything'
      end do
    end do
!
    write(77,*)ixc,iyc
    write(77,7700)((tw(i,j)/100.,i=ixc-3,ixc+3),j=iyc+3,iyc-3,-1)
 7700  format(11f5.1)
!
    iflag=0
    hmax=0.
    rmxavg=0.
    do ir=3,iimx
      rxx=float(ir)*deltar
      call calcr(rxx,rtan,xcn,ycn,yctest,dumu,dumv)
      tanavg(ir)=rtan
      if(tanavg(ir-2).lt.tanavg(ir-1).and.tanavg(ir).lt.tanavg(ir-1).and.iflag.eq.0) then
        hmax=tanavg(ir-1)
        rmxavg=rxx-deltar
        iflag=1
      end if
    end do
    print*,'found rmxavg ',rmxavg,hmax
    dxc=xcn
    dyc=ycn
    yold=xcn+xcorn
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
  END SUBROUTINE maxth