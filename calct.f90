  SUBROUTINE calct(deltar,xc,yc,yo,xf,rmxlim)
!
!  calculates the radial profile for eight azimuthal angles
!
            ! PARAMETER ( nmx=64)
      !  PARAMETER (IMX=640 , JMX=320, iimx=110)
    parameter(iimx=110)
    DIMENSION xf(imx,jmx)
    dimension idst(nmx),hmax(nmx),rmax(nmx)
!
    common /scale/rmxavg,rfind
    common /pass/xr(iimx,nmx),dist(nmx)
    COMMON  /TOTAL/ DDEL,dtha
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
!
    arad=6371.
    ! print*,'calct',deltar,xc,yc,yo,yold
    fact=cos(yo)
!
!   set the factor to determine lower limit of dist search
!
    bfact = .5
    rdistl =  bfact*rmxavg 
!
    ! print*
    ! print*,'b factor to determine rdistl: ',bfact
!
!
! assume the maximum wind is within rfavg of center (but <10.8 also)
!
    irend =int(rmxlim/deltar)
!
!       
    iravg =int(rdistl/deltar)
!
    ! print*,'lower limit and radius of dist search: ',rdistl,iravg
    ! print*,'upper limit and radius of dist search: ',rmxlim,irend
!
    DX=DDEL*(1./PI180)
    DY=DTHA/PI180
!
!  angle loop
!
    DO 10 I=1,NMX
      THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
!
      do 11 ir=1,iimx
        ro=float(ir)*deltar
        X=(RO*COS(THETA))/(fact) +XC +1.
        Y=(RO*SIN(THETA)) +YC +1.
        IX=INT(X/DX)
        IY=INT(Y/DY)
        IX1=IX+1
        IY1=IY+1
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
        XR(ir,I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY+1) &
                +(1.-Q)*P*XF(IX+1,IY) + P*Q*XF(IX+1,IY+1)
11    continue
!
! find relative max after which ro check begins 
!
      idst = 0
!
      hmax=-10.e10
      do 12 ir=iravg,irend
        hmax(i)=amax1(hmax(i),xr(ir,i))
        if(hmax(i).eq.xr(ir,i))dist(i)=float(ir)*deltar
        if(hmax(i).eq.xr(ir,i))idst(i)=ir
12    continue
!
!  if the max. value is also the endpt it maynot be a relative max
!  check backwards from irend for the last relative max
!
!
      irvgu = iravg+1
      if(irvgu.le.2) then
        irvgu = 3
      endif
!
      ! print*,'irend=',irend
      ! print*,'iravg=',iravg
      ! print*,'idst=',idst
      if(irend.eq.idst(i).or.iravg.eq.idst(i)) then
        do 13 ir=irend,irvgu,-1
          if(ir .lt. 2)then
            ! print*,'bound index i,ir,irend,irvgu,idst ', &
            !        i,ir,irend,irvgu,idst(i)
            ! print*,'xr ', xr(ir-1,i),xr(ir,i)
            ! print*,'other var ', &
            !        iravg,rdistl,rmxavg,deltar,rmxlim,iravg
          endif
          if(xr(ir-1,i).lt.0.) go to 14
          if(xr(ir-1,i).gt.xr(ir,i).and.xr(ir-1,i).ge.xr(ir-2,i)) then
            dist(i)=float(ir-1)*deltar
            ! print*,'readjusting dist',dist(i),hmax(i),xr(ir-1,i),rmxlim
            go to 14
          endif
13      continue
14      continue
      end if
!
!  
      if(idst(i).lt.iravg) then
        ! print*,'lower limit check, dist changed to rmxlim, i is: ',i
        dist(i) = rmxlim
      end if
!
!
10  CONTINUE
!
    do 450 iaa = 1 , nmx
      WRITE(94,*)iaa,dist(iaa) 
450 continue
!
    write(94,*)rdistl,rmxlim
!
!
    dist=dist*1.1
    ! print*,'relative max found '
    print 4400,dist
4400  format(25f4.1)
    write(4,400) (float(ir)*deltar,(xr(ir,i)/100.,i=1,1),ir=1,iimx)
400 format(25f5.1)
    RETURN
  END SUBROUTINE calct