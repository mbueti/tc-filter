  SUBROUTINE findra(dxc,dyc,yc,rmxavg,rfind,tanuv)
!  finds rfind from azimuthally averaged radial profile of tang. wind
!   parameter(imx=640,jmx=320,nmx=64,iimx=110)
    parameter(iimx=110)
    dimension tanuv(iimx)
! common /scale/rmxavg,rfind

    DR=0.1

    dist= rmxavg*1.5
    X1 = 0.0
    RTAN1 = 100000.
    R = 1.0
    r=dist

! only come back to 666 if rtan > 6m/s

666 CONTINUE
    rtan1=100000.
!
!  return to 777 if gradient, dist, or 3m/s are unmet
!
777 continue

    R = R + DR
    irad= int(r/dr)
    rtan= tanuv(irad)

    WRITE(56,*)R,RTAN
    RTAN2 = RTAN
    IF(RTAN.GT.600.)GO TO 666
! PRINT*,'R AND RTAN:  ',R,RTAN

    IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5)GO TO 999
    IF(RTAN2.GE.RTAN1.AND.R.GT.DIST)THEN
      X1 = 1.0
    ENDIF

    IF(RTAN.LT.300..AND.R.GT.DIST)GO TO 999
    WRITE(56,*)R,RTAN

    RTAN1 = RTAN - 4.0
    IF(R.LT.10.8)GO TO 777
999 CONTINUE

    rfind=r

    RETURN
  END SUBROUTINE findra
