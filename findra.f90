subroutine findra(dxc,dyc,yc,rmxavg,rfind,tanuv)
  implicit none
  
!  finds rfind from azimuthally averaged radial profile of tang. wind

  integer :: irad
  real :: dist, dr, r, rtan, rtan1, rtan2, x1
  real, intent(in) :: rmxavg, dxc, dyc, yc
  real, dimension(iimx), intent(in) :: tanuv
  real, intent(out) :: rfind

  DR=180.0/320.0

  dist= rmxavg*1.5
  X1 = 0.0
  RTAN1 = 100000.
  R = 10.0
  r=dist
  
! only come back to 666 if rtan > 6m/s
666 CONTINUE
  rtan1=100000.
  
!  return to 777 if gradient, dist, or 3m/s are unmet
777 continue

  R = R + DR
  irad= int(r/dr)
  WRITE(64,*) 'irad=',irad
  print *, 'irad=', irad
  rtan= tanuv(irad)
  
  WRITE(56,*)R,RTAN
  RTAN2 = RTAN
  IF(RTAN.GT.600.)GO TO 666
  
  IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5) then
    rfind = r
    return
  end if 

  IF(RTAN2.GE.RTAN1.AND.R.GT.DIST) THEN
    X1 = 1.0
  END IF
  
  IF(RTAN.LT.300..AND.R.GT.DIST) then
    rfind = r
    return
  end if 

  WRITE(56,*)R,RTAN
  
  RTAN1 = RTAN - 4.0
  IF(R.LT.10.8)GO TO 777
!
  rfind=r
!
  return
end subroutine findra
