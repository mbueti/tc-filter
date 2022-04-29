REAL FUNCTION mergefrac(lon,lat,i,j,rcls,xc,yc)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i,j
  REAL, INTENT(IN) :: xc, yc, rcls
  REAL, DIMENSION(imx), INTENT(IN) :: lon
  REAL, DIMENSION(jmx), INTENT(IN) :: lat
  REAL :: dlon, dlat, dx, dy, dr, b

  dlon = lon(i) - xc
  dlat = lat(j) - yc
  dx = dlon*111.19393
  dy = dlat*111.19393
  dr = SQRT(dx**2+dy**2)
  ! IF (dr.GE.(0.75*rcls)) THEN
    ! mergefrac = MIN(dr/(10.0*rcls), 1.0)
  mergefrac = MAX(MIN(1-(2.5*rcls-dr)/(rcls), 1.0), 0.0)**0.25
  ! mergefrac = MAX(MIN(1.0-(1.1*rcls-dr)/(1.1*rcls), 1.0), 0.0)
  ! ELSE
  !   mergefrac = 1.0
  ! END IF
END FUNCTION mergefrac
