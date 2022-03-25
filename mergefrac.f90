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
  IF (dr.GE.(0.5*rcls)) THEN
    mergefrac = MIN(dr/(10.0*rcls), 1.0)
    mergefrac = MIN(1.0-(5.0*rcls-dr)/(5.0*rcls), 1.0)
  ELSE
    mergefrac = 1.0
  END IF
END FUNCTION mergefrac
