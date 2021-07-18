MODULE filter_routines
  INTEGER :: IMX,JMX,nmx
  REAL :: PI, PI180
  PARAMETER (IMX=640, JMX=320, nmx=64)
  PARAMETER  (PI = 4.*ATAN(1.0), PI180 = PI/180.0)
    ! pi180=pi/180.
  CONTAINS
    INCLUDE 'bound.f90'
    INCLUDE 'bound2.f90'
    INCLUDE 'center.f90'
    INCLUDE 'calcr.f90'
    INCLUDE 'calcra.f90'
    INCLUDE 'maxth.f90'
    INCLUDE 'phase.f90'
END MODULE filter_routines