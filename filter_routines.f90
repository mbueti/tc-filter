MODULE filter_routines
  INTEGER, PARAMETER :: IMX=360, JMX=180, nmx=64
  REAL, PARAMETER :: PI=4.*ATAN(1.0), PI180=PI/180.0

  CONTAINS
    INCLUDE 'mergefrac.f90'
    INCLUDE 'phase.f90'
    INCLUDE 'separ.f90'
    INCLUDE 'bound.f90'
    INCLUDE 'bound2.f90'
    INCLUDE 'center.f90'
    INCLUDE 'maxth.f90'
    INCLUDE 'calcr.f90'
    INCLUDE 'findra.f90'
END MODULE filter_routines
