MODULE filter_routines
  INTEGER, PARAMETER :: IMX=360, JMX=180, nmx=64, &
                        KMAX=18, LGI=20, iimx=110
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
    INCLUDE 'calcra.f90'
    INCLUDE 'calct.f90'
    INCLUDE 'amatrix.f90'
    INCLUDE 'rodist.f90'
    INCLUDE 'findra.f90'
    INCLUDE 'wnlsm.f90'
    INCLUDE 'wnnls.f90'
    INCLUDE 'xerror.f90'
    INCLUDE 'xerrwv.f90'
    INCLUDE 'xerprt.f90'
    INCLUDE 'xersav.f90'
    INCLUDE 'xgetua.f90'
    INCLUDE 'xerabt.f90'
    INCLUDE 'xerctl.f90'
    INCLUDE 'fdump.f90'
    INCLUDE 'j4save.f90'
    INCLUDE 'i1mach.f90'
    INCLUDE 'handle_err.f90'
END MODULE filter_routines
