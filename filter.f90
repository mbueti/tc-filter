  PROGRAM FILTER
    USE datetime_module, ONLY: datetime, timedelta
    USE netcdf
    use filter_routines

    implicit none

      INTEGER :: UNCID, VNCID, STATUS, UVARID, VVARID, LATVARID, LONVARID, TIMEVARID, NLON, NLAT, NTIME, TRACKREAD, &
                 TCDATE, TCHOUR, STORMDAYS, STORMHOUR, STORMMIN, TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN, I, &
                 ITC, T, NTC, m, tcdum1, tcdum2, tcdum3, tcdum4, rcls, IBLKMX, JBLKMX, dftx, dfty, iang, &
                 icx, icy, ie, ienv, ifl, ihx, ihy, imax, iimax, irdex, iw, ix, iy, j, jj, jmax, jjmax, jn, js, ngd, ngr, nn1, &
                 nn2, nn3, nn4, nnn, NSTFLG, NTCID, ntr
      REAL :: LONC, LATC, rscale, xold, yold, xv, yv, xcorn, ycorn, xcg, ycg, xcp, ycp, xc, yc, rzr, rzrn, ddel, dtha, dist, &
              dr, dt, dx, dy, dxc, dyc, factr, rtan1, rtan2, r, rad, rrdd, rb, rmxlim, x1, xcgnew, ycgnew, rtan 
      LOGICAL :: END
      REAL, DIMENSION(2) :: TCLAT, TCLON
      INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimIDs
      REAL, DIMENSION(:, :), ALLOCATABLE :: U, V
      REAL, DIMENSION(:), ALLOCATABLE :: LAT, LON, TIME
      CHARACTER(len=220) :: UFILENAME, VFILENAME, TRACKFILENAME
      CHARACTER(len=7) :: UFIELDNAME, VFIELDNAME
      CHARACTER(len=10) :: TCNAME
      CHARACTER(len=4) :: TCORG, TIMESTRING
      CHARACTER(len=3) :: TCID, STORMID
      CHARACTER(len=8) :: DATESTRING
      CHARACTER(len=1), DIMENSION(2) :: TCNS, TCEW
      CHARACTER(len=3), DIMENSION(200) :: TCIDS

      TYPE(datetime), DIMENSION(2) :: TRACKTIME
      TYPE(datetime) :: BASETIME, TCTIME
      TYPE(timedelta) :: DTIME, shdelta
      REAL, DIMENSION(2) :: DTIMES
      real :: sixhours

      PARAMETER (IBLKMX=LGI*IMX+4*KMAX*IMX)
      PARAMETER (JBLKMX=IMX+14*KMAX*IMX)
!

      real, dimension(imx, jmx) :: del, tha, tang, tanp, ds,xf
      real, dimension(imx, jmx, 2) :: dmmm
      real, dimension(nmx) :: disti, rnot
      real, dimension(iimx,nmx) :: rtani
      integer, dimension(imx, lgi) :: typ1c
      integer, dimension(kmax, imx, 4) :: typ2c
      COMMON /FILEC/  TYP1C,TYP2C
      COMMON /WINDS/ DMMM,TANG,DEL,THA,TANP,DS
!
!
      COMMON  /GDINF/ NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX, JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
      COMMON /pass/rtani,disti
      COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
      COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
      COMMON /TOTAL/ DDEL,DTHA
      COMMON /IFACT/NNN,RNOT,RB,IENV
      COMMON /XXX/  XF,XC,YC,DX,DY
      REAL, DIMENSION(IMX) :: FACG1,FACG2,FACT1,FACT2
      REAL, DIMENSION(IBLKMX) :: FILC
      REAL, DIMENSION(JBLKMX) :: DIAG
      REAL, DIMENSION(IMX,JMX) :: XXD,US,VS,UALL,VALL,UP,VP, UFILS,VFILS,UFIL,VFIL,UFILP,VFILP
!
       CHARACTER*4 IBLOCK
!
!      READ IN PART OF GFDL HISTORY TAPE.........>
!
       EQUIVALENCE  (FILC(1),TYP1C(1,1) )

       CALL GETARG(1, UFILENAME)
       CALL GETARG(2, UFIELDNAME)
       CALL GETARG(3, VFILENAME)
       CALL GETARG(4, VFIELDNAME)
       CALL GETARG(5, TRACKFILENAME)

       BASETIME = datetime(1900, 1, 1, 0, 0, 0)
       shdelta = timedelta(0, 6, 0, 0, 0)
       sixhours = shdelta % total_seconds()

       STATUS = NF90_OPEN(path=UFILENAME, mode=NF90_WRITE, ncid=UNCID)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 69)
       STATUS = NF90_INQ_VARID(UNCID, UFIELDNAME, UVARID)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 71)
       STATUS = NF90_INQUIRE_VARIABLE(UNCID, UVARID, dimids = dimIDs)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 73)
       STATUS = NF90_INQUIRE_DIMENSION(UNCID, dimIDs(1), len = NLON)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 75)
       STATUS = NF90_INQUIRE_DIMENSION(UNCID, dimIDs(2), len = NLAT)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 77)
       STATUS = NF90_INQUIRE_DIMENSION(UNCID, dimIDs(3), len = NTIME)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 79)
       STATUS = NF90_OPEN(path=VFILENAME, mode=NF90_WRITE, ncid=VNCID)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 146)
       STATUS = NF90_INQ_VARID(VNCID, VFIELDNAME, VVARID)
       if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 148)

       ALLOCATE(U(NLON, NLAT))
       ALLOCATE(V(NLON, NLAT))
       ALLOCATE(LON(NLON))
       ALLOCATE(LAT(NLAT))
       ALLOCATE(TIME(NTIME))

17    FORMAT(A4, 1x, A3, 1x, A10, I2, I2, I2, 1x, I2, I2, 1x, F3.0, A1, 1x, F4.0, A1, 1x, I4, 1x, &
             I3, 1x, I4, 1x, I4, 1x, I4, 1x, I2, 1x, I3, 1x, I4, 1x, I4, 1x, I4, 1x, I7, 1x, I4, 1x, I4, 1x, I4)
      OPEN(TRACKREAD, file=TRACKFILENAME, status='old')
      END = .FALSE.
      NTCID = 0
      DO WHILE (.NOT.END)
        END = .TRUE.
        READ(TRACKREAD, 17, END=13) TCORG, TCID
        IF (NTCID.EQ.0) THEN
          GO TO 11
        END IF
        DO I=1, NTCID
          IF (TCIDS(I).EQ.TCID) THEN
            GO TO 12
          END IF
        END DO
11      TCIDS(NTCID + 1) = TCID
        NTCID = NTCID + 1
12      END = .FALSE.
13      CONTINUE
      END DO
      CLOSE(TRACKREAD)

      DO ITC=1, NTCID
      STORMID = TCIDS(ITC)
      print *, 'stormid=', stormid
!      STORMID = '09L'
      DO T=1, NTIME
        i=0
        m=0
        tcdum1=0
        tcdum2=0
        tcdum3=0
        tcdum4=0
        rcls=0
        dftx=0
        dfty=0
        iang=0
        icx=0
        icy=0
        ie=0
        ienv=0
        ifl=0
        ihx=0
        ihy=0
        imax=0
        iimax=0
        irdex=0
        iw=0
        ix=0
        iy=0
        j=0
        jj=0
        jmax=0
        jjmax=0
        jn=0
        js=0
        ngd=0
        ngr=0
        nn1=0
        nn2=0
        nn3=0
        nn4=0
        nnn=0
        nstflg=0
        ntr=0

        LONC=0
        LATC=0
        rscale=0
        xold=0
        yold=0
        xv=0
        yv=0
        xcorn=0
        ycorn=0
        xcg=0
        ycg=0
        xcp=0
        ycp=0
        xc=0
        yc=0
        rzr=0
        rzrn=0
        ddel=0
        dtha=0
        dist=0
        dr=0
        dt=0
        dx=0
        dy=0
        dxc=0
        dyc=0
        factr=0
        rtan1=0
        rtan2=0
        r=0
        rad=0
        rrdd=0
        rb=0
        rmxlim=0
        x1=0
        xcgnew=0
        ycgnew=0
        rtan=0

        U = 0
        V = 0
        TCLON = 0
        TCLAT = 0
        del = 0
        tha = 0
        tang = 0
        tanp = 0
        ds = 0
        xf = 0
        dmmm = 0
        disti = 0
        rnot = 0
        rtani = 0
        typ1c = 0
        typ2c = 0

         STATUS = NF90_INQ_VARID(UNCID, UFIELDNAME, UVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 118)
         STATUS = NF90_GET_VAR(UNCID, UVARID, U, START = (/ 1, 1, T /), COUNT = (/ NLON, NLAT, 1 /))
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 122)
         STATUS = NF90_INQ_VARID(UNCID, "lat", LATVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 124)
         STATUS = NF90_GET_VAR(UNCID, LATVARID, LAT)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 126)
         STATUS = NF90_INQ_VARID(UNCID, "lon", LONVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 128)
         STATUS = NF90_GET_VAR(UNCID, LONVARID, LON)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 130)
         STATUS = NF90_INQ_VARID(UNCID, "time", TIMEVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 132)
         STATUS = NF90_GET_VAR(UNCID, TIMEVARID, TIME)

         STORMDAYS = FLOOR(TIME(T))+1
         STORMHOUR = FLOOR(12 * (TIME(T) - STORMDAYS))
         STORMMIN = FLOOR(60 * (12 * (TIME(T) - STORMDAYS) - STORMHOUR))

         DTIME = timedelta(STORMDAYS, STORMHOUR, STORMMIN, 0, 0)
         TCTIME = BASETIME + DTIME


         STATUS = NF90_GET_VAR(VNCID, VVARID, V, START = (/ 1, 1, T /), COUNT = (/ NLON, NLAT, 1 /))
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 152)

!      READ (10)NSTEP,NNEST,IBLOCK
!
!      READ(10)NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,JMAX,
!    *         NSTFLG,ICKX,ICKY,IHX,IHY,DFTX,DFTY

!
!######################################################################
!
!   THE FOLLOWING IS THE CODE TO REMOVE THE HURRICANE COMPONENT OF THE
!   DISTURBANCE FILED FROM A GIVEN WIND FIELD.
!   THE RESULTING FIELD IS CALLED THE "ENVIRONMENTAL FIELD"
!   IN THE GFDL BOGUS SYSTEM,  THE SPECIFIED VORTEX IS ADDED TO THIS
!   RESULTING WIND FIELD TO OPTAIN THE FINAL INITIAL FIELD.
!
!   FIRST, INPUT THE CENTER OF THE STORM CALLED XV ,YV WHICH IS
!   DEFINED IN LINES 86-87 OF THE CODE.
!   THIS WILL CENTER AN 11 X 11 DEGREE GRIDPOINT SQUARE WITHIN WHICH
!   THE GLOBAL CENTER WILL BE DEFINED (IN ROUTINE "CENTER").
!   THIS ROUTINE WILL DETERMINE THE GLOBAL CENTER, FIRST BASED ON THE
!   CENTROID CENTER, AS DEFINED IN:  KURIHARA, BENDER, AND ROSS (1993)
!   MWR, PAGE 2039, EQUATION (6.1).
!
!   NEXT, THE GLOBAL CENTER IS REDEFINED AS THE GRIDPOINT WITH THE
!   LARGEST AZIMUTHALLY-AVERAGED WIND MAXIMUM.
!
!   NEXT, THE PARAMETER DIST(IR) IS DEFINED, AS THE STARTING LOCATION
!   OF THE SEARCH FOR R0 (THE FILTER DOMAIN) IN THE VARIOUS AXIMUTHAL
!   DIRECTIONS.
!   LASTLY, THE NON-HURRICANE WIND WITHIN RO IS DETERMINED FROM AN
!   OPTIMUM INTERPOLATION APPROACH.
!
!
!   IMPORTANT: WE ASSUME THAT THE GRID SPACING OF THE GRID POINTS IS
!   ONE DEGREE (DLON,DLAT).
!
!
!   IN THIS PROGRAM THE INPUT  FIELD WILL BE    fort.10
!   IN THIS PROGRAM THE OUTPUT FIELD WILL BE    fort.46
!
!
!
!#######################################################################
!
!             THE FOLLOWING ARE THE INPUT PARAMETERS:
!             THIS CODE ASSUMES THAT THE LAT,LON GRID SPACING IS
!             ONE DEGREE.................
!
!******************************************************************
!
!
!      DEFINE THE CENTER OF THE STORM:   XV,YV
!                                      (LON,LAT)
!
        OPEN(TRACKREAD, file=TRACKFILENAME, status='old')

        TRACKTIME(1) = BASETIME
        TRACKTIME(2) = BASETIME
        TCID = ''
        END = .FALSE.
        DO WHILE ((STORMID.NE.TCID.OR.(TRACKTIME(2)).LT.TCTIME.OR.(TRACKTIME(1)).EQ.BASETIME).AND..NOT.END)
          END = .TRUE.
          IF (TCID.EQ.STORMID) THEN
            TCLAT(1) = TCLAT(2)
            TCLON(1) = TCLON(2)
            TCNS(1) = TCNS(2)
            TCEW(1) = TCEW(2)
            TRACKTIME(1) = TRACKTIME(2)
          END IF

          READ(TRACKREAD, 17, END=16) TCORG, TCID, TCNAME, TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN, &
                                      TCLAT(2), TCNS(2), TCLON(2), TCEW(2), tcdum1, tcdum2, tcdum3, tcdum4, rcls
          TRACKTIME(2) = datetime(2000+TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN, 0)
          END=.FALSE.
16        CONTINUE
        END DO
        CLOSE(TRACKREAD)

        if (TCLON(1).lt.0) TCLON(1) = TCLON(1) + 360
        if (TCLON(2).lt.0) TCLON(2) = TCLON(2) + 360

        DTIMES(1) = TCTIME % secondsSinceEpoch()-TRACKTIME(1) % secondsSinceEpoch()
        DTIMES(2) = TRACKTIME(2) % secondsSinceEpoch()- TCTIME % secondsSinceEpoch()

!        if (abs(dtimes(2)-dtimes(1)).gt.SIXHOURS) then
!          cycle
!        end if
!     UNCOMMENT WHEN ITERATING TRACK
        IF (TCID.NE.STORMID.OR.TRACKTIME(1).GT.TCTIME.OR.TRACKTIME(2).LT.TCTIME &
            .OR.dtimes(2).gt.SIXHOURS.OR.dtimes(1).gt.SIXHOURS) THEN
           CYCLE
        END IF

        PRINT *, TCID
        PRINT *, STORMID
        print *, TCLON
        print *, TCLAT
        print *, TCTIME % isoformat()
        print *, TRACKTIME(1) % isoformat()
        print *, TRACKTIME(2) % isoformat()

        DO I=1, 2
          IF (TCEW(I).EQ.'W') THEN
            TCLON(I) = -1 * TCLON(I) / 10.0
          ELSE
            TCLON(I) = TCLON(I) / 10.0
          END IF

          IF (TCNS(I).EQ.'S') THEN
            TCLAT(I) = -1 * TCLAT(I) / 10.0
          ELSE
            TCLAT(I) = TCLAT(I) / 10.0
          END IF
        END DO

!       print *, 'dt1=', DTIMES(1) % total_seconds()
!       print *, 'dt2=', DTIMES(2) % total_seconds()

!        print *,'dt1=',(1 - ((DTIMES(1) % total_seconds())
!     *           / (SIXHOURS % total_seconds())))
!        print *,'dt2=',(1 - (DTIMES(2) % total_seconds())
!     *           / (SIXHOURS % total_seconds()))

        print *, 'T=',T
!        print *, 'sixhours=',SIXHOURS%total_seconds()
!        print *, 'tclon=',tclon(1)
!        print *, 'tclon=',tclon(2)
!        print *, 'tclat=',tclat(1)
!        print *, 'tclat=',tclat(2)
!        print *, (DTIMES(1)%total_seconds())/(SIXHOURS%total_seconds())
!        print *, 'TCTIME=', TCTIME%isoformat()
!        print *, 'TRACKTIME(1)=', TRACKTIME(1)%isoformat()
        XV = TCLON(1) * (1 - DTIMES(1)/SIXHOURS)+TCLON(2) * (1 - DTIMES(2)/SIXHOURS)

        YV = TCLAT(1) * (1 - DTIMES(1)/SIXHOURS) +TCLAT(2) * (1 - DTIMES(2)/SIXHOURS)

        print *, 'tclon=', tclon(1), tclon(2)
        print *, 'tclat=', tclat(1), tclat(2)
        PRINT *, 'XV=', XV
        PRINT *, 'YV=', YV
!
        IF (XV.LT.0) THEN
          XV = 360. + XV
        ENDIF
!
!****************SET UP THE FILTER STRENGTH YOU WANT***************
!
!   THIS IS THE FIRST FILTER, WHICH SEPERATES THE DISTURBANCE WIND
!   FIELD FROM THE BASIC FLOW.  THE BASIC FLOW WILL BE DEFINED AS
!   (US, VS).
!
!   SEE THE SUBROUTINE PHASE FOR DETAILS.
!
!    IFL=1   IS THE WEAK    FILTER
!    IFL=2   IS THE REGULAR FILTER *** CURRENTLY IN USE
!    IFL=3   IS THE STRONG  FILTER
!    IFL=4   IS VERY STRONG FILTER
!
!
!    FILTER IS DEFINED IN MWR PAPER OF KURIHARA, ET.ALL, 1990:
!
         IFL = 2
!
!
!**********************************************************
!
!
!      DEFINE THE VERTICAL LEVEL OF YOUR HISTORY TAPE WHICH
!      YOU WANT TO FILTER.
!      LEVEL 10 OF THE GFDL MODEL IS AT ABOUT 500 hPa HEIGHT.
!
!      KFIL = 1
!
!
!      DEFINE THE VERTICAL LEVEL OF THE GFDL HISTORY TAPE WHICH
!      IS USED TO DEFINE THE FILTER DOMAIN RO (KTOP).
!      THIS SHOULD BE THE MODEL LEVEL NEAR 850 hPa. FOR THE GFDL
!      MODEL THIS IS LEVEL 14.
!
!      KTOP = 1
!
!**************************************************************
!
!
!      INPUT THE HISTORY TAPE................>>>>>
!
!      THESE PARAMETERS WILL HAVE TO BE SPECIFIED FOR THE PROGRAM TO
!      WORK !!!!!!!!!
!
!
!      DDEL, DTHA IS THE GRID SPACING (SHOULD BE 1 DEGREE FOR THE FILTER
!      CHARACTERISTICS TO PROPERLY WORK)
!      THESE ARE INPUTED IN RADIANS !!!!
!
!                  NEXT:
!      DEL : IS THE 2-D ARRAY CONTAINING LONGTITUDE IN RADIANS OF THE GRID
!      THA : IS THE 2-D ARRAY CONTAINING LATITUDE   IN RADIANS OF THE GRID
!
!      DS :  IS THE AREA CONTAINED WITH ONE GRID POINT
!            THIS IS USED IN ROUTINE ....CENTER.... FOR THE FIRST GUESS
!            OF THE GLOBAL CENTER (WHICH IS A CENTROID CALCULATION).
!
!
!       U :  IS THE U COMPONENT OF THE WIND USED TO DEFINE RO
!       V :  IS THE V COMPONENT OF THE WIND USED TO DEFINE RO
!       UFIL :      U COMPONENT OF THE WIND USED TO BE FILTERED
!       VFIL :      V COMPONENT OF THE WIND USED TO BE FILTERED
!
!**************************************************************
!
        RAD = 6.371E3
!
        DR = 180.0/320.0
!
        XC = XV*PI180
        YC = YV*PI180
!
!
!
        J = 0
        DO 15 JJ = 1, JMX
          J = J + 1
!
!
!      READ IN THE JMAX ROWS OF THE GFDL HISTORY TAPE TO OBTAIN:
!      DS, THA, U AND V...............>


!
!
!      READ(10)(DIAG(JQ)  , JQ = 1,JBLKMX)
!      READ(10)(FILC(I)  ,   I = 1,IBLKMX)
!
!
!      DDEL AND DTHA WILL BE SET TO ONE DEGREE  (IN RADIANS)
!      (.017453.....)
!
!      IF(J.EQ.5)THEN
!      DDEL = TYP1C(5,2)
!      DTHA = TYP1C(5,3)
!
!      ENDIF
!
          DDEL = PI180 * (LON(2) - LON(1))
          DTHA = PI180 * (LAT(2) - LAT(1))
          DO 20 I = 1, IMX
            DS(I,J) = (PI180**2)*(LON(2) - LON(1))*(LAT(2) - LAT(1))
            DEL(I,J) = PI180 * LON(I)
            THA(I,J) = PI180 * LAT(J)
            UALL(I,J)   = U(I, J)
            VALL(I,J)   = V(I, J)
            UFIL(I,J)   = U(I, J)
            VFIL(I,J)   = V(I, J)
20        CONTINUE
15      CONTINUE
!
!      REWIND 10
!
      XCORN = DEL(1,1) / PI180
      YCORN = THA(1,1) / PI180
!
        PRINT*
        PRINT*
        PRINT*
        PRINT*
        PRINT*, 'THIS IS (XCORN,YCORN): ', XCORN,YCORN
        PRINT*
        PRINT*
        PRINT*
!
!
!       OBTAIN THE BASIC FLOW FIELD FOR THE WIND (U,V)
!
        CALL PHASE(IFL,U,V,IMX,JMX,US,VS)
        CALL PHASE(IFL,UFIL,VFIL,IMX,JMX,UFILS,VFILS)
!
!#############################################################
!
!
! NEXT FIND THE CENTER POSITION OF THE GLOBAL VORTEX AND RNOT
!
!
!       OBTAIN THE TOTAL DISTURBANCE FIELD:
!
        DO 995 J = 1, JMX
          DO 995 I = 1, IMX
            UP(I,J) = UALL(I,J) - US(I,J)
            VP(I,J) = VALL(I,J) - VS(I,J)
            UFILP(I,J) = UFIL(I,J) - UFILS(I,J)
            VFILP(I,J) = VFIL(I,J) - VFILS(I,J)
995     CONTINUE
!
!
!    FIRST FIND THE CENTER POSITION OF THE GLOBAL VORTEX
!    THIS IS THE CENTROID CALCULATION...........
!
        CALL CENTER(UP,VP,XCG,YCG)
!
        XOLD = XCG
        YOLD = YCG
!
!
!
!  adjust the center position of the global vortex
!
        PRINT*,'before maxth', xold,yold
        write(67,*)xold/pi180-360,yold/pi180
        CALL MAXTH(up,vp,xcgnew,ycgnew,rmxlim,tanp)
        xold = xcgnew+xcorn
        yold = ycgnew+ycorn
        dist = rmxlim
!
        print *, 'after maxth', xold, yold
!        print *, 'dist=', dist
        xcp = xold*pi180
        ycp = yold*pi180
        write(4,*)xold,yold,xcp,ycp
        write(66,*)xold-360,yold,xcp,ycp
!
!
!
!    NOW DETERMINE THE RNOT OF THE GLOBAL VORTEX
!
!
!  loop over nmx azimuthal directions
!   first compute the radial profiles of tangential wind for the
!   NMX  azimuthal angles
!
        dxc=xold -xcorn
        dyc=yold -ycorn
        print*,'before calct',dxc,dyc,ycp
!
!
!    CALCULATE THE RADIAL PROFILE OF TANGENTIAL WIND FOR 24 AXUMUTHAL
!    ANGELS
!
        call calct(dr,dxc,dyc,ycp,tanp,rmxlim)
!
!
!        DO iang=1,nmx
!          X1 = 0.0
!          RTAN1 = 100000.
!          R = 1.0
!          dist=disti(iang)
!        end do
        DO 10 iang=1,nmx
          X1 = 0.0
          RTAN1 = 100000.
          R = 1.0
          dist=disti(iang)
!          print *, 'dist=',dist
!
!  only return to 666 if rtan > 6m/s
!
666       CONTINUE
          Rtan1=100000.
!
!  return to 777 if dist or grad condition not met
!
777       continue
!
         CALL CALCRa(R,RTAN,iang,dist)
          irdex=int(r/dr)
          rtan = rtani(irdex,iang)
          R = R + DR
!         WRITE(56,*)R,RTAN
          RTAN2 = RTAN
          IF(RTAN.GT.600.)GO TO 666
          IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5)GO TO 999
          IF(RTAN2.GE.RTAN1.AND.R.GT.DIST)THEN
            X1 = 1.0
          ENDIF
!
          IF(RTAN.LT.300..AND.R.GT.DIST)GO TO 999
!         WRITE(56,*)R,RTAN
          RTAN1 = RTAN -4.0
          IF(R.LT.10.8)GO TO 777
999       CONTINUE
!         PRINT*
!
!
          IF (X1.EQ.1.0) THEN
            RNOT(iang) = (R-.1)/.1
          ELSE
            RNOT(iang) =  R/.1
          ENDIF
!
           rscale = 6-5*tanh(exp(1.0)*(rcls+200)/1200)
           RZR = rscale*float(rcls)/111.19393
!          rzr=dist
          m = 1
1999      CONTINUE
!          if (m.gt.300) THEN
!            goto 6666
!          end if
          m = m+1
!          print *, 'm=',m
!          print *, 'dist=',dist
!          print *, 'iang=',iang
          RZR = RZR + DR
!          print *, 'rzr=', rzr
!          print *, 'rtan=', rtan
          CALL CALCRa(RZR,RTAN,iang,dist)
!          print *, 'dist=',dist
!          print *, 'rzr=',rzr
!          print *, 'dr=', dr
          irdex=int(rzr/dr)
          rtan = rtani(irdex,iang)
          IF (RTAN.GT.0.0) GO TO 1999
6666      CONTINUE
!
        !   RZR = RZR*1.5
          RZRN = RZR
          RZR = RZR*111.19493
!
!         PRINT*,'RO WILL BE RTAN TIMES SOME CONSTANT'
!
!
!
          RRDD = RNOT(iang)*111.19493
!
          IF(RRDD.GT.RZR)THEN
            PRINT*
            PRINT*,'RNOT HAS A NEGATIVE TANGENTIAL COMPONENT'
            PRINT*,'RNOT WAS DEFINED AS:  ',RRDD
            PRINT*,'RNOT WILL BE MODIFIED'
            RRDD = RZR
            RNOT(iang) = RZRN
          ENDIF
!
          PRINT*,'THIS IS RO IN KM:  ',RRDD
10      CONTINUE
!
!
!
!
!
!################################################################
!
!
!      DO THE OPTIMUM INTERPOLATION......>>>>>
!
        print *, 'entering rodist'
        call rodist
        print *, 'exiting rodist'
!
!
!     CREATE MATRIX  [A]  CONTAINING THE DISTANCE-RELATED CORRELATIONS
!     BETWEEN THE 24 BOUNDARY POINTS JUST OUTSIDE OF THE FILTER DOMAIN
!     RNOT......
!
!
        call amatrix
!
!
!     SEPERATE THE DISTURBANCE INTO THE HURRICNAE AND NON-HURRICNAE
!     COMPONENTS
!
!
        print *, TCTIME % isoformat()
        print *, "RZR=",rzr
        LONC = XC/PI180
        LATC = YC/PI180
         DO 880 J = 1, JMX
           DO 880 I = 1, IMX
             XXD(I,J) = UFIL(I,J) - UFILS(I,J)
880      CONTINUE
         print *, 'storm=', STORMID
         print *, 'T=', T
         CALL SEPAR(XXD)
         DO 890 J = 1, JMX
           DO 890 I = 1, IMX
             UFIL(I,J) = UFILS(I,J)+mergefrac(lon,lat,i,j,rzr,lonc,latc)*XXD(I,J)
890      CONTINUE
         DO 980 J = 1 , JMX
           DO 980 I = 1 , IMX
             XXD(I,J) = VFIL(I,J) - VFILS(I,J)
980      CONTINUE
         CALL SEPAR(XXD)
         DO 990 J = 1 , JMX
           DO 990 I = 1 , IMX
             VFIL(I,J) = VFILS(I,J)+mergefrac(lon,lat,i,j,rzr,lonc,latc)*XXD(I,J)
990      CONTINUE

!
!
!      PUT THE ENVIRONMENTAL WINDS INTO THE GFDL HISTROY TAPE
!
        STATUS = NF90_PUT_VAR(UNCID, UVARID, UFIL, START = (/ 1, 1, T /), COUNT = (/ NLON, NLAT, 1 /))
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 611)
        STATUS = NF90_PUT_VAR(VNCID, VVARID, VFIL, START = (/ 1, 1, T /), COUNT = (/ NLON, NLAT, 1 /))
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 615)
       END DO
       END DO
        STATUS = NF90_CLOSE(UNCID)
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 619)
        STATUS = NF90_CLOSE(VNCID)
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 621)
        STOP
END PROGRAM FILTER
