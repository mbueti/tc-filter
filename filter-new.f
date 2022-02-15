      PROGRAM FILTER
      USE netcdf
      USE datetime_module, ONLY: datetime, timedelta

        implicit none

      INTEGER :: UNCID, VNCID, STATUS, UVARID, VVARID,
     *           LATVARID, LONVARID, TIMEVARID,
     *           NLON, NLAT, NTIME, IMX, JMX, TRACKREAD,
     *           TCDATE, TCHOUR, STORMDAYS, STORMHOUR, STORMMIN,
     *           TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR,
     *           TRACKMIN, I, ITC, T, NTC, m, tcdum1, tcdum2,
     *           tcdum3, tcdum4, rcls, nmx, kmax, lgi, iimx,
     *           IBLKMX, JBLKMX, dftx, dfty, iang, icx, icy,
     *           ie, ienv, ifl, ihx, ihy, imax, iimax, irdex,
     *           iw, ix, iy, j, jj, jmax, jjmax, jn, js,
     *           ngd, ngr, nn1, nn2, nn3, nn4, nnn, NSTFLG,
     *           NTCID, ntr
      REAL :: LONC, LATC, rscale, xold, yold, xv, yv,
     *        xcorn, ycorn, xcg, ycg, xcp, ycp, xc, yc,
     *        rzr, rzrn, ddel, dtha, dist, dr, dt, dx, dy,
     *        dxc, dyc, factr, rtan1, rtan2, r, rad, rrdd,
     *        pi, pi180, rb, rmxlim, x1, xcgnew, ycgnew,
     *        rtan
      LOGICAL :: END
      REAL, DIMENSION(2) :: TCLAT, TCLON
      INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimIDs
      REAL, DIMENSION(:, :), ALLOCATABLE :: U, V
      REAL, DIMENSION(:), ALLOCATABLE :: LAT, LON, TIME
C      real, dimension(:,:) :: del
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
      TYPE(timedelta) :: DTIME, SIXHOURS
      TYPE(timedelta), DIMENSION(2) :: DTIMES

      PARAMETER (IMX=360, JMX=180, nmx=64)
      PARAMETER (KMAX=18,LGI=20 ,iimx=100)
      PARAMETER (IBLKMX=LGI*IMX+4*KMAX*IMX)
      PARAMETER (JBLKMX=IMX+14*KMAX*IMX)
C

      real, dimension(imx, jmx) :: del, tha, tang, tanp, ds,
     *                             xf
      real, dimension(imx, jmx, 2) :: dmmm
      real, dimension(nmx) :: disti, rnot
      real, dimension(iimx,nmx) :: rtani
      integer, dimension(imx, lgi) :: typ1c
      integer, dimension(kmax, imx, 4) :: typ2c
      COMMON /FILEC/  TYP1C,TYP2C
      COMMON /WINDS/ DMMM,TANG,
     *               DEL,THA,TANP,DS
C
C
      COMMON  /GDINF/ NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,
     *                JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
      COMMON /pass/rtani,disti
      COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
      COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
      COMMON /TOTAL/ DDEL,DTHA
      COMMON /IFACT/NNN,RNOT,RB,IENV
      COMMON /XXX/  XF,XC,YC,DX,DY
      REAL, DIMENSION(IMX) :: FACG1,FACG2,FACT1,FACT2
      REAL, DIMENSION(IBLKMX) :: FILC
      REAL, DIMENSION(JBLKMX) :: DIAG
      REAL, DIMENSION(IMX,JMX) :: XXD,US,VS,UALL,VALL,UP,VP,
     *                            UFILS,VFILS,UFIL,VFIL,UFILP,VFILP
C
       CHARACTER*4 IBLOCK
C
C      READ IN PART OF GFDL HISTORY TAPE.........>
C
       EQUIVALENCE  (FILC(1),TYP1C(1,1) )

       CALL GETARG(1, UFILENAME)
       CALL GETARG(2, UFIELDNAME)
       CALL GETARG(3, VFILENAME)
       CALL GETARG(4, VFIELDNAME)
       CALL GETARG(5, TRACKFILENAME)

       BASETIME = datetime(1900, 1, 1, 0, 0, 0)
       SIXHOURS = timedelta(0, 6, 0, 0, 0)

C       PRINT *, "FILENAME=", UFILENAME
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

       ALLOCATE(U(NLON, NLAT))
       ALLOCATE(V(NLON, NLAT))
       ALLOCATE(LON(NLON))
       ALLOCATE(LAT(NLAT))
       ALLOCATE(TIME(NTIME))

C17     FORMAT(A4,A4,A10,1x,I4,I2,I2,1x,I2,I2,1x,F3.0,A1,1x,
C     *        F4.0,A1,1x,I3,1x,I3,3I5,
C     *        1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,
C     *        1x,I4,1x,I4,1x,I4)
17    FORMAT(A4, 1x, A3, 1x, A10, I2, I2, I2, 1x, I2, I2,
     *      1x, F3.0, A1, 1x, F4.0, A1, 1x, I4, 1x I3, 1x, I4,
     *      1x, I4, 1x, I4, 1x, I2, 1x, I3, 1x, I4, 1x, I4,
     *      1x, I4, 1x, I7, 1x, I4, 1x, I4, 1x, I4)
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

C     DO ITC=1, NTCID
c     STORMID = TCIDS(ITC)
      STORMID = '09L'
      DO T=1, NTIME
         STATUS = NF90_INQ_VARID(UNCID, UFIELDNAME, UVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 118)
         STATUS = NF90_GET_VAR(UNCID, UVARID, U,
     *                       START = (/ 1, 1, T /),
     *                       COUNT = (/ NLON, NLAT, 1 /))
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
         STORMMIN = FLOOR(60 * (12 * (TIME(T) - STORMDAYS)
     *                     - STORMHOUR))

         DTIME = timedelta(STORMDAYS, STORMHOUR,
     *                   STORMMIN, 0, 0)
         TCTIME = BASETIME + DTIME


         STATUS = NF90_open(VFILENAME, NF90_WRITE, VNCID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 146)
         STATUS = NF90_INQ_VARID(VNCID, VFIELDNAME, VVARID)
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 148)
         STATUS = NF90_GET_VAR(VNCID, VVARID, V,
     *                       START = (/ 1, 1, T /),
     *                       COUNT = (/ NLON, NLAT, 1 /))
         if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 152)

C      READ (10)NSTEP,NNEST,IBLOCK
C
C      READ(10)NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,JMAX,
C    *         NSTFLG,ICKX,ICKY,IHX,IHY,DFTX,DFTY

C
C######################################################################
C
C   THE FOLLOWING IS THE CODE TO REMOVE THE HURRICANE COMPONENT OF THE
C   DISTURBANCE FILED FROM A GIVEN WIND FIELD.
C   THE RESULTING FIELD IS CALLED THE "ENVIRONMENTAL FIELD"
C   IN THE GFDL BOGUS SYSTEM,  THE SPECIFIED VORTEX IS ADDED TO THIS
C   RESULTING WIND FIELD TO OPTAIN THE FINAL INITIAL FIELD.
C
C   FIRST, INPUT THE CENTER OF THE STORM CALLED XV ,YV WHICH IS
C   DEFINED IN LINES 86-87 OF THE CODE.
C   THIS WILL CENTER AN 11 X 11 DEGREE GRIDPOINT SQUARE WITHIN WHICH
C   THE GLOBAL CENTER WILL BE DEFINED (IN ROUTINE "CENTER").
C   THIS ROUTINE WILL DETERMINE THE GLOBAL CENTER, FIRST BASED ON THE
C   CENTROID CENTER, AS DEFINED IN:  KURIHARA, BENDER, AND ROSS (1993)
C   MWR, PAGE 2039, EQUATION (6.1).
C
C   NEXT, THE GLOBAL CENTER IS REDEFINED AS THE GRIDPOINT WITH THE
C   LARGEST AZIMUTHALLY-AVERAGED WIND MAXIMUM.
C
C   NEXT, THE PARAMETER DIST(IR) IS DEFINED, AS THE STARTING LOCATION
C   OF THE SEARCH FOR R0 (THE FILTER DOMAIN) IN THE VARIOUS AXIMUTHAL
C   DIRECTIONS.
C   LASTLY, THE NON-HURRICANE WIND WITHIN RO IS DETERMINED FROM AN
C   OPTIMUM INTERPOLATION APPROACH.
C
C
C   IMPORTANT: WE ASSUME THAT THE GRID SPACING OF THE GRID POINTS IS
C   ONE DEGREE (DLON,DLAT).
C
C
C   IN THIS PROGRAM THE INPUT  FIELD WILL BE    fort.10
C   IN THIS PROGRAM THE OUTPUT FIELD WILL BE    fort.46
C
C
C
C#######################################################################
C
C             THE FOLLOWING ARE THE INPUT PARAMETERS:
C             THIS CODE ASSUMES THAT THE LAT,LON GRID SPACING IS
C             ONE DEGREE.................
C
C******************************************************************
C
C
C      DEFINE THE CENTER OF THE STORM:   XV,YV
C                                      (LON,LAT)
C
        OPEN(TRACKREAD, file=TRACKFILENAME, status='old')

        TRACKTIME(1) = BASETIME
        TRACKTIME(2) = BASETIME
        TCID = ''
        END = .FALSE.
        DO WHILE (
     *    (STORMID.NE.TCID.OR.
     *    (TRACKTIME(2)).LT.TCTIME.OR.
     *    (TRACKTIME(1)).EQ.BASETIME).AND.
     *    .NOT.END
     *  )
          END = .TRUE.
          IF (TCID.EQ.STORMID) THEN
            TCLAT(1) = TCLAT(2)
            TCLON(1) = TCLON(2)
            TCNS(1) = TCNS(2)
            TCEW(1) = TCEW(2)
            TRACKTIME(1) = TRACKTIME(2)
          END IF

          READ(TRACKREAD, 17, END=16) TCORG, TCID, TCNAME, TRACKYEAR,
     *                      TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN,
     *                      TCLAT(2), TCNS(2), TCLON(2), TCEW(2),
     *                      tcdum1, tcdum2, tcdum3, tcdum4, rcls
          TRACKTIME(2) = datetime(2000+TRACKYEAR, TRACKMONTH,
     *                       TRACKDAY, TRACKHOUR, TRACKMIN, 0)
          END=.FALSE.
16        CONTINUE
        END DO
        CLOSE(TRACKREAD)

        if (TCLON(1).lt.0) TCLON(1) = TCLON(1) + 360
        if (TCLON(2).lt.0) TCLON(2) = TCLON(2) + 360

        DTIMES(1) = TCTIME - TRACKTIME(1)
        DTIMES(2) = TRACKTIME(2) - TCTIME

C     UNCOMMENT WHEN ITERATING TRACK
        IF (TCID.NE.STORMID.OR.
     *      TRACKTIME(1).GT.TCTIME.OR.
     *      TRACKTIME(2).LT.TCTIME
     *  ) THEN
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

C       print *, 'dt1=', DTIMES(1) % total_seconds()
C       print *, 'dt2=', DTIMES(2) % total_seconds()

C        print *,'dt1=',(1 - ((DTIMES(1) % total_seconds())
C     *           / (SIXHOURS % total_seconds())))
C        print *,'dt2=',(1 - (DTIMES(2) % total_seconds())
C     *           / (SIXHOURS % total_seconds()))

        print *, 'T=',T
C        print *, 'sixhours=',SIXHOURS%total_seconds()
C        print *, 'tclon=',tclon(1)
C        print *, 'tclon=',tclon(2)
C        print *, 'tclat=',tclat(1)
C        print *, 'tclat=',tclat(2)
C        print *, (DTIMES(1)%total_seconds())/(SIXHOURS%total_seconds())
C        print *, 'TCTIME=', TCTIME%isoformat()
C        print *, 'TRACKTIME(1)=', TRACKTIME(1)%isoformat()
        XV = TCLON(1) * (1 - (DTIMES(1) % total_seconds())
     *                     / (SIXHOURS % total_seconds())) +
     *       TCLON(2) * (1 - (DTIMES(2) % total_seconds())
     *                     / (SIXHOURS % total_seconds()))

        YV = TCLAT(1) * (1 - (DTIMES(1) % total_seconds())
     *                     / (SIXHOURS % total_seconds())) +
     *       TCLAT(2) * (1 - (DTIMES(2) % total_seconds())
     *                     / (SIXHOURS % total_seconds()))

        print *, 'dtime=', dtimes(1)%total_seconds(),
     *   dtimes(2)%total_seconds()
        print *, 'tclon=', tclon(1), tclon(2)
        print *, 'tclat=', tclat(1), tclat(2)
        PRINT *, 'XV=', XV
        PRINT *, 'YV=', YV
C
        IF (XV.LT.0) THEN
          XV = 360. + XV
        ENDIF
C
CC****************SET UP THE FILTER STRENGTH YOU WANT***************
C
C   THIS IS THE FIRST FILTER, WHICH SEPERATES THE DISTURBANCE WIND
C   FIELD FROM THE BASIC FLOW.  THE BASIC FLOW WILL BE DEFINED AS
C   (US, VS).
C
CC   SEE THE SUBROUTINE PHASE FOR DETAILS.
CC
CC    IFL=1   IS THE WEAK    FILTER
CC    IFL=2   IS THE REGULAR FILTER *** CURRENTLY IN USE
CC    IFL=3   IS THE STRONG  FILTER
CC    IFL=4   IS VERY STRONG FILTER
CC
CC
CC    FILTER IS DEFINED IN MWR PAPER OF KURIHARA, ET.ALL, 1990:
CC
         IFL = 1
CC
CC
CC**********************************************************
CC
C
C      DEFINE THE VERTICAL LEVEL OF YOUR HISTORY TAPE WHICH
C      YOU WANT TO FILTER.
C      LEVEL 10 OF THE GFDL MODEL IS AT ABOUT 500 hPa HEIGHT.
C
C      KFIL = 1
C
C
C      DEFINE THE VERTICAL LEVEL OF THE GFDL HISTORY TAPE WHICH
C      IS USED TO DEFINE THE FILTER DOMAIN RO (KTOP).
C      THIS SHOULD BE THE MODEL LEVEL NEAR 850 hPa. FOR THE GFDL
C      MODEL THIS IS LEVEL 14.
C
C      KTOP = 1
C
C**************************************************************
C
C
C      INPUT THE HISTORY TAPE................>>>>>
C
C      THESE PARAMETERS WILL HAVE TO BE SPECIFIED FOR THE PROGRAM TO
C      WORK !!!!!!!!!
C
C
C      DDEL, DTHA IS THE GRID SPACING (SHOULD BE 1 DEGREE FOR THE FILTER
C      CHARACTERISTICS TO PROPERLY WORK)
C      THESE ARE INPUTED IN RADIANS !!!!
C
C                  NEXT:
C      DEL : IS THE 2-D ARRAY CONTAINING LONGTITUDE IN RADIANS OF THE GRID
C      THA : IS THE 2-D ARRAY CONTAINING LATITUDE   IN RADIANS OF THE GRID
C
C      DS :  IS THE AREA CONTAINED WITH ONE GRID POINT
C            THIS IS USED IN ROUTINE ....CENTER.... FOR THE FIRST GUESS
C            OF THE GLOBAL CENTER (WHICH IS A CENTROID CALCULATION).
C
C
C       U :  IS THE U COMPONENT OF THE WIND USED TO DEFINE RO
C       V :  IS THE V COMPONENT OF THE WIND USED TO DEFINE RO
C       UFIL :      U COMPONENT OF THE WIND USED TO BE FILTERED
C       VFIL :      V COMPONENT OF THE WIND USED TO BE FILTERED
C
C**************************************************************
C
        RAD = 6.371E3
        PI = 4.*ATAN(1.0)
        PI180 = PI/180.
C
        DR = 1.0
C
        XC = XV*PI180
        YC = YV*PI180
C
C
C
        J = 0
        DO 15 JJ = 1, JMX
          J = J + 1
CC
C
C      READ IN THE JMAX ROWS OF THE GFDL HISTORY TAPE TO OBTAIN:
C      DS, THA, U AND V...............>


C
C
C      READ(10)(DIAG(JQ)  , JQ = 1,JBLKMX)
C      READ(10)(FILC(I)  ,   I = 1,IBLKMX)
C
C
C      DDEL AND DTHA WILL BE SET TO ONE DEGREE  (IN RADIANS)
C      (.017453.....)
C
C      IF(J.EQ.5)THEN
C      DDEL = TYP1C(5,2)
C      DTHA = TYP1C(5,3)
C
C      ENDIF
C
          DDEL = PI180 * (LON(2) - LON(1))
          DTHA = PI180 * (LAT(2) - LAT(1))
          DO 20 I = 1, IMX
            DS(I,J) = (PI180 ** 2) * (LON(2) - LON(1))
     *                * (LAT(2) - LAT(1))
            DEL(I,J) = PI180 * LON(I)
            THA(I,J) = PI180 * LAT(J)
            UALL(I,J)   = U(I, J)
            VALL(I,J)   = V(I, J)
            UFIL(I,J)   = U(I, J)
            VFIL(I,J)   = V(I, J)
20        CONTINUE
15      CONTINUE
C
C      REWIND 10
C
      XCORN = DEL(1,1) / PI180
      YCORN = THA(1,1) / PI180
CC
        PRINT*
        PRINT*
        PRINT*
        PRINT*
        PRINT*, 'THIS IS (XCORN,YCORN): ', XCORN,YCORN
        PRINT*
        PRINT*
        PRINT*
C
C
C       OBTAIN THE BASIC FLOW FIELD FOR THE WIND (U,V)
C
        CALL PHASE(IFL,U,V,IMX,JMX,US,VS)
        CALL PHASE(IFL,UFIL,VFIL,IMX,JMX,UFILS,VFILS)
C
C#############################################################
C
C
C NEXT FIND THE CENTER POSITION OF THE GLOBAL VORTEX AND RNOT
C
C
C       OBTAIN THE TOTAL DISTURBANCE FIELD:
C
        DO 995 J = 1, JMX
          DO 995 I = 1, IMX
            UP(I,J) = UALL(I,J) - US(I,J)
            VP(I,J) = VALL(I,J) - VS(I,J)
            UFILP(I,J) = UFIL(I,J) - UFILS(I,J)
            VFILP(I,J) = VFIL(I,J) - VFILS(I,J)
995     CONTINUE
C
C
C    FIRST FIND THE CENTER POSITION OF THE GLOBAL VORTEX
C    THIS IS THE CENTROID CALCULATION...........
C
        CALL CENTER(UP,VP,XCG,YCG)
C
        XOLD = XCG
        YOLD = YCG
C
C
c
c  adjust the center position of the global vortex
c
        PRINT*,'before maxth', xold,yold
        write(67,*)xold/pi180-360,yold/pi180
        CALL MAXTH(up,vp,xcgnew,ycgnew,rmxlim,tanp)
        xold = xcgnew+xcorn
        yold = ycgnew+ycorn
        dist = rmxlim
c
        print *, 'after maxth', xold, yold
C        print *, 'dist=', dist
        xcp = xold*pi180
        ycp = yold*pi180
        write(4,*)xold,yold,xcp,ycp
        write(66,*)xold-360,yold,xcp,ycp
C
C
C
C    NOW DETERMINE THE RNOT OF THE GLOBAL VORTEX
C
C
c  loop over nmx azimuthal directions
c   first compute the radial profiles of tangential wind for the
c   NMX  azimuthal angles
c
        dxc=xold -xcorn
        dyc=yold -ycorn
        print*,'before calct',dxc,dyc,ycp
C
C
C    CALCULATE THE RADIAL PROFILE OF TANGENTIAL WIND FOR 24 AXUMUTHAL
C    ANGELS
C
        call calct(dr,dxc,dyc,ycp,tanp,rmxlim)
C
C
C        DO iang=1,nmx
C          X1 = 0.0
C          RTAN1 = 100000.
C          R = 1.0
C          dist=disti(iang)
C        end do
        DO 10 iang=1,nmx
          X1 = 0.0
          RTAN1 = 100000.
          R = 1.0
          dist=disti(iang)
C          print *, 'dist=',dist
c
c  only return to 666 if rtan > 6m/s
c
666       CONTINUE
          Rtan1=100000.
c
c  return to 777 if dist or grad condition not met
c
777       continue
c
         CALL CALCRa(R,RTAN,iang,dist)
          irdex=int(r/dr)
          rtan = rtani(irdex,iang)
          R = R + DR
c         WRITE(56,*)R,RTAN
          RTAN2 = RTAN
          IF(RTAN.GT.600.)GO TO 666
          IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5)GO TO 999
          IF(RTAN2.GE.RTAN1.AND.R.GT.DIST)THEN
            X1 = 1.0
          ENDIF
C
          IF(RTAN.LT.300..AND.R.GT.DIST)GO TO 999
c         WRITE(56,*)R,RTAN
          RTAN1 = RTAN -4.0
          IF(R.LT.10.8)GO TO 777
999       CONTINUE
c         PRINT*
c
C
          IF (X1.EQ.1.0) THEN
            RNOT(iang) = (R-.1)/.1
          ELSE
            RNOT(iang) =  R/.1
          ENDIF
C
           rscale = 2.0 - tanh(exp(1.0)*(rcls-100)/1500)
           RZR = rscale*float(rcls)/111.19393
C          rzr=dist
          m = 1
1999      CONTINUE
C          if (m.gt.300) THEN
C            goto 6666
C          end if
          m = m+1
C          print *, 'm=',m
C          print *, 'dist=',dist
C          print *, 'iang=',iang
          RZR = RZR + DR
C          print *, 'rzr=', rzr
C          print *, 'rtan=', rtan
          CALL CALCRa(RZR,RTAN,iang,dist)
C          print *, 'dist=',dist
C          print *, 'rzr=',rzr
C          print *, 'dr=', dr
          irdex=int(rzr/dr)
          rtan = rtani(irdex,iang)
          IF (RTAN.GT.0.0) GO TO 1999
6666      CONTINUE
C
          RZRN = RZR
          RZR = RZR*111.19493
C
c         PRINT*,'RO WILL BE RTAN TIMES SOME CONSTANT'
C
C
C
          RRDD = RNOT(iang)*111.19493
C
          IF(RRDD.GT.RZR)THEN
            PRINT*
            PRINT*,'RNOT HAS A NEGATIVE TANGENTIAL COMPONENT'
            PRINT*,'RNOT WAS DEFINED AS:  ',RRDD
            PRINT*,'RNOT WILL BE MODIFIED'
            RRDD = RZR
            RNOT(iang) = RZRN
          ENDIF
C
          PRINT*,'THIS IS RO IN KM:  ',RRDD
10      CONTINUE
C
C
C
C
C
C################################################################
C
C
C      DO THE OPTIMUM INTERPOLATION......>>>>>
C
        print *, 'entering rodist'
        call rodist
        print *, 'exiting rodist'
C
C
C     CREATE MATRIX  [A]  CONTAINING THE DISTANCE-RELATED CORRELATIONS
C     BETWEEN THE 24 BOUNDARY POINTS JUST OUTSIDE OF THE FILTER DOMAIN
C     RNOT......
C
C
        call amatrix
C
C
C     SEPERATE THE DISTURBANCE INTO THE HURRICNAE AND NON-HURRICNAE
C     COMPONENTS
C
C
        print *, TCTIME % isoformat()
        LONC = XC/PI180
        LATC = YC/PI180
         DO 880 J = 1, JMX
           DO 880 I = 1, IMX
             XXD(I,J) = UFIL(I,J) - UFILS(I,J)
880      CONTINUE
         print *, 'T=', T
         CALL SEPAR(XXD)
         DO 890 J = 1, JMX
           DO 890 I = 1, IMX
             UFIL(I,J)  =  UFILS(I,J) +  XXD(I,J)
890      CONTINUE
         DO 980 J = 1 , JMX
           DO 980 I = 1 , IMX
             XXD(I,J) = VFIL(I,J) - VFILS(I,J)
980      CONTINUE
         CALL SEPAR(XXD)
         DO 990 J = 1 , JMX
           DO 990 I = 1 , IMX
             VFIL(I,J)  =  VFILS(I,J) +  XXD(I,J)
990      CONTINUE

C
C
C      PUT THE ENVIRONMENTAL WINDS INTO THE GFDL HISTROY TAPE
C
        STATUS = NF90_PUT_VAR(UNCID, UVARID, UFIL,
     *                      START = (/ 1, 1, T /),
     *                      COUNT = (/ NLON, NLAT, 1 /))
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 611)
        STATUS = NF90_PUT_VAR(VNCID, VVARID, VFIL,
     *                      START = (/ 1, 1, T /),
     *                      COUNT = (/ NLON, NLAT, 1 /))
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 615)
c       END DO
       END DO
        STATUS = NF90_CLOSE(UNCID)
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 619)
        STATUS = NF90_CLOSE(VNCID)
        if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 621)
        STOP
      END

      SUBROUTINE PHASE(IFL,U,V,IMX,JMX,US,VS)
        implicit none
        integer, intent(in) :: ifl, imx, jmx
        integer :: nx
      PARAMETER  (NX=25)
CC************************************************************************
CC                                                                       *
CC     THIS SUBROUTINE CREATES  FILTERED  FIELDS OF (U,V) WIND           *
CC                                                                       *
CC                                                                       *
CC************************************************************************
CC
CC***********************************************************************
CC                   IMPORTANT!!!                                        *
CC   WE ASSUME THAT THE SPACING OF ALL THE POINTS IS ONE DEGREE          *
CC   LATITUDE AND LONGITUDE.                                             *
CC                                                                       *
CC                                                                       *
CC     IFL = THE STRENGTH OF THE FILTER VARYING FROM 1 (WEAK DAMPING) TO *
CC           4 (VERY STRONG DAMPING). WE ARE CURRENTLY USING IFL=2.      *
CC           THUS THERE ARE 4 CHOICES FOR THE TYPE OF FILTER DESIRED,    *
CC           IFL = 1, 2, 3, OR 4.                                        *
CC                                                                       *
CC                                                                       *
CC     U,V      =  INPUT OF THE UNSMOOTHED FIELDS                      *
CC                                                                       *
CC     IMX = NUMBER OF INPUT AND OUTPUT POINTS IN X-DIRECTION            *
CC     JMX = NUMBER OF INPUT AND OUTPUT POINTS IN Y-DIRECTION            *
CC                                                                       *
CC     US,VS   =  OUTPUT OF THE SMOOTHED FIELDS                       *
CC                                                                       *
CC************************************************************************
CC
CC
      real, DIMENSION(imx, jmx), intent(in) :: U, V
      real, dimension(imx, jmx), intent(out) :: US, VS
      real, dimension(nx) :: tk
      real, dimension(100) :: ampf
      real, DIMENSION(imx, nx) :: XTU,XTV
      real, DIMENSION(jmx, nx) :: YTU,YTV
      real :: amm, amp, amp1, chg, ckg, cosf, fact,
     *        pi, tn, tnn, tff, tfr, zz
      integer :: i,  imxm, irt, ismth, k, kt, kty,
     *           j, jmxm, na, nn, nty, ntym, nz
CC
CC
      IMXM  = IMX-1
      JMXM  = JMX-1
CC
CC

CC
      TN = FLOAT(NX)
CCC
      PI = 4.*ATAN(1.0)
      COSF = COS(2.*PI/TN) - 1.
CC
CC *************************************************************
CC
CC    ...IFL...  WILL CONTROL THE EXTENT OF DAMPING REQUESTED
CC
CC     IFL IS DETERMINED IN THE PROGRAM FILTER AT THE BEGINNING
CC
CC
CC    ...NTY...  IS THE NUMBER OF PASSES THROUGH THE SMOOTHING OPERATOR
CC
       IF(IFL.EQ.1)NTY = 8
       IF(IFL.EQ.2)NTY = 11
       IF(IFL.EQ.3)NTY = 17
       IF(IFL.EQ.4)NTY = 24
CC
CC**************************************************************
CC
CC  ISMTH: IS THE PARAMETER TO TURN ON DES-SMOOTHING. WE WILL ALWAYS ASSUME
CC         DESMOOTHING IS UNNECESSARY. HOWEVER IT IS STILL IN THE CODE FOR
CC         THE PURPOSE OF GENERALIZATION.
CC
       ISMTH = 0
CC
CC****************************************************************************
CC
CC
CC
CC  NEXT WE WILL DETERMINE THE SMOOTHING PARAMETER K TO BE USED
CC  DURING EACH OF N PASSES THROUGH THE SMOOTHING OPERATION.
CC
CC
CC
CC  (SEE THE APPENDIX OF KURIHARA ET AL., FROM THE MONTHLY WEATHER
CC   REVIEW ARTICLE, 1990 .....EQUATION A2).
CC
CC
       CHG = 0.0
       KT = 0
CC
CC
CC
         DO 802 KTY = 1 , NTY
CC
CC
CCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CC
CC
CC     FILTER 1...WEAK FILTER.....
CC
CC     N = 8 .... AND m VARIES AS 2,3,4,2,5,6,7,2
CC
CC
       IF((KTY.EQ.4.OR.KTY.EQ.8).AND.IFL.EQ.1)CHG = 1.0
CC
CC
CCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CC
CC     FILTER 2....REGULAR FILTER.....CURRENTLY IN USE
CC
CC     N = 11 .... AND m VARIES AS 2,3,4,2,5,6,7,2,8,9,2
CC
CC
       IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.GE.11)
     *  .AND.IFL.EQ.2)CHG = 1.0
CC
CC
CCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CC
CC     FILTER 3....STRONG FILTER.....EFFECTIVE FOR HURRICANE GILBERT
CC
CC     N = 17 .... AND m VARIES AS 2,3,4,2,5,6,7,2,8,9,10,2,11,2,2,2,2
CC
CC
       IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.EQ.12.OR.KTY.GE.14.)
     *  .AND.IFL.EQ.3)CHG = 1.0
CC
CC
CCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CC
CC
CC     FILTER 4.....VERY STRONG FILTER....THE PATTERN WILL START TO BECOME ZONAL....
CC
CC
CC     N = 24.......AND m VARIES AS :
CC                  2,3,4,2,5,6,7,2,8,9,10,2,11,12,13,2,2,2,2,2,2,2,2,2
CC
CC
CC
       IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.EQ.12
     *  .OR.KTY.GE.16).AND.IFL.EQ.4)CHG = 1.0
CC
CC
CC
CC***************************************************************************************8
CC
CC
       IF(CHG.EQ.0)KT = KT + 1
       IF(CHG.EQ.1.0)TK(KTY) = .25
       IF(CHG.EQ.1.0)GO TO 801
       FACT = 2.0*PI/(FLOAT(KT) + 1.0)
       TK(KTY) = -.5/(COS(FACT) - 1.0)
      DO 679 NA = 2 , 25
      AMPF(NA) = 1 + 2.*TK(KTY)*(COS(2.*PI/FLOAT(NA)) - 1.0)
679   CONTINUE
CCC
CCC
801    CONTINUE
       CHG = 0.0
802    CONTINUE
CC
CCC       WRITE(6,815) (TK(KK),KK = 1 , NTY)
CCC815    FORMAT(2X,'THIS IS TK:',E12.6)
CC
CC
CC ......THIS IS THE END OF THE KTY LOOP.......................
CC
CC
CC**********DESMOOTHING IS SET UP IF NEEDED********************
CC
CC
       IF(ISMTH.EQ.1)THEN
         NTYM = NTY - 1
         TFF = 1.0
         DO 61 K = 1 , NTYM
         TFF = TFF*(1. + 2.*TK(K)*COSF)
61       CONTINUE
         TFR = 1./TFF
         TK(NTY) = (TFR - 1.0)/(2.*COSF)
         WRITE(6,816) TK(NTY)
816      FORMAT(2X,'THE DESMOOTHING CONSTANT',E12.6)
       ENDIF
CC
C***********PRINT OUT THE DAMPING CHARACTERISTICS***************
CC
      IRT = KT+1
CC
      DO 610 NZ = 2 , 40
      AMP = 1.0
      TNN = FLOAT(NZ)
      CKG = 0.0
      IF(NZ.GT.IRT)CKG = 1.0
      DO 617 KT = 1 , NTY
      AMP1 = (1. + 2*TK(KT)*(COS(2.*PI/TNN)-1.0))
      IF(CKG.EQ.0.0)GO TO 619
      AMP =  AMP1*AMP
619   IF(ABS(AMP1).LT..01)CKG = 1.0
617   CONTINUE
      AMM = AMP
      IF(NZ.LE.IRT)AMM=0.0
CC
      ZZ = FLOAT(NZ)
      WRITE(11,455)ZZ,AMM
455   FORMAT(F8.3,F8.3)
CC
CC
CC
CC  THE FOLLOWING WRITE STATEMENT WILL LET YOU KNOW THE AMOUNT OF THE
CC  WAVE THAT HAS REMAINED AFTER THE FILTERING,
CC  FOR THE WAVE OF A GIVEN LENGTH D (WHICH IS CURRENTLY ONE DEGREE)
CC
      IF(NZ.EQ.20.OR.NZ.EQ.30.OR.NZ.EQ.40)THEN
      WRITE(6,677)NZ,AMM
677   FORMAT(2X,'WAVE NUMBER',I5,2X,'PERCENT WAVE REMAINING',E12.6)
      ENDIF
CC
610   CONTINUE
CC
CC*******************************************************************
CC
CC        DO THE SMOOTHING IN THE LATITUDINAL DIRECTION:
CC                         (EQUATION A1)
CC
        DO 600 J = 1 , JMX
CC
        DO 58 NN = 1 , NTY
            XTU(1,NN)   = U(1,J)
            XTU(IMX,NN) = U(IMX,J)
            XTV(1,NN)   = V(1,J)
            XTV(IMX,NN) = V(IMX,J)
58      CONTINUE
CCC
CC
CC
        DO 60 I = 2 , IMXM
      XTU(I,1) = U(I,J)   + TK(1)*(U(I-1,J) +
     *             U(I+1,J) - 2.*U(I,J))
      XTV(I,1) = V(I,J)   + TK(1)*(V(I-1,J) +
     *             V(I+1,J) - 2.*V(I,J))
60      CONTINUE
CC
CC
        DO 65 NN = 2 , NTY
        DO 62  I = 2 , IMXM
      XTU(I,NN) = XTU(I,NN-1) + TK(NN)*(XTU(I-1,NN-1) +
     *XTU(I+1,NN-1) - 2.*XTU(I,NN-1))
      XTV(I,NN) = XTV(I,NN-1) + TK(NN)*(XTV(I-1,NN-1) +
     *XTV(I+1,NN-1) - 2.*XTV(I,NN-1))
62    CONTINUE
65    CONTINUE
CC
       DO 70 I = 1 , IMX
      US(I,J)   = XTU(I,NTY)
      VS(I,J)   = XTV(I,NTY)
70    CONTINUE
CC
600   CONTINUE
C
C
CC********************************************************************
CC
CC    NOW DO THE SMOOTHING IN THE MERIDIONAL DIRECTION:
CC                         (EQUATION A3)
CC
CC
        DO 700   I = 1 , IMX
CC
        DO 80 NN = 1 , NTY
      YTU(1,NN)   = US(I,1)
      YTU(JMX,NN) = US(I,JMX)
      YTV(1,NN)   = VS(I,1)
      YTV(JMX,NN) = VS(I,JMX)
80      CONTINUE
CC
CC
CC
CC
        DO 90 J = 2 , JMXM
      YTU(J,1) = US(I,J) + TK(1)*(US(I,J-1) + US(I,J+1)
     *                          -2.*US(I,J))
      YTV(J,1) = VS(I,J) + TK(1)*(VS(I,J-1) + VS(I,J+1)
     *                          -2.*VS(I,J))
90      CONTINUE
CC
CC
        DO 95 NN = 2 , NTY
        DO 95 J  = 2 , JMXM
      YTU(J,NN) = YTU(J,NN-1) + TK(NN)*(YTU(J-1,NN-1)  +
     *              YTU(J+1,NN-1) - 2.*YTU(J,NN-1))
      YTV(J,NN) = YTV(J,NN-1) + TK(NN)*(YTV(J-1,NN-1)  +
     *              YTV(J+1,NN-1) - 2.*YTV(J,NN-1))
95    CONTINUE
CC
CC
CC   STORE THE FILTERED FIELDS IN US,VS AND GS
CC
CC
      DO 99 J = 1 , JMX
      US(I,J)   =  YTU(J,NTY)
      VS(I,J)   =  YTV(J,NTY)
99    CONTINUE
CC
CC
700   CONTINUE
CC
        RETURN
CC
CC
        END
      SUBROUTINE SEPAR(XD)
C
C  SEPERATES A FIELD INTO HURRICANE COMPONENT AND REMAINDER
C

       integer :: ix, iy
       real :: xv, yv, xold, yold, xcorn, ycorn
       PARAMETER( NMX=64,nmx1=nmx+1,nmx2=nmx*2,nmx6=nmx*6)
       PARAMETER (IMX=360 , JMX=180)
       DIMENSION XR(NMX),XD(IMX,JMX)
CC
       COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),
     *      DEL(IMX,JMX),THA(IMX,JMX),TANP(IMX,JMX),DS(IMX,JMX)
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
       COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
       COMMON /XXX/  XF(IMX,JMX),XC,YC,DX,DY
       COMMON /TOTAL/ DDEL,DTHA
C
c new arrays
        dimension b(nmx),w(nmx),ab(nmx,nmx1),ipvt(nmx)
     1       ,wrk(121),iwrk(25)
c     1       ,wrk(nmx6),iwrk(nmx2)
        common /matrix/ a(nmx,nmx),capd2
        common /vect/xvect(nmx),yvect(nmx)
c
        DATA XR/64*0./
C
C  XC,YC ARE HURRICANE COORDINATES
C  RO  IS RADIUS AT WHICH HURRICANE COMPONENT OF FIELD GOES TO ZERO
C  XR ARRAY CONTAINS THE FIELD VALUES OF 12 EQUALLY SPACED POINTS
C     ON CIRCLE OF RADIUS RO CENTERED AT XC,YC
C
c  set ro to be max value of rovect
c
c
        ro=0.
        do 22 i=1,nmx
        ro=amax1(ro,rovect(i))
22       continue
C        print*,'rovect=',rovect
C        print*,'ro=',ro,capd2,a(1,1),a(2,1)
        ro = ro*1.5
          PI = 4.*ATAN(1.0)
       PI180 = 4.*ATAN(1.0)/180.
       FACT =  COS(YOLD*PI180)
CC
CC   XC IS THE I POSITION OF THE CENTER OF THE OLD VORTEX
CC   YC IS THE J POSITION OF THE CENTER OF THE OLD VORTEX
CC   DDEL IS THE LONG. IN RADIANS OF THE OUTER NEST
CC   DTHA IS THE LAT.  IN RADIANS OF THE OUTER NEST
CC
c no fact here
c      DX=FACT*DDEL/PI180
c
       dx=ddel/pi180
       DY=DTHA/PI180
CC
cc
       XC = (XOLD-XCORN)*DX
       YC = (YOLD-YCORN)*DY
       IS=NINT((XC-RO/fact)/DX) +1.
       IE=NINT((XC+RO/fact)/DX + 1.)
       JS=NINT((YC-RO)/DY) +1.
       JE=NINT((YC+RO)/DY + 1.)
C
        DO 1 J = 1 , JMX
        DO 1 I = 1 , IMX
          XF(I,J)  = XD(I,J)
1       CONTINUE
C
C  SUBROUTINE BOUND COMPUTES FIELD VALUES OF ARRAY XR USING
C         BILINEAR INTERPOLATION
C
c
        Print*, 'calling BOUND from SEPAR '
C        print*,'here is xr ',xr
        CALL BOUND(NMX,XR,rovect)
C        print*,'here is rovect ',rovect
C
c  xrop(nmx) are the interpolated values of the disturbance
c   field at the rovect pts
c
c romax is the maximum value in rovect(nmx). Within the loop a local
c ro is computed for use in the separation. At the start of the loop
c ro is again set to romax to define the domain.
c
c
c
        w=0.
        romax=ro
C
C        print *, 'ro=',ro
C        print *, 'xold=', xold, 'yold=', yold
C        print *, 'xc=', xc, 'yc=', yc
C        print *, 'js=', JS, ' je=', JE
        DO 10 IX=IS,IE
        DO 11 JY=JS,JE
             ro=romax
c            X=XC-RO +DX*(IX-IS)
c            Y=YC-RO +DY*(JY-JS)
c            X= DX*float(IX)
c            Y= DY*float(JY)
             x=del(ix,jy)/pi180 -xcorn
             y=tha(ix,jy)/pi180 -ycorn
              delx=(x-xc)*fact
              dely=(y-yc)
             DR=SQRT((delx)**2 +(dely)**2)
             IF(DR.GT.RO)GOTO11
             IF(delx.ne.0.) THETA=ATAN((dely)/(delx))
             if(delx.eq.0..and.dely.lt.0.)theta=270.*pi180
             if(delx.eq.0..and.dely.gt.0.)theta=90. *pi180
             IF(delx.LT.0.)THETA=THETA+PI
             IF(THETA.LT.0.)THETA=2.*PI+THETA
             N1=INT(THETA*NMX/(2.*PI))
             IF(N1.GT.nmx)PRINT*,N1,THETA*57.296
             IF(N1.LT.0)PRINT*,N1,THETA*57.296
             N2=N1+2
             IF(N2.GT.NMX)N2=N2-NMX
             DELTH=THETA- 2.*PI*FLOAT(N1)/FLOAT(NMX)
c
             ro=delth*float(nmx)/(2.*pi)*(rovect(n2)-rovect(n1+1))
     1             +rovect(n1+1)
             IF(DR.GT.ro)GOTO11
             XRO=DELTH*FLOAT(NMX)/(2.*PI)*(XR(N2)-XR(N1+1)) +XR(N1+1)
CC
c Now add new code to compute distance from each gridpt. to rovect pts
c
             do 12 ip=1,nmx
             dpij= (fact*(x-xvect(ip)))**2 +(y-yvect(ip))**2
             b(ip)=exp(-dpij/capd2)
12           continue
c
c
             do 44 ip=1,nmx
               do 43 jp=1,nmx
43               ab(ip,jp)=a(ip,jp)
               ab(ip,nmx1)=b(ip)
44           continue
c
c  solve system using constrained least squares method
c
             call wnnls(ab,nmx,0,nmx,nmx,0,[1.],w,rnm,md,iwrk,wrk)
c
             temp=0.
             do 20 ip=1,nmx
             temp=temp +w(ip)*xr(ip)
20           continue
c            xh(ix,jy)=xf(ix,jy)-temp
             xd(ix,jy)=temp
11     CONTINUE
10     CONTINUE
       RETURN
       END
        SUBROUTINE BOUND(NMX,XR,ro)
C
        real :: xv, yv
        PARAMETER (IMX=360 , JMX=180)
C
        DIMENSION XR(NMX),ro(nmx)
c       COMMON  /IFACT/NNN,RO(nmx),RB,IENV
        COMMON  /XXX/  XF(IMX,JMX),XC,YC,DX,DY
        COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
C       COMMON  /COOR/ XV,YV,XOLD,YOLD
c       COMMON  /POSIT/ XOLD,YOLD
        PI = 4.*ATAN(1.0)
        fact=cos(yold*pi/180.)
        DO 10 I=1,NMX
        THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
        X=RO(i)/fact*COS(THETA)+XC +1.
C        print *, 'ro=', ro, ' x=', x
        Y=RO(i)*SIN(THETA)+YC +1.
        IX=NINT(X/DX)
        IY=NINT(Y/DY)
        IX1=IX+1
        IY1=IY+1
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
       XR(I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY+1)
     1      +  (1.-Q)*P*XF(IX+1,IY) + P*Q*XF(IX+1,IY+1)
10     CONTINUE
         RETURN
         END
       SUBROUTINE CENTER(UP,VP,DELG,THAG)
CC
       real :: xv, yv
       PARAMETER (IMX=360 , JMX=180, nmx=64)
CC     PARAMETER  (IMX=75, JMX=75)
       PARAMETER  ( KMAX=18,  LGI=20 )
       PARAMETER  (IGL = 500)
       COMMON  /GDINF/  NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,
     *                  JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
CC
      COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
       COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),
     *  DEL(IMX,JMX),THA(IMX,JMX),TANP(IMX,JMX),DS(IMX,JMX)
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
       COMMON /TOTAL/ DDEL,DTHA
       DIMENSION  UP(IMX,JMX),VP(IMX,JMX)
       DIMENSION  CMSUM(2,6),DLL(IGL),THH(IGL),WIND(IGL)
       DIMENSION  XM(IGL),RM(IGL)
       DIMENSION  TABPR(IGL),TFOUR(IGL),TFIVE(IGL),TSIX(IGL)
       DIMENSION  TPRES(IGL),DSS(IGL),tanw(imx,jmx)

C
      AFCT = 150.
CCCCC      AFCT = 1.0E10
      DFCT = 2.0 * AFCT
      RAD = 6.371E3
      PI = 4.*ATAN(1.0)
      PI180 = PI / 180.
      XCC = XV
      YCC = YV
      DX = PI180 * (XCC - XCORN) / DDEL
      DY = PI180 * (YCC - YCORN) / DTHA
      IX = NINT(DX) + 1
      IY = NINT(DY) + 1
      PRINT*
      PRINT*,'(x,y) OF Corn:  ',xcorn,ycorn
      PRINT*,'(x,y) OF CENTER:  ',xcc,ycc
      PRINT*,'(I,J) OF CENTER:  ',IX,IY
      PRINT*
C
C
C
       DDD = DEL(IX,IY) / PI180
C      DDD1 = DDD - 360.0
       TTT = THA(IX,IY) / PI180
       PRINT*,'(LON,LAT) OF CENTER: ',DDD,TTT
C
       DO 2 J = 1 , 6
       DO 2 I = 1 , 2
       CMSUM(I,J) = 0.0
2      CONTINUE
C
       ALPHA = .125
C
       IRANG = 5
C
       IB = IX - IRANG
       IE = IX + IRANG
       JB = IY - IRANG
       JE = IY + IRANG
       ITOT = (IE-IB+1)*(JE-JB+1)
C
C
       II = 0
       DO 10 J = JB, JE
       DO 10 I = IB, IE
C
       II = II + 1
C
       DSS(II) = DS(I,J)
       DLL(II) = DEL(I,J)
       THH(II) = THA(I,J)
       WIND(II) = SQRT(UP(I,J)*UP(I,J)+VP(I,J)*VP(I,J) )
10     CONTINUE
C
C
20     CONTINUE
       PRINT*
c      WRITE(6,3403)AMN100,AMX100
3403   FORMAT(5X,'MIN. AND MAX. WIND (M/S): ',F6.2,F9.2)
C
       DO 700 I = 1 , ITOT
       ANGL= .5*(THA(IX,IY) + THH(I) )
       COSF = COS(ANGL)
       DX = COSF*RAD*ABS(DLL(I) -  DEL(IX,IY) )
       DY = RAD*ABS(THH(I) -  THA(IX,IY) )
       RM(I) = SQRT(DX*DX+DY*DY)
700    CONTINUE
C
       DO 701 I = 1 , ITOT
       IF(RM(I).LT.AFCT)THEN
         XM(I) = 1.0
       ELSE
         XM(I) = EXP(-( (RM(I) - AFCT)/DFCT)**2)
       ENDIF
701    CONTINUE
C
C
       DO 60 I = 1 , ITOT
       CMSUM(1,1) =  CMSUM(1,1)+WIND(I)*DLL(I)*DSS(I)*XM(I)
       CMSUM(1,2) =  CMSUM(1,2)+WIND(I)*DSS(I)*XM(I)
       CMSUM(2,1) =  CMSUM(2,1)+WIND(I)*THH(I)*DSS(I)*XM(I)
       CMSUM(2,2) =  CMSUM(2,2)+WIND(I)*DSS(I)*XM(I)
60     CONTINUE
      DELG=  CMSUM(1,1)/CMSUM(1,2)
      THAG=  CMSUM(2,1)/CMSUM(2,2)
c
c  print the global position from set2 computation
c
       write(6,445) delg/pi180, thag/pi180
445    format(2x,'global position from windspeed',2f9.3)
C
C
432    format(11f7.1)
C
      PRINT*,'DISTANCE For MAX WIND Check (DEGREES):  ',DIST
C
      PRINT*
      RETURN
      END
      SUBROUTINE maxth(dumu,dumv,dxc,dyc,rmxlim,tw)
        real :: xv, yv
        parameter(nmx=64,imx=360,jmx=180,lgth=60,iimx=110)
        dimension dumu(imx,jmx),dumv(imx,jmx),tw(imx,jmx)
        dimension th(imx,jmx),tanmx(imx,jmx),tprof(7,7,lgth)
     1      ,itpos(7,7),tmax(7,7),tanavg(iimx)
       COMMON  /TOTAL/ DDEL,dtha
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
       COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),
     *      DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
c
        common /scale/rmxavg,rfind
c
        pi=4.*atan(1.0)
        PI180 = pi/180.
        fact=cos(yold)
        deltar=0.1
        dxc=xold/pi180-xcorn
        dyc=yold/pi180-ycorn
        ixc=nint(dxc)+1
        iyc=nint(dyc)+1
        ist=ixc-3
        jst=iyc-3
        iend=ixc+3
        jend=iyc+3
        npts=7
        print*,'ist,iend',ist,jst,iend,jend
c
c  compute radial profile of azimuthal avg. tang. wind at each pt
c
        do 51 i=ist,iend
        do 51 j=jst,jend
       xcen=(del(i,j)-del(1,1))/pi180 +1.
       ycen=(tha(i,j)-tha(1,1))/pi180 +1.
C         xcen=del(i,j)/pi180 +1.
C         ycen=tha(i,j)/pi180 +1.
         yyo=tha(i,j)
           do 52 ir=1,lgth
           rbd=float(ir)*0.2
c           print *, 'del=', del(i,j)
c           print *, 'tha=', tha(i,j)
c           print *, 'tha=', tha(1,1)
           call bound2(dumu,dumv,tanp,rbd,xcen,ycen,yyo)
           tprof(i-ist+1,j-jst+1,ir)=tanp
52         continue
51      continue
C        print *, 'let us see the value of tprof'
C        print 333,((tprof(i,1,k),i=4,7),(tprof(i,2,k),i=4,7),
C     *    (tprof(i,3,k),i=4,7),(tprof(i,4,k),i=4,7),k=1,lgth)
C333      format(16f7.1)
C        write(220,333) ((tprof(i,1,k),i=4,7),(tprof(i,2,k),i=4,7),
C     *    (tprof(i,3,k),i=4,7),(tprof(i,4,k),i=4,7),k=1,lgth)
c
c  find the first relative maximum along each azimuthal direction
c  find the position of the largest relative maximum
c
         hmax=0.
         do 53 i=1,npts
         do 53 j=1,npts
          do 54 ir=2,lgth-1
          if(tprof(i,j,ir).gt.tprof(i,j,ir-1).and.tprof(i,j,ir)
     1     .gt.tprof(i,j,ir+1))then
            tmax(i,j)=tprof(i,j,ir)
cc             itpos(i,j)=100*(ist+i)+j+jst
cc
cc       fixed the bug found april 22, 1994.......>>>
cc
        itpos(i,j)=jmx*(ist+i-1)+j+jst-1
        if (tmax(i,j).gt.hmax) THEN
          hmax = tmax(i,j)
          ipos = itpos(i,j)
          rmxpos = float(if)*0.2
        end if
C        hmax=amax1(tmax(i,j),hmax)
C        if(hmax.eq.tmax(i,j))ipos=itpos(i,j)
C        if(hmax.eq.tmax(i,j)) rmxpos=float(ir)*0.2
         goto53
          endif
54       continue
         tmax(i,j)=tprof(i,j,1)
         itpos(i,j)=101
         hmax=amax1(hmax,tmax(i,j))
        if(hmax.eq.tmax(i,j)) ipos=itpos(i,j)
53       continue
        print *, 'hmax, rmxpos, 5x rmxpos are'
        print*,hmax,rmxpos,rmxpos/0.2
        print *, 'tmax is'
        print *,((tmax(i,j),i=1,npts),j=1,npts)
        print *, 'itpos is'
        print *,((itpos(i,j),i=1,npts),j=1,npts)
        print *, 'ipos is ', ipos, mod(ipos,100)
c
c
c
c  use position of the largest relative maximum as the adjusted
c  center location
c
C       ycn=float(mod(ipos,100))-1.
C       xcn=float(ipos/100)-1.
       ycn=float(mod(ipos,jmx))-1.
       xcn=float(ipos/jmx)-1.
       ixc=int(xcn)+1
       iyc=int(ycn)+1
       print *, 'xcn=', xcn
       print *, 'ycn=', ycn
         xctest=(xcn+xcorn)*pi180
         yctest=(ycn+ycorn)*pi180
       print *, 'xctext=', xctest/pi180
       print *, 'yctext=', yctest/pi180
c
c  recompute the tangential wind component based on new center
c
       fact = cos(tha(1,iyc))
       print *,'in maxth',ycn,xcn,fact,xctest/pi180,yctest/pi180
      print*, 'ixc=',ixc
      print*, 'iyc=',iyc
      print*, shape(del)
      print*, 'del=', del(ixc+1,iyc+1)
       print*,ixc,iyc,del(ixc+1,iyc+1)/pi180,tha(ixc+1,iyc+1)/pi180
       do 334 j=1,jmx
       do 334 i=1,imx
        dx=(del(i,j)-xctest)*fact
        dy=(tha(i,j)-yctest)
       if(dx.ne.0.)theta =atan2(dy,dx)
       if(dy.gt.0..and.dx.eq.0.)theta =90.*pi180
       if(dy.lt.0..and.dx.eq.0.)theta =270.*pi180
       tw(i,j)=-dumu(i,j)*sin(theta) +dumv(i,j)*cos(theta)
         if(i.eq.ixc.and.j.eq.iyc)print*,i,j,dumu(i,j),dumv(i,j)
     1     ,theta,tw(i,j),dx,dy,'check everything'
334    continue
c
         write(77,*)ixc,iyc
        !write(77,7700)((tw(i,j)/100.,i=ixc-5,ixc+5),j=iyc+5,iyc-5,-1)
 7700   format(11f5.1)
c
        iflag=0
        hmax=0.
        rmxavg=0.
        do 50 ir=3,iimx
         rxx=float(ir)*deltar
         call calcr( rxx,rtan,xcn,ycn,yctest,dumu,dumv )
         tanavg(ir)=rtan
         if(tanavg(ir-2).lt.tanavg(ir-1).and.tanavg(ir).lt.tanavg(ir-1)
     1     .and.iflag.eq.0)then
           hmax=tanavg(ir-1)
           rmxavg=rxx-deltar
           iflag=1
         endif
50     continue
        print*,'found rmxavg ',rmxavg,hmax
        dxc=xcn
        dyc=ycn
C      xold=xcn+xcorn
700     format(10f6.1)
         print 700,tanavg
       call findra( dxc,dyc,yctest,rmxavg,rfavg,tanavg)
c
        alim = .75
        print*
        print*,'a factor to determine rmxlim: ',alim
c
        rmxlim = alim*rmxavg + (1.-alim)*rfavg
        print*,'found rfavg ',rfavg,rmxlim,dxc,dyc
        return
        end
       SUBROUTINE CALCRa(RO,RTAN,iang,dist)
       real :: xv, yv
       PARAMETER ( NMX=64)
       PARAMETER (IMX=360 , JMX=180)
       COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),
     *      DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
C
       COMMON  /TOTAL/ DDEL,DTHA
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
C
          PI = 4.*ATAN(1.0)
       PI180 = 4.*ATAN(1.0)/180.
       FACT =  COS(YOLD)
C
       DX=DDEL/PI180
       DY=DTHA/PI180
       XC = (XOLD-XCORN)*DX
       YC = (YOLD-YCORN)*DY

c
        THETA= 2.*PI*FLOAT(iang-1)/FLOAT(NMX)
        X=RO*COS(THETA)+XC
C        print *, cos(theta)
C        print *, 'x=', x, ' ro=', ro, ' xc=', xc, ' theta=', theta
        x = max(mod(x, 361.0),1.0)
        if (NINT(X).gt.360) X = X - 360
        if (NINT(X).le.0) X = X + 360
        Y=RO*SIN(THETA)+YC
        X1=X+DX
        x1 = max(mod(x1, 361.0),1.0)
        if (NINT(X1).gt.360) X1 = X1 - 360
        if (NINT(X1).le.0) X1 = X1 + 360
        Y1=Y+DY
        IX=nint(X/DX)
        IY=nint(Y/DY)
        IX1=NINT(X1/dx)
        IY1=nint(Y1/dy)
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
C        print *, 'ix=', ix, ' iy=', iy
       rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY1)
     1      +  (1.-Q)*P*XF(IX1,IY) + P*Q*XF(IX1,IY1)
10     CONTINUE
c
c
         RETURN
         END
       SUBROUTINE INTERP(UO,VO,SIGP,ILEV,UN,VN)
C
       PARAMETER ( KMAX = 18, ILT=20 )
C
       DIMENSION  UO(ILT),VO(ILT)
       DIMENSION  SIGP(ILT),Q(KMAX)
       DIMENSION  UN(KMAX),VN(KMAX)
       DIMENSION  UF(KMAX),VF(KMAX)
C
       DATA Q/ .9949968,.9814907,.9604809,.9204018,.8563145,.7772229,
     *           .6881255,.5935378,.4974484,.4248250,.3748014,.3247711,
     *           .2747291,.2246687,.1745733,.1244004,.0739862,.0207469/
C
       TMASK = 1.0e20
C
       ILEVM = ILEV -1
       IF(ILEV.LT.2)THEN
       UN = TMASK
       VN = TMASK
       RETURN
       ENDIF
       DO 10 KK = 1 , ILEVM
C
       SGO1 = SIGP(KK)
       SG02 = SIGP(KK+1)
C
       DO 20 K  = 1 ,  KMAX
C
       SGN = Q(K)
       IF(SGN.GT.SIGP(1))THEN
       UF(K) = TMASK
       VF(K) = TMASK
       ENDIF
       IF(SGN.LT.SIGP(ILEV))THEN
       UF(K) = TMASK
       VF(K) = TMASK
       ENDIF
C
       IF(SGN.LE.SGO1.AND.SGN.GE.SG02)THEN
       DX  = SGO1  - SG02
       DX1 = SGO1  - SGN
       DX2 = SGN  - SGO2
       X1 = DX1/DX
       X2 = DX2/DX
       UF(K) = (1.-X1)*UO(KK) + X1*UO(KK+1)
       VF(K) = (1.-X1)*VO(KK) + X1*VO(KK+1)
       ENDIF
C
20     CONTINUE
10     CONTINUE
       DO 30 K = 1 , KMAX
       KK = KMAX + 1 - K
       UN(K) = UF(KK)
       VN(K) = VF(KK)
30     CONTINUE
       RETURN
       END
        subroutine rodist
        real :: xv, yv
        parameter(nmx=64)
        common /vect/xvect(nmx),yvect(nmx)
       COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
c
        pi=4.0*atan(1.0)
        PI180 = 4.*ATAN(1.0)/180.
        yo=yold*pi180
        fact=cos(yo)
        xc=xold-xcorn
        yc=yold-ycorn
c
        do 10 ip=1,nmx
c
        theta=float(ip-1)/float(nmx)*2.*pi
        r=rovect(ip)
c
        xvect(ip)=r*cos(theta)/fact +xc
        yvect(ip)=r*sin(theta) +yc
10      continue
c
        return
        end
        subroutine amatrix
          real :: xv, yv
        parameter(nmx=64)
        common /matrix/ a(nmx,nmx),capd2
        common /vect/xvect(nmx),yvect(nmx)
       COMMON  /IFACT/NNN,rovect(nmx),RB,IENV
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
c
        PI180 = 4.*ATAN(1.0)/180.
        yo=yold*pi180
        fact=cos(yo)
        capd2=(2.25)*(2.25)
        do 10 ip=1,nmx
        do 10 jp=ip,nmx
          dpij=(fact*(xvect(ip)-xvect(jp)))**2 +(yvect(ip)-yvect(jp))**2
          a(ip,jp)= exp(-dpij/capd2)
          a(jp,ip)= a(ip,jp)
10      continue
100     format(5f8.4)
        return
        end
      SUBROUTINE WNNLS(W,MDW,ME,MA,N,L,PRGOPT,X,RNORM,MODE,IWORK,WORK)
C***BEGIN PROLOGUE  WNNLS
C***DATE WRITTEN   790701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  K1A2A
C***KEYWORDS  CONSTRAINED LEAST SQUARES,CURVE FITTING,DATA FITTING,
C             EQUALITY CONSTRAINTS,INEQUALITY CONSTRAINTS,
C             NONNEGATIVITY CONSTRAINTS,QUADRATIC PROGRAMMING
C***AUTHOR  HANSON, R. J., (SNLA)
C           HASKELL, K. H., (SNLA)
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality constraints and nonnegativity constraints on
C            selected variables.
C***DESCRIPTION
C
C     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),IWORK(M+N),WORK(M+5*N)
C
C     Written by Karen H. Haskell, Sandia Laboratories,
C     and R.J. Hanson, Sandia Laboratories.
C
C     Abstract
C
C     This subfoo solves a linearly constrained least squares
C     problem.  Suppose there are given matrices E and A of
C     respective dimensions ME by N and MA by N, and vectors F
C     and B of respective lengths ME and MA.  This subroutine
C     solves the problem
C
C               EX = F, (equations to be exactly satisfied)
C
C               AX = B, (equations to be approximately satisfied,
C                        in the least squares sense)
C
C               subject to components L+1,...,N nonnegative
C
C     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
C
C     The problem is reposed as problem WNNLS
C
C               (WT*E)X = (WT*F)
C               (   A)    (   B), (least squares)
C               subject to components L+1,...,N nonnegative.
C
C     The subfoo chooses the heavy weight (or penalty parameter) WT.
C
C     The parameters for WNNLS are
C
C     INPUT..
C
C     W(*,*),MDW,  The array W(*,*) is double subscripted with first
C     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
C                  discussion let us call M = ME + MA.  Then MDW
C                  must satisfy MDW.GE.M.  The condition MDW.LT.M
C                  is an error.
C
C                  The array W(*,*) contains the matrices and vectors
C
C                       (E  F)
C                       (A  B)
C
C                  in rows and columns 1,...,M and 1,...,N+1
C                  respectively.  Columns 1,...,L correspond to
C                  unconstrained variables X(1),...,X(L).  The
C                  remaining variables are constrained to be
C                  nonnegative. The condition L.LT.0 or L.GT.N is
C                  an error.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subfoo features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case LINK=1 and the values KEY and DATA SET
C                  are not referenced. The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=LINK1 (link to first entry of next group)
C               .  PRGOPT(2)=KEY1 (key to the option change)
C               .  PRGOPT(3)=DATA VALUE (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2)=DATA VALUE
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK.GT.NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the foo limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000 an error
C                  message is printed and the subfoo returns.
C
C                  OPTIONS..
C
C                  KEY=6
C                         Scale the nonzero columns of the
C                  entire data matrix
C                  (E)
C                  (A)
C                  to have length one. The DATA SET for
C                  this option is a single value.  It must
C                  be nonzero if unit length column scaling is
C                  desired.
C
C                  KEY=7
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  with a user-provided diagonal matrix.
C                  The DATA SET for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=8
C                         Change the rank determination tolerance from
C                  the nominal value of SQRT(SRELPR).  This quantity
C                  can be no smaller than SRELPR, The arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The DATA SET for this option
C                  is the new tolerance.
C
C                  KEY=9
C                         Change the blow-up parameter from the
C                  nominal value of SQRT(SRELPR).  The reciprocal of
C                  this parameter is used in rejecting solution
C                  components as too large when a variable is
C                  first brought into the active set.  Too large
C                  means that the proposed component times the
C                  reciprocal of the parameter is not less than
C                  the ratio of the norms of the right-side
C                  vector and the data matrix.
C                  This parameter can be no smaller than SRELPR,
C                  the arithmetic-storage precision.
C
C                  For example, suppose we want to provide
C                  a diagonal matrix to scale the problem
C                  matrix and change the tolerance used for
C                  determining linear dependence of dropped col
C                  vectors.  For these options the dimensions of
C                  PRGOPT(*) must be at least N+6.  The FORTRAN
C                  statements defining these options would
C                  be as follows.
C
C                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
C                  PRGOPT(2)=7 (user-provided scaling key)
C
C                  CALL SCOPY(N,D,1,PRGOPT(3),1) (copy the N
C                  scaling factors from a user array called D(*)
C                  into PRGOPT(3)-PRGOPT(N+2))
C
C                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
C                  PRGOPT(N+4)=8 (linear dependence tolerance key)
C                  PRGOPT(N+5)=... (new value of the tolerance)
C
C                  PRGOPT(N+6)=1 (no more options to change)
C
C
C     IWORK(1),    The amounts of working storage actually allocated
C     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
C                  respectively.  These quantities are compared with
C                  the actual amounts of storage needed for WNNLS( ).
C                  Insufficient storage allocated for either WORK(*)
C                  or IWORK(*) is considered an error.  This feature
C                  was included in WNNLS( ) because miscalculating
C                  the storage formulas for WORK(*) and IWORK(*)
C                  might very well lead to subtle and hard-to-find
C                  execution errors.
C
C                  The length of WORK(*) must be at least
C
C                  LW = ME+MA+5*N
C                  This test will not be made if IWORK(1).LE.0.
C
C                  The length of IWORK(*) must be at least
C
C                  LIW = ME+MA+N
C                  This test will not be made if IWORK(2).LE.0.
C
C     OUTPUT..
C
C     X(*)         An array dimensioned at least N, which will
C                  contain the N components of the solution vector
C                  on output.
C
C     RNORM        The residual norm of the solution.  The value of
C                  RNORM contains the residual vector length of the
C                  equality constraints and least squares equations.
C
C     MODE         The value of MODE indicates the success or failure
C                  of the subfoo.
C
C                  MODE = 0  Subfoo completed successfully.
C
C                       = 1  Max. number of iterations (equal to
C                            3*(N-L)) exceeded. Nearly all problems
C                            should complete in fewer than this
C                            number of iterations. An approximate
C                            solution and its corresponding residual
C                            vector length are in X(*) and RNORM.
C
C                       = 2  Usage error occurred.  The offending
C                            condition is noted with the error
C                            processing subfoo, XERROR( ).
C
C     User-designated
C     Working arrays..
C
C     WORK(*)      A real-valued working array of length at least
C                  M + 5*N.
C
C     IWORK(*)     An integer-valued working array of length at least
C                  M+N.
C***REFERENCES  K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, SAND77-0552, JUNE 1978.
C               K.H. HASKELL AND R.J. HANSON, *SELECTED ALGORITHMS FOR
C                 THE LINEARLY CONSTRAINED LEAST SQUARES PROBLEM--
C                 A USERS GUIDE*, SAND78-1290, AUGUST 1979.
C               K.H. HASKELL AND R.H. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, MATH. PROG. 21 (1981),
C                 PP. 98-118.
C               R.J. HANSON AND K.H. HASKELL, *TWO ALGORITHMS FOR THE
C                 LINEARLY CONSTRAINED LEAST SQUARES PROBLEM*, ACM
C                 TRANS. ON MATH. SOFTWARE, SEPT. 1982.
C***ROUTINES CALLED  WNLSM,XERROR,XERRWV
C***END PROLOGUE  WNNLS
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (START AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/, DUMMY/,SNGL(DUMMY)/
C
C     WRITTEN BY KAREN H. HASKELL, SANDIA LABORATORIES,
C     AND R.J. HANSON, SANDIA LABORATORIES.
C     REVISED FEB.25, 1982.
C
C     SUBROUTINES CALLED BY WNNLS( )
C
C++
C     WNLSM         COMPANION SUBROUTINE TO WNNLS( ), WHERE
C                   MOST OF THE COMPUTATION TAKES PLACE.
C
C     XERROR,XERRWV FROM SLATEC ERROR PROCESSING PACKAGE.
C                   THIS IS DOCUMENTED IN SANDIA TECH. REPT.,
C                   SAND78-1189.
C
C     REFERENCES
C
C     1. SOLVING LEAST SQUARES PROBLEMS, BY C.L. LAWSON
C        AND R.J. HANSON.  PRENTICE-HALL, INC. (1974).
C
C     2. BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE, BY
C        C.L. LAWSON, R.J. HANSON, D.R. KINCAID, AND F.T. KROGH.
C        TOMS, V. 5, NO. 3, P. 308.  ALSO AVAILABLE AS
C        SANDIA TECHNICAL REPORT NO. SAND77-0898.
C
C     3. AN ALGORITHM FOR LINEAR LEAST SQUARES WITH EQUALITY
C        AND NONNEGATIVITY CONSTRAINTS, BY K.H. HASKELL AND
C        R.J. HANSON.  AVAILABLE AS SANDIA TECHNICAL REPORT NO.
C        SAND77-0552, AND MATH. PROGRAMMING, VOL. 21, (1981), P. 98-118.
C
C     4. SLATEC COMMON MATH. LIBRARY ERROR HANDLING
C        PACKAGE.  BY R. E. JONES.  AVAILABLE AS SANDIA
C        TECHNICAL REPORT SAND78-1189.
C
      REAL              DUMMY, W(MDW,1), PRGOPT(1), X(1),  RNORM
      REAL WORK(ME+MA+5*N)
      INTEGER IWORK(N+ME+MA)
C
C
C***FIRST EXECUTABLE STATEMENT  WNNLS
      MODE = 0
       iwork(1)=mdw*6
       iwork(2)=mdw*2
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
      IF (.NOT.(IWORK(1).GT.0)) GO TO 20
      LW = ME + MA + 5*N
      IF (.NOT.(IWORK(1).LT.LW)) GO TO 10
      NERR = 2
      IOPT = 1
      print*,'work array',iwork(1),lw
      CALL XERRWV( 'WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR WORK(*)
     1, NEED LW=I1 BELOW', 70, NERR, IOPT, 1, LW, 0, 0, DUMMY, DUMMY)
      MODE = 2
      RETURN
   10 CONTINUE
   20 IF (.NOT.(IWORK(2).GT.0)) GO TO 40
      LIW = ME + MA + N
      IF (.NOT.(IWORK(2).LT.LIW)) GO TO 30
      NERR = 2
      IOPT = 1
      CALL XERRWV( 'WNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR IWORK(*
     1), NEED LIW=I1 BELOW', 72, NERR, IOPT, 1, LIW, 0, 0, DUMMY, DUMMY)
      MODE = 2
      RETURN
   30 CONTINUE
   40 IF (.NOT.(MDW.LT.ME+MA)) GO TO 50
      NERR = 1
      IOPT = 1
      CALL XERROR( 'WNNLS( ), THE VALUE MDW.LT.ME+MA IS AN ERROR', 44,
     1 NERR, IOPT)
      MODE = 2
      RETURN
   50 IF (0.LE.L .AND. L.LE.N) GO TO 60
      NERR = 2
      IOPT = 1
      CALL XERROR( 'WNNLS( ), L.LE.0.AND.L.LE.N IS REQUIRED', 39, NERR,
     1 IOPT)
      MODE = 2
      RETURN
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
C     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
C     REQUIRED BY THE MAIN SUBROUTINE WNLSM( ).
C
   60 L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
C
      CALL WNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,
     1 IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),
     2 WORK(L4),WORK(L5))
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subfoo nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE WNLSM(W,MDW,MME,MA,N,L,PRGOPT,X,RNORM,MODE,IPIVOT,
     1   ITYPE,WD,H,SCALE,Z,TEMP,D)
C***BEGIN PROLOGUE  WNLSM
C***REFER TO  WNNLS
C
C     This is a companion subfoo to WNNLS( ).
C     The documentation for WNNLS( ) has more complete
C     usage instructions.
C
C     Written by Karen H. Haskell, Sandia Laboratories,
C     with the help of R.J. Hanson, Sandia Laboratories,
C     December 1976 - January 1978.
C     Revised March 4, 1982.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine WNNLS, the following work arrays are used in
C     subroutine WNLSM  (they are passed through the calling
C     sequence from WNNLS for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C***ROUTINES CALLED  H12,ISAMAX,SASUM,SAXPY,SCOPY,SNRM2,SROTM,SROTMG,
C                    SSCAL,SSWAP,WNLIT,XERROR
C***END PROLOGUE  WNLSM
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (BEGIN CHANGES AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SROTMG/DROTMG/,
C     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/SROTM/DROTM/,/AMAX1/DMAX1/,
C     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/E0/D0/,/SSWAP/DSWAP/,
C     /ISAMAX/IDAMAX/,/SRELPR/DRELPR/
C
C     SUBROUTINE WNLSM (W,MDW,MME,MA,N,L,PRGOPT,X,RNORM,MODE,
C    1                  IPIVOT,ITYPE,WD,H,SCALE,Z,TEMP,D)
C++
      REAL             W(MDW,N), X(1), WD(1), H(1), SCALE(1), DOPE(4)
      REAL             Z(1), TEMP(1), PRGOPT(1), D(N), SPARAM(5)
      REAL             ALAMDA, ALPHA, ALSQ, AMAX, BNORM, EANORM
      REAL             SRELPR, FAC, ONE, BLOWUP
      REAL             RNORM, SM, T, TAU, TWO, WMAX, ZERO, ZZ, Z2
      REAL             AMAX1, SQRT, SNRM2, SASUM
      INTEGER IPIVOT(25), ITYPE(1), ISAMAX, IDOPE(8)
      LOGICAL HITCON, FEASBL, DONE, POS
      DATA ZERO /0.E0/, ONE /1.E0/, TWO /2.E0/, SRELPR /0.E0/
C
C     INITIALIZE-VARIABLES
C***FIRST EXECUTABLE STATEMENT  WNLSM
      IGO998 = 10
      GO TO 180
C
C     PERFORM INITIAL TRIANGULARIZATION IN THE SUBMATRIX
C     CORRESPONDING TO THE UNCONSTRAINED VARIABLES USING
C     THE PROCEDURE INITIALLY-TRIANGULARIZE.
   10 IGO995 = 20
      GO TO 280
C
C     PERFORM WNNLS ALGORITHM USING THE FOLLOWING STEPS.
C
C     UNTIL(DONE)
C
C        COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
C        WHEN (HITCON) ADD-CONSTRAINTS
C
C        ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
C
C        FIN
C
C     COMPUTE-FINAL-SOLUTION
C
   20 IF (DONE) GO TO 80
C
      IGO991 = 30
      GO TO 300
C
C     COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
   30 IF (.NOT.(HITCON)) GO TO 50
      IGO986 = 40
      GO TO 370
   40 GO TO 70
C
C     WHEN (HITCON) ADD-CONSTRAINTS
C
   50 IGO983 = 60
      GO TO 640
   60 CONTINUE
C
C     ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
C
   70 GO TO 20
C
   80 IGO980 = 90
      GO TO 1000
C
C     COMPUTE-FINAL-SOLUTION
C
   90 RETURN
  100 CONTINUE
C
C     TO PROCESS-OPTION-VECTOR
      FAC = 1.E-4
C
C     THE NOMINAL TOLERANCE USED IN THE CODE,
      TAU = SQRT(SRELPR)
C
C     THE NOMINAL BLOW-UP FACTOR USED IN THE CODE.
      BLOWUP = TAU
C
C     THE NOMINAL COLUMN SCALING USED IN THE CODE IS
C     THE IDENTITY SCALING.
      D(1) = ONE
      CALL SCOPY(N, D, 0, D, 1)
C
C     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE.
      NOPT = 1000
C
C     DEFINE BOUND FOR POSITIVE VALUE OF LINK.
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (.NOT.(LINK.LE.0 .OR. LINK.GT.NLINK)) GO TO 110
      NERR = 3
      IOPT = 1
      CALL XERROR( 'WNNLS( ) THE OPTION VECTOR IS UNDEFINED', 39, NERR,
     1 IOPT)
      MODE = 2
      RETURN
  110 IF (.NOT.(LINK.GT.1)) GO TO 160
      NTIMES = NTIMES + 1
      IF (.NOT.(NTIMES.GT.NOPT)) GO TO 120
      NERR = 3
      IOPT = 1
      CALL XERROR( 'WNNLS( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING
     1.', 53,     NERR, IOPT)
      MODE = 2
      RETURN
  120 KEY = PRGOPT(LAST+1)
      IF (.NOT.(KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.ZERO)) GO TO 140
      DO 130 J=1,N
        T = SNRM2(M,W(1,J),1)
        IF (T.NE.ZERO) T = ONE/T
        D(J) = T
  130 CONTINUE
  140 IF (KEY.EQ.7) CALL SCOPY(N, PRGOPT(LAST+2), 1, D, 1)
      IF (KEY.EQ.8) TAU = AMAX1(SRELPR,PRGOPT(LAST+2))
      IF (KEY.EQ.9) BLOWUP = AMAX1(SRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.NLINK)) GO TO 150
      NERR = 3
      IOPT = 1
      CALL XERROR( 'WNNLS( ) THE OPTION VECTOR IS UNDEFINED', 39, NERR,
     1 IOPT)
      MODE = 2
      RETURN
  150 LAST = LINK
      LINK = NEXT
      GO TO 110
  160 DO 170 J=1,N
        CALL SSCAL(M, D(J), W(1,J), 1)
  170 CONTINUE
      GO TO 1260
  180 CONTINUE
C
C     TO INITIALIZE-VARIABLES
C
C     SRELPR IS THE PRECISION FOR THE PARTICULAR MACHINE
C     BEING USED.  THIS LOGIC AVOIDS RECOMPUTING IT EVERY ENTRY.
      IF (.NOT.(SRELPR.EQ.ZERO)) GO TO 210
c*** changed back by BROSS
c*** changed by RF Boisvert, 19-Feb-92  (fails on HP 9000 Series 300)
cross      srelpr = r1mach(4)
       SRELPR = ONE
  190 IF (ONE+SRELPR.EQ.ONE) GO TO 200
       SRELPR = SRELPR/TWO
       GO TO 190
  200 SRELPR = SRELPR*TWO
cross
  210 M = MA + MME
      ME = MME
      MEP1 = ME + 1
      IGO977 = 220
      GO TO 100
C
C     PROCESS-OPTION-VECTOR
  220 DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      LP1 = L + 1
      NSOLN = L
      NSP1 = NSOLN + 1
      NP1 = N + 1
      NM1 = N - 1
      L1 = MIN0(M,L)
C
C     COMPUTE SCALE FACTOR TO APPLY TO EQUAL. CONSTRAINT EQUAS.
      DO 230 J=1,N
        WD(J) = SASUM(M,W(1,J),1)
  230 CONTINUE
      IMAX = ISAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = SASUM(M,W(1,NP1),1)
      ALAMDA = EANORM/(SRELPR*FAC)
C
C     DEFINE SCALING DIAG MATRIX FOR MOD GIVENS USAGE AND
C     CLASSIFY EQUATION TYPES.
      ALSQ = ALAMDA**2
      DO 260 I=1,M
C
C     WHEN EQU I IS HEAVILY WEIGHTED ITYPE(I)=0, ELSE ITYPE(I)=1.
        IF (.NOT.(I.LE.ME)) GO TO 240
        T = ALSQ
        ITEMP = 0
        GO TO 250
  240   T = ONE
        ITEMP = 1
  250   SCALE(I) = T
        ITYPE(I) = ITEMP
  260 CONTINUE
C
C     SET THE SOLN VECTOR X(*) TO ZERO AND THE COL INTERCHANGE
C     MATRIX TO THE IDENTITY.
      X(1) = ZERO
      CALL SCOPY(N, X, 0, X, 1)
      DO 270 I=1,N
        IPIVOT(I) = I
  270 CONTINUE
      GO TO 1230
  280 CONTINUE
C
C     TO INITIALLY-TRIANGULARIZE
C
C     SET FIRST L COMPS. OF DUAL VECTOR TO ZERO BECAUSE
C     THESE CORRESPOND TO THE UNCONSTRAINED VARIABLES.
      IF (.NOT.(L.GT.0)) GO TO 290
      WD(1) = ZERO
      CALL SCOPY(L, WD, 0, WD, 1)
C
C     THE ARRAYS IDOPE(*) AND DOPE(*) ARE USED TO PASS
C     INFORMATION TO WNLIT().  THIS WAS DONE TO AVOID
C     A LONG CALLING SEQUENCE OR THE USE OF COMMON.
  290 IDOPE(1) = ME
      IDOPE(2) = MEP1
      IDOPE(3) = 0
      IDOPE(4) = 1
      IDOPE(5) = NSOLN
      IDOPE(6) = 0
      IDOPE(7) = 1
      IDOPE(8) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = FAC
      DOPE(4) = TAU
      CALL WNLIT(W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     1 IDOPE, DOPE, DONE)
      ME = IDOPE(1)
      MEP1 = IDOPE(2)
      KRANK = IDOPE(3)
      KRP1 = IDOPE(4)
      NSOLN = IDOPE(5)
      NIV = IDOPE(6)
      NIV1 = IDOPE(7)
      L1 = IDOPE(8)
      GO TO 1240
  300 CONTINUE
C
C     TO COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
C     SOLVE THE TRIANGULAR SYSTEM OF CURRENTLY NON-ACTIVE
C     VARIABLES AND STORE THE SOLUTION IN Z(*).
C
C     SOLVE-SYSTEM
      IGO958 = 310
      GO TO 1110
C
C     INCREMENT ITERATION COUNTER AND CHECK AGAINST MAX. NUMBER
C     OF ITERATIONS.
  310 ITER = ITER + 1
      IF (.NOT.(ITER.GT.ITMAX)) GO TO 320
      MODE = 1
      DONE = .TRUE.
C
C     CHECK TO SEE IF ANY CONSTRAINTS HAVE BECOME ACTIVE.
C     IF SO, CALCULATE AN INTERPOLATION FACTOR SO THAT ALL
C     ACTIVE CONSTRAINTS ARE REMOVED FROM THE BASIS.
  320 ALPHA = TWO
      HITCON = .FALSE.
      IF (.NOT.(L.LT.NSOLN)) GO TO 360
      DO 350 J=LP1,NSOLN
        ZZ = Z(J)
        IF (.NOT.(ZZ.LE.ZERO)) GO TO 340
        T = X(J)/(X(J)-ZZ)
        IF (.NOT.(T.LT.ALPHA)) GO TO 330
        ALPHA = T
        JCON = J
  330   HITCON = .TRUE.
  340   CONTINUE
  350 CONTINUE
  360 GO TO 1220
  370 CONTINUE
C
C     TO ADD-CONSTRAINTS
C
C     USE COMPUTED ALPHA TO INTERPOLATE BETWEEN LAST
C     FEASIBLE SOLUTION X(*) AND CURRENT UNCONSTRAINED
C     (AND INFEASIBLE) SOLUTION Z(*).
      IF (.NOT.(LP1.LE.NSOLN)) GO TO 390
      DO 380 J=LP1,NSOLN
        X(J) = X(J) + ALPHA*(Z(J)-X(J))
  380 CONTINUE
  390 FEASBL = .FALSE.
      GO TO 410
  400 IF (FEASBL) GO TO 610
C
C     REMOVE COL JCON AND SHIFT COLS JCON+1 THROUGH N TO THE
C     LEFT. SWAP COL JCON INTO THE N-TH POSITION.  THIS ACHIEVES
C     UPPER HESSENBERG FORM FOR THE NONACTIVE CONSTRAINTS AND
C     LEAVES AN UPPER HESSENBERG MATRIX TO RETRIANGULARIZE.
  410 DO 420 I=1,M
        T = W(I,JCON)
        CALL SCOPY(N-JCON, W(I,JCON+1), MDW, W(I,JCON), MDW)
        W(I,N) = T
  420 CONTINUE
C
C     UPDATE PERMUTED INDEX VECTOR TO REFLECT THIS SHIFT AND SWAP.
      ITEMP = IPIVOT(JCON)
      IF (.NOT.(JCON.LT.N)) GO TO 440
      DO 430 I=JCON,NM1
        IPIVOT(I) = IPIVOT(I+1)
  430 CONTINUE
  440 IPIVOT(N) = ITEMP
C
C     SIMILARLY REPERMUTE X(*) VECTOR.
      CALL SCOPY(N-JCON, X(JCON+1), 1, X(JCON), 1)
      X(N) = ZERO
      NSP1 = NSOLN
      NSOLN = NSOLN - 1
      NIV1 = NIV
      NIV = NIV - 1
C
C     RETRIANGULARIZE UPPER HESSENBERG MATRIX AFTER ADDING CONSTRAINTS.
      J = JCON
      I = KRANK + JCON - L
  450 IF (.NOT.(J.LE.NSOLN)) GO TO 570
      IF (.NOT.(ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0)) GO TO 470
      IGO938 = 460
      GO TO 620
C
C     (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) ZERO-IP1-TO-I-IN-COL-J
  460 GO TO 560
  470 IF (.NOT.(ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1)) GO TO 490
      IGO938 = 480
      GO TO 620
C
C     (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) ZERO-IP1-TO-I-IN-COL-J
  480 GO TO 560
  490 IF (.NOT.(ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0)) GO TO 510
      CALL SSWAP(NP1, W(I,1), MDW, W(I+1,1), MDW)
      CALL SSWAP(1, SCALE(I), 1, SCALE(I+1), 1)
      ITEMP = ITYPE(I+1)
      ITYPE(I+1) = ITYPE(I)
      ITYPE(I) = ITEMP
C
C     SWAPPED ROW WAS FORMERLY A PIVOT ELT., SO IT WILL
C     BE LARGE ENOUGH TO PERFORM ELIM.
      IGO938 = 500
      GO TO 620
C
C     ZERO-IP1-TO-I-IN-COL-J
  500 GO TO 560
  510 IF (.NOT.(ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1)) GO TO 550
      T = SCALE(I)*W(I,J)**2/ALSQ
      IF (.NOT.(T.GT.TAU**2*EANORM**2)) GO TO 530
      IGO938 = 520
      GO TO 620
  520 GO TO 540
  530 CALL SSWAP(NP1, W(I,1), MDW, W(I+1,1), MDW)
      CALL SSWAP(1, SCALE(I), 1, SCALE(I+1), 1)
      ITEMP = ITYPE(I+1)
      ITYPE(I+1) = ITYPE(I)
      ITYPE(I) = ITEMP
      W(I+1,J) = ZERO
  540 CONTINUE
  550 CONTINUE
  560 I = I + 1
      J = J + 1
      GO TO 450
C
C     SEE IF THE REMAINING COEFFS IN THE SOLN SET ARE FEASIBLE.  THEY
C     SHOULD BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.  IF ANY ARE
C     INFEASIBLE IT IS DUE TO ROUNDOFF ERROR.  ANY THAT ARE NON-
C     POSITIVE WILL BE SET TO ZERO AND REMOVED FROM THE SOLN SET.
  570 IF (.NOT.(LP1.LE.NSOLN)) GO TO 590
      DO 580 JCON=LP1,NSOLN
        IF (X(JCON).LE.ZERO) GO TO 600
  580 CONTINUE
  590 FEASBL = .TRUE.
  600 CONTINUE
      GO TO 400
  610 GO TO 1200
  620 CONTINUE
C
C     TO ZERO-IP1-TO-I-IN-COL-J
      IF (.NOT.(W(I+1,J).NE.ZERO)) GO TO 630
      CALL SROTMG(SCALE(I), SCALE(I+1), W(I,J), W(I+1,J), SPARAM)
      W(I+1,J) = ZERO
      CALL SROTM(NP1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, SPARAM)
  630 GO TO 1290
  640 CONTINUE
C
C     TO PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
      CALL SCOPY(NSOLN, Z, 1, X, 1)
      IF (.NOT.(NSOLN.LT.N)) GO TO 650
      X(NSP1) = ZERO
      CALL SCOPY(N-NSOLN, X(NSP1), 0, X(NSP1), 1)
  650 I = NIV1
  660 IF (.NOT.(I.LE.ME)) GO TO 690
C
C     RECLASSIFY LEAST SQUARES EQATIONS AS EQUALITIES AS
C     NECESSARY.
      IF (.NOT.(ITYPE(I).EQ.0)) GO TO 670
      I = I + 1
      GO TO 680
  670 CALL SSWAP(NP1, W(I,1), MDW, W(ME,1), MDW)
      CALL SSWAP(1, SCALE(I), 1, SCALE(ME), 1)
      ITEMP = ITYPE(I)
      ITYPE(I) = ITYPE(ME)
      ITYPE(ME) = ITEMP
      MEP1 = ME
      ME = ME - 1
  680 GO TO 660
C
C     FORM INNER PRODUCT VECTOR WD(*) OF DUAL COEFFS.
  690 IF (.NOT.(NSP1.LE.N)) GO TO 730
      DO 720 J=NSP1,N
        SM = ZERO
        IF (.NOT.(NSOLN.LT.M)) GO TO 710
        DO 700 I=NSP1,M
          SM = SM + SCALE(I)*W(I,J)*W(I,NP1)
  700   CONTINUE
  710   WD(J) = SM
  720 CONTINUE
  730 GO TO 750
  740 IF (POS .OR. DONE) GO TO 970
C
C     FIND J SUCH THAT WD(J)=WMAX IS MAXIMUM.  THIS DETERMINES
C     THAT THE INCOMING COL J WILL REDUCE THE RESIDUAL VECTOR
C     AND BE POSITIVE.
  750 WMAX = ZERO
      IWMAX = NSP1
      IF (.NOT.(NSP1.LE.N)) GO TO 780
      DO 770 J=NSP1,N
        IF (.NOT.(WD(J).GT.WMAX)) GO TO 760
        WMAX = WD(J)
        IWMAX = J
  760   CONTINUE
  770 CONTINUE
  780 IF (.NOT.(WMAX.LE.ZERO)) GO TO 790
      DONE = .TRUE.
      GO TO 960
C
C     SET DUAL COEFF TO ZERO FOR INCOMING COL.
  790 WD(IWMAX) = ZERO
C
C     WMAX .GT. ZERO, SO OKAY TO MOVE COL IWMAX TO SOLN SET.
C     PERFORM TRANSFORMATION TO RETRIANGULARIZE, AND TEST
C     FOR NEAR LINEAR DEPENDENCE.
C     SWAP COL IWMAX INTO NSOLN-TH POSITION TO MAINTAIN UPPER
C     HESSENBERG FORM OF ADJACENT COLS, AND ADD NEW COL TO
C     TRIANGULAR DECOMPOSITION.
      NSOLN = NSP1
      NSP1 = NSOLN + 1
      NIV = NIV1
      NIV1 = NIV + 1
      IF (.NOT.(NSOLN.NE.IWMAX)) GO TO 800
      CALL SSWAP(M, W(1,NSOLN), 1, W(1,IWMAX), 1)
      WD(IWMAX) = WD(NSOLN)
      WD(NSOLN) = ZERO
      ITEMP = IPIVOT(NSOLN)
      IPIVOT(NSOLN) = IPIVOT(IWMAX)
      IPIVOT(IWMAX) = ITEMP
C
C     REDUCE COL NSOLN SO THAT THE MATRIX OF NONACTIVE
C     CONSTRAINTS VARIABLES IS TRIANGULAR.
  800 J = M
  810 IF (.NOT.(J.GT.NIV)) GO TO 870
      JM1 = J - 1
      JP = JM1
C
C     WHEN OPERATING NEAR THE ME LINE, TEST TO SEE IF THE PIVOT ELT.
C     IS NEAR ZERO.  IF SO, USE THE LARGEST ELT. ABOVE IT AS THE PIVOT.
C     THIS IS TO MAINTAIN THE SHARP INTERFACE BETWEEN WEIGHTED AND
C     NON-WEIGHTED ROWS IN ALL CASES.
      IF (.NOT.(J.EQ.MEP1)) GO TO 850
      IMAX = ME
      AMAX = SCALE(ME)*W(ME,NSOLN)**2
  820 IF (.NOT.(JP.GE.NIV)) GO TO 840
      T = SCALE(JP)*W(JP,NSOLN)**2
      IF (.NOT.(T.GT.AMAX)) GO TO 830
      IMAX = JP
      AMAX = T
  830 JP = JP - 1
      GO TO 820
  840 JP = IMAX
  850 IF (.NOT.(W(J,NSOLN).NE.ZERO)) GO TO 860
      CALL SROTMG(SCALE(JP), SCALE(J), W(JP,NSOLN), W(J,NSOLN), SPARAM)
      W(J,NSOLN) = ZERO
      CALL SROTM(NP1-NSOLN, W(JP,NSP1), MDW, W(J,NSP1), MDW, SPARAM)
  860 J = JM1
      GO TO 810
C
C     SOLVE FOR Z(NSOLN)=PROPOSED NEW VALUE FOR X(NSOLN).
C     TEST IF THIS IS NONPOSITIVE OR TOO LARGE.
C     IF THIS WAS TRUE OR IF THE PIVOT TERM WAS ZERO REJECT
C     THE COL AS DEPENDENT.
  870 IF (.NOT.(W(NIV,NSOLN).NE.ZERO)) GO TO 890
      ISOL = NIV
      IGO897 = 880
      GO TO 980
C
C     TEST-PROPOSED-NEW-COMPONENT
  880 GO TO 940
  890 IF (.NOT.(NIV.LE.ME .AND. W(MEP1,NSOLN).NE.ZERO)) GO TO 920
C
C     TRY TO ADD ROW MEP1 AS AN ADDITIONAL EQUALITY CONSTRAINT.
C     CHECK SIZE OF PROPOSED NEW SOLN COMPONENT.
C     REJECT IT IF IT IS TOO LARGE.
      ISOL = MEP1
      IGO897 = 900
      GO TO 980
C
C     TEST-PROPOSED-NEW-COMPONENT
  900 IF (.NOT.(POS)) GO TO 910
C
C     SWAP ROWS MEP1 AND NIV, AND SCALE FACTORS FOR THESE ROWS.
      CALL SSWAP(NP1, W(MEP1,1), MDW, W(NIV,1), MDW)
      CALL SSWAP(1, SCALE(MEP1), 1, SCALE(NIV), 1)
      ITEMP = ITYPE(MEP1)
      ITYPE(MEP1) = ITYPE(NIV)
      ITYPE(NIV) = ITEMP
      ME = MEP1
      MEP1 = ME + 1
  910 GO TO 930
  920 POS = .FALSE.
  930 CONTINUE
  940 IF (POS) GO TO 950
      NSP1 = NSOLN
      NSOLN = NSOLN - 1
      NIV1 = NIV
      NIV = NIV - 1
  950 CONTINUE
  960 GO TO 740
  970 GO TO 1250
  980 CONTINUE
C
C     TO TEST-PROPOSED-NEW-COMPONENT
      Z2 = W(ISOL,NP1)/W(ISOL,NSOLN)
      Z(NSOLN) = Z2
      POS = Z2.GT.ZERO
      IF (.NOT.(Z2*EANORM.GE.BNORM .AND. POS)) GO TO 990
      POS = .NOT.(BLOWUP*Z2*EANORM.GE.BNORM)
  990 GO TO 1280
 1000 CONTINUE
C     TO COMPUTE-FINAL-SOLUTION
C
C     SOLVE SYSTEM, STORE RESULTS IN X(*).
C
      IGO958 = 1010
      GO TO 1110
C     SOLVE-SYSTEM
 1010 CALL SCOPY(NSOLN, Z, 1, X, 1)
C
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO X(*) IF KRANK.LT.L
      IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.L)) GO TO 1030
      DO 1020 I=1,KRANK
        CALL H12(2, I, KRP1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
 1020 CONTINUE
C
C     FILL IN TRAILING ZEROES FOR CONSTRAINED VARIABLES NOT IN SOLN.
 1030 IF (.NOT.(NSOLN.LT.N)) GO TO 1040
      X(NSP1) = ZERO
      CALL SCOPY(N-NSOLN, X(NSP1), 0, X(NSP1), 1)
C
C     REPERMUTE SOLN VECTOR TO NATURAL ORDER.
 1040 DO 1070 I=1,N
        J = I
 1050   IF (IPIVOT(J).EQ.I) GO TO 1060
        J = J + 1
        GO TO 1050
 1060   IPIVOT(J) = IPIVOT(I)
        IPIVOT(I) = J
        CALL SSWAP(1, X(J), 1, X(I), 1)
 1070 CONTINUE
C
C     RESCALE THE SOLN USING THE COL SCALING.
      DO 1080 J=1,N
        X(J) = X(J)*D(J)
 1080 CONTINUE
      IF (.NOT.(NSOLN.LT.M)) GO TO 1100
      DO 1090 I=NSP1,M
        T = W(I,NP1)
        IF (I.LE.ME) T = T/ALAMDA
        T = (SCALE(I)*T)*T
        RNORM = RNORM + T
 1090 CONTINUE
 1100 RNORM = SQRT(RNORM)
      GO TO 1210
C
C     TO SOLVE-SYSTEM
C
 1110 CONTINUE
      IF (.NOT.(DONE)) GO TO 1120
      ISOL = 1
      GO TO 1130
 1120 ISOL = LP1
 1130 IF (.NOT.(NSOLN.GE.ISOL)) GO TO 1190
C
C     COPY RT. HAND SIDE INTO TEMP VECTOR TO USE OVERWRITING METHOD.
      CALL SCOPY(NIV, W(1,NP1), 1, TEMP, 1)
      DO 1180 JJ=ISOL,NSOLN
        J = NSOLN - JJ + ISOL
        IF (.NOT.(J.GT.KRANK)) GO TO 1140
        I = NIV - JJ + ISOL
        GO TO 1150
 1140   I = J
 1150   IF (.NOT.(J.GT.KRANK .AND. J.LE.L)) GO TO 1160
        Z(J) = ZERO
        GO TO 1170
 1160   Z(J) = TEMP(I)/W(I,J)
        CALL SAXPY(I-1, -Z(J), W(1,J), 1, TEMP, 1)
 1170   CONTINUE
 1180 CONTINUE
 1190 GO TO 1270
 1200 GO TO (40) IGO986
 1210 GO TO (90) IGO980
 1220 GO TO (30) IGO991
 1230 GO TO (10) IGO998
 1240 GO TO (20) IGO995
 1250 GO TO (60) IGO983
 1260 GO TO (220) IGO977
 1270 GO TO (310, 1010) IGO958
 1280 GO TO (880, 900) IGO897
 1290 GO TO (460, 480, 500, 520) IGO938
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts foo execution and prints error message.
C***DESCRIPTION
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the foo.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C***BEGIN PROLOGUE  H12
C***REFER TO  HFTI,LSEI,WNNLS
C
C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C
C     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
C     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
C
C     Modified at SANDIA LABS, May 1977, to --
C
C     1)  Remove double precision accumulation, and
C     2)  Include usage of the Basic Linear Algebra Package for
C         vectors longer than a particular threshold.
C
C     Construction and/or application of a single
C     Householder transformation..     Q = I + U*(U**T)/B
C
C     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
C     LPIVOT is the index of the pivot element.
C     L1,M   If L1 .LE. M   the transformation will be constructed to
C            zero elements indexed from L1 through M.   If L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
C                   IUE is the storage increment between elements.
C                                       On exit from H1 U() and UP
C                   contain quantities defining the vector U of the
C                   Householder transformation.   On entry to H2 U()
C                   and UP should contain quantities previously computed
C                   by H1.  These will not be modified by H2.
C     C()    On entry to H1 or H2 C() contains a matrix which will be
C            regarded as a set of vectors to which the Householder
C            transformation is to be applied.  On exit C() contains the
C            set of transformed vectors.
C     ICE    Storage increment between elements of vectors in C().
C     ICV    Storage increment between vectors in C().
C     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
C            no operations will be done on C().
C***ROUTINES CALLED  SAXPY,SDOT,SSWAP
C***END PROLOGUE  H12
      DIMENSION U(IUE,M), C(1)
C***FIRST EXECUTABLE STATEMENT  H12
      ONE=1.
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GO TO 60
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M
   10     CL=AMAX1(ABS(U(1,J)),CL)
      IF (CL) 130,130,20
   20 CLINV=ONE/CL
      SM=(U(1,LPIVOT)*CLINV)**2
          DO 30 J=L1,M
   30     SM=SM+(U(1,J)*CLINV)**2
      CL=CL*SQRT(SM)
      IF (U(1,LPIVOT)) 50,50,40
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GO TO 70
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN
      B=UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
      IF (B) 80,130,130
   80 B=ONE/B
      MML1P2=M-L1+2
      IF (MML1P2.GT.20) GO TO 140
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
          DO 120 J=1,NCV
          I2=I2+ICV
          I3=I2+INCR
          I4=I3
          SM=C(I2)*UP
              DO 90 I=L1,M
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE
          IF (SM) 100,120,100
  100     SM=SM*B
          C(I2)=C(I2)+SM*UP
              DO 110 I=L1,M
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE
  120     CONTINUE
  130 RETURN
  140 CONTINUE
      L1M1=L1-1
      KL1=1+(L1M1-1)*ICE
      KL2=KL1
      KLP=1+(LPIVOT-1)*ICE
      UL1M1=U(1,L1M1)
      U(1,L1M1)=UP
      IF (LPIVOT.EQ.L1M1) GO TO 150
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
  150 CONTINUE
          DO 160 J=1,NCV
          SM=SDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
          SM=SM*B
          CALL SAXPY (MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
          KL1=KL1+ICV
  160 CONTINUE
      U(1,L1M1)=UL1M1
      IF (LPIVOT.EQ.L1M1) RETURN
      KL1=KL2
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
      RETURN
      END
      SUBROUTINE WNLIT(W,MDW,M,N,L,IPIVOT,ITYPE,H,SCALE,RNORM,IDOPE,
     1   DOPE,DONE)
C***BEGIN PROLOGUE  WNLIT
C***REFER TO  WNNLS
C
C     This is a companion subfoo to WNNLS( ).
C     The documentation for WNNLS( ) has more complete
C     usage instructions.
C
C     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
C           B as the (N+1)st col.
C
C     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
C     col interchanges.
C     Revised March 4, 1982.
C***ROUTINES CALLED  H12,ISAMAX,SCOPY,SROTM,SROTMG,SSCAL,SSWAP
C***END PROLOGUE  WNLIT
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (BEGIN CHANGES AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/,/SROTM/DROTM/,
C     /SSCAL/DSCAL/,
C     /SSWAP/DSWAP/,/AMAX1/DMAX1/,/ISAMAX/IDAMAX/,/.E-/.D-/,/E0/D0/
C
C++
      REAL             W(MDW,1), H(1), SCALE(1), DOPE(4), SPARAM(5)
      REAL             ALSQ, AMAX, EANORM, FAC, FACTOR, HBAR, ONE, RN
      REAL             RNORM, SN, T, TAU, TENM3, ZERO
      REAL             AMAX1
      INTEGER ITYPE(1), IPIVOT(1), IDOPE(8)
      INTEGER ISAMAX
      LOGICAL INDEP, DONE, RECALC
      DATA TENM3 /1.E-3/, ZERO /0.E0/, ONE /1.E0/
C
C***FIRST EXECUTABLE STATEMENT  WNLIT
      ME = IDOPE(1)
      MEP1 = IDOPE(2)
      KRANK = IDOPE(3)
      KRP1 = IDOPE(4)
      NSOLN = IDOPE(5)
      NIV = IDOPE(6)
      NIV1 = IDOPE(7)
      L1 = IDOPE(8)
C
      ALSQ = DOPE(1)
      EANORM = DOPE(2)
      FAC = DOPE(3)
      TAU = DOPE(4)
      NP1 = N + 1
      LB = MIN0(M-1,L)
      RECALC = .TRUE.
      RNORM = ZERO
      KRANK = 0
C     WE SET FACTOR=1.E0 SO THAT THE HEAVY WEIGHT ALAMDA WILL BE
C     INCLUDED IN THE TEST FOR COL INDEPENDENCE.
      FACTOR = 1.E0
      I = 1
      IP1 = 2
      LEND = L
   10 IF (.NOT.(I.LE.LB)) GO TO 150
C
C     SET IR TO POINT TO THE I-TH ROW.
      IR = I
      MEND = M
      IGO996 = 20
      GO TO 460
C
C     UPDATE-COL-SS-AND-FIND-PIVOT-COL
   20 IGO993 = 30
      GO TO 560
C
C     PERFORM-COL-INTERCHANGE
C
C     SET IC TO POINT TO I-TH COL.
   30 IC = I
      IGO990 = 40
      GO TO 520
C
C     TEST-INDEP-OF-INCOMING-COL
   40 IF (.NOT.(INDEP)) GO TO 110
C
C     ELIMINATE I-TH COL BELOW DIAG. USING MOD. GIVENS TRANSFORMATIONS
C     APPLIED TO (A B).
      J = M
      DO 100 JJ=IP1,M
        JM1 = J - 1
        JP = JM1
C     WHEN OPERATING NEAR THE ME LINE, USE THE LARGEST ELT.
C     ABOVE IT AS THE PIVOT.
        IF (.NOT.(J.EQ.MEP1)) GO TO 80
        IMAX = ME
        AMAX = SCALE(ME)*W(ME,I)**2
   50   IF (.NOT.(JP.GE.I)) GO TO 70
        T = SCALE(JP)*W(JP,I)**2
        IF (.NOT.(T.GT.AMAX)) GO TO 60
        IMAX = JP
        AMAX = T
   60   JP = JP - 1
        GO TO 50
   70   JP = IMAX
   80   IF (.NOT.(W(J,I).NE.ZERO)) GO TO 90
        CALL SROTMG(SCALE(JP), SCALE(J), W(JP,I), W(J,I), SPARAM)
        W(J,I) = ZERO
        CALL SROTM(NP1-I, W(JP,IP1), MDW, W(J,IP1), MDW, SPARAM)
   90   J = JM1
  100 CONTINUE
      GO TO 140
  110 CONTINUE
      IF (.NOT.(LEND.GT.I)) GO TO 130
C
C     COL I IS DEPENDENT. SWAP WITH COL LEND.
      MAX = LEND
C
C     PERFORM-COL-INTERCHANGE
      IGO993 = 120
      GO TO 560
  120 CONTINUE
      LEND = LEND - 1
C
C     FIND COL IN REMAINING SET WITH LARGEST SS.
      MAX = ISAMAX(LEND-I+1,H(I),1) + I - 1
      HBAR = H(MAX)
      GO TO 30
  130 CONTINUE
      KRANK = I - 1
      GO TO 160
  140 I = IP1
      IP1 = IP1 + 1
      GO TO 10
  150 KRANK = L1
  160 CONTINUE
      KRP1 = KRANK + 1
      IF (.NOT.(KRANK.LT.ME)) GO TO 290
      FACTOR = ALSQ
      DO 170 I=KRP1,ME
        IF (L.GT.0) W(I,1) = ZERO
        CALL SCOPY(L, W(I,1), 0, W(I,1), MDW)
  170 CONTINUE
C
C     DETERMINE THE RANK OF THE REMAINING EQUALITY CONSTRAINT
C     EQUATIONS BY ELIMINATING WITHIN THE BLOCK OF CONSTRAINED
C     VARIABLES.  REMOVE ANY REDUNDANT CONSTRAINTS.
      LP1 = L + 1
      RECALC = .TRUE.
      LB = MIN0(L+ME-KRANK,N)
      I = LP1
      IP1 = I + 1
  180 IF (.NOT.(I.LE.LB)) GO TO 280
      IR = KRANK + I - L
      LEND = N
      MEND = ME
      IGO996 = 190
      GO TO 460
C
C     UPDATE-COL-SS-AND-FIND-PIVOT-COL
  190 IGO993 = 200
      GO TO 560
C
C     PERFORM-COL-INTERCHANGE
C
C     ELIMINATE ELEMENTS IN THE I-TH COL.
  200 J = ME
  210 IF (.NOT.(J.GT.IR)) GO TO 230
      JM1 = J - 1
      IF (.NOT.(W(J,I).NE.ZERO)) GO TO 220
      CALL SROTMG(SCALE(JM1), SCALE(J), W(JM1,I), W(J,I), SPARAM)
      W(J,I) = ZERO
      CALL SROTM(NP1-I, W(JM1,IP1), MDW, W(J,IP1), MDW, SPARAM)
  220 J = JM1
      GO TO 210
C
C     SET IC=I=COL BEING ELIMINATED
  230 IC = I
      IGO990 = 240
      GO TO 520
C
C     TEST-INDEP-OF-INCOMING-COL
  240 IF (INDEP) GO TO 270
C
C     REMOVE ANY REDUNDANT OR DEPENDENT EQUALITY CONSTRAINTS.
      JJ = IR
  250 IF (.NOT.(IR.LE.ME)) GO TO 260
      W(IR,1) = ZERO
      CALL SCOPY(N, W(IR,1), 0, W(IR,1), MDW)
      RNORM = RNORM + (SCALE(IR)*W(IR,NP1)/ALSQ)*W(IR,NP1)
      W(IR,NP1) = ZERO
      SCALE(IR) = ONE
C     RECLASSIFY THE ZEROED ROW AS A LEAST SQUARES EQUATION.
      ITYPE(IR) = 1
      IR = IR + 1
      GO TO 250
C
C     REDUCE ME TO REFLECT ANY DISCOVERED DEPENDENT EQUALITY
C     CONSTRAINTS.
  260 CONTINUE
      ME = JJ - 1
      MEP1 = ME + 1
      GO TO 300
  270 I = IP1
      IP1 = IP1 + 1
      GO TO 180
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE
      IF (.NOT.(KRANK.LT.L1)) GO TO 420
C
C     TRY TO DETERMINE THE VARIABLES KRANK+1 THROUGH L1 FROM THE
C     LEAST SQUARES EQUATIONS.  CONTINUE THE TRIANGULARIZATION WITH
C     PIVOT ELEMENT W(MEP1,I).
C
      RECALC = .TRUE.
C
C     SET FACTOR=ALSQ TO REMOVE EFFECT OF HEAVY WEIGHT FROM
C     TEST FOR COL INDEPENDENCE.
      FACTOR = ALSQ
      KK = KRP1
      I = KK
      IP1 = I + 1
  310 IF (.NOT.(I.LE.L1)) GO TO 410
C
C     SET IR TO POINT TO THE MEP1-ST ROW.
      IR = MEP1
      LEND = L
      MEND = M
      IGO996 = 320
      GO TO 460
C
C     UPDATE-COL-SS-AND-FIND-PIVOT-COL
  320 IGO993 = 330
      GO TO 560
C
C     PERFORM-COL-INTERCHANGE
C
C     ELIMINATE I-TH COL BELOW THE IR-TH ELEMENT.
  330 IRP1 = IR + 1
      J = M
      DO 350 JJ=IRP1,M
        JM1 = J - 1
        IF (.NOT.(W(J,I).NE.ZERO)) GO TO 340
        CALL SROTMG(SCALE(JM1), SCALE(J), W(JM1,I), W(J,I), SPARAM)
        W(J,I) = ZERO
        CALL SROTM(NP1-I, W(JM1,IP1), MDW, W(J,IP1), MDW, SPARAM)
  340   J = JM1
  350 CONTINUE
C
C     TEST IF NEW PIVOT ELEMENT IS NEAR ZERO. IF SO, THE COL IS
C     DEPENDENT.
      T = SCALE(IR)*W(IR,I)**2
      INDEP = T.GT.TAU**2*EANORM**2
      IF (.NOT.INDEP) GO TO 380
C
C     COL TEST PASSED. NOW MUST PASS ROW NORM TEST TO BE CLASSIFIED
C     AS INDEPENDENT.
      RN = ZERO
      DO 370 I1=IR,M
        DO 360 J1=IP1,N
          RN = AMAX1(RN,SCALE(I1)*W(I1,J1)**2)
  360   CONTINUE
  370 CONTINUE
      INDEP = T.GT.TAU**2*RN
C
C     IF INDEPENDENT, SWAP THE IR-TH AND KRP1-ST ROWS TO MAINTAIN THE
C     TRIANGULAR FORM.  UPDATE THE RANK INDICATOR KRANK AND THE
C     EQUALITY CONSTRAINT POINTER ME.
  380 IF (.NOT.(INDEP)) GO TO 390
      CALL SSWAP(NP1, W(KRP1,1), MDW, W(IR,1), MDW)
      CALL SSWAP(1, SCALE(KRP1), 1, SCALE(IR), 1)
C     RECLASSIFY THE LEAST SQ. EQUATION AS AN EQUALITY CONSTRAINT AND
C     RESCALE IT.
      ITYPE(IR) = 0
      T = SQRT(SCALE(KRP1))
      CALL SSCAL(NP1, T, W(KRP1,1), MDW)
      SCALE(KRP1) = ALSQ
      ME = MEP1
      MEP1 = ME + 1
      KRANK = KRP1
      KRP1 = KRANK + 1
      GO TO 400
  390 GO TO 430
  400 I = IP1
      IP1 = IP1 + 1
      GO TO 310
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
C
C     IF PSEUDORANK IS LESS THAN L, APPLY HOUSEHOLDER TRANS.
C     FROM RIGHT.
      IF (.NOT.(KRANK.LT.L)) GO TO 450
      DO 440 I=1,KRANK
        J = KRP1 - I
        CALL H12(1, J, KRP1, L, W(J,1), MDW, H(J), W, MDW, 1, J-1)
  440 CONTINUE
  450 NIV = KRANK + NSOLN - L
      NIV1 = NIV + 1
      IF (L.EQ.N) DONE = .TRUE.
C
C  END OF INITIAL TRIANGULARIZATION.
      IDOPE(1) = ME
      IDOPE(2) = MEP1
      IDOPE(3) = KRANK
      IDOPE(4) = KRP1
      IDOPE(5) = NSOLN
      IDOPE(6) = NIV
      IDOPE(7) = NIV1
      IDOPE(8) = L1
      RETURN
  460 CONTINUE
C
C     TO UPDATE-COL-SS-AND-FIND-PIVOT-COL
C
C     THE COL SS VECTOR WILL BE UPDATED AT EACH STEP. WHEN
C     NUMERICALLY NECESSARY, THESE VALUES WILL BE RECOMPUTED.
C
      IF (.NOT.(IR.NE.1 .AND. (.NOT.RECALC))) GO TO 480
C     UPDATE COL SS =SUM OF SQUARES.
      DO 470 J=I,LEND
        H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
  470 CONTINUE
C
C     TEST FOR NUMERICAL ACCURACY.
      MAX = ISAMAX(LEND-I+1,H(I),1) + I - 1
      RECALC = HBAR + TENM3*H(MAX).EQ.HBAR
C
C     IF REQUIRED, RECALCULATE COL SS, USING ROWS IR THROUGH MEND.
  480 IF (.NOT.(RECALC)) GO TO 510
      DO 500 J=I,LEND
        H(J) = ZERO
        DO 490 K=IR,MEND
          H(J) = H(J) + SCALE(K)*W(K,J)**2
  490   CONTINUE
  500 CONTINUE
C
C     FIND COL WITH LARGEST SS.
      MAX = ISAMAX(LEND-I+1,H(I),1) + I - 1
      HBAR = H(MAX)
  510 GO TO 600
  520 CONTINUE
C
C     TO TEST-INDEP-OF-INCOMING-COL
C
C     TEST THE COL IC TO DETERMINE IF IT IS LINEARLY INDEPENDENT
C     OF THE COLS ALREADY IN THE BASIS.  IN THE INIT TRI
C     STEP, WE USUALLY WANT THE HEAVY WEIGHT ALAMDA TO
C     BE INCLUDED IN THE TEST FOR INDEPENDENCE.  IN THIS CASE THE
C     VALUE OF FACTOR WILL HAVE BEEN SET TO 1.E0 BEFORE THIS
C     PROCEDURE IS INVOKED.  IN THE POTENTIALLY RANK DEFICIENT
C     PROBLEM, THE VALUE OF FACTOR WILL HAVE BEEN
C     SET TO ALSQ=ALAMDA**2 TO REMOVE THE EFFECT OF THE HEAVY WEIGHT
C     FROM THE TEST FOR INDEPENDENCE.
C
C     WRITE NEW COL AS PARTITIONED VECTOR
C             (A1)  NUMBER OF COMPONENTS IN SOLN SO FAR = NIV
C             (A2)  M-NIV COMPONENTS
C     AND COMPUTE  SN = INVERSE WEIGHTED LENGTH OF A1
C                  RN = INVERSE WEIGHTED LENGTH OF A2
C     CALL THE COL INDEPENDENT WHEN RN .GT. TAU*SN
      SN = ZERO
      RN = ZERO
      DO 550 J=1,MEND
        T = SCALE(J)
        IF (J.LE.ME) T = T/FACTOR
        T = T*W(J,IC)**2
        IF (.NOT.(J.LT.IR)) GO TO 530
        SN = SN + T
        GO TO 540
  530   RN = RN + T
  540   CONTINUE
  550 CONTINUE
      INDEP = RN.GT.TAU**2*SN
      GO TO 590
  560 CONTINUE
C
C     TO PERFORM-COL-INTERCHANGE
C
      IF (.NOT.(MAX.NE.I)) GO TO 570
C     EXCHANGE ELEMENTS OF PERMUTED INDEX VECTOR AND PERFORM COL
C     INTERCHANGES.
      ITEMP = IPIVOT(I)
      IPIVOT(I) = IPIVOT(MAX)
      IPIVOT(MAX) = ITEMP
      CALL SSWAP(M, W(1,MAX), 1, W(1,I), 1)
      T = H(MAX)
      H(MAX) = H(I)
      H(I) = T
  570 GO TO 580
  580 GO TO (30, 200, 330, 120) IGO993
  590 GO TO (40, 240) IGO990
  600 GO TO (20, 190, 320) IGO996
      END
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  910131   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     This is the CMLIB version of I1MACH, the integer machine
C     constants subroutine originally developed for the PORT library.
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C === MACHINE = CRAY.64-BIT-INTEGER
       DATA IMACH( 1) /     5 /
       DATA IMACH( 2) /     6 /
       DATA IMACH( 3) /   102 /
       DATA IMACH( 4) /     6 /
       DATA IMACH( 5) /    64 /
       DATA IMACH( 6) /     8 /
       DATA IMACH( 7) /     2 /
       DATA IMACH( 8) /    63 /
      !  DATA IMACH( 9) /  777777777777777777777B /
       DATA IMACH( 9) /  777777777777777777777 /
       DATA IMACH(10) /     2 /
       DATA IMACH(11) /    47 /
       DATA IMACH(12) / -8189 /
       DATA IMACH(13) /  8190 /
       DATA IMACH(14) /    94 /
       DATA IMACH(15) / -8099 /
       DATA IMACH(16) /  8190 /
c
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16)
     1   CALL XERROR ( 'I1MACH -- I OUT OF BOUNDS',25,1,2)
C
      I1MACH=IMACH(I)
      RETURN
C
      END
       SUBROUTINE calct(deltar,xc,yc,yo,xf,rmxlim)
c
c  calculates the radial profile for eight azimuthal angles
c
       real :: xv, yv
       PARAMETER ( nmx=64)
       PARAMETER (IMX=360 , JMX=180, iimx=100)
       DIMENSION xf(imx,jmx)
       dimension idst(nmx),hmax(nmx),rmax(nmx)
C
       common /scale/rmxavg,rfind
       common /pass/xr(iimx,nmx),dist(nmx)
       COMMON  /TOTAL/ DDEL,dtha
       COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
C
       PI = 4.*ATAN(1.0)
       PI180 = 4.*ATAN(1.0)/180.
       arad=6371.
       print*,'calct',deltar,xc,yc,yo,yold
       fact=cos(yo)
c
c      set the factor to determine lower limit of dist search
c
       bfact = .5
       rdistl =  bfact*rmxavg
c
       print*
       print*,'b factor to determine rdistl: ',bfact
c
c
c assume the maximum wind is within rfavg of center (but <10.8 also)
c
        irend =int(rmxlim/deltar)
c
c
        iravg =int(rdistl/deltar)
        print *, 'irend=', irend, ' iravg=', iravg
        if (irend.eq.0.or.iravg.eq.0) THEN
          return
        end if
c
        print*,'lower limit and radius of dist search: ',rdistl,iravg
        print*,'upper limit and radius of dist search: ',rmxlim,irend
C
       DX=DDEL*(1./PI180)
       DY=DTHA/PI180
c
c  angle loop
c
        DO 10 I=1,NMX
        THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
c
        do 11 ir=1,iimx
        ro=float(ir)*deltar
        X=(RO*COS(THETA)) +XC
        if (NINT(X).gt.360) X = X - 360
        Y=(RO*SIN(THETA)) +YC + 91.0
        X1=X+DX
        if (NINT(X1).gt.360) X1 = X1 - 360
        Y1=Y+DY
        IX=NINT(X/DX)
        IY=NINT(Y/DY)
        IX1=NINT(X1/DX)
        IY1=NINT(Y1/DY)
C        IX1=IX+1
C        IY1=IY+1
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
        if (IY.gt.jmx.or.iy1.gt.jmx) THEN
          CYCLE
        end if
       XR(ir,I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY1)
     1      +  (1.-Q)*P*XF(IX1,IY) + P*Q*XF(IX1,IY1)
11     continue
c
c find relative max after which ro check begins
c
          idst = 0
c
        hmax=-10.e10
        do 12 ir=iravg,irend
          hmax(i)=amax1(hmax(i),xr(ir,i))
          if(hmax(i).eq.xr(ir,i))dist(i)=float(ir)*deltar
          if(hmax(i).eq.xr(ir,i))idst(i)=ir
12      continue
c
c  if the max. value is also the endpt it maynot be a relative max
c  check backwards from irend for the last relative max
c
c
       irvgu = iravg+1
       if(irvgu.le.2)then
       irvgu = 3
       endif
c
       if(irend.eq.idst(i).or.iravg.eq.idst(i))then
         do 13 ir=irend,irvgu,-1
            if(xr(ir-1,i).lt.0.)goto14
            if(xr(ir-1,i).gt.xr(ir,i).and.xr(ir-1,i).ge.xr(ir-2,i))then
               dist(i)=float(ir-1)*deltar
               print*,'readjusting dist'
     1                ,dist(i),hmax(i),xr(ir-1,i),rmxlim
               go to 14
c
           endif
c
13      continue
c
14      continue
       endif
c
c
         if(idst(i).lt.iravg)then
         print*,'lower limit check, dist changed to rmxlim, i is: ',i
         dist(i) = rmxlim
         endif
c
c
10       CONTINUE
c
       do 450 iaa = 1 , nmx
       WRITE(94,*)iaa,dist(iaa)
450    continue
c
       write(94,*)rdistl,rmxlim
c
c
        dist=dist*1.1
       print*,'relative max found '
        print 4400,dist
4400    format(25f4.1)
          write(4,400)
     1        (float(ir)*deltar,(xr(ir,i)/100.,i=1,1),ir=1,iimx)
400     format(25f5.1)
         RETURN
         END
       SUBROUTINE calcr(RO,RTAN,xc,yc,yold,u,v)
       PARAMETER (nmx=64)
       PARAMETER (IMX=360 , JMX=180)
       DIMENSION XR(NMX),u(imx,jmx),v(imx,jmx)
C
       COMMON  /TOTAL/ DDEL,dtha
c      COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
C
          PI = 4.*ATAN(1.0)
       PI180 = pi/180.
c      FACT =  COS(YOLD*PI180)
c      yo=yold*pi180
       fact=cos(yold)
C
       DX=DDEL/PI180
       DY=DTHA/PI180
c      XC = (XOLD-XCORN)*DX
c      YC = (YOLD-YCORN)*DY
        DO 10 I=1,NMX
        THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
        X=RO*COS(THETA)/fact +XC +1.
        if (NINT(X).gt.360) X = X - 360
        Y=RO*SIN(THETA)+YC +1.
C        X=max(mod((RO*COS(THETA)/fact) +XC, 361.0),1.0)
C        Y=max(mod((RO*SIN(THETA)) +YC + 90.0, 181.0), 1.0) - 90.0
C        X1=max(mod(X+XC, 361.0),1.0)
C        Y1=max(mod(Y+YC, 181.0), 1.0) - 90.0
        X1=X+DX
        if (NINT(X1).gt.360) X1 = X1 - 360
        Y1=Y+DY
        IX=NINT(X/DX)
        IY=NINT(Y/DY)
        IX1=NINT(X1/DX)
        IY1=NINT(Y1/DY)
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
c      XR(I)=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY+1)
c    1      +  (1.-Q)*P*XF(IX+1,IY) + P*Q*XF(IX+1,IY+1)
       xr(i)=-sin(theta)*
     1    ((1.-P)*(1.-Q)*u(IX,IY) +(1.-P)*Q*u(IX,IY1)
     1      +  (1.-Q)*P*u(IX1,IY) + P*Q*u(IX1,IY1))
     1         +cos(theta)*
     1   ((1.-P)*(1.-Q)*v(IX,IY) +(1.-P)*Q*v(IX,IY1)
     1      +  (1.-Q)*P*v(IX1,IY) + P*Q*v(IX1,IY1))
10     CONTINUE
       RTAN = 0.0
c
c calculate azimuthally averaged tangential wind at radius ro
c
       DO 40 I = 1 , NMX
       RTAN = RTAN + XR(I)
40     CONTINUE
       RTAN = RTAN/FLOAT(NMX)
         RETURN
         END
        subroutine findra(dxc,dyc,yc,rmxavg,rfind,tanuv)
c
c  finds rfind from azimuthally averaged radial profile of tang. wind
c
        parameter(imx=360,jmx=180,nmx=64,iimx=110)
        dimension tanuv(iimx)
ccc       common /scale/rmxavg,rfind
c
        DR=1.0
c
c
        dist= rmxavg*1.5
        X1 = 0.0
        RTAN1 = 100000.
        R = 10.0
        r=dist
c
c only come back to 666 if rtan > 6m/s
c
666     CONTINUE
        rtan1=100000.
c
c  return to 777 if gradient, dist, or 3m/s are unmet
c
777     continue
c
        R = R + DR
        irad= int(r/dr)
        WRITE(64,*) 'irad=',irad
        print *, 'irad=', irad
        rtan= tanuv(irad)
c
        WRITE(56,*)R,RTAN
        RTAN2 = RTAN
        IF(RTAN.GT.600.)GO TO 666
c       PRINT*,'R AND RTAN:  ',R,RTAN
c
        IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5)GO TO 999
        IF(RTAN2.GE.RTAN1.AND.R.GT.DIST)THEN
        X1 = 1.0
c
        ENDIF
C
        IF(RTAN.LT.300..AND.R.GT.DIST)GO TO 999
        WRITE(56,*)R,RTAN
c
        RTAN1 = RTAN - 4.0
        IF(R.LT.10.8)GO TO 777
999     CONTINUE
C
        rfind=r
c
        return
        end
       subroutine bound2(u,v,tanuv,r0,xc,yc,yyo)
       PARAMETER(IMX=360,JMX=180,nmx=64)
       DIMENSION u(imx,jmx),v(imx,jmx),tani(nmx)
c       COMMON  /POSIT/ XOLD,YOLD,XCORN,YCORN,Rxx,XV,YV
       COMMON  /TOTAL/ DDEL,dtha
c       COMMON  /XXX/    XC,YC,DX,DY
        PI = 4.*ATAN(1.0)
        pi180=pi/180.
        dx=ddel/pi180
        dy=dtha/pi180
        fact=cos(yyo)
        arad=6371.
         r0=r0*111.19
        DO 10 I=1,NMX
        THETA= 2.*PI*FLOAT(I-1)/FLOAT(NMX)
        X=(R0*COS(THETA))/(arad*fact*pi180)+XC
C        print *, yc
C        print *, (R0*SIN(THETA))/(arad*fact *pi180)
        Y=(R0*SIN(THETA))/(arad*pi180)+YC
C        Y=MAX(mod((R0*SIN(THETA))/(arad*pi180)+YC+90.0, 181.0), 1.0) - 90.0
C        YN=MAX(mod((R0*SIN(THETA))/(arad*pi180)+YC+90.0, 181.0), 1.0)
        IX=NINT(X/DX)
        IY=NINT(Y/DY)
        IX1=IX+1
        IY1=IY+1
        P=X/DX-FLOAT(IX)
        Q=Y/DY-FLOAT(IY)
c      ui=(1.-P)*(1.-Q)*u(IX,IY) +(1.-P)*Q*u(IX,IY+1)
c    1      +  (1.-Q)*P*u(IX+1,IY) + P*Q*u(IX+1,IY+1)
c      vi=(1.-P)*(1.-Q)*v(IX,IY) +(1.-P)*Q*v(IX,IY+1)
c    1      +  (1.-Q)*P*v(IX+1,IY) + P*Q*v(IX+1,IY+1)
c      tani(i)=-sin(theta)*ui +cos(theta)*vi
       tani(i)=-sin(theta)*
     1    ((1.-P)*(1.-Q)*u(IX,IY) +(1.-P)*Q*u(IX,IY1)
     1      +  (1.-Q)*P*u(IX1,IY) + P*Q*u(IX1,IY1))
     1         +cos(theta)*
     1   ((1.-P)*(1.-Q)*v(IX,IY) +(1.-P)*Q*v(IX,IY1)
     1      +  (1.-Q)*P*v(IX1,IY) + P*Q*v(IX1,IY1))
10     CONTINUE
c
        tanuv=0.
        do 20 i=1,nmx
        tanuv=tanuv +tani(i)/float(nmx)
20      continue
         RETURN
         END
      SUBROUTINE HANDLE_ERR(STATUS, LINE)
        INTEGER, INTENT(IN) :: STATUS, LINE

        PRINT *, "ERROR: EXITING WITH STATUS ", STATUS, "AT LINE ", LINE
        STOP STATUS
      END SUBROUTINE HANDLE_ERR
