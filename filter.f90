  PROGRAM FILTER
    USE datetime_module, ONLY: datetime, timedelta
    USE netcdf
    use filter_routines

    INTEGER :: UNCID, VNCID, STATUS, UVARID, VVARID, &
               LATVARID, LONVARID, TIMEVARID, &
               NLON, NLAT, NTIME, TRACKREAD, &
               TCDATE, TCHOUR, STORMDAYS, STORMHOUR, STORMMIN, &
               TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR, &
               TRACKMIN, I, ITC, T, NTC, trackcount
    LOGICAL :: END
    REAL, DIMENSION(2) :: TCLAT, TCLON
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimIDs
    REAL, DIMENSION(:, :), ALLOCATABLE :: U, V
    REAL, DIMENSION(:), ALLOCATABLE :: LAT, LON, TIME
    CHARACTER(len=200) :: UFILENAME, VFILENAME, TRACKFILENAME
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

    PARAMETER (KMAX=18,LGI=20 ,iimx=110)
    PARAMETER (IBLKMX=LGI*IMX+4*KMAX*IMX)
    PARAMETER (JBLKMX=IMX+14*KMAX*IMX)
!
    COMMON /FILEC/  TYP1C(IMX,LGI),TYP2C(KMAX,IMX,4)
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX), &
                   DEL(IMX,JMX),THA(IMX,JMX),TANP(IMX,JMX),DS(IMX,JMX) 
!
!
    COMMON  /GDINF/ NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX, &
                    JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
    COMMON /pass/rtani(iimx,nmx),disti(nmx)
    COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
    COMMON /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
    COMMON /TOTAL/ DDEL,DTHA
    COMMON /IFACT/NNN,RNOT(nmx),RB,IENV
    REAL, DIMENSION(IMX) :: FACG1,FACG2,FACT1,FACT2
    REAL, DIMENSION(IBLKMX) :: FILC
    REAL, DIMENSION(JBLKMX) :: DIAG
    REAL, DIMENSION(IMX,JMX) :: XXD,YYD,US,VS,UALL,VALL,UP,VP, &
                                UFILS,VFILS,UFIL,VFIL,UFILP,VFILP
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

    trackcount = 0

    BASETIME = datetime(1900, 1, 1, 0, 0, 0)
    SIXHOURS = timedelta(0, 6, 0, 0, 0)
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

17    FORMAT(A4, 1x, A3, 1x, A10, I2, I2, I2, 1x, I2, I2, &
             1x, F3.0, A1, 1x, F4.0, A1, 1x, I3, 1x I3, 1x, I4, &
             1x, I4, 1x, I4, 1x, I2, 1x, I3, 1x, I4, 1x, I4, &
             1x, I4, 1x, I7, 1x, I4, 1x, I4, 1x, I4)
       
! 17     FORMAT(A4,A4,A10,1x,I4,I2,I2,1x,I2,I2,1x,F3.0,A1,1x,
!      *        F4.0,A1,1x,I3,1x,I3,3I5,
!      *        1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,
!      *        1x,I4,1x,I4,1x,I4)
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
11    TCIDS(NTCID + 1) = TCID
      NTCID = NTCID + 1
12    END = .FALSE.
13    CONTINUE
    END DO
    CLOSE(TRACKREAD)


!   DO ITC=1, NTCID
!   STORMID = TCIDS(ITC)
    STORMID = '08L'
    DO T=1, NTIME
    !  T = 2381
      STATUS = NF90_INQ_VARID(UNCID, UFIELDNAME, UVARID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 125)
      STATUS = NF90_GET_VAR(UNCID, UVARID, U, &
                            START = (/ 1, 1, T /), &
                            COUNT = (/ NLON, NLAT, 1 /))
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 129)
      STATUS = NF90_INQ_VARID(UNCID, "lat", LATVARID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 131)
      STATUS = NF90_GET_VAR(UNCID, LATVARID, LAT)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 133)
      STATUS = NF90_INQ_VARID(UNCID, "lon", LONVARID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 135)
      STATUS = NF90_GET_VAR(UNCID, LONVARID, LON)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 137)
      STATUS = NF90_INQ_VARID(UNCID, "time", TIMEVARID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 139)
      STATUS = NF90_GET_VAR(UNCID, TIMEVARID, TIME)

      STORMDAYS = FLOOR(TIME(T))
      STORMHOUR = FLOOR(12 * (TIME(T) - STORMDAYS))
      STORMMIN = FLOOR(60 * (12 * (TIME(T) - STORMDAYS) - STORMHOUR))

      DTIME = timedelta(STORMDAYS, STORMHOUR, STORMMIN, 0, 0)
      TCTIME = BASETIME + DTIME

      STATUS = NF90_open(VFILENAME, NF90_WRITE, VNCID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 153)
      STATUS = NF90_INQ_VARID(VNCID, VFIELDNAME, VVARID)
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 155)
      STATUS = NF90_GET_VAR(VNCID, VVARID, V, &
                            START = (/ 1, 1, T /), &
                            COUNT = (/ NLON, NLAT, 1 /))
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 159)

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
      DO WHILE ( &
!      *    (STORMID.NE.TCID.OR.
                ((TRACKTIME(2)).LT.TCTIME.OR. &
                 (TRACKTIME(1)).EQ.BASETIME).AND. &
                .NOT.END &
               ) 
        END = .TRUE.
      !     print * , 'STORMID=',STORMID
        IF (TCID.EQ.STORMID) THEN
          TCLAT(1) = TCLAT(2)
          TCLON(1) = TCLON(2)
          TCNS(1) = TCNS(2)
          TCEW(1) = TCEW(2)
          TRACKTIME(1) = TRACKTIME(2)
        END IF

        READ(TRACKREAD, 17, END=16) TCORG, TCID, TCNAME, TRACKYEAR, &
                                    TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN, &
                                    TCLAT(2), TCNS(2), TCLON(2), TCEW(2)
        if (TRACKYEAR < 50) then
          TRACKYEAR = TRACKYEAR + 2000
        else
          TRACKYEAR = TRACKYEAR + 1900
        end if
        TRACKTIME(2) = datetime(TRACKYEAR, TRACKMONTH, TRACKDAY, TRACKHOUR, TRACKMIN, 0)
        END=.FALSE.
16      CONTINUE
      END DO
      CLOSE(TRACKREAD)

      DTIMES(1) = TCTIME - TRACKTIME(1)
      DTIMES(2) = TRACKTIME(2) - TCTIME
      if ((DTIMES(1)%total_seconds()).lt.0 &
          .or.(DTIMES(2)%total_seconds()).lt.0) cycle

      trackcount = trackcount + 1
      if (trackcount < 5) cycle
   
!       IF (TCID.NE.STORMID.OR.
!    *      TRACKTIME(1).GT.TCTIME.OR.
!    *      TRACKTIME(2).LT.TCTIME
!    *  ) THEN
!          CYCLE
!       END IF

      !   PRINT *, TCID
      !   PRINT *, STORMID
      !   print *, TCLON
      !   print *, TCLAT
      !   print *, TCTIME % isoformat()
      !   print *, TRACKTIME(1) % isoformat()
      !   print *, TRACKTIME(2) % isoformat()

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

      print *, 'dt1=', DTIMES(1) % total_seconds()
      print *, 'dt2=', DTIMES(2) % total_seconds()

      print *,'dt1=',(1 - ((DTIMES(1) % total_seconds())/(SIXHOURS % total_seconds())))
      print *,'dt2=',(1 - (DTIMES(2) % total_seconds())/(SIXHOURS % total_seconds()))

      XV = TCLON(1) * (1 - (DTIMES(1) % total_seconds()) / (SIXHOURS % total_seconds())) + &
           TCLON(2) * (1 - (DTIMES(2) % total_seconds()) / (SIXHOURS % total_seconds())) 

      YV = TCLAT(1) * (1 - (DTIMES(1) % total_seconds()) / (SIXHOURS % total_seconds())) + &
           TCLAT(2) * (1 - (DTIMES(2) % total_seconds()) / (SIXHOURS % total_seconds()))

      PRINT *, 'XV=', XV
      PRINT *, 'YV=', YV
!   
      IF (XV.LT.0) THEN
        XV = 360. + XV
      ENDIF
!
!C****************SET UP THE FILTER STRENGTH YOU WANT***************
!
!   THIS IS THE FIRST FILTER, WHICH SEPERATES THE DISTURBANCE WIND
!   FIELD FROM THE BASIC FLOW.  THE BASIC FLOW WILL BE DEFINED AS
!   (US, VS). 
!    
!C   SEE THE SUBROUTINE PHASE FOR DETAILS.
!C
!C    IFL=1   IS THE WEAK    FILTER
!C    IFL=2   IS THE REGULAR FILTER *** CURRENTLY IN USE
!C    IFL=3   IS THE STRONG  FILTER
!C    IFL=4   IS VERY STRONG FILTER
!C
!C
!C    FILTER IS DEFINED IN MWR PAPER OF KURIHARA, ET.ALL, 1990:
!C
      IFL = 2
!C
!C
!C**********************************************************
!C
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
      DR = .1
!
      XC = XV*PI180
      YC = YV*PI180
! 
!
!
      J = 0
      DO 15 JJ = 1, JMX
        J = J + 1
!C
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
          DS(I,J) = (PI180 ** 2) * (LON(2) - LON(1))*(LAT(2) - LAT(1))
          DEL(I,J) = PI180 * LON(I)
          THA(I,J) = PI180 * LAT(J)
          UALL(I,J)   = U(I, J)
          VALL(I,J)   = V(I, J)
          UFIL(I,J)   = U(I, J)
          VFIL(I,J)   = V(I, J)
20      CONTINUE
15    CONTINUE
!
!      REWIND 10
!
      XCORN = DEL(1,1) / PI180
      YCORN = THA(1,1) / PI180
!C  
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
      DO J = 1, JMX
        DO I = 1, IMX
          UP(I,J) = UALL(I,J) - US(I,J)
          VP(I,J) = VALL(I,J) - VS(I,J)
          UFILP(I,J) = UFIL(I,J) - UFILS(I,J)
          VFILP(I,J) = VFIL(I,J) - VFILS(I,J)
        END DO
      END DO
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
      CALL MAXTH(up,vp,xcgnew,ycgnew,rmxlim,tanp)
      xold = xcgnew+xcorn
      yold = ycgnew+ycorn
      dist = rmxlim
!   
      xcp = xold*pi180
      ycp = yold*pi180
      write(4,*)xold,yold,xcp,ycp
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
      print*,'berore calct',dxc,dyc,ycp
!
!
!    CALCULATE THE RADIAL PROFILE OF TANGENTIAL WIND FOR 24 AXUMUTHAL 
!    ANGELS
!
      call calct(dr,dxc,dyc,ycp,tanp,rmxlim)
!
!          
      print *, 'nmx=', nmx
      DO 10 iang=1, nmx
        X1 = 0.0
        RTAN1 = 100000.
        R = 1.0
        dist=disti(iang)
        print *, 'disti=',disti
!
!  only return to 666 if rtan > 6m/s
!
666     CONTINUE
        Rtan1=100000.
!
!  return to 777 if dist or grad condition not met
!
777     continue
!   
!         CALL CALCRa(R,RTAN,iang)
        irdex=int(r/dr)
        rtan = rtani(irdex,iang)
        R = R + DR 
!       WRITE(56,*)R,RTAN
        RTAN2 = RTAN
        IF(RTAN.GT.600.)GO TO 666
!       PRINT*,'R AND RTAN:  ',R,RTAN 
        IF(RTAN2.GE.RTAN1.AND.R.GT.DIST.AND.X1.GT..5)GO TO 999
        IF(RTAN2.GE.RTAN1.AND.R.GT.DIST) THEN
          X1 = 1.0
        ENDIF
!   
        IF(RTAN.LT.300..AND.R.GT.DIST)GO TO 999
!       WRITE(56,*)R,RTAN
        RTAN1 = RTAN -4.0
        IF(R.LT.10.8)GO TO 777
999     CONTINUE 
        print *, 999
!       PRINT*
!   
!   
        IF (X1.EQ.1.0) THEN
          RNOT(iang) = (R-.1)/.8
        ELSE
          RNOT(iang) =  R/.8
        ENDIF
!   
!       RZR = R  
        rzr=dist
1999    CONTINUE
        RZR = RZR + DR
        CALL CALCRa(RZR,RTAN,iang)
        irdex=min(floor(rzr/dr), 110)
        rtan = rtani(irdex,iang)
        print *, 'rtan=', rtan
        IF (RTAN.GT.0.0) GO TO 1999
!   
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
10    CONTINUE
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
      call rodist
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
!          DO 880 J = 1, JMX
!            DO 880 I = 1, IMX
!              XXD(I,J) = UFIL(I,J) - UFILS(I,J)
! 880      CONTINUE
      XXD = UFIL - UFILS
      CALL SEPAR(XXD)
      UFIL = UFILS + XXD
!          DO 890 J = 1, JMX
!            DO 890 I = 1, IMX
!              UFIL(I,J)  =  UFILS(I,J) +  XXD(I,J)
! 890      CONTINUE
      YYD = VFIL - VFILS
!          DO 980 J = 1 , JMX
!            DO 980 I = 1 , IMX
!              XXD(I,J) = VFIL(I,J) - VFILS(I,J)
! 980      CONTINUE
      CALL SEPAR(YYD)
      VFIL = VFILS + YYD
!          DO 990 J = 1 , JMX
!            DO 990 I = 1 , IMX
!              VFIL(I,J)  =  VFILS(I,J) +  XXD(I,J)
! 990      CONTINUE
      !    print *, 'XXD=', XXD
!
!
!      PUT THE ENVIRONMENTAL WINDS INTO THE GFDL HISTROY TAPE 
!   
      print *, 'T=', T
      STATUS = NF90_PUT_VAR(UNCID, UVARID, UFIL, &
                            START = (/ 1, 1, T /), &
                            COUNT = (/ NLON, NLAT, 1 /)) 
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 618)
      STATUS = NF90_PUT_VAR(VNCID, VVARID, VFIL, &
                            START = (/ 1, 1, T /), &
                            COUNT = (/ NLON, NLAT, 1 /))
      if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 622)
    END DO
!       END DO
    STATUS = NF90_CLOSE(UNCID)
    if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 626)
    STATUS = NF90_CLOSE(VNCID)
    if(STATUS /= NF90_NOERR) call HANDLE_ERR(STATUS, 628)
      !   STOP 'stopped'
  END PROGRAM FILTER
