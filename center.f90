  SUBROUTINE CENTER(UP,VP,DELG,THAG)
!
    ! PARAMETER (IMX=640 , JMX=320, nmx=64)
!   PARAMETER  (IMX=75, JMX=75)
    PARAMETER  ( KMAX=18,  LGI=20 )
    PARAMETER  (IGL = 500)
    COMMON  /GDINF/  NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
!
    COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),DEL(IMX,JMX),THA(IMX,JMX),TANP(IMX,JMX),DS(IMX,JMX)
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
    COMMON /TOTAL/ DDEL,DTHA
    DIMENSION UP(IMX,JMX),VP(IMX,JMX)
    DIMENSION CMSUM(2,6),DLL(IGL),THH(IGL),WIND(IGL)
    DIMENSION XM(IGL),RM(IGL)
    DIMENSION DSS(IGL)
!
    AFCT = 150.
    ! AFCT = 1.0E10
    DFCT = 2.0 * AFCT
    RAD = 6.371E3
    XCC = XV
    YCC = YV
    DX = PI180 * (XCC - XCORN) / DDEL
    DY = PI180 * (YCC - YCORN) / DTHA
    IX = IFIX(DX)
    IY = IFIX(DY)
    PRINT*
    PRINT*,'(I,J) OF CENTER:  ',IX,IY
    PRINT*
!
!
!
    DDD = DEL(IX,IY) / PI180
!   DDD1 = DDD - 360.0
    TTT = THA(IX,IY) / PI180
    PRINT*,'(LON,LAT) OF CENTER: ',DDD,TTT
!
    DO J = 1,6
      DO I = 1,2
        CMSUM(I,J) = 0.0
      END DO
    END DO
!   
    ALPHA = .125
!
    IRANG = 5
!
    IB = IX - IRANG
    IE = IX + IRANG
    JB = IY - IRANG
    JE = IY + IRANG
    ITOT = (IE-IB+1)*(JE-JB+1)
!
!
    II = 0
    DO J = JB, JE
      DO I = IB, IE
        II = II + 1
!
        DSS(II) = DS(I,J)
        DLL(II) = DEL(I,J)
        THH(II) = THA(I,J)
        WIND(II) = SQRT(UP(I,J)*UP(I,J)+VP(I,J)*VP(I,J) )
      END DO
    END DO
!
!
! 20     CONTINUE
       PRINT*
!      WRITE(6,3403)AMN100,AMX100
! 3403   FORMAT(5X,'MIN. AND MAX. WIND (M/S): ',F6.2,F9.2)
!      
    DO 700 I = 1 , ITOT
      ANGL= .5*(THA(IX,IY) + THH(I) )
      COSF = COS(ANGL)
      DX = COSF*RAD*ABS(DLL(I) -  DEL(IX,IY) ) 
      DY = RAD*ABS(THH(I) -  THA(IX,IY) ) 
      RM(I) = SQRT(DX*DX+DY*DY) 
700 CONTINUE
!
    DO 701 I = 1 , ITOT
      IF(RM(I).LT.AFCT) THEN
        XM(I) = 1.0
      ELSE
        XM(I) = EXP(-( (RM(I) - AFCT)/DFCT)**2)
      ENDIF
701 CONTINUE
! 
!
    DO 60 I = 1,ITOT
      CMSUM(1,1) =  CMSUM(1,1)+WIND(I)*DLL(I)*DSS(I)*XM(I)             
      CMSUM(1,2) =  CMSUM(1,2)+WIND(I)*DSS(I)*XM(I)                    
      CMSUM(2,1) =  CMSUM(2,1)+WIND(I)*THH(I)*DSS(I)*XM(I)            
      CMSUM(2,2) =  CMSUM(2,2)+WIND(I)*DSS(I)*XM(I)                    
60  CONTINUE
    DELG=  CMSUM(1,1)/CMSUM(1,2)                 
    THAG=  CMSUM(2,1)/CMSUM(2,2)
!
!  print the global position from set2 computation
!
    write(6,445) delg/pi180, thag/pi180
445 format(2x,'global position from windspeed',2f9.3)
!
!
! 432 format(11f7.1)
!
    PRINT*,'DISTANCE For MAX WIND Check (DEGREES):  ',DIST
!
    PRINT*
    RETURN
  END SUBROUTINE CENTER