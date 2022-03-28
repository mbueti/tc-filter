SUBROUTINE CENTER(UP,VP,DELG,THAG)
  IMPLICIT NONE
  
   integer, PARAMETER :: IGL=500
   real, dimension(imx,jmx), intent(in) :: up, vp
   real, intent(inout) :: delg, thag
   integer :: dftx, dfty, ntr, nstflg, nn1, nn2, nn3, nn4,i, ib, ie, j, jb, je, iw, ix, iy, jjmax, ngd,ngr, jmax, itot, IRANG, ii, &
              ihx, ihy, icx, icy,ifl, im, imax, iimax, js, jn
   real :: xv, yv, afct, alpha, angl, cosf, ddd, pi, pi180,ddel, dtha, dfct, dist, xcc, ycc, xold, yold,xcorn, ycorn, ttt, rad, &
           dx, dy, dt, factr
   real, dimension(2,6) :: cmsum
   real, dimension(IGL) :: dll, thh, wind, xm, rm, tabpr, tfour,tfive, tsix, tpres, dss
   real, dimension(imx, jmx) :: tanw, del, tanp, ds, tang, tha
   real, dimension(imx,jmx,2) :: dmmm
   COMMON  /GDINF/  NGD,NGR,NTR,DT,JS,JN,IE,IW,IIMAX,IMAX,JJMAX,JMAX,NSTFLG,ICX,ICY,IHX,IHY,DFTX,DFTY
!
  COMMON /VAR/  DIST,NN1,NN2,NN3,NN4,IFL
  COMMON /WINDS/ DMMM,TANG,DEL,THA,TANP,DS
  COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,IX,IY
  COMMON /TOTAL/ DDEL,DTHA
   
  AFCT = 150.
  DFCT = 2.0 * AFCT
  RAD = 6.371E3
  PI = 4.*ATAN(1.0)
  PI180 = PI / 180.
  XCC = XV
  YCC = YV
  DX = PI180 * (XCC - XCORN) / DDEL
  DY = PI180 * (YCC - YCORN) / DTHA
  IX = NINT(DX) + 1
  IX = MODULO(IX-1, IMX)+1
  IY = NINT(DY) + 1
  PRINT*
  PRINT*,'(x,y) OF Corn:  ',xcorn,ycorn
  PRINT*,'(x,y) OF CENTER:  ',xcc,ycc
  PRINT*,'(I,J) OF CENTER:  ',IX,IY
  PRINT*
  
  DDD = DEL(IX,IY) / PI180
  TTT = THA(IX,IY) / PI180
  PRINT*,'(LON,LAT) OF CENTER: ',DDD,TTT
  
  DO J = 1 , 6
    DO I = 1 , 2
      CMSUM(I,J) = 0.0
    END DO
  END DO
  
  ALPHA = .125
  IRANG = 5
  
  IB = IX - IRANG
  IE = IX + IRANG
  JB = IY - IRANG
  JE = IY + IRANG
  ITOT = (IE-IB+1)*(JE-JB+1)
  
  II = 0
  DO J = JB, JE
    DO I = IB, IE
      II = II + 1
      IM = MODULO(I-1,IMX)+1
      DSS(II) = DS(IM,J)
      DLL(II) = DEL(IM,J)
      THH(II) = THA(IM,J)
      WIND(II) = SQRT(UP(IM,J)*UP(IM,J)+VP(IM,J)*VP(IM,J) )
    END DO
  END DO
  PRINT*
  
3403 FORMAT(5X,'MIN. AND MAX. WIND (M/S): ',F6.2,F9.2)

  DO I = 1 , ITOT
    ANGL= .5*(THA(IX,IY) + THH(I) )
    COSF = COS(ANGL)
    DX = COSF*RAD*ABS(DLL(I) -  DEL(IX,IY) )
    DY = RAD*ABS(THH(I) -  THA(IX,IY) )
    RM(I) = SQRT(DX*DX+DY*DY)
  END DO

  DO I = 1 , ITOT
    IF(RM(I).LT.AFCT)THEN
      XM(I) = 1.0
    ELSE
      XM(I) = EXP(-( (RM(I) - AFCT)/DFCT)**2)
    ENDIF
  END DO
!
!
  DO I = 1 , ITOT
   CMSUM(1,1) =  CMSUM(1,1)+WIND(I)*DLL(I)*DSS(I)*XM(I)
   CMSUM(1,2) =  CMSUM(1,2)+WIND(I)*DSS(I)*XM(I)
   CMSUM(2,1) =  CMSUM(2,1)+WIND(I)*THH(I)*DSS(I)*XM(I)
   CMSUM(2,2) =  CMSUM(2,2)+WIND(I)*DSS(I)*XM(I)
  END DO
  DELG = CMSUM(1,1)/CMSUM(1,2)
  THAG = CMSUM(2,1)/CMSUM(2,2)
!
!  print the global position from set2 computation
!
  write(6,445) delg/pi180, thag/pi180
445 format(2x,'global position from windspeed',2f9.3)
432 format(11f7.1)
  PRINT*,'DISTANCE For MAX WIND Check (DEGREES):  ',DIST
  
  PRINT*
  RETURN
END SUBROUTINE CENTER
