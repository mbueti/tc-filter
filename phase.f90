  SUBROUTINE PHASE(IFL,U,V,IMX,JMX,US,VS)
    PARAMETER  (NX=25)
!************************************************************************
!                                                                       *
!     THIS SUBROUTINE CREATES  FILTERED  FIELDS OF (U,V) WIND           *
!                                                                       *
!                                                                       *
!************************************************************************
!                                                                        
!***********************************************************************
!                   IMPORTANT!!!                                        *
!   WE ASSUME THAT THE SPACING OF ALL THE POINTS IS ONE DEGREE          *
!   LATITUDE AND LONGITUDE.                                             *
!                                                                       * 
!                                                                       *
!     IFL = THE STRENGTH OF THE FILTER VARYING FROM 1 (WEAK DAMPING) TO *
!           4 (VERY STRONG DAMPING). WE ARE CURRENTLY USING IFL=2.      *
!           THUS THERE ARE 4 CHOICES FOR THE TYPE OF FILTER DESIRED,    * 
!           IFL = 1, 2, 3, OR 4.                                        *
!                                                                       *
!                                                                       *
!     U,V      =  INPUT OF THE UNSMOOTHED FIELDS                      *
!                                                                       *
!     IMX = NUMBER OF INPUT AND OUTPUT POINTS IN X-DIRECTION            *
!     JMX = NUMBER OF INPUT AND OUTPUT POINTS IN Y-DIRECTION            *
!                                                                       *
!     US,VS   =  OUTPUT OF THE SMOOTHED FIELDS                       *
!                                                                       *
!************************************************************************
!                                                                        
!   
    DIMENSION    U(IMX,JMX), V(IMX,JMX) 
    DIMENSION    US(IMX,JMX),VS(IMX,JMX)
    DIMENSION    TK(NX),AMPF(100)  
    DIMENSION    XTU(IMX,NX),XTV(IMX,NX)           
    DIMENSION    YTU(JMX,NX),YTV(JMX,NX)           
!                                                                        
! 
    IMXM  = IMX-1                                                                       
    JMXM  = JMX-1                                                                       
!                                              
!                                                                      
                                                  
!                                                                      
    TN = FLOAT(NX)                                                       
!                                                                     
    ! PI = 4.*ATAN(1.0)                                                   
    COSF = COS(2.*PI/TN) - 1. 
!                                          
! *************************************************************
!
!    ...IFL...  WILL CONTROL THE EXTENT OF DAMPING REQUESTED
!
!     IFL IS DETERMINED IN THE PROGRAM FILTER AT THE BEGINNING
!
!
!    ...NTY...  IS THE NUMBER OF PASSES THROUGH THE SMOOTHING OPERATOR
!                                                                       
    IF(IFL.EQ.1)NTY = 8
    IF(IFL.EQ.2)NTY = 11
    IF(IFL.EQ.3)NTY = 17
    IF(IFL.EQ.4)NTY = 24
!
!**************************************************************
!
!  ISMTH: IS THE PARAMETER TO TURN ON DES-SMOOTHING. WE WILL ALWAYS ASSUME
!         DESMOOTHING IS UNNECESSARY. HOWEVER IT IS STILL IN THE CODE FOR
!         THE PURPOSE OF GENERALIZATION.
!
    ISMTH = 0
!
!****************************************************************************
! 
!
!                                                                
!  NEXT WE WILL DETERMINE THE SMOOTHING PARAMETER K TO BE USED
!  DURING EACH OF N PASSES THROUGH THE SMOOTHING OPERATION.
!
!
!
!  (SEE THE APPENDIX OF KURIHARA ET AL., FROM THE MONTHLY WEATHER
!   REVIEW ARTICLE, 1990 .....EQUATION A2).
!
!                                                          
    CHG = 0.0                                                          
    KT = 0
!
!
!
    DO 802 KTY = 1,NTY  
! 
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!
!     FILTER 1...WEAK FILTER.....
!     
!     N = 8 .... AND m VARIES AS 2,3,4,2,5,6,7,2
!
!                                            
      IF((KTY.EQ.4.OR.KTY.EQ.8).AND.IFL.EQ.1)CHG = 1.0  
! 
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     FILTER 2....REGULAR FILTER.....CURRENTLY IN USE
!
!     N = 11 .... AND m VARIES AS 2,3,4,2,5,6,7,2,8,9,2
!
!                  
      IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.GE.11).AND.IFL.EQ.2)CHG = 1.0 
! 
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!     FILTER 3....STRONG FILTER.....EFFECTIVE FOR HURRICANE GILBERT
!
!     N = 17 .... AND m VARIES AS 2,3,4,2,5,6,7,2,8,9,10,2,11,2,2,2,2
!
!                  
      IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.EQ.12.OR.KTY.GE.14.).AND.IFL.EQ.3)CHG = 1.0 
!
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!
!     FILTER 4.....VERY STRONG FILTER....THE PATTERN WILL START TO BECOME ZONAL....
!
!
!     N = 24.......AND m VARIES AS : 
!                  2,3,4,2,5,6,7,2,8,9,10,2,11,12,13,2,2,2,2,2,2,2,2,2
!
!
!                    
      IF((KTY.EQ.4.OR.KTY.EQ.8.OR.KTY.EQ.12.OR.KTY.GE.16).AND.IFL.EQ.4)CHG = 1.0  
!  
!
!
!***************************************************************************************8
!
!                 
      IF(CHG.EQ.0)KT = KT + 1                                            
      IF(CHG.EQ.1.0)TK(KTY) = .25                                        
      IF(CHG.EQ.1.0) GO TO 801                                            
      FACT = 2.0*PI/(FLOAT(KT) + 1.0)                                    
      TK(KTY) = -.5/(COS(FACT) - 1.0)                                    
      DO 679 NA = 2 , 25                                                  
        AMPF(NA) = 1 + 2.*TK(KTY)*(COS(2.*PI/FLOAT(NA)) - 1.0)              
679   CONTINUE                                                            
!                                                                       
!                                                                       
801   CONTINUE                                                           
      CHG = 0.0                                                          
802 CONTINUE 
!
!       WRITE(6,815) (TK(KK),KK = 1 , NTY)                                 
!815    FORMAT(2X,'THIS IS TK:',E12.6)                                     
! 
!
! ......THIS IS THE END OF THE KTY LOOP....................... 
!                                                        
! 
!**********DESMOOTHING IS SET UP IF NEEDED******************** 
!                                               
!                                                                        
    IF(ISMTH.EQ.1) THEN                                                 
      NTYM = NTY - 1                                                   
      TFF = 1.0                                                        
      DO 61 K = 1 , NTYM                                               
        TFF = TFF*(1. + 2.*TK(K)*COSF)                                   
61    CONTINUE                                                         
      TFR = 1./TFF                                                     
      TK(NTY) = (TFR - 1.0)/(2.*COSF)                                  
      WRITE(6,816) TK(NTY)                                             
816   FORMAT(2X,'THE DESMOOTHING CONSTANT',E12.6)                      
    END IF                                                              
!                                                                        
!***********PRINT OUT THE DAMPING CHARACTERISTICS*************** 
! 
    IRT = KT+1                                                                   
!                                                                        
    DO 610 NZ = 2 , 40                                                  
      AMP = 1.0                                                           
      TNN = FLOAT(NZ)                                                     
      CKG = 0.0                                                           
      IF(NZ.GT.IRT)CKG = 1.0                                                
      DO 617 KT = 1,NTY                                                 
        AMP1 = (1. + 2*TK(KT)*(COS(2.*PI/TNN)-1.0))                         
        IF(CKG.EQ.0.0) GO TO 619                                             
        AMP =  AMP1*AMP                                                     
619     IF(ABS(AMP1).LT..01) CKG = 1.0                                       
617   CONTINUE
      AMM = AMP
      IF(NZ.LE.IRT) AMM=0.0 
!
      ZZ = FLOAT(NZ) 
      WRITE(11,455)ZZ,AMM   
455   FORMAT(F8.3,F8.3) 
!
!
!
!  THE FOLLOWING WRITE STATEMENT WILL LET YOU KNOW THE AMOUNT OF THE
!  WAVE THAT HAS REMAINED AFTER THE FILTERING,
!  FOR THE WAVE OF A GIVEN LENGTH D (WHICH IS CURRENTLY ONE DEGREE)
!
      IF(NZ.EQ.20.OR.NZ.EQ.30.OR.NZ.EQ.40) THEN
        WRITE(6,677)NZ,AMM                                                  
677     FORMAT(2X,'WAVE NUMBER',I5,2X,'PERCENT WAVE REMAINING',E12.6)                    
      END IF   
!                                                                    
610 CONTINUE                                                            
! 
!*******************************************************************
!                                                                       
!        DO THE SMOOTHING IN THE LATITUDINAL DIRECTION:
!                         (EQUATION A1)
!
    DO 600 J = 1 , JMX                                                                       
      DO 58 NN = 1 , NTY                                                
        XTU(1,NN)   = U(1,J)                                 
        XTU(IMX,NN) = U(IMX,J)                              
        XTV(1,NN)   = V(1,J)                                 
        XTV(IMX,NN) = V(IMX,J)                              
58    CONTINUE                                                          
!                                                                     
      DO 60 I = 2,IMXM
        XTU(I,1) = U(I,J)   + TK(1)*(U(I-1,J) +U(I+1,J) - 2.*U(I,J))                      
        XTV(I,1) = V(I,J)   + TK(1)*(V(I-1,J) +V(I+1,J) - 2.*V(I,J))
60    CONTINUE                                                          
!                                                                      
      DO 65 NN = 2 , NTY                                                
        DO 62  I = 2 , IMXM                                              
          XTU(I,NN) = XTU(I,NN-1) + TK(NN)*(XTU(I-1,NN-1) +XTU(I+1,NN-1) - 2.*XTU(I,NN-1)) 
          XTV(I,NN) = XTV(I,NN-1) + TK(NN)*(XTV(I-1,NN-1) +XTV(I+1,NN-1) - 2.*XTV(I,NN-1)) 
62      CONTINUE                                                            
65    CONTINUE                                                            
!                                                                      
      DO 70 I = 1,IMX                                                 
        US(I,J)   = XTU(I,NTY)                                           
        VS(I,J)   = XTV(I,NTY)                                           
70    CONTINUE                                                           
600 CONTINUE                                                                       
!                                                             
!                                                                        
!********************************************************************    
!                                                                        
!    NOW DO THE SMOOTHING IN THE MERIDIONAL DIRECTION:                   
!                         (EQUATION A3)     
!   
!
    DO 700   I = 1 , IMX   
      DO 80 NN = 1 , NTY                                                
        YTU(1,NN)   = US(I,1)                                           
        YTU(JMX,NN) = US(I,JMX) 
        YTV(1,NN)   = VS(I,1)                                           
        YTV(JMX,NN) = VS(I,JMX)   
80    CONTINUE                                                          
!                                                                      
      DO 90 J = 2 , JMXM                                               
        YTU(J,1) = US(I,J) + TK(1)*(US(I,J-1) + US(I,J+1)-2.*US(I,J)) 
        YTV(J,1) = VS(I,J) + TK(1)*(VS(I,J-1) + VS(I,J+1)-2.*VS(I,J)) 
90    CONTINUE                                                          
!                                                                      
      DO NN = 2 , NTY                                                
        DO J  = 2 , JMXM                                              
          YTU(J,NN) = YTU(J,NN-1) + TK(NN)*(YTU(J-1,NN-1)+YTU(J+1,NN-1) - 2.*YTU(J,NN-1))
          YTV(J,NN) = YTV(J,NN-1) + TK(NN)*(YTV(J-1,NN-1)+YTV(J+1,NN-1) - 2.*YTV(J,NN-1)) 
        END DO
      END DO
!
!
!   STORE THE FILTERED FIELDS IN US,VS AND GS 
!
!                                                                        
      DO 99 J = 1 , JMX                                                
        US(I,J)   =  YTU(J,NTY)                                          
        VS(I,J)   =  YTV(J,NTY)                                          
99    CONTINUE                                                          
!                                                                        
!                                                                        
700 CONTINUE                                                          
    RETURN                                               
  END SUBROUTINE PHASE