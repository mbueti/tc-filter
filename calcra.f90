  SUBROUTINE CALCRa(RO,RTAN,iang)
    COMMON /WINDS/ DMMM(IMX,JMX,2),TANG(IMX,JMX),DEL(IMX,JMX),THA(IMX,JMX),XF(IMX,JMX),DS(IMX,JMX)
!
    COMMON  /TOTAL/ DDEL,DTHA
    COMMON  /COOR/ XV,YV,XOLD,YOLD,XCORN,YCORN,FACTR,id1,id2
!
    FACT =  COS(YOLD)
!
    DX=DDEL/PI180
    DY=DTHA/PI180
    XC = (XOLD-XCORN)*DX
    YC = (YOLD-YCORN)*DY

!      
    THETA= 2.*PI*FLOAT(iang)/FLOAT(NMX)
        ! X=RO*fact*COS(THETA)+XC
    X=RO*COS(THETA)+XC
    Y=RO*SIN(THETA)+YC
    IX=mod(FLOOR(X/DX),imx)
    if (IX.eq.0) IX=imx
    IY=mod(FLOOR(Y/DY),jmx)
    if (IY.eq.0) IY=jmx
    IX1=mod(IX+1,imx)
    if (IX1.eq.0) IX1=imx
    IY1=mod(IY+1,jmx)
    if (IY1.eq.0) IY1=jmx
    P=X/DX-FLOAT(IX)
    Q=Y/DY-FLOAT(IY)
!        rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,MODULO(IY+1, IY))
!      1      +  (1.-Q)*P*XF(max(MODULO(IX+1, IX),1),IY) +
!      1      P*Q*XF(max(MODULO(IX+1,IX),1),max(MODULO(IY+1,IY),1))
      !   print *, XF(max(MODULO(IX+1,IX),1),max(MODULO(IY+1,IY),1))
    rtan=(1.-P)*(1.-Q)*XF(IX,IY) +(1.-P)*Q*XF(IX,IY1) +  (1.-Q)*P*XF(IX1,IY)+P*Q*XF(IX1,IY1)
    RETURN
  END SUBROUTINE CALCRa