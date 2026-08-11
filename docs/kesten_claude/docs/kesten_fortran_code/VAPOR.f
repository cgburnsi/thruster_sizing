		SUBROUTINE VAPOR(TEMP,ZV,LL,Q1,Q,H)									   0
		REAL KP,K															  10
		REAL MBSS 															  20
		INTEGER PRINT														  30
		COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  40
	   1        SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  50
		COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  50
	   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2	  60
	   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 							  70
		COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								  80
		COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			 100
		COMMON /MUVST/VISVST(30)											 110
		COMMON /FLAGS/MFLAG,KFLAG,PRINT 									 120
		COMMON /IFCE00/IFC,GATZ0											 130
		COMMON /BBBB/DP3,A,KC3,K0,XOA,CPS,CI3,GAMMA,BETA					 140	
		COMMON /CCC/H4TBL(40),H3TBL(40)										 150	
		COMMON /DDD/CFTBL4(34),CFTBL3(34),CFTBL2(34),CFTBL1(34)				 160
		JZ=LL 																 170
		Z(LL)=ZV 															 180
		TEMP=TVAP 															 190
		PFRC3D=0. 															 200
		ZBOUND=ZEND															 210
		KFLAG=0																 220
		JFLAG=0																 230
		KOUNT=0																 240
		INJECT=0															 250
		N=7																	 260
		THREE=3. 															 270
		WRITE(6,100)														 280
	100 FORMAT(109H1 ******************************************** ENTE       290
	   1RING VAPOR REGION    *********************************  //) 		 300
		CALL CONC (C1,C2,C3,C4,Z(LL),H,TEMP)								 310
		SUM=C1/WM1+C2/WM2+C3/WM3+C4/WM4										 320
		FRAC1=C1/(WM1*SUM)													 330
		FRAC2=C2/(WM2*SUM)													 340
		FRAC3=C3/(WM3*SUM)													 350
		FRAC4=C4/(WM4*SUM)													 360
		FRAC3D=(FRAC1/FRAC2-1.0)/(3.-FRAC1/FRAC2)							 370
		WRITE(6,4059)														 380
		WRITE(6,4050) Z(LL),TEMP,PRES,H,C1,C2,C3,C4 						 390
		WRITE(6,37)															 400
		WRITE(6,38) FRAC1,FRAC2,FRAC3,FRAC4,FRAC3D							 410
   3990 CALL UNBAR (SHTBL1(1),1,TEMP,0.,CP4,KK) 							 420   (This needs checked.  It says SHTBL1, but uses CP4)
        CALL UNBAR (SHTBL2(1),1,TEMP,0.,CP3,KK) 							 430   (This needs checked.  It says SHTBL2, but uses CP3)
        CALL UNBAR (SHTBL3(1),1,TEMP,0.,CP2,KK) 							 440   (This needs checked.  It says SHTBL3, but uses CP2)
        CALL UNBAR (SHTBL4(1),1,TEMP,0.,CP1,KK) 							 450   (This needs checked.  It says SHTBL4, but uses CP1)
		CAV=(C4*CP4+C3*CP3+C2*CP2+C1*CP1)/(C4+C3+C2+C1)						 460
		WMAV=(C1+C2+C3+C4)/(C1/WM1+C2/WM2+C3/WM3+C4/WM4)					 470
		CALL UNBAR (TBLH4(1),1,TEMP,0.,H4,KK)								 480
		CALL UNBAR (ZTBLD(1),1,Z(LL),0.,DELA,KK)							 490
		CALL UNBAR (ZTBLAP(1),1,Z(LL),0.,AP,KK)								 500
		CALL UNBAR (ZTBLA(1),1,Z(LL),0.,A,KK)								 510
		CALL PARAM(TEMP,Z(LL),1,C4,H4,0,G,GMMA,K,DPA,BETA)					 520
		IF(C4)7,29,7														 530
      7 DIFN=DIF4*((TEMP/492.)**1.823)*14.7/PRES							 540
		CALL UNBAR (VISVST(1),1,TEMP,0.,VIS,KK)								 550
		RHO=PRES*WMAV/(R*TEMP)												 560
		AKC=.61*G/RHO*((VIS/(RHO*DIFN))**-.667)*((G/(AP*VIS))**-.41)		 570
		DERIV(LL)=AKC*C4/DPA												 580
	  6 T4=AP*DPA*DERIV(LL)													 590
		GO TO 31															 600
	 29 T4=0. 																 610
     31 CALL UNBAR (TBLH3(1),1,TEMP,0.,H3,KK)								 620
		CALL PARAM(TEMP,Z(LL),1,C3,H3,0,G,GMMA,K,DPA,BETA)					 630
		IF(C3)13,30,13														 640
	 13 CALL SGRAD (GRAD,TGRAD) 											 650
		DERIV(LL)=GRAD/DPA													 660
		T3=AP*DPA*DERIV(LL)													 670
		GO TO 32 															 680
     30 T3=0. 																 690
	 32 RHOM = ALPHA3*C4*EXP(-CGM/TEMP)										 700
		T1=PRES*WMAV/(R*TEMP*G) 											 710
		T2=RHOM*DELA														 720
		DHDZ(LL)=-H4/G*(T2+T4)-H3/G*T3-FC/G*(H-HF)							 730
C		KFLAG IS SLOPE INDICATOR											
C		=0 --- SLOPE MOVING TOWARDS FIRST PEAK
C		=1 --- SLOPE HAS REACHED FIRST PEAK
		IF(KFLAG.EQ.1)GO TO 4000											 740
		Z1=-H4/(ENMX3*DHDZ(LL)) 											 750
		Z2=.05*(ZEND-Z(LL))													 760
		IF(DHDZ(LL)*(1.-Z1/Z2))4060,4060,15									 770
     15 DZ=Z1 																 780
   4000 DTDZ=DHDZ(LL)/CAV													 790
		W1=C1/RHO															 800
		W2=C2/RHO															 810
		W3=C3/RHO															 820
		W4=C4/RHO															 830
		S1=1./G																 840   (I think this is S1=1./G)
		S5=FC/(G*RHO)														 850  (I think this is correct.)
		DW4DZ=S1*(FC-T2-T4)-C4*S5											 860
		DW3DZ=S1*(T2*WM3/WM4+T4*WM3/WM4-T3)-C3*S5							 870
		DW2DZ=S1*(.5*T2*WM2/WM4+.5*T4*WM2/WM4+.5*T3*WM2/WM3)-C2*S5			 880
		DW1DZ=S1*(.5*T2*WM1/WM4+.5*T4*WM1/WM4+1.5*T3*WM1/WM3)-C1*S5			 890
		SUMWM=W1/WM1+W2/M2+WM3+WM4/WM4 										 900
		SUMDWDZ=DW1DZ/WM1+DW2DZ/WM2+DW3DZ/WM3+DW4DZ/WM4 		 			 910
		DMDZ=-WMAV/SUMWM*SMDWDZ												 920  (The .pdf scanned file is hard to read, but I think this is correct)
		DPDZ = (DELA-1.)/DELA**3*(1.75+75.*VIS*(1.-DELA)/(A*G))*G**2/ 		 930
	   1      (64.4*A*RHO) 													 940
		DPDZ = DPDZ/144.													 950
		DRDZR = DMDZ/WMAV-DTDZ/TEMP+DPDZ/PRES								 960
		T5=FC/G-DRDZR 														 970
		DC4DZ=T1*(FC-T2-T4)-C4*T5 											 980
		DC3DZ=T1*(T2*WM3/WM4+T4*WM3/WM4-T3)-C3*T5							 990
		DC2DZ=T1*(.5*T2*WM2/WM4+.5*T4*WM2/WM4+.5*T3*WM2/WM3)-C2*T5			1000
		DC1DZ=T1*(.5*T2*WM1/WM4+.5*T4*WM1/WM4+1.5*T3*WM1/WM3)-C1*T5			1010
		IF(KFLAG.EQ.0)GO TO 16 												1020
C		JFLAG IS DZ INDICATOR FOR NON-ZERO FEED RATE CASES
C		=0 --- DZ INCREMENT O.K.
C		=1 --- INCREMENT INITIALLY TOO SMALL


		IF(JFLAG.EQ.1)GO TO 93												1030
	 90 IF(KOUNT.EQ.4.OR.KOUNT.EQ.6.OR.KOUNT.EQ.8.OR.KOUNT.EQ.10.OR.KOUNT.  1040
       X         EQ.12.OR.KOUNT.EQ.14)N=N-1									1050
		KOUNT=KOUNT+1														1060
		DZ=DELTAZ/(THREE**N)												1070


		IF(FC)98,98,93														1080
     98 IF(IFC.EQ.0)GO TO 16 												1090
C		IF FEED RATE IS NON-ZERO WE MUST MAKE ADDITIONAL CHECKS ON STEP
C		SIZE OF Z
	 93 IF(ABS(DTDZ)*DZ.GT..01*TEMP)GO TO 19 								1100
C		CHECK IF WE HAVE REACHED THE END OF THE INJECTOR
		IF((1.+DZ/Z(LL)-Z0)+.01*Z0/ABS(Z(LL))-Z0)).GT.0.)GO TO 16 			1110
	 19 DZ1=DZ																1120
		CALL REDIVD (DZ1,DTDZ,NINT,JFLAG,I,LL)								1130
		WRITE (6,91) KOUNT, NINT											1140
     91 FORMAT (//7H KOUNT=I2 ,37H --- THIS INTERVAL HAS BEEN REDIVIDEDI4,  1150
       X6H TIMES)															1160
     16 H=H+DHDZ(LL)*DZ														1170
		IF (H.LT.HV) GO TO 106												1180
   4051 TEMP=TEMP+DTDZ*DZ													1190
		PRES = PRES+DPDZ*DZ													1200
		C4=C4+DC4DZ*DZ														1210
		C3=C3+DC3DZ*DZ														1220
		C2=C2+DC2DZ*DZ														1230
		C1=C1+DC1DZ*DZ														1240
	400 SUM=C1/WM1+C2/WM2+C3/WM3+C4/WM4										1250
		IF(C4.LT.0.)SUM=SUM-C4/WM4											1260
		IF(C3.LT.0.)SUM=SUM-C3/WM3											1270
		FRAC1=C1/(WM1*SUM)													1280
		FRAC2=C2/(WM2*SUM)													1290
		FRAC3=C3/(WM3*SUM)													1300
		FRAC4=C4/(WM4*SUM)													1310
		FRAC3D=(FRAC1/FRAC2-1.)/(3.-FRAC1/FRAC2)							1320
		IF(KFLAG.EQ.1)GO TO 500												1330
C		IF RELATIVE DIFFERENCE OF SUCCESSIVE FRAC3D'S IS GREATER THAN 5 
C		PERCENT WE RECALCULATE WITH SMALLER DZ INCREMENT
		IF((FRAC3D-PFRC3D).LT..05)GO TO 500									1340
	 17 H=H-DHDZ(LL)*DZ														1350
	 18 TEMP=TEMP-DPDZ*DZ													1360
		PRES = PRES-DPDZ*DZ													1370
		C4=C4-DC4DZ*DZ														1380
		C3=C3-DC3DZ*DZ														1390
		C2=C2-DC2DZ*DZ														1400
		C1=C1-DC1DZ*DZ														1410
		DZ=DZ/2.															1420
		GO TO 16															1430
	500 PFRC3D=FRAC3D														1440
		IF(C4)70,71,71														1450
     70 C4=0.																1460
		FRAC4=0.															1470
	 71 IF(C3)72,73,73														1480
	 72 C3=0.																1490
		FRAC3=0.															1500
	 73 LL=LL+1																1510
		Z(LL)=Z(LL-1)+DZ													1520
		IF(JFLAG.EQ.1)NINT=NINT-1											1530
		IF(NINT.EQ.0)JFLAG=0												1540
		WRITE (6,4059)														1550
   4059 FORMAT (//121H   Z   TEMP   PRES									1560
       X    H   C1    C2    C3												1570
       X    C4    )															1580
		WRITE (6,4050) Z(LL),TEMP,PRES,H,C1,C2,C3,C4						1590
   4050 FORMAT (1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E  1600
       X15.8,1X,E15.8///)													1610
		WRITE (6,37)														1620
	 37 FORMAT (98H    MFRAC1   MFRAC2 										1630
	   1  MFRAC3   MFRAC4   FRAC3D ) 										1640
		WRITE (6,38) FRAC1,FRAC2,FRAC3,FRAC4,FRAC3D							1650
	 38 FORMAT (20X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8,1X,E15.8 ///) 			1660
		IF(Z(LL).GT.ZBOUND)GO TO 1											1670
		IF(KOUNT.GT.15.AND.JFLAG.EQ.0) GO TO 1								1680
		INJECT=0															1690
		GO TO 3990															1700
C		SLOPE HAS TURNED NEGATIVE --- NEXT INCREMENTS ARE (1.3)**7,
C		(1/3)**6,...(1/3)**1 OF ZEND-Z(LL)
C		15 INCREMENTS IN TOTAL
   4060 DELTAZ=ZEND-Z(LL)													1710
		ZBOUND=ZEND+DELTAZ/3.												1720
		KFLAG=1																1730
		KOUNT=KOUNT+1														1740
		GO TO 4000															1750
	  1 IF(IFC.EQ.1)GO TO 22												1760
C		FOR ZERO FEED RATE CASES WE ADD IN 6 ADDITIONAL Z'S FOR USE IN
C		TRANSIENT MODEL
		ZGAP=Z(LL)-Z(LL-2)													1770
		DZZ=ZGAP/8.															1780
		JJ1=LL-1															1790
		JJ8=LL+8															1800  (I think this is correct, but it's hard to tell if these are 6 or 8 here)
		DO 27 I=JJ1,JJ8														1810
		Z(I)=Z(I-1)+DZZ														1820
     27 CONTINUE															1830
		LL=LL+6																1840
C		PRINT OUT VAPOR REGION Z VALUES FOR USE IN TRANSIENT MODEL
	 22 WRITE (6,62)														1850
	 62 FORMAT (1H1,54X,22H Z'S FROM VAPOR REGION)							1860
		WRITE (6,63) (Z(I),I=JZ,LL) 										1870  (That looks like I=JZ,LL), but I don't recall a JZ anywhere)
	 63 FORMAT (1X,10E13.7)													1880
		MBSS = (C1+C2+C3+C4)/(C1/WM1+C2/WM2+C3/WM3+C4/WM4)					1890  (Oh boy this is hard to see.  I think it is MBSS, but I can't tell)  I think it means MBAR at Steady State
		WRITE (6,8001) MBSS,G												1900   (yeah.  still not sure if that is MBSS or not)
   8001 FORMAT (///42X'STEADY STATE VALUES FOR MBAR AND G AT END OF BED'/   1910
       X           47X'MBAR ='E12.5,5X,'G ='E12.5)							1920
		GO TO 999															1930
	106 WRITE (6,107)														1940
	107 FORMAT (////13X'THERE IS A PUDDLE OF COLD HYDRAZINE AT THE LIQUID-  1950
	   XVAPOR/VAPOR INTERFACE --- TRY USING A LARGER VALUE FOR G0') 		1960
	999 RETURN																1970
		END																	1980

		













		

		

		
