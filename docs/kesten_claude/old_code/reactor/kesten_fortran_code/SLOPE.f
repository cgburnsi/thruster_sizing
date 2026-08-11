		SUBROUTINE SLOPE (CG,GMMA,K,BETA,EN12,RATE,DPA,A,DIFF)				   0
		REAL KP,K															  10
		INTEGER PRINT														  20
		DIMENSION K(1),     GMMA(1),BETA(1)DPA(1),CG(1),S(210),FST(210), 	  30
	   1SEC(210),D(210),CPA(210),CPB(21),E(210),A(1),TERM1(210), 			  40
       2TERM2(210),C(210),XXX(210),XOA(210) 								  50
		COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  60
	   1        SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  70
		COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  80
	   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2	  90
	   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 							 100
		COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								 110
		COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			 120
		COMMON /MUVST/VISVST(30)											 130
		COMMON /FLAGS/MFLAG,KFLAG,PRINT 									 140
		IFLAG=0 															 150
		JFLAG=0																 160
		FRAC=.99															 170
	 12 B=0.0 																 180
	702 I=1 																 190
		N0T=0 																 200
		ADIV=100.															 210
		BDIV=100.															 220
		TOL=.01 															 230
		STORE=1.0 															 240
		KJ=0 																 250    (This needs checked.  it might be K0=0)
		HOLD=0.																 260
		NI=0																 270
		MI=0																 280
		IL=0 																 290
		IK=0 																 300
		MM=0																 310
		IM=0																 320
		F=0.5																 330
	 20 MM=MM+1 															 340
		DR=B/ADIV 															 350
		BDR=(A(I)-B)/BDIV 													 360
		JINT=BDIV															 370
		INIT=ADIV															 380   
		IF  (MM.EQ.1) GO TO 15 												 390
		GO TO 16															 400
	 15 MMAX=JINT+1 														 410
		MAX=MMAX-1 															 420
		MINT=0																 430
		GO TO 17 															 440
	 16 MMAX=INIT+JINT+1 													 450
		MINT=INIT+1															 460
		MAX=MMAX-1															 470
	 17 X=0.0 																 480
		DO 40 IJ=1,MMAX														 490
		IF(IJ.GT.MINT) GO TO 49 											 500
		X=FLOAT(IJ-1)*DR													 510
		GO TO 23															 520
	 49 KMJ=IJ-MINT															 530
		X=B+FLOAT(KMJ)*BDR 													 540
	 23 CP=CG(I)*X/(A(I)-B)-(B/(A(I)-B))*CG(I) 								 550
		IF(CP.LT.0.0) CP=0.0												 560
C		THIS IS THE GENERAL EQUATION FOR A LINE WITH NEGATIVE Y-INTERCEPT
	 24 FACT1=K(I)*(CG(I))**(1.-EN12)*(CP)**EN12							 570
		FACT2=EXP(GMMA(I)*BETA(I)*(1.-CP/CG(I))/(1.+BETA(I)*(1.-CP/CG(I)))   580
	   1)                 													 590
		IF(X) 37,40,37 														 600
	 37 RCP=FACT1*FACT2 													 610
     35 S(IJ)=(1./X-1./A(I))*X*X*RCP 										 620    (I think this is line number 35, but it's hard to read)
	 40 CONTINUE 															 630
		GO TO 50															 640
	 50 SUM=0.0 															 650
		SUMA=0.0 															 660
		DO 60 IL=2,MAX 														 670
		IF(IL.GT.MINT) GO TO 61 											 680
		SUM=SUM+S(IL)														 690
		GO TO 60 															 700
	 61 SUMA=SUMA+S(IL) 													 710
	 60 CONTINUE 															 720
		S(1)=0. 															 730
		SGMA=(S(1)+2.0*SUM)*(DR/(2.0*DPA(I)))+(S(MMAX)+2.0*SUMA)*(BDR/(2. 	 740
	   10*DPA(I))) 															 750
		D2=SGMA-FRAC*CG(I)													 760
		TOT=D2+TOL*CG(I) 													 770  
		IF(MM-1) 110,110,112												 780
	112 IF(D2) 1,1,2 														 790
	  1 IF(TOT) 2,230,230 													 800
	  2 GO TO 115 															 810
	110 IF(IFLAG.EQ.1) GO TO 115 											 820
		IF(D2) 150,115,115													 830
	115 IF(MM-1) 20,120,140 												 840
	120 XOLD=B																 850
		IF (PRINT.EQ.1) WRITE(6,125)										 860
	125 FORMAT(41H0INITIAL CHOICE THRU ORIGIN IS TOO LARGE )				 870
		D1=D2																 880
C		CHANGE THE EQUATION OF THE LINE
		B=.999999    *A(I) 													 890
		GO TO 20															 900
C   	USE PRECCEDING RESULTS TO ESTIMATE B FOR THE NEW LINE
	140	TJMP=B 																 910   
		MI=MI+1 															 920
		B=B+(D2/(D2-D1))*(XOLD-B)											 930
   1741 XOLD=TJMP															 940
		D1=D2 																 950
	143 IF(MI-20) 145,145,147												 960
	145 GO TO 20 															 970
	147 IF (JFLAG.EQ.1) GO TO 230 											 980
		B=.9*A(I) 															 990  
		JFLAG=1																1000
		GO TO 702															1010











C		INITIAL CHOICE THRU ORIGIN IS SATISFACTORY
	150 KJ=1																1020
		IFLAG=1																1030
		X=0.																1040
C		CALCULATE THE VALUE OF THE TWO INTEGRALS FOR ALL DR (101 POINTS)
		DO 170 II=1,MMAX													1050
		X=FLOAT(II-1)*BDR													1060
		CP=CG(I)*X/(A(I)-B)-(B/(A(I)-B))*CG(I)								1070
		TERM1(II)=K(I)*(CG(I))**(1.-EN12)*(CP)**EN12						1080
		TERM2(II)=EXP(GMMA(I)*BETA(I)*(1.-CP/CG(I))/(1.+BETA(I)*(1.-CP/CG(  1090
       1I))) 																1100
		XXX(II)=X															1110
		XOA(II)=XXX(II)/A(I)												1120
		RCP=TERM1(II)*TERM2(II)												1130
		FST(II)=X*X*RCP														1140
		IF(X) 165,170,165													1150
	165 SEC(II)=(1./X-1./A(I))*X*X*RCP										1160
	170 CONTINUE															1170





C		THE TRAPEZOIDAL RULE IS USED TO EVALUATE BOTH INTEGRALS
	172 C(1)=0.																1180 (I think this is C(1) and not C(I))
	173 DO 175 JJ=2,MMAX 													1190
		IF(JJ.GT.MINT) GO TO 176 											1200
		C(JJ)=C(JJ-1)+(FST(JJ)+FST(JJ-1))*(DR/(2.0*DPA(I)))					1210
		GO TO 175															1220
	176	C(JJ)=C(JJ-1)+(FST(JJ)+FST(JJ-1))*(BDR/(2.0*DPA(I)))				1230
	175 CONTINUE															1240
		SEC(1)=0.															1250
		D(1)=0.																1260
		DO 180 KK=2,MMAX													1270
		IF(KK.GT.MINT) GO TO 179 											1280
		D(KK)=D(KK-1)+(SEC(KK-1)+SEC(KK))*(DR/(2.0*DPA(I)))					1290
		GO TO 180															1300
	179	D(KK)=D(KK-1)+(SEC(KK-1)+SEC(KK))*(BDR/(2.0*DPA(I)))				1310
	180 CONTINUE															1320
C		THE VALUE OF CP AT X=0 IS CG -D(101)
		E(1)=D(MMAX)														1330
	186 CPA(1)=0.0 															1340
C		NEGATIVE VALUES OF CPA(1)=CP(0) CANNOT BE USED
C 		STORE THE CPA(1) VALUE IN CASE A NEW F FACTOR MUST BE USED
	184 STORE=CPA(1)														1350
		IF(KJ.EQ.1)  GO TO 185 												1360
		IF(IM.EQ.0) GO TO 185 												1370
		IF(NI.EQ.1) GO TO 368 												1380
	185 SAM=0.0																1390
		DO 190 LL=2,MMAX													1400
		IF(LL.GT.MINT) GO TO 188 											1410
		E(LL)=E(LL-1)-(SEC(LL)+SEC(LL-1))*(DR/(2.0*DPA(I)))					1420
		SAM=FLOAT(LL-1)*DR													1430
		GO TO 189 															1440
	188 E(LL)=E(LL-1)-(SEC(LL)+SEC(LL-1))*(BDR/(2.0*DPA(I)))				1450
		SAM=SAM+BDR															1460
	189 CPA(LL)=CG(I)-((1./SAM)-(1./A(I)))*C(LL)-E(LL) 						1470
		IF(LL.LT.MINT) CPA(LL)=0.0 											1480
		IF(CPA(LL).LT.0.0) CPA(LL)=0.0										1490
	190 CONTINUE															1500
		IF(KJ.EQ.1) GO TO 280 												1510
		X=0.																1520
		IF(IM.EQ.0) GO TO 250												1530
		IK=1																1540
		GO TO 280															1550
C		THE NEXT ITERATION USES THE VALUES OF CP JUST CALCULATED
	192 DO 200 LI=1,MMAX													1560
		IF(LI.GT.MINT) GO TO 195											1570
		X=FLOAT(LI-1)*DR													1580
		GO TO 199 															1590
	195 KLK=LI-MINT															1600
		X=B+FLOAT(KLK)*BDR													1610
	199 CONTINUE															1620
		TERM1(LI)=K(I)*(CG(I))**(1.-EN12)*(CPA(LI))**EN12					1630
		TERM2(LI)=EXP(GMMA(I)*BETA(I)*(1.-CPA(LI)/CG(I))/(1.+BETA(I)*(1.-   1640
	   1CPA(LI)/CG(I)))) 													1650
		XXX(LI)=X															1660
		XOA(LI)=XXX(LI)/A(I)												1670
		RCP=TERM1(LI)*TERM2(LI)												1680
		FST(LI)=X*X*RCP														1690
		IF(X) 200,200,193													1700
	193 SEC(LI)=1./X-1./A(I))*X*X*RCP										1710
	200 CONTINUE															1720
		GO TO 172 															1730
C		THIS IS FOR THE CASE WHERE THE INITIAL GUESS WAS TOO LARGE
	230 IF(PRINT.EQ.1) WRITE(6,235)MI,B										1740
	235 FORMAT(41H SATISFACTORY STARTING CURVE FOUND AFTER I2,31H TRIALS T  1750
	   1HE VALUE OF B (X0) IS E12.7)										1760
		IF(B.GT..998*A(I))DERIF=2.6*CG(I)/(A(I)-B)							1770
		IF(B.GT..998*A(I)) GO TO 322										1780


C		THE RANGE OF X FROM ZERO TO A IS USED
	237 X=0.																1790
		NI=1																1800
		DO 240 NN=1,MMAX													1810
		IF(NN.GT.MINT) GO TO 246 											1820
		X=FLOAT(NN-1)*DR													1830
		GO TO 354 															1840
	246 KJK=NN-MINT 														1850
		X=B+FLOAT(KJK)*BDR													1860
	354 IF(IM.GE.1) GO TO 353 												1870
	238 CPB(NN) = CG(I)*X/(A(I)-B) - (B/A(A(I)-B))*CG(I)					1880
		IF(CPB(NN).LT.0.0) CPB(NN)=0.0										1890
	353 TERM1(NN)=K(I)*(CG(I))**(1.-EN12)*(CPB(NN))**EN12					1900
		TERM2(NN)= EXP(GMMA(I)*BETA(I)*(1.-CPB(NN)/CG(I))/(1.+BETA(I)*(1.-  1910
	   1CPB(NN)/CG(I))))													1920
		XXX(NN)=x															1930
		XOA(NN)=XXX(NN)/A(I)												1940
		RCP = TERM1(NN)*TERM2(NN)											1950
		FST(NN)=X*X*RCP														1960
		IF(X) 240,240,247													1970
	247 SEC(NN) = (1./X-1./A(I))*X*X*RCP									1980
	240 CONTINUE															1990
		GO TO 172															2000


	368 DO 370 NL=1,MMAX													2010
		CPB(NL)=CPA(NL)														2020
	370 CONTINUE															2030
		GO TO 185 															2040
	250 DO 260 IL=1,MMAX													2050
	252 CPA(IL)=F*CPA(IL)+(1.-F)*CPB(IL)									2060
	260 CONTINUE															2070
		IK=0																2080
C		THE VALUES OF X AT A AND NEAREST A ARE USED IN FINDING THE DERIVATIVE
	280 DERIF=(CG(I)-CCPA(MAX))/BDR											2090
		IM=IM+1																2100
		IF(IM.GT.99.AND.IFLAG.EQ.1) GO TO 701								2110
		IF(IM.GT.99) GO TO 328												2120
	310 IF(KJ.EQ.1) IK=1													2130
		IF(ABS(DERIF-HOLD)- 0.05*DERIF) 3,3,321 							2140
	  3 IF(IK.EQ.1) GO TO 322 												2150 
	321 HOLD=DERIF															2160
		IF(KJ.EQ.1) GO TO 192 												2170 
		IF(IK.EQ.1) GO TO 250												2180
		GO TO 192 															2190
	322 RATE=DERIF															2200
		IF(C4-CG(I)) 777,4,777												2210
	  4 IF(H-HL) 777,777,5 													2220 (I think this says (H-HL)
      5 DIFN=DIFF*((TEMP/492.)**1.823*14.7/PRES								2230
		CALL UNBAR (VISVST(1),1,TEMP,0.,VIS,KK)								2240
		RHO=PRES*WMAV/(R*TEMP)												2250
		AKC=(.61*G)/(RHO)*((VIS/(RHO*DIFN))**-.667*((G/(AP*VIS))**-.41) 	2260
		RAT=AKC*CG(I)/DPA(I)												2270
		IF(RATE-RAT) 778,778,6												2280
	  6 MFLAG=1																2290
		RATE=RAT															2300
		GO TO 777 															2310
	778 MFLAG=0																2320
	777 IF(B.GT..998*A(I)) GO TO 888										2330
		IF (PRINT.EQ.1) WRITE(6,182)IM										2340
	182 FORMAT (11H ITERATION=,I3)											2350
		IF (PRINT.EQ.1) WRITE (6,336)										2360
	336 FORMAT(13X,105M X/A       CPA       X/A     CPA						2370
	   1     X/A     CPA    X/A    CPA      ) 								2380
		IF (PRINT.EQ.1) WRITE(6,183)(XOA(I),CPA(I),I=1,MMAX,4)				2390
	183 FORMAT(9X,L12.7,1X,L12.7,1X,E12.7,1I,E12.7,1X,E12.7,1X,E12.7,1X,E1  2400
	   12.7,1X,E12.7) 														2410
	888 IF (PRINT.EQ.1) WRITE(6,323)RATE									2420
	323 FORMAT(/24J THE SLOPE CONVERGES TO E12.7)							2430
	327 RETURN																2440
	328 PERC=(ABS(DERIF-HOLD)/DERIF)*100.									2450
		ALL = 5.															2460
		GO TO 322															2470
C		SIMPLIFIED VERSION DOES NOT ''CONVERGE'' IN 99 ITERATIONS
C		SET B=.000001*A AND START OVER
	701	B=.000001*A(I)														2480
		IF (PRINT.EQ.1) WRITE(6,700)										2490
	700 FORMAT (/ 84H INITIAL CHOICE THRU ORIGIN SEEMINGLY OK, BUT RESULTS  2500
       X ROTTEN AFTER 99 ITERATIONS .../48H SET X0=.000001*A AND USE MORE   2510
	   XREFINED TECHNIQUE /) 												2520
		DR=B/ADIV															2530
		BDR=(A(I)-B)/BDIV													2540
		MMAX=INIT+JINT+1													2550
		MINT=INIT+1															2560
		MAX=MMAX-1															2570
		IM=0																2580
		KJ=0																2590
		GO TO 237 															2600
		END																	2610




























