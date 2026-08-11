	SUBROUTINE SGRAD (GRAD,TGRAD)										       0
	REAL K0,KP,KC3,KC4,MU												  	  10
	INTEGER PRINT														  	  20
	COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  	  30
   1        SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  	  40
	COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  	  50
   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2		  60
   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 							      70
	COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								      80
	COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			  	  90
	COMMON /BBBB/DP3,A,KC3,K0,X0A,CPS,CI3,GAMMA,BETA						 100	
	COMMON /DDD/CFTBL4(34),CFTBL3(34),CFTBL2(34),CFTBL1(34)					 110
	COMMON /MUVST/VISVST(30)											  	 120
	COMMON /FLAGS/MFLAG,KFLAG,PRINT 									     130
	COMMON /CCC/H4TBL(40),H3TBL(40)											 140
	DIMENSION  CPOX(101),PCPOX(101),DX(101),CPX(101),RHET(101) 				 150
C	DEFINE DP FUNCTION
	DP3F(X,Y,Z) = 14.7*Y/Z*(X/492.)**1.823*(1.-EXP(-.0672*Z*492./(			 160
   X             14.7*X)))													 170
C   DEFINE KC FUNCTION
    KCF(A,B,C,D,E) = .61*A/B*(C/(B*D))**-.667*(A/(E*C))**-.41				 180
C   DEFINE ANALYTIC INTEGRATION FUNCTIONS FROM INTEGRAL EQUATION
    EVAL1(A,B) = B**3/3.-A**3/3.											 190
    EVAL2(A,B) = B**2/2.-A**2/2.											 200
	WAF1 = .8																 210
	WAF2 = .2 																 220
  1 LTFLG=0																	 230
	P=PRES																	 240
	T=TEMP																	 250
	ALPH2=ALPHA2															 260
	CI1=C1																	 270
	CI2=C2																	 280
	CI3=C3																	 290
	CI4=C4																	 300
	D03=DIF3																 310
	D04=DIF4																 320
	N=EN1 																	 330   (I think this is correct.  It needs checked for what LN1 is)
	NPART = 50																 340
	LP1  = 1																 350
	TPSP = 0.																 360
	RHO = CI1+CI2+CI3+CI4													 370
	DI3 = D03*14.7/P*(T/492.)**1.823 										 380
	DI4 = D04*14.7/P*(T/492.)**1.823 										 390
	CALL UNBAR (VISVST,1,T,0.,MU,KK) 										 400
	CALL UNBAR (CFTBL4,1,T,0.,CF4,KK) 										 410
	CALL UNBAR (CFTBL3,1,T,0.,CF3,KK) 										 420
	CALL UNBAR (CFTBL2,1,T,0.,CF2,KK) 										 430
	CALL UNBAR (CFTBL1,1,T,0.,CF1,KK) 										 440
	KC3 = KCF(G,RHO,MU,DI3,AP) 												 450
	KC4 = KCF(G,RHO,MU,DI4,AP) 												 460
	CFBAR = (CI1*CF1+CI2*CF2+CI3*CF3+CI4*CF4)/RHO							 470
	HC = .74*G*CFBAR*(G/(AP*MU))**-.41 										 480
C
C	LOCATE SUITABLE X0
C
	DP3 = DP3F(T,D03,P) 													 490
C	CHOOSE STARTING VALUE FOR CPS TO BE = CI3/2.
	CPS = CI3/2. 															 500
	CMCPN = CI3-CPS 														 510
	DCPDX = KC3/DP3*(CI3-CPS)												 520
C	H4 CONSTANT FOR EACH ENTRY TO THIS ROUTINE
C	H3 VARIES WITH TEMP AT HEAT ITERATION
	CALL UNBAR (H4TBL,1,T,0.,H4,KK)											 530
	CALL UNBAR (H3TBL,1,T,0.,H3,KK)											 540
	IF (LP1.EQ.1) GO TO 6 													 550
 40 TPSPP = TPSP 															 560
	TPSP = TPS 																 570
  6	TPS = T-(H4*KC4*CI4+H3*DP3*DCPDX)/HC									 580
	IF (TPS.LT.0) TPS=1.													 590
	CALL UNBAR (H3TBL,1,TPS,0.,H3,KK)										 600
	DP3 = DP3F(TPS,D03,P)													 610
	DP3P = DP3 															     620
	H3P = H3																 630
	TMTPN = T-TPS 															 640
 61 GAMMA = BGM/TPS 														 650
	BETA = -CPS*H3*DP3/(KP*TPS)												 660
	KO = ALPH2*EXP(-GAMMA)*CI1**EN3											 670
C   LINEAR EXTRAPOLATION USED TO GUESS AT X0								 680
	X0P = X0 																 690
	X0 = A-CPS/DCPDX														 700
	IF (X0) 11,12,12														 710
 11 X0 = 0.																	 720
	X0A = 0. 																 730
	CPS = CI3/(DP3/(A*KC3)+1.)												 740
	DCPDX = CI3/A															 750
	TPS = T-(H4*KC4*CI4+H3*DP3*DCPDX)/HC									 760
	IF (TPS.LT.0.) TPS=1. 													 770
	WRITE (6,132) 															 780
132 FORMAT (// 37X,'WE HAVE CALCULATED A NEGATIVE X0 DURING ITERATION')		 790
C	INTEGRATE FOR CP EQUATION
 12 CALL TRAPP (X0A,1.,NPART,RIESUM)										 830
C	CALCULATE NEW CPS ...
	CPSP = CPS 																 840
	CMCPO = CMCPN 															 850
	CPS = CI3-A*RIESUM/KC3													 860
	IF (LTFLG-1) 80,84,80 													 870
 80 IF (CPS-(.25*CI3)) 81,81,130											 880
 81 LTFLG=1 																 890
	GO TO 22																 900
 84 LTFLG=0																	 910
	IF (CPS) 89,130,130														 920
 89 CPS=0.																	 930
	GO TO 46 																 940
130 CMCPN = CI3-CPS 														 950
C   CALCULATE NEW TP										
 13 DCPDX = KC3/DP3*(CI3-CPS) 												 960 (I think this is line number 13)
	GRAD = DCPDX*DP3 														 970
	TGRAD = HC*(T-TPS) 														 980
	TPSPP = TPSP 															 990
	TPSP = TPS 																1000
	TMTPO = TMTPN 															1010
 51 TPS = T-(H4*KC4*CI4+H3*DP3*DCPDX)/HC									1020
	IF (TPS.LT.0.) TPS=1. 													1030
	CALL UNBAR (H3TBL,1,TPS,0.,H3,KK)										1040
	DP3 = DP3F(TPS,D03,P) 													1050
	TMTPN = T-TPS 															1060
	GAMMA = BGM/TPS 														1070
	BETA = -CPS*H3*DP3/(KP*TPS)												1080
	KO = ALPH3*EXP(-GAMMA)*CI1**EN3											1090
C   TEST TEMP, CONCENTRATION FOR 5% LIMIT
	IF (ABS(TMTPO-TMTPN)/TMTPN - .05) 41,41,43								1100
 41 IF (ABS(CMCPO-CMCPN)/CMCPN - .05) 70,70,43								1110
C   TEST FOR TEMPERATURE LOOP .. COMPARE LAST 3 TEMPS
 43 IF (AMIN1(TPS,TPSP,TPSPP) - TPSP) 60,71,60 								1120
 60 IF (AMAX1(TPS,TPSP,TPSPP) - TPSP) 46,71,46 								1130
C	TEMPERATURE HAS FLUCTUATED .. TASKE AVERAGE AND RECALCULATE CPS
 71 TPSPP = TPSP 															1140
	TPSP = TPS 																1150
	TMTPO = TMTPN 															1160
	TPS = (TPSP+TPSPP)/2.													1170
	CALL UNBAR (H3TBL,1.TPS,0.,H3,KK) 										1180
	DP3 = DP3F(TPS,D03,P)													1190
	DP3P = DP3 																1200
	TMTPN = T-TPS															1210
	DCPDX = (HC*(T-TPS)-H4*KC4*CI4)/(H3*DP3)								1220
	CPSP = CPS																1230
	CMCPO = CMCPN 															1240
	CPS = CI3-DP3/KC3*DCPDX													1250
	IF (CPS.LT.0.) CPS=0. 													1260
	CMCPN = CI3-CPS 														1270
	LP1 = LP1 + 1															1280
	IF (LP1-50) 61,61,44													1290
C	NO CONVERGENCE YET ... AVERAGE THE CPS'S FOR LAST TWO CALC'S AND REPEAT
 46 CPS = .2*CPS+.8*CPSP 													1300
	GO TO 53 																1310
 22 X00 = WAF1*X0P+WAF2*X0 													1320
	CPS = CI3/(1.+DP3/(KC3*A-KC3*X00)) 										1330
 53	DCPDX = KC3/DP3P*(CI3-CPS)												1340  (I think this is line number 53.  it might be 58)
	CMCPN = CI3-CPS 														1350
	H3 = H3P 																1360
 42 LP1 = LP1+1 															1370
	IF (LP1-25) 40,40,44 													1380
 44 WAF1 = WAF1+.05															1390
	IF (WAF1.GT.0.95) GO TO 99 												1400
	WAF2 = 1.-WAF1 															1410
C	NO CONVERGENCE WITH PRESENT WEIGHTED AVERAGE FACTORS FOR X0
C	REPEAT ITERATION PROCEDURE WITH NEW FACTORS
	GO TO 1																	1420
 99 WRITE (6,98) 															1430
 98 FORMAT (///20X,'UNABLE TO FIND SUITABLE X0 AFTER FOUR TRIES OF 25		1440
   X		57X,5H X0 -,E12.5 )												1450
	CALL EXIT																1460
C 	SATISFACTORY X0 HAS BEEN FOUND
 70 IF (PRINT.EQ.1) WRITE(6,16)LP1,X0										1470
 16 FORMAT (/// 46X,27HSATISFACTORY X0 FOUND AFTER,I3,7H TRIES /			1480
   X 		57X,5H X0=,E12.5 )												1490









C
C	CALCULATE GRADIENT
C
131 LP2 = 1																	1500
	NX = 24 																1510
	NX1 = NX+1																1520
	NXM1 = NX-1																1530
291 X0A = X0/A																1540
	VNU = -KC3/DP3 															1550 
	INT1 = 1 																1560
	R = 2																	1570
	R1 = 0. 																1580
	R2 = 0. 																1590
	PS1 = 0. 																1600
	PS2 = 0. 																1610
	DELX0A = (1.-X0A)/FLOAT(NX)												1620
C 	CALCULATE PROFILE CURVES FOR INTEGRAND FUNCTIONS
	XA = X0A 																1630
	DO 770 I=1,NX1 															1640
C	CP(X/A) IS A LINEAR PROFILE DURING FIRST APPROXIMATION
	IF (LP2,GT.1) GO TO 664 												1650
	CPX(I) = (XA-X0A)/(1.-X0A)*CPS											1660
664 RHET(I) = K0*CI3**(1-N)*CPX(I)**N*EXP(GAMMA*BETA*(1.-CPX(I)/CI3)/ 		1670
   X          (1.+BETA*(1.-CPX(I)/CI3)) 									1680
	DX(I) = XA 																1690  (I think this is DX(I), but it might be DX(1))
	XA = XA+DELX0A															1700
770 CONTINUE																1710
C	TAKE INTERVAL FUNCT'N MIDPTS AS CONSTANT VALUE FOR CP(X/A) AND RHET 	
	DO 771 I=1,NX															1720
	CPX(I) = (CPX(I)+CPX(I+1))/2.											1730
	RHET(I) = (RHET(I)+RHET(I+1))/2.										1740
771	CONTINUE 																1750
	XA = X0A+DELX0A															1760
	CTRM = (A*VNU+1.)/(A*VNU)												1770
C	INTEGRAL EQUATION FOLLOWS	
C	CPOX(1) IS SPECIAL CASE ... X=X0												(I think this comment says CPOX(1), but it also might be CPOX(I).)
	DXL = X0A 																1780
	DXU = DXL+DELX0A 														1790
	RR1 = 0.																1800
	DO 377 I=1,NX															1810
	RR1 = RR1+RHET(I)*(EVAL2(DXL,DXU)-CTRM*EVAL1(DXL,DXU))					1820
	DXL = DXU 																1830
	DXU = DXU+DELX0A														1840
377 CONTINUE																1850
	CPOX(1) = CI3-A*A/DP3*RR1												1860 (Damn.  I still can't tell if this is CPOX(1) or CPOX(I).)
	IF (CPOX(1).LT.0.)  CPOX(1) = 0. 										1870 (Arg. It's still not clear if it is a 1 or an I).
C	SOLVE GENERAL EQUATION OF TWO INTEGRALS FOR CP(X/A)
769 DO 772 I=1,INT1															1880
	R1 = R1+RHET(I)*EVAL1(X0A,XA)											1890
	X0A = XA 																1900
	XA = XA+DELX0A 															1910
772 CONTINUE																1920
	R1 = R1*(1./X0A-CTRM)													1930
	XAD = XA																1940
	XA = XA-DELX0A															1950
	DO 773 I=INT1,NXM1														1960
	PS1 = PS1+RHET(I+1)*EVAL2(XA,XAD)										1970
	PS2 = PS2+RHET(I+1)*EVAL1(XA,XAD)										1980
	XA = XAD 																1990
	XAD = XAD+DELX0A														2000
773 CONTINUE																2010
	R2 = PS1-CTRM*PS2														2020
	INT1 = INT1+1															2030
	CPOX(K) = CI3-A*A/DP3*(R1+R2)											2040
	IF (CPOX(K).LT.0.) CPOX(K)=0.											2050
	X0A = X0/A																2060
	XA = X0A+DELX0A															2070
	K = K+1																	2080
	R1 = 0. 																2090
	R2 = 0. 																2100
	PS1 = 0. 																2110
	PS2 = 0. 																2120
	IF (K.LE.NX) GO TO 769 													2130
C	CPOX(NX1) IS SPECIAL CASE ... X=A 										
	DXL = X0A 																2140
	DXU = DXL+DELX0A 														2150
	RR2 = 0. 																2160
	DO 378 I=1,NX															2170
	RR2 = RR2 + RHET(I)*EVAL1(DXL,DXU)										2180
	DXL = DXU 																2190
	DXU = DXU+DELX0A 														2200
378 CONTINUE 																2210
	CPOX(NX1) = CI3-A*A/DP3*(1.-CTRM)*RR2 									2220
	IF (CPOX(NX1).LT.0.) CPOX(NX1)=0. 										2230
C	CALCULATE A NEW TPS
	DCPDX = KC3/DP3*(CI3-CPOX(NX1))											2240
	H3P = H3 																2250
	DP3P = DP3																2260
	TPS = T-(H4*KC4*CI4+H3*DP3*DCPDX)/HC 									2270
	CALL UNBAR (H3TBL,1,TPS,0.,H3,KK) 										2280
	DP3 = DP3F(TPS,D03,P) 													2290
	TMTPO = TMTPN 															2300
    TMTPN = T-TPS															2310
C	TWO PASSES NEEDED BEFORE CHECK ON TEMP,CONC CAN BE MADE
 33 IF (LP2.EQ.1) GO TO 27 													2320
	CMCPO = CMCPN 															2330
	CMCPN = CI3-CPOX(NX1) 													2340
	IF (ABS(TMTPO-TMTPN)/TMTPN - .05) 26,26,27 								2350
 26 IF (ABS(CMCPO-CMCPN)/CMCPN - .05) 88,88,27 								2360
C
C	CALCULATE NEW CPX(I) PROFILE FOR NEXT PASS
C
 27 DO 55 I=1,NX1															2370
	IF (MOD(LP2,5)) 34,57,34 												2380   (I think this is MOD(LP2,5))
C	CALCULATE WEIGHTED AVERAGE OF OLD AVERAGED AND CALCULATED PROFILES 		
 34 CPX(I) = .8*CPX(I)+.2*CPOX(I)											2390
	GO TO 56																2400
C	AVERAGE PRESENT AND PAST CALCULATED PROFILES EVERY 5TH PASS TO SMOOTH   
 57 CPX(I) = (CPOX(I)+PCPOX(I))/2. 											2410  (I think this is CPOX(I)+PCPOX(I).  It's hard to read in the scanned document.
C	STORE PRESENT CALCULATED PROFILE
 56 PCPOX(I) = CPOX(I) 														2420
 55 CONTINUE 																2430
	CMCPN = CI3-CPX(NX1) 													2440
	DCPDX = KC3/DP3P*(CI3-CPX(NX1))											2450
	TPS = T-(H4*KC4*CI4+H3P*DP3P*DCPDX)/HC 									2460
	IF (TPS.LT.0.) TPS=1. 													2470
	CALL UNBAR (H3TBL,1,TPS,0.,H3,KK)										2480
	DP3 = DP3F(TPS,D03,P) 													2490
	TMTPO = TMTPN 															2500
	TMTPN = T-TPS 															2510
	LP2 = LP2+1																2520
	IF (LP2-50) 29,29,30 													2530
 30 WRITE (6,18) CPOX(NX1) 													2540
 18 FORMAT (///31X,52HUNABLE TO CONVERGE ON CPS IN 50 TRIES ... CP(X/       2550
   XA) =,E12.5) 															2560
	WRITE (6,522) GRAD,TGRAD												2570
522 FORMAT (51X,'KC3*(CI3-CPS) =',E12.5 / 54X,'HC*(T-TPS) =',E12.5) 		2580
	GO TO 28																2590
 29 GAMMA = BGM/TPS 														2600
	BETA = -CPX(NX1)*H3*DP3/(KP*TPS) 										2610
	K0 = ALPH2*EXP(-GAMMA)*CI1**EN3											2620
	GO TO 291 																2630
 88 IF (PRINT.EQ.1) WRITE(6,19)(DX(I),CPOX(I),I=1,NX1)						2640
 19 FORMAT (///11X,114H/A   CP(X/A)   X/A  CP(X/A)							2650
   X X/A   CP(X/A)   X/A   CP(X/A)   X/A   CP(X/A) 							2660
   X     /(5X,10E12.6)) 													2670
	IF (PRINT.EQ.1) WRITE(6,82)LP2,CPOX(NX1)								2680
 82 FORMAT (///41X,34HCONCENTRATION GRADIENT FOUND AFTER,I3,6H TRIES / 		2690
   X      45X,27HCP(X) AT PARTICLE SURFACE =,E12.5) 						2700
	GRAD = DCPDX*DP3 														2710
	TGRAD = HC*(T-TPS) 														2720
	IF (PRINT.EQ.1) WRITE(6,83)GRAD,TGRAD 									2730
 83 FORMAT (51X,'KC3*CI3-CPS) =',E12.5 / 54X,'HC*T-TPS) =',E12.5) 			2740
 28 RETURN 																	2750
	END																		2760

























