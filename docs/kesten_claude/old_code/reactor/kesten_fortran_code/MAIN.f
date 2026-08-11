		REAL KP,K 															   0
		INTEGER OPTION,PRINT												  10
		COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  20
	   1        SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  30
		COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  40
	   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2	  50
	   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 							  60
		COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								  70
		COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			  80
		COMMON /MUVST/VISVST(30)											  90
		COMMON /FLAGS/MFLAG,KFLAG,PRINT 									 100
		COMMON /IFCE00/IFC,GATZ0											 110
		COMMON /LIZTBL/DHVST(18),DHLVST(18)									 120
		COMMON /DAVTBL/VPTBL(44)											 130
		DIMENSION  TITLE(14) 												 140
		READ (5,700) NCASE 													 150
	700 FORMAT (I3)															 160
		KOUNT=1																 170
	705 READ (5,608) TITLE 													 180
	608 FORMAT (14A6)														 190
		WRITE (6,609) TITLE													 200
	609 FORMAT (1H1,14A6//)													 210
		IFC=1																 220
		READ (5,809) OPTION,PRINT,NOFZ										 230
	809 FORMAT (2I2,I3)														 240
		READ (5,800) Z0,G0,FC,ALPHA3,HF,R,WM4,WM3,WM2,WM1,ALPHA1,ALPHA2,	 250
	   X      AGM,BGM,KP,CGM,TF,CFL,ENMX1,ENMX2,ENMX3,DIF3,DIF4,PRES,ZEND,   260
	   X      EN1,EN2,EN3,													 270
	800 FORMAT (8E10.5)														 280
		NZTBL = 2*NOFZ+4													 290
		NOFZ4 = NOFZ+4														 300
		NOFZ5 = NOFZ4+1 													 310
		CALL UNBAR (VPTBL(1),1,PRES,0.,TVAP,KK)								 320
		CALL UNBAR (DHVST(1),1,TVAP,0.,DELHV,KK)							 330
		CALL UNBAR (DHLVST(1),1,TVAP,0.,DELHL,KK)							 340
		HL=(TVAP-TF)*CFL													 350
		HV=HL+DELHV-DELHL													 360
		GATZ0=G0+FC*Z0														 370
		IF(FC.GT.0.)GO TO 837												 380
		IFC=0																 390
	637 WRITE (6,600)														 400
	600 FORMAT (52X,16H INPUT CONSTANTS/7X,102H HF      HL        HV  	     410
	   X       TF      TVAP    CFL     PRESSURE    KP    FC                  420
	   X       G0) 															 430
		WRITE (6,601) HF,HL,HV,TF,TVAP,CFL,PRES,KP,FC,G0					 440
	601 FORMAT (3X,10E11.6//)												 450
		WRITE (6,602)														 460
	602 FORMAT (7X,103H  R    ALPHA3   CGM   DIF3    DIF4 					 470
	   X        WM4     WM3     WM2     WM1   ZEND)							 480
		WRITE (6,601) R,ALPHA3,CGM,DIF3,DIF4,WM4,WM3,WM1,WM1,ZEND			 490
		WRITE (6,603)														 500
	603	FORMAT (6X,113H  AGM   BGM   ALPHA1    ALPHA2   N1 					 510
	   X   N2    N3   ENMX1   ENMX2   ENMX3     )							 520
		WRITE (6,601) AGM,BGM,ALPHA1,ALPHA2,EN1,EN2,EN3,ENMX1,ENMX2,ENMX3    530
		WRITE (6,617) Z0													 540
	617 FORMAT (// 8X,'Z0' / 3X,E11.6)										 550
		READ (5,20) (ZTBLA(I),I=1,4)										 560
	 20 FORMAT (4E8.4)														 570
		READ (5,21) (ZTBLA(I),I=5,NOFZ4)    								 580
     21 FORMAT (10E8.4) 													 590
		READ (5,21) (ZTBLA(I),I=NOFZ5,NZTBL)								 600
		READ (5,20) (ZTBLAP(I),I=1,4)										 610
		READ (5,21) (ZTBLAP(I),I=5,NOFZ4)									 620
		READ (5,21) (ZTBLAP(I),I=NOFZ5,NZTBL)								 630
		READ (5,20) (ZTBLD(I),I=1,4)										 640
		READ (5,21) (ZTBLD(I),I=5,NOFZ4)									 650
		READ (5,21) (ZTBLD(I),I=NOFZ5,NZTBL)								 660
		WRITE (6,604)														 670
	604	FORMAT (///55X,13H Z VS A TABLE)									 680
		WRITE (6,22) (ZTBLA(I),I=1,4										 690
	 22 FORMAT (40X,4E13.5)												     700
		WRITE (6,23) (ZTBLA(I),I=5,NOFZ4)									 710
	 23 FORMAT (1X,10E13.5)													 720
		WRITE (6,25)														 730
	 25 FORMAT ( / ) 														 740
		WRITE (6,23) (ZTBLA(I),I=NOFZ5,NZTBL)								 750
		WRITE (6,24) 														 760
	 24 FORMAT ( // )														 770
		WRITE (6,606)														 780
	606 FORMAT (54X,14H Z VS AP TABLE)										 790
		WRITE (6,22) (ZTBLAP(I),I=1,4										 800
		WRITE (6,23) (ZTBLAP(I),I=5,NOFZ4)									 810
		WRITE (6,25)														 820
		WRITE (6,23) (ZTBLAP(I),I=NOFZ5,NZTBL)								 830
		WRITE (6,24) 														 840
		WRITE (6,607)														 850
	607 FORMAT (54X,14H Z VS DELTA TABLE)									 860
		WRITE (6,22) (ZTBLD(I),I=1,4										 870
		WRITE (6,23) (ZTBLD(I),I=5,NOFZ4)									 880
		WRITE (6,25)														 890
		WRITE (6,23) (ZTBLD(I),I=NOFZ5,NZTBL)								 900
		WRITE (6,613) 														 910
	613 FORMAT (18X, ******************************** ENTERING LIQUID        920
	   X REGION     ********************************) 						 930
		MFLAG=0																 940
		DZ=0.0																 950
		Z(1)=0.0															 960
		H=HF																 970
		II=2																 980
	850 Z(II)=Z(II-1)+DZ 													 990
		TEMP=TF+(H-HF)/CFL													1000
		CALL UNBAR (TBLVP(I),1,TEMP,0.,VP,KK)								1010
		CN2H4=(VP*WM4)/(R*TEMP)												1020
		CALL UNBAR (TBLH4(I),1,TEMP,0.,H4,KK)								1030
		CALL UNBAR (ZTBLAP(I),1,Z(II),0.,AP,KK)								1040
		CALL UNBAR (ZTBLA(I),1,Z(II),0.,A,KK)								1050
		CALL PARAM(TEMP,Z(II),1,CN2H4,H4,0,G,GMMA,K,DPA,BETA)				1060
		CALL SLOPE (CN2H4,GMMA,K,BETA,EN1,DERIV(II),DPA,A,DIF4)			    1070
		IF(H-HL)777,776,777													1080
	776 IF(MI.GT.20)DERIV(II)=DERIV(II-1)									1090
	777 DHDZ(II)=-(H4*DPA*AP*DERIV(II)+FC*(H-HF))/G							1100
		DZ=-H4/(ENMX1*DHDZ(II))												1110
		WRITE(6,820)														1120
	820 FORMAT (/39X,48H  Z    TEMP    H DHDZ)								1130
		WRITE(6,860) Z(II),TEMP,H,DHDZ(II)									1140
	860 FORMAT (/30X,4E15.6)												1150
		IF(H-HL) 874,1020,874												1160
	874 H=H+DHDZ(II)*DZ														1170
		IF(H-HL) 875,1020,1000												1180
	875 II=II+1																1190
		GO TO 850															1200
C		BACKSTEP TO L-L-V-BOUNDARY 											
   1000 DZ=(HL-H)/DHDZ(II)+DZ												1210
		H=HL																1220
		II=II+1																1230
		GO TO 850															1240
   1020 IF(OPTION.EQ.2) CALL LQV2(H,Z(II),DERIV(II),II,DHDZ(II),TEMP,CN2H4) 1250
		IF(OPTION.EQ.2) GO TO 1021											1260
		CALL LQVP(H,Z(II),DERIV(II),II,DHDZ(II),TEMP)						1270
C		START VAPOR REGION
   1021 DZ=-H4/(ENMX2*DHDZ(II))												1280
		CALL VAPOR(TEMP,Z(II),II,DHDZ(II),DERIV(II),H)						1290
		KOUNT=KOUNT+1														1300
		IF(KOUNT.LE.NCASE) GO TO 705 										1310
		WRITE(6,102)														1320
	102 FORMAT (////41X,36H *****   OPERATIONS COMPLETE *****)				1330
		STOP																1340
		END																	1350

		












