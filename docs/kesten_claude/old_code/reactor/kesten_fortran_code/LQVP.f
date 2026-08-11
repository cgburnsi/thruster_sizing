		SUBROUTINE LQVP (H,ZLV,Q,JJ,Q1,TEMP)	 							   0
		REAL KP,K															  10
		INTEGER PRINT														  20								  50
		COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  30
	   1SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  		  40
		COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  50
	   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2	  60
	   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 						      70
		COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								  80
		COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			  90
		COMMON /MUVST/VISVST(30)											 100
		COMMON /FLAGS/MFLAG,KFLAG,PRINT 									 110
		WRITE(6,100) 														 120
 100	FORMAT109H1 **** ENTERING L 										 130
	   1IQUID-VAPOR REGION **** //) 										 140
		DERIV(JJ) = Q														 150
		DHDZ(JJ)=Q1 														 160
		WRITE(6,1750) 														 170
   1750 FORMAT(//30X,  3H Z, 11X,5H  TEMP 11X, 3H H 12X, 3HWFM) 			 180
		Z(JJ)=ZLV 															 190
		JJ=JJ+1																 200
   1800 CONTINUE 															 210
   1820 DERIV(JJ)=DERIV(JJ-1) 												 220
		Z(JJ)=Z(JJ-1)+DZ													 230
		TEMP=TVAP															 240
		CALL UNBAR (TBLH4(1),1,TEMP,0.,H4,KK)								 250
		CALL UNBAR (ZTBLAP(1),1,Z(JJ),0.,AP,KK) 							 260
		CALL PARAM(TEMP,Z(JJ),1,0.0,0.0,0,G,GMMA,K,DPA,BETA) 				 270
		DHDZ(JJ)=-(H4*DPA*AP*DERIV(JJ)+FC*(H-HF))/G 						 280
		DZ=-H4/(ENMX2*DHDZ(JJ)) 											 290
		IF(H-HV)82,1850,82 													 300
	82  H=H+DHDZ(JJ)*DZ 													 310
	    IF(H-HV) 1850,1850,2000												 320
   1850 WFV=(H-HL)/(HV-HL)													 330
		WRITE (6,1900) Z(JJ),TEMP,H,WFV 									 340
   1900 FORMAT(22X,E14.5,1X,E14.5,1X,E14.5,1X,E14.5)						 350
		IF(H-HV)83,1950,83													 360
	 83 JJ=JJ+1 															 370
		GO TO 1800															 380
   2000 DZ=(HV-H)/DHDZ(JJ)+DZ												 390
		H+HV																 400
		JJ=HH+1																 410
		GO TO 1820 															 420
   1950 RETURN																 430
		END																	 440
