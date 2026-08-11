		SUBROUTINE PARAM (T,ZA,LOP,CC,HR,LVOP,G,GMMA,K,DPA,BETA)				   0
		REAL KP,K 															  10
		INTEGER PRINT	         											  20
		COMMON /FTZ/TBLVP(70),TBLH4(42),TBLH3(42),SHTBL1(34),SHTBL2(34), 	  30
	   1        SHTBL3(34),SHTBL4(34),ZTBLD(46),ZTBLAP(46),ZTBLA(46) 		  40
		COMMON /CO/HL,HV,FC,TF,CFL,CGM,ENMX1,AGM,DIF3,DIF4,KP,PRES,G0,		  50
	   1        WM4,WM3,WM2,WM1,ALPHA3,R,TVAP,ZEND,BGM,HF,DZ,ALPHA1,ALPHA2	  60
	   2        ,ENMX2,ENMX3,EN1,EN2,EN3,H,RAT,MI 							  70
		COMMON /VAR/DERIV(250),DHDZ(250),Z(250)								  80
		COMMON /TOLL/ALIM,OPTION,C1,C2,C3,C4,CAV,G,TEMP,AP,WMAV,Z0,			  90
		COMMON /MUVST/VISVST(30)											 100
		COMMON /FLAGS/MFLAG,KFLAG,PRINT 									 110
		COMMON /IFCE00/IFC,GATZ0											 120
		IF(ZA-Z0)43,48,48													 130
	 43 G=G0+FC*ZA 															 140
		GO TO 52															 150
C		Z HAS EXCEEDED HYDRAZINE INJECTOR TUBE LENGTH ??? CONSTANT FROM
C       HERE TO END OF BED
	 48 G=GATZ0																 160
		FC=0. 																 170
	 52 IF(LVOP.EQ.1)GO TO 1004 											 180
		GMMA=AGM/T															 190
C		CALCULATE K,DPA FOR N2H4
		K=ALPHA1*EXP(-GMMA)													 200
   1001 DPL=DIF4*(T/492.)**1.832 *(14.7/PRES)*(1.-EXP(-.0672*(PRES*492.)/(   210
	   114.7*T))) 															 220
		DPA=DPL																 230
   1002 BETA=-(CC+HR*DPA)/(KP*T)											 240   (This really needs to be checked.  It might be CC+HR not CC*HR
		GO TO 1003 															 250
   1004 GMMA=BGM/T															 260
C		CALCULATE K,DPA FOR NH3										
		K=ALPHA2*EXP(-GMMA)	/C1**1.6    									 270
   1001 DPV=DIF3*(T/492.)**1.832 *(14.7/PRES)*(1.-EXP(-.0672*(PRES*492.)/(   280
	   114.7*T))) 															 290
		DPA=DPV																 300
		GO TO 1002															 310
   1003 RETURN																 320
		END																	 330
