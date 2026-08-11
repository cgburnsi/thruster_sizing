	   SUBROUTINE UNBAR(T,IK,XIN,YIN,ZZ,KK)                               0
      DIMENSION T(1),X(6),Y(6),A(6)                                      10
C                                                                  UNBAR003
C     ----------------- MARCH 4, 1961 ----------------             UNBAR004
C     ----------------- MODIFIED 7/62 ----------------             UNBAR005
C     --------TO DO QUADRATIC AND LINEAR INTERPOLATION ALSO        UNBAR006
C                                                                  UNBAR007
      II = IK+1                                                          20
      N = 3                                                              30
      N2 = 2                                                             40
         IF (T(II)-3.) 700,701,702                                       50
  700    IF (T(II)+0.) 60,701,704                                        60
  704    IF (T(II)-2.) 705,706,701                                       70
  705    N = 1                                                           80
         GO TO 707                                                       90
  706    N = 2                                                           100
  707    N2 = 1                                                          110
  701    II = II+1                                                       120
  702    N1 = N+1                                                        130
         DO 50 L = II,II                                                 140
         IF ( T(L) + 0. ) 60,60,51                                       150
   60    KK = -1                                                         160
         ZZ = 0.                                                         170
         GO TO 9999                                                      180
   51    NX = T(L)                                                       190
         IF (T(L+1) + 0. ) 60,52,50                                      200
   52    NY = 0                                                          210
         GO TO 53                                                        220
   50    NY = T(L+1)                                                     230
   53    CONTINUE                                                        240
         KK  = 0                                                         250
         KY  = 0                                                         260
         XX  = XIN                                                       270
         YY  = YIN                                                       280
         J1  = II+2                                                      290
         J2  = NX+II+1                                                   300
         IF(XX-T(J1))301,306,400                                         310
  400    DO 302 J=J1,J2                                                  320
         IF (XX-T(J)) 304,304,302                                        330
  302    CONTINUE                                                        340
  309    KK = 2                                                          350
		 XX = T(J2)                                                      360
  308    JX1 = J2-N                                                      370
         GO TO 305                                                       380
  301    KK = 1                                                          390
         XX = T(J1)                                                      400
  306    JX1 = J1                                                        410
         GO TO 305                                                       420
  304    IF (J-J1-1) 301,306,307                                         430
  307    IF (J-J2)   303,308,309                                         440
  303    JX1 = J-N2                                                      450
  305    CONTINUE                                                        460
         XINT = XX                                                       470
         IF (NY) 1500,1500,3000                                          480
 1500    DO 1599 L=1,N1                                                  490
         X(L) = T(JX1)                                                   500
         LY = JX1 + NX                                                   510
         Y(L) = T(LY)                                                    520
 1599    JX1 = JX1+1                                                     530
         I = 1                                                           540
         GO TO 54                                                        550
 3000    J1 = J1+NX                                                      560
         J2 = J2+NY                                                      570
         IF(YY-T(J1))311,316,401                                         580
  401    DO 312 J=J1,J2                                                  590
         IF (YY-T(J)) 314,314,312                                        600
  312    CONTINUE                                                        610
  319    KY = 6                                                          620
         YY = T(J2)                                                      630
  318    JY1 = J2-N                                                      640
         GO TO 315                                                       650
  311    KY = 3                                                          660
         YY = T(J1)                                                      670
  316    JY1 = J1                                                        680
         GO TO 315                                                       690
  314    IF (J-J1-1) 311,316,317                                         700
  317    IF (J-J2)   313,318,319                                         710
  313    JY1 = J-N2                                                      720
  315    CONTINUE                                                        730
         JX2  = JX1                                                      740
         LY   = JY1 + NY*(JX2-II-1)                                      750
         LY1  = LY                                                       760
         DO 3099 L=1,N1                                                  770
         X(L) = T(JX2)                                                   780
         Y(L) = T(LY1)                                                   790
         LY1  = LY1+NY                                                   800
 3099    JX2  = JX2+1                                                    810
         I    = 0                                                        820
         GO TO 54                                                        830
 3090    T(1) = ZZ                                                       840
         DO 4400 I=1,N                                                   850
         LY1  = LY+I                                                     860
         Y(I+1) = 0.                                                     870
         DO 4050 MM=1,N1                                                 880
         Y(I+1) = Y(I+1) + T(LY1)*X(MM)                                  890
 4050    LY1  = LY1+NY                                                   900
 4400    CONTINUE                                                        910
         DO 4199 L=1,N1                                                  920
         X(L) = T(JY1)                                                   930
 4199    JY1  = JY1+1                                                    940
         XINI = YY                                                       950
         I = 1                                                           960
   54    D = 1.                                                          970
         X(N+2) = X(1)                                                   980
         X(N+3) = X(2)                                                   990
         DO 55 J=1,N1                                                   1000
         A(J+1) = X(J+1) - X(J)                                         1010
         TPAL1  = XINT - X(J)                                           1020
         IF ( TPAL1 ) 57,58,57                                          1030
   58    ZZ = Y (J)                                                     1040
         X(1) = 0.                                                      1050
         X(2) = 0.                                                      1060
         X(3) = 0.                                                      1070
         X(4) = 0.                                                      1080
         X(J) = 1.0                                                     1090
         GO TO 59                                                       1100
   57    D = D * TPAL1                                                  1110
         GO TO (711,712,713) ,N                                         1120
  711    X(J) = TPAL1/A(J+1)                                            1130
         GO TO 55                                                       1140
  712    X(J) = -TPAL1                                                  1150
         GO TO 55                                                       1160
  713    X(J) = (X(J+2)-X(J))*TPAL1                                     1170
   55    CONTINUE                                                       1180
         A(1) = A(N+2)                                                  1190
         ZZ = 0.                                                        1200
         DO 56 J=1,N1                                                   1210
         X(J) = D/(A(J)*A(J+1)* X(J))                                   1220
         ZZ = ZZ + Y(J)* X(J)                                           1230
   56    CONTINUE                                                       1240
   59    IF (I) 3098,3098,9999                                          1250
 9999    KK = KK+KY                                                     1260
         RETURN                                                         1270



