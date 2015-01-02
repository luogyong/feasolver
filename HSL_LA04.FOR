C *******************************************************************
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 25 March 1999
C Modified 24 May 2001 to make LA04 threadsafe (LA05 -> LA15)
C 26 June 2001. Arrays W and Y of LA04E/ED made into arguments to
C     conform with Fortran 77. Workspace WS has to be enlarged by M.
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 1st March 2005 Version 1.1.0. FD05 and ZA02 dependence changed
C     to FD15 and ZA12.

      SUBROUTINE LA04AD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,X,
     +                  Z,G,RINFO,WS,LWS,IWS,LIWS)
      INTEGER LA,KB,LB,M,N,LWS,LIWS
      DOUBLE PRECISION
     +                A(LA),B(M),BND(2,KB),C(N),CNTL(15),G(N),
     +                RINFO(40),WS(LWS),X(N+M),Z(N)
      INTEGER IP(N+1),IRN(LA),IWS(LIWS),IX(M),JOB,JX(KB)
      DOUBLE PRECISION DUM
      INTEGER I,IB,IISD,II15,IK15,IC15,IPB,IR,IROW,IWB,J,JB,JIN,
     +        JV,JY,JYY,K,K1,K2,L,LP,MP,NPM,NZ
      DOUBLE PRECISION
     +                ONE,SIGMA,SIGN,T,ZERO
      PARAMETER (ONE=1,ZERO=0)
      DOUBLE PRECISION FD15AD
      DOUBLE PRECISION ZA12AD
      EXTERNAL FD15AD,ZA12AD
      EXTERNAL LA04BD,LA04CD,LA04DD,LA04SD,LA15BD
      INTRINSIC ABS,MAX
      RINFO(9) = ZA12AD(DUM)
      LP = CNTL(6)
      MP = CNTL(7)
      SIGMA = CNTL(1)
      IF (JOB.GT.7 .OR. JOB.LT.0) THEN
        IF (LP.GE.0) WRITE (LP,FMT='(/,/,A,I5)')
     +      ' LA04 error return: JOB has value',JOB
        JOB = -9
        RETURN
      END IF
      IB = MIN(LWS-4*M-4, (LIWS-10*M-11)/2)
      JV = 1 + M
      JY = JV + M
      JYY = JY + M
      JB = JYY + M
      II15 = LIWS - 11
      IK15 = II15 + 3
      IISD = LIWS
      IWB = 2*IB + 1
      IPB = IWB + 8*M
      IR = IWB + 4*M
      IROW = IR + M
      IC15 = LWS - 3
      IF (JOB.GE.6) THEN
        IF (M.GT.0) THEN
          T = ZA12AD(DUM)
          CALL LA15BD(WS(JB),IWS,IB,M,IWS(IPB),IWS(IWB),WS(JV),RINFO(3),
     +                WS,JOB.EQ.7,IWS(II15+1),IWS(IK15+1))
          RINFO(6) = RINFO(6) + ZA12AD(DUM) - T
        END IF
        GO TO 150
      END IF
      IF (JOB.GT.2) GO TO 140
      RINFO(4) = ZERO
      RINFO(5) = ZERO
      RINFO(6) = ZERO
      RINFO(7) = ZERO
      RINFO(29) = 0
      RINFO(30) = 0
      RINFO(33) = 0
      RINFO(27) = M
      CNTL(2) = MAX(CNTL(2), FD15AD('E')*100.0)
      CNTL(11) = MAX(CNTL(11), FD15AD('E')*100.0)
      NZ = IP(N+1) - 1
      IF (MP.GE.0) THEN
        WRITE (MP,FMT='(/,3(A,I7)/6(A,I7))') ' Entry to LA04 with M=',M,
     +    ' N=',N,' No. of non-zeros=',NZ,' LA=',LA,' KB=',KB,' LB=',LB,
     +    ' JOB=',JOB,' LWS=',LWS,' LIWS=',LIWS
        WRITE (MP,'(A,1P,3E11.3,0P,4F5.0,1P,/6X,4E11.3)') ' CNTL=',
     +    (CNTL(K),K=1,11)
      END IF
      IF (JOB.EQ.0) THEN
        L = 4*M + 4*N + 7 + 3
        IF (LWS.LT.L) THEN
          JOB = -13
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: LWS must be at least',L
          RINFO(35) = L
          RETURN
        END IF
        IF (LIWS.LT.LA+10) THEN
          JOB = -14
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: LIWS must be at least',LA+10
          RINFO(35) = LA
          RETURN
        END IF
        DO 10 I = 1,M
          IF (B(I).NE.ZERO) NZ = NZ + 1
   10   CONTINUE
        DO 20 J = 1,N
          IF (C(J).NE.ZERO) NZ = NZ + 1
   20   CONTINUE
      ELSE
        IF (IB.LT.M) THEN
          JOB = -4
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: IB must be at least',M
          RINFO(35) = M
          RETURN
        END IF
      END IF
      IF (LA.LT.NZ) THEN
        JOB = -10
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +      ' LA04 error return: LA must be at least',NZ
        RINFO(35) = NZ
        RETURN
      END IF
      IF (M.LT.0 .OR. KB.LT.0 .OR. LB.LE.KB .OR. N+1.LT.LB) THEN
        JOB = -3
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,/5(A,I7))')
     +      ' LA04 error return: M, N, KB, or LB is faulty',' M=',M,
     +      ' N=',N,' KB=',KB,' LB=',LB
        RETURN
      END IF
      DO 30 J = 1,KB
        IF (BND(1,J).LT.-SIGMA .OR. BND(1,J).GT.BND(2,J) .OR.
     +      BND(2,J).GT.SIGMA) THEN
          JOB = -5
          RINFO(35) = J
          IF (LP.GE.0) WRITE (LP,FMT='(/A/A,I6,A,1P,2E20.12)')
     +        ' LA04 error return: BND is faulty',' BND(:,',J,') =',
     +        BND(1,J),BND(2,J)
          RETURN
        END IF
   30 CONTINUE
      DO 40 I = 1,M
        IWS(I) = 0
   40 CONTINUE
      DO 60 J = 1,N
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2+1) THEN
          JOB = -6
          RINFO(35) = J + 1
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,/2(A,I6,A,I6))')
     +        ' LA04 error return: IP is faulty','  IP(',J,' )=',K1,
     +        '  IP(',J + 1,' )=',IP(J+1)
          RETURN
        END IF
        DO 50 K = K1,K2
          I = IRN(K)
          IF (I.LE.0 .OR. I.GT.M) THEN
            JOB = -7
            RINFO(35) = K
            IF (LP.GE.0) WRITE (LP,FMT='(/,A,/,A,I6,A,I6)')
     +          ' LA04 error return: IRN is faulty',' IRN(',K,' )=',I
            RETURN
          ELSE IF (IWS(I).EQ.J) THEN
            JOB = -8
            RINFO(35) = K
            IF (LP.GE.0) WRITE (LP,FMT='(/,A,I6,A,I6)')
     +          ' LA04 error return: two entries for row',I,
     +          ' and column',J
            RETURN
          END IF
          IWS(I) = J
   50   CONTINUE
   60 CONTINUE
      IF (JOB.EQ.0) THEN
        IF (M.GT.0 .AND. N.GT.0) CALL LA04SD(A,LA,IRN,IP,M,N,B,C,BND,KB,
     +                                SIGMA,WS(N+2),WS,WS(M+N+3),IWS)
        JOB = 1
        RETURN
      END IF
      CALL LA15ID(IWS(II15),WS(IC15),IWS(IK15))
      IWS(II15+1) = 0
      DO 70 I = 1,N + M
        X(I) = ZERO
   70 CONTINUE
      IF (CNTL(12).EQ.ZERO) CALL FA14ID(IWS(IISD))
      IF (JOB.EQ.2) THEN
        DO 80 J = 1,N
          G(J) = ONE
   80   CONTINUE
        NPM = N + M
        DO 90 I = 1,M
          IF (IX(I).LT.0 .OR. IX(I).GT.3*NPM) THEN
            IF (LP.GE.0) WRITE (LP,'(/,A,/,A,I6,A,I11)')
     +          ' LA04 error return: IX is faulty',' IX(',I,') =',IX(I)
            RINFO(35) = I
            JOB = -11
            RETURN
          END IF
          IF (IX(I).GT.NPM) THEN
            IX(I) = IX(I) - NPM
            IF (IX(I).GT.NPM) IX(I) = IX(I) - NPM
          END IF
          J = IX(I)
          IF (J.GE.1 .AND. J.LE.N) THEN
            G(J) = -I
            IF (J.LE.KB) THEN
              JX(J) = -1
              IF (BND(1,J).EQ.BND(2,J)) THEN
                G(J) = ZERO
                JX(J) = 1
                IX(I) = 0
              END IF
            END IF
          END IF
   90   CONTINUE
        DO 100 I = 1,KB
          IF (JX(I).EQ.-1) THEN
          ELSE IF (BND(1,I).EQ.BND(2,I)) THEN
            G(I) = ZERO
            JX(I) = 1
          ELSE IF (JX(I).EQ.1) THEN
            IF (BND(1,I).EQ.-SIGMA) GO TO 110
          ELSE IF (JX(I).EQ.2) THEN
            IF (BND(2,I).EQ.SIGMA) GO TO 110
          ELSE IF (JX(I).EQ.0) THEN
            IF (BND(1,I).GE.ZERO .OR. BND(2,I).LE.ZERO) GO TO 110
          ELSE
            GO TO 110
          END IF
  100   CONTINUE
  110   IF (I.LE.KB) THEN
          JOB = -12
          RINFO(35) = I
          IF (LP.GE.0) WRITE (LP,'(/,2A,I6,/,A,I6,A,1P,2E9.1)')
     +        ' LA04 error return: JX is faulty',' for variable I =',I,
     +        ' JX(I) = ',JX(I),' BND(:,I) = ',BND(1,I),BND(2,I)
          RETURN
        END IF
      END IF
      RINFO(23) = ZERO
      DO 120 I = 1,IP(N+1) - 1
        RINFO(23) = MAX(RINFO(23),ABS(A(I)))
  120 CONTINUE
      IF (JOB.EQ.1) THEN
        CALL LA04BD(A,LA,IRN,IP,M,N,BND,KB,CNTL,IX,G,RINFO,IWS(IWB),
     +              IWS(IROW),WS(JV),WS(JY))
        DO 130 J = 1,KB
          IF (G(J).GE.ZERO) THEN
            IF (BND(1,J).GE.ZERO) THEN
              JX(J) = 1
            ELSE IF (BND(2,J).LE.ZERO) THEN
              JX(J) = 2
            ELSE
              JX(J) = 0
            END IF
          ELSE
            JX(J) = -1
          END IF
  130   CONTINUE
      END IF
      JOB = 4
  140 IF (JOB.NE.3) THEN
        CALL LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +              RINFO,X,Z,WS(JB),IB,IWS,IWS(IPB),IWS(IWB),WS(JV),WS,
     +              WS(JY),WS(JYY),IWS(II15+1),WS(IC15+1),IWS(IISD))
      END IF
      IF (JOB.GE.3) THEN
        CALL LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,Z,G,
     +              RINFO,WS(JB),IWS,IWS(IPB),IWS(IWB),WS,WS(JV),
     +              IWS(IR),WS(JY),WS(JYY),IB,IWS(II15+1),WS(IC15+1))
      END IF
      IF (JOB.GE.5) THEN
        CALL LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +              RINFO,X,Z,WS(JB),IB,IWS,IWS(IPB),IWS(IWB),WS(JV),WS,
     +              WS(JY),WS(JYY),IWS(II15+1),WS(IC15+1),IWS(IISD))
        IF (JOB.GE.5) THEN
          CALL LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,Z,
     +                G,RINFO,WS(JB),IWS,IWS(IPB),IWS(IWB),WS,WS(JV),
     +                IWS(IR),WS(JY),WS(JYY),IB,IWS(II15+1),WS(IC15+1))
        END IF
      END IF
  150 RINFO(4) = ZA12AD(DUM) - RINFO(9) + RINFO(4)
      IF (JOB.GT.0) THEN
        IF(JOB.EQ.8) JOB = 4
      ELSE IF (JOB.EQ.0) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,A,I6,A,/,A,1P,E14.6)')
     +      ' LA04: Optimal solution found after',INT(RINFO(30)),
     +      ' basis factorizations',' Objective function = ',RINFO(1)
      ELSE IF (JOB.EQ.-2) THEN
        IF (RINFO(31).GT.ONE) THEN
          IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the objective function is unbounded'
        ELSE
          JOB = -1
          IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the feasible region is void'
        END IF
        DO 160 I = M,1,-1
          WS(N+I) = WS(I)
  160   CONTINUE
        DO 170 J = 1,N
          WS(J) = ZERO
  170   CONTINUE
        JIN = RINFO(26)
        SIGN = ONE
        IF (Z(JIN).LT.ZERO) SIGN = -ONE
        WS(JIN) = -SIGN
        DO 180 I = 1,M
          J = IX(I)
          IF (J.LE.N) WS(J) = WS(N+I)*SIGN
  180   CONTINUE
      ELSE IF (JOB.EQ.-1) THEN
        IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the feasible region is void'
      ELSE IF (JOB.EQ.-4) THEN
        L = RINFO(35)
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +      ' LA04 error return: IB must be at least',L
      END IF
      END
      SUBROUTINE LA04ID(CNTL)
      DOUBLE PRECISION CNTL(15)
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
      CNTL(1) = FD15AD('H')
      CNTL(2) = FD15AD('E')** (2.0D0/3.0D0)
      CNTL(3) = 0
      CNTL(4) = 1
      CNTL(5) = 0
      CNTL(6) = 6
      CNTL(7) = 6
      CNTL(8) = 0.1D0
      CNTL(9) = 1.0D0/SQRT(FD15AD('E'))
      CNTL(10) = FD15AD('E')
      CNTL(11) = CNTL(2)
      CNTL(12) = 0
      END
      SUBROUTINE LA04BD(A,LA,IRN,IP,M,N,BND,KB,CNTL,IX,G,RINFO,IW,IROW,
     +                  PV,Y)
      INTEGER LA,M,N,KB
      DOUBLE PRECISION
     +                A(LA),BND(2,KB),CNTL(15),G(N),RINFO(40),PV(M),Y(M)
      INTEGER IP(N+1),IRN(LA),IROW(M),IW(M),IX(M)
      DOUBLE PRECISION
     +               AMAX,CMAX,ONE,XI,XX,ZERO
      PARAMETER (ONE=1,ZERO=0)
      INTEGER I,I2,ISWP,J,JJ,K,K1,K2,KK,KMAX,MINZ,MP,MZ,NART,NMINZ,NP,NZ
      INTRINSIC ABS,MAX
      MP = CNTL(7)
      IF (MP.GE.0) WRITE (MP,FMT='(/,A/A)')
     +    ' LA04: Sweep      No.       No.     Min col   No. cols',
     +    '              variables artificials length   of min len'
      DO 10 J = 1,N
        G(J) = ONE
   10 CONTINUE
      DO 20 J = 1,KB
        IF (BND(1,J).EQ.BND(2,J)) G(J) = ZERO
   20 CONTINUE
      DO 30 I = 1,M
        Y(I) = ZERO
        IW(I) = 0
   30 CONTINUE
      NART = 0
      NP = 0
      MZ = 1
      DO 70 ISWP = 1,M
        MINZ = M + 1
        NMINZ = 0
        DO 60 J = 1,N
          IF (G(J).LE.ZERO) GO TO 60
          K1 = IP(J)
          K2 = IP(J+1) - 1
          IF (K2.LT.K1) GO TO 60
          AMAX = ZERO
          CMAX = ZERO
          NZ = 0
          DO 40 K = K1,K2
            CMAX = MAX(CMAX,ABS(A(K)))
            KK = IRN(K)
            IF (IW(KK).GT.0) GO TO 40
            IF (A(K).NE.ZERO) NZ = NZ + 1
            IF (AMAX.GE.ABS(A(K))) GO TO 40
            AMAX = ABS(A(K))
            KMAX = K
   40     CONTINUE
          IF (AMAX.LT.CNTL(8)*CMAX) GO TO 60
          IF (NZ.EQ.0) GO TO 60
          IF (NZ.LE.MZ) THEN
            I = IRN(KMAX)
            IW(I) = 1
            NP = NP + 1
            IX(NP) = J
            G(J) = -NP
            PV(NP) = A(KMAX)
            IROW(NP) = I
            DO 50 K = K1,K2
              IF (A(K).EQ.ZERO) GO TO 50
              KK = IRN(K)
              IF (IW(KK).GT.0) GO TO 50
              IW(KK) = 1
              NART = NART + 1
              NP = NP + 1
              IROW(NP) = KK
              PV(NP) = ONE
              IX(NP) = N + KK
   50       CONTINUE
          END IF
          IF (NZ.LE.MINZ) THEN
            IF (NZ.LT.MINZ) THEN
              NMINZ = 1
              MINZ = NZ
            ELSE
              NMINZ = NMINZ + 1
            END IF
          END IF
   60   CONTINUE
        IF (MP.GE.0) WRITE (MP,FMT='(I10,4I10)') ISWP,NP,NART,MINZ,NMINZ
        IF (NP.EQ.M) GO TO 100
        IF (NMINZ.EQ.0) GO TO 80
        MZ = MINZ
        IF (MINZ*NMINZ* (M-ISWP).LE.M-NP) MZ = MZ + 1
   70 CONTINUE
      ISWP = ISWP - 1
   80 DO 90 I = 1,M
        IF (IW(I).GT.0) GO TO 90
        NP = NP + 1
        PV(NP) = ONE
        IROW(NP) = I
        IX(NP) = N + I
   90 CONTINUE
  100 RINFO(34) = ISWP
      IF (CNTL(4).EQ.ONE) THEN
        DO 150 J = 1,N
          IF (G(J).LE.ZERO) GO TO 150
          XX = ONE
          DO 110 K = IP(J),IP(J+1) - 1
            KK = IRN(K)
            Y(KK) = A(K)
  110     CONTINUE
          DO 140 I2 = M,1,-1
            I = IROW(I2)
            IF (Y(I).EQ.ZERO) GO TO 140
            XI = Y(I)/PV(I2)
            XX = XX + XI*XI
            JJ = IX(I2)
            IF (JJ.GT.N) GO TO 130
            DO 120 K = IP(JJ),IP(JJ+1) - 1
              KK = IRN(K)
              Y(KK) = Y(KK) - XI*A(K)
  120       CONTINUE
  130       Y(I) = ZERO
  140     CONTINUE
          G(J) = XX
  150   CONTINUE
      END IF
      END
      SUBROUTINE LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +                  RINFO,X,Z,BA,IB,INDB,IPB,IWB,V,W,Y,YY,I15,C15,
     +                  ISEED)
      INTEGER LA,IB,JOB,KB,LB,M,N,ISEED
      DOUBLE PRECISION
     +                A(LA),B(M),BA(IB),BND(2,KB),C(N),CNTL(15),G(N),
     +                RINFO(40),V(M),W(M),X(N+M),Y(M),YY(M),Z(N),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IRN(LA),IWB(M,8),IX(M),JX(KB),
     +        I15(10)
      DOUBLE PRECISION
     +                AW,BL,BU,DV,DW,DZI,E,F,G5,ONE,SIGMA,T,
     +                TIME,XJ,ZGMAX,ZI,ZERO
      PARAMETER (ONE=1,ZERO=0)
      DOUBLE PRECISION DUM
      LOGICAL FEAS
      INTEGER I,IS,IU,J,JIN,K,K1,K15,K2,L,MP,NINF,NPM
      DOUBLE PRECISION FA14AD
      DOUBLE PRECISION ZA12AD
      EXTERNAL FA14AD,ZA12AD
      EXTERNAL LA15AD,LA15BD,LA04ED
      INTRINSIC ABS,MIN,MAX
      K15 = 3
      MP = CNTL(7)
      SIGMA = CNTL(1)
      NPM = N + M
      DO 10 I = 1,M
        J = IX(I)
        IF (J.GT.NPM) THEN
          J = J - NPM
          IF (J.GT.NPM) J = J - NPM
          IX(I) = J
          IF (J.LE.KB) JX(J) = -1
        END IF
        IF (J.GT.0) X(J) = ZERO
   10 CONTINUE
      IF (JOB.EQ.5) GO TO 110
      TIME = ZA12AD(DUM)
      RINFO(28) = 0
      RINFO(2) = CNTL(8)
      IF (M.EQ.0) THEN
        I15(K15+1) = 0
        I15(K15+2) = 0
        I15(K15+3) = 0
        RINFO(3) = ZERO
      ELSE
        DO 90 IU = 1,10
          L = 1
          DO 30 I = 1,M
            J = IX(I)
            IF (J.LE.0) THEN
            ELSE IF (J.LE.N) THEN
              G(J) = -I
              K1 = IP(J)
              K2 = IP(J+1) - 1
              IF (L+K2-K1.LT.IB) THEN
                DO 20 K = K1,K2
                  BA(L) = A(K)
                  INDB(L,1) = IRN(K)
                  INDB(L,2) = I
                  L = L + 1
   20           CONTINUE
              ELSE
                L = L + K2 - K1 + 1
              END IF
            ELSE
              IF (L.LE.IB) THEN
                BA(L) = ONE
                INDB(L,1) = J - N
                INDB(L,2) = I
              END IF
              L = L + 1
            END IF
   30     CONTINUE
          IF (L.GT.IB) THEN
            JOB = -4
            RINFO(35) = L
            RETURN
          END IF
          T = ZA12AD(DUM)
          C15(1) = CNTL(10)*RINFO(23)
          CALL LA15AD(BA,INDB,L-1,IB,M,IPB,IWB,Y,G5,RINFO(2),I15,
     +                C15,I15(K15+1))
          RINFO(3) = G5
          RINFO(32) = I15(K15+1) + I15(K15+2)
          RINFO(33) = I15(K15+3)
          RINFO(5) = RINFO(5) + ZA12AD(DUM) - T
          RINFO(30) = RINFO(30) + 1
          IF (G5.LT.-6.5) THEN
            JOB = -4
            RINFO(35) = NINT(IB*1.5)
            RETURN
          ELSE IF (G5.LT.-4.5) THEN
            IF (MP.GE.0) WRITE (MP,FMT='(A)')
     +          ' New basis constructed by LA04C'
            I = 0
            K = 0
            DO 80 L = 1,M
              IS = K + 1
              DO 40 K = IS,M
                IF (IWB(K,6).GE.0) GO TO 50
   40         CONTINUE
              GO TO 90
   50         J = IX(K)
              IS = I + 1
              DO 60 I = IS,M
                IF (IWB(I,5).GE.0) GO TO 70
   60         CONTINUE
              GO TO 90
   70         IF (J.GE.1 .AND. J.LE.N) THEN
                G(J) = ONE
                IF (J.LE.KB) THEN
                  IF (BND(1,J).GE.ZERO) THEN
                    JX(J) = 1
                  ELSE IF (BND(2,J).LE.ZERO) THEN
                    JX(J) = 2
                  ELSE
                    JX(J) = 0
                  END IF
                END IF
              END IF
              IX(K) = I + N
   80       CONTINUE
          ELSE
            IF (G5.LT.CNTL(9)*RINFO(23) .OR.
     +          RINFO(2).EQ.ONE) GO TO 100
            RINFO(2) = MIN(ONE,RINFO(2)+RINFO(2))
          END IF
   90   CONTINUE
  100   RINFO(24) = MAX(G5*10.,RINFO(23)*CNTL(9)**1.5)
        IF (MP.GE.0) WRITE (MP,'(/,A,F9.4,/,A,3F11.3,/,A,I9)')
     +     ' LA04: Factorization of basis. Value of U used by LA15A was'
     +      ,RINFO(2),' Total times so far for LA15A/B/C are',RINFO(5),
     +      RINFO(6),RINFO(7),' Size of LA15 arrays is',IB
      END IF
  110 DO 120 I = 1,M
        W(I) = B(I)*FA14AD(ISEED,-1)*CNTL(2)
        V(I) = B(I)
  120 CONTINUE
      DO 140 J = 1,N
        XJ = X(J)
        IF (J.LE.KB) THEN
          IF (JX(J).EQ.1) THEN
            XJ = XJ + BND(1,J)
          ELSE IF (JX(J).EQ.2) THEN
            XJ = XJ + BND(2,J)
          END IF
        END IF
        IF (XJ.NE.ZERO) THEN
          DO 130 K = IP(J),IP(J+1) - 1
            I = IRN(K)
            W(I) = W(I) + A(K)*XJ*FA14AD(ISEED,-1)*CNTL(2)
            V(I) = V(I) - A(K)*XJ
  130     CONTINUE
        END IF
  140 CONTINUE
      DO 150 J = 1,M
        V(J) = V(J) - X(N+J)
  150 CONTINUE
      IF (M.GT.0) THEN
        T = ZA12AD(DUM)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.FALSE.,I15,
     +              I15(K15+1))
        RINFO(6) = RINFO(6) + ZA12AD(DUM) - T
      END IF
      DO 160 I = 1,M
        CNTL(3) = MAX(CNTL(3),ABS(W(I)))
  160 CONTINUE
      F = SIGMA
      DO 170 I = 1,KB
        IF (BND(1,I).NE.ZERO) F = MIN(ABS(BND(1,I)),F)
        IF (BND(2,I).NE.ZERO) F = MIN(ABS(BND(2,I)),F)
  170 CONTINUE
      IF (F.NE.SIGMA) CNTL(3) = MAX(F*CNTL(2),CNTL(3))
      DO 180 I = 1,M
        W(I) = V(I)
  180 CONTINUE
      IF (M.GT.0) THEN
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),V,.FALSE.,I15,
     +              I15(K15+1))
        CALL LA04ED(A,LA,IRN,IP,M,N,W,V,DV,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.FALSE.,I15,C15)
        IF (MP.GE.0) WRITE (MP,'(A,1P,E12.4)')
     +      ' Largest change, iterative refinement of X:',DV
        CNTL(3) = MAX(DV,CNTL(3))
      END IF
      IF (MP.GE.0) WRITE (MP,'(/,A,1P,E9.1)')
     +    ' LA04: tolerance on the bounds is',CNTL(3)
      NINF = 0
      F = ZERO
      DO 190 I = 1,M
        E = V(I)
        W(I) = ZERO
        J = IX(I)
        X(J) = E
        IF (J.GT.N) THEN
          IF (E.LT.-CNTL(3)) THEN
            W(I) = -ONE
            IX(I) = J + NPM
            NINF = NINF + 1
            F = F - E
          ELSE IF (E.GT.CNTL(3)) THEN
            W(I) = ONE
            IX(I) = J + NPM + NPM
            NINF = NINF + 1
            F = F + E
          END IF
        ELSE IF (J.LE.KB) THEN
          BL = BND(1,J)
          BU = BND(2,J)
          IF (BL.NE.-SIGMA) THEN
            IF (E.LT.BL-CNTL(3)) THEN
              E = E - BL
              X(J) = E
              W(I) = -ONE
              IX(I) = J + NPM
              NINF = NINF + 1
              F = F - E
              JX(J) = 1
              GO TO 190
            ELSE IF (E.LT.BL) THEN
              GO TO 190
            END IF
          END IF
          IF (BU.NE.SIGMA) THEN
            IF (E.GT.BU+CNTL(3)) THEN
              W(I) = ONE
              E = E - BU
              X(J) = E
              IX(I) = J + NPM + NPM
              NINF = NINF + 1
              F = F + E
              JX(J) = 2
            END IF
          END IF
        ELSE IF (J.GE.LB) THEN
          IF (E.LT.-CNTL(3)) THEN
            IX(I) = J + NPM
            W(I) = -ONE
            NINF = NINF + 1
            F = F - E
          END IF
        END IF
  190 CONTINUE
      IF (NINF.GT.RINFO(27)) THEN
        IF (CNTL(5).LT.1) THEN
          CNTL(5) = 1
          IF (MP.GE.0) WRITE (MP,'(/,A)')
     +        ' LA04: starting iterative refinement'
        END IF
      END IF
      FEAS = NINF .EQ. 0
      IF (FEAS) THEN
        RINFO(31) = 2
        DO 200 J = 1,N
          F = F + C(J)*X(J)
          IF (G(J).LT.ZERO) THEN
            Z(J) = ZERO
            I = -G(J)
            W(I) = C(J)
          ELSE
            Z(J) = C(J)
          END IF
  200   CONTINUE
        DO 210 J = 1,KB
          IF (JX(J).EQ.1) THEN
            F = F + C(J)*BND(1,J)
          ELSE IF (JX(J).EQ.2) THEN
            F = F + C(J)*BND(2,J)
          END IF
  210   CONTINUE
      ELSE
        RINFO(31) = 1
        DO 220 J = 1,N
          Z(J) = ZERO
  220   CONTINUE
      END IF
      RINFO(1) = F
      RINFO(27) = NINF
      IF (MP.GE.0) THEN
        WRITE (MP,'(/,A,I5)') ' LA04: number of infeasibilites:',NINF
        WRITE (MP,'(A,1P,E14.6)') ' Objective function:',F
      END IF
      IF (M.GT.0) THEN
        DO 230 I = 1,M
          V(I) = W(I)
  230   CONTINUE
        T = ZA12AD(DUM)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.TRUE.,I15,
     +              I15(K15+1))
        RINFO(6) = RINFO(6) + ZA12AD(DUM) - T
        CALL LA04ED(A,LA,IRN,IP,M,N,V,W,DW,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.TRUE.,I15,C15)
        IF (MP.GE.0) WRITE (MP,'(A,1P,E12.4)')
     +      ' Largest change, iterative refinement of reduced costs:',DW
        DO 240 I = 1,M
          IF (ABS(W(I)).LE.DW) W(I) = ZERO
  240   CONTINUE
      END IF
      ZGMAX = ZERO
      JIN = 0
      DO 260 I = 1,N
        IF (G(I).LE.ZERO) THEN
          Z(I) = ZERO
          GO TO 260
        END IF
        ZI = Z(I)
        DZI = ABS(ZI)
        DO 250 K = IP(I),IP(I+1) - 1
          AW = A(K)*W(IRN(K))
          DZI = DZI + ABS(AW)
          ZI = ZI - AW
  250   CONTINUE
        DZI = DZI*CNTL(2)
        IF (I.GE.LB) THEN
          ZI = ZI + DZI
          Z(I) = ZI
          IF (ZI.GE.ZERO) GO TO 260
        ELSE IF (I.GT.KB) THEN
          Z(I) = ZI
          IF (ABS(ZI).LE.DZI) THEN
            Z(I) = ZERO
            GO TO 260
          END IF
        ELSE IF (JX(I).EQ.0) THEN
          Z(I) = ZI
          IF (ABS(ZI).LE.DZI) THEN
            Z(I) = ZERO
            GO TO 260
          END IF
        ELSE IF (JX(I).EQ.1) THEN
          ZI = ZI + DZI
          Z(I) = ZI
          IF (ZI.GE.ZERO) GO TO 260
        ELSE IF (JX(I).EQ.2) THEN
          Z(I) = ZI - DZI
          ZI = -Z(I)
          IF (ZI.GE.ZERO) GO TO 260
        END IF
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  260 CONTINUE
      RINFO(26) = JIN
      RINFO(10) = ZA12AD(DUM)
      IF (JOB.NE.5) RINFO(8) = RINFO(10) - TIME
      RINFO(10) = RINFO(10) - RINFO(9) + RINFO(4)
      IF (ZGMAX.EQ.ZERO) THEN
        IF (FEAS) THEN
          JOB = 0
        ELSE
          JOB = -1
        END IF
      END IF
      END
      SUBROUTINE LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,
     +                Z,G,RINFO,BA,INDB,IPB,IWB,W,V,IR,Y,YY,IB,I15,C15)
      INTEGER LA,IB,JOB,KB,LB,M,N
      DOUBLE PRECISION
     +                A(LA),BA(IB),BND(2,KB),B(M),C(N),CNTL(15),G(N),
     +                RINFO(40),V(M),W(M),X(N+M),Y(M),YY(M),Z(N),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IR(M),IRN(LA),IWB(M,4),IX(M),
     +        JX(KB),I15(10)
      DOUBLE PRECISION
     +                ALPHA,BETA,DV,DW,ERR,GQ,G5,PIVOT,PVMIN,R,RMIN,RP,
     +                RPMIN,SIGMA,SIGN,SIGNX,T,WW,XJ,ZGMAX,ZI,ZQ,
     +                ONE,ZERO
      PARAMETER (ONE=1,ZERO=0)
      DOUBLE PRECISION DUM
      LOGICAL FEAS
      INTEGER I,IK,IRP,ITER,J,JIN,JOUT,JP,K,K1,K15,K2,KK,MP,NINF,NPM
      DOUBLE PRECISION ZA12AD
      EXTERNAL ZA12AD
      EXTERNAL LA04ED,LA15BD,LA15CD
      INTRINSIC ABS,MIN,SQRT,MAX,MOD
      K15 = 3
      MP = CNTL(7)
      SIGMA = CNTL(1)
      JIN = RINFO(26)
      NINF = RINFO(27)
      IK = RINFO(28)
      ITER = RINFO(29)
      FEAS = RINFO(31) .EQ. 2
      JOUT = 0
      NPM = N + M
      IF (JIN.EQ.0) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +      ' LA04: no candidates are available for entry into basis'
        JOB = 5
        GO TO 280
      END IF
      SIGN = ONE
      IF (Z(JIN).GT.ZERO) SIGN = -ONE
      DO 10 I = 1,M
        W(I) = ZERO
   10 CONTINUE
      DO 20 K = IP(JIN),IP(JIN+1) - 1
        KK = IRN(K)
        W(KK) = A(K)
   20 CONTINUE
      IF (M.GT.0) THEN
        T = ZA12AD(DUM)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.FALSE.,I15,
     +              I15(K15+1))
        RINFO(6) = RINFO(6) + ZA12AD(DUM) - T
        IF (CNTL(5).GE.1) THEN
          CALL LA04ED(A,LA,IRN,IP,M,N,B,W,DW,IX,JIN,RINFO,BA,INDB,IPB,
     +                IWB,YY,V,IB,.FALSE.,I15,C15)
        END IF
      ELSE
        DW = ZERO
      END IF
      ZQ = ZERO
      IF (.NOT.FEAS) THEN
        DO 30 I = 1,M
          IF (IX(I).GT.NPM) THEN
            IF (IX(I).GT.NPM+NPM) THEN
              ZQ = ZQ - W(I)
            ELSE
              ZQ = ZQ + W(I)
            END IF
          END IF
   30   CONTINUE
      ELSE
        ZQ = C(JIN)
        DO 40 I = 1,M
          J = IX(I)
          IF (J.LE.N) THEN
            ZQ = ZQ - C(J)*W(I)
          END IF
   40   CONTINUE
      END IF
      ZI = Z(JIN)
      IF (ABS(ZQ-ZI).GE.MAX(ABS(ZQ),ABS(ZI))*0.1) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,1P,2E12.4)')
     +      ' Approximate and newly computed reduced costs are',Z(JIN),
     +      ZQ
        IF (JOB.EQ.3) THEN
          JOB = 5
        ELSE
          JOB = -1
          IF (FEAS) JOB = 0
        END IF
        GO TO 280
      END IF
Compute G(IQ) afresh.
      GQ = ONE
      DO 50 I = 1,M
        GQ = GQ + W(I)**2
   50 CONTINUE
      RPMIN = SIGMA
      RMIN = SIGMA
      PVMIN = CNTL(11)*SQRT(GQ)
      IF (CNTL(5).GE.1) PVMIN = MAX(DW,PVMIN)
      IRP = 0
      DO 60 I = 1,M
        WW = W(I)*SIGN
        J = IX(I)
        IF (J.GT.NPM) THEN
          IF (J.GT.NPM+NPM) THEN
            XJ = X(J-NPM-NPM)
          ELSE
            XJ = X(J-NPM)
          END IF
        ELSE
          XJ = X(J)
        END IF
        RP = ZERO
        IF (WW.GT.ZERO) THEN
          IF (J.LE.KB) THEN
            IF (BND(1,J).NE.-SIGMA) THEN
              R = (XJ-BND(1,J))/WW
              RP = R + CNTL(3)/WW
            END IF
          ELSE IF (J.GT.NPM+NPM) THEN
            R = XJ/WW
            RP = R + CNTL(3)/WW
          ELSE IF (J.GE.LB) THEN
            R = XJ/WW
            RP = R + CNTL(3)/WW
          END IF
        ELSE IF (WW.LT.ZERO) THEN
          IF (J.LE.KB) THEN
            IF (BND(2,J).NE.SIGMA) THEN
              R = (XJ-BND(2,J))/WW
              RP = R - CNTL(3)/WW
            END IF
          ELSE IF (J.GT.N .AND. J.LE.NPM+NPM) THEN
            R = XJ/WW
            RP = R - CNTL(3)/WW
          END IF
        END IF
        IF (RP.GT.ZERO) THEN
          IF (R.LE.RPMIN) THEN
            IF (RP.LT.RMIN) IRP = 0
            RPMIN = MIN(RP,RPMIN)
            RMIN = MIN(R,RMIN)
            IF (ABS(WW).GT.PVMIN) THEN
              IRP = IRP + 1
              IR(IRP) = I
            END IF
          END IF
        END IF
   60 CONTINUE
      IF (JIN.LE.KB) THEN
        R = ZERO
        IF (SIGN.GT.0) THEN
          IF (BND(2,JIN).NE.SIGMA) THEN
            IF (JX(JIN).EQ.1) THEN
              R = BND(2,JIN) - X(JIN) - BND(1,JIN)
            ELSE IF (JX(JIN).EQ.0) THEN
              R = BND(2,JIN) - X(JIN)
            END IF
          END IF
        ELSE
          IF (BND(1,JIN).NE.-SIGMA) THEN
            IF (JX(JIN).EQ.2) THEN
              R = X(JIN) + BND(2,JIN) - BND(1,JIN)
            ELSE IF (JX(JIN).EQ.0) THEN
              R = X(JIN) - BND(1,JIN)
            END IF
          END IF
        END IF
        IF (R.GT.ZERO) THEN
          IF (R.LE.RPMIN) THEN
            RMIN = R
            JOUT = JIN
            PIVOT = ZERO
            GO TO 90
          END IF
        END IF
      END IF
      PIVOT = ZERO
      DO 70 K = 1,IRP
        I = IR(K)
        WW = W(I)*SIGN
        J = IX(I)
        IF (J.GT.NPM) THEN
          IF (J.GT.NPM+NPM) THEN
            XJ = X(J-NPM-NPM)
          ELSE
            XJ = X(J-NPM)
          END IF
        ELSE
          XJ = X(J)
        END IF
        IF (J.GT.N) THEN
          R = XJ/WW
          IF (.NOT.FEAS) WW = WW*10000.0
        ELSE IF (WW.LE.ZERO) THEN
          R = (XJ-BND(2,J))/WW
        ELSE IF (J.GE.LB) THEN
          R = XJ/WW
        ELSE
          R = (XJ-BND(1,J))/WW
        END IF
        IF (R.GT.RPMIN) GO TO 70
        IF (ABS(WW).LE.ABS(PIVOT)) GO TO 70
        PIVOT = WW
        RMIN = R
        JP = I
        JOUT = J
   70 CONTINUE
      IF (PIVOT.EQ.ZERO) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,1P,2E12.4)')
     +      ' LA04: approximate and newly computed reduced costs are',
     +      Z(JIN),ZQ
        IF (JOB.NE.3) THEN
          JOB = -2
          DO 80 I = 1,M
            IF (ABS(W(I)).LE.PVMIN) W(I) = ZERO
   80     CONTINUE
        ELSE
          IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +        ' LA04: no row is acceptable as the pivot row'
          JOB = 5
        END IF
        GO TO 280
      END IF
      PIVOT = W(JP)
   90 SIGNX = 0
      IF (JOUT.GT.NPM) THEN
        SIGNX = -1
        JOUT = JOUT - NPM
        IF (JOUT.GT.NPM) THEN
          SIGNX = 1
          JOUT = JOUT - NPM
        END IF
        NINF = NINF - 1
      END IF
      RMIN = MAX(ZERO,RMIN)*SIGN
      DO 100 I = 1,M
        J = IX(I)
        IF (J.GT.NPM) THEN
          J = J - NPM
          IF (J.GT.NPM) J = J - NPM
        END IF
        X(J) = X(J) - W(I)*RMIN
  100 CONTINUE
      IF (JIN.NE.JOUT) THEN
        IF (M.GT.0) THEN
          T = ZA12AD(DUM)
          IF (CNTL(4).EQ.ONE) CALL LA15BD(BA,INDB,IB,M,IPB,IWB,V,
     +                             RINFO(3),W,.TRUE.,I15,I15(K15+1))
          RINFO(6) = ZA12AD(DUM) - T + RINFO(6)
          T = ZA12AD(DUM)
          G5 = RINFO(3)
          C15(1) = CNTL(10)*RINFO(23)
          CALL LA15CD(BA,INDB,IB,M,IPB,IWB,Y,G5,RINFO(2),JP,I15,
     +                C15,I15(K15+1))
          RINFO(3) = G5
          RINFO(7) = ZA12AD(DUM) - T + RINFO(7)
          IF (G5.GE.RINFO(24)) THEN
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +          ' LA04: refactorization because of LA15 growth'
            JOB = 8
            GO TO 280
          ELSE IF (G5.LT.-6.5) THEN
            IF (IK.LE.5) THEN
              JOB = -4
              RINFO(35) = NINT(IB*1.5)
              GO TO 280
            END IF
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,F7.1)')
     +          ' LA04: refactorization because of LA15 fill-in'
            JOB = 8
            GO TO 280
          ELSE IF (G5.LT.ZERO) THEN
            CNTL(11) = ABS(PIVOT)*1.1/SQRT(GQ)
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A/,A,1PE10.2)')
     +          ' LA04: refactorization to avoid singular basis.',
     +          ' CNTL(11) increased to ', CNTL(11)
            JOB = 8
            GO TO 280
          END IF
        END IF
        X(JIN) = X(JIN) + RMIN
        IF (JIN.LE.KB) THEN
          IF (JX(JIN).EQ.1) THEN
            X(JIN) = X(JIN) + BND(1,JIN)
          ELSE IF (JX(JIN).EQ.2) THEN
            X(JIN) = X(JIN) + BND(2,JIN)
          END IF
        END IF
      END IF
      RINFO(1) = RINFO(1) + ZQ*RMIN
      ITER = ITER + 1
      T = ZA12AD(DUM) - RINFO(9) + RINFO(4)
      IF (MP.GE.0) THEN
        IF (MOD(IK,10).EQ.0) WRITE (MP,FMT='(/,A/,A,A,/,A,A)') ' LA04:',
     +      '   Iter LENL   JIN     NINF       G        Z      ',
     +      '   GAMA       Time',
     +      '        LENU  JOUT  Obj. Fun     NCP    Aprox Z   ',
     +      'Aprox GAMA    Pivot'
        WRITE (MP,FMT='(3I6,I7,1P,E12.1,2E11.3,0P,F11.3)') ITER,
     +    I15(K15+1),JIN,NINF,RINFO(3),ZQ,GQ,T
        WRITE (MP,FMT='(I12,I6,1P,E14.6,I4,E12.3,E11.3,E11.1)')
     +    I15(K15+2),JOUT,RINFO(1),I15(K15+3),Z(JIN),G(JIN),PIVOT
      END IF
      IF (JIN.EQ.JOUT) THEN
        IF (SIGN.LT.ZERO) THEN
          JX(JOUT) = 1
        ELSE
          JX(JOUT) = 2
        END IF
        X(JOUT) = ZERO
        GO TO 170
      END IF
      IF (JOUT.LE.KB) THEN
        IF (SIGNX.EQ.ZERO) THEN
          IF (PIVOT*SIGN.GE.ZERO) THEN
            JX(JOUT) = 1
            X(JOUT) = X(JOUT) - BND(1,JOUT)
          ELSE
            JX(JOUT) = 2
            X(JOUT) = X(JOUT) - BND(2,JOUT)
          END IF
        END IF
      END IF
      IF (JIN.LE.KB) JX(JIN) = -1
      IX(JP) = JIN
      G(JIN) = -JP
      DO 110 I = 1,M
        Y(I) = ZERO
  110 CONTINUE
      Y(JP) = ONE
      IF (M.GT.0) THEN
        T = ZA12AD(DUM)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,V,RINFO(3),Y,.TRUE.,
     +              I15,I15(K15+1))
        RINFO(6) = ZA12AD(DUM) - T + RINFO(6)
      END IF
      IF (CNTL(4).EQ.1) THEN
        DO 140 I = 1,N
          IF (G(I).LE.ZERO) GO TO 140
          ALPHA = ZERO
          K1 = IP(I)
          K2 = IP(I+1) - 1
          DO 120 K = K1,K2
            KK = IRN(K)
            ALPHA = ALPHA + A(K)*Y(KK)
  120     CONTINUE
          IF (ALPHA.EQ.ZERO) GO TO 140
          BETA = ZERO
          DO 130 K = K1,K2
            KK = IRN(K)
            BETA = BETA + A(K)*W(KK)
  130     CONTINUE
          G(I) = MAX(G(I)+ALPHA* (ALPHA*GQ-2.*BETA),1.+ALPHA*ALPHA)
          Z(I) = Z(I) - ALPHA*ZQ
  140   CONTINUE
      ELSE
        DO 160 I = 1,N
          IF (G(I).LE.ZERO) GO TO 160
          ALPHA = ZERO
          K1 = IP(I)
          K2 = IP(I+1) - 1
          DO 150 K = K1,K2
            KK = IRN(K)
            ALPHA = ALPHA + A(K)*Y(KK)
  150     CONTINUE
          Z(I) = Z(I) - ALPHA*ZQ
  160   CONTINUE
      END IF
      IF (JOUT.LE.N) THEN
        ZQ = -ZQ/PIVOT - SIGNX
        Z(JOUT) = ZQ
        G(JOUT) = MAX(GQ/ (PIVOT*PIVOT),ONE)
      END IF
  170 IK = IK + 1
      T = ZA12AD(DUM) - RINFO(9) + RINFO(4)
      J = MOD(IK,10) + 1
      RINFO(10+J) = (RINFO(8)+T-RINFO(10))/IK
      IF (IK.EQ.1) RINFO(21) = RINFO(10+J)
      RINFO(21) = MIN(RINFO(21),RINFO(10+J))
      IF (IK.GE.10) THEN
        ERR = ZERO
        DO 180 I = 2,10
          IF (J.NE.I-1) ERR = MAX(ERR,RINFO(10+I-1)-RINFO(10+I))
  180   CONTINUE
        IF (J.NE.10) ERR = MAX(ERR,RINFO(10+10)-RINFO(10+1))
        IF (RINFO(10+J).GT. (RINFO(21)+ERR)*1.05) THEN
          IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +        ' LA04: refactorization to improve average iteration time'
          JOB = 4
          GO TO 280
        END IF
      END IF
      IF (CNTL(5).GE.1) THEN
        DO 190 I = 1,M
          Y(I) = X(N+I)
  190   CONTINUE
        DO 200 I = 1,M
          W(I) = B(I)
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          V(I) = X(J)
          IF (J.GT.N) Y(J-N) = ZERO
  200   CONTINUE
        DO 220 J = 1,N
          XJ = ZERO
          IF (G(J).GT.ZERO) XJ = X(J)
          IF (J.LE.KB) THEN
            IF (JX(J).EQ.1) THEN
              XJ = XJ + BND(1,J)
            ELSE IF (JX(J).EQ.2) THEN
              XJ = XJ + BND(2,J)
            END IF
          END IF
          IF (XJ.NE.ZERO) THEN
            DO 210 K = IP(J),IP(J+1) - 1
              I = IRN(K)
              W(I) = W(I) - A(K)*XJ
  210       CONTINUE
          END IF
  220   CONTINUE
        DO 230 J = 1,M
          W(J) = W(J) - Y(J)
  230   CONTINUE
        CALL LA04ED(A,LA,IRN,IP,M,N,W,V,DV,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.FALSE.,I15,C15)
        DO 240 I = 1,M
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          X(J) = V(I)
  240   CONTINUE
      END IF
      JIN = 0
      JOB = 3
      IF (.NOT.FEAS) THEN
        IF (RINFO(1).LT.CNTL(3)) GO TO 280
      END IF
      ZGMAX = ZERO
      DO 250 I = 1,KB
        IF (G(I).LE.ZERO) GO TO 250
        ZI = Z(I)
        IF (JX(I).EQ.1) THEN
          IF (ZI.GT.ZERO) GO TO 250
        ELSE IF (JX(I).EQ.2) THEN
          IF (ZI.LT.ZERO) GO TO 250
        END IF
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  250 CONTINUE
      DO 260 I = KB + 1,LB - 1
        IF (G(I).LE.ZERO) GO TO 260
        ZI = Z(I)
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  260 CONTINUE
      DO 270 I = LB,N
        IF (G(I).LE.ZERO) GO TO 270
        ZI = Z(I)
        IF (ZI.GT.ZERO) GO TO 270
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  270 CONTINUE
  280 RINFO(25) = JOUT
      RINFO(26) = JIN
      RINFO(27) = NINF
      RINFO(28) = IK
      RINFO(29) = ITER
      RINFO(32) = I15(K15+1) + I15(K15+2)
      RINFO(33) = I15(K15+3)
      END
      SUBROUTINE LA04ED(A,LA,IRN,IP,M,N,B,X,DX,IX,JIN,RINFO,BA,INDB,IPB,
     +                  IWB,W,Y,IB,TRANS,I15,C15)
      INTEGER LA,IB,M,N
      LOGICAL TRANS
      DOUBLE PRECISION DX
      DOUBLE PRECISION
     +                A(LA),BA(IB),B(M),RINFO(40),W(M),X(M),Y(M),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IRN(LA),IWB(M,4),IX(M),JIN,
     +        I15(10)
      DOUBLE PRECISION ZA12AD
      EXTERNAL ZA12AD
      EXTERNAL LA15BD
      INTEGER I,J,K,NPM,K15
      REAL SA,SX,T
      DOUBLE PRECISION
     +                DUM,OLDSUM,RA,RX,SUM,SUM2,TERM,ZERO
      PARAMETER (ZERO=0)
      K15 = 3
      NPM = N + M
      IF (JIN.EQ.0) THEN
        DO 10 I = 1,M
          W(I) = B(I)
   10   CONTINUE
      ELSE
        DO 20 I = 1,M
          W(I) = ZERO
   20   CONTINUE
        DO 30 K = IP(JIN),IP(JIN+1) - 1
          W(IRN(K)) = A(K)
   30   CONTINUE
      END IF
      IF (TRANS) THEN
        DO 50 I = 1,M
          J = IX(I)
          IF (J.GT.N+M) THEN
            J = J - N - M
            IF (J.GT.N+M) J = J - N - M
          END IF
          IF (J.LE.N) THEN
            SA = W(I)
            SUM = SA
            SUM2 = W(I) - SUM
            DO 40 K = IP(J),IP(J+1) - 1
              RX = -X(IRN(K))
              SX = RX
              RX = RX - SX
              RA = A(K)
              SA = RA
              RA = RA - SA
              TERM = DBLE(SX)*DBLE(SA)
              OLDSUM = SUM
              SUM = SUM + TERM
              TERM = (OLDSUM-SUM) + TERM
              SUM2 = SUM2 + RX*SA + SX*RA + RX*RA  + TERM
   40       CONTINUE
            W(I) = SUM + SUM2
          ELSE
            W(I) = W(I) - X(J-N)
          END IF
   50   CONTINUE
      ELSE
        DO 60 I = 1,M
          Y(I) = ZERO
   60   CONTINUE
        DO 80 I = 1,M
          IF (X(I).EQ.ZERO) GO TO 80
          RX = X(I)
          SX = RX
          RX = RX - SX
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          IF (J.LE.N) THEN
            DO 70 K = IP(J),IP(J+1) - 1
              RA = -A(K)
              SA = RA
              RA = RA - SA
              TERM = DBLE(SX)*DBLE(SA)
              OLDSUM = W(IRN(K))
              W(IRN(K)) = OLDSUM + TERM
              TERM = (OLDSUM-W(IRN(K))) + TERM
              Y(IRN(K)) = Y(IRN(K)) + RX*SA + SX*RA + RX*RA + TERM
   70       CONTINUE
          ELSE
            TERM = -X(I)
            OLDSUM = W(J-N)
            W(J-N) = OLDSUM + TERM
            TERM = (OLDSUM-W(J-N)) + TERM
            Y(J-N) = Y(J-N) + TERM
          END IF
   80   CONTINUE
        DO 90 I = 1,M
          W(I) = W(I) + Y(I)
   90   CONTINUE
      END IF
      T = ZA12AD(DUM)
      CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,TRANS,I15,
     +            I15(K15+1))
      RINFO(6) = RINFO(6) + ZA12AD(DUM) - T
      DX = ZERO
      DO 100 I = 1,M
        DX = MAX(DX,ABS(W(I)))
        X(I) = X(I) + W(I)
  100 CONTINUE
      END
      SUBROUTINE LA04SD(A,LA,IRN,IP,M,N,B,C,BND,KB,SIGMA,RS,CS,WS,JCN)
      INTEGER LA,M,N,KB
      DOUBLE PRECISION
     +                SIGMA
      DOUBLE PRECISION
     +                 A(LA),B(M),BND(2,KB),C(N),CS(N+1),RS(M+1),
     +                 WS(2*M+3*N+5)
      INTEGER IP(N+1),IRN(LA),JCN(LA)
      INTEGER I,J,K,NZ
      DOUBLE PRECISION
     +                ONE,ZERO
      PARAMETER (ONE=1,ZERO=0)
      INTRINSIC ABS,EXP
      EXTERNAL MC29AD
      NZ = IP(N+1) - 1
      DO 20 J = 1,N
        DO 10 K = IP(J),IP(J+1) - 1
          JCN(K) = J
   10   CONTINUE
   20 CONTINUE
      DO 30 I = 1,M
        IF (B(I).NE.ZERO) THEN
          NZ = NZ + 1
          A(NZ) = B(I)
          IRN(NZ) = I
          JCN(NZ) = N + 1
        END IF
   30 CONTINUE
      DO 40 J = 1,N
        IF (C(J).NE.ZERO) THEN
          NZ = NZ + 1
          A(NZ) = C(J)
          IRN(NZ) = M + 1
          JCN(NZ) = J
        END IF
   40 CONTINUE
      IF (NZ.GT.0) THEN
        CALL MC29AD(M+1,N+1,NZ,A,IRN,JCN,RS,CS,WS,0,I)
        DO 50 I = 1,M + 1
          RS(I) = EXP(RS(I))
   50   CONTINUE
        DO 60 J = 1,N + 1
          CS(J) = EXP(CS(J))
   60   CONTINUE
        RS(M+1) = ONE/CS(N+1)
      ELSE
        DO 70 I = 1,M + 1
          RS(I) = ONE
   70   CONTINUE
        DO 80 J = 1,N + 1
          CS(J) = ONE
   80   CONTINUE
      END IF
      DO 90 I = 1,M
        B(I) = B(I)*RS(I)*CS(N+1)
   90 CONTINUE
      DO 100 J = 1,N
        C(J) = C(J)*CS(J)*RS(M+1)
  100 CONTINUE
      DO 110 I = 1,KB
        IF (ABS(BND(1,I)).NE.SIGMA) THEN
          BND(1,I) = BND(1,I)*CS(N+1)/CS(I)
        END IF
        IF (ABS(BND(2,I)).NE.SIGMA) THEN
          BND(2,I) = BND(2,I)*CS(N+1)/CS(I)
        END IF
  110 CONTINUE
      DO 120 K = 1,NZ
        A(K) = A(K)*RS(IRN(K))*CS(JCN(K))
  120 CONTINUE
      END


C LA04_DEPENDENCES
C Original date 15 May 2001
C  May 2001: threadsafe version of LA05
C  Jan 2002: replaced MC20 by MC59
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 19th November 2004 Version 1.1.0. Code added to allow for MC59 not
C        handling an empty matrix.
C 21st February 2005 Version 1.2.0. FD05 dependence changed to FD15.

      SUBROUTINE LA15ID(ICNTL,CNTL,KEEP)
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      CNTL(1) = 0.0
      CNTL(2) = 0.0
      CNTL(3) = 0.0
      DO 10 I = 1,7
      KEEP(I) = 0
   10 CONTINUE
      RETURN
      END
      SUBROUTINE LA15AD(A,IND,NZERO,IA,N,IP,IW,W,G,U,ICNTL,CNTL,KEEP)
      DOUBLE PRECISION G,U
      INTEGER IA,N,NZERO
      DOUBLE PRECISION A(IA),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,8)
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
      DOUBLE PRECISION AM,AMAX,AU,EPS,JCOST,KCOST,NZ1
      INTEGER I,IDUMMY,II,IL,IN,IPP,IPV,IR,J,JP,K,K1,K2,KC,KJ,KK,KL,KLC,
     +        KN,KNP,KP,KPC,KPL,KQ,KR,KRL,KS,L,MCP,NC,NZ,NZC
      INTEGER LP,LENL,LENU,NCP,LROW,LCOL,ICNT59(10),INFO59(10)
      REAL SMALL
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
      EXTERNAL LA15ED,MC59AD
      INTRINSIC ABS,MAX
      LP    = ICNTL(1)
      SMALL = CNTL(1)
      LENL  = KEEP(1)
      LENU  = KEEP(2)
      NCP   = KEEP(3)
      LROW  = KEEP(4)
      LCOL  = KEEP(5)
      EPS = FD15AD('E')
      IF (U.GT.1.0D0) U = 1.0D0
      IF (U.LT.EPS) U = EPS
      IF (N.LT.1) GO TO 640
      IF (NZERO.LT.0) GO TO 700
      G = 0.
      DO 20 I = 1,N
        W(I) = 0.
        DO 10 J = 1,5
          IW(I,J) = 0
   10   CONTINUE
   20 CONTINUE
      L = 1
      LENU = NZERO
      DO 50 IDUMMY = 1,NZERO
        IF (L.GT.LENU) GO TO 60
        DO 30 K = L,LENU
          IF (ABS(A(K)).LE.SMALL) GO TO 40
          I = IND(K,1)
          J = IND(K,2)
          G = MAX(ABS(A(K)),G)
          IF (I.LT.1 .OR. I.GT.N) GO TO 650
          IF (J.LT.1 .OR. J.GT.N) GO TO 650
          IW(I,1) = IW(I,1) + 1
          IW(J,2) = IW(J,2) + 1
   30   CONTINUE
        GO TO 60
   40   L = K
        A(L) = A(LENU)
        IND(L,1) = IND(LENU,1)
        IND(L,2) = IND(LENU,2)
        LENU = LENU - 1
   50 CONTINUE
   60 LENL = 0
      LROW = LENU
      LCOL = LROW
      MCP = MAX(N/10,20)
      NCP = 0
      IF(LENU.GT.0)THEN
        ICNT59(1) = 1
        ICNT59(2) = 0
        ICNT59(3) = 0
        ICNT59(4) = ICNTL(1)
        IF ( ICNT59(4).EQ.0 ) ICNT59(4) = -1
        ICNT59(5) = ICNT59(4)
        ICNT59(6) = 0
        CALL MC59AD(ICNT59,N,N,LENU,IND(1,2),IA,IND(1,1),IA,A,N+1,IP,
     +            N+1,IW(1,7),INFO59)
      ELSE
        DO 70 IR = 1,N
          IP(IR,1) = 1
   70   CONTINUE
      END IF
      K = 1
      DO 80 IR = 1,N
        K = K + IW(IR,2)
        IP(IR,2) = K
   80 CONTINUE
      KL = LENU
      DO 100 II = 1,N
        IR = N + 1 - II
        KP = IP(IR,1)
        DO 90 K = KP,KL
          J = IND(K,2)
          IF (IW(J,5).EQ.IR) GO TO 630
          IW(J,5) = IR
          KR = IP(J,2) - 1
          IP(J,2) = KR
          IND(KR,1) = IR
   90   CONTINUE
        KL = KP - 1
  100 CONTINUE
      DO 120 L = 1,2
        DO 110 I = 1,N
          NZ = IW(I,L)
          IF(NZ.GT.0)THEN
            IN = IW(NZ,L+2)
            IW(NZ,L+2) = I
            IW(I,L+6) = IN
            IF (IN.NE.0) IW(IN,L+4) = I
          END IF
          IW(I,L+4) = 0
  110   CONTINUE
  120 CONTINUE
      DO 560 IPV = 1,N
        JCOST = N
        JCOST = JCOST*JCOST
        DO 210 NZ = 1,N
          NZ1 = NZ - 1
          IF (JCOST.LE.NZ1**2) GO TO 220
          J = IW(NZ,4)
          DO 160 IDUMMY = 1,N
            IF (J.LE.0) GO TO 170
            KP = IP(J,2)
            KL = KP + IW(J,2) - 1
            DO 150 K = KP,KL
              I = IND(K,1)
              KCOST = NZ1* (IW(I,1)-1)
              IF (KCOST.GE.JCOST) GO TO 150
              IF (NZ.EQ.1) GO TO 140
              AMAX = 0.D0
              K1 = IP(I,1)
              K2 = IW(I,1) + K1 - 1
              DO 130 KK = K1,K2
                AMAX = MAX(AMAX,ABS(A(KK)))
                IF (IND(KK,2).EQ.J) KJ = KK
  130         CONTINUE
              IF (ABS(A(KJ)).LT.AMAX*U) GO TO 150
  140         JCOST = KCOST
              IPP = I
              JP = J
              IF (JCOST.LE.NZ1**2) GO TO 220
  150       CONTINUE
            J = IW(J,8)
  160     CONTINUE
  170     I = IW(NZ,3)
          DO 200 IDUMMY = 1,N
            IF (I.LE.0) GO TO 210
            AMAX = 0.D0
            KP = IP(I,1)
            KL = KP + IW(I,1) - 1
            DO 180 K = KP,KL
              AMAX = MAX(ABS(A(K)),AMAX)
  180       CONTINUE
            AU = AMAX*U
            DO 190 K = KP,KL
              IF (ABS(A(K)).LT.AU) GO TO 190
              J = IND(K,2)
              KCOST = NZ1* (IW(J,2)-1)
              IF (KCOST.GE.JCOST) GO TO 190
              JCOST = KCOST
              IPP = I
              JP = J
              IF (JCOST.LE.NZ1**2) GO TO 220
  190       CONTINUE
            I = IW(I,7)
  200     CONTINUE
  210   CONTINUE
        GO TO 680
  220   KP = IP(JP,2)
        KL = IW(JP,2) + KP - 1
        DO 260 L = 1,2
          DO 250 K = KP,KL
            I = IND(K,L)
            IL = IW(I,L+4)
            IN = IW(I,L+6)
            IF (IL.EQ.0) GO TO 230
            IW(IL,L+6) = IN
            GO TO 240
  230       NZ = IW(I,L)
            IW(NZ,L+2) = IN
  240       IF (IN.GT.0) IW(IN,L+4) = IL
  250     CONTINUE
          KP = IP(IPP,1)
          KL = KP + IW(IPP,1) - 1
  260   CONTINUE
        IW(IPP,5) = -IPV
        IW(JP,6) = -IPV
        DO 290 K = KP,KL
          J = IND(K,2)
          KPC = IP(J,2)
          IW(J,2) = IW(J,2) - 1
          KLC = KPC + IW(J,2)
          DO 270 KC = KPC,KLC
            IF (IPP.EQ.IND(KC,1)) GO TO 280
  270     CONTINUE
  280     IND(KC,1) = IND(KLC,1)
          IND(KLC,1) = 0
          IF (J.EQ.JP) KR = K
  290   CONTINUE
        AU = A(KR)
        A(KR) = A(KP)
        A(KP) = AU
        IND(KR,2) = IND(KP,2)
        IND(KP,2) = JP
        NZC = IW(JP,2)
        IF (NZC.EQ.0) GO TO 520
        DO 510 NC = 1,NZC
          KC = IP(JP,2) + NC - 1
          IR = IND(KC,1)
          KR = IP(IR,1)
          KRL = KR + IW(IR,1) - 1
          DO 300 KNP = KR,KRL
            IF (JP.EQ.IND(KNP,2)) GO TO 310
  300     CONTINUE
  310     AM = A(KNP)
          A(KNP) = A(KR)
          A(KR) = AM
          IND(KNP,2) = IND(KR,2)
          IND(KR,2) = JP
          AM = -A(KR)/A(KP)
          IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 320
          IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.
     +        IA) GO TO 670
          CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
          KP = IP(IPP,1)
          KR = IP(IR,1)
  320     KRL = KR + IW(IR,1) - 1
          KQ = KP + 1
          KPL = KP + IW(IPP,1) - 1
          IF (KQ.GT.KPL) GO TO 340
          DO 330 K = KQ,KPL
            J = IND(K,2)
            W(J) = A(K)
  330     CONTINUE
  340     IP(IR,1) = LROW + 1
          IND(KR,2) = 0
          KR = KR + 1
          IF (KR.GT.KRL) GO TO 400
          DO 390 KS = KR,KRL
            J = IND(KS,2)
            AU = A(KS) + AM*W(J)
            IND(KS,2) = 0
            IF (ABS(AU).LE.SMALL) GO TO 350
            G = MAX(G,ABS(AU))
            LROW = LROW + 1
            A(LROW) = AU
            IND(LROW,2) = J
            GO TO 380
  350       LENU = LENU - 1
            K = IP(J,2)
            KL = K + IW(J,2) - 1
            IW(J,2) = KL - K
            DO 360 KK = K,KL
              IF (IND(KK,1).EQ.IR) GO TO 370
  360       CONTINUE
  370       IND(KK,1) = IND(KL,1)
            IND(KL,1) = 0
  380       W(J) = 0.D0
  390     CONTINUE
  400     IF (KQ.GT.KPL) GO TO 490
          DO 480 KS = KQ,KPL
            J = IND(KS,2)
            AU = AM*W(J)
            IF (ABS(AU).LE.SMALL) GO TO 470
            LROW = LROW + 1
            A(LROW) = AU
            IND(LROW,2) = J
            LENU = LENU + 1
            NZ = IW(J,2)
            K = IP(J,2)
            KL = K + NZ - 1
            IF (KL.NE.LCOL) GO TO 410
            IF (LCOL+LENL.GE.IA) GO TO 430
            LCOL = LCOL + 1
            GO TO 420
  410       IF (IND(KL+1,1).NE.0) GO TO 430
  420       IND(KL+1,1) = IR
            GO TO 460
  430       IF (LCOL+LENL+NZ+1.LT.IA) GO TO 440
            IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 670
            CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                  LROW,LCOL)
            K = IP(J,2)
            KL = K + NZ - 1
  440       IP(J,2) = LCOL + 1
            DO 450 KK = K,KL
              LCOL = LCOL + 1
              IND(LCOL,1) = IND(KK,1)
              IND(KK,1) = 0
  450       CONTINUE
            LCOL = LCOL + 1
            IND(LCOL,1) = IR
  460       G = MAX(G,ABS(AU))
            IW(J,2) = NZ + 1
  470       W(J) = 0.D0
  480     CONTINUE
  490     IW(IR,1) = LROW + 1 - IP(IR,1)
          IF (LENL+LCOL+1.LE.IA) GO TO 500
          IF (NCP.GE.MCP) GO TO 670
          CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                LROW,LCOL)
  500     K = IA - LENL
          LENL = LENL + 1
          A(K) = AM
          IND(K,1) = IPP
          IND(K,2) = IR
          LENU = LENU - 1
  510   CONTINUE
  520   K1 = IP(JP,2)
        K2 = IW(JP,2) + K1 - 1
        IW(JP,2) = 0
        DO 550 L = 1,2
          IF (K2.LT.K1) GO TO 540
          DO 530 K = K1,K2
            IR = IND(K,L)
            IF (L.EQ.1) IND(K,L) = 0
            NZ = IW(IR,L)
            IF (NZ.GT.0) THEN
              IN = IW(NZ,L+2)
              IW(IR,L+6) = IN
              IW(NZ,L+2) = IR
              IF (IN.NE.0) IW(IN,L+4) = IR
            END IF
            IW(IR,L+4) = 0
  530     CONTINUE
  540     K1 = IP(IPP,1) + 1
          K2 = IW(IPP,1) + K1 - 2
  550   CONTINUE
  560 CONTINUE
      DO 570 I = 1,N
        J = -IW(I,5)
        IW(J,3) = I
        J = -IW(I,6)
        IW(J,4) = I
        IW(I,2) = 0
  570 CONTINUE
      DO 590 I = 1,N
        KP = IP(I,1)
        KL = IW(I,1) + KP - 1
        DO 580 K = KP,KL
          J = IND(K,2)
          IW(J,2) = IW(J,2) + 1
  580   CONTINUE
  590 CONTINUE
      K = 1
      DO 600 I = 1,N
        K = K + IW(I,2)
        IP(I,2) = K
  600 CONTINUE
      LCOL = K - 1
      DO 620 II = 1,N
        I = IW(II,3)
        KP = IP(I,1)
        KL = IW(I,1) + KP - 1
        DO 610 K = KP,KL
          J = IND(K,2)
          KN = IP(J,2) - 1
          IP(J,2) = KN
          IND(KN,1) = I
  610   CONTINUE
  620 CONTINUE
      GO TO 710
  630 IF (LP.GT.0) WRITE (LP,FMT=9000) IR,J
 9000 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE',/,' THERE IS MOR',
     +       'E THAN ONE ENTRY IN ROW',I5,' AND COLUMN',I5)
      G = -4.
      GO TO 710
  640 IF (LP.GT.0) WRITE (LP,FMT=9010)
 9010 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE N IS NOT POSITIVE')
      G = -1.0D0
      GO TO 710
  650 IF (LP.GT.0) WRITE (LP,FMT=9020) K,I,J
 9020 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE',/,' ELEMENT',I7,
     +       ' IS IN ROW',I5,' AND COLUMN',I5)
      G = -3.D0
      GO TO 710
  670 IF (LP.GT.0) WRITE (LP,FMT=9040)
 9040 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE IA IS TOO SMALL')
      G = -7.D0
      GO TO 710
  680 IPV = IPV - 1
      IF (LP.GT.0) WRITE (LP,FMT=9050) IPV
 9050 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE THE MATRIX IS ',
     +   'SINGULAR WITH RANK ',I7)
      G = -5.D0
      GO TO 710
  700 IF (LP.GT.0) WRITE (LP,FMT=9060)
 9060 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE NZ IS NEGATIVE')
      G = -8.
  710 CONTINUE
      KEEP(1)  = LENL
      KEEP(2)  = LENU
      KEEP(3)  = NCP
      KEEP(4)  = LROW
      KEEP(5)  = LCOL
      RETURN
      END
      SUBROUTINE LA15BD(A,IND,IA,N,IP,IW,W,G,B,TRANS,ICNTL,KEEP)
      DOUBLE PRECISION G
      INTEGER IA,N
      LOGICAL TRANS
      DOUBLE PRECISION A(IA),B(N),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,4)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
      DOUBLE PRECISION AM
      INTEGER I,II,J,K,K2,KK,KL,KLL,KP,KPC,L1,N1,NZ,LP,LENL
      INTRINSIC ABS
      LP    = ICNTL(1)
      LENL  = KEEP(1)
      IF (G.LT.0.D0) GO TO 130
      KLL = IA - LENL + 1
      IF (TRANS) GO TO 80
      IF (LENL.LE.0) GO TO 20
      L1 = IA + 1
      DO 10 KK = 1,LENL
        K = L1 - KK
        I = IND(K,1)
        IF (B(I).EQ.0.D0) GO TO 10
        J = IND(K,2)
        B(J) = B(J) + A(K)*B(I)
   10 CONTINUE
   20 DO 30 I = 1,N
        W(I) = B(I)
        B(I) = 0.D0
   30 CONTINUE
      N1 = N + 1
      DO 70 II = 1,N
        I = N1 - II
        I = IW(I,3)
        AM = W(I)
        KP = IP(I,1)
        IF (KP.GT.0) GO TO 50
        KP = -KP
        IP(I,1) = KP
        NZ = IW(I,1)
        KL = KP - 1 + NZ
        K2 = KP + 1
        DO 40 K = K2,KL
          J = IND(K,2)
          AM = AM - A(K)*B(J)
   40   CONTINUE
   50   IF (AM.EQ.0.D0) GO TO 70
        J = IND(KP,2)
        B(J) = AM/A(KP)
        KPC = IP(J,2)
        KL = IW(J,2) + KPC - 1
        IF (KL.EQ.KPC) GO TO 70
        K2 = KPC + 1
        DO 60 K = K2,KL
          I = IND(K,1)
          IP(I,1) = -ABS(IP(I,1))
   60   CONTINUE
   70 CONTINUE
      GO TO 140
   80 DO 90 I = 1,N
        W(I) = B(I)
        B(I) = 0.D0
   90 CONTINUE
      DO 110 II = 1,N
        I = IW(II,4)
        AM = W(I)
        IF (AM.EQ.0.D0) GO TO 110
        J = IW(II,3)
        KP = IP(J,1)
        AM = AM/A(KP)
        B(J) = AM
        KL = IW(J,1) + KP - 1
        IF (KP.EQ.KL) GO TO 110
        K2 = KP + 1
        DO 100 K = K2,KL
          I = IND(K,2)
          W(I) = W(I) - AM*A(K)
  100   CONTINUE
  110 CONTINUE
      IF (KLL.GT.IA) RETURN
      DO 120 K = KLL,IA
        J = IND(K,2)
        IF (B(J).EQ.0.D0) GO TO 120
        I = IND(K,1)
        B(I) = B(I) + A(K)*B(J)
  120 CONTINUE
      GO TO 140
  130 IF (LP.GT.0) WRITE (LP,FMT=9000)
 9000 FORMAT (/,/,' ERROR RETURN FROM LA15BD BECAUSE EARLIER ENTRY',
     +       ' GAVE ERROR RETURN')
  140 RETURN
      END
      SUBROUTINE LA15CD(A,IND,IA,N,IP,IW,W,G,U,MM,ICNTL,CNTL,KEEP)
      DOUBLE PRECISION G,U
      INTEGER IA,MM,N
      DOUBLE PRECISION A(IA),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,4)
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
      DOUBLE PRECISION AM,AU
      INTEGER I,II,IJ,IM,IN,INS,IPP,IR,IS,J,JM,JNS,JP,K,KJ,KK,KL,KM,KNP,
     +        KP,KPL,KQ,KR,KRL,KS,L,LAST,LAST1,LAST2,M,M1,MCP,NZ
      INTEGER LP,LENL,LENU,NCP,LROW,LCOL
      REAL SMALL
      EXTERNAL LA15ED
      INTRINSIC ABS,MAX
      LP    = ICNTL(1)
      SMALL = CNTL(1)
      LENL  = KEEP(1)
      LENU  = KEEP(2)
      NCP   = KEEP(3)
      LROW  = KEEP(4)
      LCOL  = KEEP(5)
      IF (G.LT.0.D0) GO TO 620
      IF (MM.LT.1 .OR. MM.GT.N) GO TO 630
      JM = MM
      MCP = NCP + 20
      LENU = LENU - IW(JM,2)
      KP = IP(JM,2)
      IM = IND(KP,1)
      KL = KP + IW(JM,2) - 1
      IW(JM,2) = 0
      DO 30 K = KP,KL
        I = IND(K,1)
        IND(K,1) = 0
        KR = IP(I,1)
        NZ = IW(I,1) - 1
        IW(I,1) = NZ
        KRL = KR + NZ
        DO 10 KM = KR,KRL
          IF (IND(KM,2).EQ.JM) GO TO 20
   10   CONTINUE
   20   A(KM) = A(KRL)
        IND(KM,2) = IND(KRL,2)
        IND(KRL,2) = 0
   30 CONTINUE
      DO 110 II = 1,N
        I = IW(II,3)
        IF (I.EQ.IM) M = II
        IF (ABS(W(I)).LE.SMALL) GO TO 100
        LENU = LENU + 1
        LAST = II
        IF (LCOL+LENL.LT.IA) GO TO 40
        IF (NCP.GE.MCP .OR. LENL+LENU.GE.IA) GO TO 610
        CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,LROW,
     +              LCOL)
   40   LCOL = LCOL + 1
        NZ = IW(JM,2)
        IF (NZ.EQ.0) IP(JM,2) = LCOL
        IW(JM,2) = NZ + 1
        IND(LCOL,1) = I
        NZ = IW(I,1)
        KPL = IP(I,1) + NZ
        IF (KPL.GT.LROW) GO TO 50
        IF (IND(KPL,2).EQ.0) GO TO 90
   50   IF (LENL+LROW+NZ.LT.IA) GO TO 60
        IF (NCP.GE.MCP .OR. LENL+LENU+NZ.GE.IA) GO TO 610
        CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
   60   KP = IP(I,1)
        IP(I,1) = LROW + 1
        IF (NZ.EQ.0) GO TO 80
        KPL = KP + NZ - 1
        DO 70 K = KP,KPL
          LROW = LROW + 1
          A(LROW) = A(K)
          IND(LROW,2) = IND(K,2)
          IND(K,2) = 0
   70   CONTINUE
   80   LROW = LROW + 1
        KPL = LROW
   90   IW(I,1) = NZ + 1
        A(KPL) = W(I)
        IND(KPL,2) = JM
  100   W(I) = 0.D0
  110 CONTINUE
      IF (IW(IM,1).EQ.0 .OR. IW(JM,2).EQ.0 .OR. M.GT.LAST) GO TO 600
      INS = M
      M1 = M
      W(JM) = 1.0D0
      DO 140 II = M,LAST
        I = IW(II,3)
        J = IW(II,4)
        IF (W(J).EQ.0.D0) GO TO 130
        KP = IP(I,1)
        KL = KP + IW(I,1) - 1
        DO 120 K = KP,KL
          J = IND(K,2)
          W(J) = 1.0D0
  120   CONTINUE
        IW(INS,4) = I
        INS = INS + 1
        GO TO 140
  130   IW(M1,3) = I
        M1 = M1 + 1
  140 CONTINUE
      IJ = M + 1
      IF (M1.GE.LAST) GO TO 160
      LAST1 = LAST - 1
      DO 150 II = M1,LAST1
        IW(II,3) = IW(IJ,4)
        IJ = IJ + 1
  150 CONTINUE
  160 IW(LAST,3) = IM
      LAST1 = LAST
      JNS = LAST
      W(IM) = 2.D0
      J = JM
      DO 190 IJ = M1,LAST
        II = LAST + M1 - IJ
        I = IW(II,3)
        IF (W(I).NE.2.D0) GO TO 180
        K = IP(I,1)
        IF (II.NE.LAST) J = IND(K,2)
        KP = IP(J,2)
        KL = KP + IW(J,2) - 1
        IW(JNS,4) = I
        JNS = JNS - 1
        DO 170 K = KP,KL
          I = IND(K,1)
          W(I) = 2.D0
  170   CONTINUE
        GO TO 190
  180   IW(LAST1,3) = I
        LAST1 = LAST1 - 1
  190 CONTINUE
      DO 200 II = M1,LAST1
        JNS = JNS + 1
        I = IW(JNS,4)
        W(I) = 3.D0
        IW(II,3) = I
  200 CONTINUE
      DO 240 II = M1,LAST1
        KP = IP(JM,2)
        KL = KP + IW(JM,2) - 1
        IS = 0
        DO 210 K = KP,KL
          L = IND(K,1)
          IF (W(L).NE.3.D0) GO TO 210
          IF (IS.NE.0) GO TO 250
          I = L
          KNP = K
          IS = 1
  210   CONTINUE
        IF (IS.EQ.0) GO TO 600
        IND(KNP,1) = IND(KP,1)
        IND(KP,1) = I
        KP = IP(I,1)
        DO 220 K = KP,IA
          IF (IND(K,2).EQ.JM) GO TO 230
  220   CONTINUE
  230   AM = A(KP)
        A(KP) = A(K)
        A(K) = AM
        IND(K,2) = IND(KP,2)
        IND(KP,2) = JM
        JM = IND(K,2)
        IW(II,4) = I
        W(I) = 2.D0
  240 CONTINUE
      II = LAST1
      GO TO 270
  250 IN = M1
      DO 260 IJ = II,LAST1
        IW(IJ,4) = IW(IN,3)
        IN = IN + 1
  260 CONTINUE
  270 LAST2 = LAST1 - 1
      IF (M1.EQ.LAST1) GO TO 580
      DO 280 I = M1,LAST2
        IW(I,3) = IW(I,4)
  280 CONTINUE
      M1 = II
      IF (M1.EQ.LAST1) GO TO 580
      DO 290 I = 1,N
        W(I) = 0.D0
  290 CONTINUE
      IR = IW(LAST1,3)
      DO 570 II = M1,LAST1
        IPP = IW(II,3)
        KP = IP(IPP,1)
        KR = IP(IR,1)
        JP = IND(KP,2)
        IF (II.EQ.LAST1) JP = JM
        KRL = KR + IW(IR,1) - 1
        DO 300 KNP = KR,KRL
          IF (JP.EQ.IND(KNP,2)) GO TO 310
  300   CONTINUE
        IF(II-LAST1.NE.0) GO TO 570
        IF(II-LAST1.EQ.0) GO TO 600
  310   AM = A(KNP)
        A(KNP) = A(KR)
        A(KR) = AM
        IND(KNP,2) = IND(KR,2)
        IND(KR,2) = JP
        IF (II.EQ.LAST1) GO TO 320
        IF (ABS(A(KP)).LT.U*ABS(AM)) GO TO 320
        IF (ABS(AM).LT.U*ABS(A(KP))) GO TO 350
        IF (IW(IPP,1).LE.IW(IR,1)) GO TO 350
  320   IW(LAST1,3) = IPP
        IW(II,3) = IR
        IR = IPP
        IPP = IW(II,3)
        K = KR
        KR = KP
        KP = K
        KJ = IP(JP,2)
        DO 330 K = KJ,IA
          IF (IND(K,1).EQ.IPP) GO TO 340
  330   CONTINUE
  340   IND(K,1) = IND(KJ,1)
        IND(KJ,1) = IPP
  350   IF (A(KP).EQ.0.D0) GO TO 600
        IF (II.EQ.LAST1) GO TO 570
        AM = -A(KR)/A(KP)
        IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 360
        IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.
     +      IA) GO TO 610
        CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
        KP = IP(IPP,1)
        KR = IP(IR,1)
  360   KRL = KR + IW(IR,1) - 1
        KQ = KP + 1
        KPL = KP + IW(IPP,1) - 1
        IF (KQ.GT.KPL) GO TO 380
        DO 370 K = KQ,KPL
          J = IND(K,2)
          W(J) = A(K)
  370   CONTINUE
  380   IP(IR,1) = LROW + 1
        IND(KR,2) = 0
        KR = KR + 1
        IF (KR.GT.KRL) GO TO 440
        DO 430 KS = KR,KRL
          J = IND(KS,2)
          AU = A(KS) + AM*W(J)
          IND(KS,2) = 0
          IF (ABS(AU).LE.SMALL) GO TO 390
          G = MAX(G,ABS(AU))
          LROW = LROW + 1
          A(LROW) = AU
          IND(LROW,2) = J
          GO TO 420
  390     LENU = LENU - 1
          K = IP(J,2)
          KL = K + IW(J,2) - 1
          IW(J,2) = KL - K
          DO 400 KK = K,KL
            IF (IND(KK,1).EQ.IR) GO TO 410
  400     CONTINUE
  410     IND(KK,1) = IND(KL,1)
          IND(KL,1) = 0
  420     W(J) = 0.D0
  430   CONTINUE
  440   IF (KQ.GT.KPL) GO TO 530
        DO 520 KS = KQ,KPL
          J = IND(KS,2)
          AU = AM*W(J)
          IF (ABS(AU).LE.SMALL) GO TO 510
          LROW = LROW + 1
          A(LROW) = AU
          IND(LROW,2) = J
          LENU = LENU + 1
          NZ = IW(J,2)
          K = IP(J,2)
          KL = K + NZ - 1
          IF (KL.NE.LCOL) GO TO 450
          IF (LCOL+LENL.GE.IA) GO TO 470
          LCOL = LCOL + 1
          GO TO 460
  450     IF (IND(KL+1,1).NE.0) GO TO 470
  460     IND(KL+1,1) = IR
          GO TO 500
  470     IF (LCOL+LENL+NZ+1.LT.IA) GO TO 480
          IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 610
          CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                LROW,LCOL)
          K = IP(J,2)
          KL = K + NZ - 1
  480     IP(J,2) = LCOL + 1
          DO 490 KK = K,KL
            LCOL = LCOL + 1
            IND(LCOL,1) = IND(KK,1)
            IND(KK,1) = 0
  490     CONTINUE
          LCOL = LCOL + 1
          IND(LCOL,1) = IR
  500     G = MAX(G,ABS(AU))
          IW(J,2) = NZ + 1
  510     W(J) = 0.D0
  520   CONTINUE
  530   IW(IR,1) = LROW + 1 - IP(IR,1)
        IF (LENL+LCOL+1.LE.IA) GO TO 540
        IF (NCP.GE.MCP) GO TO 610
        CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,LROW,
     +              LCOL)
  540   K = IA - LENL
        LENL = LENL + 1
        A(K) = AM
        IND(K,1) = IPP
        IND(K,2) = IR
        KP = IP(JP,2)
        NZ = IW(JP,2) - 1
        KL = KP + NZ
        DO 550 K = KP,KL
          IF (IND(K,1).EQ.IR) GO TO 560
  550   CONTINUE
  560   IND(K,1) = IND(KL,1)
        IW(JP,2) = NZ
        IND(KL,1) = 0
        LENU = LENU - 1
  570 CONTINUE
  580 DO 590 II = M,LAST
        I = IW(II,3)
        K = IP(I,1)
        J = IND(K,2)
        IW(II,4) = J
  590 CONTINUE
      GO TO 640
  600 IF (LP.NE.0) WRITE (LP,FMT=9000) MM
 9000 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE',/,' SINGULAR MAT',
     +       'RIX CREATED BY REPLACEMENT OF COL',I5)
      G = -6.D0
      GO TO 640
  610 IF (LP.NE.0) WRITE (LP,FMT=9010)
 9010 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE IA IS TOO SMALL')
      G = -7.D0
      GO TO 640
  620 IF (LP.NE.0) WRITE (LP,FMT=9020)
 9020 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE EARLIER ENTRY',
     +       ' GAVE ERROR RETURN')
      GO TO 640
  630 IF (LP.NE.0) WRITE (LP,FMT=9030) MM
 9030 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE M HAS THE VALUE',
     +       I8)
      G = -9.
  640 CONTINUE
      KEEP(1)  = LENL
      KEEP(2)  = LENU
      KEEP(3)  = NCP
      KEEP(4)  = LROW
      KEEP(5)  = LCOL
      RETURN
      END
      SUBROUTINE LA15ED(A,IRN,IP,N,IW,IA,REALS,NCP,LROW,LCOL)
      INTEGER IA,N,NCP,LROW,LCOL
      LOGICAL REALS
      DOUBLE PRECISION A(IA)
      INTEGER IP(N),IRN(IA),IW(N)
      INTEGER IPI,J,K,KL,KN,NZ
      NCP = NCP + 1
      DO 10 J = 1,N
        NZ = IW(J)
        IF (NZ.LE.0) GO TO 10
        K = IP(J) + NZ - 1
        IW(J) = IRN(K)
        IRN(K) = -J
   10 CONTINUE
      KN = 0
      IPI = 0
      KL = LCOL
      IF (REALS) KL = LROW
      DO 30 K = 1,KL
        IF (IRN(K).EQ.0) GO TO 30
        KN = KN + 1
        IF (REALS) A(KN) = A(K)
        IF (IRN(K).GE.0) GO TO 20
        J = -IRN(K)
        IRN(K) = IW(J)
        IP(J) = IPI + 1
        IW(J) = KN - IPI
        IPI = KN
   20   IRN(KN) = IRN(K)
   30 CONTINUE
      IF (REALS) LROW = KN
      IF (.NOT.REALS) LCOL = KN
      RETURN
      END
* COPYRIGHT (c) 1975 AEA Technology
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       Array lengths given explicitly eg A(MAXA)
C
      SUBROUTINE MC20AD(NC,MAXA,A,INUM,JPTR,JNUM,JDISP)
C
C     .. Scalar Arguments ..
      INTEGER JDISP,MAXA,NC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(MAXA)
      INTEGER INUM(MAXA),JNUM(MAXA),JPTR(NC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JA,JB,JCE,JCEP,K,KR,LOC,NULL
C     ..
C     .. Executable Statements ..
      NULL = -JDISP
C**      CLEAR JPTR
      DO 10 J = 1,NC
        JPTR(J) = 0
   10 CONTINUE
C**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
      DO 20 K = 1,MAXA
        J = JNUM(K) + JDISP
        JPTR(J) = JPTR(J) + 1
   20 CONTINUE
C**      SET THE JPTR ARRAY
      K = 1
      DO 30 J = 1,NC
        KR = K + JPTR(J)
        JPTR(J) = K
        K = KR
   30 CONTINUE
C
C**      REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM IS AN
C        IN-PLACE SORT AND IS OF ORDER MAXA.
      DO 50 I = 1,MAXA
C        ESTABLISH THE CURRENT ENTRY.
        JCE = JNUM(I) + JDISP
        IF (JCE.EQ.0) GO TO 50
        ACE = A(I)
        ICE = INUM(I)
C        CLEAR THE LOCATION VACATED.
        JNUM(I) = NULL
C        CHAIN FROM CURRENT ENTRY TO STORE ITEMS.
        DO 40 J = 1,MAXA
C        CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT
C        POSITION TO STORE ENTRY.
          LOC = JPTR(JCE)
          JPTR(JCE) = JPTR(JCE) + 1
C        SAVE CONTENTS OF THAT LOCATION.
          ACEP = A(LOC)
          ICEP = INUM(LOC)
          JCEP = JNUM(LOC)
C        STORE CURRENT ENTRY.
          A(LOC) = ACE
          INUM(LOC) = ICE
          JNUM(LOC) = NULL
C        CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED.
          IF (JCEP.EQ.NULL) GO TO 50
C        IT DOES.  COPY INTO CURRENT ENTRY.
          ACE = ACEP
          ICE = ICEP
          JCE = JCEP + JDISP
   40   CONTINUE
C
   50 CONTINUE
C
C**      RESET JPTR VECTOR.
      JA = 1
      DO 60 J = 1,NC
        JB = JPTR(J)
        JPTR(J) = JA
        JA = JB
   60 CONTINUE
      RETURN

      END
      SUBROUTINE MC20BD(NC,MAXA,A,INUM,JPTR)
C     . .
C     .. Scalar Arguments ..
      INTEGER MAXA,NC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(MAXA)
      INTEGER INUM(MAXA),JPTR(NC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      KMAX = MAXA
      DO 50 JJ = 1,NC
        J = NC + 1 - JJ
        KLO = JPTR(J) + 1
        IF (KLO.GT.KMAX) GO TO 40
        KOR = KMAX
        DO 30 KDUMMY = KLO,KMAX
C ITEMS KOR, KOR+1, .... ,KMAX ARE IN ORDER
          ACE = A(KOR-1)
          ICE = INUM(KOR-1)
          DO 10 K = KOR,KMAX
            IK = INUM(K)
            IF (IABS(ICE).LE.IABS(IK)) GO TO 20
            INUM(K-1) = IK
            A(K-1) = A(K)
   10     CONTINUE
          K = KMAX + 1
   20     INUM(K-1) = ICE
          A(K-1) = ACE
          KOR = KOR - 1
   30   CONTINUE
C        NEXT COLUMN
   40   KMAX = KLO - 2
   50 CONTINUE
      RETURN

      END
* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      INTRINSIC MAX
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF
      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF
      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF
      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF
      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF
      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0
      PART = 0
      IF (ICNTL6.NE.0) PART = 1
      IF (ICNTL2.EQ.0) THEN
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)
      ELSE IF (ICNTL2.EQ.1) THEN
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)
      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)
      END IF
      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70
   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70
   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN
 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130
      ELSE
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF
      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE
      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (LA.GT.1) THEN
        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE
        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE
        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************
      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER
      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Hyprotech UK and
* Council for the Central Laboratory of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
C Original date March 1993
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MIN
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE,SMIN,ZERO
      PARAMETER (ONE=1D0,SMIN=0.1,ZERO=0D0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3
      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE
      DO 60 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (S.LE.SM) GO TO 160
      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (S.LE.SM) E = ZERO
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (S.LE.SM) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (S.LE.SM) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (S.LE.SM) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29AD ****',' IFAIL =',IFAIL
      END
* *******************************************************************
* COPYRIGHT (c) 1979 Hyprotech UK and
* Council for the Central Laboratory of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
C Original date 23 March 2001
C  March 2001: threadsafe version of FA04
C 20/2/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      DOUBLE PRECISION FUNCTION FA14AD(IX,I)
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
      INTEGER IX,I
      DOUBLE PRECISION X
      INTEGER FHI,K,LEFTLO,XALO,XHI
      INTRINSIC FLOAT
      XHI = IX/B16
      XALO = (IX-XHI*B16)*A
      LEFTLO = XALO/B16
      FHI = XHI*A + LEFTLO
      K = FHI/B15
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
      IF (IX.LT.0) IX = IX + P
      XHI = IX/B16
      X = (FLOAT(XHI)*65536.0D0) + FLOAT(IX-XHI*B16)
      IF (I.GE.0) FA14AD = X*4.6566128752457969241D-10
      IF (I.LT.0) FA14AD = X*9.3132257504915938482D-10 - 1.0D0
      RETURN
      END
      SUBROUTINE FA14BD(IX,MAX,NRAND)
      INTEGER A,B15,B16,P
      PARAMETER (A=16807,B15=32768,B16=65536,P=2147483647)
      INTEGER B30,Q
      PARAMETER (B30=1073741824,Q=1073741823)
      INTEGER IX,MAX,NRAND
      INTEGER BE1,BE2,C,D,F,FHI,G,K,LEFTLO,MHI,MLO,MU,NU,XALO,XHI,XLO
      XHI = IX/B16
      XALO = (IX-XHI*B16)*A
      LEFTLO = XALO/B16
      FHI = XHI*A + LEFTLO
      K = FHI/B15
      IX = (((XALO-LEFTLO*B16)-P)+ (FHI-K*B15)*B16) + K
      IF (IX.LT.0) IX = IX + P
      XHI = IX/B15
      XLO = IX - B15*XHI
      MHI = MAX/B15
      MLO = MAX - B15*MHI
      F = (XHI*MLO-P) + XLO*MHI
      IF (F.GT.0) GO TO 1
      F = F + P
      BE1 = F/B15
      BE2 = F - BE1*B15
      GO TO 2
    1 F = F - 1
      BE1 = F/B15
      BE2 = F - BE1*B15
      BE1 = BE1 + B16
    2 G = B15*BE2 + XLO*MLO
      D = G/B30
      C = XHI/2
      F = ((2* (C*MHI-Q)-1)+MHI* (XHI-2*C)) + D + BE1
      IF (F.GT.0) GO TO 3
      F = F + P
      NU = F/2
      MU = F - NU*2
      GO TO 4
    3 F = F - 1
      NU = F/2
      MU = F - 2*NU
      NU = NU + B30
    4 F = (B30*MU-P) + NU + (G-B30*D)
      NRAND = NU + 1
      IF (F.GE.0) NRAND = NRAND + 1
      RETURN
      END
      SUBROUTINE FA14CD(IX,IGEN)
      INTEGER IX,IGEN
      IGEN = IX
      RETURN
      END
      SUBROUTINE FA14DD(IX,IGEN)
      INTEGER IX,IGEN
      IX = IGEN
      RETURN
      END
      SUBROUTINE FA14ID(IX)
      INTEGER IX
      IX = 1
      RETURN
      END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END
* COPYRIGHT (c) 1989 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 28th feb 2005

C 28th February 2005 Version 1.0.0. New routine to replace za02.

      DOUBLE PRECISION FUNCTION ZA12AD(X)
      DOUBLE PRECISION X

C  Returns the current CPU time in seconds when compiled
C  with any Fortran 95 compiler.

      CALL CPU_TIME(ZA12AD)

      END

