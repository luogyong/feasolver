* *******************************************************************
* COPYRIGHT (c) 1988 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 21 Sept 1994
C
C
      SUBROUTINE MC40AD(ITYPE,N,NNZ,IRN,JCN,ICPTR,IPERM,IW,IPROF,IFLAG)
      INTEGER IFLAG,ITYPE,N,NNZ
      INTEGER ICPTR(N+1),IPERM(N),IRN(2*NNZ),IW(3*N+2),JCN(*)
      REAL(8)::IPROF(2)
      INTEGER LP,MP
      INTEGER I,IOUT,J,J1,JFLAG,KZ,L,NDIAG,OFDIAG
      LOGICAL YESA
      DOUBLE PRECISION A(1)
      EXTERNAL MC34AD,MC49AD,MC40BD
      EXTERNAL MC40KD
      COMMON /MC40ID/LP,MP
      SAVE /MC40ID/
      YESA = .FALSE.
      IFLAG = 0
      IF (N.LT.1) THEN
          IFLAG = -1
          IF (LP.GT.0) WRITE (LP,FMT=9000) IFLAG
          RETURN
      END IF
      IF (NNZ.LT.0) THEN
          IFLAG = -2
          IF (LP.GT.0) WRITE (LP,FMT=9000) IFLAG
          RETURN
      END IF
      IF (ITYPE.NE.0 .AND. ITYPE.NE.1) THEN
          IFLAG = -3
          IF (LP.GT.0) WRITE (LP,FMT=9000) IFLAG
          RETURN
      END IF
      IF (N.EQ.1) THEN
          IPROF(1) = 1
          IPROF(2) = 1
          IPERM(1) = 1
          RETURN
      END IF
      IF (NNZ.EQ.0) THEN
          IPROF(1) = N
          IPROF(2) = N
          RETURN
      END IF
      DO 10 I = 1,3*N + 2
          IW(I) = 0
   10 CONTINUE
      KZ = 0
      IOUT = 0
      NDIAG = 0
      OFDIAG = 0
      IF (ITYPE.EQ.0) THEN
          J1 = ICPTR(1)
          ICPTR(1) = 1
          DO 30 L = 1,N
              DO 20 J = J1,ICPTR(L+1) - 1
                  I = IRN(J)
                  IF (I.GT.N .OR. I.LT.1) IOUT = IOUT + 1
                  IF (I.EQ.L) NDIAG = NDIAG + 1
                  IF (I.LT.L) OFDIAG = OFDIAG + 1
                  IF (I.GT.L .AND. I.LE.N) THEN
                      KZ = KZ + 1
                      IRN(KZ) = I
                  END IF
   20         CONTINUE
              J1 = ICPTR(L+1)
              ICPTR(L+1) = KZ + 1
   30     CONTINUE
      ELSE
          DO 40 L = 1,NNZ
              I = IRN(L)
              J = JCN(L)
              IF (I.GT.N .OR. I.LT.1) IOUT = IOUT + 1
              IF (J.GT.N .OR. J.LT.1) IOUT = IOUT + 1
              IF (I.EQ.J) NDIAG = NDIAG + 1
              IF (I.LT.J) OFDIAG = OFDIAG + 1
              IF (J.GE.1 .AND. I.GT.J .AND. I.LE.N) THEN
                  KZ = KZ + 1
                  IRN(KZ) = I
                  JCN(KZ) = J
              END IF
   40     CONTINUE
      END IF
      IF (IOUT.GT.0) THEN
          IFLAG = 3
      ELSE IF (OFDIAG.GT.0) THEN
          IFLAG = 2
      ELSE IF (NDIAG.GT.0) THEN
          IFLAG = 1
      END IF
      IF (IFLAG.GT.0 .AND. MP.GT.0) WRITE (MP,FMT=9010) IFLAG
      IF (ITYPE.EQ.1) THEN
          CALL MC49AD(-1,N,N,KZ,IRN,JCN,YESA,1,A,N+1,ICPTR,N+1,
     +               IW,JFLAG)
          DO 50 I = 1,N
              IW(I) = 0
   50     CONTINUE
      END IF
      CALL MC34AD(N,IRN,ICPTR,YESA,A,IW)
      DO 60 I = 1,N
          IW(I) = 0
   60 CONTINUE
      CALL MC40BD(N,NNZ,IRN,ICPTR,IPERM,IW,IPROF)
      RETURN
 9000 FORMAT (/,3X,'ERROR MESSAGE: IFLAG = ',I2)
 9010 FORMAT (/,3X,'WARNING MESSAGE : IFLAG =',I2)
      END
      SUBROUTINE MC40BD(N,NNZ,IRN,ICPTR,IPERM,IW,IPROF)
      INTEGER N,NNZ
      INTEGER ICPTR(N+1),IPERM(N),IRN(2*NNZ),IW(3*N+2)
      REAL(8)::IPROF(2)
      INTEGER DEGREE,I,LSTNUM,N1,N2,N3,NODES,NSTOP,NSTRT
      EXTERNAL MC40CD,MC40DD,MC40ED,MC40FD
      DO 10 I = 1,N
          IPERM(I) = 1
   10 CONTINUE
      LSTNUM = 0
      DO 20 I = 1,N
          DEGREE = ICPTR(I+1) - ICPTR(I)
          IF (DEGREE.EQ.0) THEN
              LSTNUM = LSTNUM + 1
              IPERM(I) = -LSTNUM
          END IF
   20 CONTINUE
   30 IF (LSTNUM.LT.N) THEN
          N1 = 1
          N2 = N1 + N - LSTNUM
          N3 = N2 + N + 1
          CALL MC40CD(N,NNZ,IRN,ICPTR,IPERM,IW(N1),IW(N2),IW(N2),IW(N3),
     +               NSTRT,NSTOP,NODES)
          N2 = N1 + NODES
          N3 = N2 + NODES + 1
          CALL MC40DD(N,NNZ,NSTOP,NODES,IRN,ICPTR,IPERM,IW(N1),IW(N2))
          N3 = N2 + N
          CALL MC40ED(N,NNZ,NODES,NSTRT,LSTNUM,IRN,ICPTR,IPERM,IW(N1),
     +               IW(N1),IW(N2))
          GO TO 30
      END IF
      DO 40 I = 1,N
          IPERM(I) = -IPERM(I)
   40 CONTINUE
      CALL MC40FD(N,NNZ,IPERM,IRN,ICPTR,IPROF)
      IF (IPROF(2).GE.IPROF(1)) THEN
          DO 50 I = 1,N
              IPERM(I) = I
   50     CONTINUE
          IPROF(2) = IPROF(1)
      END IF
      RETURN
      END
      SUBROUTINE MC40CD(N,NNZ,IRN,ICPTR,MASK,LS,XLS,DEG,LIST,NSTRT,
     +                 NSTOP,NODES)
      INTEGER N,NNZ,NODES,NSTOP,NSTRT
      INTEGER DEG(N),ICPTR(N+1),IRN(2*NNZ),LIST(N),LS(N),MASK(N),
     +        XLS(N+1)
      INTEGER DEGREE,I,ID1,ID2,ISTOP,ISTRT,LSIZE,LWIDTH,MAXDEP,MINDEG,
     +        MINWID,NLSIZE,NLVL,NODE
      EXTERNAL MC40GD,MC40HD
      MINDEG = N
      DO 10 I = 1,N
          IF (MASK(I).EQ.1) THEN
              DEGREE = ICPTR(I+1) - ICPTR(I)
              IF (DEGREE.LT.MINDEG) THEN
                  NSTRT = I
                  MINDEG = DEGREE
              END IF
          END IF
   10 CONTINUE
      CALL MC40GD(NSTRT,N,N,NNZ,IRN,ICPTR,MASK,LS,XLS,MAXDEP,LWIDTH)
      NODES = XLS(MAXDEP+1) - 1
   20 CONTINUE
      ISTRT = XLS(MAXDEP)
      ISTOP = XLS(MAXDEP+1) - 1
      DO 30 I = ISTRT,ISTOP
          NODE = LS(I)
          DEG(NODE) = ICPTR(NODE+1) - ICPTR(NODE)
   30 CONTINUE
      LSIZE = ISTOP - ISTRT + 1
      IF (LSIZE.GT.1) CALL MC40HD(N,LSIZE,LS(ISTRT),DEG)
C*******
      ISTRT = ISTRT - 1
      DO 40 I = 1,LSIZE
          LIST(I) = LS(ISTRT+I)
   40 CONTINUE
C*******
      ID1 = DEG(LIST(1))
      NLSIZE = 1
      DO 50 I = 2,LSIZE
          ID2 = DEG(LIST(I))
          IF (ID2.NE.ID1) THEN
              NLSIZE = NLSIZE + 1
              LIST(NLSIZE) = LIST(I)
          END IF
          ID1 = ID2
   50 CONTINUE
C**********
      MINWID = NODES
      DO 60 I = 1,NLSIZE
          NODE = LIST(I)
          CALL MC40GD(NODE,MINWID,N,NNZ,IRN,ICPTR,MASK,LS,XLS,NLVL,
     +               LWIDTH)
          IF (NLVL.GT.MAXDEP .AND. LWIDTH.LT.MINWID) THEN
              NSTRT = NODE
              MAXDEP = NLVL
              GO TO 20
          ELSE IF (LWIDTH.LT.MINWID) THEN
              NSTOP = NODE
              MINWID = LWIDTH
          END IF
   60 CONTINUE
      RETURN
      END
      SUBROUTINE MC40DD(N,NNZ,ROOT,NODES,IRN,ICPTR,MASK,LS,XLS)
      INTEGER N,NNZ,NODES,ROOT
      INTEGER ICPTR(N+1),IRN(2*NNZ),LS(NODES),MASK(N),XLS(NODES+1)
      INTEGER I,J,JSTOP,JSTRT,LWIDTH,NLVL
      EXTERNAL MC40GD
      CALL MC40GD(ROOT,NODES,N,NNZ,IRN,ICPTR,MASK,LS,XLS,NLVL,LWIDTH)
      DO 20 I = 1,NLVL
          JSTRT = XLS(I)
          JSTOP = XLS(I+1) - 1
          DO 10 J = JSTRT,JSTOP
              MASK(LS(J)) = I - 1
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE MC40ED(N,NNZ,NODES,NSTRT,LSTNUM,IRN,ICPTR,STATUS,NLIST,
     +                 QUEUE,PRIOR)
      INTEGER W1,W2
      PARAMETER (W1=2,W2=1)
      INTEGER LSTNUM,N,NNZ,NODES,NSTRT
      INTEGER ICPTR(N+1),IRN(2*NNZ),NLIST(NODES),PRIOR(N),QUEUE(NODES),
     +        STATUS(N)
      INTEGER ADDRES,DEGREE,I,ISTOP,ISTRT,J,JSTOP,JSTRT,MAXDEG,MAXPRT,
     +        NABOR,NBR,NEXT,NODE,NQ,PRTY,STANAB
      MAXDEG = NODES
      DO 10 I = 1,NODES
          NODE = NLIST(I)
          DEGREE = ICPTR(NODE+1) - ICPTR(NODE) + 1
          PRIOR(NODE) = W1* (MAXDEG-DEGREE) + W2*STATUS(NODE)
          STATUS(NODE) = 2
   10 CONTINUE
      NQ = 1
      QUEUE(NQ) = NSTRT
      STATUS(NSTRT) = 1
   20 IF (NQ.GT.0) THEN
          MAXPRT = -1
          DO 30 I = 1,NQ
              PRTY = PRIOR(QUEUE(I))
              IF (PRTY.GT.MAXPRT) THEN
                  ADDRES = I
                  MAXPRT = PRTY
              END IF
   30     CONTINUE
          NEXT = QUEUE(ADDRES)
          QUEUE(ADDRES) = QUEUE(NQ)
          NQ = NQ - 1
          ISTRT = ICPTR(NEXT)
          ISTOP = ICPTR(NEXT+1) - 1
          IF (STATUS(NEXT).EQ.1) THEN
              DO 40 I = ISTRT,ISTOP
                  NBR = IRN(I)
                  PRIOR(NBR) = PRIOR(NBR) + W1
                  IF (STATUS(NBR).EQ.2) THEN
                      NQ = NQ + 1
                      QUEUE(NQ) = NBR
                      STATUS(NBR) = 1
                  END IF
   40         CONTINUE
          END IF
          LSTNUM = LSTNUM + 1
          STATUS(NEXT) = -LSTNUM
          DO 60 I = ISTRT,ISTOP
              NBR = IRN(I)
              IF (STATUS(NBR).EQ.1) THEN
                  PRIOR(NBR) = PRIOR(NBR) + W1
                  STATUS(NBR) = 0
                  JSTRT = ICPTR(NBR)
                  JSTOP = ICPTR(NBR+1) - 1
                  DO 50 J = JSTRT,JSTOP
                      NABOR = IRN(J)
                      STANAB = STATUS(NABOR)
                      IF (STANAB.GE.0) THEN
                          PRIOR(NABOR) = PRIOR(NABOR) + W1
                          IF (STANAB.EQ.2) THEN
                              NQ = NQ + 1
                              QUEUE(NQ) = NABOR
                              STATUS(NABOR) = 1
                          END IF
                      END IF
   50             CONTINUE
              END IF
   60     CONTINUE
          GO TO 20
      END IF
      RETURN
      END
      SUBROUTINE MC40FD(N,NNZ,IPERM,IRN,ICPTR,IPROF)
      INTEGER N,NNZ
      INTEGER ICPTR(N+1),IPERM(N),IRN(2*NNZ)
      REAL(8)::IPROF(2)
      INTEGER I,J,JSTOP,JSTRT,K1,NEWMIN,OLDMIN
      INTRINSIC DIM,MIN
      DO 10 I = 1,2
          IPROF(I) = 0
   10 CONTINUE
      DO 30 I = 1,N
          JSTRT = ICPTR(I)
          JSTOP = ICPTR(I+1) - 1
          OLDMIN = I
          NEWMIN = IPERM(I)
          DO 20 J = JSTRT,JSTOP
              K1 = IRN(J)
              OLDMIN = MIN(OLDMIN,K1)
              NEWMIN = MIN(NEWMIN,IPERM(K1))
   20     CONTINUE
          IPROF(1) = IPROF(1) + DIM(I,OLDMIN)
          IPROF(2) = IPROF(2) + DIM(IPERM(I),NEWMIN)
   30 CONTINUE
      DO 40 I = 1,2
          IPROF(I) = IPROF(I) + N
   40 CONTINUE
      RETURN
      END
      SUBROUTINE MC40GD(ROOT,MAXWID,N,NNZ,IRN,ICPTR,MASK,LS,XLS,NLVL,
     +                 LWIDTH)
      INTEGER LWIDTH,MAXWID,N,NLVL,NNZ,ROOT
      INTEGER ICPTR(N+1),IRN(2*NNZ),LS(N),MASK(N),XLS(N+1)
      INTEGER I,J,JSTOP,JSTRT,LBEGIN,LNBR,LVLEND,LVSIZE,NBR,NODE
      INTRINSIC MAX
      MASK(ROOT) = 0
      LS(1) = ROOT
      NLVL = 0
      LVLEND = 0
      LNBR = 1
      LWIDTH = 1
      LVSIZE = 1
   10 IF (LVSIZE.GT.0) THEN
          LWIDTH = MAX(LVSIZE,LWIDTH)
          IF (LWIDTH.GE.MAXWID) THEN
              GO TO 40
          END IF
          LBEGIN = LVLEND + 1
          LVLEND = LNBR
          NLVL = NLVL + 1
          XLS(NLVL) = LBEGIN
          DO 30 I = LBEGIN,LVLEND
              NODE = LS(I)
              JSTRT = ICPTR(NODE)
              JSTOP = ICPTR(NODE+1) - 1
              DO 20 J = JSTRT,JSTOP
                  NBR = IRN(J)
                  IF (MASK(NBR).EQ.1) THEN
                      LNBR = LNBR + 1
                      LS(LNBR) = NBR
                      MASK(NBR) = 0
                  END IF
   20         CONTINUE
   30     CONTINUE
          LVSIZE = LNBR - LVLEND
          GO TO 10
      END IF
      XLS(NLVL+1) = LVLEND + 1
   40 CONTINUE
      DO 50 I = 1,LNBR
          MASK(LS(I)) = 1
   50 CONTINUE
      RETURN
      END
      SUBROUTINE MC40HD(N,NL,LIST,KEY)
      INTEGER N,NL
      INTEGER KEY(N),LIST(NL)
      INTEGER I,J,K,VALUE
      J = 1
      K = LIST(1)
      VALUE = KEY(K)
      DO 10 I = 2,NL
          IF (KEY(LIST(I)).LT.VALUE) THEN
              J = I
              VALUE = KEY(LIST(I))
          END IF
   10 CONTINUE
      LIST(1) = LIST(J)
      LIST(J) = K
      DO 30 I = 2,NL
          J = I
          K = LIST(I)
          VALUE = KEY(K)
   20     IF (VALUE.LT.KEY(LIST(J-1))) THEN
              LIST(J) = LIST(J-1)
              J = J - 1
              GO TO 20
          END IF
          LIST(J) = K
   30 CONTINUE
      END
      BLOCK DATA MC40KD
      INTEGER LP,MP
      COMMON /MC40ID/LP,MP
      SAVE /MC40ID/
      DATA LP/6/,MP/6/
      END

