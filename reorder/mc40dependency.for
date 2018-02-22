* COPYRIGHT (c) 1993 AEA Technology
*######DATE 10 Feb 1993
C       Toolpack tool decs employed.
C
      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
C THIS SUBROUTINE ACCEPTS AS INPUT THE STANDARD DATA STRUCTURE FOR
C     A SYMMETRIC MATRIX STORED AS A LOWER TRIANGLE AND PRODUCES
C     AS OUTPUT THE SYMMETRIC MATRIX HELD IN THE SAME DATA
C     STRUCTURE AS A GENERAL MATRIX.
C N IS AN INTEGER VARIABLE THAT MUST BE SET BY THE USER TO THE
C     ORDER OF THE MATRIX. NOT ALTERED BY THE ROUTINE
C     RESTRICTION (IBM VERSION ONLY): N LE 32767.
C IRN IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY THAT
C     MUST BE SET BY THE USER TO HOLD THE ROW INDICES OF THE LOWER
C     TRIANGULAR PART OF THE SYMMETRIC MATRIX.  THE ENTRIES OF A
C     SINGLE COLUMN ARE CONTIGUOUS. THE ENTRIES OF COLUMN J
C     PRECEDE THOSE OF COLUMN J+1 (J_=_1, ..., N-1), AND THERE IS
C     NO WASTED SPACE BETWEEN COLUMNS. ROW INDICES WITHIN A COLUMN
C     MAY BE IN ANY ORDER.  ON EXIT IT WILL HAVE THE SAME MEANING
C     BUT WILL BE CHANGED TO HOLD THE ROW INDICES OF ENTRIES IN
C     THE EXPANDED STRUCTURE.  DIAGONAL ENTRIES NEED NOT BE
C     PRESENT. THE NEW ROW INDICES ADDED IN THE UPPER TRIANGULAR
C     PART WILL BE IN ORDER FOR EACH COLUMN AND WILL PRECEDE THE
C     ROW INDICES FOR THE LOWER TRIANGULAR PART WHICH WILL REMAIN
C     IN THE INPUT ORDER.
C JCOLST IS AN INTEGER ARRAY OF LENGTH N+1 THAT MUST BE SET BY
C     THE USER SO THAT JCOLST(J) IS THE POSITION IN ARRAYS IRN AND
C     A OF THE FIRST ENTRY IN COLUMN J (J_=_1, ..., N).
C     JCOLST(N+1) MUST BE SET TO ONE MORE THAN THE TOTAL NUMBER OF
C     ENTRIES.  ON EXIT, JCOLST(J) WILL HAVE THE SAME MEANING BUT
C     WILL BE CHANGED TO POINT TO THE POSITION OF THE FIRST ENTRY
C     OF COLUMN J IN THE EXPANDED STRUCTURE. THE NEW VALUE OF
C     JCOLST(N+1) WILL BE ONE GREATER THAN THE NUMBER OF ENTRIES
C     IN THE EXPANDED STRUCTURE.
C YESA IS A LOGICAL VARIABLE THAT MUST BE SET TO .TRUE. IF THE
C     USER DESIRES TO GENERATE THE EXPANDED FORM FOR THE VALUES ALSO.
C     IF YESA IS .FALSE., THE ARRAY A WILL NOT BE REFERENCED.  IT IS
C     NOT ALTERED BY THE ROUTINE.
C A IS A REAL (DOUBLE PRECISION IN THE D VERSION) ARRAY THAT
C     CAN BE SET BY THE USER SO THAT A(K) HOLDS THE VALUE OF THE
C     ENTRY IN POSITION K OF IRN, {K = 1, _..._ JCOLST(N+1)-1}.
C     ON EXIT, IF YESA IS .TRUE., THE ARRAY WILL HOLD THE VALUES
C     OF THE ENTRIES IN THE EXPANDED STRUCTURE CORRESPONDING TO
C     THE OUTPUT VALUES OF IRN.   IF YESA IS .FALSE., THE ARRAY IS
C     NOT ACCESSED BY THE SUBROUTINE.
C IW IS AN INTEGER (INTEGER*2 IN IBM VERSION) ARRAY OF LENGTH
C     N THAT WILL BE USED AS WORKSPACE.
C
C CKP1 IS A LOCAL VARIABLE USED AS A RUNNING POINTER.
C OLDTAU IS NUMBER OF ENTRIES IN SYMMETRIC STORAGE.
C     .. Scalar Arguments ..
      INTEGER N
      LOGICAL YESA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
C     ..
C     .. Local Scalars ..
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
C     ..
C     .. Executable Statements ..
C
      OLDTAU = JCOLST(N+1) - 1
C INITIALIZE WORK ARRAY
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
C
C IW(J) IS SET EQUAL TO THE TOTAL NUMBER OF ENTRIES IN COLUMN J
C     OF THE EXPANDED SYMMETRIC MATRIX.
C NDIAG COUNTS NUMBER OF DIAGONAL ENTRIES PRESENT
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1

          ELSE
            NDIAG = NDIAG + 1
          END IF

   10   CONTINUE
   20 CONTINUE
C
C NEWTAU IS NUMBER OF ENTRIES IN EXPANDED STORAGE.
      NEWTAU = 2*OLDTAU - NDIAG
C IPKP1 POINTS TO POSITION AFTER END OF COLUMN BEING CURRENTLY
C     PROCESSED
      IPKP1 = OLDTAU + 1
C CKP1 POINTS TO POSITION AFTER END OF SAME COLUMN IN EXPANDED
C     STRUCTURE
      CKP1 = NEWTAU + 1
C GO THROUGH THE ARRAY IN THE REVERSE ORDER PLACING LOWER TRIANGULAR
C     ELEMENTS IN THE APPROPRIATE SLOTS.
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
C LENK IS NUMBER OF ENTRIES IN COLUMN J OF ORIGINAL STRUCTURE
        LENK = I2 - I1
C JSTART IS RUNNING POINTER TO POSITION IN NEW STRUCTURE
        JSTART = CKP1
C SET IKP1 FOR NEXT COLUMN
        IPKP1 = I1
        I2 = I2 - 1
C RUN THROUGH COLUMNS IN REVERSE ORDER
C LOWER TRIANGULAR PART OF COLUMN MOVED TO END OF SAME COLUMN IN
C     EXPANDED FORM
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
C JCOLST IS SET TO POSITION OF FIRST ENTRY IN LOWER TRIANGULAR PART OF
C     COLUMN J IN EXPANDED FORM
        JCOLST(J) = JSTART
C SET CKP1 FOR NEXT COLUMN
        CKP1 = CKP1 - IW(J)
C RESET IW(J) TO NUMBER OF ENTRIES IN LOWER TRIANGLE OF COLUMN.
        IW(J) = LENK
   40 CONTINUE
C
C AGAIN SWEEP THROUGH THE COLUMNS IN THE REVERSE ORDER, THIS
C     TIME WHEN ONE IS HANDLING COLUMN J THE UPPER TRIANGULAR
C     ELEMENTS A(J,I) ARE PUT IN POSITION.
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
C RUN DOWN COLUMN IN ORDER
C NOTE THAT I IS ALWAYS GREATER THAN OR EQUAL TO J
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN

      END
* *******************************************************************
* COPYRIGHT (c) 1967 Hyprotech UK
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
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       SAVE statement for COMMON FA01ED added.
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C
C
      DOUBLE PRECISION FUNCTION FA01AD(I)
      INTEGER I
      DOUBLE PRECISION R,S
      INTRINSIC DINT,MOD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      EXTERNAL FA01FD
      SAVE /FA01ED/
      R = GR*9228907D0/65536D0
      S = DINT(R)
      GL = MOD(S+GL*9228907D0,65536D0)
      GR = R - S
      IF (I.GE.0) FA01AD = (GL+GR)/65536D0
      IF (I.LT.0) FA01AD = (GL+GR)/32768D0 - 1.D0
      GR = GR*65536D0
      RETURN
      END
      SUBROUTINE FA01BD(MAX,NRAND)
      INTEGER MAX,NRAND
      DOUBLE PRECISION FA01AD
      EXTERNAL FA01AD
      INTRINSIC DBLE,INT
      NRAND = INT(FA01AD(1)*DBLE(MAX)) + 1
      RETURN
      END
      SUBROUTINE FA01CD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      IL = GL
      IR = GR
      RETURN
      END
      SUBROUTINE FA01DD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      GL = IL
      GR = IR
      RETURN
      END
      BLOCK DATA FA01FD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      DATA GL/21845D0/
      DATA GR/21845D0/
      END

