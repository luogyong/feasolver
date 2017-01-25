! Copyright (C) 2013 Intel Corporation. All Rights Reserved. 
!
! The source code contained or described herein and all documents related to the source code 
! ("Material") are owned by Intel Corporation or its suppliers or licensors. Title to the 
! Material remains with Intel Corporation or its suppliers and licensors.  The Material is 
! protected by worldwide copyright laws and treaty provisions. No part of the Material may be 
! used, copied, reproduced, modified, published, uploaded, posted, transmitted, distributed, 
! or disclosed in any way except as expressly provided in the license provided with the 
! Materials.  No license under any patent, copyright, trade secret or other intellectual 
! property right is granted to or conferred upon you by disclosure or delivery of the 
! Materials, either expressly, by implication, inducement, estoppel or otherwise, except as 
! expressly provided in the license provided with the Materials.

!**********************************************************************
!
! SciGraph -- Scientific Graphs for Intel Visual Fortran
!
! 7/2013 Restructure if statements in SciGetTextLength
! 6/2001 Replace calls to setgtextvector to calls to setgtextrotation
! 6/2001 Add EOL check in loops removing leading blanks 

MODULE SGLOWLVL

USE IFQWIN
USE SGDATA
IMPLICIT NONE

CONTAINS

! This function makes sure we are in the correct graphics mode, and then
!    sets all the graphics library settings to values that are good as
!    defaults

   LOGICAL FUNCTION SciGraphicsMode(graph)
          logical lretval

          RECORD /GraphSettings/ graph
          INTEGER retVal
          RECORD /windowconfig/ wc
	  INTEGER solidi
    INTEGER*1, PARAMETER :: SOLID(8) = (/ (#FF,solidi=1,8) /)   ! solid fill

          RECORD /PrivateGraphSettings/ pgraph   ! private graph
          RECORD /GraphSettings/ cgraph          ! copy of private graph
          EQUIVALENCE(cgraph,pgraph)

          cgraph=graph
          IF (graph.setGraphMode .AND. pgraph.didPlotGraph .EQ. 0) THEN
              lretval = SGSetMaxRes()
              IF (.not. lretVal) THEN
                  SciGraphicsMode=.FALSE.
                  RETURN
              END IF
          END IF

          lretval = GETWINDOWCONFIG(wc)

          IF (graph.x1 < 0  .OR.  graph.y1 < 0  &
                 .OR.  graph.x2 > wc.numxpixels   &
                 .OR.  graph.y2 > wc.numypixels) THEN
              SciGraphicsMode=.FALSE.
              RETURN
          END IF

          CALL SETVIEWPORT(graph.x1,graph.y1,graph.x2,graph.y2)   ! clip region
          CALL CLEARSCREEN($GVIEWPORT)
          CALL SETLINESTYLE($LVDEF)
          CALL SETFILLMASK(SOLID)
          retVal=SETCOLOR(graph.graphBgColor)
          retVal=RECTANGLE($GFILLINTERIOR,graph.x1,graph.y1,  &
                                  graph.x2,graph.y2)              ! "background"
          retVal=SETCOLOR(graph.graphColor)
          retVal=SETWINDOW(LOGICAL(.false.,2),$WXMIN,$WYMIN,$WXMAX,$WYMAX)   ! rescale it

          SciGraphicsMode=.TRUE.
          RETURN
      END FUNCTION

!   Does a MOVETO_W, but allows us to not clutter up the code with lots of temp wxycoords

      SUBROUTINE SciMoveto_w(x,y)
          REAL*8 x,y
          RECORD /wxycoord/ lastp
          CALL MOVETO_W(x,y,lastp)
      END SUBROUTINE

!   Draws a line from one (x,y) to another.  Draws in the specified line type
!     from the library.  Here so we don't have to keep doing MOVETO_W, LINETO_W
!     in the code

      SUBROUTINE SciLine_w(lt,x1,y1,x2,y2)

          INTEGER(2) lt
          REAL*8 x1,y1,x2,y2

          RECORD /wxycoord/ lastp
          INTEGER retv
          LOGICAL thick
          REAL*8 dx,dy,adx,ady

          thick=SciSetLineType(lt)

          CALL MOVETO_W(x1,y1,lastp)
          retv=LINETO_W(x2,y2)

          IF (thick) THEN
              IF (y2 .EQ. y1) THEN
                  CALL MOVETO_W(x1,y1+16,lastp)
                  retv=LINETO_W(x2,y2+16)
                  CALL MOVETO_W(x1,y1-16,lastp)
                  retv=LINETO_W(x2,y2-16)
              ELSE IF (x2 .EQ. x1) THEN
                  CALL MOVETO_W(x1+16,y1,lastp)
                  retv=LINETO_W(x2+16,y2)
                  CALL MOVETO_W(x1-16,y1,lastp)
                  retv=LINETO_W(x2-16,y2)
              ELSE
                  dy=(y2-y1)/16         ! number of pixels rise
                  dx=(x2-x1)/16         !  and run

                  IF (abs(dy) > abs(dx)) THEN
                      adx=SIGN(1.0D0,dy)
                      ady=0
                  ELSE
                      ady=-SIGN(1.0D0,dx)
                      adx=0
                  END IF

                  CALL MOVETO_W(x1+adx*16,y1+ady*16,lastp)
                  retv=LINETO_W(x2+adx*16,y2+ady*16)
                  CALL MOVETO_W(x1-adx*16,y1-ady*16,lastp)
                  retv=LINETO_W(x2-adx*16,y2-ady*16)
              END IF
          END IF

          RETURN
      END SUBROUTINE


!   Draws an circle with a given center and radius.  Draws with the specified line type
!     from the library and should be filled in if desired

      SUBROUTINE SciCircle_w(x,y,r,lt,ft)

          REAL*8 x,y,r
          INTEGER(2) lt
          INTEGER*2 ft

          LOGICAL thick
          REAL*8 x1,y1,x2,y2
          INTEGER retv

          x1=x-r
          y1=y-r
          x2=x+r
          y2=y+r

          thick=SciSetLineType(lt)
          retv=ELLIPSE_W(ft,x1,y1,x2,y2)

          IF (thick) THEN
              retv=ELLIPSE_W(ft,x1-16,y1-16,x2+16,y2+16)
          END IF

          RETURN
      END SUBROUTINE

!   Sets the line type to a library specified value and returns whether it is
!     thick or not

      LOGICAL FUNCTION SciSetLineType(lt)

          INTEGER(2) lt

          LOGICAL thick

          thick=.FALSE.

          SELECT CASE (lt)
          CASE ($LTNONE)
              ! do nothing if line type is NONE
          CASE ($LTSOLID)
              CALL SETLINESTYLE($LVSOLID)
          CASE ($LTDASH)
              CALL SETLINESTYLE($LVDASH)
          CASE ($LTDASHDOT)
              CALL SETLINESTYLE($LVDASHDOT)
          CASE ($LTDASHDOTDOT)
              CALL SETLINESTYLE($LVDASHDOTDOT)
          CASE ($LTDOT)
              CALL SETLINESTYLE($LVDOT)
          CASE ($LTTHICKSOLID)
              CALL SETLINESTYLE($LVSOLID)
              thick=.TRUE.
          CASE ($LTTHICKDASH)
              CALL SETLINESTYLE($LVDASH)
              thick=.TRUE.
          CASE ($LTTHICKDASHDOT)
              CALL SETLINESTYLE($LVDASHDOT)
              thick=.TRUE.
          CASE ($LTTHICKDASHDOTDOT)
              CALL SETLINESTYLE($LVDASHDOTDOT)
              thick=.TRUE.
          CASE ($LTTHICKDOT)
              CALL SETLINESTYLE($LVDOT)
              thick=.TRUE.
          END SELECT

          SciSetLineType=thick
          RETURN
      END FUNCTION



!   Sets the font to one of the library defined values

      INTEGER FUNCTION SciSetFont(fontType)

          INTEGER(2) fontType

          INTEGER retv /-1/

          SELECT CASE(fontType)
          CASE ($FTCOUR)
                retv=SETFONT($FVCOUR)
          CASE ($FTCOURSM)
                retv=SETFONT($FVCOURSM)
          CASE ($FTTROMAN)
                retv=SETFONT($FVTROMAN)
          CASE ($FTTROMANSM)
                retv=SETFONT($FVTROMANSM)
          CASE ($FTSANSSERIF)
                retv=SETFONT($FVSANSSERIF)
          CASE ($FTSANSSERIFSM)
                retv=SETFONT($FVSANSSERIFSM)
          CASE DEFAULT
                retv=SETFONT($FVDEF)
          END SELECT

          SciSetFont=retv
          RETURN
      END FUNCTION


!   Get font heigth
      INTEGER FUNCTION SciGetTextHeight()

          RECORD /fontinfo/ fi
          INTEGER retv,th

          retv=GETFONTINFO(fi)

          ! maybe can only take height for bitmap fonts
          IF (fi.type .NE. 0 .AND. fi.pixheight .LE. 0) THEN
              th=16
          ELSE
              th=fi.pixheight
          END IF

          SciGetTextHeight=th
          RETURN
      END FUNCTION


!   Get text length
      INTEGER FUNCTION SciGetTextLength(text,column)

          CHARACTER*(*) text
          LOGICAL column            ! if text in a column

          RECORD /fontinfo/ fi
          INTEGER retv,tl

 
          IF (column) THEN
              retv=GETFONTINFO(fi)
              tl=fi.pixheight*LEN_TRIM(text)
              SciGetTextLength=tl
              RETURN
          ELSE                                        !normal text
              tl=GETGTEXTEXTENT(text)
              SciGetTextLength=tl
              RETURN
          END IF
 
      END FUNCTION

!   Centers a given string between two x and between two y

      SUBROUTINE SciCenter(xwpr,ywpr,xp1,xp2,yp1,yp2,text)

          REAL(8) xwpr,ywpr
          REAL(8) xp1,xp2,yp1,yp2
          CHARACTER*(*) text

          INTEGER nsend,nsstart
          RECORD /WXYCOORD/ lp

          REAL*8 xsp,ysp
          INTEGER tl,th

          nsstart=1          ! trash leading blanks
          DO WHILE (text(nsstart:nsstart) .EQ. ' ')
              nsstart=nsstart+1
	      if (nsstart > len(text)) exit
          END DO

          nsend=nsstart+len_trim(text(nsstart:))-1
          tl=GETGTEXTEXTENT(text(nsstart:nsend))
          xsp=xp1+(xp2-xp1-tl*xwpr)/2.0D0

          th=SciGetTextHeight()
          ysp=yp1+(yp2-yp1-th*ywpr)/2.0D0 ! text draw starts from upper left

          CALL MOVETO_W(xsp,ysp,lp)
          CALL OUTGTEXT(text(nsstart:nsend))
      END SUBROUTINE


!   Centers a given string starting or ending at an x and between two y

      SUBROUTINE SciEndCenter(xwpr,ywpr,xp1,xEnd,yp1,yp2,text)

          REAL(8) xwpr,ywpr
          REAL*8 xp1,yp1,yp2
          LOGICAL xEnd
          CHARACTER*(*) text

          INTEGER nsend,nsstart
          RECORD /WXYCOORD/ lp

          REAL*8 xsp,ysp
          INTEGER tl,th

          nsstart=1          ! trash leading blanks
          DO WHILE (text(nsstart:nsstart) .EQ. ' ')
              nsstart=nsstart+1
	      if (nsstart > len(text)) exit
          END DO

          nsend=nsstart+len_trim(text(nsstart:))-1
          tl=GETGTEXTEXTENT(text(nsstart:nsend))
          IF (xEnd) THEN
              xsp=xp1-tl*xwpr
          ELSE
              xsp=xp1
          END IF

          th=SciGetTextHeight()
          ysp=yp1+(yp2-yp1-th*ywpr)/2.0D0 ! text draw starts from upper left

          CALL MOVETO_W(xsp,ysp,lp)
          CALL OUTGTEXT(text(nsstart:nsend))
      END SUBROUTINE

!   Centers a given line of rotated text at a given x and between two y

      SUBROUTINE SciRotCenter(xwpr,ywpr,xp1,xp2,yp1,yp2,text)

          REAL(8) xwpr,ywpr
          REAL*8 xp1,xp2,yp1,yp2
          CHARACTER*(*) text

          INTEGER nsstart,nsend
          RECORD /wxycoord/ lp
          REAL*8 ysp,xsp
          INTEGER tl,th
          LOGICAL column

          nsstart=1          ! trash leading blanks
          DO WHILE (text(nsstart:nsstart) .EQ. ' ')
              nsstart=nsstart+1
	      if (nsstart > len(text)) exit
          END DO
          nsend=nsstart+len_trim(text(nsstart:))-1

!          CALL SETGTEXTVECTOR(0,1)        ! rotate 90d C-clockwise
          CALL SETGTEXTROTATION(900)       ! rotate 90d C-clockwise
          column= GRSTATUS() .NE. $GROK

          tl=SciGetTextLength(text(nsstart:nsend),column)
          ysp=yp1+(yp2-yp1+tl*ywpr)/2.0D0

          th=SciGetTextHeight()
          xsp=xp1+(xp2-xp1-th*xwpr)/2.0D0

          CALL MOVETO_W(xsp,ysp,lp)
          IF (column) THEN
              ysp=ysp-tl*ywpr
              CALL MOVETO_W(xsp,ysp,lp)
              CALL SciVertText(text(nsstart:nsend))
          ELSE
              CALL OUTGTEXT(text(nsstart:nsend))
          END IF

!          CALL SETGTEXTVECTOR(1,0)  ! set it back to the default
          CALL SETGTEXTROTATION(0)   ! set it back to the default
      END SUBROUTINE

!   Centers a given line of rotated text at a given x and between two y

      SUBROUTINE SciRotEndCenter(xwpr,ywpr,xp1,xBottom,yp1,yp2,text)

          REAL(8) xwpr,ywpr
          REAL*8 xp1,yp1,yp2
          LOGICAL xBottom
          CHARACTER*(*) text

          INTEGER nsstart,nsend
          RECORD /wxycoord/ lp
          REAL*8 ysp,xsp
          INTEGER tl,th
          LOGICAL column

          nsstart=1          ! trash leading blanks
          DO WHILE (text(nsstart:nsstart) .EQ. ' ')
              nsstart=nsstart+1
	      if (nsstart > len(text)) exit
          END DO
          nsend=nsstart+len_trim(text(nsstart:))-1

!          CALL SETGTEXTVECTOR(0,1)        ! rotate 90d C-clockwise
          CALL SETGTEXTROTATION(900)       ! rotate 90d C-clockwise
          column= GRSTATUS() .NE. $GROK

          tl=SciGetTextLength(text(nsstart:nsend),column)
          ysp=yp1+(yp2-yp1+tl*ywpr)/2.0D0

          th=SciGetTextHeight()
          IF (xBottom) THEN
              xsp=xp1-th*xwpr
          ELSE
              xsp=xp1
          END IF

          CALL MOVETO_W(xsp,ysp,lp)
          IF (column) THEN
              ysp=ysp-tl*ywpr
              CALL MOVETO_W(xsp,ysp,lp)
              CALL SciVertText(text(nsstart:nsend))
          ELSE
              CALL OUTGTEXT(text(nsstart:nsend))
          END IF

!          CALL SETGTEXTVECTOR(1,0)  ! set it back to the default
          CALL SETGTEXTROTATION(0)   ! set it back to the default
      END SUBROUTINE

      Logical Function SGSetMaxRes
      record /windowconfig/ wc
      wc.numxpixels = -1
      wc.numypixels = -1
      wc.numtextrows = -1
      wc.numtextcols = -1
      wc.numcolors = -1
      wc.fontsize = -1
      SGSetMaxRes = setwindowconfig(wc)
      end function

!
!  Vertical graphics output of a text string (from the current position)
!

      SUBROUTINE SciVertText (text)

      CHARACTER*(*)  text
      CHARACTER*1    ch
      INTEGER        status, loop, ipos, jpos, indent
      RECORD         / xycoord / xy, xyorig
      RECORD         / fontinfo / fi

!     Get the beginning location for output
      CALL getcurrentposition(xyorig)
      ipos = xyorig.xcoord
      jpos = xyorig.ycoord

!     The font info is used to center the letters, and determine how
!       far each letter should be below the previous one
      status = getfontinfo(fi)

!     Output one char of the string at a time, in an descending column
      do loop = 1, len(text)
        ch = text(loop:loop)
        indent = ipos + (fi.pixwidth-getgtextextent(ch))/2
        CALL moveto(int2(indent), int2(jpos), xy)
        CALL outgtext(ch)
        jpos = jpos + fi.pixheight
      end do

!     Set new current position to the end of the vertical text column
      CALL moveto(int2(ipos), int2(jpos), xy)
      END SUBROUTINE

END MODULE SGLOWLVL
