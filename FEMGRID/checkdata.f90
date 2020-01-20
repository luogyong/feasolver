

   subroutine checkdata()
      use meshds
	  use ds_t
	  use dflib
	  implicit none
      LOGICAL  statusmode,tof2,tof1
      real(8)::xmin1,ymin1,xmax1,ymax1
	  common /xy/ xmin1,ymin1,xmax1,ymax1
	  real(8)::t1,t2,size
	  integer(4)    iunit, ievent, ikeystate, ixpos, iypos
	  !character(1) butter
	  !common butter
	 ! TYPE (windowconfig) thescreen
	  type(wxycoord) wxy
	  type(qwinfo) winfo
      !COMMON              thescreen
	  EXTERNAL            checkgraph,showns,shownn
	  integer(4)          oldcolor,event,ix,iy,res,result
	  !common i4
	  INTEGER(2)          status
      integer(4)   keystate
	  real(8)::scale
	
	
	!open(unit=20,file='user',title='madluo')
	i4=-1
	thescreen.numypixels   = -1
	thescreen.numxpixels   = -1
    thescreen.numtextcols  = -1
    thescreen.numtextrows  = -1
    thescreen.numcolors    = -1
	thescreen.fontsize = -1
   ! thescreen.title = "luo "C
	
	xmin1=winxul
	ymin1=winyul
	xmax1=winxlr
    ymax1=winylr
	winfo%TYPE = QWIN$MAX
    result =SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	open(90,file='user',title='checkdata')
    statusmode = SETWINDOWCONFIG(thescreen)
    iF(.NOT. statusmode) statusmode = SETWINDOWCONFIG(thescreen)
    statusmode = GETWINDOWCONFIG( thescreen )

	result = INSERTMENUQQ (6, 0, $MENUENABLED, 'Check'c, nul)
	result = INSERTMENUQQ (6, 1, $MENUENABLED, 'Point_Num'c, Shownn)
	result = INSERTMENUQQ (6, 2, $MENUENABLED, 'Point_Size'c, Showns)		
	oldcolor=setbkcolorrgb(#00)
	oldcolor=setcolorrgb(#0000ff)
    CALL CLEARSCREEN( $GCLEARSCREEN )

	!CALL SETVIEWPORT( INT2(10), INT2(10), int2(940), int2(590) )
    thescreen.numxpixels=thescreen.numxpixels*0.98
	thescreen.numypixels=thescreen.numypixels*0.88

	CALL SETVIEWPORT( INT2(10), INT2(10), int2(thescreen.numxpixels), int2(thescreen.numypixels) )
    status=RECTANGLE( $GBORDER,INT2(0),INT2(0), int2(thescreen.numxpixels-10), int2(thescreen.numypixels-10))
	status=RECTANGLE( $GBORDER,INT2(1), INT2(1),int2(thescreen.numxpixels-11), int2(thescreen.numypixels-11) )
	
	CALL SETVIEWPORT( INT2(15), INT2(15), int2(thescreen.numxpixels-5), int2(thescreen.numypixels-5) )
	CALL CLEARSCREEN($GVIEWPORT)

	if(ymax1-ymin1>xmax1-xmin1) then	
	    t1=xmax1
		ymax1=ymax1+(ymax1-ymin1)/20
        ymin1=ymin1-(ymax1-ymin1)/20
		xmax1=(ymax1-ymin1)*(thescreen.numxpixels-31)/(thescreen.numypixels-31)+xmin1
		xmin1=xmin1-(xmax1-t1)/2
		xmax1=xmax1-(xmax1-t1)/2
    else
	    t1=ymax1
		xmax1=xmax1+(xmax1-xmin1)/20
        xmin1=xmin1-(xmax1-xmin1)/20		
		ymax1=(xmax1-xmin1)*(thescreen.numypixels-31)/(thescreen.numxpixels-31)+ymin1
		ymin1=ymin1-(ymax1-t1)/2
		ymax1=ymax1-(ymax1-t1)/2
	end if
	
	status=SETWINDOW( .true., xmin1,ymin1,xmax1,ymax1)
	oldcolor=setcolorrgb(#ff00)
	status=RECTANGLE_w( $GBORDER,xmin1,ymin1,xmax1,ymax1)


	call checkgraph(iunit, ievent, ikeystate, ixpos, iypos)

	event = MOUSE$LBUTTONDOWN
	event=  ior(event,MOUSE$LBUTTONUP)
    event = IOR (event, MOUSE$RBUTTONDOWN)
    event = IOR (event, MOUSE$LBUTTONDBLCLK)
    event = IOR (event, MOUSE$RBUTTONDBLCLK)
    event = IOR (event, MOUSE$MOVE)
	i4 = registermouseevent(90, event, checkgraph)

	result =  SETWSIZEQQ(90, winfo)

	do while(.true.)

	 !res = waitonmouseevent(event, keystate, ix, iy)
	 !if((MOUSE$KS_SHIFT .AND. res) == MOUSE$KS_SHIFT)   return
	 !i4 = registermouseevent(90, event, checkgraph)
	 !res = waitonmouseevent(MOUSE$move, keystate, ix, iy)

	 butter = GETCHARQQ()
	 tof1=(ichar(butter)==ichar('q').OR.ichar(butter)==ichar('Q').OR.ichar(butter)==27)
	 tof2=(ichar(butter)==ichar('c').OR.ichar(butter)==ichar('C').OR.ichar(butter)==13)

	 if(tof1) stop

	 if(tof2) then
	   !winfo%TYPE = QWIN$MIN
       !result =  SETWSIZEQQ(90, winfo)
	   close(90)
	   exit
	 end  if
	

	end do
	
	result = DELETEMENUQQ(6,0)	

	

   end subroutine

   subroutine checkgraph(iunit, ievent, ikeystate, ixpos, iypos)
      use meshDS
	  use ds_t
      USE DFLIB
	  implicit none
      integer(4)    oldcolor,iunit, ievent, ikeystate, ixpos, iypos
	  INTEGER(2)    status,result
      INTEGER(4)     res
	  character(256)  msg
	  character(256) nnum,nnum1
	  !character(1) butter
	  !common butter
	  integer::i,j,k,n1,n2
      logical::tof,tof1,tof2
      real(8)::x1,y1,x2,y2,size,scale
	  real(8)::xmin1,ymin1,xmax1,ymax1
	  common /xy/ xmin1,ymin1,xmax1,ymax1
	  !TYPE (windowconfig):: thescreen
	  !common thescreen
	
	  TYPE (wxycoord)       wxy
	  scale=1

	  if(allocated(b)) deallocate(b)
	  allocate(b(inpn)) !利用b去判断该点是否已经输出
	  b=0
	  if(i4>=0) then
         call getwindowcoord(ixpos-15, iypos-15,wxy)
		
		 if(ievent==MOUSE$MOVE) then
		    write(msg,*) 'Put C to close window and continue meshgen. Put Q to terminate program.'
            CALL SETMESSAGEQQ (msg, QWIN$MSG_INPUTPEND)
		    return
	     end if

		 if((MOUSE$KS_lbutton .AND. ikeystate) == MOUSE$KS_lbutton) scale=scale/2
	          	
		 if((MOUSE$KS_Rbutton .AND. ikeystate) == MOUSE$KS_Rbutton) scale=scale*2


		 xmin1=wxy.wx-(wxy.wx-xmin1)*scale
		 ymin1=wxy.wy-(wxy.wy-ymin1)*scale
		 xmax1=wxy.wx+(xmax1-wxy.wx)*scale
		 ymax1=wxy.wy+(ymax1-wxy.wy)*scale
         CALL SETVIEWPORT( INT2(15), INT2(15), int2(thescreen.numxpixels-5), int2(thescreen.numypixels-5) )
	     CALL CLEARSCREEN($GVIEWPORT)
		 oldcolor=setcolorrgb(#0000ff)
	     status=SETWINDOW( .true., xmin1,ymin1,xmax1,ymax1)
	     oldcolor=setcolorrgb(#ff00)
	     status=RECTANGLE_w( $GBORDER,xmin1,ymin1,xmax1,ymax1)
	
	  end if
  	
	  result = INITIALIZEFONTS()
      result = SETFONT('t''Arial''h15w8p')

	 bp=>bnhead
	 do while(.true..and.associated(bp))
		call moveto_w(bp.npt.x,bp.npt.y,wxy)
		oldcolor=setcolorrgb(#00ff00)
		status=lineto_w(bp.next.npt.x,bp.next.npt.y)

		do i=1,inpn
		   if(((arr_t(i).x-bp.next.npt.x)**2+(arr_t(i).y-bp.next.npt.y)**2)<1e-10) exit
		end do
		if(b(i)==0) then
		   call  value(nnum,i)
		   oldcolor=setcolorrgb(#ffff00)
		   call outgtext(trim(adjustL(nnum)))
		   b(i)=1
		end if
      	if(associated(bp.next,bnhead)) exit
		bp=>bp.next
	 end do

	  !各子域的从属子域的个数之和
	  n1=0

	  do i=1,cln-ccln-n1
      
	     do j=2,csl(i).num
		    call moveto_w(csl(i).conpoint(1,j-1),csl(i).conpoint(2,j-1),wxy)
		    !oldcolor=setcolorrgb(#00ffff)
        select case(mod(i,15))
		    case(1)
                oldcolor=setcolorrgb($HIRED)
		    case(2)
                oldcolor=setcolorrgb($HIGREEN)
		    case(3)
                oldcolor=setcolorrgb($HIBLUE)
		    case(4)
                oldcolor=setcolorrgb($HIYELLOW)
		    case(5)
                oldcolor=setcolorrgb($HICYAN)
		    case(6)
                oldcolor=setcolorrgb($HIMAGENTA)
		    case(7)
                oldcolor=setcolorrgb($HIWHITE)
		    case(8)
			    oldcolor=setcolorrgb($Logray)
		    case(9)
                oldcolor=setcolorrgb($LORED)
		    case(10)
                oldcolor=setcolorrgb($LOGREEN)
		    case(11)
                oldcolor=setcolorrgb($LOBLUE)
		    case(12)
                oldcolor=setcolorrgb($LOBROWN)
		    case(13)
                oldcolor=setcolorrgb($LOCYAN)
		    case(14)
                oldcolor=setcolorrgb($LOMAGENTA)
		    case(0)
                oldcolor=setcolorrgb($Lowhite)			
    !		   case default
    !              oldcolor=setcolorrgb(#fff00)
		end select                 
		     
            
		x1=csl(i).conpoint(1,j)
		y1=csl(i).conpoint(2,j)
		status=lineto_w(x1,y1)
		if(csl(i).flag==1) then
		    call moveto_w(csl(i).conpoint(1,1),csl(i).conpoint(2,1),wxy)
            !oldcolor=setcolorrgb(#00ffff)
		    x1=csl(i).conpoint(1,csl(i).num)
		    y1=csl(i).conpoint(2,csl(i).num)
		    status=lineto_w(x1,y1)
		end if
        
        do k=1,inpn
		       if(((arr_t(k).x-csl(i).conpoint(1,j))**2+(arr_t(k).y-csl(i).conpoint(2,j))**2)<1e-10) exit
		    end do
			if(b(k)==0) then
			   call  value(nnum,k)
			   oldcolor=setcolorrgb(#ffff00)
			   call outgtext(trim(adjustL(nnum)))
			   b(k)=1
			 end if
		 end do
		
		 do k=1,inpn
		    if(((arr_t(k).x-csl(i).conpoint(1,1))**2+(arr_t(k).y-csl(i).conpoint(2,1))**2)<1e-10) exit
		 end do
		
		 if(b(k)==0) then
		    call  value(nnum,k)
		    oldcolor=setcolorrgb(#ffff00)
            call moveto_w(csl(i).conpoint(1,1),csl(i).conpoint(2,1),wxy)
		    call outgtext(trim(adjustL(nnum)))
			b(k)=1
		 end if
		

	  end do


	  do i=1,aupn
	  	 size=aupoint(i).s/10
		 oldcolor=setcolorrgb(#ff00f0)
!		 do k=1,inpn
!		    if(arr_t(k).num==aupoint(i).point) exit
!		 end do
		 k=aupoint(i).point
		 status=ellipse_W($GBORDER,arr_t(k).x-size,arr_t(k).y-size,arr_t(k).x+size,arr_t(k).y+size)
         if(b(k)==0) then
			call  value(nnum,k)
			oldcolor=setcolorrgb(#ffff00)
			call moveto_w(arr_t(k).x,arr_t(k).y,wxy)
			call outgtext(trim(adjustL(nnum)))
			b(k)=1		
		 end if	
	  end do

	  do i=1,ccln
	     size=ccl(i).r
         oldcolor=setcolorrgb(#fff0f0)
!		 do k=1,inpn
!		    if(arr_t(k).num==ccl(i).point) exit
!		 end do
		 k=ccl(i).point
         status=ellipse_W($GBORDER,arr_t(k).x-size,arr_t(k).y-size,arr_t(k).x+size,arr_t(k).y+size)
         if(b(k)==0) then
			call  value(nnum,k)
			oldcolor=setcolorrgb(#ffff00)
			call moveto_w(arr_t(k).x,arr_t(k).y,wxy)
			call outgtext(trim(adjustL(nnum)))
			b(k)=1		
		 end if	
	  end do

	  do i=1,inpn
	     size=arr_t(i).s/20.0
         oldcolor=setcolorrgb(#ffff00)
		 k=i
         status=ellipse_W($GBORDER,arr_t(k).x-size,arr_t(k).y-size,arr_t(k).x+size,arr_t(k).y+size)
         if(b(k)==0) then
			call  value(nnum,k)
			oldcolor=setcolorrgb(#ffff00)
			call moveto_w(arr_t(k).x,arr_t(k).y,wxy)
			call outgtext(trim(adjustL(nnum)))
			b(k)=1		
		 end if	
	  end do


   end subroutine

	subroutine showns()
		use ds_t
		implicit none
		showvalue2=2
	end subroutine

	subroutine shownn()
		use ds_t
		implicit none
		showvalue2=1
	end subroutine

	subroutine value(nnum,n1)
		use ds_t
		implicit none
		integer::n1
		character(8)::nnum
		select case(showvalue2)
		   case(1)
		   write(nnum,'(i6)') arr_t(n1).num
		   case(2)
		   write(nnum,'(f8.3)') arr_t(n1).s
		end select
	end subroutine