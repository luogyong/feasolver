subroutine graph
	USE DFLIB
	use meshDS
	use ds_t
    use geomodel    
	implicit none
	LOGICAL statusmode
	TYPE (windowconfig) myscreen
	type(wxycoord) wxy
	type(qwinfo) winfo
	COMMON  myscreen
	EXTERNAL elementgraph,displayall,displaysubstructure, &
					 displaymainmesh,meshoutput,UpdateMaterial, &
					 tri_element,contact_element,all_element, &
					 cut_off_wall_element,EQ,outplot_checkdata,sple_element, &
					 check_BcType,check_BcValue,check_number,check_subnum, &
					 check_material,check_null,check_bandwidth,check_kcd, &
					 check_enum,check_EZONE,penaltyelementshow, &
					 DISCONTINUITIES,offset,gen_6_noded_triangle, &
					 gen_15_noded_triangle,Generate_3D_MODEL, &
					 Check_Edge,gen_dem_particle					 
	integer(4)::oldcolor,event,ix,iy,res,result
	integer(4)::iunit, ievent, ikeystate, ixpos, iypos
	!common i4
	INTEGER(2)::status,keystate
	real(8)::scale,t1
	
	
	!open(unit=20,file='user',title='madluo')
	i4=-1
	myscreen.numypixels   = -1
	myscreen.numxpixels   = -1
	myscreen.numtextcols  = -1
	myscreen.numtextrows  = -1
	myscreen.numcolors    = -1
	myscreen.fontsize = -1
	!  myscreen.title = " "C
	open(100,file='user',title='mesh')
	winfo%TYPE = QWIN$MAX
	result = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	!open(100,file='user',title='mesh')
	statusmode = SETWINDOWCONFIG(myscreen)
	iF(.NOT. statusmode) statusmode = SETWINDOWCONFIG(myscreen)
	statusmode = GETWINDOWCONFIG( myscreen )

	result = INSERTMENUQQ (4, 0, $MENUCHECKED, 'Tools'c, nul)
	result = INSERTMENUQQ (4, 1, $MENUENABLED, 'Generate 6-noded element'c, gen_6_noded_triangle)
	result = INSERTMENUQQ (4, 2, $MENUENABLED, 'Generate 15-noded element'c, gen_15_noded_triangle)
	result = INSERTMENUQQ (4, 3, $MENUENABLED, 'Generate 3D model'c, Generate_3D_MODEL)
    !result = INSERTMENUQQ (4, 4, $MENUENABLED, 'Generate gmsh volume model'c, out_volume_to_gmsh)
    result = INSERTMENUQQ (4, 4, $MENUENABLED, 'Generate gmsh volume model'c, out_volume_to_gmsh2)
    result = INSERTMENUQQ (4, 5, $MENUENABLED, 'Output Elevation at nodes'c, out_z_at_node)
    result = INSERTMENUQQ (4, 6, $MENUENABLED, 'Output Elevation at element center'c, out_z_at_center)
    result = INSERTMENUQQ (4, 7, $MENUENABLED, 'Generate 2D DEM particle'c, gen_dem_particle)
    !result = INSERTMENUQQ (4, 5, $MENUENABLED, 'Generate tet element'c, Generate_tet_element)
	result = INSERTMENUQQ (6, 0, $MENUCHECKED, 'Control'c, nul)
	result = INSERTMENUQQ (6, 1, $MENUENABLED, 'display&all'c, displayall)
	result = INSERTMENUQQ (6, 2, $MENUENABLED, 'display&substructure'c, displaysubstructure)
	result = INSERTMENUQQ (6, 3, $MENUENABLED, 'display&mainmesh'c, displaymainmesh)

	result = INSERTMENUQQ (7, 0, $MENUCHECKED, 'Element'c, nul)
	result = INSERTMENUQQ (7, 1, $MENUENABLED, 'TRIELEMENTS'c, tri_element)
	result = INSERTMENUQQ (7, 2, $MENUENABLED, 'DISCONTINUITIES'c, DISCONTINUITIES)
	result = INSERTMENUQQ (7, 3, $MENUENABLED, '防渗墙单元'c, cut_off_wall_element)
	result = INSERTMENUQQ (7, 4, $MENUENABLED, '越流单元'c, nul)
	result = INSERTMENUQQ (7, 5, $MENUENABLED, '接触单元'c, contact_element)
	result = INSERTMENUQQ (7, 6, $MENUENABLED, 'et7'c, sple_element)
	result = INSERTMENUQQ (7, 7, $MENUENABLED, 'penaltyelement'c, penaltyelementshow)
	result = INSERTMENUQQ (7, 8, $MENUENABLED, '所有单元'c, all_element)

	result = INSERTMENUQQ (8, 0, $MENUENABLED, 'Update'c, nul)
	!	result = INSERTMENUQQ (8, 1, $MENUENABLED, 'Boundarycondition'c, Updatebc)
	result = INSERTMENUQQ (8, 2, $MENUENABLED, 'Material'c, UpdateMaterial)
	result = INSERTMENUQQ (9, 0, $MENUENABLED, 'Check'c, nul)
	result = INSERTMENUQQ (9, 1, $MENUENABLED, 'BC_Type'c, check_BcType)
	result = INSERTMENUQQ (9, 2, $MENUENABLED, 'BC_Value'c, check_BcValue)
	result = INSERTMENUQQ (9, 3, $MENUENABLED, 'Node_Number'c, check_number)
	result = INSERTMENUQQ (9, 4, $MENUENABLED, 'Node_Subnum'c, check_subnum)
	result = INSERTMENUQQ (9, 5, $MENUENABLED, 'Element_Bandwidth'c, check_bandwidth)	
	!	result = INSERTMENUQQ (9, 6, $MENUENABLED, 'Element_Material'c, check_material)
	result = INSERTMENUQQ (9, 6, $MENUENABLED, 'Element_KCD'c, check_kcd)
	result = INSERTMENUQQ (9, 7, $MENUENABLED, 'Element_QUALITY'c, EQ)
	result = INSERTMENUQQ (9, 8, $MENUENABLED, 'Element_NUMBER'c, check_ENUM)
	result = INSERTMENUQQ (9, 9, $MENUENABLED, 'Element_ZONE'c, check_EZONE)
	result = INSERTMENUQQ (9, 10, $MENUENABLED, 'Edge_Number'c, check_edge)
	result = INSERTMENUQQ (9, 11, $MENUENABLED, 'Check_Cancel'c, check_null)
	result = INSERTMENUQQ (10, 0, $MENUCHECKED, 'SaveData'c, nul)
	result = INSERTMENUQQ (10, 1, $MENUENABLED, 'To Solver'c, meshoutput)
	result= INSERTMENUQQ (10, 2, $MENUGRAYED, 'Offset'c, offset)
	result = INSERTMENUQQ (10, 3, $MENUGRAYED, 'To Check'c, outplot_checkdata)

	oldcolor=setbkcolorrgb(#00)
	oldcolor=setcolorrgb(#0000ff)
  CALL CLEARSCREEN( $GCLEARSCREEN )


  myscreen.numxpixels=myscreen.numxpixels*0.98
	myscreen.numypixels=myscreen.numypixels*0.88

	CALL SETVIEWPORT( INT2(10), INT2(10), int2(myscreen.numxpixels), int2(myscreen.numypixels) )
  status=RECTANGLE( $GBORDER,INT2(0),INT2(0), int2(myscreen.numxpixels-10), int2(myscreen.numypixels-10))
	status=RECTANGLE( $GBORDER,INT2(1), INT2(1),int2(myscreen.numxpixels-11), int2(myscreen.numypixels-11) )
	
	CALL SETVIEWPORT( INT2(15), INT2(15), int2(myscreen.numxpixels-5), int2(myscreen.numypixels-5) )
	CALL CLEARSCREEN($GVIEWPORT)

	!winxlr=(winylr-winyul)*(myscreen.numxpixels-31)/(myscreen.numypixels-31)+winxul

	if(winylr-winyul>winxlr-winxul) then	
		t1=winxlr
		winylr=winylr+(winylr-winyul)/20
		winyul=winyul-(winylr-winyul)/20
		winxlr=(winylr-winyul)*(myscreen.numxpixels-31)/(myscreen.numypixels-31)+winxul
		winxul=winxul-(winxlr-t1)/2
		winxlr=winxlr-(winxlr-t1)/2
	else
		t1=winylr
		winxlr=winxlr+(winxlr-winxul)/20
		winxul=winxul-(winxlr-winxul)/20		
		winylr=(winxlr-winxul)*(myscreen.numypixels-31)/(myscreen.numxpixels-31)+winyul
		winyul=winyul-(winylr-t1)/2
		winylr=winylr-(winylr-t1)/2
	end if

	
	status=SETWINDOW( .true., winxul, winyul,winxlr, winylr)
	oldcolor=setcolorrgb(#ff00)
	status=RECTANGLE_w( $GBORDER,winxul, winyul,winxlr, winylr)
	

	CALL elementgraph(iunit, ievent, ikeystate, ixpos, iypos)
	

	event = MOUSE$LBUTTONDOWN
	event=  ior(event,MOUSE$LBUTTONUP)
  event = IOR (event, MOUSE$RBUTTONDOWN)
  event = IOR (event, MOUSE$LBUTTONDBLCLK)
  event = IOR (event, MOUSE$RBUTTONDBLCLK)
  event = IOR (event, MOUSE$MOVE)
	i4 = registermouseevent(100, event, elementgraph)
 ! Maximize child window
  result =   SETWSIZEQQ(100, winfo)
	
	
	
  !	if((MOUSE$KS_Lbutton .AND. ikeystate) == MOUSE$KS_Lbutton) then
!	  call getwindowcoord(ixpos, iypos,wxy)
!	  scale=scale/2
 !     CALL CLEARSCREEN($GVIEWPORT)
!	end if

	!event = MOUSE$RBUTTONDOWN .OR. MOUSE$LBUTTONDOWN
	do while(.true.)
    ! res = waitonmouseevent(event, keystate, ix, iy)
	 !i4 = registermouseevent(100, event, elementgraph)
	 res = waitonmouseevent(MOUSE$move, res, ix, iy)
	! if((MOUSE$KS_control .AND. res) == MOUSE$KS_control)   exit
	end do

    END

    
    
    
   SUBROUTINE elementgraph(iunit, ievent, ikeystate, ixpos, iypos)
      use meshDS
      use ds_t
	  USE DFLIB
      use packgen2d
	  implicit none
      integer(4)          oldcolor,iunit, ievent, ikeystate, ixpos, iypos
	  !external,real displayinfo
	  INTEGER(2)          status,result
      INTEGER(4)          res
	  character(48)  msg
	  character(256) nnum
	  !character(1) butter
	  integer::i,n1,nt,j,k
      logical::tof,tof1,tof2
      real(8)::x1,y1,x2,y2,size,scale
	  TYPE (windowconfig):: myscreen
	  common myscreen

	
	  TYPE (wxycoord)       wxy
	  scale=1

		nt=1


	  if(i4>=0) then
        call getwindowcoord(ixpos-15, iypos-15,wxy)
		
		if(ievent==MOUSE$MOVE) then
		   write(msg,'("(x,y)=",f15.7,",",f15.7)') wxy.wx*xyscale+xmin,wxy.wy*xyscale+ymin
           CALL SETMESSAGEQQ (trim(msg), QWIN$MSG_MOUSEINPUTPEND)
		   return
		end if
		
		if((mouse$ks_control.and.ikeystate)==mouse$ks_control) then
		  !call displayinfo(wxy.wx,wxy.wy)
		  !return
		  if(ievent==MOUSE$LBUTTONDOWN) then
		     x1=wxy.wx
			 y1=wxy.wy
		  end if

		  if(ievent==MOUSE$LBUTTONUP) then
             x2=wxy.wx
			 y2=wxy.wy
		     oldcolor=setcolorrgb(#ff00)
		     if(y1<y2) then
		       status=RECTANGLE_w( $GBORDER,x1, y1,x2,y2)
		     else
		       status=RECTANGLE_w( $GBORDER,x2, y2,x1,y1)
		     end if
		!	 res = MESSAGEBOXQQ('Type N for number,H for head,V for velocity,M for modifying node coordinate'C, &
         !    'Show Messages'C,   &
         !    MB$ICONINFORMATION .OR.MB$YESNO.OR.MB$DEFBUTTON1)
	!		 if(res==MB$IDNO) then
	!		   return
	!		 else
		 !  write(msg,*) 'Type N for number,H for head,V for velocity,M for modifying node coordinate.'

			   do while(.true.)
                  !write(msg,*) 'Type N for number,H for head,V for velocity,M for modifying node coordinate.'
				  !CALL SETMESSAGEQQ (msg, QWIN$MSG_MOUSEINPUTPEND)
			      butter = GETCHARQQ()
				  tof2=ichar(butter)==ichar('n')
                  tof2=tof2.or.ichar(butter)==ichar('s')
				  !tof2=tof2.or.ichar(butter)==ichar('v')
				  !tof2=tof2.or.ichar(butter)==ichar('m')
				  if(tof2) exit
               end do
	!		 end if

			 oldcolor=setcolorrgb(#0000ff)
			 result = INITIALIZEFONTS()
             result = SETFONT('t''Arial''h18w10p')
		     do i=1,nnode
		       tof=node(i).x<min(x1,x2)
			   if(tof) cycle
               tof=node(i).x>max(x1,x2)
               if(tof) cycle
			   tof=node(i).y<min(y1,y2)
               if(tof) cycle
               tof=node(i).y>max(y1,y2)
               if(tof) cycle
			   size=node(i).s/10
			   oldcolor=setcolorrgb(#0000ff)
			   status=ellipse_W($GBORDER,node(i).x-size,node(i).y-size,node(i).x+size,node(i).y+size)
			   select case(ichar(butter))
			      case(ichar('n'))
			         write(nnum,'(i6)') node(i).number
				  case(ichar('s'))

				  case(ichar('v'))

				  case(ichar('m'))
                     !call mc(node(i).x,node(i).y)

			   end select
			
			   call moveto_w(node(i).x-size,node(i).y-size,wxy)
			   oldcolor=setcolorrgb(#FFFFFF)
			   call outgtext(trim(adjustL(nnum)))		

		     end do

		   end if
           return
		end if
		
		if((MOUSE$KS_lbutton .AND. ikeystate) == MOUSE$KS_lbutton) scale=scale/2
	          	
		if((MOUSE$KS_Rbutton .AND. ikeystate) == MOUSE$KS_Rbutton) scale=scale*2
		
	
		
	      winxul=wxy.wx-(wxy.wx-winxul)*scale
		  winyul=wxy.wy-(wxy.wy-winyul)*scale
		  winxlr=wxy.wx+(winxlr-wxy.wx)*scale
		  winylr=wxy.wy+(winylr-wxy.wy)*scale
          CALL SETVIEWPORT( INT2(15), INT2(15), int2(myscreen.numxpixels-5), int2(myscreen.numypixels-5) )
		  CALL CLEARSCREEN($GVIEWPORT)
		  oldcolor=setcolorrgb(#0000ff)
		  ! winxlr=(winylr-winyul)*770/570+winxul
	      status=SETWINDOW( .true., winxul, winyul,winxlr , winylr)
	      oldcolor=setcolorrgb(#ff00)
	      status=RECTANGLE_w( $GBORDER,winxul, winyul,winxlr , winylr)
		
      end if

      call  particle2D.plot()
      
      
 	 oldcolor=setcolorrgb(#0000ff)
	 result = INITIALIZEFONTS()
     result = SETFONT('t''Arial Narrow''h18w10p')
	  do i=1,1	
		 select case(controldisplay)
		    case(1)
           ept=1
			case(2)
			   Ept=1
		    case(0)
	           Ept=1		   
		
		 end select

		
100 do while(ept<=nelt)
      if(elt(ept).isdel)then
     		ept=ept+1
     		cycle     
      end if
			select case(element_type)
               case(-1)
                 if(elt(ept).et/=-1) then
                    ept=ept+1
                    cycle
                 end if
               case(0)
                 if(elt(ept).et/=0) then
                    ept=ept+1
                    cycle
                 end if
               case(1)
                 if(elt(ept).et/=1) then
                    ept=ept+1
                    cycle
                 end if
               case(2)
                 if(elt(ept).et/=2) then
                    ept=ept+1
                    cycle
                 end if
               case(3)
                 if(elt(ept).et/=3) then
                    ept=ept+1
                    cycle
                 end if
               case(4)
                 if(elt(ept).et/=4) then
                    ept=ept+1
                    cycle
                 end if
               case(5)
                 if(elt(ept).et/=5) then
                    ept=ept+1
                    cycle
                 end if
               case(7)
                 if(elt(ept).et/=7) then
                    ept=ept+1
                    cycle
                 end if
               case(8)
                if(elt(ept).et/=8) then
                    ept=ept+1
                    cycle
                 end if
				case default					
					
            end select

	
		 tof1=max(node(elt(ept).node(1)).x,node(elt(ept).node(2)).x,node(elt(ept).node(3)).x)<min(winxlr,winxul)
		 if(tof1) then
		   ept=ept+1
		   cycle
		 end if
		 tof1=min(node(elt(ept).node(1)).x,node(elt(ept).node(2)).x,node(elt(ept).node(3)).x)>max(winxlr,winxul)
         if(tof1) then
		   ept=ept+1
		   cycle
		 end if
		 tof1=max(node(elt(ept).node(1)).y,node(elt(ept).node(2)).y,node(elt(ept).node(3)).y)<min(winylr,winyul)
		 if(tof1) then
		   ept=ept+1
		   cycle
		 end if
		 tof1=min(node(elt(ept).node(1)).y,node(elt(ept).node(2)).y,node(elt(ept).node(3)).y)>max(winylr,winyul)
         if(tof1) then
		   ept=ept+1
		   cycle
		 end if

		 select case(mod(elt(ept).zn,15))
		   case(1)
              oldcolor=setcolorrgb($HIWHITE)
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
              oldcolor=setcolorrgb($HIRED)
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
		
  		 call moveto_w(node(elt(ept).node(1)).x,node(elt(ept).node(1)).y,wxy)
		 status=lineto_w(node(elt(ept).node(2)).x,node(elt(ept).node(2)).y)
		 call moveto_w(node(elt(ept).node(2)).x,node(elt(ept).node(2)).y,wxy)
         status=lineto_w(node(elt(ept).node(3)).x,node(elt(ept).node(3)).y)
         call moveto_w(node(elt(ept).node(3)).x,node(elt(ept).node(3)).y,wxy)
		 if(elt(ept).et==-1) then
			status=lineto_w(node(elt(ept).node(4)).x,node(elt(ept).node(4)).y)
			call moveto_w(node(elt(ept).node(4)).x,node(elt(ept).node(4)).y,wxy)
		 end if
         status=lineto_w(node(elt(ept).node(1)).x,node(elt(ept).node(1)).y)

		 oldcolor=setcolorrgb(#FFFF00)
		  result = INITIALIZEFONTS()
      result = SETFONT('t''Arial''h15w8p')

		 select case(showvalue)
			
			case(1)
			   
			   do k=1,3
				   write(nnum,'(i6)') node(elt(ept).node(k)).number
				   call moveto_w(node(elt(ept).node(k)).x,node(elt(ept).node(k)).y,wxy)			
				   call outgtext(trim(adjustL(nnum)))
			   end do
			
				if(elt(ept).et==6) then
					do k=4,6
					   write(nnum,'(i6)') node(elt(ept).node(k)).number
					   call moveto_w(node(elt(ept).node(k)).x,node(elt(ept).node(k)).y,wxy)
					   call outgtext(trim(adjustL(nnum)))
				   end do
				end if

				if(elt(ept).et==15) then
					do k=4,15
					   write(nnum,'(i6)') node(elt(ept).node(k)).number
					   call moveto_w(node(elt(ept).node(k)).x,node(elt(ept).node(k)).y,wxy)
					   call outgtext(trim(adjustL(nnum)))
				   end do
				end if

			case(2)

			case(3)
			   
			case(4)
			

			case(5)

			case(6)
			   !oldcolor=setcolorrgb(#FFFFFF)
			   write(nnum,'(i6)') node(elt(ept).node(1)).bw
			   call moveto_w(node(elt(ept).node(1)).x,node(elt(ept).node(1)).y,wxy)			
			   call outgtext(trim(adjustL(nnum)))
			
			   write(nnum,'(i6)') node(elt(ept).node(2)).bw
			   call moveto_w(node(elt(ept).node(2)).x,node(elt(ept).node(2)).y,wxy)
			   call outgtext(trim(adjustL(nnum)))

			   write(nnum,'(i6)') node(elt(ept).node(3)).bw
			   call moveto_w(node(elt(ept).node(3)).x,node(elt(ept).node(3)).y,wxy)
			   call outgtext(trim(adjustL(nnum)))
			
			 case(7)
			   msg=''
			   write(nnum,'(i6)') elt(ept).kcd
			   msg=trim(msg)//'\'//trim(adjustL(nnum))
			   write(nnum,'(i6)') elt(ept).maxedge
			   msg=trim(msg)//'\'//trim(adjustL(nnum))
			
			   call moveto_w((node(elt(ept).node(1)).x+node(elt(ept).node(2)).x+node(elt(ept).node(3)).x)/3,(node(elt(ept).node(1)).y+node(elt(ept).node(2)).y+node(elt(ept).node(3)).y)/3,wxy)			
			   call outgtext(trim(msg))
			
			 case(8)
			   msg=''
			   write(nnum,'(i6)') ept
			   msg=trim(adjustL(nnum))
               write(nnum,'(i6)') elt(ept).NUMBER 
               msg=trim(adjustl(msg))//'N'//trim(adjustL(nnum))
			   call moveto_w(sum(node(elt(ept).node(1:elt(ept).nnum)).x)/elt(ept).nnum,sum(node(elt(ept).node(1:elt(ept).nnum)).y)/elt(ept).nnum,wxy)			
			   call outgtext(trim(msg))	
			   			 					
			 case(9)
			   msg=''
			   write(nnum,'(i6)') elt(ept).zn
			   msg=trim(adjustL(nnum))			
			   call moveto_w((node(elt(ept).node(1)).x+node(elt(ept).node(2)).x+node(elt(ept).node(3)).x)/3,(node(elt(ept).node(1)).y+node(elt(ept).node(2)).y+node(elt(ept).node(3)).y)/3,wxy)			
			   call outgtext(trim(msg))
			 case(10)
			   do k=1,elt(ept).nnum
                   !IF(EDGE(elt(ept).edge(k)).ISZONEBC/=-1.AND.EDGE(elt(ept).edge(k)).ISCEDGE==0) THEN
				       write(nnum,'(i6)') elt(ept).edge(k)
                       msg=trim(adjustL(nnum))
                       write(nnum,'(i6)') edge(elt(ept).edge(k)).num
                       msg=trim(adjustl(msg))//'N'//trim(adjustL(nnum))
				       x1=(node(edge(elt(ept).edge(k)).v(1)).x+ &
				   	    node(edge(elt(ept).edge(k)).v(2)).x)/2
				       y1=(node(edge(elt(ept).edge(k)).v(1)).y+ &
				   	    node(edge(elt(ept).edge(k)).v(2)).y)/2				   	
				       call moveto_w(x1,y1,wxy)			
				       call outgtext(trim(adjustL(msg)))
                   !ENDIF
			   end do
		 end select		


		 ept=ept+1
	
		 end do


      end do

      

    end	subroutine

    subroutine gen_dem_particle()
        use packgen2d
        implicit none       
        

        call  particle2D.gen_particle()
        
        

    endsubroutine
    
    
    subroutine displaysubstructure()
	   use ds_t
	   implicit none	
	   controldisplay=1
	end subroutine

    subroutine displaymainmesh()
	   use ds_t
	   implicit none	
	   controldisplay=2
	end subroutine

	subroutine displayall()
	   use ds_t
	   implicit none	
	   controldisplay=0
	end subroutine

    subroutine tri_element()
       use ds_t
       implicit none
       element_type=0
    end subroutine

     subroutine DISCONTINUITIES()
       use ds_t
       implicit none
       element_type=-1
    end subroutine
    subroutine contact_element()
       use ds_t
       implicit none
       element_type=5
    end subroutine
    subroutine cut_off_wall_element()
       use ds_t
       implicit none
       element_type=-1
    end subroutine
    subroutine sple_element()
       use ds_t
       implicit none
       element_type=7
    end subroutine
    
    subroutine penaltyelementshow()
       use ds_t
       implicit none
       element_type=8
    end subroutine



    subroutine all_element()
       use ds_t
       implicit none
       element_type=-99
    end subroutine


   subroutine UpdateMaterial()
	  use meshds
	  use dflib
	  implicit none
	  integer::i,msg
	  character(128)::term
!	  call element_material()
	  term="Update material completely"
	  term=trim(term)
	  msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$OK.OR.MB$DEFBUTTON1)
	  if(msg==MB$IDOK) return
   end subroutine

   subroutine check_null()
	  use ds_T
	  implicit none
	  showvalue=-1
   end subroutine

   subroutine check_bctype()
	  use ds_T
	  implicit none
	  showvalue=3
   end subroutine

   subroutine check_bcvalue()
	  use ds_T
	  implicit none
	  showvalue=4
   end subroutine

   subroutine check_number()
	  use ds_T
	  implicit none
	  showvalue=1
   end subroutine

   subroutine check_subnum()
	  use ds_T
	  implicit none
	  showvalue=2
   end subroutine

   subroutine check_material()
	  use ds_T
	  implicit none
	  showvalue=5
   end subroutine

   subroutine check_bandwidth()
	  use ds_T
	  implicit none
	  showvalue=6
   end subroutine

   subroutine check_kcd()
	  use ds_T
	  implicit none
	  showvalue=7
   end subroutine
   
   subroutine check_edge()
   	use Ds_T
   	implicit none
   	showvalue=10
   end subroutine

   subroutine check_ENUM()
	  use ds_T
	  implicit none
	  showvalue=8
   end subroutine
   subroutine check_EZONE()
	  use ds_T
	  implicit none
	  showvalue=9
   end subroutine

   subroutine EQ()
	  use meshds
	  implicit none
	  real(8)::coe
	  integer::n1,ar(0:9)=0,sum

	  ar=0
	  do ept=1,nelt
	  	if(elt(ept).isdel.or.elt(ept).et/=0) cycle
		 call Rr(ept,coe)
		 n1=int(coe/0.1)
		 ar(n1)=ar(n1)+1
	  end do
	
	  sum=0
	  do n1=0,9
		 sum=sum+ar(n1)
	  end do

	  print *, 'Element Quality statistics:'
	  print *, '    区间    ','       单元数     ','     %     '
	  write(*,'(a12,i14,f9.2)') '[0,0.1)',ar(0),ar(0)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.1,0.2)',ar(1),ar(1)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.2,0.3)',ar(2),ar(2)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.3,0.4)',ar(3),ar(3)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.4,0.5)',ar(4),ar(4)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.5,0.6)',ar(5),ar(5)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.6,0.7)',ar(6),ar(6)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.7,0.8)',ar(7),ar(7)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.8,0.9)',ar(8),ar(8)/real(sum)*100
	  write(*,'(a12,i14,f9.2)') '[0.9,1)',ar(9),ar(9)/real(sum)*100
	  		
   end subroutine


	!给出单元EL的形状系数coe：2Ri/Rc,Rc,Ri分别三角形单元的外接圆和内接圆半径。
	subroutine Rr(el,coe)
	   use meshds
	   implicit none
	   real(8)::Ri,Rc,L1,L2,L3
	   integer,intent(in)::EL
	   real(8),intent(out)::coe
	
	   L1=((node(elt(el).node(2)).x-node(elt(el).node(3)).x)**2+(node(elt(el).node(2)).y-node(elt(el).node(3)).y)**2)**0.5
	   L2=((node(elt(el).node(1)).x-node(elt(el).node(3)).x)**2+(node(elt(el).node(1)).y-node(elt(el).node(3)).y)**2)**0.5
	   L3=((node(elt(el).node(2)).x-node(elt(el).node(1)).x)**2+(node(elt(el).node(2)).y-node(elt(el).node(1)).y)**2)**0.5

	   Rc=(L1*L2*L3)/((L1+L2+L3)*(L1+L2-L3)*(L1-L2+L3)*(L2+L3-L1))**0.5
	   Ri=2*elt(el).property(1)/(L1+L2+L3)
	   coe=2*Ri/Rc
	end subroutine

subroutine offset()
	use meshds
	implicit none
	real(8),allocatable::ec(:,:) 
	integer::i
	real(8)::minx,miny,xc,yc
	
	allocate(ec(enumber,2)) 
	minx=minval(node(1:nnode).x)
	miny=minval(node(1:nnode).y)

	ec=0.0D0
	do ept=1,nelt
		if(elt(ept).isdel) cycle	
		if(elt(ept).et==0) then
			!element centroids
			xc=(node(elt(ept).node(1)).x+node(elt(ept).node(2)).x+node(elt(ept).node(3)).x)/3
			yc=(node(elt(ept).node(1)).y+node(elt(ept).node(2)).y+node(elt(ept).node(3)).y)/3
			!let new element centroid scale by 1.1 times
			ec(elt(ept).number,1)=(xc-minx)*1.15+minx
			ec(elt(ept).number,2)=(yc-miny)*1.15+miny
			!element centroids offset distance
			ec(elt(ept).number,1)=ec(elt(ept).number,1)-xc !dx,dy
			ec(elt(ept).number,2)=ec(elt(ept).number,2)-yc
		end if
	end do

	do ept=1,nelt
		if(elt(ept).isdel) cycle	
		if(elt(ept).et==0) then
			node(elt(ept).node(1)).x=node(elt(ept).node(1)).x+ec(elt(ept).number,1)
			node(elt(ept).node(1)).y=node(elt(ept).node(1)).y+ec(elt(ept).number,2)
			node(elt(ept).node(2)).x=node(elt(ept).node(2)).x+ec(elt(ept).number,1)
			node(elt(ept).node(2)).y=node(elt(ept).node(2)).y+ec(elt(ept).number,2)
			node(elt(ept).node(3)).x=node(elt(ept).node(3)).x+ec(elt(ept).number,1)
			node(elt(ept).node(3)).y=node(elt(ept).node(3)).y+ec(elt(ept).number,2)
		end if			
	end do

	end subroutine

