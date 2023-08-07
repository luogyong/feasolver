  
  program  mainmesh
		use meshds
		USE DFPORT
        use CutoffWall
        use triangle_io
		character*8 char_time
		integer::i,n1,nc,iflag
		REAL time_begin, time_end

		CALL CPU_TIME ( time_begin )

		Print *, 'MESH.LGY WORKS.'
		call TIME(char_time )
		print *, 'Begin to read in data...', char_time
		allocate(node(-3:Maxnnode),adjlist(Maxnnode),elt(maxnelement))
		call meshinp
		call checkdata() 
		call TIME(char_time)  	 
		print *, 'Read in data completed.Begin to allocate space...', char_time

        if(libtriangle.method>=0) then
            call libtriangle.exe(path_name)
            call libtriangle.getdata(element=elt,node=node,edge=edge,adjlist=adjlist,cedge=cedge)
            nelt=size(elt);nnode=size(node);nedge=size(edge);ncedge=size(cedge);            
        else
            
		    !对全区变量用区域变量的相应值进行初始化
		    !call initialize_m(i)	  
		    call TIME(char_time)
		    print *, 'Begin to insert nodes in control lines... ', char_time	
		    !call BPInsert(i)
		    call CLinsert
		    call TIME(char_time)
		    print *, 'Control lines have been divided. Begin to generate the initial mesh, please wait…: ', char_time	

		    if(keypn>0) then
    !			call BPmesh
			    call kpmesh
			    call TIME(char_time)
			    print *, ''
			    write(*,*), 'Initial mesh has been generated. Begin to repair control line, please wait…: ', char_time		 

			    iflag=0
			    call RCL_4th(iflag)
			    !call Removesupertri()			
			    call RemoveT()
			    call TIME(char_time)
			    if(isnorefined<1) then
				    print *, 'Control lines have been repaired in the initial mesh. Begin to refine the initial mesh, please warit…: ', char_time
				    call InsertPoint()
				    call TIME(char_time)
				    print *, 'Refine mesh completed. Begin to repair control lines again, please wait…: ', char_time
				    iflag=1
				    call RCL_4th(iflag)
				    call time(char_time)
			    endif
			    print *, 'Delaunay triangular Mesh completed. Begin to  handle other tasks, please wait…: ', char_time	 	 	 	 
		    end if
		endif
		!call seg_initialize()
		
		if(nsm>0) call structuremesh()
		!call group()
        CALL tri_elt_group2d()
        !do i=1,ncow
        !    call cowall(i).gen_element()
        !enddo        
		call elementgroup()
		call limitanalysisgrid()
		call time(char_time)	 
		print *, '主域网格处理完成，请检查网格。如果要保存数据,请按SaveData…: ', char_time
		CALL CPU_TIME ( time_end )
		print *, 'nnode=',nnode,'ENUMBER=',ENUMBER
		PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds'
		!call meshoutput()
		
        
		CALL CLEAR_EDGE_ELEMENT()
        !call debug()
		call graph

  end program

  
  SUBROUTINE CLEAR_EDGE_ELEMENT()
	USE MESHDS
    USE geomodel,ONLY:CHECK_ORIENT
	IMPLICIT NONE
	INTEGER::I,N1=0
    REAL(8)::VEC1(2,2)
	N1=0
	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		N1=N1+1
		EDGE(I).NUM=N1
	end do
	N1=0
    mesharea=0.d0
	do i=1,NELT
		if(ELT(i).ISDEL) cycle
		N1=N1+1
		ELT(I).NUMBER=N1
        VEC1(1,:)=NODE(ELT(I).NODE([2,3])).X-NODE(ELT(I).NODE(1)).X
        VEC1(2,:)=NODE(ELT(I).NODE([2,3])).Y-NODE(ELT(I).NODE(1)).Y
        ELT(I).PROPERTY(1)=0.5*ABS(VEC1(1,1)*VEC1(2,2)-VEC1(1,2)*VEC1(2,1)) !AREA
        mesharea=mesharea+ELT(I).PROPERTY(1)
	end do
	ENUMBER=N1
    
    CALL CHECK_ORIENT()
	
ENDSUBROUTINE
    
  subroutine group()
		use meshds
		implicit none
		integer::i
		real(8)::xc,yc
		logical::tof

		do ept=1,nelt
			if(elt(ept).isdel) cycle
			if(znum>1) then
				xc=(node(elt(ept).node(1)).x+node(elt(ept).node(2)).x+node(elt(ept).node(3)).x)/3
				yc=(node(elt(ept).node(1)).y+node(elt(ept).node(2)).y+node(elt(ept).node(3)).y)/3
				do i=1,znum
                    if(zone(i).outgmshtype==2) cycle
					call judgezone(xc,yc,zone(i),tof)                   
					if(tof) then
                        elt(ept).zn=i 
                        exit
					 end if
				end do
				do i=1,znum
                    if(zone(i).outgmshtype==2) then
					    call judgezone(xc,yc,zone(i),tof)                   
					    if(tof) then
                            elt(ept).zn2=i 
                            exit
					     end if
                     endif
				end do                
			end if

		end do

    end subroutine
    
    subroutine tri_elt_group2d()
		use meshds
		implicit none
		integer::i,j,K,IELT1,N1,N2,iseg1,ielt2
		real(8)::xc,yc,PT(2)
		logical::tof
        integer::stack1(10000),EDGE1(NEDGE),EDGE2(NEDGE),ISCHK1(NELT),IA1(NNODE)
        
        ELT(1:NELT).ZN=0
        
        call SETUP_SEARCH_ZONE_2D(SEARCHZONE,NSZONE,(Xmax-Xmin)/xyscale,0.d0,(Ymax-Ymin)/xyscale,0.d0,NODE(1:NNODE),elt(1:nelt))
        EDGE1=0;
        DO I=1,ZNUM
            ISCHK1=0
            IF(ZONE(I).FORMAT==1) THEN
                EDGE2=0
                DO K=1,ZONE(I).NUM
                    IELT1=POINTlOC_2D(ZONE(I).POINT(:,K),-1,SEARCHZONE,NODE(1:NNODE),ELT)
                    IF(IELT1<1) THEN
                        PRINT *, 'NO TRIANGLE ELEMENT IN ZONE(I).I=',I
                        CYCLE
                    ENDIF
                    STACK1(1)=IELT1
                    ELT(IELT1).ZN=I
                    if(zone(i).outgmshtype>=2) ELT(IELT1).ZN2=I
                    N2=1
                    DO WHILE(N2>0)
                        IELT1=STACK1(N2)                
                        N2=N2-1
                        IF(ISCHK1(IELT1)==1) CYCLE
                        ISCHK1(IELT1)=1
                        DO J=1,ELT(IELT1).NNUM
                            N1=ELT(IELT1).EDGE(J)
                            EDGE2(N1)=EDGE2(N1)+1
                            IF(EDGE(N1).ISCEDGE==0.AND.EDGE1(N1)/=I) THEN
                                EDGE1(N1)=I                            
                                IELT2=ELT(IELT1).ADJ(J)
                                IF(IELT2>0) THEN
                                    IF(ISCHK1(IELT2)==0) THEN
                                        N2=N2+1
                                        STACK1(N2)=IELT2
                                        ELT(IELT2).ZN=I
                                        if(zone(i).outgmshtype>=2) ELT(IELT2).ZN2=I
                                    ENDIF
                                ENDIF
                            ENDIF
                        
                        ENDDO
            
                    ENDDO
                ENDDO
                ZONE(I).NBE=COUNT(EDGE2==1)
                ZONE(I).BEDGE=PACK([1:NEDGE],EDGE2==1)
                DO J=1,ZONE(I).NBE
                    IA1(EDGE(ZONE(I).BEDGE(J)).V)=1
                ENDDO
                ZONE(I).BNODE=PACK([1:NNODE],IA1(1:NNODE)==1) !NOT SORTED
                zone(i).nbn=size(zone(i).bnode)
            ELSE
                if(zone(i).nbe==0) then
                    do j=1,zone(i).num
                        iseg1=segindex(zone(i).cp(j),zone(i).cp(mod(j,zone(i).num)+1))
                        if(iseg1<1) cycle
                        zone(i).bedge=[zone(i).bedge,seg(iseg1).get_edge(zone(i).point(:,j))]
                        if(.not.allocated(zone(i).bnode)) then
                            zone(i).bnode=[zone(i).bnode,seg(iseg1).get_node(zone(i).point(:,j))]
                        else                            
                            zone(i).bnode=[zone(i).bnode(:size(zone(i).bnode)-1),seg(iseg1).get_node(zone(i).point(:,j))]
                        endif
                    enddo
                    zone(i).bedge=abs(zone(i).bedge)
                    zone(i).nbe=size(zone(i).bedge)
                    zone(i).nbn=size(zone(i).bnode)-1 !ZONE BC IS CLOSE. EXCLUDE THE LAST NODE.
                endif                
                
                PT=Find_Point_Inside_Polygon_2D(zone(i).point)
                IELT1=POINTlOC_2D(PT,-1,SEARCHZONE,NODE(1:NNODE),ELT)
                IF(IELT1<1) THEN
                    PRINT *, 'NO TRIANGLE ELEMENT IN ZONE(I).I=',I
                    CYCLE
                ENDIF
                STACK1(1)=IELT1
                IF(ELT(IELT1).ZN==0) ELT(IELT1).ZN=I !这种格式要区分大区套小区的情况
                if(zone(i).outgmshtype>=2.and.ELT(IELT1).ZN2==-1) ELT(IELT1).ZN2=I
                N2=1
                
                edge2=0
                edge2(zone(i).bedge)=1
                
                DO WHILE(N2>0)
                    IELT1=STACK1(N2)                     
                    N2=N2-1
                    IF(ISCHK1(IELT1)==1) CYCLE
                    ISCHK1(IELT1)=1
                    DO J=1,ELT(IELT1).NNUM
                        N1=ELT(IELT1).EDGE(J)
                        IF(EDGE2(N1)==0.AND.EDGE1(N1)/=I) THEN
                            EDGE1(N1)=I                            
                            IELT2=ELT(IELT1).ADJ(J)
                            IF(IELT2>0) THEN
                                IF(ISCHK1(IELT2)==0) THEN
                                    N2=N2+1
                                    STACK1(N2)=IELT2
                                    IF(ELT(IELT2).ZN==0) ELT(IELT2).ZN=I !这种格式要区分大区套小区的情况
                                    if(zone(i).outgmshtype>=2.and.ELT(IELT2).ZN2==-1) ELT(IELT2).ZN2=I
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        

    end subroutine  
    
    


subroutine debug()
	use meshds
	implicit none
	integer::i,n1

	open(3,file=trim(path_name)//'_bug.dat',status='replace')
	write(3,'(a5)') 'NODES'
    WRITE(3,30) 
	do i=1,nnode
        N1=ADJLIST(I).COUNT
		write(3,31) i,node(i).x,node(i).y,ADJLIST(I).NODE(1:N1),ADJLIST(I).EDGE(1:N1)
	end do
	write(3,'(a5)') 'EDGES'
    WRITE(3,20) 
	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		write(3,'(<3>(I6,X))') i,edge(i).v 
	end do
	write(3,'(a8)') 'ELEMENTS'
    WRITE(3,10) 
	do ept=1,nelt
		if(elt(ept).isdel) cycle
        n1=ELT(EPT).NNUM
		write(3,11) elt(ept).number,elt(ept).node(1:N1),elt(ept).edge(1:N1),elt(ept).adj(1:N1)
	end do
	close(3)

10 format(5X,'NO',3('NODEI',X),3(3X,'EDGEI'),3(3X,'ADJEI'))
11 FORMAT(<1+3*N1>(I6,X))
20 FORMAT(5X,'NO',5X,'V1',5X,'V2')
30 FORMAT(5X,'NO',15X,'X',15X,'Y',(2X,'ADJNI'),(2X,'ADJNE')) 
31 FORMAT(I6,X,2(F15.6,X),<N1>(I6,X),<N1>(I6,X))
end subroutine

  



