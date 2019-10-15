  
  program  mainmesh
		use meshds
		USE DFPORT
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
			call Removesupertri()			
			call RemoveT()
			call TIME(char_time)
			if(isnorefined==0) then
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
		
		call seg_initialize()
		
		if(nsm>0) call structuremesh()
		call group()
        
		call elementgroup()
		call limitanalysisgrid()
		call time(char_time)	 
		print *, '主域网格处理完成，请检查网格。如果要保存数据,请按SaveData…: ', char_time
		CALL CPU_TIME ( time_end )
		print *, 'nnode=',nnode,'ENUMBER=',ENUMBER
		PRINT *, 'Time of operation was ', time_end - time_begin, ' seconds'
		!call meshoutput()
		!call debug()
		CALL CLEAR_EDGE_ELEMENT()
		call graph

  end program

  
  SUBROUTINE CLEAR_EDGE_ELEMENT()
	USE MESHDS
	IMPLICIT NONE
	INTEGER::I,N1=0
	N1=0
	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		N1=N1+1
		EDGE(I).NUM=N1
	end do
	N1=0
	do i=1,NELT
		if(ELT(i).ISDEL) cycle
		N1=N1+1
		ELT(I).NUMBER=N1
	end do
	ENUMBER=N1
	
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

subroutine debug()
	use meshds
	implicit none
	integer::i

	open(3,file=trim(path_name)//'_bug.dat',status='replace')
	write(3,'(a5)') 'NODES'
	do i=1,nnode
		write(3,'(i6,2f15.6)') i,node(i).x,node(i).y
	end do
	write(3,'(a5)') 'EDGES'
	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		write(3,'(3I6)') i,edge(i).v 
	end do
	write(3,'(a8)') 'ELEMENTS'
	do ept=1,nelt
		if(elt(ept).isdel) cycle
		write(3,'(7I6)') elt(ept).number,elt(ept).node(1:3),elt(ept).edge(1:3)
	end do
	close(3)

end subroutine

  



