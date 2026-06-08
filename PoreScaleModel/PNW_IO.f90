
program PNWmain
	use meshDS
	use PoreScaleModel
	use dflib
    use ifport
    use omp_lib
    
	implicit none


	integer::i,j,k,k1,n1,n2
	!integer::a1
	integer::ef,unit,itype
	character(256) term,keyword
	character(1)::ch
	real(8)::r_t !  外包三角形的内接圆半径，r_t=5*（包含边界关键点的最小圆半径）
	real(8)::x_t,y_t
	real(8)::t1,t2,t3
	integer(4)::length,result,msg,status,cnt
	character(256)::nme,filename
	CHARACTER(3)        drive
	CHARACTER(256)      dir
	CHARACTER(256)      name
	CHARACTER(256)      ext
	!type(constrainline_tydef) ,allocatable::csl_t(:)
	!type(qwinfo) winfo
 !
	!term="Input Files for Mesh(*.mes),*.mes;Data Files(*.dat),*.dat;All Files(*.*),*.*;"
	!term=trim(term)
	!call setmessageqq(term,QWIN$MSG_FILEOPENDLG)
	!winfo%TYPE = QWIN$MAX
	!result = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	!result=SETWSIZEQQ(0, winfo) 
	!term=''
    cnt = command_argument_count ()
    if(cnt==0) then
        call psmodel.help(psmodel.helpstring)
        stop
    endif
    call get_command_argument (1, filename, length, status)
    if(trim(adjustl(filename))=='-h'.or.trim(adjustl(filename))=='-help') then
        call psmodel.help(psmodel.helpstring)
        stop
    endif
	open(1,file=filename(1:length),status='old' )
	inquire(1,name=nme)
	length = SPLITPATHQQ(nme, drive, dir, name, ext)
	msg = CHDIR(trim(drive)//trim(dir))
	path_name=trim(drive)//trim(dir)//trim(name)
	title=trim(name)
	unit=1
	itype=0
    
    call read_execute(unit,itype,keyword)
	
	
	resultfile=trim(path_name)//'.SINP'		
 !
	!checkfile=trim(path_name)//'_check.plot'
	

	close(1)


print *,'Reading data COMPLETED.'
!if(znum==0) znum=1




end 



subroutine command(term,unit)
 use dflib
 use meshds
 use PoreScaleModel
 implicit none
 integer::i,j,k,j1

 integer::nt1,nt2,n1,n2,ef,unit,ni,n3
 integer(4)::msg,oldcolor
 character(len=*) term
 character(len=1024)::string,str1
 character(1)::ch
 CHARACTER(3)::SUPNO
 character(16)::legalC
 character(16)::noun
 integer::dnmax,dn     
 real(8)::t1,t2
 real(8),allocatable::ar(:),ar1(:)
 !type(arr_tydef)::ts1
 logical::isall1,iserror1
 
!INTERFACE
!	SUBROUTINE strtoint(unit,ar,nmax,n1,num_read,isall)
!		INTEGER:: unit, nmax, n1,num_read
!		real(8)::ar(nmax)
!		logical::isall
!		OPTIONAL::isall
!	END SUBROUTINE
!END INTERFACE
 
 
 dnmax=300
 dn=0
 if(.not.allocated(ar)) allocate(ar(dnmax))
 ar=0.d0
 !unit=1

 legalc='.0123456789+-eE*'
	 
 term=trim(term)

 select case(term)
 case('psm')
	call psmodel.handle(unit)
 
	 case default
	 
		term="No such key word!  "//trim(term)
		msg = MESSAGEBOXQQ(term,'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
		if(msg==MB$IDOK) then
		  stop
		end if

  end select


end subroutine

!
!
!subroutine indataerror(n1,noun,unit)
!   use dflib
!   implicit none
!   integer::n1,ef,unit
!   integer(4)::msg
!   character(256)::str1
!   character(16)::noun
!   character(64)::str2
!   character(16)::legalC
!   character(5)::ch
! 
!   legalc='.0123456789+-eE*'
!   write(ch,'(i5)') n1-1
!   write(str2,'(a64)') trim(noun)//' number is wrong.have read in numbers is'//trim(ch)//'.Please Check.'
!   do while(.true.)				
!	  read(unit,'(a256)',iostat=ef) str1					
!	  if(ef<0) then
!		 msg = MESSAGEBOXQQ(trim(str2),'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
!		 if(msg==MB$IDOK) then
!		   stop
!		 end if
!	  end if
!	  str1=adjustL(str1)
!	  if(len_trim(str1)==0) cycle
!	  if(str1(1:1)/='/') then
!		 if(index(legalC,str1(1:1))==0) then
!			msg = MESSAGEBOXQQ(trim(str2),'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
!			if(msg==MB$IDOK) then
!			   stop
!			end if
!		 end if
!		 backspace(unit)
!		 exit
!	 end if
!   end do
!
!end subroutine
!




!elarge Array Node size by increment of 10000
!Set Array Adjlist size equel to the size of Node
!Set Array Edge size equel to 3*Nnode+1+cln-2,according to
!Euler equation node+element-edge=2-boundary(including innner boundary)
!Subroutine EnlargeNodeRelative()
! use meshds
! implicit none
! integer::i,nc1
! type(point_tydef),allocatable::node1(:)
! type(adjlist_tydef),allocatable::adjlist1(:)
! type(edge_tydef),allocatable::edge1(:)
! real(8),allocatable::elevation1(:,:),nlayer1(:,:)
! 
! allocate(node1(-3:nnode))
! node1=node
! deallocate(node)	
! maxnnode=maxnnode+10000
! allocate(node(-3:maxnnode))
! node(-3:nnode)=node1
! deallocate(node1)
! 
! !enlarge adjlist as it must keep in the same size with node
! allocate(adjlist1(nnode))
! do i=1,nnode
!	 nc1=size(adjlist(i).node)
!	 allocate(adjlist1(i).node(nc1), &
!		 adjlist1(i).edge(nc1))
!	 adjlist1(i).count=adjlist(i).count
!	 adjlist1(i).node=adjlist(i).node
!	 adjlist1(i).edge=adjlist(i).edge
! end do
! 
! deallocate(adjlist)
! allocate(adjlist(maxnnode))
! do i=1,nnode
!	 nc1=size(adjlist1(i).node)
!	 allocate(adjlist(i).node(nc1), &
!		 adjlist(i).edge(nc1))
!	 adjlist(i).count=adjlist1(i).count
!	 adjlist(i).node=adjlist1(i).node
!	 adjlist(i).edge=adjlist1(i).edge
! end do	
! deallocate(adjlist1)
! 
! !if(allocated(elevation))then
! !	nc1=size(elevation,dim=1)
! !	allocate(elevation1(nc1,nnode),nlayer1(nc1,nnode))
! !	elevation1=elevation(:,1:nnode)
! !	nlayer1=nlayer(:,1:nnode)
! !	deallocate(elevation,nlayer)
! !	allocate(elevation(nc1,maxnnode),nlayer(nc1,maxnnode))
! !	elevation(:,1:nnode)=elevation1
! !	nlayer(:,1:nnode)=nlayer1
! !!allocate
! !end if
! 
!!	allocate(edge(3*nnode+cln-2))	
! 
!End subroutine
