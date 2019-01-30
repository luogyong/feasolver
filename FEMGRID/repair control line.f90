  
!return the first locaton of the entry whose value equals to VAL in a integer array A of length ia
! if no entry is found, then iloc=-999999 is returned.
integer function iloc(a,ia,val)
	implicit none
	integer,intent(in)::val,ia
	integer,intent(in)::a(ia)
	integer::i
	
	iloc=-999999
	do i=1,ia
		if(a(i)==val) then
			iloc=i
			return
		end if
	end do
end function


!iflag==0, fix all the constraint lines including the outer boundary
subroutine RCL_4TH(iflag)

	use meshDS
	implicit none
	integer::i,n1,iflag,ni,n2,nc1,iloc1,n3
	intent(in)::iflag
	integer,external::iloc
	real(8)::xa,ya,xb,yb,sa,sb
	logical::tof

	tof=.true.
	nc1=0
	do while(tof)
		nc1=nc1+1
		if(nc1==10) then
			print *, 'Failed to repair Contrainted Lines'
			exit
		end if
		tof=.false.
		do i=1,ncedge
			!the edge around a hole will be kept after first fix
			!cedge().cl=999999 is an indicative that the contraint edge is generated
			! from the Prism element generation. 
			if(cedge(i).cl/=999999.and.(iflag/=0.and.csl(cedge(i).cl).hole==1)) cycle
			!the boundary will be kept after first fix.
			if(iflag/=0.and.cedge(i).cl==0) cycle

			n1=cedge(i).v(1)
			n2=cedge(i).v(2)
			n3=cedge(i).edge
			if(n3/=-1) then
				!quick judge
				if(edge(n3).v(1)==n1.and.edge(n3).v(2)==n2) cycle
				if(edge(n3).v(2)==n1.and.edge(n3).v(1)==n2) cycle
			end if
			if(.not.any(adjlist(n1).node(1:adjlist(n1).count)==n2)) then
				xa=node(n1).x
				ya=node(n1).y
				sa=node(n1).s
				xb=node(n2).x
				yb=node(n2).y
				sb=node(n2).s
				if(nnode+1>maxnnode) call EnlargeNodeRelative()
				nnode=nnode+1
				node(nnode).number=nnode
				node(nnode).x=(xa+xb)/2
				node(nnode).y=(ya+yb)/2
				node(nnode).s=(((xa-xb)**2+(ya-yb)**2)**0.5)/2
				!if Prism generation, Calculate elevation
				!if(cedge(i).cl==999999) call RCL_Prism(n1,n2)
				!add constraint edge
				cedge(i).v(2)=nnode
				cedge(i).edge=-1
				ncedge=ncedge+1
				if(ncedge>maxncedge) then
					call enlargecedge()
				end if
				cedge(ncedge).v(1)=nnode
				cedge(ncedge).v(2)=n2
				cedge(ncedge).cl=cedge(i).cl

				call TRILOC(node(nnode).x,node(nnode).y)
				call GNM_Sloan(nnode)
				tof=.true.
			else
				cedge(i).edge=adjlist(n1).edge( &
				iloc(adjlist(n1).node,adjlist(n1).count,n2))
			end if	
		end do

		if(.not.tof) print *, 'Completed successfully in & 
		repairing Constrainted Lines,nc1=',nc1
	end do

end subroutine

!!when insert a new node in the Prism generation mode
!!calculate array elevation() and nlayer() 
!subroutine RCL_Prism(ni,nj)
!	use meshds
!	implicit none
!	integer,intent(in)::ni,nj
!	
!	elevation(:,nnode)=(elevation(:,ni)+elevation(:,nj))/2.0
!	call Dsort(elevation(:,nnode),nlayer(:,nnode),nelevation)
!	
!end subroutine

  
