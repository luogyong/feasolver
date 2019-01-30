

subroutine gen_6_noded_triangle
	use meshds
	use dflib
	implicit none
	integer::i,j,k,v1,v2
	integer(4)::msg
	integer,allocatable::ninedge(:)
	
!	call gen_element_edge_set
	
	allocate(ninedge(nedge)) !mid-point number in node() of the nedgeth edge.
	ninedge=0

	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		if(nnode+1>maxnnode) call EnlargeNodeRelative()
		nnode=nnode+1
		node(nnode).number=nnode
		v1=edge(i).v(1)
		v2=edge(i).v(2)
		node(nnode).x=(node(v1).x+node(v2).x)/2.0
		node(nnode).y=(node(v1).y+node(v2).y)/2.0
		ninedge(i)=nnode
	end do
	
	do i=1,nelt
		if(elt(i).isdel.or.elt(i).et/=0) then
			cycle
		end if
		elt(i).et=6
		elt(i).node(4:6)=ninedge(elt(i).edge(1:3))
	end do

	deallocate(ninedge)
	
	msg = MODIFYMENUFLAGSQQ (4, 1, $MENUGRAYED )
	msg = MODIFYMENUFLAGSQQ (4, 2, $MENUGRAYED )

end subroutine

subroutine gen_15_noded_triangle
	use meshds
	use dflib
	implicit none
	integer::i,j,k,v1,v2
	integer(4)::msg
	integer,allocatable::ninedge(:,:)
	
!	call gen_element_edge_set

	allocate(ninedge(nedge,3)) !inner-point number in node() of the nedgeth edge.
	ninedge=0

	do i=1,nedge
		if(edge(i).v(1)*edge(i).v(2)==0) cycle
		v1=edge(i).v(1)
		v2=edge(i).v(2)
		do j=1,3
			if(nnode+1>maxnnode) call EnlargeNodeRelative()
			nnode=nnode+1
			node(nnode).number=nnode
			node(nnode).x=((4-j)*node(v1).x+j*node(v2).x)/4.0
			node(nnode).y=((4-j)*node(v1).y+j*node(v2).y)/4.0 
			ninedge(i,j)=nnode
		end do
	end do
	
	do j=1,nelt
		if(elt(j).isdel.or.elt(j).et/=0) then
			cycle
		end if
		elt(j).et=15
		!boundary nodes
		do i=1,3
			elt(j).node(3+i)=ninedge(elt(j).edge(i),2)
			if(elt(j).node(i)==edge(elt(j).edge(i)).v(1)) then				
				elt(j).node(6+2*(i-1)+1)=ninedge(elt(j).edge(i),1)
				elt(j).node(6+2*(i-1)+2)=ninedge(elt(j).edge(i),3)
			else
				elt(j).node(6+2*(i-1)+1)=ninedge(elt(j).edge(i),3)
				elt(j).node(6+2*(i-1)+2)=ninedge(elt(j).edge(i),1)
			end if
		end do
		!inner nodes
		do i=1,3
			if(nnode+1>maxnnode) call EnlargeNodeRelative()
			nnode=nnode+1
			node(nnode).x=(node(elt(j).node(3+i)).x+node(elt(j).node(4+mod(i,3))).x)/2.0
			node(nnode).y=(node(elt(j).node(3+i)).y+node(elt(j).node(4+mod(i,3))).y)/2.0
			node(nnode).number=nnode
			elt(j).node(12+i)=nnode
		end do
	end do

	deallocate(ninedge)

	msg = MODIFYMENUFLAGSQQ (4, 1, $MENUGRAYED )
	msg = MODIFYMENUFLAGSQQ (4, 2, $MENUGRAYED )
end subroutine


