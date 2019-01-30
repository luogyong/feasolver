!generate a discontinuity element between to adjacent trielements
subroutine limitanalysisgrid()
		use meshds
		implicit none
		integer::i,issept(3),p1,p2,nnum1,nnum2,n1,n2
		real(8)::x1,y1,x2,y2,x3,y3
		integer::ept1
		logical,external::isacw
		
		if(nislm==0) return
		
!		open(2,file='debug.dat',status='replace')
				
		
		do ept=1,nelt
			if(elt(ept).isdel) cycle
!			write(2,'(3i15)') node(elt(ept).node(1)).number,node(elt(ept).node(2)).number,node(elt(ept).node(3)).number 
			if(elt(ept).et==0) then
				elt(ept).maxedge=0  !borrow it
				if(any(islm==elt(ept).zn)) then
					do i=1,3
						if(nnode+1>maxnnode) call EnlargeNodeRelative()
						nnode=nnode+1
						node(nnode)=node(elt(ept).node(i))
						node(nnode).number=nnode
						elt(ept).node(i)=nnode
					end do						
					if(elt(ept).node(4)>0) then
						if(nnode+1>maxnnode) call EnlargeNodeRelative()
						nnode=nnode+1
						node(nnode)=node(elt(ept).node(4))
						node(nnode).number=nnode
						elt(ept).node(4)=nnode		
					end if
				end if
			end if
!			write(2,'(3i15)') node(elt(ept).node(1)).number,node(elt(ept).node(2)).number,node(elt(ept).node(3)).number

		end do

!		close(2)

		do ept=1,nelt
			if(elt(ept).isdel) cycle
			if(elt(ept).et/=-1) then
				if(any(islm==elt(ept).zn)) then
					do i=1,3
						ept1=elt(ept).adj(i)
						if(ept1==-1) cycle						
						p1=node(elt(ept).node(i)).number
						p2=node(elt(ept).node(mod(i,3)+1)).number																

						if(elt(ept1).maxedge/=-1) then
							nelt=nelt+1
							elt(nelt).et=-1
							elt(nelt).zn=elt(ept).zn
							x1=node(p1).x
							y1=node(p1).y
							x2=node(p2).x
							y2=node(p2).y
							x3=(node(elt(ept1).node(1)).x+node(elt(ept1).node(2)).x+node(elt(ept1).node(3)).x)/3
							y3=(node(elt(ept1).node(1)).y+node(elt(ept1).node(2)).y+node(elt(ept1).node(3)).y)/3
							call fen(p1,ept1,nnum1)
							call fen(p2,ept1,nnum2)
							if(.not. isacw(x1,y1,x2,y2,x3,y3)) then
								elt(nelt).node(1)=p1
								elt(nelt).node(2)=nnum1
								elt(nelt).node(3)=nnum2
								elt(nelt).node(4)=p2																
							else
								elt(nelt).node(1)=p2
								elt(nelt).node(2)=nnum2
								elt(nelt).node(3)=nnum1
								elt(nelt).node(4)=p1
							end if
							elt(ept).maxedge=-1
						end if
					end do

				end if				
			end if
		end do		
					
end subroutine

!return the element EPT1 nodal index,NNUM, whose coordinate is node(p) .
subroutine fen(p,ept1,nnum)
	use meshds
	implicit none
	integer::i
	real(8)::x,y,t1
	integer,intent(in)::p,ept1
	integer,intent(out)::nnum
	
	x=node(p).x
	y=node(p).y
	do i=1,2
		t1=(node(elt(ept1).node(i)).x-x)**2+(node(elt(ept1).node(i)).y-y)**2
		if(abs(t1)<1e-7) then
			nnum=node(elt(ept1).node(i)).number
			return
		end if
	end do
	nnum=node(elt(ept1).node(3)).number
end subroutine

!judge whether the points (x1,y1),(x2,y2),(x3,y3) are arranged in a anticlockwise order
! yes, iscaw=.true. 
logical function isacw(x1,y1,x2,y2,x3,y3)
	implicit none
	real(8)::x1,y1,x2,y2,x3,y3,t1
	
	isacw=.false.
	y2=y2-y1
	x2=x2-x1
	y3=y3-y1
	x3=x3-x1
	t1=x2*y3-y2*x3
	if(t1>0) isacw=.true.

end function
