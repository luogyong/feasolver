		!ATTENTION, EDGE() AND ADJLIST() HAS NOT MAINTENED IN HERE
		subroutine structuremesh()
			use meshds
			implicit none
			integer::i,j,k,icsl,v(3)=0,irow,& 
				jcol,k1,k2,nv(5)=0,n1,n2,n3,n4,n5,adn(3)=0
			real(8)::t1,t2,shafun(4)=0.0D0
			real(8),allocatable::Lx(:),Ly(:)
			integer,allocatable::elist(:)

			do i=1, nsm
				icsl=structmesh(i).csl
				do j=1,csl(icsl).num
					structmesh(i).v(j,:)=csl(icsl).conpoint(1:2,j)
				end do
				k=0
				cpp=>cpphead(icsl)
				do while(.not.associated(cpp.next, cpphead(icsl)))
					k=k+1
					do j=1,3
						if(v(j)/=0) cycle
						t1=abs(cpp.npt.x-structmesh(i).v(j,1))
						if(t1>1e-7) cycle
						t1=abs(cpp.npt.y-structmesh(i).v(j,2))
						if(t1>1e-7) cycle
						v(j)=k
					end do
					cpp=>cpp.next
				end do
				jcol=v(2)-v(1)+1
				irow=v(3)-v(2)+1
				if(allocated(lx)) deallocate(lx)
				if(allocated(ly)) deallocate(ly)
				allocate(lx(jcol),ly(irow))
				structmesh(i).nnum=irow*jcol+(irow-1)*(jcol-1)  
				!the second term is added due to that every quadrilateral is divided into four 
				!triangular element, a node must be added in the center of the quadrilateral element.

				allocate(structmesh(i).node(structmesh(i).nnum,2))
				
				k=0
				k1=1
				k2=1
				t1=0.0D0
				t2=0.0D0
				lx(1)=0.0D0
				ly(1)=0.0D0
				cpp=>cpphead(icsl)
				do while(.not.associated(cpp.next, cpphead(icsl)))
					k=k+1
					if(k>=v(1).and.k<v(2)) then
						k1=k1+1
						lx(k1)=lx(k1-1)+((cpp.npt.x-cpp.next.npt.x)**2+ (cpp.npt.y-cpp.next.npt.y)**2)**0.5
					end if
					if(k>=v(2).and.k<v(3)) then
						k2=k2+1
						ly(k2)=ly(k2-1)+((cpp.npt.x-cpp.next.npt.x)**2+ (cpp.npt.y-cpp.next.npt.y)**2)**0.5
					end if
					cpp=>cpp.next
				end do
				!scale the region into a isoparameter 4-node element. the four vertexes of the element are
				!(-1,-1),(1,-1),(1,1) and (-1,1)
				lx=lx/lx(jcol)*2-1.0
				ly=ly/ly(irow)*2-1.0
				
				!local coordinates of nodes
				k1=0
				do j=1,irow
					do k=1,jcol
						k1=k1+1
						structmesh(i).node(k1,1)=lx(k)
						structmesh(i).node(k1,2)=ly(j)											
					end do
				end do
				!the center node of each quadrilateral element
				do j=1,jcol-1
					lx(j)=(lx(j)+lx(j+1))/2.0
				end do
				do j=1,irow-1
					ly(j)=(ly(j)+ly(j+1))/2.0
				end do
				
				do j=1,irow-1
					do k=1,jcol-1
						k1=k1+1
						structmesh(i).node(k1,1)=lx(k)
						structmesh(i).node(k1,2)=ly(j)	
					end do
				end do	
				
				!globe coordinates of nodes
				do j=1,structmesh(i).nnum
					shafun(1)=(1-structmesh(i).node(j,1))*(1-structmesh(i).node(j,2))/4.0
					shafun(2)=(1+structmesh(i).node(j,1))*(1-structmesh(i).node(j,2))/4.0
					shafun(3)=(1+structmesh(i).node(j,1))*(1+structmesh(i).node(j,2))/4.0
					shafun(4)=(1-structmesh(i).node(j,1))*(1+structmesh(i).node(j,2))/4.0
					structmesh(i).node(j,:)=0.0D0
					do k=1,4
						structmesh(i).node(j,1)=structmesh(i).node(j,1)+shafun(k)*structmesh(i).v(k,1)
						structmesh(i).node(j,2)=structmesh(i).node(j,2)+shafun(k)*structmesh(i).v(k,2)
					end do
				end do
				!generate triangular element
				structmesh(i).enum=4*(irow-1)*(jcol-1)
				allocate(structmesh(i).element(structmesh(i).enum,3))
				allocate(structmesh(i).adjelem(structmesh(i).enum,3))
				n1=structmesh(i).nnum
				n2=irow*jcol
				n3=0
				n5=4*(jcol-1) !elements number in each row
				
				do j=n2+1,n1
					nv(5)=j
					nv(1)=nv(5)-n2+int(n3/n5)
					nv(2)=nv(1)+1
					nv(3)=nv(2)+jcol
					nv(4)=nv(3)-1
					n4=n3
					do k=1,4
						n3=n3+1
						structmesh(i).element(n3,1)=nv(5)
						structmesh(i).element(n3,2)=nv(k)
						structmesh(i).element(n3,3)=nv(mod(k,4)+1)
						!adajency table
						select case(k)
							case(1)
								adn(1)=n4+4
								adn(3)=n4+2
								adn(2)=n4+3-n5
								if(adn(2)<0) adn(2)=0
							case(2)
								adn(1)=n4+1
								adn(3)=n4+3
								adn(2)=n4+4+4
								if(mod(n3+2,n5)==0) adn(2)=0
							case(3)
								adn(1)=n4+2
								adn(3)=n4+4
								adn(2)=n4+1+n5
								if(adn(2)>structmesh(i).enum) adn(2)=0
							case(4)
								adn(1)=n4+3
								adn(3)=n4+1
								adn(2)=n4+2-4
								if(mod(n3-3,n5)==1) adn(2)=0																
						end select
						structmesh(i).adjelem(n3,1)=adn(1)
						structmesh(i).adjelem(n3,2)=adn(2)
						structmesh(i).adjelem(n3,3)=adn(3)
					end do	
				end do
				!add nodes and elements to the globe mesh
				if(allocated(elist)) deallocate(elist)
				allocate(elist(structmesh(i).enum))
				n1=nnode				
				do j=1,structmesh(i).nnum
					if(nnode+1>maxnnode) call EnlargeNodeRelative()
					nnode=nnode+1
					node(nnode).x=structmesh(i).node(j,1)
					node(nnode).y=structmesh(i).node(j,2)
					node(nnode).number=nnode
				end do
				ept=nelt
				do j=1,structmesh(i).enum
					if(nelt+1>maxnelement) call enlargeelement()
					nelt=nelt+1
					elist(j)=nelt
					elt(nelt).et=0
					elt(nelt).zn=structmesh(i).csl
					elt(nelt).node(1:3)=structmesh(i).element(j,1:3)+n1
				end do
				
				do j=1,structmesh(i).enum
					do k=1,3
						if(structmesh(i).adjelem(j,k)>0) &
							elt(ept+j).adj(k)=elist(structmesh(i).adjelem(j,k))
					end do
				end do
				
				deallocate(structmesh(i).node,structmesh(i).element,structmesh(i).adjelem,elist)

			end do
		end subroutine
