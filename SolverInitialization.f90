

subroutine Initialization()

!**************************************************************************************************************
!For each node, calculate the freedom number
!For each element, calcuate the element positioning vector, formulate the elastic element matrix
!Input variables: None, use those defined in solverds
!Ouput variables: None, use those defined in solverds
!Modulus Used:
!Solverds
!Subroutines Called:
!Programmer: LUO Guanyong
!Last update: 2009,10,04 
!**************************************************************************************************************
	use solverds
    USE SolverMath
    !USE MESHGEO
	implicit none
	integer::i,j,k,nj,p,j1,j2,iset1
	integer::n1,n2,n3,n4
	real(kind=DPN)::t1=0,vcos=0,vsin=0,rpi,coord1(3,4)=0,trans1(12,12)=0,c1(3,3)=0,b2(3),c2(3),R1,R2,R3,R4,G1
	real(kind=DPN)::km1(6,6)=0.d0
	integer::dof1(MNDOF)
	character(64)::ermsg=''
	
!	byte,allocatable::bar(:,:)

	rpi=pi()

	!open(2,file='fea_dug.dat',status='replace')
    
    do i=1,enum
        do k=1,ndimension
		    element(i).bbox(1,k)=MINVAL(NODE(element(i).node).COORD(K))
            element(i).bbox(2,k)=MAXVAL(NODE(element(i).node).COORD(K))
        enddo
    ENDDO
	
	do i=1,enum

        
		select case(element(i).et)
			case(CONDUCT1D) !activated dof is 4
				node(element(i).node(1:element(i).nnum)).dof(4)=0
				!initialize element stiffness matrix
				allocate(element(i).km(element(i).ndof,element(i).ndof))
				t1=0
				do j=1,3
					t1=t1+(node(element(i).node(1)).coord(j)-node(element(i).node(2)).coord(j))**2
				end do
				t1=t1**0.5
				if(abs(t1)<1e-6) then
					print *, 'element lenghth is zero. the element number is' ,i,'sub Initialization'
					stop
				end if
				element(i).km(1,1)=1/t1
				element(i).km(1,2)=-1/t1
				element(i).km(2,1)=-1/t1
				element(i).km(2,2)=1/t1
			case(UB3)
				 do j=1,2
					node(element(i).node(1:element(i).nnum)).dof(j)=0
				 end do
				 allocate(element(i).km(3,element(i).ndof))
				 p=int(material(element(i).mat).property(1))
				 allocate(element(i).a12(3,p))
				 element(i).km=0.0D0
				 do j=1,3
					element(i).km(1,2*j-1)=node(element(i).node(mod(j,3)+1)).coord(2)- &
															node(element(i).node(mod(j+1,3)+1)).coord(2)
					element(i).km(2,2*j)=node(element(i).node(mod(j+1,3)+1)).coord(1)- &
															node(element(i).node(mod(j,3)+1)).coord(1)
					element(i).km(3,2*j)=element(i).km(1,2*j-1)
					element(i).km(3,2*j-1)=element(i).km(2,2*j)
				 end do
				 element(i).property(1)=abs(element(i).km(1,1)*element(i).km(2,4)-element(i).km(2,2)*element(i).km(1,3))
				 element(i).property(3)=(node(element(i).node(1)).coord(2)+node(element(i).node(2)).coord(2)+ &
														node(element(i).node(3)).coord(2))/3.0D0								
				 element(i).km=element(i).km/element(i).property(1)
				 element(i).a12=0.0D0
				do j=1,p
					t1=2*rpi*j/p
					element(i).a12(1,j)=dcos(t1)+dsin(material(element(i).mat).property(2)/180*rpi)
					element(i).a12(2,j)=-dcos(t1)+dsin(material(element(i).mat).property(2)/180*rpi)
					element(i).a12(3,j)=2*dsin(t1)
				end do
				element(i).a12=-element(i).a12
			case(LB3)
				do j=1,3
					node(element(i).node(1:element(i).nnum)).dof(j)=0
				end do
				allocate(element(i).km(2,element(i).ndof))
				element(i).km=0.0D0
				do j=1,3
					element(i).km(1,3*j-2)=node(element(i).node(mod(j,3)+1)).coord(2)- &
											node(element(i).node(mod(j+1,3)+1)).coord(2)
					element(i).km(1,3*j)=node(element(i).node(mod(j+1,3)+1)).coord(1)- &
											node(element(i).node(mod(j,3)+1)).coord(1)
					element(i).km(2,3*j)=element(i).km(1,3*j-2)
					element(i).km(2,3*j-1)=element(i).km(1,3*j)					
				end do				
				element(i).property(1)=abs(element(i).km(1,1)*element(i).km(1,9)-element(i).km(1,7)*element(i).km(1,3))
				element(i).km=element(i).km/element(i).property(1)
				p=int(material(element(i).mat).property(1))
				allocate(element(i).a12(p,3))
				t1=material(element(i).mat).property(2)/180.0*rpi
				do j=1,p
					element(i).a12(p,1)=dcos(2*rpi*j/p)+dsin(t1)*dcos(rpi/p)
					element(i).a12(p,2)=-dcos(2*rpi*j/p)+dsin(t1)*dcos(rpi/p)
					element(i).a12(p,3)=2*dsin(2*rpi*j/p)
				end do
				
				
			case(UBZT4)
				do j=1,2
					node(element(i).node(1:element(i).nnum)).dof(j)=0
				end do
				allocate(element(i).km(4,element(i).ndof))
				allocate(element(i).a12(4,4))
				!calculate the angle between the element and the x-axial.
				t1=0
				do j=1,3
					t1=t1+(node(element(i).node(1)).coord(j)-node(element(i).node(4)).coord(j))**2
				end do
				t1=t1**0.5
				vcos=(node(element(i).node(1)).coord(1)-node(element(i).node(4)).coord(1))/t1
				vsin=(node(element(i).node(1)).coord(2)-node(element(i).node(4)).coord(2))/t1
				element(i).property(1)=t1
				 element(i).property(3)=(node(element(i).node(1)).coord(2)+node(element(i).node(4)).coord(2))/2.0D0
				element(i).km=0.0D0
				element(i).km(1,1)=-vcos
				element(i).km(1,2)=-vsin
				element(i).km(1,3)=vcos
				element(i).km(1,4)=vsin
				element(i).km(2,1)=vsin
				!element(i).km(2,1)=-vsin
				element(i).km(2,2)=-vcos
				element(i).km(2,3)=-vsin
				!element(i).km(2,3)=vsin
				element(i).km(2,4)=vcos
				element(i).km(3,7)=-vcos
				element(i).km(3,8)=-vsin
				element(i).km(3,5)=vcos
				element(i).km(3,6)=vsin
				element(i).km(4,7)=vsin
				!element(i).km(4,7)=-vsin
				element(i).km(4,8)=-vcos
				element(i).km(4,5)=-vsin
				!element(i).km(4,5)=vsin
				element(i).km(4,6)=vcos
				
				element(i).a12=0.0D0
				element(i).a12(1,1)=1
				element(i).a12(1,2)=-1
				element(i).a12(2,1)=dtan(material(element(i).mat).property(2)/180*rpi)
				element(i).a12(2,2)=element(i).a12(2,1)
				element(i).a12(3,3)=1 
				element(i).a12(3,4)=-1
				element(i).a12(4,3)=element(i).a12(2,1)
				element(i).a12(4,4)=element(i).a12(2,1)
				element(i).a12=-element(i).a12
			case(LBZT4)
				do j=1,3
					node(element(i).node(1:element(i).nnum)).dof(j)=0
				end do
				allocate(element(i).km(4,element(i).ndof))
				!calculate the angle between the element and the x-axial.
				t1=0
				do j=1,3
					t1=t1+(node(element(i).node(1)).coord(j)-node(element(i).node(4)).coord(j))**2
				end do
				t1=t1**0.5
				vcos=(node(element(i).node(4)).coord(1)-node(element(i).node(1)).coord(1))/t1
				vsin=(node(element(i).node(4)).coord(2)-node(element(i).node(1)).coord(2))/t1
				element(i).property(1)=t1
				
				element(i).km=0.0D0
				element(i).km(1,1)=vsin**2
				element(i).km(1,2)=vcos**2
				element(i).km(1,3)=-2*vsin*vcos
				element(i).km(2,1)=-0.5*element(i).km(1,3)
				element(i).km(2,2)=-element(i).km(2,1)
				element(i).km(2,3)=2*vcos**2-1
				element(i).km(1:2,4:6)=-element(i).km(1:2,1:3)
				element(i).km(3:4,7:9)=element(i).km(1:2,1:3)
				element(i).km(3:4,10:12)=-element(i).km(1:2,1:3)
			case(bar,bar2D)
				do j=1,ndimension
					node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
				end do
				!km

		
				allocate(element(i).km(2*ndimension,2*ndimension),element(i).g2l(3,3),element(i).gforce(2*ndimension),&
				element(i).gforceILS(2),element(i).Dgforce(2*ndimension))
				element(i).gforce=0.d0
				element(i).gforceILS=0.0D0
				element(i).Dgforce=0.D0

				element(i).property(1)=0.d0

				do j=1,ndimension
					t1=node(element(i).node(2)).coord(j)-node(element(i).node(1)).coord(j)															  
					element(i).property(1)=element(i).property(1)+t1**2						
				end do
				element(i).property(1)=element(i).property(1)**0.5
				
				!local system for bar elements
				! x'-axis, along the bar and the positive directioin along the increasing one of xi 
				!y'-axis, 如果x'不平行与Z轴，则y'=x'*Z, 否则 y'=x'*Y
				!z'-axis, z'=x'*y'
				! it direction cosine is stored in g2l(:,:)				
				
				!determine the direction of x'.
				do j=1,ndimension
					t1=node(element(i).node(2)).coord(j)-node(element(i).node(1)).coord(j)		
					if(abs(t1)>1e-7) then
						if(t1<0) then
							n1=1
							n2=2
						else
							n1=2
							n2=1
						end if
						exit
					end if
				end do
				
				!调整单元节点顺序，使局部坐标的x'由节点1指向节点2
				if(n1==1) then
					n3=element(i).node(1)
					element(i).node(1)=element(i).node(2)
					element(i).node(2)=n3
				end if
				
				c1=0.d0
				
				do j=1,ndimension
					c1(1,j)=(node(element(i).node(2)).coord(j)-node(element(i).node(1)).coord(j))/element(i).property(1)				
				end do
				c1(2,1:3)=0.0d0
				if(1.d0-abs(c1(1,3))>1e-7) then
					c1(2,3)=1.d0
				else
					c1(2,2)=1.d0
				end if
				do j=1,2				
					c1(3,1)=c1(1,2)*c1(2,3)-c1(1,3)*c1(2,2)
					c1(3,2)=-(c1(1,1)*c1(2,3)-c1(1,3)*c1(2,1))
					c1(3,3)=c1(1,1)*c1(2,2)-c1(1,2)*c1(2,1)
					t1=(c1(3,1)**2+c1(3,2)**2+c1(3,3)**2)**0.5
					c1(3,:)=c1(3,:)/t1
					if(j==1) c1(2,:)=c1(3,:)
				end do
				element(i).g2l=c1
				
                element(i).property(2)=element(i).property(1)				
				R1=material(element(i).mat).property(1)*material(element(i).mat).property(2)/ &
														element(i).property(1)
				element(i).property(1)=1.0D0
				
				element(i).km=0.0D0
				do j=1,ndimension	
					do k=1,ndimension	
						element(i).km(j,k)=R1*element(i).g2l(1,j)*element(i).g2l(1,k)
						element(i).km(j,k+ndimension)=-element(i).km(j,k)
						element(i).km(j+ndimension,k)=-element(i).km(j,k)
						element(i).km(j+ndimension,k+ndimension)=element(i).km(j,k)				
                    end do
                    !if(element(i).km(j,j)==0) element(i).km(j,j)=1.0D0
                end do
                
                CALL BARFAMILY_EXTREMEVALUE(I)
				
				!write(2,'(6e15.7)') (element(i).property(1)*element(i).km(j,1:6),j=1,6) passed	
			case(soilspringx,soilspringy,soilspringz,springx,springy,springz,springmx,springmy,springmz)

				
				allocate(element(i).gforceILS(1),element(i).Dgforce(1))
				
				if(.not.allocated(element(i).gforce)) allocate(element(i).gforce(1))
				if(.not.allocated(element(i).km)) allocate(element(i).km(1,1))
                
				element(i).gforce=0.d0
				element(i).gforceILS=0.0D0
				element(i).Dgforce=0.D0
				if(element(i).ec==soilspring) then
					node(element(i).node(1)).dof(element(i).et-soilspringx+1)=0
					element(i).km(1,1)=element(i).property(4)
					
				else
					if(element(i).et>=springmx) then
						node(element(i).node(1)).dof(element(i).et-springx+2)=0
					else
					    node(element(i).node(1)).dof(element(i).et-springx+1)=0
                    endif
					IF(ELEMENT(I).MAT>0) THEN
						element(i).km(1,1)=material(element(i).mat).property(1)
                        element(i).gforce(1)=material(element(i).mat).property(4)
					ELSE
						element(i).km(1,1)=element(i).property(4)
					ENDIF
				endif
				

				
			case(beam,beam2d,ssp2d)
			
				do j=1,ndimension
					node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
				end do
				if(element(i).et==beam) then
					do j=5,7
						node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
					end do
					allocate(element(i).km(12,12),element(i).gforce(12),element(i).gforceILS(12),element(i).Dgforce(12))
				else
					node(element(i).node(1:element(i).nnum)).dof(7)=0 !beam2d Mz/=0
					allocate(element(i).km(6,6),element(i).gforce(6),element(i).Dgforce(6),element(i).gforceILS(6))
				end if
				
				if(element(i).ifreedof>0) then
                    n1=freedof(element(i).ifreedof).newnode
					node(n1).dof=inactive
					node(n1).dof(freedof(element(i).ifreedof).dof)=0
                endif
                
				element(i).gforce=0.D0
				element(i).Dgforce=0.D0
				element(i).gforceILS=0.0D0
				
				element(i).property(1)=0.d0
				do j=1,ndimension
					element(i).property(1)=element(i).property(1)+(node(element(i).node(1)).coord(j)-node(element(i).node(2)).coord(j))**2 
				end do
				element(i).property(1)=element(i).property(1)**0.5	
				
				element(i).km=0.0D0	
				R1=material(element(i).mat).property(1)*material(element(i).mat).property(2)/element(i).property(1)
				R2=material(element(i).mat).property(1)*material(element(i).mat).property(5)/element(i).property(1)**3
				R3=material(element(i).mat).property(1)*material(element(i).mat).property(6)/element(i).property(1)**3
				R4=material(element(i).mat).property(1)/(2*(material(element(i).mat).property(3)+1))*material(element(i).mat).property(4)/element(i).property(1)
				!G=E/(2*(1+v))

				if(element(i).et==beam) then
					element(i).km(1,1)=R1
					element(i).km(2,2)=12*R2
					element(i).km(3,3)=12*R3
					element(i).km(4,4)=R4														
					element(i).km(5,5)=4*R3*(element(i).property(1))**2									
					element(i).km(6,6)=4*R2*(element(i).property(1))**2
					element(i).km(2,6)=6*R2*element(i).property(1)
					element(i).km(3,5)=-6*R3*element(i).property(1)
					element(i).km(1,7)=-element(i).km(1,1)
					element(i).km(2,8)=-element(i).km(2,2)
					element(i).km(2,12)=element(i).km(2,6)
					element(i).km(3,9)=-element(i).km(3,3)
					element(i).km(3,11)=element(i).km(3,5)
					element(i).km(4,10)=-element(i).km(4,4)									
					element(i).km(5,9)=-element(i).km(3,5)
					element(i).km(5,11)=2*R3*(element(i).property(1))**2
					element(i).km(6,8)=-element(i).km(2,6)
					element(i).km(6,12)=2*R2*(element(i).property(1))**2
					do j=1,6
						element(i).km(j+6,j+6)=element(i).km(j,j)					
					end do
					element(i).km(8,12)=-element(i).km(2,6)
					element(i).km(9,11)=-element(i).km(3,5)
					n1=12
				else
					element(i).km(1,1)=R1
					element(i).km(2,2)=12*R2
					element(i).km(3,3)=4*R2*element(i).property(1)**2
					do j=1,3
						element(i).km(j+3,j+3)=element(i).km(j,j)					
					end do
					element(i).km(2,3)=6*R2*element(i).property(1)
					element(i).km(1,4)=-element(i).km(1,1)
					element(i).km(2,5)=-element(i).km(2,2)
					element(i).km(3,5)=-element(i).km(2,3)
					element(i).km(2,6)=element(i).km(2,3)
					element(i).km(3,6)=2*R2*element(i).property(1)**2
					element(i).km(5,6)=-element(i).km(2,3)
					
					n1=6
				end if
				
				if(element(i).et==ssp2d) then
					!!轴偏移引起的变化项，Ref. 傅永华, 关于偏心梁单元刚度矩阵的几点说明. 力学与实践, 1996(02): 65-67.
					element(i).km(3,3)=element(i).km(3,3)+R1*material(element(i).mat).property(19)**2
					element(i).km(3,6)=element(i).km(3,6)-R1*material(element(i).mat).property(19)**2
					element(i).km(6,6)=element(i).km(3,3)
					element(i).km(1,6)=R1*material(element(i).mat).property(19)
					element(i).km(3,4)=element(i).km(1,6)
					element(i).km(1,3)=-element(i).km(1,6)
					element(i).km(4,6)=element(i).km(1,3)
				end if
				
				do j=1,n1
					do k=1,j-1
						element(i).km(j,k)=element(i).km(k,j)
					end do
				end do				
				!the local system for beam element:
				!x',along the elementand the positive directioin along the increasing one of xi 
				!for beam, y',provided by the user , obtained by element.system(2,1:3)
				!z',determined by the right-hand rule.
				!To complete the transformation maxtrix from the global system to the local system 
				!
				
				!determine the direction of x'.
				do j=1,ndimension
					t1=node(element(i).node(2)).coord(j)-node(element(i).node(1)).coord(j)		
					if(abs(t1)>1e-7) then
						if(t1<0) then
							n1=1
							n2=2
						else
							n1=2
							n2=1
						end if
						exit
					end if
				end do
				
				!调整单元节点顺序，使局部坐标的x'由节点1指向节点2
				if(n1==1) then
					n3=element(i).node(1)
					element(i).node(1)=element(i).node(2)
					element(i).node(2)=n3
				end if
				
				c1=0.d0
				c1(1,1)=(node(element(i).node(2)).coord(1)-node(element(i).node(1)).coord(1))/element(i).property(1)
				c1(1,2)=(node(element(i).node(2)).coord(2)-node(element(i).node(1)).coord(2))/element(i).property(1)
				if(element(i).et==beam) then
					c1(1,3)=(node(element(i).node(2)).coord(3)-node(element(i).node(1)).coord(3))/element(i).property(1)
					c1(2,:)=coordinate(element(i).system).c(2,:) !For BEAM, the y' direction is used only.
					c1(3,1)=c1(1,2)*c1(2,3)-c1(1,3)*c1(2,2)
					c1(3,2)=-(c1(1,1)*c1(2,3)-c1(1,3)*c1(2,1))
					c1(3,3)=c1(1,1)*c1(2,2)-c1(1,2)*c1(2,1)
					t1=(c1(3,1)**2+c1(3,2)**2+c1(3,3)**2)**0.5
					c1(3,:)=c1(3,:)/t1
					n1=4
					n2=12
				else
					c1(2,1)=-c1(1,2)
					c1(2,2)=c1(1,1)
					c1(3,3)=1.0d0
					n1=2
					n2=6
				end if
				
				!complete the local system
				allocate(element(i).g2l(3,3))
				element(i).g2l=c1
				trans1=0.0d0
				do j=1,n1					
					trans1((j-1)*3+1:j*3,(j-1)*3+1:j*3)=c1
				end do

				element(i).km=matmul(transpose(trans1(1:n2,1:n2)),matmul(element(i).km,trans1(1:n2,1:n2)))
                element(i).property(2)=element(i).property(1)
				element(i).property(1)=1.0d0
                !if(i==33) then
                !    pause
                !endif    
				!write(2,'(i7,X,<element(i).ndof>e15.7)') (i,element(i).property(1)*element(i).km(j,1:element(i).ndof),j=1,element(i).ndof) 
                
                CALL BARFAMILY_EXTREMEVALUE(I)
			
			!case(lme2d)
			!	do j=1,ndimension
			!		node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
			!	end do
			!	node(element(i).node(3)).dof(8:9)=0 !接触单元局部坐标下的接触力Nx,Ny,即拉氏乘子。
			!	allocate(element(i).km(6,6))
			!	element(i).km=0.d0
			case(pe_ssp2d) 
				node(element(i).node(1:2)).dof(1)=0
				allocate(element(i).km(2,2),element(i).gforce(2),element(i).Dgforce(2),element(i).gforceILS(2))
				element(i).gforce=0.D0
				element(i).Dgforce=0.D0
				element(i).gforceILS=0.0D0
				element(i).property(1)=1.0d0
				element(i).km(1,1)=um
				element(i).km(2,2)=um
				element(i).km(1,2)=-um
				element(i).km(2,1)=-um
				
				!每个点对点的罚单元，对应一个slave-master节点对
				!假定slave-master节点对是一一对应的，两节点之间是唯一对应的。
				do j=1, nsmnp
					if(smnp(j).master==element(i).node(1).or.smnp(j).master==element(i).node(2)) then
						smnp(j).pe=i
						element(i).ngp=j !借用ngp
						cycle
					end if				
				end do
				
				
			case(ssp2d1)  !master-slaver method，假定局部坐标与整体坐标一样，
				do j=1,ndimension
					node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
				end do
				node(element(i).node(1:element(i).nnum)).dof(7)=0 ! Mz/=0
				node(element(i).node(1:element(i).nnum)).dof(8)=0	!x' for ssp2d
				node(element(i).node(1:element(i).nnum)).dof(9)=0	!theta' for ssp2d
				allocate(element(i).km(10,10),element(i).gforce(10),element(i).Dgforce(10),element(i).gforceILS(10))

				element(i).gforce=0.D0
				element(i).Dgforce=0.D0
				element(i).gforceILS=0.0D0
				
				element(i).property(1)=0.d0
				do j=1,ndimension
					element(i).property(1)=element(i).property(1)+(node(element(i).node(1)).coord(j)-node(element(i).node(2)).coord(j))**2 
				end do
				element(i).property(1)=element(i).property(1)**0.5	
				
				element(i).km=0.0D0	
				km1=0.d0
				R1=material(element(i).mat).property(1)*material(element(i).mat).property(2)/element(i).property(1)
				R2=material(element(i).mat).property(1)*material(element(i).mat).property(5)/element(i).property(1)**3

				km1(1,1)=R1
				km1(2,2)=12*R2
				km1(3,3)=4*R2*element(i).property(1)**2
				do j=1,3
					km1(j+3,j+3)=km1(j,j)					
				end do
				km1(2,3)=6*R2*element(i).property(1)
				km1(1,4)=-km1(1,1)
				km1(2,5)=-km1(2,2)
				km1(3,5)=-km1(2,3)
				km1(2,6)=km1(2,3)
				km1(3,6)=2*R2*element(i).property(1)**2
				km1(5,6)=-km1(2,3)
				
				n1=6

				

				!!轴偏移引起的变化项，Ref. 傅永华, 关于偏心梁单元刚度矩阵的几点说明. 力学与实践, 1996(02): 65-67.
				km1(3,3)=km1(3,3)+R1*material(element(i).mat).property(19)**2
				km1(3,6)=km1(3,6)-R1*material(element(i).mat).property(19)**2
				km1(6,6)=km1(3,3)
				
				
			
				
				do j=1,n1
					do k=1,j-1
						km1(j,k)=km1(k,j)
					end do
				end do	
				
				
				do j=1,2
					if(j==1) then
						dof1(1:6)=(/1,2,3,6,7,8/)
						km1(1,6)=R1*material(element(i).mat).property(19)
					else
						dof1(1:6)=(/4,2,5,9,7,10/)
						km1(1,6)=-km1(1,6)
					end if
					km1(3,4)=km1(1,6)
					km1(1,3)=-km1(1,6)
					km1(4,6)=km1(1,3)
					!对称项
					km1(6,1)=km1(1,6)
					km1(4,3)=km1(3,4)
					km1(3,1)=km1(1,3)
					km1(6,4)=km1(4,6)
					
					
					do k=1,6
						do j1=1,6
							element(i).km(dof1(k),dof1(j1))=element(i).km(dof1(k),dof1(j1))+km1(k,j1)
						end do
					end do
					
				end do
				
				!!初始化为两梁固定
				!ELEMENT(I).KM(1,1)=ELEMENT(I).KM(1,1)+UM
				!ELEMENT(I).KM(4,4)=ELEMENT(I).KM(4,4)+UM				
				!ELEMENT(I).KM(1,4)=ELEMENT(I).KM(1,4)+UM
				!ELEMENT(I).KM(4,1)=ELEMENT(I).KM(1,4)
				
				!ELEMENT(I).KM(6,6)=ELEMENT(I).KM(6,6)+UM
				!ELEMENT(I).KM(9,9)=ELEMENT(I).KM(9,9)+UM				
				!ELEMENT(I).KM(6,9)=ELEMENT(I).KM(6,9)+UM
				!ELEMENT(I).KM(9,6)=ELEMENT(I).KM(6,9)				
				
				!the local system for beam element:
				!x',along the elementand the positive directioin along the increasing one of xi 
				!for beam, y',provided by the user , obtained by element.system(2,1:3)
				!z',determined by the right-hand rule.
				!To complete the transformation maxtrix from the global system to the local system 
				!
				
				!determine the direction of x'.
				do j=1,ndimension
					t1=node(element(i).node(2)).coord(j)-node(element(i).node(1)).coord(j)		
					if(abs(t1)>1e-7) then
						if(t1<0) then
							n1=1
							n2=2
						else
							n1=2
							n2=1
						end if
						exit
					end if
				end do
				
				!调整单元节点顺序，使局部坐标的x'由节点1指向节点2
				if(n1==1) then
					n3=element(i).node(1)
					element(i).node(1)=element(i).node(2)
					element(i).node(2)=n3
				end if
				
				c1=0.d0
				c1(1,1)=(node(element(i).node(2)).coord(1)-node(element(i).node(1)).coord(1))/element(i).property(1)
				c1(1,2)=(node(element(i).node(2)).coord(2)-node(element(i).node(1)).coord(2))/element(i).property(1)
				c1(2,1)=-c1(1,2)
				c1(2,2)=c1(1,1)
				c1(3,3)=1.0d0
				n1=2
				n2=6

				
				!complete the local system
				allocate(element(i).g2l(3,3))
				element(i).g2l=c1
				!trans1=0.0d0
				!do j=1,n1					
				!	trans1((j-1)*3+1:j*3,(j-1)*3+1:j*3)=c1
				!end do

				!element(i).km=matmul(transpose(trans1(1:n2,1:n2)),matmul(element(i).km,trans1(1:n2,1:n2)))
                !element(i).property(2)=element(i).property(1)
				!element(i).property(1)=1.0d0
				!!write(2,'(12e15.7)') (element(i).property(1)*element(i).km(j,1:12),j=1,12) 
                
                CALL BARFAMILY_EXTREMEVALUE(I)					
				
				
				
				
			case(shell3)
				do j=1,3
					node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
				end do
				do j=5,7
					node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
				end do
				!h, shell thickness
				allocate(element(i).km(18,18),element(i).gforce(18),element(i).Dgforce(18),element(i).gforceILS(18), &
							element(i).g2l(3,3))

				
				element(i).km=0.0d0
				element(i).gforce=0.0d0
				element(i).Dgforce=0.0d0
				element(i).gforceILS=0.0d0
				element(i).g2l=0.0d0
				call SM3SHELL (node(element(i).node(1:3)).coord(1), &
				 node(element(i).node(1:3)).coord(2), &
				 node(element(i).node(1:3)).coord(3), element(i).property(3), &
				 material(element(i).mat).property(1), &
				 material(element(i).mat).property(2), &
				 element(i).km,element(i).ndof,element(i).g2l,ermsg, &
				 element(i).property(2))
				if(len_trim(ermsg)/=0) then
					print *, ermsg(1:len_trim(ermsg))//'Errors in shell element(i), i=  ',i					 
					stop
				end if
				 element(i).property(1)=1.0

				!redefine the local system as follows:(Referred to ABAQUS)
				!The x-axis is projected on to the surface in question.This line represents the local 1-direction. 
				!If the global x-axis is within 0.1 degrees of the normal to the surface, 
				!the local 1-direction is the projection of the global z-axis onto the surface. 
				!The local 2-direction is the at right angles to the local 1-direction, 
				!so that the local 1-direction, local 2-direction and the positive normal to the surface 
				!form a right-handed set. 
				!The positive normal direction is defined in an element 
				!by the right-handed rotation rule going around the nodes of the element 
				!(in the anticlockwise direction).
				!
				t1=0.1
				if(abs(element(i).g2l(3,1))<Dcosd(t1)) then
					element(i).g2l(2,1)=0
					element(i).g2l(2,2)=element(i).g2l(3,3)
					element(i).g2l(2,3)=-element(i).g2l(3,2)	
				else
					element(i).g2l(2,1)=element(i).g2l(3,2)
					element(i).g2l(2,2)=-element(i).g2l(3,1)
					element(i).g2l(2,3)=0					
				end if
				t1=(element(i).g2l(2,1)**2+element(i).g2l(2,2)**2+element(i).g2l(2,3)**2)**0.5
				element(i).g2l(2,:)=element(i).g2l(2,:)/t1
				element(i).g2l(1,1)=element(i).g2l(2,2)*element(i).g2l(3,3)-element(i).g2l(3,2)*element(i).g2l(2,3)
				element(i).g2l(1,2)=-(element(i).g2l(2,1)*element(i).g2l(3,3)-element(i).g2l(3,1)*element(i).g2l(2,3))
				element(i).g2l(1,3)=element(i).g2l(2,1)*element(i).g2l(3,2)-element(i).g2l(3,1)*element(i).g2l(2,2)
				t1=(element(i).g2l(1,1)**2+element(i).g2l(1,2)**2+element(i).g2l(1,3)**2)**0.5
				element(i).g2l(1,:)=element(i).g2l(1,:)/t1 
                
                CALL BARFAMILY_EXTREMEVALUE(I)
			case(shell3_kjb)
			
			case(pipe2,wellbore,WELLBORE_SPGFACE)
				node(element(i).node(1:element(i).nnum)).dof(4)=0
				allocate(element(i).km(element(i).nnum,element(i).nnum))
				element(i).km=0.0D0
                IF(element(i).et/=WELLBORE_SPGFACE) THEN
                    t1=vk(material(element(i).mat).property(2))/9.8
                    if(abs(material(element(i).mat).property(4)-1.d0)>1.d-7) then
                        t1=t1/24./3600. !UNIT,M.DAY
                    endif
                    element(i).property(1)=(2*material(element(i).mat).property(1))**2/t1/32
                    T1=RPI*material(element(i).mat).property(1)**2/ &
                        NORM2(NODE(ELEMENT(I).NODE(1)).COORD-NODE(ELEMENT(I).NODE(2)).COORD)
                               
                    ELEMENT(I).PROPERTY(1)=1./(T1*ELEMENT(I).PROPERTY(1))
                ELSE
                    ELEMENT(I).PROPERTY(1)=1.D10 !模拟不透水,大阻力
                ENDIF
                
                if(element(i).et==WELLBORE.or.element(i).et==WELLBORE_SPGFACE) then
                    CALL INI_WELLBORE(I)                   
                ELSE
                    !ALLOCATE(ELEMENT(I).KM(2,2))
                    ELEMENT(I).KM=1.D0
                    ELEMENT(I).KM(1,2)=-1.D0;ELEMENT(I).KM(2,1)=-1.D0
                    ELEMENT(I).KM=ELEMENT(I).KM/ELEMENT(I).PROPERTY(1)
                endif
                ISOUT_WELL_FILE=.TRUE.
                
            CASE(SEMI_SPHFLOW,SPHFLOW)
                node(element(i).node(1:element(i).nnum)).dof(4)=0
				allocate(element(i).km(element(i).nnum,element(i).nnum))
				element(i).km=0.0D0
                CALL INI_SPHFLOW(I)
                ISOUT_WELL_FILE=.TRUE.
                
			CASE(ZT4_SPG,ZT6_SPG) !assume no flow occurs in the diections parallel to element faces.
				node(element(i).node(1:element(i).nnum)).dof(4)=0
                
               
                CALL ZT_SPG_INI2(I)
				
				
			case default
				select case(element(i).ec)
					case(spg2d,spg,CAX_SPG)
						node(element(i).node(1:element(i).nnum)).dof(4)=0
					case(cpl,CAX_CPL)
						do j=1,NDIMENSION
							node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
						end do
						node(element(i).node(1:element(i).nnum)).dof(4)=0
					case(c3d,cpe,cps,cax)
						do j=1,NDIMENSION
							node(element(i).node(1:element(i).nnum)).dof(j)=0 !ACTIVATIVE THE NODAL DOF
						end do						
																	
                end select
				
                IF(ELEMENT(I).ET==ZT4_SPG2.OR.ELEMENT(I).ET==ZT6_SPG2) CALL ZT_SPG_INI2(I)    
                    
               
               
                
				!call Calangle(i)
				!allocate room for km,b,d.
				call el_alloc_room(i)
				!according to the material,initialize d
				call el_ini_D(i)
				!calculate the globe co-ordinates for the gassian point
				call xygp(i)				
				!n1=element(i).ngp+element(i).nnum				
				n1=element(i).ngp
				call JACOB2(i,element(i).b,element(i).nd,element(i).ndof,n1, &
									element(i).detjac,element(i).ngp)
				call BTDB(element(i).b,element(i).nd,element(i).ndof,element(i).ngp, &
								element(i).d,element(i).nd,element(i).nd, &
								element(i).km,element(i).ndof,element(i).ndof, &
								ecp(element(i).et).weight,ecp(element(i).et).ngp, &
								element(i).detjac,element(i).ngp,i)
                !if(isnan(element(i).km))
				
!				write(99,*) 'element,', i
!				write(99,999) ((element(i).km(j,k),k=1,element(i).ndof),j=1,element(i).ndof)
!				coord(:,1)=node(element(i).node).coord(1)				
!				coord(:,2)=node(element(i).node).coord(2)
!				element(i).km=0.0D0
!				call stiff4(element(i).km,coord,material(element(i).mat).property(1),material(element(i).mat).property(2))											
!				write(2,'(/,a)') ''
!				write(2,999) ((-element(i).km(j,k),k=1,element(i).ndof),j=1,element(i).ndof)

!				999 FORMAT(<ELEMENT(I).NDOF>E15.7)
		
				


		end select
	end do

	! calculate the dof number of each node,and number it
	! set ndof value
	
	do i=1,nSMNP
		node(smnp(i).slave).dof(smnp(i).sdof)=-i
	end do
	
	n1=size(node(1).dof)
	do i=1,nnum
		do j=1,n1
			if(node(i).dof(j)==inactive) cycle
			
			if(node(i).dof(j)==0) then
				node(i).ndof=node(i).ndof+1
				ndof=ndof+1
				node(i).dof(j)=ndof
			end if
			if(node(i).dof(j)<0) then
				n2=-node(i).dof(j)
				if(node(smnp(n2).master).dof(smnp(n2).mdof)>0) then !如果主节点mdof已经编号
					node(i).ndof=node(i).ndof+1
					node(i).dof(j)=node(smnp(n2).master).dof(smnp(n2).mdof)
				else
					node(i).ndof=node(i).ndof+1
					ndof=ndof+1
					node(i).dof(j)=ndof
					
					node(smnp(n2).master).dof(smnp(n2).mdof)=ndof
					node(smnp(n2).master).ndof=node(smnp(n2).master).ndof+1
				end if
			end if
		end do
	end do
	nuvar=ndof
	
	do i=1,nfreedof
		dof1(freedof(i).dof)=node(freedof(i).newnode).dof(freedof(i).dof)
		node(freedof(i).newnode).dof=node(freedof(i).node).dof
        node(freedof(i).newnode).ndof=node(freedof(i).node).ndof
		node(freedof(i).newnode).dof(freedof(i).dof)=dof1(freedof(i).dof)
	enddo
	
	!Initialize the element positioning vector
	!initialize the band width of each dof
	allocate(bw(ndof),Lmre(ndof))
	!bw=0
	!LDLUBW=0 
	bw=10000000 !location of the most LEFT entry
	Lmre=0 !location of the most RIGHT entry
!	if(.not.solver_control.issym) then
!		allocate(ubw(ndof))
!		ubw=0
!	end if

	!initialize the element position vector.	
	do i=1,enum
		select case(element(i).et)
			case(CONDUCT1D)				
				dof1=0
				dof1(1)=4
				call fepv(i,dof1)
				call dofbw(i)

			case(UB3)
				p=int(material(element(i).mat).property(1))
				element(i).ndof=element(i).ndof+p
				ncons=ncons+3
				nnz=nnz+12+3*p
				!allocate(element(i).g(element(i).ndof))
				dof1=0
				dof1(1)=1
				dof1(2)=2
				call fepv(i,dof1)
				n1=element(i).ndof-p					
				do j=1,p
					ndof=ndof+1
					n1=n1+1
					element(i).g(n1)=ndof
				end do
			case(UBZT4)
				element(i).ndof=element(i).ndof+4
				ncons=ncons+4
				nnz=nnz+24
				!allocate(element(i).g(element(i).ndof))
				dof1=0
				dof1(1)=1
				dof1(2)=2
				call fepv(i,dof1)
				n1=element(i).ndof-4				
				do j=1,4
					ndof=ndof+1
					n1=n1+1
					element(i).g(n1)=ndof
				end do
			case(CPE3,CPE6,CPE4,CPE8,CPE4R,CPE8R,CPE15, &
						CPS3,CPS4,CPS4R,CPS8,CPS8R,CPS6,CPS15, &
						CAX3,CAX4,CAX4R,CAX6,CAX15,CAX8,CAX8R,BAR2D)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
			case(PRM6,PRM15,BAR)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				dof1(3)=3
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)				
			case(CPE3_SPG,CPE6_SPG,CPE4_SPG,CPE8_SPG,CPE4R_SPG,CPE8R_SPG,CPE15_SPG, &
					 CPS4_SPG,CPS4R_SPG,CPS8_SPG,CPS8R_SPG,CPS6_SPG,CPS15_SPG, &	
					 CAX3_SPG,CAX4_SPG,CAX4R_SPG,CAX6_SPG,CAX15_SPG,CAX8_SPG,CAX8R_SPG, &
					 PRM6_SPG,PRM15_SPG,TET4_SPG,TET10_SPG,ZT4_SPG2,ZT6_SPG2)
				dof1=0
				dof1(4)=4
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)

				call Cal_epsilon(i)
!				!epslong one(low)
!				element(i).property(3)=minval(element(i).xygp(ndimension,:))- &
!						minval(node(element(i).node).coord(ndimension))
!				!epslong two(up)
!				element(i).property(2)=maxval(node(element(i).node).coord(ndimension))- &
!										maxval(element(i).xygp(ndimension,:))
						

			case(CPE3_CPL,CPE6_CPL,CPE4_CPL,CPE8_CPL,CPE4R_CPL,CPE8R_CPL,CPE15_CPL, &
					 CPS4_CPL,CPS4R_CPL,CPS8_CPL,CPS8R_CPL,CPS6_CPL,CPS15_CPL, &	
					 CAX4_CPL,CAX4R_CPL,CAX6_CPL,CAX15_CPL,CAX8_CPL,CAX8R_CPL)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				dof1(4)=4
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
			case(PRM6_CPL,PRM15_CPL,TET4_CPL,TET10_CPL)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				dof1(3)=3
				dof1(4)=4
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
			case(BEAM,shell3,shell3_kjb)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				dof1(3)=3
				dof1(5)=5
				dof1(6)=6
				dof1(7)=7
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
			case(BEAM2D,SSP2D)
				dof1=0
				dof1(1)=1
				dof1(2)=2
				dof1(7)=7
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
                !write(2,"(i5,<element(i).ndof>i5)") element(i).g
			case(ssp2d1)
				dof1=0
				dof1(1:2)=1
				dof1(7)=7
				dof1(8)=8
				dof1(9)=9
				call fepv(i,dof1)
				call dofbw(i)
				
			case(dkt3)
				dof1=0
				dof1(3)=3
				dof1(5)=5
				dof1(6)=6
				!allocate(element(i).g(element(i).ndof))
				call fepv(i,dof1)
				call dofbw(i)
				
			case(pipe2,wellbore,WELLBORE_SPGFACE,SEMI_SPHFLOW,SPHFLOW,ZT4_SPG,ZT6_SPG)
				dof1=0
				dof1(4)=4
				call fepv(i,dof1)
				call dofbw(i)
			case(soilspringx,soilspringy,soilspringz,springx,springy,springz,springmx,springmy,springmz)
				dof1=0
				if(element(i).ec==soilspring) then
					dof1(element(i).et-soilspringx+1)=element(i).et-soilspringx+1	
				else
					if(element(i).et>=springmx) then
						dof1(element(i).et-springx+2)=element(i).et-springx+2
					else
						dof1(element(i).et-springx+1)=element(i).et-springx+1	
					endif
				endif
				call fepv(i,dof1)
				call dofbw(i)
		end select
		
	end do
	
	!if(.not.solver_control.issym) then
	bw=lmre-bw+1 !BAND WIDTH 
	!end if
	allocate(load(ndof),tdisp(ndof),adof(ndof))
	adof=0
	load=0.0D0
	tdisp=0.0D0
	if(.not.allocated(Tstepdis)) then
		allocate(Tstepdis(ndof,0:nstep))
		Tstepdis=0.D0
	end if
	!默认读入的初值是0步的结果。
	do i=1,NiniV
		Tstepdis(node(inivalue(i).node).dof(inivalue(i).dof),0)=inivalue(i).value
	end do


!	if(.not.stepinfo(1).issteady) call cmm_spg_cal() 
	
	
	if(solver_control.solver==LA04 &
		.or.solver_control.solver==LPSOLVER &
		.or.solver_control.solver==MOSEK)  return
	
	bwmax=maxval(bw)
	!diaglkmloc(1)=1
	do i=2,ndof
		bw(i)=bw(i)+bw(i-1) !bw(i) is NOW the location of the most right entry in the i row in the total matrix.
		!if(.not.solver_control.ismkl) diaglkmloc(i)=bw(i) !!!1!对角元素的在总刚中位置,这时以总刚以下三角的形式存储
	end do

	! allocate space for load() and km
	!nnz=bw(ndof)
    !
	!if(.not.solver_control.issym) nnz=UBW(NDOF)

	load=0.0D0


	call irowjcol_SEC()
	if(.not.allocated(km)) allocate(km(nnz))
	km=0.d0


	
	!convert body force to nodal force
	if(bfnum>0) call cbf2nf()
	
	!if(solver_control.i2ncal==spr) then
	!	!setting environments for nodal stress calculation using SPR method.
	!	call spr_elist()
	!	!use SPR method to calculate the polynomials for each patch.
	!	call spr_polynomial_cal()
	!end if
	
	
	!form element list
	!call FormElementList()

	
    
	!INITIALIZED DOFHEAD
	CALL NDOF_HEAD()
	
	!generate brickelement information for bar and beam family element to post-processing.
	do i=1,neset
        iset1=esetid(i)
		if(eset(iset1).ec==stru) then
			if(eset(iset1).et==bar.or.eset(iset1).et==bar2d.or.eset(iset1).et==beam2d.or.eset(iset1).et==beam.or.eset(iset1).et==ssp2d)	call generate_brickelement_bar_beam(iset1)
		end if
    end do
	!对于临空的渗流薄元节点，设成出溢边界
    CALL ZT_SPFACE_BC()
	
end subroutine




! form the element(ienum) position vector
subroutine fepv(ienum,dof1)
	use solverds
	implicit none
	integer::ienum,dof1(*)
	integer::j,n1,n2
	
	allocate(element(ienum).g(element(ienum).ndof))
	n1=1
	do j=1,element(ienum).nnum
		n2=1
		do while(n2<=MNDOF)
			if(dof1(n2)>0) then
			 	element(ienum).g(n1)=node(element(ienum).node(j)).dof(dof1(n2))
				n1=n1+1			
			end if
			n2=n2+1	
		end do	
	end do
end subroutine

!according element connection, calculate the bandwidth.
subroutine dofbw(ienum)
	use solverds
    USE quicksort
	implicit none
	integer::ienum,n1,n2,n3,n4,j,DOF1(200),I
    

	!n1=minval(element(ienum).g)
	!n3=maxval(element(ienum).g)
    DOF1(1:element(ienum).ndof)=element(ienum).g
    
    CALL quick_sort(DOF1(1:element(ienum).ndof))
    n1=dof1(1);n3=dof1(element(ienum).ndof)
	do j=1,element(ienum).ndof
		
		
		if(solver_control.issym) then
			
			if(.not.solver_control.ismkl) then
				!for default solver,存下三角,bw(i)为第i行下角形带宽。
!				n2=dof1(j)-n1+1
!				if(bw(dof1(j))<n2) bw(dof1(j))=n2
				
				if(bw(dof1(j))>n1) bw(dof1(j))=n1 !location of the most Left entry
				Lmre(dof1(j))=dof1(j) !location of the most right entry
				
                DO I=1,J
                    CALL SETUP_DOFADJL(DOF1(I),DOF1(J))
                ENDDO
                
			else
				!for mkl solver, 存上三角，bw(i)为第i行上角形带宽
!				n2=n3-dof1(j)+1
!				if(bw(dof1(j))<n2) bw(dof1(j))=n2
				
				bw(dof1(j))=dof1(j) !location of the most Left entry
				if(Lmre(dof1(j))<n3) Lmre(dof1(j))=n3 !location of the most right entry			
				 DO I=J,element(ienum).ndof
                    CALL SETUP_DOFADJL(DOF1(I),DOF1(J))
                 ENDDO
			end if

		else
			if(bw(dof1(j))>n1) bw(dof1(j))=n1 !location of the most Left entry
			if(Lmre(dof1(j))<n3) Lmre(dof1(j))=n3 !location of the most right entry	
            DO I=1,element(ienum).ndof
                CALL SETUP_DOFADJL(DOF1(I),DOF1(J))                
            ENDDO
		end if
		
		
		

	end do

end subroutine


SUBROUTINE SETUP_DOFADJL(ITEM,IDOF)
    USE solverds
    USE MESHGEO,ONLY:I_ENLARGE_AR
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ITEM,IDOF
    INTEGER::I,N1
    
    IF(.NOT.ALLOCATED(DOFADJL)) THEN
        ALLOCATE(DOFADJL(NDOF))
        DO I=1,NDOF
            ALLOCATE(DOFADJL(I).DOF(DOFADJL(I).DOF_SIZE))    
            DOFADJL(I).NDOF=1
            DOFADJL(I).DOF(1)=I
        ENDDO    
    ENDIF

    IF(.NOT.ANY(DOFADJL(IDOF).DOF(1:DOFADJL(IDOF).NDOF)-ITEM==0)) THEN
            
        IF(DOFADJL(IDOF).NDOF+1>DOFADJL(IDOF).DOF_SIZE) THEN
            CALL I_ENLARGE_AR(DOFADJL(IDOF).DOF,10)
            DOFADJL(IDOF).DOF_SIZE=DOFADJL(IDOF).DOF_SIZE+10
        ENDIF
        IF(ITEM>DOFADJL(IDOF).DOF(DOFADJL(IDOF).NDOF)) THEN
            N1=DOFADJL(IDOF).NDOF+1
        ELSEIF(ITEM<DOFADJL(IDOF).DOF(1)) THEN
            N1=1
            DOFADJL(IDOF).DOF(DOFADJL(IDOF).NDOF+1:2:-1)=DOFADJL(IDOF).DOF(DOFADJL(IDOF).NDOF:1:-1)
        ELSE
            N1=MINVAL(PACK([1:DOFADJL(IDOF).NDOF],DOFADJL(IDOF).DOF(:DOFADJL(IDOF).NDOF)-ITEM>0))
            DOFADJL(IDOF).DOF(DOFADJL(IDOF).NDOF+1:N1+1:-1)=DOFADJL(IDOF).DOF(DOFADJL(IDOF).NDOF:N1:-1)
        ENDIF
        DOFADJL(IDOF).DOF(N1)=ITEM
        DOFADJL(IDOF).NDOF=DOFADJL(IDOF).NDOF+1
    ENDIF
    
    

    ENDSUBROUTINE
    
    

! allocate room for element.km,b,d.
subroutine el_alloc_room(ienum)
	use solverds
	implicit none
	integer::i,ienum
	
	allocate(element(ienum).km(element(ienum).ndof,element(ienum).ndof))
	element(ienum).km=0.0D0
	allocate(element(ienum).b(element(ienum).nd,element(ienum).ndof,element(ienum).ngp))
	element(ienum).b=0.0D0
	allocate(element(ienum).d(element(ienum).nd,element(ienum).nd))
	allocate(element(ienum).detjac(element(ienum).ngp))
	element(ienum).detjac=0.0D0
	allocate(element(ienum).xygp(3,element(ienum).ngp))
	element(ienum).xygp=0.0
	select case(element(ienum).ec)
		case(spg2d,spg,CAX_spg)
			allocate(element(ienum).igrad(element(ienum).nd,element(ienum).ngp+element(ienum).nnum))
			element(ienum).igrad=0.0D0
			allocate(element(ienum).velocity(element(ienum).nd,element(ienum).ngp+element(ienum).nnum))
			element(ienum).velocity=0.0D0
			allocate(element(ienum).flux(element(ienum).nnum))
			element(ienum).flux=0.0D0
			allocate(element(ienum).Kr(element(ienum).ngp+element(ienum).nnum),element(ienum).Mw(element(ienum).ngp+element(ienum).nnum))
			element(ienum).kr=0.0D0
			element(ienum).Mw=0.0D0
			allocate(element(ienum).sita_ini(element(ienum).ngp),element(ienum).sita_fin(element(ienum).ngp))
			element(ienum).sita_ini=0.0d0
			element(ienum).sita_fin=0.0d0
            
			!allocate(element(ienum).lamda(element(ienum).ngp))
			!element(ienum).lamda=1.0
			do i=1,element(ienum).nnum
				if(.not.allocated(node(element(ienum).node(i)).igrad)) then
					allocate(node(element(ienum).node(i)).igrad(ndimension))
					allocate(node(element(ienum).node(i)).velocity(ndimension))
				end if				
			end do
			!if(.not.stepinfo(1).issteady) then
			!	allocate(element(ienum).cmm(element(ienum).nnum,element(ienum).nnum))
			!	element(ienum).cmm=0.0D0
			!end if
		case(CPE,CPS,c3d,CPL,CAX,CAX_CPL)
			allocate(element(ienum).stress(6,element(ienum).ngp+element(ienum).nnum))
			allocate(element(ienum).Dstress(6,element(ienum).ngp+element(ienum).nnum))
			element(ienum).stress=0.0D0
			element(ienum).Dstress=0.0D0
			allocate(element(ienum).strain(6,element(ienum).ngp+element(ienum).nnum))
			allocate(element(ienum).Dstrain(6,element(ienum).ngp+element(ienum).nnum))
			element(ienum).strain=0.0D0
			element(ienum).Dstrain=0.0D0
			allocate(element(ienum).pstrain(6,element(ienum).ngp+element(ienum).nnum))
			allocate(element(ienum).evp(6,element(ienum).ngp+element(ienum).nnum))	
			allocate(element(ienum).ev(element(ienum).ngp), &
			element(ienum).e(element(ienum).ngp),element(ienum).uw(element(ienum).ngp))
			allocate(element(ienum).sfr(8,element(ienum).ngp+element(ienum).nnum))	
			element(ienum).pstrain=0.0D0
			element(ienum).evp=0.0D0	
			element(ienum).ev=0.0
			element(ienum).e=0.0
            element(ienum).uw=0.0
			element(ienum).sfr=0.d0
			do i=1,element(ienum).nnum
				if(.not.allocated(node(element(ienum).node(i)).stress)) then
					allocate(node(element(ienum).node(i)).stress(6))
					allocate(node(element(ienum).node(i)).strain(6))
					allocate(node(element(ienum).node(i)).pstrain(6))
					IF(OUTVAR(SFR).VALUE>0) then
                        allocate(node(element(ienum).node(i)).sfr(11))
                        !node(element(ienum).node(i)).sfr(9)=SOLVER_CONTROL.slidedirection
                    endif
					IF(OUTVAR(PSIGMA).VALUE>0) allocate(node(element(ienum).node(i)).PSIGMA(4))
				end if				
			end do			
	end select
	

end subroutine

!according the element material,initialize the elemental strain-stress matrix
subroutine el_ini_D(ienum)
	use solverds
    use SolverMath
	implicit none
	integer::ienum
	integer::nd1,mat1,isys1,i
	real(8)::t1=0,t2=0,VSIN1,VCOS1,V1(3,3)
	
	nd1=element(ienum).nd
	mat1=element(ienum).mat
	select case(element(ienum).ec)
		case(cps)
			element(ienum).d=0.0D0
			element(ienum).d(1,1)=1.0
			element(ienum).d(2,2)=1.0
			t1=material(mat1).property(2)
			t2=material(mat1).property(1)
			element(ienum).d(4,4)=(1-t1)/2.0
			element(ienum).d(1,2)=t1
			element(ienum).d(2,1)=element(ienum).d(1,2)
			element(ienum).d=element(ienum).d*material(mat1).property(1) &
				/(1-material(mat1).property(2)**2)
		case(cpe,cax,c3d)	
			element(ienum).d=0.0D0
			element(ienum).d(1,1)=1.0
			element(ienum).d(2,2)=1.0
			element(ienum).d(3,3)=1.0
			t1=material(mat1).property(2)
			t2=material(mat1).property(1)
			do i=4,2*ndimension
				element(ienum).d(i,i)=(1.0-2*t1)/(2*(1.0-t1))
			end	do
			t1=t1/(1-t1)
			element(ienum).d(2,1)=t1
			element(ienum).d(1,2)=t1
			element(ienum).d(1,3)=t1
			element(ienum).d(2,3)=t1
			element(ienum).d(3,1)=t1
			element(ienum).d(3,2)=t1
			t1=material(mat1).property(1)* (1-material(mat1).property(2))/ &
				((1+material(mat1).property(2))* (1-2*material(mat1).property(2)))
			element(ienum).d=element(ienum).d*t1
		case(spg2d,CAX_spg)
			element(ienum).d=0.0D0
			element(ienum).d(1,1)=material(mat1).property(1)
			element(ienum).d(2,2)=material(mat1).property(2)
			isys1=int(material(mat1).property(4))            
            IF(element(ienum).ET==ZT4_SPG2.AND.(element(ienum).d(1,1)-element(ienum).d(2,2))>.1D-7) THEN
                ISYS1=-1
                T1=NORM2(GNODE(:,ELEMENT(IENUM).NODE2(4))-GNODE(:,ELEMENT(IENUM).NODE2(1)))
                VCOS1=GNODE(1,ELEMENT(IENUM).NODE2(4))-GNODE(1,ELEMENT(IENUM).NODE2(1))/T1
                VSIN1=GNODE(2,ELEMENT(IENUM).NODE2(4))-GNODE(2,ELEMENT(IENUM).NODE2(1))/T1
                coordinate(isys1).c(1,1)=VCOS1;coordinate(isys1).c(2,2)=VCOS1;
                coordinate(isys1).c(2,1)=VSIN1;coordinate(isys1).c(1,2)=-VSIN1;
            ENDIF
			if(isys1/=0.AND.ABS(element(ienum).d(1,1)-element(ienum).d(2,2))>.1D-7) then
                
				element(ienum).d(1:nd1,1:nd1)=matmul(matmul(coordinate(isys1).c(1:nd1,1:nd1), &
													element(ienum).d(1:nd1,1:nd1)),transpose(coordinate(isys1).c(1:nd1,1:nd1)))
			end if
		case(spg)
			element(ienum).d=0.0D0
			element(ienum).d(1,1)=material(mat1).property(1)
			element(ienum).d(2,2)=material(mat1).property(2)
			element(ienum).d(3,3)=material(mat1).property(3)
			isys1=int(material(mat1).property(4))
            IF(element(ienum).ET==ZT6_SPG2.AND.(element(ienum).d(1,1)-element(ienum).d(3,3))>.1D-7) THEN
                ISYS1=-1
                T1=NORM2(GNODE(:,ELEMENT(IENUM).NODE2(4))-GNODE(:,ELEMENT(IENUM).NODE2(1)))
                coordinate(isys1).c(1,:)=GNODE(:,ELEMENT(IENUM).NODE2(4))-GNODE(:,ELEMENT(IENUM).NODE2(1))/T1
                T1=NORM2(GNODE(:,ELEMENT(IENUM).NODE2(2))-GNODE(:,ELEMENT(IENUM).NODE2(1)))
                coordinate(isys1).c(2,:)=GNODE(:,ELEMENT(IENUM).NODE2(2))-GNODE(:,ELEMENT(IENUM).NODE2(1))/T1
                v1(:,1)=coordinate(isys1).c(1,:);v1(:,2)=coordinate(isys1).c(2,:);v1(:,3)=0.d0
                coordinate(isys1).c(3,:)=NORMAL_TRIFACE(v1)
                coordinate(isys1).c(3,:)=coordinate(isys1).c(3,:)/norm2(coordinate(isys1).c(3,:))
                coordinate(isys1).c=transpose(coordinate(isys1).c)
            ENDIF            
			if(isys1/=0) then
				element(ienum).d(1:nd1,1:nd1)=matmul(matmul(coordinate(isys1).c(1:nd1,1:nd1), &
													element(ienum).d(1:nd1,1:nd1)),transpose(coordinate(isys1).c(1:nd1,1:nd1)))
			end if
	end select
end subroutine


!convert bodfy force to nodal force
subroutine cbf2nf()
	use solverds
	implicit none
	integer::i,j,k,el1,nnum1
	real(8)::bf1=0,t1=0,r1=0
	type(bc_tydef),allocatable::bc_load1(:),bc_load2(:)
	integer::bfnum1
	
	!allocate(bfload(ndof))
	!bfload=0.0
	!cal the number of load
	bfnum1=0
	do i=1,bfnum
		bfnum1=bfnum1+element(bf(i).node).nnum
	end do
	allocate(bc_load1(bfnum1))

	!t1=sum(element.property(2))
	
	bfnum1=0
	t1=0	
	do i=1, bfnum
		el1=bf(i).node
		nnum1=element(el1).nnum
		do j=1,nnum1
			bf1=0.0
			select case(element(el1).et)
			case(shell3) !shell3 
				bf1=element(el1).property(2)*bf(i).value/3.0 
			case default				
				do k=1,element(el1).ngp
					
					if(element(el1).ec==CAX.or.element(el1).ec==CAX_SPG) THEN
						r1=abs(element(EL1).xygp(1,K))
						if(abs(r1)<1e-7) r1=1.0e-7
					ELSE
						R1=1.0
					END IF
					bf1=bf1+R1*ecp(element(el1).et).Lshape(j,k)* &
						ecp(element(el1).et).weight(k)* &
						element(el1).detjac(k)*bf(i).value 
				end do				
			end select
			bfnum1=bfnum1+1
			bc_load1(bfnum1).value=bf1
			bc_load1(bfnum1).node=element(el1).node(j)
			bc_load1(bfnum1).dof=bf(i).dof
			bc_load1(bfnum1).sf=bf(i).sf
!			t1=t1+bf1
			!write(10,*) bf1
		end do
	end do
	allocate(bc_load2(bfnum1+bl_num))
	if(allocated(bc_load)) then
		bc_load2(1:bl_num)=bc_load(1:bl_num)
		deallocate(bc_load)
	end if
	bc_load2(bl_num+1:bl_num+bfnum1)=bc_load1(1:bfnum1)
	deallocate(bc_load1)

	bl_num=bl_num+bfnum1
	allocate(bc_load(bl_num))
	bc_load=bc_load2
	write(99,*) bc_load.value 
	deallocate(bc_load2)
	
	
	!t1=sum(bc_load.value)
	!print *, t1


end subroutine




!calculate the globe co-ordinates of gausian point for each element.
subroutine xygp(ienum)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::j,k,k1
	

	do j=1,element(ienum).ngp
		
		do k=1,ndimension
			element(ienum).xygp(k,j)=dot_product(ecp(element(ienum).et).lshape(:,j),node(element(ienum).node).coord(k))
			IF(ELEMENT(IENUM).ET==ZT4_SPG2.OR.ELEMENT(IENUM).ET==ZT6_SPG2) element(ienum).xygp(k,j)=dot_product(ecp(element(ienum).et).lshape(:,j),Gnode(K,element(ienum).node2))
		end do
	end do

end subroutine

!calculate the cut-off values of epsilon1 and epsilon2
subroutine cal_epsilon(ienum)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::nlow1,nup2,nlowgp1,nupgp2
	real(8)::hg1,lg1,coord_low1(3),coord_up2(3)

	nlow1=minloc(node(element(ienum).node).coord(ndimension),1)
	nup2=maxloc(node(element(ienum).node).coord(ndimension),1)
	coord_low1=node(element(ienum).node(nlow1)).coord
	coord_up2=node(element(ienum).node(nup2)).coord
	IF(ELEMENT(IENUM).ET==ZT4_SPG2.OR.ELEMENT(IENUM).ET==ZT6_SPG2) THEN
		nlow1=minloc(Gnode(NDIMENSION,element(ienum).node2),1)
		nup2=maxloc(Gnode(NDIMENSION,element(ienum).node2),1)
		coord_low1=Gnode(:,element(ienum).node2(nlow1))
		coord_up2=Gnode(:,element(ienum).node2(nup2))
	ENDIF
	
	nlowgp1=minloc(element(ienum).xygp(ndimension,:),1)
	nupgp2=maxloc(element(ienum).xygp(ndimension,:),1)

	hg1=element(ienum).xygp(ndimension,nlowgp1)-coord_low1(ndimension)
	lg1=abs(element(ienum).xygp(1,nlowgp1)-coord_low1(1))
	if(element(ienum).ec==spg) then
		lg1=lg1**2
		lg1=lg1+(element(ienum).xygp(2,nlowgp1)-coord_low1(2))**2
		lg1=lg1**0.5
	end if
	element(ienum).property(3)=hg1/2.0*((1+(lg1/hg1)**2)**0.5+1)
	if(element(ienum).property(3)<eps1) eps1=element(ienum).property(3)


	hg1=coord_up2(ndimension)-element(ienum).xygp(ndimension,nupgp2)
	lg1=abs(element(ienum).xygp(1,nupgp2)-coord_up2(1))
	if(element(ienum).ec==spg) then
		lg1=lg1**2
		lg1=lg1+(element(ienum).xygp(2,nupgp2)-coord_up2(2))**2
		lg1=lg1**0.5
	end if
	element(ienum).property(2)=hg1/2.0*((1+(lg1/hg1)**2)**0.5+1)
	if(element(ienum).property(2)<eps2) eps2=element(ienum).property(2)
	

	!hg1=node(element(ienum).node(nup2)).coord(ndimension)- &
	!		node(element(ienum).node(nlow1)).coord(ndimension)
	!lg1=node(element(ienum).node(nup2)).coord(1)- &
	!		node(element(ienum).node(nlow1)).coord(1)
	!lg1=abs(lg1)		
	!element(ienum).property(1)=-hg1/2.0*((1+(lg1/hg1)**2)**0.5+1)
	!if(element(ienum).property(1)<minNPH) minNPH=element(ienum).property(1)
	
	
	!element(ienum).property(2:3)=0.1

	!!epslong one(low)
	!element(ienum).property(3)=minval(element(ienum).xygp(ndimension,:))- &
	!		minval(node(element(ienum).node).coord(ndimension))
	!!epslong two(up)
	!element(ienum).property(2)=maxval(node(element(ienum).node).coord(ndimension))- &
	!						maxval(element(ienum).xygp(ndimension,:))

end subroutine 

subroutine BARFAMILY_EXTREMEVALUE(IENUM)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::i,j
	
	do i=1,element(ienum).nnum
		do j=1,ndimension
			if(node(element(ienum).node(i)).coord(j)<barfamily_minxyz(j)) barfamily_minxyz(j)=node(element(ienum).node(i)).coord(j)
			if(node(element(ienum).node(i)).coord(j)>barfamily_maxxyz(j)) barfamily_maxxyz(j)=node(element(ienum).node(i)).coord(j)
		end do
	end do
	
end subroutine

!COUNT AND INDEX THE HEAD DOFS.
SUBROUTINE NDOF_HEAD()
	USE SOLVERDS
	IMPLICIT NONE
    !INTEGER,INTENT(IN)::ISTEP
	INTEGER::I
	
	NDOFHEAD=0
	DO I=1,NNUM
		IF(NODE(I).DOF(4)>0) NDOFHEAD=NDOFHEAD+1
	ENDDO
	
	IF(NDOFHEAD<1) RETURN
	
	IF(.NOT.ALLOCATED(DOFHEAD)) ALLOCATE(DOFHEAD(NDOFHEAD))
	NDOFHEAD=0
	DO I=1,NNUM
		IF(NODE(I).DOF(4)>0) THEN
			NDOFHEAD=NDOFHEAD+1
			DOFHEAD(NDOFHEAD)=NODE(I).DOF(4)
		ENDIF
	ENDDO
	

ENDSUBROUTINE


!SUBROUTINE ZT_SPG_INI(IELT)
!    USE solverds
!    USE SolverMath
!    IMPLICIT NONE
!    INTEGER,INTENT(IN)::IELT
!    integer::i,j,K
!	integer::n1,n2,n3,n4,AELT1(2),IN1(10)
!	real(kind=DPN)::t1=0,vcos=0,vsin=0,rpi,coord1(3,4)=0,cent1(3,2)=0,c1(3,3)=0,VEC1(3),t2
!    
!    
!    if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
!    !找出zt单元的由face1指向face2的向量VEC1
!    IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
!        AELT1=ELEMENT(IELT).ADJELT([1,3])
!    ELSE
!        AELT1=ELEMENT(IELT).ADJELT([1,2])
!    ENDIF
!    where(AELT1<1)   AELT1=IELT
!
!    DO I=1,2
!        DO J=1,NDIMENSION
!            CENT1(J,I)=SUM(NODE(ELEMENT(AELT1(I)).NODE).COORD(J))/ELEMENT(AELT1(I)).NNUM            
!        ENDDO
!    ENDDO
!    VEC1=(CENT1(:,2)-CENT1(:,1))
!    t1=NORM2(VEC1)
!    
!    IF(ABS(t1)<1.E-6) THEN
!        PRINT *, 'ZT单元I的指向向量无法确定.I=',ielt
!    ELSE
!        vec1=vec1/t1
!    ENDIF
!    
!    T2=1.0D0
!	IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
!        
!		!GENERATE THE GHOST NODE
!		T1=((NODE(ELEMENT(IELT).NODE(1)).COORD(1)-NODE(ELEMENT(IELT).NODE(2)).COORD(1))**2+ &
!			(NODE(ELEMENT(IELT).NODE(1)).COORD(2)-NODE(ELEMENT(IELT).NODE(2)).COORD(2))**2)**0.5D0
!		!LOCAL TO GLOBAL 
!		VCOS=(NODE(ELEMENT(IELT).NODE(2)).COORD(1)-NODE(ELEMENT(IELT).NODE(1)).COORD(1))/T1
!		VSIN=(NODE(ELEMENT(IELT).NODE(2)).COORD(2)-NODE(ELEMENT(IELT).NODE(1)).COORD(2))/T1
!		C1(1,1)=VCOS
!		C1(1,2)=-VSIN
!		C1(2,2)=VCOS
!		C1(2,1)=VSIN
!        t2=sign(1.d0,DOT_PRODUCT(VEC1(1:2),C1(1:2,2)))
!        C1(1:2,2)=T2*C1(1:2,2)        
!        !IF(DOT_PRODUCT(VEC1(1:2),C1(1:2,2))<0.D0) C1(1:2,2)=-C1(1:2,2)
!        !4#节点的局部坐标
!		!COORD1=0.D0;COORD1(2,1)=Matproperty(ELEMENT(IELT).MAT,14,1)
!  !      IF(ABS(COORD1(2,1))<1E-6) COORD1(2,1)=1.0D0
!		!COORD1(1,2)=0.D0;COORD1(2,2)=COORD1(2,1)
!		!IF(.NOT.ALLOCATED(GNODE)) ALLOCATE(GNODE(3,1000))
!		!IF(NGNODE+4>SIZE(GNODE,DIM=2)) CALL ENLARGE_GNODE(1000)
!		!NGNODE=NGNODE+1;GNODE(:,NGNODE)=NODE(ELEMENT(IELT).NODE(1)).COORD
!		!NGNODE=NGNODE+1;GNODE(:,NGNODE)=NODE(ELEMENT(IELT).NODE(2)).COORD
!		!COORD1(1:2,1)=MATMUL(C1(1:2,1:2),COORD1(1:2,1))
!		!NGNODE=NGNODE+1;GNODE(:,NGNODE)=COORD1(:,1)+NODE(ELEMENT(IELT).NODE(2)).COORD
!		!NGNODE=NGNODE+1;GNODE(:,NGNODE)=COORD1(:,1)+NODE(ELEMENT(IELT).NODE(1)).COORD
!		!ALLOCATE(ELEMENT(IELT).NODE2(4))
!		!ELEMENT(IELT).NODE2=[NGNODE-4+1:NGNODE:1]
!        N2=2
!        !IF(T2<0) THEN
!        !    IN1()=[3,4,0]
!        !ELSE
!        !    IN1=[1,2,0]            
!        !ENDIF
!                
!	ENDIF
!	IF(ELEMENT(IELT).ET==ZT6_SPG) THEN
!	    !LOCAL TO GLOBAL
!                    
!        do j=1,3
!            C1(:,j)=NODE(ELEMENT(IELT).NODE(j)).COORD;
!        enddo
!        C1(3,:)=NORMAL_TRIFACE(C1)
!        C1(3,:)=C1(3,:)/norm2(C1(3,:))
!        
!        t2=sign(1.d0,dot_product(vec1,c1(3,:)))
!        C1(3,:)=t2*C1(3,:)
!
!        T1=NORM2(NODE(ELEMENT(IELT).NODE(1)).COORD-NODE(ELEMENT(IELT).NODE(2)).COORD)
!        C1(1,:)=(NODE(ELEMENT(IELT).NODE(2)).COORD-NODE(ELEMENT(IELT).NODE(1)).COORD)/T1
!        C1(2,:)=0.d0
!        C1(2,:)=NORMAL_TRIFACE(RESHAPE([C1(2,:),C1(1,:),C1(3,:)],([3,3])))
!        C1(2,:)=C1(2,:)/norm2(C1(2,:))
!        C1=TRANSPOSE(C1)
!        !4#节点的坐标
!        IF(T2<0) THEN
!            IN1(1:6)=ELEMENT(IELT).NODE([4,6,5,1,3,2])
!            ELEMENT(IELT).NODE=IN1(1:6)
!            IN1(1:9)=-ELEMENT(IELT).EDGE([6,5,4,3,2,1,7,9,8])
!            ELEMENT(IELT).EDGE=IN1(1:9)
!            IN1(1:5)=ELEMENT(IELT).FACE([2,1,5,4,3])
!            ELEMENT(IELT).FACE=IN1(1:5)
!            IN1(1:5)=ELEMENT(IELT).ADJELT([2,1,5,4,3])
!            ELEMENT(IELT).ADJELT=IN1(1:5)
!        ENDIF
!        N2=3
!        
!    ENDIF 
!    
!    
!    
!	IF(.NOT.ALLOCATED(GNODE)) ALLOCATE(GNODE(3,1000))
!	IF(NGNODE+2*N2>SIZE(GNODE,DIM=2)) CALL ENLARGE_GNODE(1000)
!    COORD1=0.D0;COORD1(N2,1)=Matproperty(ELEMENT(IELT).MAT,14,1)
!    IF(ABS(COORD1(N2,1))<1E-6) COORD1(N2,1)=1.0D0
!    
!	COORD1(1:N2,1)=MATMUL(C1(1:N2,1:N2),COORD1(1:N2,1))
!   
!    DO I=1,N2    
!	    !NGNODE=NGNODE+1;
!        GNODE(:,NGNODE+I)=NODE(ELEMENT(IELT).NODE(I)).COORD
!        IF(ELEMENT(IELT).ET==ZT6_SPG) THEN
!            N3=I
!        ELSE
!            N3=MOD(I,N2)+1
!        ENDIF
!        GNODE(:,NGNODE+N2+N3)=COORD1(:,1)+NODE(ELEMENT(IELT).NODE(I)).COORD
!    ENDDO
!
!	ALLOCATE(ELEMENT(IELT).NODE2(2*N2))
!	ELEMENT(IELT).NODE2=[NGNODE+1:NGNODE+2*N2]
!    
!    NGNODE=NGNODE+2*N2
!
!ENDSUBROUTINE

