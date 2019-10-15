!Assumption
!1) the most bottom layer never intersects with upper layer.
!
!Generate 3D 6-node prism element.
Subroutine Generate_3D_MODEL()
	use meshds
	implicit none
	logical::istet
	integer::et1,i,j,k,n1,n2
    real(8)::t1

	
	!Output the node and read in the elevation data
    call st_membrance()
    
	call First_Order_PrismElement_Gen()
    
    call meshoutput()
end Subroutine
    
!Subroutine Generate_tet_element()
!	use meshds
!	implicit none
!	logical::istet
!	integer::et1
!	
!	!Output the node and read in the elevation data
!	call st_membrance()
!	
!	istet=.true.
!    
!	call First_Order_PrismElement_Gen(istet)
!    !call out_elt_gmsh()
!    
!end Subroutine

    

	subroutine First_Order_PrismElement_Gen()
	    use meshds
	    implicit none
		!logical,intent(in)::istet
	    integer::i,j,k,n1,n2,n3,n4(2)=0,iat(3)=0,Vtet(4)=0,shift1=1
	    real(8)::t1,t2
	    type(element_tydef),allocatable::ELT1(:),PRMELT(:),TetELT(:)
	    integer::NPRMELT=0,NTetelt=0
	    
	    !if(et1==153) call gen_6_noded_triangle
        
        
	    
	    allocate(PRMELT(soillayer*nelt))
        
	    allocate(Tetelt(3*soillayer*nelt))
        
	    elt.ismodel=0
	    do i=1,nelt
            
	      if(elt(i).isdel) cycle         
	         do j=1,soillayer
				 if(zone(elt(i).zn).ismodel(j)/=1) cycle
                 if(index(zone(elt(i).zn).solver_et(j),'cpe')>0) then
                    elt(i).ismodel=1
                    cycle
                 endif
                 if(index(zone(elt(i).zn).solver_et(j),'prm')<1.and. &
                 index(zone(elt(i).zn).solver_et(j),'tet')<1) cycle
                 nprmelt=nprmelt+1
                 prmelt(nprmelt).mother=i
                 !prmelt(nprmelt).matlayer=j
                 prmelt(nprmelt).et=63
                 prmelt(nprmelt).ilayer=j
                 !prmelt(nprmelt).nlayer(2)=j
                 prmelt(nprmelt).nnum=6
                 prmelt(nprmelt).zn=elt(i).zn
                 prmelt(nprmelt).mat=zone(elt(i).zn).mat(j)
                 prmelt(nprmelt).ismodel=1
                 IF(.NOT.ALLOCATED(zone(elt(i).zn).nprm)) then
                    allocate(zone(elt(i).zn).nprm(soillayer))
                    zone(elt(i).zn).nprm=0
                 endif
                 zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)+1
                 
                 iat=0
                 
                 do k=1,3
						prmelt(nprmelt).node(k)=elt(i).node(k)+(j-1)*nnode
						!�ؽڵ㣨�߳���ȵĽڵ㣩
						if(node(prmelt(nprmelt).node(k)).subbw/=0) &
							prmelt(nprmelt).node(k)=node(prmelt(nprmelt).node(k)).subbw
							
						prmelt(nprmelt).node(3+k)=elt(i).node(k)+(j)*nnode
						!�ؽڵ㣨�߳���ȵĽڵ㣩
						if(node(prmelt(nprmelt).node(3+k)).subbw/=0) & 
							prmelt(nprmelt).node(3+k)=node(prmelt(nprmelt).node(3+k)).subbw
						
						if(prmelt(nprmelt).node(k)==prmelt(nprmelt).node(3+k)) then
							iat(k)=1
						end if						
                 end do
				 
                 !if(any(prmelt(nprmelt).node(1:6)==0)) then
                 !   pause
                 !endif
                 
                 IF(COUNT(IAT==1)==3) THEN
                    !deallocate(prmelt(nprmelt).node)
					nprmelt=nprmelt-1
                    zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)-1
                    CYCLE
                 ENDIF
                 
                 if(index(zone(elt(i).zn).solver_et(j),'tet')>0) then	
                    IF(.NOT.ALLOCATED(zone(elt(i).zn).ntet) ) then
                        allocate(zone(elt(i).zn).ntet(soillayer))
                        zone(elt(i).zn).ntet=0
                     endif                                      
					 select case(count(iat==1)) 
							case(1) !������������嵥Ԫ
								prmelt(nprmelt).nnum=4
								prmelt(nprmelt).et=43                                
                                n1=maxloc(iat,dim=1)
								Vtet(1)=mod(n1,3)+1
								Vtet(2)=mod(Vtet(1),3)+1
                                Vtet(3)=Vtet(2)+3
                                Vtet(4)=Vtet(1)+3
                                
                                do k=1,2
                                    ntetelt=ntetelt+1  
									tetelt(ntetelt)=prmelt(nprmelt)	
                                    if(k==1) tetelt(ntetelt).node(1:3)=prmelt(nprmelt).node(vtet([1,3,2]))
	                                if(k==2) tetelt(ntetelt).node(1:3)=prmelt(nprmelt).node(vtet([1,4,3]))
									tetelt(ntetelt).node(4)=prmelt(nprmelt).node(n1)
                                    zone(elt(i).zn).ntet(j)=zone(elt(i).zn).ntet(j)+1
                                end do    
                                
								!deallocate(prmelt(nprmelt).node)
								nprmelt=nprmelt-1
								zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)-1
							case(2) !���һ�������嵥Ԫ
								
								prmelt(nprmelt).nnum=4
								prmelt(nprmelt).et=43
								n1=prmelt(nprmelt).node(3+minloc(iat,dim=1))
								prmelt(nprmelt).node(1:3)=prmelt(nprmelt).node(1:3)
								prmelt(nprmelt).node(4)=n1
								prmelt(nprmelt).node(5:6)=0
								ntetelt=ntetelt+1
								tetelt(ntetelt)=prmelt(nprmelt)
								zone(elt(i).zn).ntet(j)=zone(elt(i).zn).ntet(j)+1
								
								!deallocate(prmelt(nprmelt).node)
								nprmelt=nprmelt-1
                                zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)-1
							case(3) !dead,delect
							
								!deallocate(prmelt(nprmelt).node)
								nprmelt=nprmelt-1
                                zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)-1
							case default    !�ֽ�Ϊ3������
                                prmelt(nprmelt).nnum=4
								prmelt(nprmelt).et=43
                                
                                do k=1,3
                                    ntetelt=ntetelt+1
								    tetelt(ntetelt)=prmelt(nprmelt)
                                    if(k==1) tetelt(ntetelt).node(1:4)=prmelt(nprmelt).node(1:4)
                                    if(k==2) tetelt(ntetelt).node(1:4)=prmelt(nprmelt).node([2,5,6,4])
                                    if(k==3) tetelt(ntetelt).node(1:4)=prmelt(nprmelt).node([2,6,3,4])
                                enddo
                                zone(elt(i).zn).ntet(j)=zone(elt(i).zn).ntet(j)+3
                                !deallocate(prmelt(nprmelt).node)
								nprmelt=nprmelt-1
                                zone(elt(i).zn).nprm(j)=zone(elt(i).zn).nprm(j)-1                                
					 end select 
                 end if                    
	         end do
	    end do
	    
        
        allocate(elt1(nelt+nprmelt+ntetelt))
        elt1(1:nelt)=elt
        if(nprmelt>0) elt1(nelt+1:nelt+nprmelt)=prmelt
        if(ntetelt>0) elt1(nelt+nprmelt+1:nelt+nprmelt+ntetelt)=tetelt
        deallocate(elt)
        allocate(elt,source=elt1)
        deallocate(elt1,prmelt,tetelt)
        nelt=nelt+nprmelt+ntetelt
        nprmelt=0;ntetelt=0
	    
	end subroutine
	

    

	subroutine st_membrance()
		use meshds
		implicit none
		integer::i,j,k,ng,maxbw1,a1,iflag,n1,j1,j2,n2,k1,pv1(3),err
		logical::isinp=.true.
		real(8)::t1,b1(3),c1(3),xy(2,3),t2
		real(8),allocatable::Tm1(:),tm1b(:),load1(:),at1(:)	
		INTEGER,ALLOCATABLE::IPERM(:),bw1(:)
	    logical,allocatable::Lt1(:)
        type(point_tydef),allocatable::node1(:)
 
		
	 
		!if(nmeminp2>0) call seg_initialize()
        
        if(.not.ismeminpdone.or.inpmethod==0.or.any(node(1:).havesoildata<2)) then
	
		    node.subbw=-999
		
		    do i=1,nelt
			    if(elt(i).isdel) cycle
			    n1=3
			    if(elt(i).et==6) n1=6
			    if(elt(i).et==15) n1=15
			    node(elt(i).node(1:n1)).subbw=0
		    end do

		    ng=0

		    do i=1,nnode
			    if(node(i).subbw/=0) then
				    cycle
			    end if
			    ng=ng+1
			    node(i).number=ng  !!!!
		    end do

		    ALLOCATE(IPERM(nnode),Noutputorder(nnode),bw1(nnode),stat=err)
		    call reorder_nodal_number(IPERM,nnode,node(1:),nelt,elt)
            IPERM=IPERM-(NNODE-NG) !ȥ������ģ���еĽڵ�
		    do i=1,nnode
			    if(node(i).subbw==0) then
				    noutputorder(IPERM(node(i).number))=i	
				    node(i).number=IPERM(node(i).number)
			    end if			
		    end do		
	    
		    call csb(nnode,node(1:),nelt,elt)
		
		    bw1=0
		    bw1(1)=node(noutputorder(1)).bw
		    maxbw1=bw1(1)
		    do i=2,ng
			    n1=noutputorder(i)
			    if(maxbw1<node(n1).bw) maxbw1=node(n1).bw
			    bw1(i)=node(n1).bw+bw1(i-1)								
		    end do
		
		    !assemble the total km
		    allocate(tm1(bw1(ng)),tm1b(bw1(ng)),load1(ng))
		    tm1=0.0
		
		    do i=1,nelt
			    if(elt(i).isdel) cycle
			    pv1(1:3)=node(elt(i).node(1:3)).number
			    select case(elt(i).et)
				    case(0)

					    xy(1,1:3)=node(elt(i).node(1:3)).x
					    xy(2,1:3)=node(elt(i).node(1:3)).y
				
					    do j=1,3
						    j1=mod(j,3)+1
						    j2=mod(j1,3)+1
						    b1(j)=xy(2,j1)-xy(2,j2)
						    c1(j)=xy(1,j2)-xy(1,j1)
					    end do
		      
					    t1=(b1(1)*c1(2)-b1(2)*c1(1))/2
					    do j=1,3
					      do k=1,3				    
							    if(pv1(j)<pv1(k))  cycle   
							    t2=(b1(j)*b1(k)+c1(j)*c1(k))/(4*t1)
							    a1=bw1(pv1(j))-(pv1(j)-pv1(k))
							    tm1(a1)=tm1(a1)+t2							  
						    end do
					    end do
					
					
				    case default
					    Print *, 'Error! Membrance Interpolation. Only 3-noded element is considered.'
					    stop
				
			    end select
	
			
		    end do
		
		    !�γɱ߽�
		
		    call bc_meminp()	
		
		


		    !���
		
		    do i=0,soillayer
			    tm1b=tm1
			    load1=0.D0
			    !��״�߽�
			    do j=1,nmeminp
				    do k=1,meminp(j).nvb
					    if(abs(meminp(j).vbc(k,i)+999)<1e-6) cycle
                        n1=node(meminp(j).nbc(k)).number
					    tm1b(bw1(n1))=um
					    load1(n1)=meminp(j).vbc(k,i)*um
				    end do
			    end do
			    do j=1,nmeminp2
				    do k=1,meminp2(j).nvb
                    
					    if(abs(meminp2(j).vbc(k,i)+999)<1e-6) cycle
                        n1=node(meminp2(j).nbc(k)).number
					    tm1b(bw1(n1))=um
					    load1(n1)=meminp2(j).vbc(k,i)*um
				    end do
			    end do
			    !��״�߽�
			    do j=1,ngeo
				    if(geology(j).isini==0) cycle
				    if(abs(geology(j).elevation(i)+999)<1e-6) cycle
				    tm1b(bw1(node(geology(j).node).number))=um
				    load1(node(geology(j).node).number)=geology(j).elevation(i)*um	  
			    end do
			
			    call chodec(tm1b,bw1,ng)
			    call chosol(tm1b,bw1,load1,ng,maxbw1)
			
			
			    !if(.not.allocated(elevation)) allocate(elevation(0:soillayer,nnode))
			    !n1=nnode*i
			    do j=1,ng
                    n2=noutputorder(j)
				    !node(n1+n2).z=load1(j)
                    if(.not.allocated(node(n2).elevation)) then
                        allocate(node(n2).elevation(0:soillayer))                    
                    endif
                    node(n2).elevation(i)=load1(j)
				    !elevation(i,noutputorder(j))=load1(j)  !!!
			    end do
		    end do
		
		    !		
		    DEALLOCATE(noutputorder,bw1,tm1,tm1b,load1,stat=err)
		
            deallocate(IPERM,stat=err)

	    endif
        
        ismeminpdone=.true.
        node(1:nnode).havesoildata=2        
            	!enlarge node space
        allocate(node1(nnode))
	    node1=node(1:nnode)
	    deallocate(node)
	    allocate(node(nnode*(soillayer+1)))
	    do i=0,soillayer
		    n1=nnode*i
            if(i==0) then
		        node((n1+1):(n1+nnode))=node1(1:nnode)
            else
                node((n1+1):(n1+nnode)).x=node1(1:nnode).x;node((n1+1):(n1+nnode)).y=node1(1:nnode).y;
            endif
            node((n1+1):(n1+nnode)).z=node1(1:nnode).elevation(i);
		    node((n1+1):(n1+nnode)).layer=i			
	    end do
	    deallocate(node1)
    
    	    !check elevation
        
	    do j=0,soillayer
		    n1=nnode*j
		    do i=1,nnode
			    !recover
			    node(n1+i).number=n1+i
			    node(n1+i).subbw=0
		    end do
	    end do
        
	    allocate(at1(0:soillayer),Lt1(0:soillayer))	
    
	    do i=1,nnode
		    Lt1=.false.
		    do j=0, soillayer-1
			    if(Lt1(j)) cycle
			    n1=(j)*nnode+i
			    do k=j+1,soillayer
				    if(Lt1(k)) cycle
				    n2=(k)*nnode+i
				    t1=node(n1).z-node(n2).z
				    if(abs(t1)<1e-4) then
					    node(n2).subbw=n1 !<>0,dead
					    Lt1(k)=.true.					
                    end if
                end do
		    end do
	    end do	
		

        

        
	end subroutine
	
	
	subroutine bc_meminp()
		use ds_t
		use meshds
		implicit none
		integer::i,j,k,k1,n1,n2,j1
        integer,allocatable::node1(:)
		real(8)::t1,t2,ar1(1000),x0,y0
		
		
		allocate(node1(nnode))
        
		do i=1,nmeminp2
		
			!ͳ�Ƹõ��������ܽڵ�����
			node1=0
			do j=1,meminp2(i).nnum-1
				!n1=n1+seg(segindex(meminp2(i).cp(j),meminp2(i).cp(j+1))).nnum-1 !���˵��ݲ����룬�������
                node1(seg(segindex(meminp2(i).cp(j),meminp2(i).cp(j+1))).node)=j
			end do
			!������պϣ������һ���˽ڵ���롣
			!if(meminp2(i).cp(1)/=meminp2(i).cp(meminp2(i).nnum)) n1=n1+1 
			
			n1=count(node1>0)
            meminp2(i).nvb=n1
			allocate(meminp2(i).nbc(n1),meminp2(i).vbc(n1,0:soillayer),meminp2(i).niseg(n1))
			meminp2(i).nbc=pack([1:nnode],node1>0)
            meminp2(i).niseg=pack(node1,node1>0)
            
		
		
			!�߲�ϵ��
			t1=0
			t2=0
			do j=1,meminp2(i).nnum-1
				k=j+1
				n1=meminp2(i).cp(j)
				n2=meminp2(i).cp(k)
				t2=((arr_t(n1).x-arr_t(n2).x)**2+ &
						 (arr_t(n1).y-arr_t(n2).y)**2)**0.5
				if(abs(t2)<1e-7) then
					print *, 'Error in st_membrance'
					pause
				end if
				
				if(j==1) allocate(meminp2(i).lincof(2,meminp2(i).nnum-1,0:soillayer))
				do k1=0,soillayer
					meminp2(i).lincof(1,j,k1)=(meminp2(i).elevation(k,k1)-meminp2(i).elevation(j,k1))/t2
					meminp2(i).lincof(2,j,k1)=meminp2(i).elevation(j,k1)
				end do
				t1=t2
			end do		
		
			
			
			do j=1,meminp2(i).nvb
                n2=meminp2(i).niseg(j)
                n1=meminp2(i).cp(n2)
                j1=meminp2(i).nbc(j)	
				x0=arr_t(n1).x;y0=arr_t(n1).y
				t1=((node(j1).x-x0)**2+ &
						 (node(j1).y-y0)**2)**0.5                
				do k=0,soillayer
					meminp2(i).vbc(j,k)=meminp2(i).lincof(1,n2,k)*t1+meminp2(i).lincof(2,n2,k)
				end do				

			end do			
			

		end do		
		
        deallocate(node1)
		
	end subroutine
	
	
	!subroutine seg_initialize()
	!	use ds_t
	!	use meshds
	!	implicit none
	!	integer::i,nnode1(1000)
	!	real(8)::t1,t2
	!	logical::tof1
	!	
	!	do i=1,nseg
	!		cpp=>cpphead(seg(i).icl)
	!		tof1=.false.			
	!		do while(.true.)
	!			t1=((arr_t(seg(i).sv).x-cpp.npt.x)**2+(arr_t(seg(i).sv).y-cpp.npt.y)**2)**0.5
	!			if(abs(t1)<precision) then
	!				seg(i).svp=>cpp
	!				if(.not.associated(seg(i).evp)) tof1=.true.
	!			end if
	!			t2=((arr_t(seg(i).ev).x-cpp.npt.x)**2+(arr_t(seg(i).ev).y-cpp.npt.y)**2)**0.5
	!			if(abs(t2)<precision) then
	!				seg(i).evp=>cpp
	!				if(.not.associated(seg(i).svp)) then
	!					seg(i).isa2z=0
	!					tof1=.true.
	!				end if
	!			end if
	!			if(tof1) then
	!				seg(i).nnum=seg(i).nnum+1
	!				nnode1(seg(i).nnum)=cpp.npt.number
	!			end if
	!			if(associated(seg(i).evp).and.associated(seg(i).svp)) exit
	!			cpp=>cpp.next
	!			if(.not.associated(cpp).or.associated(cpp,cpphead(seg(i).icl)).or.associated(cpp,bnhead)) exit 
	!		end do
	!		!�����β�������Σ������ת,����ͳ�ƽڵ�����
	!		tof1=.false.		
	!		if(seg(i).ist2h==1) then
	!			if(seg(i).isa2z==0) then
	!				seg(i).isa2z=1
	!				cpp=>seg(i).svp
	!				cpptail=>seg(i).evp
	!			else
	!				tof1=.true. !mark
	!				seg(i).isa2z=0
	!				cpp=>seg(i).evp
	!				cpptail=>seg(i).svp
	!			end if
 !
	!			seg(i).nnum=0
	!			do while(.true.)
	!				seg(i).nnum=seg(i).nnum+1
	!				nnode1(seg(i).nnum)=cpp.npt.number
	!				cpp=>cpp.next
	!				if(cpp.npt.number==cpptail.npt.number) then !!!!!why using if(associated(cpp,cpptail)) is wrong? 
	!					seg(i).nnum=seg(i).nnum+1
	!					nnode1(seg(i).nnum)=cpp.npt.number
	!					exit
	!				end if
	!			end do
	!			if(tof1) nnode1(1:seg(i).nnum)=nnode1(seg(i).nnum:1:-1)
	!		end if
	!		
	!		allocate(seg(i).node(seg(i).nnum))
	!		if(seg(i).isa2z==1) then
	!			seg(i).node=nnode1(1:seg(i).nnum)
	!		else
	!			seg(i).node=nnode1(seg(i).nnum:1:-1)
	!		end if
 !
 !           
	!	end do
	!
	!end subroutine
	
	subroutine seg_initialize()
		use ds_t
		use meshds
		implicit none
		integer::i,j,nnode1(1000),n1
		real(8)::t1,t2,xi,yi,xj,yj,x1,y1
		logical::tof1
		
        node.subbw=0!
        
		do i=1,nseg
            
            xi=arr_t(seg(i).sv).x;yi=arr_t(seg(i).sv).y
            xj=arr_t(seg(i).ev).x;yj=arr_t(seg(i).ev).y
            n1=0;nnode1=0            
            do j=1,ncedge
                if(cedge(j).cl/=seg(i).icl) cycle
                x1=sum(node(cedge(j).v).x)/2.0
                y1=sum(node(cedge(j).v).y)/2.0
                call vins2d(xi,yi,xj,yj,x1,y1,tof1)
                if(tof1) then
                    n1=n1+1
                    nnode1(n1)=cedge(j).edge
                    node(cedge(j).v).subbw=i
                endif
			enddo
            
            if(n1>0) then
			    seg(i).nnum=count(node(1:).subbw==i)
			    allocate(seg(i).node(seg(i).nnum))
                seg(i).node=pack([1:nnode],node(1:).subbw==i)
                seg(i).nedge=n1
			    allocate(seg(i).edge(n1))
                seg(i).edge=nnode1(1:n1)
            endif
            
		end do
        
        node.subbw=0
	
	end subroutine
    
    
    
	
subroutine chodec(tm_t1,sbw,nnum_t) !cholesky decomposition
	implicit none
	integer::tmn,nnum_t,sbw(nnum_t)
	real(8)::tm_t1(sbw(nnum_t))
	integer::i,j
	integer::K1,i1,j1
	real(8)::sum,art(25)
	real(8),allocatable::T(:)
	allocate(T(nnum_t))
	t=0

	tmn=sbw(nnum_t)

	do i=2,nnum_t
		i1=i-sbw(i)+sbw(i-1)+1  !�̨�iDD�̨���???��?��??a??��?��Do?
		do j=i1,i-1 
			if(j==1) then
				j1=1
			else
				j1=j-sbw(j)+sbw(j-1)+1
			end if
			if(i1>j1) j1=i1
				sum=0
			do K1=j1,j-1
				sum=sum+T(K1)*tm_t1(sbw(j)-j+K1)
			end do
				T(j)=tm_t1(sbw(i)-i+j)-sum
				if(abs(tm_t1(sbw(j)))<1e-15) tm_t1(sbw(j))=1e-15
				tm_t1(sbw(i)-i+j)=T(j)/tm_t1(sbw(j))
				!if(isnan(tm_t1(sbw(i)-i+j))) pause
				tm_t1(sbw(i))=tm_t1(sbw(i))-T(j)*tm_t1(sbw(i)-i+j)
		end do     
	end do 
	deallocate(T)
end subroutine

subroutine chosol(tm_t1,sbw,value,nnum_t,maxbdw_t)

	implicit none
	integer::tmn,nnum_t,maxbdw_t,sbw(nnum_t)
	real(8)::tm_t1(sbw(nnum_t)),value(nnum_t)
	integer::i,j
	integer::i1,j1

	tmn=sbw(nnum_t)

	do i=2,nnum_t
		i1=i-sbw(i)+sbw(i-1)+1
		do j=i1,i-1
			value(i)=value(i)-tm_t1(sbw(i)-i+j)*value(j)
		end do
	end do

	do i=1,nnum_t
		value(i)=value(i)/tm_t1(sbw(i))
	end do
	do i=nnum_t-1,1,-1
		i1=i+Maxbdw_t
		if(i1>nnum_t) i1=nnum_t
		do j=i+1,i1
			j1=j-sbw(j)+sbw(j-1)+1
			if(i>=j1) value(i)=value(i)-tm_t1(sbw(j)-j+i)*value(j)
		end do  
	end do

end subroutine
