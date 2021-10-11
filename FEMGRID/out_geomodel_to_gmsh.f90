module geomodel
	USE IFQWIN 
    use meshDS, only : node,nnode,elt,nelt,edge,nedge,adjlist,zone,znum, &
						bnode,nbnode,soillayer,&
						path_name,ENUMBER,xyscale,xmin,ymin,write_help,strtoint,&
                        ENLARGE_AR,adjlist_tydef,addadjlist,removeadjlist
    use ds_t
    use CutoffWall,only:xzone,nxzone
    implicit none
    private
    public out_volume_to_gmsh,BLO,NBLO,CHECK_ORIENT,Gen_SubElement,out_volume_to_gmsh2,NGNODE
    
    INTERFACE INCOUNT
        MODULE PROCEDURE INCOUNT,INCOUNT2
    END INTERFACE

    INTERFACE GM_ENLARGE_AR
        MODULE PROCEDURE FACE_ENLARGE_AR,SELT_ENLARGE_AR,GEDGE_ENLARGE_AR
    END INTERFACE
    
    INTEGER::ngnode=0
    
    INTEGER,ALLOCATABLE::BEDGE(:),ZBEDGE(:,:) !区边界inode对应的每个节点生成的edge的整体编号, ZBEGGE(ISOILLAYER,INODE)=1表示边区边界inode对应的第isoillayer的厚度（竖向线段长度）为0;
    TYPE OFF_MODEL_TYDEF
        INTEGER::NNODE=0,NFACE=0
        REAL(8),ALLOCATABLE::NODE(:,:)
        INTEGER,ALLOCATABLE::FACE(:,:)
    ENDTYPE
    TYPE(OFF_MODEL_TYDEF)::OFF_MATHER_MODEL
    type lineloop_tydef
        integer::nnum=0
        integer,allocatable::node(:)
        real(8),allocatable::elevation(:)
    endtype
    type blo_tydef
        INTEGER::IBLO=1
        integer::nloop=0,ndim=2,itype=2
        type(lineloop_tydef),allocatable::bloop(:)

        
        character(1024)::helpstring= &
            & "\n blo的输入格式为:\n &
            & 1)nblo //结构体数  \n &	
            & 2) nloop,itype //分别为结构体的边界环数及类型。\n &
            &   2.1) nPt   //环的边界点数,如为闭环，则令首尾节点相同 \n&
            &   2.2）P1,P2,...Pnpt. //边界点号\n &
            &   2.3) E1,E2,...Enpt,[E12,E22,...Enpt2]. //边界点高程，[E12...]仅当itype=4时输入.\n  &
            & notes: \n &
            & 1）itype =0 x-y平面上的多边形，高程为E1; \n &
            &          =1,棱柱体，棱平行于z轴，上下底面平行于x-y平面，上下底面的高程为E1和E2.这时nloop只能为1; \n &
            &          =2,3D空间多边形，各点高程为E1,E2,...,Enpt; \n & 
            &          =3,3D空间直线，各点高程为E1,E2,...,Enpt; \n & 
            &          =4,平行z轴面，模拟防渗墙，各点底高程为E1,E2,...,Enpt,顶高程为E12,E22,...Enpt2; \n &     
            & "C
    contains
        procedure,nopass::help=>write_help
        procedure::readin=>blo_readin
        procedure::write=>blo_output
    endtype
    TYPE(blo_tydef),ALLOCATABLE::BLO(:)
    INTEGER::NBLO=0
    
    type subelt_tydef
        logical::isdel=.true.
        integer::et=63,nnum=6,ilayer,zn,nface,isze=0 !et=83 6面体单元;isze，是否为0厚度的土层单元
        integer,allocatable::node(:)
        integer,allocatable::adj(:),face(:) !下、上表面，四周
        integer::cutface=0,subelt(2)=0
        integer::isb=0
    endtype
    type(subelt_tydef),allocatable::selt(:)
    INTEGER::NSELT=0
    
    type face_tydef
        logical::isdel=.true.
        integer::nnum=3,e(2)=-1,isb=0,isv=0,iptr !isb=1,zone boundary face;=2,model boundary face;=0,zone inside face;=3,same zone but differen soillayer 
        !isv=1,vertical face;=0,horizontal face.
        integer::ilayer
        integer,allocatable::node(:)
        integer,allocatable::edge(:)
        integer::cutedge=0,subface(2)=0 
    contains
        procedure::getface=>get_child_face
    endtype
    type(face_tydef),allocatable::face(:)
    integer::nface=0
    
    type gedge_tydef
        logical::isdel=.false.
        integer::v(2)=0
        integer::isoverlap=0,isb=0,nface=0,iptr=0 !!isb=1,zone boundary edge;=2,model boundary edge;=0,zone inside edge
        integer::cutpoint=0,subedge(2)=0
        integer,allocatable::face(:),subid(:)
    contains
        procedure::getedge=>get_child_edge
    endtype
    type(gedge_tydef),allocatable::gedge(:)
    integer::ngedge=0
    
    type(adjlist_tydef),allocatable::gadjlist(:)
    
    contains
    
    recursive subroutine get_child_edge(this,iedge,children,nchd) 
        class(gedge_tydef)::this
        integer,intent(in)::iedge
        integer,intent(inout)::nchd
        integer,allocatable::children(:)
        integer::i
        
        if(.not.allocated(children)) then
            allocate(children(1))
            nchd=0      
        endif
        
        if(this.isdel==.false.) then
            children=[children(:nchd),iedge]
            nchd=nchd+1
        elseif(this.cutpoint>0) then
            do i=1,2
                if(gedge(this.subedge(i)).isdel==.false.) then
                    children=[children(:nchd),this.subedge(i)]
                    nchd=nchd+1
                else
                    call get_child_edge(gedge(this.subedge(i)),this.subedge(i),children,nchd)
                endif
            enddo            
        endif
        
    endsubroutine

    recursive subroutine get_child_face(this,iface,children,nchd) 
        class(face_tydef)::this
        integer,intent(in)::iface
        integer,intent(inout)::nchd
        integer,allocatable::children(:)
        integer::i
        
        if(.not.allocated(children)) then
            allocate(children(1))
            nchd=0      
        endif
        
        if(this.isdel==.false.) then
            children=[children(:nchd),iface]
            nchd=nchd+1
        elseif(this.cutedge>0) then
            do i=1,2
                if(face(this.subface(i)).isdel==.false.) then
                    children=[children(:nchd),this.subface(i)]
                    nchd=nchd+1
                else
                    call get_child_face(face(this.subface(i)),this.subface(i),children,nchd)
                endif
            enddo            
        endif
        
    endsubroutine    
    
    function node2vgedge(inode) result(vedge)
    !Given: inode, the node id in the mother 2d model
    !Return: the first vertical edge id in the gedge,extruded by the node
        implicit none
        integer,intent(in)::inode
        integer,allocatable::vedge(:)        
        integer::i,n1,n2=0,nvedge1=0
        
        n1=nedge*(soillayer+1)+(inode-1)*soillayer+1
        n2=0
        do i=n1,n1+soillayer-1            
            call gedge(i).getedge(i,vedge,nvedge1)                
        enddo
        
        vedge=vedge(:nvedge1)
        
    end function
    
    function edge2vface(iedge) result(vface)
    !Given: iedge, the edge id in the mother 2d model
    !Return: the first vertical face id in the face,extruded by the edge
        implicit none
        integer,intent(in)::iedge
        integer,allocatable::vface(:)        
        integer::i,n1,n2=0,nvface1=0        
       
        n1=nelt*(soillayer+1)+(iedge-1)*soillayer+1
        n2=0
        do i=n1,n1+soillayer-1            
            call face(i).getface(i,vface,nvface1)                
        enddo
        
        vface=vface(:nvface1)        
        
    end function    
    
    subroutine Gen_SubElement()
       
	    implicit none
		!logical,intent(in)::istet
	    integer::i,j,k,n1,n2,n3,n4=0,iat(4)=0,e1(2),N5,v1(2)
	    real(8)::t1,t2
	    !type(element_tydef),allocatable::ELT1(:),PRMELT(:),TetELT(:)
	    integer::ielt1,IA1(100)
        integer,allocatable::iw1(:)
	    !if(et1==153) call gen_6_noded_triangle
        !n2=count(elt(:nelt).isdel==.false.)
        nface=nelt*(soillayer+1)+nedge*soillayer
        allocate(selt(nelt*soillayer),face(nface))
        ALLOCATE(gadjlist(ngnode))
        !element的顺序，按elt的顺序由下往上，上下单元号相差1
        !face的顺序，先水平层，后按edge的顺序生成竖向面,相邻的竖向面相差1.
        !node的顺序，按层生成，上下相邻节点号相差nnode
        !edge的顺序，按层生成，上下相邻节点号相差nedge        
        !gen gedge
        ngedge=nedge*(soillayer+1)+nnode*soillayer
        allocate(gedge(ngedge))
        n1=0;
        !horizontal edge
        do i=0,soillayer
            do j=1,nedge
                n1=n1+1
                if(edge(j).v(1)*edge(j).v(2)==0) cycle
                gedge(n1).V=node(edge(j).v+nnode*i).iptr
                
                !不输出无厚度单元的侧边
                if(edge(j).isxyoverlap/=0) then
                    gedge(n1).isoverlap=1
                    gedge(n1).isdel=.true.
                endif
                
                !v1=node(gedge(n1).v).iptr
                !where(v1>0) gedge(n1).v=v1
                call addadjlist(gadjlist,gedge(n1).v(1),gedge(n1).v(2),n1,gedge(n1).iptr)
                if(gedge(n1).iptr/=n1) gedge(n1).isdel=.true.
                !if(all(v1>0)) then
                !    if(node(v1(1)).layer==node(v1(2)).layer) then
                !        n2=node(v1(1)).layer-node(gedge(n1).v(1)).layer
                !        gedge(n1).iptr=n1+n2*nedge
                !        gedge(n1).isdel=.true.
                !    endif    
                !endif
                
                
            enddo
        enddo
        !vertical edge
        do i=1,nnode
            do j=0,soillayer-1
                n1=n1+1
                gedge(n1).v=node([i+j*nnode,i+(j+1)*nnode]).iptr
                !where(node(gedge(n1).v).iptr>0) gedge(n1).v=node(gedge(n1).v).iptr
                if(gedge(n1).v(1)-gedge(n1).v(2)==0) then
                    gedge(n1).isoverlap=1
                    gedge(n1).isdel=.true.
                endif
                call addadjlist(gadjlist,gedge(n1).v(1),gedge(n1).v(2),n1,gedge(n1).iptr)
                if(gedge(n1).iptr/=n1) gedge(n1).isdel=.true.
                
            enddo
        enddo
        
        

        !gen face
        n1=0;N3=0
        !horizontal face
        do i=0,soillayer
            do j=1,nelt
                n1=n1+1
                if(elt(j).isdel)  cycle
                n2=elt(j).nnum
                face(n1).nnum=n2
                if(size(face(n1).node)<face(n1).nnum) then
                    call enlarge_ar(face(n1).node,10)
                    call enlarge_ar(face(n1).edge,10)
                endif
                face(n1).node(1:n2)=node(elt(j).node(1:n2)+nnode*i).iptr
                !where(node(face(n1).node(1:n2)).iptr>0) face(n1).node(1:n2)=node(face(n1).node(1:n2)).iptr
                face(n1).edge(1:n2)=gedge((elt(j).edge(1:n2)+nedge*i)).iptr
                !where(gedge(abs(face(n1).edge(1:n2))).iptr>0) face(n1).edge(1:n2)=gedge(abs(face(n1).edge(1:n2))).iptr
                face(n1).edge(1:n2)=face(n1).edge(1:n2)*elt(j).orient(1:n2)
                face(n1).ilayer=i
                face(n1).isdel=.false.
                face(n1).iptr=n1
                if(elt(j).et==-1) face(n1).isdel=.true. !不输出无厚度单元的上下面
                call edge_adjlist_add_face(n1)
                
                N3=N3+1
            enddo
        enddo
        !vertical face
        do i=1,nedge
            do j=1,soillayer
                n1=n1+1
                face(n1).isv=1;face(n1).ilayer=j
                face(n1).iptr=n1
                if(edge(i).v(1)*edge(i).v(2)==0) cycle
                if(size(face(n1).node)<face(n1).nnum) then
                    call enlarge_ar(face(n1).node,10)
                    call enlarge_ar(face(n1).edge,10)
                endif
                face(n1).node(1:4)=node([edge(i).v(1:2)+(j-1)*nnode,edge(i).v(2:1:-1)+j*nnode]).iptr
                
                face(n1).edge(1)=i+(j-1)*nedge
                face(n1).edge(3)=-(i+(j)*nedge)
                face(n1).edge(2)=(edge(i).v(2)-1)*soillayer+j+nedge*(soillayer+1)
                face(n1).edge(4)=-((edge(i).v(1)-1)*soillayer+j+nedge*(soillayer+1))
                face(n1).edge(:4)=sign(gedge(abs(face(n1).edge(:4))).iptr,face(n1).edge(:4))
                !where(node(face(n1).node(:face(n1).nnum)).iptr>0) &
                !    face(n1).node(:face(n1).nnum)=node(face(n1).node(:face(n1).nnum)).iptr
                n2=4
                if(face(n1).node(1)==face(n1).node(4)) then
                    n2=3                    
                endif
                if(face(n1).node(2)==face(n1).node(3)) then
                    n2=n2-1
                    if(n2==3) then
                        face(n1).node(3)=face(n1).node(4)                        
                    endif
                    face(n1).edge(2:3)=face(n1).edge(3:4)
                endif
                face(n1).nnum=n2
                
                if(n2>2) face(n1).isdel=.false.
                !不输出无厚度单元的侧边的拉伸面
                if(edge(i).isxyoverlap/=0) face(n1).isdel=.true.
                call edge_adjlist_add_face(n1)
                N3=N3+1
            enddo
        enddo
        !gen element
        nselt=0;N5=0
	    do i=1,nelt
            if(elt(i).isdel) cycle 
            allocate(elt(i).selt(soillayer))
            elt(i).selt=0
            ielt1=(I-1)*soillayer
	        do j=1,soillayer
                nselt=ielt1+j
                selt(nselt).nnum=2*elt(i).nnum
                if(size(selt(nselt).node)<selt(nselt).nnum) call enlarge_ar(selt(nselt).node,10)
                if(elt(i).et==0) then
                    selt(nselt).et=63
                    selt(nselt).nface=5
                endif
                if(elt(i).et==-1) then
                    selt(nselt).et=83
                    selt(nselt).nface=6
                endif
                if(size(selt(nselt).face)<selt(nselt).nface) then
                    call enlarge_ar(selt(nselt).face,10)
                    call enlarge_ar(selt(nselt).adj,10)
                endif
                
                selt(nselt).ilayer=j
                selt(nselt).zn=elt(i).zn
                selt(nselt).isdel=.false.
                
                !单元下、上底面的邻接单元 
                if(j==1) then
                    selt(nselt).adj(1)=-1;
                else
                    selt(nselt).adj(1)=ielt1+j-1;
                endif
                if(j==soillayer) then
                    selt(nselt).adj(2)=-1
                else
                    selt(nselt).adj(2)=ielt1+j+1
                endif
                
                do k=1,2
                    n2=(j+k-2)*nelt+i
                    selt(nselt).face(k)=n2
                    if(face(n2).e(1)==-1) then
                        face(n2).e(1)=nselt
                    else
                        face(n2).e(2)=nselt
                    endif
                enddo
                !单元四周邻接单元 
                n1=elt(i).nnum
                do k=1,n1
                    n2=elt(i).adj(k)
                    n3=elt(i).edge(k)
                    if(n2>0) then
                        selt(nselt).adj(2+k)=(n2-1)*soillayer+j
                        !IF(NSELT==3491) THEN
                        !    PRINT *, 'ADJ(J)',2+K,selt(nselt).adj(2+k)
                        !ENDIF
                    else
                        selt(nselt).adj(2+k)=-1 
                    endif
                    n4=(n3-1)*soillayer+j+nelt*(soillayer+1)
                    selt(nselt).face(2+k)=n4
                    if(face(n4).e(1)==-1) then
                        face(n4).e(1)=nselt
                    else
                        face(n4).e(2)=nselt
                    endif
                enddo  
                 
                iat=0
                do k=1,n1
				    selt(nselt).node(k)=node(elt(i).node(k)+(j-1)*nnode).iptr
				    !重节点（高程相等的节点）
				    !if(node(selt(nselt).node(k)).iptr/=0) &
					   ! selt(nselt).node(k)=node(selt(nselt).node(k)).iptr
							
				    selt(nselt).node(n1+k)=node(elt(i).node(k)+(j)*nnode).iptr
				    !重节点（高程相等的节点）
				    !if(node(selt(nselt).node(n1+k)).iptr/=0) & 
					   ! selt(nselt).node(n1+k)=node(selt(nselt).node(n1+k)).iptr
						
				    if(selt(nselt).node(k)==selt(nselt).node(n1+k)) then
					    iat(k)=1
				    end if						
                end do
				 
                 
                elt(i).selt(j)=nselt
                 
                if(elt(i).et==0) then
                                    
                    select case(count(iat==1)) 
                    case(1) !变成1个五面体单元
	                    selt(nselt).nnum=5
	                    selt(nselt).et=53
                        n1=maxloc(iat,dim=1)
                        n2=mod(n1,3)+1
                        n3=mod(n2,3)+1
                        selt(nselt).node(1:5)=selt(nselt).node([n3,n2,n2+3,n3+3,n1])
                        selt(nselt).node(6)=-1
                                
                    case(2) !变成一个四面体单元
								
	                    selt(nselt).nnum=4
	                    selt(nselt).et=43
                        n1=minloc(iat,dim=1)
	                    selt(nselt).node(1:3)=selt(nselt).node(1:3)
	                    selt(nselt).node(4)=selt(nselt).node(3+n1)
                        selt(nselt).node(5:6)=-1
                        n2=mod(n1,3)+1
                        !n3=mod(n2,3)+1
                        
                        !selt(nselt).adj(2+n2)=selt(nselt).adj(2+n3)
                        !selt(nselt).adj(2+n3)=-1
                        selt(nselt).face(2+n2)=-1
                        !selt(nselt).face(2+n3)=-1
                        iw1=pack(selt(nselt).adj(:selt(nselt).nface),selt(nselt).face(:selt(nselt).nface)>0)
                        selt(nselt).adj(:4)=iw1
                        iw1=pack(selt(nselt).face(:selt(nselt).nface),selt(nselt).face(:selt(nselt).nface)>0)                        
                        selt(nselt).face(:4)=iw1             
                                                
                        selt(nselt).nface=4
                        
                        
                    case(3) !dead,delete
                        
                        !selt(nselt).face(2)=face(selt(nselt).face(1)).iptr
                        !face(selt(nselt).face(2)).iptr=face(selt(nselt).face(1)).iptr
                        !
                        !if(selt(nselt).adj(2)>0) selt(selt(nselt).adj(2)).adj(1)=selt(nselt).adj(1)
                        !if(selt(nselt).adj(1)>0) selt(selt(nselt).adj(1)).adj(2)=selt(nselt).adj(2)
                        call remove_selt(nselt)
                        selt(nselt).isze=1
                        
                             
                    end select
                !else
                    
                    !print *, 'the function for et=i is to be improved.sub=Gen_SubElement(). i=', elt(i).et
                    
                end if

                 
	         end do
        end do
        
        
        
	    !假定土体厚度为0的单元。        
        do i=1,nelt
            if(elt(i).isdel) cycle
            
            
            N2=count(selt(elt(i).selt).isze==0)
            if(n2==soillayer) cycle
            
            ia1(1:n2)=pack(elt(i).selt,selt(elt(i).selt).isze==0)
            
            do j=1,n2
                n3=ia1(j)
                !n4=selt(n3).adj(1)
                if(j==1)then
                    selt(n3).adj(1)=-1
                    where(face(selt(n3).face(1)).e/=n3) face(selt(n3).face(1)).e=-1
                else  
                    if(selt(n3).adj(1)/=ia1(j-1)) then
                        selt(n3).adj(1)=ia1(j-1)
                        N4=selt(n3).face(1)
                        face(selt(n3).face(1)).isdel=.true.
                        selt(n3).face(1)=selt(ia1(j-1)).face(2)
                        face(selt(n3).face(1)).isdel=.FALSE.
                        where(face(selt(n3).face(1)).e/=ia1(j-1)) face(selt(n3).face(1)).e=n3
                    endif
                    
                endif
                !!更新face(1)的相邻单元
                !if(n4>0) then
                !    where(face(selt(n3).face(1)).e==n4) face(selt(n3).face(1)).e=selt(n3).adj(1)
                !endif 
                
                    
                !n4=selt(n3).adj(2)    
                if(j==n2)then
                    selt(n3).adj(2)=-1
                    where(face(selt(n3).face(2)).e/=n3) face(selt(n3).face(2)).e=-1
                else
                    if(selt(n3).adj(2)/=ia1(j+1)) then
                        selt(n3).adj(2)=ia1(j+1)
                        N4=selt(n3).face(2)
                        face(selt(n3).face(2)).isdel=.true.
                        selt(n3).face(2)=selt(ia1(j+1)).face(1)
                        face(selt(n3).face(2)).isdel=.FALSE.
                        where(face(selt(n3).face(2)).e/=ia1(j+1)) face(selt(n3).face(2)).e=n3
                    endif
                endif
                !!更新face(2)的相邻单元
                !if(n4>0) then
                !    where(face(selt(n3).face(2)).e==n4) face(selt(n3).face(2)).e=selt(n3).adj(2)
                !endif
            enddo
        enddo
        
        call xzone_cut_selt()
        
        
        !find faces/edges on boundary
        face.isdel=.true.
        do i=1,nselt
            if(selt(i).isdel)cycle
            face(selt(i).face(:selt(i).nface)).isdel=.false.
        enddo

        do i=1,nface
            if(face(i).isdel)  cycle
            e1=face(i).e
            
            if(all(e1<1)) THEN
                FACE(I).ISDEL=.TRUE. !无相邻单元的面，不输出，自身也死。
                call edge_adjlist_remove_face(i)
                CYCLE
            ENDIF
            IF(ALL(E1>0)) THEN
                IF(ALL(SELT(E1).ISDEL)) THEN
                    FACE(I).ISDEL=.TRUE. !相邻两边单元都死的面，自身也死。
                    call edge_adjlist_remove_face(i)
                    CYCLE
                ENDIF
            ENDIF
            
            
            if(any(e1<0)) then
                face(i).isb=1
            elseif(selt(e1(1)).zn/=selt(e1(2)).zn) then
                face(i).isb=2 
            elseif(selt(e1(1)).zn==selt(e1(2)).zn.and.(selt(e1(1)).ilayer/=selt(e1(2)).ilayer)) then
                face(i).isb=3
            endif
        enddo
        
        
        DO j=1,NXZONE
        !MARK THE FACES ON THE CUT FACE
            IF(XZONE(J).ISX/=2) THEN
                WHERE(FACE(XZONE(j).CUTFACE).ISB==0) FACE(XZONE(j).CUTFACE).ISB=4
            ELSE
                WHERE(GEDGE(XZONE(j).CUTFACE).ISB==0) GEDGE(XZONE(j).CUTFACE).ISB=4
            ENDIF                
        ENDDO
        
        do i=1,nface
            if(face(i).isdel)  cycle
            n1=face(i).nnum
            if(face(i).isb>0) then
                gedge(abs(face(i).edge(1:n1))).isb=face(i).isb
                node(face(i).node(1:n1)).isb=face(i).isb
                call bface2zonebface(i)
            endif
        end do
        
        DO I=1,ZNUM
            DO J=1,SOILLAYER
                CALL VOLUME_NUMBER_INSIDE_ZONE(I,J)
            ENDDO
        ENDDO
        

        
    contains
    
        subroutine bface2zonebface(iface)
            implicit none
            integer,intent(in)::iface
            integer::i,elt1(2),n1,n2
            
            elt1=face(iface).e
            do i=1,2
                if(elt1(i)>0) then
                    n1=selt(elt1(i)).zn
                    n2=selt(elt1(i)).ilayer
                    if(.not.allocated(zone(n1).bface)) then
                        allocate(zone(n1).bface(soillayer))
                        do j=1,soillayer
                            allocate(zone(n1).bface(j).node(100))
                        enddo
                    endif
                    zone(n1).bface(n2).nnum=zone(n1).bface(n2).nnum+1
                    if(zone(n1).bface(n2).nnum>size(zone(n1).bface(n2).node)) &
                        call ENLARGE_AR(zone(n1).bface(n2).node,100)
                    
                    zone(n1).bface(n2).node(zone(n1).bface(n2).nnum)=iface
                endif
            enddo
        end subroutine
   
        
    end subroutine 
    
    subroutine edge_adjlist_add_face(iface)
        implicit none
        integer,intent(in)::iface
        integer::i,ie1
        
        if(face(iface).isdel) return
        
        do i=1,face(iface).nnum
            ie1=abs(face(iface).edge(i))
            if(ie1<1) cycle

            if(any(gedge(ie1).face(1:gedge(ie1).nface)==iface)) then
                return
            else
                gedge(ie1).nface=gedge(ie1).nface+1
                if(gedge(ie1).nface>size(gedge(ie1).face)) then
                    call enlarge_AR(gedge(ie1).face,10)
                    call enlarge_AR(gedge(ie1).subid,10)
                endif
                gedge(ie1).face(gedge(ie1).nface)=iface
                gedge(ie1).subid(gedge(ie1).nface)=i
            endif                
            
        enddo
        
    end subroutine
    
    subroutine edge_adjlist_remove_face(iface)
        implicit none
        integer,intent(in)::iface
        integer::i,j,ie1
        
        if(.not.face(iface).isdel) return
        
        do i=1,face(iface).nnum
            ie1=abs(face(iface).edge(i))
            if(ie1<1) cycle
            do j=1,gedge(ie1).nface
                if(gedge(ie1).face(j)==iface) then
                    gedge(ie1).face(j:gedge(ie1).nface-1)=gedge(ie1).face(j+1:gedge(ie1).nface)
                    gedge(ie1).subid(j:gedge(ie1).nface-1)=gedge(ie1).subid(j+1:gedge(ie1).nface)
                    gedge(ie1).nface=gedge(ie1).nface-1
                    exit
                endif
            enddo
        enddo
        
    end subroutine    
    
    subroutine xzone_cut_selt()
    
        implicit none
        integer::i,j,k,k1,k2,j1,n1,n2,n3,n4,n5,v1(2),nc1,nvf1
        integer::ia2(8),staterr
        real(8)::t1,t2,tol1=1.d-6,xi1(3)
        integer,allocatable::viz1(:),eiz1(:)
        
        
        do i=1,nxzone
            !InsPt1=0 !InsPt1(i)=0,need to check;-1,no intercecting point;>0,intercepting, the intersecting point is node(gedge(i).cutpoint)
            !cutedge1=0  !InsPt1(i)=0,need to check;any is -1,no intercecting point;any >0,intercepting, the intersecting point is node(gedge(i)).cutpoint              
            t1=xzone(i).te
            !if(allocated(iw1)) deallocate(iw1)
            !allocate(iw1(ngnode))
            !借用bw
            node.bw=0
            gedge.isb=0
            face.isb=0
            selt.isb=0
            !only  insert zone boundary line into the gmodel
            if(xzone(i).isx==2) then
                do j=1,zone(xzone(i).izone).nbn
                    eiz1=node2vgedge(zone(xzone(i).izone).bnode(j))
                    
                    do k=1,size(eiz1)
                        !if(gedge(k).isdel) cycle
                        call cut_edge_by_z(eiz1(k),t1)
                    enddo
                
                enddo
                
                do j=1,zone(xzone(i).izone).nbe
                    eiz1=edge2vface(zone(xzone(i).izone).bedge(j))
                    
                    do k=1,size(eiz1)
                        !if(face(k).isdel) cycle
                        call cut_face2(eiz1(k),t1)
                    enddo
                
                enddo
                
                !edge around the cut face
                xzone(i).cutface=pack([1:ngedge],gedge(1:ngedge).isb==2)
                xzone(i).ncutf=size(xzone(i).cutface)                
                
                cycle
            endif
            
            do j=1,nselt
                if(selt(j).isdel) cycle
                !if(j>old_nselt1) cycle
                if(selt(j).isdel.or.selt(j).zn/=xzone(i).izone) cycle
                n1=selt(j).nnum
                if(maxval(node(selt(j).node(1:n1)).z)<=t1-tol1) cycle                
                if(xzone(i).isx==1.and.minval(node(selt(j).node(1:n1)).z)>=t1-tol1) then
                    selt(j).isdel=.true.
                    call remove_selt(j)
                    cycle
                endif
                
                selt(j).isb=1
            enddo
            
            viz1=pack([1:nselt],selt(1:nselt).isb>0)             
            do j=1,size(viz1)
                !mark zone node                
                node(selt(viz1(j)).node(:selt(viz1(j)).nnum)).bw=1
            enddo
            
            eiz1=pack([1:ngnode],node.bw>0)
            where(abs(node(eiz1).z-t1)<tol1.and.node(eiz1).iptr==eiz1) node(eiz1).bw=2 !bw=2 node on the cut face
            
            do j1=1,size(viz1)
                !mark zone edge 
                j=viz1(j1)
                
                do k=1,selt(j).nface
                    n1=selt(j).face(k)
                    if(face(n1).isb/=0) cycle
                    face(n1).isb=-1 !no need to check again
                    n2=face(n1).nnum
                    n4=0
                    do k1=1,n2
                        n3=abs(face(n1).edge(k1))
                        if(gedge(n3).isb/=0) then
                            if(gedge(n3).isb==1) then
                                face(n1).isb=1
                            elseif(gedge(n3).isb==2) then
                                n4=n4+1
                            endif
                        else
                            if(all(node(gedge(n3).v).bw==2))then
                                gedge(n3).isb=2 !on the cut face
                                n4=n4+1
                            else
                                if((node(gedge(n3).v(1)).z-t1)*(node(gedge(n3).v(2)).z-t1)<0)then
                                    gedge(n3).isb=1 !cross the cut face
                                    face(n1).isb=1 !cross the cut face
                                else
                                    gedge(n3).isb=-1 !no need to check again
                                endif
                            endif
                        endif
                    enddo
                    if(n4==n2) face(n1).isb=2 !on the cut face
                    
                    !if(n4==n2) then
                    !    face(n1).isb=2 !on the cut face
                    !else   
                    !    !if((maxval(node(face(n1).node(:n2)).z)-t1)*(minval(node(face(n1).node(:n2)).z)-t1)<0) then
                    !    !    face(n1).isb=1
                    !    !else
                    !    !    face(n1).isb=-1 !no need to check again
                    !    !endif
                    !endif
                enddo
            enddo
            

            
            viz1=pack([1:ngedge],gedge(:ngedge).isb==1)
            !cut edge
            do j=1,size(viz1)
                n1=viz1(j)
                call cut_edge_by_z(n1,t1)
                !v1=gedge(n1).v
                !!if((node(v1(1)).z-t1)*(node(v1(2)).z-t1)<0.d0) then
                !t2=(t1-node(v1(1)).z)/(node(v1(2)).z-node(v1(1)).z)
                !
                !xi1(1)=node(v1(1)).x+t2*(node(v1(2)).x-node(v1(1)).x)
                !xi1(2)=node(v1(1)).y+t2*(node(v1(2)).y-node(v1(1)).y)
                !xi1(3)=t1
                !!如果交点与节点距离很小，则令节点处于切割面上。
                !do k=1,2
                !    if(norm2([node(v1(k)).x,node(v1(k)).y,node(v1(k)).z]-xi1)*xyscale<1.0d-3) then
                !        node(v1(k)).bw=2 !mark
                !        node(v1(k)).x=xi1(1)
                !        node(v1(k)).y=xi1(2)
                !        node(v1(k)).z=xi1(3)
                !        exit
                !    endif
                !enddo
                !if(k>2) then
                !    ngnode=ngnode+1
                !    if(size(node)<ngnode) call enlarge_ar(node,100)
                !    node(ngnode).x=xi1(1);node(ngnode).y=xi1(2);node(ngnode).z=xi1(3);
                !    node(ngnode).bw=2 !mark
                !    node(ngnode).iptr=ngnode
                !    gedge(n1).cutpoint=ngnode
                !    call cut_edge(n1,ngnode)
                !endif
                !!endif
            enddo
            !cut face
            viz1=pack([1:nface],face(1:nface).isb==1) 
            
            do j=1,size(viz1)
                n1=viz1(j)
                call cut_face2(n1,t1)
            enddo            
            
            
            viz1=pack([1:nselt],selt(1:nselt).isb>0) 
             do j=1,size(viz1)
                n1=viz1(j)
                call cut_selt2(N1,t1,xzone(i).isx)
             enddo
             
             !face on the cut face
             xzone(i).cutface=pack([1:nface],face(1:nface).isb==2)
             xzone(i).ncutf=size(xzone(i).cutface)
            
        enddo
        
        !还
        node.bw=0
        gedge.isb=0
        face.isb=0
        selt.isb=0
        
        deallocate(viz1,eiz1,stat=staterr)
    
    endsubroutine    
    
    
    subroutine cut_edge_by_z(iedge,z)
       
        implicit none
        integer,intent(in)::iedge
        real(8),intent(in)::z
         real(8)::v1(2),t2,xi1(3)
        
        if(gedge(iedge).isdel) return
        v1=gedge(iedge).v
        
        !平行xy平面的edge
        if(abs(node(v1(1)).z-node(v1(2)).z)<1.d-6) then
            if(abs(node(v1(1)).z-z)<1.d-6) then
                node(v1).bw=2
                !if(gedge(iedge).cutpoint>0) node(gedge(iedge).cutpoint).bw=2
            endif
            return
        endif
        
        t2=(z-node(v1(1)).z)/(node(v1(2)).z-node(v1(1)).z)
        
        if(t2<0.d0.or.t2>1.d0) return !no cross
        
        !if(gedge(iedge).cutpoint>0) then
        !               
        !    t2=(z-node(v1(1)).z)/(node(gedge(iedge).cutpoint).z-node(v1(1)).z)
        !    if(t2>1) then
        !        iedge1=gedge(iedge).subedge(2)
        !    else
        !        iedge1=gedge(iedge).subedge(1)
        !    endif
        !    v1=gedge(iedge1).v 
        !    t2=(z-node(v1(1)).z)/(node(v1(2)).z-node(v1(1)).z)
        !    
        !else
        !    iedge1=iedge
        !endif
            
        
        
        
        
        !if(t2<0.d0.or.t2>1.d0) return !no cross
        !cross at end nodes
        if(t2<1.d-6) then
            node(v1(1)).bw=2  
            return
        elseif(1.0d0-t2<1.d-6) then
            node(v1(2)).bw=2
            return
        endif
        
        xi1(1)=node(v1(1)).x+t2*(node(v1(2)).x-node(v1(1)).x)
        xi1(2)=node(v1(1)).y+t2*(node(v1(2)).y-node(v1(1)).y)
        xi1(3)=z

        ngnode=ngnode+1
        if(size(node)<ngnode) call enlarge_ar(node,100)
        node(ngnode).x=xi1(1);node(ngnode).y=xi1(2);node(ngnode).z=xi1(3);
        node(ngnode).bw=2 !mark
        node(ngnode).iptr=ngnode
        gedge(iedge).cutpoint=ngnode
        call cut_edge(iedge,ngnode)
   
    
    
    endsubroutine
    
        
        subroutine cut_selt2(ielt,be,isx)
            implicit none
            integer,intent(in)::ielt,isx
            real(8),intent(in)::be
            integer::edge1(10),nf1,i,j,node1(11),n1,ischeck1(10),if1,cf1,cf2,n2
            integer,allocatable::iw1(:),iw2(:),eface1(:)
            
            nf1=selt(ielt).nface
            if(any(face(selt(ielt).face(:nf1)).isb==2)) return
            
            n1=0
            do i=1,nf1
                if1=selt(ielt).face(i)                
                do j=1,face(if1).nnum
                    n2=abs(face(if1).edge(j))
                    if(gedge(n2).isb==2) then
                        if(any(edge1(:n1)==n2)) cycle
                        n1=n1+1
                        edge1(n1)=gedge(n2).iptr
                    endif
                enddo
            enddo
           
            
            if(n1>2) then            
           
                !sorted the edges
                
                nface=nface+1
                if(nface>size(face)) call GM_ENLARGE_AR(face,100)
                if(size(face(nface).node)<face(nface).nnum) then
                    call enlarge_ar(face(nface).node,10)
                    call enlarge_ar(face(nface).edge,10)
                endif
                ischeck1=0
                ischeck1(1)=1
                face(nface).edge(1)=edge1(1)
                node1(1:2)=gedge(edge1(1)).v
                do i=2,n1
                    do j=2,n1
                        if(ischeck1(j)/=0) cycle
                        if(node1(i)==gedge(edge1(j)).v(1)) then
                            node1(i+1)=gedge(edge1(j)).v(2)
                            ischeck1(j)=1
                            face(nface).edge(i)=edge1(j)
                            exit
                        elseif(node1(i)==gedge(edge1(j)).v(2)) then
                            node1(i+1)=gedge(edge1(j)).v(1)
                            ischeck1(j)=1
                            face(nface).edge(i)=-edge1(j)
                            exit
                        endif
                    enddo
                enddo
            
                if(node1(1)/=node1(n1+1)) then
                    print *, 'edges are not closed. sub=cut_selt'
                    error stop
                endif

                face(nface).isdel=.false.
                face(nface).nnum=n1
                face(nface).isb=2 !on the cut face
                if(size(face(nface).node)<face(nface).nnum) then
                    call enlarge_ar(face(nface).node,10)
                    call enlarge_ar(face(nface).edge,10)
                endif
                face(nface).node(:n1)=node1(1:n1)
                face(nface).iptr=nface
                !face(nface).edge(:n1)=sign(edge1(abs(ischeck1(:n1))),ischeck1(:n1))
                call edge_adjlist_add_face(nface)
                
                selt(ielt).cutface=nface
                selt(ielt).subelt=[nselt+1:nselt+2]
                
                if(nselt+2>size(selt)) call GM_ENLARGE_AR(selt,100)
                selt(nselt+1:nselt+2)=selt(ielt)
                selt(nselt+1:nselt+2).cutface=0
                selt(nselt+1:nselt+2).nnum=0
                selt(nselt+1:nselt+2).nface=0
                !all element faces
                cf1=0;cf2=0
                do i=1,nf1
                    if1=selt(ielt).face(i)

                    if(any(node(face(if1).node(:face(if1).nnum)).z<be-1.d-6)) then
                        cf2=cf2+1
                        if(size(selt(nselt+2).face)<cf2+1) then
                            call enlarge_ar(selt(nselt+2).face,5)
                            call enlarge_ar(selt(nselt+2).adj,5)
                        endif
                        selt(nselt+2).face(cf2)=if1
                        n2=selt(ielt).adj(i)
                        selt(nselt+2).adj(cf2)=n2
                        if(n2>0) then
                            where(selt(n2).face(:selt(n2).nface)==if1) selt(n2).adj(:selt(n2).nface)=nselt+2
                        endif
                    else
                        cf1=cf1+1
                        if(size(selt(nselt+1).face)<cf1+1) then
                            call enlarge_ar(selt(nselt+1).face,5)
                            call enlarge_ar(selt(nselt+1).adj,5)
                        endif                        
                        selt(nselt+1).face(cf1)=if1
                        n2=selt(ielt).adj(i)
                        selt(nselt+1).adj(cf1)=n2
                        if(n2>0) then
                            where(selt(n2).face(:selt(n2).nface)==if1) selt(n2).adj(:selt(n2).nface)=nselt+1
                        endif
                    endif
                    
                 
                enddo
                
                selt(nselt+1).face(cf1+1)=nface
                selt(nselt+2).face(cf2+1)=nface
                selt(nselt+1).nface=cf1+1
                selt(nselt+2).nface=cf2+1
                
                selt(nselt+1).adj(cf1+1)=nselt+2
                selt(nselt+2).adj(cf2+1)=nselt+1
                
                !子单元的节点顺序及单元类型在此还处于凌乱状态。
                iw1=pack(selt(ielt).node(:selt(ielt).nnum),node(selt(ielt).node(:selt(ielt).nnum)).z>be+1.d-6)
                selt(nselt+1).nnum=size(iw1)+face(nface).nnum
                if(size(selt(nselt+1).node)<selt(nselt+1).nnum) call enlarge_ar(selt(nselt+1).node,10)
                selt(nselt+1).node(:selt(nselt+1).nnum)=[iw1,face(nface).node(:face(nface).nnum)]
                selt(nselt+1).et=631
                iw1=pack(selt(ielt).node(:selt(ielt).nnum),node(selt(ielt).node(:selt(ielt).nnum)).z<be-1.d-6)
                selt(nselt+2).nnum=size(iw1)+face(nface).nnum
                if(size(selt(nselt+2).node)<selt(nselt+2).nnum) call enlarge_ar(selt(nselt+2).node,10)
                selt(nselt+2).node(:selt(nselt+2).nnum)=[iw1,face(nface).node(:face(nface).nnum)]
                selt(nselt+2).et=632                
                
                !update face 邻接关系
                
                do i=1,2
                    where(face(selt(nselt+1).face(:cf1)).e(i)==ielt) face(selt(nselt+1).face(:cf1)).e(i)=nselt+1
                    where(face(selt(nselt+2).face(:cf2)).e(i)==ielt) face(selt(nselt+2).face(:cf2)).e(i)=nselt+2
                enddo
                face(nface).e=[nselt+1,nselt+2]                
                if(isx==1) call remove_selt(nselt+1) 
                call remove_selt(ielt)
                nselt=nselt+2
                
            endif
            
            
        endsubroutine         
        

        subroutine cut_face2(iface,be)
            implicit none
            integer,intent(in)::iface
            real(8),intent(in)::be
            integer::j,ip1,ip2,eface1(10),nv1,ie1,i,vface1(10),nvf1,n1
            logical::tof1
            integer,allocatable::iw1(:),iw2(:)
            
            
            if(face(iface).isdel) return 
            
            nv1=face(iface).nnum
            vface1(:nv1)=face(iface).node(:nv1)
            
            iw1=pack(vface1(:nv1),node(vface1(:nv1)).bw==2)
            n1=size(iw1)
            if(n1/=2) then
                !print *,'the number of intersecting points between two faces should be 2.but it was ',n1
                !stop
            !else
                !if(n1>2) return
                return
            endif
            
            
            !the cut line is one edge,assuming the face is convex.
            tof1=.false.
            do i=1,nv1
                n1=abs(face(iface).edge(i))
                if(all(node(gedge(n1).v).bw==2)) then
                    gedge(gedge(n1).iptr).isb=2
                    tof1=.true.
                endif
            enddo
            if(tof1) return
            !if((minval(node(vface1(:nv1)).z)-be)*(maxval(node(vface1(:nv1)).z)-be)>=0.d0) return
            
            ngedge=ngedge+1
            if(ngedge>size(gedge)) call gm_enlarge_ar(gedge,100)
            face(iface).cutedge=ngedge
            gedge(ngedge).v=iw1
            gedge(ngedge).nface=0
            gedge(ngedge).isb=2 !on the cut face
            call addadjlist(gadjlist,gedge(ngedge).v(1),gedge(ngedge).v(2),ngedge,gedge(ngedge).iptr)
            if(gedge(ngedge).iptr/=ngedge) gedge(ngedge).isdel=.true.
            if(.not.allocated(gedge(ngedge).face)) allocate(gedge(ngedge).face(10),gedge(ngedge).subid(10))
            gedge(ngedge).face=[nface+1,nface+2]
            gedge(ngedge).nface=2
            
            !把一个面分成2面
             
            if(nface+3>size(face)) call GM_ENLARGE_AR(face,100)
            face(nface+1:nface+2)=face(iface)
            face(nface+1).iptr=nface+1
            face(nface+2).iptr=nface+2
            face(nface+1:nface+2).isdel=.false.
            ip1=minloc(abs(vface1(:nv1)-iw1(1)),dim=1)
            ip2=minloc(abs(vface1(1:nv1)-iw1(2)),dim=1)
            if(ip1>ip2) then
                iw1=[ip1:nv1,1:ip2]
                iw2=[ip2:ip1]
            else
                iw1=[ip1:ip2]
                iw2=[ip2:nv1,1:ip1]
            endif
            
            eface1(:nv1)=face(iface).edge(:nv1)                    
            face(nface+1).nnum=size(iw1)
            if(size(face(nface+1).node)<face(nface+1).nnum) then
                    call enlarge_ar(face(nface+1).node,10)
                    call enlarge_ar(face(nface+1).edge,10)
            endif
            face(nface+1).node(1:size(iw1))=vface1(iw1)
            face(nface+1).edge(1:size(iw1))=eface1(iw1)
            face(nface+1).edge(face(nface+1).nnum)=-gedge(ngedge).iptr
            gedge(ngedge).subid(1)=face(nface+1).nnum
                                
            face(nface+2).nnum=size(iw2)
            if(size(face(nface+2).node)<face(nface+2).nnum) then
                call enlarge_ar(face(nface+2).node,10)
                call enlarge_ar(face(nface+2).edge,10)
            endif            
            face(nface+2).node(1:size(iw2))=vface1(iw2)
            face(nface+2).edge(1:size(iw2))=eface1(iw2)
            face(nface+2).edge(face(nface+2).nnum)=gedge(ngedge).iptr
            gedge(ngedge).subid(2)=face(nface+2).nnum
            
            !令第一个子面nface+1为上子面。
            !if(any(node(vface1(iw1)).z<be)) then
            !    face(nface+3)=face(nface+1)
            !    face(nface+1)=face(nface+2)
            !    face(nface+2)=face(nface+3)
            !endif
            face(nface+1:nface+3).cutedge=0                 
            call edge_adjlist_add_face(nface+1)  
            call edge_adjlist_add_face(nface+2)
            !update face element
            do i=1,2
                ip1=face(iface).e(i)
                if(ip1>0) then
                    do j=1,selt(ip1).nface
                        if(selt(ip1).face(j)==iface) then
                            selt(ip1).face(j)=nface+1
                            selt(ip1).nface=selt(ip1).nface+1
                            if(size(selt(ip1).face)<selt(ip1).nface) then
                                call enlarge_ar(selt(ip1).face,5)
                                call enlarge_ar(selt(ip1).adj,5)
                            endif
                            selt(ip1).face(selt(ip1).nface)=nface+2
                            selt(ip1).adj(selt(ip1).nface)=selt(ip1).adj(j)
                            exit
                        endif
                    enddo
                endif
            enddo  
            
            face(iface).subface=[nface+1:nface+2]
            nface=nface+2
            face(iface).isdel=.true.
            call edge_adjlist_remove_face(iface)
            
        endsubroutine       
        
    
    subroutine remove_selt(ielt)
        implicit none
        integer,intent(in)::ielt
        integer::i,if1
        
        selt(ielt).isdel=.true.
        do i=1,selt(ielt).nface
            if1=selt(ielt).face(i)
            if(if1>0) where(face(if1).e==ielt) face(if1).e=-1
            if1=selt(ielt).adj(i)
            if(if1>0) where(selt(if1).adj==ielt) selt(if1).adj=-1
        enddo
    
    end subroutine
    
    subroutine cut_edge(iedge,inode)
        implicit none
        integer,intent(in)::iedge,inode
        integer::ev1(2)
        
        gedge(iedge).cutpoint=inode
        ev1=gedge(iedge).v
        ngedge=ngedge+1
        if(ngedge+2>size(gedge)) call GM_ENLARGE_AR(gedge,100)
        gedge(ngedge:ngedge+1)=gedge(iedge)
        gedge(ngedge:ngedge+1).cutpoint=0
        gedge(ngedge).v=[ev1(1),inode]
        call addadjlist(gadjlist,gedge(ngedge).v(1),gedge(ngedge).v(2),ngedge,gedge(ngedge).iptr)
        if(gedge(ngedge).iptr/=ngedge) gedge(ngedge).isdel=.true.
        gedge(ngedge+1).v=[inode,ev1(2)]        
        gedge(iedge).subedge(1:2)=[ngedge,ngedge+1]
        ngedge=ngedge+1
        call addadjlist(gadjlist,gedge(ngedge).v(1),gedge(ngedge).v(2),ngedge,gedge(ngedge).iptr)
        if(gedge(ngedge).iptr/=ngedge) gedge(ngedge).isdel=.true.
        call update_faces_sharing_the_edge(iedge,inode)
        
    endsubroutine
    
    subroutine update_faces_sharing_the_edge(iedge,inode)
        implicit none
        integer,intent(in)::iedge,inode
        integer::i,j,if1,ie1,ie2,nc1,elt1(20),nelt1=0,ielt1

        elt1=0;nelt1=0
        do i=1,gedge(iedge).nface
            if1=gedge(iedge).face(i)
            ie1=gedge(iedge).subid(i)
            !update face node and edge
            nc1=face(if1).nnum+1
            if(nc1>size(face(if1).node)) then
                call enlarge_ar(face(if1).node,5)
                call enlarge_ar(face(if1).edge,5)
            endif
            face(if1).node(:nc1)=[face(if1).node(1:ie1),inode,face(if1).node(ie1+1:)]
            if(face(if1).edge(ie1)>0) then
                face(if1).edge(:nc1)=[face(if1).edge(1:ie1-1),gedge(iedge).subedge,face(if1).edge(ie1+1:face(if1).nnum)]
            else
                face(if1).edge(:nc1)=[face(if1).edge(1:ie1-1),-gedge(iedge).subedge(2:1:-1),face(if1).edge(ie1+1:face(if1).nnum)] 
            endif
            face(if1).nnum=nc1
            !update subid
            do j=ie1,nc1
                ie2=abs(face(if1).edge(j))
                if(j<ie1+2) then
                    where(gedge(ie2).face(:gedge(ie2).nface)==if1) & 
                        gedge(ie2).subid(:gedge(ie2).nface)=j                    
                else 
                    where(gedge(ie2).face(:gedge(ie2).nface)==if1) & 
                        gedge(ie2).subid(:gedge(ie2).nface)=gedge(ie2).subid(:gedge(ie2).nface)+1
                endif
                
            enddo
            !fine element sharing the edge
            do j=1,2
                ielt1=face(if1).e(j)
                if(ielt1>0) then
                    if(any(elt1(:nelt1)==ielt1)) cycle
                    nelt1=nelt1+1
                    elt1(nelt1)=ielt1
                endif
            enddo
        enddo
        !update element sharing the edge
        do i=1,nelt1
            selt(elt1(i)).nnum=selt(elt1(i)).nnum+1
            if(size(selt(elt1(i)).node)<selt(elt1(i)).nnum) call enlarge_ar(selt(elt1(i)).node,10)
            selt(elt1(i)).node(selt(elt1(i)).nnum)=inode
        enddo
        
        !gedge(iedge).nface=0
        gedge(iedge).isdel=.true.
        call removeadjlist(gadjlist,gedge(iedge).v(1),gedge(iedge).v(2))    
    
    endsubroutine
    
    
   
    
    
    subroutine out_volume_to_gmsh2()
		character*256 term
		integer(4)::msg,LEN1
		character(512)::outfile,OUTFILE2
		
        call st_membrance()
       
		outfile=trim(path_name)//'_geomodel.geo'
        call WRITE2GMSH_GEOFILEi2(outfile)
        
  !      OUTFILE2=trim(path_name)//'_offmodel.off'
		!CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE VOLUME GEO FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	LEN1=LEN_trim(term)
     	msg = MESSAGEBOXQQ(TRIM(term(1:LEN1)),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
		
    endsubroutine    
    
    
    SUBROUTINE WRITE2GMSH_GEOFILEi2(FILE_L)
		
		IMPLICIT NONE

		CHARACTER(LEN=*),INTENT(IN)::FILE_L
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::V1(2),allostat,METHOD1
		REAL(8)::T1,AT1(3)
		REAL(8),ALLOCATABLE::MINDIS1(:)
		INTEGER,ALLOCATABLE::IA1(:),IA2(:) 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3,CH4,CH5,CH6		
		CHARACTER(8192*3)::STR1
        
	
        
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,10)
        
        !do i=1,nnode
        !    if(node(i).havesoildata==0) cycle
        !    do j=soillayer,1,-1
        !        if(node(i).elevation(j)-node(i).elevation(j-1)<0.01) node(i).elevation(j-1)=node(i).elevation(j)-0.01          
        !    enddo
        !enddo
		!MINIMAL DISTANCE BETWEEN NODES
		!ALLOCATE(MINDIS1(NNODE))
		!MINDIS1=1.D20
		!DO I=1,NNODE			
		!	DO J=I+1,NNODE				
		!		AT1=[(NODE(I).X-NODE(J).X)*XYSCALE,(NODE(I).Y-NODE(J).Y)*XYSCALE,0.d0]
		!		T1=MIN(XYSCALE/20,MAX(NORM2(AT1),5.0)) !SET MINDIS>0.1
		!		IF (MINDIS1(I)>T1) MINDIS1(I)=T1
		!		IF (MINDIS1(J)>T1) MINDIS1(J)=T1
		!	ENDDO
		!ENDDO
		
		!DO J=0,SOILLAYER
		!	N1=NNODE*J
		!	DO I=1,NNODE
		!		N2=N1+I
  !              IF(NODE(N2).IPTR>0.OR.NODE(N2).ISB==0) CYCLE !重节点和内节点不输出
		!		INC=INCOUNT(N2)
		!		WRITE(50,20) N2,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).elevation(j),MINDIS1(I)
		!	ENDDO
		!ENDDO
        T1=XYSCALE/20
		DO I=1,NGNODE			
            IF(NODE(I).IPTR/=I.OR.NODE(I).ISB==0) CYCLE !重节点和内节点不输出
			INC=INCOUNT(I)
			WRITE(50,20) I,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).Z,T1
		ENDDO
		
		
		DO I=1,NGEDGE
			IF(Gedge(i).v(1)*Gedge(i).v(2)==0.OR.GEDGE(I).ISOVERLAP==1.OR.GEDGE(I).ISB==0.OR.GEDGE(I).ISDEL) cycle
			INC=INCOUNT(I)
            V1=GEDGE(I).V
            !DO J=1,2
            !    IF(NODE(GEDGE(I).V(J)).IPTR>0) V1(J)=NODE(GEDGE(I).V(J)).IPTR
            !ENDDO
			WRITE(50,30), I,V1            
		ENDDO			
		

		DO I=1,NFACE
			IF (FACE(I).ISDEL.OR.FACE(I).ISB==0) CYCLE 
            N1=FACE(I).NNUM
			INC=INCOUNT(I)            
			NV1=N1-1
			WRITE(50,40) I,FACE(I).EDGE(1:N1)
			WRITE(50,50) I,I
        ENDDO	
        
        METHOD1=1
        IF(METHOD1==1) THEN
		    !VOLUME 
		    N1=0;
		    DO I=1,ZNUM
            
                IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
            
			    DO J=1,SOILLAYER
                    WRITE(CH1,'(I4)') I
				    WRITE(CH2,'(I4)') J
                    CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))                
                    CH4=TRIM(ADJUSTL(CH3))//'_HFACE'
                    CH5=TRIM(ADJUSTL(CH3))//'_VFACE'
                    NV1=N1
                    DO K=1,ZONE(I).VBFACE(J).NVOL
                        IF(ZONE(I).VBFACE(J).BFACE(K).NNUM<4) CYCLE
                        IF(COUNT((SELT.ISDEL==.FALSE.).AND.(SELT.ZN==I).AND.(SELT.ILAYER==J))==0) CYCLE
				        N1=N1+1
				        INC=INCOUNT(N1)
  
                        !WRITE(CH6,'(I4)') K
                        !CH6=TRIM(ADJUSTL(CH3))//'_V'//TRIM(ADJUSTL(CH6))
                    
				    	
			            !n4=ZONE(I).NTRIE3N
                        INC3=0
                        STR1=""
				        DO K1=1,ZONE(I).VBFACE(J).BFACE(K).NNUM
                            N2=ZONE(I).VBFACE(J).BFACE(K).NODE(K1)
					        INC2=INCOUNT(N2)
					        CH1=""
					        WRITE(CH1,90) N2                                            
					        STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                            INC3=INC3+1    
                            IF(MOD(INC3,100)==0) THEN
                                STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                            ENDIF
				        ENDDO
			
				        !WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
				        !               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
                        INC2=LEN_TRIM(STR1)
                        IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				        IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
				        IF(ZONE(I).VBFACE(J).BFACE(K).NNUM>0) THEN
					        WRITE(50,61) N1,STR1(1:INC2)
					        WRITE(50,70) N1,N1
                        ENDIF
                    
                    ENDDO
                
                    IF(N1-NV1==0) CYCLE
                
                    LEN1=LEN(TRIM(ADJUSTL(CH3)))
                    IF(N1-NV1==1) THEN
                        WRITE(50,82) TRIM(CH3),N1
                    ELSEIF(N1-NV1>1) THEN
                        INC=INCOUNT(NV1+1);INC2=INCOUNT(N1)
                        WRITE(50,81) TRIM(CH3),NV1+1,N1
                    ENDIF
                
                    !volume outer horizontal faces
                    N3=ZONE(I).BFACE(J).NNUM
                    IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                            .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==0 &
                            .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ILAYER/=0)
                    INC3=0
                    STR1=""
                    IF(SIZE(IA1)>0) THEN
				        DO K=1,SIZE(IA1)
                            N2=IA1(K)
					        INC2=INCOUNT(N2)
					        CH1=""
					        WRITE(CH1,90) N2                                            
					        STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                            INC3=INC3+1    
                            IF(MOD(INC3,100)==0) THEN
                                STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                            ENDIF
                        ENDDO
                        INC2=LEN_TRIM(STR1)
                        IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				        IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                        LEN1=LEN(TRIM(ADJUSTL(CH4)))
				        IF(SIZE(IA1)>0) THEN
					        WRITE(50,100) TRIM(CH4),STR1(1:INC2)
                        ENDIF
                    ENDIF
                
                    !volume outer vertical faces
                    IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                            .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==1)
                    INC3=0
                    STR1=""
                    IF(SIZE(IA1)>0) THEN
				        DO K=1,SIZE(IA1)
                            N2=IA1(K)
					        INC2=INCOUNT(N2)
					        CH1=""
					        WRITE(CH1,90) N2                                            
					        STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                            INC3=INC3+1    
                            IF(MOD(INC3,100)==0) THEN
                                STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                            ENDIF
                        ENDDO
                        INC2=LEN_TRIM(STR1)                    
                        IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				        IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                        LEN1=LEN(TRIM(ADJUSTL(CH5)))
				        IF(SIZE(IA1)>0) THEN
					        WRITE(50,100) TRIM(CH5),STR1(1:INC2)
                        ENDIF
                    ENDIF                
			    
                ENDDO
                
            ENDDO        
            
            DO I=1,NXZONE
                IF(XZONE(I).NCUTF<1) CYCLE
                STR1=""
                DO K=1,XZONE(I).NCUTF
                    N2=XZONE(I).CUTFACE(K)
					INC2=INCOUNT(N2)
					CH1=""
					WRITE(CH1,90) N2                                            
					STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                    INC3=INC3+1    
                    IF(MOD(INC3,100)==0) THEN
                        STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                    ENDIF
                ENDDO
                INC2=INCOUNT(I)
                CH1=""
				WRITE(CH1,91) I
                IF(XZONE(I).ISX/=2) THEN
                    CH5='XZONE_'//TRIM(ADJUSTL(CH1))//'_CUTFACE'
                ELSE
                    CH5='XZONE_'//TRIM(ADJUSTL(CH1))//'_CUTEDGE'
                ENDIF
                INC2=LEN_TRIM(STR1)                    
                IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                LEN1=LEN(TRIM(ADJUSTL(CH5)))
				IF(XZONE(I).NCUTF>0) THEN
                    IF(XZONE(I).ISX/=2) THEN
					    WRITE(50,100) TRIM(CH5),STR1(1:INC2)
                    ELSE
                        WRITE(50,101) TRIM(CH5),STR1(1:INC2)
                    ENDIF
                    
                ENDIF                
                
            ENDDO            
                
            
        ELSE
        
		
		    !VOLUME 
		    N1=0;
		    DO I=1,ZNUM
            
                IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
            
			    DO J=1,SOILLAYER
  
                    IF(ZONE(I).BFACE(J).NNUM<3) CYCLE
                    IF(COUNT((SELT.ISDEL==.FALSE.).AND.(SELT.ZN==I).AND.(SELT.ILAYER==J))==0) CYCLE
				    N1=N1+1
				    INC=INCOUNT(N1)
				    WRITE(CH1,'(I4)') I
				    WRITE(CH2,'(I4)') J
				    CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))
                    CH4=TRIM(ADJUSTL(CH3))//'_HFACE'
                    CH5=TRIM(ADJUSTL(CH3))//'_VFACE'
				    LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			        !n4=ZONE(I).NTRIE3N
                    INC3=0
                    STR1=""
				    DO K=1,ZONE(I).BFACE(J).NNUM
                        N2=ZONE(I).BFACE(J).NODE(K)
					    INC2=INCOUNT(N2)
					    CH1=""
					    WRITE(CH1,90) N2                                            
					    STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                        INC3=INC3+1    
                        IF(MOD(INC3,100)==0) THEN
                            STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                        ENDIF
				    ENDDO
			
				    !WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
				    !               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
                    INC2=LEN_TRIM(STR1)
                    IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				    IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
				    IF(ZONE(I).BFACE(J).NNUM>0) THEN
					    WRITE(50,61) N1,STR1(1:INC2)
					    WRITE(50,70) N1,N1
					    WRITE(50,80) TRIM(CH3),N1,N1
                    ENDIF
                
                    !volume outer horizontal faces
                    N3=ZONE(I).BFACE(J).NNUM
                    IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                         .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==0)
                    INC3=0
                    STR1=""
                    IF(SIZE(IA1)>0) THEN
				        DO K=1,SIZE(IA1)
                            N2=IA1(K)
					        INC2=INCOUNT(N2)
					        CH1=""
					        WRITE(CH1,90) N2                                            
					        STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                            INC3=INC3+1    
                            IF(MOD(INC3,100)==0) THEN
                                STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                            ENDIF
                        ENDDO
                        INC2=LEN_TRIM(STR1)
                        IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				        IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                        LEN1=LEN(TRIM(ADJUSTL(CH4)))
				        IF(SIZE(IA1)>0) THEN
					        WRITE(50,100) TRIM(CH4),STR1(1:INC2)
                        ENDIF
                    ENDIF
                
                    !volume outer horizontal faces
                    IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                         .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==1)
                    INC3=0
                    STR1=""
                    IF(SIZE(IA1)>0) THEN
				        DO K=1,SIZE(IA1)
                            N2=IA1(K)
					        INC2=INCOUNT(N2)
					        CH1=""
					        WRITE(CH1,90) N2                                            
					        STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                            INC3=INC3+1    
                            IF(MOD(INC3,100)==0) THEN
                                STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                            ENDIF
                        ENDDO
                        INC2=LEN_TRIM(STR1)                    
                        IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				        IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                        LEN1=LEN(TRIM(ADJUSTL(CH5)))
				        IF(SIZE(IA1)>0) THEN
					        WRITE(50,100) TRIM(CH5),STR1(1:INC2)
                        ENDIF
                    ENDIF                
			
			    ENDDO
		    ENDDO
		ENDIF
		!N1=0;INC3=0
		!DO I=1,ZNUM
  !          IF(ZONE(I).OUTGMSHTYPE<2) CYCLE
		!	J=ZONE(I).IELEVATION
		!	N1=N1+1
		!	INC=INCOUNT(N1)
		!	WRITE(CH1,'(I4)') I
		!	WRITE(CH2,'(I4)') J
		!	CH3='Z'//TRIM(ADJUSTL(CH1))//'_E'//TRIM(ADJUSTL(CH2))
		!	LEN1=LEN(TRIM(ADJUSTL(CH3)))	
		!	!n4=ZONE(I).NTRIE3N
  !          STR1=""
		!	DO K=1,ZONE(I).NTRIE3N                    
		!		N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
		!		N2=ENUMBER*(J-1)+N3;
  !              
  !              IF(IA1(N2)==0) CYCLE
  !              
		!		INC2=INCOUNT(N2)
		!		CH1=""
		!		WRITE(CH1,90) N2
		!		STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
  !              INC3=INC3+1    
  !              IF(MOD(INC3,100)==0) THEN
  !                  STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
  !              ENDIF                
  !          ENDDO
  !          INC2=LEN_TRIM(STR1)
  !          IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
		!	IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','		
		!	WRITE(50,100) TRIM(CH3),N1,STR1(1:INC2)				
  !      ENDDO        
  !      
        DO I=1,NBLO
            CALL BLO(i).write(50)
        ENDDO
        
		CLOSE(50)
        
        DEALLOCATE(MINDIS1,IA1,stat=allostat)

	10  FORMAT('SetFactory("OpenCASCADE");Geometry.OCCAutoFix = 0;', /, 'meshscale=1;')
	20	FORMAT("Point(",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")
	30	FORMAT("Line(",I<INC>,")={",I7,",",I7,"};")
	40  FORMAT("Line Loop(",I<INC>,")={",<NV1>(I7,","),I7,"};")
	50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
	60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I7,","),I7,"};")
	61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
	70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
    80  FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
    81  FORMAT('Physical Volume("',A<LEN1>,'")={',I<INC>,":",I<INC2>"};")
    82	FORMAT('Physical Volume("',A<LEN1>,'")={',I<INC>,"};")
    90  FORMAT(I<INC2>,",")
    91 	FORMAT(I<INC2>)
	100 FORMAT('Physical Surface("',A<LEN1>,'"',")={",A<INC2>"};")
	101 FORMAT('Physical Line("',A<LEN1>,'"',")={",A<INC2>"};")   
	ENDSUBROUTINE	    
    
    
    subroutine blo_output(this,unit)
        class(blo_tydef)::this
        integer,intent(in)::unit
        
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::NEDGE1=0,NBCEDGE1=0
        INTEGER,ALLOCATABLE::IW1(:),IA1(:)
		REAL(8)::T1,AT1(3)
		!INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3		
		CHARACTER(8192*3)::STR1
		


		SELECT CASE(THIS.ITYPE)
            
        CASE(0,2,3) !X-Y平面内的多边形
            !point
            WRITE(UNIT,11)
            N2=0;N3=0;N4=0 
            DO I=1,THIS.NLOOP               
                DO J=1,THIS.BLOOP(I).NNUM
                    N1=THIS.BLOOP(I).NODE(J)
                    IF(J>2.AND.THIS.BLOOP(I).NNUM==J.AND.N1==THIS.BLOOP(I).NODE(1)) EXIT !跳过闭环的重复节点
                    N2=N2+1
                    INC=INCOUNT(N2)
                    !point1
                    if(this.itype==0) THEN
                        t1=THIS.BLOOP(I).ELEVATION(1)
                    ELSE
                        T1=THIS.BLOOP(I).ELEVATION(J)
                    ENDIF
                    WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,T1,10.0
                ENDDO
                DO J=1,THIS.BLOOP(i).NNUM-1 
                    N3=N3+1                    
                    IA1=[N3,N3,MOD(N3,(THIS.BLOOP(i).NNUM-1))+1+(THIS.BLOOP(i).NNUM-1)*(I-1)]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,30) IA1
                ENDDO                 
                !LINELOOP
                IA1=[I,N4+1,N4+THIS.BLOOP(I).NNUM-1]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,40) IA1
                N4=N4+THIS.BLOOP(I).NNUM-1
            ENDDO
            if(this.itype/=3) then 
                !SURFACE
                IF(THIS.NLOOP>1)    THEN
                    !SURFACE
                    IA1=[1,1,THIS.NLOOP]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,51) 1,1,THIS.NLOOP
                ELSE
                    INC=1
                    WRITE(UNIT,50) 1,1
                ENDIF
                !PHYSICAL Surface
                inc2=INCOUNT(this.iblo)
                IW1=[INC2,1]
                WRITE(UNIT,100) this.iblo,1 
            else
                !Physical line
                IA1=[this.iblo,1,sum(this.bloop.nnum) -this.nloop]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,110) IA1 
            endif
        
                    
        CASE(1,4)
        
            IF(THIS.NLOOP/=1) THEN
                PRINT *,'NLOOP SHOULD BE 1 WHEN ITYPE=1.'
                STOP
            ENDIF
            
            WRITE(UNIT,11)
            !POINT
            N2=0;N3=0;N4=0 
            DO I=1,2             
                DO J=1,THIS.BLOOP(1).NNUM 
                    !point1
                    N1=THIS.BLOOP(1).NODE(J)
                    IF(.NOT.(J>2.AND.THIS.BLOOP(1).NNUM==J.AND.N1==THIS.BLOOP(1).NODE(1))) THEN
                        !跳过闭环的重复节点
                        N2=N2+1
                        INC=INCOUNT(N2)
                        if(this.itype==1) then
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION(I),10.0
                        else
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION((I-1)*THIS.BLOOP(1).NNUM+J),10.0
                        endif 
                    ENDIF
                ENDDO
                DO J=1,THIS.BLOOP(1).NNUM-1 
                    N3=N3+1                    
                    IA1=[N3,N3,MOD(N3,(THIS.BLOOP(1).NNUM-1))+1+(THIS.BLOOP(1).NNUM-1)*(I-1)]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,30) IA1
                ENDDO                
                
                IF(THIS.ITYPE==1) THEN
                    !LINELOOP
                    IA1=[I,N4+1,N4+THIS.BLOOP(1).NNUM-1]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,40) IA1
                    N4=N4+THIS.BLOOP(1).NNUM-1
                    !SURFACE
                    INC=IW1(1)
                    WRITE(UNIT,50) I,I
                ENDIF
            ENDDO
            !VERTICAL EDGE
            N2=2*(THIS.BLOOP(1).NNUM-1)
            N4=THIS.BLOOP(1).NNUM
            IF(THIS.BLOOP(1).NNUM>2.AND.THIS.BLOOP(1).NODE(1)==THIS.BLOOP(1).NODE(THIS.BLOOP(1).NNUM))  N4=N4-1
            
            DO J=1,THIS.BLOOP(1).NNUM
                N1=THIS.BLOOP(1).NODE(J)
                IF(J>2.AND.THIS.BLOOP(1).NNUM==J.AND.N1==THIS.BLOOP(1).NODE(1)) EXIT !跳过闭环的重复节点
                N2=N2+1
                IA1=[N2,J,J+N4]
                IW1=INCOUNT(IA1)                              
                WRITE(UNIT,30) IA1    
            ENDDO
            N2=2;
            IF(THIS.ITYPE==4) N2=0
            N3=(THIS.BLOOP(1).NNUM-1)
            !VERTICAL SURFACE
            DO J=1,THIS.BLOOP(1).NNUM-1
                N2=N2+1
                IA1=[N2,J,2*N3+MOD(J,N3)+1,(N3+J),(2*N3+J)]
                IW1=INCOUNT(IA1)
                !LINE LOOP
                WRITE(UNIT,41) IA1
                !SURFACE 
                INC=IW1(1)
                WRITE(UNIT,50) N2,N2
            ENDDO 
            IF(THIS.ITYPE==1) THEN
                !VOLUME
                !SURFACE LOOP
                IA1=[1,1,N2]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,60) IA1
                !VOLUME
                INC=IW1(1)
                WRITE(UNIT,70) 1,1
                !PHYSICAL VOLUME                
                inc2=INCOUNT(this.iblo)
                IW1=[INC2,1]
                WRITE(UNIT,80) this.iblo,1
            ELSE
                !PHYSICAL Surface
                IA1=[this.iblo,1,n2]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,101)  IA1
            ENDIF
            
        
            
        END SELECT
        


10      FORMAT('SetFactory("OpenCASCADE");', /, 'meshscale=1;')
11      FORMAT('p1=newp-1;e1=newl-1;lp1=newll-1;s1=news-1;sp1=newsl-1;v1=newv;')           
20      FORMAT("Point(p1+",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")       
30      FORMAT("Line(e1+",I<IW1(1)>,")={p1+",I<IW1(2)>,",p1+",I<IW1(3)>,"};")       
40      FORMAT("Line Loop(lp1+",I<IW1(1)>,")={e1+",I<IW1(2)>,":e1+",I<IW1(3)>"};")
41      FORMAT("Line Loop(lp1+",I<IW1(1)>,")={","e1+",I<IW1(2)>,",e1+",I<IW1(3)>,",-(e1+",I<IW1(4)>,"),-(e1+",I<IW1(5)>")};")
50      FORMAT("Plane Surface(s1+",I<INC>,")={lp1+",I<INC>,"};")
51      FORMAT("Plane Surface(s1+",I<IW1(1)>,")={lp1+",I<IW1(2)>,":lp1+",I<IW1(3)>,"};")       
60      FORMAT("Surface Loop(sp1+",I<IW1(1)>,")={s1+",I<IW1(2)>,":s1+",I<IW1(3)>,"};")
61	    FORMAT("Surface Loop(sp1+",I<IW1(1)>,")={s1+",A<IW1(2)>"};")
70	    FORMAT("Volume(v1+",I<INC>,")={sp1+",I<INC>,"};")
80	    FORMAT('Physical Volume("BLO_VOL_',I<IW1(1)>,'"',")={v1+",I<IW1(2)>,"};")
90 	    FORMAT(I<INC2>,",")
100     FORMAT('Physical Surface("BLO_FACE',I<IW1(1)>,'"',")={s1+",I<IW1(2)>"};")
101     FORMAT('Physical Surface("BLO_FACE_',I<IW1(1)>,'"',")={s1+",I<IW1(2)>,":s1+",I<IW1(3)>,"};")        
110     FORMAT('Physical Line("BLO_LINE_',I<IW1(1)>,'"',")={e1+",I<IW1(2)>,":e1+",I<IW1(3)>,"};")
  
 
	ENDSUBROUTINE
    
    subroutine blo_readin(this,unit,ID)
        class(blo_tydef)::this
        integer,intent(in)::unit,ID
        INTEGER::I,J,N1,DN=0,NINC1
        INTEGER,PARAMETER::DNMAX=1000
        REAL(8)::AR(DNMAX) 
        
 	    call strtoint(unit,ar,dnmax,dn,dnmax)
	    !n1=I
        THIS.nloop=int(ar(1)) 
        !THIS.ndim=int(ar(2))
        THIS.itype=int(ar(2))
        THIS.IBLO=ID
        allocate(this.bloop(this.nloop))
        do i=1,this.nloop
            call strtoint(unit,ar,dnmax,dn,dnmax)
            n1=int(ar(1))
            this.bloop(i).nnum=n1
            call strtoint(unit,ar,dnmax,dn,dnmax)
            this.bloop(i).node=int(ar(1:n1))
            call strtoint(unit,ar,dnmax,dn,dnmax)
            if(this.itype/=4) then
                this.bloop(i).elevation=(ar(1:n1)) 
            else
                this.bloop(i).elevation=(ar(1:2*n1)) 
            endif
            
        enddo
        
     
    endsubroutine

   
    
    subroutine out_volume_to_gmsh()
		character*256 term
		integer(4)::msg
		character(512)::outfile,OUTFILE2
		
        call st_membrance()
        call find_zone_boudary()
        call node_on_boundary()        
        CALL CHECK_ORIENT()
		outfile=trim(path_name)//'_geomodel.geo'
        call WRITE2GMSH_GEOFILEi(outfile)
        
        OUTFILE2=trim(path_name)//'_offmodel.off'
		CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE VOLUME GEO FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	term=trim(term)
     	msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
		
    endsubroutine
    
	subroutine node_on_boundary()
		integer::i
		do i=1,nedge
            IF(edge(i).v(1)*edge(i).v(2)==0) cycle
            if(edge(i).iszonebc<0) cycle            
			node(edge(I).v).onbdy=1		
		enddo
        
		nbnode=count(node.onbdy==1)
		allocate(bnode(nbnode),BEDGE(NNODE),ZBEDGE(SOILLAYER,NNODE))
        ZBEDGE=0
		bnode=pack(node(1:nnode).number,node(1:nnode).onbdy==1)
	end subroutine
	
    subroutine find_zone_boudary()
		integer::izone
        integer::i,j,k,ielt,iedge,n1
		integer::bedge1(5000),nbe=0
		
        do izone=1,znum
        
            if(zone(izone).outgmshtype==2) cycle
            
		    nbe=0
            n1=zone(izone).ntrie3n
		    do i=1,zone(izone).ntrie3n
			    ielt=zone(izone).trie3n(i)			
			    do j=1,3
				    iedge=elt(ielt).edge(j)
				    if(edge(iedge).e(1)+1==0.or.edge(iedge).e(2)+1==0) then
					    nbe=nbe+1
					    bedge1(nbe)=iedge
				    else
					    if(elt(edge(iedge).e(1)).zn/=elt(edge(iedge).e(2)).zn) then
						    nbe=nbe+1
						    bedge1(nbe)=iedge
					    endif
				    endif
			    enddo
		    enddo
		    if(nbe>0) then
                IF(ALLOCATED(zone(izone).bedge)) DEALLOCATE(zone(izone).bedge)
			    allocate(zone(izone).bedge,source=bedge1(1:nbe))
                EDGE(bedge1(1:nbe)).ISZONEBC=IZONE
                ZONE(IZONE).NBE=NBE
                !DO J=1,NBE
                !    IF(EDGE(bedge1(J)).ISCEDGE==0) THEN                        
                !        PRINT *, 'SOME WRONG.'
                !    ENDIF
                !ENDDO

		    endif
	    enddo
    endsubroutine
    
    SUBROUTINE OUT_OFF_MATHER_MODEL(FILE_L)
        IMPLICIT NONE
        
        CHARACTER(LEN=*),INTENT(IN)::FILE_L
        INTEGER::I
        
        CALL GEN_OFF_MATHER_MODEL()
        
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,'(A3)') 'OFF'
        WRITE(50,'(3I7)') OFF_MATHER_MODEL.NNODE,OFF_MATHER_MODEL.NFACE,0
        DO I=1,OFF_MATHER_MODEL.NNODE
            WRITE(50,'(3(E15.7,1X))') OFF_MATHER_MODEL.NODE(:,I)
        ENDDO
        DO I=1,OFF_MATHER_MODEL.NFACE
            WRITE(50,'(4(I7,1X))') 3,OFF_MATHER_MODEL.FACE(:,I)-1
        ENDDO
        CLOSE(50)
        
    ENDSUBROUTINE
    
    SUBROUTINE GEN_OFF_MATHER_MODEL()
        IMPLICIT NONE
        
        INTEGER::I,J,N1,N2,N3
        
         
        !NNODE
        OFF_MATHER_MODEL.NNODE=NNODE*(SOILLAYER+1)
        ALLOCATE(OFF_MATHER_MODEL.NODE(3,OFF_MATHER_MODEL.NNODE))
        
		DO J=0,SOILLAYER
			N1=NNODE*J
			DO I=1,NNODE
				N2=N1+I
				OFF_MATHER_MODEL.NODE(:,N2)=[NODE(N2).X*XYSCALE+XMIN,NODE(N2).Y*XYSCALE+YMIN,NODE(N2).Z]
			ENDDO
		ENDDO
        !NFACE
        N1=COUNT(EDGE.ISZONEBC>0)
        OFF_MATHER_MODEL.NFACE=ENUMBER*(SOILLAYER+1)+N1*SOILLAYER*2
        ALLOCATE(OFF_MATHER_MODEL.FACE(3,OFF_MATHER_MODEL.NFACE))
        		!honrizontal triangle
		N1=0
		DO J=0,SOILLAYER
			N2=NNODE*J;
			DO I=1,NELT
				IF (ELT(I).ISDEL) CYCLE
				N1=N1+1			
                OFF_MATHER_MODEL.FACE(:,N1)=ELT(I).NODE(1:3)+N2
			ENDDO		
		ENDDO
		!vertial boundary rectangle
		DO J=1,SOILLAYER
			N2=NNODE*(J-1);N3=NNODE*J
			DO I=1,NEDGE
                IF(edge(i).v(1)*edge(i).v(2)==0) cycle
                IF(EDGE(I).ISZONEBC<0) CYCLE
				N1=N1+1			
				OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V+N2,EDGE(I).V(2)+N3]
                N1=N1+1			
				OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V(2:1:-1)+N3,EDGE(I).V(1)+N2]
			ENDDO		
		ENDDO	
    
	END SUBROUTINE
    
	SUBROUTINE WRITE2GMSH_GEOFILEi(FILE_L)
		
		IMPLICIT NONE

		CHARACTER(LEN=*),INTENT(IN)::FILE_L
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::NEDGE1=0,NBCEDGE1=0
		REAL(8)::T1,AT1(3)
		REAL(8),ALLOCATABLE::MINDIS1(:)
		INTEGER,ALLOCATABLE::IA1(:),IA2(:) 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3		
		CHARACTER(8192*3)::STR1
        
		

		
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,10)
        
        do i=1,nnode
            if(node(i).havesoildata==0) cycle
            do j=soillayer,1,-1
                if(node(i).elevation(j)-node(i).elevation(j-1)<0.01) node(i).elevation(j-1)=node(i).elevation(j)-0.01          
            enddo
        enddo
		!MINIMAL DISTANCE BETWEEN NODES
		ALLOCATE(MINDIS1(NNODE))
		MINDIS1=1.D20
		DO I=1,NNODE			
			DO J=I+1,NNODE				
				AT1=[(NODE(I).X-NODE(J).X)*XYSCALE,(NODE(I).Y-NODE(J).Y)*XYSCALE,0.d0]
				T1=MIN(XYSCALE/20,MAX(NORM2(AT1),5.0)) !SET MINDIS>0.1
				IF (MINDIS1(I)>T1) MINDIS1(I)=T1
				IF (MINDIS1(J)>T1) MINDIS1(J)=T1
			ENDDO
		ENDDO
		
		DO J=0,SOILLAYER
			N1=NNODE*J
			DO I=1,NNODE
				N2=N1+I
                IF(NODE(N2).IPTR/=N2) CYCLE !重节点不输出
				INC=INCOUNT(N2)
				WRITE(50,20) N2,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).elevation(j),MINDIS1(I)
			ENDDO
		ENDDO
		
		N1=0;NEDGE1=0
        NBCEDGE1=COUNT(EDGE.ISZONEBC>0)
		DO J=0,SOILLAYER
			N2=NNODE*J
			DO I=1,NEDGE
				IF(edge(i).v(1)*edge(i).v(2)==0) cycle
				IF(J==0) THEN
					NEDGE1=NEDGE1+1					
                ENDIF
				N1=N1+1
                !XY坐标重合的边也不输出，无厚度防渗墙不输出。
                IF(EDGE(I).ISXYOVERLAP==0) THEN
				    INC=INCOUNT(N1)
				    WRITE(50,30), N1,EDGE(I).V+N2
                ENDIF
			ENDDO			
		ENDDO
		!VERTICAL BOUNDARY EDGE 
		DO J=1,SOILLAYER
			N3=NNODE*(J-1);N4=NNODE*J
			DO I=1,NBNODE
				N1=N1+1
				IF(J==1) BEDGE(BNODE(I))=N1
                !0长度的边不输出
				IF((NODE(BNODE(I)+N3).IPTR==NODE(BNODE(I)+N4).IPTR)) THEN
                    ZBEDGE(J,BNODE(I))=1
                ELSE
                    INC=INCOUNT(N1)
				    WRITE(50,30), N1,BNODE(I)+N3,BNODE(I)+N4
                ENDIF
			ENDDO
		ENDDO
		
		!honrizontal triangle
		N1=0
        ALLOCATE(IA1(NELT*(SOILLAYER+1)+NBCEDGE1*SOILLAYER))
        IA1=0
		DO J=0,SOILLAYER
			N2=NEDGE1*J;
			DO I=1,NELT
				IF (ELT(I).ISDEL) CYCLE 
				N1=N1+1
                IF (ELT(I).ET==-1) CYCLE !无厚度防渗墙单元也不输出
				INC=INCOUNT(N1)
				NV1=2
				WRITE(50,40) N1,(EDGE(ELT(I).EDGE(1:ELT(I).NNUM)).NUM+N2)*ELT(I).ORIENT(1:ELT(I).NNUM)
				WRITE(50,50) N1,N1
                IA1(N1)=1
			ENDDO		
		ENDDO
		!vertial boundary rectangle
		DO J=1,SOILLAYER
			N2=NEDGE1*(J-1);N3=NBNODE*(J-1)
			DO I=1,NEDGE
                IF(edge(i).v(1)*edge(i).v(2)==0) cycle
                IF(EDGE(I).ISZONEBC<0) CYCLE                                
				N1=N1+1
                IF(J==1) EDGE(I).BRECT=N1
				
                AN1(1)=EDGE(I).NUM+N2
                
                IF(ZBEDGE(J,EDGE(I).V(2))==0) THEN
                    AN1(2)=BEDGE(EDGE(I).V(2))+N3
                    AN1(3)=-(EDGE(I).NUM+N2+NEDGE1)
                    N4=3
                ELSE
                    AN1(2)=-(EDGE(I).NUM+N2+NEDGE1)
                    N4=2
                ENDIF
                IF(ZBEDGE(J,EDGE(I).V(1))==0) THEN
                    N4=N4+1
                    AN1(N4)=-(BEDGE(EDGE(I).V(1))+N3)
                ENDIF
                IF(N4>2) THEN
                    NV1=N4-1
                    INC=INCOUNT(N1) 
				    WRITE(50,40) N1,AN1(1:N4)
				    WRITE(50,50) N1,N1
                    IA1(N1)=1 !MARK OUTPUT FACES
                ENDIF
				
			ENDDO		
		ENDDO		
		
		!VOLUME 
		N1=0;
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
			DO J=0,SOILLAYER-1
				N1=N1+1
				INC=INCOUNT(N1)
				WRITE(CH1,'(I4)') I
				WRITE(CH2,'(I4)') J+1
				CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))
				LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			    !n4=ZONE(I).NTRIE3N
                INC3=0
                STR1=""
				DO K=1,ZONE(I).NTRIE3N
                    
					N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
					DO K1=1,2
						IF(K1==1) THEN
							N2=ENUMBER*J+N3; !ENUMBER INCLUDING THE ZERO-THICKNESS ELEMENT.
						ELSE
							N2=ENUMBER*(J+1)+N3;
                        ENDIF
                        
                        IF(IA1(N2)==0) CYCLE !SKIP THE FACE NOT OUTPUT
                        
						INC2=INCOUNT(N2)
						CH1=""
						WRITE(CH1,90) N2                                            
						STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                        INC3=INC3+1    
                        IF(MOD(INC3,100)==0) THEN
                            STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                        ENDIF
					ENDDO
				ENDDO
				DO K=1,ZONE(I).NBE
					N2=NBCEDGE1*(J)+EDGE(ZONE(I).BEDGE(K)).BRECT
                    
                    IF(IA1(N2)==0) CYCLE
                    
					INC2=INCOUNT(N2)
					CH1=""
					WRITE(CH1,90) N2
					STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
                    INC3=INC3+1    
                    IF(MOD(INC3,100)==0) THEN
                        STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                    ENDIF
                ENDDO			
				!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
				!               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
                INC2=LEN_TRIM(STR1)
                IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
				IF(ZONE(I).NBE>0) THEN
					WRITE(50,61) N1,STR1(1:INC2)
					WRITE(50,70) N1,N1
					WRITE(50,80) TRIM(CH3),N1,N1
				ENDIF
			
			ENDDO
		ENDDO
		
		N1=0;INC3=0
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE<2) CYCLE
			J=ZONE(I).IELEVATION
			N1=N1+1
			INC=INCOUNT(N1)
			WRITE(CH1,'(I4)') I
			WRITE(CH2,'(I4)') J
			CH3='Z'//TRIM(ADJUSTL(CH1))//'_E'//TRIM(ADJUSTL(CH2))
			LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			!n4=ZONE(I).NTRIE3N
            STR1=""
			DO K=1,ZONE(I).NTRIE3N                    
				N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
				N2=ENUMBER*(J-1)+N3;
                
                IF(IA1(N2)==0) CYCLE
                
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                INC3=INC3+1    
                IF(MOD(INC3,100)==0) THEN
                    STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                ENDIF                
            ENDDO
            INC2=LEN_TRIM(STR1)
            IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
			IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','		
			WRITE(50,100) TRIM(CH3),N1,STR1(1:INC2)				
        ENDDO        
        
        DO I=1,NBLO
            CALL BLO(i).write(50)
        ENDDO
        
		CLOSE(50)
        
        DEALLOCATE(MINDIS1)

	10  FORMAT('SetFactory("OpenCASCADE");', /, 'meshscale=1;')
	20	FORMAT("Point(",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")
	30	FORMAT("Line(",I<INC>,")={",I7,",",I7,"};")
	40  FORMAT("Line Loop(",I<INC>,")={",<NV1>(I7,","),I7,"};")
	50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
	60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I7,","),I7,"};")
	61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
	70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
	80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
	90 	FORMAT(I<INC2>,",")
	100 FORMAT('Physical Surface("',A<LEN1>,'",',I<INC>,")={",A<INC2>"};")
    
    CONTAINS
        !FUNCTION ORIENT_EDGE()
        FUNCTION EDGE_NODE(IEDGE,ILAYER) RESULT(V)
            INTEGER,INTENT(IN)::IEDGE,ILAYER !O<=ILAYER<=SOILLAYER
            INTEGER::V(2)           
            
            V=EDGE(IEDGE).V+NNODE*(ILAYER)
        
    
        ENDFUNCTION
        FUNCTION VEDGE_NODE(INODE,ISOIL) RESULT(V)
            INTEGER,INTENT(IN)::INODE,ISOIL !ISOIL>=1
            INTEGER::V(2)           
            
            V=[INODE+NNODE*(ISOIL-1),INODE+NNODE*(ISOIL)]
        
    
        ENDFUNCTION
	ENDSUBROUTINE	 
    
    INTEGER FUNCTION INCOUNT(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL(8)::T1
    INTEGER::I,N1


    IF(N==0) THEN
        INCOUNT=1
        RETURN
    ENDIF
    T1=ABS(N)
    
    INCOUNT=INT(LOG10(T1))+1
	IF(N<0) INCOUNT=INCOUNT+1
    

    END FUNCTION

    FUNCTION INCOUNT2(N) RESULT(IW)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N(:)
    INTEGER,ALLOCATABLE::IW(:)
    REAL(8)::T1
    INTEGER::I,N1

    N1=SIZE(N)
    ALLOCATE(IW(N1))
    DO I=1,N1
        IF(N(I)==0) THEN
            IW(I)=1
            CYCLE
        ENDIF
        T1=ABS(N(I))
    
        IW(I)=INT(LOG10(T1))+1
	    IF(N(I)<0) IW(I)=IW(I)+1
    ENDDO

	END FUNCTION
    
    !SUBROUTINE CHECK_ORIENT()
    !    INTEGER::I,J,V1,V2
    !    
    !    DO I=1,NELT
    !        IF(ELT(I).ISDEL) CYCLE
    !        V2=EDGE(ELT(I).EDGE(1)).V(2)
    !        DO J=2,ELT(I).NNUM
    !            V1=EDGE(ELT(I).EDGE(J)).V(1)
    !            IF(V2/=V1)THEN
    !                ELT(I).ORIENT(J)=-1
    !                V2=EDGE(ELT(I).EDGE(J)).V(1)
    !            ELSE
    !                ELT(I).ORIENT(J)=1;V2=EDGE(ELT(I).EDGE(J)).V(2)
    !            ENDIF
    !        ENDDO
    !    ENDDO
    !    
    !ENDSUBROUTINE
    
    
    SUBROUTINE CHECK_ORIENT()
        INTEGER::I,J        
        DO I=1,NELT
            IF(ELT(I).ISDEL) CYCLE
            
            DO J=1,ELT(I).NNUM
                IF(ELT(I).NODE(J)/=EDGE(ELT(I).EDGE(J)).V(1)) THEN
                    ELT(I).ORIENT(J)=-1
                ELSE
                    ELT(I).ORIENT(J)=1                
                ENDIF
            ENDDO
        ENDDO
        
    ENDSUBROUTINE    
    
    SUBROUTINE SELT_ENLARGE_AR(AVAL,DSTEP)
        TYPE(subelt_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(subelt_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE 
    
   SUBROUTINE FACE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(FACE_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(FACE_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
   END SUBROUTINE 

   SUBROUTINE GEDGE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(GEDGE_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(GEDGE_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
   END SUBROUTINE 
   
    SUBROUTINE VOLUME_NUMBER_INSIDE_ZONE(IZONE,ILAYER)
        IMPLICIT NONE
        INTEGER,INTENT(IN)::IZONE,ILAYER
        INTEGER::ALLOSTAT,I,J,NGELT1,IELT1,STACK1(10000),N1,N2, &
            EDGE1(NFACE),EDGE2(NFACE),IELT2,NVBF1,ISCHK1(NSELT)
        INTEGER,ALLOCATABLE::GELT1(:)
        
        
        IF(.NOT.ALLOCATED(ZONE(IZONE).VBFACE)) ALLOCATE(ZONE(IZONE).VBFACE(SOILLAYER))
        ALLOCATE(ZONE(IZONE).VBFACE(ILAYER).BFACE(5))
        NVBF1=0
        GELT1=PACK([1:NSELT],SELT(:NSELT).ZN==IZONE.AND.SELT(:NSELT).ILAYER==ILAYER.AND.SELT(:NSELT).ISDEL==.FALSE.)
        NGELT1=SIZE(GELT1)
        !ALLOCATE(ISCHK1(NGELT1))
        ISCHK1(GELT1)=0
        
        DO WHILE(ANY(ISCHK1(GELT1)==0))    
            IELT1=MINLOC(ISCHK1(GELT1),DIM=1,MASK=(ISCHK1(GELT1)==0))
            STACK1(1)=GELT1(IELT1)
            N2=1
            EDGE2=0;EDGE1=0
            NVBF1=NVBF1+1
            IF(NVBF1>SIZE(ZONE(IZONE).VBFACE(ILAYER).BFACE)) CALL ENLARGE_AR(ZONE(IZONE).VBFACE(ILAYER).BFACE,5)
            DO WHILE(N2>0)
                IELT1=STACK1(N2)                
                N2=N2-1
                IF(ISCHK1(IELT1)==1) CYCLE
                ISCHK1(IELT1)=1
                DO J=1,SELT(IELT1).NFACE
                    N1=SELT(IELT1).FACE(J)
                    EDGE2(N1)=EDGE2(N1)+1
                    IF(FACE(N1).ISB==0.AND.EDGE1(N1)/=1) THEN
                        EDGE1(N1)=1                            
                        IELT2=SELT(IELT1).ADJ(J)
                        IF(IELT2>0) THEN
                            IF(SELT(IELT2).ISDEL==.FALSE..AND.ISCHK1(IELT2)/=1) THEN
                                N2=N2+1
                                STACK1(N2)=IELT2 
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
            
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).NODE=PACK([1:NFACE],EDGE2==1)
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).NNUM=SIZE(ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).NODE)
        ENDDO
        ZONE(IZONE).VBFACE(ILAYER).NVOL=NVBF1
   
        DEALLOCATE(GELT1,STAT=ALLOSTAT)
   
    ENDSUBROUTINE
    
    
end module geomodel



