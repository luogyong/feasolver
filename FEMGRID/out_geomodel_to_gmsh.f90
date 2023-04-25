module geomodel
	USE IFQWIN 
    use meshDS, only : node,nnode,elt,nelt,edge,nedge,adjlist,zone,znum, &
						bnode,nbnode,soillayer,poly3d,ismerged,iscompounded,&
						path_name,ENUMBER,xyscale,xmin,ymin,write_help,strtoint,&
                        ENLARGE_AR,adjlist_tydef,addadjlist,removeadjlist,ar2d_tydef2,&
                        iar2str,segindex,seg,zmin,ISMESHSIZE
    use ds_t
    use quicksort
    use CutoffWall,only:xzone,nxzone,cowall,ncow,vseg_pg,nvseg_pg,vface_pg,nvface_pg,node_pg,nnode_pg,vol_pg,nvol_pg
    implicit none
    private
    public BLO,NBLO,CHECK_ORIENT,Gen_SubElement,out_volume_to_gmsh2,NGNODE,out_z_at_center,out_z_at_node
    real(8),parameter::Pi=3.141592653589793
    
    INTERFACE INCOUNT
        MODULE PROCEDURE INCOUNT,INCOUNT2
    END INTERFACE

    INTERFACE GM_ENLARGE_AR
        MODULE PROCEDURE FACE_ENLARGE_AR,SELT_ENLARGE_AR,GEDGE_ENLARGE_AR
    END INTERFACE
    
    INTEGER::ngnode=0
    
    INTEGER,ALLOCATABLE::BEDGE(:),ZBEDGE(:,:) !区边界inode对应的每个节点生成的edge的整体编号, ZBEDGE(ISOILLAYER,INODE)=1表示边区边界inode对应的第isoillayer的厚度（竖向线段长度）为0;
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
        integer::et=63,nnum=6,ilayer,zn,nface,isze=0,isb=0,IGMSHVOL=0 !et=83 6面体单元;isze，是否为0厚度的土层单元
        !IGMSHVOL=单元所在的几何体编号
        integer::cutface=0,subelt(2)=0
        real(8)::vol
        integer,allocatable::node(:)
        integer,allocatable::adj(:),face(:) !下、上表面，四周
    contains
        procedure::calvol=>cal_selt_volume
    endtype
    type(subelt_tydef),allocatable::selt(:)
    INTEGER::NSELT=0,NGMSHVOL=0 !NGMSHVOL !模型几何体的数量
    
    type face_tydef
        logical::isdel=.true.
        integer::nnum=3,e(2)=-1,isb=0,isv=0,iptr,MARKER=0 
        !isb=1,zone boundary face;=2,model boundary face;=0,zone inside face;=3,same zone but differen soillayer;=4,on the cut face;=5,on the cowall; =-1,be merged.
        !isv=1,vertical face;=0,horizontal face.
        integer::ilayer,ismerged=0
        integer,allocatable::node(:)
        integer,allocatable::edge(:)
        integer::cutedge=0,subface(2)=0 
        real(8)::normal(3)=0.d0
    contains
        procedure::getface=>get_child_face
    endtype
    type(face_tydef),allocatable::face(:)
    integer::nface=0
    
    type gedge_tydef
        logical::isdel=.false.
        integer::v(2)=0
        integer::isoverlap=0,isb=0,nface=0,iptr=0,iscedge=0 !!isb=1,zone boundary edge;=2,model boundary edge;=0,zone inside edge
        integer::cutpoint=0,subedge(2)=0,mface(2)=0,it1=0
        integer,allocatable::face(:),subid(:)
    contains
        procedure::getedge=>get_child_edge 
    endtype
    type(gedge_tydef),allocatable::gedge(:)
    integer::ngedge=0

    type general_poly_face_tydef
        integer::id=0
        integer::nedge=0,npoly=0,nhole=0,oneface=0,nfp !face(oneface) is one of faces be merged into this facets.
        integer,allocatable::edge(:),floatp(:)
        real(8),allocatable::hole(:,:)
        type(ar2d_tydef2),allocatable::poly(:) !all edges of the poly facets, the edges is not in order
    contains
        procedure::set_facet=>poly_set_facet
        procedure::facet2gmsh=>tetgen_facet2gmsh
    endtype
    type(general_poly_face_tydef),allocatable::GPFACET(:)
    integer::NGPF=0
    
    type(adjlist_tydef),allocatable::gadjlist(:)
    
    contains
    
    subroutine out_z_at_node()

        implicit none
        character*256 term
        INTEGER::I,J,NC1
		integer(4)::msg,LEN1
		character(512)::outfile,OUTFILE2
		character(8)::CH1,ZNAME1(0:SOILLAYER)
        
        call st_membrance()
       
		outfile=trim(path_name)//'_nodal_z.dat'
        open(unit=60,file=outfile,status='replace')
        DO I=0,SOILLAYER
            NC1=INCOUNT(I)
            WRITE(ch1,10) I 
            ZNAME1(I)='Z'//trim(adjustl(ch1))
        ENDDO
        WRITE(60,20) ZNAME1(0:SOILLAYER)
        DO I=1,NNODE
            
            write(60,30) i,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).ELEVATION(:SOILLAYER)*XYSCALE+ZMIN
            
        ENDDO
        CLOSE(60)
  !      OUTFILE2=trim(path_name)//'_offmodel.off'
		!CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE NODE Z FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	LEN1=LEN_trim(term)
     	msg = MESSAGEBOXQQ(TRIM(term(1:LEN1)),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))	
10      FORMAT(I<NC1>)
20 FORMAT(14X,'N',X,14X,'X',X,14X,'Y',X,<1+SOILLAYER>(12X,A3,X))        
30 FORMAT(I15,X,<3+soillayer>(F15.3,X))
        
    endsubroutine  

    subroutine out_z_at_center()

        implicit none
        character*256 term
        INTEGER::I,J,NC1,k
		integer(4)::msg,LEN1
		character(512)::outfile,OUTFILE2
		character(8)::CH1,ZNAME1(0:SOILLAYER)
        real(8)::x1,y1,z1(0:soillayer)
        call st_membrance()
       
		outfile=trim(path_name)//'_center_z.dat'
        open(unit=60,file=outfile,status='replace')
        DO I=0,SOILLAYER
            NC1=INCOUNT(I)
            WRITE(ch1,10) I 
            ZNAME1(I)='Z'//trim(adjustl(ch1))
        ENDDO
        WRITE(60,20) ZNAME1(0:SOILLAYER)
        DO I=1,NELT
            if(elt(i).isdel) cycle 
            x1=sum(NODE(elt(i).node(1:3)).X)/3
            y1=sum(NODE(elt(i).node(1:3)).y)/3
            do j=0,soillayer
                z1(j)=0
                do k=1,3
                    z1(j)=z1(j)+NODE(elt(i).node(k)).ELEVATION(j)
                enddo
            enddo
            write(60,30) i,elt(i).zn,x1*XYSCALE+XMIN,Y1*XYSCALE+YMIN,z1(:SOILLAYER)*XYSCALE+ZMIN
            
        ENDDO
        CLOSE(60)
  !      OUTFILE2=trim(path_name)//'_offmodel.off'
		!CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE ELEMENTAL CENTER Z FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	LEN1=LEN_trim(term)
     	msg = MESSAGEBOXQQ(TRIM(term(1:LEN1)),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))	
10      FORMAT(I<NC1>)
20 FORMAT(9X,'IELT',X,9X,'IZONE',X,14X,'X',X,14X,'Y',X,<1+SOILLAYER>(12X,A3,X))        
30 FORMAT(I15,X,<3+soillayer>(F15.3,X))
        
    endsubroutine  
        
    
    subroutine tetgen_facet2gmsh(this,facebase,unit)
        implicit none
        class(general_poly_face_tydef)::this
        integer,intent(in)::facebase,unit
        integer::i,j,len1,n1,inc,lloop1,nseg1,edge1(100),W1
        character(len=:),allocatable::str
        !line loop
        write(unit,30)
        lloop1=0;nseg1=0
        do i=1,this.npoly
            if(this.poly(i).nnum<3) then
                nseg1=nseg1+1
                edge1(nseg1)=this.poly(i).edge(1)                
            else            
                str=iar2str(this.poly(i).edge(1:this.poly(i).nnum))
                inc=incount(lloop1)
                len1=len_trim(str)
                write(unit,40) lloop1,trim(str) 
                lloop1=lloop1+1
            endif
        enddo
        n1=THIS.ID+facebase
        W1=INCOUNT(n1)
        if(lloop1>1) then
            inc=incount(lloop1-1)            
            write(unit,50) n1,lloop1-1 
        elseif(lloop1==1) then
            write(unit,51) n1
        endif
        if(nseg1>0) then
            str=iar2str(edge1(1:nseg1))
            len1=len_trim(str)
            write(unit,60) trim(str),n1
        endif
        if(THIS.NFP>0) then
            str=iar2str(THIS.FLOATP)
            len1=len_trim(str)
            write(unit,70) trim(str),n1
        endif        
        
    30  FORMAT("lloop1=newll;")    
	40  FORMAT("Line Loop(lloop1+",I<INC>,")={",A<len1>,"};")
	50	FORMAT("Plane Surface(",I<W1>,")={lloop1:(lloop1+",I<INC>,")};")        
    51  FORMAT("Plane Surface(",I<W1>,")={lloop1};")
    60  FORMAT("Line {",A<len1>,"} In Surface {",I<W1>"};")  
    70  FORMAT("Point {",A<len1>,"} In Surface {",I<W1>"};")    
    endsubroutine
    
    subroutine vface_pg_handle()
        implicit none

        integer::iseg1,nnum1,i,j,k,err_info
        integer,allocatable::edge1(:),face1(:)
       
        do k=1,nvface_pg
        
            do i=1,vface_pg(k).ncp-1
                iseg1=segindex(vface_pg(k).cp(i),vface_pg(k).cp(i+1))
                   
                edge1=seg(iseg1).get_edge([arr_t(vface_pg(k).cp(i)).x,arr_t(vface_pg(k).cp(i)).y])
                do j=1,size(edge1)
                    face1=edge2vface(edge1(j))
                    vface_pg(k).face=[vface_pg(k).face,face1]
                enddo
            
            enddo
            vface_pg(k).nface=size(vface_pg(k).face)
            
            call vface_pg(k).getnode()
            
            do i=1,2
                vface_pg(k).bvedge=[vface_pg(k).bvedge,node2vgedge(vface_pg(k).node(i))]
            enddo
            
            gedge(vface_pg(k).bvedge).iscedge=-10 !mark it cannot be merged
        
        enddo
        
        deallocate(edge1,face1,STAT=err_info)
        
            
    end subroutine
    
    subroutine vseg_pg_handle()
    !把竖向线物理的控制点插入模型，并找出其所包含的线段
        implicit none
        integer::i,j,k,n1,n2,n3,n4
        integer,allocatable::vedge1(:),n2uz1(:),o2n1(:)
        real(8),allocatable::uz1(:)
        real(8)::zmin1,zmax1
        do i=1,nvseg_pg
            call vseg_pg(i).getnode()
            
            zmin1=node(vseg_pg(i).ip).elevation(0)
            zmax1=node(vseg_pg(i).ip).elevation(soillayer)
            
            if(vseg_pg(i).z(1)<zmin1.or.vseg_pg(i).z(vseg_pg(i).nnode)>zmax1) then
                if(vseg_pg(i).z(1)<zmin1) vseg_pg(i).z=[vseg_pg(i).z,zmin1]
                if(vseg_pg(i).z(vseg_pg(i).nnode)>zmax1) vseg_pg(i).z=[vseg_pg(i).z,zmax1]
                if(allocated(n2uz1)) deallocate(n2uz1)
                allocate(n2uz1(size(vseg_pg(i).z)))
                call quick_sort(vseg_pg(i).z,uz1,n2uz1)
                vseg_pg(i).z=uz1
                vseg_pg(i).nnode=size(uz1)
                deallocate(vseg_pg(i).node)
                allocate(vseg_pg(i).node(vseg_pg(i).nnode))
                vseg_pg(i).node=0
            endif
            
            do j=1,vseg_pg(i).nnode
                if(vseg_pg(i).z(j)<zmin1.or.vseg_pg(i).z(j)>zmax1 ) then
                    ngnode=ngnode+1
                    if(size(node)<ngnode) call enlarge_ar(node,100)
                    node(ngnode).x=node(vseg_pg(i).ip).x;node(ngnode).y=node(vseg_pg(i).ip).y;node(ngnode).z=vseg_pg(i).z(j);        
                    node(ngnode).s=0.d0
                    node(ngnode).iptr=ngnode
                    vseg_pg(i).node(j)=ngnode
                    cycle                    
                endif
                
                vedge1=node2vgedge(vseg_pg(i).ip)
                
                do k=1,size(vedge1)
                    call cut_edge_by_z(vedge1(k),vseg_pg(i).z(j),vseg_pg(i).node(j))                    
                    if(vseg_pg(i).node(j)>0) exit
                enddo
                                
            enddo
                    
            
            vedge1=node2vgedge(vseg_pg(i).ip)
            n1=0;n2=0
            if(vseg_pg(i).z(1)<zmin1) then
                n3=minloc(abs(vseg_pg(i).z-zmin1),dim=1)
            else
                n3=1
            endif
            if(vseg_pg(i).z(vseg_pg(i).nnode)>zmax1) then
                n4=minloc(abs(vseg_pg(i).z-zmax1),dim=1)
            else
                n4=vseg_pg(i).nnode
            endif          
            
            do k=1,size(vedge1)
                if(gedge(vedge1(k)).v(1)==vseg_pg(i).node(n3)) n1=k
                if(gedge(vedge1(k)).v(2)==vseg_pg(i).node(n4)) n2=k
                if(n1*n2/=0) then
                    vseg_pg(i).edge=vedge1(n1:n2)                    
                    exit
                endif
            enddo
            
            !add edges outside the model
            do k=1,n3-1
                ngedge=ngedge+1
                if(ngedge+1>size(gedge)) call GM_ENLARGE_AR(gedge,100)        
                gedge(ngedge).cutpoint=0
                gedge(ngedge).v=[vseg_pg(i).node(k),vseg_pg(i).node(k+1)]
                call addadjlist(gadjlist,gedge(ngedge).v(1),gedge(ngedge).v(2),ngedge,gedge(ngedge).iptr)
                if(gedge(ngedge).iptr/=ngedge) gedge(ngedge).isdel=.true.
                
            enddo
            if(n3>1) vseg_pg(i).edge=[ngedge-(n3-1)+1:ngedge,vseg_pg(i).edge]
            do k=n4,vseg_pg(i).nnode-1
                ngedge=ngedge+1
                if(ngedge+1>size(gedge)) call GM_ENLARGE_AR(gedge,100)        
                gedge(ngedge).cutpoint=0
                gedge(ngedge).v=[vseg_pg(i).node(k),vseg_pg(i).node(k+1)]
                call addadjlist(gadjlist,gedge(ngedge).v(1),gedge(ngedge).v(2),ngedge,gedge(ngedge).iptr)
                if(gedge(ngedge).iptr/=ngedge) gedge(ngedge).isdel=.true.
            enddo 
            if(n4<vseg_pg(i).nnode) vseg_pg(i).edge=[vseg_pg(i).edge,ngedge-(vseg_pg(i).nnode-n4)+1:ngedge]
            allocate(vseg_pg(i).igmshvol(size(vseg_pg(i).edge)))
            vseg_pg(i).igmshvol=0
            gedge(vseg_pg(i).edge).iscedge=-100 !cannot be merged
        enddo
        if(allocated(vedge1)) deallocate(vedge1)
        if(allocated(n2uz1)) deallocate(n2uz1)
        if(allocated(uz1)) deallocate(uz1)
        if(allocated(o2n1)) deallocate(o2n1)
    endsubroutine
    subroutine node_pg_handle(iflag)
        !把控制点插入模型
        implicit none
        integer,optional::iflag 
        !如果iflag非零,则处理z=-999的情况。此时,找处于地表且处于"输出"状态的节点
        integer::i,j,k,n1,n2,iflag1
        integer,allocatable::vedge1(:) 

        iflag1=0
        if(present(iflag)) iflag1=iflag          
        do i=1,nnode_pg
            if(iflag1==0) then
                call node_pg(i).getnode()            
                do j=1,node_pg(i).nnode
                    vedge1=node2vgedge(node_pg(i).node(j))
                    n1=0
                    if(abs(node_pg(i).z+999.d0)>1.d-6) then                    
                        do k=1,size(vedge1)                        
                            call cut_edge_by_z(vedge1(k),node_pg(i).getz(node_pg(i).node(j)),n1)
                            if(n1>0) then
                                node_pg(i).node(j)=n1
                                exit
                            endif
                        enddo
                    endif
                enddo
            elseif(abs(node_pg(i).z+999.d0)<1.d-6) then
                !此时这些控制点必然是存在的，不必新插入
                call node_pg(i).getnode()            
                do j=1,node_pg(i).nnode
                    vedge1=node2vgedge(node_pg(i).node(j))
                    n1=0
                    n1=maxloc(node(gedge(vedge1).v(2)).z,mask=node(gedge(vedge1).v(2)).isb>0,dim=1)
                    if(n1>0) node_pg(i).node(j)=gedge(vedge1(n1)).v(2)
                enddo              
                
            endif

            
        enddo
        if(allocated(vedge1)) deallocate(vedge1)   
    endsubroutine    
    function cal_selt_volume(this) result(vol)
        implicit none
        class(subelt_tydef)::this
        real(8)::vol
        integer::i,j,k,n1,n2,n3,nid1(ngnode),morder1,order1(100), &
            e1,v2,e2f1(2,ngedge),if1,sf1,iso1(this.nface)
        integer,allocatable::fnode1(:,:)
        real(8)::coord1(3,100)
        
        nid1(this.node(:this.nnum))=[1:this.nnum]
        coord1(1,1:this.nnum)=node(this.node(:this.nnum)).x
        coord1(2,1:this.nnum)=node(this.node(:this.nnum)).y
        coord1(3,1:this.nnum)=node(this.node(:this.nnum)).z
        order1(1:this.nface)=face(this.face(1:this.nface)).nnum
        morder1=maxval(order1(1:this.nface))
        if(allocated(fnode1)) deallocate(fnode1)
        allocate(fnode1(this.nface,morder1))
        sf1=0;
        iso1=0;
        do j=1,this.nface
            n1=this.face(j)
            if(face(n1).isdel) then
                iso1(j)=-1
                cycle
            endif
            if(sf1==0) then
                sf1=j;iso1(sf1)=1
            endif
            fnode1(j,1:face(n1).nnum)=nid1(face(n1).node(1:face(n1).nnum))
            do i=1,face(n1).nnum
                n2=abs(face(n1).edge(i))
                if(e2f1(1,n2)>=-this.nface) then
                    e2f1(2,n2)=sign(j,face(n1).edge(i))
                else
                    e2f1(1,n2)=sign(j,face(n1).edge(i))
                endif
            enddo            
        enddo
        !get consistent nodal order
        do while(any(iso1==0))
            j=minloc(abs(iso1-1),mask=iso1==1,dim=1)
            n1=this.face(j)
            do i=1,face(n1).nnum
                n2=abs(face(n1).edge(i))
                if1=abs(e2f1(1,n2))
                if(iso1(if1)>0) then
                    if1=abs(e2f1(2,n2))
                    if(iso1(if1)>0) then
                        if(e2f1(1,n2)*e2f1(2,n2)>0) then
                            error stop 'error in sub=cal_selt_volume'
                        else
                            cycle
                        endif
                    endif
                endif
                iso1(if1)=1
                if(e2f1(1,n2)*e2f1(2,n2)>0) then
                    n3=this.face(if1)
                    fnode1(if1,1:face(n3).nnum)=nid1(face(n3).node(face(n3).nnum:1:-1))
                    do k=1,face(n3).nnum
                        e1=abs(face(n3).edge(k))
                        if(abs(e2f1(1,e1))==if1) then
                            e2f1(1,e1)=sign(if1,-face(n3).edge(k))
                        else
                            e2f1(2,e1)=sign(if1,-face(n3).edge(k))
                        endif
                    enddo
                endif
            enddo 
            iso1(j)=2
        enddo
        vol=polyhedron_volume_3d (coord1(:,1:this.nnum) ,morder1 , &
            this.nface, fnode1, this.nnum, order1 )        
    endfunction
    
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
       
        n1=nelt*(soillayer+1)+(abs(iedge)-1)*soillayer+1
        n2=0
        do i=n1,n1+soillayer-1            
            call face(i).getface(i,vface,nvface1)                
        enddo
        
        vface=vface(:nvface1)        
        
    end function    
    
    function edge2hgedge(iedge) result(hedge)
    !Given: iedge, the edge id in the mother 2d model
    !Return: the first vertical face id in the face,extruded by the edge
        implicit none
        integer,intent(in)::iedge
        integer,allocatable::hedge(:)        
        integer::i,n1,nhedge1=0        

        do i=0,soillayer
            n1=iedge+i*nedge
            call gedge(n1).getedge(n1,hedge,nhedge1)                
        enddo
        
        hedge=hedge(:nhedge1)        
        
    end function        
    
    subroutine Gen_SubElement()
       
	    implicit none
		!logical,intent(in)::istet
	    integer::i,j,k,n1,n2,n3,n4=0,iat(4)=0,e1(2),N5,v1(2)
	    real(8)::t1,t2,ver1(3,3)
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
        !edge的顺序，按层生成，上下相邻横线段号相差nedge        
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
                gedge(n1).iscedge=edge(j).iscedge
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
                !cal the normal of the prototype face.scale back.
                ver1(1,1:3)=node(face(n1).node(1:3)).x
                ver1(2,1:3)=node(face(n1).node(1:3)).y
                ver1(3,1:3)=node(face(n1).node(1:3)).z 
                face(n1).normal=NORMAL_TRIFACE(ver1,.true.)
            enddo
        enddo
        !vertical face
        do i=1,nedge
            if(edge(i).v(1)*edge(i).v(2)/=0) then
                ver1(1,1)=node(edge(i).v(2)).x-node(edge(i).v(1)).x
                ver1(2,1)=node(edge(i).v(2)).y-node(edge(i).v(1)).y
                ver1(:,2)=0.d0
                if(abs(ver1(1,1))<1.d-7) then
                    if(ver1(2,1)>0) then
                        ver1(1,2)=-1.d0
                    else
                        ver1(1,2)=1.d0
                    endif
                else
                    t1=atan(ver1(2,1)/ver1(1,1))+Pi/2.d0
                    ver1(1,2)=cos(t1);ver1(2,2)=sin(t1)
                endif
            endif
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
                n2=4;
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
                
                face(n1).normal=ver1(:,2)

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
                selt(nselt).vol=abs(ELT(I).PROPERTY(1)*SUM(node(selt(nselt).node(4:6)).z-node(selt(nselt).node(1:3)).z)/3.D0)
                !t1=selt(nselt).calvol()
                !if(abs(selt(nselt).vol-t1)>1.d-7) then
                !    print *, 'v1,v2=',selt(nselt).vol,t1
                !endif
                
                !if(selt(nselt).vol<1.d-7) write(*,*) 'the volume of element i is less than 1.d-7. i=', nselt,selt(nselt).vol                

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
        
        call set_cow_face_edge()
        
        call vseg_pg_handle()
        
        call vface_pg_handle()

        call node_pg_handle()
        !find faces/edges on boundary
        face.isdel=.true.
        face.isb=0;gedge.isb=0;node.isb=0
        do i=1,nselt
            if(selt(i).isdel)cycle
            !if(selt(i).vol<1.d-7) write(*,*) 'the volume of element i is less than 1.d-7. i=', i,selt(i).vol 
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
                WHERE(GEDGE(XZONE(j).CUTFACE).ISB==0) 
                    GEDGE(XZONE(j).CUTFACE).ISB=4
                    GEDGE(XZONE(j).CUTFACE).ISCEDGE=-J !MARK IT AND CANNOT BE MERGED
                ENDWHERE
            ENDIF                
        ENDDO
        
        
        
        DO i=1,NCOW
            !MARK THE FACES ON THE WALL
            WHERE(FACE(COWALL(I).FACE(1).NODE).ISB==0) FACE(COWALL(I).FACE(1).NODE).ISB=5
            GEDGE(COWALL(I).FACE(4).NODE).ISCEDGE=-I*100 !MARK IT AND CANNOT BE MERGED
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
        
       
        
        
        !call mergeface()
        
        DO I=1,ZNUM
            IF(.NOT.ALLOCATED(ZONE(I).HFACE)) ALLOCATE(ZONE(I).HFACE(SOILLAYER))
            IF(.NOT.ALLOCATED(ZONE(I).VFACE)) ALLOCATE(ZONE(I).VFACE(SOILLAYER))
            IF(.NOT.ALLOCATED(ZONE(I).BFACE)) CYCLE
            DO J=1,SOILLAYER
                
                CALL VOLUME_NUMBER_INSIDE_ZONE(I,J)
                
                !outer HORIZONTAL faces
                N3=ZONE(I).BFACE(J).NNUM
                ZONE(I).HFACE(J).NODE=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==0 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ILAYER/=0)
                ZONE(I).HFACE(J).NNUM=SIZE(ZONE(I).HFACE(J).NODE)
                
                !outer vertical faces
                ZONE(I).VFACE(J).NODE=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB==1 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==1)
                ZONE(I).VFACE(J).NNUM=SIZE(ZONE(I).VFACE(J).NODE)
            ENDDO
        ENDDO

        !vertical lines
        do i=1,nvseg_pg            
            do j=1,size(vseg_pg(i).edge)
                if(gedge(vseg_pg(i).edge(j)).isb==0)then
                    gedge(vseg_pg(i).edge(j)).isb=6
                    node(gedge(vseg_pg(i).edge(j)).v).isb=77777
                    !find containing volume
                    iw1=gedge2selt(vseg_pg(i).edge(j))
                    if(size(iw1)>0) then
                        if(all(selt(iw1).igmshvol-selt(iw1(1)).igmshvol==0)) then
                            vseg_pg(i).igmshvol(J)=selt(iw1(1)).igmshvol
                        else
                            print *, 'Unexpected error. All elts should be in the same volume. loc=nvseg_pg.'
                        endif
                        
                    endif                   
                else
                    vseg_pg(i).igmshvol(J)=-1 !no need to embed in a volume.
                endif
            enddo
        enddo 

        call node_pg_handle(iflag=1)

        do i=1,nnode_pg
            node(node_pg(i).node).isb=77777 !cann't be merged
        enddo       
        
        
        
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
    
    
    subroutine merge_face()
        implicit none
        integer::i,j,k,n1,n2,f1(2),f2(2),e1(2),mid1=0,if1=0,v1(3)
        integer,allocatable::edge1(:)
        integer::stack1(1000),nst1=0,nedge1,iw1(ngnode)
        integer,allocatable::iw2(:)
        real(8)::t1
        do i=1,ngedge
            if(gedge(i).isdel.or.gedge(i).isb<1) cycle
            n1=gedge(i).nface
            if(n1<1) cycle
            n2=count(face(gedge(i).face(:n1)).isb>0)
            if(n2==2) then                
                f1=pack([1:n1],face(gedge(i).face(:n1)).isb>0)                
                f2=gedge(i).face(f1)
                gedge(i).mface=f2
                if(gedge(i).iscedge/=0) cycle !the constrained edges cannot be merged.
                t1=dot_product(face(f2(1)).normal,face(f2(2)).normal)
                if(abs(abs(t1)-1.0d0)<1.d-10) then
                    gedge(i).isb=-1 !mark the edges can be  merged
                    !face(f2).isb=-1                    
                    !node(gedge(i).v).isb=-1
                endif

            endif
            
        enddo
        
        NGPF=0;iw1=0
        do while(any(gedge.isb==-1))
            NGPF=NGPF+1
            IF(NGPF>SIZE(GPFACET)) CALL GPFACET_ENLARGE_AR(GPFACET,10)  
            i=minloc([1:ngedge],mask=gedge.isb==-1,dim=1)
            stack1(1:2)=gedge(i).mface
            gedge(i).isb=-2 !is handled,no need to check again.
            nst1=2
            GPFACET(NGPF).ONEFACE=stack1(1)
            !if(allocated(edge1)) deallocate(edge1)
            nedge1=0
            
            do while(nst1>0)
                if1=stack1(nst1)
                nst1=nst1-1
                if(face(if1).isb==-NGPF) CYCLE
                face(if1).isb=-NGPF
                do j=1,face(if1).nnum
                    iw1(face(if1).node(j))=ngpf !mark facet node
                    
                    n1=abs(face(if1).edge(j))
                    if(gedge(n1).isb==-2) cycle !is handled
                    if(gedge(n1).isb==-1) then
                        !face(gedge(n1).mface).isb=-NGPF
                        gedge(n1).isb=-2 !marked
                        nst1=nst1+1
                        if(gedge(n1).mface(1)==if1) then                            
                            stack1(nst1)=gedge(n1).mface(2)
                        else
                            stack1(nst1)=gedge(n1).mface(1)
                        endif
                    else
                        nedge1=nedge1+1
                        if(nedge1>size(edge1)) call enlarge_ar(edge1,50)
                        edge1(nedge1)=n1
                    endif
                
                enddo
            
            enddo
            GPFACET(NGPF).EDGE=EDGE1(1:nedge1)            
            GPFACET(NGPF).NEDGE=nedge1
            GPFACET(NGPF).ID=NGPF
            !find the floating points on the facet
            iw1(GEDGE(EDGE1(:nedge1)).V(1))=0
            iw1(GEDGE(EDGE1(:nedge1)).V(2))=0
            iw2=pack([1:ngnode],iw1==ngpf)
            do j=1,size(iw2)
                n1=iw2(j)
                !如果这点连接着一个输出边或为node_pg上的点，这个点必须输出
                if(any(gedge(gadjlist(n1).edge(1:gadjlist(n1).count)).isb>0)) then
                    node(n1).isb=77777
                endif               
            enddo
            GPFACET(NGPF).FLOATP=pack(iw2,node(iw2).isb==77777) 
            !GPFACET(NGPF).FLOATP=pack([1:ngnode],iw1==ngpf.and.node.isb==77777) 
            GPFACET(NGPF).NFP=SIZE(GPFACET(NGPF).FLOATP)
            IW1(GPFACET(NGPF).FLOATP)=0
            WHERE(IW1==NGPF) NODE(1:ngnode).ISB=-1 !MERGED POINTS
            
            IF(ISMERGED>1) CALL GPFACET(NGPF).SET_FACET()
        enddo
        !ISMERGED==1,MERGED FACES AND EDGES BOTH;==2,MERGED FACES ONLY;==0, NOT BOTH.
        IF(ISMERGED==1) THEN
        ! merge edges
            do i=1,ngnode
                if(node(i).isb<1.or.node(i).isb==77777) cycle
                n1=gadjlist(i).count
                n2=count(gedge(gadjlist(i).edge(1:n1)).isb>0)
                if(n2==2) then
                    e1=pack([1:n1],gedge(gadjlist(i).edge(1:n1)).isb>0)
                    v1=[gadjlist(i).node(e1(1)),i,gadjlist(i).node(e1(2))]
                    if(iscoline([node(v1(1)).x,node(v1(1)).y,node(v1(1)).z],&
                        [node(v1(2)).x,node(v1(2)).y,node(v1(2)).z],&
                        [node(v1(3)).x,node(v1(3)).y,node(v1(3)).z])) then
                        
                        node(i).isb=-1
                        j=gadjlist(i).edge(e1(1))
                        if(gedge(j).v(1)==i) then
                            gedge(j).v(1)=v1(3)
                        else
                            gedge(j).v(2)=v1(3)
                        endif

                        gedge(gadjlist(i).edge(e1(2))).isb=-1

                        !update adjlist
                        do k=1,3,2
                            n1=gadjlist(v1(k)).count
                            n2=K+2
                            if(n2>3) n2=1
                            where(gadjlist(v1(k)).node(1:n1)==i)
                                gadjlist(v1(k)).node=v1(n2)
                                gadjlist(v1(k)).edge=J
                            ENDWHERE
                        enddo

                        gadjlist(i).edge(e1)=-1;gadjlist(i).node(e1)=-1;gadjlist(i).count=gadjlist(i).count-2

                    endif

                endif
            enddo
            
            !update face
            do i=1,nface
                if(face(i).isb<1) cycle
                n1=face(i).nnum
                face(i).edge=pack(face(i).edge(1:n1),gedge(abs(face(i).edge(1:n1))).isb>0)
                face(i).node=pack(face(i).node(1:n1),node(face(i).node(1:n1)).isb>0)
                face(i).nnum=size(face(i).edge)
            enddo
            !update gpf
            do i=1,ngpf
                n1=GPFACET(i).NEDGE
                GPFACET(i).EDGE=pack(GPFACET(i).EDGE(1:n1),gedge(GPFACET(i).EDGE(1:n1)).isb>0)
                GPFACET(i).NEDGE=size(GPFACET(i).EDGE)
                GPFACET(I).FLOATP=pack(GPFACET(I).FLOATP,node(GPFACET(I).FLOATP).isb>0) 
                GPFACET(i).NFP=SIZE(GPFACET(I).FLOATP)
                CALL GPFACET(i).SET_FACET()
            enddo
        ENDIF
    endsubroutine
    
    logical function iscoline(v1,v2,v3)
        implicit none 
        real(8),intent(in)::v1(3),v2(3),v3(3)
        real(8)::v4(3),v5(3),t1,t2,tol1

        iscoline=.false.
        !tol1=100.d0*epsilon(tol1)
        tol1=1.0d-10
        v4=v2-v1;v5=v3-v1
        t1=norm2(v4);t2=norm2(v5)
        if(abs(t1)<tol1.or.abs(t2)<tol1) then
            iscoline=.true.
        else
            v4=v4/t1;v5=v5/t2
            t1=abs(dot_product(v4,v5))
            if(abs(t1-1.d0)<tol1) iscoline=.true.
        endif

    endfunction

    SUBROUTINE poly_set_facet(THIS)
        implicit none
        class(general_poly_face_tydef)::this
		integer::i,j,k,n1,n2,n3,e1,sign1,norder1(ngnode),g2ln(ngnode),ia1(2),stat_err,NLOOP1
        integer,allocatable::NODE1(:),edge1(:),l2gn(:)
        real(8),allocatable::coord1(:,:)
		real(8)::t1,t2,xi,yi,xj,yj,x1,y1,HOLE1(3,100)        
		logical::tof1
        
        call FACET_INNER_EDGES(this)
        
        norder1=0        
        do j=1,this.nedge              
            norder1(gedge(this.edge(j)).v)=norder1(gedge(this.edge(j)).v)+1            
        enddo
        !find loop
        n1=maxloc(abs(face(this.oneface).normal),dim=1)
        n1=mod(n1,3)+1
        n2=mod(n1,3)+1
        if(n1>n2) then
            n3=n1;n1=n2;n2=n3
        endif
        ia1=[n1,n2]
        !n3=count(norder1>0)
        l2gn=pack([1:ngnode],norder1>0)        
        n3=size(l2gn)
        do i=1,n3
            g2ln(l2gn(i))=i
        enddo
        allocate(coord1(2,n3),node1(n3+1),edge1(this.nedge))
        
        do i=1,2            
            select case(ia1(i))
            case(1)
                coord1(i,:)=node(l2gn).x
            case(2)
                coord1(i,:)=node(l2gn).y
            case(3)
                coord1(i,:)=node(l2gn).z
            endselect
        enddo

      
        
        GEDGE(THIS.EDGE).IT1=THIS.ID
        nloop1=0
        do while(any(GEDGE(this.edge).IT1==THIS.ID))
            nloop1=nloop1+1
            node.subbw=0!借用
            node.layer=0
            do j=1,this.nedge
                if(GEDGE(this.edge(j)).IT1/=THIS.ID) cycle
                do k=1,2
                    n2=gedge(this.edge(j)).v(k)
                    if(k==1) then
                        sign1=1 !first vetex
                    else
                        sign1=-1 !second vetex
                    endif
                    if(node(n2).subbw==0) then
                        node(n2).subbw=this.edge(j)*sign1
                    else
                        node(n2).layer=this.edge(j)*sign1
                    endif
                enddo
            enddo
            
            n1=minloc(coord1(1,:),mask=norder1(l2gn)>0,dim=1)            
            NODE1(1)=L2GN(N1)
            E1=OUTER_EDGE(N1,PI/2.0D0)
            N1=2
            IF(E1<0) N1=1
            NODE1(2)=GEDGE(ABS(E1)).V(N1)
            J=2
            edge1(1)=e1
            N3=ABS(E1)
            NORDER1(NODE1(1:2))=NORDER1(NODE1(1:2))-1
            GEDGE(ABS(E1)).IT1=0
            DO WHILE(NODE1(J)/=NODE1(1))                
                N2=NODE1(J)
                if(norder1(n2)==1) then
                    if(abs(node(n2).subbw)==n3) then
                        n3=node(n2).layer
                    else
                        n3=node(n2).subbw
                    endif
                elseif(norder1(n2)>1) then
                    T1=ATAN2(COORD1(2,G2LN(NODE1(J-1)))-COORD1(2,G2LN(N2)),&
                        COORD1(1,G2LN(NODE1(J-1)))-COORD1(1,G2LN(N2)))
                    if(nloop1==1) then
                        N3=OUTER_EDGE(G2LN(NODE1(J)),T1) !outer loop 
                    else
                        N3=OUTER_EDGE(G2LN(NODE1(J)),T1,.FALSE.)    !inner loop                    
                    endif
                else
                    error stop 'nodal order should be >1. sub=INNERLOOP'
                endif
            
                j=j+1    
                if(n3>0) then
                    node1(j)=gedge(n3).v(2)
                    edge1(j-1)=n3 
                else
                    n3=abs(n3)
                    node1(j)=gedge(n3).v(1)
                    edge1(j-1)=-n3
                endif
                NORDER1(GEDGE(N3).V)=NORDER1(GEDGE(N3).V)-1
                GEDGE(N3).IT1=0
                
            ENDDO
        
            THIS.NPOLY=THIS.NPOLY+1
            IF(THIS.NPOLY>SIZE(THIS.POLY)) CALL ENLARGE_AR(THIS.POLY,2)
            THIS.POLY(THIS.NPOLY).NNUM=J-1
            THIS.POLY(THIS.NPOLY).NODE=NODE1(1:J-1)
            THIS.POLY(THIS.NPOLY).EDGE=EDGE1(1:J-1)
            IF(NLOOP1>1) THEN
                THIS.NHOLE=THIS.NHOLE+1
                IF(THIS.POLY(THIS.NPOLY).NNUM<5) THEN
                    HOLE1(1,THIS.NHOLE)=SUM(NODE(NODE1(1:J-1)).X)/(J-1)
                    HOLE1(2,THIS.NHOLE)=SUM(NODE(NODE1(1:J-1)).Y)/(J-1)
                    HOLE1(3,THIS.NHOLE)=SUM(NODE(NODE1(1:J-1)).Z)/(J-1)
                ELSE
                    DO I=1,J-1                        
                        N2=MOD(I,J-1)+1
                        N3=MOD(N2,J-1)+1
                        N1=G2LN(NODE1(I));N2=G2LN(NODE1(N2));N3=G2LN(NODE1(N3))
                        !ASSUME ALL NODES OF THE LOOPS ARE ORDERED IN CCW
                        K=triangle_orientation_2d(coord1(:,[N1,N2,N3]))
                        IF(K==0) THEN
                            HOLE1(1,THIS.NHOLE)=SUM(NODE(L2GN([N1,N2,N3])).X)/3.D0
                            HOLE1(2,THIS.NHOLE)=SUM(NODE(L2GN([N1,N2,N3])).Y)/3.D0
                            HOLE1(3,THIS.NHOLE)=SUM(NODE(L2GN([N1,N2,N3])).Z)/3.D0
                        ENDIF
                        
                    ENDDO
                ENDIF
            ENDIF
            
        ENDDO
        THIS.HOLE=HOLE1(:,1:THIS.NHOLE) 
        
        deallocate(coord1,node1,edge1,l2gn,stat=stat_err)
             
        
    CONTAINS
        INTEGER FUNCTION OUTER_EDGE(IN1,REF_PHI,ISOUT)
            IMPLICIT NONE
            INTEGER,INTENT(IN)::IN1
            REAL(8),INTENT(IN)::REF_PHI !已知上一边界(参考边)的atan2的值
            LOGICAL,INTENT(IN),OPTIONAL::ISOUT !是否找外边界
            INTEGER::I,N1,N2,N3
            REAL(8)::DX1,DY1,PHI1,PHI2,PN1
            LOGICAL::ISOUT1
            
            ISOUT1=.TRUE.
            IF(PRESENT(ISOUT)) ISOUT1=ISOUT
            
            PN1=PI-REF_PHI !ROTATE THE REFERENCE EDGE TO -X-AXIS(PI) 
            IF(ISOUT1) THEN
                PHI2=PI
            ELSE
                PHI2=-PI
            ENDIF
            
            do i=1,gadjlist(l2gn(IN1)).count
                n2=gadjlist(l2gn(IN1)).edge(i)
                IF(GEDGE(N2).IT1==THIS.ID) THEN
                    N3=GADJLIST(l2gn(IN1)).NODE(i)
                    DX1=coord1(1,g2ln(N3))-COORD1(1,IN1)
                    DY1=coord1(2,g2ln(N3))-COORD1(2,IN1)
                    PHI1=ATAN2(DY1,DX1)+PN1 
                    !LET THE ANGLE IN -PI<=PHI1<=PI
                    IF(PHI1>PI) THEN
                        PHI1=PHI1-2*PI
                    ELSEIF(PHI1<-PI) THEN
                        PHI1=PHI1+2*PI
                    ENDIF
                    IF(ISOUT1) THEN
                        IF(PHI1<PHI2) THEN
                            PHI2=PHI1  
                            OUTER_EDGE=N2
                            IF(GEDGE(N2).V(1)==N3) THEN
                                OUTER_EDGE=-N2
                            ENDIF
                        ENDIF
                    ELSE
                        IF(ABS(PHI1-PI)>1.D-7.AND.PHI1>PHI2) THEN
                            PHI2=PHI1  
                            OUTER_EDGE=N2
                            IF(GEDGE(N2).V(1)==N3) THEN
                                OUTER_EDGE=-N2
                            ENDIF
                        ENDIF                        
                    ENDIF
                ENDIF            
            ENDDO       
        
        ENDFUNCTION
    ENDSUBROUTINE    
    
    SUBROUTINE FACET_INNER_EDGES(THIS)
        implicit none
        class(general_poly_face_tydef)::this 
        INTEGER::I,J,IW1(NGEDGE)
        
        IW1=0
        !当某线段为不可合并的非环线段(不成为某个闭环的线段)时，其编号将在this.edge中会重复出现2次
        IW1(THIS.EDGE)=IW1(THIS.EDGE)+1
        
        DO I=1,THIS.NEDGE
            !IF(ANY(GEDGE(THIS.EDGE(I)).MFACE<1)) CYCLE
            !IF(FACE(GEDGE(THIS.EDGE(I)).MFACE(1)).ISB==FACE(GEDGE(THIS.EDGE(I)).MFACE(2)).ISB) THEN
            IF(IW1(THIS.EDGE(I))==2) THEN
                THIS.NPOLY=THIS.NPOLY+1
                IF(THIS.NPOLY>SIZE(THIS.POLY)) CALL ENLARGE_AR(THIS.POLY,2)
                THIS.POLY(THIS.NPOLY).NNUM=2
                THIS.POLY(THIS.NPOLY).NODE=GEDGE(THIS.EDGE(I)).V
                THIS.POLY(THIS.NPOLY).EDGE=[THIS.EDGE(I)]
                IW1(THIS.EDGE(I))=-2
            ENDIF 
            
        ENDDO
        
        THIS.EDGE=PACK(THIS.EDGE,IW1(THIS.EDGE)==1)
        THIS.NEDGE=SIZE(THIS.EDGE)
                
    ENDSUBROUTINE
    
    function gedge2selt(iedge) result(seltf)
        implicit none
        integer,intent(in)::iedge
        integer,allocatable::seltf(:)
        integer::i,j,elt2(2),elt1(100),n1
        
        n1=0
        do i=1,gedge(iedge).nface
            elt2=face(gedge(iedge).face(i)).e
            do j=1,2
                if(elt2(j)>0.and.all(elt1(:n1)-elt2(j)/=0)) then
                    n1=n1+1
                    elt1(n1)=elt2(j)
                endif
            enddo            
        enddo        
        seltf=elt1(1:n1)
        
    endfunction
    
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
    
    !face and edge needed to define the cowall
    subroutine set_cow_face_edge()
        implicit none
        integer::i,j,k,staterr,n1,n2,n3,n4
        real(8)::t1,t2,pt1(3,2),be1(2),te1(2),ite1,ibe1,dx1,ia1(2)
        integer,allocatable::eiz1(:)
        
        do i=1,ncow
            !faces on the wall
            n1=size(cowall(i).node) 
            n2=size(cowall(i).edge)
            gedge.isb=0
            face.isb=0

            do j=1,n2
                n3=j+1
                if(n3>n1) n3=1
                pt1(1,:)=node(cowall(i).node([j,n3])).x
                pt1(2,:)=node(cowall(i).node([j,n3])).y
                !pt1(3,:)=node(node1([j,n3])).z
                ia1=[j,n3]
                do k=1,2                    
                    te1(k)=node(cowall(i).node(ia1(k))).we(1)
                    be1(k)=node(cowall(i).node(ia1(k))).we(2)
                enddo
                eiz1=edge2vface(abs(cowall(i).edge(j)))
                do k=1,size(eiz1)
                    n4=eiz1(k)
                    if(face(n4).isdel) cycle
                    t1=sum(node(face(n4).node(1:face(n4).nnum)).z)/face(n4).nnum
                    
                    dx1=pt1(1,2)- pt1(1,1)
                    if(abs(dx1)>1.d-7) then                        
                        t2=sum(node(face(n4).node(1:face(n4).nnum)).x)/face(n4).nnum
                        t2=(t2- pt1(1,1))/dx1
                    else
                        t2=sum(node(face(n4).node(1:face(n4).nnum)).y)/face(n4).nnum
                        t2=(t2- pt1(2,1))/( pt1(2,2)- pt1(2,1))
                    endif                    
                    
                    ite1=te1(1)+t2*(te1(2)-te1(1))
                    ibe1=be1(1)+t2*(be1(2)-be1(1))
                    if(ite1>t1.and.ibe1<t1) then
                        face(n4).isb=2
                    else
                        face(n4).isb=1
                    endif
                enddo
            enddo 
            cowall(i).face(1).node=pack([1:nface],face(1:nface).isb==2.or.face(1:nface).isb==1) !all faces streteched by the wall line
            cowall(i).face(1).nnum=size(cowall(i).face(1).node)
            cowall(i).face(2).node=pack([1:nface],face(1:nface).isb==2) !just the faces on the real wall
            cowall(i).face(2).nnum=size(cowall(i).face(2).node)            
            !get boundary edges
            gedge.isb=0
            do j=1,cowall(i).face(1).nnum
                n3=cowall(i).face(1).node(j)
                gedge(abs(face(n3).edge(1:face(n3).nnum))).isb=gedge(abs(face(n3).edge(1:face(n3).nnum))).isb+1
            enddo
            cowall(i).face(3).node=pack([1:ngedge],gedge(1:ngedge).isb==1)
            cowall(i).face(3).nnum=size(cowall(i).face(3).node)
            gedge.isb=0
            do j=1,cowall(i).face(2).nnum
                n3=cowall(i).face(2).node(j)
                gedge(abs(face(n3).edge(1:face(n3).nnum))).isb=gedge(abs(face(n3).edge(1:face(n3).nnum))).isb+1
            enddo
            cowall(i).face(4).node=pack([1:ngedge],gedge(1:ngedge).isb==1)
            cowall(i).face(4).nnum=size(cowall(i).face(4).node)            
        enddo

        deallocate(eiz1,stat=staterr)
        
    endsubroutine
    
    subroutine cut_face_by_cow()
        implicit none
        integer::i,j,k,k1,n1,n2,n3,staterr
        real(8)::pt1(3,2)
        integer,allocatable::eiz1(:),node1(:)
        
        do i=1,ncow
            if(.not.allocated(cowall(i).node)) call cowall(i).setboundary()
            n1=size(cowall(i).node) 
            allocate(node1(n1))
            
            do k1=1,2
                node.bw=0
                gedge.isb=0
                face.isb=0
                node1=0
                do j=1,n1
                    eiz1=node2vgedge(cowall(i).node(j))
                    
                    do k=1,size(eiz1)
                        !if(gedge(k).isdel) cycle
                        call cut_edge_by_z(eiz1(k),node(cowall(i).node(j)).we(k1),node1(j))
                        if(node1(j)>0) exit
                    enddo
                    
                enddo
                !cut bedge by seg in the same vertical face
                n2=size(cowall(i).edge)
                do j=1,n2
                    eiz1=edge2hgedge(abs(cowall(i).edge(j)))
                    n3=j+1
                    if(n3>n1) n3=1
                    pt1(1,:)=node(node1([j,n3])).x
                    pt1(2,:)=node(node1([j,n3])).y
                    pt1(3,:)=node(node1([j,n3])).z
                    
                    do k=1,size(eiz1)
                        !if(gedge(k).isdel) cycle
                        call cut_edge_by_seg_in_same_vertical_face(eiz1(k),pt1)
                    enddo
                    
                    deallocate(eiz1)
                    eiz1=edge2vface(abs(cowall(i).edge(j)))                    
                    do k=1,size(eiz1)
                        !if(face(k).isdel) cycle
                        call cut_face2(eiz1(k))
                    enddo
                    
                enddo      
                !edge around the cut face
                !cowall(i).wedge(k1).node=pack([1:ngedge],gedge(1:ngedge).isb==2)
                !cowall(i).wedge(k1).nnum=size(cowall(i).wedge(k1).node)
                
            enddo            

        enddo
        
        deallocate(eiz1,node1,stat=staterr)
    end subroutine
    
    subroutine xzone_cut_selt()
    
        implicit none
        integer::i,j,k,k1,k2,j1,n1,n2,n3,n4,n5,v1(2),nc1,nvf1,mbw1(2)
        integer::ia2(8),staterr,instate
        real(8)::t1,t2,tol1=1.d-8,xi1(3),p1(3),p2(3),p3(3)
        integer,allocatable::viz1(:),eiz1(:)
        
        call cut_face_by_cow()

        
        do i=1,nxzone
            !InsPt1=0 !InsPt1(i)=0,need to check;-1,no intercecting point;>0,intercepting, the intersecting point is node(gedge(i).cutpoint)
            !cutedge1=0  !InsPt1(i)=0,need to check;any is -1,no intercecting point;any >0,intercepting, the intersecting point is node(gedge(i)).cutpoint              
            !t1=xzone(i).te
            !if(allocated(iw1)) deallocate(iw1)
            !allocate(iw1(ngnode))
            !借用bw
            node.bw=0
            gedge.isb=0
            face.isb=0
            selt.isb=0
            node.s=0
            !借用s
            do j=1,ngnode
                node(j).s=node(j).z-xzone(i).getz([node(j).x,node(j).y])
            enddo
            
            !only  insert zone boundary line into the gmodel
            if(xzone(i).isx==2) then
                do j=1,zone(xzone(i).izone).nbn
                    n1=zone(xzone(i).izone).bnode(j)
                    eiz1=node2vgedge(n1)
                    t1=xzone(i).getz([node(n1).x,node(n1).y])
                    do k=1,size(eiz1)
                        !if(gedge(k).isdel) cycle
                        call cut_edge_by_z(eiz1(k),t1)
                    enddo
                
                enddo
                
                do j=1,zone(xzone(i).izone).nbe
                    eiz1=edge2hgedge(zone(xzone(i).izone).bedge(j))
                    
                    do k=1,size(eiz1)
                        call edge_Intersect_cutface(i,eiz1(k))                        
                    enddo
                
                enddo                
                
                do j=1,zone(xzone(i).izone).nbe
                    eiz1=edge2vface(zone(xzone(i).izone).bedge(j))
                    
                    do k=1,size(eiz1)
                        !if(face(k).isdel) cycle
                        call cut_face2(eiz1(k))
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
                if(maxval(node(selt(j).node(1:n1)).s)<=-tol1) cycle                
                if(xzone(i).isx==1.and.minval(node(selt(j).node(1:n1)).s)>-tol1) then
                    selt(j).isdel=.true.
                    call remove_selt(j)
                    cycle
                endif
                
                selt(j).isb=1
            enddo
            
            viz1=pack([1:nselt],selt(1:nselt).isb>0)             
            do j=1,size(viz1)
                !mark zone node 
                n1=viz1(j)
                do j1=1,selt(n1).nnum
                    n2=selt(n1).node(j1)
                    t2=node(n2).s
                    if(abs(t2)<tol1) then
                        node(n2).bw=2 !on cut face
                    elseif(t2>0.d0) then
                        node(n2).bw=1 !over cut face
                    else
                        node(n2).bw=-1 !below cut face
                    endif                        
                enddo
            enddo
            
            
            do j1=1,size(viz1)
                !mark zone edge 
                j=viz1(j1)
                
                do k=1,selt(j).nface
                    n1=selt(j).face(k)
                    if(face(n1).isb/=0) cycle
                    face(n1).isb=-1 !no need to check again
                    n2=face(n1).nnum
                    n4=0;mbw1=0
                    do k1=1,n2
                        n3=abs(face(n1).edge(k1))
                        if(gedge(n3).isb/=0) then
                            if(gedge(n3).isb==1) then
                                face(n1).isb=1
                            elseif(gedge(n3).isb==2) then
                                n4=n4+1
                            endif
                        else
                            v1=node(gedge(n3).v).bw                            
                            if(all(v1==2))then
                                gedge(n3).isb=2 !on the cut face
                                n4=n4+1
                            else
                                
                                if(mbw1(1)==0.and.any(v1==1)) mbw1(1)=1 !cut plane just through the nodes, not cutting the edge.
                                if(mbw1(2)==0.and.any(v1==-1)) mbw1(2)=-1
                                
                                if((v1(1)*v1(2))<0)then
                                    gedge(n3).isb=1 !cross the cut face
                                    face(n1).isb=1 !cross the cut face
                                else
                                    gedge(n3).isb=-1 !no need to check again
                                endif
                            endif
                        endif
                    enddo
                    if(n4==n2) face(n1).isb=2 !on the cut face
                    
                    !the cut plane is just through the nodes, not cutting any edge.
                    if(face(n1).isb==-1.and.mbw1(1)*mbw1(2)<0) then
                         face(n1).isb=1 !cross the cut face
                    endif
                    

                enddo
            enddo
            

            
            viz1=pack([1:ngedge],gedge(:ngedge).isb==1)
            !cut edge
            do j=1,size(viz1)
                n1=viz1(j)
                !call cut_edge_by_z(n1,t1)
                call edge_Intersect_cutface(i,n1)
            enddo
            !cut face
            viz1=pack([1:nface],face(1:nface).isb==1) 
            
            do j=1,size(viz1)
                n1=viz1(j)
                call cut_face2(n1)
            enddo            
            
            
            viz1=pack([1:nselt],selt(1:nselt).isb>0) 
             do j=1,size(viz1)
                n1=viz1(j)
                call cut_selt2(N1,xzone(i).isx)
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
        node.s=0
        deallocate(viz1,eiz1,stat=staterr)
    
    endsubroutine    
    
    subroutine edge_Intersect_cutface(ixz,iedge)
        implicit none
        integer,intent(in)::ixz,iedge
        integer::n1,n2,instate
        real(8)::p1(3),p2(3),p3(3),t1
        
        if(gedge(iedge).isdel) return
        
        if(xzone(ixz).iplane==0) then
            t1=xzone(ixz).getz([0.d0,0.0d0])
            call cut_edge_by_z(iedge,t1)
            
        elseif(xzone(ixz).iplane==1) then    
            n2=gedge(iedge).v(1)
            p1=[node(n2).x,node(n2).y,node(n2).z]
            n2=gedge(iedge).v(2)
            p2=[node(n2).x,node(n2).y,node(n2).z]                            
            call xzone(ixz).edgexface(p1,p2,instate,p3,t1)
            if(instate==1) then
                
                !call cut_edge_by_z(iedge,p3(3))
                 if(abs(t1)<1.d-8) then
                    node(gedge(iedge).v(1)).bw=2
                    node(gedge(iedge).v(1)).z=p3(3)
                elseif(abs(1.0d0-t1)<1.d-8) then
                    node(gedge(iedge).v(2)).bw=2
                    node(gedge(iedge).v(2)).z=p3(3)
                else
                    ngnode=ngnode+1
                    if(size(node)<ngnode) call enlarge_ar(node,100)
                    node(ngnode).x=p3(1);node(ngnode).y=p3(2);node(ngnode).z=p3(3);
                    node(ngnode).bw=2 !mark
                    node(ngnode).s=0.d0
                    node(ngnode).iptr=ngnode
                    gedge(iedge).cutpoint=ngnode
                    call cut_edge(iedge,ngnode) 
                endif                
            elseif(instate==2) then                
                node(gedge(iedge).v).bw=2
                node(gedge(iedge).v).z=p3(3) 
                gedge(iedge).isb=2
            endif
        else
            print *, 'No such cut plane  type=', xzone(ixz).iplane
            error stop
        endif
    end subroutine
    
    
    subroutine cut_edge_by_z(iedge,z,inode)
       
        implicit none
        integer,intent(in)::iedge
        integer,intent(out),optional::inode
        real(8),intent(in)::z
        real(8)::t2,xi1(3)
        integer::v1(2),inode1
        
        inode1=0
        if(present(inode)) inode=inode1
        if(gedge(iedge).isdel) return
        v1=gedge(iedge).v
        
        !平行xy平面的edge
        !下面平行判断只使用xy平面
        if(abs(node(v1(1)).z-node(v1(2)).z)<1.d-8) then
            if(abs(node(v1(1)).z-z)<1.d-8) then
                node(v1).bw=2
                inode1=v1(1)
                node(v1).z=z
                !if(gedge(iedge).cutpoint>0) node(gedge(iedge).cutpoint).bw=2
            endif
            if(present(inode)) inode=inode1
            return
        endif
        
        t2=(z-node(v1(1)).z)/(node(v1(2)).z-node(v1(1)).z)
        
        
       
        !if(t2<0.d0.or.t2>1.d0) return !no cross
        !cross at end nodes
        if(abs(t2)<1.d-8) then
            node(v1(1)).bw=2
            node(v1(1)).z=z
            inode1=v1(1)
            if(present(inode)) inode=inode1
            return
        elseif(abs(1.0d0-t2)<1.d-8) then
            node(v1(2)).bw=2
            node(v1(2)).z=z
            inode1=v1(2)
            if(present(inode)) inode=inode1            
            return
        else
            if(t2<0.d0.or.t2>1.d0) return !no cross
        endif
        
        
        
        xi1(1)=node(v1(1)).x+t2*(node(v1(2)).x-node(v1(1)).x)
        xi1(2)=node(v1(1)).y+t2*(node(v1(2)).y-node(v1(1)).y)
        xi1(3)=z

        ngnode=ngnode+1
        if(size(node)<ngnode) call enlarge_ar(node,100)
        node(ngnode).x=xi1(1);node(ngnode).y=xi1(2);node(ngnode).z=xi1(3);
        node(ngnode).bw=2 !mark
        node(ngnode).s=0.d0
        node(ngnode).iptr=ngnode
        gedge(iedge).cutpoint=ngnode
        call cut_edge(iedge,ngnode)
        inode1=ngnode
        if(present(inode)) inode=inode1
    
    
    endsubroutine
    
    
    subroutine cut_edge_by_seg_in_same_vertical_face(iedge,pt)
        implicit none
        integer,intent(in)::iedge
        real(8),intent(in)::pt(3,2)
        real(8)::v1(3,2),t1(2),t2,v2(3)
        integer::iv1(2)
        
        v1(1,1:2)=node(gedge(iedge).v).x
        v1(2,1:2)=node(gedge(iedge).v).y
        v1(3,1:2)=node(gedge(iedge).v).z
        
        if(abs(v1(1,1)-pt(1,1))>1.d-7) then
            if(abs(v1(1,2)-pt(1,1))>1.d-7) then
                error stop 'the two edges are not in a same vertical face1.'
            else
                if(abs(v1(1,1)-pt(1,2))>1.d-7) then
                    error stop 'the two edges are not in a same vertical face2.'
                else
                    iv1=[2,1]
                endif
            endif
        else
             if(abs(v1(1,2)-pt(1,2))>1.d-7) then
                error stop 'the two edges are not in a same vertical face3.'
             else
                 iv1=[1,2]                
            endif           
        endif
        
        t1=v1(3,iv1)-pt(3,:)
        if(t1(1)*t1(2)<0.d0) then
            t2=abs(t1(1))/(abs(t1(1))+abs(t1(2)))
            v2=pt(:,1)+t2*(pt(:,2)-pt(:,1))
            ngnode=ngnode+1
            if(size(node)<ngnode) call enlarge_ar(node,100)
            node(ngnode).x=v2(1);node(ngnode).y=v2(2);node(ngnode).z=v2(3);
            node(ngnode).bw=2 !mark
            node(ngnode).s=0.d0
            node(ngnode).iptr=ngnode
            gedge(iedge).cutpoint=ngnode
            call cut_edge(iedge,ngnode)
        endif
                    
            
        
        
    end subroutine
    
        
        subroutine cut_selt2(ielt,isx)
            implicit none
            integer,intent(in)::ielt,isx
            real(8)::ver1(3,3)
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
                
                !face(nface).normal=[0.d0,0.d0,1.d0] !assume the cut face is horizontal
                ver1(1,1:3)=node(node1(1:3)).x
                ver1(2,1:3)=node(node1(1:3)).y
                ver1(3,1:3)=node(node1(1:3)).z               
                face(nface).normal=NORMAL_TRIFACE(ver1,.true.)  
                
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

                    if(any(node(face(if1).node(:face(if1).nnum)).s<-1.d-8)) then
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
                iw1=pack(selt(ielt).node(:selt(ielt).nnum),node(selt(ielt).node(:selt(ielt).nnum)).s>1.d-8)
                selt(nselt+1).nnum=size(iw1)+face(nface).nnum
                if(size(selt(nselt+1).node)<selt(nselt+1).nnum) call enlarge_ar(selt(nselt+1).node,10)
                selt(nselt+1).node(:selt(nselt+1).nnum)=[iw1,face(nface).node(:face(nface).nnum)]
                selt(nselt+1).et=631
                iw1=pack(selt(ielt).node(:selt(ielt).nnum),node(selt(ielt).node(:selt(ielt).nnum)).s<-1.d-8)
                selt(nselt+2).nnum=size(iw1)+face(nface).nnum
                if(size(selt(nselt+2).node)<selt(nselt+2).nnum) call enlarge_ar(selt(nselt+2).node,10)
                selt(nselt+2).node(:selt(nselt+2).nnum)=[iw1,face(nface).node(:face(nface).nnum)]
                selt(nselt+2).et=632
                
                !check volume


                !selt(nselt+1).vol=selt(nselt+1).calvol()
                !selt(nselt+2).vol=selt(ielt).vol-selt(nselt+1).vol
                !if(selt(nselt+1).vol<1.d-7) write(*,*) '1the volume of element i is less than 1.d-7. i=', nselt+1
                !if(selt(nselt+2).vol<1.d-7) then
                !    write(*,*) '2the volume of element i is less than 1.d-7. i=', nselt+2
                !endif
                !polyhedron_volume_3d ( coord, order_max, face_num, node, node_num, order )
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
        

        subroutine cut_face2(iface)
            implicit none
            integer,intent(in)::iface
            !real(8),intent(in)::be
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
        if(ismerged>0) CALL merge_face()
        call cal_meshsize()
		outfile=trim(path_name)//'_geomodel.geo'
        call WRITE2GMSH_GEOFILEi2(outfile)
        outfile=trim(path_name)//'_geomodel.poly'
        CALL WRITE_3DPOLY_FILE(OUTFILE)
  !      OUTFILE2=trim(path_name)//'_offmodel.off'
		!CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE VOLUME GEO/POLY FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	LEN1=LEN_trim(term)
     	msg = MESSAGEBOXQQ(TRIM(term(1:LEN1)),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
		
    endsubroutine    
    
    
    SUBROUTINE WRITE2GMSH_GEOFILEi2(FILE_L)
		
		IMPLICIT NONE

		CHARACTER(LEN=*),INTENT(IN)::FILE_L
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::V1(2),allostat,METHOD1
		REAL(8)::T1,AT1(3),T2		
		INTEGER,ALLOCATABLE::IA1(:),IA2(:) 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3,CH4,CH5,CH6		
		CHARACTER(LEN=:),ALLOCATABLE::STR1
        
	         
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,10)
        
        T1=XYSCALE/20
		DO I=1,NGNODE
            
            IF(NODE(I).IPTR/=I.OR.NODE(I).ISB<1) CYCLE !重节点和内节点不输出            
			INC=INCOUNT(I)
            IF(ISMESHSIZE>0) THEN
                T2=NODE(I).S
            ELSE
                T2=T1
            ENDIF
			WRITE(50,20) I,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).Z*XYSCALE+ZMIN,T2
		ENDDO
		
		
		DO I=1,NGEDGE
            
			IF(Gedge(i).v(1)*Gedge(i).v(2)==0.OR.GEDGE(I).ISOVERLAP==1.OR.GEDGE(I).ISB<1.OR.GEDGE(I).ISDEL) cycle
			INC=INCOUNT(I)
            V1=GEDGE(I).V
            !DO J=1,2
            !    IF(NODE(GEDGE(I).V(J)).IPTR>0) V1(J)=NODE(GEDGE(I).V(J)).IPTR
            !ENDDO
			WRITE(50,30), I,V1            
		ENDDO			
		

		DO I=1,NFACE
            !IF(FACE(I).ISB/=9999) CYCLE
			IF (FACE(I).ISDEL.OR.FACE(I).ISB<1) CYCLE 
            N1=FACE(I).NNUM
			INC=INCOUNT(I)            
			NV1=N1-1
			WRITE(50,40) I,FACE(I).EDGE(1:N1)
			WRITE(50,50) I,I
        ENDDO
        
        !merged faces
        DO I=1,NGPF
            !IF(.NOT.ANY(IA1==NFACE+I)) CYCLE
            CALL GPFACET(I).facet2gmsh(nface,50)            
        ENDDO
        IF(NGPF>0) ALLOCATE(IA2(NGPF))

		!VOLUME 
		N1=0;
		DO I=1,ZNUM
            
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
            IF(.NOT.ALLOCATED(ZONE(I).VBFACE)) CYCLE
            
			DO J=1,SOILLAYER              
                
                WRITE(CH1,'(I4)') I
				WRITE(CH2,'(I4)') J
                CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))                
                CH4=TRIM(ADJUSTL(CH3))//'_HFACE'
                CH5=TRIM(ADJUSTL(CH3))//'_VFACE'
                NV1=N1
                
                DO K=1,ZONE(I).VBFACE(J).NVOL
                    !IF(ZONE(I).VBFACE(J).BFACE(K).NNUM<4) CYCLE
                    !IF(COUNT((SELT.ISDEL==.FALSE.).AND.(SELT.ZN==I).AND.(SELT.ILAYER==J))==0) CYCLE
				    N1=N1+1
				    INC=INCOUNT(N1)
                    !ZONE(I).VBFACE(J).BFACE(K).IVOL=N1 !VOLUME ID
                    !IF(N1/=19) CYCLE
                    !UPDATE TEH BFACES AFTER MERGE
                    IF(NGPF>0) THEN
                        CALL UPDATE_FACES_AFTER_MERGED(ZONE(I).VBFACE(J).BFACE(K).NODE,&
                                ZONE(I).VBFACE(J).BFACE(K).NNUM)                        
                    ENDIF
                    STR1=IAR2STR(ZONE(I).VBFACE(J).BFACE(K).NODE)
                    INC2=LEN_TRIM(STR1)
                    WRITE(50,61) N1,STR1(1:INC2)
					WRITE(50,70) N1,N1
			        
                ENDDO
                
                !CYCLE
                
                IF(N1-NV1==0) CYCLE
                
                LEN1=LEN(TRIM(ADJUSTL(CH3)))
                IF(N1-NV1==1) THEN
                    WRITE(50,82) TRIM(CH3),N1
                ELSEIF(N1-NV1>1) THEN
                    INC=INCOUNT(NV1+1);INC2=INCOUNT(N1)
                    WRITE(50,81) TRIM(CH3),NV1+1,N1
                ENDIF
                
                 
                
                IF(NGPF>0) THEN
                    !UPDATE outer horizontal faces AFTER MERGE
                    IA2=0
                    N2=ZONE(I).HFACE(J).NNUM
                    IF(N2>0) THEN
                        CALL UPDATE_FACES_AFTER_MERGED(ZONE(I).HFACE(J).NODE,&
                                ZONE(I).HFACE(J).NNUM)                        
                    ENDIF
                    !UPDATE outer VERTICAL faces AFTER MERGE
                    IA2=0
                    N2=ZONE(I).VFACE(J).NNUM
                    IF(N2>0) THEN
                        CALL UPDATE_FACES_AFTER_MERGED(ZONE(I).VFACE(J).NODE,&
                                ZONE(I).VFACE(J).NNUM)                        
                    ENDIF
                ENDIF
                IF(ZONE(I).HFACE(J).NNUM>0) THEN
                    STR1=IAR2STR(ZONE(I).HFACE(J).NODE)
                    LEN1=LEN(TRIM(ADJUSTL(CH4)))
                    INC2=LEN_TRIM(STR1)
                    WRITE(50,100) TRIM(CH4),STR1(1:INC2)
                    !COMPOUND SURFACE
                    IF(J==SOILLAYER.and.iscompounded==1) WRITE(50,120) STR1(1:INC2)
                        
                ENDIF
                IF(ZONE(I).VFACE(J).NNUM>0) THEN
                    STR1=IAR2STR(ZONE(I).VFACE(J).NODE)
                    LEN1=LEN(TRIM(ADJUSTL(CH5)))
                    INC2=LEN_TRIM(STR1)
                    WRITE(50,100) TRIM(CH5),STR1(1:INC2)
                ENDIF   
			    
            ENDDO
                            
        ENDDO        
        
        !RETURN
        
        DO I=1,NXZONE
            IF(XZONE(I).NCUTF<1) CYCLE    
            
            INC2=INCOUNT(I)
            CH1=""
			WRITE(CH1,91) I
            IF(XZONE(I).ISX/=2) THEN
                CH5='XZONE_'//TRIM(ADJUSTL(CH1))//'_CUTFACE'
                CALL UPDATE_FACES_AFTER_MERGED(XZONE(I).CUTFACE,XZONE(I).NCUTF)
            ELSE
                CH5='XZONE_'//TRIM(ADJUSTL(CH1))//'_CUTEDGE'
                !UPDATE EDGE AFTER MERGED EDGE
                IF(ISMERGED==1) THEN
                    XZONE(I).CUTFACE=PACK(XZONE(I).CUTFACE(1:XZONE(I).NCUTF),GEDGE(XZONE(I).CUTFACE).ISB>0)
                    XZONE(I).NCUTF=SIZE(XZONE(I).CUTFACE)
                ENDIF
            ENDIF
            
            STR1=IAR2STR(XZONE(I).CUTFACE)
            INC2=LEN_TRIM(STR1)                    
   !         IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
			!IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
            LEN1=LEN(TRIM(ADJUSTL(CH5)))
			IF(XZONE(I).NCUTF>0) THEN
                IF(XZONE(I).ISX/=2) THEN
					WRITE(50,100) TRIM(CH5),STR1(1:INC2)
                ELSE
                    WRITE(50,101) TRIM(CH5),STR1(1:INC2)
                ENDIF
                    
            ENDIF                
                
        ENDDO
            
        DO I=1,NCOW
            DO J=1,4
                IF(COWALL(I).FACE(J).NNUM<1) CYCLE

                INC2=INCOUNT(I)
                CH1=""
				WRITE(CH1,91) I
                IF(J==1) THEN
                    CH5='COW_'//TRIM(ADJUSTL(CH1))//'_CRACKFACE'
                ELSEIF(J==2) THEN
                    CH5='COW_'//TRIM(ADJUSTL(CH1))//'_WALL'
                ELSEIF(J==4) THEN
                    CH5='COW_'//TRIM(ADJUSTL(CH1))//'_CRACKEDGE'
                ELSE 
                    CH5='COW_'//TRIM(ADJUSTL(CH1))//'_BOUNDARY_EDGES'
                ENDIF
                IF(J<3) THEN
                    CALL UPDATE_FACES_AFTER_MERGED(COWALL(I).FACE(J).NODE, &
                    COWALL(I).FACE(J).NNUM)
                ELSEIF(ISMERGED==1) THEN
                    COWALL(I).FACE(J).NODE=PACK(COWALL(I).FACE(J).NODE,GEDGE(COWALL(I).FACE(J).NODE).ISB>0)
                    COWALL(I).FACE(J).NNUM=SIZE(COWALL(I).FACE(J).NODE)
                ENDIF
                STR1=IAR2STR(COWALL(I).FACE(J).NODE)
                INC2=LEN_TRIM(STR1)                    
    !            IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				!IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
                LEN1=LEN(TRIM(ADJUSTL(CH5)))
				!IF(COWALL(I).FACE(J).NNUM>0) THEN
                IF(J<3) THEN
                    WRITE(50,100) TRIM(CH5),STR1(1:INC2)  
                ELSE
                    WRITE(50,101) TRIM(CH5),STR1(1:INC2) 
                ENDIF
                !ENDIF              
            ENDDO
                
        ENDDO
            
        DO I=1,NVSEG_PG
            if(SIZE(VSEG_PG(I).EDGE)<1) cycle
            IF(ISMERGED==1) THEN
                VSEG_PG(I).IGMSHVOL=PACK(VSEG_PG(I).IGMSHVOL,GEDGE(VSEG_PG(I).EDGE).ISB>0)
                VSEG_PG(I).EDGE=PACK(VSEG_PG(I).EDGE,GEDGE(VSEG_PG(I).EDGE).ISB>0)                
            ENDIF

            DO J=1,SIZE(VSEG_PG(I).EDGE)
                INC2=INCOUNT(VSEG_PG(I).EDGE(J))
                IF(VSEG_PG(I).IGMSHVOL(J)>0) THEN                    
                    INC=INCOUNT(VSEG_PG(I).IGMSHVOL(J))
                    WRITE(50,110) VSEG_PG(I).EDGE(J),VSEG_PG(I).IGMSHVOL(J)
                ENDIF
            ENDDO            
            INC2=INCOUNT(I)
            WRITE(CH1,91) I
            CH5='VSEG_PG_'//TRIM(ADJUSTL(CH1)) 
            STR1=IAR2STR(VSEG_PG(I).EDGE)
            LEN1=LEN_TRIM(CH5)
            INC2=LEN_TRIM(STR1)                    
   !         IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
			!IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','            
            WRITE(50,101) CH5,STR1
        ENDDO

        DO I=1,NNODE_PG 
            
            !STR1=""
                      
            INC2=INCOUNT(I)
            WRITE(CH1,91) I
            CH5='NODE_PG_'//TRIM(ADJUSTL(CH1)) 
            STR1=IAR2STR(NODE_PG(I).NODE)
            LEN1=LEN_TRIM(CH5)
            INC2=LEN_TRIM(STR1)                    
   !         IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
			!IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','            
            WRITE(50,130) CH5,STR1
        ENDDO
        
        DO I=1,NVFACE_PG
            IF(VFACE_PG(I).NFACE<1) CYCLE    
            
            INC2=INCOUNT(I)
            CH1=""
            WRITE(CH1,91) I
            
            CH5='VFACE_PG_'//TRIM(ADJUSTL(CH1))
            CALL UPDATE_FACES_AFTER_MERGED(VFACE_PG(I).FACE,VFACE_PG(I).NFACE)
            
            
            STR1=IAR2STR(VFACE_PG(I).FACE)
            INC2=LEN_TRIM(STR1)                    
   
			LEN1=LEN(TRIM(ADJUSTL(CH5)))
			IF(VFACE_PG(I).NFACE>0) THEN
				WRITE(50,100) TRIM(CH5),STR1(1:INC2)
            ENDIF                
                
        ENDDO 

       DO I=1,NVOL_PG 
            CALL VOL_PG(I).GETVOL()
            
            IF(VOL_PG(I).ngvol<1) CYCLE    
            
            INC2=INCOUNT(I)
            CH1=""
            WRITE(CH1,91) I
            
            CH5='USER_VOL_PG'//TRIM(ADJUSTL(CH1))
            
            STR1=IAR2STR(VOL_PG(I).GVOL)
            INC2=LEN_TRIM(STR1)                    
   
			LEN1=LEN(TRIM(ADJUSTL(CH5)))
			IF(VOL_PG(I).NGVOL>0) THEN
				WRITE(50,140) TRIM(CH5),STR1(1:INC2)
            ENDIF                
                
        ENDDO      
        
        DO I=1,NBLO
            CALL BLO(i).write(50)
        ENDDO
        
		CLOSE(50)
        
        DEALLOCATE(IA1,IA2,stat=allostat)

	10  FORMAT('SetFactory("Built-in");Geometry.OCCAutoFix = 0;Mesh.AngleToleranceFacetOverlap=1.e-6;Mesh.MshFileVersion=2.2;', /, 'ms=1;')
	20	FORMAT("Point(",I<INC>,")={",3(E24.15,","),E24.15,"*ms};")
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
    102 FORMAT('Physical Line("',A<LEN1>,'"',")={",I<INC2>"};") 
110     FORMAT('Curve {',I<INC2>,'} In Volume {',I<INC>,'};') 
120     FORMAT('Compound Surface{',A<INC2>,'};')
130     FORMAT('Physical Point("',A<LEN1>,'"',")={",A<INC2>"};") 
140     FORMAT('Physical Volume("',A<LEN1>,'"',")={",A<INC2>"};")       
    ENDSUBROUTINE	    
    
    SUBROUTINE UPDATE_FACES_AFTER_MERGED(IA,NIA)
        IMPLICIT NONE
        INTEGER,DIMENSION(:),ALLOCATABLE::IA,IA1,IA2
        INTEGER::NIA,ALLOSTAT
        
        
        IF(NIA>0) THEN
            ALLOCATE(IA2(NGPF))
            IA2=0
            IA1=IA(1:NIA)
            WHERE(FACE(IA1).ISB<0) IA2(-FACE(IA1).ISB)=NFACE-FACE(IA1).ISB
            IA=[PACK(IA1,FACE(IA1).ISB>0),PACK(IA2,IA2>0)]
            NIA=SIZE(IA)  
        ENDIF
    
        DEALLOCATE(IA1,IA2,stat=allostat)
        
    END SUBROUTINE
    
    SUBROUTINE WRITE_3DPOLY_FILE(FILE_L)
		
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
        
	    if(poly3d==0) return
        
        
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		N1=0
		DO I=1,NGNODE
            IF(NODE(I).IPTR/=I.OR.NODE(I).ISB<1) CYCLE !重节点和内节点不输出
            N1=N1+1
        ENDDO
        WRITE(50,20) N1,0
        N1=0
		DO I=1,NGNODE
            IF(NODE(I).IPTR/=I.OR.NODE(I).ISB<1) CYCLE !重节点和内节点不输出
            N1=N1+1
            NODE(I).BW=N1
			WRITE(50,30) N1,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).Z*XYSCALE+ZMIN
        ENDDO
        
		!FACE.MARKER 
		DO I=1,ZNUM
            
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
            IF(.NOT.ALLOCATED(ZONE(I).BFACE)) CYCLE
			DO J=1,SOILLAYER
   
                
                !volume outer horizontal faces
                N3=ZONE(I).BFACE(J).NNUM
                IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB/=0 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==0 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ILAYER/=0)
              
                IF(SIZE(IA1)>0) THEN
                    FACE(IA1).MARKER=I*1000+J*10+0 !IZONE*1000+ILAYER*10+0/1
                ENDIF
                
                !volume outer vertical faces
                IA1=PACK(ZONE(I).BFACE(J).NODE(:N3),FACE(ZONE(I).BFACE(J).NODE(:N3)).ISB/=0 &
                        .AND.FACE(ZONE(I).BFACE(J).NODE(:N3)).ISV==1)
               
                IF(SIZE(IA1)>0) THEN
				    FACE(IA1).MARKER=I*1000+J*10+1 
                ENDIF                
			    
            ENDDO
                
        ENDDO      
        DO I=1,NCOW
            FACE(COWALL(I).FACE(1).NODE(:COWALL(I).FACE(1).NNUM)).MARKER=-I            
        ENDDO        
        
        WRITE(50,'(A)') '#FACET DATA.MARKER=imarker'
        WRITE(50,'(A)') '#when imarker>0,imarker=IZONE*1000+ELAYER*10+a(a=0,zone horizontal boundaries;a=1,vertical boundaries.)'
        WRITE(50,'(A)') '#when imarker<0,imarker=-a(the a-th cut off wall boundaries.)'
        N1=COUNT((.NOT.FACE(:NFACE).ISDEL).AND.FACE(:NFACE).ISB>0)+NGPF
        WRITE(50,40) N1
		
        N3=0
		DO I=1,NFACE
			IF (FACE(I).ISDEL.OR.FACE(I).ISB<1) CYCLE
            N3=N3+1
            WRITE(50,42) 1,0,FACE(I).MARKER,N3            
            N1=FACE(I).NNUM
			WRITE(50,41) N1,NODE(FACE(I).NODE(1:N1)).BW
        ENDDO
        !merged facets
        DO I=1,NGPF            
            WRITE(50,421) GPFACET(I).NPOLY+GPFACET(I).NFP,GPFACET(I).NHOLE,FACE(GPFACET(I).ONEFACE).MARKER,I
            DO J=1,GPFACET(I).NPOLY
                N1=GPFACET(I).POLY(J).NNUM
			    WRITE(50,41) N1,NODE(GPFACET(I).POLY(J).NODE(1:N1)).BW
            ENDDO
            DO J=1,GPFACET(I).NFP
                N1=1
			    WRITE(50,41) N1,NODE(GPFACET(I).FLOATP(J)).BW                
            ENDDO
            DO J=1,GPFACET(I).NHOLE                
			    WRITE(50,43) J,GPFACET(I).HOLE(1,J)*XYSCALE+XMIN, GPFACET(I).HOLE(2,J)*XYSCALE+YMIN,GPFACET(I).HOLE(3,J)*XYSCALE+ZMIN
            ENDDO            
        ENDDO
        
        
        !hole
        WRITE(50,'(A)') '#HOLE DATA'
        WRITE(50,'(I3)') 0
        
        !REGION        

		N1=0;
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
			DO J=1,SOILLAYER
                IF(COUNT((SELT.ISDEL==.FALSE.).AND.(SELT.ZN==I).AND.(SELT.ILAYER==J))==0) CYCLE
                
                DO K=1,ZONE(I).VBFACE(J).NVOL
                    IF(ZONE(I).VBFACE(J).BFACE(K).NNUM<4) CYCLE
                    N1=N1+1
                ENDDO
            ENDDO
        ENDDO
        
        WRITE(50,'(A)') '#REGION DATA,REGION NUMBER=IZONE*1000+ELAYER*10'
        WRITE(50,'(I3)') N1
        N1=0
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
			DO J=1,SOILLAYER
                IF(COUNT((SELT.ISDEL==.FALSE.).AND.(SELT.ZN==I).AND.(SELT.ILAYER==J))==0) CYCLE
                
                DO K=1,ZONE(I).VBFACE(J).NVOL
                    IF(ZONE(I).VBFACE(J).BFACE(K).NNUM<4) CYCLE
                    N1=N1+1
                    WRITE(50,'(I3,1X,3(F15.7,1X),I7)') N1,ZONE(I).VBFACE(J).BFACE(K).INSIDEPT(1)*xyscale+xmin,&
                        ZONE(I).VBFACE(J).BFACE(K).INSIDEPT(2)*xyscale+ymin,ZONE(I).VBFACE(J).BFACE(K).INSIDEPT(3)*XYSCALE+ZMIN,I*1000+J*10 
                ENDDO
            ENDDO
        ENDDO         
         
        
            

        
		CLOSE(50)
        


10      FORMAT(A<len1>)  
20      FORMAT(I7,1X,'3',1X,I3,1X,'0')
30      FORMAT(I7,1X,<3>(F15.7,1X))
31      FORMAT(I7,1X,<2>(F15.7,1X),<N1>('-999.0',1X),2(I3,1X)) 
32      FORMAT(I7,1X,<2+N1+1>(F15.7,1X),2(I3,1X))
33      FORMAT(I7,1X,<3>(F15.7,1X),<N1>('-999.0',1X),2(I3,1X))         
40      FORMAT(I7,1X,'1')  
41      FORMAT(<MIN(N1+1,50)>(I7,1X))
42      FORMAT(3(I7,1X),'#FACET ',I4)
421     FORMAT(3(I7,1X),'#Merged FACET ',I4)        
43      FORMAT(I7,1X,3(F15.7,1X))         
50      FORMAT(I7,2(F15.7,1X))
60      FORMAT(I7,2(F15.7,1X),I) 
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
                    WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,T1*XYSCALE+ZMIN,10.0
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
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION(I)*XYSCALE+ZMIN,10.0
                        else
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION((I-1)*THIS.BLOOP(1).NNUM+J)*XYSCALE+ZMIN,10.0
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
                this.bloop(i).elevation=(ar(1:n1)-Zmin)/xyscale 
            else
                this.bloop(i).elevation=(ar(1:2*n1)-Zmin)/xyscale 
            endif
            
        enddo
        
     
    endsubroutine


    
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
            NGMSHVOL=NGMSHVOL+1
            IF(NVBF1>SIZE(ZONE(IZONE).VBFACE(ILAYER).BFACE)) CALL ENLARGE_AR(ZONE(IZONE).VBFACE(ILAYER).BFACE,5)
            !该体内的一点，取该体内一单元的形心
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).INSIDEPT(1)=SUM(NODE(SELT(STACK1(1)).NODE(:SELT(STACK1(1)).NNUM)).X)/SELT(STACK1(1)).NNUM
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).INSIDEPT(2)=SUM(NODE(SELT(STACK1(1)).NODE(:SELT(STACK1(1)).NNUM)).Y)/SELT(STACK1(1)).NNUM
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).INSIDEPT(3)=SUM(NODE(SELT(STACK1(1)).NODE(:SELT(STACK1(1)).NNUM)).Z)/SELT(STACK1(1)).NNUM
            !蚕食法
            DO WHILE(N2>0)
                IELT1=STACK1(N2)                
                N2=N2-1
                IF(ISCHK1(IELT1)==1) CYCLE
                ISCHK1(IELT1)=1
                SELT(IELT1).IGMSHVOL=NGMSHVOL
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
            ZONE(IZONE).VBFACE(ILAYER).BFACE(NVBF1).IVOL=NGMSHVOL
        ENDDO
        ZONE(IZONE).VBFACE(ILAYER).NVOL=NVBF1
   
        DEALLOCATE(GELT1,STAT=ALLOSTAT)
   
    ENDSUBROUTINE
    

    
function polyhedron_volume_3d ( coord, order_max, face_num, node, &
  node_num, order ) result(volume)

!*****************************************************************************80
!
!! POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COORD(3,NODE_NUM), the coordinates of 
!    the vertices.  The vertices may be listed in any order.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer ( kind = 4 ) NODE(FACE_NUM,ORDER_MAX).  Face I is defined by
!    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ( kind = 4 ) ORDER(FACE_NUM), the number of vertices making 
!    up each face.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the polyhedron.
!https://stackoverflow.com/questions/1838401/general-formula-to-calculate-polyhedron-volume
!1.Take the polygons and break them into triangles.
!2.Consider the tetrahedron formed by each triangle and an arbitrary point (the origin).
!3.Sum the signed volumes of these tetrahedra.
!Notes:
!
!(1)This will only work if you can keep a consistent CW or CCW order to the triangles as viewed from the outside.
!(2)The signed volume of the tetrahedron is equal to 1/6 the determinant of the following matrix:
![ x1 x2 x3 x4 ]
![ y1 y2 y3 y4 ]
![ z1 z2 z3 z4 ]
![ 1 1 1 1 ]
!
!where the columns are the homogeneous coordinates of the verticies (x,y,z,1).
!
!It works even if the shape does not enclose the origin by subracting off that volume as well as adding it in, 
!but that depends on having a consistent ordering.
!
!If you can't preserve the order you can still find some way to break it into tetrahedrons and sum 1/6 absolute v
!alue of the determinant of each one.
!
!Edit: I'd like to add that for triangle mesh where one vertex (say V4) of the tetrahedron is (0,0,0) the determinante 
!of the 4x4 matrix can be simplified to the upper left 3x3 (expansion along the 0,0,0,1 column) and that can be 
!simplified to Vol = V1xV2.V3 where "x" is cross product and "." is dot product. So compute that expression for 
!every triangle, sum those volumes and divide by 6.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) coord(dim_num,node_num)
  integer ( kind = 4 ) face
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) node(face_num,order_max)
  integer ( kind = 4 ) order(face_num)
  integer ( kind = 4 ) v
  real ( kind = 8 ) volume

  volume = 0.0D+00
!
!  Triangulate each face.
!
  do face = 1, face_num

    n3 = node(face,order(face))

    do v = 1, order(face) - 2

      n1 = node(face,v)
      n2 = node(face,v+1)

      volume = volume &
        + coord(1,n1) &
        * ( coord(2,n2) * coord(3,n3) - coord(2,n3) * coord(3,n2) ) &
        + coord(1,n2) &
        * ( coord(2,n3) * coord(3,n1) - coord(2,n1) * coord(3,n3) ) &
        + coord(1,n3) &
        * ( coord(2,n1) * coord(3,n2) - coord(2,n2) * coord(3,n1) )

    end do

  end do

  volume = abs(volume / 6.0D+00)
  
  return
end function    
    
function polyhedron_volume_3d_2 ( coord, order_max, face_num, node, &
  node_num, order ) result(volume)

!*****************************************************************************80
!
!! POLYHEDRON_VOLUME_3D_2 computes the volume of a polyhedron in 3D.
!
!  Discussion:
!
!    The computation is not valid unless the faces of the polyhedron
!    are planar polygons.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    AP Professional, 1995, T385.G6975.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COORD(3,NODE_NUM), the vertices.
!    The vertices may be listed in any order.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer ( kind = 4 ) NODE(FACE_NUM,ORDER_MAX).  Face I is defined 
!    by the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ( kind = 4 ) ORDER(FACE_NUM), the number of vertices making 
!    up each face.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the polyhedron.
!
  implicit none

  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) coord(dim_num,node_num)
  integer ( kind = 4 ) face
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) node(face_num,order_max)
  real ( kind = 8 ) normal(dim_num)
  integer ( kind = 4 ) order(face_num)
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) volume

  volume = 0.0D+00

  do face = 1, face_num

    v(1:dim_num) = 0.0D+00
!
!  Compute the area vector for this face.
!
    do j = 1, order(face)

      k1 = node(face,j)

      if ( j < order(face) ) then
        k2 = node(face,j+1)
      else
        k2 = node(face,1)
      end if
!
!  Compute the cross product.
!
      normal(1) = coord(2,k1) * coord(3,k2) - coord(3,k1) * coord(2,k2)
      normal(2) = coord(3,k1) * coord(1,k2) - coord(1,k1) * coord(3,k2)
      normal(3) = coord(1,k1) * coord(2,k2) - coord(2,k1) * coord(1,k2)

      v(1:dim_num) = v(1:dim_num) + normal(1:dim_num)

    end do
!
!  Area vector dot any vertex.
!
    k = node(face,1)
    volume = volume + dot_product ( v(1:dim_num), coord(1:dim_num,k) )

  end do

  volume =abs( volume / 6.0D+00)

  return
  end    
  
  
    function NORMAL_TRIFACE(V,isunify) result (Normal)

    !*****************************************************************************80
    !V, XY OF THE FACET.
    !CALCULATE THE NORMAL VECTOR OF A TRI-FACET. 
    !

    !
      implicit none

      REAL(8),INTENT(IN)::V(:,:) !3*3
      logical,optional::isunify
      real ( kind = 8 ) v1(SIZE(V,DIM=1))
      real ( kind = 8 ) v2(SIZE(V,DIM=1))
      real ( kind = 8 ) normal(SIZE(V,DIM=1))
      logical::isunify1=.true.
      
      if(.not.present(isunify)) then
          isunify1=.true.
      else
          isunify1=isunify
      endif
      
      V1=V(:,2)-V(:,1);V2=V(:,3)-V(:,1);

      normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
      normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
      normal(3) = v1(1) * v2(2) - v1(2) * v2(1)
      if(isunify1) normal=normal/norm2(normal)
      return
  
    end function
    
    SUBROUTINE GPFACET_ENLARGE_AR(AVAL,DSTEP)
        TYPE(general_poly_face_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(general_poly_face_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE  
    
function triangle_orientation_2d ( t )

!*****************************************************************************80
!
!! TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
!
!  Discussion:
!
!    Three distinct non-colinear points in the plane define a circle.
!    If the points are visited in the order P1, P2, and then
!    P3, this motion defines a clockwise or counter clockwise
!    rotation along the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, integer ( kind = 4 ) TRIANGLE_ORIENTATION_2D, reports if the 
!    three points lie clockwise on the circle that passes through them.  
!    The possible return values are:
!    0, the points are distinct, noncolinear, and lie counter clockwise
!    on their circle.
!    1, the points are distinct, noncolinear, and lie clockwise
!    on their circle.
!    2, the points are distinct and colinear.
!    3, at least two of the points are identical.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) det
  integer ( kind = 4 ) triangle_orientation_2d
  real ( kind = 8 ) t(dim_num,3)

  if ( all ( t(1:dim_num,1) == t(1:dim_num,2) ) .or. &
       all ( t(1:dim_num,2) == t(1:dim_num,3) ) .or. &
       all ( t(1:dim_num,3) == t(1:dim_num,1) ) ) then
    triangle_orientation_2d = 3
    return
  end if

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  if ( det == 0.0D+00 ) then
    triangle_orientation_2d = 2
  else if ( det < 0.0D+00 ) then
    triangle_orientation_2d = 1
  else if ( 0.0D+00 < det ) then
    triangle_orientation_2d = 0
  end if

  return
end

subroutine cal_meshsize()
    implicit None
    integer::I,J,V1(2)
    REAL(8)::T1,AR1(3,2)

    NODE(1:NGNODE).S=XYSCALE/20
    DO I=1,NGEDGE
        IF(Gedge(i).v(1)*Gedge(i).v(2)==0.OR.GEDGE(I).ISOVERLAP==1.OR.GEDGE(I).ISB<1.OR.GEDGE(I).ISDEL) cycle
        
        V1=GEDGE(I).V
        DO J=1,2
            AR1(1,J)=NODE(V1(J)).X*XYSCALE+XMIN
            AR1(2,J)=NODE(V1(J)).Y*XYSCALE+YMIN
            AR1(3,J)=(NODE(V1(J)).Z*XYSCALE+ZMIN)*10
        ENDDO
        T1=MAX(NORM2(AR1(:,1)-AR1(:,2)),xyscale/200.)
        WHERE(NODE(V1).S>T1) NODE(V1).S=T1
        
    ENDDO	   
    
end subroutine

   
 !   
 !   subroutine out_volume_to_gmsh()
	!	character*256 term
	!	integer(4)::msg
	!	character(512)::outfile,OUTFILE2
	!	
 !       call st_membrance()
 !       call find_zone_boudary()
 !       call node_on_boundary()        
 !       CALL CHECK_ORIENT()
	!	outfile=trim(path_name)//'_geomodel.geo'
 !       call WRITE2GMSH_GEOFILEi(outfile)
 !       
 !       OUTFILE2=trim(path_name)//'_offmodel.off'
	!	CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
 !       
	!	term="THE VOLUME GEO FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	! 	term=trim(term)
 !    	msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
 !    	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
	!	
 !   endsubroutine
 !   
	!subroutine node_on_boundary()
	!	integer::i
	!	do i=1,nedge
 !           IF(edge(i).v(1)*edge(i).v(2)==0) cycle
 !           if(edge(i).iszonebc<0) cycle            
	!		node(edge(I).v).onbdy=1		
	!	enddo
 !       
	!	nbnode=count(node.onbdy==1)
	!	allocate(bnode(nbnode),BEDGE(NNODE),ZBEDGE(SOILLAYER,NNODE))
 !       ZBEDGE=0
	!	bnode=pack(node(1:nnode).number,node(1:nnode).onbdy==1)
	!end subroutine
	!
 !   subroutine find_zone_boudary()
	!	integer::izone
 !       integer::i,j,k,ielt,iedge,n1
	!	integer::bedge1(5000),nbe=0
	!	
 !       do izone=1,znum
 !       
 !           if(zone(izone).outgmshtype==2) cycle
 !           
	!	    nbe=0
 !           n1=zone(izone).ntrie3n
	!	    do i=1,zone(izone).ntrie3n
	!		    ielt=zone(izone).trie3n(i)			
	!		    do j=1,3
	!			    iedge=elt(ielt).edge(j)
	!			    if(edge(iedge).e(1)+1==0.or.edge(iedge).e(2)+1==0) then
	!				    nbe=nbe+1
	!				    bedge1(nbe)=iedge
	!			    else
	!				    if(elt(edge(iedge).e(1)).zn/=elt(edge(iedge).e(2)).zn) then
	!					    nbe=nbe+1
	!					    bedge1(nbe)=iedge
	!				    endif
	!			    endif
	!		    enddo
	!	    enddo
	!	    if(nbe>0) then
 !               IF(ALLOCATED(zone(izone).bedge)) DEALLOCATE(zone(izone).bedge)
	!		    allocate(zone(izone).bedge,source=bedge1(1:nbe))
 !               EDGE(bedge1(1:nbe)).ISZONEBC=IZONE
 !               ZONE(IZONE).NBE=NBE
 !               !DO J=1,NBE
 !               !    IF(EDGE(bedge1(J)).ISCEDGE==0) THEN                        
 !               !        PRINT *, 'SOME WRONG.'
 !               !    ENDIF
 !               !ENDDO
 !
	!	    endif
	!    enddo
 !   endsubroutine
 !   
 !   SUBROUTINE OUT_OFF_MATHER_MODEL(FILE_L)
 !       IMPLICIT NONE
 !       
 !       CHARACTER(LEN=*),INTENT(IN)::FILE_L
 !       INTEGER::I
 !       
 !       CALL GEN_OFF_MATHER_MODEL()
 !       
	!	OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
	!	WRITE(50,'(A3)') 'OFF'
 !       WRITE(50,'(3I7)') OFF_MATHER_MODEL.NNODE,OFF_MATHER_MODEL.NFACE,0
 !       DO I=1,OFF_MATHER_MODEL.NNODE
 !           WRITE(50,'(3(E15.7,1X))') OFF_MATHER_MODEL.NODE(:,I)
 !       ENDDO
 !       DO I=1,OFF_MATHER_MODEL.NFACE
 !           WRITE(50,'(4(I7,1X))') 3,OFF_MATHER_MODEL.FACE(:,I)-1
 !       ENDDO
 !       CLOSE(50)
 !       
 !   ENDSUBROUTINE
 !   
 !   SUBROUTINE GEN_OFF_MATHER_MODEL()
 !       IMPLICIT NONE
 !       
 !       INTEGER::I,J,N1,N2,N3
 !       
 !        
 !       !NNODE
 !       OFF_MATHER_MODEL.NNODE=NNODE*(SOILLAYER+1)
 !       ALLOCATE(OFF_MATHER_MODEL.NODE(3,OFF_MATHER_MODEL.NNODE))
 !       
	!	DO J=0,SOILLAYER
	!		N1=NNODE*J
	!		DO I=1,NNODE
	!			N2=N1+I
	!			OFF_MATHER_MODEL.NODE(:,N2)=[NODE(N2).X*XYSCALE+XMIN,NODE(N2).Y*XYSCALE+YMIN,NODE(N2).Z]
	!		ENDDO
	!	ENDDO
 !       !NFACE
 !       N1=COUNT(EDGE.ISZONEBC>0)
 !       OFF_MATHER_MODEL.NFACE=ENUMBER*(SOILLAYER+1)+N1*SOILLAYER*2
 !       ALLOCATE(OFF_MATHER_MODEL.FACE(3,OFF_MATHER_MODEL.NFACE))
 !       		!honrizontal triangle
	!	N1=0
	!	DO J=0,SOILLAYER
	!		N2=NNODE*J;
	!		DO I=1,NELT
	!			IF (ELT(I).ISDEL) CYCLE
	!			N1=N1+1			
 !               OFF_MATHER_MODEL.FACE(:,N1)=ELT(I).NODE(1:3)+N2
	!		ENDDO		
	!	ENDDO
	!	!vertial boundary rectangle
	!	DO J=1,SOILLAYER
	!		N2=NNODE*(J-1);N3=NNODE*J
	!		DO I=1,NEDGE
 !               IF(edge(i).v(1)*edge(i).v(2)==0) cycle
 !               IF(EDGE(I).ISZONEBC<0) CYCLE
	!			N1=N1+1			
	!			OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V+N2,EDGE(I).V(2)+N3]
 !               N1=N1+1			
	!			OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V(2:1:-1)+N3,EDGE(I).V(1)+N2]
	!		ENDDO		
	!	ENDDO	
 !   
	!END SUBROUTINE
    
	!SUBROUTINE WRITE2GMSH_GEOFILEi(FILE_L)
	!	
	!	IMPLICIT NONE
 !
	!	CHARACTER(LEN=*),INTENT(IN)::FILE_L
	!	
	!	INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
	!	INTEGER::NEDGE1=0,NBCEDGE1=0
	!	REAL(8)::T1,AT1(3)
	!	REAL(8),ALLOCATABLE::MINDIS1(:)
	!	INTEGER,ALLOCATABLE::IA1(:),IA2(:) 
	!	CHARACTER(512)::GEOFILE1
	!	CHARACTER(32)::CH1,CH2,CH3		
	!	CHARACTER(8192*3)::STR1
 !       
	!	
 !
	!	
	!	OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
	!	WRITE(50,10)
 !       
 !       do i=1,nnode
 !           if(node(i).havesoildata==0) cycle
 !           do j=soillayer,1,-1
 !               if(node(i).elevation(j)-node(i).elevation(j-1)<0.01) node(i).elevation(j-1)=node(i).elevation(j)-0.01          
 !           enddo
 !       enddo
	!	!MINIMAL DISTANCE BETWEEN NODES
	!	ALLOCATE(MINDIS1(NNODE))
	!	MINDIS1=1.D20
	!	DO I=1,NNODE			
	!		DO J=I+1,NNODE				
	!			AT1=[(NODE(I).X-NODE(J).X)*XYSCALE,(NODE(I).Y-NODE(J).Y)*XYSCALE,0.d0]
	!			T1=MIN(XYSCALE/20,MAX(NORM2(AT1),5.0)) !SET MINDIS>0.1
	!			IF (MINDIS1(I)>T1) MINDIS1(I)=T1
	!			IF (MINDIS1(J)>T1) MINDIS1(J)=T1
	!		ENDDO
	!	ENDDO
	!	
	!	DO J=0,SOILLAYER
	!		N1=NNODE*J
	!		DO I=1,NNODE
	!			N2=N1+I
 !               IF(NODE(N2).IPTR/=N2) CYCLE !重节点不输出
	!			INC=INCOUNT(N2)
	!			WRITE(50,20) N2,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).elevation(j),MINDIS1(I)
	!		ENDDO
	!	ENDDO
	!	
	!	N1=0;NEDGE1=0
 !       NBCEDGE1=COUNT(EDGE.ISZONEBC>0)
	!	DO J=0,SOILLAYER
	!		N2=NNODE*J
	!		DO I=1,NEDGE
	!			IF(edge(i).v(1)*edge(i).v(2)==0) cycle
	!			IF(J==0) THEN
	!				NEDGE1=NEDGE1+1					
 !               ENDIF
	!			N1=N1+1
 !               !XY坐标重合的边也不输出，无厚度防渗墙不输出。
 !               IF(EDGE(I).ISXYOVERLAP==0) THEN
	!			    INC=INCOUNT(N1)
	!			    WRITE(50,30), N1,EDGE(I).V+N2
 !               ENDIF
	!		ENDDO			
	!	ENDDO
	!	!VERTICAL BOUNDARY EDGE 
	!	DO J=1,SOILLAYER
	!		N3=NNODE*(J-1);N4=NNODE*J
	!		DO I=1,NBNODE
	!			N1=N1+1
	!			IF(J==1) BEDGE(BNODE(I))=N1
 !               !0长度的边不输出
	!			IF((NODE(BNODE(I)+N3).IPTR==NODE(BNODE(I)+N4).IPTR)) THEN
 !                   ZBEDGE(J,BNODE(I))=1
 !               ELSE
 !                   INC=INCOUNT(N1)
	!			    WRITE(50,30), N1,BNODE(I)+N3,BNODE(I)+N4
 !               ENDIF
	!		ENDDO
	!	ENDDO
	!	
	!	!honrizontal triangle
	!	N1=0
 !       ALLOCATE(IA1(NELT*(SOILLAYER+1)+NBCEDGE1*SOILLAYER))
 !       IA1=0
	!	DO J=0,SOILLAYER
	!		N2=NEDGE1*J;
	!		DO I=1,NELT
	!			IF (ELT(I).ISDEL) CYCLE 
	!			N1=N1+1
 !               IF (ELT(I).ET==-1) CYCLE !无厚度防渗墙单元也不输出
	!			INC=INCOUNT(N1)
	!			NV1=2
	!			WRITE(50,40) N1,(EDGE(ELT(I).EDGE(1:ELT(I).NNUM)).NUM+N2)*ELT(I).ORIENT(1:ELT(I).NNUM)
	!			WRITE(50,50) N1,N1
 !               IA1(N1)=1
	!		ENDDO		
	!	ENDDO
	!	!vertial boundary rectangle
	!	DO J=1,SOILLAYER
	!		N2=NEDGE1*(J-1);N3=NBNODE*(J-1)
	!		DO I=1,NEDGE
 !               IF(edge(i).v(1)*edge(i).v(2)==0) cycle
 !               IF(EDGE(I).ISZONEBC<0) CYCLE                                
	!			N1=N1+1
 !               IF(J==1) EDGE(I).BRECT=N1
	!			
 !               AN1(1)=EDGE(I).NUM+N2
 !               
 !               IF(ZBEDGE(J,EDGE(I).V(2))==0) THEN
 !                   AN1(2)=BEDGE(EDGE(I).V(2))+N3
 !                   AN1(3)=-(EDGE(I).NUM+N2+NEDGE1)
 !                   N4=3
 !               ELSE
 !                   AN1(2)=-(EDGE(I).NUM+N2+NEDGE1)
 !                   N4=2
 !               ENDIF
 !               IF(ZBEDGE(J,EDGE(I).V(1))==0) THEN
 !                   N4=N4+1
 !                   AN1(N4)=-(BEDGE(EDGE(I).V(1))+N3)
 !               ENDIF
 !               IF(N4>2) THEN
 !                   NV1=N4-1
 !                   INC=INCOUNT(N1) 
	!			    WRITE(50,40) N1,AN1(1:N4)
	!			    WRITE(50,50) N1,N1
 !                   IA1(N1)=1 !MARK OUTPUT FACES
 !               ENDIF
	!			
	!		ENDDO		
	!	ENDDO		
	!	
	!	!VOLUME 
	!	N1=0;
	!	DO I=1,ZNUM
 !           IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
	!		DO J=0,SOILLAYER-1
	!			N1=N1+1
	!			INC=INCOUNT(N1)
	!			WRITE(CH1,'(I4)') I
	!			WRITE(CH2,'(I4)') J+1
	!			CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))
	!			LEN1=LEN(TRIM(ADJUSTL(CH3)))	
	!		    !n4=ZONE(I).NTRIE3N
 !               INC3=0
 !               STR1=""
	!			DO K=1,ZONE(I).NTRIE3N
 !                   
	!				N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
	!				DO K1=1,2
	!					IF(K1==1) THEN
	!						N2=ENUMBER*J+N3; !ENUMBER INCLUDING THE ZERO-THICKNESS ELEMENT.
	!					ELSE
	!						N2=ENUMBER*(J+1)+N3;
 !                       ENDIF
 !                       
 !                       IF(IA1(N2)==0) CYCLE !SKIP THE FACE NOT OUTPUT
 !                       
	!					INC2=INCOUNT(N2)
	!					CH1=""
	!					WRITE(CH1,90) N2                                            
	!					STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
 !                       INC3=INC3+1    
 !                       IF(MOD(INC3,100)==0) THEN
 !                           STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
 !                       ENDIF
	!				ENDDO
	!			ENDDO
	!			DO K=1,ZONE(I).NBE
	!				N2=NBCEDGE1*(J)+EDGE(ZONE(I).BEDGE(K)).BRECT
 !                   
 !                   IF(IA1(N2)==0) CYCLE
 !                   
	!				INC2=INCOUNT(N2)
	!				CH1=""
	!				WRITE(CH1,90) N2
	!				STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
 !                   INC3=INC3+1    
 !                   IF(MOD(INC3,100)==0) THEN
 !                       STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
 !                   ENDIF
 !               ENDDO			
	!			!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
	!			!               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
 !               INC2=LEN_TRIM(STR1)
 !               IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
	!			IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
	!			IF(ZONE(I).NBE>0) THEN
	!				WRITE(50,61) N1,STR1(1:INC2)
	!				WRITE(50,70) N1,N1
	!				WRITE(50,80) TRIM(CH3),N1,N1
	!			ENDIF
	!		
	!		ENDDO
	!	ENDDO
	!	
	!	N1=0;INC3=0
	!	DO I=1,ZNUM
 !           IF(ZONE(I).OUTGMSHTYPE<2) CYCLE
	!		J=ZONE(I).IELEVATION
	!		N1=N1+1
	!		INC=INCOUNT(N1)
	!		WRITE(CH1,'(I4)') I
	!		WRITE(CH2,'(I4)') J
	!		CH3='Z'//TRIM(ADJUSTL(CH1))//'_E'//TRIM(ADJUSTL(CH2))
	!		LEN1=LEN(TRIM(ADJUSTL(CH3)))	
	!		!n4=ZONE(I).NTRIE3N
 !           STR1=""
	!		DO K=1,ZONE(I).NTRIE3N                    
	!			N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
	!			N2=ENUMBER*(J-1)+N3;
 !               
 !               IF(IA1(N2)==0) CYCLE
 !               
	!			INC2=INCOUNT(N2)
	!			CH1=""
	!			WRITE(CH1,90) N2
	!			STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
 !               INC3=INC3+1    
 !               IF(MOD(INC3,100)==0) THEN
 !                   STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
 !               ENDIF                
 !           ENDDO
 !           INC2=LEN_TRIM(STR1)
 !           IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
	!		IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','		
	!		WRITE(50,100) TRIM(CH3),N1,STR1(1:INC2)				
 !       ENDDO        
 !       
 !       DO I=1,NBLO
 !           CALL BLO(i).write(50)
 !       ENDDO
 !       
	!	CLOSE(50)
 !       
 !       DEALLOCATE(MINDIS1)
 !
	!10  FORMAT('SetFactory("OpenCASCADE");', /, 'meshscale=1;')
	!20	FORMAT("Point(",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")
	!30	FORMAT("Line(",I<INC>,")={",I7,",",I7,"};")
	!40  FORMAT("Line Loop(",I<INC>,")={",<NV1>(I7,","),I7,"};")
	!50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
	!60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I7,","),I7,"};")
	!61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
	!70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
	!80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
	!90 	FORMAT(I<INC2>,",")
	!100 FORMAT('Physical Surface("',A<LEN1>,'",',I<INC>,")={",A<INC2>"};")
 !   
 !   CONTAINS
 !       !FUNCTION ORIENT_EDGE()
 !       FUNCTION EDGE_NODE(IEDGE,ILAYER) RESULT(V)
 !           INTEGER,INTENT(IN)::IEDGE,ILAYER !O<=ILAYER<=SOILLAYER
 !           INTEGER::V(2)           
 !           
 !           V=EDGE(IEDGE).V+NNODE*(ILAYER)
 !       
 !   
 !       ENDFUNCTION
 !       FUNCTION VEDGE_NODE(INODE,ISOIL) RESULT(V)
 !           INTEGER,INTENT(IN)::INODE,ISOIL !ISOIL>=1
 !           INTEGER::V(2)           
 !           
 !           V=[INODE+NNODE*(ISOIL-1),INODE+NNODE*(ISOIL)]
 !       
 !   
 !       ENDFUNCTION
	!ENDSUBROUTINE	 


end module geomodel



